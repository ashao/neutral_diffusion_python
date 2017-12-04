module neutral_diffusion_driver

use MOM_EOS,                   only : EOS_type, EOS_manual_init, calculate_density_derivs
use MOM_EOS,                   only : EOS_LINEAR, EOS_TEOS10, EOS_WRIGHT
use MOM_neutral_diffusion_aux, only : mark_unstable_cells, mark_unstable_cells_i, refine_nondim_position
use MOM_neutral_diffusion_aux, only : calc_delta_rho, check_neutral_positions
use MOM_remapping,             only : remapping_CS, initialize_remapping
use MOM_remapping,             only : extract_member_remapping_CS, build_reconstructions_1d, average_value_ppoly
use polynomial_functions,      only : evaluation_polynomial, first_derivative_polynomial
use MOM_error_handler,         only : MOM_error, WARNING, FATAL
use neutral_aux,               only : find_neutral_surface_positions_discontinuous

implicit none

public neutral_diffusion_init
public find_neutral_surfaces

type neutral_diffusion_CS
  integer :: ppoly_deg = 2 ! Degree of polynomial used for reconstructions
  logical :: continuous_reconstruction = .true.   ! True if using continuous PPM reconstruction at interfaces
  logical :: boundary_extrap = .true.
  logical :: refine_pos = .false.
  integer :: max_iter ! Maximum number of iterations if refine_position is defined
  double precision :: tolerance   ! Convergence criterion representing difference from true neutrality
  double precision :: ref_pres    ! Reference pressure, negative if using locally referenced neutral density
  character(len=80) :: EOS_name
  double precision :: dRho_dT
  double precision :: dRho_dS

  character(len=80) :: remapping_scheme
end type neutral_diffusion_CS

contains

function neutral_diffusion_init(ref_pres, remapping_scheme, boundary_extrap, refine_position, tolerance, &
                                max_iter, EOS_name, drho_dT, drho_dS) result (CS)
  double precision,                       intent(in   ) :: ref_pres
  character(len=80),          intent(in   ) :: remapping_scheme
  logical,                    intent(in   ) :: boundary_extrap
  logical,                    intent(in   ) :: refine_position
  double precision,                       intent(in   ) :: tolerance
  integer,                    intent(in   ) :: max_iter
  character(len=80),          intent(in   ) :: EOS_name
  double precision, optional,             intent(in   ) :: drho_dT
  double precision, optional,             intent(in   ) :: drho_dS

  type(neutral_diffusion_CS) :: CS

  CS%continuous_reconstruction = .false.
  ! Remapping
  CS%boundary_extrap = boundary_extrap
  ! Related to finding neutral surface
  CS%ref_pres = ref_pres
  CS%refine_pos = refine_position
  CS%tolerance = tolerance
  CS%max_iter = max_iter
  CS%EOS_name = EOS_name
  if (present(drho_dT)) CS%drho_dT = drho_DT
  if (present(drho_dS)) CS%drho_dS = drho_DS
  CS%remapping_scheme = remapping_scheme

end function neutral_diffusion_init

subroutine find_neutral_surfaces(CS, nk, Tl, Sl, hl, Tr, Sr, hr, PoL, PoR, KoL, KoR, hEff)
  type(neutral_diffusion_CS), intent(inout) :: CS
  integer,                    intent(in   ) :: nk
  double precision, dimension(nk),        intent(in   ) :: Tl
  double precision, dimension(nk),        intent(in   ) :: Sl
  double precision, dimension(nk),        intent(in   ) :: hl
  double precision, dimension(nk),        intent(in   ) :: Tr
  double precision, dimension(nk),        intent(in   ) :: Sr
  double precision, dimension(nk),        intent(in   ) :: hr
  double precision, dimension(4*nk),      intent(  out) :: PoL
  double precision, dimension(4*nk),      intent(  out) :: PoR
  integer, dimension(4*nk),   intent(  out) :: KoL
  integer, dimension(4*nk),   intent(  out) :: KoR
  double precision, dimension(4*nk-1),    intent(  out) :: hEff

  ! Local variables
  type(EOS_type), target  :: EOS_main
  type(EOS_type), pointer :: EOS
  type(remapping_CS)      :: remap_CS
  integer :: k
  integer :: ns_l, ns_r
  integer :: iMethod
  double precision, dimension(nk,2)  :: Ti_l, Si_l, dRdT_il, dRdS_il
  double precision, dimension(nk,2)  :: Ti_r, Si_r, dRdT_ir, dRdS_ir
  double precision, dimension(nk,2)  :: ppoly_Tl, ppoly_Sl
  double precision, dimension(nk,2)  :: ppoly_Tr, ppoly_Sr
  double precision, dimension(nk,2)  :: ppoly_r_S
  double precision, dimension(nk+1)  :: Pres_l, Pres_r
  double precision, dimension(nk)    :: Play_l, Play_r, ref_pres
  double precision, dimension(nk)    :: dRdT_ll, dRdS_ll, dRdT_lr, dRdS_lr
  logical, dimension(nk) :: stable_cell_l, stable_cell_r
  print *, "Initializing"
  ! Initialize equation of state
  EOS => EOS_main
  if (trim(CS%EOS_name) == 'LINEAR') then
    call EOS_manual_init(EOS, form_of_EOS = EOS_LINEAR, drho_dT = CS%drho_dT, dRho_dS = CS%drho_dS)
  elseif (trim(CS%EOS_name) == 'WRIGHT') then
    call EOS_manual_init(EOS, form_of_EOS = EOS_WRIGHT)
  elseif (trim(CS%EOS_name) == 'TEOS10') then
    call EOS_manual_init(EOS, form_of_EOS = EOS_TEOS10)
  endif

  ! Initialize remapping
  call initialize_remapping( remap_CS, trim(CS%remapping_scheme), boundary_extrapolation = CS%boundary_extrap )
  call extract_member_remapping_CS(remap_CS, degree=CS%ppoly_deg)

  ! Calculate pressure at interfaces and mid-layer
  Pres_l(:) = 0.
  Pres_r(:) = 0.
  do k = 2, nk
    Pres_l(k) = Pres_l(k-1) + hl(k)
    Play_l(k-1) = 0.5*(Pres_l(k-1) + Pres_l(k))
    Pres_r(k) = Pres_r(k-1) + hr(k)
    Play_r(k-1) = 0.5*(Pres_r(k-1) + Pres_r(k))
  enddo
  print *, "Calculated pressures"
  if (CS%ref_pres < 0.) then
    call calculate_density_derivs(Tl, Sl, Play_l, dRdT_ll, dRdS_ll, 1, nk, EOS)
    call calculate_density_derivs(Tr, Sr, Play_r, dRdT_lr, dRdS_lr, 1, nk, EOS)
  else
    ref_pres(:) = CS%ref_pres
    call calculate_density_derivs(Tl, Sl, ref_pres, dRdT_ll, dRdS_ll, 1, nk, EOS)
    call calculate_density_derivs(Tr, Sr, ref_pres, dRdT_lr, dRdS_lr, 1, nk, EOS)
  endif
  call mark_unstable_cells( nk, dRdT_ll, dRdS_ll, Tl, Sl, stable_cell_l, ns_l )
  call mark_unstable_cells( nk, dRdT_lr, dRdS_lr, Tr, Sr, stable_cell_r, ns_r )
  print *, "Marked unstable cells"

  call build_reconstructions_1d(remap_CS, nk, hl, Tl, ppoly_Tl, Ti_l, ppoly_r_S, iMethod)
  call build_reconstructions_1d(remap_CS, nk, hl, Sl, ppoly_Sl, Si_l, ppoly_r_S, iMethod)
  call build_reconstructions_1d(remap_CS, nk, hr, Tr, ppoly_Tr, Ti_r, ppoly_r_S, iMethod)
  call build_reconstructions_1d(remap_CS, nk, hr, Sr, ppoly_Sr, Si_r, ppoly_r_S, iMethod)
  print *, "Built reconstrutions"

  if (CS%ref_pres < 0.) then
    ! Left column
    do k = 1,nk ; ref_pres(k) = Pres_l(k) ; enddo
    call calculate_density_derivs(Ti_l(:,1), Si_l(:,1), ref_pres, dRdT_il(:,1), dRdS_il(:,1), 1, nk, EOS)
    do k = 1,nk ; ref_pres(k) = Pres_l(k+1) ; enddo
    call calculate_density_derivs(Ti_l(:,2), Si_l(:,2), ref_pres, dRdT_il(:,2), dRdS_il(:,2), 1, nk, EOS)
    ! Right column
    do k = 1,nk ; ref_pres(k) = Pres_l(k) ; enddo
    call calculate_density_derivs(Ti_r(:,1), Si_r(:,1), ref_pres, dRdT_ir(:,1), dRdS_ir(:,1), 1, nk, EOS)
    do k = 1,nk ; ref_pres(k) = Pres_r(k+1) ; enddo
    call calculate_density_derivs(Ti_r(:,2), Si_r(:,2), ref_pres, dRdT_ir(:,2), dRdS_ir(:,2), 1, nk, EOS)
  else
    ref_pres(:) = CS%ref_pres
    call calculate_density_derivs(Ti_l(:,1), Si_l(:,1), ref_pres, dRdT_il(:,1), dRdS_il(:,1), 1, nk, EOS)
    call calculate_density_derivs(Ti_l(:,2), Si_l(:,2), ref_pres, dRdT_il(:,2), dRdS_il(:,2), 1, nk, EOS)
    call calculate_density_derivs(Ti_r(:,1), Si_r(:,1), ref_pres, dRdT_ir(:,1), dRdS_ir(:,1), 1, nk, EOS)
    call calculate_density_derivs(Ti_r(:,2), Si_r(:,2), ref_pres, dRdT_ir(:,2), dRdS_ir(:,2), 1, nk, EOS)
  endif

  KoL(:) = 0  ; KoR(:) = 0
  PoL(:) = 0. ; PoR(:) = 0.
  hEff(:) = 0.
  print *, "Finding positions", nk, ns_l, ns_r
  call find_neutral_surface_positions_discontinuous(nk, ns_l + ns_r, CS%ppoly_deg,                                    &
                                                    Pres_l, hl, Ti_l, Si_l, dRdT_il, dRdS_il, stable_cell_l,          &
                                                    Pres_r, hr, Ti_r, Si_r, dRdT_ir, dRdS_ir, stable_cell_r,          &
                                                    PoL, PoR, KoL, KoR, hEff, CS%refine_pos, ppoly_Tl, ppoly_Sl,      &
                                                    ppoly_Tr, ppoly_Sr, EOS, CS%max_iter, CS%tolerance, CS%ref_pres)
  EOS => NULL()
  do k=1,4*nk-1
    print *, PoL(k), PoR(k), KoL(k), KoR(k), hEff(k)
  enddo
  print *, "Found positions"
end subroutine find_neutral_surfaces

end module neutral_diffusion_driver
