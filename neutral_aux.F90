module neutral_aux

use MOM_EOS,                   only : EOS_type, EOS_manual_init, calculate_density_derivs
use MOM_EOS,                   only : EOS_LINEAR, EOS_TEOS10, EOS_WRIGHT
use MOM_neutral_diffusion_aux, only : mark_unstable_cells, mark_unstable_cells_i, refine_nondim_position
use MOM_neutral_diffusion_aux, only : calc_delta_rho, check_neutral_positions
use MOM_remapping,             only : remapping_CS, initialize_remapping
use MOM_remapping,             only : extract_member_remapping_CS, build_reconstructions_1d, average_value_ppoly
use polynomial_functions,      only : evaluation_polynomial, first_derivative_polynomial
use MOM_error_handler,         only : MOM_error, WARNING, FATAL

implicit none

contains

!> Higher order version of find_neutral_surface_positions. Returns positions within left/right columns
!! of combined interfaces using intracell reconstructions of T/S
subroutine find_neutral_surface_positions_discontinuous(nk, ns, deg,                                                   &
                Pres_l, hcol_l, Tl, Sl, dRdT_l, dRdS_l, stable_l, Pres_r, hcol_r, Tr, Sr, dRdT_r, dRdS_r, stable_r,    &
                PoL, PoR, KoL, KoR, hEff,                                                                              &
                refine_pos_in, ppoly_T_l, ppoly_S_l, ppoly_T_r, ppoly_S_r, EOS, max_iter, tolerance, ref_pres)
  integer,                    intent(in)    :: nk        !< Number of levels
  integer,                    intent(in)    :: ns        !< Number of neutral surfaces
  integer,                    intent(in)    :: deg       !< Degree of polynomial used for reconstructions
  real, dimension(nk+1),      intent(in)    :: Pres_l    !< Left-column interface pressure (Pa)
  real, dimension(nk),        intent(in)    :: hcol_l    !< Left-column layer thicknesses
  real, dimension(nk,2),      intent(in)    :: Tl        !< Left-column top interface potential temperature (degC)
  real, dimension(nk,2),      intent(in)    :: Sl        !< Left-column top interface salinity (ppt)
  real, dimension(nk,2),      intent(in)    :: dRdT_l    !< Left-column, top interface dRho/dT (kg/m3/degC)
  real, dimension(nk,2),      intent(in)    :: dRdS_l    !< Left-column, top interface dRho/dS (kg/m3/ppt)
  logical, dimension(nk),     intent(in)    :: stable_l  !< Left-column, top interface dRho/dS (kg/m3/ppt)
  real, dimension(nk+1),      intent(in)    :: Pres_r    !< Right-column interface pressure (Pa)
  real, dimension(nk),        intent(in)    :: hcol_r    !< Left-column layer thicknesses
  real, dimension(nk,2),      intent(in)    :: Tr        !< Right-column top interface potential temperature (degC)
  real, dimension(nk,2),      intent(in)    :: Sr        !< Right-column top interface salinity (ppt)
  real, dimension(nk,2),      intent(in)    :: dRdT_r    !< Right-column, top interface dRho/dT (kg/m3/degC)
  real, dimension(nk,2),      intent(in)    :: dRdS_r    !< Right-column, top interface dRho/dS (kg/m3/ppt)
  logical, dimension(nk),     intent(in)    :: stable_r  !< Left-column, top interface dRho/dS (kg/m3/ppt)
  real, dimension(4*nk),      intent(inout) :: PoL       !< Fractional position of neutral surface within
                                                         !! layer KoL of left column
  real, dimension(4*nk),      intent(inout) :: PoR       !< Fractional position of neutral surface within
                                                         !! layer KoR of right column
  integer, dimension(4*nk),   intent(inout) :: KoL       !< Index of first left interface above neutral surface
  integer, dimension(4*nk),   intent(inout) :: KoR       !< Index of first right interface above neutral surface
  real, dimension(4*nk-1),    intent(inout) :: hEff      !< Effective thickness between two neutral surfaces (Pa)
  logical,                    optional, intent(in) :: refine_pos_in !< True if rootfinding is used for position
  real, dimension(nk,deg+1),  optional, intent(in) :: ppoly_T_l !< Left-column coefficients of T reconstruction
  real, dimension(nk,deg+1),  optional, intent(in) :: ppoly_S_l !< Left-column coefficients of S reconstruction
  real, dimension(nk,deg+1),  optional, intent(in) :: ppoly_T_r !< Right-column coefficients of T reconstruction
  real, dimension(nk,deg+1),  optional, intent(in) :: ppoly_S_r !< Right-column coefficients of S reconstruction
  type(EOS_type),             optional, pointer    :: EOS       !< Equation of state structure
  integer,                    optional, intent(in) :: max_iter  !< Maximum number of iterations in refine_position
  real,                       optional, intent(in) :: tolerance !< Convergence criterion for refine_position
  real,                       optional, intent(in) :: ref_pres  !< Reference pressure to use for deriviative calculation

  ! Local variables
  integer :: k_surface              ! Index of neutral surface
  integer :: kl_left, kl_right      ! Index of layers on the left/right
  integer :: ki_left, ki_right      ! Index of interfaces on the left/right
  logical :: searching_left_column  ! True if searching for the position of a right interface in the left column
  logical :: searching_right_column ! True if searching for the position of a left interface in the right column
  logical :: reached_bottom         ! True if one of the bottom-most interfaces has been used as the target
  logical :: refine_pos             ! Use rootfinding to find the true neutral surface position
  integer :: k, kl_left_0, kl_right_0
  real    :: dRho, dRhoTop, dRhoBot, dRhoTopm1, hL, hR
  integer :: lastK_left, lastK_right
  real    :: lastP_left, lastP_right
  real    :: min_bound
  logical, dimension(nk) :: top_connected_l, top_connected_r
  logical, dimension(nk) :: bot_connected_l, bot_connected_r
  logical :: debug_this_module = .true.
  top_connected_l(:) = .false. ; top_connected_r(:) = .false.
  bot_connected_l(:) = .false. ; bot_connected_r(:) = .false.
  ! Vectors with all the values of the discontinuous reconstruction.
  ! Dimensions are [number of layers x number of interfaces]. Second dimension = 1 for top interface, = 2 for bottom
!  real, dimension(nk,2) :: Sl, Sr, Tl, Tr, dRdT_l, dRdS_l, dRdT_r, dRdS_r

! Check to make sure that polynomial reconstructions were passed if refine_pos defined)
  refine_pos = .false.
  if (present(refine_pos_in)) then
    refine_pos = refine_pos_in
    if (refine_pos .and. (.not. ( present(ppoly_T_l) .and. present(ppoly_S_l) .and.  &
                                  present(ppoly_T_r) .and. present(ppoly_S_r) .and.  &
                                  present(tolerance) .and. present(max_iter) .and. present(ref_pres) ) )) &
        call MOM_error(FATAL, "fine_neutral_surface_positions_discontinuous: refine_pos is requested, but polynomial"// &
                              "coefficients not available for T and S")
  endif

  do k = 1,nk
    if (stable_l(k)) then
      kl_left = k
      kl_left_0 = k
      exit
    endif
  enddo
  do k = 1,nk
    if (stable_r(k)) then
      kl_right = k
      kl_right_0 = k
      exit
    endif
  enddo

  ! Initialize variables for the search
  ki_right = 1 ; lastK_right = 1 ; lastP_right = -1.
  ki_left = 1  ; lastK_left = 1  ; lastP_left = -1.

  reached_bottom = .false.
  searching_left_column = .false.
  searching_right_column = .false.

  ! Loop over each neutral surface, working from top to bottom
  neutral_surfaces: do k_surface = 1, ns
    ! Potential density difference, rho(kr) - rho(kl)
    dRho = 0.5 * &
      ( ( dRdT_r(kl_right,ki_right) + dRdT_l(kl_left,ki_left) ) * ( Tr(kl_right,ki_right) - Tl(kl_left,ki_left) ) &
      + ( dRdS_r(kl_right,ki_right) + dRdS_l(kl_left,ki_left) ) * ( Sr(kl_right,ki_right) - Sl(kl_left,ki_left) ) )
    if (debug_this_module)  write(*,'(A,I2,A,E12.4,A,I2,A,I2,A,I2,A,I2)') "k_surface=",k_surface,"  dRho=",dRho,"  &
        kl_left=",kl_left, "  ki_left=",ki_left,"  kl_right=",kl_right, "  ki_right=",ki_right
    ! Which column has the lighter surface for the current indexes, kr and kl
    if (.not. reached_bottom) then
      if (dRho < 0.) then
        searching_left_column = .true.
        searching_right_column = .false.
      elseif (dRho > 0.) then
        searching_right_column = .true.
        searching_left_column = .false.
      else ! dRho == 0.
        if (  ( kl_left == kl_left_0) .and. ( kl_right == kl_right_0 ) .and. (ki_left + ki_right == 2) ) then ! Still at surface
          searching_left_column = .true.
          searching_right_column = .false.
        else ! Not the surface so we simply change direction
          searching_left_column = .not. searching_left_column
          searching_right_column = .not. searching_right_column
        endif
      endif
    endif

    if (searching_left_column) then
      ! Determine differences between right column interface and potentially three different parts of the left
      ! Potential density difference, rho(kl-1) - rho(kr) (should be negative)
      dRhoTop = 0.5 * &
        ( ( dRdT_l(kl_left,1) + dRdT_r(kl_right,ki_right) ) * ( Tl(kl_left,1) - Tr(kl_right,ki_right) ) &
        + ( dRdS_l(kl_left,1) + dRdS_r(kl_right,ki_right) ) * ( Sl(kl_left,1) - Sr(kl_right,ki_right) ) )
      ! Potential density difference, rho(kl) - rho(kl_right,ki_right) (will be positive)
      dRhoBot = 0.5 * &
        ( ( dRdT_l(kl_left,2) + dRdT_r(kl_right,ki_right) ) * ( Tl(kl_left,2) - Tr(kl_right,ki_right) ) &
        + ( dRdS_l(kl_left,2) + dRdS_r(kl_right,ki_right) ) * ( Sl(kl_left,2) - Sr(kl_right,ki_right) ) )
      if (lastK_left /= kl_left .and. kl_left>kl_left_0) then
        if (stable_l(kl_left-1) ) then ! Calculate the density difference at top of discontinuity
          dRhoTopm1 = 0.5 * &
            ( ( dRdT_l(kl_left-1,2) + dRdT_r(kl_right,ki_right) ) * ( Tl(kl_left-1,2) - Tr(kl_right,ki_right) ) &
            + ( dRdS_l(kl_left-1,2) + dRdS_r(kl_right,ki_right) ) * ( Sl(kl_left-1,2) - Sr(kl_right,ki_right) ) )
        endif
      else
        dRhoTopm1 = dRhoTop
      endif
      if (debug_this_module) then
        write(*,'(A,I2,A,E12.4,A,E12.4,A,E12.4)') "Searching left layer ", kl_left, ":  dRhoTopm1=", dRhoTopm1, &
                                                  "  dRhoTop=", dRhoTop, "  dRhoBot=", dRhoBot
        write(*,'(A,I2,X,I2)') "Searching from right: ", kl_right, ki_right
        write(*,*) "Temp/Salt Reference: ", Tr(kl_right,ki_right), Sr(kl_right,ki_right)
        write(*,*) "Temp/Salt Top L: ", Tl(kl_left,1), Sl(kl_left,1)
        write(*,*) "Temp/Salt Bot L: ", Tl(kl_left,2), Sl(kl_left,2)
      endif

      ! Set the position within the starting column
      PoR(k_surface) = REAL(ki_right-1)
      KoR(k_surface) = kl_right

      ! Set position within the searched column
      call search_other_column_discontinuous(dRhoTopm1, dRhoTop, dRhoBot, Pres_l(kl_left), Pres_l(kl_left+1), &
        lastP_left, lastK_left, kl_left, kl_left_0, ki_left, tolerance, top_connected_l, bot_connected_l, &
        PoL(k_surface), KoL(k_surface))
      if ( refine_pos .and. (PoL(k_surface) > 0.) .and. (PoL(k_surface) < 1.) ) then
        min_bound = 0.
        if (k_surface > 1) then
          if ( KoL(k_surface) == KoL(k_surface-1) ) min_bound = PoL(k_surface-1)
        endif
        PoL(k_surface) = refine_nondim_position(max_iter, tolerance, Tr(kl_right,ki_right), Sr(kl_right,ki_right),     &
                    dRdT_r(kl_right,ki_right), dRdS_r(kl_right,ki_right), Pres_l(kl_left), Pres_l(kl_left+1),          &
                    deg, ppoly_T_l(kl_left,:), ppoly_S_l(kl_left,:), EOS, min_bound, dRhoTop, dRhoBot, min_bound, &
                    ref_pres)
      endif
      if (PoL(k_surface) == 0.) top_connected_l(KoL(k_surface)) = .true.
      if (PoL(k_surface) == 1.) bot_connected_l(KoL(k_surface)) = .true.
      call increment_interface(nk, kl_right, ki_right, stable_r, reached_bottom, searching_right_column, searching_left_column)
      lastK_left = KoL(k_surface)  ; lastP_left = PoL(k_surface)

    elseif (searching_right_column) then
      ! Interpolate for the neutral surface position within the right column, layer krm1
      ! Potential density difference, rho(kr-1) - rho(kl) (should be negative)
      dRhoTop = 0.5 * &
        ( ( dRdT_r(kl_right,1) + dRdT_l(kl_left,ki_left) ) * ( Tr(kl_right,1) - Tl(kl_left,ki_left) ) &
        + ( dRdS_r(kl_right,1) + dRdS_l(kl_left,ki_left) ) * ( Sr(kl_right,1) - Sl(kl_left,ki_left) ) )
      dRhoBot = 0.5 * &
        ( ( dRdT_r(kl_right,2) + dRdT_l(kl_left,ki_left) ) * ( Tr(kl_right,2) - Tl(kl_left,ki_left) ) &
        + ( dRdS_r(kl_right,2) + dRdS_l(kl_left,ki_left) ) * ( Sr(kl_right,2) - Sl(kl_left,ki_left) ) )
      if (lastK_right /= kl_right .and. kl_right>kl_right_0) then
        if(stable_r(kl_right-1)) then
          dRhoTopm1 = 0.5 * &
            ( ( dRdT_r(kl_right-1,2) + dRdT_l(kl_left,ki_left) ) * ( Tr(kl_right-1,2) - Tl(kl_left,ki_left) ) &
            + ( dRdS_r(kl_right-1,2) + dRdS_l(kl_left,ki_left) ) * ( Sr(kl_right-1,2) - Sl(kl_left,ki_left) ) )
        endif
      else
        dRhoTopm1 = dRhoTop
      endif
      if (debug_this_module) then
        write(*,'(A,I2,A,E12.4,A,E12.4,A,E12.4)') "Searching right layer ", kl_right, ":  dRhoTopm1=", dRhoTopm1, &
                                                  "  dRhoTop=", dRhoTop, "  dRhoBot=", dRhoBot
        write(*,'(A,I2,X,I2)') "Searching from left: ", kl_left, ki_left
        write(*,*) "Temp/Salt Reference: ", Tl(kl_left,ki_left), Sl(kl_left,ki_left)
        write(*,*) "Temp/Salt Top R: ", Tr(kl_right,1), Sr(kl_right,1)
        write(*,*) "Temp/Salt Bot R: ", Tr(kl_right,2), Sr(kl_right,2)
      endif
      ! Set the position within the starting column
      PoL(k_surface) = REAL(ki_left-1)
      KoL(k_surface) = kl_left

      ! Set position within the searched column
      call search_other_column_discontinuous(dRhoTopm1, dRhoTop, dRhoBot, Pres_r(kl_right), Pres_r(kl_right+1),        &
         lastP_right, lastK_right, kl_right, kl_right_0, ki_right, tolerance, top_connected_r, bot_connected_r,                   &
         PoR(k_surface), KoR(k_surface))
      if ( refine_pos .and. (PoR(k_surface) > 0. .and. PoR(k_surface) < 1.) ) then
        min_bound = 0.
        if (k_surface > 1) then
          if ( KoR(k_surface) == KoR(k_surface-1) )  min_bound = PoR(k_surface-1)
        endif
        PoR(k_surface) = refine_nondim_position(max_iter, tolerance, Tl(kl_left,ki_left), Sl(kl_left,ki_left),         &
                  dRdT_l(kl_left,ki_left), dRdS_l(kl_left,ki_left), Pres_r(kl_right), Pres_r(kl_right+1),              &
                  deg, ppoly_T_r(kl_right,:), ppoly_S_r(kl_right,:), EOS, min_bound, dRhoTop, dRhoBot, min_bound, &
                  ref_pres)
      endif
      if (PoR(k_surface) == 0.) top_connected_r(KoR(k_surface)) = .true.
      if (PoR(k_surface) == 1.) bot_connected_r(KoR(k_surface)) = .true.
      call increment_interface(nk, kl_left, ki_left, stable_l, reached_bottom, searching_left_column, searching_right_column)
      lastK_right = KoR(k_surface) ; lastP_right = PoR(k_surface)

    else
      stop 'Else what?'
    endif
    lastK_left = KoL(k_surface)  ; lastP_left = PoL(k_surface)
    lastK_right = KoR(k_surface) ; lastP_right = PoR(k_surface)

    if (debug_this_module)  write(*,'(A,I3,A,ES16.6,A,I2,A,ES16.6)') "KoL:", KoL(k_surface), " PoL:", PoL(k_surface), "     KoR:", &
      KoR(k_surface), " PoR:", PoR(k_surface)
    ! Effective thickness
    if (k_surface>1) then
      ! This is useful as a check to make sure that positions are monotonically increasing
      hL = absolute_position(nk,ns,Pres_l,KoL,PoL,k_surface) - absolute_position(nk,ns,Pres_l,KoL,PoL,k_surface-1)
      hR = absolute_position(nk,ns,Pres_r,KoR,PoR,k_surface) - absolute_position(nk,ns,Pres_r,KoR,PoR,k_surface-1)
      print *, hL, hR, hEff(k_surface-1)
      ! In the case of a layer being unstably stratified, may get a negative thickness. Set the previous position
      ! to the current location
      if ( hL<0. .or. hR<0. ) then
        hEff(k_surface-1) = 0.
        call MOM_error(WARNING, "hL or hR is negative")
      elseif ( hL > 0. .and. hR > 0.) then
        hL = (PoL(k_surface) - PoL(k_surface-1))*hcol_l(KoL(k_surface))
        hR = (PoR(k_surface) - PoR(k_surface-1))*hcol_r(KoR(k_surface))
        hEff(k_surface-1) = 2. * ( (hL * hR) / ( hL + hR ) )! Harmonic mean
      else
        hEff(k_surface-1) = 0.
      endif
      print *, hL, hR, hEff(k_surface-1)
    endif
  enddo neutral_surfaces
  ! Check to make sure that neutral surfaces are truly neutral
  if (debug_this_module) then
    do k_surface = 1,ns-1
      if (hEff(k_surface)>0.) then
        kl_left = KoL(k_surface)
        kl_right = KoR(k_surface)
        if (refine_pos) then
          if ( check_neutral_positions(deg, EOS, &
                PoL(k_surface), ppoly_T_l(kl_left,:),  ppoly_S_l(kl_left,:),  (/Pres_l(kl_left),Pres_l(kl_left+1)/),  &
                PoR(k_surface), ppoly_T_r(kl_right,:), ppoly_S_r(kl_right,:), (/Pres_r(kl_right),Pres_r(kl_right+1)/),&
                tolerance, ref_pres) ) then
            call MOM_error(WARNING,"Endpoints of neutral surfaces have different densities")
          endif
        endif
      endif
    enddo
  endif

end subroutine find_neutral_surface_positions_discontinuous

!> Converts non-dimensional position within a layer to absolute position (for debugging)
real function absolute_position(n,ns,Pint,Karr,NParr,k_surface)
  integer, intent(in) :: n            !< Number of levels
  integer, intent(in) :: ns           !< Number of neutral surfaces
  real,    intent(in) :: Pint(n+1)    !< Position of interfaces (Pa)
  integer, intent(in) :: Karr(ns)     !< Index of interface above position
  real,    intent(in) :: NParr(ns)    !< Non-dimensional position within layer Karr(:)

  ! Local variables
  integer :: k_surface, k

  k = Karr(k_surface)
  if (k>n) stop 'absolute_position: k>nk is out of bounds!'
  absolute_position = Pint(k) + NParr(k_surface) * ( Pint(k+1) - Pint(k) )

end function absolute_position

!> Searches the "other" (searched) column for the position of the neutral surface
subroutine search_other_column_discontinuous(dRhoTopm1, dRhoTop, dRhoBot, Ptop, Pbot, lastP, lastK, kl, kl_0, ki, &
                                             tolerance, top_connected, bot_connected, out_P, out_K)
  real,                  intent(in   ) :: dRhoTopm1      !< Density difference across previous interface
  real,                  intent(in   ) :: dRhoTop        !< Density difference across top interface
  real,                  intent(in   ) :: dRhoBot        !< Density difference across top interface
  real,                  intent(in   ) :: Ptop           !< Pressure at top interface
  real,                  intent(in   ) :: Pbot           !< Pressure at bottom interface
  real,                  intent(in   ) :: lastP          !< Last position connected in the searched column
  integer,               intent(in   ) :: lastK          !< Last layer connected in the searched column
  integer,               intent(in   ) :: kl             !< Layer in the searched column
  integer,               intent(in   ) :: kl_0           !< Layer in the searched column
  integer,               intent(in   ) :: ki             !< Interface of the searched column
  real,                  intent(in   ) :: tolerance      !< How close to 0 "neutral" is defined
  logical, dimension(:), intent(inout) :: top_connected  !< True if the top interface was pointed to
  logical, dimension(:), intent(inout) :: bot_connected  !< True if the top interface was pointed to
  real,                  intent(  out) :: out_P          !< Position within searched column
  integer,               intent(  out) :: out_K          !< Layer within searched column

  if (kl > kl_0) then ! Away from top cell
    if (kl == lastK) then ! Searching in the same layer
      if (dRhoTop > tolerance) then
        out_P = max(0.,lastP) ; out_K = kl
      elseif ( ABS(dRhoTop - dRhoBot)<tolerance ) then
        if (top_connected(kl)) then
          out_P = 1. ; out_K = kl
        else
          out_P = max(0.,lastP) ; out_K = kl
        endif
      elseif (dRhoTop >= (dRhoBot+tolerance)) then
        out_P = 1. ; out_K = kl
      else
        out_K = kl
        out_P = max(interpolate_for_nondim_position( dRhoTop, Ptop, dRhoBot, Pbot ),lastP)
      endif
    else ! Searching across the interface
      if (.not. bot_connected(kl-1) ) then
        out_K = kl-1
        out_P = 1.
      else
        out_K = kl
        out_P = 0.
      endif
    endif
  else ! At the top cell
    if (ki == 1) then
      out_P = 0. ; out_K = kl
    elseif (dRhoTop > tolerance) then
      out_P = max(0.,lastP) ; out_K = kl
    elseif ( (dRhoTop - dRhoBot)<tolerance) then
      if (top_connected(kl)) then
        out_P = 1. ; out_K = kl
      else
        out_P = max(0.,lastP) ; out_K = kl
      endif
    elseif (dRhoTop >= (dRhoBot+tolerance)) then
      out_P = 1. ; out_K = kl
    else
      out_K = kl
      out_P = max(interpolate_for_nondim_position( dRhoTop, Ptop, dRhoBot, Pbot ),lastP)
    endif
  endif
end subroutine search_other_column_discontinuous

!> Returns the non-dimensional position between Pneg and Ppos where the
!! interpolated density difference equals zero.
!! The result is always bounded to be between 0 and 1.
real function interpolate_for_nondim_position(dRhoNeg, Pneg, dRhoPos, Ppos)
  real, intent(in) :: dRhoNeg !< Negative density difference
  real, intent(in) :: Pneg    !< Position of negative density difference
  real, intent(in) :: dRhoPos !< Positive density difference
  real, intent(in) :: Ppos    !< Position of positive density difference

  if (Ppos<Pneg) then
    stop 'interpolate_for_nondim_position: Houston, we have a problem! Ppos<Pneg'
  elseif (dRhoNeg>dRhoPos) then
    write(0,*) 'dRhoNeg, Pneg, dRhoPos, Ppos=',dRhoNeg, Pneg, dRhoPos, Ppos
  elseif (dRhoNeg>dRhoPos) then
    stop 'interpolate_for_nondim_position: Houston, we have a problem! dRhoNeg>dRhoPos'
  endif
  if (Ppos<=Pneg) then ! Handle vanished or inverted layers
    interpolate_for_nondim_position = 0.5
  elseif ( dRhoPos - dRhoNeg > 0. ) then
    interpolate_for_nondim_position = min( 1., max( 0., -dRhoNeg / ( dRhoPos - dRhoNeg ) ) )
  elseif ( dRhoPos - dRhoNeg == 0) then
    if (dRhoNeg>0.) then
      interpolate_for_nondim_position = 0.
    elseif (dRhoNeg<0.) then
      interpolate_for_nondim_position = 1.
    else ! dRhoPos = dRhoNeg = 0
      interpolate_for_nondim_position = 0.5
    endif
  else ! dRhoPos - dRhoNeg < 0
    interpolate_for_nondim_position = 0.5
  endif
  if ( interpolate_for_nondim_position < 0. ) stop 'interpolate_for_nondim_position: Houston, we have a problem! Pint < Pneg'
  if ( interpolate_for_nondim_position > 1. ) stop 'interpolate_for_nondim_position: Houston, we have a problem! Pint > Ppos'
end function interpolate_for_nondim_position

!> Increments the interface which was just connected and also set flags if the bottom is reached
subroutine increment_interface(nk, kl, ki, stable, reached_bottom, searching_this_column, searching_other_column)
  integer, intent(in   )                :: nk                     !< Number of vertical levels
  integer, intent(inout)                :: kl                     !< Current layer (potentially updated)
  integer, intent(inout)                :: ki                     !< Current interface
  logical, dimension(nk), intent(in   ) :: stable                 !< True if the cell is stably stratified
  logical, intent(inout)                :: reached_bottom         !< Updated when kl == nk and ki == 2
  logical, intent(inout)                :: searching_this_column  !< Updated when kl == nk and ki == 2
  logical, intent(inout)                :: searching_other_column !< Updated when kl == nk and ki == 2
  integer :: k

  if (ki == 1) then
    ki = 2
  elseif ((ki == 2) .and. (kl < nk) ) then
    do k = kl+1,nk
      if (stable(kl)) then
        kl = k
        ki = 1
        exit
      endif
      ! If we did not find another stable cell, then the current cell is essentially the bottom
      ki = 2
      reached_bottom = .true.
      searching_this_column = .true.
      searching_other_column = .false.
    enddo
  elseif ((kl == nk) .and. (ki==2)) then
    reached_bottom = .true.
    searching_this_column = .true.
    searching_other_column = .false.
  else
    call MOM_error(FATAL,"Unanticipated eventuality in increment_interface")
  endif
end subroutine increment_interface
end module neutral_aux
