#!/bin/bash

# patch MOM_remapping.F90 to expose everything as public (otherwise
# the interface can't access anything)
sed -e 's/private//g' MOM6/src/ALE/MOM_remapping.F90 > MOM_remapping.F90

# build the original fortran source to modules
# this is a pretty ugly way of doing it, but it works
# we've hardcoded the source files in an order that satisfies dependencies
# ideally we'd use a proper build system (ninja, make)
# but hopefully we don't touch the actual remapping core

core_files="MOM_types.F90"

remap_files="MOM_error_handler.F90 MOM6/src/framework/MOM_string_functions.F90 MOM6/src/ALE/regrid_solvers.F90 MOM6/src/ALE/polynomial_functions.F90 MOM6/src/ALE/regrid_edge_slopes.F90 MOM6/src/ALE/regrid_edge_values.F90 MOM6/src/ALE/PCM_functions.F90 MOM6/src/ALE/PLM_functions.F90 MOM6/src/ALE/PPM_functions.F90 MOM6/src/ALE/PQM_functions.F90 MOM_remapping.F90"

eos_files="MOM6/src/equation_of_state/TEOS10/gsw_mod_toolbox.f90 MOM6/src/equation_of_state/TEOS10/gsw_mod_kinds.f90 MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_exact.f90 MOM6/src/equation_of_state/TEOS10/gsw_specvol_second_derivatives.f90 MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_poly.f90 MOM6/src/equation_of_state/TEOS10/gsw_ct_from_t.f90 MOM6/src/equation_of_state/TEOS10/gsw_specvol.f90 MOM6/src/equation_of_state/TEOS10/gsw_chem_potential_water_t_exact.f90 MOM6/src/equation_of_state/TEOS10/gsw_pt_from_ct.f90 MOM6/src/equation_of_state/TEOS10/gsw_pt_from_t.f90 MOM6/src/equation_of_state/TEOS10/gsw_sp_from_sr.f90 MOM6/src/equation_of_state/TEOS10/gsw_gibbs.f90 MOM6/src/equation_of_state/TEOS10/gsw_entropy_part.f90 MOM6/src/equation_of_state/TEOS10/gsw_mod_specvol_coefficients.f90 MOM6/src/equation_of_state/TEOS10/gsw_ct_from_pt.f90 MOM6/src/equation_of_state/TEOS10/gsw_rho.f90 MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_poly.f90 MOM6/src/equation_of_state/TEOS10/gsw_t_from_ct.f90 MOM6/src/equation_of_state/TEOS10/gsw_rho_first_derivatives.f90 MOM6/src/equation_of_state/TEOS10/gsw_pt0_from_t.f90 MOM6/src/equation_of_state/TEOS10/gsw_specvol_first_derivatives.f90  MOM6/src/equation_of_state/TEOS10/gsw_sr_from_sp.f90 MOM6/src/equation_of_state/TEOS10/gsw_gibbs_ice.f90 MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_exact.f90 MOM6/src/equation_of_state/TEOS10/gsw_mod_freezing_poly_coefficients.f90 MOM6/src/equation_of_state/TEOS10/gsw_rho_second_derivatives.f90 MOM6/src/equation_of_state/TEOS10/gsw_t_deriv_chem_potential_water_t_exact.f90 MOM6/src/equation_of_state/TEOS10/gsw_gibbs_pt0_pt0.f90 MOM6/src/equation_of_state/TEOS10/gsw_entropy_part_zerop.f90 MOM6/src/equation_of_state/TEOS10/gsw_mod_gibbs_ice_coefficients.f90 MOM6/src/equation_of_state/TEOS10/gsw_mod_teos10_constants.f90 MOM6/src/equation_of_state/MOM_TFreeze.F90 MOM6/src/equation_of_state/MOM_EOS_linear.F90 MOM6/src/equation_of_state/MOM_EOS_NEMO.F90 MOM6/src/equation_of_state/MOM_EOS_TEOS10.F90 MOM6/src/equation_of_state/MOM_EOS_UNESCO.F90 MOM6/src/equation_of_state/MOM_EOS_Wright.F90 mod_mom6/MOM_EOS.F90"
ndiff_files="MOM6/src/tracer/MOM_neutral_diffusion_aux.F90"

fflags="-fdefault-double-8 -fdefault-real-8 -Waliasing -ffree-line-length-none -fno-range-check -O3 -fPIC"
#fflags="-fdefault-real-8 -ffree-line-length-none -O0 -fPIC -g -fbacktrace"
#fflags="-fdefault-real-8 -ffree-line-length-none -fPIC"
gfortran ${fflags} -I. -c ${core_files} ${remap_files} ${eos_files} ${ndiff_files} neutral_diffusion_driver.F90
gfortran ${fflags} -I. -c ${core_files} ${remap_files} ${eos_files} ${ndiff_files} neutral_diffusion_driver.F90 neutral_aux.F90

# if the API has changed, regenerate (note that changes are required!)
f90wrap -m neutral_diffusion neutral_diffusion_driver.F90 
## build the interface
f2py --f90flags="${fflags}" -c -m _neutral_diffusion f90wrap_*.f90 *.o
#f2py --f90flags="${fflags}" -c -m neutral_diffusion neutral_diffusion_driver.F90 *.o
