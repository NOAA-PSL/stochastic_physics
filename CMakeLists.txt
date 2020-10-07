if(32BIT)
  message ("Force 64 bits in stochastic_physics")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    if(REPRO)
      string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    else()
      string (REPLACE "-i4 -real-size 32" "-i4 -real-size 64 -no-prec-div -no-prec-sqrt" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    endif()
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")
  endif()
endif()

list(APPEND _stoch_phys_srcs
  kinddef.F90
  mpi_wrapper.F90
  halo_exchange.fv3.F90
  plumes.f90

  stochy_gg_def.f
  stochy_resol_def.f
  stochy_layout_lag.f
  four_to_grid_stochy.F
  glats_stochy.f
  sumfln_stochy.f
  gozrineo_stochy.f
  num_parthds_stochy.f
  get_ls_node_stochy.f
  get_lats_node_a_stochy.f
  setlats_a_stochy.f
  setlats_lag_stochy.f
  epslon_stochy.f
  getcon_lag_stochy.f
  pln2eo_stochy.f
  dozeuv_stochy.f
  dezouv_stochy.f
  mersenne_twister.F

  spectral_layout.F90
  getcon_spectral.F90
  stochy_namelist_def.F90
  compns_stochy.F90
  stochy_internal_state_mod.F90
  stochastic_physics.F90
  stochy_patterngenerator.F90
  stochy_data_mod.F90
  get_stochy_pattern.F90
  initialize_spectral_mod.F90
  cellular_automata_global.F90
  cellular_automata_sgs.F90
  update_ca.F90
  lndp_apply_perts.F90
)

add_library(stochastic_physics ${_stoch_phys_srcs})
add_dependencies(stochastic_physics fms)

target_compile_definitions(stochastic_physics PRIVATE INTERNAL_FILE_NML)

target_link_libraries(stochastic_physics PUBLIC sp::sp_d)
target_link_libraries(stochastic_physics PUBLIC fms)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(stochastic_physics PUBLIC OpenMP::OpenMP_Fortran)
endif()
