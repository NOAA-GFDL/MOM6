# MOM6 Core Library
#
# This file builds the main MOM6 library (mom6_core) containing all
# the physics, dynamics, and framework code.

#------------------------------------------------------------------------------
# Memory and infrastructure configuration
#------------------------------------------------------------------------------

# Select memory configuration directory
if(MOM6_ASYMMETRIC)
  set(MEMORY_DIR "${CMAKE_SOURCE_DIR}/config_src/memory/dynamic_nonsymmetric")
else()
  set(MEMORY_DIR "${CMAKE_SOURCE_DIR}/config_src/memory/dynamic_symmetric")
endif()

# Select FMS infrastructure directory
if(MOM6_FMS1)
  set(INFRA_DIR "${CMAKE_SOURCE_DIR}/config_src/infra/FMS1")
else()
  set(INFRA_DIR "${CMAKE_SOURCE_DIR}/config_src/infra/FMS2")
endif()

# External stubs directory
set(EXTERNAL_DIR "${CMAKE_SOURCE_DIR}/config_src/external")

#------------------------------------------------------------------------------
# Create the core library
#------------------------------------------------------------------------------

add_library(mom6_core STATIC)
add_library(MOM6::core ALIAS mom6_core)

# Module output directory
set(MOM6_MOD_DIR "${CMAKE_CURRENT_BINARY_DIR}/mod")
set_target_properties(mom6_core PROPERTIES
  Fortran_MODULE_DIRECTORY "${MOM6_MOD_DIR}"
  POSITION_INDEPENDENT_CODE ON
)

# Include directories for header files (.h) and Fortran modules (.mod)
target_include_directories(mom6_core PUBLIC
  $<BUILD_INTERFACE:${MEMORY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/framework>
  $<BUILD_INTERFACE:${MOM6_MOD_DIR}>
  $<INSTALL_INTERFACE:include/mom6>
)

# CVMix module directory (CVMix doesn't properly export BUILD_INTERFACE includes)
target_include_directories(mom6_core PRIVATE
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
)

# Link dependencies
# Note: FMS::fms_r8 is an alias to FMS::fms created in FetchFMS.cmake
target_link_libraries(mom6_core PUBLIC
  FMS::fms_r8
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran
)

target_link_libraries(mom6_core PRIVATE
  cvmix_static
  gsw
)

if(MOM6_OPENMP)
  target_link_libraries(mom6_core PUBLIC OpenMP::OpenMP_Fortran)
endif()

#------------------------------------------------------------------------------
# Source files organized by subsystem
#------------------------------------------------------------------------------

add_subdirectory(ALE)
add_subdirectory(core)
add_subdirectory(diagnostics)
add_subdirectory(equation_of_state)
add_subdirectory(framework)
add_subdirectory(ice_shelf)
add_subdirectory(initialization)
add_subdirectory(ocean_data_assim)
add_subdirectory(parameterizations)
add_subdirectory(tracer)
add_subdirectory(user)

#------------------------------------------------------------------------------
# Infrastructure sources (FMS1 or FMS2)
#------------------------------------------------------------------------------

target_sources(mom6_core PRIVATE
  ${INFRA_DIR}/MOM_coms_infra.F90
  ${INFRA_DIR}/MOM_constants.F90
  ${INFRA_DIR}/MOM_couplertype_infra.F90
  ${INFRA_DIR}/MOM_cpu_clock_infra.F90
  ${INFRA_DIR}/MOM_data_override_infra.F90
  ${INFRA_DIR}/MOM_diag_manager_infra.F90
  ${INFRA_DIR}/MOM_domain_infra.F90
  ${INFRA_DIR}/MOM_ensemble_manager_infra.F90
  ${INFRA_DIR}/MOM_error_infra.F90
  ${INFRA_DIR}/MOM_interp_infra.F90
  ${INFRA_DIR}/MOM_io_infra.F90
  ${INFRA_DIR}/MOM_time_manager.F90
)

#------------------------------------------------------------------------------
# External stubs (MARBL, drifters, etc.)
#------------------------------------------------------------------------------

target_sources(mom6_core PRIVATE
  # MARBL stubs
  ${EXTERNAL_DIR}/MARBL/marbl_constants_mod.F90
  ${EXTERNAL_DIR}/MARBL/marbl_interface.F90
  ${EXTERNAL_DIR}/MARBL/marbl_interface_public_types.F90
  ${EXTERNAL_DIR}/MARBL/marbl_logging.F90

  # Drifter stubs
  ${EXTERNAL_DIR}/drifters/MOM_particles.F90
  ${EXTERNAL_DIR}/drifters/MOM_particles_types.F90

  # GFDL BGC stub
  ${EXTERNAL_DIR}/GFDL_ocean_BGC/MOM_generic_tracer.F90

  # Database communications stub
  ${EXTERNAL_DIR}/database_comms/database_client_interface.F90
  ${EXTERNAL_DIR}/database_comms/MOM_database_comms.F90

  # Stochastic physics stubs
  ${EXTERNAL_DIR}/stochastic_physics/get_stochy_pattern.F90
  ${EXTERNAL_DIR}/stochastic_physics/stochastic_physics.F90

  # ODA hooks stubs
  ${EXTERNAL_DIR}/ODA_hooks/kdtree.f90
  ${EXTERNAL_DIR}/ODA_hooks/ocean_da_core.F90
  ${EXTERNAL_DIR}/ODA_hooks/ocean_da_types.F90
  ${EXTERNAL_DIR}/ODA_hooks/write_ocean_obs.F90
)

#------------------------------------------------------------------------------
# Installation
#------------------------------------------------------------------------------

include(GNUInstallDirs)

install(TARGETS mom6_core
  EXPORT MOM6Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY ${MOM6_MOD_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mom6
  FILES_MATCHING PATTERN "*.mod"
)
