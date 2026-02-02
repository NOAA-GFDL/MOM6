# MOM6 Solo Driver
#
# This builds the standalone ocean-only MOM6 executable.
# It provides surface forcing routines for idealized and forced simulations.

add_executable(MOM6)

target_link_libraries(MOM6 PRIVATE MOM6::core)

target_sources(MOM6 PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_driver.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_surface_forcing.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MESO_surface_forcing.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/user_surface_forcing.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/atmos_ocean_fluxes.F90
)

# Installation
include(GNUInstallDirs)

install(TARGETS MOM6
  EXPORT MOM6Targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
