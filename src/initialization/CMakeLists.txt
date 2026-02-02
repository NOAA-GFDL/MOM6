target_sources(mom6_core PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_coord_initialization.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_fixed_initialization.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_grid_initialize.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_shared_initialization.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_state_initialization.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_tracer_initialization_from_Z.F90
)
