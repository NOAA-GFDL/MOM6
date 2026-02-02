target_sources(mom6_core PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_ice_shelf.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_ice_shelf_diag_mediator.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_ice_shelf_dynamics.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_ice_shelf_initialize.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_ice_shelf_state.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/MOM_marine_ice.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/user_shelf_init.F90
)
