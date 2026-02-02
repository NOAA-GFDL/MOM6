# MOM6 Drivers
#
# This file controls which driver executables are built.
# - solo_driver: Standalone MOM6 ocean-only executable

# Solo driver (ocean-only standalone executable)
if(MOM6_SOLO)
  add_subdirectory(solo_driver)
endif()
