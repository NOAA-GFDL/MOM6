# FetchFMS.cmake - FetchContent-based FMS acquisition
#
# This module provides a macro to download and build FMS as part of the MOM6 build.
# It follows the MetalQuicha pattern for fetching dependencies.
#
# FMS_Fortran_FLAGS can be set to pass additional compiler flags to FMS.

include(FetchContent)

# Macro to fetch FMS from a Git repository
#
# Usage:
#   fetch_fms(URL <git_url> TAG <git_tag>)
#
# Arguments:
#   URL - Git repository URL (default: https://github.com/NOAA-GFDL/FMS.git)
#   TAG - Git tag, branch, or commit hash (default: 2026.01)
#
# This macro will:
#   1. Configure FMS build options
#   2. Fetch FMS source via FetchContent
#   3. Build FMS as part of the main build
#   4. Provide FMS::fms_r8 target for linking
#
macro(fetch_fms)
  cmake_parse_arguments(FMS "" "URL;TAG" "" ${ARGN})

  # Set defaults if not specified
  if(NOT FMS_URL)
    set(FMS_URL "https://github.com/NOAA-GFDL/FMS.git")
  endif()
  if(NOT FMS_TAG)
    set(FMS_TAG "2026.01")
  endif()

  message(STATUS "Fetching FMS from ${FMS_URL} (tag: ${FMS_TAG})")

  # Apply FMS-specific Fortran flags if provided
  if(FMS_Fortran_FLAGS)
    message(STATUS "FMS additional Fortran flags: ${FMS_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FMS_Fortran_FLAGS}")
  endif()

  # Configure FMS build options before fetching
  # These must be set as cache variables to affect the FMS build
  #
  # NOTE: We do NOT set 64BIT or 32BIT, which causes FMS to build its default
  # double-precision library as FMS::fms. This is required because FMS's own
  # tests link to FMS::fms, which doesn't exist when explicit kind options are set.
  # The default FMS::fms is double precision (r8) with mixed precision support.
  set(OPENMP ${MOM6_OPENMP} CACHE BOOL "Enable OpenMP in FMS" FORCE)

  # Disable FMS install
  set(FMS_INSTALL OFF CACHE BOOL "" FORCE)

  # Declare FMS as a dependency
  FetchContent_Declare(fms
    GIT_REPOSITORY ${FMS_URL}
    GIT_TAG        ${FMS_TAG}
    GIT_SHALLOW    TRUE
    GIT_PROGRESS   TRUE
  )

  # Make FMS available
  FetchContent_MakeAvailable(fms)

  # Create FMS::fms_r8 alias for compatibility with code expecting the r8 target
  if(NOT TARGET FMS::fms_r8)
    add_library(FMS::fms_r8 ALIAS fms)
  endif()

  # FMS creates the following targets when built with 64BIT=ON:
  #   - FMS::fms_r8      (64-bit real precision library)
  #   - FMS::fms_r4      (32-bit real precision library, if built)
  #
  # For MOM6, we use FMS::fms_r8 exclusively.

  message(STATUS "FMS fetched and configured successfully")
endmacro()
