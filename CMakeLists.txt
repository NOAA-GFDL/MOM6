cmake_minimum_required(VERSION 3.22)

project(MOM6
  VERSION 6.0.0
  LANGUAGES C Fortran
  DESCRIPTION "GFDL's Modular Ocean Model version 6"
)

# Prevent in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not allowed. Please create a separate build directory.")
endif()

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()

#------------------------------------------------------------------------------
# Build Options
#------------------------------------------------------------------------------

option(MOM6_SOLO              "Build solo driver executable"       ON)
option(MOM6_ASYMMETRIC        "Use asymmetric memory layout"       OFF)
option(MOM6_FMS1              "Use FMS1 I/O infrastructure"        OFF)
option(MOM6_OPENMP            "Enable OpenMP threading"            OFF)
option(MOM6_ACCESS            "Use find_package for FMS (spack)"   OFF)

# FMS configuration (only used when MOM6_ACCESS=OFF)
set(FMS_GIT_URL "https://github.com/NOAA-GFDL/FMS.git"
    CACHE STRING "FMS Git repository URL")
set(FMS_GIT_TAG "2026.01"
    CACHE STRING "FMS Git tag/branch/commit")
set(FMS_Fortran_FLAGS ""
    CACHE STRING "Additional Fortran flags for FMS compilation")

#------------------------------------------------------------------------------
# CMake modules
#------------------------------------------------------------------------------

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Apply compiler flags (MetalQuicha-style)
include(MOM6CompilerFlags)

#------------------------------------------------------------------------------
# System dependencies
#------------------------------------------------------------------------------

# MPI is required
find_package(MPI REQUIRED COMPONENTS C Fortran)

# NetCDF is required for I/O
find_package(NetCDF REQUIRED COMPONENTS C Fortran)

# OpenMP is optional
if(MOM6_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

#------------------------------------------------------------------------------
# Fetch and build FMS
#------------------------------------------------------------------------------

include(FetchFMS)
fetch_fms(URL ${FMS_GIT_URL} TAG ${FMS_GIT_TAG})

#------------------------------------------------------------------------------
# Submodules
#------------------------------------------------------------------------------

set(CVMIX_BUILD_DRIVER OFF CACHE BOOL "" FORCE)
set(CVMIX_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
add_subdirectory(pkg/CVMix-src)

# Generate CMakeLists.txt if submodule lacks CMake support (temporary)
# TODO: udpate the gitsubmodule to use CMake
include(GenerateGSWCMake)
add_subdirectory(pkg/GSW-Fortran)

#------------------------------------------------------------------------------
# MOM6 core library
#------------------------------------------------------------------------------

add_subdirectory(src)

#------------------------------------------------------------------------------
# Drivers
#------------------------------------------------------------------------------

add_subdirectory(config_src/drivers)

#------------------------------------------------------------------------------
# Installation
#------------------------------------------------------------------------------

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Generate package config files
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/MOM6Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/MOM6Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOM6
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/MOM6ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# Install cmake config files
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/MOM6Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/MOM6ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOM6
)

# Export targets
install(EXPORT MOM6Targets
  FILE MOM6Targets.cmake
  NAMESPACE MOM6::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MOM6
)

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------

message(STATUS "")
message(STATUS "=== MOM6 Build Configuration ===")
message(STATUS "  Version:            ${PROJECT_VERSION}")
message(STATUS "  Build type:         ${CMAKE_BUILD_TYPE}")
message(STATUS "  Fortran compiler:   ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  Memory layout:      ${MOM6_ASYMMETRIC}")
message(STATUS "  FMS infrastructure: ${MOM6_FMS1}")
message(STATUS "  OpenMP:             ${MOM6_OPENMP}")
message(STATUS "  Solo driver:        ${MOM6_SOLO}")
message(STATUS "  ACCESS mode:        ${MOM6_ACCESS}")
if(NOT MOM6_ACCESS)
  message(STATUS "  FMS URL:            ${FMS_GIT_URL}")
  message(STATUS "  FMS tag:            ${FMS_GIT_TAG}")
endif()
message(STATUS "================================")
message(STATUS "")
