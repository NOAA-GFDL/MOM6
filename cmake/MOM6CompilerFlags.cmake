# MOM6 Compiler Flags - MetalQuicha-style category-based flags
#
# This module defines global compiler flags for MOM6 based on the compiler vendor.
# Flags are organized into categories:
#   FLAGS_STANDARD  - Always applied (line length, endianness, real8, etc.)
#   FLAGS_DEBUG     - Debug builds (trapping, bounds checking)
#   FLAGS_RELEASE   - Optimized builds
#   FLAGS_WARNINGS  - Warning flags

# Helper function to join list elements into a string
function(list_to_string LIST_VAR OUTPUT_VAR)
  string(REPLACE ";" " " _result "${${LIST_VAR}}")
  set(${OUTPUT_VAR} "${_result}" PARENT_SCOPE)
endfunction()

# GNU Fortran (gfortran)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(FLAGS_STANDARD
    -ffree-line-length-none
    -fconvert=big-endian
    -fbacktrace
    -fdefault-real-8
    -fdefault-double-8
    -fcray-pointer
  )

  set(FLAGS_DEBUG
    -O0
    -g
    -fcheck=bounds
  )
  # Note: -ffpe-trap=invalid,zero,overflow is too aggressive for FMS
  # Use -fcheck=bounds instead of -fcheck=all for compatibility

  set(FLAGS_RELEASE
    -O2
    -funroll-loops
  )

  set(FLAGS_WARNINGS
    -Wall
    -Wextra
    -Wno-compare-reals
    -Wno-unused-dummy-argument
    -Wno-unused-variable
    -Wno-unused-parameter
  )

  # GCC 10+ requires -fallow-argument-mismatch for legacy MPI interfaces
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0")
    list(APPEND FLAGS_STANDARD -fallow-argument-mismatch -fallow-invalid-boz)
  endif()

# Intel Classic Fortran (ifort)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(FLAGS_STANDARD
    -convert big_endian
    -assume byterecl
    -traceback
    -r8
    -ftz
    -fpp
  )

  set(FLAGS_DEBUG
    -O0
    -g
    -check bounds
  )
  # Note: -fpe0 causes FPE trapping which breaks FMS initialization

  set(FLAGS_RELEASE
    -O2
    -debug minimal
  )

  set(FLAGS_WARNINGS
    -warn all
    -warn nointerfaces
    -warn nounused
  )

# Intel LLVM Fortran (ifx)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  set(FLAGS_STANDARD
    -convert big_endian
    -assume byterecl
    -traceback
    -r8
    -ftz
    -fpp
  )

  set(FLAGS_DEBUG
    -O0
    -g
    -check bounds
  )
  # Note: -fpe0 causes FPE trapping which breaks FMS initialization

  set(FLAGS_RELEASE
    -O2
    -debug minimal
  )

  set(FLAGS_WARNINGS
    -warn all
    -warn nointerfaces
    -warn nounused
  )

# NVIDIA HPC Fortran (nvfortran)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  set(FLAGS_STANDARD
    -Mfree
    -byteswapio
    -traceback
    -r8
    -Mpreprocess
  )

  set(FLAGS_DEBUG
    -O0
    -g
    -Mbounds
    -Mchkptr
    -Ktrap=fp
  )

  set(FLAGS_RELEASE
    -O2
    -fast
  )

  set(FLAGS_WARNINGS
    -Minform=warn
  )

# Cray Fortran (crayftn)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(FLAGS_STANDARD
    -f free
    -h byteswapio
    -s real64
    -eF
  )

  set(FLAGS_DEBUG
    -O0
    -g
    -K trap=fp
    -R bps
  )

  set(FLAGS_RELEASE
    -O2
  )

  set(FLAGS_WARNINGS "")

else()
  message(WARNING "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}. Using default flags.")
  set(FLAGS_STANDARD "")
  set(FLAGS_DEBUG "-O0 -g")
  set(FLAGS_RELEASE "-O2")
  set(FLAGS_WARNINGS "")
endif()

# Convert lists to space-separated strings
list_to_string(FLAGS_STANDARD FLAGS_STANDARD_STR)
list_to_string(FLAGS_DEBUG FLAGS_DEBUG_STR)
list_to_string(FLAGS_RELEASE FLAGS_RELEASE_STR)
list_to_string(FLAGS_WARNINGS FLAGS_WARNINGS_STR)

# Apply flags globally to CMAKE_Fortran_FLAGS
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FLAGS_STANDARD_STR} ${FLAGS_WARNINGS_STR}"
    CACHE STRING "Fortran compiler flags" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG "${FLAGS_DEBUG_STR}"
    CACHE STRING "Fortran compiler flags for Debug builds" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "${FLAGS_RELEASE_STR}"
    CACHE STRING "Fortran compiler flags for Release builds" FORCE)
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${FLAGS_RELEASE_STR} -g"
    CACHE STRING "Fortran compiler flags for RelWithDebInfo builds" FORCE)

# Report the flags being used
message(STATUS "MOM6 Fortran Compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "MOM6 Fortran Standard Flags: ${FLAGS_STANDARD_STR}")
message(STATUS "MOM6 Fortran Warning Flags: ${FLAGS_WARNINGS_STR}")
