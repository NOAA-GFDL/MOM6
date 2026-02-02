# MOM6 CMake Build

## Quick Start

```bash
git submodule init 
git submodule update 
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j
```

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `MOM6_SOLO` | ON | Build solo driver executable |
| `MOM6_ASYMMETRIC` | OFF | Use asymmetric memory layout |
| `MOM6_OPENMP` | OFF | Enable OpenMP |
| `MOM6_FMS1` | OFF | Use FMS1 I/O (legacy) |
| `FMS_GIT_TAG` | 2026.01 | FMS version to fetch |
| `FMS_Fortran_FLAGS` | "" | Extra flags for FMS compilation |

## Dependencies

- MPI (required)
- NetCDF C and Fortran (required)
- FMS (fetched automatically)

If NetCDF isn't found automatically:
```bash
cmake .. -DCMAKE_PREFIX_PATH=/path/to/netcdf
```

## With .testing Makefile

```bash
cd .testing

# Build and test with CMake
make test USE_CMAKE=1

# Use existing CMake build
make test CMAKE_BUILD_DIR=../build

# Pass NetCDF path
make test USE_CMAKE=1 CMAKE_FLAGS="-DCMAKE_PREFIX_PATH=/path/to/netcdf"
```

## Structure

- `cmake/` - CMake modules (compiler flags, FMS fetching, etc.)
- `src/*/CMakeLists.txt` - Source files per subdirectory
- `config_src/drivers/solo_driver/` - Standalone executable

FMS is fetched via FetchContent and built as part of the project.
