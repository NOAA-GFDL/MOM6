# Makefile template for FMS
#
# Compiler flags are configured by autoconf's configure script.
#
# Source code dependencies are configured by makedep and saved to Makefile.dep.

CC = mpicc
FC = mpifort
LD = mpifort
AR = ar
PYTHON = /opt/homebrew/Caskroom/miniforge/base/envs/mommy2/bin/python
MAKEDEP = /Users/aakash/mom6_files/ML_diffusivity_development/MOM6-examples/src/MOM6/ac/makedep

DEFS = -DPACKAGE_NAME=\"FMS\" -DPACKAGE_TARNAME=\"fms\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"FMS\" -DPACKAGE_BUGREPORT=\"https://github.com/NOAA-GFDL/FMS/issues\" -DPACKAGE_URL=\"\" -DHAVE_MPI=1 -D__NR_gettid=SYS_gettid -D__APPLE__=1 -DHAVE_MPI=1 -Duse_libMPI=1 -DHAVE_STDIO_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_STRINGS_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_UNISTD_H=1 -DSTDC_HEADERS=1 -DHAVE_NETCDF_H=1 -Duse_netCDF=1
CPPFLAGS = -D_FORTIFY_SOURCE=2 -isystem /opt/homebrew/Caskroom/miniforge/base/envs/mommy2/include
FCFLAGS = -g -O0 -fcray-pointer -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fallow-invalid-boz 
LDFLAGS = -Wl,-pie -Wl,-headerpad_max_install_names -Wl,-dead_strip_dylibs -Wl,-rpath,/opt/homebrew/Caskroom/miniforge/base/envs/mommy2/lib -L/opt/homebrew/Caskroom/miniforge/base/envs/mommy2/lib
LIBS = 
ARFLAGS = rv

-include Makefile.dep

.PHONY: depend
depend: Makefile.dep
Makefile.dep:
	$(PYTHON) $(MAKEDEP) $(DEFS) -o Makefile.dep -e -x libFMS.a -s ../src/test_fms ../src
