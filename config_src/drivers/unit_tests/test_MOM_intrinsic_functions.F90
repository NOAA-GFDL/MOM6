! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Driver for MOM_intrinsic_functions unit tests
program test_MOM_intrinsic_functions

use MOM_error_handler, only : set_skip_mpi
use MOM_intrinsic_functions, only : intrinsic_functions_unit_tests
use MOM_intrinsic_functions_tests, only : run_intrinsic_functions_tests

call set_skip_mpi(.true.)

if (intrinsic_functions_unit_tests(.true.)) stop 1

call run_intrinsic_functions_tests

end program test_MOM_intrinsic_functions
