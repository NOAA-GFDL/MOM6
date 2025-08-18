program test_MOM_ANN

use MOM_ANN, only : ANN_unit_tests
use MOM_error_handler, only : set_skip_mpi

call set_skip_mpi(.true.) ! This unit tests is not expecting MPI to be used

if ( ANN_unit_tests(.true.) ) stop 1

end program test_MOM_ANN
