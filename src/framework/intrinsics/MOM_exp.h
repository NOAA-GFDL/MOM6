!// MOM_exp.h - Preprocessor macros for reproducible exp()
!//
!// Usage:
!//   #include "MOM_exp.h"
!//
!// Required:
!//   use ieee_arithmetic, only : ieee_rint
!//
!// Build flags:
!//   -DENABLE_FAST_RINT : Use fast_rint() bit manipulation
!//   (default)          : Use ieee_rint() fallback

#ifdef ENABLE_FAST_RINT
!// Fast round-to-nearest-integer using bit shifts
#define NEAREST_INT(x) fast_rint(x)
#else
!// Use IEEE's rint() method
#define NEAREST_INT(x) ieee_rint(x)
#endif
