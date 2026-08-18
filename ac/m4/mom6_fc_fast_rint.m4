dnl Determine the flag required to honor parentheses in floating-point
dnl expressions (i.e., prevent reassociation of (a + b) - b into a).
dnl
dnl This is required for the fast round-to-nearest-integer trick:
dnl   K = (x + round_bias) - round_bias
dnl where round_bias = 1.5 * 2^52.
dnl
dnl Compiler flags that enable this:
dnl   GCC gfortran:    -fprotect-parens (default with -std=f2008+)
dnl   Intel ifx/ifort: -assume protect_parens
dnl   NVHPC nvfortran: (appears to be default)
dnl
AC_DEFUN([MOM6_FC_FAST_RINT], [
  AC_ARG_ENABLE([fast-rint],
    [AS_HELP_STRING([--disable-fast-rint],
      [do not use fast round-to-nearest-integer])])
  AS_IF([test "x$enable_fast_rint" != xno], [
    AC_CACHE_CHECK([for $FC option to honor parentheses],
        [mom6_cv_fc_fast_rint], [
          mom6_cv_fc_fast_rint="unsupported"
          ac_fc_fr_FCFLAGS_save=${FCFLAGS}
          AC_LANG_PUSH([Fortran])
          for ac_flag in none \
            -fprotect-parens \
            "-assume protect_parens"
          do
            test "$ac_flag" != none \
              && FCFLAGS="$ac_fc_fr_FCFLAGS_save $ac_flag"
            AC_RUN_IFELSE([
              AC_LANG_PROGRAM([], [
dnl ---
      real, parameter :: bias = 1.5 * 2.**(digits(1.) - 1)
      if (fast_rint(1.23) == 1.) stop 0
      stop 1
      contains
      function fast_rint(x) result(y)
        real, intent(in) :: x
        real :: y
        y = (x + bias) - bias
      end function
dnl ---
              ])
            ], [
              mom6_cv_fc_fast_rint="$ac_flag"
              break
            ], [], [
              mom6_cv_fc_fast_rint="unsupported"
              break
            ])
          done
          AC_LANG_POP([Fortran])
          FCFLAGS=$ac_fc_fr_FCFLAGS_save
        ]
      )
      AS_CASE([$mom6_cv_fc_fast_rint],
        [none], [
          mom6_cv_fc_fast_rint="none needed"
          FAST_RINT_FCFLAGS=""
          AC_DEFINE([ENABLE_FAST_RINT], [1],
            [Define to 1 if parentheses are protected in FP expressions])
        ],
        [unsupported], [
          AC_MSG_WARN(
            [No known flag found to protect parentheses; using ieee_rint fallback])
          FAST_RINT_FCFLAGS=""
        ],
        [
          FAST_RINT_FCFLAGS=$mom6_cv_fc_fast_rint
          AC_DEFINE([ENABLE_FAST_RINT], [1],
            [Define to 1 if parentheses are protected in FP expressions])
        ]
      )
  ])
  AC_SUBST([FAST_RINT_FCFLAGS])
])
