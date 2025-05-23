!> Provides functions used with the Piecewise-Parabolic-Method in the vertical ALE algorithm.
module PPM_functions

! This file is part of MOM6. See LICENSE.md for the license.

! First version was created by Laurent White, June 2008.
! Substantially re-factored January 2016.

!! @todo Re-factor PPM_boundary_extrapolation to give round-off safe and
!!       optimization independent results.

use regrid_edge_values, only : bound_edge_values, check_discontinuous_edge_values

implicit none ; private

public PPM_reconstruction, PPM_boundary_extrapolation, PPM_monotonicity

contains

!> Builds quadratic polynomials coefficients from cell mean and edge values.
subroutine PPM_reconstruction( N, h, u, edge_values, ppoly_coef, h_neglect, answer_date)
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: h !< Cell widths [H]
  real, dimension(N),   intent(in)    :: u !< Cell averages in arbitrary coordinates [A]
  real, dimension(N,2), intent(inout) :: edge_values !< Edge values [A]
  real, dimension(N,3), intent(inout) :: ppoly_coef !< Polynomial coefficients, mainly [A]
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables
  integer   :: k              ! Loop index
  real      :: edge_l, edge_r ! Edge values (left and right) [A]

  ! PPM limiter
  call PPM_limiter_standard( N, h, u, edge_values, h_neglect, answer_date=answer_date )

  ! Loop over all cells
  do k = 1,N

    edge_l = edge_values(k,1)
    edge_r = edge_values(k,2)

    ! Store polynomial coefficients
    ppoly_coef(k,1) = edge_l
    ppoly_coef(k,2) = 4.0 * ( u(k) - edge_l ) + 2.0 * ( u(k) - edge_r )
    ppoly_coef(k,3) = 3.0 * ( ( edge_r - u(k) ) + ( edge_l - u(k) ) )

  enddo

end subroutine PPM_reconstruction

!> Adjusts edge values using the standard PPM limiter (Colella & Woodward, JCP 1984)
!! after first checking that the edge values are bounded by neighbors cell averages
!! and that the edge values are monotonic between cell averages.
subroutine PPM_limiter_standard( N, h, u, edge_values, h_neglect, answer_date )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Potentially modified edge values [A]
  real,                 intent(in)    :: h_neglect !< A negligibly small width [H]
  integer,    optional, intent(in)    :: answer_date  !< The vintage of the expressions to use

  ! Local variables
  integer   :: k              ! Loop index
  real      :: u_l, u_c, u_r  ! Cell averages (left, center and right) [A]
  real      :: edge_l, edge_r ! Edge values (left and right) [A]
  real      :: expr1, expr2   ! Temporary expressions [A2]

  ! Bound edge values
  call bound_edge_values( N, h, u, edge_values, h_neglect, answer_date=answer_date )

  ! Make discontinuous edge values monotonic
  call check_discontinuous_edge_values( N, u, edge_values )

  ! Loop on interior cells to apply the standard
  ! PPM limiter (Colella & Woodward, JCP 84)
  do k = 2,N-1

    ! Get cell averages
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)

    edge_l = edge_values(k,1)
    edge_r = edge_values(k,2)

    if ( (u_r - u_c)*(u_c - u_l) <= 0.0) then
      ! Flatten extremum
      edge_l = u_c
      edge_r = u_c
    else
      expr1 = 3.0 * (edge_r - edge_l) * ( (u_c - edge_l) + (u_c - edge_r))
      expr2 = (edge_r - edge_l) * (edge_r - edge_l)
      if ( expr1 > expr2 ) then
        ! Place extremum at right edge of cell by adjusting left edge value
        edge_l = u_c + 2.0 * ( u_c - edge_r )
        edge_l = max( min( edge_l, max(u_l, u_c) ), min(u_l, u_c) ) ! In case of round off
      elseif ( expr1 < -expr2 ) then
        ! Place extremum at left edge of cell by adjusting right edge value
        edge_r = u_c + 2.0 * ( u_c - edge_l )
        edge_r = max( min( edge_r, max(u_r, u_c) ), min(u_r, u_c) ) ! In case of round off
      endif
    endif
    ! This checks that the difference in edge values is representable
    ! and avoids overshoot problems due to round off.
    !### The 1.e-60 needs to have units of [A], so this dimensionally inconsistent.
    if ( abs( edge_r - edge_l )<max(1.e-60,epsilon(u_c)*abs(u_c)) ) then
      edge_l = u_c
      edge_r = u_c
    endif

    edge_values(k,1) = edge_l
    edge_values(k,2) = edge_r

  enddo ! end loop on interior cells

  ! PCM within boundary cells
  edge_values(1,:) = u(1)
  edge_values(N,:) = u(N)

end subroutine PPM_limiter_standard

!> Adjusts edge values using the original monotonicity constraint (Colella & Woodward, JCP 1984)
!! Based on hybgen_ppm_coefs
subroutine PPM_monotonicity( N, u, edge_values )
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: u !< cell average properties (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values !< Potentially modified edge values [A]

  ! Local variables
  integer   :: k      ! Loop index
  real      :: a6, da ! Normalized scalar curvature and slope [A]

  ! Loop on interior cells to impose monotonicity
  ! Eq. 1.10 of (Colella & Woodward, JCP 84)
  do k = 2,N-1
    if (((u(k+1)-u(k))*(u(k)-u(k-1)) <= 0.)) then !local extremum
      edge_values(k,1) = u(k)
      edge_values(k,2) = u(k)
    else
      da = edge_values(k,2)-edge_values(k,1)
      a6 = 6.0*u(k) - 3.0*(edge_values(k,1)+edge_values(k,2))
      if (da*a6 > da*da) then !peak in right half of zone
        edge_values(k,1) = 3.0*u(k) - 2.0*edge_values(k,2)
      elseif (da*a6 < -da*da) then !peak in left half of zone
        edge_values(k,2) = 3.0*u(k) - 2.0*edge_values(k,1)
      endif
    endif
  enddo ! end loop on interior cells

end subroutine PPM_monotonicity

!------------------------------------------------------------------------------
!> Reconstruction by parabolas within boundary cells
subroutine PPM_boundary_extrapolation( N, h, u, edge_values, ppoly_coef, h_neglect)
!------------------------------------------------------------------------------
! Reconstruction by parabolas within boundary cells.
!
! The following explanations apply to the left boundary cell. The same
! reasoning holds for the right boundary cell.
!
! A parabola needs to be built in the cell and requires three degrees of
! freedom, which are the right edge value and slope and the cell average.
! The right edge values and slopes are taken to be that of the neighboring
! cell (i.e., the left edge value and slope of the neighboring cell).
! The resulting parabola is not necessarily monotonic and the traditional
! PPM limiter is used to modify one of the edge values in order to yield
! a monotonic parabola.
!
! N:     number of cells in grid
! h:     thicknesses of grid cells
! u:     cell averages to use in constructing piecewise polynomials
! edge_values : edge values of piecewise polynomials
! ppoly_coef : coefficients of piecewise polynomials
!
! It is assumed that the size of the array 'u' is equal to the number of cells
! defining 'grid' and 'ppoly'. No consistency check is performed here.
!------------------------------------------------------------------------------

  ! Arguments
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(:),   intent(in)    :: h !< cell widths (size N) [H]
  real, dimension(:),   intent(in)    :: u !< cell averages (size N) [A]
  real, dimension(:,:), intent(inout) :: edge_values    !< edge values of piecewise polynomials [A]
  real, dimension(:,:), intent(inout) :: ppoly_coef !< coefficients of piecewise polynomials, mainly [A]
  real,                 intent(in)    :: h_neglect  !< A negligibly small width for
                                           !! the purpose of cell reconstructions [H]

  ! Local variables
  integer :: i0, i1
  real    :: u0, u1     ! Average concentrations in the two neighboring cells [A]
  real    :: h0, h1     ! Thicknesses of the two neighboring cells [H]
  real    :: a, b, c    ! An edge value, normalized slope and normalized curvature
                        ! of a reconstructed distribution [A]
  real    :: u0_l, u0_r ! Edge values of a neighboring cell [A]
  real    :: u1_l, u1_r ! Neighboring cell slopes renormalized by the thickness of
                        ! the cell being worked on [A]
  real    :: slope      ! The normalized slope [A]
  real    :: exp1, exp2 ! Temporary expressions [A2]

  ! ----- Left boundary -----
  i0 = 1
  i1 = 2
  h0 = h(i0)
  h1 = h(i1)
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the left edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly_coef(i1,2)
  u1_r = b *((h0+h_neglect)/(h1+h_neglect))     ! derivative evaluated at xi = 0.0,
                        ! expressed w.r.t. xi (local coord. system)

  ! Limit the right slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 )
  if ( abs(u1_r) > abs(slope) ) then
    u1_r = slope
  endif

  ! The right edge value in the boundary cell is taken to be the left
  ! edge value in the neighboring cell
  u0_r = edge_values(i1,1)

  ! Given the right edge value and slope, we determine the left
  ! edge value and slope by computing the parabola as determined by
  ! the right edge value and slope and the boundary cell average
  u0_l = 3.0 * u0 + 0.5 * u1_r - 2.0 * u0_r

  ! Apply the traditional PPM limiter
  exp1 = (u0_r - u0_l) * (u0 - 0.5*(u0_l+u0_r))
  exp2 = (u0_r - u0_l) * (u0_r - u0_l) / 6.0

  if ( exp1 > exp2 ) then
    u0_l = 3.0 * u0 - 2.0 * u0_r
  endif

  if ( exp1 < -exp2 ) then
    u0_r = 3.0 * u0 - 2.0 * u0_l
  endif

  edge_values(i0,1) = u0_l
  edge_values(i0,2) = u0_r

  a = u0_l
  b = 6.0 * u0 - 4.0 * u0_l - 2.0 * u0_r
  c = 3.0 * ( u0_r + u0_l - 2.0 * u0 )

  ppoly_coef(i0,1) = a
  ppoly_coef(i0,2) = b
  ppoly_coef(i0,3) = c

  ! ----- Right boundary -----
  i0 = N-1
  i1 = N
  h0 = h(i0)
  h1 = h(i1)
  u0 = u(i0)
  u1 = u(i1)

  ! Compute the right edge slope in neighboring cell and express it in
  ! the global coordinate system
  b = ppoly_coef(i0,2)
  c = ppoly_coef(i0,3)
  u1_l = (b + 2*c)                  ! derivative evaluated at xi = 1.0
  u1_l = u1_l * ((h1+h_neglect)/(h0+h_neglect))

  ! Limit the left slope by the PLM limited slope
  slope = 2.0 * ( u1 - u0 )
  if ( abs(u1_l) > abs(slope) ) then
    u1_l = slope
  endif

  ! The left edge value in the boundary cell is taken to be the right
  ! edge value in the neighboring cell
  u0_l = edge_values(i0,2)

  ! Given the left edge value and slope, we determine the right
  ! edge value and slope by computing the parabola as determined by
  ! the left edge value and slope and the boundary cell average
  u0_r = 3.0 * u1 - 0.5 * u1_l - 2.0 * u0_l

  ! Apply the traditional PPM limiter
  exp1 = (u0_r - u0_l) * (u1 - 0.5*(u0_l+u0_r))
  exp2 = (u0_r - u0_l) * (u0_r - u0_l) / 6.0

  if ( exp1 > exp2 ) then
    u0_l = 3.0 * u1 - 2.0 * u0_r
  endif

  if ( exp1 < -exp2 ) then
    u0_r = 3.0 * u1 - 2.0 * u0_l
  endif

  edge_values(i1,1) = u0_l
  edge_values(i1,2) = u0_r

  a = u0_l
  b = 6.0 * u1 - 4.0 * u0_l - 2.0 * u0_r
  c = 3.0 * ( u0_r + u0_l - 2.0 * u1 )

  ppoly_coef(i1,1) = a
  ppoly_coef(i1,2) = b
  ppoly_coef(i1,3) = c

end subroutine PPM_boundary_extrapolation

end module PPM_functions
