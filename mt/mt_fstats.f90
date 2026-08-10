! ---------------------------------------------------------------------------
! The array exports below take n-sized explicit-shape dummies: every trail
! array is described by the sample count `n` that travels with it, so the
! caller passes its own buffer by first element and nothing here is bounded by
! a constant the caller has to agree with. `n` is therefore declared ahead of
! the arrays it sizes -- a specification expression may only use entities
! already declared. There is no upper bound on n at this boundary: the arrays
! are exactly as long as the caller says, and any resource cap (the quadratic
! Theil-Sen estimator excepted, which carries its own in mt_robust.f90) is the
! caller's policy, not an ABI constraint.
! ---------------------------------------------------------------------------
module MT

implicit none

contains

subroutine ComputeLowessFit(x, y, ys, n, fsmooth, nstps) bind(C, name="ComputeLowessFit")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use fstats, only : lowess

    ! The number of points in the input arrays.
    integer(c_int), intent(in) :: n
    ! An N-element array containing the independent variable data.  This
    ! array must be monotonically increasing.
    real(c_double), intent(in), dimension(n) :: x
    ! An N-element array containing the dependent variable data.
    real(c_double), intent(in), dimension(n) :: y
    ! An N-element array where the smoothed results will be written.
    real(c_double), intent(out), dimension(n) :: ys
    ! An optional input that specifies the amount of smoothing.
    ! Specifically, this value is the fraction of points used to compute
    ! each value.  As this value increases, the output becomes smoother.
    ! Choosing a value in the range of 0.2 to 0.8 typically results in a
    ! good fit.  The default value is 0.2.
    real(c_double), intent(in) :: fsmooth
    ! An optional input that specifies the number of iterations.  If set to
    ! zero, a non-robust fit is returned.  The default value is set to 2.
    integer(c_int), intent(in) :: nstps

    call lowess(x, y, ys, fsmooth, nstps)

end subroutine ComputeLowessFit

subroutine ComputeExpFitLBFGSB(&
    x, y, n, seed, lower, upper, maxit, lmm, factr, pgtol, coef, fval, convergence, status &
) bind(C, name="ComputeExpFitLBFGSB")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use MTExpFit, only : ExpFitLBFGSB

    ! The number of samples in x and y.
    integer(c_int), intent(in) :: n
    ! An N-element array of the independent variable (e.g. time).
    real(c_double), intent(in), dimension(n) :: x
    ! An N-element array of the dependent variable (e.g. distance).
    real(c_double), intent(in), dimension(n) :: y
    ! Initial guess for (a0, aa, bb, cc).
    real(c_double), intent(in), dimension(1:4) :: seed
    ! Lower box bounds for (a0, aa, bb, cc).
    real(c_double), intent(in), dimension(1:4) :: lower
    ! Upper box bounds for (a0, aa, bb, cc).
    real(c_double), intent(in), dimension(1:4) :: upper
    ! Maximum number of L-BFGS-B iterations.
    integer(c_int), intent(in) :: maxit
    ! Number of BFGS corrections retained (the "m" of L-BFGS-B).
    integer(c_int), intent(in) :: lmm
    ! Relative-reduction convergence control (factr*epsmch).
    real(c_double), intent(in) :: factr
    ! Projected-gradient convergence tolerance.
    real(c_double), intent(in) :: pgtol
    ! Fitted (a0, aa, bb, cc) on return.
    real(c_double), intent(out), dimension(1:4) :: coef
    ! Sum of squared residuals at the fitted parameters.
    real(c_double), intent(out) :: fval
    ! Convergence code: 0 converged, 1 iter limit, 51 warning,
    ! 52 abnormal/other.
    integer(c_int), intent(out) :: convergence
    ! 0 if the optimizer ran; 1 if it could not be started (invalid
    ! sizes / allocation failure).
    integer(c_int), intent(out) :: status

    coef = seed
    fval = 0.0d0
    convergence = 52
    status = 0

    ! An empty series has nothing to fit. ExpFitLBFGSB rejects it too; rejecting
    ! it here as well keeps every export in this file answering an unusable n the
    ! same way, whatever the routine behind it does.
    if (n < 1) then
        status = 1
        return
    end if

    call ExpFitLBFGSB(&
        n, x, y, seed, lower, upper, maxit, lmm, factr, pgtol, coef, fval, &
        convergence, status &
    )

end subroutine ComputeExpFitLBFGSB

subroutine ComputeExpFitGN(&
    x, y, n, seed, maxit, tol, minFactor, coef, se, sigma, convergence, status &
) bind(C, name="ComputeExpFitGN")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use MTExpFit, only : ExpFitGaussNewton

    ! The number of samples in x and y.
    integer(c_int), intent(in) :: n
    ! An N-element array of the independent variable (e.g. time).
    real(c_double), intent(in), dimension(n) :: x
    ! An N-element array of the dependent variable (e.g. distance).
    real(c_double), intent(in), dimension(n) :: y
    ! Initial guess for (a0, aa, bb, cc).
    real(c_double), intent(in), dimension(1:4) :: seed
    ! Maximum number of Gauss-Newton iterations.
    integer(c_int), intent(in) :: maxit
    ! Relative-offset convergence tolerance.
    real(c_double), intent(in) :: tol
    ! Minimum step-halving factor before declaring non-convergence.
    real(c_double), intent(in) :: minFactor
    ! Fitted (a0, aa, bb, cc) on return (the last iterate on non-convergence).
    real(c_double), intent(out), dimension(1:4) :: coef
    ! Per-coefficient standard errors sqrt(sigma^2 * (J^T J)^-1_ii); zero on
    ! non-convergence.
    real(c_double), intent(out), dimension(1:4) :: se
    ! Residual standard error sqrt(RSS / (n - 4)).
    real(c_double), intent(out) :: sigma
    ! Convergence code: 0 converged, 1 iteration limit, 2 step factor below
    ! minFactor, 3 singular Jacobian / linear-algebra failure.
    integer(c_int), intent(out) :: convergence
    ! 0 if the fit ran; 1 if it could not be started (n < 5, maxit < 1, or
    ! non-positive tol/minFactor).
    integer(c_int), intent(out) :: status

    coef = 0.0d0
    se = 0.0d0
    sigma = 0.0d0
    convergence = 1
    status = 0

    if (n < 1) then
        status = 1
        return
    end if

    call ExpFitGaussNewton(&
        n, x, y, seed, maxit, tol, minFactor, coef, se, sigma, convergence, status &
    )

end subroutine ComputeExpFitGN

! -----------------------------------------------------------------------
! ComputeExpFitPort -- Bounded exp-deceleration model fit, the native
! bounded exp-deceleration fit: the same analytic-Jacobian Gauss-Newton as
! ComputeExpFitGN, with every trial step clamped into the box
! [lower, upper]. Stateless.
! -----------------------------------------------------------------------
subroutine ComputeExpFitPort(&
    x, y, n, seed, lower, upper, maxit, tol, coef, se, sigma, convergence, status &
) bind(C, name="ComputeExpFitPort")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use iso_fortran_env, only : real64
    use MTExpFit, only : ExpFitGaussNewtonBounded

    ! The number of samples in x and y.
    integer(c_int), intent(in) :: n
    ! An N-element array of the independent variable (e.g. time).
    real(c_double), intent(in), dimension(n) :: x
    ! An N-element array of the dependent variable (e.g. distance).
    real(c_double), intent(in), dimension(n) :: y
    ! Initial guess for (a0, aa, bb, cc); projected onto the box if outside it.
    real(c_double), intent(in), dimension(1:4) :: seed
    ! Lower box bounds for (a0, aa, bb, cc).
    real(c_double), intent(in), dimension(1:4) :: lower
    ! Upper box bounds for (a0, aa, bb, cc).
    real(c_double), intent(in), dimension(1:4) :: upper
    ! Maximum number of Gauss-Newton iterations.
    integer(c_int), intent(in) :: maxit
    ! Relative-offset convergence tolerance.
    real(c_double), intent(in) :: tol
    ! Fitted (a0, aa, bb, cc) on return (the last iterate on non-convergence).
    real(c_double), intent(out), dimension(1:4) :: coef
    ! Per-coefficient standard errors sqrt(sigma^2 * (J^T J)^-1_ii); zero on
    ! non-convergence.
    real(c_double), intent(out), dimension(1:4) :: se
    ! Residual standard error sqrt(RSS / (n - 4)).
    real(c_double), intent(out) :: sigma
    ! Convergence code: 0 converged, 1 iteration limit, 2 step factor below
    ! minFactor, 3 singular Jacobian / linear-algebra failure.
    integer(c_int), intent(out) :: convergence
    ! 0 if the fit ran; 1 if it could not be started (n < 5, maxit < 1,
    ! non-positive tol, or an empty box lower > upper).
    integer(c_int), intent(out) :: status

    ! Baked in because this wrapper does not expose minFactor as a knob.
    real(real64), parameter :: MINFACTOR = 1.0d0 / 1024.0d0

    coef = 0.0d0
    se = 0.0d0
    sigma = 0.0d0
    convergence = 1
    status = 0

    if (n < 1) then
        status = 1
        return
    end if

    call ExpFitGaussNewtonBounded(&
        n, x, y, seed, lower, upper, maxit, tol, MINFACTOR, coef, se, sigma, &
        convergence, status &
    )

end subroutine ComputeExpFitPort

subroutine ComputeTheilSenFit(&
    x, y, n, slope, intercept, sigma, seSlope, seIntercept, status &
) bind(C, name="ComputeTheilSenFit")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use MTRobust, only : TheilSenFit, MAX_THEILSEN_LENGTH

    ! The number of samples in x and y.
    integer(c_int), intent(in) :: n
    ! An N-element array of the independent variable (e.g. time).
    real(c_double), intent(in), dimension(n) :: x
    ! An N-element array of the dependent variable (e.g. distance).
    real(c_double), intent(in), dimension(n) :: y
    ! Theil-Sen slope estimate.
    real(c_double), intent(out) :: slope
    ! Intercept: median(y - slope*x).
    real(c_double), intent(out) :: intercept
    ! Residual standard error sqrt(sum(resid^2)/(n-2)).
    real(c_double), intent(out) :: sigma
    ! MAD of the pairwise slopes (the coefficient error).
    real(c_double), intent(out) :: seSlope
    ! MAD of the per-point intercepts (the coefficient error).
    real(c_double), intent(out) :: seIntercept
    ! 0 if the fit was computed; 1 if it could not be started (n < 3,
    ! n > MAX_THEILSEN_LENGTH, or a working-buffer allocation failure); 2 if the
    ! data are degenerate (all x equal, so no slope is defined).
    integer(c_int), intent(out) :: status

    slope = 0.0d0
    intercept = 0.0d0
    sigma = 0.0d0
    seSlope = 0.0d0
    seIntercept = 0.0d0
    status = 0

    ! n >= 3 keeps the sigma denominator (n-2) positive. The upper bound is
    ! MAX_THEILSEN_LENGTH: this estimator is defined over all n(n-1)/2 pairwise
    ! slopes and has to materialise them, so its cost is quadratic in memory and
    ! time. It is the one size cap this file still imposes -- it is the
    ! estimator's own, not a marshalling limit. See mt_robust.f90.
    if (n < 3 .or. n > MAX_THEILSEN_LENGTH) then
        status = 1
        return
    end if

    ! With every x identical there is no distinct-x pair and no slope; report
    ! it rather than divide by zero.
    if (all(x == x(1))) then
        status = 2
        return
    end if

    call TheilSenFit(n, x, y, slope, intercept, sigma, seSlope, seIntercept, status)

end subroutine ComputeTheilSenFit

subroutine ComputeSiegelFit(&
    x, y, n, slope, intercept, sigma, seSlope, seIntercept, status &
) bind(C, name="ComputeSiegelFit")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use MTRobust, only : SiegelFit

    ! The number of samples in x and y.
    integer(c_int), intent(in) :: n
    ! An N-element array of the independent variable (e.g. time).
    real(c_double), intent(in), dimension(n) :: x
    ! An N-element array of the dependent variable (e.g. distance).
    real(c_double), intent(in), dimension(n) :: y
    ! Siegel repeated-median slope estimate.
    real(c_double), intent(out) :: slope
    ! Siegel repeated-median intercept estimate.
    real(c_double), intent(out) :: intercept
    ! Residual standard error sqrt(sum(resid^2)/(n-2)).
    real(c_double), intent(out) :: sigma
    ! MAD of the per-point median slopes (the coefficient error).
    real(c_double), intent(out) :: seSlope
    ! MAD of the per-point median intercepts (the coefficient error).
    real(c_double), intent(out) :: seIntercept
    ! 0 if the fit was computed; 1 if it could not be started (n < 3 or a
    ! working-buffer allocation failure); 2 if the data are degenerate (all x
    ! equal, so no slope is defined).
    integer(c_int), intent(out) :: status

    slope = 0.0d0
    intercept = 0.0d0
    sigma = 0.0d0
    seSlope = 0.0d0
    seIntercept = 0.0d0
    status = 0

    if (n < 3) then
        status = 1
        return
    end if

    if (all(x == x(1))) then
        status = 2
        return
    end if

    call SiegelFit(n, x, y, slope, intercept, sigma, seSlope, seIntercept, status)

end subroutine ComputeSiegelFit

! -----------------------------------------------------------------------
! ComputeOlsFit -- Ordinary least-squares polynomial fit with
! per-coefficient standard
! errors. The design matrix is
!   nCoef = 2 : [1, x]        (a plain linear fit)
!   nCoef = 3 : [1, x, x^2/2] (the half-squared-times convention)
! so the third column is x^2/2 -- NOT x^2. That factor of one half changes
! the third coefficient and its standard error, so this routine builds the
! design matrix explicitly rather than calling fstats' design_matrix (which
! only ever produces the plain-power columns [1, x, x^2, ...]). The
! coefficients, residuals and standard errors are then formed exactly as
! fstats' linear_least_squares / calculate_regression_statistics do:
! coeffs = (X^T X)^-1 X^T y (via fstats' SVD-based covariance_matrix),
! sigma = sqrt(SSR / (n - nCoef)), se_i = sqrt(sigma^2 * (X^T X)^-1_ii).
! Stateless: all working storage is local to this call.
!
! RANK: covariance_matrix is an SVD PSEUDO-inverse, so it does not fail on a
! rank-deficient design -- it quietly returns the minimum-norm solution,
! whose fitted values are right but whose split between coefficients is
! arbitrary and whose standard errors are meaningless. This export returns
! a fixed-width coefficient vector and has no way to signal "column
! dropped", so instead the design is rank-tested up front from its
! singular values and a deficient one is rejected with status = 2.
! -----------------------------------------------------------------------
subroutine ComputeOlsFit(x, y, n, nCoef, coef, se, sigma, status) bind(C, name="ComputeOlsFit")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use iso_fortran_env, only : real64, int32
    use fstats_regression, only : covariance_matrix
    use ferror, only : errors

    ! The number of samples in x and y.
    integer(c_int), intent(in) :: n
    ! An N-element array of the independent variable (e.g. time).
    real(c_double), intent(in), dimension(n) :: x
    ! An N-element array of the dependent variable (e.g. distance).
    real(c_double), intent(in), dimension(n) :: y
    ! The number of coefficients: 2 for [1, x], 3 for [1, x, x^2/2].
    integer(c_int), intent(in) :: nCoef
    ! The fitted coefficients (only the leading nCoef are meaningful).
    real(c_double), intent(out), dimension(1:3) :: coef
    ! The per-coefficient standard errors (only the leading nCoef meaningful).
    real(c_double), intent(out), dimension(1:3) :: se
    ! The residual standard error sqrt(SSR / (n - nCoef)).
    real(c_double), intent(out) :: sigma
    ! 0 if the fit was computed; 1 if it could not be started (nCoef not in
    ! {2, 3}, or n <= nCoef so the residual degrees of freedom would be <= 0);
    ! 2 if the design matrix is numerically rank-deficient (see RANK_TOL below),
    ! or the covariance solve failed (an allocation failure inside
    ! covariance_matrix).
    integer(c_int), intent(out) :: status

    ! LAPACK SVD, used only for the singular values (resolved via the linked
    ! liblapack).
    external :: dgesvd

    real(real64), parameter :: zero = 0.0d0
    real(real64), parameter :: half = 0.5d0
    ! Relative singular-value floor for the rank test: a design is rank-deficient
    ! when s(nCoef) <= RANK_TOL * s(1) -- see the header.
    real(real64), parameter :: RANK_TOL = 1.0d-7

    real(real64), allocatable :: a(:,:), cov(:,:), xty(:), coeffs(:)
    real(real64), allocatable :: ymod(:), resid(:)
    real(real64), allocatable :: asv(:,:), sval(:), svwork(:)
    real(real64) :: ssr, var, svq(1), udum(1,1), vdum(1,1)
    integer(int32) :: i, dof, svinfo, svlwork
    ! A private error handler with termination disabled, so a failure inside
    ! covariance_matrix is reported back as status 2 rather than aborting the
    ! whole host process (ferror's default is exit-on-error).
    type(errors) :: errmgr

    coef = zero
    se = zero
    sigma = zero
    status = 0

    ! nCoef selects the design; only the linear and quadratic forms are supported.
    if (nCoef < 2 .or. nCoef > 3) then
        status = 1
        return
    end if
    ! n > nCoef keeps the residual degrees of freedom (n - nCoef) positive.
    if (n <= nCoef) then
        status = 1
        return
    end if

    allocate(&
        a(n, nCoef), cov(nCoef, nCoef), xty(nCoef), coeffs(nCoef), ymod(n), resid(n) &
    )

    ! Design matrix columns: [1, x] and, for nCoef = 3, x^2/2 (see the header).
    a(:, 1) = 1.0d0
    a(:, 2) = x
    if (nCoef == 3) a(:, 3) = half * x * x

    ! Numerical rank of the design from its singular values. dgesvd destroys its
    ! input and needs no singular vectors here ('N', 'N'), so it runs on a scratch
    ! copy. A failed SVD is itself a reason not to trust the fit -> status 2.
    allocate(asv(n, nCoef), sval(min(n, nCoef)))
    asv = a
    call dgesvd('N', 'N', n, nCoef, asv, n, sval, udum, 1, vdum, 1, svq, -1, svinfo)
    svlwork = max(int(svq(1)), 1)
    allocate(svwork(svlwork))
    call dgesvd('N', 'N', n, nCoef, asv, n, sval, udum, 1, vdum, 1, svwork, svlwork, svinfo)
    if (svinfo /= 0 .or. sval(1) <= zero .or. &
            sval(nCoef) <= RANK_TOL * sval(1)) then
        status = 2
        deallocate(a, cov, xty, coeffs, ymod, resid)
        deallocate(asv, sval, svwork)
        return
    end if
    deallocate(asv, sval, svwork)

    ! cov = (X^T X)^-1 via fstats' covariance_matrix (DGEMM + SVD pseudo-inverse).
    call errmgr%set_exit_on_error(.false.)
    call covariance_matrix(a, cov, errmgr)
    if (errmgr%has_error_occurred()) then
        status = 2
        deallocate(a, cov, xty, coeffs, ymod, resid)
        return
    end if

    ! coeffs = cov * (X^T y); model values and residuals follow directly.
    xty = matmul(transpose(a), y)
    coeffs = matmul(cov, xty)
    ymod = matmul(a, coeffs)
    resid = ymod - y

    ! Residual standard error and coefficient standard errors, matching
    ! var = SSR/(n - nCoef), se_i = sqrt(var * cov_ii).
    dof = n - nCoef
    ssr = sum(resid * resid)
    var = ssr / real(dof, real64)
    sigma = sqrt(var)
    do i = 1, nCoef
        coef(i) = coeffs(i)
        se(i) = sqrt(var * cov(i, i))
    end do

    deallocate(a, cov, xty, coeffs, ymod, resid)

end subroutine ComputeOlsFit

! -----------------------------------------------------------------------
! ComputeFTestPValue -- Nested-model F-test p-value, replacing
! Given the residual sums of squares and residual
! degrees of freedom of a null (restricted) and an alternative (full)
! model, returns the upper-tail F probability
!   p = pf(F, dfNull-dfAlt, dfAlt, lower.tail=FALSE)
! with F = ((rssNull-rssAlt)/(dfNull-dfAlt)) / (rssAlt/dfAlt), computed as a
! DIRECT upper tail from the regularized incomplete beta (see below), not as
! 1 - CDF. Edge cases:
!   * dfAlt <= 0 or dfNull <= dfAlt -> status 1 (caller returns NA/NaN):
!     the F statistic would have non-positive degrees of freedom.
!   * non-finite rssAlt -> status 1: a rank-deficient alternative fit
!     cannot be assessed (NA coefficients yield an NA RSS).
!   * rssAlt <= 0 with rssNull > 0 -> p = 0: a perfect alternative fit is
!     infinitely significant. When the null fit is ALSO perfect (or rssNull
!     is NaN), F is 0/0 = NaN, so that corner is undefined (status 1):
!     two exact fits carry no evidence for the extra term,
!     and p = 0 would report infinite significance on none.
!   * F <= 0 -> p = 1: the alternative model did not reduce the residual sum
!     of squares (possible from roundoff when the extra term is worthless).
!     The beta argument below would be out of range, so this is handled
!     up front.
! rssNull is deliberately NOT guarded for finiteness: an NaN rssNull gives
! an NaN F and hence NaN, and an infinite rssNull gives F = +Inf and hence
! p = 0.
! Stateless.
! -----------------------------------------------------------------------
subroutine ComputeFTestPValue(&
    rssNull, dfNull, rssAlt, dfAlt, p, status &
) bind(C, name="ComputeFTestPValue")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use iso_fortran_env, only : real64
    use fstats_special_functions, only : regularized_beta

    ! Residual sum of squares of the null (restricted) model.
    real(c_double), intent(in) :: rssNull
    ! Residual degrees of freedom of the null model.
    integer(c_int), intent(in) :: dfNull
    ! Residual sum of squares of the alternative (full) model.
    real(c_double), intent(in) :: rssAlt
    ! Residual degrees of freedom of the alternative model.
    integer(c_int), intent(in) :: dfAlt
    ! The upper-tail p-value (meaningful only when status = 0).
    real(c_double), intent(out) :: p
    ! 0 if p was computed (including the perfect-alternative p = 0 branch);
    ! 1 for the undefined cases (dfAlt <= 0, dfNull <= dfAlt, non-finite rssAlt, or
    ! a perfect alternative fit without a worse null fit to compare against).
    integer(c_int), intent(out) :: status

    real(real64) :: fstat, ddf, adf, z

    p = 0.0d0
    status = 0

    ! Non-positive F degrees of freedom -> NA.
    if (dfAlt <= 0 .or. dfNull <= dfAlt) then
        status = 1
        return
    end if

    ! Non-finite alternative RSS -> NA.
    if (.not. ieee_is_finite(rssAlt)) then
        status = 1
        return
    end if

    ! Perfect (or, out of contract, negative-RSS) alternative fit. R forms
    ! F = ((rssNull - rssAlt)/ddf) / (rssAlt/adf): with a positive rssNull that
    ! is +Inf and p is 0 -- infinitely significant. With rssNull 0 (both
    ! fits exact) it is 0/0 = NaN -> undefined, status 1. The
    ! test is written NaN-rejecting (.not. >) so an NaN rssNull lands on the
    ! same NA path here that it reaches through fstat in the normal branch.
    if (rssAlt <= 0.0d0) then
        if (.not. (rssNull > 0.0d0)) then
            status = 1
        else
            p = 0.0d0
        end if
        return
    end if

    ddf = real(dfNull - dfAlt, real64)
    adf = real(dfAlt, real64)
    fstat = ((rssNull - rssAlt) / ddf) / (rssAlt / adf)

    ! A non-positive F carries no evidence for the alternative: p = 1. Guarding
    ! here also keeps the beta argument below inside [0, 1] -- 1 - CDF instead
    ! returned values ABOVE 1 (1.67 at F = -0.5), an out-of-range "probability"
    ! nothing downstream could detect.
    if (fstat <= 0.0d0) then
        p = 1.0d0
        return
    end if

    ! The lower tail is pf(x, d1, d2) = I_z(d1/2, d2/2) with z = d1*x/(d1*x+d2)
    ! (fstats' f_distribution%cdf). Taking 1 - that loses ALL relative precision
    ! on a significant test, because the CDF is itself formed as 1 - small_tail:
    ! p collapsed to exactly 0 from about 1e-16 down, so an ordinary F = 9490 on
    ! (2, 19) d.o.f. reported 0 where R reports 3.16e-29. The beta symmetry
    ! I_z(a,b) = 1 - I_{1-z}(b,a) gives the upper tail DIRECTLY, evaluated at
    ! 1-z = d2/(d1*x+d2) with the arguments swapped -- no cancellation.
    z = adf / (ddf * fstat + adf)
    p = regularized_beta(0.5d0 * adf, 0.5d0 * ddf, z)

    ! The continued fraction is a fixed 20-term truncation, so pin the result to
    ! the valid range. NaN (a non-finite rssNull) fails both tests and survives
    ! as NaN, which is the documented R-NA behaviour.
    if (p < 0.0d0) p = 0.0d0
    if (p > 1.0d0) p = 1.0d0

end subroutine ComputeFTestPValue

! -----------------------------------------------------------------------
! ComputeLoessFit -- Native loess (locally-weighted regression) matching
! the plot overlays and the light-mass/speed loess columns.
!
! Each fitted value comes from an exact local polynomial fit at its own
! abscissa rather than from interpolating a surface -- see mt_loess.f90.
! Stateless: all working storage is local to this call.
!
! For each input abscissa x_i, take q = min(n, floor(span*n + 1e-5)) nearest
! neighbours by |x - x_i|, weight them with the tricube kernel
! (1 - (d/dmax)^3)^3 where dmax is the distance to the q-th neighbour, and
! fit a weighted polynomial of the given degree; ys(i) is its value at x_i.
! -----------------------------------------------------------------------
subroutine ComputeLoessFit(x, y, n, span, degree, ys, status) bind(C, name="ComputeLoessFit")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use MTLoess, only : LoessFit

    ! The number of samples in x, y and ys.
    integer(c_int), intent(in) :: n
    ! An N-element array containing the independent variable data.
    real(c_double), intent(in), dimension(n) :: x
    ! An N-element array containing the dependent variable data.
    real(c_double), intent(in), dimension(n) :: y
    ! Neighbourhood fraction; must be > 0.
    real(c_double), intent(in) :: span
    ! Local polynomial degree: 1 (locally linear) or 2 (locally quadratic).
    integer(c_int), intent(in) :: degree
    ! An N-element array where the fitted values are written. Meaningful only
    ! when status = 0; on either rejection below it is left untouched, and on a
    ! degenerate fit it holds whatever LoessFit reached before giving up.
    real(c_double), intent(out), dimension(n) :: ys
    ! 0 if the fit was computed; 1 if it could not be started (n < 1, degree not
    ! in {1, 2}, or span <= 0); 2 if a local neighbourhood was degenerate (too
    ! few positive-weight points to determine the local polynomial, coincident
    ! abscissae, or a singular normal matrix).
    integer(c_int), intent(out) :: status

    status = 0

    ! Only the linear and quadratic surfaces are supported; degree selects the
    ! local polynomial. n > 0 and a positive span are required for a meaningful
    ! neighbourhood.
    if (n < 1) then
        status = 1
        return
    end if
    if (degree < 1 .or. degree > 2) then
        status = 1
        return
    end if
    if (span <= 0.0d0) then
        status = 1
        return
    end if

    call LoessFit(n, x, y, span, degree, ys, status)

end subroutine ComputeLoessFit

! -----------------------------------------------------------------------
! ComputeBorovickaFit -- Native per-sample Borovicka straight-line
! radiant solver. Minimises the weighted sum of squared
! line-to-line distances between a straight trajectory and every sight, over
! a 4-parameter pivot chart, using the vendored L-BFGS-B 3.0 optimizer with
! an analytic gradient. See mt_borovicka.f90 for the full algorithm.
!
! Non-convergence (a failed fit, a still-bound-active re-pivot, or a
! degenerate radiant) is signalled by NaN outputs with converged = 0 -- the
! DOCUMENTED signal, not an error. All quantities are in Mm (megametres);
! the Pascal wrapper scales the point to metres.
! -----------------------------------------------------------------------
subroutine ComputeBorovickaFit(&
    stations, nStations, stIdx, xi, eta, zeta, w, nPoints, seedPoint, seedRadiant, maxit, lmm, &
    factr, pgtol, outRadiant, outPoint, converged, status &
) bind(C, name="ComputeBorovickaFit")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use MTBorovicka, only : BorovickaMinimize4

    ! Station coordinates, column-major 3 x nStations (Mm): stations(:, j) is
    ! the (x, y, z) of station j.
    integer(c_int), intent(in) :: nStations
    real(c_double), intent(in) :: stations(3, nStations)
    ! The number of sights (rows of the observation set).
    integer(c_int), intent(in) :: nPoints
    ! 1-based station index of each sight.
    integer(c_int), intent(in) :: stIdx(nPoints)
    ! Per-sight line-of-sight unit-direction components and weight.
    real(c_double), intent(in) :: xi(nPoints), eta(nPoints), zeta(nPoints)
    real(c_double), intent(in) :: w(nPoints)
    ! Seed trajectory point (Mm) and seed radiant (unitless direction).
    real(c_double), intent(in) :: seedPoint(3), seedRadiant(3)
    ! L-BFGS-B controls. Pass maxit = 300 -- see mt_borovicka.f90.
    integer(c_int), intent(in) :: maxit, lmm
    real(c_double), intent(in) :: factr, pgtol
    ! Unit radiant and trajectory point (Mm) on return; NaN when converged = 0.
    real(c_double), intent(out) :: outRadiant(3), outPoint(3)
    ! 1 if the fit converged, 0 otherwise (NaN outputs).
    integer(c_int), intent(out) :: converged
    ! 0 if the solver ran; 1 if it could not be started (bad sizes, an
    ! out-of-range station index, or a workspace allocation failure).
    integer(c_int), intent(out) :: status

    call BorovickaMinimize4(&
        nStations, stations, nPoints, stIdx, xi, eta, zeta, w, seedPoint, seedRadiant, maxit, &
        lmm, factr, pgtol, outRadiant, outPoint, converged, status &
    )

end subroutine ComputeBorovickaFit

! -----------------------------------------------------------------------
! ComputeBorovickaKernel -- Signed line-to-line distance kernel.
! TEST-SUPPORT export: it exposes the internal kernel so the Pascal suite
! can assert the gauge invariances (R -> kR, P -> P + sR) and the finite
! parallel-sentinel path. Stateless, pure.
! -----------------------------------------------------------------------
subroutine ComputeBorovickaKernel(P, R, S, u, dist) bind(C, name="ComputeBorovickaKernel")
    use, intrinsic :: iso_c_binding, only: c_double
    use MTBorovicka, only : BorovickaLineDistance

    ! Trajectory point P and direction R, station S and sight unit direction u.
    real(c_double), intent(in) :: P(3), R(3), S(3), u(3)
    ! The signed distance (a large finite sentinel when the lines are parallel).
    real(c_double), intent(out) :: dist

    dist = BorovickaLineDistance(P, R, S, u)

end subroutine ComputeBorovickaKernel

! -----------------------------------------------------------------------
! ComputeBorovickaCostGrad -- FF4 cost and its analytic gradient at a
! given 4-parameter iterate q under the pivot chart. TEST-SUPPORT export:
! it lets the Pascal suite check the analytic gradient against central
! differences and the near-parallel finiteness of both cost and gradient.
! Stateless.
! -----------------------------------------------------------------------
subroutine ComputeBorovickaCostGrad(&
    stations, nStations, stIdx, xi, eta, zeta, w, nPoints, q, k, pointFixed, cost, grad &
) bind(C, name="ComputeBorovickaCostGrad")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use MTBorovicka, only : BorovickaFF4Cost, BorovickaFF4Grad

    ! Station coordinates, column-major 3 x nStations (Mm).
    integer(c_int), intent(in) :: nStations
    real(c_double), intent(in) :: stations(3, nStations)
    ! Sights: count, 1-based station index, LOS components and weights.
    integer(c_int), intent(in) :: nPoints
    integer(c_int), intent(in) :: stIdx(nPoints)
    real(c_double), intent(in) :: xi(nPoints), eta(nPoints), zeta(nPoints)
    real(c_double), intent(in) :: w(nPoints)
    ! The 4 free parameters, the pivot index k (1-based) and the pinned
    ! coordinate pointFixed defining the chart.
    real(c_double), intent(in) :: q(4)
    integer(c_int), intent(in) :: k
    real(c_double), intent(in) :: pointFixed
    ! The cost and its 4-element gradient at q. Both are NaN if k is not a valid
    ! pivot index (1, 2 or 3) or a station index is out of range -- this export
    ! has no status argument, so NaN is how it reports unusable input.
    real(c_double), intent(out) :: cost
    real(c_double), intent(out) :: grad(4)

    ! Neither of these is a theoretical concern only. k indexes 3-element arrays
    ! inside Reconstruct (P(k) = pointFixed, R(k) = 1) with no internal check, and
    ! FreeIndices maps any k outside {1, 2, 3} to the k = 3 case rather than
    ! failing, so a bad k is an out-of-bounds WRITE; an out-of-range stIdx is an
    ! out-of-bounds read of stations(:, st). ComputeBorovickaFit validates both
    ! before it charts, and this export -- test support, but an exported C symbol
    ! all the same -- must not be the softer way in.
    if (k < 1 .or. k > 3 .or. nStations < 1 .or. nPoints < 1) then
        cost = ieee_value(0.0_c_double, ieee_quiet_nan)
        grad = ieee_value(0.0_c_double, ieee_quiet_nan)
        return
    end if
    if (any(stIdx < 1) .or. any(stIdx > nStations)) then
        cost = ieee_value(0.0_c_double, ieee_quiet_nan)
        grad = ieee_value(0.0_c_double, ieee_quiet_nan)
        return
    end if

    call BorovickaFF4Cost(&
        nStations, stations, nPoints, stIdx, xi, eta, zeta, w, q, k, pointFixed, cost &
    )
    call BorovickaFF4Grad(&
        nStations, stations, nPoints, stIdx, xi, eta, zeta, w, q, k, pointFixed, grad &
    )

end subroutine ComputeBorovickaCostGrad

end module MT
