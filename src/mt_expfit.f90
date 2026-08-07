! =======================================================================
! MTExpFit -- Exponential-deceleration model fit.
!
! The model is  s(x) = a0 + aa*x - bb*exp(cc*x)  with parameter vector
! p = (a0, aa, bb, cc), fitted by minimising the sum of squared residuals
!   f(p) = sum( (y - (p1 + p2*x - p3*exp(p4*x)))^2 )
!
! Three fits are offered, all with an analytic Jacobian or gradient:
!   ExpFitLBFGSB              box-bounded, via the vendored L-BFGS-B 3.0
!   ExpFitGaussNewton         unbounded Gauss-Newton
!   ExpFitGaussNewtonBounded  Gauss-Newton with each step projected onto
!                             the box
!
! Stateless: all working storage is local to each call.
! =======================================================================
module MTExpFit

use, intrinsic :: iso_fortran_env, only : real64, int32
implicit none
private

public :: ExpModelCost, ExpModelGrad, ExpFitLBFGSB
public :: ExpFitGaussNewton, ExpFitGaussNewtonBounded

! Number of model parameters (a0, aa, bb, cc). The optimization runs in
! this dimension; it is distinct from the number of data samples.
integer(int32), parameter :: NPARAM = 4

contains

! -----------------------------------------------------------------------
! ExpModelValues -- Fitted values of the exp-deceleration model s(x) = a0
! + aa*x - bb*exp(cc*x) with p = (a0, aa, bb, cc). The single site for
! the model formula; the cost, residual and Jacobian below all derive
! from it.
! -----------------------------------------------------------------------
pure subroutine ExpModelValues(n, x, p, f)
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: x(n)
    real(real64), intent(in) :: p(NPARAM)
    real(real64), intent(out) :: f(n)

    f = p(1) + p(2)*x - p(3)*exp(p(4)*x)
end subroutine ExpModelValues

! -----------------------------------------------------------------------
! ExpModelCost -- Sum of squared residuals of the exp-deceleration model.
! -----------------------------------------------------------------------
pure subroutine ExpModelCost(n, x, y, p, f)
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: x(n), y(n)
    real(real64), intent(in) :: p(NPARAM)
    real(real64), intent(out) :: f

    real(real64) :: fvals(n), e(n)

    call ExpModelValues(n, x, p, fvals)
    e = y - fvals
    f = sum(e*e)
end subroutine ExpModelCost

! -----------------------------------------------------------------------
! ExpModelJacobian -- Analytic Jacobian of the model FITTED VALUES w.r.t.
! (a0, aa, bb, cc):
!   d/da0 = 1,  d/daa = x,  d/dbb = -exp(cc*x),  d/dcc = -bb*x*exp(cc*x).
! These are the columns of the fitted-value Jacobian used by the
! Gauss-Newton step, NOT the cost gradient in ExpModelGrad (d(RSS)/dp).
! -----------------------------------------------------------------------
pure subroutine ExpModelJacobian(n, x, p, jac)
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: x(n)
    real(real64), intent(in) :: p(NPARAM)
    real(real64), intent(out) :: jac(n, NPARAM)

    real(real64) :: ecx(n)

    ecx = exp(p(4)*x)
    jac(:, 1) = 1.0_real64
    jac(:, 2) = x
    jac(:, 3) = -ecx
    jac(:, 4) = -p(3) * x * ecx
end subroutine ExpModelJacobian

! -----------------------------------------------------------------------
! ExpModelGrad -- Analytic gradient of `ExpModelCost` w.r.t. (a0, aa, bb, cc).
!
!   g1 = -2 sum(e)             ( d/da0 )
!   g2 = -2 sum(e*x)           ( d/daa )
!   g3 = +2 sum(e*exp(cc*x))   ( d/dbb )
!   g4 = +2 sum(e*bb*x*exp(cc*x)) ( d/dcc )
! with e = y - (a0 + aa*x - bb*exp(cc*x)).
! -----------------------------------------------------------------------
pure subroutine ExpModelGrad(n, x, y, p, g)
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: x(n), y(n)
    real(real64), intent(in) :: p(NPARAM)
    real(real64), intent(out) :: g(NPARAM)

    real(real64) :: e(n), ecx(n)

    ecx = exp(p(4)*x)
    e = y - (p(1) + p(2)*x - p(3)*ecx)
    g(1) = -2.0_real64 * sum(e)
    g(2) = -2.0_real64 * sum(e*x)
    g(3) =  2.0_real64 * sum(e*ecx)
    g(4) =  2.0_real64 * sum(e * p(3) * x * ecx)
end subroutine ExpModelGrad

! -----------------------------------------------------------------------
! ExpFitLBFGSB -- Drive the vendored setulb reverse-communication loop to
! minimise `ExpModelCost` under box bounds. Stateless: all L-BFGS-B
! working storage is local to this call.
!
! On return `coef` holds the final iterate, `fval` its cost, and
! `convergence` the termination code:
!   0  -> converged (setulb task 'CONVERGENCE...')
!   1  -> iteration limit reached (iter > maxit)
!   51 -> warning termination (setulb task 'WARNING...')
!   52 -> abnormal/other termination (incl. 'ABNORMAL...'/'ERROR...')
! `status` is 0 when the optimizer ran, or 1 when it could not be
! started at all (invalid sizes / workspace allocation failure), in
! which case `coef` = `seed`, `fval` = 0 and `convergence` = 52.
! -----------------------------------------------------------------------
subroutine ExpFitLBFGSB(&
    ndata, x, y, seed, lower, upper, maxit, lmm, factr, pgtol, coef, fval, convergence, status &
)
    integer(int32), intent(in) :: ndata, maxit, lmm
    real(real64), intent(in) :: x(ndata), y(ndata)
    real(real64), intent(in) :: seed(NPARAM), lower(NPARAM), upper(NPARAM)
    real(real64), intent(in) :: factr, pgtol
    real(real64), intent(out) :: coef(NPARAM)
    real(real64), intent(out) :: fval
    integer(int32), intent(out) :: convergence, status

    ! setulb reverse-communication state (written by setulb, opaque here)
    character(len=60) :: task, csave
    logical :: lsave(4)
    integer(int32) :: isave(44)
    real(real64) :: dsave(29)

    ! problem description and working storage for setulb
    integer(int32) :: nbd(NPARAM)
    integer(int32) :: iwa(3*NPARAM)
    real(real64), allocatable :: wa(:)
    real(real64) :: g(NPARAM), f
    integer(int32) :: iter, ierr

    ! iprint < 0 suppresses the L-BFGS-B iteration report and the iterate.dat
    ! file. It does NOT cover two unguarded write(6,*) calls in the reference
    ! sources; those are commented out in the vendored copy (see README.vendored),
    ! without which this would print to stdout from inside the library.
    integer(int32), parameter :: IPRINT = -1

    convergence = 52
    status = 0
    fval = 0.0_real64
    coef = seed
    f = 0.0_real64
    g = 0.0_real64

    if (ndata < 1 .or. lmm < 1 .or. maxit < 0) then
        status = 1
        return
    end if

    ! setulb workspace: wa(2*m*n + 5*n + 11*m*m + 8*m), iwa(3*n), n=NPARAM,
    ! m=lmm (see the setulb header in vendored/lbfgsb.f).
    allocate(wa(2*lmm*NPARAM + 5*NPARAM + 11*lmm*lmm + 8*lmm), stat=ierr)
    if (ierr /= 0) then
        status = 1
        return
    end if

    ! nbd(i) = 2: every parameter has both a lower and an upper bound.
    nbd = 2
    iter = 0
    task = 'START'

    ! Reverse-communication loop.
    do
        call setulb(&
            NPARAM, lmm, coef, lower, upper, nbd, f, g, factr, pgtol, wa, iwa, task, IPRINT, &
            csave, lsave, isave, dsave &
        )

        if (task(1:2) == 'FG') then
            ! setulb requests a function + gradient evaluation at `coef`.
            call ExpModelCost(ndata, x, y, coef, f)
            call ExpModelGrad(ndata, x, y, coef, g)
        else if (task(1:5) == 'NEW_X') then
            ! A new outer iterate is available; count it and enforce maxit.
            ! Callers should pass maxit = 300. The factr*epsmch stop uses
            ! epsilon(one) = 2.22e-16, which needs ~101 outer iterations on
            ! representative data; a limit of 100 stops one step short.
            iter = iter + 1
            if (iter > maxit) then
                convergence = 1
                exit
            end if
        else if (task(1:4) == 'WARN') then
            convergence = 51
            exit
        else if (task(1:4) == 'CONV') then
            convergence = 0
            exit
        else
            ! 'ERROR...', 'ABNORMAL...', or any unexpected terminal task.
            convergence = 52
            exit
        end if
    end do

    fval = f
    deallocate(wa)
end subroutine ExpFitLBFGSB

! -----------------------------------------------------------------------
! ExpFitGaussNewtonBounded -- Projected (box-constrained) Gauss-Newton
! fit of the exp-deceleration model. It runs the same Gauss-Newton
! iteration as the unbounded ExpFitGaussNewton below, but every trial
! iterate is PROJECTED onto the box [lower, upper]: each candidate p +
! factor*delta is clamped component-wise into the bounds. Where the
! optimum is interior to the box -- the meteor speed-fit case -- that is
! exact; it can differ from an ideal bounded solver only when a bound is
! active at the optimum.
!
! Each iteration builds the analytic Jacobian J = d(model)/dp (columns
! 1, x, -exp(cc*x), -bb*x*exp(cc*x); see ExpModelJacobian) and the residual
! r = y - model(p), then QR-factorises J. The rotated residual rr = Q^T r
! gives the relative-offset convergence criterion, evaluated BEFORE
! stepping:
!     conv = sqrt( sum(rr(1:p)^2) / sum(rr(p+1:n)^2) ) < tol
! rr(1:p) is the residual
! projected onto the Jacobian column space -- the improvement a Gauss-Newton
! step could still achieve -- and rr(p+1:n) is the orthogonal complement, the
! residual that remains at the linearised solution. A zero denominator is
! split on the numerator: 0/0 is an exact interpolating fit and converges,
! while a positive numerator over a zero denominator is +Inf and does not.
! The step solves
! R*delta = rr(1:p) (i.e. delta = (J^T J)^-1 J^T r) and is halved (starting
! from the running factor) until the residual sum of squares does not
! increase; if the factor falls below minFactor the fit is declared
! non-converged.
!
! Convergence codes: 0 converged; 1 iteration limit reached; 2 step factor
! below minFactor; 3 singular Jacobian / LAPACK failure. On any
! non-zero code `coef` is the last real iterate (never garbage) and `se` is
! zero. `sigma` is sqrt(RSS/(n-4)) at that iterate for codes 1 and 2; for
! code 3 it is zero when the failure happened DURING the iteration (no
! trustworthy iterate was reached), but non-zero in the one case where code 3
! is raised AFTER convergence -- an inversion failure while forming the
! standard errors, where the coefficients and sigma are sound and only `se` is
! unavailable. A caller that needs to tell those apart should read sigma:
! code 3 with sigma > 0 is "fit is good, no standard errors".
! On success se(i) = sqrt(sigma^2 * (J^T J)^-1_ii) with
! sigma = sqrt(RSS/(n-4)).
! -----------------------------------------------------------------------
subroutine ExpFitGaussNewtonBounded(&
    ndata, x, y, seed, lower, upper, maxit, tol, minFactor, coef, se, sigma, convergence, status &
)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    integer(int32), intent(in) :: ndata, maxit
    real(real64), intent(in) :: x(ndata), y(ndata)
    real(real64), intent(in) :: seed(NPARAM)
    real(real64), intent(in) :: lower(NPARAM), upper(NPARAM)
    real(real64), intent(in) :: tol, minFactor
    real(real64), intent(out) :: coef(NPARAM), se(NPARAM), sigma
    integer(int32), intent(out) :: convergence, status

    ! LAPACK QR primitives (resolved via the linked liblapack).
    external :: dgeqrf, dormqr, dtrtrs, dtrtri

    real(real64), allocatable :: jac(:,:), rr(:), resid(:), fvals(:), work(:)
    real(real64) :: p(NPARAM), delta(NPARAM), newp(NPARAM), tau(NPARAM)
    real(real64) :: rinv(NPARAM, NPARAM)
    real(real64) :: dev, newdev, factor, projss, orthss, conv, var, qwork(1)
    integer(int32) :: iter, i, j, info, lwork, dof, ierr
    logical :: converged

    convergence = 1          ! iteration limit unless a branch below says otherwise
    status = 0
    coef = seed
    se = 0.0_real64
    sigma = 0.0_real64
    converged = .false.

    ! A 4-parameter fit needs n > 4 so the (n-4) residual d.o.f. of both sigma and
    ! the relative-offset denominator stay positive; controls must be usable.
    if (ndata < NPARAM + 1 .or. maxit < 1 .or. tol <= 0.0_real64 &
            .or. minFactor <= 0.0_real64) then
        status = 1
        return
    end if
    ! An empty box (any lower bound above its upper bound) has no feasible point.
    if (any(lower > upper)) then
        status = 1
        return
    end if

    allocate(jac(ndata, NPARAM), rr(ndata), resid(ndata), fvals(ndata), stat=ierr)
    if (ierr /= 0) then
        status = 1
        return
    end if

    dof = ndata - NPARAM
    ! Project the seed onto the box so the iteration starts feasible.
    p = min(max(seed, lower), upper)
    coef = p
    call ExpModelCost(ndata, x, y, p, dev)

    ! One workspace query sized for both dgeqrf and dormqr (dimensions are fixed
    ! across iterations); contents of jac/tau/rr are irrelevant to the query.
    tau = 0.0_real64
    jac = 0.0_real64
    rr = 0.0_real64
    call dgeqrf(ndata, NPARAM, jac, ndata, tau, qwork, -1, info)
    lwork = int(qwork(1))
    call dormqr('L', 'T', ndata, 1, NPARAM, jac, ndata, tau, rr, ndata, qwork, -1, info)
    lwork = max(lwork, int(qwork(1)))
    lwork = max(lwork, 1)
    allocate(work(lwork), stat=ierr)
    if (ierr /= 0) then
        deallocate(jac, rr, resid, fvals)
        status = 1
        return
    end if

    factor = 1.0_real64

    do iter = 1, maxit
        ! Jacobian and residual at the current iterate.
        call ExpModelJacobian(ndata, x, p, jac)
        call ExpModelValues(ndata, x, p, fvals)
        resid = y - fvals

        ! QR of J (jac is overwritten with R in its upper triangle and the
        ! Householder vectors below); rr := Q^T * resid.
        call dgeqrf(ndata, NPARAM, jac, ndata, tau, work, lwork, info)
        if (info /= 0) then
            convergence = 3
            coef = p
            exit
        end if
        rr = resid
        call dormqr('L', 'T', ndata, 1, NPARAM, jac, ndata, tau, rr, ndata, work, lwork, info)
        if (info /= 0) then
            convergence = 3
            coef = p
            exit
        end if

        ! Relative-offset criterion sqrt(sum(rr(1:p)^2)/sum(rr(p+1:n)^2)) < tol
        ! on rr = Q^T*resid, checked BEFORE stepping.
        projss = sum(rr(1:NPARAM)**2)
        orthss = sum(rr(NPARAM+1:ndata)**2)
        if (orthss <= 0.0_real64) then
            if (projss <= 0.0_real64) then
                ! 0/0: the residual is identically zero, so the model interpolates
                ! the data exactly and there is genuinely nothing left to gain.
                conv = 0.0_real64
            else
                ! projss > 0 with no orthogonal residual: the residual lies wholly
                ! in the Jacobian column space, so a Gauss-Newton step can still
                ! reduce it, so this does NOT converge -- reporting success
                ! here would claim a solution that is not one.
                conv = huge(1.0_real64)
            end if
        else
            conv = sqrt(projss / orthss)
        end if
        if (conv < tol) then
            converged = .true.
            convergence = 0
            coef = p
            exit
        end if

        ! Gauss-Newton step: solve R*delta = rr(1:p) (R = upper triangle of jac).
        delta = rr(1:NPARAM)
        call dtrtrs('U', 'N', 'N', NPARAM, 1, jac, ndata, delta, NPARAM, info)
        if (info /= 0) then
            convergence = 3
            coef = p
            exit
        end if

        ! Step-halving from the running factor until the RSS does not increase.
        ! Each trial step is PROJECTED onto the box: the Gauss-Newton point
        ! p + factor*delta is clamped component-wise into [lower, upper].
        do
            newp = min(max(p + factor * delta, lower), upper)
            call ExpModelCost(ndata, x, y, newp, newdev)
            if (ieee_is_finite(newdev) .and. newdev <= dev) exit
            factor = factor * 0.5_real64
            if (factor < minFactor) then
                convergence = 2
                coef = p
                sigma = sqrt(dev / real(dof, real64))
                deallocate(jac, rr, resid, fvals, work)
                return
            end if
        end do

        p = newp
        dev = newdev
        ! Relax the step cap back toward a full step after a success.
        factor = min(2.0_real64 * factor, 1.0_real64)
    end do

    if (.not. converged) then
        ! Iteration limit (code 1) or a code-3 exit above: report a finite sigma
        ! at the last iterate; coef was already set on the failure path.
        if (convergence == 1) then
            coef = p
            sigma = sqrt(dev / real(dof, real64))
        end if
        deallocate(jac, rr, resid, fvals, work)
        return
    end if

    ! Converged: sigma and per-coefficient SEs from sigma^2 * (J^T J)^-1.
    ! At convergence jac still holds the QR of J at `coef`; (J^T J)^-1 = (R^T R)^-1
    ! = Rinv * Rinv^T, so (J^T J)^-1_ii = sum_j Rinv(i,j)^2 (Rinv upper triangular).
    var = dev / real(dof, real64)
    sigma = sqrt(var)

    rinv = 0.0_real64
    do j = 1, NPARAM
        do i = 1, j
            rinv(i, j) = jac(i, j)
        end do
    end do
    call dtrtri('U', 'N', NPARAM, rinv, NPARAM, info)
    if (info /= 0) then
        ! R (the QR factor) is singular: SEs are undefined. Keep the converged
        ! coefficients and sigma but flag the linear-algebra failure.
        se = 0.0_real64
        convergence = 3
        deallocate(jac, rr, resid, fvals, work)
        return
    end if
    do i = 1, NPARAM
        se(i) = sqrt(var * sum(rinv(i, i:NPARAM)**2))
    end do

    deallocate(jac, rr, resid, fvals, work)
end subroutine ExpFitGaussNewtonBounded

! -----------------------------------------------------------------------
! ExpFitGaussNewton -- Unbounded Gauss-Newton fit of the exp-deceleration
! model. Stateless: all working storage is local to this call.
!
! Each iteration builds the analytic Jacobian J = d(model)/dp (columns
! 1, x, -exp(cc*x), -bb*x*exp(cc*x); see ExpModelJacobian) and the residual
! r = y - model(p), then QR-factorises J. The rotated residual rr = Q^T r
! gives the relative-offset convergence criterion, evaluated BEFORE
! stepping:
!     conv = sqrt( sum(rr(1:p)^2) / sum(rr(p+1:n)^2) ) < tol
! The step solves
! R*delta = rr(1:p) (i.e. delta = (J^T J)^-1 J^T r) and is halved (starting
! from the running factor) until the residual sum of squares does not
! increase; if the factor falls below minFactor the fit is declared
! non-converged.
!
! Convergence codes: 0 converged; 1 iteration limit reached; 2 step factor
! below minFactor; 3 singular Jacobian / LAPACK failure. On any
! non-zero code `coef` is the last real iterate (never garbage) and `se` is
! zero. `sigma` is sqrt(RSS/(n-4)) at that iterate for codes 1 and 2; for
! code 3 it is zero when the failure happened DURING the iteration (no
! trustworthy iterate was reached), but non-zero in the one case where code 3
! is raised AFTER convergence -- an inversion failure while forming the
! standard errors, where the coefficients and sigma are sound and only `se` is
! unavailable. A caller that needs to tell those apart should read sigma:
! code 3 with sigma > 0 is "fit is good, no standard errors".
! On success se(i) = sqrt(sigma^2 * (J^T J)^-1_ii) with
! sigma = sqrt(RSS/(n-4)).
!
! Implemented by delegating to ExpFitGaussNewtonBounded with an unbounded
! (-inf, +inf) box, so there is one source of truth for the Gauss-Newton
! machinery. The projection is then a genuine no-op: for every finite
! iterate min(max(v, -inf), +inf) = v, and an overflowing trial step stays
! +/-inf and is rejected by the same finiteness guard the unbounded loop
! used -- so this path is numerically identical to a hand-written unbounded
! Gauss-Newton.
! -----------------------------------------------------------------------
subroutine ExpFitGaussNewton(&
    ndata, x, y, seed, maxit, tol, minFactor, coef, se, sigma, convergence, status &
)
    use, intrinsic :: ieee_arithmetic, only : ieee_value, &
        ieee_negative_inf, ieee_positive_inf
    integer(int32), intent(in) :: ndata, maxit
    real(real64), intent(in) :: x(ndata), y(ndata)
    real(real64), intent(in) :: seed(NPARAM)
    real(real64), intent(in) :: tol, minFactor
    real(real64), intent(out) :: coef(NPARAM), se(NPARAM), sigma
    integer(int32), intent(out) :: convergence, status

    real(real64) :: lower(NPARAM), upper(NPARAM)

    lower = ieee_value(0.0_real64, ieee_negative_inf)
    upper = ieee_value(0.0_real64, ieee_positive_inf)
    call ExpFitGaussNewtonBounded(&
        ndata, x, y, seed, lower, upper, maxit, tol, minFactor, coef, se, sigma, convergence, &
        status &
    )
end subroutine ExpFitGaussNewton

end module MTExpFit
