*DECK DVODPK
      SUBROUTINE DVODPK (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1  ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, PSOL, MF, RPAR, IPAR)
      EXTERNAL F, JAC, PSOL
      DOUBLE PRECISION Y, T, TOUT, RTOL, ATOL, RWORK, RPAR
      INTEGER NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, MF, IPAR
      DIMENSION Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW),
     1          RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! This is the 26 April 2002 version of
! DVODPK: Variable-coefficient Ordinary Differential equation solver
!         with the Preconditioned Krylov method GMRES for the solution
!         of linear systems.
!
! This version is in double precision.
!
! DVODPK solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
! DVODPK is a package based on the VODE and LSODPK packages, and on
! the October 23, 1978 version of the ODEPACK user interface standard,
! with minor modifications.
!-----------------------------------------------------------------------
! Authors:
!               Alan C. Hindmarsh and Peter N. Brown
!               Center for Applied Scientific Computing, L-561
!               Lawrence Livermore National Laboratory
!               Livermore, CA 94551
! and
!               George D. Byrne
!               Illinois Institute of Technology
!               Chicago, IL 60616
!-----------------------------------------------------------------------
! References:
! 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE, A Variable-
!    Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10  (1989),
!    pp., 1038-1051.  Also LLNL report UCRL-98412, June 1988.
! 2. P. N. Brown and A. C. Hindmarsh, "Reduced Storage Matrix Methods
!    in Stiff ODE Systems," J. Appl. Math. & Comp., 31 (1989), pp.40-91.
!    Also LLNL report UCRL-95088, Rev. 1, June 1987.
! 3. G. D. Byrne, "Pragmatic Experiments with Krylov Methods in the
!    Stiff ODE Setting," Computational Ordinary Differential Equations,
!    J. Cash and I. Gladwell, eds., Oxford Univ. Press, Oxford, 1992,
!    pp. 323-356.
!-----------------------------------------------------------------------
! Introduction.
!
! This is a modification of the VODE package which incorporates
! the preconditioned Krylov subspace iterative method SPIGMR for the
! linear algebraic systems that arise in the case of stiff systems.
! SPIGMR denotes a scaled preconditioned incomplete version of the
! GMRES (Generalized Minimum Residual) method.
!
! The linear systems that are solved have the form
!   A * x  = b ,  where  A = I - hrl1 * (df/dy) .
! here hrl1 is a scalar, I is the identity matrix, and df/dy is the
! Jacobian matrix of partial derivatives of f with respect to y
! (an NEQ by NEQ matrix).
!
! The particular Krylov method is chosen by setting the second digit,
! MITER, in the method flag MF.
! Currently, the values of MITER have the following meanings:
!
!          1 means SPIGMR, a scaled, preconditioned, incomplete version
!            of GMRES, a generalized minimum residual method.
!            This is the best choice in general.
!
!          9 means that only a user-supplied matrix P (approximating A)
!            will be used, with no Krylov iteration done internally to
!            DVODPK.  This option allows the user to provide the
!            complete linear system solution algorithm, if desired.
!
! The user can apply preconditioning to the linear system A*x = b,
! by means of arbitrary matrices (the preconditioners).
!
!     In the case of SPIGMR, one can apply left and right
! preconditioners P1 and P2, and the basic iterative method is then
! applied to the matrix (P1-inverse)*A*(P2-inverse) instead of to the
! matrix A.  The product P1*P2 should be an approximation to A
! such that linear systems with P1 or P2 are easier to solve than with
! A alone.  Preconditioning from the left only or right only means using
! P2 = I  or  P1 = I, respectively.
!
!     If the Jacobian  J = df/dy  splits in a natural way into a sum
! J = J1 + J2, then one possible choice of preconditioners is
!            P1 = I - hrl1 * J1  and  P2 = I - hrl1 * J2
! provided each of these is easy to solve (or to approximately solve).
!
! NOTE:  To achieve an efficient solution, the preconditioned Krylov
! methods in DVODPK generally require a thoughtful choice of
! preconditioners.  If the ODE system produces linear systems that are
! not amenable to solution by such iterative methods, the cost can be
! higher than with a solver that uses sparse direct methods.  However,
! for many systems, careful use of DVODPK can be highly effective.
!
! See Ref. 2 for more details on the methods and applications.
!-----------------------------------------------------------------------
! Summary of usage.
!
! Communication between the user and the DVODPK package, for normal
! situations, is summarized here.  This summary describes only a subset
! of the full set of options available.  See full description (below)
! for details, including optional communication, nonstandard options,
! and instructions for special situations.  See also the example
! program embedded in the comments below.
!
! A. First provide a subroutine of the form
!     SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
!     DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), RPAR(*)
!     INTEGER IPAR(*)
! which supplies the vector function f by loading YDOT(i) with f(i).
!
! B. Next determine (or guess) whether or not the problem is stiff.
! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest.  If the problem is nonstiff,
! use method flag MF = 10.  If it is stiff, MF should be 21.
!
! The following four parameters must also be set.
!  IWORK(1) = LWP  = length of real array WP for preconditioning.
!  IWORK(2) = LIWP = length of integer array IWP for preconditioning.
!  IWORK(3) = JPRE = preconditioner type flag:
!                  = 0 for no preconditioning (P1 = P2 = I)
!                  = 1 for left-only preconditioning (P2 = I)
!                  = 2 for right-only preconditioning (P1 = I)
!                  = 3 for two-sided preconditioning
!  IWORK(4) = JACFLG = flag for whether JAC is called.
!                    = 0 if JAC is not to be called,
!                    = 1 if JAC is to be called.
!  Use JACFLG = 1 if JAC computes any nonconstant data for use in
!  preconditioning, such as Jacobian elements.  See next paragraph.
!  The arrays WP and IWP are work arrays under the user's control,
!  for use in the routines that perform preconditioning operations.
!
! C. If the problem is stiff, you must supply two routines that deal
! with the preconditioning of the linear systems to be solved.
! These are as follows:
!
!     SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY, V, HRL1, WP, IWP,
!    1                IER, RPAR, IPAR)
!     DOUBLE PRECISION T, Y(NEQ), YSV(NEQ), REWT(NEQ), FTY(NEQ), V(NEQ),
!    1                 HRL1, WP(*), RPAR(*)
!     INTEGER IWP(*), IPAR(*)
!
!        This routine is optional, and is to evaluate and preprocess
!     any parts of the Jacobian matrix df/dy involved in the
!     preconditioners P1 and P2.
!     The Y and FTY arrays contain the current values of y and f(t,y),
!     respectively, and YSV also contains the current value of y.
!     The array V is work space of length NEQ.
!     JAC must multiply all computed Jacobian elements by the scalar
!     -hrl1, add the identity matrix I, and do any factorization
!     operations called for, in preparation for solving linear systems
!     with a coefficient matrix of P1 or P2.  The matrix P1*P2 should
!     be an approximation to  I - hrl1 * (df/dy), where hrl1 is a
!     scalar stored in HRL1.
!     JAC should return IER = 0 if successful, and IER .ne. 0 if not.
!     (If IER .ne. 0, a smaller time step will be tried.)
!
!     SUBROUTINE PSOL (NEQ, T, Y, FTY, WK, HRL1, WP, IWP, B, LR, IER,
!    1                 RPAR, IPAR)
!     DOUBLE PRECISION T, Y(NEQ), FTY(NEQ), WK(NEQ), HRL1, WP(*),
!    1                 B(NEQ), RPAR(*)
!     INTEGER IWP(*), IPAR(*)
!
!        This routine must solve a linear system with b (stored in B)
!     as right-hand side and one of the preconditioning matrices, P1 or
!     P2, as coefficient matrix, and return the solution vector in B.
!     LR is a flag concerning left vs. right preconditioning, input
!     to PSOL.  PSOL is to use P1 if LR = 1, and P2 if LR = 2.
!
!        PSOL can use data generated in the JAC routine and stored in
!     WP and IWP.  WK is a work array of length NEQ.
!     The argument HRL1 is the current value of the scalar appearing
!     in the linear system.  If the old value, at the time of the last
!     JAC call, is needed, it must have been saved by JAC in WP.
!     on return, PSOL should set the error flag  IER as follows:
!        IER = 0 if PSOL was successful,
!        IER .gt. 0 if a recoverable error occurred, meaning that the
!              time step will be retried,
!        IER .lt. 0 if an unrecoverable error occurred, meaning that the
!              solver is to stop immediately.
!
! D. Write a main program which calls subroutine DVODPK once for
! each point at which answers are desired.  This should also provide
! for possible use of logical unit 6 for output of error messages
! by DVODPK.  on the first call to DVODPK, supply arguments as follows:
! F      = name of subroutine for right-hand side vector f.
!          This name must be declared EXTERNAL in calling program.
! NEQ    = number of first order ODEs.
! Y      = array of initial values, of length NEQ.
! T      = the initial value of the independent variable.
! TOUT   = first point where output is desired (.ne. T).
! ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
! RTOL   = relative tolerance parameter (scalar).
! ATOL   = absolute tolerance parameter (scalar or array).
!          The estimated local error in Y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL (or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control.  Caution: Actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of Y at t = TOUT.
! ISTATE = integer flag (input and output).  Set ISTATE = 1.
! IOPT   = 0 to indicate no optional input used.
! RWORK  = real work array of length at least:
!             20 + 16*NEQ           for MF = 10,
!             61 + 17*NEQ + LWP     for MF = 21.
! LRW    = declared length of RWORK (in user's DIMENSION statement).
! IWORK  = integer work array of length at least:
!             30            for MF = 10,
!             30 + LIWP     for MF = 21.
! LIW    = declared length of IWORK (in user's DIMENSION statement).
! JAC,PSOL = names of subroutines for preconditioning.  These names
!            must be declared EXTERNAL in the user's calling program.
! MF     = method flag.  Standard values are:
!          10 for nonstiff (Adams) method.
!          21 for stiff (BDF) method, with SPIGMR.
!
! RPAR, IPAR  User-specified arrays used to communicate real and integer
!             parameters (respectively) to user-supplied subroutines.
!             to user-supplied subroutines.  If RPAR is a vector, then
!             it must be dimensioned in the user's main program.  If it
!             is unused or a scalar, then it need not be dimensioned.
!
! IPAR     User-specified array used to communicate integer parameter
!          to user-supplied subroutines.  The comments on dimensioning
!          RPAR apply to IPAR.
!
! Note that the user's main (calling) program must declare arrays
! Y, RWORK, IWORK, and possibly ATOL, RPAR, and IPAR.
!
! E. The output from the first call (or any call) is:
!      Y = array of computed values of y(t) vector.
!      T = corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if DVODPK was successful, negative otherwise.
!         -1 means excess work done on this call (perhaps wrong MF).
!         -2 means excess accuracy requested (tolerances too small).
!         -3 means illegal input detected (see printed message).
!         -4 means repeated error test failures (check all input).
!         -5 means repeated convergence failures (perhaps bad JAC
!            or PSOL routine supplied or wrong choice of MF or
!            tolerances, or this solver is inappropriate).
!         -6 means error weight became zero during problem. (Solution
!            component i vanished, and ATOL or ATOL(i) = 0.)
!         -7 means an unrecoverable error occurred in JAC or PSOL.
!
! F. To continue the integration after a successful return, simply
! reset TOUT and call DVODPK again.  No other parameters need be reset.
!
!-----------------------------------------------------------------------
! Example problem.
! An ODE system is generated from the following 2-species diurnal
! kinetics advection-diffusion PDE system in 2 space dimensions:
!
! dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz)
!                 + Ri(c1,c2,t)      for i = 1,2,   where
!   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
!   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
!   Kv(z) = Kv0*exp(z/5) ,
! Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
! vary diurnally.   The problem is posed on the square
!   0 .le. x .le. 20,    30 .le. z .le. 50   (all in km),
! with homogeneous Neumann boundary conditions, and for time t in
!   0 .le. t .le. 86400 sec (1 day).
! The PDE system is treated by central differences on a uniform
! 10 x 10 mesh, with simple polynomial initial profiles.
! The problem is solved with DVODPK, with the BDF/GMRES method and
! the block-diagonal part of the Jacobian as a left preconditioner.
!-----------------------------------------------------------------------
!      EXTERNAL FEX, JACBD, SOLBD
!      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3,A4, OM, C3, DZ, HDCO,VDCO,HACO
!      COMMON /PCOM/ Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO,MX,MZ,MM
!      DOUBLE PRECISION ATOL, AVDIM, CX, CZ, DKH, DKV0, DX, FLOOR,
!     1     HALFDA, PI, RPAR, RTOL, RWORK, T, TOUT, TWOHR, VEL, X, Y, Z
!      DIMENSION Y(2,10,10), RWORK(3861), IWORK(230)
!      DATA DKH/4.0D-6/, VEL/0.001D0/, DKV0/1.0D-8/, HALFDA/4.32D4/,
!     1  PI/3.1415926535898D0/, TWOHR/7200.0D0/, RTOL/1.0D-5/,
!     2  FLOOR/100.0D0/, LRW/3861/, LIW/230/, MF/21/, JPRE/1/, JACFLG/1/
!
! Load Common block of problem parameters.
!      MX = 10
!      MZ = 10
!      MM = MX*MZ
!      Q1 = 1.63D-16
!      Q2 = 4.66D-16
!      A3 = 22.62D0
!      A4 = 7.601D0
!      OM = PI/HALFDA
!      C3 = 3.7D16
!      DX = 20.0D0/(MX - 1.0D0)
!      DZ = 20.0D0/(MZ - 1.0D0)
!      HDCO = DKH/DX**2
!      HACO = VEL/(2.0D0*DX)
!      VDCO = (1.0D0/DZ**2)*DKV0
! Set other input arguments.
!      ATOL = RTOL*FLOOR
!      NEQ = 2*MX*MZ
!      IWORK(1) = 4*MX*MZ
!      IWORK(2) = NEQ
!      IWORK(3) = JPRE
!      IWORK(4) = JACFLG
!      T = 0.0D0
!      TOUT = TWOHR
!      ISTATE = 1
! Set initial profiles.
!      DO 20 JZ = 1,MZ
!        Z = 30.0D0 + (JZ - 1.0D0)*DZ
!        CZ = (0.1D0*(Z - 40.0D0))**2
!        CZ = 1.0D0 - CZ + 0.5D0*CZ**2
!        DO 10 JX = 1,MX
!          X = (JX - 1.0D0)*DX
!          CX = (0.1D0*(X - 10.0D0))**2
!          CX = 1.0D0 - CX + 0.5D0*CX**2
!          Y(1,JX,JZ) = 1.0D6*CX*CZ
!          Y(2,JX,JZ) = 1.0D12*CX*CZ
! 10       CONTINUE
! 20     CONTINUE
!
! Loop over output points, call DVODPK, print sample solution values.
!      DO 70 IOUT = 1,12
!        CALL DVODPK (FEX, NEQ, Y, T, TOUT, 1, RTOL, ATOL, 1, ISTATE, 0,
!     1            RWORK, LRW, IWORK, LIW, JACBD, SOLBD, MF, RPAR, IPAR)
!        WRITE(6,50) T,IWORK(11),IWORK(14),RWORK(11)
! 50     FORMAT(/' t =',D10.2,5X,'no. steps =',I5,
!     1                      '   order =',I3,'   stepsize =',D10.2)
!        WRITE(6,60) Y(1,1,1), Y(1,5,5), Y(1,10,10),
!     1              Y(2,1,1), Y(2,5,5), Y(2,10,10)
! 60     FORMAT('  c1 (bot.left/middle/top rt.) =',3D12.3/
!     1         '  c2 (bot.left/middle/top rt.) =',3D12.3)
!        IF (ISTATE .NE. 2) STOP
!        TOUT = TOUT + TWOHR
! 70     CONTINUE
!
! Print final statistics.
!      LENRW = IWORK(17)
!      LENIW = IWORK(18)
!      NST = IWORK(11)
!      NFE = IWORK(12)
!      NPE = IWORK(13)
!      NPS = IWORK(24)
!      NNI = IWORK(20)
!      NLI = IWORK(23)
!      AVDIM = REAL(NLI)/REAL(NNI)
!      NCFN = IWORK(21)
!      NCFL = IWORK(25)
!      WRITE (6,80) LENRW,LENIW,NST,NFE,NPE,NPS,NNI,NLI,AVDIM,NCFN,NCFL
! 80   FORMAT(//' Final statistics:'/
!     1 ' RWORK size =',I5,5X,' IWORK size =',I4/
!     2 ' Number of steps        =',I5,5X,'Number of f evals.     =',I5/
!     3 ' Number of prec. evals. =',I5,5X,'Number of prec. solves =',I5/
!     4 ' Number of nonl. iters. =',I5,5X,'Number of lin. iters.  =',I5/
!     5 ' Average Krylov subspace dimension (NLI/NNI)  =',F8.4/
!     6 ' Number of conv. failures:  nonlinear =',I3,'  linear =',I3)
!      STOP
!      END
!
!      SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
!      DOUBLE PRECISION T, Y, YDOT, RPAR
!      DIMENSION Y(2,*), YDOT(2,*)
!      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3,A4, OM, C3, DZ, HDCO,VDCO,HACO
!      COMMON /PCOM/ Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO,MX,MZ,MM
!      DOUBLE PRECISION C1, C2, C1DN, C2DN, C1UP, C2UP, C1LT, C2LT,
!     1    C1RT, C2RT, CZDN, CZUP, HORD1, HORD2, HORAD1, HORAD2,
!     2    QQ1, QQ2, QQ3, QQ4, RKIN1, RKIN2, S, VERTD1, VERTD2, ZDN, ZUP
!
! Set diurnal rate coefficients.
!      S = SIN(OM*T)
!      IF (S .GT. 0.0D0) THEN
!        Q3 = EXP(-A3/S)
!        Q4 = EXP(-A4/S)
!      ELSE
!        Q3 = 0.0D0
!        Q4 = 0.0D0
!      ENDIF
! Loop over all grid points.
!      DO 20 JZ = 1,MZ
!        ZDN = 30.0D0 + (JZ - 1.5D0)*DZ
!        ZUP = ZDN + DZ
!        CZDN = VDCO*EXP(0.2D0*ZDN)
!        CZUP = VDCO*EXP(0.2D0*ZUP)
!        IBLOK0 = (JZ-1)*MX
!        IDN = -MX
!        IF (JZ .EQ. 1) IDN = MX
!        IUP = MX
!        IF (JZ .EQ. MZ) IUP = -MX
!        DO 10 JX = 1,MX
!          IBLOK = IBLOK0 + JX
!          C1 = Y(1,IBLOK)
!          C2 = Y(2,IBLOK)
! Set kinetic rate terms.
!          QQ1 = Q1*C1*C3
!          QQ2 = Q2*C1*C2
!          QQ3 = Q3*C3
!          QQ4 = Q4*C2
!          RKIN1 = -QQ1 - QQ2 + 2.0D0*QQ3 + QQ4
!          RKIN2 = QQ1 - QQ2 - QQ4
! Set vertical diffusion terms.
!          C1DN = Y(1,IBLOK+IDN)
!          C2DN = Y(2,IBLOK+IDN)
!          C1UP = Y(1,IBLOK+IUP)
!          C2UP = Y(2,IBLOK+IUP)
!          VERTD1 = CZUP*(C1UP - C1) - CZDN*(C1 - C1DN)
!          VERTD2 = CZUP*(C2UP - C2) - CZDN*(C2 - C2DN)
! Set horizontal diffusion and advection terms.
!          ILEFT = -1
!          IF (JX .EQ. 1) ILEFT = 1
!          IRIGHT = 1
!          IF (JX .EQ. MX) IRIGHT = -1
!          C1LT = Y(1,IBLOK+ILEFT)
!          C2LT = Y(2,IBLOK+ILEFT)
!          C1RT = Y(1,IBLOK+IRIGHT)
!          C2RT = Y(2,IBLOK+IRIGHT)
!          HORD1 = HDCO*(C1RT - 2.0D0*C1 + C1LT)
!          HORD2 = HDCO*(C2RT - 2.0D0*C2 + C2LT)
!          HORAD1 = HACO*(C1RT - C1LT)
!          HORAD2 = HACO*(C2RT - C2LT)
! Load all terms into YDOT.
!          YDOT(1,IBLOK) = VERTD1 + HORD1 + HORAD1 + RKIN1
!          YDOT(2,IBLOK) = VERTD2 + HORD2 + HORAD2 + RKIN2
! 10       CONTINUE
! 20     CONTINUE
!      RETURN
!      END
!
!      SUBROUTINE JACBD (F, NEQ, T, Y, YSV, REWT, F0, F1, HRL1,
!     1                  BD, IPBD, IER, RPAR, IPAR)
!      EXTERNAL F
!      DOUBLE PRECISION T, Y, YSV, REWT, F0, F1, HRL1, BD, RPAR
!      DIMENSION Y(2, *), YSV(*), REWT(*), F0(*), F1(*), BD(2, 2, *),
!     1          IPBD(2, *)
!      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3,A4, OM, C3, DZ, HDCO,VDCO,HACO
!      COMMON /PCOM/ Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO,MX,MZ,MM
!      DOUBLE PRECISION C1, C2, CZDN, CZUP, DIAG, ZDN, ZUP
!
! Compute diagonal Jacobian blocks, multiplied by -HRL1
!   (using q3 and q4 values computed on last F call).
!      DO 20 JZ = 1,MZ
!        ZDN = 30.0D0 + (JZ - 1.5D0)*DZ
!        ZUP = ZDN + DZ
!        CZDN = VDCO*EXP(0.2D0*ZDN)
!        CZUP = VDCO*EXP(0.2D0*ZUP)
!        DIAG = -(CZDN + CZUP + 2.0D0*HDCO)
!        IBLOK0 = (JZ-1)*MX
!        DO 10 JX = 1,MX
!          IBLOK = IBLOK0 + JX
!          C1 = Y(1,IBLOK)
!          C2 = Y(2,IBLOK)
!          BD(1,1,IBLOK) = -HRL1*( (-Q1*C3 - Q2*C2) + DIAG )
!          BD(1,2,IBLOK) = -HRL1*( -Q2*C1 + Q4 )
!          BD(2,1,IBLOK) = -HRL1*( Q1*C3 - Q2*C2 )
!          BD(2,2,IBLOK) = -HRL1*( (-Q2*C1 - Q4) + DIAG )
! 10       CONTINUE
! 20     CONTINUE
! Add identity matrix and do LU decompositions on blocks.
!      DO 40 IBLOK = 1,MM
!        BD(1,1,IBLOK) = BD(1,1,IBLOK) + 1.0D0
!        BD(2,2,IBLOK) = BD(2,2,IBLOK) + 1.0D0
!        CALL DGEFA (BD(1,1,IBLOK), 2, 2, IPBD(1,IBLOK), IER)
!        IF (IER .NE. 0) RETURN
! 40     CONTINUE
!      RETURN
!      END
!
!      SUBROUTINE SOLBD (NEQ, T, Y, F0, WK, HRL1, BD, IPBD, V, LR, IER,
!     1                  RPAR, IPAR)
!      DOUBLE PRECISION T, Y, F0, WK, HRL1, BD, V, RPAR
!      DIMENSION BD(2,2,*), IPBD(2,*), V(2,*)
!      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3,A4, OM, C3, DZ, HDCO,VDCO,HACO
!      COMMON /PCOM/ Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO,MX,MZ,MM
!
! Solve the block-diagonal system Px = v using LU factors stored in BD
! and pivot data in IPBD, and return the solution in V.
!      IER = 0
!      DO 10 I = 1,MM
!        CALL DGESL (BD(1,1,I), 2, 2, IPBD(1,I), V(1,I), 0)
! 10     CONTINUE
!      RETURN
!      END
!
! The output of this program, on a Cray-1 in single precision,
! is as follows:
!
! t =  7.20e+03     no. steps =  194   order =  5   stepsize =  1.17e+02
!  c1 (bot.left/middle/top rt.) =   1.047e+04   2.964e+04   1.119e+04
!  c2 (bot.left/middle/top rt.) =   2.527e+11   7.154e+11   2.700e+11
!
! t =  1.44e+04     no. steps =  227   order =  5   stepsize =  2.73e+02
!  c1 (bot.left/middle/top rt.) =   6.659e+06   5.316e+06   7.301e+06
!  c2 (bot.left/middle/top rt.) =   2.582e+11   2.057e+11   2.833e+11
!
! t =  2.16e+04     no. steps =  252   order =  5   stepsize =  4.21e+02
!  c1 (bot.left/middle/top rt.) =   2.665e+07   1.036e+07   2.931e+07
!  c2 (bot.left/middle/top rt.) =   2.993e+11   1.028e+11   3.313e+11
!
! t =  2.88e+04     no. steps =  291   order =  4   stepsize =  2.13e+02
!  c1 (bot.left/middle/top rt.) =   8.702e+06   1.292e+07   9.650e+06
!  c2 (bot.left/middle/top rt.) =   3.380e+11   5.029e+11   3.751e+11
!
! t =  3.60e+04     no. steps =  321   order =  5   stepsize =  9.90e+01
!  c1 (bot.left/middle/top rt.) =   1.404e+04   2.029e+04   1.561e+04
!  c2 (bot.left/middle/top rt.) =   3.387e+11   4.894e+11   3.765e+11
!
! t =  4.32e+04     no. steps =  374   order =  4   stepsize =  4.44e+02
!  c1 (bot.left/middle/top rt.) =  -5.457e-09  -4.365e-09  -6.182e-09
!  c2 (bot.left/middle/top rt.) =   3.382e+11   1.355e+11   3.804e+11
!
! t =  5.04e+04     no. steps =  393   order =  5   stepsize =  5.22e+02
!  c1 (bot.left/middle/top rt.) =   3.396e-12   2.798e-12   3.789e-12
!  c2 (bot.left/middle/top rt.) =   3.358e+11   4.930e+11   3.864e+11
!
! t =  5.76e+04     no. steps =  407   order =  5   stepsize =  3.54e+02
!  c1 (bot.left/middle/top rt.) =   7.738e-12   6.455e-12   8.598e-12
!  c2 (bot.left/middle/top rt.) =   3.320e+11   9.650e+11   3.909e+11
!
! t =  6.48e+04     no. steps =  419   order =  5   stepsize =  5.90e+02
!  c1 (bot.left/middle/top rt.) =  -2.018e-11  -1.680e-11  -2.243e-11
!  c2 (bot.left/middle/top rt.) =   3.313e+11   8.922e+11   3.963e+11
!
! t =  7.20e+04     no. steps =  432   order =  5   stepsize =  5.90e+02
!  c1 (bot.left/middle/top rt.) =  -2.837e-11  -2.345e-11  -3.166e-11
!  c2 (bot.left/middle/top rt.) =   3.330e+11   6.186e+11   4.039e+11
!
! t =  7.92e+04     no. steps =  444   order =  5   stepsize =  5.90e+02
!  c1 (bot.left/middle/top rt.) =  -4.861e-14  -4.433e-14  -5.162e-14
!  c2 (bot.left/middle/top rt.) =   3.334e+11   6.669e+11   4.120e+11
!
! t =  8.64e+04     no. steps =  456   order =  5   stepsize =  5.90e+02
!  c1 (bot.left/middle/top rt.) =   2.511e-15   2.071e-15   2.802e-15
!  c2 (bot.left/middle/top rt.) =   3.352e+11   9.107e+11   4.163e+11
!
!
! Final statistics:
! RWORK size = 3861      IWORK size = 230
! Number of steps        =  456     Number of f evals.     = 1317
! Number of prec. evals. =   82     Number of prec. solves = 1226
! Number of nonl. iters. =  571     Number of lin. iters.  =  743
! Average Krylov subspace dimension (NLI/NNI)  =  1.3012
! Number of conv. failures:  nonlinear =  0  linear =  0
!-----------------------------------------------------------------------
! Full description of user interface to DVODPK.
!
! The user interface to DVODPK consists of the following parts.
!
! i.   The call sequence to subroutine DVODPK, which is a driver
!      routine for the solver.  This includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      Following these descriptions are
!        * a description of optional input available through the
!          call sequence,
!        * a description of optional output (in the work arrays), and
!        * instructions for interrupting and restarting a solution.
!
! ii.  Descriptions of other routines in the DVODPK package that may be
!      (optionally) called by the user.  These provide the ability to
!      alter error message handling, save and restore the internal
!      COMMON, and obtain specified derivatives of the solution y(t).
!
! iii. Descriptions of COMMON blocks to be declared in overlay
!      or similar environments.
!
! iv.  Description of two routines in the DVODPK package, either of
!      which the user may replace with the user's own version, if
!      desired.  These relate to the measurement of errors.
!
!-----------------------------------------------------------------------
! Part i.  Call Sequence.
!
! The call sequence parameters used for input only are
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW,
!     JAC, PSOL, MF,
! and those used for both input and output are
!     Y, T, ISTATE.
! The work arrays RWORK and IWORK are also used for conditional and
! optional input and optional output.  (The term output here refers
! to the return from subroutine DVODPK to the user's calling program.)
!
! The legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by ISTATE = 3 in the input.
!
! The descriptions of the call arguments are as follows.
!
! F      = The name of the user-supplied subroutine defining the
!          ODE system.  The system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  Subroutine F is to
!          compute the function f.  It is to have the form
!               SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR)
!               DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), RPAR(*)
!               INTEGER IPAR(*)
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!          is output.  Y and YDOT are arrays of length NEQ.
!          (In the DIMENSION statement above, NEQ  can be replaced by
!          *  to make  Y  and  YDOT  assumed size arrays.)
!          Subroutine F should not alter Y or T.
!          F must be declared EXTERNAL in the calling program.
!
!          Subroutine F may access user-defined real and integer
!          work arrays RPAR and IPAR, which are to be dimensioned
!          in the user's calling (main) program.
!
!          If quantities computed in the F routine are needed
!          externally to DVODPK, an extra call to F should be made
!          for this purpose, for consistent and accurate results.
!          If only the derivative dy/dt is needed, use DVINDY instead.
!
! NEQ    = The size of the ODE system (number of first order
!          ordinary differential equations).  Used only for input.
!          NEQ may not be increased during the problem, but
!          can be decreased (with ISTATE = 3 in the input).
!
! Y      = A real array for the vector of dependent variables, of
!          length NEQ or more.  Used for both input and output on the
!          first call (ISTATE = 1), and only for output on other calls.
!          On the first call, Y must contain the vector of initial
!          values.  In the output, Y contains the computed solution
!          evaluated at T.  If desired, the Y array may be used
!          for other purposes between calls to the solver.
!
!          This array is passed as the Y argument in all calls to
!          F, JAC, and PSOL.
!
! T      = The independent variable.  In the input, T is used only on
!          the first call, as the initial point of the integration.
!          In the output, after each call, T is the value at which a
!          computed solution Y is evaluated (usually the same as TOUT).
!          On an error return, T is the farthest point reached.
!
! TOUT   = The next value of t at which a computed solution is desired.
!          Used only for input.
!
!          When starting the problem (ISTATE = 1), TOUT may be equal
!          to T for one call, then should .ne. T for the next call.
!          For the initial T, an input value of TOUT .ne. T is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem.  Integration in either direction
!          (forward or backward in t) is permitted.
!
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!          the first call (i.e. the first call with TOUT .ne. T).
!          Otherwise, TOUT is required on every call.
!
!          If ITASK = 1, 3, or 4, the values of TOUT need not be
!          monotone, but a value of TOUT which backs up is limited
!          to the current internal t interval, whose endpoints are
!          TCUR - HU and TCUR.  (See optional output, below, for
!          TCUR and HU.)
!
! ITOL   = An indicator for the type of error control.  See
!          description below under ATOL.  Used only for input.
!
! RTOL   = A relative error tolerance parameter, either a scalar or
!          an array of length NEQ.  See description below under ATOL.
!          Input only.
!
! ATOL   = An absolute error tolerance parameter, either a scalar or
!          an array of length NEQ.  Input only.
!
!          The input parameters ITOL, RTOL, and ATOL determine
!          the error control performed by the solver.  The solver will
!          control the vector e = (e(i)) of estimated local errors
!          in Y, according to an inequality of the form
!                      rms-norm of ( e(i)/EWT(i) )   .le.   1,
!          where       EWT(i) = RTOL(i)*abs(Y(i)) + ATOL(i),
!          and the rms-norm (root-mean-square norm) here is
!          rms-norm(v) = sqrt(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
!          is a vector of weights which must always be positive, and
!          the values of RTOL and ATOL should all be non-negative.
!          The following table gives the types (scalar/array) of
!          RTOL and ATOL, and the corresponding form of EWT(i).
!
!             ITOL    RTOL       ATOL          EWT(i)
!              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
!              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
!              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
!              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
!
!          When either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!
!          If none of the above choices (with ITOL, RTOL, and ATOL
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of EWT and/or for
!          the norm calculation.  See Part iv below.
!
!          If global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of RTOL and ATOL (i.e. of EWT) should be scaled
!          down uniformly.
!
! ITASK  = An index specifying the task to be performed.
!          Input only.  ITASK has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = TOUT (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = TOUT and return.
!          4  means normal computation of output values of y(t) at
!             t = TOUT but without overshooting t = TCRIT.
!             TCRIT must be input as RWORK(1).  TCRIT may be equal to
!             or beyond TOUT, but not behind it in the direction of
!             integration.  This option is useful if the problem
!             has a singularity at or beyond t = TCRIT.
!          5  means take one step, without passing TCRIT, and return.
!             TCRIT must be input as RWORK(1).
!
!          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!          (within roundoff), it will return T = TCRIT (exactly) to
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!          in which case answers at T = TOUT are returned first).
!
! ISTATE = an index used for input and output to specify the
!          the state of the calculation.
!
!          In the input, the values of ISTATE are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  See note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly TOUT and ITASK.
!             (If ITOL, RTOL, and/or ATOL are changed between calls
!             with ISTATE = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             TOUT and ITASK.  Changes are allowed in
!             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
!             and any of the optional input except H0.
!
!          Note:  A preliminary call with TOUT = T is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (Such a call is sometimes useful to include
!          the initial conditions in the output.)
!          Thus the first call for which TOUT .ne. T requires
!          ISTATE = 1 in the input.
!
!          In the output, ISTATE has the following values and meanings.
!           1  means nothing was done, as TOUT was equal to T with
!              ISTATE = 1 in the input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than MXSTEP
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as T.  (MXSTEP is an optional input
!              and is normally 500.)  To continue, the user may
!              simply reset ISTATE to a value .gt. 1 and call again.
!              (The excess work step counter will be reset to 0.)
!              In addition, the user may increase MXSTEP to avoid
!              this error return.  (See optional input below.)
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  This was detected before
!              completing the requested task, but the integration
!              was successful as far as T.  To continue, the tolerance
!              parameters must be reset, and ISTATE must be set
!              to 3.  The optional output TOLSF may be used for this
!              purpose.  (Note: If this condition is detected before
!              taking any steps, then an illegal input return
!              (ISTATE = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  See written message for details.
!              Note:  If the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              The problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              This may be caused by a poor preconditioner matrix.
!          -6  means EWT(i) became zero for some i during the
!              integration.  Pure relative error control (ATOL(i)=0.0)
!              was requested on a variable which has now vanished.
!              The integration was successful as far as T.
!          -7  means an unrecoverable error occurred in JAC or PSOL.
!              Either JAC returned IER .ne. 0, or PSOL returned
!              IER .lt. 0.
!
!          Note:  Since the normal output value of ISTATE is 2,
!          it does not need to be reset for normal continuation.
!          Also, since a negative input value of ISTATE will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other input, before
!          calling the solver again.
!
! IOPT   = An integer flag to specify whether or not any optional
!          input is being used on this call.  Input only.
!          The optional input is listed separately below.
!          IOPT = 0 means no optional input is being used.
!                   Default values will be used in all cases.
!          IOPT = 1 means optional input is being used.
!
! RWORK  = A real working array (double precision).
!          The length of RWORK must be at least
!             20 + NYH*(MAXORD + 1) + 3*NEQ + LENK + LWP   where
!          NYH    = the initial value of NEQ,
!          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!                   smaller value is given as an optional input),
!          LENK = length of work space for Krylov-related data:
!          LENK = 0                                 if MITER = 0,
!          LENK = NEQ*(MAXL+3+MIN(1,MAXL-KMP))
!                  + (MAXL+3)*MAXL + 1              if MITER = 1,
!          LENK = 3*NEQ                             if MITER = 9.
!          LWP = length of real user work space for preconditioning.
!          (See JAC/PSOL.)
!          (See the MF description for METH and MITER.)
!          Thus if MAXORD etc. have default values and NEQ is constant,
!          this length is:
!             20 + 16*NEQ                    for MF = 10,
!             61 + 24*NEQ + LWP              for MF = 11,
!             20 + 19*NEQ + LWP              for MF = 19,
!             20 + 9*NEQ                     for MF = 20,
!             61 + 17*NEQ + LWP              for MF = 21,
!             20 + 12*NEQ + LWP              for MF = 29
!          The first 20 words of RWORK are reserved for conditional
!          and optional input and optional output.
!
!          The following word in RWORK is a conditional input:
!            RWORK(1) = TCRIT = critical value of t which the solver
!                       is not to overshoot.  Required if ITASK is
!                       4 or 5, and ignored otherwise.  (See ITASK.)
!
! LRW    = The length of the array RWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! IWORK  = An integer work array.  The length of IWORK must be at least
!             30        if MITER = 0  (MF = 10, 20), or
!             30 + LIWP  otherwise (MF = 11, 21, 19, 29).
!          LIWP = length of integer user work space for preconditioning.
!          (See conditional input list following).
!
!          The first 30 words of IWORK are reserved for conditional and
!          optional input and optional output.
!
!          The following 4 words in IWORK are conditional input,
!          required if MITER .ge. 1:
!
!          IWORK(1) = LWP  = length of real array WP for use in
!                     preconditioning (part of RWORK array).
!          IWORK(2) = LIWP = length of integer array IWP for use in
!                     preconditioning (part of IWORK array).
!                     The arrays WP and IWP are work arrays under the
!                     user's control, for use in the routines that
!                     perform preconditioning operations (JAC and PSOL).
!          IWORK(3) = JPRE = preconditioner type flag:
!                   = 0 for no preconditioning (P1 = P2 = I
!                   = 1 for left-only preconditioning (P2 = I)
!                   = 2 for right-only preconditioning (P1 = I)
!                   = 3 for two-sided preconditioning
!          IWORK(4) = JACFLG = flag for whether JAC is called.
!                   = 0 if JAC is not to be called,
!                   = 1 if JAC is to be called.
!                     Use JACFLG = 1 if JAC computes any nonconstant
!                     data needed in preconditioning operations,
!                     such as some of the Jacobian elements.
!
!
! LIW    = the length of the array IWORK, as declared by the user.
!          (This will be checked by the solver.)
!
! Note:  The work arrays must not be altered between calls to DVODPK
! for the same problem, except possibly for the conditional and
! optional input, and except for the last 3*NEQ words of RWORK.
! The latter space is used for internal scratch space, and so is
! available for use by the user outside DVODPK between calls, if
! desired (but not for use by F or JAC).
!
! JAC    = The name of the user-supplied routine (MITER = 1 or 9) to
!          compute the Jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y.  It is to have the form
!             SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY, V, HRL1,
!            1                WP, IWP, IER, RPAR, IPAR)
!             EXTERNAL F
!             DOUBLE PRECISION T, Y(NEQ), YSV(NEQ), REWT(NEQ), FTY(NEQ),
!            1                 V(NEQ), HRL1, WP(*), RPAR(*)
!             INTEGER IWP(*), IPAR(*)
!          This routine must evaluate and preprocess any parts of the
!          Jacobian matrix df/dy used in the preconditioners P1, P2 .
!          The Y and FTY arrays contain the current values of y and
!          f(t,y), respectively, and YSV also contains the current
!          value of y.  The array V is work space of length
!          NEQ for use by JAC.  REWT is the array of reciprocal error
!          weights (1/ewt).  JAC must multiply all computed Jacobian
!          elements by the scalar -hrl1, add the identity matrix I and
!          do any factorization operations called for, in preparation
!          for solving linear systems with a coefficient matrix of
!          P1 or P2.  The matrix P1*P2 should be an approximation to
!          I - hrl1 * (df/dy), where hrl1 is stored in HRL1.  JAC should
!          return IER = 0 if successful, and IER .ne. 0 if not.
!          (If IER .ne. 0, a smaller time step will be tried.)
!          The arrays WP (of length LWP) and IWP (of length LIWP)
!          are for use by JAC and PSOL for work space and for storage
!          of data needed for the solution of the preconditioner
!          linear systems.  Their lengths and contents are under the
!          user's control.
!          The JAC routine may save relevant Jacobian elements (or
!          approximations) used in the preconditioners, along with the
!          value of hrl1, and use these to reconstruct preconditioner
!          matrices later without reevaluationg those elements.
!          This may be cost-effective if JAC is called with hrl1
!          considerably different from its earlier value, indicating
!          that a corrector convergence failure has occurred because
!          of the change in hrl1, not because of changes in the
!          value of the Jacobian.  In doing this, use the saved and
!          current values of hrl1 to decide whether to use saved
!          or reevaluated elements.
!          JAC may alter V, but not Y, YSV, REWT, FTY, or HRL1.
!          JAC must be declared external in the calling program.
!
! PSOL   = the name of the user-supplied routine for the
!          solution of preconditioner linear systems.
!          It is to have the form
!             SUBROUTINE PSOL (NEQ, T, Y, FTY, WK, HRL1, WP, IWP, B, LR,
!            1                 IER, RPAR, IPAR)
!             DOUBLE PRECISION T, Y(NEQ), FTY(NEQ), WK(NEQ), HRL1,
!            1                 WP(*), B(NEQ), RPAR(*)
!             INTEGER  IWP(*), IPAR(*)
!          This routine must solve a linear system with b (stored in B)
!          as right-hand side and one of the preconditioning matrices,
!          P1 or P2, as coefficient matrix, and return the solution
!          vector in B.  LR is a flag concerning left vs. right
!          preconditioning, input to PSOL.  PSOL is to use P1 if LR = 1
!          and P2 if LR = 2.  In the case MITER = 9 (no Krylov
!          iteration), LR will be 1 and then 2, according to JPRE, and
!          PSOL is to return in B the desired approximate solution to 
!          A * x = b, where A = I - hrl1 * (df/dy).  (hrl1 is stored in
!          HRL1.)  PSOL can use data generated in the JAC routine and
!          stored in WP and IWP.  The Y and FTY arrays contain the 
!          current values of y and f(t,y), respectively.
!          The array WK is work space of length NEQ for use by PSOL.
!          The argument HRL1 is the current value of the scalar appear-
!          ing in the linear system.  If the old value, as of the last
!          JAC call, is needed, it must have been saved by JAC in WP.
!          On return, PSOL should set the error flag IER as follows:
!            IER = 0 if PSOL was successful,
!            IER .gt. 0 on a recoverable error, meaning that the
!                   time step will be retried,
!            IER .lt. 0 on an unrecoverable error, meaning that the
!                   solver is to stop immediately.
!          PSOL may not alter Y, FTY, or HRL1.
!          PSOL must be declared external in the calling program.
!
! MF     = The method flag.  Used only for input.  The legal values of
!          MF are 10, 11, 19, 20, 21, 29 .
!          MF is a two-digit integer, MF = 10*METH + MITER .
!          METH indicates the basic linear multistep method:
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on backward
!                     differentiation formulas (BDF-s).
!          MITER indicates the corrector iteration method.  Currently,
!            the values of MITER have the following meanings:
!
!          0 means functional iteration is used (no Jacobian matrix
!            is involved).
!
!          1 means SPIGMR, a scaled, preconditioned, incomplete version
!            of GMRES, a generalized minimum residual method, is used.
!            This is the best choice in general.
!
!          9 means that only a user-supplied matrix P (approximating A)
!            will be used, with no Krylov iteration done internally to
!            DVODPK.  This option allows the user to provide the
!            complete linear system solution algorithm, if desired.
!
! The user can apply preconditioning to the linear system A*x = b,
! by means of arbitrary matrices (the preconditioners).
!
! RPAR     User-specified array used to communicate real parameters
!          to user-supplied subroutines.  If RPAR is a vector, then
!          it must be dimensioned in the user's main program.  If it
!          is unused or a scalar, then it need not be dimensioned.
!
! IPAR     User-specified array used to communicate integer parameter
!          to user-supplied subroutines.  The comments on dimensioning
!          RPAR apply to IPAR.
!-----------------------------------------------------------------------
! Optional Inputs.
!
! The following is a list of the optional input provided for in the
! call sequence.  (See also Part ii.)  For each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! The use of any of this input requires IOPT = 1, and in that
! case all of this input is examined.  A value of zero for any
! of these optional input variables will cause the default value to be
! used.  Thus to use a subset of the optional input, simply preload
! locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! NAME    LOCATION      MEANING AND DEFAULT VALUE
!
! H0      RWORK(5)  The step size to be attempted on the first step.
!                   The default value is determined by the solver.
!
! HMAX    RWORK(6)  The maximum absolute step size allowed.
!                   The default value is infinite.
!
! HMIN    RWORK(7)  The minimum absolute step size allowed.
!                   The default value is 0.  (This lower bound is not
!                   enforced on the final step before reaching TCRIT
!                   when ITASK = 4 or 5.)
!
! DELT    RWORK(8)  Convergence test constant used in Krylov iteration
!                   algorithm.  The default value is 0.05.
!
! MAXORD  IWORK(5)  The maximum order to be allowed.  The default
!                   value is 12 if METH = 1, and 5 if METH = 2.
!                   If MAXORD exceeds the default value, it will
!                   be reduced to the default value.
!                   If MAXORD is changed during the problem, it may
!                   cause the current order to be reduced.
!
! MXSTEP  IWORK(6)  Maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   The default value is 500.
!
! MXHNIL  IWORK(7)  Maximum number of messages printed (per problem)
!                   warning that T + H = T on a step (H = step size).
!                   This must be positive to result in a non-default
!                   value.  The default value is 10.
!
! MAXL    IWORK(8)  maximum number of iterations in the SPIGMR
!                   algorithm (.le. NEQ).  The default is
!                   MAXL = min(5, NEQ).
!
! KMP     IWORK(9)  number of vectors on which orthogonalization
!                   is done in the SPIGMR algorithm (.le. MAXL).
!                   The default is KMP = MAXL (complete GMRES method).
!                   See Ref. 2 for details on incomplete GMRES.
!                   Note:  When KMP .lt. MAXL and MITER = 1, the length
!                   of RWORK must be set accordingly.  See RWORK above.
!-----------------------------------------------------------------------
! Optional Outputs.
!
! As optional additional output from DVODPK, the variables listed
! below are quantities related to the performance of DVODPK
! which are available to the user.  These are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! Except where stated otherwise, all of this output is defined
! on any successful return from DVODPK, and on any return with
! ISTATE = -1, -2, -4, -5, -6, or -7.  On an illegal input return
! (ISTATE = -3), they will be unchanged from their existing values
! (if any), except possibly for TOLSF, LENRW, and LENIW.
! On any error return, output relevant to the error will be defined,
! as noted below.
!
! NAME    LOCATION      MEANING
!
! HU      RWORK(11) The step size in t last used (successfully).
!
! HCUR    RWORK(12) The step size to be attempted on the next step.
!
! TCUR    RWORK(13) The current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  In the output,
!                   TCUR will always be at least as far from the
!                   initial value of t as the current argument T,
!                   but may be farther (if interpolation was done).
!
! TOLSF   RWORK(14) A tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise).  If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
!
! NST     IWORK(11) The number of steps taken for the problem so far.
!
! NFE     IWORK(12) The number of f evaluations for the problem so far.
!
! NPE     IWORK(13) The number of preconditioner evaluations (JAC calls)
!                   so far.
!
! NQU     IWORK(14) The method order last used (successfully).
!
! NQCUR   IWORK(15) The order to be attempted on the next step.
!
! IMXER   IWORK(16) The index of the component of largest magnitude in
!                   the weighted local error vector ( e(i)/EWT(i) ),
!                   on an error return with ISTATE = -4 or -5.
!
! LENRW   IWORK(17) The length of RWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! LENIW   IWORK(18) The length of IWORK actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! NNI     IWORK(20) The number of nonlinear iterations so far (each of
!                   which calls the Krylov iterative linear solver).
!
! NCFN    IWORK(21) The number of convergence failures of the nonlinear
!                   (Newton) iteration so far.
!                   Note: A measure of success is the overall rate of
!                   nonlinear convergence failures, NCFN/NST.
!
! NETF    IWORK(22) The number of error test failures of the integrator
!                   so far.
!
! NLI     IWORK(23) The number of linear iterations so far.
!                   Note: a measure of the success of SPIGMR algorithm
!                   is the average number of linear iterations per
!                   nonlinear iteration, given by NLI/NNI.
!                   If this is close to MAXL, MAXL may be too small.
!
! NPS     IWORK(24) The number of preconditioning solve operations
!                   (PSOL calls) so far.
!
! NCFL    IWORK(25) The number of convergence failures of the linear
!                   iteration so far.
!                   Note: A measure of success is the overall rate of
!                   linear convergence failures, NCFL/NNI.
!
! The following two arrays are segments of the RWORK array which
! may also be of interest to the user as optional output.
! For each array, the table below gives its internal name,
! its base address in RWORK, and its description.
!
! NAME    BASE ADDRESS      DESCRIPTION
!
! YH      21             The Nordsieck history array, of size NYH by
!                        (NQCUR + 1), where NYH is the initial value
!                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
!                        of YH contains HCUR**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the
!                        solution, evaluated at t = TCUR.
!
! ACOR     LENRW-NEQ+1   Array of size NEQ used for the accumulated
!                        corrections on each step, scaled in the output
!                        to represent the estimated local error in Y on
!                        the last step.  This is the vector e in the
!                        description of the error control.  Defined
!                        only on a successful return from DVODPK.
!
!-----------------------------------------------------------------------
! Interrupting and Restarting
!
! If the integration of a given problem by DVODPK is to be
! interrrupted and then later continued, such as when restarting
! an interrupted run or alternating between two or more ODE problems,
! the user should save, following the return from the last DVODPK call
! prior to the interruption, the contents of the call sequence
! variables and internal COMMON blocks, and later restore these
! values before the next DVODPK call for that problem.  To save
! and restore the COMMON blocks, use subroutine DVKSRC, as
! described below in Part ii.
!
! In addition, if non-default values for either LUN or MFLAG are
! desired, an extra call to XSETUN and/or XSETF should be made just
! before continuing the integration.  See Part ii below for details.
!
!-----------------------------------------------------------------------
! Part ii.  Other Routines Callable.
!
! The following are optional calls which the user may make to
! gain additional capabilities in conjunction with DVODPK.
! (The routines XSETUN and XSETF are designed to conform to the
! SLATEC error handling package.)
!
!     FORM OF CALL                  FUNCTION
!
!  CALL XSETUN(LUN)           Set the logical unit number, LUN, for
!                             output of messages from DVODPK, if
!                             the default is not desired.
!                             The default value of LUN is 6.
!
!  CALL XSETF(MFLAG)          Set a flag to control the printing of
!                             messages by DVODPK.
!                             MFLAG = 0 means do not print. (Danger:
!                             This risks losing valuable information.)
!                             MFLAG = 1 means print (the default).
!
!                             Either of the above calls may be made at
!                             any time and will take effect immediately.
!
!  CALL DVKSRC(RSAV,ISAV,JOB) Saves and restores the contents of
!                             the internal COMMON blocks used by
!                             DVODPK. (See Part iii below.)
!                             RSAV must be a real array of length 52
!                             or more, and ISAV must be an integer
!                             array of length 52 or more.
!                             JOB=1 means save COMMON into RSAV/ISAV.
!                             JOB=2 means restore COMMON from RSAV/ISAV.
!
!                                DVKSRC is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with DVODPK.
!
!  CALL DVINDY(,,,,,)         Provide derivatives of y, of various
!        (See below.)         orders, at a specified point T, if
!                             desired.  It may be called only after
!                             a successful return from DVODPK.
!
! The detailed instructions for using DVINDY are as follows.
! The form of the call is:
!
!  CALL DVINDY (T, K, RWORK(21), NYH, DKY, IFLAG)
!
! The input parameters are:
!
! T         = Value of independent variable where answers are desired
!             (normally the same as the T last returned by DVODPK).
!             For valid results, T must lie between TCUR - HU and TCUR.
!             (See optional output for TCUR and HU.)
! K         = Integer order of the derivative desired.  K must satisfy
!             0 .le. K .le. NQCUR, where NQCUR is the current order
!             (see optional output).  The capability corresponding
!             to K = 0, i.e. computing y(T), is already provided
!             by DVODPK directly.  Since NQCUR .ge. 1, the first
!             derivative dy/dt is always available with DVINDY.
! RWORK(21) = The base address of the history array YH.
! NYH       = Column length of YH, equal to the initial value of NEQ.
!
! The output parameters are:
!
! DKY       = A real array of length NEQ containing the computed value
!             of the K-th derivative of y(t).
! IFLAG     = Integer flag, returned as 0 if K and T were legal,
!             -1 if K was illegal, and -2 if T was illegal.
!             On an error return, a message is also written.
!-----------------------------------------------------------------------
! Part iii.  COMMON Blocks.
! If DVODPK is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in:
!   (1) the call sequence to DVODPK,
!   (2) the three internal COMMON blocks
!         /DVOD01/  of length  81  (48 double precision words
!                         followed by 33 integer words),
!         /DVOD02/  of length  9  (1 double precision word
!                         followed by 8 integer words),
!         /DVPK01/  of length 14 (3 double precision words
!                         followed by 11 integer words)
!
! If DVODPK is used on a system in which the contents of internal
! COMMON blocks are not preserved between calls, the user should
! declare the above three COMMON blocks in the calling (main) program
! to insure that their contents are preserved.
!
!-----------------------------------------------------------------------
! Part iv.  Optionally Replaceable Solver Routines.
!
! Below are descriptions of two routines in the DVODPK package which
! relate to the measurement of errors.  Either routine can be
! replaced by a user-supplied version, if desired.  However, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (Note: The means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
!
! (a) DEWSET.
! The following subroutine is called just before each internal
! integration step, and sets the array of error weights, EWT, as
! described under ITOL/RTOL/ATOL above:
!     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
! where NEQ, ITOL, RTOL, and ATOL are as in the DVODPK call sequence,
! YCUR contains the current dependent variable vector, and
! EWT is the array of weights set by DEWSET.
!
! If the user supplies this subroutine, it must return in EWT(i)
! (i = 1,...,NEQ) a positive quantity suitable for comparison with
! errors in Y(i).  The EWT array returned by DEWSET is passed to the
! DVNORM routine (see below), and also used by DVODPK in the computation
! of the optional output IMXER, the diagonal Jacobian approximation,
! and the increments for difference quotient Jacobians.
!
! In the user-supplied version of DEWSET, it may be desirable to use
! the current values of derivatives of y.  Derivatives up to order NQ
! are available from the history array YH, described above under
! Optional Output.  In DEWSET, YH is identical to the YCUR array,
! extended to NQ + 1 columns with a column length of NYH and scale
! factors of h**j/factorial(j).  On the first call for the problem,
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
! NYH is the initial value of NEQ.  The quantities NQ, H, and NST
! can be obtained by including in DEWSET the statements
!     COMMON /DVOD01/ RVOD(48), IVOD(33)
!     COMMON /DVOD02/ HU, NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
!     NQ = IVOD(28)
!     H = RVOD(21)
! Thus, for example, the current value of dy/dt can be obtained as
! YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
! unnecessary when NST = 0).
!
! (b) DVNORM.
! The following is a real function routine which computes the weighted
! root-mean-square norm of a vector v:
!     D = DVNORM (N, V, W)
! where:
!   N = the length of the vector,
!   V = real array of length N containing the vector,
!   W = real array of length N containing weights,
!   D = sqrt( (1/N) * sum(V(i)*W(i))**2 ).
! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
! EWT is as set by subroutine DEWSET.
!
! If the user supplies this function, it should return a non-negative
! value of DVNORM suitable for use in the error control in DVODPK.
! None of the arguments should be altered by DVNORM.
! For example, a user-supplied DVNORM routine might:
!   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
!   -ignore some components of V in the norm, with the effect of
!    suppressing the error control on those components of Y.
!-----------------------------------------------------------------------
!
! Revision History (YYYYMMDD)
! 19910315  DATE WRITTEN
! 19910415  Minor revisions to VODPK prologue.
! 19920715  In demo, corrected name R1MACH to D1MACH.
! 19921106  In VSTEP, added ETAQ and ETAQM1 to SAVE statement.
! 19930701  IN VNLSK, moved line setting HRL1 below statement 220.
! 19940502  Minor revisions to VODPK prologue and internal comments.
!           In VODPK, set JACFLG = 0 if 0 < MITER < 9 and JPRE = 0.
!           In VNLSK, add conditions on rescaling of correction vector.
! 19940504  In demo programs, fixed logic in SOLSBG involving LR.
! 19970515  Minor revisions to VODPK prologue and internal comments.
!           In VHIN, attached sign to H in second derivative estimation.
! 19981111  In VODPK, at end of Block B, when ISTATE = 3, jump to 200.
! 20020423  Major upgrade: Added *DECK lines.  Renamed all routines and
!           Common blocks for uniqueness across single/double prec.
!           versions and for sharing of routines with VODE and ODEPACK.
!           Changed names R1MACH/D1MACH to RUMACH/DUMACH.
!           Converted intrinsic names to generic form.
!           Numerous revisions to main prologue.
!           Revisions to demo program - formats, intrinsics, comments.
! 20020426  Converted upgraded single precision version to double prec.
!
!-----------------------------------------------------------------------
! Other Routines in the DVODPK Package.
!
! In addition to subroutine DVODPK, the DVODPK package includes the
! following subroutines and function routines:
!  DVHIN    computes an approximate step size for the initial step.
!  DVINDY   computes an interpolated value of the y vector at t = TOUT.
!  DVSTEP   is the core integrator, which does one step of the
!           integration and the associated error control.
!  DVSET    sets all method coefficients and test constants.
!  DVJUST   adjusts the history array on a change of order.
!  DVNLSK   solves the underlying nonlinear system -- the corrector.
!  DVSLPK   manages solution of linear system in chord iteration.
!  DVSPIG   performs the SPIGMR algorithm.
!  DVATV    computes a scaled, preconditioned product (I-hrl1*J)*v.
!  DORTHOG  orthogonalizes a vector against previous basis vectors.
!  DHEQR    generates a QR factorization of a Hessenberg matrix.
!  DHELS    finds the least squares solution of a Hessenberg system.
!  DVUSOL   interfaces to the user's PSOL routine (MITER = 9).
!  DEWSET   sets the error weight vector EWT before each step.
!  DVNORM   computes the weighted r.m.s. norm of a vector.
!  DVKSRC   is a user-callable routine to save and restore
!           the contents of the internal COMMON blocks.
!  DAXPY, DCOPY, DDOT, DNRM2, and DSCAL are basic linear
!           algebra modules (BLAS) used by this package.
!  DUMACH   computes the unit roundoff in a machine-independent manner.
!  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH handle the printing of all
!           error messages and warnings.  XERRWD is machine-dependent.
! Note:  DVNORM, DDOT, DNRM2, DUMACH, IXSAV, and IUMACH are function
! routines.  All the others are subroutines.
!
!-----------------------------------------------------------------------
!
! Declarations for external routines and function subroutines called ---
      EXTERNAL DVNLSK
      DOUBLE PRECISION DUMACH, DVNORM
!
! Declarations for local variables -------------------------------------
!
      LOGICAL IHIT, LAVD, LCFN, LCFL, LWARN
      DOUBLE PRECISION ATOLI, AVDIM, BIG, EWTI, FOUR, H0, HMAX, HMX,
     1   HUN, ONE, PT05, PT2, PT9, RCFL, RCFN, RH, RTOLI, SIZE,
     2   TCRIT, TNEXT, TOLSF, TP, TWO, ZERO
      INTEGER I, IER, IFLAG, IMXER, KGO, LENIW, LENIWK, LENRW, LENWK,
     1   LENWM, LF0, LIWP, LWP, MORD, MXHNL0, MXSTP0, NCFL0, NCFN0,
     2   NITER, NLI0, NNI0, NNID, NSTD, NSLAST, NWARN
      CHARACTER*80 MSG
      DIMENSION MORD(2)
!-----------------------------------------------------------------------
! The following Fortran-77 declarations are to cause the values of the
! listed (local) variables to be saved between calls to DVODPK.
!-----------------------------------------------------------------------
      SAVE MORD, MXHNL0, MXSTP0
      SAVE ZERO, ONE, TWO, FOUR, HUN, PT05, PT2, PT9
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
!
! Type declarations for labeled COMMON block DVPK01 --------------------
!
      DOUBLE PRECISION DELT, SQRTN, RSQRTN
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LVSAV, KMP, MAXL, MNEWT,
     1      NLI, NPS, NCFL
!
!-----------------------------------------------------------------------
! The following internal COMMON blocks contain variables which are
! communicated between subroutines in the DVODPK package, or which are
! to be saved between calls to DVODPK.
! In each block, real variables precede integers.
! The block /DVOD01/ appears in subroutines DVODPK, DVINDY, DVSTEP,
! DVSET, DVJUST, DVNLSK, DVSLPK, DVATV, and DVKSRC.
! The block /DVOD02/ appears in subroutines DVODPK, DVINDY, DVSTEP,
! DVNLSK, DVSLPK, DVATV, and DVKSRC.
! The block /DVPK01/ appears in subroutines DVODPK, DVNLSK, DVSLPK,
! and DVKSRC.
!
! The variables stored in the internal COMMON blocks are as follows:
!
! ACNRM  = Weighted r.m.s. norm of accumulated correction vectors.
! CCMXJ  = Threshhold on DRC for updating the Jacobian. (See DRC.)
! CONP   = The saved value of TQ(5).
! CRATE  = Estimated corrector convergence rate constant.
! DRC    = Relative change in H*RL1 since last VJAC call.
! EL     = Real array of integration coefficients.  See DVSET.
! ETA    = Saved tentative ratio of new to old H.
! ETAMAX = Saved maximum value of ETA to be allowed.
! H      = The step size.
! HMIN   = The minimum absolute value of the step size H to be used.
! HMXI   = Inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
! HNEW   = The step size to be attempted on the next step.
! HSCAL  = Stepsize in scaling of YH array.
! PRL1   = The saved value of RL1.
! RC     = Ratio of current H*RL1 to value on last VJAC call.
! RL1    = The reciprocal of the coefficient EL(1).
! TAU    = Real vector of past NQ step sizes, length 13.
! TQ     = A real vector of length 5 in which DVSET stores constants
!          used for the convergence test, the error test, and the
!          selection of H at a new order.
! TN     = The independent variable, updated on each step taken.
! UROUND = The machine unit roundoff.  The smallest positive real number
!          such that  1.0 + UROUND .ne. 1.0
! ICF    = Integer flag for convergence failure in DVNLSK:
!            0 means no failures.
!            1 means convergence failure with out of date Jacobian
!                   (recoverable error).
!            2 means convergence failure with current Jacobian or
!                   singular matrix (unrecoverable error).
! INIT   = Saved integer flag indicating whether initialization of the
!          problem has been done (INIT = 1) or not.
! IPUP   = Saved flag to signal updating of Newton matrix.
! JCUR   = Output flag from VJAC showing Jacobian status:
!            JCUR = 0 means J is not current.
!            JCUR = 1 means J is current.
! JSTART = Integer flag used as input to DVSTEP:
!            0  means perform the first step.
!            1  means take a new step continuing from the last.
!            -1 means take the next step with a new value of MAXORD,
!                  HMIN, HMXI, N, METH, MITER, and/or matrix parameters.
!          On return, DVSTEP sets JSTART = 1.
! JSV    = Integer flag for Jacobian saving, = sign(MF).
! KFLAG  = A completion code from DVSTEP with the following meanings:
!               0      the step was succesful.
!              -1      the requested error could not be achieved.
!              -2      corrector convergence could not be achieved.
!              -3, -4  fatal error in VNLS.
! KUTH   = Input flag to DVSTEP showing whether H was reduced by the
!          driver.  KUTH = 1 if H was reduced, = 0 otherwise.
! L      = Integer variable, NQ + 1, current order plus one.
! LMAX   = MAXORD + 1 (used for dimensioning).
! LOCJS  = A pointer to the saved Jacobian, whose storage starts at
!          WM(LOCJS), if JSV = 1.
! LYH, LEWT, LACOR, LSAVF, LWM, LIWM = Saved integer pointers
!          to segments of RWORK and IWORK.
! MAXORD = The maximum order of integration method to be allowed.
! METH/MITER = The method flags.  See MF.
! MSBJ   = The maximum number of steps between J evaluations, = 50.
! MXHNIL = Saved value of optional input MXHNIL.
! MXSTEP = Saved value of optional input MXSTEP.
! N      = The number of first-order ODEs, = NEQ.
! NEWH   = Saved integer to flag change of H.
! NEWQ   = The method order to be used on the next step.
! NHNIL  = Saved counter for occurrences of T + H = T.
! NQ     = Integer variable, the current integration method order.
! NQNYH  = Saved value of NQ*NYH.
! NQWAIT = A counter controlling the frequency of order changes.
!          An order change is about to be considered if NQWAIT = 1.
! NSLJ   = The number of steps taken as of the last Jacobian update.
! NSLP   = Saved value of NST as of last Newton matrix update.
! NYH    = Saved value of the initial value of NEQ.
!
! HU     = The step size in t last used.
! NCFN   = Number of nonlinear convergence failures so far.
! NETF   = The number of error test failures of the integrator so far.
! NFE    = The number of f evaluations for the problem so far.
! NPE    = The number of preconditioner evaluations (JAC calls) so far.
! NLU    = The number of matrix LU decompositions so far.
! NNI    = Number of nonlinear iterations so far.
! NQU    = The method order last used.
! NST    = The number of steps taken for the problem so far.
!
! DELT   = Convergence test constant in Krylov iterations.
! SQRTN  = SQRT(NEQ), for use in weights in Krylov convergence tests.
! RSQRTN = 1.0/SQRTN, also for use in convergence weights.
! JPRE   = Preconditioner type flag.
! JACFLG = Indicator for presence of user-supplied JAC routine.
! LOCWP  = Location of start of user's WP array in WM work array.
! LOCIWP = Location of start of user's IWP array in IWM work array.
! LVSAV  = Saved pointer to VSAV array in RWORK.
! KMP    = Number of vectors on which orthogonalization is done in
!          Krylov iteration.
! MAXL   = Maximum dimension of Krylov subspace used.
! MNEWT  = Newton iteration index.
! NLI    = Number of linear (Krylov) iterations done.
! NPS    = Number of preconditioner solvers (PSOL calls) done.
! NCFL   = Number of convergence failures in Krylov iteration.
!-----------------------------------------------------------------------
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
      COMMON /DVPK01/ DELT, SQRTN, RSQRTN, JPRE, JACFLG, LOCIWP,
     1                LOCWP, LVSAV, KMP, MAXL, MNEWT, NLI, NPS, NCFL
!
      DATA  MORD(1),MORD(2)/12,5/, MXSTP0/500/, MXHNL0/10/
      DATA ZERO/0.0D0/, ONE/1.0D0/, TWO/2.0D0/, FOUR/4.0D0/,
     1   PT05/0.05D0/, PT2/0.2D0/, PT9/0.9D0/, HUN/100.0D0/
!-----------------------------------------------------------------------
! Block A.
! This code block is executed on every call.
! It tests ISTATE and ITASK for legality and branches appropriately.
! If ISTATE .gt. 1 but the flag INIT shows that initialization has
! not yet been done, an error return occurs.
! If ISTATE = 1 and TOUT = T, jump to Block G and return immediately.
!-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
!-----------------------------------------------------------------------
! Block B.
! The next code block is executed for the initial call (ISTATE = 1),
! or for a continuation call with parameter changes (ISTATE = 3).
! It contains checking of all inputs and various initializations.
!
! First check legality of the non-optional inputs NEQ, ITOL, IOPT, MF.
!-----------------------------------------------------------------------
 20   IF (NEQ .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ .GT. N) GO TO 605
 25   N = NEQ
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      JSV = 0
      METH = MF/10
      MITER = MF - 10*METH
      IF (METH .LT. 1 .OR. METH .GT. 2) GO TO 608
      IF (MITER .LT. 0) GO TO 608
      IF (MITER .GT. 1 .AND. MITER .LT. 9) GO TO 608
      IF (MITER .GE. 1) JPRE = IWORK(3)
      JACFLG = 0
      IF (MITER .GE. 1) JACFLG = IWORK(4)
      IF (MITER .GE. 1 .AND. MITER .NE. 9 .AND. JPRE .EQ. 0) JACFLG = 0
! Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      MAXORD = MORD(METH)
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      IF (ISTATE .EQ. 1) H0 = ZERO
      HMXI = ZERO
      HMIN = ZERO
      MAXL = MIN(5,N)
      KMP = MAXL
      DELT = PT05
      GO TO 60
 40   MAXORD = IWORK(5)
      IF (MAXORD .LT. 0) GO TO 611
      IF (MAXORD .EQ. 0) MAXORD = 100
      MAXORD = MIN(MAXORD,MORD(METH))
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      IF ((TOUT - T)*H0 .LT. ZERO) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. ZERO) GO TO 615
      HMXI = ZERO
      IF (HMAX .GT. ZERO) HMXI = ONE/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. ZERO) GO TO 616
      MAXL = IWORK(8)
      IF (MAXL .EQ. 0) MAXL = 5
      MAXL = MIN(MAXL,N)
      KMP = IWORK(9)
      IF (KMP .EQ. 0 .OR. KMP .GT. MAXL) KMP = MAXL
      DELT = RWORK(8)
      IF (DELT .EQ. 0.0D0) DELT = PT05
!-----------------------------------------------------------------------
! Set work array pointers and check lengths lrw and liw.
! Pointers to segments of RWORK and iwork are named by prefixing l to
! the name of the segment.  e.g., the segment YH starts at RWORK(LYH).
! Segments of RWORK (in order) are  YH, WM, EWT, SAVF, VSAV, ACOR.
! Within WM, LOCWP is the location of the WP work array,
! and within IWM, LOCIWP is the location of the IWP work array.
!-----------------------------------------------------------------------
 60   LYH = 21
      IF (ISTATE .EQ. 1) NYH = N
      LWM = LYH + (MAXORD + 1)*NYH
      IF (MITER .EQ. 0) LENWK = 0
      IF (MITER .EQ. 1)
     1   LENWK = N*(MAXL+2+MIN(1,MAXL-KMP)) + (MAXL+3)*MAXL + 1
      IF (MITER .EQ. 9) LENWK = 2*N
      LWP = 0
      IF (MITER .GE. 1) LWP = IWORK(1)
      LENWM = LENWK + LWP
      LOCWP = LENWK + 1
      LEWT = LWM + LENWM
      LSAVF = LEWT + N
      LVSAV = LSAVF + N
      LACOR = LVSAV + N
      IF (MITER .EQ. 0) LACOR = LVSAV
      LENRW = LACOR + N - 1
      IWORK(17) = LENRW
      LIWM = 31
      LENIWK = 0
      LIWP = 0
      IF (MITER .GE. 1) LIWP = IWORK(2)
      LENIW = 30 + LENIWK + LIWP
      LOCIWP = LENIWK + 1
      IWORK(18) = LENIW
      IF (LENRW .GT. LRW) GO TO 617
      IF (LENIW .GT. LIW) GO TO 618
! Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 70 I = 1, N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. ZERO) GO TO 619
        IF (ATOLI .LT. ZERO) GO TO 620
 70     CONTINUE
! Load SQRT(N) and its reciprocal in common. ---------------------------
      SQRTN = SQRT(DBLE(N))
      RSQRTN = ONE/SQRTN
      IF (ISTATE .EQ. 1) GO TO 100
! If ISTATE = 3, set flag to signal parameter changes to DVSTEP. -------
      JSTART = -1
      IF (NQ .LE. MAXORD) GO TO 200
! MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. ---------
      CALL DCOPY (N, RWORK(LWM), 1, RWORK(LSAVF), 1)
      GO TO 200
!-----------------------------------------------------------------------
! Block C.
! The next block is for the initial call only (ISTATE = 1).
! It contains all remaining initializations, the initial call to F,
! and the calculation of the initial step size.
! The error weights in EWT are inverted after being loaded.
!-----------------------------------------------------------------------
 100  UROUND = DUMACH()
      TN = T
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. ZERO) GO TO 625
      IF (H0 .NE. ZERO .AND. (T + H0 - TCRIT)*H0 .GT. ZERO)
     1   H0 = TCRIT - T
 110  JSTART = 0
      CCMXJ = PT2
      MSBJ = 50
      NHNIL = 0
      NST = 0
      NSLAST = 0
      HU = ZERO
      NQU = 0
      NPE = 0
      NLI0 = 0
      NNI0 = 0
      NCFN0 = 0
      NCFL0 = 0
      NWARN = 0
      NNI = 0
      NLI = 0
      NPS = 0
      NETF = 0
      NCFN = 0
      NCFL = 0
! Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (N, T, Y, RWORK(LF0), RPAR, IPAR)
      NFE = 1
! Load the initial value vector in YH. ---------------------------------
      CALL DCOPY (N, Y, 1, RWORK(LYH), 1)
! Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = ONE
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1, N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 621
 120    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
      IF (H0 .NE. ZERO) GO TO 180
! Call DVHIN to set initial step size H0 to be attempted. --------------
      CALL DVHIN (N, T, RWORK(LYH), RWORK(LF0), F, RPAR, IPAR, TOUT,
     1   UROUND, RWORK(LEWT), ITOL, ATOL, Y, RWORK(LACOR), H0,
     2   NITER, IER)
      NFE = NFE + NITER
      IF (IER .NE. 0) GO TO 622
! Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. ONE) H0 = H0/RH
! Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      CALL DSCAL (N, H0, RWORK(LF0), 1)
      GO TO 270
!-----------------------------------------------------------------------
! Block D.
! The next code block is for continuation calls only (ISTATE = 2 or 3)
! and is to check stop conditions before taking a step.
!-----------------------------------------------------------------------
 200  NSLAST = NST
      KUTH = 0
      NLI0 = NLI
      NNI0 = NNI
      NCFN0 = NCFN
      NCFL0 = NCFL
      NWARN = 0
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(ONE + HUN*UROUND)
      IF ((TP - TOUT)*H .GT. ZERO) GO TO 623
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. ZERO) GO TO 625
      IF ((TN - TOUT)*H .LT. ZERO) GO TO 245
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. ZERO) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + HNEW*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
!-----------------------------------------------------------------------
! Block E.
! The next block is normally executed for all calls and contains
! the call to the one-step core integrator DVSTEP.
!
! This is a looping point for the integration steps.
!
! First check for too many steps being taken,
! check for poor Newton/Krylov performance, update EWT (if not at
! start of problem), check for too much accuracy being requested, and
! check for H below the roundoff level in T.
!-----------------------------------------------------------------------
 250  CONTINUE
      IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      NSTD = NST - NSLAST
      NNID = NNI - NNI0
      IF (NSTD .LT. 10 .OR. NNID .EQ. 0) GO TO 255
      AVDIM = REAL(NLI - NLI0)/REAL(NNID)
      RCFN = REAL(NCFN - NCFN0)/REAL(NSTD)
      RCFL = REAL(NCFL - NCFL0)/REAL(NNID)
      LAVD = AVDIM .GT. (MAXL - PT05)
      LCFN = RCFN .GT. PT9
      LCFL = RCFL .GT. PT9
      LWARN = LAVD .OR. LCFN .OR. LCFL
      IF (.NOT.LWARN) GO TO 255
      NWARN = NWARN + 1
      IF (NWARN .GT. 10) GO TO 255
      IF (LAVD) THEN
        MSG = 'DVODPK- Warning. Poor iterative algorithm performance   '
        CALL XERRWD (MSG, 56, 111, 0, 0, 0, 0, 0, ZERO, ZERO)
        MSG = '      at T = R1. Average no. of linear iterations = R2  '
        CALL XERRWD (MSG, 56, 111, 0, 0, 0, 0, 2, TN, AVDIM)
        ENDIF
      IF (LCFN) THEN
        MSG = 'DVODPK- Warning. Poor iterative algorithm performance   '
        CALL XERRWD (MSG, 56, 112, 0, 0, 0, 0, 0, ZERO, ZERO)
        MSG = '      at T = R1. Nonlinear convergence failure rate = R2'
        CALL XERRWD (MSG, 56, 112, 0, 0, 0, 0, 2, TN, RCFN)
        ENDIF
      IF (LCFL) THEN
        MSG = 'DVODPK- Warning. Poor iterative algorithm performance   '
        CALL XERRWD (MSG, 56, 113, 0, 0, 0, 0, 0, ZERO, ZERO)
        MSG = '      at T = R1. Linear convergence failure rate = R2   '
        CALL XERRWD (MSG, 56, 113, 0, 0, 0, 0, 2, TN, RCFL)
        ENDIF
 255  CONTINUE
      DO 260 I = 1, N
        IF (RWORK(I+LEWT-1) .LE. ZERO) GO TO 510
 260    RWORK(I+LEWT-1) = ONE/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*DVNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. ONE) GO TO 280
      TOLSF = TOLSF*TWO
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'DVODPK-  Warning: internal T (=R1) and H (=R2) are'
      CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      (H = step size). solver will continue anyway'
      CALL XERRWD (MSG, 50, 101, 1, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'DVODPK-  Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 50, 102, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      it will not be issued again for this problem'
      CALL XERRWD (MSG, 50, 102, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
 290  CONTINUE
!-----------------------------------------------------------------------
!  CALL DVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR, WM, IWM,
!                                     F, JAC, PSOL, DVNLSK, RPAR, IPAR)
!-----------------------------------------------------------------------
      CALL DVSTEP (Y, RWORK(LYH), NYH, RWORK(LYH), RWORK(LEWT),
     1   RWORK(LSAVF), RWORK(LVSAV), RWORK(LACOR), RWORK(LWM),
     2   IWORK(LIWM), F, JAC, PSOL, DVNLSK, RPAR, IPAR)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540, 550, 555), KGO
!-----------------------------------------------------------------------
! Block F.
! The following block handles the case of a successful return from the
! core integrator (KFLAG = 0).  Test for stop conditions.
!-----------------------------------------------------------------------
 300  INIT = 1
      KUTH = 0
      GO TO (310, 400, 330, 340, 350), ITASK
! ITASK = 1.  if TOUT has been reached, interpolate. -------------------
 310  IF ((TN - TOUT)*H .LT. ZERO) GO TO 250
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
! ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. ZERO) GO TO 400
      GO TO 250
! ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. ZERO) GO TO 345
      CALL DVINDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(ONE + FOUR*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. ZERO) GO TO 250
      H = (TCRIT - TN)*(ONE - FOUR*UROUND)
      KUTH = 1
      GO TO 250
! ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. HUN*UROUND*HMX
!-----------------------------------------------------------------------
! Block G.
! The following block handles all successful returns from DVODPK.
! If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
! ISTATE is set to 2, and the optional outputs are loaded into the
! work arrays before returning.
!-----------------------------------------------------------------------
 400  CONTINUE
      CALL DCOPY (N, RWORK(LYH), 1, Y, 1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NPE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      IWORK(23) = NLI
      IWORK(24) = NPS
      IWORK(25) = NCFL
      RETURN
!-----------------------------------------------------------------------
! Block H.
! The following block handles all unsuccessful returns other than
! those for illegal input.  First the error message routine is called.
! if there was an error test or convergence test failure, IMXER is set.
! Then Y is loaded from YH, and T is set to TN.  The optional outputs
! are loaded into the work arrays before returning.
!-----------------------------------------------------------------------
! The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'DVODPK-  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 50, 201, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 1, 1, MXSTEP, 0, 1, TN, ZERO)
      ISTATE = -1
      GO TO 580
! EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'DVODPK-  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWD (MSG, 50, 202, 1, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
! Too much accuracy requested for machine precision. -------------------
 520  MSG = 'DVODPK-  At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      for precision of machine:  See TOLSF (=R2)  '
      CALL XERRWD (MSG, 50, 203, 1, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
! KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'DVODPK-  At T(=R1) and step size H(=R2), the error'
      CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      test failed repeatedly or with abs(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
! KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'DVODPK-  At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      or with abs(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 1, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 560
! KFLAG = -3.  Unrecoverable error from JAC. ---------------------------
 550  MSG = 'DVODPK-  at T (=R1) an unrecoverable error return '
      CALL XERRWD(MSG, 50, 206, 0, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      was made from subroutine JAC      '
      CALL XERRWD(MSG, 40, 206, 0, 0, 0, 0, 1, TN, ZERO)
      ISTATE = -7
      GO TO 580
! KFLAG = -4.  Unrecoverable error from PSOL. --------------------------
 555  MSG = 'DVODPK-  at T (=R1) an unrecoverable error return '
      CALL XERRWD(MSG, 50, 207, 0, 0, 0, 0, 0, ZERO, ZERO)
      MSG = '      was made from subroutine PSOL     '
      CALL XERRWD(MSG, 40, 207, 0, 0, 0, 0, 1, TN, ZERO)
      ISTATE = -7
      GO TO 580
! Compute IMXER if relevant. -------------------------------------------
 560  BIG = ZERO
      IMXER = 1
      DO 570 I = 1, N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
! Set Y vector, T, and optional outputs. -------------------------------
 580  CONTINUE
      CALL DCOPY (N, RWORK(LYH), 1, Y, 1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NPE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(20) = NNI
      IWORK(21) = NCFN
      IWORK(22) = NETF
      IWORK(23) = NLI
      IWORK(24) = NPS
      IWORK(25) = NCFL
      RETURN
!-----------------------------------------------------------------------
! Block I.
! The following block handles all error returns due to illegal input
! (ISTATE = -3), as detected before calling the core integrator.
! Call the error message routine and then return.
!-----------------------------------------------------------------------
 601  MSG = 'DVODPK-  ISTATE (=I1) illegal '
      CALL XERRWD (MSG, 30, 1, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'DVODPK-  ITASK (=I1) illegal  '
      CALL XERRWD (MSG, 30, 2, 1, 1, ITASK, 0, 0, ZERO, ZERO)
      GO TO 700
 603  MSG='DVODPK-   ISTATE (=I1) .gt. 1 but DVODPK not initialized    '
      CALL XERRWD (MSG, 60, 3, 1, 1, ISTATE, 0, 0, ZERO, ZERO)
      GO TO 700
 604  MSG = 'DVODPK-  NEQ (=I1) .lt. 1     '
      CALL XERRWD (MSG, 30, 4, 1, 1, NEQ, 0, 0, ZERO, ZERO)
      GO TO 700
 605  MSG = 'DVODPK-  ISTATE = 3 and NEQ increased (I1 to I2)  '
      CALL XERRWD (MSG, 50, 5, 1, 2, N, NEQ, 0, ZERO, ZERO)
      GO TO 700
 606  MSG = 'DVODPK-  ITOL (=I1) illegal   '
      CALL XERRWD (MSG, 30, 6, 1, 1, ITOL, 0, 0, ZERO, ZERO)
      GO TO 700
 607  MSG = 'DVODPK-  IOPT (=I1) illegal   '
      CALL XERRWD (MSG, 30, 7, 1, 1, IOPT, 0, 0, ZERO, ZERO)
      GO TO 700
 608  MSG = 'DVODPK-  MF (=I1) illegal     '
      CALL XERRWD (MSG, 30, 8, 1, 1, MF, 0, 0, ZERO, ZERO)
      GO TO 700
 611  MSG = 'DVODPK-  MAXORD (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 11, 1, 1, MAXORD, 0, 0, ZERO, ZERO)
      GO TO 700
 612  MSG = 'DVODPK-  MXSTEP (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 12, 1, 1, MXSTEP, 0, 0, ZERO, ZERO)
      GO TO 700
 613  MSG = 'DVODPK-  MXHNIL (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 13, 1, 1, MXHNIL, 0, 0, ZERO, ZERO)
      GO TO 700
 614  MSG = 'DVODPK-  TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 40, 14, 1, 0, 0, 0, 2, TOUT, T)
      MSG = '      integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 1, 0, 0, 0, 1, H0, ZERO)
      GO TO 700
 615  MSG = 'DVODPK-  HMAX (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 15, 1, 0, 0, 0, 1, HMAX, ZERO)
      GO TO 700
 616  MSG = 'DVODPK-  HMIN (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 16, 1, 0, 0, 0, 1, HMIN, ZERO)
      GO TO 700
 617  CONTINUE
      MSG='DVODPK-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 17, 1, 2, LENRW, LRW, 0, ZERO, ZERO)
      GO TO 700
 618  CONTINUE
      MSG='DVODPK-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWD (MSG, 60, 18, 1, 2, LENIW, LIW, 0, ZERO, ZERO)
      GO TO 700
 619  MSG = 'DVODPK-  RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 19, 1, 1, I, 0, 1, RTOLI, ZERO)
      GO TO 700
 620  MSG = 'DVODPK-  ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 20, 1, 1, I, 0, 1, ATOLI, ZERO)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'DVODPK-  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWD (MSG, 40, 21, 1, 1, I, 0, 1, EWTI, ZERO)
      GO TO 700
 622  CONTINUE
      MSG='DVODPK-  TOUT (=R1) too close to T(=R2) to start integration'
      CALL XERRWD (MSG, 60, 22, 1, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  CONTINUE
      MSG='DVODPK-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 60, 23, 1, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  CONTINUE
      MSG='DVODPK-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD (MSG, 60, 24, 1, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  CONTINUE
      MSG='DVODPK-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD (MSG, 60, 25, 1, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'DVODPK-  At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 50, 26, 1, 0, 0, 0, 0, ZERO, ZERO)
      MSG='      requested for precision of machine:  See TOLSF (=R1)  '
      CALL XERRWD (MSG, 60, 26, 1, 0, 0, 0, 1, TOLSF, ZERO)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG='DVODPK-  Trouble from DVINDY. ITASK = I1, TOUT = R1         '
      CALL XERRWD (MSG, 60, 27, 1, 1, ITASK, 0, 1, TOUT, ZERO)
!
 700  CONTINUE
      ISTATE = -3
      RETURN
!
 800  MSG = 'DVODPK-  Run aborted: apparent infinite loop      '
      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, ZERO, ZERO)
      RETURN
!----------------------- End of Subroutine DVODPK ----------------------
      END
*DECK DVNLSK
      SUBROUTINE DVNLSK (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM,
     1                   F, JAC, PSOL, NFLAG, RPAR, IPAR)
!
      EXTERNAL F, JAC, PSOL
      DOUBLE PRECISION Y, YH, VSAV, SAVF, EWT, ACOR, WM, RPAR
      INTEGER IWM, LDYH, NFLAG, IPAR
      DIMENSION Y(*), YH(LDYH, *),  SAVF(*), VSAV(*), EWT(*), ACOR(*),
     1          IWM(*), WM(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- YH, LDYH, F, JAC, EWT, ACOR, IWM, WM,
!                        NFLAG, RPAR, IPAR
! Call sequence output -- Y, YH, VSAV, SAVF, ACOR, IWM, WM, NFLAG
! COMMON block variables accessed:
!        /DVOD01/  ACNRM, CRATE, DRC, H, ICF, IPUP, JCUR, JSTART,
!                  METH, MITER, N, NSLP, RC, RL1, TN, TQ
!        /DVOD02/  NFE, NNI, NPE, NST
!        /DVPK01/  JACFLG, LOCIWP, LOCWP, MNEWT
! Subroutines called: F, JAC, PSOL, DAXPY, DCOPY, DSCAL, DVSLPK
! Function subroutines called: DVNORM
!-----------------------------------------------------------------------
! Subroutine DVNLSK is a nonlinear system solver, which uses either
! functional iteration (MITER = 0), or a combination of an inexact
! Newton method and preconditioned Krylov iteration (MITER .gt. 0)
! to solve the implicit system for the corrector y vector.
! It calls Subroutine JAC (user-supplied) for preprocessing the
! preconditioner, and Subroutine DVSLPK for the Krylov iteration.
!
! In addition to variables described elsewhere, communication with
! DVNLSK uses the following variables:
!
! Y          = The dependent variable, a vector of length N, input.
! YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
!              and output.  On input, it contains predicted values.
! LDYH       = A constant .ge. N, the first dimension of YH, input.
! VSAV       = A work array of length N.
! SAVF       = A work array of length N.
! EWT        = An error weight vector of length N, input.
! ACOR       = A work array of length N, used for the accumulated
!              corrections to the predicted y vector.
! WM,IWM     = Real and integer work arrays associated with matrix
!              operations in Newton iteration (MITER .ne. 0).
! F          = Dummy name for user-supplied routine for f.
! JAC        = Dummy name for user-supplied routine for Jacobian data
!              and associated preconditioner matrix.
! PSOL       = Dummy name for user-supplied subroutine to solve
!              preconditioner linear system.
! NFLAG      = Input/output flag, with values and meanings as follows:
!              INPUT
!                  0 first call for this time step.
!                 -1 convergence failure in previous call to DVNLSK.
!                 -2 error test failure in DVSTEP.
!              OUTPUT
!                  0 successful completion of nonlinear solver.
!                 -1 convergence failure or failure in JAC.
!                 -2 unrecoverable error in matrix preprocessing
!                    (cannot occur here).
!                 -3 unrecoverable error in PSOL.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
!
! IPUP       = Own variable flag with values and meanings as follows:
!              0,            do not update preconditioner.
!              MITER .ne. 0, update the preconditioner, because it is
!                            the initial step, user input changed,
!                            there was an error test failure, or an
!                            update is indicated by a change in the
!                            scalar RC or step counter NST.
!
! For more details, see comments in driver subroutine.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
!
! Type declarations for labeled COMMON block DVPK01 --------------------
!
      DOUBLE PRECISION DELT, SQRTN, RSQRTN
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LVSAV, KMP, MAXL, MNEWT,
     1      NLI, NPS, NCFL
!
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
      COMMON /DVPK01/ DELT, SQRTN, RSQRTN, JPRE, JACFLG, LOCIWP,
     1                LOCWP, LVSAV, KMP, MAXL, MNEWT, NLI, NPS, NCFL
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION CCMAX, CRDOWN, CSCALE, DEL, DCON, DELP, HRL1,RDIV
      DOUBLE PRECISION ONE, TWO, ZERO
      INTEGER I, IERPJ, IERSL, M, MAXCOR, MSBP
!
! Type declaration for function subroutines called ---------------------
!
      DOUBLE PRECISION DVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declarations are to cause the values of the
! listed (local) variables to be saved between calls to DVODPK.
      SAVE CCMAX, CRDOWN, MAXCOR, MSBP, RDIV
      SAVE ONE, TWO, ZERO
!-----------------------------------------------------------------------
      DATA CCMAX /0.3D0/, CRDOWN /0.3D0/, MAXCOR /3/, MSBP /20/,
     1     RDIV   /2.0D0/
      DATA ONE /1.0D0/, TWO /2.0D0/, ZERO /0.0D0/
!-----------------------------------------------------------------------
! Up to MAXCOR corrector iterations are taken.  A convergence test is
! made on the RMS norm of each correction, weighted by the error
! weight vector EWT.  The sum of the corrections is accumulated in the
! vector ACOR(*).  The YH array is not altered in the corrector loop.
!-----------------------------------------------------------------------
      IF (JSTART .EQ. 0) NSLP = 0
      IF (NFLAG .EQ. 0) ICF = 0
      IF (NFLAG .EQ. -2) IPUP = MITER
      IF ( (JSTART .EQ. 0) .OR. (JSTART .EQ. -1) ) IPUP = MITER
      IF (JACFLG .EQ. 0) THEN
        IPUP = 0
        CRATE = ONE
        GO TO 220
      ENDIF
      DRC = ABS(RC-ONE)
      IF (DRC .GT. CCMAX .OR. NST .GE. NSLP+MSBP) IPUP = MITER
 220  M = 0
      HRL1 = H*RL1
      DELP = ZERO
      MNEWT = 0
      CALL DCOPY (N, YH(1,1), 1, Y, 1 )
      CALL F (N, TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
!-----------------------------------------------------------------------
! If indicated, the preconditioner matrix is reevaluated and
! preprocessed before starting the corrector iteration.  IPUP is set
! to 0 as an indicator that this has been done.
!-----------------------------------------------------------------------
      JCUR = 1
      IERPJ = 0
      CALL JAC (F, N, TN, Y, YH, EWT, SAVF, ACOR, HRL1,
     1   WM(LOCWP), IWM(LOCIWP), IERPJ, RPAR, IPAR)
      NPE = NPE + 1
      IPUP = 0
      RC = ONE
      DRC = ZERO
      CRATE = ONE
      NSLP = NST
      IF (IERPJ .NE. 0) GO TO 420
 250  DO 260 I = 1, N
 260    ACOR(I) = ZERO
 270  IF (MITER .NE. 0) GO TO 350
!-----------------------------------------------------------------------
! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.
!-----------------------------------------------------------------------
      DO 290 I = 1, N
        SAVF(I) = RL1*(H*SAVF(I) - YH(I,2))
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = DVNORM (N, Y, EWT)
      DO 300 I = 1, N
 300    Y(I) = YH(I,1) + SAVF(I)
      CALL DCOPY (N, SAVF, 1, ACOR, 1)
      GO TO 400
!-----------------------------------------------------------------------
! In the case of the Newton method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! A as coefficient matrix.  In the case of Modified Newton iteration
! with BDF, the correction is scaled by the factor 2/(1+RC) to
! account for changes in H*RL1 since the last JAC call.
!-----------------------------------------------------------------------
 350  DO 360 I = 1, N
 360    VSAV(I) = HRL1*SAVF(I) - (RL1*YH(I,2) + ACOR(I))
      CALL DVSLPK (Y, SAVF, VSAV, EWT, WM, IWM, F, PSOL, IERSL,
     1             RPAR, IPAR)
      NNI = NNI + 1
      IF (METH .EQ. 2 .AND. JACFLG .EQ. 1 
     1    .AND. MITER .EQ. 9 .AND. RC .NE. ONE) THEN
        CSCALE = TWO/(ONE + RC)
        CALL DSCAL (N, CSCALE, VSAV, 1)
      ENDIF
      IF (IERSL .LT. 0) GO TO 440
      IF (IERSL .GT. 0) GO TO 410
      DEL = DVNORM (N, VSAV, EWT)
      CALL DAXPY (N, ONE, VSAV, 1, ACOR, 1)
      DO 380 I = 1, N
 380    Y(I) = YH(I,1) + ACOR(I)
!-----------------------------------------------------------------------
! Test for convergence.  If M.gt.0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.
!-----------------------------------------------------------------------
 400  IF (M .NE. 0) CRATE = MAX(CRDOWN*CRATE,DEL/DELP)
      DCON = DEL*MIN(ONE,CRATE)/TQ(4)
      IF (DCON .LE. ONE) GO TO 450
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. RDIV*DELP) GO TO 410
      MNEWT = M
      DELP = DEL
      CALL F (N, TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      GO TO 270
!
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1 .OR. JACFLG .EQ. 0) GO TO 420
      ICF = 1
      IPUP = MITER
      GO TO 220
!
 420  CONTINUE
      ICF = 2
      NFLAG = -1
      RETURN
 440  CONTINUE
      NFLAG = -3
      RETURN
! Return for successful step. ------------------------------------------
 450  NFLAG = 0
      JCUR = 0
      ICF = 0
      IF (M .EQ. 0) ACNRM = DEL
      IF (M .GT. 0) ACNRM = DVNORM (N, ACOR, EWT)
      RETURN
!----------------------- End of Subroutine DVNLSK ----------------------
      END
*DECK DVSLPK
      SUBROUTINE DVSLPK (Y, SAVF, X, EWT, WM, IWM, F, PSOL, IERSL,
     1                   RPAR, IPAR)
      EXTERNAL F, PSOL
      DOUBLE PRECISION Y, SAVF, X, EWT, WM, RPAR
      INTEGER IWM, IERSL, IPAR
      DIMENSION Y(*), SAVF(*), X(*), EWT(*), WM(*), IWM(*),
     1   RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, SAVF, X, EWT, F, PSOL, RPAR, IPAR
! Call sequence output -- Y, SAVF, X, WM, IWM, IERSL
! COMMON block variables accessed:
!        /DVOD01/  H, RL1, TQ, TN, MITER, N
!        /DVPK01/  DELT, SQRTN, RSQRTN, JPRE, LOCIWP, LOCWP,
!                  KMP, MAXL, MNEWT, NLI, NPS, NCFL
! Subroutines called: F, PSOL, DCOPY, DSCAL, DVSPIG, DVUSOL
!-----------------------------------------------------------------------
! This routine interfaces with  DVSPIG  or  DVUSOL  for the solution of
! the linear system arising from a Newton iteration (MITER .ne. 0).
!
! In addition to variables described elsewhere, communication with
! DVSLPK uses the following variables:
! WM    = real work space containing data for the algorithm
!         (Krylov basis vectors, Hessenberg matrix, etc.)
! IWM   = integer work space containing data for the algorithm
! X     = the right-hand side vector on input, and the solution vector
!         on output, of length N.
! IERSL = output flag (in COMMON):
!         IERSL =  0 means no trouble occurred.
!         IERSL =  1 means the iterative method failed to converge.
!                    If the preconditioner is out of date, the step
!                    is repeated with a new preconditioner.  Otherwise,
!                    the stepsize is reduced (forcing a new evalua-
!                    tion of the preconditioner) and the step is
!                    repeated.
!         IERSL = -1 means there was a nonrecoverable error in the
!                    iterative solver.  The stepsize is reduced in
!                    DVSTEP and the step is repeated.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
!
! Type declarations for labeled COMMON block DVPK01 --------------------
!
      DOUBLE PRECISION DELT, SQRTN, RSQRTN
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LVSAV, KMP, MAXL, MNEWT,
     1      NLI, NPS, NCFL
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION DELTA, HRL1
      INTEGER IFLAG, LB, LDL, LGMR, LHES, LQ, LV, LWK, MAXLP1, NPSL
!-----------------------------------------------------------------------
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVPK01/ DELT, SQRTN, RSQRTN, JPRE, JACFLG, LOCIWP,
     1                LOCWP, LVSAV, KMP, MAXL, MNEWT, NLI, NPS, NCFL
!-----------------------------------------------------------------------
!
      IERSL = 0
      HRL1 = H*RL1
      DELTA = DELT*TQ(4)
      IF (MITER .EQ. 1) THEN
!-----------------------------------------------------------------------
! Use the SPIGMR algorithm to solve the linear system A*x = -f.
!-----------------------------------------------------------------------
        MAXLP1 = MAXL + 1
        LV = 1
        LB = LV + N*MAXL
        LHES = LB + N + 1
        LQ = LHES + MAXL*MAXLP1
        LWK = LQ + 2*MAXL
        LDL = LWK + MIN(1,MAXL-KMP)*N
        CALL DCOPY (N, X, 1, WM(LB), 1)
        CALL DSCAL (N, RSQRTN, EWT, 1)
        CALL DVSPIG (TN, Y, SAVF, WM(LB), EWT, N, MAXL, MAXLP1, KMP,
     1     DELTA, HRL1, JPRE, MNEWT, F, PSOL, NPSL, X, WM(LV), WM(LHES),
     2     WM(LQ), LGMR, WM(LOCWP), IWM(LOCIWP), WM(LWK), WM(LDL),
     3     RPAR, IPAR, IFLAG)
        NLI = NLI + LGMR
        NPS = NPS + NPSL
        CALL DSCAL (N, SQRTN, EWT, 1)
        IF (IFLAG .NE. 0) NCFL = NCFL + 1
        IF (IFLAG .GE. 2) IERSL = 1
        IF (IFLAG .LT. 0) IERSL = -1
        RETURN
      ELSE IF (MITER .EQ. 9) THEN
!-----------------------------------------------------------------------
! Use DVUSOL, which interfaces to PSOL, to solve the linear system
! (No Krylov iteration).
!-----------------------------------------------------------------------
        LB = 1
        LWK = LB + N
        CALL DCOPY (N, X, 1, WM(LB), 1)
        CALL DVUSOL (N, TN, Y, SAVF, WM(LB), EWT, DELTA, HRL1, JPRE,
     1     MNEWT, PSOL, NPSL, X, WM(LOCWP), IWM(LOCIWP), WM(LWK),
     2     RPAR, IPAR, IFLAG)
        NPS = NPS + NPSL
        IF (IFLAG .NE. 0) NCFL = NCFL + 1
        IF (IFLAG .EQ. 3) IERSL = 1
        IF (IFLAG .LT. 0) IERSL = -1
        RETURN
       ENDIF
!----------------------- End of Subroutine DVSLPK ----------------------
      END
*DECK DVSPIG
      SUBROUTINE DVSPIG (TN, Y, SAVF, B, WGHT, N, MAXL, MAXLP1,
     1  KMP, DELTA, HB0, JPRE, MNEWT, F, PSOL, NPSL, X, V, HES, Q,
     2  LGMR, WP, IWP, WK, DL, RPAR, IPAR, IFLAG)
      EXTERNAL F, PSOL
      DOUBLE PRECISION TN, Y, SAVF, B, WGHT, DELTA, HB0, X, V, HES,
     1   Q, WP, WK, DL, RPAR
      INTEGER N, MAXL, MAXLP1, KMP, JPRE, MNEWT, NPSL, LGMR, IWP,
     1   IFLAG, IPAR
      DIMENSION Y(*), SAVF(*), B(*), WGHT(*), X(*), V(N,*),
     1   HES(MAXLP1,*), Q(*), WP(*), IWP(*), WK(*), DL(*),
     2   RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input --  TN, Y, SAVF, B, WGHT, N, MAXL, MAXLP1, DELTA,
!                         HB0, JPRE, MNEWT, F, PSOL, RPAR, IPAR
! Call sequence output -- B, KMP, DELTA, NPSL, X, V, HES, Q, LGMR, WP,
!                         IWP, WK, DL, RPAR, IPAR, IFLAG
! COMMON block variables accessed: None
! Subroutines called: F, DORTHOG, PSOL, DAXPY, DCOPY, DHELS, DHEQR,
!                      DSCAL, DVATV
! Function subroutines called: DNRM2
!-----------------------------------------------------------------------
! This routine solves the linear system  A * x = b using SPIGMR,
! a scaled preconditioned incomplete version of the generalized
! minimum residual method GMRES.
! An initial guess of x = 0 is assumed.
!-----------------------------------------------------------------------
!
!      On entry
!
!           TN = current value of t.
!
!            Y = array containing current dependent variable vector.
!
!         SAVF = array containing current value of f(t,y).
!
!            B = the right hand side of the system A*x = b.
!                B is also used as work space when computing
!                the final approximation.
!                (B is the same as V(*,MAXL+1) in the call to DVSPIG.)
!
!         WGHT = the vector of length N containing the nonzero
!                elements of the diagonal scaling matrix.
!
!            N = the order of the matrix A, and the lengths
!                of the vectors WGHT, B and X.
!
!         MAXL = the maximum allowable order of the matrix HES.
!
!       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
!
!          KMP = the number of previous vectors the new vector VNEW
!                must be made orthogonal to.  KMP .le. MAXL.
!
!        DELTA = tolerance on residuals  b - A*x  in weighted RMS norm.
!
!          HB0 = current value of (step size h) * (coefficient beta0).
!
!         JPRE = preconditioner type flag.
!
!        MNEWT = Newton iteration counter (.ge. 0).
!
!           WK = real work array used by routine DVATV and PSOL.
!
!           DL = real work array used for calculation of the residual
!                norm rho when the method is incomplete (KMP.lt.MAXL).
!
!           WP = real work array used by preconditioner PSOL.
!
!          IWP = integer work array used by preconditioner PSOL.
!
!      On return
!
!         X    = the final computed approximation to the solution
!                of the system A*x = b.
!
!         LGMR = the number of iterations performed and the current
!                order of the upper Hessenberg matrix HES.
!
!         NPSL = the number of calls to PSOL.
!
!         V    = the N by (LGMR+1) array containing the LGMR
!                orthogonal vectors V(*,1) to V(*,LGMR).
!
!         HES  = the upper triangular factor of the QR decomposition
!                of the (LGMR+1) by LGMR upper Hessenberg matrix whose
!                entries are the scaled inner-products of A*V(*,i)
!                and V(*,k).
!
!         Q    = real array of length 2*MAXL containing the components
!                of the Givens rotations used in the QR decomposition
!                of HES.  It is loaded in DHEQR and used in DHELS.
!
!        IFLAG = integer error flag:
!                0 means convergence in LGMR iterations, LGMR.le.MAXL.
!                1 means the convergence test did not pass in MAXL
!                  iterations, but the residual norm is .lt. 1,
!                  or .lt. norm(b) if MNEWT = 0, and so x is computed.
!                2 means the convergence test did not pass in MAXL
!                  iterations, residual .gt. 1, and x is undefined.
!                3 means there was a recoverable error in PSOL
!                  caused by the preconditioner being out of date.
!               -1 means there was a nonrecoverable error in PSOL.
!
!-----------------------------------------------------------------------
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION BNRM, BNRM0, C, DLNRM, PROD, RHO, S, SNORMW, TEM
      INTEGER I, IER, INFO, IP1, I2, J, K, LL, LLP1
!
! Type declaration for function subroutines called ---------------------
!
      DOUBLE PRECISION DNRM2
!
      IFLAG = 0
      LGMR = 0
      NPSL = 0
!-----------------------------------------------------------------------
! The initial residual is the vector b.  Apply scaling to b, and test
! for an immediate return with x = 0 or x = b.
!-----------------------------------------------------------------------
      DO 10 I = 1, N
 10     V(I,1) = B(I)*WGHT(I)
      BNRM0 = DNRM2 (N, V, 1)
      BNRM = BNRM0
      IF (BNRM0 .GT. DELTA) GO TO 30
      IF (MNEWT .GT. 0) GO TO 20
      CALL DCOPY (N, B, 1, X, 1)
      RETURN
 20   DO 25 I = 1, N
 25     X(I) = 0.0D0
      RETURN
 30   CONTINUE
! Apply inverse of left preconditioner to vector b. --------------------
      IER = 0
      IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GO TO 55
      CALL PSOL (N, TN, Y, SAVF, WK, HB0, WP, IWP, B, 1,
     1           IER, RPAR, IPAR)
      NPSL = 1
      IF (IER .NE. 0) GO TO 300
! Calculate norm of scaled vector V(*, 1) and normalize it. ------------
      DO 50 I = 1, N
 50     V(I,1) = B(I)*WGHT(I)
      BNRM = DNRM2 (N, V, 1)
      DELTA = DELTA*(BNRM/BNRM0)
 55   TEM = 1.0D0/BNRM
      CALL DSCAL (N, TEM, V(1,1), 1)
! Zero out the HES array. ----------------------------------------------
      DO 65 J = 1, MAXL
        DO 60 I = 1, MAXLP1
 60        HES (I,J) = 0.0D0
 65     CONTINUE
!-----------------------------------------------------------------------
! Main loop to compute the vectors V(*,2) to V(*,MAXL).
! The running product PROD is needed for the convergence test.
!-----------------------------------------------------------------------
      PROD = 1.0D0
      DO 90 LL = 1, MAXL
        LGMR = LL
!-----------------------------------------------------------------------
! Call routine DVATV to compute VNEW = Abar*v(ll), where Abar is
! the matrix A with scaling and inverse preconditioner factors applied.
! Call routine DORTHOG to orthogonalize the new vector VNEW = V(*,LL+1).
! Call routine DHEQR to update the factors of HES.
!-----------------------------------------------------------------------
        CALL DVATV (Y, SAVF, V(1,LL), WGHT, X, F, PSOL, RPAR, IPAR,
     1              V(1,LL+1), WK, WP, IWP, HB0, JPRE, IER, NPSL)
        IF (IER .NE. 0) GO TO 300
        CALL DORTHOG (V(1,LL+1), V,  HES , N, LL, MAXLP1, KMP, SNORMW)
        HES (LL+1,LL) = SNORMW
        CALL DHEQR (HES, MAXLP1, LL, Q, INFO, LL)
        IF (INFO .EQ. LL) GO TO 120
!-----------------------------------------------------------------------
! Update RHO, the estimate of the norm of the residual b - A*xl.
! If KMP .lt. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
! necessarily orthogonal for LL .gt. KMP.  The vector DL must then
! be computed, and its norm used in the calculation of RHO.
!-----------------------------------------------------------------------
        PROD = PROD*Q(2*LL)
        RHO = ABS(PROD*BNRM)
        IF (LL.GT.KMP .AND. KMP.LT.MAXL) THEN
          IF (LL .EQ. KMP+1) THEN
            CALL DCOPY (N, V(1,1), 1, DL, 1)
            DO 75 I = 1, KMP
              IP1 = I + 1
              I2 = I*2
              S = Q(I2)
              C = Q(I2-1)
              DO 70 K = 1, N
 70             DL(K) = S*DL(K) + C*V(K,IP1)
 75           CONTINUE
            ENDIF
          S = Q(2*LL)
          C = Q(2*LL-1)/SNORMW
          LLP1 = LL + 1
          DO 80 K = 1, N
 80         DL(K) = S*DL(K) + C*V(K,LLP1)
          DLNRM = DNRM2 (N, DL, 1)
          RHO = RHO*DLNRM
          ENDIF
!-----------------------------------------------------------------------
! Test for convergence.  If passed, compute approximation xl.
! If failed and LL .lt. MAXL, then continue iterating.
!-----------------------------------------------------------------------
        IF (RHO .LE. DELTA) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
!-----------------------------------------------------------------------
! Rescale so that the norm of V(1,LL+1) is one.
!-----------------------------------------------------------------------
        TEM = 1.0D0/SNORMW
        CALL DSCAL (N, TEM, V(1,LL+1), 1)
 90     CONTINUE
 100  CONTINUE
      IF (RHO .LE. 1.0D0) GO TO 150
      IF (RHO .LE. BNRM .AND. MNEWT .EQ. 0) GO TO 150
 120  CONTINUE
      IFLAG = 2
      RETURN
 150  IFLAG = 1
!-----------------------------------------------------------------------
! Compute the approximation xl to the solution.
! Since the vector X was used as work space, and the initial guess
! of the Newton correction is zero, X must be reset to zero.
!-----------------------------------------------------------------------
 200  CONTINUE
      LL = LGMR
      LLP1 = LL + 1
      DO 210 K = 1, LLP1
 210    B(K) = 0.0D0
      B(1) = BNRM
      CALL DHELS (HES , MAXLP1, LL, Q, B)
      DO 220 K = 1, N
 220    X(K) = 0.0D0
      DO 230 I = 1, LL
        CALL DAXPY (N, B(I), V(1,I), 1, X, 1)
 230    CONTINUE
      DO 240 I = 1, N
 240    X(I) = X(I)/WGHT(I)
      IF (JPRE .LE. 1) RETURN
      CALL PSOL (N, TN, Y, SAVF, WK, HB0, WP, IWP, X, 2,
     1          IER, RPAR, IPAR)
      NPSL = NPSL + 1
      IF (IER .NE. 0) GO TO 300
      RETURN
!-----------------------------------------------------------------------
! This block handles error returns forced by routine PSOL.
!-----------------------------------------------------------------------
 300  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
!
      RETURN
!----------------------- End of Subroutine DVSPIG ----------------------
      END
*DECK DVATV
      SUBROUTINE DVATV (Y, SAVF, V, WGHT, FTEM, F, PSOL, RPAR, IPAR,
     1                Z, VTEM, WP, IWP, HB0, JPRE, IER, NPSL)
      EXTERNAL F, PSOL
      DOUBLE PRECISION Y, SAVF, V, WGHT, FTEM, RPAR, Z, VTEM, WP, HB0
      INTEGER IPAR, IWP, JPRE, IER, NPSL
      DIMENSION Y(*), SAVF(*), V(*), WGHT(*), FTEM(*), Z(*),
     1   VTEM(*), WP(*), IWP(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, SAVF, V, WGHT, F, PSOL, RPAR, IPAR,
!                        WP, IWP, HB0, NPSL
! Call sequence output --Z, IER, NPSL
! COMMON block variables accessed:
!        /DVOD01/  TN, N
!        /DVOD02/  NFE
! Subroutines called: F, PSOL, DCOPY
! Function subroutines called: DNRM2
!-----------------------------------------------------------------------
! This routine computes the product
!
!   (D-inverse)*(P1-inverse)*(I - hb0*df/dy)*(P2-inverse)*(D*v),
!
! where D is a diagonal scaling matrix, and P1 and P2 are the
! left and right preconditioning matrices, respectively.
! v is assumed to have L2 norm equal to 1.
! The product is stored in Z.  This is computed by a
! difference quotient, a call to F, and two calls to PSOL.
!-----------------------------------------------------------------------
!
!      On entry
!
!            Y = array containing current dependent variable vector.
!
!         SAVF = array containing current value of f(t,y).
!
!            V = real array of length N (can be the same array as Z).
!
!         WGHT = array of length N containing scale factors.
!                1/WGHT(i) are the diagonal elements of the matrix D.
!
!         FTEM = work array of length N.
!
!         VTEM = work array of length N used to store the
!                unscaled version of v.
!
!           WP = real work array used by preconditioner PSOL.
!
!          IWP = integer work array used by preconditioner PSOL.
!
!          HB0 = current value of (step size h) * (coefficient beta0).
!
!         JPRE = preconditioner type flag.
!
!
!      On return
!
!            Z = array of length N containing desired scaled
!                matrix-vector product.
!
!          IER = error flag from PSOL.
!
!         NPSL = the number of calls to PSOL.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION FAC, RNORM, TEMPN
      INTEGER I
!
! Type declaration for function subroutines called ---------------------
!
      DOUBLE PRECISION DNRM2
!
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
!-----------------------------------------------------------------------
! Set vtem = D * v. ----------------------------------------------------
      DO 10 I = 1, N
 10     VTEM(I) = V(I)/WGHT(I)
      IER = 0
      IF (JPRE .GE. 2) GO TO 30
!
! JPRE = 0 or 1.  Save y in Z and increment Y by VTEM. -----------------
      CALL DCOPY (N, Y, 1, Z, 1)
      DO 20 I = 1, N
 20     Y(I) = Z(I) + VTEM(I)
      FAC = HB0
      GO TO 60
!
! JPRE = 2 or 3.  Apply inverse of right preconditioner to VTEM. -------
 30   CONTINUE
      CALL PSOL (N, TN, Y, SAVF, FTEM, HB0, WP, IWP, VTEM, 2,
     1          IER, RPAR, IPAR)
      NPSL = NPSL + 1
      IF (IER .NE. 0) RETURN
! Calculate l-2 norm of (D-inverse) * VTEM. ----------------------------
      DO 40 I = 1, N
 40     Z(I) = VTEM(I)*WGHT(I)
      TEMPN = DNRM2 (N, Z, 1)
      RNORM = 1.0D0/TEMPN
! Save y in Z and increment Y by VTEM/norm. ----------------------------
      CALL DCOPY (N, Y, 1, Z, 1)
      DO 50 I = 1, N
 50     Y(I) = Z(I) + VTEM(I)*RNORM
      FAC = HB0*TEMPN
!
! For all JPRE, call F with incremented Y argument, and restore Y. -----
 60   CONTINUE
      CALL F (N, TN, Y, FTEM, RPAR, IPAR)
      NFE = NFE + 1
      CALL DCOPY (N, Z, 1, Y, 1)
! Set Z = (I - HB0*Jacobian) * VTEM, using difference quotient. --------
      DO 70 I = 1, N
 70     Z(I) = FTEM(I) - SAVF(I)
      DO 80 I = 1, N
 80     Z(I) = VTEM(I) - FAC*Z(I)
! Apply inverse of left preconditioner to Z, if nontrivial. ------------
      IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GO TO 85
      CALL PSOL (N, TN, Y, SAVF, FTEM, HB0, WP, IWP, Z, 1,
     1           IER, RPAR, IPAR)
      NPSL = NPSL + 1
      IF (IER .NE. 0) RETURN
 85   CONTINUE
! Apply D-inverse to Z and return. -------------------------------------
      DO 90 I = 1, N
 90     Z(I) = Z(I)*WGHT(I)
      RETURN
!----------------------- End of Subroutine DVATV -----------------------
      END
*DECK DVUSOL
      SUBROUTINE DVUSOL (N, TN, Y, SAVF, B, WGHT, DELTA, HB0, JPRE,
     1   MNEWT, PSOL, NPSL, X, WP, IWP, WK, RPAR, IPAR, IFLAG)
      EXTERNAL PSOL
      DOUBLE PRECISION TN, Y, SAVF, B, WGHT, DELTA, HB0, X, WP, WK, RPAR
      INTEGER N, JPRE, MNEWT, NPSL, IWP, IPAR, IFLAG
      DIMENSION Y(*), SAVF(*), B(*), WGHT(*), X(*),
     1   WP(*), IWP(*), WK(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! This routine solves the linear system A * x = b using only
! calls to the user-supplied routine PSOL (no Krylov iteration).
! If the norm of the right-hand side vector b is smaller than DELTA,
! the vector x returned is x = b (if MNEWT = 0) or x = 0 otherwise.
! PSOL is called with an LR argument of 1 (if JPRE = 1 or 3),
! then 2 (if JPRE = 2 or 3).
!-----------------------------------------------------------------------
!
!      On entry
!
!          NEQ = problem size, passed to F and PSOL (NEQ(1) = N).
!
!           TN = current value of t.
!
!            Y = array containing current dependent variable vector.
!
!         SAVF = array containing current value of f(t,y).
!
!            B = the right hand side of the system A*x = b.
!
!         WGHT = the vector of length N containing the nonzero
!                elements of the diagonal scaling matrix.
!
!            N = the order of the matrix A, and the lengths
!                of the vectors WGHT, b and x.
!
!        DELTA = tolerance on residuals  b - A*x  in weighted RMS norm.
!
!          HB0 = current value of (step size h) * (coefficient beta0).
!
!         JPRE = preconditioner type flag.
!
!        MNEWT = Newton iteration counter (.ge. 0).
!
!           WK = real work array used by PSOL.
!
!           WP = real work array used by preconditioner PSOL.
!
!          IWP = integer work array used by preconditioner PSOL.
!
!      On return
!
!         X    = the final computed approximation to the solution
!                of the system A*x = b.
!
!         NPSL = the number of calls to PSOL.
!
!        IFLAG = integer error flag:
!                0 means no trouble occurred.
!                3 means there was a recoverable error in PSOL
!                  caused by the preconditioner being out of date.
!               -1 means there was a nonrecoverable error in PSOL.
!
!-----------------------------------------------------------------------
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION BNRM
      INTEGER I, IER
!
! Type declaration for function subroutines called ---------------------
!
      DOUBLE PRECISION DVNORM
!
      IFLAG = 0
      NPSL = 0
!-----------------------------------------------------------------------
! Test for an immediate return with x = 0 or x = b.
!-----------------------------------------------------------------------
      BNRM = DVNORM (N, B, WGHT)
      IF (BNRM .GT. DELTA) GO TO 30
      IF (MNEWT .GT. 0) GO TO 10
      CALL DCOPY (N, B, 1, X, 1)
      RETURN
 10   DO 20 I = 1, N
 20     X(I) = 0.0D0
      RETURN
! Apply inverse of left preconditioner to vector b. --------------------
 30   IER = 0
      IF (JPRE .EQ. 0 .OR. JPRE .EQ. 2) GO TO 40
      CALL PSOL (N, TN, Y, SAVF, WK, HB0, WP, IWP, B, 1,
     1           IER, RPAR,IPAR)
      NPSL = 1
      IF (IER .NE. 0) GO TO 100
! Apply inverse of right preconditioner to result, and copy to X. ------
 40   IF (JPRE .LE. 1) GO TO 50
      CALL PSOL (N, TN, Y, SAVF, WK, HB0, WP, IWP, B, 2,
     1           IER, RPAR, IPAR)
      NPSL = NPSL + 1
      IF (IER .NE. 0) GO TO 100
 50   CALL DCOPY (N, B, 1, X, 1)
      RETURN
!-----------------------------------------------------------------------
! This block handles error returns forced by routine PSOL.
!-----------------------------------------------------------------------
 100  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
      RETURN
!----------------------- End of Subroutine DVUSOL ----------------------
      END
*DECK DVKSRC
      SUBROUTINE DVKSRC (RSAV, ISAV, JOB)
      DOUBLE PRECISION RSAV
      INTEGER ISAV, JOB
      DIMENSION RSAV(*), ISAV(*)
!-----------------------------------------------------------------------
! Call sequence input -- RSAV, ISAV, JOB
! Call sequence output -- RSAV, ISAV
! COMMON block variables accessed: all of /DVOD01/, /DVOD02/, /DVPK01/
!
! Subroutines/functions called by DVKSRC: None
!-----------------------------------------------------------------------
! This routine saves or restores (depending on JOB) the contents of the
! COMMON blocks DVOD01, DVOD02, DVPK01, used internally by DVODPK.
!
! RSAV = real array of length 52 or more.
! ISAV = integer array of length 52 or more.
! JOB  = flag indicating to save or restore the COMMON blocks:
!        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV).
!        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV).
!        A call with JOB = 2 presumes a prior call with JOB = 1.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION RVOD1
      INTEGER IVOD1
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION RVOD2
      INTEGER IVOD2
!
! Type declarations for labeled COMMON block DVPK01 --------------------
!
      DOUBLE PRECISION RVPK1
      INTEGER IVPK1
!
! Type declarations for local variables --------------------------------
!
      INTEGER I, IOFF, LENIV1, LENIV2, LENRV1, LENRV2, LRVK1, LIVK1
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE LENRV1, LENIV1, LENRV2, LENIV2, LRVK1, LIVK1
!-----------------------------------------------------------------------
      COMMON /DVOD01/ RVOD1(48), IVOD1(33)
      COMMON /DVOD02/ RVOD2(1), IVOD2(8)
      COMMON /DVPK01/ RVPK1(3), IVPK1(11)
      DATA LENRV1 /48/, LENIV1 /33/, LENRV2 /1/, LENIV2 /8/,
     1   LRVK1 /3/, LIVK1 /11/
!
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1, LENRV1
 10     RSAV(I) = RVOD1(I)
      DO 12 I = 1, LENRV2
 12     RSAV(LENRV1+I) = RVOD2(I)
      IOFF = LENRV1 + LENRV2
      DO 14 I = 1, LRVK1
 14     RSAV(IOFF+I) = RVPK1(I)
!
      DO 20 I = 1, LENIV1
 20     ISAV(I) = IVOD1(I)
      DO 22 I = 1, LENIV2
 22     ISAV(LENIV1+I) = IVOD2(I)
      IOFF = LENIV1 + LENIV2
      DO 24 I = 1, LIVK1
 24     ISAV(IOFF+I) = IVPK1(I)
!
      RETURN
!
 100  CONTINUE
      DO 110 I = 1, LENRV1
 110     RVOD1(I) = RSAV(I)
      DO 112 I = 1, LENRV2
 112     RVOD2(I) = RSAV(LENRV1+I)
      IOFF = LENRV1 + LENRV2
      DO 114 I = 1, LRVK1
 114    RVPK1(I) = RSAV(IOFF+I)
!
      DO 120 I = 1, LENIV1
 120     IVOD1(I) = ISAV(I)
      DO 122 I = 1, LENIV2
 122     IVOD2(I) = ISAV(LENIV1+I)
      IOFF = LENIV1 + LENIV2
      DO 124 I = 1, LIVK1
 124    IVPK1(I) = ISAV(IOFF+I)
!
      RETURN
!----------------------- End of Subroutine DVKSRC ----------------------
      END
*DECK DVHIN
      SUBROUTINE DVHIN (N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
     1   EWT, ITOL, ATOL, Y, TEMP, H0, NITER, IER)
      EXTERNAL F
      DOUBLE PRECISION T0, Y0, YDOT, RPAR, TOUT, UROUND, EWT, ATOL, Y,
     1   TEMP, H0
      INTEGER N, IPAR, ITOL, NITER, IER
      DIMENSION Y0(*), YDOT(*), EWT(*), ATOL(*), Y(*),
     1   TEMP(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND,
!                        EWT, ITOL, ATOL, Y, TEMP
! Call sequence output -- H0, NITER, IER
! COMMON block variables accessed -- None
!
! Subroutines called by DVHIN:  F
! Function routines called by DVHI: DVNORM
!-----------------------------------------------------------------------
! This routine computes the step size, H0, to be attempted on the
! first step, when the user has not supplied a value for this.
!
! First we check that TOUT - T0 differs significantly from zero.  Then
! an iteration is done to approximate the initial second derivative
! and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
! A bias factor of 1/2 is applied to the resulting h.
! The sign of H0 is inferred from the initial values of TOUT and T0.
!
! Communication with DVHIN is done with the following variables:
!
! N      = Size of ODE system, input.
! T0     = Initial value of independent variable, input.
! Y0     = Vector of initial conditions, input.
! YDOT   = Vector of initial first derivatives, input.
! F      = Name of subroutine for right-hand side f(t,y), input.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
! TOUT   = First output value of independent variable
! UROUND = Machine unit roundoff
! EWT, ITOL, ATOL = Error weights and tolerance parameters
!                   as described in the driver routine, input.
! Y, TEMP = Work arrays of length N.
! H0     = Step size to be attempted, output.
! NITER  = Number of iterations (and of f evaluations) to compute H0,
!          output.
! IER    = The error flag, returned with the value
!          IER = 0  if no trouble occurred, or
!          IER = -1 if TOUT and T0 are considered too close to proceed.
!-----------------------------------------------------------------------
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION AFI, ATOLI, DELYI, H, HALF, HG, HLB, HNEW, HRAT,
     1     HUB, HUN, PT1, T1, TDIST, TROUND, TWO, YDDNRM
      INTEGER I, ITER
!
! Type declaration for function subroutines called ---------------------
!
      DOUBLE PRECISION DVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE HALF, HUN, PT1, TWO
      DATA HALF /0.5D0/, HUN /100.0D0/, PT1 /0.1D0/, TWO /2.0D0/
!
      NITER = 0
      TDIST = ABS(TOUT - T0)
      TROUND = UROUND*MAX(ABS(T0),ABS(TOUT))
      IF (TDIST .LT. TWO*TROUND) GO TO 100
!
! Set a lower bound on h based on the roundoff level in T0 and TOUT. ---
      HLB = HUN*TROUND
! Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. -
      HUB = PT1*TDIST
      ATOLI = ATOL(1)
      DO 10 I = 1, N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        DELYI = PT1*ABS(Y0(I)) + ATOLI
        AFI = ABS(YDOT(I))
        IF (AFI*HUB .GT. DELYI) HUB = DELYI/AFI
 10     CONTINUE
!
! Set initial guess for h as geometric mean of upper and lower bounds. -
      ITER = 0
      HG = SQRT(HLB*HUB)
! If the bounds have crossed, exit with the mean value. ----------------
      IF (HUB .LT. HLB) THEN
        H0 = HG
        GO TO 90
      ENDIF
!
! Looping point for iteration. -----------------------------------------
 50   CONTINUE
! Estimate the second derivative as a difference quotient in f. --------
      H = SIGN (HG, TOUT - T0)
      T1 = T0 + H
      DO 60 I = 1, N
 60     Y(I) = Y0(I) + H*YDOT(I)
      CALL F (N, T1, Y, TEMP, RPAR, IPAR)
      DO 70 I = 1, N
 70     TEMP(I) = (TEMP(I) - YDOT(I))/H
      YDDNRM = DVNORM (N, TEMP, EWT)
! Get the corresponding new value of h. --------------------------------
      IF (YDDNRM*HUB*HUB .GT. TWO) THEN
        HNEW = SQRT(TWO/YDDNRM)
      ELSE
        HNEW = SQRT(HG*HUB)
      ENDIF
      ITER = ITER + 1
!-----------------------------------------------------------------------
! Test the stopping conditions.
! Stop if the new and previous h values differ by a factor of .lt. 2.
! Stop if four iterations have been done.  Also, stop with previous h
! if HNEW/HG .gt. 2 after first iteration, as this probably means that
! the second derivative value is bad because of cancellation error.
!-----------------------------------------------------------------------
      IF (ITER .GE. 4) GO TO 80
      HRAT = HNEW/HG
      IF ( (HRAT .GT. HALF) .AND. (HRAT .LT. TWO) ) GO TO 80
      IF ( (ITER .GE. 2) .AND. (HNEW .GT. TWO*HG) ) THEN
        HNEW = HG
        GO TO 80
      ENDIF
      HG = HNEW
      GO TO 50
!
! Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ----
 80   H0 = HNEW*HALF
      IF (H0 .LT. HLB) H0 = HLB
      IF (H0 .GT. HUB) H0 = HUB
 90   H0 = SIGN(H0, TOUT - T0)
      NITER = ITER
      IER = 0
      RETURN
! Error return for TOUT - T0 too small. --------------------------------
 100  IER = -1
      RETURN
!----------------------- End of Subroutine DVHIN -----------------------
      END
*DECK DVINDY
      SUBROUTINE DVINDY (T, K, YH, LDYH, DKY, IFLAG)
      DOUBLE PRECISION T, YH, DKY
      INTEGER K, LDYH, IFLAG
      DIMENSION YH(LDYH,*), DKY(*)
!-----------------------------------------------------------------------
! Call sequence input -- T, K, YH, LDYH
! Call sequence output -- DKY, IFLAG
! COMMON block variables accessed:
!     /DVOD01/ --  H, TN, UROUND, L, N, NQ
!     /DVOD02/ --  HU
!
! Subroutines called by DVINDY: DSCAL, XERRWD
! Function routines called by DVINDY: None
!-----------------------------------------------------------------------
! DVINDY computes interpolated values of the K-th derivative of the
! dependent variable vector y, and stores it in DKY.  This routine
! is called within the package with K = 0 and T = TOUT, but may
! also be called by the user for any K up to the current order.
! (See detailed instructions in the usage documentation.)
!-----------------------------------------------------------------------
! The computed values in DKY are gotten by interpolation using the
! Nordsieck history array YH.  This array corresponds uniquely to a
! vector-valued polynomial of degree NQCUR or less, and DKY is set
! to the K-th derivative of this polynomial at T.
! The formula for DKY is:
!              q
!  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
!             j=K
! where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
! The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
! communicated by COMMON.  The above sum is done in reverse order.
! IFLAG is returned negative if either K or T is out of bounds.
!
! Discussion above and comments in driver explain all variables.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION C, HUN, R, S, TFUZZ, TN1, TP, ZERO
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
      CHARACTER*80 MSG
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE HUN, ZERO
!
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
      DATA HUN /100.0D0/, ZERO /0.0D0/
!
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TFUZZ = HUN*UROUND*SIGN(ABS(TN) + ABS(HU), HU)
      TP = TN - HU - TFUZZ
      TN1 = TN + TFUZZ
      IF ((T-TP)*(T-TN1) .GT. ZERO) GO TO 90
!
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1, NQ
 10     IC = IC*JJ
 15   C = REAL(IC)
      DO 20 I = 1, N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1, JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1, J
 30       IC = IC*JJ
 35     C = REAL(IC)
        DO 40 I = 1, N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      CALL DSCAL (N, R, DKY, 1)
      RETURN
!
 80   MSG = 'DVINDY-- K (=I1) illegal      '
      CALL XERRWD (MSG, 30, 51, 1, 1, K, 0, 0, ZERO, ZERO)
      IFLAG = -1
      RETURN
 90   MSG = 'DVINDY-- T (=R1) illegal      '
      CALL XERRWD (MSG, 30, 52, 1, 0, 0, 0, 1, T, ZERO)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWD (MSG, 60, 52, 1, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      RETURN
!----------------------- End of Subroutine DVINDY ----------------------
      END
*DECK DVSTEP
      SUBROUTINE DVSTEP (Y, YH, LDYH, YH1, EWT, SAVF, VSAV, ACOR,
     1                  WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR)
      EXTERNAL F, JAC, PSOL, VNLS
      DOUBLE PRECISION Y, YH, YH1, EWT, SAVF, VSAV, ACOR, WM, RPAR
      INTEGER LDYH, IWM, IPAR
      DIMENSION Y(*), YH(LDYH,*), YH1(*), EWT(*), SAVF(*), VSAV(*),
     1   ACOR(*), WM(*), IWM(*), RPAR(*), IPAR(*)
!-----------------------------------------------------------------------
! Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV,
!                        ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR
! Call sequence output -- YH, ACOR, WM, IWM
! COMMON block variables accessed:
!     /DVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13),
!               TQ(5), TN, JCUR, JSTART, KFLAG, KUTH,
!               L, LMAX, MAXORD, N, NEWQ, NQ, NQWAIT
!     /DVOD02/  HU, NCFN, NETF, NFE, NQU, NST
!
! Subroutines called by DVSTEP: F, DAXPY, DCOPY, DSCAL,
!                               DVJUST, VNLS, DVSET
! Function routines called by DVSTEP: DVNORM
!-----------------------------------------------------------------------
! DVSTEP performs one step of the integration of an initial value
! problem for a system of ordinary differential equations.
! DVSTEP calls subroutine VNLS for the solution of the nonlinear system
! arising in the time step.  Thus it is independent of the problem
! Jacobian structure and the type of nonlinear system solution method.
! DVSTEP returns a completion flag KFLAG (in COMMON).
! A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10
! consecutive failures occurred.  On a return with KFLAG negative,
! the values of TN and the YH array are as of the beginning of the last
! step, and H is the last step size attempted.
!
! Communication with DVSTEP is done with the following variables:
!
! Y      = An array of length N used for the dependent variable vector.
! YH     = An LDYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ).  On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! LDYH   = A constant integer .ge. N, the first dimension of YH.
!          N is the number of ODEs in the system.
! YH1    = A one-dimensional array occupying the same space as YH.
! EWT    = An array of length N containing multiplicative weights
!          for local error measurements.  Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = An array of working storage, of length N.
!          also used for input of YH(*,MAXORD+2) when JSTART = -1
!          and MAXORD .lt. the current order NQ.
! VSAV   = A work array of length N passed to subroutine VNLS.
! ACOR   = A work array of length N, used for the accumulated
!          corrections.  On a successful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = Real and integer work arrays associated with matrix
!          operations in VNLS.
! F      = Dummy name for the user supplied subroutine for f.
! JAC    = Dummy name for the user supplied Jacobian subroutine.
! PSOL   = Dummy name for the subroutine passed to VNLS, for
!          possible use there.
! VNLS   = Dummy name for the nonlinear system solving subroutine,
!          whose real name is dependent on the method used.
! RPAR, IPAR = Dummy names for user's real and integer work arrays.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
!
! Type declarations for labeled COMMON block DVOD02 --------------------
!
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION ADDON, BIAS1,BIAS2,BIAS3, CNQUOT, DDN, DSM, DUP,
     1     ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF,
     2     ETAQ, ETAQM1, ETAQP1, FLOTL, ONE, ONEPSM,
     3     R, THRESH, TOLD, ZERO
      INTEGER I, I1, I2, IBACK, J, JB, KFC, KFH, MXNCF, NCF, NFLAG
!
! Type declaration for function subroutines called ---------------------
!
      DOUBLE PRECISION DVNORM
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE ADDON, BIAS1, BIAS2, BIAS3,
     1     ETACF, ETAMIN, ETAMX1, ETAMX2, ETAMX3, ETAMXF, ETAQ, ETAQM1,
     2     KFC, KFH, MXNCF, ONEPSM, THRESH, ONE, ZERO
!-----------------------------------------------------------------------
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
!
      DATA KFC/-3/, KFH/-7/, MXNCF/10/
      DATA ADDON  /1.0D-6/,    BIAS1  /6.0D0/,     BIAS2  /6.0D0/,
     1     BIAS3  /10.0D0/,    ETACF  /0.25D0/,    ETAMIN /0.1D0/,
     2     ETAMXF /0.2D0/,     ETAMX1 /1.0D4/,     ETAMX2 /10.0D0/,
     3     ETAMX3 /10.0D0/,    ONEPSM /1.00001D0/, THRESH /1.5D0/
      DATA ONE/1.0D0/, ZERO/0.0D0/
!
      KFLAG = 0
      TOLD = TN
      NCF = 0
      JCUR = 0
      NFLAG = 0
      IF (JSTART .GT. 0) GO TO 20
      IF (JSTART .EQ. -1) GO TO 100
!-----------------------------------------------------------------------
! On the first call, the order is set to 1, and other variables are
! initialized.  ETAMAX is the maximum ratio by which H can be increased
! in a single step.  It is normally 10, but is larger during the
! first step to compensate for the small initial H.  If a failure
! occurs (in corrector convergence or error test), ETAMAX is set to 1
! for the next increase.
!-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      NQNYH = NQ*LDYH
      TAU(1) = H
      PRL1 = ONE
      RC = ZERO
      ETAMAX = ETAMX1
      NQWAIT = 2
      HSCAL = H
      GO TO 200
!-----------------------------------------------------------------------
! Take preliminary actions on a normal continuation step (JSTART.GT.0).
! If the driver changed H, then ETA must be reset and NEWH set to 1.
! If a change of order was dictated on the previous step, then
! it is done here and appropriate adjustments in the history are made.
! On an order decrease, the history array is adjusted by DVJUST.
! On an order increase, the history array is augmented by a column.
! On a change of step size H, the history array YH is rescaled.
!-----------------------------------------------------------------------
 20   CONTINUE
      IF (KUTH .EQ. 1) THEN
        ETA = MIN(ETA,H/HSCAL)
        NEWH = 1
        ENDIF
 50   IF (NEWH .EQ. 0) GO TO 200
      IF (NEWQ .EQ. NQ) GO TO 150
      IF (NEWQ .LT. NQ) THEN
        CALL DVJUST (YH, LDYH, -1)
        NQ = NEWQ
        L = NQ + 1
        NQWAIT = L
        GO TO 150
        ENDIF
      IF (NEWQ .GT. NQ) THEN
        CALL DVJUST (YH, LDYH, 1)
        NQ = NEWQ
        L = NQ + 1
        NQWAIT = L
        GO TO 150
      ENDIF
!-----------------------------------------------------------------------
! The following block handles preliminaries needed when JSTART = -1.
! If N was reduced, zero out part of YH to avoid undefined references.
! If MAXORD was reduced to a value less than the tentative order NEWQ,
! then NQ is set to MAXORD, and a new H ratio ETA is chosen.
! Otherwise, we take the same preliminary actions as for JSTART .gt. 0.
! In any case, NQWAIT is reset to L = NQ + 1 to prevent further
! changes in order for that many steps.
! The new H ratio ETA is limited by the input H if KUTH = 1,
! by HMIN if KUTH = 0, and by HMXI in any case.
! Finally, the history array YH is rescaled.
!-----------------------------------------------------------------------
 100  CONTINUE
      LMAX = MAXORD + 1
      IF (N .EQ. LDYH) GO TO 120
      I1 = 1 + (NEWQ + 1)*LDYH
      I2 = (MAXORD + 1)*LDYH
      IF (I1 .GT. I2) GO TO 120
      DO 110 I = I1, I2
 110    YH1(I) = ZERO
 120  IF (NEWQ .LE. MAXORD) GO TO 140
      FLOTL = REAL(LMAX)
      IF (MAXORD .LT. NQ-1) THEN
        DDN = DVNORM (N, SAVF, EWT)/TQ(1)
        ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
        ENDIF
      IF (MAXORD .EQ. NQ .AND. NEWQ .EQ. NQ+1) ETA = ETAQ
      IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ+1) THEN
        ETA = ETAQM1
        CALL DVJUST (YH, LDYH, -1)
        ENDIF
      IF (MAXORD .EQ. NQ-1 .AND. NEWQ .EQ. NQ) THEN
        DDN = DVNORM (N, SAVF, EWT)/TQ(1)
        ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL) + ADDON)
        CALL DVJUST (YH, LDYH, -1)
        ENDIF
      ETA = MIN(ETA,ONE)
      NQ = MAXORD
      L = LMAX
 140  IF (KUTH .EQ. 1) ETA = MIN(ETA,ABS(H/HSCAL))
      IF (KUTH .EQ. 0) ETA = MAX(ETA,HMIN/ABS(HSCAL))
      ETA = ETA/MAX(ONE,ABS(HSCAL)*HMXI*ETA)
      NEWH = 1
      NQWAIT = L
      IF (NEWQ .LE. MAXORD) GO TO 50
! Rescale the history array for a change in H by a factor of ETA. ------
 150  R = ONE
      DO 180 J = 2, L
        R = R*ETA
        CALL DSCAL (N, R, YH(1,J), 1 )
 180    CONTINUE
      H = HSCAL*ETA
      HSCAL = H
      RC = RC*ETA
      NQNYH = NQ*LDYH
!-----------------------------------------------------------------------
! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! DVSET is called to calculate all integration coefficients.
! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
!-----------------------------------------------------------------------
 200  TN = TN + H
      I1 = NQNYH + 1
      DO 220 JB = 1, NQ
        I1 = I1 - LDYH
        DO 210 I = I1, NQNYH
 210      YH1(I) = YH1(I) + YH1(I+LDYH)
 220  CONTINUE
      CALL DVSET
      RL1 = ONE/EL(2)
      RC = RC*(RL1/PRL1)
      PRL1 = RL1
!
! Call the nonlinear system solver. ------------------------------------
!
      CALL VNLS (Y, YH, LDYH, VSAV, SAVF, EWT, ACOR, IWM, WM,
     1           F, JAC, PSOL, NFLAG, RPAR, IPAR)
!
      IF (NFLAG .EQ. 0) GO TO 450
!-----------------------------------------------------------------------
! The VNLS routine failed to achieve convergence (NFLAG .NE. 0).
! The YH array is retracted to its values before prediction.
! The step size H is reduced and the step is retried, if possible.
! Otherwise, an error exit is taken.
!-----------------------------------------------------------------------
        NCF = NCF + 1
        NCFN = NCFN + 1
        ETAMAX = ONE
        TN = TOLD
        I1 = NQNYH + 1
        DO 430 JB = 1, NQ
          I1 = I1 - LDYH
          DO 420 I = I1, NQNYH
 420        YH1(I) = YH1(I) - YH1(I+LDYH)
 430      CONTINUE
        IF (NFLAG .LT. -1) GO TO 680
        IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 670
        IF (NCF .EQ. MXNCF) GO TO 670
        ETA = ETACF
        ETA = MAX(ETA,HMIN/ABS(H))
        NFLAG = -1
        GO TO 150
!-----------------------------------------------------------------------
! The corrector has converged (NFLAG = 0).  The local error test is
! made and control passes to statement 500 if it fails.
!-----------------------------------------------------------------------
 450  CONTINUE
      DSM = ACNRM/TQ(2)
      IF (DSM .GT. ONE) GO TO 500
!-----------------------------------------------------------------------
! After a successful step, update the YH and TAU arrays and decrement
! NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved
! for use in a possible order increase on the next step.
! If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2.
!-----------------------------------------------------------------------
      KFLAG = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO 470 IBACK = 1, NQ
        I = L - IBACK
 470    TAU(I+1) = TAU(I)
      TAU(1) = H
      DO 480 J = 1, L
        CALL DAXPY (N, EL(J), ACOR, 1, YH(1,J), 1 )
 480    CONTINUE
      NQWAIT = NQWAIT - 1
      IF ((L .EQ. LMAX) .OR. (NQWAIT .NE. 1)) GO TO 490
      CALL DCOPY (N, ACOR, 1, YH(1,LMAX), 1 )
      CONP = TQ(5)
 490  IF (ETAMAX .NE. ONE) GO TO 560
      IF (NQWAIT .LT. 2) NQWAIT = 2
      NEWQ = NQ
      NEWH = 0
      ETA = ONE
      HNEW = H
      GO TO 690
!-----------------------------------------------------------------------
! The error test failed.  KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again.  Compute the optimum step size for the
! same order.  After repeated failures, H is forced to decrease
! more rapidly.
!-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      NETF = NETF + 1
      NFLAG = -2
      TN = TOLD
      I1 = NQNYH + 1
      DO 520 JB = 1, NQ
        I1 = I1 - LDYH
        DO 510 I = I1, NQNYH
 510      YH1(I) = YH1(I) - YH1(I+LDYH)
 520  CONTINUE
      IF (ABS(H) .LE. HMIN*ONEPSM) GO TO 660
      ETAMAX = ONE
      IF (KFLAG .LE. KFC) GO TO 530
! Compute ratio of new H to current H at the current order. ------------
      FLOTL = REAL(L)
      ETA = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
      ETA = MAX(ETA,HMIN/ABS(H),ETAMIN)
      IF ((KFLAG .LE. -2) .AND. (ETA .GT. ETAMXF)) ETA = ETAMXF
      GO TO 150
!-----------------------------------------------------------------------
! Control reaches this section if 3 or more consecutive failures
! have occurred.  It is assumed that the elements of the YH array
! have accumulated errors of the wrong order.  The order is reduced
! by one, if possible.  Then H is reduced by a factor of 0.1 and
! the step is retried.  After a total of 7 consecutive failures,
! an exit is taken with KFLAG = -1.
!-----------------------------------------------------------------------
 530  IF (KFLAG .EQ. KFH) GO TO 660
      IF (NQ .EQ. 1) GO TO 540
      ETA = MAX(ETAMIN,HMIN/ABS(H))
      CALL DVJUST (YH, LDYH, -1)
      L = NQ
      NQ = NQ - 1
      NQWAIT = L
      GO TO 150
 540  ETA = MAX(ETAMIN,HMIN/ABS(H))
      H = H*ETA
      HSCAL = H
      TAU(1) = H
      CALL F (N, TN, Y, SAVF, RPAR, IPAR)
      NFE = NFE + 1
      DO 550 I = 1, N
 550    YH(I,2) = H*SAVF(I)
      NQWAIT = 10
      GO TO 200
!-----------------------------------------------------------------------
! If NQWAIT = 0, an increase or decrease in order by one is considered.
! Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could
! be multiplied at order q, q-1, or q+1, respectively.
! The largest of these is determined, and the new order and
! step size set accordingly.
! A change of H or NQ is made only if H increases by at least a
! factor of THRESH.  If an order change is considered and rejected,
! then NQWAIT is set to 2 (reconsider it after 2 steps).
!-----------------------------------------------------------------------
! Compute ratio of new H to current H at the current order. ------------
 560  FLOTL = REAL(L)
      ETAQ = ONE/((BIAS2*DSM)**(ONE/FLOTL) + ADDON)
      IF (NQWAIT .NE. 0) GO TO 600
      NQWAIT = 2
      ETAQM1 = ZERO
      IF (NQ .EQ. 1) GO TO 570
! Compute ratio of new H to current H at the current order less one. ---
      DDN = DVNORM (N, YH(1,L), EWT)/TQ(1)
      ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL - ONE)) + ADDON)
 570  ETAQP1 = ZERO
      IF (L .EQ. LMAX) GO TO 580
! Compute ratio of new H to current H at current order plus one. -------
      CNQUOT = (TQ(5)/CONP)*(H/TAU(2))**L
      DO 575 I = 1, N
 575    SAVF(I) = ACOR(I) - CNQUOT*YH(I,LMAX)
      DUP = DVNORM (N, SAVF, EWT)/TQ(3)
      ETAQP1 = ONE/((BIAS3*DUP)**(ONE/(FLOTL + ONE)) + ADDON)
 580  IF (ETAQ .GE. ETAQP1) GO TO 590
      IF (ETAQP1 .GT. ETAQM1) GO TO 620
      GO TO 610
 590  IF (ETAQ .LT. ETAQM1) GO TO 610
 600  ETA = ETAQ
      NEWQ = NQ
      GO TO 630
 610  ETA = ETAQM1
      NEWQ = NQ - 1
      GO TO 630
 620  ETA = ETAQP1
      NEWQ = NQ + 1
      CALL DCOPY (N, ACOR, 1, YH(1,LMAX), 1)
! Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ----
 630  IF (ETA .LT. THRESH .OR. ETAMAX .EQ. ONE) GO TO 640
      ETA = MIN(ETA,ETAMAX)
      ETA = ETA/MAX(ONE,ABS(H)*HMXI*ETA)
      NEWH = 1
      HNEW = H*ETA
      GO TO 690
 640  NEWQ = NQ
      NEWH = 0
      ETA = ONE
      HNEW = H
      GO TO 690
!-----------------------------------------------------------------------
! All returns are made through this section.
! On a successful return, ETAMAX is reset and ACOR is scaled.
!-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  IF (NFLAG .EQ. -2) KFLAG = -3
      IF (NFLAG .EQ. -3) KFLAG = -4
      GO TO 720
 690  ETAMAX = ETAMX3
      IF (NST .LE. 10) ETAMAX = ETAMX2
 700  R = ONE/TQ(2)
      CALL DSCAL (N, R, ACOR, 1)
 720  JSTART = 1
      RETURN
!----------------------- End of Subroutine DVSTEP ----------------------
      END
*DECK DVSET
      SUBROUTINE DVSET
!-----------------------------------------------------------------------
! Call sequence communication: None
! COMMON block variables accessed:
!     /DVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1),
!                 METH, NQ, NQWAIT
!
! Subroutines called by DVSET: None
! Function routines called by DVSET: None
!-----------------------------------------------------------------------
! DVSET is called by DVSTEP and sets coefficients for use there.
!
! For each order NQ, the coefficients in EL are calculated by use of
!  the generating polynomial lambda(x), with coefficients EL(i).
!      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
! For the backward differentiation formulas,
!                                     NQ-1
!      lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) .
!                                     i = 1
! For the Adams formulas,
!                              NQ-1
!      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) ,
!                              i = 1
!      lambda(-1) = 0,    lambda(0) = 1,
! where c is a normalization constant.
! In both cases, xi(i) is defined by
!      H*xi(i) = t sub n  -  t sub (n-i)
!              = H + TAU(1) + TAU(2) + ... TAU(i-1).
!
!
! In addition to variables described previously, communication
! with DVSET uses the following:
!   TAU    = A vector of length 13 containing the past NQ values
!            of H.
!   EL     = A vector of length 13 in which vset stores the
!            coefficients for the corrector formula.
!   TQ     = A vector of length 5 in which vset stores constants
!            used for the convergence test, the error test, and the
!            selection of H at a new order.
!   METH   = The basic method indicator.
!   NQ     = The current order.
!   L      = NQ + 1, the length of the vector stored in EL, and
!            the number of columns of the YH array being used.
!   NQWAIT = A counter controlling the frequency of order changes.
!            An order change is about to be considered if NQWAIT = 1.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION AHATN0, ALPH0, CNQM1, CORTES, CSUM, ELP, EM,
     1     EM0, FLOTI, FLOTL, FLOTNQ, HSUM, ONE, RXI, RXIS, S, SIX,
     2     T1, T2, T3, T4, T5, T6, TWO, XI, ZERO
      INTEGER I, IBACK, J, JP1, NQM1, NQM2
!
      DIMENSION EM(13)
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE CORTES, ONE, SIX, TWO, ZERO
!
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
!
      DATA CORTES /0.1D0/
      DATA ONE  /1.0D0/, SIX /6.0D0/, TWO /2.0D0/, ZERO /0.0D0/
!
      FLOTL = REAL(L)
      NQM1 = NQ - 1
      NQM2 = NQ - 2
      GO TO (100, 200), METH
!
! Set coefficients for Adams methods. ----------------------------------
 100  IF (NQ .NE. 1) GO TO 110
      EL(1) = ONE
      EL(2) = ONE
      TQ(1) = ONE
      TQ(2) = TWO
      TQ(3) = SIX*TQ(2)
      TQ(5) = ONE
      GO TO 300
 110  HSUM = H
      EM(1) = ONE
      FLOTNQ = FLOTL - ONE
      DO 115 I = 2, L
 115    EM(I) = ZERO
      DO 150 J = 1, NQM1
        IF ((J .NE. NQM1) .OR. (NQWAIT .NE. 1)) GO TO 130
        S = ONE
        CSUM = ZERO
        DO 120 I = 1, NQM1
          CSUM = CSUM + S*EM(I)/REAL(I+1)
 120      S = -S
        TQ(1) = EM(NQM1)/(FLOTNQ*CSUM)
 130    RXI = H/HSUM
        DO 140 IBACK = 1, J
          I = (J + 2) - IBACK
 140      EM(I) = EM(I) + EM(I-1)*RXI
        HSUM = HSUM + TAU(J)
 150    CONTINUE
! Compute integral from -1 to 0 of polynomial and of x times it. -------
      S = ONE
      EM0 = ZERO
      CSUM = ZERO
      DO 160 I = 1, NQ
        FLOTI = REAL(I)
        EM0 = EM0 + S*EM(I)/FLOTI
        CSUM = CSUM + S*EM(I)/(FLOTI+ONE)
 160    S = -S
! In EL, form coefficients of normalized integrated polynomial. --------
      S = ONE/EM0
      EL(1) = ONE
      DO 170 I = 1, NQ
 170    EL(I+1) = S*EM(I)/REAL(I)
      XI = HSUM/H
      TQ(2) = XI*EM0/CSUM
      TQ(5) = XI/EL(L)
      IF (NQWAIT .NE. 1) GO TO 300
! For higher order control constant, multiply polynomial by 1+x/xi(q). -
      RXI = ONE/XI
      DO 180 IBACK = 1, NQ
        I = (L + 1) - IBACK
 180    EM(I) = EM(I) + EM(I-1)*RXI
! Compute integral of polynomial. --------------------------------------
      S = ONE
      CSUM = ZERO
      DO 190 I = 1, L
        CSUM = CSUM + S*EM(I)/REAL(I+1)
 190    S = -S
      TQ(3) = FLOTL*EM0/CSUM
      GO TO 300
!
! Set coefficients for BDF methods. ------------------------------------
 200  DO 210 I = 3, L
 210    EL(I) = ZERO
      EL(1) = ONE
      EL(2) = ONE
      ALPH0 = -ONE
      AHATN0 = -ONE
      HSUM = H
      RXI = ONE
      RXIS = ONE
      IF (NQ .EQ. 1) GO TO 240
      DO 230 J = 1, NQM2
! In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------
        HSUM = HSUM + TAU(J)
        RXI = H/HSUM
        JP1 = J + 1
        ALPH0 = ALPH0 - ONE/REAL(JP1)
        DO 220 IBACK = 1, JP1
          I = (J + 3) - IBACK
 220      EL(I) = EL(I) + EL(I-1)*RXI
 230    CONTINUE
      ALPH0 = ALPH0 - ONE/REAL(NQ)
      RXIS = -EL(2) - ALPH0
      HSUM = HSUM + TAU(NQM1)
      RXI = H/HSUM
      AHATN0 = -EL(2) - RXI
      DO 235 IBACK = 1, NQ
        I = (NQ + 2) - IBACK
 235    EL(I) = EL(I) + EL(I-1)*RXIS
 240  T1 = ONE - AHATN0 + ALPH0
      T2 = ONE + REAL(NQ)*T1
      TQ(2) = ABS(ALPH0*T2/T1)
      TQ(5) = ABS(T2/(EL(L)*RXI/RXIS))
      IF (NQWAIT .NE. 1) GO TO 300
      CNQM1 = RXIS/EL(L)
      T3 = ALPH0 + ONE/REAL(NQ)
      T4 = AHATN0 + RXI
      ELP = T3/(ONE - T4 + T3)
      TQ(1) = ABS(ELP/CNQM1)
      HSUM = HSUM + TAU(NQ)
      RXI = H/HSUM
      T5 = ALPH0 - ONE/REAL(NQ+1)
      T6 = AHATN0 - RXI
      ELP = T2/(ONE - T6 + T5)
      TQ(3) = ABS(ELP*RXI*(FLOTL + ONE)*T5)
 300  TQ(4) = CORTES*TQ(2)
      RETURN
!----------------------- End of Subroutine DVSET -----------------------
      END
*DECK DVJUST
      SUBROUTINE DVJUST (YH, LDYH, IORD)
      DOUBLE PRECISION YH
      INTEGER LDYH, IORD
      DIMENSION YH(LDYH,*)
!-----------------------------------------------------------------------
! Call sequence input -- YH, LDYH, IORD
! Call sequence output -- YH
! COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N
! COMMON block variables accessed:
!     /DVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ,
!
! Subroutines called by DVJUST: DAXPY
! Function routines called by DVJUST: None
!-----------------------------------------------------------------------
! This subroutine adjusts the YH array on reduction of order,
! and also when the order is increased for the stiff option (METH = 2).
! Communication with DVJUST uses the following:
! IORD  = An integer flag used when METH = 2 to indicate an order
!         increase (IORD = +1) or an order decrease (IORD = -1).
! HSCAL = Step size H used in scaling of Nordsieck array YH.
!         (If IORD = +1, DVJUST assumes that HSCAL = TAU(1).)
! See References 1 and 2 for details.
!-----------------------------------------------------------------------
!
! Type declarations for labeled COMMON block DVOD01 --------------------
!
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
!
! Type declarations for local variables --------------------------------
!
      DOUBLE PRECISION ALPH0, ALPH1, HSUM, ONE, PROD, T1, XI,XIOLD, ZERO
      INTEGER I, IBACK, J, JP1, LP1, NQM1, NQM2, NQP1
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE ONE, ZERO
!
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
!
      DATA ONE /1.0D0/, ZERO /0.0D0/
!
      IF ((NQ .EQ. 2) .AND. (IORD .NE. 1)) RETURN
      NQM1 = NQ - 1
      NQM2 = NQ - 2
      GO TO (100, 200), METH
!-----------------------------------------------------------------------
! Nonstiff option...
! Check to see if the order is being increased or decreased.
!-----------------------------------------------------------------------
 100  CONTINUE
      IF (IORD .EQ. 1) GO TO 180
! Order decrease. ------------------------------------------------------
      DO 110 J = 1, LMAX
 110    EL(J) = ZERO
      EL(2) = ONE
      HSUM = ZERO
      DO 130 J = 1, NQM2
! Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). -----------------
        HSUM = HSUM + TAU(J)
        XI = HSUM/HSCAL
        JP1 = J + 1
        DO 120 IBACK = 1, JP1
          I = (J + 3) - IBACK
 120      EL(I) = EL(I)*XI + EL(I-1)
 130    CONTINUE
! Construct coefficients of integrated polynomial. ---------------------
      DO 140 J = 2, NQM1
 140    EL(J+1) = REAL(NQ)*EL(J)/REAL(J)
! Subtract correction terms from YH array. -----------------------------
      DO 170 J = 3, NQ
        DO 160 I = 1, N
 160      YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
 170    CONTINUE
      RETURN
! Order increase. ------------------------------------------------------
! Zero out next column in YH array. ------------------------------------
 180  CONTINUE
      LP1 = L + 1
      DO 190 I = 1, N
 190    YH(I,LP1) = ZERO
      RETURN
!-----------------------------------------------------------------------
! Stiff option...
! Check to see if the order is being increased or decreased.
!-----------------------------------------------------------------------
 200  CONTINUE
      IF (IORD .EQ. 1) GO TO 300
! Order decrease. ------------------------------------------------------
      DO 210 J = 1, LMAX
 210    EL(J) = ZERO
      EL(3) = ONE
      HSUM = ZERO
      DO 230 J = 1,NQM2
! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
        HSUM = HSUM + TAU(J)
        XI = HSUM/HSCAL
        JP1 = J + 1
        DO 220 IBACK = 1, JP1
          I = (J + 4) - IBACK
 220      EL(I) = EL(I)*XI + EL(I-1)
 230    CONTINUE
! Subtract correction terms from YH array. -----------------------------
      DO 250 J = 3,NQ
        DO 240 I = 1, N
 240      YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
 250    CONTINUE
      RETURN
! Order increase. ------------------------------------------------------
 300  DO 310 J = 1, LMAX
 310    EL(J) = ZERO
      EL(3) = ONE
      ALPH0 = -ONE
      ALPH1 = ONE
      PROD = ONE
      XIOLD = ONE
      HSUM = HSCAL
      IF (NQ .EQ. 1) GO TO 340
      DO 330 J = 1, NQM1
! Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). ---------------
        JP1 = J + 1
        HSUM = HSUM + TAU(JP1)
        XI = HSUM/HSCAL
        PROD = PROD*XI
        ALPH0 = ALPH0 - ONE/REAL(JP1)
        ALPH1 = ALPH1 + ONE/XI
        DO 320 IBACK = 1, JP1
          I = (J + 4) - IBACK
 320      EL(I) = EL(I)*XIOLD + EL(I-1)
        XIOLD = XI
 330    CONTINUE
 340  CONTINUE
      T1 = (-ALPH0 - ALPH1)/PROD
! Load column L + 1 in YH array. ---------------------------------------
      LP1 = L + 1
      DO 350 I = 1, N
 350    YH(I,LP1) = T1*YH(I,LMAX)
! Add correction terms to YH array. ------------------------------------
      NQP1 = NQ + 1
      DO 370 J = 3, NQP1
        CALL DAXPY (N, EL(J), YH(1,LP1), 1, YH(1,J), 1 )
 370  CONTINUE
      RETURN
!----------------------- End of Subroutine DVJUST ----------------------
      END
*DECK DEWSET
      SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
!***BEGIN PROLOGUE  DEWSET
!***SUBSIDIARY
!***PURPOSE  Set error weight vector.
!***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This subroutine sets the error weight vector EWT according to
!      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
!  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
!  depending on the value of ITOL.
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DEWSET
!**End
      INTEGER N, ITOL
      INTEGER I
      DOUBLE PRECISION RTOL, ATOL, YCUR, EWT
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
!
!***FIRST EXECUTABLE STATEMENT  DEWSET
      GO TO (10, 20, 30, 40), ITOL
 10   CONTINUE
      DO 15 I = 1,N
 15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 20   CONTINUE
      DO 25 I = 1,N
 25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
      RETURN
 30   CONTINUE
      DO 35 I = 1,N
 35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 40   CONTINUE
      DO 45 I = 1,N
 45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
      RETURN
!----------------------- END OF SUBROUTINE DEWSET ----------------------
      END
*DECK DVNORM
      DOUBLE PRECISION FUNCTION DVNORM (N, V, W)
!***BEGIN PROLOGUE  DVNORM
!***SUBSIDIARY
!***PURPOSE  Weighted root-mean-square vector norm.
!***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  This function routine computes the weighted root-mean-square norm
!  of the vector of length N contained in the array V, with weights
!  contained in the array W of length N:
!    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
!
!***SEE ALSO  DLSODE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791129  DATE WRITTEN
!   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
!   890503  Minor cosmetic changes.  (FNF)
!   930809  Renamed to allow single/double precision versions. (ACH)
!***END PROLOGUE  DVNORM
!**End
      INTEGER N,   I
      DOUBLE PRECISION V, W,   SUM
      DIMENSION V(N), W(N)
!
!***FIRST EXECUTABLE STATEMENT  DVNORM
      SUM = 0.0D0
      DO 10 I = 1,N
 10     SUM = SUM + (V(I)*W(I))**2
      DVNORM = SQRT(SUM/N)
      RETURN
!----------------------- END OF FUNCTION DVNORM ------------------------
      END
*DECK DORTHOG
      SUBROUTINE DORTHOG (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
      INTEGER N, LL, LDHES, KMP
      DOUBLE PRECISION VNEW, V, HES, SNORMW
      DIMENSION VNEW(*), V(N,*), HES(LDHES,*)
!-----------------------------------------------------------------------
! This routine orthogonalizes the vector VNEW against the previous
! KMP vectors in the V array.  It uses a modified Gram-Schmidt
! orthogonalization procedure with conditional reorthogonalization.
! This is the version of 28 may 1986.
!-----------------------------------------------------------------------
!
!      On entry
!
!         VNEW = the vector of length N containing a scaled product
!                of the Jacobian and the vector V(*,LL).
!
!         V    = the N x l array containing the previous LL
!                orthogonal vectors v(*,1) to v(*,LL).
!
!         HES  = an LL x LL upper Hessenberg matrix containing,
!                in HES(i,k), k.lt.LL, scaled inner products of
!                A*V(*,k) and V(*,i).
!
!        LDHES = the leading dimension of the HES array.
!
!         N    = the order of the matrix A, and the length of VNEW.
!
!         LL   = the current order of the matrix HES.
!
!          KMP = the number of previous vectors the new vector VNEW
!                must be made orthogonal to (KMP .le. MAXL).
!
!
!      On return
!
!         VNEW = the new vector orthogonal to V(*,i0) to V(*,LL),
!                where i0 = MAX(1, LL-KMP+1).
!
!         HES  = upper Hessenberg matrix with column LL filled in with
!                scaled inner products of A*V(*,LL) and V(*,i).
!
!       SNORMW = L-2 norm of VNEW.
!
!-----------------------------------------------------------------------
      INTEGER I, I0
      DOUBLE PRECISION ARG, DDOT, DNRM2, SUMDSQ, TEM, VNRM
!
! Get norm of unaltered VNEW for later use. ----------------------------
      VNRM = DNRM2 (N, VNEW, 1)
!-----------------------------------------------------------------------
! Do modified Gram-Schmidt on VNEW = A*v(LL).
! Scaled inner products give new column of HES.
! Projections of earlier vectors are subtracted from VNEW.
!-----------------------------------------------------------------------
      I0 = MAX(1,LL-KMP+1)
      DO 10 I = I0,LL
        HES(I,LL) = DDOT (N, V(1,I), 1, VNEW, 1)
        TEM = -HES(I,LL)
        CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
 10     CONTINUE
!-----------------------------------------------------------------------
! Compute SNORMW = norm of VNEW.
! If VNEW is small compared to its input value (in norm), then
! reorthogonalize VNEW to V(*,1) through V(*,LL).
! Correct if relative correction exceeds 1000*(unit roundoff).
! finally, correct SNORMW using the dot products involved.
!-----------------------------------------------------------------------
      SNORMW = DNRM2 (N, VNEW, 1)
      IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0.0D0
      DO 30 I = I0,LL
        TEM = -DDOT (N, V(1,I), 1, VNEW, 1)
        IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GO TO 30
        HES(I,LL) = HES(I,LL) - TEM
        CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
        SUMDSQ = SUMDSQ + TEM**2
 30     CONTINUE
      IF (SUMDSQ .EQ. 0.0D0) RETURN
      ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
!
      RETURN
!----------------------- End of Subroutine DORTHOG ---------------------
      END
*DECK DHEQR
      SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)
      INTEGER LDA, N, INFO, IJOB
      DOUBLE PRECISION A(LDA,*), Q(*)
!-----------------------------------------------------------------------
!     This routine performs a QR decomposition of an upper
!     Hessenberg matrix A.  There are two options available:
!
!          (1)  performing a fresh decomposition
!          (2)  updating the QR factors by adding a row and a
!               column to the matrix A.
!-----------------------------------------------------------------------
!     DHEQR decomposes an upper Hessenberg matrix by using Givens
!     rotations.
!
!     On entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the matrix to be decomposed.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                A is an (N+1) by N Hessenberg matrix.
!
!        IJOB    INTEGER
!                = 1     means that a fresh decomposition of the
!                        matrix A is desired.
!                .ge. 2  means that the current decomposition of A
!                        will be updated by the addition of a row
!                        and a column.
!     On return
!
!        A       the upper triangular matrix R.
!                The factorization can be written Q*A = R, where
!                Q is a product of Givens rotations and R is upper
!                triangular.
!
!        Q       DOUBLE PRECISION(2*N)
!                the factors c and s of each Givens rotation used
!                in decomposing A.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = k  if  A(k,k) .eq. 0.0 .  This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DHELS will divide by zero
!                     if called.
!
!     Modification of LINPACK, by Peter Brown, LLNL.
!     Written 1/13/86.  This version dated 6/20/01.
!-----------------------------------------------------------------------
      INTEGER I, IQ, J, K, KM1, KP1, NM1
      DOUBLE PRECISION C, S, T, T1, T2
!
      IF (IJOB .GT. 1) GO TO 70
!
! A new facorization is desired.
!
!     QR decomposition without pivoting
!
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
!
!           Compute kth column of R.
!           First, multiply the kth column of A by the previous
!           k-1 Givens rotations.
!
            IF (KM1 .LT. 1) GO TO 20
            DO 10 J = 1, KM1
              I = 2*(J-1) + 1
              T1 = A(J,K)
              T2 = A(J+1,K)
              C = Q(I)
              S = Q(I+1)
              A(J,K) = C*T1 - S*T2
              A(J+1,K) = S*T1 + C*T2
   10         CONTINUE
!
!           Compute Givens components c and s
!
   20       CONTINUE
            IQ = 2*KM1 + 1
            T1 = A(K,K)
            T2 = A(KP1,K)
            IF (T2 .NE. 0.0D0) GO TO 30
              C = 1.0D0
              S = 0.0D0
              GO TO 50
   30       CONTINUE
            IF (ABS(T2) .LT. ABS(T1)) GO TO 40
              T = T1/T2
              S = -1.0D0/SQRT(1.0D0+T*T)
              C = -S*T
              GO TO 50
   40       CONTINUE
              T = T2/T1
              C = 1.0D0/SQRT(1.0D0+T*T)
              S = -C*T
   50       CONTINUE
            Q(IQ) = C
            Q(IQ+1) = S
            A(K,K) = C*T1 - S*T2
            IF (A(K,K) .EQ. 0.0D0) INFO = K
   60 CONTINUE
      RETURN
!
! The old factorization of A will be updated.  A row and a column
! has been added to the matrix A.
! N by N-1 is now the old size of the matrix.
!
  70  CONTINUE
      NM1 = N - 1
!
! Multiply the new column by the N previous Givens rotations.
!
      DO 100 K = 1,NM1
        I = 2*(K-1) + 1
        T1 = A(K,N)
        T2 = A(K+1,N)
        C = Q(I)
        S = Q(I+1)
        A(K,N) = C*T1 - S*T2
        A(K+1,N) = S*T1 + C*T2
 100    CONTINUE
!
! Complete update of decomposition by forming last Givens rotation,
! and multiplying it times the column vector (A(N,N), A(N+1,N)).
!
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF (T2 .NE. 0.0D0) GO TO 110
        C = 1.0D0
        S = 0.0D0
        GO TO 130
 110  CONTINUE
      IF (ABS(T2) .LT. ABS(T1)) GO TO 120
        T = T1/T2
        S = -1.0D0/SQRT(1.0D0+T*T)
        C = -S*T
        GO TO 130
 120  CONTINUE
        T = T2/T1
        C = 1.0D0/SQRT(1.0D0+T*T)
        S = -C*T
 130  CONTINUE
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
!----------------------- End of Subroutine DHEQR -----------------------
      END
*DECK DHELS
      SUBROUTINE DHELS (A, LDA, N, Q, B)
      INTEGER LDA, N
      DOUBLE PRECISION A(LDA,*), B(*), Q(*)
!-----------------------------------------------------------------------
! This is part of the LINPACK routine DGESL with changes
! due to the fact that A is an upper Hessenberg matrix.
!-----------------------------------------------------------------------
!     DHELS solves the least squares problem
!
!           min (b-A*x, b-A*x)
!
!     using the factors computed by DHEQR.
!
!     On entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the output from DHEQR which contains the upper
!                triangular factor R in the QR decomposition of A.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                A is originally an (N+1) by N matrix.
!
!        Q       DOUBLE PRECISION(2*N)
!                The coefficients of the N givens rotations
!                used in the QR factorization of A.
!
!        B       DOUBLE PRECISION(N+1)
!                the right hand side vector.
!
!     On return
!
!        B       the solution vector  x .
!
!     Modification of LINPACK, by Peter Brown, LLNL.
!     Written 1/13/86.  This version dated 6/20/01.
!
!     BLAS called: DAXPY
!-----------------------------------------------------------------------
      INTEGER IQ, K, KB, KP1
      DOUBLE PRECISION C, S, T, T1, T2
!
!        Minimize (b-A*x, b-A*x)
!        First form Q*b.
!
         DO 20 K = 1, N
            KP1 = K + 1
            IQ = 2*(K-1) + 1
            C = Q(IQ)
            S = Q(IQ+1)
            T1 = B(K)
            T2 = B(KP1)
            B(K) = C*T1 - S*T2
            B(KP1) = S*T1 + C*T2
   20    CONTINUE
!
!        Now solve  R*x = Q*b.
!
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY (K-1, T, A(1,K), 1, B(1), 1)
   40    CONTINUE
      RETURN
!----------------------- End of Subroutine DHELS -----------------------
      END
*DECK XERRWD
      SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
!***BEGIN PROLOGUE  XERRWD
!***SUBSIDIARY
!***PURPOSE  Write error message with values.
!***CATEGORY  R3C
!***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
!  as given here, constitute a simplified version of the SLATEC error
!  handling package.
!
!  All arguments are input arguments.
!
!  MSG    = The message (character array).
!  NMES   = The length of MSG (number of characters).
!  NERR   = The error number (not used).
!  LEVEL  = The error level..
!           0 or 1 means recoverable (control returns to caller).
!           2 means fatal (run is aborted--see note below).
!  NI     = Number of integers (0, 1, or 2) to be printed with message.
!  I1,I2  = Integers to be printed, depending on NI.
!  NR     = Number of reals (0, 1, or 2) to be printed with message.
!  R1,R2  = Reals to be printed, depending on NR.
!
!  Note..  this routine is machine-dependent and specialized for use
!  in limited context, in the following ways..
!  1. The argument MSG is assumed to be of type CHARACTER, and
!     the message is printed with a format of (1X,A).
!  2. The message is assumed to take only one line.
!     Multi-line messages are generated by repeated calls.
!  3. If LEVEL = 2, control passes to the statement   STOP
!     to abort the run.  This statement may be machine-dependent.
!  4. R1 and R2 are assumed to be in double precision and are printed
!     in D21.13 format.
!
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   920831  DATE WRITTEN
!   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
!   930329  Modified prologue to SLATEC format. (FNF)
!   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
!   930922  Minor cosmetic change. (FNF)
!***END PROLOGUE  XERRWD
!
!*Internal Notes:
!
! For a different default logical unit number, IXSAV (or a subsidiary
! routine that it calls) will need to be modified.
! For a different run-abort command, change the statement following
! statement 100 at the end.
!-----------------------------------------------------------------------
! Subroutines called by XERRWD.. None
! Function routine called by XERRWD.. IXSAV
!-----------------------------------------------------------------------
!**End
!
!  Declare arguments.
!
      DOUBLE PRECISION R1, R2
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
      CHARACTER*(*) MSG
!
!  Declare local variables.
!
      INTEGER LUNIT, IXSAV, MESFLG
!
!  Get logical unit number and message print flag.
!
!***FIRST EXECUTABLE STATEMENT  XERRWD
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
!
!  Write the message.
!
      WRITE (LUNIT,10)  MSG
 10   FORMAT(1X,A)
      IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
 20   FORMAT(6X,'In above message,  I1 =',I10)
      IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
 30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
      IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
 40   FORMAT(6X,'In above message,  R1 =',D21.13)
      IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
 50   FORMAT(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)
!
!  Abort the run if LEVEL = 2.
!
 100  IF (LEVEL .NE. 2) RETURN
      STOP
!----------------------- End of Subroutine XERRWD ----------------------
      END
*DECK XSETF
      SUBROUTINE XSETF (MFLAG)
!***BEGIN PROLOGUE  XSETF
!***PURPOSE  Reset the error print control flag.
!***CATEGORY  R3A
!***TYPE      ALL (XSETF-A)
!***KEYWORDS  ERROR CONTROL
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!   XSETF sets the error print control flag to MFLAG:
!      MFLAG=1 means print all messages (the default).
!      MFLAG=0 means no printing.
!
!***SEE ALSO  XERRWD, XERRWV
!***REFERENCES  (NONE)
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Added SLATEC format prologue. (FNF)
!   930407  Corrected SEE ALSO section. (FNF)
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  XSETF
!
! Subroutines called by XSETF.. None
! Function routine called by XSETF.. IXSAV
!-----------------------------------------------------------------------
!**End
      INTEGER MFLAG, JUNK, IXSAV
!
!***FIRST EXECUTABLE STATEMENT  XSETF
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV (2,MFLAG,.TRUE.)
      RETURN
!----------------------- End of Subroutine XSETF -----------------------
      END
*DECK XSETUN
      SUBROUTINE XSETUN (LUN)
!***BEGIN PROLOGUE  XSETUN
!***PURPOSE  Reset the logical unit number for error messages.
!***CATEGORY  R3B
!***TYPE      ALL (XSETUN-A)
!***KEYWORDS  ERROR CONTROL
!***DESCRIPTION
!
!   XSETUN sets the logical unit number for error messages to LUN.
!
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***SEE ALSO  XERRWD, XERRWV
!***REFERENCES  (NONE)
!***ROUTINES CALLED  IXSAV
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Added SLATEC format prologue. (FNF)
!   930407  Corrected SEE ALSO section. (FNF)
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  XSETUN
!
! Subroutines called by XSETUN.. None
! Function routine called by XSETUN.. IXSAV
!-----------------------------------------------------------------------
!**End
      INTEGER LUN, JUNK, IXSAV
!
!***FIRST EXECUTABLE STATEMENT  XSETUN
      IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
      RETURN
!----------------------- End of Subroutine XSETUN ----------------------
      END
*DECK IXSAV
      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
!***BEGIN PROLOGUE  IXSAV
!***SUBSIDIARY
!***PURPOSE  Save and recall error message control parameters.
!***CATEGORY  R3C
!***TYPE      ALL (IXSAV-A)
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
!
!  IXSAV saves and recalls one of two error message parameters:
!    LUNIT, the logical unit number to which messages are printed, and
!    MESFLG, the message print flag.
!  This is a modification of the SLATEC library routine J4SAVE.
!
!  Saved local variables..
!   LUNIT  = Logical unit number for messages.  The default is obtained
!            by a call to IUMACH (may be machine-dependent).
!   MESFLG = Print control flag..
!            1 means print all messages (the default).
!            0 means no printing.
!
!  On input..
!    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
!    IVALUE = The value to be set for the parameter, if ISET = .TRUE.
!    ISET   = Logical flag to indicate whether to read or write.
!             If ISET = .TRUE., the parameter will be given
!             the value IVALUE.  If ISET = .FALSE., the parameter
!             will be unchanged, and IVALUE is a dummy argument.
!
!  On return..
!    IXSAV = The (old) value of the parameter.
!
!***SEE ALSO  XERRWD, XERRWV
!***ROUTINES CALLED  IUMACH
!***REVISION HISTORY  (YYMMDD)
!   921118  DATE WRITTEN
!   930329  Modified prologue to SLATEC format. (FNF)
!   930915  Added IUMACH call to get default output unit.  (ACH)
!   930922  Minor cosmetic changes. (FNF)
!   010425  Type declaration for IUMACH added. (ACH)
!***END PROLOGUE  IXSAV
!
! Subroutines called by IXSAV.. None
! Function routine called by IXSAV.. IUMACH
!-----------------------------------------------------------------------
!**End
      LOGICAL ISET
      INTEGER IPAR, IVALUE
!-----------------------------------------------------------------------
      INTEGER IUMACH, LUNIT, MESFLG
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this routine.
!-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/-1/, MESFLG/1/
!
!***FIRST EXECUTABLE STATEMENT  IXSAV
      IF (IPAR .EQ. 1) THEN
        IF (LUNIT .EQ. -1) LUNIT = IUMACH()
        IXSAV = LUNIT
        IF (ISET) LUNIT = IVALUE
        ENDIF
!
      IF (IPAR .EQ. 2) THEN
        IXSAV = MESFLG
        IF (ISET) MESFLG = IVALUE
        ENDIF
!
      RETURN
!----------------------- End of Function IXSAV -------------------------
      END
*DECK IUMACH
      INTEGER FUNCTION IUMACH()
!***BEGIN PROLOGUE  IUMACH
!***PURPOSE  Provide standard output unit number.
!***CATEGORY  R1
!***TYPE      INTEGER (IUMACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        INTEGER  LOUT, IUMACH
!        LOUT = IUMACH()
!
! *Function Return Values:
!     LOUT : the standard logical unit for Fortran output.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   930915  DATE WRITTEN
!   930922  Made user-callable, and other cosmetic changes. (FNF)
!***END PROLOGUE  IUMACH
!
!*Internal Notes:
!  The built-in value of 6 is standard on a wide range of Fortran
!  systems.  This may be machine-dependent.
!**End
!***FIRST EXECUTABLE STATEMENT  IUMACH
      IUMACH = 6
!
      RETURN
!----------------------- End of Function IUMACH ------------------------
      END
*DECK DUMACH
      DOUBLE PRECISION FUNCTION DUMACH ()
!***BEGIN PROLOGUE  DUMACH
!***PURPOSE  Compute the unit roundoff of the machine.
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Hindmarsh, Alan C., (LLNL)
!***DESCRIPTION
! *Usage:
!        DOUBLE PRECISION  A, DUMACH
!        A = DUMACH()
!
! *Function Return Values:
!     A : the unit roundoff of the machine.
!
! *Description:
!     The unit roundoff is defined as the smallest positive machine
!     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
!     in a machine-independent manner.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DUMSUM
!***REVISION HISTORY  (YYYYMMDD)
!   19930216  DATE WRITTEN
!   19930818  Added SLATEC-format prologue.  (FNF)
!   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
!***END PROLOGUE  DUMACH
!
      DOUBLE PRECISION U, COMP
!***FIRST EXECUTABLE STATEMENT  DUMACH
      U = 1.0D0
 10   U = U*0.5D0
      CALL DUMSUM(1.0D0, U, COMP)
      IF (COMP .NE. 1.0D0) GO TO 10
      DUMACH = U*2.0D0
      RETURN
!----------------------- End of Function DUMACH ------------------------
      END
      SUBROUTINE DUMSUM(A,B,C)
!     Routine to force normal storing of A + B, for DUMACH.
      DOUBLE PRECISION A, B, C
      C = A + B
      RETURN
      END

