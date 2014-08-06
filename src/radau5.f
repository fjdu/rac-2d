      subroutine radau5(n,fcn,x,y,xend,h,
     &                  rtol,atol,itol,
     &                  jac ,ijac,mljac,mujac,
     &                  mas ,imas,mlmas,mumas,
     &                  solout,iout,
     &                  work,lwork,iwork,liwork,rpar,ipar,idid)
! ----------------------------------------------------------
!     numerical solution of a stiff (or differential algebraic)
!     system of first 0rder ordinary differential equations
!                     m*y'=f(x,y).
!     the system can be (linearly) implicit (mass-matrix m .ne. i)
!     or explicit (m=i).
!     the method used is an implicit runge-kutta method (radau iia)
!     of order 5 with step size control and continuous output.
!     cf. section iv.8
!
!     authors: e. hairer and g. wanner
!              universite de geneve, dept. de mathematiques
!              ch-1211 geneve 24, switzerland 
!              e-mail:  ernst.hairer@math.unige.ch
!                       gerhard.wanner@math.unige.ch
!     
!     this code is part of the book:
!         e. hairer and g. wanner, solving ordinary differential
!         equations ii. stiff and differential-algebraic problems.
!         springer series in computational mathematics 14,
!         springer-verlag 1991, second edition 1996.
!      
!     version of july 9, 1996
!     (latest small correction: january 18, 2002)
!
!     input parameters  
!     ----------------  
!     n           dimension of the system 
!
!     fcn         name (external) of subroutine computing the
!                 value of f(x,y):
!                    subroutine fcn(n,x,y,f,rpar,ipar)
!                    double precision x,y(n),f(n)
!                    f(1)=...   etc.
!                 rpar, ipar (see below)
!
!     x           initial x-value
!
!     y(n)        initial values for y
!
!     xend        final x-value (xend-x may be positive or negative)
!
!     h           initial step size guess;
!                 for stiff equations with initial transient, 
!                 h=1.d0/(norm of f'), usually 1.d-3 or 1.d-5, is good.
!                 this choice is not very important, the step size is
!                 quickly adapted. (if h=0.d0, the code puts h=1.d-6).
!
!     rtol,atol   relative and absolute error tolerances. they
!                 can be both scalars or else both vectors of length n.
!
!     itol        switch for rtol and atol:
!                   itol=0: both rtol and atol are scalars.
!                     the code keeps, roughly, the local error of
!                     y(i) below rtol*abs(y(i))+atol
!                   itol=1: both rtol and atol are vectors.
!                     the code keeps the local error of y(i) below
!                     rtol(i)*abs(y(i))+atol(i).
!
!     jac         name (external) of the subroutine which computes
!                 the partial derivatives of f(x,y) with respect to y
!                 (this routine is only called if ijac=1; supply
!                 a dummy subroutine in the case ijac=0).
!                 for ijac=1, this subroutine must have the form
!                    subroutine jac(n,x,y,dfy,ldfy,rpar,ipar)
!                    double precision x,y(n),dfy(ldfy,n)
!                    dfy(1,1)= ...
!                 ldfy, the column-length of the array, is
!                 furnished by the calling program.
!                 if (mljac.eq.n) the jacobian is supposed to
!                    be full and the partial derivatives are
!                    stored in dfy as
!                       dfy(i,j) = partial f(i) / partial y(j)
!                 else, the jacobian is taken as banded and
!                    the partial derivatives are stored
!                    diagonal-wise as
!                       dfy(i-j+mujac+1,j) = partial f(i) / partial y(j).
!
!     ijac        switch for the computation of the jacobian:
!                    ijac=0: jacobian is computed internally by finite
!                       differences, subroutine "jac" is never called.
!                    ijac=1: jacobian is supplied by subroutine jac.
!
!     mljac       switch for the banded structure of the jacobian:
!                    mljac=n: jacobian is a full matrix. the linear
!                       algebra is done by full-matrix gauss-elimination.
!                    0<=mljac<n: mljac is the lower bandwith of jacobian 
!                       matrix (>= number of non-zero diagonals below
!                       the main diagonal).
!
!     mujac       upper bandwith of jacobian  matrix (>= number of non-
!                 zero diagonals above the main diagonal).
!                 need not be defined if mljac=n.
!
!     ----   mas,imas,mlmas, and mumas have analog meanings      -----
!     ----   for the "mass matrix" (the matrix "m" of section iv.8): -
!
!     mas         name (external) of subroutine computing the mass-
!                 matrix m.
!                 if imas=0, this matrix is assumed to be the identity
!                 matrix and needs not to be defined;
!                 supply a dummy subroutine in this case.
!                 if imas=1, the subroutine mas is of the form
!                    subroutine mas(n,am,lmas,rpar,ipar)
!                    double precision am(lmas,n)
!                    am(1,1)= ....
!                    if (mlmas.eq.n) the mass-matrix is stored
!                    as full matrix like
!                         am(i,j) = m(i,j)
!                    else, the matrix is taken as banded and stored
!                    diagonal-wise as
!                         am(i-j+mumas+1,j) = m(i,j).
!
!     imas       gives information on the mass-matrix:
!                    imas=0: m is supposed to be the identity
!                       matrix, mas is never called.
!                    imas=1: mass-matrix  is supplied.
!
!     mlmas       switch for the banded structure of the mass-matrix:
!                    mlmas=n: the full matrix case. the linear
!                       algebra is done by full-matrix gauss-elimination.
!                    0<=mlmas<n: mlmas is the lower bandwith of the
!                       matrix (>= number of non-zero diagonals below
!                       the main diagonal).
!                 mlmas is supposed to be .le. mljac.
!
!     mumas       upper bandwith of mass-matrix (>= number of non-
!                 zero diagonals above the main diagonal).
!                 need not be defined if mlmas=n.
!                 mumas is supposed to be .le. mujac.
!
!     solout      name (external) of subroutine providing the
!                 numerical solution during integration. 
!                 if iout=1, it is called after every successful step.
!                 supply a dummy subroutine if iout=0. 
!                 it must have the form
!                    subroutine solout (nr,xold,x,y,cont,lrc,n,
!                                       rpar,ipar,irtrn)
!                    double precision x,y(n),cont(lrc)
!                    ....  
!                 solout furnishes the solution "y" at the nr-th
!                    grid-point "x" (thereby the initial value is
!                    the first grid-point).
!                 "xold" is the preceeding grid-point.
!                 "irtrn" serves to interrupt the integration. if irtrn
!                    is set <0, radau5 returns to the calling program.
!           
!          -----  continuous output: -----
!                 during calls to "solout", a continuous solution
!                 for the interval [xold,x] is available through
!                 the function
!                        >>>   contr5(i,s,cont,lrc)   <<<
!                 which provides an approximation to the i-th
!                 component of the solution at the point s. the value
!                 s should lie in the interval [xold,x].
!                 do not change the entries of cont(lrc), if the
!                 dense output function is used.
!
!     iout        switch for calling the subroutine solout:
!                    iout=0: subroutine is never called
!                    iout=1: subroutine is available for output.
!
!     work        array of working space of length "lwork".
!                 work(1), work(2),.., work(20) serve as parameters
!                 for the code. for standard use of the code
!                 work(1),..,work(20) must be set to zero before
!                 calling. see below for a more sophisticated use.
!                 work(21),..,work(lwork) serve as working space
!                 for all vectors and matrices.
!                 "lwork" must be at least
!                             n*(ljac+lmas+3*le+12)+20
!                 where
!                    ljac=n              if mljac=n (full jacobian)
!                    ljac=mljac+mujac+1  if mljac<n (banded jac.)
!                 and                  
!                    lmas=0              if imas=0
!                    lmas=n              if imas=1 and mlmas=n (full)
!                    lmas=mlmas+mumas+1  if mlmas<n (banded mass-m.)
!                 and
!                    le=n               if mljac=n (full jacobian)
!                    le=2*mljac+mujac+1 if mljac<n (banded jac.)
!
!                 in the usual case where the jacobian is full and the
!                 mass-matrix is the indentity (imas=0), the minimum
!                 storage requirement is 
!                             lwork = 4*n*n+12*n+20.
!                 if iwork(9)=m1>0 then "lwork" must be at least
!                          n*(ljac+12)+(n-m1)*(lmas+3*le)+20
!                 where in the definitions of ljac, lmas and le the
!                 number n can be replaced by n-m1.
!
!     lwork       declared length of array "work".
!
!     iwork       integer working space of length "liwork".
!                 iwork(1),iwork(2),...,iwork(20) serve as parameters
!                 for the code. for standard use, set iwork(1),..,
!                 iwork(20) to zero before calling.
!                 iwork(21),...,iwork(liwork) serve as working area.
!                 "liwork" must be at least 3*n+20.
!
!     liwork      declared length of array "iwork".
!
!     rpar, ipar  real and integer parameters (or parameter arrays) which  
!                 can be used for communication between your calling
!                 program and the fcn, jac, mas, solout subroutines. 
!
! ----------------------------------------------------------------------
! 
!     sophisticated setting of parameters
!     -----------------------------------
!              several parameters of the code are tuned to make it work 
!              well. they may be defined by setting work(1),...
!              as well as iwork(1),... different from zero.
!              for zero input, the code chooses default values:
!
!    iwork(1)  if iwork(1).ne.0, the code transforms the jacobian
!              matrix to hessenberg form. this is particularly
!              advantageous for large systems with full jacobian.
!              it does not work for banded jacobian (mljac<n)
!              and not for implicit systems (imas=1).
!
!    iwork(2)  this is the maximal number of allowed steps.
!              the default value (for iwork(2)=0) is 100000.
!
!    iwork(3)  the maximum number of newton iterations for the
!              solution of the implicit system in each step.
!              the default value (for iwork(3)=0) is 7.
!
!    iwork(4)  if iwork(4).eq.0 the extrapolated collocation solution
!              is taken as starting value for newton's method.
!              if iwork(4).ne.0 zero starting values are used.
!              the latter is recommended if newton's method has
!              difficulties with convergence (this is the case when
!              nstep is larger than naccpt + nrejct; see output param.).
!              default is iwork(4)=0.
!
!       the following 3 parameters are important for
!       differential-algebraic systems of index > 1.
!       the function-subroutine should be written such that
!       the index 1,2,3 variables appear in this order. 
!       in estimating the error the index 2 variables are
!       multiplied by h, the index 3 variables by h**2.
!
!    iwork(5)  dimension of the index 1 variables (must be > 0). for 
!              ode's this equals the dimension of the system.
!              default iwork(5)=n.
!
!    iwork(6)  dimension of the index 2 variables. default iwork(6)=0.
!
!    iwork(7)  dimension of the index 3 variables. default iwork(7)=0.
!
!    iwork(8)  switch for step size strategy
!              if iwork(8).eq.1  mod. predictive controller (gustafsson)
!              if iwork(8).eq.2  classical step size control
!              the default value (for iwork(8)=0) is iwork(8)=1.
!              the choice iwork(8).eq.1 seems to produce safer results;
!              for simple problems, the choice iwork(8).eq.2 produces
!              often slightly faster runs
!
!       if the differential system has the special structure that
!            y(i)' = y(i+m2)   for  i=1,...,m1,
!       with m1 a multiple of m2, a substantial gain in computertime
!       can be achieved by setting the parameters iwork(9) and iwork(10).
!       e.g., for second order systems p'=v, v'=g(p,v), where p and v are 
!       vectors of dimension n/2, one has to put m1=m2=n/2.
!       for m1>0 some of the input parameters have different meanings:
!       - jac: only the elements of the non-trivial part of the
!              jacobian have to be stored
!              if (mljac.eq.n-m1) the jacobian is supposed to be full
!                 dfy(i,j) = partial f(i+m1) / partial y(j)
!                for i=1,n-m1 and j=1,n.
!              else, the jacobian is banded ( m1 = m2 * mm )
!                 dfy(i-j+mujac+1,j+k*m2) = partial f(i+m1) / partial y(j+k*m2)
!                for i=1,mljac+mujac+1 and j=1,m2 and k=0,mm.
!       - mljac: mljac=n-m1: if the non-trivial part of the jacobian is full
!                0<=mljac<n-m1: if the (mm+1) submatrices (for k=0,mm)
!                     partial f(i+m1) / partial y(j+k*m2),  i,j=1,m2
!                    are banded, mljac is the maximal lower bandwidth
!                    of these mm+1 submatrices
!       - mujac: maximal upper bandwidth of these mm+1 submatrices
!                need not be defined if mljac=n-m1
!       - mas: if imas=0 this matrix is assumed to be the identity and
!              need not be defined. supply a dummy subroutine in this case.
!              it is assumed that only the elements of right lower block of
!              dimension n-m1 differ from that of the identity matrix.
!              if (mlmas.eq.n-m1) this submatrix is supposed to be full
!                 am(i,j) = m(i+m1,j+m1)     for i=1,n-m1 and j=1,n-m1.
!              else, the mass matrix is banded
!                 am(i-j+mumas+1,j) = m(i+m1,j+m1)
!       - mlmas: mlmas=n-m1: if the non-trivial part of m is full
!                0<=mlmas<n-m1: lower bandwidth of the mass matrix
!       - mumas: upper bandwidth of the mass matrix
!                need not be defined if mlmas=n-m1
!
!    iwork(9)  the value of m1.  default m1=0.
!
!    iwork(10) the value of m2.  default m2=m1.
!
! ----------
!
!    work(1)   uround, the rounding unit, default 1.d-16.
!
!    work(2)   the safety factor in step size prediction,
!              default 0.9d0.
!
!    work(3)   decides whether the jacobian should be recomputed;
!              increase work(3), to 0.1 say, when jacobian evaluations
!              are costly. for small systems work(3) should be smaller 
!              (0.001d0, say). negativ work(3) forces the code to
!              compute the jacobian after every accepted step.     
!              default 0.001d0.
!
!    work(4)   stopping criterion for newton's method, usually chosen <1.
!              smaller values of work(4) make the code slower, but safer.
!              default min(0.03d0,rtol(1)**0.5d0)
!
!    work(5) and work(6) : if work(5) < hnew/hold < work(6), then the
!              step size is not changed. this saves, together with a
!              large work(3), lu-decompositions and computing time for
!              large systems. for small systems one may have
!              work(5)=1.d0, work(6)=1.2d0, for large full systems
!              work(5)=0.99d0, work(6)=2.d0 might be good.
!              defaults work(5)=1.d0, work(6)=1.2d0 .
!
!    work(7)   maximal step size, default xend-x.
!
!    work(8), work(9)   parameters for step size selection
!              the new step size is chosen subject to the restriction
!                 work(8) <= hnew/hold <= work(9)
!              default values: work(8)=0.2d0, work(9)=8.d0
!
!-----------------------------------------------------------------------
!
!     output parameters 
!     ----------------- 
!     x           x-value for which the solution has been computed
!                 (after successful return x=xend).
!
!     y(n)        numerical solution at x
! 
!     h           predicted step size of the last accepted step
!
!     idid        reports on successfulness upon return:
!                   idid= 1  computation successful,
!                   idid= 2  comput. successful (interrupted by solout)
!                   idid=-1  input is not consistent,
!                   idid=-2  larger nmax is needed,
!                   idid=-3  step size becomes too small,
!                   idid=-4  matrix is repeatedly singular.
!
!   iwork(14)  nfcn    number of function evaluations (those for numerical
!                      evaluation of the jacobian are not counted)  
!   iwork(15)  njac    number of jacobian evaluations (either analytically
!                      or numerically)
!   iwork(16)  nstep   number of computed steps
!   iwork(17)  naccpt  number of accepted steps
!   iwork(18)  nrejct  number of rejected steps (due to error test),
!                      (step rejections in the first step are not counted)
!   iwork(19)  ndec    number of lu-decompositions of both matrices
!   iwork(20)  nsol    number of forward-backward substitutions, of both
!                      systems; the nstep forward-backward substitutions,
!                      needed for step size selection, are not counted
!-----------------------------------------------------------------------
! *** *** *** *** *** *** *** *** *** *** *** *** ***
!          declarations 
! *** *** *** *** *** *** *** *** *** *** *** *** ***
      implicit double precision (a-h,o-z)
      dimension y(n),atol(*),rtol(*),work(lwork),iwork(liwork)
      dimension rpar(*),ipar(*)
      logical implct,jband,arret,startn,pred
      external fcn,jac,mas,solout
! *** *** *** *** *** *** ***
!        setting the parameters 
! *** *** *** *** *** *** ***
       nfcn=0
       njac=0
       nstep=0
       naccpt=0
       nrejct=0
       ndec=0
       nsol=0
       arret=.false.
! -------- uround   smallest number satisfying 1.0d0+uround>1.0d0  
      if (work(1).eq.0.0d0) then
         uround=1.0d-16
      else
         uround=work(1)
         if (uround.le.1.0d-19.or.uround.ge.1.0d0) then
            write(6,*)' coefficients have 20 digits, uround=',work(1)
            arret=.true.
         end if
      end if
! -------- check and change the tolerances
      expm=2.0d0/3.0d0
      if (itol.eq.0) then
          if (atol(1).le.0.d0.or.rtol(1).le.10.d0*uround) then
              write (6,*) ' tolerances are too small'
              arret=.true.
          else
              quot=atol(1)/rtol(1)
              rtol(1)=0.1d0*rtol(1)**expm
              atol(1)=rtol(1)*quot
          end if
      else
          do i=1,n
          if (atol(i).le.0.d0.or.rtol(i).le.10.d0*uround) then
              write (6,*) ' tolerances(',i,') are too small'
              arret=.true.
          else
              quot=atol(i)/rtol(i)
              rtol(i)=0.1d0*rtol(i)**expm
              atol(i)=rtol(i)*quot
          end if
          end do
      end if
! -------- nmax , the maximal number of steps -----
      if (iwork(2).eq.0) then
         nmax=100000
      else
         nmax=iwork(2)
         if (nmax.le.0) then
            write(6,*)' wrong input iwork(2)=',iwork(2)
            arret=.true.
         end if
      end if
! -------- nit    maximal number of newton iterations
      if (iwork(3).eq.0) then
         nit=7
      else
         nit=iwork(3)
         if (nit.le.0) then
            write(6,*)' curious input iwork(3)=',iwork(3)
            arret=.true.
         end if
      end if
! -------- startn  switch for starting values of newton iterations
      if(iwork(4).eq.0)then
         startn=.false.
      else
         startn=.true.
      end if
! -------- parameter for differential-algebraic components
      nind1=iwork(5)
      nind2=iwork(6)
      nind3=iwork(7)
      if (nind1.eq.0) nind1=n
      if (nind1+nind2+nind3.ne.n) then
       write(6,*)' curious input for iwork(5,6,7)=',nind1,nind2,nind3
       arret=.true.
      end if
! -------- pred   step size control
      if(iwork(8).le.1)then
         pred=.true.
      else
         pred=.false.
      end if
! -------- parameter for second order equations
      m1=iwork(9)
      m2=iwork(10)
      nm1=n-m1
      if (m1.eq.0) m2=n
      if (m2.eq.0) m2=m1
      if (m1.lt.0.or.m2.lt.0.or.m1+m2.gt.n) then
       write(6,*)' curious input for iwork(9,10)=',m1,m2
       arret=.true.
      end if
! --------- safe     safety factor in step size prediction
      if (work(2).eq.0.0d0) then
         safe=0.9d0
      else
         safe=work(2)
         if (safe.le.0.001d0.or.safe.ge.1.0d0) then
            write(6,*)' curious input for work(2)=',work(2)
            arret=.true.
         end if
      end if
! ------ thet     decides whether the jacobian should be recomputed;
      if (work(3).eq.0.d0) then
         thet=0.001d0
      else
         thet=work(3)
         if (thet.ge.1.0d0) then
            write(6,*)' curious input for work(3)=',work(3)
            arret=.true.
         end if
      end if
! --- fnewt   stopping criterion for newton's method, usually chosen <1.
      tolst=rtol(1)
      if (work(4).eq.0.d0) then
         fnewt=max(10*uround/tolst,min(0.03d0,tolst**0.5d0))
      else
         fnewt=work(4)
         if (fnewt.le.uround/tolst) then
            write(6,*)' curious input for work(4)=',work(4)
            arret=.true.
         end if
      end if
! --- quot1 and quot2: if quot1 < hnew/hold < quot2, step size = const.
      if (work(5).eq.0.d0) then
         quot1=1.d0
      else
         quot1=work(5)
      end if
      if (work(6).eq.0.d0) then
         quot2=1.2d0
      else
         quot2=work(6)
      end if
      if (quot1.gt.1.0d0.or.quot2.lt.1.0d0) then
         write(6,*)' curious input for work(5,6)=',quot1,quot2
         arret=.true.
      end if
! -------- maximal step size
      if (work(7).eq.0.d0) then
         hmax=xend-x
      else
         hmax=work(7)
      end if 
! -------  facl,facr     parameters for step size selection
      if(work(8).eq.0.d0)then
         facl=5.d0
      else
         facl=1.d0/work(8)
      end if
      if(work(9).eq.0.d0)then
         facr=1.d0/8.0d0
      else
         facr=1.d0/work(9)
      end if
      if (facl.lt.1.0d0.or.facr.gt.1.0d0) then
            write(6,*)' curious input work(8,9)=',work(8),work(9)
            arret=.true.
         end if
! *** *** *** *** *** *** *** *** *** *** *** *** ***
!         computation of array entries
! *** *** *** *** *** *** *** *** *** *** *** *** ***
! ---- implicit, banded or not ?
      implct=imas.ne.0
      jband=mljac.lt.nm1
! -------- computation of the row-dimensions of the 2-arrays ---
! -- jacobian  and  matrices e1, e2
      if (jband) then
         ldjac=mljac+mujac+1
         lde1=mljac+ldjac
      else
         mljac=nm1
         mujac=nm1
         ldjac=nm1
         lde1=nm1
      end if
! -- mass matrix
      if (implct) then
          if (mlmas.ne.nm1) then
              ldmas=mlmas+mumas+1
              if (jband) then
                 ijob=4
              else
                 ijob=3
              end if
          else
              mumas=nm1
              ldmas=nm1
              ijob=5
          end if
! ------ bandwith of "mas" not smaller than bandwith of "jac"
          if (mlmas.gt.mljac.or.mumas.gt.mujac) then
             write (6,*) 'bandwith of "mas" not smaller than bandwith of
     & "jac"'
            arret=.true.
          end if
      else
          ldmas=0
          if (jband) then
             ijob=2
          else
             ijob=1
             if (n.gt.2.and.iwork(1).ne.0) ijob=7
          end if
      end if
      ldmas2=max(1,ldmas)
! ------ hessenberg option only for explicit equ. with full jacobian
      if ((implct.or.jband).and.ijob.eq.7) then
         write(6,*)' hessenberg option only for explicit equations with 
     &full jacobian'
         arret=.true.
      end if
! ------- prepare the entry-points for the arrays in work -----
      iez1=21
      iez2=iez1+n
      iez3=iez2+n
      iey0=iez3+n
      iescal=iey0+n
      ief1=iescal+n
      ief2=ief1+n
      ief3=ief2+n
      iecon=ief3+n
      iejac=iecon+4*n
      iemas=iejac+n*ldjac
      iee1=iemas+nm1*ldmas
      iee2r=iee1+nm1*lde1
      iee2i=iee2r+nm1*lde1
! ------ total storage requirement -----------
      istore=iee2i+nm1*lde1-1
      if(istore.gt.lwork)then
         write(6,*)' insufficient storage for work, min. lwork=',istore
         arret=.true.
      end if
! ------- entry points for integer workspace -----
      ieip1=21
      ieip2=ieip1+nm1
      ieiph=ieip2+nm1
! --------- total requirement ---------------
      istore=ieiph+nm1-1
      if (istore.gt.liwork) then
         write(6,*)' insuff. storage for iwork, min. liwork=',istore
         arret=.true.
      end if
! ------ when a fail has occured, we return with idid=-1
      if (arret) then
         idid=-1
         return
      end if
! -------- call to core integrator ------------
      call radcor(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,
     &   jac,ijac,mljac,mujac,mas,mlmas,mumas,solout,iout,idid,
     &   nmax,uround,safe,thet,fnewt,quot1,quot2,nit,ijob,startn,
     &   nind1,nind2,nind3,pred,facl,facr,m1,m2,nm1,
     &   implct,jband,ldjac,lde1,ldmas2,work(iez1),work(iez2),
     &   work(iez3),work(iey0),work(iescal),work(ief1),work(ief2),
     &   work(ief3),work(iejac),work(iee1),work(iee2r),work(iee2i),
     &   work(iemas),iwork(ieip1),iwork(ieip2),iwork(ieiph),
     &   work(iecon),nfcn,njac,nstep,naccpt,nrejct,ndec,nsol,rpar,ipar)
      iwork(14)=nfcn
      iwork(15)=njac
      iwork(16)=nstep
      iwork(17)=naccpt
      iwork(18)=nrejct
      iwork(19)=ndec
      iwork(20)=nsol
! -------- restore tolerances
      expm=1.0d0/expm
      if (itol.eq.0) then
              quot=atol(1)/rtol(1)
              rtol(1)=(10.0d0*rtol(1))**expm
              atol(1)=rtol(1)*quot
      else
          do i=1,n
              quot=atol(i)/rtol(i)
              rtol(i)=(10.0d0*rtol(i))**expm
              atol(i)=rtol(i)*quot
          end do
      end if
! ----------- return -----------
      return
      end
!
!     end of subroutine radau5
!
! ***********************************************************
!
      subroutine radcor(n,fcn,x,y,xend,hmax,h,rtol,atol,itol,
     &   jac,ijac,mljac,mujac,mas,mlmas,mumas,solout,iout,idid,
     &   nmax,uround,safe,thet,fnewt,quot1,quot2,nit,ijob,startn,
     &   nind1,nind2,nind3,pred,facl,facr,m1,m2,nm1,
     &   implct,banded,ldjac,lde1,ldmas,z1,z2,z3,
     &   y0,scal,f1,f2,f3,fjac,e1,e2r,e2i,fmas,ip1,ip2,iphes,
     &   cont,nfcn,njac,nstep,naccpt,nrejct,ndec,nsol,rpar,ipar)
! ----------------------------------------------------------
!     core integrator for radau5
!     parameters same as in radau5 with workspace added 
! ---------------------------------------------------------- 
!         declarations 
! ---------------------------------------------------------- 
      implicit double precision (a-h,o-z)
      dimension y(n),z1(n),z2(n),z3(n),y0(n),scal(n),f1(n),f2(n),f3(n)
      dimension fjac(ldjac,n),fmas(ldmas,nm1),cont(4*n)
      dimension e1(lde1,nm1),e2r(lde1,nm1),e2i(lde1,nm1)
      dimension atol(*),rtol(*),rpar(*),ipar(*)
      integer ip1(nm1),ip2(nm1),iphes(nm1)
      common /conra5/nn,nn2,nn3,nn4,xsol,hsol,c2m1,c1m1
      common/linal/mle,mue,mbjac,mbb,mdiag,mdiff,mbdiag
      logical reject,first,implct,banded,caljac,startn,calhes
      logical index1,index2,index3,last,pred
      external fcn
! *** *** *** *** *** *** ***
!  initialisations
! *** *** *** *** *** *** ***
! --------- duplify n for common block cont -----
      nn=n
      nn2=2*n
      nn3=3*n 
      lrc=4*n
! -------- check the index of the problem ----- 
      index1=nind1.ne.0
      index2=nind2.ne.0
      index3=nind3.ne.0
! ------- compute mass matrix for implicit case ----------
      if (implct) call mas(nm1,fmas,ldmas,rpar,ipar)
! ---------- constants ---------
      sq6=dsqrt(6.d0)
      c1=(4.d0-sq6)/10.d0
      c2=(4.d0+sq6)/10.d0
      c1m1=c1-1.d0
      c2m1=c2-1.d0
      c1mc2=c1-c2
      dd1=-(13.d0+7.d0*sq6)/3.d0
      dd2=(-13.d0+7.d0*sq6)/3.d0
      dd3=-1.d0/3.d0
      u1=(6.d0+81.d0**(1.d0/3.d0)-9.d0**(1.d0/3.d0))/30.d0
      alph=(12.d0-81.d0**(1.d0/3.d0)+9.d0**(1.d0/3.d0))/60.d0
      beta=(81.d0**(1.d0/3.d0)+9.d0**(1.d0/3.d0))*dsqrt(3.d0)/60.d0
      cno=alph**2+beta**2
      u1=1.0d0/u1
      alph=alph/cno
      beta=beta/cno
      t11=9.1232394870892942792d-02
      t12=-0.14125529502095420843d0
      t13=-3.0029194105147424492d-02
      t21=0.24171793270710701896d0
      t22=0.20412935229379993199d0
      t23=0.38294211275726193779d0
      t31=0.96604818261509293619d0
      ti11=4.3255798900631553510d0
      ti12=0.33919925181580986954d0
      ti13=0.54177053993587487119d0
      ti21=-4.1787185915519047273d0
      ti22=-0.32768282076106238708d0
      ti23=0.47662355450055045196d0
      ti31=-0.50287263494578687595d0
      ti32=2.5719269498556054292d0
      ti33=-0.59603920482822492497d0
      if (m1.gt.0) ijob=ijob+10
      posneg=sign(1.d0,xend-x)
      hmaxn=min(abs(hmax),abs(xend-x)) 
      if (abs(h).le.10.d0*uround) h=1.0d-6
      h=min(abs(h),hmaxn)
      h=sign(h,posneg)
      hold=h
      reject=.false.
      first=.true.
      last=.false.
      if ((x+h*1.0001d0-xend)*posneg.ge.0.d0) then
         h=xend-x
         last=.true.
      end if
      hopt=h
      faccon=1.d0
      cfac=safe*(1+2*nit)
      nsing=0
      xold=x
      if (iout.ne.0) then
          irtrn=1
          nrsol=1
          xosol=xold
          xsol=x
          do i=1,n
             cont(i)=y(i)
          end do
          nsolu=n
          hsol=hold
          call solout(nrsol,xosol,xsol,y,cont,lrc,nsolu,
     &                rpar,ipar,irtrn)
          if (irtrn.lt.0) goto 179
      end if
      mle=mljac
      mue=mujac
      mbjac=mljac+mujac+1
      mbb=mlmas+mumas+1
      mdiag=mle+mue+1
      mdiff=mle+mue-mumas
      mbdiag=mumas+1
      n2=2*n
      n3=3*n
      if (itol.eq.0) then
          do i=1,n
             scal(i)=atol(1)+rtol(1)*abs(y(i))
          end do
      else
          do i=1,n
             scal(i)=atol(i)+rtol(i)*abs(y(i))
          end do
      end if
      hhfac=h
      call fcn(n,x,y,y0,rpar,ipar)
      nfcn=nfcn+1
! --- basic integration step  
  10  continue
! *** *** *** *** *** *** ***
!  computation of the jacobian
! *** *** *** *** *** *** ***
      njac=njac+1
      if (ijac.eq.0) then
! --- compute jacobian matrix numerically
         if (banded) then
! --- jacobian is banded
            mujacp=mujac+1
            md=min(mbjac,m2)
            do mm=1,m1/m2+1
               do k=1,md
                  j=k+(mm-1)*m2
 12               f1(j)=y(j)
                  f2(j)=dsqrt(uround*max(1.d-5,abs(y(j))))
                  y(j)=y(j)+f2(j)
                  j=j+md
                  if (j.le.mm*m2) goto 12 
                  call fcn(n,x,y,cont,rpar,ipar)
                  j=k+(mm-1)*m2
                  j1=k
                  lbeg=max(1,j1-mujac)+m1
 14               lend=min(m2,j1+mljac)+m1
                  y(j)=f1(j)
                  mujacj=mujacp-j1-m1
                  do l=lbeg,lend
                     fjac(l+mujacj,j)=(cont(l)-y0(l))/f2(j) 
                  end do
                  j=j+md
                  j1=j1+md
                  lbeg=lend+1
                  if (j.le.mm*m2) goto 14
               end do
            end do
         else
! --- jacobian is full
            do i=1,n
               ysafe=y(i)
               delt=dsqrt(uround*max(1.d-5,abs(ysafe)))
               y(i)=ysafe+delt
               call fcn(n,x,y,cont,rpar,ipar)
               do j=m1+1,n
                 fjac(j-m1,i)=(cont(j)-y0(j))/delt
               end do
               y(i)=ysafe
            end do
         end if
      else
! --- compute jacobian matrix analytically
         call jac(n,x,y,fjac,ldjac,rpar,ipar)
      end if
      caljac=.true.
      calhes=.true.
  20  continue
! --- compute the matrices e1 and e2 and their decompositions
      fac1=u1/h
      alphn=alph/h
      betan=beta/h
      call decomr(n,fjac,ldjac,fmas,ldmas,mlmas,mumas,
     &            m1,m2,nm1,fac1,e1,lde1,ip1,ier,ijob,calhes,iphes)
      if (ier.ne.0) goto 78
      call decomc(n,fjac,ldjac,fmas,ldmas,mlmas,mumas,
     &            m1,m2,nm1,alphn,betan,e2r,e2i,lde1,ip2,ier,ijob)
      if (ier.ne.0) goto 78
      ndec=ndec+1
  30  continue
      nstep=nstep+1
      if (nstep.gt.nmax) goto 178
      if (0.1d0*abs(h).le.abs(x)*uround) goto 177
          if (index2) then
             do i=nind1+1,nind1+nind2
                scal(i)=scal(i)/hhfac
             end do
          end if
          if (index3) then
             do i=nind1+nind2+1,nind1+nind2+nind3
                scal(i)=scal(i)/(hhfac*hhfac)
             end do
          end if
      xph=x+h
! *** *** *** *** *** *** ***
!  starting values for newton iteration
! *** *** *** *** *** *** ***
      if (first.or.startn) then
         do i=1,n
            z1(i)=0.d0
            z2(i)=0.d0
            z3(i)=0.d0
            f1(i)=0.d0
            f2(i)=0.d0
            f3(i)=0.d0
         end do
      else
         c3q=h/hold
         c1q=c1*c3q
         c2q=c2*c3q
         do i=1,n
            ak1=cont(i+n)
            ak2=cont(i+n2)
            ak3=cont(i+n3)
            z1i=c1q*(ak1+(c1q-c2m1)*(ak2+(c1q-c1m1)*ak3))
            z2i=c2q*(ak1+(c2q-c2m1)*(ak2+(c2q-c1m1)*ak3))
            z3i=c3q*(ak1+(c3q-c2m1)*(ak2+(c3q-c1m1)*ak3))
            z1(i)=z1i
            z2(i)=z2i
            z3(i)=z3i
            f1(i)=ti11*z1i+ti12*z2i+ti13*z3i
            f2(i)=ti21*z1i+ti22*z2i+ti23*z3i
            f3(i)=ti31*z1i+ti32*z2i+ti33*z3i
         end do
      end if
! *** *** *** *** *** *** ***
!  loop for the simplified newton iteration
! *** *** *** *** *** *** ***
            newt=0
            faccon=max(faccon,uround)**0.8d0
            theta=abs(thet)
  40        continue
            if (newt.ge.nit) goto 78
! ---     compute the right-hand side
            do i=1,n
               cont(i)=y(i)+z1(i)
            end do
            call fcn(n,x+c1*h,cont,z1,rpar,ipar)
            do i=1,n
               cont(i)=y(i)+z2(i)
            end do
            call fcn(n,x+c2*h,cont,z2,rpar,ipar)
            do i=1,n
               cont(i)=y(i)+z3(i)
            end do
            call fcn(n,xph,cont,z3,rpar,ipar)
            nfcn=nfcn+3
! ---     solve the linear systems
           do i=1,n
              a1=z1(i)
              a2=z2(i)
              a3=z3(i)
              z1(i)=ti11*a1+ti12*a2+ti13*a3
              z2(i)=ti21*a1+ti22*a2+ti23*a3
              z3(i)=ti31*a1+ti32*a2+ti33*a3
           end do
        call slvrad(n,fjac,ldjac,mljac,mujac,fmas,ldmas,mlmas,mumas,
     &          m1,m2,nm1,fac1,alphn,betan,e1,e2r,e2i,lde1,z1,z2,z3,
     &          f1,f2,f3,cont,ip1,ip2,iphes,ier,ijob)
            nsol=nsol+1
            newt=newt+1
            dyno=0.d0
            do i=1,n
               denom=scal(i)
               dyno=dyno+(z1(i)/denom)**2+(z2(i)/denom)**2
     &          +(z3(i)/denom)**2
            end do
            dyno=dsqrt(dyno/n3)
! ---     bad convergence or number of iterations to large
            if (newt.gt.1.and.newt.lt.nit) then
                thq=dyno/dynold
                if (newt.eq.2) then
                   theta=thq
                else
                   theta=sqrt(thq*thqold)
                end if
                thqold=thq
                if (theta.lt.0.99d0) then
                    faccon=theta/(1.0d0-theta)
                    dyth=faccon*dyno*theta**(nit-1-newt)/fnewt
                    if (dyth.ge.1.0d0) then
                         qnewt=dmax1(1.0d-4,dmin1(20.0d0,dyth))
                         hhfac=.8d0*qnewt**(-1.0d0/(4.0d0+nit-1-newt))
                         h=hhfac*h
                         reject=.true.
                         last=.false.
                         if (caljac) goto 20
                         goto 10
                    end if
                else
                    goto 78
                end if
            end if
            dynold=max(dyno,uround)
            do i=1,n
               f1i=f1(i)+z1(i)
               f2i=f2(i)+z2(i)
               f3i=f3(i)+z3(i)
               f1(i)=f1i
               f2(i)=f2i
               f3(i)=f3i
               z1(i)=t11*f1i+t12*f2i+t13*f3i
               z2(i)=t21*f1i+t22*f2i+t23*f3i
               z3(i)=t31*f1i+    f2i
            end do
            if (faccon*dyno.gt.fnewt) goto 40
! --- error estimation  
      call estrad (n,fjac,ldjac,mljac,mujac,fmas,ldmas,mlmas,mumas,
     &          h,dd1,dd2,dd3,fcn,nfcn,y0,y,ijob,x,m1,m2,nm1,
     &          e1,lde1,z1,z2,z3,cont,f1,f2,ip1,iphes,scal,err,
     &          first,reject,fac1,rpar,ipar)
! --- computation of hnew
! --- we require .2<=hnew/h<=8.
      fac=min(safe,cfac/(newt+2*nit))
      quot=max(facr,min(facl,err**.25d0/fac))
      hnew=h/quot
! *** *** *** *** *** *** ***
!  is the error small enough ?
! *** *** *** *** *** *** ***
      if (err.lt.1.d0) then
! --- step is accepted  
         first=.false.
         naccpt=naccpt+1
         if (pred) then
!       --- predictive controller of gustafsson
            if (naccpt.gt.1) then
               facgus=(hacc/h)*(err**2/erracc)**0.25d0/safe
               facgus=max(facr,min(facl,facgus))
               quot=max(quot,facgus)
               hnew=h/quot
            end if
            hacc=h
            erracc=max(1.0d-2,err)
         end if
         xold=x
         hold=h
         x=xph 
         do i=1,n
            y(i)=y(i)+z3(i)  
            z2i=z2(i)
            z1i=z1(i)
            cont(i+n)=(z2i-z3(i))/c2m1
            ak=(z1i-z2i)/c1mc2
            acont3=z1i/c1
            acont3=(ak-acont3)/c2
            cont(i+n2)=(ak-cont(i+n))/c1m1
            cont(i+n3)=cont(i+n2)-acont3
         end do
         if (itol.eq.0) then
             do i=1,n
                scal(i)=atol(1)+rtol(1)*abs(y(i))
             end do
         else
             do i=1,n
                scal(i)=atol(i)+rtol(i)*abs(y(i))
             end do
         end if
         if (iout.ne.0) then
             nrsol=naccpt+1
             xsol=x
             xosol=xold
             do i=1,n
                cont(i)=y(i)
             end do
             nsolu=n
             hsol=hold
             call solout(nrsol,xosol,xsol,y,cont,lrc,nsolu,
     &                   rpar,ipar,irtrn)
             if (irtrn.lt.0) goto 179
         end if
         caljac=.false.
         if (last) then
            h=hopt
            idid=1
            return
         end if
         call fcn(n,x,y,y0,rpar,ipar)
         nfcn=nfcn+1
         hnew=posneg*min(abs(hnew),hmaxn)
         hopt=hnew
         hopt=min(h,hnew)
         if (reject) hnew=posneg*min(abs(hnew),abs(h)) 
         reject=.false.
         if ((x+hnew/quot1-xend)*posneg.ge.0.d0) then
            h=xend-x
            last=.true.
         else
            qt=hnew/h 
            hhfac=h
            if (theta.le.thet.and.qt.ge.quot1.and.qt.le.quot2) goto 30
            h=hnew 
         end if
         hhfac=h
         if (theta.le.thet) goto 20
         goto 10
      else
! --- step is rejected  
         reject=.true.
         last=.false.
         if (first) then
             h=h*0.1d0
             hhfac=0.1d0
         else 
             hhfac=hnew/h
             h=hnew
         end if
         if (naccpt.ge.1) nrejct=nrejct+1
         if (caljac) goto 20
         goto 10
      end if
! --- unexpected step-rejection
  78  continue
      if (ier.ne.0) then
          nsing=nsing+1
          if (nsing.ge.5) goto 176
      end if
      h=h*0.5d0 
      hhfac=0.5d0
      reject=.true.
      last=.false.
      if (caljac) goto 20
      goto 10
! --- fail exit
 176  continue
      write(6,979)x   
      write(6,*) ' matrix is repeatedly singular, ier=',ier
      idid=-4
      return
 177  continue
      write(6,979)x   
      write(6,*) ' step size t0o small, h=',h
      idid=-3
      return
 178  continue
      write(6,979)x   
      write(6,*) ' more than nmax =',nmax,'steps are needed' 
      idid=-2
      return
! --- exit caused by solout
 179  continue
      write(6,979)x
 979  format(' exit of radau5 at x=',e18.4) 
      idid=2
      return
      end
!
!     end of subroutine radcor
!
! ***********************************************************
!
      double precision function contr5(i,x,cont,lrc) 
! ----------------------------------------------------------
!     this function can be used for coninuous output. it provides an
!     approximation to the i-th component of the solution at x.
!     it gives the value of the collocation polynomial, defined for
!     the last successfully computed step (by radau5).
! ----------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension cont(lrc)
      common /conra5/nn,nn2,nn3,nn4,xsol,hsol,c2m1,c1m1
      s=(x-xsol)/hsol
      contr5=cont(i)+s*(cont(i+nn)+(s-c2m1)*(cont(i+nn2)
     &     +(s-c1m1)*cont(i+nn3)))
      return
      end
!
!     end of function contr5
!
! ***********************************************************
