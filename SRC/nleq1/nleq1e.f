      SUBROUTINE NLEQ1E(N,X,RTOL,IERR)
C*    Begin Prologue NLEQ1E
      INTEGER N
      DOUBLE PRECISION X(N), RTOL
      INTEGER IERR
C     ------------------------------------------------------------
C
C*  Title
C
C     Numerical solution of nonlinear (NL) equations (EQ)
C     especially designed for numerically sensitive problems.
C     (E)asy-to-use driver routine for NLEQ1.
C
C*  Written by        U. Nowak, L. Weimann 
C*  Purpose           Solution of systems of highly nonlinear equations
C*  Method            Damped affine invariant Newton method
C                     (see references below)
C*  Category          F2a. - Systems of nonlinear equations
C*  Keywords          Nonlinear equations, Newton methods
C*  Version           2.3
C*  Revision          November 1991
C*  Latest Change     November 1991
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C*  Copyright     (c) Konrad-Zuse-Zentrum fuer
C                     Informationstechnik Berlin (ZIB)
C                     Takustrasse 7, D-14195 Berlin-Dahlem
C                     phone : + 49/30/84185-0
C                     fax   : + 49/30/84185-125
C*  Contact           Bodo Erdmann
C                     ZIB, Division Scientific Computing, 
C                          Department Numerical Analysis and Modelling
C                     phone : + 49/30/84185-185
C                     fax   : + 49/30/84185-107
C                     e-mail: erdmann@zib.de
C
C*    References:
C
C     /1/ P. Deuflhard:
C         Newton Methods for Nonlinear Problems. -
C         Affine Invariance and Adaptive Algorithms.
C         Series Computational Mathematics 35, Springer (2004)
C
C     /2/ U. Nowak, L. Weimann:
C         A Family of Newton Codes for Systems of Highly Nonlinear
C         Equations - Algorithm, Implementation, Application.
C         ZIB, Technical Report TR 90-10 (December 1990)
C
C  ---------------------------------------------------------------
C
C* Licence
C    You may use or modify this code for your own non commercial
C    purposes for an unlimited time. 
C    In any case you should not deliver this code without a special 
C    permission of ZIB.
C    In case you intend to use the code commercially, we oblige you
C    to sign an according licence agreement with ZIB.
C
C* Warranty 
C    This code has been tested up to a certain level. Defects and
C    weaknesses, which may be included in the code, do not establish
C    any warranties by ZIB. ZIB does not take over any liabilities
C    which may follow from aquisition or application of this code.
C
C* Software status 
C    This code is under care of ZIB and belongs to ZIB software class 1.
C
C     ------------------------------------------------------------
C
C*    Summary:
C     ========
C     Damped Newton-algorithm for systems of highly nonlinear
C     equations - damping strategy due to Ref. (1).
C
C     (The iteration is done by subroutine N1INT actually. NLEQ1E
C      calls the standard interface driver NLEQ1, which itself does 
C      some house keeping and builds up workspace.)
C
C     Jacobian approximation by numerical differences.
C
C     The numerical solution of the arising linear equations is
C     done by means of the subroutines *GEFA and *GESL ( Gauss-
C     algorithm with column-pivoting and row-interchange; 
C     replace '*' by 'S' or 'D' for single or double precision
C     version respectively).
C
C     ------------------------------------------------------------
C
C*    External subroutine to be supplied by the user
C     ==============================================
C 
C     FCN(N,X,F,IFAIL) Ext    Function subroutine - the problem function
C                             (The routine must be named exactly FCN !)
C       N              Int    Number of vector components (input)
C                             Must not be altered!
C       X(N)           Dble   Vector of unknowns (input)
C                             Must not be altered!
C       F(N)           Dble   Vector of function values (output)
C       IFAIL          Int    FCN evaluation-failure indicator. (output)
C                             On input:  Has always value 0 (zero).
C                             On output: Indicates failure of FCN eval-
C                                uation, if having a nonzero value.
C                             If <0: NLEQ1E will be terminated with 
C                                    IFAIL returned via IERR.
C                             If =1: A new trial Newton iterate will be
C                                    computed, with the damping factor
C                                    reduced to it's half.
C                             If =2: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced by a reduction factor,
C                                    which must be output through F(1)
C                                    by FCN, and it's value must be
C                                    >0 and < 1.
C
C*    Parameters list description
C     ===========================
C
C     Name   Type    In/Out Meaning
C
C     N      Int     In     Number of unknowns ( N .LE. 50 )
C     X(N)   Dble    In     Initial estimate of the solution
C                    Out    Solution values ( or final values,
C                           respectively )
C     RTOL   Dble    In     Required relative precision of
C                           solution components -
C                           RTOL.GE.EPMACH*TEN*N
C                    Out    Finally achieved accuracy
C     IERR   Int     Out    Return value parameter
C                           < 0 Termination forced by user function FCN
C                               due to IFAIL<0 on output, IERR is set
C                               to IFAIL 
C                           = 0 successfull completion of the iteration,
C                               solution has been computed
C                           > 0 see list of error/warning messages below
C
C*   Error and warning messages:
C    ===========================
C
C      1    Termination, since Jacobian matrix became singular
C      2    Termination after 100 iterations
C      3    Termination, since damping factor became to small
C      4    Warning: Superlinear or quadratic convergence slowed down
C           near the solution.
C           Iteration has been stopped therefore with an approximation
C           of the solution not such accurate as requested by RTOL,
C           because possibly the RTOL requirement may be too stringent
C           (i.e. the nonlinear problem is ill-conditioned)
C      5    Warning: Iteration stopped with termination criterion 
C           (using RTOL as requested precision) satisfied, but no 
C           superlinear or quadratic convergence has been indicated yet.
C           Therefore, possibly the error estimate for the solution may
C           not match good enough the really achieved accuracy.
C     20    Bad input value to parameter N; 1 .LE. N .LE. 50 required
C     21    Nonpositive value for RTOL supplied
C     82    Termination, because user routine FCN returned with IFAIL>0
C
C     Note   : in case of failure:
C        -    use better initial guess
C        -    or refine model
C        -    or use non-standard options and/or analytical Jacobian
C             via the standard interface NLEQ1
C
C*    Machine dependent constants used:
C     =================================
C
C     DOUBLE PRECISION EPMACH  in  N1PCHK, N1INT
C     DOUBLE PRECISION GREAT   in  N1PCHK
C     DOUBLE PRECISION SMALL   in  N1PCHK, N1INT, N1SCAL
C
C*    Subroutines called: NLEQ1
C
C     ------------------------------------------------------------
C*    End Prologue
      INTEGER NMAX, LIOPT, LIWK, LRWK, LUPRT
      PARAMETER ( NMAX=50 )
      PARAMETER ( LIOPT=50, LIWK=NMAX+50, LRWK=(NMAX+13)*NMAX+60 )
      PARAMETER ( LUPRT=6 )
      EXTERNAL FCN, NLEQ1
      DOUBLE PRECISION XSCAL(NMAX)
      INTEGER IOPT(LIOPT), IWK(LIWK)
      DOUBLE PRECISION RWK(LRWK)
      INTEGER I, NIW, NRW
      CHARACTER CHGDAT*20, PRODCT*8
C
C     Version: 2.3               Latest change:
C     -----------------------------------------
C
      DATA      CHGDAT      /'November 15, 1991   '/
      DATA      PRODCT      /'NLEQ1E  '/
C*    Begin
      NIW = N+50
      NRW = (N+13)*N+60
C     Checking dimensional parameter N
      IF ( N.LT.1 .OR. N.GT.NMAX ) THEN
        WRITE(LUPRT,1001) NMAX, N
1001    FORMAT(/,' Error: Bad input to parameter N supplied',
     $         /,8X,'choose 1 .LE. N .LE. ',I3,
     $               ' , your input is: N = ',I5)
        IERR = 20
        RETURN
      ENDIF
      DO 10 I=1,LIOPT
        IOPT(I) = 0
10    CONTINUE
      DO 20 I=1,NIW
        IWK(I) = 0
20    CONTINUE
      DO 30 I=1,NRW
        RWK(I) = 0.0D0
30    CONTINUE
      DO 40 I=1,N
        XSCAL(I) = 0.0D0
40    CONTINUE
C     Print errors, warnings, monitor and time monitor
C     to standard output
      IOPT(11) = 3
      IOPT(12) = LUPRT
      IOPT(13) = 3
      IOPT(14) = LUPRT
      IOPT(19) = 1
      IOPT(20) = LUPRT
C     Maximum number of Newton iterations
      IWK(31) = 100
C
      CALL NLEQ1(N,FCN,DUMMY,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
      IF (IERR.EQ.82 .AND. IWK(23).LT.0) IERR=IWK(23)
C     End of subroutine NLEQ1E
      RETURN
      END
