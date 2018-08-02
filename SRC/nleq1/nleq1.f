      SUBROUTINE NLEQ1(N,FCN,JAC,X,XSCAL,RTOL,IOPT,IERR,
     $LIWK,IWK,LRWK,RWK)
C*    Begin Prologue NLEQ1
      INTEGER N
      EXTERNAL FCN,JAC
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*  Title
C
C     Numerical solution of nonlinear (NL) equations (EQ)
C     especially designed for numerically sensitive problems.
C
C*  Written by        U. Nowak, L. Weimann 
C*  Purpose           Solution of systems of highly nonlinear equations
C*  Method            Damped affine invariant Newton method
C                     (see references below)
C*  Category          F2a. - Systems of nonlinear equations
C*  Keywords          Nonlinear equations, Newton methods
C*  Version           2.4
C*  Revision          May 2009
C*  Latest Change     May 2009
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
C    which may follow from acquisition or application of this code.
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
C     (The iteration is done by subroutine N1INT currently. NLEQ1
C      itself does some house keeping and builds up workspace.)
C
C     Jacobian approximation by numerical differences or user
C     supplied subroutine JAC.
C
C     The numerical solution of the arising linear equations is
C     done by means of the subroutines *GETRF and *GETRS ( Gauss-
C     algorithm with column-pivoting and row-interchange ) in the
C     dense matrix case, or by the subroutines *GBTRF and *GBTRS in
C     the band matrix case from LAPACK (replace '*' by 'S' or 'D'
C     for single or double precision version respectively).
C     For special purposes these routines may be substituted.
C
C     This is a driver routine for the core solver N1INT.
C
C     ------------------------------------------------------------
C
C*    Parameters list description (* marks inout parameters)
C     ======================================================
C
C*    External subroutines (to be supplied by the user)
C     =================================================
C 
C     (Caution: Arguments declared as (input) must not
C               be altered by the user subroutines ! )
C
C     FCN(N,X,F,IFAIL) Ext    Function subroutine
C       N              Int    Number of vector components (input)
C       X(N)           Dble   Vector of unknowns (input)
C       F(N)           Dble   Vector of function values (output)
C       IFAIL          Int    FCN evaluation-failure indicator. (output)
C                             On input:  Has always value 0 (zero).
C                             On output: Indicates failure of FCN eval-
C                                uation, if having a value <= 2.
C                             If <0: NLEQ1 will be terminated with 
C                                    error code = 82, and IFAIL stored
C                                    to IWK(23).
C                             If =1: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced to it's half.
C                             If =2: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced by a reduct. factor, which
C                                    must be output through F(1) by FCN,
C                                    and it's value must be >0 and < 1.
C                             Note, that if IFAIL = 1 or 2, additional
C                             conditions concerning the damping factor,
C                             e.g. the minimum damping factor or the
C                             bounded damping strategy may also influ-
C                             ence the value of the reduced damping 
C                             factor.
C
C     JAC(N,LDJAC,X,DFDX,IFAIL) 
C                       Ext    Jacobian matrix subroutine
C       N                 Int    Number of vector components (input)
C       LDJAC             Int    Leading dimension of Jacobian array 
C                                (input)
C       X(N)              Dble   Vector of unknowns (input)
C       DFDX(LDJAC,N)     Dble   DFDX(i,k): partial derivative of
C                                I-th component of FCN with respect
C                                to X(k) (output)
C       IFAIL             Int    JAC evaluation-failure indicator.
C                                (output)
C                                Has always value 0 (zero) on input.
C                                Indicates failure of JAC evaluation
C                                and causes termination of NLEQ1,
C                                if set to a negative value on output
C
C
C*    Input parameters of NLEQ1
C     =========================
C
C     N              Int    Number of unknowns
C   * X(N)           Dble   Initial estimate of the solution
C   * XSCAL(N)       Dble   User scaling (lower threshold) of the 
C                           iteration vector X(N)
C   * RTOL           Dble   Required relative precision of
C                           solution components -
C                           RTOL.GE.EPMACH*TEN*N
C   * IOPT(50)       Int    Array of run-time options. Set to zero
C                           to get default values (details see below)
C
C*    Output parameters of NLEQ1
C     ==========================
C
C   * X(N)           Dble   Solution values ( or final values,
C                           respectively )
C   * XSCAL(N)       Dble   After return with IERR.GE.0, it contains
C                           the latest internal scaling vector used
C                           After return with IERR.EQ.-1 in onestep-
C                           mode it contains a possibly adapted 
C                           (as described below) user scaling vector:
C                           If (XSCAL(I).LT. SMALL) XSCAL(I) = SMALL ,
C                           If (XSCAL(I).GT. GREAT) XSCAL(I) = GREAT .
C                           For SMALL and GREAT, see section machine
C                           constants below  and regard note 1.
C   * RTOL           Dble   Finally achieved (relative) accuracy.
C                           The estimated absolute error of component i
C                           of x_out is approximately given by
C                             abs_err(i) = RTOL * XSCAL_out(i) ,
C                           where (approximately)
C                             XSCAL_out(i) = 
C                                max(abs(X_out(i)),XSCAL_in(i)).
C     IERR           Int    Return value parameter
C                           =-1 sucessfull completion of one iteration
C                               step, subsequent iterations are needed 
C                               to get a solution. (stepwise mode only) 
C                           = 0 successfull completion of the iteration,
C                               solution has been computed
C                           > 0 see list of error messages below
C
C     Note 1.
C        The machine dependent values SMALL and EPMACH are
C        gained from calls of the ZIBCONST.
C
C*    Workspace parameters of NLEQ1
C     =============================
C
C     LIWK           Int    Declared dimension of integer workspace.
C                           Required minimum (for standard linear system
C                           solver) : N+50
C  *  IWK(LIWK)      Int    Integer Workspace
C     LRWK           Int    Declared dimension of real workspace.
C                           Required minimum (for standard linear system
C                           solver and Jacobian computed by numerical
C                           approximation - if the Jacobian is computed
C                           by a user subroutine JAC, decrease the 
C                           expressions noted below by N):
C                           for full case Jacobian: (N+NBROY+13)*N+61
C                           for a band-matrix Jacobian:
C                              (2*ML+MU+NBROY+14)*N+61 with
C                           ML = lower bandwidth , MU = upper bandwidth
C                           NBROY = Maximum number of Broyden steps
C                           (Default: if Broyden steps are enabled, e.g.
C                                                IOPT(32)=1            -
C                                       NBROY=N (full Jacobian), 
C                                            =ML+MU+1 (band Jacobian),
C                                       but at least NBROY=10 
C                                     else (if IOPT(32)=0) - 
C                                       NBROY=0 ;
C                            see equally named IOPT and IWK-fields below)
C   * RWK(LRWK)      Dble   Real Workspace
C
C     Note 2a.  A test on sufficient workspace is made. If this
C               test fails, IERR is set to 10 and an error-message
C               is issued from which the minimum of required
C               workspace size can be obtained.
C
C     Note 2b.  The first 50 elements of IWK and RWK are partially 
C               used as input for internal algorithm parameters (for
C               details, see below). In order to set the default values
C               of these parameters, the fields must be set to zero.
C               Therefore, it's recommanded always to initialize the
C               first 50 elements of both workspaces to zero.
C
C*   Options IOPT:
C    =============
C
C     Pos. Name   Default  Meaning
C
C       1  QSUCC  0        =0 (.FALSE.) initial call:
C                             NLEQ1 is not yet initialized, i.e. this is
C                             the first call for this nonlinear system.
C                             At successfull return with MODE=1,
C                             QSUCC is set to 1.
C                          =1 (.TRUE.) successive call:
C                             NLEQ1 is initialized already and is now
C                             called to perform one or more following
C                             Newton-iteration steps.
C                             ATTENTION:
C                                Don't destroy the contents of
C                                IOPT(i) for 1 <= i <= 50 ,
C                                IWK(j)  for 1 <= j < NIWKFR and
C                                RWK(k)  for 1 <= k < NRWKFR.
C                                (Nevertheless, some of the options, e.g.
C                                 FCMIN, SIGMA, MPR..., can be modified
C                                 before successive calls.)
C       2  MODE   0        =0 Standard mode initial call:
C                             Return when the required accuracy for the
C                             iteration vector is reached. User defined
C                             parameters are evaluated and checked.
C                             Standard mode successive call:
C                             If NLEQ1 was called previously with MODE=1,
C                             it performs all remaining iteration steps.
C                          =1 Stepwise mode:
C                             Return after one Newton iteration step.
C       3  JACGEN 0        Method of Jacobian generation
C                          =0 Standard method is JACGEN=2
C                          =1 User supplied subroutine JAC will be 
C                             called to generate Jacobian matrix
C                          =2 Jacobian approximation by numerical
C                             differentation (see subroutines N1JAC
C                             and N1JACB)
C                          =3 Jacobian approximation by numerical
C                             differentation with feedback control
C                             (see subroutines N1JCF and N1JCFB)
C       4  MSTOR  0        =0 The Jacobian A is a dense matrix
C                          =1 A is a band matrix
C       5                  Reserved
C       6  ML     0        Lower bandwidth of A (excluding the
C                          diagonal);
C                          IOPT(6) ignored, if IOPT(4).NE.1
C       7  MU     0        Upper bandwidth of A (excluding the
C                          diagonal);
C                          IOPT(7) ignored, if IOPT(4).NE.1
C       8                  Reserved
C       9  ISCAL  0        Determines how to scale the iterate-vector:
C                          =0 The user supplied scaling vector XSCAL is
C                             used as a (componentwise) lower threshold
C                             of the current scaling vector
C                          =1 The vector XSCAL is always used as the
C                             current scaling vector
C      10                  Reserved
C      11  MPRERR 0        Print error messages
C                          =0 No output
C                          =1 Error messages
C                          =2 Warnings additionally
C                          =3 Informal messages additionally
C      12  LUERR  6        Logical unit number for error messages
C      13  MPRMON 0        Print iteration Monitor
C                          =0 No output
C                          =1 Standard output
C                          =2 Summary iteration monitor additionally
C                          =3 Detailed iteration monitor additionally
C                          =4,5,6 Outputs with increasing level addi-
C                             tional increasing information for code
C                             testing purposes. Level 6 produces
C                             in general extremely large output!
C      14  LUMON  6        Logical unit number for iteration monitor
C      15  MPRSOL 0        Print solutions
C                          =0 No output
C                          =1 Initial values and solution values
C                          =2 Intermediate iterates additionally
C      16  LUSOL  6        Logical unit number for solutions
C      17..18              Reserved
C      19  MPRTIM 0        Output level for the time monitor
C                          = 0 : no time measurement and no output
C                          = 1 : time measurement will be done and
C                                summary output will be written -
C                                regard note 5a.
C      20  LUTIM  6        Logical output unit for time monitor
C      21..30              Reserved
C      31  NONLIN 3        Problem type specification
C                          =1 Linear problem
C                             Warning: If specified, no check will be
C                             done, if the problem is really linear, and
C                             NLEQ1 terminates unconditionally after one
C                             Newton-iteration step.
C                          =2 Mildly nonlinear problem
C                          =3 Highly nonlinear problem
C                          =4 Extremely nonlinear problem
C      32  QRANK1 0        =0 (.FALSE.) Rank-1 updates by Broyden-
C                             approximation are inhibited.
C                          =1 (.TRUE.) Rank-1 updates by Broyden-
C                             approximation are allowed.
C      33  QORDI  0        =0 (.FALSE.) Standard program mode 
C                          =1 (.TRUE.)  Special program mode:
C                             Ordinary Newton iteration is done, e.g.:
C                             No damping strategy and no monotonicity
C                             test is applied
C      34  QSIMPL 0        =0 (.FALSE.) Standard program mode
C                          =1 (.TRUE.)  Special program mode:
C                             Simplified Newton iteration is done, e.g.:
C                             The Jacobian computed at the starting
C                             point is fixed for all subsequent
C                             iteration steps, and
C                             no damping strategy and no monotonicity
C                             test is applied.
C      35  QNSCAL 0        Inhibit automatic row scaling: 
C                          =0 (.FALSE.) Automatic row scaling of
C                             the linear system is activ: 
C                             Rows i=1,...,N will be divided by
C                             max j=1,...,N (abs(a(i,j))) 
C                          =1 (.TRUE.) No row scaling of the linear
C                             system. Recommended only for well row-
C                             scaled nonlinear systems.
C      36..37              Reserved
C      38  IBDAMP          Bounded damping strategy switch:
C                          =0 The default switch takes place, dependent
C                             on the setting of NONLIN (=IOPT(31)):
C                             NONLIN = 0,1,2,3 -> IBDAMP = off ,
C                             NONLIN = 4 -> IBDAMP = on
C                          =1 means always IBDAMP = on 
C                          =2 means always IBDAMP = off
C      39  IORMON          Convergence order monitor 
C                          =0 Standard option is IORMON=2 
C                          =1 Convergence order is not checked,
C                             the iteration will be always proceeded
C                             until the solution has the required 
C                             precision RTOL (or some error condition
C                             occured)
C                          =2 Use additional 'weak stop' criterion:
C                             Convergence order is monitored
C                             and termination due to slowdown of the
C                             convergence may occur.
C                          =3 Use additional 'hard stop' criterion:
C                             Convergence order is monitored
C                             and termination due to superlinear 
C                             convergence slowdown may occur. 
C                          In case of termination due to convergence
C                          slowdown, the warning code IERR=4 will be
C                          set.
C                          In cases, where the Newton iteration con-
C                          verges but superlinear convergence order has
C                          never been detected, the warning code IERR=5 
C                          is returned.
C      40..45              Reserved
C      46..50              User options (see note 5b)
C
C     Note 3:
C         If NLEQ1 terminates with IERR=2 (maximum iterations)
C         or  IERR=3 (small damping factor), you may try to continue
C         the iteration by increasing NITMAX or decreasing FCMIN
C         (see RWK) and setting QSUCC to 1.
C
C     Note 4 : Storage of user supplied banded Jacobian
C        In the band matrix case, the following lines may build
C        up the analytic Jacobian A;
C        Here AFL denotes the quadratic matrix A in dense form,
C        and ABD the rectangular matrix A in banded form :
C
C                   ML = IOPT(6)
C                   MU = IOPT(7)
C                   MH = MU+1
C                   DO 20 J = 1,N
C                     I1 = MAX0(1,J-MU)
C                     I2 = MIN0(N,J+ML)
C                     DO 10 I = I1,I2
C                       K = I-J+MH
C                       ABD(K,J) = AFL(I,J)
C           10        CONTINUE
C           20      CONTINUE
C
C         The total number of rows needed in  ABD  is  ML+MU+1 .
C         The  MU by MU  upper left triangle and the
C         ML by ML  lower right triangle are not referenced.
C
C     Note 5a:
C        The integrated time monitor calls the machine dependent
C        subroutine SECOND to get the current time stamp in form
C        of a real number (Single precision). As delivered, this
C        subroutine always return 0.0 as time stamp value. Refer
C        to the compiler- or library manual of the FORTRAN compiler
C        which you currently use to find out how to get the current
C        time stamp on your machine.
C
C     Note 5b:
C         The user options may be interpreted by the user replacable
C         routines N1SOUT, N1FACT, N1SOLV - the distributed version
C         of N1SOUT currently uses IOPT(46) as follows:
C         0 = standard plotdata output (may be postprocessed by a user-
C             written graphical program)
C         1 = plotdata output is suitable as input to the graphical
C             package GRAZIL (based on GKS), which has been developed
C             at ZIB. 
C
C
C*   Optional INTEGER input/output in IWK:
C    =======================================
C
C     Pos. Name          Meaning
C
C      1   NITER  IN/OUT Number of Newton-iterations
C      2                 reserved
C      3   NCORR  IN/OUT Number of corrector steps
C      4   NFCN   IN/OUT Number of FCN-evaluations
C      5   NJAC   IN/OUT Number of Jacobian generations or
C                        JAC-calls
C      6                 reserved
C      7                 reserved
C      8   NFCNJ  IN/OUT Number of FCN-evaluations for Jacobian
C                        approximation
C      9   NREJR1 IN/OUT Number of rejected Newton iteration steps
C                        done with a rank-1 approximated Jacobian
C     10..11             Reserved
C     12   IDCODE IN/OUT Output: The 8 decimal digits program identi-
C                        fication number ppppvvvv, consisting of the
C                        program code pppp and the version code vvvv.
C                        Input: If containing a negative number,
C                        it will only be overwritten by the identi-
C                        fication number, immediately followed by
C                        a return to the calling program.      
C     13..15             Reserved
C     16   NIWKFR OUT    First element of IWK which is free to be used
C                        as workspace between Newton iteration steps
C                        for standard linear solvers: 51
C     17   NRWKFR OUT    First element of RWK which is free to be used
C                        as workspace between Newton iteration steps.
C                        For standard linear solvers and numerically 
C                        approximated Jacobian computed by one of the 
C                        expressions:
C                        (N+7+NBROY)*N+61        for a full Jacobian
C                        (2*ML+MU+8+NBROY)*N+61  for a banded Jacobian
C                        If the Jacobian is computed by a user routine
C                        JAC, subtract N in both expressions.
C     18   LIWKA  OUT    Length of IWK currently required
C     19   LRWKA  OUT    Length of RWK currently required
C     20..22             Reserved
C     23   IFAIL  OUT    Set in case of failure of N1FACT (IERR=80),
C                        N1SOLV (IERR=81), FCN (IERR=82) or JAC(IERR=83)
C                        to the nonzero IFAIL value returned by the 
C                        routine indicating the failure .
C     24   ICONV  OUT    Current status of of the convergence monitor
C                        (only if convergence order monitor is on - 
C                         see IORMON(=IOPT(39)))
C                        =0: No convergence indicated yet
C                        =1: Damping factor is 1.0d0
C                        =2: Superlinear convergence in progress
C                        =3: Quadratic convergence in progress
C     25..30             Reserved
C     31   NITMAX IN     Maximum number of permitted iteration
C                        steps (default: 50)
C     32                 Reserved
C     33   NEW    IN/OUT Count of consecutive rank-1 updates
C     34..35             Reserved
C     36   NBROY  IN     Maximum number of possible consecutive 
C                        iterative Broyden steps. The total real 
C                        workspace needed (RWK) depends on this value
C                        (see LRWK above).
C                        Default is N (see parameter N),
C                        if MSTOR=0 (=IOPT(4)), 
C                        and ML+MU+1 (=IOPT(6)+IOPT(7)+1), if MSTOR=1
C                        (but minimum is always 10) - 
C                        provided that Broyden is allowed. 
C                        If Broyden is inhibited, NBROY is always set to
C                        zero.
C     37..50             Reserved
C
C*   Optional REAL input/output in RWK:
C    ====================================
C
C     Pos. Name          Meaning
C
C      1..16             Reserved
C     17   CONV   OUT    The achieved relative accuracy after the  
C                        current step
C     18   SUMX   OUT    Natural level (not Normx of printouts)
C                        of the current iterate, i.e. Sum(DX(i)**2),
C                        where DX = scaled Newton correction.
C     19   DLEVF  OUT    Standard level (not Normf of printouts)
C                        of the current iterate, i.e. Norm2(F(X)),
C                        where F =  nonlinear problem function.
C     20   FCBND  IN     Bounded damping strategy restriction factor
C                        (Default is 10)
C     21   FCSTRT IN     Damping factor for first Newton iteration -
C                        overrides option NONLIN, if set (see note 6)
C     22   FCMIN  IN     Minimal allowed damping factor (see note 6)
C     23   SIGMA  IN     Broyden-approximation decision parameter
C                        Required choice: SIGMA.GE.1. Increasing this
C                        parameter make it less probable that the algo-
C                        rith performs rank-1 updates.
C                        Rank-1 updates are inhibited, if 
C                        SIGMA.GT.1/FCMIN is set. (see note 6)
C     24   SIGMA2 IN     Decision parameter about increasing damping
C                        factor to corrector if predictor is small.
C                        Required choice: SIGMA2.GE.1. Increasing this
C                        parameter make it less probable that the algo-
C                        rith performs rank-1 updates.
C     25                 Reserved
C     26   AJDEL  IN     Jacobian approximation without feedback:
C                        Relative pertubation for components
C                        (Default: sqrt(epmach*10), epmach: relative
C                         machine precision) 
C     27   AJMIN  IN     Jacobian approximation without feedback:
C                        Threshold value (Default: 0.0d0)
C                          The absolute pertubation for component k is
C                          computed by 
C                          DELX := AJDEL*max(abs(Xk),AJMIN)
C     28  ETADIF  IN     Jacobian approximation with feedback:
C                        Target value for relative pertubation ETA of X
C                        (Default: 1.0d-6)
C     29  ETAINI  IN     Jacobian approximation with feedback:
C                        Initial value for denominator differences
C                        (Default: 1.0d-6)
C     30..50             Reserved
C
C     Note 6:
C       The default values of the internal parameters may be obtained
C       from the monitor output with at least IOPT field MPRMON set to 2
C       and by initializing the corresponding RWK-fields to zero. 
C
C*   Error and warning messages:
C    ===========================
C
C      1    Termination, since jacobian matrix became singular
C      2    Termination after NITMAX iterations ( as indicated by
C           input parameter NITMAX=IWK(31) )
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
C     10    Integer or real workspace too small
C     20    Bad input to dimensional parameter N
C     21    Nonpositive value for RTOL supplied
C     22    Negative scaling value via vector XSCAL supplied
C     30    One or more fields specified in IOPT are invalid
C           (for more information, see error-printout)
C     80    Error signalled by linear solver routine N1FACT,
C           for more detailed information see IFAIL-value
C           stored to IWK(23)
C     81    Error signalled by linear solver routine N1SOLV,
C           for more detailed information see IFAIL-value
C           stored to IWK(23)
C           (not used by standard routine N1SOLV)
C     82    Error signalled by user routine FCN (Nonzero value
C           returned via IFAIL-flag; stored to IWK(23) )
C     83    Error signalled by user routine JAC (Nonzero value
C           returned via IFAIL-flag; stored to IWK(23) )
C
C     Note 7 : in case of failure:
C        -    use non-standard options
C        -    or turn to Newton-algorithm with rank strategy
C        -    use another initial guess
C        -    or reformulate model
C        -    or apply continuation techniques
C
C*    Machine dependent constants used:
C     =================================
C
C     DOUBLE PRECISION EPMACH  in  N1PCHK, N1INT
C     DOUBLE PRECISION GREAT   in  N1PCHK
C     DOUBLE PRECISION SMALL   in  N1PCHK, N1INT, N1SCAL
C
C*    Subroutines called: N1PCHK, N1INT
C
C     ------------------------------------------------------------
C*    End Prologue
C
C*    Summary of changes:
C     ===================
C      
C     2.2.1  91, May 30     Time monitor included
C     2.2.2  91, May 30     Bounded damping strategy implemented
C     2.2.3  91, June 19    AJDEL, AJMIN as RWK-options for JACGEN.EQ.2,
C                           ETADIF, ETAINI as RWK-opts. for JACGEN.EQ.3
C                           FCN-count changed for anal. Jacobian
C     2.2.4  91, August  9  Convergence order monitor included
C     2.2.5  91, August 13  Standard Broyden updates replaced by
C                           iterative Broyden
C     2.2.6  91, Sept.  16  Damping factor reduction by FCN-fail imple-
C                           mented
C     2.3    91, Dec.   20  New Release for CodeLib
C            92, March  11  Level of accepted simplified correction
C                           stored to RWK(IRWKI+4)
C            00, July   12  RTOL output-value bug fixed
C            06, Jan.   24  IERR=5 no longer returned if residuum of
C                           final iterate is exactly zero
C     2.4    09, May    29  Changed routines N1FACT and N1SOLV to use
C                           LAPACK Routines DGETRF, DGETRS, DGBTRF and
C                           DGBTRS instead of LINPACK routines to solve 
C                           linear systems.
C            10, July   26  Subroutine N1INT: Initialization of unitialized
C                           Variable FCMON fixed.
C   
C     ------------------------------------------------------------
C
C     PARAMETER (IRWKI=xx, LRWKI=yy)  
C     IRWKI: Start position of internally used RWK part
C     LRWKI: Length of internally used RWK part
C     (current values see parameter statement below)
C
C     INTEGER L4,L5,L51,L6,L61,L62,L63,L7,L71,L8,L9,L10,L11,L12,L121,
C             L13,L14,L20
C     Starting positions in RWK of formal array parameters of internal
C     routine N1INT (dynamically determined in driver routine NLEQ1,
C     dependent on N and options setting)
C
C     Further RWK positions (only internally used)
C
C     Position  Name     Meaning
C
C     IRWKI     FCKEEP   Damping factor of previous successfull iter.
C     IRWKI+1   FCA      Previous damping factor
C     IRWKI+2   FCPRI    A priori estimate of damping factor
C     IRWKI+3   DMYCOR   Number My of latest corrector damping factor
C                        (kept for use in rank-1 decision criterium)
C     IRWKI+4   SUMXS    natural level of accepted simplified correction
C     IRWKI+(5..LRWKI-1) Free
C
C     Internal arrays stored in RWK (see routine N1INT for descriptions)
C
C     Position  Array         Type   Remarks
C
C     L4        A(M1,N)       Perm   M1=N (full mode) or 
C                                    M1=2*IOPT(6)+IOPT(7)+1 (band mode)
C     L41       DXSAVE(N,NBROY)
C                             Perm   NBROY=IWK(36) (Default: N or 0)
C     L5        DX(N)         Perm  
C     L51       DXQ(N)        Perm 
C     L6        XA(N)         Perm
C     L61       F(N)          Perm
C     L62       FW(N)         Perm
C     L63       XWA(N)        Perm
C     L7        FA(N)         Perm
C     L71       ETA(N)        Perm   Only used for JACGEN=IOPT(3)=3
C     L9        XW(N)         Temp
C     L11       DXQA(N)       Temp
C     L12       T1(N)         Temp
C     L121      T2(N)         Temp
C     L13       T3(N)         Temp
C     L14                     Temp   Start position of array workspace 
C                                    needed for linear solver  
C     
C
      EXTERNAL N1INT
      INTRINSIC DBLE
      INTEGER IRWKI, LRWKI
      PARAMETER (IRWKI=51, LRWKI=10)  
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER NITMAX,LUERR,LUMON,LUSOL,MPRERR,MPRMON,MPRSOL,
     $MSTOR,M1,M2,NRWKFR,NRFRIN,NRW,NIWKFR,NIFRIN,NIW,NONLIN,JACGEN
      INTEGER L4,L41,L5,L51,L6,L61,L62,L63,L7,L71,L8,L9,L11,L12,L121,
     $L13,L14,L20
      DOUBLE PRECISION FC,FCMIN,PERCI,PERCR
      LOGICAL QINIMO,QRANK1,QFCSTR,QSUCC,QBDAMP,QSIMPL
      CHARACTER CHGDAT*20, PRODCT*8
C     Which version ?
      LOGICAL QVCHK
      INTEGER IVER
      PARAMETER( IVER=21112401 )
C
C     Version: 2.4.0.1           Latest change:
C     -----------------------------------------
C
      DATA      CHGDAT      /'July 26, 2010       '/
      DATA      PRODCT      /'NLEQ1   '/
C*    Begin
      IERR = 0
      QVCHK = IWK(12).LT.0
      IWK(12) = IVER
      IF (QVCHK) RETURN
C        Print error messages?
      MPRERR = IOPT(11)
      LUERR = IOPT(12)
      IF (LUERR .EQ. 0) THEN
        LUERR = 6
        IOPT(12)=LUERR
      ENDIF
C        Print iteration monitor?
      MPRMON = IOPT(13)
      LUMON = IOPT(14)
      IF (LUMON .LE. 0 .OR. LUMON .GT. 99) THEN
        LUMON = 6
        IOPT(14)=LUMON
      ENDIF
C        Print intermediate solutions?
      MPRSOL = IOPT(15)
      LUSOL = IOPT(16)
      IF (LUSOL .EQ. 0) THEN
        LUSOL = 6
        IOPT(16)=LUSOL
      ENDIF
C        Print time summary statistics?
      MPRTIM = IOPT(19)
      LUTIM = IOPT(20)
      IF (LUTIM .EQ. 0) THEN
        LUTIM = 6
        IOPT(20)=LUTIM
      ENDIF
      QSUCC = IOPT(1).EQ.1
      QINIMO = MPRMON.GE.1.AND..NOT.QSUCC
C     Print NLEQ1 heading lines
      IF(QINIMO)THEN
10000   FORMAT('   N L E Q 1  *****  V e r s i o n  ',
     $         '2 . 4 ***',//,1X,'Newton-Method ',
     $         'for the solution of nonlinear systems',//)
        WRITE(LUMON,10000)
      ENDIF
C     Check input parameters and options
      CALL N1PCHK(N,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
C     Exit, if any parameter error was detected till here
      IF (IERR.NE.0) RETURN 
C
      MSTOR=IOPT(4)
      IF (MSTOR.EQ.0) THEN
        M1=N
        M2=N
      ELSE IF (MSTOR.EQ.1) THEN
        ML=IOPT(6)
        MU=IOPT(7)
        M1=2*ML+MU+1
        M2=ML+MU+1
      ENDIF
      JACGEN=IOPT(3)
      IF (JACGEN.EQ.0) JACGEN=2
      IOPT(3)=JACGEN
      QRANK1=IOPT(32).EQ.1
      QSIMPL=IOPT(34).EQ.1
      IF (QRANK1) THEN
        NBROY=IWK(36)
        IF (NBROY.EQ.0) NBROY=MAX(M2,10)
        IWK(36)=NBROY
      ELSE
        NBROY=0
      ENDIF
C     WorkSpace: RWK
      L4=IRWKI+LRWKI
      L41=L4+M1*N
      L5=L41+NBROY*N
      L51=L5+N
      L6=L51+N
      L61=L6+N
      L62=L61+N
      L63=L62+N
      L7=L63+N
      L71=L7+N
      IF (JACGEN.NE.3) THEN
        L8=L71
      ELSE
        L8=L71+N
      ENDIF
      NRWKFR = L8
      L9=L8
      L11=L9+N
      L12=L11+N
      L121=L12+N
      L13=L121+N
      L14=L13+N
      NRW=L14-1
C     End WorkSpace at NRW
C     WorkSpace: IWK
      L20=51
      NIWKFR = L20
      IF (QRANK1.OR.QSIMPL) NIWKFR = NIWKFR+N
      NIW=L20-1
C     End WorkSpace at NIW
      IWK(16) = NIW+1
      IWK(17) = NRW+1
      NIFRIN = NIW+1
      NRFRIN = NRW+1
C
      IF(NRW.GT.LRWK.OR.NIW.GT.LIWK)THEN
        IERR=10
      ELSE
        IF(QINIMO)THEN
          PERCR = DBLE(NRW)/DBLE(LRWK)*100.0D0
          PERCI = DBLE(NIW)/DBLE(LIWK)*100.0D0
C         Print statistics concerning workspace usage
10050     FORMAT(' Real    Workspace declared as ',I9,
     $    ' is used up to ',I9,' (',F5.1,' percent)',//,
     $    ' Integer Workspace declared as ',I9,
     $    ' is used up to ',I9,' (',F5.1,' percent)',//)
          WRITE(LUMON,10050)LRWK,NRW,PERCR,LIWK,NIW,PERCI
        ENDIF
        IF(QINIMO)THEN
10051     FORMAT(/,' N =',I4,//,' Prescribed relative ',
     $    'precision',D10.2,/)
          WRITE(LUMON,10051)N,RTOL
10052     FORMAT(' The Jacobian is supplied by ',A)
          IF (JACGEN.EQ.1) THEN
            WRITE(LUMON,10052) 'a user subroutine'
          ELSE IF (JACGEN.EQ.2) THEN
             WRITE(LUMON,10052) 
     $        'numerical differentiation (without feedback strategy)'
          ELSE IF (JACGEN.EQ.3) THEN
             WRITE(LUMON,10052) 
     $        'numerical differentiation (feedback strategy included)'
          ENDIF
10055     FORMAT(' The Jacobian will be stored in ',A,' mode')
          IF (MSTOR.EQ.0) THEN
            WRITE(LUMON,10055) 'full'
          ELSE IF (MSTOR.EQ.1) THEN
            WRITE(LUMON,10055) 'banded'
10056       FORMAT(' Lower bandwidth : ',I3,'   Upper bandwidth : ',I3)
            WRITE(LUMON,10056) ML,MU
          ENDIF
10057     FORMAT(' Automatic row scaling of the Jacobian is ',A,/)
          IF (IOPT(35).EQ.1) THEN
            WRITE(LUMON,10057) 'inhibited'
          ELSE
            WRITE(LUMON,10057) 'allowed'
          ENDIF
        ENDIF
        NONLIN=IOPT(31)
        IF (IOPT(38).EQ.0) QBDAMP = NONLIN.EQ.4
        IF (IOPT(38).EQ.1) QBDAMP = .TRUE.
        IF (IOPT(38).EQ.2) QBDAMP = .FALSE.
        IF (QBDAMP) THEN
          IF (RWK(20).LT.ONE) RWK(20) = TEN
        ENDIF
10064   FORMAT(' Rank-1 updates are ',A)
        IF (QINIMO) THEN
          IF (QRANK1) THEN
            WRITE(LUMON,10064) 'allowed'
          ELSE
            WRITE(LUMON,10064) 'inhibited'
          ENDIF
10065     FORMAT(' Problem is specified as being ',A)
          IF (NONLIN.EQ.1) THEN
            WRITE(LUMON,10065) 'linear'
          ELSE IF (NONLIN.EQ.2) THEN
            WRITE(LUMON,10065) 'mildly nonlinear'
          ELSE IF (NONLIN.EQ.3) THEN
            WRITE(LUMON,10065) 'highly nonlinear'
          ELSE IF (NONLIN.EQ.4) THEN
            WRITE(LUMON,10065) 'extremely nonlinear'
          ENDIF
10066     FORMAT(' Bounded damping strategy is ',A,:,/, 
     $           ' Bounding factor is ',D10.3)
          IF (QBDAMP) THEN
            WRITE(LUMON,10066) 'active', RWK(20)
          ELSE
            WRITE(LUMON,10066) 'off'
          ENDIF
10067     FORMAT(' Special mode: ',A,' Newton iteration will be done')
          IF (IOPT(33).EQ.1) WRITE(LUMON,10067) 'Ordinary'
          IF (IOPT(34).EQ.1) WRITE(LUMON,10067) 'Simplified'
        ENDIF
C       Maximum permitted number of iteration steps
        NITMAX=IWK(31)
        IF (NITMAX.LE.0) NITMAX=50
        IWK(31)=NITMAX
10068   FORMAT(' Maximum permitted number of iteration steps : ',
     $         I6)
        IF (QINIMO) WRITE(LUMON,10068) NITMAX
C       Initial damping factor for highly nonlinear problems
        QFCSTR=RWK(21).GT.ZERO
        IF (.NOT.QFCSTR) THEN
          RWK(21)=1.0D-2
          IF (NONLIN.EQ.4) RWK(21)=1.0D-4
        ENDIF
C       Minimal permitted damping factor
        IF (RWK(22).LE.ZERO) THEN
          RWK(22)=1.0D-4
          IF (NONLIN.EQ.4) RWK(22)=1.0D-8
        ENDIF
        FCMIN=RWK(22)
C       Rank1 decision parameter SIGMA
        IF (RWK(23).LT.ONE) RWK(23)=3.0D0
        IF (.NOT.QRANK1) RWK(23)=10.0D0/FCMIN
C       Decision parameter about increasing too small predictor
C       to greater corrector value
        IF (RWK(24).LT.ONE) RWK(24)=10.0D0/FCMIN       
C       Starting value of damping factor (FCMIN.LE.FC.LE.1.0)
        IF(NONLIN.LE.2.AND..NOT.QFCSTR)THEN
C         for linear or mildly nonlinear problems
          FC = ONE
        ELSE
C         for highly or extremely nonlinear problems
          FC = RWK(21)
        ENDIF
C       Simplied Newton iteration implies ordinary Newton it. mode
        IF (IOPT(34).EQ.1) IOPT(33)=1
C       If ordinary Newton iteration, factor is always one
        IF (IOPT(33).EQ.1) FC=1.0D0
        RWK(21)=FC
        IF (MPRMON.GE.2.AND..NOT.QSUCC) THEN
10069     FORMAT(//,' Internal parameters:',//,
     $      ' Starting value for damping factor FCSTART = ',D9.2,/,
     $      ' Minimum allowed damping factor FCMIN = ',D9.2,/,
     $      ' Rank-1 updates decision parameter SIGMA = ',D9.2)
          WRITE(LUMON,10069) RWK(21),FCMIN,RWK(23)
        ENDIF
C       Store lengths of currently required workspaces
        IWK(18) = NIFRIN-1
        IWK(19) = NRFRIN-1
C
C       Initialize and start time measurements monitor
C
        IF ( IOPT(1).EQ.0 .AND. MPRTIM.NE.0 ) THEN
          CALL MONINI (' NLEQ1',LUTIM)
          CALL MONDEF (0,'NLEQ1')
          CALL MONDEF (1,'FCN')
          CALL MONDEF (2,'Jacobi')
          CALL MONDEF (3,'Lin-Fact')
          CALL MONDEF (4,'Lin-Sol')
          CALL MONDEF (5,'Output')
          CALL MONSTR (IERR)
        ENDIF
C
C
C
        IERR=-1
C       If IERR is unmodified on exit, successive steps are required
C       to complete the Newton iteration
        IF (NBROY.EQ.0) NBROY=1
        CALL N1INT(N,FCN,JAC,X,XSCAL,RTOL,NITMAX,NONLIN,IOPT,IERR,
     $  LRWK,RWK,NRFRIN,LIWK,IWK,NIFRIN,
     $  M1,M2,NBROY,
     $  RWK(L4),RWK(L41),RWK(L5),RWK(L51),RWK(L6),RWK(L63),RWK(L61),
     $  RWK(L7),
     $  RWK(L71),RWK(L9),RWK(L62),RWK(L11),RWK(L12),RWK(L121),RWK(L13),
     $  RWK(21),RWK(22),RWK(23),RWK(24),RWK(IRWKI+1),RWK(IRWKI),
     $  RWK(IRWKI+2),RWK(IRWKI+3),RWK(17),RWK(18),RWK(IRWKI+4),RWK(19),
     $  MSTOR,MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,IWK(1),IWK(3),
     $  IWK(4),IWK(5),IWK(8),IWK(9),IWK(33),IWK(24),QBDAMP)
C
        IF (MPRTIM.NE.0.AND.IERR.NE.-1.AND.IERR.NE.10) THEN
          CALL MONHLT
          CALL MONPRT
        ENDIF
C
C       Free workspaces, so far not used between steps
        IWK(16) = NIWKFR
        IWK(17) = NRWKFR
      ENDIF
C     Print statistics
      IF (MPRMON.GE.1.AND.IERR.NE.-1.AND.IERR.NE.10) THEN
10080   FORMAT(/, '   ******  Statistics * ', A8, ' *******', /,
     $            '   ***  Newton iterations : ', I7,'  ***', /,
     $            '   ***  Corrector steps   : ', I7,'  ***', /,
     $            '   ***  Rejected rk-1 st. : ', I7,'  ***', /,
     $            '   ***  Jacobian eval.    : ', I7,'  ***', /,
     $            '   ***  Function eval.    : ', I7,'  ***', /,
     $            '   ***  ...  for Jacobian : ', I7,'  ***', /,
     $            '   *************************************', /)
        WRITE (LUMON,10080) PRODCT,IWK(1),IWK(3),IWK(9),IWK(5),
     $  IWK(4),IWK(8)
      ENDIF
C     Print workspace requirements, if insufficient
      IF (IERR.EQ.10) THEN
10090   FORMAT(//,' ',20('*'),'Workspace Error',20('*'))
        IF (MPRERR.GE.1) WRITE(LUERR,10090)
        IF(NRW.GT.LRWK)THEN
10091     FORMAT(/,' Real Workspace dimensioned as',1X,I9,
     $    1X,'must be enlarged at least up to ',
     $    I9,//)
          IF (MPRERR.GE.1) WRITE(LUERR,10091)LRWK,NRFRIN-1
        ENDIF
        IF(NIW.GT.LIWK)THEN
10092     FORMAT(/,' Integer Workspace dimensioned as ',
     $    I9,' must be enlarged at least up ',
     $    'to ',I9,//)
          IF (MPRERR.GE.1) WRITE(LUERR,10092)LIWK,NIFRIN-1
        ENDIF
      ENDIF
C     End of subroutine NLEQ1
      RETURN
      END
C
      SUBROUTINE N1PCHK(N,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
C*    Begin Prologue N1PCHK
      INTEGER N
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 1 P C H K : Checking of input parameters and options
C                   for NLEQ1.
C
C*    Parameters:
C     ===========
C
C     See parameter description in driver routine.
C
C*    Subroutines called: ZIBCONST
C
C*    Machine dependent constants used:
C     =================================
C
C     EPMACH = relative machine precision
C     GREAT = squareroot of maxreal divided by 10
C     SMALL = squareroot of "smallest positive machine number
C             divided by relative machine precision"
      DOUBLE PRECISION EPMACH,GREAT,SMALL
C
C     ------------------------------------------------------------
C*    End Prologue
C
      INTRINSIC DBLE
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
      PARAMETER (NUMOPT=50)
      INTEGER MSTOR,IOPTL(NUMOPT),IOPTU(NUMOPT)
      DOUBLE PRECISION TOLMIN,TOLMAX,DEFSCL
C
      DATA IOPTL /0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1,
     $            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     $            0,0,0,0,0,0,0,0,0,0,
     $            -9999999,-9999999,-9999999,-9999999,-9999999/
      DATA IOPTU /1,1,3,1,0,9999999,9999999,0,1,0,3,99,6,99,3,99,0,0,1,
     $            99,0,0,0,0,0,0,0,0,0,0,4,1,1,1,1,
     $            0,0,2,3,0,0,0,0,0,0,
     $            9999999,9999999,9999999,9999999,9999999/
C
      CALL ZIBCONST(EPMACH,SMALL)
      GREAT  = 1.0D0/SMALL
      IERR = 0
C        Print error messages?
      MPRERR = IOPT(11)
      LUERR = IOPT(12)
      IF (LUERR .LE. 0 .OR. LUERR .GT. 99) THEN
        LUERR = 6
        IOPT(12)=LUERR
      ENDIF
C
C     Checking dimensional parameter N
      IF ( N.LE.0 ) THEN
        IF (MPRERR.GE.1)  WRITE(LUERR,10011) N
10011   FORMAT(/,' Error: Bad input to dimensional parameter N supplied'
     $         ,/,8X,'choose N positive, your input is: N = ',I5)
        IERR = 20
      ENDIF
C
C     Problem type specification by user
      NONLIN=IOPT(31)
      IF (NONLIN.EQ.0) NONLIN=3
      IOPT(31)=NONLIN
C
C     Checking and conditional adaption of the user-prescribed RTOL
      IF (RTOL.LE.ZERO) THEN
        IF (MPRERR.GE.1) 
     $      WRITE(LUERR,'(/,A)') ' Error: Nonpositive RTOL supplied'
        IERR = 21
      ELSE
        TOLMIN = EPMACH*TEN*DBLE(N)
        IF(RTOL.LT.TOLMIN) THEN
          RTOL = TOLMIN
          IF (MPRERR.GE.2) 
     $      WRITE(LUERR,10012) 'increased ','smallest',RTOL
        ENDIF
        TOLMAX = 1.0D-1
        IF(RTOL.GT.TOLMAX) THEN
          RTOL = TOLMAX
          IF (MPRERR.GE.2) 
     $      WRITE(LUERR,10012) 'decreased ','largest',RTOL
        ENDIF
10012   FORMAT(/,' Warning: User prescribed RTOL ',A,'to ',
     $         'reasonable ',A,' value RTOL = ',D11.2)
      ENDIF
C     
C     Test user prescribed accuracy and scaling on proper values
      IF (N.LE.0) RETURN 
      IF (NONLIN.GE.3) THEN
        DEFSCL = RTOL
      ELSE
        DEFSCL = ONE
      ENDIF
      DO 10 I=1,N
        IF (XSCAL(I).LT.ZERO) THEN
          IF (MPRERR.GE.1) THEN 
            WRITE(LUERR,10013) I
10013       FORMAT(/,' Error: Negative value in XSCAL(',I5,') supplied')
          ENDIF
          IERR = 22
        ENDIF
        IF (XSCAL(I).EQ.ZERO) XSCAL(I) = DEFSCL
        IF ( XSCAL(I).GT.ZERO .AND. XSCAL(I).LT.SMALL ) THEN
          IF (MPRERR.GE.2) THEN
            WRITE(LUERR,10014) I,XSCAL(I),SMALL
10014       FORMAT(/,' Warning: XSCAL(',I5,') = ',D9.2,' too small, ',
     $             'increased to',D9.2)
          ENDIF
          XSCAL(I) = SMALL
        ENDIF
        IF (XSCAL(I).GT.GREAT) THEN
          IF (MPRERR.GE.2) THEN
            WRITE(LUERR,10015) I,XSCAL(I),GREAT
10015       FORMAT(/,' Warning: XSCAL(',I5,') = ',D9.2,' too big, ',
     $             'decreased to',D9.2)
          ENDIF
          XSCAL(I) = GREAT
        ENDIF
10    CONTINUE
C     Special dependence on Jacobian's storage mode for ML and MU
      MSTOR = IOPT(4)
      IF (MSTOR.EQ.0) THEN
        IOPTU(6)=0
        IOPTU(7)=0
      ELSE IF (MSTOR.EQ.1) THEN
        IOPTU(6)=N-1
        IOPTU(7)=N-1
      ENDIF
C     Checks options
      DO 20 I=1,30
        IF (IOPT(I).LT.IOPTL(I) .OR. IOPT(I).GT.IOPTU(I)) THEN
          IERR=30
          IF (MPRERR.GE.1) THEN
            WRITE(LUERR,20001) I,IOPT(I),IOPTL(I),IOPTU(I)
20001       FORMAT(' Invalid option specified: IOPT(',I2,')=',I12,';',
     $             /,3X,'range of permitted values is ',I8,' to ',I8)
          ENDIF
        ENDIF
20    CONTINUE
C     End of subroutine N1PCHK
      RETURN
      END
C
      SUBROUTINE N1INT(N,FCN,JAC,X,XSCAL,RTOL,NITMAX,NONLIN,IOPT,IERR,
     $LRWK,RWK,NRWKFR,LIWK,IWK,NIWKFR,M1,M2,NBROY,
     $A,DXSAVE,DX,DXQ,XA,XWA,F,FA,ETA,XW,FW,DXQA,T1,T2,T3,FC,FCMIN,
     $SIGMA,SIGMA2,FCA,FCKEEP,FCPRI,DMYCOR,CONV,SUMX,SUMXS,DLEVF,MSTOR,
     $MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,NITER,NCORR,NFCN,NJAC,
     $NFCNJ,NREJR1,NEW,ICONV,QBDAMP)
C*    Begin Prologue N1INT
      INTEGER N
      EXTERNAL FCN,JAC
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER NITMAX,NONLIN
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER NRWKFR,LIWK
      INTEGER IWK(LIWK)
      INTEGER NIWKFR,M1,M2,NBROY
      DOUBLE PRECISION A(M1,N),DXSAVE(N,NBROY)
      DOUBLE PRECISION DX(N),DXQ(N),XA(N),XWA(N),F(N),FA(N),ETA(N)
      DOUBLE PRECISION XW(N),FW(N),DXQA(N),T1(N),T2(N),T3(N)
      DOUBLE PRECISION FC,FCMIN,SIGMA,SIGMA2,FCA,FCKEEP,CONV,SUMX,SUMXS,
     $                 DLEVF,FCPRI,DMYCOR
      INTEGER MSTOR,MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,NITER,
     $NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW,ICONV
      LOGICAL QBDAMP
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 1 I N T : Core routine for NLEQ1 .
C     Damped Newton-algorithm for systems of highly nonlinear
C     equations especially designed for numerically sensitive
C     problems.
C
C*    Parameters:
C     ===========
C
C       N,FCN,JAC,X,XSCAL,RTOL   
C                         See parameter description in driver routine
C
C       NITMAX      Int    Maximum number of allowed iterations
C       NONLIN      Int    Problem type specification
C                          (see IOPT-field NONLIN)
C       IOPT        Int    See parameter description in driver routine
C       IERR        Int    See parameter description in driver routine
C       LRWK        Int    Length of real workspace
C       RWK(LRWK)   Dble   Real workspace array
C       NRWKFR      Int    First free position of RWK on exit 
C       LIWK        Int    Length of integer workspace
C       IWK(LIWK)   Int    Integer workspace array
C       NIWKFR      Int    First free position of IWK on exit 
C       M1          Int    Leading dimension of Jacobian array A
C                          for full case Jacobian: N
C                          for banded Jacobian: 2*ML+MU+1;
C                          ML, MU see IOPT-description in driver routine
C       M2          Int    for full case Jacobian: N
C                          for banded Jacobian: ML+MU+1
C       NBROY       Int    Maximum number of possible consecutive
C                          iterative Broyden steps. (See IWK(36))
C       A(M1,N)     Dble   Holds the Jacobian matrix (decomposed form
C                          after call of linear decomposition
C                          routine)
C       DXSAVE(X,NBROY)
C                   Dble   Used to save the quasi Newton corrections of
C                          all previously done consecutive Broyden
C                          steps.
C       DX(N)       Dble   Current Newton correction
C       DXQ(N)      Dble   Simplified Newton correction J(k-1)*X(k)
C       XA(N)       Dble   Previous Newton iterate
C       XWA(N)      Dble   Scaling factors used for latest decomposed
C                          Jacobian for column scaling - may differ
C                          from XW, if Broyden updates are performed
C       F(N)        Dble   Function (FCN) value of current iterate
C       FA(N)       Dble   Function (FCN) value of previous iterate
C       ETA(N)      Dble   Jacobian approximation: updated scaled
C                          denominators
C       XW(N)       Dble   Scaling factors for iteration vector
C       FW(N)       Dble   Scaling factors for rows of the system
C       DXQA(N)     Dble   Previous Newton correction
C       T1(N)       Dble   Workspace for linear solvers and internal
C                          subroutines
C       T2(N)       Dble   Workspace array for internal subroutines
C       T3(N)       Dble   Workspace array for internal subroutines
C       FC          Dble   Current Newton iteration damping factor.
C       FCMIN       Dble   Minimum permitted damping factor. If
C                          FC becomes smaller than this value, one
C                          of the following may occur:
C                          a.    Recomputation of the Jacobian
C                                matrix by means of difference
C                                approximation (instead of Rank1
C                                update), if Rank1 - update
C                                previously was used
C                          b.    Fail exit otherwise
C       SIGMA       Dble   Decision parameter for rank1-updates
C       SIGMA2      Dble   Decision parameter for damping factor
C                          increasing to corrector value
C       FCA         Dble   Previous Newton iteration damping factor.
C       FCKEEP      Dble   Keeps the damping factor as it is at start
C                          of iteration step.
C       CONV        Dble   Scaled maximum norm of the Newton-
C                          correction. Passed to RWK-field on output.
C       SUMX        Dble   Square of the natural level (see equal-
C                          named RWK-output field)
C       SUMXS       Dble   Square of the "simplified" natural level
C                          (see equal-named RWK-internal field)
C       DLEVF       Dble   The standard level (see equal-
C                          named RWK-output field)
C       MSTOR       Int    see description of IOPT-field MSTOR
C       MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,
C       NITER,NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW :
C                          See description of equal named IWK-fields
C                          in the driver subroutine
C       QBDAMP      Logic  Flag, that indicates, whether bounded damping
C                          strategy is active:
C                          .true.  = bounded damping strategy is active
C                          .false. = normal damping strategy is active
C
C*    Internal double variables
C     =========================
C
C       AJDEL    See RWK(26) (num. diff. without feedback)
C       AJMIN    See RWK(27) (num. diff. without feedback)
C       CONVA    Holds the previous value of CONV .
C       DMUE     Temporary value used during computation of damping 
C                factors predictor.
C       EPDIFF   sqrt(10*epmach) (num. diff. with feedback)
C       ETADIF   See description of RWK(28) (num. diff. with feedback)
C       ETAINI   Initial value for all ETA-components (num. diff. fb.)
C       ETAMAX   Maximum allowed pertubation (num. diff. with feedback)
C       ETAMIN   Minimum allowed pertubation (num. diff. with feedback)
C       FCDNM    Used to compute the denominator of the damping 
C                factor FC during computation of it's predictor,
C                corrector and aposteriori estimate (in the case of
C                performing a Rank1 update) .
C       FCK2     Aposteriori estimate of FC.
C       FCH      Temporary storage of new corrector FC.
C       FCMIN2   FCMIN**2 . Used for FC-predictor computation.
C       FCNUMP   Gets the numerator of the predictor formula for FC.
C       FCNMP2   Temporary used for predictor numerator computation.
C       FCNUMK   Gets the numerator of the corrector computation 
C                of FC .
C       SUMXA    Natural level of the previous iterate.
C       TH       Temporary variable used during corrector- and 
C                aposteriori computations of FC.
C
C*    Internal integer variables
C     ==========================
C
C     IFAIL      Gets the return value from subroutines called from
C                N1INT (N1FACT, N1SOLV, FCN, JAC) 
C     ISCAL      Holds the scaling option from the IOPT-field ISCAL      
C     MODE       Matrix storage mode (see IOPT-field MODE) 
C     NRED       Count of successive corrector steps
C     NILUSE     Gets the amount of IWK used by the linear solver
C     NRLUSE     Gets the amount of RWK used by the linear solver
C     NIWLA      Index of first element of IWK provided to the
C                linear solver
C     NRWLA      Index of first element of RWK provided to the
C                linear solver
C     LIWL       Holds the maximum amount of integer workspace
C                available to the linear solver
C     LRWL       Holds the maximum amount of real workspace
C                available to the linear solver
C
C
C*    Internal logical variables
C     ==========================
C
C     QGENJ      Jacobian updating technique flag:
C                =.TRUE.  : Call of analytical subroutine JAC or
C                           numerical differentiation
C                =.FALSE. : rank1- (Broyden-) update
C     QINISC     Iterate initial-scaling flag:
C                =.TRUE.  : at first call of N1SCAL
C                =.FALSE. : at successive calls of N1SCAL
C     QSUCC      See description of IOPT-field QSUCC.
C     QJCRFR     Jacobian refresh flag:
C                set to .TRUE. if damping factor gets too small
C                and Jacobian was computed by rank1-update. 
C                Indicates, that the Jacobian needs to be recomputed
C                by subroutine JAC or numerical differentation.
C     QLINIT     Initialization state of linear solver workspace:
C                =.FALSE. : Not yet initialized
C                =.TRUE.  : Initialized - N1FACT has been called at
C                           least one time.
C     QSCALE     Holds the value of .NOT.QNSCAL. See description
C                of IOPT-field QNSCAL.
C
C*    Subroutines called:
C     ===================
C
C       N1FACT, N1SOLV, N1JAC,  N1JACB, N1JCF,  N1JCFB, N1LVLS, 
C       N1SCRF, N1SCRB, N1SOUT, N1PRV1, N1PRV2, N1SCAL,
C       MONON,  MONOFF
C
C*    Functions called:
C     =================
C
C       ZIBCONST, WNORM
C
C*    Machine constants used
C     ======================
C
      DOUBLE PRECISION EPMACH,SMALL
C 
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL N1FACT, N1SOLV, N1JAC,  N1JACB, N1JCF,  N1JCFB, N1LVLS,
     $         N1SCRF, N1SCRB, N1SOUT, N1PRV1, N1PRV2, N1SCAL,
     $         MONON,  MONOFF, WNORM
      INTRINSIC DSQRT,DMIN1,MAX0,MIN0
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=10.0D0)
      INTEGER IFAIL,ILOOP,ISCAL,K,MODE,NRED,NILUSE,NRLUSE,NIWLA,
     $NRWLA,LIWL,LRWL,L1,L2,JACGEN
      DOUBLE PRECISION AJDEL,AJMIN,ALFA1,ALFA2,ALFA,BETA,CONVA,
     $DLEVXA,DMYPRI,DXANRM,DXNRM,WNORM,EPDIFF,ETAMIN,ETAMAX,
     $ETAINI,ETADIF,FCDNM,FCMIN2,FCNUMP,FCCOR,FCNMP2,FCH,FCBND,
     $FCBH,FCK2,FCNUMK,FCREDU,DLEVFN,SUMXA,SUM1,SUM2,S1,TH,RSMALL,
     $APREC
      LOGICAL QGENJ,QINISC,QSUCC,QJCRFR,QLINIT,QSCALE,QNEXT,QREP,
     $        QRANK1,QMIXIO,QORDI,QSIMPL,QLU
CWEI
      INTRINSIC DLOG
      DOUBLE PRECISION CLIN0,CLIN1,CALPHA,CALPHK,ALPHAE,ALPHAK,ALPHAA,
     $                 SUMXA0,SUMXA1,SUMXA2,SUMXTE,FCMON,DLOG
      INTEGER IORMON
      LOGICAL QMSTOP
      SAVE CLIN0,CLIN1,CALPHA,ALPHAE,ALPHAK,ALPHAA,SUMXA0,SUMXA1,SUMXA2,
     $     QMSTOP,FCMON
C
      CALL ZIBCONST(EPMACH,SMALL)
C*    Begin
C       ----------------------------------------------------------
C       1 Initialization
C       ----------------------------------------------------------
C       1.1 Control-flags and -integers
        QSUCC = IOPT(1).EQ.1
        QSCALE = .NOT. IOPT(35).EQ.1
        QORDI  = IOPT(33).EQ.1
        QSIMPL = IOPT(34).EQ.1
        QRANK1 = IOPT(32).EQ.1
        IORMON = IOPT(39)
        IF (IORMON.EQ.0) IORMON=2
        ISCAL = IOPT(9)
        MODE = IOPT(2)
        JACGEN = IOPT(3)
        QMIXIO = LUMON.EQ.LUSOL .AND. MPRMON.NE.0 .AND. MPRSOL.NE.0
        QLU    = .NOT. QSIMPL
        MPRTIM = IOPT(19)
C       ----------------------------------------------------------
C       1.2 Derivated dimensional parameters
        IF (MSTOR.EQ.0) THEN
          ML=0
        ELSE IF (MSTOR.EQ.1) THEN
          ML=M1-M2
          MU=M2-1-ML
        ENDIF
C       ----------------------------------------------------------
C       1.3 Derivated internal parameters
        FCMIN2 = FCMIN*FCMIN
        RSMALL = DSQRT(TEN*RTOL)
C       ----------------------------------------------------------
C       1.4 Adaption of input parameters, if necessary
        IF(FC.LT.FCMIN) FC = FCMIN
        IF(FC.GT.ONE) FC = ONE
C       ----------------------------------------------------------
C       1.5 Initial preparations
        QJCRFR = .FALSE.
        QLINIT = .FALSE.
        IFAIL = 0
        FCBND = ZERO
        IF (QBDAMP) FCBND = RWK(20)
C       ----------------------------------------------------------
C       1.5.1 Numerical differentiation related initializations
        IF (JACGEN.EQ.2) THEN
          AJDEL = RWK(26)
          IF (AJDEL.LE.SMALL) AJDEL = DSQRT(EPMACH*TEN)
          AJMIN = RWK(27)
        ELSE IF (JACGEN.EQ.3) THEN
          ETADIF = RWK(28)
          IF (ETADIF .LE. SMALL) ETADIF = 1.0D-6
          ETAINI = RWK(29)
          IF (ETAINI .LE. SMALL) ETAINI = 1.0D-6
          EPDIFF = DSQRT(EPMACH*TEN)
          ETAMAX = DSQRT(EPDIFF)
          ETAMIN = EPDIFF*ETAMAX
        ENDIF
C       ----------------------------------------------------------
C       1.5.2 Miscellaneous preparations of first iteration step
        IF (.NOT.QSUCC) THEN
          NITER = 0
          NCORR = 0
          NREJR1 = 0
          NFCN = 0
          NJAC = 0
          NFCNJ = 0
          QGENJ = .TRUE.
          QINISC = .TRUE.
          FCKEEP = FC
          FCA = FC
          FCPRI = FC
          FCK2 = FC
          FCMON = FC
          CONV = ZERO
          IF (JACGEN.EQ.3) THEN
            DO 1520 L1=1,N
              ETA(L1)=ETAINI
1520        CONTINUE
          ENDIF
          DO 1521 L1=1,N
            XA(L1)=X(L1)
1521      CONTINUE
CWEI      
          ICONV = 0
          ALPHAE = ZERO
          SUMXA1 = ZERO
          SUMXA0 = ZERO
          CLIN0  = ZERO
          QMSTOP = .FALSE.
C         ------------------------------------------------------
C         1.6 Print monitor header
          IF(MPRMON.GE.2 .AND. .NOT.QMIXIO)THEN
16003       FORMAT(///,2X,66('*'))
            WRITE(LUMON,16003)
16004       FORMAT(/,8X,'It',7X,'Normf ',10X,'Normx ',8X,
     $             'Damp.Fct.',3X,'New')
            WRITE(LUMON,16004)
          ENDIF
C         --------------------------------------------------------
C         1.7 Startup step
C         --------------------------------------------------------
C         1.7.1 Computation of the residual vector
          IF (MPRTIM.NE.0) CALL MONON(1)
          CALL FCN(N,X,F,IFAIL)
          IF (MPRTIM.NE.0) CALL MONOFF(1)
          NFCN = NFCN+1
C     Exit, if ...
          IF (IFAIL.NE.0) THEN
            IERR = 82
            GOTO 4299
          ENDIF
        ELSE
          QINISC = .FALSE.
        ENDIF
C
C       Main iteration loop
C       ===================
C
C       Repeat
2       CONTINUE
C         --------------------------------------------------------
C         2 Startup of iteration step
          IF (.NOT.QJCRFR) THEN
C           ------------------------------------------------------
C           2.1 Scaling of variables X(N)
            CALL N1SCAL(N,X,XA,XSCAL,XW,ISCAL,QINISC,IOPT,LRWK,RWK)
            QINISC = .FALSE.
            IF(NITER.NE.0)THEN
C             ----------------------------------------------------
C             2.2 Aposteriori estimate of damping factor
              DO 2200 L1=1,N
                DXQA(L1)=DXQ(L1)
2200          CONTINUE
              IF (.NOT.QORDI) THEN
                FCNUMP = ZERO
                DO 2201 L1=1,N
                  FCNUMP=FCNUMP+(DX(L1)/XW(L1))**2
2201            CONTINUE
                TH = FC-ONE
                FCDNM = ZERO
                DO 2202 L1=1,N
                  FCDNM=FCDNM+((DXQA(L1)+TH*DX(L1))/XW(L1))**2
2202            CONTINUE
C               --------------------------------------------------
C               2.2.2 Decision criterion for Jacobian updating
C                     technique:
C                     QGENJ.EQ..TRUE. numerical differentation,
C                     QGENJ.EQ..FALSE. rank1 updating
                QGENJ = .TRUE.
                IF (FC.EQ.FCPRI) THEN
                  QGENJ = FC.LT.ONE.OR.FCA.LT.ONE.OR.DMYCOR.LE.FC*SIGMA
     $                    .OR. .NOT.QRANK1 .OR. NEW+2.GT.NBROY 
                  FCA = FC
                ELSE
                  DMYCOR = FCA*FCA*HALF*DSQRT(FCNUMP/FCDNM)
                  IF (NONLIN.LE.3) THEN
                    FCCOR = DMIN1(ONE,DMYCOR)
                  ELSE
                    FCCOR = DMIN1(ONE,HALF*DMYCOR)
                  ENDIF
                  FCA = DMAX1(DMIN1(FC,FCCOR),FCMIN)
C$Test-begin
                  IF (MPRMON.GE.5) THEN
                    WRITE(LUMON,22201) FCCOR, FC, DMYCOR, FCNUMP,
     $                                 FCDNM
22201               FORMAT (/, ' +++ aposteriori estimate +++', /,
     $                    ' FCCOR  = ', D18.10, '  FC     = ', D18.10, /,
     $                    ' DMYCOR = ', D18.10, '  FCNUMP = ', D18.10, /,
     $                    ' FCDNM  = ', D18.10, /,
     $                       ' ++++++++++++++++++++++++++++', /)
                  ENDIF
C$Test-end 
                ENDIF
                FCK2 = FCA
C               ------------------------------------------------------
C               2.2.1 Computation of the numerator of damping
C                     factor predictor
                FCNMP2 = ZERO
                DO 221 L1=1,N
                  FCNMP2=FCNMP2+(DXQA(L1)/XW(L1))**2
221             CONTINUE
                FCNUMP = FCNUMP*FCNMP2
              ENDIF
            ENDIF
          ENDIF
          QJCRFR =.FALSE.
C         --------------------------------------------------------
C         2.3 Jacobian matrix (stored to array A(M1,N))
C         --------------------------------------------------------
C         2.3.1 Jacobian generation by routine JAC or
C               difference approximation (If QGENJ.EQ..TRUE.)
C               - or -
C               Rank-1 update of Jacobian (If QGENJ.EQ..FALSE.)
          IF (QGENJ .AND. (.NOT.QSIMPL .OR. NITER.EQ.0)) THEN
            NEW = 0
            IF (JACGEN.EQ.1) THEN
               IF (MPRTIM.NE.0) CALL MONON(2)
               CALL JAC(N,M1,X,A,IFAIL)
               IF (MPRTIM.NE.0) CALL MONOFF(2)
            ELSE
              IF (MSTOR.EQ.0) THEN
                IF (MPRTIM.NE.0) CALL MONON(2)
                IF (JACGEN.EQ.3) 
     $            CALL N1JCF(FCN,N,N,X,F,A,XW,ETA,ETAMIN,ETAMAX,
     $                       ETADIF,CONV,NFCNJ,T1,IFAIL)
                IF (JACGEN.EQ.2) 
     $            CALL N1JAC(FCN, N, N, X, F, A, XW, AJDEL, AJMIN,
     $                       NFCNJ, T1, IFAIL)
                IF (MPRTIM.NE.0) CALL MONOFF(2)
              ELSE IF (MSTOR.EQ.1) THEN
                IF (MPRTIM.NE.0) CALL MONON(2)
                IF (JACGEN.EQ.3) 
     $            CALL N1JCFB(FCN,N,M1,ML,X,F,A,XW,ETA,ETAMIN,
     $                        ETAMAX,ETADIF,CONV,NFCNJ,T1,T2,T3,IFAIL)
                IF (JACGEN.EQ.2) 
     $            CALL N1JACB (FCN, N, M1, ML, X, F, A, XW,  AJDEL,
     $                         AJMIN, NFCNJ, T1, T2, T3, IFAIL)
                IF (MPRTIM.NE.0) CALL MONOFF(2)
              ENDIF
            ENDIF
            NJAC = NJAC + 1
C     Exit, If ...
            IF (JACGEN.EQ.1 .AND. IFAIL.LT.0) THEN
              IERR = 83
              GOTO 4299
            ENDIF
            IF (JACGEN.NE.1 .AND. IFAIL.NE.0) THEN
              IERR = 82
              GOTO 4299
            ENDIF
          ELSE IF (.NOT.QSIMPL) THEN
            NEW = NEW+1
          ENDIF
          IF ( NEW.EQ.0 .AND. (QLU.OR.NITER.EQ.0) ) THEN
C             ------------------------------------------------------
C             2.3.2.1 Save scaling values
              DO 2321 L1=1,N
                XWA(L1) = XW(L1)
2321          CONTINUE
C             ------------------------------------------------------
C             2.3.2.2 Prepare Jac. for use by band-solver DGBFA/DGBSL
              IF (MSTOR.EQ.1) THEN
                DO 2322 L1=1,N
                  DO 2323 L2=M2,1,-1
                    A(L2+ML,L1)=A(L2,L1)
2323              CONTINUE
2322            CONTINUE
              ENDIF
C             ------------------------------------------------------
C             2.4 Prepare solution of the linear system
C             ------------------------------------------------------
C             2.4.1 internal column scaling of matrix A
              IF (MSTOR.EQ.0) THEN
                DO 2410 K=1,N
                  S1 =-XW(K)
                  DO 2412 L1=1,N
                    A(L1,K)=A(L1,K)*S1
2412              CONTINUE
2410            CONTINUE
              ELSE IF (MSTOR.EQ.1) THEN
                DO 2413 K=1,N
                  L2=MAX0(1+M2-K,ML+1)
                  L3=MIN0(N+M2-K,M1)
                  S1 =-XW(K)
                  DO 2415 L1=L2,L3
                    A(L1,K)=A(L1,K)*S1
2415              CONTINUE
2413            CONTINUE
              ENDIF
C             ----------------------------------------------------
C             2.4.2 Row scaling of matrix A
              IF (QSCALE) THEN
                IF (MSTOR.EQ.0) THEN
                  CALL N1SCRF(N,N,A,FW)
                ELSE IF (MSTOR.EQ.1) THEN
                  CALL N1SCRB(N,M1,ML,MU,A,FW)
                ENDIF
              ELSE
                DO 242 K=1,N
                  FW(K)=ONE
242             CONTINUE
              ENDIF
          ENDIF
C         --------------------------------------------------------
C         2.4.3 Save and scale values of F(N)
          DO 243 L1=1,N
            FA(L1)=F(L1)
            T1(L1)=F(L1)*FW(L1)
243       CONTINUE
C         --------------------------------------------------------
C         3 Central part of iteration step
C         --------------------------------------------------------
C         3.1 Solution of the linear system
C         --------------------------------------------------------
C         3.1.1 Decomposition of (N,N)-matrix A
          IF (.NOT.QLINIT) THEN
            NIWLA = IWK(18)+1
            NRWLA = IWK(19)+1
            LIWL = LIWK-NIWLA+1
            LRWL = LRWK-NRWLA+1
          ENDIF
          IF ( NEW.EQ.0 .AND. (QLU.OR.NITER.EQ.0)) THEN
            IF (MPRTIM.NE.0) CALL MONON(3)
            CALL N1FACT(N,M1,ML,MU,A,IOPT,IFAIL,LIWL,IWK(NIWLA),
     $                  NILUSE,LRWL,RWK(NRWLA),NRLUSE)
            IF (MPRTIM.NE.0) CALL MONOFF(3)
            IF (.NOT.QLINIT) THEN
              NIWKFR = NIWKFR+NILUSE
              NRWKFR = NRWKFR+NRLUSE
C             Store lengths of currently required workspaces
              IWK(18) = NIWKFR-1
              IWK(19) = NRWKFR-1
            ENDIF
C       Exit Repeat If ...
            IF(IFAIL.NE.0) THEN
              IF (IFAIL.EQ.1) THEN
                IERR = 1
              ELSE
                IERR = 80
              ENDIF
              GOTO 4299
            ENDIF
          ENDIF
          QLINIT = .TRUE.
C         --------------------------------------------------------
C         3.1.2 Solution of linear (N,N)-system
          IF(NEW.EQ.0) THEN 
            IF (MPRTIM.NE.0) CALL MONON(4)
            CALL N1SOLV(N,M1,ML,MU,A,T1,IOPT,IFAIL,LIWL,IWK(NIWLA),
     $                 IDUMMY,LRWL,RWK(NRWLA),IDUMMY)
            IF (MPRTIM.NE.0) CALL MONOFF(4)
C     Exit Repeat If ...
            IF(IFAIL.NE.0)  THEN
              IERR = 81
              GOTO 4299
            ENDIF
          ELSE  
            ALFA1=ZERO
            ALFA2=ZERO
            DO 3121 I=1,N
              ALFA1=ALFA1+DX(I)*DXQ(I)/XW(I)**2
              ALFA2=ALFA2+DX(I)**2/XW(I)**2
3121        CONTINUE
            ALFA=ALFA1/ALFA2
            BETA=ONE-ALFA
            DO 3122 I=1,N
              T1(I)=(DXQ(I)+(FCA-ONE)*ALFA*DX(I))/BETA
3122        CONTINUE
            IF(NEW.EQ.1) THEN
              DO 3123 I=1,N
                DXSAVE(I,1)=DX(I)
3123          CONTINUE
            ENDIF
            DO 3124 I=1,N
              DXSAVE(I,NEW+1)=T1(I)
              DX(I)=T1(I)
              T1(I)=T1(I)/XW(I)
3124        CONTINUE
          ENDIF
C         --------------------------------------------------------
C         3.2 Evaluation of scaled natural level function SUMX
C             scaled maximum error norm CONV
C             evaluation of (scaled) standard level function
C             DLEVF ( DLEVF only, if MPRMON.GE.2 )
C             and computation of ordinary Newton corrections 
C             DX(N)
          IF (.NOT. QSIMPL) THEN
            CALL N1LVLS(N,T1,XW,F,DX,CONV,SUMX,DLEVF,MPRMON,NEW.EQ.0)
          ELSE
            CALL N1LVLS(N,T1,XWA,F,DX,CONV,SUMX,DLEVF,MPRMON,NEW.EQ.0)
          ENDIF
          DO 32 L1=1,N
            XA(L1)=X(L1)
32        CONTINUE
          SUMXA = SUMX
          DLEVXA = DSQRT(SUMXA/DBLE(FLOAT(N)))
          CONVA = CONV
          DXANRM = WNORM(N,DX,XW)
C         --------------------------------------------------------
C         3.3 A - priori estimate of damping factor FC
          IF(NITER.NE.0.AND.NONLIN.NE.1.AND.NEW.EQ.0.AND.
     $       .NOT. QORDI)THEN
C           ------------------------------------------------------
C           3.3.1 Computation of the denominator of a-priori
C                 estimate
            FCDNM = ZERO
            DO 331 L1=1,N
              FCDNM=FCDNM+((DX(L1)-DXQA(L1))/XW(L1))**2
331         CONTINUE
            FCDNM = FCDNM*SUMX
C           ------------------------------------------------------
C           3.3.2 New damping factor
            IF(FCDNM.GT.FCNUMP*FCMIN2 .OR.
     $        (NONLIN.EQ.4 .AND. FCA**2*FCNUMP .LT. 4.0D0*FCDNM)) THEN
              DMYPRI = FCA*DSQRT(FCNUMP/FCDNM)
              FCPRI = DMIN1(DMYPRI,ONE)
              IF (NONLIN.EQ.4) FCPRI = DMIN1(HALF*DMYPRI,ONE)
            ELSE
              FCPRI = ONE
C$Test-begin
              DMYPRI = -1.0D0
C$Test-end
            ENDIF
C$Test-begin
            IF (MPRMON.GE.5) THEN
              WRITE(LUMON,33201) FCPRI, FC, FCA, DMYPRI, FCNUMP,
     $                           FCDNM
33201         FORMAT (/, ' +++ apriori estimate +++', /,
     $                ' FCPRI  = ', D18.10, '  FC     = ', D18.10, /,
     $                ' FCA    = ', D18.10, '  DMYPRI = ', D18.10, /,
     $                ' FCNUMP = ', D18.10, '  FCDNM  = ', D18.10, /,
     $                   ' ++++++++++++++++++++++++', /)
            ENDIF
C$Test-end 
            FC = DMAX1(FCPRI,FCMIN)
            IF (QBDAMP) THEN
              FCBH = FCA*FCBND
              IF (FC.GT.FCBH) THEN
                FC = FCBH
                IF (MPRMON.GE.4)
     $            WRITE(LUMON,*) ' *** incr. rest. act. (a prio) ***'
              ENDIF
              FCBH = FCA/FCBND
              IF (FC.LT.FCBH) THEN
                FC = FCBH
                IF (MPRMON.GE.4)
     $            WRITE(LUMON,*) ' *** decr. rest. act. (a prio) ***'
              ENDIF
            ENDIF
          ENDIF
CWEI
          IF (IORMON.GE.2) THEN
            SUMXA2=SUMXA1
            SUMXA1=SUMXA0
            SUMXA0=DLEVXA
            IF (SUMXA0.EQ.ZERO) SUMXA0=SMALL
C           Check convergence rates (linear and superlinear)
C           ICONV : Convergence indicator
C                   =0: No convergence indicated yet
C                   =1: Damping factor is 1.0d0
C                   =2: Superlinear convergence detected (alpha >=1.2)
C                   =3: Quadratic convergence detected (alpha > 1.8)
            FCMON = DMIN1(FC,FCMON)
            IF (FCMON.LT.ONE) THEN
              ICONV = 0
              ALPHAE = ZERO
            ENDIF
            IF (FCMON.EQ.ONE .AND. ICONV.EQ.0) ICONV=1
            IF (NITER.GE.1) THEN
              CLIN1 = CLIN0
              CLIN0 = SUMXA0/SUMXA1
            ENDIF
            IF (ICONV.GE.1.AND.NITER.GE.2) THEN
              ALPHAK = ALPHAE
              ALPHAE = ZERO
              IF (CLIN1.LE.0.95D0) ALPHAE = DLOG(CLIN0)/DLOG(CLIN1)
              IF (ALPHAK.NE.ZERO) ALPHAK =0.5D0*(ALPHAE+ALPHAK)
              ALPHAA = DMIN1(ALPHAK,ALPHAE)
              CALPHK = CALPHA
              CALPHA = ZERO
              IF (ALPHAE.NE.ZERO) CALPHA = SUMXA1/SUMXA2**ALPHAE
              SUMXTE = DSQRT(CALPHA*CALPHK)*SUMXA1**ALPHAK-SUMXA0
              IF (ALPHAA.GE.1.2D0 .AND. ICONV.EQ.1) ICONV = 2
              IF (ALPHAA.GT.1.8D0) ICONV = 3
              IF (MPRMON.GE.4)  WRITE(LUMON,32001) ICONV, ALPHAE, 
     $                            CALPHA, CLIN0, ALPHAK, SUMXTE
32001         FORMAT(' ** ICONV: ',I1,'  ALPHA: ',D9.2,
     $               '  CONST-ALPHA: ',D9.2,'  CONST-LIN: ',D9.2,' **',
     $               /,' **',11X,'ALPHA-POST: ',D9.2,' CHECK: ',D9.2,
     $               25X,'**')
              IF ( ICONV.GE.2 .AND. ALPHAA.LT.0.9D0 ) THEN
                 IF (IORMON.EQ.3) THEN
                   IERR = 4
                   GOTO 4299
                 ELSE
                   QMSTOP = .TRUE.
                 ENDIF 
              ENDIF
            ENDIF
          ENDIF
          FCMON = FC
C
C         --------------------------------------------------------
C         3.4 Save natural level for later computations of
C             corrector and print iterate
          FCNUMK = SUMX
          IF (MPRMON.GE.2) THEN
            IF (MPRTIM.NE.0) CALL MONON(5)
            CALL N1PRV1(DLEVF,DLEVXA,FCKEEP,NITER,NEW,MPRMON,LUMON,
     $                  QMIXIO)
            IF (MPRTIM.NE.0) CALL MONOFF(5)
          ENDIF
          NRED = 0
          QNEXT = .FALSE.
          QREP  = .FALSE.   
C         QREP = ITER .GT. ITMAX   or  QREP = ITER .GT. 0
C
C         Damping-factor reduction loop
C         ================================
C         DO (Until)
34        CONTINUE
C           ------------------------------------------------------
C           3.5 Preliminary new iterate
            DO 35 L1=1,N
              X(L1)=XA(L1)+DX(L1)*FC
35          CONTINUE
C           -----------------------------------------------------
C           3.5.2 Exit, if problem is specified as being linear
C     Exit Repeat If ...
            IF( NONLIN.EQ.1 )THEN
              IERR = 0
              GOTO 4299
            ENDIF
C           ------------------------------------------------------
C           3.6.1 Computation of the residual vector
            IF (MPRTIM.NE.0) CALL MONON(1)
            CALL FCN(N,X,F,IFAIL)
            IF (MPRTIM.NE.0) CALL MONOFF(1)
            NFCN = NFCN+1
C     Exit, If ...
            IF(IFAIL.LT.0)THEN
              IERR = 82
              GOTO 4299
            ENDIF
            IF(IFAIL.EQ.1 .OR. IFAIL.EQ.2) THEN
              IF (QORDI) THEN
                IERR = 82
                GOTO 4299
              ENDIF
              IF (IFAIL.EQ.1) THEN
                FCREDU = HALF
              ELSE
                FCREDU = F(1)
C     Exit, If ...
                IF (FCREDU.LE.0 .OR. FCREDU.GE.1) THEN
                  IERR = 82
                  GOTO 4299
                ENDIF
              ENDIF
              IF (MPRMON.GE.2) THEN
36101           FORMAT(8X,I2,' FCN could not be evaluated  ',
     $                 8X,F7.5,4X,I2)
                WRITE(LUMON,36101)NITER,FC,NEW
              ENDIF
              FCH = FC
              FC = FCREDU*FC
              IF (FCH.GT.FCMIN) FC = DMAX1(FC,FCMIN)
              IF (QBDAMP) THEN
                FCBH = FCH/FCBND
                IF (FC.LT.FCBH) THEN
                  FC = FCBH
                  IF (MPRMON.GE.4) WRITE(LUMON,*)
     $               ' *** decr. rest. act. (FCN redu.) ***'
                ENDIF
              ENDIF
              IF (FC.LT.FCMIN) THEN
                IERR = 3
                GOTO 4299
              ENDIF  
C     Break DO (Until) ...
              GOTO 3109
            ENDIF
            IF (QORDI) THEN
C             -------------------------------------------------------
C             3.6.2 Convergence test for ordinary Newton iteration
C     Exit Repeat If ...
              IF( DXANRM.LE.RTOL )THEN
                IERR = 0
                GOTO 4299
              ENDIF
            ELSE
              DO 361 L1=1,N
               T1(L1)=F(L1)*FW(L1)
361           CONTINUE
C             --------------------------------------------------
C             3.6.3 Solution of linear (N,N)-system
              IF (MPRTIM.NE.0) CALL MONON(4)
              CALL N1SOLV(N,M1,ML,MU,A,T1,IOPT,IFAIL,LIWL,IWK(NIWLA),
     $                   IDUMMY,LRWL,RWK(NRWLA),IDUMMY)
              IF (MPRTIM.NE.0) CALL MONOFF(4)
C     Exit Repeat If ...
              IF(IFAIL.NE.0)  THEN
                IERR = 81
                GOTO 4299
              ENDIF
              IF(NEW.GT.0) THEN 
                DO 3630 I=1,N
                  DXQ(I) = T1(I)*XWA(I)
3630            CONTINUE                   
                DO 363 ILOOP=1,NEW 
                  SUM1=ZERO
                  SUM2=ZERO
                  DO 3631 I=1,N
                    SUM1=SUM1+(DXQ(I)*DXSAVE(I,ILOOP))/ XW(I)**2
                    SUM2=SUM2+(DXSAVE(I,ILOOP)/XW(I))**2
3631              CONTINUE
                  BETA=SUM1/SUM2
                  DO 3632 I=1,N
                    DXQ(I)=DXQ(I)+BETA*DXSAVE(I,ILOOP+1)
                    T1(I) = DXQ(I)/XW(I)
3632              CONTINUE
363             CONTINUE
              ENDIF
C             ----------------------------------------------------
C             3.6.4 Evaluation of scaled natural level function
C                   SUMX
C                   scaled maximum error norm CONV and evaluation
C                   of (scaled) standard level function DLEVFN
              IF (.NOT. QSIMPL) THEN
                CALL N1LVLS(N,T1,XW,F,DXQ,CONV,SUMX,DLEVFN,MPRMON,
     $                      NEW.EQ.0)
              ELSE
                CALL N1LVLS(N,T1,XWA,F,DXQ,CONV,SUMX,DLEVFN,MPRMON,
     $                      NEW.EQ.0)
              ENDIF  
              DXNRM = WNORM(N,DXQ,XW)
C             -----------------------------------------------------
C             3.6.5 Convergence test
C     Exit Repeat If ...
              IF ( DXNRM.LE.RTOL .AND. DXANRM.LE.RSMALL .AND. 
     $            FC.EQ.ONE ) THEN
                IERR = 0
                GOTO 4299
              ENDIF
C           
              FCA = FC
C             ----------------------------------------------------
C             3.6.6 Evaluation of reduced damping factor
              TH = FCA-ONE
              FCDNM = ZERO
              DO 39 L1=1,N
                FCDNM=FCDNM+((DXQ(L1)+TH*DX(L1))/XW(L1))**2
39            CONTINUE
              IF (FCDNM.NE.ZERO) THEN
                DMYCOR = FCA*FCA*HALF*DSQRT(FCNUMK/FCDNM)
              ELSE
                DMYCOR = 1.0D+35
              ENDIF
              IF (NONLIN.LE.3) THEN
                FCCOR = DMIN1(ONE,DMYCOR)
              ELSE
                FCCOR = DMIN1(ONE,HALF*DMYCOR)
              ENDIF
C$Test-begin
              IF (MPRMON.GE.5) THEN
                WRITE(LUMON,39001) FCCOR, FC, DMYCOR, FCNUMK,
     $                             FCDNM, FCA
39001           FORMAT (/, ' +++ corrector computation +++', /,
     $                ' FCCOR  = ', D18.10, '  FC     = ', D18.10, /,
     $                ' DMYCOR = ', D18.10, '  FCNUMK = ', D18.10, /,
     $                ' FCDNM  = ', D18.10, '  FCA    = ', D18.10, /,
     $                   ' +++++++++++++++++++++++++++++', /)
              ENDIF
C$Test-end 
            ENDIF
C           ------------------------------------------------------
C           3.7 Natural monotonicity test
            IF(SUMX.GT.SUMXA .AND. .NOT.QORDI)THEN
C             ----------------------------------------------------
C             3.8 Output of iterate
              IF(MPRMON.GE.3) THEN
                IF (MPRTIM.NE.0) CALL MONON(5)
                CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,
     $                      NITER,MPRMON,LUMON,QMIXIO,'*')
                IF (MPRTIM.NE.0) CALL MONOFF(5)
              ENDIF
              IF (QMSTOP) THEN
                IERR = 4
                GOTO 4299
              ENDIF
              FCH = DMIN1(FCCOR,HALF*FC)
              IF (FC.GT.FCMIN) THEN
                FC=DMAX1(FCH,FCMIN)
              ELSE
                FC=FCH
              ENDIF
              IF (QBDAMP) THEN
                FCBH = FCA/FCBND
                IF (FC.LT.FCBH) THEN
                  FC = FCBH
                  IF (MPRMON.GE.4)
     $              WRITE(LUMON,*) ' *** decr. rest. act. (a post) ***'
                ENDIF
              ENDIF
CWEI
              FCMON = FC
C
C$Test-begin
                IF (MPRMON.GE.5) THEN
                  WRITE(LUMON,39002) FC
39002             FORMAT (/, ' +++ corrector setting 1 +++', /,
     $                    ' FC     = ', D18.10, /,
     $                       ' +++++++++++++++++++++++++++', /)
                ENDIF
C$Test-end 
              QREP = .TRUE.
              NCORR = NCORR+1
              NRED = NRED+1
C             ----------------------------------------------------
C             3.10 If damping factor is too small:
C                  Refresh Jacobian,if current Jacobian was computed
C                  by a Rank1-update, else fail exit
              QJCRFR  = FC.LT.FCMIN.OR.NEW.GT.0.AND.NRED.GT.1
C     Exit Repeat If ...
              IF(QJCRFR.AND.NEW.EQ.0)THEN
                IERR = 3
                GOTO 4299
              ENDIF
            ELSE
              IF (.NOT.QORDI.AND..NOT.QREP .AND. FCCOR.GT.SIGMA2*FC)
     $        THEN
                IF(MPRMON.GE.3) THEN
                  IF (MPRTIM.NE.0) CALL MONON(5)
                  CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,
     $                        NITER,MPRMON,LUMON,QMIXIO,'+')
                  IF (MPRTIM.NE.0) CALL MONOFF(5)
                ENDIF
                FC = FCCOR
C$Test-begin
                IF (MPRMON.GE.5) THEN
                  WRITE(LUMON,39003) FC
39003             FORMAT (/, ' +++ corrector setting 2 +++', /,
     $                    ' FC     = ', D18.10, /,
     $                       ' +++++++++++++++++++++++++++', /)
                ENDIF
C$Test-end 
                QREP = .TRUE.
              ELSE
                QNEXT = .TRUE.
              ENDIF
            ENDIF
3109      CONTINUE
          IF(.NOT.(QNEXT.OR.QJCRFR)) GOTO  34
C         UNTIL ( expression - negated above)
C         End of damping-factor reduction loop
C         =======================================
          IF(QJCRFR)THEN
C           ------------------------------------------------------
C           3.11 Restore former values for repeting iteration
C                step
            NREJR1 = NREJR1+1
            DO 3111 L1=1,N
              X(L1)=XA(L1)
3111        CONTINUE
            DO 3112 L1=1,N
              F(L1)=FA(L1)
3112        CONTINUE
            IF(MPRMON.GE.2)THEN
31130           FORMAT(8X,I2,' Not accepted damping factor ',
     $                 8X,F7.5,4X,I2)
                WRITE(LUMON,31130)NITER,FC,NEW
            ENDIF
            FC = FCKEEP
            FCA = FCK2
            IF(NITER.EQ.0)THEN
              FC = FCMIN
            ENDIF
            QGENJ = .TRUE.
          ELSE
C           ------------------------------------------------------
C           4 Preparations to start the following iteration step
C           ------------------------------------------------------
C           4.1 Print values
            IF(MPRMON.GE.3 .AND. .NOT.QORDI) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,NITER+1,
     $                    MPRMON,LUMON,QMIXIO,'*')
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
C           Print the natural level of the current iterate and return
C           it in one-step mode
            SUMXS = SUMX
            SUMX = SUMXA
            IF(MPRSOL.GE.2.AND.NITER.NE.0) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N1SOUT(N,XA,2,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ELSE IF(MPRSOL.GE.1.AND.NITER.EQ.0)THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N1SOUT(N,XA,1,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
            NITER = NITER+1
            DLEVF = DLEVFN
C     Exit Repeat If ...
            IF(NITER.GE.NITMAX)THEN
              IERR = 2
              GOTO 4299
            ENDIF
            FCKEEP = FC
C           ------------------------------------------------------
C           4.2 Return, if in one-step mode
C
C Exit Subroutine If ...
            IF (MODE.EQ.1) THEN
              IWK(18)=NIWLA-1
              IWK(19)=NRWLA-1
              IOPT(1)=1
              RETURN
            ENDIF
          ENDIF
        GOTO 2
C       End Repeat
4299    CONTINUE
C       End of main iteration loop
C       ==========================
C       ----------------------------------------------------------
C       9 Exits
C       ----------------------------------------------------------
C       9.1 Solution exit
        APREC = -1.0D0
C 
        IF(IERR.EQ.0 .OR. IERR.EQ.4)THEN
          IF (NONLIN.NE.1) THEN
            IF (.NOT.QORDI) THEN
              IF ( IERR.EQ.0 ) THEN
                APREC = DSQRT(SUMX/DBLE(FLOAT(N)))
                DO 91 L1=1,N
                  X(L1)=X(L1)+DXQ(L1)
91              CONTINUE
              ELSE 
                APREC = DSQRT(SUMXA/DBLE(FLOAT(N)))
                IF (ALPHAA.GT.ZERO .AND. IORMON.EQ.3) THEN
                  DO 92 L1=1,N
                    X(L1)=X(L1)+DX(L1)
92                CONTINUE
                ENDIF
              ENDIF
C             Print final monitor output
              IF(MPRMON.GE.2) THEN
                IF (IERR.EQ.0) THEN
                  IF (MPRTIM.NE.0) CALL MONON(5)
                  CALL N1PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,
     $                        NITER+1,MPRMON,LUMON,QMIXIO,'*')
                  IF (MPRTIM.NE.0) CALL MONOFF(5)
                ELSE IF (IORMON.EQ.3) THEN
                  IF (MPRTIM.NE.0) CALL MONON(5)
                  CALL N1PRV1(DLEVFN,DSQRT(SUMXA/DBLE(FLOAT(N))),FC,
     $                        NITER,NEW,MPRMON,LUMON,QMIXIO)
                  IF (MPRTIM.NE.0) CALL MONOFF(5)
                ENDIF
              ENDIF
              IF (  IORMON.GE.2 ) THEN
                IF ( ICONV.LE.1 .AND. ALPHAE .NE. ZERO 
     $                          .AND. ALPHAK .NE. ZERO ) IERR = 5
              ENDIF
            ELSE
C           IF (QORDI) THEN
              APREC = DSQRT(SUMXA/DBLE(FLOAT(N)))
            ENDIF
            IF(MPRMON.GE.1) THEN
91001         FORMAT(///' Solution of nonlinear system ',
     $        'of equations obtained within ',I3,
     $        ' iteration steps',//,' Achieved relative accuracy',D10.3)
              IF (QORDI .OR. IERR.EQ.4) THEN
                WRITE(LUMON,91001) NITER,APREC
              ELSE
                WRITE(LUMON,91001) NITER+1,APREC
              ENDIF 
            ENDIF
          ELSE
            IF(MPRMON.GE.1) THEN
91002         FORMAT(///' Solution of linear system ',
     $        'of equations obtained by NLEQ1',//,' No estimate ',
     $        'available for the achieved relative accuracy')
                WRITE(LUMON,91002)
            ENDIF
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       9.2 Fail exit messages
C       ----------------------------------------------------------
C       9.2.1 Termination, since jacobian matrix became singular
        IF(IERR.EQ.1.AND.MPRERR.GE.1)THEN
92101     FORMAT(/,' Iteration terminates due to ',
     $    'singular jacobian matrix',/)
          WRITE(LUERR,92101)
        ENDIF
C       ----------------------------------------------------------
C       9.2.2 Termination after more than NITMAX iterations
        IF(IERR.EQ.2.AND.MPRERR.GE.1)THEN
92201     FORMAT(/,' Iteration terminates after NITMAX ',
     $    '=',I3,'  Iteration steps')
          WRITE(LUERR,92201)NITMAX
        ENDIF
C       ----------------------------------------------------------
C       9.2.3 Damping factor FC became too small
        IF(IERR.EQ.3.AND.MPRERR.GE.1)THEN
92301     FORMAT(/,' Damping factor has become too ',
     $    'small: lambda =',D10.3,2X,/)
          WRITE(LUERR,92301)FC
        ENDIF
CWEI
C       ----------------------------------------------------------
C       9.2.4.1 Superlinear convergence slowed down
        IF(IERR.EQ.4.AND.MPRERR.GE.1)THEN
92401     FORMAT(/,' Warning: Monotonicity test failed after ',A,
     $           ' convergence was already checked;',/,
     $    ' RTOL requirement may be too stringent',/)
92402     FORMAT(/,' Warning: ',A,' convergence slowed down;',/,
     $    ' RTOL requirement may be too stringent',/)
          IF (QMSTOP) THEN
            IF (ICONV.EQ.2) WRITE(LUERR,92401) 'superlinear'
            IF (ICONV.EQ.3) WRITE(LUERR,92401) 'quadratic'
          ELSE
            IF (ICONV.EQ.2) WRITE(LUERR,92402) 'superlinear'
            IF (ICONV.EQ.3) WRITE(LUERR,92402) 'quadratic'
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       9.2.4.2 Convergence criterion satisfied before superlinear
C               convergence has been established
        IF (IERR.EQ.5.AND.DLEVFN.EQ.ZERO) IERR=0
        IF(IERR.EQ.5.AND.MPRERR.GE.1)THEN
92410     FORMAT(/,' Warning: No quadratic or superlinear convergence ',
     $           'established yet',/,
     $           10X,'your solution may perhaps may be less accurate ',
     $           /,10X,'as indicated by the standard error estimate')
          WRITE(LUERR,92410)
        ENDIF
C       ----------------------------------------------------------
C       9.2.5 Error exit due to linear solver routine N1FACT
        IF(IERR.EQ.80.AND.MPRERR.GE.1)THEN
92501     FORMAT(/,' Error ',I5,' signalled by linear solver N1FACT')
          WRITE(LUERR,92501) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.6 Error exit due to linear solver routine N1SOLV
        IF(IERR.EQ.81.AND.MPRERR.GE.1)THEN
92601     FORMAT(/,' Error ',I5,' signalled by linear solver N1SOLV')
          WRITE(LUERR,92601) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.7 Error exit due to fail of user function FCN
        IF(IERR.EQ.82.AND.MPRERR.GE.1)THEN
92701     FORMAT(/,' Error ',I5,' signalled by user function FCN')
          WRITE(LUERR,92701) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.8 Error exit due to fail of user function JAC
        IF(IERR.EQ.83.AND.MPRERR.GE.1)THEN
92801     FORMAT(/,' Error ',I5,' signalled by user function JAC')
          WRITE(LUERR,92801) IFAIL
        ENDIF
        IF(IERR.GE.80.AND.IERR.LE.83) IWK(23) = IFAIL
        IF ((IERR.EQ.82.OR.IERR.EQ.83).AND.NITER.LE.1.AND.MPRERR.GE.1)
     $  THEN
          WRITE (LUERR,92810)
92810     FORMAT(' Try to find a better initial guess for the solution')
        ENDIF
C       ----------------------------------------------------------
C       9.3 Common exit
        IF (MPRERR.GE.3.AND.IERR.NE.0.AND.IERR.NE.4.AND.NONLIN.NE.1)
     $    THEN
93100     FORMAT(/,'    Achieved relative accuracy',D10.3,2X)
          WRITE(LUERR,93100)CONVA
          APREC = CONVA
        ENDIF
        RTOL = APREC
        SUMX = SUMXA
        IF(MPRSOL.GE.2.AND.NITER.NE.0) THEN
           MODE=2
           IF (QORDI) MODE=3
           IF (MPRTIM.NE.0) CALL MONON(5)
           CALL N1SOUT(N,XA,MODE,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
           IF (MPRTIM.NE.0) CALL MONOFF(5)
        ELSE IF(MPRSOL.GE.1.AND.NITER.EQ.0)THEN
           IF (MPRTIM.NE.0) CALL MONON(5)
           CALL N1SOUT(N,XA,1,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
           IF (MPRTIM.NE.0) CALL MONOFF(5)
        ENDIF
        IF (.NOT.QORDI) THEN
          IF (IERR.NE.4) NITER = NITER+1
          DLEVF = DLEVFN
          IF(MPRSOL.GE.1)THEN
C           Print Solution or final iteration vector
            IF(IERR.EQ.0)THEN
               MODEFI = 3
            ELSE
               MODEFI = 4
            ENDIF
            IF (MPRTIM.NE.0) CALL MONON(5)
            CALL N1SOUT(N,X,MODEFI,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
            IF (MPRTIM.NE.0) CALL MONOFF(5)
          ENDIF
        ENDIF
C       Return the latest internal scaling to XSCAL
        DO 93 I=1,N
          XSCAL(I)=XW(I)
93      CONTINUE
C       End of exits
C       End of subroutine N1INT
      RETURN
      END
C
      SUBROUTINE N1SCAL(N,X,XA,XSCAL,XW,ISCAL,QINISC,IOPT,LRWK,RWK)
C*    Begin Prologue SCALE
      INTEGER N
      DOUBLE PRECISION X(N),XSCAL(N),XA(N),XW(N)
      INTEGER ISCAL
      LOGICAL QINISC
      INTEGER IOPT(50),LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*    Summary :
C    
C     S C A L E : To be used in connection with NLEQ1 .
C       Computation of the internal scaling vector XW used for the
C       Jacobian matrix, the iterate vector and it's related
C       vectors - especially for the solution of the linear system
C       and the computations of norms to avoid numerical overflow.
C
C*    Input parameters
C     ================
C
C     N         Int     Number of unknowns
C     X(N)      Dble    Current iterate
C     XA(N)     Dble    Previous iterate
C     XSCAL(N)  Dble    User scaling passed from parameter XSCAL
C                       of interface routine NLEQ1
C     ISCAL     Int     Option ISCAL passed from IOPT-field
C                       (for details see description of IOPT-fields)
C     QINISC    Logical = .TRUE.  : Initial scaling
C                       = .FALSE. : Subsequent scaling
C     IOPT(50)  Int     Options array passed from NLEQ1 parameter list
C     LRWK      Int     Length of real workspace
C     RWK(LRWK) Dble    Real workspace (see description above)
C
C*    Output parameters
C     =================
C
C     XW(N)     Dble   Scaling vector computed by this routine
C                      All components must be positive. The follow-
C                      ing relationship between the original vector
C                      X and the scaled vector XSCAL holds:
C                      XSCAL(I) = X(I)/XW(I) for I=1,...N
C
C*    Subroutines called: ZIBCONST
C
C*    Machine constants used
C     ======================
C
      DOUBLE PRECISION EPMACH, SMALL
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS,DMAX1
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      INTEGER MPRMON,LUMON
      CALL ZIBCONST(EPMACH, SMALL)
C*    Begin
      DO 1 L1=1,N
        IF (ISCAL.EQ.1) THEN
          XW(L1) = XSCAL(L1)
        ELSE
          XW(L1)=DMAX1(XSCAL(L1),(DABS(X(L1))+DABS(XA(L1)))*HALF,SMALL)
        ENDIF
1     CONTINUE
C$Test-Begin
      MPRMON = IOPT(13)
      IF (MPRMON.GE.6) THEN
        LUMON = IOPT(14)
        WRITE(LUMON,*) ' '
        WRITE(LUMON,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(LUMON,*) '      X-components   Scaling-components    '
        WRITE(LUMON,10) (X(L1), XW(L1), L1=1,N)
10      FORMAT('  ',D18.10,'  ',D18.10)
        WRITE(LUMON,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(LUMON,*) ' '
      ENDIF
C$Test-End
C     End of subroutine N1SCAL
      RETURN
      END
C
      SUBROUTINE N1SCRF(M,N,A,FW)
C*    Begin Prologue SCROWF
      INTEGER M,N
      DOUBLE PRECISION A(M,N),FW(M)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S C R O W F : Row Scaling of a (M,N)-matrix in full storage
C                   mode
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       M           Int    Number of rows of the matrix
C       N           Int    Number of columns of the matrix
C     * A(M,N)      Dble   Matrix to be scaled
C
C*    Output parameters
C     =================
C
C       FW(M)       Dble   Row scaling factors - FW(i) contains
C                          the factor by which the i-th row of A
C                          has been multiplied
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER J,K
      DOUBLE PRECISION S1,S2
C*    Begin
      DO 1 K=1,M
        S1=ZERO
        DO 2 J=1,N
          S2=DABS(A(K,J))
          IF (S2.GT.S1) S1=S2
2       CONTINUE
        IF (S1.GT.ZERO) THEN
          S1=ONE/S1
          FW(K)=S1
          DO 3 J=1,N
            A(K,J)=A(K,J)*S1
3         CONTINUE
        ELSE
          FW(K)=ONE
        ENDIF
1     CONTINUE
C     End of subroutine N1SCRF
      RETURN
      END
C
      SUBROUTINE N1SCRB(N,LDA,ML,MU,A,FW)
C*    Begin Prologue SCROWB
      INTEGER N,LDA,ML,MU
      DOUBLE PRECISION A(LDA,N),FW(N)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S C R O W B : Row Scaling of a (N,N)-matrix in band storage
C                   mode
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       N           Int    Number of rows and columns of the matrix
C       LDA         Int    Leading dimension of the matrix array
C       ML          Int    Lower bandwidth of the matrix (see IOPT)
C       MU          Int    Upper bandwidth of the matrix (see IOPT)
C     * A(LDA,N)    Dble   Matrix to be scaled (see Note in routine
C                          DGBFA for details on storage technique)
C
C*    Output parameters
C     =================
C
C       FW(N)       Dble   Row scaling factors - FW(i) contains
C                          the factor by which the i-th row of
C                          the matrix has been multiplied
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS,MAX0,MIN0
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER K,K1,L1,L2,L3,M2
      DOUBLE PRECISION S1,S2
C*    Begin
      M2=ML+MU+1
      DO 1 K=1,N
        S1=ZERO
        L2=MAX0(1,K-ML)
        L3=MIN0(N,K+MU)
        K1=M2+K
        DO 2 L1=L2,L3
          S2=DABS(A(K1-L1,L1))
          IF (S2.GT.S1) S1=S2
2       CONTINUE
        IF (S1.GT.ZERO) THEN
          S1=ONE/S1
          FW(K)=S1
          DO 3 L1=L2,L3
            A(K1-L1,L1)=A(K1-L1,L1)*S1
3         CONTINUE
        ELSE
          FW(K)=ONE
        ENDIF
1     CONTINUE
C     End of subroutine N1SCRB
      RETURN
      END
C
      SUBROUTINE N1FACT(N,LDA,ML,MU,A,IOPT,IFAIL,LIWK,IWK,LAIWK,LRWK,
     $RWK,LARWK)
C*    Begin Prologue FACT
      INTEGER N,LDA,ML,MU
      DOUBLE PRECISION A(LDA,N)
      INTEGER IOPT(50)
      INTEGER IFAIL
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LAIWK,LRWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER LARWK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     F A C T : Call linear algebra subprogram for factorization of
C               a (N,N)-matrix
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C     N             Int    Order of the linear system
C     LDA           Int    Leading dimension of the matrix array A
C     ML            Int    Lower bandwidth of the matrix (only for
C                          banded systems)
C     MU            Int    Upper bandwidth of the matrix (only for
C                          banded systems)
C   * A(LDA,N)      Dble   Matrix storage. See main routine NLEQ1,
C                          note 4 for how to store a banded Jacobian
C     IOPT(50)      Int    Option vector passed from NLEQ1
C
C*    Output parameters
C     =================
C
C     IFAIL         Int    Error indicator returned by this routine:
C                          = 0 matrix decomposition successfull
C                          = 1 decomposition failed - 
C                              matrix is numerically singular
C                          =10 supplied (integer) workspace too small
C
C*    Workspace parameters
C     ====================
C
C     LIWK          Int    Length of integer workspace passed to this
C                          routine (In)
C     IWK(LIWK)     Int    Integer Workspace supplied for this routine
C     LAIWK         Int    Length of integer Workspace used by this 
C                          routine (out)       
C     LRWK          Int    Length of real workspace passed to this
C                          routine (In)                  
C     RWK(LRWK)     Dble   Real Workspace supplied for this routine
C     LARWK         Int    Length of real Workspace used by this 
C                          routine (out)
C
C*    Subroutines called:  DGETRF, DGBTRF
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL DGETRF, DGBTRF
C*    Begin
      LAIWK = N
      LARWK = 0
C     ( Real workspace starts at RWK(NRWKFR) )
C     IF (LIWK.GE.LAIWK.AND.LRWK.GE.LARWK) THEN
      IF (LIWK.GE.LAIWK) THEN
        MSTOR = IOPT(4)
        IF (MSTOR.EQ.0) THEN
          CALL DGETRF(N,N,A,LDA,IWK,IFAIL)
        ELSE IF (MSTOR.EQ.1) THEN
          CALL DGBTRF(N,N,ML,MU,A,LDA,IWK,IFAIL)
        ENDIF
        IF(IFAIL.NE.0)THEN
          IFAIL = 1
        ENDIF
      ELSE
        IFAIL = 10
        MPRERR=IOPT(11)
        LUERR=IOPT(12)
10001   FORMAT(/,' Insuffient workspace for linear solver,',
     $         ' at least needed more needed : ',/,
     $         ' ',A,' workspace : ',I4)
        IF (LIWK.LT.LAIWK.AND.MPRERR.GT.0) 
     $    WRITE(LUERR,10001) 'Integer',LAIWK-LIWK
      ENDIF
      RETURN
      END
C
      SUBROUTINE N1SOLV(N,LDA,ML,MU,A,B,IOPT,IFAIL,LIWK,IWK,LAIWK,
     $LRWK,RWK,LARWK)
C*    Begin Prologue SOLVE
      INTEGER N,LDA,ML,MU
      DOUBLE PRECISION A(LDA,N)
      DOUBLE PRECISION B(N)
      INTEGER IOPT(50)
      INTEGER IFAIL
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK,LAIWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER LARWK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S O L V E : Call linear algebra subprogram for solution of
C                 the linear system A*Z = B
C
C*    Parameters
C     ==========
C
C     N,LDA,ML,MU,A,IOPT,IFAIL,LIWK,IWK,LAIWK,LRWK,RWK,LARWK :
C                        See description for subroutine N1FACT.          
C     B          Dble    In:  Right hand side of the linear system
C                        Out: Solution of the linear system
C
C     Subroutines called: DGETRS, DGBTRS
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL DGETRS, DGBTRS
C*    Begin
      MSTOR = IOPT(4)
      IF (MSTOR.EQ.0) THEN
        CALL DGETRS('N',N,1,A,LDA,IWK,B,N,IFAIL)
      ELSE IF (MSTOR.EQ.1) THEN
        CALL DGBTRS('N',N,ML,MU,1,A,LDA,IWK,B,N,IFAIL)
      ENDIF
      RETURN
      END
C
      SUBROUTINE N1LVLS(N,DX1,XW,F,DXQ,CONV,SUMX,DLEVF,MPRMON,QDSCAL)
C*    Begin Prologue LEVELS
      INTEGER N,MPRMON
      DOUBLE PRECISION DX1(N),XW(N),F(N),DXQ(N)
      DOUBLE PRECISION CONV,SUMX,DLEVF
      LOGICAL QDSCAL
C
C     ------------------------------------------------------------
C
C*    Summary :
C
C     L E V E L S : To be used in connection with NLEQ1 .
C     provides descaled solution, error norm and level functions
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       N              Int    Number of parameters to be estimated
C       DX1(N)         Dble   array containing the scaled Newton
C                             correction
C       XW(N)          Dble   Array containing the scaling values
C       F(N)           Dble   Array containing the residuum
C
C*    Output parameters
C     =================
C
C       DXQ(N)         Dble   Array containing the descaled Newton
C                             correction
C       CONV           Dble   Scaled maximum norm of the Newton
C                             correction
C       SUMX           Dble   Scaled natural level function value
C       DLEVF          Dble   Standard level function value (only
C                             if needed for print)
C       MPRMON         Int    Print information parameter (see
C                             driver routine NLEQ1 )
C       QDSCAL         Logic  .TRUE., if descaling of DX1 required,
C                             else .FALSE.
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1
      DOUBLE PRECISION S1
C*    Begin
      IF (QDSCAL) THEN
C       ------------------------------------------------------------
C       1.2 Descaling of solution DX1 ( stored to DXQ )
        DO 12 L1=1,N
          DXQ(L1)=DX1(L1)*XW(L1)
12      CONTINUE
      ENDIF
C     ------------------------------------------------------------
C     2 Evaluation of scaled natural level function SUMX and
C       scaled maximum error norm CONV
      CONV = ZERO
      DO 20 L1=1,N
        S1 = DABS(DX1(L1))
        IF(S1.GT.CONV) CONV=S1
20    CONTINUE
      SUMX = ZERO
      DO 21 L1=1,N
        SUMX = SUMX+DX1(L1)**2
21    CONTINUE
C     ------------------------------------------------------------
C     3 Evaluation of (scaled) standard level function DLEVF
      DLEVF = ZERO
      DO 3 L1=1,N
        DLEVF = DLEVF+F(L1)**2
3     CONTINUE
      DLEVF = DSQRT(DLEVF/DBLE(FLOAT(N)))
C     End of subroutine N1LVLS
      RETURN
      END
C
      SUBROUTINE N1JAC (FCN, N, LDA, X, FX, A, YSCAL, AJDEL, AJMIN,
     $                  NFCN, FU, IFAIL)
C* Begin Prologue N1JAC
      EXTERNAL FCN
      INTEGER N, LDA
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), AJDEL, AJMIN
      INTEGER NFCN
      DOUBLE PRECISION FU(N)
      INTEGER IFAIL
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Evaluation of a dense Jacobian matrix using finite difference
C  approximation adapted for use in nonlinear systems solver NLEQ1
C
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   May 1990
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C  FCN        Ext     FCN (N, X, FX, IFAIL)
C                     Subroutine in order to provide right-hand
C                     side of first-order differential equations
C    N        Int     Number of rows and columns of the Jacobian
C    X(N)     Dble    The current scaled iterates
C    FX(N)    Dble    Array containing FCN(X)
C    IFAIL    Int     Return code
C                     Whenever a negative value is returned by FCN
C                     routine N1JAC is terminated immediately.
C
C
C* Input parameters (* marks inout parameters)
C  ----------------
C
C  N          Int     Number of rows and columns of the Jacobian
C  LDA        Int     Leading Dimension of array A
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(N)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  AJDEL      Dble    Perturbation of component k: abs(Y(k))*AJDEL
C  AJMIN      Dble    Minimum perturbation is AJMIN*AJDEL
C  NFCN       Int  *  FCN - evaluation count
C
C* Output parameters (* marks inout parameters)
C  -----------------
C
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C
C* Workspace parameters
C  --------------------
C
C  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C
C* Called
C  ------
C
      INTRINSIC DABS, DMAX1, DSIGN
C  ---------------------------------------------------------------------
C
C* End Prologue
C
C* Local variables
C  ---------------
C
      INTEGER I, K
      DOUBLE PRECISION U, W
C
C* Begin
C
      IFAIL = 0
      DO 1 K = 1,N
         W = X(K)
         U = DSIGN(DMAX1(DABS(X(K)),AJMIN,YSCAL(K))*AJDEL, X(K))
         X(K) = W + U
C
         CALL FCN (N, X, FU, IFAIL)
         NFCN = NFCN + 1
         IF (IFAIL .NE. 0) GOTO 99
C
         X(K) = W
         DO 11 I = 1,N
            A(I,K) = (FU(I) - FX(I)) / U  
 11      CONTINUE
 1    CONTINUE
C
99    CONTINUE
      RETURN
C
C
C* End of N1JAC
C
      END
      SUBROUTINE N1JACB (FCN, N, LDA, ML, X, FX, A, YSCAL, AJDEL, AJMIN,
     $     NFCN, FU, U, W, IFAIL)
C* Begin Prologue N1JACB
      EXTERNAL FCN
      INTEGER N, LDA, ML
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), AJDEL, AJMIN
      INTEGER NFCN
      DOUBLE PRECISION FU(N), U(N), W(N)
      INTEGER IFAIL
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Evaluation of a banded Jacobian matrix using finite difference
C  approximation adapted for use in nonlinear systems solver NLEQ1
C
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   May 1990
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C  FCN        Ext     FCN (N, X, FX, IFAIL)
C                     Subroutine in order to provide right-hand
C                     side of first-order differential equations
C    N        Int     Number of rows and columns of the Jacobian
C    X(N)     Dble    The current scaled iterates
C    FX(N)    Dble    Array containing FCN(X)
C    IFAIL    Int     Return code
C                     Whenever a negative value is returned by FCN
C                     routine N1JACB is terminated immediately.
C
C
C* Input parameters (* marks inout parameters)
C  ----------------
C
C  N          Int     Number of rows and columns of the Jacobian
C  LDA        Int     Leading Dimension of array A
C  ML         Int     Lower bandwidth of the Jacobian matrix
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(N)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  AJDEL      Dble    Perturbation of component k: abs(Y(k))*AJDEL
C  AJMIN      Dble    Minimum perturbation is AJMIN*AJDEL
C  NFCN       Int  *  FCN - evaluation count
C
C* Output parameters (* marks inout parameters)
C  -----------------
C
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C
C* Workspace parameters
C  --------------------
C
C  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C  U(N)       Dble    Array to contain dx(i)
C  W(N)       Dble    Array to save original components of X
C
C* Called
C  ------
C
      INTRINSIC DABS, DMAX1, MAX0, MIN0
C
C* Constants
C  ---------
C
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.0 D0)
C
C  ---------------------------------------------------------------------
C
C* End Prologue
C
C* Local variables
C  ---------------
C
      INTEGER I, I1, I2, JJ, K, MH, MU
C
C* Begin
C
      IFAIL = 0
      MU = LDA - 2*ML - 1
      LDAB = ML + MU + 1
      DO 1 I = 1,LDAB
         DO 11 K = 1,N
            A(I,K) = ZERO
 11      CONTINUE
 1    CONTINUE
C
      DO 2 JJ = 1,LDAB
         DO 21 K = JJ,N,LDAB
            W(K) = X(K)
            U(K) = DSIGN(DMAX1(DABS(X(K)),AJMIN,YSCAL(K))*AJDEL, X(K))
            X(K) = W(K) + U(K)
 21      CONTINUE
C
         CALL FCN (N, X, FU, IFAIL)
         NFCN = NFCN + 1
         IF (IFAIL .NE. 0) GOTO 99
C
         DO 22 K = JJ,N,LDAB 
            X(K) = W(K)
            I1 = MAX0 (1, K - MU)
            I2 = MIN0 (N, K + ML)
            MH = MU + 1 - K
            DO 221 I = I1,I2
               A(MH+I,K) = (FU(I) - FX(I)) / U(K)
 221        CONTINUE
 22      CONTINUE
C
 2    CONTINUE
C
99    CONTINUE
      RETURN
C
C
C* End of N1JACB
C
      END
      SUBROUTINE N1JCF (FCN, N, LDA, X, FX, A, YSCAL, ETA, ETAMIN,
     $     ETAMAX, ETADIF, CONV, NFCN, FU, IFAIL)
C* Begin Prologue N1JCF
      EXTERNAL FCN
      INTEGER N, LDA
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), ETA(N),
     $     ETAMIN, ETAMAX, ETADIF, CONV
      INTEGER NFCN
      DOUBLE PRECISION FU(N)
      INTEGER IFAIL
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Approximation of dense Jacobian matrix for nonlinear systems solver
C  NLEQ1 with feed-back control of discretization and rounding errors
C
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   May 1990
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C  FCN        Ext     FCN (N, X, FX, IFAIL)
C                     Subroutine in order to provide right-hand
C                     side of first-order differential equations
C    N        Int     Number of rows and columns of the Jacobian
C    X(N)     Dble    The current scaled iterates
C    FX(N)    Dble    Array containing FCN(X)
C    IFAIL    Int     Return code
C                     Whenever a negative value is returned by FCN
C                     routine N1JCF is terminated immediately.
C
C
C* Input parameters (* marks inout parameters)
C  ----------------
C
C  N          Int     Number of rows and columns of the Jacobian
C  LDA        Int     Leading dimension of A (LDA .GE. N)
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(N)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  ETA(N)     Dble *  Array containing the scaled denominator
C                     differences
C  ETAMIN     Dble    Minimum allowed scaled denominator
C  ETAMAX     Dble    Maximum allowed scaled denominator
C  ETADIF     Dble    DSQRT (1.1*EPMACH)
C                     EPMACH = machine precision
C  CONV       Dble    Maximum norm of last (unrelaxed) Newton correction
C  NFCN       Int  *  FCN - evaluation count
C
C* Output parameters (* marks inout parameters)
C  -----------------
C
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  ETA(N)     Dble *  Scaled denominator differences adjusted
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C
C* Workspace parameters
C  --------------------
C
C  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C
C* Called
C  ------
C
      INTRINSIC DABS, DMAX1, DMIN1, DSIGN, DSQRT
C
C* Constants
C  ---------
C
      DOUBLE PRECISION SMALL2, ZERO
      PARAMETER (SMALL2 = 0.1D0,
     $           ZERO   = 0.0D0)
C
C  ---------------------------------------------------------------------
C
C* End Prologue
C
C* Local variables
C  ---------------
C
      INTEGER I, K, IS
      DOUBLE PRECISION FHI, HG, U, SUMD, W
      LOGICAL QFINE
C
C* Begin
C
      DO 1 K = 1,N
         IS = 0
C        DO (Until)
 11         CONTINUE
            W = X(K)
            U = DSIGN (ETA(K)*YSCAL(K), X(K))
            X(K) = W + U
            CALL FCN (N, X, FU, IFAIL)
            NFCN = NFCN + 1
C           Exit, If ...
            IF (IFAIL .NE. 0) GOTO 99
            X(K) = W
            SUMD = ZERO
            DO 111 I = 1,N
               HG = DMAX1 (DABS (FX(I)), DABS (FU(I)))
               FHI = FU(I) - FX(I)
               IF (HG .NE. ZERO) SUMD = SUMD + (FHI/HG)**2
               A(I,K) = FHI / U
 111        CONTINUE
            SUMD = DSQRT (SUMD / DBLE(N))
            QFINE = .TRUE.
            IF (SUMD .NE. ZERO .AND. IS .EQ. 0)THEN
               ETA(K) = DMIN1 (ETAMAX,
     $              DMAX1 (ETAMIN, DSQRT (ETADIF / SUMD)*ETA(K)))
               IS = 1
               QFINE = CONV .LT. SMALL2 .OR. SUMD .GE. ETAMIN
            ENDIF
            IF (.NOT.(QFINE)) GOTO  11
C        UNTIL ( expression - negated above)
 1    CONTINUE
C
C     Exit from DO-loop
 99   CONTINUE
C
      RETURN
C
C* End of subroutine N1JCF
C
      END
      SUBROUTINE N1JCFB (FCN, N, LDA, ML, X, FX, A, YSCAL, ETA,
     $     ETAMIN, ETAMAX, ETADIF, CONV, NFCN, FU, U, W, IFAIL)
C* Begin Prologue N1JCFB
      EXTERNAL FCN
      INTEGER N, LDA, ML
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), ETA(N),
     $     ETAMIN, ETAMAX, ETADIF, CONV
      INTEGER NFCN
      DOUBLE PRECISION FU(N), U(N), W(N)
      INTEGER IFAIL
C     
C     ---------------------------------------------------------------------
C     
C* Title
C  
C  Approximation of banded Jacobian matrix for nonlinear systems solver
C  NLEQ1 with feed-back control of discretization and rounding errors
C  
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   May 1990
C  
C  
C* Parameter list description
C  --------------------------
C  
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C  
C  FCN        Ext     FCN (N, X, FX, IFAIL)
C                     Subroutine in order to provide right-hand
C                     side of first-order differential equations
C    N        Int     Number of rows and columns of the Jacobian
C    X(N)     Dble    The current scaled iterates
C    FX(N)    Dble    Array containing FCN(X)
C    IFAIL    Int     Return code
C                     Whenever a negative value is returned by FCN
C                     routine N1JCFB is terminated immediately.
C  
C  
C* Input parameters (* marks inout parameters)
C  ----------------
C  
C  N          Int     Number of rows and columns of the Jacobian
C  LDA        Int     Leading dimension of A (LDA .GE. ML+MU+1)
C  ML         Int     Lower bandwidth of Jacobian matrix
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(N)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  ETA(N)     Dble *  Array containing the scaled denominator
C                     differences
C  ETAMIN     Dble    Minimum allowed scaled denominator
C  ETAMAX     Dble    Maximum allowed scaled denominator
C  ETADIF     Dble    DSQRT (1.1*EPMACH)
C                     EPMACH = machine precision
C  CONV       Dble    Maximum norm of last (unrelaxed) Newton correction
C  NFCN       Int  *  FCN - evaluation count
C  
C* Output parameters (* marks inout parameters)
C  -----------------
C  
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  ETA(N)     Dble *  Scaled denominator differences adjusted
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C  
C* Workspace parameters
C  --------------------
C  
C  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C  U(N)       Dble    Array to contain dx(i)
C  W(N)       Dble    Array to save original components of X
C  
C* Called
C  ------
C  
      INTRINSIC DABS, DMAX1, DMIN1, DSIGN, DSQRT, MAX0, MIN0
C  
C* Constants
C  ---------
C  
      DOUBLE PRECISION SMALL2, ZERO
      PARAMETER (SMALL2 = 0.1D0,
     $           ZERO   = 0.0D0)
C  
C  ---------------------------------------------------------------------
C  
C* End Prologue
C  
C* Local variables
C  ---------------
C  
      INTEGER I, IS, I1, I2, JJ, K, MH, MU, LDAB
      DOUBLE PRECISION FHI, HG, SUMD
      LOGICAL QFINE
C  
C* Begin
C  
      MU = LDA - 2*ML - 1
      LDAB = ML + MU + 1
      DO 1 I = 1,LDAB
         DO 11 K = 1,N
            A(I,K) = ZERO
 11      CONTINUE
 1    CONTINUE
C  
      DO 2 JJ = 1,LDAB
         IS = 0
C        DO (Until)
 21         CONTINUE
            DO 211  K = JJ,N,LDAB
               W(K) = X(K)
               U(K) = DSIGN (ETA(K)*YSCAL(K), X(K))
               X(K) = W(K) + U(K)
 211        CONTINUE
C     
            CALL FCN (N, X, FU, IFAIL)
            NFCN = NFCN + 1
C           Exit, If ...
            IF (IFAIL .NE. 0) GOTO 99
C     
            DO 212 K = JJ,N,LDAB
               X(K) = W(K)
               SUMD = ZERO
               I1 = MAX0 (1, K - MU)
               I2 = MIN0 (N, K + ML)
               MH = MU + 1 - K
               DO 213 I = I1,I2
                  HG = DMAX1 (DABS(FX(I)), DABS(FU(I)))
                  FHI = FU(I) - FX(I)
                  IF (HG .NE. ZERO) SUMD = SUMD + (FHI/HG)**2
                  A(MH+I, K) = FHI / U(K)
 213           CONTINUE
               SUMD = DSQRT (SUMD / DBLE(N))
               QFINE = .TRUE.
               IF (SUMD .NE. ZERO .AND. IS .EQ. 0) THEN
                  ETA(K) = DMIN1 (ETAMAX,
     $                 DMAX1 (ETAMIN, DSQRT (ETADIF / SUMD)*ETA(K)))
                  IS = 1
                  QFINE = CONV .LT. SMALL2 .OR. SUMD .GE. ETAMIN
               ENDIF
 212        CONTINUE
            IF (.NOT.(QFINE)) GOTO 21
C        UNTIL ( expression - negated above)
 2    CONTINUE
C
C     Exit from DO-loop
 99   CONTINUE
C
      RETURN
C
C* End of subroutine N1JCFB
C
      END
C
      SUBROUTINE N1PRV1(DLEVF,DLEVX,FC,NITER,NEW,MPRMON,LUMON,QMIXIO)
C*    Begin Prologue N1PRV1
      DOUBLE PRECISION DLEVF,DLEVX,FC
      INTEGER NITER,NEW,MPRMON,LUMON
      LOGICAL QMIXIO
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 1 P R V 1 : Printing of intermediate values (Type 1 routine)
C
C*    Parameters
C     ==========
C
C     DLEVF, DLEVX   See descr. of internal double variables of N1INT
C     FC,NITER,NEW,MPRMON,LUMON
C                  See parameter descr. of subroutine N1INT
C     QMIXIO Logical  = .TRUE.  , if LUMON.EQ.LUSOL
C                     = .FALSE. , if LUMON.NE.LUSOL
C
C     ------------------------------------------------------------
C*    End Prologue
C     Print Standard - and natural level
      IF(QMIXIO)THEN
1       FORMAT(2X,66('*'))
        WRITE(LUMON,1)
2       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',20X,'New')
        IF (MPRMON.GE.3) WRITE(LUMON,2)
3       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',8X,'Damp.Fct.',3X,'New')
        IF (MPRMON.EQ.2) WRITE(LUMON,3)
      ENDIF
4     FORMAT(6X,I4,5X,D10.3,2X,4X,D10.3,17X,I2)
      IF (MPRMON.GE.3.OR.NITER.EQ.0) 
     $  WRITE(LUMON,4) NITER,DLEVF,DLEVX,NEW
5     FORMAT(6X,I4,5X,D10.3,6X,D10.3,6X,F7.5,4X,I2)
      IF (MPRMON.EQ.2.AND.NITER.NE.0) 
     $  WRITE(LUMON,5) NITER,DLEVF,DLEVX,FC,NEW
      IF(QMIXIO)THEN
6       FORMAT(2X,66('*'))
        WRITE(LUMON,6)
      ENDIF
C     End of subroutine N1PRV1
      RETURN
      END
C
      SUBROUTINE N1PRV2(DLEVF,DLEVX,FC,NITER,MPRMON,LUMON,QMIXIO,
     $                  CMARK)
C*    Begin Prologue N1PRV2
      DOUBLE PRECISION DLEVF,DLEVX,FC
      INTEGER NITER,MPRMON,LUMON
      LOGICAL QMIXIO
      CHARACTER*1 CMARK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 1 P R V 2 : Printing of intermediate values (Type 2 routine)
C
C*    Parameters
C     ==========
C
C     DLEVF, DLEVX   See descr. of internal double variables of N2INT
C     FC,NITER,MPRMON,LUMON
C                  See parameter descr. of subroutine N2INT
C     QMIXIO Logical  = .TRUE.  , if LUMON.EQ.LUSOL
C                     = .FALSE. , if LUMON.NE.LUSOL
C     CMARK Char*1    Marker character to be printed before DLEVX
C
C     ------------------------------------------------------------
C*    End Prologue
C     Print Standard - and natural level, and damping
C     factor
      IF(QMIXIO)THEN
1       FORMAT(2X,66('*'))
        WRITE(LUMON,1)
2       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',8X,'Damp.Fct.')
        WRITE(LUMON,2)
      ENDIF
3     FORMAT(6X,I4,5X,D10.3,4X,A1,1X,D10.3,2X,4X,F7.5)
      WRITE(LUMON,3)NITER,DLEVF,CMARK,DLEVX,FC
      IF(QMIXIO)THEN
4       FORMAT(2X,66('*'))
        WRITE(LUMON,4)
      ENDIF
C     End of subroutine N1PRV2
      RETURN
      END
C
      SUBROUTINE N1SOUT(N,X,MODE,IOPT,RWK,NRW,IWK,NIW,MPRINT,LUOUT)
C*    Begin Prologue SOLOUT
      INTEGER N
      DOUBLE PRECISION X(N)
      INTEGER NRW
      INTEGER MODE
      INTEGER IOPT(50)
      DOUBLE PRECISION RWK(NRW)
      INTEGER NIW
      INTEGER IWK(NIW)
      INTEGER MPRINT,LUOUT
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S O L O U T : Printing of iterate (user customizable routine)
C
C*    Input parameters
C     ================
C
C     N         Int Number of equations/unknowns
C     X(N)   Double iterate vector
C     MODE          =1 This routine is called before the first
C                      Newton iteration step
C                   =2 This routine is called with an intermedi-
C                      ate iterate X(N)
C                   =3 This is the last call with the solution
C                      vector X(N)
C                   =4 This is the last call with the final, but
C                      not solution vector X(N)
C     IOPT(50)  Int The option array as passed to the driver
C                   routine(elements 46 to 50 may be used
C                   for user options)
C     MPRINT    Int Solution print level 
C                   (see description of IOPT-field MPRINT)
C     LUOUT     Int the solution print unit 
C                   (see description of see IOPT-field LUSOL)
C
C
C*    Workspace parameters
C     ====================
C
C     NRW, RWK, NIW, IWK    see description in driver routine
C
C*    Use of IOPT by this routine
C     ===========================
C
C     Field 46:       =0 Standard output
C                     =1 GRAZIL suitable output
C
C     ------------------------------------------------------------
C*    End Prologue
      LOGICAL QGRAZ,QNORM
C*    Begin
      QNORM = IOPT(46).EQ.0
      QGRAZ = IOPT(46).EQ.1
      IF(QNORM) THEN
1        FORMAT('  ',A,' data:',/)
         IF (MODE.EQ.1) THEN
101        FORMAT('  Start data:',/,'  N =',I5,//,
     $            '  Format: iteration-number, (x(i),i=1,...N), ',
     $            'Normf , Normx ',/)
           WRITE(LUOUT,101) N
           WRITE(LUOUT,1) 'Initial'
         ELSE IF (MODE.EQ.3) THEN
           WRITE(LUOUT,1) 'Solution'
         ELSE IF (MODE.EQ.4) THEN
           WRITE(LUOUT,1) 'Final'
         ENDIF
2        FORMAT(' ',I5)
C        WRITE          NITER
         WRITE(LUOUT,2) IWK(1)
3        FORMAT((12X,3(D18.10,1X)))
         WRITE(LUOUT,3)(X(L1),L1=1,N)
C        WRITE          DLEVF,   DLEVX
         WRITE(LUOUT,3) RWK(19),DSQRT(RWK(18)/DBLE(FLOAT(N)))
         IF(MODE.EQ.1.AND.MPRINT.GE.2) THEN
           WRITE(LUOUT,1) 'Intermediate'
         ELSE IF(MODE.GE.3) THEN
           WRITE(LUOUT,1) 'End'
         ENDIF
      ENDIF
      IF(QGRAZ) THEN
        IF(MODE.EQ.1) THEN
10        FORMAT('&name com',I3.3,:,255(7(', com',I3.3,:),/))
          WRITE(LUOUT,10)(I,I=1,N+2)
15        FORMAT('&def  com',I3.3,:,255(7(', com',I3.3,:),/))
          WRITE(LUOUT,15)(I,I=1,N+2)
16        FORMAT(6X,': X=1, Y=',I3)
          WRITE(LUOUT,16) N+2
        ENDIF
20      FORMAT('&data ',I5)
C        WRITE          NITER
        WRITE(LUOUT,20) IWK(1) 
21      FORMAT((6X,4(D18.10)))
        WRITE(LUOUT,21)(X(L1),L1=1,N)
C        WRITE          DLEVF,   DLEVX
        WRITE(LUOUT,21) RWK(19),DSQRT(RWK(18)/DBLE(FLOAT(N)))
        IF(MODE.GE.3) THEN
30        FORMAT('&wktype 3111',/,'&atext x ''iter''')
          WRITE(LUOUT,30)
35        FORMAT('&vars = com',I3.3,/,'&atext y ''x',I3,'''',
     $           /,'&run')
          WRITE(LUOUT,35) (I,I,I=1,N)
36        FORMAT('&vars = com',I3.3,/,'&atext y ''',A,'''',
     $           /,'&run')
          WRITE(LUOUT,36) N+1,'Normf ',N+2,'Normx '
C39       FORMAT('&stop')
C         WRITE(LUOUT,39)
        ENDIF
      ENDIF
C     End of subroutine N1SOUT
      RETURN
      END
C*    End package
