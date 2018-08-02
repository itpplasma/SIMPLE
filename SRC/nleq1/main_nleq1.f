      PROGRAM TNLEQ1
      IMPLICIT DOUBLEPRECISION(S)
C
C     ____________________________________________________________
C
C     Testexample for NLEQ1: Computation of the Tchebychev
C     Polynom's evaluation points for quadrature-formulas.
C
C*  Written by        L. Weimann 
C*  Purpose           Testexample for code NLEQ1
C*  Version           2.3
C*  Revision          January 1992
C*  Latest Change     June 1992
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C
C     ____________________________________________________________
C
      INTEGER IRW
      PARAMETER (IRW=400)
      INTEGER IIW
      PARAMETER (IIW=60)
      INTEGER NN
      PARAMETER (NN=13)
      INTEGER NMAXP,N,I,N1
      DOUBLE PRECISION EPS
      INTEGER IOPT(50)
      INTEGER IERR, IFAIL
      DOUBLE PRECISION X(NN),XSCAL(NN),RW(IRW)
      INTEGER IW(IIW)
      REAL  STIME, ETIME, CPTIME
      EXTERNAL F
      EXTERNAL DF
C:    Begin
      OPEN(2,FILE='nleq1.dat')
      OPEN(9,FILE='nleq1.out')
      WRITE(6,*) ' monitor: nleq1.out , data: nleq1.dat '
      NMAXP = 9
      N = 2
C:    While (expression)
70    IF(N.LE.NMAXP)THEN
        EPS = 1.0D-5
        N1 = N+1
        DO 710 I=1,50
          IOPT(I)=0
710      CONTINUE
        DO 711 I=1,IIW
          IW(I)=0
711     CONTINUE
        DO 712 I=1,IRW
          RW(I)=0.0D0
712     CONTINUE
C       Execution mode: 0=Standard Mode, 1=Stepwise mode
        IOPT(2)=1
C       Jacobian: 0=(same as value 3)
C                 1=supplied by user routine JAC
C                 2=computed by numerical differentation (no feedback) 
C                 3=computed by numerical differentation (with feedback)
        IOPT(3)=1
C       Jacobian storage mode: 0=full matrix, 1=band matrix
        IOPT(4)=0
C       For a band matrix Jacobian: lower bandwidth
        IOPT(6)=0
C       For a band matrix Jacobian: upper bandwidth
        IOPT(7)=0
C       Broyden updates: 0 = inhibit, 1=allow
C       IOPT(32)=1
C     1 = linear , 2 = mildly nonlinear , 3 = highly nonlinear
C     4 = extremely nonlinear
C       IOPT(31)=3
C       Set MPRERR, MPRMON, MPRSOL, MPRTIM
        IOPT(11)=3
        IOPT(13)=3
        IOPT(15)=2
        IOPT(19)=1
C       Set print units LUERR, LUMON, LUSOL, LUTIM
        IOPT(12)=9
        IOPT(14)=9
        IOPT(16)=2
        IOPT(20)=9
C       Solution output format:
C       0=standard format, 1= GRAZIL readable output
        IOPT(46)=0
C       Override maximum allowed number of iterations:
        IW(31)=10
C       Override starting damping factor:
C       RW(21)=1.0D0
C       Override minimal allowed damping factor:
C       RW(22)=1.0D-2
C       Override rank1-decision parameter SIGMA:
C       RW(23)=3.0D0
C       Override 'take corrector instead of small predictor'-
C       decision parameter SIGMA2:
C       RW(24)=3.0D0
C       Starting vector for the iteration
        DO 72 I=1,N
          X(I)=DBLE(I)/DBLE(N1)
72      CONTINUE
        DO 75 I=1,N
          XSCAL(I) = 0.0
75      CONTINUE
        IERR=-1
        I=0
        CALL ZIBSEC(STIME,IFAIL)
31      IF (IERR.EQ.-1) THEN
          CALL NLEQ1(N,F,DF,X,XSCAL,EPS,IOPT,IERR,IIW,IW,IRW,RW)
C         Clear workspace declared not to be used
          NIFREE=IW(16)
          DO 311 K=NIFREE,IIW
            IW(K)=0
311       CONTINUE
          NRFREE=IW(17)
          DO 312 K=NRFREE,IRW
            RW(K)=0.0D0
312       CONTINUE
          I=I+1
32        FORMAT(' Returned from call ',I4,' of NLEQ1')
          WRITE(9,32)I
C         IOPT(2)=0
          GOTO 31
        ENDIF
        CALL ZIBSEC(ETIME,IFAIL)
        CPTIME = ETIME-STIME
73      FORMAT(//,1X,'Time ','used ','=',F9.3,1X,'Sec',//,66('*'),
     *    /)
        WRITE(9,73)CPTIME
80      N = N+1
      GOTO 70
      ENDIF
C.    EndWhile
      END
      SUBROUTINE F(N,X,FX,IFLAG)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N
      DOUBLE PRECISION X(N),FX(N)
C:    End Parameter
      INTEGER I,I1,L
      DOUBLE PRECISION TI2,TI1,TI,FACTT
C:    Begin
      DO 73 I=2,N,2
        I1 = I-1
        FX(I1)=0.0
        FX(I)=DBLE(N)/DBLE(I*I-1)
73    CONTINUE
      IF(MOD(N,2).EQ.1)THEN
        FX(N)=0.0
      ENDIF
      DO 74 L=1,N
        FACTT = 4.0*X(L)-2.0
        TI2 = 1.0
        TI1 = 0.5*FACTT
        FX(1)=TI1+FX(1)
        DO 75 I=2,N
          TI = FACTT*TI1-TI2
          FX(I)=TI+FX(I)
          TI2 = TI1
          TI1 = TI
75      CONTINUE
74    CONTINUE
      RETURN
      END
      SUBROUTINE DF(N,LDDFX,X,DFX,IFLAG)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N,LDDFX
      DOUBLE PRECISION X(N)
      DOUBLE PRECISION DFX(LDDFX,N)  
C:    End Parameter
      INTEGER I,J
      DOUBLE PRECISION TI2,TI1,TI,FACTT,TABLI2,TABLI1,TABLI
C:    Begin
      DO 76 J=1,N
        FACTT = 4.0*X(J)-2.0
        TI2 = 1.0
        TI1 = 0.5*FACTT
        TABLI2 = 0.0
        TABLI1 = 2.0
        DFX(1,J)=TABLI1
        DO 77 I=2,N
          TI = FACTT*TI1-TI2
          TABLI = 4.0*TI1+FACTT*TABLI1-TABLI2
          DFX(I,J)=TABLI
          TI2 = TI1
          TI1 = TI
          TABLI2 = TABLI1
          TABLI1 = TABLI
77      CONTINUE
76    CONTINUE
      RETURN
      END
