      PROGRAM TNLE1E
C
C     ____________________________________________________________
C
C     Testexample for NLEQ1E: Computation of the Tchebychev
C     Polynom's evaluation points for quadrature-formulas.
C
C*  Written by        L. Weimann 
C*  Purpose           Testexample for easy to use driver NLEQ1E
C*  Version           2.3
C*  Revision          November 1991
C*  Latest Change     November 1991
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C
C     ____________________________________________________________
C
      INTEGER NN
      PARAMETER ( NN=9 )
      INTEGER N,I,N1
      DOUBLE PRECISION RTOL
      INTEGER IERR, IFAIL
      DOUBLE PRECISION X(NN)
      EXTERNAL NLEQ1E, SECOND
      REAL STIME,ETIME,CPTIME
      N = 2
10    IF(N.LE.NN)THEN
        RTOL = 1.0D-5
        N1 = N+1
        DO 20 I=1,N
          X(I)=DBLE(I)/DBLE(N1)
20      CONTINUE
        CALL ZIBSEC(STIME,IFAIL)
        CALL NLEQ1E(N,X,RTOL,IERR)
        CALL ZIBSEC(ETIME,IFAIL)
        CPTIME = ETIME-STIME
1001    FORMAT(//,1X,'Time ','used ','=',F9.3,1X,'Sec',//,66('*'),
     $    /)
        WRITE(*,1001)CPTIME
        N = N+1
      GOTO 10
      ENDIF
      END
C
      SUBROUTINE FCN(N,X,F,IFAIL)
      INTEGER N
      DOUBLE PRECISION X(N),F(N)
C:    End Parameter
      INTEGER I,I1,L
      DOUBLE PRECISION TI2,TI1,TI,FACTT
C:    Begin
      DO 73 I=2,N,2
        I1 = I-1
        F(I1)=0.0
        F(I)=DBLE(N)/DBLE(I*I-1)
73    CONTINUE
      IF(MOD(N,2).EQ.1)THEN
        F(N)=0.0
      ENDIF
      DO 74 L=1,N
        FACTT = 4.0*X(L)-2.0
        TI2 = 1.0
        TI1 = 0.5*FACTT
        F(1)=TI1+F(1)
        DO 75 I=2,N
          TI = FACTT*TI1-TI2
          F(I)=TI+F(I)
          TI2 = TI1
          TI1 = TI
75      CONTINUE
74    CONTINUE
      RETURN
      END
