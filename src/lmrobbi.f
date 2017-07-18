C=======================================================================
C     ROBUST BOUNDED INFLUENCE ESTIMATION   { "Robust Library" =: "RL" }
C            =------ =--------                   names  "RL...BI" -> re-capitalized
C     MATHSOFT, INC.
C     10/25/99
C=======================================================================

      FUNCTION rlUPCVbi(S,IUCV,A,B)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA GAM,DSPI/1.E-6,2.506628274631001D0/
C-----------------------------------------------------------------------
C     COMPUTES THE U'-FUNCTION FOR AFFINE INVARIANT COVARIANCES
C-----------------------------------------------------------------------
C     IUCV: I, CLASSICAL IF IUCV=0,
C              HUBER     IF IUCV=1,
C              H-K       IF IUCV=2,
C              K-W       IF IUCV=3,
C              MALLOWRLU IF IUCV=4.
C-----------------------------------------------------------------------
      rlUPCVbi=0.D0
      IF (IUCV.EQ.0) RETURN
      ZED=S
      IF (IUCV.EQ.1) GOTO 100
      IF (IUCV.EQ.2) GOTO 200
      IF (IUCV.EQ.3) GOTO 300
      IF (IUCV.EQ.4) GOTO 400
 100  A2 = A*A
      B2 = B*B
      IF (S*S .GE. A2 .AND. S .GE. 0.D0) GOTO 110
      IF (S .GT. GAM) GOTO 110
      ZED=GAM
 110  Z2=ZED*ZED
      IF (Z2 .GT. B2) rlUPCVbi=-2.D0*B2/Z2/ZED
      IF (Z2 .LT. A2) rlUPCVbi=-2.D0*A2/Z2/ZED
      RETURN
 200  IF (S .LE. 0.D0) RETURN
      IF (S .LE. GAM) ZED=GAM
      Z2=ZED*ZED
      Q=A/ZED
      Q2=Q*Q
      PD=DEXP(-Q2/2.D0)/DSPI
      rlUPCVbi=2.D0*PD*(-A/Z2)
      RETURN
 300  IF (S .LE. 0.D0) RETURN
      IF (S .LE. GAM) ZED=GAM
      Q=A/ZED
      CALL rlGAUSbi(Q,PC)
      rlUPCVbi=-4.D0*(Q*Q)*(1.D0-PC)/ZED
      RETURN
 400  IF (S .LT. A) RETURN
      IF (S .GT. GAM) GOTO 410
      ZED=GAM
 410  IF (S .GT. A) rlUPCVbi=-A/(ZED*ZED)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlBET0bi(WGT,N,ITYPE,ISQW,TOL,BT0)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N)
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     COMPUTE THE CONSTANT BETA0 FOR NORMAL ERRORS
C-----------------------------------------------------------------------
      P=0.75D0
      CALL rlQUNTbi(P,BT0)
      IF (ITYPE .NE. 2) RETURN
      XN=DBLE(N)
      IF (ISQW .EQ. 0) GOTO 10
      E=2.D0
      IF (ISQW .EQ. 1) E=0.5D0
      DO 5 I=1,N
         IF (WGT(I) .LE. ZERO) GOTO 5
         WGT(I)=WGT(I)**E
 5    CONTINUE
 10   BT0=ZERO
 20   SMF=ZERO
      SMFP=ZERO
      DO 30 I=1,N
         IF (WGT(I) .LE. ZERO) GOTO 30
         X=BT0/WGT(I)
         CALL rlGAUSbi(X,A)
         CALL rlXERFbi(2,X,B)
         SMF=SMF+A
         SMFP=SMFP+B/WGT(I)
 30   CONTINUE
      F=SMF/XN-P
      FP=SMFP/XN
      BT0=BT0-F/FP
      IF (DABS(F) .LT. TOL) GOTO 40
      GOTO 20
 40   IF (ISQW .EQ. 0) RETURN
      E=1.D0/E
      DO 45 I=1,N
         IF (WGT(I) .LE. ZERO) GOTO 45
         WGT(I)=WGT(I)**E
 45   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlNLGMbi(N,GL)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA PI/3.1415926535898D0/
C-----------------------------------------------------------------------
C     LOGARITHM OF THE GAMMA-FUNCTION AT N/2 FOR INTEGER N
C-----------------------------------------------------------------------
      GL2=DLOG(2.D0)
      GL=0.D0
      K=N-2
 20   IF (K.LE.1) GOTO 30
      GL=GL+DLOG(DBLE(K))-GL2
      K=K-2
      GOTO 20
 30   IF (K.EQ.1) GL=GL+DLOG(DSQRT(PI))-GL2
      IF (N.EQ.1) GL=DLOG(DSQRT(PI))
      RETURN
      END
C=======================================================================
      SUBROUTINE rlXERPbi(IP,XLCNST,S,F)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA CMIN/-0.2257913526D0/
C-----------------------------------------------------------------------
      S2=-S*S/2.D0
      PP=DBLE(IP)
      IF (XLCNST.GT.CMIN .OR. XLCNST.EQ.0.D0) GOTO 30
      CALL rlNLGMbi(IP,XLGM)
      XLCNST=(1.D0-PP/2.D0)*DLOG(2.D0)-XLGM
 30   F=0.D0
      IF (S .LE. 0.D0) RETURN
      XLP=(PP-1.D0)*DLOG(S)+S2+XLCNST
      F=DEXP(XLP)
      RETURN
      END
C=======================================================================
      FUNCTION rl2PHIbi(S,WGT,N,SIGM,INDEX,TUNINGC,XLCNST,FPSI)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N),INDEX(7),TUNINGC(9)
      EXTERNAL FPSI
C
C avoid compiler warnings
C

      dummy = WGT(1)
      dummy = SIGM
      dummy = XLCNST

C-----------------------------------------------------------------------
C compute  Psi^2 (s)  or  Psi^2 (s / wgt[i]) -- psi(.) specified by FPSI(.)
c     used as integrand passed to rlIGRDbi(.)
C-----------------------------------------------------------------------
      R=S
      I=INDEX(6)
      CALL rlXERFbi(2,R,PHI)
      IF (INDEX(5) .EQ. 3) R=R/WGT(I)
      rl2PHIbi=FPSI(R,INDEX(4),TUNINGC(5))**2.D0*PHI
      RETURN
      END
C=======================================================================
      SUBROUTINE rlGAUSbi(X,P)
C.......................................................................
      DOUBLE PRECISION X, ROOT2, ROBLIBERF, ROBLIBERFC,P
      DATA ROOT2 /1.4142135623730950488D0/
C-----------------------------------------------------------------------
C     COMPUTES THE GAUSSIAN DISTRIBUTION FUNCTION
C-----------------------------------------------------------------------
C     X: I, QUANTILE
C     P: O, PROBABILITY.
C-----------------------------------------------------------------------
      IF (X .EQ. 0.D0) THEN
         P = 0.5D0
      ELSEIF(X .GT. 0.D0) THEN
         P = (1.D0+ROBLIBERF(X/ROOT2))/2.D0
      ELSE
         P = ROBLIBERFC(-X/ROOT2)/2.D0
      ENDIF
      RETURN
      END
C=======================================================================
      FUNCTION rlUCVbi(S,IUCV,A,B)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA GAM,DSPI/1.E-6,2.506628274631001D0/
C-----------------------------------------------------------------------
C     COMPUTES THE U-FUNCTION FOR AFFINE INVARIANT COVARIANCES
C-----------------------------------------------------------------------
C     IUCV=0: CLASSICAL
C     IUCV=1: HUBER MINMAX
C     IUCV=2: HAMPEL-KRASKER
C     IUCV=3: KRASKER-WELCH
C     IUCV=4: MALLOW UNSTANDARDIZED
C     IUCV=5: USERFD AS ON PAGE 139 OF MARAZZI (1993) WITH A2=4.0
C-----------------------------------------------------------------------
      rlUCVbi=1.D0
      IF (IUCV.EQ.0) RETURN
      ZED=S
      IF (IUCV.EQ.1) GOTO 100
      IF (IUCV.EQ.2) GOTO 200
      IF (IUCV.EQ.3) GOTO 300
      IF (IUCV.EQ.4) GOTO 400
      IF (IUCV.EQ.5) GOTO 500
 100  A2 = A*A
      B2 = B*B
      IF (S*S .GE. A2) GOTO 110
      IF (S .GT. GAM) GOTO 110
      ZED=GAM
 110  Z2=ZED*ZED
      IF (Z2 .GT. B2) rlUCVbi=B2/Z2
      IF (Z2 .LT. A2) rlUCVbi=A2/Z2
      RETURN
 200  IF (S .LE. 0.D0) RETURN
      IF (S .LE. GAM) ZED=GAM
      Q=A/ZED
      CALL rlGAUSbi(Q,PC)
      rlUCVbi=2.D0*PC-1.D0
      RETURN
 300  IF (S .LE. 0.D0) RETURN
      IF (S .LE. GAM) ZED=GAM
      Q=A/ZED
      Q2=Q*Q
      CALL rlGAUSbi(Q,PC)
      PD=DEXP(-Q2/2.D0)/DSPI
      rlUCVbi=Q2+(1.D0-Q2)*(2.D0*PC-1.D0)-2.D0*Q*PD
      RETURN
 400  IF (S .LE. A) RETURN
      IF (S .GT. GAM) GOTO 410
      ZED=GAM
 410  rlUCVbi=A/ZED
      RETURN
 500  rlUCVbi=A*1.D12
      IF (S .GT. 1.E-6) rlUCVbi=A/(S*S)
      RETURN
      END
C=======================================================================
      FUNCTION rlWWWbi(Z,IWWW,IUCV,A2,B2)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA GAM/1.D-6/
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF THE W BAR-FUNCTION
C-----------------------------------------------------------------------
      rlWWWbi=1.D0
      IF (IWWW .EQ. 0) RETURN
      IF (IWWW .EQ. 1) GOTO 100
      IF (IWWW .EQ. 2) GOTO 200
      IF (IWWW .EQ. 3) GOTO 300
 100  IF (Z .GT. GAM) GOTO 110
      Z=GAM
 110  rlWWWbi=1.D0/Z
      RETURN
 200  rlWWWbi=rlUCVbi(Z,IUCV,A2,B2)
      RETURN
 300  rlWWWbi=DSQRT(rlUCVbi(Z,IUCV,A2,B2))
      RETURN
      END
C=======================================================================
      SUBROUTINE rlXERFbi(KODE,X,P)
C.......................................................................
      DOUBLE PRECISION X,X2,P,SPI
      DATA SPI/2.506628274631D0/
C-----------------------------------------------------------------------
C     GAUSSIAN DENSITY FUNCTION
C-----------------------------------------------------------------------
      X2=-X*X/2.D0
      P=DEXP(X2)
      IF (KODE .EQ. 2) P=P/SPI
      RETURN
      END
C=======================================================================
      SUBROUTINE rlEPSHbi(C,EPSI2,EPSIP)
C.......................................................................
      DOUBLE PRECISION C,EPSI2,EPSIP,PC,PD,C2
C-----------------------------------------------------------------------
C     EXPECTED VALUE OF PSI(X,C)**2 AND PSI'(X,C)  for Huber Psi(.)
C-----------------------------------------------------------------------
c PC := pnorm(c)
      CALL rlGAUSbi(C,PC)
c PD := dnorm(c)
      CALL rlXERFbi(2,C,PD)
      C2=C*C
      EPSI2=C2+(1.D0-C2)*(2.D0*PC-1.D0)-2.D0*C*PD
      EPSIP=(2.D0*PC-1.D0)
      RETURN
      END
C=======================================================================
      FUNCTION rlPPHIbi(S,WGT,N,SIGM,INDEX,TUNINGC,XLCNST,FPSI)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N),INDEX(7),TUNINGC(9)
      EXTERNAL FPSI
C
C avoid compiler warnings
C

      dummy = WGT(1)
      dummy = SIGM
      dummy = XLCNST

C-----------------------------------------------------------------------
C compute  Psi(s)*phi(s)*s  or ... (s/wgt[i]) -- psi(.) specified by FPSI(.)
c     phi(.) = dnorm(.);  used to compute  integral[..  psi'(.) ]
c     used as integrand passed to rlIGRDbi(.)
C-----------------------------------------------------------------------
      R=S
      CALL rlXERFbi(2,R,PHI)
      PHI=R*PHI
      IF (INDEX(5) .EQ. 3) R=R/WGT(INDEX(6))
      rlPPHIbi=FPSI(R,INDEX(4),TUNINGC(5))*PHI
      RETURN
      END
C=======================================================================
      SUBROUTINE rlQK15bi(F,FARR,N,FEXT,A,RESULT,ABSERR,RESABS,
     +     RESASC,SIGM,INDEX,TUNINGC,XLCNST)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C     BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C     CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
C-----------------------------------------------------------------------
      DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8),FARR(N)
      DIMENSION INDEX(7),TUNINGC(9)
      EXTERNAL F,FEXT
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8)/
     *     9.914553711208126D-01,   9.491079123427585D-01,
     *     8.648644233597691D-01,   7.415311855993944D-01,
     *     5.860872354676911D-01,   4.058451513773972D-01,
     *     2.077849550078985D-01,   0.0D+00              /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8)/
     *     2.293532201052922D-02,   6.309209262997855D-02,
     *     1.047900103222502D-01,   1.406532597155259D-01,
     *     1.690047266392679D-01,   1.903505780647854D-01,
     *     2.044329400752989D-01,   2.094821410847278D-01/
      DATA WG(1),WG(2),WG(3),WG(4)/
     *     1.294849661688697D-01,   2.797053914892767D-01,
     *     3.818300505051189D-01,   4.179591836734694D-01/
C-----------------------------------------------------------------------
C     LIST OF MAJOR VARIABLES
C     -----------------------
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C     MACHINE DEPENDENT CONSTANTS
C     ---------------------------
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C-----------------------------------------------------------------------
      B = TUNINGC(8)
      CENTR = (A+B)/2.D0
      HLGTH = (B-A)/2.D0
      DHLGTH = DABS(HLGTH)
      EPMACH = TUNINGC(6)
      UFLOW = TUNINGC(7)
C-----------------------------------------------------------------------
C     COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C     THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C-----------------------------------------------------------------------
      FC = F(CENTR,FARR,N,SIGM,INDEX,TUNINGC,XLCNST,FEXT)
      RESG = FC*WG(4)
      RESK = FC*WGK(8)
      RESABS = DABS(RESK)
      DO 10 J=1,3
         JTW=J*2
         ABSC=HLGTH*XGK(JTW)
         FVAL1=F(CENTR-ABSC,FARR,N,SIGM,INDEX,TUNINGC,XLCNST,FEXT)
         FVAL2=F(CENTR+ABSC,FARR,N,SIGM,INDEX,TUNINGC,XLCNST,FEXT)
         FV1(JTW) = FVAL1
         FV2(JTW) = FVAL2
         FSUM = FVAL1+FVAL2
         RESG = RESG+WG(J)*FSUM
         RESK = RESK+WGK(JTW)*FSUM
         RESABS = RESABS+WGK(JTW)*(DABS(FVAL1)+DABS(FVAL2))
 10   CONTINUE
      DO 15 J = 1,4
         JTWM1 = J*2-1
         ABSC = HLGTH*XGK(JTWM1)
         FVAL1=F(CENTR-ABSC,FARR,N,SIGM,INDEX,TUNINGC,XLCNST,FEXT)
         FVAL2=F(CENTR+ABSC,FARR,N,SIGM,INDEX,TUNINGC,XLCNST,FEXT)
         FV1(JTWM1) = FVAL1
         FV2(JTWM1) = FVAL2
         FSUM = FVAL1+FVAL2
         RESK = RESK+WGK(JTWM1)*FSUM
         RESABS = RESABS+WGK(JTWM1)*(DABS(FVAL1)+DABS(FVAL2))
 15   CONTINUE
      RESKH = RESK/2.D0
      RESASC = WGK(8)*DABS(FC-RESKH)
      DO 20 J=1,7
         RESASC = RESASC+WGK(J)*(DABS(FV1(J)-RESKH)+
     +        DABS(FV2(J)-RESKH))
 20   CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = DABS((RESK-RESG)*HLGTH)
      IF(RESASC .NE. 0.D0 .AND. ABSERR .NE. 0.D0)
     +     ABSERR = RESASC*DMIN1(1.D0,((2.0D+02)*ABSERR/RESASC)**1.5D0)
      IF(RESABS .GT. UFLOW/((5.0D+01)*EPMACH))
     +     ABSERR = DMAX1((EPMACH*(5.0D+01))*RESABS,ABSERR)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlQSRTbi(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ELIST(LAST),IORD(LAST)
C-----------------------------------------------------------------------
C     CHECK WHETHER THE LIST CONTAINS MORE THAN TWO ERROR ESTIMATES.
C-----------------------------------------------------------------------
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
C-----------------------------------------------------------------------
C     THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
C     DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
C     ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
C     START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
C-----------------------------------------------------------------------
 10   ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
         ISUCC = IORD(NRMAX-1)
         IF(ERRMAX .LE. ELIST(ISUCC)) GO TO 30
         IORD(NRMAX) = ISUCC
         NRMAX = NRMAX-1
 20   CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE
C     MAINTAINED IN DESCENDING ORDER. THIS NUMBER
C     DEPENDS ON THE NUMBER OF SUBDIVISIONS STILL ALLOWED.
C-----------------------------------------------------------------------
 30   JUPBN = LAST
      IF(LAST .GT. (LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
C-----------------------------------------------------------------------
C     INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
C     STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
C-----------------------------------------------------------------------
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG .GT. JBND) GO TO 50
      DO 40 I=IBEG,JBND
         ISUCC = IORD(I)
         IF(ERRMAX .GE. ELIST(ISUCC)) GO TO 60
         IORD(I-1) = ISUCC
 40   CONTINUE
 50   IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
C-----------------------------------------------------------------------
C     INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
C-----------------------------------------------------------------------
 60   IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
         ISUCC = IORD(K)
         IF(ERRMIN .LT. ELIST(ISUCC)) GO TO 80
         IORD(K+1) = ISUCC
         K = K-1
 70   CONTINUE
      IORD(I) = LAST
      GO TO 90
 80   IORD(K+1) = LAST
C-----------------------------------------------------------------------
C     SET MAXERR AND ERMAX.
C-----------------------------------------------------------------------
 90   MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlIGRDbi(F,FARR,N,FEXT,A,EPSREL,LIMIT,RESULT,
     +     ABSERR,NEVAL,IER,SIGM,INDEX,TUNINGC,XLCNST,
     +     ALIST,BLIST,RLIST,ELIST,IORD)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IORD(LIMIT),INDEX(7)
      DIMENSION FARR(N),TUNINGC(9)
      DIMENSION ALIST(LIMIT),BLIST(LIMIT),RLIST(LIMIT),ELIST(LIMIT)
      EXTERNAL F,FEXT
c              ^  to be integrated:  Integral_A^B  F(x,...) dx
      DATA ZERO/0.D0/
c-----------------------------------------------------------------------
c Numerical Integration for  E[ .. ] computations ... unfortunately not further documented
C-----------------------------------------------------------------------
C     INDEX: (IPP,IWWW,IUCV,IPSI,ITYPE,I,IER1)
C     TUNINGC: (ZBAR2,BET2,A2,B2,C,EPMACH,UFLOW)
C-----------------------------------------------------------------------
C     LIST OF MAJOR VARIABLES
C     -----------------------
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                      (ALIST(I),BLIST(I))
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
C                       ERROR ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C
C     MACHINE DEPENDENT CONSTANTS
C     ---------------------------
C           EPMACH  IS THE LARGEST RELATIVE SPACING.
C           UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE.
C-----------------------------------------------------------------------
C     TEST ON VALIDITY OF PARAMETERS
C-----------------------------------------------------------------------
      NEVAL = 0
      LAST = 0
      B = TUNINGC(8)
      EPSABS = TUNINGC(9)
      RESULT = ZERO
      ABSERR = ZERO
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = ZERO
      ELIST(1) = ZERO
      IORD(1) = 0
      IER=0
      EPMACH = TUNINGC(6)
      UFLOW = TUNINGC(7)
C-----------------------------------------------------------------------
C     FIRST APPROXIMATION TO THE INTEGRAL
C-----------------------------------------------------------------------
      CALL rlQK15bi(F,FARR,N,FEXT,A,RESULT,ABSERR,DEFABS,RESABS,
     +     SIGM,INDEX,TUNINGC,XLCNST)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
C-----------------------------------------------------------------------
C     TEST ON ACCURACY
C-----------------------------------------------------------------------
      ERRBND = DMAX1(EPSABS,EPSREL*DABS(RESULT))
      IF(ABSERR.LE.(5.0D+01)*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0 .OR. (ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS)
     +     .OR.ABSERR.EQ.ZERO) GO TO 60
C-----------------------------------------------------------------------
C     INITIALIZATION
C-----------------------------------------------------------------------
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      NRMAX = 1
      IROFF1 = 0
      IROFF2 = 0
C-----------------------------------------------------------------------
C     MAIN DO-LOOP
C-----------------------------------------------------------------------
      DO 30 LAST = 2,LIMIT
C-----------------------------------------------------------------------
C     BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
C-----------------------------------------------------------------------
         A1 = ALIST(MAXERR)
         B1 = (5.0D-01)*(ALIST(MAXERR)+BLIST(MAXERR))
         A2 = B1
         B2 = BLIST(MAXERR)
         TUNINGC(8) = B1
         CALL rlQK15bi(F,FARR,N,FEXT,A1,AREA1,ERROR1,RESABS,
     +        DEFAB1,SIGM,INDEX,TUNINGC,XLCNST)
         TUNINGC(8) = B2
         CALL rlQK15bi(F,FARR,N,FEXT,A2,AREA2,ERROR2,RESABS,
     +        DEFAB2,SIGM,INDEX,TUNINGC,XLCNST)
C-----------------------------------------------------------------------
C     IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C     AND ERROR AND TEST FOR ACCURACY.
C-----------------------------------------------------------------------
         NEVAL = NEVAL+1
         AREA12 = AREA1+AREA2
         ERRO12 = ERROR1+ERROR2
         ERRSUM = ERRSUM+ERRO12-ERRMAX
         AREA = AREA+AREA12-RLIST(MAXERR)
         IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 5
         IF(DABS(RLIST(MAXERR)-AREA12).LE.(1.0D-05)*DABS(AREA12)
     +        .AND.ERRO12.GE.(9.9D-01)*ERRMAX) IROFF1 = IROFF1+1
         IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF2 = IROFF2+1
 5       RLIST(MAXERR) = AREA1
         RLIST(LAST) = AREA2
         ERRBND = DMAX1(EPSABS,EPSREL*DABS(AREA))
         IF(ERRSUM.LE.ERRBND) GO TO 8
C-----------------------------------------------------------------------
C     TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C-----------------------------------------------------------------------
         IF(IROFF1.GE.6.OR.IROFF2.GE.20) IER = 2
C-----------------------------------------------------------------------
C     SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
C     SUBINTERVALS EQUALS LIMIT.
C-----------------------------------------------------------------------
         IF(LAST.EQ.LIMIT) IER = 1
C-----------------------------------------------------------------------
C     SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C     AT A POINT OF THE INTEGRATION RANGE.
C-----------------------------------------------------------------------
         IF(DMAX1(DABS(A1),DABS(B2)).LE.((1.D0)+(1.0D+03)*
     +        EPMACH)*(DABS(A2)+(1.0D+04)*UFLOW)) IER = 3
C-----------------------------------------------------------------------
C     APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C-----------------------------------------------------------------------
 8       IF(ERROR2.GT.ERROR1) GO TO 10
         ALIST(LAST) = A2
         BLIST(MAXERR) = B1
         BLIST(LAST) = B2
         ELIST(MAXERR) = ERROR1
         ELIST(LAST) = ERROR2
         GO TO 20
 10      ALIST(MAXERR) = A2
         ALIST(LAST) = A1
         BLIST(LAST) = B1
         RLIST(MAXERR) = AREA2
         RLIST(LAST) = AREA1
         ELIST(MAXERR) = ERROR2
         ELIST(LAST) = ERROR1
C-----------------------------------------------------------------------
C     CALL SUBROUTINE rlQSRTbi TO MAINTAIN THE DESCENDING ORDERING
C     IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C     WITH THE LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C-----------------------------------------------------------------------
 20      CALL rlQSRTbi(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
         IF(IER.NE.0.OR.ERRSUM.LE.ERRBND) GO TO 40
 30   CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE FINAL RESULT.
C-----------------------------------------------------------------------
 40   TUNINGC(8) = B
      RESULT = 0.D0
      DO 50 K=1,LAST
         RESULT = RESULT+RLIST(K)
 50   CONTINUE
      ABSERR = ERRSUM
 60   NEVAL = 30*NEVAL+15
      RETURN
      END
C=======================================================================
      SUBROUTINE rlEPSUbi(EXPSI,ERREST,EPSI2,EPSIP,SIGM,INDEX,
     +     TUNINGC,XLCNST)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(1),INDEX(7),TUNINGC(9),IWRK(40)
      DIMENSION WRK1(40),WRK2(40),WRK3(40),WRK4(40)
      EXTERNAL rl2PHIbi,rlPPHIbi,EXPSI
C-----------------------------------------------------------------------
C     EXPECTED VALUES OF PSI(X)**2 AND PSP(X) = psi'(X)
C     	(WHERE X IS A STANDARD NORMAL DEVIATE)
c      for "non-Huber" psi functions, EXPSI(.) , further specified
c      via INDEX(.) & TUNINGC(.)
C-----------------------------------------------------------------------
      LIMIT=40
      CALL rlIGRDbi(rl2PHIbi,WGT,1,EXPSI,0.D0,TUNINGC(9),LIMIT,
     +     EPSI2,ERRST1,NEVAL1,IER1,SIGM,INDEX,TUNINGC,XLCNST,
     +     WRK1,WRK2,WRK3,WRK4,IWRK)
      EPSI2=2.D0*EPSI2
      CALL rlIGRDbi(rlPPHIbi,WGT,1,EXPSI,0.D0,TUNINGC(9),LIMIT,
     +     EPSIP,ERRST2,NEVAL2,IER2,SIGM,INDEX,TUNINGC,XLCNST,
     +     WRK1,WRK2,WRK3,WRK4,IWRK)
      EPSIP=2.D0*EPSIP
      NEVAL=NEVAL1+NEVAL2
      ERREST=DMAX1(ERRST1,ERRST2)
      IER=MAX0(IER1,IER2)
      RETURN
      END
C=======================================================================
      FUNCTION rlINS1bi(S,WGT,N,SIGM,INDEX,TUNINGC,XLCNST,EXPSI)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N),IWRK(20),WRK1(20),WRK2(20),WRK3(20),WRK4(20)
      DIMENSION INDEX(7),TUNINGC(9)
      EXTERNAL rlPPHIbi,EXPSI
C-----------------------------------------------------------------------
C     rlINS1bi(S)=E[ETA'(|ZI|)]*dG(S)
C     SERVES TO COMPUTE S1 WHEN OF THE FORM B1
C-----------------------------------------------------------------------
      IPP =INDEX(1)
      IPSI=INDEX(4)
      ITP =INDEX(5)
      I   =INDEX(6)
      IER1=INDEX(7)
      ZBAR2=TUNINGC(1)
      BET2 =TUNINGC(2)
      C    =TUNINGC(5)
C-----------------------------------------------------------------------
      ANS=1.D0
      Z=DSQRT(ZBAR2+BET2*S*S)
      WGT(I)=rlWWWbi(Z,INDEX(2),INDEX(3),TUNINGC(3),TUNINGC(4))
C-----------------------------------------------------------------------
C     WEIGHTS FOR SCHWEPPE ESTIMATORS WHEN FUNCTION PSI IS NOT
C     FROM HUBER TYP
C-----------------------------------------------------------------------
      IF (IPSI.EQ.3) GOTO 15
      LIMIT=20
      ERRST1=0.D0
      CALL rlIGRDbi(rlPPHIbi,WGT,N,EXPSI,0.D0,0.D0,LIMIT,RES,
     +     ERREST,NEVAL,IER,SIGM,INDEX,TUNINGC,XLCNST,WRK1,WRK2,
     +     WRK3,WRK4,IWRK)
      IER1=MAX0(IER1,IER)
      INDEX(7)=IER1
      ERRST1=DMAX1(ERREST,ERRST1)
      RES1=2.D0*RES*WGT(I)
      GOTO 10
C-----------------------------------------------------------------------
C     WEIGHTS FOR SCHWEPPE ESTIMATORS WHEN FUNCTION PSI IS
C     FROM HUBER TYP
C-----------------------------------------------------------------------
 15   C0=C*WGT(I)
      CALL rlEPSHbi(C0,EPSI2,RES1)
 10   IF(IPP.GT.0) THEN
         SBAR=S/SIGM
         CALL rlXERPbi(IPP,XLCNST,SBAR,ANS)
         ANS=ANS/SIGM
      ENDIF
C-----------------------------------------------------------------------
C     FOR MALLOWS AND SCHWEPPE CASES, E[ETA'(|ZI|)]dG(S)
C-----------------------------------------------------------------------
      IF (ITP.LT.3) rlINS1bi=WGT(I)*ANS
      IF (ITP.EQ.3) rlINS1bi=RES1*ANS
      WGT(I)=ZBAR2
      RETURN
      END
C=======================================================================
      FUNCTION rlINS2bi(S,WGT,N,SIGM,INDEX,TUNINGC,XLCNST,EXPSI)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N),IWRK(20),WRK1(20),WRK2(20),WRK3(20),WRK4(20)
      DIMENSION INDEX(7),TUNINGC(9)
      EXTERNAL rl2PHIbi,EXPSI
C-----------------------------------------------------------------------
C     rlINS2bi(S)=E[ETA**2(|ZI|)]*dG(S)
C     SERVES TO COMPUTE S2 WHEN OF THE FORM B1
C-----------------------------------------------------------------------
      IPP =INDEX(1)
      IPSI=INDEX(4)
      ITP =INDEX(5)
      I   =INDEX(6)
      IER1=INDEX(7)
      ZBAR2=TUNINGC(1)
      BET2 =TUNINGC(2)
      C    =TUNINGC(5)
C-----------------------------------------------------------------------
      ANS=1.D0
      Z=DSQRT(ZBAR2+BET2*S*S)
      WGT(I)=rlWWWbi(Z,INDEX(2),INDEX(3),TUNINGC(3),TUNINGC(4))
C-----------------------------------------------------------------------
C     WEIGHTS FOR SCHWEPPE ESTIMATORS
C     (FUNCTION PSI IS NOT FROM HUBER TYPE)
C-----------------------------------------------------------------------
      IF (IPSI.EQ.3) GOTO 15
      LIMIT=20
      ERRST1=0.D0
      CALL rlIGRDbi(rl2PHIbi,WGT,N,EXPSI,0.D0,0.D0,LIMIT,RES,
     +     ERREST,NEVAL,IER,SIGM,INDEX,TUNINGC,XLCNST,WRK1,WRK2,
     +     WRK3,WRK4,IWRK)
      IER1=MAX0(IER1,IER)
      INDEX(7)=IER1
      ERRST1=DMAX1(ERRST1,ERREST)
      RES1=2.D0*RES*WGT(I)*WGT(I)
      GOTO 10
C-----------------------------------------------------------------------
C     WEIGHTS FOR SCHWEPPE ESTIMATORS (FUNCTION PSI IS FROM HUBER TYPE)
C-----------------------------------------------------------------------
 15   C1=C*WGT(I)
      CALL rlEPSHbi(C1,RES1,EPSIP)
C-----------------------------------------------------------------------
C     MALLOWS et SCHWEPPE E[ETA**2]*dG(S)
C-----------------------------------------------------------------------
 10   IF (IPP.GT.0) THEN
         SBAR=S/SIGM
         CALL rlXERPbi(IPP,XLCNST,SBAR,ANS)
         ANS=ANS/SIGM
      ENDIF
      IF (ITP.LT.3) rlINS2bi=WGT(I)*WGT(I)*ANS
      IF (ITP.EQ.3) rlINS2bi=RES1*ANS
      WGT(I)=ZBAR2
      RETURN
      END
C=======================================================================
      FUNCTION rlINS3bi(S,WGT,N,SIGM,INDEX,TUNINGC,XLCNST,EXPSI)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N),INDEX(7),TUNINGC(9)
      EXTERNAL EXPSI,rlINS1bi
C-----------------------------------------------------------------------
C     rlINS3bi(S)=1/n*SUM{E[ETA'(|ZI|)]}*dG(S)
C     TO COMPUTE S1 WHEN OF THE FORM B2
C-----------------------------------------------------------------------
      SUM=0.D0
      DO 10 J=1,N
         I=J
         INDEX(6)=I
         TUNINGC(1)=WGT(I)
         SUM=SUM+rlINS1bi(S,WGT,N,SIGM,INDEX,TUNINGC,XLCNST,EXPSI)
 10   CONTINUE
      rlINS3bi=SUM*S*S/DBLE(N)
      RETURN
      END
C=======================================================================
      FUNCTION rlINS4bi(S,WGT,N,SIGM,INDEX,TUNINGC,XLCNST,EXPSI)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N),INDEX(7),TUNINGC(9)
      EXTERNAL EXPSI,rlINS2bi
C-----------------------------------------------------------------------
C     rlINS4bi(S)=1/n*SUM{E[ETA**2(|ZI|)]}*dG(S)
C     TO COMPUTE S2 WHEN OF THE FORM B2
C-----------------------------------------------------------------------
      SUM=0.D0
      DO 10 J=1,N
         I=J
         INDEX(6)=I
         TUNINGC(1)=WGT(I)
         SUM=SUM+rlINS2bi(S,WGT,N,SIGM,INDEX,TUNINGC,XLCNST,EXPSI)
 10   CONTINUE
      rlINS4bi=SUM*S*S/DBLE(N)
      RETURN
      END
C======================================================================
      FUNCTION rlUZEDbi(S,WGT,N,SIGM,INDEX,TUNINGC,XLCNST,EXU)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N),INDEX(7),TUNINGC(9)
      EXTERNAL EXU
C
C
C

      dummy = WGT(1)

C-----------------------------------------------------------------------
C     rlUZEDbi(S)=U(SQRT(ZBAR2+BET2*S^2))*dG(S)
C-----------------------------------------------------------------------
      IPP=INDEX(1)
      ZBAR2=TUNINGC(1)
      BET2=TUNINGC(2)
C-----------------------------------------------------------------------
      IF (IPP.GT.0) GOTO 15
      Z=DSQRT(ZBAR2)
      ANS=1.D0
      GOTO 20
 15   SBAR=S/SIGM
      CALL rlXERPbi(IPP,XLCNST,SBAR,ANS)
      ANS=ANS/SIGM
      Z=DSQRT(ZBAR2+BET2*S*S)
 20   rlUZEDbi=EXU(Z,INDEX(3),TUNINGC(3),TUNINGC(4))*ANS
      RETURN
      END
C=======================================================================
      FUNCTION rlUZD2bi(S,WGT,N,SIGM,INDEX,TUNINGC,XLCNST,EXU)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N),INDEX(7),TUNINGC(9)
      EXTERNAL EXU
C-----------------------------------------------------------------------
C     rlUZD2bi(S)=AVE{U(SQRT(ZBAR2+BET2*S**2))*S**2*dG(S)}
C-----------------------------------------------------------------------
      IPP=INDEX(1)
      BET2=TUNINGC(2)
C-----------------------------------------------------------------------
      U=0.D0
      DO 10 L=1,N
         ZBAR2=WGT(L)*WGT(L)
         Z=DSQRT(ZBAR2+BET2*S*S)
         U=U+EXU(Z,INDEX(3),TUNINGC(3),TUNINGC(4))
 10   CONTINUE
      TUNINGC(1)=ZBAR2
      SBAR=S/SIGM
      CALL rlXERPbi(IPP,XLCNST,SBAR,ANS)
      XN=DBLE(N)*SIGM
      rlUZD2bi=(U/XN)*S*S*ANS
      RETURN
      END
C=======================================================================
      SUBROUTINE rlREF0bi(INDEX,TUNINGC,XLCNST,IALFA,SIGM,MAXIT,
     +     TOL,NIT,ALFA,BETA,REFF)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(1),INDEX(7),TUNINGC(9)
      DIMENSION IWORK(20),WORK1(20),WORK2(20),WORK3(20),WORK4(20)
      EXTERNAL rlINS1bi,rlINS2bi,rlINS3bi,rlINS4bi,RLPSIM2,rlUZEDbi,
     +     rlUZD2bi,rlUCVbi
      DATA ZERO,ONE /0.D0,1.D0/
C-----------------------------------------------------------------------
C     INDEX: (IPP,IWWW,IUCV,IPSI,ITYPE,I,IER1)
C     TUNINGC: (ZBAR2,BET2,A2,B2,C,EPMACH,UFLOW,UPPER,TIL)
C-----------------------------------------------------------------------
C     ASYMPTOTIC RELATIVE EFFICIENCY OF A GM-ESTIMATOR
C-----------------------------------------------------------------------
      IPP  = INDEX(1)
      IPSI = INDEX(4)
      ITYP = INDEX(5)
      C    = TUNINGC(5)
C-----------------------------------------------------------------------
C     SPECIAL FOR  ITYPE.EQ.1 .OR. (MU.EQ.0 .AND. ITYPE.EQ.2):
C     Compute Epsi**2/Epsi'
C-----------------------------------------------------------------------
      IF (ITYP.EQ.3) GOTO 5
      IF (IPSI.EQ.3) THEN
         IF (C .LE. ZERO) C=1.345D0
         CALL rlEPSHbi(C,G1,G0)
      ELSE
         CALL rlEPSUbi(RLPSIM2,ERREST,G1,G0,SIGM,INDEX,TUNINGC,XLCNST)
      ENDIF
      IF (ITYP.NE.1 .AND. (IPP.GT.0 .OR. ITYP.NE.2)) GOTO 5
      ALFA=G1
      BETA=G0
      REFF=(G0**2)/G1
      RETURN
 5    CONTINUE
      P=DBLE(IPP)
      Q=DBLE(IALFA)
      LIMIT=20
C-----------------------------------------------------------------------
C     COVARIANCE LS-ESTIMATOR
C-----------------------------------------------------------------------
      TRCVLS=Q+P/(SIGM**2)
C-----------------------------------------------------------------------
C     COMPUTATION OF ALFA AND BETA (ITERATIVE ALGORITHM)
C-----------------------------------------------------------------------
      DS=ONE
      DALF=ZERO
      DBET=ZERO
C-----------------------------------------------------------------------
C     STEP 0: INITIALIZATION
C-----------------------------------------------------------------------
      NIT=1
      TUNINGC(1)=ONE
      TUNINGC(2)=ONE
      IF (IALFA.EQ.0) TUNINGC(1)=ZERO
      IF (IPP  .EQ.0) TUNINGC(2)=ZERO
 10   IF (IALFA.EQ.0) GOTO 20
C-----------------------------------------------------------------------
C     STEP 1: SOLVE FOR ALFA
C-----------------------------------------------------------------------
      IF (IPP .EQ. 0)
     +     ANS1=rlUZEDbi(DS,WGT,1,SIGM,INDEX,TUNINGC,XLCNST,rlUCVbi)
      IF (IPP.GT.0)
     +     CALL rlIGRDbi(rlUZEDbi,WGT,1,rlUCVbi,ZERO,ZERO,
     +     LIMIT,ANS1,ERRSTD,NEVAL,IER,SIGM,INDEX,TUNINGC,XLCNST,
     +     WORK1,WORK2,WORK3,WORK4,IWORK)
      TUNINGC(1)=ONE/ANS1
      ALF=DSQRT(TUNINGC(1))
      DALF=DABS(ALFA-ALF)
      ALFA=ALF
 20   IF (IPP .EQ. 0) GOTO 30
C-----------------------------------------------------------------------
C     STEP 2: SOLVE FOR BETA
C-----------------------------------------------------------------------
      WGT(1)=DSQRT(TUNINGC(1))
      CALL rlIGRDbi(rlUZD2bi,WGT,1,rlUCVbi,ZERO,ZERO,LIMIT,
     +     ANS2,ERRSTD,NEVAL,IER,SIGM,INDEX,TUNINGC,XLCNST,
     +     WORK1,WORK2,WORK3,WORK4,IWORK)
      TUNINGC(2)=P/ANS2
      BET=DSQRT(TUNINGC(2))
      DBET=DABS(BETA-BET)
      BETA=BET
C-----------------------------------------------------------------------
C     STEP 3: CHECK CONVERGENCE
C-----------------------------------------------------------------------
 30   IF ((DALF.LT.TOL.AND.DBET.LT.TOL).OR.(NIT.EQ.MAXIT)) GOTO 40
      NIT=NIT+1
      GOTO 10
C-----------------------------------------------------------------------
C     COVARIANCE M-ESTIMATORS
C-----------------------------------------------------------------------
 40   TRCOV=ZERO
      IF (IALFA.EQ.0 .AND. IPP.GT.0) GOTO 50
      IF (IALFA.EQ.1 .AND. IPP.EQ.0) GOTO 60
C-----------------------------------------------------------------------
C  STEP 4: COMPUTE S1 AND S2 FOR THE QUALITATIVE COVARIATE
C-----------------------------------------------------------------------
      CALL rlIGRDbi(rlINS2bi,WGT,1,RLPSIM2,ZERO,ZERO,LIMIT,
     +     ANS2,ERRSTD,NEVAL,IER,SIGM,INDEX,TUNINGC,XLCNST,
     +     WORK1,WORK2,WORK3,WORK4,IWORK)
      INDEX(7)=0
      CALL rlIGRDbi(rlINS1bi,WGT,1,RLPSIM2,ZERO,ZERO,LIMIT,
     +     ANS1,ERRSTD,NEVAL,IER,SIGM,INDEX,TUNINGC,XLCNST,
     +     WORK1,WORK2,WORK3,WORK4,IWORK)
      INDEX(7)=0
      TRCOV=ANS2/ANS1**2
C-----------------------------------------------------------------------
C     COMPUTE S1 AND S2 FOR THE QUANTITATIVE COVARIATES
C-----------------------------------------------------------------------
 50   CALL rlIGRDbi(rlINS3bi,WGT,1,RLPSIM2,ZERO,ZERO,LIMIT,
     +     ANS3,ERRSTD,NEVAL,IER,SIGM,INDEX,TUNINGC,XLCNST,
     +     WORK1,WORK2,WORK3,WORK4,IWORK)
      INDEX(7)=0
      CALL rlIGRDbi(rlINS4bi,WGT,1,RLPSIM2,ZERO,ZERO,LIMIT,
     +     ANS4,ERRSTD,NEVAL,IER,SIGM,INDEX,TUNINGC,XLCNST,
     +     WORK1,WORK2,WORK3,WORK4,IWORK)
      INDEX(7)=0
C-----------------------------------------------------------------------
C     STEP 5: COMPUTE ARE1
C-----------------------------------------------------------------------
      FONCT=TRCOV+P*P*(ANS4/ANS3**2.D0)
      REFF=(TRCVLS)/FONCT
      GOTO 70
 60   ANS2=rlINS2bi(DS,WGT,1,SIGM,INDEX,TUNINGC,XLCNST,RLPSIM2)
      ANS1=rlINS1bi(DS,WGT,1,SIGM,INDEX,TUNINGC,XLCNST,RLPSIM2)
      FONCT=ANS2/ANS1**2
      REFF=ONE/FONCT
 70   IF (ITYP.EQ.3) RETURN
      REFF=REFF*(G0**2)/G1
      RETURN
      END
C=======================================================================
      SUBROUTINE rlCOLbi(V1,V2,MLT,M,IOUT)
C.......................................................................
      DOUBLE PRECISION V1(M),V2(M),MLT
C-----------------------------------------------------------------------
C     AUXILIARY ROUTINE FOR rlLARSbi
C-----------------------------------------------------------------------
      DO 220 I=1,M
         IF (I .EQ. IOUT) GOTO 220
         V1(I)=V1(I)-V2(I)*MLT
 220  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlICHGbi(A,B)
C.......................................................................
C     AUXILIARY ROUTINE FOR rlLARSbi
C-----------------------------------------------------------------------
      DOUBLE PRECISION A,B,C
      C=A
      A=B
      B=C
      RETURN
      END
C=======================================================================
      SUBROUTINE rlLARSbi(X,Y,N,NP,MDX,MDT,TOL,NIT,K,
     +     KODE,SIGMA,THETA,RS,SC1,SC2,SC3,SC4,BET0)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(MDT),RS(N),SC1(N),SC2(NP),
     +     SC3(NP),SC4(NP)
      INTEGER OUT
      LOGICAL STAGE,TEST
      DATA ZERO,TWO,EPS,BIG/0.D0,2.D0,1.0D-10,3.401D38/
C      DATA ZERO,TWO,EPS,BIG/0.D0,2.D0,2.22D-16,1.796D308/
C-----------------------------------------------------------------------
C     LEAST ABSOLUTE RESIDUALS -- aka  L_1 - Regression
C      --> Result in THETA[1:NP]
C-----------------------------------------------------------------------
      SUM=ZERO
      DO 10 J=1,NP
         SC4(J)=DBLE(J)
         SC2(J)=ZERO
 10   CONTINUE
      DO 40 I=1,N
         SC1(I)=DBLE(NP+I)
         THETA(I)=Y(I)
         IF (Y(I) .GE. ZERO) GOTO 30
         DO 20 J=1,NP
            X(I,J)=-X(I,J)
 20      CONTINUE
         THETA(I)=-THETA(I)
         SC1(I)=-SC1(I)
 30      SUM=SUM+THETA(I)
 40   CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE THE MARGINAL COSTS.
C-----------------------------------------------------------------------
      SUMIN=SUM
      DO 60 J=1,NP
         SUM=ZERO
         DO 50 I=1,N
            SUM=SUM+X(I,J)
 50      CONTINUE
         SC3(J)=SUM
 60   CONTINUE
C-----------------------------------------------------------------------
C     STAGE I. DETERMINE THE VECTOR TO ENTER THE BASIS.
C-----------------------------------------------------------------------
      STAGE=.TRUE.
      KOUNT=0
      KR=1
      KL=1
 70   VMAX=-1.D0
      DNP=DBLE(NP)
      DO 80 J=KR,NP
         IF (DABS(SC4(J)) .GT. DNP) GOTO 80
         D=DABS(SC3(J))
         IF (D-VMAX .LE. ZERO) GOTO 80
         IF (D-VMAX .LE. EPS)  GOTO 80
         VMAX=D
         IN=J
 80   CONTINUE
      IF (SC3(IN) .GE. ZERO) GOTO 100
      DO 90 I=1,N
         X(I,IN)=-X(I,IN)
 90   CONTINUE
      SC3(IN)=-SC3(IN)
      SC4(IN)=-SC4(IN)
C-----------------------------------------------------------------------
C     DETERMINE THE VECTOR TO LEAVE THE BASIS.
C-----------------------------------------------------------------------
 100  K=0
      DO 110 I=KL,N
         D=X(I,IN)
         IF (D .LE. TOL) GOTO 110
         K=K+1
         Y(K)=THETA(I)/D
         RS(K)=DBLE(I)
         TEST=.TRUE.
 110  CONTINUE
 120  IF (K .GT. 0) GOTO 130
      TEST=.FALSE.
      GOTO 150
 130  VMIN=BIG
      DO 140 I=1,K
         IF (Y(I)-VMIN .GE. ZERO) GOTO 140
         IF (VMIN-Y(I) .LE. EPS)  GOTO 140
         J=I
         VMIN=Y(I)
         OUT=INT(RS(I))
 140  CONTINUE
      Y(J)=Y(K)
      RS(J)=RS(K)
      K=K-1
C-----------------------------------------------------------------------
C     CHECK FOR LINEAR DEPENDENCE IN STAGE I.
C-----------------------------------------------------------------------
 150  IF (TEST .OR. .NOT.STAGE) GOTO 170
      DO 160 I=1,N
         CALL rlICHGbi(X(I,KR),X(I,IN))
 160  CONTINUE
      CALL rlICHGbi(SC3(KR),SC3(IN))
      CALL rlICHGbi(SC4(KR),SC4(IN))
      KR=KR+1
      GOTO 260
 170  IF (TEST) GOTO 180
      KODE=2
      GOTO 350
 180  PIVOT=X(OUT,IN)
      IF (SC3(IN)-PIVOT-PIVOT .LE. TOL) GOTO 200
      DO 190 J=KR,NP
         D=X(OUT,J)
         SC3(J)=SC3(J)-D-D
         X(OUT,J)=-D
 190  CONTINUE
      D=THETA(OUT)
      SUMIN=SUMIN-D-D
      THETA(OUT)=-D
      SC1(OUT)=-SC1(OUT)
      GOTO 120
C-----------------------------------------------------------------------
C     PIVOT ON X(OUT,IN).
C-----------------------------------------------------------------------
 200  DO 210 J=KR,NP
         IF (J.EQ.IN) GOTO 210
         X(OUT,J)=X(OUT,J)/PIVOT
 210  CONTINUE
      THETA(OUT)=THETA(OUT)/PIVOT
      DO 230 J=KR,NP
         IF (J .EQ. IN) GOTO 230
         D=X(OUT,J)
         SC3(J)=SC3(J)-D*SC3(IN)
         CALL rlCOLbi(X(1,J),X(1,IN),D,N,OUT)
 230  CONTINUE
      SUMIN=SUMIN-SC3(IN)*THETA(OUT)
      DO 240 I=1,N
         IF (I .EQ. OUT) GOTO 240
         D=X(I,IN)
         THETA(I)=THETA(I)-D*THETA(OUT)
         X(I,IN)=-D/PIVOT
 240  CONTINUE
      SC3(IN)=-SC3(IN)/PIVOT
      X(OUT,IN)=1.D0/PIVOT
      CALL rlICHGbi(SC1(OUT),SC4(IN))
      KOUNT=KOUNT+1
      IF (.NOT.STAGE) GOTO 270
C-----------------------------------------------------------------------
C     INTERCHANGE ROWS IN STAGE I.
C-----------------------------------------------------------------------
      KL=KL+1
      DO 250 J=KR,NP
         CALL rlICHGbi(X(OUT,J),X(KOUNT,J))
 250  CONTINUE
      CALL rlICHGbi(THETA(OUT),THETA(KOUNT))
      CALL rlICHGbi(SC1(OUT),SC1(KOUNT))
 260  IF (KOUNT+KR .NE. NP+1) GOTO 70
C-----------------------------------------------------------------------
C     STAGE II. DETERMINE THE VECTOR TO ENTER THE BASIS.
C-----------------------------------------------------------------------
      STAGE=.FALSE.
 270  VMAX=-BIG
      DO 290 J=KR,NP
         D=SC3(J)
         IF (D .GE. ZERO) GOTO 280
         IF (D+TWO .GT. ZERO) GOTO 290
         D=-D-TWO
 280     IF (D-VMAX .LE. ZERO) GOTO 290
         IF (D-VMAX .LE. EPS)  GOTO 290
         VMAX=D
         IN=J
 290  CONTINUE
      IF (VMAX .LE. TOL) GOTO 310
      IF (SC3(IN) .GT. ZERO) GOTO 100
      DO 300 I=1,N
         X(I,IN)=-X(I,IN)
 300  CONTINUE
      SC3(IN)=-SC3(IN)-2.D0
      SC4(IN)=-SC4(IN)
      GOTO 100
C-----------------------------------------------------------------------
C     PREPARE OUTPUT
C-----------------------------------------------------------------------
 310  L=KL-1
      DO 330 I=1,N
         RS(I)=ZERO
         IF (I .GT. L .OR. THETA(I) .GE. ZERO) GOTO 330
         DO 320 J=KR,NP
            X(I,J)=-X(I,J)
 320     CONTINUE
         THETA(I)=-THETA(I)
         SC1(I)=-SC1(I)
 330  CONTINUE
      KODE=0
      IF (KR .NE. 1) GOTO 350
      DO 340 J=1,NP
         D=DABS(SC3(J))
         IF (D .LE. TOL .OR. TWO-D .LE. TOL) GOTO 350
 340  CONTINUE
      KODE=1
 350  DO 380 I=1,N
         K=INT(SC1(I))
         D=THETA(I)
         IF (K .GT. 0) GOTO 360
         K=-K
         D=-D
 360     IF (I .GE. KL) GOTO 370
         SC2(K)=D
         GOTO 380
 370     K=K-NP
         RS(K)=D
 380  CONTINUE
      K=NP+1-KR
      SUM=ZERO
      DO 390 I=KL,N
         SUM=SUM+THETA(I)
 390  CONTINUE
      SUMIN=SUM
      NIT=KOUNT
      DO 400 J=1,NP
         THETA(J)=SC2(J)
 400  CONTINUE
      DO 500 I=1,N
         Y(I)=DABS(RS(I))
 500  CONTINUE
      N2=N/2+1
      CALL RLSTORM2(Y,N,N2,SIGMA)
      SIGMA=SIGMA/BET0
      RETURN
      END
C=======================================================================
      SUBROUTINE rlSRT1bi(A,N,K1,K2)
C.......................................................................
      DOUBLE PRECISION A(N),X
C----------------------------------------------------------------------
C     SORTS THE VECTOR A[K1:K2] IN ASCENDING ORDER
C----------------------------------------------------------------------
      N1=K2-K1+1
      I=1
 10   I=I+I
      IF (I .LE. N1) GOTO 10
      M=I-1
 20   M=M/2
      IF (M .EQ. 0) GOTO 90
      K=N1-M
      DO 40 J=1,K
         L=J
 50      IF (L .LT. 1) GOTO 40
         LPM=L+M
         LPM1=LPM+K1-1
         L1=L+K1-1
         IF (A(LPM1) .GE. A(L1)) GOTO 40
         X=A(LPM1)
         A(LPM1)=A(L1)
         A(L1)=X
         L=L-M
         GOTO 50
 40   CONTINUE
      GOTO 20
 90   CONTINUE
      END
C=======================================================================
      SUBROUTINE rlKEDHbi(WGT,N,C,ITYPE,D,E)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N),E(N),D(N)
C-----------------------------------------------------------------------
C     COMPUTES THE DIAGONAL MATRIX FOR HUBER FUNCTION
C-----------------------------------------------------------------------
      IF (ITYPE.EQ.3) GOTO 30
C-----------------------------------------------------------------------
C     MALLOWS CASE
C-----------------------------------------------------------------------
      C2=C*C
      CALL rlGAUSbi(C,PC)
      CALL rlXERFbi(2,C,PD)
      G1=C2+(1.D0-C2)*(2.D0*PC-1.D0)-2.D0*C*PD
      F1=2.D0*PC-1.D0
      DO 20 I=1,N
         D(I)=F1*WGT(I)
         E(I)=G1*WGT(I)*WGT(I)
 20   CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     SCHWEPPE CASE
C-----------------------------------------------------------------------
 30   DO 40 I=1,N
         Z=C*WGT(I)
         Z2=Z*Z
         CALL rlGAUSbi(Z,PC)
         CALL rlXERFbi(2,Z,PD)
         E(I)=Z2+(1.D0-Z2)*(2.D0*PC-1.D0)-2.D0*Z*PD
         D(I)=2.D0*PC-1.D0
 40   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlKTASbi(X,D,E,N,NP,MDX,MDSC,NCOV,TAU,IA,F,F1,
     +     IAINV,A,S1INV,S2,AINV,COV,SC)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),D(N),E(N),S1INV(NCOV),S2(NCOV),
     +     A(NCOV),AINV(NCOV),COV(NCOV),SC(MDSC,NP)
      DATA DZERO/0.D0/
C-----------------------------------------------------------------------
C     COMPUTES THE COVARIANCE MATRIX OF COEFFICIENT ESTIMATES
C-----------------------------------------------------------------------
      NN=NP*(NP+1)/2
      XN1=DBLE(N)
C-----------------------------------------------------------------------
C     IF IA.EQ.-1 SET S1INV=F1*A
C-----------------------------------------------------------------------
      IF (IA .NE. -1) GOTO 40
      DO 35 L=1,NCOV
         S1INV(L)=A(L)
 35   CONTINUE
      IF (F1 .GT. DZERO) CALL RLSCALM2(S1INV,F1,NN,1,NN)
 40   CONTINUE
C-----------------------------------------------------------------------
C     IF IA.EQ.0 SET S1INV=F1*(A**T)*A
C-----------------------------------------------------------------------
      IF (IA .NE. 0) GOTO 45
      CALL RLMTT1M2(A,S1INV,NP,NN)
      IF (F1 .GT. DZERO) CALL RLSCALM2(S1INV,F1,NN,1,NN)
 45   CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE S2=X**T*E*X/N (AND STORE IT IN S2).
C     IF IA.EQ.1 COMPUTE S1=X**T*D*X/N
C     AND STORE IT TEMPORARILY IN COV)
C-----------------------------------------------------------------------
      L=0
      DO 62 I=1,NP
         DO 60 J=1,I
            L=L+1
            SM2=DZERO
            SM1=DZERO
            DO 50 K=1,N
               DXX=X(K,I)*X(K,J)
               SM2=SM2+DXX*E(K)
               IF (IA .EQ. 1) SM1=SM1+DXX*D(K)
 50         CONTINUE
            S2(L)=SM2/XN1
            IF (IA.EQ.1) COV(L)=SM1/XN1
 60      CONTINUE
 62   CONTINUE
C-----------------------------------------------------------------------
C     IF IA .EQ.1 COMPUTE A LOWER TRIANGULAR MATRIX A (AND ITS INVERSE)
C     SUCH THAT S1INV=S1**(-1)=A**T*A (IF IAINV.EQ.1 STORE THE INVERSE
C     OF A IN AINV)
C-----------------------------------------------------------------------
      IF (IA .EQ. -1 .OR. IA .EQ. 0) GOTO 80
      CALL RLMCHLM2(COV,NP,NN,INFO)
      IF (INFO .EQ. 0) GOTO 65
      RETURN
 65   CONTINUE
      DO 70 L=1,NN
         IF (IAINV .EQ. 1) AINV(L)=COV(L)
         A(L)=COV(L)
 70   CONTINUE
      CALL RLMINVM2(A,NP,NN,TAU,ISING)
      IF (ISING.EQ.0) GOTO 75
      RETURN
 75   CONTINUE
      CALL RLMTT1M2(A,S1INV,NP,NN)
 80   CALL rlMSSDbi(S2,S1INV,SC,NP,NN,MDSC)
C-----------------------------------------------------------------------
C     COMPUTE COV=F*S1**(-1)*S2*S1**(-1)
C-----------------------------------------------------------------------
      CALL rlMSF1bi(S1INV,SC,COV,NP,NN,MDSC)
      IF (F .GT. DZERO) CALL RLSCALM2(COV,F,NN,1,NN)
C-----------------------------------------------------------------------
C     IF IAINV.EQ.1 (AND IA.NE.1) COMPUTE THE INVERSE
C     OF A AND STORE IT IN AINV
C-----------------------------------------------------------------------
      IF (IA .EQ. 1 .OR. IAINV .EQ. 0) RETURN
      DO 90 L=1,NN
         AINV(L)=A(L)
 90   CONTINUE
      CALL RLMINVM2(AINV,NP,NN,TAU,ISING)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlMLYDbi(A,Y,N,NN,NY,IYE)
C.......................................................................
      DOUBLE PRECISION A(NN),Y(NY),SM,DZERO
      DATA DZERO/0.D0/
C-----------------------------------------------------------------------
C     MULTIPLIES A LOWER TRIANGULAR MATRIX A BY A VECTOR Y
C-----------------------------------------------------------------------
C     A: I, LOWER TRIANGULAR MATRIX WITH DIMENSION NN=N*(N+1)/2
C     Y: I/O, VECTOR Y ON INPUT, AND VECTOR A*Y ON OUTPUT.
C-----------------------------------------------------------------------
      IA=NN
      IY1=N*IYE+1
      DO 20 J1=1,N
         J=N-J1+1
         IY1=IY1-IYE
         IY=IY1
         SM=DZERO
         DO 10 I=1,J
            SM=SM+A(IA)*Y(IY)
            IA=IA-1
            IY=IY-IYE
 10      CONTINUE
         Y(IY1)=SM
 20   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlNRM2bi(X,N,INCX,MDX,XNRM)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX)
      DATA ZERO,ONE,CUTLO,CUTHI/0.D0,1.D0,4.441D-16,1.304D19/
C-----------------------------------------------------------------------
C     FORMS THE EUCLIDEAN NORM OF A VECTOR
C-----------------------------------------------------------------------
      IF (N .GT. 0) GOTO 10
      XNRM=ZERO
      GOTO 300
c10   ASSIGN 30 TO NEXT
 10   NEXT=30
      SUM=ZERO
      NN=N*INCX
C-----------------------------------------------------------------------
C     BEGIN MAIN LOOP
C-----------------------------------------------------------------------
      I=1
 20   DXI=X(I)
c     GOTO NEXT,(30,50,70,110)
      if(next.eq.30) goto 30
      if(next.eq.50) goto 50
      if(next.eq.70) goto 70
      if(next.eq.110) goto 110
 30   IF (DABS(X(I)) .GT. CUTLO) GOTO 85
c     ASSIGN 50 TO NEXT
      NEXT=50
      XMAX=ZERO
C-----------------------------------------------------------------------
C     PHASE1.  SUM IS ZERO
C-----------------------------------------------------------------------
 50   IF (X(I) .EQ. ZERO) GOTO 200
      IF (DABS(X(I)) .GT. CUTLO) GOTO 85
c     ASSIGN 70 TO NEXT
      next=70
      GOTO 105
 100  I=J
      DXI=X(I)
c     ASSIGN 110 TO NEXT
      next=110
      SUM=(SUM/DXI)/DXI
 105  XMAX=DABS(DXI)
      GOTO 115
C-----------------------------------------------------------------------
C     PHASE 2.  SUM IS SMALL. SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C-----------------------------------------------------------------------
 70   IF (DABS(X(I)) .GT. CUTLO) GOTO 75
C-----------------------------------------------------------------------
C     COMMON CODE FOR PHASE 2 AND 4.
C     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C-----------------------------------------------------------------------
 110  IF (DABS(DXI) .LE. XMAX) GOTO 115
      SUM=ONE+SUM*(XMAX/DXI)**2
      XMAX=DABS(DXI)
      GOTO 200
 115  SUM=SUM+(DXI/XMAX)**2
      GOTO 200
 75   SUM=(SUM*XMAX)*XMAX
 85   HITEST=CUTHI/DBLE(N)
C-----------------------------------------------------------------------
C     PHASE3.  SUM IS MID-RANGE.  NO SCALING.
C-----------------------------------------------------------------------
      DO 95 J=I,NN,INCX
         IF (DABS(X(J)) .GE. HITEST) GOTO 100
         SUM=SUM+X(J)*X(J)
 95   CONTINUE
      XNRM=DSQRT(SUM)
      GOTO 300
 200  CONTINUE
      I=I+INCX
      IF (I .LE. NN) GOTO 20
C-----------------------------------------------------------------------
C     END MAIN LOOP
C-----------------------------------------------------------------------
      XNRM=XMAX*DSQRT(SUM)
 300  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlWEDVbi(X,NVAR,NCOV,MDX,ITYPW,INIT,NFIRST,A,SC)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NCOV),X(MDX,NVAR),SC(NFIRST)
      DATA ZERO,ONE,TL/0.D0,1.D0,1.D-10/
C-----------------------------------------------------------------------
C     COMPUTE THE INITIAL VALUE OF SCATTER MATRIX A
C     ITYPW=1: STANDARDIZED CASE, A LOWER TRIANGULAR;
C     ITYPW=2: UNSTANDARDIZED CASE, A SYMMETRIC.
C-----------------------------------------------------------------------
C     X: I, DESIGN MATRIX
C     A: O, LOWER TRIANGLE OF SCATTER MATRIX A
C-----------------------------------------------------------------------
      DO 10 J=1,NCOV
         A(J)=ZERO
 10   CONTINUE
      DO 15 J=1,NVAR
         JJ=J*(J+1)/2
         A(JJ)=ONE
 15   CONTINUE
      IF (INIT.EQ.1) RETURN
      IF (ITYPW.EQ.2) GOTO 100
      DO 50 J=1,NVAR
         CALL rlLMDDbi(X(1,J),SC,NFIRST,1,XME,XMD,XSD)
         SQDEV2=DSQRT(XSD**2+XME**2)
         JJ=(J*J+J)/2
         IF (SQDEV2 .GT. TL) GOTO 40
         A(JJ)=9999.D0
         GOTO 50
 40      A(JJ)=ONE/SQDEV2
 50   CONTINUE
      RETURN
 100  DO 150 J=1,NVAR
         JJ=J*(J+1)/2
         CALL rlLMDDbi(X(1,J),SC,NFIRST,1,XME,XMD,XSD)
         DEV2=XSD**2+XME**2
         IF (DEV2 .GT. TL) GOTO 145
         A(JJ)=9999.D0
         GOTO 150
 145     A(JJ)=ONE/DEV2
 150  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlGRADbi(X,HBRS,N,NP,MDX,GRAD)
C.......................................................................
      DOUBLE PRECISION X(MDX,NP),HBRS(N),GRAD(NP),SUM
C-----------------------------------------------------------------------
      DO 20 J=1,NP
         SUM=0.D0
         DO 10 I=1,N
            SUM=SUM+X(I,J)*HBRS(I)
 10      CONTINUE
         GRAD(J)=SUM
 20   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlMSFDbi(A,B,C,N,NN,M,MDB,MDC)
C.......................................................................
      DOUBLE PRECISION A(NN),B(MDB,M),C(MDC,M),DZERO,SM
      DATA DZERO/0.D0/
C-----------------------------------------------------------------------
C     MULTIPLIES A SYMMETRIC MATRIX "A" BY A FULL MATRIX "B"
C-----------------------------------------------------------------------
C     A: I, SYMMETRIC MATRIX WITH DIMENSION NN=N*(N+1)/2,
C     B: I, A FULL MATRIX,
C     C: O, A * B.
C-----------------------------------------------------------------------
      DO 20 J=1,M
         IH=1
         DO 15 I=1,N
            L=IH
            KK=1
            SM=DZERO
            DO 10 K=1,N
               SM=SM+A(L)*B(K,J)
               IF (K .GE. I) KK=K
               L=L+KK
 10         CONTINUE
            C(I,J)=SM
            IH=IH+I
 15      CONTINUE
 20   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlMTT3bi(A,B,C,N,NN)
C.......................................................................
      DOUBLE PRECISION A(NN),B(NN),C(NN),SM,DZERO
      DATA DZERO/0.D0/
C-----------------------------------------------------------------------
C     MULTIPLIES A TRIANGULAR MATRIX BY A TRIANGULAR MATRIX
C-----------------------------------------------------------------------
      IC=0
      JJ=0
      DO 30 J=1,N
         II=0
         DO 20 I=1,J
            II=II+I
            IL=II
            IC=IC+1
            SM=DZERO
            DO 10 L=I,J
               JL=JJ+L
               SM=SM+A(IL)*B(JL)
               IL=IL+L
 10         CONTINUE
            C(IC)=SM
 20      CONTINUE
         JJ=JJ+J
 30   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlPRSFbi(SU1,NP,NCOV,TAU,INFO)
C.......................................................................
      DOUBLE PRECISION SU1(NCOV), TAU
C-----------------------------------------------------------------------
C     PRESCRIPTION F0 (TO BE USED WITH IALG=1)
C-----------------------------------------------------------------------
      CALL RLMCHLM2(SU1,NP,NCOV,INFO)
      IF (INFO.EQ.0) GOTO 100
      INFO=1
      RETURN
 100  CALL RLMINVM2(SU1,NP,NCOV,TAU,INFO)
      IF (INFO .NE. 0) INFO=2
      RETURN
      END
C=======================================================================
      SUBROUTINE rlWFAGbi(X,A,GWT,NOBS,NVAR,NVARQ,NCOV,MDX,TAU,MAXIT,
     +     ICNV,ITYPW,IGWT,TOL,NIT,DIST,SU,SA,ST,SD,SZ,IUCV,A2,B2)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NCOV),X(MDX,NVAR),DIST(NOBS),GWT(NOBS),SA(NCOV),
     +     ST(NCOV),SD(NVAR),SZ(NVAR),SU(NOBS),WUP(1)
      INTEGER rlICNVbi
      DATA ZERO,ONE/0.D0,1.D0/
C-----------------------------------------------------------------------
C     Fixed point algorithm for the computation of the matrix  A
C
C     ITYPW = 1: Standardized case,   A lower triangular;
C     ITYPW = 2: Unstandardized case, A symmetric.
c
C-----------------------------------------------------------------------
C     STEP 0 : INITIALIZATION
C-----------------------------------------------------------------------
      IALG=1
      NU=1
      NIT=0
      XN=DBLE(NOBS)
      IF (ICNV.EQ.1) THEN
         L=0
         DO 20 I=1,NVAR
            DO 10 J=1,I
               L=L+1
               SA(L)=ZERO
               IF (I.EQ.J) SA(L)=-ONE
 10         CONTINUE
 20      CONTINUE
      ENDIF
      DO 30 L=1,NOBS
         DIST(L)=-ONE
 30   CONTINUE
C-----------------------------------------------------------------------
C     STEP 1: COMPUTE WEIGHTED COVARIANCE (ST) AND AUXILIARY VALUES
C-----------------------------------------------------------------------
 100  IF (ITYPW.EQ.1) THEN
         CALL rlUCOWbi(X,A,ST,NOBS,NVAR,NVARQ,NCOV,MDX,MDX,NU,IALG,ICNV,
     +        IGWT,NIT,GWT,DELTA,DIST,SU,WUP,X,SD,IUCV,A2,B2)
      ELSE
         DO 130 I=1,NCOV
            ST(I)=ZERO
 130     CONTINUE
         DELTA=ZERO
         DO 180 L=1,NOBS
            DO 150 J=1,NVAR
               SD(J)=X(L,J)
 150        CONTINUE
            CALL rlMSFDbi(A,SD,SZ,NVAR,NCOV,1,NVAR,NVAR)
            CALL rlNRM2bi(SZ,NVAR,1,NVAR,ZNR)
            DISTL=ZNR
            IF (ICNV.NE.1) DELTA=DMAX1(DELTA,DABS(DISTL-DIST(L)))
            DIST(L)=DISTL
            U=rlUCVbi(DISTL,IUCV,A2,B2)
            SU(L)=U
            IF (IGWT.EQ.1) U=U*GWT(L)
            IJ=0
            DO 175 I=1,NVAR
               DO 170 J=1,I
                  IJ=IJ+1
                  ST(IJ)=ST(IJ)+(SD(I)*U)*SD(J)
 170           CONTINUE
 175        CONTINUE
 180     CONTINUE
         DO 190 IJ=1,NCOV
            ST(IJ)=ST(IJ)/XN
 190     CONTINUE
      ENDIF
C-----------------------------------------------------------------------
C     STEP 2: CHECK CONVERGENCE
C-----------------------------------------------------------------------
      ITMP = rlICNVbi(NCOV,DELTA,A,SA,TOL,ICNV)
      IF (NIT .EQ. MAXIT .OR. ITMP .EQ. 1) GOTO 500
C-----------------------------------------------------------------------
C     STEP 3: FIND IMPROVEMENT MATRIX ST FOR A
C-----------------------------------------------------------------------
      INFO=0
      CALL rlPRSFbi(ST,NVAR,NCOV,TAU,INFO)
C-----------------------------------------------------------------------
C     STEP 4: SET SA:=A AND A:=(ST)*SA IF ITYPW.EQ.1 ELSE A:=(ST**T)*ST
C-----------------------------------------------------------------------
      DO 410 IJ=1,NCOV
         SA(IJ)=A(IJ)
 410  CONTINUE
      IF (ITYPW.EQ.1) THEN
         CALL rlMTT3bi(SA,ST,A,NVAR,NCOV)
      ELSE
         CALL RLMTT1M2(ST,A,NVAR,NCOV)
      ENDIF
      NIT=NIT+1
      GOTO 100
 500  RETURN
      END
C=======================================================================
      SUBROUTINE rlQUNTbi(P,X)
C.......................................................................
      DOUBLE PRECISION C(6),P,P1,T,X,XN,XZ
      DATA C(1),C(2),C(3),C(4),C(5),C(6)/
     +     2.515517D0,0.802853D0,0.010328D0,
     +     1.432788D0,0.189269D0,0.001308D0/
C-----------------------------------------------------------------------
C     INVERSE OF GAUSSIAN DISTRIBUTION FUNCTION --  X := qnorm(P)
C-----------------------------------------------------------------------
C     P: I, PROBABILITY,
C     X: O, QUANTILE.
C-----------------------------------------------------------------------
      P1=P
      IF (P .GT. 0.5D0) P1=1.D0-P
      T=DSQRT(-2.D0*DLOG(P1))
      XZ=(C(3)*T+C(2))*T+C(1)
      XN=((C(6)*T+C(5))*T+C(4))*T+1.D0
      X=T-XZ/XN
      IF (P .LT. 0.5D0) X=-X
      RETURN
      END
C=======================================================================
      SUBROUTINE rlMSSDbi(A,B,C,N,NN,MDC)
C.......................................................................
      DOUBLE PRECISION A(NN),B(NN),C(MDC,N),SM,DZERO
      DATA DZERO/0.D0/
C-----------------------------------------------------------------------
C     MULTIPLIES A SYMMETRIC MATRIX BY A SYMMETRIC MATRIX -  C <- A %*% B
C-----------------------------------------------------------------------
      LI=1
      DO 60 IR=1,N
         LJ=1
         DO 50 JC=1,N
            I=LI
            J=LJ
            SM=DZERO
            DO 40 K=1,N
               SM=A(I)*B(J)+SM
               IF (K.LT.IR) GOTO 10
               I=I+K
               GOTO 20
 10            I=I+1
 20            IF (K.LT.JC) GOTO 30
               J=J+K
               GOTO 40
 30            J=J+1
 40         CONTINUE
            C(IR,JC)=SM
            LJ=LJ+JC
 50      CONTINUE
         LI=LI+IR
 60   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlMSF1bi(A,B,C,N,NN,MDB)
C.......................................................................
      DOUBLE PRECISION A(NN),B(MDB,N),C(NN),SM,DZERO
      DATA DZERO/0.D0/
C-----------------------------------------------------------------------
C     MULTIPLIES A SYMMETRIC MATRIX BY A FULL MATRIX WHEN THE RESULT
C     IS A SYMMETRIC MATRIX
C-----------------------------------------------------------------------
      LC=1
      DO 15 J=1,N
         IBEG=1
         DO 10 I=1,J
            L=IBEG
            KK=1
            SM=DZERO
            DO 5 K=1,N
               SM=SM+A(L)*B(K,J)
               IF (K.GE.I) KK=K
               L=L+KK
 5          CONTINUE
            C(LC)=SM
            LC=LC+1
            IBEG=IBEG+I
 10      CONTINUE
 15   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlLMDDbi(X,Y,N,ISORT,XME,XMD,XSD)
C.......................................................................
      DOUBLE PRECISION XME,X1,X2,XMD,XSD,X(N),Y(N),ZERO
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     COMPUTES THE MEDIAN AND MEDIAN ABSOLUTE DEVIATION (MAD):
C     IF ISORT=0, THEN X IS ALREADY SORTED IN ASCENDING ORDER.
C-----------------------------------------------------------------------
C     X  : I, ORIGINAL OBSERVATIONS
C     Y  : O, ORDERED OBSERVATIONS
C     XME: O, MEDIAN.
C     XMD: O, MAD.
C     XSD: O, MAD/0.6745.
C-----------------------------------------------------------------------
      KM=(N+1)/2
      DO 20 I=1,N
         Y(I)=X(I)
 20   CONTINUE
      IF (ISORT .NE. 0) CALL rlSRT1bi(Y,N,1,N)
      XME=Y(KM)
      IF (KM*2 .EQ. N) XME=(XME+Y(KM+1))/2.D0
      K=0
      K1=KM
      K2=KM
      X1=ZERO
      X2=ZERO
 30   IF (K .GE. KM) GOTO 50
      K=K+1
      IF (X1 .GT. X2) GOTO 40
      K1=K1-1
      IF (K1 .EQ. 0) GOTO 50
      X1=XME-Y(K1)
      GOTO 30
 40   K2=K2+1
      IF (K2 .GT. N) GOTO 50
      X2=Y(K2)-XME
      GOTO 30
 50   XMD=DMIN1(X1,X2)
      XSD=XMD/0.6745D0
      RETURN
      END
C=======================================================================
      SUBROUTINE rlUCOWbi(X,SA,ST,N,NP,NQ,NCOV,MDX,MDZ,NU,
     +     IALG,ICNV,IGWT,NIT,GWT,ZMAX,DIST,SU,SUP,SZ,SD,IUCV,A2,B2)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),DIST(N),SA(NCOV),SD(NP),ST(NCOV),SU(N),
     +     SUP(NU),SZ(MDZ,NP),GWT(N)
      DATA ZERO,XN,SQPMQ,NQP1/0.D0,2*0.D0,0/
C-----------------------------------------------------------------------
C     COMPUTE WEIGHTED COVARIANCE MATRIX F(A)
C     IF IALG=2, U' IS CALCULATED AND STORED IN SUP.
C     IF IALG=3, X (ACTUALLY Z=A*X) IS SAVED IN SZ.
C-----------------------------------------------------------------------
C     ST: O, LOWER TRIANGLE OF MATRIX F(A)
C-----------------------------------------------------------------------
      IF (NIT .GT. 1) GOTO 10
      XN=DBLE(N)
      SQPMQ=DSQRT(DBLE(NP-NQ))
      NQP1=NQ+1
 10   ZMAX=ZERO
      DO 20 IJ=1,NCOV
         ST(IJ)=ZERO
 20   CONTINUE
      DO 100 L=1,N
         DO 50 J=1,NP
            SD(J)=X(L,J)
 50      CONTINUE
         CALL rlMLYDbi(SA,SD,NP,NCOV,NP,1)
         CALL rlNRM2bi(SD(NQP1),NP-NQ,1,NP-NQ,ZNR)
         DISTL = ZNR
         IF (NQ .NE. 0) DISTL=DISTL/SQPMQ
         IF (ICNV .EQ. 2) ZMAX=DMAX1(ZMAX,DABS(DISTL-DIST(L)))
         DIST(L)=DISTL
         U=rlUCVbi(DISTL,IUCV,A2,B2)
         SU(L)=U
         IF (IGWT .EQ. 1) U=U*GWT(L)
         IF (IALG .EQ. 1) GOTO 80
         UP=rlUPCVbi(DISTL,IUCV,A2,B2)
         IF (NQ .NE. 0) UP=UP/SQPMQ
         SUP(L)=UP
         IF (IALG .EQ. 2) GOTO 80
         DO 70 I=1,NP
            SZ(L,I)=SD(I)
 70      CONTINUE
 80      IJ=0
         DO 95 I=1,NP
            DO 90 J=1,I
               IJ=IJ+1
               ST(IJ)=ST(IJ)+(SD(I)*U)*SD(J)
 90         CONTINUE
 95      CONTINUE
 100  CONTINUE
      DO 110 IJ=1,NCOV
         ST(IJ)=ST(IJ)/XN
 110  CONTINUE
      RETURN
      END
C=======================================================================
      FUNCTION rlICNVbi(NCOV,DELTA,SA,SA0,TOL,ICNV)
C.......................................................................
      DOUBLE PRECISION SA(NCOV),SA0(NCOV),SDMAX,DELTA,TOL
      INTEGER rlICNVbi
C-----------------------------------------------------------------------
C     CHECK CONVERGENCE
C-----------------------------------------------------------------------
      rlICNVbi=0
      IF (ICNV .EQ. 1) THEN
         DO 10 IJ=1,NCOV
            SA0(IJ)=SA(IJ)-SA0(IJ)
 10      CONTINUE
         CALL rlNRM2bi(SA0,NCOV,1,NCOV,SDMAX)
         DELTA=SDMAX
      ENDIF
      IF (DELTA .LT. TOL) rlICNVbi=1
      RETURN
      END
C=======================================================================
      SUBROUTINE rlC0HKbi(X,N,NP,MDX,CONST)
C.......................................................................
      DOUBLE PRECISION SUMNRM,CONST,XNRM,X(MDX,NP)
C-----------------------------------------------------------------------
C     INITIALIZE THE TUNING CONSTANT FOR HAMPEL-KRASKER ESTIMATOR
C-----------------------------------------------------------------------
      SUMNRM=0.D0
      DO 10 I=1,N
         CALL rlNRM2bi(X(I,1),NP,MDX,MDX*NP-I+1,XNRM)
         SUMNRM=SUMNRM+XNRM
 10   CONTINUE
      CONST=DBLE(NP)*DSQRT(1.5707963D0)/(SUMNRM/DBLE(N))
      RETURN
      END
C=======================================================================
      SUBROUTINE rlC0MUbi(X,N,NP,MDX,CONST)
C.......................................................................
      DOUBLE PRECISION SUMNRM,CONST,XNRM,X(MDX,NP)
C-----------------------------------------------------------------------
C     INITIALIZE THE TUNING CONSTANT FOR MALLOWS(U) ESTIMATOR
C-----------------------------------------------------------------------
      SUMNRM=0.D0
      DO 10 I=1,N
         CALL rlNRM2bi(X(I,1),NP,MDX,MDX*NP-I+1,XNRM)
         SUMNRM=SUMNRM+XNRM
 10   CONTINUE
      CONST=DBLE(NP)/(SUMNRM/DBLE(N))
      RETURN
      END
C=======================================================================
      SUBROUTINE rlFUDGbi(SS,NP,NCOV,XKAP,GAMMA)
C.......................................................................
      DOUBLE PRECISION SS(NCOV),E,XKAP,GAMMA,SII
C-----------------------------------------------------------------------
C     COMPUTE FUDGE FACTOR GAMMA
C-----------------------------------------------------------------------
      E=0.D0
      DO 10 I=1,NP
         II=I*(I+1)/2
         SII=SS(II)
         E=DMAX1(E,DABS(SII))
 10   CONTINUE
      GAMMA=1.D0/DMAX1(1.D0,XKAP*E)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlHUBbi(RS,WGT,WGT2,SIGMB,N,ITYPE,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RS(N),WGT(N),WGT2(N)
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     COMPUTES HUBERIZED RESIDUALS IN THE HUBER'S ALGORITHM.
C-----------------------------------------------------------------------
      IF (ITYPE .NE. 1) GOTO 20
C-----------------------------------------------------------------------
C     HUBER-TYPE
C-----------------------------------------------------------------------
      DO 10 I=1,N
         S=RS(I)/SIGMB
         RS(I)=RLPSIM2(S,IPS,XK)*SIGMB
 10   CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     MALLOWS-TYPE (WGT2 is the square root of the WEIGHTS)
C-----------------------------------------------------------------------
 20   IF (ITYPE .NE. 2) GOTO 40
      DO 30 I=1,N
         SW=WGT2(I)*SIGMB
         IF (SW .GT. ZERO) GOTO 25
         RS(I)=ZERO
         GOTO 30
 25      S=RS(I)/SIGMB
         RS(I)=RLPSIM2(S,IPS,XK)*SW
 30   CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     SCHWEPPE-TYPE
C-----------------------------------------------------------------------
 40   DO 50 I=1,N
         SW=WGT(I)*SIGMB
         IF (SW .GT. ZERO .AND. WGT(I) .GT. ZERO) GOTO 45
         RS(I)=ZERO
         GOTO 50
 45      S=RS(I)/SW
         RS(I)=RLPSIM2(S,IPS,XK)*SW
 50   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlKEDCbi(WGT,RS,N,SIGMA,ITYPE,D,E,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N),RS(N),D(N),E(N)
      DATA ZERO,ONE/0.D0,1.D0/
C-----------------------------------------------------------------------
C     COMPUTES THE DIAGONAL MATRIX
C-----------------------------------------------------------------------
      IF (ITYPE .EQ. 3) GOTO 20
C-----------------------------------------------------------------------
C     MALLOWS CASE
C-----------------------------------------------------------------------
      DO 10 I=1,N
         IF (WGT(I) .GT. ZERO) GOTO 5
         D(I)=-ONE
         E(I)=ZERO
         GOTO 10
 5       X=RS(I)/SIGMA
         D(I)=RLPSPM2(X,IPS,XK)*WGT(I)
         E(I)=(RLPSIM2(X,IPS,XK)*WGT(I))**2
 10   CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     SCHWEPPE CASE
C-----------------------------------------------------------------------
 20   DO 40 I=1,N
         IF (WGT(I) .GT. ZERO) GOTO 30
         D(I)=-ONE
         E(I)=ZERO
         GOTO 40
 30      X=RS(I)/SIGMA/WGT(I)
         D(I)=RLPSPM2(X,IPS,XK)
         E(I)=(RLPSIM2(X,IPS,XK)*WGT(I))**2
 40   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlPRSHbi(SU1,SS,DIST,SU,SUP,SV,SVPZ,N,NP,NCOV)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DIST(N),SU1(NCOV),SS(NCOV),SU(N),SUP(N)
      DATA TL/1.D-10/
C-----------------------------------------------------------------------
C     MODIFIED NEWTON PRESCRIPTION (FOR IALG=2 FOR PRESCRIPTION)
C-----------------------------------------------------------------------
      XN=DBLE(N)
      XP=DBLE(NP)
      C1=SVPZ/XN
      D1=SV/XN
      A=0.D0
      B=0.D0
      GOTO 20
 10   XPDEN=1.D0
      BTMD=-D1
      GOTO 40
 20   DO 30 L=1,N
         ZNR=DIST(L)
         ZNR2=ZNR*ZNR
         A=A+SU(L)*ZNR2
         B=B+SUP(L)*ZNR*ZNR2
 30   CONTINUE
      A1=A/XN
      B1=B/XN/(XP+2.D0)
      DEN=A1+B1
      IF (DABS(DEN) .LT. TL) GOTO 10
      TDEN=2.D0*DEN+XP*(B1-C1)
      IF (DABS(TDEN) .LT. TL) GOTO 10
      T=(D1*XP-A1)/TDEN
      XPDEN=XP/DEN
      BTMD=(B1-C1)*T-D1
 40   IJ=0
      DO 70 I=1,NP
         IF (I.EQ.1) GOTO 60
         I1=I-1
         DO 50 J=1,I1
            IJ=IJ+1
            SS(IJ)=SU1(IJ)*XPDEN
 50      CONTINUE
 60      IJ=IJ+1
         SS(IJ)=(SU1(IJ)+BTMD)*(XPDEN/2.D0)
 70   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlBETHbi(WGT,N,D,ITYPE,BETA)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION WGT(N)
C-----------------------------------------------------------------------
C     COMPUTES THE CONSTANT BETA FOR HUBER FUNCTION
C-----------------------------------------------------------------------
      SM=0.D0
      XN=DBLE(N)
      C2=D*D
      IF (ITYPE.EQ.3) GOTO 30
C-----------------------------------------------------------------------
C     HUBER CASE
C-----------------------------------------------------------------------
      CALL rlGAUSbi(D,PC)
      CALL rlXERFbi(2,D,DC)
      BETA=-D*DC+PC-0.5D0+C2*(1.D0-PC)
      IF (ITYPE.EQ.1) RETURN
C-----------------------------------------------------------------------
C     MALLOWS CASE
C-----------------------------------------------------------------------
      DO 20 I=1,N
         SM=SM+WGT(I)
 20   CONTINUE
      BETA=SM*BETA/XN
      RETURN
C-----------------------------------------------------------------------
C     SCHWEPPE CASE
C-----------------------------------------------------------------------
 30   DO 40 I=1,N
         W2=WGT(I)*WGT(I)
         CW=D*WGT(I)
         CALL rlGAUSbi(CW,PC)
         CALL rlXERFbi(2,CW,DC)
         B=C2*(1.D0-PC)+(-CW*DC+PC-0.5D0)/W2
         SM=SM+B*W2/XN
 40   CONTINUE
      BETA=SM
      RETURN
      END

C=======================================================================
      SUBROUTINE rlRNAGbi(X,Y,THETA,WGT,COV,SIGMAI,N,NP,MDX,MDT,NCOV,
     +     GAM,TOL,TAU,ITYPE,IOPT,ISIGMA,ICNV,MAXIT,MAXIS,NIT,SIGMAF,
     +     QS1,RS,DELTA,GRAD,HESSNV,SD,SW,SF,SG,SH,IP,SX,IPS,XK,
     +     BETA,BET0)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(MDT),WGT(N),COV(NCOV),RS(N),
     +     DELTA(NP),GRAD(NP),SD(N),HESSNV(NCOV),
     +     SF(NP),SG(NP),SH(NP),SW(N),SX(MDX,NP)
      INTEGER IP(NP),RLICTHM2,RLISIGM2
      LOGICAL FIRST
      DATA ZERO,TL/0.D0,1.D-10/
C----------------------------------------------------------------------
C     MODIFIED NEWTON ALGORITHMS FOR ROBUST AND BOUNDED INFLUENCE
C     LINEAR REGRESSION
C----------------------------------------------------------------------
      NN=NP*(NP+1)/2
      XN=DBLE(N)
      SIGMA=SIGMAI
      SIGMB=SIGMA
      IASG=IABS(ISIGMA)
      PSP0=RLPSPM2(ZERO,IPS,XK)
      INTCH=1
      ITYP=ITYPE
C----------------------------------------------------------------------
      IF (ITYP .EQ. 1) GOTO 15
      N0=N
      E=2.D0
      IF (ITYP .EQ. 2) E=0.5D0
      DO 10 I=1,N
         IF (WGT(I) .LE. ZERO) THEN
            SW(I)=-1.D0
            N0=N0-1
         ELSE
            SW(I)=WGT(I)**E
         ENDIF
 10   CONTINUE
      IF (N0 .EQ. 0) ITYP=1
 15   IF (IASG .EQ. 0) CONST=ZERO
      IF (IASG .EQ. 1) CONST=BETA*DBLE(N-NP)
      IF (IASG .EQ. 2) CONST=BET0*DBLE(N-NP)
C----------------------------------------------------------------------
C     STEP 1. SET NIT:=1
C----------------------------------------------------------------------
      NIT=1
C----------------------------------------------------------------------
C     STEP 2. COMPUTE RESIDUALS
C----------------------------------------------------------------------
 200  CALL RLRESDM2(X,Y,THETA,N,NP,MDX,RS)
C----------------------------------------------------------------------
C     STEP 3. COMPUTE A NEW VALUE SIGMB FOR SIGMA
C----------------------------------------------------------------------
      IF (ISIGMA .LT. 0 .AND. NIT .EQ. 1) GOTO 300
      IF (ISIGMA .EQ. 0) GOTO 300
      SIGMA=SIGMB
      CALL RLRSIGM2(RS,WGT,SIGMA,N,NP,TOL,ITYP,ISIGMA,MAXIS,
     +     NIS,SIGMB,SW,SD,IPS,XK,BETA,BET0)
      IF (SIGMB .GT. TL) GOTO 300
      RETURN
 300  CALL RLQRSSM2(RS,WGT,SW,N,ITYP,SIGMB,CONST,QS0,IPS,XK)
C----------------------------------------------------------------------
C     STEP 4. COMPUTE THE (UNSCALED) NEGATIVE GRADIENT
C----------------------------------------------------------------------
      DO 410 I=1,N
         SD(I)=RS(I)
 410  CONTINUE
      CALL rlHUBbi(SD,WGT,WGT,SIGMB,N,ITYP,IPS,XK)
      CALL rlGRADbi(X,SD,N,NP,MDX,GRAD)
C----------------------------------------------------------------------
C     STEP 5. COMPUTE WEIGHTS AND APPLY THEM TO X; STORE RESULT IN SX
C----------------------------------------------------------------------
      DO 550 I=1,N
         T=RS(I)/SIGMB
         IF (ITYP .EQ. 1) GOTO 510
         SDI=0.D0
         IF (WGT(I) .LE. ZERO) GOTO 520
         IF (ITYP .EQ. 3) T=T/WGT(I)
 510     SDI=RLPSPM2(T,IPS,XK)
 520     SQD=DSQRT(SDI)
         IF (ITYP .EQ. 2) SQD=SQD*SW(I)
         DO 530 J=1,NP
            SX(I,J)=X(I,J)*SQD
 530     CONTINUE
 550  CONTINUE
      FIRST=.TRUE.
C----------------------------------------------------------------------
C     STEP 6. COMPUTE GENERALIZED INVERSE OF UNSCALED HESSIAN MATRIX
C----------------------------------------------------------------------
 600  CALL RLRMTRM2(SX,N,NP,MDX,INTCH,TAU,K,SF,SG,SH,IP)
      IF (K.EQ.0) RETURN
      CALL RLKIASM2(SX,K,NP,MDX,NN,1.D0,1.D0,HESSNV)
      CALL RLKFASM2(SX,HESSNV,K,NP,MDX,NN,1.D0,DELTA,SG,IP)
C----------------------------------------------------------------------
C     STEP 7. COMPUTE THE INCREMENT VECTOR
C----------------------------------------------------------------------
      CALL rlMSFDbi(HESSNV,GRAD,DELTA,NP,NN,1,NP,NP)
      DO 710 J=1,NP
         DELTA(J)=GAM*DELTA(J)
         IF (FIRST) THEN
            SD(J)=THETA(J)
            THETA(J)=THETA(J)+DELTA(J)
         ELSE
            SF(J)=THETA(J)
            THETA(J)=SD(J)+DELTA(J)
         ENDIF
 710  CONTINUE
C----------------------------------------------------------------------
C     STEP 8. IF NECESSARY DETERMINE ANOTHER STEP LENGTH.
C----------------------------------------------------------------------
      CALL RLRESDM2(X,Y,THETA,N,NP,MDX,RS)
      IF (.NOT.FIRST) QSF=QS1
      CALL RLQRSSM2(RS,WGT,SW,N,ITYP,SIGMB,CONST,QS1,IPS,XK)
      IF (QS1 .LE. QS0) GOTO 900
      IF (.NOT.FIRST) GOTO 880
      FIRST=.FALSE.
      IF (IOPT .EQ. 1) GOTO 800
      IF (IOPT .EQ. 2) GOTO 830
C----------------------------------------------------------------------
C     H-ALGORITHM OPTION
C----------------------------------------------------------------------
 800  IF (ITYP .EQ. 2) THEN
         DO 815 J=1,NP
            DO 810 I=1,N
               SX(I,J)=X(I,J)*SW(I)
 810        CONTINUE
 815     CONTINUE
      ELSE
         DO 825 J=1,NP
            DO 820 I=1,N
               SX(I,J)=X(I,J)
 820        CONTINUE
 825     CONTINUE
      ENDIF
      GOTO 600
C----------------------------------------------------------------------
C     W-ALGORITHM OPTION
C----------------------------------------------------------------------
 830  DO 850 I=1,N
         SDI=PSP0
         IF (RS(I) .EQ. ZERO) GOTO 840
         T=RS(I)/SIGMB
         IF (ITYP .EQ. 1) GOTO 835
         SDI=ZERO
         IF (WGT(I) .LE. ZERO) GOTO 840
         IF (ITYP .EQ. 2) GOTO 835
         T=T/WGT(I)
 835     SDI=RLPSIM2(T,IPS,XK)/T
 840     PI=DSQRT(SDI)
         IF (ITYP .EQ. 2) PI=PI*SW(I)
         DO 845 J=1,NP
            SX(I,J)=PI*X(I,J)
 845     CONTINUE
 850  CONTINUE
      GOTO 600
 880  IF (QSF .LT. QS1) THEN
         DO 890 J=1,NP
            THETA(J)=SF(J)
 890     CONTINUE
      ENDIF
C----------------------------------------------------------------------
C     STEP 9. STOP ITERATIONS IF DESIRED PRECISION HAS BEEN REACHED
C----------------------------------------------------------------------
 900  IF (ISIGMA .LT. 0 .AND. NIT .EQ. 1) GOTO 950
      IF (NIT .EQ. MAXIT) GOTO 990
      IF (RLICTHM2(NP,NCOV,DELTA,SIGMA,COV,TOL,ICNV) .EQ. 1
     +     .AND. RLISIGM2(SIGMA,SIGMB,TOL) .EQ. 1) GOTO 990
 950  NIT=NIT+1
      GOTO 200
 990  SIGMAF=SIGMB
      QS1=QS1/XN
      RETURN
      END
C=======================================================================
      SUBROUTINE rlUDATbi(SS,SA0,SA,GAMMA,NP,NCOV)
C.......................................................................
      DOUBLE PRECISION SS(NCOV),SA0(NCOV),SA(NCOV),GAMMA,GAMD
      IJ=0
      GAMD=-GAMMA
      DO 20 I=1,NP
         DO 10 J=1,I
            IJ=IJ+1
            SA(IJ)=GAMD*SS(IJ)
            IF (I.EQ.J) SA(IJ)=1.D0+SA(IJ)
 10      CONTINUE
 20   CONTINUE
      CALL rlMTT3bi(SA0,SA,SA,NP,NCOV)
      RETURN
      END

C=======================================================================
      SUBROUTINE rlWNAGbi(X,A,NOBS,NVAR,NCOV,MDX,MAXIT,ICNV,TOL,
     +     XFUD,NIT,DIST,SA,SS,SU,SUP,ST,SD,IUCV,A2,B2)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NVAR),DIST(NOBS),A(NCOV),SA(NCOV),SS(NCOV),
     +     ST(NCOV),SU(NOBS),SUP(NOBS),SD(NVAR)
      INTEGER rlICNVbi
      DATA ZERO,ONE/0.D0,1.D0/
C-----------------------------------------------------------------------
C     NEWTON-HUBER ALGORITHM FOR THE COMPUTATION OF SCATTER MATRIX A
C     (STANDARDIZED CASE, A LOWER TRIANGULAR)
C-----------------------------------------------------------------------
C     STEP 0 : INITIALIZATION
C-----------------------------------------------------------------------
      IALG=2
      NU=NOBS
      XN=DBLE(NOBS)
      NIT=0
      NVARQ=0
      IF (ICNV.EQ.1) THEN
         L=0
         DO 20 I=1,NVAR
            DO 10 J=1,I
               L=L+1
               SA(L)=ZERO
               IF (I.EQ.J) SA(L)=-ONE
 10         CONTINUE
 20      CONTINUE
      ENDIF
      DO 30 L=1,NOBS
         DIST(L)=-ONE
 30   CONTINUE
C-----------------------------------------------------------------------
C     STEP 1: COMPUTE WEIGHTED COVARIANCE (ST) AND AUXILIARY VALUES
C-----------------------------------------------------------------------
 100  CALL rlUCOWbi(X,A,ST,NOBS,NVAR,NVARQ,NCOV,MDX,MDX,
     +     NU,IALG,ICNV,0,NIT,DIST,DELTA,DIST,SU,SUP,X,SD,IUCV,A2,B2)
C-----------------------------------------------------------------------
C     STEP 2: CHECK CONVERGENCE
C-----------------------------------------------------------------------
      IF (NIT .EQ. MAXIT .OR. rlICNVbi(NCOV,DELTA,A,SA,TOL,ICNV) .EQ. 1)
     +     GOTO 500
C-----------------------------------------------------------------------
C     STEP 3: FIND IMPROVEMENT MATRIX SS FOR A
C-----------------------------------------------------------------------
      CALL rlPRSHbi(ST,SS,DIST,SU,SUP,XN,0.D0,NOBS,NVAR,NCOV)
C-----------------------------------------------------------------------
C     STEP 4: COMPUTE GAM0 AND SET A:=(I-GAM0*SS)*SA
C-----------------------------------------------------------------------
      DO 410 IJ=1,NCOV
         SA(IJ)=A(IJ)
 410  CONTINUE
      CALL rlFUDGbi(SS,NVAR,NCOV,XFUD,GAM0)
      CALL rlUDATbi(SS,SA,A,GAM0,NVAR,NCOV)
      NIT=NIT+1
      GOTO 100
  500 RETURN
      END

C======================================================================
      SUBROUTINE rlWWWAbi(N,SVALS,FVALS,IWWW,IUCV,A2,B2)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SVALS(N),FVALS(N)
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF THE W BAR-FUNCTION
C-----------------------------------------------------------------------
      DO 50 I=1,N
         FVALS(I)=rlWWWbi(SVALS(I),IWWW,IUCV,A2,B2)
 50   CONTINUE
      RETURN
      END
C=======================================================================







