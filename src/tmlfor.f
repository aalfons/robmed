C
C=======================================================================
      SUBROUTINE RLHSESM3(X,Y,N,NP,NQ,NCOV,MDX,MDW,MDI,IOPT,INTCH,NREP,
     +     TOLS,TOLR,TAU,GAM,MAXIT,MAXS1,MAXS2,ISEED,IERR,SMIN,
     +     THETA,RS,IT1,COV,WORK,IWORK,IPS,XK,BETA,BET0,ITRACE)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(N),RS(N),IT1(NQ),COV(NCOV)
      DIMENSION WORK(MDW),IWORK(MDI)
C-----------------------------------------------------------------------
C     Resampling algorithm for the computation of S-estimates
C-----------------------------------------------------------------------
      NP1=NP+1
      N0=NP*NQ+1
      N1=N0+NQ
      N2=N1+NQ
      N3=N2+NP
      N4=N3+NP
      N5=N4+NP
      N6=N5+MDX*NP
      CALL RLHSE2M3(X,Y,N,NP,NQ,NCOV,MDX,IOPT,INTCH,NREP,TOLS,TOLR,TAU,
     +  GAM,MAXIT,MAXS1,MAXS2,ISEED,IERR,SMIN,THETA,RS,IT1,COV,
     +  WORK(1),WORK(N0),WORK(N1),WORK(N2),WORK(N3),WORK(N4),WORK(N5),
     +  WORK(N6),WORK(1),IWORK(1),IWORK(NP1),IPS,XK,BETA,BET0,ITRACE)
      RETURN
      END
c
      SUBROUTINE RLHSE2M3(X,Y,N,NP,NQ,NCOV,MDX,IOPT,INTCH,NREP,
     +   TOLS,TOLR,TAU,GAM,MAXIT,MAXS1,MAXS2,ISEED,IERR,SMIN,THETA,
     +   RS,IT1,COV,XX,YY,XTHETA,SF,SG,SH,SX,SZ,SA,SP,IT,
     +   IPS,XK,BETA,BET0,ITRACE)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(N),RS(N),COV(NCOV),XX(NQ,NP),
     * YY(NQ),XTHETA(NQ),SF(NP),SG(NP),SH(NP),SX(MDX,NP),SZ(N),SA(NCOV)
      INTEGER IT1(NQ),SP(NP),IT(NQ)
      EXTERNAL RLPSPM2,RLPSIM2,RLCHIM2,ICNREP
      LOGICAL ALLZERO
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     Resampling algorithm for the computation of S-estimates
C-----------------------------------------------------------------------
C     STEP 0: INITIALIZATIONS
C-----------------------------------------------------------------------
      K1=N/2+1
      CONST=BETA*DFLOAT(N-NP)
      IF (ITRACE .EQ. 1) THEN
         call mminitclk(iclock)
         ielapse=0
         if (np .le. 10) then
            ninc=15000
         else if (np .le. 15) then
            ninc=10000
         else
            ninc=5000
         endif
      endif
      NIT=1
      IERR=2
      SMIN=ZERO
      ITYPE=1
      IF (IOPT.NE.2) NREP=ICNREP(N,NQ,IOPT,0)
      PSP0=RLPSPM2(ZERO,IPS,XK)
C      DO 10 I=1,NP
C  10  SP(I)=I
      DO I=1,NP
        SP(I)=I
      END DO
C-----------------------------------------------------------------------
C     STEP 1: DRAW A SUBSAMPLE
C-----------------------------------------------------------------------
 100  IF (IOPT.NE.3) THEN
         DO 130 K=1,NQ
 110        CALL RLRNDM2(ISEED,RND)
            ITK=idint(RND*N)+1
            DO 120 KK=1,K-1
               IF (ITK.EQ.IT(KK)) GOTO 110
 120        CONTINUE
            IT(K)=ITK
 130     CONTINUE
      ELSE
         IF (NIT.EQ.1) THEN
            DO 140 K=1,NQ
               IT(K)=K
 140        CONTINUE
         ELSE
            CALL RLNCOMM2(N,NQ,IT)
         ENDIF
      ENDIF
      DO 160 K=1,NQ
         ITK=IT(K)
         DO 150 J=1,NP
            XX(K,J)=X(ITK,J)
 150     CONTINUE
         YY(K)=Y(ITK)
 160  CONTINUE
C-----------------------------------------------------------------------
C     STEP 2: DECOMPOSE SAMPLE MATRIX
C-----------------------------------------------------------------------
      ALLZERO = .TRUE.
      K = 1
      DO WHILE (ALLZERO .AND. K .LE. NQ)
         IF (YY(K) .NE. ZERO) ALLZERO = .FALSE.
         K = K + 1
      ENDDO
      IF (ALLZERO) GOTO 700
      CALL RLRMTRM2(XX,NQ,NP,NQ,INTCH,TAU,KK,SF,SG,SH,SP)
      IF(KK.NE.NP) GOTO 700
C-----------------------------------------------------------------------
C     STEP 3: SOLVE SYSTEM OF LINEAR EQUATIONS
C-----------------------------------------------------------------------
      CALL RLRICLM2(XX,YY,NQ,NP,NQ,XTHETA,SH,SP)
C-----------------------------------------------------------------------
C     STEP 4: COMPUTE RESIDUALS
C-----------------------------------------------------------------------
      DO 420 I=1,N
         S=Y(I)
         DO 410 J=1,NP
            S=S-XTHETA(J)*X(I,J)
 410     CONTINUE
         RS(I)=S
 420  CONTINUE
      IF (SMIN .EQ. ZERO) THEN
         S=1.0D7
         DO 430 I=1,N
            ARI=DABS(RS(I))
            SZ(I)=ARI
            IF (ARI .NE. ZERO) S=DMIN1(S,ARI)
 430     CONTINUE
         IF (S .EQ. 1.0D7) GOTO 915
         CALL RLSTORM2(SZ,N,K1,S0)
         S0=2.D0*S0
         IF (S0 .EQ. ZERO) S0=S
         SRES=S0
      ENDIF
 435  D=ZERO
      DO 440 I=1,N
         D=D+RLCHIM2(RS(I)/SRES,IPS,XK)
 440  CONTINUE
      IF (SMIN .NE. ZERO .AND. D .GT. CONST) GOTO 700
      IF (D .LE. CONST) GOTO 500
      S0=1.5D0*S0
      SRES=S0
      GOTO 435
C-----------------------------------------------------------------------
C     STEP 5: SOLVE FOR SRES
C-----------------------------------------------------------------------
 500  CALL RLRSIGM2(RS,SZ,S0,N,NP,TOLR,ITYPE,1,MAXS1,NIS,SRES,SZ,
     +     SZ,IPS,XK,BETA,BET0)
C-----------------------------------------------------------------------
C     STEP 6: UPDATE BEST FIT
C-----------------------------------------------------------------------
      IERR=0
      SMIN=SRES
      S0=SMIN
      DO 610 K=1,NP
         THETA(K)=XTHETA(K)
 610  CONTINUE
      DO 620 K=1,NQ
         IT1(K)=IT(K)
 620  CONTINUE
      IF (SRES .LE. TOLS) THEN
         IERR=1
         GOTO 800
      ENDIF
C-----------------------------------------------------------------------
C     STEP 7: END OF MAIN LOOP
C-----------------------------------------------------------------------
 700  IF (NIT.EQ.NREP) GOTO 800
      IF (ITRACE .EQ. 1) THEN
         itmp = nit/ninc
         if (itmp .gt. 0 .and. nit-ninc*itmp .eq. 0) then
            call mmprint(nrep,itmp,iclock,ielapse,ninc)
         endif
      ENDIF
      NIT=NIT+1
      GOTO 100
C-----------------------------------------------------------------------
C     STEP 8: EXIT
C-----------------------------------------------------------------------
 800  IF (IERR.EQ.2) RETURN
      DO 820 I=1,N
         S=Y(I)
         DO 810 J=1,NP
            S=S-THETA(J)*X(I,J)
 810     CONTINUE
         RS(I)=S
 820  CONTINUE
      K=1
      MAXIW=1
      ISIGMA=-1
  830 SWI=0.D0
C      DO 860 I=1,N
      DO I=1,N
        WI=0.D0
        IF (RS(I).EQ.0.D0) GOTO 840
        T=RS(I)/SMIN
        WI=RLPSIM2(T,IPS,XK)/T
        SWI=SWI+WI
        WI=SQRT(WI)
C  840   DO 850 J=1,NP
C  850   SX(I,J)=WI*X(I,J)
  840   DO J=1,NP
          SX(I,J)=WI*X(I,J)
        END DO
      END DO
C  860 CONTINUE
      CALL RLKFFAM2(RS,N,NP,SMIN,FH,IPS,XK)
      FACT=FH*SWI/DFLOAT(N)
      IF (K.EQ.0) FACT=FACT*SMIN*SMIN
      CALL RLKTASM2(SX,N,NP,MDX,NCOV,TAU,FACT,SA,COV)
c        call dblepr('cov',3,cov,ncov)
      IF (K.EQ.0) RETURN
      SRES=SMIN
      ICNV=1
C      DO 870 J=1,NP
C  870 XTHETA(J)=THETA(J)
      DO J=1,NP
        XTHETA(J)=THETA(J)
      END DO
      IF (MAXIW.EQ.1) CALL RLQRSHM2(RS,N,NP,SRES,QR0,IPS,XK)
  880  CALL RLRWAGM2(X,Y,THETA,SZ,COV,PSP0,SRES,N,NP,MDX,NCOV,
     *      TOLR,GAM,TAU,ITYPE,ISIGMA,ICNV,MAXIW,MAXS2,NIT8,
     *      SMIN,RS,YY,SZ,SF,SG,SH,SP,SZ,SX,IPS,XK,BETA,BET0)

C
C STEP 9: EXIT
C -------
      IF (MAXIW.EQ.1) THEN
        CALL RLQRSHM2(RS,N,NP,SRES,QR1,IPS,XK)
        IF (QR0.LE.QR1) GOTO 910
        ISIGMA=1
        MAXIW=MAXIW+MAXIT
        GOTO 880
      ENDIF
c
      IF (SMIN.LT.SRES) THEN
        K=0
        GOTO 830
      ENDIF
  910 CONTINUE
      SMIN=SRES
      FACT=SMIN*SMIN
      CALL RLSCALM2(COV,FACT,NCOV,1,NCOV)
C  915 DO 920 J=1,NP
C  920 THETA(J)=XTHETA(J)
  915 DO J=1,NP
        THETA(J)=XTHETA(J)
      END DO
C      DO 940 I=1,N
C      S=Y(I)
C      DO 930 J=1,NP
C  930 S=S-THETA(J)*X(I,J)
C  940 RS(I)=S
      DO I=1,N
        S=Y(I)
        DO J=1,NP
          S=S-THETA(J)*X(I,J)
        END DO
        RS(I)=S
      END DO
      RETURN
      END
c
      FUNCTION ICNREP(N,NP,IOPT,IMODE)
C
C  COMPUTE NUMBER OF REPETITIONS IN RLHSESM3
C  M  NUMBER OF OBSERVATIONS
C  NP NUMBER OF PARAMETERS
C  IOPT  0  QUICK VERSION
C        1  EXTENDED VERSION
C        2  NOT USED
C        3  ALL COMBINATIONS
C
      DIMENSION NREPQ(8),NREPE(5)
      DATA NREPQ/150,300,400,500,600,700,850,1250/
      DATA NREPE/500,1000,1500,2000,2500/
C      GOTO (1,2,3,4) IOPT+1
      if (IOPT + 1 == 2) then
        goto 2
      else if (IOPT + 1 == 3) then
        goto 3
      else if (IOPT + 1 == 4) then
        goto 4
      else
        goto 1
      end if
    1 IF(NP .GE. 9) THEN
         ICNREP=1500
      ELSE
         ICNREP=NREPQ(NP)
      ENDIF
      RETURN
    2 IF(NP .GE. 6) THEN
         ICNREP=3000
      ELSE
         ICNREP=NREPE(NP)
      ENDIF
    3 RETURN
    4 NN=N
      NR=1
C      DO 10 I=1,NP
C         NR=(NR*NN)/I
C   10    NN=NN-1
      DO I=1,NP
         NR=(NR*NN)/I
         NN=NN-1
      END DO
      IF (IMODE.GE.3) NR=NR*2**(NP-1)
      ICNREP=NR
      RETURN
      END
C
C-----------------------------------------------------------------------
C
C      FUNCTION RLICTHM2(NP,NCOV,DELTA,SIGMA,S,TOL,ICNV)
C.......................................................................
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      INTEGER RLICTHM2
C      DIMENSION DELTA(NP),S(NCOV)
C-----------------------------------------------------------------------
C      RLICTHM2=0
C      TOL1=TOL*SIGMA
C      IF (ICNV.EQ.2) GOTO 200
C      IF (ICNV.EQ.3) GOTO 300
C      L=0
C      DO 100 J=1,NP
C         L=L+J
C         TOL2=TOL1*DSQRT(S(L))
C         IF (TOL2 .LT. DABS(DELTA(J))) RETURN
C 100  CONTINUE
C      GOTO 500
C 200  CALL RLXSYM2(DELTA,DELTA,S,NP,NCOV,TOL2)
C      TOL2=DSQRT(TOL2)
C      IF (TOL1 .GE. TOL2) RLICTHM2=1
C      RETURN
C 300  L=0
C      DO 350 J=1,NP
C         L=L+J
C         TOL2=DABS(DELTA(J))*DSQRT(S(L))
C         IF (TOL1 .LT. TOL2) RETURN
C 350  CONTINUE
C 500  RLICTHM2=1
C      RETURN
C      END
C
C-----------------------------------------------------------------------
C
C      SUBROUTINE RLINTGRT(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,
C     1           LIMIT,RESULT,ABSERR,NEVAL,IER,WORK,IWORK,NPR,PARAM)
C
C.......................................................................
C
C   AUTHOR :     QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
C      implicit double precision(a-h, o-z)
C      DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,RESULT,FEXT,WORK
C      INTEGER IER,KEY,LAST,LIMIT,NEVAL,ALIST,BLIST,ELIST,RLIST
C
C      DIMENSION FARR(N),WORK(4*LIMIT),IWORK(LIMIT),PARAM(NPR)
C
C      EXTERNAL F,FEXT,GEXT
C
C         LIMIT IS THE MAXIMUM NUMBER OF SUBINTERVALS ALLOWED IN THE
C         SUBDIVISION PROCESS OF QAGE1D. TAKE CARE THAT LIMIT.GE.1.
C
C**** DATA LIMIT/500/
C
C      IF ((EPSABS.LT.0..AND.EPSREL.LT.0.).OR.LIMIT.LE.1
C     1   .OR.LIMIT.GT.500)  CALL MESSGE(500,'INTGRD',1)
C      ALIST=1
C      BLIST=ALIST+LIMIT
C      RLIST=BLIST+LIMIT
C      ELIST=RLIST+LIMIT
C      CALL RLQAGE1T(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,LIMIT,
C     1     RESULT,ABSERR,NEVAL,IER,
C     2     WORK,WORK(BLIST),WORK(RLIST),WORK(ELIST),IWORK,LAST,
C     2	   NPR,PARAM)
C
C      RETURN
C      END
C
C-----------------------------------------------------------------------
C
C      SUBROUTINE RLQAGE1T(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,
C     *  LIMIT,
C     *  RESULT,ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST,
C     * 	NPR,PARAM)
C.......................................................................
C
C   AUTHOR :     QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
C      implicit double precision(a-h, o-z)
C      DOUBLE PRECISION A,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,
C     *  AA1,AA2,B,
C     *  BLIST,BB1,BB2,C,DABS,DEFABS,DEFAB1,DEFAB2,DMAX1,ELIST,EPMACH,
C     *  EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,F,OFLOW,
C     *  RESABS,RESULT,RLIST,UFLOW,FEXT
C      INTEGER IER,IORD,IROFF1,IROFF2,K,KEY,KEYF,LAST,LIMIT,MAXERR,NEVAL,
C     *  NRMAX
C
C      DIMENSION ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),IORD(LIMIT),
C     *  RLIST(LIMIT),FARR(N),param(npr)
C
C      EXTERNAL F,FEXT,GEXT
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
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
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH  IS THE LARGEST RELATIVE SPACING.
C           UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE.
C           OFLOW  IS THE LARGEST MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENTS
C      CALL RLMACHD(7,EPMACH)
C      CALL RLMACHD(4,UFLOW)
C      CALL RLMACHD(6,OFLOW)
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
C      NEVAL = 0
C      LAST = 0
C      RESULT = 0.0D+00
C      ABSERR = 0.0D+00
C      ALIST(1) = A
C      BLIST(1) = B
C      RLIST(1) = 0.0D+00
C      ELIST(1) = 0.0D+00
C      IORD(1) = 0
C      IER=6
C      IF ((EPSABS.LT.0..AND.EPSREL.LT.0.).OR.LIMIT.LE.1
C     1   .OR.LIMIT.GT.500)  CALL MESSGE(500,'QAGE1D',1)
C      IER = 0
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
C      KEYF = KEY
C      IF(KEY.LE.0) KEYF = 1
C      IF(KEY.GE.7) KEYF = 6
C      C = dble(FLOAT(KEYF))
C      NEVAL = 0
C      IF (KEYF.EQ.1)
C     *  CALL RLQ1K15T(F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,DEFABS,
C     *  RESABS,NPR,PARAM)
C
C      LAST = 1
C      RLIST(1) = RESULT
C      ELIST(1) = ABSERR
C      IORD(1) = 1
C
C           TEST ON ACCURACY.
C
C      ERRBND = DMAX1(EPSABS,EPSREL*DABS(RESULT))
C      IF(ABSERR.LE.5.0D+01*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
C LIMIT==1 MEANS THAT ONLY 1 SUBINTERVAL IS ALLOWED..
C      IF(LIMIT.EQ.1) IER = 1
C      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS)
C     *  .OR.ABSERR.EQ.0.0D+00) GO TO 60
C
C           INITIALIZATION
C           --------------
C
C
C      ERRMAX = ABSERR
C      MAXERR = 1
C      AREA = RESULT
C      ERRSUM = ABSERR
C      NRMAX = 1
C      IROFF1 = 0
C      IROFF2 = 0
C
C           MAIN DO-LOOP
C           ------------
C
C      DO 30 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
C
C        AA1 = ALIST(MAXERR)
C        BB1 = 5.0D-01*(ALIST(MAXERR)+BLIST(MAXERR))
C        AA2 = BB1
C        BB2 = BLIST(MAXERR)
C        IF (KEYF.EQ.1)
C     *  CALL RLQ1K15T(F,FARR,N,FEXT,GEXT,AA1,BB1,AREA1,ERROR1,
C     *          RESABS,DEFAB1,NPR,PARAM)
C        IF (KEYF.EQ.1)
C     *  CALL RLQ1K15T(F,FARR,N,FEXT,GEXT,AA2,BB2,AREA2,ERROR2,
C     *          RESABS,DEFAB2,NPR,PARAM)
C        NEVAL = NEVAL+1
C        AREA12 = AREA1+AREA2
C        ERRO12 = ERROR1+ERROR2
C        ERRSUM = ERRSUM+ERRO12-ERRMAX
C        AREA = AREA+AREA12-RLIST(MAXERR)
C        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 5
C        IF(DABS(RLIST(MAXERR)-AREA12).LE.1.0D-05*DABS(AREA12)
C     *  .AND.ERRO12.GE.9.9D-01*ERRMAX) IROFF1 = IROFF1+1
C        IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF2 = IROFF2+1
C    5   RLIST(MAXERR) = AREA1
C        RLIST(LAST) = AREA2
C        ERRBND = DMAX1(EPSABS,EPSREL*DABS(AREA))
C        IF(ERRSUM.LE.ERRBND) GO TO 8
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
C        IF(IROFF1.GE.6.OR.IROFF2.GE.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
C           SUBINTERVALS EQUALS LIMIT.
C
C        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE.
C
C        IF(DMAX1(DABS(AA1),DABS(BB2)).LE.(1.0D+00+C*1.0D+03*
C     *  EPMACH)*(DABS(AA2)+1.0D+04*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
C    8   IF(ERROR2.GT.ERROR1) GO TO 10
C        ALIST(LAST) = AA2
C        BLIST(MAXERR) = BB1
C        BLIST(LAST) = BB2
C        ELIST(MAXERR) = ERROR1
C        ELIST(LAST) = ERROR2
C        GO TO 20
C   10   ALIST(MAXERR) = AA2
C        ALIST(LAST) = AA1
C        BLIST(LAST) = BB1
C        RLIST(MAXERR) = AREA2
C        RLIST(LAST) = AREA1
C        ELIST(MAXERR) = ERROR2
C        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE QSORTD TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH THE LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
C   20   CALL RLQSORTD(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C***JUMP OUT OF DO-LOOP
C        IF(IER.NE.0.OR.ERRSUM.LE.ERRBND) GO TO 40
C   30 CONTINUE
C
C           COMPUTE FINAL RESULT.
C           ---------------------
C
C   40 RESULT = 0.0D+00
C      DO 50 K=1,LAST
C        RESULT = RESULT+RLIST(K)
C   50 CONTINUE
C      ABSERR = ERRSUM
C   60 IF(KEYF.NE.1) NEVAL = (10*KEYF+1)*(2*NEVAL+1)
C      IF(KEYF.EQ.1) NEVAL = 30*NEVAL+15
C  999 IF (IER.NE.0) CALL MESSGE(400+IER,'QAGE1D',0)
C 999 CONTINUE
C      RETURN
C      END
C
C-----------------------------------------------------------------------
C
C      SUBROUTINE RLQ1K15T
C     *  (F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,RESABS,RESASC,NPR,PARAM)
C.......................................................................
C
C   AUTHOR :     QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
C      implicit double precision(a-h, o-z)
C      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DABS,DHLGTH,DMAX1,DMIN1,
C     *  EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,OFLOW,RESABS,RESASC,
C     *  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK,FEXT
C      INTEGER J,JTW,JTWM1
C      EXTERNAL F,FEXT,GEXT
C
C      DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8),FARR(N),param(npr)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
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
C
C      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8)/
C     *     9.914553711208126D-01,   9.491079123427585D-01,
C     *     8.648644233597691D-01,   7.415311855993944D-01,
C     *     5.860872354676911D-01,   4.058451513773972D-01,
C     *     2.077849550078985D-01,   0.0D+00              /
C      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8)/
C     *     2.293532201052922D-02,   6.309209262997855D-02,
C     *     1.047900103222502D-01,   1.406532597155259D-01,
C     *     1.690047266392679D-01,   1.903505780647854D-01,
C     *     2.044329400752989D-01,   2.094821410847278D-01/
C      DATA WG(1),WG(2),WG(3),WG(4)/
C     *     1.294849661688697D-01,   2.797053914892767D-01,
C     *     3.818300505051189D-01,   4.179591836734694D-01/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C           OFLOW IS THE LARGEST MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENTS
C      CALL RLMACHD(7,EPMACH)
C      CALL RLMACHD(4,UFLOW)
C      CALL RLMACHD(6,OFLOW)
C
C      CENTR = 5.0D-01*(A+B)
C      HLGTH = 5.0D-01*(B-A)
C      DHLGTH = DABS(HLGTH)
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
C      FC = F(CENTR,FARR,N,FEXT,GEXT,NPR,PARAM)
C
C      RESG = FC*WG(4)
C      RESK = FC*WGK(8)
C      RESABS = DABS(RESK)
C      DO 10 J=1,3
C        JTW = J*2
C        ABSC = HLGTH*XGK(JTW)
C        FVAL1 = F(CENTR-ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)
C
C        FVAL2 = F(CENTR+ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)
C
C        FV1(JTW) = FVAL1
C        FV2(JTW) = FVAL2
C        FSUM = FVAL1+FVAL2
C        RESG = RESG+WG(J)*FSUM
C        RESK = RESK+WGK(JTW)*FSUM
C        RESABS = RESABS+WGK(JTW)*(DABS(FVAL1)+DABS(FVAL2))
C   10 CONTINUE
C      DO 15 J = 1,4
C        JTWM1 = J*2-1
C        ABSC = HLGTH*XGK(JTWM1)
C        FVAL1 = F(CENTR-ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)
C        FVAL2 = F(CENTR+ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)
C        FV1(JTWM1) = FVAL1
C        FV2(JTWM1) = FVAL2
C        FSUM = FVAL1+FVAL2
C        RESK = RESK+WGK(JTWM1)*FSUM
C        RESABS = RESABS+WGK(JTWM1)*(DABS(FVAL1)+DABS(FVAL2))
C   15 CONTINUE
C      RESKH = RESK*5.0D-01
C      RESASC = WGK(8)*DABS(FC-RESKH)
C      DO 20 J=1,7
C        RESASC = RESASC+WGK(J)*(DABS(FV1(J)-RESKH)+DABS(FV2(J)-RESKH))
C   20 CONTINUE
C      RESULT = RESK*HLGTH
C      RESABS = RESABS*DHLGTH
C      RESASC = RESASC*DHLGTH
C      ABSERR = DABS((RESK-RESG)*HLGTH)
C      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
C     *  ABSERR = RESASC*DMIN1(1.0D+00,(2.0D+02*ABSERR/RESASC)**1.5D+00)
C      IF(RESABS.GT.UFLOW/(5.0D+01*EPMACH)) ABSERR = DMAX1
C     *  ((EPMACH*5.0D+01)*RESABS,ABSERR)
C      RETURN
C      END
C
C      SUBROUTINE RLQSORTD(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
C      DOUBLE PRECISION ELIST,ERMAX,ERRMAX,ERRMIN
C      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
C     *  NRMAX
C      DIMENSION ELIST(LAST),IORD(LAST)
C
C           CHECK WHETHER THE LIST CONTAINS MORE THAN
C           TWO ERROR ESTIMATES.
C
C***FIRST EXECUTABLE STATEMENT
C      IF(LAST.GT.2) GO TO 10
C      IORD(1) = 1
C      IORD(2) = 2
C      GO TO 90
C
C           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
C           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
C           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
C           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
C
C   10 ERRMAX = ELIST(MAXERR)
C      IF(NRMAX.EQ.1) GO TO 30
C      IDO = NRMAX-1
C      DO 20 I = 1,IDO
C        ISUCC = IORD(NRMAX-1)
C***JUMP OUT OF DO-LOOP
C        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
C        IORD(NRMAX) = ISUCC
C        NRMAX = NRMAX-1
C   20    CONTINUE
C
C           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE
C           MAINTAINED IN DESCENDING ORDER. THIS NUMBER
C           DEPENDS ON THE NUMBER OF SUBDIVISIONS STILL ALLOWED.
C
C   30 JUPBN = LAST
C      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
C      ERRMIN = ELIST(LAST)
C
C           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
C           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
C
C      JBND = JUPBN-1
C      IBEG = NRMAX+1
C      IF(IBEG.GT.JBND) GO TO 50
C      DO 40 I=IBEG,JBND
C        ISUCC = IORD(I)
C***JUMP OUT OF DO-LOOP
C        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
C        IORD(I-1) = ISUCC
C   40 CONTINUE
C   50 IORD(JBND) = MAXERR
C      IORD(JUPBN) = LAST
C      GO TO 90
C
C           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
C
C   60 IORD(I-1) = MAXERR
C      K = JBND
C      DO 70 J=I,JBND
C        ISUCC = IORD(K)
C***JUMP OUT OF DO-LOOP
C        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
C        IORD(K+1) = ISUCC
C        K = K-1
C   70 CONTINUE
C      IORD(I) = LAST
C      GO TO 90
C   80 IORD(K+1) = LAST
C
C           SET MAXERR AND ERMAX.
C
C   90 MAXERR = IORD(NRMAX)
C      ERMAX = ELIST(MAXERR)
C      RETURN
C      END
C
C-----------------------------------------------------------------------
C
C      SUBROUTINE RLXERF(KODE,X,P)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C      implicit double precision(a-h,o-z)
C      EXTERNAL RLXEXPD
C      DATA SPI/2.506628274631D0/
C      IF (KODE.NE.1.AND.KODE.NE.2) CALL MESSGE(500,'XERF  ',1)
C      X2=-X*X/2.D0
C      P=RLXEXPD(X2)
C      IF (KODE.EQ.2) P=P/SPI
C      RETURN
C      END
C
C-----------------------------------------------------------------------
C
C      SUBROUTINE RLGAUSSD (KODE,X,P)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C
C.......................................................................
C
C      DOUBLE PRECISION   P,X,SQR1D2,CD
C      DATA               SQR1D2/.7071067811865475D0/
C
C      IF (KODE.NE.1.AND.KODE.NE.2) CALL MESSGE(500,'GAUSSD',1)
C      CALL RLCERFD(-X*SQR1D2,CD)
C      P = .5D0 * CD
C      IF (KODE.EQ.2) P=1.D0-P
C      RETURN
C      END
C
C-----------------------------------------------------------------------
C
C      SUBROUTINE RLCERFD(X,F)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C
C.......................................................................
C
C      DOUBLE PRECISION   F,X,RLXEXPD
C      DIMENSION          P(5),Q(4),P1(9),Q1(8),P2(6),Q2(5)
C      DOUBLE PRECISION   P,Q,P1,Q1,P2,Q2,XMIN,XLARGE,SQRPI,XX,
C     *                   RES,XSQ,XNUM,XDEN,XI,XBIG
C      INTEGER            ISW,I
C      EXTERNAL RLXEXPD
C                                  COEFFICIENTS FOR 0.0 .LE. Y .LT.
C                                  .477
C      DATA               P(1)/113.8641541510502D0/,
C     *                   P(2)/377.4852376853020D0/,
C     *                   P(3)/3209.377589138469D0/,
C     *                   P(4)/.1857777061846032D0/,
C     *                   P(5)/3.161123743870566D0/
C      DATA               Q(1)/244.0246379344442D0/,
C     *                   Q(2)/1282.616526077372D0/,
C     *                   Q(3)/2844.236833439171D0/,
C     *                   Q(4)/23.60129095234412D0/
C                                  COEFFICIENTS FOR .477 .LE. Y
C                                  .LE. 4.0
C      DATA               P1(1)/8.883149794388376D0/,
C     *                   P1(2)/66.11919063714163D0/,
C     *                   P1(3)/298.6351381974001D0/,
C     *                   P1(4)/881.9522212417691D0/,
C     *                   P1(5)/1712.047612634071D0/,
C     *                   P1(6)/2051.078377826071D0/,
C     *                   P1(7)/1230.339354797997D0/,
C     *                   P1(8)/2.153115354744038D-8/,
C     *                   P1(9)/.5641884969886701D0/
C      DATA               Q1(1)/117.6939508913125D0/,
C     *                   Q1(2)/537.1811018620099D0/,
C     *                   Q1(3)/1621.389574566690D0/,
C     *                   Q1(4)/3290.799235733460D0/,
C     *                   Q1(5)/4362.619090143247D0/,
C     *                   Q1(6)/3439.367674143722D0/,
C     *                   Q1(7)/1230.339354803749D0/,
C     *                   Q1(8)/15.74492611070983D0/
C                                  COEFFICIENTS FOR 4.0 .LT. Y
C      DATA               P2(1)/-3.603448999498044D-01/,
C     *                   P2(2)/-1.257817261112292D-01/,
C     *                   P2(3)/-1.608378514874228D-02/,
C     *                   P2(4)/-6.587491615298378D-04/,
C     *                   P2(5)/-1.631538713730210D-02/,
C     *                   P2(6)/-3.053266349612323D-01/
C      DATA               Q2(1)/1.872952849923460D0/,
C     *                   Q2(2)/5.279051029514284D-01/,
C     *                   Q2(3)/6.051834131244132D-02/,
C     *                   Q2(4)/2.335204976268692D-03/,
C     *                   Q2(5)/2.568520192289822D0/
C                                  CONSTANTS
C      DATA               XMIN/1.0D-10/,XLARGE/6.375D0/
C                                  CERFD(XBIG) .APPROX. DETAP
C      DATA               XBIG/13.3D0/
C      DATA               SQRPI/.5641895835477563D0/
C
C      Y=X
C      XX = Y
C      ISW = 1
C      IF (XX.GE.0.0D0) GO TO 5
C      ISW = -1
C      XX = -XX
C    5 IF (XX.LT..477D0) GO TO 10
C      IF (XX.LE.4.0D0) GO TO 30
C      IF (ISW .GT. 0) GO TO 40
C      IF (XX.LT.XLARGE) GO TO 45
C      RES = 2.0D0
C      GO TO 70
C                                  ABS(Y) .LT. .477, EVALUATE
C                                  APPROXIMATION FOR CERFD
C   10 IF (XX.LT.XMIN) GO TO 20
C      XSQ = XX*XX
C      XNUM = P(4)*XSQ+P(5)
C      XDEN = XSQ+Q(4)
C      DO 15 I = 1,3
C         XNUM = XNUM*XSQ+P(I)
C         XDEN = XDEN*XSQ+Q(I)
C   15 CONTINUE
C      RES = XX*XNUM/XDEN
C      GO TO 25
C   20 RES = XX*P(3)/Q(3)
C   25 IF (ISW.EQ.-1) RES = -RES
C      RES = 1.0D0-RES
C      GO TO 70
C                                  .477 .LE. ABS(Y) .LE. 4.0
C                                  EVALUATE APPROXIMATION FOR CERFD
C   30 XSQ = XX*XX
C      XNUM = P1(8)*XX+P1(9)
C      XDEN = XX+Q1(8)
C      DO 35 I=1,7
C         XNUM = XNUM*XX+P1(I)
C         XDEN = XDEN*XX+Q1(I)
C   35 CONTINUE
C      RES = XNUM/XDEN
C      GO TO 60
C                                  4.0 .LT. ABS(Y), EVALUATE
C                                  MINIMAX APPROXIMATION FOR CERFD
C   40 IF (XX.GT.XBIG) GO TO 65
C   45 XSQ = XX*XX
C      XI = 1.0D0/XSQ
C      XNUM= P2(5)*XI+P2(6)
C      XDEN = XI+Q2(5)
C      DO 50 I = 1,4
C         XNUM = XNUM*XI+P2(I)
C         XDEN = XDEN*XI+Q2(I)
C   50 CONTINUE
C      RES = (SQRPI+XI*XNUM/XDEN)/XX
C   60 RES = RES*RLXEXPD(-XSQ)
C      IF (ISW.EQ.-1) RES = 2.0D0-RES
C      GO TO 70
C   65 RES = 0.0D0
C   70 F = RES
C      RETURN
C      END
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION DNORM0(X)
      DOUBLE PRECISION TMP,SPI,X2,X,EXMIN
      DATA NCALL,EXMIN,SPI/0,0.D0,2.506628274631D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
      ENDIF
      X2=-X*X/2.D0
      TMP=0.D0
      IF (X2.GT.EXMIN) TMP=DEXP(X2)
      TMP=TMP/SPI
      DNORM0=TMP
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION PNORM0(Z)
      DOUBLE PRECISION Z,TMP
      CALL RLGAUSSD(1,Z,TMP)
      PNORM0=TMP
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION RLChisk(S,K0)
      DOUBLE PRECISION S,S2,ABST,TMP,K0
c     COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
c
c     K0=DBLE(XK)
      TMP=1.D0
      ABST=DABS(S)
      IF (ABST.GE.K0) GOTO 400
      S2=(S/K0)**2
      TMP=(S2*(S2-3.D0)+3.D0)*S2
  400 RLCHISK=TMP-0.5D0
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION RLPSI1N(Z,U)
      DOUBLE PRECISION Z,U
      RLPSI1N=0.D0
      IF (Z.LT.-U.OR.Z.GT.U) RETURN
      RLPSI1N=Z
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION RLPSI2N(Z,U)
      DOUBLE PRECISION Z,U
      RLPSI2N=0.D0
      IF (Z.LT.-U.OR.Z.GT.U) RETURN
      RLPSI2N=Z*Z
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION RLBETAN(u)
      double precision Alfa,u,sum,dnorm0,pnorm0,pnrm0
      external dnorm0,pnorm0
      pnrm0=pnorm0(u)
      Alfa=2.d0*pnrm0-1.d0
      sum=2.d0*(-u*dnorm0(u)+pnrm0-0.5d0)
      RLBETAN=SUM/Alfa
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION IALPHAN(Z0,U,SIGMA,IS0)
      DOUBLE PRECISION Z0,U,SIGMA,IS0,ETA,RHO,TMP,DNORM0,pnorm0,
     *       XLGMN,YLGMN
      EXTERNAL DNORM0,pnorm0
      DATA NCALL,XLGMN,YLGMN/0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
        CALL RLMACHD(4,XLGMN)
        CALL RLMACHD(5,YLGMN)
        NCALL=1
      ENDIF
      ETA=dnorm0(U)
      TMP=-YLGMN
      IF (ETA.GT.XLGMN) TMP=-dlog(ETA)
      ETA=TMP
      RHO=dnorm0(z0)
      TMP=-YLGMN
      IF (RHO.GT.XLGMN) TMP=-dlog(RHO)
      RHO=TMP
      TMP=2.d0*U*DNORM0(U)*IS0/SIGMA-(2.d0*pnorm0(U)-1.d0)
      IF (ETA.GE.RHO) TMP=TMP+1.D0
      IALPHAN=TMP
      RETURN
      END
c
      SUBROUTINE RLD1N(U,SIGMA,IT0,XtX,NP,VAL)
      DOUBLE PRECISION U,L,SIGMA,IT0(NP),XtX(NP,NP),TMP1,
     +       TMP, DNORM0,EZU,VAL(NP)
      EXTERNAL DNORM0
      L=-U
c     TMP2=(U*U-L*L)*IS0=0.D0
      TMP1=U-L
      EZU=DNORM0(U)
      DO 200 I=1,NP
      TMP=0.D0
      DO 100 J=1,NP
      TMP=TMP+ XtX(I,J)*IT0(J)
  100 CONTINUE
      VAL(I)=TMP1*TMP
      VAL(I)=EZU*(VAL(I))/SIGMA
  200 CONTINUE
      RETURN
      END
c
      SUBROUTINE RLD2N(U,SIGMA,IS0,VAL)
      DOUBLE PRECISION L,U,SIGMA,IS0,TMP2,DNORM0,EZU,U2,VAL
      EXTERNAL DNORM0
      L=-U
      U2=U*U
      TMP2=(U*U2-L*U2)*IS0
c      TMP1=0.D0
      EZU=DNORM0(U)
c      TMP=0.D0
c      DO 100 J=1,NP
c      TMP=TMP+XBAR(J)*IT0(J)
c  100 CONTINUE
c      VAL=TMP*TMP1
      VAL=EZU*TMP2/SIGMA
      RETURN
      END
c
      SUBROUTINE rlavtmlnf(X,y,n,np,ncov,u,k0,theta,sigma,invm0,
     +           invm1,avts0,avts,xbar,XtX,sa,sc1,x0,its0,its)
      implicit double precision(A-H,O-Z)
      double precision X(n,np),k0,XtX(np,np),xbar(np),y(n),theta(np),
     +       x0(np),avts0(np+1,np+1),avts(np+1,np+1),invm0(np+1,np+1),
     +       invm1(np+1,np+1),its0(np+1),its(np+1),is0,ialf,IALPHAN,
     +       sa(ncov),sc1(ncov)
      external pnorm0,IALPHAN,RLPSI1N,RLPSI2N, rlpsim2, rlchisk
c     COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
c
c     IPSI=4
c     K0=1.548D0
c
C      do 10 i=1,np+1
C      do 10 j=1,np+1
C        avts0(i,j)=0.d0
C        avts(i,j)=0.d0
C   10 continue
      do i=1,np+1
        do j=1,np+1
          avts0(i,j)=0.d0
          avts(i,j)=0.d0
        end do
      end do
c      avs0=0.d0
c      avs=0.d0
      tmp1=xbar(1)
      en2=dfloat(n)*dfloat(n-np)
      pnrm0=pnorm0(u)
      alfa=2.d0*pnrm0-1.d0
      beta=RLBETAN(u)

      do 500 k=1,n
        y0=y(k)
        z0=y0
        do 220 j=1,np
        x0(j)=X(k,j)
        z0=z0-x0(j)*theta(j)
  220   continue
        z0=z0/sigma
        tmp1=rlpsim2(z0,2,k0)
C        do 235 i=1,np
C  235   sc1(i)=tmp1*x0(i)
        do i=1,np
          sc1(i)=tmp1*x0(i)
        end do
        sc1(np+1)=rlchisk(z0,k0)
        do 240 i=1,np+1
        its0(i)=0.d0
        do 230 j=1,np+1
        its0(i)=its0(i)+invm0(i,j)*sc1(j)
  230   continue
  240   continue
        is0=its0(np+1)
c
c        z0=y0
c        do 260 j=1,np
c        x0(j)=X(k,j)
c        z0=z0-x0(j)*theta1(j)
c  260   continue
c        z0=z0/sigma1
        ialf=IALPHAN(Z0,U,SIGMA,IS0)
        tmp1=RLPSI1N(z0,u)
        call RLD1N(U,SIGMA,ITS0,XtX,NP,SA)
        call RLD2N(U,SIGMA,IS0,TMP2)
        tmp2=tmp2+RLPSI2N(z0,u)-alfa*beta-beta*ialf
C        do 265 i=1,np
C  265   sc1(i)=tmp1*x0(i)+sa(i)
        do i=1,np
          sc1(i)=tmp1*x0(i)+sa(i)
        end do
        sc1(np+1)=tmp2
        do 280 i=1,np+1
        its(i)=0.d0
        do 270 j=1,np+1
        its(i)=its(i)+invm1(i,j)*sc1(j)
  270   continue
  280   continue
c
c        do 300 i=1,np+1
c        do 300 j=1,i
c        avts0(i,j)=avts0(i,j)+its0(i)*its0(j)/en2
c        if (i.ne.j) avts0(j,i)=avts0(i,j)
c        avts(i,j)=avts(i,j)+its(i)*its(j)/en2
c        if (i.ne.j) avts(j,i)=avts(i,j)
c  300   continue
      do i=1,np+1
        do j=1,i
          avts0(i,j)=avts0(i,j)+its0(i)*its0(j)/en2
          if (i.ne.j) avts0(j,i)=avts0(i,j)
          avts(i,j)=avts(i,j)+its(i)*its(j)/en2
          if (i.ne.j) avts(j,i)=avts(i,j)
        end do
      end do
  500 continue
      return
      end
C
C********************************************************************************
C

      DOUBLE PRECISION FUNCTION RLRHOW(Z,CONST)
      implicit double precision(a-h,o-z)
      EXTERNAL RLXEXPD
c     COMMON/ZEZCOM/CONST
      RLRHOW=RLXEXPD(Z)-CONST-Z
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLRGFL2(F,CONST,Y,A,B,TOL,MAXIT,X,ITERM)
      implicit double precision(a-h,o-z)
      EXTERNAL F
c     COMMON/ZEZCOM/CONST
c     DATA TL/1.D-16/
C
C  INITIALIZE
C
      TL=DMIN1(1.D-10,0.1D0*TOL)
      ITR=1
      MESS=0
   10 FA=F(A,CONST)-Y
      FB=F(B,CONST)-Y
C
C  REGULA FALSI ITERATION
C
   20 IF (DABS(FA-FB).GT.TL) GOTO 30
      MESS=MESS+1
      IF (MESS.LE.2) THEN
        A=A/10.D0
        GOTO 10
      ENDIF
C     CALL MESSGE(401,'RGFL2 ',0)
      RETURN
   30 XN=(A*FB-B*FA)/(FB-FA)
      FN=F(XN,CONST)-Y
C
C  TEST TO SEE IF MAXIMUM NUMBER OF ITERATIONS HAS BEEN EXECUTED
C
      IF (ITR.GE.MAXIT) GOTO 60
C
C  TEST TO SEE IF ROOT HAS BEEN FOUND
C
      IF (DABS(FN).LT.TOL) GOTO 70
      IF (FA*FN.LE.0.D0) GOTO 40
      A=XN
      FA=FN
      GOTO 50
   40 B=XN
      FB=FN
C
C  INCREMENT ITERATION COUNTER
C
   50 ITR=ITR+1
      GOTO 20
C
   60 ITERM=2
      X=XN
      RETURN
   70 ITERM=1
      X=XN
      RETURN
      END
C
      SUBROUTINE RLF0W(U,TOL,MAXIT,P)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION LOW
c     COMMON/ZEZCOM/CONST
      EXTERNAL RLXEXPD,RLRHOW
      P=0.D0
      IF (U.LE.1.D0) RETURN
      P=1.D0
      IF (U.GT.16.D0) RETURN
      CONST=U
      IF (U.GT.1.5D0) THEN
        LOW=-U
        UP=-U+1.5D0
        CALL RLRGFL2(RLRHOW,CONST,0.D0,LOW,UP,TOL,MAXIT,TL,ITRM)
      ELSE
        TLO=TOL
        IF (U-1.D0.LT.1.D-3) TLO=DMIN1(TOL,1.D-8)
        LOW=-U
        UP=0.D0
        CALL RLRGFL2(RLRHOW,CONST,0.D0,LOW,UP,TLO,MAXIT,TL,ITRM)
      ENDIF
      ALOGU=DLOG(U)
      CALL RLRGFL2(RLRHOW,CONST,0.D0,ALOGU,U,TOL,MAXIT,TU,ITRM)
      XU=RLXEXPD(TU)
      CALL RLPWEIBL(1.D0,1.D0,XU,PU)
      XL=RLXEXPD(TL)
      CALL RLPWEIBL(1.D0,1.D0,XL,PL)
      P=PU-PL
      RETURN
      END
C
C      SUBROUTINE RLPWEIBL(ALPHA,SIGMA,X,P)
C      implicit double precision(a-h,o-z)
C      DATA NCALL,EXMIN,XLGMN,YLGMN/0,0.D0,0.D0,0.D0/
C      IF (NCALL.EQ.0) THEN
C        NCALL=1
C        CALL RLMACHD(3,EXMIN)
C        CALL RLMACHD(4,XLGMN)
C        CALL RLMACHD(5,YLGMN)
C      ENDIF
C     IF (ALPHA.LE.0..OR.SIGMA.LE.0.) CALL MESSGE(500,'PWEIBL',1)
C      P=0.D0
C      IF (X.LE.0.D0) RETURN
C      ALGXS=YLGMN
C      XS=X/SIGMA
C      IF (XS.GT.XLGMN) ALGXS=DLOG(XS)
C      T=ALPHA*ALGXS
C      EXPT=0.D0
C      IF (T.GT.EXMIN) EXPT=DEXP(T)
C      XXPT=0.D0
C      IF (-EXPT.GT.EXMIN) XXPT=DEXP(-EXPT)
C      P=1.D0-XXPT
C      RETURN
C      END
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION CHIS1WP(DX,WGT,N,EXU,EXV)
      implicit double precision(a-h,o-z)
      DIMENSION WGT(N)
      EXTERNAL EXU,EXV
      I= idint(WGT(1))
      B1=WGT(2)
      ANS=EXU(DX)
      IF (I.GE.4) THEN
        VV=B1
        ZV=DX/VV
        XZV=EXV(ZV)
      ENDIF
c      GOTO (10,20,30,40,50) I
      if(I == 2) then
        goto 20
      else if (I == 3) then 
        goto 30
      else if (I == 4) then
        goto 40
      else if (I == 5) then
        goto 50
      else
        goto 10
      end if
   10 CHIS1WP=(EXV(DX-B1)-1.D0)*ANS
      RETURN
   20 CHIS1WP=EXV(DX-B1)*ANS
      RETURN
   30 CHIS1WP=DX*(EXV(DX)-1.D0)*ANS
      RETURN
   40 CHIS1WP=ZV*(XZV-1.D0)*ANS
      RETURN
   50 CHIS1WP=-ZV*(XZV-1.D0+XZV*ZV)*ANS/VV
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE RLINTMW(IWGT,TL,TU,B1,TIL,SUM)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION rlezez,ERRSTD,WORK,CHIS1WP,LO,HI,
     +       TL,TU,TIL,SUM
      DIMENSION WGT(2),IWORK(80),WORK(320)
      EXTERNAL rlezez,CHIS1WP,RLXEXPD
C
C     INITIALISATION
C
      LIMIT=80
      KEY=1
      LO=TL
      HI=TU
      WGT(1)=DFLOAT(IWGT)
      WGT(2)=B1
      CALL RLINTGRT(CHIS1WP,WGT,2,rlezez,RLXEXPD,LO,HI,TIL,TIL,
     1            KEY,LIMIT,SUM,ERRSTD,NEVAL,IER,WORK,IWORK,2,WGT)
      RETURN
      END

C======================================================================

      DOUBLE PRECISION FUNCTION rlezez(Z)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION Z,EXMIN,TMP,VAL
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
      ENDIF
      TMP=Z
      IF (Z.GE.EXMIN) TMP=Z-DEXP(Z)
      VAL=0.D0
      IF (TMP.GT.EXMIN) VAL=DEXP(TMP)
      rlezez=VAL
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION rlpezez(Z)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION Z,EXMIN,TMP,VAL
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
      ENDIF
      TMP=0.D0
      IF (Z.GT.EXMIN) TMP=-DEXP(Z)
      VAL=0.D0
      IF (TMP.GT.EXMIN) VAL=DEXP(TMP)
      rlpezez=1.D0-VAL
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION RLPSI1W(Z,L,U)
      DOUBLE PRECISION Z,L,U,EXMIN,TMP
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
      ENDIF
      RLPSI1W=0.D0
      IF (Z.LT.L.OR.Z.GT.U) RETURN
      TMP=-1.D0
      IF (Z.GT.EXMIN) TMP=DEXP(Z)-1.D0
      RLPSI1W=TMP
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION RLPSI2W(Z,L,U)
      DOUBLE PRECISION Z,L,U,EXMIN,TMP
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
      ENDIF
      RLPSI2W=0.D0
      IF (Z.LT.L.OR.Z.GT.U) RETURN
      TMP=-Z
      IF (Z.GT.EXMIN) TMP=Z*(DEXP(Z)-1.D0)
      RLPSI2W=TMP
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION RLBetaw(l,u)
      implicit double precision(a-h,o-z)
      double precision Alfa,l,u,tild,sum,rlpezez
      external rlpezez
      Alfa=rlpezez(u)-rlpezez(l)
      tild=1.D-4
      CALL RLINTMW(3,L,U,0.D0,TILD,SUM)
      RLbetaw=SUM/Alfa
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION RLIALFAW(Z0,L,U,SIGMA,IS0)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION Z0,L,U,SIGMA,IS0,ETA,RHO,TMP,RLEZEZ,RLPEZEZ,EXMIN
      EXTERNAL RLEZEZ,RLPEZEZ
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
      ENDIF
      ETA=DEXP(U)-U
      RHO=-Z0
      IF (Z0.GT.EXMIN) RHO=DEXP(Z0)-Z0
      TMP=(U*RLEZEZ(U)-L*RLEZEZ(L))*IS0/SIGMA-(RLPEZEZ(U)-RLPEZEZ(L))
      IF (ETA.GE.RHO) TMP=TMP+1.D0
      RLIALFAW=TMP
      RETURN
      END

c
      SUBROUTINE RLD1W(L,U,SIGMA,IT0,IS0,XtX,XBAR,NP,VAL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION L,U,SIGMA,IT0(NP),IS0,XtX(NP,NP),TMP,TMP1,TMP2,
     +       RLEZEZ,EZU,VAL(NP),EXMIN,DXPL,XBAR(NP)
      EXTERNAL RLEZEZ
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
      ENDIF
      DXPL=0.D0
      IF (L.GT.EXMIN) DXPL=DEXP(L)
      TMP1=(DEXP(U)-DXPL)
      TMP2=(U*DEXP(U)-U-L*DXPL+L)*IS0
      EZU=RLEZEZ(U)
      DO 200 I=1,NP
      TMP=0.D0
      DO 100 J=1,NP
      TMP=TMP+XtX(I,J)*IT0(J)
  100 CONTINUE
      VAL(I)=TMP1*TMP+TMP2*XBAR(I)
      VAL(I)=EZU*(VAL(I))/SIGMA
  200 CONTINUE
      RETURN
      END
c
      SUBROUTINE RLD2W(L,U,SIGMA,IT0,IS0,XBAR,NP,VAL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION L,U,SIGMA,IT0(NP),IS0,XBAR(NP),TMP,TMP1,TMP2,
     +       RLEZEZ,EZU,VAL,EXMIN,DXPL
      EXTERNAL RLEZEZ
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
      ENDIF
      DXPL=0.D0
      IF (L.GT.EXMIN) DXPL=DEXP(L)
      TMP2=(U*U*(DEXP(U)-1.D0)-L*L*(DXPL-1.D0))*IS0
      TMP=U*(DEXP(U)-1.D0)-L*(DXPL-1.D0)
      EZU=RLEZEZ(U)
      TMP1=0.D0
      DO 100 J=1,NP
      TMP1=TMP1+XBAR(J)*IT0(J)
  100 CONTINUE
      VAL=TMP*TMP1
      VAL=EZU*(VAL+TMP2)/SIGMA
      RETURN
      END
c
      SUBROUTINE rlavtmlwf(X,y,n,np,ncov,l,u,xk,theta,sigma,invm0,
     +           invm1,avts0,avts,xbar,XtX,sa,sc1,x0,its0,its)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision X(n,np),l,u,xbar(np),y(n),theta(np),XtX(np,np),
     +       invm0(np+1,np+1),invm1(np+1,np+1),x0(np),is0,sa(ncov),
     +       sc1(ncov),avts0(np+1,np+1),avts(np+1,np+1),ialf,
     +       its0(np+1),its(np+1)
      external rlpezez,RLIALFAW,RLPsi1w,RLPsi2w,rlchisk
c     COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
c
c
c      do 10 i=1,np+1
c      do 10 j=1,np+1
c        avts0(i,j)=0.d0
c        avts(i,j)=0.d0
c   10 continue
      do i=1,np+1
        do j=1,np+1
          avts0(i,j)=0.d0
          avts(i,j)=0.d0
        end do
      end do
      en2=dfloat(n)*dfloat(n-np)

      alfa=rlpezez(u)-rlpezez(l)
      beta=RLBetaw(l,u)
      do 500 k=1,n
        y0=y(k)
        z0=y0
        do 220 j=1,np
           x0(j)=X(k,j)
           z0=z0-x0(j)*theta(j)
  220   continue
        z0=z0/sigma
        tmp1=rlpsim2(z0,2,xk)
c        do 235 i=1,np
c  235   sc1(i)=tmp1*x0(i)
        do i=1,np
          sc1(i)=tmp1*x0(i)
        end do
        sc1(np+1)=rlchisk(z0,xk)
        do 240 i=1,np+1
        its0(i)=0.d0
        do 230 j=1,np+1
        its0(i)=its0(i)+invm0(i,j)*sc1(j)
  230   continue
  240   continue
        is0=its0(np+1)
        its0(1)=its0(1)+0.1352D0*is0
c
c        z0=y0
c        do 260 j=1,np
c        x0(j)=X(k,j)
c        z0=z0-x0(j)*theta1(j)
c  260   continue
c        z0=z0/sigma1
        ialf=RLIALFAW(Z0,L,U,SIGMA,IS0)
        tmp1=RLPsi1w(z0,l,u)
        call RLD1W(L,U,SIGMA,ITS0,IS0,XtX,XBAR,NP,SA)
        call RLD2W(L,U,SIGMA,ITS0,IS0,XBAR,NP,TMP2)
        tmp2=tmp2+RLPsi2w(z0,l,u)-alfa*beta-beta*ialf
c        do 265 i=1,np
c  265   sc1(i)=tmp1*x0(i)+sa(i)
        do i=1,np
          sc1(i)=tmp1*x0(i)+sa(i)
        end do
        sc1(np+1)=tmp2
        do 280 i=1,np+1
        its(i)=0.d0
        do 270 j=1,np+1
        its(i)=its(i)+invm1(i,j)*sc1(j)
  270   continue
  280   continue
c
c        do 300 i=1,np+1
c        do 300 j=1,i
c        avts0(i,j)=avts0(i,j)+its0(i)*its0(j)/en2
c        if (i.ne.j) avts0(j,i)=avts0(i,j)
c        avts(i,j)=avts(i,j)+its(i)*its(j)/en2
c        if (i.ne.j) avts(j,i)=avts(i,j)
c  300   continue
      do i=1,np+1
        do j=1,i
          avts0(i,j)=avts0(i,j)+its0(i)*its0(j)/en2
          if (i.ne.j) avts0(j,i)=avts0(i,j)
          avts(i,j)=avts(i,j)+its(i)*its(j)/en2
          if (i.ne.j) avts(j,i)=avts(i,j)
        end do
      end do
  500 continue
      return
      end


