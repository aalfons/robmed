C
C       file: glmrob.f
C
C       Matias - Jeff
C
C
      SUBROUTINE RLGINTAC(X,Y,NI,OI,MDX,MDT,N,NP,NCOV,ICASE,MAXTT,
     +     MAXTA,TAU,TOLT,TOLA,B,C,NITT,NITA,SIGMA,A,THETA,CI,DIST,
     +     RW1,RW2,IW1,DW1,IPS,XK)
C.......................................................................
C
C  INITIAL VALUES FOR THE LOGISTIC REGRESSION
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION X(MDX,NP),Y(N),THETA(MDT),OI(N),CI(N),DIST(N),
     +          RW1(5*NCOV+3*N),RW2(MDX,NP)
      DOUBLE PRECISION A(NCOV),DW1(2*NCOV+NP+N)
      INTEGER NI(N),IW1(NP)

      IWG=NCOV+1
      ISC=IWG+N
      ISF=ISC+NCOV
      ISG=ISF+NCOV
      ISH=ISG+NCOV
      ISY=ISH+NCOV
      ISW=ISY+N
      IST=NCOV+1
      ISD=IST+NCOV
      ISU=ISD+NP
      CALL RLGITAC2(X,Y,NI,OI,MDX,MDT,N,NP,NCOV,ICASE,MAXTT,MAXTA,
     +     TAU,TOLT,TOLA,B,C,NITT,NITA,SIGMA,A,THETA,CI,DIST,RW1,
     +     RW1(IWG),RW1(ISC),RW1(ISF),RW1(ISG),RW1(ISH),RW1(ISY),
     +     RW1(ISW),RW2,IW1,DW1,DW1(IST),DW1(ISD),DW1(ISU),IPS,XK)
      RETURN
      END

C-----------------------------------------------------------------------
C
      SUBROUTINE RLGITAC2(X,Y,NI,OI,MDX,MDT,N,NP,NCOV,ICASE,MAXTT,
     +     MAXTA,TAU,TOLT,TOLA,B,C,NITT,NITA,SIGMA,A,THETA,CI,
     +     DIST,COV,WGT,SC,SF,SG,SH,SY,SW,SX,SP,SA,ST,SD,SU,IPS,XK)
C.......................................................................
C
C  INITIAL VALUES FOR THE LOGISTIC REGRESSION
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),COV(NCOV),CI(N),THETA(MDT),WGT(N),
     +     DIST(N),OI(N),SC(NCOV),SF(NCOV),SG(NCOV),SH(NCOV),
     +     SX(MDX,NP),SY(N),SW(N)
      DIMENSION A(NCOV),SA(NCOV),SD(NP),ST(NCOV),SU(N)
      INTEGER NI(N),SP(NP)
      DATA ZERO,ONE/0.D0,1.D0/
C
C these are for the function ucv_glmrob (ex ucv) in rlwyfalg
C
      IUCV=1
      A2=ZERO
      B2=B*B

      NFIRST=N
C      IERR=0
      INIT=1
C convergence criteria for wyfalg
      ICNV=2
      ITYPW=1
C
C INITIAL VALUE FOR A
C
      CALL RLWEDVBI(X,NP,NCOV,MDX,1,INIT,NFIRST,A,SY)

C
C COMPUTES A (FIXED POINT ALGORITHM)
C       (I REMOVED THE ARGUMENT "NITMON")
C
      CALL RLWFAGBI(X,A,SW,N,NP,0,NCOV,MDX,TAU,MAXTA,ICNV,ITYPW,
     +     0,TOLA,NITA,WGT,SU,SA,ST,SD,SD,IUCV,A2,B2)

      IPS=3
      XK=B

      DO 5 I=1,N
         S=WGT(I)
         IF (S.LE.1.D-3) S=1.D-3
         WGT(I)=RLPSIM2(S,IPS,XK)/S
 5    CONTINUE
      DO 30 L=1,N
         IF (ICASE.EQ.3) GOTO 10
         IF (ICASE.EQ.1) THEN
            ENP1=2.D0
         ELSE
            ENP1=DBLE(NI(L))+ONE
         ENDIF
         S=(Y(L)+0.5D0)/ENP1
         SY(L)=DLOG(S/(ONE-S))-OI(L)
         GOTO 15
 10      S=Y(L)
         IF (S.LE.ZERO) S=0.5D0
         SY(L)=DLOG(S)-OI(L)
 15      SU(L)=SY(L)
         DO 20 K=1,NP
            SX(L,K)=X(L,K)
 20      CONTINUE
 30   CONTINUE

C
C SOLVE L1 PROBLEM
C
C INITIAL BETA -> BET0
      CALL RLBET0BI(WGT,N,2,1,TOLT,BET0)
      CALL RLLARSBI(SX,SY,N,NP,MDX,MDT,TAU,NIT,K,KODE,
     +     SIG0,THETA,SW,CI,SF,SG,SH,BET0)

c
c     sig0 == 0 !!!!
c
c     maybe this is making sigma=0 in RLrwagM2
c

      IF (SIG0 .LE. TAU) SIG0=ONE
      CPSI=C

      IA=1
      IAINV=0
      F1=ONE/DBLE(N)
      F0=ZERO

C
C     SY <- double(SU)
C
      DO 40 L=1,N
         SY(L)=SU(L)
 40   CONTINUE

C
C     SU <- double(Y)
C
      DO 50 L=1,N
         SU(L)=Y(L)
 50   CONTINUE
      ITYP=2
      CALL RLKEDHBI(WGT,N,CPSI,ITYP,CI,SW)
      CALL RLKTASBI(X,CI,SW,N,NP,MDX,MDX,NCOV,TAU,IA,F1,F0,IAINV,
     +     SC,SF,SG,SH,COV,SX)

      ISIGMA=2
      PSP0=ONE
      GAMT=ONE
      ICNV=1
      MAXIS=1

C
C     should the first bet0 be beta in the following?
C
      CALL RLRWAGM2(X,SY,THETA,WGT,COV,PSP0,SIG0,N,NP,MDX,
     +     NCOV,TOLT,GAMT,TAU,ITYP,ISIGMA,ICNV,MAXTT,MAXIS,
     +     NITT,SIGMA,Y,SC,CI,SF,SG,SH,SP,SW,SX,IPS,XK,BET0,BET0)

C     IF (NITT.EQ.MAXTT) IERR=1

      DO 60 L=1,N
         CI(L)=ZERO
 60   CONTINUE
      DO 70 L=1,N
         Y(L)=SU(L)
 70   CONTINUE
      DO 90 I=1,N
         DO 80 J=1,NP
            SD(J)=X(I,J)
 80      CONTINUE
         CALL RLMLYDBI(A,SD,NP,NCOV,NP,1)
         CALL RLNRM2BI(SD,NP,1,NP,Z)
         DIST(I)=Z
 90   CONTINUE
      RETURN
      END

C-----------------------------------------------------------------------
      SUBROUTINE RLMACHD(I,X)
C.......................................................................
C
C
C                   DOUBLE PRECISION VERSION
C                   ************************
C
C  MACHINE PARAMETERS : TO  ALTER  THIS  SUBROUTINE  FOR  A PARTICULAR
C  ++++++++++++++++++   ENVIRONMENT, THE DESIRED SET OF DATA STATEMENT
C  SHOULD BE ACTIVATED BY REMOVING THE "C" FROM COLUMN ONE  AND ADDING
C  THE "C" FOR THE TWO LINES AFTER "... VAX FORTRAN (V5) compiler".
C
C  RADIX IS ALWAYS EQUAL TO 2.D0
C  PREC CAN BE FOUND BY CALLING THE ROBETH SUBROUTINE "PRECD"
C  EPMACH IS APPROXIMATELY EQUAL TO THE EXPONENT PART OF PREC
C  EXMIN, XLGMN, YLGMN AND XBIG CAN BE FOUND BY TRIAL AND ERROR
C
      DOUBLE PRECISION X,RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C  for VAX FORTRAN (V5) compiler (VMS)
C       DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C     *  /2.D0,1.38778D-17,-88.722D0,2.939D-39,-88.7227D0,1.7D38,1.D-17/
C  for IBM-PC F77L compiler
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *  /2.D0,5.422D-20,-708.D0,1.D-307,-706.591D0,1.D308,1.D-19/
C  for IBM-PC MICROSOFT FORTRAN compiler
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *  /2.D0,5.47522D-18,-745.133D0,0.9D-48,-110.629D0,1.D308,1.D-17/
C  for WATCOM F77 Compiler (32 bits version)
      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
     +     /2.,0.1121D-15,-709.782D0,9.74D-290,-718.433D0,1.797D308,
     +     1.0D-17/
C    */2.,0.1121D-15,-709.782D0,0.974D-312,-718.433D0,1.797D308,1.0D-17/
C  for ULTRIX DEC FORTRAN-77 compiler
C     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,1.12133D-16,-744.44D0,2.226D-308,-708.396D0,1.797D308,1.0D-17/
C  for SUN FORTRAN compiler
C     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,1.12133D-16,-745.13D0,0.494D-323,-744.44D0,1.797D308,1.0D-17/
C  for SILICON GRAPHICS MIPS FORTRAN 77 Compiler
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,1.12133D-16,-744.04D0,0.758D-323,-743.75D0,1.797D308,1.0D-17/
C  for HP-UX FORTRAN 77 Compiler
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,1.12133D-16,-708.396D0,0.1D-308,-709.09D0,1.797D308,1.0D-17/
C  for DEC-ALPHA FORTRAN Compiler (OpenVMS)
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,1.1102D-16,-709.782D0,0.1057D-45,-105.863D0,0.898D307,1.0D-17/
C
      IF (I.EQ.1) X=RADIX
      IF (I.EQ.2) X=PREC
      IF (I.EQ.3) X=EXMIN
      IF (I.EQ.4) X=XLGMN
      IF (I.EQ.5) X=YLGMN
      IF (I.EQ.6) X=XBIG
      IF (I.EQ.7) X=EPMACH
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGFEDCA(VTHETA,CI,WA,NI,OI,N,ICASE,D,E)
C.......................................................................
C
C  D AND E MATRICES FOR THE COV. MATRIX OF THE COEFF. ESTIMATES
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VTHETA(N),CI(N),OI(N),WA(N),NI(N),D(N),E(N)
      DATA ZERO,ONE/0.D0,1.D0/
      DATA PREC/0.D0/

      IF (PREC.EQ.ZERO) CALL RLMACHD(2,PREC)
      DO 500 I=1,N
         GI=VTHETA(I)+OI(I)
         CC=CI(I)
         AA=WA(I)
         PI=RLGFUN(ICASE,1,GI)
         IF (ICASE.EQ.1) THEN
C==>      LOGISTIC BERNOULLI
            TEMP=-PI-CC
            T0=DMIN1(AA,DABS(TEMP))
            IF (TEMP.LT.ZERO) T0=-T0
            TEMP=ONE-PI-CC
            A = dble(1)
            T1=DMIN1(A,ABS(TEMP))
            IF (TEMP.LT.0.D0) T1=-T1
            D(I)=T1*(ONE-PI)*PI-T0*PI*(ONE-PI)
            E(I)=T1**2*PI+T0**2*(ONE-PI)
            GOTO 500
         ELSEIF (ICASE.EQ.2) THEN
C==>      LOGISTIC BINOMIAL
            DI=ZERO
            EI=ZERO
            MI=NI(I)
            DNI=DBLE(MI)
            MED=1+MI/2
            DO 200 J=0,MI
               CALL RLPROBIN(J,MI,PI,PJ)
               FJME=DBLE(J)-DNI*PI
               TMP=FJME-CC
               TT=DMIN1(AA,DABS(TMP))
               IF (TMP.LT.ZERO) TT=-TT
               TPJ=TT*PJ
               TMP=TPJ*FJME
               DI=DI+TMP
               TE=TPJ*TT
               EI=EI+TE
               IF (J.LE.MED) GOTO 200
               IF (DABS(TT).GT.ZERO .AND. DABS(TMP).LE.PREC
     +              .AND.TE.LE.PREC) GOTO 350
 200        CONTINUE

         ELSEIF (ICASE.EQ.3) THEN
C==>      LOGISTIC POISSON
            DI=ZERO
            EI=ZERO
            MI=INT(100.*PI)
            DO 300 J=0,MI
               CALL RLPRPOIS(PI,J,PJ)
               FJME=DBLE(J)-PI
               TMP=FJME-CC
               TT=DMIN1(AA,DABS(TMP))
               IF (TMP.LT.ZERO) TT=-TT
               TPJ=TT*PJ
               TMP=TPJ*FJME
               DI=DI+TMP
               TE=TPJ*TT
               EI=EI+TE
               IF (DABS(TT).GT.ZERO .AND. DABS(TMP).LE.PREC
     +              .AND.TE.LE.PREC) GOTO 350
 300        CONTINUE
         ENDIF
 350     D(I)=DI
         E(I)=EI
 500  CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLPROBIN(K,N,P,PK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION LPL
      DATA NCALL,KL,LPL,EMIN,SML,ALSML/0,0,0.D0,0.D0,0.D0,0.D0/
      DATA ALP,ALQ/0.D0,0.D0/

      PK=0.D0
      IF (NCALL.EQ.0) THEN
         CALL RLMACHD(3,EMIN)
         CALL RLMACHD(4,SML)
         CALL RLMACHD(5,ALSML)
         NCALL=1
      ENDIF
      IF (P .NE. 0.D0) GO TO 15
      PK = 1.D0
      IF (K .NE. 0) PK=0.D0
      GOTO 900
 15   IF (P .NE. 1.) GO TO 20
      PK = 1.D0
      IF (K .NE. N) PK=0.D0
      GO TO 900
 20   IF (K.EQ.0.OR.KL+1.NE.K) THEN
         Q1=1.D0-P
         ALQ=ALSML
         IF (Q1.GT.SML) ALQ=DLOG(Q1)
         ALP=ALSML
         IF (P.GT.SML) ALP=DLOG(P)
      ENDIF
      IF (K.EQ.0) THEN
         PK=0.D0
         LPL=DBLE(N)*ALQ
         IF (LPL.GT.EMIN) PK=DEXP(LPL)
      ELSEIF (KL+1.NE.K.OR.LPL.LE.ALSML) THEN
         CALL RLBINPRD(K,N,P,S1,PK)
         GOTO 900
      ELSE
         LPL=LPL+DLOG(DBLE(N-K+1))+ALP-DLOG(DBLE(K))-ALQ
         PK=0.D0
         IF (LPL.GT.EMIN) PK=DEXP(LPL)
      ENDIF
      GOTO 950
 900  LPL=ALSML
      IF (PK.GT.SML) LPL=DLOG(PK)
 950  KL=K
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLPRPOIS(E,K,PK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION LPL,LE
      DATA NCALL,KL,LPL,LE,ESML,XLMN,YLMN/0,0,0.D0,0.D0,0.D0,0.D0,0.D0/

      PK=0.D0
      IF (NCALL.EQ.0) THEN
         CALL RLMACHD(3,ESML)
         CALL RLMACHD(4,XLMN)
         CALL RLMACHD(5,YLMN)
         NCALL=1
      ENDIF
      IF (K.GT.1100000) THEN
C       For E.LE.1E6 and K.GT.1100000 the probability
C       PK is smaller than 1E-2000.
         LPL = YLMN
         PK=0.D0
         GOTO 950
      ELSEIF (E.LT.DSQRT(XLMN)) THEN
C       This case is treated here in order to reduce underflow
C       problems
         PK = 0.D0
         IF (K.EQ.0) PK = 1.D0
         IF (K.EQ.1) PK = E
         GOTO 900
      ENDIF
      IF (K.EQ.0.OR.KL+1.NE.K) THEN
         LE=YLMN
         IF (E.GT.XLMN) LE=DLOG(E)
      ENDIF
      IF (K.EQ.0) THEN
         LPL = -E
      ELSEIF (KL+1.NE.K.OR.LPL.LE.YLMN) THEN
         CALL RLPOISSN(E,K,S1,PK)
         GOTO 900
      ELSE
         LPL=LPL+LE-DLOG(DBLE(K))
      ENDIF
      PK=0.D0
      IF (LPL.GT.ESML) PK = DEXP(LPL)
      GOTO 950
 900  LPL=YLMN
      IF (PK.GT.XLMN) LPL=DLOG(PK)
 950  KL=K
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLPOISSN(LAMBDA,K,PS,PK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION LAMBDA,IAX,JAX,LPK,MEDIAN

      PS = 0.D0
      PK = 0.D0
C      PGT = 0.D0

      CALL RLMACHD(3,EXMIN)
      CALL RLMACHD(4,XLGMN)
C
C     Returns the lower tail and point probabilities
C     associated with a Poisson distribution.
C
C     Let Z denote a random variable having a Poisson distribution with
C     parameters N and P. The routine computes for given LAMBDA and K:
C
      IF (K.GT.1100000) THEN
C       For LAMBDA.LE.1E6 and K.GT.1100000 the probability
C       PK is smaller than 1E-2000.
         PS = 1.D0
         PK = 0.D0
         RETURN
      ELSEIF (LAMBDA.LT.DSQRT(XLGMN)) THEN
C       This case is treated here in order to reduce underflow
C       problems
         PS = 1.D0
         PK = 0.D0
         IF (K.EQ.0) PK = 1.D0
         IF (K.EQ.1) PK = LAMBDA
         RETURN
      ENDIF
      A = DBLE(K+1)
      I2A=2*(K+1)
      X=LAMBDA
C
C     Computation of PK
C
      IF (A.EQ.1.D0) THEN
         LPK = -X
      ELSE
         CALL RLNLGMBI(I2A,GL)
         LPK=-X+(A-1.D0)*DLOG(X)-GL
      ENDIF
      PK = RLXEXPD(LPK)
      MEDIAN = A - 0.33D0
      IF (X.LE.MEDIAN) THEN
C
C       Compute PS = 1 - PK*IAX
C
         IF (LPK.GE.EXMIN) THEN
            CALL RLINTGM0(X,A,IAX)
            PS = 1.D0 - PK*IAX
         ELSE IF (2*X.LE.A) THEN
            PS = 1.D0
         ELSE
            Q = X/A
            ARG = LPK + DLOG(Q/(1.D0-Q))
            IF (ARG.LE.EXMIN) THEN
               PS = 1.D0
            ELSE
               CALL RLINTGM0(X,A,IAX)
               ARG = LPK + DLOG(IAX)
               PS = 1.D0 - RLXEXPD(ARG)
            ENDIF
         ENDIF
      ELSE
C
C       Compute PS = PK*JAX
C
         IF (LPK.GE.EXMIN) THEN
            CALL RLINTGM1(X,A,JAX)
            PS = PK*JAX
         ELSE
            Q = (A-1.D0)/X
            ARG = LPK - DLOG(1.D0-Q)
            IF (ARG.LE.EXMIN) THEN
               PS = 0.D0
            ELSE
               CALL RLINTGM1(X,A,JAX)
               ARG = LPK + DLOG(JAX)
               PS = RLXEXPD(ARG)
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLINTGM0(X,A,IAX)
C
C     Computes I(A,X) = Integral from t=0 to t=X of
C                        (t/X)**(A-1) * EXP(X-t)
C     0 .LT. X .LE. A .LE. 1E6
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION IAX
      DATA EPS/0.D-6/
C
C     Determine NSTEP, the number of steps in the backward recursion:
C
      NSTEP = 0
      FAC = 1.D0
      B = A
 20   CONTINUE
      FAC = FAC*X/B
      B = B + 1.D0
      NSTEP = NSTEP + 1
      IF (FAC.GT.EPS) GO TO 20
C
C     Backward recursion with NSTEP iteration steps:
C
      IAX = 0.D0
      DO 40 I = 1, NSTEP
         B = B - 1.D0
         IAX = (1.D0+IAX)*X/B
 40   CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLINTGM1(X,A,JAX)
C
C     Computes J(A,X) = Integral from t=X to t=infinity of
C                        (t/X)**(A-1) * EXP(X-t)
C     0 .LT. A .LE. X .LE. 1E6 (A VALUE IS INTEGER)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION JAX
      DATA EPS/0.5D-6/
C
C     Determine NSTEP, the number of steps in the forward recursion:
C
      NSTEP = 0
      FAC = 1.D0
      B = A
 20   CONTINUE
      B = B - 1.D0
      FAC = FAC*B/X
      NSTEP = NSTEP + 1
      IF (FAC.GT.EPS) GO TO 20
C
C     Forward recursion with NSTEP iteration steps:
C
      JAX = 1.D0
      DO 40 I = 2, NSTEP
         B = B + 1.D0
         JAX = 1.D0 + JAX*B/X
   40 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE RLBINPRD(K,N,P,PS,PK)
C.......................................................................
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PK = 0.D0
      PS = 0.D0
      CALL RLMACHD(4,SML)
      CALL RLMACHD(5,ALSML)
      IF (P .NE. 0.D0) GO TO 15
      PS = 1.D0
      IF (K .NE. 0) GO TO 900
      PK = 1.D0
      GO TO 900
 15   IF (P .NE. 1.D0) GO TO 20
      IF (K .NE. N) GO TO 900
      PK = 1.D0
      PS = 1.D0
      GO TO 900
 20   P1 = P
      Q1 = 1.D0-P
      K1 = K
      XN = N
      XX = XN*P
      IF (K .LE. XX) GO TO 25
      P1 = Q1
      Q1 = P
      K1 = N-K
 25   ALQN = XN*DLOG(Q1)
      ICNT = idint(ALQN/ALSML)
      ALQN = ALQN-ICNT*ALSML
      PK = RLXEXPD(ALQN)
      IF (K1 .EQ. 0) GO TO 35
      QP = P1/Q1
      XJ = 0.D0
      XN = XN+1.D0
      DO 30 J = 1,K1
         IF (ICNT .EQ. 0) PS = PS+PK
         XJ = XJ+1.D0
         PK = PK*(QP*(XN-XJ))
         IF (PK .LT. XJ) GO TO 30
         PK = PK*SML
         ICNT = ICNT-1
         PK = PK/XJ
 30   CONTINUE
 35   IF (ICNT .NE. 0) PK = 0.D0
      IF (K .GT. XX) GO TO 40
      PS = PS+PK
      GO TO 900
 40   PS = 1.D0-PS
 900  RETURN
      END
C
C-----------------------------------------------------------------------
C
      double precision FUNCTION RLGFUN(ICASE,NI,GI)
C.......................................................................
C
C  G-FUNCTION FOR LOGISTIC REGRESSION
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DATA NCALL,DMIN,DMAX,XBIG/0,0.D0,0.D0,0.D0/

      IF (NCALL.EQ.1) GOTO 10
      CALL RLMACHD(3,DMIN)
      CALL RLMACHD(6,XBIG)
      XBIG=XBIG/10.D0
      DMAX=DLOG(XBIG)
      NCALL=1
 10   IF (ICASE.LE.2) THEN
C==>    LOGISTIC BERNOUILLI OR BINOMIAL
         IF (GI.LE.DMIN) THEN
            RLGFUN=0.D0
         ELSEIF (GI.GE.DMAX) THEN
            RLGFUN=DBLE(NI)
         ELSE
            EXGI=DEXP(GI)
            RLGFUN=DBLE(NI)*EXGI/(1.D0+EXGI)
         ENDIF
      ELSE
C==>    LOGISTIC POISSON (ICASE=3)
         RLGFUN=RLXEXPD(GI)
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      double precision FUNCTION RLXEXPD(X)
C.......................................................................
C
C     Extended exp() function
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DATA NCALL,DMIN,DMAX,XBIG/0,0.D0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
         CALL RLMACHD(3,DMIN)
         CALL RLMACHD(6,XBIG)
         XBIG=XBIG/10.D0
         DMAX=DLOG(XBIG)
         NCALL=1
      ENDIF
      IF (X.LE.DMIN) THEN
         RLXEXPD=0.D0
      ELSEIF (X.GE.DMAX) THEN
         RLXEXPD=XBIG
      ELSE
         RLXEXPD=DEXP(X)
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGYMAIN(X,Y,NI,COV,A,THETA,OI,MDX,N,NP,NCOV,B,GAM,TAU,
     +     ICASE,IUGL,IOPT,IALG,ICNVT,ICNVA,MAXIT,MAXTT,MAXTA,MAXTC,
     +     TOL,TOLT,TOLA,TOLC,NIT,CI,WA,VTHETA,
     +     DELTA,GRAD,HESSNV,RW1,RW2,IW1,DW1, TRACE)
C.......................................................................
C
C  MAIN ALGORITHM FOR THE LOGISTIC REGRESSION
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),COV(NCOV),THETA(NP),CI(N),VTHETA(N),
     +     WA(N),DELTA(NP),GRAD(NP),HESSNV(NCOV),RW2(MDX,NP),
     +     RW1(5*NCOV+3*N),OI(N),A(NCOV),DW1(2*NCOV+NP+N)
      INTEGER IW1(NP),NI(N), TRACE
      DATA ZMIN/1D-6/

      IF1=N+1
      IF2=IF1+N
      ISC=IF2+N
      ISE=ISC+NCOV
      ISF=ISE+NCOV
      ISG=ISF+NP
      ISH=ISG+NP
      IST=NCOV+1
      ISD=IST+NCOV
      ISU=ISD+NP
      CALL RLGMAIN2(X,Y,NI,COV,A,THETA,OI,MDX,N,NP,NCOV,B,GAM,TAU,
     +     ICASE,IUGL,IOPT,IALG,ICNVT,ICNVA,MAXIT,MAXTT,MAXTA,MAXTC,
C
C     NITMNT,NITMNA,
C
     +     TOL,TOLT,TOLA,TOLC,ZMIN,NIT,CI,WA,VTHETA,
     +     DELTA,GRAD,HESSNV,RW1,RW1(IF1),RW1(IF2),RW1(ISC),RW1(ISE),
     +     RW1(ISF),RW1(ISG),RW1(ISH),RW2,IW1,DW1,DW1(IST),DW1(ISD),
     +     DW1(ISU), TRACE)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGMAIN2(X,Y,NI,COV,A,THETA,OI,MDX,N,NP,NCOV,B,
     1     GAM,
     1     TAU,
     1     ICASE,IUGL,IOPT,IALG,ICNVT,ICNVA,MAXIT,MAXTT,MAXTA,
     1     MAXTC,
     1     TOL,TOLT,TOLA,TOLC,ZMIN,NIT,CI,WA,VTHETA,
     1     DELTA,GRAD,HESSNV,F0,F1,F2,SC,SE,SF,SG,SH,SX,SP,SA,ST,SD,
     1     SU, TRACE)
C.......................................................................
C
C  MAIN ALGORITHM FOR THE LOGISTIC REGRESSION
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),WA(N),F0(N),F1(N),F2(N),CI(N),VTHETA(N),
     +     THETA(NP),DELTA(NP),GRAD(NP),COV(NCOV),HESSNV(NCOV),
     +     SC(NCOV),SE(NCOV),SF(NP),SG(NP),SH(NP),SX(MDX,NP),OI(N),
     +     A(NCOV),SA(NCOV),SD(NP),ST(NCOV),SU(N)
      INTEGER SP(NP),NI(N),RLICTHM2, TRACE, OURNIT(2)
      EXTERNAL RLICTHM2
      DATA ZERO/0.D0/
C
C avoid compiler warning
C

      idummy = TRACE

C
C  STEP 0 : INITIALIZATIONS
C  ------

      OURNIT(2) = MAXIT
      NIT=1
      DO 10 I=1,N
         CI(I)=ZERO
 10   CONTINUE
      DO 40 I=1,N
         DO 20 J=1,NP
            SD(J)=X(I,J)
 20      CONTINUE
         CALL RLMLYDBI(A,SD,NP,NCOV,NP,1)
         CALL RLNRM2BI(SD,NP,1,NP,ZNR)
         IF (ZNR.GT.ZMIN) GOTO 30
         ZNR=ZMIN
 30      WA(I)=B/ZNR
 40   CONTINUE
C
C  STEP 1 : COMPUTE THETA
C  ------
 100  DO 150 I=1,NP
         SD(I)=THETA(I)
 150  CONTINUE
C
C     is this theta step working? it doesn't look like it is...
C

C     ICASE==1 -> Binomial
C     ICASE==2 -> Bernulli
C     ICASE==3 -> Poisson

C
C     It looks like IOPT==1
C
      CALL RLGYTST2(X,Y,CI,THETA,WA,COV,NI,OI,N,NP,MDX,NCOV,GAM,
     +     TOLT,TAU,1D-6,15,IOPT,ICASE,ICNVT,MAXTT,
     +     NITT,Q0,DELTA,
     +     F0,F1,F2,VTHETA,GRAD,HESSNV,SE,SF,SG,SH,SC,SX,SP)

C
C  STEP 2 : CHECK CONVERGENCE
C  ------
      IF (NIT.EQ.MAXIT) GOTO 500

C
C     SD contains the old THETA (check line 100 above)
C

      DO 200 I=1,NP
         DELTA(I)=THETA(I)-SD(I)
 200  CONTINUE

C
C     RLICTHM2 checks if the vector DELTA is small enough
C

      IF (RLICTHM2(NP,NCOV,DELTA,1.D0,COV,TOL,ICNVT).EQ.1)
     +     RETURN

C
C  STEP 3 : COMPUTE THE A MATRIX AND THE ai's
C  ------
      CALL RLGYASTP(X,Y,NI,VTHETA,CI,A,OI,B,IUGL,ICASE,N,NP,NCOV,
     +     MDX,TAU,MAXTA,
     +     ICNVA,TOLA,NITA,WA,SU,SA,ST,SD)

      DO 340 I=1,N
         ZNR=WA(I)
         IF (ZNR.GT.ZMIN) GOTO 320
         ZNR=ZMIN
 320     WA(I)=B/ZNR
 340  CONTINUE

C
C  STEP 4 : COMPUTE THE ci's
C  ------
      CALL RLGICSTP(ICASE,IALG,NI,VTHETA,WA,OI,N,TOLC,MAXTC,CI)
C
C  STEP 5 : SET NIT:=NIT+1 AND GOTO STEP 1
C  ------
      NIT=NIT+1
      OURNIT(1) = NIT
      GOTO 100
 500  RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGYTST2(X,Y,CI,THETA,WA,COV,NI,OI,N,NP,MDX,NCOV,GAM,
     +     TOL,TAU,ZETA,IQ,IOPT,ICASE,ICNV,MAXIT,
     +     NIT,Q0,DELTA,F0,F1,F2,VTHETA,GRAD,HESSNV,SE,SF,SG,SH,ST,
     +     SX,IP)
C.......................................................................
C
C  NEWTON ALGORITHM FOR ROBUST LOGISTIC REGRESSION : THETA-STEP
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(NP),WA(N),COV(NCOV),DELTA(NP),
     +     F0(N),F1(N),F2(N),CI(N),VTHETA(N),GRAD(NP),HESSNV(NCOV),
     +     SE(NP),SF(NP),SG(NP),SH(NP),ST(NP),SX(MDX,NP),OI(N)
      INTEGER NI(N),IP(NP),RLICTHM2
      LOGICAL FIRST,NULF2
      EXTERNAL RLICTHM2

      INTCH=1
      SE(1)=0.D0

C
C  STEP 1.   SET NIT=1
C  ------
C
      NIT=1
C
C  STEP 2.   COMPUTE CURRENT OBJECTIVE FUNCTION VALUE AND (-)DERIVATIVES
C  ------
 200  continue
      CALL RLMFYD(X,THETA,VTHETA,N,NP,MDX,NP,1,N,1)
      CALL RLLRFNCT(ICASE,Y,CI,VTHETA,OI,WA,NI,N,1,1,1,F0,F1,F2,Q0)
C
C     we use ICASE in the line above
C

C
C  STEP 3.   COMPUTE THE NEGATIVE GRADIENT
C  ------
C
      CALL RLGRADBI(X,F1,N,NP,MDX,GRAD)
C
C  STEP 4.  COMPUTE THE GENERALIZED INVERSE OF THE NEGATIVE HESSIAN MAT.
C  ------
C
      FIRST=.TRUE.
 400  NULF2=.TRUE.
c      DO 410 I=1,N
c         SQF2=DSQRT(F2(I))
c         IF (SQF2.GT.1.D-5) NULF2=.FALSE.
c         DO 410 J=1,NP
c            SX(I,J)=X(I,J)*SQF2
c 410     CONTINUE
      DO I=1,N
        SQF2=DSQRT(F2(I))
        IF (SQF2.GT.1.D-5) NULF2=.FALSE.
        DO J=1,NP
          SX(I,J)=X(I,J)*SQF2
        END DO
      END DO
         IF (NULF2) THEN
            K=0
         ELSE
            CALL RLRMTRM2(SX,N,NP,MDX,INTCH,TAU,K,SF,SG,SH,IP)
         ENDIF
C
         L=0
         DO 425 I=1,NP
            DO 420 J=1,I
               L=L+1
               HESSNV(L)=0.D0
               IF (J.EQ.I) HESSNV(L)=1.D0
 420        CONTINUE
 425     CONTINUE
         IF (K.EQ.0) GOTO 500
         CALL RLKIASM2(SX,K,NP,MDX,NCOV,1.D0,1.D0,HESSNV)
         CALL RLKFASM2(SX,HESSNV,K,NP,MDX,NCOV,1.D0,
     +        DELTA,SG,IP)
C
C  STEP 5.   COMPUTE THE INCREMENT VECTOR
C  ------
C
 500     continue
         CALL RLMSFDBI(HESSNV,GRAD,DELTA,NP,NCOV,1,NP,NP)
         GAM0=GAM
         IF (K.LT.NP.AND.IOPT.NE.2) THEN
            CALL RLDOTPM2(DELTA,DELTA,NP,1,1,NP,NP,GAM0)
            IF (GAM0.GT.1.D0) GAM0=1.D0/DSQRT(GAM0)
         ENDIF
         DO 510 J=1,NP
            DELTA(J)=-DELTA(J)*GAM0
            IF (FIRST) THEN
               ST(J)=THETA(J)
               THETA(J)=THETA(J)+DELTA(J)
            ELSE
               SF(J)=THETA(J)
               THETA(J)=ST(J)+DELTA(J)
            ENDIF
 510     CONTINUE
C
C  STEP 6.   DETERMINE THE STEP-LENGTH
C  ------
C
         CALL RLMFYD(X,THETA,VTHETA,N,NP,MDX,NP,1,N,1)
         IF (.NOT.FIRST) Q01=Q0L
         CALL RLLRFNCT(ICASE,Y,CI,VTHETA,OI,WA,NI,N,1,0,0,F0,F1,
     +        F2,Q0L)
C
C     the above line only computes F0
C
C     we use ICASE in the line above
C
         IF (Q0L.LE.Q0) GOTO 700
         IF (.NOT.FIRST) GOTO 650
         FIRST=.FALSE.
C
C     It looks like IOPT==1
C
         IF (IOPT.EQ.1) THEN
            IF (ICASE.LE.2) THEN
               CALL RLDBINOM(Y,CI,VTHETA,WA,NI,F0,OI,N,1.D-6,F2)
            endif
            IF (ICASE.EQ.3) THEN
                 CALL RLDPOISS(Y,CI,VTHETA,WA,F0,OI,N,1.D-6,F2)
            endif
            GOTO 400
         ELSE
            CALL RLSTPLRG(ICASE,X,Y,CI,OI,ZETA,IQ,THETA,DELTA,WA,
     +           NI,GRAD,N,NP,MDX,Q0,Q01,GAM0,ST,F0,VTHETA)
            GOTO 700
         ENDIF
 650     IF (Q01.LT.Q0L) THEN
            DO 670 J=1,NP
               THETA(J)=SF(J)
 670        CONTINUE
         ENDIF

C
C  STEP 7. STOP ITERATIONS IF DESIRED PRECISION HAS BEEN REACHED
C  -------
 700     IF (NIT.EQ.MAXIT) GOTO 730
         IF (RLICTHM2(NP,NCOV,DELTA,1.D0,COV,TOL,ICNV).EQ.1) GOTO 730
         NIT=NIT+1
         GOTO 200
 730     continue
         CALL RLMFYD(X,THETA,VTHETA,N,NP,MDX,NP,1,N,1)
         CALL RLLRFNCT(ICASE,Y,CI,VTHETA,OI,WA,NI,N,1,1,1,F0,F1,F2,Q0)
         RETURN
         END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLMFYD(A,Y,Z,M,N,MDA,NY,IYE,NZ,IZE)
C.......................................................................
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(MDA,N),Y(NY),Z(NZ)

      NA1=(N-1)*MDA+1
      IZ=-IZE+1
      DO 20 I=1,M
C      DO 20 I=1,NZ
         IZ=IZ+IZE
         CALL RLDOTPM2(A(I,1),Y,N,MDA,IYE,NA1,NY,R)
         Z(IZ)=R
 20   CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RLDXLOG(X,XMIN,YMIN)
C.......................................................................
C
      DOUBLE PRECISION X,XMIN,YMIN
C
C  EXTENDED NATURAL LOGARITHM FUNCTION
C
      IF (X.LE.0.D0) THEN
         RLDXLOG=0.D0
      ELSEIF (X.LE.XMIN) THEN
         RLDXLOG=YMIN
      ELSE
        RLDXLOG=DLOG(X)
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLDBINOM(Y,CI,VTHETA,WA,NI,F0,OI,N,KAP,D)
C.......................................................................
C
C  APPROXIMATION FOR D(I) IN THE THETA STEP. BINOMIAL CASE.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER NI(N)
      DIMENSION Y(N),CI(N),VTHETA(N),OI(N),WA(N),F0(N),D(N)
      DOUBLE PRECISION KAP
      DATA NCALL,DMIN,XMIN,YMIN,DMAX/0,0.D0,0.D0,0.D0,0.D0/

      IF (NCALL.EQ.1) GOTO 10
      CALL RLMACHD(3,DMIN)
      CALL RLMACHD(4,XMIN)
      CALL RLMACHD(5,YMIN)
      CALL RLMACHD(6,XBIG)
      XBIG=XBIG/10.D0
      DMAX=DLOG(XBIG)
      NCALL=1
 10   CONTINUE
      DO 500 I=1,N
         GG=VTHETA(I)
         OO=OI(I)
         YY=Y(I)-CI(I)
         AI=WA(I)
         NII=NI(I)
         ENI=DBLE(NII)
         ENO=ENI
         IF (YY.GT.AI) THEN
            IF (-YY+ENI.LE.-AI) THEN
               D(I)=KAP
               GOTO 500
            ELSEIF (-YY+ENI.LE.AI) THEN
               IF (ENI.GT.YY) THEN
                  CALL RLTS12BI(YY,AI,ENO,OO,XMIN,YMIN,T01,S01,T2,S2)
                  IF (GG.GE.T01) GOTO 300
 100              ENI=ENI+1
                  IF (-YY+ENI.LE.AI) GOTO 100
                  CALL RLTS12BI(YY,AI,ENI,OO,XMIN,YMIN,T1,S1,T2,S2)
                  T01=GG+(S1-S01)/AI
                  DNI=ENI
                  DI=(S1-S2)/(2.D0*AI) - T01
                  CALL RLBIGGBI(T01,DNI,DMIN,DMAX,S2)
                  ELG0=(-YY*T01+ENI*S2)
                  CALL RLBIGGBI(T1,DNI,DMIN,DMAX,S2)
                  ELG1=(-YY*T1+ENI*S2)
                  CALL RLBIGGBI(T2,DNI,DMIN,DMAX,S2)
                  ELG2=(-YY*T2+ENI*S2)
                  IF (ELG0.LE.DMAX1(ELG1,ELG2)) GOTO 300
                  D(I)=DABS(AI/DI)
                  GOTO 500
               ELSE
                  GOTO 300
               ENDIF
            ENDIF
         ELSEIF (YY.LE.AI.AND.-YY.LT.AI) THEN
            S2=ENI/4.D0
            IF (-YY+ENI.GT.AI) THEN
               T2=(YY+AI)/(ENI-YY-AI)
               T2=RLDXLOG(T2,XMIN,YMIN)-OO
               IF (T2.LT.0.D0) THEN
                  S2=RLGFUN(2,1,T2+OO)
                  S2=S2*(1.D0-S2)*ENI
               ENDIF
            ENDIF
            D(I)=S2
            GOTO 500
         ELSEIF (-YY.GE.AI) THEN
            D(I)=KAP
            GOTO 500
         ENDIF
         CALL RLTS12BI(YY,AI,ENI,OO,XMIN,YMIN,T1,S1,T2,S2)
         ELG0=F0(I)
         ELG1=(-AI*T1+S1)
         ELG2=( AI*T2+S2)
         IF (ELG0.GE.DMAX1(ELG1,ELG2)) THEN
            DI=(S1-S2)/(2.D0*AI)-GG
            D(I)=DABS(AI/DI)
         ELSE
            IF (T2.LE.0.D0) THEN
               S2=RLGFUN(2,1,T2+OO)
               D(I)=S2*(1.D0-S2)*ENI
            ELSEIF (T1.LE.0.D0.AND.T2.GT.0.D0) THEN
               D(I)=ENI/4.D0
            ELSE
               S1=RLGFUN(2,1,T1+OO)
               D(I)=S1*(1.D0-S1)*ENI
            ENDIF
         ENDIF
         GOTO 500
 300     T1=(YY-AI)/(ENI-YY+AI)
         T1=RLDXLOG(T1,XMIN,YMIN)-OO
         S2=ENI/4.D0
         IF (T1.GT.0.D0) THEN
            S2=RLGFUN(2,1,T1+OO)
            S2=S2*(1.D0-S2)*ENI
         ENDIF
         D(I)=S2
 500  CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLBIGGBI(X,DNI,DMIN,DMAX,Y)
C.......................................................................
C
      DOUBLE PRECISION DNI,DMIN,DMAX,X,Y
C
C  AUXILIARY SUBROUTINE FOR RLDBINOM (BINOMIAL CASE).
C
      IF (X.LE.DMIN) THEN
         Y=0.D0
      ELSEIF (X.GE.DMAX) THEN
         Y=DNI*X
      ELSE
        Y=DNI*DLOG(1.D0+DEXP(X))
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLTS12BI(YI,AI,DNI,OI,XMIN,YMIN,T1,S1,T2,S2)
C.......................................................................
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C  AUXILIARY SUBROUTINE FOR RLDBINOM (BINOMIAL CASE).
C
      T1=(YI-AI)/(DNI-YI+AI)
      T1=RLDXLOG(T1,XMIN,YMIN)-OI
      S1=DNI/(DNI-YI+AI)
      S1=RLDXLOG(S1,XMIN,YMIN)
      S1=-(YI-AI)*T1+DNI*S1
      T2=(YI+AI)/(DNI-YI-AI)
      T2=RLDXLOG(T2,XMIN,YMIN)-OI
      S2=DNI/(DNI-YI-AI)
      S2=RLDXLOG(S2,XMIN,YMIN)
      S2=-(YI+AI)*T2+DNI*S2
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLDPOISS(Y,CI,VTHETA,WA,F0,OI,N,KAP,D)
C.......................................................................
C
C  APPROXIMATION FOR D(I) IN THE THETA STEP. POISSON CASE.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(N),CI(N),VTHETA(N),OI(N),WA(N),F0(N),D(N)
      DOUBLE PRECISION KAP
      DATA NCALL,DMIN,XMIN,YMIN,DMAX/0,0.D0,0.D0,0.D0,0.D0/
C
C     They've just set up NCALL==0!!
C
      IF (NCALL.EQ.1) GOTO 10
      CALL RLMACHD(3,DMIN)
      CALL RLMACHD(4,XMIN)
      CALL RLMACHD(5,YMIN)
      CALL RLMACHD(6,XBIG)
      XBIG=XBIG/10.D0
      DMAX=DLOG(XBIG)
      NCALL=1
 10   CONTINUE
      DO 500 I=1,N
         GG=VTHETA(I)
         OO=OI(I)
         YY=Y(I)-CI(I)
         AI=WA(I)
         IF (YY.GT.AI) THEN
            CALL RLTS12PO(YY,AI,OO,XMIN,YMIN,T1,S1,T2,S2)
         ELSEIF (YY.GT.-AI) THEN
            T2=RLDXLOG(YY+AI,XMIN,YMIN)-OO
            IF (GG.LE.T2.OR.YY.LE.0.D0) GOTO 300
            ENO=RLDXLOG(YY,XMIN,YMIN)
            IF (F0(I).LT.AI*(ENO-GG)/2.D0+YY*(1.D0-ENO)) GOTO 300
            D(I)=AI/DABS(ENO-GG)
            GOTO 500
         ELSE
            D(I)=KAP
            GOTO 500
         ENDIF
         ELG0=F0(I)
         ELG1=(-AI*T1+S1)
         ELG2=( AI*T2+S2)
         IF (ELG0.GE.DMAX1(ELG1,ELG2)) THEN
            DI=(S1-S2)/(2*AI)-GG
            D(I)=DABS(AI/DI)
            GOTO 500
         ENDIF
 300     IF (T2+OO.LE.DMIN) THEN
            D(I)=KAP
         ELSEIF (T2+OO.GT.DMAX) THEN
            D(I)=XBIG
         ELSE
            D(I)=DEXP(T2+OO)
         ENDIF
 500  CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLTS12PO(YI,AI,OI,XMIN,YMIN,T1,S1,T2,S2)
C.......................................................................
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C  AUXILIARY SUBROUTINE FOR DPOISS (POISSON CASE).
C
      T1=RLDXLOG(YI-AI,XMIN,YMIN)-OI
      S1=-(YI-AI)*T1+(YI-AI)
      T2=RLDXLOG(YI+AI,XMIN,YMIN)-OI
      S2=-(YI+AI)*T2+(YI+AI)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLLRFNCT(ICASE,Y,C,VTHETA,OI,WA,NN,N,I0,I1,I2,F0,F1,
     +     F2,SF0)
C.......................................................................
C
C  FUNCTIONS L, L' and L" FOR THE THETA STEP
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(N),C(N),VTHETA(N),OI(N),WA(N),NN(N),F0(N),F1(N),
     +     F2(N)
      DATA NCALL,DMIN,XMIN,YMIN,DMAX/0,0.D0,0.D0,0.D0,0.D0/
C
C     They've just set NCALL==0!!
C
      IF (NCALL.EQ.1) GOTO 10
      CALL RLMACHD(3,DMIN)
      CALL RLMACHD(4,XMIN)
      CALL RLMACHD(5,YMIN)
      CALL RLMACHD(6,XBIG)
      XBIG=XBIG/10.D0
      DMAX=DLOG(XBIG)
      NCALL=1
 10   SUM=0.D0
      DO 500 I=1,N
         GI=VTHETA(I)
         OF=OI(I)
         GO=GI+OF
         YI=Y(I)-C(I)
         AI=WA(I)
         NI=1
         IF (ICASE.EQ.2) NI=NN(I)
         ENI=DBLE(NI)
         IF (-YI.GE.AI) THEN
            S2=0.D0
            GOTO 200
         ENDIF
         IF (ICASE.EQ.3) GOTO 30
C==>    BERNOUILLI OR BINOMIAL CASE
         IF (-YI.LT.-AI) THEN
            IF (-YI+ENI.LE.-AI) THEN
               S1=0.D0
               GOTO 100
            ENDIF
            T1=(YI-AI)/(ENI-YI+AI)
            T1=RLDXLOG(T1,XMIN,YMIN)-OF
            IF (GI.LT.T1) THEN
               S1=ENI/(ENI-YI+AI)
               S1=RLDXLOG(S1,XMIN,YMIN)
               S1=-(YI-AI)*T1+ENI*S1
               GOTO 100
            ENDIF
            GOTO 50
         ELSE
            GOTO 50
         ENDIF
 30      IF (-YI.LT.-AI) THEN
C==>    POISSON CASE
            T1=YI-AI
            T2=YI+AI
            T1=RLDXLOG(T1,XMIN,YMIN)-OF
            T2=RLDXLOG(T2,XMIN,YMIN)-OF
            IF (GI.LT.T1) THEN
               S1=-(YI-AI)*T1+(YI-AI)
               GOTO 100
            ELSEIF (GI.GT.T2) THEN
               S2=-(YI+AI)*T2+(YI+AI)
               GOTO 200
            ELSE
               GOTO 400
            ENDIF
         ELSE
            T2=YI+AI
            T2=RLDXLOG(T2,XMIN,YMIN)-OF
            IF (GI.LE.T2) GOTO 400
            S2=-(YI+AI)*T2+(YI+AI)
            GOTO 200
         ENDIF
 50      IF (-YI+ENI.LE.AI) GOTO 300
         T2=(YI+AI)/(ENI-YI-AI)
         T2=RLDXLOG(T2,XMIN,YMIN)-OF
         IF (GI.LE.T2) GOTO 300
         S2=ENI/(ENI-YI-AI)
         S2=RLDXLOG(S2,XMIN,YMIN)
         S2=-(YI+AI)*T2+ENI*S2
         GOTO 200
 100     F0I=-AI*GI+S1
         IF (I0.NE.0) F0(I)=F0I
         IF (I1.NE.0) F1(I)=-AI
         IF (I2.NE.0) F2(I)=0.D0
         GOTO 450
 200     F0I=AI*GI+S2
         IF (I0.NE.0) F0(I)=F0I
         IF (I1.NE.0) F1(I)=AI
         IF (I2.NE.0) F2(I)=0.D0
         GOTO 450
 300     IF (GI+OF.LE.DMIN) THEN
            S1=0.D0
         ELSEIF (GI+OF.GE.DMAX) THEN
            S1=ENI*(GI+OF)
         ELSE
            S1=ENI*DLOG(1.D0+DEXP(GI+OF))
         ENDIF
         S2=RLGFUN(ICASE,1,GO)
         F0I=-YI*GI+S1
         IF (I0.NE.0) F0(I)=F0I
         IF (I1.NE.0) F1(I)=-YI+ENI*S2
         IF (I2.NE.0) F2(I)=ENI*S2*(1.D0-S2)
         GOTO 450
 400     EXPGI=RLGFUN(ICASE,1,GO)
         F0I=-YI*GI+EXPGI
         IF (I0.NE.0) F0(I)=F0I
         IF (I1.NE.0) F1(I)=-YI+EXPGI
         IF (I2.NE.0) F2(I)=EXPGI
 450     SUM=SUM+F0I
 500  CONTINUE
      SF0=SUM
      RETURN
      END
C
C-----------------------------------------------------------------------
      SUBROUTINE RLSTPLRG(ICASE,X,Y,C,OI,ZETA,IQ,THETA,DELTA,WA,NI,
     +     GRAD,N,NP,MDX,SF0,SF1,GAM,ST,F0,VTHETA)
C.......................................................................
C
C  ARMIJO-GOLDSTEIN STEP LENGTH ALGORITHM
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),C(N),THETA(NP),DELTA(NP),WA(N),NI(N),
     +     GRAD(NP),ST(NP),F0(N),VTHETA(N),OI(N)

      CALL RLDOTPM2(DELTA,GRAD,NP,1,1,NP,NP,S0)
      IF (DABS(S0).GT.1.D-5) GOTO 10
      ETA=1.D0
      DO 450 NLOOP=1,IQ
         ETA=ETA*0.5D0
         DO 400 J=1,NP
            ST(J)=THETA(J)+ETA*DELTA(J)
 400     CONTINUE
         CALL RLMFYD(X,ST,VTHETA,N,NP,MDX,NP,1,N,1)
         CALL RLLRFNCT(ICASE,Y,C,VTHETA,OI,WA,NI,N,1,0,0,F0,F0,F0,SF1)
         IF (SF1.LT.SF0) GOTO 500
 450  CONTINUE
      GOTO 500
 10   IP=0
      ETA=2.0
 100  IF (IP.EQ.IQ) GOTO 500
      ETA=0.5D0**IP
      DO 200 J=1,NP
         ST(J)=THETA(J)+ETA*DELTA(J)
 200  CONTINUE
      CALL RLMFYD(X,ST,VTHETA,N,NP,MDX,NP,1,N,1)
      CALL RLLRFNCT(ICASE,Y,C,VTHETA,OI,WA,NI,N,1,0,0,F0,F0,F0,SF1)
      IF (SF1.LT.SF0) GOTO 500
      IF ((SF1-SF0)/ETA/S0.GT.ZETA) GOTO 500
      IP=IP+1
      GOTO 100
 500  GAM=ETA
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGYASTP(X,Y,NI,VTHETA,CI,A,OI,B,IUGL,ICASE,NOBS,NVAR,
     +     NCOV,MDX,TAU,MAXIT,
     +     ICNV,TOL,NIT,DIST,SU,SA,ST,SD)
C.......................................................................
C
C  FIXED POINT ALGORITHM FOR THE COMPUTATION OF THE MATRIX A
C  (STANDARDIZED CASE, A LOWER TRIANGULAR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(MDX,NVAR),Y(NOBS),OI(NOBS),VTHETA(NOBS),CI(NOBS),
     +     DIST(NOBS),A(NCOV),SA(NCOV),ST(NCOV),SD(NVAR),SU(NOBS)
      INTEGER NI(NOBS),RLICNVBI
C
C  STEP 0 : INITIALIZATION
C  ------

      NIT=0
      IF (ICNV.EQ.1) THEN
         L=0
         DO 20 I=1,NVAR
            DO 10 J=1,I
               L=L+1
               SA(L)=0.D0
               IF (I.EQ.J) SA(L)=-1.D0
 10         CONTINUE
 20      CONTINUE
      ENDIF
      DO 30 L=1,NOBS
         DIST(L)=-1.D0
 30   CONTINUE
C
C  STEP 1: COMPUTE WEIGHTED COVARIANCE (ST) AND AUXILIARY VALUES
C  ------
      CALL RLUCOWJ(X,Y,NI,VTHETA,OI,CI,A,ST,NOBS,NVAR,NCOV,MDX,
     +     ICNV,NIT,DELTA,DIST,SU,SD,IUGL,B,ICASE)

C
C  STEP 2: CHECK CONVERGENCE
C  ------
      IF (NIT.EQ.MAXIT.OR.RLICNVBI(NCOV,DELTA,A,SA,TOL,ICNV).EQ.1)
     +     GOTO 500
C
C  STEP 3: FIND IMPROVEMENT MATRIX SS=I-ST FOR A
C  ------
      INFO=0
      CALL RLPRSFBI(ST,NVAR,NCOV,TAU,INFO)

C
C  STEP 4: SET SA:=A AND A:=(I-SS)*SA
C  -------
      DO 410 IJ=1,NCOV
         SA(IJ)=A(IJ)
 410  CONTINUE
      CALL RLMTT3BI(SA,ST,A,NVAR,NCOV)
      NIT=NIT+1

C
C  STEP 5: EXIT
C  ------
 500  RETURN
      END
C
C-----------------------------------------------------------------------
C
      double precision FUNCTION RLUGL(UPAR,S,IUGL,ICASE,B)
C.......................................................................
C
C  PURPOSE
C  -------
C  WEIGHT FUNCTION FOR THE A-STEP
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      DIMENSION UPAR(NPAR)
      DIMENSION UPAR(4)
      DATA STOL,PREC/1.D-6,0.D0/

      IF (PREC.EQ.0.D0) CALL RLMACHD(2,PREC)
      YI=UPAR(1)
      ENI=UPAR(2)
      GI=UPAR(3)
      CI=UPAR(4)
      NI=INT(ENI+0.001)
      SI=S
      IF (S.LE.STOL) SI=STOL
      A=B/SI
      IF (IUGL.EQ.1) THEN
C
C     it seems IUGL==1 (always)
C
C====>  OPTION 1
         PP=RLGFUN(ICASE,1,GI)
         IF (ICASE.EQ.1) THEN
C==>      LOGISTIC BERNOULLI
            TEMP=DABS(1.D0-PP-CI)
            T1=A**2
            IF (TEMP.LT.A) T1=TEMP**2
            TEMP=DABS(-PP-CI)
            T2=A**2
            IF (TEMP.LT.A) T2=TEMP**2
            RLUGL=T1*PP+T2*(1.D0-PP)
         ELSEIF (ICASE.EQ.2) THEN
C==>      LOGISTIC BINOMIAL
            T2=0.D0
            DO 200 J=0,NI
               CALL RLPROBIN(J,NI,PP,PJ)
               TEMP=DABS(DBLE(J)-ENI*PP-CI)
               T1=A**2
               IF (TEMP.LT.A) T1=TEMP**2
               T2=T2+T1*PJ
 200        CONTINUE
            RLUGL=T2
         ELSEIF (ICASE.EQ.3) THEN
C==>      LOGISTIC POISSON
            MI=INT(100.D0*PP)
            T2=0.D0
            DO 300 J=0,MI
               CALL RLPRPOIS(PP,J,PJ)
               FJME=DBLE(J)-PP
               TEMP=DABS(FJME-CI)
               T1=A**2
               IF (TEMP.LT.A) T1=TEMP**2
               IF (FJME.GT.0.D0 .AND. T1*PJ.LT.PREC) GOTO 350
               T2=T2+T1*PJ
 300        CONTINUE
 350        RLUGL=T2
         ENDIF
      ELSE
C====>  OPTION 2 (IUGL=2)
         PP=RLGFUN(ICASE,NI,GI)
         TEMP=DABS(YI-PP-CI)
         RLUGL=A**2
         IF (TEMP.LT.A) RLUGL=TEMP**2
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLUCOWJ(X,Y,NI,VTHETA,OI,CI,SA,ST,N,NP,NCOV,MDX,
     +     ICNT,NIT,ZMAX,DIST,SU,SD,IUGL,B,ICASE)
C.......................................................................
C
C  COMPUTE WEIGHTED COVARIANCE MATRIX FOR LOGISTIC REGRESSION;
C  (STORE EXUL VALUES IN SU)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),CI(N),VTHETA(N),OI(N),DIST(N),UARR(4),
     +     SA(NCOV),ST(NCOV),SU(N),SD(NP)
      INTEGER NI(N)
      DATA XN/0.D0/
      IF (NIT.GT.1) GOTO 10
      XN=DBLE(N)
 10   ZMAX=0.D0
      DO 20 IJ=1,NCOV
         ST(IJ)=0.D0
 20   CONTINUE
      YL=1.D0
      NL=1
      DO 100 L=1,N
         DO 50 J=1,NP
            SD(J)=X(L,J)
 50      CONTINUE
         CALL RLMLYDBI(SA,SD,NP,NCOV,NP,1)
         CALL RLNRM2BI(SD,NP,1,NP,DISTL)
         IF (ICNT.EQ.2) ZMAX=DMAX1(ZMAX,DABS(DISTL-DIST(L)))
         DIST(L)=DISTL
         GL=VTHETA(L)+OI(L)
         CL=CI(L)
         IF (IUGL.EQ.1) YL=Y(L)
         IF (ICASE.EQ.2) NL=NI(L)
         UARR(1)=YL
         UARR(2)=DBLE(NL)
         UARR(3)=GL
         UARR(4)=CL
         U=RLUGL(UARR,DISTL,IUGL,ICASE,B)
         SU(L)=U
         IJ=0
         DO 90 I=1,NP
            DO 85 J=1,I
               IJ=IJ+1
               ST(IJ)=ST(IJ)+(SD(I)*U)*SD(J)
 85         CONTINUE
 90      CONTINUE
 100  CONTINUE
      DO 110 IJ=1,NCOV
         ST(IJ)=ST(IJ)/XN
 110  CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGICSTP(ICASE,IALG,NN,VTHETA,WA,OI,N,TOL,MAXIT,C)
C.......................................................................
C
C  H-,W-,S-ALGORITHM FOR ROBUST LOGISTIC REGRESSION : C-STEP solution
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION VTHETA(N),OI(N),WA(N),C(N)
      INTEGER NN(N)


      DO 500 I=1,N
         GI=VTHETA(I)+OI(I)
         A=WA(I)
         NI=1
         IF (ICASE.EQ.2) NI=NN(I)
         PP=RLGFUN(ICASE,NI,GI)
         T=C(I)+PP
         CALL RLGYCSTP(ICASE,IALG,NI,A,PP,TOL,MAXIT,T)
         C(I)=T-PP
 500  CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGYCSTP(ICASE,IALG,NI,A,E,TOL,MAXIT,T)
C.......................................................................
C
C  NEWTON-TYPE ALGORITHM FOR THE C-STEP
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER H
      DATA INICA,INIAL,PREC/0,0,0.D0/

      IF (INICA.EQ.ICASE.AND.INIAL.EQ.IALG) GOTO 10
      INICA=ICASE
      INIAL=IALG
      CALL RLMACHD(2,PREC)

C
C  STEP 0.   SET NIT=1
C  ------

 10   NIT=1
C
C  STEP 1.   COMPUTE CURRENT OBJECTIVE FUNCTION VALUE
C  ------
      IF (ICASE.EQ.1) THEN
C==>    LOGISTIC BERNOULLI : SOLVE EXPLICITLY
         P1=E
         P0=1.D0-P1
         T=0.D0
         IF (P1.LT.0.5D0 .AND. A.LT.P0) THEN
            T=A*P1/P0 - P1
         ELSEIF (P1.GT.0.5D0 .AND. A.LT.P1) THEN
            T=P0-A*P0/P1
         ENDIF
         T=T+E
         RETURN
      ENDIF
      IF (ICASE.EQ.2) PI=E/DBLE(NI)
      M=NI
      IF (ICASE.EQ.3) M=MAX0(INT(100.*E),5000)
 100  H=-1
      K=-1
      J1=0
      J2=M
      IF (IALG.EQ.1) GOTO 200
      IF (IALG.EQ.2) GOTO 200
C
C  STEP 2.  COMPUTE AUXILIARY QUANTITIES FOR STEP 2, 4 AND 5.
C  _______
      AUX=T-A
      H=INT(AUX)
      IF (AUX.LT.0.D0 .AND. AUX.NE.DBLE(H)) H=H-1
      H=MAX0(H,-1)
      AUX=T+A
      K=INT(AUX)
      IF (AUX.LT.0.D0 .AND. AUX.NE.DBLE(K)) K=K-1
      K=MIN0(K,NI)
      J1=0
      J2=H
 200  SH=0.D0
      EH=0.D0
      SK=0.D0
      EK=0.D0
      S1=0.D0
      E1=0.D0
      T1=0.D0
      DJBAR=0.D0
 210  IF (ICASE.EQ.2) THEN
C==>    LOGISTIC BINOMIAL
         DO 220 J=J1,J2
            CALL RLPROBIN(J,NI,PI,PJ)
            DEN=(DBLE(J)-T)
            TEMP=DMIN1(A,DABS(DEN))
            IF (DEN.LT.0.D0) TEMP=-TEMP
            T1=T1+TEMP*PJ
            E1=E1+DBLE(J)*PJ
            IF (IABS(IALG).NE.2) GOTO 220
            DJ=PJ
            IF (DABS(DEN).GT.1.D-6) DJ=TEMP*PJ/DEN
            DJBAR=DJBAR+DJ
 220     CONTINUE
      ELSEIF (ICASE.EQ.3) THEN
C==>    LOGISTIC POISSON
         DO 230 J=J1,J2
            CALL RLPRPOIS(E,J,PJ)
            DEN=(DBLE(J)-T)
            TEMP=DMIN1(A,DABS(DEN))
            IF (DEN.LT.0.D0) TEMP=-TEMP
            TMP=TEMP*PJ
            IF (DABS(TMP).LT.PREC) TMP=0.D0
            T1=T1+TMP
            TMPJ=DBLE(J)*PJ
            IF (DABS(TMPJ).LT.PREC) TMPJ=0.D0
            E1=E1+TMPJ
            IF (ABS(IALG).NE.2) GOTO 225
            DJ=PJ
            IF (DABS(DEN).GT.1.D-6) DJ=TMP/DEN
            DJBAR=DJBAR+DJ
 225        IF (TMP.EQ.0.D0.AND.TMPJ.EQ.0.D0) GOTO 240
 230     CONTINUE
      ENDIF
 240  IF (IALG.EQ.1) GOTO 400
      IF (IALG.EQ.2) GOTO 500
      IF (J1.EQ.0.AND.J2.EQ.H) THEN
         IF (H.NE.-1) SH=S1
         IF (H.NE.-1) EH=E1
         J1=H+1
         J2=K
         GOTO 210
      ELSEIF (J1.EQ.H+1.AND.J2.EQ.K) THEN
         SK=S1
         EK=E1
         J1=K+1
         J2=NI
         GOTO 210
      ENDIF
      IF (K.EQ.NI) THEN
         SK=1.D0
         EK=E
      ENDIF
      TSTAR=A*(1.D0-SH-SK)+(EK-EH)
      DEN=SK-SH
      IF (DABS(DEN).LE.PREC) DEN=PREC
      TSTAR=TSTAR/DEN
C
C  STEP 3.  VERIFY IF H AND K ARE GOOD
C  _______
      TH=DBLE(H)-TSTAR
      IF (H+1.EQ.0) TH=-A-1.D0
      THP1=TH+1.D0
      TK=DBLE(K)-TSTAR
      TKP1=TK+1
      IF (K.EQ.NI) TKP1=A+1.D0
      IF (TH.LE.-A .AND. -A.LT.THP1 .AND. TK.LE.A.AND.A.LT.TKP1)
     +     THEN
         T=TSTAR
         RETURN
      ENDIF
C
C  STEP 4.  H-ALGORITHM : SET T:=T0 + DELTA
C  _______
 400  DELTA=T1
      IF (ABS(IALG).EQ.1) T=T+DELTA
C
C   STEP 5.  W-ALGORITHM : SET T:=T0 + DELTA/DBAR
C   _______
 500  IF (ABS(IALG).EQ.2) THEN
         IF (DABS(DJBAR).LE.1.D-5) DJBAR=DSIGN(1.D0,DJBAR)
         DELTA=T1/DJBAR
         T=T+DELTA
      ENDIF
C
C   STEP 6.   CHECK CONVERGENCE
C   ------
C
      IF (DABS(DELTA).LT.TOL.OR.NIT.EQ.MAXIT) GOTO 600
      NIT=NIT+1
      GOTO 100
 600  RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGLMDEV(Y,NI,CI,WA,VTHETA,OI,N,ICASE,DEV,THETAS,LI,SC)
C.......................................................................
C
C     NEW VERSION (10/99)
C.......................................................................
C
C     GLM DEVIANCE COMPUTATION
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(N),CI(N),WA(N),VTHETA(N),OI(N),THETAS(N),SC(N)
      DOUBLE PRECISION LI(N)
      INTEGER NI(N)

      CALL RLLRFNCT(ICASE,Y,CI,VTHETA,OI,WA,NI,N,1,0,0,LI,WA,WA,Q)
      DO 700 I=1,N
         TMP=(Y(I)-CI(I))/DBLE(NI(I))
         THETAS(I)=RLFLINK(ICASE,TMP)-OI(I)
 700  CONTINUE
      QS=0.D0
      DO 800 I=1,N
         ENI=DBLE(NI(I))
         YI=Y(I)
         IF (ICASE.LE.2) THEN
            TMP=ENI*DLOG(ENI)
            IF (YI.GT.0.D0) TMP=TMP-YI*DLOG(YI)
            ENYI=ENI-YI
            IF (ENYI.GT.0.D0) TMP=TMP-ENYI*DLOG(ENYI)
         ELSE
            TMP=YI
            IF (YI.GT.0.D0) TMP=TMP-YI*DLOG(YI)
         ENDIF
         QS=QS+TMP
         SC(I)=TMP
 800  CONTINUE
      DEV=2.D0*DABS(Q-QS)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      double precision FUNCTION RLFLINK(ICASE,EM)
C.......................................................................
C
C     NEW VERSION (10/99)
C.......................................................................
C
C     LINK FUNCTION FOR LOGISTIC REGRESSION
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DATA NCALL,XMIN,YMIN/0,0.D0,0.D0/

      IF (NCALL.EQ.1) GOTO 10
      CALL RLMACHD(4,XMIN)
      CALL RLMACHD(5,YMIN)
      NCALL=1
 10   RLFLINK=-999.D0
      IF (EM.LE.0.D0) RETURN
      TT=YMIN
      IF (EM.GT.XMIN) TT=DLOG(EM)
      TMP=0.D0
      IF (ICASE.GT.2) GOTO 20
      IF (1.D0-EM.LE.0.D0) RETURN
      TMP=YMIN
      IF (1.D0-EM.GT.XMIN) TMP=DLOG(1.D0-EM)
 20   RLFLINK=TT-TMP
      RETURN
      END


C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGYTSTP(X,Y,CI,THETA,WA,COV,NI,OI,N,NP,MDX,NCOV,GAM,
     1           TOL,TAU,IOPT,ICASE,ICNV,MAXIT,NITMON,NIT,Q0,DELTA,
     2           F0,F1,F2,VTHETA,GRAD,HESSNV,RW1,RW2,IW1)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. RANDRIAMIHARISOA / A. MARAZZI
C.......................................................................
C
C  NEWTON ALGORITHM FOR ROBUST LOGISTIC REGRESSION : THETA-STEP
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(NP),OI(N),WA(N),COV(NCOV),
     1          F0(N),F1(N),F2(N),CI(N),VTHETA(N),GRAD(NP),DELTA(NP),
     2          RW1(5*NP),RW2(MDX,NP),HESSNV(NCOV)
      INTEGER NI(N),IW1(NP)
      DATA ZETA,IQ/0.D001,10/
C
C avoid compiler warnings
C

      idummy = nitmon

C
C  PARAMETER CHECK
C
C      NN=NP*(NP+1)/2
C      NPRCHK=NP.GT.0.AND.NP.LE.N.AND.MDX.GE.N.AND.NN.EQ.NCOV
C     1       .AND.GAM.GT.0..AND.GAM.LT.2..AND.TAU.GE.0..AND.
C     2       TOL.GT.0..AND.(ICASE.EQ.1.OR.ICASE.EQ.2.OR.ICASE.EQ.3)
C     3       .AND.MAXIT.GT.0.AND.(IOPT.EQ.1.OR.IOPT.EQ.2).AND.
C     4       (ICNV.EQ.1.OR.ICNV.EQ.2.OR.ICNV.EQ.3)
C      IF (.NOT.NPRCHK) CALL MESSGE(500,'GYTSTP',1)
      ISF=NP+1
      ISG=ISF+NP
      ISH=ISG+NP
      IST=ISH+NP



      CALL RLGYTST2(X,Y,CI,THETA,WA,COV,NI,OI,N,NP,MDX,NCOV,GAM,TOL,
     +TAU,ZETA,IQ,IOPT,ICASE,ICNV,MAXIT,NIT,Q0,DELTA,F0,
     +F1,F2,VTHETA,GRAD,HESSNV,RW1,RW1(ISF),RW1(ISG),RW1(ISH),
     +RW1(IST),RW2,IW1)

      RETURN
      END
C
C-----------------------------------------------------------------------
C




