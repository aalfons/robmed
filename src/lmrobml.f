C=======================================================================
C     Robust Efficient Weighted Least Squares (REWLS)
C     AUTHOR: JEFFREY WANG
C     MATHSOFT, INC.
C     10/06/99
C=======================================================================
      SUBROUTINE RLCHI1ML(X,P)
C.......................................................................
      DOUBLE PRECISION X,P,TMP,ZERO
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     CALCULATE THE CDF FOR ABS(NORMAL R.V.)
C     NOTE: X >= 0
C-----------------------------------------------------------------------
      IF (X .LT. ZERO) X=ZERO
      IF (X .EQ. ZERO) THEN
         P=ZERO
         RETURN
      ENDIF
      CALL RLGAUSBI( X,P)
      CALL RLGAUSBI(-X,TMP)
      P = P-TMP
      RETURN
      END
C=======================================================================
      SUBROUTINE RLRWETML(X,P)
C.......................................................................
      DOUBLE PRECISION X,AX,P,COEF,ZERO
      DIMENSION COEF(4)
      DATA ZERO/0.D0/
      DATA COEF/-19.7187928669416D0,
     +           82.3045267489739D0,
     +         -105.4526748971229D0,
     +           42.8669410150906D0/
C-----------------------------------------------------------------------
C     CALCULATE WEIGHT FUNCTION FOR REWEIGHTING
C     NOTE: X >= 0
C-----------------------------------------------------------------------
      AX = DABS(X)
      IF (AX .GE. 1.D0) THEN
         P = ZERO
      ELSE IF (AX .LE. 0.8D0) THEN
         P = 1.D0
      ELSE
         P = COEF(1)+COEF(2)*AX**2+COEF(3)*AX**4+COEF(4)*AX**6
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE RLRWEPML(X,P)
C.......................................................................
      DOUBLE PRECISION X,AX,P,COEF,ZERO
      DIMENSION COEF(4)
      DATA ZERO/0.D0/
      DATA COEF/-19.7187928669416D0,
     +           82.3045267489739D0,
     +         -105.4526748971229D0,
     +           42.8669410150906D0/
C-----------------------------------------------------------------------
C     CALCULATE THE DERIVATIVE OF WEIGHT FUNCTION
C-----------------------------------------------------------------------
      AX=DABS(X)
      IF (AX .GE. 1.D0 .OR. AX .LE. 0.8D0) THEN
         P=ZERO
      ELSE
         P=2.D0*COEF(2)*X+4.D0*COEF(3)*X**3+6.D0*COEF(4)*X**5
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE RLFINLML(X,Y,WGT,RS,N,NP,MDX,THETA,SCAL,SF,SG,
     +     SH,IP,SX,SY,TAU,ETA,IERR,IPS,XK,FAC,U)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),SX(MDX,NP),SY(N),WGT(N),RS(N)
      DIMENSION THETA(N),SF(NP),SG(NP),SH(NP),U(N)
      INTEGER IP(NP)
      DATA ZERO,ONE/0.D0,1.D0/
C-----------------------------------------------------------------------
C     COMPUTES THE FINAL EFFICIENT WEIGHTED LS
C-----------------------------------------------------------------------
C     STEP 1. COMPUTE THE CUTOFF VALUE T
C-----------------------------------------------------------------------
      DO 10 I=1,N
         DM = RS(I)/SCAL
         U(I) = DM
         WGT(I)=DABS(DM)
 10   CONTINUE
      CALL RLSRT1BI(WGT,N,1,N)
      DN = DBLE(N)
      DM = ZERO
      DO 20 I=1,N
         CALL RLCHI1ML(WGT(I),TMP)
         TMP = DMAX1(TMP-(DBLE(I)-ONE)/DN,ZERO)
         IF (TMP .GT. DM) DM = TMP
 20   CONTINUE
      T = DMAX1(ETA,WGT(N-INT(N*DM)))
C-----------------------------------------------------------------------
C     STEP 2. COMPUTE WEIGHTS AND OTHER VALUES FOR FAC
C-----------------------------------------------------------------------
      PP = ZERO
      HT = ZERO
      H1T = ZERO
      DO 25 I=1,N
         PP = PP + RLPSPM2(U(I),IPS,XK)
         RT = U(I)/T
         CALL RLRWEPML(RT,TMP)
         HT = HT + TMP*RT
         CALL RLRWETML(RT,TMP)
         H1T = H1T + TMP
         WGT(I) = DSQRT(TMP)
 25   CONTINUE
      PP = PP/DN
      HT = -HT/DN
      H1T = H1T/DN
      FAC = ZERO
      DO 35 I=1,N
         FAC = FAC + (WGT(I)*WGT(I)*U(I)+HT/PP*RLPSIM2(U(I),
     +        IPS,XK))**2
 35   CONTINUE
      FAC = FAC/DN/(H1T*H1T)
C-----------------------------------------------------------------------
C     STEP 3. COMPUTE THE WEIGHTED LEAST SQUARES
C-----------------------------------------------------------------------
      DO 50 I=1,N
         DO 40 J=1,NP
            SX(I,J) = X(I,J)*WGT(I)
 40      CONTINUE
         SY(I) = Y(I)*WGT(I)
 50   CONTINUE
      INTCH = 1
      IERR = 0
      CALL RLRMTRM2(SX,N,NP,MDX,INTCH,TAU,KK,SF,SG,SH,IP)
      IF (KK .NE. NP) THEN
         IERR = 1
         GOTO 80
      ENDIF
      CALL RLRICLM2(SX,SY,N,NP,MDX,THETA,SH,IP)
      DO 70 I=1,N
         TMP = Y(I)
         DO 60 J=1,NP
            TMP = TMP-X(I,J)*THETA(J)
 60      CONTINUE
         RS(I) = TMP
 70   CONTINUE
 80   RETURN
      END
C=======================================================================
      SUBROUTINE RLFRSTML(X1,X2,Y,N,NP1,NP2,NQ,MDX,T1,T2,X1C,
     +     XTHETA1,XTHETA2,XX,YY,IOPT,INTCH,NREP,TOLS,TOLR,TAU,MAXS1,
     +     ISEED,IERR,SMIN,RS,RS0,SFGH1,SFGH2,IP2,SZ,IT,IPS,XK,BETA,
     +     BET0,ITRACE)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X1(MDX,NP1),X2(MDX,NP2),Y(N),RS(N),RS0(N)
      DIMENSION XX(NQ,NP2),YY(NQ),T1(NP1),T2(NP2),X1C(MDX,NP1)
      DIMENSION XTHETA1(N),XTHETA2(NQ),SFGH1(NP1,3),SFGH2(NP2,3),SZ(N)
      INTEGER IP2(NP2),IT(NQ)
      LOGICAL ALLZERO
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     COMPUTE INITIAL APPROXIMATION FOR M-S ESTIMATE
C-----------------------------------------------------------------------
C     STEP 0: INITIALIZATIONS
C-----------------------------------------------------------------------
      K1=N/2+1
      CONST=BETA*DBLE(N-NP2)
      IF (ITRACE .EQ. 1) THEN
         call mminitclk(iclock)
         ielapse=0
         if (np2 .le. 10) then
            ninc=15000
         else if (np2 .le. 15) then
            ninc=10000
         else
            ninc=5000
         endif
      endif
      NIT=1
      IERR=2
      SMIN=ZERO
      ITYPE=1
      ISIGMA=1

c --- Big Loop --- for(it in 1:nrep) { .... }

C-----------------------------------------------------------------------
C     STEP 1: Draw a subsample of size 'nq' from (X2,Y)
C-----------------------------------------------------------------------
 100  IF (IOPT.NE.3) THEN
         DO 130 K=1,NQ
c     random sample (w/o replacement): rlrndm2(*, U) <=> U <- runif(1)
c     --> ./lmrobmm.f
 110        CALL RLRNDM2(ISEED,RND)
            ITK= idint(RND*N)+1
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
         DO 150 J=1,NP2
            XX(K,J)=X2(ITK,J)
 150     CONTINUE
         YY(K)=Y(ITK)
 160  CONTINUE
C-----------------------------------------------------------------------
C     STEP 2: DECOMPOSE SAMPLE MATRIX AND SOLVE EQUATIONS
C-----------------------------------------------------------------------
      ALLZERO = .TRUE.
      K = 1
      DO WHILE (ALLZERO .AND. K .LE. NQ)
         IF (YY(K) .NE. ZERO) ALLZERO = .FALSE.
         K = K + 1
      ENDDO
      IF (ALLZERO) GOTO 700
      CALL RLRMTRM2(XX,NQ,NP2,NQ,INTCH,TAU,KK,SFGH2(1,1),
     +     SFGH2(1,2),SFGH2(1,3),IP2)
      IF(KK.NE.NP2) GOTO 700
      CALL RLRICLM2(XX,YY,NQ,NP2,NQ,XTHETA2,SFGH2(1,3),IP2)
      DO 210 I=1,N
         S=Y(I)
         DO 200 J=1,NP2
            S=S-XTHETA2(J)*X2(I,J)
 200     CONTINUE
         RS(I)=S
 210  CONTINUE
C-----------------------------------------------------------------------
C     STEP 3: OBTAIN M-ESTIMATE OF B1
C-----------------------------------------------------------------------
      DO 310 I=1,N
         DO 300 J=1,NP1
            X1C(I,J) = X1(I,J)
 300     CONTINUE
 310  CONTINUE
      CALL RLLARSBI(X1C,RS,N,NP1,MDX,N,TAU,NIS,KK,KODE,
     +     S,XTHETA1,RS0,SZ,SFGH1(1,1),SFGH1(1,2),SFGH1(1,3),BET0)
C-----------------------------------------------------------------------
C     STEP 4: COMPUTE THE SCALE ESTIMATE
C-----------------------------------------------------------------------
      IF (SMIN .EQ. ZERO) THEN
         S=1.0D7
         DO 430 I=1,N
            ARI=DABS(RS0(I))
            SZ(I)=ARI
            IF (ARI .NE. ZERO) S=DMIN1(S,ARI)
 430     CONTINUE
         IF (S .EQ. 1.0D7) GOTO 800
         CALL RLSTORM2(SZ,N,K1,S0)
         S0=2.D0*S0
         IF (S0 .EQ. ZERO) S0=S
         SRES=S0
      ENDIF
 435  D=ZERO
      DO 440 I=1,N
         D=D+RLCHIM2(RS0(I)/SRES,IPS,XK)
 440  CONTINUE
      IF (SMIN .NE. ZERO .AND. D .GT. CONST) GOTO 700
      IF (D .LE. CONST) GOTO 500
      S0=1.5D0*S0
      SRES=S0
      GOTO 435
C-----------------------------------------------------------------------
C     STEP 5: SOLVE FOR SRES
C-----------------------------------------------------------------------
 500  CALL RLRSIGM2(RS0,SZ,S0,N,NP2,TOLR,ITYPE,ISIGMA,MAXS1,
     +     NIS,SRES,SZ,SZ,IPS,XK,BETA,BET0)
C-----------------------------------------------------------------------
C     STEP 6: UPDATE BEST FIT
C-----------------------------------------------------------------------
      IERR=0
      SMIN=SRES
      S0=SMIN
      DO 610 K=1,NP1
         T1(K)=XTHETA1(K)
 610  CONTINUE
      DO 620 K=1,NP2
         T2(K)=XTHETA2(K)
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
C     STEP 8: SOLVE SYSTEM OF EQUATIONS FOR THETA AND SRES
C-----------------------------------------------------------------------
 800  IF (IERR.EQ.2) THEN
         DO 810 K=1,NP1
            T1(K)=XTHETA1(K)
 810     CONTINUE
         DO 820 K=1,NP2
            T2(K)=XTHETA2(K)
 820     CONTINUE
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE RLLUDCM2(X,N,IDX,W,IFAIL)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N,N),W(N),IDX(N)
      DATA ZERO,ONE,EPS/0.D0,1.D0,2.220446e-16/
C-----------------------------------------------------------------------
C     LU DECOMPOSITION OF A SQUARE MATRIX
C-----------------------------------------------------------------------
      DO I=1,N
         TMP=ZERO
         DO J=1,N
            TEMP=DABS(X(I,J))
            IF (TEMP .GT. TMP) TMP=TEMP
         ENDDO
         IF (TMP .EQ. ZERO) THEN
            IFAIL=1
            RETURN
         ENDIF
         W(I)=ONE/TMP
      ENDDO
      DO J=1,N
         DO I=1,J-1
            TMP=X(I,J)
            DO K=1,I-1
               TMP=TMP-X(I,K)*X(K,J)
            ENDDO
            X(I,J)=TMP
         ENDDO
         TMP=ZERO
         DO I=J,N
            TEMP=X(I,J)
            DO K=1,J-1
               TEMP=TEMP-X(I,K)*X(K,J)
            ENDDO
            X(I,J)=TEMP
            TEMP=W(I)*DABS(TEMP)
            IF (TEMP .GE. TMP) THEN
               TMP=TEMP
               IMAX=I
            ENDIF
         ENDDO
         IF (J .NE. IMAX) THEN
            DO K=1,N
               TEMP=X(IMAX,K)
               X(IMAX,K)=X(J,K)
               X(J,K)=TEMP
            ENDDO
            W(IMAX)=W(J)
         ENDIF
         IDX(J)=IMAX
         IF (DABS(X(J,J)) .LE. EPS) THEN
            IFAIL=1
            RETURN
         ENDIF
         IF (J .NE. N) THEN
            TMP=ONE/X(J,J)
            DO I=J+1,N
               X(I,J)=X(I,J)*TMP
            ENDDO
         ENDIF
      ENDDO
      RETURN
      END
C=======================================================================
      SUBROUTINE RLLUSLM2(X,N,IDX,Y)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N,N),Y(N),IDX(N)
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     SOLVE A LU DECOMPOSED SQUARE MATRIX
C-----------------------------------------------------------------------
      II=-1
      DO I=1,N
         IP=IDX(I)
         TMP=Y(IP)
         Y(IP)=Y(I)
         IF (II .GE. 0) THEN
            DO J=II,I-1
               TMP=TMP-X(I,J)*Y(J)
            ENDDO
         ELSE IF(TMP .GT. ZERO) THEN
            II=I
         ENDIF
         Y(I)=TMP
      ENDDO
      DO I=N,1,-1
         TMP=Y(I)
         DO J=I+1,N
            TMP=TMP-X(I,J)*Y(J)
         ENDDO
         Y(I)=TMP/X(I,I)
      ENDDO
      RETURN
      END
C=======================================================================
      SUBROUTINE RLLUINM2(X,X1,N,IDX,W,IFAIL)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N,N),X1(N,N),W(N),IDX(N)
      DATA ZERO,ONE/0.D0,1.D0/
C-----------------------------------------------------------------------
C     MATRIX INVERSION USING LU DECOMPOSITION
C-----------------------------------------------------------------------
      DO I=1,N
         DO J=1,N
            X1(I,J)=X(I,J)
         ENDDO
      ENDDO
      IFAIL=0
      CALL RLLUDCM2(X1,N,IDX,W,IFAIL)
      IF (IFAIL .EQ. 1) RETURN
      DO J=1,N
         DO I=1,N
            W(I)=ZERO
         ENDDO
         W(J)=ONE
         CALL RLLUSLM2(X1,N,IDX,W)
         DO I=1,N
            X(I,J)=W(I)
         ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
      SUBROUTINE RLINVSM2(X,N,IDX)
C.......................................................................
C     INVERSE A SYMMETRIC MATRIX USING CHOLESKY DECOMPOSITION
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N,N)
      DATA ZERO,ONE/0.D0,1.D0/
      DO I=1,N
         DO J=1,I
            SUM=X(J,I)
            DO K=1,J-1
               SUM=SUM-X(J,K)*X(I,K)
            ENDDO
            IF (I .EQ. J) THEN
               IF (SUM .LE. ZERO) THEN
                  IDX=1
                  RETURN
               ENDIF
               X(I,J)=DSQRT(SUM)
            ELSE
               X(I,J)=SUM/X(J,J)
            ENDIF
         ENDDO
      ENDDO
      DO I=1,N
         X(I,I)=ONE/X(I,I)
         DO J=I+1,N
            SUM=ZERO
            DO K=I,J-1
               SUM=SUM-X(J,K)*X(K,I)
            ENDDO
            X(J,I)=SUM/X(J,J)
         ENDDO
      ENDDO
      DO I=1,N
         DO J=I,N
            SUM=ZERO
            DO K=J,N
               SUM=SUM+X(K,I)*X(K,J)
            ENDDO
            X(J,I)=SUM
         ENDDO
         DO J=1,I-1
            X(J,I)=X(I,J)
         ENDDO
      ENDDO
      RETURN
      END
C=======================================================================
      FUNCTION RLCOVGM2(X,MDX,N,NP,DELTA,SIGMA)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),DELTA(NP)
      DATA ZERO/0.D0/
      RLCOVGM2=ZERO
      DO I=1,N
         TMP=ZERO
         DO J=1,NP
            TMP=TMP+X(I,J)*DELTA(J)
         ENDDO
         TMP=TMP/SIGMA
         IF (TMP .GT. RLCOVGM2) RLCOVGM2=TMP
      ENDDO
      RETURN
      END
C=======================================================================
      SUBROUTINE RLYWAGM2(X,Y,THETA,SIGMA,N,NP,MDX,
     +     TOL,GAM,TAU,MAXIT,NIT,RS,DELTA,SC,SF,SG,SH,IP,SX)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(NP),RS(N)
      DIMENSION DELTA(NP),SC(N),SF(NP),SG(NP),SH(NP),SX(MDX,NP)
      INTEGER IP(NP)
      DATA ZERO,ONE/0.D0,1.D0/
C-----------------------------------------------------------------------
C     W-ALGORITHM FOR ROBUST AND BOUNDED INFLUENCE LINEAR REGRESSION
C-----------------------------------------------------------------------
      MDXP1=MDX+1
      LDIAG=MIN0(N,NP)
      XK=ONE
      INTCH=1
C-----------------------------------------------------------------------
C     STEP 1. SET NIT := 1
C-----------------------------------------------------------------------
      NIT=1
C-----------------------------------------------------------------------
C     STEP 2. COMPUTE RLRESDM2 AS R=Y-X*THETA
C-----------------------------------------------------------------------
 200  CALL RLRESDM2(X,Y,THETA,N,NP,MDX,RS)
C-----------------------------------------------------------------------
C     STEP 3. COMPUTE WEIGHTS AND APPLY THEM TO X; STORE RESULT IN SX
C-----------------------------------------------------------------------
      DO 430 I=1,N
         SC(I)=ONE
         IF (RS(I) .EQ. ZERO) GOTO 410
         T=RS(I)/SIGMA
         SC(I)=RLPSIM2(T,4,XK)/T
 410     PI=DSQRT(SC(I))
         RS(I)=PI*RS(I)
         DO 420 J=1,NP
            SX(I,J)=PI*X(I,J)
 420     CONTINUE
 430  CONTINUE
C-----------------------------------------------------------------------
C     STEP 4. SOLVE FOR DELTA
C-----------------------------------------------------------------------
      CALL RLRMTRM2(SX,N,NP,MDX,INTCH,TAU,K,SF,SG,SH,IP)
      IF (K.EQ.0) RETURN
      KK=MDX*(K-1)+K
      IF (K.NE.NP) CALL RLSWAPM2(SX,SF,K,MDXP1,1,KK,K)
      DO 500 JJ=1,LDIAG
         J=JJ
         CALL RLH12M2(2,J,J+1,N,SX(1,J),1,SH(J),RS,1,N,1,N)
 500  CONTINUE
      IF (K.NE.NP) CALL RLSWAPM2(SX,SF,K,MDXP1,1,KK,K)
      CALL RLSOLVM2(SX,RS,NP,K,MDX,N)
      IF (K.EQ.NP) GOTO 530
      KP1=K+1
      DO 510 J=KP1,NP
         RS(J)=ZERO
 510  CONTINUE
      DO 520 J=1,K
         I=J
         CALL RLH12M2(2,I,KP1,NP,SX(I,1),MDX,SG(I),RS,1,N,1,NP)
 520  CONTINUE
 530  DO 540 J=1,NP
         DELTA(J)=GAM*RS(J)
 540  CONTINUE
      CALL RLPERMM2(DELTA,IP,LDIAG,NP)
C-----------------------------------------------------------------------
C     STEP 6. COMPUTE NEW SOLUTION
C-----------------------------------------------------------------------
      DO 600 J=1,NP
         THETA(J)=THETA(J)+DELTA(J)
 600  CONTINUE
C-----------------------------------------------------------------------
C     STEP 7. STOP ITERATIONS IF DESIRED PRECISION IS REACHED
C-----------------------------------------------------------------------
      IF (NIT.EQ.MAXIT) GOTO 800
      IF(RLCOVGM2(X,MDX,N,NP,DELTA,SIGMA) .LE. TOL) GOTO 800
      NIT=NIT+1
      GOTO 200
 800  CALL RLRESDM2(X,Y,THETA,N,NP,MDX,RS)
      RETURN
      END
C=======================================================================
      SUBROUTINE RLBETAM2(X1,X2,Y,N,NP1,NP2,MDX,S0,S1,B1,B2,T1,T2,RS,
     +     RSTMP,TOLR,TAU,MAXIT,MAXS1,SFGH,IPS,XK,BETA,BET0,
     +     IFAIL,UV,A,B,CC,C2,D,BD,H,TC,X1C,IP,IDX)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X1(MDX,NP1),X2(MDX,NP2),X1C(MDX,NP1),Y(N)
      DIMENSION B1(NP1),T1(NP1)
      DIMENSION B2(NP2),T2(NP2)
      DIMENSION RS(N),RSTMP(N),SFGH(NP1,3)
      DIMENSION UV(N),A(NP1,NP2),B(NP1,NP1),BD(NP1,NP2),D(N,NP2)
      DIMENSION CC(NP2,NP2),C2(NP2,NP2),H(NP2)
      INTEGER IP(NP1),IDX(NP2)
      DATA ZERO,ONE/0.D0,1.D0/
C-----------------------------------------------------------------------
C     STEP 1. COMPUTES WEIGHTS U AND MATRICES A, B AND D
C-----------------------------------------------------------------------
      DO I=1,N
         TMP=RS(I)/(TC*S0)
         RS(I)=TMP
         UV(I)=RLPSPM2(TMP,4,ONE)
         TMP=Y(I)
         DO J=1,NP1
            TMP=TMP-X1(I,J)*B1(J)
         ENDDO
         RSTMP(I)=TMP
      ENDDO
      DO I=1,NP1
         DO J=1,NP2
            TMP=ZERO
            DO K=1,N
               TMP=TMP+UV(K)*X1(K,I)*X2(K,J)
            ENDDO
            A(I,J)=TMP
         ENDDO
         DO J=1,NP1
            TMP=ZERO
            DO K=1,N
               TMP=TMP+UV(K)*X1(K,I)*X1(K,J)
            ENDDO
            B(I,J)=TMP
         ENDDO
      ENDDO
      IFAIL=0
      CALL RLINVSM2(B,NP1,IFAIL)
      IF (IFAIL .EQ. 1) RETURN
      DO I=1,NP1
         DO J=1,NP2
            TMP=ZERO
            DO K=1,NP1
               TMP=TMP+B(I,K)*A(K,J)
            ENDDO
            BD(I,J)=TMP
         ENDDO
      ENDDO
      DO I=1,N
         DO J=1,NP2
            TMP=X2(I,J)
            DO K=1,NP1
               TMP=TMP-X1(I,K)*BD(K,J)
            ENDDO
            D(I,J)=TMP
         ENDDO
      ENDDO
C-----------------------------------------------------------------------
C     STEP 2. COMPUTES WEIGHTS V AND MATRIX C AND H
C-----------------------------------------------------------------------
      DO 100 I=1,N
         UV(I)=ONE
         TMP=RS(I)
         IF (TMP .EQ. ZERO) GOTO 100
         UV(I)=RLPSIM2(TMP,IPS,XK)/TMP
 100  CONTINUE
      DO I=1,NP2
         DO J=1,NP2
            TMP=ZERO
            DO K=1,N
               TMP=TMP+UV(K)*D(K,I)*X2(K,J)
            ENDDO
            CC(I,J)=TMP
         ENDDO
         TMP=ZERO
         DO K=1,N
            TMP=TMP+D(K,I)*UV(K)*RSTMP(K)
         ENDDO
         H(I)=TMP
      ENDDO
C-----------------------------------------------------------------------
C     STEP 3. UPDATES B1 AND B2 TO T1 AND T2
C-----------------------------------------------------------------------
      CALL RLLUINM2(CC,C2,NP2,IDX,B2,IFAIL)
      IF (IFAIL .EQ. 1) RETURN
      DO I=1,NP2
         TMP=ZERO
         DO J=1,NP2
            TMP=TMP+CC(I,J)*H(J)
         ENDDO
         T2(I)=TMP
      ENDDO
      CALL RLRESDM2(X2,Y,T2,N,NP2,MDX,RSTMP)
      DO I=1,NP1
         T1(I)=B1(I)
      ENDDO
      CALL RLYWAGM2(X1,RSTMP,T1,S0,N,NP1,MDX,TOLR,ONE,TAU,MAXIT,
     +     NIS,RS,B1,UV,SFGH(1,1),SFGH(1,2),SFGH(1,3),IP,X1C)
C-----------------------------------------------------------------------
C     STEP 4. UPDATES RESIDUALS AND SCALE
C-----------------------------------------------------------------------
      CALL RLRESDM2(X1,RSTMP,T1,N,NP1,MDX,RS)
      NP=NP1+NP2
      ITYPE=1
      ISIGMA=1
      CALL RLRSIGM2(RS,UV,S0,N,NP,TOLR,ITYPE,ISIGMA,MAXS1,
     +     NIS,S1,UV,UV,IPS,XK,BETA,BET0)
      RETURN
      END
C=======================================================================
      SUBROUTINE RLDSCNM2(X1,X2,Y,N,NP1,NP2,MDX,S0,S1,B1,B2,T1,
     +     T2,RS,RSTMP,TOLR,TAU,MAXIT,MAXS1,SFGH,IPS,XK,BETA,BET0,
     +     IFAIL,UV,A,B,CC,C2,D,BD,H,TC,X1C,IP,IDX,WP1,WP2,NIT,MAXK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X1(MDX,NP1),X2(MDX,NP2),X1C(MDX,NP1),Y(N)
      DIMENSION B1(NP1),T1(NP1)
      DIMENSION B2(NP2),T2(NP2)
      DIMENSION RS(N),RSTMP(N),SFGH(NP1,3),WP1(NP1),WP2(NP2)
      DIMENSION UV(N),A(NP1,NP2),B(NP1,NP1),BD(NP1,NP2),D(N,NP2)
      DIMENSION CC(NP2,NP2),C2(NP2,NP2),H(NP2)
      INTEGER IP(NP1),IDX(NP2)
      DATA ONE,TWO/1.D0,2.D0/
C-----------------------------------------------------------------------
C     STEP 1. SET NIT := 1
C-----------------------------------------------------------------------
      NIT=1
      NP=NP1+NP2
      ITYPE=1
      ISIGMA=1
      ERR=20.D0
C-----------------------------------------------------------------------
C     STEP 2. SET NIT := 1
C-----------------------------------------------------------------------
 200  DO I=1,NP1
         WP1(I)=B1(I)
      ENDDO
      DO I=1,NP2
         WP2(I)=B2(I)
      ENDDO
      CALL RLBETAM2(X1,X2,Y,N,NP1,NP2,MDX,S0,S1,WP1,WP2,T1,T2,RS,
     +     RSTMP,TOLR,TAU,MAXIT,MAXS1,SFGH,IPS,XK,BETA,BET0,
     +     IFAIL,UV,A,B,CC,C2,D,BD,H,TC,X1C,IP,IDX)
      IF (IFAIL .EQ. 1) GOTO 800
      IF (NIT .GE. MAXIT) GOTO 800
      IF (ERR .LE. TOLR) GOTO 800
      IF (S1 .GT. S0) THEN
         DO I=1,NP2
            H(I)=T2(I)-B2(I)
         ENDDO
         K=1
 300     DO I=1,NP2
            H(I)=H(I)/TWO
            T2(I)=B2(I)+H(I)
         ENDDO
         CALL RLRESDM2(X2,Y,T2,N,NP2,MDX,RS)
         DO I=1,NP1
            T1(I)=B1(I)
         ENDDO
         CALL RLYWAGM2(X1,RS,T1,S0,N,NP1,MDX,TOLR,ONE,TAU,MAXIT,
     +        NIS,RSTMP,WP1,UV,SFGH(1,1),SFGH(1,2),SFGH(1,3),IP,X1C)
         CALL RLRESDM2(X1,RSTMP,T1,N,NP1,MDX,RS)
         CALL RLRSIGM2(RS,UV,S0,N,NP,TOLR,ITYPE,ISIGMA,MAXS1,
     +        NIS,S1,UV,UV,IPS,XK,BETA,BET0)
         IF (S1 .LT. S0) GOTO 600
         IF (K .GE. MAXK) GOTO 600
         K=K+1
         GOTO 300
      ENDIF
 600  ERR=S0/S1-ONE
      DO I=1,NP1
         B1(I)=T1(I)
      ENDDO
      DO I=1,NP2
         B2(I)=T2(I)
      ENDDO
      S0=S1
      GOTO 200
 800  RETURN
      END

