C=======================================================================
C     ROBUST MM-ESTIMATION
C     MATHSOFT, INC.
C     08/09/00
C=======================================================================
c
c All subroutine/function names originally were
c "RL....M2" : [R]obust [L]ibrary ... "MM"-Estimation
c  ~~    ~~
c           change ----> only capitalize the information-content part
c
      SUBROUTINE rlEXCHm2(S,N,NN,H,K)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(NN)
      INTEGER H
C-----------------------------------------------------------------------
c  EXCHange  S(..H ..) with S(.. K ..)  :
c  -------   auxiliary only called from  rlKFASm2()  below
C-----------------------------------------------------------------------

      LH=H*(H+1)/2
      LK=K*(K+1)/2
      T=S(LH)
      S(LH)=S(LK)
      S(LK)=T
      LH=LH-H
      LK=LK-K
      M=H-1
      IF (M.EQ.0) GOTO 15
      DO 10 I=1,M
         LH=LH+1
         LK=LK+1
         T=S(LH)
         S(LH)=S(LK)
         S(LK)=T
 10   CONTINUE
 15   LH=LH+1
      LK=LK+1
      M=K-H-1
      IF (M.EQ.0) GOTO 30
      DO 20 I=1,M
         LH=LH+H-1+I
         LK=LK+1
         T=S(LH)
         S(LH)=S(LK)
         S(LK)=T
 20   CONTINUE
 30   LH=LH+K-1
      LK=LK+1
      M=N-K
      IF (M.EQ.0) GOTO 45
      DO 40 I=1,M
         LH=LH+K+I-1
         LK=LK+K+I-1
         T=S(LH)
         S(LH)=S(LK)
         S(LK)=T
 40   CONTINUE
 45   RETURN
      END
C=======================================================================
      SUBROUTINE rlKFASm2(XT,COV,K,NP,MDX,NCOV,F,SE,SG,IP)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION COV(NCOV),XT(MDX,NP),SE(NP),SG(NP)
      INTEGER IP(NP)
C-----------------------------------------------------------------------
      KP1=K+1
      LDIAG=MIN0(MDX,NP)
C-----------------------------------------------------------------------
C     TRANSFORM UNSCALED COVARIANCE MATRIX TO COMPENSATE
C     HOUSEHOLDER TRANSFORMATIONS AND PERMUTATIONS
c
c Called from S /R 's lmRob.ucovcoef()  {and from [g]lmRob's Fortran}
C-----------------------------------------------------------------------
      IF (K.EQ.NP) GOTO 130
      DO 120 II=1,K
         I=II
         CALL rlVSVm2(I,KP1,NP,XT(I,1),MDX,SG(I),COV,NCOV,SE)
 120  CONTINUE
 130  CONTINUE
      DO 150 JJ=1,LDIAG
         J=LDIAG-JJ+1
         IF (IP(J).EQ.J) GOTO 150
         L=IP(J)
         CALL rlEXCHm2(COV,NP,NCOV,J,L)
 150  CONTINUE
C-----------------------------------------------------------------------
C     MULTIPLY COV BY THE SCALE FACTOR F
C-----------------------------------------------------------------------
      IF (F .GT. 0.D0) CALL rlSCALm2(COV,F,NCOV,1,NCOV)
      RETURN
      END

C=======================================================================
      SUBROUTINE rlKIASm2(XT,K,NP,MDX,NCOV,FU,FB,COV)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XT(MDX,NP),COV(NCOV)
      DATA DZERO/0.D0/
C-----------------------------------------------------------------------
      LDIAG=MIN0(MDX,NP)
      KP1=K+1
      L=0
      DO 20 J=1,K
         DO 10 I=1,J
            L=L+1
            COV(L)=XT(I,J)
 10      CONTINUE
 20   CONTINUE
C-----------------------------------------------------------------------
C     INVERT U UPON ITSELF
C-----------------------------------------------------------------------
      DO 30 J=1,K
         XT(J,J) = 1.D0/XT(J,J)
 30   CONTINUE
      IF (K.EQ.1) GOTO 60
      KM1=K-1
      DO 55 I=1,KM1
         IP1=I+1
         DO 50 J=IP1,K
            JM1=J-1
            SM=DZERO
            DO 40 L=I,JM1
               SM=SM+XT(I,L)*XT(L,J)
 40         CONTINUE
            XT(I,J)=-SM*XT(J,J)
 50      CONTINUE
 55   CONTINUE
C-----------------------------------------------------------------------
C     REPLACE U**(-1) BY UPPER TRIANG.PART OF (U*U**T)**(-1)
C-----------------------------------------------------------------------
 60   CONTINUE
      DO 85 I=1,K
         DO 80 J=I,K
            SM=DZERO
            DO 70 L=J,K
               SM=SM+XT(I,L)*XT(J,L)
 70         CONTINUE
            XT(I,J)=SM
 80      CONTINUE
 85   CONTINUE
C-----------------------------------------------------------------------
C     INTERCH. (U*U**T)**(-1) WITH COV(1)...COV(K*(K+1)/2)
C-----------------------------------------------------------------------
      L=0
      DO 100 J=1,K
         DO 90 I=1,J
            L=L+1
            AIJ=XT(I,J)
            XT(I,J)=COV(L)
            COV(L)=AIJ
 90      CONTINUE
 100  CONTINUE
C-----------------------------------------------------------------------
C     MULTIPLY COV BY THE SCALE FACTOR FU
C-----------------------------------------------------------------------
      IF (FU .GT. DZERO) CALL rlSCALm2(COV,FU,NCOV,1,NCOV)
C-----------------------------------------------------------------------
C     COMPLETE COV
C-----------------------------------------------------------------------
      IF (K.EQ.NP) RETURN
      II=K*(K+1)/2+1
      J=KP1
      L=II+K
      DO 160 I=II,NCOV
         COV(I)=DZERO
         IF (I.NE.L) GOTO 160
         COV(I)=FB
         J=J+1
         L=L+J
 160  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlMCHLm2(A,N,NN,INFO)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NN)
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     CHOLESKY DECOMPOSITION OF A SYMMETRIC MATRIX
C-----------------------------------------------------------------------
      JJ=0
      DO 30 J=1,N
         INFO=J
         S=ZERO
         JM1=J-1
         KJ=JJ
         KK=0
         IF (JM1.LT.1) GOTO 20
         DO 10 K=1,JM1
            KJ=KJ+1
            CALL rlDOTPm2(A(KK+1),A(JJ+1),K-1,1,1,NN-KK,NN-JJ,DTP)
            T=A(KJ)-DTP
            KK=KK+K
            T=T/A(KK)
            A(KJ)=T
            S=S+T*T
 10      CONTINUE
 20      CONTINUE
         JJ=JJ+J
         S=A(JJ)-S
         IF (S .LE. ZERO) GOTO 40
         A(JJ)=DSQRT(S)
 30   CONTINUE
      INFO=0
 40   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlSCALm2(X,SA,N,INCX,MDX)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX)
C-----------------------------------------------------------------------
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1) GOTO 20
C-----------------------------------------------------------------------
C     CODE FOR INCREMENT NOT EQUAL TO 1
C-----------------------------------------------------------------------
      NINCX=N*INCX
      DO 10 I=1,NINCX,INCX
         X(I)=SA*X(I)
 10   CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     CODE FOR INCREMENT EQUAL TO 1
C-----------------------------------------------------------------------
 20   M=MOD(N,5)
      IF (M.EQ.0) GOTO 40
      DO 30 I=1,M
         X(I)=SA*X(I)
 30   CONTINUE
      IF (N.LT.5) RETURN
 40   MP1=M+1
      DO 50 I=MP1,N,5
         X(I)=SA*X(I)
         X(I+1)=SA*X(I+1)
         X(I+2)=SA*X(I+2)
         X(I+3)=SA*X(I+3)
         X(I+4)=SA*X(I+4)
 50   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlSWAPm2(X,Y,N,INCX,INCY,MDX,MDY)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX),Y(MDY)
C-----------------------------------------------------------------------
      IF (N.EQ.0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1) GOTO 20
C-----------------------------------------------------------------------
C     CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT
C     EQUAL TO 1
C-----------------------------------------------------------------------
      IX=1
      IY=1
      IF (INCX.LT.0) IX=(-N+1)*INCX+1
      IF (INCY.LT.0) IY=(-N+1)*INCY+1
      DO 10 I=1,N
         TEMP=X(IX)
         X(IX)=Y(IY)
         Y(IY)=TEMP
         IX=IX+INCX
         IY=IY+INCY
 10   CONTINUE
      RETURN
C-----------------------------------------------------------------------
C     CODE FOR BOTH INCREMENTS EQUAL TO 1
C-----------------------------------------------------------------------
 20   M=MOD(N,3)
      IF (M.EQ.0) GOTO 40
      DO 30 I=1,M
         TEMP=X(I)
         X(I)=Y(I)
         Y(I)=TEMP
 30   CONTINUE
      IF (N.LT.3) RETURN
 40   MP1=M+1
      DO 50 I=MP1,N,3
         TEMP=X(I)
         X(I)=Y(I)
         Y(I)=TEMP
         TEMP=X(I+1)
         X(I+1)=Y(I+1)
         Y(I+1)=TEMP
         TEMP=X(I+2)
         X(I+2)=Y(I+2)
         Y(I+2)=TEMP
 50   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlMINVm2(R,N,NN,TAU,ISING)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(NN)
      DATA DZERO,ONE/0.D0,1.D0/
C-----------------------------------------------------------------------
C     INVERTS A TRIANGULAR MATRIX
C-----------------------------------------------------------------------
      ISING=0
      I1=0
      DO 10 I=1,N
         I1=I1+I
         IF (DABS(R(I1)) .LE. TAU) GOTO 900
         R(I1)=ONE/R(I1)
 10   CONTINUE
      IF (N.EQ.1) RETURN
      I1=0
      NM1=N-1
      DO 40 I=1,NM1
         I1=I1+I
         J1=I1+I
         IP1=I+1
         DO 30 J=IP1,N
            SM=DZERO
            IL=I1
            LJ=J1
            JM1=J-1
            DO 20 L=I,JM1
               SM=SM+R(IL)*R(LJ)
               LJ=LJ+1
               IL=IL+L
 20         CONTINUE
            R(J1)=-R(LJ)*SM
            J1=J1+J
 30      CONTINUE
 40   CONTINUE
      RETURN
 900  ISING=1
      RETURN
      END
C=======================================================================
      SUBROUTINE rlMTT1m2(A,B,N,NN)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NN),B(NN)
      DATA DZERO/0.D0/
C-----------------------------------------------------------------------
C     MULTIPLIES AN UPPER TRIANGULAR MATRIX BY ITS TRANSPOSE
C-----------------------------------------------------------------------
      IJ=0
      JJ=0
      DO 30 J=1,N
         DO 20 I=1,J
            IJ=IJ+1
            SM=DZERO
            IL=JJ+I
            JL=JJ+J
            DO 10 L=J,N
               SM=SM+A(IL)*A(JL)
               IL=IL+L
               JL=JL+L
 10         CONTINUE
            B(IJ)=SM
 20      CONTINUE
         JJ=JJ+J
 30   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlSTORm2(Y,N,J,YJ)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N)
C-----------------------------------------------------------------------
C     rlSTORm2 SEARCHES THE J-TH VALUE IN ORDER OF MAGNITUDE IN
C     A VECTOR OF LENGTH N.
C-----------------------------------------------------------------------
      L=1
      LR=N
 20   IF (L.GE.LR) GOTO 90
      AX=Y(J)
      JNC=L
      JJ=LR
 30   IF(JNC.GT.JJ) GOTO 80
 40   IF (Y(JNC).GE.AX) GOTO 50
      JNC=JNC+1
      GOTO 40
 50   IF(Y(JJ).LE.AX) GOTO 60
      JJ=JJ-1
      GOTO 50
 60   IF(JNC.GT.JJ) GOTO 70
      WA=Y(JNC)
      Y(JNC)=Y(JJ)
      Y(JJ)=WA
      JNC=JNC+1
      JJ=JJ-1
 70   GOTO 30
 80   IF(JJ.LT.J) L=JNC
      IF(J.LT.JNC) LR=JJ
      GOTO 20
 90   YJ=Y(J)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlH12m2(MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV,MDC)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(IUE,M),C(MDC)
      DATA ONE/1.D0/
C-----------------------------------------------------------------------
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN
      CL=DABS(U(1,LPIVOT))
      IF (MODE.EQ.2) GOTO 60
C-----------------------------------------------------------------------
C     CONSTRUCT THE TRANSFORMATION
C-----------------------------------------------------------------------
      DO 10 J=L1,M
         CL=DMAX1(DABS(U(1,J)),CL)
 10   CONTINUE
C      IF (CL) 130,130,20
      if(CL > 0.D0) then
        goto 20
      else
        goto 130
      end if
 20   CLINV=ONE/CL
      SM=(U(1,LPIVOT)*CLINV)**2
      DO 30 J=L1,M
         SM=SM+(U(1,J)*CLINV)**2
 30   CONTINUE
C-----------------------------------------------------------------------
C     CONVERT DBLE. PRE. SM TO SNGL. PREC. SM1
C-----------------------------------------------------------------------
      SM1=SM
      CL=CL*DSQRT(SM1)
C      IF (U(1,LPIVOT)) 50,50,40
      if(U(1,LPIVOT) > 0.D0) then
        goto 40
      else
        goto 50
      end if
 40   CL=-CL
 50   UP=U(1,LPIVOT)-CL
      U(1,LPIVOT)=CL
      GOTO 70
C-----------------------------------------------------------------------
C     APPLY THE TRANSFORMATION I+U*(U**T)/B TO C
C-----------------------------------------------------------------------
c 60   IF (CL) 130,130,70
 60   if(CL > 0.D0) then
        goto 70
      else
        goto 130
      end if
 70   IF (NCV.LE.0) RETURN
      B=UP*U(1,LPIVOT)
C-----------------------------------------------------------------------
C     B MUST BE NONPOSITIVE HERE. IF B=0., RETURN.
C-----------------------------------------------------------------------
c      IF (B) 80,130,130
      if(B < 0.D0) then
        goto 80
      else
        goto 130
      end if
 80   B=ONE/B
      I2=1-ICV+ICE*(LPIVOT-1)
      INCR=ICE*(L1-LPIVOT)
      DO 120 J=1,NCV
         I2=I2+ICV
         I3=I2+INCR
         I4=I3
         SM=C(I2)*UP
         DO 90 I=L1,M
            SM=SM+C(I3)*U(1,I)
            I3=I3+ICE
 90      CONTINUE
C         IF (SM) 100,120,100
      if(SM > 0.D0 .or. SM < 0.D0) then
        goto 100
      else
        goto 120
      end if
 100     SM=SM*B
         C(I2)=C(I2)+SM*UP
         DO 110 I=L1,M
            C(I4)=C(I4)+SM*U(1,I)
            I4=I4+ICE
 110     CONTINUE
 120  CONTINUE
 130  RETURN
      END
C=======================================================================
      SUBROUTINE rlXSYm2(X,Y,S,N,NN,RESULT)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N),S(NN)
C-----------------------------------------------------------------------
      SM=0.D0
      L=0
      DO 25 I=1,N
         L=L+I
         L1=L-I+1
         K=0
         DO 20 J=L1,L
            K=K+1
            IF (J.EQ.L) GOTO 10
            SM=SM+S(J)*(X(I)*Y(K)+X(K)*Y(I))
            GOTO 20
 10         SM=SM+S(J)*X(I)*Y(I)
 20      CONTINUE
 25   CONTINUE
      RESULT=SM
      RETURN
      END
C=======================================================================
      SUBROUTINE rlRNDm2(ISEED,RN)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION T(128)
      DATA INIT,T/0,128*0.D0/
C-----------------------------------------------------------------------
C     RANDOM NUMBER GENERATOR ACCORDING TO THE LINEAR CONGRUENT SCHEME
C                  ISEED=ISEED*5761+999 MODULO 65536
C     IMPROVED AFTER MACLAREN-MARSAGLIA
C-----------------------------------------------------------------------
      IF (INIT.EQ. 0 .OR. INIT.NE. ISEED) THEN
         ISEED=ISEED-(ISEED/65536)*65336
         DO 100 I=1,128
            ISEED=ISEED*5761+999
            ISEED=ISEED-(ISEED/65536)*65536
            T(I)=DBLE(ISEED)/65536.0D0
 100     CONTINUE
      ENDIF
      ISEED=ISEED*5761+999
      ISEED=ISEED-(ISEED/65536)*65536
      I=128*ISEED/65536
      RN=T(I+1)
      ISEED=ISEED*5761+999
      ISEED=ISEED-(ISEED/65536)*65536
      T(I+1)=DBLE(ISEED)/65536.0D0
      INIT=ISEED
      RETURN
      END
C=======================================================================
      SUBROUTINE rlNCOMm2(N,NP,IT)
C.......................................................................
      DIMENSION IT(NP)
C-----------------------------------------------------------------------
C     COMPUTE ALL COMBINATIONS FOR RESAMPLING ALGORITHM
c     i.e., given previous  it[1:np] in {1,..,n} -- compute *next* one
C-----------------------------------------------------------------------
      IN=NP
 10   IT(IN)=IT(IN)+1
      IF(IT(IN).GT.N-NP+IN) THEN
         IN=IN-1
         GOTO 10
      ENDIF
      IF(IN.NE.NP) THEN
         DO 20 I=IN+1,NP
            IT(I)=IT(I-1)+1
 20      CONTINUE
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE rlRICLm2(XT,Y,N,NP,MDXT,THETA,SH,SP)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XT(MDXT,NP),Y(N),THETA(MDXT),SH(NP)
      INTEGER SP(NP)
C-----------------------------------------------------------------------
      DO 20 JJ=1,NP
         J=JJ
         CALL rlH12m2(2,J,J+1,N,XT(1,J),1,SH(J),Y,1,N,1,N)
 20   CONTINUE
C-----------------------------------------------------------------------
C     SOLVE THE SYSTEM
C-----------------------------------------------------------------------
      DO 30 I=1,N
         THETA(I)=Y(I)
 30   CONTINUE
      CALL rlSOLVm2(XT,THETA,NP,NP,MDXT,N)
C-----------------------------------------------------------------------
C     TRANSFORM THE SOLUTION VECTOR FOR OUTPUT
C-----------------------------------------------------------------------
      CALL rlPERMm2(THETA,SP,NP,NP)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlNSIGm2(RS,WGT,WGT2,SIGMA,SIGMB,N,ITYPE,IPS,XK,CONST)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RS(N),WGT(N),WGT2(N)
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     COMPUTES A NEW VALUE SIGMB FOR THE ROBUST ESTIMATE OF THE
C     ERROR STANDARD DEVIATION IN THE HUBER'S ALGORITHM FOR REGRESSION.
C     rlNSIGm2 CALLS THE FUNCTION rlCHIm2.
C-----------------------------------------------------------------------
      TMP=ZERO
      IF (ITYPE.NE.1) GOTO 20
C-----------------------------------------------------------------------
C     HUBER-TYPE
C-----------------------------------------------------------------------
      DO 10 I=1,N
         S=RS(I)/SIGMA
         TMP=TMP+rlCHIm2(S,IPS,XK)
 10   CONTINUE
      GOTO 90
C-----------------------------------------------------------------------
C     MALLOWS-TYPE
C-----------------------------------------------------------------------
 20   IF (ITYPE.NE.2) GOTO 40
      DO 30 I=1,N
         S=RS(I)/SIGMA
         IF (WGT(I) .LE. ZERO) GOTO 30
         TMP=TMP+rlCHIm2(S,IPS,XK)*WGT(I)
 30   CONTINUE
      GOTO 90
C-----------------------------------------------------------------------
C     SCHWEPPE-TYPE
C-----------------------------------------------------------------------
 40   DO 50 I=1,N
         SW=SIGMA*WGT(I)
         IF (SW .EQ. ZERO .OR. WGT(I) .LE. ZERO) GOTO 50
         S=RS(I)/SW
         TMP=TMP+rlCHIm2(S,IPS,XK)*WGT2(I)
 50   CONTINUE
 90   SIGMB=DSQRT(TMP/CONST)*SIGMA
      RETURN
      END
C=======================================================================
      FUNCTION rlISIGm2(SIGMA,SIGMB,TOL)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER rlISIGm2
C-----------------------------------------------------------------------
      DS=DABS(SIGMA-SIGMB)/DMAX1(1.D0,SIGMA)
      rlISIGm2=0
      IF (TOL.GE.DS) rlISIGm2=1
      RETURN
      END
C=======================================================================
      SUBROUTINE rlRESDm2(X,Y,THETA,N,NP,MDX,RS)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(NP),RS(N)
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     COMPUTES RESIDUALS: RS(I)=Y(I)-SUM X(I,J)*THETA(J)
C-----------------------------------------------------------------------
      DO 200 I=1,N
         SUM=ZERO
         DO 100 J=1,NP
            SUM=SUM+X(I,J)*THETA(J)
 100     CONTINUE
         RS(I)=Y(I)-SUM
 200  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlSOLVm2(X,THETA,NP,K,MDX,MDT)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),THETA(MDT)
      DATA DZERO/0.D0/
C-----------------------------------------------------------------------
C     LET X BE THE K BY K UPPER TRIANGULAR MATRIX WITH ELEMENTS
C     X(1,1)...X(1,K),X(2,2)...X(2,K),...X(K,K).
C     SOLV SOLVES THE TRIANGULAR SYSTEM U*THETA=Y (BY BACK SUBSTI-
C     TUTION). ON INPUT Y IS CONTAINED IN THETA.  ON OUTPUT
C     THETA(1)...THETA(K) CONTAIN THE DESIRED SOLUTION.
C-----------------------------------------------------------------------
C     ERRORS
C     1   AN ELEMENT OF THE PRINCIPAL DIAGONAL OF X IS =0.
C-----------------------------------------------------------------------
      KP1=K+1
      DO 90 L=1,K
         SM=DZERO
         I=KP1-L
         IF (I.EQ.K) GOTO 60
         IP1=I+1
         DO 50 J=IP1,K
            SM=SM+X(I,J)*THETA(J)
 50      CONTINUE
 60      SM1=SM
C         IF (X(I,I)) 80,70,80
         if(X(I,I) > 0.D0 .or. X(I,I) < 0.D0) then
           goto 80
         else
           goto 70
         end if
 70      CALL XERROR('Singular matrix',15,10,-1)
 80      THETA(I)=(THETA(I)-SM1)/X(I,I)
 90   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlPERMm2(X,SP,N,NDIM)
C.......................................................................
      DOUBLE PRECISION TMP,X(NDIM)
      INTEGER SP(NDIM)
C-----------------------------------------------------------------------
C     PERMUTE COMPONENTS OF X TO COMPENSATE COLUMN INTERCH.
C-----------------------------------------------------------------------
      DO 10 JJ=1,N
         J=N-JJ+1
         IF (SP(J).EQ.J) GOTO 10
         L=SP(J)
         TMP=X(L)
         X(L)=X(J)
         X(J)=TMP
 10   CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlQRSSm2(RS,WGT,WGT2,N,ITYPE,SIGMA,CONST,QR,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RS(N),WGT(N),WGT2(N)
      DATA ZERO,ONE/0.D0,1.D0/
C-----------------------------------------------------------------------
      TMP=ZERO
      IF (ITYPE.NE.1) GOTO 15
C-----------------------------------------------------------------------
C     HUBER-TYPE
C-----------------------------------------------------------------------
      DO 10 I=1,N
        S=RS(I)/SIGMA
        TMP=TMP+rlRHOm2(S,IPS,XK)
 10   CONTINUE
      GOTO 50
C-----------------------------------------------------------------------
C     MALLOWS-TYPE
C-----------------------------------------------------------------------
 15   IF (ITYPE.NE.2) GOTO 30
      DO 20 I=1,N
         IF (WGT(I) .EQ. ZERO .OR. WGT(I) .EQ. -ONE) GOTO 20
         S=RS(I)/SIGMA
         TMP=TMP+rlRHOm2(S,IPS,XK)*WGT(I)
 20   CONTINUE
      GOTO 50
C-----------------------------------------------------------------------
C     SCHWEPPE-TYPE
C-----------------------------------------------------------------------
 30   DO 40 I=1,N
         IF (WGT(I) .EQ. ONE .OR. WGT(I) .EQ. -ONE) GOTO 40
         S=RS(I)/(SIGMA*WGT(I))
         TMP=TMP+rlRHOm2(S,IPS,XK)*WGT2(I)
 40   CONTINUE
 50   QR=(TMP+CONST)*SIGMA
      RETURN
      END
C=======================================================================
      SUBROUTINE rlQRSHm2(RS,N,NP,SIGMA,QR,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RS(N)
C-----------------------------------------------------------------------
      TMP=0.D0
      DO 10 I=1,N
        S=RS(I)/SIGMA
        TMP=TMP+rlRHOm2(S,IPS,XK)
 10   CONTINUE
      QR=TMP/DBLE(N-NP)
      RETURN
      END
C=======================================================================
      FUNCTION rlICTHm2(NP,NCOV,DELTA,SIGMA,S,TOL,ICNV)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER rlICTHm2
      DIMENSION DELTA(NP),S(NCOV)
C-----------------------------------------------------------------------
      rlICTHm2=0
      TOL1=TOL*SIGMA
      IF (ICNV.EQ.2) GOTO 200
      IF (ICNV.EQ.3) GOTO 300
      L=0
      DO 100 J=1,NP
         L=L+J
         TOL2=TOL1*DSQRT(S(L))
         IF (TOL2 .LT. DABS(DELTA(J))) RETURN
 100  CONTINUE
      GOTO 500
 200  CALL rlXSYm2(DELTA,DELTA,S,NP,NCOV,TOL2)
      TOL2=DSQRT(TOL2)
      IF (TOL1 .GE. TOL2) rlICTHm2=1
      RETURN
 300  L=0
      DO 350 J=1,NP
         L=L+J
         TOL2=DABS(DELTA(J))*DSQRT(S(L))
         IF (TOL1 .LT. TOL2) RETURN
 350  CONTINUE
 500  rlICTHm2=1
      RETURN
      END
C=======================================================================
      SUBROUTINE rlFACSm2(RS,N,K,SIGMA,TL,XKAPPA,SUM2,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RS(N)
C-----------------------------------------------------------------------
C     COMPUTES CORRECTION FACTORS XKAPPA AND SUM2 FOR
C     THE COVARIANCE MATRIX.
C     rlFACSm2 CALLS THE FUNCTIONS rlPSIm2 AND rlPSPm2
C-----------------------------------------------------------------------
      TMP1=0.D0
      TMP2=0.D0
      DN=DBLE(N)
      DO 10 J=1,N
         S=RS(J)/SIGMA
         TMP1=TMP1+rlPSPm2(S,IPS,XK)
         PS=rlPSIm2(S,IPS,XK)
         TMP2=TMP2+PS*PS
 10   CONTINUE
      XMU=TMP1/DN
      SUM2=TMP2
      VAR=0.D0
      DO 20 J=1,N
         S=RS(J)/SIGMA
         VAR=VAR+(rlPSPm2(S,IPS,XK)-XMU)**2
 20   CONTINUE
      VAR=VAR/DN
      XKAPPA=0.D0
      IF (XMU.LE.TL) RETURN
      XMU2=XMU*XMU
      XKAPPA=1.D0+DBLE(K)*VAR/DN/XMU2
      SUM2=SUM2/XMU2/DBLE(N-K)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlVSVm2(LPIVOT,L1,M,U,IUE,UP,S,NCOV,SB)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(NCOV),U(IUE,M),SB(M)
      INTEGER H,HP1
C-----------------------------------------------------------------------
C     LET V BE THE ELEMENTARY HOUSEHOLDER TRANSFORMATION DEFINED
C     BY THE VECTOR U=(U(1)...U(N)) (WITH U(1)=...=U(LPIVOT-1)=0,
C     U(LPIVOT)=UP,U(LPIVOT+1)=...=U(K)=0 AND U(K+1)...U(N) POSSIBLY
C     DIFFERENT FROM 0) AND S A SYMMETRIC MATRIX STORED COLUMNWISE
C     IN THE ARRAY S OF LENGTH NCOV=N*(N+1)/2.  rlVSVm2 COMPUTES THE
C     SYMMETRIC MATRIX V*S*V AND STORES IT IN S.  IUE IS THE STORAGE
C     INCREMENT BETWEEN ELEMENTS OF THE VECTOR U.  THE LPIVOT-TH COMPO-
C     NENT OF U IS STORED IN UP WHEREAS THE STORAGE LOCATION U(LPIVOT)
C     CONTAINS THE NORM S (S.L6,P.55,FORMULA 10.5).
C-----------------------------------------------------------------------
C     ERRORS
C     1.  NCOV.NE.N*(N+1)/2.  NO COMPUTATION DONE BY rlVSVm2.
C-----------------------------------------------------------------------
      IF (L1.GT.M) RETURN
      ONE=1.D0
      K=L1-1
      B=U(1,LPIVOT)*UP
C      IF (B) 15,999,999
      if(B < 0.D0) then
        goto 15
      else
        goto 999
      end if
 15   B=ONE/B
C-----------------------------------------------------------------------
C     COMPUTE THE SCALAR PRODUCTS OF U WITH THE H-TH COLUMN OF S FOR H=1..M
C-----------------------------------------------------------------------
      L=0
      DO 85 H=1,M
         L=L+H
         L0=L-H
         IF (H.GE.LPIVOT) GOTO 20
         I=(LPIVOT-1)*LPIVOT/2+H
         GOTO 30
 20      I=L0+LPIVOT
 30      SM=UP*S(I)
         IF (H.LE.K) GOTO 60
         L0=L0+K
         DO 40 I=L1,H
            L0=L0+1
            SM=SM+S(L0)*U(1,I)
 40      CONTINUE
         HP1=H+1
         IF (H.EQ.M) GOTO 80
         DO 50 J=HP1,M
            L0=L0+J-1
            SM=SM+U(1,J)*S(L0)
 50      CONTINUE
         GOTO 80
 60      L0=(K-1)*K/2+H
         DO 70 J=L1,M
            L0=L0+J-1
            SM=SM+U(1,J)*S(L0)
 70      CONTINUE
 80      SB(H)=SM*B
 85   CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE THE QUADRATIC FORM U**T*S*U
C-----------------------------------------------------------------------
      SM=UP*SB(LPIVOT)
      DO 90 J=L1,M
         SM=SM+SB(J)*U(1,J)
 90   CONTINUE
      S1=SM*B
C-----------------------------------------------------------------------
C     SET U(LPIVOT)=UP
C-----------------------------------------------------------------------
      CSC=U(1,LPIVOT)
      U(1,LPIVOT)=UP
C-----------------------------------------------------------------------
C     COMPUTE S(1,LPIVOT)...S(LPIVOT-1,LPIVOT)
C-----------------------------------------------------------------------
      LPM1=LPIVOT-1
      L0=LPIVOT*LPM1/2
      IF (LPM1.LT.1) GOTO 105
      DO 100 I=1,LPM1
         L0=L0+1
         S(L0)=S(L0)+SB(I)*U(1,LPIVOT)
 100  CONTINUE
 105  CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE S(LPIVOT,LPIVOT)
C-----------------------------------------------------------------------
      L0=L0+1
      S(L0)=S(L0)+U(1,LPIVOT)*(S1*U(1,LPIVOT)+2.D0*SB(LPIVOT))
C-----------------------------------------------------------------------
C     COMPUTE S(LPIVOT,LPIVOT+1)...S(LPIVOT,L1-1)
C-----------------------------------------------------------------------
      LPP1=LPIVOT+1
      IF (LPP1.GT.K) GOTO 115
      DO 110 J=LPP1,K
         L0=L0+J-1
         S(L0)=S(L0)+SB(J)*U(1,LPIVOT)
 110  CONTINUE
 115  CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE S(I,J),J=L1...N,I=1...K (J.GE.L1,I.LT.L1)
C-----------------------------------------------------------------------
      L0=(K+1)*K/2-K
      DO 135 J=L1,M
         L0=L0+J-1
         DO 130 I=1,K
            L=L0+I
            S(L)=S(L)+SB(I)*U(1,J)
 130     CONTINUE
         I=LPIVOT
         L=L0+I
         S(L)=S(L)+SB(J)*U(1,I)+U(1,I)*S1*U(1,J)
 135  CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE S(I,J),I=L1...M,J=I...M
C-----------------------------------------------------------------------
      L0=K*(K+1)/2-K
      DO 145 J=L1,M
         L0=L0+J-1
         DO 140 I=L1,J
            L=L0+I
            S(L)=S(L)+(S1*U(1,I)*U(1,J)+(U(1,I)*SB(J)+U(1,J)*SB(I)))
 140     CONTINUE
 145  CONTINUE
      U(1,LPIVOT)=CSC
 999  RETURN
      END
C=======================================================================
      SUBROUTINE rlDOTPm2(X,Y,N,INCX,INCY,NX,NY,RESULT)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NX),Y(NY)
C-----------------------------------------------------------------------
      DTEMP=0.D0
      RESULT=0.D0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1) GOTO 20
C-----------------------------------------------------------------------
C     CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1
C-----------------------------------------------------------------------
      IX=1
      IY=1
      IF (INCX.LT.0) IX=(-N+1)*INCX+1
      IF (INCY.LT.0) IY=(-N+1)*INCY+1
      DO 10 I=1,N
         DTEMP=DTEMP+X(IX)*Y(IY)
         IX=IX+INCX
         IY=IY+INCY
 10   CONTINUE
      RESULT=DTEMP
      RETURN
C-----------------------------------------------------------------------
C     CODE FOR BOTH INCREMENTS EQUAL TO 1
C-----------------------------------------------------------------------
 20   M=MOD(N,5)
      IF (M.EQ.0) GOTO 40
      DO 30 I=1,M
         DTEMP=DTEMP+X(I)*Y(I)
 30   CONTINUE
      IF (N.LT.5) GOTO 60
 40   MP1=M+1
      DO 50 I=MP1,N,5
         DTEMP=DTEMP+X(I)*Y(I)+X(I+1)*Y(I+1)+X(I+2)*Y(I+2)
     +        +X(I+3)*Y(I+3)+X(I+4)*Y(I+4)
 50   CONTINUE
 60   RESULT=DTEMP
      RETURN
      END
C=======================================================================
      SUBROUTINE rlRMTRm2(X,N,NP,MDX,INTCH,TAU,K,SF,SG,SH,IP)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),SF(NP),SG(NP),SH(NP)
      INTEGER IP(NP)
C-----------------------------------------------------------------------
      FACTOR=0.001D0
      LDIAG=MIN0(N,NP)
      DO 80 JJ=1,LDIAG
         J=JJ
         IF (INTCH.EQ.0) IP(J)=J
         IF (INTCH.EQ.0) GOTO 70
         IF (J.EQ.1) GOTO 20
C-----------------------------------------------------------------------
C     UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
C-----------------------------------------------------------------------
         LMAX=J
         DO 10 L=J,NP
            SH(L)=SH(L)-X(J-1,L)**2
            IF(SH(L).GT.SH(LMAX)) LMAX=L
 10      CONTINUE
C         IF (HMAX+FACTOR*SH(LMAX)-HMAX) 20,20,50
      if((HMAX+FACTOR*SH(LMAX)-HMAX) > 0.D0) then
        goto 50
      else
        goto 20
      end if
C-----------------------------------------------------------------------
C     COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
C-----------------------------------------------------------------------
 20      LMAX=J
         DO 40 L=J,NP
            SH(L)=0.D0
            DO 30 I=J,N
               SH(L)=SH(L)+X(I,L)**2
 30         CONTINUE
            IF (SH(L).GT.SH(LMAX)) LMAX=L
 40      CONTINUE
         HMAX=SH(LMAX)
C-----------------------------------------------------------------------
C     LMAX HAS BEEN DETERMINED: INTERCHANGE COLUMNS IF NEEDED
C-----------------------------------------------------------------------
 50      CONTINUE
         IP(J)=LMAX
         IF (IP(J).EQ.J) GOTO 70
         DO 60 I=1,N
            TMP=X(I,J)
            X(I,J)=X(I,LMAX)
            X(I,LMAX)=TMP
 60      CONTINUE
         SH(LMAX)=SH(J)
C-----------------------------------------------------------------------
C     COMPUTE THE HOUSEHOLDER TRANSF. Q AND APPLY IT TO X
C-----------------------------------------------------------------------
 70      MDC=NP-J
         IF (MDC.GT.0)
     +        CALL rlH12m2(1,J,J+1,N,X(1,J),1,SH(J),X(1,J+1),1,MDX,
     +        MDC,MDX*MDC)
         IF (MDC.EQ.0)
     +        CALL rlH12m2(1,J,J+1,N,X(1,J),1,SH(J),SF,1,1,0,1)
 80   CONTINUE
C-----------------------------------------------------------------------
C     X CONTAINS NOW THE TRANSFORMED DESIGN MATRIX Q*X.
C     DETERMINE THE PSEUDORANK K USING THE TOLERANCE TAU
C-----------------------------------------------------------------------
      DO 100 J=1,LDIAG
         IF (DABS(X(J,J)).LE.TAU) GOTO 110
 100  CONTINUE
      K=LDIAG
      GOTO 120
 110  K=J-1
 120  KP1=K+1
C-----------------------------------------------------------------------
C     IF THE PSEUDORANK IS LESS THAN NP STORE THE FIRST K
C     DIAG.ELEMENTS OF X FOR FURTHER APPLICATIONS OF Q
C-----------------------------------------------------------------------
      IF (K.EQ.NP) GOTO 130
      DO 125 I=1,K
         SF(I)=X(I,I)
 125  CONTINUE
 130  CONTINUE
C-----------------------------------------------------------------------
C     SPECIAL FOR PSEUDORANK=0
C-----------------------------------------------------------------------
      IF (K.GT.0) GOTO 140
C     CALL XERROR('Pseudorank is zero',18,10,2)
      RETURN
C-----------------------------------------------------------------------
C     IF THE PSEUDORANK IS LESS THAN NP COMPUTE Q*X*V
C-----------------------------------------------------------------------
 140  IF (K.EQ.NP) GOTO 160
      MDC=MDX*(NP-1)
      DO 150 II=1,K
         I=KP1-II
         CALL rlH12m2(1,I,KP1,NP,X(I,1),MDX,SG(I),X,MDX,1,I-1,
     +        MDC+I-1)
 150  CONTINUE
 160  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlRWAGm2(X,Y,THETA,WGT,COV,PSP0,SIGMAI,N,NP,
     +     MDX,NCOV,TOL,GAM,TAU,ITYPE,ISIGMA,ICNV,MAXIT,MAXIS,NIT,
     +     SIGMAF,RS,DELTA,SC,SF,SG,SH,IP,SW,SX,IPS,XK,BETA,BET0)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(NP),WGT(N),COV(NCOV),RS(N),
     +     DELTA(NP),SC(N),SF(NP),SG(NP),SH(NP),SW(N),SX(MDX,NP)
      INTEGER IP(NP),rlICTHm2,rlISIGm2
      DATA ZERO,ONE,TL/0.D0,1.D0,1.D-10/
C-----------------------------------------------------------------------
C     W-ALGORITHM FOR ROBUST AND BOUNDED INFLUENCE LINEAR REGRESSION
C-----------------------------------------------------------------------
      MDXP1=MDX+1
      LDIAG=MIN0(N,NP)
      SIGMA=SIGMAI
      SIGMB=SIGMA
      IASG=IABS(ISIGMA)
      INTCH=1
      ITYP=ITYPE
      IF (ITYP.EQ.1) GOTO 15
      N0=N
      E=2.D0
      IF (ITYP.EQ.2) E=0.5D0
      DO 10 I=1,N
        IF (WGT(I).LE.ZERO) THEN
          SW(I)=-ONE
          N0=N0-1
        ELSE
          SW(I)=WGT(I)**E
        ENDIF
 10   CONTINUE
      IF (N0.EQ.0) ITYP=1
 15   IF (IASG.EQ.0) CONST=ZERO
      IF (IASG.EQ.1) CONST=BETA*DBLE(N-NP)
      IF (IASG.EQ.2) CONST=BET0*DBLE(N-NP)
C-----------------------------------------------------------------------
C     STEP 1. SET NIT := 1
C-----------------------------------------------------------------------
      NIT=1
C-----------------------------------------------------------------------
C     STEP 2. COMPUTE rlRESDm2 AS R=Y-X*THETA
C-----------------------------------------------------------------------
 200  CALL rlRESDm2(X,Y,THETA,N,NP,MDX,RS)
C-----------------------------------------------------------------------
C     STEP 3. COMPUTE A NEW VALUE SIGMB FOR SIGMA.
C-----------------------------------------------------------------------
      IF (ISIGMA .LT. 0 .AND. NIT .EQ. 1) GOTO 300
      IF (ISIGMA .EQ. 0) GOTO 300
      SIGMA=SIGMB
      CALL rlRSIGm2(RS,WGT,SIGMA,N,NP,TOL,ITYP,ISIGMA,MAXIS,
     +     NIS,SIGMB,SW,SC,IPS,XK,BETA,BET0)
      IF (SIGMB.LE.TL) RETURN
C-----------------------------------------------------------------------
C     STEP 4. COMPUTE WEIGHTS AND APPLY THEM TO X; STORE RESULT IN SX
C-----------------------------------------------------------------------
 300  DO 430 I=1,N
         SC(I)=PSP0
         IF (RS(I) .EQ. ZERO) GOTO 410
         T=RS(I)/SIGMB
         IF (ITYP.EQ.1) GOTO 400
         SC(I)=ZERO
         IF (WGT(I) .LE. ZERO) GOTO 410
         IF (ITYP.EQ.2) GOTO 400
         T=T/WGT(I)
 400     SC(I)=rlPSIm2(T,IPS,XK)/T
 410     PI=DSQRT(SC(I))
         IF (ITYP.EQ.2) PI=PI*SW(I)
         RS(I)=PI*RS(I)
         DO 420 J=1,NP
            SX(I,J)=PI*X(I,J)
 420     CONTINUE
 430  CONTINUE
C-----------------------------------------------------------------------
C     STEP 5. SOLVE FOR DELTA
C-----------------------------------------------------------------------
      CALL rlRMTRm2(SX,N,NP,MDX,INTCH,TAU,K,SF,SG,SH,IP)
      IF (K.EQ.0) RETURN
      KK=MDX*(K-1)+K
      IF (K.NE.NP) CALL rlSWAPm2(SX,SF,K,MDXP1,1,KK,K)
      DO 500 JJ=1,LDIAG
         J=JJ
         CALL rlH12m2(2,J,J+1,N,SX(1,J),1,SH(J),RS,1,N,1,N)
 500  CONTINUE
      IF (K.NE.NP) CALL rlSWAPm2(SX,SF,K,MDXP1,1,KK,K)
C-----------------------------------------------------------------------
C     SOLVE FOR DELTA
C-----------------------------------------------------------------------
      CALL rlSOLVm2(SX,RS,NP,K,MDX,N)
      IF (K.EQ.NP) GOTO 530
      KP1=K+1
      DO 510 J=KP1,NP
         RS(J)=ZERO
 510  CONTINUE
      DO 520 J=1,K
         I=J
         CALL rlH12m2(2,I,KP1,NP,SX(I,1),MDX,SG(I),RS,1,N,1,NP)
 520  CONTINUE
 530  DO 540 J=1,NP
         DELTA(J)=GAM*RS(J)
 540  CONTINUE
      CALL rlPERMm2(DELTA,IP,LDIAG,NP)
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
      IF (ISIGMA.LT.0.AND.NIT.EQ.1) GOTO 700
      IF(rlICTHm2(NP,NCOV,DELTA,SIGMA,COV,TOL,ICNV).EQ.1
     +     .AND. rlISIGm2(SIGMA,SIGMB,TOL).EQ.1) GOTO 800
 700  NIT=NIT+1
      GOTO 200
 800  SIGMAF=SIGMB
      CALL rlRESDm2(X,Y,THETA,N,NP,MDX,RS)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlRSIGm2(RS,WGT,SIGMAI,N,NP,TOL,ITYPE,ISIGMA,MAXIS,
     +     NIT,SIGMAF,SW,SC,IPS,XK,BETA,BET0)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RS(N),WGT(N),SW(N),SC(N)
      INTEGER rlISIGm2
      DATA ONE,ZERO,TL/1.D0,0.D0,1.D-10/
C-----------------------------------------------------------------------
C     UPDATE THE SCALE PARAMETER OF AN M-ESTIMATE
C-----------------------------------------------------------------------
      N0=N
      SIGMB=SIGMAI
      IASG=IABS(ISIGMA)
      ITYP=ITYPE
      IF (ITYP.EQ.1) GOTO 20
      IF (SIGMAI.EQ.SIGMAF) GOTO 20
      E=2.D0
      IF (ITYP.EQ.2) E=0.5D0
      DO 10 I=1,N
         IF (WGT(I) .LE. ZERO) THEN
            SW(I)=-ONE
            N0=N0-1
         ELSE
            SW(I)=WGT(I)**E
         ENDIF
 10   CONTINUE
      IF (N0.EQ.0) ITYP=1
 20   CONTINUE
      IF (IASG.EQ.2) GOTO 500
      CONST=BETA*DBLE(N-NP)
C-----------------------------------------------------------------------
C     STEP 1. SET NIT := 1
C-----------------------------------------------------------------------
      NIT=1
C-----------------------------------------------------------------------
C     STEP 2. COMPUTE A NEW VALUE SIGMB FOR SIGMA
C-----------------------------------------------------------------------
 100  SIGMA=SIGMB
      CALL rlNSIGm2(RS,WGT,SW,SIGMA,SIGMB,N,ITYP,IPS,XK,CONST)
      IF (SIGMB.GT.TL) GOTO 300
      RETURN
C-----------------------------------------------------------------------
C     STEP 3. STOP ITERATIONS IF DESIRED PRECISION IS REACHED
C-----------------------------------------------------------------------
 300  IF (rlISIGm2(SIGMA,SIGMB,TOL).EQ.1.OR.NIT.EQ.MAXIS) GOTO 400
      NIT=NIT+1
      GOTO 100
 400  SIGMAF=SIGMB
      RETURN
C-----------------------------------------------------------------------
C     COMPUTE SIGMA USING MEDIAN
C-----------------------------------------------------------------------
 500  IF (ITYPE.NE.1) GOTO 650
C-----------------------------------------------------------------------
C     HUBER-TYPE
C-----------------------------------------------------------------------
      DO 600 I=1,N
         SC(I)=DABS(RS(I))
 600  CONTINUE
      N0=N
      GOTO 900
C-----------------------------------------------------------------------
C     MALLOWS
C-----------------------------------------------------------------------
 650  IF (ITYPE.NE.2) GOTO 750
      N0=0
      DO 700 I=1,N
         IF (SW(I) .LE. 0.D0) GOTO 700
         N0=N0+1
         SC(N0)=DABS(RS(I))*SW(I)
 700  CONTINUE
      GOTO 900
C-----------------------------------------------------------------------
C     SCHWEPPE-TYPE
C-----------------------------------------------------------------------
 750  N0=0
      DO 800 I=1,N
         IF (WGT(I) .EQ. ZERO) GOTO 800
         N0=N0+1
         SC(N0)=DABS(RS(I))
 800  CONTINUE
 900  MED=(N0/2)+1
      CALL rlSTORm2(SC,N0,MED,SIGMAF)
      SIGMAF=SIGMAF/BET0
      RETURN
      END
C=======================================================================
      SUBROUTINE rlHSESm2(X,Y,N,NP,NQ,MDX,MDW,MDI,IOPT,INTCH,
     +     NREP,TOLS,TOLR,TAU,MAXS1,ISEED,IERR,SMIN,
     +     THETA,RS,IT1,WORK,IWORK,IPS,XK,BETA,BET0,ITRACE)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(N),RS(N),IT1(NQ)
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
      CALL rlHSE2m2(X,Y,N,NP,NQ,MDX,IOPT,INTCH,NREP,TOLS,TOLR,TAU,
     +     MAXS1,ISEED,IERR,SMIN,THETA,RS,IT1,WORK(1),WORK(N0),
     +     WORK(N1),WORK(N2),WORK(N3),WORK(N4),WORK(N5),IWORK(1),
     +     IWORK(NP1),IPS,XK,BETA,BET0,ITRACE)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlHSE2m2(X,Y,N,NP,NQ,MDX,IOPT,INTCH,NREP,TOLS,TOLR,
     +     TAU,MAXS1,ISEED,IERR,SMIN,THETA,RS,IT1,XX,YY,XTHETA,
     +     SF,SG,SH,SZ,SP,IT,IPS,XK,BETA,BET0,ITRACE)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(N),RS(N),XX(NQ,NP),YY(NQ)
      DIMENSION XTHETA(NQ),SF(NP),SG(NP),SH(NP),SZ(N)
      INTEGER IT1(NQ),SP(NP),IT(NQ)
      LOGICAL ALLZERO
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C     Resampling algorithm for the computation of S-estimates
C-----------------------------------------------------------------------
C     STEP 0: INITIALIZATIONS
C-----------------------------------------------------------------------
      K1=N/2+1
      CONST=BETA*DBLE(N-NP)
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
      PSP0=rlPSPm2(ZERO,IPS,XK)
C-----------------------------------------------------------------------
C     STEP 1: DRAW A SUBSAMPLE
C-----------------------------------------------------------------------
 100  IF (IOPT.NE.3) THEN
         DO 130 K=1,NQ
 110        CALL rlRNDm2(ISEED,RND)
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
            CALL rlNCOMm2(N,NQ,IT)
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
      CALL rlRMTRm2(XX,NQ,NP,NQ,INTCH,TAU,KK,SF,SG,SH,SP)
      IF(KK.NE.NP) GOTO 700
C-----------------------------------------------------------------------
C     STEP 3: SOLVE SYSTEM OF LINEAR EQUATIONS
C-----------------------------------------------------------------------
      CALL rlRICLm2(XX,YY,NQ,NP,NQ,XTHETA,SH,SP)
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
         IF (S .EQ. 1.0D7) GOTO 800
         CALL rlSTORm2(SZ,N,K1,S0)
         S0=2.D0*S0
         IF (S0 .EQ. ZERO) S0=S
         SRES=S0
      ENDIF
 435  D=ZERO
      DO 440 I=1,N
         D=D+rlCHIm2(RS(I)/SRES,IPS,XK)
 440  CONTINUE
      IF (SMIN .NE. ZERO .AND. D .GT. CONST) GOTO 700
      IF (D .LE. CONST) GOTO 500
      S0=1.5D0*S0
      SRES=S0
      GOTO 435
C-----------------------------------------------------------------------
C     STEP 5: SOLVE FOR SRES
C-----------------------------------------------------------------------
 500  CALL rlRSIGm2(RS,SZ,S0,N,NP,TOLR,ITYPE,1,MAXS1,NIS,SRES,SZ,
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
      RETURN
      END
C=======================================================================
      SUBROUTINE rlKTASm2(X,N,NP,MDX,NCOV,TAU,F,A,COV)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),A(NCOV),COV(NCOV)
      DATA DZERO/0.D0/
C-----------------------------------------------------------------------
C     COVARIANCE MATRIX OF THE FORM f*(X'X)^(-1)
C-----------------------------------------------------------------------
      NN=NP*(NP+1)/2
C-----------------------------------------------------------------------
C     COMPUTE X'X AND STORE IT TEMPORARILY IN COV
C-----------------------------------------------------------------------
      L=0
      DO 65 I=1,NP
         DO 60 J=1,I
            L=L+1
            SM1=DZERO
            DO 50 K=1,N
               SM1=SM1+X(K,I)*X(K,J)
 50         CONTINUE
            COV(L)=SM1
 60      CONTINUE
 65   CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE A LOWER TRIANGULAR MATRIX A SUCH THAT
C     (X'X)^(-1)=A'A; SET COV=A'A.
C-----------------------------------------------------------------------
      CALL rlMCHLm2(COV,NP,NN,INFO)
      IF (INFO.EQ.0) GOTO 68
      RETURN
 68   DO 70 L=1,NN
         A(L)=COV(L)
 70   CONTINUE
      CALL rlMINVm2(A,NP,NN,TAU,ISING)
      IF (ISING.EQ.0) GOTO 75
      RETURN
 75   CALL rlMTT1m2(A,COV,NP,NN)
      IF (F .GT. 0.D0) CALL rlSCALm2(COV,F,NCOV,1,NCOV)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlKFFAm2(RS,N,NP,SIGMA,FH,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RS(N)
      DATA TL/1.D-10/
C-----------------------------------------------------------------------
C     HUBER'S CORRECTION FACTOR FOR AS. COVARIANCE
C     MATRIX OF PARAMETER ESTIMATES
C-----------------------------------------------------------------------
      FH=1.D0
      IF (NP.EQ.N) RETURN
      CALL rlFACSm2(RS,N,NP,SIGMA,TL,XKAPPA,SUM2,IPS,XK)
      IF (XKAPPA .EQ. 0.D0) RETURN
      FH=(XKAPPA*XKAPPA)*SUM2
      RETURN
      END
C=======================================================================
      subroutine rlGENEm2(x,y,n,np,npopsize,probmutate,initgen,
     +     nbirths,nstock,maxslen,objvec,ntable,nstocklen,noldstock,
     +     stockprob,intch,tolr,tau,maxs1,smin,theta,rs,
     +     sz,sp,sg,sf,xtheta,yy,sh,xx,ntind,ips,xk,beta,bet0)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      dimension x(n,np),y(n),objvec(npopsize),theta(np),rs(n)
      dimension xx(maxslen,np),yy(maxslen),xtheta(maxslen),probmutate(4)
      dimension sf(np),sg(np),sh(np),sz(n),stockprob(npopsize)
      integer npins(2),nstock(maxslen,npopsize)
      integer nstocklen(npopsize),ntind(maxslen)
      integer ntable(2*maxslen),sp(np)
C-----------------------------------------------------------------------
      psp0 = rlPSPm2(0.D0,ips,xk)
      itype = 1
      call fseedi()
      if (noldstock .GT. 0) then
        goto 9
      else
        goto 11
      endif
 9    DO 10 I = 1, NOLDSTOCK
         call rlGEN2m2(x,y,nstock(1,i),nstocklen(i),n,np,
     +        maxslen,xx,yy,xtheta,rs,sres,sf,sg,
     +        sh,sp,sz,intch,tolr,tau,maxs1,itype,ips,xk,
     +        beta,bet0)
         objvec(i) = sres
 10   CONTINUE
 11   if (noldstock .LT. npopsize) then
         goto 12
      else
         goto 31
      endif
 12   do 30 i=noldstock+1,npopsize
         call getrandind(n,np,maxslen,ntind,ni)
         call rlGEN2m2(x,y,ntind,ni,n,np,maxslen,xx,yy,
     +        xtheta,rs,sres,sf,sg,sh,sp,sz,
     +        intch,tolr,tau,maxs1,itype,ips,xk,
     +        beta,bet0)
         objvec(i) = sres
         do 20 j=1,ni
            nstock(j,i) = ntind(j)
 20      CONTINUE
         nstocklen(i) = ni
 30   CONTINUE
 31   call rlGMAXm2(stmax,indmax,npopsize,objvec)
C-----------------------------------------------------------------------
C     Take random samples
C-----------------------------------------------------------------------
      if (initgen .GT. 0) then
         i = 1
 40      if (i .LE. initgen) then
            call getrandind(n,np,maxslen,ntind,ni)
            call rlGEN2m2(x,y,ntind,ni,n,np,maxslen,xx,yy,
     +           xtheta,rs,sres,sf,sg,sh,sp,sz,
     +           intch,tolr,tau,maxs1,itype,ips,xk,
     +           beta,bet0)
            tmp = sres
            if (tmp .LT. stmax) then
               objvec(indmax) = tmp
               do 50 j=1,ni
                  nstock(j,indmax) = ntind(j)
 50            CONTINUE
               nstocklen(indmax) = ni
               call rlGMAXm2(stmax,indmax,npopsize,objvec)
            endif
            i = i + 1
            goto 40
         endif
      endif
C-----------------------------------------------------------------------
C     GENETIC PART
C-----------------------------------------------------------------------
      i = 1
 60   if (i .LE. nbirths) then
         call marriage(nstock,maxslen,npopsize,stockprob,
     +        nstocklen,probmutate,ntind,ni,n,np,npins,ntable)
         call rlGEN2m2(x,y,ntind,ni,n,np,maxslen,xx,yy,xtheta,rs,
     +        sres,sf,sg,sh,sp,sz,intch,tolr,tau,maxs1,itype,ips,xk,
     +        beta,bet0)
         tmp = sres
         ntestm = 0
         ntestf = 0
         if (tmp .GT. objvec(npins(1))) ntestm = 1
         if (tmp .GT. objvec(npins(2))) ntestf = 1
         if (ntestm .EQ. 1 .OR. ntestf .EQ. 1) then
            if (objvec(npins(1)) .LT. objvec(npins(2))) then
               itmp = npins(2)
            else
               itmp = npins(1)
            endif
            objvec(itmp) = tmp
            nstocklen(itmp) = ni
            do 70 j=1,ni
               nstock(j,itmp) = ntind(j)
 70         CONTINUE
         endif
         i = i + 1
         goto 60
      endif
      call fseedo()

      tmp = 1.0d36
      do 80 j=1,npopsize
         if (objvec(j) .LT. tmp) then
            tmp = objvec(j)
            ind = j
         endif
 80   continue
      ni = nstocklen(ind)
      do 90 i=1,ni
         ntind(i) = nstock(i,ind)
 90   CONTINUE
      call rlGEN2m2(x,y,ntind,ni,n,np,maxslen,xx,yy,xtheta,rs,sres,
     +     sf,sg,sh,sp,sz,intch,tolr,tau,maxs1,itype,ips,xk,
     +     beta,bet0)
      smin = sres
      do 95 k=1,np
         theta(k) = xtheta(k)
 95   CONTINUE
C-----------------------------------------------------------------------
C     Compute the values to be returned
C-----------------------------------------------------------------------
      do 110 i=1,n
         s = y(i)
         do 100 j=1,np
            s = s - theta(j)*x(i,j)
 100     CONTINUE
         rs(i) = s
 110  CONTINUE
      RETURN
      END
C=======================================================================
      subroutine rlGEN2m2(x,y,ntind,ni,n,np,maxslen,xx,yy,xtheta,rs,
     +     sres,sf,sg,sh,sp,sz,intch,tolr,tau,maxs1,itype,ips,xk,
     +     beta,bet0)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      integer ntind(ni),sp(np)
      dimension x(n,np),y(n),xx(maxslen,np),yy(maxslen)
      dimension sf(np),sg(np),sh(np),xtheta(maxslen)
      dimension rs(n),sz(n)
C-----------------------------------------------------------------------
      const = beta*dble(n-np)
      big = 1.0d36
      k1 = n/2 + 1
      do 100 i=1,ni
         itmp = ntind(i)
         do 90 j=1,np
            xx(i,j) = x(itmp,j)
 90      CONTINUE
         yy(i) = y(itmp)
 100  CONTINUE
C-----------------------------------------------------------------------
C     Decompose Sampled Matrix: tmat
C-----------------------------------------------------------------------
      call rlRMTRm2(xx,ni,np,maxslen,intch,tau,kk,sf,sg,sh,sp)
      if (kk .NE. np) then
         sres = big
         return
      endif
      call rlRICLm2(xx,yy,ni,np,maxslen,xtheta,sh,sp)
      do 120 i=1,n
         s = y(i)
         do 110 j=1,np
            s = s - xtheta(j)*x(i,j)
 110     CONTINUE
         rs(i) = s
 120  CONTINUE
      s = dble(1.0E7)
      do 130 i=1,n
         ari = dabs(rs(i))
         sz(i) = ari
         if (ari .NE. 0.D0) s = dmin1(s,ari)
 130  CONTINUE
      CALL rlSTORm2(SZ,N,K1,S0)
      S0 = 2.D0*S0
      IF (S0 .EQ. 0.D0) S0 = S
      SRES = S0
 135  D = 0.D0
      DO 140 I=1,N
         D = D + rlCHIm2(RS(I)/SRES,IPS,XK)
 140  CONTINUE
      IF (D .LE. CONST) GOTO 150
      S0 = 1.5D0*S0
      SRES = S0
      GOTO 135
 150  CALL rlRSIGm2(rs,sz,s0,n,np,tolr,itype,1,maxs1,nis,
     +     sres,sz,sz,ips,xk,beta,bet0)
      RETURN
      END
C=======================================================================
      subroutine rlGMAXm2(stmax,indmax,npopsize,objvec)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      dimension objvec(npopsize)
C-----------------------------------------------------------------------
      stmax = objvec(1)
      indmax = 1
      if (npopsize .GT. 1) then
         j = 2
 5       if (j .LE. npopsize) then
            if (objvec(j) .GT. stmax) then
               stmax = objvec(j)
               INDMAX = J
            ENDIF
            J = J +1
            GOTO 5
         ENDIF
      ENDIF
      RETURN
      END
C=======================================================================
      FUNCTION rlPSIm2(S,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Compute  psi(s, xk) -- vectorized as rlPSIAm2(N,SVALS,FVALS,IPS,XK)
c
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C     IPS = 4: SMOOTHED HUBER FUNCTION
C-----------------------------------------------------------------------
      ABST=DABS(S)
      IF (IPS .EQ. 1) GOTO 100
      IF (IPS .EQ. 2) GOTO 200
      IF (IPS .EQ. 3) GOTO 300
      IF (IPS .EQ. 4) GOTO 400
 100  R1= -1.944D0
      R2=  1.728D0
      R3= -0.312D0
      R4=  0.016D0
      AX=ABST/XK
      IF (AX .GT. 3.D0) THEN
         rlPSIm2=0.D0
      ELSE IF(AX .GT. 2.D0) THEN
         AX=S/XK
         IF (AX .GT. 0.D0) THEN
            rlPSIm2=DMAX1(0.D0,XK*(R4*AX**7+R3*AX**5+R2*AX**3+R1*AX))
         ELSE
            rlPSIm2=-DABS(XK*(R4*AX**7+R3*AX**5+R2*AX**3+R1*AX))
         ENDIF
      ELSE
         rlPSIm2=S
      ENDIF
      RETURN
 200  rlPSIm2=0.D0
      IF (ABST .LT. XK) THEN
         SK=S/XK
         rlPSIm2=(6.D0*SK/XK)*(1.D0-SK*SK)*(1.D0-SK*SK)
      ENDIF
      RETURN
 300  rlPSIm2=DMIN1(XK,ABST)
      IF (S .LT. 0.D0) rlPSIm2=-rlPSIm2
      RETURN
 400  IF (ABST .LE. XK) THEN
         rlPSIm2=S
      ELSE
         rlPSIm2=S/ABST*XK*(1.D0+(1.D0-(ABST/XK)**(-3.D0))/3.D0)
      ENDIF
      RETURN
      END
C=======================================================================
      FUNCTION rlRHOm2(S,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Compute  rho(s, xk) -- vectorized as rlRHOAm2(N,SVALS,FVALS,IPS,XK)
c
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C-----------------------------------------------------------------------
      ABST=DABS(S)
      S2=S*S
      IF (IPS .EQ. 1) GOTO 100
      IF (IPS .EQ. 2) GOTO 200
      IF (IPS .EQ. 3 .OR. IPS .EQ. 4) GOTO 300
 100  R1=-1.944D0/2.0D0
      R2= 1.728D0/4.0D0
      R3=-0.312D0/6.0D0
      R4= 0.016D0/8.0D0
      AX=ABST/XK
      IF (AX .GT. 3.D0) THEN
         rlRHOm2=3.25D0*XK*XK
      ELSE IF (AX .GT. 2.D0) THEN
         rlRHOm2=XK*XK*(R1*AX**2+R2*AX**4+R3*AX**6+R4*AX**8+1.792D0)
      ELSE
         rlRHOm2=S2/2.D0
      ENDIF
      RETURN
 200  rlRHOm2=1.D0
      IF (ABST .LT. XK) THEN
         S2=S2/(XK**2)
         rlRHOm2=(S2*(S2-3.D0)+3.D0)*S2
      ENDIF
      RETURN
 300  rlRHOm2=S2/2.D0
      IF (ABST .GT. XK) rlRHOm2=XK*(ABST-XK/2.D0)
      RETURN
      END
C=======================================================================
      FUNCTION rlPSPm2(S,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     Compute  psi'(s, xk)  = psp(s,.) -- vectorized as rlPSPAm2(N,SVALS,FVALS,IPS,XK)
c
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C     IPS = 4: SMOOTHED HUBER FUNCTION
C-----------------------------------------------------------------------
      ABST=DABS(S)
      IF (IPS .EQ. 1) GOTO 100
      IF (IPS .EQ. 2) GOTO 200
      IF (IPS .EQ. 3) GOTO 300
      IF (IPS .EQ. 4) GOTO 400
 100  R1= -1.944D0
      R2=  1.728D0
      R3= -0.312D0
      R4=  0.016D0
      AX=ABST/XK
      IF (AX .GT. 3.D0) THEN
         rlPSPm2=0.D0
      ELSE IF (AX .GT. 2.D0) THEN
         rlPSPm2=7.D0*R4*AX**6+5.D0*R3*AX**4+3.D0*R2*AX**2+R1
      ELSE
         rlPSPm2=1.D0
      ENDIF
      RETURN
 200  rlPSPm2=0.D0
      IF (ABST .LT. XK) THEN
         S2=(S/XK)**2
         rlPSPm2=(6.D0/XK)*(1.D0-S2)*(1.D0-5.D0*S2)/XK
      ENDIF
      RETURN
 300  rlPSPm2=0.D0
      IF (ABST .LE. XK) rlPSPm2=1.D0
      RETURN
 400  rlPSPm2=1.D0
      IF (ABST .GT. XK) rlPSPm2=(ABST/XK)**(-3.D0)
      RETURN
      END
C=======================================================================
      FUNCTION rlCHIm2(S,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
c Compute  chi(s, xk)  -- vectorized as rlCHIAm2(N,SVALS,FVALS,IPS,XK)
c
C     COMPUTES THE VALUE OF CHI FUNCTION
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C-----------------------------------------------------------------------
      IF (IPS .EQ. 1) GOTO 100
      IF (IPS .EQ. 2) GOTO 200
      IF (IPS .EQ. 3 .OR. IPS .EQ. 4) GOTO 300
c                    ^^^^^^^^^^^^^^^ why?  FIXME?
 100  R1=-1.944D0/2.0D0
      R2= 1.728D0/4.0D0
      R3=-0.312D0/6.0D0
      R4= 0.016D0/8.0D0
      AX=DABS(S/XK)
      IF (AX .GT. 3.D0) THEN
         rlCHIm2=3.25D0*XK*XK
      ELSE IF (AX .GT. 2.D0) THEN
         rlCHIm2=XK*XK*(R1*AX**2+R2*AX**4+R3*AX**6+R4*AX**8+1.792D0)
      ELSE
         rlCHIm2=S*S/2.D0
      ENDIF
      RETURN
 200  rlCHIm2=1.D0
      ABST=DABS(S)
      IF (ABST .LT. XK) THEN
         S2=(S/XK)**2
         rlCHIm2=(S2*(S2-3.D0)+3.D0)*S2
      ENDIF
      RETURN
 300  ABST=DABS(S)
      PS=DMIN1(XK,ABST)
      rlCHIm2=PS*PS/2.D0
      RETURN
      END
c_______________________________________________________________________
c
c--> The following *vectorized* versions of
c-->     psi(), psi'() = psp(), rho(), chi()
c--> are  called from S / R  as   psi.weight(), etc :
c
C=======================================================================
      SUBROUTINE rlPSIAm2(N,SVALS,FVALS,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SVALS(N),FVALS(N)
C-----------------------------------------------------------------------
c computes  psi(s[i], .) , i = 1:n
c
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C     IPS = 4: SMOOTHED HUBER FUNCTION
C-----------------------------------------------------------------------
      DO 150 I=1,N
         FVALS(I)=rlPSIm2(SVALS(I),IPS,XK)
 150  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlRHOAm2(N,SVALS,FVALS,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SVALS(N),FVALS(N)
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF RHO FUNCTION
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C-----------------------------------------------------------------------
      DO 150 I=1,N
         FVALS(I)=rlRHOm2(SVALS(I),IPS,XK)
 150  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlPSPAm2(N,SVALS,FVALS,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SVALS(N),FVALS(N)
C-----------------------------------------------------------------------
c computes  psi'(s[i], *) = psp(s[i], *), i = 1:n
c
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C     IPS = 4: SMOOTHED HUBER FUNCTION
C-----------------------------------------------------------------------
      DO 150 I=1,N
         FVALS(I)=rlPSPm2(SVALS(I),IPS,XK)
 150  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlCHIAm2(N,SVALS,FVALS,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SVALS(N),FVALS(N)
C-----------------------------------------------------------------------
c computes  chi(s[i], *) , i = 1:n

C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C-----------------------------------------------------------------------
      DO 150 I=1,N
         FVALS(I)=rlCHIm2(SVALS(I),IPS,XK)
 150  CONTINUE
      RETURN
      END
