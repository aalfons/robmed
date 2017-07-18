C wblrob.f
C 26/04/00
C
C 11/08/00 # fixed RLWLA123


      double precision FUNCTION RLWEIBLN(TAU,V,X)
      implicit double precision(a-h,o-z)
      double precision lower

      CALL RLWEILIM(TAU,V,LOWER,UPPER)
      RLWEIBLN=0.D0
      IF (X.LE.LOWER) RETURN
      IF (X.GE.UPPER) RETURN
      Z=(X-TAU)/V
      TMP=Z-DEXP(Z)
      TMP=DEXP(TMP)/V
      RLWEIBLN=TMP
      RETURN
      END
C
C-----------------------------------------------------------------------
C

      DOUBLE PRECISION FUNCTION RLXLOGD(X)
C.......................................................................
      implicit double precision(a-h,o-z)
      DATA NCALL,XMIN,YMIN/0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
        CALL RLMACHD(4,XMIN)
        CALL RLMACHD(5,YMIN)
        NCALL=1
      ENDIF 
C
C  EXTENDED NATURAL LOGARITHM FUNCTION
C
      IF (X.LE.0.D0) THEN
        RLXLOGD=0.D0
      ELSEIF (X.LE.XMIN) THEN
        RLXLOGD=YMIN
      ELSE
        RLXLOGD=DLOG(X)
      ENDIF
      RETURN
      END 

C
C==========================================================================
C
      DOUBLE PRECISION FUNCTION RLWZANS(DX,WGT,N,EXU,EXWLN,TAU,V,
     *            A11,A21,A22,B1,B2,C1,C2,U12X11,BETA,yb)
      implicit double precision(a-h,o-z)
      DIMENSION WGT(N),yb(8,2)
      EXTERNAL EXWLN,EXU,RLXEXPD 
C
C avoid compiler warnings
C
      dummy = EXU(1.D0,1,1.D0,1.D0)
C
C  Initializations
C
      dummy = u12x11
      dummy = yb(1,1)
      ANS=EXWLN(TAU,V,DX)
      RLWZANS=0.D0
      IF (ANS.EQ.0.D0) RETURN 
      BB1=DBLE(B1)
      BB2=DBLE(B2)
      G=RLXEXPD(DX)
      X1=G-1.D0-C1
      Z1=A11*X1
      U1=1.D0
      TMP=DABS(Z1)
      IF (TMP.GT.BB1) U1=BB1/TMP
      IF (WGT(1).EQ.4.D0) GOTO 40
      X2=DX*G-DX-1.D0-C2
      Z2=A22*X2+A21*X1
      U2=1.D0
      TMP=DABS(Z2)
      IF (TMP.GT.BB2) U2=BB2/TMP
      I=int(WGT(1))
c      GOTO (10,20,30,40,50,60,70,80) I
      if (I == 2) then
        goto 20
      else if (I == 3) then
        goto 30
      else if (I == 4) then
        goto 40
      else if (I == 5) then
        goto 50
      else if (I == 6) then
        goto 60
      else if (I == 7) then
        goto 70
      else if (I == 8) then
        goto 80
      else
        goto 10
      end if
   10 RLWZANS=U1*U2*X1*X2*DBLE(ANS)
      RETURN
   20 RLWZANS=U1*U2*X1*X1*DBLE(ANS)
      RETURN
   30 RLWZANS=(U2*(BETA*X1+X2))**2*DBLE(ANS)
      RETURN
   40 RLWZANS=(U1*X1)**2*DBLE(ANS)
      RETURN
   50 RLWZANS=U2*X2*DBLE(ANS)
      RETURN
   60 RLWZANS=U2*X1*DBLE(ANS)
      RETURN
   70 RLWZANS=U1*Z1*U2*Z2*DBLE(ANS)
      RETURN
   80 RLWZANS=U2*Z2*DBLE(ANS)
      RETURN
      END
C
C========================================================================
C
      SUBROUTINE RLINTUXW(X,NREP,IOPT,TIL,SUM,WLO,WHI,
     *     TAU,V,A11,A21,A22,B1,B2,C1,C2,U12X11,BETA,YB)

      implicit double precision(a-h,o-z)
      dimension ss(7),X(NREP),YB(8,2),WGT(1),work(320),iwork(80)
      double precision lo
      EXTERNAL RLEXU,RLWZANS,RLWEIBLN
C
C     INITIALISATION
C
      LIMIT=80
      KEY=1
      HI=WLO
      WGT(1)=dble(IOPT)
      X(NREP)=WHI
      SUM=0.D0
      TILD=DBLE(TIL)
      I=0
  100 I=I+1
      LO=HI
      IF (I.GT.NREP) GOTO 200
      HI=DMIN1(X(I),WHI)
      IF (LO.GE.HI) THEN
        SS(I)=0.D0
        GOTO 100
      ENDIF
      CALL RLINTGRW(RLWZANS,WGT,1,RLEXU,RLWEIBLN,LO,HI,TILD,TILD,
     1            KEY,LIMIT,SS(I),ERRSTD,NEVAL,IER,WORK,IWORK,TAU,V,
     2              A11,A21,A22,B1,B2,C1,c2,U12X11,
     3              BETA,YB)
      SUM=SUM+SS(I)
      IF (DABS(HI-WHI).LT.1.D-6) GOTO 200
      GOTO 100
  200 RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLWEQTA1(AA,FA,A11,a21,a22,B1,b2,C1,c2,TOLD,wlo,whi,
     *    tau,v,u12x11,beta,yb)

      implicit double precision(a-h,o-z)
      dimension x(3),yb(8,2)
      EXTERNAL RLXLOGD 

      XL=-B1/A11+1.D0+C1
      NSOL=0
      IF (XL.GT.0.D0) THEN
        NSOL=NSOL+1
        X(1)=RLXLOGD(XL)
      ENDIF
      XU=B1/A11+1.D0+C1
      IF (XU.GT.0.D0) THEN
        NSOL=NSOL+1
        X(NSOL)=RLXLOGD(XU)
      ENDIF
      AA=0.D0
      NSOL=NSOL+1
      CALL RLINTUXW(X,NSOL,4,dble(TOLD),SUM1,wlo,whi,
     *    tau,v,a11,a21,a22,b1,b2,c1,c2,u12x11,beta,yb)

      XL=SUM1
      IF (DABS(SUM1).LT.1.D-10) XL=DSIGN(1.D-10,SUM1)
      AA=1.D0/DSQRT(XL)
      FA=A11*A11*SUM1-1.D0
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLWEQTA2(AA,Fa,A11,a21,A22,B1,b2,C1,c2,U12X11,beta,
     *        yb,wlo,whi,tau,v,X2,NSOL,told)

      implicit double precision (a-h, o-z)
      dimension SC(7),Z0(7),x(7),x2(4),yb(8,2)
      EXTERNAL RLXLOGD 
      DATA Z0/7*0.D0/

      DO 10 I=1,NSOL
      X(I)=X2(I)
   10 CONTINUE
      XL=-B1/A11+1.D0+C1
      N1=0
      IF (XL.GT.0.D0) THEN
        N1=1
        X(NSOL+1)=RLXLOGD(XL)
      ENDIF
      XU=B1/A11+1.D0+C1
      IF (XU.GT.0.D0) THEN
        N1=N1+1
        X(NSOL+N1)=RLXLOGD(XU)
      ENDIF
      N2=NSOL+N1
c      DO 20 I=1,N2
c   20 SC(I)=X(I) 
      DO I=1,N2
        SC(I)=X(I) 
      END DO
      CALL RLSRT2(SC,Z0,7,1,N2)
c      DO 30 I=1,N2
c   30 X(I)=SC(I)
      DO I=1,N2
        X(I)=SC(I)
      END DO
      IOPT=1
   50 continue
      CALL RLINTUXW(X,N2+1,IOPT,TOLD,SUMI,wlo,whi,tau,v,
     *        a11,a21,a22,b1,b2,c1,c2,u12x11,beta,yb)
      IF (IOPT.EQ.1) THEN
        U12X12=SUMI
        IOPT=2
        GOTO 50
      ENDIF
      U12X11=SUMI
      IF (U12X11.LT.1.D-6) U12X11=1.D-6
      BETA=-U12X12/U12X11
      IF (NSOL.GT.0) CALL RLSRT2(X2,Z0,NSOL,1,NSOL)   
c      DO 150 I=1,NSOL
c  150 X(I)=X2(I)
      DO I=1,NSOL
        X(I)=X2(I)
      END DO
      CALL RLINTUXW(X,NSOL+1,3,TOLD,SUM,wlo,whi,tau,v,
     *        a11,a21,a22,b1,b2,c1,c2,u12x11,beta,yb)
      XL=SUM
      IF (SUM.LT.1.D-10) XL=DSIGN(1.D-10,SUM)
      AA=1.D0/DSQRT(XL)
      FA=A22*A22*SUM-1.D0
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLWEQTA3(AA,FA,A11,A21,A22,U12X11,BETA)
      implicit double precision(a-h,o-z)
      AA=BETA*AA
      FA=A11*U12X11*(A21-A22*BETA)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLWLA123(MAXIT,TOL,IOPT,A,FA,NIT,
     *                A11,A21,A22,B1,B2,
     *                C1,C2,WLO,WHI,tau,v,NSOL,X2,Y2,
     *                U12X11,BETA,yb)
      implicit double precision(a-h,o-z) 
      dimension A(3),FA(3),X2(4),Y2(0:4), yb(8,2)
C
C  STEP 0 : INITIALIZATION
C  ------
      NIT=1
      TLD=TOL*TOL
      A11=A(1)
      A21=A(2)
      A22=A(3)
      IF (IOPT.EQ.1) FA(1)=0.D0
      FA(2)=0.D0
      FA(3)=0.D0
C
C  STEP 1: Compute new value of (Aij) and check convergence
C  ------

  100 IF (IOPT.EQ.1) CALL RLWEQTA1(a1,FA(1),A11,a21,a22,B1,
     *                b2,c1,c2,ToL,wlo,whi,tau,v,
     *                u12x11,beta,yb)
      CALL RLSOLWX(B2,TOL,NSOL,X2,Y2,A21,A22,C1,C2,WLO,WHI)
      CALL RLWEQTA2(a2,FA(2),A11,a21,A22,B1,b2,C1,c2,
     *        U12X11,beta,yb,wlo,whi,tau,v,X2,NSOL,tld)
      a3=a2
      CALL RLWEQTA3(a3,FA(3),A11,A21,A22,U12X11,BETA)
      S=(fa(3))**2+(fa(2))**2+(fa(1))**2
      IF (IOPT.EQ.1) A11=a1
      A21=a3
      A22=a2
      IF (S.LT.TLD) GOTO 200
      if(nit.ge.maxit)goto 200
      nit=nit+1
      goto 100
C
C  STEP 2: EXIT
C  ------
  200 A(1)=A11
      A(2)=A21
      A(3)=A22
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLCRETABW(SB1,SB2,A,MAXTA,MAXTC,MAXIT,TIL,TOL,
     +           MONIT,TAB,TPAR)

      implicit double precision(a-h,o-z)
      DOUBLE PRECISION ZERO(2),TAB(5),ver(9),tpar(6),a(3), 
     +     fa(3),calf(2),ac(5), 
     +     veca(24),tabc(3,8),x2(4), y2(0:4),yb(8,2)
C
      DATA VECA/-0.3857397D0, 3.2524269D0, 4.701294D0,-0.3805570D0,
     +  2.5858480D0,3.742037D0,-0.3745761D0,2.1899542D0,3.175076D0,
     + -0.3681503D0,1.9200350D0,2.806144D0,-0.3613750D0,1.718788D0,
     +  2.550337D0,-0.3542827D0,1.5640367D0,2.366607D0,-0.346828D0,
     +  1.441594D0, 2.231880D0,0.D0,0.D0,0.D0/
c
c avoid compiler warnings
c
      idummy = MONIT
c

      IJ=0
c      DO 5 J=1,8
c      DO 5 I=1,3
c      IJ=IJ+1
c      TABC(I,J)=VECA(IJ)
c    5 CONTINUE
      DO J=1,8
        DO I=1,3
          IJ=IJ+1
          TABC(I,J)=VECA(IJ)
        END DO
      END DO
      B1 = SB1
      B2 = SB2
      N=1
      FA(1)=0.D0
      FA(2)=0.D0
      FA(3)=0.D0
      FA1=0.D0
      FA2=0.D0
      FA3=0.D0
      FC1=0.D0
      FC2=0.D0
      TAU=0.D0
      TILD=TIL
      TOLD=TOL
      V=1.D0
      delta=0.D0
      tpar(1)=B1
      tpar(2)=B2
      tpar(3)=1.D0
      tpar(4)=1.D0
      tpar(5)=1.D0
      tpar(6)=delta
      MXTA=MAXTA
      MXTC=MAXTC
      suma=dabs(a(1))+dabs(a(2))+dabs(a(3))
      DO 100 IA=0,1
      ALFA=tpar(3)+dble(IA)*tpar(6)
      A11=1.D0
      A22=0.779697D0 
      A21=-0.32964D0
      c1=0.D0
      c2=0.D0
      if (suma.ne.0.d0) then
        A11=A(1)
        A21=A(2)
        A22=A(3)
      endif
      CALL RLWEILIM(0.D0,1.D0,WLO,WHI)
      CALF(1)=c1
      CALF(2)=c2
      A(1)=A11
      A(2)=A21
      A(3)=A22
      CALL RLSOLWX(B2,TOL,NSOL,X2,Y2,A21,A22,C1,C2,wlo,whi)
      NIT=1
      MAXCUR=MAXIT-3
   10 continue
      IF (NIT.GT.MAXCUR) THEN
        CALL RLWLNAC1(MAXIT,TOL,AC,FA1,FC1,AA,CALF,FA,ZERO,NIT1,
     *     A11,a21,a22,B1,b2,C1,c2,WLO,WHI,tau,v,ux12,beta,yb,til)
        FC1=ZERO(1)
        FA1=FA(1)
        CALL RLWLNAC2(MAXIT,TOL,AC,FA2,FA3,AA,CALF,FA,ZERO,NIT2,
     *    A11,A21,A22,b1,B2,C1,C2,UX12,BETA,x2,y2,
     *    wlo,whi,tau,v,yb,nsol,til)
        FC2=ZERO(2)
        FA2=FA(2)
        FA3=FA(3)
        A(1)=A11
        A(2)=A21
        A(3)=A22
        NIT=MAXIT
        GOTO 40
      ENDIF
      CALL RLWLNC12(MXTC,TOL,1,CALF,ZERO,NITC,A11,A21,
     * A22,B1,B2,C1,C2,tau,v,nsol,wlo,whi,x2,y2,ux12,beta,yb,tild)
      CALL RLWLA123(MXTA,TOL,1,A,FA,NITA,A11,A21,A22,B1,B2,
     *    C1,C2,WLO,WHI,tau,v,NSOL,X2,Y2,UX12,BETA,yb) 
      CALL RLSOLWX(B2,TOL,NSOL,X2,Y2,A21,A22,C1,C2,WLO,WHI)
      call RLWEQTC2(FC2,F2A,F2b,x2,y2,a11,A21,A22,b1,b2,c1,c2,
     *       ux12,beta,yb,TAU,V,NSOL,wlo,whi,tild)
      IF (NIT.NE.MAXCUR.AND.DABS(FC2).GT.TOLD) GOTO 20
      CALL RLWEQTC1(FC1,F1A,F1B,A11,B1,C1,TAU,V)
      IF (NIT.NE.MAXCUR.AND.DABS(FC1).GT.TOLD) GOTO 20
      CALL RLWEQTA1(AA,FA1,A11,a21,a22,B1,b2,C1,c2,TiLD,wlo,whi,
     *        tau,v,ux12,beta,yb)
      IF (NIT.NE.MAXCUR.AND.DABS(FA1).GT.TOLD) GOTO 20
      CALL RLWEQTA2(aa,FA2,A11,a21,A22,B1,b2,C1,c2,
     *        UX12,beta,yb,wlo,whi,tau,v,X2,NSOL,tild)
      IF (NIT.NE.MAXCUR.AND.DABS(FA2).GT.TOLD) GOTO 20
      CALL RLWEQTA3(AA,FA3,A11,A21,A22,UX12,BETA)
      IF (NIT.NE.MAXCUR.AND.DABS(FA3).GT.TOLD) GOTO 20
      IF (NIT.LT.MAXCUR) GOTO 40
   20 IF (NIT.EQ.MAXCUR) GOTO 30
      NIT=NIT+1
      GOTO 10
   30 CONTINUE
      IF (MAXCUR.EQ.MAXIT-3) THEN
        BMIN=DMIN1(B1,B2)
        IF (BMIN.LT.1.0875D0) THEN
          IB=1
        ELSEIF (BMIN.LT.1.1125D0) THEN
          IB=2
        ELSEIF (BMIN.LT.1.1375D0) THEN
          IB=3 
        ELSEIF (BMIN.LT.1.1625D0) THEN
          IB=4
        ELSEIF (BMIN.LT.1.1875D0) THEN
          IB=5
        ELSEIF (BMIN.LT.1.2125D0) THEN
          IB=6
        ELSEIF (BMIN.LT.1.2375D0) THEN
          IB=7
        ELSE
          IB=8
          TABC(1,8)=AC(5)
          TABC(2,8)=AC(3)
          TABC(3,8)=AC(2)
        ENDIF
        AC(1)=A11
        AC(2)=TABC(3,IB)
        AC(3)=TABC(2,IB)
        AC(4)=C1
        AC(5)=TABC(1,IB)
        NIT=NIT+1
        MAXCUR=MAXIT-2
        GOTO 10
      ELSEIF (MAXCUR.EQ.MAXIT-2) THEN
        NIT=MAXCUR+1
        GOTO 10
      ENDIF
   40 tab(1)=calf(1)
      tab(2)=calf(2)
      tab(3)=a(1)
      tab(4)=a(2)
      tab(5)=a(3)
        ver(1)=alfa
        ver(2)=dble(nit)
        ver(3)=fc1
        ver(4)=fc2
        ver(5)=fa1
        ver(6)=fa2
        ver(7)=fa3
        ver(1)=dble(nsol)
        ver(2)=x2(1)
        ver(3)=x2(2)
        ver(4)=x2(3)
        ver(5)=x2(4)
      if (delta.eq.0.D0) RETURN 
  100 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLWEQTC2(F,FP1,FP2,x,y,a11,A21,A22,b1,b2,c1,c2,
     *        u12x11,beta,yb,TAU,V,NSOL,wlo,whi,told)
      implicit double precision(a-h,o-z)
      dimension sum(10), xx(5), x(4), y(0:4), yb(8,2)

      XL=X(1)
      XU=X(2)
      F=Y(0)
      FP1=0.D0
      FP2=0.D0
      IF (NSOL.GE.2) THEN
c       DO 10 I=1,NSOL
c   10  XX(I)=X(I)
       DO I=1,NSOL
         XX(I)=X(I)
       END DO
       CALL RLINTUXW(XX,NSOL+1,8,TOLD,TMP,wlo,whi,tau,v,
     *        a11,a21,a22,b1,b2,c1,c2,u12x11,beta,yb)
       CALL RLSUMWLN(XL,TAU,V,SUM(1))
       CALL RLSUMWLN(XU,TAU,V,SUM(5))
       F=TMP
       FP1=-A21*(SUM(5)-SUM(1))
       FP2=-A22*(SUM(5)-SUM(1))
      ENDIF
      IF (NSOL.EQ.4) THEN
       XL=X(3)
       XU=X(4)
       CALL RLSUMWLN(XL,TAU,V,SUM(6))
       CALL RLSUMWLN(XU,TAU,V,SUM(10))
       FP1=FP1-A21*(SUM(10)-SUM(6))
       FP2=FP2-A22*(SUM(10)-SUM(6))
      ENDIF
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE RLWEQTC1(F,FP1,FP2,A11,B1,C1,TAU,V)
      implicit double precision(a-h,o-z)
      dimension sum(3)
      EXTERNAL RLXLOGD 

      XL=-B1/A11+1.D0+C1
      NSOL=0
      IF (XL.GT.0.D0) THEN
        NSOL=NSOL+1
        XL=RLXLOGD(XL)
      ENDIF
      XU=B1/A11+1.D0+C1
      IF (XU.GT.0.D0) THEN
        XU=RLXLOGD(XU)
        NSOL=NSOL+1
      ENDIF
      IF (NSOL.EQ.0) THEN
        TMP=B1
        FP1=0.D0
      ELSEIF (NSOL.EQ.1) THEN
        CALL RLEXPWLN(XU,TAU,V,SUM(2))
        CALL RLSUMWLN(XU,TAU,V,SUM(3))
        TMP=A11*SUM(2)-A11*(1.D0+C1)*SUM(3)+B1*(1.D0-SUM(3))
        FP1=-A11*SUM(3)
      ELSE
        CALL RLSUMWLN(XL,TAU,V,SUM(1))
        CALL RLEXPWLN(XL,TAU,V,TMP)
        CALL RLEXPWLN(XU,TAU,V,SUM(2))
        SUM(2)=SUM(2)-TMP
        CALL RLSUMWLN(XU,TAU,V,SUM(3))
        TMP=-B1*SUM(1)+A11*SUM(2)-A11*(1.D0+C1)*(SUM(3)-SUM(1))+
     *       B1*(1.D0-SUM(3))
        FP1=-A11*(SUM(3)-SUM(1))
      ENDIF
      F=TMP
      FP2=0.D0
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE RLWLNAC1(MAXIT,TOL,AC,FA1,FC1,AA,CALF,FA,ZERO,NIT,
     * A11,a21,a22,B1,b2,C1,c2,WLO,WHI,TAU,V,u12x11,beta,yb,tild)
      implicit double precision(a-h,o-z)
      dimension ZERO(2),calf(2),fa(3),ac(5),yb(8,2)
      DATA TL/1.D-6/
C
C  STEP 0.   INITIALIZATIONS
C  ------
C
      IF (DABS(FA1).LT.TOL.AND.DABS(FC1).LT.TOL) THEN
         FA(1)=FA1
         ZERO(1)=FC1
         CALF(1)=C1
         RETURN
      ENDIF
      IF (FA1.NE.1.D0.AND.DABS(FA1).GE.TOL) 
     + A11=A11+DSIGN(1.5D0,A11-AC(1))        
      C1=CALF(1)
      IF (FC1.NE.1.D0.AND.DABS(FC1).GE.TOL) 
     + C1=C1+2.D0*DSIGN(TOL,C1-AC(4))        
      NIT=0
      MAXTA=0
      MAXTC=1
      F10=FC1
      GOTO 130
C
C  STEP 1.   COMPUTE NEW VALUE of (A1,C1) AND CHECK CONVERGENCE
C  ------
C
  100 IF (MOD(NIT,20).EQ.0) THEN
        MAXTA=MAXTA+1
        MAXTC=MAXTC+1
      ENDIF
      DO 110 I=1,MAXTA
      CALL RLWEQTA1(AA,FA(1),A11,a21,a22,B1,b2,C1,c2,Tild,wlo,
     *        whi,tau,v,u12x11,beta,yb)
      A11=AA
  110 CONTINUE
      DO 120 I=1,MAXTC
      CALL RLWEQTC1(F10,F1A,F1B,A11,B1,C1,TAU,V)
      IF (DABS(F1A).LE.TL) F1A=DSIGN(TL,F1A)
      C1=C1-F10/F1A
  120 CONTINUE
  130 CALL RLWEQTA1(AA,FA(1),A11,a21,a22,b1,B2,C1,c2,Tild,
     *            wlo,whi,tau,v,u12x11,beta,yb)
      A11=AA
      IF (DABS(F10).LT.TOL.AND.DABS(FA(1)).LT.TOL) GOTO 200
      IF (NIT.EQ.MAXIT) GOTO 200
      NIT=NIT+1
      GOTO 100
C
C  STEP 2.   EXIT
C  ------
  200 CALF(1)=C1
      ZERO(1)=F10
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE RLWLNAC2(MAXIT,TOL,AC,FA2,FA3,AA,CALF,FA,ZERO,NIT,
     *        A11,A21,A22,b1,B2,C1,C2,UX12,BETA,x2,y2,WLO,WHI,tau,v,
     *        yb,NSOL,tild)
      implicit double precision(a-h, o-z)
      dimension ZERO(2),calf(2),fa(3),ac(5),x2(4),y2(0:4),yb(8,2)
C
C  STEP 0.   INITIALIZATIONS
C  ------
C
      NIT=1
      TIL=TOL
      C2=AC(5)
      A21=AC(3)
      A22=AC(2)
      CALL RLSOLWX(B2,TOL,NSOL,X2,Y2,A21,A22,C1,C2,WLO,WHI)
      CALL RLWEQTA2(aa,FA(2),A11,a21,A22,B1,b2,C1,c2,
     *        UX12,beta,yb,wlo,whi,tau,v,X2,NSOL,tild)
      FA(3)=A11*UX12*(A21-A22*BETA)
      call RLWEQTC2(F20,F2a,F2b,x2,y2,a11,A21,A22,b1,b2,c1,c2,
     *            ux12,beta,yb,TAU,V,NSOL,wlo,whi,tild)
      IF (DABS(F20).LT.TOL.AND.DABS(FA(2)).LT.TOL.AND.
     +  DABS(FA(3)).LT.TOL) GOTO 200
C
C  STEP 1.   COMPUTE NEW VALUE of (A22,A21,C2) AND CHECK CONVERGENCE
C  -----
  100 DO 120 II=1,1
      CALL RLSOLWX(B2,TOL,NSOL,X2,Y2,A21,A22,C1,C2,WLO,WHI)
      CALL RLWEQTA2(aa,FA(2),A11,a21,A22,B1,b2,C1,c2,
     *        UX12,beta,yb,wlo,whi,tau,v,X2,NSOL,tild)
        FA(3)=A11*UX12*(A21-A22*BETA)
        A22=AA
        A21=BETA*AA
  120 CONTINUE
      DO 130 II=1,1
        CALL RLSOLWX(B2,TOL,NSOL,X2,Y2,A21,A22,C1,C2,WLO,WHI)
        call RLWEQTC2(F20,F2a,F2b,x2,y2,a11,A21,A22,b1,b2,c1,c2,
     *            ux12,beta,yb,TAU,V,NSOL,wlo,whi,tild)
        IF (DABS(F2B).LE.1.D-6) F2B=DSIGN(1.D0,F2B)
        TMP=F20/F2B
        C2=C2-TMP
  130 CONTINUE
      IF (DABS(F20).LT.TOL.AND.DABS(FA(2)).LT.TOL.AND.
     +  DABS(FA(3)).LT.TOL) GOTO 200
      IF (NIT.EQ.MAXIT) GOTO 200
      NIT=NIT+1
      GOTO 100
C
C  STEP 2.   EXIT
C  ------
  200 CALF(2)=C2
      ZERO(2)=F20
      FA2=FA(2)
      FA3=FA(3)
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE RLWLNC12(MAXIT,TOL,IOPT,CALF,ZERO,NIT,A11,A21,A22,
     *             B1,B2,C1,C2,tau,v,nsol,wlo,whi,x,y,ux12,beta,yb,tild)
      implicit double precision(a-h,o-z)
      dimension ZERO(2), calf(2),x(4),y(0:4),yb(8,2)
      double precision LAMDA
      DATA TL/1.D-6/
C
C  STEP 0.   INITIALIZATIONS
C  ------
C
      NIT=1
      TOL2=TOL*TOL
      C1=CALF(1)
      C2=CALF(2)
      CALL RLWEQTC1(F10,F1A,F1B,A11,B1,C1,tau,v)
      call RLWEQTC2(F20,F2a,F2b,x,y,a11,A21,A22,b1,b2,c1,c2,
     *       ux12,beta,yb,TAU,V,NSOL,wlo,whi,tild)
      S=F10**2+F20**2
      IF (S.LT.TOL2) GOTO 200
C
C  STEP 1.   COMPUTE NEW VALUE of (C1,C2) AND CHECK CONVERGENCE
C  ------
C
  100 C01=C1
      C02=C2
      LAMDA=0.D0
      S0=S
  110 F1A=F1A+LAMDA
      F2B=F2B+LAMDA
      D=F1A*F2B-F1B*F2A
      IF (DABS(D).LT.TL) THEN
        LAMDA=LAMDA+DSQRT(TL)
        GOTO 110
      ENDIF
      DLT1=( F2B*F10-F1B*F20)/D
      DLT2=(-F2A*F10+F1A*F20)/D
      GAM1=1.D0
      GAM2=1.D0
      NSTP=0
  150 IF (IOPT.EQ.1) C1=C01-GAM1*DLT1
      C2=C02-GAM2*DLT2
      CALL RLSOLWX(B2,TOL,NSOL,X,Y,A21,A22,C1,C2,wlo,whi)
      CALL RLWEQTC1(F10,F1A,F1B,A11,B1,C1,tau,v)
      call RLWEQTC2(F20,F2a,F2b,x,y,a11,A21,A22,b1,b2,c1,c2,
     *            ux12,beta,yb,TAU,V,NSOL,wlo,whi,tild)
      S=F10**2+F20**2
      IF (S.LT.TOL2) GOTO 200
      IF (S.GT.S0.AND.NSTP.LT.10) THEN
        NSTP=NSTP+1
        GAM1=GAM1/2.D0
        GAM2=GAM2/2.D0
        GOTO 150
      ENDIF
      IF (NIT.EQ.MAXIT) GOTO 200
      NIT=NIT+1
      GOTO 100
C
C  STEP 2.   EXIT
C  ------
  200 CALF(1)=C1
      CALF(2)=C2
      ZERO(1)=F10
      ZERO(2)=F20
      RETURN
      END
C
C------------------------------------------------------------------------
C
C
      SUBROUTINE RLWEILIM(TAU,V,LOWER,UPPER)

      implicit double precision(a-h,o-z)
      double precision lower
      DATA NCALL,EXMIN,ZUP,ZLOW/0,0.D0,0.D0,0.D0/

      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
        ZUP=4.2D0
  100   ZUP=ZUP+0.01D0
        TMP=ZUP-DEXP(ZUP)
        IF (TMP.GT.EXMIN) GOTO 100
        ZUP=ZUP-0.05D0
        ZLOW=EXMIN+0.05D0
      ENDIF

      LOWER=ZLOW*V + TAU
      UPPER=ZUP*V + TAU
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLSUMWLN(HI,TAU,V,WL)
C
C  INT  Weibulln(x,tau,v)dx  for  x=-Infty to HI
C
      implicit double precision(a-h, o-z)
      double precision lower

      CALL RLWEILIM(TAU,V,LOWER,UPPER)
      WL=0.D0
      IF (HI.LT.LOWER) RETURN
      WL=1.D0
      IF (HI.GT.UPPER) RETURN
      ZHI=(HI-TAU)/V
      ZHI=-DEXP(ZHI)
      WL=1.D0-DEXP(ZHI)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLEXPWLN(HI,TAU,V,WL)
C
C  INT  exp((x-tau)/v)*Weibulln(x,tau,v) dx  for  x=-Infty to HI
C
      implicit double precision(a-h,o-z)
      double precision lower

      CALL RLWEILIM(TAU,V,LOWER,UPPER)
      WL=0.D0
      IF (HI.LT.LOWER) RETURN
      WL=1.D0
      IF (HI.GT.UPPER) RETURN
      ZHI=(HI-TAU)/V
      TMP=-DEXP(ZHI)
      WL=1.D0-DEXP(ZHI+TMP)-DEXP(TMP)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLZDERIV(X0,Y0,YPP0,A21,A22,C1,C2)
C
C  Solution of  y' = f'(z) = a22*(exp(z)+z*exp(z)-1)+a21*exp(z)=0
C        or     y' = a22*[1+z-exp(-z)] + a21 = 0
C  and value of y  = f(z0)
C      and      y''= f''(z0) = exp(z0)*[a22*(2+z0) + a21]               
C
      implicit double precision(a-h, o-z)
      EXTERNAL RLXEXPD 
C
C  NEWTON ALGORITHM FOR F'(Z)=0
C
      Z0=0.D0
      NIT=0           
      G =A22*(1.D0+Z0-RLXEXPD(-Z0)) + A21
  100 GP=A22*(1+RLXEXPD(-Z0))
      IF (DABS(GP).LT.1.D-6) GP=DSIGN(1.D-6,GP)
      Z=Z0-G/GP
      G=A22*(1.D0+Z-RLXEXPD(-Z)) + A21
      IF (DABS(G).GT.1.D-4) THEN
        Z0=Z
        NIT=NIT+1
        IF (NIT.LT.100) GOTO 100
      ENDIF  
      Z0=Z
      G=RLXEXPD(Z0)     
      FZ0=A21*(G-1.D0-C1)+A22*(Z0*G-Z0-1.D0-C2)
      FPPZ0=A22*(2.D0+Z0)+A21
      X0=Z0
      Y0=FZ0
      YPP0=FPPZ0
      RETURN
      END    
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLSOLWX(B,TOL,NSOL,XB,YB,A21,A22,C1,C2,wlo,whi)
C
C  Solutions of y=a22*(z*exp(z)-z-1-v*c2)+a21*(exp(z)-1-v*c1)=+/-v*B
C               a22 <> 0 and v = 1
C
      implicit double precision(a-h, o-z)
      DIMENSION XB(4),YB(0:4)
      do 50 j=1,4
      xb(j)=0.D0
  50  continue
      CALL RLZDERIV(X0,Y0,YPP0,A21,A22,C1,C2)
      BS=B
      IF (YPP0.GT.0.D0) BS=-B
      IF ((YPP0.GT.0.D0.AND.Y0.LT.-B).OR.
     *    (YPP0.LT.0.D0.AND.Y0.GT.B)) THEN
C       2 solutions for Y=-BS
        YB(0)=-BS
        CALL RLSOLWX0(-BS,TOL,x0,y0,1,XB(1),A21,A22,C1,C2,wlo,whi)
        YB(1)=0.D0
        CALL RLSOLWX0(-BS,TOL,x0,y0,2,XB(4),A21,A22,C1,C2,wlo,whi)
        YB(4)=-BS
C       2 solutions for Y=BS
        YB(2)=BS
        CALL RLSOLWX0(BS,TOL,x0,y0,1,XB(2),A21,A22,C1,C2,wlo,whi)
        YB(3)=0.D0
        CALL RLSOLWX0(BS,TOL,x0,y0,2,XB(3),A21,A22,C1,C2,wlo,whi)
        NSOL=4
      ELSEIF ((YPP0.LT.0.D0.AND.Y0.GT.-B).OR.
     +                       (YPP0.GT.0.D0.AND.Y0.LT.B)) THEN
C       0 solution for Y=BS & 2 solutions for Y=-BS
        YB(0)=-BS
        CALL RLSOLWX0(-BS,TOL,x0,y0,1,XB(1),A21,A22,C1,C2,wlo,whi)
        YB(1)=0.D0
        CALL RLSOLWX0(-BS,TOL,x0,y0,2,XB(2),A21,A22,C1,C2,wlo,whi)
        YB(2)=-BS
        NSOL=2
      ELSE
C       0 solution for Y=BS & 0 solution for y=-BS
        XB(1)=0.D0
        YB(0)=-BS
        NSOL=0
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RLFZY(Z,W,NOBS,param)
      implicit double precision(a-h, o-z)
      DIMENSION W(NOBS), param(2) 
      EXTERNAL RLXEXPD 
      DATA NCALL,XBIG/0,0.D0/

      A21 = param(1)
      A22 = param(2)

      IF (NCALL.EQ.0) THEN
        NCALL=1
        W(1)=1.D0
        NONE=1 
        CALL RLMACHD(6,XBIG)
      ENDIF
C
C  Initializations
C
      G=RLXEXPD(Z)     
      TMP=DABS(A21)+DABS(A22*Z)
      IF (TMP.LT.1.D0) GOTO 10 
      IF (G.GE.XBIG/TMP) G=XBIG/TMP  
  10  RLFZY=A21*G+A22*Z*(G-1.D0)
      RETURN      
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLSOLWX0(B,TOL,x0,y0,istep,X,A21,A22,C1,C2,wlo,whi)
      implicit double precision(a-h,o-z)    
      DIMENSION ZZ(1), param(2)
      EXTERNAL RLFZY
      DATA MAXIT/300/
C
C  Initializations
C
      Y=B+A22+A22*C2+A21+A21*C1
      TOLD=TOL
      Z=X0      
      IF(ISTEP.EQ.1)THEN
C
C  1 SOLUTION LESSER THAN X0
C
         param(1) = a21
         param(2) = a22
        AA=RLFZY(WLO,ZZ,1,param)
        X=WLO 
        IF ((AA.LE.B.AND.AA.GT.Y0).OR.(AA.GE.B.AND.AA.LT.Y0)) RETURN 
        AA=Z-2.D0
        BB=Z
        CALL RLRGFLD(RLFZY,ZZ,Y,AA,BB,TOLD,MAXIT,XX,ITERM,1,param)
        IF (ITERM.EQ.2) XX=WLO 
      ELSE
C
C  1 SOLUTION GREATER THAN X0
C
         param(1) = a21
         param(2) = a22
        AA=RLFZY(WHI,ZZ,1,param)
        X=WHI 
        IF ((AA.LE.B.AND.AA.GT.Y0).OR.(AA.GE.B.AND.AA.LT.Y0)) RETURN 
        AA=Z
        BB=Z+2.D0
        CALL RLRGFLD(RLFZY,ZZ,Y,AA,BB,TOLD,MAXIT,XX,ITERM,1,param)
        IF (ITERM.EQ.2) XX=WHI 
      ENDIF
      X=XX
      RETURN
      END
C
c
      SUBROUTINE RLWESTIM(Y,NOBS,TAB,MAXIT,TPAR,TOL,ALFA1,ALFA2,
     +           ALPHA,SIGMA,NIT,c1c2,a123)
C
C   INPUT: Y(NOBS)    Log(Observations)
C          TAB(5)     Table of A, c1, c2.
C          ALFA1      Search starts at ALFA1
C          ALFA2      Search ends   at ALFA2
C          TOL        Tolerance
C
C   OUTPUT: NIT = 0   No solution found for ALPHA in [ALFA1,ALFA2]
C           NIT = 1   Solution = (ALPHA, SIGMA, C1C2, A123)
C
      implicit double precision(a-h,o-z)
      DIMENSION Y(NOBS),TAB(5),c1c2(2),a123(3),tpar(6)
      dimension A(3),sign9(2),param(6)
      EXTERNAL RLWEQTN9,RLWEQTN10,RLXLOGD,RLXEXPD

      B1=TPAR(1)
      B2=TPAR(2)
      N=NOBS
      LA=2
      SIG0=SIGMA
      SIGM10=SIG0
      IF (SIG0.NE.0.D0) THEN
        TAU0=RLXLOGD(SIG0)
        LA=1
      ENDIF
      ITERMY=1
      TOLD=TOL
      DELTA=ALFA2-ALFA1
      IF (DELTA.GE.30.D0) THEN
        DLTA=0.2D0
      ELSEIF (DELTA.GE.10.D0) THEN
        DLTA=0.1D0
      ELSE
        DLTA=0.05D0
      ENDIF
      ALF1=ALFA1
      IF (ALF1.EQ.0.D0) ALF1=0.2D0
      ALF2=ALFA2
      IF (ALF2.EQ.0.D0) ALF2=25.2D0
      NIT=0
      NRFN=0
      SIGN9(1)=-9.D0
      SIGN9(2)=-9.D0
      I9=2
      ALF9=0.D0
      SIG9=0.D0
      TAU9=0.D0
      NK=-1
  100 NK=NK+1
      ALFA=ALF1+DBLE(NK)*DLTA
      IF (ALFA.GT.ALF2) GOTO 500
      IF (TAB(3).EQ.0.D0) THEN
        MAXTA=1
        MAXTC=1
        CALL RLCRETABW(B1,B2,A,MAXTA,MAXTC,MAXIT,TOL,TOL,0,TAB,TPAR)
      ENDIF
        C1=TAB(1)
        C2=TAB(2)
        A11=TAB(3)
        A21=TAB(4)
        A22=TAB(5)
      V=1.D0/ALFA
C
C  Special case: Solve only the 2nd equation for alpha
C
      IF (SIG0.NE.0.D0) THEN
        IF (NIT.EQ.1) THEN
          ALF9=ALFA
          NRFN=NRFN+1
          param(1) = v
          param(2) = a21
          param(3) = a22
          param(4) = b2
          param(5) = c1 
          param(6) = c2
          SIGN9(I9)=RLWEQTN9(TAU0,Y,NOBS,param)
          GOTO 300
        ENDIF
        I9=3-I9
        param(1) = v
        param(2) = a21
        param(3) = a22
        param(4) = b2
        param(5) = c1 
        param(6) = c2
        SIGN9(I9)=RLWEQTN9(TAU0,Y,NOBS,param)
        GOTO 120
      ENDIF
C
C  Solve equation 10 for sigma
C
      TMIN=Y(1)-V*RLXLOGD(1.D0+C1)
      TMAX=Y(NOBS)-V*RLXLOGD(1.D0+C1)
      param(1) = v
      param(2) = a11
      param(3) = b1
      param(4) = c1
      FMAX=RLWEQTN10(TMIN,Y,NOBS,param)
      FMIN=RLWEQTN10(TMAX,Y,NOBS,param)
      IF     (FMAX.LT.0.D0) THEN
             ITERMy=4
             SIGM10=RLXEXPD(TMIN)
      ELSEIF (FMIN.GT.0.D0) THEN
             ITERMy=3
             SIGM10=RLXEXPD(TMAX)
      ELSE
        CALL RLRGFLD(RLWEQTN10,Y,0.D0,TMIN,TMAX,TOLD,MAXIT,TAU0,
     *             ITERMY,NOBS,param)
        SIGM10=RLXEXPD(TAU0)
      ENDIF
      IF (NIT.EQ.1) GOTO 200
      IF (ITERMY.NE.1) GOTO 100
      I9=3-I9
      param(1) = v
      param(2) = a21
      param(3) = a22
      param(4) = b2
      param(5) = c1 
      param(6) = c2
      SIGN9(I9)=RLWEQTN9(TAU0,Y,NOBS,param)
  120 TMP=DABS(SIGN9(I9))
      IF (SIGN9(2).EQ.-9.D0) THEN
        FMIN9=DABS(SIGN9(1))
        SIGN9(2)=SIGN9(1)
        ALF9=ALFA
        SIG9=SIGM10
        GOTO 100
      ENDIF
      IF (TMP.LT.FMIN9) THEN
        ALF9=ALFA
        SIG9=SIGM10
        TAU9=TAU0
        FMIN9=TMP
      ENDIF
      IF (SIGN9(1)*SIGN9(2).GT.0.D0) GOTO 100
      NIT=1
C
C The solution is between ALFA-DLTA and ALFA. Refinement step
C
  150 NK=2*NK-2
      DLTA=DLTA/2.D0
      IF (10.D0*DLTA.LT.TOL) GOTO 500
      GOTO 100
  200 NRFN=NRFN+1
      IF (SIG0.EQ.0.D0) then
        param(1) = v
        param(2) = a21
        param(3) = a22
        param(4) = b2
        param(5) = c1 
        param(6) = c2
        SIGN9(I9)=RLWEQTN9(TAU0,Y,NOBS,param)
      endif
      TMP=DABS(SIGN9(I9))
      IF (TMP.LT.FMIN9) THEN
        ALF9=ALFA
        SIG9=SIGM10
        TAU9=TAU0
        FMIN9=TMP
      ENDIF
  300 IF (SIGN9(1)*SIGN9(2).GT.0.D0) THEN
        NK=2*NK
        DLTA=DLTA/2.D0
        IF (10.D0*DLTA.LT.TOL) GOTO 500
        I9=3-I9
        IF (NRFN.LE.10) GOTO 100
      ELSE
        IF (NRFN.LE.10) GOTO 150
      ENDIF
  500 ALPHA=ALF9
      SIGMA=SIG9
      TAU=TAU9
      V=1.D0/ALPHA
      c1c2(1)=c1/v
      c1c2(2)=c2/v
      a123(1)=a11*v
      a123(2)=a21*v
      a123(3)=a22*v
      RETURN
      END 
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLWEILIK(SY,N,MAXIT,TOL,ALPHA,SIGMA,ZERO,NIT)
C
C Solution for (@,beta) of n*beta - sum(y_i^@)=0   and
C n/@ + sum(log(y_i)) - n*sum(y_i^@*log(y_i))/sum(y_i^@)=0
C Then sigma=exp(log(beta)/@)
C
      implicit double precision (a-h,o-z)
      DIMENSION SY(N)
      double precision logx
      DATA NCALL,EXMIN,XLGMN,YLGMN,TIL/0,0.D0,0.D0,0.D0,1.D-6/
C
C Compute sample mean and variance; set initial values
C
      NIT=1
      TOLD=TOL
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
        CALL RLMACHD(4,XLGMN)
        CALL RLMACHD(5,YLGMN)
      ENDIF

      SUML=0.D0
      EN=DBLE(N)

      DO 10 I=1,N
      XI=SY(I)
      LOGX=YLGMN
      IF (XI.GT.XLGMN) LOGX=DLOG(XI)
      SUML=SUML+LOGX
      SY(I)=LOGX
   10 CONTINUE
      SUML=SUML/EN
      TEST=0.D0
      ALF0=1.D0
C
C Compute the objective function and its derivative
C
  100 ALFA=ALF0
      SUMYA=0.D0
      SUMYAL=0.D0
      SYAL2=0.D0
c      DO 110 I=1,N
c      TA=0.D0
c      TMP=ALFA*SY(I)
c      IF (TMP.GT.EXMIN) TA=DEXP(TMP)
c      SUMYA=SUMYA+TA
c      SUMYAL=SUMYAL+TA*SY(I)
c  110 SYAL2=SYAL2+TA*(SY(I)**2)
      DO I=1,N
        TA=0.D0
        TMP=ALFA*SY(I)
        IF (TMP.GT.EXMIN) TA=DEXP(TMP)
        SUMYA=SUMYA+TA
        SUMYAL=SUMYAL+TA*SY(I)
        SYAL2=SYAL2+TA*(SY(I)**2)
      END DO
      DEN=ALFA
      IF (DEN.LT.TIL) DEN=TIL
      F =1.D0/DEN + SUML - SUMYAL/SUMYA
      IF (TEST.GT.0.) GOTO 250
  120 ALF2=(ALFA)**2
      IF (ALF2.LT.TIL) ALF2=TIL
      FP=-1.D0/ALF2 - (SYAL2*SUMYA-SUMYAL**2)/(SUMYA**2)
C
C Compute new values
C
      NSTP=1
      IF (DABS(FP).LT.TIL) FP=DSIGN(TIL,FP)
      ALF0=ALFA-(F/FP)
  200 IF (ALF0.LE.0.d0) THEN
        NSTP=NSTP+1
        FP=2.D0*FP
        ALF0=ALFA-(F/FP)
        GOTO 200
      ENDIF
      TEST=1.D0
      ALF1=ALFA
      GOTO 100
  250 TEST=0.d0
      ZERO=(EN*F)
      IF (dABS(ZERO).LT.TOL) GOTO 300
      IF (dABS(ALF1-ALF0) .LT. DMIN1(1.D0,DABS(ALF0))*TOL 
     1     .AND. NSTP.LE.2) GOTO 300
      IF (NIT.EQ.MAXIT) GOTO 300
      NIT=NIT+1
      ALFA=ALF0
      GOTO 120
C
C Exit
C
  300 ALPHA=ALFA
      SIGMA=DEXP(DLOG(SUMYA/EN)/ALFA)
      RETURN
      END


      DOUBLE PRECISION FUNCTION RLWEQTN9(TAU,Y,NOBS,param)
      implicit double precision (a-h,o-z)
      DIMENSION Y(NOBS),param(6)
      EXTERNAL RLPSI2,RLXEXPD

      v = param(1)
      a21 = param(2)
      a22 = param(3)
      b2 = param(4)
      c1 = param(5)
      c2 = param(6)
      SUM=0.D0
      EN=dble(NOBS)
      DO 120 I=1,NOBS
        S=(Y(I)-TAU)/V
        TMP=RLXEXPD(S)
        S=A22*(S*(TMP-1.D0)-1.D0-C2)+A21*(TMP-1.D0-C1) 
        SUM=SUM+RLPSI2(S,b2)
  120 CONTINUE
      RLWEQTN9=SUM/EN
      RETURN
      END
C
C------------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RLWEQTN10(TAU,Y,NOBS,param)
      implicit double precision (a-h,o-z)
      DIMENSION Y(NOBS),param(4)
      EXTERNAL RLPSI1,RLXEXPD
     
      v = param(1)
      a11 = param(2)
      b1 = param(3)
      c1 = param(4)
      SUM=0.D0
      EN=Dble(NOBS)
      DO 120 I=1,NOBS
        S=(Y(I)-TAU)/V
        TMP=RLXEXPD(S)
        S=A11*(TMP-1.D0-C1) 
        SUM=SUM+RLPSI1(S,b1)
  120 CONTINUE
      RLWEQTN10=SUM/EN
      RETURN
      END


      FUNCTION RLWSCORC(X,IS,C1,C2)
c
c     v*(SCORES CORRIGES)
c     IS=1   v*df/dtau=v*SC1=EXP((X-TAU)/V)-1-C1
c     IS=2   v*df/dv  =v*SC2=((X-TAU)/V)*(EXP((X-TAU)/V)-1)-1-C2
c     N.B: Integrations are computed with the variable z=(x-tau)/v
C
      implicit double precision (a-h,o-z)
      EXTERNAL RLXEXPD
      S=X
      TMP=RLXEXPD(S)
c      GOTO (10,20) IS
      if (IS == 2) then
        goto 20
      else
        goto 10
      end if
   10 RLWSCORC=TMP-1.D0-C1
      RETURN
   20 RLWSCORC=S*(TMP-1.D0)-1.D0-C2
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION RLWSCOR(X,IS)
c
c     IS=1   v*df/dtau=v*SC1=EXP((S-TAU)/V)-1
c     IS=2   v*df/dv  =v*SC2=((S-TAU)/V)*(EXP((S-TAU)/V)-1)-1
c     N.B: Integrations are computed with the variable z=(x-tau)/v
c
      implicit double precision (a-h,o-z)
      EXTERNAL RLXEXPD
      S=X
      TMP=RLXEXPD(S)
c      GOTO (10,20) IS
      if (IS == 2) then
        goto 20
      else
        goto 10
      end if
   10 RLWSCOR=TMP-1.D0
      RETURN
   20 RLWSCOR=S*(TMP-1.D0)-1.D0
      RETURN
      END
C
C=============================================================
C
      FUNCTION RLWZSCOR(X,IZ,A11,A21,A22,C1,C2)
C
C     combined scores
C     iz=1  z1=a11*sc1
C     iz=2  z2=a21*sc1+a22*sc2
c     N.B: Integrations are computed with the variable z=(x-tau)/v
C
      implicit double precision (a-h,o-z)
      EXTERNAL RLWSCORC
      is=1
      sc1=RLwscorc(x,IS,C1,C2)
c      GOTO (10,20) IZ
      if (IZ == 2) then
        goto 20
      else
        goto 10
      end if
   10 RLWZSCOR=a11*sc1
      RETURN
   20 is=2
      sc2=RLwscorc(x,IS,C1,C2)
      RLWZSCOR=a21*sc1+a22*sc2
      RETURN
      END


      FUNCTION RLWDPSI(X,JPSI,JPS0,A11,A21,A22,C1,C2,B1,B2)
C
C     jpsi=1  psi=psi1(z1)
C     jpsi=2  psi=psi2(z2)
C     jps0=-1 psi=-b   jps0=0 psi=z  jps0=1 psi=b
c     N.B: Integrations are computed with the variable z=(x-tau)/v
C
      implicit double precision (a-h,o-z)
      EXTERNAL RLWZSCOR
      IZ=JPSI
      ZS=RLWZSCOR(X,IZ,A11,A21,A22,C1,C2)
      BS=B2
      IF (JPSI.EQ.1) BS=B1
      IF (JPS0.EQ.0) THEN
        RLWDPSI=ZS
      ELSEIF (JPS0.EQ.-1) THEN
        RLWDPSI=-BS
      ELSEIF (JPS0.EQ.1) THEN
        RLWDPSI=BS
      ENDIF
      RETURN
      END
C
C--------------------------------------------------------------
C
      FUNCTION RLWPSIS(DX,WGT,N,EXPSI,EXWLN,TAU,V,
     2         A11,A21,A22,B1,B2,C1,C2,UX12,BETA,YB)

C
C     psi*score
c     IOPT=1  psis=psi1*s1  IOPT=3 psis=psi1*s2
c     IOPT=2  psis=psi2*s1  IOPT=4 psis=psi2*s2
c     N.B: Integrations are computed with the variable z=(s-tau)/v
C
      implicit double precision (a-h,o-z)
      dimension yb(8,2)
      EXTERNAL RLWSCOR,EXPSI,EXWLN
      DIMENSION WGT(N)

      dummy = tau
      dummy = a11
      dummy = a21
      dummy = a22 
      dummy = c1
      dummy = c2
      dummy = beta
      dummy = ux12
      ans=exwln(0.d0,1.d0,dx)
      iopt=int(wgt(1))
      i=int(wgt(2))
      if(iopt.eq.1.or.iopt.eq.3)then
       jpsi=1
       jps0=int(yb(i,1))
       ps1=expsi(dx,jpsi,jps0,a11,a21,a22,c1,c2,b1,b2)
      else
       jpsi=2
       jps0=int(yb(i,2))
       ps2=expsi(dx,jpsi,jps0,a11,a21,a22,c1,c2,b1,b2)
      endif
      is=1
      S1=rlwscor(dx,is)/v
      is=2
      S2=rlwscor(dx,is)/v
c      goto (10,20,30,40) IOPT
      if (IOPT == 2) then
        goto 20
      else if (IOPT == 3) then
        goto 30  
      else if (IOPT == 4) then
        goto 40
      else
        goto 10
      end if
  10  rlwpsis=ps1*s1*ans
      return
  20  rlwpsis=ps2*s1*ans
      return
  30  rlwpsis=ps1*s2*ans
      return
  40  rlwpsis=ps2*s2*ans
      return
      END
C
C--------------------------------------------------------------
C
      FUNCTION RLWPSIPS(DX,WGT,N,EXPSI,EXWLN,TAU,V,
     2     A11,A21,A22,B1,B2,C1,C2,UX12,BETA,YB)
C
C     IOPT=1   WPSIPS=PSI1*PSI1     IOPT=3 PSI2*PSI1
C     N.B: Integrations are computed with the variable z=(s-tau)/v
C

      implicit double precision (a-h,o-z)
      DIMENSION WGT(N),YB(8,2)
      EXTERNAL EXPSI,EXWLN 

      dummy = tau
      dummy = v
      dummy = ux12
      dummy = beta
      ans=exwln(0.d0,1.d0,dx)
      iopt=int(wgt(1))
      i=int(wgt(2))
      jpsi=1
      jps0=int(yb(i,1))
      ps1=expsi(dx,jpsi,jps0,a11,a21,a22,c1,c2,b1,b2)
      jpsi=2
      jps0=int(yb(i,2))
      ps2=expsi(dx,jpsi,jps0,a11,a21,a22,c1,c2,b1,b2)
c      goto (10,20,30,40) iopt
      if (iopt == 2) then
        goto 20
      else if (iopt == 3) then
        goto 30  
      else if (iopt == 4) then
        goto 40
      else
        goto 10
      end if
  10  rlwpsips=ps1*ps1*ans
      return
  20  rlwpsips=ps1*ps2*ans
      return
  30  rlwpsips=ps2*ps1*ans
      return
  40  rlwpsips=ps2*ps2*ans
      return
      end
c
c----------------------------------------------------------------------
c
      SUBROUTINE RLWBRKPT(XLOWER,UPPER,xb,yb,ns,A11,A21,A22,
     *    B1,B2,C1,C2)

      implicit double precision(a-h,o-z)
      DIMENSION XB(8),YB(8,2),Z0(8)
      EXTERNAL RLWZSCOR
      DATA Z0/8*0.D0/

      wlo = xlower
      whi = upper

      XX=-B1/A11+C1+1.D0
      IF (XX.GT.0.D0) THEN
        XB(1)=DLOG(XX)
      ELSE
        XB(1)=xlower-1.D0
      ENDIF
      xx=B1/A11+C1+1.D0
      IF (XX.GT.0.D0) THEN
        xb(2)=dlog(XX)
      ELSE
        xb(2)=xlower-1.D0
      ENDIF
      xb(3)=xlower
      xb(4)=upper 
      call rlsolwx(b2,0.001D0,nsol,xb(5),yb,A21,A22,C1,C2,wlo,whi)
      ns=8
      if (xb(7).eq.0.d0) ns=6
      call rlsrt2(xb,z0,8,1,ns)
      jup=ns
      jl0=0
      do 20 j=1,ns
      if (xb(j).le.xlower) then
        xb(j)=xlower
        jl0=j
      endif
      if (xb(j).ge.upper) then
        xb(j)=upper
        jup=min(jup,j)
      endif
  20  continue
      j0=0
      do 30 j=jl0,jup
      j0=j0+1      
      xb(j0)=xb(j)
  30  continue
      ns=j0
      do 40 j=1,ns-1
      x=(xb(j)+xb(j+1))/2.D0
      xx=x
      iz=1
      zs1=rlwzscor(xx,IZ,A11,A21,A22,C1,C2)
      iz=2
      zs2=rlwzscor(xx,IZ,A11,A21,A22,C1,C2)
      tmp=zs1
      yb(j,1)=0.D0
      yb(j,2)=0.D0
      if (tmp.gt.b1) yb(j,1)=zs1/tmp
      tmp=dabs(zs2)
      if (tmp.gt.b2) yb(j,2)=zs2/tmp
  40  continue
      return
      end
c
c=====================================================================
c
      SUBROUTINE RLVARWEB(ALPHA,SIGMA,TAB,TPAR,TIL,M,Q,MI,V,VSIGA,VMOY)
      implicit double precision(a-h,o-z)
      double precision TAB(5),B(4),M(4),Q(4),MI(4),MIT(4),TPAR(6),
     +     V(4),VSIGA(4),THETA(2),xb(8),yb(8,2)
      EXTERNAL RLGAMDIGAMA 

        W=1.D0/ALPHA
        TAU=DLOG(SIGMA)
        b1=tpar(1)
        b2=tpar(2)
        C1=TAB(1)
        C2=TAB(2)
        A11=TAB(3)
        A21=TAB(4)
        A22=TAB(5)

      TL=TIL 
      ALF1=1.D0+W 
      CALL RLLGAMAD(ALF1,TMP)
      WMU=DEXP(TMP+TAU)
      THETA(1)=WMU/DEXP(TAU) 
      THETA(2)=-DEXP(TAU)*RLGAMDIGAMA(ALF1)*W*W 
      B(1)=DEXP(TAU)
      B(2)=0.D0
      B(3)=0.D0
      B(4)=-1.D0/(W*W)
      CALL RLWEILIM(0.D0,1.D0,WLOW,WHI)
      CALL RLWBRKPT(WLOW,WHI,XB,YB,NS,A11,A21,A22,B1,B2,C1,C2)
C
C                  M    =-int (dpsi/dtheta)=int(psi*s)
c                  Q    = int(psi * psi^T)
C
      CALL RLAUXWAS(TIL,M,Q,tau,w,a11,a21,a22,b1,b2,c1,c2,xb,yb,ns)
C
C  VARIANCE ASS :
c     V(tau,v)     : V    =M^-1 * Q * M^-T
c     V(sigma,alfa): Vsiga=B^T * V * B
c
      CALL RLINVERS(M,MI)
      CALL RLTRNSPO(MI,MIT)
      CALL RLMULTIP(MI,Q,MIT,V)
      CALL RLMULTIP(B,V,B,VSIGA)
      CALL RLMULT2(THETA,VSIGA,THETA,VMOY)
      RETURN
      END

      SUBROUTINE RLAUXWAS(TIL,M,Q,tau,v,a11,a21,a22,b1,b2,c1,c2,
     *           xb,yb,ns)
c
c     N.B: Integrations are computed with the variable z=(x-tau)/v
c          => dx = v * dz 
c
      implicit double precision(a-h,o-z)
      double precision M(4),Q(4),WGT(2),xb(8),yb(8,2)
      dimension iwork(80),work(320)
      EXTERNAL RLWPSIPS,RLWEIBLN,RLWDPSI,RLWPSIS
      u12x11=1.d0
      beta=1.d0
      tild=til
      key=1
      limit=80
C
C  COMPUTE M
C
      DO 100 IOPT=1,4
       WGT(1)=dble(IOPT)
       TS=0.D0
       DO 10 I=1,NS-1
        WGT(2)=I
        CALL RLINTGRW(RLWPSIS,WGT,2,RLWDPSI,RLWEIBLN,XB(I),XB(I+1),
     1    TILD,0.D0,KEY,LIMIT,T,ERRST,NEVAL,IER,WORK,IWORK,
     2      tau,v,a11,a21,a22,b1,b2,c1,c2,u12x11,beta,yb)
        TS=T+TS
  10   CONTINUE
       M(IOPT)=TS
  100 CONTINUE

C
C   COMPUTE Q
C

      DO 200 IOPT=1,4
       WGT(1)=dble(IOPT)
       TS=0.D0
       DO 20 I=1,NS-1
        WGT(2)=I
        CALL RLINTGRW(RLWPSIPS,WGT,2,RLWDPSI,RLWEIBLN,XB(I),XB(I+1),
     1       TILD,0.D0,KEY,LIMIT,T,ERRST,NEVAL,IER,WORK,IWORK,
     2       tau,v,a11,a21,a22,b1,b2,c1,c2,u12x11,beta,yb)
        TS=TS+T
   20  CONTINUE
       Q(IOPT)=TS
  200 CONTINUE
      RETURN
      END

c
c======================================================================
C
      SUBROUTINE RLLIMWBL(SIGMA,ALFA,LOWER,UPPER)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION LOWER
      DATA NCALL,EXMIN,XLGMN,YLGMN,WBLIM/0,0.D0,0.D0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
        CALL RLMACHD(4,XLGMN)
        CALL RLMACHD(5,YLGMN)
        WBLIM=ALOG(1.E-7)
      ENDIF
      LOWER=0.D0
      UPPER=2000.D0
      IF (ALFA.LE.0.2D0) RETURN 
      ALM1=1.D0/ALFA
      CALL RLLGAMAD(1.D0+ALM1,GL)
      ALSIG=DLOG(SIGMA)
      WMU=0.D0
      TMP=GL+ALSIG
      IF (TMP.GT.EXMIN) WMU=DEXP(TMP)
      ALA=DLOG(ALFA)
      CONST=ALA-ALSIG
      IF (ALFA.LE.5.D0) GOTO 250
      X=WMU/2.D0
  100 X=X-0.1D0 
      IF (X.LE.0.D0) GOTO 250
      S=X/SIGMA
      ALOGS=YLGMN
      IF (S.GT.XLGMN) ALOGS=DLOG(S)
      TMP=ALFA*ALOGS
      ALOWEI=CONST+(ALFA-1.D0)*ALOGS
      IF (TMP.GT.EXMIN) ALOWEI=ALOWEI-DEXP(TMP)
      IF (ALOWEI.GT.WBLIM) GOTO 100
      LOWER=X
  250 X=WMU*2.D0
  300 X=X+1.D0
      S=X/SIGMA
      ALOGS=YLGMN
      IF (S.GT.XLGMN) ALOGS=DLOG(S)
      TMP=ALFA*ALOGS
      ALOWEI=CONST+(ALFA-1.D0)*ALOGS
      IF (TMP.GT.EXMIN) ALOWEI=ALOWEI-DEXP(TMP)       
      IF (ALOWEI.GT.WBLIM) GOTO 300
      UPPER=X
      RETURN
      END
C
c
C-----------------------------------------------------------------------
C
      SUBROUTINE RLINTGRW(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,
     *             LIMIT,
     1           RESULT,ABSERR,NEVAL,IER,WORK,IWORK,TAU,V,
     2              A11,A21,A22,B1,B2,C1,c2,U12X11,
     3              BETA,YB)

      implicit double precisioN(a-h, o-z)
      INTEGER ALIST,BLIST,ELIST,RLIST
      DIMENSION FARR(N),WORK(4*LIMIT),IWORK(LIMIT),YB(8,2)
      EXTERNAL F,FEXT,GEXT
C
C         LIMIT IS THE MAXIMUM NUMBER OF SUBINTERVALS ALLOWED IN THE
C         SUBDIVISION PROCESS OF QAGE1D. TAKE CARE THAT LIMIT.GE.1.
C
C**** DATA LIMIT/500/
C
      ALIST=1
      BLIST=ALIST+LIMIT
      RLIST=BLIST+LIMIT
      ELIST=RLIST+LIMIT
      CALL RLQAGE1W(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,LIMIT,
     1     RESULT,ABSERR,NEVAL,IER,
     2     WORK,WORK(BLIST),WORK(RLIST),WORK(ELIST),IWORK,LAST,
     2       TAU,V,
     2       A11,A21,A22,B1,B2,C1,C2,U12X11,
     3       BETA,YB)
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLQAGE1W(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,
     *  LIMIT,
     *  RESULT,ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST,
     *     TAU,V, A11,A21,A22,B1,B2,C1,C2,U12X11,BETA,YB)

      implicit double precision(a-h, o-z)
      dimension yb(8,2)
      DIMENSION ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),IORD(LIMIT),
     *  RLIST(LIMIT),FARR(N)
      EXTERNAL F,FEXT,GEXT

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
      CALL RLMACHD(7,EPMACH)
      CALL RLMACHD(4,UFLOW)
      CALL RLMACHD(6,OFLOW)
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      IER=6
      IER = 0
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      KEYF = KEY
      IF(KEY.LE.0) KEYF = 1
      IF(KEY.GE.7) KEYF = 6
      C = dble(FLOAT(KEYF))
      NEVAL = 0
      IF (KEYF.EQ.1)
     *  CALL RLQ1K15W(F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,DE FABS,
     *  RESABS,TAU,V,
     2              A11,A21,A22,B1,B2,C1,C2,U12X11,
     3              BETA,yb)

      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
C
C           TEST ON ACCURACY.
C
      ERRBND = DMAX1(EPSABS,EPSREL*DABS(RESULT))
      IF(ABSERR.LE.5.0D+01*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
C LIMIT==1 MEANS THAT ONLY 1 SUBINTERVAL IS ALLOWED..
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS)
     *  .OR.ABSERR.EQ.0.0D+00) GO TO 60
C
C           INITIALIZATION
C           --------------
C
C
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      NRMAX = 1
      IROFF1 = 0
      IROFF2 = 0
C
C           MAIN DO-LOOP
C           ------------
C
      DO 30 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
C
        AA1 = ALIST(MAXERR)
        BB1 = 5.0D-01*(ALIST(MAXERR)+BLIST(MAXERR))
        AA2 = BB1
        BB2 = BLIST(MAXERR)
        IF (KEYF.EQ.1)
     *  CALL RLQ1K15W(F,FARR,N,FEXT,GEXT,AA1,BB1,AREA1,ERROR1,
     *          RESABS,DEFAB1,TAU,V,
     2              A11,A21,A22,B1,B2,C1,C2,U12X11,
     3              BETA,yb)
        IF (KEYF.EQ.1)
     *  CALL RLQ1K15W(F,FARR,N,FEXT,GEXT,AA2,BB2,AREA2,ERROR2,
     *          RESABS,DEFAB2,TAU,V,
     2              A11,A21,A22,B1,B2,C1,C2,U12X11,
     3              BETA,yb)
        NEVAL = NEVAL+1
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 5
        IF(DABS(RLIST(MAXERR)-AREA12).LE.1.0D-05*DABS(AREA12)
     *  .AND.ERRO12.GE.9.9D-01*ERRMAX) IROFF1 = IROFF1+1
        IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF2 = IROFF2+1
    5   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = DMAX1(EPSABS,EPSREL*DABS(AREA))
        IF(ERRSUM.LE.ERRBND) GO TO 8
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
        IF(IROFF1.GE.6.OR.IROFF2.GE.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
C           SUBINTERVALS EQUALS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE.
C
        IF(DMAX1(DABS(AA1),DABS(BB2)).LE.(1.0D+00+C*1.0D+03*
     *  EPMACH)*(DABS(AA2)+1.0D+04*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
    8   IF(ERROR2.GT.ERROR1) GO TO 10
        ALIST(LAST) = AA2
        BLIST(MAXERR) = BB1
        BLIST(LAST) = BB2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 20
   10   ALIST(MAXERR) = AA2
        ALIST(LAST) = AA1
        BLIST(LAST) = BB1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE QSORTD TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH THE LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
   20   CALL RLQSORTD(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C***JUMP OUT OF DO-LOOP
        IF(IER.NE.0.OR.ERRSUM.LE.ERRBND) GO TO 40
   30 CONTINUE
C
C           COMPUTE FINAL RESULT.
C           ---------------------
C
   40 RESULT = 0.0D+00
      DO 50 K=1,LAST
        RESULT = RESULT+RLIST(K)
   50 CONTINUE
      ABSERR = ERRSUM
   60 IF(KEYF.NE.1) NEVAL = (10*KEYF+1)*(2*NEVAL+1)
      IF(KEYF.EQ.1) NEVAL = 30*NEVAL+15
      CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLQ1K15W
     *  (F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,RESABS,RESASC,TAU,V,
     2              A11,A21,A22,B1,B2,C1,C2,U12X11,
     3              BETA,yb)

      implicit double precision(a-h, o-z)
      EXTERNAL F,FEXT,GEXT
      DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8),FARR(N),yb(8,2)
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
      CALL RLMACHD(7,EPMACH)
      CALL RLMACHD(4,UFLOW)
      CALL RLMACHD(6,OFLOW)

      CENTR = 5.0D-01*(A+B)
      HLGTH = 5.0D-01*(B-A)
      DHLGTH = DABS(HLGTH)
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR,FARR,N,FEXT,GEXT,TAU,V,
     2              A11,A21,A22,B1,B2,C1,C2,U12X11,
     3              BETA,yb)

      RESG = FC*WG(4)
      RESK = FC*WGK(8)
      RESABS = DABS(RESK)
      DO 10 J=1,3
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC,FARR,N,FEXT,GEXT,TAU,V,
     2              A11,A21,A22,B1,B2,C1,C2,U12X11,
     3              BETA,yb)

        FVAL2 = F(CENTR+ABSC,FARR,N,FEXT,GEXT,TAU,V,
     2              A11,A21,A22,B1,B2,C1,C2,U12X11,
     3              BETA,yb)

        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(DABS(FVAL1)+DABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,4
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC,FARR,N,FEXT,GEXT,TAU,V,
     2              A11,A21,A22,B1,B2,C1,C2,U12X11,
     3              BETA,yb)
        FVAL2 = F(CENTR+ABSC,FARR,N,FEXT,GEXT,TAU,V,
     2              A11,A21,A22,B1,B2,C1,C2,U12X11,
     3              BETA,yb)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(DABS(FVAL1)+DABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*5.0D-01
      RESASC = WGK(8)*DABS(FC-RESKH)
      DO 20 J=1,7
        RESASC = RESASC+WGK(J)*(DABS(FV1(J)-RESKH)+DABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = DABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     *  ABSERR = RESASC*DMIN1(1.0D+00,(2.0D+02*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(5.0D+01*EPMACH)) ABSERR = DMAX1
     *  ((EPMACH*5.0D+01)*RESABS,ABSERR)
      RETURN
      END




 









