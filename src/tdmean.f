c$$$      SUBROUTINE RLMACHD(I,X)
c$$$C.......................................................................
c$$$C
c$$$C
c$$$C                   DOUBLE PRECISION VERSION
c$$$C                   ************************
c$$$C
c$$$C  MACHINE PARAMETERS : TO  ALTER  THIS  SUBROUTINE  FOR  A PARTICULAR
c$$$C  ++++++++++++++++++   ENVIRONMENT, THE DESIRED SET OF DATA STATEMENT
c$$$C  SHOULD BE ACTIVATED BY REMOVING THE "C" FROM COLUMN ONE  AND ADDING
c$$$C  THE "C" FOR THE TWO LINES AFTER "... VAX FORTRAN (V5) compiler".
c$$$C
c$$$C  EPMACH IS APPROXIMATELY EQUAL TO THE EXPONENT PART OF PREC
c$$$C  EXMIN, XLGMN, YLGMN AND XBIG CAN BE FOUND BY TRIAL AND ERROR
c$$$C
c$$$C  RADIX IS ALWAYS EQUAL TO 2.D0
c$$$C  PREC CAN BE FOUND BY CALLING THE ROBETH SUBROUTINE "PRECD"
c$$$      DOUBLE PRECISION X,RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c$$$C  for VAX FORTRAN (V5) compiler (VMS)
c$$$C       DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c$$$C     *  /2.D0,1.38778D-17,-88.722D0,2.939D-39,-88.7227D0,1.7D38,1.D-17/
c$$$C  for IBM-PC F77L compiler
c$$$C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c$$$C    *  /2.D0,5.422D-20,-708.D0,1.D-307,-706.591D0,1.D308,1.D-19/
c$$$C  for IBM-PC MICROSOFT FORTRAN compiler
c$$$C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c$$$C    *  /2.D0,5.47522D-18,-745.133D0,0.9D-48,-110.629D0,1.D308,1.D-17/
c$$$C  for WATCOM F77 Compiler (32 bits version)
c$$$      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c$$$     +     /2.,0.1121D-15,-709.782D0,9.74D-290,-718.433D0,1.797D308,
c$$$     +     1.0D-17/
c$$$C    */2.,0.1121D-15,-709.782D0,0.974D-312,-718.433D0,1.797D308,1.0D-17/
c$$$C  for ULTRIX DEC FORTRAN-77 compiler
c$$$C     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c$$$C    */2.,1.12133D-16,-744.44D0,2.226D-308,-708.396D0,1.797D308,1.0D-17/
c$$$C  for SUN FORTRAN compiler
c$$$C     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c$$$C    */2.,1.12133D-16,-745.13D0,0.494D-323,-744.44D0,1.797D308,1.0D-17/
c$$$C  for SILICON GRAPHICS MIPS FORTRAN 77 Compiler
c$$$C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c$$$C    */2.,1.12133D-16,-744.04D0,0.758D-323,-743.75D0,1.797D308,1.0D-17/
c$$$C  for HP-UX FORTRAN 77 Compiler
c$$$C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c$$$C    */2.,1.12133D-16,-708.396D0,0.1D-308,-709.09D0,1.797D308,1.0D-17/
c$$$C  for DEC-ALPHA FORTRAN Compiler (OpenVMS)
c$$$C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c$$$C    */2.,1.1102D-16,-709.782D0,0.1057D-45,-105.863D0,0.898D307,1.0D-17/
c$$$C
c$$$      IF (I.EQ.1) X=RADIX
c$$$      IF (I.EQ.2) X=PREC
c$$$      IF (I.EQ.3) X=EXMIN
c$$$      IF (I.EQ.4) X=XLGMN
c$$$      IF (I.EQ.5) X=YLGMN
c$$$      IF (I.EQ.6) X=XBIG
c$$$      IF (I.EQ.7) X=EPMACH
c$$$      RETURN
c$$$      END
c$$$C
c$$$C==========================================================================
C
c$$$      DOUBLE PRECISION FUNCTION RLEXU(X)
c$$$      DOUBLE PRECISION X
c$$$      RLEXU=X
c$$$      RETURN
c$$$      END
C
C-----------------------------------------------------------------------
C
c$$$      DOUBLE PRECISION FUNCTION RLXLOGD(X)
c$$$C.......................................................................
c$$$C
c$$$C   COPYRIGHT 1992 Alfio Marazzi
c$$$C
c$$$C   AUTHOR : A. RANDRIAMIHARISOA
c$$$C.......................................................................
c$$$      implicit double precision(a-h,o-z)
c$$$      DATA NCALL,XMIN,YMIN/0,0.D0,0.D0/
c$$$      IF (NCALL.EQ.0) THEN
c$$$        CALL RLMACHD(4,XMIN)
c$$$        CALL RLMACHD(5,YMIN)
c$$$        NCALL=1
c$$$      ENDIF
c$$$C
c$$$C  EXTENDED NATURAL LOGARITHM FUNCTION
c$$$C
c$$$      IF (X.LE.0.D0) THEN
c$$$        RLXLOGD=0.D0
c$$$      ELSEIF (X.LE.XMIN) THEN
c$$$        RLXLOGD=YMIN
c$$$      ELSE
c$$$        RLXLOGD=DLOG(X)
c$$$      ENDIF
c$$$      RETURN
c$$$      END
c$$$C
c$$$C-----------------------------------------------------------------------
c$$$C
c$$$      DOUBLE PRECISION FUNCTION RLXEXPD(X)
c$$$C.......................................................................
c$$$C
c$$$C   COPYRIGHT 1992 Alfio Marazzi
c$$$C
c$$$C   AUTHOR : A. RANDRIAMIHARISOA
c$$$C.......................................................................
c$$$      DOUBLE PRECISION X,DMIN,DMAX,XBIG
c$$$      DATA NCALL,DMIN,DMAX,XBIG/0,0.D0,0.D0,0.D0/
c$$$      IF (NCALL.EQ.0) THEN
c$$$        CALL RLMACHD(3,DMIN)
c$$$        CALL RLMACHD(6,XBIG)
c$$$        XBIG=XBIG/10.D0
c$$$        DMAX=DLOG(XBIG)
c$$$        NCALL=1
c$$$      ENDIF
c$$$C
c$$$C  EXTENDED EXPONENTIAL FUNCTION
c$$$C
c$$$      IF (X.LE.DMIN) THEN
c$$$        RLXEXPD=0.D0
c$$$      ELSEIF (X.GE.DMAX) THEN
c$$$        RLXEXPD=XBIG
c$$$      ELSE
c$$$        RLXEXPD=DEXP(X)
c$$$      ENDIF
c$$$      RETURN
c$$$      END
c$$$C
      SUBROUTINE RLXERF(KODE,X,P)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      implicit double precision(a-h,o-z)
      EXTERNAL RLXEXPD
      DATA SPI/2.506628274631D0/
c      IF (KODE.NE.1.AND.KODE.NE.2) CALL MESSGE(500,'XERF  ',1)
      X2=-X*X/2.D0
      P=RLXEXPD(X2)
      IF (KODE.EQ.2) P=P/SPI
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLGAUSSD (KODE,X,P)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C
C.......................................................................
C
      DOUBLE PRECISION   P,X,SQR1D2,CD
      DATA               SQR1D2/.7071067811865475D0/
C
c      IF (KODE.NE.1.AND.KODE.NE.2) CALL MESSGE(500,'GAUSSD',1)
      CALL RLCERFD(-X*SQR1D2,CD)
      P = .5D0 * CD
      IF (KODE.EQ.2) P=1.D0-P
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLCERFD(X,F)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C
C.......................................................................
C
      DOUBLE PRECISION   F,X,RLXEXPD
      DIMENSION          P(5),Q(4),P1(9),Q1(8),P2(6),Q2(5)
      DOUBLE PRECISION   P,Q,P1,Q1,P2,Q2,XMIN,XLARGE,SQRPI,XX,
     *                   RES,XSQ,XNUM,XDEN,XI,XBIG
      INTEGER            ISW,I
      EXTERNAL RLXEXPD
C                                  COEFFICIENTS FOR 0.0 .LE. Y .LT.
C                                  .477
      DATA               P(1)/113.8641541510502D0/,
     *                   P(2)/377.4852376853020D0/,
     *                   P(3)/3209.377589138469D0/,
     *                   P(4)/.1857777061846032D0/,
     *                   P(5)/3.161123743870566D0/
      DATA               Q(1)/244.0246379344442D0/,
     *                   Q(2)/1282.616526077372D0/,
     *                   Q(3)/2844.236833439171D0/,
     *                   Q(4)/23.60129095234412D0/
C                                  COEFFICIENTS FOR .477 .LE. Y
C                                  .LE. 4.0
      DATA               P1(1)/8.883149794388376D0/,
     *                   P1(2)/66.11919063714163D0/,
     *                   P1(3)/298.6351381974001D0/,
     *                   P1(4)/881.9522212417691D0/,
     *                   P1(5)/1712.047612634071D0/,
     *                   P1(6)/2051.078377826071D0/,
     *                   P1(7)/1230.339354797997D0/,
     *                   P1(8)/2.153115354744038D-8/,
     *                   P1(9)/.5641884969886701D0/
      DATA               Q1(1)/117.6939508913125D0/,
     *                   Q1(2)/537.1811018620099D0/,
     *                   Q1(3)/1621.389574566690D0/,
     *                   Q1(4)/3290.799235733460D0/,
     *                   Q1(5)/4362.619090143247D0/,
     *                   Q1(6)/3439.367674143722D0/,
     *                   Q1(7)/1230.339354803749D0/,
     *                   Q1(8)/15.74492611070983D0/
C                                  COEFFICIENTS FOR 4.0 .LT. Y
      DATA               P2(1)/-3.603448999498044D-01/,
     *                   P2(2)/-1.257817261112292D-01/,
     *                   P2(3)/-1.608378514874228D-02/,
     *                   P2(4)/-6.587491615298378D-04/,
     *                   P2(5)/-1.631538713730210D-02/,
     *                   P2(6)/-3.053266349612323D-01/
      DATA               Q2(1)/1.872952849923460D0/,
     *                   Q2(2)/5.279051029514284D-01/,
     *                   Q2(3)/6.051834131244132D-02/,
     *                   Q2(4)/2.335204976268692D-03/,
     *                   Q2(5)/2.568520192289822D0/
C                                  CONSTANTS
      DATA               XMIN/1.0D-10/,XLARGE/6.375D0/
C                                  CERFD(XBIG) .APPROX. DETAP
      DATA               XBIG/13.3D0/
      DATA               SQRPI/.5641895835477563D0/
C
      Y=real(X)
      XX = Y
      ISW = 1
      IF (XX.GE.0.0D0) GO TO 5
      ISW = -1
      XX = -XX
    5 IF (XX.LT..477D0) GO TO 10
      IF (XX.LE.4.0D0) GO TO 30
      IF (ISW .GT. 0) GO TO 40
      IF (XX.LT.XLARGE) GO TO 45
      RES = 2.0D0
      GO TO 70
C                                  ABS(Y) .LT. .477, EVALUATE
C                                  APPROXIMATION FOR CERFD
   10 IF (XX.LT.XMIN) GO TO 20
      XSQ = XX*XX
      XNUM = P(4)*XSQ+P(5)
      XDEN = XSQ+Q(4)
      DO 15 I = 1,3
         XNUM = XNUM*XSQ+P(I)
         XDEN = XDEN*XSQ+Q(I)
   15 CONTINUE
      RES = XX*XNUM/XDEN
      GO TO 25
   20 RES = XX*P(3)/Q(3)
   25 IF (ISW.EQ.-1) RES = -RES
      RES = 1.0D0-RES
      GO TO 70
C                                  .477 .LE. ABS(Y) .LE. 4.0
C                                  EVALUATE APPROXIMATION FOR CERFD
   30 XSQ = XX*XX
      XNUM = P1(8)*XX+P1(9)
      XDEN = XX+Q1(8)
      DO 35 I=1,7
         XNUM = XNUM*XX+P1(I)
         XDEN = XDEN*XX+Q1(I)
   35 CONTINUE
      RES = XNUM/XDEN
      GO TO 60
C                                  4.0 .LT. ABS(Y), EVALUATE
C                                  MINIMAX APPROXIMATION FOR CERFD
   40 IF (XX.GT.XBIG) GO TO 65
   45 XSQ = XX*XX
      XI = 1.0D0/XSQ
      XNUM= P2(5)*XI+P2(6)
      XDEN = XI+Q2(5)
      DO 50 I = 1,4
         XNUM = XNUM*XI+P2(I)
         XDEN = XDEN*XI+Q2(I)
   50 CONTINUE
      RES = (SQRPI+XI*XNUM/XDEN)/XX
   60 RES = RES*RLXEXPD(-XSQ)
      IF (ISW.EQ.-1) RES = 2.0D0-RES
      GO TO 70
   65 RES = 0.0D0
   70 F = RES
      RETURN
      END
c
C-----------------------------------------------------------------------
C
      SUBROUTINE RLINTGRT(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,
     1           LIMIT,RESULT,ABSERR,NEVAL,IER,WORK,IWORK,NPR,PARAM)

C.......................................................................
C
C   AUTHOR :     QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
      implicit double precision(a-h, o-z)
      DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,RESULT,FEXT,WORK
      INTEGER IER,KEY,LAST,LIMIT,NEVAL,ALIST,BLIST,ELIST,RLIST
C
      DIMENSION FARR(N),WORK(4*LIMIT),IWORK(LIMIT),PARAM(NPR)
C
      EXTERNAL F,FEXT,GEXT
C
C         LIMIT IS THE MAXIMUM NUMBER OF SUBINTERVALS ALLOWED IN THE
C         SUBDIVISION PROCESS OF QAGE1D. TAKE CARE THAT LIMIT.GE.1.
C
C**** DATA LIMIT/500/
C
C      IF ((EPSABS.LT.0..AND.EPSREL.LT.0.).OR.LIMIT.LE.1
C     1   .OR.LIMIT.GT.500)  CALL MESSGE(500,'INTGRD',1)
      ALIST=1
      BLIST=ALIST+LIMIT
      RLIST=BLIST+LIMIT
      ELIST=RLIST+LIMIT
      CALL RLQAGE1T(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,LIMIT,
     1     RESULT,ABSERR,NEVAL,IER,
     2     WORK,WORK(BLIST),WORK(RLIST),WORK(ELIST),IWORK,LAST,
     3     NPR,PARAM)
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLQAGE1T(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,
     *     LIMIT,
     *     RESULT,ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST,
     *     NPR,PARAM)
C.......................................................................
C
C   AUTHOR :     QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
      implicit double precision(a-h, o-z)
C      DOUBLE PRECISION A,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,
C     *  AA1,AA2,B,
C     *  BLIST,BB1,BB2,C,DABS,DEFABS,DEFAB1,DEFAB2,DMAX1,ELIST,EPMACH,
C     *  EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,F,OFLOW,
C     *  RESABS,RESULT,RLIST,UFLOW,FEXT
      INTEGER IER,IORD,IROFF1,IROFF2,K,KEY,KEYF,LAST,LIMIT,MAXERR,NEVAL,
     *  NRMAX
C
      DIMENSION ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),IORD(LIMIT),
     *  RLIST(LIMIT),FARR(N),param(npr)
C
      EXTERNAL F,FEXT,GEXT
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
C      IF ((EPSABS.LT.0..AND.EPSREL.LT.0.).OR.LIMIT.LE.1
C     1   .OR.LIMIT.GT.500)  CALL MESSGE(500,'QAGE1D',1)
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
     *  CALL RLQ1K15T(F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,DEFABS,
     *  RESABS,NPR,PARAM)

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
     *  CALL RLQ1K15T(F,FARR,N,FEXT,GEXT,AA1,BB1,AREA1,ERROR1,
     *          RESABS,DEFAB1,NPR,PARAM)
        IF (KEYF.EQ.1)
     *  CALL RLQ1K15T(F,FARR,N,FEXT,GEXT,AA2,BB2,AREA2,ERROR2,
     *          RESABS,DEFAB2,NPR,PARAM)
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
C  999 IF (IER.NE.0) CALL MESSGE(400+IER,'QAGE1D',0)
      CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLQ1K15T
     *  (F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,RESABS,RESASC,NPR,PARAM)
C.......................................................................
C
C   AUTHOR :     QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
      implicit double precision(a-h, o-z)
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DABS,DHLGTH,DMAX1,DMIN1,
     *  EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,OFLOW,RESABS,RESASC,
     *  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK,FEXT
      INTEGER J,JTW,JTWM1
      EXTERNAL F,FEXT,GEXT
C
      DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8),FARR(N),param(npr)
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
C
      CENTR = 5.0D-01*(A+B)
      HLGTH = 5.0D-01*(B-A)
      DHLGTH = DABS(HLGTH)
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR,FARR,N,FEXT,GEXT,NPR,PARAM)

      RESG = FC*WG(4)
      RESK = FC*WGK(8)
      RESABS = DABS(RESK)
      DO 10 J=1,3
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)

        FVAL2 = F(CENTR+ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)

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
        FVAL1 = F(CENTR-ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)
        FVAL2 = F(CENTR+ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)
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
c$$$c
c$$$c***********************************************************************
c$$$c
c$$$      SUBROUTINE RLQSORTD(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
c$$$C.......................................................................
c$$$C
c$$$C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
c$$$C
c$$$C   PROGRAMMER : QUADPACK
c$$$C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
c$$$C.......................................................................
c$$$C
c$$$      DOUBLE PRECISION ELIST,ERMAX,ERRMAX,ERRMIN
c$$$      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
c$$$     *  NRMAX
c$$$      DIMENSION ELIST(LAST),IORD(LAST)
c$$$C
c$$$C           CHECK WHETHER THE LIST CONTAINS MORE THAN
c$$$C           TWO ERROR ESTIMATES.
c$$$C
c$$$C***FIRST EXECUTABLE STATEMENT
c$$$      IF(LAST.GT.2) GO TO 10
c$$$      IORD(1) = 1
c$$$      IORD(2) = 2
c$$$      GO TO 90
c$$$C
c$$$C           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
c$$$C           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
c$$$C           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
c$$$C           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
c$$$C
c$$$   10 ERRMAX = ELIST(MAXERR)
c$$$      IF(NRMAX.EQ.1) GO TO 30
c$$$      IDO = NRMAX-1
c$$$      DO 20 I = 1,IDO
c$$$        ISUCC = IORD(NRMAX-1)
c$$$C***JUMP OUT OF DO-LOOP
c$$$        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
c$$$        IORD(NRMAX) = ISUCC
c$$$        NRMAX = NRMAX-1
c$$$   20    CONTINUE
c$$$C
c$$$C           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE
c$$$C           MAINTAINED IN DESCENDING ORDER. THIS NUMBER
c$$$C           DEPENDS ON THE NUMBER OF SUBDIVISIONS STILL ALLOWED.
c$$$C
c$$$   30 JUPBN = LAST
c$$$      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
c$$$      ERRMIN = ELIST(LAST)
c$$$C
c$$$C           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
c$$$C           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
c$$$C
c$$$      JBND = JUPBN-1
c$$$      IBEG = NRMAX+1
c$$$      IF(IBEG.GT.JBND) GO TO 50
c$$$      DO 40 I=IBEG,JBND
c$$$        ISUCC = IORD(I)
c$$$C***JUMP OUT OF DO-LOOP
c$$$        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
c$$$        IORD(I-1) = ISUCC
c$$$   40 CONTINUE
c$$$   50 IORD(JBND) = MAXERR
c$$$      IORD(JUPBN) = LAST
c$$$      GO TO 90
c$$$C
c$$$C           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
c$$$C
c$$$   60 IORD(I-1) = MAXERR
c$$$      K = JBND
c$$$      DO 70 J=I,JBND
c$$$        ISUCC = IORD(K)
c$$$C***JUMP OUT OF DO-LOOP
c$$$        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
c$$$        IORD(K+1) = ISUCC
c$$$        K = K-1
c$$$   70 CONTINUE
c$$$      IORD(I) = LAST
c$$$      GO TO 90
c$$$   80 IORD(K+1) = LAST
c$$$C
c$$$C           SET MAXERR AND ERMAX.
c$$$C
c$$$   90 MAXERR = IORD(NRMAX)
c$$$      ERMAX = ELIST(MAXERR)
c$$$      RETURN
c$$$      END
c$$$C
c$$$C-----------------------------------------------------------------------
c$$$C
c$$$      SUBROUTINE RLLGAMAD(X,GL)
c$$$C.......................................................................
c$$$C
c$$$C   AUTHORS :     M.C. PIKE AND I.D. HILL (1966)
c$$$C                 ALGORITHM 291: LOGARITHM OF GAMMA FUNCTION.
c$$$C                 COMMUNICATIONS OF THE ACM, VOL.9, P 684.
c$$$C                 ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
c$$$C.......................................................................
c$$$C
c$$$      DOUBLE PRECISION X,GL,V,F,Z
c$$$C      IF (X.LE.0.D0) CALL MESSGE(500,'LGAMAD',1)
c$$$      V=X
c$$$      F=0.D0
c$$$      IF (X.GE.7.D0) GOTO 300
c$$$      F=1.D0
c$$$      Z=X-1.D0
c$$$  100 Z=Z+1.D0
c$$$      IF (Z.GE.7.D0) GOTO 200
c$$$      V=Z
c$$$      F=F*Z
c$$$      GOTO 100
c$$$  200 V=V+1.D0
c$$$      F=-DLOG(F)
c$$$  300 Z=1.D0/V**2
c$$$      GL=F+(V-0.5D0)*DLOG(V)-V+.9189385332D0+(((-.000595238D0*Z+
c$$$     +   .0007936507D0)*Z - .0027777778D0)*Z+.0833333333D0)/V
c$$$      RETURN
c$$$      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLINGAMA(X,P,G)
C.......................................................................
C
C   AUTHOR :     G. P. BHATTACHARJEE (1970)
C                ALGORITHM AS 32 "THE INCOMPLETE GAMA INTEGRAL"
C                APPLIED STATISTICS, VOL.19.
C                REPRINT FROM PP.285-287 WITH THE PERMISSION OF
C                BLACKWELL PUBLISHERS.
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
      implicit double precision(a-h,o-z)
      DIMENSION PN(6)
      EXTERNAL RLXEXPD
      DATA TOL/1.0D-8/
C
      G=0.D0
      IF (X.EQ.0.D0) RETURN
C      IF (X.LT.0..OR.P.LE.0.) CALL MESSGE(500,'INGAMA',1)
      CALL RLMACHD(6,OFLO)
      OFLO=OFLO*1.D-15
      CALL RLLGAMAD(P,GP)
      GIN=0.D0
      FACTOR=RLXEXPD(P*DLOG(X)-X-GP)
      IF (X.GT.1.D0.AND.X.GE.P) GOTO 30
C
C  CALCULATION BY SERIES EXPANSION
C
      GIN=1.D0
      TERM=1.D0
      RN=P
   20 RN=RN+1.D0
      TERM=TERM*X/RN
      GIN=GIN+TERM
      IF (TERM.GT.TOL) GOTO 20
      GIN=GIN*FACTOR/P
      GOTO 50
C
C  CALCULATION BY CONTINUED FRACTION
C
   30 A=1.D0-P
      B=A+X+1.D0
      TERM=0.D0
      PN(1)=1.D0
      PN(2)=X
      PN(3)=X+1.D0
      PN(4)=X*B
      GIN=PN(3)/PN(4)
   32 A=A+1.D0
      B=B+2.D0
      TERM=TERM+1.D0
      AN=A*TERM
C      DO 33 I=1,2
C   33 PN(I+4)=B*PN(I+2)-AN*PN(I)
      DO I=1,2
        PN(I+4)=B*PN(I+2)-AN*PN(I)
      END DO
      IF (PN(6).EQ.0.D0) GOTO 35
      RN=PN(5)/PN(6)
      DIF=DABS(GIN-RN)
      IF (DIF.GT.TOL) GOTO 34
      IF (DIF.LE.TOL*RN) GOTO 42
   34 GIN=RN
C   35 DO 36 I=1,4
C   36 PN(I)=PN(I+2)
   35 DO I=1,4
        PN(I)=PN(I+2)
      END DO
      IF (DABS(PN(5)).LT.OFLO) GOTO 32
C      DO 41 I=1,4
C   41 PN(I)=PN(I)/OFLO
      DO I=1,4
        PN(I)=PN(I)/OFLO
      END DO
      GOTO 32
   42 GIN=1.D0-FACTOR*GIN
   50 G=GIN
      RETURN
      END
C
c$$$      DOUBLE PRECISION FUNCTION RLGAMDIGAMA(X)
c$$$C.......................................................................
c$$$C
c$$$C   DI-GAMMA(X) FOR X.GT.0.
c$$$C   METHOD :
c$$$C   SEE FORMULAS 6.3.5 AND 6.3.18, ABRAMOWITZ AND STEGUN, 1970
c$$$C   AUTHOR : G. VAN MELLE
c$$$C.......................................................................
c$$$C
c$$$      implicit double precision(a-h,o-z)
c$$$c      DOUBLE PRECISION T,Z,Z2,S
c$$$c      IF (X.LE.0.) CALL MESSGE(500,'DIGAMA',1)
c$$$      T=0.D0
c$$$      Z=X
c$$$      IF (X.LT.5.D0) THEN
c$$$        N=5-INT(X+1.D-10)
c$$$        DO 100 K=1,N
c$$$          T=T+1.D0/Z
c$$$          Z=Z+1.D0
c$$$  100   CONTINUE
c$$$      ENDIF
c$$$      Z2=1.D0/(Z*Z)
c$$$      S=DLOG(Z) - 1.D0/(2.D0*Z) + Z2*(- 1.D0/12.D0 + Z2*(1.D0/120.D0 +
c$$$     + Z2*(-1.D0/252.D0 + Z2*(1.D0/240.D0 + Z2*(- 1.D0/132.D0 + Z2*(
c$$$     + 691.D0/(12.D0*2730.D0) - Z2/12.D0)))))) - T
c$$$      RLGAMDIGAMA=S
c$$$      RETURN
c$$$      END
C
      SUBROUTINE RLDIGAMA(X,RESULT)
      DOUBLE PRECISION X,RESULT,RLGAMDIGAMA
      EXTERNAL RLGAMDIGAMA
      RESULT=RLGAMDIGAMA(X)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLSRT1(A,N,K1,K2)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      DOUBLE PRECISION A(N),X
c      LOGICAL NPRCHK
C
c      NPRCHK=K1.GE.1.AND.K2.GT.K1.AND.K2.LE.N
c      IF (.NOT.NPRCHK) CALL MESSGE(500,'SRT1  ',1)
      N1=K2-K1+1
c      I=1
c   10 I=I+I
c      IF (I.LE.N1) GOTO 10
      M=N1
   20 M=M/2
      IF (M.EQ.0) GOTO 90
      K=N1-M
      DO 40 J=1,K
      L=J
   50 IF (L.LT.1) GOTO 40
      LPM=L+M
      LPM1=LPM+K1-1
      L1=L+K1-1
      IF (A(LPM1).GE.A(L1)) GOTO 40
      X=A(LPM1)
      A(LPM1)=A(L1)
      A(L1)=X
      L=L-M
      GOTO 50
   40 CONTINUE
      GOTO 20
   90 CONTINUE
      END
c$$$C
c$$$C-----------------------------------------------------------------------
c$$$C
c$$$      SUBROUTINE RLSRT2(A,B,N,K1,K2)
c$$$C.......................................................................
c$$$C
c$$$C   COPYRIGHT 1992 Alfio Marazzi
c$$$C
c$$$C   AUTHOR : A. MARAZZI
c$$$C.......................................................................
c$$$C
c$$$      DOUBLE PRECISION A(N),B(N),X,Y
c$$$c      LOGICAL NPRCHK
c$$$C
c$$$c      NPRCHK=N.GT.0.AND.K1.GE.1.AND.K2.GE.K1.AND.K2.LE.N
c$$$c      IF (.NOT.NPRCHK) CALL MESSGE(500,'SRT2  ',1)
c$$$      N1=K2-K1+1
c$$$c      I=1
c$$$c   10 I=I+I
c$$$c      IF (I.LE.N1) GOTO 10
c$$$      M=N1
c$$$   20 M=M/2
c$$$      IF (M.EQ.0) GOTO 90
c$$$      K=N1-M
c$$$      DO 40 J=1,K
c$$$      L=J
c$$$   50 IF (L.LT.1) GOTO 40
c$$$      LPM=L+M
c$$$      LPM1=LPM+K1-1
c$$$      L1=L+K1-1
c$$$      IF (A(LPM1).GE.A(L1)) GOTO 40
c$$$      X=A(LPM1)
c$$$      Y=B(LPM1)
c$$$      A(LPM1)=A(L1)
c$$$      B(LPM1)=B(L1)
c$$$      A(L1)=X
c$$$      B(L1)=Y
c$$$      L=L-M
c$$$      GOTO 50
c$$$   40 CONTINUE
c$$$      GOTO 20
c$$$   90 CONTINUE
c$$$      END
C***
C======================================================================
C***
      SUBROUTINE RLQUANTD(P,XX)
C**********************************************************************
C     Normal distribution inverse
C
C     The  rational   function   on  page 95    of Kennedy  and  Gentle,
C     Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
C     value for the Newton method of finding roots.
C
C     If P or 1-P .lt. machine EPS returns +/- xx(EPS)
C**********************************************************************
      DOUBLE PRECISION p,q,eps,xx,strtx,xcur,cum,pp,dx,
     +       rlstvaln,rldennor
      LOGICAL qporq
      EXTERNAL rlstvaln,rldennor
      DATA NCALL,MAXIT,EPS/0,100,0.D0/
      IF (NCALL.EQ.0) THEN
        CALL RLMACHD(2,EPS)
        NCALL=1
      ENDIF
C
C     FIND MINIMUM OF P AND Q
C
      q=1.D0-p
      qporq = p .LE. q
      IF (.NOT. (qporq)) GO TO 10
      pp = p
      GO TO 20

   10 pp = q
C
C     INITIALIZATION STEP
C
   20 strtx = rlstvaln(pp)
      xcur = strtx
C
C     NEWTON INTERATIONS
C
      DO 30,i = 1,maxit
          call RLGAUSSD(1,xcur,cum)
          dx = (cum-pp)/rldennor(xcur)
          xcur = xcur - dx
          IF (DABS(dx/xcur).LT.eps) GO TO 40
   30 CONTINUE
      xx = strtx
C
C     IF WE GET HERE, NEWTON HAS FAILED
C
      IF (.NOT.qporq) xx = -xx
      RETURN
C
C     IF WE GET HERE, NEWTON HAS SUCCEDED
C
   40 xx = xcur
      IF (.NOT.qporq) xx = -xx
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION rldennor(x)
      DOUBLE PRECISION x,r2pi,nhalf
      DATA  r2pi,nhalf/0.3989422804014326D0,-0.5D0/
      rldennor = r2pi*dexp(nhalf*x*x)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION rlstvaln(p)
C**********************************************************************
C                    STarting VALue for Neton-Raphon
C                calculation of Normal distribution Inverse
C**********************************************************************
      DOUBLE PRECISION p,sign,y,z,xden(5),xnum(5),rldevlpl
      EXTERNAL rldevlpl
      DATA xnum/-0.322232431088D0,-1.000000000000D0,-0.342242088547D0,
     +     -0.204231210245D-1,-0.453642210148D-4/
      DATA xden/0.993484626060D-1,0.588581570495D0,0.531103462366D0,
     +     0.103537752850D0,0.38560700634D-2/
C     ..
      IF (.NOT. (p.LE.0.5D0)) GO TO 10
      sign = -1.0D0
      z = p
      GO TO 20

   10 sign = 1.0D0
      z = 1.0D0 - p
   20 y = dsqrt(-2.0D0*dlog(z))
      rlstvaln = y + rldevlpl(xnum,5,y)/rldevlpl(xden,5,y)
      rlstvaln = sign*rlstvaln
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION rldevlpl(a,n,x)
C**********************************************************************
C     returns  A(1) + A(2)*X + ... + A(N)*X**(N-1)
C**********************************************************************
      DOUBLE PRECISION x,a(n),term
      term = a(n)
      DO 10,i=n-1,1,-1
        term = a(i) + term*x
   10 CONTINUE
      rldevlpl = term
      RETURN
      END

C ----------------------------------------------------------------------
      SUBROUTINE RLQGAMMA (P, A, X)
C ----------------------------------------------------------------------
C            INVERSE INCOMPLETE GAMMA RATIO FUNCTION
C
C     GIVEN A>0 , AND NONEGATIVE P THEN X IS COMPUTED WHERE P(A,X) = P.
C     SCHRODER ITERATION IS EMPLOYED.
C                      ------------
C
C     X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0,
C     OTHERWISE, RLQGAMMA ATTEMPTS TO OBTAIN A SOLUTION FOR P(A,X) = P.
C     IF THE ROUTINE IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X.
C
C     MESSAGE THAT REPORTS THE ERROR STATUS OF THE RESULTS.
C
C       MESSAGE 403   NO SOLUTION WAS OBTAINED. THE RATIO (1-P)/A
C                     IS TOO LARGE.
C       MESSAGE 407   ITERATION FAILED. NO VALUE IS GIVEN FOR X.
C                     THIS MAY OCCUR WHEN X IS APPROXIMATELY 0.
C       MESSAGE 408   A VALUE FOR X HAS BEEN OBTAINED, BUT THE
C                     ROUTINE IS NOT CERTAIN OF ITS ACCURACY.
C                     ITERATION CANNOT BE PERFORMED IN THIS
C                     CASE. THIS CAN OCCUR IF P OR 1-P IS
C                     APPROXIMATELY 0.
C ----------------------------------------------------------------------
C     WRITTEN BY ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WEAPONS CENTER
C        DAHLGREN, VIRGINIA
C     ADAPTED FOR ROBETH BY A.RANDRIAMIHARISOA
C     -------------------
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION LN10, EPS0(2),AMIN(2),BMIN(2),DMIN(2),EMIN(2)
      EXTERNAL RLXEXPD,RLXLOGD
C     -------------------
C     LN10 = LN(10)
C     C = EULER CONSTANT
C     -------------------
      DATA LN10 /2.302585D0/
      DATA C  /.577215664901533D0/
C     -------------------
      DATA A0 /3.31125922108741D0/, A1 /11.6616720288968D0/,
     *     A2 /4.28342155967104D0/, A3 /.213623493715853D0/
      DATA B1 /6.61053765625462D0/, B2 /6.40691597760039D0/,
     *     B3 /1.27364489782223D0/, B4 /.036117081018842D0/
C     -------------------
      DATA EPS0(1) /1.D-10/, EPS0(2) /1.D-08/
      DATA AMIN(1) / 500.D0/, AMIN(2) / 100.D0/
      DATA BMIN(1) /1.D-28/, BMIN(2) /1.D-13/
      DATA DMIN(1) /1.D-06/, DMIN(2) /1.D-04/
      DATA EMIN(1) /2.D-03/, EMIN(2) /6.D-03/
C     -------------------
      DATA NCALL,TOL,E,EXMIN,XMIN,XMAX,XL,YL/0,1.D-5,6*0.D0/
C     -------------------
C     ****** E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS.
C            E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0.
C            XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE
C            LARGEST POSITIVE NUMBER.
C
       IF (NCALL.EQ.0) THEN
         NCALL=1
         CALL RLMACHD(2,E)
         CALL RLMACHD(3,EXMIN)
         XMIN=DEXP(EXMIN)
         CALL RLMACHD(4,XL)
         CALL RLMACHD(5,YL)
         CALL RLMACHD(6,XMAX)
       ENDIF
      X0=0.D0
      Q=1.D0-P
      X = 0.D0
c      IF (A.LE.0..OR.P.LT.0..OR.P.GT.1.) CALL MESSGE(500,'RLQGAMMA',1)
C     T = DBLE(P) + DBLE(Q) - 1.D0
C     IF (ABS(T) .GT. E) GO TO 520
C
      IERR = 0
      IF (P .EQ. 0.D0) RETURN
      IF (Q .EQ. 0.D0) GO TO 400
      IF (DABS(A-1.D0) .LE. 1.D-6) GO TO 410
C
      E2 = 2.D0*E
      AMAX = 0.4D-10/(E*E)
      IOP = 1
      IF (E .GT. 1.D-10) IOP = 2
      EPS = EPS0(IOP)
      XN = X0
      IF (X0 .GT. 0.D0) GO TO 100
C
C        SELECTION OF THE INITIAL APPROXIMATION XN OF X
C                       WHEN A .LT. 1
C
      IF (A .GT. 1.D0) GO TO 30
      CALL RLLGAMAD(A+1.D0,GL)
      G = RLXEXPD(GL)
C     G = GAMMA(A + 1.0)
      QG = Q*G
      IF (QG .EQ. 0.D0) GO TO 560
      B = QG/A
      IF (QG .GT. 0.6D0*A) GO TO 20
      IF (A .GE. 0.30D0 .OR. B .LT. 0.35D0) GO TO 10
         T = RLXEXPD(-(B + C))
         U = T*RLXEXPD(T)
         XN = T*RLXEXPD(U)
         GO TO 100
C
   10 IF (B .GE. 0.45D0) GO TO 20
      IF (B .EQ. 0.D0) GO TO 560
C      Y = -RLXLOGD(B,XL,YL)
      Y = -RLXLOGD(B)
      S = 0.5D0 + (0.5D0 - A)
C      Z = RLXLOGD(Y,XL,YL)
      Z = RLXLOGD(Y)
      T = Y - S*Z
      IF (B .LT. 0.15D0) GO TO 11
C         XN = Y - S*RLXLOGD(T,XL,YL) - RLXLOGD(1.D0+S/(T+1.D0),XL,YL)
         XN = Y - S*RLXLOGD(T) - RLXLOGD(1.D0+S/(T+1.D0))
         GO TO 200
   11 IF (B .LE. 1.D-2) GO TO 12
         U = ((T + 2.D0*(3.D0 - A))*T + (2.D0 - A)*(3.D0 - A))/
     *           ((T + (5.D0 - A))*T + 2.D0)
C         XN = Y - S*RLXLOGD(T,XL,YL) - RLXLOGD(U,XL,YL)
         XN = Y - S*RLXLOGD(T) - RLXLOGD(U)
         GO TO 200
   12 C1 = -S*Z
      C2 = -S*(1.D0 + C1)
      C3 =  S*((0.5D0*C1 + (2.D0 - A))*C1 + (2.5D0 - 1.5D0*A))
      C4 = -S*(((C1/3.D0 + (2.5D0 - 1.5D0*A))*C1 + ((A - 6.D0)*A +
     *       7.D0))*C1 + ((11.D0*A - 46.D0)*A + 47.D0)/6.D0)
      C5 = -S*((((-C1/4.D0 + (11.D0*A - 17.D0)/6.D0)*C1
     *           + ((-3.D0*A + 13.D0)*A - 13.D0))*C1
     *           + 0.5D0*(((2.D0*A - 25.D0)*A + 72.D0)*A - 61.D0))*C1
     *           + (((25.D0*A - 195.D0)*A + 477.D0)*A -379.D0)/12.D0)
      XN = ((((C5/Y + C4)/Y + C3)/Y + C2)/Y + C1) + Y
      IF (A .GT. 1.D0) GO TO 200
      IF (B .GT. BMIN(IOP)) GO TO 200
      X = XN
      RETURN
C
   20 IF (B*Q .GT. 1.D-8) GO TO 21
         XN = RLXEXPD(-(Q/A + C))
         GO TO 23
   21 IF (P .LE. 0.9D0) GO TO 22
         CALL RLLGAMAD(1.D0+A,GL)
C         XN = RLXEXPD((RLXLOGD(1.D0-Q,XL,YL) + GL)/A)
         XN = RLXEXPD((RLXLOGD(1.D0-Q) + GL)/A)
         GO TO 23
   22 XN = DEXP(DLOG(P*G)/A)
   23 IF (XN .EQ. 0.D0) GO TO 510
      T = 0.5D0 + (0.5D0 - XN/(A + 1.D0))
      XN = XN/T
      GO TO 100
C
C        SELECTION OF THE INITIAL APPROXIMATION XN OF X
C                       WHEN A .GT. 1
C
   30 IF (Q .LE. 0.5D0) GO TO 31
C         W = RLXLOGD(P,XL,YL)
         W = RLXLOGD(P)
         GO TO 32
C   31 W = RLXLOGD(Q,XL,YL)
   31 W = RLXLOGD(Q)
   32 T = SQRT(-2.D0*W)
      S = T - (((A3*T + A2)*T + A1)*T + A0)/((((B4*T + B3)*T
     *                + B2)*T + B1)*T + 1.D0)
      IF (Q .GT. 0.5D0) S = -S
C
      RTA = SQRT(A)
      S2 = S*S
      XN = A + S*RTA + (S2 - 1.D0)/3.D0 + S*(S2 - 7.D0)/(36.D0*RTA)
     1       - ((3.D0*S2 + 7.D0)*S2 - 16.D0)/(810.D0*A)
     2       + S*((9.D0*S2 + 256.D0)*S2 - 433.D0)/(38880.D0*A*RTA)
      XN = DMAX1(XN, 0.D0)
      IF (A .LT. AMIN(IOP)) GO TO 40
      X = XN
      D = 0.5D0 + (0.5D0 - X/A)
      IF (DABS(D) .LE. DMIN(IOP)) RETURN
C
   40 IF (P .LE. 0.5D0) GO TO 50
      IF (XN .LT. 3.D0*A) GO TO 200
      CALL RLLGAMAD(A,GL)
      Y = -(W + GL)
      D = DMAX1(2.D0, A*(A - 1.D0))
      IF (Y .LT. LN10*D) GO TO 41
         S = 1.D0 - A
         Z = DLOG(Y)
         GO TO 12
   41 T = A - 1.D0
C      XN = Y + T*RLXLOGD(XN,XL,YL) - RLXLOGD(1.D0-T/(XN + 1.D0),XL,YL)
      XN = Y + T*RLXLOGD(XN) - RLXLOGD(1.D0-T/(XN + 1.D0))
C      XN = Y + T*RLXLOGD(XN,XL,YL) - RLXLOGD(1.D0-T/(XN + 1.D0),XL,YL)
      XN = Y + T*RLXLOGD(XN) - RLXLOGD(1.D0-T/(XN + 1.D0))
      GO TO 200
C
   50 AP1 = A + 1.D0
      IF (XN .GT. 0.7D0*AP1) GO TO 101
      CALL RLLGAMAD(AP1,GL)
      W = W + GL
      IF (XN .GT. 0.15D0*AP1) GO TO 60
         AP2 = A + 2.D0
         AP3 = A + 3.D0
         X = RLXEXPD((W + X)/A)
C         TX = RLXLOGD(1.D0 + (X/AP1)*(1.D0 + X/AP2),XL,YL)
         TX = RLXLOGD(1.D0 + (X/AP1)*(1.D0 + X/AP2))
         X = RLXEXPD((W + X - TX)/A)
C         TX = RLXLOGD(1.D0 + (X/AP1)*(1.D0 + X/AP2),XL,YL)
         TX = RLXLOGD(1.D0 + (X/AP1)*(1.D0 + X/AP2))
         X = RLXEXPD((W + X - TX)/A)
C         TX = RLXLOGD(1.D0+(X/AP1)*(1.D0+(X/AP2)*(1.D0+X/AP3)),XL,YL)
         TX = RLXLOGD(1.D0+(X/AP1)*(1.D0+(X/AP2)*(1.D0+X/AP3)))
         X = RLXEXPD((W + X - TX)/A)
         XN = X
         IF (XN .GT. 1.D-2*AP1) GO TO 60
         IF (XN .LE. EMIN(IOP)*AP1) RETURN
         GO TO 101
C
   60 APN = AP1
      T = XN/APN
      SUM = 1.D0 + T
   61    APN = APN + 1.D0
         T = T*(XN/APN)
         SUM = SUM + T
         IF (T .GT. 1.D-4) GO TO 61
C      T = W - RLXLOGD(SUM,XL,YL)
      T = W - RLXLOGD(SUM)
      XN = RLXEXPD((XN + T)/A)
C      XN = XN*(1.D0 - (A*RLXLOGD(XN,XL,YL) - XN -T)/(A - XN))
      XN = XN*(1.D0 - (A*RLXLOGD(XN) - XN -T)/(A - XN))
      GO TO 101
C
C                 SCHRODER ITERATION USING P
C
  100 IF (P .GT. 0.5D0) GO TO 200
  101 IF (P .LE. 1.D10*XMIN) GO TO 550
      AM1 = (A - 0.5D0) - 0.5D0
  102 IF (A .LE. AMAX) GO TO 110
      D = 0.5D0 + (0.5D0 - XN/A)
      IF (DABS(D) .LE. E2) GO TO 550
C
  110 IF (IERR .GE. 20) GO TO 530
      IERR = IERR + 1
      CALL RLINGAMA (XN, A, PN)
      QN=1.D0-PN
      IF (PN .EQ. 0.D0 .OR. QN .EQ. 0.D0) GO TO 550
      R = RLRCOMP(A,XN)
      IF (R .EQ. 0.D0) GO TO 550
      T = (PN - P)/R
      W = 0.5D0*(AM1 - XN)
      IF (DABS(T) .LE. 0.1D0 .AND. DABS(W*T) .LE. 0.1D0) GO TO 120
         X = XN*(1.D0 - T)
         IF (X .LE. 0.D0) GO TO 540
         D = DABS(T)
         GO TO 121
C
  120 H = T*(1.D0 + W*T)
      X = XN*(1.D0 - H)
      IF (X .LE. 0.D0) GO TO 540
      IF (DABS(W) .GE. 1.D0 .AND. DABS(W)*T*T .LE. EPS) RETURN
      D = DABS(H)
  121 XN = X
      IF (D .GT. TOL) GO TO 102
      IF (D .LE. EPS) RETURN
      IF (DABS(P - PN) .LE. TOL*P) RETURN
      GO TO 102
C
C                 SCHRODER ITERATION USING Q
C
  200 IF (Q .LE. 1.D10*XMIN) GO TO 550
      AM1 = (A - 0.5D0) - 0.5D0
  201 IF (A .LE. AMAX) GO TO 210
      D = 0.5D0 + (0.5D0 - XN/A)
      IF (DABS(D) .LE. E2) GO TO 550
C
  210 IF (IERR .GE. 20) GO TO 530
      IERR = IERR + 1
      CALL RLINGAMA (XN, A, PN)
      QN=1.D0-PN
      IF (PN .EQ. 0.D0 .OR. QN .EQ. 0.D0) GO TO 550
      R = RLRCOMP(A,XN)
      IF (R .EQ. 0.D0) GO TO 550
      T = (Q - QN)/R
      W = 0.5D0*(AM1 - XN)
      IF (DABS(T) .LE. 0.1D0 .AND. DABS(W*T) .LE. 0.1D0) GO TO 220
         X = XN*(1.D0 - T)
         IF (X .LE. 0.D0) GO TO 540
         D = DABS(T)
         GO TO 221
C
  220 H = T*(1.D0 + W*T)
      X = XN*(1.D0 - H)
      IF (X .LE. 0.D0) GO TO 540
      IF (DABS(W) .GE. 1.D0 .AND. DABS(W)*T*T .LE. EPS) RETURN
      D = DABS(H)
  221 XN = X
      IF (D .GT. TOL) GO TO 201
      IF (D .LE. EPS) RETURN
      IF (DABS(Q - QN) .LE. TOL*Q) RETURN
      GO TO 201
C
C                       SPECIAL CASES
C
  400 X = XMAX
      RETURN
C
C  410 X = -RLXLOGD(1.-P,XL,YL)
  410 X = -RLXLOGD(1.-P)
      RETURN
C
  510 X=-403
c     CALL MESSGE(403,'RLQGAMMA',0)
      RETURN
C
  530 X=-406
c     CALL MESSGE(406,'RLQGAMMA',0)
      RETURN
C
  540 X=-407
c     CALL MESSGE(407,'RLQGAMMA',0)
      RETURN
C
  550 X = XN
c      CALL MESSGE(408,'RLQGAMMA',0)
      RETURN
C
  560 X = XMAX
c      CALL MESSGE(408,'RLQGAMMA',0)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLREXP(X)
      implicit double precision(a-h,o-z)
C     ------------------------------------------------------------------
C     COMPUTATION OF EXP(X) - 1
C     ------------------------------------------------------------------
      DATA P1/ .914041914819518D-09/, P2/ .238082361044469D-01/,
     *     Q1/-.499999999085958D+00/, Q2/ .107141568980644D+00/,
     *     Q3/-.119041179760821D-01/, Q4/ .595130811860248D-03/
C     ------------------
      IF (DABS(X) .GT. 0.15D0) GO TO 10
      RLREXP = X*(((P2*X + P1)*X + 1.D0)/((((Q4*X + Q3)*X + Q2)*X
     *                 + Q1)*X + 1.D0))
      RETURN
C
   10 W = DEXP(X)
      IF (X .GT. 0.D0) GO TO 20
         RLREXP = (W - 0.5D0) - 0.5D0
         RETURN
   20 RLREXP = W*(0.5D0 + (0.5D0 - 1.D0/W))
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLRLOG(X)
      implicit double precision(a-h,o-z)
C     -------------------
C     COMPUTATION OF  X - 1 - LN(X)
C     -------------------
      DATA A/.566749439387324D-01/
      DATA B/.456512608815524D-01/
C     -------------------
      DATA P0/ .333333333333333D+00/, P1/-.224696413112536D+00/,
     *     P2/ .620886815375787D-02/
      DATA Q1/-.127408923933623D+01/, Q2/ .354508718369557D+00/
C     -------------------
      IF (X .LT. 0.61D0 .OR. X .GT. 1.57D0) GO TO 100
      IF (X .LT. 0.82D0) GO TO 10
      IF (X .GT. 1.18D0) GO TO 20
C
C              ARGUMENT REDUCTION
C
      U = (X - 0.5D0) - 0.5D0
      W1 = 0.D0
      GO TO 30
C
   10 U = DBLE(X) - 0.7D0
      U = U/0.7D0
      W1 = A - U*0.3D0
      GO TO 30
C
   20 U = 0.75D0*DBLE(X) - 1.D0
      W1 = B + U/3.D0
C
C               SERIES EXPANSION
C
   30 R = U/(U + 2.D0)
      T = R*R
      W = ((P2*T + P1)*T + P0)/((Q2*T + Q1)*T + 1.D0)
      RLRLOG = 2.D0*T*(1.D0/(1.D0 - R) - R*W) + W1
      RETURN
C
C
  100 R = (X - 0.5D0) - 0.5D0
      RLRLOG = R - DLOG(X)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLGAM1(A)
C     ------------------------------------------------------------------
C     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5
C     ------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION P(7), Q(5), R(9)
C     -------------------
      DATA P(1)/ .577215664901533D+00/, P(2)/-.409078193005776D+00/,
     *     P(3)/-.230975380857675D+00/, P(4)/ .597275330452234D-01/,
     *     P(5)/ .766968181649490D-02/, P(6)/-.514889771323592D-02/,
     *     P(7)/ .589597428611429D-03/
C     -------------------
      DATA Q(1)/ .100000000000000D+01/, Q(2)/ .427569613095214D+00/,
     *     Q(3)/ .158451672430138D+00/, Q(4)/ .261132021441447D-01/,
     *     Q(5)/ .423244297896961D-02/
C     -------------------
      DATA R(1)/-.422784335098468D+00/, R(2)/-.771330383816272D+00/,
     *     R(3)/-.244757765222226D+00/, R(4)/ .118378989872749D+00/,
     *     R(5)/ .930357293360349D-03/, R(6)/-.118290993445146D-01/,
     *     R(7)/ .223047661158249D-02/, R(8)/ .266505979058923D-03/,
     *     R(9)/-.132674909766242D-03/
C     -------------------
      DATA S1  / .273076135303957D+00/, S2  / .559398236957378D-01/
C     -------------------
      T = A
      D = A - 0.5D0
      IF (D .GT. 0.D0) T = D - 0.5D0
C      IF (T) 30,10,20
      if(T < 0.D0) then
        goto 30
      else if (T > 0.D0) then
        goto 20
      else
        goto 10
      end if
C
   10 RLGAM1 = 0.D0
      RETURN
C
   20 TOP = (((((P(7)*T + P(6))*T + P(5))*T + P(4))*T + P(3))*T
     *                  + P(2))*T + P(1)
      BOT = (((Q(5)*T + Q(4))*T + Q(3))*T + Q(2))*T + 1.D0
      W = TOP/BOT
      IF (D .GT. 0.D0) GO TO 21
         RLGAM1 = A*W
         RETURN
   21 RLGAM1 = (T/A)*((W - 0.5D0) - 0.5D0)
      RETURN
C
   30 TOP = (((((((R(9)*T + R(8))*T + R(7))*T + R(6))*T + R(5))*T
     *                    + R(4))*T + R(3))*T + R(2))*T + R(1)
      BOT = (S2*T + S1)*T + 1.D0
      W = TOP/BOT
      IF (D .GT. 0.D0) GO TO 31
         RLGAM1 = A*((W + 0.5D0) + 0.5D0)
         RETURN
   31 RLGAM1 = T*W/A
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION RLRCOMP(A,X)
      implicit double precision(a-h,o-z)
C     -------------------
C     EVALUATION OF EXP(-X)*X**A/GAMMA(A)
C     -------------------
C     RT2PIN = 1/SQRT(2*PI)
C     -------------------
      DATA RT2PIN/.398942280401433D0/
C     -------------------
      RLRCOMP = 0.D0
      IF (A .GE. 20.D0) GO TO 20
      T = A*DLOG(X) - X
      IF (A .GE. 1.D0) GO TO 10
         RLRCOMP = (A*DEXP(T))*(1.D0 + RLGAM1(A))
         RETURN
   10 CALL RLLGAMAD(A,AL)
      RLRCOMP = RLXEXPD(T-AL)  !/GAMMA(A)
      RETURN
C
   20 U = X/A
      IF (U .EQ. 0.D0) RETURN
      T = (1.D0/A)**2
      T1 = (((0.75D0*T - 1.D0)*T + 3.5D0)*T - 105.D0)/(A*1260.D0)
      T1 = T1 - A*RLRLOG(U)
      RLRCOMP = RT2PIN*DSQRT(A)*RLXEXPD(T1)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE RLPWEIBL(ALPHA,SIGMA,X,P)
      implicit double precision(a-h,o-z)
      DATA NCALL,EXMIN,XLGMN,YLGMN/0,3*0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
        CALL RLMACHD(4,XLGMN)
        CALL RLMACHD(5,YLGMN)
      ENDIF
c     IF (ALPHA.LE.0..OR.SIGMA.LE.0.) CALL MESSGE(500,'PWEIBL',1)
      P=0.D0
      IF (X.LE.0.D0) RETURN
      ALGXS=YLGMN
      XS=X/SIGMA
      IF (XS.GT.XLGMN) ALGXS=DLOG(XS)
      T=ALPHA*ALGXS
      EXPT=0.D0
      IF (T.GT.EXMIN) EXPT=DEXP(T)
      XXPT=0.D0
      IF (-EXPT.GT.EXMIN) XXPT=DEXP(-EXPT)
      P=1.D0-XXPT
      RETURN
      END
C
      SUBROUTINE RLQWEIBL(ALPHA,SIGMA,P,X)
      implicit double precision(a-h,o-z)
      DATA NCALL,EXMIN,XLGMN,YLGMN,XMAX/0,4*0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
        CALL RLMACHD(4,XLGMN)
        CALL RLMACHD(5,YLGMN)
        CALL RLMACHD(6,XMAX)
      ENDIF
c      IF (ALPHA.LE.0.D0.OR.SIGMA.LE.0..OR.P.LT.0..OR.P.GT.1.)
c     +CALL MESSGE(500,'QWEIBL',1)
      X=0.D0
      IF (P.LE.0.D0) RETURN
      X=XMAX
      IF (P.GE.1.D0) RETURN
      T=1.D0-P
      ALOGT=YLGMN
      IF (T.GT.XLGMN) ALOGT=DLOG(T)
      ALLGT=YLGMN
      IF (-ALOGT.GT.XLGMN) ALLGT=DLOG(-ALOGT)
      T=ALLGT/ALPHA
      XT=0.D0
      IF (XT.GT.EXMIN) XT=DEXP(T)
      X=XT*SIGMA
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE RLPLNORM(ALPHA,SIGMA,X,P)
      implicit double precision(a-h,o-z)
      DATA NCALL,XLGMN,YLGMN/0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(4,XLGMN)
        CALL RLMACHD(5,YLGMN)
      ENDIF
c      IF (SIGMA.LE.0.D0.OR.X.LT.0.D0) CALL MESSGE(500,'PLNORM',1)
      P=0.D0
      IF (X.LE.0.D0) RETURN
      ALOGX=YLGMN
      IF (X.GT.YLGMN) ALOGX=DLOG(X)
      Z=(ALOGX-ALPHA)/SIGMA
      CALL RLGAUSSD(1,Z,P)
      RETURN
      END
C
      SUBROUTINE RLQLNORM(ALPHA,SIGMA,P,X)
      implicit double precision(a-h,o-z)
c      DOUBLE PRECISION QQ
      EXTERNAL RLXEXPD
      DATA NCALL,EXMIN,XMAX/0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
        CALL RLMACHD(6,XMAX)
      ENDIF
c      IF (SIGMA.LE.0.D0.OR.P.LT.0.D0.OR.P.GT.1.) CALL MESSGE(500,'QLNORM',1)
      X=0.D0
      IF (P.LE.0.D0) RETURN
      X=XMAX
      IF (P.GE.1.D0) RETURN
      CALL RLQUANTD(DBLE(P),QQ)
      T=ALPHA+DBLE(QQ*SIGMA)
      X=RLXEXPD(T)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE RLTMEANE(X,N,BETA,POS)
      implicit double precision(a-h,o-z)
      DIMENSION X(N)
      CALL RLSRT1(X,N,1,N)
      EN=DBLE(FLOAT(N))
      IF (DABS(BETA-0.5D0).LT.1.D-5) THEN
        MED=INT(EN/2.)
        POS=X(MED+1)
        IF (2*MED.EQ.N) POS=(X(MED)+POS)/2.D0
      ELSEIF (BETA.LT.1.D-5) THEN
        SUM=X(1)
C        DO 100 I=2,N
C  100   SUM=SUM+X(I)
        DO I=2,N
          SUM=SUM+X(I)
        END DO
        POS=SUM/EN
      ELSE
        IU=INT(EN*(1.D0-BETA))
        IL=INT(EN*BETA)
        SU=EN*(1.D0-BETA) - FLOAT(IU)
        SL=EN*BETA - FLOAT(IL)
        I=0
        SUMU=0.D0
        SUML=0.D0
        IUL=MAX0(IU,IL)
  200   I=I+1
        IF (I.LE.IU) SUMU=SUMU+X(I)
        IF (I.LE.IL) SUML=SUML+X(I)
        IF (I.LT.IUL) GOTO 200
        AREA=SUMU+SU*X(IU+1)-SUML-SL*X(IL+1)
        UNITS=SU-SL+FLOAT(IU-IL)
        POS=AREA/UNITS
      ENDIF
      RETURN
      END
C
      SUBROUTINE RLTMADVE(X,N,BETA,GAM,POS,SCAL,SX)
      implicit double precision(a-h,o-z)
      DIMENSION X(N),SX(N)
      CALL RLTMEANE(X,N,BETA,POS)
C      DO 100 I=1,N
C  100 SX(I)=DABS(X(I)-POS)
      DO I=1,N
        SX(I)=DABS(X(I)-POS)
      END DO
      CALL RLTMEANE(SX,N,GAM,SCAL)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE RLTRMNG(ALPHA,SIGMA,BETA,mF)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION mF
      IF (DABS(BETA-0.5D0).LT.1.D-5) THEN
        CALL RLQGAMMA(0.5D0,ALPHA,Q)
        mF=SIGMA*Q
        RETURN
      ENDIF
      IF (BETA.LT.1.D-5) THEN
        mF=SIGMA*ALPHA
      ELSE
        CALL RLQGAMMA(BETA,ALPHA,QL)
        CALL RLQGAMMA(1.D0-BETA,ALPHA,QU)
        CALL RLINGAMA(QL,ALPHA+1.D0,GGL)
        CALL RLINGAMA(QU,ALPHA+1.D0,GGU)
        mF=SIGMA*(GGU-GGL)*ALPHA/(1.D0-2.D0*BETA)
      ENDIF
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLEQAD1G(DD,V,NV,PARAM)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION DD,V(NV),M,MMD,PARAM(1)
      ALPHA=V(1)
      GAM  =V(3)
      M    =V(4)
      D=DBLE(DD)
      DUMMY=PARAM(1)
      CALL RLINGAMA(M+D,ALPHA,P1)
      MMD=M-D
      IF (MMD.LT.0.D0) MMD=0.D0
      CALL RLINGAMA(MMD,ALPHA,P2)
      TMP=P1-P2-(1.D0-GAM)
      RLEQAD1G=DBLE(TMP)
      RETURN
      END
C
      SUBROUTINE RLQAD1DG(ALPHA,BETA,GAM,TOL,QAD1,ISOL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION RLEQAD1G,QAD,SMIN,SMAX,PARAM(1)
      DIMENSION Q(4)
      EXTERNAL RLEQAD1G
c      COMMON/DATI/NQ
      QAD1=0.D0
      ISOL=0
      NQ=4
      Q(1)=ALPHA
      Q(2)=BETA
      Q(3)=GAM
      CALL RLTRMNG(ALPHA,1.D0,BETA,Q(4))
      CALL RLQGAMMA(0.999D0,ALPHA,UPPER)
      SMIN=0.D0
      SMAX=DBLE(UPPER)
c      TOLD=DBLE(TOL)
      CALL RLRGFLD(RLEQAD1G,Q,0.D0,SMIN,SMAX,TOL,100,QAD,ITERM,NQ,PARAM)
c     IF (ITERM.NE.1) RETURN
      QAD1=DBLE(QAD)
      ISOL=1
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLASLVDG(AA,V)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION V(5),M
      POS =V(1)
      SCAL=V(2)
      BETA=V(3)
      GAM =V(4)
      TOL =V(5)
      ALPHA=DBLE(AA)
      CALL RLTRMADG(ALPHA,1.D0,BETA,GAM,TOL,M,D)
      TMP=M/D - POS/SCAL
      RLASLVDG=DBLE(TMP)
      RETURN
      END
C
      SUBROUTINE RLSOLVDG(POS,SCAL,BETA,GAM,ALFA1,ALFA2,TOL,ALPHA,ISOL)
      implicit double precision(a-h,o-z)
c      DOUBLE PRECISION ASLVDG,A,B,FA,FB,TL,TOLD,X,XN,FN
      DIMENSION V(5)
      EXTERNAL RLASLVDG
      DATA TL/1.D-10/
      ALPHA=0.D0
      ISOL=0
      V(1)=POS
      V(2)=SCAL
      V(3)=BETA
      V(4)=GAM
      V(5)=TOL
      A=DBLE(ALFA1)
      B=DBLE(ALFA2)
      TOLD=DBLE(TOL)
      IF (A.GE.100.D0) TOLD=DMIN1(5.D-3,TOLD)
      MAXIT=100
C     WARNING:
C     A CALL TO RGFLD HERE CAUSES RECURSION AND WRONG RESULTS!!
C     CALL RLRGFLD(ASLVDG,V,0.D0,A,B,TOLD,MAXIT,X,ITERM)
      ITR=1
      FA=RLASLVDG(A,V)
      FB=RLASLVDG(B,V)
C
C  REGULA FALSI ITERATION
C
   20 IF (DABS(FA-FB).GT.TL) GOTO 30
c      CALL MESSGE(401,'RGFLD ',0)
      RETURN
   30 XN=(A*FB-B*FA)/(FB-FA)
      TOLD=DBLE(TOL)
      IF (XN.GE.100.D0) TOLD=DMIN1(5.D-3,TOLD)
      FN=RLASLVDG(XN,V)
      IF (ITR.GE.MAXIT) GOTO 60
      IF (DABS(FN).LT.TOLD) GOTO 70
      IF (FA*FN.LE.0.D0) GOTO 40
      A=XN
      FA=FN
      GOTO 50
   40 B=XN
      FB=FN
   50 ITR=ITR+1
      GOTO 20
C
   60 ITERM=2
      X=XN
      ALPHA=DBLE(X)
      IF (DABS(FN).LT.10.D0*TOLD) GOTO 70
      RETURN
   70 ITERM=1
      X=XN
      ALPHA=DBLE(X)
      ISOL=1
      RETURN
      END
C
      SUBROUTINE RLTRMADG(ALPHA,SIGMA,BETA,GAM,TOL,mF,sF)
      implicit double precision(a-h,o-z)
c      DOUBLE PRECISION AA,A0,PA0,B0,A1,PA1,B1,A2,PA2,B2,SS,SUP,SLOW
      DOUBLE PRECISION mF,MMQ
      CALL RLTRMNG(ALPHA,SIGMA,BETA,mF)
      IF (DABS(GAM-0.5D0).LT.1.D-5) THEN
        CALL RLQAD1DG(ALPHA,BETA,GAM,TOL,QAD1,ISOL)
        sF=SIGMA*QAD1
      ELSE
        GAM0=GAM
        AA=DBLE(ALPHA)
        SS=DBLE(SIGMA)
        CALL RLINGAMA(DBLE(mF)/SS,AA+1.D0,PA0)
        A0=PA0*AA*SS
        CALL RLINGAMA(DBLE(mF)/SS,AA,B0)
 100    CALL RLQAD1DG(ALPHA,BETA,GAM0,TOL,QAD1,ISOL)
        QADF=SIGMA*QAD1
        CALL RLINGAMA(DBLE(mF + QADF)/SS,AA+1.D0,PA1)
        A1=PA1*AA*SS
        MMQ=(mF-QADF)/SIGMA
        IF (MMQ.LT.0.D0) MMQ=0.D0
        CALL RLINGAMA(DBLE(MMQ),AA+1.D0,PA2)
        A2=PA2*AA*SS
        CALL RLINGAMA(DBLE(mF + QADF)/SS,AA,B1)
        CALL RLINGAMA(DBLE(MMQ),AA,B2)
        IF (DABS(GAM-GAM0).LT.1.D-6) THEN
          SUP=(A1 + A2 - 2.D0*A0) - DBLE(mF) * (B1 + B2 - 2.D0*B0)
          GAM0=1.D0-GAM
          GOTO 100
        ELSE
          SLOW=(A1 + A2 - 2.D0*A0) - DBLE(mF) * (B1 + B2 - 2.D0*B0)
        ENDIF
        sF=DBLE(SUP-SLOW)/(1.D0-2.D0*GAM)
      ENDIF
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION RLLEQNG(LL,V,NV,PARAM)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION LL, V(NV),PARAM(1)
      DATA U0,ALF0,F0/3*0.D0/
      U=V(1)
      ALPHA=V(2)
      EL=DBLE(LL)
      DUMMY=PARAM(1)
      IF (DABS(EL-U).LT.0.0001D0) EL=EL-0.0002D0
      IF (DABS(U-U0).LT.1.D-5.AND.DABS(ALPHA-ALF0).LT.1.D-5) THEN
        F2=F0
        X=EL
        CALL RLQGAMMA(X,ALPHA,QX)
      ELSE
        X=U
        U0=U
        ALF0=ALPHA
        QX=V(3)
      ENDIF
  100 IF (X.LT.1D-5) THEN
        FX=0.D0
      ELSE
        IF (1.D0-X.LT.1D-5) THEN
          FX=1.D0
        ELSE
          CALL RLINGAMA(QX,ALPHA+1.D0,FX)
        ENDIF
      ENDIF
      IF (DABS(X-U).LT.1.D-6) THEN
        F0=FX
        F2=FX
        X=EL
        CALL RLQGAMMA(X,ALPHA,QX)
        GOTO 100
      ELSE
        V(3)=QX
        F1=FX
      ENDIF
      TMP=(F2-F1)/(U-EL)-1.D0
      RLLEQNG=DBLE(TMP)
      RETURN
      END
C
      SUBROUTINE RLQUQLDG(U,ALPHA,SIGMA,TOL,QL,QU,ISOL)
      implicit double precision(a-h,o-z)
c      DOUBLE PRECISION RLLEQNG,QQL,SMIN,SMAX,TOLD
      DIMENSION Q(3),PARAM(1)
      EXTERNAL RLLEQNG
c      COMMON/DATI/NQ
      QL=0.D0
      ISOL=0
      NQ=3
      Q(1)=U
      Q(2)=ALPHA
      CALL RLQGAMMA(U,ALPHA,QU)
      Q(3)=QU
      PARAM(1)=QU
      QU=QU*SIGMA
      SMIN=0.D0
      SMAX=1.D0
      TOLD=DBLE(TOL)
      CALL RLRGFLD(RLLEQNG,Q,0.D0,SMIN,SMAX,TOLD,100,QQL,ITERM,NQ,PARAM)
      IF (ITERM.NE.1) RETURN
      QL=Q(3)*SIGMA
      ISOL=1
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLGAMMAD(SIGMA,ALFA,X)
C.......................................................................
C
C   COPYRIGHT 1993 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      implicit double precision(a-h,o-z)
c      DOUBLE PRECISION X,XLGMN,YLGMN,GALIM,S,ALOGS,ALOGAM,ALF,F,V,Z,GL
      DATA NCALL,XLGMN,YLGMN,GALIM/0,0.D0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,GALIM)
        CALL RLMACHD(4,XLGMN)
        CALL RLMACHD(5,YLGMN)
c        GALIM=DLOG(1.D-10)
      ENDIF
c      IF (SIGMA.LE.0.D0.OR.ALFA.LE.0.D0.OR.X.LT.0.D0)
c     1 CALL MESSGE(500,'RLGAMMAD',1)
      RLGAMMAD=0.D0
      IF (X.EQ.0.D0) RETURN
      S=X/DBLE(SIGMA)
      ALOGS=YLGMN
      IF (S.GT.XLGMN) ALOGS=DLOG(S)
      ALF=DBLE(ALFA)
      V=ALF
      F=0.D0
      IF (ALF.GE.7.D0) GOTO 300
      F=1.D0
      Z=ALF-1.D0
  100 Z=Z+1.D0
      IF (Z.GE.7.D0) GOTO 200
      V=Z
      F=F*Z
      GOTO 100
  200 V=V+1.D0
      F=-DLOG(F)
  300 Z=1.D0/(V*V)
      GL=F+(V-0.5D0)*DLOG(V)-V+.9189385332D0+(((-.000595238D0*Z+
     +   .0007936507D0)*Z - .0027777778D0)*Z+.0833333333D0)/V
      ALOGAM=(ALF-1.D0)*ALOGS-S-DLOG(DBLE(SIGMA))-GL
      IF (ALOGAM.LE.GALIM) RETURN
      RLGAMMAD=DEXP(ALOGAM)
      RETURN
      END
C
      SUBROUTINE RLMEDMAD(X,NT,TETA,TMEANF,TMADF)
c      implicit double precision(a-h,o-z)
      DOUBLE PRECISION BETA,GAM,mF,fm,uF,lF,WF,TRMNF,TRMDF,GX,
     + FQGLOW,FQGUP,D1,D1F,D2,D2F,QGUP2F,quF,fquF,A1F,B1,B1F,QGLOW2F
      DOUBLE PRECISION A2F,B2,B2F,QGUP1F,A3F,B3,B3F,QGLOW1F,A4F,B4,
     + B4F,FA,FB,Y,X,TMEANF,TMADF,TETA(NT)
c      COMMON/TETAPR/TETA(60)
      BETA=DBLE(TETA(4))
      GAM=DBLE(TETA(5))
      Y=DBLE(X)
c IF.trmean.f
      mF=DBLE(TETA(6))
      fm=DBLE(TETA(27))
      TRMNF=0.D0
      uF=DBLE(TETA(12))
      lF=DBLE(TETA(13))
      WF=DBLE(TETA(14))
      IF (Y.LT.lF) TRMNF=TRMNF+lF-WF
      IF (lF.LE.Y.AND.Y.LE.uF) TRMNF=TRMNF+Y-WF
      IF (uF.LT.Y) TRMNF=TRMNF+(uF-WF)
      TRMNF=TRMNF/(1.D0-2.D0*BETA)
      TMEANF=TRMNF
c IF.D1.f
      FQGLOW=DBLE(TETA(29))
      FQGUP=DBLE(TETA(28))
      D1=DBLE(TETA(8))
      GX=-GAM
      IF (mF+D1.GE.Y) GX=GX+1.D0
      IF (mF-D1.GE.Y) GX=GX-1.D0
      D1F=(-GX+TRMNF*(FQGLOW-FQGUP))/(FQGLOW+FQGUP)
c IF.D2.f
      FQGLOW=DBLE(TETA(31))
      FQGUP=DBLE(TETA(30))
      D2=DBLE(TETA(9))
      GX=-(1.D0-GAM)
      IF (mF+D2.GE.Y) GX=GX+1.D0
      IF (mF-D2.GE.Y) GX=GX-1.D0
      D2F=(-GX+TRMNF*(FQGLOW-FQGUP))/(FQGLOW+FQGUP)
c IF.QGup2.f
      QGUP2F=TRMNF+D2F
c IF.A1.f
      quF=DBLE(TETA(25))
      fquF=DBLE(TETA(30))
      A1F=-DBLE(TETA(15))
      IF (quF.GE.Y) A1F=A1F+Y
      A1F=A1F+fquF*quF*QGUP2F
c IF.B1.f
      B1=DBLE(TETA(19))
      B1F=-B1
      IF (quF.GE.Y) B1F=B1F+1.D0
      B1F=B1F+fquF*QGUP2F
c IF.QGlow2.f
      QGLOW2F=TRMNF-D2F
c IF.A2.f
      quF=DBLE(TETA(26))
      fquF=DBLE(TETA(31))
      A2F=-DBLE(TETA(16))
      IF (Y.LE.quf) A2F=A2F+Y
      A2F=A2F+fquF*quF*QGLOW2F
c IF.B2.f
      B2=DBLE(TETA(20))
      B2F=-B2
      IF (Y.LE.quf) B2F=B2F+1.D0
      B2F=B2F+fquF*QGLOW2F
c IF.QGup1.f
      QGUP1F=TRMNF+D1F
c IF.A3.f
      quF=DBLE(TETA(23))
      fquF=DBLE(TETA(28))
      A3F=-DBLE(TETA(17))
      IF (quF.GE.Y) A3F=A3F+Y
      A3F=A3F+fquF*quF*QGUP1F
c IF.B3.f
      B3=DBLE(TETA(21))
      B3F=-B3
      IF (quF.GE.Y) B3F=B3F+1.D0
      B3F=B3F+fquF*QGUP1F
c IF.QGlow1.f
      QGLOW1F=TRMNF-D1F
c IF.A4.f
      quF=DBLE(TETA(24))
      fquF=DBLE(TETA(29))
      A4F=-DBLE(TETA(18))
      IF (Y.LE.quF) A4F=A4F+Y
      A4F=A4F+fquF*quF*QGLOW1F
c IF.B4.f
      B4=DBLE(TETA(22))
      B4F=-B4
      IF (Y.LE.quF) B4F=B4F+1.D0
      B4F=B4F+fquF*QGLOW1F
c IF.trmadv.f
      FA=A1F+A2F-A3F-A4F
      FB=B1F+B2F-B3F-B4F
      TRMDF=(FA-FB*mF-TRMNF*(B1+B2-B3-B4))/(1.D0-2.D0*GAM)
      TMADF=TRMDF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RLIFGANS(DX,WGT,N,EXGAM,EXPSI,NT,TETA)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION ANS,EXGAM,EXPSI,DX,IFG(1),WGT(N),TETA(NT),XX(1)
      EXTERNAL EXGAM,EXPSI
C     COMMON/TETAPR/TETA(60)
C
C avoid compiler warnings
C
      dummy = expsi(1.D0,1,1,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0)
C
      RLIFGANS=0.D0
      ANS=EXGAM(WGT(2),WGT(1),DX)
      IF (ANS.LE.1.D-15) RETURN
      ALPHA=WGT(1)
      SIGMA=WGT(2)
      ITC=INT(SNGL(WGT(3)))
      XX(1)=DX
      CALL RLIFGAMA(XX,TETA,1,NT,ALPHA,SIGMA,ITC,0,IFG)
      IF (ITC.GE.0) IFG(1)=IFG(1)*IFG(1)
      RLIFGANS=IFG(1)*ANS
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RLIFGAMA(DX,TETA,NX,NT,ALPHA,SIGMA,ITC,ITT,IFG)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION  mF,MpF,madF,medF,J1,J0G,J1G,IFG(NX),X
      DIMENSION TETA(NT),DX(NX)
C     COMMON/TETAPR/TETA(60)
C
c      IF (ITT.EQ.1) THEN
c        DO 100 I=1,NT
c        TETA(I)=THETA(I)
c  100   CONTINUE
c      ENDIF
      DUMMY=ITT
      DO 500 IX=1,NX
      X=DX(IX)
      IF (TETA(1).EQ.2) GOTO 200
      CALL RLMEDMAD(X,NT,TETA,TRMNF,TRMDF)
      IF (ITC.NE.1.AND.ITC.NE.2) GOTO 103
      IF (ITC.EQ.1) XIF=TRMNF
      IF (ITC.EQ.2) XIF=TRMDF
      GOTO 300
c IF.alpha.Dsm.g
  103 mF=TETA(6)
      sF=TETA(10)
      MpF=TETA(32)
      SpF=TETA(33)
      A=(TRMNF-mF*TRMDF/sF)/sF
      B=(MpF-SpF*mF/sF)*SIGMA/sF
      ALDSMG=A/B
      IF (ITC.NE.3) GOTO 104
      XIF=ALDSMG
      GOTO 300
c IF.sigma.Dsm.g
  104 SIDSMG=(TRMDF-SpF*ALDSMG*SIGMA)*SIGMA/sF
      IF (IABS(ITC).NE.4) GOTO 105
      XIF=SIDSMG
      IF (ITC.EQ.-4) XIF=XIF*ALDSMG
      GOTO 300
c IF.mu.g
  105 XIF=SIGMA*ALDSMG+ALPHA*SIDSMG
      IF (ITC.EQ.5) GOTO 300
c IF.qu1.g
      QU1G=-ALDSMG*TETA(51)/TETA(50)
      IF (ITC.NE.6) GOTO 107
      XIF=QU1G
      GOTO 300
c IF.qu.Dsm.g
  107 QUDSMG=SIDSMG*TETA(37)+QU1G*SIGMA
      IF (ITC.NE.7) GOTO 108
      XIF=QUDSMG
      GOTO 300
c IF.H0.g
  108 quF=TETA(38)
      fquF=TETA(41)
      H0=TETA(44)
      H0G=-H0
      IF (quF.GE.X) H0G=H0G+1.
      H0G=H0G+fquF*QUDSMG
      IF (ITC.NE.8) GOTO 109
      XIF=H0G
      GOTO 300
cc IF.H1.g
  109 H1=TETA(45)
      H1G=-H1
      IF (quF.GE.X) H1G=H1G+X
      H1G=H1G+fquF*quF*QUDSMG
      IF (ITC.NE.9) GOTO 110
      XIF=H1G
      GOTO 300
c IF.ql1.g
  110 u=TETA(34)
      B=(u-TETA(48))*ALDSMG+(TETA(56)-TETA(55)-TETA(52)*TETA(35))*ALDSMG
     +   -TETA(53)*QU1G
      A=TETA(49)*TETA(35)-TETA(54)
      QL1G=B/A
      IF (ITC.NE.10) GOTO 111
      XIF=QL1G
      GOTO 300
c IF.ql.Dsm.g
  111 QLDSMG=SIDSMG*TETA(39)+QL1G*SIGMA
      IF (ITC.NE.11) GOTO 112
      XIF=QLDSMG
      GOTO 300
c IF.J0.g
  112 qlF=TETA(40)
      fqlF=TETA(42)
      J0G=-TETA(46)
      IF (qlF.GE.X) J0G=J0G+1.
      J0G=J0G+fqlF*QLDSMG
      IF (ITC.NE.12) GOTO 113
      XIF=J0G
      GOTO 300
c IF.J1.g
  113 J1=TETA(47)
      J1G=-J1
      IF (qlF.GE.X) J1G=J1G+X
      J1G=J1G+fqlF*qlF*QLDSMG
      IF (ITC.NE.13) GOTO 114
      XIF=J1G
      GOTO 300
c IF.tcmean.g
  114 H0mJ0=TETA(44)-TETA(46)
      XIF=(-(H1-J1)*(H0G-J0G)/H0mJ0 + H1G-J1G)/H0mJ0
      GOTO 300
c IF.med.f
  200 medF=TETA(31)
      emdF=0.5
      IF (medF.GE.X) emdF=-0.5
      emdF=emdF/TETA(33)
      IF (ITC.NE.1) GOTO 202
      XIF=emdF
      GOTO 300
c IF.mad.f
  202 madF=TETA(32)
      fMmD=TETA(35)
      fMpD=TETA(34)
      amdF=emdF*(fMmD-fMpD)
      SGN=DABS(X-medF)-madF
      IF (DABS(SGN).GT.1D-6) amdF=amdF+DSIGN(0.5D0,SGN)
      amdF=amdF/(fMmD+fMpD)
      IF (ITC.NE.2) GOTO 203
      XIF=amdF
      GOTO 300
c IF.alpha.D.g
  203 A=(emdF-medF*amdF/madF)/madF
      D1F=TETA(28)
      DpF=TETA(30)
      B=(TETA(29)-DpF*TETA(27)/D1F)/D1F
      ALFDG=A/B
      IF (ITC.NE.3) GOTO 204
      XIF=ALFDG
      GOTO 300
c IF.sigma.D.g
  204 SIGDG=(amdF-madF*DpF*ALFDG/D1F)/D1F
      IF (IABS(ITC).NE.4) GOTO 205
      XIF=SIGDG
      IF (ITC.EQ.-4) XIF=XIF*ALFDG
      GOTO 300
c IF.mu.g
  205 XIF=SIGMA*ALFDG+ALPHA*SIGDG
      IF (ITC.EQ.5) GOTO 300
c IF.qu1.g
      QU1G=-ALFDG*TETA(19)/TETA(18)
      IF (ITC.NE.6) GOTO 207
      XIF=QU1G
      GOTO 300
c IF.qu.D.g
  207 QUDG=SIGDG*TETA(5)+QU1G*SIGMA
      IF (ITC.NE.7) GOTO 208
      XIF=QUDG
      GOTO 300
c IF.H0.g
  208 quF=TETA(6)
      fquF=TETA(9)
      u=TETA(2)
      H0=TETA(12)
      H0G=-H0
      IF (quF.GE.X) H0G=H0G+1.
      H0G=H0G+fquF*QUDG
      IF (ITC.NE.8) GOTO 209
      XIF=H0G
      GOTO 300
c IF.H1.g
  209 H1=TETA(13)
      H1G=-H1
      IF (quF.GE.X) H1G=H1G+X
      H1G=H1G+fquF*quF*QUDG
      IF (ITC.NE.9) GOTO 210
      XIF=H1G
      GOTO 300
c IF.ql1.g
  210 B=(u-TETA(16))*ALFDG+(TETA(24)-TETA(23)-TETA(20)*TETA(3))*ALFDG
     +   - TETA(21)*QU1G
      A=TETA(17)*TETA(3)-TETA(22)
      QL1G=B/A
      IF (ITC.NE.10) GOTO 211
      XIF=QL1G
      GOTO 300
c IF.ql.D.g
  211 QLDG=SIGDG*TETA(7)+QL1G*SIGMA
      IF (ITC.NE.11) GOTO 212
      XIF=QLDG
      GOTO 300
c IF.J0.g
  212 qlF=TETA(8)
      fqlF=TETA(10)
      J0G=-TETA(14)
      IF (qlF.GE.X) J0G=J0G+1.
      J0G=J0G+fqlF*QLDG
      IF (ITC.NE.12) GOTO 213
      XIF=J0G
      GOTO 300
c IF.J1.g
  213 J1=TETA(15)
      J1G=-J1
      IF (qlF.GE.X) J1G=J1G+X
      J1G=J1G+fqlF*qlF*QLDG
      IF (ITC.NE.13) GOTO 214
      XIF=J1G
      GOTO 300
c IF.tcmean.g
  214 H0mJ0=TETA(12)-TETA(14)
      XIF=(-(H1-J1)*(H0G-J0G)/H0mJ0 + H1G-J1G)/H0mJ0
  300 IFG(IX)=DBLE(XIF)
  500 CONTINUE
      RETURN
      END
C
      SUBROUTINE RLAVTCMG(TETA,NT,ALPHA,SIGMA,ITC,UPPER,TIL,SUM,
     +           IWORK,WORK)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION RLIFGANS,TILD,ERRSTD,WORK,RLGAMMAD,LO,HI,SUM
      DIMENSION TETA(NT),WGT(3),IWORK(80),WORK(320)
c      COMMON/INTEGN/AINTEG(4),IWORK(80),WORK(320)
c      COMMON/TETAPR/TETA(60)
      EXTERNAL RLIFGANS,RLGAMMAD,rlexu
C
C     INITIALISATION
C
      LIMIT=80
      KEY=1
C
      WGT(1)=ALPHA
      WGT(2)=SIGMA
      WGT(3)=DBLE(FLOAT(ITC))
      LO=0.D0
      HI=DBLE(UPPER)
c      DO 100 I=1,NT
c      TETA(I)=THETA(I)
c  100 CONTINUE
      SUM=0.D0
      TILD=DBLE(TIL)
      CALL RLINTGRT(RLIFGANS,WGT,3,RLGAMMAD,rlexu,LO,HI,TILD,TILD,
     1            KEY,LIMIT,SUM,ERRSTD,NEVAL,IER,WORK,IWORK,NT,TETA)
      RETURN
      END
C
C=======================================================================
      SUBROUTINE RLTRMNLW(ALPHA,SIGMA,BETA,mF)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION GL,GU,mF
      IF (DABS(BETA-0.5D0).LT.1.D-5) THEN
        mF=-0.3665129D0
        GOTO 100
      ENDIF
      IF (BETA.LT.1.D-5) THEN
        mF=-0.5772157D0
      ELSE
        CALL RLQWEIBL(1.D0,1.D0,BETA,QL)
        CALL RLQWEIBL(1.D0,1.D0,1.D0-BETA,QU)
        CALL RLSUMLGM(DBLE(QL),1.D0,GL)
        CALL RLSUMLGM(DBLE(QU),1.D0,GU)
        mF=DBLE(GU-GL)/(1.D0-2.D0*BETA)
      ENDIF
  100 mF=mF/ALPHA + DLOG(SIGMA)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLEQADLW(DD,V,NV,PARAM)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION DD,V(NV),M,PARAM(1)
      BETA =V(1)
      GAM  =V(2)
      DUMMY=PARAM(1)
      D=DBLE(DD)
      CALL RLTRMNLW(1.D0,1.D0,BETA,M)
      CALL RLPWEIBL(1.D0,1.D0,DEXP(M+D),P1)
      CALL RLPWEIBL(1.D0,1.D0,DEXP(M-D),P2)
      TMP=P1-P2-(1.D0-GAM)
      RLEQADLW=DBLE(TMP)
      RETURN
      END
C
      SUBROUTINE RLQAD1LW(BETA,GAM,TOL,QAD1,ISOL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION RLEQADLW,QAD,SMIN,SMAX,TOLD,PARAM(1),Q(2)
      EXTERNAL RLEQADLW
C     COMMON/DATI/NQ
      QAD1=0.D0
      ISOL=0
      NQ=2
      Q(1)=BETA
      Q(2)=GAM
      PARAM(1)=0.D0
      SMIN=0.D0
      SMAX=10.D0
      TOLD=DBLE(TOL)
      CALL RLRGFLD(RLEQADLW,Q,0.D0,SMIN,SMAX,TOLD,100,QAD,ITRM,NQ,PARAM)
      IF (ITRM.NE.1) RETURN
      QAD1=DBLE(QAD)
      ISOL=1
      RETURN
      END
C
      SUBROUTINE RLTRMADLW(ALPHA,BETA,GAM,TOL,mF,sF,ISOL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION AA0,AA1,AA2,mF
      ISOL=1
      CALL RLTRMNLW(1.D0,1.D0,BETA,mF)
      IF (DABS(GAM-0.5D0).LT.1.D-5) THEN
        CALL RLQAD1LW(BETA,GAM,TOL,sF,ISOL)
      ELSEIF(DABS(GAM-0.4D0).LT.1.D-5.AND.DABS(BETA-0.4D0).LT.1.D-5)THEN
        sF=0.7707968D0
      ELSE
        GAM0=GAM
        EmF=DEXP(mF)
        CALL RLSUMLGM(DBLE(EmF),1.D0,AA0)
        A0=DBLE(AA0)
        CALL RLPWEIBL(1.D0,1.D0,EmF,B0)
 100    CALL RLQAD1LW(BETA,GAM0,TOL,QADF,JSOL)
        EmQ=DEXP(mF+QADF)
        CALL RLSUMLGM(DBLE(EmQ),1.D0,AA1)
        A1=DBLE(AA1)
        CALL RLPWEIBL(1.D0,1.D0,EmQ,B1)
        EmQ=DEXP(mF-QADF)
        CALL RLSUMLGM(DBLE(EmQ),1.D0,AA2)
        A2=DBLE(AA2)
        CALL RLPWEIBL(1.D0,1.D0,EmQ,B2)
        IF (DABS(GAM-GAM0).LT.1.D-6) THEN
          SUP=(A1 + A2 - 2.D0 * A0) - mF * (B1 + B2 - 2.D0 * B0)
          GAM0=1.-GAM
          ISOL=JSOL
          GOTO 100
        ELSE
          SLOW=(A1 + A2 - 2.D0 * A0) - mF * (B1 + B2 - 2.D0 * B0)
          IF (JSOL.EQ.0) ISOL=0
        ENDIF
        sF=(SUP-SLOW)/(1.D0-2.D0*GAM)
      ENDIF
      sF=sF/ALPHA
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION RLLEQNW(LL,V,NV,PARAM)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION LL,LX,V(NV),PARAM(1)
      U=V(1)
      ALPHA=V(2)
      EL=DBLE(LL)
      DUMMY=PARAM(1)
      IF (DABS(EL-U).LT.0.0001D0) EL=EL-0.0002D0
      X=U
  100 IF (X.LT.1D-5) THEN
        FX=0.D0
      ELSE
        IF (1.D0-X.LT.1D-5) THEN
          FX=1.D0
        ELSE
          LX=DLOG(1.D0/(1.D0-X))
          CALL RLINGAMA(LX,1.D0/ALPHA+1.D0,FX)
        ENDIF
      ENDIF
      IF (DABS(X-U).LT.1D-6) THEN
        F2=FX
        X=EL
        GOTO 100
      ELSE
        F1=FX
      ENDIF
      TMP=(F2-F1)/(U-EL)-1.D0
      RLLEQNW=DBLE(TMP)
      RETURN
      END
C
      SUBROUTINE RLQUQLDW(U,ALPHA,SIGMA,TOL,QL,QU,ISOL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION RLLEQNW,QQL,SMIN,SMAX,TOLD,PARAM(1),Q(2)
      EXTERNAL RLLEQNW
c     COMMON/DATI/NQ
      QL=0.D0
      ISOL=0
      NQ=2
      CALL RLQWEIBL(ALPHA,SIGMA,U,QU)
      Q(1)=U
      Q(2)=ALPHA
      PARAM(1)=0.D0
      SMIN=0.D0
      SMAX=0.5D0
      TOLD=DBLE(TOL)
      CALL RLRGFLD(RLLEQNW,Q,0.D0,SMIN,SMAX,TOLD,100,QQL,ITERM,NQ,PARAM)
      IF (ITERM.NE.1) RETURN
      QQ=DBLE(QQL)
      CALL RLQWEIBL(ALPHA,SIGMA,QQ,QL)
      ISOL=1
      RETURN
      END
C
      SUBROUTINE RLTRMNW(ALPHA,SIGMA,BETA,mF)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION mF,LPL,LPU
      EXTERNAL RLXEXPD
      IF (DABS(BETA-0.5D0).LT.1.D-5) THEN
        CALL RLQWEIBL(ALPHA,SIGMA,0.5D0,mF)
        RETURN
      ENDIF
      ALF1=1.D0+1.D0/ALPHA
      CALL RLLGAMAD(ALF1,GA)
      IF (BETA.LT.1.D-5) THEN
        mF=SIGMA*RLXEXPD(GA)
      ELSE
c       CALL RLQWEIBL(ALPHA,1.D0,BETA,QL)
c       CALL RLQWEIBL(ALPHA,1.D0,1.D0-BETA,QU)
c       CALL RLPWEIBL(ALPHA,1.D0,QL,PL)
c       CALL RLPWEIBL(ALPHA,1.D0,QU,PU)
        PU=1.D0-BETA
        GU=PU*RLXEXPD(GA)
        IF (DABS(PU-1.D0).LT.1D-5) GOTO 20
        LPU=DLOG(1.D0/(1.D0-PU))
        CALL RLINGAMA(LPU,ALF1,GU)
        GU=GU*RLXEXPD(GA)
   20   PL=BETA
        GL=PL*RLXEXPD(GA)
        IF (DABS(PL-1.D0).LT.1D-5) GOTO 40
        LPL=DLOG(1.D0/(1.D0-PL))
        CALL RLINGAMA(LPL,ALF1,GL)
        GL=GL*RLXEXPD(GA)
   40   mF=SIGMA*(GU-GL)/(1.D0-2.D0*BETA)
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION RLEQADW(DD,V)
      DOUBLE PRECISION DD,V(3),M,D,ALPHA,BETA,GAM,P1,P2
      ALPHA=V(1)
      BETA =V(2)
      GAM  =V(3)
      D=DBLE(DD)
      CALL RLTRMNW(ALPHA,1.D0,BETA,M)
      CALL RLPWEIBL(ALPHA,1.D0,M+D,P1)
      CALL RLPWEIBL(ALPHA,1.D0,M-D,P2)
      RLEQADW= P1- P2 - (1.D0-GAM)
      RETURN
      END
C
      SUBROUTINE RLQAD1W(ALPHA,BETA,GAM,TOL,QAD1,ISOL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION RLEQADW,QAD,SMIN,SMAX,TOLD,PARAM(1),Q(3)
      EXTERNAL RLEQADW
C     COMMON/DATI/NQ
      QAD1=0.D0
      ISOL=0
      NQ=3
      Q(1)=ALPHA
      Q(2)=BETA
      Q(3)=GAM
      PARAM(1)=0.D0
      SMIN=0.D0
      CALL RLQWEIBL(ALPHA,1.D0,0.999D0,Q9)
      SMAX=DBLE(Q9)
      TOLD=DBLE(TOL)
      CALL RLRGFLD(RLEQADW,Q,0.D0,SMIN,SMAX,TOLD,100,QAD,ITERM,NQ,PARAM)
      IF (ITERM.NE.1) RETURN
      QAD1=DBLE(QAD)
      ISOL=1
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION FUNCTION RLWEIBUD(SIGMA,ALPHA,X)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION X,XLGMN,YLGMN,EXMIN,T,TA,TMP,ALF,SS
      DATA NCALL,EXMIN,XLGMN,YLGMN/0,0.D0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
        CALL RLMACHD(4,XLGMN)
        CALL RLMACHD(5,YLGMN)
      ENDIF
c      IF (SIGMA.LE.0..OR.ALPHA.LE.0.0.OR.X.LT.0.D0)
c     + CALL MESSGE(500,'RLWEIBUD',1)
      RLWEIBUD=0.D0
      IF (X.LE.0.D0) RETURN
      SS=DBLE(SIGMA)
      T=X/SS
      TMP=YLGMN
      IF (T.GT.XLGMN) TMP=DLOG(T)
      TA=0.D0
      ALF=DBLE(ALPHA)
      IF (ALF*TMP.GT.EXMIN) TA=DEXP(ALF*TMP)
      TMP=DLOG(ALF)-DLOG(SS)+(ALF-1.D0)*TMP-TA
      IF (TMP.LE.EXMIN) RETURN
      RLWEIBUD=DEXP(TMP)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLIFWANS(DX,WGT,N,EXWEB,EXPSI,NT,TETA)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION ANS,EXWEB,DX(1),IFW(1),WGT(N),TETA(NT)
      EXTERNAL EXWEB,EXPSI
C     COMMON/TETAPR/TETA(60)
C
C avoid compiler warnings
C
      dummy = expsi(1.D0,1,1,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0)
C
C
C  Initializations
C
      RLIFWANS=0.D0
      ANS=EXWEB(WGT(2),WGT(1),DX)
      IF (ANS.LE.1.D-15) RETURN
      ALPHA=WGT(1)
      SIGMA=WGT(2)
      ITC=INT(SNGL(WGT(3)))
      CALL RLIFWEIB(DX,TETA,1,NT,ALPHA,SIGMA,ITC,0,IFW)
      IF (ITC.GE.0) IFW(1)=IFW(1)*IFW(1)
      RLIFWANS=IFW(1)*ANS
      RETURN
      END
C
      SUBROUTINE RLIFWEIB(DX,TETA,NX,NT,ALPHA,SIGMA,ITC,ITT,IFW)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION TETA(NT),madF,medF,MU1W,J1,J0W,J1W,DX(NX),IFW(NX)
      EXTERNAL RLXEXPD,RLGAMDIGAMA
c      COMMON/TETAPR/TETA(60)
      DATA NCALL,XLGMN,YLGMN,ALF0,GA,DGA/0,5*0.D0/
C
C  Initializations
C
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(4,XLGMN)
        CALL RLMACHD(5,YLGMN)
      ENDIF
      IF (DABS(ALPHA-ALF0).GT.1D-6) THEN
        TA=1.D0+1.D0/ALPHA
        ALF0=ALPHA
        CALL RLLGAMAD(TA,GA)
        DGA=RLGAMDIGAMA(TA)
      ENDIF
c      IF (ITT.EQ.1) THEN
c        DO 100 I=1,NT
c        TETA(I)=THETA(I)
c  100   CONTINUE
c      ENDIF
      DUMMY=ITT
      DO 500 IX=1,NX
      X=DBLE(DX(IX))
      IF (TETA(1).EQ.2) GOTO 200
      Y=YLGMN
      IF (X.GT.XLGMN) Y=DLOG(X)
      CALL RLMEDMAD(Y,NT,TETA,TRMNF,TRMDF)
      IF (ITC.EQ.1) XIF=TRMNF
      IF (ITC.EQ.2) XIF=TRMDF
      IF (ITC.EQ.1.OR.ITC.EQ.2) GOTO 300
c IF.alpha.Dsm.w
      s1F=TETA(11)
      ALDSMW=-TRMDF*ALPHA**2/s1F
      IF (ITC.NE.3) GOTO 104
      XIF=ALDSMW
      GOTO 300
c IF.tau.Dsm.w
  104 TAUDSMW=TRMNF-TRMDF*TETA(7)/s1F
c IF.sigma.Dsm.w
      SIDSMW=TAUDSMW*SIGMA
      IF (IABS(ITC).NE.4) GOTO 105
      XIF=SIDSMW
      IF (ITC.EQ.-4) XIF=XIF*ALDSMW
      GOTO 300
c IF.mu.w
  105 IFW(IX)=DEXP(DBLE(GA))*DBLE(SIDSMW-SIGMA*ALDSMW*DGA/ALPHA**2)
      IF (ITC.EQ.5) GOTO 500
c IF.qu1.w
      QU1W=-ALDSMW*TETA(51)/TETA(50)
      IF (ITC.NE.6) GOTO 107
      XIF=QU1W
      GOTO 300
c IF.qu.Dsm.w
  107 QUDSMW=SIDSMW*TETA(37)+QU1W*SIGMA
      IF (ITC.NE.7) GOTO 108
      XIF=QUDSMW
      GOTO 300
c IF.H0.w
  108 quF=TETA(38)
      fquF=TETA(41)
      H0=TETA(44)
      u=TETA(34)
      H0W=-H0
      IF (quF.GE.X) H0W=H0W+1.D0
      H0W=H0W+fquF*QUDSMW
      IF (ITC.NE.8) GOTO 109
      XIF=H0W
      GOTO 300
c IF.H1.w
  109 H1=TETA(45)
      H1W=-H1
      IF (quF.GE.X) H1W=H1W+X
      H1W=H1W+fquF*quF*QUDSMW
      IF (ITC.NE.9) GOTO 110
      XIF=H1W
      GOTO 300
c IF.mu1.w
  110 MU1W=-DGA*ALDSMW*RLXEXPD(GA)/ALPHA**2
c IF.ql1.w
      B=(u-TETA(48))*MU1W+(TETA(56)-TETA(55)-TETA(52)*TETA(35))*ALDSMW-
     +   TETA(53)*QU1W
      A=TETA(49)*TETA(35)-TETA(54)
      QL1W=B/A
      IF (ITC.NE.10) GOTO 111
      XIF=QL1W
      GOTO 300
c IF.ql.Dsm.w
  111 QLDSMW=SIDSMW*TETA(39)+QL1W*SIGMA
      IF (ITC.NE.11) GOTO 112
      XIF=QLDSMW
      GOTO 300
c IF.J0.w
  112 qlF=TETA(40)
      fqlF=TETA(42)
      J0W=-TETA(46)
      IF (qlF.GE.X) J0W=J0W+1.D0
      J0W=J0W+fqlF*QLDSMW
      IF (ITC.NE.12) GOTO 113
      XIF=J0W
      GOTO 300
c IF.J1.w
  113 J1=TETA(47)
      J1W=-J1
      IF (qlF.GE.X) J1W=J1W+X
      J1W=J1W+fqlF*qlF*QLDSMW
      IF (ITC.NE.13) GOTO 114
      XIF=J1W
      GOTO 300
c IF.tcmean.w
  114 H0mJ0=TETA(44)-TETA(46)
      XIF=(-(H1-J1)*(H0W-J0W)/H0mJ0 + H1W-J1W)/H0mJ0
      GOTO 300
c IF.med.f
  200 medF=TETA(31)
      emdF=0.5D0
      IF (medF.GE.X) emdF=-0.5D0
      emdF=emdF/TETA(33)
      IF (ITC.NE.1) GOTO 202
      XIF=emdF
      GOTO 300
c IF.mad.f
  202 madF=TETA(32)
      fMmD=TETA(35)
      fMpD=TETA(34)
      amdF=emdF*(fMmD-fMpD)
      SGN=DABS(X-medF)-madF
      IF (DABS(SGN).GT.1D-6) amdF=amdF+DSIGN(0.5D0,SGN)
      amdF=amdF/(fMmD+fMpD)
      IF (ITC.NE.2) GOTO 203
      XIF=amdF
      GOTO 300
c IF.alpha.D.w
  203 A=(emdF-medF*amdF/madF)/madF
      D1F=TETA(28)
      DpF=TETA(30)
      B=(TETA(29)-DpF*TETA(27)/D1F)/D1F
      ALFDW=A/B
      IF (ITC.NE.3) GOTO 204
      XIF=ALFDW
      GOTO 300
c IF.sigma.D.w
  204 SIGDW=(amdF-madF*DpF*ALFDW/D1F)/D1F
      IF (IABS(ITC).NE.4) GOTO 205
      XIF=SIGDW
      IF (ITC.EQ.-4) XIF=XIF*ALFDW
      GOTO 300
c IF.mu.w
  205 IFW(IX)=DEXP(DBLE(GA))*DBLE(SIGDW-SIGMA*ALFDW*DGA/ALPHA**2)
      IF (ITC.EQ.5) GOTO 500
c IF.qu1.w
      QU1W=-ALFDW*TETA(19)/TETA(18)
      IF (ITC.NE.6) GOTO 207
      XIF=QU1W
      GOTO 300
c IF.qu.D.w
  207 QUDW=SIGDW*TETA(5)+QU1W*SIGMA
      IF (ITC.NE.7) GOTO 208
      XIF=QUDW
      GOTO 300
c IF.H0.w
  208 quF=TETA(6)
      fquF=TETA(9)
      H0=TETA(12)
      u=TETA(2)
      H0W=-H0
      IF (quF.GE.X) H0W=H0W+1.D0
      H0W=H0W+fquF*QUDW
      IF (ITC.NE.8) GOTO 209
      XIF=H0W
      GOTO 300
c IF.H1.w
  209 H1=TETA(13)
      H1W=-H1
      IF (quF.GE.X) H1W=H1W+X
      H1W=H1W+fquF*quF*QUDW
      IF (ITC.NE.9) GOTO 210
      XIF=H1W
      GOTO 300
c IF.mu1.w
  210 MU1W=-DGA*ALFDW*RLXEXPD(GA)/ALPHA**2
c IF.ql1.w
      B=(u-TETA(16))*MU1W+(TETA(24)-TETA(23)-TETA(20)*TETA(3))*ALFDW
     +   - TETA(21)*QU1W
      A=TETA(17)*TETA(3)-TETA(22)
      QL1W=B/A
      IF (ITC.NE.10) GOTO 211
      XIF=QL1W
      GOTO 300
c IF.ql.D.w
  211 QLDW=SIGDW*TETA(7)+QL1W*SIGMA
      IF (ITC.NE.11) GOTO 212
      XIF=QLDW
      GOTO 300
c IF.J0.w
  212 qlF=TETA(8)
      fqlF=TETA(10)
      J0W=-TETA(14)
      IF (qlF.GE.X) J0W=J0W+1.D0
      J0W=J0W+fqlF*QLDW
      IF (ITC.NE.12) GOTO 213
      XIF=J0W
      GOTO 300
c IF.J1.w
  213 J1=TETA(15)
      J1W=-J1
      IF (qlF.GE.X) J1W=J1W+X
      J1W=J1W+fqlF*qlF*QLDW
      IF (ITC.NE.13) GOTO 214
      XIF=J1W
      GOTO 300
c IF.tcmean.w
  214 H0mJ0=TETA(12)-TETA(14)
      XIF=(-(H1-J1)*(H0W-J0W)/H0mJ0 + H1W-J1W)/H0mJ0
  300 IFW(IX)=DBLE(XIF)
  500 CONTINUE
      RETURN
      END
C
      SUBROUTINE RLAVTCMW(TETA,NT,ALPHA,SIGMA,ITC,UPPER,TIL,SUM,
     +           IWORK,WORK)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION RLIFWANS,TILD,ERRSTD,WORK,RLWEIBUD,LO,HI,SUM
      DIMENSION TETA(NT),WGT(3),IWORK(80),WORK(320)
C     COMMON/TETAPR/TETA(60)
      EXTERNAL RLIFWANS,RLWEIBUD,rlexu
C
C     INITIALISATION
C
      LIMIT=80
      KEY=1
C
      WGT(1)=ALPHA
      WGT(2)=SIGMA
      WGT(3)=DBLE(FLOAT(ITC))
      LO=0.D0
      HI=DBLE(UPPER)
c      DO 100 I=1,NT
c      TETA(I)=THETA(I)
c  100 CONTINUE
      SUM=0.D0
      TILD=DBLE(TIL)
      CALL RLINTGRT(RLIFWANS,WGT,3,RLWEIBUD,rlexu,LO,HI,TILD,TILD,
     1            KEY,LIMIT,SUM,ERRSTD,NEVAL,IER,WORK,IWORK,NT,TETA)
      RETURN
      END
C
C=======================================================================
C
      SUBROUTINE RLTRMNN(ALPHA,BETA,mF)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION QL,QU,BB,DL,DU,GL,GU,PL,PU,rldennor,mF
      EXTERNAL rldennor
      IF (DABS(BETA-0.5).LT.1.D-5) THEN
        mF=ALPHA
        RETURN
      ENDIF
      IF (BETA.LT.1.D-5) THEN
        mF=ALPHA
      ELSE
        BB=DBLE(BETA)
        CALL RLQUANTD(BB,QL)
        CALL RLQUANTD(1.D0-BB,QU)
        DL=rldennor(QL)
        CALL RLGAUSSD(1,QL,PL)
        GL=-DL+DBLE(ALPHA)*PL
        DU=rldennor(QU)
        CALL RLGAUSSD(1,QU,PU)
        GU=-DU+DBLE(ALPHA)*PU
        mF=DBLE(GU-GL)/(1.D0-2.D0*BETA)
      ENDIF
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLEQADN(DD,V,NV,PARAM)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION DD,MM,P1,P2,V(NV),M,PARAM(1)
      BETA =V(1)
      GAM  =V(2)
      DUMMY=PARAM(1)
      CALL RLTRMNN(0.D0,BETA,M)
      MM=DBLE(M)
      CALL RLGAUSSD(1,MM+DD,P1)
      CALL RLGAUSSD(1,MM-DD,P2)
      RLEQADN=P1-P2-DBLE(1.D0-GAM)
      RETURN
      END
C
      SUBROUTINE RLQAD1N(BETA,GAM,TOL,QAD1,ISOL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION RLEQADN,QAD,SMIN,SMAX,TOLD,PARAM(1),Q(2)
      EXTERNAL RLEQADN
C     COMMON/DATI/NQ
      QAD1=0.D0
      ISOL=0
      NQ=2
      Q(1)=BETA
      Q(2)=GAM
      PARAM(1)=0.D0
      CALL RLQUANTD(0.999D0,SMAX)
      SMIN=-SMAX
      TOLD=DBLE(TOL)
      CALL RLRGFLD(RLEQADN,Q,0.D0,SMIN,SMAX,TOLD,100,QAD,ITERM,NQ,PARAM)
      IF (ITERM.NE.1) RETURN
      QAD1=DBLE(QAD)
      ISOL=1
      RETURN
      END
C
      SUBROUTINE RLTRMADN(SIGMA,BETA,GAM,TOL,sF,ISOL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION mF
      IF (DABS(GAM-0.5D0).LT.1.D-5) THEN
        CALL RLQAD1N(BETA,GAM,TOL,sF,ISOL)
      ELSE
        CALL RLTRMNN(0.D0,BETA,mF)
        GAM0=GAM
        CALL RLXERF(2,mF,A0)
        A0=-A0
        CALL RLGAUSSD(1,mF,B0)
 100    CALL RLQAD1N(BETA,GAM0,TOL,QADF,JSOL)
        CALL RLXERF(2,mF+QADF,A1)
        A1=-A1
        CALL RLXERF(2,mF-QADF,A2)
        A2=-A2
        CALL RLGAUSSD(1,mF+QADF,B1)
        CALL RLGAUSSD(1,mF-QADF,B2)
        IF (DABS(GAM-GAM0).LT.1.D-6) THEN
          SUP=(A1 + A2 - 2.D0 * A0) - mF * (B1 + B2 - 2.D0 * B0)
          ISOL=JSOL
          GAM0=1.D0-GAM
          GOTO 100
        ELSE
          SLOW=(A1 + A2 - 2.D0 * A0) - mF * (B1 + B2 - 2.D0 * B0)
          IF (JSOL.EQ.0) ISOL=0
        ENDIF
        sF=(SUP-SLOW)/(1.D0-2.D0*GAM)
      ENDIF
      sF=sF*SIGMA
      RETURN
      END
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION RLLEQNL(L,V,NV,PARAM)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION L,X,X0,XMS,F1,F2,FX,UU,LL,V(NV),PARAM(1)
      U=V(1)
      SIGMA=V(2)
      DUMMY=PARAM(1)
      LL=L
      UU=DBLE(U)
      IF (DABS(LL-UU).LT.0.0001D0) LL=LL-0.0002D0
      X=UU
  100 CALL RLQUANTD(X,X0)
      XMS=X0-DBLE(SIGMA)
      CALL RLGAUSSD(1,XMS,FX)
      IF (DABS(X-UU).LT.1D-6) THEN
        F2=FX
        X=LL
        GOTO 100
      ELSE
        F1=FX
      ENDIF
      XMS=(F2-F1)/(UU-LL)-1.D0
      RLLEQNL=XMS
      RETURN
      END
C
      SUBROUTINE RLQUQLDL(U,ALPHA,SIGMA,TOL,QL,QU,ISOL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION RLLEQNL,QQL,QQU,SMIN,SMAX,TOLD,PARAM(1),Q(2)
      EXTERNAL RLLEQNL,RLXEXPD
C     COMMON/DATI/NQ
      QL=0.D0
      ISOL=0
      NQ=2
      CALL RLQUANTD(DBLE(U),QQU)
      QQU=DBLE(SIGMA)*QQU+DBLE(ALPHA)
      QU=RLXEXPD(DBLE(QQU))
      Q(1)=U
      Q(2)=SIGMA
      PARAM(1)=0.D0
      SMIN=1.D-4
      SMAX=0.98D0
      TOLD=DBLE(TOL)
      CALL RLRGFLD(RLLEQNL,Q,0.D0,SMIN,SMAX,TOLD,100,QQL,ITERM,NQ,PARAM)
      IF (ITERM.NE.1) RETURN
      QQ=DBLE(QQL)
      CALL RLQLNORM(ALPHA,SIGMA,QQ,QL)
      ISOL=1
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLGAUSDD(SIGMA,MU,X)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION X,EXMIN,SPI,SS,T,X2,MU
      DATA NCALL,EXMIN,SPI/0,0.D0,2.506628274631D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
      ENDIF
c      IF (SIGMA.LE.0.) CALL MESSGE(500,'GAUSDD',1)
      RLGAUSDD=0.D0
      SS=SIGMA
      T=(X-MU)/SS
      X2=-T*T/2.D0
      IF (X2.LE.EXMIN) RETURN
      RLGAUSDD=DEXP(X2)/(SS*SPI)
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION RLIFLANS(DX,WGT,N,EXGAU,EXPSI,NT,TETA)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION ANS,EXGAU,EXMIN,DX,IFL(1),Z(1),WGT(N),TETA(NT)
      EXTERNAL EXGAU,EXPSI
c     COMMON/TETAPR/TETA(60)
      DATA NCALL,EXMIN/0,0.D0/
C
C avoid compiler warnings
C
      dummy = expsi(1.D0,1,1,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0,1.D0)
C
C
C  Initializations
C
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(3,EXMIN)
      ENDIF
      Z=0.D0
      IF (Z(1).GT.EXMIN) Z=DEXP(DX)
      RLIFLANS=0.D0
      ANS=EXGAU(WGT(2),WGT(1),DX)
      IF (ANS.LE.1.D-20) RETURN
      ALPHA=WGT(1)
      SIGMA=WGT(2)
      ITC=INT(SNGL(WGT(3)))
      CALL RLIFLOGN(Z,TETA,1,NT,ALPHA,SIGMA,ITC,0,IFL)
      IF (ITC.GE.0) IFL(1)=IFL(1)*IFL(1)
      RLIFLANS=IFL(1)*ANS
      RETURN
      END
C
      SUBROUTINE RLIFLOGN(DX,TETA,NX,NT,ALPHA,SIGMA,ITC,ITT,IFL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION TETA(NT),J1,J0L,J1L,MU1L,DX(NX),IFL(NX)
c     COMMON/TETAPR/TETA(60)
      DATA NCALL,XLGMN,YLGMN/0,0.D0,0.D0/
C
C  Initializations
C
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL RLMACHD(4,XLGMN)
        CALL RLMACHD(5,YLGMN)
      ENDIF
c      IF (ITT.EQ.1) THEN
c        DO 100 I=1,NT
c        TETA(I)=THETA(I)
c  100   CONTINUE
c      ENDIF
      DUMMY=ITT
      DO 500 IX=1,NX
      X=DBLE(DX(IX))
      Y=YLGMN
      IF (X.GT.XLGMN) Y=DLOG(X)
      CALL RLMEDMAD(Y,NT,TETA,TRMNF,TRMDF)
      IF (ITC.EQ.1) XIF=TRMNF
      IF (ITC.EQ.2) XIF=TRMDF
      IF (ITC.EQ.1.OR.ITC.EQ.2) GOTO 300
c IF.alpha.Dsm.l
      ALDSML=TRMNF
      IF (ITC.NE.3) GOTO 104
      XIF=ALDSML
      GOTO 300
c IF.sigma.Dsm.l
  104 SIDSML=TRMDF/TETA(11)
      IF (IABS(ITC).NE.4) GOTO 105
      XIF=SIDSML
      IF (ITC.EQ.-4) XIF=XIF*ALDSML
      GOTO 300
c IF.mu.l
  105 IF (ITC.EQ.5) THEN
        Z=DBLE(ALPHA+0.5*SIGMA**2)
        IFL(IX)=DEXP(Z)*DBLE(ALDSML+SIDSML*SIGMA)
        GOTO 500
      ENDIF
c IF.mu1.l
      MU1L=SIGMA*SIDSML*TETA(35)
c IF.qu1.l
      Z=DBLE(SIDSML)*DBLE(TETA(51))/DBLE(TETA(50))
      QU1L=-DBLE(Z)
      IF (ITC.NE.6) GOTO 107
      XIF=QU1L
      GOTO 300
c IF.qu.Dsm.l
  107 Z=(DBLE(ALDSML)*DBLE(TETA(37))-Z)*DEXP(DBLE(ALPHA))
      QUDSML=DBLE(Z)
      IF (ITC.NE.7) GOTO 108
      XIF=QUDSML
      GOTO 300
c IF.H0.l
  108 quF=TETA(38)
      fquF=TETA(41)
      H0=TETA(44)
      u=TETA(34)
      H0L=-H0
      IF (quF.GE.X) H0L=H0L+1.
      H0L=H0L+fquF*QUDSML
      IF (ITC.NE.8) GOTO 109
      XIF=H0L
      GOTO 300
c IF.H1.l
  109 H1=TETA(45)
      H1L=-H1
      IF (quF.GE.X) H1L=H1L+X
      H1L=H1L+fquF*quF*QUDSML
      IF (ITC.NE.9) GOTO 110
      XIF=H1L
      GOTO 300
c IF.ql1.l
  110 B=(u-TETA(48))*MU1L+(TETA(56)-TETA(55)-TETA(52)*TETA(35))*SIDSML-
     +   TETA(53)*QU1L
      A=TETA(49)*TETA(35)-TETA(54)
      Z=DBLE(B)/DBLE(A)
      QL1L=DBLE(Z)
      IF (ITC.NE.10) GOTO 111
      XIF=QL1L
      GOTO 300
c IF.ql.Dsm.l
  111 Z=(DBLE(ALDSML)*DBLE(TETA(39))+Z)*DEXP(DBLE(ALPHA))
      QLDSML=DBLE(Z)
      IF (ITC.NE.11) GOTO 112
      XIF=QLDSML
      GOTO 300
c IF.J0.l
  112 qlF=TETA(40)
      fqlF=TETA(42)
      J0L=-TETA(46)
      IF (qlF.GE.X) J0L=J0L+1.
      J0L=J0L+fqlF*QLDSML
      IF (ITC.NE.12) GOTO 113
      XIF=J0L
      GOTO 300
c IF.J1.l
  113 J1=TETA(47)
      J1L=-J1
      IF (qlF.GE.X) J1L=J1L+X
      J1L=J1L+fqlF*qlF*QLDSML
      IF (ITC.NE.13) GOTO 114
      XIF=J1L
      GOTO 300
c IF.tcmean.l
  114 H0mJ0=TETA(44)-TETA(46)
      XIF=(-(H1-J1)*(H0L-J0L)/H0mJ0 + H1L-J1L)/H0mJ0
  300 IFL(IX)=DBLE(XIF)
  500 CONTINUE
      RETURN
      END
C
      SUBROUTINE RLAVTCML(TETA,NT,ALPHA,SIGMA,ITC,LOWER,UPPER,TIL,SUM,
     +           IWORK,WORK)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION RLIFLANS,ERRSTD,WORK,RLGAUSDD,LO,HI,SUM,LOWER
      DIMENSION TETA(NT),WGT(3),IWORK(80),WORK(320)
c     COMMON/TETAPR/TETA(60)
      EXTERNAL RLIFLANS,RLGAUSDD,RLEXU
C
C     INITIALISATION
C
      LIMIT=80
      KEY=1
C
      WGT(1)=ALPHA
      WGT(2)=SIGMA
      WGT(3)=DBLE(FLOAT(ITC))
      LO=DBLE(LOWER)
      HI=DBLE(UPPER)
c      DO 100 I=1,NT
c      TETA(I)=THETA(I)
c  100 CONTINUE
      SUM=0.D0
      TILD=DBLE(TIL)
      CALL RLINTGRT(RLIFLANS,WGT,3,RLGAUSDD,RLEXU,LO,HI,TILD,TILD,
     1         KEY,LIMIT,SUM,ERRSTD,NEVAL,IER,WORK,IWORK,NT,TETA)
      RETURN
      END
C
C      SUBROUTINE IF2TST(IOPT,THETA,NT,X,N,Y)
C      DOUBLE PRECISION IFGANS,IFWANS,IFLANS,GAUSDD,RLWEIBUD,RLGAMMAD,XX,Z
C      REAL THETA(NT),X(N),Y(N),WGT(3)
C      COMMON/TETAPR/TETA(60)
C      EXTERNAL IFGANS,IFWANS,IFLANS,GAUSDD,RLWEIBUD,RLGAMMAD,rlexu
CC
C      DO 100 I=1,NT
C      TETA(I)=THETA(I)
C  100 CONTINUE
C      WGT(1)=TETA(2)
C      WGT(2)=TETA(3)
C      WGT(3)=1.
C      IF (IOPT.LT.0) WGT(3)=0
C      if (teta(1).eq.2.) then
C       wgt(1)=teta(25)
C       wgt(2)=teta(26)
C      endif
C      JOPT=IABS(IOPT)
C      DO 200 I=1,N
C      XX=DBLE(X(I))
C      IF (JOPT.EQ.1) Z=IFGANS(XX,WGT,3,RLGAMMAD,rlexu)
C      IF (JOPT.EQ.2) Z=IFWANS(XX,WGT,3,RLWEIBUD,rlexu)
C      IF (JOPT.EQ.3) Z=IFLANS(XX,WGT,3,GAUSDD,rlexu)
C      Y(I)=DBLE(Z)
C  200 CONTINUE
C      RETURN
C      END
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c$$$C
c$$$      SUBROUTINE RLRGFLD(F,V,Y,A,B,TOL,MAXIT,X,ITERM,nv,param)
c$$$C.......................................................................
c$$$      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c$$$      DIMENSION V(nv),param(1)
c$$$      EXTERNAL F
c$$$      DATA TL/1.D-10/
c$$$C-----------------------------------------------------------------------
c$$$C     SOLUTION OF THE SCALAR EQUATION F(X)=Y USING REGULA-FALSI METHOD
c$$$C     USES EXTERNAL FUNCTION F WITH ARGUMENTS F(A,V) WHERE A IS A SCALAR
c$$$C     AND V IS A 3 ELEMENT ARRAY
c$$$C     MF401 = 1 IF MESSGE 401 IS CALLED AND PROGRAM STOPS
c$$$C-----------------------------------------------------------------------
c$$$      MF401=0
c$$$      ITR=1
c$$$      FA=F(A,V,nv,param)-Y
c$$$      FB=F(B,V,nv,param)-Y
c$$$ 20   IF (DABS(FA-FB).GT.TL) GOTO 30
c$$$      MF401=1
c$$$      RETURN
c$$$ 30   XN=(A*FB-B*FA)/(FB-FA)
c$$$      FN=F(XN,V,nv,param)-Y
c$$$      IF (ITR.GE.MAXIT) GOTO 60
c$$$C-----------------------------------------------------------------------
c$$$C     TEST TO SEE IF ROOT HAS BEEN FOUND
c$$$C-----------------------------------------------------------------------
c$$$      IF (DABS(FN).LT.TOL) GOTO 70
c$$$      IF (FA*FN.LE.0.D0) GOTO 40
c$$$      A=XN
c$$$      FA=FN
c$$$      GOTO 50
c$$$ 40   B=XN
c$$$      FB=FN
c$$$C-----------------------------------------------------------------------
c$$$C     INCREMENT ITERATION COUNTER
c$$$C-----------------------------------------------------------------------
c$$$ 50   ITR=ITR+1
c$$$      GOTO 20
c$$$ 60   ITERM=2
c$$$      X=XN
c$$$      RETURN
c$$$ 70   ITERM=1
c$$$      X=XN
c$$$      RETURN
c$$$      END
c
c-----------------------------------------------------------------------
c
c$$$      SUBROUTINE RLSUMLGM(HI,ALPHA,GL)
c$$$C
c$$$C  INT  Log(x)*gamma(x,alpha) dx for x=0 to HI
c$$$C
c$$$      DOUBLE PRECISION ALF,ALPHA,HI,LOGHI,SUM,GA1,GG,GL,GX,
c$$$     1       PREC,SS,TERM
c$$$      DATA NCALL,PREC/0,0.D0/
c$$$      IF (NCALL.EQ.0) THEN
c$$$        NCALL=1
c$$$        CALL RLMACHD(2,PREC)
c$$$      ENDIF
c$$$      GL=0.D0
c$$$      IF (HI.LE.0.D0) RETURN
c$$$      ALF=ALPHA
c$$$      LOGHI=DLOG(HI)
c$$$      CALL RLLGAMAD(ALF+1.D0,GA1)
c$$$      GG=ALF*LOGHI-HI-GA1
c$$$      SS=1.D0/ALF
c$$$      SUM=DEXP(GG+DLOG(SS))
c$$$  100 ALF=ALF+1.D0
c$$$      GG=GG+LOGHI-DLOG(ALF)
c$$$      SS=SS+1.D0/ALF
c$$$      TERM=DEXP(GG+DLOG(SS))
c$$$      SUM=SUM+TERM
c$$$      IF (TERM.GT.PREC) GOTO 100
c$$$      CALL RLINGAMA(HI,ALPHA,GX)
c$$$      GL=LOGHI*GX-SUM
c$$$      RETURN
c$$$      END
c$$$
