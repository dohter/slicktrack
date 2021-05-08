C   06/10/83 510231710  MEMBER NAME  PERTOP   (SEPT95.S)    FORTRAN
      SUBROUTINE PERTOP(E0,IE0,U0,CIR,TM6A,TUN)
C
C==========WITH NEW PERTURBED CLOSED ORBIT GET PERTURBED MACHINE========
C                 PARAMETERS & DISPERSION FUNCTIONS
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      DIMENSION TREV6(6,6),TREV66(6,6),TM6A(6,6),TM6B(6,6),TM6(6,6)
      DIMENSION TRIN6(6,6),RR6(6),RI6(6),VR6(6,6),VI6(6,6)
C      DIMENSION INTGE6(6),WW6(6,6),ZZ(6,6)
      PARAMETER( LWORK=64*6 )   ! now needed for  F02EBF: MV/09.10.2007
      DIMENSION WORK(LWORK)     ! now needed for  F02EBF: MV/09.10.2007
      DIMENSION WW6(6,6),ZZ(6,6)
      DIMENSION WR(6),WI(6),TUN(6),VC6A(6),VC6B(6),EIG6(6,6)
      DIMENSION TRIN66(6,6)
      DIMENSION TM4A(4,4),VC4A(4),VC4B(4),LL(4),MM(4)
      DIMENSION A(3,3),B(3,3),P(3,3),PDIF(3,3),PTEMP(3,3),PTOT(3,3)

C======Variables with ``NC'' are for getting at objects for no cavities.
      DIMENSION ZZNC(6,6),TM6NCA(6,6),RRNC6(6),RINC6(6),VRNC6(6,6),
     +        VINC6(6,6),WRNC(6),WINC(6),TUNNC(6),TM6NCB(6,6)

CDB      REAL*8 DNAME/'D'/,IP/'IP'/
C
C
C
      INCLUDE "cloorb.for"
      INCLUDE "cnlist.for"
      INCLUDE "clatic.for"
      INCLUDE "csol.for"
      INCLUDE "cdisp.for"
C
C
C
      PI=3.1415926535897932D0
      PI2=2.D0*PI
C
C      
C
      IT=ITYPE(1)
      DO 10 I=1,6
      DO 10 J=1,6
   10 TREV6(I,J)=TMAT(I,J,IT)
C      CALL UCOPY(TREV6,TREV66,72)
      TREV66 = TREV6
      DO 555 II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      CALL MX66(TMAT(1,1,ITY),ITY,IID,II,  XX(ITY),X2(ITY),YY(ITY),TM6)

C      WRITE(53,913)
C      DO 45 I=1,6
C   45 WRITE(53,'(I10,3X,A8,6F12.5)')II,NAME(ITY),(TM6(I,J),J=1,6)

      CALL JAM666(TREV6,TM6,TREV6)
      IF(IID.EQ.5) GOTO 555               ! Get matrix without cavities --> get
      CALL JAM666(TREV66,TM6,TREV66)      ! betas in uncoupled regions.


      IF(II.EQ.NELEM)THEN 
CCC      WRITE(53,913)
CCC      DO 55 I=1,6
CCC   55 WRITE(53,'(I10,3X,A8,6F12.5)')II,NAME(ITY),(TREV6(I,J),J=1,6)
      ENDIF 

C      IF(II.EQ.1700)STOP

  555 CONTINUE

C      IF(1.EQ.1)STOP

C      CALL UCOPY(TREV6,TM6B,72)
      TM6B = TREV6    

C      IF(1.EQ.1)STOP

C     CALL EIV6(TREV6,TM6A,WR,WI,IERRO)
C     IF(IERRO .NE. 0) GOTO 9999
C=====USE THE NAG ROUTINE. THIS IS MORE ROBUST THAN EIV6.
C      CALL UCOPY(TREV6,TRIN6,72)
      TRIN6 = TREV6
      IFAIL=0
C      CALL F02AGF(TRIN6,6,6,RR6,RI6,VR6,6,VI6,6,INTGE6,IFAIL)
      CALL F02EBF('V',6,TRIN6,6,RR6,RI6,VR6,6,VI6,6,WORK,LWORK,IFAIL) ! <<replmnt
      IF(IFAIL.NE.0)GO TO 9999
C      CALL UCOPY(TREV66,TRIN66,72)
      TRIN66 = TREV66

c======Repeat for case of no cavities.
      IFAIL=0
C      CALL F02AGF(TRIN66,6,6,RRNC6,RINC6,VRNC6,6,VINC6,6,INTGE6,IFAIL)
      CALL F02EBF('V',6,TRIN66,6,RRNC6,RINC6,VRNC6,6,VINC6,6,
     + WORK,LWORK,IFAIL)
      IF(IFAIL.NE.0)GO TO 9999
C=====WRITE OUT F02AGF RESULTS
C     WRITE(53,9266)
C9266 FORMAT(' ','NOW THE F02AGF RESULTS: 8X8')
      DO 2822 I=1,3
      IT=2*I-1
      DO 2833 J=1,6
      ZZNC(J,IT  )=VRNC6(J,IT)
      ZZNC(J,IT+1)=VINC6(J,IT)
      ZZ(J,IT  )  =VR6(J,IT)
 2833 ZZ(J,IT+1)  =VI6(J,IT)
 2822 CONTINUE
      DO 2844 I=1,6
      WRNC(I)=RRNC6(I)
      WINC(I)=RINC6(I)
      WR(I)  =RR6(I)
 2844 WI(I)  =RI6(I)
C      CALL UCOPY(ZZ,TM6A,72)
      TM6A = ZZ 
C      CALL UCOPY(ZZNC,TM6NCA,72)
      TM6NCA = ZZNC
C
      DO 318 I=1,6
C     TUN(I)=DACOS(WR(I))/PI2
      TUN(I)=DATAN2(WI(I),WR(I))/PI2
      DAMP=-DLOG(DSQRT(WR(I)*WR(I)+WI(I)*WI(I)))
  318 WRITE(53,916)WR(I),WI(I),(TM6A(J,I),J=1,6),TUN(I),DAMP




C======Get re-ordered eigenvectors and tunes.
      CALL RENORM(TM6A,  6,WR,  WI)
C      CALL RENORM(TM6NCA,6,WRNC,WINC)
CC
CCC      WRITE(53,915)
CCC      DO 388 I=1,6
CCCC     TUN(I)=DACOS(WR(I))/PI2
CCC      TUN(I)=DATAN2(WI(I),WR(I))/PI2
CCC      DAMP=-DLOG(DSQRT(WR(I)*WR(I)+WI(I)*WI(I)))
CCC 388  WRITE(53,916)WR(I),WI(I),(TM6A(J,I),J=1,6),TUN(I),DAMP
C=====TEST FOR SYMPLECTICITY
CCC      CALL SYMP6(TM6B)

C      IF(1.EQ.1)STOP


C
C=====WRITE THE DISPERSION FUNCTION
      IF(IE0.EQ.1.AND.IDISP.EQ.1)WRITE(53,910)
      DO 558 I=1,4
      DO 559 J=1,4
      TM4A(I,J)=TREV66(I,J)
  559 IF(I.EQ.J)TM4A(I,J)=TREV66(I,J)-1.D0
  558 VC4B(I)=-TREV66(I,6)
      CALL SIMQ(TM4A,VC4B,4,IERRO)
C=====TEST USING MINV TO SOLVE FOR DISPERSION INSTEAD.
CDB   CALL DMINV(TM4A,4,DET,LL,MM)
CDB   DO 5599 I=1,4
CDB   VC4B(I)=0.D0
CDB   DO 5599 J=1,4
CDB   VC4B(I)=-TM4A(I,J)*TREV66(J,6)+VC4B(I)
CDB 5599 CONTINUE
C
C
      IF(IERRO.NE.0)STOP
      S=0.
      IF(IE0.EQ.1)WRITE(53,911)NAME(ITYPE(1)),S,VC4B
      ALPHA1=-TREV66(5,6)/CIR
      DO 569 I=1,4
  569 ALPHA1=ALPHA1-TREV66(5,I)*VC4B(I)/CIR
      EXMAX=DABS(VC4B(1))
      EYMAX=DABS(VC4B(3))
      EXRMS=VC4B(1)**2
      EYRMS=VC4B(3)**2
      D1(1)=VC4B(1)
      D2(1)=VC4B(2)
      D3(1)=VC4B(3)
      D4(1)=VC4B(4)
C
C
C

C      IF(1.EQ.1)STOP

C
C
C
C=====NOW PROPAGATE AROUND THE LATTICE: GET DISP.AT END OF LAST ELEMENT.
C=====                                    I.E. AT START OF ELEMENT II
      S=0.
      DO 564 II=2,NELEM
      JTY=ITYPE(II-1)
      IID=ID(JTY)
      IF(IID.EQ.5.OR.IID.EQ.6.OR.IID.EQ.7.OR.IID.EQ.17)GO TO 699
      IF( (NAME(JTY)(1:2).EQ.'CQ'.AND.IID.EQ.3) !Ignore artificial lengths. 
     +     .OR.(NAME(JTY)(1:2).EQ.'RQ'.AND.IID.EQ.4))GO TO 699 
       S=S+YY(JTY)
  699 CONTINUE  

      CALL MX66(TMAT(1,1,JTY),JTY,IID,II-1,XX(JTY),X2(JTY),YY(JTY),TM6)
      IF(IID.EQ.5)GO TO 574
      CALL JAM666(TM6NCA,TM6,TM6NCA)    !propagate the eigenvectors to get the
                                        !betas even with coupling. 
      BETACX(II) = 2.D0*     (TM6NCA(1,1)**2 + 1.D0*TM6NCA(1,2)**2)
      BETACY(II) = 2.D0*     (TM6NCA(3,3)**2 + 1.D0*TM6NCA(3,4)**2)
      DO 566 I=1,4
      VC4A(I)=TM6(I,6)
      DO 566 J=1,4
  566 VC4A(I)=VC4A(I)+TM6(I,J)*VC4B(J)
      DO 565 I=1,4
  565 VC4B(I)=VC4A(I)
      IF(DABS(VC4B(1)).GT.EXMAX)SMAXX=S
      IF(DABS(VC4B(1)).GT.EXMAX)EXMAX=DABS(VC4B(1))
      IF(DABS(VC4B(3)).GT.EYMAX)SMAXY=S
      IF(DABS(VC4B(3)).GT.EYMAX)EYMAX=DABS(VC4B(3))
C  574 IF(IE0.EQ.1.AND.IDISP.EQ.1.AND.(S.LT. 250..OR.S.GT.6000.)  ! 27/8/03
  574 IF(IE0.EQ.1.AND.IDISP.EQ.1                               
     +.AND.NAME(JTY).NE.'DRIFT')WRITE(53,911)NAME(JTY),S,VC4B,
     +                                            BETACX(II), BETACY(II)
  911 FORMAT(4X,A8,T16,7F11.5)
C=====STORE DISPERSION:  USE OLD DISPERSION IF LAST ELEMENT WAS CAVITY.
      D1(II)=VC4B(1)
      D2(II)=VC4B(2)
      D3(II)=VC4B(3)
      D4(II)=VC4B(4)
      EXRMS=EXRMS+VC4B(1)**2
      EYRMS=EYRMS+VC4B(3)**2
  564 CONTINUE
C
C
C
C      IF(1.EQ.1)STOP
C
C
C=====SUMMARISE THE PERTURBED MACHINE FUNCTIONS
      WRITE(53,912)ALPHA1,U0
      WRITE(53,913)
      DO 35 I=1,6
   35 WRITE(53,914)(TM6B(I,J),J=1,6)
  914 FORMAT(6F12.5)

      WRITE(53,9133)
      DO 355 I=1,6
  355 WRITE(53,914)(TREV66(I,J),J=1,6)


      WRITE(53,915)
      DO 38 I=1,6
C     TUN(I)=DACOS(WR(I))/PI2
      TUN(I)=DATAN2(WI(I),WR(I))/PI2
      DAMP=-DLOG(DSQRT(WR(I)*WR(I)+WI(I)*WI(I)))
   38 WRITE(53,916)WR(I),WI(I),(TM6A(J,I),J=1,6),TUN(I),DAMP
  916 FORMAT(T6,F12.5,T21,F12.5,T35,'(',6F11.5,')',F15.8,E15.5)
C=====TEST FOR SYMPLECTICITY
      CALL SYMP6(TM6B)

C======2/2004: Kill writing this: the normalisation in RENORM goes crazy and it's of no 
C      real interest as far as I can see right now.
C      WRITE(53,9155)
C      DO 388 I=1,6
C      TUNNC(I)=DATAN2(WINC(I),WRNC(I))/PI2
C      DAMPNC=-DLOG(DSQRT(WRNC(I)*WRNC(I)+WINC(I)*WINC(I)))
C  388 WRITE(53,916)WRNC(I),WINC(I),(TM6NCB(J,I),J=1,6),TUNNC(I),DAMPNC




      EXRMS=DSQRT(EXRMS/NELEM)
      EYRMS=DSQRT(EYRMS/NELEM)
      WRITE(53,909)EXMAX,SMAXX,EYMAX,SMAXY,EXRMS,EYRMS
  561 CONTINUE
C
C
C
C
C
C      IF(1.EQ.1)STOP
C
C
C
C
C=====GO ONCE MORE  ROUND  THE  RING  TO  EXTRACT  PHASE  ADVANCES  FROM
C=====6-ELEMENT EIGEN-VECTORS----MIGHT BE USEFUL.  EVEN  IF  THE  1-TURN
C=====MATRIX IS COUPLED,THE FRACTIONAL ACCUMULATED PHASE FROM THE CHOSEN
C=====COLUMN ELEMENT SHOULD AGREE WITH  THE  EIGENVALUE  SINCE  DEFINING
C=====START PHASE AS ZERO IS SAME AS DEFINING CHOSEN COLUMN  ELEMENT  AS
C=====REAL.FROM EACH EIGEN-VECTOR TYPE EXTRACT 3  PHASE QUANTITIES  FROM
C=====ROWS 1,3,5. PRESUME THESE PHASES ARE RELATED TO 'PETROS TYPE' 4
C=====PHASE IDEAS.
C
C
C      CALL UCOPY(TM6A,EIG6,72)
      EIG6 = TM6A
      DO 105 I=1,3
      DO 105 J=1,3
      A(J,I)=EIG6(2*J-1,2*I-1)
      B(J,I)=EIG6(2*J-1,2*I)
      P(J,I)=0.D0
      IF(A(J,I).NE.0.D0)P(J,I)=DATAN2(B(J,I),A(J,I))
      IF(B(J,I).LT.0.D0.AND.A(J,I).NE.0.D0)P(J,I)=PI2+P(J,I)
      PTEMP(J,I)=P(J,I)
  105 CONTINUE
C
C
      IT=ITYPE(1)
C      CALL VZERO(PTOT,18)
      PTOT = 0.D0
C     DO 101 I=1,6
C     DO 101 J=1,6
C 101 TREV6(I,J)=TMAT(I,J,IT)
      DO 104 II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      CALL MX66(TMAT(1,1,ITY),ITY,IID,II,  XX(ITY),X2(ITY),YY(ITY),TM6)
      CALL JAM666(EIG6,TM6,EIG6)
C
  110 DO 116 I=1,3
      DO 106 J=1,3
      A(J,I)=EIG6(2*J-1,2*I-1)
      B(J,I)=EIG6(2*J-1,2*I)
      P(J,I)=0.D0
      IF(A(J,I).NE.0.D0)P(J,I)=DATAN2(B(J,I),A(J,I))
      IF(P(J,I).LT.0.D0)P(J,I)=PI2+P(J,I)
      PDIF(J,I)=P(J,I)-PTEMP(J,I)
      IF(PTEMP(J,I).GT.P(J,I).AND.DABS(PTEMP(J,I)-P(J,I)).GT.PI)
     +                                           PDIF(J,I)=PDIF(J,I)+PI2
      IF(PTEMP(J,I).LT.P(J,I).AND.DABS(PTEMP(J,I)-P(J,I)).GT.PI)
     +                                           PDIF(J,I)=PDIF(J,I)-PI2
      PTOT(J,I)=PTOT(J,I)+PDIF(J,I)/PI2
      PTEMP(J,I)=P(J,I)
  106 CONTINUE
  116 CONTINUE
C
  104 CONTINUE
C
      WRITE(53,1107)
 1107 FORMAT('0',//////,' ACCUMULATED PHASE ADVANCES:')
      DO 111 JJ=1,3
  111 WRITE(53,1106)JJ,(PTOT(I,JJ),I=1,3)
 1106 FORMAT('0',4X,'EIGENVECTOR',I2,4X,3F12.6)
C
C
C
C
C      IF(1.EQ.1)STOP
C
C
C
      RETURN
 9999 WRITE(53,921)
  921 FORMAT(' ERROR IN EIGEN VALUE ROUTINE')
      STOP
  909 FORMAT(///,' MAX.DISP.FUNCTIONS: ETAXMAX =',F10.4,' M.    AT',
     +   F10.4,' M.',/,
     +   T22,'ETAYMAX =',F10.4,' M.    AT',F10.4,' M.',
     +   /,' RMS DISP.FUNCTIONS: ETAXRMS =',
     +    F10.4,' M.',/,T22,'ETAYRMS =',F10.4,' M.')
  910 FORMAT(1H1,1X,'PERTURBED DISPERSION FUNCTIONS AND APPROX. BETAS AT 
     + END OF ELEMENT:'
     +                                               ,///,3X,' NAME--'
     +    ,T20,'POS--',T31,'ETAX--',T42,'ETAXP--',T53,
     +                  'ETAY--',T64,'ETAYP--',T75,'BETAX--',T86,
     +                                                   'BETAY--'/)
  912 FORMAT(1H1,' PERTURBED MACHINE PARAMETERS:',/,T6,'MOM.COMP.',
     >    'FACTOR =',F13.6,/,T5,'SYN.RAD.LOSS =',F13.6,'  MEV')
  913 FORMAT(///,' PERTURBED TOTAL TRANSF. AROUND THE 1-ST ELEMENT')
 9133 FORMAT(///,' PERTURBED TOTAL TRANSF. AROUND THE 1-ST ELEMENT', 
     +             '--WITHOUT CAVITIES')



  915 FORMAT(///,' PERTURBED EIGEN-TUNES AND EIGENVECTORS:',/,T6,
     >    'COS(2*PI*NU)',T21,'SIN(2*PI*NU)',T36,'EIGENVECTORS',T106,'
     >TUNES          DAMPING?')

 9155 FORMAT(///,' PERTURBED EIGEN-TUNES AND EIGENVECTORS--NO CAVITIES:'
     +,/,T6,
     +    'COS(2*PI*NU)',T21,'SIN(2*PI*NU)',T36,'EIGENVECTORS',T106,'
     +TUNES          DAMPING?')



      END
