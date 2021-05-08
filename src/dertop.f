C   06/10/83 911021123  MEMBER NAME  DERTOP   (N3000.S)     FORTRAN
      SUBROUTINE DERTOP(E0,IE0,U0,CIR,TM6A,TUN)

C
C==========WITH NEW PERTURBED CLOSED ORBIT GET PERTURBED MACHINE========
C          PARAMETERS IN UNCOUPLED DISPERSION FORMALISM.
C          ASSUMES NO HOR./VERT COUPLING.
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      DIMENSION TREV6(6,6),TM6A(6,6),TM6B(6,6),TM6(6,6),TM6I(6,6)
      DIMENSION TRIN6(6,6),RR6(6),RI6(6),VR6(6,6),VI6(6,6)
C      DIMENSION INTGE6(6),WW6(6,6),ZZ(6,6)
      PARAMETER( LWORK=64*6 )   ! now needed for  F02EBF: MV/09.10.2007
      DIMENSION WORK(LWORK)     ! now needed for  F02EBF: MV/09.10.2007
      DIMENSION WW6(6,6),ZZ(6,6)
      DIMENSION WR(6),WI(6),TUN(6),VC6A(6),VC6B(6),EIG6(6,6)
      DIMENSION TM4A(4,4),VC4A(4),VC4B(4),LL(4),MM(4)
      DIMENSION A(3,3),B(3,3),P(3,3),PDIF(3,3),PTEMP(3,3),PTOT(3,3)
      COMPLEX*16 VP,VM,V11N0,V1M1N0,XPH,SPH,EXPN0,ALX
CDP      REAL*8 DNAME/'D'/,IP/'IP'/
      DIMENSION SYM(6,6)
      REAL*8 MUX,MUY,MUS
      DATA N0/47/,QX/47.184/,QS/-0.0537/
C
C
C
      INCLUDE "cloorb.for"
      INCLUDE "cnlist.for"
      INCLUDE "clatic.for"
      INCLUDE "csol.for"
      INCLUDE "cdisp.for"
C
C=====6X6 SYMPLECTIC UNIT MATRIX
      DATA SYM/36*0.D0/
      SYM(1,2)=-1.D0
      SYM(2,1)= 1.D0
      SYM(3,4)=-1.D0
      SYM(4,3)= 1.D0
      SYM(5,6)=-1.D0
      SYM(6,5)= 1.D0
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
      DO 555 II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      CALL DX66(TMAT(1,1,ITY),ITY,IID,II,  XX(ITY),X2(ITY),YY(ITY),TM6,
     +                                                             TM6I)
      CALL JAM666(TREV6,TM6,TREV6)
  555 CONTINUE
C      CALL UCOPY(TREV6,TM6B,72)
      TM6B = TREV6 
C      CALL EIV6(TREV6,TM6A,WR,WI,IERRO)
C      CALL UCOPY(TREV6,TRIN6,72)
      TRIN6 = TREV6
      IFAIL=0
C      CALL F02AGF(TRIN6,6,6,RR6,RI6,VR6,6,VI6,6,INTGE6,IFAIL)
      CALL F02EBF('V',6,TRIN6,6,RR6,RI6,VR6,6,VI6,6,WORK,LWORK,IFAIL) ! <<replmnt
      IF(IFAIL.NE.0)GO TO 9999
      DO 2822 I=1,3
      IT=2*I-1
      DO 2833 J=1,6
      ZZ(J,IT  )  =VR6(J,IT)
 2833 ZZ(J,IT+1)  =VI6(J,IT)
 2822 CONTINUE
      DO 2844 I=1,6
      WR(I)  =RR6(I)
 2844 WI(I)  =RI6(I)
C      CALL UCOPY(ZZ,TM6A,72)
      TM6A = ZZ

C======Get re-ordered eigenvectors and tunes.
      CALL RENORM(TM6A,  6,WR,  WI)


C
C=====GET 3 LOTS OF TWISSES: ASSUME NO COUPLING.
      MUX=DACOS(0.5D0*(TREV6(1,1)+TREV6(2,2)))
      IF(TREV6(1,2) .LT. 0.)MUX=2.D0*PI-MUX
      MUY=DACOS(0.5D0*(TREV6(3,3)+TREV6(4,4)))
      IF(TREV6(3,4) .LT. 0.)MUY=2.D0*PI-MUY
      MUS=DACOS(0.5D0*(TREV6(5,5)+TREV6(6,6)))
      IF(TREV6(5,6) .LT. 0.)MUS=2.D0*PI-MUS
      PSIX=0.D0
      PSIY=0.D0
      PSIS=0.D0
      BETAX=TREV6(1,2)/DSIN(MUX)
      BETAY=TREV6(3,4)/DSIN(MUY)
      BETAS=TREV6(5,6)/DSIN(MUS)
      ALPHAX=(TREV6(1,1)-TREV6(2,2))/(2.D0*DSIN(MUX))
      ALPHAY=(TREV6(3,3)-TREV6(4,4))/(2.D0*DSIN(MUY))
      ALPHAS=(TREV6(5,5)-TREV6(6,6))/(2.D0*DSIN(MUS))
C
C
C
      WRITE(53,103)
      WRITE(53,103)
  103 FORMAT(/,'  ')


      WRITE(53,917)
  917 FORMAT('1',1X,
     + 'EIGENVECTORS AND TUNES FROM THE BETATRON-DISPERSION VERSION',
     + ' OF THE 6X6 FORMALISM')
      DO 38 I=1,6
C      TUN(I)=DACOS(WR(I))/PI2
      TUN(I)=DATAN2(WI(I),WR(I))/PI2
   38 WRITE(53,916)WR(I),WI(I),(TM6A(J,I),J=1,6),TUN(I),DAMP
  916 FORMAT(T6,F12.5,T21,F12.5,T35,'(',6F11.5,')',F15.8,E15.5)
C
C
      WRITE(53,101)BETAX,BETAY,BETAS,ALPHAX,ALPHAY,ALPHAS,MUX,MUY,MUS
  101 FORMAT(' ',9F10.4)
C
C
C
C
C
C
C=====NOW PROPAGATE AROUND THE LATTICE:
C=====GET TWISSES AT END OF LAST ELEMENT. I.E. AT START OF ELEMENT II
      WRITE(53,902)
      S=0.D0
      BXTEMP=BETAX
      BYTEMP=BETAY
      BSTEMP=BETAS
      AXTEMP=ALPHAX
      AYTEMP=ALPHAY
      ASTEMP=ALPHAS
      V11N0 =DCMPLX(0.D0,0.D0)
      V1M1N0=DCMPLX(0.D0,0.D0)
      DO 564 II=2,NELEM
      JTY=ITYPE(II-1)
      IID=ID(JTY)
      IF(IID.EQ.5.OR.IID.EQ.6.OR.IID.EQ.7)GO TO 700
      IF((NAME(JTY)(1:2).EQ.'CQ'.AND.IID.EQ.3) !Ignore artificial lengths. 
     +  .OR.(NAME(JTY)(1:2).EQ.'RQ'.AND.IID.EQ.4).OR.IID.EQ.17)GO TO 700 
       S=S+YY(JTY)
  700 CONTINUE 

      CALL DX66(TMAT(1,1,JTY),JTY,IID,II-1,XX(JTY),X2(JTY),YY(JTY),TM6,
     +                                                             TM6I)
C=====USE SIMILARITY TRANS. TO GET 1-TURN MATRIX AFTER THIS ELEMENT.
C=====FIRST GET THE INVERSE MATRIX USING THE SYMP. UNIT MATRIX.
      CALL JAM666(TREV6,TREV6,TM6I)
      CALL JAM666(TREV6,TM6,TREV6)
C=====WE NOW HAVE THE 1-TURN 6X6 MATRIX AT THE END OF THE ELEMENT II-1
C
C
C
C=====GET 3 LOTS OF TWISSES AGAIN: ASSUME NO COUPLING.
      BETAX=TREV6(1,2)/DSIN(MUX)
      BETAY=TREV6(3,4)/DSIN(MUY)
      BETAS=TREV6(5,6)/DSIN(MUS)
      ALPHAX=(TREV6(1,1)-TREV6(2,2))/(2.D0*DSIN(MUX))
      ALPHAY=(TREV6(3,3)-TREV6(4,4))/(2.D0*DSIN(MUY))
      ALPHAS=(TREV6(5,5)-TREV6(6,6))/(2.D0*DSIN(MUS))
C
C=====GET PHASES:ASSUME THAT THERE ARE NO HOR/VERT COUPLING ELEMENTS.
      SINPSI=TM6(1,2)/DSQRT(BETAX*BXTEMP)
      COSPSI=TM6(1,1)*DSQRT(BXTEMP/BETAX)-AXTEMP*SINPSI
      DPSI=DATAN2(SINPSI,COSPSI)
C     WRITE(53,102)DPSI
C 102 FORMAT(' ','DPSI ',E20.8)
      IF(SINPSI.LT.-1.D-20)DPSI=2.D0*PI+DPSI
      PSIX=PSIX+DPSI/2.D0/PI
C
      SINPSI=TM6(3,4)/DSQRT(BETAY*BYTEMP)
      COSPSI=TM6(3,3)*DSQRT(BYTEMP/BETAY)-AYTEMP*SINPSI
      DPSI=DATAN2(SINPSI,COSPSI)
C     WRITE(53,102)DPSI
      IF(SINPSI.LT.-1.D-20)DPSI=2.D0*PI+DPSI
      PSIY=PSIY+DPSI/2.D0/PI
C=====SYNC.PHASE STEPS ARE ALWAYS SMALL AND CAN BE NEGATIVE.
      SINPSI=TM6(5,6)/DSQRT(BETAS*BSTEMP)
      COSPSI=TM6(5,5)*DSQRT(BSTEMP/BETAS)-ASTEMP*SINPSI
      DPSI=DATAN2(SINPSI,COSPSI)
C     WRITE(53,102)DPSI
C     IF(SINPSI.LT.-1.D-8)DPSI=2.D0*PI+DPSI
      PSIS=PSIS+DPSI/2.D0/PI
C
      BXTEMP=BETAX
      BYTEMP=BETAY
      BSTEMP=BETAS
      AXTEMP=ALPHAX
      AYTEMP=ALPHAY
      ASTEMP=ALPHAS
      IF(IID.NE.5)GO TO 563
C
C=====CALCULATE  'LINEAR' SYNCHRO-BETA STOP BANDS ALA G.RIPKEN
      A1010=XX(JTY)*D1(II-1)
      A0110=XX(JTY)*D2(II-1)
C     WRITE(53,5577)A1010,A0110,XX(JTY),DXP(II-1),DX(II-1)
C5577 FORMAT(' ',5F20.10)
      CHIX=2.D0*PI*(PSIX-QX*S/CIR)
      CHIS=2.D0*PI*(PSIS-QS*S/CIR)
      EN0=N0*2.D0*PI*S/CIR
      XPH=  DCMPLX(DCOS(CHIX), DSIN(CHIX))
      SPH=  DCMPLX(DCOS(CHIS), DSIN(CHIS))
      SPHM= DCMPLX(DCOS(CHIS),-DSIN(CHIS))
      EXPN0=DCMPLX(DCOS(EN0),  DSIN(EN0))
      ALX=DCMPLX(ALPHAX,-1D0)
      VP=0.5D0*XPH*SPH*EXPN0
      VP=VP*(A1010*DSQRT(BETAX*BETAS)-A0110*DSQRT(BETAS/BETAX)*ALX)
      V11N0=V11N0+VP
      VM=0.5D0*XPH*SPHM*EXPN0
      VM=VM*(A1010*DSQRT(BETAX*BETAS)-A0110*DSQRT(BETAS/BETAX)*ALX)
      V1M1N0=V1M1N0+VM
C
  563 IF(S.LT.10000.D0.OR.S.GT.CIR-10.)
     +WRITE(53,93)NAME(ITYPE(II-1)),S,BETAX,ALPHAX,BETAY,ALPHAY,BETAS,
     +    ALPHAS,PSIX,PSIY,PSIS,NAME(ITYPE(II-1))
   93 FORMAT(1X,A8,10F10.6,1X,A4)

      WRITE(55,'(I10,F20.2,2E20.4)')II,S,PSIS


C
C
C
C
C
  564 CONTINUE
C
C
C=====GET STOPBAND WIDTHS.
  565 FORMAT('0','STOP BAND WIDTHS--',2E14.6)
      DEL0P=CDABS(V11N0)/PI
      DEL0M=CDABS(V1M1N0)/PI
      WRITE(53,565)DEL0P,DEL0M
C
C
C
C
C
C
C
C
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
      CALL DX66(TMAT(1,1,ITY),ITY,IID,II,  XX(ITY),X2(ITY),YY(ITY),TM6,
     +                                                             TM6I)
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
C
C
C
C
      RETURN
 9999 WRITE(53,921)
  902 FORMAT('0',1X,'HOR.,VERT. & LONG. TWISS PARAMETERS AT END OF EA',
     +    'CH ELEMENT',///,1X,'NAME',T13,'S______',T23,'BETAX___',T34,
     +    'ALPHAX__',T43,'BETAY___',T54,'ALPHAY__',T64,'BETAS___'
     +    ,T74,'ALPHAS__',T84,'PSIX____',T94,'PSIY____',T104,'PSIS----',
     +    /)
  921 FORMAT(' ERROR IN EIGEN VALUE ROUTINE')
      STOP
  909 FORMAT(///,' MAX.DISP.FUNCTIONS: ETAXMAX =',F10.4,' M.    AT',
     +   F10.4,' M.',/,
     +   T22,'ETAYMAX =',F10.4,' M.    AT',F10.4,' M.',
     +   /,' RMS DISP.FUNCTIONS: ETAXRMS =',
     +    F10.4,' M.',/,T22,'ETAYRMS =',F10.4,' M.')
  910 FORMAT(1H1,1X,'PERTURBED DISPERSION FUNCTIONS AT END OF ELEMENT:'
     +                                               ,///,3X,'NAME--'
     +    ,T15,'POS--',T26,'ETAX--',T37,'ETAXP--',T48,
     +                                         'ETAY--',T59,'ETAYP--',/)
  912 FORMAT(1H1,' PERTURBED MACHINE PARAMETERS:',/,T6,'MOM.COMP.',
     >    'FACTOR =',F13.6,/,T5,'SYN.RAD.LOSS =',F13.6,'  MEV')
  913 FORMAT(///,' PERTURBED TOTAL TRANSF. AROUND THE 1-ST ELEMENT')
  915 FORMAT(///,' PERTURBED EIGEN-TUNES AND EIGENVECTORS:',/,T6,
     >    'COS(2*PI*NU)',T21,'SIN(2*PI*NU)',T36,'EIGENVECTORS',T106,'
     >TUNES          DAMPING?')
      END