
C   06/10/83 603081617  MEMBER NAME  LINOPT   (SEPT95.S) M  FORTRAN
      SUBROUTINE LINOPENOPT(E0,IE0,U0,CIR)
C
C
C===========LINEAR MACHINE FUNCTIONS AROUND THE RING.=================
C
C 
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER *1 TEXT(80)
C
      DIMENSION TREV(7,7),TREV77(7,7),TEMP(7,7),TM7(7,7)
      DIMENSION TM6A(6,6),TM6B(6,6),WR(6),WI(6), DVEC(6)
      DIMENSION TRIN6(6,6),RR6(6),RI6(6),VR6(6,6),VI6(6,6)
      DIMENSION INTGE6(6),WW6(6,6),ZZ(6,6)
      COMPLEX*16 TVEC(6) 
      REAL*8 MUX,MUX0,MUY, MUY0
      DOUBLE PRECISION BETAX, ALPHAX, BETAY, ALPHAY, ALPHA0 

      DATA IH/0/,IV/0/,IPU/0/
C
C
      INCLUDE "cdesy2.for"
      INCLUDE "cnlist.for"
      INCLUDE "clatic.for"
      INCLUDE "csol.for"
      INCLUDE "csex.for"
      INCLUDE "cintpt.for"

C
      PI=3.1415926535897932D0
C
       I131=1

      WRITE(53,103)
      WRITE(53,103)
  103 FORMAT(/,'  ')
      WRITE(53,929)
  929 FORMAT('1','Entering Subroutine LINOPENOPT') 
      WRITE(53,103)
      WRITE(53,103)

C=====================Input Courant-Snyder Parameters=====

      IT=ITYPE(1)


      READ(5,'(80A1)')TEXT 
      WRITE(53,'(80A1)')TEXT 
C      WRITE(53,'(A)')' Externally supplied initial beam parameters:'
      READ(5,'(30A1,F20.14)') (TEXT(KT),KT=1,30),BETAX
      WRITE(53,'(A,30A1,F20.14)')' ',(TEXT(KT),KT=1,30),BETAX
      READ(5,'(30A1,F20.14)') (TEXT(KT),KT=1,30),ALPHAX
      WRITE(53,'(A,30A1,F20.14)')' ',(TEXT(KT),KT=1,30),ALPHAX
      READ(5,'(30A1,F20.14)') (TEXT(KT),KT=1,30),MUX
      WRITE(53,'(A,30A1,F20.14)')' ',(TEXT(KT),KT=1,30),MUX
      READ(5,'(30A1,F20.14)') (TEXT(KT),KT=1,30),BETAY
      WRITE(53,'(A,30A1,F20.14)')' ',(TEXT(KT),KT=1,30),BETAY
      READ(5,'(30A1,F20.14)') (TEXT(KT),KT=1,30),ALPHAY
      WRITE(53,'(A,30A1,F20.14)')' ',(TEXT(KT),KT=1,30),ALPHAY
      READ(5,'(30A1,F20.14)') (TEXT(KT),KT=1,30),MUY
      WRITE(53,'(A,30A1,F20.14)')' ',(TEXT(KT),KT=1,30),MUY

      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,902)



C==================================
      TVEX1R=DSQRT(BETAX*0.5D0)
      TVEY1R=DSQRT(BETAY*0.5D0) 
      TVEX2R=-ALPHAX/DSQRT(BETAX*2.D0)
      TVEX2I=-1.D0/DSQRT(BETAX*2.D0)
      TVEY2R=-ALPHAY/DSQRT(BETAY*2.D0)
      TVEY2I=-1.D0/DSQRT(BETAY*2.D0)

      TVEC(1)=DCMPLX(TVEX1R,0.D0)
      TVEC(2)=DCMPLX(TVEX2R,TVEX2I) 
      TVEC(3)=DCMPLX(TVEY1R,0.D0)
      TVEC(4)=DCMPLX(TVEY2R,TVEY2I) 
      TVEC(5)=DCMPLX(0.D0,0.D0)
      TVEC(6)=DCMPLX(0.D0,0.D0)

C=====Set initial dispersion to ZERO for now==========   
                                                                              
      Do 3 I=1,5
 3    DVEC(I)=0.0
      DVEC(6)=1.D0   
      DVEC(1)=-0.029153
      DVEC(2)= 0.000680
C=========== MAIN CYCLE=(131)===========================
      IMET=0
      IMET2=0
      IMET1=0
      S = 0.D0

   15 DO 131 I=2,NELEM
      NELEM1=NELEM-5
      ITY=ITYPE(I)
      IID=ID(ITY) 
      IF (IID.EQ.3) IMET=IMET+1
      IF (IID.EQ.2.OR.IID.EQ.15) IMET2=IMET2+1 
      IF (IID.EQ.1.AND.I131.GT.379) IMET1=IMET1+1 


      IF(IID.EQ.5.OR.IID.EQ.6.OR.IID.EQ.7.OR.
     +  (IID.EQ.3.AND.NAME(ITY)(1:2).EQ.'CQ').OR.
     +  (IID.EQ.4.AND.NAME(ITY)(1:2).EQ.'RQ').OR.
     +  (IID.EQ.3.AND.NAME(ITY)(1:1).EQ.'E'))GO TO 777
      S = S + YY(ITY)
  777 CONTINUE      
  
      ALPHA0=ALPHAX
      GO TO(130,130,130,132,130,1322,1322,132,130, 132,130,130,130,130,
     +                                                 130,130,130),IID
      IF(IID.EQ.99)GO TO 131
C  130 CALL UCOPY(TMAT(1,1,ITY),TEMP,98)
  130 TEMP  =  TMAT(1:7,1:7,ITY)
      GOTO 133
  132 CALL UNIT(7,TEMP)
      IF(NAME(ITY)(1:2).NE.'RQ')THEN
      TEMP(1,2)=YY(ITY)
      TEMP(3,4)=YY(ITY)
      ENDIF
      GOTO 133
 1322 CALL UNIT(7,TEMP)
  133 CONTINUE


C========FIND PHASE ADVANCE==========
       DO 174 I74=1,6
       DO 174 J74=1,6
       TM6B(I74,J74)=TEMP(I74,J74)
 174   CONTINUE  
        TVEC=MATMUL(TM6B,TVEC)
        DVEC=MATMUL(TM6B,DVEC)
        I131=I131+1

C===== Find phase advance using old alpha and beta(see p.49 handbk)

        CTM6B=BETAX*TM6B(1,1)-ALPHAX*TM6B(1,2)
        MUX=DATAN2(TM6B(1,2),CTM6B)/PI/2.D0
        MUX0=MUX0+MUX
        CTM6C=BETAY*TM6B(1,1)-ALPHAY*TM6B(1,2)
        MUY=DATAN2(TM6B(1,2),CTM6C)/PI/2.D0
        MUY0=MUY0+MUY
        

C============ NEW Courant-Snyders parameters=======

        WR(1)=DBLE(TVEC(2)/TVEC(1))
        WR(2)=DIMAG(TVEC(2)/TVEC(1))
        ALPHAX=WR(1)/WR(2)
        BETAX=-1.D0/WR(2)
        WR(3)=DBLE(TVEC(4)/TVEC(3))
        WR(4)=DIMAG(TVEC(4)/TVEC(3))
        ALPHAY=WR(3)/WR(4)
        BETAY=-1.D0/WR(4)
        YD=ALPHAX-ALPHA0


      ETAX =  DVEC(1)
      ETAXP = DVEC(2)
      ETAY =  DVEC(3)
      ETAYP = DVEC(4)
      PSIX  = MUX0
      PSIY  = MUY0

C      IF(I.EQ.2.OR.I.EQ.NELEM)THEN
      IF(NAME(ITYPE(I-0))(1:2).EQ.'IP')THEN
      WRITE(53,93)NAME(ITYPE(I-0)),S,BETAX,BETAY,ALPHAX,
     +    ALPHAY,ETAX,ETAXP,ETAY,ETAYP,PSIX,PSIY,NAME(ITYPE(I-0))
      ELSEIF(NAME(ITYPE(I-1))(1:5).NE.'DRIFT') THEN
      IF(IPTBET.EQ.1)WRITE(53,93)NAME(ITYPE(I-1)),S,BETAX,BETAY,ALPHAX,
     +    ALPHAY,ETAX,ETAXP,ETAY,ETAYP,PSIX,PSIY,NAME(ITYPE(I-1))
      ENDIF


  131 CONTINUE

      Do 456 I=1,6
      RR6(I)=DBLE(TVEC(I))
      RI6(I)=DIMAG(TVEC(I))
 456  WRITE (53,193)RR6(I),RI6(I)                   
C===13/03/07--LINAC!
       WRITE (53,202) ALPHAX,BETAX,MUX0,
     + DVEC(1),
     + DVEC(2), ALPHAY,BETAY,MUY0,DVEC(3),DVEC(4),I131                                                                                      
                                     

 
C      CALL SYMP6(TM6B)

      RETURN
 9999 WRITE(53,92)
   92 FORMAT(' ERROR IN EIGEN VALUE ROUTINE')
      STOP
C
C
C
C
C
 193   FORMAT (//, 'REAL=',ES17.9,TR5,'IM=',ES17.9)  

 202   FORMAT (/,ES16.9,ES16.9,ES16.9,ES16.9,ES16.9,
     + /,ES16.9,ES16.9,ES16.9,ES16.9,ES16.9,/,'N=',I6)

 902   FORMAT(' UNPERTURBED MACHINE FUNCTIONS AT END OF EA',
     >    'CH ELEMENT',///,1X,'NAME',T14,'S______',T25,'BETAX___',T35,
     >    'BETAY___',T45,'ALPHAX__',T55,'ALPHAY__',T65,'ETAX____'
     >    ,T75,'ETAXP___',T85,'ETAY____',T95,'ETAYP___',T105,'PSIX----',
     >    T115,'PSIY----',/)
   93 FORMAT(1X,A8,11F10.4,1X,A4)
C  904 FORMAT(1H1,T4,'UNPERTURBED MACHINE PARAMETERS:',/,T5,'TOTAL ',
C     +    '# OF BEAM-LINE ELEMENTS=',I6,/,T5,'NUX=',F9.5,T23,'NUY=',
C     +    F9.5,/,T5,'BETAX*=',F9.4,' M',/,T5,
C     +    'BETAY*=',F9.4,' M',/,T5,'BEAM ENERGY=',F9.4,' GEV',/,T5,
C     +    'CIRCUMFERENCE=',F10.4,T30,'M',/,T5,'MOM.COM.FACTOR=',
C     +   F9.6,/,T5,'SYN.RAD.LOSS=',F13.6,'  MEV',/,T5,'CHROMATICITIES'
C    +    ,':',T25,'NATURAL',T36,'SEX.CONTRIBUTION',T57,'TOTAL',
C     +    /,T20,'X',T23,F10.4,T37,F10.4,T54,F10.4,/,T20,'Y',T23,F10.4,
C     +    T37,F10.4,T54,F10.4,/,'    REQUIRED 2-FAM SEXT STRENGTHS:',
C    +    F14.6,F17.6,' M-2',10X,'SEXT. SUMS:',4F9.3)
C
C
      END
