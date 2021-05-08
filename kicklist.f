C
C
       PROGRAM KICKLIST 
 
C
C
C                 MAIN ROUTINE FOR RENAMING KICKERS ETC BEFORE RUNNING SLICK
C                 ----------------------------------------------------------
C
C                    This is only OK for FORMAT 4.
C
C   This reads a standard SLICK input file, identifies the kickers in the geometry list and
C   then renames them with increasing numerical identifiers. A type list for the kickers
C   and a new geometry list is written to output.
C   
C   Only a shortened, modified version of LATTIN is needed. This is integral to this
C   main routine. The input is read in the usual way but ignored. The focus is on the
C   geometry list.
C
C   NTYPE=(NUMBER OF DIFFERENT MAGNET-TYPES) <= 1000
C   WARNING ==> WHEN CHANGING NO. OF TYPES, TRANS1 MUST BE CHANGED
C   BEAM-LINE MUST CONTAIN ALL ELEMENTS IN STORAGE RING.
C   2 <= (TOTAL NUMBER OF ELEMENTS IN THE BEAM-LINE) <= (NO. IN CLOORB)
C
C
C   EACH MAGNET-TYPE IS SPECIFIED BY (ID,NAME,XX,YY ------).
C   ID=  ( 1 , 2 , 3 ,    4   , 5 ,   6    ,   7    ,  8 , 9 ,  10)
C   MEANS(DRF,BBX,QUA,SKEW-QUA,CAV,HOR-KICK,VER-KICK,SEXT,BBY,SOLENOID)
C   ID=  (11         ,12           ,13           )
C   MEANS(ROTAT.QUAD  FOC.HOR.BEND  FOC.VERT.BEND)
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*80 RING,OUTPUT
      NAMELIST/FILES/RING,OUTPUT
C======Ring     = input file.
C======Output   = printout
C
C
      INCLUDE "clatic.for"
      INCLUDE "csol.for"
  
      CHARACTER*4 KAME,JAME,NAM1,NAM2
      CHARACTER*1 LAME(4)
      DIMENSION KD(3000),KAME(3000),XXK(3000),X2K(3000),YYK(3000),             
     +        NUNTK(3000),TWISTK(3000),NSOLK(3000),MAME(3000)
      EQUIVALENCE (JAME,LAME(1))
C
C
C
      DIMENSION PSUP(100,4)
      DATA PSUP/400*1.D40/
      CHARACTER *4 TEXT(20)
      CHARACTER STRING *80

      CHARACTER *8 NAM(6)
      DIMENSION IPOS(6)
C
C
      PI=3.1415926535897932D0
      RING=    'he7basml.inp'
      OUTPUT=  'kicklist.out' 

      READ(5,FILES)                                                
      OPEN(52,FILE=RING,    STATUS='OLD') 
      OPEN(53,FILE=OUTPUT,  STATUS='UNKNOWN')
C

C
C    FORMAT=4:  4 MAGNET POSITIONS/LINE & 1/10 MM. UNITS
C    FORMAT=5:  5 MAGNET POSITIONS/LINE & 1/10 MM. UNITS
C    FORMAT=6:  6 MAGNET POSITIONS/LINE &    1 MM. UNITS
C

C
      CALL DATIK
C
C
C
C
C=====READ THE DATA FILE HEADING
      READ (52,901)TEXT
      WRITE(53,901)TEXT
      READ (52,901)TEXT
      WRITE(53,901)TEXT
      READ (52,901)TEXT
      WRITE(53,901)TEXT
C
C=====CHECK FORMAT AND RANDOM NUMBER GENERATOR SEEDS.
      READ(52,801)IFORM
      READ(52,805)KSEED
      READ(52,805)NINS
      IF(IFORM.EQ.4.OR.IFORM.EQ.5.OR.IFORM.EQ.6)GO TO 804
      WRITE(53,803)
      STOP


  804 CONTINUE


      SCL=0.001D+0
      SEP=-0.001001D+0
      SEP=-0.0001001D+0

      IF(IFORM.EQ.5.OR.IFORM.EQ.4)SCL=0.0001D+0
      IF(IFORM.EQ.5.OR.IFORM.EQ.4)SEP=-0.001101D+0
C
C
C
C
      NTY=1
C=====LOOP TO READ THE MAGNET PARAMETERS & CHECK FOR COMMENTS.
      WRITE(53,952)
   26 READ(52,'(A)',END=999)STRING
      IF(STRING(:1).NE.' ')THEN
      WRITE(53,'(1X,A,A)')STRING, '  Comment'
      GO TO 26
      ELSE
      NAME(NTY)=' '
      IF(IFORM.EQ.4)READ(STRING,99)KOM,ID(NTY),NAME(NTY),XX(NTY),
     +                   X2(NTY),YY(NTY),NUNT(NTY),TWIST(NTY),NSOL(NTY)
      IF(IFORM.NE.4)READ(STRING,89)KOM,ID(NTY),NAME(NTY),XX(NTY),
     +                   X2(NTY),YY(NTY),NUNT(NTY),TWIST(NTY),NSOL(NTY)
      NTWIST(NTY)=TWIST(NTY)
      WRITE(53,89)KOM,ID(NTY),NAME(NTY),XX(NTY),X2(NTY),YY(NTY),             
     +        NUNT(NTY),TWIST(NTY),NSOL(NTY)
      ENDIF


      IF(NAME(NTY)(1:3).EQ.'END')GO TO 3

CCC   IF(ID(NTY).NE.15)GO TO 192
CC      if(nunt(nty).eq.1)go to 192                
CCC   XX(NTY)=XX(NTY)/16.D+0                                                    
CCC   X2(NTY)=X2(NTY)/16.D+0
CCC   NUNT(NTY)=16                                                    
CC      XX(NTY)=XX(NTY)/2.D+0                                                    
CC      X2(NTY)=X2(NTY)/2.D+0
CC      NUNT(NTY)=2          
  192 CONTINUE     

C
C
C=====Check for legal type.
      J=ID(NTY)
      GO TO(105,105,105,105,105,105,105,105,105,105,105,105,105,105,
C             D   B   Q  QK  RF  HK  VK  SX  BY SOL  RQ FBH FBV ROT
     +105,105),J
C     CFH CFV
      WRITE(53,96) J
      STOP
  999 WRITE (53,950)
      STOP
  105 CONTINUE

C
C=====SET UP A FIRST DRIFT OF ZERO LENGTH----THIS IS NOT THE I.P.!
C=====THE I.P. WILL BE THE 3RD ELEMENT!---AFTER THIS DUMMY ELEMENT & THE
C=====FOLLOWING ZERO LENGTH DRIFT.
 1044 NTY=NTY+1
      IF(NTY.LE.3000)GO TO 26
      WRITE(53,94)
      STOP
    3 CONTINUE
      ND1=NTY
      NAME(NTY)='DRIFT'
      ID(NTY)=1
      XX(NTY)=0.D0
      YY(NTY)=0.D0
      NT1=NTY-1

      NELEM=1
      ITYPE(1)=NTY
      IREAD=0
      NEXTRA = 0
      NVC = 0
      NHC = 0
      NVQ = 0
      NHQ = 0
      NVD = 0
      NRQ = 0
      NCQ = 0
C
C
C
C
C=====READ THE STRUCTURE OF THE RING INTO STORAGE.
C=====CAVITIES & KICKERS GET ZERO LENGTH (FOR CAVITIES YY=0 ON INPUT)
C=====EVERYTHING ELSE GETS ITS ACTUAL LENGTH.
C=====WRITE OUT THE STRUCTURE WITH NEW KICKER NAMES.
      WRITE(53,933)
    4 CONTINUE
      IF(IFORM.EQ.6)READ(52,92  ,END=30) (NAM(KK),IPOS(KK),KK=1,6)
      IF(IFORM.EQ.5)READ(52,922 ,END=30) (NAM(KK),IPOS(KK),KK=1,5)
      IF(IFORM.EQ.4)READ(52,9222,END=30) (NAM(KK),IPOS(KK),KK=1,4)
      IF(IREAD.EQ.0)TOTL=IPOS(1)*SCL
      IREAD=1
C      IF(IFORM.EQ.6)WRITE(53,95)  (NAM(KK),IPOS(KK),KK=1,6)
C      IF(IFORM.EQ.5)WRITE(53,955) (NAM(KK),IPOS(KK),KK=1,5)
C      IF(IFORM.EQ.4)WRITE(53,9555)(NAM(KK),IPOS(KK),KK=1,4)
      DO 13 KK=1,IFORM
      IF(NAM(KK).EQ.'        ')GO TO 13
      IF(NAM(KK).EQ.'END')GO TO 30
      IREM=0
C=====SEARCH FOR NAME IN TYPE LIST.
      DO 5 IN=1,NT1
      IF(NAM(KK).NE.NAME(IN))GO TO 5
      IREM=ID(IN)
      YYREM=YY(IN)
      GO TO 6
    5 CONTINUE
      WRITE(53,91)NAM(KK)
      STOP

    6 CONTINUE   

C=======Set up new kicker names.
      IF(NAM(KK).NE.'VC      '.AND.
     +   NAM(KK).NE.'HC      '.AND.
     +   NAM(KK).NE.'VQ      '.AND.
     +   NAM(KK).NE.'HQ      '.AND.
     +   NAM(KK).NE.'VD      '.AND.
     +   NAM(KK).NE.'RQ      '.AND.
     +   NAM(KK).NE.'CQ      ')GO TO 66
       
      NEXTRA = NEXTRA + 1 
      MAME(NEXTRA) = 0
 
      IF(NAM(KK).EQ.'VC      ')THEN
      KD(NEXTRA)=7
      NVC = NVC + 1
      LAME(1)='V'
      LAME(2)='C'
      LAME(3)='0'
      LAME(4)='0'
      MAME(NEXTRA)=NVC 
      KAME(NEXTRA)=JAME
      NAM1=JAME       
      WRITE(NAM2,'(I4)')NVC
      NAM(KK)=NAM1//NAM2
      ENDIF

      IF(NAM(KK).EQ.'HC      ')THEN
      KD(NEXTRA)=6
      NHC = NHC + 1
      LAME(1)='H'
      LAME(2)='C'
      LAME(3)='0'
      LAME(4)='0'
      MAME(NEXTRA)=NHC 
      KAME(NEXTRA)=JAME
      NAM1=JAME       
      WRITE(NAM2,'(I4)')NHC
      NAM(KK)=NAM1//NAM2
      ENDIF

      IF(NAM(KK).EQ.'VQ      ')THEN
      KD(NEXTRA)=7
      NVQ = NVQ + 1
      LAME(1)='V'
      LAME(2)='Q'
      LAME(3)='0'
      LAME(4)='0'
      MAME(NEXTRA)=NVQ 
      KAME(NEXTRA)=JAME
      NAM1=JAME       
      WRITE(NAM2,'(I4)')NVQ
      NAM(KK)=NAM1//NAM2
      ENDIF

      IF(NAM(KK).EQ.'HQ      ')THEN
      KD(NEXTRA)=6
      NHQ = NHQ + 1
      LAME(1)='H'
      LAME(2)='Q'
      LAME(3)='0'
      LAME(4)='0'
      MAME(NEXTRA)=NHQ
      KAME(NEXTRA)=JAME
      NAM1=JAME       
      WRITE(NAM2,'(I4)')NHQ
      NAM(KK)=NAM1//NAM2
      ENDIF

      IF(NAM(KK).EQ.'VD      ')THEN
      KD(NEXTRA)=7
      NVD = NVD + 1
      LAME(1)='V'
      LAME(2)='D'
      LAME(3)='0'
      LAME(4)='0'
      MAME(NEXTRA)=NVD
      KAME(NEXTRA)=JAME
      NAM1=JAME       
      WRITE(NAM2,'(I4)')NVD
      NAM(KK)=NAM1//NAM2
      ENDIF

      IF(NAM(KK).EQ.'RQ      ')THEN
      KD(NEXTRA)=4
      NRQ = NRQ + 1
      LAME(1)='R'
      LAME(2)='Q'
      LAME(3)='0'
      LAME(4)='0'
      MAME(NEXTRA)=NRQ
      KAME(NEXTRA)=JAME
      NAM1=JAME       
      WRITE(NAM2,'(I4)')NRQ
      NAM(KK)=NAM1//NAM2
      ENDIF

      IF(NAM(KK).EQ.'CQ      ')THEN
      KD(NEXTRA)=3
      NCQ = NCQ + 1
      LAME(1)='C'
      LAME(2)='Q'
      LAME(3)='0'
      LAME(4)='0'
      MAME(NEXTRA)=NCQ
      KAME(NEXTRA)=JAME
      NAM1=JAME       
      WRITE(NAM2,'(I4)')NCQ
      NAM(KK)=NAM1//NAM2
      ENDIF

      XXK(NEXTRA) = 0.D0                            ! Set zero, fill 
      X2K(NEXTRA) = 0.D0                            ! during running.
      YYK(NEXTRA) = YYREM
      NUNTK(NEXTRA) = 1
      TWISTK(NEXTRA)=0.D0
      NSOLK(NEXTRA) =0
      IF(NAM(KK).EQ.'VC      ')YYK(NEXTRA) = 0.1D0  ! A first shot.
      IF(NAM(KK).EQ.'HC      ')YYK(NEXTRA) = 0.1D0  !Prevents problems 
      IF(NAM(KK).EQ.'VQ      ')YYK(NEXTRA) = 0.1D0  ! with radiation.
      IF(NAM(KK).EQ.'HQ      ')YYK(NEXTRA) = 0.1D0                   
      IF(NAM(KK).EQ.'VD      ')YYK(NEXTRA) = 0.1D0   
      IF(NAM(KK).EQ.'RQ      ')YYK(NEXTRA) = 0.0D0   
      IF(NAM(KK).EQ.'CQ      ')YYK(NEXTRA) = 0.0D0      


   66 CONTINUE  

      NU=NUNT(IN)
      NUS=NU
      PLT=IPOS(KK)*SCL
      PL=PLT-YY(IN)*NU*0.5D0-TOTL

      IF(IREM.EQ.10.OR.IREM.EQ.14)PL=PLT-TOTL
C
C=====HOR. & VERT. KICKERS WHICH ARE GIVEN ZERO LENGTH IN THE LATTICE.
      IF(IREM.EQ. 6.OR.IREM.EQ. 7)PL=PLT-TOTL

C
C
C=====CONSTRUCT THE DRIFT SPACE UP TO THE NEXT THICK LENS.
   11 DO 7 I=ND1,NTY
      IF(DABS(PL*1000.D0-XX(I)*1000.D0).LT.0.1000001D0)GO TO 10
    7 CONTINUE
      NTY=NTY+1
      IF(NTY.LE.3000)GO TO 8
      WRITE(53,94)
      STOP
    8 CONTINUE
      SNAME(NTY)=NAME(IN)
      NAME(NTY)= 'DRIFT'
      ID(NTY)=1
C=====TAKE CARE OF DIGITISATION ERRORS ON INPUT FILE--LIMIT GIVEN BY SEP
      IF(PL.GE.SEP.AND.PL.LT.0.0D0)PL=0.D0
C     IF(PL.GE.-0.0006.AND.PL.LT.0.0D0)PL=0.D0
      XX(NTY)=PL
      YY(NTY)=PL
      I=NTY
   10 NELEM=NELEM+1
      IF(NELEM.LE.50000)GO TO 12
      WRITE(53,907)TOTL,NELEM
      STOP
   12 ITYPE(NELEM)=I

      DO 212 INU=1,NU
      NELEM=NELEM+1
  212 CONTINUE  


      TOTL=YY(IN)*NUS*0.5D0+PLT
      IF(IREM.EQ.10.OR.IREM.EQ.14)TOTL=YY(IN)*NUS+PLT
      IF(IREM.EQ. 6.OR.IREM.EQ. 7)TOTL=PLT
   13 CONTINUE

      IF(IFORM.EQ.6)WRITE(53,905)  (NAM(KK),IPOS(KK),KK=1,6)
      IF(IFORM.EQ.5)WRITE(53,9055) (NAM(KK),IPOS(KK),KK=1,5)
      IF(IFORM.EQ.4)WRITE(53,90555)(NAM(KK),IPOS(KK),KK=1,4)

      GO TO 4


   30 CONTINUE
      WRITE(53,902)NELEM,NTY
C
C
C
C
C======PRINT OUT TYPE LIST.
      WRITE(53,952)
      DO 14 I=1,NTY
      WRITE(53,89)KOM,ID(I),NAME(I),XX(I),X2(I),YY(I),             
     +        NUNT(I),TWIST(I),NSOL(I)
C      WRITE(53,97)I,ID(I),NAME(I),XX(I),X2(I),YY(I),NUNT(I),TWIST(I),
C     +SNAME(I) ,NSOL(I),NTWIST(I)
      IF(ID(I).EQ.1.AND.XX(I).LT.0.D0)WRITE(53,98)
      IF(ID(I).EQ.1.AND.XX(I).LT.0.D0)STOP
   14 CONTINUE
C

C======PRINT OUT KICKET LIST WITH NEW NAMES.
      WRITE(53,933)
      WRITE(53,952)
      DO 144 I=1,NEXTRA
      WRITE(53,899)KOM,KD(I),KAME(I),MAME(I),XXK(I),X2K(I),YYK(I),             
     +        NUNTK(I),TWISTK(I),NSOLK(I)
 144  CONTINUE




C
C
C
   60 FORMAT(I6,A8,I4,4(F14.5))
   80 FORMAT(80A1)
   91 FORMAT(' NAME ',A8,' NOT FOUND IN THE TYPE LIST')
   92 FORMAT(6(A4,I7,1X))
  922 FORMAT(5(A4,1X,I8,1X))
 9222 FORMAT(4(A8,1X,I8,1X))
C  95 FORMAT(5(4A1,1X,I8,1X))
   93 FORMAT('1',' THE LATTICE AS READ FROM THE INPUT DATA---')
  933 FORMAT('1',' THE LATTICE WITH NEW KICKER NAMES---') 
 9333 FORMAT('1',' NEW KICKER TYPE LIST---') 
   94 FORMAT(' SIZE OF TYPE LIST EXCEEDED !!!')
   95 FORMAT(6(1X,A4,I7,1X))
  955 FORMAT(5(1X,A4,1X,I8,1X))
 9555 FORMAT(4(1X,A8,1X,I8,1X))
  905 FORMAT(6(1X,A4,I7,1X))
 9055 FORMAT(5(1X,A4,1X,I8,1X))
90555 FORMAT(4(A8,1X,I8,1X))   
   96 FORMAT('  FALSE MAGNET TYPE  ',I6)
   97 FORMAT(' ',2I5,1X,A8,3F12.8,I5,2X,F11.6,4X,A8,I4,I10)
   98 FORMAT(' ','A NEGATIVE DRIFT LENGTH HAS BEEN FOUND--SO STOP')
   99 FORMAT(  A1,I4,1X,A8,3E12.8,I5   ,F11.6,I5)
  199 FORMAT(8X,I2,1A1,3X)
   89 FORMAT(  A1,I4,1X,A8,3F12.8,I5   ,F11.6,I5)
  899 FORMAT(  A1,I4,1X,A4,I4,3F12.8,I5   ,F11.6,I5)
  801 FORMAT(9X,I1)
  803 FORMAT(' ',' STOP--BAD LATTICE FORMAT ')
  805 FORMAT(8X,I2)
  901 FORMAT(20A4)
  902 FORMAT(' LATTICE FILLED UP TO ',I7,', TYPE LIST FILLED UP TO ',I7)
C
C
  907 FORMAT(' SIZE OF LATTICE LIST EXCEEDED !!!',/' LTOT =',F11.4,
     +' ELEMENT=',I6)
  950 FORMAT(' SOMETHING IS WRONG WITH THE INPUT DATA-FILE ')
  951 FORMAT(1X,80A1)
  952 FORMAT('0',//,'  TYPE NAME      STRENGTH1  STRENGTH2   LENGTH',
     +  '   NO.SLICES  SFLAG GMATRIX FLAG')

C
C
      CLOSE(53)  
      CLOSE(52) 
C
C
      STOP
C
C
      END
