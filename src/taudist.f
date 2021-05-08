C   19/05/82 603081607  MEMBER NAME  SLICK    (SEPT95.S)    FORTRAN
C
C 
       PROGRAM TAUDIST  
 
C
C
C       MAIN ROUTINE FOR CHECKING DISTRIBUTIONS OF TAUD
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*80 TAUDISTDAT,OUTPUT
      NAMELIST/FILES/TAUDISTDAT,OUTPUT
      DIMENSION KGBIN1(-1000:1000),KGBIN2(-1000:1000)
C
C
C
C
      PI=3.1415926535897932D0
      OUTPUT       = 'taudist.out'
      TAUDISTDAT   = 'taudist_dat.out'
      READ(5,FILES)                                                        
      OPEN(53,FILE=OUTPUT,      STATUS='UNKNOWN')
      OPEN(51,FILE=TAUDISTDAT,  STATUS='UNKNOWN')


      CALL G05CBF(4)

      NG     = 0
      SUMG   = 0.D0
      SUMG2  = 0.D0 
      OFFSET = 1.D0
      WIDTH  = 0.1D0
      GBIN1  = WIDTH*0.1D0
      GBIN2  = 1.D0/WIDTH *0.001D0

      OFFSET = 0.1D0
      WIDTH  = 1.0D0
      GBIN1  = OFFSET*0.5D0
      GBIN2  = 1/WIDTH *0.1D0

      KGBIN1 = 0
      KGBIN2 = 0

 1001 GAUS1 = G05DDF(OFFSET, WIDTH)
      SUMG2 = SUMG2 + GAUS1*GAUS1
      SUMG  = SUMG  + GAUS1
      NG    = NG + 1
      NGAUS1 = GAUS1/GBIN1
      NGAUS2 = 1.D0/GAUS1/GBIN2
      IF(IABS(NGAUS1).LE.1000)KGBIN1(NGAUS1) = KGBIN1(NGAUS1) + 1
      IF(IABS(NGAUS2).LE.1000)KGBIN2(NGAUS2) = KGBIN2(NGAUS2) + 1
      IF(NG.LT.100000)GO TO 1001 

C======Write out statistics on the Gaussians.
      IF(NG.GT.0)THEN
      SUMG2 = SUMG2/NG
      SUMG  = SUMG/NG
      SUMG2 = DSQRT(SUMG2)
      WRITE(53,'(1X,A,I10,2(1X,F15.6))')
     + 'Check the Gaussian: calls, mean, r.m.s. ',NG,SUMG,SUMG2
      ENDIF
      DO 556 JG = -1000,1000
  556 WRITE(51,'(1X,I10,1X,I10,1X,I10)')JG,KGBIN1(JG),KGBIN2(JG)

      CLOSE(51) 
      CLOSE(53)
C
C
C
C
C
      STOP
      END