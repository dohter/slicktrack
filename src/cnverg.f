C   18/12/95 601231727  MEMBER NAME  CNVERG   (SEPT95.S)    FVS
C
      PROGRAM CNVERG
C
C
C
C  MAIN ROUTINE TO CHECK FOR VANISHING   OF SUMS OF OSCILLATING TERMS.
C  -------------------------------------------------------------------
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      PI=3.1415926535897932D0
C     TUNE=SQRT(7.0D0)
C     TUNE=3.D0
      TUNE=2.D0/3.D0
      TWOPI=2.D0*PI
C
C
      CSUM=0.D0
      SSUM=0.D0
      DO 10 N=1,200000
      C=COS(TWOPI*TUNE*N)
      S=SIN(TWOPI*TUNE*N)
      CSUM=CSUM+C
      SSUM=SSUM+S
      CSUMAV=CSUM/N
      SSUMAV=SSUM/N
      IF(MOD(N,100).EQ.0)WRITE(6,1000)N,C,S,CSUM,SSUM,CSUMAV,SSUMAV
 1000 FORMAT(' ','N,C,S,CSUM,SSUM,CSUMAV,SSUMAV: ',I10,6F10.5)
   10 CONTINUE
C
C
      STOP
      END