C   05/09/85 504111748  MEMBER NAME  CHROTN   (MAR95.S)     FORTRAN
      SUBROUTINE CHROTN(E0,NTY,ICALL)
C
C
C=====ROUTINE TO RESET THE E3 ROTATORS FOR USE AT A DIFFERENT ENERGY.
C=====PARAMETERS FROM J.BUON'S 'ROTATOR' FILE .
C=====SEE END OF MEMBER 'ROTATOR' FOR DETAILS. ANGLES ARE READ IN RADS.
C
C     ICALL=0: READ FILE
C     ICALL=1: CHANGE ENERGY,FIND NEARBY DATA LINES.
C     ICALL=2: SET MAGNETS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CROT/H4(1500),V3(1500),H3(1500),V2(1500),H2(1500),V1(1500),
     +            H1(1500),EROT(1500),FC,NL,JL,JL1
      REAL*8 TNAME
      LOGICAL*1 ANAME2(4)
      REAL*4 ANAME(2)
      EQUIVALENCE (TNAME,ANAME(1))
      EQUIVALENCE (ANAME2(1),ANAME(2))
#include "cnlist.for"
#include "clatic.for"
#include "csol.for"
C
C
      REAL*4  BA05/'BA05'/
      REAL*4  BB01/'BB01'/
      REAL*4  BC05/'BC05'/
      REAL*4  BD05/'BD05'/
      REAL*4  BF06/'BF06'/
      REAL*4  BG06/'BG06'/
      LOGICAL*1 P/'P'/
      LOGICAL*1 M/'M'/
      LOGICAL*1 OST/'O'/
      REAL*4    EL/'OL  '/
      REAL*4    ER/'OR  '/
C
C
      ESET=E0
C     ESET=29.23
      IF(ICALL.GT.0)GO TO 15
C
C=====READ THE ROTATOR DATA FILE.
      NL=1
    9 READ(5,102)V3(NL),H3(NL),V2(NL),H2(NL),V1(NL),H1(NL)
      IF(V3(NL).GT.1.0D5)GO TO 12
      READ(5,102)H4(NL),DNU,EROT(NL)
      NL=NL+1
      GO TO 9
   12 CONTINUE
  102 FORMAT(6E13.7)
      RETURN
C
C=====USE THE ENERGY ESET TO SELECT THE CORRECT SETTINGS FOR EACH NEW
C=====ROTATOR ENERGY: ROTATOR SETTINGS COME IN ABOUT 20MEV STEPS.
C=====TRY INTERPOLATING TO SEE HOW WELL IT PERFORMS.
   15 IF(ICALL.EQ.2)GO TO 16
C     JL=(ESET+0.001D0-EROT(1))/0.02D0+1
      JL1=1
      DO 20 K=1,NL
      IF(ESET.GT.EROT(K))GO TO 20
      JL1=K
      GO TO 21
   20 CONTINUE
   22 WRITE(6,103)ESET
  103 FORMAT(' ',' ENERGY',F6.2,' OUT OF DESIGN RANGE OF ROTATOR')
      STOP
   21 CONTINUE
      JL=JL1-1
      IF(JL1.EQ.1)GO TO 22
      FC=(ESET-EROT(JL))/(EROT(JL1)-EROT(JL))
      RETURN
C
C
C
C=====FILL IN THE DIPOLE STRENGTHS.
   16 NN=NUNT(NTY)
C     NN=1
C     NUNT(NTY)=1
      TNAME=NAME(NTY)

C
C     IF(ANAME(1).EQ.BF06)WRITE(6,111)ANAME(1),ANAME(2),ANAME2
  111 FORMAT(' ',A4,2X,A4,2X,A1,2X,A1,2X,A1,2X,A1)
      IF(ANAME(1).EQ.BF06)XX(NTY)= ((V1(JL1)-V1(JL))*FC+V1(JL))/NN
      IF(ANAME(1).EQ.BF06.AND.ANAME2(4).EQ.M)XX(NTY)=-XX(NTY)
      IF(ANAME(1).EQ.BG06)XX(NTY)= ((V3(JL1)-V3(JL))*FC+V3(JL))/NN
      IF(ANAME(1).EQ.BG06.AND.ANAME2(4).EQ.M)XX(NTY)=-XX(NTY)
      IF(ANAME(1).EQ.BA05)XX(NTY)= ((H1(JL1)-H1(JL))*FC+H1(JL))/NN
C     IF(ANAME(1).EQ.BA05)XX(NTY)= XX(NTY)-0.0034D0/NN
      IF(ANAME(1).EQ.BB01)XX(NTY)= ((H2(JL1)-H2(JL))*FC+H2(JL))/NN/2.D0
      IF(ANAME(1).EQ.BC05)XX(NTY)= ((H3(JL1)-H3(JL))*FC+H3(JL))/NN/2.D0
      IF(ANAME(1).EQ.BD05)XX(NTY)= ((H4(JL1)-H4(JL))*FC+H4(JL))/NN/2.D0
C=====TURN OF EAST VERTICAL BENDS.
C     IF(ANAME(1).EQ.BF06.AND.(ANAME2(1).EQ.EL.OR.ANAME2(1).EQ.OST))
C    +                                                XX(NTY)=0.D0
C     IF(ANAME(1).EQ.BG06.AND.(ANAME2(1).EQ.EL.OR.ANAME2(1).EQ.OST))
C    +                                                XX(NTY)=0.D0
C
C     WRITE(6,112)E0,ANAME(1),NAME(NTY),XX(NTY)
C 112 FORMAT(' ',F14.8,2X,A4,2X,A8,2X,F14.8)
C
C
C
      RETURN
      END
