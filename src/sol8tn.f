C   15/11/79 402121404  MEMBER NAME  SOL8TN   (S)           FORTRAN
      SUBROUTINE SOL8TN(I,X,Y,S,ZW,A,B,NSOL)
C
C
C     VERSION TO BE USED IF IT IS REQUIRED TO SPLIT THE SOLENOID INTO
C     MANY SLICES FOR CLOSED ORBIT TESTS.
C     IF SINGLE SLICE ANALYTIC METHOD IS NEEDED USE SOL8AN
C
C=====BARBER AUGUST 1982 REPAIR THE ERROR IN T4 (SOLENOID MID SECTION)
C     T4(7,6)= T1(3,3)*X*(1+A)
C     T4(8,6)=-T1(3,2)*X*(1+A)
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DATA ICALL/0/
      DIMENSION ZW(3,3),A(3,3),B(8,8),T1(3,3),T2(3,3),T3(3,3),T4(8,8),
     +T5(8,8),T6(4,4)
C
C
#include "cloorb.for"
C
C
      XX1=X*(1.+S)
C     IF(ICALL.EQ.1)GO TO 103
C     WRITE(6,104)
C 104 FORMAT(' ','SOL88 INPUT===================')
C     WRITE(6,105)I,X,Y,S,DX(I),DY(I),DXP(I),DYP(I),DEL(I)
C 105 FORMAT(' ',I6,8E12.5,/'  ZW===================')
C     WRITE(6,106)ZW
C 106 FORMAT(' ',3F12.5)
C 103 CONTINUE
      CALL RSPIN(XX1*DX(I)*.5/Y,XX1*DY(I)*.5/Y,0.D0,T1)
      CALL JAM333(T2,T1,ZW)
      CALL AVER(ZW,T2,3,T3)
      CALL UNIT(8,B)
      B(2,3)=.5*X/Y
      B(4,1)=-B(2,3)
      B(7,1)=.5*XX1*T3(1,3)/Y
      B(7,3)=.5*XX1*T3(2,3)/Y
      B(8,1)=-.5*XX1*T3(1,2)/Y
      B(8,3)=-.5*XX1*T3(2,2)/Y
      IF(NSOL.EQ.0)GO TO 15
      B(7,1)=0.D0
      B(7,3)=0.D0
      B(8,1)=0.D0
      B(8,3)=0.D0
   15 CONTINUE
      CALL SOLXYP(I,X/Y,DXPA,DYPA)
      XDELA=X*(DEL(I)+DEL(I+1))*.5
      CALL RSPIN(X*S*DXPA,X*S*DYPA,-X+XDELA,T1)
C     CALL RSPIN(X*S*DXPA,X*S*DYPA,-X/2.+XDELA,T1)
      CALL JAM333(T3,T1,T2)
      CALL AVER(T2,T3,3,T1)
      CALL UNIT(8,T4)
      CALL SOL6(X,Y,T6)
      DO 10 III=1,4
      DO 10 JJJ=1,4
 10   T4(III,JJJ)=T6(III,JJJ)
      DS=DSIN(X)
      DC1=1.D0-DCOS(X)
C     IF(ICALL.EQ.0)WRITE(6,108)
C 108 FORMAT(' ','T1 ARRAY================')
C     IF(ICALL.EQ.0)WRITE(6,107)T1
C 107 FORMAT(' ',3F12.5)
      T4(7,2)=S*T1(1,3)*DS-S*T1(2,3)*DC1
      T4(7,4)=S*T1(2,3)*DS+S*T1(1,3)*DC1
      T4(7,6)=T1(3,3)*(1.+0.0011596)*X
      T4(8,2)=-S*T1(1,2)*DS+S*T1(2,2)*DC1
      T4(8,4)=-S*T1(2,2)*DS-S*T1(1,2)*DC1
      T4(8,6)=-T1(3,2)*(1.+0.0011596)*X
      IF(NSOL.EQ.0)GO TO 16
      T4(7,2)=0.D0
      T4(7,4)=0.D0
      T4(7,6)=0.D0
      T4(8,2)=0.D0
      T4(8,4)=0.D0
      T4(8,6)=0.D0
   16 CONTINUE
C     WRITE(6,109)
C 109 FORMAT(' ','T4---THE G(2X6) MATRIX')
C     WRITE(6,110)T4(7,2),T4(7,4),T4(7,6),T4(8,2),T4(8,4),T4(8,6)
C 110 FORMAT(' ',3E20.6)
      CALL JAM888(T5,T4,B)
      CALL RSPIN(-XX1*DX(I+1)*.5/Y,-XX1*DY(I+1)*.5/Y,0.D0,T1)
      CALL JAM333(A,T1,T3)
      CALL AVER(T3,A,3,T1)
      CALL UNIT(8,T4)
      T4(2,3)=-.5*X/Y
      T4(4,1)=-T4(2,3)
      T4(7,1)=-.5*T1(1,3)*XX1/Y
      T4(7,3)=-.5*XX1*T1(2,3)/Y
      T4(8,1)=.5*XX1*T1(1,2)/Y
      T4(8,3)=.5*XX1*T1(2,2)/Y
      IF(NSOL.EQ.0)GO TO 17
      T4(7,1)=0.D0
      T4(7,3)=0.D0
      T4(8,1)=0.D0
      T4(8,3)=0.D0
   17 CONTINUE
      CALL JAM888(B,T4,T5)
C     IF(ICALL.EQ.1)RETURN
C     WRITE(6,101)
C 101 FORMAT(' ','SOL88 OUTPUT=========================')
C     WRITE(6,102)((B(LL,MM),MM=1,8),LL=1,8)
C 102 FORMAT(' ',8E12.5)
      ICALL=1
      RETURN
      END
