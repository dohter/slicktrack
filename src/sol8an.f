C   15/11/79 510241814  MEMBER NAME  SOL8AN   (S)           FORTRAN
      SUBROUTINE SOL8AN(I,X,Y,SS,ZW,A,B,NSOL)
C

C     1/6/2007: changed from .5 to 0.5D0 in many places. Old version in save_sol8an.f
C
C
C
C
C=====NEED TO PUT IN EXACT G IN RSPIN AND HERE.
C     USE SOL8AN IN ROTS.
C
C     'THICK' VERSION IN WHICH INTEGRALS ARE DONE ANALYTICALLY.
C     NEED CARE IF HAVE CLOSED ORBIT DISTORTIONS.
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
      DIMENSION RSOL(3)
C
C
      INCLUDE "cloorb.for"
C
C=====MAKE A SPIN TRANSPARENT THIN ROTATOR.
C     IF(1.EQ.1)GO TO 177
      NSOLM6=NSOL-6
      IF(NSOLM6.GT.0)GO TO 177
C
      S=SS
      XX1=X*(1.D0+S)  
      XX2=X*(1.D0+S)                          *1.D0  !Switch off C.O. effects Sept 2003.
      CALL RSPIN(XX2*DX(I)*0.5D0/Y,XX2*DY(I)*0.5D0/Y,0.D0,T1)
      CALL JAM333(T2,T1,ZW)
      CALL AVER(ZW,T2,3,T3)
      CALL UNIT(8,B)
      B(2,3)=0.5D0*X/Y
      B(4,1)=-B(2,3)
      B(7,1)=0.5D0*XX1*T3(1,3)/Y
      B(7,3)=0.5D0*XX1*T3(2,3)/Y
      B(8,1)=-0.5D0*XX1*T3(1,2)/Y
      B(8,3)=-0.5D0*XX1*T3(2,2)/Y
      IF(NSOL.EQ.0)GO TO 15
      B(7,1)=0.D0
      B(7,3)=0.D0
      B(8,1)=0.D0
      B(8,3)=0.D0
   15 CONTINUE
      CALL SOLXYP(I,X/Y,DXPA,DYPA)
      XDELA=X*(DEL(I)+DEL(I+1))*0.5D0           *1.D0 ! Switch of C.O. effects Sept 2003
      DXPA = DXPA                               *1.D0
      DYPA = DYPA                               *1.D0
C      CALL RSPIN(X*S*DXPA,X*S*DYPA,-X+XDELA,T1)
      CALL RSPIN(X*S*DXPA,X*S*DYPA,-X*(1.D0 +0.001159652D0)+XDELA,T1)    !  corrected (1+a) on Oct 6 2020
      CALL JAM333(T3,T1,T2)
      CALL AVER(T2,T3,3,T1)
      CALL UNIT(8,T4)
      CALL SOL6(X,Y,T6)
      DO 10 III=1,4
      DO 10 JJJ=1,4
 10   T4(III,JJJ)=T6(III,JJJ)
C     DS=DSIN(X)
C     DC1=1.D0-DCOS(X)
C     T4(7,2)=S*T1(1,3)*DS-S*T1(2,3)*DC1
C     T4(7,4)=S*T1(2,3)*DS+S*T1(1,3)*DC1
C     T4(8,2)=-S*T1(1,2)*DS+S*T1(2,2)*DC1
C     T4(8,4)=-S*T1(2,2)*DS-S*T1(1,2)*DC1
C     T4(7,6)= T1(3,3)*(1.+0.0011596)*X
C     T4(8,6)=-T1(3,2)*(1.+0.0011596)*X
C
      T4(7,6)= T2(3,3)*(1.D0 +0.0011596)*X
      T4(8,6)=-T2(3,2)*(1.D0 +0.0011596)*X
      T4(7,2)= S*X*T2(1,3)
      T4(7,4)= S*X*T2(2,3)
      T4(8,2)=-S*X*T2(1,2)
      T4(8,4)=-S*X*T2(2,2)
      IF(NSOL.EQ.0)GO TO 16
      T4(7,2)=0.D0
      T4(7,4)=0.D0
      T4(7,6)=0.D0
      T4(8,2)=0.D0
      T4(8,4)=0.D0
      T4(8,6)=0.D0
   16 CONTINUE
      CALL JAM888(T5,T4,B)
      CALL RSPIN(-XX2*DX(I+1)*0.5D0/Y,-XX2*DY(I+1)*0.5D0/Y,0.D0,T1)
      CALL JAM333(A,T1,T3)
      CALL AVER(T3,A,3,T1)
      CALL UNIT(8,T4)
      T4(2,3)=-0.5D0*X/Y
      T4(4,1)=-T4(2,3)
      T4(7,1)=-0.5D0*T1(1,3)*XX1/Y
      T4(7,3)=-0.5D0*XX1*T1(2,3)/Y
      T4(8,1)=0.5D0*XX1*T1(1,2)/Y
      T4(8,3)=0.5D0*XX1*T1(2,2)/Y
      IF(NSOL.EQ.0)GO TO 17
      T4(7,1)=0.D0
      T4(7,3)=0.D0
      T4(8,1)=0.D0
      T4(8,3)=0.D0
   17 CONTINUE
      CALL JAM888(B,T4,T5)
      ICALL=1
C      DO 273 III=1,8
C  273 WRITE(53,925)(B(III,JJJ),JJJ=1,8)
C  925 FORMAT(' ','Solenoid  ',8D16.8)
      RETURN
C
C
C=====MAKE ZERO LENGTH TRANSPARENT ROTATOR.
  177 RSOL=0.D0
      RSOL(NSOLM6)=-X
      CALL RSPIN(RSOL(1),RSOL(2),RSOL(3),T1)
      CALL JAM333(A,T1,ZW)
      CALL UNIT(8,B)
      RETURN
      END
