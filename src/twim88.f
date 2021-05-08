C   05/06/84 406052129  MEMBER NAME  TWIM88   (S)           FORTRAN
      SUBROUTINE TWIM88(ID,ITY,XX,YYY,DX,DY,D1,D2,D3,D4,S,ZW,A)
C
C
C
C
C=====ROUTINE FOR CONSTRUCTING 8X8 TRANSPORT MATRICES USING THE
C     TWISS PARAMETER BASIS INSTEAD OF THE (X,X',Z,Z',DS,DE/E) BASIS.
C     ACCORDING TO THE FORMALISM IN DESY M 84/?
C     THE ROUTINE IS CALLED IN THE SAME WAY AS THE CONVENTIONAL MX88.
C     AS WITH MX88 ALL MATRIX ELEMENTS ARE GENERATED AFRESH FROM THE
C     STORED CLOSED ORBIT AND DISPERSION VALUES IN ORDER TO SAVE TIME.
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 DX,DY
      DIMENSION ZW(3,3),A(8,8)
C
#include "csol.for"
C
C
      CALL VZERO(A,128)
      A(1,1)=1.D0
      A(2,2)=1.D0
      A(3,3)=1.D0
      A(4,4)=1.D0
      A(5,5)=1.D0
      A(6,6)=1.D0
      A(7,7)=1.D0
      A(8,8)=1.D0
C
C
      NNSOL=NSOL(ITY)
      GOTO(1,2,3,4,5,6,7,8,9,10),ID
    1 RETURN
C
C
C=====QUADRUPOLE
    3 CONTINUE
      TEMP=(1.+S)*XX
      A(2,1)=-XX
      A(4,3)=+XX
      A(2,6)= XX*DX
      A(4,6)=-XX*DY
      A(5,1)=-XX*DX
      A(5,3)= XX*DY
      A(7,1)=0.D0
      A(7,3)=0.D0
      A(8,1)=0.D0
      A(8,3)=0.D0
      IF(NNSOL.EQ.1)RETURN
      A(7,1)=-TEMP*ZW(2,3)
      A(7,3)=-TEMP*ZW(1,3)
      A(8,1)= TEMP*ZW(2,2)
      A(8,3)= TEMP*ZW(1,2)
      RETURN
C
C
C=====SKEW QUADRUPOLE
    4 TEMP=(1.+S)*XX
      A(2,3)=-XX
      A(4,1)=-XX
      A(2,6)= XX*DY
      A(4,6)= XX*DX
      A(5,1)=-XX*DY
      A(5,3)=-XX*DX
      A(7,1)=0.D0
      A(7,3)=0.D0
      A(8,1)=0.D0
      A(8,3)=0.D0
      IF(NNSOL.EQ.1)RETURN
      A(7,1)= TEMP*ZW(1,3)
      A(7,3)=-TEMP*ZW(2,3)
      A(8,1)=-TEMP*ZW(1,2)
      A(8,3)= TEMP*ZW(2,2)
      RETURN
C
C
C=====HORIZONTAL DIPOLES WITH 1/RHO**2 FOCUSSING
    2 CONTINUE
      A(2,6)= XX
      A(5,1)=-XX
      A(2,1)=-XX*XX/YYY
      A(7,1)= 0.D0
      A(7,4)= 0.D0
      A(7,6)= 0.D0
      A(8,1)= 0.D0
      A(8,4)= 0.D0
      A(8,6)= 0.D0
      IF(NNSOL.EQ.1)RETURN
      XXS=XX*XX/YYY
      A(7,1)=-(1.+S)*XXS*ZW(2,3)
      A(7,4)=     S *XX *ZW(3,3)
      A(7,6)=        XX *ZW(2,3)
      A(8,1)=+(1.+S)*XXS*ZW(2,2)
      A(8,4)=    -S *XX *ZW(3,2)
      A(8,6)=       -XX *ZW(2,2)
      RETURN
C
C
C=====CAVITIES
    5 TEMP=YYY*(1.+S)
      A(6,5)= XX
      A(7,2)=-0.D0
      A(7,4)= 0.D0
      A(8,2)= 0.D0
      A(8,4)=-0.D0
      IF(NNSOL.EQ.1)RETURN
      A(7,2)=-TEMP*ZW(2,3)
      A(7,4)= TEMP*ZW(1,3)
      A(8,2)= TEMP*ZW(2,2)
      A(8,4)=-TEMP*ZW(1,2)
      RETURN
C
C
C=====VERTICAL DIPOLES WITH 1/RHO**2 FOCUSSING
    9 CONTINUE
      A(4,6)=-XX
      A(5,3)= XX
      A(4,3)=-XX*XX/YYY
      A(7,2)= 0.D0
      A(7,3)= 0.D0
      A(7,6)= 0.D0
      A(8,2)= 0.D0
      A(8,3)= 0.D0
      A(8,6)= 0.D0
      IF(NNSOL.EQ.1)RETURN
      XXS=XX*XX/YYY
      A(7,2)=     S *XX *ZW(3,3)
      A(7,3)=+(1.+S)*XXS*ZW(1,3)
      A(7,6)=        XX *ZW(1,3)
      A(8,2)=    -S *XX *ZW(3,2)
      A(8,3)=-(1.+S)*XXS*ZW(1,2)
      A(8,6)=       -XX *ZW(1,2)
      RETURN
C
C
C=====HORIZONTAL KICKER.
    6 CONTINUE
      A(2,6)= XX
      A(5,1)=-XX
      A(7,4)= 0.D0
      A(7,6)= 0.D0
      A(8,4)=-0.D0
      A(8,6)= 0.D0
      IF(NNSOL.EQ.1)RETURN
      A(7,4)= S*XX*ZW(3,3)
      A(7,6)=   XX*ZW(2,3)
      A(8,4)=-S*XX*ZW(3,2)
      A(8,6)=  -XX*ZW(2,2)
      RETURN
C
C
C=====VERTICAL KICKER.
    7 CONTINUE
      A(4,6)=-XX
      A(5,3)= XX
      A(7,2)= 0.D0
      A(8,2)=-0.D0
      A(7,6)= 0.D0
      A(8,6)=-0.D0
      IF(NNSOL.EQ.1)RETURN
      A(7,2)= XX*S*ZW(3,3)
      A(8,2)=-XX*S*ZW(3,2)
      A(7,6)= XX*ZW(1,3)
      A(8,6)=-XX*ZW(1,2)
C
C
C=====SEXTUPOLE
    8 A(2,1)=-XX*DX
      A(2,3)=XX*DY
      A(2,6)=XX*(DX*DX-DY*DY)/2.
      A(4,1)=XX*DY
      A(4,3)=XX*DX
      A(4,6)=-XX*DX*DY
      A(5,1)=-A(2,6)
      A(5,3)=-A(4,6)
      TEMP=(1.+S)*XX
      A(7,1)=0.D0
      A(7,3)=0.D0
      A(8,1)=0.D0
      A(8,3)=0.D0
      IF(NNSOL.EQ.1)RETURN
      A(7,1)=-TEMP*(DY*ZW(1,3)+DX*ZW(2,3))
      A(7,3)=-TEMP*(DX*ZW(1,3)-DY*ZW(2,3))
      A(8,1)= TEMP*(DY*ZW(1,2)+DX*ZW(2,2))
      A(8,3)= TEMP*(DX*ZW(1,2)-DY*ZW(2,2))
      RETURN
C
C
C=====SOLENOIDS ARE CALCULATED SEPARATELY
   10 RETURN
      END
