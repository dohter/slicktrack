C   15/11/79 309261523  MEMBER NAME  SEX88    (S)           FORTRAN
      SUBROUTINE SEX88(XX,S,DX,DY,ZW,A,NSOL)
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 DX,DY
      DIMENSION A(8,8),ZW(3,3)
C
C
C      CALL VZERO(A,128)
      A     = 0.D0 
      A(1,1)=1
      A(2,2)=1
      A(3,3)=1
      A(4,4)=1
      A(5,5)=1
      A(6,6)=1
      A(7,7)=1
      A(8,8)=1
      A(2,1)=-XX*DX
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
      IF(NSOL.NE.0)RETURN
      NX=1
      NY=1
C     NX=0
      A(7,1)=-TEMP*(DY*ZW(1,3)+DX*ZW(2,3))*NX
      A(7,3)=-TEMP*(DX*ZW(1,3)-DY*ZW(2,3))*NY
      A(8,1)= TEMP*(DY*ZW(1,2)+DX*ZW(2,2))*NX
      A(8,3)= TEMP*(DX*ZW(1,2)-DY*ZW(2,2))*NY
      RETURN
      END
