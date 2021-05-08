C   15/11/79 311071918  MEMBER NAME  SOL66    (S)           FORTRAN
      SUBROUTINE SOL66(XX,YY,A)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(4,4),B(4,4),C(4,4)
      CALL UNITM(4,B)
      B(2,3)=0.5*XX/YY
      B(4,1)=-B(2,3)
      CALL SOL6(XX,YY,A)
      CALL JAM444(C,A,B)
      B(2,3)=-B(2,3)
      B(4,1)=-B(4,1)
      CALL JAM444(A,B,C)
      RETURN
      END
