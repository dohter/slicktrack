C   15/11/79 510241909  MEMBER NAME  UNIT     (S)           FORTRAN
      SUBROUTINE UNIT(N,U)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION U(N,N)
      DO 10 I=1,N
      DO 10 J=1,N
      U(I,J)=0.D0
   10 IF(I.EQ.J)U(I,J)=1.D0
      RETURN
      END
