C   15/11/79 406051926  MEMBER NAME  TRANS1   (S)           FORTRAN
      SUBROUTINE TRANS1(T,N,ITY,A,M)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION T(N,N),A(M,M)
      DO 10 I=1,M
      DO 10 J=1,M
   10 A(I,J)=T(I,J)
      RETURN
      END
