C   06/06/84 406061310  MEMBER NAME  TRANS1   (S)           FORTRAN
C   15/11/79 205191641  MEMBER NAME  TRANS1   (S)           FORTRAN
        SUBROUTINE TRANS1(T,N,ITY,A,M)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION T(N,N,400),A(M,M)
        DO 10 I=1,M
        DO 10 J=1,M
 10     A(I,J)=T(I,J,ITY)
        RETURN
        END
