C   15/11/79            MEMBER NAME  TRANS    (SLIMS)       FORTRAN
        SUBROUTINE TRANS(A,N,B)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION A(N,N),B(N,N)
        DO 10 I=1,N
        DO 10 J=1,N
 10     B(I,J)=A(I,J)
        RETURN
        END
