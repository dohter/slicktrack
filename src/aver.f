C   15/11/79            MEMBER NAME  AVER     (SLIMS)       FORTRAN
        SUBROUTINE AVER(A,B,N,C)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION A(N,N),B(N,N),C(N,N)
        DO 10 I=1,N
        DO 10 J=1,N
 10     C(I,J)=(A(I,J)+B(I,J))/2.
        RETURN
        END
