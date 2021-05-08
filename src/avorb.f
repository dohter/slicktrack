C   15/11/79 503081829  MEMBER NAME  AVORB    (MAR95.S)     FORTRAN
        SUBROUTINE AVORB(A,B,N,C)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION A(N),B(N),C(N)
        DO 10 I=1,N
        DO 10 J=1,N
 10     C(I)=(A(I)+B(I))/2.
        RETURN
        END
