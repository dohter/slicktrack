C   15/11/79            MEMBER NAME  AVER     (SLIMS)       FORTRAN
        SUBROUTINE AVER(A,B,N,C)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION A(N,N),B(N,N),C(N,N)
        DO I=1,N
          DO J=1,N
            C(I,J)=(A(I,J)+B(I,J))/2.
          ENDDO
        ENDDO

        RETURN
        END
