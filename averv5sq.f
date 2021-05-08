C   15/11/79            MEMBER NAME  AVER     (SLIMS)       FORTRAN
        SUBROUTINE AVERV5SQ(A,B,V5SQ1,V5SQ3,V5SQ5)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION A(6,6),B(6,6)

C======Average the squares of the fifth components.

        V5SQ1 = 0.5D0*(A(5,1)**2 + A(5,2)**2 + B(5,1)**2 + B(5,2)**2)
        V5SQ3 = 0.5D0*(A(5,3)**2 + A(5,4)**2 + B(5,3)**2 + B(5,4)**2)
        V5SQ5 = 0.5D0*(A(5,5)**2 + A(5,6)**2 + B(5,5)**2 + B(5,6)**2)


        RETURN
        END
