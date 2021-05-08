C   22/01/90 510031959  MEMBER NAME  NORM2    (SEPT95.S) M  FORTRAN
      SUBROUTINE NORM2(A,N,C)
C  *****************************************************************
C  *       SMILE version of NORM:  Modified to Fortran IV.         *
C  *       calculates eigenvector norm C(J), j=1,3,5 such that     *
C  *       F*E = -i where F is conj. transp. of E                  *
C  *       Even if N=8, only look at the orbital part.             *
C  *****************************************************************

      IMPLICIT REAL*8 (A-H,O-Z)

CDB      COMMON/JAMJAM/DJAM(20,20)

      DIMENSION A(N,N),C(N)

C===== eigenvectors are stored in A as column vectors, real and imaginary part alternating
      DO 10 J=1,6,2
        C(J) = 0.D0
        C(J+1) = 0.D0
        DO 10 I=1,6,2
          C(J)=C(J) + A(I,J)*A(I+1,J+1) - A(I+1,J)*A(I,J+1)
   10 CONTINUE

      DO 20 I=1,6
   20 C(I) = C(I)*2.D0
C
      WRITE(53,30) C(1), C(2), C(3), C(4), C(5), C(6)
   30 FORMAT('0','--NORM2-',6E15.5)
      RETURN
      END
