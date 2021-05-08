C   15/11/79 411071514  MEMBER NAME  SIG      (S)           FORTRAN
      SUBROUTINE SIG(AA1,AA3,AA5,A,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(6,6),B(6,6)
C
C
C
      DO 10 I=1,6
      DO 10 J=1,6
      B(I,J)=
     +       2.*AA1*(A(I,1)*A(J,1)+A(I,2)*A(J,2))
     +      +2.*AA3*(A(I,3)*A(J,3)+A(I,4)*A(J,4))
     +      +2.*AA5*(A(I,5)*A(J,5)+A(I,6)*A(J,6))
   10 CONTINUE
C
C
      RETURN
      END
