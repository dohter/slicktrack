C   15/11/79 411062115  MEMBER NAME  NORM     (S)           FORTRAN
      SUBROUTINE NORM(A,N,AB)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,N)
      DIMENSION AB(6)
C
      AB(1)=-A(1,1)*A(2,2)+A(2,1)*A(1,2)
     +      -A(3,1)*A(4,2)+A(4,1)*A(3,2)
     +      -A(5,1)*A(6,2)+A(6,1)*A(5,2)
C
      AB(3)=-A(1,3)*A(2,4)+A(2,3)*A(1,4)
     +      -A(3,3)*A(4,4)+A(4,3)*A(3,4)
     +      -A(5,3)*A(6,4)+A(6,3)*A(5,4)
C
      AB(5)=-A(1,5)*A(2,6)+A(2,5)*A(1,6)
     +      -A(3,5)*A(4,6)+A(4,5)*A(3,6)
     +      -A(5,5)*A(6,6)+A(6,5)*A(5,6)
C
      AB(2)=AB(1)
      AB(4)=AB(3)
      AB(6)=AB(5)
C
C
      WRITE(53,10)AB(1),AB(3),AB(5)
   10 FORMAT('0','--NORM--',3E20.5)
      RETURN
      END
