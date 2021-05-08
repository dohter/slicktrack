C   15/12/83 503061851  MEMBER NAME  ORTCHK   (MAR95.S)     FORTRAN
      SUBROUTINE ORTCHK(ROT)
C
C
C=====ROUTINE TO CHECK ORTHOGONALITY OF SPIN ROTATION MATRIX.
C=====BY MULTIPLYING WIT ITS INVERSE.

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ROT(3,3),ROTROT(3,3)
C
C
C
C
      DO 2 J=1,3
      DO 2 K=1,3
      ROTROT(J,K)=0.D0
      DO 1 L=1,3
    1 ROTROT(J,K)=ROTROT(J,K)+ROT(J,L)*ROT(K,L)
    2 CONTINUE
      ROTROT(1,1)=ROTROT(1,1)-1.D0
      ROTROT(2,2)=ROTROT(2,2)-1.D0
      ROTROT(3,3)=ROTROT(3,3)-1.D0
      WRITE(53,10)ROTROT
   10 FORMAT('0',' ORTHOG.TEST---',9E11.2)
C
C
      RETURN
      END
