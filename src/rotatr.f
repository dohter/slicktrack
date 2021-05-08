C   10/02/84 406061632  MEMBER NAME  ROTATR   (S)           FORTRAN
      SUBROUTINE ROTATR(NTY)
C
C
C=====ROUTINE TO READ IN AND ASSIGN SPECIAL TMATRIX ELEMENTS TO
C=====A COUPLED SECTION OF THE RING WHICH OVERALL IS UNCOUPLED.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE "clatic.for"
      INCLUDE "csol.for"
C
C
C=====READ 4 ROWS OF MATRIX ELEMENTS.
    1 FORMAT(A1,I4,1X,A4,2F12.8,I5)
      READ(52,1)KOM,ID(NTY),NAME(NTY),TMAT(1,1,NTY),TMAT(1,2,NTY),
     +         NUNT(NTY)
      READ(52,1)KOM,ID(NTY),NAME(NTY),TMAT(2,1,NTY),TMAT(2,2,NTY),
     +         NUNT(NTY)
      READ(52,1)KOM,ID(NTY),NAME(NTY),TMAT(3,3,NTY),TMAT(3,4,NTY),
     +         NUNT(NTY)
      READ(52,1)KOM,ID(NTY),NAME(NTY),TMAT(4,3,NTY),TMAT(4,4,NTY),
     +         NUNT(NTY)
C
C
      WRITE(53,40)NTY
   40 FORMAT(' ', 'ELEMENT ',I4,' HAS THE FOLLOWING MATRIX')
      DO 20 J=1,7
   20 WRITE(53,30)(TMAT(J,K,NTY),K=1,7)
   30 FORMAT(' ',30X,7F12.5)
C
C
      RETURN
      END