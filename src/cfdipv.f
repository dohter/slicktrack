C   06/05/84 602092042  MEMBER NAME  CFDIPV   (S)           FORTRAN
      SUBROUTINE CFDIPV(X,X2,YY,T,ECHROM)
C
C
C
C=====ROUTINE TO SET UP AN 7X7 THICK VERTICAL COMBINED FUNCTION
C                        DIPOLE MATRIX
C
C======X IS THE TOTAL BEND ANGLE,X2 IS THE K * LENGTH (I.E. FOCUSSING)
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 T(7,7)
C
      WRITE(6,1)
    1 FORMAT (' ',' VERTICAL COMBINED FUNCTION MAGNETS NOT YET CODED--
     +  SO STOP')
C
C
C
C
      STOP
      END
