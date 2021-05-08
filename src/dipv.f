C   06/05/84 603061441  MEMBER NAME  DIPV     (S)           FORTRAN
      SUBROUTINE DIPV(X,YY,T)
C
C
C
C=====ROUTINE TO SET UP AN 7X7 THICK VERTICAL DIPOLE MATRIX
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 T(7,7)


C
C
C
C
      XX=X/YY
C=====INITIALISE MATRICES
      T(1,2)=YY
      T(3,4)=YY
      IF(XX.EQ.0.D0)RETURN
C
C
C
C=====CONSTRUCT A THICK VERTICAL DIPOLE: XX IS THE CURVATURE.
      ANG=X
      C=DCOS(ANG)
      S=DSIN(ANG)
      T(3,3)=C
      T(3,4)=S/XX
      T(4,3)=-XX*S
      T(4,4)=C
      T(3,6)=(1.D0-C)/XX
      T(4,6)=S
      T(5,3)=-S
      T(5,4)=-(1.D0-C)/XX
      T(5,6)=-(YY-S/XX)
C
C
C
C
C
C
      RETURN
      END
