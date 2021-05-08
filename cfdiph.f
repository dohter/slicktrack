C   06/05/84 602092037  MEMBER NAME  CFDIPH   (S)           FORTRAN
      SUBROUTINE CFDIPH(X,X2,YY,T,ECHROM)
C
C
C
C=====ROUTINE TO SET UP AN 7X7 THICK HORIZONTAL COMBINED FUNCTION
C                        DIPOLE MATRIX
C
C======X IS THE TOTAL BEND ANGLE,X2 IS THE K * LENGTH (I.E. FOCUSSING)
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 T(7,7)
      REAL*8 ST(6,6),KX
C
C
C
C
C
C=====REDEFINE THE STRENGTH 'X2'
      X2=X2/(1.D0+ECHROM)
      KX=X/YY
      G=X2/YY
C=====INITIALISE MATRICES
      T(1,2)=YY
      T(3,4)=YY
      IF(KX.EQ.0.D0.AND.G.EQ.0D0)RETURN
C
C
C
C=====CONSTRUCT A THICK C.F.HORIZONTAL DIPOLE: KX IS THE CURVATURE.
C                                               G IS THE FOCUSSING.
C
C=====GET VARIOUS CASES DEPENDING ON THE SIGN AND SIZE OF G WRT KX.


      G1=KX*KX+G
      IF(G1.LT.0.D0)GO TO 1
C=====+VE OR ZERO G1
      RG1=DSQRT(G1)
      C=DCOS(RG1*YY)
      S=DSIN(RG1*YY)
      T(1,1)=C
      T(1,2)=S/RG1
      T(2,1)=-RG1*S
      T(2,2)=C
      T(1,6)=(1.D0-C)*KX/G1
      T(2,6)=S*KX/RG1
      T(5,1)=-T(2,6)
      T(5,2)=-T(1,6)
      T(5,6)=-KX*KX/G1*(YY-S/RG1)
      GO TO 2
    1 CONTINUE
C=====-VE G1
      RG1=DSQRT(-G1)
      CH=DCOSH(RG1*YY)
      SH=DSINH(RG1*YY)
      T(1,1)=CH
      T(1,2)=SH/RG1
      T(2,1)=RG1*SH
      T(2,2)=CH
      T(1,6)=-(1.D0-CH)*KX/(-G1)
      T(2,6)=SH*KX/RG1
C=====SIGNS OF THE NEXT 3 ELEMENTS CHANGED 27/06/85
      T(5,1)=-T(2,6)
      T(5,2)=-T(1,6)
      T(5,6)=+KX*KX/(-G1)*(YY-SH/RG1)
C
C
C=====NOW VERTICAL MOTION
C
    2 IF(G.EQ.0.D0)GO TO 4
      IF(G.LT.0.D0)GO TO 3
C=====+VE OR ZERO G
      G2=G
      RG2=DSQRT(G2)
      CH=DCOSH(RG2*YY)
      SH=DSINH(RG2*YY)
      T(3,3)=CH
      T(3,4)=SH/RG2
      T(4,3)=SH*RG2
      T(4,4)=CH
      GO TO 4
    3 CONTINUE
C=====-VE G
      G2=-G
      RG2=DSQRT(G2)
      C=DCOS(RG2*YY)
      S=DSIN(RG2*YY)
      T(3,3)=C
      T(3,4)=S/RG2
      T(4,3)=-S*RG2
      T(4,4)=C
C
C
    4 CONTINUE
C


C     DO 100 J=1,6
C     DO 100 K=1,6
C 100 ST(J,K)=T(J,K)
C     CALL SYMP(ST)
C
C
C
C
C
C
C
C
C     WRITE(6,12)
C  12 FORMAT(///)
C     DO 11 I=1,7
C  11 WRITE(6,10)(T(I,J),J=1,7)
C  10 FORMAT(' ',7E14.7)
C
C
C
C
      RETURN
      END
