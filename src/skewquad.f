C   06/05/84 201221453  MEMBER NAME  SKEWQD   (N3000.S)     FORTRAN
      SUBROUTINE SKEWQUAD(X,YY,IKICK,NTY,T,TILT)
C
C
C
C=====ROUTINE TO SET UP A 7X7 THICK NORMAL-SKEW-QUAD MATRIX.
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE "csol.for"

      REAL*8 T(7,7),T1(7,7),T2(7,7)
C
C=====SET AS A UNIT MATRIX TO BEGIN and on subsequent entries so that
C     the method of construction doesn't make a mess.
      CALL UNIT(7,T)
C
      IF(YY.EQ.0.D0)RETURN
C=====REDEFINE THE STRENGTH 'X'
      XX=X/YY
C=====INITIALISE MATRICES
      T(1,2)=YY
      T(3,4)=YY
      IF(XX.EQ.0.D0)RETURN
C
C
C
C=====INITIALISE MATRICES
      DO 2 J=1,7
      DO 2 K=1,7
        T1(J,K)=0.D0
        T2(J,K)=0.D0
    2 CONTINUE

      CTILT = DCOS(TILT)
      STILT = DSIN(TILT)
      DO 3 J=1,7
        T1(J,J)= CTILT
        T2(J,J)= CTILT
    3 CONTINUE
      T1(1,3)= STILT
      T1(2,4)= STILT
      T1(3,1)=-STILT
      T1(4,2)=-STILT
      T2(1,3)=-STILT
      T2(2,4)=-STILT
      T2(3,1)= STILT
      T2(4,2)= STILT
C
C
C
C=====CONSTRUCT A THICK SKEW-QUAD BY ROTATING A NORMAL QUAD.
      RXX=DSQRT(DABS(XX))
      FQ1=0.D0
      FQ2=1.D0
      IF(XX.GT.0.D0)FQ1=1.D0
      IF(XX.GT.0.D0)FQ2=0.D0
      ARG=RXX*YY
      DC=DCOS(ARG)
      DS=DSIN(ARG)
      DCH=DCOSH(ARG)
      DSH=DSINH(ARG)
      T(1,1)= FQ1*DC     + FQ2*DCH
      T(2,2)= FQ1*DC     + FQ2*DCH
      T(1,2)= FQ1*DS/RXX + FQ2*DSH/RXX
      T(2,1)=-FQ1*DS*RXX + FQ2*DSH*RXX
      T(3,3)= FQ2*DC     + FQ1*DCH
      T(4,4)= FQ2*DC     + FQ1*DCH
      T(3,4)= FQ2*DS/RXX + FQ1*DSH/RXX
      T(4,3)=-FQ2*DS*RXX + FQ1*DSH*RXX
C
      CALL JAM777(T,T,T1)
      CALL JAM777(T,T2,T)
C
C
C     DO 273 I=1,7
C 273 WRITE(53,925)(T(I,J),J=1,7)
C 925 FORMAT(7D16.8)
C     IF(1.EQ.1)STOP
C
C
C
      RETURN
      END
