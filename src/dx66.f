C   15/11/79 603171550  MEMBER NAME  DX66     (S)           FORTRAN
      SUBROUTINE DX66(T,ITY,ID,II,XX,X2,YY,A,AI)

C
C
C=====GENERATE THE 6X6 USING THE UNCOUPLED DISPERSION FORMALISM.
C=====THIS IS A STRIPPED DOWN VERSION OF DX88
C=====THICK LENS VERSION: PICK UP MATRICES FROM TMAT
C                         C.O. ERRORS NOT YET ALLOWED


      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(6,6),B(6,6),AI(6,6)
      DIMENSION T(7,7)
      REAL*8 KX
C
C
      INCLUDE "cnlist.for"
      INCLUDE "csol.for"
      INCLUDE "cloorb.for"
      INCLUDE "cdisp.for"
C
C
C
C
C=====SET UP INITIAL OPTICAL PART (4X4) OF THE 6X6 MATRIX & INVERSE.
C=====2X2 SYNCH. PART IS FILLED LATER.
C      CALL VZERO(A,72)
C      CALL VZERO(AI,72)
      A  = 0.D0
      AI = 0.D0
      DO 1 J=1,4
      DO 1 K=1,4
      AI(J,K)=T(J,K)
    1 A(J,K) =T(J,K)
      A(5,5)= 1.D0
      A(6,6)= 1.D0
      AI(5,5)=1.D0
      AI(6,6)=1.D0
C
C
C=====USE AVERAGE CLOSED ORBIT SHIFTS.
      DELX =(DX(II) +DX(II+1))/2.
      DELY =(DY(II) +DY(II+1))/2.
      DELE =(DEL(II)+DEL(II+1))/2.
C
      NNSOL=NSOL(ITY)
C=====GET THE INPUT DISPERSIONS: FOR ELEMENT(5,6)
      D1I=D1(II)
      D2I=D2(II)
      D3I=D3(II)
      D4I=D4(II)
      D1D=D1I
      D2D=D2I
      D3D=D3I
      D4D=D4I
C     WRITE(53,41)ID,D1I,D2I,D3I,D4I
C  41 FORMAT(' ',I4,4F15.5)
C     IF(NNSOL.EQ.1)D1I=0.D0
C     IF(NNSOL.EQ.1)D2I=0.D0
C     NNSOL=0
C
C
C
C
      GOTO(11,12,13,14,15,18,19,20,17,10,11,11,11,11,21,22),ID
C=====DRIFT
   11 CONTINUE
      A(1,2)=YY
      A(3,4)=YY
      AI(1,2)=-YY
      AI(3,4)=-YY
      RETURN
C
C=====HORIZONTAL DIPOLES
   12 CONTINUE
      AI(1,2)=-A(1,2)
      AI(2,1)=-A(2,1)
      AI(3,4)=-A(3,4)
      AI(4,3)=-A(4,3)
      A(5,6)=0.D0
      IF(XX.NE.0.D0)
     +      A(5,6)=-D1D*XX/YY*A(1,2)+D2D*(A(1,1)-1.D0)*YY/XX-(YY-A(1,2))
      AI(5,6)=-A(5,6)
      RETURN
C
C=====QUADRUPOLES.
   13 CONTINUE
      AI(1,2)=-A(1,2)
      AI(2,1)=-A(2,1)
      AI(3,4)=-A(3,4)
      AI(4,3)=-A(4,3)
      RETURN
C
C=====SKEW QUADS.
   14 CONTINUE
C=====INVERSE NOT YET CODED
      RETURN
C
C=====CAVITIES
   15 CONTINUE
      A(6,5)=  T(6,5)
      AI(6,5)=-A(6,5)
      RETURN
C
C=====VERTICAL DIPOLES
   17 CONTINUE
      AI(1,2)=-A(1,2)
      AI(2,1)=-A(2,1)
      AI(3,4)=-A(3,4)
      AI(4,3)=-A(4,3)
      A(5,6)=0.D0
      IF(XX.NE.0.D0)
     +    A(5,6)=-D3D*XX/YY*A(3,4)+D4D*(A(3,3)-1.D0)*YY/XX-(YY-A(3,4))
      AI(5,6)=-A(5,6)
      RETURN
C
C=====HORIZONTAL KICKER.
   18 CONTINUE
      RETURN
C
C=====VERTICAL KICKER.
   19 CONTINUE
      RETURN
C
C=====SEXTUPOLE:CODE AS THIN LENS BETWEEN 2 HALF DRIFTS.
   20 CONTINUE
C=====NOT PROPERLY CODED YET:PUT C.O. TO ZERO.
C
      A(1,2)=YY
      A(3,4)=YY
      AI(1,2)=-YY
      AI(3,4)=-YY
      RETURN
C=====SOLENOID: TREATED A A DRIFT
   10 CONTINUE
      A(1,2)=YY
      A(3,4)=YY
      AI(1,2)=-YY
      AI(3,4)=-YY
      RETURN
C
C=====HORIZONTAL C.F. DIPOLE.
   21 CONTINUE
      AI(1,2)=-A(1,2)
      AI(2,1)=-A(2,1)
      AI(3,4)=-A(3,4)
      AI(4,3)=-A(4,3)
      KX=XX/YY
      G=X2/YY
      G1=KX*KX+G
      A(5,6)=0.D0
      IF(XX.EQ.0.D0.AND.X2.EQ.0.D0)RETURN
      A(5,6)=-D1D*KX*A(1,2)+D2D*KX/G1*(A(1,1)-1.D0)-(YY-A(1,2))*KX*KX/G1
      AI(5,6)=-A(5,6)
C
      RETURN
C=====VERTICAL C.F. DIPOLE.
   22 CONTINUE
      STOP
C
C
      END
