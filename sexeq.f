
      SUBROUTINE SEXEQ(IRAD,S,II,ITY,A)
C
C=====CALCULATES SEXTUPOLE TERMS IN PRESENCE OF A C.O. SHIFT USING
C=====A THIN LENS SANDWICHED BETWEEN 2 HALF LENGTH DRIFTS SO THAT IT IS
C=====SYMPLECTIC.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(7,7),B(7,7),C(7,7)
C
      INCLUDE 'cloorb.for'
      INCLUDE 'clatic.for'
C
C======DX etc are at the ENTRANCE to an element. Average from entrance to end.
      DELX=(DX(II)+DX(II+1))/2.D0
      DELY=(DY(II)+DY(II+1))/2.D0
C
C
C=====JUNE97 KILL CURVATURE TERMS FOR C.O GENERATION.
C
C=====SET UP THIN LENS TERMS
      CALL UNIT(7,B)
      CALL UNIT(7,C)
      C(1,2)=YY(ITY)*0.5D0
      C(3,4)=YY(ITY)*0.5D0
      B(1,2)= 0.D0
      B(3,4)= 0.D0
C     B(2,6)=XX(ITY)*(DELX**2-DELY**2)/2.D0
C     B(5,1)=-B(2,6)
C     B(2,7)=-B(2,6)
      B(2,7)=-XX(ITY)*(DELX**2-DELY**2)/2.D0
C     B(4,6)=-XX(ITY)*DELX*DELY
C     B(5,3)=-B(4,6)
C     B(4,7)=-B(4,6)
      B(4,7)=+XX(ITY)*DELX*DELY
C=====ADD IN THE QUAD/SKEW QUAD TERMS:BARBER JULY 1982
C=====REMOVE COMMENTS AGAIN 5/11/84 TO MAKE CONSISTENT WITH MX66 ETC.
      B(2,1)=-XX(ITY)*DELX
      B(2,3)= XX(ITY)*DELY
      B(4,1)= B(2,3)
      B(4,3)=-B(2,1)
      IF(IRAD.EQ.0)GO TO 10
      B(6,7)=-S*XX(ITY)**2*(DELX**2+DELY**2)**2/(4.D0*YY(ITY))
C=====DO THIS IN LINE LATER WHEN IT IS WORKING
   10 CALL JAM777(A,B,C)
      CALL JAM777(A,C,A)
      RETURN
      END
