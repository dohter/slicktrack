C   15/11/79 411161900  MEMBER NAME  MX66     (MAY92.S)     FORTRAN
      SUBROUTINE MX66(T,ITY,ID,II,XX,X2,YY,A)
C
C=====GET 6X6 TERMS IN PRESENCE OF C.O. SHIFTS.
C
C
C=====PUT IN SEXTUPOLE  CURVATURE TERMS AS THIN LENS
C=====AT THE LENS CENTRE SANDWICHED BETWEEN HALF LENGTH DRIFTS?
C     SEXEQ,MX66,MX88,MXDAMP
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(6,6),T(7,7),B(6,6),C(6,6)
C
C
      INCLUDE "cnlist.for"
      INCLUDE "cloorb.for"
      INCLUDE "csol.for"
C
C
      I2ND=0
C     IF(NFREQ.GT.0)I2ND=1
C
C======DX etc are at the entrance to an element. Average from entrance to end.
      DELX =(DX(II) +DX(II+1))/2.D0
      DELY =(DY(II) +DY(II+1))/2.D0
      DELE =(DEL(II)+DEL(II+1))/2.D0
C
C
C      DO 10 I=1,6
C      DO 10 J=1,6
C   10 A(I,J)=T(I,J)
    
      A = T(1:6,1:6)

C      IF(ID.EQ.17)THEN                           !Ignore beam-beam kicks.
C      A(2,1)=0.D0
C      A(4,3)=0.D0
C      ENDIF

      IF(ID.EQ.8)GO TO 11
      IF(ID.EQ.3)GO TO 12
      IF(ID.EQ.4)GO TO 13
      IF(ID.EQ.15)GO TO 15
      RETURN
C
C=====QUADRUPOLE: WITH C.O. SHIFT FUDGE THE CURVATURE STUFF BY PUTTING
C=====IN DISPERSION GENERATING KICKERS OF 1/2 STRENGTH FORE & AFT.
C=====SINCE WE ALREADY HAVE SLICES THAT SHOULD NOT BE TOO BAD.
   12 CONTINUE
      CALL UNIT(6,B)
      B(2,6)= XX*DELX/2.D0     *(1.D0-DELE*I2ND)  *1.D0
      B(4,6)=-XX*DELY/2.D0     *(1.D0-DELE*I2ND)  *1.D0
      B(5,1)=-XX*DELX/2.D0     *(1.D0-DELE*I2ND)  *1.D0
      B(5,3)= XX*DELY/2.D0     *(1.D0-DELE*I2ND)  *1.D0
      CALL JAM666(A,A,B)
      CALL JAM666(A,B,A)
      RETURN
C
C=====SKEW QUADRUPOLE
C=====THESE ALSO NEED THE DISPERSION STUFF.
   13 CONTINUE
C     A(2,6)= XX*DELY
C     A(4,6)= XX*DELX
C     A(5,1)=-XX*DELY
C     A(5,3)=-XX*DELX
      RETURN
   15 CONTINUE
C=====Horizontal combined function: add in vertical curvature effect
C     due to vertical closed orbit shift.
C=====NOV 94. NOT YET CHECKED!!!!  FIX MX88 ALSO.
      DEX=DX(II)
      DEY=DY(II)
      A(4,6)=-X2*DEY     *0.D0      !Kill Sept 2003, it breaks symp.
      A(5,3)= X2*DEY     *0.D0
      RETURN
C
C=====SEXTUPOLE: THIN LENS SANDWICHED BETWEEN 2 HALF LENGTH DRIFTS.
   11 CONTINUE
      CALL UNIT(6,B)
      CALL UNIT(6,C)
      C(1,2)=YY*0.5D0
      C(3,4)=YY*0.5D0
      XS=XX*(1.D0-DELE*I2ND)
      B(2,1)=-XS*DELX
      B(2,3)= XS*DELY
      B(2,6)= XS*(DELX**2-DELY**2)/2.
      B(4,1)= B(2,3)
      B(4,3)=-B(2,1)
      B(4,6)=-XS*DELX*DELY
      B(5,1)=-B(2,6)
      B(5,3)=-B(4,6)
      CALL JAM666(A,B,C)
      CALL JAM666(A,C,A)
      RETURN
      END
