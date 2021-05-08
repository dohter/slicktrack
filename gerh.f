C   15/11/79 501201723  MEMBER NAME  GERH     (S)           FORTRAN
      SUBROUTINE MX88(ID,II,ITY,T,XX,X2,YY,S,ZWM,ZWI,A)
C
C
C=====GENERATE THE 8X8 MATRIX USING THE SPIN BASIS.
C=====THICK LENS VERSION: PICK UP MATRICES FROM TMAT
C


      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ZWM(3,3),ZWI(3,3),A(8,8)
      DIMENSION T(7,7)
      REAL*8 LXM,LZM,LSM,MXM,MZM,MSM
      REAL*8 LXI,LZI,LSI,MXI,MZI,MSI
      REAL*8 KX
C
C
#include "csol.for"
#include "cloorb.for"
C
C
C
C
C=====SET UP INITIAL OPTICAL PART OF 8X8 MATRIX: SPIN PART REMAINS EMPTY
C=====UNLESS FILLED LATER.
      CALL VZERO(A,128)
      DO 1 J=1,6
      DO 1 K=1,6
    1 A(J,K)=T(J,K)
      A(7,7)=1.D0
      A(8,8)=1.D0
C
C
C=====USE AVERAGE CLOSED ORBIT SHIFTS.
      DELX=(DX(II)+DX(II+1))/2.
      DELY=(DY(II)+DY(II+1))/2.
      DELE =(DEL(II)+DEL(II+1))/2.
C
C
C=====SET SPIN VECTORS---AVERAGED INPUT & OUTPUT VALUES.
      LXM=ZWM(1,3)
      LZM=ZWM(2,3)
      LSM=ZWM(3,3)
      MXM=ZWM(1,2)
      MZM=ZWM(2,2)
      MSM=ZWM(3,2)
C=====SET SPIN VECTORS---INPUT VALUES.
      LXI=ZWI(1,3)
      LZI=ZWI(2,3)
      LSI=ZWI(3,3)
      MXI=ZWI(1,2)
      MZI=ZWI(2,2)
      MSI=ZWI(3,2)
C
C
C
C
C
C
      NNSOL=NSOL(ITY)
      GOTO(11,12,13,14,15,18,19,20,17,11,11,11,11,11,21,17),ID
   11 RETURN
C
C=====HORIZONTAL DIPOLES
   12 CONTINUE
      IF(NNSOL.EQ.1)RETURN
      SO=DSIN(-S*XX)
      CO=DCOS(-S*XX)
      AGP1=(1.D0+S)
      A(7,1)= AGP1* A(2,1)      *LZI
      A(8,1)=-AGP1* A(2,1)      *MZI
      A(7,2)= AGP1*(A(2,2)-1.D0)*LZI
      A(8,2)=-AGP1*(A(2,2)-1.D0)*MZI
      A(7,4)=-SO*LSI-(CO-1.D0)*LXI
      A(8,4)= SO*MSI+(CO-1.D0)*MXI
      A(7,6)=-AGP1*(XX-A(2,6))*LZI+XX*LZI
      A(8,6)= AGP1*(XX-A(2,6))*MZI-XX*MZI
      RETURN
C
C=====QUADRUPOLES.
   13 CONTINUE
      AGP1=(1.D0+S)
      A(2,6)= XX*DELX
      A(4,6)=-XX*DELY
      A(5,1)=-XX*DELX
      A(5,3)= XX*DELY
      IF(NNSOL.EQ.1)RETURN
      A(7,1)= AGP1* A(2,1)      *LZM
      A(8,1)=-AGP1* A(2,1)      *MZM
      A(7,2)= AGP1*(A(2,2)-1.D0)*LZM
      A(8,2)=-AGP1*(A(2,2)-1.D0)*MZM
      A(7,3)=-AGP1* A(4,3)      *LXM
      A(8,3)= AGP1* A(4,3)      *MXM
      A(7,4)=-AGP1*(A(4,4)-1.D0)*LXM
      A(8,4)= AGP1*(A(4,4)-1.D0)*MXM
      RETURN
C
C=====SKEW QUADS.
   14 AGP1=(1.+S)
      A(2,6)= XX*DELY
      A(4,6)= XX*DELX
      A(5,1)=-XX*DELY
      A(5,3)=-XX*DELX
      IF(NNSOL.EQ.1)RETURN
      A(7,1)= AGP1*(A(2,1)*LZM-A(2,3)*LXM)
      A(8,1)=-AGP1*(A(2,1)*MZM-A(2,3)*MXM)
      A(7,2)= AGP1*((A(2,2)-1.D0)*LZM-A(2,4)*LXM)
      A(8,2)=-AGP1*((A(2,2)-1.D0)*MZM-A(2,4)*MXM)
      A(7,3)= AGP1*(A(2,3)*LZM-A(2,1)*LXM)
      A(8,3)=-AGP1*(A(2,3)*MZM-A(2,1)*MXM)
      A(7,4)= AGP1*(A(2,4)*LZM-(A(2,2)-1.D0)*LXM)
      A(8,4)=-AGP1*(A(2,4)*MZM-(A(2,2)-1.D0)*MXM)
      RETURN
C
C=====CAVITIES
C=====12/SEPT /83 ----MAKE IT THE SAME AS A.CHAO'S PAPERS
   15 TEMP=YY*(1.+S)
      IF(NNSOL.EQ.1)RETURN
      A(7,2)=-TEMP*LZM
      A(7,4)= TEMP*LXM
      A(8,2)= TEMP*MZM
      A(8,4)=-TEMP*MXM
      RETURN
C
C=====VERTICAL DIPOLES
   17 CONTINUE
      IF(NNSOL.EQ.1)RETURN
      SO=DSIN(-S*XX)
      CO=DCOS(-S*XX)
      AGP1=(1.D0+S)
      A(7,2)= SO*LSI-(CO-1.D0)*LZI
      A(8,2)=-SO*MSI+(CO-1.D0)*MZI
      A(7,3)=-AGP1* A(4,3)      *LXI
      A(8,3)= AGP1* A(4,3)      *MXI
      A(7,4)=-AGP1*(A(4,4)-1.D0)*LXI
      A(8,4)= AGP1*(A(4,4)-1.D0)*MXI
      A(7,6)= AGP1*(XX-A(4,6))*LXI-XX*LXI
      A(8,6)=-AGP1*(XX-A(4,6))*MXI+XX*MXI
      RETURN
C
C=====HORIZONTAL KICKER.
   18 CONTINUE
      IF(NNSOL.EQ.1)RETURN
      A(7,4)= S*XX*LSM
      A(7,6)=   XX*LZM
      A(8,4)=-S*XX*MSM
      A(8,6)=  -XX*MZM
      RETURN
C
C=====VERTICAL KICKER.
   19 CONTINUE
      IF(NNSOL.EQ.1)RETURN
      A(7,2)= XX*S*LSM
      A(7,6)=   XX*LXM
      A(8,2)=-XX*S*MSM
      A(8,6)=  -XX*MXM
C
C=====SEXTUPOLE: AS THICK/THIN LENS MIXTURE.
   20 A(2,1)=-XX*DELX
      A(2,3)= XX*DELY
      A(2,6)= XX*(DELX*DELX-DELY*DELY)/2.
      A(4,1)= XX*DELY
      A(4,3)= XX*DELX
      A(4,6)=-XX*DELX*DELY
      A(5,1)=-A(2,6)
      A(5,3)=-A(4,6)
      TEMP=(1.+S)*XX
      IF(NNSOL.EQ.1)RETURN
      A(7,1)=-TEMP*(DELY*LXM+DELX*LZM)
      A(7,3)=-TEMP*(DELX*LXM-DELY*LZM)
      A(8,1)= TEMP*(DELY*MXM+DELX*MZM)
      A(8,3)= TEMP*(DELX*MXM-DELY*MZM)
C
C=====HORIZONTAL C.F. DIPOLE.
   21 CONTINUE
      IF(NNSOL.EQ.1)RETURN
      AGP1=(1.D0+S)
      A(7,1)= AGP1* A(2,1)      *LZI
      A(8,1)=-AGP1* A(2,1)      *MZI
      A(7,2)= AGP1*(A(2,2)-1.D0)*LZI
      A(8,2)=-AGP1*(A(2,2)-1.D0)*MZI
C
C
      A(7,6)=-AGP1*(XX-A(2,6))*LZI+XX*LZI
      A(8,6)= AGP1*(XX-A(2,6))*MZI-XX*MZI
C
C     A(7,4)=-SO*LSI-(CO-1.D0)*LXI
C     A(8,4)= SO*LZI-(CO-1.D0)*MXI
      IF(XX.EQ.0.D0.AND.X2.EQ.0.D0)RETURN
      A3=A(4,3)
      A4=A(4,4)
      SO=DSIN(-S*XX)
      CO=DCOS(-S*XX)
      KX=XX/YY
      OM=-S*KX
      G=X2/YY
C
C=====+VE OR ZERO G
      IF(G.LT.0.D0)GO TO 31
      RG=DSQRT(G)
      A(7,3)=-AGP1*G*LSI*(  A3*SO-OM*A4*CO+OM)
     +       -AGP1*G*LXI*(  A3*CO+OM*A4*SO   )
     +       +  S*KX*LSI*(G*A4*CO+OM*A3*SO- G)
     +       -  S*KX*LXI*(G*A4*SO-OM*A3*CO   )
      A(7,3)=A(7,3)/(OM*OM+G)
      A(8,3)=-AGP1*G*MSI*(  A3*SO-OM*A4*CO+OM)
     +       -AGP1*G*MXI*(  A3*CO+OM*A4*SO   )
     +       +  S*KX*MSI*(G*A4*CO+OM*A3*SO- G)
     +       -  S*KX*MXI*(G*A4*SO-OM*A3*CO   )
      A(8,3)=-A(8,3)/(OM*OM+G)
      A(7,4)=-AGP1*LSI*  (G*A4*SO-OM*A3*CO   )
     +       -AGP1*LXI*  (G*A4*CO+OM*A3*SO- G)
     +       +S*KX*LSI*  (  A3*CO+OM*A4*SO   )
     +       -S*KX*LXI*  (  A3*SO-OM*A4*CO+OM)
      A(7,4)=A(7,4)/(OM*OM+G)
      A(8,4)=-AGP1*MSI*  (G*A4*SO-OM*A3*CO   )
     +       -AGP1*MXI*  (G*A4*CO+OM*A3*SO- G)
     +       +S*KX*MSI*  (  A3*CO+OM*A4*SO   )
     +       -S*KX*MXI*  (  A3*SO-OM*A4*CO+OM)
      A(8,4)=-A(8,4)/(OM*OM+G)
C     A(7,4)=-SO*LSI-(CO-1.D0)*LXI
C     A(8,4)= SO*MSI+(CO-1.D0)*MXI
      RETURN
   31 CONTINUE
C=====G<0
      G=-X2/YY
      RG=DSQRT(G)
      A(7,3)=-AGP1*G*LSI*(  A3*SO-OM*A4*CO+OM)
C
     +       -AGP1*G*LXI*(  A3*CO+OM*A4*SO   )
C
     +       +  S*KX*LSI*(G*A4*CO-OM*A3*SO- G)
C
     +       -  S*KX*LXI*(G*A4*SO+OM*A3*CO   )
C
      A(7,3)=A(7,3)/(G-OM*OM)
C
C
C
      A(7,4)=-AGP1*LSI*  (G*A4*SO+OM*A3*CO   )
C
     +       -AGP1*LXI*  (G*A4*CO-OM*A3*SO- G)
C
     +       +S*KX*LSI*  (  A3*CO+OM*A4*SO   )
C
     +       -S*KX*LXI*  (  A3*SO-OM*A4*CO+OM)
C
      A(7,4)=A(7,4)/(G-OM*OM)
C
C
      RETURN
C
C
      END
