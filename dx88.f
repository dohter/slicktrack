C   15/11/79 509051800  MEMBER NAME  DX88     (SEPT95.S) M  FORTRAN
      SUBROUTINE DX88(ID,II,ITY,T,XX,X2,YY,S,ZWM,ZWI,A)

C
C
C=====GENERATE THE 8X8 MATRIX USING THE SPIN BASIS.
C     AND THE DISPERSION FORMALISM.


C     INCLUDE BEAM-BEAM CENTRED ON THE C.O.
C!!!!!MX66 DOES NOT INCLUDED BEAM-BEAM. SO THE STORED DISPERSIONS CAN BE OFF!!!!!
C
C=====THICK LENS VERSION: PICK UP MATRICES FROM TMAT
C                         C.O. ERRORS NOT YET ALLOWED


      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ZWM(3,3),ZWI(3,3),A(8,8),B(8,8)
      DIMENSION T(7,7)
      REAL*8 LXM,LZM,LSM,MXM,MZM,MSM
      REAL*8 LXI,LZI,LSI,MXI,MZI,MSI
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
C=====SET UP INITIAL OPTICAL PART (4X4) OF THE 8X8 MATRIX.
C=====THE SPIN PART AND THE 2X2 SYNCH. PART ARE FILLED LATER.
C      CALL VZERO(A,128)
      A = 0.D0
      DO 1 J=1,4
      DO 1 K=1,4
    1 A(J,K)=T(J,K)
      A(5,5)=1.D0
      A(6,6)=1.D0
      A(7,7)=1.D0
      A(8,8)=1.D0
C
C
C=====USE AVERAGE CLOSED ORBIT SHIFTS.
      DELX =(DX(II) +DX(II+1))/2.
      DELY =(DY(II) +DY(II+1))/2.
      DELE =(DEL(II)+DEL(II+1))/2.
C
      NNSOL=NSOL(ITY)
C=====GET THE INPUT DISPERSIONS: FOR ELEMENT(5,6) & FOR SPIN.
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
      GOTO(11,12,13,14,15,18,19,20,17,11,11,11,11,11,21,17,22),ID
   11 RETURN
C
C=====HORIZONTAL DIPOLES
   12 CONTINUE
      A(5,6)=0.D0
      IF(XX.NE.0.D0)
     +      A(5,6)=-D1D*XX/YY*A(1,2)+D2D*(A(1,1)-1.D0)*YY/XX-(YY-A(1,2))
      IF(NNSOL.EQ.1)RETURN
      SO=DSIN(-S*XX)
      CO=DCOS(-S*XX)
      AGP1=(1.D0+S)
      A(7,1)= AGP1* A(2,1)      *LZI *1.D0
      A(8,1)=-AGP1* A(2,1)      *MZI *1.D0
      A(7,2)= AGP1*(A(2,2)-1.D0)*LZI *1.D0
      A(8,2)=-AGP1*(A(2,2)-1.D0)*MZI *1.D0
      A(7,4)=-SO*LSI-(CO-1.D0)*LXI
      A(8,4)= SO*MSI+(CO-1.D0)*MXI
C     IF(NNSOL.EQ.3)D1I=0.D0
C     IF(NNSOL.EQ.3)D2I=0.D0
      A(7,6)= A(7,1)*D1I+A(7,2)*D2I+A(7,3)*D3I+A(7,4)*D4I+XX*LZI *1.D0
      A(8,6)= A(8,1)*D1I+A(8,2)*D2I+A(8,3)*D3I+A(8,4)*D4I-XX*MZI *1.D0
      RETURN
C
C=====QUADRUPOLES.
   13 CONTINUE
C     IF(NNSOL.EQ.3)D1I=0.D0
C     IF(NNSOL.EQ.3)D2I=0.D0
C
      AGP1=(1.D0+S*(1.D0+ECHROM))
C     AGP1=(1.D0+S*(1.D0       ))
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
      IF(NNSOL.EQ.2)RETURN
      A(7,6)= A(7,1)*D1I+A(7,2)*D2I
      A(8,6)= A(8,1)*D1I+A(8,2)*D2I
      A(7,6)= A(7,6)+A(7,3)*D3I+A(7,4)*D4I
      A(8,6)= A(8,6)+A(8,3)*D3I+A(8,4)*D4I
      RETURN
C
C=====SKEW QUADS.
   14 AGP1=(1.D0+S*(1.D0+ECHROM))
C  14 AGP1=(1.D0+S*(1.D0       ))
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
      A(7,6)= A(7,1)*D1I+A(7,2)*D2I+A(7,3)*D3I+A(7,4)*D4I
      A(8,6)= A(8,1)*D1I+A(8,2)*D2I+A(8,3)*D3I+A(8,4)*D4I
      RETURN
C
C=====CAVITIES
C=====12/SEPT /83 ----MAKE IT THE SAME AS A.CHAO'S PAPERS
C=====5/Sept/95 Add in terms detailed in HERA 94-02.
   15 TEMP=YY*(1.+S)
      A65=T(6,5)
      A(1,1)= D2I*D1I*A65    +1.D0
      A(2,1)= D2I*D2I*A65
      A(3,1)= D2I*D3I*A65
      A(4,1)= D2I*D4I*A65
C
      A(1,2)=-D1I*D1I*A65
      A(2,2)=-D1I*D2I*A65    +1.D0
      A(3,2)=-D1I*D3I*A65
      A(4,2)=-D1I*D4I*A65
C
      A(1,3)= D4I*D1I*A65
      A(2,3)= D4I*D2I*A65
      A(3,3)= D4I*D3I*A65    +1.D0
      A(4,3)= D4I*D4I*A65
C
      A(1,4)=-D3I*D1I*A65
      A(2,4)=-D3I*D2I*A65
      A(3,4)=-D3I*D3I*A65
      A(4,4)=-D3I*D4I*A65    +1.D0
C
      A(1,5)=    -D1I*A65
      A(2,5)=    -D2I*A65
      A(3,5)=    -D3I*A65
      A(4,5)=    -D4I*A65
C
      A(6,1)=    -D2I*A65
      A(6,2)=     D1I*A65
      A(6,3)=    -D4I*A65
      A(6,4)=     D3I*A65
      A(6,5)=T(6,5)
      IF(NNSOL.EQ.1)RETURN
      A(7,2)=-TEMP*LZM
      A(7,4)= TEMP*LXM
      A(8,2)= TEMP*MZM
      A(8,4)=-TEMP*MXM
C     IF(NNSOL.EQ.3)RETURN
      A(7,6)= A(7,2)*D2I+A(7,4)*D4I
      A(8,6)= A(8,2)*D2I+A(8,4)*D4I
C     A(7,2)=0.D0
C     A(8,2)=0.D0
      RETURN
C
C=====VERTICAL DIPOLES
   17 CONTINUE
      A(5,6)=0.D0
      IF(XX.NE.0.D0)
     +    A(5,6)=-D3D*XX/YY*A(3,4)+D4D*(A(3,3)-1.D0)*YY/XX-(YY-A(3,4))
      IF(NNSOL.EQ.1)RETURN
      SO=DSIN(+S*XX)
      CO=DCOS(+S*XX)
      AGP1=(1.D0+S)
      A(7,2)=-SO*LSI+(CO-1.D0)*LZI
      A(8,2)=+SO*MSI-(CO-1.D0)*MZI
      A(7,3)=-AGP1* A(4,3)      *LXI
      A(8,3)= AGP1* A(4,3)      *MXI
      A(7,4)=-AGP1*(A(4,4)-1.D0)*LXI
      A(8,4)= AGP1*(A(4,4)-1.D0)*MXI
      A(7,6)= A(7,1)*D1I+A(7,2)*D2I+A(7,3)*D3I+A(7,4)*D4I-XX*LXI
      A(8,6)= A(8,1)*D1I+A(8,2)*D2I+A(8,3)*D3I+A(8,4)*D4I+XX*MXI
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
      RETURN
C
C=====SEXTUPOLE:CODE AS THIN LENS BETWEEN 2 HALF DRIFTS.
   20 CONTINUE
C=====TEST: REPLACE C.O. WITH DISPERSION. VERT DISP IS USUALLY ZERO.
      DELX=ECHROM*D1I
      DELY=ECHROM*D3I
C
      CALL UNIT(8,A)
      CALL UNIT(8,B)
      B(1,2)=YY*0.5D0
      B(3,4)=YY*0.5D0
      A(2,1)=-XX*DELX
      A(2,3)= XX*DELY
C     A(2,6)= XX*(DELX*DELX-DELY*DELY)/2.
      A(4,1)= XX*DELY
      A(4,3)= XX*DELX
C     A(4,6)=-XX*DELX*DELY
C     A(5,1)=-A(2,6)
C     A(5,3)=-A(4,6)
      TEMP=(1.+S)*XX
C     IF(NNSOL.EQ.1)RETURN  TOO MESSY TO USE THIS WITH DRIFT MULTS.
      A(7,1)=-TEMP*(DELY*LXM+DELX*LZM)
      A(7,3)=-TEMP*(DELX*LXM-DELY*LZM)
      A(8,1)= TEMP*(DELY*MXM+DELX*MZM)
      A(8,3)= TEMP*(DELX*MXM-DELY*MZM)
      A(7,6)= A(7,1)*D1I+A(7,2)*D2I+A(7,3)*D3I+A(7,4)*D4I
      A(8,6)= A(8,1)*D1I+A(8,2)*D2I+A(8,3)*D3I+A(8,4)*D4I
      CALL JAM888(A,A,B)
      CALL JAM888(A,B,A)
      RETURN
C
C=====HORIZONTAL C.F. DIPOLE.
   21 CONTINUE
      KX=XX/YY
      G=X2/YY
      G1=KX*KX+G
      A(5,6)=0.D0
      IF(XX.EQ.0.D0.AND.X2.EQ.0.D0)RETURN
      A(5,6)=-D1D*KX*A(1,2)+D2D*KX/G1*(A(1,1)-1.D0)-(YY-A(1,2))*KX*KX/G100004280
      IF(NNSOL.EQ.1)RETURN
      AGP1=(1.D0+S*(1.D0+ECHROM))
C     AGP1=(1.D0+S*(1.D0       ))
      A(7,1)= AGP1* A(2,1)      *LZI
      A(8,1)=-AGP1* A(2,1)      *MZI
      A(7,2)= AGP1*(A(2,2)-1.D0)*LZI
      A(8,2)=-AGP1*(A(2,2)-1.D0)*MZI
      A3=A(4,3)
      A4=A(4,4)
      SO=DSIN(-S*XX)
      CO=DCOS(-S*XX)
      OM=-S*KX
C
C=====+VE OR ZERO G
      IF(G.LT.0.D0)GO TO 31
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
      GO TO 32
   31 CONTINUE
C=====G<0
      G=-X2/YY
      A(7,3)=-AGP1*G*LSI*(  A3*SO-OM*A4*CO+OM)
     +       -AGP1*G*LXI*(  A3*CO+OM*A4*SO   )
     +       +  S*KX*LSI*(G*A4*CO-OM*A3*SO- G)
     +       -  S*KX*LXI*(G*A4*SO+OM*A3*CO   )
      A(7,3)=A(7,3)/(G-OM*OM)
      A(8,3)=-AGP1*G*MSI*(  A3*SO-OM*A4*CO+OM)
     +       -AGP1*G*MXI*(  A3*CO+OM*A4*SO   )
     +       +  S*KX*MSI*(G*A4*CO-OM*A3*SO- G)
     +       -  S*KX*MXI*(G*A4*SO+OM*A3*CO   )
      A(8,3)=-A(8,3)/(G-OM*OM)
      A(7,4)=-AGP1*LSI*  (G*A4*SO+OM*A3*CO   )
     +       -AGP1*LXI*  (G*A4*CO-OM*A3*SO- G)
     +       -S*KX*LSI*  (  A3*CO+OM*A4*SO   )
     +       +S*KX*LXI*  (  A3*SO-OM*A4*CO+OM)
      A(7,4)=A(7,4)/(G-OM*OM)
      A(8,4)=-AGP1*MSI*  (G*A4*SO+OM*A3*CO   )
     +       -AGP1*MXI*  (G*A4*CO-OM*A3*SO- G)
     +       -S*KX*MSI*  (  A3*CO+OM*A4*SO   )
     +       +S*KX*MXI*  (  A3*SO-OM*A4*CO+OM)
      A(8,4)=-A(8,4)/(G-OM*OM)
   32 CONTINUE
      A(7,6)= A(7,1)*D1I+A(7,2)*D2I+A(7,3)*D3I+A(7,4)*D4I+XX*LZI
      A(8,6)= A(8,1)*D1I+A(8,2)*D2I+A(8,3)*D3I+A(8,4)*D4I-XX*MZI
C
C     A(7,1)= 0.D0
C     A(8,1)= 0.D0
C     A(7,2)= 0.D0
C     A(8,2)= 0.D0
C
      RETURN

C====BEAM-BEAM
   22 CONTINUE   
      AGP1 = (1.D0 + S)
      IF(NNSOL.EQ.1)RETURN
      A(7,1)= AGP1* A(2,1)      *LZM
      A(8,1)=-AGP1* A(2,1)      *MZM
      A(7,3)=-AGP1* A(4,3)      *LXM
      A(8,3)= AGP1* A(4,3)      *MXM
      IF(NNSOL.EQ.2)RETURN
      A(7,6)= A(7,1)*D1I+A(7,2)*D2I
      A(8,6)= A(8,1)*D1I+A(8,2)*D2I
      A(7,6)= A(7,6)+A(7,3)*D3I+A(7,4)*D4I
      A(8,6)= A(8,6)+A(8,3)*D3I+A(8,4)*D4I
      RETURN
C
C
      END
