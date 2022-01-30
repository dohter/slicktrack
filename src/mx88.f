C   15/11/79 401171758  MEMBER NAME  MX88     (MAY92.S)     FORTRAN
      SUBROUTINE MX88(ID,II,ITY,T,XX,X2,YY,S,ZWM,ZWI,A)
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C=====GENERATE THE 8X8 MATRIX USING THE SPIN BASIS.
C=====THICK LENS VERSION: PICK UP MATRICES FROM TMAT
C


      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ZWM(3,3),ZWI(3,3),A(8,8),B(8,8)
      DIMENSION T(7,7)
      REAL*8 LXM,LZM,LSM,MXM,MZM,MSM
      REAL*8 LXI,LZI,LSI,MXI,MZI,MSI
      REAL*8 KX
C
C
      INCLUDE "csol.for"
      INCLUDE "cloorb.for"
C
C
C
C
C=====SET UP INITIAL OPTICAL PART OF 8X8 MATRIX: SPIN PART REMAINS EMPTY
C=====UNLESS FILLED LATER. INCLUDE BEAM-BEAM CENTRED ON THE C.O.
C     CALL VZERO(A,128)
C     DO 1 J=1,6
C     DO 1 K=1,6
C   1 A(J,K)=T(J,K)

      AP1 = 1.159652D-3 + 1.0D0

      A(1:6,1:6) = T(1:6,1:6)
      A(7:8,1:8) = 0.D0
      A(1:6,7:8) = 0.D0
      A(7,7)=1.D0
      A(8,8)=1.D0
C 9277 WRITE(53,'(A,1X,I5,F16.6)')' In MX88 8X8 TRANSF.MATRIX at',ID,YY
C      DO 2733 I=1,8
C 2733 WRITE(53,925)(A(I,J),J=1,8)
C  925 FORMAT(8D16.8)


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
      GOTO(11,12,13,14,15,18,19,20,17,11,11,11,11,11,21,17,22),ID
   11 RETURN
C
C=====HORIZONTAL DIPOLES
   12 CONTINUE
      DMUL=1.D0
C     IF(NNSOL.EQ.1)DMUL=0.D0
      A(1,6)= A(1,6)     *DMUL
      A(2,6)= A(2,6)     *DMUL
      A(5,1)= A(5,1)     *DMUL
      A(5,2)= A(5,2)     *DMUL
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
C     IF(NNSOL.EQ.1)RETURN
      A(7,6)=-AGP1*(XX-A(2,6))*LZI+XX*LZI
      A(8,6)= AGP1*(XX-A(2,6))*MZI-XX*MZI
C 9277 WRITE(53,'(A,1X,I5,F16.6)')' In MX88 8X8 TRANSF.MATRIX at',ID,YY
C      DO 273 III=1,8
C  273 WRITE(53,925)(A(III,JJJ),JJJ=1,8)
C  925 FORMAT(' ',8D16.8)
      RETURN
C
C=====QUADRUPOLES: FUDGE IN ROW 5 & COL 6 TERMS AS IN MX66.
C=====ALEX NEVER INCLUDED SPIN TERMS TO CORRESPOND TO THE COL.6 TERMS.
C=====BUT WHICH SPIN BASIS SHOULD THEN BE TAKEN AS UNPERTURBED?
C=====ANYWAY THESE PURE DIPOLE TERMS WOULD BE VERY SMALL
C=====NEED HALF THE EFFECT AT FRONT AND BACK: CORRECTED 8/10/86.
   13 CONTINUE
      CALL UNIT(8,B)
      B(2,6)= XX*DELX*0.5D0
      B(4,6)=-XX*DELY*0.5D0
      B(5,1)=-XX*DELX*0.5D0
      B(5,3)= XX*DELY*0.5D0
      CALL JAM888(A,A,B)
      CALL JAM888(A,B,A)
      AGP1=(1.D0+S)
      FAC=1.0D0
C      WRITE(53,331)LZI,MZI
C  331 FORMAT(' ','lz,mz:  ',2F10.6)
C     IF(NNSOL.EQ.1)FAC=0.9
C     IF(NNSOL.NE.2)RETURN
      A(7,1)= AGP1* A(2,1)      *LZM*FAC
      A(8,1)=-AGP1* A(2,1)      *MZM*FAC


      A(7,2)= AGP1*(A(2,2)-1.D0)*LZM*FAC
      A(8,2)=-AGP1*(A(2,2)-1.D0)*MZM*FAC

      A(7,3)=-AGP1* A(4,3)      *LXM           *1.D0
      A(8,3)= AGP1* A(4,3)      *MXM           *1.D0
      A(7,4)=-AGP1*(A(4,4)-1.D0)*LXM           *1.D0
      A(8,4)= AGP1*(A(4,4)-1.D0)*MXM           *1.D0
      IF(NNSOL.EQ.999)THEN  ! HORIZONTAL EDGE FIELD
C       WRITE(*,*) X2, XX, A(2,1), AP1, LZM
        A(7,3)=-AGP1* A(4,3)*LXM + AP1*X2*LSM
        A(8,3)= AGP1* A(4,3)*MXM - AP1*X2*MSM
C       WRITE(*,*) "THIS IS THE EDGE"
      ELSEIF(NNSOL.EQ.998)THEN  ! VERTICAL EDGE FIELD
C       WRITE(*,*) X2, XX, A(2,1), AP1, LZM
        A(7,1)= AGP1* A(2,1)*LZM - AP1*X2*LSM
        A(8,1)=-AGP1* A(2,1)*MZM + AP1*X2*MSM
C       WRITE(*,*) "THIS IS THE EDGE"
      ELSE  ! REGULAR QUADRUPOLE
C     WRITE(*,*) "THIS IS THE QUAD"
      ENDIF
C 9278 WRITE(53,'(A,1X,I5,F16.6)')' In MX88 8X8 TRANSF.MATRIX at',ID,YY
C      DO 274 III=1,8
C  274 WRITE(53,926)(A(III,JJJ),JJJ=1,8)
C  926 FORMAT(' ',8D16.8)
      RETURN
C
C=====SKEW QUADS.
   14 AGP1=(1.+S)
C=====ALEX NEVER INCLUDED SPIN TERMS TO CORRESPOND TO THE COL.6 TERMS.
C     A(2,6)= XX*DELY
C     A(4,6)= XX*DELX
C     A(5,1)=-XX*DELY
C     A(5,3)=-XX*DELX
      IF(NNSOL.EQ.1)RETURN
C     IF(NNSOL.EQ.1)GO TO 2273
      A(7,1)= AGP1*(A(2,1)*LZM-A(2,3)*LXM)
      A(8,1)=-AGP1*(A(2,1)*MZM-A(2,3)*MXM)
      A(7,2)= AGP1*((A(2,2)-1.D0)*LZM-A(2,4)*LXM)
      A(8,2)=-AGP1*((A(2,2)-1.D0)*MZM-A(2,4)*MXM)
      A(7,3)= AGP1*(A(2,3)*LZM-A(2,1)*LXM)
      A(8,3)=-AGP1*(A(2,3)*MZM-A(2,1)*MXM)
      A(7,4)= AGP1*(A(2,4)*LZM-(A(2,2)-1.D0)*LXM)
      A(8,4)=-AGP1*(A(2,4)*MZM-(A(2,2)-1.D0)*MXM)
C2273 DO 273 III=1,8
C 273 WRITE(53,925)(A(III,JJJ),JJJ=1,6)
C 925 FORMAT(' ','Skew Quad  ',6D16.8)
      RETURN
C
C=====CAVITIES
C=====12/SEPT /83 ----MAKE IT THE SAME AS A.CHAO'S PAPERS
   15 TEMP=YY*(1.+S)
C     A(5,5)=1.D0
C     A(6,6)=1.D0
C     A(6,5)=0.D0
      IF(NNSOL.EQ.1)RETURN
      A(7,2)=-TEMP*LZM
      A(7,4)= TEMP*LXM
      A(8,2)= TEMP*MZM
      A(8,4)=-TEMP*MXM
      RETURN
C
C=====VERTICAL DIPOLES
   17 CONTINUE
C     MXI=1.D0
C     LXI=0.D0
C     MZI=0.D0
C     LZI=-1.D0
C     MSI=0.D0
C     LSI=0.D0
C
C     MXI=0.D0
C     LXI=1.D0
C     MZI=-1.D0
C     LZI=0.D0
C     MSI=0.D0
C     LSI=0.D0
C
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
      A(7,6)= AGP1*(XX-A(4,6))*LXI-XX*LXI
      A(8,6)=-AGP1*(XX-A(4,6))*MXI+XX*MXI
C     IF(NNSOL.NE.30)RETURN
C     WRITE(53,926)
C 926 FORMAT(' ','vert bend')
C     DO 273 III=1,8
C 273 WRITE(53,925)(A(III,JJJ),JJJ=1,8)
C 925 FORMAT(' ',8D16.8)
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
C=====ALEX NEVER INCLUDED SPIN TERMS TO CORRESPOND TO THE COL.6 TERMS.
   20 CONTINUE
      CALL UNIT(8,A)
      CALL UNIT(8,B)
      B(1,2)=YY*0.5D0
      B(3,4)=YY*0.5D0
      A(2,1)=-XX*DELX
      A(2,3)= XX*DELY
      A(2,6)= XX*(DELX*DELX-DELY*DELY)/2.
      A(4,1)= XX*DELY
      A(4,3)= XX*DELX
      A(4,6)=-XX*DELX*DELY
      A(5,1)=-A(2,6)
      A(5,3)=-A(4,6)
      TEMP=(1.+S)*XX
C     IF(NNSOL.EQ.1)RETURN  TOO MESSY TO USE THIS WITH DRIFT MULTS.
      A(7,1)=-TEMP*(DELY*LXM+DELX*LZM)
      A(7,3)=-TEMP*(DELX*LXM-DELY*LZM)
      A(8,1)= TEMP*(DELY*MXM+DELX*MZM)
      A(8,3)= TEMP*(DELX*MXM-DELY*MZM)
      CALL JAM888(A,A,B)
      CALL JAM888(A,B,A)
      RETURN
C
C=====HORIZONTAL C.F. DIPOLE.---ASSUMES A SECTOR SHAPE.
C=====N.B. IN TSPIN THE GRADIENT WAS NEGLECTED IN GETTING THE PERTURBED
C=====SPIN BASIS.
   21 CONTINUE
      IF(XX.EQ.0.D0.AND.X2.EQ.0.D0)RETURN
      IF(NNSOL.EQ.1)RETURN
      AGP1=(1.D0+S)
      A(7,1)= AGP1* A(2,1)      *LZI
      A(8,1)=-AGP1* A(2,1)      *MZI
      A(7,2)= AGP1*(A(2,2)-1.D0)*LZI
      A(8,2)=-AGP1*(A(2,2)-1.D0)*MZI
      A(7,6)=-AGP1*(XX-A(2,6))*LZI+XX*LZI
      A(8,6)= AGP1*(XX-A(2,6))*MZI-XX*MZI
C
      A3=A(4,3)
      A4=A(4,4)
      SO=DSIN(-S*XX)
      CO=DCOS(-S*XX)
      KX=XX/YY
      OM=-S*KX
      G=X2/YY
C
C=====TURN OFF VERT. EFFECT.
C     RETURN
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
      RETURN
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
      RETURN

C====BEAM-BEAM
   22 CONTINUE
      AGP1=(1.D0+S)
      IF(NNSOL.EQ.1)RETURN
      A(7,1)= AGP1* A(2,1)      *LZM
      A(8,1)=-AGP1* A(2,1)      *MZM
      A(7,3)=-AGP1* A(4,3)      *LXM
      A(8,3)= AGP1* A(4,3)      *MXM
C
C
      END
