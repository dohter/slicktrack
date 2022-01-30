C   15/11/79 401171758  MEMBER NAME  MX88     (MAY92.S)     FORTRAN
      SUBROUTINE MX99DAMP(A,ITY,ID,II,XX,X2,YY,S,ZWM,ZWI,CRAD,
     +                                                NM,IDAMPFLG)



C=====GENERATE THE THICK LENS 9X9 MATRIX USING THE n0,m,l SPIN BASIS.
C     The bottom 3 rows are the normal G(2x6) matrix augmented
C     with an extra row for needed to get higher order resonances
C     with M-C tracking.
C     Since n0 obeys the T-BMT just in the same way as m and l
C     the 9th row trivially mirrors the 7th and 8th rows.
C     No need to be careful with the sign of row 9. Just know what it is!
C     Choose the sign to be the same as for l.
C     The signs of rows 7 and 8 are(!) important for the G matrix.
C=====PICK UP SYMPLECTIC MATRICES FROM TMAT. INCLUDE DAMPING AFTER EACH ELEMENT.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ZWM(3,3),ZWI(3,3),A(9,9),B(9,9)
      REAL*8 LXM,LZM,LSM, NXM,NZM,NSM, MXM,MZM,MSM
      REAL*8 LXI,LZI,LSI, NXI,NZI,NSI, MXI,MZI,MSI
      REAL*8 KX
      CHARACTER *8 NM
C
C
      INCLUDE "cnlist.for"
      INCLUDE "csol.for"
      INCLUDE "cloorb.for"
C
C
      AP1 = 1.159652D-3 + 1.0D0

      A(7:9,1:9) = 0.D0
      A(1:6,7:9) = 0.D0
      A(7,7)=1.D0
      A(8,8)=1.D0
      A(9,9)=1.D0
C
C
C=====USE AVERAGE CLOSED ORBIT SHIFTS.
      DELX =(DX(II) +DX(II+1))/2.
      DELY =(DY(II) +DY(II+1))/2.
      DELE =(DEL(II)+DEL(II+1))/2.
C
C
C=====SET SPIN BASIS VECTORS --- AVERAGED INPUT & OUTPUT VALUES.
      LXM=ZWM(1,3)
      LZM=ZWM(2,3)
      LSM=ZWM(3,3)

      NXM=ZWM(1,1)
      NZM=ZWM(2,1)
      NSM=ZWM(3,1)

      MXM=ZWM(1,2)
      MZM=ZWM(2,2)
      MSM=ZWM(3,2)
C=====SET SPIN BASIS VECTORS---INPUT VALUES.
      LXI=ZWI(1,3)
      LZI=ZWI(2,3)
      LSI=ZWI(3,3)

      NXI=ZWI(1,1)
      NZI=ZWI(2,1)
      NSI=ZWI(3,1)

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
      A(1,6)= A(1,6)
      A(2,6)= A(2,6)
      A(5,1)= A(5,1)
      A(5,2)= A(5,2)
      IF(NNSOL.EQ.1)GO TO 112
      SO=DSIN(-S*XX)
      CO=DCOS(-S*XX)
      AGP1=(1.D0+S)
      A(7,1)= AGP1* A(2,1)      *LZI
      A(8,1)=-AGP1* A(2,1)      *MZI
      A(9,1)= AGP1* A(2,1)      *NZI

      A(7,2)= AGP1*(A(2,2)-1.D0)*LZI
      A(8,2)=-AGP1*(A(2,2)-1.D0)*MZI
      A(9,2)= AGP1*(A(2,2)-1.D0)*NZI

      A(7,4)=-SO*LSI-(CO-1.D0)*LXI
      A(8,4)= SO*MSI+(CO-1.D0)*MXI
      A(9,4)=-SO*MSI+(CO-1.D0)*NXI
C     IF(NNSOL.EQ.1)RETURN
      A(7,6)=-AGP1*(XX-A(2,6))*LZI+XX*LZI
      A(8,6)= AGP1*(XX-A(2,6))*MZI-XX*MZI
      A(9,6)=-AGP1*(XX-A(2,6))*NZI+XX*NZI
  112 CONTINUE
C======Damping: ignoring 1/RHO**2 terms
      A(6,6)=A(6,6)-CRAD*XX*XX*2.D0/YY            *NDAMP3*IDAMPFLG
      RETURN
C
C=====QUADRUPOLES: FUDGE IN ROW 5 & COL 6 TERMS AS IN MX66.
C=====ALEX NEVER INCLUDED SPIN TERMS TO CORRESPOND TO THE COL.6 TERMS.
C=====BUT WHICH SPIN BASIS SHOULD THEN BE TAKEN AS UNPERTURBED?
C=====ANYWAY THESE PURE DIPOLE TERMS WOULD BE VERY SMALL
C=====NEED HALF THE EFFECT AT FRONT AND BACK: CORRECTED 8/10/86.
   13 CONTINUE
      CALL UNIT(9,B)
      B(2,6)= XX*DELX*0.5D0
      B(4,6)=-XX*DELY*0.5D0
      B(5,1)=-XX*DELX*0.5D0
      B(5,3)= XX*DELY*0.5D0
      A = MATMUL(A,B)
      A = MATMUL(B,A)
      AGP1=(1.D0+S)
      IF(NNSOL.EQ.1)GO TO 113
C     IF(NNSOL.EQ.1)FAC=0.9
      A(7,1)= AGP1* A(2,1)      *LZM
      A(8,1)=-AGP1* A(2,1)      *MZM
      A(9,1)= AGP1* A(2,1)      *NZM

      A(7,2)= AGP1*(A(2,2)-1.D0)*LZM
      A(8,2)=-AGP1*(A(2,2)-1.D0)*MZM
      A(9,2)= AGP1*(A(2,2)-1.D0)*NZM
C     IF(NNSOL.NE.2)RETURN
      A(7,3)=-AGP1* A(4,3)      *LXM
      A(8,3)= AGP1* A(4,3)      *MXM
      A(9,3)=-AGP1* A(4,3)      *NXM

      A(7,4)=-AGP1*(A(4,4)-1.D0)*LXM
      A(8,4)= AGP1*(A(4,4)-1.D0)*MXM
      A(9,4)=-AGP1*(A(4,4)-1.D0)*NXM
  113 CONTINUE
      IF(NNSOL.EQ.999)THEN  ! HORIZONTAL EDGE FIELD
C       WRITE(*,*) X2, XX, A(2,1), AP1, LZM
        A(7,3)=-AGP1* A(4,3)*LXM + AP1*X2*LSM
        A(8,3)= AGP1* A(4,3)*MXM - AP1*X2*MSM
        A(9,3)=-AGP1* A(4,3)*NXM + AP1*X2*NSM
C       WRITE(*,*) "THIS IS THE EDGE"
        RETURN
      ELSEIF(NNSOL.EQ.998)THEN  ! VERTICAL EDGE FIELD
C       WRITE(*,*) X2, XX, A(2,1), AP1, LZM
        A(7,1)= AGP1* A(2,1)*LZM - AP1*X2*LSM
        A(8,1)=-AGP1* A(2,1)*MZM + AP1*X2*MSM
        A(9,1)= AGP1* A(2,1)*NZM - AP1*X2*NSM
C       WRITE(*,*) "THIS IS THE EDGE"
        RETURN
      ELSE  ! REGULAR QUADRUPOLE
C     WRITE(*,*) "THIS IS THE QUAD"
      ENDIF
      XY7=CRAD*2.D0*XX*XX/YY                      *NDAMP3*IDAMPFLG
      A(6,6)=A(6,6)-XY7*(DELX*DELX+DELY*DELY)
      A(6,1)=-XY7*DELX
      A(6,3)=-XY7*DELY
      RETURN
C
C=====SKEW QUADS.
   14 AGP1=(1.+S)
C=====ALEX NEVER INCLUDED SPIN TERMS TO CORRESPOND TO THE COL.6 TERMS.
C     A(2,6)= XX*DELY
C     A(4,6)= XX*DELX
C     A(5,1)=-XX*DELY
C     A(5,3)=-XX*DELX
      IF(NNSOL.EQ.1)GO TO 114
      A(7,1)= AGP1*(A(2,1)*LZM-A(2,3)*LXM)
      A(8,1)=-AGP1*(A(2,1)*MZM-A(2,3)*MXM)
      A(9,1)= AGP1*(A(2,1)*NZM-A(2,3)*NXM)

      A(7,2)= AGP1*((A(2,2)-1.D0)*LZM-A(2,4)*LXM)
      A(8,2)=-AGP1*((A(2,2)-1.D0)*MZM-A(2,4)*MXM)
      A(9,2)= AGP1*((A(2,2)-1.D0)*NZM-A(2,4)*NXM)

      A(7,3)= AGP1*(A(2,3)*LZM-A(2,1)*LXM)
      A(8,3)=-AGP1*(A(2,3)*MZM-A(2,1)*MXM)
      A(9,3)= AGP1*(A(2,3)*NZM-A(2,1)*NXM)

      A(7,4)= AGP1*(A(2,4)*LZM-(A(2,2)-1.D0)*LXM)
      A(8,4)=-AGP1*(A(2,4)*MZM-(A(2,2)-1.D0)*MXM)
      A(9,4)= AGP1*(A(2,4)*NZM-(A(2,2)-1.D0)*NXM)
  114 CONTINUE
      XY7=CRAD*2.D0*XX*XX/YY                      *NDAMP3*IDAMPFLG
      A(6,6)=A(6,6)-XY7*(DELX*DELX+DELY*DELY)
      A(6,1)=-XY7*DELX
      A(6,3)=-XY7*DELY
      RETURN
C
C=====CAVITIES
C=====12/SEPT /83 ----MAKE IT THE SAME AS A.CHAO'S PAPERS
   15 TEMP=YY*(1.D0+S)
C     A(5,5)=1.D0
C     A(6,6)=1.D0
C     A(6,5)=0.D0
      IF(NNSOL.EQ.1)GO TO 115
      A(7,2)=-TEMP*LZM
      A(8,2)= TEMP*MZM
      A(9,2)=-TEMP*NZM

      A(7,4)= TEMP*LXM
      A(8,4)=-TEMP*MXM
      A(9,4)= TEMP*NXM
  115 CONTINUE
      A(2,2)=A(2,2)-(YY+XX*DELE)                  *NDAMP3*IDAMPFLG
      A(4,4)=A(2,2)
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
      IF(NNSOL.EQ.1)GO TO 117
      SO=DSIN(+S*XX)
      CO=DCOS(+S*XX)
      AGP1=(1.D0+S)
      A(7,2)=-SO*LSI+(CO-1.D0)*LZI
      A(8,2)=+SO*MSI-(CO-1.D0)*MZI
      A(9,2)=-SO*LSI+(CO-1.D0)*NZI

      A(7,3)=-AGP1* A(4,3)      *LXI
      A(8,3)= AGP1* A(4,3)      *MXI
      A(9,3)=-AGP1* A(4,3)      *NXI

      A(7,4)=-AGP1*(A(4,4)-1.D0)*LXI
      A(8,4)= AGP1*(A(4,4)-1.D0)*MXI
      A(9,4)=-AGP1*(A(4,4)-1.D0)*NXI

      A(7,6)= AGP1*(XX-A(4,6))*LXI-XX*LXI
      A(8,6)=-AGP1*(XX-A(4,6))*MXI+XX*MXI
      A(9,6)= AGP1*(XX-A(4,6))*NXI-XX*NXI
  117 CONTINUE
C======Damping: ignoring 1/RHO**2 terms
      A(6,6)=A(6,6)-CRAD*XX*XX*2.D0/YY            *NDAMP3*IDAMPFLG
      RETURN
C
C=====HORIZONTAL KICKER.
   18 CONTINUE
      IF(NNSOL.EQ.1)GO TO 118
      A(7,4)= S*XX*LSM
      A(8,4)=-S*XX*MSM
      A(9,4)= S*XX*NSM

      A(7,6)=   XX*LZM
      A(8,6)=  -XX*MZM
      A(9,6)=   XX*NZM
  118 CONTINUE
C======Damping:
      A(6,6)=A(6,6)-CRAD*XX*XX*2.D0/YY            *NDAMP3*IDAMPFLG
      RETURN
C
C=====VERTICAL KICKER.
   19 CONTINUE
      IF(NNSOL.EQ.1)GO TO 119
      A(7,2)= XX*S*LSM
      A(8,2)=-XX*S*MSM
      A(9,2)= XX*S*NSM

      A(7,6)=   XX*LXM
      A(8,6)=  -XX*MXM
      A(9,6)=   XX*NXM
  119 CONTINUE
C======Damping:
      A(6,6)=A(6,6)-CRAD*XX*XX*2.D0/YY            *NDAMP3*IDAMPFLG
      RETURN
C
C=====SEXTUPOLE:CODE AS THIN LENS BETWEEN 2 HALF DRIFTS.
C=====ALEX NEVER INCLUDED SPIN TERMS TO CORRESPOND TO THE COL.6 TERMS.
   20 CONTINUE
      CALL UNIT(9,A)
      CALL UNIT(9,B)
      B(1,2)= YY*0.5D0
      B(3,4)= YY*0.5D0
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
      A(8,1)= TEMP*(DELY*MXM+DELX*MZM)
      A(9,1)=-TEMP*(DELY*NXM+DELX*NZM)

      A(7,3)=-TEMP*(DELX*LXM-DELY*LZM)
      A(8,3)= TEMP*(DELX*MXM-DELY*MZM)
      A(9,3)=-TEMP*(DELX*NXM-DELY*NZM)

      A = MATMUL(A,B)
      A = MATMUL(B,A)
      XY7=CRAD*(DELX*DELX+DELY*DELY)*XX*XX/YY     *NDAMP3*IDAMPFLG
      A(6,6)=A(6,6)-XY7*0.5D0*(DELX*DELX+DELY*DELY)
      A(6,1)=-XY7*DELX
      A(6,3)=-XY7*DELY
      RETURN
C
C=====HORIZONTAL C.F. DIPOLE.---ASSUMES A SECTOR SHAPE.
C=====N.B. IN TSPIN THE GRADIENT WAS NEGLECTED IN GETTING THE PERTURBED
C=====SPIN BASIS.
   21 CONTINUE
      IF(XX.EQ.0.D0.AND.X2.EQ.0.D0)RETURN
      IF(NNSOL.EQ.1)GO TO 221
      AGP1=(1.D0+S)
      A(7,1)= AGP1* A(2,1)      *LZI
      A(8,1)=-AGP1* A(2,1)      *MZI
      A(9,1)= AGP1* A(2,1)      *NZI

      A(7,2)= AGP1*(A(2,2)-1.D0)*LZI
      A(8,2)=-AGP1*(A(2,2)-1.D0)*MZI
      A(9,2)= AGP1*(A(2,2)-1.D0)*NZI

      A(7,6)=-AGP1*(XX-A(2,6))*LZI+XX*LZI
      A(8,6)= AGP1*(XX-A(2,6))*MZI-XX*MZI
      A(9,6)=-AGP1*(XX-A(2,6))*NZI+XX*NZI
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
      IF(G.GE.0.D0)THEN
      A(7,3)=-AGP1*G*LSI*(  A3*SO-OM*A4*CO+OM)
     +       -AGP1*G*LXI*(  A3*CO+OM*A4*SO   )
     +       +  S*KX*LSI*(G*A4*CO+OM*A3*SO- G)
     +       -  S*KX*LXI*(G*A4*SO-OM*A3*CO   )
      A(7,3)= A(7,3)/(OM*OM+G)

      A(8,3)=-AGP1*G*MSI*(  A3*SO-OM*A4*CO+OM)
     +       -AGP1*G*MXI*(  A3*CO+OM*A4*SO   )
     +       +  S*KX*MSI*(G*A4*CO+OM*A3*SO- G)
     +       -  S*KX*MXI*(G*A4*SO-OM*A3*CO   )
      A(8,3)=-A(8,3)/(OM*OM+G)

      A(9,3)=-AGP1*G*NSI*(  A3*SO-OM*A4*CO+OM)
     +       -AGP1*G*NXI*(  A3*CO+OM*A4*SO   )
     +       +  S*KX*NSI*(G*A4*CO+OM*A3*SO- G)
     +       -  S*KX*NXI*(G*A4*SO-OM*A3*CO   )
      A(9,3)= A(9,3)/(OM*OM+G)

      A(7,4)=-AGP1*LSI*  (G*A4*SO-OM*A3*CO   )
     +       -AGP1*LXI*  (G*A4*CO+OM*A3*SO- G)
     +       +S*KX*LSI*  (  A3*CO+OM*A4*SO   )
     +       -S*KX*LXI*  (  A3*SO-OM*A4*CO+OM)
      A(7,4)= A(7,4)/(OM*OM+G)

      A(8,4)=-AGP1*MSI*  (G*A4*SO-OM*A3*CO   )
     +       -AGP1*MXI*  (G*A4*CO+OM*A3*SO- G)
     +       +S*KX*MSI*  (  A3*CO+OM*A4*SO   )
     +       -S*KX*MXI*  (  A3*SO-OM*A4*CO+OM)
      A(8,4)=-A(8,4)/(OM*OM+G)

      A(9,4)=-AGP1*NSI*  (G*A4*SO-OM*A3*CO   )
     +       -AGP1*NXI*  (G*A4*CO+OM*A3*SO- G)
     +       +S*KX*NSI*  (  A3*CO+OM*A4*SO   )
     +       -S*KX*NXI*  (  A3*SO-OM*A4*CO+OM)
      A(9,4)= A(9,4)/(OM*OM+G)

      ELSE
C=====G<0
      G=-X2/YY
      A(7,3)=-AGP1*G*LSI*(  A3*SO-OM*A4*CO+OM)
     +       -AGP1*G*LXI*(  A3*CO+OM*A4*SO   )
     +       +  S*KX*LSI*(G*A4*CO-OM*A3*SO- G)
     +       -  S*KX*LXI*(G*A4*SO+OM*A3*CO   )
      A(7,3)= A(7,3)/(G-OM*OM)

      A(8,3)=-AGP1*G*MSI*(  A3*SO-OM*A4*CO+OM)
     +       -AGP1*G*MXI*(  A3*CO+OM*A4*SO   )
     +       +  S*KX*MSI*(G*A4*CO-OM*A3*SO- G)
     +       -  S*KX*MXI*(G*A4*SO+OM*A3*CO   )
      A(8,3)=-A(8,3)/(G-OM*OM)

      A(9,3)=-AGP1*G*NSI*(  A3*SO-OM*A4*CO+OM)
     +       -AGP1*G*NXI*(  A3*CO+OM*A4*SO   )
     +       +  S*KX*NSI*(G*A4*CO-OM*A3*SO- G)
     +       -  S*KX*NXI*(G*A4*SO+OM*A3*CO   )
      A(9,3)= A(9,3)/(G-OM*OM)


      A(7,4)=-AGP1*LSI*  (G*A4*SO+OM*A3*CO   )
     +       -AGP1*LXI*  (G*A4*CO-OM*A3*SO- G)
     +       -S*KX*LSI*  (  A3*CO+OM*A4*SO   )
     +       +S*KX*LXI*  (  A3*SO-OM*A4*CO+OM)
      A(7,4)= A(7,4)/(G-OM*OM)

      A(8,4)=-AGP1*MSI*  (G*A4*SO+OM*A3*CO   )
     +       -AGP1*MXI*  (G*A4*CO-OM*A3*SO- G)
     +       -S*KX*MSI*  (  A3*CO+OM*A4*SO   )
     +       +S*KX*MXI*  (  A3*SO-OM*A4*CO+OM)
      A(8,4)=-A(8,4)/(G-OM*OM)

      A(9,4)=-AGP1*NSI*  (G*A4*SO+OM*A3*CO   )
     +       -AGP1*NXI*  (G*A4*CO-OM*A3*SO- G)
     +       -S*KX*NSI*  (  A3*CO+OM*A4*SO   )
     +       +S*KX*NXI*  (  A3*SO-OM*A4*CO+OM)
      A(9,4)= A(9,4)/(G-OM*OM)

      ENDIF
C
  221 CONTINUE
C=====NOTE X*X2 CROSS TERMS WHICH CAN GIVE ANTIDAMPING.

      A(6,6)=A(6,6)-CRAD*2.D0*XX*XX/YY             *NDAMP3*IDAMPFLG
      A(6,1)=(-CRAD*2.D0*XX*X2/YY-CRAD*XX**3/YY**2)*NDAMP3*IDAMPFLG
      A(6,3)= 0.D0                           !No ``vertical'' cross term.
      RETURN
C
C
      END
