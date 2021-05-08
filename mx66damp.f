C   15/11/79 411161900  MEMBER NAME  MX66     (MAY92.S)     FORTRAN
      SUBROUTINE MX66DAMP(A,ITY,IID,II,XX,X2,YY,CRAD,NM,IDAMPFLG)
C
C=====GET 6X6 TERMS IN PRESENCE OF C.O. SHIFTS. 
C     To gain speed (just one CALL), include damping from mxdamp.f,  
C     i.e. using thin lens terms: 
C     Probably OK if we divide the elements into slices.
C     The test is in whether the damping comes out correct during tracking.
C     Feed in the matrix A and send it out again, modified. 
C
C=====PUT IN SEXTUPOLE  CURVATURE TERMS AS THIN LENS
C=====AT THE LENS CENTRE SANDWICHED BETWEEN HALF LENGTH DRIFTS?
C     SEXEQ,MX66,MX88,MXDAMP
C
C     Only include LINEAR beam-beam here. To check effect of
C     nonlinear kicks on the beam, use SCRURITA2(3).
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(6,6),B(6,6),C(6,6)
      CHARACTER *8 NM
C
C
      INCLUDE "cnlist.for"
      INCLUDE "cloorb.for"
      INCLUDE "csol.for"
C
C
C      I2ND=0
C      WRITE(53,'(A,A,F20.8)')' ','CRAD',CRAD

C
      DELX= (DX(II) +DX(II+1))/2.D0
      DELY= (DY(II) +DY(II+1))/2.D0
      DELE =(DEL(II)+DEL(II+1))/2.D0
      DELL =(DL(II) +DL(II+1))/2.D0
C
C
      GO TO (1,2,3,4,5,2,2,8,2,10,1,1,1,1,11,11,1),IID

C      Drift, beam-beam and various exotics.
    1 RETURN

C=====VARIOUS DIPOLES: IGNORING 1/RHO**2 TERMS
    2 A(6,6)=A(6,6)-CRAD*XX*XX*2.D0/YY            *NDAMP1*IDAMPFLG
      RETURN
C


C=====QUADRUPOLE: WITH C.O. SHIFT FUDGE THE CURVATURE STUFF BY PUTTING
C=====IN DISPERSION GENERATING KICKERS OF 1/2 STRENGTH FORE & AFT.
C=====SINCE WE ALREADY HAVE SLICES THAT SHOULD NOT BE TOO BAD.
    3 CONTINUE
      CALL UNIT(6,B)
      B(2,6)= XX*DELX/2.D0    ! *(1.D0-DELE*I2ND)  *1.D0
      B(4,6)=-XX*DELY/2.D0    ! *(1.D0-DELE*I2ND)  *1.D0
      B(5,1)=-XX*DELX/2.D0    ! *(1.D0-DELE*I2ND)  *1.D0
      B(5,3)= XX*DELY/2.D0    ! *(1.D0-DELE*I2ND)  *1.D0
C      CALL JAM666(A,A,B)
C      CALL JAM666(A,B,A)
      A = MATMUL(A,B)
      A = MATMUL(B,A)
      IF(NM(1:1).EQ.'E')RETURN     !Kill radiation in edge fields. 
      XY7=CRAD*2.*XX*XX/YY                        *NDAMP1*IDAMPFLG 
      A(6,6)=A(6,6)-XY7*(DELX*DELX+DELY*DELY)    
      A(6,1)=-XY7*DELX                           
      A(6,3)=-XY7*DELY   
      RETURN
C
C=====SKEW QUADRUPOLE
C=====THESE ALSO NEED THE DISPERSION STUFF.
  4   CONTINUE
C     A(2,6)= XX*DELY
C     A(4,6)= XX*DELX
C     A(5,1)=-XX*DELY
C     A(5,3)=-XX*DELX
      XY7=CRAD*2.*XX*XX/YY                        *NDAMP1*IDAMPFLG 
      A(6,6)=A(6,6)-XY7*(DELX*DELX+DELY*DELY)    
      A(6,1)=-XY7*DELX                           
      A(6,3)=-XY7*DELY                           
      RETURN

C=====CAVITIES
    5 A(2,2)=A(2,2)-(YY+XX*DELL)                  *NDAMP1*IDAMPFLG
      A(4,4)=A(2,2)     
      RETURN


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
    8 CONTINUE
      CALL UNIT(6,B)
      CALL UNIT(6,C)
      C(1,2)=YY*0.5D0
      C(3,4)=YY*0.5D0
      XS=XX*(1.D0)                      !-DELE*I2ND
      B(2,1)=-XS*DELX
      B(2,3)= XS*DELY
      B(2,6)= XS*(DELX**2-DELY**2)/2.D0
      B(4,1)= B(2,3)
      B(4,3)=-B(2,1)
      B(4,6)=-XS*DELX*DELY
      B(5,1)=-B(2,6)
      B(5,3)=-B(4,6)
C      CALL JAM666(A,B,C)
C      CALL JAM666(A,C,A)
      A = MATMUL(B,C) 
      A = MATMUL(C,A)
      XY7=CRAD*(DELX*DELX+DELY*DELY)*XX*XX/YY     *NDAMP1*IDAMPFLG
      A(6,6)=A(6,6)-XY7*0.5D0*(DELX*DELX+DELY*DELY)
      A(6,1)=-XY7*DELX
      A(6,3)=-XY7*DELY
      RETURN

   10 CONTINUE  
C     Solenoids: skip the damping for now.
C     CALL SOLXYP(I,X/Y,DXPA,DYPA)
C     XY7=CRAD*2.*X*X/Y
C     A(6,6)=A(6,6)-XY7*(DXPA*DXPA+DYPA*DYPA)
C     A(6,2)=-XY7*DXPA
C     A(6,4)=-XY7*DYPA
C     A(6,6)=0.D0
      A(6,2)=0.D0
      A(6,4)=0.D0
      RETURN

C
C=====C.F. DIPOLES.----NO CLOSED ORBIT STUFF YET
C=====NOTE X*X2 CROSS TERMS WHICH CAN GIVE ANTIDAMPING.
   11 A(6,6)=A(6,6)-CRAD*2.D0*XX*XX/YY             *NDAMP1*IDAMPFLG
      A(6,1)=(-CRAD*2.D0*XX*X2/YY-CRAD*XX**3/YY**2)*NDAMP1*IDAMPFLG
      A(6,3)= 0.D0
C
      END


