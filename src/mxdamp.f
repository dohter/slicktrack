C   15/11/79 510201616  MEMBER NAME  MXDAMP   (SEPT95.S)    FORTRAN
      SUBROUTINE MXDAMP(I,IID,ITY,A,CRAD,X,X2,Y,NM)
C
C==========ROUTINE TO CALCULATE DAMPING MATRIX ELEMENTS.================
C
C
C=====USING THIN LENS TERMS: PROBABLY O.K. SINCE IN PRACTICE WE
C=====DIVIDE THE ELEMENTS INTO THICK SLICES.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(6,6)
      CHARACTER *8 NM
C
C
      INCLUDE "cloorb.for"
      INCLUDE "csol.for"
      INCLUDE "cshift.for"
C
C
      DELX=(DX(I)+DX(I+1))/2.D0
      DELY=(DY(I)+DY(I+1))/2.D0
C

C      IF(1.EQ.1)RETURN

C
      GO TO (1,2,3,4,5,2,2,8,2,10,1,1,1,1,11,11,1),IID
    1 RETURN
C
C=====VARIOUS DIPOLES: IGNORING 1/RHO**2 TERMS
    2 CONTINUE
      A(6,6)=A(6,6)-CRAD*X*X*2.D0/Y                    *1.D0
      RETURN
C
C=====QUADRUPOLES
    3 CONTINUE
      IF(NM(1:1).EQ.'E')RETURN
      XY7=CRAD*2.D0*X*X/Y                              *1.D0 
C     IF(NQSH(I).EQ.3)DELY=DY(I)-QSH(I)
C     IF(NQSH(I).EQ.4)DELX=DX(I)-QSH(I)
      A(6,6)=A(6,6)-XY7*(DELX*DELX+DELY*DELY)    
      A(6,1)=-XY7*DELX                           
      A(6,3)=-XY7*DELY                   
      RETURN
C
C
C=====SKEW QUADRUPOLES.
    4 CONTINUE
      XY7=CRAD*2.D0*X*X/Y                              *1.D0 
      A(6,6)=A(6,6)-XY7*(DELX*DELX+DELY*DELY)    
      A(6,1)=-XY7*DELX                           
      A(6,3)=-XY7*DELY                           
      RETURN
C
C=====CAVITIES
    5 A(2,2)=A(2,2)-(Y+X*DL(I))                        *1.D0
      A(4,4)=A(2,2)     
      RETURN
C
C=====SEXTUPOLES 
    8 XY7=CRAD*(DELX*DELX+DELY*DELY)*X*X/Y             *1.D0
      A(6,6)=A(6,6)-XY7*0.5D0*(DELX*DELX+DELY*DELY)
      A(6,1)=-XY7*DELX
      A(6,3)=-XY7*DELY
      RETURN
C
C=====SOLENOIDS
C=====AS IN DECEMBER 1995: AS IN EMITNC TURN OFF RADIATION IN SOLENOIDS.
   10 CONTINUE 
C  10 CALL SOLXYP(I,X/Y,DXPA,DYPA)
C     XY7=CRAD*2.D0*X*X/Y
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
   11 CONTINUE
      A(6,6)=A(6,6)-CRAD*2.D0*X*X/Y                    *1.D0
      A(6,1)=-CRAD*2.D0*X*X2/Y-CRAD*X**3/Y**2          *1.D0
      A(6,3)= 0.D0
C
      RETURN
      END
