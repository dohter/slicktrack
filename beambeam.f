C   15/11/79 401171758  MEMBER NAME  MX88     (MAY92.S)     FORTRAN
C
      SUBROUTINE BEAMBEAM(A,ITY,ID,II,XX,X2,YY,S,ZWI)


C=====GENERATE THE THIN LENS BEAM-BEAM KICKS USING THE SPIN BASIS.
C     PICK UP SYMPLECTIC MATRICES FROM TMAT. 
C     If |IBMBM| = 1, just use the prcalculated orbit matrix.


C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ZWI(3,3),A(8,8)
      REAL*8 LXI,LZI,LSI,MXI,MZI,MSI
C
C
      INCLUDE "cnlist.for" 
      INCLUDE "csol.for"
      INCLUDE "cloorb.for"
C
C
      A(7:8,1:8) = 0.D0
      A(1:6,7:8) = 0.D0
      A(7,7)=1.D0
      A(8,8)=1.D0
C

C=====Assume that the beams are relatively centred.
C=====USE AVERAGE CLOSED ORBIT SHIFTS.
C      DELX =(DX(II) +DX(II+1))/2.
C      DELY =(DY(II) +DY(II+1))/2.
C      DELE =(DEL(II)+DEL(II+1))/2.
C
C

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
      AGP1=(1.D0+S)
      NNSOL=NSOL(ITY)
      IF(NNSOL.EQ.1)RETURN
      A(7,1)= AGP1* A(2,1)      *LZI
      A(8,1)=-AGP1* A(2,1)      *MZI
      A(7,3)=-AGP1* A(4,3)      *LXI     
      A(8,3)= AGP1* A(4,3)      *MXI   
      RETURN
      END
