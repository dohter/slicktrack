C   06/10/83 510142106  MEMBER NAME  DAMPER   (SEPT95.S)    FORTRAN
      SUBROUTINE PSEQUIL(E0,IE0,U0,CRAD,CIR,TDAMP,TUN)
C
C
C======Set up an equilibrium phase space distribution.
C      This is an modification, to include ``big'' photon emission, of the subroutine DAMPER
C      Track NPART particles, all staring on the C.O. for convenience.
C      We need the correlation matrix but the eigenanalysis is not needed.
C      Emit photons after each dipole or CF and scale the strength according to 
C      the curvature and length. 
C      No radiation in the correction coils.
C
C=====USING THICK LENS TRANSPORT MATRICES BUT THIN LENS DAMPING MATRIX
C=====ELEMENTS.
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C
      DIMENSION VC6A(6),VC6B(6),VC6C(6),DAMP(6),TDAMP(6),TUN(6),DTUN(6)
      DIMENSION TREV6(6,6)
      DIMENSION TRIN6(6,6)
      DIMENSION INTGE6(6)
      DIMENSION TM6A(6,6),TM6B(6,6)
C
C
C
      INCLUDE "cloorb.for"
      INCLUDE "cnlist.for"
      INCLUDE "clatic.for"
      INCLUDE "csol.for"
C
C
C
C
C
C
C
      ANTIDP=1.
      PI=3.1415926535897932D0
C
C
C
      CALL UNIT(6,TREV6)
      DO 20 I=1,NELEM
      ITY=ITYPE(I)
      IID=ID(ITY)
      CALL MX66(TMAT(1,1,ITY),ITY,IID,I,XX(ITY),X2(ITY),YY(ITY),TM6B)
      CALL MXDAMP(I,IID,ITY,TM6B,CRAD,XX(ITY),X2(ITY),YY(ITY))
      CALL JAM666(TREV6,TM6B,TREV6)
   20 CONTINUE
C
C      CALL UCOPY(TREV6,TM6A,72)
       TM6A = TREV6

   94 FORMAT(6F12.5)
C
C
C

C
C
C
C
C
C
C
      RETURN
C
C

C
C
      STOP
      END
