C   06/05/84 405191653  MEMBER NAME  NQMA8    (DMINS)       FORTRAN
      SUBROUTINE NQMA8(STRE,EL,SPIBAS,GA1,B)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     CALCULATION OF THE  QUADRUPOL 8X8 TRANSFERMATRIX.
C     NO CLOSED ORBIT DEVIATION INCLUDED.
C     THICK LENS VERSION.
C     DERIVED FROM DES BARBER'S NORMQ.
C     PARAMETER STRE=E/PC*(DBX/DZ), EL=LENGTH,SPIBAS=SPIN BASE MATRIX
C                   GA1=SPIN TUNE(GAMMA*A)+1,B=8X8 MATRIX
C     B IS CHANGED ON RETURN.
C
C
C      JOERG KEWISCH MAY 1984
C
C
C
      REAL*8 B(8,8)
      REAL*8 SPIBAS(3,3)
      REAL*8 GA1


C
C
C
C
C
      CALL DMXUTY(B,8)
      B(1,2)=EL
      B(3,4)=EL
      IF(STRE.EQ.0.D0)RETURN
C
C
C
C
      RSTRE=DSQRT(DABS(STRE))
      RSTREL=RSTRE*EL
      DCH=DCOSH(RSTREL)
      DSH=DSINH(RSTREL)
      DC =DCOS (RSTREL)
      DS =DSIN (RSTREL)
      IF(STRE) 10,11,12
   10 B(1,1)= DCH
      B(2,2)= DCH
      B(1,2)= DSH/RSTRE
      B(2,1)= DSH*RSTRE
      B(3,3)= DC
      B(4,4)= DC
      B(3,4)= DS/RSTRE
      B(4,3)=-DS*RSTRE
      GOTO 13
C
   12 B(1,1)= DC
      B(2,2)= DC
      B(1,2)= DS/RSTRE
      B(2,1)=-DS*RSTRE
      B(3,3)= DCH
      B(4,4)= DCH
      B(3,4)= DSH/RSTRE
      B(4,3)= DSH*RSTRE
C
   13 CONTINUE
C
C     SET UP 2X6 PART USING GENERALISED 'COUPLED' FORM.
      B(7,1)=  GA1*SPIBAS(2,3)* B(2,1)
      B(7,2)=  GA1*SPIBAS(2,3)*(B(2,2)-1.D0)
      B(7,3)= -GA1*SPIBAS(1,3)* B(4,3)
      B(7,4)= -GA1*SPIBAS(1,3)*(B(4,4)-1.D0)
      B(8,1)= -GA1*SPIBAS(2,2)* B(2,1)
      B(8,2)= -GA1*SPIBAS(2,2)*(B(2,2)-1.D0)
      B(8,3)=  GA1*SPIBAS(1,2)* B(4,3)
      B(8,4)=  GA1*SPIBAS(1,2)*(B(4,4)-1.D0)
C
C
C
   11 RETURN
      END
