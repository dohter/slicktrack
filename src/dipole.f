C   25/07/82 406051707  MEMBER NAME  DIPOLE   (S)           FORTRAN
      SUBROUTINE DIPOLE(IN,NTY,TOTL,NU,PL)
C
C    *GENERATES THE DETAILED STRUCTURE OF A FOCUSSING DIPOLE. *
C    *                                                        *
C    *                                                        *
C    * A FOCUSSING DIPOLE IS GIVEN THE STRUCTURE:             *
C    *                                                        *
C    * HALFQUAD--DIPOLE--QUAD--DIPOLE--QUAD--DIPOLE--HALFQUAD *
C    *                                                        *
C    * AND ALL DRIFT SPACES ETC ARE SET UP.                   *
C    **********************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER DNAME/'D'/
C
C
#include "clatic.for"
#include "csol.for"
C
C
C
C=====SET UP LAST DRIFT SPACE AND SINCE IT IS NEW,STORE IT.
      PL=PL-YY(IN)/2.
      TOTL=TOTL+PL
      NTY=NTY+1
      IF(NTY.GT.400) GO TO 9
      SNAME(NTY)=NAME(IN)
      NAME(NTY)= DNAME
      ID(NTY)=1
      XX(NTY)=PL
      YY(NTY)=0.
      TMAT(1,2,NTY)= PL
      TMAT(3,4,NTY)= PL
      NELEM=NELEM+1
      IF(NELEM.GT.19099)GO TO 12
      ITYPE(NELEM)=NTY
C
C=====STORE DETAILS OF FIRST QUAD.--DETAILS TO BE FOUND 2 ENTRIES AFTER
C=====THE DIPOLE IN THE TYPE LIST.
      NELEM=NELEM+1
      IF(NELEM.GT.19099)GO TO 12
      ITYPE(NELEM)=IN+2

C=====SET UP QUAD-DIP.DRIFT SPACE AND SINCE IT IS NEW,STORE IT.
      PL=YY(IN)/2.
      NTY=NTY+1
      IF(NTY.GT.400) GO TO 9
      SNAME(NTY)=NAME(IN)
      NAME(NTY)= DNAME
      ID(NTY)=1
      XX(NTY)=PL
      YY(NTY)=0.
      TMAT(1,2,NTY)= PL
      TMAT(3,4,NTY)= PL
      NELEM=NELEM+1
      IF(NELEM.GT.19099)GO TO 12
      ITYPE(NELEM)=NTY
C
C=====LOOP TO SET UP THE DIPOLE-QUAD SANDWICH IN THE TYPE LIST.
      DO 10 J=1,NU
C=====FIRST THE DIPOLE SLICE
      NELEM=NELEM+1
      IF(NELEM.GT.19099)GO TO 12
      ITYPE(NELEM)=IN
C
C=====THEN THE DRIFT
      NELEM=NELEM+1
      IF(NELEM.GT.19099)GO TO 12
      ITYPE(NELEM)=NTY
C
C=====THEN THE QUAD
      NELEM=NELEM+1
      IF(NELEM.GT.19099)GO TO 12
      ITYPE(NELEM)=IN+1
      IF(J.EQ.NU)ITYPE(NELEM)=IN+2
C
C=====THEN THE DRIFT
      IF(J.EQ.NU)GO TO 10
      NELEM=NELEM+1
      IF(NELEM.GT.19099)GO TO 12
      ITYPE(NELEM)=NTY
   10 CONTINUE
C
      TOTL=TOTL+YY(IN)*NU
C
C
      RETURN
   12 WRITE(6,95)  TOTL
   95 FORMAT(' SIZE OF LATTICE LIST EXCEEDED !!!',/' LTOT =',F18.4)
      STOP
    9 WRITE(6,94)
   94 FORMAT(' SIZE OF TYPELISTE EXCEEDED !!!')
      STOP
      END
