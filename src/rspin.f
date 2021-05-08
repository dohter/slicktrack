C   15/11/79 309261514  MEMBER NAME  RSPIN    (S)           FORTRAN
       SUBROUTINE RSPIN(ROTX,ROTY,ROTZ,ROT)
C
C
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION ROT(3,3)
C
C
       U=DSQRT(ROTX**2+ROTY**2+ROTZ**2)
       IF(U .LT. 1.D-15) GO TO 10
       UX=ROTX/U
       UY=ROTY/U
       UZ=ROTZ/U
       C=DCOS(U)
       S=DSIN(U)
       TEMP=1.D0-C
       ROT(1,1)=UX*UX*TEMP+C
       ROT(1,2)=UX*UY*TEMP-UZ*S
       ROT(1,3)=UX*UZ*TEMP+UY*S
       ROT(2,1)=UX*UY*TEMP+UZ*S
       ROT(2,2)=UY*UY*TEMP+C
       ROT(2,3)=UY*UZ*TEMP-UX*S
       ROT(3,1)=UX*UZ*TEMP-UY*S
       ROT(3,2)=UY*UZ*TEMP+UX*S
       ROT(3,3)=UZ*UZ*TEMP+C
       RETURN
 10    CALL UNIT(3,ROT)
       RETURN
       END
