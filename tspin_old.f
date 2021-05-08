C   15/11/79 511212038  MEMBER NAME  TSPIN    (SEPT95.S) M  FORTRAN
      SUBROUTINE TSPIN(II,IID,ITY,XXX,XX2,YYY,S,C)
C
C
C=====ROUTINE TO ROTATE THE SPIN BASIS AS IT TRAVERSES A MAGNET.
C
C=====IN C.F. MAGNETS WITH C.O. SHIFT,
C=====NO 'C.O.*GRADIENT' PART IS ADDED YET
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(3,3),A(3,3),B(3,3)
C
C
      INCLUDE "csol.for"
      INCLUDE "cloorb.for"
      INCLUDE "clatic.for"
C
C
C
      PI=3.1415926535897932D0
C
C
      IF(IID.EQ.1.OR.IID.EQ.99)GO TO 1
      DDX =(DX(II)+DX(II+1))/2.D0       *1.D0     !Sept 2003 kill effect of C.O. distn.
      DDY =(DY(II)+DY(II+1))/2.D0       *1.D0     !Also SOL8AN
      XXX1=XXX*(S+1.D0)                 *1.D0     !Kill too for kickers.
      XXX2=XXX*S
      GG=XXX*(DEL(II)+DEL(II+1))/2.D0   *1.D0
      GGX=(DXP(II)+DXP(II+1))/2.D0      *1.D0
      GGY=(DYP(II)+DYP(II+1))/2.D0      *1.D0
      GOTO (1,2,3,4,5,6,7,8,9,10,1,1,1,11,15,9),IID
C
    1 CALL UNIT(3,C)
      RETURN
    2 CALL RSPIN(0.D0,-XXX2+GG,XXX2*GGY,C)
      RETURN
    3 CONTINUE
C     IF(NTWIST(ITY).EQ.3)DDY=DY(II)-TWIST(ITY)
C     IF(NTWIST(ITY).EQ.4)DDX=DX(II)-TWIST(ITY)
      CALL RSPIN(-XXX1*DDY,-XXX1*DDX,0.D0,C)
      RETURN
    4 CALL RSPIN( XXX1*DDX,-XXX1*DDY,0.D0,C)
      RETURN
    5 GGCAV=YYY+XXX*DL(II)             
      CALL RSPIN((1.D0+S)*GGCAV*GGY,-(1.D0+S)*GGCAV*GGX,0.D0,C)
      RETURN
    6 CALL RSPIN(0.D0,-XXX1+GG,XXX2*GGY,C)
      RETURN
    7 CALL RSPIN(+XXX1-GG,0.D0,-XXX2*GGX,C)
      RETURN
    8 CALL RSPIN(-XXX1*DDX*DDY,-XXX1*(DDX**2-DDY**2)/2.,0.D0,C)
      RETURN
    9 CALL RSPIN(+XXX2-GG,0.D0,-XXX2*GGX,C)
      RETURN
   10 CONTINUE
C      Sept 2003: fix this so that the C.O. effects can be killed.
      FDX = DX(II)                                    *1.D0
      FDY = DY(II)                                    *1.D0
      CALL RSPIN( XXX1*FDX*0.5D0/YYY,XXX1*FDY*0.5D0/YYY,0.D0,C)
      CALL SOLXYP(II,XXX/YYY,DXPA,DYPA)
      DXPA = DXPA                                     *1.D0 ! Switch off C.O. Effect. 
      DYPA = DYPA                                     *1.D0
      CALL RSPIN(XXX2*DXPA,XXX2*DYPA,-XXX+GG,B)
      CALL JAM333(A,B,C)
      FDX = DX(II+1)                                  *1.D0 
      FDY = DY(II+1)                                  *1.D0 
      CALL RSPIN(-XXX1*FDX*0.5D0/YYY,-XXX1*FDY*0.5D0/YYY,0.D0,B)
      CALL JAM333(C,B,A)
      RETURN
   11 CONTINUE
C=====Full snakes: Types 1,2,3 rotate around the 3 machine axes.
C     1:x, 2:s,3:z.
C     To get combinations use more types
C     Type 4: +22.5 deg wrt s.
C     Type 5: -22.5 deg wrt s.

      IF(NSOL(ITY).EQ.1)CALL RSPIN(PI  ,0.D0,0.D0,C)
      IF(NSOL(ITY).EQ.2)CALL RSPIN(0.D0,0.D0,PI  ,C)
      IF(NSOL(ITY).EQ.3)CALL RSPIN(0.D0,PI  ,0.D0,C)
      IF(NSOL(ITY).EQ.4)ANG=PI/8.D0
      IF(NSOL(ITY).EQ.4)CALL RSPIN(PI*DSIN(ANG),0.D0,PI*DCOS(ANG),C)
      IF(NSOL(ITY).EQ.5)ANG=-PI/8.D0
      IF(NSOL(ITY).EQ.5)CALL RSPIN(PI*DSIN(ANG),0.D0,PI*DCOS(ANG),C)
C
   15 CONTINUE
      CALL RSPIN(-XXX2*DDY,-XXXS+GG -XXX2*DDX,XXXS*GGY,C)
    2 CALL RSPIN(0.D0,-XXX2+GG,XXX2*GGY,C)
      CALL RSPIN(-XXX1*DDY,-XXX1*DDX,0.D0,C)
C
C     WRITE(53,1000)NAME(ITY),NSOL(ITY),C(1,1),
C    +                                 C(1,2),
C    +                                 C(1,3),
C    +                                 C(2,1),
C    +                                 C(2,2),
C    +                                 C(2,3),
C    +                                 C(3,1),
C    +                                 C(3,2),
C    +                                 C(3,3)
1000  FORMAT(' ',A8,I7,9F10.4)
C
C
C
C
      RETURN
      END
