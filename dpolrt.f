C   17/11/75 102091120  MEMBER NAME  DPOLRT   (JJS)         FORTRAN
      SUBROUTINEDPOLRT(XCOF,COF,M,ROOTR,ROOTI,IER)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XCOF(*),COF(*),ROOTR(*),ROOTI(*)
      IFIT=0
      N=M
      IER=0
      IF(XCOF(N+1))10,25,10
   10 IF(N) 15,15,32
   15 IER=1
   20 RETURN
   25 IER=4
      GO TO 20
   30 IER=2
      GO TO 20
   32 IF(N-36) 35,35,30
   35 NX=N
      NXX=N+1
      N2=1
      KJ1 = N+1
      DO 40 L=1,KJ1
      MT=KJ1-L+1
   40 COF(MT)=XCOF(L)
   45 XO=.00500101
      YO=0.01000101
      IN=0
   50 X=XO
      XO = - 1D1 * YO
      YO = - 1D1 * X
      X=XO
      Y=YO
      IN=IN+1
      GO TO 59
   55 IFIT=1
      XPR=X
      YPR=Y
   59 ICT=0
   60 UX = 0D0
      UY = 0D0
      V = 0D0
      YT = 0D0
      XT = 1D0
      U=COF(N+1)
      IF(U) 65,130,65
   65 DO 70 I=1,N
      L =N-I+1
      TEMP=COF(L)
      XT2=X*XT-Y*YT
      YT2=X*YT+Y*XT
      U=U+TEMP*XT2
      V=V+TEMP*YT2
      FI=I
      UX=UX+FI*XT*TEMP
      UY=UY-FI*YT*TEMP
      XT=XT2
   70 YT=YT2
      SUMSQ=UX*UX+UY*UY
      IF(SUMSQ) 75,110,75
   75 DX=(V*UY-U*UX)/SUMSQ
      X=X+DX
      DY=-(U*UY+V*UX)/SUMSQ
      Y=Y+DY
   78 IF(DABS(DY)+DABS(DX)-1D-12)100,80,80
   80 ICT=ICT+1
      IF(ICT-500) 60,85,85
   85 IF(IFIT)100,90,100
   90 IF(IN-5) 50,95,95
   95 IER=3
      GO TO 20
  100 DO 105 L=1,NXX
      MT=KJ1-L+1
      TEMP=XCOF(MT)
      XCOF(MT)=COF(L)
  105 COF(L)=TEMP
      ITEMP=N
      N=NX
      NX=ITEMP
      IF(IFIT) 120,55,120
  110 IF(IFIT) 115,50,115
  115 X=XPR
      Y=YPR
  120 IFIT=0
  122 IF(DABS(Y/X)-1D-10)135,125,125
  125 ALPHA=X+X
      SUMSQ=X*X+Y*Y
      N=N-2
      GO TO 140
  130 X=0D0
      NX=NX-1
      NXX=NXX-1
  135 Y=0D0
      SUMSQ =0D0
      ALPHA=X
      N=N-1
  140 COF(2)=COF(2)+ALPHA*COF(1)
  145 DO 150 L=2,N
  150 COF(L+1)=COF(L+1)+ALPHA*COF(L)-SUMSQ*COF(L-1)
  155 ROOTI(N2)=Y
      ROOTR(N2)=X
      N2=N2+1
      IF(SUMSQ) 160,165,160
  160 Y=-Y
      SUMSQ = 0D0
      GO TO 155
  165 IF(N) 20,20,45
      END
C   17/11/75            MEMBER NAME  DGMPRD   (S)           FORTRAN
      SUBROUTINE DGMPRD(A,B,R,N,M,L)
      REAL*8 A(1),B(1),R(1)
      IR=0
      IK=-M
      DO10 K=1,L
      IK=IK+M
      DO10J=1,N
      IR=IR+1
      JI=J-N
      IB=IK
      R(IR)=0.
      DO10 I=1,M
      JI=JI+N
      IB=IB+1
   10 R(IR)=R(IR)+A(JI)*B(IB)
      RETURN
      END
C   17/11/75            MEMBER NAME  DMCPY    (S)           FORTRAN
      SUBROUTINE DMCPY(A,R,N,M,MS)
      REAL*8 A(1),R(1)
      CALL LOC(N,M,IT,N,M,MS)
      DO1 I=1,IT
    1 R(I)=A(I)
      RETURN
      END
