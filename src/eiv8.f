C   16/01/80 311092020  MEMBER NAME  EIV8     (MAY92.S)     FORTRAN
      SUBROUTINE EIV8(A,EVECT,LAMR,LAMI,IERRO)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      REAL*8 C
      REAL*8  LAMR,LAMI,NN
      DIMENSION CLAM(8),CM(8),CN(8),
     *          CX(8,8),CNX(8,8)
      DIMENSION A(8,8),Z(8,9),ZI(8,8),LL(8),MM(8),P(9),NN(9),LAMR(8),
     *          LAMI(8),R(8),CA(8,8,8),CD(8,8),
     *          RCNX(2,8,8),EVECT(8,8)
      DIMENSION ASP(8),AMAX(8),CV(8),CP(8,8,8),CB(8,8,8)
      EQUIVALENCE (RCNX(1,1,1),CNX(1,1))
C
C
      Z(1,1)=1.2D0
      Z(2,1)=-0.9D0
      Z(3,1)=3.4D0
      Z(4,1)=-7.2D0
      Z(5,1)=5.1D0
      Z(6,1)=-9.9D0
      Z(7,1)=2.3D0
      Z(8,1)=-6.7D0
C
C
C     WRITE(6,906) ((A(I,K),K=1,8),I=1,8)
  906 FORMAT(' DIE MATRIX A(I,K):'/(' ',8G16.8))
C
C     BERECHNUNG DER VEKTOREN Z(I,K)
      DO 70 I=1,8
   70 CALL DGMPRD (A,Z(1,I),Z(1,I+1),8,8,1)
C
C     BERECHNUNG DER INVERSEN MATRIX ZI
      CALL DMCPY(Z,ZI,8,8,0)
      CALL DMINV(ZI,8,DET,LL,MM)
C     WRITE(6,4) ((ZI(I,K),K=1,8),I=1,8)
    4 FORMAT(' DIE INVERSE MATRIX ZI(M,N):'/(' ',8G12.5))
C
C     BERECHNUNG DES KOEFFIZIENTENVEKTORS P:
      CALL DGMPRD(ZI,Z(1,9),P,8,8,1)
      P(9)=-1.D0
C     WRITE(6,5) DET
    5 FORMAT(' DET=',D16.8)
C
C     BERECHNUNG DER EIGENWERTE
      CALL DPOLRT(P,NN,8,LAMR,LAMI,IERRO)
      IF(IERRO .EQ. 0) GOTO 30
      WRITE(6,109) IERRO
  109 FORMAT(' ERROR IN POLYNOMIAL NULLPOINTS--EIV8')
      RETURN

C  30 WRITE (6,6)  LAMR,LAMI
    6 FORMAT(' DIE EIGENWERTE LAM=LAMR+LAMI:'/(' ',8G16.8))
   30 DO 20 I=1,8
   20 CLAM(I)=DCMPLX(LAMR(I),LAMI(I))
C
C     BERECHNUNG DER EIGENVEKTOREN
C
      DO 21 N=1,8
      DO 22 I=1,8
      DO 23 K=1,8
   23 CA(K,I,N)=A(K,I)
   22 CONTINUE
   21 CONTINUE
      DO 24 N=1,8
      DO 25 I=1,8
   25 CA(I,I,N) = CA(I,I,N) - CLAM(N)
   24 CONTINUE
C
C
      CALL DCGMPR(CA(1,1,2),CA(1,1,3),CP(1,1,2),8,8,8)
      CALL DCGMPR(CP(1,1,2),CA(1,1,4),CP(1,1,3),8,8,8)
      CALL DCGMPR(CP(1,1,3),CA(1,1,5),CP(1,1,4),8,8,8)
      CALL DCGMPR(CP(1,1,4),CA(1,1,6),CP(1,1,5),8,8,8)
      CALL DCGMPR(CP(1,1,5),CA(1,1,7),CP(1,1,6),8,8,8)
      CALL DCGMPR(CP(1,1,6),CA(1,1,8),CB(1,1,1),8,8,8)
C
      CALL DCGMPR(CA(1,1,3),CA(1,1,4),CP(1,1,2),8,8,8)
      CALL DCGMPR(CP(1,1,2),CA(1,1,5),CP(1,1,3),8,8,8)
      CALL DCGMPR(CP(1,1,3),CA(1,1,6),CP(1,1,4),8,8,8)
      CALL DCGMPR(CP(1,1,4),CA(1,1,7),CP(1,1,5),8,8,8)
      CALL DCGMPR(CP(1,1,5),CA(1,1,8),CP(1,1,6),8,8,8)
      CALL DCGMPR(CP(1,1,6),CA(1,1,1),CB(1,1,2),8,8,8)
C
      CALL DCGMPR(CA(1,1,4),CA(1,1,5),CP(1,1,2),8,8,8)
      CALL DCGMPR(CP(1,1,2),CA(1,1,6),CP(1,1,3),8,8,8)
      CALL DCGMPR(CP(1,1,3),CA(1,1,7),CP(1,1,4),8,8,8)
      CALL DCGMPR(CP(1,1,4),CA(1,1,8),CP(1,1,5),8,8,8)
      CALL DCGMPR(CP(1,1,5),CA(1,1,1),CP(1,1,6),8,8,8)
      CALL DCGMPR(CP(1,1,6),CA(1,1,2),CB(1,1,3),8,8,8)
C
      CALL DCGMPR(CA(1,1,5),CA(1,1,6),CP(1,1,2),8,8,8)
      CALL DCGMPR(CP(1,1,2),CA(1,1,7),CP(1,1,3),8,8,8)
      CALL DCGMPR(CP(1,1,3),CA(1,1,8),CP(1,1,4),8,8,8)
      CALL DCGMPR(CP(1,1,4),CA(1,1,1),CP(1,1,5),8,8,8)
      CALL DCGMPR(CP(1,1,5),CA(1,1,2),CP(1,1,6),8,8,8)
      CALL DCGMPR(CP(1,1,6),CA(1,1,3),CB(1,1,4),8,8,8)
C
      CALL DCGMPR(CA(1,1,6),CA(1,1,7),CP(1,1,2),8,8,8)
      CALL DCGMPR(CP(1,1,2),CA(1,1,8),CP(1,1,3),8,8,8)
      CALL DCGMPR(CP(1,1,3),CA(1,1,1),CP(1,1,4),8,8,8)
      CALL DCGMPR(CP(1,1,4),CA(1,1,2),CP(1,1,5),8,8,8)
      CALL DCGMPR(CP(1,1,5),CA(1,1,3),CP(1,1,6),8,8,8)
      CALL DCGMPR(CP(1,1,6),CA(1,1,4),CB(1,1,5),8,8,8)
C
      CALL DCGMPR(CA(1,1,7),CA(1,1,8),CP(1,1,2),8,8,8)
      CALL DCGMPR(CP(1,1,2),CA(1,1,1),CP(1,1,3),8,8,8)
      CALL DCGMPR(CP(1,1,3),CA(1,1,2),CP(1,1,4),8,8,8)
      CALL DCGMPR(CP(1,1,4),CA(1,1,3),CP(1,1,5),8,8,8)
      CALL DCGMPR(CP(1,1,5),CA(1,1,4),CP(1,1,6),8,8,8)
      CALL DCGMPR(CP(1,1,6),CA(1,1,5),CB(1,1,6),8,8,8)
C
      CALL DCGMPR(CA(1,1,8),CA(1,1,1),CP(1,1,2),8,8,8)
      CALL DCGMPR(CP(1,1,2),CA(1,1,2),CP(1,1,3),8,8,8)
      CALL DCGMPR(CP(1,1,3),CA(1,1,3),CP(1,1,4),8,8,8)
      CALL DCGMPR(CP(1,1,4),CA(1,1,4),CP(1,1,5),8,8,8)
      CALL DCGMPR(CP(1,1,5),CA(1,1,5),CP(1,1,6),8,8,8)
      CALL DCGMPR(CP(1,1,6),CA(1,1,6),CB(1,1,7),8,8,8)
C
      CALL DCGMPR(CA(1,1,1),CA(1,1,2),CP(1,1,2),8,8,8)
      CALL DCGMPR(CP(1,1,2),CA(1,1,3),CP(1,1,3),8,8,8)
      CALL DCGMPR(CP(1,1,3),CA(1,1,4),CP(1,1,4),8,8,8)
      CALL DCGMPR(CP(1,1,4),CA(1,1,5),CP(1,1,5),8,8,8)
      CALL DCGMPR(CP(1,1,5),CA(1,1,6),CP(1,1,6),8,8,8)
      CALL DCGMPR(CP(1,1,6),CA(1,1,7),CB(1,1,8),8,8,8)
C
C     WRITE(6,100) ((CB(I,K,1),K=1,4),I=1,8)
  100 FORMAT(' DIE MATRIX CB(I,K,1):'/(' ',4G12.5))
C
      DO 203 N=1,8
      DO 201 K=1,8
      ASP(K)=0.D0
      DO 200 I=1,8
  200 ASP(K)=ASP(K)+CB(I,K,N)*DCONJG(CB(I,K,N))
  201 CONTINUE
      IZW=1
      AMAX(1)=ASP(1)
      DO 202 J=2,8
      IF(ASP(J) .LE. AMAX(1)) GO TO 202
      AMAX(1)=ASP(J)
      IZW=J
  202 CONTINUE
C
      CV(1)=0.D0
      CV(2)=0.D0
      CV(3)=0.D0
      CV(4)=0.D0
      CV(5)=0.D0
      CV(6)=0.D0
      CV(7)=0.D0
      CV(8)=0.D0
      CV(IZW)=1.D0
C
      CALL DCGMPR(CB(1,1,N),CV,CX(1,N),8,8,1)
  203 CONTINUE
C
      DO 13 I=1,8
      R(I)=0.D0
      DO 14 K=1,8
   14 R(I)=R(I)+CX(K,I)*DCONJG(CX(K,I))
      R(I)=DSQRT(R(I))
      DO 16 K=1,8
   16 CNX(K,I)=CX(K,I)/R(I)
   13 CONTINUE
C
C
      DO 31 I=1,8
      DO 31 J=1,8,2
      EVE         =RCNX(1,I,J)
      IF(DABS(EVE) .LT. 1.E-14) EVE=0
      EVECT(I,J  )=EVE
      EVE         =RCNX(2,I,J)
      IF(DABS(EVE) .LT. 1.E-14) EVE=0
      EVECT(I,J+1)=EVE
   31 CONTINUE
C     WRITE(6,108) ((EVECT(I,J),J=1,8),I=1,8)
  108 FORMAT('0DIE REELLEN EIGENVECTOREN:',8(/1X,8G16.8))
C
C
C
C      WRITE(6,40) CNX(1,1),CNX(2,1),CNX(3,1),CNX(4,1),
C    *            CNX(5,1),CNX(6,1),CNX(7,1),CNX(8,1)
   40 FORMAT(' CNX(1,1)=',2G12.5/
     2       ' CNX(2,1)=',2G12.5/
     3       ' CNX(3,1)=',2G12.5/
     4       ' CNX(4,1)=',2G12.5/
     5       ' CNX(5,1)=',2G12.5/
     6       ' CNX(6,1)=',2G12.5/
     7       ' CNX(7,1)=',2G12.5/
     8       ' CNX(8,1)=',2G12.5/)
C
C     WRITE(6,41) CNX(1,2),CNX(2,2),CNX(3,2),CNX(4,2),
C    *            CNX(5,2),CNX(6,2),CNX(7,2),CNX(8,2)
   41 FORMAT(' CNX(1,2)=',2G12.5/
     2       ' CNX(2,2)=',2G12.5/
     3       ' CNX(3,2)=',2G12.5/
     4       ' CNX(4,2)=',2G12.5/
     5       ' CNX(5,2)=',2G12.5/
     6       ' CNX(6,2)=',2G12.5/
     7       ' CNX(7,2)=',2G12.5/
     8       ' CNX(8,2)=',2G12.5/)
C
C     WRITE(6,42) CNX(1,3),CNX(2,3),CNX(3,3),CNX(4,3),
C    *            CNX(5,3),CNX(6,3),CNX(7,3),CNX(8,3)
   42 FORMAT(' CNX(1,3)=',2G12.5/
     2       ' CNX(2,3)=',2G12.5/
     3       ' CNX(3,3)=',2G12.5/
     4       ' CNX(4,3)=',2G12.5/
     5       ' CNX(5,3)=',2G12.5/
     6       ' CNX(6,3)=',2G12.5/
     7       ' CNX(7,3)=',2G12.5/
     8       ' CNX(8,3)=',2G12.5/)
C
C     WRITE(6,43) CNX(1,4),CNX(2,4),CNX(3,4),CNX(4,4),
C    *            CNX(5,4),CNX(6,4),CNX(7,4),CNX(8,4)
   43 FORMAT(' CNX(1,4)=',2G12.5/
     2       ' CNX(2,4)=',2G12.5/
     3       ' CNX(3,4)=',2G12.5/
     4       ' CNX(4,4)=',2G12.5/
     5       ' CNX(5,4)=',2G12.5/
     6       ' CNX(6,4)=',2G12.5/
     7       ' CNX(7,4)=',2G12.5/
     8       ' CNX(8,4)=',2G12.5/)
C
C     WRITE(6,44) CNX(1,5),CNX(2,5),CNX(3,5),CNX(4,5),
C    *            CNX(5,5),CNX(6,5),CNX(7,5),CNX(8,5)
   44 FORMAT(' CNX(1,5)=',2G12.5/
     2       ' CNX(2,5)=',2G12.5/
     3       ' CNX(3,5)=',2G12.5/
     4       ' CNX(4,5)=',2G12.5/
     5       ' CNX(5,5)=',2G12.5/
     6       ' CNX(6,5)=',2G12.5/
     7       ' CNX(7,5)=',2G12.5/
     8       ' CNX(8,5)=',2G12.5/)
C
C     WRITE(6,45) CNX(1,6),CNX(2,6),CNX(3,6),CNX(4,6),
C    *            CNX(5,6),CNX(6,6),CNX(7,6),CNX(8,6)
   45 FORMAT(' CNX(1,6)=',2G12.5/
     2       ' CNX(2,6)=',2G12.5/
     3       ' CNX(3,6)=',2G12.5/
     4       ' CNX(4,6)=',2G12.5/
     5       ' CNX(5,6)=',2G12.5/
     6       ' CNX(6,6)=',2G12.5/
     7       ' CNX(7,6)=',2G12.5/
     8       ' CNX(8,6)=',2G12.5/)
C
C     WRITE(6,46) CNX(1,7),CNX(2,7),CNX(3,7),CNX(4,7),
C    *            CNX(5,7),CNX(6,7),CNX(7,7),CNX(8,7)
   46 FORMAT(' CNX(1,7)=',2G12.5/
     2       ' CNX(2,7)=',2G12.5/
     3       ' CNX(3,7)=',2G12.5/
     4       ' CNX(4,7)=',2G12.5/
     5       ' CNX(5,7)=',2G12.5/
     6       ' CNX(6,7)=',2G12.5/
     7       ' CNX(7,7)=',2G12.5/
     8       ' CNX(8,7)=',2G12.5/)
C
C     WRITE(6,48) CNX(1,8),CNX(2,8),CNX(3,8),CNX(4,8),
C    *            CNX(5,8),CNX(6,8),CNX(7,8),CNX(8,8)
   48 FORMAT(' CNX(1,8)=',2G12.5/
     2       ' CNX(2,8)=',2G12.5/
     3       ' CNX(3,8)=',2G12.5/
     4       ' CNX(4,8)=',2G12.5/
     5       ' CNX(5,8)=',2G12.5/
     6       ' CNX(6,8)=',2G12.5/
     7       ' CNX(7,8)=',2G12.5/
     8       ' CNX(8,8)=',2G12.5/)
C
      CALL DCGMPR(CA(1,1,1),CNX(1,1),CD(1,1),8,8,1)
      CALL DCGMPR(CA(1,1,2),CNX(1,2),CD(1,2),8,8,1)
      CALL DCGMPR(CA(1,1,3),CNX(1,3),CD(1,3),8,8,1)
      CALL DCGMPR(CA(1,1,4),CNX(1,4),CD(1,4),8,8,1)
      CALL DCGMPR(CA(1,1,5),CNX(1,5),CD(1,5),8,8,1)
      CALL DCGMPR(CA(1,1,6),CNX(1,6),CD(1,6),8,8,1)
      CALL DCGMPR(CA(1,1,7),CNX(1,7),CD(1,7),8,8,1)
      CALL DCGMPR(CA(1,1,8),CNX(1,8),CD(1,8),8,8,1)
C
C     WRITE (6,50) CD(1,1),CD(2,1),CD(3,1),CD(4,1),
C    *             CD(5,1),CD(6,1),CD(7,1),CD(8,1)
   50 FORMAT(' CD(1,1)=',2G12.5/
     2       ' CD(2,1)=',2G12.5/
     3       ' CD(3,1)=',2G12.5/
     4       ' CD(4,1)=',2G12.5/
     5       ' CD(5,1)=',2G12.5/
     6       ' CD(6,1)=',2G12.5/
     7       ' CD(7,1)=',2G12.5/
     8       ' CD(8,1)=',2G12.5/)
C
C     WRITE (6,51) CD(1,2),CD(2,2),CD(3,2),CD(4,2),
C    *             CD(5,2),CD(6,2),CD(7,2),CD(8,2)
   51 FORMAT(' CD(1,2)=',2G12.5/
     2       ' CD(2,2)=',2G12.5/
     3       ' CD(3,2)=',2G12.5/
     4       ' CD(4,2)=',2G12.5/
     5       ' CD(5,2)=',2G12.5/
     6       ' CD(6,2)=',2G12.5/
     7       ' CD(7,2)=',2G12.5/
     8       ' CD(8,2)=',2G12.5/)
C
C     WRITE (6,52) CD(1,3),CD(2,3),CD(3,3),CD(4,3),
C    *             CD(5,3),CD(6,3),CD(7,3),CD(8,3)
   52 FORMAT(' CD(1,3)=',2G12.5/
     2       ' CD(2,3)=',2G12.5/
     3       ' CD(3,3)=',2G12.5/
     4       ' CD(4,3)=',2G12.5/
     5       ' CD(5,3)=',2G12.5/
     6       ' CD(6,3)=',2G12.5/
     7       ' CD(7,3)=',2G12.5/
     8       ' CD(8,3)=',2G12.5/)
C
C     WRITE (6,53) CD(1,4),CD(2,4),CD(3,4),CD(4,4),
C    *             CD(5,4),CD(6,4),CD(7,4),CD(8,4)
   53 FORMAT(' CD(1,4)=',2G12.5/
     2       ' CD(2,4)=',2G12.5/
     3       ' CD(3,4)=',2G12.5/
     4       ' CD(4,4)=',2G12.5/
     5       ' CD(5,4)=',2G12.5/
     6       ' CD(6,4)=',2G12.5/
     7       ' CD(7,4)=',2G12.5/
     8       ' CD(8,4)=',2G12.5/)
C
C     WRITE (6,54) CD(1,5),CD(2,5),CD(3,5),CD(4,5),
C    *             CD(5,5),CD(6,5),CD(7,5),CD(8,5)
   54 FORMAT(' CD(1,5)=',2G12.5/
     2       ' CD(2,5)=',2G12.5/
     3       ' CD(3,5)=',2G12.5/
     4       ' CD(4,5)=',2G12.5/
     5       ' CD(5,5)=',2G12.5/
     6       ' CD(6,5)=',2G12.5/
     7       ' CD(7,5)=',2G12.5/
     8       ' CD(8,5)=',2G12.5/)
C
C     WRITE (6,55) CD(1,6),CD(2,6),CD(3,6),CD(4,6),
C    *             CD(5,6),CD(6,6),CD(7,6),CD(8,6)
   55 FORMAT(' CD(1,6)=',2G12.5/
     2       ' CD(2,6)=',2G12.5/
     3       ' CD(3,6)=',2G12.5/
     4       ' CD(4,6)=',2G12.5/
     5       ' CD(5,6)=',2G12.5/
     6       ' CD(6,6)=',2G12.5/
     7       ' CD(7,6)=',2G12.5/
     8       ' CD(8,6)=',2G12.5/)
C
C     WRITE (6,56) CD(1,7),CD(2,7),CD(3,7),CD(4,7),
C    *             CD(5,7),CD(6,7),CD(7,7),CD(8,7)
   56 FORMAT(' CD(1,7)=',2G12.5/
     2       ' CD(2,7)=',2G12.5/
     3       ' CD(3,7)=',2G12.5/
     4       ' CD(4,7)=',2G12.5/
     5       ' CD(5,7)=',2G12.5/
     6       ' CD(6,7)=',2G12.5/
     7       ' CD(7,7)=',2G12.5/
     8       ' CD(8,7)=',2G12.5/)
C
C     WRITE (6,57) CD(1,8),CD(2,8),CD(3,8),CD(4,8),
C    *             CD(5,8),CD(6,8),CD(7,8),CD(8,8)
   57 FORMAT(' CD(1,8)=',2G12.5/
     2       ' CD(2,8)=',2G12.5/
     3       ' CD(3,8)=',2G12.5/
     4       ' CD(4,8)=',2G12.5/
     5       ' CD(5,8)=',2G12.5/
     6       ' CD(6,8)=',2G12.5/
     7       ' CD(7,8)=',2G12.5/
     8       ' CD(8,8)=',2G12.5/)
C
      RETURN
      END
