C   10/01/80 311092019  MEMBER NAME  EIV6     (MAY92.S)     FORTRAN
      SUBROUTINE EIV6(A,EVECT,LAMR,LAMI,IERRO)
      IMPLICIT REAL*8(A-B,D-H,O-Z),COMPLEX*16(C)
      REAL*8 C
      REAL*8  LAMR,LAMI,NN
      DIMENSION CLAM(6),CM(6),CN(6),CO1(6),CO2(6),CO3(6),
     *          CX(6,6),CNX(6,6),RCNX(2,6,6),EVECT(6,6)
      DIMENSION A(6,6),Z(6,7),ZI(6,6),LL(6),MM(6),P(7),NN(7),LAMR(6),
     *          LAMI(6),R(6),CA(6,6,6),CD(6,6)
      EQUIVALENCE (RCNX(1,1,1),CNX(1,1))
C
C
      Z(1,1)=1.2D0
      Z(2,1)=-0.9D0
      Z(3,1)=3.4D0
      Z(4,1)=-7.2D0
      Z(5,1)=5.1D0
      Z(6,1)=-9.9D0
C
C
C
C     WRITE(6,100) ((A(I,K),K=1,6),I=1,6)
  100 FORMAT(' DIE MATRIX A(I,K):'/(' ',6G16.8))
C
C     BERECHNUNG DER VEKTOREN Z(I,K)
      DO 70 I=1,6
   70 CALL DGMPRD (A,Z(1,I),Z(1,I+1),6,6,1)
C
C     BERECHNUNG DER INVERSEN MATRIX ZI
      CALL DMCPY(Z,ZI,6,6,0)
      CALL DMINV(ZI,6,DET,LL,MM)
C     WRITE(6,101)   ((ZI(I,K),K=1,6),I=1,6)
  101 FORMAT('0DIE INVERSE MATRIX ZI(M,N):'/(' ',6G12.5))
C
C     BERECHNUNG DES KOEFFIZIENTENVEKTORS P:
      CALL DGMPRD(ZI,Z(1,7),P,6,6,1)
      P(7)=-1.D0
C     WRITE(6,102) DET
  102 FORMAT(' DET=',D16.8)
C
C     BERECHNUNG DER EIGENWERTE
      CALL DPOLRT(P,NN,6,LAMR,LAMI,IERRO)
      IF(IERRO .EQ. 0) GOTO 30
      WRITE(6,109) IERRO
  109 FORMAT(' ERROR IN POLYNOMIAL NULLPOINTS--EIV6')
      RETURN
C
   30 CONTINUE
C  30 WRITE(6,103)  LAMR,LAMI
  103 FORMAT('0DIE EIGENWERTE LAM=LAMR+I*LAMI:'/(' ',6G16.8))
C
C     BERECHNUNG DER KOEFFIZIENTEN FUER DIE EIGENVEKTOREN
      DO 20 I=1,6
   20 CLAM(I)=DCMPLX(LAMR(I),LAMI(I))
      CMA=CLAM(1)*CLAM(2)*CLAM(3)*CLAM(4)*CLAM(5)*CLAM(6)
      DO 11 I=1,6
   11 CM(I)=CMA/CLAM(I)
C
      CNA=CLAM(1)+CLAM(2)+CLAM(3)+CLAM(4)+CLAM(5)+CLAM(6)
      DO 12 I=1,6
   12 CN(I)=CNA-CLAM(I)
C
      CO1(1)=CLAM(3)*CLAM(4)*CLAM(5)*CLAM(6)
     2      +CLAM(4)*CLAM(5)*CLAM(6)*CLAM(2)
     3      +CLAM(5)*CLAM(6)*CLAM(2)*CLAM(3)
     4      +CLAM(6)*CLAM(2)*CLAM(3)*CLAM(4)
     5      +CLAM(2)*CLAM(3)*CLAM(4)*CLAM(5)
C
      CO1(2)=CLAM(3)*CLAM(4)*CLAM(5)*CLAM(6)
     2      +CLAM(4)*CLAM(5)*CLAM(6)*CLAM(1)
     3      +CLAM(5)*CLAM(6)*CLAM(1)*CLAM(3)
     4      +CLAM(6)*CLAM(1)*CLAM(3)*CLAM(4)
     5      +CLAM(1)*CLAM(3)*CLAM(4)*CLAM(5)
C
      CO1(3)=CLAM(1)*CLAM(4)*CLAM(5)*CLAM(6)
     2      +CLAM(4)*CLAM(5)*CLAM(6)*CLAM(2)
     3      +CLAM(5)*CLAM(6)*CLAM(2)*CLAM(1)
     4      +CLAM(6)*CLAM(2)*CLAM(1)*CLAM(4)
     5      +CLAM(2)*CLAM(1)*CLAM(4)*CLAM(5)
C
      CO1(4)=CLAM(3)*CLAM(1)*CLAM(5)*CLAM(6)
     2      +CLAM(1)*CLAM(5)*CLAM(6)*CLAM(2)
     3      +CLAM(5)*CLAM(6)*CLAM(2)*CLAM(3)
     4      +CLAM(6)*CLAM(2)*CLAM(3)*CLAM(1)
     5      +CLAM(2)*CLAM(3)*CLAM(1)*CLAM(5)
C
      CO1(5)=CLAM(3)*CLAM(4)*CLAM(1)*CLAM(6)
     2      +CLAM(4)*CLAM(1)*CLAM(6)*CLAM(2)
     3      +CLAM(1)*CLAM(6)*CLAM(2)*CLAM(3)
     4      +CLAM(6)*CLAM(2)*CLAM(3)*CLAM(4)
     5      +CLAM(2)*CLAM(3)*CLAM(4)*CLAM(1)
C
      CO1(6)=CLAM(3)*CLAM(4)*CLAM(5)*CLAM(1)
     2      +CLAM(4)*CLAM(5)*CLAM(1)*CLAM(2)
     3      +CLAM(5)*CLAM(1)*CLAM(2)*CLAM(3)
     4      +CLAM(1)*CLAM(2)*CLAM(3)*CLAM(4)
     5      +CLAM(2)*CLAM(3)*CLAM(4)*CLAM(5)
C
      CO2(1)=CLAM(4)*CLAM(5)*CLAM(6)
     2      +CLAM(3)*CLAM(5)*CLAM(6)
     3      +CLAM(3)*CLAM(4)*CLAM(6)
     4      +CLAM(3)*CLAM(4)*CLAM(5)
     5      +CLAM(2)*CLAM(5)*CLAM(6)
     6      +CLAM(2)*CLAM(4)*CLAM(6)
     7      +CLAM(2)*CLAM(4)*CLAM(5)
     8      +CLAM(2)*CLAM(3)*CLAM(6)
     9      +CLAM(2)*CLAM(3)*CLAM(5)
     A      +CLAM(2)*CLAM(3)*CLAM(4)
C
      CO2(2)=CLAM(4)*CLAM(5)*CLAM(6)
     2      +CLAM(3)*CLAM(5)*CLAM(6)
     3      +CLAM(3)*CLAM(4)*CLAM(6)
     4      +CLAM(3)*CLAM(4)*CLAM(5)
     5      +CLAM(1)*CLAM(5)*CLAM(6)
     6      +CLAM(1)*CLAM(4)*CLAM(6)
     7      +CLAM(1)*CLAM(4)*CLAM(5)
     8      +CLAM(1)*CLAM(3)*CLAM(6)
     9      +CLAM(1)*CLAM(3)*CLAM(5)
     A      +CLAM(1)*CLAM(3)*CLAM(4)
C
      CO2(3)=CLAM(4)*CLAM(5)*CLAM(6)
     2      +CLAM(1)*CLAM(5)*CLAM(6)
     3      +CLAM(1)*CLAM(4)*CLAM(6)
     4      +CLAM(1)*CLAM(4)*CLAM(5)
     5      +CLAM(2)*CLAM(5)*CLAM(6)
     6      +CLAM(2)*CLAM(4)*CLAM(6)
     7      +CLAM(2)*CLAM(4)*CLAM(5)
     8      +CLAM(2)*CLAM(1)*CLAM(6)
     9      +CLAM(2)*CLAM(1)*CLAM(5)
     A      +CLAM(2)*CLAM(1)*CLAM(4)
C
      CO2(4)=CLAM(1)*CLAM(5)*CLAM(6)
     2      +CLAM(3)*CLAM(5)*CLAM(6)
     3      +CLAM(3)*CLAM(1)*CLAM(6)
     4      +CLAM(3)*CLAM(1)*CLAM(5)
     5      +CLAM(2)*CLAM(5)*CLAM(6)
     6      +CLAM(2)*CLAM(1)*CLAM(6)
     7      +CLAM(2)*CLAM(1)*CLAM(5)
     8      +CLAM(2)*CLAM(3)*CLAM(6)
     9      +CLAM(2)*CLAM(3)*CLAM(5)
     A      +CLAM(2)*CLAM(3)*CLAM(1)
C
      CO2(5)=CLAM(4)*CLAM(1)*CLAM(6)
     2      +CLAM(3)*CLAM(1)*CLAM(6)
     3      +CLAM(3)*CLAM(4)*CLAM(6)
     4      +CLAM(3)*CLAM(4)*CLAM(1)
     5      +CLAM(2)*CLAM(1)*CLAM(6)
     6      +CLAM(2)*CLAM(4)*CLAM(6)
     7      +CLAM(2)*CLAM(4)*CLAM(1)
     8      +CLAM(2)*CLAM(3)*CLAM(6)
     9      +CLAM(2)*CLAM(3)*CLAM(1)
     A      +CLAM(2)*CLAM(3)*CLAM(4)
C
      CO2(6)=CLAM(4)*CLAM(5)*CLAM(1)
     2      +CLAM(3)*CLAM(5)*CLAM(1)
     3      +CLAM(3)*CLAM(4)*CLAM(1)
     4      +CLAM(3)*CLAM(4)*CLAM(5)
     5      +CLAM(2)*CLAM(5)*CLAM(1)
     6      +CLAM(2)*CLAM(4)*CLAM(1)
     7      +CLAM(2)*CLAM(4)*CLAM(5)
     8      +CLAM(2)*CLAM(3)*CLAM(1)
     9      +CLAM(2)*CLAM(3)*CLAM(5)
     A      +CLAM(2)*CLAM(3)*CLAM(4)
C
      CO3(1)=CLAM(2)*CLAM(3)+CLAM(2)*CLAM(4)+CLAM(2)*CLAM(5)
     3      +CLAM(2)*CLAM(6)+CLAM(3)*CLAM(4)+CLAM(3)*CLAM(5)
     4      +CLAM(3)*CLAM(6)+CLAM(4)*CLAM(5)+CLAM(4)*CLAM(6)
     5      +CLAM(5)*CLAM(6)
C
      CO3(2)=CLAM(1)*CLAM(3)+CLAM(1)*CLAM(4)+CLAM(1)*CLAM(5)
     2      +CLAM(1)*CLAM(6)+CLAM(3)*CLAM(4)+CLAM(3)*CLAM(5)
     3      +CLAM(3)*CLAM(6)+CLAM(4)*CLAM(5)+CLAM(4)*CLAM(6)
     4      +CLAM(5)*CLAM(6)
C
      CO3(3)=CLAM(2)*CLAM(1)+CLAM(2)*CLAM(4)+CLAM(2)*CLAM(5)
     2      +CLAM(2)*CLAM(6)+CLAM(1)*CLAM(4)+CLAM(1)*CLAM(5)
     3      +CLAM(1)*CLAM(6)+CLAM(4)*CLAM(5)+CLAM(4)*CLAM(6)
     4      +CLAM(5)*CLAM(6)
C
      CO3(4)=CLAM(2)*CLAM(3)+CLAM(2)*CLAM(1)+CLAM(2)*CLAM(5)
     2      +CLAM(2)*CLAM(6)+CLAM(3)*CLAM(1)+CLAM(3)*CLAM(5)
     3      +CLAM(3)*CLAM(6)+CLAM(1)*CLAM(5)+CLAM(1)*CLAM(6)
     4      +CLAM(5)*CLAM(6)
C
      CO3(5)=CLAM(2)*CLAM(3)+CLAM(2)*CLAM(4)+CLAM(2)*CLAM(1)
     2      +CLAM(2)*CLAM(6)+CLAM(3)*CLAM(4)+CLAM(3)*CLAM(1)
     3      +CLAM(3)*CLAM(6)+CLAM(4)*CLAM(1)+CLAM(4)*CLAM(6)
     4      +CLAM(1)*CLAM(6)
C
      CO3(6)=CLAM(2)*CLAM(3)+CLAM(2)*CLAM(4)+CLAM(2)*CLAM(5)
     2      +CLAM(2)*CLAM(1)+CLAM(3)*CLAM(4)+CLAM(3)*CLAM(5)
     3      +CLAM(3)*CLAM(1)+CLAM(4)*CLAM(5)+CLAM(4)*CLAM(1)
     4      +CLAM(5)*CLAM(1)
C
C     BERECHNUNG DER EIGENVEKTOREN
      DO 13 I=1,6
      DO 15 K=1,6
   15 CX(K,I)=-CM(I)*Z(K,1)+CO1(I)*Z(K,2)-CO2(I)*Z(K,3)
     *        +CO3(I)*Z(K,4)-CN(I)*Z(K,5)+Z(K,6)
      R(I)=0.
      DO 14 K=1,6
   14 R(I)=R(I)+CX(K,I)*DCONJG(CX(K,I))
      R(I)=DSQRT(R(I))
      DO 16 K=1,6
   16 CNX(K,I)=CX(K,I)/R(I)
   13 CONTINUE
C
C
C     WRITE(6,104) CNX
  104 FORMAT('0DIE EIGENVEKTOREN CNX(K,I):',6(/6(/' ',2G12.5)))
      DO 31 I=1,6
      DO 31 J=1,6,2
      EVE         =RCNX(1,I,J)
      IF(DABS(EVE) .LT. 1.E-14) EVE=0
      EVECT(I,J  )=EVE
      EVE         =RCNX(2,I,J)
      IF(DABS(EVE) .LT. 1.E-14) EVE=0
      EVECT(I,J+1)=EVE
   31 CONTINUE
C     WRITE(6,108) ((EVECT(I,J),J=1,6),I=1,6)
  108 FORMAT('0DIE REELLEN EIGENVECTOREN:',6(/1X,6G16.8))
C
      RETURN
C
C
C     PROBE MACHEN
C
      DO 21 N=1,6
      DO 22 I=1,6
      DO 23 K=1,6
   23 CA(K,I,N)=A(K,I)
   22 CONTINUE
   21 CONTINUE
      DO 24 N=1,6
      DO 25 I=1,6
   25 CA(I,I,N) = CA(I,I,N) - CLAM(N)
   24 CONTINUE
      DO 26 N=1,6
      CALL DCGMPR(CA(1,1,N),CNX(1,N),CD(1,N),6,6,1)
   26 WRITE(6,105)(CD(K,N),K=1,6)
  105 FORMAT ('0PROBEVEKTOR CD(K,N):'(' ',2G12.5))
      RETURN
      END
