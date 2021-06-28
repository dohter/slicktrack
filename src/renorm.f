C   22/01/90 510032056  MEMBER NAME  RENORM   (SEPT95.S) M  FORTRAN
      SUBROUTINE RENORM(A,N,WR,WI)
C  *****************************************************************
C  *        From SMILE and modified to Fortran IV.                 *
C  *        renormalizes a 6x6 or 8x8 matrix                       *
C  *        Uses NORM2 (SMILE version of NORM).                    *
C  *****************************************************************

C  (1) normalize A(I,J) to AB(J)=1 (see subr.norm).
C  (2) identify the mode (horizontal, vertical and longitudinal) and
C      reorder A(I,J) and WR,WI into this order. the spin is assumed
C      to be the 7-th and 8-th component. if the coupling is very
C      strong, the identification may fail but it does not affect the
C      final polarization.
C  (3) normalize again so that A(1,1)+I*A(1,2), A(3,3)+I*A(3,4) and
C      A(5,5)+I*A(5,6) are real positive. then, A(1,1)=SQRT(BETA) etc.
C  if the process (2) does not fail, then
C   ATAN2(WI(I),WR(I))/TWOPI (I=1,2,3) (times -1 for I=3) gives the tune
C   with correct signs.
C  only A(J,I) (J=1..N, I=1..6) are changed.

      IMPLICIT REAL*8(A-H,O-Z)

CDB      COMMON/JAMJAM/DJAM(20,20)

      DIMENSION A(N,N),WR(N),WI(N),B(8,8),C(8),C2(8),KK(3)

      CALL NORM2(A,N,C)
      IF(C(1).EQ.0.D0.OR.C(3).EQ.0.D0.OR.C(5).EQ.0.D0)THEN
      WRITE(53,'(A,A)')' Normalisation of 8-eigenvectors in RENORM',
     +                  ' crazy --  STOP.'
      STOP
      ENDIF


      DO 10 K=1,3
   10 KK(K)=0

      DO 60 II=1,3
        I=2*II-1
C-- Normalize the orbit vectors so that norm=1

C        IF (DABS(C(I)).LT.1.D-6)  C(I) = 1.D0   !Sept 2003: Not too sure what the problem is here
        IF (DABS(C(I)).LT.1.D-10)  C(I) = 1.D0  ! but the vectors sometimes get reversed.
C       C(I)=1.D0
        IF(C(I).GE.0.D0) GOTO 30
        C(I)=-C(I)

        DO 20 J=1,N
   20   A(J,I+1)=-A(J,I+1)

        WI(I)=-WI(I)
        WI(I+1)=-WI(I+1)

   30   C1=1.D0/DSQRT(C(I))

        DO 40 J=1,N
          A(J,I)=C1*A(J,I)
   40   A(J,I+1)=C1*A(J,I+1)

C-- identify x,y,z
        D = -1.D0
       DO 50 JJ=1,3
       D1=A(2*JJ-1,I)**2+A(2*JJ,I)**2+A(2*JJ-1,I+1)**2+A(2*JJ,I+1)**2
          IF(D1.LE.D) GOTO 50
          D=D1
          K=JJ
   50   CONTINUE

        KK(K)=II
   60 CONTINUE



C     IF(KK(1)*KK(2)*KK(3).NE.0) THEN
      IF(KK(1)*KK(2)*KK(3).EQ.0) GO TO 120
C-- reorder in x,y,z
        DO 80 II=1,3
          I=2*II-1
          K=2*KK(II)-1
          C(I)=WR(K)
          C(I+1)=WR(K+1)
          C2(I)=WI(K)
          C2(I+1)=WI(K+1)
          DO 70 J=1,N
            B(J,I)=A(J,K)
            B(J,I+1)=A(J,K+1)
   70     CONTINUE
   80   CONTINUE

        DO 90 I=1,6
          WR(I)=C(I)
   90   WI(I)=C2(I)

C-- phase rotation. diagonal=real.
        DO 110 II=1,3
          I=2*II-1
          R=DSQRT(B(I,I)**2+B(I,I+1)**2)
          CO=B(I,I)/R
          SI=-B(I,I+1)/R
          DO 100 J=1,N
            A(J,I)=CO*B(J,I)-SI*B(J,I+1)
  100     A(J,I+1)=SI*B(J,I)+CO*B(J,I+1)

  110   CONTINUE
C     ELSE
      RETURN

  120 CONTINUE
        WRITE(53,10000) (KK(I),I=1,3)
C     ENDIF
10000 FORMAT(' (SUBR.RENORM) failed to identify eigenvectors as X, Y ',
     +'and Z.  KK=',3I2)
      STOP
      END
