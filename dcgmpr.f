C   07/01/80            MEMBER NAME  DCGMPR   (EVS)         FORTRAN
      SUBROUTINE DCGMPR(A,B,R,N,M,L)
C     COMPLEX*16  VERSION VON GMPRD
      COMPLEX*16  A(*),B(*),R(*)
      IR=0
      IK=-M
      DO 10 K=1,L
      IK=IK+M
      DO 10 J=1,N
      IR=IR+1
      JI=J-N
      IB=IK
      R(IR)=0
      DO 10 I=1,M
      JI=JI+N
      IB=IB+1
   10 R(IR)=R(IR)+A(JI)*B(IB)
      RETURN
      END
