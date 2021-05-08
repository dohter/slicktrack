C   15/11/79 501131548  MEMBER NAME  DAMINT   (S)           FORTRAN
      SUBROUTINE DAMINT(I,IID,ITY,CRAD,X,X2,Y,V,DAMMR,NM)
C
C==========ROUTINE TO CALCULATE DAMPING CONSTANTS USING ================
C                        M & R INTEGRALS.
C
C=====THESE SHOULD BE COMPARED WITH THE NAIVE THIN LENS RESULTS.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(6,6),V(6,6),DAMMR(6)
      CHARACTER *8 NM
C
C
C
C
C
C      CALL VZERO(A,72)
      A = 0.D0
      CALL MXDAMP(I,IID,ITY,A,CRAD,X,X2,Y,NM)
C
C
C
C      IF(IID.EQ.5)WRITE(53,'(A,2F20.10)')'CAVITY ', A(2,2),A(4,4)


      DO 10 K=1,5,2
      DAMMR(K)=DAMMR(K)
     +        -A(2,2)*(V(1,K+1)*V(2,K)-V(1,K)*V(2,K+1))
     +        -A(4,4)*(V(3,K+1)*V(4,K)-V(3,K)*V(4,K+1))
     +        -A(6,1)*(V(5,K+1)*V(1,K)-V(5,K)*V(1,K+1))
     +        -A(6,2)*(V(5,K+1)*V(2,K)-V(5,K)*V(2,K+1))
     +        -A(6,3)*(V(5,K+1)*V(3,K)-V(5,K)*V(3,K+1))
     +        -A(6,4)*(V(5,K+1)*V(4,K)-V(5,K)*V(4,K+1))
     +        -A(6,6)*(V(5,K+1)*V(6,K)-V(5,K)*V(6,K+1))
      DAMMR(K+1)=DAMMR(K)
C
C      IF(IID.EQ.5)WRITE(53,'(A,4F20.10)')'CAVITY ', 
C     +         V(1,K+1),V(2,K),V(1,K),V(2,K+1)
C
   10 CONTINUE
C
C
      RETURN
      END
