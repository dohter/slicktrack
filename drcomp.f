      SUBROUTINE DRCOMP
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      COMMON/DRIFTS/DR1(10000),DR2(10000),IDR1,IDR2 
     DIMENSION DRDIFF(10000)

C
      DO 10,I=1,7000
      DRDIFF(I)=DR1(I)-DR2(I)
      WRITE(53,100)I,DR1(I),DR2(I),DRDIFF(I)
      WRITE(6,100)I,DR1(I),DR2(I),DRDIFF(I)
  100 FORMAT(' ','I,DR1,DR2,DIFF :',I5,3E20.10) 
   10 CONTINUE  
      RETURN 
      END
