      S U B R O U T I N E JAM888(AUS,EIN1,EIN2)
C     aus=ein1*ein2

      IMPLICIT REAL*8 (A,C-H,O-Z)
      REAL*8 AUS(8,8),EIN1(8,8),EIN2(8,8)
      REAL*8 S(8,8)

C     matrizenmultiplikation

      DO  1   I=1,8
         DO  2   K=1,8
            S(I,K)=0.D+0
            DO  3   J=1,8
               S(I,K)=S(I,K)+EIN1(I,J)*EIN2(J,K)
 3          CONTINUE
 2       CONTINUE
 1    CONTINUE
      DO  5   I=1,8
         DO  4   K=1,8
            AUS(I,K)=S(I,K)
 4       CONTINUE
 5    CONTINUE

      RETURN

      E N D 
