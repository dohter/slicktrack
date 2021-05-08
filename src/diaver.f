C   15/11/79 205261903  MEMBER NAME  DIAVER   (MAY92.S)  M  FORTRAN
      SUBROUTINE DIAVER(A,B,C,II)

C     SUBROUTINE DIAVER(A,B,C,II,NAME)
C
C=====AVERAGEING ROUTINE FOR USE WITH 'DISPERSION' VERSION OF SLICK:
C=====CALCULATES THE SCALAR PRODUCT OF THE EIGENVECTORS AND DISPERSION
C=====VECTOR AND THEN AVERAGES THEM.
C     'A' IS INPUT,'B' IS OUTPUT,'C' IS AVERAGE.
C
C=====26 May 92, change sign of C(5,J) to change sign of numerator
C=====gdn/dg term to get agreement with Chao formalism. No
C=====guarentee that this is the correct proceedure.
C
      IMPLICIT REAL*8(A-H,O-Z)
CDB      REAL*8 NAME
      DIMENSION A(8,8),B(8,8),C(8,8)
C
      INCLUDE "cdisp.for"
C
C
      D1I=D1(II)
      D2I=D2(II)
      D3I=D3(II)
      D4I=D4(II)
      D1O=D1(II+1)
      D2O=D2(II+1)
      D3O=D3(II+1)
      D4O=D4(II+1)
C     WRITE(53,41)NAME,D1I,D2I,D3I,D4I
C  41 FORMAT(' ',A8,4F15.5)
C     WRITE(53,41)D1I,D2I,D3I,D4I
   41 FORMAT(' ',4F15.5)
C
C
C
      DO 9 J=1,8
      VA=A(1,J)*D2I-A(2,J)*D1I+A(3,J)*D4I-A(4,J)*D3I-A(5,J)
      VB=B(1,J)*D2O-B(2,J)*D1O+B(3,J)*D4O-B(4,J)*D3O-B(5,J)
C      VA= -A(5,J)
C      VB= -B(5,J)
      C(5,J)=-(VA+VB)*0.5D0
      C(7,J)=(A(7,J)+B(7,J))*0.5D0
      C(8,J)=(A(8,J)+B(8,J))*0.5D0
    9 CONTINUE
C
C
C
      RETURN
      END
