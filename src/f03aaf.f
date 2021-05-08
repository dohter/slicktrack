      SUBROUTINE F03AAF(A,LDA,N,DET,WKSPCE,IFAIL)
C     MARK 15 RE-ISSUE. NAG COPYRIGHT 1991.
C
C     Determinant of real matrix.
C     1st August 1971
C
C     Rewritten to call F07ADGN, a modified version of LAPACK routine
C     SGETRF/F07ADFN; new IFAIL exit inserted for illegal input
C     parameters; error messages inserted. February 1991.
C
C     .. Parameters ..
      CHARACTER*6      SRNAME
      DOUBLE PRECISION ONE, ZERO
      PARAMETER        (ONE=1.0D+0,ZERO=0.0D+0)
C     .. Scalar Arguments ..
      DOUBLE PRECISION DET
      INTEGER          IFAIL, LDA, N
      INTEGER          IPIV(N)
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), WKSPCE(N)
C     Old: IFAIL = 0
C     CALL f03aaf(A,LDA,N,DET,WKSPCE,IFAIL)
C     New: INTEGER IPIV(N)
C     ...

      CALL dgetrf(N,N,A,LDA,IPIV,INFO)
      IFAIL = 0
      DET = 1.d0
      DO I = 1,N
        if (ipiv(i).ne.i) then
          DET = -DET * A(I,I)
        else
          DET =  DET * A(I,I)
        endif
      ENDDO
      END
