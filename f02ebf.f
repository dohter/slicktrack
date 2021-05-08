      SUBROUTINE F02EBF(JOB,N,A,LDA,WR,WI,VR,LDVR,VI,LDVI,WORK,LWORK,
     *                  IFAIL)
C     MARK 16 RELEASE. NAG COPYRIGHT 1992.
C     MARK 20 REVISED. IER-3133 (JAN 2001).
C
C     F02EBF computes all the eigenvalues, and optionally all the
C     eigenvectors, of a real general matrix A.
C
C     F02EBF is a driver routine which calls computational routines
C     from LAPACK in Chapter F08.
C
C     .. Parameters ..
      INTEGER          IFAIL, LDA, LDVI, LDVR, LWORK, N
      CHARACTER        JOB
C     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), VI(LDVI,*), VR(LDVR,*), WI(N),
     *                 WORK(LWORK), WR(N)
      DOUBLE PRECISION VV(LDVR)
      integer i, j , ii, ic
      DOUBLE PRECISION ZERO
      PARAMETER        ( ZERO = 0.0d0 )
      CHARACTER JOBVR, INTEGER ICOUNT
      IF (JOB.EQ.'N') THEN
        JOBVR = 'N'
      ELSE
        JOBVR = 'V'
      END IF
      CALL dgeev('N',JOBVR,N,A,LDA,WR,WI,VI,LDVI,VR,LDVR,
     +  WORK,LWORK,IFAIL)
      DO I = 1, N
         J = 1
         DO WHILE( J.LE.N )
            IF( WI( J ).EQ.ZERO) THEN
               VI(I,J) = 0.d0
               J = J + 1
            ELSE
               VI(I,J) =    VR(I,J+1)
               VI(I,J+1) = -VR(I,J+1)
               VR(I,J+1) = VR(I,J)
               J = J + 2
            END IF
         END DO
      END DO
      END
