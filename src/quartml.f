       SUBROUTINE QUARTM(qc,qb,qa)

C        Routine for multiplying quarternions: e.g. page 15 of MV's thesis.
C
       
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION qa(0:3),qb(0:3),qc(0:3),qd(0:3)

      qd(0) = qb(0)*qa(0) - qb(1)*qa(1) - qb(2)*qa(2) - qb(3)*qa(3)
      
      qd(1) = qb(0)*qa(1) + qa(0)*qb(1) + qb(2)*qa(3) - qb(3)*qa(2)

      qd(2) = qb(0)*qa(2) + qa(0)*qb(2) - qb(1)*qa(3) + qb(3)*qa(1)

      qd(3) = qb(0)*qa(3) + qa(0)*qb(3) + qb(1)*qa(2) - qb(2)*qa(1)

      qc = qd

      RETURN
      END


C     qd(0) = qb(0)*qa(0) 
C     
C     qd(1) = qb(0)*qa(1) + qa(0)*qb(1)
C
C     qd(2) = qb(0)*qa(2) + qa(0)*qb(2)
C
C     qd(3) = qb(0)*qa(3) + qa(0)*qb(3)
C
C     qc = qd
