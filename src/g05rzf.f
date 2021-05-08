
      SUBROUTINE G05RZF (MODE, N, M, XMU, C, LDC, R, LR, STATE, X, LDX
     &  , IFAIL)
      INTEGER  MODE, N, M, LDC, LR, STATE(*), LDX, IFAIL
      DOUBLE PRECISION  XMU(M), C(LDC,M), R(LR), X(LDX,M), CTEMP(LDC,M)
      DOUBLE PRECISION WORK(M), PARM(M*(M+3)/2+1)
      DOUBLE PRECISION PART(M)
      INTEGER I

      CTEMP = C
      call setgmn ( XMU, CTEMP, M, PARM)
      DO I = 1,N
        call genmn(PARM, PART, WORK)
        X(I,:) = PART
      ENDDO
      IFAIL = 0

      END SUBROUTINE G05rzf
