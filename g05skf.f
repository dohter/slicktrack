      SUBROUTINE G05SKF (N, XMU, VAR, STATE, X, IFAIL)
        INTEGER            N, STATE(*), IFAIL
        double precision XMU, VAR, X(N)
        double precision gennor

        DO I = 1,N
          X(I) = gennor(XMU, SQRT(VAR))
        ENDDO

      END SUBROUTINE G05SKF
