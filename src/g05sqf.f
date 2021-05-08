      SUBROUTINE G05SQF (N, A, B, STATE, X, IFAIL)
        INTEGER            N, STATE(*), IFAIL, I
        double precision A, B, X(N)
        double precision r8_uni_01

        do I=1,N
          X(I) = (B-A) * r8_uni_01() + A
        enddo

      END SUBROUTINE G05SQF
