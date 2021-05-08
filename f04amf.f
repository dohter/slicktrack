      Subroutine f04amf (a, lda, x, ldx, b, ldb, m, n, ir, eps, qr, ldqr
     +  , alpha, e, y, z, r, ipiv, ifail)
      Integer, Intent (In) ::lda, ldx, ldb, m, n, ir, ldqr
      Integer, Intent (Inout) :: ifail
      Integer, Intent (Out):: ipiv(n)
      double precision, Intent (In):: a(lda,n), b(ldb,ir), eps
      double precision, Intent (Inout) :: x(ldx,ir), qr(ldqr,n)
      double precision, Intent (Out) :: alpha(n), e(n), y(n), z(n), r(m)

      end subroutine
