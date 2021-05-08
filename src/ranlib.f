      subroutine genmn ( parm, x, work )

c*********************************************************************72
c
cc GENMN generates a multivariate normal deviate.
c
c  Discussion:
c
c    The method is:
c    1) Generate P independent standard normal deviates - Ei ~ N(0,1)
c    2) Using Cholesky decomposition find A so that A'*A = COVM
c    3) A' * E + MEANV ~ N(MEANV,COVM)
c
c    Note that PARM contains information needed to generate the
c    deviates, and is set up by SETGMN.
c
c    PARM(1) contains the size of the deviates, P
c    PARM(2:P+1) contains the mean vector.
c    PARM(P+2:P*(P+3)/2+1) contains the upper half of the Cholesky
c    decomposition of the covariance matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    Modifications by John Burkardt.
c
c  Parameters:
c
c    Input, real PARM(P*(P+3)/2+1), parameters set by SETGMN.
c
c    Output, real X(P), a random deviate from the distribution.
c
c    Workspace, real WORK(P).
c
      implicit none

      double precision ae
      integer i
      integer icount
      integer j
      integer p
      double precision parm(*)
      double precision snorm
      double precision work(*)
      double precision x(*)

      p = int ( parm(1) )
c
c  Generate P independent normal deviates.
c
      do i = 1, p
        work(i) = snorm ( )
      end do
c
c  Compute X = MEANV + A' * WORK
c
      do i = 1, p
        icount = 0
        ae = 0.0E+00
        do j = 1, i
          icount = icount + j - 1
          ae = ae + parm(i+(j-1)*p-icount+p+1) * work(j)
        end do

        x(i) = ae + parm(i+1)

      end do

      return
      end

      function gennor ( av, sd )

c*********************************************************************72
c
cc GENNOR generates a normal random deviate.
c
c  Discussion:
c
c    This procedure generates a single random deviate from a normal distribution
c    with mean AV, and standard deviation SD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Joachim Ahrens, Ulrich Dieter,
c    Extensions of Forsythe's Method for Random
c    Sampling from the Normal Distribution,
c    Mathematics of Computation,
c    Volume 27, Number 124, October 1973, page 927-937.
c
c  Parameters:
c
c    Input, real AV, the mean of the normal distribution.
c
c    Input, real SD, the standard deviation of the normal distribution.
c
c    Output, real GENNOR, a random deviate from the distribution.
c
      implicit none

      double precision av
      double precision gennor
      double precision sd
      double precision snorm

      gennor = sd * snorm ( ) + av

      return
      end
      function genunf ( low, high )

c*********************************************************************72
c
cc GENUNF generates a uniform random deviate.
c
c  Discussion:
c
c    This procedure generates a real deviate uniformly distributed between
c    LOW and HIGH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    Modifications by John Burkardt.
c
c  Parameters:
c
c    Input, real LOW, HIGH, the lower and upper bounds.
c
c    Output, real GENUNF, a random deviate from the distribution.
c
      implicit none

      double precision genunf
      double precision high
      double precision low
      double precision r8_uni_01

      genunf = low + ( high - low ) * r8_uni_01 ( )

      return
      end
      function sdot ( n, sx, incx, sy, incy )

c*********************************************************************72
c
cc SDOT forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, real X(*), one of the vectors to be multiplied.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input, real Y(*), one of the vectors to be multiplied.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
c    Output, real SDOT, the dot product of X and Y.
c
      implicit none

      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n
      double precision sdot
      double precision stemp
      double precision sx(*)
      double precision sy(*)

      sdot = 0.0d0

      if ( n .le. 0 ) then
        return
      end if

      stemp = 0.0d0
c
c  Code for unequal increments or equal increments not equal to 1.
c
      if ( incx .ne. 1 .or. incy .ne. 1 ) then

        if ( incx .lt. 0 ) then
          ix = ( - n + 1 ) * incx + 1
        else
          ix = 1
        end if

        if ( incy .lt. 0 ) then
          iy = ( - n + 1 ) * incy + 1
        else
          iy = 1
        end if

        do i = 1, n
          stemp = stemp + sx(ix) * sy(iy)
          ix = ix + incx
          iy = iy + incy
        end do
c
c  Code for both increments equal to 1.
c
      else

        m = mod ( n, 5 )

        do i = 1, m
          stemp = stemp + sx(i) * sy(i)
        end do

        do i = m + 1, n, 5
          stemp = stemp
     &      + sx(i)     * sy(i)
     &      + sx(i + 1) * sy(i + 1)
     &      + sx(i + 2) * sy(i + 2)
     &      + sx(i + 3) * sy(i + 3)
     &      + sx(i + 4) * sy(i + 4)
        end do

      end if

      sdot = stemp

      return
      end

      subroutine setgmn ( meanv, covm, p, parm )

c*********************************************************************72
c
cc SETGMN sets data for the generation of multivariate normal deviates.
c
c  Discussion:
c
c    This procedure places P, MEANV, and the Cholesky factorization of
c    COVM in GENMN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    Modifications by John Burkardt.
c
c  Parameters:
c
c    Input, real MEANV(P), the means of the multivariate normal distribution.
c
c    Input/output, real COVM(P,P).  On input, the covariance matrix of the
c    multivariate distribution.  On output, the information in COVM has been
c    overwritten.
c
c    Input, integer P, the number of dimensions.
c
c    Output, real PARM(P*(P+3)/2+1), parameters needed to generate
c    multivariate normal deviates.
c
      implicit none

      integer p

      double precision covm(p,p)
      integer i
      integer icount
      integer info
      integer j
      double precision meanv(p)
      double precision parm(p*(p+3)/2+1)

      if ( p .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETGMN - Fatal error!'
        write ( *, '(a)' ) '  P was not positive.'
        stop 1
      end if
c
c  Store P.
c
      parm(1) = p
c
c  Store MEANV.
c
      do i = 2, p + 1
        parm(i) = meanv(i-1)
      end do
c
c  Compute the Cholesky decomposition.
c
      call spofa ( covm, p, p, info )

      if ( info .ne. 0) then
        write ( *, '(a)' ) ' INFO'
        write ( *, * ) info
        write ( *, '(a)' ) 'SETGMN - Fatal error!'
        write ( *, '(a)' ) '  SPOFA finds COVM not positive definite.'
        stop 1
      end if
c
c  Store the upper half of the Cholesky factor.
c
      icount = p + 1

      do i = 1, p
        do j = i, p
          icount = icount + 1
          parm(icount) = covm(i,j)
        end do
      end do

      return
      end
      function snorm ( )

c*********************************************************************72
c
cc SNORM samples the standard normal distribution.
c
c  Discussion:
c
c    This procedure corresponds to algorithm FL, with M = 5, in the reference.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 March 2013
c
c  Author:
c
c    Original FORTRAN77 version by Barry Brown, James Lovato.
c    Modifications by John Burkardt.
c
c  Reference:
c
c    Joachim Ahrens, Ulrich Dieter,
c    Extensions of Forsythe's Method for Random
c    Sampling from the Normal Distribution,
c    Mathematics of Computation,
c    Volume 27, Number 124, October 1973, page 927-937.
c
c  Parameters:
c
c    Output, real SNORM, a random deviate from the distribution.
c
      implicit none

      double precision a(32)
      double precision aa
      double precision d(31)
      double precision h(31)
      integer i
      double precision r8_uni_01
      double precision s
      double precision snorm
      double precision t(31)
      double precision tt
      double precision u
      double precision ustar
      double precision w
      double precision y

      save a
      save d
      save h
      save t

      data a /
     &  0.0000000d+00, 0.3917609d-01, 0.7841241d-01, 0.1177699d+00,
     &  0.1573107d+00, 0.1970991d+00, 0.2372021d+00, 0.2776904d+00,
     &  0.3186394d+00, 0.3601299d+00, 0.4022501d+00, 0.4450965d+00,
     &  0.4887764d+00, 0.5334097d+00, 0.5791322d+00, 0.6260990d+00,
     &  0.6744898d+00, 0.7245144d+00, 0.7764218d+00, 0.8305109d+00,
     &  0.8871466d+00, 0.9467818d+00, 1.009990d+00,  1.077516d+00,
     &  1.150349d+00,  1.229859d+00,  1.318011d+00,  1.417797d+00,
     &  1.534121d+00,  1.675940d+00,  1.862732d+00,  2.153875d+00 /

      data d /
     &  0.0000000d+00, 0.0000000d+00, 0.0000000d+00, 0.0000000d+00,
     &  0.0000000d+00, 0.2636843d+00, 0.2425085d+00, 0.2255674d+00,
     &  0.2116342d+00, 0.1999243d+00, 0.1899108d+00, 0.1812252d+00,
     &  0.1736014d+00, 0.1668419d+00, 0.1607967d+00, 0.1553497d+00,
     &  0.1504094d+00, 0.1459026d+00, 0.1417700d+00, 0.1379632d+00,
     &  0.1344418d+00, 0.1311722d+00, 0.1281260d+00, 0.1252791d+00,
     &  0.1226109d+00, 0.1201036d+00, 0.1177417d+00, 0.1155119d+00,
     &  0.1134023d+00, 0.1114027d+00, 0.1095039d+00 /

      data h /
     &  0.3920617d-01, 0.3932705d-01, 0.3950999d-01, 0.3975703d-01,
     &  0.4007093d-01, 0.4045533d-01, 0.4091481d-01, 0.4145507d-01,
     &  0.4208311d-01, 0.4280748d-01, 0.4363863d-01, 0.4458932d-01,
     &  0.4567523d-01, 0.4691571d-01, 0.4833487d-01, 0.4996298d-01,
     &  0.5183859d-01, 0.5401138d-01, 0.5654656d-01, 0.5953130d-01,
     &  0.6308489d-01, 0.6737503d-01, 0.7264544d-01, 0.7926471d-01,
     &  0.8781922d-01, 0.9930398d-01, 0.1155599d+00, 0.1404344d+00,
     &  0.1836142d+00, 0.2790016d+00, 0.7010474d+00 /

      data t /
     &  0.7673828d-03, 0.2306870d-02, 0.3860618d-02, 0.5438454d-02,
     &  0.7050699d-02, 0.8708396d-02, 0.1042357d-01, 0.1220953d-01,
     &  0.1408125d-01, 0.1605579d-01, 0.1815290d-01, 0.2039573d-01,
     &  0.2281177d-01, 0.2543407d-01, 0.2830296d-01, 0.3146822d-01,
     &  0.3499233d-01, 0.3895483d-01, 0.4345878d-01, 0.4864035d-01,
     &  0.5468334d-01, 0.6184222d-01, 0.7047983d-01, 0.8113195d-01,
     &  0.9462444d-01, 0.1123001d+00, 0.1364980d+00, 0.1716886d+00,
     &  0.2276241d+00, 0.3304980d+00, 0.5847031d+00 /

      u = r8_uni_01 ( )
      if ( u .le. 0.5d0 ) then
        s = 0.0d0
      else
        s = 1.0d0
      end if
      u = 2.0d0 * u - s
      u = 32.0d0 * u
      i = int ( u )
      if ( i .eq. 32 ) then
        i = 31
      end if
c
c  Center
c
      if ( i .ne. 0 ) then

        ustar = u - dble ( i )
        aa = a(i)

10      continue

        if ( t(i) .lt. ustar ) then

          w = ( ustar - t(i) ) * h(i)

          y = aa + w

          if ( s .ne. 1.0d0 ) then
            snorm = y
          else
            snorm = -y
          end if

          return

        end if

        u = r8_uni_01 ( )
        w = u * ( a(i+1) - aa )
        tt = ( 0.5d0 * w + aa ) * w

20      continue

        if ( tt .lt. ustar ) then
          y = aa + w
          if ( s .ne. 1.0d0 ) then
            snorm = y
          else
            snorm = -y
          end if
          return
        end if

        u = r8_uni_01 ( )

        if ( u .le. ustar ) then
          tt = u
          ustar = r8_uni_01 ( )
          go to 20
        end if

        ustar = r8_uni_01 ( )
        go to 10
c
c  Tail
c
      else

        i = 6
        aa = a(32)

30      continue

        u = u + u

        if ( u .lt. 1.0d0 ) then
          aa = aa + d(i)
          i = i + 1
          go to 30
        end if

        u = u - 1.0d0
        w = u * d(i)
        tt = ( 0.5d0 * w + aa ) * w

40      continue

        ustar = r8_uni_01 ( )

        if ( tt .lt. ustar ) then
          y = aa + w
          if ( s .ne. 1.0d0 ) then
            snorm = y
          else
            snorm = -y
          end if
          return
        end if

        u = r8_uni_01 ( )

        if ( u .le. ustar ) then
          tt = u
        else
          u = r8_uni_01 ( )
          w = u * d(i)
          tt = ( 0.5d0 * w + aa ) * w
        end if

        go to 40

      end if

      end
      subroutine spofa ( a, lda, n, info )

c*********************************************************************72
c
cc SPOFA factors a real symmetric positive definite matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2013
c
c  Author:
c
c    Cleve Moler
c
c  Parameters:
c
c    Input/output, real A(LDA,N).  On input, the symmetric matrix to be factored.
c    Only the diagonal and upper triangle are accessed.  On output, the strict
c    lower triangle has not been changed.  The diagonal and upper triangle contain
c    an upper triangular matrix R such that A = R' * R.  If INFO is nonzero,
c    the factorization was not completed.
c
c    Input, integer LDA, the leading dimension of the array A.
c    N <= LDA.
c
c    Input, integer N, the order of the matrix.
c
c    Output, integer INFO, error flag.
c    0, no error was detected.
c    K, the leading minor of order K is not positive definite.
c
      implicit none

      integer lda
      integer n

      double precision a(lda,n)
      integer info
      integer j
      integer jm1
      integer k
      double precision s
      double precision sdot
      double precision t

      info = 0

      do j = 1, n
        info = j
        s = 0.d0
        jm1 = j - 1
        do k = 1, jm1
          t = a(k,j) - sdot ( k-1, a(1,k), 1, a(1,j), 1 )
          t = t / a(k,k)
          a(k,j) = t
          s = s + t * t
        end do
        s = a(j,j) - s
        if ( s .le. 0.0d0 ) then
          info = j
          return
        end if
        a(j,j) = sqrt ( s )
      end do

      info = 0
      return
      end
