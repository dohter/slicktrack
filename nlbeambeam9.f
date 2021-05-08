C   15/11/79 401171758  MEMBER NAME  MX88     (MAY92.S)     FORTRAN
C
      SUBROUTINE NLBEAMBEAM9(II,X,X2,YY,S,ZWBB,SORBVEC,
     +                                       SPINKICKA,NSOLBB,MODIBMBM)

C                             II,XXBB(II),X2BB(II),YYBB(II),NU,ZWBB,SORBVEC,
C     +                                   SPINKICKA,NSOLBB(II),MODIBMBM)

C     X  = sigma_x, X2 = sigma_y, 
C     YY = current (amps) of oncoming beam.


C=====GENERATE THE THIN LENS NONLINEAR BEAM-BEAM KICKS USING THE SPIN BASIS.
C     USE EITHER ELLIPTICAL  GAUSSIAN BEAMS (IBMBM = +/-2) WITH EXACT CALC OR
C                ROUND       GAUSSIAN BEAMS (IBMBM = +/-3) OR 
C                ELLIPTICAL  GAUSSIAN BEAMS (IBMBM = +/-4) WITH A QUICK + DIRTY
C                                                          CALC. 
C

      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE "cnlist.for" 

      REAL*8 LXI,LZI,LSI,MXI,MZI,MSI,NXI,NZI,NSI

      PARAMETER (NPART   = 10000) ! Maximum allowed particles. Use NPART3 of them.
      PARAMETER (LIMSECT = 3000)  ! Maximum allowed sections.
      PARAMETER (C =299792458.D0)            ! Speed of light, metres/second
      PARAMETER (QE=1.602176D-19)            ! Elementary electric charge, Coulombs
      PARAMETER (RE=2.817940D-15)            ! Classical electron radius, metres.
      PARAMETER (ENERGE=0.510999E-3)         ! Energy equivalent of electron, GeV. 
      PARAMETER (PI=3.1415926535897932D0)    ! Pi.

      DIMENSION ZWBB(3,3,LIMSECT)
      DIMENSION SORBVEC(8,NPART)
      DIMENSION SPINKICKA(3,NPART)

C===Complex error function stuff.
      COMPLEX*16 Z1,Z2,W1,W2,DELTA
      COMPLEX*16 S15DDF
      EXTERNAL   S15DDF
C
C
      IF(YY.EQ.0.D0)RETURN
      IF(X*X2.EQ.0.D0)THEN
      WRITE(53,103)
      WRITE(53,103)
  103 FORMAT(/,'  ')
      WRITE(53,'(A,A,A)') 
     +  ' In NLBEAMBEAM9: ',
     +  ' the oncoming beam size is zero in at least one plane.',      
     +  '  So STOP.'
      STOP
      ENDIF

C=====Assume that the beams are relatively centred
C=====so that CLOSED ORBIT SHIFTS can be ignored.

C=====SET SPIN VECTORS---INPUT VALUES.
      LXI=ZWBB(1,3,II)
      LZI=ZWBB(2,3,II)
      LSI=ZWBB(3,3,II)
      MXI=ZWBB(1,2,II)
      MZI=ZWBB(2,2,II)
      MSI=ZWBB(3,2,II)
      NXI=ZWBB(1,1,II)
      NZI=ZWBB(2,1,II)
      NSI=ZWBB(3,1,II)

C    Bunch charge: Coulombs 
      QB = YY*CIRCUM/C
C    Number of oncoming particles.
      ENB = QB/QE                 

C    Set gamma as in SETBB so that it's consistent.
      GAMMA = ECAV/ENERGE

      CONST = 2.D0*RE*ENB/GAMMA*IBMBM/MODIBMBM
      AGP1=(1.D0+S)



      IF(MODIBMBM.EQ.2)THEN
C===Get spin-orbit transformation for elliptical Gaussian beams.
C   Use the formulae for the orbit deflections in J. Kewisch's thesis 
C   Even with the mixture of magnetic and electric fields, the spin kick can 
C   be derived from the orbit kick in the usual way at high energy and 
C   with the trajectory along the longitudinal axis.
C   It is assumed that the horizontal sigma is greater than the vertical sigma.

      IF(X.LE.X2)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,A,A)')
     +  ' In NLBEAMBEAM9: ',
     +  'the horizontal size of the oncoming beam is less than or the', 
     +  ' same as the vertical size. So STOP.'
      STOP
      ENDIF


      DEN = 1.D0/DSQRT(2.D0*(X**2 - X2**2))
      R1 = X2/X
      R2 = X/X2
      S1 = 0.5D0/X**2
      S2 = 0.5D0/X2**2
      RPI = DSQRT(PI) 
      CONSTE = CONST*RPI*DEN

      IFAIL = 0
      DO 9 IP = 1,NPART3 
      RSQ1     = SORBVEC(1,IP)**2 *1.D0                             ! *0 to switch off hor. position.
      RSQ2     = SORBVEC(3,IP)**2 
      RSQ      = RSQ1 + RSQ2
      IF(RSQ.EQ.0.D0)GO TO 9

      Z1 = DEN * DCMPLX(SORBVEC(1,IP)*1.D0,   SORBVEC(3,IP)   )     ! *0 to switch off hor. position.
      Z2 = DEN * DCMPLX(SORBVEC(1,IP)*R1*1.D0,SORBVEC(3,IP)*R2)     ! *0 to switch off hor. position.
      EXPFAC =  DEXP(-RSQ1*S1 - RSQ2*S2)

      W1 = S15DDF(Z1, IFAIL)
      W2 = S15DDF(Z2, IFAIL)    

      DELTA  = (-1.D0) * CONSTE * (EXPFAC * W2 - W1) 
      COEFF1 = DIMAG(DELTA)         ! Get the imaginary part
      COEFF2 = DELTA                ! Get the real part
 

      ACOEFF1 = AGP1*COEFF1
      ACOEFF2 = AGP1*COEFF2

      SORBVEC(2,IP) = SORBVEC(2,IP) + COEFF1           !*10.D0
      SORBVEC(4,IP) = SORBVEC(4,IP) + COEFF2           !*10.D0
C
      IF(NSOLBB.GT.0)GO TO 9
      SPINKICKA(1,IP)=   ACOEFF1*LZI 
     +                 - ACOEFF2*LXI
      SPINKICKA(2,IP)= - ACOEFF1*MZI 
     +                 + ACOEFF2*MXI
      SPINKICKA(3,IP)=   ACOEFF1*NZI 
     +                 - ACOEFF2*NXI

    9 CONTINUE  

C   Check that the complex error function stuff was clean for all particles.   
      IF(IFAIL.NE.0)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,A,A,I3)')
     +  ' In NLBEAMBEAM: ',
     +  'the complex error function flagged a problem for a particle.',
     +  ' So STOP. IFAIL= ',IFAIL
      STOP
      ENDIF


      ELSEIF(MODIBMBM.EQ.3)THEN
C=====Round Gaussian beams.     
      TXSQ = 2.D0*X**2       !The size of the oncoming beam is given by the hor. sigma.

      DO 10 IP = 1,NPART3 
      RSQ           = SORBVEC(1,IP)**2 + SORBVEC(3,IP)**2
      IF(RSQ.EQ.0.D0)GO TO 10
      COEFF         = CONST*(1.D0 - DEXP(-RSQ/TXSQ))/RSQ
      ACOEFF = AGP1*COEFF
      SORBVEC(2,IP) = SORBVEC(2,IP) + COEFF*SORBVEC(1,IP)
      SORBVEC(4,IP) = SORBVEC(4,IP) + COEFF*SORBVEC(3,IP)
C
      IF(NSOLBB.GT.0)GO TO 10
      SPINKICKA(1,IP)=   ACOEFF*( SORBVEC(1,IP)*LZI - SORBVEC(3,IP)*LXI)
      SPINKICKA(2,IP)=   ACOEFF*(-SORBVEC(1,IP)*MZI + SORBVEC(3,IP)*MXI)
      SPINKICKA(3,IP)=   ACOEFF*( SORBVEC(1,IP)*NZI - SORBVEC(3,IP)*NXI)
   10 CONTINUE  
C
      ELSEIF(MODIBMBM.EQ.4)THEN
C=====A quick + dirty treatment of elliptical Gaussian beams.
C     The vertical   kick is taken to be indep. of the x position and
C     the horizontal kick is taken to be indep. of the y position. 
C     The distance scales are set by the parameters of the usual line linear 
C     approximation.
C     So there is no indirect coupling between planes and the effect on
C     each plane is exaggerated because the fall off, in the other plane,
C     of the field at the edge of the beam is ignored.  


      TXSQ1 = X *(X + X2)
      TXSQ2 = X2*(X + X2)

      DO 11 IP = 1,NPART3 
      RSQ1          = SORBVEC(1,IP)**2
      RSQ2          = SORBVEC(3,IP)**2
      RSQ           = RSQ1 + RSQ2

      IF(RSQ.EQ.0.D0)GO TO 11
      COEFF1         = CONST*(1.D0 - DEXP(-RSQ1/TXSQ1))/RSQ1
      COEFF2         = CONST*(1.D0 - DEXP(-RSQ2/TXSQ2))/RSQ2
 
C      COEFF1         = CONST/TXSQ1
C      COEFF2         = CONST/TXSQ2

      ACOEFF1 = AGP1*COEFF1
      ACOEFF2 = AGP1*COEFF2

      SORBVEC(2,IP) = SORBVEC(2,IP) + COEFF1*SORBVEC(1,IP)    !*1000.D0
      SORBVEC(4,IP) = SORBVEC(4,IP) + COEFF2*SORBVEC(3,IP)    !*1000.D0
C


      IF(NSOLBB.GT.0)GO TO 11
      SPINKICKA(1,IP)=   ACOEFF1*SORBVEC(1,IP)*LZI 
     +                 - ACOEFF2*SORBVEC(3,IP)*LXI
      SPINKICKA(2,IP)= - ACOEFF1*SORBVEC(1,IP)*MZI 
     +                 + ACOEFF2*SORBVEC(3,IP)*MXI
      SPINKICKA(3,IP)=   ACOEFF1*SORBVEC(1,IP)*NZI 
     +                 - ACOEFF2*SORBVEC(3,IP)*NXI
   11 CONTINUE  


      ENDIF


      RETURN
      END


C======Guide for geting spin kicks.
C      IF(NNSOL.EQ.1)RETURN
C      A(7,1)= AGP1* A(2,1)      *LZI
C      A(8,1)=-AGP1* A(2,1)      *MZI
C      A(7,3)=-AGP1* A(4,3)      *LXI     
C      A(8,3)= AGP1* A(4,3)      *MXI   
