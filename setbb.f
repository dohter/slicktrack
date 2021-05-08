C   06/05/84 602092036  MEMBER NAME  THIQAD   (S)           FORTRAN
      SUBROUTINE SETBB(X,X2,YY,T)
C
C=====Routine to set up a 7X7 THIN matrix and other things for   
C     beam-beam effects.
C     It will be assumed that the oncoming beam is centred on the 
C     electron(positron) beam. So the full 7x7 technology is not really needed. 

C     X = sigma_x, X2 = sigma_y, 
C     YY = current (amps) of oncoming protons

      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (C =299792458.D0)            ! Speed of light, metres/second
      PARAMETER (QE=1.602176D-19)            ! Elementary electric charge, Coulombs
      PARAMETER (RE=2.817940D-15)            ! Classical electron radius, metres.
      PARAMETER (ENERGE=0.510999E-3)         ! Energy equivalent of electron, GeV. 
      PARAMETER (PI=3.1415926535897932D0)    ! Pi.

      INCLUDE "csol.for"
      INCLUDE "cnlist.for"
      REAL*8 T(7,7)
C
      IF(IBMBM.EQ.0)RETURN
      IF(YY.EQ.0.D0)RETURN
      IF(X*X2.EQ.0.D0)THEN
      WRITE(53,'(A)') 
     +          ' Oncoming beam size is zero in at least one plane.'      
      STOP
      ENDIF

C    Bunch charge: Coulombs 
      QB = YY*CIRCUM/C

C    Number of oncoming particles.
      ENB = QB/QE                 ! Too big for an integer.

C    Set gamma for the central energy of the scan -- to avoid updating at
C    each energy step. 
C    NOOOOOO! that makes the tunes a bit dependent on DE0 and NSTEP. 
C    The contribution of the b-b to the total tune is small and the
C    fractional energy ranges are smallish.
C    So kill that piece of cleverness!!!
C    Instead, set gamma using ECAV so that it is independent of the 
C    starting energy too.
C      GAMMA = 0.5D0*(E00 + E00 + DE0*NSTEP)/ENERGE
C      GAMMA = (E00                )/ENERGE
      GAMMA = (ECAV                )/ENERGE
C    Derive the (gradient*length) parameters for the linear thin lens.     
      FDSIGN = IBMBM/IABS(IBMBM) 
      SX = 2.D0*RE*ENB/GAMMA/X /(X+X2) * FDSIGN
      SY = 2.D0*RE*ENB/GAMMA/X2/(X+X2) * FDSIGN
C
C    SX is -ve for unlike colliding charges (  focusing). 
C    SX is +ve for   like colliding charges (defocusing). 
C    These signs automatically give the correct signs for T(2,1) and T(4,3). 
      T(2,1) = SX
      T(4,3) = SY


      WRITE(53,'(A,6E15.6)')
     + ' Beam-beam: number of particles, gamma, T(2,1),T(4.3) ', 
     +                                    ENB,GAMMA,SX,SY
C
C
      RETURN
      END
