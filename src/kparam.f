C   13/09/83 205221650  MEMBER NAME  IPARAM   (MAY92.S)     FORTRAN
      SUBROUTINE IPARAM
C
C
C=====ROUTINE TO PRINT NAMELIST PARAMETERS IN A CIVILISED FORMAT
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE "cnlist.for"
C
C
C!!!!!!!!!!!!!!!!!!!!!!!!!!Some nonsense: this NINS is not the one read in!!!!!!!!!!!!
      WRITE(53,1)
    1 FORMAT(///,' INPUT NAMELIST PARAMETERS')
      WRITE(53,2)E00,
     +          DE0,
     +          NSTEP,
     +          NINS,
     +          IRAD,
     +          ISPIN,
     +          IPTCO,
     +          IDISP,
     +          IPTNML,
     +          IPTD,
     +          IPTBET
      WRITE(53,3)ICLBP,
     +          ITWP,
     +          ISIG,
     +          MRDISP,
     +          ISECT,
     +          IVROT,
     +          KICK,
     +          SCLKIK,
     +          ECHROM,
     +          IHARM,
     +          ECAV,
     +          ITER,
     +          SFREQ,
     +          DFREQ,
     +          NFREQ
    2 FORMAT(/,' STARTING BEAM ENERGY      E00    = ',1X,F9.5,/,
     +         ' ENERGY STEP SIZE          DE0    = ',F9.4,/,
     +         ' NUMBER OF ENERGY STEPS    NSTEP  = ',I4,/,
     +         ' NO. OF SUPERPERIODS       NINS   = ',I4,/,
     +         '                           IRAD   = ',I4,/,
     +         '                           ISPIN  = ',I4,/,
     +         '           F               IPTCO  = ',I4,/,
     +         '                           IDISP  = ',I4,/,
     +         '           L               IPTNML = ',I4,/,
     +         '                           IPTD   = ',I4,/,
     +         '           A               IPTBET = ',I4)
    3 FORMAT(  '                           ICLBP  = ',I4,/,
     +         '           G               ITWP   = ',I4,/,
     +         '                           ISIG   = ',I4,/,
     +         '           S               MRDISP = ',I4,/,
     +         '                           ISECT  = ',I4,/,
     +         '                           IVROT  = ',I4,/,
     +         '                           KICK   = ',I4,/,
     +         ' KICKER STRENGTH SCALE     SCLKIK = ',F9.4,/,
     +         ' CHROMATIC ENERGY ERROR    ECHROM = ',F9.4,/,
     +         ' HARMONIC NUMBER           IHARM  = ',I7,/,
     +         ' FIDUCIAL ENERGY(CAVITIES) ECAV   = ',1X,F9.5,/,
     +         ' NO. OF C.O. ITERATIONS    ITER   = ',I4,/,
     +         ' IST. FRAC FREQUENCY SHIFT SFREQ  = ',2X,F10.7,/,
     +         ' FRACT. RF FREQUENCY STEP  DFREQ  = ',2X,F10.7,/,
     +         ' NO. OF RF FREQUENCY STEPS NFREQ  = ',I4//)
C
C
C
C
      RETURN
      END
