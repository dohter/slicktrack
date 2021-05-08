C   13/09/83 205221650  MEMBER NAME  IPARAM   (MAY92.S)     FORTRAN
      SUBROUTINE SFPARAM
C
C
C=====ROUTINE TO PRINT NAMELIST PARAMETERS IN A CIVILISED FORMAT
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE "csfnlist.for"
C
C
      WRITE(53,1)
    1 FORMAT(///,' INPUT NAMELIST PARAMETERS')

      WRITE(53,2)E00,
     +              DE0,
     +              NSTEP,
     +              QS, 
     +              DCONST,
     +              CFAC,
     +              CIRC, 
     +              CURVX,
     +              IEXTSIZE,
     +              NSECT
      WRITE(53,3)NPART3,
     +              NTURN3,
     +              NDAMP3,
     +              KSEED15,
     +              KSEED16,
     +              RADFAC

    2 FORMAT(/,' STARTING BEAM ENERGY      E00         =  ',1X,F9.5,/,
     +         ' ENERGY STEP SIZE          DE0         =  ',F9.4,/,
     +         ' NUMBER OF ENERGY STEPS    NSTEP       =  ',I4,/,
     +         ' SYNCHROTRON TUNE          QS          =  ',F9.4,/,
     +         ' DAMPING CONSTANT          DCONST      =  ',F9.4,/,
     +         ' COMPACTION FACTOR         CFAC        =  ',F9.4,/,
     +         ' CIRCUMFERENCE             CIRC        =  ',F9.4,/,
     +         ' COMMON DIPOLE CURVATURE   CURVX       =  ',F9.4,/,
     +         ' EXTERNAL BEAM SIZE        IEXTSIZE    =  ',I4,/,
     +         ' NUMBER OF SECTIONS        NSECT       =  ',I4,/)

    3 FORMAT(  ' TRACKING FLAGS/PARAMS' ,/,
     +         '                           NPART3 = ',I6,'      PARTIC'
     +                                                   ,'LES' 
     +                                                               ,/,  
     +         '                           NTURN3 = ',I6,'      M-C TUR'
     +                                                   ,'NS'
     +                                                               ,/,
     +         '                           NDAMP3 = ',I6,'      M-C DAM'
     +                                                   ,'PING FACTOR'
     +                                                               ,/,
     +         '                           KSEED15= ',I6,'      CHOOSE '
     +                                                   ,'INITIAL S-O',
     +                                                   ' DISTRIBUTION'
     +                                                               ,/,
     +         '                           KSEED16= ',I6,'      CHOOSE '
     +                                                   ,'BIG PHOTON ',
     +                                                    'SERIES',/,
     +         '                           RADFAC =   ',F9.4 //) 
C
C
C
C
      RETURN
      END
