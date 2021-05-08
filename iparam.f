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
     +          CIRCUM,
     +          NINS,
     +          IBMBM,
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
     +          ICORR,
     +          IKMIN,
     +          IEDGE,
     +          ISEXTU,
     +          ILIN, 
     +          IEXTSIZE,
     +          SCLKIK,
     +          ECHROM,
     +          IHARM,
     +          ECAV,
     +          NITER,
     +          SFREQ,
     +          DFREQ,
     +          NFREQ,
     +          NDIFFDUMP,
     +          NSKIP  
      WRITE(53,4)KSEED1,
     +           WID1,
     +           SCUT1,
     +           KSEED2,
     +           WID2,
     +           SCUT2,
     +           KSEED3,
     +           WID3,
     +           SCUT3,
     +           KSEED4,
     +           WID4,
     +           SCUT4,
     +           KSEED5,
     +           WID5,
     +           SCUT5,
     +           KSEED6,
     +           WID6,
     +           SCUT6,
     +           KSEED7,
     +           WID7,
     +           SCUT7,
     +           KSEED8,
     +           WID8,
     +           SCUT8,
     +           KSEED9,
     +           WID9,
C     +           SCUT9,
     +           KSEED10,
     +           WID10,
     +           SCUT10
      WRITE(53,5)IBEQUIL,
     +              NPART1,
     +              NTURN1,
     +              NDAMP1,
     +              KSEED11,
     +              KSEED12,
C 
     +              IDEPLIN,
     +              NPART2,
     +              NTURN2,
     +              NDAMP2,
     +              KSEED13,
     +              KSEED14,
C
     +              IDEPNON,
     +              NPART3,
     +              NTURN3,
     +              NDAMP3,
     +              KSEED15,
     +              KSEED16,
     +              MODES,
     +              RADFAC
    2 FORMAT(/,' STARTING BEAM ENERGY      E00      =  ',1X,F9.5,/,
     +         ' ENERGY STEP SIZE          DE0      =  ',F9.4,/,
     +         ' NUMBER OF ENERGY STEPS    NSTEP    =  ',I4,/,
     +         ' APPROX. OPPOSITE CIRCUM.  CIRCUM   =  ',F9.4,/,
     +         ' NO. OF SUPERPERIODS       NINS     =  ',I4,/,
     +         '----------------------------------------------------',/,
     +         ' BEAM-BEAM                 IBMBM    =  ',I4,/, 
     +         '           F               IRAD     =  ',I4,/,
     +         '                           ISPIN    =  ',I4,/,
     +         '           L               IPTCO    =  ',I4,/,
     +         '                           IDISP    =  ',I4,/,
     +         '           A               IPTNML   =  ',I4,/,
     +         '                           IPTD     =  ',I4,/,
     +         '           G               IPTBET   =  ',I4)
    3 FORMAT( '                           ICLBP    =  ',I4,/,
     +         '           S               ITWP     =  ',I4,/,
     +         '                           ISIG     =  ',I4,/,
     +         '                           MRDISP   =  ',I4,/,
     +         '           F               ISECT    =  ',I4,/,
     +         '                           IVROT    =  ',I4,/,
     +         '           L               KICK     =  ',I4,/,
     +         '                           ICORR    =  ',I4,/,
     +         '           A               IKMIN    =  ',I4,/,
     +         '                           IEDGE    =  ',I4,/,  
     +         '           G               ISEXTU   =  ',I4,/,  
     +         '                           ILIN     =  ',I4,/,
     +         '           S               IEXTSIZE =  ',I4,/,
     +         '----------------------------------------------------',/,
     +         ' DISTORTION STRENGTH SCALE SCLKIK   =  ',F9.4,/,
     +         ' CHROMATIC ENERGY ERROR    ECHROM   =  ',F9.4,/,
     +         ' HARMONIC NUMBER           IHARM    = ',I5,/,
     +         ' FIDUCIAL ENERGY(CAVITIES) ECAV     =  ',1X,F9.5,/,
     +         ' NO. OF C.O. ITERATIONS    NITER    =  ',I4,/,
     +         ' IST. FRAC FREQUENCY SHIFT SFREQ    =  ',2X,F10.7,/,
     +         ' FRACT. RF FREQUENCY STEP  DFREQ    =  ',2X,F10.7,/,
     +         ' NO. OF RF FREQUENCY STEPS NFREQ    =  ',I4,/,
     +         ' E STEP FOR DIFFN PLOTS    NDIFFDUMP=  ',I4,/,
     +         ' SKIPPED DAMPING TIMES     NSKIP    =  ',I4,//)
    4 FORMAT(  ' RANDOM GENERATOR SEEDS    KSEED1 = ',I4,'      V. QUAD' 
     +                                                    ,' SHIFT' ,/,  
     +         ' WIDTHS                    WID1   = ',F9.4,' WIDTH' ,/,
     +         ' CUTS                      SCUT1  = ',F9.4,' CUT',/,
     +         '                           --------------------------',
     +                                                     '-------',/,
     +         '                           KSEED2 = ',I4,'      H. QUAD'
     +                                                    ,' SHIFT' ,/,
     +         '                           WID2   = ',F9.4,' WIDTH' ,/,
     +         '                           SCUT2  = ',F9.4,' CUT',/,
     +         '                           --------------------------',
     +                                                     '-------',/,
     +         '                           KSEED3 = ',I4,'      QUAD.'
     +                                                    ,' ROLE'  ,/,
     +         '                           WID3   = ',F9.4,' WIDTH' ,/,
     +         '                           SCUT3  = ',F9.4,' CUT',/,
     +         '                           --------------------------',
     +                                                     '-------',/,
     +         '                           KSEED4 = ',I4,'      QUAD.'
     +                                                   ,' STRENGTH',/,
     +         '                           WID4   = ',F9.4,' WIDTH' ,/,
     +         '                           SCUT4  = ',F9.4,' CUT',/,
     +         '                           --------------------------',
     +                                                     '-------',/,
     +         '                           KSEED5 = ',I4,'      DIPOLE'
     +                                                    ,' ROLE'  ,/,
     +         '                           WID5   = ',F9.4,' WIDTH' ,/,
     +         '                           SCUT5  = ',F9.4,' CUT',/,
     +         '                           --------------------------',
     +                                                     '-------',/,
     +         '                           KSEED6 = ',I4,'      MONITOR' 
     +                                                    ,' SHIFT' ,/,
     +         '                           WID6   = ',F9.4,' WIDTH' ,/,
     +         '                           SCUT6  = ',F9.4,' CUT',/,
     +         '                           --------------------------',
     +                                                     '-------',/,
     +         '                           KSEED7 = ',I4,'      MONITOR'
     +                                                    ,' SCALE' ,/,   
     +         '                           WID7   = ',F9.4,' WIDTH' ,/,
     +         '                           SCUT7  = ',F9.4,' CUT',/,
     +         '                           --------------------------',
     +                                                     '-------',/,
     +         '                           KSEED8 = ',I4,'      Q-M '
     +                                                    ,'OFFSET' ,/,   
     +         '                           WID8   = ',F9.4,' WIDTH' ,/,
     +         '                           SCUT8  = ',F9.4,' CUT',/,
     +         '                           --------------------------',
     +                                                     '-------',/,
     +        '                           KSEED9 = ',I4,'      MONITORS'
     +                                                   ,' ON-OFF',/,   
     +        '                           WID9   = ',F9.4,' FRACTION',/,
     +        '                           SCUT9  = ', '   ------ ',/,
     +         '                           --------------------------',
     +                                                     '-------',/,
     +         '                           KSEED10= ',I4,'             '
     +                                                    ,'      ' ,/,   
     +         '                           WID10  = ',F9.4,' WIDTH' ,/,
     +         '                           SCUT10 = ',F9.4,' CUT' //) 

   5  FORMAT(  ' TRACKING FLAGS/PARAMS     IBEQUIL= ',I6,'      M-C' 
     +                                                   ,' PHASE SPACE'
     +                                                               ,/,  
     +         '                           NPART1 = ',I6,'      PARTIC'
     +                                                   ,'LES' 
     +                                                               ,/,
     +         '                           NTURN1 = ',I6,'      M-C TUR'
     +                                                   ,'NS'
     +                                                               ,/,
     +         '                           NDAMP1 = ',I6,'      M-C DAM'
     +                                                 ,'PING FACTOR',/,
     +         '                           KSEED11= ',I6,'      CHOOSE '
     +                                                 ,'INITIAL ORBIT',
     +                                                   ' DISTRIBUTION'
     +                                                               ,/,
     +         '                           KSEED12= ',I6,'      CHOOSE '
     +                                                   ,'BIG PHOTON ',
     +                                                    'SERIES',/,
     +         '                           --------------------------',
     +                              '------------------------------',/,
     +         '                           IDEPLIN= ',I6,'      M-C'
     +                                                   ,' DEPOLARISAT'
     +                                                ,'ION:LINEAR SPIN'
     +                                                               ,/,  
     +         '                           NPART2 = ',I6,'      PARTIC'
     +                                                   ,'LES' 
     +                                                               ,/,  
     +         '                           NTURN2 = ',I6,'      M-C TUR'
     +                                                   ,'NS'
     +                                                               ,/,
     +         '                           NDAMP2 = ',I6,'      M-C DAM'
     +                                                   ,'PING FACTOR'
     +                                                               ,/,
     +         '                           KSEED13= ',I6,'      CHOOSE '
     +                                                   ,'INITIAL S-O',
     +                                                   ' DISTRIBUTION'
     +                                                               ,/,
     +         '                           KSEED14= ',I6,'      CHOOSE '
     +                                                   ,'BIG PHOTON ',
     +                                                    'SERIES',/,
     +         '                           --------------------------',
     +                              '------------------------------',/,
     +         '                           IDEPNON= ',I6,'      M-C'
     +                                                   ,' DEPOLARISAT'
     +                                                ,'ION: 3-D SPIN'
     +                                                               ,/,  
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
     +         '                           --------------------------',
     +                              '------------------------------',/,
     +         '                           MODES  = ',I6,'      DEPOLN:'
     +                                        ,'MODES SEPARATE(TOTAL) ',
     +                                                    ' 1(0)',/,
     +         '                           RADFAC =   ',F9.4 //) 
C
C
C
C
      RETURN
      END
