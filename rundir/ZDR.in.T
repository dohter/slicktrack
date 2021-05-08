************************* EXPLANATION ****************************************
*       ROTATOR      => INPUT FILE FOR THE HERA ROTATOR DATA                 *
*       RING         => INPUT DATA FOR THE LATTICE                           *
*       OUTPUT       => FILE FOR TOTAL OUTPUTK                               *
*       POLARI       => LIST OF POLARIZATION AND ENERGY                      *
*       AZIFUNCS     => OPTICAL FUNCTIONS AROUND THE RING                    *
*       GAUSSES      => CHECK OF GAUSSIAN DISTRIBUTIONS                      *
*       ELLIPSE1     => BEAM ELLIPSE FROM THEORETICAL COVARIANCE MATRIX      *
*       ELLIPSE1     => BEAM ELLIPSE FROM TRACKED ENSEMBLE                   *
*       BEAMDIST     => DISTRIBUTION OF TRACKED ENSEMBLE                     *
*       SPINDIST     => DISTRIBUTION OF SPIN ALPHA & BETA                    *
*       SPINDIFF     => TIME DEPENDENCE OF THE SPIN SPREAD                   *
*       ISPINSPREAD  => EQUILIBRIUM SPIN  ALPHA & BETA COVARIANCES           *
*       MCTDEP       => M-C DEPOLARISATION TIMES                             *
*       ABRATIO      => RATIO OF 1ST TURN SPIN ALPHA & BETA COVARIANCES      *
*       SPINELLIPSE1 => SLICK CONSTANT DENSITY CONTOUR FOR ALPHA AND BETA    *
*       SPINELLIPSE2 => M_C   CONSTANT DENSITY CONTOUR FOR ALPHA AND BETA    *
*       FFT          => FAST FOURIER TRANSFORM OF ORBITAL MOTION             *
******************************************************************************

 &FILES ROTATOR     ='../slick_rotator.inp',
        RING        ='ZDR.slick',
        OUTPUT      ='ZDR_results/slicktrack_PPPP.out',
        POLARI      ='ZDR_results/pol_PPPP.out',
        AZIFUNCS    ='ZDR_results/funcs_PPPP.out',
        GAUSSES     ='ZDR_results/gausses_PPPP.out',
        ELLIPSE1    ='ZDR_results/ellipse1_PPPP.out',
        ELLIPSE2    ='ZDR_results/ellipse2_PPPP.out',
        BEAMDIST    ='ZDR_results/beamdist_PPPP.out',
        SPINDIST    ='ZDR_results/spindist_PPPP.out',
        SPINDIFF    ='ZDR_results/spindiff_PPPP.out',
        MCTDEP      ='ZDR_results/mctdep_PPPP.out',
        ABRATIO     ='ZDR_results/ABratio_PPPP.out',
        SPINELLIPSE1='ZDR_results/spinellipse1_PPPP.out',
        SPINELLIPSE2='ZDR_results/spinellipse2_PPPP.out',
        FFT         ='ZDR_results/fft_PPPP.out',
        ISPINSPREAD ='ZDR_results/ispinspread_PPPP.out'

 &END
*********************** EXPLANATION *************************************
*       E00    => STARTING BEAM ENERGY                                  *
*       ECAV   => FIDUCIAL ENERGY(CAVITIES)                             *
*       DE0    => ENERGY STEP SIZE                                      *
*       NSTEP  => >0 NUMBER OF ENERGY STEPS                             *
*************************************************************************
*       CIRCUM => APPROX.CIRCUMFERENCE OF OTHER RING TO CONVERT mA to   *
*                 PARTICLES/BUNCH in B-B CALCS.                         *
*       IBMBM  => +/-  0/1/2/3/4  NO                     BEAM-BEAM      *
*                                 LINEAR                 BEAM-BEAM      *
*                                 ELLIPTICAL NONLINEAR   BEAM-BEAM      *
*                                 ROUND B-B              BEAM-BEAM      *
*                                 CHEEP-DIRTY ELLIPTICAL BEAM-BEAM      *
*                                                                       *
*************************************************************************
*       IPTCO  => 0/1/2 SKIPS/PRINTS/ONLY ONE CLOSED ORBIT              *
*       KICK   => Applies distortions                                   *
*       SCLKIK => DISTORTION STRENGTH SCALE FOR SENSITIVITY CHECKS      *
*       NITER  => NO. OF C.O. ITERATIONS                                *
*       IRAD   => 1/0 INCLUDES/NEGLECTS ENERGY JUMPS IN ELEMENTS        *
*       ICORR  => 1/0 APPLIES/NEGLECTS CONVENTIONAL ORBIT CORRECTION    *
*       IKMIN  => 1/0 APPLIES/NEGLECTS KICK MINIMISATION                *
*************************************************************************
*       ISPIN  => 1 INCLUDES POLARIZATION CALCULATIONS                  *
*       IDISP  => 1 CALCULATES THE PERTURBED DISPERSION FUNCTIONS       *
*       IPTNML => 1 PRINTS OUT SPIN BASE VECTORS                        *
*       IPTD   => 1 PRINTS OUT "D" VECTORS AT DIPOLES WITH NTWIST=1     *
*       IPTBET =>                                                       *
*       ICLBP  => 1 DUMP OUT THE 6-DIM CLOSED ORBIT AND STOP            *
*       ITWP   => 1 DUMP OUT THE TWISS PARAMETERS AND STOP              *
*       ISIG   => 1/2 BEAM SIZES AROUND THE RING/AT INTERSECTIONS       *
*       ISECT  => 1 GET 8X8 MATRIX FOR A SECTION OF RING                *
*       IEDGE  => 1 Apply edge focussing to dipoles                     *
*       ISEXTU => 1 Switch on sextupoles                                *
*       ILIN   => 1 Run LINOPT                                          *
*       IEXTEM => 1 Use external initial emittances for tracking M-C    *
*************************************************************************
*       MRDISP => USE MAIS & RIPKEN 'DISPERSION FORMALISM' IN SPIN      *
*       IVROT  => 0/1 ROTATORS FLAT/ACTIVE                              *
*       IHARM  => HARMONIC NUMBER                                       *
*************************************************************************
*       SFREQ  => IST. FRAC FREQUENCY SHIFT                             *
*       DFREQ  => FRACT. RF FREQUENCY STEP                              *
*       NFREQ  => NO. OF RF FREQUENCY STEPS                             *
*       NDIFFDUMP => ENERGY STEP FOR DIFFUSION PLOTS                    *
*       NSKIP  => NO. OF SKIPPED DAMPING TIMES BEFORE FITTING A SLOPE   *
*************************************************************************
*       NINS   => NO. OF SUPERPERIODS                                   *
*************************************************************************
*       KSEED1,WID1,SCUT1=> VERTICAL QUAD SHIFTS, WIDTH, GAUSSIAN CUT   *
*       KSEED2,WID2,SCUT2=> HORIZONTAL QUAD SHIFTS, WIDTH, GAUSSIAN CUT *
*       KSEED3,WID3,SCUT3=> QUAD ROLE, WIDTH, GAUSSIAN CUT              *
*       KSEED4,WID4,SCUT4=> FRACTIONAL QUAD GRADIENT ERROR, WIDTH,G.CUT *
*       KSEED5,WID5,SCUT5=> DIPOLE ROLE, WIDTH, GAUSSIAN CUT            *
*       KSEED6,WID6,SCUT6=> MONITOR SHIFT, WIDTH, GAUSSIAN CUT          *
*       KSEED7,WID7,SCUT7=> MONITOR SCALE, WIDTH, GAUSSIAN CUT          *
*       KSEED8,WID8,SCUT8=> QUAD-MONITOR OFFSET FOR KICK MINIMISATION   *
*       KSEED9,WID9,SCUT9=> MONITORS ON/OFF, FRACTION                   *
*       KSEED10,WID10,SCUT10=>                                          *
*************************************************************************
*       IBEQUIL=> TRACK TO GET EQUILIBRIUM PHASE SPACE                  *
*       NPART1 => NUMBER OF PARTICLES                                   *
*       NTURN1 => NUMBER OF (ARTIFICIAL) DAMPING TIMES                  *
*       NDAMP1 => ARTIFICIAL DAMPING FOR TRACKING                       *
*       KSEED11=> SEED FOR ORBIT EQUILIBRIUM. FOR 0 START ON C.O.       *
*       KSEED12=> SEED FOR BIG PHOTONS                                  *
*************************************************************************
*       IDEPLIN=> LINEAR S-O (8X8)TRACKING TO MEASURE DEPOLARISATION    *
*       NPART2 => NUMBER OF PARTICLES                                   *
*       NTURN2 => NUMBER OF (ARTIFICIAL) DAMPING TIMES.                 *
*                 IF NTURN2 = 1, SET THE LOOPS IN SCRURITA2 TO DO JUST  *
*                 1 TURN TO CHECK AGAINST 1-TURN SLICK DIFFUSION.       *
*       NDAMP2 => ARTIFICIAL DAMPING FOR TRACKING                       *
*       KSEED13=> SEED FOR S-O EQUILIBRIUM. FOR 0 START ON C.O.         *
*       KSEED14=> SEED FOR BIG PHOTONS                                  *
*************************************************************************
*       IDEPNON=> LINEAR ORBIT + 3-D SPIN TRACKING TO MEASURE DEPOLN    *
*       NPART3 => NUMBER OF PARTICLES                                   *
*       NTURN3 => NUMBER OF (ARTIFICIAL) DAMPING TIMES.                 *
*                 IF NTURN2 = 1, SET THE LOOPS IN SCRURITA2 TO DO JUST  *
*                 1 TURN TO CHECK AGAINST 1-TURN SLICK DIFFUSION.       *
*       NDAMP3 => ARTIFICIAL DAMPING FOR TRACKING                       *
*       KSEED15=> SEED FOR S-O EQUILIBRIUM. FOR 0 START ON C.O.         *
*       KSEED16=> SEED FOR BIG PHOTONS                                  *
*************************************************************************
*       MODES  => 1/0: DEPN TIME for SEPARATE/TOTAL ORBITAL MODES       *
*       RADFAC         SCALE FOR SIZE OF BIG PHOTONS                    *
*************************************************************************
&RUNSPEC E00=EEEE,ECAV=9.714380556603,DE0=0.002,CIRCUM=3800.,IBMBM=-0,NSTEP=NNN,
         IPTCO=1,KICK=1,SCLKIK=1.0,NITER=3,ICORR=1,IKMIN=0,IRAD=0,
         IPTNML=1,IPTD=0,ISPIN=1,IDISP=1,ICLBP=1,ISECT=0,ISIG=2,ITWP=0,
         IHARM=2040,IVROT=0,MRDISP=0,IEDGE=0,ISEXTU=0,ILIN=0,NINS=1,
         NDIFFDUMP=1,NSKIP=5,IEXTSIZE=0,
         KSEED1=8,WID1=0.3, SCUT1=3.0,
         KSEED2=0,WID2=0.0, SCUT2=0.0,
         KSEED3=6,WID3=0.3, SCUT3=3.0,
         KSEED4=0,WID4=1.0, SCUT4=1.0,
         KSEED5=1,WID5=0.3, SCUT5=3.0,
         KSEED6=4,WID6=0.3, SCUT6=1.0,
         KSEED7=5,WID7=0.0, SCUT7=1.0,
         KSEED8=7,WID8=0.1, SCUT8=1.0,
         KSEED9=2,WID9=0.0, SCUT9=0.0,
         IBEQUIL=0, NPART1=1000, NTURN1=20,  NDAMP1=5,
         IDEPLIN=0, NPART2=100, NTURN2=200 NDAMP2=5, KSEED13=3, KSEED14=3,
         IDEPNON=1, NPART3=100, NTURN3=300, NDAMP3=5, KSEED15=3, KSEED16=4,
         MODES=0,RADFAC=1.D0,
&END


*****&RUNSPEC E00=9.2500,ECAV=9.714380556603,DE0=0.058,CIRCUM=3800.,IBMBM=-0,NSTEP=2,

















****NSTEP=1330
****&SL E00=9.250000000000,ECAV=9.714380556603,DE0=0.001,NSTEP=1,

****E00=9.2820
****E00=9.320000000000
****E00=9.354000000000
****E00=9.250000000000
****E00=9.315


****&SL E00=9.714380556603,ECAV=9.714380556603,DE0=0.01,NSTEP=70,
    IPTCO=1,KICK=0,SCLKIK=1.0,ITER=3,ICORR=1,IRAD=0,
    IPTNML=1,IPTD=1,ISPIN=1,IDISP=1,ICLBP=1,ISECT=0,ISIG=2,
    IHARM=1170,IVROT=0,MRDISP=0,
    KSEED1=8,WID1=0.3,SCUT1=3.0,
    KSEED2=0,WID2=0.0,SCUT2=0.0,
    KSEED3=6,WID3=0.3,SCUT3=3.0,
    KSEED4=7,WID4=1.0,SCUT4=1.0,
    KSEED5=1,WID5=1.0,SCUT5=3.0,
    KSEED6=4,WID6=0.1,SCUT6=1.0,
    KSEED7=5,WID7=0.5,SCUT7=1.0,
    IBEQUIL=0, NPART1=1, NTURN1=1000,  NDAMP1=5,
    IDEPLIN=1, NPART2=1, NTURN2=1000,  NDAMP2=1,
****&END






*KSEED1=3 gives big corrector values.
*&SL E00=9.714380556603,DE0=0.010000,NSTEP=70,ITER=3,IHARM=1170,
*  IPTCO=1,KICK=1,SCLKIK=+1.000000,IPTNML=1,IPTD=1,IRAD=0,ISIG=2,ISPIN=1,
*  ECAV=9.714380556603,ISECT=0,IDISP=1,MRDISP=0,IVROT=0,  &END

*&SL E00=14.5715708349045,DE0=0.010000,NSTEP=100,ITER=3,IHARM=1170,
*  IPTCO=1,KICK=1,SCLKIK=+1.000000,IPTNML=1,IPTD=1,IRAD=0,ISIG=2,ISPIN=1,
*  ECAV=9.714380556603,ISECT=0,IDISP=1,MRDISP=0,IVROT=0,  &END

*SL E00=9.514380556603,DE0=0.005000,NSTEP=200,ITER=3,IHARM=1170,
* IPTCO=1,KICK=0,SCLKIK=+1.00000,IPTNML=1,IPTD=1,IRAD=0,ISIG=2,ISPIN=1,
*  ECAV=9.714380556603,ISECT=0,IDISP=1,MRDISP=0,IVROT=0,  &END



*&SL E00=5.014380556603,DE0=0.010000,NSTEP=1200,ITER=3,IHARM=1170,
*  IPTCO=1,KICK=0,SCLKIK=+0.0450,IPTNML=1,IPTD=1,IRAD=0,ISIG=2,ISPIN=1,
*  ECAV=9.714380556603,ISECT=0,IDISP=1,MRDISP=0,IVROT=0,  &END
