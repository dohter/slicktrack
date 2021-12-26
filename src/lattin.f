C   16/07/81 509051828  MEMBER NAME  LATTIN   (SEPT95.S)    FORTRANg1(
      SUBROUTINE LATTIN(RLENGE,IROT)
C
C                     READS THE SLIM INPUT DATA
C
C                     D.BARBER 1982 -->
C
C
C    DATA LAYOUT IS THE SAME AS FOR THE THIN LENS PROGRAM BUT IS
C    USED FOR A SLIGHTLY DIFFERENT LATTICE LIST ORGANISATION.
C                      D.BARBER   JUNE /84
C
C    FORMAT=4:  4 MAGNET POSITIONS/LINE & 1/10 MM. UNITS
C    FORMAT=5:  5 MAGNET POSITIONS/LINE & 1/10 MM. UNITS
C    FORMAT=6:  6 MAGNET POSITIONS/LINE &    1 MM. UNITS
C
C-------------------------------------------------------------------------------
C    E0            Initial energy (GeV)
C    DE0           Energy steps
C    ECAV          ``Fiducial energy'' for fixing the synch. tune.
C    CIRCUM        Approx. circum. of other ring
C                  to convert mA to particles/bunch in b-b calcs.
C    IHARM         Harmonic number.
C    SFREQ         First frac. rf frequency shift --> change damping (thin lens:SLIM)
C    DFREQ         Frac. rf frequency step      s --> change damping (thin lens:SLIM)
C    NFREQ         Number of rf frequency steps.
C    NSTEP         Number of energy steps.
C    ECHROM        Obsolescent: adds chromatic effects to quads.
C    NINS          Obsolescent: number of superperiods for time saving in LINOPT.
C    NSKIP         Number of skipped damping times before fitting slope.
C    NDIFFDUMP     Energy step for dumping diffusion plots.
C    NITER         Number of iterations to get the closed orbit.
C    SCLKIK        Scale strengths of kicks which simulate misalignments.
C------------------------------------------------------------------------
C   (IBMBM.EQ.0    ) Ignore beam-beam effects when b-b (type 17) elements are present.
C   ( "   .EQ.+/- 1) Linear b-b effects:   +/- for   like/unlike charges.
C   ( "   .EQ.+/- 2) Nonlinear b-b effects for M-C, elliptical oncoming beam.
C   ( "   .EQ.+/- 3) Nonlinear b-b effects for M-C, round      oncoming beam.

C   ( "  .EQ.0)    IGNORES ENERGY JUMPS AND COMPUTES THE CLOSED ORBIT
C                  ONLY ONCE BEFORE CALLING CPSIN IN A SEPARATE LOOP.
C   (IPTCO.EQ.0)   SKIPS CLOSED ORBIT CALCULATION AND PUTS IN ZERO DEV'N
C   (  "      1)   GETS CLOSED-ORBIT 6-VECTOR AND PRINTS IT.
C   (  "      2)   CALCULATES LINEAR OPTICS TO GET SOME PARAMETERS & THEN
C   (  "      2)   SETS C.O.TO ZERO & JUMPS DIRECTLY TO THE SPIN
C   (  "      2)   WITHOUT GETTING THE PERTURBED 6x6 OPTICS AND THE EMITTANCE.
C   (ISIG.EQ.1)    CALCULATES THE BEAM SIZES AROUND THE RING.
C   (ISIG.EQ.2)    CALCULATES THE BEAM SIZES AT INTERSECTION POINTS.
C   (ISPIN .EQ.1)  INCLUDES POLARIZATION CALCULATIONS.
C   (IDISP .EQ.1)  CALCULATES THE PERTURBED DISPERSION FUNCTIONS.
C   (IPTNML.EQ.1)  PRINTS OUT SPIN BASE VECTORS.
C   (IPTD  .EQ.1)  PRINTS OUT "D" VECTORS AT DIPOLES WITH NTWIST=1
C   (KICK  .NE.0)  APPLIES DISTORTIONS.
C   (ICLBP.EQ.1)   DUMP OUT THE 6-DIM CLOSED ORBIT.
C   (ICLBP.EQ.2)   DUMP OUT THE 6-DIM CLOSED ORBIT & STOP
C   (ITWP .EQ.1)   DUMP OUT THE TWISS PARAMETERS AND STOP
C   (ISECT.EQ.1)   TREAT THE LATICE AS AN OPEN BEAM LINE.
C   (MRDISP.EQ.1)  USE MAIS & RIPKEN 'DISPERSION FORMALISM' IN CSPIN.
C   (IVROT.EQ.1)   VARY ROTATOR SETTING WITH ENERGY TO KEEP N VERTICAL.
C   (ICORR.NE.0)   CORRECT THE CLOSED ORBIT.
C   (IKMIN.NE.0)   Apply kick minimisation to closed orbit.
C    IIEDGE.NE.0)  Apply edge focussing to dipoles.
C   (ISEXTU.NE.0)  Switch on sextupoles.
C   (ILIN.NE.0)    Run LINOPT
C   (IEXTSIZE.NE.0)USE EXERNAL INITIAL EMITTANCES FOR TRACKING M-C.
C   (IRAD.EQ.1)    INCLUDES ENERGY LOSS IN ELEMENTS FOR GETTING THE CO.
C--------------------------------------------------------------------.
C   (KSEED1,  WID1, SCUT1)  VERTICAL QUAD SHIFTS, WIDTH, GAUSSIAN CUT
C   (KSEED2,  WID2, SCUT2)  HORIZONTAL QUAD SHIFTS, WIDTH, GAUSSIAN CUT
C   (KSEED3,  WID4, SCUT3)  QUAD ROLE, WIDTH, GAUSSIAN CUT
C   (KSEED4,  WID5, SCUT4)  QUAD FRACTIONAL STRENGTH ERROR
C   (KSEED5,  WID3, SCUT5)  DIPOLE ROLE, WIDTH, GAUSSIAN CUT
C   (KSEED6,  WID6, SCUT6)  MONITOR SHIFT, WIDTH, GAUSSIAN CUT
C   (KSEED7,  WID7, SCUT7)  MONITOR SHIFT, WIDTH, GAUSSIAN CUT
C   (KSEED8,  WID8, SCUT8)  QUAD-MONITOR OFFSET FOR KICK MINIMISATION
C   (KSEED9,  WID9, SCUT9)  KILL MONITORS, FRACTION.
C   (KSEED10, WID10 SCUT10)
C--------------------------------------------------------------------
C   (IBEQUIL.NE.0)    TRACK TO GET EQUILIBRIUM PHASE SPACE
C   (NPART1)          NUMBER OF PARTICLES
C   (NTURN1)          NUMBER OF TURNS
C   (NDAMP1)          ARTIFICIAL DAMPING FOR TRACKING
C   (KSEED11)         INITIAL ORBIT DISTRIBUTION. FOR 0 START ON C.O.
C   (KSEED12)         BIG PHOTON SERIES
C--------------------------------------------------------------------
C   (IDEPLIN.NE.0)    TRACK TO MEASURE DEPOLARISATION: LINEARISED SPIN
C   (NPART2)          NUMBER OF PARTICLES
C   (NTURN2)          NUMBER OF TURNS
C   (NDAMP2)          ARTIFICIAL DAMPING FOR TRACKING
C   (KSEED13)         INITIAL S-O DISTRIBUTION. FOR 0 START ON C.O.
C   (KSEED14)         BIG PHOTON SERIES
C--------------------------------------------------------------------
C   (IDEPNON.NE.0)    TRACK TO MEASURE DEPOLARISATION: 3_D SPIN
C   (NPART3)          NUMBER OF PARTICLES
C   (NTURN3)          NUMBER OF TURNS
C   (NDAMP3)          ARTIFICIAL DAMPING FOR TRACKING
C   (KSEED15)         INITIAL S-O DISTRIBUTION. FOR 0 START ON C.O.
C   (KSEED16)         BIG PHOTON SERIES
C   (MODES)           1/0: DEPN TIME for SEPARATE/TOTAL ORBITAL MODES
C    RADFAC        Scale for size of big photons.
C   =================================================================

C
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE "cnlist.for"
      INCLUDE "clatic.for"
      INCLUDE "csol.for"
      INCLUDE "cintpt.for"
      INCLUDE "csex.for"
      INCLUDE "cglist.for"
      INCLUDE "cmags.for"
      INCLUDE "ccentr.for"
      INCLUDE "cmontc.for"
C
C
C
      DIMENSION KGBIN1(1000),KGBIN2(1000)
      DATA     KGBIN1/1000*0/,KGBIN2/1000*0/
      DIMENSION XXTEMP(4)
      DIMENSION KPS(2)
      DIMENSION PSUP(100,4)
      DATA PSUP/400*1.D40/
      CHARACTER *4 TEXT(20)
      CHARACTER STRING *80
      CHARACTER *8 NAM
      REAL*8 RPOS
      CHARACTER *8 NAMEREM
      CHARACTER *1 KOM

      character(100) str
      integer i, ix

      INTEGER END_COUNT ! FIX OF BMAD PRODUCING TWO END ELEMENTS

C===========Initialization for NAG library G05KFF=================
      INTEGER      NNAG
      PARAMETER    (NNAG=20000)
      INTEGER      MSTATE, MSEED
      PARAMETER    (MSTATE=633, MSEED=1)
      INTEGER      GENID, SUBID, LSEED, LSTATE, IFAIL
      DOUBLE PRECISION XNAG(NNAG), ANAG, BNAG
      INTEGER      SEED(MSEED), STATE(MSTATE)
C=================================================================

C=====NAMELIST READS THE INPUT PARAMETERS: ENERGY,PRINT-COMMANDS ETC
       NAMELIST/RUNSPEC/NTY,NINS,IRAD,ISPIN,IPTCO,IDISP,IPTNML,ISIG,
     + NSTEP,E00,DE0,CIRCUM,IBMBM,NITER,IHARM,KICK,SCLKIK,ICLBP,ITWP,
     + IPTBET,ISECT,ECAV,SFREQ,DFREQ,NFREQ,IPTD,MRDISP,IVROT,ECHROM,
     + ICORR,IKMIN,IEDGE,ISEXTU,ILIN,NDIFFDUMP,NSKIP,RADFAC,IEXTSIZE,
     + KSEED1,  WID1, SCUT1,
     + KSEED2,  WID2, SCUT2,
     + KSEED3,  WID3, SCUT3,
     + KSEED4,  WID4, SCUT4,
     + KSEED5,  WID5, SCUT5,
     + KSEED6,  WID6, SCUT6,
     + KSEED7,  WID7, SCUT7,
     + KSEED8,  WID8, SCUT8,
     + KSEED9,  WID9, SCUT9,
     + KSEED10, WID10,SCUT10,
     + KSEED11, WID11,SCUT11,
     + KSEED12, WID12,SCUT12,
     + KSEED13, WID13,SCUT13,
     + KSEED14, WID14,SCUT14,
     + KSEED15, WID15,SCUT15,
     + KSEED16, WID16,SCUT16,
     + KSEED17, WID17,SCUT17,
     + KSEED18, WID18,SCUT18,
     + KSEED19, WID19,SCUT19,
     + KSEED20, WID20,SCUT20,
     + IBEQUIL, NPART1, NTURN1, NDAMP1,
     + IDEPLIN, NPART2, NTURN2, NDAMP2,
     + IDEPNON, NPART3, NTURN3, NDAMP3, MODES

C=====Sept. 2003: replace the routine RANNOR for Gaussians with a NAG routine.
C     Comment out RANNOR code.
C      DIMENSION ISEED(20)
C      DATA ISEED/    12345, 51233, 45123, 34515, 23451,
C     +              112345,151233,145123,134515,123451,
C     +              212345,251233,245123,234515,223451,
C     +              312345,351233,345123,334515,323451/
C
C
C
C
C
C
C
C
C
C
C
C
C=====SET DEFAULT VALUES
      PI=3.1415926535897932D0
      NINS=1
      MODES=0
      IBMBM=0
      IRAD=0
      ISPIN=1
      IPTCO=1
      NITER=3
      IDISP=0
      IPTBET=1
      IPTNML=0
      IHARM=3840
      KICK=0
      SCLKIK=1.D0
      ICORR=0
      IKMIN=0
      ECHROM=0.D0
      ICLBP=1
      ITWP=0
      ISIG=0
      NSTEP=10
      ISECT=0
      DE0=0.008D+0
      SIGYBB=0.00006D+0
      ECAV=0.D+0
      SFREQ=0.D+0
      DFREQ=0.D+0
      NFREQ=0
      IROT=0
      IOPEN=0
      IEDGE=0
      ISEXTU=0
      ILIN=1
      NDIFFDUMP=1
      NSKIP=1
      RADFAC=1.D0
      IFORM=0
      BETAX0=0.D+0
      BETAY0=0.D+0
      NSEX=0
      KSEED1 =0;KSEED2 =0;KSEED3 =0;KSEED4 =0;KSEED5 =0
      KSEED6 =0;KSEED7 =0;KSEED8 =0;KSEED9 =0;KSEED10=0
      KSEED11=0;KSEED12=0;KSEED13=0;KSEED14=0;KSEED15=0
      KSEED16=0;KSEED17=0;KSEED18=0;KSEED19=0;KSEED20=0
      SCUT1 =0;SCUT2 =0;SCUT3 =0;SCUT4 =0;SCUT5 =0
      SCUT6 =0;SCUT7 =0;SCUT8 =0;SCUT9 =0;SCUT10=0
      SCUT11=0;SCUT12=0;SCUT13=0;SCUT14=0;SCUT15=0
      SCUT16=0;SCUT17=0;SCUT18=0;SCUT19=0;SCUT20=0
      WID1 =0.D0;WID2 =0.D0;WID3 =0.D0;WID4 =0.D0;WID5 =0.D0
      WID6 =0.D0;WID7 =0.D0;WID8 =0.D0;WID9 =0.D0;WID10=0.D0
      WID11=0.D0;WID12=0.D0;WID13=0.D0;WID14=0.D0;WID15=0.D0
      WID16=0.D0;WID17=0.D0;WID18=0.D0;WID19=0.D0;WID20=0.D0
      IBEQUIL=0; NPART1=1000; NTURN1=100000;NDAMP1=1
      IDEPLIN=0; NPART2=1000; NTURN2=100000;NDAMP2=1
      IDEPNON=0; NPART3=1000; NTURN3=100000;NDAMP2=1
      KG1 =0;KG2 =0;KG3 =0;KG4 =0;KG5 =0
      KG6 =0;KG7 =0;KG8 =0;KG9 =0;KG10=0
      KG11=0;KG12=0;KG13=0;KG14=0;KG15=0
      KG16=0;KG17=0;KG18=0;KG19=0;KG20=0
      IBB = 0
C
      CENPOS = -1000.D0 !Element centres. They remain  large and -ve if not set later.
      RADSCALE = 0.D0   !Radiation scale for the kicks at the starts of sections.
C
C
C
      CALL DATET
C
C=====READ & WRITE THE NAMELISTS,MANIPULATE PARAMETERS.
      READ (5,RUNSPEC)
      CALL IPARAM

      IF(IBEQUIL+IDEPLIN+IDEPNON.GT.1)THEN
      WRITE(53,103)
  103 FORMAT(/,'  ')
      WRITE(53,'(3A)')' STOP: attempt to use more than 1 one kind of',
     +                 '  Monte-Carlo algorithm.'
      STOP
      ENDIF

      IF((IDEPNON.NE.0. AND.NTURN3.LE.NSKIP).OR.
     +    (IDEPLIN.NE.0. AND.NTURN2.LE.NSKIP))THEN
      WRITE(53,103)
      WRITE(53,'(3A)')' STOP: NDAMP3 is too small for an estimate of',
     + ' the depolarisation rate.'
      STOP
      ENDIF

C
C
C======Set up independent Gaussian random number series within the cuts.
C      G9 uses a uniform distribution to get on/off flags for the monitors.
C      It is set to 1.D0 as a default, i.e. by default all monitors are on.
C      Test the Gaussians by plotting a G1 distribution.
      G1=0.D0;G2=0.D0;G3=0.D0;G4=0.D0;G5=0.D0;
      G6=0.D0;G7=0.D0;G8=0.D0;G9=0.D0;G10=0.D0
      IF(KICK.GT.0)THEN

      GENID   = 1
      SUBID   = 1
      LSTATE  = MSTATE
      LSEED   = MSEED
      IFAIL   = 1
C------Vertical quad. offsets.
      IF(KSEED1.NE.0)THEN
      SEED(1) = KSEED1
C      CALL G05CBF(KSEED1)   G05CBF is replaced by G05KFF
      CALL G05KFF(GENID,SUBID,SEED,LSEED,STATE,LSTATE,IFAIL)
      IF(IFAIL .NE. 0)THEN
      WRITE(53,'(I10)') IFAIL
      ENDIF

      NG1   = 0
      NG    = 0
      GBIN  = 0.1D0
      SUMG  = 0.D0
      SUMG2 = 0.D0
C 1001 GAUS1 = G05DDF(0.D0,WID1)   G05DDF is replaced by G05SKF
      CALL G05SKF(NNAG,0.D0,WID1,STATE,XNAG,IFAIL)
      IF(IFAIL .NE. 0)THEN
        WRITE(53,'(I10)') IFAIL
      ENDIF

      DO JNAG = 1, NNAG
        SUMG2 = SUMG2 + XNAG(JNAG)**2
        SUMG  = SUMG  + XNAG(JNAG)
        NG    = NG + 1
        NGAUS1 = XNAG(JNAG)/GBIN + 50
        IF(NGAUS1.GT.0)KGBIN1(NGAUS1) = KGBIN1(NGAUS1) + 1
C        WRITE(31,'(1X,I10,1X,F10.5)')JNAG,XNAG(JNAG)
      ENDDO

      DO JNAG = 1, NNAG
        IF(DABS(XNAG(JNAG)).LT.SCUT1*WID1)THEN
        NG1 = NG1 + 1
        G1(NG1) = XNAG(JNAG)*SCLKIK
        NGAUS2 = G1(NG1)/GBIN + 50
        IF(NGAUS2.GT.0)KGBIN2(NGAUS2) = KGBIN2(NGAUS2) + 1
C        WRITE(32,'(1X,I10,1X,F10.5)')NG1,G1(NG1)
        ENDIF
      ENDDO
      ENDIF
C======Write out statistics on the Gaussians.
      IF(NG.GT.0)THEN
      SUMG2 = SUMG2/NG
      SUMG  = SUMG/NG
      SUMG2 = DSQRT(SUMG2)
      WRITE(53,'(1X,A,I10,2(1X,F15.6))')
     + 'Check the Gaussian: calls, mean, r.m.s. ',NG,SUMG,SUMG2
      ENDIF
      DO JG = 1,1000
        WRITE(56,'(1X,I10,1X,I10,1X,I10)')JG,KGBIN1(JG),KGBIN2(JG)
      ENDDO

C      IF(1 .EQ. 1)STOP

C------Horizontal quad. offsets.
      IF(KSEED2.NE.0)THEN
      SEED(1) = KSEED2
C      CALL G05CBF(KSEED2)   G05CBF is replaced by G05KFF
      CALL G05KFF(GENID,SUBID,SEED,LSEED,STATE,LSTATE,IFAIL)
      NG2 = 0
C 1002 GAUS1 = G05DDF(0.D0,WID2)   G05DDF is replaced by G05SKF
      CALL G05SKF(NNAG,0.D0,WID2,STATE,XNAG,IFAIL)

      DO 7002 JNAG = 1, NNAG
      IF(DABS(XNAG(JNAG)).LT.SCUT2*WID2)THEN
      NG2 = NG2 + 1
      G2(NG2) = XNAG(JNAG)*SCLKIK
      ENDIF
7002  CONTINUE
      ENDIF
C------Quadrupole role.
      IF(KSEED3.NE.0)THEN
      SEED(1) = KSEED3
C      CALL G05CBF(KSEED3)   G05CBF is replaced by G05KFF
      CALL G05KFF(GENID,SUBID,SEED,LSEED,STATE,LSTATE,IFAIL)
      NG3 = 0
C 1003 GAUS1 = G05DDF(0.D0,WID3)   G05DDF is replaced by G05SKF
      CALL G05SKF(NNAG,0.D0,WID3,STATE,XNAG,IFAIL)

      DO 7003 JNAG = 1, NNAG
      IF(DABS(XNAG(JNAG)).LT.SCUT3*WID3)THEN
      NG3 = NG3 + 1
      G3(NG3) = XNAG(JNAG)*SCLKIK
      ENDIF
7003  CONTINUE
      ENDIF
C------Quadrupole fractional strength error.
      IF(KSEED4.NE.0)THEN
      SEED(1) = KSEED4
C      CALL G05CBF(KSEED4)   G05CBF is replaced by G05KFF
      CALL G05KFF(GENID,SUBID,SEED,LSEED,STATE,LSTATE,IFAIL)
      NG4 = 0
C 1004 GAUS1 = G05DDF(0.D0,WID4)   G05DDF is replaced by G05SKF
      CALL G05SKF(NNAG,0.D0,WID4,STATE,XNAG,IFAIL)

      DO 7004 JNAG = 1, NNAG
      IF(DABS(XNAG(JNAG)).LT.SCUT4*WID4)THEN
      NG4 = NG4 + 1
      G4(NG4) = XNAG(JNAG)*SCLKIK
      ENDIF
7004  CONTINUE
      ENDIF
C------Dipole role.
      IF(KSEED5.NE.0)THEN
      SEED(1) = KSEED5
C      CALL G05CBF(KSEED5)   G05CBF is replaced by G05KFF
      CALL G05KFF(GENID,SUBID,SEED,LSEED,STATE,LSTATE,IFAIL)
      NG5 = 0
C 1005 GAUS1 = G05DDF(0.D0,WID5)   G05DDF is replaced by G05SKF
      CALL G05SKF(NNAG,0.D0,WID5,STATE,XNAG,IFAIL)

      DO 7005 JNAG = 1, NNAG
      IF(DABS(XNAG(JNAG)).LT.SCUT5*WID5)THEN
      NG5 = NG5 + 1
      G5(NG5) = XNAG(JNAG)*SCLKIK
      ENDIF
7005  CONTINUE
      ENDIF
C------Monitor offset.
      IF(KSEED6.NE.0)THEN
      SEED(1) = KSEED6
C      CALL G05CBF(KSEED6)   G05CBF is replaced by G05KFF
      CALL G05KFF(GENID,SUBID,SEED,LSEED,STATE,LSTATE,IFAIL)
      NG6 = 0
C 1006 GAUS1 = G05DDF(0.D0,WID6)   G05DDF is replaced by G05SKF
      CALL G05SKF(NNAG,0.D0,WID6,STATE,XNAG,IFAIL)

      DO 7006 JNAG = 1, NNAG
      IF(DABS(XNAG(JNAG)).LT.SCUT6*WID6)THEN
      NG6 = NG6 + 1
      G6(NG6) = XNAG(JNAG)*SCLKIK
      ENDIF
7006  CONTINUE
      ENDIF
C------Monitor scale.
      IF(KSEED7.NE.0)THEN
      SEED(1) = KSEED7
C      CALL G05CBF(KSEED7)   G05CBF is replaced by G05KFF
      CALL G05KFF(GENID,SUBID,SEED,LSEED,STATE,LSTATE,IFAIL)
      NG7 = 0
C 1007 GAUS1 = G05DDF(0.D0,WID7)   G05DDF is replaced by G05SKF
      CALL G05SKF(NNAG,0.D0,WID7,STATE,XNAG,IFAIL)

      DO 7007 JNAG = 1, NNAG
      IF(DABS(XNAG(JNAG)).LT.SCUT7*WID7)THEN
      NG7 = NG7 + 1
      G7(NG7) = XNAG(JNAG)*SCLKIK
      ENDIF
7007  CONTINUE
      ENDIF
C------Quad-monitor offset for kick minimisation.
      IF(KSEED8.NE.0)THEN
      SEED(1) = KSEED8
C      CALL G05CBF(KSEED8)   G05CBF is replaced by G05KFF
      CALL G05KFF(GENID,SUBID,SEED,LSEED,STATE,LSTATE,IFAIL)
      NG8 = 0
C 1008 GAUS1 = G05DDF(0.D0,WID8)   G05DDF is replaced by G05SKF
      CALL G05SKF(NNAG,0.D0,WID8,STATE,XNAG,IFAIL)

      DO 7008 JNAG = 1, NNAG
      IF(DABS(XNAG(JNAG)).LT.SCUT8*WID8)THEN
      NG8 = NG8 + 1
      G8(NG8) = XNAG(JNAG)*SCLKIK
      ENDIF
7008  CONTINUE
      ENDIF
C------Monitor on/off: WID8 is the fraction of dead monitors.
      IF(KSEED9.NE.0)THEN
      SEED(1) = KSEED9
C      CALL G05CBF(KSEED9)   G05CBF is replaced by G05KFF
      CALL G05KFF(GENID,SUBID,SEED,LSEED,STATE,LSTATE,IFAIL)
      NG9 = 0
      ANAG = 0.D0
      BNAG = 1.D0
C 1009 CALL G05FAF(0.D0,1.D0,1,ONOFF)   G05FAF is replaced by G05SQF
      CALL G05SQF(NNAG,ANAG,BNAG,STATE,XNAG,IFAIL)

      DO 7009 JNAG = 1, NNAG
      NG9 = NG9 + 1
      G9(NG9) = 1.D0
      IF(XNAG(JNAG).GE.(1.D0 - WID9))G9(NG9) = 0.D0
7009  CONTINUE
      ENDIF
      ENDIF
C======Avoid using the same seeds for the monitor shift and the quad-monitor offset
      IF((KSEED6.EQ.KSEED8).AND.KSEED6.NE.0)THEN
      WRITE(53,'(A)')' M-C seeds 8 and 9 are the same, so STOP'
      STOP
      ENDIF

C=====SET UP UNIT MATRICES
      DO I=1,20000
      DO J=1,7
      DO K=1,7
      TMAT(J,K,I)=0.D+0
      IF(J.EQ.K)TMAT(J,K,I)=1.D+0
      ENDDO
      ENDDO
      ENDDO
C
      NUNT = 1
C
C
C
C
      IF(NSTEP.LE.0)NSTEP=1
      NKICK=(KICK/10)*10
      IKICK=KICK-NKICK
C
C
C
C
C
C=====READ THE DATA FILE HEADING
C     READ (52,901)TEXT
C     WRITE(53,901)TEXT
C     READ (52,901)TEXT
C     WRITE(53,901)TEXT
C     READ (52,901)TEXT
C     WRITE(53,901)TEXT
C
C=====CHECK FORMAT AND RANDOM NUMBER GENERATOR SEEDS.
C      CALL RDMOUT(JJSEED)
C     READ(52,801)IFORM
C     WRITE(53,801)IFORM
C     READ(52,805)KSEED                ! Old random generator seed.
C     WRITE(53,805)KSEED
C     READ(52,805)NINSOLD              ! Ignore: history.
C     WRITE(53,805)NINSOLD

C      JSEED=ISEED(KSEED)
C      WRITE(53,802)IFORM,JSEED,JJSEED
C     IF(IFORM.EQ.4.OR.IFORM.EQ.5.OR.IFORM.EQ.6)GO TO 804
C     WRITE(53,803)
C     STOP


  804 CONTINUE
C      CALL G05CBF(KSEED1)   G05CBF is replaced by G05KFF
C      SEED(1) = KSEED1
C      CALL G05KFF(GENID,SUBID,SEED,LSEED,STATE,LSTATE,IFAIL)

C     SCL= 0.001D+0
C    SEP=-0.001001D+0
     SEP=-0.0001001D+0

C     IF(IFORM.EQ.5.OR.IFORM.EQ.4)SCL=0.0001D+0
C     IF(IFORM.EQ.5.OR.CFORM.EQ.4)SEP=-0.001101D+0
C
C
C
      IF(IVROT.NE.0)CALL CHROTN(E00,0,1)
C      CALL RDMIN(JSEED)   ! Kill March 09 when dumping CERNLIB
C
C
C
C
      NTY=1
C=====LOOP TO READ THE MAGNET PARAMETERS & CHECK FOR COMMENTS.
      WRITE(53,952)
      END_COUNT = 0
      DO WHILE(.true.)
        READ(52,'(A)')STR
C       write(*,*) STR
        IF(STR(:1).NE.' ')THEN
          WRITE(53,'(1X,A,A)')STR, '  Comment'
          CYCLE
        ENDIF
        ix = 0
        call string_trim(str(ix+1:), str, ix)
        if (ix == 0) CYCLE
        read (str(1:ix), *) ID(NTY)
        call string_trim(str(ix+1:), str, ix)
        NAME(NTY) = str(1:ix)
        IF(NAME(NTY).EQ.'END'.and.END_COUNT==1) EXIT
        IF(NAME(NTY).EQ.'END'.and.END_COUNT==0) THEN
          END_COUNT = END_COUNT+1
          CYCLE
        ENDIF
        call string_trim(str(ix+1:), str, ix)
        read (str(1:ix), *) XX(NTY)
        call string_trim(str(ix+1:), str, ix)
        read (str(1:ix), *) X2(NTY)
        call string_trim(str(ix+1:), str, ix)
        read (str(1:ix), *) YY(NTY)
        call string_trim(str(ix+1:), str, ix)
        read (str(1:ix), *) NUNT(NTY)
        call string_trim(str(ix+1:), str, ix)
        read (str(1:ix), *) TWIST(NTY)
        call string_trim(str(ix+1:), str, ix)
        read (str(1:ix), *) NSOL(NTY)
        SNAME(NTY)= '--------'

C======Switch sextupoles on/off.
        IF(ID(NTY).EQ.8.AND.ISEXTU.EQ.0)XX(NTY)=0.D0

C
C=====MODIFY LENGTHS ACCORDING TO MAGNET SUBDIVISION REQUESTED.
C     YY(NTY)=IDINT(YY(NTY)*1000.D0+0.5)/1000.D0
C=====14/10/85 KILL NEXT LINE. WHAT WAS IT FOR?
C     YY(NTY)=IDINT(YY(NTY)*1000.D0+0.5*0.D0)/1000.D0

C=====If the slicing is patched into the optic files, then XX and X2 must be fixed too.
        IF(NUNT(NTY).GT.1)XX(NTY)=XX(NTY)/NUNT(NTY)     ! 30/6/09 testing with NUNT = 100.
        IF(NUNT(NTY).GT.1)X2(NTY)=X2(NTY)/NUNT(NTY)     ! 30/6/09 testing with NUNT = 100.

        IF(NUNT(NTY).GT.1)YY(NTY)=YY(NTY)/NUNT(NTY)
C=====INPUT SOLENOID STRENGTH REFERS TO THE WHOLE SOLENOID--SO DIVIDE
C=====IT UP.NOTE CONTRAST WITH OTHER CASES WHERE THE DIVISION WAS DONE
C=====AS THE PETROS FILE WAS READ AND CONVERTED.
        IF(NUNT(NTY).GT.1.AND.ID(NTY).EQ.10)XX(NTY)=XX(NTY)/NUNT(NTY)
C
C=====SET UP 6X6 TRANSFER MATRICES.
        J=ID(NTY)

        IF(J.eq.4.and.X2(NTY).lt.1.D-16) THEN
          ID(NTY) = 3
          J = 3
        ENDIF
        SELECT CASE(J)
          CASE(99)
C           DO NOTHING
C=======D
C=====DRIFT: STRENGTH=LENGTH: GET IT (EXCEPTIONALLY) FROM YY. ALSO STORE IN XX
          CASE(1)
            TMAT(1,2,NTY)= YY(NTY)
            TMAT(3,4,NTY)= YY(NTY)
            XX(NTY)=YY(NTY)
C=======B
          CASE(2)
C=====HORIZONTAL BENDING + EDGE FIELD.
C=====TO SET UP EDGE FIELD MATRIX USE THIN LENS QUAD. SET LENGTH TO SMALL
C=====TO MAINTAIN ORBIT LENGTH & STORE XX AS ZERO TO AVOID RESULTING
C=====LARGE RADIATION EFFECTS.SLOT INTO THE TYPE LIST BEHIND THE DIPOLE
C=====THIS IS A MESS. IT NEEDS FIXING.
         IF(IVROT.NE.0)CALL CHROTN(E00,NTY,2)
C        YY(NTY)=YY(NTY)*XX(NTY)/SIN(XX(NTY))
         !  modified on Sep 29, 2020 as David Sagan pointed out
         CALL DIPH(XX(NTY),YY(NTY),TMAT(1,1,NTY))

         IF(NAME(NTY)(1:1).EQ.'W')THEN
         TMAT(5:6,1:6,NTY) = 0.D0
         TMAT(1:6,5:6,NTY) = 0.D0
         TMAT(5,5,NTY)     = 1.D0
         TMAT(6,6,NTY)     = 1.D0
         ID(NTY) = 3
         ENDIF
C        NTY=NTY+1
C        CALL HEDGE(XX(NTY-1),YY(NTY-1),TMAT(1,1,NTY),NUNT(NTY-1),
C    +                               XX(NTY),X2(NTY),NSOL(NTY-1))
C        NAME(NTY)='EDGE'
C        SNAME(NTY)= '--------'
C        NUNT(NTY)=1
C        ID(NTY)=3
C        XX(NTY)=0.D0                        !Give zero strength to handle radiation.
C        X2(NTY)=0.D0
C        YY(NTY)=1.D-15                      !And give edge fields some length too.
C        IF(IEDGE.EQ.0)THEN            !Turn off edge focussing.
C          YY(NTY)=0.D0
C          XX(NTY)=0.D0
C          ID(NTY)=1
C        ENDIF
C=======Q
          CASE(3)
C=====NORMAL QUAD WITH SPECIAL TRICK FOR B-B SUBSTITUTE
           IF(NAME(NTY)(1:2).NE.'BB')THEN
           CALL THIQUAD(XX(NTY),YY(NTY),IKICK,NTY,TMAT(1,1,NTY),0.D0)
     C
           ELSE
           TMAT(2,1,NTY)= -XX(NTY)             !Thin lens b-b. Use XX and X2 to
           TMAT(4,3,NTY)= -X2(NTY)             !handle elliptical oncoming beams.
           XX(NTY) = XX(NTY)*0.000000D0        !Then kill XX to avoid problems with
           X2(NTY) = X2(NTY)*0.000000D0        !radiation. Kill X2 too for good measure.
           ENDIF                              !This facility is essential if worries
                                              !arise about sectioning: with this trick
                                               !the b-b is absorbed into an existing
                                               !section.
C=======QK
          CASE(4)
C=====SKEW QUAD
            CALL SKEWQUAD(XX(NTY),YY(NTY),IKICK,
     +        NTY,TMAT(1,1,NTY),X2(NTY))
C=======RF
          CASE(5)
C=====CAVITY
            TMAT(6,5,NTY)= XX(NTY)
            ! This is overwritten with a sensible value in cavity.f.
C=======HK
          CASE(6)
C=====HOR. KICKER-----USE THIN LENS TYPE MATRIX ELEMENTS.
            TMAT(2,6,NTY)= XX(NTY)
            TMAT(5,1,NTY)=-XX(NTY)
            TMAT(2,7,NTY)=-XX(NTY)
C=======VK
          CASE(7)
C=====VERT KICKER-----USE THIN LENS TYPE MATRIX ELEMENTS:M&R SIGNS
            TMAT(4,6,NTY)=+XX(NTY)
            TMAT(5,3,NTY)=-XX(NTY)
            TMAT(4,7,NTY)=-XX(NTY)
C=======BY
          CASE(9)
C=====VERTICAL BEND
            IF(IVROT.NE.0)CALL CHROTN(E00,NTY,2)
            CALL DIPV(XX(NTY),YY(NTY),TMAT(1,1,NTY))
C=======SOL
          CASE(10)
C=====USE G.RIPKEN'S 7X7 MATRICES FOR SHIFTED SOLENOIDS.
C=====CHECK PREVIOUS ELEMENT TO SEE IF IT IS A SPLIT SOLENOID.
            ISPLIT=0
            IF(ID(NTY-1).EQ.10.AND.XX(NTY-1)*XX(NTY).GT.0.D0)ISPLIT=1
            CALL SOL77(XX(NTY),YY(NTY),TMAT(1,1,NTY),IKICK,ISPLIT,
     +        TWIST(NTY),NSOL(NTY))
C=======RQ
          CASE(11)
C=====THIN ROTATED QUADS.
            CALL ROQUAD(NTY)
C=======FBH,FBV
          CASE(12)
            WRITE(53,953)
          CASE(13)
            WRITE(53,953)
C=====UNCOUPLED SECTION OF ROTATOR
C 116 CALL ROTATR(NTY)
C=======ROT
          CASE(14)
            CALL SNAKE(NTY)
            IROT=1
C=======SX
          CASE(8)
C=====SEXTUPOLE:SET UP AS A DRIFT TO START WITH.
C=====FOR 2-FAMILY CORRECTION,LABEL ACCORDING TO HERA LENGTHS.
            TMAT(1,2,NTY)= YY(NTY)
            TMAT(3,4,NTY)= YY(NTY)
            IF(YY(NTY).GT.0.61.AND.YY(NTY).LT.0.63)ISEXFM(NTY)=1
            IF(YY(NTY).GT.0.27.AND.YY(NTY).LT.0.29)ISEXFM(NTY)=2
C=======CFH
          CASE(15)
C=====HORIZONTAL COMBINED FUNCTION DIPOLE
            CALL CFDIPH(XX(NTY),X2(NTY),YY(NTY),TMAT(1,1,NTY),0.D0)
C=====Put in an edge field just as for hor dipoles.
C           NTY=NTY+1
C           XIN=XX(NTY-1)
C           CALL HEDGE((XIN),YY(NTY-1),TMAT(1,1,NTY),NUNT(NTY-1),
C    +                               XX(NTY),X2(NTY),NSOL(NTY-1))
C           NAME(NTY)='EDGE'
C           SNAME(NTY)= '--------'
C           NUNT(NTY)=1
C           ID(NTY)=3
C           XX(NTY)=0.D0                     !Give zero strength to handle radiation.
C           X2(NTY)=0.D0
C           YY(NTY)=1.D-10                   !And give edge fields some length too.
C           IF(IEDGE.EQ.0)THEN        ! Turn off edge focussing.
C           YY(NTY)=0.D0
C           XX(NTY)=0.D0
C           ID(NTY)=1
C           ENDIF
C=======CFV
          CASE(16)
C=====VERTICAL COMBINED FUNCTION DIPOLE
            CALL CFDIPV(XX(NTY),X2(NTY),YY(NTY),TMAT(1,1,NTY),0.D0)
C=======BB
          CASE(17)
C=====BEAM-BEAM: need width, height and bunch current. Also need
C     circumference and energy.
            CALL SETBB(XX(NTY),X2(NTY),YY(NTY),TMAT(1,1,NTY))
C=======EDGE FIELDS for SPRINT PROGRAMS
          CASE(97)
C=====Edge focussing from input file. Redefine the type to 3.
C           XX(NTY) = XX(NTY)
C           XX(NTY) = 0.d0
            TMAT(2,1,NTY) =  XX(NTY)  * 1.0D0     !Fanglei modified on May 28, 2020, horizontal is defocusing, vertical is focusing
            TMAT(4,3,NTY) = -XX(NTY)  * 1.0D0     !SY's book, page 66, exercise 2.2.2

            ID(NTY)=3                      !Change the type number from 97 to 3.
            XX(NTY)=0.D0                   !Give zero strength.
            YY(NTY)=1.D-16                 !Don't give edge fields zero length.
            NSOL(NTY) = 999                ! MARK AS EDGE FIELD
C           YY(NTY)=5.D-2                  !Give edge fields some length too.
C=====SET UP A FIRST DRIFT OF ZERO LENGTH----THIS IS NOT THE I.P.!
C=====THE I.P. WILL BE THE 3RD ELEMENT!---AFTER THIS DUMMY ELEMENT & THE
C=====FOLLOWING ZERO LENGTH DRIFT.
C=======EDGE FIELDS for SPRINT PROGRAMS
          CASE(98)
C=====Edge focussing from input file. Redefine the type to 3.
C           XX(NTY) = XX(NTY)
C           XX(NTY) = 0.d0
            TMAT(2,1,NTY) = -XX(NTY)  * 1.0D0
            TMAT(4,3,NTY) =  XX(NTY)  * 1.0D0

            ID(NTY)=3                      !Change the type number from 97 to 3.
            XX(NTY)=0.D0                   !Give zero strength.
            YY(NTY)=1.D-16                 !Don't give edge fields zero length.
            NSOL(NTY) = 998                ! MARK AS EDGE FIELD
C           YY(NTY)=5.D-2                  !Give edge fields some length too.
C=====SET UP A FIRST DRIFT OF ZERO LENGTH----THIS IS NOT THE I.P.!
C=====THE I.P. WILL BE THE 3RD ELEMENT!---AFTER THIS DUMMY ELEMENT & THE
C=====FOLLOWING ZERO LENGTH DRIFT.
          CASE DEFAULT
            WRITE(*,*) 'UNDEFINED TYPE J=',J
            WRITE(53,96) J
            STOP
        END SELECT

C
C
        NTY=NTY+1
        IF(NTY.GE.20000) THEN
          WRITE(53,94)
          STOP
        ENDIF
      ENDDO

      ND1=NTY
      NAME(NTY)='DRIFT'
      SNAME(NTY)= '--------'
      ID(NTY)=1
      XX(NTY)=0.D0
      YY(NTY)=0.D0
      NT1=NTY-1

      NELEM=1
      ITYPE(1)=NTY
      IREAD=0
      TOTDIP=0.
C
C
C
C
C=====READ THE STRUCTURE OF THE RING INTO STORAGE.
C=====CAVITIES,KICKERS AND B-B GET ZERO LENGTH (FOR CAVITIES YY=0 ON INPUT)
C=====EVERYTHING ELSE INCLUDING SEXTUPOLES GETS ITS ACTUAL LENGTH.
C=====Initialise the memory for quadrupole and dipole strengths to apply to various error
C     thin lenses.
      QREM  =  0.D0
      QYREM = 10.D0
      DREM  =  0.D0
      DYREM = 10.D0

      KK = 0
      WRITE(53,93)
      DO WHILE(.true.)
        IF(KK.eq.1) EXIT
        READ(52,'(A)')STR
        IF(STR(:2).EQ.'--') CYCLE
        ix = 0
        DO WHILE(.true.)
          call string_trim(str(ix+1:), str, ix)
          if (ix == 0) EXIT
          ! Read name
          NAM = str(1:ix)
C         IF(NAM.eq.'IP') write(*,*)' spotted ip'
          call string_trim(str(ix+1:), str, ix)
          ! Read posiztion
          read (str(1:ix), *) RPOS
          IF(NAM.EQ.'END') THEN
            KK=1
            EXIT
          ENDIF
          IF(NAM.EQ.'OPEN')THEN
            IOPEN=1
            KK=1
            EXIT
          ENDIF
          WRITE(53,9555) NAM,RPOS
          IF(IREAD.EQ.0)TOTL=RPOS
          IREAD=1
          IREM=0
C=====    SEARCH FOR NAME IN TYPE LIST.
          DO INN=1,NT1
            IN = INN
            IF(NAM.EQ.NAME(IN))THEN
            ! FOUND
              IREM=ID(IN)
              EXIT
            ENDIF
          ENDDO
          IF(IREM .eq. 0) THEN
            ! NOT FOUND
            WRITE(53,91)NAM
            STOP
          ENDIF
          NU=NUNT(IN)
          NUS=NU
          IF(ID(IN).EQ.2.OR.ID(IN).EQ.15)TOTDIP=TOTDIP+YY(IN)*NUNT(IN)
          PLT=RPOS
          PL=PLT-YY(IN)*NU*0.5D0-TOTL

C=====    SOL.AND ERSATZ CASE WHERE ENTRANCE COORDS(NOT CENTRES) ARE GIVEN.
C=====    =Can comment this out for files from MAD via SPRINT. But it doesn't work!????
          IF(IREM.EQ.10.OR.IREM.EQ.14)PL=PLT-TOTL
C
C=====    HOR. & VERT. KICKERS AND B-B WHICH ARE GIVEN ZERO LENGTH IN THE LATTICE.
          IF(IREM.EQ. 6.OR.IREM.EQ. 7.OR.IREM.EQ.17.OR.IREM.EQ.5)
     +        PL=PLT-TOTL
          IF(IREM.EQ.17)WRITE(*,'(A,E20.6)')' PL  ', PL
C
C=====    =If a quad has been reached, save the strength for use with a VQ,HQ,RQ,CQ later.
C          Need care to avoid picking up the length of a CQ
C          Note that the lengths of RQ and CQ are set to zero in the input file
C          to allow the lattice to be constructed. So use this to distinguish.
          IF(IREM.EQ.3.AND.YY(IN).GT.0.D0)THEN
          QREM  = XX(IN)   *2.D0    ! *2 to get the strength of the full quad.
          QYREM = YY(IN)   *2.D0    ! *2 to get the   length of the full quad.
          INREM = IN
          ENDIF
C
C=====    =If a CF has been reached, save the strength for use with a VQ,HQ,RQ,CQ later.
C          Need care to avoid picking up the length of a CQ
C          Note that the lengths of RQ and CQ are set to zero in the input file
C          to allow the lattice to be constructed. So use this to distinguish.
          IF(IREM.EQ.15.AND.YY(IN).GT.0.D0)THEN
          QREM  = X2(IN)   *2.D0    ! *2 to get the strength of the full CF.
          QYREM = YY(IN)   *2.D0    ! *2 to get the   length of the full CF.
          ENDIF
C
C=====    =If a dipole or CF has been reached, save the strength for use with a VD.
C          By construction, any VD is in the middle of a magnet of finite length
          IF(IREM.EQ.2.OR.IREM.EQ.15)THEN
          DREM  = XX(IN)   *2.D0    ! *2 to get the strength of the full dipole/CF.
          DYREM = YY(IN)   *2.D0    ! *2 to get the   length of the full dipole/CF.
          NAMEREM = NAME(IN)
          ENDIF
C
C=====    =If a skew quad has been reached, save the strength for use with a VQ,HQ,RQ,CQ later.
C          Need care to avoid picking up the length of a CQ
C          Note that the lengths of RQ and CQ are set to zero in the input file
C          to allow the lattice to be constructed. So use this to distinguish.
          IF(IREM.EQ.4.AND.YY(IN).GT.0.D0)THEN
          QREM  = XX(IN)   *2.D0    ! *2 to get the strength of the full skew quad.
          QYREM = YY(IN)   *2.D0    ! *2 to get the   length of the full skew quad.
          ENDIF
C
C         If there is no drift, and allowing some initial zero length drifts to be found,
C         then don't introduce any more zero length drifts!
C         One can't easily check that superimposed elements or elements right next to others
C         are thin lenses because thin lenses might have been given a finite length to handle
C         radiation. The information ``thin lens'' is in the matrices.
C          WRITE(53,'(A,4I10)')' IPOS etc ',NELEM,IPOSOLD,IPOSNEW
C          IF(NELEM.GT.10.AND.IPOSNEW.EQ.IPOSOLD)GO TO 100
C          PL is in metres. So this cut of 0.00011D0 is very tight (0.11mm)and might
C          need relaxing to take care of rounding errors so that drifts are not excluded.
          IF(NELEM.GT.10.AND.PL.LT.0.00011D0)GO TO 100

C
C=====CONSTRUCT THE DRIFT SPACE UP TO THE NEXT THICK LENS.
C     Check if the drift has already been constructed.
   11     DO 7 III=ND1,NTY
          I = III
          IF(DABS(PL*1000.D0-XX(III)*1000.D0).LT.0.1000001D0)GO TO 10
    7     CONTINUE
          NTY=NTY+1
          IF(NTY.GT.20000)THEN
          WRITE(53,94)
          STOP
          ENDIF

          SNAME(NTY)=NAME(IN)
          NAME(NTY)= 'DRIFT'
          ID(NTY)=1
C====     TAKE CARE OF DIGITISATION ERRORS ON INPUT FILE--LIMIT GIVEN BY SEP
          IF(PL.GE.SEP.AND.PL.LT.0.0D0)PL=0.D0
          XX(NTY)=PL
          YY(NTY)=PL
          TMAT(1,2,NTY)= PL
          TMAT(3,4,NTY)= PL
          I=NTY

   10     CONTINUE
          NELEM=NELEM+1
C          LL=1
C          WRITE(53,908)LL,TOTL,NELEM
          IF(NELEM.GT.50000)THEN
            WRITE(53,907)TOTL,NELEM
            STOP
          ENDIF
          ITYPE(NELEM)=I

  100     CONTINUE
          RPOSOLD=RPOSNEW
C=====    FOR HORIZONTAL DIPOLES or HORIZONTAL CF SLOT IN THE EDGE FIELD QUADS: F/B
C=====    For split dipoles/CF, only need EDGE fields at the front of th efirst and behind
C=====    the second. Leave it for now: 26/12/03. Need a new flag in the NAMELIST to do
C=====    a clean switch off instead of killing in HEDGE.
C=====    ALSO SLOT IN A VERTICAL KICKER AT EVERY SECOND DIPOLE SLICE AFTER THE
C=====    FIRST-- UP TO A TOTAL OF 4 KICKERS.
C         IF(IREM.EQ.2.OR.IREM.EQ.15)THEN
C           NELEM=NELEM+1
C            LL=2
C            WRITE(53,908)LL,TOTL,NELEM
C           ITYPE(NELEM)=IN+1
C         ENDIF
          NKI=0
C=====    =August 2003: Kill kicks inserted for twist.
          DO INU=1,NU
            NELEM=NELEM+1
C            LL=3
C            WRITE(53,908)NU,TOTL,NELEM
            ITYPE(NELEM)=IN
            CENPOS(NELEM) = RPOS     ! Store the positions, but just at real elements.
            IF(1.EQ.1) CYCLE
            IF(KICK.EQ.0) CYCLE
            IF(NU.NE.2.OR.INU.NE.1) CYCLE
            IF(IREM.NE.2.OR.IREM.EQ.15) CYCLE
            NKI=NKI+1
            NTY=NTY+1
            NELEM=NELEM+1
            ITYPE(NELEM)=NTY
            NAME(NTY)='VD'
            SNAME(NTY) = NAME(IN)
            XX(NTY)=0.D0
            YY(NTY)=2.D0*YY(IN)            ! Give the kickers the same length as the dipole.
            ID(NTY)=7                      ! That also avoids radiation problems.
C  771       CALL RANNOR(GAUS1,GAUS2)
            KG5 = KG5 + 1
            XX(NTY)=XX(IN)*G5(KG5)*1.D-3 *2.D0 ! *2 to get the bend of the full
            TMAT(4,6,NTY)=+XX(NTY)             ! dipole.
            TMAT(5,3,NTY)=-XX(NTY)
            TMAT(4,7,NTY)=-XX(NTY)
          ENDDO
C         IF(IREM.EQ.2.OR.IREM.EQ.15)THEN
C           NELEM=NELEM+1
C            LL=4
C            WRITE(53,908)LL,TOTL,NELEM
C           ITYPE(NELEM)=IN+1
C         ENDIF

C=====    RE-ANCHOR THE LATTICE ON THIS ELEMENT BEFORE PROCEEDING.
          TOTL=YY(IN)*NUS*0.5D0+PLT
          IF(IREM.EQ.10.OR.IREM.EQ.14)TOTL=YY(IN)*NUS+PLT
          IF(IREM.EQ. 6.OR.IREM.EQ. 7.OR.IREM.EQ.17)TOTL=PLT
          IF((IREM.EQ.3.AND.NAME(IN)(1:2).EQ.'CQ').OR.
     +       (IREM.EQ.4.AND.NAME(IN)(1:2).EQ.'RQ'))TOTL=PLT
C    +       (IREM.EQ.4.AND.NAME(IN)(1:2).EQ.'RQ').OR.
C    +       (IREM.EQ.3.AND.NAME(IN)(1:2).EQ.'EE'))TOTL=PLT
C
C=====    =Set the strength of v. kickers used to shift quads.
C          Allow the quads a maximum misplacement of WID1*SCUT1*1.D-3 m.
C          Pick the strengths from the pre-calculated list.
C          Need the most recent quad strength.
C          Update length with the quad length.
C          Set safe parameters at start in case it has no length.
C          A +ve quad strength implies vertical defocusing. Thus the vertical kick
C          to a betatron trajectory is (quad strength)*(betatron trajectory position).
C          Then for a quad shifted up by G1, the change to the vertical kick is -(quad strength)*G1.
C          However, with the Mais-Ripken sign convention used here, a negative
C          vertical kick implies a +ve kicker strength.
C          A horizontal offset leads to a - sign -- see the HQ stuff.
          IF(IREM.EQ.7.AND.NAME(IN)(1:2).EQ.'VQ')THEN
            XX(IN)=0.D0
            YY(IN)=QYREM
            TMAT(4,6,IN)=0.D0
            TMAT(5,3,IN)=0.D0
            TMAT(4,7,IN)=0.D0
            IF(KICK.EQ.1)THEN
C  777         CALL RANNOR(GAUS1,GAUS2)
              KG1 = KG1 + 1
              XX(IN)=QREM*G1(KG1)*1.D-3    ! Note the plus sign.
              YY(IN)=QYREM
              WRITE(*,*) NAME(INREM), G1(KG1)*1.D-3
              TMAT(4,6,IN)=+XX(IN)
              TMAT(5,3,IN)=-XX(IN)
              TMAT(4,7,IN)=-XX(IN)
            ENDIF
          ENDIF
C=====    =Set the strength of h. kickers used to shift quads.
C          Set safe parameters at start,in case it has no length.
          IF(IREM.EQ.6.AND.NAME(IN)(1:2).EQ.'HQ')THEN
            XX(IN)=0.D0
            YY(IN)=QYREM
            TMAT(2,6,IN)=0.D0
            TMAT(5,1,IN)=0.D0
            TMAT(2,7,IN)=0.D0
            IF(KICK.EQ.1)THEN
              KG2 = KG2 + 1
              XX(IN)= -QREM*G2(KG2)*1.D-3      ! Note the minus sign.
              YY(IN)=  QYREM
              TMAT(2,6,IN)= XX(IN)
              TMAT(5,1,IN)=-XX(IN)
              TMAT(2,7,IN)=-XX(IN)
            ENDIF
          ENDIF
C=====    =Set the strength of skew quad used to rotate quads. Set safe parameters at start,
C          in case it has no length.
          IF(IREM.EQ.4.AND.NAME(IN)(1:2).EQ.'RQ')THEN
            XX(IN)=0.D0
            YY(IN)=QYREM
            TMAT(2,3,IN)=0.D0
            TMAT(4,1,IN)=0.D0
            IF(KICK.EQ.1)THEN
              KG3 = KG3 + 1
              XX(IN)=QREM*G3(KG3)*1.D-3 *2.D0      ! Approx. sin(2phi)
              YY(IN)=QYREM
              TMAT(2,3,IN)=-XX(IN)
              TMAT(4,1,IN)=-XX(IN)
            ENDIF
          ENDIF
C=====    =Set the fraction quad strength error.Set safe parameters at start,in case it has no length.
          IF(IREM.EQ.3.AND.NAME(IN)(1:2).EQ.'CQ')THEN
            XX(IN)=0.D0
            YY(IN)=QYREM
            TMAT(2,1,IN)=0.D0
            TMAT(4,3,IN)=0.D0
            IF(KICK.EQ.1)THEN
              KG4 = KG4 + 1
              XX(IN)=QREM*G4(KG4)*1.D-3            !   * TWIST(IN)
              YY(IN)=QYREM
              TMAT(2,1,IN)=-XX(IN)
              TMAT(4,3,IN)= XX(IN)
            ENDIF
          ENDIF
C=====    =Set the dipole rotation.Set safe parameters at start,in case it has no length.
          IF(IREM.EQ.7.AND.NAME(IN)(1:2).EQ.'VD')THEN
            XX(IN)=0.D0
            YY(IN)=DYREM
            TMAT(4,6,IN)=0.D0
            TMAT(5,3,IN)=0.D0
            TMAT(4,7,IN)=0.D0
            IF(KICK.EQ.1)THEN
              KG5 = KG5 + 1
              XX(IN)=DREM*G5(KG5)*1.D-3
              YY(IN)=DYREM
              TMAT(4,6,IN)=+XX(IN)
              TMAT(5,3,IN)=-XX(IN)
              TMAT(4,7,IN)=-XX(IN)
            ENDIF
C=====    Set the kick scales for the tracking. Associate kicks with VDs/sections
            RADSCALE(IN) = DABS(DREM)**3/DYREM**2 !Used this for eRHIC but without the sqrt.
          ENDIF

C          CENPOS(NELEM) = RPOS(KK)   ! Explicitly store the positions.
C          IF(NELEM.LT.1000)write(*,'(A,A,I10,F15.6)')' ',
C         +                     NAME(IN),NELEM,CENPOS(NELEM)
        END DO
      END DO
C======Jump to here when the structure has been read (EOF found).
      NTYPE=NTY
      WRITE(53,902)NELEM,NTY
C
C
C
C=====PRINT OUT TYPE LIST: EXCLUDE VERTICAL KICKERS FOR NOW? No: kill by excluding type 77.
      WRITE(53,'(A,A)')'1','Summarise the final type list.'
      WRITE(53,954)
      DO I=1,NTY
        WRITE(53,97)I,ID(I),NAME(I),XX(I),X2(I),YY(I),NUNT(I),TWIST(I),
     +              SNAME(I) ,NSOL(I),NTWIST(I)
        IF(ID(I).EQ.1.AND.XX(I).LT.0.D0)THEN
        WRITE(53,98)
        STOP
        ENDIF
      ENDDO
C
C
C
C
C
      WRITE(53,'(A,A)')'1','Sample of the final structure.'

C=====NUMBER OF MAGNETS AND LENGTH OF RING
      IBEND =0
      ICOMF =0
      IQUAD =0
      IQUADC=0
      ISKQ  =0
      ISEXT =0
      IHCOR =0
      IHD   =0
      IHQ   =0
      IHC   =0
      IVCOR =0
      IVD   =0
      IVQ   =0
      IVC   =0
      ICAV=0
      RLENGE=0.D0
      DO 69 I=1,NELEM
C     IF(I.LE.5000)THEN
        WRITE(53,60)I,NAME(ITYPE(I)),ID(ITYPE(I)),XX(ITYPE(I)),
     +           X2(ITYPE(I)),YY(ITYPE(I)),RADSCALE(ITYPE(I)),RLENGE
C     ENDIF


      JID=ID(ITYPE(I))
      GO TO(61,62,62,62,71,67,67,62,62,72,62,62,62,72,62,62,62),JID
      GO TO 69
   61 RLENGE=RLENGE+XX(ITYPE(I))*1.D0
      GO TO 69
   62 IF((NAME(ITYPE(I))(1:2).EQ.'CQ'.AND.JID.EQ.3) !Ignore artificial lengths.
     +.OR.(NAME(ITYPE(I))(1:2).EQ.'RQ'.AND.JID.EQ.4).OR.JID.EQ.17)
     +                                                       GO TO 699
      RLENGE=RLENGE+YY(ITYPE(I))*1.D0
  699 CONTINUE
      IF(JID.EQ. 2.OR.JID.EQ. 9)IBEND=IBEND+1
      IF(JID.EQ.15.OR.JID.EQ.16)ICOMF=ICOMF+1
      IF(JID.EQ. 3.AND.NAME(ITYPE(I))(1:2).NE.'CQ')IQUAD=IQUAD+1
      IF(JID.EQ. 3.AND.NAME(ITYPE(I))(1:2).EQ.'CQ')IQUADC=IQUADC+1
      IF(JID.EQ. 4             )ISKQ =ISKQ +1
      IF(JID.EQ. 8             )ISEXT=ISEXT+1
      IF(JID.EQ.17.OR.NAME(ITYPE(I))(1:2).EQ.'BB')then
        IBB=IBB+1
        write(*,*) JID, ITYPE(I)
      endif
      GO TO 69
   67 RLENGE=RLENGE+YY(ITYPE(I))*0.D0       !Kickers are given zero length.
      IF(JID.EQ.7)IVCOR=IVCOR+1
      IF(JID.EQ.7.AND.NAME(ITYPE(I))(1:2).EQ.'VD')IVD=IVD+1
      IF(JID.EQ.7.AND.NAME(ITYPE(I))(1:2).EQ.'VQ')IVQ=IVQ+1
      IF(JID.EQ.7.AND.NAME(ITYPE(I))(1:2).EQ.'VC')IVC=IVC+1
      IF(JID.EQ.6)IHCOR=IHCOR+1
      IF(JID.EQ.6.AND.NAME(ITYPE(I))(1:2).EQ.'HD')IHD=IHD+1
      IF(JID.EQ.6.AND.NAME(ITYPE(I))(1:2).EQ.'HQ')IHQ=IHQ+1
      IF(JID.EQ.6.AND.NAME(ITYPE(I))(1:2).EQ.'HC')IHC=IHC+1
      GO TO 69
   71 ICAV=ICAV+1
      GO TO 69
   72 RLENGE=RLENGE+YY(ITYPE(I))*1.D0

C   69 WRITE(53,'(A,A,2I10,F15.5)')' RLENGE ',NAME(ITYPE(I)),I,JID,RLENGE
   69 CONTINUE
      WRITE(53,70)IBEND,ICOMF,IQUAD,IQUADC,ISKQ,ISEXT,IBB,
     +IHCOR,IHD,IHQ,IHC,
     +IVCOR,IVD,IVQ,IVC,
     +ICAV, RLENGE,TOTL,TOTDIP
C
C      IF(1.EQ.1)STOP
C
C
C
C
C=====TEST TOTAL BEND ANGLES TO SEE IF THEY GIVE ZERO OR 2*PI
      AH = 0.D0
      AV = 0.D0
      DO 23 II=1,NELEM
      ITY=ITYPE(II)
      IF(ID(ITY).EQ.2.OR.ID(ITY).EQ.15)AH=AH+XX(ITY)
      IF(ID(ITY).EQ.9.OR.ID(ITY).EQ.16)AV=AV+XX(ITY)
   23 CONTINUE
      PIE2=2.*PI
      WRITE(53,905) AH,AV,PIE2
C
C
C
C=====FOR RESTARTING RANDOM NUMBER GENERATOR
C      CALL RDMOUT(JSEED)
C      WRITE(53,956)JSEED
C
C
C     IF(1.EQ.1)STOP
      RETURN
C
C
C
C
   60 FORMAT(I6,2X,A8,I4,3(F14.5),E15.6,F14.5)
C
   70 FORMAT(/,' NUMBER OF BENDING MAGNETS :          ',I4,/,
     +' NUMBER OF COMB. FUNC.(CF) MAGNETS  : ',I4,/,
     +' NUMBER OF QUADRUPOLES              : ',I4,/,
     +' NUMBER OF QUAD AND CF PERTURBATIONS: ',I4,/,
     +' NUMBER OF SKEW QUADRUPOLES         : ',I4,/,
     +' NUMBER OF SEXTUPOLES               : ',I4,/,
     +' NUMBER OF BEAM-BEAM ELEMENTS       : ',I4,/,
     +' ',/,
     +' NUMBER OF HORIZONTAL KICKERS       : ',I4, ' and of these: ',/,
     +' NUMBER FOR HOR. DIPOLE ERRORS      : ',I4,/,
     +' NUMBER FOR HOR. QUAD SHIFTS        : ',I4,/,
     +' NUMBER FOR HOR. CORRECTION         : ',I4,/,
     +' ',/,
     +' NUMBER OF VERTICAL KICKERS         : ',I4, ' and of these: ',/,
     +' NUMBER FOR VERT. DIPOLE ROLE       : ',I4,/,
     +' NUMBER FOR VERT. QUAD SHIFTS       : ',I4,/,
     +' NUMBER FOR VERT. CORRECTION        : ',I4,/,
     +' ',/,
     +' NUMBER OF CAVITIES                 : ',I4,/,
     +' ',/,
     +' LENGTHS OF THE RING                : ',2F12.4,/,
     +' LENGTH OF HOR. DIPOLES             : ',1F10.4,/)
C
   80 FORMAT(80A1)
   91 FORMAT(' NAME ',A8,' NOT FOUND IN THE TYPE LIST')
   92 FORMAT(6(A4,I7,1X))
  922 FORMAT(5(A4,1X,I8,1X))
 9222 FORMAT(4(A8,F16.8,1X))
C  95 FORMAT(5(4A1,1X,I8,1X))
   93 FORMAT('1',' THE LATTICE AS READ FROM THE INPUT DATA---')
   94 FORMAT(' SIZE OF TYPE LIST EXCEEDED !!!')
   95 FORMAT(6(1X,A4,I7,1X))
  955 FORMAT(5(1X,A4,1X,I9,1X))
 9555 FORMAT(1X,A8,1X,F16.8,1X)
   96 FORMAT('  FALSE MAGNET TYPE  ',I6)
   97 FORMAT(' ',2I5,1X,A8,3F12.8,I5,2X,F11.6,4X,A8,I4,I10)
   98 FORMAT(' ','A NEGATIVE DRIFT LENGTH HAS BEEN FOUND--SO STOP')
C   99 FORMAT(  A1,I4,1X,A8,3F12.8,I5   ,F11.6,I5)
   99 FORMAT(  A1,I4,1X,A8,3E12.8,I5   ,F11.6,I5)
C 199 FORMAT(  A1,I4,1X,2X,I2,1A1,3X,3F12.8,I5   ,F11.6,I5)
  199 FORMAT(8X,I2,1A1,3X)
   89 FORMAT(  A1,I4,1X,A8,3F12.8,I5   ,F11.6,I5)
  801 FORMAT(9X,I1)
  802 FORMAT(' ','LATTICE FORMAT= ',I1,/,
     +           'NEW RANDOM SEED= ',I20,/,
     +           'DEF.RANDOM SEED= ',I20)
  803 FORMAT(' ',' STOP--BAD LATTICE FORMAT ')
  805 FORMAT(8X,I2)
  901 FORMAT(20A4)
  902 FORMAT(' LATTICE FILLED UP TO ',I7,', TYPE LIST FILLED UP TO ',I7,
     +       /,/,/)
C
  905 FORMAT(' TOTAL BENDING ANGLE',
     + /' HORIZONTAL ',F15.10,
     + /' VERTICAL   ',F15.10,
     + /' 2.*PI      ',F15.10,/)
C
  907 FORMAT(' SIZE OF LATTICE LIST EXCEEDED !!!',/' LTOT =',F11.4,
     +' ELEMENT=',I6)
  908 FORMAT(' NELEM LL ', I10, ' LTOT =', F11.4,' ELEMENT=',I6)
  950 FORMAT(' SOMETHING IS WRONG WITH THE INPUT DATA-FILE ')
  951 FORMAT(1X,80A1)
  952 FORMAT('0',//,'  TYPE NAME      STRENGTH1  STRENGTH2   LENGTH',
     +  '   NO.SLICES  SFLAG GMATRIX FLAG')
  953 FORMAT('0','SANDWICH STRUCTURE FOCUSSING DIPOLES NOT ALLOWED')
  954 FORMAT('0',//,'       TYPE NAME      STRENGTH1   STRENGTH2  ',
     +' LENGTH      SLICES  SFLAG                  GMATRIX ')
  956 FORMAT('0','LAST RANDOM GENERATOR SEED=',I20)
      END




C      IF(NAMEREM(1:2).EQ.'BW')RADSCALE(IN)= RADSCALE(IN)*0.0
C      RADSCALE(IN) = DREM**2/DYREM
C      IF(NAMEREM(1:2).EQ.'BM')RADSCALE(IN)= RADSCALE(IN)*2.5
C      IF(NAMEREM(1:2).EQ.'BM')RADSCALE(IN)= RADSCALE(IN)*2.5
C      IF(NAMEREM(1:2).EQ.'BA')RADSCALE(IN)= RADSCALE(IN)*2.0
C      IF(NAMEREM(1:2).EQ.'BA')RADSCALE(IN)= RADSCALE(IN)*2.0
C      IF(NAMEREM(1:2).EQ.'BM')
C     + write(53,'(A,I10,F15.6)')' Lattin:reached here ',IN,RADSCALE(IN)

C      WRITE(53,'(A,A8,1X,2I10,F15.6)')
C     +                        ' VD on?',NAME(IN),IREM,KG5,DREM
C      WRITE(53,'(A,I10,4A,3E16.5)')
C     +  ' Lattin',IN,'  ',NAME(IN),'  ',NAMEREM,RADSCALE(IN),DREM,DYREM
