C
C     August 2010. At last recompiling with the special 32-bit tricks in the makefile.
C     So now tweek the alignments as in EpsSLICK.
      COMMON/CNLIST/NINS,IRAD,ISPIN,IPTCO,IDISP,IPTNML,ISIG,NSTEP,
     +              NITER,IHARM,KICK,IPTBET,ICLBP,ITWP,NFREQ,
     +              ISECT,IPTD,MRDISP,IVROT,ICORR,IKMIN,IEDGE,ISEXTU,
     +              ILIN,IEXTSIZE,NDIFFDUMP,NSKIP,IBMBM,
C=====First the integers. Need an even number of 4 bytes.
     +              KSEED1,  KPAD01, 
     +              KSEED2,  KPAD02, 
     +              KSEED3,  KPAD03, 
     +              KSEED4,  KPAD04, 
     +              KSEED5,  KPAD05, 
     +              KSEED6,  KPAD06, 
     +              KSEED7,  KPAD07, 
     +              KSEED8,  KPAD08, 
     +              KSEED9,  KPAD09, 
     +              KSEED10, KPAD10, 
     +              KSEED11, KPAD11, 
     +              KSEED12, KPAD12, 
     +              KSEED13, KPAD13, 
     +              KSEED14, KPAD14, 
     +              KSEED15, KPAD15, 
     +              KSEED16, KPAD16, 
     +              KSEED17, KPAD17, 
     +              KSEED18, KPAD18, 
     +              KSEED19, KPAD19, 
     +              KSEED20, KPAD20, 
     +              IBEQUIL, NPART1, NTURN1, NDAMP1,
     +              IDEPLIN, NPART2, NTURN2, NDAMP2,
     +              IDEPNON, NPART3, NTURN3, NDAMP3, MODES,KPAD00,
C===== Now the reals
     +              WID1, SCUT1,
     +              WID2, SCUT2,
     +              WID3, SCUT3,
     +              WID4, SCUT4,
     +              WID5, SCUT5,
     +              WID6, SCUT6,
     +              WID7, SCUT7,
     +              WID8, SCUT8,
     +              WID9, SCUT9,
     +              WID10,SCUT10,
     +              WID11,SCUT11,
     +              WID12,SCUT12,
     +              WID13,SCUT13,
     +              WID14,SCUT14,
     +              WID15,SCUT15,
     +              WID16,SCUT16,
     +              WID17,SCUT17,
     +              WID18,SCUT18,
     +              WID19,SCUT19,
     +              WID20,SCUT20,
     +              E00,DE0,ECAV,CIRCUM,RADFAC,SFREQ,DFREQ,SCLKIK,ECHROM
C
C
C
C     +              KSEED1,  WID1, SCUT1,
C     +              KSEED2,  WID2, SCUT2,
C     +              KSEED3,  WID3, SCUT3,
C     +              KSEED4,  WID4, SCUT4,
C     +              KSEED5,  WID5, SCUT5,
C     +              KSEED6,  WID6, SCUT6,
C     +              KSEED7,  WID7, SCUT7,
C     +              KSEED8,  WID8, SCUT8,
C     +              KSEED9,  WID9, SCUT9,
C     +              KSEED10 ,WID10,SCUT10,
C     +              KSEED11 ,WID11,SCUT11,
C     +              KSEED12 ,WID12,SCUT12,
C     +              KSEED13 ,WID13,SCUT13,
C     +              KSEED14 ,WID14,SCUT14,
C     +              KSEED15 ,WID15,SCUT15,
C     +              KSEED16 ,WID16,SCUT16,
C     +              KSEED17 ,WID17,SCUT17,
C     +              KSEED18 ,WID18,SCUT18,
C     +              KSEED19 ,WID19,SCUT19,
C     +              KSEED20 ,WID20,SCUT20,
C     +              IBEQUIL, NPART1, NTURN1, NDAMP1,
C     +              IDEPLIN, NPART2, NTURN2, NDAMP2,
C     +              IDEPNON, NPART3, NTURN3, NDAMP3,MODES

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
C   (  "      1)   PRINTS OUT CLOSED-ORBIT 6-VECTOR.
C   (  "      2)   CALCULATES LINEAR OPTICS TO GET SOME PARAMETERS & THEN
C   (  "      2)   SETS C.O.TO ZERO & JUMPS DIRECTLY TO SPIN CALCULATION
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
C   (ISECT.EQ.1)   GET 8X8 MATRIX FOR A SECTION OF RING.
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
C   (KSEED4,  WID4, SCUT3)  QUAD ROLE, WIDTH, GAUSSIAN CUT
C   (KSEED5,  WID5, SCUT4)  QUAD FRACTIONAL STRENGTH ERROR
C   (KSEED3,  WID3, SCUT5)  DIPOLE ROLE, WIDTH, GAUSSIAN CUT
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
C    RADFAC           SCALE FOR SIZE OF BIG PHOTONS.
C   =================================================================
