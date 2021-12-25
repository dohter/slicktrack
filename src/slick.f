C   19/05/82 603081607  MEMBER NAME  SLICK    (SEPT95.S)    FORTRAN
C
C
       PROGRAM SLICK

C
C
C                 MAIN ROUTINE FOR RUNNING SLICKTRACK
C                 -----------------------------------
C
C  **************************************************************************************
C  *                                                                                    *
C  *      The history of SLICKTRACK, a request to the user, and a disclaimer.           *
C  *      ------------------------------------------------------------------            *
C  *                                                                                    *
C  *                                                                                    *
C  *                                                                                    *
C  *                                   History                                          *
C  *                                   -------                                          *
C  *                                                                                    *
C  *    SLICKTRACK has it roots in the code SLIM, for estimating the                    *
C  *    equilibrium spin polarisation in electron/positron storage rings.               *
C  *    SLIM was written A.W. Chao in the early 1980's and is based on the.             *
C  *    algorithm described in A.W. Chao, Nucl.Inst.Meth.  180 (1981) 29.               *
C  *                                                                                    *
C  *    The original SLIM code, which uses thin-lens optics, was quickly                *
C  *    adapted for use at DESY by R. Schmidt and J. Kewisch.                           *
C  *                                                                                    *
C  *    In 1983, a major upgrade was then undertaken by D.P. Barber of DESY             *
C  *    in order to to base it on thick-lens optics and to give it a more               *
C  *    convenient structure with respect to the layout of the subroutines.             *
C  *    The new version was called                                                      *
C  *                                                                                    *
C  *                                     SLICK.                                         *
C  *                                                                                    *
C  *    References to the thick lens formalism can be found in:                         *
C  *    D.P. Barber and G. Ripken,                                                      *
C  *    ``Handbook of Accelerator Physics and Engineering'',                            *
C  *    Eds: A.W. Chao and M. Tigner, 3rd edition, World Scientific (2006).             *
C  *    and refernces therein.                                                          *
C  *                                                                                    *
C  *    Both SLIM and SLICK linearise the orbital and spin motion so that               *
C  *    only first order spin-orbit resonances are apparent. This restriction           *
C  *    to first order resonances is a serious limitation at high energy.               *
C  *    Neither code can handle beam-beam forces beyond trivial linearisation.          *
C  *                                                                                    *
C  *    SLICK has therefore been extended by the inclusion of spin-orbit                *
C  *    tracking routines with inclusion of Monte-Carlo algorithms to simulate          *
C  *    photon emission and estimate the rate of spin diffusion and thence              *
C  *    the rate of depolarisation. The resulting code is called                        *
C  *                                                                                    *
C  *                                  SLICKTRACK                                        *
C  *                                                                                    *
C  *    The Monte-Carlo algorithms are purely classical and in the spirit of the        *
C  *    stochastic analysis of H. Mais and G. Ripken, as for example                    *
C  *    in DESY Report 83-62 (1983) and subsequent papers by these                      *
C  *    with D.P. Barber and K. Heinemann.                                              *
C  *                                                                                    *
C  *    At the time of writing this introduction (April 2009), the orbital              *
C  *    motion in the Monte-Carlo simulations is linearised but the                     *
C  *    spin-orbit tracking can be undertaken either with linearised                    *
C  *    spin motion (within the spirit of the above mentioned stochastic                *
C  *    analysis), or with full 3-D spin motion to be described in                      *
C  *    DESY Report 09-15 in preparation. It is the full 3-D algorithm that             *
C  *    delivers the the higher order spin-orbit resonances.                            *
C  *                                                                                    *
C  *    The Monte-Carlo tracking algorithms can include kicks to the electrons          *
C  *    from an oncoming bunch. It is assumed that the oncoming bunch is                *
C  *    unaffected by the collisions but the beam-beam forces can be handled in         *
C  *    various approximations including the full non-linear force on an electron       *
C  *    from oncoming bunches of elliptical cross section.                              *
C  *                                                                                    *
C  *    It is hoped that in the future, non-linear orbital motion in the lattice        *
C  *    can be included but there is no fixed time scale for this.                      *
C  *                                                                                    *
C  *                                                                                    *
C  *                                                                                    *
C  *                               Request to the user                                  *
C  *                               -------------------                                  *
C  *                                                                                    *
C  *    Until April 2009, SLICK and now SLICKTRACK, have been personal codes, developed *
C  *    at DESY by D.P. Barber. They have been carefully tested and the development     *
C  *    has entailed much work and deep knowledge of the subject.                       *
C  *                                                                                    *
C  *    D.P. Barber therefore requests that presentations and published work            *
C  *    based on this code acknowledge the authorship and mention the name SLICKTRACK.  *
C  *                                                                                    *
C  *                                                                                    *
C  *    D.P. Barber also requests to be informed of problems and bugs and requests that *
C  *    these paragraphs of material on the history etc are not removed.                *
C  *                                                                                    *
C  *                                                                                    *
C  *                                                                                    *
C  *                                   Disclaimer                                       *
C  *                                   ----------                                       *
C  *                                                                                    *
C  *    This code has been built and maintained by  D.P. Barber for doing  physics.     *
C  *    There is no manual but there is some guidance in the ...flags.inp file and in   *
C  *    this main program. The code has been kept as clean as possible in order to      *
C  *    minimise errors and to allow the physics to be easily followed. Beyond that,    *
C  *    no attempt has been made to make the code super elegant. Until now, people      *
C  *    at DESY have been the only users.                                               *
C  *                                                                                    *
C  *    Neither D.P. Barber  nor DESY accepts any responsibilty for loss or damage,     *
C  *    whether material, financial  or of any other kind resulting from the use of     *
C  *    SLICKTRACK                                                                      *
C  *                                                                                    *
C  *                                                                                    *
C  **************************************************************************************
C
C
C
C
C        'THICK SLICE' VERSION WITH EXACT MATRICES FOR THICK
C         CONTIGUOUS SLICES SO THAT INTEGRALS CAN BE HANDLED.
C
C         FAST VERSION USING SHORT CUTS FOR DRIFT SPACES.
C
C
C   NTYPE=(NUMBER OF DIFFERENT MAGNET-TYPES) <= 3000
C   WARNING ==> WHEN CHANGING NO. OF TYPES, TRANS1 MUST BE CHANGED
C   BEAM-LINE MUST CONTAIN ALL ELEMENTS IN STORAGE RING.
C   2 <= (TOTAL NUMBER OF ELEMENTS IN THE BEAM-LINE) <= (NO. IN CLOORB)
C
C
C   EACH MAGNET-TYPE IS SPECIFIED BY (ID,NAME,XX,YY ------).
C   ID=  ( 1 , 2 , 3 ,    4   , 5 ,   6    ,   7    ,  8 , 9 ,  10)
C   MEANS(DRF,BBX,QUA,SKEW-QUA,CAV,HOR-KICK,VER-KICK,SEXT,BBY,SOLENOID)
C   ID=  (11         ,12           ,13           )
C   MEANS(ROTAT.QUAD  FOC.HOR.BEND  FOC.VERT.BEND)
C   ID=  (14        ,15          ,16            ,17)
C   MEANS(ROTATOR  HOR.CF      VERT.CF      BEAM-BEAM)
C
C
C   (ID .EQ. 99) DEFINES HALF SUPERPERIOD IN THE UNPERTURBED MACHINE.
C   NINS=(NUMBER OF HALF SUPERPERIODS IN THE BEAM-LINE).
C   (ID .NE. 1,2,3,9) ARE TREATED AS PERTURBATIONS TO THE UNCOUPLED LINEAR OPTICS.
C
C
C   FOR THICK LENS VERSION INPUT(!!!) STRENGTHS XX ARE:
C   XX=(ELEMENT STRENGTH). FOR DIFFERENT ID NUMBERS, XX ARE GIVEN BY:
C      (L, L*BY/BRHO, L*(DBY/DX)/BRHO, L*(DBY/DY)/BRHO,
C     V*COS(PHIS)*H/(R*E0),L*BY/BRHO,L*BX/BRHO,L*(D**2*BY/DX**2)/BRHO,
C     L*BX/BRHO,L*BZ/BRHO)
C   YY=0 ---- IF ID=1,5,99.
C     =(ELEMENT LENGTH) ---- OTHERWISE.
C
C
C
C   T(J,K,I)=(JK-TH ELEMENT OF THE TRANSFORMATION MATRIX FOR
C     THE I-TH MAGNET-TYPE) ---- THICK-LENS ASSUMED.
C
C
C
C   FOR MEANING OF INPUT FLAG SETTINGS (IN NAMELIST) SEE SUBROUTINE
C   LATTIN.
C
C   THE 7-DIMENSIONAL VECTOR IS (X,X',Y,Y',DL,DE/E,1).
C   THE FIRST AND LAST BEAM-LINE ELEMENTS MUST BE DRIFTS WITH XX=0.
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*120 ROTATOR,RING,OUTPUT,GMT,POLARI,AZIFUNCS,GAUSSES
      CHARACTER*120 ELLIPSE1,ELLIPSE2,BEAMDIST,SPINDIST,SPINDIFF
      CHARACTER*120 SPINELLIPSE1,SPINELLIPSE2,MCTDEP,ABRATIO,FFT
      CHARACTER*120 ISPINSPREAD

      DIMENSION TM6A(6,6),TDAMP(6),TUN(6),COVMAT(6,6),SCALEMAT(6,6)
      NAMELIST/FILES/ROTATOR,RING,OUTPUT,GMT,POLARI,AZIFUNCS,GAUSSES,
     +               ELLIPSE1,ELLIPSE2,BEAMDIST,SPINDIST,SPINDIFF,
     +               SPINELLIPSE1,SPINELLIPSE2,MCTDEP,ABRATIO,FFT,
     +               ISPINSPREAD

C======Rotator      = Look-up table of setting for HERA rotators.
C======Ring         = Input file.
C======Output       = Printout
C======Polari       = Polarization (etc) vs. energy.
C======Azifuncs     = Optic and polarization functions plotted vs. azimuth.
C======Gausses      = Gaussian distributions.
C======Ellipse1     = Beam ellipse from theoretical covariance matrix.
C======Ellipse2     = Beam ellipse from covariance matrix of a tracked ensemble.
C======Beamdist     = Distributions of a tracked ensemble.
C======Spindist     = Distribution of the spin angles alpha and beta.
C======Ispinspread  = Equilibrium distribution of the spin angles alpha and beta.
C======Spindiff     = Time dependence of the spread of the spin angles alpha and beta.
C======
C======MCtdep       = The depolarising times from the Monte-Carlo simulation.
C======ABratio      = Ratio of first turn alph-beta covariances from SPIN and SCRURITA2.
C======Spinellipse1 = SLICK constant density contour for alpha and beta.
C======Spinellipse2 = M-C   constant density contour for alpha and beta.
C======FFT          = FFT of the horizontal and vertical transverse motion for
C                     1 particle at the first energy step. Includes damping.
C
C
      INCLUDE "cloorb.for"
      INCLUDE "cnlist.for"
      INCLUDE "clatic.for"
      INCLUDE "csol.for"
C
C
C
C
      PI=3.1415926535897932D0
      ROTATOR       = 'slick_rotator.inp'
      RING          = 'he7basml.inp'
      OUTPUT        = 'slick.out'
      GMT           = 'slick.gmt'
      POLARI        = 'slick_pol.out'
      AZIFUNCS      = 'slick_dndelta_funcs.out'
      GAUSSES       = 'slick_gausses.out'
      ELLIPSE1      = 'slick_ellipse1.out'
      ELLIPSE2      = 'slick_ellipse2.out'
      BEAMDIST      = 'slick_beamdist.out'
      SPINDIST      = 'slick_spindist.out'
      SPINDIFF      = 'slick_spindiff.out'
      MCTDEP        = 'slick_mctdep.out'
      ABRATIO       = 'slick_abratio.out'
      SPINELLIPSE1  = 'slick_spinellipse1.out'
      SPINELLIPSE2  = 'slick_spinellipse2.out'
      FFT           = 'slick_fft.out'
      ISPINSPREAD   = 'slick_ispinspread.out'

      READ(5,FILES)
      OPEN(51,FILE=ROTATOR,      STATUS='OLD')
      OPEN(52,FILE=RING,         STATUS='OLD')
      OPEN(53,FILE=OUTPUT,       STATUS='UNKNOWN')

C    *************** Fanglei comments some of the following open files on 20Dec12 ********
C    *************** Fanglei uncomments on 12May20 ********
      OPEN(533,FILE=GMT,         STATUS='UNKNOWN')
      OPEN(54,FILE=POLARI,       STATUS='UNKNOWN')
      WRITE(54,9467)"#", "E0", "AGAMMA","NU0",
     +              "PBKS","PTOT","PX","PY","PS",
     +              "TAUBKS","TAUD0","TAUDX","TAUDY","TAUDS",
     +              "TAUTOT"
      WRITE(54,9468)"#",[(I, I=1,14)]
     +
     +
     +
9467  FORMAT(' ',(A1,A14,3X,13(A15,3X)))
9468  FORMAT(' ',(A1,I14,3X,13(I15,3X)))
      OPEN(55,FILE=AZIFUNCS,     STATUS='UNKNOWN')
C      OPEN(56,FILE=GAUSSES,      STATUS='UNKNOWN')
      OPEN(57,FILE=ELLIPSE1,     STATUS='UNKNOWN')
C      OPEN(58,FILE=ELLIPSE2,     STATUS='UNKNOWN')
C      OPEN(59,FILE=BEAMDIST,     STATUS='UNKNOWN')
      OPEN(60,FILE=SPINDIST,     STATUS='UNKNOWN')
      OPEN(61,FILE=SPINDIFF,     STATUS='UNKNOWN')
C      OPEN(62,FILE=SPINELLIPSE,  STATUS='UNKNOWN')
      OPEN(63,FILE=MCTDEP,       STATUS='UNKNOWN')
C      OPEN(64,FILE=ABRATIO,      STATUS='UNKNOWN')
C      OPEN(65,FILE=SPINELLIPSE1, STATUS='UNKNOWN')
C      OPEN(66,FILE=SPINELLIPSE2, STATUS='UNKNOWN')
C      OPEN(67,FILE=FFT,          STATUS='UNKNOWN')
      OPEN(68,FILE=ISPINSPREAD,  STATUS='UNKNOWN')
C
C
C
C
C
C


C    ***************  CHROTN is reading the HERA rotator look-up table *****
      CALL CHROTN(0.D0,0,0)

      CALL LATTIN(RLENGE,IROT)


C      IF(1.EQ.1)STOP
C
C
C
C
C
C
C
C     ***************************************************************
C     *  LOOP OVER THE ENERGY STEPS.                                *
C     *  TWO MODES ARE RUNNING.                                     *
C     *                                                             *
C     *  MODE 1 :                                                   *
C     *                                                             *
C     *  1)FIRST CALCULATE THE UNPERTURBED LINEAR OPTICS.           *
C     *  2)THEN CALCULATE THE PERTURBED 6X6 OPTICS.                 *
C     *  3)THEN CALCULATE PERTURBED DAMPING AND EXCITATION.         *
C     *  4)THEN,IF IRAD=1 DO THE SPIN CALCULATION IN SPIN.          *
C     *  5)REPEAT FROM 1).                                          *
C     *                                                             *
C     *  MODE 2 :                                                   *
c     *                                                             *
C     *  OTHERWISE, DO 1) TO 4) JUST ONECE, I.E. TO GET THE         *
C     *  PERTURBED 6X6 OPTICS AND THE EMITTANCE BUT NO SPIN STUFF.  *
C     *  THEN, FOR THE SPIN DO A SEPARATE LOOP OVER ENERGY (SINCE   *
C     *  THE CLOSED ORBIT DOESN'T CHANGE IF ENERGY JUMPS IN THE     *
C     *  MAGNETS ARE NEGLECTED). HOWEVER THE CAVITY VOLTAGE MUST    *
C     *  BE RECALCULATED IN THE NEW SPIN LOOP.                      *
C     *  ALSO, IF IPTCO=2, EVEN SKIP THE PERTURBED 6X6 OPTICS AND   *
C     *  THE EMITTANCE.                                             *
C     *                                                             *
C     *  MODE 3 :                                                   *
C     *                                                             *
C     *  IF IPTCO=2, THEN JUST LOOP THE SPIN.                       *
C     ***************************************************************
C
C
C
      IF(IRAD.EQ.1.AND.IPTCO.EQ.2)THEN
      WRITE(53,'(A)')' Incompatible values of IRAD and IPTCO. So STOP'
      STOP
      ENDIF
C
C
C=====For crazy normalisations in SPIN.F
      NCRAZY = 0

C=====ENERGY LOOP---NUMBER OF STEPS =1  UNLESS WE LOOP OVER CLOSED ORBITS.
      ISTEP=1
      IF(IRAD.EQ.1)ISTEP=NSTEP
C
C
      DO 50 IE0=1,ISTEP
      E0=E00+(IE0-1)*DE0
C=====CHANGE ROTATOR SETTINGS AT EACH ENERGY.
      IF(IVROT.EQ.0)GO TO 611
      CALL CHROTN(E0,0,1)
      DO 622 I=1,NTYPE
      IID=ID(I)
      IF(IID.EQ.2.OR.IID.EQ.9)CALL CHROTN(E0,I,2)
      IF(IID.EQ.2)CALL DIPH(XX(I),YY(I),TMAT(1,1,I))
      IF(IID.EQ.9)CALL DIPV(XX(I),YY(I),TMAT(1,1,I))
  622 CONTINUE
  611 CONTINUE

      CALL CAVITY(E0,U0,RLENGE,CRAD)
      U0KEEP=U0
C
C      IF(1.EQ.1)STOP

C
      IF(ISECT.EQ.0.AND.ILIN.NE.0)
     +                              CALL LINOPT(E0,IE0,U0,CIR)

C      IF(1.EQ.1)STOP

      IF(ISECT.EQ.1.AND.ILIN.NE.0)
     +                              CALL LINOPENOPT(E0,IE0,U0,CIR)
      CIR=RLENGE

C      IF(1.EQ.1)STOP
C
      IF(ITWP.NE.0)STOP
      IF(IPTCO.EQ.2)GO TO 50    !Have done 1 energy step with just LINOPT.
C                          !Now go to the special SPIN loop.
C                          !IPTCO = 2 makes no sense if IRAD = 1.
C                          !KEEP GOING FOR JUST ONE ENERGY STEP IF IPTCO = 1 or 0.
C
C
C
C
C
C
C=====INCLUDE ENERGY JUMP IN EACH MAGNET: CORRECTED 25/8/83 TO FORM
C=====TMAT(6,7) BEFORE TMAT(2,7)&TMAT(4,7). AT 2ND & LATER ENERGY STEPS
C=====REDEFINE THE (2,7) & (4,7) ELEMENTS BACK TO THEIR ENERGY INDEP.
C=====VALUES BEFORE ADDING ON THE NEW (6,7) TERM.
      IF(IRAD.EQ.0)GO TO 55
      DO 140 I=1,NTYPE
      IID=ID(I)
      GO TO(140,141,140,140,143,144,145,140,142,140,140,140,140,140,
     +                                               141,142,140),IID
      GO TO 140
C=====H. BENDS: M & R CONVENTIONS---ALSO INCLUDE C.F MAGNETS TEMPORARILY
  141 TMAT(1,7,I)=-CRAD*(XX(I)-DSIN(XX(I)))
      TMAT(2,7,I)=-CRAD*XX(I)/YY(I)*(1.D0-DCOS(XX(I)))
      TMAT(5,7,I)=-CRAD*(XX(I)*XX(I)/2.D0-1.D0+DCOS(XX(I)))
      TMAT(6,7,I)=-CRAD*XX(I)**2/YY(I)
      GOTO 140
C=====V. BENDS
  142 TMAT(3,7,I)=-CRAD*(XX(I)-DSIN(XX(I)))
      TMAT(4,7,I)=-CRAD*XX(I)/YY(I)*(1.D0-DCOS(XX(I)))
      TMAT(5,7,I)=-CRAD*(XX(I)*XX(I)/2.D0-1.D0+DCOS(XX(I)))
      TMAT(6,7,I)=-CRAD*XX(I)**2/YY(I)
      GOTO 140
C=====H. KICKER
  144 IF(IE0.GT.1)TMAT(2,7,I)=-XX(I)*TMAT(6,7,I)/2.D0+TMAT(2,7,I)
      TMAT(6,7,I)=-CRAD*XX(I)**2/YY(I)
      TMAT(2,7,I)= XX(I)*TMAT(6,7,I)/2.D0+TMAT(2,7,I)
      GOTO 140
C=====V. KICKER
  145 IF(IE0.GT.1)TMAT(4,7,I)=-XX(I)*TMAT(6,7,I)/2.D0+TMAT(4,7,I)
      TMAT(6,7,I)=-CRAD*XX(I)**2/YY(I)
      TMAT(4,7,I)=+XX(I)*TMAT(6,7,I)/2.D0+TMAT(4,7,I)
      GOTO 140
C=====CAVITY TERM: IF DIPOLE RADIATION EFFECTS ARE SWITCHED OFF THIS
C=====SHOULD NOT BE USED SINCE IT REPRESENTS A SUPPLY OF ENERGY WHICH IS
C=====NOT DISSIPATED BY RADIATION. IT IS NOT THEN POSSIBLE TO FIND
C=====A 6-DIM CLOSED ORBIT SINCE THE ENERGY WOULD INCREASE TURN BY TURN.
  143 TMAT(6,7,I)=YY(I)
  140 CONTINUE
C
C
C
C
C
C
C
C
C
   55 CONTINUE

      IF(IRAD.eq.0) U0=U0KEEP   ! Oleksii NOW U0 is updated only if
C                                 IRAD=0
C
C
C      IF(ISECT.NE.0)GO TO 12345

C
C       write(*,'(A,I10,F15.6)')' Slick: reached here 2',ICORR
C
C=======Get and correct the closed orbit. But only if
      IF(IPTCO.EQ.1.AND.ICORR*KICK*(KSEED1+KSEED2+KSEED5).NE.0)THEN
      CALL FIXORB1(E0,IE0,U0,CRAD,FREQ,KKICK)   ! Get roughly corrected orbit
C       write(*,'(A,I10,F15.6)')' Slick: reached here 3'
C      CALL FIXORB2(E0,IE0,U0,CRAD,FREQ,KKICK)   ! Do a fine correction
C      CALL FIXORB5(E0,IE0,U0,CRAD,FREQ,KKICK)   ! Test the outcome.
      ENDIF

C      IF(1.EQ.1)STOP


      IF(IPTCO.EQ.1.AND.IKMIN*KICK*(KSEED1+KSEED2+KSEED5).NE.0)THEN
      CALL FIXORB1(E0,IE0,U0,CRAD,FREQ,KKICK)   ! Get roughly corrected orbit
C      CALL FIXORB5(E0,IE0,U0,CRAD,FREQ,KKICK)
      CALL FIXORB3(E0,IE0,U0,CRAD,FREQ,KKICK)
      CALL FIXORB5(E0,IE0,U0,CRAD,FREQ,KKICK)   ! Test the outcome.
      ENDIF
C=======Construct the corrected closed orbit to be used for the rest.
C      IF(1.EQ.1)STOP
      IF(IPTCO.EQ.1)CALL ORBIT(E0,IE0,U0,CRAD,FREQ,KKICK)

      IF(ICLBP.EQ.2)STOP

C12345 CONTINUE

      IF(IPTCO.EQ.1)GO TO 267
C======Set the C.O. to zero if IPTCO = 0.
      DO 266 I=1,NELEM
      DX(I)=0.D0
      DY(I)=0.D0
      DXP(I)=0.D0
      DYP(I)=0.D0
      DL(I)=0.D0
      DEL(I)=0.D0
  266 CONTINUE
  267 CONTINUE
C
C
C
C      IF(1.EQ.1)STOP
C

C      IF(ISECT.NE.0)GO TO 54321

      CALL PERTOP(E0,IE0,U0,CIR,TM6A,TUN)
      IF(MRDISP.EQ.1)CALL DERTOP(E0,IE0,U0,CIR,TM6A,TUN)

C
C      IF(1.EQ.1)STOP

      CALL DAMPER(E0,IE0,U0,CRAD,CIR,TDAMP,TUN,ANTIDP)
      IF(ANTIDP.LT.0)THEN
      WRITE(53,1063)
      STOP
      ENDIF
C
C
C      IF(1.EQ.1)STOP
      CALL EMITNC(E0,IE0,CIR,CRAD,TDAMP,TUN,TM6A,COVMAT,SCALEMAT,
     +                                                OFFSETL,OFFSETE)

C      IF(1.EQ.1)STOP
C=====IF IRAD =1, GET THE POLARIZATION IN THIS LOOP.
      TAUY=1.
C      IF(1.EQ.1)STOP

      IF(ISPIN.EQ.1.AND.IRAD.EQ.1)
     +                 CALL SPIN(IE0,E0,CIR,TAUY,BETAY0,TUN,NCRAZY)

      IF(ISPIN.EQ.1.AND.IRAD.EQ.1.AND.IBEQUIL.NE.0)
     +        CALL SCRURITA1(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT)
C      IF(ISPIN.EQ.1.AND.IRAD.EQ.1.AND.IBEQUIL.NE.0)
C     +      CALL SCRURITA2(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT)
C      IF(ISPIN.EQ.1.AND.IRAD.EQ.1.AND.IDEPNON.NE.0)
C     +        CALL SCRURITA3(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT)
C
C
C
C=====END OF ENERGY LOOP
   50 CONTINUE

C54321 CONTINUE



      IF(IRAD.EQ.1.OR.ISPIN.EQ.0)THEN    ! Jan 2006: I don't understand this bit now.
      WRITE(53,103)
  103 FORMAT(/,/,/,/,/,'  ')
      IF(NCRAZY.NE.0)WRITE(53, '(A)')' Crazy normalisations in SPIN'
      STOP
      ENDIF



C
C
C
C
C
C=====CAN ONLY USE DISPERSION VERSION IN SPIN IF PERTOP IS CALLED FIRST.
C=====BUT FOR IPTCO = 2  PERTOP WAS NOT CALLED.
C=====TO SAVE TIME GET THE DISPERSION IN PERTOP FOR THE FIRST ENERGY AND HOPE
C=====THAT IT DOES NOT CHANGE IT MUCH WITH ENERGY DUE (SAY) TO TUNING A ROTATOR
C=====FOR EACH ENERGY.
C=====N.B. IN DERTOP AND DX88 WE ARE FREE TO CHOOSE THE DISPERSION.
C=====GET NORMAL DISPERSION SINCE NO CHROMATIC SHIFT HAS BEEN MADE YET.
      IF(MRDISP.NE.0.AND.IPTCO.GE.2)WRITE(53,1059)
      IF(MRDISP.NE.0.AND.IPTCO.GE.2.AND.ISECT.EQ.0)
     +                              CALL PERTOP(E0,IE0,U0,CIR,TM6A,TUN)
C
C
C
C   JAN 2006: WHY DOES THIS NEXT STUFF DEPEND ON MRDISP? KILL IT!
      IF(1.EQ.1)GO TO 761
      IF(MRDISP.EQ.0)GO TO 761
C=====IF RUNNING OFF ENERGY (ECHROM.NE.0) RE-SETUP SOME MATRICES BY
C=====REPEATING PART OF LATTIN FOR QUADS  + C.F.
      DO 72 I=1,NTYPE
      IID=ID(I)
      IF(IID.EQ.99)GO TO 72
      GO TO(72,72,107,108,72,72,72,72,72, 72,72, 72, 72, 72,118,119),IID
C            D  B   Q  QK RF HK VK SX BY SOL RQ FBH FBV ROT CFH CFV
C=====NORMAL QUAD
  107 CALL THIQUAD(XX(I),YY(I),IKICK,I,TMAT(1,1,I),ECHROM)
      GO TO 72
C=====SKEW QUAD
C     CALL UNITM(7,TMAT(1,1,I)):  ZERO IT INSIDE THE SKEWQD ITSELF.
  108 CALL SKEWQUAD(XX(I),YY(I),IKICK,I,TMAT(1,1,I),ECHROM)
      GO TO 72
C=====HORIZONTAL COMBINED FUNCTION DIPOLE
  118 CALL CFDIPH(XX(I),X2(I),YY(I),TMAT(1,1,I),ECHROM)
      GO TO 72
C=====VERTICAL COMBINED FUNCTION DIPOLE
  119 CALL CFDIPV(XX(I),X2(I),YY(I),TMAT(1,1,I),ECHROM)
      GO TO 72
   72 CONTINUE
  761 CONTINUE
C
C
C
C
C=====IF IRAD=0 CALCULATE THE SPIN STUFF IN A LOOP OVER ENERGY HERE.
C=====IF IPTCO=2 SET THE C.O. TO ZERO FOR THE SPIN STUFF.
C     IF IPTCO=0 THE C.O. WAS NEVER CALCULATED AND IS THUS ZERO ANTWAY.
C=====IF IVROT=1 RESET THE ROTATORS & MATRICES AT EACH NEW ENERGY.
      IF(IPTCO.NE.2)GO TO 265
      DO 264 I=1,NELEM
      DX(I)=0.D0
      DY(I)=0.D0
      DXP(I)=0.D0
      DYP(I)=0.D0
      DL(I)=0.D0
      DEL(I)=0.D0
  264 CONTINUE
C
C
  265 CONTINUE

      IF(1.EQ.1)GO TO 778
C     IF(1.EQ.1)GO TO 777
C===NOW THE ENERGY LOOP JUST OVER SPIN STUFF
  778 DO 60 IE0=1,NSTEP
      E0=E00+(IE0-1)*DE0
      CALL CAVITY(E0,U0,RLENGE,CRAD)
      IF(IVROT.EQ.0)GO TO 61
      CALL CHROTN(E0,0,1)
      DO 62 I=1,NTYPE
      IID=ID(I)
      IF(IID.EQ.2.OR.IID.EQ.9)CALL CHROTN(E0,I,2)
      IF(IID.EQ.2)CALL DIPH(XX(I),YY(I),TMAT(1,1,I))
      IF(IID.EQ.9)CALL DIPV(XX(I),YY(I),TMAT(1,1,I))
   62 CONTINUE
   61 CONTINUE
C
C     IF(ISECT.EQ.0)CALL SPIN(IE0,E0,CIR,TAUY,BETAY0,TUN)
C     IF(ISECT.EQ.0)CALL SODOM(IE0,E0,CIR,TUN)
      IF(ISECT.EQ.0)CALL SPIN(IE0,E0,CIR,TAUY,BETAY0,TUN,NCRAZY)
C      IF(1.EQ.1)STOP

      IF(IBEQUIL.NE.0.AND.ISECT.EQ.0)
     +        CALL SCRURITA1(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT)
      IF(IDEPLIN.NE.0.AND.ISECT.EQ.0)
     +        CALL SCRURITA2(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT)
      IF(IDEPNON.NE.0.AND.ISECT.EQ.0)
     +        CALL SCRURITA3(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT,
     +                                                 OFFSETL,OFFSETE)
      IF(IDEPNON.NE.0.AND.ISECT.EQ.1)
     +        CALL SCRURITA4(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT)

      IF(ISECT.GT.0)CALL SSECT(IE0,E0)
   60 CONTINUE
      GO TO 779
C
C
C===NOW THE SPIN TUNE LOOP JUST OVER SPIN STUFF
  777 DO 700 IE0=1,NSTEP
      E0=E00
      DO 702 I=1,NTYPE
      IID=ID(I)
      NSOLM6=NSOL(I)-6
      NNSOL=NSOL(I)
      XY=XX(I)
      IF(IID.EQ.10.AND.NSOLM6.GT.0.AND.XY.GT.0.D0)
     +   XX(I)=XX(I) + 2.D0*3.1415926359*DE0/(1.D0 + 0.00115965218073)
      IF(IID.EQ.10.AND.NSOLM6.GT.0.AND.XY.LT.0.D0)
     +   XX(I)=XX(I) - 2.D0*3.1415926359*DE0/(1.D0 + 0.00115965218073)
  702 CONTINUE
C
      IF(ISECT.EQ.0)CALL SPIN(IE0,E0,CIR,TAUY,BETAY0,TUN,NCRAZY)

      IF(IBEQUIL.NE.0.AND.ISECT.EQ.0)
     +        CALL SCRURITA1(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT)
      IF(IDEPLIN.NE.0.AND.ISECT.EQ.0)
     +        CALL SCRURITA2(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT)
      IF(IDEPNON.NE.0.AND.ISECT.EQ.0)
     +        CALL SCRURITA3(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT,
     +                                                 OFFSETL,OFFSETE)
      IF(IDEPNON.NE.0.AND.ISECT.EQ.1)
     +        CALL SCRURITA4(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT)

      IF(ISECT.GT.0)CALL SSECT(IE0,E0)
  700 CONTINUE


  779 WRITE(53,103)
      IF(NCRAZY.NE.0)WRITE(53, '(A)')
     +            ' WARNING: At least one crazy normalisation in SPIN'
      WRITE(53,103)



      CALL DATET



      CLOSE(68)
      CLOSE(67)
      CLOSE(66)
      CLOSE(65)
      CLOSE(64)
      CLOSE(63)
      CLOSE(62)
      CLOSE(61)
      CLOSE(60)
      CLOSE(59)
      CLOSE(58)
      CLOSE(57)
      CLOSE(56)
      CLOSE(55)
      CLOSE(54)
      CLOSE(53)
      CLOSE(52)
      CLOSE(51)
C
C
C
C
C
      STOP
C
C
 1059 FORMAT(' ',//,'IF USING DISPERSION FORMALISM  CALL PERTOP',
     +' IF IT WAS NOT ALREADY CALLED')
 1060 WRITE(53,1061) NREM
 1061 FORMAT(' ','RUN TERMINATED WITH ',I5,' SECONDS LEFT')
 1062 FORMAT(' ','TIME REMAINING:     ',I5,' SECONDS ')
      STOP
 1063 FORMAT(' ','MOTION IS ANTIDAMPED---SO STOP')
C
C
C
C
C
C
      END
