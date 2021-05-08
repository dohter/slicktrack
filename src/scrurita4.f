      SUBROUTINE SCRURITA4(IE0,E0,U0,CIR,COVMATIN,CRAD,TDAMP,SCALEMAT)
C
C
C
C
C   Routine to handle spin-orbit tracking in a linac and estimate the loss of 
C   polarisation (new version, was scrurita5.f)
C   with full 3-D spin motion derived from the extended G-matrix.
C   The acceleration means that a*gamma is changing from cavity to cavity.
C   We have the increase in energy in the cavities. 
C   In addition  the acceleration causes the damping of the beam sizes. Use MXDamp!
C
C   In contrast to the other scrurita's, no eigen-analysis is needed.
C   So this is a stripped down version of scrurita3.

C  

C
      IMPLICIT REAL*8(A-H,O-Z)
C 
      CHARACTER *1 TEXT(80)
      INCLUDE "cnlist.for"
      INCLUDE "clatic.for"
      INCLUDE "cloorb.for"
      INCLUDE "csol.for"
      INCLUDE "cemit.for"
      INCLUDE "cdisp.for"
      INCLUDE "cglist.for"
      INCLUDE "cmags.for"
      INCLUDE "cmontc.for"
      INCLUDE "cspindiff.for"
      INCLUDE "cpol.for"
     
      REAL*8   TOTALKICK(9),NEWNU, DELTA, E0,AAEXT(8)
      PARAMETER (NPART   = 10000)  ! Maximum allowed particles. Use NPART3 of them.
      PARAMETER (LIMSECT = 30000)  ! Maximum allowed sections.
      PARAMETER (MAXT    = 50000)  ! Maximum number of turns.
      PARAMETER (MAXE    = 2000)  ! Maximum number of energy steps: redundant
      REAL*8     NU, N0END(3)
      DIMENSION ZR3(3,3),ZI3(3,3),A(3,3),B(3,3)
      DIMENSION ROT(3,3),TM3A(3,3),ROTT(3,3)
      DIMENSION TM3B(3,3),WR3(3),WI3(3),ZW(3,3),ZWSAVE(3,3),TM3C(3,3)
      DIMENSION TRIN3(3,3),RR3(3),RI3(3),VR3(3,3),VI3(3,3),INTGE3(3)
      DIMENSION TW8A(8,8),TW2A(2,2)
      DIMENSION DISTMEAN(6),COVVEC(28),COVMATIN(6,6),COVMAT(6,6)
      DIMENSION SCALEMAT(6,6)
      DIMENSION COVSIM(8,8),COVTRACK(8,8,MAXT)
      DIMENSION COVSPINA(2,2,MAXT)
      DIMENSION SORBVEC(8,NPART),SORBVET(NPART,8) 
      DIMENSION SPINVECA(2,NPART)
C,JEFFVEC(9)
      DIMENSION SPINKICKA(3,NPART),SPINKICK1(3,NPART)
      DIMENSION ENVEC(2,NPART),REDSPIN(2,NPART),REDSPINT(NPART,2)
      DIMENSION SPINVETA(NPART,2)

      DIMENSION ONE(NPART),SORBAVE(8),SORBSIG(8)
      DIMENSION CC(8,8),AAA(8,8)
      DIMENSION CC9(9,9),AAA9(9,9),SECTMAP(9,9,LIMSECT),PSCALE(LIMSECT)
      DIMENSION IBBSECT(LIMSECT)
      DIMENSION XXBB(LIMSECT),X2BB(LIMSECT),YYBB(LIMSECT) 
      DIMENSION TREV8(8,8),TREV8D(8,8),TREV8WOWB(8,8)
      DIMENSION BIGPHOT(NPART),TDAMP(6)
      DIMENSION SDIFFUSION(MAXT),    DFIT(2),    DVEC(3),   DIFFMAT(2,2)
      DIMENSION SIGDFIT(2)
      DIMENSION SDIFFUSIONA(MAXT),   DFITA(2),   DVECA(3)
      DIMENSION SDIFFUSIONB(MAXT),   DFITB(2),   DVECB(3)
      DIMENSION SDIFFUSIONAB(MAXT),  DFITAB(2),  DVECAB(3)

      DIMENSION TAUDMC(MAXE),TAUDMC1(MAXE),TAUDMC2(MAXE),TAUDMC3(MAXE)
      DIMENSION TAUDMC13(MAXE),TAUDMC23(MAXE),TAUDMC12(MAXE)


      DIMENSION ODIFFUSION(6,MAXT)

      DIMENSION KGBINX(-500:500),KGBINY(-500:500),KBIN(100)
      DIMENSION NRB1(0:100),NRB2(0:100),NRB3(0:100)
      DIMENSION SORBDUM(8)
C=====8x8 eigenvector stuff.
      DIMENSION TRIN8(8,8),RR8(8),RI8(8),VR8(8,8),VI8(8,8)
      DIMENSION INTGE8(8),WW8(8,8)
      DIMENSION ZZ(8,8),ZV(8,8),WR8(8),WI8(8),ZZSAVE(8,8),USYMP(6,6)
      DIMENSION TM8A(8,8),TM8B(8,8)
      DIMENSION TN(8),PTUN(6)
      DIMENSION AB(6)
      DIMENSION TEMP(6,NPART),ORBAMP(6,NPART),GMAT(2,6),SIXVEC(6,NPART)
      DIMENSION GMAT9(3,6)
      DIMENSION TESTVEC(6)
      DIMENSION ZZSECT(8,8,LIMSECT),ZZCONJ(6,6)
      DIMENSION ZZWOWBCONJ(6,6)
      DIMENSION RA(0:3),RB(0:3),RC(0:3),TWQ(0:3)
      DIMENSION RQA(0:3,NPART),RQ1(0:3,NPART)
C
C===Set initial unit quarternion.
      RQA(0,:) =1.D0; RQA(1,:) =0.D0; RQA(2,:) =0.D0; RQA(3,:)  = 0.D0 

      PI=3.1415926535897932D0
      PI2=2.D0*PI
      NPART3=1000
     
C====Sync. excitation constant from Compton wavelength and electron radius
C    See physics 9709025 by Klaus Heinemann. 
      SYNCON = 2.818E-15 * 3.86E-13 * 55.D0/24.D0/DSQRT(3.D0)
      SYNCON = SYNCON * (E0/0.000511D0)**5

C====Naive spin tune.
      NU=E0/0.440652D0     ! Can scale a gamma indep. of the energy.

C===22/02/07 New lines for including acceleration and energy gain in cavities.
      NEWNU=NU
      E00=E0
      DELTAE=32.634*0.001D0
 202  FORMAT (ES16.9,ES16.9,ES16.9) 

      WRITE(53,103)
      WRITE(53,103)
  103 FORMAT(/,'  ')


C     MODES = 1                 !Switch on/off single and combined mode calcs.
      WRITE(53,929)MODES,IE0
  929 FORMAT('1','Entering subroutine SCRURITA4 for M-C tracking to ``me
     +asure'' the amount of depolarisation.', '   MODES = ',I2,
     +                                   F10.4)



      IF(NPART3.GT.NPART.OR.NSTEP.GT.MAXE)THEN
      WRITE(53,103)
      WRITE(53,'(3A)')' STOP: attempt to work with more than NPART',
     +                 '  particles or more than MAXE energy steps.'
      STOP
      ENDIF
C
C   *************************************************************
C   * SPIN ROTATION MATRIX AND THE ORTHONORMAL SPIN BASE VECTORS*
C   * AT THE FIRST BEAM-LINE ELEMENT                            *
C   *************************************************************

C
C=====CALCULATE THE SPIN TRANSPORT MATRIX=================================
      CALL UNIT(3,ROT)
      ANGH=0.D0
      ANGV=0.D0
      DO 231 II=1,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      IF(IID.EQ.1)GO TO 231
      XY=XX(ITY)
      IF(IID.EQ.2.OR.IID.EQ.15)ANGH=ANGH+XY
      IF(IID.EQ.9.OR.IID.EQ.16)ANGV=ANGV+XY 
      XX2=X2(ITY)
      YYY=YY(ITY)

      IF(IID.EQ.5.AND.NAME(ITY)(1:2).EQ.'ML')THEN
      XY=0.D0
      XX2=X2(ITY)
      YYY=DELTAE/E0
      ENDIF

      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NEWNU,TM3A)

      CALL JAM333(ROT,TM3A,ROT)
C------- In case of RF cavity--- energy update
      IF(IID.EQ.5.AND.NAME(ITY)(1:2).EQ.'ML')THEN
      E0=E0+DELTAE
      NEWNU=E0/0.440652D0
      ENDIF
     
  231 CONTINUE

C=====CHECK TOTAL ORBIT DEFLECTION.
C=====IS O.K.
      DANGH=DABS(PI2-ANGH)
      WRITE(53,2312)ANGH,ANGV
 2312 FORMAT(' ','Total bending angles', 2F15.10)
 2311 CONTINUE
C
C=====1)GET SPIN TUNE FROM TRACE OF MATRIX. 2)EFFECT OF ARGUS SOLENOID ON TUNE.
      STUNE=(ROT(1,1)+ROT(2,2)+ROT(3,3)-1.D0)*0.5D0
      STUNE=DACOS(STUNE)/(2.D0*PI)
C=====3)CHECK ORTHOGONALITY OF THE ROTATION MATRIX
      CALL ORTCHK(ROT)
C
C
C
C===Get spin basis at start by inverting from a longitudinal n_0 at the end.

      ROTT = TRANSPOSE(ROT)

      N0END(1)= 0.0D0
      N0END(2)= 0.0D0
      N0END(3)= 1.0D0

      ZW(:,1) =MATMUL(ROTT,N0END)

      ZW(1,2)=0.D0
      ZW(2,2)=1.D0
      ZW(3,2)=0.D0
      ZW(1,3)=ZW(2,1)*ZW(3,2)-ZW(3,1)*ZW(2,2)
      ZW(2,3)=ZW(3,1)*ZW(1,2)-ZW(1,1)*ZW(3,2)
      ZW(3,3)=ZW(1,1)*ZW(2,2)-ZW(2,1)*ZW(1,2)

      ZWSAVE = ZW
C

      WRITE(53,933)NU,STUNE,RNU
      DO 250 I=1,3
  250 WRITE(53,932)(ROT(I,J),J=1,3),WR3(I),WI3(I),(ZW(J,I),J=1,3)

 933  FORMAT(' ',' E0/.440652=',T14,F20.8,3X,F9.6,3X,F9.6,
     +                               //,' SPIN ROT. MATRIX AROUN',
     +   'D',T50,'            ',T79,'CHOSEN SPIN BASIS:',/,
     +   T5,'THE 1-ST BEAM-LINE ELEMENT:',T52,'REAL    IMAG')
 932  FORMAT(3F10.5,T50,2F9.5,T71,'--->  (',T79,3F9.5,T107,')')

      WRITE(53,98)
C
C
C-------------------------------------------------------------------

      COVMATIN      = 0.0D0


C============March07: read data from input file instead======
C============ READ the data for covariance matrix elements========= 

      READ(5,'(80A1)')TEXT 
      WRITE(53,'(80A1)')TEXT
C Externally supplied initial beam parameters:'
      READ(5,'(30A1,F12.8)') (TEXT(KT),KT=1,30), AAEXT(1)
      WRITE(53,'(A,30A1,F12.8)')' ',(TEXT(KT),KT=1,30),AAEXT(1)
      READ(5,'(30A1,F12.8)') (TEXT(KT),KT=1,30),AAEXT(2)
      WRITE(53,'(A,30A1,F12.8)')' ',(TEXT(KT),KT=1,30),AAEXT(2)
      READ(5,'(30A1,F12.8)') (TEXT(KT),KT=1,30),AAEXT(3)
      WRITE(53,'(A,30A1,F12.8)')' ',(TEXT(KT),KT=1,30),AAEXT(3)
      READ(5,'(30A1,F12.8)') (TEXT(KT),KT=1,30),AAEXT(4)
      WRITE(53,'(A,30A1,F12.8)')' ',(TEXT(KT),KT=1,30),AAEXT(4)
      READ(5,'(30A1,F12.8)') (TEXT(KT),KT=1,30),AAEXT(5)
      WRITE(53,'(A,30A1,F12.8)')' ',(TEXT(KT),KT=1,30),AAEXT(5)
      READ(5,'(30A1,F12.8)') (TEXT(KT),KT=1,30),AAEXT(6)
      WRITE(53,'(A,30A1,F12.8)')' ',(TEXT(KT),KT=1,30),AAEXT(6)
      READ(5,'(30A1,F12.8)') (TEXT(KT),KT=1,30),AAEXT(7)
      WRITE(53,'(A,30A1,F12.8)')' ',(TEXT(KT),KT=1,30),AAEXT(7)
      READ(5,'(30A1,F12.8)') (TEXT(KT),KT=1,30),AAEXT(8)
      WRITE(53,'(A,30A1,F12.8)')' ',(TEXT(KT),KT=1,30),AAEXT(8)

C====Get real emittances 
      COVMATIN(1,1)=AAEXT(1)**2 
      COVMATIN(2,2)=AAEXT(2)**2
      COVMATIN(3,3)=AAEXT(3)**2 
      COVMATIN(4,4)=AAEXT(4)**2 
      COVMATIN(5,5)=AAEXT(5)**2
      COVMATIN(6,6)=AAEXT(6)**2
      COVMATIN(1,2) = AAEXT(7)*DSQRT(COVMATIN(1,1)*COVMATIN(2,2))
      COVMATIN(2,1) = COVMATIN(1,2)
      COVMATIN(3,4) =AAEXT(8)*DSQRT(COVMATIN(3,3)*COVMATIN(4,4))
      COVMATIN(4,3) = COVMATIN(3,4)


C======Set up a 6-D Gaussian distribution according to the covariance matrix
C      from emitnc.f. Use NPART3 particles.  

 1111   CONTINUE
        WRITE(53,98)


        WRITE(53,951)COVMATIN
   98 FORMAT(//,' The theoretical beam covariance matrix <Xi*Xj>'
     + ,'(mm*mm,mm*mrad...) from Linac_covmat:')
  951 FORMAT(T6,6F20.12)


C      IF(IE0.EQ.1)
C      CALL G05CBF(KSEED15)       !Seed for initial s-o distribution. Original = 15 
                                    
      DISTMEAN = 0.D0
      EPS = 1.D-12
      IFAIL1 = 0
C      CALL G05EAF(DISTMEAN,6,COVMATIN,6,EPS,COVVEC,28,IFAIL1)
      IF(IFAIL1.NE.0)THEN
      WRITE(53,'(A,A,I10,A)')' ', 'IFAIL1 = ',IFAIL1,' STOP'  
      ENDIF

      WRITE(53,103)
      WRITE(53,103)


      SORBVEC   = 0.D0     ! Set to zero in case G05EZF is switched off.
      SPINVECA  = 0.D0


      IF(KSEED15.GT.0)THEN
      DO 1  I = 1, NPART3
      IFAIL2 = 0 
C      CALL G05EZF(SORBVEC(1:6,I),6,COVVEC,28,IFAIL2)
      IF(IFAIL2.NE.0)THEN
      WRITE(53,'(A,A,I10,A)')' ', 'IFAIL2 = ',IFAIL2,' STOP'  
      ENDIF
    1 CONTINUE
      ENDIF


C-----------------------------------------------------------------------
C======Get the 1-pass spin-(symplectic) orbit matrix again.
      IDAMPFLG = 0
      S=0.D0
      E0=E00
      NEWNU=NU

      CALL UNIT(8,TREV8)
      TM3C = ZWSAVE
      DO 8  II=2,NELEM
      ITY =ITYPE(II)
      IID =ID(ITY)
      XY  =XX(ITY)
      XX2 =X2(ITY)
      YYY =YY(ITY)
C===9/03/07-- New  parameters in  cavity!
      IF(IID.EQ.5.AND.NAME(ITY)(1:2).EQ.'ML')THEN
      XX(ITY)=0.0
      YY(ITY)=DELTAE/E0
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY) 
      ENDIF
      IF(IID.EQ.1.AND.XY.LT.0.00011D0)GO TO 8   !Skip zero length drifts.
      IF(IID.EQ.1)THEN                                !Treat drifts separately.
      TREV8(1,:) = TREV8(1,:) + XY * TREV8(2,:) 
      TREV8(3,:) = TREV8(3,:) + XY * TREV8(4,:) 
      GO TO 8
      ELSEIF(IID.EQ.10)THEN
      ZW = TM3C
      CALL SOL8AN(II,XY,YYY,NEWNU,ZW,TM3C,CC,NSOL(ITY)) 
      ELSE
      CC(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
C      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NEWNUNU,TM3A)
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NEWNU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 
      CALL MX88DAMP(CC,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NEWNU,TM3B,ZW,
     +                                          CRAD,NAME(ITY),IDAMPFLG)
     
      ENDIF
      TREV8 = MATMUL(CC,TREV8)
C ===== In case of cavity The energy update!
      IF(IID.EQ.5.AND.NAME(ITY)(1:2).EQ.'ML')THEN
      E0=E0+DELTAE
      NEWNU=E0/0.440652D0
      ENDIF
    8 CONTINUE

      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')
     +   '1-pass spin-(symplectic)orbit matrix:IDMPF=0,generated afresh'

      DO 36 I=1,8
   36 WRITE(53,914)(TREV8(I,J),J=1,8)
  914 FORMAT(8F12.5)

      CALL SYMP8(TREV8)

      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,934)
  934 FORMAT(' ',' Spin basis at the exit')
  935 FORMAT(' ',T79,3F9.5)
      DO 251 I=1,3
  251 WRITE(53,935)(ZW(J,I),J=1,3)


C

C===Set up the starting spin basis again==========================================

      WRITE(53,103)
      DO 66 INP = 1,NPART3
      IF(INP.LE.10)THEN
      WRITE(53,'(A)')
     + ' A sample of starting spin-orbit vectors (mm,mrad..)'
      WRITE(53,'(A,A,I5,8F12.6)')' ', 
     +     'Spin-orbit vector:    all',
     +                          INP,SORBVEC(1:6,INP),SPINVECA(:,INP)
      ENDIF
   66 CONTINUE  

C----------------------------------------------------------------------
C=====Get the covariance matrices for this sample with clever use of MATMUL
      SORBVET(1:NPART3,:) = TRANSPOSE(SORBVEC(:,1:NPART3))
      COVSIM = MATMUL(SORBVEC(:,1:NPART3),SORBVET(1:NPART3,:))/NPART3
      WRITE(53,899)

  899 FORMAT(//,' The generated spin-orbit covariance matrix',
     +     ' <Xi*Xj>(mm*mm,mm*mrad...) at the beginning of linac:')
      WRITE(53,955)COVSIM
  955 FORMAT(T6,8F16.8)

  923 FORMAT(2F12.5,2X,'--->  (',8F10.5,' )',F15.8)
  105 FORMAT(' ',2I10, 6E16.6)

C====Now set all particles on the C.O. with spins along n_0
C    so that a ``natural'' spin distribution can develop even if there
C    are strong orbital nonlinearities.
C    This also prevents potentially large oscillations of the two
C    SPINVEC covariances (although their sum doesn't oscillate even 
C    when the oscillations of the two parts are large).

C      SORBVEC  = SORBVEC  *0.D0
      SPINVECA = SPINVECA *0.D0

C===Put some starting tilt on the whole ensemble.
C      SPINVECA(1,:) = 100.D0      ! 100 mrad of initial tilt wrt n_0


C----------------------------------------------------------------------
C======Get the 1-pass spin-(damped) orbit matrix.

      IDAMPFLG = 1
      S=0.
      CALL UNIT(8,TREV8D)
      TM3C = ZWSAVE
      E0=E00
      NEWNU=NU
      DO 888  II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
C====9/03/07 New Parameters in cavities!
      IF(IID.EQ.5.AND.NAME(ITY)(1:2).EQ.'ML')THEN
      XX(ITY)=0.0
      YY(ITY)=DELTAE/E0
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
      ENDIF
      IF(IID.EQ.1.AND.XY.LT.0.00011D0)GO TO 888 !Skip zero length drifts 
      IF(IID.EQ.1)THEN
      TREV8D(1,:) = TREV8D(1,:) + XY * TREV8D(2,:) 
      TREV8D(3,:) = TREV8D(3,:) + XY * TREV8D(4,:) 
      GO TO 888

      ELSEIF(IID.EQ.10)THEN
      ZW = TM3C
      CALL SOL8AN(II,XY,YYY,NEWNU,ZW,TM3C,CC,NSOL(ITY))
      ELSE
      CC(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NEWNU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 
      CALL MX88DAMP(CC,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NEWNU,TM3B,ZW,
     +                                         CRAD,NAME(ITY),IDAMPFLG)
      ENDIF
C 26/02/07===== In case of cavity The energy update!
      IF(IID.EQ.5.AND.NAME(ITY)(1:2).EQ.'ML')THEN
      E0=E0+DELTAE
      NEWNU=E0/0.440652D0
      ENDIF   
      TREV8D = MATMUL(CC,TREV8D)

  888 CONTINUE
     
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,/,A,I1,A)')'1-pass spin-orbit (damped) matrix',
     +                               ' NDAMP3 = ',NDAMP3, ' !!!!!!!'
      DO 336 I=1,8
  336 WRITE(53,914)(TREV8D(I,J),J=1,8)



C----------------------------------------------------------------------------------
C      NOW 9x9 
C======Make 1 pass along the linac with  acceleration  to generate damped
C      9x9 spin-orbit matrices for sections between the dipole and cavity centres using
C      the fact that the C.O. is fixed.  
C      So don't have to recalc all fixed stuff for each particle and can track 
C      using the sections.
C      Get the radiation scale factors which account for the various dipole  
C      strengths.
C     
C      
      NSECT = 0                                  !Count sections  
      JSECT = 0                                  !Count non-dipole sections 
C===13/03/07     
      IDAMPFLG = 1                               !Switch on/off radn.
      CALL UNIT(9,AAA9)
      IBBSECT = 0 
      TM3C = ZWSAVE
      RADTOT = 0.D0

      E0=E00
      NEWNU=NU
      DO 7  II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)

C===26/02/07-- New  parameters in  cavity!
      IF(IID.EQ.5.AND.NAME(ITY)(1:2).EQ.'ML')THEN 
      XX(ITY)=0.0
      YY(ITY)=DELTAE/E0
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY) 
      ENDIF
C==========
      IF(IID.EQ.10)THEN
      ZW = TM3C
      CALL SOL9AN(II,XY,YYY,NU,ZW,TM3C,CC9,NSOL(ITY))
      ELSE
      CC9(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NEWNU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 

      CALL MX99DAMP(CC9,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NEWNU,TM3B,
     + ZW,CRAD,NAME(ITY),IDAMPFLG)
   
      ENDIF

      AAA9 = MATMUL(CC9,AAA9)
C=======================================
      IF((IID.EQ.7.AND.NAME(ITY)(1:2).EQ.'VD').OR.(II.EQ.NELEM).
     +   OR.(IID.EQ.5.AND.NAME(ITY)(1:2).EQ.'ML').            
     + OR.(II.NE.NELEM.AND.ID(ITYPE(II+1)).EQ.17))THEN  
      IF(II.EQ.NELEM)JSECT = JSECT + 1
      NSECT = NSECT + 1
     
      IBBSECT(NSECT) = 0
      IF(NSECT.GT.LIMSECT)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')' STOP: the number of sections exceeds LIMSECT'
      STOP
      ENDIF
     
C 26/02/07===== In case of cavity:  Energy update!
      IF(IID.EQ.5.AND.NAME(ITY)(1:2).EQ.'ML')THEN
      E0=E0+DELTAE
      NEWNU=E0/0.440652D0 
      ENDIF   
      SECTMAP(:,:,NSECT)   = AAA9

      PSCALE(NSECT) = RADSCALE(ITY)

C      IF(IE0.EQ.1) WRITE(54,'(A,4I10,2A,E16.5)')
C     +      ' ',II,IID,NSECT,JSECT,'  ',NAME(ITY),PSCALE(NSECT)

      RADTOT = RADTOT + PSCALE(NSECT) *SYNCON

      CALL UNIT(9,AAA9)      
      ENDIF

    7 CONTINUE 
      WRITE(53,'(A)') 'NSECT'
      WRITE(53,'(A)') ' Total Kick Angle in the lattice' 
      WRITE(53,9144) (TOTALKICK(IJJ),IJJ=1,4)
      WRITE(53,'(A,E16.5,A,E16.5)')' Total radiation excitation: ',
     +                                                       RADTOT


      RADAVE = RADTOT/(NSECT-JSECT)

C======Check the 1-pass matrix using the sections.
      CALL UNIT(9,AAA9)
      DO 3 II = 1,NSECT
      AAA9 = MATMUL(SECTMAP(:,:,II),AAA9)

 9149 FORMAT(9F12.5)
 9144 FORMAT(4F12.5)
    3 CONTINUE


      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,A,I6)')' ','NSECT', NSECT
      WRITE(53,103)
      WRITE(53,103)

      WRITE(53,103)
      WRITE(53,'(A,A,/,A,I1,A)')' 1-pass spin-(damped)orbit matrix',
     +           ' using sections.', ' NDAMP3 = ',NDAMP3, ' !!!!!!!'


      DO 35 I=1,9
   35 WRITE(53,9149)(AAA9(I,J),J=1,9)
C      WRITE(53,'(A)')' Element (9,6) should be about 2pi.a.gamma.'


C-------------------------------------------------------------------------------------
C======Now track. 
C======After each dipole or C-F, radiate ``big'' photons according to a 
C      centred top-hat distribution. 
C      In contrast to a Gaussian, there are no long tails so that no cut 
C      against too large (particle expelling kicks) is needed.
C      For linear motion this should also be enough to get Gaussian phase space
C      by the central limit theorem. 
C      Use the seed, KSEED(16).
C      Instead of the Boege formula, use the explicit formula for
C      the rms size of ``big photons'' in a dipople. The rms ``big'' photon 
C      energy is very much smaller than the energy spread of the beam.
C      This can be important for dipoles where n_0 is horizontal and there 
C      is no spin match.
C      No radiation in the correction coils. 
C      Set the damping  factor NDAMP3=1. 
C
C      For linear tracking, no need to scale the s-orbit vector to metres etc. 
C      But do it anyway to prepare for the future when nonlinear orbit stuff is
C      used. Then use metre,radian,seconds for all the remaining calcs but 
C      print in mm, mrad.,msecs.
C      NSECT includes the starting point (IP). But that give no radiation
C      So to get the size of big photons use NSECT-1.
C
C     
C
C
      SORBVEC   = SORBVEC*1.D-3
      SPINVECA  = SPINVECA*1.D-3



C   Set up unit quarternions to represent elements (7,8) of SORBVEC
C   i.e. the spin field at the starting particle coordinates. 
C   Be careful with correct assignment of axes and signs.
C   The SLICK dreibein n0,m,l can be read as l,n0,m and assigned the M.Vogt 
C   quarternion indices 1,2,3. 
C   The angle SPINVECA(1,:) (``alpha'' = kick along m) is a +ve rotation around l.
C   The angle SPINVECA(2,:) (``beta '' = kick along l) is a -ve rotation around m.
C   The angle SPINVECA(3,:) (no-name                 ) is a +ve rotation around n_0. 
C   SPINVECA(3,:) is not assigned at the start since the other 2 provide a good
C   enough starting representation of the spin field.  
C   See further comments below.
C   All modes together.


      DO 666 IQ = 1,NPART3
      RB(0) =  1.0D0
      RB(1) =  0.5D0*SPINVECA(1,IQ)
      RB(2) =  0.D0
      RB(3) = -0.5D0*SPINVECA(2,IQ)
C   
C   Renormalise as in SPRINT. RB(2) is always zero here. So ignore it.
      RBNORM = DSQRT(RB(0)*RB(0) + RB(1)*RB(1) + RB(3)*RB(3))  
      RB(0)  = RB(0)/RBNORM
      RB(1)  = RB(1)/RBNORM
      RB(3)  = RB(3)/RBNORM
C   Now get the quarternion for the spin field vectors starting with the initial
C   RQ which represents n0.
      CALL QUARTM(RQA(:,IQ),RB,RQA(:,IQ))
C


  666 CONTINUE 

     

      COVMAT  = COVMATIN *1.D-6
      SCALEMAT= COVMAT



      WRITE(53,103)  
      WRITE(53,103)
 

      NDAMPTURNS = 1


      IF(NTURN3*NDAMPTURNS.GT.MAXT)THEN
      WRITE(53,103)
      WRITE(53,'(3A)')' STOP: attempt to work with more than MAXT turns'
      STOP
      ENDIF

      IDAMPFLG = 1                              !Switch on/off radn.
      TH = DSQRT(3.D0)                          !Gives a top hat of unit variance. 
      VARTH = TH*TH/3.D0                        !The variance for a top hat distn.


      WRITE(53,103)
      WRITE(53,103)

      WRITE(53,'(A,F15.7,/,A,F15.7,/,A,F15.7)')
     +  ' Variance of generator             ', VARTH

      WRITE(53,103)


C      IF(IE0.EQ.1)
C      CALL G05CBF(KSEED16)             !Set the seed for the big photons.

C===Special action if NTURN3 = 1
      IF(NTURN3.EQ.1)NDAMPTURNS = 1


      SPK1MAX = 0.D0
      SPK2MAX = 0.D0
      SPK3MAX = 0.D0
      SIGK1 = 0.D0
      SIGK2 = 0.D0
      SIGK3 = 0.D0
      AVEK1 = 0.D0
      AVEK2 = 0.D0
      AVEK3 = 0.D0
      NRB1  = 0
      NRB2  = 0
      NRB3  = 0
      NSIGK = 0

      IDNT = 0
      NSPINPLOT = 0 
      KBIN     = 0
      SPINAVE1 = 0.D0
      SPINSIG  = 0.D0
      NSPIN    = 0
      ANGBIN   = 0.01D0

      NTND = 1
      DO 9 INT = 1,1              !Loop over damping periods.
      DO 4 IDT = 1,1              !Loop over turns in a damping period. 

      IDNT = IDNT + 1
      IF(IDNT.GT.MAXT)THEN
      WRITE(53,'(A)')' Number of turns exceeds MAXT. SO STOP'
      STOP
      ENDIF

      DO 5 II  = 1,NSECT                !Loop over sections in a turn.
      GMAT9= SECTMAP(7:9,1:6,II)

      SPINKICKA(:,1:NPART3)  = MATMUL(GMAT9,SORBVEC(1:6,1:NPART3))

C     
C   Set up the unit quarternions and use them to transport spins. 
C   Be careful with correct assignment of axes and signs.
C   The SLICK dreibein n0,m,l can be read as l,n0,m and assigned the M.Vogt 
C   quarternion indices 1,2,3. 
C   The angle SPINKICK(1,:) (``alpha'') is a +ve rotation around l.
C   The angle SPINKICK(2,:) ( no name ) is a +ve rotation around n_0.
C   The angle SPINKICK(3,:) (``beta '') is a -ve rotation around m. 
C
C   Alternatively and more straightforwardly, begin with eq. 2.48
C   in M. Vogt's thesis and use the first order solution obtained 
C   simply by integrating the Omegas which, in this case, are the 
C   omega.l, omega.n_0 and omega.m . These integrals are given by GMAT9*ORBVEC
C   except that row 2 of GMAT9 has the SLIM/SLICK sign so that the sign of 
C   SPINKICK(3, ) must be reversed in order to get the correct sign for the 
C   3rd quarternion component. One could also directly reverse the sign in 
C   MX99DAMP. But that could cause confusion under comparison with MX88DAMP.  
C
C   Eq. 3.32 in M. Berglund's thesis also indicate the the omega.l and omega.n_0 
C   terms have the same sign but that the sign on the omega.m term should be 
C   reversed.


      DO 777 IQQ = 1,NPART3
      RB(0) =  1.0D0
      RB(1) =  0.5D0*SPINKICKA(1,IQQ)
Cmore test10/01/07
      RB(2) =  0.5D0*SPINKICKA(3,IQQ) *(1.D0) 
      RB(3) = -0.5D0*SPINKICKA(2,IQQ)

C    Check the size of kicks. Search for the largest on the last turn   
      IF(IE0.EQ.1.AND.INT.EQ.NTURN3.AND.IDT.EQ.NDAMPTURNS)THEN
      IF(DABS(RB(1)).GE.SPK1MAX)SPK1MAX=DABS(RB(1))
      IF(DABS(RB(2)).GE.SPK2MAX)SPK2MAX=DABS(RB(2))
      IF(DABS(RB(3)).GE.SPK3MAX)SPK3MAX=DABS(RB(3))
      SIGK1 = SIGK1 + DABS(RB(1))**2
      SIGK2 = SIGK2 + DABS(RB(2))**2
      SIGK3 = SIGK3 + DABS(RB(3))**2
      AVEK1 = AVEK1 + DABS(RB(1))
      AVEK2 = AVEK2 + DABS(RB(2))
      AVEK3 = AVEK3 + DABS(RB(3))
      NSIGK= NSIGK + 1

C====Bin the RB(:) in 1 mrad bins.
      KRB1 = DABS(RB(1))*1.D3  
      KRB2 = DABS(RB(2))*1.D3  
      KRB3 = DABS(RB(3))*1.D3  
      NRB1(KRB1) = NRB1(KRB1) + 1
      NRB2(KRB2) = NRB2(KRB2) + 1
      NRB3(KRB3) = NRB3(KRB3) + 1

      ENDIF


      CALL QUARTM(RQA(:,IQQ),RB,RQA(:,IQQ))


  777 CONTINUE 



C    Update SORBVEC(1-6) only after it has been used for other things.
C  

      IF(IBBSECT(II).NE.1.OR.(IBBSECT(II).EQ.1.AND.MODIBMBM.LT.2))
     +          SORBVEC(1:6,1:NPART3) = 
     +                MATMUL(SECTMAP(1:6,1:6,II),SORBVEC(1:6,1:NPART3))

C=====Now radiate. But only at dipoles, not at the IP. 
C     Weight according to the calculated big photon size.
C     Also apply the rescale factor PHOTSCALE. This should be close to sqrt(NDAMP3).
C     Also apply the scale factor RADFAC to enhance the size of big photons 
C     artificially.
C     Apply the energy kick to the vector at the END of this section using the strength
C     of the radiator marking the end of this section. Then we are ready for the next
C     section.


C      IF(II.NE.NSECT)THEN
C      WRITE(53,'(A,I10,E16.5)')' PSCALE:  ',II,PSCALE(II)
      IF(PSCALE(II).NE.0.D0)THEN
C      CALL G05FAF(-TH,TH,NPART3,BIGPHOT)
      DO 55 IPHS = 1,NPART3
      
      SORBVEC(6,IPHS) =  SORBVEC(6,IPHS) + BIGPHOT(IPHS)
     +                              * DSQRT(PSCALE(II)*SYNCON)*RADFAC

   55 CONTINUE

      ENDIF

    5 CONTINUE





C===Get the covariance matrices and the 3-D spin info.
C   Before setting the rotations which describe the initial spins along n,
C   the starting spin vector was along n0. This is (0,1,0) in M.Vogt
C   coordinates. So the projections (beta and alpha) on l and m are R(1,2) and R(3,2)
C   We could also just use the n0 component directly but we would be burying 
C   information that could be useful for diagnostics. 
C   Now, at last, apply the renormalisation factors.



      DO 91 IQQQ = 1,NPART3
C   Renormalise as in SPRINT. 
      RQNORM = DSQRT(RQA(0,IQQQ)**2 + RQA(1,IQQQ)**2
     +              +RQA(2,IQQQ)**2 + RQA(3,IQQQ)**2)
      RQA(0,IQQQ)  = RQA(0,IQQQ)/RQNORM
      RQA(1,IQQQ)  = RQA(1,IQQQ)/RQNORM
      RQA(2,IQQQ)  = RQA(2,IQQQ)/RQNORM
      RQA(3,IQQQ)  = RQA(3,IQQQ)/RQNORM

C

      SPINVECA(2,IQQQ) 
     +      = 2.D0*(RQA(1,IQQQ)*RQA(2,IQQQ) - RQA(0,IQQQ)*RQA(3,IQQQ)) !beta along l
      SPINVECA(1,IQQQ) 
     +      = 2.D0*(RQA(3,IQQQ)*RQA(2,IQQQ) + RQA(0,IQQQ)*RQA(1,IQQQ)) !alpha along m

  91  CONTINUE 


C====Covariance matrices.
      SORBVET(1:NPART3,:)     = TRANSPOSE(SORBVEC(:,1:NPART3))
      SPINVETA(1:NPART3,:)    = TRANSPOSE(SPINVECA(:,1:NPART3)) 

      IF(IDNT.EQ.1)
     +        REDSPINT(1:NPART3,:)   = TRANSPOSE(REDSPIN(:,1:NPART3))
      COVTRACK(:,:,IDNT) = 0.D0   !Cleans out spin--which is not calculated.
      COVTRACK(1:6,1:6,IDNT) 
     +     = MATMUL(SORBVEC(1:6,1:NPART3),SORBVET(1:NPART3,1:6))/NPART3
      COVSPINA(:,:,IDNT)
     +     = MATMUL(SPINVECA(:,1:NPART3),SPINVETA(1:NPART3,:))/NPART3


      SDIFFUSION(IDNT)      = COVSPINA(1,1,IDNT)+COVSPINA(2,2,IDNT)
      SDIFFUSIONA(IDNT)     = COVSPINA(1,1,IDNT)
      SDIFFUSIONB(IDNT)     = COVSPINA(2,2,IDNT)
      SDIFFUSIONAB(IDNT )   = COVSPINA(1,2,IDNT)


      DO 90 IK = 1,6
      ODIFFUSION(IK,IDNT) = 1.D0
      IF(SCALEMAT(IK,IK).NE.0.D0)
     +       ODIFFUSION(IK,IDNT) = COVTRACK(IK,IK,IDNT)/SCALEMAT(IK,IK)
   90 CONTINUE




    4 CONTINUE 


C====At the end of the line write out sample info

      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,85)
   85 FORMAT(//,' The  transported s-o covariance matrix',
     +          ' <Xi*Xj>(mm*mm,mm*mrad...)',
     +          ' at the end of the linac ')
      WRITE(53,957)COVTRACK(:,:,IDNT)*1.D6
  957 FORMAT(T6,8F16.12)

      WRITE(53,86)
   86 FORMAT(/,' The  transported spin covariance matrix',
     +          ' mode-A <Si*Sj>(mrad*mrad...)',
     +          ' at the end of the linac')
      WRITE(53,956)COVSPINA(:,:,IDNT)*1.D6
  956 FORMAT(T6,2F16.3)


    9 CONTINUE


      IF(1.EQ.1)RETURN



C======Check the variance of the last big photons. 
      WRITE(53,103)
      WRITE(53,103)
      VARPHOT = 0
      IF(PHOTSCALE.GT.0.D0)THEN
      DO 16 IPH = 1,NPART3
   16 VARPHOT = VARPHOT + BIGPHOT(IPH)**2
      VARPHOT = VARPHOT/NPART3/PHOTSCALE**2
      ENDIF
      WRITE(53,'(A,A,F10.6)')' ','Variance of the top-hat samples ',
     +                           VARPHOT


C-------------------------------------------------------------------------
C      Print some statistical stuff. Rescale to mm. etc.

      ONE = 1.D0
      SORBAVE = MATMUL(SORBVEC(:,1:NPART3),ONE(1:NPART3))/NPART3
      WRITE(53,96)
      WRITE(53,95)SORBAVE*1.D3
   96 FORMAT(//,' The transported spin-beam average <Xi>(mm,mrad...)',
     +           ' at the 1-st beam line element:')
   95 FORMAT(T6,8F13.5)
 
      DO 94 IJ=1,8
   94 SORBSIG(IJ)=DSQRT(COVTRACK(IJ,IJ,IDNT))
      WRITE(53,93)
      WRITE(53,95)SORBSIG*1.D3
   93 FORMAT(//,' The transported rms spin-beam sizes(mm,mrad...)',
     +           ' at the 1-st beam line element:')


C======Plot histograms of projections: 
C      XYBIN  = 0.01D0
C      DO 101 IO = 1,NPART3
C      NGAUSX = SORBVEC(1,IO)*1.D3/XYBIN
C      KGBINX(NGAUSX) = KGBINX(NGAUSX) + 1
C      NGAUSY = SORBVEC(3,IO)*1.D3/XYBIN
C      KGBINY(NGAUSY) = KGBINY(NGAUSY) + 1
C  101 CONTINUE 
C      DO 102 JO = -500,500
C  102 WRITE(59,'(1X,I10,1X,I10,1X,I10)')JO,KGBINX(JO),KGBINY(JO)

C======Plot histogram of spin angles. 
      DO 102 JS = 1,100
  102 WRITE(59,'(1X,I10,1X,I10,1X,I10)')JS,KBIN(JS)
      SPINAVE      =  SPINAVE/NSPIN
      SPINSIG      =  SPINSIG/NSPIN
      SPINSPREAD   =  DSQRT(SPINSIG - SPINAVE**2)
      SPINSPREAD   =  SPINSPREAD * 180.D0/pi
      SPINAVE      =  SPINAVE * 180.D0/pi
      DECOHANG     =  NU * DSQRT(SCALEMAT(6,6))/DABS(TN(5)) * 180.D0/pi
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(2A,I6,2F10.4,A)')
     + ' The number of spins, the spin spread in the alpha-beta plane',
     + ' and a model prediction: ',
     + NSPIN,SPINSPREAD,DECOHANG,' degrees'
C----------------------------------------------------------------------------------------
C======Plot histograms of horizontal and vertical projections: 
      IF(IBEQUIL.EQ.0)THEN
      XYBIN  = 0.01D0
      DO 101 IO = 1,NPART3
      NGAUSX = SORBVEC(1,IO)*1.D3/XYBIN
      KGBINX(NGAUSX) = KGBINX(NGAUSX) + 1
      NGAUSY = SORBVEC(3,IO)*1.D3/XYBIN
      KGBINY(NGAUSY) = KGBINY(NGAUSY) + 1
  101 CONTINUE 
C      DO 102 JO = -500,500
C  102 WRITE(59,'(1X,I10,1X,I10,1X,I10)')JO,KGBINX(JO),KGBINY(JO)
C  102 CONTINUE    
      ENDIF
C======Plot scatter of the spin angles alpha and beta at first energy point.
      IF(IE0.EQ.1)THEN
      DO 104 IOO = 1, NPART3 
C  104 WRITE(60,105)IOO,SORBVEC(7,IOO),SORBVEC(8,IOO)
  104 CONTINUE
      ENDIF



      RETURN

 9999 WRITE(53,92)
   92 FORMAT(' ERROR IN EIGEN VALUE ROUTINE')
      STOP
      END
