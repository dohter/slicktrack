      SUBROUTINE SCRURITA2(IE0,E0,U0,CIR,COVMATIN,CRAD,TDAMP,SCALEMAT)
C
C
C
C
C   Routine to handle spin-orbit tracking and estimate the rate of depolarisaiton
C   ------------------------------------------------------------------------------
C   with linearised spin motion as in SPIN.
C

C   NOTE!!!!!!!!!!!!!!!!!!!! -- the covariance matrix and the damping constants 
C   are usually (for convenience) those derived at the starting energy.
C   So, the details of the M-C results can change w.r.t. a run starting at low  
C   energy, if the run is restarted at a higher energy within the scan range 
C   of the earlier run. 
C
C

C
      IMPLICIT REAL*8(A-H,O-Z)
C
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
C 
      PARAMETER (NPART   = 10000)   ! Maximum allowed particles. Use NPART2 of them.
      PARAMETER (LIMSECT = 3000)    ! Maximum allowed sections.
      PARAMETER (MAXT    = 150000)   ! Maximum number of turns.
      PARAMETER (MAXE    = 2000)    ! Maximum number of energy steps.
      PARAMETER (NFFT    = 10000)   ! Maximum number of turns for the FFT of the orbit.
      REAL*8     NU
      DIMENSION ZR3(3,3),ZI3(3,3),A(3,3),B(3,3)
      DIMENSION ROT(3,3),TM3A(3,3)
      DIMENSION TM3B(3,3),WR3(3),WI3(3),ZW(3,3),ZWSAVE(3,3),TM3C(3,3)
      DIMENSION TRIN3(3,3),RR3(3),RI3(3),VR3(3,3),VI3(3,3),INTGE3(3)
      DIMENSION TW8A(8,8),TW2A(2,2)
      DIMENSION DISTMEAN(6),COVVEC(28),COVMATIN(6,6),COVMAT(6,6)
      DIMENSION SCALEMAT(6,6)
      DIMENSION COVSIM(8,8),COVTRACK(8,8,MAXT)
      DIMENSION COVSPIN1(2,2,MAXT),COVSPIN2(2,2,MAXT),COVSPIN3(2,2,MAXT)
      DIMENSION COVREDSPIN(2,2,MAXT)
      DIMENSION SORBVEC(8,NPART),SORBVET(NPART,8)
      DIMENSION ORBVEC1(6,NPART),ORBVEC2(6,NPART),ORBVEC3(6,NPART)
      DIMENSION SPINVEC1(2,NPART),SPINVEC2(2,NPART),SPINVEC3(2,NPART)
      DIMENSION ENVEC1(2,NPART),ENVEC2(2,NPART),ENVEC3(2,NPART)
      DIMENSION ENVEC(2,NPART),REDSPIN(2,NPART),REDSPINT(NPART,2)
      DIMENSION SPINVET1(NPART,2),SPINVET2(NPART,2),SPINVET3(NPART,2)
      DIMENSION ONE(NPART),SORBAVE(8),SORBSIG(8)
      DIMENSION CC(8,8),AAA(8,8),SECTMAP(8,8,LIMSECT)
      DIMENSION IBBSECT(LIMSECT),PSCALE(LIMSECT)
      DIMENSION ZWBB(3,3,LIMSECT),NSOLBB(LIMSECT)
      DIMENSION XXBB(LIMSECT),X2BB(LIMSECT),YYBB(LIMSECT) 
      DIMENSION TREV8(8,8),TREV8D(8,8),TREV8WOWB(8,8)
      DIMENSION BIGPHOT(NPART),TDAMP(6)
      DIMENSION SDIFFUSION(MAXT),    DFIT(2),    DVEC(3),   DIFFMAT(2,2)
      DIMENSION SIGDFIT(2)
      DIMENSION SDIFFUSIONA(MAXT),   DFITA(2),   DVECA(3)
      DIMENSION SDIFFUSIONB(MAXT),   DFITB(2),   DVECB(3)
      DIMENSION SDIFFUSIONAB(MAXT),  DFITAB(2),  DVECAB(3)

      DIMENSION SDIFFUSION1(MAXT),   DFIT1(2),   DVEC1(3)
      DIMENSION SIGDFIT1(2)
      DIMENSION SDIFFUSIONA1(MAXT),  DFITA1(2),  DVECA1(3)
      DIMENSION SDIFFUSIONB1(MAXT),  DFITB1(2),  DVECB1(3)
      DIMENSION SDIFFUSIONAB1(MAXT), DFITAB1(2), DVECAB1(3)

      DIMENSION SDIFFUSION2(MAXT),   DFIT2(2),   DVEC2(3)
      DIMENSION SIGDFIT2(2)
      DIMENSION SDIFFUSIONA2(MAXT),  DFITA2(2),  DVECA2(3)
      DIMENSION SDIFFUSIONB2(MAXT),  DFITB2(2),  DVECB2(3)
      DIMENSION SDIFFUSIONAB2(MAXT), DFITAB2(2), DVECAB2(3)

      DIMENSION SDIFFUSION3(MAXT),   DFIT3(2),   DVEC3(3)
      DIMENSION SIGDFIT3(2)
      DIMENSION SDIFFUSIONA3(MAXT),  DFITA3(2),  DVECA3(3)
      DIMENSION SDIFFUSIONB3(MAXT),  DFITB3(2),  DVECB3(3)
      DIMENSION SDIFFUSIONAB3(MAXT), DFITAB3(2), DVECAB3(3)

      DIMENSION TAUDMC(MAXE),TAUDMC1(MAXE),TAUDMC2(MAXE),TAUDMC3(MAXE)

      DIMENSION ODIFFUSION(6,MAXT)

      DIMENSION KGBINX(-500:500),KGBINY(-500:500)
      DIMENSION SORBDUM(8)
C=====8x8 eigenvector stuff.
      DIMENSION TRIN8(8,8),RR8(8),RI8(8),VR8(8,8),VI8(8,8)
C      DIMENSION INTGE8(8),WW8(8,8)
      PARAMETER( LWORK=64*8 )   ! now needed for  F02EBF: MV/09.10.2007
      DIMENSION WORK(LWORK)     ! now needed for  F02EBF: MV/09.10.2007
      DIMENSION WW8(8,8)
      DIMENSION ZZ(8,8),ZV(8,8),WR8(8),WI8(8),ZZSAVE(8,8),USYMP(6,6)
      DIMENSION ZZWOWB(8,8)
      DIMENSION TM8A(8,8),TM8B(8,8)
      DIMENSION TN(8),PTUN(6)
      DIMENSION AB(6)
      DIMENSION TEMP(6,NPART),ORBAMP(6,NPART),GMAT(2,6),SIXVEC(6,NPART)
      DIMENSION TESTVEC(6)
      DIMENSION ZZSECT(8,8,LIMSECT),ZZCONJ(6,6)
      DIMENSION ZZWOWBCONJ(6,6)

      DIMENSION  XFFT(0:NFFT-1),WXFFT(NFFT)
      DIMENSION  YFFT(0:NFFT-1),WYFFT(NFFT)
      DIMENSION AXFFT(0:NFFT-1),BXFFT(0:NFFT-1)
      DIMENSION AYFFT(0:NFFT-1),BYFFT(0:NFFT-1)
      DIMENSION XXFFT(0:NFFT-1),YYFFT(0:NFFT-1)
C
C=====Set up the  SYMPLECTIC UNIT MATRIX: S
      USYMP     = 0.D0
      USYMP(1,2)=-1.D0
      USYMP(2,1)= 1.D0
      USYMP(3,4)=-1.D0
      USYMP(4,3)= 1.D0
      USYMP(5,6)=-1.D0
      USYMP(6,5)= 1.D0


      PI=3.1415926535897932D0
      PI2=2.D0*PI

C====Sync. excitation constant from Compton wavelength and electron radius
C    See physics 9709025 by Klaus Heinemann. 
      SYNCON = 2.818E-15 * 3.86E-13 * 55.D0/24.D0/DSQRT(3.D0)
      SYNCON = SYNCON * (E0/0.000511D0)**5

      NU=E0/0.440652D0 * 1.0D0       ! Can scale a gamma indep. of the energy.



      WRITE(53,103)
      WRITE(53,103)
  103 FORMAT(/,'  ')


      WRITE(53,929)MODES,IE0
  929 FORMAT('1','Entering subroutine SCRURITA2 for M-C tracking to ``me
     +asure'' the rate of depolarisation.', '   MODES = ',I2,
     +                                      '  Energy step =  ',I5)
C      write(*,929)

C      IF(IE0.EQ.1)NPART2 = 10 
C      IF(IE0.EQ.2)NPART2 = 1000


      IF(NPART2.GT.NPART.OR.NSTEP.GT.MAXE)THEN
      WRITE(53,103)
      WRITE(53,'(4A,I6,A,I6)')' STOP: attempt to work with more than',
     +           '  NPART particles or more than MAXE energy steps.',
     =                 'NPART2= ', NPART2,' NSTEP= ',NSTEP 
      STOP
      ENDIF
C
C   *************************************************************
C   * SPIN ROTATION MATRIX AND THE ORTHONORMAL SPIN BASE VECTORS*
C   * AT THE FIRST BEAM-LINE ELEMENT                            *
C   *************************************************************

C
C=====CALCULATE THE SPIN REVOL. MATRIX=================================
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
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      CALL JAM333(ROT,TM3A,ROT)
  231 CONTINUE
C
C=====CHECK TOTAL ORBIT DEFLECTION TO ENSURE THAT ROTATOR INTERPOLATION
C=====IS O.K.
      DANGH=DABS(PI2-ANGH)
      WRITE(53,2312)ANGH,ANGV
 2312 FORMAT(' ','Total bending angles', 2F15.10)
 2311 CONTINUE
C
C=====1)GET SPIN TUNE FROM TRACE OF MATRIX. 2)EFFECT OF ARGUS SOLENOID ON TUNE.
      STUNE=(ROT(1,1)+ROT(2,2)+ROT(3,3)-1.D0)*0.5D0
      STUNE=DACOS(STUNE)/(2.D0*PI)
      RNU=DACOS(DCOS(0.01D0)*DCOS(PI*NU))/PI-1.D0
C=====3)CHECK ORTHOGONALITY OF THE ROTATION MATRIX
      CALL ORTCHK(ROT)
C
C
C
C=====GET EIGENVECTORS & ORDER THEM=====================================
      TRIN3 = ROT
      IFAIL=0
C      CALL F02AGF(TRIN3,3,3,RR3,RI3,VR3,3,VI3,3,INTGE3,IFAIL)
      CALL F02EBF('V',3,TRIN3,3,RR3,RI3,VR3,3,VI3,3,WORK,LWORK,IFAIL) ! <<replmnt
      IF(IFAIL.NE.0)GO TO 9999
  926 FORMAT(' ','NOW THE F02AGF RESULTS: 3X3')

      WR3 = RR3      !real parts of eigenvalues
      WI3 = RI3      !imag parts of eigenvalues
      ZR3 = VR3      !real parts of eigenvectors
      ZI3 = VI3      !imag parts of eigenvectors 


      TUN=999999.
  922 FORMAT(2F12.5,2X,'--->  (',3F10.5,' )',F15.8)
      NN=4                                     !Set to a too big value.
      ICOUNT=0
      DO 255 I=1,3
      TM3A(I,1)=WR3(I)
      TM3A(I,2)=WI3(I)
C=====Locate the real unit eigenvalue:
  255 IF(DABS(WI3(I)).LT.1.D-9.AND.WR3(I).GT.0.9999D0)NN=I
      IF(NN .NE. 4)GO TO 241
      WRITE(53,242)
  242 FORMAT(' ','No real unit spin eigenvalue found----so STOP')
      STOP
C
  241 CONTINUE
      MM=MOD(NN,3)+1
      LL=MOD(MM,3)+1
      AZ1=DSQRT(ZR3(1,NN)**2+ZR3(2,NN)**2+ZR3(3,NN)**2)
      AZ2=DSQRT(ZR3(1,MM)**2+ZR3(2,MM)**2+ZR3(3,MM)**2)
C======We just need the normed real parts of each eigenvector. 
      DO 256 I=1,3
      ZW(I,1)=ZR3(I,NN)/AZ1
  256 ZW(I,2)=ZR3(I,MM)/AZ2
C      Force the m vector to be vertical.
C      ZW(1,2)=0.D0
C      ZW(2,2)=1.D0
C      ZW(3,2)=0.D0
      ZW(1,3)=ZW(2,1)*ZW(3,2)-ZW(3,1)*ZW(2,2)
      ZW(2,3)=ZW(3,1)*ZW(1,2)-ZW(1,1)*ZW(3,2)
      ZW(3,3)=ZW(1,1)*ZW(2,2)-ZW(2,1)*ZW(1,2)

C======Set signs of imaginary parts of the 2nd and 3rd eigenvectors.
      SGN=ZW(1,3)*ZI3(1,MM)+ZW(2,3)*ZI3(2,MM)+ZW(3,3)*ZI3(3,MM)
C=====Careful at half integer spin tune.
C     SGN=SGN/DABS(SGN)
      IF (SGN.GE.0.D0) SGN= 1.D0
      IF (SGN.LT.0.D0) SGN=-1.D0
      WR3(1)=TM3A(NN,1)
      WR3(2)=TM3A(MM,1)
      WR3(3)=TM3A(LL,1)
      WI3(1)=TM3A(NN,2)
      WI3(2)=TM3A(MM,2)*SGN
      WI3(3)=TM3A(LL,2)*SGN
C
      ZWSAVE = ZW
C

      WRITE(53,933)NU,STUNE,RNU
      DO 250 I=1,3
  250 WRITE(53,932)(ROT(I,J),J=1,3),WR3(I),WI3(I),(ZW(J,I),J=1,3)

  933 FORMAT(' ',' E0/.440652=',T14,F9.6,3X,F9.6,3X,F9.6,
     +                               //,' SPIN ROT. MATRIX AROUN',
     +   'D',T50,'EIGENVALUES:',T79,'SPIN BASIS FROM EIGENVECTORS:',/,
     +   T5,'THE 1-ST BEAM-LINE ELEMENT:',T52,'REAL    IMAG')
  932 FORMAT(3F10.5,T50,2F9.5,T71,'--->  (',T79,3F9.5,T107,')')
C
C
C----------------------------------------------------------------------------------
C======Set up a 6-D Gaussian distribution according to the covariance matrix
C      from emitnc.f. Use NPART2 particles.  

      WRITE(53,98)
C      write(*,98)

      WRITE(53,951)COVMATIN
   98 FORMAT(//,' The theoretical beam covariance matrix <Xi*Xj>'
     + ,'(mm*mm,mm*mrad...) at the 1-st beam line element:')
  951 FORMAT(T6,6F13.5)

C      IF(IE0.EQ.1)
C      CALL G05CBF(KSEED13)       !Seed for initial s-o distribution. Original = 15 
                                    
      DISTMEAN = 0.D0
      EPS = 1.D-12
      IFAIL1 = 0
C      CALL G05EAF(DISTMEAN,6,COVMATIN,6,EPS,COVVEC,28,IFAIL1)
      IF(IFAIL1.NE.0)THEN
      WRITE(53,'(A,A,I10,A)')' ', 'IFAIL1 = ',IFAIL1,' STOP'  
      ENDIF

      WRITE(53,103)
      WRITE(53,103)



      SORBVEC  = 0.D0     ! Set to zero in case G05EZF is switched off.
      SPINVEC1 = 0.D0
      SPINVEC2 = 0.D0      
      SPINVEC3 = 0.D0

      IF(KSEED13.GT.0)THEN
      DO 1  I = 1, NPART2
      IFAIL2 = 0 
C      CALL G05EZF(SORBVEC(1:6,I),6,COVVEC,28,IFAIL2)
      IF(IFAIL2.NE.0)THEN
      WRITE(53,'(A,A,I10,A)')' ', 'IFAIL2 = ',IFAIL2,' STOP'  
      ENDIF
    1 CONTINUE
      ENDIF

C======Plot histograms of projections: 
C      XYBIN  = 0.01D0
C      DO 101 IO = 1,NPART2
C      NGAUSX = SORBVEC(1,IO)*1.D3/XYBIN
C      KGBINX(NGAUSX) = KGBINX(NGAUSX) + 1
C      NGAUSY = SORBVEC(3,IO)*1.D3/XYBIN
C      KGBINY(NGAUSY) = KGBINY(NGAUSY) + 1
C  101 CONTINUE 
C      DO 102 JO = -500,500
C  102 WRITE(59,'(1X,I10,1X,I10,1X,I10)')JO,KGBINX(JO),KGBINY(JO)

C-----------------------------------------------------------------------
C======Get the 1-turn spin-(symplectic) orbit matrix again.
      IDAMPFLG = 0
      S=0.D0
      CALL UNIT(8,TREV8)
      TM3C = ZWSAVE
      DO 8  II=2,NELEM
      ITY =ITYPE(II)
      IID =ID(ITY)
      XY  =XX(ITY)
      XX2 =X2(ITY)
      YYY =YY(ITY)
      IF(IID.EQ.1.AND.XY.LT.0.00011D0)GO TO 8   !Skip zero length drifts.
      IF(IID.EQ.1)THEN                                !Treat drifts separately.
      TREV8(1,:) = TREV8(1,:) + XY * TREV8(2,:) 
      TREV8(3,:) = TREV8(3,:) + XY * TREV8(4,:) 
      GO TO 8
      ELSEIF(IID.EQ.10)THEN
      ZW = TM3C
      CALL SOL8AN(II,XY,YYY,NU,ZW,TM3C,CC,NSOL(ITY)) 
      ELSEIF(IID.EQ.17)THEN
      ZW = TM3C
      CC(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      CALL LBEAMBEAM(ITY,CC,ZW,NU)
      ELSE
      CC(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 
      CALL MX88DAMP(CC,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NU,TM3B,ZW,
     +                                         CRAD,NAME(ITY),IDAMPFLG)
      ENDIF
      TREV8 = MATMUL(CC,TREV8)
    8 CONTINUE

C===Store the 1-turn map without the wind-back.
      TREV8WOWB = TREV8 

C======WIND BACK THE SPIN BASIS==========================================
      CALL UNIT(8,TW8A)
      TW8A(7,7)= WR3(2)
      TW8A(7,8)= WI3(2)
      TW8A(8,7)=-WI3(2)
      TW8A(8,8)= WR3(2)

      TREV8    = MATMUL(TW8A,TREV8)

      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')
     +   ' 1-turn spin-(symplectic)orbit matrix, generated afresh.'

      DO 36 I=1,8
   36 WRITE(53,914)(TREV8(I,J),J=1,8)

C=====GET 8X8 EIGENVECTORS FOR ONE REVOLUTION & ORDER THEM==============
      TRIN8 = TREV8
      IFAIL=0
C      CALL F02AGF(TRIN8,8,8,RR8,RI8,VR8,8,VI8,8,INTGE8,IFAIL)
      CALL F02EBF('V',8,TRIN8,8,RR8,RI8,VR8,8,VI8,8,WORK,LWORK,IFAIL) ! <<replmnt
      IF(IFAIL.NE.0)GO TO 9999
C=====WRITE OUT F02AGF RESULTS
C     WRITE(53,9266)
 9266 FORMAT(' ','NOW THE F02AGF RESULTS: 8x8')
      DO 2822 I=1,4
      IT=2*I-1
      DO 2833 J=1,8
      ZZ(J,IT  )=VR8(J,IT)
 2833 ZZ(J,IT+1)=VI8(J,IT)
 2822 CONTINUE
      DO 2844 I=1,8
      WR8(I)=RR8(I)
 2844 WI8(I)=RI8(I)
      TUN=999999
      IF(DABS(RR8(IT)) .LE. 1.D0)TUN=DACOS(RR8(IT))/(2.*PI)
C
      DO 290 I=1,8
      SUM=0.
      DO 291 J=1,6
  291 SUM=SUM+ZZ(J,I)**2
      IF(SUM.LT.1.D-12)IS=I
      TM8B(I,1)=WR8(I)
  290 TM8B(I,2)=WI8(I)
C
      TM8A = ZZ
      DO 292 I=1,8
      IS=MOD(IS,8)+1
      WR8(I)=TM8B(IS,1)
      WI8(I)=TM8B(IS,2)
      DO 292 J=1,8
  292 ZZ(J,I)=TM8A(J,IS)
C
C
C=====Renormalize vectors as in SMILE but keep the usual NORM routine.
C=====Note: The next few lines redefine the tunes from RENORM.
      CALL RENORM(ZZ,8,WR8,WI8)
C
      WRITE(53,924)
  924 FORMAT(////,' Eigenvalues and eigenvectors:',/,T8,
     +                                          'REAL',
     +                                          T18,'IMAG',T122,'TUNES')
      DO 281 I=1,8
      TUN=999999.
C=====GET TUNES: THIS DOES NOT DIFFERENTIATE BETWEEN +/- ANGLES
C=====WITHOUT USING THE SINE.
C     IF(DABS(WR8(I)) .LE. 1.D0)TUN=DACOS(WR8(I))/(2.*PI)
      TUN=DATAN2(WI8(I),WR8(I))/(2.*PI)
      TN(I)=TUN
  281 WRITE(53,923)WR8(I),WI8(I),(ZZ(J,I),J=1,8),TUN
  923 FORMAT(2F12.5,2X,'--->  (',8F10.5,' )',F15.8)
      CALL NORM(ZZ,8,AB)

      IF(AB(1).EQ.0.D0.OR.AB(3).EQ.0.D0.OR.AB(5).EQ.0.D0)THEN 
      WRITE(53,'(A)')' Normalisation of 8-eigenvectors crazy.'
      NCRAZY = 1
      RETURN
      ENDIF 

C=====Save the starting (damping-free) eigenvectors
      ZZSAVE = ZZ
C=====Get the eigenvectors at the end of the ring but without wind-back.
      ZZWOWB = MATMUL(TREV8WOWB,ZZ)
C-----------------------------------------------------------------------
C=====Now that we have the s-o eigenvectors we can project out the 3-modes
C     and, in particular, get n-n_0 for each particle. Then we can start the 
C     tracking as close as possible to s-o equilibrium. 

C     Use SORBVEC to get the vectors corresponding to the three orbital modes.
C     Do this using the ``real'' forms for the eigenvectors, i.e., don't
C     convert to complex form. Then must be careful with identifying the  
C     correct pairs of real vectors and with signs. 
      ZZCONJ           =  TRANSPOSE(ZZSAVE(1:6,1:6))
      TEMP(1,1:NPART2) = -SORBVEC(2,1:NPART2)
      TEMP(2,1:NPART2) =  SORBVEC(1,1:NPART2)
      TEMP(3,1:NPART2) = -SORBVEC(4,1:NPART2)
      TEMP(4,1:NPART2) =  SORBVEC(3,1:NPART2)
      TEMP(5,1:NPART2) = -SORBVEC(6,1:NPART2)
      TEMP(6,1:NPART2) =  SORBVEC(5,1:NPART2)
      ORBAMP(:,1:NPART2) = MATMUL(ZZCONJ,TEMP(:,1:NPART2))

C=====Now set the (spin) alphas and betas for the 3 modes. 

      WRITE(53,103)
      DO 66 INP = 1,NPART2

      ORBVEC1(:,INP) = 2.D0*
     + (ORBAMP(2,INP)*ZZSAVE(1:6,1) - ORBAMP(1,INP)*ZZSAVE(1:6,2))
      ORBVEC2(:,INP) = 2.D0*
     + (ORBAMP(4,INP)*ZZSAVE(1:6,3) - ORBAMP(3,INP)*ZZSAVE(1:6,4))
      ORBVEC3(:,INP) = 2.D0*
     + (ORBAMP(6,INP)*ZZSAVE(1:6,5) - ORBAMP(5,INP)*ZZSAVE(1:6,6))

      SPINVEC1(:,INP) = 2.D0*
     + (ORBAMP(2,INP)*ZZSAVE(7:8,1) - ORBAMP(1,INP)*ZZSAVE(7:8,2))
      SPINVEC2(:,INP) = 2.D0*
     + (ORBAMP(4,INP)*ZZSAVE(7:8,3) - ORBAMP(3,INP)*ZZSAVE(7:8,4))
      SPINVEC3(:,INP) = 2.D0*
     + (ORBAMP(6,INP)*ZZSAVE(7:8,5) - ORBAMP(5,INP)*ZZSAVE(7:8,6))

      SORBVEC(7:8,INP) = SPINVEC1(:,INP)+SPINVEC2(:,INP)+SPINVEC3(:,INP)


C====Plot scatter of the starting spin angles alpha and beta at required energy point.
C    Differential scaling to help display.
      NSPINPLOT = 0
      IF(IE0.EQ.NDIFFDUMP)THEN
      WRITE(60,105)NSPINPLOT,INP,SORBVEC(7,INP)*1.D+1,
     +                            SORBVEC(8,INP)*1.D-1
  105 FORMAT(' ',2I10, 6E16.6) 
      ENDIF


      IF(INP.LE.5)THEN
      WRITE(53,'(A)')
     + ' A sample of starting spin-orbit vectors (mm,mrad..)'
      WRITE(53,'(A,A,I5,8F12.6)')' ', 
     +     'Spin-orbit vector:    all',INP,SORBVEC(1:8,INP)
      WRITE(53,'(A,A,I5,8F12.6)')' ',
     +     'Spin-orbit vector: mode-1',
     +                          INP,ORBVEC1(1:6,INP),SPINVEC1(:,INP)
      WRITE(53,'(A,A,I5,8F12.6)')' ',
     +     'Spin-orbit vector: mode-2',
     +                          INP,ORBVEC2(1:6,INP),SPINVEC2(:,INP)
      WRITE(53,'(A,A,I5,8F12.6)')' ',
     +     'Spin-orbit vector: mode-3',
     +                          INP,ORBVEC3(1:6,INP),SPINVEC3(:,INP)
      ENDIF
   66 CONTINUE  

C----------------------------------------------------------------------
C=====Get the covariance matrices for this sample with clever use of MATMUL
      SORBVET(1:NPART2,:) = TRANSPOSE(SORBVEC(:,1:NPART2))
      COVSIM = MATMUL(SORBVEC(:,1:NPART2),SORBVET(1:NPART2,:))/NPART2
      WRITE(53,899)

  899 FORMAT(//,' The generated spin-orbit covariance matrix',
     +     ' <Xi*Xj>(mm*mm,mm*mrad...) at the 1-st beam line element:')
      WRITE(53,955)COVSIM
  955 FORMAT(T6,8F13.5)


      SPINVET1(1:NPART2,:)= TRANSPOSE(SPINVEC1(:,1:NPART2))
      SPINVET2(1:NPART2,:)= TRANSPOSE(SPINVEC2(:,1:NPART2))
      SPINVET3(1:NPART2,:)= TRANSPOSE(SPINVEC3(:,1:NPART2))
      COVSPIN1(:,:,1)
     +     = MATMUL(SPINVEC1(:,1:NPART2),SPINVET1(1:NPART2,:))/NPART2
      COVSPIN2(:,:,1)
     +     = MATMUL(SPINVEC2(:,1:NPART2),SPINVET2(1:NPART2,:))/NPART2
      COVSPIN3(:,:,1)
     +     = MATMUL(SPINVEC3(:,1:NPART2),SPINVET3(1:NPART2,:))/NPART2

      WRITE(53,877)
  877 FORMAT(/,' The generated spin covariance matrix',
     +          ' mode-1 <Si*Sj>(mrad*mrad...)',
     +          ' at the 1-st beam line element:')
      WRITE(53,956)COVSPIN1(:,:,1)
      WRITE(53,878)
  878 FORMAT(/,' The generated spin covariance matrix',
     +          ' mode-2 <Si*Sj>(mrad*mrad...)',
     +          ' at the 1-st beam line element:')
      WRITE(53,956)COVSPIN2(:,:,1)
      WRITE(53,879)
  879 FORMAT(/,' The generated spin covariance matrix',
     +          ' mode-3 <Si*Sj>(mrad*mrad...)',
     +          ' at the 1-st beam line element:')
      WRITE(53,956)COVSPIN3(:,:,1)

C====Put the equilibrium alpha and beta variances on the ``slope plot''.
      KOL1 = 0
      IF(IE0.EQ.NDIFFDUMP)THEN
      WRITE(68,'(I10,40E13.4)')KOL1,
     +             COVSIM(7,7)+COVSIM(8,8),COVSIM(7,7),COVSIM(8,8)
      ENDIF

C      IF(1.EQ.1)STOP


C====Now set all particles on the C.O. with spins along n_0
C    so that a ``natural'' spin distribution can develop even if there
C    are strong orbital nonlinearities.
C    This also prevents potentially large oscillations of the two
C    SPINVEC covariances (although their sum doesn't oscillate even 
C    when the oscillations of the two parts are large).

      SORBVEC         = SORBVEC  *0.D0
C      SORBVEC(7:8,:)  = SORBVEC(7:8,:)  *3.D0

      SPINVEC1 = SPINVEC1 *0.D0
      SPINVEC2 = SPINVEC2 *0.D0
      SPINVEC3 = SPINVEC3 *0.D0




C----------------------------------------------------------------------
C======Get the 1-turn spin-(damped) orbit matrix.
      IDAMPFLG = 1
      S=0.
      CALL UNIT(8,TREV8D)
      TM3C = ZWSAVE
      DO 888  II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
      IF(IID.EQ.1.AND.XY.LT.0.00011D0)GO TO 888 !Skip zero length drifts 
      IF(IID.EQ.1)THEN
      TREV8D(1,:) = TREV8D(1,:) + XY * TREV8D(2,:) 
      TREV8D(3,:) = TREV8D(3,:) + XY * TREV8D(4,:) 
      GO TO 888
      ELSEIF(IID.EQ.10)THEN
      ZW = TM3C
      CALL SOL8AN(II,XY,YYY,NU,ZW,TM3C,CC,NSOL(ITY))
      ELSEIF(IID.EQ.17)THEN
      ZW = TM3C
      CC(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      CALL LBEAMBEAM(ITY,CC,ZW,NU)
      ELSE
      CC(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 
      CALL MX88DAMP(CC,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NU,TM3B,ZW,
     +                                         CRAD,NAME(ITY),IDAMPFLG)
      ENDIF
      TREV8D = MATMUL(CC,TREV8D)
  888 CONTINUE

C======WIND BACK THE SPIN BASIS==========================================
      CALL UNIT(8,TW8A)
      TW8A(7,7)= WR3(2)
      TW8A(7,8)= WI3(2)
      TW8A(8,7)=-WI3(2)
      TW8A(8,8)= WR3(2)
      TREV8D   = MATMUL(TW8A,TREV8D)


      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,/,A,I1,A)')' 1-turn spin-(damped)orbit matrix',
     +                               ' NDAMP2 = ',NDAMP2, ' !!!!!!!'
      DO 336 I=1,8
  336 WRITE(53,914)(TREV8D(I,J),J=1,8)



C----------------------------------------------------------------------------------
C======Make 1 pass around the ring with damping to generate damped
C      spin-orbit matrices for sections between the dipole centres using
C      the fact that the C.O. is fixed.  
C      So don't have to recalc all fixed stuff on each turn and can track 
C      using the sections.
C      Get the radiation scale factors which account for the various dipole  
C      strengths.
C      Since the beam-beam kicks can be nonlinear, it is simplest to treat
C      them separately and introduce additional sections. 
C      But for debugging by comparing the 1-turn maps, just use the LINEAR 
C      beam-beam maps at this stage. 
C      
      NSECT = 0                                  !Count sections  
      JSECT = 0                                  !Count non-dipole sections 
      IDAMPFLG = 1                               !Switch on/off radn.
      CALL UNIT(8,AAA)
      IBBSECT = 0
      TM3C = ZWSAVE
      RADTOT = 0.D0

      DO 7  II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
      CALL UNIT(8,CC)


C====If it's a beam-beam element, treat this as a section.
      IF(IID.EQ.17)THEN
      ZW = TM3C
      CC(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      CALL LBEAMBEAM(ITY,CC,ZW,NU)
      AAA = CC
      NSECT = NSECT + 1
      JSECT = JSECT + 1
      IBBSECT(NSECT) = 1                        ! Label it as a b-b section
      ZWBB(:,:,NSECT) = ZW                         ! Store the spin basis.
      XXBB(NSECT) = XY
      X2BB(NSECT) = XX2
      YYBB(NSECT) = YYY
      NSOLBB(NSECT) =NSOL(ITY) 
      IF(NSECT.GT.LIMSECT)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')' STOP: at b-b number of sections exceeds LIMSECT'
      STOP
      ENDIF
      SECTMAP(:,:,NSECT)   = AAA
      PSCALE(NSECT) = RADSCALE(ITY)
      IF(IE0.EQ.1) WRITE(53,'(A,3I10,2A,E16.5)')
     +                     ' ',II,IID,NSECT,'  ',NAME(ITY),PSCALE(NSECT)
      RADTOT = RADTOT + PSCALE(NSECT) * SYNCON
      CALL UNIT(8,AAA)      
      GO TO 7
      ENDIF

      IF(IID.EQ.10)THEN
      ZW = TM3C
      CALL SOL8AN(II,XY,YYY,NU,ZW,TM3C,CC,NSOL(ITY))
      ELSE
      CC(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 
      CALL MX88DAMP(CC,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NU,TM3B,ZW,
     +                                         CRAD,NAME(ITY),IDAMPFLG)
      ENDIF
      AAA = MATMUL(CC,AAA)
      IF((IID.EQ.7.AND.NAME(ITY)(1:2).EQ.'VD').
     +                OR.II.EQ.NELEM.
     +                OR.(II.NE.NELEM.AND.ID(ITYPE(II+1)).EQ.17))THEN  
      IF(II.EQ.NELEM)JSECT = JSECT + 1 
      NSECT = NSECT + 1
      IBBSECT(NSECT) = 0
      IF(NSECT.GT.LIMSECT)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')' STOP: the number of sections exceeds LIMSECT'
      STOP
      ENDIF

      SECTMAP(:,:,NSECT)   = AAA
      PSCALE(NSECT) = RADSCALE(ITY)
      IF(IE0.EQ.1) WRITE(53,'(A,3I10,2A,E16.5)')
     +                    ' ',II,IID,NSECT,'  ',NAME(ITY),PSCALE(NSECT)

      RADTOT = RADTOT + PSCALE(NSECT) * SYNCON
      CALL UNIT(8,AAA)      
      ENDIF

    7 CONTINUE

      WRITE(53,'(A,E16.5,A,E16.5)')' Total radiation excitation: ',
     +                                                       RADTOT


      RADAVE = RADTOT/(NSECT-JSECT)



C======Check the 1-turn matrix using the sections.
      CALL UNIT(8,AAA)
      DO 3 II = 1,NSECT
      AAA = MATMUL(SECTMAP(:,:,II),AAA)
      IF(II.LE.-20)THEN
      WRITE(53,103)
      WRITE(53,'(A,I5,A)')' ',II,' section matrix.'
      DO 337 I=1,8
  337 WRITE(53,914)(SECTMAP(I,J,II),J=1,8)
      ENDIF
    3 CONTINUE
C======WIND BACK THE SPIN BASIS==========================================
      AAA      = MATMUL(TW8A,AAA)

      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,A,I6)')' ','NSECT', NSECT
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,A,/,A,I1,A)')' 1-turn spin-(damped)orbit matrix',
     +           ' using sections.', ' NDAMP2 = ',NDAMP2, ' !!!!!!!'
      DO 35 I=1,8
   35 WRITE(53,914)(AAA(I,J),J=1,8)
  914 FORMAT(8F12.5)


C      IF(1.EQ.1)STOP


C-------------------------------------------------------------------------------------
C======Make another pass around the ring with NO damping to generate s-o eigenvectors
C      at the start of each section.
C      So don't have to recalc all fixed eigenvectors on each turn.
C      This is basically a repeat of the above stuff on sections above but 
C      without the damping.
C      Just use the linear b-b.
      MSECT = 0
      IDAMPFLG = 0                               !Switch on/off radn.
      CALL UNIT(8,AAA)
      TM3C = ZWSAVE
      ZZSECT(:,:,1) = ZZSAVE


      DO 77  II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
C====If it's a beam-beam element, treat this as a section.
      IF(IID.EQ.17)THEN
      ZW = TM3C
      CC(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      CALL LBEAMBEAM(ITY,CC,ZW,NU)
      AAA = CC
      MSECT = MSECT + 1
      IF(NSECT.GT.LIMSECT)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')' STOP: at b-b number of sections exceeds LIMSECT'
      STOP
      ENDIF
      ZZSECT(:,:,MSECT+1)  = MATMUL(AAA,ZZSECT(:,:,MSECT))
      CALL UNIT(8,AAA)      
      GO TO 77
      ENDIF

      IF(IID.EQ.10)THEN
      ZW = TM3C
      CALL SOL8AN(II,XY,YYY,NU,ZW,TM3C,CC,NSOL(ITY))
      ELSE
      CC(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 
      CALL MX88DAMP(CC,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NU,TM3B,ZW,
     +                                         CRAD,NAME(ITY),IDAMPFLG)
      ENDIF
      AAA = MATMUL(CC,AAA)
      IF((IID.EQ.7.AND.NAME(ITY)(1:2).EQ.'VD').
     +                OR.II.EQ.NELEM.
     +                OR.(II.NE.NELEM.AND.ID(ITYPE(II+1)).EQ.17))THEN
      MSECT = MSECT + 1
      IF(MSECT.GT.LIMSECT)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')' STOP: the number of sections exceeds LIMSECT'
      STOP
      ENDIF

      ZZSECT(:,:,MSECT+1)  = MATMUL(AAA,ZZSECT(:,:,MSECT))
      CALL UNIT(8,AAA)      
      ENDIF

   77 CONTINUE


C-------------------------------------------------------------------------------- 
      WRITE(53,103)
      WRITE(53,103)

      WRITE(53,9244)
 9244 FORMAT(////,' Eigenvalues and eigenvectors --',
     +               ' no damping, original, i.e., at the start:',/,T8,
     +                                          'REAL',
     +                                          T18,'IMAG',T122,'TUNES')
      DO 2811 I=1,8
      TUN=999999.
C=====GET TUNES: THIS DOES NOT DIFFERENTIATE BETWEEN +/- ANGLES
C=====WITHOUT USING THE SINE.
C     IF(DABS(WR8(I)) .LE. 1.D0)TUN=DACOS(WR8(I))/(2.*PI)
      TUN=DATAN2(WI8(I),WR8(I))/(2.*PI)
      TN(I)=TUN
 2811 WRITE(53,923)WR8(I),WI8(I),(ZZSAVE(J,I),J=1,8),TUN


      WRITE(53,9245)
 9245 FORMAT(////,' Eigenvalues and eigenvectors --', 
     +                 ' no damping, after 1 turn with sections:',/,T8,
     +                                          'REAL',
     +                                          T18,'IMAG',T122,'TUNES')

      DO 2812 I=1,8
      TUN=999999.
C=====GET TUNES: THIS DOES NOT DIFFERENTIATE BETWEEN +/- ANGLES
C=====WITHOUT USING THE SINE.
C     IF(DABS(WR8(I)) .LE. 1.D0)TUN=DACOS(WR8(I))/(2.*PI)
      TUN=DATAN2(WI8(I),WR8(I))/(2.*PI)
      TN(I)=TUN
 2812 WRITE(53,923)WR8(I),WI8(I),(ZZSECT(J,I,MSECT+1),J=1,8),TUN









C      WRITE(53,'(A,A,F20.8)')' ','CRAD',CRADIN

C======Now track to get the equilibrium beam and the non-equilibrium spin 
C      components.
C======After each dipole or C-F, radiate ``big'' photons according to a 
C      centred top-hat distribution. 
C      In contrast to a Gaussian, there are no long tails so that no cut 
C      against too large (particle expelling kicks) is needed.
C      For linear motion this should also be enough to get Gaussian phase space
C      by the central limit theorem. 
C      Use the seed, KSEED(14).
C      Use the Boege formula for the strength. The rms ``big'' photon energy is
C      very much smaller than the energy spread of the beam.
C      Scale the `big'' photons according to length/|rho|**3
C      This can be important for dipoles where n_0 is horizontal and there 
C      is no spin match.
C      No radiation in the correction coils. 
C      Scale the damping up with a factor NDAMP2 to decrease the damping time  
C      and thus get answers quicker. 
C
C      For linear tracking, no need to scale the s-orbit vector to metres etc. 
C      But do it anyway to prepare for the future when nonlinear stuff is
C      used. Then use metre,radian,seconds for all the remaining calcs but 
C      print in mm, mrad.,msecs.
C      NSECT includes the starting point (IP). But that give no radiation
C      So to get the size of big photons use NSECT-1.
C
C      We work in the periodic dreibein n_0, m, l. Then we can save computing 
C      time by setting up the 8x8 matrices for each section just once. 
C      If we were to use the nonperiodic dreibein n_0, m_0, l_0 , we would
C      need new 8x8 matrices for each turn.
C      Then we would see that the rms diffusions of alpha and beta along the  
C      m_0, l_0 axes would change from turn to turn as the G matrix changed 
C      and would differ from that given by SLICK, which is for just the first turn. 
C      In fact, one would find that the rms diffusions of alpha and beta 
C      would become equal after enough turns, as the rows of the G matrix   
C      interchanged many times. 
C      Since we work with the periodic dreibein n_0, m, l, we have a fixed
C      G matrix, but we wind back the alphas and betas at the end of each turn. 
C      This also has the effect of equalising the rms diffusions of alpha and 
C      beta, as evidenced by the numerical results. 
C      Thus to compare with the (first turn) predictions of SLICK, one must
C      just look at the diffusion away from n at the end of the first turn and
C      do that using a huge number of particles. Then, to avoid wasting time,
C      use the setting ``NTURN2 = 1'' to set all loops to just one pass. 
C
C
      SORBVEC  = SORBVEC*1.D-3
      SPINVEC1 = SPINVEC1*1.D-3
      SPINVEC2 = SPINVEC2*1.D-3
      SPINVEC3 = SPINVEC3*1.D-3


      MODIBMBM = IABS(IBMBM)


      IF(IE0.EQ.1)THEN
      COVMAT  = COVMATIN *1.D-6
      SCALEMAT= SCALEMAT *1.D-6
      TDAMP   = TDAMP*1.D-3
      ENDIF      

      CTIME =  CIR/3.0D8
      DAMPTURNS = TDAMP(5)/CTIME/NDAMP2


      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,F8.1,A,I3)')' Turns per sync. damping time', 
     +                       DAMPTURNS,'  NDAMP2=', NDAMP2
      NDAMPTURNS = DAMPTURNS   


      IF(NTURN2*NDAMPTURNS.GT.MAXT)THEN
      WRITE(53,103)
      WRITE(53,'(3A)')' STOP: attempt to work with more than MAXT turns'
      STOP
      ENDIF

      IDAMPFLG = 1                              !Switch on/off radn.
      TH = DSQRT(3.D0)                          !Gives a top hat of unit variance. 
      VARTH = TH*TH/3.D0                        !The variance for a top hat distn.


C    Get the scale factor needed to tweak the big photons. 
C    PHOTSCALE should be close to the sqrt of NDAMP2.
      PHOTSCALE = CTIME*SCALEMAT(6,6)/(TDAMP(5)*RADTOT)
      PHOTSCALE = 2.D0*DSQRT(PHOTSCALE*NDAMP2)*IDAMPFLG
      TH = TH * PHOTSCALE

      WRITE(53,103)
      WRITE(53,103)

      WRITE(53,'(A,F15.7,/,A,F15.7,/,A,F15.7)')
     +  ' Variance of generator             ', VARTH, 
     +  ' Required rescaling factor         ', PHOTSCALE,
     +  ' (r.m.s. kick)/(rel energy spread) ', 
     +       DSQRT(RADAVE)/DSQRT(SCALEMAT(6,6))
      WRITE(53,103)

      IF(IE0.EQ.1)WRITE(53,'(A,A,A,I6,3E15.7)')' ','NSECT, PHOTSCALE,', 
     +                         ' (r.m.s. kick)/(rel energy spread) ', 
     +     NSECT,PHOTSCALE,TH/DSQRT(3.D0)/DSQRT(COVMAT(6,6))

      WRITE(53,103)


      CALL UNIT(2,TW2A)
      TW2A(1,1)= WR3(2)
      TW2A(1,2)= WI3(2)
      TW2A(2,1)=-WI3(2)
      TW2A(2,2)= WR3(2)


C      IF(IE0.EQ.1)
C      CALL G05CBF(KSEED14)             !Set the seed for the big photons.

C===Special action if NTURN2 = 1
      IF(NTURN2.EQ.1)NDAMPTURNS = 1

      IFFT = 0
      XFFT = 0.D0
      YFFT = 0.D0
      IDNT = 0
      NSPINPLOT = 0 

      DO 9 INT = 1,NTURN2               !Loop over damping periods.
      DO 4 IDT = 1,NDAMPTURNS           !Loop over turns in a damping period. 

      IDNT = IDNT + 1
      IF(IDNT.GT.MAXT)THEN
      WRITE(53,'(A)')' Number of turns exceeds MAXT. SO STOP'
      STOP
      ENDIF

C===Do a FFT on the transverse motion to discover the effect
C   of the nonlinearity of the beam-beam force.
C   Use just 1 particle and sample after the first 5 damping times.   
C   Sample at the IP. 
CCCC      IF(IE0.EQ.1.AND.INT.GT.5)THEN 
CCCC      IF(IFFT.LE.NFFT-1)THEN 
CCCC      XFFT(IFFT) = SORBVEC(1,1)
CCCCC      XFFT(IFFT) = DSIN(2.D0*PI*0.40D0*IFFT) + DSIN(2.D0*PI*0.25D0*IFFT)
CCCC      YFFT(IFFT) = SORBVEC(3,1)    
CCCCC      WRITE(53,'(A,3I10,2F16.6)')' Reached here 2: ',
CCCCC     +                       IE0,INT,IFFT,XFFT(IFFT),YFFT(IFFT)

CCCC      IFFT = IFFT + 1
CCCC      KIFFT = IFFT                   !Final number of data points.
CCCC      ENDIF
CCCC      ENDIF


      DO 5 II  = 1,NSECT                !Loop over sections in a turn.

C=====Use SORBVEC to get the vectors corresponding to the three orbital modes.
C     Do this using the ``real'' forms for the eigenvectors, i.e., don't
C     convert to complex form. Then must be careful with identifying the  
C     correct pairs of real vectors and with signs. 
C     Check by rebuilding the original vector.
      ZZCONJ(1:6,1:6)  =  ZZSECT(1:6,1:6,II)
      ZZCONJ           =  TRANSPOSE(ZZCONJ)
      IF(MODES.EQ.1)THEN
      TEMP(1,1:NPART2) = -SORBVEC(2,1:NPART2)
      TEMP(2,1:NPART2) =  SORBVEC(1,1:NPART2)
      TEMP(3,1:NPART2) = -SORBVEC(4,1:NPART2)
      TEMP(4,1:NPART2) =  SORBVEC(3,1:NPART2)
      TEMP(5,1:NPART2) = -SORBVEC(6,1:NPART2)
      TEMP(6,1:NPART2) =  SORBVEC(5,1:NPART2)
      ORBAMP(:,1:NPART2) = MATMUL(ZZCONJ,TEMP(:,1:NPART2))
      ENDIF

C=====Now reconstruct the 3 real orbit vectors. 
      IF(MODES.EQ.1)THEN
      DO 6 INP = 1,NPART2
      ORBVEC1(:,INP) = 2.D0*
     + (ORBAMP(2,INP)*ZZSECT(1:6,1,II) - ORBAMP(1,INP)*ZZSECT(1:6,2,II))
      ORBVEC2(:,INP) = 2.D0*
     + (ORBAMP(4,INP)*ZZSECT(1:6,3,II) - ORBAMP(3,INP)*ZZSECT(1:6,4,II))
      ORBVEC3(:,INP) = 2.D0*
     + (ORBAMP(6,INP)*ZZSECT(1:6,5,II) - ORBAMP(5,INP)*ZZSECT(1:6,6,II))
   6  CONTINUE  
      ENDIF 
C=====Sum the 3 vectors to check that one gets the original.
C      IF(II.EQ.77.AND.IDT.EQ.1.AND.INT.EQ.10)THEN
C      WRITE(53,'(6F11.7)')ZZCONJ
C      WRITE(53,'(6F11.7)')TEMP(:,5)
C      WRITE(53,'(6F11.7)')USYMP
C      WRITE(53,'(A,A,8F11.7)')' ',
C     +                'Current spin-orbit vector       ',SORBVEC(1:6,5)
C      WRITE(53,'(A,A,6F11.7)')' ',
C     +                'Orbit amplitudes                ',ORBAMP(:,5)
C      WRITE(53,'(A,A,6F11.7)')' ',
C     +                'orbit vector_1                  ',ORBVEC1(:,5)
C      WRITE(53,'(A,A,6F11.7)')' ',
C     +                'orbit vector_2                  ',ORBVEC2(:,5)
C      WRITE(53,'(A,A,6F11.7)')' ',
C     +                'orbit vector_3                  ',ORBVEC3(:,5)
C      TESTVEC(:) = ORBVEC1(:,5) + ORBVEC2(:,5) + ORBVEC3(:,5)
C      WRITE(53,'(A,A,6F11.7)')' ',
C     +                'Check vector                    ',TESTVEC
C
C      IF(1.EQ.1)STOP
C
C      ENDIF 

C===Do a FFT on the transverse motion to discover the effect
C   of the nonlinearity of the beam-beam force.
C   Use just 1 particle and sample after the first 5 damping times.   
C   Sample at the start of just 1 chosen section in the arc. 
      IF(IE0.EQ.1.AND.INT.GT.5.AND.II.EQ.40)THEN 
      IF(IFFT.LE.NFFT-1)THEN 
      XFFT(IFFT) = SORBVEC(1,1)
      YFFT(IFFT) = SORBVEC(3,1)    

      IFFT = IFFT + 1
      KIFFT = IFFT                   !Final number of data points.
      ENDIF
      ENDIF





C=====Now transport the various vectors.
C     For linear beam-beam effect, just use the linear section maps. 
C     Otherwise use the round or elliptical beam-beam maps.

      IF(IBBSECT(II).EQ.1)THEN    !If the section is a type 17 element. 
      IF(MODIBMBM.GE.2)THEN       !If we want full nonlinear b-b.
                                        !If IBMBM=0, this is skipped. 

      CALL NLBEAMBEAM(II,XXBB(II),X2BB(II),YYBB(II),NU,ZWBB,SORBVEC,
     +                                             NSOLBB(II),MODIBMBM)
C    Or test that we get exactly the same results as when using a type 3 
C    b-b quad, by inserting the linear b-b SECTMAP that is always set up. 
C    even if it's a type 17 element. If IBMBM=0, NLBEAMBEAM is skipped.
C    That's fine since the SECTMAP = identity.
      ELSE

      SORBVEC(:,1:NPART2) = MATMUL(SECTMAP(:,:,II),SORBVEC(:,1:NPART2))

      IF(MODES.EQ.1)THEN
      GMAT = SECTMAP(7:8,1:6,II)
      SPINVEC1(:,1:NPART2) = MATMUL(GMAT,ORBVEC1(:,1:NPART2))
     +                                      + SPINVEC1(:,1:NPART2)
      SPINVEC2(:,1:NPART2) = MATMUL(GMAT,ORBVEC2(:,1:NPART2))
     +                                      + SPINVEC2(:,1:NPART2)
      SPINVEC3(:,1:NPART2) = MATMUL(GMAT,ORBVEC3(:,1:NPART2))
     +                                      + SPINVEC3(:,1:NPART2)      
      ENDIF
      ENDIF

      ELSE

      DO 55553 IPP = 1, NPART2
      SORBVEC(:,IPP) = MATMUL(SECTMAP(:,:,II),SORBVEC(:,IPP))
55553 CONTINUE

9144  FORMAT(8F20.12)

      GMAT = SECTMAP(7:8,1:6,II)

      IF(MODES.EQ.1)THEN
      GMAT = SECTMAP(7:8,1:6,II)
      SPINVEC1(:,1:NPART2) = MATMUL(GMAT,ORBVEC1(:,1:NPART2))
     +                                      + SPINVEC1(:,1:NPART2)
      SPINVEC2(:,1:NPART2) = MATMUL(GMAT,ORBVEC2(:,1:NPART2))
     +                                      + SPINVEC2(:,1:NPART2)
      SPINVEC3(:,1:NPART2) = MATMUL(GMAT,ORBVEC3(:,1:NPART2))
     +                                      + SPINVEC3(:,1:NPART2)      
      ENDIF

      ENDIF



C=====These individual components are now computed earlier. 
C     ORBVEC1(:,1:NPART2)
C    +           = MATMUL(SECTMAP(1:6,1:6,II),ORBVEC1(:,1:NPART2))
C     ORBVEC2(:,1:NPART2)
C    +           = MATMUL(SECTMAP(1:6,1:6,II),ORBVEC2(:,1:NPART2))
C     ORBVEC3(:,1:NPART2)
C    +           = MATMUL(SECTMAP(1:6,1:6,II),ORBVEC3(:,1:NPART2))


C      IF(II.EQ.180)
C      WRITE(53,'(A,3I10,4E24.12)')'SECT,IDNT,IBBSECT ',
C     +                             II,IDNT,IBBSECT(II),
C     +                             SORBVEC(1,16),SORBVEC(2,16),
C     +                             SORBVEC(3,16),SORBVEC(4,16)

C      IF(IDNT.EQ.2)STOP

C=====Now radiate. But only at dipoles, not at the IP or b-b. Weight according to the
C     dipole strength. 
C     Apply the energy kick to the vector at the END of this section using the strength  
C     of the radiator marking the end of this section. Then we are ready for the next
C     section.


      IF(PSCALE(II).NE.0.D0)THEN
C      CALL G05FAF(-TH,TH,NPART2,BIGPHOT)
      DO 55 IPHS = 1,NPART2

      SORBVEC(6,IPHS) =  SORBVEC(6,IPHS) + BIGPHOT(IPHS)
     +                              * DSQRT(PSCALE(II)*SYNCON)*RADFAC
C      WRITE(53,'(A,I10,E20.8);')'Bigphot ',IPHS,BIGPHOT(IPHS)     
C      IF(IPHS.EQ.16)WRITE(53,'(A,E24.12)')'Photon ', 
C     +                 BIGPHOT(IPHS) * DSQRT(PSCALE(II)*SYNCON)*RADFAC  
   55 CONTINUE
      ENDIF

    5 CONTINUE



C=====Get the n vector for each new particle position. Need the eigenvectors
C     at the end of a turn without wind-back: ZZWOWB   
      ZZWOWBCONJ       =  TRANSPOSE(ZZWOWB(1:6,1:6))
      TEMP(1,1:NPART2) = -SORBVEC(2,1:NPART2)
      TEMP(2,1:NPART2) =  SORBVEC(1,1:NPART2)
      TEMP(3,1:NPART2) = -SORBVEC(4,1:NPART2)
      TEMP(4,1:NPART2) =  SORBVEC(3,1:NPART2)
      TEMP(5,1:NPART2) = -SORBVEC(6,1:NPART2)
      TEMP(6,1:NPART2) =  SORBVEC(5,1:NPART2)
      ORBAMP(:,1:NPART2) = MATMUL(ZZWOWBCONJ,TEMP(:,1:NPART2))

      IF(IDNT.EQ.1)THEN
      DO 67 INPR = 1,NPART2
      ENVEC1(:,INPR) = 2.D0*
     + (ORBAMP(2,INPR)*ZZWOWB(7:8,1) - ORBAMP(1,INPR)*ZZWOWB(7:8,2))
      ENVEC2(:,INPR) = 2.D0*
     + (ORBAMP(4,INPR)*ZZWOWB(7:8,3) - ORBAMP(3,INPR)*ZZWOWB(7:8,4))
      ENVEC3(:,INPR) = 2.D0*
     + (ORBAMP(6,INPR)*ZZWOWB(7:8,5) - ORBAMP(5,INPR)*ZZWOWB(7:8,6))
      ENVEC(:,INPR) = ENVEC1(:,INPR) + ENVEC2(:,INPR) + ENVEC3(:,INPR)
C===Subtract the n-vector components from the spin components
C   in order to measure diffusion away from n,  
      REDSPIN(:,INPR) = SORBVEC(7:8,INPR) - ENVEC(:,INPR)
   67 CONTINUE  
      ENDIF


C=====Wind back spin basis at the end of a turn but to reflect 
C     the calc in SPIN, only after logging the spin coordinates
C     for the ``spin ellipses''.
      SORBVEC(7:8,1:NPART2)  = MATMUL(TW2A,SORBVEC(7:8,1:NPART2))
      IF(MODES.EQ.1)THEN      
      SPINVEC1(:,1:NPART2)   = MATMUL(TW2A,SPINVEC1(:,1:NPART2))
      SPINVEC2(:,1:NPART2)   = MATMUL(TW2A,SPINVEC2(:,1:NPART2))
      SPINVEC3(:,1:NPART2)   = MATMUL(TW2A,SPINVEC3(:,1:NPART2))
      ENDIF



C===Get the covariance matrices.
      SORBVET(1:NPART2,:)    = TRANSPOSE(SORBVEC(:,1:NPART2))
      IF(MODES.EQ.1)THEN
      SPINVET1(1:NPART2,:)   = TRANSPOSE(SPINVEC1(:,1:NPART2))
      SPINVET2(1:NPART2,:)   = TRANSPOSE(SPINVEC2(:,1:NPART2))
      SPINVET3(1:NPART2,:)   = TRANSPOSE(SPINVEC3(:,1:NPART2))
      ENDIF
      IF(IDNT.EQ.1)
     +        REDSPINT(1:NPART2,:)   = TRANSPOSE(REDSPIN(:,1:NPART2))
      COVTRACK(:,:,IDNT) 
     +       = MATMUL(SORBVEC(:,1:NPART2),SORBVET(1:NPART2,:))/NPART2
      IF(MODES.EQ.1)THEN
      COVSPIN1(:,:,IDNT)
     +     = MATMUL(SPINVEC1(:,1:NPART2),SPINVET1(1:NPART2,:))/NPART2
      COVSPIN2(:,:,IDNT)
     +     = MATMUL(SPINVEC2(:,1:NPART2),SPINVET2(1:NPART2,:))/NPART2
      COVSPIN3(:,:,IDNT)
     +     = MATMUL(SPINVEC3(:,1:NPART2),SPINVET3(1:NPART2,:))/NPART2
      ENDIF
      IF(IDNT.EQ.1)COVREDSPIN(:,:,IDNT)
     +      = MATMUL(REDSPIN(:,1:NPART2),REDSPINT(1:NPART2,:))/NPART2

      SDIFFUSION(IDNT)     = COVTRACK(7,7,IDNT)+COVTRACK(8,8,IDNT)
      SDIFFUSIONA(IDNT)    = COVTRACK(7,7,IDNT)
      SDIFFUSIONB(IDNT)    = COVTRACK(8,8,IDNT)
      SDIFFUSIONAB(IDNT)   = COVTRACK(7,8,IDNT)

      IF(MODES.EQ.1)THEN
      SDIFFUSION1(IDNT)    = COVSPIN1(1,1,IDNT)+COVSPIN1(2,2,IDNT)
      SDIFFUSIONA1(IDNT)   = COVSPIN1(1,1,IDNT)
      SDIFFUSIONB1(IDNT)   = COVSPIN1(2,2,IDNT)
      SDIFFUSIONAB1(IDNT)  = COVSPIN1(1,2,IDNT)

      SDIFFUSION2(IDNT)    = COVSPIN2(1,1,IDNT)+COVSPIN2(2,2,IDNT)
      SDIFFUSIONA2(IDNT)   = COVSPIN2(1,1,IDNT)
      SDIFFUSIONB2(IDNT)   = COVSPIN2(2,2,IDNT)
      SDIFFUSIONAB2(IDNT)  = COVSPIN2(1,2,IDNT)

      SDIFFUSION3(IDNT)    = COVSPIN3(1,1,IDNT)+COVSPIN3(2,2,IDNT)
      SDIFFUSIONA3(IDNT)   = COVSPIN3(1,1,IDNT)
      SDIFFUSIONB3(IDNT)   = COVSPIN3(2,2,IDNT)
      SDIFFUSIONAB3(IDNT)  = COVSPIN3(1,2,IDNT)
      ENDIF

      DO 90 IK = 1,6
   90 ODIFFUSION(IK,IDNT) = COVTRACK(IK,IK,IDNT)/SCALEMAT(IK,IK)


C=====Wind back spin basis at the end of a turn but to reflect 
C     the calc in SPIN, only after logging the spin coordinates.
C      SORBVEC(7:8,1:NPART2)  = MATMUL(TW2A,SORBVEC(7:8,1:NPART2))
C      IF(MODES.EQ.1)THEN      
C      SPINVEC1(:,1:NPART2)   = MATMUL(TW2A,SPINVEC1(:,1:NPART2))
C      SPINVEC2(:,1:NPART2)   = MATMUL(TW2A,SPINVEC2(:,1:NPART2))
C      SPINVEC3(:,1:NPART2)   = MATMUL(TW2A,SPINVEC3(:,1:NPART2))
C      ENDIF


C=====Plot scatter of the spin angles alpha and beta at chosen energy step.
C     Use the wound back angles.
C     Differential scaling to help display and sideways shifts to aid viewing.
      IF(MOD(IDNT,1).EQ.0.AND.IDNT.GE.3980)NSPINPLOT = NSPINPLOT + 1 
C      IF(IE0.EQ.NDIFFDUMP.AND.MOD(IDNT,2).EQ.0.AND.IDNT.GE.3980)THEN
      IF(IE0.EQ.NDIFFDUMP.AND.MOD(IDNT,1).EQ.0.AND.IDNT.GE.NTND-10)THEN
C      DO 104 IOO = 1, NPART2 
C  104 WRITE(60,105)NSPINPLOT,IOO,SORBVEC(7,IOO)*1.D4+1.D3*NSPINPLOT,
C     +                            SORBVEC(8,IOO)*1.D2
      ENDIF





C===Dump out info on 1st turn diffusion.
C===This will only make sense for a huge (say 10000) number of particles.
C   Then use SPECIAL RUNS with just the first turn



      IF(NTURN2.EQ.1)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,3E12.3,A,E12.3)')
     +                    ' alpha-beta covariances after first turn:  ', 
     +   COVREDSPIN(1,1,IDNT),COVREDSPIN(2,2,IDNT),COVREDSPIN(1,2,IDNT),
     + '  alpha-beta ratio ',  COVREDSPIN(1,1,IDNT)/COVREDSPIN(2,2,IDNT)  !Careful if all radn is off

      WRITE(53,'(A,3E12.3,A,E12.3)')
     +   ' From SPIN: alpha-beta ellipse parameters: ',WD0A,WD0B,WD0AB,
     +'  alpha-beta ratio ',  WD0ABRATIO

      WRITE(64,'(I10,12E13.5)')IE0,NU,WD0A,WD0B,WD0AB,WD0ABRATIO,
     +  COVREDSPIN(1,1,IDNT),COVREDSPIN(2,2,IDNT),COVREDSPIN(1,2,IDNT), 
     +  COVREDSPIN(1,1,IDNT)/COVREDSPIN(2,2,IDNT)                         !Careful if all radn is off

C===Get major and minor axes of M-C 1-sigma spin ellipse from COVREDSPIN. 
      IF(IE0.EQ.1)THEN 
      STILT=DATAN2(2.D0*COVREDSPIN(1,2,IDNT),(COVREDSPIN(1,1,IDNT)        !Careful if all radn is off
     +                                    - COVREDSPIN(2,2,IDNT)))        ! 2 theta
      CSTILT=DCOS(STILT)                                                  !cos (2 theta)
C      CSTILT = 0.5                                 !A fix if all radn is off 
      STILT=STILT*90.D0/PI
      AB3=COVREDSPIN(1,1,IDNT)*COVREDSPIN(2,2,IDNT)
     +                                      - COVREDSPIN(1,2,IDNT)**2
      IF(AB3.LT.0.D0.AND.AB3.GT.-1.D-8)AB3=0.D0                !Are these gymnastics needed?
      SAXISMAJ   = DABS(AB3/ 
     +            ( 0.5D0*((COVREDSPIN(2,2,IDNT)-COVREDSPIN(1,1,IDNT))
     +            /CSTILT + COVREDSPIN(2,2,IDNT)+COVREDSPIN(1,1,IDNT))))
      SAXISMAJ   = DSQRT(SAXISMAJ)
      SAXISMIN   = DABS(AB3/
     +            (-0.5D0*((COVREDSPIN(2,2,IDNT)-COVREDSPIN(1,1,IDNT))
     +            /CSTILT - COVREDSPIN(2,2,IDNT)-COVREDSPIN(1,1,IDNT))))
      SAXISMIN   = DSQRT(SAXISMIN)
C===Plot the ellipse at the starting point. Take 1000 points. Also 
C   principle axes.
      SDLTA  = STILT*PI/180.D0                     
      DO 71 IPHS = 1,1000
      SPH    = 0.001D0*2.D0*PI*(IPHS-1)
      SHR    = SAXISMAJ*DCOS(SPH)    
      SVR    = SAXISMIN*DSIN(SPH)    
      SHOR   = SHR*DCOS(SDLTA       ) - SVR*DSIN(SDLTA       )         !Turn off tilt
      SVER   = SHR*DSIN(SDLTA       ) + SVR*DCOS(SDLTA       )         !Turn off tilt
      SAXMAJ =    (-1.D0 + 0.002D0*(IPHS-1))*DTAN(SDLTA)
      SAXMIN =    (-1.D0 + 0.002D0*(IPHS-1))*DTAN(SDLTA+PI/2.D0)
      WRITE(66,105)IPHS,SHOR,SVER,-1.D0+0.002D0*(IPHS-1),SAXMAJ,SAXMIN
   71 CONTINUE  

C   Skip the following trick. Now do things properly with WD0AB in SPIN. 
C===Get major and minor axes of THEORETICAL 1-sigma spin ellipse from SPIN.
C   Transform the info from SPIN to simulate the rotation STILT of the m,l  
C   vectors, so that the ellipse has no tilt and can be compared with the M-C.
C   First transform the raw SPIN data. 
      IF(IE0.LT.-1)THEN
      AXISMAJ  =  (DCOS(SDLTA))**2 * WD0A + (DSIN(SDLTA))**2 * WD0B
     +                                  +  WD0AB * SIN(2.D0*SDLTA)
      AXISMIN  =  (DSIN(SDLTA))**2 * WD0A + (DCOS(SDLTA))**2 * WD0B
     +                                  -  WD0AB * SIN(2.D0*SDLTA)
      AXISMAJ  =  DSQRT(AXISMAJ)
      AXISMIN  =  DSQRT(AXISMIN)

C===Plot the SPIN ellipse after one turn. Take 1000 points. Also principle axes.
C   Now plot, using the tilt from the M-C!
      DO 70 IPH = 1,1000
      PH    = 0.001D0*2.D0*PI*(IPH-1)
      HR    = AXISMAJ*DCOS(PH)    
      VR    = AXISMIN*DSIN(PH)    
      HOR   = HR*DCOS(SDLTA) - VR*DSIN(SDLTA)
      VER   = HR*DSIN(SDLTA) + VR*DCOS(SDLTA) 
      AXMAJ =    (-1.D0 + 0.002D0*(IPH-1))*DTAN(SDLTA)
      AXMIN =    (-1.D0 + 0.002D0*(IPH-1))*DTAN(SDLTA+PI/2.D0)
C      WRITE(65,1105)IPH,HOR,VER,-1.D0+0.002D0*(IPH-1),AXMAJ,AXMIN
 1105 FORMAT(' ',I10, 6E16.6) 
   70 CONTINUE  
      ENDIF


      ENDIF
      ENDIF 


    4 CONTINUE 


C====At the end of each damping time, write out sample info

      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,86)INT,IDNT
   86 FORMAT(//,' The ',I3,'th  transported s-o covariance matrix',
     +          ' <Xi*Xj>(mm*mm,mm*mrad...)',
     +          ' at the 1-st beam line element after ',I6,' turns')
      WRITE(53,955)COVTRACK(:,:,IDNT)*1.D6

C      IF(1.EQ.1)STOP


      IF(MODES.EQ.1)THEN
      WRITE(53,87)INT,IDNT
   87 FORMAT(/,' The ',I3,'th  transported spin covariance matrix',
     +          ' mode-1 <Si*Sj>(mrad*mrad...)',
     +          ' at the 1-st beam line element after ',I6,' turns')
      WRITE(53,956)COVSPIN1(:,:,IDNT)*1.D6
  956 FORMAT(T6,2F16.3)
      WRITE(53,88)INT,IDNT
   88 FORMAT(/,' The ',I3,'th  transported spin covariance matrix',
     +          ' mode-2 <Si*Sj>(mrad*mrad...)',
     +          ' at the 1-st beam line element after ',I6,' turns')
      WRITE(53,956)COVSPIN2(:,:,IDNT)*1.D6
      WRITE(53,89)INT,IDNT
   89 FORMAT(/,' The ',I3,'th  transported spin covariance matrix',
     +          ' mode-3 <Si*Sj>(mrad*mrad...)',
     +          ' at the 1-st beam line element after ',I6,' turns')
      WRITE(53,956)COVSPIN3(:,:,IDNT)*1.D6
      ENDIF

      WRITE(53,103)
C      WRITE(53,'(A,A,F13.5,A,F13.5)')' ','(7,7) + (8,8): ',
C     +                             SDIFFUSION(IDNT),'    (7,7)/(8,8): ',
C     +                         COVTRACK(7,7,IDNT)/COVTRACK(8,8,IDNT)



    9 CONTINUE

C======Check the variance of the last big photons. 
      WRITE(53,103)
      WRITE(53,103)
      VARPHOT = 0
      IF(PHOTSCALE.GT.0.D0)THEN
      DO 16 IPH = 1,NPART2
   16 VARPHOT = VARPHOT + BIGPHOT(IPH)**2
      VARPHOT = VARPHOT/NPART2/PHOTSCALE**2
      ENDIF
      WRITE(53,'(A,A,F10.6)')' ','Variance of the top-hat samples ',
     +                           VARPHOT



C      IF(1.EQ.1)STOP


C-------------------------------------------------------------------------
C      Print some statistical stuff. Rescale to mm. etc.

      ONE = 1.D0
      SORBAVE = MATMUL(SORBVEC(:,1:NPART2),ONE(1:NPART2))/NPART2
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
C      DO 101 IO = 1,NPART2
C      NGAUSX = SORBVEC(1,IO)*1.D3/XYBIN
C      KGBINX(NGAUSX) = KGBINX(NGAUSX) + 1
C      NGAUSY = SORBVEC(3,IO)*1.D3/XYBIN
C      KGBINY(NGAUSY) = KGBINY(NGAUSY) + 1
C  101 CONTINUE 
C      DO 102 JO = -500,500
C  102 WRITE(59,'(1X,I10,1X,I10,1X,I10)')JO,KGBINX(JO),KGBINY(JO)

C----------------------------------------------------------------------------------------
C======Plot histograms of horizontal and vertical projections: 
      IF(IBEQUIL.EQ.0)THEN
      XYBIN  = 0.01D0
      DO 101 IO = 1,NPART2
      NGAUSX = SORBVEC(1,IO)*1.D3/XYBIN
      KGBINX(NGAUSX) = KGBINX(NGAUSX) + 1
      NGAUSY = SORBVEC(3,IO)*1.D3/XYBIN
      KGBINY(NGAUSY) = KGBINY(NGAUSY) + 1
  101 CONTINUE 
      DO 102 JO = -500,500
  102 WRITE(59,'(1X,I10,1X,I10,1X,I10)')JO,KGBINX(JO),KGBINY(JO)
      ENDIF
C======Plot scatter of the spin angles alpha and beta at first energy point.
C      IF(IE0.EQ.1)THEN
C      DO 104 IOO = 1, NPART2 
C  104 WRITE(60,105)IOO,SORBVEC(7,IOO),SORBVEC(8,IOO)
C  105 FORMAT(' ',I10, 6E16.6) 
C      ENDIF


C======Plot the rms diffusion angle vs. turns and do a LSQ fit.
C      If the requested NTURN2 is too small, skip this. 
      IF(NTURN2.GT.3)THEN
      DIFFMAT  = 0.D0  
      DVEC     = 0.D0
      DVECA    = 0.D0
      DVECB    = 0.D0
      DVECAB   = 0.D0
      DVEC1    = 0.D0
      DVECA1   = 0.D0
      DVECB1   = 0.D0
      DVECAB1  = 0.D0
      DVEC2    = 0.D0
      DVECA2   = 0.D0
      DVECB2   = 0.D0
      DVECAB2  = 0.D0
      DVEC3    = 0.D0
      DVECA3   = 0.D0
      DVECB3   = 0.D0
      DVECAB3  = 0.D0


C      IF(1.EQ.1)STOP


C====Now do fits. Skip data from the first few damping times.

      NSKIPT = NSKIP * NDAMPTURNS
      DO 106 JF = NSKIPT+1,IDNT
      TJF = JF
      DIFFMAT(1,1) = DIFFMAT(1,1) + TJF*TJF   ! This can overflow if real*4 integers are used.
      DIFFMAT(1,2) = DIFFMAT(1,2) + TJF
      DIFFMAT(2,2) = DIFFMAT(2,2) + 1.D0
      DVEC(1)   = DVEC(1)   + JF*SDIFFUSION(JF)
      DVEC(2)   = DVEC(2)   +    SDIFFUSION(JF)
      DVEC(3)   = DVEC(3)   +    SDIFFUSION(JF) * SDIFFUSION(JF)      
      DVECA(1)  = DVECA(1)  + JF*SDIFFUSIONA(JF)
      DVECA(2)  = DVECA(2)  +    SDIFFUSIONA(JF)
      DVECB(1)  = DVECB(1)  + JF*SDIFFUSIONB(JF)
      DVECB(2)  = DVECB(2)  +    SDIFFUSIONB(JF)
      DVECAB(1) = DVECAB(1) + JF*SDIFFUSIONAB(JF)
      DVECAB(2) = DVECAB(2) +    SDIFFUSIONAB(JF)

      DVEC1(1)   = DVEC1(1)   + JF*SDIFFUSION1(JF)
      DVEC1(2)   = DVEC1(2)   +    SDIFFUSION1(JF)
      DVEC1(3)   = DVEC1(3)   +    SDIFFUSION1(JF) * SDIFFUSION1(JF)
      DVECA1(1)  = DVECA1(1)  + JF*SDIFFUSIONA1(JF)
      DVECA1(2)  = DVECA1(2)  +    SDIFFUSIONA1(JF)
      DVECB1(1)  = DVECB1(1)  + JF*SDIFFUSIONB1(JF)
      DVECB1(2)  = DVECB1(2)  +    SDIFFUSIONB1(JF)
      DVECAB1(1) = DVECAB1(1) + JF*SDIFFUSIONAB1(JF)
      DVECAB1(2) = DVECAB1(2) +    SDIFFUSIONAB1(JF)

      DVEC2(1)   = DVEC2(1)   + JF*SDIFFUSION2(JF)
      DVEC2(2)   = DVEC2(2)   +    SDIFFUSION2(JF)
      DVEC2(3)   = DVEC2(3)   +    SDIFFUSION2(JF) * SDIFFUSION2(JF)
      DVECA2(1)  = DVECA2(1)  + JF*SDIFFUSIONA2(JF)
      DVECA2(2)  = DVECA2(2)  +    SDIFFUSIONA2(JF)
      DVECB2(1)  = DVECB2(1)  + JF*SDIFFUSIONB2(JF)
      DVECB2(2)  = DVECB2(2)  +    SDIFFUSIONB2(JF)
      DVECAB2(1) = DVECAB2(1) + JF*SDIFFUSIONAB2(JF)
      DVECAB2(2) = DVECAB2(2) +    SDIFFUSIONAB2(JF)

      DVEC3(1)   = DVEC3(1)   + JF*SDIFFUSION3(JF)
      DVEC3(2)   = DVEC3(2)   +    SDIFFUSION3(JF)
      DVEC3(3)   = DVEC3(3)   +    SDIFFUSION3(JF) * SDIFFUSION3(JF)
      DVECA3(1)  = DVECA3(1)  + JF*SDIFFUSIONA3(JF)
      DVECA3(2)  = DVECA3(2)  +    SDIFFUSIONA3(JF)
      DVECB3(1)  = DVECB3(1)  + JF*SDIFFUSIONB3(JF)
      DVECB3(2)  = DVECB3(2)  +    SDIFFUSIONB3(JF)
      DVECAB3(1) = DVECAB3(1) + JF*SDIFFUSIONAB3(JF)
      DVECAB3(2) = DVECAB3(2) +    SDIFFUSIONAB3(JF)
  106 CONTINUE 

C      IF(1.EQ.1)STOP


      DIFFMAT(2,1) = DIFFMAT(1,2)
      DET = DIFFMAT(1,1)*DIFFMAT(2,2) - DIFFMAT(1,2)**2
      SXX = DIFFMAT(1,1) - (DIFFMAT(1,2)**2/DIFFMAT(2,2))
      DFIT(1)  =  ( DIFFMAT(2,2)*DVEC(1)   - DIFFMAT(1,2)*DVEC(2))/DET
      DFIT(2)  =  (-DIFFMAT(2,1)*DVEC(1)   + DIFFMAT(1,1)*DVEC(2))/DET
      SYY = DVEC(3) - DVEC(2)**2/DIFFMAT(2,2)
      SXY = DVEC(1) - DVEC(2)*DIFFMAT(1,2)/DIFFMAT(2,2)
      SEST   = (SYY - SXY**2/SXX)/(DIFFMAT(2,2) - 2.D0)
      SEST   = DSQRT(SEST)
      SIGDFIT(1) = SEST/DSQRT(SXX)    !Error on the slope
      DFITA(1) =  ( DIFFMAT(2,2)*DVECA(1)  - DIFFMAT(1,2)*DVECA(2))/DET
      DFITA(2) =  (-DIFFMAT(2,1)*DVECA(1)  + DIFFMAT(1,1)*DVECA(2))/DET
      DFITB(1) =  ( DIFFMAT(2,2)*DVECB(1)  - DIFFMAT(1,2)*DVECB(2))/DET
      DFITB(2) =  (-DIFFMAT(2,1)*DVECB(1)  + DIFFMAT(1,1)*DVECB(2))/DET
      DFITAB(1)=  ( DIFFMAT(2,2)*DVECAB(1) - DIFFMAT(1,2)*DVECAB(2))/DET
      DFITAB(2)=  (-DIFFMAT(2,1)*DVECAB(1) + DIFFMAT(1,1)*DVECAB(2))/DET

      DFIT1(1)  =( DIFFMAT(2,2)*DVEC1(1)  -DIFFMAT(1,2)*DVEC1(2))/DET
      DFIT1(2)  =(-DIFFMAT(2,1)*DVEC1(1)  +DIFFMAT(1,1)*DVEC1(2))/DET
      SYY = DVEC1(3) - DVEC1(2)**2/DIFFMAT(2,2)
      SXY = DVEC1(1) - DVEC1(2)*DIFFMAT(1,2)/DIFFMAT(2,2)
      SEST   = (SYY - SXY**2/SXX)/(DIFFMAT(2,2) - 2.D0)
      SEST1  = DSQRT(SEST)
      SIGDFIT1(1) = SEST/DSQRT(SXX)
      DFITA1(1) =( DIFFMAT(2,2)*DVECA1(1) -DIFFMAT(1,2)*DVECA1(2))/DET
      DFITA1(2) =(-DIFFMAT(2,1)*DVECA1(1) +DIFFMAT(1,1)*DVECA1(2))/DET
      DFITB1(1) =( DIFFMAT(2,2)*DVECB1(1) -DIFFMAT(1,2)*DVECB1(2))/DET
      DFITB1(2) =(-DIFFMAT(2,1)*DVECB1(1) +DIFFMAT(1,1)*DVECB1(2))/DET
      DFITAB1(1)=(DIFFMAT(2,2)*DVECAB1(1)-DIFFMAT(1,2)*DVECAB1(2))/DET
      DFITAB1(2)=
     +         (-DIFFMAT(2,1)*DVECAB1(1)+DIFFMAT(1,1)*DVECAB1(2))/DET

      DFIT2(1)  =( DIFFMAT(2,2)*DVEC2(1)  -DIFFMAT(1,2)*DVEC2(2))/DET
      DFIT2(2)  =(-DIFFMAT(2,1)*DVEC2(1)  +DIFFMAT(1,1)*DVEC2(2))/DET
      SYY = DVEC2(3) - DVEC2(2)**2/DIFFMAT(2,2)
      SXY = DVEC2(1) - DVEC2(2)*DIFFMAT(1,2)/DIFFMAT(2,2)
      SEST   = (SYY - SXY**2/SXX)/(DIFFMAT(2,2) - 2.D0)
      SEST2  = DSQRT(SEST)
      SIGDFIT2(1) = SEST/DSQRT(SXX)
      DFITA2(1) =( DIFFMAT(2,2)*DVECA2(1) -DIFFMAT(1,2)*DVECA2(2))/DET
      DFITA2(2) =(-DIFFMAT(2,1)*DVECA2(1) +DIFFMAT(1,1)*DVECA2(2))/DET
      DFITB2(1) =( DIFFMAT(2,2)*DVECB2(1) -DIFFMAT(1,2)*DVECB2(2))/DET
      DFITB2(2) =(-DIFFMAT(2,1)*DVECB2(1) +DIFFMAT(1,1)*DVECB2(2))/DET
      DFITAB2(1)=(DIFFMAT(2,2)*DVECAB2(1)-DIFFMAT(1,2)*DVECAB2(2))/DET
      DFITAB2(2)=
     +         (-DIFFMAT(2,1)*DVECAB2(1)+DIFFMAT(1,1)*DVECAB2(2))/DET

      DFIT3(1)  =( DIFFMAT(2,2)*DVEC3(1)  -DIFFMAT(1,2)*DVEC3(2))/DET
      DFIT3(2)  =(-DIFFMAT(2,1)*DVEC3(1)  +DIFFMAT(1,1)*DVEC3(2))/DET
      SYY = DVEC3(3) - DVEC3(2)**2/DIFFMAT(2,2)
      SXY = DVEC3(1) - DVEC3(2)*DIFFMAT(1,2)/DIFFMAT(2,2)
      SEST   = (SYY - SXY**2/SXX)/(DIFFMAT(2,2) - 2.D0)
      SEST3  = DSQRT(SEST)
      SIGDFIT3(1) = SEST/DSQRT(SXX)
      DFITA3(1) =( DIFFMAT(2,2)*DVECA3(1) -DIFFMAT(1,2)*DVECA3(2))/DET
      DFITA3(2) =(-DIFFMAT(2,1)*DVECA3(1) +DIFFMAT(1,1)*DVECA3(2))/DET
      DFITB3(1) =( DIFFMAT(2,2)*DVECB3(1) -DIFFMAT(1,2)*DVECB3(2))/DET
      DFITB3(2) =(-DIFFMAT(2,1)*DVECB3(1) +DIFFMAT(1,1)*DVECB3(2))/DET
      DFITAB3(1)=(DIFFMAT(2,2)*DVECAB3(1)-DIFFMAT(1,2)*DVECAB3(2))/DET
      DFITAB3(2)=
     +          (-DIFFMAT(2,1)*DVECAB3(1)+DIFFMAT(1,1)*DVECAB3(2))/DET


C======Plot the diffusion data including any transients.

C      IF(E0.EQ.9.45D0)THEN
      IF(IE0.EQ.NDIFFDUMP)THEN
      DO 107 KF = 1,IDNT
      IF(KF.LT.1000.OR.MOD(KF,10).EQ.0)THEN !Dump 1st 1000, then every 10th.
      WRITE(61,'(I10,100E13.4)')KF,
     +                      SDIFFUSION(KF),   DFIT(1)*KF    +DFIT(2),
     +                      SDIFFUSIONA(KF),  DFITA(1)*KF   +DFITA(2),
     +                      SDIFFUSIONB(KF) , DFITB(1)*KF   +DFITB(2),
     +                      SDIFFUSIONAB(KF), DFITAB(1)*KF  +DFITAB(2),
     +                      SDIFFUSION1(KF),  DFIT1(1)*KF   +DFIT1(2),
     +                      SDIFFUSIONA1(KF), DFITA1(1)*KF  +DFITA1(2),
     +                      SDIFFUSIONB1(KF), DFITB1(1)*KF  +DFITB1(2),
     +                      SDIFFUSIONAB1(KF),DFITAB1(1)*KF +DFITAB1(2),
     +                      SDIFFUSION2(KF),  DFIT2(1)*KF   +DFIT2(2),
     +                      SDIFFUSIONA2(KF), DFITA2(1)*KF  +DFITA2(2),
     +                      SDIFFUSIONB2(KF), DFITB2(1)*KF  +DFITB2(2),
     +                      SDIFFUSIONAB2(KF),DFITAB2(1)*KF +DFITAB2(2),
     +                      SDIFFUSION3(KF),  DFIT3(1)*KF   +DFIT3(2),
     +                      SDIFFUSIONA3(KF), DFITA3(1)*KF  +DFITA3(2),
     +                      SDIFFUSIONB3(KF), DFITB3(1)*KF  +DFITB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),   
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),      !Padding
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      SDIFFUSIONAB3(KF),DFITAB3(1)*KF +DFITAB3(2),
     +                      (ODIFFUSION(IORB,KF),IORB=1,6)
      ENDIF

  107 CONTINUE
      ENDIF

C      IF(1.EQ.1)STOP




C=====Convert DFIT(1) to a depolarisation time. Correct for using enhanced 
C     radiation/damping: the rates are NDAMP2 times smaller than ``measured'' 
C     rates obtained with the ``used'' damping time.


      WRITE(53,103)
      WRITE(53,'(A,A,F13.5)')' ','Slope ratios (7,7)/(8,8): ',
     +                                               DFITA(1)/DFITB(1)

      IF(IE0.EQ.1)THEN       !For -ve slopes use the last sensible time. 
      TAUDMC(1)  = 0.D0        !Begin with zero.
      TAUDMC1(1) = 0.D0
      TAUDMC2(1) = 0.D0
      TAUDMC3(1) = 0.D0
      ENDIF

      RATEMC    = 0.5D0*   DFIT(1)/CTIME/NDAMP2
      SIGRATEMC = 0.5D0*SIGDFIT(1)/CTIME/NDAMP2
      IF(IE0.GT.1)TAUDMC(IE0) = TAUDMC(IE0-1)
      IF(RATEMC.GT.0.D0)TAUDMC(IE0) = 1.D0/RATEMC
      IF(RATEMC.LE.0.D0)RATEMC = 1.D-30
      WRITE(53,103)
      WRITE(53,'(A,A,4D15.4)')' Total:  fitted slope (+error), offset',
     +      ' and spread around fit  ', DFIT(1),SIGDFIT(1),DFIT(2),SEST
      WRITE(53,'(A,A,D15.4,A)')' ', 'M-C TAUD   ', TAUDMC(IE0),' SEC'

      IF(MODES.EQ.1)THEN
      RATEMC1    = 0.5D0*   DFIT1(1)/CTIME/NDAMP2
      SIGRATEMC1 = 0.5D0*SIGDFIT1(1)/CTIME/NDAMP2
      IF(IE0.GT.1)TAUDMC1(IE0) = TAUDMC1(IE0-1)
      IF(RATEMC1.GT.0.D0)TAUDMC1(IE0) = 1.D0/RATEMC1
      IF(RATEMC1.LE.0.D0)RATEMC1 = 1.D-30
      WRITE(53,'(A,A,4D15.4)')' Mode-1: fitted slope (+error), offset',
     + ' and spread around fit: ', DFIT1(1),SIGDFIT1(1),DFIT1(2),SEST1
      WRITE(53,'(A,A,D15.4,A)')' ', 'M-C TAUD1  ', TAUDMC1(IE0),' SEC'

      RATEMC2    = 0.5D0   *DFIT2(1)/CTIME/NDAMP2
      SIGRATEMC2 = 0.5D0*SIGDFIT2(1)/CTIME/NDAMP2
      IF(IE0.GT.1)TAUDMC2(IE0) = TAUDMC2(IE0-1)
      IF(RATEMC2.GT.0.D0)TAUDMC2(IE0) = 1.D0/RATEMC2
      IF(RATEMC2.LE.0.D0)RATEMC2 = 1.D-30
      WRITE(53,'(A,A,4D15.4)')' Mode-2: fitted slope (+error), offset',
     + ' and spread around fit: ', DFIT2(1),SIGDFIT2(1),DFIT2(2),SEST2 
      WRITE(53,'(A,A,D15.4,A)')' ', 'M-C TAUD2  ', TAUDMC2(IE0),' SEC'

      RATEMC3    = 0.5D0*   DFIT3(1)/CTIME/NDAMP2
      SIGRATEMC3 = 0.5D0*SIGDFIT3(1)/CTIME/NDAMP2
      IF(IE0.GT.1)TAUDMC3(IE0) = TAUDMC3(IE0-1)
      IF(RATEMC3.GT.0.D0)TAUDMC3(IE0) = 1.D0/RATEMC3
      IF(RATEMC3.LE.0.D0)RATEMC3 = 1.D-30
      WRITE(53,'(A,A,4D15.4)')' Mode-3: fitted slope (+error), offset',
     + ' and spread around fit: ', DFIT3(1),SIGDFIT3(1),DFIT3(2),SEST3  
      WRITE(53,'(A,A,D15.4,A)')' ', 'M-C TAUD3  ', TAUDMC3(IE0),' SEC'

      ENDIF

      F = 60.D0       !Convert to minutes
      WRITE(63,'(I10,20E16.6)')IE0,NU,
     +                          TAUDMC(IE0) /F,RATEMC *F,SIGRATEMC *F,
     +                          TAUDMC1(IE0)/F,RATEMC1*F,SIGRATEMC1*F, 
     +                          TAUDMC2(IE0)/F,RATEMC2*F,SIGRATEMC2*F, 
     +                          TAUDMC3(IE0)/F,RATEMC3*F,SIGRATEMC3*F,
     +                          E0 


      ENDIF

C======Get major and minor axes of 1-sigma transverse ellipse from COVTRACK. 
      IF(IBEQUIL.EQ.0.AND.IE0.EQ.0)THEN
      TILT=DATAN2(2.D0*COVTRACK(1,3,IDNT),COVTRACK(1,1,IDNT)
     +                                    -COVTRACK(3,3,IDNT)) ! 2 theta
      CTILT=DCOS(TILT)                                           !cos (2 theta)
      TILT=TILT*90.D0/PI
      XY3=COVTRACK(1,1,IDNT)*COVTRACK(3,3,IDNT)
     +                                      -COVTRACK(1,3,IDNT)**2
      IF(XY3.LT.0..AND.XY3.GT.-1.D-8)XY3=0.
      AXISMAJ   = DSQRT(XY3/ 
     +            ( 0.5D0*((COVTRACK(3,3,IDNT)-COVTRACK(1,1,IDNT))
     +             /CTILT + COVTRACK(3,3,IDNT)+COVTRACK(1,1,IDNT))))
      AXISMIN   = DSQRT(XY3/
     +            (-0.5D0*((COVTRACK(3,3,IDNT)-COVTRACK(1,1,IDNT))
     +             /CTILT - COVTRACK(3,3,IDNT)-COVTRACK(1,1,IDNT))))
C======Now transform the sigmas 
      DLTA  = TILT*PI/180.D0 
C======Plot the ellipse at the starting point. Take 1000 points. Also 
C      principle axes.
      DO 72 IPH = 1,1000
      PH    = 0.001D0*2.D0*PI*(IPH-1)
      HR    = AXISMAJ*DCOS(PH)    
      VR    = AXISMIN*DSIN(PH)    
      HOR   = HR*DCOS(DLTA) - VR*DSIN(DLTA)
      VER   = HR*DSIN(DLTA) + VR*DCOS(DLTA) 
      AXMAJ =    (-1.D0 + 0.002D0*(IPH-1))*DTAN(DLTA)
      AXMIN =    (-1.D0 + 0.002D0*(IPH-1))*DTAN(DLTA+PI/2.D0)
      WRITE(58,105)IPH,HOR,VER,-1.D0+0.002D0*(IPH-1),AXMAJ,AXMIN
   72 CONTINUE  
      ENDIF

C
C
C----------------------------------------------------------------------------
C======Get the approximate damping times: assumes no coupling, ignores 1st two 
C      damping times
      IF(NTURN2.GT.NSKIP)THEN
      DEC1 = COVTRACK(1,1,IDNT)/COVTRACK(1,1,3)
      DEC1 = -0.5D0*DLOG(DEC1)/(NDAMPTURNS*(IDNT - 2))
      DEC1 = CIR/(3.D8*DEC1)*1.D3
      DPER1 = DEC1*NDAMP2/TDAMP(1)
      DEC2 = COVTRACK(3,3,IDNT)/COVTRACK(3,3,3)
      DEC2 = -0.5D0*DLOG(DEC2)/(NDAMPTURNS*(IDNT - 2))
      DEC2 = CIR/(3.D8*DEC2)*1.D3
      DPER2 = DEC2*NDAMP2/TDAMP(3)
      DEC3 = COVTRACK(5,5,IDNT)/COVTRACK(5,5,3)
      DEC3 = -0.5D0*DLOG(DEC3)/(NDAMPTURNS*(IDNT - 2))
      DEC3 = CIR/(3.D8*DEC3)*1.D3
      DPER3 = DEC3*NDAMP2/TDAMP(5)
      ENDIF
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,3(1X,F11.1),A,3(1X,F11.1))')
     +                           ' Damping times (msec) ',
     +                            DEC1,DEC2,DEC3,  
     +                           '  Damping periods with NDAMP2 ',
     +                            DPER1,DPER2,DPER3



C===Do the FFT's
      IF(IE0.EQ.1)THEN
      WRITE(53,'(A,I6)')' KIFFT= ', KIFFT

      IFAIL = 0
C      CALL C06FAF(XFFT,KIFFT,XWFFT,IFAIL)
      IF(IFAIL.NE.0)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,A,A,I3)')
     +  ' In SCRURITA2: ',
     +  'the Fourier routine flagged a problem for the X motion.',
     +  ' So STOP. IFAIL= ',IFAIL
      STOP      
      ENDIF
      IFAIL = 0
C      CALL C06FAF(YFFT,KIFFT,YWFFT,IFAIL)
      IF(IFAIL.NE.0)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,A,A,I3)')
     +  ' In SCRURITA2: ',
     +  'the Fourier routine flagged a problem for the Y motion.',
     +  ' So STOP. IFAIL= ',IFAIL
      STOP      
      ENDIF

      AXFFT(0) = XFFT(0)
      BXFFT(0) = 0.D0
      XXFFT(0) = DSQRT(AXFFT(0)**2 + BXFFT(0)**2)
      AYFFT(0) = YFFT(0)
      BYFFT(0) = 0.D0
      YYFFT(0) = DSQRT(AYFFT(0)**2 + BYFFT(0)**2)
      HARM = 0.D0
      WRITE(67,'(3E16.6)') HARM/KIFFT,1.D-30,1.D-30
      WRITE(67,'(3E16.6)') HARM/KIFFT,XXFFT(0),YYFFT(0)
      WRITE(67,'(3E16.6)') HARM/KIFFT,1.D-30,1.D-30

      NFFT2 = (KIFFT-1)/2
      DO 60 J = 1, NFFT2
      NJ = KIFFT - J
      AXFFT(J)  =  XFFT(J)
      BXFFT(J)  =  XFFT(NJ)
      AYFFT(J)  =  YFFT(J)
      BYFFT(J)  =  YFFT(NJ)
      XXFFT(J)  = DSQRT(AXFFT(J)**2  + BXFFT(J)**2)
      YYFFT(J)  = DSQRT(AYFFT(J)**2  + BYFFT(J)**2)
      HARM = J
      WRITE(67,'(3E16.6)') HARM/KIFFT,1.D-30,1.D-30
      WRITE(67,'(3E16.6)') HARM/KIFFT,XXFFT(J),YYFFT(J)
      WRITE(67,'(3E16.6)') HARM/KIFFT,1.D-30,1.D-30
      AXFFT(NJ) =  XFFT(J)  
      BXFFT(NJ) = -XFFT(NJ) 
      AYFFT(NJ) =  YFFT(J)  
      BYFFT(NJ) = -YFFT(NJ) 
      HARM = NJ
      XXFFT(NJ) = DSQRT(AXFFT(NJ)**2 + BXFFT(NJ)**2)
      YYFFT(NJ) = DSQRT(AYFFT(NJ)**2 + BYFFT(NJ)**2)
C      WRITE(67,'(3E16.6)') HARM/KIFFT,1.D-30,1.D-30
C      WRITE(67,'(3E16.6)') HARM/KIFFT,XXFFT(NJ),YYFFT(NJ)
C      WRITE(67,'(3E16.6)') HARM/KIFFT,1.D-30,1.D-30
  60  CONTINUE
      IF(MOD(NFFT2,2).EQ.0) THEN
      AXFFT(NFFT2+1) = XFFT(NFFT2+1)
      BXFFT(NFFT2+1) = 0.D0
      XXFFT(NFFT2+1) = DSQRT(AXFFT(NFFT2+1)**2 + BXFFT(NFFT2+1)**2)
      AYFFT(NFFT2+1) = YFFT(NFFT2+1)
      BYFFT(NFFT2+1) = 0.D0
      YYFFT(NFFT2+1) = DSQRT(AYFFT(NFFT2+1)**2 + BYFFT(NFFT2+1)**2)
      HARM = NFFT2+1
      WRITE(67,'(3E16.6)') HARM/KIFFT,1.D-30,1.D-30      
      WRITE(67,'(3E16.6)') HARM/KIFFT,XXFFT(NFFT2+1),YYFFT(NFFT2+1)
      WRITE(67,'(3E16.6)') HARM/KIFFT,1.D-30,1.D-30      
      ENDIF

      ENDIF

      RETURN

 9999 WRITE(53,92)
   92 FORMAT(' ERROR IN EIGEN VALUE ROUTINE')
      STOP
      END

