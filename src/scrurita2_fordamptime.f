      SUBROUTINE SCRURITA2(IE0,E0,U0,CIR,COVMATIN,CRAD,TDAMP)
C
C
C
C
C   ROUTINE TO HANDLE SPIN-ORBIT TRACKING AND ESTIMATE THE RATE OF DEPOLARISATION.
C   ------------------------------------------------------------------------------
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

C
      PARAMETER (NPART   = 1000)  ! Maximum allowed particles. Use NPART2 of them.
      PARAMETER (LIMSECT = 1000)  ! Maximum allowed sections.
      PARAMETER (MAXD = 700)      ! Maximum number of damping times.
      PARAMETER (MAXE = 2000)     ! Maximum number of energy steps.
      REAL*8     NU
      DIMENSION ZR3(3,3),ZI3(3,3),A(3,3),B(3,3)
      DIMENSION ROT(3,3),TM3A(3,3)
      DIMENSION TM3B(3,3),WR3(3),WI3(3),ZW(3,3),ZWSAVE(3,3),TM3C(3,3)
      DIMENSION TRIN3(3,3),RR3(3),RI3(3),VR3(3,3),VI3(3,3),INTGE3(3)
      DIMENSION TW8A(8,8),TW2A(2,2)
      DIMENSION DISTMEAN(6),COVVEC(28),COVMATIN(6,6),COVMAT(6,6)
      DIMENSION COVSIM(8,8),COVTRACK(8,8,MAXD)
      DIMENSION COVSPIN1(2,2,MAXD),COVSPIN2(2,2,MAXD),COVSPIN3(2,2,MAXD)
      DIMENSION SORBVEC(8,NPART),SORBVET(NPART,8)
      DIMENSION ORBVEC1(6,NPART),ORBVEC2(6,NPART),ORBVEC3(6,NPART)
      DIMENSION SPINVEC1(2,NPART),SPINVEC2(2,NPART),SPINVEC3(2,NPART)
      DIMENSION SPINVET1(NPART,2),SPINVET2(NPART,2),SPINVET3(NPART,2)
      DIMENSION ONE(NPART),SORBAVE(8),SORBSIG(8)
      DIMENSION BB(8,8),AAA(8,8),SECTMAP(8,8,LIMSECT),PSCALE(LIMSECT)
      DIMENSION TREV8(8,8),TREV8D(8,8)
      DIMENSION BIGPHOT(NPART),TDAMP(6)
      DIMENSION SDIFFUSION(MAXD),    DFIT(2),    DVEC(2),   DIFFMAT(2,2)
      DIMENSION SDIFFUSIONA(MAXD),   DFITA(2),   DVECA(2)
      DIMENSION SDIFFUSIONB(MAXD),   DFITB(2),   DVECB(2)
      DIMENSION SDIFFUSIONAB(MAXD),  DFITAB(2),  DVECAB(2)

      DIMENSION SDIFFUSION1(MAXD),   DFIT1(2),   DVEC1(2)
      DIMENSION SDIFFUSIONA1(MAXD),  DFITA1(2),  DVECA1(2)
      DIMENSION SDIFFUSIONB1(MAXD),  DFITB1(2),  DVECB1(2)
      DIMENSION SDIFFUSIONAB1(MAXD), DFITAB1(2), DVECAB1(2)

      DIMENSION SDIFFUSION2(MAXD),   DFIT2(2),   DVEC2(2)
      DIMENSION SDIFFUSIONA2(MAXD),  DFITA2(2),  DVECA2(2)
      DIMENSION SDIFFUSIONB2(MAXD),  DFITB2(2),  DVECB2(2)
      DIMENSION SDIFFUSIONAB2(MAXD), DFITAB2(2), DVECAB2(2)

      DIMENSION SDIFFUSION3(MAXD),   DFIT3(2),   DVEC3(2)
      DIMENSION SDIFFUSIONA3(MAXD),  DFITA3(2),  DVECA3(2)
      DIMENSION SDIFFUSIONB3(MAXD),  DFITB3(2),  DVECB3(2)
      DIMENSION SDIFFUSIONAB3(MAXD), DFITAB3(2), DVECAB3(2)

      DIMENSION TAUDMC(MAXE),TAUDMC1(MAXE),TAUDMC2(MAXE),TAUDMC3(MAXE)



      DIMENSION ODIFFUSION(6,MAXD)

      DIMENSION KGBINX(-500:500),KGBINY(-500:500)
      DIMENSION SORBDUM(8)
C=====8x8 eigenvector stuff.
      DIMENSION TRIN8(8,8),RR8(8),RI8(8),VR8(8,8),VI8(8,8)
      DIMENSION INTGE8(8),WW8(8,8)
      DIMENSION ZZ(8,8),ZV(8,8),WR8(8),WI8(8),ZZSAVE(8,8),USYMP(6,6)
      DIMENSION TM8A(8,8),TM8B(8,8)
      DIMENSION TN(8),PTUN(6)
      DIMENSION AB(6)
      DIMENSION TEMP(6,NPART),ORBAMP(6,NPART),GMAT(2,6),SIXVEC(6,NPART)
      DIMENSION TESTVEC(6)
      DIMENSION ZZSECT(8,8,LIMSECT),ZZCONJ(6,6),ORBAMP(6,NPART)


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
      NU=E0/0.440652D0 * 1.0D0       ! Can scale a gamma indep. of the energy.


      WRITE(53,103)
      WRITE(53,103)
  103 FORMAT(/,'  ')


      WRITE(53,929)
  929 FORMAT('1','Entering subroutine SCRURITA2 for M-C tracking to ``me
     +asure'' the rate of depolarisation...')
C      write(*,929)

      IF(NPART2.GT.NPART.OR.NTURN2.GT.MAXD.OR.NSTEP.GT.MAXE)THEN
      WRITE(53,103)
      WRITE(53,'(3A)')' STOP: attempt to work with more than NPART',
     +                 '  particles or more than MAXD damping times',
     +                 '  or more than MAXE energy steps.'
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
      CALL F02AGF(TRIN3,3,3,RR3,RI3,VR3,3,VI3,3,INTGE3,IFAIL)
      IF(IFAIL.NE.0)GO TO 9999
  926 FORMAT(' ','NOW THE F02AGF RESULTS: 3X3')
C      DO 284 I=1,3
C      WR3(I)=RR3(I)
C      WI3(I)=RI3(I)
C      DO 285 J=1,3
C      ZR3(J,I)=VR3(J,I)
C      ZI3(J,I)=VI3(J,I)
C  285 CONTINUE
C  284 CONTINUE
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

      CALL G05CBF(15)              !No freedom to choose the starting seed yet.
                                    !Original = 15
      DISTMEAN = 0.D0
      EPS = 1.D-12
      IFAIL1 = 0
      CALL G05EAF(DISTMEAN,6,COVMATIN,6,EPS,COVVEC,28,IFAIL1)
      IF(IFAIL1.NE.0)THEN
      WRITE(53,'(A,A,I10,A)')' ', 'IFAIL1 = ',IFAIL1,' STOP'  
      ENDIF

      WRITE(53,103)
      WRITE(53,103)



      SORBVEC  = 0.D0     ! Set to zero in case G05EZF is switched off.
      SPINVEC1 = 0.D0
      SPINVEC2 = 0.D0      
      SPINVEC3 = 0.D0

      DO 1  I = 1, NPART2
      IFAIL2 = 0 
C      CALL G05EZF(SORBVEC(1:6,I),6,COVVEC,28,IFAIL2)
      IF(IFAIL2.NE.0)THEN
      WRITE(53,'(A,A,I10,A)')' ', 'IFAIL2 = ',IFAIL2,' STOP'  
      ENDIF
    1 CONTINUE

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
      IF(IID.EQ.1.AND.XY.LT.0.00011D0)GO TO 8   !Skip zero length drifts
      IF(IID.EQ.1)THEN                          !Treat drifts separately.
      TREV8(1,:) = TREV8(1,:) + XY * TREV8(2,:) 
      TREV8(3,:) = TREV8(3,:) + XY * TREV8(4,:) 
      GO TO 8
      ELSEIF(IID.EQ.10)THEN
      ZW = TM3C
      CALL SOL8AN(II,XY,YYY,NU,ZW,TM3C,BB,NSOL(ITY)) 
      ELSE
      BB(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 
      CALL MX88DAMP(BB,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NU,TM3B,ZW,
     +                                                 CRAD,IDAMPFLG)
      ENDIF
      TREV8 = MATMUL(BB,TREV8)
    8 CONTINUE

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
C      write(*,'(A)')
C     +   ' 1-turn spin-(symplectic)orbit matrix, generated afresh.'
      DO 36 I=1,8
   36 WRITE(53,914)(TREV8(I,J),J=1,8)

C=====GET 8X8 EIGENVECTORS FOR ONE REVOLUTION & ORDER THEM==============
      TRIN8 = TREV8
      IFAIL=0
      CALL F02AGF(TRIN8,8,8,RR8,RI8,VR8,8,VI8,8,INTGE8,IFAIL)
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
C      write(*,'(A)')' Reached here 1'
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

C      write(*,'(A)')' Reached here 2'

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
      CALL SOL8AN(II,XY,YYY,NU,ZW,TM3C,BB,NSOL(ITY))
      ELSE
      BB(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 
      CALL MX88DAMP(BB,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NU,TM3B,ZW,
     +                                                 CRAD,IDAMPFLG)
      ENDIF
      TREV8D = MATMUL(BB,TREV8D)
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
      WRITE(53,'(A)')' 1-turn spin-(damped)orbit matrix.'
C      write(*,'(A)') ' 1-turn spin-(damped)orbit matrix.'
      DO 336 I=1,8
  336 WRITE(53,914)(TREV8D(I,J),J=1,8)



C-------------------------------------------------------------------------------------
C======Make 1 pass around the ring with damping to generate damped
C      spin-orbit matrices for sections between the dipole centres using
C      the fact that the C.O. is fixed.  
C      So don't have to recalc all fixed stuff on each turn and can track 
C      using the sections.
C      Get the radiation scale factors which account for the various dipole  
C      strengths.
      NSECT = 0
      IDAMPFLG = 1                               !Switch on/off radn.
      CALL UNIT(8,AAA)
      TM3C = ZWSAVE
      RADTOT = 0.D0

      DO 7  II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
      IF(IID.EQ.10)THEN
      ZW = TM3C
      CALL SOL8AN(II,XY,YYY,NU,ZW,TM3C,BB,NSOL(ITY))
      ELSE
      BB(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 
      CALL MX88DAMP(BB,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NU,TM3B,ZW,
     +                                                 CRAD,IDAMPFLG)
      ENDIF
      AAA = MATMUL(BB,AAA)
      IF((IID.EQ.7.AND.NAME(ITY)(1:2).EQ.'VD').OR.II.EQ.NELEM)THEN  
      NSECT = NSECT + 1

      IF(NSECT.GT.LIMSECT)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')' STOP: the number of sections exceeds LIMSECT'
      STOP
      ENDIF

      SECTMAP(:,:,NSECT)   = AAA
      PSCALE(NSECT) = RADSCALE(ITY)
      IF(IE0.EQ.1) WRITE(53,'(A,2I10,2A,E16.5)')
     +                          ' ',II,IID,'  ',NAME(ITY),PSCALE(NSECT)

      RADTOT = RADTOT + PSCALE(NSECT)
      CALL UNIT(8,AAA)      
      ENDIF

    7 CONTINUE

      WRITE(53,'(A,E16.5,A,E16.5)')' Total radiation excitation: ',
     +                                                       RADTOT
C      write(*,'(A)')' Reached here 3'



      RADAVE = RADTOT/(NSECT-1)
      PSCALE(1:NSECT) = PSCALE(1:NSECT)/RADAVE  ! Scale factors for individual dipoles.

C      PSCALE(1) = 0.D0
C      PSCALE(2) = 0.D0
C      PSCALE(3) = 0.D0
C      PSCALE(176) = 0.D0
C      PSCALE(177) = 0.D0
C      PSCALE(178) = 0.D0


       
C      write(*,'(A)')' Reached here 4'

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
C      write(*,'(A,A,A,I5)')' 1-turn spin-(damped)orbit matrix',
C     +                      ' using sections.','   NDAMP2= ',NDAMP2
      WRITE(53,'(A,A,A,I5)')' 1-turn spin-(damped)orbit matrix',
     +                      ' using sections.','   NDAMP2= ',NDAMP2
      DO 35 I=1,8
   35 WRITE(53,914)(AAA(I,J),J=1,8)
  914 FORMAT(8F12.5)



C-------------------------------------------------------------------------------------
C======Make another pass around the ring with NO damping to generate s-o eigenvectors
C      at the start of each section.
C      So don't have to recalc all fixed eigenvector on each turn.
C      This is basically a repeat of the above stuff on sections above but 
C      without the damping.
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
      IF(IID.EQ.10)THEN
      ZW = TM3C
      CALL SOL8AN(II,XY,YYY,NU,ZW,TM3C,BB,NSOL(ITY))
      ELSE
      BB(1:6,1:6)  = TMAT(1:6,1:6,ITY)
      ZW = TM3C
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      TM3C = MATMUL(TM3A,ZW) 
      TM3B = 0.5D0*(ZW + TM3C) 
      CALL MX88DAMP(BB,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),NU,TM3B,ZW,
     +                                                 CRAD,IDAMPFLG)
      ENDIF
      AAA = MATMUL(BB,AAA)
      IF((IID.EQ.7.AND.NAME(ITY)(1:2).EQ.'VD').OR.II.EQ.NELEM)THEN  
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

C      WRITE(53,'(A,A,F20.8)')' ','CRAD',CRADIN

C======Now track to get the equilibrium beam and the non-equilibrium spin 
C      components.
C======After each dipole or C-F, radiate ``big'' photons according to a 
C      centred top-hat distribution. 
C      In contrast to a Gaussian, there are no long tails so that no cut 
C      against too large (particle expelling kicks) is needed.
C      For linear motion this should also be enough to get Gaussian phase space
C      by the central limit theorem. Get a vector of NPART2 samples.  
C      Use the current seed.
C      Use the Boege formula for the strength. The rms ``big'' photon energy is
C      very much smaller than the energy spread of the beam.
C      All ``big'' photons have the same energy spread.
C      Later,scale the strengths according to the curvature and length of the 
C      dipoles
C      This could be important for dipoles where n_0 is horizontal and there 
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


      SORBVEC  = SORBVEC*1.D-3
      SPINVEC1 = SPINVEC1*1.D-3
      SPINVEC2 = SPINVEC2*1.D-3
      SPINVEC3 = SPINVEC3*1.D-3

      IF(IE0.EQ.1)THEN
      COVMAT  = COVMATIN *1.D-6
      TDAMP   = TDAMP*1.D-3
      ENDIF      

      CTIME =  CIR/3.0D8
      DAMPTURNS = TDAMP(5)/CTIME/NDAMP2
      WRITE(53,'(A,F8.1,A,I3)')' Turns per sync. damping time', 
     +                       DAMPTURNS,'  NDAMP2=', NDAMP2
      NDAMPTURNS = DAMPTURNS   
      IDAMPFLG = 1                              !Switch on/off radn.
      TH = 0.5D0
      VARTH = TH*TH/3.D0
      PHOTSCALE = CTIME*COVMAT(6,6)/(TDAMP(5)*VARTH*(NSECT-1))
      PHOTSCALE = 2.D0*DSQRT(PHOTSCALE*NDAMP2)*IDAMPFLG
      TH = TH * PHOTSCALE
      WRITE(53,'(A,A,I6,E15.7)')' ','NSECT,PHOTSCALE', NSECT,PHOTSCALE
      WRITE(53,103)


      CALL UNIT(2,TW2A)
      TW2A(1,1)= WR3(2)
      TW2A(1,2)= WI3(2)
      TW2A(2,1)=-WI3(2)
      TW2A(2,2)= WR3(2)

      IDNT = 0
      DO 9 INT = 1,NTURN2               !Loop over damping periods.
      DO 4 IDT = 1,NDAMPTURNS           !Loop over turns in a damping period. 
      DO 5 II  = 1,NSECT                !Loop over sections in a turn.

C=====Use SORBVEC to get the vectors corresponding to the three orbital modes.
C     Do this using the ``real'' forms for the eigenvectors, i.e., don't
C     convert to complex form. Then must be careful with identifying the  
C     correct pairs of real vectors and with signs. 
C     Check by rebuilding the original vector.
      ZZCONJ(1:6,1:6)  =  ZZSECT(1:6,1:6,II)
      ZZCONJ           =  TRANSPOSE(ZZCONJ)
      TEMP(1,1:NPART2) = -SORBVEC(2,1:NPART2)
      TEMP(2,1:NPART2) =  SORBVEC(1,1:NPART2)
      TEMP(3,1:NPART2) = -SORBVEC(4,1:NPART2)
      TEMP(4,1:NPART2) =  SORBVEC(3,1:NPART2)
      TEMP(5,1:NPART2) = -SORBVEC(6,1:NPART2)
      TEMP(6,1:NPART2) =  SORBVEC(5,1:NPART2)
      ORBAMP(:,1:NPART2) = MATMUL(ZZCONJ,TEMP(:,1:NPART2))


C=====Now reconstruct the 3 real orbit vectors. 
      DO 6 INP = 1,NPART2
      ORBVEC1(:,INP) = 2.D0*
     + (ORBAMP(2,INP)*ZZSECT(1:6,1,II) - ORBAMP(1,INP)*ZZSECT(1:6,2,II))
      ORBVEC2(:,INP) = 2.D0*
     + (ORBAMP(4,INP)*ZZSECT(1:6,3,II) - ORBAMP(3,INP)*ZZSECT(1:6,4,II))
      ORBVEC3(:,INP) = 2.D0*
     + (ORBAMP(6,INP)*ZZSECT(1:6,5,II) - ORBAMP(5,INP)*ZZSECT(1:6,6,II))
   6  CONTINUE  
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

C=====Now transport the various vectors.
      SORBVEC(:,1:NPART2) = MATMUL(SECTMAP(:,:,II),SORBVEC(:,1:NPART2))

      GMAT = SECTMAP(7:8,1:6,II)
      SPINVEC1(:,1:NPART2) = MATMUL(GMAT,ORBVEC1(:,1:NPART2))
     +                                      + SPINVEC1(:,1:NPART2)
      SPINVEC2(:,1:NPART2) = MATMUL(GMAT,ORBVEC2(:,1:NPART2))
     +                                      + SPINVEC2(:,1:NPART2)
      SPINVEC3(:,1:NPART2) = MATMUL(GMAT,ORBVEC3(:,1:NPART2))
     +                                      + SPINVEC3(:,1:NPART2)      


C     ORBVEC1(:,1:NPART2)
C    +           = MATMUL(SECTMAP(1:6,1:6,II),ORBVEC1(:,1:NPART2))
C     ORBVEC2(:,1:NPART2)
C    +           = MATMUL(SECTMAP(1:6,1:6,II),ORBVEC2(:,1:NPART2))
C     ORBVEC3(:,1:NPART2)
C    +           = MATMUL(SECTMAP(1:6,1:6,II),ORBVEC3(:,1:NPART2))



C      write(*,*)II,IDT,INT

C=====Radiate. But only at dipoles, not at the IP. Weight according to the
C     dipole strength.
C      IF(II.NE.NSECT)THEN
C      WRITE(53,'(A,I10,E16.5)')' PSCALE:  ',II,PSCALE(II)
      CALL G05FAF(-TH,TH,NPART2,BIGPHOT)
      DO 55 IPHS = 1,NPART2
      SORBVEC(6,IPHS) = SORBVEC(6,IPHS) + BIGPHOT(IPHS)               !*PSCALE(II)
   55 CONTINUE
C      ENDIF

    5 CONTINUE

C=====WIND BACK THE SPIN BASIS AT THE END OF A TURN==========================
      SORBVEC(7:8,1:NPART2)  = MATMUL(TW2A,SORBVEC(7:8,1:NPART2))
      SPINVEC1(:,1:NPART2)   = MATMUL(TW2A,SPINVEC1(:,1:NPART2))
      SPINVEC2(:,1:NPART2)   = MATMUL(TW2A,SPINVEC2(:,1:NPART2))
      SPINVEC3(:,1:NPART2)   = MATMUL(TW2A,SPINVEC3(:,1:NPART2))



    4 CONTINUE   
C======At the end of each sync damping time, get the covariance matrix 
C      And write out some sample vectors.

C     DO 666 INP = 1,NPART2
C     IF(INP.LE.5)THEN
C     WRITE(53,'(A,A,I5,8F12.4)')' ',
C    + 'Transported spin-orbit vector (mm,mrad..):    all',
C    +                                        INP,SORBVEC(1:8,INP)*1.D3
C     WRITE(53,'(A,A,I5,8F12.4)')' ',
C    + 'Transported spin-orbit vector (mm,mrad..): mode-1',
C    +                  INP,ORBVEC1(1:6,INP)*1.D3,SPINVEC1(:,INP)*1.D3
C     WRITE(53,'(A,A,I5,8F12.4)')' ',
C    + 'Transported spin-orbit vector (mm,mrad..): mode-2',
C    +                  INP,ORBVEC2(1:6,INP)*1.D3,SPINVEC2(:,INP)*1.D3
C     WRITE(53,'(A,A,I5,8F12.4)')' ',
C    + 'Transported spin-orbit vector (mm,mrad..): mode-3',
c    +                  INP,ORBVEC3(1:6,INP)*1.D3,SPINVEC3(:,INP)*1.D3
C     ENDIF
C 666 CONTINUE  


      SORBVET(1:NPART2,:) = TRANSPOSE(SORBVEC(:,1:NPART2))
      SPINVET1(1:NPART2,:)= TRANSPOSE(SPINVEC1(:,1:NPART2))
      SPINVET2(1:NPART2,:)= TRANSPOSE(SPINVEC2(:,1:NPART2))
      SPINVET3(1:NPART2,:)= TRANSPOSE(SPINVEC3(:,1:NPART2))
      COVTRACK(:,:,INT) 
     +     = MATMUL(SORBVEC(:,1:NPART2),SORBVET(1:NPART2,:))/NPART2
      COVSPIN1(:,:,INT)
     +     = MATMUL(SPINVEC1(:,1:NPART2),SPINVET1(1:NPART2,:))/NPART2
      COVSPIN2(:,:,INT)
     +     = MATMUL(SPINVEC2(:,1:NPART2),SPINVET2(1:NPART2,:))/NPART2
      COVSPIN3(:,:,INT)
     +     = MATMUL(SPINVEC3(:,1:NPART2),SPINVET3(1:NPART2,:))/NPART2

      WRITE(53,103)
      WRITE(53,86)INT,INT*IDT
   86 FORMAT(//,' The ',I3,'th  transported s-o covariance matrix',
     +          ' <Xi*Xj>(mm*mm,mm*mrad...)',
     +          ' at the 1-st beam line element after ',I6,' turns')
      WRITE(53,955)COVTRACK(:,:,INT)*1.D6
      WRITE(53,87)INT,INT*IDT
   87 FORMAT(/,' The ',I3,'th  transported spin covariance matrix',
     +          ' mode-1 <Si*Sj>(mrad*mrad...)',
     +          ' at the 1-st beam line element after ',I6,' turns')
      WRITE(53,956)COVSPIN1(:,:,INT)*1.D6
  956 FORMAT(T6,2F16.3)
      WRITE(53,88)INT,INT*IDT
   88 FORMAT(/,' The ',I3,'th  transported spin covariance matrix',
     +          ' mode-2 <Si*Sj>(mrad*mrad...)',
     +          ' at the 1-st beam line element after ',I6,' turns')
      WRITE(53,956)COVSPIN2(:,:,INT)*1.D6
      WRITE(53,89)INT,INT*IDT
   89 FORMAT(/,' The ',I3,'th  transported spin covariance matrix',
     +          ' mode-3 <Si*Sj>(mrad*mrad...)',
     +          ' at the 1-st beam line element after ',I6,' turns')
      WRITE(53,956)COVSPIN3(:,:,INT)*1.D6

      WRITE(53,103)
      SDIFFUSION(INT)     = COVTRACK(7,7,INT)+COVTRACK(8,8,INT)
      SDIFFUSIONA(INT)    = COVTRACK(7,7,INT)
      SDIFFUSIONB(INT)    = COVTRACK(8,8,INT)
      SDIFFUSIONAB(INT)   = COVTRACK(7,8,INT)

      SDIFFUSION1(INT)    = COVSPIN1(1,1,INT)+COVSPIN1(2,2,INT)
      SDIFFUSIONA1(INT)   = COVSPIN1(1,1,INT)
      SDIFFUSIONB1(INT)   = COVSPIN1(2,2,INT)
      SDIFFUSIONAB1(INT)  = COVSPIN1(1,2,INT)

      SDIFFUSION2(INT)    = COVSPIN2(1,1,INT)+COVSPIN2(2,2,INT)
      SDIFFUSIONA2(INT)   = COVSPIN2(1,1,INT)
      SDIFFUSIONB2(INT)   = COVSPIN2(2,2,INT)
      SDIFFUSIONAB2(INT)  = COVSPIN2(1,2,INT)

      SDIFFUSION3(INT)    = COVSPIN3(1,1,INT)+COVSPIN3(2,2,INT)
      SDIFFUSIONA3(INT)   = COVSPIN3(1,1,INT)
      SDIFFUSIONB3(INT)   = COVSPIN3(2,2,INT)
      SDIFFUSIONAB3(INT)  = COVSPIN3(1,2,INT)

      DO 90 IK = 1,6
   90 ODIFFUSION(IK,INT) = COVTRACK(IK,IK,INT)


      WRITE(53,'(A,A,F13.5,A,F13.5)')' ','(7,7) + (8,8): ',
     +                             SDIFFUSION(INT),'    (7,7)/(8,8): ',
     +                         COVTRACK(7,7,INT)/COVTRACK(8,8,INT)

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
   94 SORBSIG(IJ)=DSQRT(COVTRACK(IJ,IJ,NTURN2))
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
C======Plot scatter of the spin angles alpha and beta. 
      DO 104 IOO = 1, NPART2 
  104 WRITE(60,105)IOO,SORBVEC(7,IOO),SORBVEC(8,IOO)
  105 FORMAT(' ',I10, 6E16.6) 

C======Plot the rms diffusion angle vs. damping times and do a LSQ fit.
C      Start at the 3rd data point to exclude transients.
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

      DO 106 JF = 4, NTURN2
      DIFFMAT(1,1) = DIFFMAT(1,1) + JF*JF 
      DIFFMAT(1,2) = DIFFMAT(1,2) + JF
      DIFFMAT(2,2) = DIFFMAT(2,2) + 1.D0
      DVEC(1)   = DVEC(1)   + JF*SDIFFUSION(JF)
      DVEC(2)   = DVEC(2)   +    SDIFFUSION(JF)
      DVECA(1)  = DVECA(1)  + JF*SDIFFUSIONA(JF)
      DVECA(2)  = DVECA(2)  +    SDIFFUSIONA(JF)
      DVECB(1)  = DVECB(1)  + JF*SDIFFUSIONB(JF)
      DVECB(2)  = DVECB(2)  +    SDIFFUSIONB(JF)
      DVECAB(1) = DVECAB(1) + JF*SDIFFUSIONAB(JF)
      DVECAB(2) = DVECAB(2) +    SDIFFUSIONAB(JF)

      DVEC1(1)   = DVEC1(1)   + JF*SDIFFUSION1(JF)
      DVEC1(2)   = DVEC1(2)   +    SDIFFUSION1(JF)
      DVECA1(1)  = DVECA1(1)  + JF*SDIFFUSIONA1(JF)
      DVECA1(2)  = DVECA1(2)  +    SDIFFUSIONA1(JF)
      DVECB1(1)  = DVECB1(1)  + JF*SDIFFUSIONB1(JF)
      DVECB1(2)  = DVECB1(2)  +    SDIFFUSIONB1(JF)
      DVECAB1(1) = DVECAB1(1) + JF*SDIFFUSIONAB1(JF)
      DVECAB1(2) = DVECAB1(2) +    SDIFFUSIONAB1(JF)

      DVEC2(1)   = DVEC2(1)   + JF*SDIFFUSION2(JF)
      DVEC2(2)   = DVEC2(2)   +    SDIFFUSION2(JF)
      DVECA2(1)  = DVECA2(1)  + JF*SDIFFUSIONA2(JF)
      DVECA2(2)  = DVECA2(2)  +    SDIFFUSIONA2(JF)
      DVECB2(1)  = DVECB2(1)  + JF*SDIFFUSIONB2(JF)
      DVECB2(2)  = DVECB2(2)  +    SDIFFUSIONB2(JF)
      DVECAB2(1) = DVECAB2(1) + JF*SDIFFUSIONAB2(JF)
      DVECAB2(2) = DVECAB2(2) +    SDIFFUSIONAB2(JF)

      DVEC3(1)   = DVEC3(1)   + JF*SDIFFUSION3(JF)
      DVEC3(2)   = DVEC3(2)   +    SDIFFUSION3(JF)
      DVECA3(1)  = DVECA3(1)  + JF*SDIFFUSIONA3(JF)
      DVECA3(2)  = DVECA3(2)  +    SDIFFUSIONA3(JF)
      DVECB3(1)  = DVECB3(1)  + JF*SDIFFUSIONB3(JF)
      DVECB3(2)  = DVECB3(2)  +    SDIFFUSIONB3(JF)
      DVECAB3(1) = DVECAB3(1) + JF*SDIFFUSIONAB3(JF)
      DVECAB3(2) = DVECAB3(2) +    SDIFFUSIONAB3(JF)
  106 CONTINUE 
      DIFFMAT(2,1) = DIFFMAT(1,2)
      DET = DIFFMAT(1,1)*DIFFMAT(2,2) - DIFFMAT(1,2)**2
      DFIT(1)  =  ( DIFFMAT(2,2)*DVEC(1)   - DIFFMAT(1,2)*DVEC(2))/DET
      DFIT(2)  =  (-DIFFMAT(2,1)*DVEC(1)   + DIFFMAT(1,1)*DVEC(2))/DET
      DFITA(1) =  ( DIFFMAT(2,2)*DVECA(1)  - DIFFMAT(1,2)*DVECA(2))/DET
      DFITA(2) =  (-DIFFMAT(2,1)*DVECA(1)  + DIFFMAT(1,1)*DVECA(2))/DET
      DFITB(1) =  ( DIFFMAT(2,2)*DVECB(1)  - DIFFMAT(1,2)*DVECB(2))/DET
      DFITB(2) =  (-DIFFMAT(2,1)*DVECB(1)  + DIFFMAT(1,1)*DVECB(2))/DET
      DFITAB(1)=  ( DIFFMAT(2,2)*DVECAB(1) - DIFFMAT(1,2)*DVECAB(2))/DET
      DFITAB(2)=  (-DIFFMAT(2,1)*DVECAB(1) + DIFFMAT(1,1)*DVECAB(2))/DET

      DFIT1(1)  =( DIFFMAT(2,2)*DVEC1(1)  -DIFFMAT(1,2)*DVEC1(2))/DET
      DFIT1(2)  =(-DIFFMAT(2,1)*DVEC1(1)  +DIFFMAT(1,1)*DVEC1(2))/DET
      DFITA1(1) =( DIFFMAT(2,2)*DVECA1(1) -DIFFMAT(1,2)*DVECA1(2))/DET
      DFITA1(2) =(-DIFFMAT(2,1)*DVECA1(1) +DIFFMAT(1,1)*DVECA1(2))/DET
      DFITB1(1) =( DIFFMAT(2,2)*DVECB1(1) -DIFFMAT(1,2)*DVECB1(2))/DET
      DFITB1(2) =(-DIFFMAT(2,1)*DVECB1(1) +DIFFMAT(1,1)*DVECB1(2))/DET
      DFITAB1(1)=(DIFFMAT(2,2)*DVECAB1(1)-DIFFMAT(1,2)*DVECAB1(2))/DET
      DFITAB1(2)=
     +         (-DIFFMAT(2,1)*DVECAB1(1)+DIFFMAT(1,1)*DVECAB1(2))/DET

      DFIT2(1)  =( DIFFMAT(2,2)*DVEC2(1)  -DIFFMAT(1,2)*DVEC2(2))/DET
      DFIT2(2)  =(-DIFFMAT(2,1)*DVEC2(1)  +DIFFMAT(1,1)*DVEC2(2))/DET
      DFITA2(1) =( DIFFMAT(2,2)*DVECA2(1) -DIFFMAT(1,2)*DVECA2(2))/DET
      DFITA2(2) =(-DIFFMAT(2,1)*DVECA2(1) +DIFFMAT(1,1)*DVECA2(2))/DET
      DFITB2(1) =( DIFFMAT(2,2)*DVECB2(1) -DIFFMAT(1,2)*DVECB2(2))/DET
      DFITB2(2) =(-DIFFMAT(2,1)*DVECB2(1) +DIFFMAT(1,1)*DVECB2(2))/DET
      DFITAB2(1)=(DIFFMAT(2,2)*DVECAB2(1)-DIFFMAT(1,2)*DVECAB2(2))/DET
      DFITAB2(2)=
     +         (-DIFFMAT(2,1)*DVECAB2(1)+DIFFMAT(1,1)*DVECAB2(2))/DET

      DFIT3(1)  =( DIFFMAT(2,2)*DVEC3(1)  -DIFFMAT(1,2)*DVEC3(2))/DET
      DFIT3(2)  =(-DIFFMAT(2,1)*DVEC3(1)  +DIFFMAT(1,1)*DVEC3(2))/DET
      DFITA3(1) =( DIFFMAT(2,2)*DVECA3(1) -DIFFMAT(1,2)*DVECA3(2))/DET
      DFITA3(2) =(-DIFFMAT(2,1)*DVECA3(1) +DIFFMAT(1,1)*DVECA3(2))/DET
      DFITB3(1) =( DIFFMAT(2,2)*DVECB3(1) -DIFFMAT(1,2)*DVECB3(2))/DET
      DFITB3(2) =(-DIFFMAT(2,1)*DVECB3(1) +DIFFMAT(1,1)*DVECB3(2))/DET
      DFITAB3(1)=(DIFFMAT(2,2)*DVECAB3(1)-DIFFMAT(1,2)*DVECAB3(2))/DET
      DFITAB3(2)=
     +          (-DIFFMAT(2,1)*DVECAB3(1)+DIFFMAT(1,1)*DVECAB3(2))/DET


C======Plot the diffusion data beyond any transients.

      IF(E0.EQ.9.45D0)THEN
      DO 107 KF = 1,NTURN2
      WRITE(61,'(I10,40E13.4)')KF,
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
     +                      (ODIFFUSION(IORB,KF),IORB=1,6)

  107 CONTINUE
      ENDIF

C=====Convert DFIT(1) to a depolarisation time. Correct for using enhanced 
C     radiation/damping.
      WRITE(53,103)
      WRITE(53,'(A,A,F13.5)')' ','Slope ratios (7,7)/(8,8): ',
     +                                               DFITA(1)/DFITB(1)

C=====Get the depoln rates and then times. The rates are NDAMP2 times 
C     smaller than ``measured'' rates obtained with the ``used'' damping time. So in the 
C     end one just needs to divide by the real damping time.

      IF(IE0.EQ.1)THEN       !For -ve slopes use the last sensible time. 
      TAUDMC(1)  = 0.D0        !Begin with zero.
      TAUDMC1(1) = 0.D0
      TAUDMC2(1) = 0.D0
      TAUDMC3(1) = 0.D0
      ENDIF


      RATEMC = 0.5D0*DFIT(1)/TDAMP(5)
      IF(IE0.GT.1)TAUDMC(IE0) = TAUDMC(IE0-1)
      IF(RATEMC.GT.0.D0)TAUDMC(IE0) = 1.D0/RATEMC
      WRITE(53,103)
      WRITE(53,'(A,A,2D15.8,A)')' ', ' Fitted slope and offset:      ',
     +                                                              DFIT 
      WRITE(53,'(A,A,D15.4,A)')' ', 'M-C TAUD   ', TAUDMC(IE0),' SEC'

      RATEMC1 = 0.5D0*DFIT1(1)/TDAMP(5)
      IF(IE0.GT.1)TAUDMC1(IE0) = TAUDMC1(IE0-1)
      IF(RATEMC1.GT.0.D0)TAUDMC1(IE0) = 1.D0/RATEMC1
      WRITE(53,'(A,A,2D15.8,A)')' ', 
     +                     'Mode-1: fitted slope and offset: ',DFIT1 
      WRITE(53,'(A,A,D15.4,A)')' ', 'M-C TAUD1  ', TAUDMC1(IE0),' SEC'

      RATEMC2 = 0.5D0*DFIT2(1)/TDAMP(5)
      IF(IE0.GT.1)TAUDMC2(IE0) = TAUDMC2(IE0-1)
      IF(RATEMC2.GT.0.D0)TAUDMC2(IE0) = 1.D0/RATEMC2
      WRITE(53,'(A,A,2D15.8,A)')' ', 
     +                     'Mode-2: fitted slope and offset: ',DFIT2 
      WRITE(53,'(A,A,D15.4,A)')' ', 'M-C TAUD2  ', TAUDMC2(IE0),' SEC'

      RATEMC3 = 0.5D0*DFIT3(1)/TDAMP(5)
      IF(IE0.GT.1)TAUDMC3(IE0) = TAUDMC3(IE0-1)
      IF(RATEMC3.GT.0.D0)TAUDMC3(IE0) = 1.D0/RATEMC3
      WRITE(53,'(A,A,2D15.8,A)')' ', 
     +                     'Mode-3: fitted slope and offset: ',DFIT3 
      WRITE(53,'(A,A,D15.4,A)')' ', 'M-C TAUD3  ', TAUDMC3(IE0),' SEC'


      WRITE(63,'(I10,10E16.5)')IE0,NU,TAUDMC(IE0) ,RATEMC,
     +                                 TAUDMC1(IE0),RATEMC1, 
     +                                 TAUDMC2(IE0),RATEMC2, 
     +                                 TAUDMC3(IE0),RATEMC3,E0 
      ENDIF

C======Get major and minor axes of 1-sigma transverse ellipse from COVTRACK. 
      IF(IBEQUIL.EQ.0)THEN
      TILT=DATAN2(2.D0*COVTRACK(1,3,NTURN2),COVTRACK(1,1,NTURN2)
     +                                    -COVTRACK(3,3,NTURN2)) ! 2 theta
      CTILT=DCOS(TILT)                                           !cos (2 theta)
      TILT=TILT*90.D0/PI
      XY3=COVTRACK(1,1,NTURN2)*COVTRACK(3,3,NTURN2)
     +                                      -COVTRACK(1,3,NTURN2)**2
      IF(XY3.LT.0..AND.XY3.GT.-1.D-8)XY3=0.
      AXISMAJ   = DSQRT(XY3/ 
     +            ( 0.5D0*((COVTRACK(3,3,NTURN2)-COVTRACK(1,1,NTURN2))
     +             /CTILT + COVTRACK(3,3,NTURN2)+COVTRACK(1,1,NTURN2))))
      AXISMIN   = DSQRT(XY3/
     +            (-0.5D0*((COVTRACK(3,3,NTURN2)-COVTRACK(1,1,NTURN2))
     +             /CTILT - COVTRACK(3,3,NTURN2)-COVTRACK(1,1,NTURN2))))
C======Now transform the sigmas 
      DLTA  = TILT*PI/180.D0 
C======Plot the ellipse at the starting point. Take 1000 points. Also 
C      principle axes.
      DO 70 IPH = 1,1000
      PH    = 0.001D0*2.D0*PI*(IPH-1)
      HR    = AXISMAJ*DCOS(PH)    
      VR    = AXISMIN*DSIN(PH)    
      HOR   = HR*DCOS(DLTA) - VR*DSIN(DLTA)
      VER   = HR*DSIN(DLTA) + VR*DCOS(DLTA) 
      AXMAJ =    (-1.D0 + 0.002D0*(IPH-1))*DTAN(DLTA)
      AXMIN =    (-1.D0 + 0.002D0*(IPH-1))*DTAN(DLTA+PI/2.D0)
      WRITE(58,105)IPH,HOR,VER,-1.D0+0.002D0*(IPH-1),AXMAJ,AXMIN
   70 CONTINUE  
      ENDIF

C
C
C----------------------------------------------------------------------------
C======Get the approximate damping times: assumes no coupling, ignores 1st two 
C      damping times
      IF(NTURN2.GT.3)THEN
      DEC1 = COVTRACK(1,1,NTURN2)/COVTRACK(1,1,3)
      DEC1 = -0.5D0*DLOG(DEC1)/(NDAMPTURNS*(NTURN2 - 2))
      DEC1 = CIR/(3.D8*DEC1)*1.D3
      DPER1 = DEC1*NDAMP2/TDAMP(1)
      DEC2 = COVTRACK(3,3,NTURN2)/COVTRACK(3,3,3)
      DEC2 = -0.5D0*DLOG(DEC2)/(NDAMPTURNS*(NTURN2 - 2))
      DEC2 = CIR/(3.D8*DEC2)*1.D3
      DPER2 = DEC2*NDAMP2/TDAMP(3)
      DEC3 = COVTRACK(5,5,NTURN2)/COVTRACK(5,5,3)
      DEC3 = -0.5D0*DLOG(DEC3)/(NDAMPTURNS*(NTURN2 - 2))
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


      RETURN

 9999 WRITE(53,92)
   92 FORMAT(' ERROR IN EIGEN VALUE ROUTINE')
      STOP
      END

