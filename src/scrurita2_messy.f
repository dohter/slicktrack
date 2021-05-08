      SUBROUTINE SCRURITA2(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP)
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
C
      PARAMETER (NPART   = 10000) ! Maximum allowed particles. Use NPART2 of them.
      PARAMETER (LIMSECT = 1000)  ! Maximum allowed sections.
      REAL*8   NU
      DIMENSION ZR3(3,3),ZI3(3,3),A(3,3),B(3,3)
      DIMENSION ROT(3,3),TM3A(3,3)
      DIMENSION TM3B(3,3),WR3(3),WI3(3),ZW(3,3),ZWSAVE(3,3),TM3C(3,3)
      DIMENSION TRIN3(3,3),RR3(3),RI3(3),VR3(3,3),VI3(3,3),INTGE3(3)
      DIMENSION TW8A(8,8),TW2A(2,2)
      DIMENSION DISTMEAN(6),COVVEC(28),COVMAT(6,6)
      DIMENSION COVSIM(8,8),COVTRACK(8,8,1500)
      DIMENSION SORBVEC(8,NPART),SORBVET(NPART,8)
      DIMENSION ONE(NPART),SORBAVE(8),SORBSIG(8)
      DIMENSION BB(8,8),AAA(8,8),SECTMAP(8,8,LIMSECT)
      DIMENSION TREV8(8,8),TREV8D(8,8)
      DIMENSION BIGPHOT(NPART),TDAMP(6)
      DIMENSION SDIFFUSION(1500),   DFIT(2),   DVEC(2), DIFFMAT(2,2)
      DIMENSION SDIFFUSIONA(1500),  DFITA(2),  DVECA(2)
      DIMENSION SDIFFUSIONB(1500),  DFITB(2),  DVECB(2)
      DIMENSION SDIFFUSIONAB(1500), DFITAB(2), DVECAB(2)
      DIMENSION ODIFFUSION(6,1500)

      DIMENSION KGBINX(-500:500),KGBINY(-500:500)
      DIMENSION SORBDUM(8)
C
C
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

      IF(NPART2.GT.NPART)THEN
      WRITE(53,103)
      WRITE(53,'(A)')' STOP: attempt to work with more than NPART',
     +                 '  particles'
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
C      IF(DABS(WI3(I)).LT.1.D-9.AND.WR3(I).GT.0.9999D0)ICOUNT=ICOUNT+1
  255 IF(DABS(WI3(I)).LT.1.D-9.AND.WR3(I).GT.0.9999D0)NN=I
      IF(NN .NE. 4)GO TO 241
      WRITE(53,242)
  242 FORMAT(' ','No real unit spin eigenvalue found----so STOP')
      STOP
C
  241 CONTINUE
CC      IF(ICOUNT.EQ.1)GO TO 243
C     NN=1
CC  243 CONTINUE
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
C      Exchange m and l.
C     S1 = ZW(1,2)
C     S2 = ZW(2,2)
C     S3 = ZW(3,2)
C     ZW(1,2) = ZW(1,3)
C     ZW(2,2) = ZW(2,3)
C     ZW(3,2) = ZW(3,3)
C     ZW(1,3) = S1
C     ZW(2,3) = S2
C     ZW(3,3) = S3
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
C---------------------------------------------------------------------------------------
C======Set up a 6-D Gaussian distribution according to the covariance matrix
C      from emitnc.f. Use NPART2 particles.  

      WRITE(53,98)
C      write(*,98)

      WRITE(53,951)COVMAT
   98 FORMAT(//,' The theoretical beam covariance matrix <Xi*Xj>'
     + ,'(mm*mm,mm*mrad...) at the 1-st beam line element:')
  951 FORMAT(T6,6F13.5)

      CALL G05CBF(15)                  !No freedom to chosse the starting seed yet.

      DISTMEAN = 0.D0
      EPS = 1.D-12
      IFAIL1 = 0
      CALL G05EAF(DISTMEAN,6,COVMAT,6,EPS,COVVEC,28,IFAIL1)
      IF(IFAIL1.NE.0)THEN
      WRITE(53,'(A,A,I10,A)')' ', 'IFAIL1 = ',IFAIL1,' STOP'  
      ENDIF

      WRITE(53,103)
      WRITE(53,103)



      SORBVEC = 0.D0

      DO 1  I = 1, NPART2
      IFAIL2 = 0 
C      CALL G05EZF(SORBVEC(1:6,I),6,COVVEC,28,IFAIL2)
      IF(IFAIL2.NE.0)THEN
      WRITE(53,'(A,A,I10,A)')' ', 'IFAIL2 = ',IFAIL2,' STOP'  
      ENDIF
      IF(I.LE.5)WRITE(53,'(A,A,I5,8F10.6)')' ',
     +                  'Starting spin-orbit vector ',I,SORBVEC(1:8,I)
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

C-------------------------------------------------------------------------------------
C======Get the covariance matrix for this sample with clever use of MATMUL
C      Restrict to the non-zero part for cleanliness.
      SORBVET(1:NPART2,:) = TRANSPOSE(SORBVEC(:,1:NPART2))
      COVSIM = MATMUL(SORBVEC(:,1:NPART2),SORBVET(1:NPART2,:))/NPART2
C      SORBVET = TRANSPOSE(SORBVEC)
C      COVSIM  = MATMUL(SORBVEC,SORBVET)/NPART2      
      WRITE(53,89)
C      write(*,89)
   89 FORMAT(//,' The generated spin-orbit covariance matrix',
     +     ' <Xi*Xj>(mm*mm,mm*mrad...) at the 1-st beam line element:')
      WRITE(53,955)COVSIM
  955 FORMAT(T6,8F13.5)

C------------------------------------------------------------------------------------
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

C----------------------------------------------------------------------
C======Get the 1-turn spin-(damped) orbit matrix.
      IDAMPFLG = 1
      S=0.
      CALL UNIT(8,TREV8D)
      TM3C = ZWSAVE
      DO 88  II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
      IF(IID.EQ.1.AND.XY.LT.0.00011D0)GO TO 88 !Skip zero length drifts 
      IF(IID.EQ.1)THEN
      TREV8D(1,:) = TREV8D(1,:) + XY * TREV8D(2,:) 
      TREV8D(3,:) = TREV8D(3,:) + XY * TREV8D(4,:) 
      GO TO 88
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
   88 CONTINUE

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
C======Make 1 pass around the ring to generate damped spin-orbit matrices for sections
C      between the dipole centres using the fact that the C.O. is fixed.
C      So don't have to recalc all fixed stuff on each turn and can track using the sections.
      NSECT = 0
      IDAMPFLG = 1                               !Switch on/off radn.
      CALL UNIT(8,AAA)
      TM3C = ZWSAVE
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

      SECTMAP(:,:,NSECT) = AAA
      CALL UNIT(8,AAA)      
      ENDIF

    7 CONTINUE

       
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
      CALL UNIT(8,TW8A)
      TW8A(7,7)= WR3(2)
      TW8A(7,8)= WI3(2)
      TW8A(8,7)=-WI3(2)
      TW8A(8,8)= WR3(2)
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


C------------------------------------------------------------------------------------ 
      WRITE(53,103)
      WRITE(53,103)

C      WRITE(53,'(A,A,F20.8)')' ','CRAD',CRADIN

C======Now track to get the equilibrium beam and the non-equilibrium spin components.
C======After each dipole or C-F, radiate ``big'' photons according to a centred 
C      top-hat distribution. 
C      In contrast to a Gaussian, there are no long tails so that no cut against
C      too large (particle expelling kicks) is needed.
C      For linear motion this should also be enough to get Gaussian phase space
C      by the central limit theorem. Get a vector of NPART2 samples.  
C      Use the current seed.
C      Use the Boege formula for the strength. The rms ``big'' photon energy is
C      very much smaller than th eenergy spread of the beam.
C      All ``big'' photons have the same energy spread.
C      Later, scale the strengths according to the curvature and length of the dipoles. 
C      This could be important for dipoles where n_0 is horizontal and there is no 
C      spin match.
C      No radiation in the correction coils. 
C      Scale the damping up with a factor NDAMP2 to decrease the damping time  
C      and thus get answers quicker. Track for 10 compressed sync damping times
C
C      For linear tracking, no need to scale the orbit vector to metres etc. 
C      but do it anyway as a reminder for the future when nonlinear stuff is used.

      SORBVEC = SORBVEC*1.D-3
      
      CIRCTIME =  CIR/3.0D8
      DAMPTURNS = 1.D-3*TDAMP(5)/CIRCTIME/NDAMP2
      WRITE(53,'(A,F8.1,A,I3)')' Turns per sync. damping time', 
     +                       DAMPTURNS,'  NDAMP2=', NDAMP2
      NDAMPTURNS = DAMPTURNS   
      IDAMPFLG = 1                              !Switch on/off radn.
      TH = 0.5D0
      VARTH = TH*TH/3.D0
      PHOTSCALE = CIRCTIME*COVMAT(6,6)/(1.D-3*TDAMP(5)*VARTH*NSECT)
      PHOTSCALE = 2.D0*1.D-3*DSQRT(PHOTSCALE*NDAMP2)*IDAMPFLG
      TH = TH * PHOTSCALE
      WRITE(53,'(A,A,I6,E15.7)')' ','NSECT,PHOTSCALE', NSECT,PHOTSCALE

      IDNT = 0
      DO 9 INT = 1,NTURN2               !Loop over damping periods.
      DO 4 IDT = 1,NDAMPTURNS           !Loop over turns in a damping period. 
C      CALL UNIT(8,AAA)
      DO 5 II  = 1,NSECT                !Loop over sections in a turn.
      SORBVEC(:,1:NPART2) = MATMUL(SECTMAP(:,:,II),SORBVEC(:,1:NPART2))
C     SORBVEC = MATMUL(SECTMAP(:,:,II),SORBVEC)

C      write(*,*)II,IDT,INT

C======Radiate.
      CALL G05FAF(-TH,TH,NPART2,BIGPHOT)
      SORBVEC(6,1:NPART2) = SORBVEC(6,1:NPART2) + BIGPHOT(1:NPART2)
C      SORBVEC(6,:) = SORBVEC(6,:) + BIGPHOT
C      AAA = MATMUL(SECTMAP(:,:,II),AAA)
    5 CONTINUE

C======WIND BACK THE SPIN BASIS AT THE END OF A TURN=================================
      CALL UNIT(8,TW8A)
      CALL UNIT(2,TW2A)
      TW8A(7,7)= WR3(2)
      TW8A(7,8)= WI3(2)
      TW8A(8,7)=-WI3(2)
      TW8A(8,8)= WR3(2)

C      TW8A(7,7)= WR3(2)
C      TW8A(7,8)=-WI3(2)
C      TW8A(8,7)= WI3(2)
C      TW8A(8,8)= WR3(2)

      TW2A(1,1)= WR3(2)
      TW2A(1,2)= WI3(2)
      TW2A(2,1)=-WI3(2)
      TW2A(2,2)= WR3(2)
C======WIND FORWARD THE SPIN BASIS AT THE END OF A TURN--as a test ===================
C      TM2A(1,1)= WR3(2)
C      TM2A(1,2)=-WI3(2)
C      TM2A(2,1)= WI3(2)
C      TM2A(2,2)= WR3(2)

      SORBVEC(:,1:NPART2)  = MATMUL(TW8A,SORBVEC(:,1:NPART2))
CCCCCC      SORBVEC  = MATMUL(TW8A,SORBVEC)
C      AAA      = MATMUL(TW8A,AAA)
C      SORBVEC(7:8,1:NPART2)  = MATMUL(TW2A,SORBVEC(7:8,1:NPART2))
C      IF(IT.LT.4)THEN
C      WRITE(53,103)
C      WRITE(53,'(A,I5,A)')' ',IT,' dampings OTM'
C      DO 338 I=1,8
C  338 WRITE(53,914)(AAA(I,J),J=1,8)
C      ENDIF


CC      IDNT = IDNT + 1
CCC      SORBVET = TRANSPOSE(SORBVEC)
CCC      COVTRACK(:,:,IDNT) = MATMUL(SORBVEC,SORBVET)*1.D6/NPART2      
CCC      SORBVET(1:NPART2,:) = TRANSPOSE(SORBVEC(:,1:NPART2))
CCC      COVTRACK(:,:,IDNT) 
CCC     +   = MATMUL(SORBVEC(:,1:NPART2),SORBVET(1:NPART2,:))*1.D6/NPART2
CCC      WRITE(53,103)
CCC      WRITE(53,103)
CCC      WRITE(53,86)IDNT,IDNT*IDNT
CCC      WRITE(53,955)COVTRACK(:,:,IDNT)
CCC   86 FORMAT(//,' The ',I3,'th  transported s-o covariance matrix',
CCC     +          ' <Xi*Xj>(mm*mm,mm*mrad...)',
CCC     +          ' at the 1-st beam line element after ',I6,' turns')
CCC      WRITE(53,103)
CCC      WRITE(53,103)
CCC      SDIFFUSION(IDNT)    = COVTRACK(7,7,IDNT)+COVTRACK(8,8,IDNT)
CCC      SDIFFUSIONA(IDNT)   = COVTRACK(7,7,IDNT)
CCC      SDIFFUSIONB(IDNT)   = COVTRACK(8,8,IDNT)
CCC      SDIFFUSIONAB(IDNT)  = COVTRACK(7,8,IDNT)
CCC      DO 87 IK = 1,6
CCC   87 ODIFFUSION(IK,IDNT) = COVTRACK(IK,IK,IDNT)







    4 CONTINUE   
C======At the end of each sync damping time, get the covariance matrix 
C======The orbit has not yet been rescaled.
C      SORBVET(1:NPART2,:) = TRANSPOSE(SORBVEC(:,1:NPART2))
CCC      SORBVET = TRANSPOSE(SORBVEC)
CCC      COVTRACK(:,:,INT) = MATMUL(SORBVEC,SORBVET)*1.D6/NPART2      
      SORBVET(1:NPART2,:) = TRANSPOSE(SORBVEC(:,1:NPART2))
      COVTRACK(:,:,INT) 
     +     = MATMUL(SORBVEC(:,1:NPART2),SORBVET(1:NPART2,:))*1.D6/NPART2
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,86)INT,INT*IDT
      WRITE(53,955)COVTRACK(:,:,INT)
   86 FORMAT(//,' The ',I3,'th  transported s-o covariance matrix',
     +          ' <Xi*Xj>(mm*mm,mm*mrad...)',
     +          ' at the 1-st beam line element after ',I6,' turns')
      WRITE(53,103)
      WRITE(53,103)
      SDIFFUSION(INT)    = COVTRACK(7,7,INT)+COVTRACK(8,8,INT)
      SDIFFUSIONA(INT)   = COVTRACK(7,7,INT)
      SDIFFUSIONB(INT)   = COVTRACK(8,8,INT)
      SDIFFUSIONAB(INT)  = COVTRACK(7,8,INT)
      DO 87 IK = 1,6
   87 ODIFFUSION(IK,INT) = COVTRACK(IK,IK,INT)


      WRITE(53,'(A,A,F13.5,A,F13.5)')' ','(7,7) + (8,8): ',
     +                             SDIFFUSION(INT),'    (7,7)/(8,8): ',
     +                         COVTRACK(7,7,INT)/COVTRACK(8,8,INT)

    9 CONTINUE

C======Check the variance of the last big photons. 
      WRITE(53,103)
      WRITE(53,103)
      VARPHOT = 0
      DO 6 IPH = 1,NPART2
    6 VARPHOT = VARPHOT + BIGPHOT(IPH)**2
      VARPHOT = VARPHOT/NPART2/PHOTSCALE**2
      WRITE(53,'(A,A,F10.6)')' ','Variance of the top-hat samples ',
     +                           VARPHOT


C-------------------------------------------------------------------------------------
C      Print some statistical stuff. Rescale to mm. etc.

      ONE = 1.D0
      SORBAVE = MATMUL(SORBVEC(:,1:NPART2),ONE(1:NPART2))*1.D3/NPART2
C      SORBAVE = MATMUL(SORBVEC,ONE)*1.D3/NPART2      
      WRITE(53,96)
      WRITE(53,95)SORBAVE
   96 FORMAT(//,' The transported spin-beam average <Xi>(mm,mrad...)',
     +           ' at the 1-st beam line element:')
   95 FORMAT(T6,8F13.5)
 
      DO 94 IJ=1,8
   94 SORBSIG(IJ)=DSQRT(COVTRACK(IJ,IJ,NTURN2))
      WRITE(53,93)
      WRITE(53,95)SORBSIG
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
      IF(NTURN2.GT.3)THEN
C      IF(NDAMPTURNS.GT.3)THEN
      DIFFMAT = 0.D0  
      DVEC    = 0.D0
      DVECA   = 0.D0
      DVECB   = 0.D0

      DO 106 JF = 3, NTURN2
C      DO 106 JF = 3, NDAMPTURNS
C      DO 106 JF = 3, IDNT
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
C======Plot the diffusion data beyond any transients.

      DO 107 KF = 1,NTURN2
C      DO 107 KF = 1,NDAMPTURNS
C      DO 107 KF = 1,IDNT
 107  WRITE(61,'(I10,16F13.5)')KF,
     +                        SDIFFUSION(KF),  DFIT(1)*KF   +DFIT(2),
     +                        SDIFFUSIONA(KF), DFITA(1)*KF  +DFITA(2),
     +                        SDIFFUSIONB(KF), DFITB(1)*KF  +DFITB(2),
     +                        SDIFFUSIONAB(KF),DFITAB(1)*KF +DFITAB(2),
     +                        (ODIFFUSION(IORB,KF),IORB=1,6)
C======Convert DFIT(1) to a depolarisation time. Correct for using enhanced radiation/damping.
      WRITE(53,103)
      WRITE(53,'(A,A,F13.5)')' ','Slope ratios (7,7)/(8,8): ',
     +                                               DFITA(1)/DFITB(1)

      TAUDMC = 0.5D0*DFIT(1)*1.D-6/(DAMPTURNS*TDAMP(5)*1.D-3) * NDAMP2 
      TAUDMC = 1.D0/TAUDMC
      WRITE(53,103)
      WRITE(53,'(A,A,2D15.8,A)')' ', 'Fitted slope and offset: ',DFIT 
      WRITE(53,103)
      WRITE(53,'(A,A,D15.4,A)')' ', 'M-C TAUD  ', TAUDMC,' SEC'
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
C======Plot the ellipse at the starting point. Take 1000 points. Also principle axes.
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
C----------------------------------------------------------------------------------------
C======Get the approximate damping times: assumes no coupling, ignores 1st two damping times
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

C      DO 440 IIP = 1,1000
C         DO 441 IIR = 1,8
C            SORBDUM(IIR) = 0.D+0
C            DO 442 IIC = 1,8
C               SORBDUM(IIR) = SORBDUM(IIR)
C     +              +  SECTMAP(IIR,IIC,II)*SORBVEC(IIC,IIP)
C 442        CONTINUE
C 441     CONTINUE
C         DO 443 IIR=1,8
C            SORBVEC(IIR,IIP) = SORBDUM(IIR)
C 443     CONTINUE
C 440  CONTINUE
