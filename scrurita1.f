      SUBROUTINE SCRURITA1(IE0,E0,U0,CIR,COVMAT,CRAD,TDAMP,SCALEMAT)
C
C
C
C
C   ROUTINE TO HANDLE SPIN-ORBIT TRACKING AND ESTIMATE THE RATE OF DEPOLARISATION.
C   ------------------------------------------------------------------------------
C   
C   First code written Tue 7/10/2003 17.00. Everything running a week later after
C   3 days of work in total. 
C   Almost equilibrium beam with the Boege kick size on 8/10/2003 at 22.00
C
C
C   For purely orbital tracking, 100 particles for 1000 turns of eRHIC (7500/2 elements
C   after the drifts have been made much faster) take 100 secs. 
C   With 400 sections and big photons it takes 25 secs. The big photons take the same  
C   time as the matrix multiplication.
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
      PARAMETER (NPART  = 10000)  ! Maximum allowed particles. Use NPART1 of them.
      PARAMETER (LIMSECT = 1000)  ! Maximum allowed sections.
      DIMENSION ORBVEC(6,NPART),ORBVET(NPART,6)
      DIMENSION ONE(NPART),ORBAVE(6),ORBSIG(6)
      DIMENSION DISTMEAN(6),COVVEC(28)
      DIMENSION COVMAT(6,6),SCALEMAT(6,6)
      DIMENSION COVSIM(6,6),COVSIMD(6,6),COVSIMOUT(6,6),
     +           COVSIMOUTD(6,6,100)
      DIMENSION AA(6,6),TMAP(6,6,3000),AAA(6,6),SECTMAP(6,6,LIMSECT)
      DIMENSION BB(6,6),BBB(6,6),AAAS(6,6)
      DIMENSION BIGPHOT(NPART),TDAMP(6)
      DIMENSION KGBINX(-500:500),KGBINY(-500:500)
      DIMENSION VC6A(6),VC6B(6),VC6C(6),DAMP(6),DTUN(6)
      DIMENSION TREV6(6,6),WR(6),WI(6)
      DIMENSION RR6(6),RI6(6),VR6(6,6),VI6(6,6)
C     DIMENSION INTGE6(6)
      PARAMETER( LWORK=64*6 )   ! now needed for  F02EBF: MV/09.10.2007
      DIMENSION WORK(LWORK)     ! now needed for  F02EBF: MV/09.10.2007
      DIMENSION WW6(6,6),ZZ(6,6),TM6B(6,6)

C
C
      PI=3.1415926535897932D0
      PI2=2.D0*PI

      WRITE(53,103)
      WRITE(53,103)
  103 FORMAT(/,'  ')


      WRITE(53,929)
  929 FORMAT('1','Entering subroutine SCRURITA1 for M-C tracking to esta
     +blish the correct excitation strength...')

      IF(NPART1.GT.NPART)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')' Attempt to work with more than NPART particles'
      STOP
      ENDIF
C---------------------------------------------------------------------------------------
C======Set up a 6-D Gaussian distribution according to the covariance matrix
C      from emitnc.f. Use NPART1 particles.  

      WRITE(53,98)
      WRITE(53,95)COVMAT
   98 FORMAT(//,' The theoretical beam covariance matrix <Xi*Xj>(mm*mm)'
     + ,' at the 1-st beam line element:')
   95 FORMAT(T6,6F13.5)
      

C      CALL G05CBF(15)

      DISTMEAN = 0.D0
      EPS = 1.D-12
      IFAIL1 = 0
C      CALL G05EAF(DISTMEAN,6,COVMAT,6,EPS,COVVEC,28,IFAIL1)
      IF(IFAIL1.NE.0)THEN
      WRITE(53,'(A,A,I10,A)')' ', 'IFAIL1 = ',IFAIL1,' STOP'  
      ENDIF

      WRITE(53,103)
      WRITE(53,103)


      ORBVEC = 0.D0

      DO 1  I = 1, NPART1

      IFAIL2 = 0 
C      CALL G05EZF(ORBVEC(:,I),6,COVVEC,28,IFAIL2)
      IF(IFAIL2.NE.0)THEN
      WRITE(53,'(A,A,I10,A)')' ', 'IFAIL2 = ',IFAIL2,' STOP'  
      ENDIF

      IF(I.LE.5)WRITE(53,'(A,A,6F10.6)')' ',
     +                        'Trial orbit vector ',ORBVEC(:,I)
    1 CONTINUE

C======Plot histograms of projections: 
C      XYBIN  = 0.01D0
C      DO 101 IO = 1,NPART1
C      NGAUSX = ORBVEC(1,IO)*1.D3/XYBIN
C      KGBINX(NGAUSX) = KGBINX(NGAUSX) + 1
C      NGAUSY = ORBVEC(3,IO)*1.D3/XYBIN
C      KGBINY(NGAUSY) = KGBINY(NGAUSY) + 1
C  101 CONTINUE 
C      DO 102 JO = -500,500
C  102 WRITE(59,'(1X,I10,1X,I10,1X,I10)')JO,KGBINX(JO),KGBINY(JO)

C-------------------------------------------------------------------------------------
C======Get the covariance matrix for this sample with clever use of MATMUL
C      Restrict to the non-zero part for cleanliness.
      ORBVET = TRANSPOSE(ORBVEC)
      COVSIM = MATMUL(ORBVEC(:,1:NPART1),ORBVET(1:NPART1,:))/NPART1
      
      WRITE(53,89)
      WRITE(53,95)COVSIM
   89 FORMAT(//,' The generated beam covariance matrix <Xi*Xj>(mm*mm)',
     +                             'at the 1-st beam line element:')


C------------------------------------------------------------------------------------
C======Transport the NPART1 particles NPART/NPART1 times with the symplectic
C      6x6 matrix to fill up the phase space some more so that a better 
C      estimate of the covariance matrix can be made. 
C      This is only safe if the tunes are not too rational.

      NBLOCK = NPART/NPART1
C======First get the 1-turn symplectic matrix again.
      IDAMPFLG = 0
      CALL UNIT(6,BBB)
      DO 8  II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      BB  = TMAT(1:6,1:6,ITY)
      CALL MX66DAMP(BB,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),CRAD,
     +                                        NAME(ITY),IDAMPFLG)
      BBB = MATMUL(BB,BBB)
    8 CONTINUE

      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')
     +          ' 1-turn symplectic orbit matrix, generated afresh.'
      DO 36 I=1,6
   36 WRITE(53,914)(BBB(I,J),J=1,6)
       

      DO 9 IB = 1,NBLOCK
      ORBVEC(:,IB*NPART1+1:(IB+1)*NPART1) 
     + = MATMUL(BBB,ORBVEC(:,(IB-1)*NPART1+1:IB*NPART1))
    9 CONTINUE    

C======The orbit has not yet been rescaled.
C======Get the covariance matrix with clever use of MATMUL. Now use the full ORBVEC
      ORBVET  = TRANSPOSE(ORBVEC)
      COVSIMD = MATMUL(ORBVEC,ORBVET)/NPART
      
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,99)
      WRITE(53,95)COVSIMD
   99 FORMAT(//,' The generated beam covariance matrix <Xi*Xj>(mm*mm)',
     + 'at the 1-st beam line element for a densified distribution:')

C      Check the beam average.
      ONE = 1.D0
      ORBAVE = MATMUL(ORBVEC(:,1:NPART),ONE(1:NPART))/NPART
      
      WRITE(53,88)
      WRITE(53,95)ORBAVE
   88 FORMAT(//,' The generated beam average <Xi>(mm)',
     + ' at the 1-st beam line element for a densified distribution:')
 

      
C      Put the orbit vector into metres, rad  etc
C      ORBVEC = ORBVEC*1.D-3

C-------------------------------------------------------------------------------------
C======Make 1 pass around the ring to update the element matrices and to generate
C      matrices for sections between the dipole centres using the 
C      fact that the C.O. is fixed. So don't have to recalc all fixed stuff on each turn
C      and can track using the sections.
      NSECT = 0
      IDAMPFLG = 1                               !Switch on/off radn.
      CALL UNIT(6,AAA)
      DO 3  II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      AA  = TMAT(1:6,1:6,ITY)
      CALL MX66DAMP(AA,ITY,IID,II,XX(ITY),X2(ITY),YY(ITY),CRAD,
     +                                        NAME(ITY),IDAMPFLG)

      TMAP(:,:,ITY) = AA

      AAA = MATMUL(AA,AAA)
      IF((IID.EQ.7.AND.NAME(ITY)(1:2).EQ.'VD').OR.II.EQ.NELEM)THEN  
      NSECT = NSECT + 1

      IF(NSECT.GT.LIMSECT)THEN
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A)')' STOP:the number of sections exceeds LIMSECT'
      STOP
      ENDIF

      SECTMAP(:,:,NSECT) = AAA
      CALL UNIT(6,AAA)      
      ENDIF

    3 CONTINUE

       
C======Check the 1-turn matrix using the sections.
      CALL UNIT(6,AAA)
      DO 7 II = 1,NSECT
      AAA = MATMUL(SECTMAP(:,:,II),AAA)
    7 CONTINUE


      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,A,I6)')' ','NSECT', NSECT
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,A,I5)')' 1-turn damped orbit matrix using sections.',
     +               '   NDAMP1= ',NDAMP1
      DO 35 I=1,6
   35 WRITE(53,914)(AAA(I,J),J=1,6)
  914 FORMAT(6F12.5)

C      IF(1.EQ.1)RETURN

C======Get the tunes to gauge the tune shift due to the damping.
      AAAS = AAA
      IFAIL=0
C      CALL F02AGF(AAAS,6,6,RR6,RI6,VR6,6,VI6,6,INTGE6,IFAIL)
      CALL F02EBF('V',6,AAAS,6,RR6,RI6,VR6,6,VI6,6,WORK,LWORK,IFAIL) ! <<replmnt
      IF(IFAIL.NE.0)GO TO 9999
      DO 2822 I=1,3
      IT=2*I-1
      ZZ(:,IT  )=VR6(:,IT)
      ZZ(:,IT+1)=VI6(:,IT)
 2822 CONTINUE
      WR = RR6
      WI = RI6
      TM6B = ZZ
      CALL RENORM(TM6B,6,WR,WI)
      WRITE(53,39)
   39 FORMAT(//,' DAMPED EIGEN-TUNES AND EIGENVECTORS:',/,T12,
     +    'COS(2*PI*NU)',T27,'SIN(2*PI*NU)',T42,'EIGENVECTORS',T112,
     +    'TUNES(OSC.PART)')
      DO 38 I=1,6
      DTUN(I)=DATAN2(WI(I),WR(I))/(2.*PI)
   38 WRITE(53,999)WR(I),WI(I),(TM6B(J,I),J=1,6),DTUN(I)
 999  FORMAT(T12,F12.5,T27,F12.5,T41,'(',6F11.5,')',F15.8)
C=====EXTRACT DAMPING CONSTANTS FROM COMPLEX EIGENVALUES.
      DO 23 ISX=1,6
      VC6A(ISX)=DSQRT(WR(ISX)**2+WI(ISX)**2)
      XY7=-DLOG(VC6A(ISX))
      VC6C(ISX)=XY7
      VC6B(ISX)=XY7*2.*E0*1000./U0
   23 DAMP(ISX)=CIR/(3.D8*XY7)*1000.
      ERAT=2.D0*U0/E0/1000 *NDAMP1
      WRITE(53,366) WR,WI,VC6A,VC6C,ERAT,VC6B,DAMP
  366 FORMAT(//,' EIGENVALUES WITH RADIATION DAMPING:',/,
     +    T7,'REAL',           T18,6(1X,F11.6),/,
     +    T7,'IMAG',           T18,6(1X,F11.6),/,
     +    T7,'ABS ',           T18,6(1X,F11.6),/,
     +    T3,'DAMPING CONST.', T18,6(1X,F11.6),'  2U0*NDAMP1/E= ',F9.5/,
     +    T3,'PARTITION NO.',  T18,6(1X,F11.6),/,
     +    T3,'DAMP.TIME(MSEC)',T18,6(1X,F11.6))
C
C
C


C------------------------------------------------------------------------------------ 
      WRITE(53,103)
      WRITE(53,103)

C      WRITE(53,'(A,A,F20.8)')' ','CRAD',CRADIN

C======Now track to get the equilibrium beam
C======After each dipole or C-F, radiate ``big'' photons according to a centred 
C      top-hat distribution. 
C      In contrast to a Gaussian, there are no long tails so that no cut against
C      too large (particle expelling kicks) is needed.
C      For linear motion this should also be enough to get Gaussian phase space
C      by the central limit theorem. Get a vector of NPART1 samples.  
C      Use the current seed.
C      Use the Boege formula for the strength.
C      Later, scale the strengths according to the curvature and length. 
C      No radiation in the correction coils. 
C      Scale the damping up with a factor NDAMP1 to decrease the damping time  
C      and thus get answers quicker. Track for 10 compressed sync damping times
C
C      For linear tracking, no need to scale the orbit vector to metres etc. 
C      but do it anyway.

      ORBVEC = ORBVEC*1.D-3
      
      CIRCTIME =  CIR/3.0D8
      NDAMPTURNS = 1.D-3*TDAMP(5)/CIRCTIME/NDAMP1 * 1.D0
      WRITE(53,'(A,I6,A,I6)')' Turns per sync. damping time ', 
     +                       NDAMPTURNS,'  NDAMP1=  ', NDAMP1
      IDAMPFLG = 1                              !Switch on/off radn.
      IMAG = IBEND + ICOMF
      TH = 0.5D0
      VARTH = TH*TH/3.D0
      PHOTSCALE = CIRCTIME*SCALEMAT(6,6)/(1.D-3*TDAMP(5)*VARTH*NSECT)
      PHOTSCALE = 2.D0*1.D-3*DSQRT(PHOTSCALE*NDAMP1)*IDAMPFLG

      WRITE(53,'(A,A,I6,E15.7)')' ','NSECT,PHOTSCALE', NSECT,PHOTSCALE

      IDAMPCNT = 0
      DO 4 IT = 1,NDAMPTURNS*NTURN1
      DO 5 II = 1,NSECT
      ORBVEC(:,1:NPART1) = MATMUL(SECTMAP(:,:,II),ORBVEC(:,1:NPART1))
C======Radiate.
C      CALL G05FAF(-TH,TH,NPART1,BIGPHOT)
      ORBVEC(6,1:NPART1) = ORBVEC(6,1:NPART1) 
     +                                   + BIGPHOT(1:NPART1)*PHOTSCALE
    5 CONTINUE

C======At the end of each sync damping time, get the covariance matrix 
C      for the densified distribution. 
      IF(MOD(IT,NDAMPTURNS).EQ.0)THEN
C      Clean up ORBVEC for safety.
      ORBVEC(:,NPART1+1:NPART) = 0.D0
      DO 11 IB = 1,NBLOCK
      ORBVEC(:,IB*NPART1+1:(IB+1)*NPART1) 
     + = MATMUL(BBB,ORBVEC(:,(IB-1)*NPART1+1:IB*NPART1))
   11 CONTINUE    

C======The orbit has not yet been rescaled.
      IDAMPCNT =  IDAMPCNT + 1     
      ORBVET = TRANSPOSE(ORBVEC)
      COVSIMOUTD(:,:,IDAMPCNT) = MATMUL(ORBVEC,ORBVET)*1.D6/NPART
      
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,86)IDAMPCNT,IT
      WRITE(53,95)COVSIMOUTD(:,:,IDAMPCNT)
   86 FORMAT(//,' The ',I3,'th  transported beam covariance matrix',
     + ' <Xi*Xj>(mm*mm)',
     + ' at the 1-st beam line element for a densified distribution',
     + ' after ',I6,' turns')
      ENDIF

    4 CONTINUE   

C======Check the variance of the last big photons. 
      VARPHOT = 0
      DO 6 IPH = 1,NPART1
    6 VARPHOT = VARPHOT + BIGPHOT(IPH)**2
      VARPHOT = VARPHOT/NPART1
      WRITE(53,'(A,A,F10.6)')' ','Variance of the top-hat samples ',
     +                           VARPHOT


C-------------------------------------------------------------------------------------
C      Check the covariance matrix etc after tracking: rescale to millmetres etc
      ORBVET = TRANSPOSE(ORBVEC)
      COVSIMOUT = 
     +     MATMUL(ORBVEC(:,1:NPART1),ORBVET(1:NPART1,:))*1.D6/NPART1

      ONE = 1.D0
      ORBAVE = MATMUL(ORBVEC(:,1:NPART1),ONE(1:NPART1))*1.D3/NPART1
      
      WRITE(53,96)
      WRITE(53,95)ORBAVE
   96 FORMAT(//,' The transported beam average <Xi>(mm)',
     +           ' at the 1-st beam line element:')
 
      DO 94 IJ=1,6
   94 ORBSIG(IJ)=DSQRT(COVSIMOUT(IJ,IJ))
      WRITE(53,93)
      WRITE(53,95)ORBSIG
   93 FORMAT(//,' The transported rms beam sizes(mm)',
     +           ' at the 1-st beam line element:')

      WRITE(53,97)
      WRITE(53,95)COVSIMOUT
   97 FORMAT(//,' The transported beam covariance matrix <Xi*Xj>(mm*mm)'
     +           ,' at the 1-st beam line element:')

C======Plot histograms of horizontal and vertical projections: 
      IF(IDEPLIN.EQ.0)THEN
      XYBIN  = 0.01D0
      DO 101 IO = 1,NPART1
      NGAUSX = ORBVEC(1,IO)*1.D3/XYBIN
      KGBINX(NGAUSX) = KGBINX(NGAUSX) + 1
      NGAUSY = ORBVEC(3,IO)*1.D3/XYBIN
      KGBINY(NGAUSY) = KGBINY(NGAUSY) + 1
  101 CONTINUE 
      DO 102 JO = -500,500
  102 WRITE(59,'(1X,I10,1X,I10,1X,I10)')JO,KGBINX(JO),KGBINY(JO)
      ENDIF

C----------------------------------------------------------------------------------------
C      Now densify the phase space. Clean up ORBVEC for safety.
      ORBVEC(:,NPART1+1:NPART) = 0.D0
      DO 10 IB = 1,NBLOCK
      ORBVEC(:,IB*NPART1+1:(IB+1)*NPART1) 
     + = MATMUL(BBB,ORBVEC(:,(IB-1)*NPART1+1:IB*NPART1))
  10  CONTINUE    

C======The orbit has not yet been rescaled.
C======Get the covariance matrix with clever use of MATMUL. Now use the full ORBVEC
      ORBVET = TRANSPOSE(ORBVEC)
      COVSIMOUTD(:,:,IDAMPCNT) = MATMUL(ORBVEC,ORBVET)*1.D6/NPART
      
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,87)
      WRITE(53,95)COVSIMOUTD(:,:,IDAMPCNT)
   87 FORMAT(//,' The transported beam covariance matrix <Xi*Xj>(mm*mm)'
     + ,' at the 1-st beam line element for a densified distribution:')
C-------------------------------------------------------------------------
C======Get major and minor axes of 1-sigma transverse ellipse from COVSIMOUT. 
      TILT=DATAN2(2.D0*COVSIMOUT(1,3),COVSIMOUT(1,1)-COVSIMOUT(3,3)) ! 2 theta
      CTILT=DCOS(TILT)                                           !cos (2 theta)
      TILT=TILT*90.D0/PI
      XY3=COVSIMOUT(1,1)*COVSIMOUT(3,3)-COVSIMOUT(1,3)**2
      IF(XY3.LT.0..AND.XY3.GT.-1.D-8)XY3=0.
      AXISMAJ   = DSQRT(XY3/ 
     +            ( 0.5D0*((COVSIMOUT(3,3)-COVSIMOUT(1,1))/CTILT + 
     +                COVSIMOUT(3,3)+COVSIMOUT(1,1))))
      AXISMIN   = DSQRT(XY3/
     +            (-0.5D0*((COVSIMOUT(3,3)-COVSIMOUT(1,1))/CTILT - 
     +                COVSIMOUT(3,3)-COVSIMOUT(1,1))))
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
  105 FORMAT(' ',I10, 6E16.6) 
   70 CONTINUE  
C
C
C
C----------------------------------------------------------------------------------------
C======Get the approximate damping times: assumes no coupling.
      DEC1 = COVSIMOUTD(1,1,IDAMPCNT)/COVSIMOUTD(1,1,3)
      DEC1 = -0.5D0*DLOG(DEC1)/(NDAMPTURNS*7)
      DEC1 = CIR/(3.D8*DEC1)*1.D3
      DPER1 = DEC1*NDAMP1/TDAMP(1)
      DEC2 = COVSIMOUTD(3,3,IDAMPCNT)/COVSIMOUTD(3,3,3)
      DEC2 = -0.5D0*DLOG(DEC2)/(NDAMPTURNS*7)
      DEC2 = CIR/(3.D8*DEC2)*1.D3
      DPER2 = DEC2*NDAMP1/TDAMP(3)
      DEC3 = COVSIMOUTD(5,5,IDAMPCNT)/COVSIMOUTD(5,5,3)
      DEC3 = -0.5D0*DLOG(DEC3)/(NDAMPTURNS*7)
      DEC3 = CIR/(3.D8*DEC3)*1.D3
      DPER3 = DEC3*NDAMP1/TDAMP(5)
      WRITE(53,103)
      WRITE(53,103)
      WRITE(53,'(A,3(1X,F11.1),A,3(1X,F11.1))')
     +                           ' Damping times (msec) ',
     +                            DEC1,DEC2,DEC3,  
     +                           '  Damping periods with NDAMP1 ',
     +                            DPER1,DPER2,DPER3


      RETURN

 9999 WRITE(53,92)
   92 FORMAT(' ERROR IN EIGEN VALUE ROUTINE')
      STOP
      END
