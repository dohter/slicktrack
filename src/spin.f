C   07/03/95 509071949  MEMBER NAME  SPIN     (SEPT95.S) M  FORTRAN
C   22/07/81 503061853  MEMBER NAME  SPIN     (MAR95.S)  M  FORTRAN
      SUBROUTINE SPIN(IE0,E0,CIR,TAUY,BETAY0,PTUN,NCRAZY)
C
C
C
C
C   ROUTINE TO HANDLE 8X8 SPIN EIGEN ANALYSIS & COMPUTE POLARIZATION
C   ----------------------------------------------------------------
C
C
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
      INCLUDE "cspindiff.for"
      INCLUDE "cpol.for"
C
      DIMENSION ZR3(3,3),ZI3(3,3),A(3,3),B(3,3),P(3,3),PTEMP(3,3)
      DIMENSION PDIF(3,3),PTOT(3,3)
      DIMENSION ROT(3,3),TM3A(3,3)
      DIMENSION TM3B(3,3),WR3(3),WI3(3),ZW(3,3),TM3C(3,3)
      DIMENSION TRIN3(3,3),RR3(3),RI3(3),VR3(3,3),VI3(3,3)
      DIMENSION INTGE3(3),WW3(3,3)
      DIMENSION TREV8(8,8),ZZ(8,8),ZV(8,8),WR8(8),WI8(8)
      DIMENSION TREV08(8,8),RESSTR(2,6)                  ! For resonance strengths.
      DIMENSION TRIN8(8,8),RR8(8),RI8(8),VR8(8,8),VI8(8,8)
C      DIMENSION INTGE8(8),WW8(8,8)
      PARAMETER( LWORK=64*8 )   ! now needed for  F02EBF: MV/09.10.2007
      DIMENSION WORK(LWORK)     ! now needed for  F02EBF: MV/09.10.2007
      DIMENSION WW8(8,8)
      DIMENSION TM8A(8,8),TM8B(8,8)
      DIMENSION SOL(8,8)
      DIMENSION AB(6)
      DIMENSION TN(8),PTUN(6)
      REAL*8 NU
      DATA ZV/64*0.D0/
      LOGICAL IPRIN1,IPRIN2
C
C
C
C
C=====STORAGE FOR THE F02AGF ROUTINE
C     DIMENSION TRIN8(8,8),RR(8),RI(8),VR(8,8),VI(8,8),INTGER(8),WW(8,8)
C
C
C
C
C
C
C
      PI=3.1415926535897932D0
      PI2=2.D0*PI
      NU=E0/0.440652D0 * 1.0D0       ! Can scale a gamma indep. of the energy.
      SUMNZ = 0
      SUMNZA= 0
      SUMNZM= 0
      SUMNTM= 0
      NNZ=0
      NNZA=0
      NNZM=0
      NNTM=0
C
C
C
C
C   *************************************************************
C   * SPIN ROTATION MATRIX AND THE ORTHONORMAL SPIN BASE VECTORS*
C   * AT THE FIRST BEAM-LINE ELEMENT                            *
C   *************************************************************
C
      WRITE(53,929)IE0
  929 FORMAT('1','Entering Subroutine SPIN. Energy step = ',I5)
C
C
C=====CALCULATE THE SPIN REVOL. MATRIX=================================
      CALL UNIT(3,ROT)
      ANGH=0.D0
      ANGV=0.D0
      DO 231 II=1,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      IF(IID.EQ.1)GO TO 231
C=====TURN OFF STORED EDGE FIELD STRENGTH AS A TEST.
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
C     IF(DANGH.LT.1.D-4.AND.DABS(ANGV).LT.1.D-8)GO TO 2311
      WRITE(53,2312)ANGH,ANGV
 2312 FORMAT(' ','Total bending angles', 2F15.10)
C     STOP
 2311 CONTINUE
C
C=====1)GET SPIN TUNE FROM TRACE OF MATRIX. 2)EFFECT OF ARGUS ON TUNE.
      STUNE=(ROT(1,1)+ROT(2,2)+ROT(3,3)-1.D0)*0.5D0
      STUNE=DACOS(STUNE)/(2.D0*PI)
      RNU=DACOS(DCOS(0.01D0)*DCOS(PI*NU))/PI-1.D0
C=====3)CHECK ORTHOGONALITY OF THE ROTATION MATRIX
      CALL ORTCHK(ROT)
C
C
C
C=====GET EIGENVECTORS & ORDER THEM=====================================
C=====The old G.R. routines fail for 1/2 integer spin tune. So use NAG
C     CALL EIV3(ROT,ZR3,ZI3,WR3,WI3,IERRO)
C      CALL UCOPY(ROT  ,TRIN3,18)
      TRIN3 = ROT
      IFAIL=0
C      CALL F02AGF(TRIN3,3,3,RR3,RI3,VR3,3,VI3,3,INTGE3,IFAIL)
      CALL F02EBF('V',3,TRIN3,3,RR3,RI3,VR3,3,VI3,3,WORK,LWORK,IFAIL) ! <<replmnt
      IF(IFAIL.NE.0)GO TO 9999
C=====WRITE OUT F02AGF RESULTS
C     WRITE(53,926)
  926 FORMAT(' ','NOW THE F02AGF RESULTS: 3X3')
      DO 284 I=1,3
      WR3(I)=RR3(I)
      WI3(I)=RI3(I)
      DO 285 J=1,3
      ZR3(J,I)=VR3(J,I)
      ZI3(J,I)=VI3(J,I)
  285 CONTINUE
  284 CONTINUE
      TUN=999999.
CDB      IF(DABS(RR3(IT)) .LE. 1.D0)TUN=DACOS(RR3(IT))/(2.*PI)
C=====GET TUNES: THIS DOES NOT DIFFERENTIATE BETWEEN +/- ANGLES
C=====WITHOUT USING THE SINE.
C     WRITE(53,922)RR3(1 ),RI3(1 ),(VR3(J,1 ),J=1,3),TUN
C     WRITE(53,922)RR3(1 ),RI3(1 ),(VI3(J,1 ),J=1,3),TUN
C     WRITE(53,922)RR3(2 ),RI3(2 ),(VR3(J,2 ),J=1,3),TUN
C     WRITE(53,922)RR3(2 ),RI3(2 ),(VI3(J,2 ),J=1,3),TUN
C     WRITE(53,922)RR3(3 ),RI3(3 ),(VR3(J,3 ),J=1,3),TUN
C     WRITE(53,922)RR3(3 ),RI3(3 ),(VI3(J,3 ),J=1,3),TUN
  922 FORMAT(2F12.5,2X,'--->  (',3F10.5,' )',F15.8)
      NN=4
      ICOUNT=0
      DO 255 I=1,3
      TM3A(I,1)=WR3(I)
      TM3A(I,2)=WI3(I)
C=====Locate the real unit eigenvalue:
      IF(DABS(WI3(I)).LT.1.D-9.AND.WR3(I).GT.0.9999D0)ICOUNT=ICOUNT+1
 255  IF(DABS(WI3(I)).LT.1.D-9.AND.WR3(I).GT.0.9999D0)NN=I
      IF(NN .NE. 4)GO TO 241
      WRITE(53,242)
 242  FORMAT(' ','No real unit spin eigenvalue found----so STOP')
      STOP
C
 241  CONTINUE
      IF(ICOUNT.EQ.1)GO TO 243
C     NN=1
 243  CONTINUE
      MM=MOD(NN,3)+1
      LL=MOD(MM,3)+1
      AZ1=DSQRT(ZR3(1,NN)**2+ZR3(2,NN)**2+ZR3(3,NN)**2)
      AZ2=DSQRT(ZR3(1,MM)**2+ZR3(2,MM)**2+ZR3(3,MM)**2)
      DO 256 I=1,3
      ZW(I,1)=ZR3(I,NN)/AZ1
 256  ZW(I,2)=ZR3(I,MM)/AZ2

C      Force n_0 to be vertical.
C      ZW(1,1)=0.D0
C      ZW(2,1)=1.D0
C      ZW(3,1)=0.D0
C      Force m_0 to be radial.
C      ZW(1,2)=1.D0
C      ZW(2,2)=0.D0
C      ZW(3,2)=0.D0

C      Force the m vector to be vertical.
C      ZW(1,2)=0.D0
C      ZW(2,2)=1.D0
C      ZW(3,2)=0.D0
      ZW(1,3)=ZW(2,1)*ZW(3,2)-ZW(3,1)*ZW(2,2)
      ZW(2,3)=ZW(3,1)*ZW(1,2)-ZW(1,1)*ZW(3,2)
      ZW(3,3)=ZW(1,1)*ZW(2,2)-ZW(2,1)*ZW(1,2)
C      Force the l vector to be vertical.
C      ZW(1,3)=0.D0
C      ZW(2,3)=1.D0
C      ZW(3,3)=0.D0
C      ZW(1,2)=-ZW(2,1)*ZW(3,3)+ZW(3,1)*ZW(2,3)
C      ZW(2,2)=-ZW(3,1)*ZW(1,3)+ZW(1,1)*ZW(3,3)
C      ZW(3,2)=-ZW(1,1)*ZW(2,3)+ZW(2,1)*ZW(1,3)
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
C
C
      WRITE(53,933)NU,STUNE,RNU
      DO 250 I=1,3
  250 WRITE(53,932)(ROT(I,J),J=1,3),WR3(I),WI3(I),(ZW(J,I),J=1,3)
C
C
C
C
C
C
C    ***************************************************************
C    * ORTHONORMAL SPIN BASE VECTORS AND 8X8 MATRICES FOR ONE TURN *
C    ***************************************************************
C
C
      IPRIN1=IPTNML.EQ.1.AND.IE0.EQ.1
      IF(IPRIN1)WRITE(53,931)
C
C
C
C
C
      S=0.
      ISOL=0
      ANGH=0.D0
      CALL UNIT(8,SOL)
      CALL UNIT(8,TREV8)
C      CALL UCOPY(ZW,ROT,18)
      ROT = ZW
C      CALL UCOPY(ZW,TM3C,18)
      TM3C = ZW
C=====LOOP AROUND THE LATTICE===========================================
      DO 257 II=1,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)


C=====GET THE LENGTH and bend angle:
      IF(IID.EQ.2.OR.IID.EQ.15)ANGH=ANGH+XY
      IF(IID.EQ.5.OR.IID.EQ.6.OR.IID.EQ.7)GO TO 699
      IF( (NAME(ITY)(1:2).EQ.'CQ'.AND.IID.EQ.3) !Ignore artificial lengths.
     + .OR.(NAME(ITY)(1:2).EQ.'RQ'.AND.IID.EQ.4).OR.IID.EQ.17)GO TO 699
       S=S+YY(ITY)
  699 CONTINUE

C      CALL UCOPY(TM3C,ZW,18)
      ZW = TM3C

      IF(IID.EQ.1)THEN                     !Treat drifts separately.
      IF(XX(ITY).LT.0.00011D0)GO TO 2588   !Skip zero length drifts
      TREV8(1,:) = TREV8(1,:) + XY * TREV8(2,:)
      TREV8(3,:) = TREV8(3,:) + XY * TREV8(4,:)
      GO TO 2588
      ENDIF

      IF(IID.EQ.10)GO TO 1001
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      CALL JAM333(TM3C,TM3A,ZW)
      CALL AVER(ZW,TM3C,3,TM3B)

      IF(MRDISP.EQ.0)
     +CALL MX88(IID,II,ITY,TMAT(1,1,ITY),XY,XX2,YYY,NU,TM3B,ZW,TM8A)
      IF(MRDISP.EQ.1)
     +CALL DX88(IID,II,ITY,TMAT(1,1,ITY),XY,XX2,YYY,NU,TM3B,ZW,TM8A)
      GO TO 258
 1001 CONTINUE
      CALL SOL8AN(II,XY,YY(ITY),NU,ZW,TM3C,TM8A,NSOL(ITY))
C
C
  258 CONTINUE

C      WRITE(533,2581)NAME(ITY),S
C 2581 FORMAT(///,' 8X8 TRANSF.MATRIX at ',  A8, F15.9, /)
C     DO 2735 I=1,8
C 2735 WRITE(533,925)(TM8A(I,J),J=1,8)

C      CALL JAM888(TREV8,TM8A,TREV8)
      TREV8 = MATMUL(TM8A,TREV8)

 2588 CONTINUE

C
C
C
C     IF((NAME(ITY).EQ.'IP'.OR.II.LT.1000.OR.II.GT.6930).AND.IPRIN1)
      IF(                     IPRIN1 )
C      IF((NAME(ITY).EQ.'IP'.AND.IPRIN1))
     +                               WRITE(53,930)S,NAME(ITY),IID,
     +                (ZW(J,1),J=1,3),(ZW(J,2),J=1,3),(ZW(J,3),J=1,3),II
      IF(NAME(ITY).EQ.'IP6'.or.NAME(ITY).EQ.'IP12')THEN
      WRITE(53,9277)S
      WRITE(533,9277)S
 9277 FORMAT(///,' 8X8 TRANSF.MATRIX at',  F10.4, /)
      DO 2733 I=1,8
      WRITE(533,925)(TREV8(I,J),J=1,8)
 2733 WRITE(53,925)(TREV8(I,J),J=1,8)
      ENDIF

C     WRITE(534,9278)NAME(ITY),S
C 9278 FORMAT(///,' 8X8 TRANSF.MATRIX at   ',  A8,  F15.9, /)
C      DO 2734 I=1,8
C 2734 WRITE(534,925)(TREV8(I,J),J=1,8)

C
C
C
C
C=====SUM THE n0-AXIS DEVIATIONS OVER THE NON-VERTICAL BEND REGION:
C     THE ARCS,THE MICRO BETA REGION & THE SUM OF THESE.
C
  267 ZZ2=1.D0-ZW(2,1)*ZW(2,1)
      ZZT=1.D0-ZW(3,1)*ZW(3,1)
      IF(ZZ2.LT.0.D0)ZZ2=0.D0
      IF(ZZT.LT.0.D0)ZZT=0.D0

C      KTILT=0
C     IF((S.LT.12.0).OR.(S.GT.132.0.AND.S.LT.156.0).OR.(S.GT.276.0))
C    +KTILT=1
C     IF((S.GT.20.0.AND.S.LT.124.).OR.(S.GT.164.AND.S.LT.268))
C    +KTILT=2


C=====FOR HERA
      KTILT=0
C     IF(S.LT.25.  )KTILT=1                            ! For old HERA
C     IF(S.LT.2.D0 )KTILT=1                            ! For eRHIC
C      IF(S.LT.0.5.OR. S.GT.6335.8                      ! For the LHeC
C     +    .OR.DABS(S-1583.96).LT.0.5
C     +    .OR.DABS(S-4751.86).LT.0.5 )KTILT=1

      IF(S.LT.398.23 )KTILT=1                            ! For MEIC


C      IF(S.GT.230.0.AND.S.LT.792.0)KTILT=2            ! For HERA
C      IF(S.GT.100.0.AND.S.LT.1177.0)KTILT=2           ! For eRHIC
C      KTILT = 2                                       ! For a simple flat ring.
C      IF(S.GT.2000.D0.AND.S.LT.24000.D0)KTILT=2        ! For the LHeC

      IF(S.GT.398.23.AND.S.LT.670.22)KTILT=2        ! For MEIC

      IF(KTILT.EQ.1)THEN
      SUMNZM=SUMNZM+DSQRT(ZZ2)
      NNZM=NNZM+1
      SUMNTM=SUMNTM+DSQRT(ZZT)
      NNTM=NNTM+1
      ENDIF
      IF(KTILT.EQ.2.AND.IID.EQ.3)THEN
      SUMNZA=SUMNZA+DSQRT(ZZ2)
      NNZA=NNZA+1
      ENDIF
      IF(KTILT.EQ.0)GO TO 257
      SUMNZ=SUMNZ+DSQRT(ZZ2)
      NNZ=NNZ+1



C======End of lattice loop
  257 CONTINUE
C
C
C======WIND BACK THE SPIN BASIS==========================================
      CALL UNIT(8,TM8A)
      TM8A(7,7)= WR3(2)
      TM8A(7,8)= WI3(2)
      TM8A(8,7)=-WI3(2)
      TM8A(8,8)= WR3(2)

C======WIND FORWARD THE SPIN BASIS: This gives crazy polarization
C      TM8A(7,7)= WR3(2)
C      TM8A(7,8)=-WI3(2)
C      TM8A(8,7)= WI3(2)
C      TM8A(8,8)= WR3(2)

C      CALL UCOPY(TREV8,TREV08,128)             !Save the ``unrotated'' matrix.
      TREV08 = TREV8
C      CALL JAM888(TREV8,TM8A,TREV8)
      TREV8 = MATMUL(TM8A,TREV8)
      WRITE(53,927)
  927 FORMAT(///,' 8X8 TRANSF.MATRIX FOR ONE REVOLUTION AROUND THE',
     +                                        ' 1-ST BEAM-LINE ELEM:',/)
      DO 273 I=1,8
  273 WRITE(53,925)(TREV8(I,J),J=1,8)
  925 FORMAT(8E16.8)
      COL1=SQRT(TREV8(7,1)**2+TREV8(8,1)**2)
      COL2=SQRT(TREV8(7,2)**2+TREV8(8,2)**2)
      COL3=SQRT(TREV8(7,3)**2+TREV8(8,3)**2)
      COL4=SQRT(TREV8(7,4)**2+TREV8(8,4)**2)
      COL5=SQRT(TREV8(7,5)**2+TREV8(8,5)**2)
      COL6=SQRT(TREV8(7,6)**2+TREV8(8,6)**2)
      WRITE(53,1908)COL1,COL2,COL3,COL4,COL5,COL6
 1908 FORMAT(' ','Abs. values of G-matrix columns:', /6D16.8)
C

      CALL SYMP8(TREV8)
C
C
C=====GET 8X8 EIGENVECTORS FOR ONE REVOLUTION & ORDER THEM==============
C=====The old G.R. routines fail for 1/2 integer spin tune. So use NAG
C     CALL EIV8(TREV8,ZZ,WR8,WI8,IERRO)
C     IF(IERRO .NE. 0)GO TO 9999
C      CALL UCOPY(TREV8,TRIN8,128)
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
C=====GET TUNES: THIS DOES NOT DIFFERENTIATE BETWEEN +/- ANGLES
C=====WITHOUT USING THE SINE.
C     WRITE(53,923)RR(IT),RI(IT),(VR(J,IT),J=1,8),TUN
C     WRITE(53,923)RR(IT),RI(IT),(VI(J,IT),J=1,8),TUN
C
C
C
C
C
C
C
      DO 290 I=1,8
      SUM=0.
      DO 291 J=1,6
  291 SUM=SUM+ZZ(J,I)**2
      IF(SUM.LT.1.D-12)IS=I
      TM8B(I,1)=WR8(I)
  290 TM8B(I,2)=WI8(I)
C
C      CALL UCOPY(ZZ,TM8A,128)
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
  924 FORMAT(////,' EIGENVALUES AND EIGENVECTORS:',/,T8,'REAL',
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
      CALL NORM(ZZ,8,AB)   !Needed because we don't use directly the mormalised eigenvevtors
C                          !but forms such as for XY2 below. E.g. XY2 is pure real and only
C                          !half of the size that one gets by multiplying complex eigenvector
C                          !components.
      IF(AB(1).EQ.0.D0.OR.AB(3).EQ.0.D0.OR.AB(5).EQ.0.D0)THEN
      WRITE(53,'(A)')' Normalisation of 8-eigenvectors crazy.'
      NCRAZY = 1
      RETURN
      ENDIF
C
C
C======Get scaled 1-turn integrals:  multiply the un-rotated G-matric (in TREV08) onto the
C      orbit eigenvectors.
C      Actually, the absolute values are independent of the last spin basis rotation.
C      Off resonance, the values depend on the azimuth.
C      On resonance they are indep of azimuth and they are prop. to resonance strengths.
C      If a 1-turn integral from some azimuth is zero at all energies due to spin matching
C      then it is zero at resonance at ALL azimuths.
C      Resonance strengths can be got by doing inf. turn integrals in very(!) small energy steps
C      to avoid missing the non-zero spikes. Here, the approx. resonance strengths could be got by
C      choosing the energies roughly and using their smoothness. That's how SPRINT does it.
C
C
      DO  150 IR=1,2
      DO  151 KR=1,6
      RESSTR(IR,KR)=0.D0
      DO  152 JR=1,6
  152 RESSTR(IR,KR)=RESSTR(IR,KR) + TREV08(6+IR,JR)*ZZ(JR,KR)
  151 CONTINUE
  150 CONTINUE

C======Get the real and imag parts for the plus and minus orbital tunes
C      Use the prper normalisation for th eeigenvectors.
      R1RP = (RESSTR(2,1)-RESSTR(1,2))/DSQRT(DABS(AB(1)))
      R1IP = (RESSTR(1,1)+RESSTR(2,2))/DSQRT(DABS(AB(1)))
      R1RM = (RESSTR(2,1)+RESSTR(1,2))/DSQRT(DABS(AB(1)))
      R1IM = (RESSTR(1,1)-RESSTR(2,2))/DSQRT(DABS(AB(1)))

      R2RP = (RESSTR(2,3)-RESSTR(1,4))/DSQRT(DABS(AB(3)))
      R2IP = (RESSTR(1,3)+RESSTR(2,4))/DSQRT(DABS(AB(3)))
      R2RM = (RESSTR(2,3)+RESSTR(1,4))/DSQRT(DABS(AB(3)))
      R2IM = (RESSTR(1,3)-RESSTR(2,4))/DSQRT(DABS(AB(3)))

      R3RP = (RESSTR(2,5)-RESSTR(1,6))/DSQRT(DABS(AB(5)))
      R3IP = (RESSTR(1,5)+RESSTR(2,6))/DSQRT(DABS(AB(5)))
      R3RM = (RESSTR(2,5)+RESSTR(1,6))/DSQRT(DABS(AB(5)))
      R3IM = (RESSTR(1,5)-RESSTR(2,6))/DSQRT(DABS(AB(5)))

      R1P  = DSQRT(R1RP**2 + R1IP**2)
      R1M  = DSQRT(R1RM**2 + R1IM**2)

      R2P  = DSQRT(R2RP**2 + R2IP**2)
      R2M  = DSQRT(R2RM**2 + R2IM**2)

      R3P  = DSQRT(R3RP**2 + R3IP**2)
      R3M  = DSQRT(R3RM**2 + R3IM**2)
C
C
      WRITE(53,'(A,6(1X,F15.5))')' Generalised resonance strengths',
     +                          R1P,R1M,R2P,R2M,R3P,R3M

      WRITE(53,234)
  234 FORMAT('0','Identify resonances:          N1    N2    N3    INT')
      CALL RESON(E0,TN(1),TN(3),TN(5),TN(7))
C
C
C
C
C
      IF(MRDISP.EQ.0)GO TO 765
C=====IF USING THE DISPERSION VERSION GET THE UNCOUPLED TWISSES.
      RMUX=DACOS(0.5D0*(TREV8(1,1)+TREV8(2,2)))
      IF(TREV8(1,2) .LT. 0.)RMUX=2.D0*PI-RMUX
      RMUY=DACOS(0.5D0*(TREV8(3,3)+TREV8(4,4)))
      IF(TREV8(3,4).LT.0.D0)RMUY=2.D0*PI-RMUY
      BETAX=TREV8(1,2)/DSIN(RMUX)
      BETAY=TREV8(3,4)/DSIN(RMUY)
      ALPHAX=(TREV8(1,1)-TREV8(2,2))/(2.*DSIN(RMUX))
      ALPHAY=(TREV8(3,3)-TREV8(4,4))/(2.*DSIN(RMUY))
      WRITE(53,763)
  763 FORMAT('0','Twisses from uncoupled 8X8 matrix:',/
     +'     Betax     Betay    Alphax     Alphay')
      WRITE(53,762)BETAX,BETAY,ALPHAX,ALPHAY
  762 FORMAT(' ',4F10.4)
  765 CONTINUE
C
C
C
C
C
C=====FOR GETTING D/K  N AXIS TILT:
C=====SINCE THE 6X6 & 8X8 EIGENVALUE LISTS NEED NOT COME IN
C=====THE SAME ORDER, CHOOSE THE CORRECT EMITTANCES.
C=====THE 4TH EIGEN VECTOR PAIR HAS JUST BEEN DEFINED TO BE THE
C====='PURE SPIN' EIGEN VECTOR. ONLY IF PERTOP HAS BEEN CALLED.
      IF(IPTCO.EQ.2)GO TO 26
      WRITE(53,233)
  233 FORMAT('0','Check the matching of the 6x6 and 8x8 eigenvalues')
      NTUN=0
      DO 25 JJ=1,5,2
      DO 24 J =1,5,2
C     DIF=DABS(TN(JJ)-PTUN(J))
      DIF=DABS(DABS(TN(JJ))-DABS(PTUN(J)))
      IF(DIF.LT.0.0010)AASPIN(JJ)=AAN(J)
      IF(DIF.LT.0.0010)NTUN=NTUN+1
      IF(DIF.LT.0.0010)WRITE(53,23)JJ,J
   23 FORMAT(' ','JJ,J',2I10)
   24 CONTINUE
   25 CONTINUE
      IF(NTUN.NE.3)GO TO 9988
   26 CONTINUE
C
C
C
C
C
C      IF(1.EQ.1)RETURN
C
C
C
C
C
C
C     ****************************************************
C     *      THE POLARIZATION AND DEPOLARIZATION.        *
C     ****************************************************
C
C
      IPRIN2=IPTD.EQ.1.AND.IE0.EQ.1
      IF(IPRIN2)WRITE(53,936)
C
C
C
C=====PURE SOKOLOV/TERNOV POLARISATION TERMS:DENOMINATOR,NUMERATOR RESP.
      WP0  =0.D0
      WP1  =0.D0
C=====11/18*(D**2) DEPOLARISATION TERMS.
      WD0  =0.D0
      WD0A =0.D0
      WD0B =0.D0
      WD0AB=0.D0
      WD1  =0.D0
      WD2  =0.D0
      WD3  =0.D0
      WD4  =0.D0
C=====NUMERATOR TERM INCLUDING 'D' DEPOL. TERM.
      WND  =0.D0
C
C
C
C
C
C
      NNAXIS=0
      ENAXIS=0.D0
      AXIST=0.D0
      ANGH=0.D0
C      CALL UCOPY(ZZ,TM8B,128)
      TM8B = ZZ
      S=0.
      DO 300 II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY =XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)

C=====GET the Length and bend angle:
      IF(IID.EQ.2.OR.IID.EQ.15)ANGH=ANGH+XY
      IF(IID.EQ.5.OR.IID.EQ.6.OR.IID.EQ.7)GO TO 700
      IF((NAME(ITY)(1:2).EQ.'CQ'.AND.IID.EQ.3) !Ignore artificial lengths.
     +  .OR.(NAME(ITY)(1:2).EQ.'RQ'.AND.IID.EQ.4).OR.IID.EQ.17)GO TO 700
       S=S+YY(ITY)
  700 CONTINUE

C      CALL UCOPY(TM8B,TREV8,128)
      TREV8 = TM8B
C      CALL UCOPY(ROT,TM3B,18)
      TM3B = ROT

      IF(IID.EQ.1)THEN                     !Treat drifts separately
      IF(XX(ITY).LT.0.00011D0)GO TO 300    !Skip zero length drifts
      TM8B(1,:) = TREV8(1,:) + XY * TREV8(2,:)
      TM8B(3,:) = TREV8(3,:) + XY * TREV8(4,:)
      GO TO 300
      ENDIF

      IF(IID .EQ. 10)GO TO 301
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      CALL JAM333(ROT,TM3A,TM3B)
      CALL AVER(ROT,TM3B,3,ZW)

      IF(MRDISP.EQ.0)
     +CALL MX88(IID,II,ITY,TMAT(1,1,ITY),XY,XX2,YYY,NU,ZW,TM3B,TM8A)
      IF(MRDISP.EQ.1)
     +CALL DX88(IID,II,ITY,TMAT(1,1,ITY),XY,XX2,YYY,NU,ZW,TM3B,TM8A)
      GO TO 202
  301 CALL SOL8AN(II,XY,YY(ITY),NU,TM3B,ROT,TM8A,NSOL(ITY))
C  202 CALL JAM888(TM8B,TM8A,TREV8)
  202 TM8B = MATMUL(TM8A,TREV8)

  211 CONTINUE
C
C
C
C=====GET 1/(RHO**3) INTEGRANDS INCLUDING C.O. DEVIATION EFFECTS
      IDIP=0
      DELX=(DX(II)+DX(II+1))/2.
      DELY=(DY(II)+DY(II+1))/2.
      XY1=0.D0
      SOLFLD=0.D0
      CURVAT=0.D0
      GO TO(311,302,303,304,300,302,302,308,302,310,300,300,300,300,
     +                                                      302,302),IID
      GO TO 300

  311 CONTINUE
      IF(II.NE.2)GO TO 300
      GO TO 210
C
C
  302 XY1=DABS(XY)**3/YY(ITY)**2
      CURVAT=XY/YY(ITY)*10.D0
C      IF(NSOL(ITY).EQ.20)XY1=0.D0

C=====REMOVE SYNCHROTRON RADIATION FROM SOME DIPOLES, ADD BY FANGLEI 04 JAN 2013========
C      IF(NAME(ITY)(1:2).EQ.'BS')XY1=0.D0
C      IF(NAME(ITY)(1:3).EQ.'BS1')XY1=0.D0
C      IF(NAME(ITY)(1:4).EQ.'BSR1')XY1=0.D0
C      IF(NAME(ITY)(1:3).EQ.'BS2')XY1=0.D0
C      IF(NAME(ITY)(1:4).EQ.'BSR2')XY1=0.D0

      IDIP=1
C      IF(IID.NE.6.AND.IID.NE.7)THEN
C      WRITE(53,'(A,2I10,2A,E16.5)')
C     +                        ' ',II,IID,'  ',NAME(ITY),XY**2/YY(ITY)
C      ENDIF
C      IF(IID.EQ.2.OR.IID.EQ.15)XY1 = 0.28490E-05

      GO TO 210
C
C
  303 CONTINUE
C     IF(NTWIST(ITY).EQ.3)DELY=DELY-TWIST(ITY)
C     IF(NTWIST(ITY).EQ.4)DELX=DELX-TWIST(ITY)
      XY1=0.D0
      IF(NAME(ITY)(1:1).NE.'E')          !Exclude edge fields.
     +XY1=((XY*DELX)**2+(XY*DELY)**2)**1.5/YY(ITY)**2
      GO TO 210
C
  304 CONTINUE
      XY1=((XY*DELX)**2+(XY*DELY)**2)**1.5/YY(ITY)**2
      GO TO 210
C
  308 XY1=((XY*(DELX**2-DELY**2)/2.)**2+(XY*DELX
     +                                        *DELY)**2)**1.5/YY(ITY)**2
      GO TO 210
C
C
C=====As in December 1995: as in EMITNC turn off radiation in solenoids.
  310 CONTINUE
      XY1=0.D0
      SOLFLD=XY
C 310 CALL SOLXYP(II,XY/YY(ITY),VX,VY)
C     XY1=((XY*VX)**2+(XY*VY)**2)**1.5/YY(ITY)**2
C
C
  210 CONTINUE
C
      IF(MRDISP.EQ.0)CALL AVER(TREV8,TM8B,8,ZV)
C=====PUT IN HERE (TREV8)*(Dispersion vector) TO DEFINE ZV. NEED A SPECIAL AVER!
      IF(MRDISP.EQ.1)CALL DIAVER(TREV8,TM8B,ZV,II)
C     IF(MRDISP.EQ.1)CALL DIAVER(TREV8,TM8B,ZV,II,NAME(ITY))
      XY2=(-ZV(5,2)*ZV(7,1)+ZV(5,1)*ZV(7,2))/AB(1)
      XY3=(-ZV(5,4)*ZV(7,3)+ZV(5,3)*ZV(7,4))/AB(3)
      XY4=(-ZV(5,6)*ZV(7,5)+ZV(5,5)*ZV(7,6))/AB(5)
      XY5=(-ZV(5,2)*ZV(8,1)+ZV(5,1)*ZV(8,2))/AB(1)
      XY6=(-ZV(5,4)*ZV(8,3)+ZV(5,3)*ZV(8,4))/AB(3)
      XY7=(-ZV(5,6)*ZV(8,5)+ZV(5,5)*ZV(8,6))/AB(5)

C===Get phases of the synch. eigenvector as in EMITNC
      PHASE1 = DATAN2(TREV8(5,6),TREV8(5,5))/2.D0/pi
      PHASE2 = DATAN2(TREV8(6,6),TREV8(6,5))/2.D0/pi - 0.25D0

      IF(AB(1).GT.0.D0.OR.AB(3).GT.0.D0.OR.AB(5).GT.0.D0)THEN
      WRITE(53,'(A,A,3F10.4)')' Normalisation of 8-eigenvectors crazy:',
     +                   '  AB(1), AB(2), AB(3) -- ', AB(1),AB(3),AB(5)
      NCRAZY = 1
      RETURN
      ENDIF


C
C=====Stuff for analysing kinetic polarization in AmPs and Bates.
      XY51=ZV(5,1)/DSQRT(-AB(1))
      XY52=ZV(5,2)/DSQRT(-AB(1))
      XY512=DSQRT(XY51**2 + XY52**2)
      XY53=ZV(8,1)/DSQRT(-AB(1))
      XY54=ZV(8,2)/DSQRT(-AB(1))
      XY534=DSQRT(XY53**2 + XY54**2)
      QD5=TM8A(8,6)*XY512
C
      XY71=ZV(5,5)/DSQRT(-AB(5))
      XY72=ZV(5,6)/DSQRT(-AB(5))
      XY712=DSQRT(XY71**2 + XY72**2)
      XY73=ZV(8,5)/DSQRT(-AB(5))
      XY74=ZV(8,6)/DSQRT(-AB(5))
      XY734=DSQRT(XY73**2 + XY74**2)
      QD7=TM8A(8,6)*XY712

      HD=D1(II)
C      CURVAT=XY1*100.D0
C      CURVAT=(XY1)**(1.D0/3.D0)*100.D0
C
C
C      WRITE(53,967)II,NAME(ITY),XY71,XY72,XY73,XY74,QD7,HD
  967 FORMAT(' ',I5,A8,6E15.6)
C
C
C     WRITE(53,927)
C     DO 2273 I=1,8
C2273 WRITE(53,925)(ZV(I,J),J=1,8)
C
C
C=====GET N AXIS TILT AT THIS POINT.
C=====NB: WE WOULD NORMALLY DIVIDE BY 2.ABN (SEE SLIM NOTES). BUT
C=====WE NEED 2* EMITTANCE IN THE NUMERATOR. SO THE TWO'S CANCEL.
C=====MODIFIED MAY92 TO AGREE WITH CALCULATION IN SLIM:
C=====VIZ WEIGHT BY MAGNET LENGTH, INCLUDE ALL ELEMENTS EXCEPT DRIFTS,
C=====KEEP AS A MEAN SQ DEVN AND DIVIDE BY 2.
C=====Is this formalism OK in the dispersion formalism?
      IF(IPTCO.EQ.2)GO TO 27
      ENTILT=
C            +AASPIN(1)*(ZV(7,1)*ZV(7,1)+ZV(7,2)*ZV(7,2))/DABS(AB(1))
C     +      +AASPIN(1)*(ZV(8,1)*ZV(8,1)+ZV(8,2)*ZV(8,2))/DABS(AB(1))
C     +      +AASPIN(3)*(ZV(7,3)*ZV(7,3)+ZV(7,4)*ZV(7,4))/DABS(AB(3))
C     +      +AASPIN(3)*(ZV(8,3)*ZV(8,3)+ZV(8,4)*ZV(8,4))/DABS(AB(3))
     +      +AASPIN(5)*(ZV(7,5)*ZV(7,5)+ZV(7,6)*ZV(7,6))/DABS(AB(5))
     +      +AASPIN(5)*(ZV(8,5)*ZV(8,5)+ZV(8,6)*ZV(8,6))/DABS(AB(5))
      TNTILT=ENTILT
      ENTILT=ENTILT*YY(ITY)
      NNAXIS=NNAXIS+1
      ENAXIS=ENAXIS+YY(ITY)
C     AXIST=AXIST+DSQRT(ENTILT)
      AXIST=AXIST+     (ENTILT)
      SNTILT=DSQRT(TNTILT)*1000.D0
C     IF(II.LT.4000)WRITE(53,934)NAME(ITY),SNTILT
C=====THESE IF STATEMENTS DON'T WORK IN THIS POSITION IF DRIFTS
C=====ARE SKIPPED.
C     IF(NAME(ITY).EQ.'DRIFT')WRITE(53,934)NAME(ITY),SNTILT
C     IF(II.EQ.2        )WRITE(53,934)NAME(ITY),SNTILT
C                        WRITE(53,934)NAME(ITY),SNTILT
  277 CONTINUE
   27 CONTINUE
C
C
C=====DDM/DDL ARE 'D' PERTURBATION VECTORS ALONG M/L DIRECTIONS.
C=====AT THIS POINT "ZW" IS THE AVERAGE N,M,L BASIS--EXCEPT FOR SOLS.
C      XY3 = XY3/1.5D0
C      XY6 = XY6/1.5D0
C      XY4 = XY4/2.25D0      ! To compensate for increasing the spin phase advances X 1.5.
C      XY7 = XY7/2.25D0      ! Don't want to increase the synchrobeta coupling

      DDM  =-XY2-XY3-XY4
      DDL  =-XY5-XY6-XY7
C      DDML = XY2+XY5+XY3*XY6+XY4*XY7
      DDML = DDM * DDL

C      DDM=-XY4
C      DDL=-XY7
C
C      DDM=-XY2
C      DDL=-XY5

C      DDM=-XY3
C      DDL=-XY6

      DMX=ZW(1,2)*DDM
      DMY=ZW(2,2)*DDM
      DMZ=ZW(3,2)*DDM
      DLX=ZW(1,3)*DDL
      DLY=ZW(2,3)*DDL
      DLZ=ZW(3,3)*DDL
C      DDMABS=DSQRT(DMX**2 + DMZ**2)
C      DDMABS=DSQRT(DMX**2 + DMZ**2)
      DDABS=DSQRT(DDM**2 + DDL**2)
      DDLABS=DSQRT(DLX**2 + DLZ**2)
      DDMRUF=DABS(NU*(PI-ANGH)) - NU*PI
C      IF(NTWIST(ITY).GT.0.AND.IPRIN2.AND.IDIP.EQ.1)
C      DDLABS=DSQRT(DLX**2 + DLZ**2)
C      DDLRUF=NU*(PI-ANGH)
      IF(IPRIN2)
     +         WRITE(53,937)S,NAME(ITY),IID,
     +         (ZW(J,1),J=1,3),DMX,DMY,DMZ,DDM,DLX,DLY,DLZ,DDL
      IF(IPRIN2)WRITE(55,938)S,NAME(ITY),II,IID,
     +(ZW(J,1),J=1,3),DMX,DMY,DMZ,DDM,DLX,DLY,DLZ,DDL,
     +                       DDABS,DDLABS,DDMRUF,
     +                XY71,XY72,XY712,XY73,XY74,XY734,QD7,
     +                XY51,XY52,XY512,XY53,XY54,XY534,QD5,
     +                D1(II),D3(II),CURVAT,SOLFLD,DX(II),DY(II),
     +                BETACX(II),BETACY(II),PHASE1,PHASE2
C
C
C
C=====GET SUM 1/(RHO**3) * EIGENVECTOR PRODUCTS: TOTAL & INDIVIDUAL.
C      IF(NSOL(ITY).NE.20)WD0=WD0+XY1*(DDM**2+DDL**2)
      WD0  =WD0  +XY1*(DDM**2+DDL**2)
      WD0A =WD0A +XY1* DDM**2
      WD0B =WD0B +XY1* DDL**2
      WD0AB=WD0AB+XY1* DDML
C
C      IF(NSOL(ITY).NE.20)WD1=WD1+XY1*(XY2**2+XY5**2)
      WD1=WD1+XY1*(XY2**2+XY5**2)
C
      WD2=WD2+XY1*(XY3**2+XY6**2)
C
C      IF(NSOL(ITY).NE.20)WD3=WD3+XY1*(XY4**2+XY7**2)
C      IF(NSOL(ITY).EQ.20)WD3=WD3+XY1*(XY4**2+XY7**2)*100.
      WD3=WD3+XY1*(XY4**2+XY7**2)
C
      WD4=WD4+XY1*(XY2**2+XY5**2)
     +       +XY1*(XY3**2+XY6**2)
     +       +XY1*(XY4**2+XY7**2)
C
C
C
C=====GET VELOCITY COMPONENTS & PROJECTION OF SPIN ON VELOCITY.
C     IF(IID.NE.10)VX=(DXP(II+1)+DXP(II))/2.
C     IF(IID.NE.10)VY=(DYP(II+1)+DYP(II))/2.
      VX=(DXP(II+1)+DXP(II))/2.
      VY=(DYP(II+1)+DYP(II))/2.
      VZ=DSQRT(1.-VX*VX-VY*VY)
      VN=VX*ZW(1,1)+VY*ZW(2,1)+VZ*ZW(3,1)
      VM=VX*ZW(1,2)+VY*ZW(2,2)+VZ*ZW(3,2)
      VL=VX*ZW(1,3)+VY*ZW(2,3)+VZ*ZW(3,3)
      SGN=-1.D0
      IF(XY.GT.0.D0)SGN=1.D0
C
C
C
C=====GET    N *(V * VDOT)/VDOT FACTOR---->WP2
C=====AND (N-D)*(V * VDOT)/VDOT FACTOR---->WP3
C=====IF XY2 IS SMALL IN LENSES(I.E. FIELD IS SMALL)---THEN IGNORE.
      GO TO(300,402,403,404,300,402,407,408,407,410,209,209,209,209,
     +                                               402,407,300),IID
C
C
  402 WP2=(ZW(2,1)-VY*VN)/DSQRT(VX*VX+VZ*VZ)
      WP3=(ZW(2,1)-VY*VN-DDM*(ZW(2,2)-VY*VM)-DDL*(ZW(2,3)-VY*VL))/
     +                                                DSQRT(VX*VX+VZ*VZ)
      GO TO 209
C
C
  407 WP2=(ZW(1,1)-VX*VN)/DSQRT(VY*VY+VZ*VZ)
      WP3=(ZW(1,1)-VX*VN-DDM*(ZW(1,2)-VX*VM)-DDL*(ZW(1,3)-VX*VL))/
     +                                                DSQRT(VY*VY+VZ*VZ)
C=====CHANGE SIGNS TO AGREE WITH M & R USE OF CURVATURE (NOT FIELD)
      WP2=-WP2
      WP3=-WP3
      GO TO 209
C
C
  403 CONTINUE
C     IF(NTWIST(ITY).EQ.3)DELY=DELY-TWIST(ITY)
C     IF(NTWIST(ITY).EQ.4)DELX=DELX-TWIST(ITY)
      XY2=DSQRT((DELX*VZ)**2+(DELY*VZ)**2+(VX*DELX-VY*DELY)**2)
      XY3=0.
C=====TAKE CARE OF CASE WHERE THERE IS NO CLOSED ORBIT SHIFT.
      IF(XY2.GT.1.D-30)XY3=1./XY2
      WP2=DELY*ZW(1,1)+DELX*ZW(2,1)-(VN)*(DELY*VX+DELX*VY)
      WP3=DELY*ZW(1,1)+DELX*ZW(2,1)-(VN-DDM*VM-DDL*VL)*
     +                                             (DELY*VX+DELX*VY)
     +                              -DDM*(DELY*ZW(1,2)+DELX*ZW(2,2))
     +                              -DDL*(DELY*ZW(1,3)+DELX*ZW(2,3))
      WP2=WP2*XY3
      WP3=WP3*XY3
      GO TO 209
C
C
  404 XY2=DSQRT((DELY*VZ)**2+(DELX*VZ)**2+(DELY*VX+DELX*VY)**2)
      XY3=0.
      IF(XY2.GT.1.D-30)XY3=1./XY2
      WP2=DELY*ZW(2,1)-DELX*ZW(1,1)-(VN)*(DELY*VY-DELX*VX)
      WP3=DELY*ZW(2,1)-DELX*ZW(1,1)-(VN-DDM*VM-DDL*VL)*
     +                                             (DELY*VY-DELX*VX)
     +                              +DDM*(DELX*ZW(1,2)-DELY*ZW(2,2))
     +                              +DDL*(DELX*ZW(1,3)-DELY*ZW(2,3))
      WP2=WP2*XY3
      WP3=WP3*XY3
      GO TO 209
C
C
  408 XY4=DELX*DELY
      XY5=(DELX**2-DELY**2)/2.
      XY2=DSQRT((VZ*XY5)**2+(VZ*XY4)**2+(VX*XY5-VY*XY4)**2)
      XY3=0.
      IF(XY2.GT.1.D-30)XY3=1./XY2
      WP2=XY4*ZW(1,1)+XY5*ZW(2,1)-(XY4*VX+XY5*VY)*(VN)
      WP3=XY4*ZW(1,1)+XY5*ZW(2,1)-(XY4*VX+XY5*VY)*(VN-DDM*VM-DDL*VL)
     +                           -DDM*(XY4*ZW(1,2)+XY5*ZW(2,2))
     +                           -DDL*(XY4*ZW(1,3)+XY5*ZW(2,3))
      WP2=WP2*XY3
      WP3=WP3*XY3
      GO TO 209
C
C
C=====THIS SOLENOID TREATMENT IS NOT COMPLETE: SEPARATE ENDS & CENTRE.
  410 XY2=DSQRT(VX*VX+VY*VY)
      XY3=0.
      IF(XY2.GT.1.D-30)XY3=1./XY2
      WP2=ZW(3,1)-VZ*(VN)
      WP3=ZW(3,1)-VZ*(VN-DDM*VM-DDL*VL)-DDM*ZW(3,2)-DDL*ZW(3,3)
      WP2=WP2*XY3
      WP3=WP3*XY3
C
C
C
C=====GET NUM. & DENOM. INTEGRALS--->PURE SOKOLOV/TERNOV POLARISATION.
  209 CONTINUE
C     WP0I=XY1*(1.-2.*VN*VN/9.)
C     WP1I=SGN*XY1*WP2
      WP1=WP1+SGN*XY1*WP2
      WND=WND+SGN*XY1*WP3
      WP0=WP0+XY1*(1.D0 - 2.D0*VN*VN/9.D0)
C
C     IF(NAME(ITY).EQ.EDGE)WRITE(53,155)IID,S,XY1,WP2,WP1I,WP0I,DELX,DELY
C 155 FORMAT(' ','EDGE ID,S,XY1,WP2,WP1I,WP0I,DELX,DELY',I4,F7.,6E13.3)
C
C
C
C
C=====END OF LATTICE LOOP===============================================
  300 CONTINUE
C
C
C
C
C
C
C
C
C=====AVERAGE N-AXIS TILT TO COMPARE WITH SMILE3
C     IF(NNAXIS.NE.0)AXIST=AXIST/NNAXIS
      IF(ENAXIS.NE.0D0)AXIST=AXIST/ENAXIS
C     AXIST=     (AXIST)*1000000.D0           *0.5D0
      AXIST=DSQRT(AXIST)*1000.D0
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
C
C
C
C
C
C
C=====USE THE INTEGRALS TO GET THE POLARISATION.
      SYNCON=3.534D-19*(E0/0.000511D0)**5/CIR
C      SYNCON=1.D0
      WP0=WP0*SYNCON
      WP1=WP1*SYNCON
      WND=WND*SYNCON
      TAUP=1.D0/WP0
      P0=DABS(92.38*WP1/WP0)
      WRITE(53,943)WP1,WP0,P0,TAUP
C====Store for SCRURITA3
      SPINP0   = P0
      SPINTAUP = TAUP
C
C=====GET TOTAL & INDIVIDUAL 11/18*(D**2) TERMS.
      WD0UN=WD0
      WD0=WD0*11.D0/18.D0*SYNCON
      IF(WD0.LT.1.D-30)WD0=1.D-30
      TAUD0=1./WD0
      WD0UN=WD0UN/(1.D0+WD0UN/16.D0)**2
      WD0UN=WD0UN*11.D0/18.D0*SYNCON
C
C
      WD1=WD1*11.D0/18.D0*SYNCON
      IF(WD1.LT.1.D-30)WD1=1.D-30
      TAUD1=1./WD1
C
C
      WD2=WD2*11.D0/18.D0*SYNCON
      IF(WD2.LT.1.D-30)WD2=1.D-30
      TAUD2=1./WD2
C
C
      WD3=WD3*11.D0/18.D0*SYNCON
      IF(WD3.LT.1.D-30)WD3=1.D-30
      TAUD3=1./WD3
C
C
      WD4=WD4*11.D0/18.D0*SYNCON
      IF(WD4.LT.1.D-30)WD4=1.D-30
      TAUD4=1./WD4
C
C
C=====ADD WD0 INTO SOKOLOV/TERNOV POLN. TO GET TOTAL
      POL1= DABS(92.38*WP1/(WP0+WD0))
      POLM1=DABS(92.38*WP1/(WP0+WD1))
      POLM2=DABS(92.38*WP1/(WP0+WD2))
      POLM3=DABS(92.38*WP1/(WP0+WD3))
      POLM4=DABS(92.38*WP1/(WP0+WD4))
      WRITE(53,944)POL1,E0
      POL1UN=DABS(92.38*WP1/(WP0+WD0UN))
      WRITE(53,945)POL1UN,E0
      POL2=DABS(92.38*WND/(WP0+WD0))
      WRITE(53,946)POL2,E0
      EMPOL2= POL2
C
C
      XXX0=TAUP/TAUD0
      XXX1=TAUP/TAUD1
      XXX2=TAUP/TAUD2
      XXX3=TAUP/TAUD3
      XXX4=TAUP/TAUD4
      WRITE(53,942)POL1, TAUD0,XXX0
      WRITE(53,947)POLM1,TAUD1,XXX1
      WRITE(53,948)POLM2,TAUD2,XXX2
      WRITE(53,949)POLM3,TAUD3,XXX3
      WRITE(53,950)POLM4,TAUD4,XXX4

      IF(NNZ.NE.0)SUMNZ =SUMNZ/NNZ
      SUMNZ =     (SUMNZ)*1000.D0
      IF(NNZM.NE.0)SUMNZM=SUMNZM/NNZM
      SUMNZM=     (SUMNZM)*1000.D0
      IF(NNTM.NE.0)SUMNTM=SUMNTM/NNTM
      SUMNTM=     (SUMNTM)*1000.D0
      IF(NNZA.NE.0)SUMNZA=SUMNZA/NNZA
      SUMNZA=     (SUMNZA)*1000.D0


      WRITE(54,9466)E0,NU,
     +              P0,POL1,POLM1,POLM2,POLM3,
     +              EMPOL2,DMY,
     +              R1RP,R1IP,R1RM,R1IM,
     +              R2RP,R2IP,R2RM,R2IM,
     +              R3RP,R3IP,R3RM,R3IM,
     +              TAUP,TAUD0,TAUD1,TAUD2,TAUD3,
     +              SUMNZ,SUMNZM,SUMNTM,SUMNZA,AXIST,
     +               R1P,R1M,R2P,R2M,R3P,R3M
 9466 FORMAT(' ',40(E15.7,3X))

C      WRITE(54,9466)E0,NU,TN(2),TN(4),TN(6),TN(7),
C     +              P0,POL1,TAUP,TAUD0,TAUD1,TAUD2,TAUD3
C 9466 FORMAT(' ',15(E16.6,3X))
C
C
C
C
C
C
C
C      WRITE (53,980)SUMNZ
      WRITE (53,983)SUMNTM
      WRITE (53,981)SUMNZA
C      WRITE (53,982)SUMNZM
      WRITE (53,984)AXIST

C======Get the spin ellipse.
      WRITE(53,103)
      WRITE(53,103)
  103 FORMAT(/,'  ')
      IF(WD0B.EQ.0.D0)WD0B=1.D-10
      WD0ABRATIO = WD0A/WD0B
      WRITE(53,'(A,A,3E12.3,A,E12.3)')'0',
     +        'Spin alpha-beta ellipse parameters: ',WD0A,WD0B,WD0AB,
     +'  Ratio of alpha,betas axes: ',  WD0ABRATIO

      IF(IE0.LT.-1)THEN                  ! Do this ellipse stuff in SCRURITA2
C      IF(IE0.EQ.1)THEN
C======Get major and minor axes of ``1-turn'' spin ellipse from WD0A,WD0B,WD0AB.
      TILT=DATAN2(2.D0*WD0AB,WD0A - WD0B)                 ! 2 theta
      CTILT=DCOS(TILT)                                    !cos (2 theta)
      TILT=TILT*90.D0/PI                                  !Just 1 theta
      AB3=WD0A*WD0B - WD0AB**2
      IF(AB3.LT.0.D0.AND.AB3.GT.-1.D-8)AB3=0.D0         !Why are these gymnastics needed?
      AXISMAJ  =  DABS(AB3/( 0.5D0*((WD0B - WD0A)/CTILT + WD0B + WD0A)))
      AXISMAJ  =  DSQRT(AXISMAJ)
      AXISMIN  =  DABS(AB3/(-0.5D0*((WD0B - WD0A)/CTILT - WD0B - WD0A)))
      AXISMIN  =  DSQRT(AXISMIN)
C======Plot the ellipse at the starting point. Take 1000 points. Also principle axes.
      DLTA  = TILT*PI/180.D0
      DO 70 IPH = 1,1000
      PH    = 0.001D0*2.D0*PI*(IPH-1)
      HR    = AXISMAJ*DCOS(PH)
      VR    = AXISMIN*DSIN(PH)
      HOR   = HR*DCOS(DLTA) - VR*DSIN(DLTA)
      VER   = HR*DSIN(DLTA) + VR*DCOS(DLTA)
      AXMAJ =    (-1.D0 + 0.002D0*(IPH-1))*DTAN(DLTA)
      AXMIN =    (-1.D0 + 0.002D0*(IPH-1))*DTAN(DLTA+PI/2.D0)
      WRITE(65,1105)IPH,HOR,VER,-1.D0+0.002D0*(IPH-1),AXMAJ,AXMIN
 1105 FORMAT(' ',I10, 6E16.6)
   70 CONTINUE
      ENDIF


C
C
C
C
C
C=====GO ONCE MORE  ROUND  THE  RING  TO  EXTRACT  SPIN  PHASE ADVANCES
C=====FROM THE 3-ELEMENT EIGEN-VECTORS----C.F. ORBITAL PHASE
C=====CALCULATION IN PERTOP
C
      IF(IE0.GT.1)RETURN
      DO 1055 J=1,3
      A(J,2)=ZR3(J,MM)
      B(J,2)=ZI3(J,MM)
      A(J,3)=ZR3(J,LL)
 1055 B(J,3)=ZI3(J,LL)
      DO 105 I=2,3
      DO 105 J=1,3
      P(J,I)=0.D0
      IF(A(J,I).NE.0.D0)P(J,I)=DATAN2(B(J,I),A(J,I))
      IF(B(J,I).LT.0.D0.AND.A(J,I).NE.0.D0)P(J,I)=PI2+P(J,I)
      PTEMP(J,I)=P(J,I)
  105 CONTINUE
C
C
      S=0.D0
C     WRITE(53,1107)
C      CALL VZERO(PTOT,18)
      PTOT = 0.D0
      DO 104 II=1,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
      IF(IID.NE.5.AND.IID.NE.6.AND.IID.NE.7.AND.IID.NE.17)S=S+YY(ITY)
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      CALL JAM333(A,TM3A,A)
      CALL JAM333(B,TM3A,B)
C
      DO 116 I=2,3
      DO 106 J=1,3
      P(J,I)=0.D0
      IF(A(J,I).NE.0.D0)P(J,I)=DATAN2(B(J,I),A(J,I))
      IF(P(J,I).LT.0.D0)P(J,I)=PI2+P(J,I)
      PDIF(J,I)=P(J,I)-PTEMP(J,I)
      IF(PTEMP(J,I).GT.P(J,I).AND.DABS(PTEMP(J,I)-P(J,I)).GT.PI)
     +                                           PDIF(J,I)=PDIF(J,I)+PI2
      IF(PTEMP(J,I).LT.P(J,I).AND.DABS(PTEMP(J,I)-P(J,I)).GT.PI)
     +                                           PDIF(J,I)=PDIF(J,I)-PI2
      PTOT(J,I)=PTOT(J,I)+PDIF(J,I)/PI2
      PTEMP(J,I)=P(J,I)
  106 CONTINUE
  116 CONTINUE
C     IF(S.LT.200.0)WRITE(53,1066)
C    +                S,NAME(ITY),IID,((PTOT(J1,J2),J1=1,3),J2=2,3)
 1066 FORMAT(' ',1X,F9.2,1X,A4,I4,6F12.6)
C
C
  104 CONTINUE
C
      WRITE(53,1107)
 1107 FORMAT('0',//,' ACCUMULATED SPIN PHASE ADVANCES:')
      DO 111 JJ=2,3
  111 WRITE(53,1106)JJ,(PTOT(I,JJ),I=1,3)
 1106 FORMAT('0',4X,'EIGENVECTOR',I2,4X,3F12.6)
C
C
C
C
C
C
C
      RETURN
C
C
C
C
C
 9988 WRITE(53,93)NDAMP
   93 FORMAT(' 6X6 & 8X8  EIGENTUNES MISS-MATCHED SO STOP',I4)
      STOP
 9999 WRITE(53,921)
  921 FORMAT(' ERROR IN EIGEN VALUE ROUTINE')
      STOP
C
C
C
C
C
C
C=====LIST THE FORMATS
C
  933 FORMAT(' ',' E0/.440652=',T14,F9.6,3X,F9.6,3X,F9.6,
     +                               //,' SPIN ROT. MATRIX AROUN',
     +   'D',T50,'EIGENVALUES:',T79,'SPIN BASIS FROM EIGENVECTORS:',/,
     +   T5,'THE 1-ST BEAM-LINE ELEMENT:',T52,'REAL    IMAG')
C
C
  932 FORMAT(3F10.5,T50,2F9.5,T71,'--->  (',T79,3F9.5,T107,')')
C
C
  931 FORMAT('0',' ORTHONORMAL SPIN BASE VECTORS (N,M,L):',//,4X,
     +' POS     NAME    IID ',T28,'NX',T37,'NY',T46,'NZ',T65,'MX',T74
     +,'MY',T83,'MZ',T102,'LX',T111,'LY',T120,'LZ',/)
C
  934 FORMAT(' ','ELEMENT,N AXIS TILT(MRAD) ',1X,A8,F10.4)
C
  936 FORMAT('1','   "N" VECTOR & "M" AND "L" PARTS OF "D" VECTOR',//,
     +'    POS   NAME   IID ',3X,'NX',6X,'NY',6X,'NZ',
     +                       9X,'DMX',7X,'DMY',7X,'DMZ',7X,'DDM',
     +                      11X,'DLX',7X,'DLY',7X,'DLZ',7X,'DDL'/)
C
  930 FORMAT(1X,F9.2,1X,A8,1X,I2,3F9.5,10X,3F9.5,10X,3F9.5,I8)
C
C
  937 FORMAT(' ',F6.1, 1X,A8,1X,I2, 3F8.4, 4X,4F10.4,4X,4F10.4)
  938 FORMAT(' ',E16.6,1X,A8,1X,2I6,3E16.6,4X,4E16.6,4X,42E16.6)
C
C
  943 FORMAT('0',' PURE SOKOLOV/TERNOV:  NUMERATOR    = ',E18.8,/,
     +          '                        DENOMINATOR  = ',E18.8,/,
     +          '                        POLARISATION = ',F10.4,/,
     +          '                        TIME(SECONDS)= ',F10.4)
C
C
  980 FORMAT(//,' R.M.S. TILT OF n0-AXIS VS. VERTICAL IN ARCS & MICRO',
     +' BETA REGIONS=>',F8.4,' MRADS.')
  983 FORMAT(//,
     +       ' R.M.S. TILT OF n0-AXIS VS. BEAM     IN MICRO BETA REGIONS
     +       =>',F8.4,' MRADS.')
  981 FORMAT(' R.M.S. TILT OF n0-AXIS VS. VERTICAL IN ARC QUADS
     +       =>',F8.4,' MRADS.')
  982 FORMAT(' R.M.S. TILT OF n0-AXIS VS. VERTICAL IN MICRO BETA REGIONS
     +       =>',F8.4,' MRADS.')
  984 FORMAT(' R.M.S. TILT OF DERBENEV-KONDRATENKO N-AXIS FROM N0
     +       =>',F8.4,' MRADS.')
C
C
  942 FORMAT(' TOT. POLN = ',F10.4,
     +       ' TAUD (DIFFUSION ) = ',D15.8,' SEC    TAUP/TAUD = ',D15.8)
  947 FORMAT(' POLN. 1   = ',F10.4,
     +       ' TAUD (DIFFUSION1) = ',D15.8,' SEC    TAUP/TAUD = ',D15.8)
  948 FORMAT(' POLN. 2   = ',F10.4,
     +       ' TAUD (DIFFUSION2) = ',D15.8,' SEC    TAUP/TAUD = ',D15.8)
  949 FORMAT(' POLN. 3   = ',F10.4,
     +       ' TAUD (DIFFUSION3) = ',D15.8,' SEC    TAUP/TAUD = ',D15.8)
  950 FORMAT(' POLN. 4   = ',F10.4,
     +       ' TAUD (DIFFUSION4) = ',D15.8,' SEC    TAUP/TAUD = ',D15.8)
C
C
  944 FORMAT('0',' DEGREE OF POLARISATION WITH DIFFUSION  = ',F7.3,
     +                                        ' ENERGY = ',F6.3,' GEV')
  946 FORMAT('0',' OR WITH NUMERATOR "D" TERM  ALSO       = ',F7.3,
     +                                        ' ENERGY = ',F6.3,' GEV'/)
  945 FORMAT('0',' OR WITH UNITARISED"D" TERM  ALSO       = ',F7.3,
     +                                        ' ENERGY = ',F6.3,' GEV'/)
C
C
C
C
C
C
      END
