C   21/11/95 511212042  MEMBER NAME  SODOM    (SEPT95.S) M  FORTRAN
C   07/03/95 510041423  MEMBER NAME  SODOM    (MAR95.S)  M  FORTRAN
      SUBROUTINE SODOM(IE0,E0,CIR,PTUN)
C
C
C
C
C   Routine to manage the SODOM algorithm of K.Yokoya.
C   -------------------------------------------------
C
C   The strategy is to use a loop structure similar to that in SPIN
C   and in the process use as many as possible of the facilities
C   already existing.
C   Several calculations already made in SPIN are repeated in SODOM
C   so as to make SODOM independent of SPIN.
C
C   Steps:
C
C   1. Get the SLIM spin dreibein. Set m and l so that m is perp. to the
C      longitudinal axis as in FIDO.
C
C   2. Get the (re)normalised 8x8 eigen vectors. Just need the 6x6 part.
C
C   3. Initially track one orbit around the ring using TSPIN and RSPIN.
C
C   4. Transform the 3x3 rotation matrix into the C.O. Dreibein.
C
C   5. Write the new 3x3 matrix in spinor form.
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
      IMPLICIT REAL*8(A-H,O-Z)
C
#include "cnlist.for"
#include "clatic.for"
#include "cloorb.for"
#include "csol.for"
#include "cemit.for"
C
      DIMENSION ZR3(3,3),ZI3(3,3),A(3,3),B(3,3),P(3,3),PTEMP(3,3)
      DIMENSION PDIF(3,3),PTOT(3,3)
      DIMENSION ROT(3,3),TM3A(3,3),ROTOFF(3,3),ROTDIP(3,3),ROTQUA(3,3)
      DIMENSION ROTINV(3,3),RINV(3,3)
      DIMENSION SBAS(3,3)
      DIMENSION TM3B(3,3),WR3(3),WI3(3),ZW(3,3),TM3C(3,3)
      DIMENSION TREV8(8,8),ZZ(8,8),ZV(8,8),WR8(8),WI8(8)
      DIMENSION TM8A(8,8),TM8B(8,8)
      DIMENSION SOL(8,8)
      DIMENSION AB(6)
      DIMENSION TN(8),PTUN(6)
      DIMENSION TRAJ(8),TRAJ1(8),TRAJAV(8)
      REAL*8 NU
      DATA ZV/64*0.D0/
      REAL*8 DNAME/'D'/,IP/'IP'/
      REAL*8                       K6/'K6  '/
      REAL*8                     EDGE/'EDGE'/
      LOGICAL IPRIN1,IPRIN2
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
      PI=3.1415926535897932D0
      PI2=2.D0*PI
      NU=E0/0.440652D0
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
C
      WRITE(6,929)
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
      IF(DANGH.LT.1.D-4.AND.DABS(ANGV).LT.1.D-8)GO TO 2311
      WRITE(6,2312)ANGH,ANGV
 2312 FORMAT(' ','TOTAL BENDING ANGLE CRAZY ',2F15.10)
      STOP
 2311 CONTINUE
C
C=====1)GET SPIN TUNE FROM TRACE OF MATRIX.
      STUNE=(ROT(1,1)+ROT(2,2)+ROT(3,3)-1.D0)*0.5D0
      STUNE=DARCOS(STUNE)/(2.D0*PI)
      RNU=DARCOS(DCOS(0.01D0)*DCOS(PI*NU))/PI-1.D0
C=====3)CHECK ORTHOGONALITY OF THE ROTATION MATRIX
      CALL ORTCHK(ROT)
C
C
C
C=====GET EIGENVECTORS & ORDER THEM=====================================
      CALL EIV3(ROT,ZR3,ZI3,WR3,WI3,IERRO)
      IF(IERRO .NE. 0) GO TO 9999
      NN=4
      DO 255 I=1,3
      TM3A(I,1)=WR3(I)
      TM3A(I,2)=WI3(I)
 255  IF(DABS(WI3(I)) .LT. 1.D-5)NN=I
      IF(NN .EQ. 4)STOP
      MM=MOD(NN,3)+1
C
C=====Arrange things so that independently of the 'm' eigenvector
C=====from EIV3, we always have an 'm' which is perpendicular to the
C=====longitudinal beam axis so that spin basis vector components are
C=====continuous with energy and don't depend on arbitrary phases
C=====arising in the EIV3 formalism. So, REDEFINE 'm'.
C=====The resulting 'm' must be normal to 'n'.
C
      ZR3(1,MM)= ZR3(2,NN)/DSQRT(ZR3(1,NN)**2+ZR3(2,NN)**2)
      ZR3(2,MM)=-ZR3(1,NN)/DSQRT(ZR3(1,NN)**2+ZR3(2,NN)**2)
      ZR3(3,MM)= 0.D0
C=====Or for an 'm' with no transverse component:
C     ZR3(1,MM)= 0.D0
C     ZR3(2,MM)=-ZR3(3,NN)/DSQRT(ZR3(2,NN)**2+ZR3(3,NN)**2)
C     ZR3(3,MM)= ZR3(2,NN)/DSQRT(ZR3(2,NN)**2+ZR3(3,NN)**2)
C
      AB1=DSQRT(ZR3(1,NN)**2+ZR3(2,NN)**2+ZR3(3,NN)**2)
      AB2=DSQRT(ZR3(1,MM)**2+ZR3(2,MM)**2+ZR3(3,MM)**2)
      DO 130 I=1,3
      ZW(I,1)=ZR3(I,NN)/AB1
  130 ZW(I,2)=ZR3(I,MM)/AB2
C=====Define l = n x m  -- guarantees rh coord system
      ZW(1,3)=ZW(2,1)*ZW(3,2)-ZW(3,1)*ZW(2,2)
      ZW(2,3)=ZW(3,1)*ZW(1,2)-ZW(1,1)*ZW(3,2)
      ZW(3,3)=ZW(1,1)*ZW(2,2)-ZW(2,1)*ZW(1,2)
C
C=====Must get spin tune, including tune-shift due to rotator
      CALL JAM333(ZR3,ROT,ZW)
      AB1=DSQRT(ZR3(1,2)**2+ZR3(2,2)**2+ZR3(3,2)**2)
      AB2=DSQRT(ZR3(1,3)**2+ZR3(2,3)**2+ZR3(3,3)**2)
      DO 140 I=1,3
      ZR3(I,2) = ZR3(I,2)/AB1
  140 ZR3(I,3) = ZR3(I,3)/AB2
      DOT1 = ZR3(1,2)*ZW(1,2) + ZR3(2,2)*ZW(2,2) + ZR3(3,2)*ZW(3,2)
      DOT2 = ZR3(1,2)*ZW(1,3) + ZR3(2,2)*ZW(2,3) + ZR3(3,2)*ZW(3,3)
      DOT3 = ZR3(1,3)*ZW(1,3) + ZR3(2,3)*ZW(2,3) + ZR3(3,3)*ZW(3,3)
      DOT4 = ZR3(1,3)*ZW(1,2) + ZR3(2,3)*ZW(2,2) + ZR3(3,3)*ZW(3,2)
      IF(DABS(DOT1-DOT3).LT.1.D-4.AND.DABS(DOT2+DOT4).LT.1.D-4)GOTO 150
      WRITE(6,145) DOT1,DOT2,DOT3,DOT4
  145 FORMAT(' ERROR IN ROTATION OF L & M',/1X,2F9.5,/1X,2F9.5)
      STOP
C=====Spin tune in all Mais & Ripken theory is defined by angle of
C=====rotation of the dreibein--not by the matrix eigenvalues.
  150 STUNE = DATAN2(-DOT2,DOT1)/2.D0/PI
      IF (STUNE.LT.0) STUNE=STUNE+1
      TTUNE=(ROT(1,1)+ROT(2,2)+ROT(3,3)-1.D0)*0.5D0
      TTUNE=DARCOS(TTUNE)/(2.D0*PI)
      WRITE(6,185) NU,STUNE,TTUNE
  185 FORMAT('0',' GAMMA*A     = ',F9.4,9X,'TUNE = ',2F9.4,/)
C
C
      WRITE(6,933)NU,STUNE
      DO 250 I=1,3
  250 WRITE(6,932)(ROT(I,J),J=1,3),WR3(I),WI3(I),(ZW(J,I),J=1,3)
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
      IF(IPRIN1)WRITE(6,931)
C
C
C
C
C
      S=0.
      ISOL=0
      CALL UNIT(8,SOL)
      CALL UNIT(8,TREV8)
      CALL UCOPY(ZW,ROT,18)
      CALL UCOPY(ZW,TM3C,18)
C=====LOOP AROUND THE LATTICE===========================================
      DO 257 II=1,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
C=====GET THE LENGTH:
      IF(IID.NE.5.AND.IID.NE.6.AND.IID.NE.7)S=S+YY(ITY)
      CALL UCOPY(TM3C,ZW,18)
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
      CALL JAM888(TREV8,TM8A,TREV8)
 2588 CONTINUE
C
C
C
C     IF((NAME(ITY).EQ.IP.OR.II.LT.1000.OR.II.GT.6930).AND.IPRIN1)
C     IF(                     IPRIN1 )
      IF((NAME(ITY).EQ.IP.AND.IPRIN1))
     +                               WRITE(6,930)S,NAME(ITY),IID,
     +                (ZW(J,1),J=1,3),(ZW(J,2),J=1,3),(ZW(J,3),J=1,3),II
C
C
C
C
C=====SUM THE N-AXIS DEVIATIONS OVER THE NON-VERTICAL BEND REGION:
C     THE ARCS,THE MICRO BETA REGION & THE SUM OF THESE.
C
  267 ZZ2=1.D0-ZW(2,1)*ZW(2,1)
      ZZT=1.D0-ZW(3,1)*ZW(3,1)
C     KTILT=0
C     IF((S.LT.12.0).OR.(S.GT.132.0.AND.S.LT.156.0).OR.(S.GT.276.0))
C    +KTILT=1
C     IF((S.GT.20.0.AND.S.LT.124.).OR.(S.GT.164.AND.S.LT.268))
C    +KTILT=2
C=====FOR HERA
      KTILT=0
      IF(S.LT.25.  )KTILT=1
      IF(S.GT.230.0.AND.S.LT.792.0)KTILT=2
      IF(KTILT.EQ.1)SUMNZM=SUMNZM+DSQRT(ZZ2)
      IF(KTILT.EQ.1)NNZM=NNZM+1
      IF(KTILT.EQ.1)SUMNTM=SUMNTM+DSQRT(ZZT)
      IF(KTILT.EQ.1)NNTM=NNTM+1
      IF(KTILT.EQ.2)SUMNZA=SUMNZA+DSQRT(ZZ2)
      IF(KTILT.EQ.2)NNZA=NNZA+1
      IF(KTILT.EQ.0)GO TO 257
      SUMNZ=SUMNZ+DSQRT(ZZ2)
      NNZ=NNZ+1
  257 CONTINUE
C
C
C=====WIND BACK THE SPIN BASIS==========================================
      CALL UNIT(8,TM8A)
      TM8A(7,7)=WR3(2)
      TM8A(7,8)=WI3(2)
      TM8A(8,7)=-WI3(2)
      TM8A(8,8)=WR3(2)
C=====For comparison with the 3x3 matrix in the Chao frame, kill
C     this wind back.
C     CALL JAM888(TREV8,TM8A,TREV8)
      WRITE(6,927)
  927 FORMAT(///,' 8X8 TRANSF.MATRIX FOR ONE REVOLUTION AROUND THE',
     +                                        ' 1-ST BEAM-LINE ELEM:',/)
      DO 273 I=1,8
  273 WRITE(6,925)(TREV8(I,J),J=1,8)
  925 FORMAT(8D16.8)
C
C
      IF(2.EQ.2)GO TO 30000
C
C
C=====GET 8X8 EIGENVECTORS FOR ONE REVOLUTION & ORDER THEM==============
      CALL EIV8(TREV8,ZZ,WR8,WI8,IERRO)
      IF(IERRO .NE. 0)GO TO 9999
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
      CALL UCOPY(ZZ,TM8A,128)
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
      WRITE(6,924)
  924 FORMAT(////,' EIGENVALUES AND EIGENVECTORS:',/,T8,'REAL',
     +                                          T18,'IMAG',T122,'TUNES')
      DO 281 I=1,8
      TUN=999999.
C=====GET TUNES: THIS DOES NOT DIFFERENTIATE BETWEEN +/- ANGLES
C=====WITHOUT USING THE SINE.
C     IF(DABS(WR8(I)) .LE. 1.D0)TUN=DARCOS(WR8(I))/(2.*PI)
      TUN=DATAN2(WI8(I),WR8(I))/(2.*PI)
      TN(I)=TUN
  281 WRITE(6,923)WR8(I),WI8(I),(ZZ(J,I),J=1,8),TUN
  923 FORMAT(2F12.5,2X,'--->  (',8F10.5,' )',F15.8)
      CALL NORM(ZZ,8,AB)
C
C
C
      CALL RESON(TN(1),TN(3),TN(5),TN(7))
C
C
C
C
C
      IF(MRDISP.EQ.0)GO TO 765
C=====IF USING THE DISPERSION VERSION GET THE UNCOUPLED TWISSES.
      RMUX=DARCOS(0.5D0*(TREV8(1,1)+TREV8(2,2)))
      IF(TREV8(1,2) .LT. 0.)RMUX=2.D0*PI-RMUX
      RMUY=DARCOS(0.5D0*(TREV8(3,3)+TREV8(4,4)))
      IF(TREV8(3,4).LT.0.D0)RMUY=2.D0*PI-RMUY
      BETAX=TREV8(1,2)/DSIN(RMUX)
      BETAY=TREV8(3,4)/DSIN(RMUY)
      ALPHAX=(TREV8(1,1)-TREV8(2,2))/(2.*DSIN(RMUX))
      ALPHAY=(TREV8(3,3)-TREV8(4,4))/(2.*DSIN(RMUY))
      WRITE(6,763)
  763 FORMAT('0','TWISSES FROM UNCOUPLED 8X8 MATRIX:',/
     +'     BETAX     BETAY    ALPHAX     ALPHAY')
      WRITE(6,762)BETAX,BETAY,ALPHAX,ALPHAY
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
      NTUN=0
      DO 25 JJ=1,5,2
      DO 24 J =1,5,2
C     DIF=DABS(TN(JJ)-PTUN(J))
      DIF=DABS(DABS(TN(JJ))-DABS(PTUN(J)))
      IF(DIF.LT.0.0010)AASPIN(JJ)=AAN(J)
      IF(DIF.LT.0.0010)NTUN=NTUN+1
      IF(DIF.LT.0.0010)WRITE(6,23)JJ,J
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
C
C
C
30000 CONTINUE
C
C
C
C     *****************************************************
C     *      Tracking a sample particle for one turn      *
C     *      to get its 3x3 spin rotation matrix.         *
C     *                                                   *
C     *      Do it two ways:                              *
C     *                                                   *
C     *      1) In machine coordinates using the D.O. AND *
C     *         the synchobeta terms in dipoles + C.F.    *
C     *                                                   *
C     *      2) In the Chao dreibein at each step using   *
C     *         just the synchrobeta terms at each step   *
C     *         with  a transform back into machine coords*
C     *         at the end.                               *
C     *                                                   *
C     *                                                   *
C     *                                                   *
C     *****************************************************
C
C
C=====Method 1.
C
C
C=====Store the one turn 8x8 matrix.
      CALL UCOPY(ZZ,TM8B,128)
C=====The starting Chao spin basis (dreibein) is now in ROT
C=====having been copied from ZW. Now store it in SBAS.
C=====So SBAS encodes the orientation of the Chao dreibein wrt machine
C=====coords.
      CALL UCOPY(ROT,SBAS,18)
      S=0.
C
C=====The off-D.O. rotation matrix is in ROTOFF and includes the
C=====D.O. pure dipole effect as well as the off D.O.synchrobeta terms.
      CALL UNIT(3,ROTOFF)
C=====Sample trajectory including linear spin components.
      TRAJ1(1)=0.00010D0    *0.0D0
      TRAJ1(2)=0.0000D0     *0.0D0
      TRAJ1(3)=0.000100D0   *1.0D0
      TRAJ1(4)=0.0000D0     *0.0D0
      TRAJ1(5)=0.0001D0     *0.0D0
      TRAJ1(6)=0.00001D0    *0.0D0
      TRAJ1(7)=0.0000D0
      TRAJ1(8)=0.0000D0
      WRITE(6,9511)TRAJ1
 9511 FORMAT('1','Method 1',/,' Input  TRAJ1 ',8F10.5)


      DO 300 II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY =XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
      IF(IID.NE.5.AND.IID.NE.6.AND.IID.NE.7)S=S+YY(ITY)
      CALL UCOPY(TRAJ1,TRAJ,16)
      CALL UCOPY(SBAS,TM3B,18)
C=====Standard TSPIN to rotate the dreibein: pure design orbit rotation.
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      CALL JAM333(SBAS,TM3A,TM3B)
C=====Average dreibein.
      CALL AVER(SBAS,TM3B,3,ZW)
C=====Transport the orbit (and linear spin).
      CALL MX88(IID,II,ITY,TMAT(1,1,ITY),XY,XX2,YYY,NU,ZW,TM3B,TM8A)
      CALL JAM881(TRAJ1,TM8A,TRAJ)
      CALL AVORB(TRAJ,TRAJ1,8,TRAJAV)
C=====Using the average orbit in the dipole, quadrupole or
C     combined function magnet, get the 3x3 rotation matrix.
      CALL TTSPIN(II,IID,ITY,XY,XX2,YYY,NU,TRAJAV,TM3A)
      CALL JAM333(ROTOFF,TM3A,ROTOFF)
C
C
C
C=====END OF LATTICE LOOP===============================================
  300 CONTINUE
C
C
      WRITE(6,951)TRAJ1
  951 FORMAT(' ','Method 1',/,' Output TRAJ1 ',8F10.5)
      WRITE(6,973)NU,STUNE
      WRITE(6,974)
  974 FORMAT(' ','Method 1: the matrix wrt machine coodinates.')
      DO 350 I=1,3
  350 WRITE(6,932)(ROTOFF(I,J),J=1,3),WR3(I),WI3(I),(ZW(J,I),J=1,3)


C
C
C=====Method 2.
C
C
C=====Store the one turn 8x8 matrix.
      CALL UCOPY(ZZ,TM8B,128)
C=====The starting Chao spin basis (dreibein) is now in ROT
C=====having been copied from ZW. Now store it in SBAS.
C=====So SBAS encodes the orientation of the Chao dreibein wrt machine
C=====coords.
      CALL UCOPY(ROT,SBAS,18)
      S=0.
C
C=====The rotation matrix in the Chao frame just due to the off-D.O.
C=====terms is in ROTQUA.
C=====The rotation matrix just due to D.O. dipole efffect is in ROTDIP.
      CALL UNIT(3,ROTDIP)
      CALL UNIT(3,ROTQUA)
C=====Sample trajectory including linear spin components.
      TRAJ1(1)=0.00010D0    *0.0D0
      TRAJ1(2)=0.0000D0     *0.0D0
      TRAJ1(3)=0.000100D0   *1.0D0
      TRAJ1(4)=0.0000D0     *0.0D0
      TRAJ1(5)=0.0001D0     *0.0D0
      TRAJ1(6)=0.00001D0    *0.0D0
      TRAJ1(7)=0.0000D0
      TRAJ1(8)=0.0000D0
      WRITE(6,9512)TRAJ1
 9512 FORMAT('1','Method 2',/,' Input  TRAJ1 ',8F10.5)
C
C
      DO 3000 II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY =XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
      IF(IID.NE.5.AND.IID.NE.6.AND.IID.NE.7)S=S+YY(ITY)
      CALL UCOPY(TRAJ1,TRAJ,16)
      CALL UCOPY(SBAS,TM3B,18)
C=====Standard TSPIN to rotate the dreibein: pure design orbit rotation.
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      CALL JAM333(SBAS,TM3A,TM3B)
C=====Accumulate  the total rotation of the machine coords.
      CALL JAM333(ROTDIP,TM3A,ROTDIP)
C=====Average dreibein.
      CALL AVER(SBAS,TM3B,3,ZW)
C=====Transport the orbit (and linear spin).
      CALL MX88(IID,II,ITY,TMAT(1,1,ITY),XY,XX2,YYY,NU,ZW,TM3B,TM8A)
      CALL JAM881(TRAJ1,TM8A,TRAJ)
      CALL AVORB(TRAJ,TRAJ1,8,TRAJAV)
C=====Using the average orbit in the quadrupole or combined function
C==== magnet, get the 3x3 rotation matrix wrt the Chao dreibein.
C=====Only include the synchrobeta terms.The D.O.part was done in TSPIN.
      CALL TQSPIN(II,IID,ITY,XY,XX2,YYY,NU,TRAJAV,ZW,TM3A)
      CALL JAM333(ROTQUA,TM3A,ROTQUA)
C
C
C
C=====END OF LATTICE LOOP===============================================
 3000 CONTINUE
C
      WRITE(6,9513)TRAJ1
 9513 FORMAT(' ','Method 2',/,' Output TRAJ1 ',8F10.5)
      WRITE(6,973)NU,STUNE
      WRITE(6,9741)
 9741 FORMAT(' ','Method 2: the matrix wrt Chao    coodinates.')
      DO 3501 I=1,3
 3501 WRITE(6,932)(ROTQUA(I,J),J=1,3),WR3(I),WI3(I),(ZW(J,I),J=1,3)
      RINV(1,1)=ROT(1,1)
      RINV(1,2)=ROT(2,1)
      RINV(1,3)=ROT(3,1)
      RINV(2,1)=ROT(1,2)
      RINV(2,2)=ROT(2,2)
      RINV(2,3)=ROT(3,2)
      RINV(3,1)=ROT(1,3)
      RINV(3,2)=ROT(2,3)
      RINV(3,3)=ROT(3,3)
C=====Rotate back from the Chao dreibein into the 'dipole' system.
      CALL JAM333(ROTQUA,ROTQUA,RINV)
      CALL JAM333(ROTQUA,ROT ,ROTQUA)
C=====Rotate back into the machine coords.
      CALL JAM333(ROTOFF,ROTDIP,ROTQUA)
C
      WRITE(6,9742)
 9742 FORMAT(' ','Method 2: the matrix wrt machine coodinates.')
      WRITE(6,973)NU,STUNE
      DO 3500 I=1,3
 3500 WRITE(6,932)(ROTOFF(I,J),J=1,3),WR3(I),WI3(I),(ZW(J,I),J=1,3)




C
      IF(1.EQ.1)RETURN
C=====GET 1/(RHO**3) INTEGRANDS INCLUDING C.O. DEVIATION EFFECTS
      IDIP=0
      DELX=(DX(II)+DX(II+1))/2.
      DELY=(DY(II)+DY(II+1))/2.
      XY1=0.D0
      GO TO(311,302,303,304,300,302,302,308,302,310,300,300,300,300,
     +                                                      302,302),IID
      GO TO 300
  311 IF(II.NE.2)GO TO 300
      XY1=0.D0
      GO TO 210
C
C
  302 XY1=DABS(XY)**3/YY(ITY)**2
      IDIP=1
C
C=====CHECK EFFECT OF CHANGING THE LENGTH OF SPECIAL DIPOLES.
C     IF(NAME(ITY).EQ.LON1.OR.NAME(ITY).EQ.LON2)XY1=XY1/9.
C     IF(NAME(ITY).EQ.LON3.OR.NAME(ITY).EQ.LON4)XY1=XY1/9.
C     IF(NAME(ITY).EQ.LON1.OR.NAME(ITY).EQ.LON2)XY1=XY1*100.
C     IF(NAME(ITY).EQ.LON3.OR.NAME(ITY).EQ.LON4)XY1=XY1/9.
C     IF(NAME(ITY).EQ.ZON1.OR.NAME(ITY).EQ.ZON2)XY1=XY1/9.
C     IF(NAME(ITY).EQ.BEN2.OR.NAME(ITY).EQ.BEN3)XY1=XY1*4.
C     IF(NAME(ITY).EQ.BEN5.OR.NAME(ITY).EQ.BEN6)XY1=XY1*4.
C     IF(NAME(ITY).EQ.K6  .OR.NAME(ITY).EQ.BEN6)XY1=XY1/100.
      GO TO 210
C
C
  303 CONTINUE
C     IF(NTWIST(ITY).EQ.3)DELY=DELY-TWIST(ITY)
C     IF(NTWIST(ITY).EQ.4)DELX=DELX-TWIST(ITY)
      XY1=((XY*DELX)**2+(XY*DELY)**2)**1.5/YY(ITY)**2
      GO TO 210
C
  304 CONTINUE
      XY1=((XY*DELX)**2+(XY*DELY)**2)**1.5/YY(ITY)**2
      GO TO 210
C
  308 XY1=((XY*(DELX**2-DELY**2)/2.)**2+(XY*DELX
     +                                        *DELY)**2)**1.5/YY(ITY)**200062800
      GO TO 210
C
C
  310 CALL SOLXYP(II,XY/YY(ITY),VX,VY)
      XY1=((XY*VX)**2+(XY*VY)**2)**1.5/YY(ITY)**2
C
C
  210 CONTINUE
C
      IF(MRDISP.EQ.0)CALL AVER(TREV8,TM8B,8,ZV)
C=====PUT IN HERE (TREV8)*('D'VECTOR) TO DEFINE ZV. NEED A SPECIAL AVER!
      IF(MRDISP.EQ.1)CALL DIAVER(TREV8,TM8B,ZV,II)
C     IF(MRDISP.EQ.1)CALL DIAVER(TREV8,TM8B,ZV,II,NAME(ITY))
      XY2=(-ZV(5,2)*ZV(7,1)+ZV(5,1)*ZV(7,2))/AB(1)
      XY3=(-ZV(5,4)*ZV(7,3)+ZV(5,3)*ZV(7,4))/AB(3)
      XY4=(-ZV(5,6)*ZV(7,5)+ZV(5,5)*ZV(7,6))/AB(5)
      XY5=(-ZV(5,2)*ZV(8,1)+ZV(5,1)*ZV(8,2))/AB(1)
      XY6=(-ZV(5,4)*ZV(8,3)+ZV(5,3)*ZV(8,4))/AB(3)
      XY7=(-ZV(5,6)*ZV(8,5)+ZV(5,5)*ZV(8,6))/AB(5)
C
C
C     IF(IID.EQ.2)WRITE(6,967)NAME(ITY),XY2,XY3,XY4,XY5,XY6,XY7
  967 FORMAT(' ',A8,7E15.4)
C
C
C     WRITE(6,927)
C     DO 2273 I=1,8
C2273 WRITE(6,925)(ZV(I,J),J=1,8)
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
      ENTILT=AASPIN(1)*(ZV(7,1)*ZV(7,1)+ZV(7,2)*ZV(7,2))/DABS(AB(1))
     +      +AASPIN(1)*(ZV(8,1)*ZV(8,1)+ZV(8,2)*ZV(8,2))/DABS(AB(1))
     +      +AASPIN(3)*(ZV(7,3)*ZV(7,3)+ZV(7,4)*ZV(7,4))/DABS(AB(3))
     +      +AASPIN(3)*(ZV(8,3)*ZV(8,3)+ZV(8,4)*ZV(8,4))/DABS(AB(3))
     +      +AASPIN(5)*(ZV(7,5)*ZV(7,5)+ZV(7,6)*ZV(7,6))/DABS(AB(5))
     +      +AASPIN(5)*(ZV(8,5)*ZV(8,5)+ZV(8,6)*ZV(8,6))/DABS(AB(5))
      TNTILT=ENTILT
      ENTILT=ENTILT*YY(ITY)
      NNAXIS=NNAXIS+1
      ENAXIS=ENAXIS+YY(ITY)
C     AXIST=AXIST+DSQRT(ENTILT)
      AXIST=AXIST+     (ENTILT)
      SNTILT=DSQRT(TNTILT)*1000.D0
C     IF(II.LT.4000)WRITE(6,934)NAME(ITY),SNTILT
C=====THESE IF STATEMENTS DON'T WORK IN THIS POSITION IF DRIFTS
C=====ARE SKIPPED.
      IF(NAME(ITY).EQ.DNAME)WRITE(6,934)NAME(ITY),SNTILT
C     IF(II.EQ.2        )WRITE(6,934)NAME(ITY),SNTILT
C                        WRITE(6,934)NAME(ITY),SNTILT
  277 CONTINUE
   27 CONTINUE
C
C
C=====DDM/DDL ARE 'D' PERTURBATION VECTORS ALONG M/L DIRECTIONS.
C=====AT THIS POINT "ZW" IS THE AVERAGE N,M,L BASIS--EXCEPT FOR SOLS.
      DDM=-XY2-XY3-XY4
      DDL=-XY5-XY6-XY7
      DMX=ZW(1,2)*DDM
      DMY=ZW(2,2)*DDM
      DMZ=ZW(3,2)*DDM
      DLX=ZW(1,3)*DDL
      DLY=ZW(2,3)*DDL
      DLZ=ZW(3,3)*DDL
      IF(NTWIST(ITY).GT.0.AND.IPRIN2.AND.IDIP.EQ.1)
     +                                      WRITE(6,937)S,NAME(ITY),IID,
     +                   (ZW(J,1),J=1,3),DMX,DMY,DMZ,DDM,DLX,DLY,DLZ,DDL
C
C
C
C=====GET SUM 1/(RHO**3) * EIGENVECTOR PRODUCTS: TOTAL & INDIVIDUAL.
      WD0=WD0+XY1*(DDM**2+DDL**2)
C
      WD1=WD1+XY1*(XY2**2+XY5**2)
C
      WD2=WD2+XY1*(XY3**2+XY6**2)
C
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
      SGN=-1.
      IF(XY.GT.0.)SGN=1.
C
C
C
C=====GET    N *(V * VDOT)/VDOT FACTOR---->WP2
C=====AND (N-D)*(V * VDOT)/VDOT FACTOR---->WP3
C=====IF XY2 IS SMALL IN LENSES(I.E. FIELD IS SMALL)---THEN IGNORE.
      GO TO(300,402,403,404,300,402,407,408,407,410,209,209,209,209,
     +                                                      402,407),IID
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
      WP0=WP0+XY1*(1.-2.*VN*VN/9.)
C
C     IF(NAME(ITY).EQ.EDGE)WRITE(6,155)IID,S,XY1,WP2,WP1I,WP0I,DELX,DELY
C 155 FORMAT(' ','EDGE ID,S,XY1,WP2,WP1I,WP0I,DELX,DELY',I4,F7.,6E13.3)
C
C
C
C
C=====END OF LATTICE LOOP===============================================
C 300 CONTINUE
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
      SYNCON=3.534D-19*(E0/0.000511)**5/CIR
      WP0=WP0*SYNCON
      WP1=WP1*SYNCON
      WND=WND*SYNCON
      TAUP=1.D0/WP0
      P0=92.38*WP1/WP0
      WRITE(6,943)WP1,WP0,P0,TAUP
C
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
      POL1=92.38*WP1/(WP0+WD0)
      POLM1=92.38*WP1/(WP0+WD1)
      POLM2=92.38*WP1/(WP0+WD2)
      POLM3=92.38*WP1/(WP0+WD3)
      POLM4=92.38*WP1/(WP0+WD4)
      WRITE(6,944)POL1,E0
      POL1UN=92.38*WP1/(WP0+WD0UN)
      WRITE(6,945)POL1UN,E0
      POL2=92.38*WND/(WP0+WD0)
      WRITE(6,946)POL2,E0
C
C
C
C
      XXX0=TAUP/TAUD0
      XXX1=TAUP/TAUD1
      XXX2=TAUP/TAUD2
      XXX3=TAUP/TAUD3
      XXX4=TAUP/TAUD4
      WRITE(6,942)POL1,TAUD0,XXX0
      WRITE(6,947)POLM1,TAUD1,XXX1
      WRITE(6,948)POLM2,TAUD2,XXX2
      WRITE(6,949)POLM3,TAUD3,XXX3
      WRITE(6,950)POLM4,TAUD4,XXX4
C
C
C
C
      SUMNZ =SUMNZ/NNZ
      SUMNZ =     (SUMNZ)*1000
      SUMNZM=SUMNZM/NNZM
      SUMNZM=     (SUMNZM)*1000
      SUMNTM=SUMNTM/NNTM
      SUMNTM=     (SUMNTM)*1000
      SUMNZA=SUMNZA/NNZA
      SUMNZA=     (SUMNZA)*1000
C
C
      WRITE (6,980)SUMNZ
      WRITE (6,983)SUMNTM
      WRITE (6,981)SUMNZA
      WRITE (6,982)SUMNZM
      WRITE (6,984)AXIST
C
      WRITE(20)WP1,WP0,WD0,WD1,WD2,WD3,WD4,P0,TAUP,
     +        POL1,E0, POL2,E0, POL1,TAUD0,XXX0,
     +                          POLM1,TAUD1,XXX1,
     +                          POLM2,TAUD2,XXX2,
     +                          POLM3,TAUD3,XXX3,
     +                          POLM4,TAUD4,XXX4,
     +                          SUMNZ, SUMNTM, SUMNZA,
     +                          SUMNZM, AXIST
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
C     WRITE(6,1107)
      CALL VZERO(PTOT,18)
      DO 104 II=1,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
      IF(IID.NE.5.AND.IID.NE.6.AND.IID.NE.7)S=S+YY(ITY)
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
     +                                           PDIF(J,I)=PDIF(J,I)+PI200099500
      IF(PTEMP(J,I).LT.P(J,I).AND.DABS(PTEMP(J,I)-P(J,I)).GT.PI)
     +                                           PDIF(J,I)=PDIF(J,I)-PI200099700
      PTOT(J,I)=PTOT(J,I)+PDIF(J,I)/PI2
      PTEMP(J,I)=P(J,I)
  106 CONTINUE
  116 CONTINUE
C     IF(S.LT.200.0)WRITE(6,1066)
C    +                S,NAME(ITY),IID,((PTOT(J1,J2),J1=1,3),J2=2,3)
 1066 FORMAT(' ',1X,F9.2,1X,A4,I4,6F12.6)
C
C
  104 CONTINUE
C
      WRITE(6,1107)
 1107 FORMAT('0',//,' ACCUMULATED SPIN PHASE ADVANCES:')
      DO 111 JJ=2,3
  111 WRITE(6,1106)JJ,(PTOT(I,JJ),I=1,3)
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
 9988 WRITE(6,93)NDAMP
   93 FORMAT(' 6X6 & 8X8  EIGENTUNES MISS-MATCHED SO STOP',I4)
      STOP
 9999 WRITE(6,921)
  921 FORMAT(' ERROR IN EIGEN VALUE ROUTINE')
      STOP
C
C
C
C
C
C
C=====LIST THE FORMATS
  929 FORMAT('1','Entering Subroutine SODOM...')
C
  933 FORMAT(' ',' E0/.440652=',T14,F9.6,3X,F9.6,12X,
     +                               //,' SPIN ROT. MATRIX AROUN',
     +   'D',T50,'EIGENVALUES:',T79,'SPIN BASIS FROM EIGENVECTORS:',/,
     +   T5,'THE 1-ST BEAM-LINE ELEMENT:',T52,'REAL    IMAG')
C
C
C
  973 FORMAT(' ',' E0/.440652=',T14,F9.6,3X,F9.6,12X,
     +                               //,' Off C.O.  matrix around the',
     +       T50,'Eigenvalues:',T79,'Spin basis from eigenvectors:',/,
     +   T2,'1-st beam line element:',T52,'Real    Imag')
C
C
  932 FORMAT(3F12.8,T50,2F9.5,T71,'--->  (',T79,3F9.5,T107,')')
C
C
  931 FORMAT('0',' ORTHONORMAL SPIN BASE VECTORS (N,M,L):',//,4X,
     +' POS     NAME    IID ',T25,'NX',T34,'NY',T43,'NZ',T62,'MX',T71
     +,'MY',T80,'MZ',T99,'LX',T108,'LY',T117,'LZ',/)
C
  934 FORMAT(' ','ELEMENT,N AXIS TILT(MRAD) ',1X,A8,F10.4)
C
  936 FORMAT('1','   "N" VECTOR & "M" AND "L" PARTS OF "D" VECTOR',//,
     +'    POS   NAME IID ',3X,' NX',6X,' NY',6X,' NZ',
     +                      10X,'DMX',7X,'DMY',7X,'DMZ',7X,'DDM',
     +                      11X,'DLX',7X,'DLY',7X,'DLZ',7X,'DDL'/)
C
  930 FORMAT(1X,F9.2,1X,A8,1X,I2,3F9.5,10X,3F9.5,10X,3F9.5,I8)
C
C
  937 FORMAT(' ',F8.1,1X,A8,1X,I2,3F9.4,4X,4E10.2,4X,4E10.2)
C
C
  943 FORMAT('0',' PURE SOKOLOV/TERNOV:  NUMERATOR    = ',E18.8,/,
     +          '                        DENOMINATOR  = ',E18.8,/,
     +          '                        POLARISATION = ',F10.4,/,
     +          '                        TIME(SECONDS)= ',F10.4)
C
C
  980 FORMAT(//,' R.M.S. TILT OF N-AXIS VS. VERTICAL IN ARCS & MICRO',
     +' BETA REGIONS=>',F8.4,' MRADS.')
  983 FORMAT(' R.M.S. TILT OF N-AXIS VS. BEAM     IN MICRO BETA REGIONS
     +      =>',F8.4,' MRADS.')
  981 FORMAT(' R.M.S. TILT OF N-AXIS VS. VERTICAL IN ARCS
     +      =>',F8.4,' MRADS.')
  982 FORMAT(' R.M.S. TILT OF N-AXIS VS. VERTICAL IN MICRO BETA REGIONS
     +      =>',F8.4,' MRADS.')
  984 FORMAT(' R.M.S. TILT OF DERBENEV-KONDRATENKO N-AXIS FROM N0
     +      =>',F8.4,' MRADS.')
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
      END
