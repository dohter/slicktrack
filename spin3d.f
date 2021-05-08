C   07/03/95 509071949  MEMBER NAME  SPIN     (SEPT95.S) M  FORTRAN
C   22/07/81 503061853  MEMBER NAME  SPIN     (MAR95.S)  M  FORTRAN
      SUBROUTINE SCRURITA(IE0,E0,CIR,TAUY,BETAY0,PTUN)
C
C
C
C
C   ROUTINE TO HANDLE SPIN-ORBIT TRACKING AND ESTIMATE THE RATE OF DEPOLARISATION.
C   ------------------------------------------------------------------------------
C   
C  In F90 with -fast, using the eRHIC lattice on solar10 
C  10000 turns take 140 secs for 100 particles if they all get transported en bloc with
C  MATMUL
C  44 energy steps and 10000 turns then needs 17 hours, i.e. 1 day. 
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
C
      PARAMETER (NPART = 100)
      DIMENSION ZR3(3,3),ZI3(3,3),A(3,3),B(3,3),P(3,3),PTEMP(3,3)
      DIMENSION PDIF(3,3),PTOT(3,3)
      DIMENSION ROT(3,3),TM3A(3,3)
      DIMENSION TM3B(3,3),WR3(3),WI3(3),ZW(3,3),TM3C(3,3)
      DIMENSION TRIN3(3,3),RR3(3),RI3(3),VR3(3,3),VI3(3,3)
      DIMENSION INTGE3(3),WW3(3,3)
      DIMENSION TREV8(8,8),ZZ(8,8),ZV(8,8),WR8(8),WI8(8)
      DIMENSION TREV08(8,8),RESSTR(2,6)                  ! For resonance strengths.
      DIMENSION TRIN8(8,8),RR8(8),RI8(8),VR8(8,8),VI8(8,8)
      DIMENSION INTGE8(8),WW8(8,8)
      DIMENSION TM8A(8,8),TM8B(8,8)
      DIMENSION SOL(8,8)
      DIMENSION ORBVEC(8,NPART),ORBVET(NPART,8)
      DIMENSION AB(6)
      DIMENSION TN(8),PTUN(6)
      REAL*8 NU
      DATA ZV/64*0.D0/
      DATA ORBVET/NPART*1.0D-4,
     +           NPART*1.0D-4,
     +           NPART*1.0D-4,
     +           NPART*1.0D-4,
     +           NPART*1.0D-2,
     +           NPART*1.0D-3,
     +           NPART*0.D0,
     +           NPART*0.D0/

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
      WRITE(53,929)
  929 FORMAT('1','Entering Subroutine SPIN3D...')
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
      CALL F02AGF(TRIN3,3,3,RR3,RI3,VR3,3,VI3,3,INTGE3,IFAIL)
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

 243  CONTINUE
      MM=MOD(NN,3)+1
      LL=MOD(MM,3)+1
      AZ1=DSQRT(ZR3(1,NN)**2+ZR3(2,NN)**2+ZR3(3,NN)**2)
      AZ2=DSQRT(ZR3(1,MM)**2+ZR3(2,MM)**2+ZR3(3,MM)**2)
      DO 256 I=1,3
      ZW(I,1)=ZR3(I,NN)/AZ1
 256  ZW(I,2)=ZR3(I,MM)/AZ2
C     Force the m vector to be vertical.
C      ZW(1,2)=0.D0
C      ZW(2,2)=1.D0
C      ZW(3,2)=0.D0
      ZW(1,3)=ZW(2,1)*ZW(3,2)-ZW(3,1)*ZW(2,2)
      ZW(2,3)=ZW(3,1)*ZW(1,2)-ZW(1,1)*ZW(3,2)
      ZW(3,3)=ZW(1,1)*ZW(2,2)-ZW(2,1)*ZW(1,2)
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
      ORBVEC = TRANSPOSE(ORBVET)
      
      DO 2258 MO=1,NPART
      ORBVEC(6,MO) =  ORBVEC(6,MO) - 1.D-6 * MO
 2258 ORBVEC(1,MO) =  ORBVEC(1,MO) + 1.D-6 * MO  

C======Loop over many turns

      DO 2257 MT= 1,1000

      




C=====LOOP AROUND THE LATTICE===========================================
      DO 257 II=1,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
     
      IF(IID.EQ.1)THEN
      ORBVEC(1,:) = ORBVEC(1,:) + XY * ORBVEC(2,:)
      ORBVEC(3,:) = ORBVEC(3,:) + XY * ORBVEC(4,:)
      S=S+YY(ITY)
      ELSE
C=====GET THE LENGTH and bend angle:
      IF(IID.EQ.2.OR.IID.EQ.15)ANGH=ANGH+XY
      IF(IID.NE.5.AND.IID.NE.6.AND.IID.NE.7)S=S+YY(ITY)
C      CALL UCOPY(TM3C,ZW,18)
C      IF(IID.EQ.10)GO TO 1001
C      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
C      CALL JAM333(TM3C,TM3A,ZW)
C      CALL AVER(ZW,TM3C,3,TM3B)
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
C      CALL JAM881(ORBVEC,TM8A,ORBVEC)
      ORBVEC = MATMUL(TM8A,ORBVEC)
      ENDIF



 2588 CONTINUE







  257 CONTINUE 
 2257 CONTINUE


C      WRITE(53,'(A,A,8(1X,F12.5),/)')
C     +                   ' ','Output vectors', ORBVEC




      RETURN


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
  938 FORMAT(' ',E16.6,1X,A8,1X,2I6,3E16.6,4X,4E16.6,4X,40E16.6)
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
  983 FORMAT(' R.M.S. TILT OF n0-AXIS VS. BEAM     IN MICRO BETA REGIONS
     +      =>',F8.4,' MRADS.')
  981 FORMAT(' R.M.S. TILT OF n0-AXIS VS. VERTICAL IN ARC QUADS
     +      =>',F8.4,' MRADS.')
  982 FORMAT(' R.M.S. TILT OF n0-AXIS VS. VERTICAL IN MICRO BETA REGIONS
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
C
C

      END
