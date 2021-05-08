C   22/07/81 507301043  MEMBER NAME  SSECT    (S)           FORTRAN
      SUBROUTINE SSECT(IE0,E0)
C
C
C
C
C   ROUTINE TO HANDLE 8X8 SPIN MATRIX ANALYSIS OF A SHORT SECTION OF
C   ----------------------------------------------------------------
C
C                            RING.
C                            -----
C   NEW GO FASTER VERSION---CUTTING OUT TIME WASTING DRIFT SPACE
C   MANIPULATIONS  14/4/83
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE "cnlist.for"
      INCLUDE "clatic.for"
      INCLUDE "cloorb.for"
      INCLUDE "csol.for"
C
      DIMENSION ROT(3,3),TM3A(3,3)
      DIMENSION TM3B(3,3),ZW(3,3),TM3C(3,3)
      DIMENSION TREV8(8,8),ZZ(8,8),WR8(8),WI8(8)
      DIMENSION TM8A(8,8),TM8B(8,8)
      REAL*8 NU
CDB      REAL*8 DNAME'D'/,IP/'IP'/
      DIMENSION SOL(8,8)
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
      NU=E0/0.440652D0
C
C
      WRITE(53,103)
      WRITE(53,103)
  103 FORMAT(/,'  ')


      WRITE(53,929)
  929 FORMAT('1','Entering subroutine SSECT to get the 8x8 matrix')
C
C
C
C    **************************************************
C    * ORTHONORMAL SPIN BASE VECTORS AND 8X8 MATRICES *
C    **************************************************
C
C=====PUT IN A SPECIAL SPIN BASIS FOR TESTING
C      CALL VZERO(ZW,18)
      ZW = 0.D0

C     CALL RSPIN(-1.65D0,0.1D0,0.3D0,T1)
C     CALL JAM333(ZW,T1,ZW)
C      ZW(1,1)= 0.D0
C      ZW(2,1)=-1.D0
C      ZW(3,1)= 0.D0

C
C      Force n_0 to be vertical.
C      ZW(1,1)=0.D0
C      ZW(2,1)=1.D0
C      ZW(3,1)=0.D0

C      Force n_0 to be longitudinal
      ZW(1,1)= 0.D0
      ZW(2,1)= 0.D0
      ZW(3,1)= 1.D0
C      Force m_0 to be radial.
      ZW(1,2)= 1.D0
      ZW(2,2)= 0.D0
      ZW(3,2)= 0.D0



      ZW(1,3)=ZW(2,1)*ZW(3,2)-ZW(3,1)*ZW(2,2)
      ZW(2,3)=ZW(3,1)*ZW(1,2)-ZW(1,1)*ZW(3,2)
      ZW(3,3)=ZW(1,1)*ZW(2,2)-ZW(2,1)*ZW(1,2)

      WRITE(53,931)
      WRITE(53,930)S,NAME(ITYPE(1)),ID(ITYPE(1)),
     +                   (ZW(J,1),J=1,3),(ZW(J,2),J=1,3),(ZW(J,3),J=1,3)
C
C
C
C
      S=0.
      ISOL=0
      CALL UNIT(8,TREV8)
      CALL UNIT(8,SOL)
C      CALL UCOPY(ZW,ROT,18)
      ROT = ZW
C      CALL UCOPY(ZW,TM3C,18)
      TM3C = ZW
C=====LOOP OVER THE LATTICE SECTION===================================
      DO 257 II=1,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      YYY=YY(ITY)
C=====GET THE LENGTH:
      IF(IID.NE.5)S=S+YY(ITY)
C      CALL UCOPY(TM3C,ZW,18)
      ZW = TM3C
      IF(IID.EQ.10)GO TO 1001
      CALL TSPIN(II,IID,ITY,XY,XX2,YYY,NU,TM3A)
      CALL JAM333(TM3C,TM3A,ZW)
      CALL AVER(ZW,TM3C,3,TM3B)
      CALL MX88(IID,II,ITY,TMAT(1,1,ITY),XY,XX2,YYY,NU,TM3B,ZW,TM8A)
      GOTO 258
 1001 CALL SOL8AN(II,XY,YY(ITY),NU,ZW,TM3C,TM8A,NSOL(ITY))
  258 CONTINUE
C
      CALL JAM888(TREV8,TM8A,TREV8)
 2588 CONTINUE
C
C
C
      IF(NAME(ITY).EQ.'IP')WRITE(53,930)S,NAME(ITY),IID,
     +                (ZW(J,1),J=1,3),(ZW(J,2),J=1,3),(ZW(J,3),J=1,3),II
C
C
C
C=====END OF LATTICE LOOP
  257 CONTINUE
C
C
      WRITE(53,927)
  927   FORMAT(///,' 8X8 TRANSF.MATRIX FOR ONE REVOLUTION AROUND THE'
     >    ,' 1-ST BEAM-LINE ELEM:',/)
      DO 273 I=1,8
  273 WRITE(53,925)(TREV8(I,J),J=1,8)
  925 FORMAT(8D16.8)
C
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
C=====LIST THE FORMATS
  933 FORMAT('1E0/.440652=',T14,F9.6,3X,F9.6,3X,F9.6,
     +                               //,' SPIN ROT. MATRIX AROUN',
     +   'D',T50,'EIGENVALUES:',T79,'SPIN BASIS FROM EIGENVECTORS:',/,
     +   T5,'THE 1-ST BEAM-LINE ELEMENT:',T52,'REAL    IMAG')
C
C
  932 FORMAT(3F10.5,T50,2F9.5,T71,'--->  (',T79,3F9.5,T107,')')
C
C
C
C
  931 FORMAT(' ORTHONORMAL SPIN BASE VECTORS (N,M,L):',//,
     +' POS     NAME IID ',2X,'NX',9X,'NY',9X,'NZ',12X,'MX',9X,'MY',
     +9X,'MZ',12X,'LX',9X,'LY',9X,'LZ',/)
C
C
  930 FORMAT(1X,F9.2,1X,A4,1X,I2,3F11.7,3X,3F11.7,3X,3F11.7,I6)

      END
