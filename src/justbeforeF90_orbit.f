C   06/10/83 702051841  MEMBER NAME  ORBIT    (S)           FORTRAN
      SUBROUTINE ORBIT(E0,IE0,U0,CRAD,FREQ,KKICK)
C
C
C================CALCULATE THE CLOSED-ORBIT 6-VECTOR====================
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
      DIMENSION TREV(7,7),TEMP(7,7),TM7(7,7),VC6A(6),VC6B(6),TM6(6,6)


C
C
C
      INCLUDE  "cdesy2.for"
      INCLUDE  "cloorb.for"
      INCLUDE  "cnlist.for"
      INCLUDE  "clatic.for"
      INCLUDE  "csol.for"
C
C
      IPU=0
      NQUADS = 0
C
C
C=====ITERATE TO FIND CLOSED ORBIT
      DO 23 ITORB=1,ITER
      IF(IE0.EQ.1.AND.IPTCO.EQ.1.AND.ITORB.EQ.ITER)
     +                                        WRITE(53,97)KKICK,FREQ
      IT=ITYPE(1)
      CALL UCOPY(TMAT(1,1,IT),TREV,98)
C
C=====CALCULATE 7X7 REVOLUTION MATRIX ON CLOSED ORBIT: INCLUDE RAD'N
C=====EFFECTS ONLY ON SECOND ITERATION.
      S=0.
      DO 26 II=2,NELEM
      ITY=ITYPE(II)
      IF(ID(ITY).EQ.1.AND.XX(ITY).LT.0.00000000001D0)GO TO 26 !Skip zero length drifts
      NTW=NTWIST(ITY)                                         !Sept 2003. 
      IF(ID(ITY).NE.5.AND.ID(ITY).NE.6.AND.ID(ITY).NE.7)S=S+YY(ITY)
      CALL UCOPY(TMAT(1,1,ITY),TM7,98)
      IID=ID(ITY)
C=====JUNE 97 as in:
C=====Oct 95:finally fix the kickers.Take code from SMILE3.JULY92.V2.
C=====For kickers, turn off curvature terms when getting the C.O.
       IF(IID.NE.6.AND.IID.NE.7) GO TO 27
           TM7(5,1)=0.D0
           TM7(5,3)=0.D0
           TM7(2,6)=0.D0
           TM7(4,6)=0.D0
   27 CONTINUE
      IF(ITORB.EQ.1)GO TO 71
C     IF(IRAD.EQ.0)GO TO 71
      DELX=(DX(II)+DX(II+1))/2.
      DELY=(DY(II)+DY(II+1))/2.
      GO TO (71,71,73,74,71,71,71,78,71,80,71,71,71,71,71,71),IID
   73 CONTINUE
C     IF(NTWIST(ITY).EQ.3)DELY=DY(II)-TWIST(ITY)
C     IF(NTWIST(ITY).EQ.4)DELX=DX(II)-TWIST(ITY)
      IF(IRAD.NE.0)TM7(6,7)=-CRAD*XX(ITY)**2*(DELX**2+DELY**2)/YY(ITY)
      GO TO 71
   74 IF(IRAD.NE.0)TM7(6,7)=-CRAD*XX(ITY)**2*(DELX**2+DELY**2)/YY(ITY)
      GO TO 71
   78 CALL SEXEQ(IRAD,CRAD,II,ITY,TM7)
C=====AS IN DECEMBER 1995: AS IN EMITNC TURN OFF RADIATION IN SOLENOIDS.
   80 CONTINUE
      GO TO 71
C  80 CALL SOLXYP(II,XX(ITY)/YY(ITY),DXPA,DYPA)
C     IF(IRAD.NE.0)
C    +           TM7(6,7)=-CRAD*XX(ITY)**2*(DXPA*DXPA+DYPA*DYPA)/YY(ITY)
   71 CONTINUE
C     IF(ITORB.EQ.2.AND.IID.EQ.8)WRITE(53,101)ITORB,II,ITY,IID
C     IF(ITORB.EQ.2.AND.IID.EQ.8)
C    +                          WRITE(53,102)((TM7(IL,IJ),IJ=1,7),IL=1,7)
C 101 FORMAT(' ','ITORB,II,ITY,ID', 4I7)
C 102 FORMAT(' ',7E16.6)
      CALL JAM777(TREV,TM7,TREV)
C
C
   26 CONTINUE
C
C
C
C
C
      DO 12 I=1,6
      DO 32 J=1,6
      TM6(I,J)=TREV(I,J)
   32 IF(I.EQ.J)TM6(I,J)=TREV(I,J)-1.D0
   12 VC6B(I)=-TREV(I,7)*1000.D0
      CALL SIMQ(TM6,VC6B,6,IERRO)
      IF(IERRO .NE. 0) STOP
      DX(1) =VC6B(1)/1000.D0
      DY(1) =VC6B(3)/1000.D0
      DXP(1)=VC6B(2)/1000.D0
      DYP(1)=VC6B(4)/1000.D0
      DL(1) =VC6B(5)/1000.D0
      DEL(1)=VC6B(6)/1000.D0
      DEN=DEL(1) *E0*1000.D0
      DXMAX= DABS(VC6B(1))
      DYMAX= DABS(VC6B(3))
      DSMAX= DABS(VC6B(5))
      DENMAX=DABS(DEN)
      DXRMS= VC6B(1)**2
      DYRMS= VC6B(3)**2
      DSRMS= VC6B(5)**2
      DENRMS=DEN**2
      DYPRMS=VC6B(4)**2
      IBEND = 0
      DYBEND= 0
      DYABS = 0
C
C
C
C
C=====CALCULATE NEW CLOSED ORBIT AT ENTRANCE TO ELEMENT II.
C=====DX(II-1) IS DX AT ENTRANCE TO ELEMENT(II-1) ETC.
      S=0.
      CALL VZERO(VC6A,12)
      DO 28 II=2,NELEM
      JTY=ITYPE(II-1)
C     IF(ID(JTY).EQ.1.AND.XX(JTY).LT.0.00000000001D0)GO TO 52 !Skip zero length drift 
      IF(ID(JTY).NE.5.AND.ID(JTY).NE.6.AND.ID(JTY).NE.7)S=S+YY(JTY)
      MTW=NTWIST(JTY)
      CALL UCOPY(TMAT(1,1,JTY),TEMP,98)
      IID=ID(JTY)
      IF(ITORB.EQ.1) GO TO 51
C     IF(IRAD.EQ.0)GO TO 51
      DELX=(DX(II-1)+DX(II))/2.
      DELY=(DY(II-1)+DY(II))/2.
      GO TO (51,51,53,54,51,51,51,58,51,60,51,51,51,51,51,51),IID
   53 CONTINUE
C     IF(NTWIST(JTY).EQ.3)DELY=DY(II-1)-TWIST(JTY)
C     IF(NTWIST(JTY).EQ.4)DELX=DX(II-1)-TWIST(JTY)
      IF(IRAD.NE.0)TEMP(6,7)=-CRAD*XX(JTY)**2*(DELX**2+DELY**2)/YY(JTY)
      GO TO 51
   54 IF(IRAD.NE.0)TEMP(6,7)=-CRAD*XX(JTY)**2*(DELX**2+DELY**2)/YY(JTY)
      GO TO 51
   58 CALL SEXEQ(IRAD,CRAD,II-1,JTY,TEMP)
      GO TO 51
C=====AS IN DECEMBER 1995: AS IN EMITNC TURN OFF RADIATION IN SOLENOIDS.
   60 CONTINUE
C  60 CALL SOLXYP(II-1,XX(JTY)/YY(JTY),DXPA,DYPA)
C     IF(IRAD.NE.0)
C    +          TEMP(6,7)=-CRAD*XX(JTY)**2*(DXPA*DXPA+DYPA*DYPA)/YY(JTY)
C
C
C

   51 CONTINUE 
      DO 210 I=1,6
      VC6A(I)=TEMP(I,7)*1000.D0
      DO 210 J=1,6
  210 VC6A(I)=VC6A(I)+TEMP(I,J)*VC6B(J)

   52 CONTINUE 
      DO 211 I=1,6
  211 VC6B(I)=VC6A(I)
      DX(II) =VC6B(1)/1000.D0
      DXP(II)=VC6B(2)/1000.D0
      DY(II) =VC6B(3)/1000.D0
      DYP(II)=VC6B(4)/1000.D0
      DL(II) =VC6B(5)/1000.D0
      DEL(II)=VC6B(6)/1000.D0
      DEN= DEL(II)*E0*1000.D0
      IF(DABS(VC6B(1)).GT.DXMAX)SXMAX=S
      IF(DABS(VC6B(1)).GT.DXMAX)DXMAX=DABS(VC6B(1))
      IF(DABS(VC6B(3)).GT.DYMAX)SYMAX=S
      IF(DABS(VC6B(3)).GT.DYMAX)DYMAX=DABS(VC6B(3))
      IF(DABS(VC6B(5)).GT.DSMAX)SSMAX=S
      IF(DABS(VC6B(5)).GT.DSMAX)DSMAX=DABS(VC6B(5))
      IF(DABS(DEL(II)*E0*1000.D0).GT.DENMAX)
     +                          DENMAX=DABS(DEL(II)*E0*1000.D0)
      IF(ITORB.NE.ITER) GO TO 28
      IF(IID.EQ.5)U0=U0+XX(JTY)*E0*VC6B(5)
      IF(IID.EQ.3.AND.NAME(JTY).NE.'EDGE')THEN
      DXRMS=DXRMS+VC6B(1)**2
      DYRMS=DYRMS+VC6B(3)**2
      DSRMS=DSRMS+VC6B(5)**2
      DENRMS=DENRMS+DEN**2
      DYPRMS=DYPRMS+VC6B(4)**2
      DYABS=DYABS+DABS(VC6B(3))
      NQUADS = NQUADS + 1
      ENDIF
      IF(IID.EQ.2.OR.IID.EQ.15)IBEND=IBEND+1
      IF(IID.EQ.2.OR.IID.EQ.15)DYBEND=DYBEND+VC6B(3)**2
C     IF(IE0.GT.-1.AND.ICLBP.NE.0.AND.IID.NE.1.AND.S.LT.250.)
C    +                    WRITE(53,98)II,S,NAME(JTY),VC6B,DEN,NAME(JTY)
C     IF(IE0.GT.-1.AND.ICLBP.NE.0.AND.ITORB.EQ.ITER.AND.(NAME(JTY).EQ.
C    +'IP'.OR.IID.EQ.5))WRITE(53,98)II,S,NAME(JTY),VC6B,DEN,NAME(JTY)
C     IF(IE0.GT.-1.AND.ICLBP.NE.0.AND.ITORB.EQ.ITER.AND.(NAME(JTY).EQ.
C    +'IP'            ))WRITE(53,98)II,S,NAME(JTY),VC6B,DEN,NAME(JTY)
C
      IF(IE0.EQ. 1.AND.ICLBP.NE.0.AND.ITORB.EQ.ITER.AND.(NAME(JTY).EQ.
     +'IP'.OR.IID.EQ.7.OR.IID.EQ.6.OR.NAME(JTY).EQ.'PU').AND.NAME(JTY)
     +.NE.'VKICK')
     +                    WRITE(53,98)II,S,NAME(JTY),VC6B,DEN,NAME(JTY)
C
C     IF(IE0.GT.-1.AND.ICLBP.NE.0.AND.MOD(II,10).EQ.0)
C    +           WRITE(53,98)ITORB,S,NAME(JTY),VC6B,DEN,NAME(JTY)
   98 FORMAT(1X,I4,1X,F10.4,1X,A8,2X,7F10.4,2X,A8)
C     IF(IE0.EQ.1.AND.ICLBP.NE.0 )WRITE(11)II,S,NAME(JTY),VC6B
C
C
C=====STORE ORBIT AT THE PICKUPS: USE MILLIMETERS.
      IF(ITORB.NE.ITER.OR.NAME(JTY).NE.'PU')GO TO 28
      IF(NAME(JTY).EQ.'PU')IPU=IPU+1
      LAG(2,IPU)=VC6B(1)
      LAG(1,IPU)=VC6B(3)
C
C
   28 CONTINUE
C
C
C     IF(DXMAX.LE.1.D-8.AND.DYMAX.LE.1.D-8)GO TO 63
C=====END CLOSED ORBIT ITERATION
   23 CONTINUE
C
C
C
C
C
      DXRMS=DSQRT(DXRMS/NQUADS)
      DYRMS=DSQRT(DYRMS/NQUADS)
      DSRMS=DSQRT(DSRMS/NQUADS)
      DENRMS=DSQRT(DENRMS/NQUADS)
      DYPRMS=DSQRT(DYPRMS/NQUADS)
      DYABS=DYABS/NQUADS
   63 CONTINUE
      WRITE(53,91)DXMAX,SXMAX,DYMAX,SYMAX,DSMAX,SSMAX,DENMAX,
     +           DXRMS,DYRMS,DSRMS,DENRMS,
     +DYPRMS,DYABS,DYBEND,NQUADS
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
      RETURN
   97 FORMAT('1',I4,' TH CLOSED-ORBIT 6-VECTOR WITH FRACTIONAL R.F.',
     +' FREQUENCY SHIFT OF: ',F11.7,///,1X,'NR. DISTANCE NAME',9X
     +,'DX*1000',3X,'DXP*1000',2X,'DY*1000',3X,'DYP*1000',2X,'DL*1000',
     +    2X,'DE/E*1000','  DE (MEV)'/)
C
C
   91 FORMAT(///,' MAX.ORB.DISTORTN: DXMAX =',F15.8,' MM. AT',F6.0,'M'/,
     +                          T20,'DYMAX =',F15.8,' MM. AT',F6.0,'M'/,
     +                          T20,'DSMAX =',F15.8,' MM. AT',F6.0,'M'/,
     +                          T20,'DEMAX =',F15.8,' MEV.',/,
     +           ' RMS ORB.DISTORTN: DXRMS =',F15.8,' MM. ',/,
     +           ' AT QUADS',   T20,'DYRMS =',F15.8,' MM. ',/,
     +                          T20,'DSRMS =',F15.8,' MM. ',/,
     +                          T20,'DERMS =',F15.8,' MEV.',/,
     +                          T20,'DYPRMS=',F15.8,' MRAD',/,
     +                          T20,'DYABS =',F15.8,' MM. ',/,
     +                          T20,'DYBEND=',F15.8,' MM. ',/,
     +           ' QUAD SAMPLES',T20,'NQUAD =',I15      )
      END