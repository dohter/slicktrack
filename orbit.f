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
      CHARACTER*8 TEMPNAME

C
C
C
      INCLUDE  "cdesy2.for"
      INCLUDE  "cloorb.for"
      INCLUDE  "cnlist.for"
      INCLUDE  "clatic.for"
      INCLUDE  "csol.for"
      INCLUDE  "ccentr.for"
C
C
      IPU=0
      NQUADS = 0
C
C
C=====ITERATE TO FIND CLOSED ORBIT
      WRITE(53,103)
  103 FORMAT(/,/,/,/,/,'  ')
      WRITE(53,'(A,A,A)')
     +   '1',' ORBIT: Get the final closed orbit'
      DO 23 ITORB=1,NITER
      IF(IE0.EQ.1.AND.IPTCO.EQ.1.AND.ITORB.EQ.NITER)
     +                                        WRITE(53,97)KKICK,FREQ
      IT=ITYPE(1)
C      CALL UCOPY(TMAT(1,1,IT),TREV,98)
      TREV =  TMAT(1:7,1:7,IT)
C


C
C=====CALCULATE 7X7 REVOLUTION MATRIX ON DESIGN ORBIT: INCLUDE RAD'N
C=====EFFECTS ONLY ON SECOND ITERATION.
      S=0.
      DO 26 II=2,NELEM
C      write(*,'(A,I10)')' NELEM=',II
      ITY=ITYPE(II)
      IID=ID(ITY)
      IF(ID(ITY).EQ.1.AND.XX(ITY).LT.1.D-10)GO TO 26 !Skip zero length drifts
      NTW=NTWIST(ITY)                                !Sept 2003. 
      IF(ID(ITY).NE.5.AND.ID(ITY).NE.6.AND.ID(ITY).NE.7.AND.
     + NAME(ITY)(1:2).NE.'CQ'.AND.NAME(ITY)(1:2).NE.'RQ'.AND.IID.NE.17)
     +                                                      S=S+YY(ITY)
C      write(*,'(2A,2I10)')' NELEM,NAME=',NAME(ITY),ITY,II
C      CALL UCOPY(TMAT(1,1,ITY),TM7,98)
      TM7  =  TMAT(1:7,1:7,ITY)
C      write(*,'(2A,2I10)')' NELEM,NAME=',NAME(ITY),ITY,II
C=====JUNE 97 as in:
C=====Oct 95:finally fix the kickers.Take code from SMILE3.JULY92.V2.
C=====For kickers, turn off curvature terms when getting the C.O.
      IF(IID.EQ.6.OR.IID.EQ.7)THEN
      TM7(5,1)=0.D0
      TM7(5,3)=0.D0
      TM7(2,6)=0.D0
      TM7(4,6)=0.D0
      ENDIF
      IF(IID.EQ.17.OR.NAME(ITY)(1:2).EQ.'BB')THEN !Ignore beam-beam kicks.
      TM7(2,1)=0.D0
      TM7(4,3)=0.D0
      ENDIF
C      write(*,'(2A,2I10)')' NELEM,NAME=',NAME(ITY),ITY,II
      IF(ITORB.NE.1)THEN
      DELX=(DX(II+1)+DX(II))/2.
      DELY=(DY(II+1)+DY(II))/2.
      GO TO (71,71,73,74,71,71,71,78,71,80,71,71,71,71,71,71,71),IID
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
      ENDIF

      CALL JAM777(TREV,TM7,TREV)
C
C
   26 CONTINUE
C
C
C
C
C======3/2004: clean out the 1000 factors and divisors in the calculation
c              to make it more transparent.
      DO 12 I=1,6
      DO 32 J=1,6
      TM6(I,J)=TREV(I,J)
   32 IF(I.EQ.J)TM6(I,J)=TREV(I,J)-1.D0
   12 VC6B(I)=-TREV(I,7)
      CALL SIMQ(TM6,VC6B,6,IERRO)
      IF(IERRO .NE. 0) STOP
C
C
C======Initialise statistical info. 
      DX(1) =VC6B(1)               ! Why are these DX(1) etc needed?
      DY(1) =VC6B(3)
      DXP(1)=VC6B(2)
      DYP(1)=VC6B(4)
      DL(1) =VC6B(5)
      DEL(1)=VC6B(6)
      DEN  = DABS(VC6B(6)) *E0
      DXMAX= DABS(VC6B(1))
      DYMAX= DABS(VC6B(3))
      DSMAX= DABS(VC6B(5))
      DENMAX=DEN    
      DXRMS   = 0.D0
      DYRMS   = 0.D0
      DXVQRMS = 0.D0
      DYVQRMS = 0.D0
      DXVCRMS = 0.D0
      DYVCRMS = 0.D0
      DSRMS= 0.D0
      DENRMS=0.D0
      DYPRMS=0.D0
      IBEND = 0
      DYBEND= 0
      DYABS = 0
      NQUADS   = 0
      NQUADSVQ = 0 
      NQUADSVC = 0
C
C
C
C
C======Calculate new closed orbit at ENTRANCE to element II.
C======DX(II-1) is DX at ENTRANCE to element (II-1) etc.
C======The distance S is at the end of element (II-1).
      S=0.D0
C      CALL VZERO(VC6A,12)
      VC6A = 0.D0
      DO 28 II=2,NELEM
C      IF(II.EQ.16.AND.ID(ITYPE(II)).EQ.7)XX(ITYPE(II)) = 0.D0 !With NITER = 1, kill 1st kicker
C      IF(II.EQ.16.AND.ID(ITYPE(II)).EQ.7)TMAT(4,7,ITYPE(II)) = 0.D0 ! To see its effect.
C      IF(II.EQ.29.AND.ID(ITYPE(II)).EQ.7)XX(ITYPE(II)) = 0.D0
C      IF(II.EQ.29.AND.ID(ITYPE(II)).EQ.7)TMAT(4,7,ITYPE(II)) = 0.D0

      JTY=ITYPE(II-1)
      IID=ID(JTY)
C=====IF(ID(JTY).EQ.1.AND.XX(JTY).LT.1.D-10)GO TO 52    !Skip zero length drifts 
      IF(ID(JTY).NE.5.AND.ID(JTY).NE.6.AND.ID(JTY).NE.7.AND.
     + NAME(JTY)(1:2).NE.'CQ'.AND.NAME(JTY)(1:2).NE.'RQ'.AND.IID.NE.17)
     +                                                 S=S+YY(JTY)
      MTW=NTWIST(JTY)
C      CALL UCOPY(TMAT(1,1,JTY),TEMP,98)
      TEMP =  TMAT(1:7,1:7,JTY)
      IF(IID.EQ.6.OR.IID.EQ.7)THEN          !03/04:Kill generation of dispersion.
      TEMP(5,1)=0.D0                        !That is only needed for 6x6 wrt C.O.  
      TEMP(5,3)=0.D0
      TEMP(2,6)=0.D0
      TEMP(4,6)=0.D0
      ENDIF
      IF(IID.EQ.17.OR.NAME(JTY)(1:2).EQ.'BB')THEN !Ignore beam-beam kicks.
      TM7(2,1)=0.D0
      TM7(4,3)=0.D0
      ENDIF


      IF(ITORB.NE.1)THEN
C     IF(IRAD.EQ.0)GO TO 51
      DELX=(DX(II-1)+DX(II))/2.
      DELY=(DY(II-1)+DY(II))/2.
      GO TO (51,51,53,54,51,51,51,58,51,60,51,51,51,51,51,51,51),IID
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

   51 CONTINUE 
      ENDIF

C======Get the C.O. at the end of element II-1, i.e. at the entrance to element II. 
      DO 210 I=1,6
      VC6A(I)=TEMP(I,7)
      DO 210 J=1,6
  210 VC6A(I)=VC6A(I)+TEMP(I,J)*VC6B(J)

   52 CONTINUE 
      DO 211 I=1,6
  211 VC6B(I)=VC6A(I)
      DX(II) =VC6B(1)
      DXP(II)=VC6B(2)
      DY(II) =VC6B(3)
      DYP(II)=VC6B(4)
      DL(II) =VC6B(5)
      DEL(II)=VC6B(6)
      DEN= DEL(II)*E0
C
      IF(ITORB.EQ.NITER)THEN
      IF(DABS(VC6B(1)).GT.DXMAX)THEN
      SXMAX=S
      DXMAX=DABS(VC6B(1))
      ENDIF
      IF(DABS(VC6B(3)).GT.DYMAX)THEN
      SYMAX=S
      DYMAX=DABS(VC6B(3))
      ENDIF
      IF(DABS(VC6B(5)).GT.DSMAX)THEN
      SSMAX=S
      DSMAX=DABS(VC6B(5))
      ENDIF
      IF(DABS(DEN).GT.DENMAX)DENMAX=DABS(DEN)

      IF(IID.EQ.5)U0=U0+XX(JTY)*E0*VC6B(5)
      IF(ID(ITYPE(II)).EQ.3.AND.
     + NAME(ITYPE(II)).NE.'EDGE'.AND.NAME(ITYPE(II))(1:2).NE.'CQ')THEN
      DXRMS  = DXRMS+VC6B(1)**2   
      DYRMS  = DYRMS+VC6B(3)**2   
      DSRMS  = DSRMS+VC6B(5)**2   
      DENRMS = DENRMS+DEN**2     
      DYPRMS = DYPRMS+VC6B(4)**2 
      DYABS  = DYABS+DABS(VC6B(3))
      NQUADS = NQUADS + 1
      ENDIF
      IF(ID(ITYPE(II)).EQ.7.AND.NAME(ITYPE(II))(1:2).EQ.'VQ')THEN
      DXVQRMS  = DXVQRMS+VC6B(1)**2   
      DYVQRMS  = DYVQRMS+VC6B(3)**2   
      NQUADSVQ = NQUADSVQ + 1
      ENDIF
      IF(ID(ITYPE(II)).EQ.7.AND.NAME(ITYPE(II))(1:2).EQ.'VC')THEN
      DXVCRMS  = DXVCRMS+VC6B(1)**2   
      DYVCRMS  = DYVCRMS+VC6B(3)**2   
      NQUADSVC = NQUADSVC + 1
      ENDIF

      IF(IID.EQ.2.OR.IID.EQ.15)THEN
      IBEND=IBEND+1
      DYBEND=DYBEND+VC6B(3)**2 
      ENDIF
C
      IF(IE0.EQ. 1.AND.ICLBP.NE.0.AND.ITORB.EQ.NITER.AND.(NAME(JTY).EQ.
     +'IP'.OR.IID.EQ.7.OR.IID.EQ.6.OR.NAME(JTY).EQ.'PU'.OR.
     +NAME(JTY)(1:2).EQ.'HC'.OR.NAME(JTY)(1:2).EQ.'VC'.OR.
     +NAME(JTY)(1:2).EQ.'HQ'.OR.NAME(JTY)(1:2).EQ.'VQ').AND.NAME(JTY)
     +.NE.'VKICK')
     +               WRITE(53,98)II,S,NAME(JTY),VC6B*1.D3,DEN,NAME(JTY)
C
   98 FORMAT(1X,I6,1X,F10.4,1X,A8,2X,7F10.4,2X,A8)
      QMARK = 0.D0                          ! Mark the positions of elements.
C-----------------------------------------------------------------------
C      IF((ID(ITYPE(II)).EQ.3.OR.ID(ITYPE(II)).EQ.15.OR.
C     +  NAME(ITYPE(II))(1:2).EQ.'VQ').AND.NAME(ITYPE(II))(1:2).NE.'CQ')
C     +                                                              THEN
C      PLOTPOS = CENPOS(II)-YY(ITYPE(II))/2.D0
C      IF(NAME(ITYPE(II))(1:2).EQ.'VQ')PLOTPOS = CENPOS(II)
C      WRITE(55,938)PLOTPOS,CENPOS(II),NAME(ITYPE(II)),II,ID(ITYPE(II)),
C     +             ZW1,ZW2,ZW3,DMX,DMY,DMZ,DDM,DLX,DLY,DLZ,DDL,
C     +             DDABS,DDLABS,DDMRUF,
C     +             XY71,XY72,XY712,XY73,XY74,XY734,QD7,
C     +             XY51,XY52,XY512,XY53,XY54,XY534,QD5,
C     +             D1,D3,CURVAT,SOLFLD,DX(II),DY(II),
C     +             BETACX,BETACY,QMARK
C  938 FORMAT(' ',2E16.6,1X,A8,1X,2I6,34(1X,E16.6))
      TEMPNAME = NAME(ITYPE(II))
      PLOTPOS = CENPOS(II)
      IF(NAME(ITYPE(II))(1:2).EQ.'VC')PLOTPOS = 10000.D0+CENPOS(II)
      IF(NAME(ITYPE(II))(1:2).EQ.'HC')PLOTPOS = 10000.D0+CENPOS(II)
      IF(NAME(ITYPE(II))(1:2).EQ.'VQ')PLOTPOS = 10000.D0+CENPOS(II)
      IF(NAME(ITYPE(II))(1:2).EQ.'HQ')PLOTPOS = 10000.D0+CENPOS(II)
      IF(NAME(ITYPE(II))(1:2).EQ.'RQ')PLOTPOS = 10000.D0+CENPOS(II)
      IF(NAME(ITYPE(II))(1:2).EQ.'CQ')PLOTPOS = 10000.D0+CENPOS(II)
C      WRITE(55,938)S,PLOTPOS,TEMPNAME,II,ID(ITYPE(II)),
C     +                          DX(II),DXP(II),DY(II),DYP(II),QMARK
  938 FORMAT(' ',2F16.6,1X,A8,1X,2I6,5(1X,E16.6))

C      ENDIF
C-------------------------------------------------------------------------

C
C
C=====DESY-II: STORE ORBIT AT THE PICKUPS: USE MILLIMETERS.
      IF(ITORB.EQ.NITER.AND.NAME(JTY).EQ.'PU')THEN
      IF(NAME(JTY).EQ.'PU')IPU=IPU+1
      LAG(2,IPU)=VC6B(1)
      LAG(1,IPU)=VC6B(3)
      ENDIF
C
C
      ENDIF

C======End of loop over elements.
   28 CONTINUE
C
C


C=====END CLOSED ORBIT ITERATION
   23 CONTINUE

C
C
C
C
C
      IF(NQUADS.GT.0)THEN
      DXRMS   = DSQRT(DXRMS/NQUADS)
      DYRMS   = DSQRT(DYRMS/NQUADS)
      DSRMS= DSQRT(DSRMS/NQUADS)
      DENRMS=DSQRT(DENRMS/NQUADS)
      DYPRMS=DSQRT(DYPRMS/NQUADS)
      DYABS= DYABS/NQUADS
      ENDIF
      IF(NQUADSVQ.GT.0)THEN
      DXVQRMS = DSQRT(DXVQRMS/NQUADSVQ)
      DYVQRMS = DSQRT(DYVQRMS/NQUADSVQ)
      ENDIF
      IF(NQUADSVC.GT.0)THEN
      DXVCRMS = DSQRT(DXVCRMS/NQUADSVC)
      DYVCRMS = DSQRT(DYVCRMS/NQUADSVC)
      ENDIF

      WRITE(53,91)DXMAX   *1.D3,
     +            DYMAX   *1.D3,
     +            DSMAX   *1.D3,
     +            DENMAX  *1.D3,
     +            DXRMS   *1.D3,
     +            DXVQRMS *1.D3,
     +            DXVCRMS *1.D3,
     +            DYRMS   *1.D3,
     +            DYVQRMS *1.D3,
     +            DYVCRMS *1.D3,
     +            DSRMS   *1.D3,
     +            DENRMS  *1.D3,
     +            DYPRMS  *1.D3,
     +            DYABS   *1.D3,
     +            DYBEND  *1.D6,
     +            NQUADS       ,
     +            NQUADSVQ     ,   
     +            NQUADSVC
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
   97 FORMAT('0',I4,' TH CLOSED-ORBIT 6-VECTOR WITH FRACTIONAL R.F.',
     +' FREQUENCY SHIFT OF: ',F11.7,//,5X,'NR.  DISTANCE NAME',10X
     +,'DX (mm)',3X,'DXP(mrad)',1X,'DY (mm)',3X,'DYP(mrad)',1X,
     + 'DL (mm)', 2X,'DE/E*1000','  DE (MEV)'/)
C
C
   91 FORMAT(///,' Statistics for the final orbit:',/,
     +           ' MAX.ORB.DISTORTN: DXMAX   =',F15.8,' MM. ',/,
     +                          T20,'DYMAX   =',F15.8,' MM. ',/,
     +                          T20,'DSMAX   =',F15.8,' MM. ',/,
     +                          T20,'DEMAX   =',F15.8,' MEV.',/,
     +           ' RMS ORB.DISTORTN: DXRMS   =',F15.8,' MM. ',
     +                          '    DXVQRMS =',F15.8,' MM. ',
     +                          '    DXVCRMS =',F15.8,' MM. ',/,
     +           ' AT QUADS'   ,T20,'DYRMS   =',F15.8,' MM. ',
     +                          '    DYVQRMS =',F15.8,' MM. ',
     +                          '    DYVCRMS =',F15.8,' MM. ',/,
     +                          T20,'DSRMS   =',F15.8,' MM. ',/,
     +                          T20,'DERMS   =',F15.8,' MEV.',/,
     +                          T20,'DYPRMS  =',F15.8,' MRAD',/,
     +                          T20,'DYABS   =',F15.8,' MM. ',/,
     +                          T20,'DYBEND  =',F15.8,' MM. ',/,
     +          ' QUAD SAMPLES',T20,'NQUAD   =',I15,/,
     +          ' VQ   SAMPLES',T20,'NQUADVQ =',I15,/,
     +          ' VC   SAMPLES',T20,'NQUADVC =',I15 )



      END
