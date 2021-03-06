      SUBROUTINE FIXORB2(E0,IE0,U0,CRAD,FREQ,KKICK)
C
C
C======CALCULATE THE CLOSED-ORBIT 6-VECTOR and correct to get minimum rms
C      deviation from the measurements.
C      Adapted from orbit.f
C      Simulate shift and scale errors on the monitors.
C      Constrain the size of the correction kicks.
C      This is similar to FIXORB1 but includes sextupoles on the first iteration and
C      is used once FIXORB1 has done a first sufficient orbit correction.

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C      PARAMETER (NM = 620, NC= 620)         !Number of measurements and correcters
      PARAMETER (NM = 977, NC= 977)
      PARAMETER (LWORK = NM)
      DIMENSION TREV(7,7),TEMP(7,7),TM7(7,7),VC6A(6),VC6B(6),TM6(6,6)
C======Storage for the least squares NAG routine.
      DIMENSION A(2*NM,NC),AUSE(2*NM,NC),C(NC,1),BMEAS(2*NM,1),
     +        QR(2*NM,NC), ALPHA(NC),BREAL(2*NM,1),BTEMP(2*NM,1),
     +        E(NC),Y(NC),Z(NC),R(2*NM),IPIV(NC),RESID(2*NM),AA(NM,NC),
     +        XXSAVE(NC),XXKEEP(NC),
     +        WORK(LWORK)
      DIMENSION   NVCPOS(NC)
      CHARACTER*8 NVNAME(NC)
C
C
C
      INCLUDE  "cdesy2.for"
      INCLUDE  "cloorb.for"
      INCLUDE  "cnlist.for"
      INCLUDE  "clatic.for"
      INCLUDE  "csol.for"
      INCLUDE  "cglist.for"
C
C

      WRITE(53,103)
  103 FORMAT(/,/,/,/,/,'  ')
      WRITE(53,'(A,A,A)')
     +   '1',' FIXORB2: Correct the orbit at the monitors at the',
     +       '  correction coils'
C======Loop over the VC correctors, giving them sequencially non-zero strengths to set up the
C      matrix A of perturbed closed orbits. Reset the strengths to zero at the end of the loop.

C======First save old VC settings and perhaps set the VC's to zero for tests.
      NVCCNT = 0
      DO 76 MA = 1,NELEM
      MATY=ITYPE(MA)
      IF(NAME(MATY)(1:2).EQ.'VC')THEN
      NVCCNT = NVCCNT + 1
      XXKEEP(NVCCNT) = XX(MATY)             !Save the old values.
C      XX(MATY)       = 0.D0
C      TMAT(4,7,MATY) =-XX(MATY)
C      write(*,'(A,A,I10,F16.5)')' ',NAME(MATY),MA,XX(MATY)
      ENDIF
 76   CONTINUE



      A = 0.D0
      BMEAS=0.D0
      BREAL=0.D0
      KG6 = 0
      KG7 = 0

      DELXX  = 1.D-6                  !Perturb by 1 microradian.
      NVCCNT = 0
      DO 36 NA = 0,NELEM              !NA=0 for the starting,uncorrected orbit.
      IF(NA.EQ.0)GO TO 37             !Skip the perturbations for the first pass.
      NDO = 0

      NATY=ITYPE(NA)
      IF(NAME(NATY)(1:2).EQ.'VC')THEN
      NVCCNT = NVCCNT + 1
      NVCPOS(NVCCNT) = NATY
      NVNAME(NVCCNT) = NAME(NATY)
      XXSAVE(NVCCNT) = XX(NATY)             !Save the old values.
      XX(NATY)       = XXSAVE(NVCCNT) + DELXX
      TMAT(4,7,NATY) =-XX(NATY)
      NDO = 1

C======WRITE(53,103)
C  103 FORMAT(/,/,/,/,/,'  ')
C      WRITE(53,'(A,A,1X,I6,1X,I6,1X,I6,1X,F15.9)')' ',NAME(NATY),
C     +                                   NA,NATY,NVCCNT,XX(NATY)
      ENDIF

C======Now get the new orbit with this (NVCCNT) corrector set.
      IF(NDO.EQ.0)GO TO 36

C
C
   37 CONTINUE
C======ITERATE TO FIND CLOSED ORBIT:
C      For minimising the orbit, switch off the iteration initially.
C      In contrast to FIXORB1, the sextupoles are now included since it is safe to
C      include them once the C.O. is already good.
C      Is the energy sawtooth a problem? Ignore it for now (3/2004).
C

      MITER = 1                              !To force just one iteration.
      DO 23 ITORB=1,MITER
      IT=ITYPE(1)
C      CALL UCOPY(TMAT(1,1,IT),TREV,98)
      TREV =  TMAT(1:7,1:7,IT)
C
C======CALCULATE 7X7 REVOLUTION MATRIX ON DESIGN ORBIT:
C      Perhaps include rad'n on a second iteration.

      S=0.
      DO 26 II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      IF(IID.EQ.1.AND.XX(ITY).LT.1.D-10)GO TO 26 !Skip zero length drifts
      IF(IID.NE.5.AND.IID.NE.6.AND.IID.NE.7.AND.IID.NE.17)S=S+YY(ITY)
C      CALL UCOPY(TMAT(1,1,ITY),TM7,98)
      TM7 = TMAT(1:7,1:7,ITY)
      IF(IID.EQ.6.OR.IID.EQ.7)THEN      !Kill generation of dispersion.
      TM7(5,1)=0.D0                     !That is only needed for 6x6 wrt C.O.
      TM7(5,3)=0.D0
      TM7(2,6)=0.D0
      TM7(4,6)=0.D0
      ENDIF
      IF(IID.EQ.17.OR.NAME(ITY)(1:2).EQ.'BB')THEN !Ignore beam-beam kicks.
      TM7(2,1)=0.D0
      TM7(4,3)=0.D0
      ENDIF

C======IF(ITORB.NE.1)THEN                   ! On first iteration, skip nonlinearities.
C======IF(IRAD.EQ.0)GO TO 71
      DELX=(DX(II+1)+DX(II))/2.D0      !03/04: Average from end of last to end of this
      DELY=(DY(II+1)+DY(II))/2.D0
C      IF(NA.EQ.0)THEN
C      IF(IID.EQ.7.AND.NAME(ITY)(1:2).EQ.'VC')
C     +            write(53,'(A,I10,2F20.6)')' VC  . ',II,DX(II),DY(II)
C      IF(IID.EQ.8)write(53,'(A,I10,2F20.6)')' Sext. ',II,DX(II),DY(II)
C      IF(IID.EQ.3)write(53,'(A,I10,2F20.6)')' Quad. ',II,DX(II),DY(II)
C      ENDIF
      GO TO (71,71,73,74,71,71,71,78,71,80,71,71,71,71,71,71,71),IID
   73 CONTINUE
      IF(IRAD.NE.0)TM7(6,7)=-CRAD*XX(ITY)**2*(DELX**2+DELY**2)/YY(ITY)
      GO TO 71
   74 IF(IRAD.NE.0)TM7(6,7)=-CRAD*XX(ITY)**2*(DELX**2+DELY**2)/YY(ITY)
      GO TO 71
   78 CALL SEXEQ(IRAD,CRAD,II,ITY,TM7)
C======AS IN DECEMBER 1995: AS IN EMITNC TURN OFF RADIATION IN SOLENOIDS.
   80 CONTINUE
      GO TO 71
   71 CONTINUE
C======ENDIF

      TREV = MATMUL(TM7,TREV)
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
   12 VC6B(I)=-TREV(I,7)
      CALL SIMQ(TM6,VC6B,6,IERRO)
      IF(IERRO .NE. 0) THEN
      WRITE(53,'(A)')' Problem getting C.O. in FIXORB2, so STOP'
      ENDIF
C
C======Initialise statistical info.
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
C======Calculate new closed orbit at ENTRANCE to element II.
C======DX(II-1) is DX at ENTRANCE to element (II-1) etc.
C      And store the C.O. at the monitors -- which coincide with the VC's
      S=0.D0
C      CALL VZERO(VC6A,12)
      VC6A = 0.D0
      NVM = 0

      DO 28 II=2,NELEM
      JTY=ITYPE(II-1)
      IID=ID(JTY)
C=====IF(IID.EQ.1.AND.XX(JTY).LT.1.D-10)GO TO 52 !Skip zero length drifts
      IF(IID.NE.5.AND.IID.NE.6.AND.IID.NE.7.AND.IID.NE.17)S=S+YY(JTY)
C      CALL UCOPY(TMAT(1,1,JTY),TEMP,98)
      TEMP =  TMAT(1:7,1:7,JTY)
      IF(IID.EQ.6.OR.IID.EQ.7)THEN          !Kill generation of dispersion.
      TEMP(5,1)=0.D0                        !That is only needed for 6x6 wrt C.O.
      TEMP(5,3)=0.D0
      TEMP(2,6)=0.D0
      TEMP(4,6)=0.D0
      ENDIF
      IF(IID.EQ.17.OR.NAME(JTY)(1:2).EQ.'BB')THEN !Ignore beam-beam kicks.
      TM7(2,1)=0.D0
      TM7(4,3)=0.D0
      ENDIF

C======IF(ITORB.NE.1)THEN           ! On first iteration, skip nonlinearities.
      DELX=(DX(II-1)+DX(II))/2.D0
      DELY=(DY(II-1)+DY(II))/2.D0
      GO TO (51,51,53,54,51,51,51,58,51,60,51,51,51,51,51,51,51),IID
   53 CONTINUE
      IF(IRAD.NE.0)TEMP(6,7)=-CRAD*XX(JTY)**2*(DELX**2+DELY**2)/YY(JTY)
      GO TO 51
   54 IF(IRAD.NE.0)TEMP(6,7)=-CRAD*XX(JTY)**2*(DELX**2+DELY**2)/YY(JTY)
      GO TO 51
   58 CALL SEXEQ(IRAD,CRAD,II-1,JTY,TEMP)
      GO TO 51
   60 CONTINUE
   51 CONTINUE
C======ENDIF
C
C======Get the C.O. at the end of element II-1, i.e. at the entrance to element II.
      DO 210 I=1,6
      VC6A(I)=TEMP(I,7)
      DO 210 J=1,6
  210 VC6A(I)=VC6A(I)+TEMP(I,J)*VC6B(J)



      DO 211 I=1,6
  211 VC6B(I)=VC6A(I)
      DX(II) =VC6B(1)
      DXP(II)=VC6B(2)
      DY(II) =VC6B(3)
      DYP(II)=VC6B(4)
      DL(II) =VC6B(5)
      DEL(II)=VC6B(6)
      DEN= DEL(II)*E0
      IF(DABS(VC6B(1)).GT.DXMAX)DXMAX=DABS(VC6B(1))
      IF(DABS(VC6B(3)).GT.DYMAX)DYMAX=DABS(VC6B(3))
      IF(DABS(VC6B(5)).GT.DSMAX)DSMAX=DABS(VC6B(5))
      IF(DABS(DEN)    .GT.DENMAX)DENMAX=DABS(DEN)
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

C======Set up the ``measured'' orbit (B) and SENSITIVITY matrix (A) at the monitors
C      which we take to be at the correctors VC.
C      NVM labels the measurements. NVCCNT labels the perturbing kicks.

C======Give the ``measurements'' a shift error (G6) and a scale error (G7).
C      The orbit recorded by a monitor is: (real orbit - offset)*(1 + scale error).
C      We want to find corrector settings which give a C.O. which cancels the existing C.O.

C      IF(NAME(JTY)(1:2).EQ.'VC')THEN
C======DY(II) is the V.C.O. just before the element II.
      IF(ID(ITYPE(II)).EQ.7.AND.NAME(ITYPE(II))(1:2).EQ.'VC')THEN
      NVM=NVM + 1
      IF(NA.EQ.0)THEN
      KG6 = KG6 + 1
      KG7 = KG7 + 1
      BMEAS(NVM,1) = (DY(II)- G6(KG6)*0.001D0)*(1.D0 + G7(KG7))
      BREAL(NVM,1) =  DY(II)
      ENDIF
      IF(NA.NE.0)THEN
      A(NVM,NVCCNT) = (DY(II)-BREAL(NVM,1))/DELXX
C======WRITE(53,'(A,A,1X,I10,1X,F10.5)')' ',NAME(JTY),NVCCNT,
C     +                                     A(NVM,NVCCNT)
      ENDIF
      ENDIF

C======IF(ITORB.EQ.MITER)THEN
C      IF(IID.EQ.5)U0=U0+XX(JTY)*E0*VC6B(5)
C      ENDIF

C======End of loop over elements.
   28 CONTINUE

      IF(NA.EQ.0)THEN
      DXRMS   = DSQRT(DXRMS/NQUADS)
      DYRMS   = DSQRT(DYRMS/NQUADS)
      DXVQRMS = DSQRT(DXVQRMS/NQUADSVQ)
      DYVQRMS = DSQRT(DYVQRMS/NQUADSVQ)
      DXVCRMS = DSQRT(DXVCRMS/NQUADSVC)
      DYVCRMS = DSQRT(DYVCRMS/NQUADSVC)
      DSRMS= DSQRT(DSRMS/NQUADS)
      DENRMS=DSQRT(DENRMS/NQUADS)
      DYPRMS=DSQRT(DYPRMS/NQUADS)
      DYABS= DYABS/NQUADS

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
   91 FORMAT(///,' Statistics for the uncorrected orbit:',/,
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
C
      WRITE(53,1033)
 1033 FORMAT(/,/,'  ')
      ENDIF
C
C======End of iteration loop.
   23 CONTINUE
C
C
C
C
C
C
C======Reset this correction coil.
      IF(NA.NE.0)THEN
      XX(NATY)      = XXSAVE(NVCCNT)
      TMAT(4,7,NATY)=-XX(NATY)
      ENDIF


C======End of loop over vertical correctors.
   36 CONTINUE




C
C
C======Now find the VC strengths for a least squares fit.
C      Include a constraint to limit the size of the correction kicks.
      CALL UNIT(NM,AA)
      A(NVCCNT+1:2*NVCCNT,1:NVCCNT)   = AA(1:NVCCNT,1:NVCCNT)*1.D0
      BMEAS(NVCCNT+1:2*NVCCNT,1)      = -XXKEEP(1:NVCCNT)*1.D0
C      BMEAS(NVM+1:2*NVM,1)            = -XXKEEP(1:NVM)*1.D0

C======Instead of A, use AUSE, the actual matrix to be used for fitting once some monitors
C      have been switched off. We need A later to look at the residuals at the dead
C      monitors to see if they are indeed large.

C======Kill lines of AUSE to, in effect, switch off monitors. First set a complete AUSE
C      Killing lines actually amounts to giving some monitors zero weight in the fit.
C      Should we kill the strength constraint too?
C
C      Killing lines of AUSE, leaves the columns (i.e. the correctors) untouched,
C      i.e. there are still NVCCNT correctors.

      AUSE = A
      DO 42 IG = 1,NVCCNT
  42  AUSE(IG,1:NVCCNT) = G9(IG)* A(IG,1:NVCCNT)

      NMEAS = 1
      EPS = 1.D-15
      IFAIL = 0
      BTEMP = BMEAS
      CALL DGELS('N', 2*NVCCNT, NVCCNT, NMEAS, AUSE, 2*NM, BTEMP, 2*NM,
     +               WORK, LWORK, IFAIL)
      DO I = 1, NMEAS
        C(1:NVCCNT,I) = BTEMP(1:NVCCNT,I)
      ENDDO
c     CALL F04AMF(AUSE,2*NM,C,NC,BMEAS,2*NM,2*NVCCNT,NVCCNT,NMEAS,EPS,
c    +            QR,2*NM,ALPHA,E,Y,Z,R,IPIV,IFAIL)


C======Set up the new strengths and matrices of the correctors
      WRITE(53,1033)
      WRITE(53,'(A)')' The coil strengths needed for correction: '

      DO 38 IC = 1,NVCCNT
      WRITE(53,'(A,A,1X,A,I10,1X,I10,1X,F15.6)')
     + ' ','Coil name, pos1,pos2,value',
     +      NVNAME(IC),IC,NVCPOS(IC),C(IC,1)*1.D3
      XX(NVCPOS(IC)) = XXSAVE(IC) - C(IC,1) ! Add corrections to the originals.
      IF(DABS(C(IC,1)).LT.0.005D0)THEN      ! Note the minus sign.
      TMAT(4,7,NVCPOS(IC))= -XX(NVCPOS(IC))
      TMAT(4,6,NVCPOS(IC))= +XX(NVCPOS(IC))
      TMAT(5,3,NVCPOS(IC))= -XX(NVCPOS(IC))
      ELSE
      WRITE(53,1033)
      WRITE(53,'(A,A,1X,F10.6,1X,A)')
     + ' ','FIXORB2: Required correction strength', C(IC,1),
     +     ': Too large. STOP.'
C=======Also the radiation can be too strong
      STOP
      ENDIF
  38  CONTINUE

      WRITE(53,103)
      CORRMS = 0.D0
      TOTRMS = 0.D0
      DO 41 ICCC = 1,NVCCNT
      TOTRMS = TOTRMS + XX(NVCPOS(ICCC))*XX(NVCPOS(ICCC))
  41  CORRMS = CORRMS + C(ICCC,1)*C(ICCC,1)
      CORRMS = DSQRT(CORRMS/NVCCNT)*1.D3
      TOTRMS = DSQRT(TOTRMS/NVCCNT)*1.D3
      WRITE(53,'(A,A,1X,F15.6,A)')' ',
     +  'rms correction coil strength change  ', CORRMS,' mrad'
      WRITE(53,'(A,A,1X,F15.6,A)')' ',
     +  'rms correction coil total strength   ', TOTRMS,' mrad'
C
C
      WRITE(53,1033)
      NNG8 = 0
      DO 39 IM = 1,NVCCNT
      RESID(IM) = 0.D0
      DO 40 ICC = 1,NVCCNT
      RESID(IM)= RESID(IM) + A(IM,ICC)*C(ICC,1)
   40 CONTINUE
      RESID(IM)= (RESID(IM) - BMEAS(IM,1))*1000.D0
      NG8 = G8(IM)
      IF(NG8.EQ.0)NNG8 = NNG8 + 1
      WRITE(53,'(A,1X,I10,1X,F15.6,1X,I4,1X,I4,1X,F15.6)')
     +' Mon.number, value(mm), on/off, cumulative off, residual(mm)',
     +                       IM,BMEAS(IM,1)*1.D3,NG8,NNG8,RESID(IM)
   39 CONTINUE


      RETURN


   98 FORMAT(1X,I4,1X,F10.4,1X,A8,2X,7F10.4,2X,A8)




      END
