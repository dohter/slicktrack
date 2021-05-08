      SUBROUTINE FIXORB3(E0,IE0,U0,CRAD,FREQ,KKICK)
C
C
C
C======Use Kick Minimisation to minimise the generation of vertical dispersion and the
C      tilt of n_0 in the arcs (but not necessarily in the quads!).
C
C======This is done by calculating the closed orbit 6-vector and correcting it to get minimum
C      rms combined vertical kick from the correction coil and the adjacent quadrupole.
C      Simulate shift and scale errors on the monitors.
C      The bare algorithm is very unstable, the physics of which is obvious.
C      So constrain
C      1) the size of the correction kicks.
C      2) the VCO excursion as in FIXORB2.
C      The response matrix therefore has 3 sections, each of which can have its own weight.
C      It is assumed that FIXORB1 has done a first sufficient orbit correction.
C
C
C      Care with the signs of kicks is more important than usual.
C      A quad with +ve strength XXq causes vertical defocusing and if the betatron
C      orbit is at Y, the upward kick is XXq*Y.
C      The signs for the strengths of the vertical kickers follow the Mais-Ripken
C      convention (i.e.not the Chao convention). So  +ve strength XXk causes a downward kick.
C      and TMAT(4,7,.) = -XXk
C      The total vertical kick to the C.O. of a quad-corrector pair is then
C      XXq*YCO - XXk . Note also the relative signs for VCO effects in quads in TSPIN.
C
C      As mentioned in LATTIN, if a quad with strength XXq is shifted vertically
C      by Dq, the quad in efect gives an extra kick of - XXq*Dq.
C
C      If the monitor at the corrector is shifted by Dm vertically it registers
C      a position YCO - Dm.  The estimate of the kick in the quadrupole is then
C      changed  by -XXq*Dm. If the monitor is then shifted by the same amount as the
C      quadrupole, so that Dm = Dq, the kick in the quadrupole is still given
C      correctly by XXq*(YCO - Dm), i.e. in terms of the measured(!) C.O.
C      Then, if the monitor and the quad are perfectly relatively aligned, the
C      vertical kick can be compensated exactly in the approximation that the
C      quadrupole is thin.
C
C      If the monitor is shifted up by Dmq wrt the quad, the measured C.O. is too
C      small by Dmq and the estimated kick from the quad will then be too small
C      by XXq*Dmq.
C
C      In any case, if Dmq is known exactly from beam based calibration, the kick
C      from the quad can be found exactly and the exact value of Dmq is irrelevant.
C      However, uncertainty in the quad kick can come from uncertainly in the measured
C      Dmq.
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C      PARAMETER (NM = 620, NC= 620)         !Number of measurements and correcters
      PARAMETER (NM = 977, NC= 977)
      PARAMETER (LWORK = NM)
      DIMENSION TREV(7,7),TEMP(7,7),TM7(7,7),VC6A(6),VC6B(6),TM6(6,6)
C======Storage for the response matrices.
      DIMENSION A1(NM,NC),A2(NM,NC)
      DIMENSION AUPPER(NM,NC),AMIDDLE(NM,NC),ALOWER(NM,NC)
C======Storage for the least squares NAG routine.
      DIMENSION AUSE(3*NM,NC),C(NC,1),
     +        BMEAS(3*NM,1),BQUAD(NM,1),BQUADP(NM,1),
     +        QR(3*NM,NC), ALPHA(NC),
     +        E(NC),Y(NC),Z(NC),R(3*NM),IPIV(NC),RESID(3*NM),AA(NM,NC),
     +        XXSAVE(NC),XXKEEP(NC),XXDIFF(NC),
     +        BTEMP(3*NM,1), WORK(LWORK)
      DIMENSION ORBMEAS1(NM),ORBMEAS2(NM),ORBREAL(NM)
      DIMENSION   NVCPOS(NC)
      CHARACTER*8 NVNAME(NC)
      DIMENSION NBIND(NM), IBIND(NM,10)
      DIMENSION WKSPACE(NC),ATEMP(NM,NC)
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
      WEIGHT1 = 0.D0
      WEIGHT2 = 1.D0
      WEIGHT3 = 1.D0


      WRITE(53,103)
  103 FORMAT(/,/,/,/,/,'  ')
      WRITE(53,'(A,A,A,A)')
     +   '1',' FIXORB3: Correct the combined vertical kicks of',
     +       ' correctors and their adjacent quadrupoles.'

C======First loop over all elements to set up the connection between vertical correctors
C      and nearby relevant magnets.
C      This algorithm depends on details of the layout. So it must be checked for each
C      new optic/layout.

      NVCCNT = 0
      RMSORIG = 0.D0
      DO 46 NV = 1,NELEM
      NVTY=ITYPE(NV)
      IVD =ID(NVTY)
      IF(IVD.NE.7.OR.NAME(NVTY)(1:2).NE.'VC')GO TO 46
C======A vertical corrector has been found. Now search for relevant nearby elements
      NVCCNT = NVCCNT + 1
C      WRITE(53,1033)
C      WRITE(53,'(A,A,3I10)')' Corrector found    :  ',
C     +                      NAME(NVTY),IVD,NV,NVCCNT
C======Search forward to find up to 4 bound elements. Try to be independent of known
C      ordering.

C======To test the consistency of the correction, one can set the incoming strength of the
C      VC to ZERO.
      XXKEEP(NVCCNT)  =  XX(NVTY)
      RMSORIG         =  RMSORIG + XX(NVTY)**2
C      XX(NVTY)       =  0.D0
C      TMAT(4,7,NVTY) = -XX(NVTY)

      NBIND(NVCCNT) = 0
      NQ=0
      NVQ=0
      NSP = NV + 1
      NSPP= NV +10
      IF(NSP .GT.NELEM)NSP =NELEM
      IF(NSPP.GT.NELEM)NSPP=NELEM
      DO 55 NS = NSP,NSPP
      NSTY=ITYPE(NS)
      ISD =ID(NSTY)
      IF(ISD.EQ.3.OR.ISD.EQ.15.OR.ISD.EQ.7)THEN
      NBIND(NVCCNT) = NBIND(NVCCNT) + 1
      IBIND(NVCCNT,NBIND(NVCCNT)) = NS
      IF(ISD.EQ.3.OR.ISD.EQ.15)NQ=NQ+1
      IF(ISD.EQ.7)NVQ=NVQ+1
C      WRITE(53,'(A,I10,3X,A,3I10)')' Bound element:  ', NS,
C     +     NAME(NSTY),ISD,NBIND(NVCCNT),IBIND(NVCCNT,NBIND(NVCCNT))
      ENDIF
      IF(NBIND(NVCCNT).EQ.4.AND.NQ.EQ.3.AND.NVQ.EQ.1)GO TO 461 !Bindings found.
  55  CONTINUE

C      WRITE(53,1033)
C======If no bound elements have been found, search backwards to find up to 4 bound elements.
      NBIND(NVCCNT) = 0
      NQ=0
      NVQ=0
C      WRITE(53,'(A,A,3I10)')' Corrector found    :  ',
C     +                      NAME(NVTY),IVD,NV,NVCCNT
      NSM = NV - 1
      NSMM= NV - 10
      IF(NSM .LT.1)NSM =1
      IF(NSMM.LT.1)NSMM=1
      DO 56 NS = NSM,NSMM,-1
      NSTY=ITYPE(NS)
      ISD =ID(NSTY)
      IF(ISD.EQ.3.OR.ISD.EQ.15.OR.ISD.EQ.7)THEN
      NBIND(NVCCNT) = NBIND(NVCCNT) + 1
      IBIND(NVCCNT,NBIND(NVCCNT)) = NS
      IF(ISD.EQ.3.OR.ISD.EQ.15)NQ=NQ+1
      IF(ISD.EQ.7)NVQ=NVQ+1
C      WRITE(53,'(A,I10,3X,A,3I10)')' Bound element:  ', NS,
C     +     NAME(NSTY),ISD,NBIND(NVCCNT),IBIND(NVCCNT,NBIND(NVCCNT))
      ENDIF
      IF(NBIND(NVCCNT).EQ.4.AND.NQ.EQ.3.AND.NVQ.EQ.1)GO TO 461 !Bindings found.
  56  CONTINUE

C======All bound elements for this corrector have been found.
C 461  WRITE(53,1033)
 461  CONTINUE


C=======End of element loop.
 46   CONTINUE



C======Loop over the VC correctors, giving them sequentially shifted strengths to set up the
C      response matrix A, i.e. the Jacobean of perturbed combined quad-corrector kicks.
C      Reset the perturbations to zero at the  end of the loop.


      A = 0.D0
      BMEAS=0.D0
      BQUAD=0.D0
      BQUADP=0.D0
      KG6 = 0
      KG7 = 0
      KG8 = 0

      DELXX  = 1.D-6                  !Perturb by 1 microradian.  1.D-6
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
C      Then the sextupoles are not included even if the initial orbit from FIXORB1 is OK.
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
      TM7 =  TMAT(1:7,1:7,ITY)
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

      IF(ITORB.NE.1)THEN                   ! The first iteration, skips nonlinearities.
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
      ENDIF

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
      WRITE(53,'(A)')' Problem getting C.O. in FIXORB4, so STOP'
      ENDIF
C
C======Initialise statistical info.
      DEN  = DABS(VC6B(6)) *E0
      DXMAX= DABS(VC6B(1))
      DYMAX= DABS(VC6B(3))
      DSMAX= DABS(VC6B(5))
      DENMAX=DEN
      DXRMS    = 0.D0
      DYRMS    = 0.D0
      DXVQRMS  = 0.D0
      DYVQRMS  = 0.D0
      DXVCRMS  = 0.D0
      DYVCRMS  = 0.D0
      DSRMS    = 0.D0
      DENRMS   = 0.D0
      DYPRMS   = 0.D0
      IBEND = 0
      DYBEND= 0
      DYABS = 0
      NQUADS   = 0
      NQUADSVQ = 0
      NQUADSVC = 0
      RMSKICK  = 0.D0
      RMSQUAD  = 0.D0
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

      IF(ITORB.NE.1)THEN           ! On first iteration, skip nonlinearities.
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
      ENDIF
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
C      IF(ID(ITYPE(II)).EQ.7.AND.NAME(ITYPE(II))(1:2).EQ.'VC')THEN
C      write(*,'(A,F12.6)')' DY(II)1',DY(II)*1.D3
C      ENDIF
C      IF(ID(ITYPE(II)).EQ.3.AND.NAME(ITYPE(II))(1:2).EQ.'CQ')THEN
C      write(*,'(A,F12.6)')' DY(II)3',DY(II)*1.D3
C      ENDIF
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


C======Set up the ``measured'' orbit (BMEAS) and RESPONSE matrix (A) for the combined
C      quad-corrector kicks.
C      NVM labels the measurements. NVCCNT labels the perturbing kicks.

C======Give the ``measurements'' a shift error (G6) and a scale error (G7).
C      The orbit recorded by a monitor is: (real orbit - offset)*(1 + scale error).
C      We want to find corrector settings which give a C.O. and correctors for
C      which the existing combined kicks are cancelled.

C      IF(NAME(JTY)(1:2).EQ.'VC')THEN
C======DY(II) is the V.C.O. just before the element II.
      IF(ID(ITYPE(II)).EQ.7.AND.NAME(ITYPE(II))(1:2).EQ.'VC')THEN
      NVM=NVM + 1

      IF(NA.EQ.0.AND.NVM.LE.1)write(53,1033)
      IF(NA.EQ.0.AND.NVM.EQ.1)WRITE(53,'(A)')' Specimen kick data'
      IF(NA.EQ.0.AND.NVM.LE.10)write(53,1033)
      IF(NA.EQ.0.AND.NVM.LE.10)
     +          write(53,'(A,A,3I10,F12.6)')' Corrector found    :  ',
     +         NAME(ITYPE(II)),ID(ITYPE(II)),II,NVM,XX(ITYPE(II))

      IF(NA.EQ.0)THEN
      KG6 = KG6 + 1
      KG7 = KG7 + 1
      KG8 = KG8 + 1

C======Now get the combined kick of a quad and the nearby kicker VC.
C      The shift of the quad is already embodied in the strength of the VQ and this
C      can be used to take into account the common shift of the monitor and the quad.
C      So the estimate of the total kick goes in two steps.
C      1)Simulate the offset between the monitor and the quad with the random shifts G8.
C      2)Include the effect of the overall common shift of the monitor and the quad
c        using the strength of the VQ.
C
C      Try to improve the precision by using TMAT(4,3,.) to represent the kick in the
C      quad.
      ORBMEAS1(NVM) = (DY(II)- G8(KG8)*0.001D0)
      ORBMEAS2(NVM) = (DY(II)- G6(KG6)*0.001D0)*(1.D0 + G7(KG7)*1.D0)
      ORBREAL(NVM) =  DY(II)
      BMEAS(NVM,1) = 0.D0
      BQUAD(NVM,1) = 0.D0
      QCOMB = 0.D0
      ICNT=0
C      Search for the nearby nearby quads:
      DO 66 JN = 1,NBIND(NVM)
      IK  =  IBIND(NVM,JN)
      KTY =  ITYPE(IK)
      IKD =  ID(KTY)
      IF(NVM.LE.10)
     +      write(53,'(A,I10,3X,A,3I10)')' Bound element:  ', II,
     +                                        NAME(KTY),IKD,JN,IK
      IF(IKD.EQ.3)THEN
C      QCOMB = QCOMB + XX(KTY)
      QCOMB = QCOMB + TMAT(4,3,KTY)
      IF(NVM.LE.10)write(53,'(A,F16.5)')'Strength: ',XX(KTY)
      ICNT =ICNT + 1
      ENDIF
      IF(IKD.EQ.15)THEN
C      QCOMB = QCOMB + X2(KTY)
      QCOMB = QCOMB + TMAT(4,3,KTY)
      IF(NVM.LE.10)write(53,'(A,F16.5)')'Strength: ',X2(KTY)
      ICNT =ICNT + 1
      ENDIF
      IF(IKD.EQ.7)THEN
      XXVQ = XX(KTY)
      IF(NVM.LE.10)write(53,'(A,F16.5)')'Strength: ',XX(KTY)
      ICNT =ICNT + 1
      ENDIF
 66   CONTINUE
      IF(ICNT.NE.4)THEN
      WRITE(53,'(A,I4,A,I3,A)')' Monitor ',NVM,
     +                ' Number of bound elements = :',ICNT,' So STOP.'
      STOP
      ENDIF
      BMEAS(NVM,1) = QCOMB*ORBMEAS1(NVM) - XX(ITYPE(II)) - XXVQ !Note the signs.
      BQUAD(NVM,1) = QCOMB*ORBREAL(NVM)     !For A only the quad GRADIENTS are needed.
      IF(NVM.LE.10)WRITE(53,'(A,I4,A,7F11.6)')' Monitor ',NVM,
     +              ' strengths (mrad)',ORBMEAS1(NVM)*1.D3,
     +                                  QCOMB*ORBMEAS1(NVM)*1.D3,
     +                                  XX(ITYPE(II))*1.D3,
     +                                  XXVQ*1.D3,
     +                                  BMEAS(NVM,1)*1.D3,
     +                                  BQUAD(NVM,1)*1.D3
      RMSKICK = RMSKICK + BMEAS(NVM,1)**2
      RMSQUAD = RMSQUAD + BQUAD(NVM,1)**2
      NVMKEEP = NVM
      ENDIF

      IF(NA.NE.0)THEN
      BQUADP(NVM,1) = 0.D0
      QCOMB = 0.D0
      DO 67 JN = 1,NBIND(NVM)
      IK  = IBIND(NVM,JN)
      KTY = ITYPE(IK)
      IKD = ID(KTY)
      IF(IKD.EQ.3)THEN
C      QCOMB = QCOMB + XX(KTY)
      QCOMB = QCOMB + TMAT(4,3,KTY)
      ENDIF
      IF(IKD.EQ.15)THEN
C      QCOMB = QCOMB + X2(KTY)
      QCOMB = QCOMB + TMAT(4,3,KTY)
      ENDIF
 67   CONTINUE
      BQUADP(NVM,1) = QCOMB*DY(II)

C======Response matrix for the kicks
      A1(NVM,NVCCNT) = (BQUADP(NVM,1)-BQUAD(NVM,1))/DELXX
C======Response matrix for the orbit.
      A2(NVM,NVCCNT) = (DY(II)-ORBREAL(NVM))/DELXX
C      IF(NVCCNT.LE.10)write(*,'(A,A,1X,I10,1X,I10,1X,F10.5)')' ',
C     +                           NAME(KTY),NVM,NVCCNT,A(NVM,NVCCNT)
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
      RMSKICK = DSQRT(RMSKICK/NQUADSVC)
      RMSQUAD = DSQRT(RMSQUAD/NQUADSVC)

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
      WRITE(53,'(A,A,1X,F15.6,A)')' ',
     +    'Corresponding rms combined kicks (mrad): ',RMSKICK*1.D3
      WRITE(53,'(A,A,1X,F15.6,A)')' ',
     +    'Original rms corrector strengths (mrad): ',RMSORIG*1.D3
      WRITE(53,'(A,A,1X,F15.6,A)')' ',
     +    'Original rms quad kicks          (mrad): ',RMSQUAD*1.D3

C
C
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
C      The full response matrix AUSE has 3 parts:
C      AUPPER  for the combined quad-corrector kicks,
C      AMIDDLE for the corrector strengths,
C      ALOWER  for the VCO.
C
C      AUPPER must be completed to include the contributions of the correctors themselves
C      by adding a unit matrix.
      CALL UNIT(NM,AA)
      AUSE = 0.D0
      AUPPER(1:NVM,1:NVCCNT) = A1(1:NVM,1:NVCCNT) - AA(1:NVM,1:NVCCNT)
      ATEMP = AUPPER                        !F03AAF overwrites the input matrix.
      WRITE(53,'(A,2I10)')' Going into NAG routines. NVM,NVCCNT: ',
     +                    NVM,NVCCNT
      IFAIL = 0
C======Check for pathologies.
      CALL F03AAF(ATEMP,NM,NC,DETA,WKSPACE,IFAIL)
      WRITE(53,'(A,I4)')' Get determinant in F03AAF: IFAIL= ',IFAIL
      IF(IFAIL.NE.0)THEN
      WRITE(53,'(A)')' Getting det(A) in F03AAF failed. So STOP.'
      STOP
      ENDIF
      WRITE(53,'(A,E12.6)')' Det(A)= ',DETA


      CALL UNIT(NM,AMIDDLE)
      AUSE(NVM+1:2*NVM,1:NVCCNT) =  AMIDDLE(1:NVCCNT,1:NVCCNT)*WEIGHT2
      BMEAS(NVM+1:2*NVM,1)       = -XXKEEP(1:NVM)*WEIGHT2


      ALOWER = A2
      BMEAS(2*NVM+1:3*NVM,1)       = ORBMEAS2(1:NVM)*WEIGHT3

C======Kill lines of AUSE to, in effect, switch off monitors.
C      Killing lines actually amounts to giving some monitors zero weight in the fit.
C      No need to kill the strength constraint too
C
C      Killing lines of AUSE, leaves the columns (i.e. the correctors) untouched,
C      i.e. there are still NVCCNT correctors.

      DO 42 IG = 1,NVM
  42  AUSE(IG,1:NVCCNT)       = G9(IG)* AUPPER(IG,1:NVCCNT)*WEIGHT1

      DO 43 IG = 1,NVM
  43  AUSE(2*NVM+IG,1:NVCCNT) = G9(IG)* ALOWER(IG,1:NVCCNT)*WEIGHT3


      NMEAS = 1
      EPS = 1.D-15
      IFAIL = 0
C======Get the set C needed to generate the original combined kicks.
      BTEMP = BMEAS
      CALL DGELS('N', 3*NVCCNT, NVCCNT, NMEAS, AUSE, 3*NM, BTEMP, 3*NM,
     +               WORK, LWORK, IFAIL)
      DO I = 1, NMEAS
        C(1:NVCCNT,I) = BTEMP(1:NVCCNT,I)
      ENDDO
C     CALL F04AMF(AUSE,3*NM,C,NC,BMEAS,3*NM,3*NVCCNT,NVCCNT,NMEAS,EPS,
C    +            QR,3*NM,ALPHA,E,Y,Z,R,IPIV,IFAIL)
      WRITE(53,'(A,I4)')' Fit in F04AMF: IFAIL= ',IFAIL
      IF(IFAIL.NE.0)THEN
      WRITE(53,'(A)')' Fit in F04AMF failed. So STOP.'
      STOP
      ENDIF

C======Set up the new strengths and matrices of the correctors
      WRITE(53,1033)
      WRITE(53,'(A)')' Specimen coil strengths needed for correction: '
      DO 38 IC = 1,NVCCNT
      XX(NVCPOS(IC)) = XXSAVE(IC) - C(IC,1) ! Correct the originals.
      IF(IC.LE.10)WRITE(53,'(A,1X,A,I10,1X,I10,4(1X,F12.6))')
     + ' Coil name, pos1,pos2,original,change,total new(mrad)',
     +   NVNAME(IC),IC,NVCPOS(IC),XXKEEP(IC)*1.D3,
     +   C(IC,1)*1.D3,XX(NVCPOS(IC))*1.D3
      IF(DABS(C(IC,1)).LT.0.005D0)THEN           ! wrt that from FIXORB1.
      TMAT(4,7,NVCPOS(IC))= -XX(NVCPOS(IC))
      TMAT(4,6,NVCPOS(IC))= +XX(NVCPOS(IC))
      TMAT(5,3,NVCPOS(IC))= -XX(NVCPOS(IC))
      ELSE
      WRITE(53,1033)
      WRITE(53,'(A,A,1X,F10.6,1X,A)')
     + ' ','FIXORB4: Required correction strength', C(IC,1),
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
      WRITE(53,'(A)')' Specimen fit data'
      NNG8 = 0
      DO 39 IM = 1,NVCCNT
      RESID(IM) = 0.D0
      DO 40 ICC = 1,NVCCNT
      RESID(IM)= RESID(IM) + AUSE(IM,ICC)*C(ICC,1)
   40 CONTINUE
      RESID(IM)= (RESID(IM) - BMEAS(IM,1))*1000.D0
      NG8 = G8(IM)
      IF(NG8.EQ.0)NNG8 = NNG8 + 1
      IF(IM.LE.10)THEN
      WRITE(53,'(A,A,A,1X,I10,1X,F12.6,1X,I4,1X,I4,1X,F12.6,
     +                                             1X,E12.3)')' ',
     +'Mon.number, value(mrad), on/off, cumulative off',
     +                  ' residual(twice) (mrad)',
     +       IM,BMEAS(IM,1)*1.D3,NG8,NNG8,RESID(IM),RESID(IM)
      ENDIF
   39 CONTINUE


      RETURN

      END
