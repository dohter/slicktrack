      SUBROUTINE FIXORB5(E0,IE0,U0,CRAD,FREQ,KKICK)


C       A short routine to check total kicks given by FIXORB1,2,3,4.
C       Keep it as close to FIXORBN as possible to avoid confusions.
C       Note that due to the instability of the kick minimisation algorithm
C       the combined kicks here might not well reproduce those of the the kick minimisation
C       routine.

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C      PARAMETER (NM = 620, NC= 620)         !Number of measurements and correcters
      PARAMETER (NM = 977, NC= 977)
      DIMENSION TREV(7,7),TEMP(7,7),TM7(7,7),VC6A(6),VC6B(6),TM6(6,6)
C======Storage for the least squares NAG routine.
      DIMENSION BMEAS(2*NM,1),BQUAD(2*NM,1),XXSAVE(NC) 
      DIMENSION ORBMEAS(NM),ORBREAL(NM)

      DIMENSION   NVCPOS(NC)
      CHARACTER*8 NVNAME(NC)
      DIMENSION NBIND(NM), IBIND(NM,10)
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
      WRITE(53,'(A,A,A,A)')
     +   '1',' FIXORB5: Construct the closed orbit to check the',
     + ' combined kicks of correctors and their adjacent quadrupoles.'

C======First loop over all elements to set up the connection between vertical correctors
C      and nearby relevant magnets.
C      This algorithm depends on details of the layout. So it must be checked for each 
C      new optic/layout.

      NVCCNT = 0
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
      WRITE(53,'(A,I4,A)')' Monitor ',NVCCNT,
     +                ' Too few bound elements: So STOP.'
      STOP
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
      KG8 = 0

      NVCCNT = 0

      NA = 0
      NDO = 0

C======FIND THE CLOSED ORBIT AND THE KICKS: 
C      Switch off the iteration so that sextupoles are ignored
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

      IF(ITORB.NE.1)THEN                   ! On first iteration, skip nonlinearities.
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
      WRITE(53,'(A)')' Problem getting C.O. in FIXORB5, so STOP'
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
      RMSKICK = 0.D0
      RMSQUAD = 0.D0

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

      
C======Get the combined quad-corrector kicks for this VCO.
C      NVM labels the measurements. NVCCNT labels the perturbing kicks.

C======Give the ``measurements'' a shift error (G6) and a scale error (G7).
C      The orbit recorded by a monitor is: (real orbit - offset)*(1 + scale error).  
C      We want to find corrector settings which give a C.O. and correctors for
C      which the existing combined kicks are cancelled.

C======DY(II) is the V.C.O. just before the element II.
      IF(ID(ITYPE(II)).EQ.7.AND.NAME(ITYPE(II))(1:2).EQ.'VC')THEN
      NVM=NVM + 1 

      IF(NA.EQ.0.AND.NVM.LE.1)write(53,1033)
      IF(NA.EQ.0.AND.NVM.EQ.1)WRITE(53,'(A)')' Specimen kick data' 
      IF(NA.EQ.0.AND.NVM.LE.10)write(53,1033)
      IF(NA.EQ.0.AND.NVM.LE.10)
     +                write(53,'(A,A,3I10)')' Corrector found    :  ',
     +                      NAME(ITYPE(II)),ID(ITYPE(II)),II,NVM

      IF(NA.EQ.0)THEN
      KG8 = KG8 + 1
C======Now get the combined kick of a quad and the nearby kicker VC.
C      The shift of the quad is already embodied in the strength of the VQ.
C      but include it in FIXORB5 for checking the total kick at each quad.
C      Use TSPIN conventons for the signs of quad and corrector kickes.
C      Get the nearby quads:
      ORBMEAS(NVM) = (DY(II) - G8(KG8)*0.001D0)
      ORBREAL(NVM) =  DY(II)
      BMEAS(NVM,1) = 0.D0
      BQUAD(NVM,1) = 0.D0
      QCOMB = 0.D0
      ICNT=0
      XXVQ = 0.D0
      DO 66 JN = 1,NBIND(NVM)
      IK  =  IBIND(NVM,JN)
      KTY =  ITYPE(IK) 
      IKD =  ID(KTY)
      IF(NVM.LE.10)
     +      write(53,'(A,I10,3X,A,3I10)')' Bound element:  ', II,
     +                                        NAME(KTY),IKD,JN,IK
      IF(IKD.EQ.3)THEN
C      QCOMB= QCOMB + XX(KTY)
      QCOMB = QCOMB + TMAT(4,3,KTY)
      IF(NVM.LE.10)write(53,'(A,F16.5)')'Strength: ',XX(KTY)
      ICNT =ICNT + 1
      ENDIF
      IF(IKD.EQ.15)THEN
C      QCOMB= QCOMB + X2(KTY)
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
     +                ' Number of bound elements:',ICNT,' So STOP.'
      STOP
      ENDIF
      BMEAS(NVM,1) = QCOMB*ORBMEAS(NVM) - XX(ITYPE(II)) - XXVQ
      BQUAD(NVM,1) = QCOMB*ORBREAL(NVM)
      IF(NVM.LE.10)WRITE(53,'(A,I4,A,4F11.6,3E15.6)')' Monitor ',NVM,
     + 'Various kick strengths (mrad) ',ORBMEAS(NVM)*1.D3,
     +                                  QCOMB*ORBMEAS(NVM)*1.D3,
     +                                  XX(ITYPE(II))*1.D3,
     +                                  XXVQ*1.D3,
     +                                  BMEAS(NVM,1)*1.D3,
     +                                  BQUAD(NVM,1)*1.D3

      RMSKICK = RMSKICK + BMEAS(NVM,1)**2
      RMSQUAD = RMSQUAD + BQUAD(NVM,1)**2
      NVMKEEP = NVM
      ENDIF
   
      ENDIF
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
     +     'Corresponding rms combined kicks (mrad): ',RMSKICK*1.D3
      WRITE(53,'(A,A,1X,F15.6,A)')' ',
     +     'rms quad kicks                   (mrad): ',RMSQUAD*1.D3

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


      RETURN

      END


