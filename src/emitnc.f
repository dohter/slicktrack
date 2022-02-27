C   06/10/83 205151903  MEMBER NAME  EMITNC   (MAY92.S)  M  FORTRAN
      SUBROUTINE EMITNC(E0,IE0,CIR,CRAD,TDAMP,TUN,TM6A,COVMAT,SCALEMAT,
     +                                                  OFFSETL,OFFSETE)
C
C
C
C
C==================BEAM EMITTANCE CALCULATIONS==========================
C
C   ALSO RECALCULATE DAMPING USING THE M & R INTEGRATION METHOD.
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION TREV6(6,6),TREV66(6,6),VC6A(6),TDAMP(6),TUN(6)
      DIMENSION TM6A(6,6),TM6B(6,6),TM6C(6,6),TM6D(6,6),TM6E(6,6),
     +          TM6F(6,6),COVMAT(6,6),SCALEMAT(6,6)
      DIMENSION DAMMR(6),DAMMT(6),AA(6),AB(6),AAEXT(6)
      CHARACTER *1 TEXT(80)
CDP      REAL*8 IP/'IP'/
C
C
C
      INCLUDE "cloorb.for"
      INCLUDE "cnlist.for"
      INCLUDE "clatic.for"
      INCLUDE "csol.for"
      INCLUDE "cintpt.for"
      INCLUDE "cemit.for"
C
C
      WRITE(53,929)IE0
 929  FORMAT('1','Entering subroutine EMITNC, Energy step = ', I5)
C
C
      PI=3.1415926535897932D0
C
C
C      CALL UCOPY(TM6A,TREV66,72)
C      CALL UCOPY(TM6A,  TM6B,72)
      TREV66 = TM6A
      TM6B   = TM6A


C      CALL VZERO(DAMMR,12)
C      CALL VZERO(AA,12)
C      CALL VZERO(AAEXT,12)
       DAMMR = 0.D0
       AA    = 0.D0
       AAEXT = 0.D0

      CALL NORM(TM6A,6,AB)
C

      WRITE(53,1033)
      WRITE(53,1033)
 1033 FORMAT(/,'  ')
      WRITE(53,'(3A)')'Quantities from synchrotron eigenvector', 
     +' components to check the validity of the "oscillator model"' 


      DO 64 II=2,NELEM
      ITY=ITYPE(II)
      IID=ID(ITY)
      XY=XX(ITY)
      XX2=X2(ITY)
      IF(IID.EQ.5.OR.IID.EQ.6.OR.IID.EQ.7)GO TO 700
      IF((NAME(ITY)(1:2).EQ.'CQ'.AND.IID.EQ.3) !Ignore artificial lengths. 
     +  .OR.(NAME(ITY)(1:2).EQ.'RQ'.AND.IID.EQ.4).OR.IID.EQ.17)GO TO 700 
       S=S+YY(ITY)
  700 CONTINUE



      IF(IID.NE.1)GO TO 52
      DO 51 M=1,6
      TM6B(1,M)=TM6B(1,M)+XY*TM6B(2,M)
   51 TM6B(3,M)=TM6B(3,M)+XY*TM6B(4,M)
      GO TO 64
   52 IT=ITYPE(II)
      CALL MX66(TMAT(1,1,IT),IT,IID,II,XY,XX2,YY(IT),TREV6)
      CALL JAM666(TM6C,TREV6,TM6B)
      DELX=(DX(II)+DX(II+1))/2.
      DELY=(DY(II)+DY(II+1))/2.
      IF(IID.EQ.5)WRITE(53,'(A)')' Cavity!!!!'
      GOTO (600,602,603,604,611,602,602,608,602,610,600,600,600,600,
     +                                              602,602,600),IID
      GOTO 600
  602 CONTINUE
      XY1=DABS(XY)**3 /YY(ITY)**2
      GOTO 605
  603 CONTINUE
C     IF(NTWIST(ITY).EQ.3)DELY=DY(II)-TWIST(ITY)
C     IF(NTWIST(ITY).EQ.4)DELX=DX(II)-TWIST(ITY)
      XY1=0.D0
      IF(NAME(ITY)(1:1).NE.'E')              !Exclude edge fields.
     +                XY1=((XY*DELX)**2+(XY*DELY)**2)**1.5/YY(ITY)**2 

      GOTO 605
  604 CONTINUE
      XY1=((XY*DELX)**2+(XY*DELY)**2)**1.5/YY(ITY)**2
      GOTO 605
  608 XY1=((XY*(DELX**2-DELY**2)/2.)**2+(XY*DELX
     +                                      *DELY)**2)**1.5/YY(ITY)**2
      GOTO 605
CAs inDecember 1995:Turn off SOLXYP. i.e. ignore radiation in solenoids.
  610 CONTINUE
C 610 CALL SOLXYP(II,XY/YY(ITY), XY2,XY3)
C     XY1=((XY*XY2)**2+(XY*XY3)**2)**1.5/YY(ITY)**2
      XY1=0.D0

C====New, Aug 2008: TM6D was not defined at 600 for the first element.
C    So the (old) skip to 600 for cavities gave zero cavity damping if the 
C    cavity was the first element.
  611 CONTINUE 


  605 CONTINUE
      CALL AVER(TM6B,TM6C,6,TM6D)                    ! Need for DAMINT
      AA(1)=AA(1)+XY1*(TM6D(5,1)**2+TM6D(5,2)**2)
      AA(3)=AA(3)+XY1*(TM6D(5,3)**2+TM6D(5,4)**2)
      AA(5)=AA(5)+XY1*(TM6D(5,5)**2+TM6D(5,6)**2)

C      CALL AVERV5SQ(TM6B,TM6C,V5SQ1,V5SQ3,V5SQ5)
C      AA(1)=AA(1) + XY1*V5SQ1
C      AA(3)=AA(3) + XY1*V5SQ3
C      AA(5)=AA(5) + XY1*V5SQ5

C===Look at the value of the 5th component of the synch. eigenvector,
C   the value of the putative alpha and the phase evolution
C   to see how close we cme to the ``oscillator model''.
C   According to DESY86-147, the synch phase should only advance
C   where there is both horizontal curvature and dispersion.
      PHASE1 = DATAN2(TM6D(5,6),TM6D(5,5))/2.D0/pi
      PHASE2 = DATAN2(TM6D(6,6),TM6D(6,5))/2.D0/pi - 0.25D0
      BETAS  = 2.D0 * (TM6D(5,5)**2 + TM6D(5,6)**2)
      SIZES  = (TM6D(6,5)**2 + TM6D(6,6)**2)
      ALPHAS = 2.D0 * BETAS * ((TM6D(6,5)**2 + TM6D(6,6)**2)) - 1.D0
C      ALPHAS = DSQRT(ALPHAS)


C      IF(IID.EQ.2.OR.IID.EQ.15.OR.IID.EQ.3)
CCC      IF(NAME(ITY)(1:2).EQ.'EE'.OR.NAME(ITY)(1:2).EQ.'VD'.OR. 
CC     +   NAME(ITY)(1:2).EQ.'BP'.OR.NAME(ITY)(1:2).EQ.'BM')THEN
CCC     +                               NAME(ITY)(1:2).EQ.'BA')THEN
CC      WRITE(53,'(I10,2A,6E20.6)')II,' ',NAME(ITY),V5SQ1,V5SQ3,V5SQ5
CC     +   II,' ',NAME(ITY), TM6D(5,1),TM6D(5,2),TM6D(6,1),TM6D(6,2),         
CC     +                                 (TM6D(5,1)**2+TM6D(5,2)**2)
CC      WRITE(53,'(I10,2A,6E20.6)')
CC     +   II,' ',NAME(ITY), TM6D(5,3),TM6D(5,4),TM6D(6,3),TM6D(6,4),
CC     +                                 (TM6D(5,3)**2+TM6D(5,4)**2)
C     +WRITE(53,'(I6,F8.2,2A,2E11.3,2E12.4,A,2E11.3,3E12.4,A,F6.3)')
C     +   II,S,' ',NAME(ITY), TM6D(5,5),TM6D(5,6),BETAS,PHASE1,' |',
C     +            TM6D(6,5),TM6D(6,6),SIZES,ALPHAS,PHASE2,' |',
C     +           (TM6D(5,5)**2+TM6D(5,6)**2)*(TM6D(6,5)**2+TM6D(6,6)**2)
CCCC     +          II,'  ', NAME(ITY),(TM6D(5,5)**2+TM6D(5,6)**2)
CC      WRITE(53,1033)


CCC      ENDIF   
  600 CONTINUE
      CALL DAMINT(II,IID,ITY,CRAD,XX(ITY),X2(ITY),YY(ITY),TM6D,
     +                                           DAMMR,NAME(ITY))
C      CALL UCOPY(TM6C,TM6B,72)
      TM6B = TM6C
   64 CONTINUE
C
C
C      IF(1.EQ.1)STOP
C
C
C
C
C=====GET GENERALISED EMITTANCES: CAREFUL WITH VECTOR CONJUGATES:
C     E.G. DABS(AA) & DABS(AB)ARE NEEDED BECAUSE THE NORMALISATION AB
C     CAN (AS CONSTRUCTED) BE -VE BUT THE AA (AS CONSTRUCTED WITH OR
C     WITHOUT PRIOR NORMALISATION) ARE ALWAYS +VE.
      DO 702 J=1,5,2
      AA(J)=AA(J)*2.16D-19*(E0/.000511)**5/(2.*CIR*AB(J)*1.D3/TDAMP(J))
      AAN(J)=DABS(AA(J))
  702 AA(J)=DABS(AA(J))
      WRITE(53,101)(AA(KK),KK=1,5,2),AA(1)+AA(3),
     +                          AA(1)/TDAMP(1) + AA(3)/TDAMP(3)



C
C
      IF(AB(1).NE.0.D0)AA(1)=AA(1)*1.D6/2./DABS(AB(1))
      IF(AB(3).NE.0.D0)AA(3)=AA(3)*1.D6/2./DABS(AB(3))
      IF(AB(5).NE.0.D0)AA(5)=AA(5)*1.D6/2./DABS(AB(5))


      CALL SIG(AA(1),AA(3),AA(5),TREV66,TM6E)
      HEMIT=0.D0
      VEMIT=0.D0
      COUP =0.D0
      IF(BETAX0.NE.0.D0)HEMIT=TM6E(1,1)/BETAX0/1.D6
      IF(BETAY0.NE.0.D0)VEMIT=TM6E(3,3)/BETAY0/1.D6
      IF(HEMIT.NE.0.D0)COUP=VEMIT/HEMIT
      WRITE(53,104)HEMIT,VEMIT,COUP
C
C
C
C
C=====THE FOLLOWING LOGIC IS NOT NEEDED IF RENORM IS USED AND DATAN2
C=====FOR THE TUNES.
C=====Identify (hopefully) the vector associated predominantly with
C=====sync. oscillations: search for the smallest tune.
C     TUMIN=1.D10
C     LSYN=0
C     DO 700 J=1,5,2
C     IF(TUN(J).GT.TUMIN)GO TO 700
C     TUMIN=TUN(J)
C     LSYN=J
C 700 CONTINUE
C====='ENERGY EMITTANCE'
C     SYNCEM=AA(LSYN)/(3.D8)*E0*1.D9*2.*DABS(AB(LSYN))/1.D6
      SYNCEM=AA(5)/(3.D8)*E0*1.D9*2.*DABS(AB(5))/1.D6
      KSYN=(LSYN-1)/2+1
      KSYN=(5   -1)/2+1
      WRITE(53,102)KSYN,SYNCEM
C
      WRITE(53,98)
      WRITE(53,95)TM6E
      COVMAT   = TM6E
      SCALEMAT = TM6E


      WRITE(53,106)TM6E(1,1)+TM6E(3,3)
  106 FORMAT(' ',//,'      (1,1)       +       (3,3):       ', F13.5)
      WRITE(53,107)TM6E(1,1)/0.19/TDAMP(1) + TM6E(3,3)/0.27/TDAMP(3)
 107  FORMAT(' ',//,'dampx*(1,1)/betax + dampy*(3,3)/betay: ', F13.5)
C
C
C
C=====GET THE M & R DAMPING CONSTANTS FROM DAMMR
      DO 701 J=1,6
      IF(AB(J).NE.0.D0)THEN
      DAMMR(J)=DAMMR(J)/2./AB(J)
      DAMMT(J)=CIR/(3.D8*DAMMR(J))*1000.
      ENDIF
  701 CONTINUE

      WRITE(53,103)DAMMR,DAMMT
C
C
C
C
C
C
C
C
C
C====================CALCULATION OF BEAM SIZE==========================
C
C
C
C
C
C
C
C
C
      IF(IE0.NE.1)GO TO 69
      IF(ISIG.NE.0)WRITE(53,96)
      IF(ISIG.EQ.0)WRITE(53,99)
      NESIZE=1
      IF(ISIG.NE.0)NESIZE=NELEM
      S=0.D0
      DO 66 II=1,NESIZE
      IF(II.EQ.1) GO TO 67
      JTY=ITYPE(II-1)
      IID=ID(JTY)
      XY=XX(JTY)
      IF(IID.EQ.5.OR.IID.EQ.6.OR.IID.EQ.7)GO TO 801
      IF((NAME(JTY)(1:2).EQ.'CQ'.AND.IID.EQ.3) !Ignore artificial lengths. 
     +  .OR.(NAME(JTY)(1:2).EQ.'RQ'.AND.IID.EQ.4).OR.IID.EQ.17)GO TO 801 
       S=S+YY(JTY)
 801   CONTINUE



      CALL MX66(TMAT(1,1,JTY),JTY,IID,II-1,XX(JTY),X2(JTY),YY(JTY),TM6F)
      CALL JAM666(TREV66,TM6F,TREV66)
      CALL SIG(AA(1),AA(3),AA(5),TREV66,TM6E)
   67 TILT=DATAN2(2.D0*TM6E(1,3),TM6E(1,1)-TM6E(3,3))  ! 2 theta
      CTILT=DCOS(TILT)                              !cos (2 theta)

C=====Get major and minor axes of 1-sigma transverse ellipse in millimetres 
C     and their ratio. 
C      CTILT = 1.D0
      TILT=TILT*90.D0/PI
      XY3=TM6E(1,1)*TM6E(3,3)-TM6E(1,3)**2
      IF(XY3.LT.0..AND.XY3.GT.-1.D-8)XY3=0.
      AREA=PI*DSQRT(XY3)
      IF(CTILT.NE.0.D0.AND.XY3.NE.0.D0)THEN
      AXISMAJ   = DSQRT( 
     +XY3/( 0.5D0*((TM6E(3,3)-TM6E(1,1))/CTILT + TM6E(3,3)+TM6E(1,1))))
      AXISMIN   = DSQRT(
     +XY3/(-0.5D0*((TM6E(3,3)-TM6E(1,1))/CTILT - TM6E(3,3)-TM6E(1,1))))
      AREAROT = PI*AXISMAJ*AXISMIN
      ASPECTR   = AXISMIN/AXISMAJ
      SIGYOSIGX = TM6E(3,3)/TM6E(1,1)
      ENDIF
C======Now transform the sigmas 
      DLTA  = TILT*PI/180.D0 
      SIGXROT  = DSQRT(TM6E(1,1)*DCOS(DLTA)**2 + TM6E(3,3)*DSIN(DLTA)**2
     +            + 2.D0*DCOS(DLTA)*DSIN(DLTA)*TM6E(1,3))
      SIGYROT  = DSQRT(TM6E(1,1)*DSIN(DLTA)**2 + TM6E(3,3)*DCOS(DLTA)**2
     +            - 2.D0*DCOS(DLTA)*DSIN(DLTA)*TM6E(1,3))

C      This next line shows that SIGXROT=AXISMAJ,SIGYROT=AXISMIN,AREAROT=AREA
C      IF((ISIG.EQ.2.AND.NAME(ITYPE(II)).EQ.'IP').OR.II.EQ.1)
C     +WRITE(53,97)NAME(ITYPE(II)),
C     +           AXISMAJ,AXISMIN,SIGXROT,SIGYROT,SIGYROT/SIGXROT,AREAROT

      IF(AA(1).NE.0.D0)BETAXEFF = SIGXROT**2/AA(1)  !Betas if the dispersion is zero (using Chao coords).
      IF(AA(3).NE.0.D0)BETAYEFF = SIGYROT**2/AA(3)  
C                                          
      DO 68 IJ=1,6
   68 VC6A(IJ)=DSQRT(TM6E(IJ,IJ))
      IF((ISIG.EQ.2.AND.NAME(ITYPE(II)).EQ.'IP').OR.II.EQ.1)
     +WRITE(53,97)NAME(ITYPE(II)),
     +            VC6A,TILT,AREA,ASPECTR,SIGYOSIGX,BETAXEFF,BETAYEFF
      IF(ISIG.EQ.1.AND.II.LT.5000)
     +WRITE(53,97)NAME(ITYPE(II)),
     +            VC6A,TILT,AREA,ASPECTR,SIGYOSIGX,BETAXEFF,BETAYEFF

C======Plot the ellipse at the starting point. Take 1000 points. Also principle axes.
      IF(II.EQ.1)THEN
      DO 70 IPH = 1,1000
      PH    = 0.001D0*2.D0*PI*(IPH-1)
      HR    = AXISMAJ*DCOS(PH)    
      VR    = AXISMIN*DSIN(PH)    
      HOR   = HR*DCOS(DLTA) - VR*DSIN(DLTA)
      VER   = HR*DSIN(DLTA) + VR*DCOS(DLTA) 
      AXMAJ =    (-1.D0 + 0.002D0*(IPH-1))*DTAN(DLTA)
      AXMIN =    (-1.D0 + 0.002D0*(IPH-1))*DTAN(DLTA+PI/2.D0)
      WRITE(57,105)IPH,HOR,VER,-1.D0+0.002D0*(IPH-1),AXMAJ,AXMIN
  105 FORMAT(' ',I10, 6E16.6) 
   70 CONTINUE  
      ENDIF

   66 CONTINUE
C
   69 CONTINUE
C
C
C
C====Now insert the initial covariance matrix of an injected beam.    
C    We use the transverse emittances, the energy spread and the bunch length.
C    Since we have the transverse emittances, and the longitudinal parameters
C    are essentially independent of the azimuth, we can chose the injection point
C    as we wish. 
C    The input data are too sparse to allow the full matrix to be constructed. 

      COVMAT = SCALEMAT
      OFFSETL = 0.D0
      OFFSETE = 0.D0
      IF(IEXTSIZE.EQ.0)RETURN


      WRITE(53,1033)
      WRITE(53,1033)

      READ(5,'(80A1)')TEXT 
      WRITE(53,'(A)')' Externally supplied initial beam parameters:'
      READ(5,'(30A1,F10.4)') (TEXT(KT),KT=1,30),AAEXT(1)
      WRITE(53,'(A,30A1,F10.4)')' ',(TEXT(KT),KT=1,30),AAEXT(1)
      READ(5,'(30A1,F10.4)') (TEXT(KT),KT=1,30),AAEXT(3)
      WRITE(53,'(A,30A1,F10.4)')' ',(TEXT(KT),KT=1,30),AAEXT(3)
      READ(5,'(30A1,F10.4)') (TEXT(KT),KT=1,30),BLENGTH
      WRITE(53,'(A,30A1,F10.4)')' ',(TEXT(KT),KT=1,30),BLENGTH
      READ(5,'(30A1,F10.4)') (TEXT(KT),KT=1,30),ESPREAD
      WRITE(53,'(A,30A1,F10.4)')' ',(TEXT(KT),KT=1,30),ESPREAD
      READ(5,'(30A1,F10.4)') (TEXT(KT),KT=1,30),SIZEFAC
      WRITE(53,'(A,30A1,F10.4)')' ',(TEXT(KT),KT=1,30),SIZEFAC
      READ(5,'(30A1,F10.4)') (TEXT(KT),KT=1,30),OFFSETL
      WRITE(53,'(A,30A1,F10.4)')' ',(TEXT(KT),KT=1,30),OFFSETL
      READ(5,'(30A1,F10.4)') (TEXT(KT),KT=1,30),OFFSETE
      WRITE(53,'(A,30A1,F10.4)')' ',(TEXT(KT),KT=1,30),OFFSETE
      AA(5) = 0.D0
      GAMMA =  E0/.000511
C====Get real emittances 
      AAEXT(1) = AAEXT(1)/GAMMA                          
      AAEXT(3) = AAEXT(3)/GAMMA
      AAEXT(5) = 0.D0
C====Get the COVMAT in mm^2,mrad^2 and mm.mrad for comparison with the original. 
      CALL SIG(AAEXT(1),AAEXT(3),AAEXT(5),TREV66,TM6E)
      COVMAT = 0.D0
      COVMAT(1:4,1:4) = TM6E(1:4,1:4) * SIZEFAC**2

C    Put in the ``external values for (5,5) and (6,6).
      COVMAT(5,5) = BLENGTH**2 * SIZEFAC**2
      COVMAT(6,6) = (ESPREAD/(E0*1000.D0))**2 * 1.D6 * SIZEFAC**2
      OFFSETE     =  OFFSETE/(E0*1000.D0)     * 1.D3
      WRITE(53,981)
      WRITE(53,95)COVMAT
C
C
C

C      IF(1.EQ.1)STOP
      RETURN
C
   95 FORMAT(T6,6F13.5)
   96 FORMAT(//,' BEAM SIZE SIGMAS AROUND THE RING:',//,T5,'NAME',
     +    T14,'SIGX(mm)',T24,'SIGXP(MRAD)',T36,'SIGY(mm)',T46,'SIGYP',
     +    '(MRAD)',T58,'SIGZ(mm)',T68,'SIGE/E*1000',T80,'TILT(DEG)',
     +    T90,'AREA(mm*mm)',T103,'MIN/MAJ',T113,'SIGY/SIGX',
     +    T124, 'BETAXEFF',T135,'BETAYEFF',/)
   97 FORMAT(T5,A8,T11,12F11.5)
   98 FORMAT(//,' BEAM COVARIANCE MATRIX <Xi*Xj>(mm*mm,mm*mrad...)',
     +          ' AT THE 1-ST BEAM-LINE ELEMENT FOR A STABLE BEAM:')
  981 FORMAT(//,' BEAM COVARIANCE MATRIX <Xi*Xj>(mm*mm,mm*mrad...)',
     +          ' AT THE 1-ST BEAM-LINE ELEMENT FOR THE INJECTED BEAM:')
   99 FORMAT(//,' BEAM SIZE SIGMAS IN 1. ELEMENT : ',//,T5,'NAME',
     +    T14,'SIGX(mm)',T24,'SIGXP(MRAD)',T36,'SIGY(mm)',T46,'SIGYP',
     +    '(MRAD)',T58,'SIGZ(mm)',T68,'SIGE/E*1000',T80,'TILT(DEG)',
     +    T90,'AREA(mm*mm)',T103,'MIN/MAJ',T113,'SIGY/SIGX',
     +    T124, 'BETAXEFF',T135,'BETAYEFF',/)
  101 FORMAT(' ',///,' GENERALISED EMITTANCES (METRE-RAD) ',3E12.3,
     +               ',   Sum of first two', E12.3,
     +',   Damping weighted sum of first two', E12.3)
  102 FORMAT(' ',///,' THE SYNC. VECTOR IS NUMBER ',I3,/,
     +               ' THE SYNC. PHASE SPACE IS   ',E12.3,' EV.SECS')
  103 FORMAT(///,' DAMPING CALCULATED FROM M & R INTEGRALS:',/,
     +    T3,'DAMPING CONST.', T18,6(1X,F11.5),/,
     +    T3,'DAMP.TIME(MSEC)',T18,6(1X,F11.5))
  104 FORMAT(' ',  /,' "BETATRON"  EMITTANCES (METRE RAD) ',2E12.3,/
     +              ,' "COUPLING"                         ', E12.3)
C
C
      END
