C   02/09/86 609111554  MEMBER NAME  ORCORR   (S)           FORTRAN
      SUBROUTINE ORCORR(BX,BZ,LAG,KXZ,KICK1,KICK2,WIN,WPU,KON,K32,Q)
C
C=====ROUTINE FOR DESY22 CONTROL SYSTEM TO MAKE CORRECTIONS TO THE
C=====ORBIT.
C=====THE METHOD EMPLOYS A ROW OF CLOSED BUMPS TO MAKE CORRECTIONS
C=====TO MONITORS BETWEEN KICKERS KICK1 & KICK2 IN PLANE KXZ.
C=====NCOR=(KICK2-KICK1) MUST GE.2.
C
C
      INTEGER C
      REAL LAG
      DIMENSION BX(4,24),BZ(4,24),B(4,24),C(4,24),LAG(2,24)
      DIMENSION WIN(2,24),ALAG(2,24),WPU(2,24),KON(2,24),Q(2),Y(24)
C
C
C=====COPY THE TWISS PARAMETERS INTO THE WORKING SPACE B
      DO 10 J=1,24
      DO 10 K=1,4
      B(K,J)=BX(K,J)
      IF(KXZ.EQ.1)B(K,J)=BZ(K,J)
   10 CONTINUE
C
C=====COPY MEASURED ORBIT INTO WORKING SPACE,ZERO THE KICKER ANGLE
C=====UPDATE ARRAY.CALCULATE AVERAGES
      AVELAG=0.
      RMSLAG=0.
      DO 11 J=1,24
      WIN(KXZ,J)=0.
      ALAG(KXZ,J)=LAG(KXZ,J)
      AVELAG=AVELAG+LAG(KXZ,J)
      RMSLAG=RMSLAG+LAG(KXZ,J)*LAG(KXZ,J)
   11 CONTINUE
      AVELAG=AVELAG/24.
      RMSLAG=SQRT(RMSLAG/24.)
      WRITE(6,21)AVELAG,RMSLAG
   21 FORMAT(' ','AVERAGE ORBIT ',F10.4,' MM',' RMS ORBIT ',F10.4,' MM')
C
C
C
C=====LOOP OVER THE KICKERS MAKING A SUCCESSION OF CLOSED BUMP
C=====CORRECTIONS USING THREE KICKERS AFTER SETTING UP THE PHASE
C=====DIFFERENCES ETC.
C=====USE A SINGLE PARAMETER,THE STRENGTH OF THE FIRST KICKER,TO
C=====MINIMISE THE SUM SQ. DEVIATION OF THE TWO MONITORS IN BETWEEN.
C=====IF SOME KICKERS ARE DEAD THEN THE SPAN IS LARGER AND MORE MONITORS
C=====ARE INCLUDED IN THE MINIMISATION.
C=====ITERATE 3 TIMES TO SUCCESSIVELY IMPROVE THE RESULT.
C=====IF WHOLE RING WAS REQUESTED THEN MAKE AN OVERLAP AT THE END.
      K32=0
      NCOR=KICK2-KICK1
      IF(NCOR.LE.0)NCOR=NCOR+24
      IF(NCOR.EQ.24)NCOR=26
      NCOR=NCOR-1
C
C=====ITERATION LOOP
      DO 3 JJ=1,3
C
C=====RUN THE BUMP ACROSS THE MONITORS.CODE ALLOWS TO CIRCULATE
C=====CONTINUOUSLY AROUND THE MACHINE.
      DO 1 L=1,NCOR
      J1=KICK1+L-1
      IF(J1.GT.24)J1=J1-24
      IF(KON(KXZ,J1).EQ.0)GO TO 1
      J2=J1+1
      JJ2=J2
      J2SAVE=J2
  102 IF(J2.GT.24)J2=J2-24
      IF(KON(KXZ,J2).NE.0)GO TO 101
      J2=JJ2+1
      JJ2=J2
      J2SAVE=J2
      GO TO 102
  101 J3=J2+1
      JJ3=J3
  104 IF(J3.GT.24)J3=J3-24
      IF(KON(KXZ,J3).NE.0)GO TO 103
      J3=JJ3+1
      JJ3=J3
      GO TO 104
C=====NEED A TRAP TO STOP IT LOOPING FOR EVER.
C
  103 CONTINUE
C     WRITE(6,61)ALAG(KXZ,J2),ALAG(KXZ,J3),J1,J2,J3,Q(KXZ)
   61 FORMAT(' ','STARTING ORBIT',2F10.4,3I6,F10.4)
C
C=====GET PHASES SHIFTS BETWEEN THE KICKERS.
      DK31=B(2,J3)-B(2,J1)
      IF(DK31.LT.0.)DK31=DK31+Q(KXZ)
      DK31S=DK31
      DK31=SIN(DK31)
      DK21=B(2,J2)-B(2,J1)
      IF(DK21.LT.0.)DK21=DK21+Q(KXZ)
      DK21S=DK21
      DK21=SIN(DK21)
      DK32=B(2,J3)-B(2,J2)
      IF(DK32.LT.0.)DK32=DK32+Q(KXZ)
      DK32S=DK32
      DK32=SIN(DK32)
      IF(DK32.NE.0.0)GO TO 2
      K32=1
      RETURN
    2 R21=-B(1,J1)/B(1,J2)*DK31/DK32
      R31= B(1,J1)/B(1,J3)*DK21/DK32
C
C=====LOOP OVER MONITORS BETWEEN THE WORKING KICKERS.
C=====GET PHASE SHIFTS BETWEEN THE KICKERS AND MONITORS.
      UP=0.
      DOWN=0.
      DO 5 II=1,24
      II1=J1+II
      IISAVE=II1
      IF(II1.GT.24)II1=II1-24
      IF(IISAVE.GT.J2SAVE)GO TO 8
      DPU1=B(4,II1)-B(2,J1)
      IF(DPU1.LT.0.)DPU1=DPU1+Q(KXZ)
      DPU1=SIN(DPU1)
      Y(II)=B(1,J1)*B(3,II1)*DPU1
      UP=UP-Y(II)*ALAG(KXZ,II1)*WPU(KXZ,II1)
      DOWN=DOWN+Y(II)*Y(II)*WPU(KXZ,II1)
      GO TO 5
    8 DPU2=B(4,II1)-B(2,J1)
      IF(DPU2.LT.0.)DPU2=DPU2+Q(KXZ)
      DPU2=SIN(DPU2)
      Y(II)=B(1,J1)*B(3,II1)*DPU2
      DPU3=B(4,II1)-B(2,J2)
      IF(DPU3.LT.0.)DPU3=DPU3+Q(KXZ)
      DPU3=SIN(DPU3)
      Y(II)=Y(II)+B(1,J2)*B(3,II1)*DPU3*R21
C     IF(WPU(KXZ,J2)+WPU(KXZ,J3).LT.0.001)GO TO 1
      UP=UP-Y(II)*ALAG(KXZ,II1)*WPU(KXZ,II1)
      DOWN=DOWN+Y(II)*Y(II)*WPU(KXZ,II1)
      IF(II1.EQ.J3)GO TO 6
C     A1=-ALAG(KXZ,J2)/Y(II)
C     A1=0.003
    5 CONTINUE
    6 CONTINUE
      A1=UP/DOWN
      A2=A1*R21
      A3=A1*R31
C=====ACCUMULATED CORRECTION ANGLES.
C     WRITE(6,72)WIN(KXZ,J1),WIN(KXZ,J2),WIN(KXZ,J3),L,NCOR
      WIN(KXZ,J1)=WIN(KXZ,J1)+A1
      WIN(KXZ,J2)=WIN(KXZ,J2)+A2
      WIN(KXZ,J3)=WIN(KXZ,J3)+A3
C     WRITE(6,72)WIN(KXZ,J1),WIN(KXZ,J2),WIN(KXZ,J3),L,NCOR
  72  FORMAT(' ','WIN ',3F10.4,2I4)
C=====ACCUMULATED PREDICTED NEW ORBIT.
      DO 55 II=1,24
      II1=J1+II
      IF(II1.GT.24)II1=II1-24
      ALAG(KXZ,II1)=ALAG(KXZ,II1)+Y(II)*A1
      IF(II1.EQ.J3)GO TO 56
   55 CONTINUE
   56 CONTINUE
C     WRITE(6,51)ALAG(KXZ,J2),ALAG(KXZ,J3),J1,J2,J3
   51 FORMAT(' ','ORCORR ALAG   ',2F10.4,3I6)
C     WRITE(6,52)DK31S,DK21S,DK32S,DPU1S,DPU2S,DPU3S
   52 FORMAT(' ','PHASES        ',6F10.4)
    1 CONTINUE
    3 CONTINUE
C
C=====GET AVERAGE + RMS OF CORRECTED ORBIT
      AVELAG=0.
      RMSLAG=0.
      DO 41 J=1,24
      AVELAG=AVELAG+ALAG(KXZ,J)
      RMSLAG=RMSLAG+ALAG(KXZ,J)*ALAG(KXZ,J)
   41 CONTINUE
      AVELAG=AVELAG/24.
      RMSLAG=SQRT(RMSLAG/24.)
      WRITE(6,21)AVELAG,RMSLAG
C
C
      RETURN
      END
