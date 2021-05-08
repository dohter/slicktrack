
      P R O G R A M SPLMAG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ PROGRAM    :  SPLMAG
C@ AUTHORS    :  M.Vogt (vogtm@mail.desy.de) / CODE
C-               D.Barber (mpybar.desy.de)   / SLIKCTRACK FORMAT
C@ VERSION    :  0.00.01
C@ RELEASE    :  APR-22-2005                             ! ==>> 'ctrl-s (***)'
C@ REMARKS    :  prog for slicing quads and bends in slim file
C-            :  also continuesns input reader and writer for other
C-            :  lattice-lanuages
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)


C     Variables:
C
C     MBITS=31 can be MBITS=32 once we trust signed integers
C
C     set JBNOSP : element is a ``spin-drift''
C      "  JBSNEG : scale negative
C      "  JBSPC7 : special (type-dependent) flag #7
C      "  JBSPC8 : special (type-dependent) flag #8
C      "  JBSPC9 : special (type-dependent) flag #9

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MAXVNR=99)

C%%      write(*,*)'SPLMAG'

      READ(*,*) CJOB
      READ(*,*) CLATI
      READ(*,*) CLATO
      READ(*,*) NUMOPT
      IF(NUMOPT.GT.10) CALL TOFF('SPLMAG: too many options')
      DO 1000 I=1,NUMOPT
         READ(*,*) COPT(I) 
 1000 CONTINUE
         


      IE = ICSEND(CJOB)
      IF(IE.GT.MAXNME-8) CALL TOFF('SPLMAG: jobname too long')
      CJOB(IE+1:IE+5) = '.log'
      CALL NEWVNR(CJOB,MAXVNR,CJOBV,IVNR,LOUT,IERR,*9901)
      OPEN(16, FILE=CJOBV(1:LOUT), ACCESS='SEQUENTIAL', STATUS='NEW',
     +     FORM='FORMATTED',ERR=9902)
      CJOB(IE+1:IE+5) = '    '

      OPEN(15, FILE=CLATI(1:ICSEND(CLATI)), ACCESS='SEQUENTIAL',
     +     STATUS='OLD', FORM='FORMATTED',ERR=9903)

      IE = ICSEND(CLATO)
      IF(IE.GT.MAXNME-4) CALL TOFF('SPLMAG: outfilename too long')
      CALL NEWVNR(CLATO,MAXVNR,CLATOV,IVNR,LOUT,IERR,*9992)
      OPEN(17, FILE=CLATOV(1:LOUT), ACCESS='SEQUENTIAL', STATUS='NEW',
     +     FORM='FORMATTED',ERR=9992)

      WRITE(16,'(3A)') 'LOGFILE    >>', CJOB (:ICSEND(CJOB )), '<<'
      WRITE(16,'(3A)') 'INPUTFILE  >>', CLATI(:ICSEND(CLATI)), '<<'
      WRITE(16,'(3A)') 'OUTPUTFILE >>', CLATO(:ICSEND(CLATO)), '<<'
      WRITE(16,'(10(A,I2,3A/))')
     +     ( 'OPT#',I,' >>',COPT(I)(:ICSEND(COPT(I))),'<<', I=1,NUMOPT )

      CALL LATIN

      CALL GENELE

      CALL LATOUT


      STOP
C=====error handling

 9901 CONTINUE
      IF    (IERR.EQ.1) THEN
         CALL TOFF('SPLMAG: jobname too long')
      ELSEIF(IERR.EQ.2) THEN
         CALL TOFF('SPLMAG: too many verisons of logfile')
      ENDIF
 9902 CONTINUE
      CALL TOFF('SPLMAG: error opening logfile')
 9903 CONTINUE
      CALL TOFF('SPLMAG: error opening inpfile')
 9991 CONTINUE
      IF    (IERR.EQ.1) THEN
         CALL TOFF('SPLMAG: outfilename too long')
      ELSEIF(IERR.EQ.2) THEN
         CALL TOFF('SPLMAG: too many verisons of outfile')
      ENDIF
 9992 CONTINUE
      CALL TOFF('SPLMAG: error opening outfile')

      E N D


      S U B R O U T I N E GENELE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  GENELE
C@ AUTHORS    :  M.Vogt
C@ COMMAND    :  
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  16 -> logfile
C@ DESCRIPTION:  insert special kicks into certain elements e.g. for use
C-               with SLICKTRACK + add thin quads for entr/exit pole faces
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C>>   [1b,3b,4b,6b] added to cope with edge fields / 08.10.2005 / vogtm
C
C     0: remove END in typelist and lattice list
C     1: iterate original typelist :
C       - checks if NAMES of pattern HQ*, VQ* ,RQ*, CQ*, VD* exist
C         => if TRUE then CALL TOFF(....)
C     1b: doe same with paterns EI*, EO*, EE* (EdgeIn, EdgeOut, Edge[in=out])
C     2: check corrector coil & monitor doubles :
C       - iterate through original typelist :
C       -- for each corrector coil : iterate through original lattice :
C       --- if double : 
C       ----- add # to name in typelist ; #++ 
C       ----- append name//# to typelist
C     3: iterate original typelist :
C       - splits every ITYPE in
C         {IBX, IMQ, ISQ, IBY, IRQ, IRX, IRY, ICX, ICY, ICH, ICV} into 2 halfs
C         -- thereby NAME -> NAME //'H' & BSPLIT = .TRUE.
C     3b: for every ITYPE in {IBX, IBY, IRX, IRY, ICX, ICY, ICH, ICV} check if
C         -- one of E1, E2, HGAP, FINT .NE. 0. or ITYPE in {IIRX, IRY, ICH, ICV}
C            (RBEND type) ; then
C            --- E1X = E1 + ANGLE/2 if RBEND, = E1 otherwise (same with E2X) 
C            --- change RBEND type to SBEND type
C            --- if E1X.EQ.E1X then create single thin quad EE_# and
C                    IEDGE(1,elenumber) = NELE ; IEDGE(2,elenumber)=NELE
C            ---  else create EI_* and EO_* in type list and
C                    IEDGE(1,elenumber) = NELE-1 ; IEDGE(2,elenumber)=NELE
C     4: iterate original lattice list:
C       - count ITYPE in {IBX, IBY, IRX, IRY} (plain bends) -> NBEND
C       - count ITYPE in {IMQ, ISQ, IRQ} (plain quads)      -> NQUAD
C       - count ITYPE in {ICX, ICY, ICH, ICV} (comb.func.)  -> NCMBF
C       - NEWELE = NBEND + 4*NQUAD + 4*NCMBF must be .LT. MAXELE-NELE
C         => if FALSE then CALL TOFF(...)
C     4b:(note) athough RBEND types should not exist, keep testing !
C     5: copy original lattice list to temp lattice list
C        (original will be overwritten )
C     6: iterate temp lattice list:
C        - if .NOT.BSPLIT => copy to lattice list
C          else 
C          -- if PLAIN BEND :
C             --- append VD_# to typelist
C   (6b)      --- if IEDGE(1,LATELE) .EQ. 0 then 
C                ---- replace BEND in lattice list by HalfBend VD HalfBEND'
C   (6b)      --- else (works for EdgeIn.eq.EdgeOut or not !!!)
C                ---- replace by EdgeIn HalfBend VD HalfBEND EdgeOut
C             --- increment # of VD
C          -- if PLAIN QUAD :
C             --- append HQ_#, VQ_#, RQ_#, CQ_#  to typelist
C             --- replace QUAD in lattice list by HalfQUAD HQ VQ RQ CQ halfQuad
C             --- increment # of HQ, VQ, RQ, CQ
C          -- if COMBINED FUNCTION :
C             --- append HQ_#, VD_#, RQ_#, CQ_#  to typelist
C   (6b)      --- if IEDGE(1,LATELE) .EQ. 0 then 
C              ---- replace CMBF in lattice list by HalfCMBF HQ VD RQ CQ HalfCMBF
C   (6b)      --- else (works for EdgeIn.eq.EdgeOut or not !!!)
C              ---- replace CMBF by EdgeIn HalfCMBF HQ VD RQ CQ HalfCMBF EdgeOut
C             --- increment # of HQ, VD, RQ, CQ
C     7: add END in typelist and lattice list


      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MAXNUM=MAELNM-4)
C      PARAMETER(DLRAD=.1D+0)
      DIMENSION BSPLIT(MAXELE), IEDGE(2,MAXELE)
      CHARACTER APPOVWR*(MAELNM), APPEND*(MAELNM),
     +     CFINFL*(MAXNUM), CDUM*(MAXNUM)
      DIMENSION LATEL0(MAXLIST),PREDR0(MAXLIST),LSLI0(MAXLIST) ! loc.cp /LATTICE/

C=====statementfunctions ...
C=====... for dealing with edge effects
      BSTRN(I) = STRN(1,I).NE.0.D+0 .OR. STRN(2,I).NE.0.D+0 
     +     .OR. STRN(3,I).NE.0.D+0 .OR. STRN(4,I).NE.0.D+0  ! edge treatment
      BFRIN(I) = STRN(3,I).NE.0.D+0 .OR. STRN(4,I).NE.0.D+0 ! fringe corr: x.ne.y
      BSAME(I) = STRN(1,I) .EQ. STRN(2,I)                   ! EdgeIn=EdgeOut
      DKEX(E,ANGLE,DL) = -ANGLE*DTAN(E)/DL
      DKEY(E,ANGLE,DL,HGAP,FINT) =
     +     -ANGLE*DTAN( E - ANGLE*HGAP*FINT/DL * (1.D+0+DSIN(E)**2) )/DL

C%%      write(*,*)'genele'

      CALL INFO('GENELE')

C     0: remove END in typelist and lattice list
C>>>>>IF( CNAME(NELE).EQ.'X______EOL' ) THEN <<< changed / vogtm 23.04.09
      IF( CNAME(NELE).EQ.'X______EOL' .OR. CNAME(NELE).EQ.'END') THEN
         NELE = NELE -1
      ELSE
         CALL TOFF('GENELE, last element in typelist must be END '//
     +        'or X______EOL')
      ENDIF
C>>>>>IF( CNAME(LATELE(NLAT)).EQ.'X______EOL' ) THEN <<<changed/vogtm 23.04.09
      IF( CNAME(LATELE(NLAT)).EQ.'X______EOL' .OR.
     +    CNAME(LATELE(NLAT)).EQ.'END' ) THEN
         PRESAVE = PREDRF(NLAT)
         NLAT = NLAT -1
      ELSE
         CALL TOFF('GENELE, last element in lattice must be END '//
     +        'or X______EOL')
      ENDIF

C     1: iterate original typelist :
C       - checks if NAMES of pattern HQ*, VQ* ,RQ*, CQ*, VD* exist
C         => if TRUE then CALL TOFF(....)
C     1b: doe same with paterns EI*, EO*, EE* (EdgeIn, EdgeOut, Edge[in=out])

      NBADNME = 0

      DO 1000 I=1,NELE
         IF( CNAME(I)(1:2).EQ.'HQ' .OR. CNAME(I)(1:2).EQ.'VQ'
     +        .OR. CNAME(I)(1:2).EQ.'RQ' .OR. CNAME(I)(1:2).EQ.'CQ'
     +        .OR. CNAME(I)(1:2).EQ.'VD' .OR. CNAME(I)(1:2).EQ.'EI'
     +        .OR. CNAME(I)(1:2).EQ.'EO' .OR. CNAME(I)(1:2).EQ.'EE'
     +        ) THEN
            NBADNME = NBADNME + 1
            WRITE(26,'(1X,3A,A,A,I10)')
     +           'invalid (slicktrack)-name : ', CNAME(I), ' # ', I
         ENDIF
 1000 CONTINUE

      IF( NBADNME.GT.0 ) CALL TOFF(
     +     'GENELE : invalid (slicktrack)-names encountered')

C     2: check corrector coil and monitor doubles :
C       - iterate through original typelist :
C       -- for each corrector coil : iterate through original lattice :
C       --- if double : 
C       ----- add # to name in typelist ; #++ 
C       ----- append name//# to typelist

      NENEW = NELE

      DO 2000 I=1,NELE
         IF( ITYPE(I).EQ.IKX .OR. ITYPE(I).EQ.IKY
     +        .OR. ITYPE(I).EQ.IMO ) THEN ! corrector/monitor
            NC = 0
            DO 2100 J=1,NLAT
               IF( LATELE(J).EQ.I ) NC = NC + 1
 2100       CONTINUE
            IF( NC.LE.1 ) GOTO 2000 ! zero or one found => unique
            WRITE(16,'(1X,A,A,A,I5,A)')
     +           'GENELE : corrector ', CNAME(I), ' found with ',
     +           NC, ' copies'
            IE = ICSEND( CNAME(I) )
            LF = MAELNM - IE    ! free spaces in namefield
            IF( LF.EQ.0 .OR. (LF.EQ.1 .AND. NC.GT.9) .OR.
     +           (LF.EQ.2 .AND. NC.GT.99  ) .OR.
     +           (LF.EQ.3 .AND. NC.GT.999 ) .OR.
     +           (LF.EQ.4 .AND. NC.GT.9999) ) THEN
               CALL WARN('GENELE, corrector '//CNAME(I)//
     +              ' cannot be made unique')
               GOTO 2000
            ENDIF
            LF   = MIN0(MAXNUM,LF)
            IC   = 1
            CDUM = CFINFL(IC,'0')
            CNAME(I)(IE+1:MAELNM) = CDUM(MAXNUM-LF+1:MAXNUM) 
            DO 2200 J=1,NLAT
               IF( LATELE(J).EQ.I ) THEN
                  J1 = J+1
                  GOTO 2201     ! first occurrence does not get changed
               ENDIF
 2200       CONTINUE
 2201       CONTINUE
            DO 2300 J=J1,NLAT   ! start after first occurrence
               IF( LATELE(J).EQ.I ) THEN
                  IC = IC + 1
                  CDUM = CFINFL(IC,'0')
                  NENEW = NENEW +1 
C                 if( nele+nenew .gt. maxele )  call toff...! WHY NELE+NENEW ???
                  IF( NENEW .GT. MAXELE )  CALL TOFF(
     +                 'GENELE, increase MAXELE')
                  ITYPE(NENEW) = ITYPE(I)
                  STRX (NENEW) = STRX (I)
                  RLEN (NENEW) = RLEN (I)
                  CNAME(NENEW) =CNAME(I)(1:IE)//CDUM(MAXNUM-LF+1:MAXNUM)
                  LATELE(J) = NENEW
               ENDIF
 2300       CONTINUE
         ENDIF
 2000 CONTINUE

      WRITE(16,'(1X,A,I10,A,I10,A)')
     +     'GENELE : typelist increased from ', NELE,
     +     ' to ', NENEW, ' after making correctors unique'  
      NELE = NENEW

C     3: iterate original typelist :
C       - splits every ITYPE in
C         {IBX, IMQ, ISQ, IBY, IRQ, IRX, IRY, ICX, ICY, ICH, ICV} into 2 halfs
C         -- thereby NAME -> NAME //'H' & BSPLIT = .TRUE.
C     3b: for every ITYPE in {IBX, IBY, IRX, IRY, ICX, ICY, ICH, ICV} check if
C         -- one of E1, E2, HGAP, FINT .NE. 0. or ITYPE in {IIRX, IRY, ICH, ICV}
C            (RBEND type) ; then
C            --- E1X = E1 + ANGLE/2 if RBEND, = E1 otherwise (same with E2X) 
C            --- change RBEND type to SBEND type
C            --- if E1X.EQ.E1X then create single thin quad EE_# and
C                    IEDGE(1,elenumber) = NELE ; IEDGE(2,elenumber)=NELE
C            ---  else create EI_* and EO_* in type list and
C                    IEDGE(1,elenumber) = NELE-1 ; IEDGE(2,elenumber)=NELE

      NSPLIT = 0
      NENEX  = 0
      NSAME  = 0


      DO 3000 I=1,NELE          ! as block-if to be prepared for specialities

         IF    ( ITYPE(I) .EQ. IBX) THEN
            NSPLIT = NSPLIT + 1
            BSPLIT(I) = .TRUE.
            IF( BSTRN(I) ) CALL NEWEDGE(
     +           I,NENEW,NENEX,NSAME,IEDGE(1,I))
            STRX(I) = STRX(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSEIF( ITYPE(I) .EQ. IBY) THEN
            NSPLIT = NSPLIT + 1
            BSPLIT(I) = .TRUE.
            IF( BSTRN(I) ) CALL NEWEDGE(
     +           I,NENEW,NENEX,NSAME,IEDGE(1,I))
            STRX(I) = STRX(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSEIF( ITYPE(I) .EQ. IRX) THEN
            NSPLIT = NSPLIT + 1
            BSPLIT(I) = .TRUE.
            CALL RECEDGE(
     +           I,NENEW,NENEX,NSAME,IEDGE(1,I))
            STRX(I) = STRX(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSEIF( ITYPE(I) .EQ. IRY) THEN
            NSPLIT = NSPLIT + 1
            BSPLIT(I) = .TRUE.
            CALL RECEDGE(
     +           I,NENEW,NENEX,NSAME,IEDGE(1,I))
            STRX(I) = STRX(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSEIF( ITYPE(I) .EQ. IMQ) THEN
            NSPLIT = NSPLIT + 1
            BSPLIT(I) = .TRUE.
            STRX(I) = STRX(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSEIF( ITYPE(I) .EQ. ISQ) THEN
            NSPLIT = NSPLIT + 1
            BSPLIT(I) = .TRUE.
            STRX(I) = STRX(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSEIF( ITYPE(I) .EQ. IRQ) THEN 
            NSPLIT = NSPLIT + 1
            BSPLIT(I) = .TRUE.
            STRX(I) = STRX(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSEIF( ITYPE(I) .EQ. ICX) THEN
            NSPLIT = NSPLIT + 1
            BSPLIT(I) = .TRUE.
            IF( BSTRN(I) ) CALL NEWEDGE(
     +           I,NENEW,NENEX,NSAME,IEDGE(1,I))
            STRX(I) = STRX(I) *0.5D+0
            STRY(I) = STRY(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSEIF( ITYPE(I) .EQ. ICY) THEN
            NSPLIT = NSPLIT + 1
            BSPLIT(I) = .TRUE.
            IF( BSTRN(I) ) CALL NEWEDGE(
     +           I,NENEW,NENEX,NSAME,IEDGE(1,I))
            STRX(I) = STRX(I) *0.5D+0
            STRY(I) = STRY(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSEIF( ITYPE(I) .EQ. ICH) THEN
            NSPLIT = NSPLIT + 1
            BSPLIT(I) = .TRUE.
            CALL RECEDGE(
     +           I,NENEW,NENEX,NSAME,IEDGE(1,I))
            STRX(I) = STRX(I) *0.5D+0
            STRY(I) = STRY(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSEIF( ITYPE(I) .EQ. ICV) THEN
            NSPLIT = NSPLIT + 1
            CALL RECEDGE(
     +           I,NENEW,NENEX,NSAME,IEDGE(1,I))
            BSPLIT(I) = .TRUE.
            STRX(I) = STRX(I) *0.5D+0
            STRY(I) = STRY(I) *0.5D+0
            RLEN(I) = RLEN(I) *0.5D+0
            CNAME(I) = APPEND(CNAME(I),'H')
         ELSE
            BSPLIT(I) = .FALSE.            
         ENDIF

 3000 CONTINUE
 
      WRITE(16,'(1X,A,I10,A)')
     +     'GENELE : ', NSPLIT, ' element types split in typelist'

      WRITE(16,'(1X,A,I10,A,I10,A)')
     +     'GENELE : typelist increased from ', NELE,
     +     ' to ', NENEW, ' creating thin quads for bend edges'  
      NELE = NENEW

C     4: iterate original lattice list:
C       - count ITYPE in {IBX, IBY, IRX, IRY} (plain bends) -> NBEND
C       - count ITYPE in {IMQ, ISQ, IRQ} (plain quads)      -> NQUAD
C       - count ITYPE in {ICX, ICY, ICH, ICV} (comb.func.)  -> NCMBF
C       - NEWELE = NBEND + 4*NQUAD + 4*NCMBF must be .LT. MAXELE-NELE
C         => if FALSE then CALL TOFF(...)
C     4b:(note) athough RBEND types should not exist, keep testing !
C     5: copy original lattice list to temp lattice list
C        (original will be overwritten )

      NBEND = 0
      NQUAD = 0
      NCMBF = 0

      DO 4000 I=1,NLAT
         L   = LATELE(I)
         ITP = ITYPE (L)
         IF    ( ITP.EQ.IBX .OR. ITP.EQ.IBY
     +        .OR. ITP.EQ.IRX .OR. ITP.EQ.IRY ) THEN ! bend
            NBEND = NBEND + 1
         ELSEIF( ITP.EQ.IMQ .OR. ITP.EQ.ISQ .OR. ITP.EQ.IRQ) THEN ! quad
            NQUAD = NQUAD + 1
         ELSEIF( ITP.EQ.ICX .OR. ITP.EQ.ICY
     +           .OR. ITP.EQ.ICH .OR. ITP.EQ.ICV ) THEN ! combined function
            NCMBF = NCMBF + 1
         ENDIF                  ! <- ad 4:  & ad 5: ->

         LATEL0(I) = LATELE(I)
         PREDR0(I) = PREDRF(I)
         LSLI0 (I) = LSLIC (I)
 4000 CONTINUE

      NEWELE =   NBEND + 4*NQUAD +4*NCMBF
      NEWLAT = 2*NBEND + 5*NQUAD + 5*NCMBF
      WRITE(16,'(1X,A,I10,A)')
     +     'GENELE : ',NBEND,' bends found in lattice -> VD' 
      WRITE(16,'(1X,A,I10,A)')
     +     'GENELE : ',NQUAD,' quads found in lattice -> HQ,VQ,RQ,CQ' 
      WRITE(16,'(1X,A,I10,A)')
     +     'GENELE : ',NCMBF,' cmbfs found in lattice -> HQ,VQ,RQ,CQ' 
      WRITE(16,'(1X,A,I10,A)')
     +     'GENELE : ',NEWELE,' new elements needed in typelist' 
      WRITE(16,'(1X,A,I10,A)')
     +     'GENELE : total of ',NELE+NEWELE,
     +     ' elements needed in typelist' 
      WRITE(16,'(1X,A,I10,A)')
     +     'GENELE : ',NEWLAT,' new elements needed in lattice' 
      WRITE(16,'(1X,A,I10,A)')
     +     'GENELE : total of ',NLAT+NEWLAT,
     +     ' elements needed in lattice' 

      IF( NELE+NEWELE .GT. MAXELE )  CALL TOFF(
     +     'GENELE, increase MAXELE')
      IF( NLAT+NEWLAT .GT. MAXLIST ) CALL TOFF(
     +     'GENELE, increase MAXLIST')

C     6: iterate temp lattice list:
C        - if .NOT.BSPLIT => copy to lattice list
C          else 
C          -- if PLAIN BEND :
C             --- append VD_# to typelist
C   (6b)      --- if IEDGE(1,LATELE) .EQ. 0 then 
C                ---- replace BEND in lattice list by HalfBend VD HalfBEND'
C   (6b)      --- else (works for EdgeIn.eq.EdgeOut or not !!!)
C                ---- replace by EdgeIn HalfBend VD HalfBEND EdgeOut
C             --- increment # of VD
C          -- if PLAIN QUAD :
C             --- append HQ_#, VQ_#, RQ_#, CQ_#  to typelist
C             --- replace QUAD in lattice list by HalfQUAD HQ VQ RQ CQ halfQuad
C             --- increment # of HQ, VQ, RQ, CQ
C          -- if COMBINED FUNCTION :
C             --- append HQ_#, VD_#, RQ_#, CQ_#  to typelist
C   (6b)      --- if IEDGE(1,LATELE) .EQ. 0 then 
C              ---- replace CMBF in lattice list by HalfCMBF HQ VD RQ CQ HalfCMBF
C   (6b)      --- else (works for EdgeIn.eq.EdgeOut or not !!!)
C              ---- replace CMBF by EdgeIn HalfCMBF HQ VD RQ CQ HalfCMBF EdgeOut
C             --- increment # of HQ, VD, RQ, CQ

      NVD   = 0
      NHD   = 0
      NHQ   = 0
      NVQ   = 0
      NRQ   = 0
      NCQ   = 0
      NLNEW = 0

      DO 6000 I=1,NLAT
         L   = LATEL0(I)
         
         IF( BSPLIT(L) ) THEN 
            ITP = ITYPE (L)
            IF    ( ITP.EQ.IBX .OR. ITP.EQ.IBY
     +           .OR. ITP.EQ.IRX .OR. ITP.EQ.IRY ) THEN ! bend

               IF( IEDGE(1,L).EQ.0 ) THEN ! no edge effects
                  NLNEW = NLNEW + 3 ! HalfBend,VD,HalfBend 
                  IF(NLNEW.GT.MAXLIST)CALL TOFF(
     +                 'GENELE, increase MAXLIST')
                  NELE = NELE + 1
                  IF(NELE.GT.MAXELE)CALL TOFF(
     +                 'GENELE, increase MAXELE')
                  IF( ITP.EQ.IBX .OR. ITP.EQ.IRX ) THEN ! hbend : add vbend
                     NVD  = NVD + 1
                     CNAME(NELE) = 'VD'//CFINFL(NVD,'0')
                     ITYPE(NELE) = IKY
                  ELSE          ! vbend : add hbend
                     NHD  = NHD + 1
                     CNAME(NELE) = 'HD'//CFINFL(NHD,'0')
                     ITYPE(NELE) = IKX
                  ENDIF
                  STRX (NELE) = 0.D+0
                  RLEN (NELE) = 0.D+0
C                  RLEN (NELE) = DLRAD ! this is a quick and dirty fix !
                  
                  LATELE(NLNEW-2) = L 
                  PREDRF(NLNEW-2) = PREDR0(I)
                  LSLIC (NLNEW-2) = LSLI0 (I)/2 ! 1-st half bend
                  
                  LATELE(NLNEW-1) = NELE
                  PREDRF(NLNEW-1) = 0.D+0
                  LSLIC (NLNEW-1) = 0 ! central dipol kick
                  
                  LATELE(NLNEW)   = L
                  PREDRF(NLNEW)   = 0.D+0
                  LSLIC (NLNEW)   = LSLI0 (I)/2 ! 2-nd half bend
               ELSE             ! edge effects
                  NLNEW = NLNEW + 5 ! Edge HalfBend,VD,HalfBend Edge
                  IF(NLNEW.GT.MAXLIST)CALL TOFF(
     +                 'GENELE, increase MAXLIST')
                  NELE = NELE + 1
                  IF(NELE.GT.MAXELE)CALL TOFF(
     +                 'GENELE, increase MAXELE')
                  IF( ITP.EQ.IBX .OR. ITP.EQ.IRX ) THEN ! hbend : add vbend
                     NVD  = NVD + 1
                     CNAME(NELE) = 'VD'//CFINFL(NVD,'0')
                     ITYPE(NELE) = IKY
                  ELSE          ! vbend : add hbend
                     NHD  = NHD + 1
                     CNAME(NELE) = 'HD'//CFINFL(NHD,'0')
                     ITYPE(NELE) = IKX
                  ENDIF
                  STRX (NELE) = 0.D+0
                  RLEN (NELE) = 0.D+0
C                  RLEN (NELE) = DLRAD ! this is a quick and dirty fix !
                  
                  LATELE(NLNEW-4) = IEDGE(1,L)
                  PREDRF(NLNEW-4) = PREDR0(I)
                  LSLIC (NLNEW-4) = 0 ! entrance edge
                  
                  LATELE(NLNEW-3) = L 
                  PREDRF(NLNEW-3) = 0.D+0
                  LSLIC (NLNEW-3) = LSLI0 (I)/2 ! 1-st half bend
                  
                  LATELE(NLNEW-2) = NELE
                  PREDRF(NLNEW-2) = 0.D+0
                  LSLIC (NLNEW-2) = 0 ! central dipol kick
                  
                  LATELE(NLNEW-1) = L
                  PREDRF(NLNEW-1) = 0.D+0
                  LSLIC (NLNEW-1) = LSLI0 (I)/2 ! 2-nd half bend
                  
                  LATELE(NLNEW)   = IEDGE(2,L) ! .eq. iedge(1,l) or not !
                  PREDRF(NLNEW)   = 0.D+0
                  LSLIC (NLNEW)   = 0 ! exit edge
               ENDIF

            ELSEIF( ITP.EQ.IMQ .OR. ITP.EQ.ISQ .OR. ITP.EQ.IRQ) THEN ! quad

               NLNEW = NLNEW + 6 ! HalfQuad,HQ,VQ,RQ,CQ,HalfQuad 
               IF(NLNEW.GT.MAXLIST)CALL TOFF('GENELE, increase MAXLIST')
               NHQ  = NHQ + 1
               NVQ  = NVQ + 1
               NRQ  = NRQ + 1
               NCQ  = NCQ + 1
               NELE = NELE + 4
               IF(NELE.GT.MAXELE)CALL TOFF('GENELE, increase MAXELE')
               CNAME(NELE-3) = 'HQ'//CFINFL(NHQ,'0')
               ITYPE(NELE-3) = IKX
               STRX (NELE-3) = 0.D+0
               RLEN (NELE-3) = 0.D+0
C               RLEN (NELE-3) = DLRAD ! this is a quick and dirty fix !
               CNAME(NELE-2) = 'VQ'//CFINFL(NVQ,'0')
               ITYPE(NELE-2) = IKY
               STRX (NELE-2) = 0.D+0
               RLEN (NELE-2) = 0.D+0
C               RLEN (NELE-2) = DLRAD ! this is a quick and dirty fix !
               IF( ITP.EQ.IRQ ) THEN ! rotated 
                  CNAME(NELE-1) = 'RQ'//CFINFL(NRQ,'0')
                  ITYPE(NELE-1) = IRQ
                  STRX (NELE-1) = 0.D+0
                  RLEN (NELE-1) = 0.D+0
                  TWIS (NELE-1) = TWIS(L) + .25D+0 * PI
                  CNAME(NELE  ) = 'CQ'//CFINFL(NCQ,'0')
                  ITYPE(NELE  ) = IRQ
                  STRX (NELE  ) = 0.D+0
                  RLEN (NELE  ) = 0.D+0                  
                  TWIS (NELE-1) = TWIS(L)
               ELSE             ! upright or skew
                  CNAME(NELE-1) = 'RQ'//CFINFL(NRQ,'0')
                  ITYPE(NELE-1) = ISQ
                  STRX (NELE-1) = 0.D+0
                  RLEN (NELE-1) = 0.D+0
                  CNAME(NELE  ) = 'CQ'//CFINFL(NCQ,'0')
                  ITYPE(NELE  ) = IMQ
                  STRX (NELE  ) = 0.D+0
                  RLEN (NELE  ) = 0.D+0
               ENDIF

               LATELE(NLNEW-5) = L 
               PREDRF(NLNEW-5) = PREDR0(I)
               LSLIC (NLNEW-5) = LSLI0 (I)/2 ! 1-st half quad

               LATELE(NLNEW-4) = NELE-3
               PREDRF(NLNEW-4) = 0.D+0
               LSLIC (NLNEW-4) = 0           ! central horizontal dipol kick

               LATELE(NLNEW-3) = NELE-2
               PREDRF(NLNEW-3) = 0.D+0
               LSLIC (NLNEW-3) = 0           ! central vertical dipol kick

               LATELE(NLNEW-2) = NELE-1 
               PREDRF(NLNEW-2) = 0.D+0
               LSLIC (NLNEW-2) = 0           ! central skew quadrupole kick

               LATELE(NLNEW-1) = NELE
               PREDRF(NLNEW-1) = 0.D+0
               LSLIC (NLNEW-1) = 0           ! central quadrupole kick

               LATELE(NLNEW)   = L
               PREDRF(NLNEW)   = 0.D+0
               LSLIC (NLNEW)   = LSLI0 (I)/2 ! 2-nd half quad

            ELSEIF( ITP.EQ.ICX .OR. ITP.EQ.ICY
     +              .OR. ITP.EQ.ICH .OR. ITP.EQ.ICV ) THEN ! combined function

               IF( IEDGE(1,L).EQ.0 ) THEN ! no edge effects
                  NLNEW = NLNEW + 6 ! HalfCMBF,HQ,VQ,RQ,CQ,HalfCMBF 
                  IF(NLNEW.GT.MAXLIST)CALL TOFF(
     +                 'GENELE, increase MAXLIST')
                  NRQ  = NRQ + 1
                  NCQ  = NCQ + 1
                  NELE = NELE + 4
                  IF(NELE.GT.MAXELE)CALL TOFF(
     +                 'GENELE, increase MAXELE')
                  IF( ITP.EQ.ICX .OR. ITP.EQ.ICH ) THEN ! hcmbf : add VD and HQ
                     NHQ  = NHQ + 1
                     NVD  = NVD + 1
                     CNAME(NELE-3) = 'HQ'//CFINFL(NHQ,'0')
                     ITYPE(NELE-3) = IKX
                     STRX (NELE-3) = 0.D+0
                     RLEN (NELE-3) = 0.D+0
C                     RLEN (NELE-3) = DLRAD ! this is a quick and dirty fix !
                     CNAME(NELE-2) = 'VD'//CFINFL(NVD,'0')
                     ITYPE(NELE-2) = IKY
                     STRX (NELE-2) = 0.D+0
                     RLEN (NELE-2) = 0.D+0
C                     RLEN (NELE-2) = DLRAD ! this is a quick and dirty fix !
                  ELSE          ! vcmbf : add HD and VQ
                     NHD  = NHD + 1
                     NVQ  = NVQ + 1
                     CNAME(NELE-3) = 'HD'//CFINFL(NHD,'0')
                     ITYPE(NELE-3) = IKX
                     STRX (NELE-3) = 0.D+0
                     RLEN (NELE-3) = 0.D+0
C                     RLEN (NELE-3) = DLRAD ! this is a quick and dirty fix !
                     CNAME(NELE-2) = 'VQ'//CFINFL(NVQ,'0')
                     ITYPE(NELE-2) = IKY
                     STRX (NELE-2) = 0.D+0
                     RLEN (NELE-2) = 0.D+0
C                     RLEN (NELE-2) = DLRAD ! this is a quick and dirty fix !
                  ENDIF   
                  CNAME(NELE-1) = 'RQ'//CFINFL(NRQ,'0')
                  ITYPE(NELE-1) = ISQ
                  STRX (NELE-1) = 0.D+0
                  RLEN (NELE-1) = 0.D+0
                  CNAME(NELE  ) = 'CQ'//CFINFL(NCQ,'0')
                  ITYPE(NELE  ) = IMQ
                  STRX (NELE  ) = 0.D+0
                  RLEN (NELE  ) = 0.D+0
                  
                  LATELE(NLNEW-5) = L 
                  PREDRF(NLNEW-5) = PREDR0(I)
                  LSLIC (NLNEW-5) = LSLI0 (I)/2 ! 1-st half comb.func
                  
                  LATELE(NLNEW-4) = NELE-3
                  PREDRF(NLNEW-4) = 0.D+0
                  LSLIC (NLNEW-4) = 0 ! central horizontal dipol kick
                  
                  LATELE(NLNEW-3) = NELE-2
                  PREDRF(NLNEW-3) = 0.D+0
                  LSLIC (NLNEW-3) = 0 ! central vertical dipol kick
                  
                  LATELE(NLNEW-2) = NELE-1 
                  PREDRF(NLNEW-2) = 0.D+0
                  LSLIC (NLNEW-2) = 0 ! central skew quadrupole kick
                  
                  LATELE(NLNEW-1) = NELE
                  PREDRF(NLNEW-1) = 0.D+0
                  LSLIC (NLNEW-1) = 0 ! central quadrupole kick
                  
                  LATELE(NLNEW)   = L
                  PREDRF(NLNEW)   = 0.D+0
                  LSLIC (NLNEW)   = LSLI0 (I)/2 ! 2-nd half comb.func
               ELSE             ! edge effects
                  NLNEW = NLNEW + 8 ! Edge HalfCMBF,HQ,VQ,RQ,CQ,HalfCMBF Edge
                  IF(NLNEW.GT.MAXLIST)CALL TOFF(
     +                 'GENELE, increase MAXLIST')
                  NRQ  = NRQ + 1
                  NCQ  = NCQ + 1
                  NELE = NELE + 4
                  IF(NELE.GT.MAXELE)CALL TOFF(
     +                 'GENELE, increase MAXELE')
                  IF( ITP.EQ.ICX .OR. ITP.EQ.ICH ) THEN ! hcmbf : add VD and HQ
                     NHQ  = NHQ + 1
                     NVD  = NVD + 1
                     CNAME(NELE-3) = 'HQ'//CFINFL(NHQ,'0')
                     ITYPE(NELE-3) = IKX
                     STRX (NELE-3) = 0.D+0
                     RLEN (NELE-3) = 0.D+0
C                     RLEN (NELE-3) = DLRAD ! this is a quick and dirty fix !
                     CNAME(NELE-2) = 'VD'//CFINFL(NVD,'0')
                     ITYPE(NELE-2) = IKY
                     STRX (NELE-2) = 0.D+0
                     RLEN (NELE-2) = 0.D+0
C                     RLEN (NELE-2) = DLRAD ! this is a quick and dirty fix !
                  ELSE          ! vcmbf : add HD and VQ
                     NHD  = NHD + 1
                     NVQ  = NVQ + 1
                     CNAME(NELE-3) = 'HD'//CFINFL(NHD,'0')
                     ITYPE(NELE-3) = IKX
                     STRX (NELE-3) = 0.D+0
                     RLEN (NELE-3) = 0.D+0
C                     RLEN (NELE-3) = DLRAD ! this is a quick and dirty fix !
                     CNAME(NELE-2) = 'VQ'//CFINFL(NVQ,'0')
                     ITYPE(NELE-2) = IKY
                     STRX (NELE-2) = 0.D+0
                     RLEN (NELE-2) = 0.D+0
C                     RLEN (NELE-2) = DLRAD ! this is a quick and dirty fix !
                  ENDIF   
                  CNAME(NELE-1) = 'RQ'//CFINFL(NRQ,'0')
                  ITYPE(NELE-1) = ISQ
                  STRX (NELE-1) = 0.D+0
                  RLEN (NELE-1) = 0.D+0
                  CNAME(NELE  ) = 'CQ'//CFINFL(NCQ,'0')
                  ITYPE(NELE  ) = IMQ
                  STRX (NELE  ) = 0.D+0
                  RLEN (NELE  ) = 0.D+0

                  LATELE(NLNEW-7) = IEDGE(1,L)
                  PREDRF(NLNEW-7) = PREDR0(I)
                  LSLIC (NLNEW-7) = 0 ! entrance edge

                  LATELE(NLNEW-6) = L 
                  PREDRF(NLNEW-6) = 0.D+0
                  LSLIC (NLNEW-6) = LSLI0 (I)/2 ! 1-st half comb.func
                  
                  LATELE(NLNEW-5) = NELE-3
                  PREDRF(NLNEW-5) = 0.D+0
                  LSLIC (NLNEW-5) = 0 ! central horizontal dipol kick
                  
                  LATELE(NLNEW-4) = NELE-2
                  PREDRF(NLNEW-4) = 0.D+0
                  LSLIC (NLNEW-4) = 0 ! central vertical dipol kick
                  
                  LATELE(NLNEW-3) = NELE-1 
                  PREDRF(NLNEW-3) = 0.D+0
                  LSLIC (NLNEW-3) = 0 ! central skew quadrupole kick
                  
                  LATELE(NLNEW-2) = NELE
                  PREDRF(NLNEW-2) = 0.D+0
                  LSLIC (NLNEW-2) = 0 ! central quadrupole kick
                  
                  LATELE(NLNEW-1) = L
                  PREDRF(NLNEW-1) = 0.D+0
                  LSLIC (NLNEW-1) = LSLI0 (I)/2 ! 2-nd half comb.func

                  LATELE(NLNEW)   = IEDGE(2,L)
                  PREDRF(NLNEW)   = 0.D+0
                  LSLIC (NLNEW)   = 0 ! exit edge

               ENDIF
            ENDIF            
         ELSE
            NLNEW = NLNEW + 1
            IF( NLNEW.GT.MAXLIST ) CALL TOFF('GENELE, increase MAXLIST')
            LATELE(NLNEW) = L
            PREDRF(NLNEW) = PREDR0(I)
            LSLIC (NLNEW) = LSLI0 (I)
         ENDIF
 6000 CONTINUE

      NLAT = NLNEW

C     7: add END in typelist and lattice list

      NELE  = NELE + 1
      IF( NELE.GT.MAXELE ) CALL TOFF('GENELE, increase MAXELE for END')
      CNAME(NELE) = 'X______EOL'
      ITYPE(NELE) = IMR
      NLAT  = NLAT + 1
      IF(NLAT.GT.MAXLIST) CALL TOFF('GENELE, increase MAXLIST for END')
      LATELE(NLAT) = NELE
      PREDRF(NLAT) = PRESAVE


      RETURN

      E N D 


      CHARACTER*(*) F U N C T I O N APPOVWR(CS1,CS2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  APPOVWR
C@ AUTHORS    :  M.Vogt
C@ ARGS/IN    :  CS1,CS2 :    
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  ICSEND
C@ IO_UNITS   :  NONE
C@ DESCRIPTION:  If cs1 holds enough unused space it appends cs2 after the last
C-               non-blank char of cs1. If cs1 is too short it overwrites as much
C-               of cs1's back as is needed to contain cs2. Warns on overwrite!
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CS1*(*), CS2*(*)

C%%      write(*,*)'appovwr'
      
      LN1 = LEN   (CS1)         ! length of cs1
      IE1 = ICSEND(CS1)         ! used length of cs1
      LF1 = LN1 - IE1           ! free space in cs1
      IE2 = ICSEND(CS2)         ! used length of cs2

      APPOVWR = CS1

      IF( IE2.GT.LN1 ) CALL TOFF(
     +     'APPOVWR, attempt to overwrite beyond start')

      IF( LF1.GE.IE2 ) THEN
         APPOVWR( IE1+1    :IE1+IE2 ) = CS2( 1:IE2 )
      ELSE
         CALL WARN('APPOVWR, string truncated :')
         APPOVWR( LN1-IE2+1:LN1     ) = CS2( 1:IE2 )
         WRITE(16,'(2X,A,A)') CS1(1:IE1) ,' ->', APPOVWR
      ENDIF

      RETURN
      E N D

      CHARACTER*(*) F U N C T I O N APPEND(CS1,CS2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  APPEND
C@ AUTHORS    :  M.Vogt
C@ ARGS/IN    :  CS1,CS2 :    
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  ICSEND
C@ IO_UNITS   :  NONE
C@ DESCRIPTION:  Appends cs2 after the last non-blank char of cs1. 
C-               If cs1 is too short, only the first part of cs2 is used.
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CS1*(*), CS2*(*)

C%%      write(*,*)'append'
      
      LN1 = LEN   (CS1)         ! length of cs1
      IE1 = ICSEND(CS1)         ! used length of cs1
      APPEND = CS1
      APPEND(IE1+1:LN1) = CS2

      RETURN
      E N D


      S U B R O U T I N E NEWEDGE(IELE,NENEW,NENEX,NSAME,IEDGE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  NEWEDGE
C@ AUTHORS    :  M.Vogt
C@ ARGS/IN    :  IELE    
C@ ARGS/INOUT :  NENEW,NENEX,NSAME,IEDGE
C@ ARGS/OUT   :  
C@ EXTERNAL   :  TOFF
C@ IO_UNITS   :  16
C@ DESCRIPTION:  Creates edge focusing quads for SBEND like type with rotated
C-               polefaces
C@ REMARKS    :  CHECK FOR TREATMENT OF VERTICAL TYPES !!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====local stuff
      PARAMETER(MAXNUM=MAELNM-4)
      CHARACTER  CFINFL*(MAXNUM)
      DIMENSION IEDGE(2)

C=====statementfunctions ...
C=====... for dealing with edge effects in SBENDs
      BSTRN(I) = STRN(1,I).NE.0.D+0 .OR. STRN(2,I).NE.0.D+0 
     +     .OR. STRN(3,I).NE.0.D+0 .OR. STRN(4,I).NE.0.D+0  ! edge treatment
      BFRIN(I) = STRN(3,I).NE.0.D+0 .OR. STRN(4,I).NE.0.D+0 ! fringe corr: x.ne.y
      BSAME(I) = STRN(1,I) .EQ. STRN(2,I)                   ! EdgeIn=EdgeOut
      DKEX(E,ANGLE,DL) = ANGLE*DTAN(E)/DL
      DKEY(E,ANGLE,DL,HGAP,FINT) =
     +     ANGLE*DTAN( E - ANGLE*HGAP*FINT/DL * (1.D+0+DSIN(E)**2) )/DL

C%%      write(*,*)'newedge'

      IT = ITYPE(IELE)

      IF( BSAME(IELE) ) THEN
         NENEW = NENEW +1 
         IF( NENEW .GT. MAXELE )  CALL TOFF(
     +        'GENELE, increase MAXELE')
         NSAME = NSAME +1
         CNAME(NENEW) = 'EE'//CFINFL(NSAME,'0')
         ITYPE(NENEW) = ITQ

         IF( IT.EQ.IBX .OR. IT.EQ.ICX ) THEN ! horizontal
            STRX(NENEW) = -DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))
            IF( BFRIN(IELE) ) THEN
               STRY(NENEW) =  DKEY(STRN(1,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
            ELSE
               STRY(NENEW) = 0.D+0
            ENDIF
         ELSE                                ! vertical
            IF( BFRIN(IELE) ) THEN
               STRX(NENEW) = DKEY(STRN(1,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
               STRY(NENEW) = -DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))
            ELSE
               STRX(NENEW) = DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))
               STRY(NENEW) = 0.D+0
            ENDIF
         ENDIF

         IEDGE(1) = NENEW
         IEDGE(2) = NENEW
      ELSE
         NENEW = NENEW +2 
         IF( NENEW .GT. MAXELE )  CALL TOFF(
     +        'GENELE, increase MAXELE')
         NENEX = NENEX +1
         CNAME(NENEW-1) = 'EI'//CFINFL(NENEX,'0')
         ITYPE(NENEW-1) = ITQ
         CNAME(NENEW)   = 'EO'//CFINFL(NENEX,'0')
         ITYPE(NENEW)   = ITQ

         IF( IT.EQ.IBX .OR. IT.EQ.ICX ) THEN ! horizontal
            STRX(NENEW-1) =-DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))            
            STRX(NENEW)   =-DKEX(STRN(2,IELE),STRX(IELE),RLEN(IELE))
            IF( BFRIN(IELE) ) THEN
               STRY(NENEW-1) = DKEY(STRN(1,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
               STRY(NENEW)   = DKEY(STRN(2,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
            ELSE
               STRY(NENEW-1) = 0.D+0
               STRY(NENEW)   = 0.D+0
            ENDIF
         ELSE                                ! vertical
            IF( BFRIN(IELE) ) THEN
               STRX(NENEW-1) = DKEY(STRN(1,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
               STRY(NENEW-1) =-DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))
               STRX(NENEW)   = DKEY(STRN(2,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
               STRY(NENEW)   =-DKEX(STRN(2,IELE),STRX(IELE),RLEN(IELE))
            ELSE
               STRX(NENEW-1) = DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))
               STRY(NENEW-1) = 0.D+0
               STRX(NENEW)   = DKEX(STRN(2,IELE),STRX(IELE),RLEN(IELE))
               STRY(NENEW)   = 0.D+0
            ENDIF
         ENDIF

         IEDGE(1) = NENEW-1
         IEDGE(2) = NENEW
      ENDIF

      WRITE(16,'(1X,3A,I10,A/2X,3A,I10/2X,3A,I10)')
     +     'SBEND : ',CNAME(IELE),' at TPL-index : ',IELE,' Edges >>',
     +     'EdgeIn  : ',CNAME(IEDGE(1)),' at TPL-index : ',IEDGE(1),
     +     'EdgeOut : ',CNAME(IEDGE(2)),' at TPL-index : ',IEDGE(2)


      CALL ZERO(MXSTRN,STRN(1,IELE)) ! do not apply edge effects twice! 

      RETURN
      E N D


      S U B R O U T I N E RECEDGE(IELE,NENEW,NENEX,NSAME,IEDGE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  RECEDGE
C@ AUTHORS    :  M.Vogt
C@ ARGS/IN    :  IELE    
C@ ARGS/INOUT :  NENEW,NENEX,NSAME,IEDGE
C@ ARGS/OUT   :  
C@ EXTERNAL   :  TOFF, ZERO
C@ IO_UNITS   :  16
C@ DESCRIPTION:  Creates edge focusing quads for RBEND like type with rotated
C-               polefaces
C@ REMARKS    :  CHECK FOR TREATMENT OF VERTICAL TYPES !!!!!!!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====local stuff
      PARAMETER(MAXNUM=MAELNM-4)
      CHARACTER  CFINFL*(MAXNUM)
      DIMENSION IEDGE(2)

C=====statementfunctions ...
C=====... for dealing with edge effects in RBENDs
      BSTRN(I) = STRN(1,I).NE.0.D+0 .OR. STRN(2,I).NE.0.D+0 
     +     .OR. STRN(3,I).NE.0.D+0 .OR. STRN(4,I).NE.0.D+0  ! edge treatment
      BFRIN(I) = STRN(3,I).NE.0.D+0 .OR. STRN(4,I).NE.0.D+0 ! fringe corr: x.ne.y
      BSAME(I) = STRN(1,I) .EQ. STRN(2,I)                   ! EdgeIn=EdgeOut
      DKEX(E,ANGLE,DL) = ANGLE*DTAN(E+ANGLE/2.D+0)/DL
      DKEY(E,ANGLE,DL,HGAP,FINT) =
     +     ANGLE*DTAN( E+ANGLE/2.D+0 - ANGLE*HGAP*FINT/DL *
     +     (1.D+0+DSIN(E+ANGLE/2.D+0)**2) )/DL

C%%      write(*,*)'recedge'

      IF( BSAME(IELE) ) THEN
         NENEW = NENEW +1 
         IF( NENEW .GT. MAXELE )  CALL TOFF(
     +        'GENELE, increase MAXELE')
         NSAME = NSAME + 1
         CNAME(NENEW) = 'EE'//CFINFL(NSAME,'0')
         ITYPE(NENEW) = ITQ

         IF( IT.EQ.IRX .OR. IT.EQ.ICH ) THEN ! horizontal
            STRX(NENEW) = -DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))
            IF( BFRIN(IELE) ) THEN
               STRY(NENEW) =  DKEY(STRN(1,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
            ELSE
               STRY(NENEW) = 0.D+0
            ENDIF
         ELSE                                ! vertical
            IF( BFRIN(IELE) ) THEN
               STRX(NENEW) = DKEY(STRN(1,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
               STRY(NENEW) = -DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))
            ELSE
               STRX(NENEW) = DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))
               STRY(NENEW) = 0.D+0
            ENDIF
         ENDIF

         IEDGE(1) = NENEW
         IEDGE(2) = NENEW
      ELSE
         NENEW = NENEW +2 
         IF( NENEW .GT. MAXELE )  CALL TOFF(
     +        'GENELE, increase MAXELE')
         NENEX = NENEX +1
         CNAME(NENEW-1) = 'EI'//CFINFL(NENEX,'0')
         ITYPE(NENEW-1) = ITQ
         CNAME(NENEW)   = 'EO'//CFINFL(NENEX,'0')
         ITYPE(NENEW)   = ITQ

         IF( IT.EQ.IBX .OR. IT.EQ.ICX ) THEN ! horizontal
            STRX(NENEW-1) =-DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))            
            STRX(NENEW)   =-DKEX(STRN(2,IELE),STRX(IELE),RLEN(IELE))
            IF( BFRIN(IELE) ) THEN
               STRY(NENEW-1) = DKEY(STRN(1,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
               STRY(NENEW)   = DKEY(STRN(2,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
            ELSE
               STRY(NENEW-1) = 0.D+0
               STRY(NENEW)   = 0.D+0
            ENDIF
         ELSE                                ! vertical
            IF( BFRIN(IELE) ) THEN
               STRX(NENEW-1) = DKEY(STRN(1,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
               STRY(NENEW-1) =-DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))
               STRX(NENEW)   = DKEY(STRN(2,IELE),STRX(IELE),RLEN(IELE),
     +              STRN(3,IELE),STRN(4,IELE))
               STRY(NENEW)   =-DKEX(STRN(2,IELE),STRX(IELE),RLEN(IELE))
            ELSE
               STRX(NENEW-1) = DKEX(STRN(1,IELE),STRX(IELE),RLEN(IELE))
               STRY(NENEW-1) = 0.D+0
               STRX(NENEW)   = DKEX(STRN(2,IELE),STRX(IELE),RLEN(IELE))
               STRY(NENEW)   = 0.D+0
            ENDIF
         ENDIF

         IEDGE(1) = NENEW-1
         IEDGE(2) = NENEW
      ENDIF

      WRITE(16,'(1X,3A,I10,A/2X,3A,I10/2X,3A,I10)')
     +     'RBEND : ',CNAME(IELE),' at TPL-index : ',IELE,' Edges >>',
     +     'EdgeIn  : ',CNAME(IEDGE(1)),' at TPL-index : ',IEDGE(1),
     +     'EdgeOut : ',CNAME(IEDGE(2)),' at TPL-index : ',IEDGE(2)

C=====now it's become an SBEND type
      IF    ( ITYPE(IELE).EQ.IRX ) THEN
         ITYPE(IELE) = IBX
      ELSEIF( ITYPE(IELE).EQ.IRY ) THEN
         ITYPE(IELE) = IBY
      ELSEIF( ITYPE(IELE).EQ.ICH ) THEN
         ITYPE(IELE) = ICX
      ELSEIF( ITYPE(IELE).EQ.ICV ) THEN
         ITYPE(IELE) = ICY
      ELSE
         CALL TOFF('RECEDGE, this should be impossible !')
      ENDIF

      CALL ZERO(MXSTRN,STRN(1,IELE)) ! do not apply edge effects twice! 

      RETURN
      E N D



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 4                                                                           C
C@ PROGRAM=SPRINT, MODULE=COMMAND, SUBODULE=UTIL_LAT, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      S U B R O U T I N E WRMAD0(CSEQ0,CREF0,CSXPND,CSNAM0)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  WRMAD0
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ COMMAND    :  
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     writes the  c u r r e n t  status of /ELEIN/ and /LATTICE/ into
C     CMAD0F.mad0.VNR
C     !!! write_mad0_lattice '$mylattice' '$MYSEQUENCE'
C         ${'ENTRY','CENTRE','EXIT'} $expand_superperiod
C         ${'distance','number'} ;    !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
C=====In order to compile SPRINT on 'AIn`t uniX' MBUF is just 400 here !
C=====MAD supplies a comand buffer of 4000 !!!
C=====DON'T YOU FORGET TO FORGET ABOUT IBM !!!
      PARAMETER(MM0COL=80,MAXVNR=99,MALPH=26,DM2D=10.D+0,
     +     MXNRPO=6,MBUF=400)   !>>>> MXNRPRO = 6 (normal), 8 (large)
      DIMENSION BALPH(MALPH,MALPH)
      DIMENSION NUMBER(MAXELE)
      CHARACTER CSEQ0*(*), CREF0*(*)
      CHARACTER CSEQNM*15, CREFER*6
      CHARACTER CREGMR*2, CFINFL*(MXNRPO), CSNAME*8
      CHARACTER CALPH*(MALPH), CREG*2, CRPOS*(MXNRPO),
     +     CLATEL*(MAELNM+MXNRPO), CDUM*10, CBUF*(MBUF), 
     +     CCLASS*(MAELNM), CSXPND*(*), CSNAM0*(*), CASEUP*5
      DATA CALPH /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

C%%      write(*,*)'wrmad0'

C=====preliminaries.......
      CALL INFO('WRMAD0')

      IF(ICSEND(CREF0)-ICSBEG(CREF0) .GT. MAXNME)
     +     CALL TOFF('WRMAD0, referencing mode invalid')
      CREFER = CREF0(ICSBEG(CREF0):ICSEND(CREF0))

      IF(CREFER.EQ.'ENTRY') THEN
         DREFER = 0.D+0
      ELSEIF(CREFER.EQ.'CENTRE') THEN
         DREFER = 0.5D+0
      ELSEIF(CREFER.EQ.'CENTER') THEN
C=====(The Stormtroopers Of Death (S.O.D '87): speak english or die !)
         CREFER = 'CENTRE'
         DREFER = 0.5D+0
      ELSEIF(CREFER.EQ.'EXIT') THEN
         DREFER = 1.D+0
      ELSE
         CALL TOFF('WRMAD0, invalid referencing mode: '//CREFER)
      ENDIF

      IF(ICSEND(CSNAM0)-ICSBEG(CSNAM0) .GT. MAXNME)
     +     CALL TOFF('WRMAD0, sequence naming  mode invalid')
      CSNAME = CSNAM0(ICSBEG(CSNAM0):ICSEND(CSNAM0))

      IF    (CSNAME.EQ.'distance') THEN
         BDIST = .TRUE.
         BNUMB = .FALSE.
      ELSEIF(CSNAME.EQ.'number'  ) THEN
         BDIST = .FALSE.
         BNUMB = .TRUE.
      ELSE
         CALL TOFF('WRMAD0, invalid referencing mode: '//CREFER)
      ENDIF

      IF(ICSEND(CSEQ0)-ICSBEG(CSEQ0) .GT. 15)
     +     CALL TOFF('WRMAD0, sequence-name too long')
      CSEQNM = CSEQ0(ICSBEG(CSEQ0):ICSEND(CSEQ0))
      IF(ICSEND(CSEQNM).LE.MAELNM .AND. NUMELE(CSEQNM).NE.0)
     +     CALL TOFF('WRMAD0, sequence name not unique: '//CSEQNM)

      IUMAD0 = 17


      IBXP = ICSBEG(CSXPND)
      IEXP = ICSEND(CSXPND)
      IF(IEXP-IBXP+1.GT.5) CALL TOFF('WRMAD0, logical must be'//
     +     ' ''true'' or ''false''')

      IF    (CASEUP(CSXPND(IBXP:IEXP)).EQ.'TRUE' .AND. NSUPER.GT.1) THEN
         WRITE(CDUM,'(I10)') NSUPER
         CALL WARN('WRMAD0, superperiod'//CDUM//' will be expanded')
         NSXP = NSUPER
      ELSEIF(CASEUP(CSXPND(IBXP:IEXP)).EQ.'FALSE' .OR. NSUPER.EQ.1) THEN
         NSXP = 1
      ELSE
         CALL TOFF('WRMAD0, logical must be ''true'' or ''false''')
      ENDIF

      
C=====mad0-header
      WRITE(IUMAD0,'(A,A)') '! ', CLATOV(:MIN(ICSEND(CLATOV),MM0COL))
      WRITE(IUMAD0,'(A,A)') '! job: ',
     +     CJOB(:MIN(ICSEND(CJOB),MM0COL-4))
      WRITE(IUMAD0,'(A,A)') '! original: ',
     +     CLATI(:MIN(ICSEND(CLATI),MM0COL-12))

C=====type-list
      
      DO 10 I=1,NELE
         NUMBER(I) = 1          ! for numbering mode in sequence
         CBUF = ' '
         ITY  = ITYPE(I)
         ITY2 = ITPEFF(I)
         IFH  = IFLAG(I)
         IPH  = IPS(I)
         ISH  = ISEC(I) * IGLBSL

         IF(ISH.GT.1) THEN
            DSH  = DBLE(ISH)
            RL   = RLEN(I) * DSH
            SX   = STRX(I) * DSH
            SY   = STRY(I) * DSH
         ELSE
            ISH  = 1
            RL   = RLEN(I)
            SX   = STRX(I)
            SY   = STRY(I)
         ENDIF

         IF(  ITY.EQ.IBX .OR. ITY.EQ.IBY .OR. ITY.EQ.IRX .OR.ITY.EQ.IRY
     +        .OR.
     +        ITY.EQ.ICX .OR. ItY.EQ.ICY .OR. ITY.EQ.ICH .OR.ITY.EQ.ICV
     +        ) THEN
            IF(  STRN(1,I).NE.0.D+0 .OR. STRN(2,I).NE.0.D+0 .OR.
     +           STRN(3,I).NE.0.D+0 .OR. STRN(4,I).NE.0.D+0 ) THEN
               BFRINGE = .TRUE.
            ELSE
               BFRINGE = .FALSE.
            ENDIF
         ENDIF


         IF    (ITY.EQ.IMR) THEN
            WRITE(CBUF,'(1X,2A)')
     +           CNAME(I),': MARKER'
         ELSEIF(ITY.EQ.IDL) THEN
            WRITE(CBUF,'(1X,A,1(A,E13.7))')
     +           CNAME(I),': DRIFT, L=',RL
            CALL EXTMAD(CBUF,ITY,0,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IMO) THEN
            WRITE(CBUF,'(1X,A,1(A,E13.7))')
     +           CNAME(I),': MONITOR, L=',RL
            CALL EXTMAD(CBUF,ITY,0,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IBX) THEN
            WRITE(CBUF,'(1X,A,A,E13.7,A,E22.15)')
     +           CNAME(I),': SBEND, L=',RL,', ANGLE=',SX
            IF( BFRINGE ) THEN
               IEBFF = ICSEND(CBUF)
               WRITE( CBUF(IEBFF+1:), '(4(A,E13.7))' ) 
     +              ', E1=', STRN(1,I), ', E2=', STRN(2,I),
     +              ', HGAP=', STRN(3,I), ', FINT=', STRN(4,I)
            ENDIF
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IRX) THEN
            WRITE(CBUF,'(1X,A,A,E13.7,A,E22.15)')
     +           CNAME(I),': RBEND, L=',RL,', ANGLE=',SX
            IF( BFRINGE ) THEN
               IEBFF = ICSEND(CBUF)
               WRITE( CBUF(IEBFF+1:), '(4(A,E13.7))' ) 
     +              ', E1=', STRN(1,I), ', E2=', STRN(2,I),
     +              ', HGAP=', STRN(3,I), ', FINT=', STRN(4,I)
            ENDIF
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IMQ) THEN
            WRITE(CBUF,'(1X,A,A,E13.7,A,E15.9)')
     +           CNAME(I),': QUADRUPOLE, L=',RL,', K1=',DIVI0(SX,RL)
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.ISQ) THEN
            WRITE(CBUF,'(1X,A,A,E13.7,A,E15.9)') 
     +           CNAME(I),': QUADRUPOLE, TILT, L=',RL,', K1=',
     +           DIVI0(SX,RL)
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IRF) THEN
            WRITE(CBUF,'(1X,A,3(A,E13.7))')
     +           CNAME(I),': RFCAVITY, L=',RL,', VOLT=',SX*1.D+3,
     +           ', FREQ=',SY*1.D-6
            CALL EXTMAD(CBUF,ITY,IPH,0,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IKX) THEN
            WRITE(CBUF,'(1X,A,2(A,E13.7))')
     +           CNAME(I),': HKICKER, L=',RL,', KICK=',SX
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IKY) THEN
            WRITE(CBUF,'(1X,A,2(A,E13.7))')
     +           CNAME(I),': VKICKER, L=',RL,', KICK=',SX
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IMS) THEN
            WRITE(CBUF,'(1X,A,A,E13.7,A,E15.9)')
     +           CNAME(I),': SEXTUPOLE, L=',RL,', K2=',DIVI0(SX,RL)
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IBY) THEN
            WRITE(CBUF,'(1X,A,A,E13.7,A,E22.15)')
     +           CNAME(I),': SBEND, TILT, L=',RL,', ANGLE=',-SX
            IF( BFRINGE ) THEN
               IEBFF = ICSEND(CBUF)
               WRITE( CBUF(IEBFF+1:), '(4(A,E13.7))' ) 
     +              ', E1=', STRN(1,I), ', E2=', STRN(2,I),
     +              ', HGAP=', STRN(3,I), ', FINT=', STRN(4,I)
            ENDIF
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IRY) THEN
            WRITE(CBUF,'(1X,A,A,E13.7,A,E22.15)')
     +           CNAME(I),': RBEND, TILT, L=',RL,', ANGLE=',-SX
            IF( BFRINGE ) THEN
               IEBFF = ICSEND(CBUF)
               WRITE( CBUF(IEBFF+1:), '(4(A,E13.7))' ) 
     +              ', E1=', STRN(1,I), ', E2=', STRN(2,I),
     +              ', HGAP=', STRN(3,I), ', FINT=', STRN(4,I)
            ENDIF
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.ISO) THEN
            WRITE(CBUF,'(1X,A,A,E13.7,A,E15.9)')
     +           CNAME(I),': SOLENOID, L=',RL,', KS=',SX
            CALL EXTMAD(CBUF,ITY,IPH,0,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IRQ) THEN
            WRITE(CBUF,'(1X,A,2(A,E13.7),A,E15.9)')
     +           CNAME(I),': QUADRUPOLE, TILT=',TWIS(I),', L=',RL,
     +           ', K1=',DIVI0(SX,RL)
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IRO) THEN
            CALL WARN('WRMAD0, can''t convert rotator '//
     +           'with itype=14')
         ELSEIF(ITY.EQ.ICX) THEN

            IF(ITY2.EQ.IMQ) THEN
               SY = SX
               SX = 0.D+0
            ENDIF

            WRITE(CBUF,'(1X,A,A,E13.7,2(A,E22.15))')
     +           CNAME(I),': SBEND, L=',RL,', ANGLE=',SX,
     +           ', K1=',DIVI0(SY,RL)
            IF( BFRINGE ) THEN
               IEBFF = ICSEND(CBUF)
               WRITE( CBUF(IEBFF+1:), '(4(A,E13.7))' ) 
     +              ', E1=', STRN(1,I), ', E2=', STRN(2,I),
     +              ', HGAP=', STRN(3,I), ', FINT=', STRN(4,I)
            ENDIF
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.ICY) THEN

            IF(ITY2.EQ.IMQ) THEN
               SY = SX
               SX = 0.D+0
            ENDIF

            WRITE(CBUF,'(1X,A,A,E13.7,2(A,E22.15))')
     +           CNAME(I),': SBEND, TILT, L=',RL,', ANGLE=',-SX,
     +           ', K1=',DIVI0(SY,RL)
            IF( BFRINGE ) THEN
               IEBFF = ICSEND(CBUF)
               WRITE( CBUF(IEBFF+1:), '(4(A,E13.7))' ) 
     +              ', E1=', STRN(1,I), ', E2=', STRN(2,I),
     +              ', HGAP=', STRN(3,I), ', FINT=', STRN(4,I)
            ENDIF
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.ICH) THEN

            IF(ITY2.EQ.IMQ) THEN
               SY = SX
               SX = 0.D+0
            ENDIF

            WRITE(CBUF,'(1X,A,A,E13.7,2(A,E22.15))')
     +           CNAME(I),': RBEND, L=',RL,', ANGLE=',SX,
     +           ', K1=',DIVI0(SY,RL)
            IF( BFRINGE ) THEN
               IEBFF = ICSEND(CBUF)
               WRITE( CBUF(IEBFF+1:), '(4(A,E13.7))' ) 
     +              ', E1=', STRN(1,I), ', E2=', STRN(2,I),
     +              ', HGAP=', STRN(3,I), ', FINT=', STRN(4,I)
            ENDIF
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.ICV) THEN

            IF(ITY2.EQ.IMQ) THEN
               SY = SX
               SX = 0.D+0
            ENDIF

            WRITE(CBUF,'(1X,A,A,E13.7,2(A,E22.15))')
     +           CNAME(I),': RBEND, TILT, L=',RL,', ANGLE=',-SX,
     +           ', K1=',DIVI0(SY,RL)
            IF( BFRINGE ) THEN
               IEBFF = ICSEND(CBUF)
               WRITE( CBUF(IEBFF+1:), '(4(A,E13.7))' ) 
     +              ', E1=', STRN(1,I), ', E2=', STRN(2,I),
     +              ', HGAP=', STRN(3,I), ', FINT=', STRN(4,I)
            ENDIF
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.ISS) THEN
            WRITE(CBUF,'(''!@ '',A,A,E13.7,3(A,F9.4))')
     +           CNAME(I),': SPINROTATOR, L=',RL,', PSI=',TWIS(I),
     +           ', PHI=',SX,', THETA=',SY
            CALL EXTMAD(CBUF,ITY,IPH,0,IFH,0.D+0,.FALSE.)
         ELSEIF(ITY.EQ.ISF) THEN
            WRITE(CBUF,'(''!@ '',A,A,E13.7,3(A,F9.4))')
     +           CNAME(I),': SPINROTATOR, FREE, L=',RL,
     +           ', PSI=',TWIS(I),', PHI=',SX,', THETA=',SY
            CALL EXTMAD(CBUF,ITY,0,0,IFH,0.D+0,.FALSE.)
         ELSEIF(ITY.EQ.ISM) THEN
            WRITE(CBUF,'(''!@ '',A,A,E13.7)')
     +           CNAME(I),': SPINMATCH, L=',RL
            CALL EXTMAD(CBUF,ITY,0,0,IFH,0.D+0,.FALSE.)
         ELSEIF(ITY.EQ.ITQ) THEN
            WRITE(CBUF,'(1X,A,1(A,E15.9))')
     +           CNAME(I),': MULTIPOLE, K1L=',SX
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,SY,.TRUE.)
         ELSEIF(ITY.EQ.IAS) THEN
            WRITE(CBUF,'(''!@ '',A,A,E13.7,3(A,F9.4))')
     +           CNAME(I),': SHIFTER, L=',RL,
     +           ', PHX=',SX,', PHY=',SY,', PHT=',TWIS(I)
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.IMA) THEN
            WRITE(CBUF,'(''!@ '',A,A,3(A,F9.4))')
     +           CNAME(I),': TWISSMAT',
     +           ', BETA=',SY,', ALPHA=',SX,', MU=',TWIS(I)
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.ISR) THEN
            WRITE(CBUF,
     +           '(''!@ '',A,A,E13.7,A,F9.4,4(A,F16.8),2(A,I6))')
     +           CNAME(I),': SRMCELL, L=',RL,
     +           ', ANGLE=',SY,', MUY=',SX,', EPS0=',STRZ(I),
     +           ', MUEPS=',TWIS(I),', OFFSET=',ISTR0(I),
     +           ', ORDER=',ISTR1(I)
            CALL EXTMAD(CBUF,ITY,IPH,ISH,IFH,0.D+0,.TRUE.)
         ELSEIF(ITY.EQ.I00) THEN
            CALL WARN('WRMAD0, can''t convert testelement '//
     +           'with itype=98')
         ENDIF

C========write to IUMAD0 (mad0 fmt, 80 cols, continuation, extensions...)
         CALL WRM80X(IUMAD0,CBUF)

 10   CONTINUE

      IF(BDIST) THEN
C=====find unique labels for lattice-elements
         DO 20 I=1,MALPH
            DO 120 J=1,MALPH
               CREG       = CALPH(I:I)//CALPH(J:J)
               BALPH(I,J) = .TRUE.
               
               IF(INDEX(CSEQNM(1:2),CREG).NE.0) THEN
                  BALPH(I,J) = .FALSE.
                  GOTO 120
               ENDIF
               
               DO 1120 K=1,NELE
                  IF(INDEX(CNAME(K)(1:2),CREG).NE.0) THEN
                     BALPH(I,J) = .FALSE.
                     GOTO 120
                  ENDIF
 1120          CONTINUE
               
 120        CONTINUE
 20      CONTINUE
         ENDIF

C=====sequence (lattice-list)
      IF(NSXP.EQ.NSUPER) THEN
         WRITE(IUMAD0,'(/,1X,4A)')CSEQNM,': SEQUENCE, REFER=',CREFER,
     +        ' !@, SUPER= 1'
      ELSE
         WRITE(IUMAD0,'(/,1X,4A,I3)')CSEQNM,': SEQUENCE, REFER=',CREFER,
     +        ' !@, SUPER=',NSUPER
      ENDIF

      EPOS = 0.D+0
      RPOS = 0.D+0
      NREG = 1
      IF(BDIST) THEN
         CREG = CREGMR(BALPH,NREG)
         IF(CREG.EQ.'!!') CALL TOFF(
     +        'WRMAD0, no free identifier even for first sector !')
      ENDIF


      DO 30 ISXP=1,NSXP
         WRITE(IUMAD0,'(2(A,I5))') '! SuperPeriod ',ISXP,' of ',NSUPER
         BSUP = ISXP.EQ.NSUPER
         
         DO 130 I=1,NLAT
            LE  = LATELE(I)
            DS  = DMAX1( 1.D+0 , DBLE(LSLIC(I)))
            DL  = PREDRF(I) + DREFER*RLEN(LE)*DS
            AT  = EPOS + DL
            AR  = RPOS + DL
            ITY = ITYPE(LE) 
            
            IF(ITY.EQ.IMR .AND. I.GT.1) THEN
               RPOS = 0.D+0
               AR   = 0.D+0
               NREG = NREG + 1
               IF(BDIST) THEN
                  CREG = CREGMR(BALPH,NREG)
                  IF(CREG.EQ.'  !!') CALL WARN( 'WRMAD0,'//
     +                 ' not enough free identifiers for all sectors')
               ENDIF
            ENDIF

            IF    (BDIST) THEN
               CRPOS  = CFINFL(NINT(DM2D*AR),'0')   
               CLATEL = CREG//CRPOS
            ELSEIF(BNUMB) THEN
               CCLASS = CNAME(LE)
               CLATEL = CCLASS(1:ICSEND(CCLASS))//CFINFL(NUMBER(LE),'_')
               NUMBER(LE) =  NUMBER(LE) + 1
            ENDIF

            IF(I.NE.NLAT .OR. BSUP) THEN
               IF(ITY.EQ.ISS .OR. ITY.EQ.ISF .OR. ITY.EQ.ISM .OR.
     +              ITY.EQ.IMA) THEN
                  WRITE(IUMAD0,'(''!@'',4A,F13.6,A,F10.6)')
     +                 CLATEL,': ',CNAME(LE),', AT=',AT,
     +                 ',    PRDRF=',PREDRF(I)
               ELSE
                  WRITE(IUMAD0,'(2X,4A,F13.6,A,F10.6)')
     +                 CLATEL,': ',CNAME(LE),', AT=',AT,
     +                 ' !@, PRDRF=',PREDRF(I)
               ENDIF
            ENDIF

            EPOS = EPOS + PREDRF(I) + RLEN(LE)*DS

            IF(ITY.NE.IMR .OR. I.EQ.1) THEN
               RPOS = RPOS + PREDRF(I) + RLEN(LE)*DS
            ENDIF
 130     CONTINUE
 30   CONTINUE

      WRITE(IUMAD0,'(1X,A)')'ENDSEQUENCE'

C=====done & bye
      CLOSE(IUMAD0)
      RETURN

C=====error while creating output unit
C=====from NEWVNR
 9991 CONTINUE
      IF(IERR.EQ.1) THEN
         CALL TOFF('WRMAD0, increase MAXNME')
      ELSEIF(IERR.EQ.2) THEN
         CALL TOFF('WRMAD0, increase MAXVNR')
      ENDIF
C=====from OPEN
 9993 CONTINUE
         WRITE(CDUM,'(A10)') IERR
         CALL TOFF('WRMAD0, system i/o error number : '//CDUM//
     +        ' while opening new slim-file')

      E N D

      S U B R O U T I N E WRSLI1(CSXPND)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  WRSLI1
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ COMMAND    :  
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     writes the  c u r r e n t  status of /ELEIN/ and /LATTICE/ into
C     CSLIM0.slim1.VNR
C     uses NEW  & c o n s i s t e n t  slim1 format with :
C     (general)
C     - 132 collumns
C     - no fixed form header but many (max 1000) comment lines before
C       'NSUPER = $NSUPER' ( <-- almost free form :
C             'NSUPER= n' ; 'NSUPER = n' ; 'NSUPER =n' ) 
C     - no 'IFORMAT'
C     (type list)
C     - strengths & lengths are stored with higher precision
C     - drifts have length and no strength
C     - element names are allowed char*(10) (as for internal representation)
C     - strengths are  i n t e r n a l y  divided by ISEC and don't have to be
C       scaled by hand 
C     - IFLAG is a bitarray (-> see parameter statement below !) and
C       is externally represented by a stream of 'T' and 'F' 
C     - IPS intruduced (P_ower S_upply number) for run-time scaling
C     (lattice list)
C     - position is given in meters (not tenths of mm) but with \mu m
C       precision
C     !!! write_slim1_lattice '$mylattice' $expand_superperiod; !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MSLCOL=132,MAXVNR=99)
      CHARACTER CSXPND*(*), CASEUP*5
      CHARACTER CSFIB*(MBITS)
      CHARACTER CNCHECK*(MAELNM), CN*(MAELNM), CNMHLP*(MAELNM), CDUM*10
      CHARACTER CLINE*(MSLCOL)
      DIMENSION DPOS(0:3), CN(0:3)

C%%      write(*,*)'wrsli1'

C=====preliminaries.......
      CALL INFO('WRSLI1')

      IUSLIM = 17

      IBXP = ICSBEG(CSXPND)
      IEXP = ICSEND(CSXPND)
      IF(IEXP-IBXP+1.GT.5) CALL TOFF('WRSLI1, logical must be'//
     +     ' ''true'' or ''false''')

      IF    (CASEUP(CSXPND(IBXP:IEXP)).EQ.'TRUE' .AND. NSUPER.GT.1) THEN
         WRITE(CDUM,'(I10)') NSUPER
         CALL WARN('WRSLI1, superperiod'//CDUM//' will be expanded')
         NSXP = NSUPER
      ELSEIF(CASEUP(CSXPND(IBXP:IEXP)).EQ.'FALSE' .OR. NSUPER.EQ.1) THEN
         NSXP = 1
      ELSE
         CALL TOFF('WRSLI1, logical must be ''true'' or ''false''')
      ENDIF
      
C=====slim-header
      WRITE(IUSLIM,'(A,A)') '-this file   : ',CLATOV(:ICSEND(CLATOV))
      WRITE(IUSLIM,'(A,A)') '-job name    : ',CJOB(:ICSEND(CJOB))
      WRITE(IUSLIM,'(A,A)') '-orig. latt. : ',CLATI(:ICSEND(CLATI))
      WRITE(IUSLIM,'(A,A)') '-WARNING: new (slim1) format !',
     +     ' drifts have strength 0 but nonvanishing length  !'
      WRITE(IUSLIM,'(A,A)') '-WARNING: new (slim1) format !',
     +     'strengths are  n o t  divided by slices !'
      CALL WARN(
     +     'WRSLI1, new (more consistent) ''slim1'' format :')
      CALL WARN(
     +     'WRSLI1, drifts have strength 0 but nonvanishing length')
      CALL WARN(
     +     'WRSLI1, strengths are  n o t  divided by slices')

      IF(NSXP.EQ.NSUPER) THEN
         WRITE(IUSLIM,'(A,I4)') ' NSUPER= ',1
      ELSE
         WRITE(IUSLIM,'(A,I4)') ' NSUPER= ',NSUPER
      ENDIF

C=====type-list
      
      DO 10 I=1,NELE
         ITP  = ITYPE (I)
         ITP2 = ITPEFF(I)

         IF(ITP.EQ.ISO .OR. ITP.EQ.IRF .OR.
     +        ITP.EQ.ISS .OR. ITP.EQ.ISF .OR.
     +        ITP.EQ.ISM .OR. ITP.EQ.IMR .OR. ITP.EQ.ISR) THEN
C========prevent certain elements from future slicing !!!!
            ISLHLP = 0
            RLEHLP = RLEN(I)
            STXHLP = STRX(I)
            STYHLP = STRY(I)
         ELSEIF(ISEC(I)*IGLBSL.GT.1) THEN
C========undo slicing
            ISLHLP = ISEC(I)*IGLBSL
            DSLHLP = DBLE(ISLHLP)
            RLEHLP = RLEN(I)*DSLHLP
            STXHLP = STRX(I)*DSLHLP
            STYHLP = STRY(I)*DSLHLP
         ELSE
            ISLHLP = 1
            RLEHLP = RLEN(I)
            STXHLP = STRX(I)
            STYHLP = STRY(I)
         ENDIF

         IF((ITP.EQ.ICX .OR. ITP.EQ.ICY) .AND. ITP2.EQ.IMQ) THEN
            STYHLP = STXHLP
            STXHLP = 0.D+0
         ENDIF

         IF(ITP.EQ.ISO) THEN
            STXHLP = STXHLP*RLEHLP
         ELSEIF(ITP.EQ.IRF) THEN
            STXHLP = STXHLP*1.D+3
         ENDIF

         CNMHLP = CNCHECK(CNAME(I),'ALL')

         WRITE(CLINE,
     +        '(1X,I2,1X,A10,1X,3(D16.9,1X),I3,1X,D10.3,1X,I3,1X,A)')
     +        ITP,CNMHLP,STXHLP,STYHLP,RLEHLP,
     +        ISLHLP,TWIS(I),IPS(I),CSFIB(IFLAG(I))
         WRITE(IUSLIM,'(A)') CLINE(1:ICSEND(CLINE))

         IF(ITP.EQ.ISR) WRITE(IUSLIM,
     +        '(1X,I2,1X,A10,1X,D16.9,2(1X,I16))')
     +        0,CNMHLP,STRZ(I),ISTR0(I),ISTR1(I)

 10   CONTINUE

C=====lattice-list
      IEND   = NUMELE('X______EOL')
      ENDPOS = 0.D+0

      DO 20 ISXP=1,NSXP
         WRITE(IUSLIM,'(2(A,I5))') '- SUPERPERIOD ',ISXP,' OF ',NSUPER
         BSUP = ISXP.EQ.NSUPER

         DO 120 I=1,NLAT,4
            JS = 3
            
            DO 1120 J=0,3
               LL      = I + J
               LT      = LATELE(LL)
               CN(J)   = CNCHECK(CNAME(LT),'ALL')
               DS      = DMAX1( 1.D+0 , DBLE(LSLIC(LL)) )
               RL      = RLEN(LT)*DS
               DPOS(J) = ENDPOS+PREDRF(LL)+RL*0.5D+0 
               ENDPOS  = ENDPOS+PREDRF(LL)+RL

               IF(LT.EQ.IEND) THEN
                  IF(BSUP) THEN
                     JS = J
                  ELSE
                     JS = J - 1
                  ENDIF
                  
                  GOTO 1121
               ENDIF
               
 1120       CONTINUE
 1121       CONTINUE
            
            WRITE(IUSLIM,'(X,4(A10,1X,F13.6,1X))')
     +           (CN(JW),DPOS(JW),JW=0,JS)
            
 120     CONTINUE
 20   CONTINUE

C=====done & bye
      CLOSE(IUSLIM)
      RETURN

C=====error while creating output unit
C=====from newvnr
 9991 CONTINUE
      IF(IERR.EQ.1) THEN
         CALL TOFF('WRSLI1, increase MAXNME')
      ELSEIF(IERR.EQ.2) THEN
         CALL TOFF('WRSLI1, increase MAXVNR')
      ENDIF
C=====from OPEN
 9993 CONTINUE
         WRITE(CDUM,'(A10)') IERR
         CALL TOFF('WRSLI1, system i/o error number : '//CDUM//
     +        ' while opening new slim-file')
      E N D

      S U B R O U T I N E WRSLI0(CSLFM0,CSXPND)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  WRSLI0
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ COMMAND    :  
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     writes the  c u r r e n t  status of /ELEIN/ and /LATTICE/ into
C     CSLIM0.slim0.VNR
C     THIS IS THE HISTORICALLY ORIGINAL FORMAT ! (as far as I know!)
C     It is as it is !
C     In case of questions as the inventors of SLIM and the keepers of
C     the code !
C     CSLFM0 = 'protons' or 'proton' =>
C              typelist format = '(1X,I4,1X,A8,3E14.6,I5,F11.6,I5)'
C              (D.B.Barber/protons)
C     CSLFM0 = 'electrons' or 'electron' =>
C              typelist format = '(1X,I4,1X,A8,3F12.8,I5,F11.6,I5)'
C              (D.B.Barber/original (possibly A.Chao ?)
C     !!! write_slim0_lattice '$mylattice' ${'protons','electrons'}
C                             $expand_superperiod ; !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MSLCOL=80,DM2LP = 1.D+4,MAXVNR=99)
      CHARACTER CNCHECK*(MAELNM), CN*(MAELNM), CNMHLP*(MAELNM),
     +     CDUM*10, CSLFM0*(*), CSLFMT*9, CTPLFM*32, CTPLF2*24,
     +     CSXPND*(*), CASEUP*5 
      DIMENSION RPOS(0:3), CN(0:3)

C%%      write(*,*)'wrsli0'

C=====preliminaries.......
      CALL INFO('WRSLI0')

      IUSLIM = 17

      IF(ICSEND(CSLFM0)-ICSBEG(CSLFM0).GT.9) CALL TOFF(
     +     'WRSLI0, invalid slim0-format qualifier')
      CSLFMT = CSLFM0(ICSBEG(CSLFM0):ICSEND(CSLFM0))
      
      IF    (CSLFMT.EQ.'protons' .OR. CSLFMT.EQ.'proton') THEN
C========D.P.Barber/protons :
         CTPLFM = '(1X,I4,1X,A8,3E14.6,I5,F11.6,I5)'
         CTPLF2 = '(1X,I4,1X,A8,E14.6,2I14)'
      ELSEIF(CSLFMT.EQ.'electrons' .OR. CSLFMT.EQ.'electron') THEN
C========D.P.Barber/original (possibly A.Chao ?!?!?!?) :
         CTPLFM = '(1X,I4,1X,A8,3F12.8,I5,F11.6,I5)'
         CTPLF2 = '(1X,I4,1X,A8,F12.8,2I12)'
      ELSE
         CALL TOFF('WRSLI0, invalid slim0-format qualifier')
      ENDIF

      IBXP = ICSBEG(CSXPND)
      IEXP = ICSEND(CSXPND)
      IF(IEXP-IBXP+1.GT.5) CALL TOFF('WRSLI0, logical must be'//
     +     ' ''true'' or ''false''')

      IF    (CASEUP(CSXPND(IBXP:IEXP)).EQ.'TRUE' .AND. NSUPER.GT.1) THEN
         WRITE(CDUM,'(I10)') NSUPER
         CALL WARN('WRSLI0, superperiod'//CDUM//' will be expanded')
         NSXP = NSUPER
      ELSEIF(CASEUP(CSXPND(IBXP:IEXP)).EQ.'FALSE' .OR. NSUPER.EQ.1) THEN
         NSXP = 1
      ELSE
         CALL TOFF('WRSLI0, logical must be ''true'' or ''false''')
      ENDIF

C=====slim-header
      IF(LSLIMF.GT.MSLCOL) CALL WARN('WRSLI0, slim-file-name '//
     +     'will be cut in header')
      WRITE(IUSLIM,'(A)') CLATO(:MIN(ICSEND(CLATO),MSLCOL))
      WRITE(IUSLIM,'(A,A)') '-job :',
     +     CJOB(:MIN(ICSEND(CJOB),MSLCOL-4))
      WRITE(IUSLIM,'(A,A)') '-original:',
     +     CLATI(:MIN(ICSEND(CLATI),MSLCOL-10))
      WRITE(IUSLIM,'(A)') 'IFORMAT= 4'
      WRITE(IUSLIM,'(A)') 'RNDMINDX 10'
      WRITE(IUSLIM,'(A,I3)') 'NSUPER=',NSUPER
      WRITE(IUSLIM,'(A,A)') '-type list format : ',CSLFMT
      WRITE(IUSLIM,'(A)')
     +     '-WARNING: old (slim0) format : '//
     +     'drifts have length 0 but strength = length  !'
      WRITE(IUSLIM,'(A)')
     +     '-WARNING: old (slim0 format) : '//
     +     'strengths a r e  divided by slices !'
      CALL WARN(
     +     'WRSLI0, old (historical) ''slim0'' format :') 
      CALL WARN(
     +     'WRSLI0, drifts have length 0 but strength = length')
      CALL WARN('WRSLI0, strengths  a r e  divided by slices')

C=====type-list
      
      DO 10 I=1,NELE
         ITP     = ITYPE (I)
         ITP2    = ITPEFF(I)
         IFLHLP  = 0

         IF(ISEC(I)*IGLBSL.GT.1) THEN
C========undo slicing for the length only
            ISLHLP = ISEC(I)*IGLBSL
            RLEHLP = RLEN(I)* DBLE(ISLHLP)
            STXHLP = STRX(I)
            STYHLP = STRY(I)
         ELSE
            ISLHLP = 1
            RLEHLP = RLEN(I)
            STXHLP = STRX(I)
            STYHLP = STRY(I)
         ENDIF

         IF(ITP.EQ.IDL) THEN
            STXHLP = RLEHLP
            RLEHLP = 0.D+0
         ELSEIF((ITP.EQ.ICX .OR. ITP.EQ.ICY) .AND. ITP2.EQ.IMQ) THEN
            STYHLP = STXHLP
            STXHLP = 0.D+0
         ELSEIF(ITP.EQ.ISO) THEN
            STXHLP = STXHLP*RLEHLP
         ELSEIF(ITP.EQ.IRF) THEN
            STXHLP = STXHLP*1.D+3
         ELSEIF(ITP.EQ.IMA) THEN
            IF    (BTEST(IFLAG(I),JBSPC7)) THEN
               IFLHLP = 1
            ELSEIF(BTEST(IFLAG(I),JBSPC8)) THEN
               IFLHLP = 2
            ELSEIF(BTEST(IFLAG(I),JBSPC9)) THEN
               IFLHLP = 3
            ENDIF
         ENDIF
         
         CNMHLP = CNCHECK(CNAME(I),'ALL')
         
         IF(ITP.EQ.IMR) ITP = IDL
   
         IF(ICSEND(CNMHLP).GT.8) THEN
            WRITE(IUSLIM,'(A,A10,A)')
     +           '-WARNING ',CNMHLP,' cut to 8 characters !'
            CALL WARN('WRSLI0, '//CNMHLP//' cut to 8 characters !')
         ENDIF

         WRITE(IUSLIM,CTPLFM)
     +        ITP,CNMHLP,STXHLP,STYHLP,RLEHLP,
     +        ISLHLP,TWIS(I),IFLHLP
         IF(ITP.EQ.ISR) WRITE(IUSLIM,CTPLF2)
     +        0,CNMHLP,STRZ(I),ISTR0(I),ISTR1(I)
 10   CONTINUE

C=====lattice-list
      IEND = NUMELE('X______EOL')
      ENDPOS = 0.D+0
      
      DO 20 ISXP=1,NSXP
         BSUP = ISXP.EQ.NSUPER
         
         DO 120 I=1,NLAT,4
            JS = 3
            
            DO 1120 J=0,3
               LL    = I+J
               LT    = LATELE(LL)
               CN(J) = CNCHECK(CNAME(LT),'ALL')
               DS    = DMAX1( 1.D+0 , DBLE(LSLIC(LL)) )
               IF(ICSEND(CN(J)).GT.8) CALL WARN(
     +              'SLATWR, '//CN(J)(:8)//' might be ambigous')
               RL      = RLEN(LT)*DS
               IF(ITYPE(LT).EQ.ISO) THEN 
                  RPOS(J) = ENDPOS+PREDRF(LL)
               ELSE
                  RPOS(J) = ENDPOS+PREDRF(LL)+RL*0.5D0
               ENDIF
               ENDPOS = ENDPOS+PREDRF(LL)+RL
               
               IF(LT.EQ.IEND) THEN
                  IF(BSUP) THEN
                     JS = J
                  ELSE
                     JS = J - 1
                  ENDIF
                  
                  GOTO 1121
               ENDIF
               
 1120       CONTINUE
 1121       CONTINUE
            
            WRITE(IUSLIM,'(4(A8,f16.8,1X))')
     +           (CN(JW),RPOS(JW),JW=0,JS)
            
 120     CONTINUE
 20   CONTINUE

C=====done & bye
      CLOSE(IUSLIM)
      RETURN

C=====error while creating output unit
C=====from newvnr
 9991 CONTINUE
      IF(IERR.EQ.1) THEN
         CALL TOFF('WRSLI0, increase MAXNME')
      ELSEIF(IERR.EQ.2) THEN
         CALL TOFF('WRSLI0, increase MAXVNR')
      ENDIF
C=====from OPEN
 9993 CONTINUE
         WRITE(CDUM,'(A10)') IERR
         CALL TOFF('WRSLI0, system i/o error number : '//CDUM//
     +        ' while opening new slim-file')
      E N D

      S U B R O U T I N E WRFOXY
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  WRFOXY
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ COMMAND    :  
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     writes the  c u r r e n t  status of /ELEIN/ and /LATTICE/ into
C     CFOXY.fox.VNR
C     - This file includes the COSY.fox soubroutine library and
C       can directly be run by FOXY.
C     - strengths & lengths are stored with higher precision
C     - IFLAG is a bitarray (-> see parameter statement below !) and
C       is externally represented by a stream of 'T' and 'F' 
C       it is passed to the foxy file as a variable and can be used later.
C     - IPS intruduced (P_ower S_upply number) for run-time scaling
C       (lattice list) The scaling is introduced in the propper
C       place depending on which iflag mode is used.
C     - All number are converted to SI units
C     !!! write_foxy_lattice '$mylattice' ; !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MAXVNR=99)
      CHARACTER CDUM*10

C%%      write(*,*)'wrfoxy'

C=====preliminaries.......
      CALL INFO('WRFOXY')

      IUFOXY = 17

C=====foxy-header
      WRITE(IUFOXY,'(A)') '{'
      WRITE(IUFOXY,'(A,A)') ' this file   : ',CLATO(:ICSEND(CLATO))
      WRITE(IUFOXY,'(A,A)') ' job name    : ',CJOB(:ICSEND(CJOB))
      WRITE(IUFOXY,'(A,A)') ' orig. latt. : ',CLATI(:ICSEND(CLATI))
      WRITE(IUFOXY,'(A)') '}'

      CALL FOXYBE(IUFOXY)

C=====type-list
      CALL FOXYELEOUT(IUFOXY)

C=====lattice-list
      CALL FOXYLATOUT(IUFOXY)

      CALL FOXYEN(IUFOXY)

C=====done & bye
      CLOSE(IUFOXY)
      RETURN

C=====error while creating output unit
C=====from newvnr
 9991 CONTINUE
      IF(IERR.EQ.1) THEN
         CALL TOFF('WRFOXY, increase MAXNME')
      ELSEIF(IERR.EQ.2) THEN
         CALL TOFF('WRFOXY, increase MAXVNR')
      ENDIF
C=====from OPEN
 9993 CONTINUE
         WRITE(CDUM,'(A10)') IERR
         CALL TOFF('WRFOXY, system i/o error number : '//CDUM//
     +        ' while opening new slim-file')
      E N D

      S U B R O U T I N E WRPETROS(CPAFM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  WRPETROS
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ COMMAND    :  
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     writes the  c u r r e n t  status of /ELEIN/ and /LATTICE/ into
C     CPETROS.pet.VNR
C     - This file contains the total lattice, even when a super-period
C       is used.  This is done, since PETROS officially does not support
C       super-periods.
C     - strengths & lengths are stored with higher precision
C     - All number are converted to SI units
C     CPAFM = 'protons' or 'proton' => PXCU=PSTR  etc. for dipoles
C     CPAFM = 'electrons' or 'electron' => PXCU=0 etc. for dipoles
C     !!! write_petros_lattice '$mylattice' ; !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MAXVNR=99)
      CHARACTER CPETROS1*(MAXNME),CPETROSF*(MAXNME+4),
     +     CDUM*10,CPAFM*(*)

C%%      write(*,*)'wrpetros'

C=====preliminaries.......
      CALL INFO('WRPETROS')

      IUPETROS = 17

      BEL = CPAFM.EQ.'electron'.OR.CPAFM.EQ.'electrons'
      BPR = CPAFM.EQ.'proton'.OR.CPAFM.EQ.'protons'
      IF(.NOT.(BEL.OR.BPR))
     +     CALL TOFF('WRPETROS, invalid particle format')

C=====petros-header
      CALL PETROSBE(IUPETROS)

C      WRITE(IUPETROS,'(A,A)') '# this file   : ',CPETROSF(:LPETROSF)
C      WRITE(IUPETROS,'(A,A)') '# sprint vers.: ',CSHDPR(:ICSEND(CSHDPR))
C      WRITE(IUPETROS,'(A,A)') '# inp. header : ',CSHDIP(:ICSEND(CSHDIP))
C      WRITE(IUPETROS,'(A,A)') '# sprint run  : ',CSHDDT(:ICSEND(CSHDDT))
C      WRITE(IUPETROS,'(A,A)') '# job name    : ',COUTF1(:ICSEND(COUTF1))
C      WRITE(IUPETROS,'(A,A)') '# orig. latt. : ',CLATF(:ICSEND(CLATF))

C=====type-list
      CALL PETROSELEOUT(IUPETROS,BPR)

C=====lattice-list
      CALL PETROSLATOUT(IUPETROS)

      CALL PETROSEN(IUPETROS)

C=====done & bye
      CLOSE(IUPETROS)
      RETURN

C=====error while creating output unit
C=====from newvnr
 9991 CONTINUE
      IF(IERR.EQ.1) THEN
         CALL TOFF('WRPETROS, increase MAXNME')
      ELSEIF(IERR.EQ.2) THEN
         CALL TOFF('WRPETROS, increase MAXVNR')
      ENDIF
C=====from OPEN
 9993 CONTINUE
         WRITE(CDUM,'(A10)') IERR
         CALL TOFF('WRPETROS, system i/o error number : '//CDUM//
     +        ' while opening new PETROS-file')
      E N D


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C15                                                                           C
C@ PROGRAM=SPRINT, MODULE=LATTICE, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      S U B R O U T I N E LATIN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  LATIN
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ EXTERNAL   :  DEFSLIM,LAYSLIM,INFCS,ICSEND,TOFF,DEFPET,LAYPET,RDMAD0,RDMDTW
C@ IO_UNITS   :  15   : lattice file (read), allocated by MAIN
C-                       constant unit number (15) will be replaced by
C-                       dynamical unit (-> NEWUNT) in  later releases !!!
C-               16   : logfile (write)
C@ DESCRIPTION:  Determines flavour of lattice definition language from the
C-                extention of lattice definition file (neam.extension.vnr)
C-                in CLATF.
C-                 .mad0.               MAD format w/o block structure
C-                                          (RECOMMEMENDED)
C-                 .slim0.              old SLIM format,     
C-                 .slim1.              new SLIM format      
C-                 .twiss. or .survey.  MAD output
C-                 .pet.                PETROS format (NOT RECOMMENDED),
C-               Then reads headers of lattice-file CLATF and calls appropriate
C-                file-format-reader
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MLIN=132)
      CHARACTER CLINE*(MLIN),CPOS*4,CNUM*10
      DATA CNUM /'0123456789'/

C%%      write(*,*)'latin'

      CALL INFO('LATIN')

      INSL00 = INDEX(CLATI,'.slim0')
      INSL10 = INDEX(CLATI,'.slim1')
      INPET0 = INDEX(CLATI,'.pet')
      INMAD0 = INDEX(CLATI,'.mad0')
      INTWI0 = INDEX(CLATI,'.twiss')
      INSUR0 = INDEX(CLATI,'.survey')
      INSL0V = INDEX(CLATI,'.slim0.')
      INSL1V = INDEX(CLATI,'.slim1.')
      INPETV = INDEX(CLATI,'.pet.')
      INMADV = INDEX(CLATI,'.mad0.')
      INTWIV = INDEX(CLATI,'.twiss.')
      INSURV = INDEX(CLATI,'.survey.')

      IF(INSL0V+INSL1V+INPETV+INMADV+INTWIV+INSURV.NE.0) THEN
         DO 10 I=ICSEND(CLATI),1,-1
            IF(INDEX(CNUM,CLATI(I:I)).EQ.0) THEN
               LCLATI = I
               GOTO 11
            ENDIF
 10      CONTINUE
 11      CONTINUE
      ELSE
         LCLATI = ICSEND(CLATI)
      ENDIF


      IF(INSL00.GT.0 .AND.
     +     (INSL00.EQ.LCLATI-5 .OR. INSL0V.EQ.LCLATI-6)) THEN
C========old (slim0) format
         READ(15,'(A)') CLINE
          WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
         READ(15,'(A)') CLINE
          WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
         READ(15,'(A)') CLINE
          WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
         READ(15,'(A)') CLINE
          WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
         READ(15,'(A)') CLINE
          WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
C     IFORMAT = INFCS(CPOS(CLINE,2))
         READ(15,'(A)') CLINE
          WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
         NSUPER = INFCS(CPOS(CLINE,2))
         CALL DEFSLIM(.FALSE.)
         CALL LAYSLIM(.FALSE.)
      ELSEIF(INSL10.GT.0 .AND.
     +     (INSL10.EQ.LCLATI-5 .OR. INSL1V.EQ.LCLATI-6)) THEN
C========new (slim1) format
         CALL DEFSLIM(.TRUE.)
         CALL LAYSLIM(.TRUE.)
      ELSEIF(INPET0.GT.0 .AND.
     +        (INPET0.EQ.LCLATI-3 .OR. INPETV.EQ.LCLATI-4)) THEN
C========petros format
         READ(15,'(A)') CLINE
          WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
         READ(15,'(A)') CLINE
          WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
         READ(15,'(A)') CLINE
          WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
         READ(15,'(A)') CLINE
          WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
         I = INDEX(CLINE,'NSUPER=')
C
C        NSUPER=### is here allowed in PETROS, although
C        this is not part of the PETROS standard (if there is such a
C        thing at all). NSUPER=### has to be in the fourth line of
C        the lattice file and there can not be a space between NSUPER and =.
C
         IF(I.NE.0) THEN
            J = INDEX(CLINE(I+7:),' ')
            IF(J.GE.2) THEN
               NSUPER = INFCS(CLINE(I+7:I+7+J-1))
            ELSE
               CALL TOFF('LATIN, NSUPER= without value in .pet file')
            ENDIF
         ELSE
            NSUPER = 1
         ENDIF
         CALL DEFPET
         CALL LAYPET
      ELSEIF(INMAD0.GT.0 .AND.
     +        (INMAD0.EQ.LCLATI-4 .OR. INMADV.EQ.LCLATI-5)) THEN
C========mad0 format
         CALL RDMAD0
      ELSEIF(INTWI0+INSUR0.GT.0 .AND.
     +        (INTWI0.EQ.LCLATI-5 .OR. INTWIV.EQ.LCLATI-6 .OR.
     +        INSUR0.EQ.LCLATI-6 .OR. INSURV.EQ.LCLATI-7)) THEN
C========mad-output (twiss & survay) 
         CALL RDMDTW
      ELSE
C========well folks, don't you think that 6 lattice formats suffice ?
C========everybody oposing is invited to include her/his vavourite format !
         CALL TOFF('LATIN, unknown format')
      ENDIF

      RETURN
      E N D

      S U B R O U T I N E LATOUT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  LATIN
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ EXTERNAL   :  WRMAD0, WRSLI1, WRSLI0, WRFOXY,  WRPETROS
C@ IO_UNITS   :  17   : lattice file (read), allocated by MAIN
C-                       constant unit number (15) will be replaced by
C-                       dynamical unit (-> NEWUNT) in  later releases !!!
C-               16   : logfile (write)
C@ DESCRIPTION:  Determines flavour of lattice definition language from the
C-                extention of lattice definition file (neam.extension.vnr)
C-                in CLATO.
C-                 .mad0.               MAD format w/o block structure
C-                                          (RECOMMEMENDED)
C-                 .slim0.              old SLIM format,     
C-                 .slim1.              new SLIM format      
C-                 .twiss. or .survey.  MAD output
C-                 .pet.                PETROS format (NOT RECOMMENDED),
C-               Then reads headers of lattice-file CLATF and calls appropriate
C-                file-format-reader
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MLIN=132)
      CHARACTER CLINE*(MLIN),CPOS*4,CNUM*10
      DATA CNUM /'0123456789'/

C%%      write(*,*)'latout'

      INSL00 = INDEX(CLATO,'.slim0')
      INSL10 = INDEX(CLATO,'.slim1')
      INPET0 = INDEX(CLATO,'.pet')
      INMAD0 = INDEX(CLATO,'.mad0')
      INFOX0 = INDEX(CLATO,'.fox')

      LCLATO = ICSEND(CLATO)


      IF(INSL00.GT.0 .AND. INSL00.EQ.LCLATO-5) THEN
C========old (slim0) format
         CALL WRSLI0( COPT(1),COPT(2) )
      ELSEIF(INSL10.GT.0 .AND. INSL10.EQ.LCLATO-5) THEN
C========new (slim1) format
         CALL WRSLI1( COPT(1) )
      ELSEIF(INPET0.GT.0 .AND. INPET0.EQ.LCLATO-3) THEN
C========petros format
         CALL WRPETROS( COPT(1) )
      ELSEIF(INMAD0.GT.0 .AND. INMAD0.EQ.LCLATO-4) THEN
C========mad0 format
         CALL WRMAD0( COPT(1),COPT(2),COPT(3),COPT(4) )
      ELSEIF(INFOX0+INSUR0.GT.0 .AND. INFOX0.EQ.LCLATO-3) THEN
C========cosy format 
         CALL WRFOXY
      ELSE
C========well folks, don't you think that 6 lattice formats suffice ?
C========everybody oposing is invited to include her/his vavourite format !
         CALL TOFF('LATOUT, unknown format')
      ENDIF

      RETURN
      E N D

      INTEGER F U N C T I O N NUMELE(CWORD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  NUMELE
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     returns the position of 'CWORD' in the
C     element-list

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====local stuff
      CHARACTER*(*) CWORD

C%%      write(*,*)'numele'

      DO 100 I=1,NELE
         IF(CWORD.EQ.CNAME(I)) THEN
            NUMELE = I
            RETURN
         ENDIF
 100  CONTINUE

      NUMELE = 0

      RETURN
      E N D

      CHARACTER*(*) F U N C T I O N CNCHECK(CWORD,CFORM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  CNCHECK
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Changes unpleasant names to more common ones.
C
C     Input Parameters:
C     CWORD        : Element name which should be checked and maybe changed
C                    if it does not fit in with the format in which this
C                    element should be written.
C     CFORM        : Descibes the format in which this element name should be
C                    written.
C           (FOXY) : FOXY output format
C           (ALL)  : everything else EXECPT MAD0 !!!!!!!!!!!!!!!
C
C     Output Parameters:
C     CNCHECK      : Agrees with CWORD, if this element name can be used in
C                    the format CFORM.  If not, the name will have changed.

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(NFOXY=2,NALL=1)
      CHARACTER CWORD*(*), CFORM*(*)
      CHARACTER COFOXY*10, CNFOXY*10, COALL*10, CNALL*10
      DIMENSION COFOXY(NFOXY), CNFOXY(NFOXY),
     +     COALL(NALL), CNALL(NALL)
      DATA COFOXY / 'RF','X______EOL' /
      DATA CNFOXY / 'HF','ENDOL' /
      DATA COALL / 'X______EOL' /
      DATA CNALL / 'END' /

C%%      write(*,*)'cncheck'

      IF(CFORM.EQ.'FOXY') THEN
         DO 10 I=1,NFOXY
            IF(CWORD.EQ.COFOXY(I)) THEN
               CNCHECK = CNFOXY(I)
               RETURN
            ENDIF
 10      CONTINUE
      ELSEIF(CFORM.EQ.'ALL') THEN
         DO 20 I=1,NALL
            IF(CWORD.EQ.COALL(I)) THEN
               CNCHECK = CNALL(I)
               RETURN
            ENDIF
 20      CONTINUE
      ELSE
         CALL TOFF('CNCHECK, illegal format')
      ENDIF

      CNCHECK = CWORD

      RETURN
      E N D


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C25                                                                           C
C@ PROGRAM=SPRINT, MODULE=COSY, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      S U B R O U T I N E FOXYBE(IUFOXY)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  FOXYBE
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C%%      write(*,*)'foxybe'

      WRITE(IUFOXY,'(A)')
     +     'include ''/afs/desy.de/user/h/hoff/ext/COSY'' ;'
      WRITE(IUFOXY,'(3(A/)A)')
     +     '   procedure run ;        variable isuper 1 ;',
     +     '      variable gevp0  1 ; variable gevm0  1 ;',
     +     '      variable eleco  1 ; variable aspin  1 ;',
     +     '      variable bspin  1 ; variable nsuper 1 ;'
      WRITE(IUFOXY,'(A/3(A,E25.17,A/),A,E25.17,A)')
     +     '      procedure data ;',
     +     '         gevp0  :=',GEVP0,' ;',
     +     '         gevm0  :=',GEVM0,' ;',
     +     '         eleco  :=',ELECO,' ;',
     +     '         aspin  :=',ASPIN,' ;'
      IF(BSPIN) THEN
         WRITE(IUFOXY,'(A)')
     +     '         bspin  :=                     true ;'
      ELSE
         WRITE(IUFOXY,'(A)')
     +     '         bspin  :=                not(true) ;'
      ENDIF
      WRITE(IUFOXY,'(A,I25,A)')
     +     '         nsuper :=',NSUPER,' ; endprocedure ;'

      RETURN
      E N D

      S U B R O U T I N E FOXYELEOUT(IUFOXY)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  FOXYELEOUT
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      CHARACTER CNCHECK*(MAELNM),CN*(MAELNM)

C%%      write(*,*)'foxyeleout'

      DO 100 I=1,NELE

         ITP = ITPEFF(I)

C========undo slicing
         DSL = 1.D+0
         IF(ITP.NE.ISO .AND. ITP.NE.IRF .AND. ITP.NE.ISS .AND.
     +      ITP.NE.ISF .AND. ITP.NE.ISM .AND. ITP.NE.IMR) THEN
            J = ISEC(I)*IGLBSL
            IF(J.GT.1) DSL = DBLE(J)
         ENDIF
         
         CN = CNCHECK(CNAME(I),'FOXY')

         WRITE(IUFOXY,'(6X,A10,A10,A2)') 'procedure ',CN,' ;'
         CALL COSYELE(IUFOXY,I,CN,DSL*RLEN(I),DSL*STRX(I),DSL*STRY(I))
         WRITE(IUFOXY,'(9X,A)') 'endprocedure ;'
 100  CONTINUE

      RETURN
      E N D

      S U B R O U T I N E COSYELE(IUFOXY,IELE,CNAMEI,RLENI,STRXI,STRYI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  COSYELE
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      CHARACTER CNAMEI*(*)

C%%      write(*,*)'cosyele'

      APERTURE = 0.1D+0
      ITP      = ITYPE (IELE)
      ISECI    = ISEC  (IELE)
      TWISI    = TWIS  (IELE)
      IFLAGI   = IFLAG (IELE)
      IPSI     = IPS   (IELE)
C     Warning: IFLAG and IPS are not includend in the FOXY file.
C              They can later be included if needed.

*******************************************
*     Drift                               *
*******************************************
      IF    (ITP.EQ.IDL) THEN
         WRITE(IUFOXY,'(9X,A,E17.9,A)') 'dl',RLENI,' ;'
*******************************************
*     Horizontal Wedge Dipole             *
*******************************************
      ELSEIF(ITP.EQ.IBX) THEN
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3,4I2,A)')
     +        'dil',RLENI,DABS(STRXI),'/degrad',APERTURE,0,0,0,0,' ;'
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
*******************************************
*     Quadrupole                          *
*******************************************
      ELSEIF(ITP.EQ.IMQ) THEN
         WRITE(IUFOXY,'(9X,A,2E17.9,A,SP,E15.9,SS,A,E11.3,A)')
     +        'mqk',RLENI,STRXI,'/(',RLENI,')',APERTURE,' ;'
*******************************************
*     Skew Quadrupole                     *
*******************************************
      ELSEIF(ITP.EQ.ISQ) THEN
         WRITE(IUFOXY,'(9X,A)') 'ra -45 ;'
         WRITE(IUFOXY,'(9X,A,2E17.9,A,SP,E15.9,SS,A,E11.3,A)')
     +        'mqk',RLENI,STRXI,'/(',RLENI,')',APERTURE,' ;'
         WRITE(IUFOXY,'(9X,A)') 'ra  45 ;'
*******************************************
*     Cavity                              *
*******************************************
      ELSEIF(ITP.EQ.IRF) THEN
         WRITE(IUFOXY,'(9X,A,9X,A,E17.9,A)') 'variable kvolt 1 1 1 ;',
     +        'kvolt(1,1) := ',STRXI*1.D+6,' ;'
         IF(RLENI.NE.0.D+0) WRITE(IUFOXY,'(9X,A2,E17.9,A)')
     +        'dl',RLENI,'/2 ;'
         WRITE(IUFOXY,'(9X,A,E17.9,A,E11.3,A)')
     +        'rf kvolt 0 ',STRYI,' 0 ',APERTURE,' ;'
         IF(RLENI.NE.0.D+0)
     +        WRITE(IUFOXY,'(9X,A,E17.9,A)') 'dl',RLENI,'/2 ;'
*******************************************
*     Horizontal Kicker                   *
*******************************************
      ELSEIF(ITP.EQ.IKX) THEN
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3,4I2,A)')
     +        'dil',RLENI,DABS(STRXI),'/degrad',APERTURE,0,0,0,0,' ;'
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         IF(STRXI.NE.0.D+0) THEN
            WRITE(IUFOXY,'(9X,A,E17.9,A)') 'ta',STRXI,'/degrad ;'
            WRITE(IUFOXY,'(9X,A,E15.9,A,SP,E15.9,A,E15.9,SS,A)')
     +           'sa ',RLENI,'*(1-cos(',STRXI,'))/(',STRXI,') 0 ;'
            WRITE(IUFOXY,'(9X,A,E15.9,A,SP,E15.9,A,E15.9,SS,A)')
     +           'dl',RLENI,'*(1-sin(',STRXI,')/(',STRXI,')) ;'
         ENDIF
*******************************************
*     Vertical Kicker                     *
*******************************************
      ELSEIF(ITP.EQ.IKY) THEN
         WRITE(IUFOXY,'(9X,A)') 'ra -90 ;'
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3,4I2,A)')
     +        'dil',RLENI,DABS(STRXI),'/degrad',APERTURE,0,0,0,0,' ;'
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         IF(STRXI.NE.0.D+0) THEN
            WRITE(IUFOXY,'(9X,A,E17.9,A)') 'ta',STRXI,'/degrad ;'
            WRITE(IUFOXY,'(9X,A,E15.9,A,SP,E15.9,A,E15.9,SS,A)')
     +           'sa ',RLENI,'*(1-cos(',STRXI,'))/(',STRXI,') 0 ;'
            WRITE(IUFOXY,'(9X,A,E15.9,A,SP,E15.9,A,E15.9,SS,A)')
     +           'dl',RLENI,'*(1-sin(',STRXI,')/(',STRXI,')) ;'
         ENDIF
         WRITE(IUFOXY,'(9X,A)') 'ra  90 ;'
*******************************************
*     Sextupole                           *
*******************************************
      ELSEIF(ITP.EQ.IMS) THEN
         WRITE(IUFOXY,'(9X,A,2E17.9,A,SP,E15.9,SS,A,E11.3,A)')
     +        'msk',RLENI,STRXI,'/(',RLENI,'*2)',APERTURE,' ;'
*******************************************
*     Vertical Bend                       *
*******************************************
      ELSEIF(ITP.EQ.IBY) THEN
         WRITE(IUFOXY,'(9X,A)') 'ra -90 ;'
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3,4I2,A)')
     +        'dil',RLENI,DABS(STRXI),'/degrad',APERTURE,0,0,0,0,' ;'
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         WRITE(IUFOXY,'(9X,A)') 'ra  90 ;'
*******************************************
*     Solenoid                            *
*******************************************
      ELSEIF(ITP.EQ.ISO) THEN
         WRITE(IUFOXY,'(9X,A,E17.9,A,SP,E15.9,SS,A,E11.3,E17.9,A)')
     +        'cmsk',STRXI,'/(',RLENI,')',APERTURE,RLENI,' ;'
*******************************************
*     Rotated Quadrupole                  *
*******************************************
      ELSEIF(ITP.EQ.IRQ) THEN
         WRITE(IUFOXY,'(9X,A,E17.9,A)') 'ra',-TWISI,'/degrad ;'
         WRITE(IUFOXY,'(9X,A,2E17.9,A,SP,E15.9,SS,A,E11.3,A)')
     +        'mqk',RLENI,STRXI,'/(',RLENI,')',APERTURE,' ;'
         WRITE(IUFOXY,'(9X,A,E17.9,A)') 'ra',TWISI,'/degrad ;'
*******************************************
*     Horizontal rectangular bend         * 
*******************************************
      ELSEIF(ITP.EQ.IRX) THEN
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3/10X,2(E17.9,A,I2),A)')
     +        'dil',RLENI,DABS(STRXI),'/degrad',APERTURE,
     +        DABS(STRXI),'/degrad/2',0,DABS(STRXI),'/degrad/2',0,' ;'
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
*******************************************
*     Vertical rectangular bend           *
*******************************************
      ELSEIF(ITP.EQ.IRY) THEN
         WRITE(IUFOXY,'(9X,A)') 'ra -90 ;'
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3/10X,2(E17.9,A,I2),A)')
     +        'dil',RLENI,DABS(STRXI),'/degrad',APERTURE,
     +        DABS(STRXI),'/degrad/2',0,DABS(STRXI),'/degrad/2',0,' ;'
         IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         WRITE(IUFOXY,'(9X,A)') 'ra  90 ;'
*******************************************
*     Instantanious spin match            *
*******************************************
      ELSEIF(ITP.EQ.ISM) THEN
         CALL TOFF('COSYELE, ism not implemented')
*******************************************
*     Reads Map From File                 *
*******************************************
      ELSEIF(ITP.EQ.IRO) THEN
         CALL TOFF('COSYELE, type 14 needs a stored map')
*******************************************************
*     Horizontal Combined Function Rectangular Magnet *
*******************************************************
      ELSEIF(ITP.EQ.ICH) THEN
         IF(STRXI.NE.0.D+0) THEN
            IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
            WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3/'//
     +           '9X,E17.9,A,SP,E15.9,A,E15.9,SS,A)')
     +           'mpl',RLENI,DABS(STRXI),'/degrad',APERTURE,
     +           -STRYI,'*(',RLENI,')/(',DABS(STRXI),')^2 0 0 0 0 ;'
            IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         ELSEIF(STRYI.NE.0.D+0) THEN
            WRITE(IUFOXY,'(9X,A,2E17.9,A,SP,E15.9,SS,A,E11.3,A)')
     +           'mqk',RLENI,STRYI,'/(',RLENI,')',APERTURE,' ;'
         ELSE
            WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3,E17.9,4I2,A)')
     +           'mpl',RLENI,0.d+0,'/degrad',APERTURE,
     +           0.d+0,0,0,0,0,' ;'
         ENDIF
**************************************************
*     Horizontal Combined Function Sector Magnet *
**************************************************
      ELSEIF(ITP.EQ.ICX) THEN
         IF(STRXI.NE.0.D+0) THEN
            IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
            WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3/'//
     +           '9X,E17.9,A,SP,E15.9,A,E15.9,SS,A)')
     +           'msl',RLENI,DABS(STRXI),'/degrad',APERTURE,
     +           -STRYI,'*(',RLENI,')/(',DABS(STRXI),')^2 0 0 0 0 ;'
            IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         ELSEIF(STRYI.NE.0.D+0) THEN
            WRITE(IUFOXY,'(9X,A,2E17.9,A,SP,E15.9,SS,A,E11.3,A)')
     +           'mqk',RLENI,STRYI,'/(',RLENI,')',APERTURE,' ;'
         ELSE
            WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3,E17.9,4I2,A)')
     +           'msl',RLENI,0.d+0,'/degrad',APERTURE,
     +           0.d+0,0,0,0,0,' ;'
         ENDIF
*****************************************************
*     Vertical Combined Function Rectangular Magnet *
*****************************************************
      ELSEIF(ITP.EQ.ICV) THEN
         WRITE(IUFOXY,'(9X,A)') 'ra -90 ;'
         IF(STRXI.NE.0.D+0) THEN
            IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
            WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3/'//
     +           '9X,E17.9,A,SP,E15.9,A,E15.9,SS,A)')
     +           'mpl',RLENI,DABS(STRXI),'/degrad',APERTURE,
     +           STRYI,'*(',RLENI,')/(',DABS(STRXI),')^2 0 0 0 0 ;'
            IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         ELSEIF(STRYI.NE.0.D+0) THEN
            WRITE(IUFOXY,'(9X,A,2E17.9,A,SP,E15.9,SS,A,E11.3,A)')
     +           'mqk',RLENI,-STRYI,'/(',RLENI,')',APERTURE,' ;'
         ELSE
            WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3/9X,E17.9,A)')
     +           'mpl',RLENI,DABS(STRXI),'/degrad',APERTURE,
     +           0.d+0,' 0 0 0 0 ;'
         ENDIF
         WRITE(IUFOXY,'(9X,A)') 'ra  90 ;'
************************************************
*     Vertical Combined Function Sector Magnet *
************************************************
      ELSEIF(ITP.EQ.ICY) THEN
         WRITE(IUFOXY,'(9X,A)') 'ra -90 ;'
         IF(STRXI.NE.0.D+0) THEN
            IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
            WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3/'//
     +           '9X,E17.9,A,SP,E15.9,A,E15.9,SS,A)')
     +           'msl',RLENI,DABS(STRXI),'/degrad',APERTURE,
     +           STRYI,'*(',RLENI,')/(',DABS(STRXI),')^2 0 0 0 0 ;'
            IF(STRXI.LT.0) WRITE(IUFOXY,'(9X,A)') 'cb ;'
         ELSEIF(STRYI.NE.0.D+0) THEN
            WRITE(IUFOXY,'(9X,A,2E17.9,A,SP,E15.9,SS,A,E11.3,A)')
     +           'mqk',RLENI,-STRYI,'/(',RLENI,')',APERTURE,' ;'
         ELSE
            WRITE(IUFOXY,'(9X,A,2E17.9,A,E11.3/9X,E17.9,A)')
     +           'msl',RLENI,DABS(STRXI),'/degrad',APERTURE,
     +           0.d+0,' 0 0 0 0 ;'
         ENDIF
         WRITE(IUFOXY,'(9X,A)') 'ra  90 ;'
*******************************************
*     Snake                               *
*******************************************
      ELSEIF(ITP.EQ.ISS) THEN
         WRITE(IUFOXY,'(9X,A,4E15.7,A)')
     +        'snake',RLENI,STRYI,STRXI,TWISI,' ;'
*******************************************
*     Free Snake                          *
*******************************************
      ELSEIF(ITP.EQ.ISF) THEN
         WRITE(IUFOXY,'(9X,A/9X,A/9X,A)')
     +        'variable t 1 ; variable p 1 ; variable a 1 ;',
     +        't = 0 ; p = 0 ; a = 0 ;',
     +        'snake 0 t p a ;'
*******************************************
*     Instantanious spin match            *
*******************************************
      ELSEIF(ITP.EQ.ISM) THEN
         WRITE(IUFOXY,'(9X,A,E15.7,A)')
     +        'dl',RLENI,' ; {Instantanious spin match}'
*******************************************
*     Marker                              *
*******************************************
      ELSEIF(ITP.EQ.IMR) THEN
         WRITE(IUFOXY,'(9X,A)')
     +        'dl 0 ; {Marker}'
*******************************************
*     Test element                        *
*******************************************
      ELSEIF(ITP.EQ.IMR) THEN
         WRITE(IUFOXY,'(9X,A)')
     +        'dl 0 ; {Test element}'
*******************************************
*     Unknown type of element             *
*******************************************
      ELSE
         WRITE(IUFOXY,'(9X,A,I6,A,A,A,E9.2,A)')
     +        '{unknown element type ',ITP,' - ',
     +        CNAMEI(:ICSEND(CNAMEI)),' with strenght ',STRXI,'}'
         WRITE(IUFOXY,'(9X,A)') 'dl 0 ;'
      ENDIF

      RETURN
      E N D

      S U B R O U T I N E FOXYLATOUT(IUFOXY)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  FOXYLATOUT
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      CHARACTER CNCHECK*(MAELNM)

C%%      write(*,*)'foxylatout'

      WRITE(IUFOXY,'(6X,A)') 'procedure lattice ;'

      DO 100 I=1,NLAT
         WRITE(IUFOXY,'(9X,A,E25.17,A,A,A)')
     +        'dl',PREDRF(I),' ; ',
     +        CNCHECK(CNAME(LATELE(I)),'FOXY'),' ;'
 100  CONTINUE

      WRITE(IUFOXY,'(9X,A)') 'endprocedure ;'

      RETURN
      E N D

      S U B R O U T I N E FOXYEN(IUFOXY)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  FOXYEN
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     This subroutine creates the end and execution part of the
C     FOXY program.

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C%%      write(*,*)'foxyen'

      WRITE(IUFOXY,'(6X,2A/6X,A))')
     +     'data ; ov 1 3 0 ; ',
     +     'rpm gevp0*1000 gevm0*1000/amumev eleco ;',
     +     'if aspin#0 ; spn 1 aspin ; endif ;'
      IF(NSUPER.GT.1) THEN
         WRITE(IUFOXY,'(6X,A,I4,A)')
     +        'um ; loop isuper 1 ',NSUPER,' ; lattice ; endloop ;'
      ELSE
         WRITE(IUFOXY,'(6X,A)') 'um ; lattice ; wm 6 ;'
      ENDIF
      WRITE(IUFOXY,'(6X,A)') 'endprocedure ; run ; end ;'

      RETURN
      E N D


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C26                                                                           C
C@ PROGRAM=SPRINT, MODULE=SLIM, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      S U B R O U T I N E DEFSLIM(BNEW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  DEFSLIM
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads data from 15, loads '/ELEIN/' for creation
C     of 'TMAP'-array in sub. 'CREAMAP'.
C     The element definition part has to end with the first
C     occurence of a marker named 'END'.  Internally this marker is
C     renamed to 'X______EOL'.
C     Input Parameters:
C     BNEW (.FALSE.): the lattice file has the old 'slim0' format
C          (.TRUE.) : the lattice file has the new 'slim1' format

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MLIN=132)
C>>>>>now the slim format is entering the `90s !!
      CHARACTER CLINE*(MLIN),CPOS*32,CASEUP*10
      PARAMETER(MAXINP=5000000)

C%%      write(*,*)'defslim'

      BSUPER = .FALSE.
      NELE = 0

      DO 100 I=1,MAXINP
         READ(15,'(A)',END=200) CLINE

         IF(CLINE(1:1).EQ.' ') THEN

            IF(BNEW .AND. ( CPOS(CLINE,1).EQ.'NSUPER' .OR.
     +              CPOS(CLINE,1).EQ.'NSUPER=' )) THEN
               NSUPER = INFCS(CPOS(CLINE(INDEX(CLINE,'=')+1:),1))
               BSUPER = .TRUE.
               GOTO 100
            ENDIF

            NELE = NELE + 1
            IF(NELE.GT.MAXELE) CALL TOFF(
     +           'DEFSLIM , too many elements, increase MAXELE')
            STRZ(NELE)  = 0.D+0 
            ISTR0(NELE) = 0
            ISTR1(NELE) = 0
            ITP         = INFCS(CPOS(CLINE,1))
            ITYPE(NELE) = ITP
            CNAME(NELE) = CASEUP(CPOS(CLINE,2))
            STRX(NELE)  = REFCS(CPOS(CLINE,3))
            STRY(NELE)  = REFCS(CPOS(CLINE,4))
            RLEN(NELE)  = REFCS(CPOS(CLINE,5))

            IF(ITP.EQ.IDL) THEN
C===========works for old and new !
               IF(RLEN(NELE).NE.0.D+0.AND.STRX(NELE).NE.0.D+0) THEN
                  CALL TOFF('DEFSLIM, slim error in drift')
               ELSE
                  RLEN(NELE) = STRX(NELE)+RLEN(NELE)
                  STRX(NELE) = 0.D+0
               ENDIF
            ELSEIF(ITP.EQ.IMR) THEN
               IF(RLEN(NELE).NE.0.D+0) CALL TOFF(
     +              'DEFSLIM, no length allowed for markers')
            ELSEIF(ITP.EQ.ISO) THEN
               STRX(NELE) = STRX(NELE)/RLEN(NELE)               
            ELSEIF(ITP.EQ.IRF) THEN
               STRX(NELE) = STRX(NELE)/1.D+3               
            ENDIF

            ISEC(NELE)  = INFCS(CPOS(CLINE,6))
            IF(ISEC(NELE).LT.0) CALL TOFF(
     +           'DEFSLIM, negative slicing')
            TWIS(NELE)  = REFCS(CPOS(CLINE,7))

            IF(BNEW) THEN
               IPS(NELE)   = INFCS(CPOS(CLINE,8))
               IFLAG(NELE) = IBFCS(CPOS(CLINE,9))
               IF(IPS(NELE).LT.0) THEN
                  IPS(NELE)   = -IPS(NELE)
                  IFLAG(NELE) = IBSET(IFLAG(NELE),JBSNEG)
               ENDIF
            ELSEIF(ITP.EQ.IMA) THEN
C==============twissmat also for old slim0 .......
               IF    (INFCS(CPOS(CLINE,8)).EQ.1) THEN
                  IFLAG(NELE) = IBSET(IFLAG(NELE),JBSPC7)
               ELSEIF(INFCS(CPOS(CLINE,8)).EQ.2) THEN
                  IFLAG(NELE) = IBSET(IFLAG(NELE),JBSPC8)
               ELSEIF(INFCS(CPOS(CLINE,8)).EQ.3) THEN
                  IFLAG(NELE) = IBSET(IFLAG(NELE),JBSPC9)
               ENDIF
            ELSE
               IPS(NELE)   = 0
               IFLAG(NELE) = 0
            ENDIF

            IF(ITP.EQ.ISR) THEN
               READ(15,'(A)',END=200) CLINE

Ctest> 
C      write(17,'(A/A132)')   'DEFSLIM:',CLINE
C      write(17,'(A,I6,2A12)')'        ',INFCS(CPOS(CLINE,1)),
C     +     CASEUP(CPOS(CLINE,2)),CNAME(NELE)
Ctest<

               IF(CLINE(1:1).EQ.' ' .AND. 
     +              INFCS(CPOS(CLINE,1)).EQ.0  .AND.   
     +              CASEUP(CPOS(CLINE,2)).EQ.CNAME(NELE)) THEN
                  STRZ(NELE)  = REFCS(CPOS(CLINE,3))
                  ISTR0(NELE) = INFCS(CPOS(CLINE,4))
                  ISTR1(NELE) = INFCS(CPOS(CLINE,5))
                  BSRM        = .TRUE.
               ELSE
                  CALL TOFF('DEFSLIM inconsistent continuation line')
               ENDIF
            ENDIF

            IF(IGLBSL.NE.0) THEN

               IF(ITP.EQ.ISO) THEN
                  IF(ISEC(NELE).GT.1 .AND. .NOT.BNEW) CALL TOFF(
     +                 'DEFSLIM, slicing not allowed for solenoids')
C=================in OLD format slicing might be fatal !!!
                  ISEC(NELE) = 0
C=================prevent global slicing for solenoids, spin rot., etc
               ELSEIF(ITP.EQ.ISS .OR. ITP.EQ.ISF) THEN
                  IF(ISEC(NELE).GT.1 .AND. .NOT.BNEW) CALL TOFF(
     +                 'DEFSLIM, slicing not allowed for spin rotat.')
                  ISEC(NELE) = 0
                ELSEIF(ITP.EQ.ISM) THEN
                  IF(ISEC(NELE).GT.1 .AND. .NOT.BNEW) CALL TOFF(
     +                 'DEFSLIM, slicing not allowed for spin match')
                   ISEC(NELE) = 0
               ELSEIF(ITP.EQ.IMR) THEN
                  IF(ISEC(NELE).GT.1 .AND. .NOT.BNEW) CALL TOFF(
     +                 'DEFSLIM, slicing not allowed for markers')
                   ISEC(NELE) = 0
               ELSEIF(ITP.EQ.IRF) THEN
                  IF(ISEC(NELE).GT.1 .AND. .NOT.BNEW) CALL TOFF(
     +                 'DEFSLIM, slicing not allowed for cavities')
C=================ramping!!!!!
                  ISEC(NELE) = 0
               ELSEIF(ITP.EQ.ISR) THEN
                  IF(ISEC(NELE).GT.1 .AND. .NOT.BNEW) CALL TOFF(
     +                 'DEFSLIM, slicing not allowed for srm-cell')
                  ISEC(NELE) = 0
               ELSEIF(ISEC(NELE)*IGLBSL.GT.1) THEN
                  D          = 1.0D+0/DBLE(ISEC(NELE)*IGLBSL)
                  RLEN(NELE) = RLEN(NELE) * D

                  IF(BNEW) THEN
                     STRX(NELE) = STRX(NELE) * D
                     STRY(NELE) = STRY(NELE) * D                     
                  ENDIF
               ENDIF

            ENDIF

            IF(CNAME(NELE).EQ.'END') THEN
C===========create new e-o-l
               CNAME(NELE) = 'X______EOL'
               ITYPE(NELE) = IMR
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = 0.D+0
               TWIS(NELE)  = 0.D+0
               ISEC(NELE)  = 0
               IFLAG(NELE) = 0
               IPS(NELE)   = 0

               IF(BNEW.AND. .NOT.BSUPER) THEN
                  NSUPER = 1
                  CALL WARN('DEFSLIM, no superperiod given -> 1')
               ENDIF

               RETURN
            ENDIF

         ELSE
            WRITE(*,'(A)') ' COMMENT:'
            WRITE(*,'(A)') CLINE
            WRITE(16,'(A)') CLINE(:ICSEND(CLINE))
         ENDIF

 100  CONTINUE

      CALL TOFF('DEFSLIM, file too long.'//
     +          ' increase MAXINP')
 200  CONTINUE
      CALL TOFF('DEFSLIM, no end')
C Hier sollten wir ein TOFF machen, weil n"amlich ohne END in der type-list
C die lattice-list nicht richtig interpretiert werden kann !  

      RETURN
      E N D

      S U B R O U T I N E LAYSLIM(BNEW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  LAYSLIM
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills the lattice of a particular machine with predefined
C     elements in 'SLIM' format

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MAXINP=5000000, MLIN=132)
      CHARACTER CLINE*(MLIN),CPOS*20,CDUM*20,CDUMPO*20,CDUMPR*14
      CHARACTER CSDUM*4

C%%      write(*,*)'layslim'

      IF(BNEW) THEN
         RLP2M = 1.D+0 
      ELSE
         RLP2M = 1.D-4
      ENDIF

      NLAT   = 0
      IEND   = NUMELE('X______EOL')
      ENDPOS = 0.D+0
      NCAV   = 0
      NTOTSL = 0

      DO 100 I=1,MAXINP
         READ(15,'(A)',END=200) CLINE

         IF(CLINE(1:1).NE.'-') THEN

            DO 110 IPOS=1,40,2
               CDUM = CPOS(CLINE,IPOS)
               IF(CDUM.EQ.' ') GOTO 100
               LELE = NUMELE(CDUM)
               IF(CDUM.EQ.'END') LELE = IEND
               
               IF(LELE.EQ.0) THEN
                  CALL WARN(
     +                 'LAYSLIM, unknown element: '
     +                 //CDUM//' -> predrift')
                  GOTO 110
               ELSEIF(LELE.EQ.ISO) THEN
                  IF(BNEW) THEN
                     RCORR = 0.D+0
                  ELSE
                     RCORR = RLEN(LELE)*.5D+0
                  ENDIF
               ENDIF
               
               NLAT = NLAT + 1
               IF(NLAT.GT.MAXLIST) CALL TOFF(
     +              'LAYSLIM, increase MAXLIST')
               LATELE(NLAT) = LELE

               IF(ITPEFF(LELE).EQ.ISR) THEN
                  LSLIC(NLAT) = -1
                  RLHLP = RLEN(LELE)               
               ELSEIF(ISEC(LELE)*IGLBSL.GT.1) THEN
                  LSLIC(NLAT) =  ISEC(LELE)*IGLBSL
                  RLHLP = RLEN(LELE) * DBLE(LSLIC(NLAT))
               ELSE
                  LSLIC(NLAT) = 1
                  RLHLP = RLEN(LELE)
               ENDIF
               
               NTOTSL = NTOTSL + MAX0(1,LSLIC(NLAT))
               
C===========works for both formats if we choose RLP2M, RCORR dependent
               PRETST =
     +              REFCS(CPOS(CLINE,IPOS+1))*RLP2M-ENDPOS-RLHLP*0.5D+0
     +              + RCORR
               
               IF(PRETST.LT.-1.D-10) THEN
                  CDUMPO = CPOS(CLINE,IPOS+1)
                  WRITE(CDUMPR,'(D14.6)') PRETST
                  CALL WARN('LAYSLIM, negative predrift'//
     +                 CDUMPR//' before '//
     +                 CDUM(:ICSEND(CDUM))//' at '//
     +                 CDUMPO(:ICSEND(CDUMPO)))
               ENDIF
               
               IF(DABS(PRETST).GT.1.D-9) THEN
                  PREDRF(NLAT) = PRETST
               ELSE
                  PREDRF(NLAT) = 0.D+0
               ENDIF
               
               IF(ITYPE(LELE).EQ.IRF) NCAV = NCAV + NSUPER
               
               IF(LELE.EQ.IEND) THEN 
                  
                  IF(NCAV.GT.1) THEN
                     WRITE(CSDUM,'(I4)') NCAV
                     CALL WARN('LAYSLIM, NCAV = '
     +                    //CSDUM//' might be time consuming at ramp')
                  ENDIF

                  DRINGL = DBLE(NSUPER) * (ENDPOS+PREDRF(NLAT)+RLHLP)

                  WRITE(16,'(2(A,I8/),A,F16.8)')
     +                 ' NUMBER OF DIFFERENT ELEMENTS : ',NELE,
     +                 ' TOT. NUMBER OF ELEM. IN LAT. : ',NLAT,
     +                 ' TOTAL LENGTH :         '        ,DRINGL

                  RETURN
               ENDIF
               
               ENDPOS = ENDPOS+PREDRF(NLAT)+RLHLP
 110        CONTINUE

         ELSE
            WRITE(*,'(A)') ' COMMENT:'
            WRITE(*,'(A)') CLINE
            WRITE(16,'(A)') CLINE
         ENDIF

 100  CONTINUE

      CALL TOFF('LAYSLIM, file too long.'//
     +          'increase MAXINP')
 200  CONTINUE
      CALL WARN('LAYSLIM, no ''END''')

      RETURN
      E N D


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C27                                                                           C
C@ PROGRAM=SPRINT, MODULE=PETROS, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      S U B R O U T I N E DEFPET
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  DEFPET
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads data from 15, loads '/ELEIN/' for creation
C     of 'TMAP'-array in sub. 'CREAMAP' in 'PETROS' format
C     The element definition part has to end with the first
C     occurence of the word 'END'.  Internally this creates a marker
C     named 'X______EOL'.
C
C     who volunteers to implement a half--way consistent PETROS implementation
C     of srm-cell ?????????????


      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MAXINP=5000000)
      CHARACTER CELNM*(MAELNM)
      CHARACTER CLINE*80,CPOS*20,CASEUP*10

C%%      write(*,*)'defpet'

      NELE = 0

C-----loop over all element definition lines of the 'PETROS' input file
      DO 100 I=1,MAXINP
         READ(15,'(A)',END=200) CLINE

C--------sort out lines which are commented out by a leading #
         IF(CLINE(1:1).NE.'#') THEN
            NELE = NELE + 1
            IF(NELE.GT.MAXELE) CALL TOFF(
     +           'DEFPET , too many elements.'//
     +           'increase maxele')
C==========='IPET' is a number which indicates the type of the element
            IPET        =  INFCS(CPOS(CLINE,2))
            CELNM       = CASEUP(CPOS(CLINE,1))
C==========='CNAME' is the name of the element
            CNAME(NELE) = CELNM

            READ(15,'(A)',END=200) CLINE

            PLEN =  REFCS(CPOS(CLINE,1))/1000.D+0
            PROT =  REFCS(CPOS(CLINE,2))
            PSTR = -REFCS(CPOS(CLINE,3))
            PXCU = -REFCS(CPOS(CLINE,4))
            PYCU = -REFCS(CPOS(CLINE,5))
            IFLAG(NELE) = 0
            IPS(NELE)   = 0

C======undefined
            IF(IPET.EQ.0) THEN
               IF(PLEN.NE.0.D+0)
     +              CALL WARN(
     +              'DEFPET, undefined in petros>drift: '//CELNM)
               ITYPE(NELE) = IDL
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======missing
            ELSEIF(IPET.EQ.1) THEN
               IF(PLEN.NE.0.D+0)
     +              CALL WARN(
     +              'DEFPET, missing in petros>drift: '//CELNM)
               ITYPE(NELE) = IDL
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======drit
            ELSEIF(IPET.EQ.2) THEN
               ITYPE(NELE) = IDL
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======cavity
            ELSEIF(IPET.EQ.3) THEN
               ITYPE(NELE) =  IRF
               STRX(NELE)  = -PSTR
               STRY(NELE)  =  1.D+0
               RLEN(NELE)  =  PLEN
               TWIS(NELE)  = 0.D+0
C======interaction point / marker
            ELSEIF(IPET.EQ.4) THEN
               ITYPE(NELE) = IMR
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======monitor
            ELSEIF(IPET.EQ.5) THEN
               ITYPE(NELE) = IMO
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======rbend
            ELSEIF(IPET.EQ.6) THEN
               IF(PROT.EQ.0.D+0) THEN
                  IF(PYCU.NE.0.D+0) CALL TOFF(
     +                 'DEFPET, petros error in rbend: '//CELNM)
                  ITYPE(NELE) = IRX
                  IF(PXCU.NE.PSTR) THEN
                     IF(PXCU*PSTR.NE.0.D+0) THEN
                        CALL TOFF(
     +                       'DEFPET, petros error in rbend: '//CELNM)
                     ELSE
                        CALL WARN(
     +                       'DEFPET, PXCU.ne.PSTR for rbend: '//CELNM)
                        PXCU = PXCU + PSTR
                        PSTR = PXCU
                     ENDIF
                  ENDIF
                  STRX(NELE)  = PLEN*PXCU
               ELSEIF(DABS(DABS(PROT)-.15707963E+01).LT.1E-6) THEN
                  ITYPE(NELE) = IRY
                  IF(PYCU.NE.PSTR) THEN
                     IF(PYCU*PSTR.NE.0.D+0) THEN
                        CALL TOFF(
     +                       'DEFPET, petros error in rbend: '//CELNM)
                     ELSE
                        CALL WARN(
     +                       'DEFPET, PYCU.ne.PSTR for rbend: '//CELNM)
                        PYCU = PYCU + PSTR
                        PSTR = PYCU
                     ENDIF
                  ENDIF
                  IF(DABS(PROT+.15707963E+01).LT.1E-6) THEN
                     STRX(NELE)  = -PLEN*PYCU
                  ELSEIF(DABS(PROT-.15707963E+01).LT.1E-6) THEN
                     STRX(NELE)  =  PLEN*PYCU
                  ELSE
                     CALL TOFF('DEFPET, internal error')
                  ENDIF
               ELSE
                  CALL TOFF(
     +                 'DEFPET, no rotated rbends yet: '//CELNM)
               ENDIF
               STRY(NELE) = 0.D+0
               RLEN(NELE) = PLEN
               TWIS(NELE) = 0.D+0
C======electrostatic plate
            ELSEIF(IPET.EQ.7) THEN
               CALL WARN(
     +              'DEFPET, electr. plate not defed yet>drift '
     +              //CELNM)
               ITYPE(NELE) = IDL
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======corrector dipole
            ELSEIF(IPET.EQ.8) THEN
               IF(PROT.EQ.0.D+0) THEN
                  IF(PYCU.NE.0.D+0) CALL TOFF(
     +            'DEFPET, petros error in corrector dipole: '//CELNM)
                  ITYPE(NELE) = IKX
                  IF(PXCU.NE.PSTR) THEN
                     IF(PXCU*PSTR.NE.0.D+0) THEN
                        CALL TOFF(
     +            'DEFPET, petros error in corrector dipole: '//CELNM)
                     ELSE
                        CALL WARN(
     +            'DEFPET, PXCU.ne.PSTR for corrector dipole: '//CELNM)
                        PXCU = PXCU + PSTR
                        PSTR = PXCU
                     ENDIF
                  ENDIF
                  STRX(NELE)  = PLEN*PXCU
               ELSEIF(DABS(DABS(PROT)-.15707963E+01).LT.1E-6) THEN
                  ITYPE(NELE) = IKY
                  IF(PYCU.NE.PSTR) THEN
                     IF(PYCU*PSTR.NE.0.D+0) THEN
                        CALL TOFF(
     +             'DEFPET, petros error in corrector dipole: '//CELNM)
                     ELSE
                        CALL WARN(
     +             'DEFPET, PYCU.ne.PSTR for corrector dipole: '//CELNM)
                        PYCU = PYCU + PSTR
                        PSTR = PYCU
                     ENDIF
                  ENDIF
                  IF(DABS(PROT+.15707963E+01).LT.1E-6) THEN
                     STRX(NELE)  = -PLEN*PYCU
                  ELSEIF(DABS(PROT-.15707963E+01).LT.1E-6) THEN
                     STRX(NELE)  =  PLEN*PYCU
                  ELSE
                     CALL TOFF('DEFPET, internal error')
                  ENDIF
               ELSE
                  CALL TOFF(
     +            'DEFPET, no rotated corrector dipole yet: '//CELNM)
               ENDIF
               STRY(NELE) = 0.D+0
               RLEN(NELE) = PLEN
               TWIS(NELE) = 0.D+0
C======quadrupole
            ELSEIF(IPET.EQ.9) THEN
               IF(DABS(PXCU)+DABS(PYCU).EQ.0.D+0) THEN
                  IF(PROT.EQ.0.D+0) THEN
                     ITYPE(NELE) = IMQ
                     TWIS(NELE)  = 0.D+0
                     STRX(NELE)  = PSTR*PLEN
                  ELSEIF(DABS(PROT-PI/4.D+0).LT.1.D-6) THEN
                     ITYPE(NELE) = ISQ
                     TWIS(NELE)  = 0.D+0
                     STRX(NELE)  = PSTR*PLEN
                  ELSEIF(DABS(PROT+PI/4.D+0).LT.1.D-6) THEN
                     ITYPE(NELE) = ISQ
                     TWIS(NELE)  = 0.D+0
                     STRX(NELE)  = -PSTR*PLEN
                  ELSE
                     ITYPE(NELE) = IRQ
                     TWIS(NELE)  = PROT
                     STRX(NELE)  = PSTR*PLEN
                  ENDIF

                  STRY(NELE) = 0.D+0
                  RLEN(NELE) = PLEN
               ELSE
                  CALL WARN(
     +                 'DEFPET, quadrupole has curvature>rect cbf: '
     +                 //CELNM)
                  IF(PROT.EQ.0.D+0) THEN
                     ITYPE(NELE) = ICH
                     STRX(NELE)  = PLEN*PXCU
                     STRY(NELE)  = PSTR*PLEN
                     IF(PYCU.NE.0.D+0) CALL TOFF(
     +                    'DEFPET, petros error in rect cbf: '//CELNM)
                  ELSEIF(DABS(DABS(PROT)-.15707963E+01).LT.1E-6) THEN
                     ITYPE(NELE) = ICV
                     IF(DABS(PROT+.15707963E+01).LT.1E-6) THEN
                        STRX(NELE)  = -PLEN*PYCU
                     ELSEIF(DABS(PROT-.15707963E+01).LT.1E-6) THEN
                        STRX(NELE)  =  PLEN*PYCU
                     ELSE
                        CALL TOFF('DEFPET, internal error')
                     ENDIF
                     STRY(NELE)  = PSTR*PLEN
                  ELSE
                     CALL TOFF(
     +               'DEFPET, no rotated cbfs yet: '//CELNM)
                  ENDIF
                  RLEN(NELE) = PLEN
                  TWIS(NELE) = 0.D+0
               ENDIF
C======electric quadrupole
            ELSEIF(IPET.EQ.10) THEN
               CALL WARN(
     +              'DEFPET, electric quad not defed yet>drift '
     +              //CELNM)
               ITYPE(NELE) = IDL
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======sextupole
            ELSEIF(IPET.EQ.11) THEN
               ITYPE(NELE) = IMS
               STRX(NELE)  = PSTR*PLEN
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======octupole
            ELSEIF(IPET.EQ.12) THEN
               CALL WARN(
     +              'DEFPET, octupole not defed yet>drift '//CELNM)
               ITYPE(NELE) = IDL
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======solenoid
            ELSEIF(IPET.EQ.13) THEN
               ITYPE(NELE) = ISO
               STRX(NELE)  = PSTR
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======thin lens quad
            ELSEIF(IPET.EQ.14) THEN
               IF(BTHICK) CALL WARN(
     +              'DEFPET, all elements are thick: '//CELNM)
               IF(DABS(PXCU)+DABS(PYCU).EQ.0.D+0) THEN
                  IF(PROT.EQ.0.D+0) THEN
                     ITYPE(NELE) = IMQ
                     TWIS(NELE)  = 0.D+0
                  ELSE
                     ITYPE(NELE) = IRQ
                     TWIS(NELE)  = PROT
                  ENDIF
                  STRX(NELE) = PSTR*PLEN
                  STRY(NELE) = 0.D+0
                  RLEN(NELE) = PLEN
               ELSE
                  CALL WARN(
     +                 'DEFPET, thin quad has curvature>rect cbf: '
     +                 //CELNM)
                  IF(PROT.EQ.0.D+0) THEN
                     ITYPE(NELE) = ICH
                     STRX(NELE)  = PLEN*PXCU
                     STRY(NELE)  = PSTR*PLEN
                     IF(PYCU.NE.0.D+0)
     +                    CALL TOFF(
     +                    'DEFPET, petros error in thin quad: '//CELNM)
                  ELSEIF(DABS(DABS(PROT)-.15707963E+01).LT.1E-6) THEN
                     ITYPE(NELE) = ICV
                     IF(DABS(PROT+.15707963E+01).LT.1E-6) THEN
                        STRX(NELE)  = -PLEN*PYCU
                     ELSEIF(DABS(PROT-.15707963E+01).LT.1E-6) THEN
                        STRX(NELE)  =  PLEN*PYCU
                     ELSE
                        CALL TOFF('DEFPET, internal error')
                     ENDIF
                     STRY(NELE)  = PSTR*PLEN
                  ELSE
                     CALL TOFF('DEFPET, no rotated bends yet: '
     +                    //CELNM)
                  ENDIF
                  RLEN(NELE) = PLEN
                  TWIS(NELE) = 0.D+0
               ENDIF
C======x direction monitor
            ELSEIF(IPET.EQ.17) THEN
               CALL WARN(
     +              'DEFPET, no specific monitor -> general '//CELNM)
               ITYPE(NELE) = IMO
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======y direction monitor
            ELSEIF(IPET.EQ.18) THEN
               CALL WARN(
     +              'DEFPET, no specific monitor -> general '//CELNM)
               ITYPE(NELE) = IMO
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======kicker
            ELSEIF(IPET.EQ.19) THEN
               CALL WARN(
     +              'DEFPET, kicks are treated as correction dipoles')
               IF(PROT.EQ.0.D+0) THEN
                  IF(PYCU.NE.0.D+0) CALL TOFF(
     +                 'DEFPET, petros error in kicker : '//CELNM)
                  ITYPE(NELE) = IKX
                  IF(PXCU.NE.PSTR) THEN
                     IF(PXCU*PSTR.NE.0.D+0) THEN
                        CALL TOFF(
     +                  'DEFPET, petros error in kicker: '//CELNM)
                     ELSE
                        CALL WARN(
     +                  'DEFPET, PXCU.ne.PSTR for kicker: '//CELNM)
                        PXCU = PXCU + PSTR
                        PSTR = PXCU
                     ENDIF
                  ENDIF
                  STRX(NELE)  = PLEN*PXCU
               ELSEIF(DABS(DABS(PROT)-.15707963E+01).LT.1E-6) THEN
                  ITYPE(NELE) = IKY
                  IF(PYCU.NE.PSTR) THEN
                     IF(PYCU*PSTR.NE.0.D+0) THEN
                        CALL TOFF(
     +                  'DEFPET, petros error in kicker: '//CELNM)
                     ELSE
                        CALL WARN(
     +                  'DEFPET, PYCU.ne.PSTR for kicker: '//CELNM)
                        PYCU = PYCU + PSTR
                        PSTR = PYCU
                     ENDIF
                  ENDIF
                  IF(DABS(PROT+.15707963E+01).LT.1E-6) THEN
                     STRX(NELE)  = -PLEN*PYCU
                  ELSEIF(DABS(PROT-.15707963E+01).LT.1E-6) THEN
                     STRX(NELE)  =  PLEN*PYCU
                  ELSE
                     CALL TOFF('DEFPET, internal error')
                  ENDIF
               ELSE
                  CALL TOFF(
     +                 'DEFPET, no rotated kicker yet: '//CELNM)
               ENDIF
               STRY(NELE) = 0.D+0
               RLEN(NELE) = PLEN
               TWIS(NELE) = 0.D+0
C======sbend
            ELSEIF(IPET.EQ.21) THEN
               IF(PROT.EQ.0.D+0) THEN
                  IF(PYCU.NE.0.D+0) CALL TOFF(
     +                 'DEFPET, petros error in sbend : '//CELNM)
                  ITYPE(NELE) = IBX
                  IF(PXCU.NE.PSTR) THEN
                     IF(PXCU*PSTR.NE.0.D+0) THEN
                        CALL TOFF(
     +                       'DEFPET, petros error in sbend: '//CELNM)
                     ELSE
                        CALL WARN(
     +                       'DEFPET, PXCU.ne.PSTR for sbend: '//CELNM)
                        PXCU = PXCU + PSTR
                        PSTR = PXCU
                     ENDIF
                  ENDIF
                  STRX(NELE)  = PLEN*PXCU
               ELSEIF(DABS(DABS(PROT)-.15707963E+01).LT.1E-6) THEN
                  ITYPE(NELE) = IBY
                  IF(PYCU.NE.PSTR) THEN
                     IF(PYCU*PSTR.NE.0.D+0) THEN
                        CALL TOFF(
     +                       'DEFPET, petros error in sbend: '//CELNM)
                     ELSE
                        CALL WARN(
     +                       'DEFPET, PYCU.ne.PSTR for sbend: '//CELNM)
                        PYCU = PYCU + PSTR
                        PSTR = PYCU
                     ENDIF
                  ENDIF
                  IF(DABS(PROT+.15707963E+01).LT.1E-6) THEN
                     STRX(NELE)  = -PLEN*PYCU
                  ELSEIF(DABS(PROT-.15707963E+01).LT.1E-6) THEN
                     STRX(NELE)  =  PLEN*PYCU
                  ELSE
                     CALL TOFF('DEFPET, internal error')
                  ENDIF
               ELSE
                  CALL TOFF('DEFPET, no rotated sbends yet: '
     +                 //CELNM)
               ENDIF
               STRY(NELE) = 0.D+0
               RLEN(NELE) = PLEN
               TWIS(NELE) = 0.D+0
C======combined function sbend
            ELSEIF(IPET.EQ.22) THEN
               IF(PROT.EQ.0.D+0) THEN
                  IF(PYCU.NE.0.D+0) CALL TOFF(
     +                 'DEFPET, petros error in sect cbf : '//CELNM)
                  ITYPE(NELE) = ICX
                  STRX(NELE)  = PLEN*PXCU
                  STRY(NELE)  = PLEN*PSTR
               ELSEIF(DABS(DABS(PROT)-.15707963E+01).LT.1E-6) THEN
                  ITYPE(NELE) = ICY
                  IF(DABS(PROT+.15707963E+01).LT.1E-6) THEN
                     STRX(NELE)  = -PLEN*PYCU
                  ELSEIF(DABS(PROT-.15707963E+01).LT.1E-6) THEN
                     STRX(NELE)  =  PLEN*PYCU
                  ELSE
                     CALL TOFF('DEFPET, internal error')
                  ENDIF
                  STRY(NELE)  = PLEN*PSTR
               ELSE
                  CALL TOFF('DEFPET, no rotated sect cbfs yet: '
     +                 //CELNM)
               ENDIF
               RLEN(NELE) = PLEN
               TWIS(NELE) = 0.D+0
C======snake (like 'SLIM'):
C======0, hor. azimuth, angle from hor., rot. angle, all in degree
            ELSEIF(IPET.EQ.23) THEN
               ITYPE(NELE) = ISS
               STRX(NELE)  = PSTR
               STRY(NELE)  = PXCU
               IF(PROT.NE.0.D+0) THEN
                  CALL WARN(
     +                 'DEFPET, snakes do not need rotations: '//CELNM)
                  STRY(NELE) = STRY(NELE) + PROT/DEG2RD
               ENDIF
               RLEN(NELE) = PLEN
               TWIS(NELE) = PYCU
C======free snake
            ELSEIF(IPET.EQ.24) THEN
               ITYPE(NELE) = ISF
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = PLEN
               TWIS(NELE)  = 0.D+0
C======instantanious spin match
            ELSEIF(IPET.EQ.25) THEN
               ITYPE(NELE) = ISM
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               IF(PLEN.NE.0.D+0)
     +              CALL TOFF(
     +              'DEFPET, length for instantanious spin match')
               RLEN(NELE)  = 0.D+0
               TWIS(NELE)  = 0.D+0
            ELSE
               CALL TOFF('DEFPET, undefined element: '//CELNM)
            ENDIF

            ISEC(NELE) = 1

            IF(IGLBSL.GT.0) THEN
               ITP = ITYPE(NELE)
               IF(ITP.EQ.ISO .OR. ITP.EQ.ISS .OR. ITP.EQ.ISF .OR.
     +              ITP.EQ.ISM .OR. ITP.EQ.IMR .OR. ITP.EQ.IRF) THEN
                  ISEC(NELE) = 0
C=================prevent global slicing for solenoids, spin rot., etc
               ELSE
                  D          = 1.0D+0/DBLE(IGLBSL)
                  STRX(NELE) = STRX(NELE) * D
                  STRY(NELE) = STRY(NELE) * D
                  RLEN(NELE) = RLEN(NELE) * D
               ENDIF

            ENDIF

            IFLAG(NELE) = 0

            IF(CNAME(NELE).EQ.'END') THEN
C===========create new e-o-l
               CNAME(NELE) = 'X______EOL'
               ITYPE(NELE) = IMR
               STRX(NELE)  = 0.D+0
               STRY(NELE)  = 0.D+0
               RLEN(NELE)  = 0.D+0
               TWIS(NELE)  = 0.D+0
               ISEC(NELE)  = 0
               IFLAG(NELE) = 0
               IPS(NELE)   = 0

               RETURN
            ENDIF


         ELSE
            WRITE(*,'(A)') ' COMMENT:'
            WRITE(*,'(A)') CLINE
            WRITE(16,'(A)') CLINE
         ENDIF

 100  CONTINUE

      CALL TOFF('DEFPET, file too long.'//
     +          ' increase MAXINP')
 200  CONTINUE
      CALL TOFF('DEFPET, no end')
C Hier sollten wir ein TOFF machen, weil n"amlich ohne END in der type-list
C der ANFANG der lattice-list nicht gefunden werden kann !  

      RETURN
      E N D

      S U B R O U T I N E LAYPET
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  LAYPET
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills the lattice of a particular machine with predefined
C     elements in 'PETROS' format.
C     The lattice ends with the first occurence of an element named 'END'.

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MAXINP=5000000, RLP2M = 1.D-4)
      CHARACTER CLINE*80,CPOS*30,CDUM*30, CSDUM*4

C%%      write(*,*)'laypet'

      NLAT   = 0
      IEND   = NUMELE('X______EOL')
      ENDPOS = 0.D+0
      NCAV   = 0
      NTOTSL = 0
      BNOEND = .TRUE.

C-----read all lines of the 'PETROS' file after the element definition part
      DO 100 I=1,MAXINP
         READ(15,'(A)',END=101) CLINE

C--------sort out lines which are commented out by a leading #
         IF(CLINE(1:1).NE.'#') THEN
            NLAT = NLAT + 1
            IF(NLAT.GT.MAXLIST) CALL TOFF(
     +           'LAYPET, increase MAXLIST')

            LELE = INFCS(CPOS(CLINE,1))
            IF(LELE.EQ.0.OR.LELE.GT.NELE) THEN
               WRITE(CDUM,'(I10,I10)') LELE,NLAT
               CALL TOFF('LAYPET, unknown element #'//CDUM(1:10)//
     +              ' at lattice entry #'//CDUM(11:20))
            ENDIF
            LATELE(NLAT) = LELE

            IF(ISEC(LELE)*IGLBSL.GT.1) THEN
               LSLIC(NLAT) =  ISEC(LELE)*IGLBSL
            ELSE
               LSLIC(NLAT) = 1
            ENDIF

            NTOTSL = NTOTSL + MAX0(1,LSLIC(NLAT))

            PRETST = REFCS(CPOS(CLINE,2))/1000.D+0
            IF( DABS(REFCS(CPOS(CLINE,3))/1000.D+0 -
     +              RLEN(LELE)*DMAX1(1.D+0,DBLE(LSLIC(NLAT))))
     +              .GT.1.D-6 ) THEN
               WRITE(CDUM,'(I10,I10)'), LELE, NLAT
               CALL TOFF('LAYPET, length inconsistency (A) in petros '//
     +              'at element #'//CDUM(1:10)//' at lattice entry'//
     +              CDUM(11:20))
            ENDIF
            DUM = REFCS(CPOS(CLINE,4))/1000.D+0

            IF(NLAT.GT.1) THEN
               IF(DABS(ENDPOS+PREDRF(NLAT-1)+RLEN(LATELE(NLAT-1))-DUM)
     +              .GT.1.D-2) THEN 
                  WRITE(CDUM,'(I10,I10)'), LELE, NLAT
                  CALL TOFF(
     +                 'LAYPET, length inconsistency (B) in petros '//
     +                 'at element #'//CDUM(1:10)//' at lattice entry'//
     +                 CDUM(11:20))
               ENDIF
            ENDIF
            
            IF(DABS(PRETST).GT.1.D-9) THEN
               PREDRF(NLAT) = PRETST
            ELSE
               PREDRF(NLAT) = 0.D+0
            ENDIF

            ENDPOS = DUM
            IF(ITYPE(LELE).EQ.IRF) NCAV = NCAV + NSUPER
            READ(15,'(A)',END=101) CLINE

            IF(LELE.EQ.IEND) THEN
               BNOEND = .FALSE.
               GOTO 101
            ENDIF
         ELSE            
            WRITE(*,'(A)') ' COMMENT:'
            WRITE(*,'(A)') CLINE
            WRITE(16,'(A)') CLINE
         ENDIF
 100  CONTINUE

      CALL TOFF('LAYPET, file too long.'//
     +          'increase MAXINP')

 101  CONTINUE

      IF(BNOEND) THEN
         NLAT = NLAT + 1
         IF(NLAT.GT.MAXLIST) CALL TOFF(
     +        'LAYPET, increase MAXLIST')
         LATELE(NLAT) = IEND
         PREDRF(NLAT) = 0.D+0
         LSLIC (NLAT) = 1
         NTOTSL       = NTOTSL + 1     
      ENDIF

      IF(NCAV.GT.1) THEN
         WRITE(CSDUM,'(I4)') NCAV
         CALL WARN('LAYPET, NCAV = '
     +        //CSDUM//' might be time consuming at ramp')
      ENDIF

      DRINGL = DBLE(NSUPER) * ENDPOS
      
      WRITE(16,'(2(A,I8/),A,F16.8)')
     +     ' NUMBER OF DIFFERENT ELEMENTS : ',NELE,
     +     ' TOT. NUMBER OF ELEM. IN LAT. : ',NLAT,
     +     ' TOTAL LENGTH :         '        ,DRINGL

      RETURN
      E N D

      S U B R O U T I N E PETROSBE(IUPETROS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  PETROSBE
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Writes four arbitrary lines which are by default skipped in a PETROS
C     file.
C
      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C%%      write(*,*)'petrosbe'

      WRITE(IUPETROS,'(3(A/),A)')
     +     'PETROS skips four lines at the beginning:',
     +     '2nd line','3rd line','4th line'

      RETURN
      E N D

      S U B R O U T I N E PETROSELEOUT(IUPETROS,BPR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  PETROSELEOUT
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      CHARACTER CNCHECK*(MAELNM),CN*(MAELNM)

C%%      write(*,*)'petroseleout'

      DO 100 IELE=1,NELE

         ITP = ITYPE(IELE)

C========undo slicing
         DSL = 1.D+0
         IF(ITP.NE.ISO .AND. ITP.NE.IRF .AND. ITP.NE.ISS .AND.
     +      ITP.NE.ISF .AND. ITP.NE.ISM .AND. ITP.NE.IMR) THEN
            J = ISEC(IELE)*IGLBSL
            IF(J.GT.1) DSL = DBLE(J)
         ENDIF

         CN = CNCHECK(CNAME(IELE),'ALL')
         CALL PETROSELE(IUPETROS,IELE,CN,DSL*RLEN(IELE),
     +        DSL*STRX(IELE),DSL*STRY(IELE),BPR)
 100  CONTINUE

      RETURN
      E N D

      S U B R O U T I N E PETROSELE(IUPETROS,IELE,
     +     CNAMEI,RLENI,STRXI,STRYI,BPR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  PETROSELE
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      CHARACTER CNAMEI*(*),CLINE*(40)

C%%      write(*,*)'petrosele'

      ITP      = ITYPE (IELE)
      ITP2     = ITPEFF(IELE)
      ISECI    = ISEC  (IELE)
      TWISI    = TWIS  (IELE)
      IFLAGI   = IFLAG (IELE)
      IPSI     = IPS   (IELE)

      PLEN     = 1000.D+0*RLENI

C     Warning: IFLAG and IPS are not includend in the PETROS file.
C     Including them wouldn't be worth the effort since *.pet files suck !

*******************************************
*     Drift                               *
*******************************************
      IF    (ITP.EQ.IDL) THEN
         IF(IELE.EQ.1.AND.PLEN.EQ.0.D+0) THEN
            CALL WARN('PETROSELE, first element gets type 1')
            ITY=1
         ELSE
            ITY=2
         ENDIF
         WRITE(CLINE,'(I18)') ITY
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = 0.D+0
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Monitor                             *
*******************************************
      ELSEIF(ITP.EQ.IMO) THEN
         WRITE(CLINE,'(I18)') 5
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = 0.D+0
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Horizontal Wedge Dipole  (SBEND)    *
*******************************************
      ELSEIF(ITP.EQ.IBX) THEN
         WRITE(CLINE,'(I18)') 21
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = DIVI0(STRXI,RLENI)
         IF(BPR) THEN
            PXCU = PSTR
         ELSE
            PXCU = 0.D+0
         ENDIF
         PYCU = 0.D+0
*******************************************
*     Vertical  Wedge Dipole   (SBEND)    *
*******************************************
      ELSEIF(ITP.EQ.IBY) THEN
         WRITE(CLINE,'(I18)') 21
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = PI/2.D+0
         PSTR = DIVI0(STRXI,RLENI)
         PXCU = 0.D+0
         IF(BPR) THEN
            PYCU = PSTR
         ELSE
            PYCU = 0.D+0
         ENDIF
*******************************************
*     Quadrupole                          *
*******************************************
      ELSEIF(ITP.EQ.IMQ) THEN
         WRITE(CLINE,'(I18)') 9
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = DIVI0(STRXI,RLENI)
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Skew Quadrupole                     *
*******************************************
      ELSEIF(ITP.EQ.ISQ) THEN
         WRITE(CLINE,'(I18)') 9
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
C         PROT = .657079644
         PROT = PI/4.D+0
C the sign of prot needs to be tested
         PSTR = DIVI0(STRXI,RLENI)
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Cavity                              *
*******************************************
      ELSEIF(ITP.EQ.IRF) THEN
         WRITE(CLINE,'(I18)') 3
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = -STRXI*STRYI
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Horizontal Kicker                   *
*******************************************
      ELSEIF(ITP.EQ.IKX) THEN
         WRITE(CLINE,'(I18)') 8
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = DIVI0(STRXI,RLENI)
         IF(BPR) THEN
            PXCU = PSTR
         ELSE
            PXCU = 0.D+0
         ENDIF
         PYCU = 0.D+0
*******************************************
*     Vertical Kicker                     *
*******************************************
      ELSEIF(ITP.EQ.IKY) THEN
         WRITE(CLINE,'(I18)') 8
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = PI/2.D+0
         PSTR = DIVI0(STRXI,RLENI)
         PXCU = 0.D+0
         IF(BPR) THEN
            PYCU = PSTR
         ELSE
            PYCU = 0.D+0
         ENDIF
*******************************************
*     Sextupole                           *
*******************************************
      ELSEIF(ITP.EQ.IMS) THEN
         WRITE(CLINE,'(I18)') 11
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = DIVI0(STRXI,RLENI)
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Solenoid                            *
*******************************************
      ELSEIF(ITP.EQ.ISO) THEN
         WRITE(CLINE,'(I18)') 13
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = STRXI
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Rotated Quadrupole                  *
*******************************************
      ELSEIF(ITP.EQ.IRQ) THEN
         WRITE(CLINE,'(I18)') 9
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = TWISI
         PSTR = DIVI0(STRXI,RLENI)
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Horizontal rectangular bend         *
*******************************************
      ELSEIF(ITP.EQ.IRX) THEN
         WRITE(CLINE,'(I18)') 6
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = DIVI0(STRXI,RLENI)
         IF(BPR) THEN
            PXCU = PSTR
         ELSE
            PXCU = 0.D+0
         ENDIF
         PYCU = 0.D+0
*******************************************
*     Vertical rectangular bend           *
*******************************************
      ELSEIF(ITP.EQ.IRY) THEN
         WRITE(CLINE,'(I18)') 6
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = PI/2.D+0
         PSTR = DIVI0(STRXI,RLENI)
         PXCU = 0.D+0
         IF(BPR) THEN
            PYCU = PSTR
         ELSE
            PYCU = 0.D+0
         ENDIF
*******************************************
*     Instantanious spin match            *
*******************************************
      ELSEIF(ITP.EQ.ISM) THEN
         CALL TOFF('PETROSELE, ISM not implemented')
*******************************************
*     Reads Map From File                 *
*******************************************
      ELSEIF(ITP.EQ.IRO) THEN
         CALL TOFF('PETROSELE, type 14 needs a stored map')
*******************************************************
*     Horizontal Combined Function Rectangular Magnet *
*******************************************************
      ELSEIF(ITP.EQ.ICH) THEN
         WRITE(CLINE,'(I18)') 9
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PYCU = 0.D+0

         IF(ITP2.EQ.IMQ) THEN
            PSTR = DIVI0(STRXI,RLENI)
            PXCU = 0.D+0
         ELSE
            PSTR = DIVI0(STRYI,RLENI)
            PXCU = DIVI0(STRXI,RLENI)
         ENDIF

**************************************************
*     Horizontal Combined Function Sector Magnet *
**************************************************
      ELSEIF(ITP.EQ.ICX) THEN
         WRITE(CLINE,'(I18)') 22
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PYCU = 0.D+0

         IF(ITP2.EQ.IMQ) THEN
            PSTR = DIVI0(STRXI,RLENI)
            PXCU = 0.D+0
         ELSE
            PSTR = DIVI0(STRYI,RLENI)
            PXCU = DIVI0(STRXI,RLENI)
         ENDIF

*****************************************************
*     Vertical Combined Function Rectangular Magnet *
*****************************************************
      ELSEIF(ITP.EQ.ICV) THEN
         WRITE(CLINE,'(I18)') 9
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = PI/2.D+0
         PXCU = 0.D+0

         IF(ITP2.EQ.IMQ) THEN
            PSTR = DIVI0(STRXI,RLENI)
            PYCU = 0.D+0
         ELSE
            PSTR = DIVI0(STRYI,RLENI)
            PYCU = DIVI0(STRXI,RLENI)
         ENDIF

************************************************
*     Vertical Combined Function Sector Magnet *
************************************************
      ELSEIF(ITP.EQ.ICY) THEN
         WRITE(CLINE,'(I18)') 22
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = PI/2.D+0
         PXCU = 0.D+0

         IF(ITP2.EQ.IMQ) THEN
            PSTR = DIVI0(STRXI,RLENI)
            PYCU = 0.D+0
         ELSE
            PSTR = DIVI0(STRYI,RLENI)
            PYCU = DIVI0(STRXI,RLENI)
         ENDIF

*******************************************
*     Snake                               *
*******************************************
      ELSEIF(ITP.EQ.ISS) THEN
         WRITE(CLINE,'(I18)') 23
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = STRXI
         PXCU = STRYI
         PYCU = TWISI
*******************************************
*     Free Snake                          *
*******************************************
      ELSEIF(ITP.EQ.ISF) THEN
         WRITE(CLINE,'(I18)') 24
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = 0.D+0
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Instantanious spin match            *
*******************************************
      ELSEIF(ITP.EQ.ISM) THEN
         WRITE(CLINE,'(I18)') 25
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = 0.D+0
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Marker                              *
*******************************************
      ELSEIF(ITP.EQ.IMR) THEN
         IF(CNAME(IELE).EQ.'X______EOL') THEN
            CALL WARN('PETROSELE, last element gets type 0')
            ITY=0
         ELSE
            ITY=4
         ENDIF
         WRITE(CLINE,'(I18)') ITY
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = 0.D+0
         PXCU = 0.D+0
         PYCU = 0.D+0
*******************************************
*     Unknown type of element             *
*******************************************
      ELSE
         WRITE(IUPETROS,'(9X,A,I6,A,A,A,E9.2,A)')
     +        '# unknown element type ',ITP,' - ',
     +        CNAMEI(:ICSEND(CNAMEI)),' with strenght ',STRXI,
     +        'make marker'
         WRITE(CLINE,'(I18)') 4
         CLINE(:ICSEND(CNAMEI)) = CNAMEI(:ICSEND(CNAMEI))
         WRITE(IUPETROS,'(A,10X,A,I5)') CLINE,'#',IELE
         PROT = 0.D+0
         PSTR = 0.D+0
         PXCU = 0.D+0
         PYCU = 0.D+0
      ENDIF

      WRITE(IUPETROS,'(5E16.8)') PLEN,PROT,-PSTR,-PXCU,-PYCU

      RETURN
      E N D

      S U B R O U T I N E PETROSLATOUT(IUPETROS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  PETROSLATOUT
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      CHARACTER CNCHECK*(MAELNM),CN*(MAELNM)

C%%      write(*,*)'petroslatout'

      ENDPOS = 0.D+0
      IEND   = NUMELE('X______EOL')
      PREEND = 0.D+0

      DO 100 ISUP=1,NSUPER
         BSUP = ISUP.EQ.NSUPER
C         WRITE(IUPETROS,'(A,I5,X,A,I5)')
C     +        '# Beginning of super-period',ISUP,'of',NSUPER

         DO 110 ILAT=1,NLAT
            IELE=LATELE(ILAT)

            IF( IELE.NE.IEND .OR. (BSUP.AND.(ILAT.EQ.NLAT)) ) THEN !screwpetros
C            IF( IELE.NE.IEND ) THEN ! petros dies with END marker in lattice
               CN = CNCHECK(CNAME(IELE),'ALL')
               WRITE(IUPETROS,'(I10,3E16.8,7X,A/,4E16.8)')
     +              IELE, (PREDRF(ILAT)+PREEND)*1000.D+0,
     +              RLEN(IELE)*1000.D+0, ENDPOS*1000.D+0,
     +              CN(:ICSEND(CN)),
     +              0.D+0,0.D+0,0.D+0,0.D+0
               ENDPOS = ENDPOS + PREDRF(ILAT) + PREEND + RLEN(IELE)
               PREEND = 0.D+0
            ELSE
               IF(RLEN(IELE).NE.0)
     +              CALL TOFF('PETROSLATOUT, END with length')
               PREEND = PREEND + PREDRF(ILAT)
            ENDIF

 110     CONTINUE
 100  CONTINUE

      RETURN
      E N D

      S U B R O U T I N E PETROSEN(IUPETROS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  PETROSEN
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     This subroutine creates the end and execution part of the
C     PETROS program.

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C%%      write(*,*)'petrosen'

      WRITE(IUPETROS,'(A))') '# WRPETROS finished'

      RETURN
      E N D


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C28                                                                           C
C@ PROGRAM=SPRINT, MODULE=MAD0, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      S U B R O U T I N E RDMAD0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  RDMAD0
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     driver for the mad-0 (mad with some restrictions !@ plus some
C     extensions that are considered a comment for mad -> '!@')
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====elements / internal !!!
      DIMENSION TMAP(6,6,MAXELE), VGAM(4,7,MAXELE)
      COMMON /ELEMAP/ TMAP, VGAM

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
C=====In order to compile SPRINT on 'AIn`t uniX' MBUF is just 400 here !
C=====MAD supplies a comand buffer of 4000 !!!
C=====DON'T YOU FORGET TO FORGET ABOUT IBM !!!
      PARAMETER(MAXINP=5000000,MLIN=132,MBUF=400)
      CHARACTER CLINE*(MLIN), CBUF*(MBUF), CSDUM*4

C%%      write(*,*)'rdmad0'

      BELL   = .FALSE.
      BLAT   = .FALSE.
      NELE   = 0

      NLAT   = 0
      NCAV   = 0 
      NTOTSL = 0

      DO 10 I=1,MAXINP
         CBUF = ' '
         IEB   = 0 
         BCONT = .FALSE.

 110     CONTINUE
         CLINE = ' '
         READ(15,'(A)',END=11) CLINE
         CALL STRIP0(CLINE,IE,BCONT,*10,*110)
         CALL STRIPX(CLINE,IE,BCONT,*10,*110)
         CALL STRIPC(CLINE,IE,BCONT,*10,*110)

         IEB1 = IEB + IE
         IF(IEB.GT.MBUF) CALL TOFF('RDMAD0, increase MBUF')

         IF(CLINE(IE:IE).EQ.'&') THEN
            CBUF(IEB+1:IEB1) = CLINE(1:IE-1)
            IEB   = IEB1 - 1
            BCONT = .TRUE.
            GOTO    110
         ELSE
            CBUF(IEB+1:IEB1) = CLINE(1:IE)
            IEB   = IEB1
         ENDIF

         CALL PRSGEN(CBUF,IEB,BELL,BLAT)
         IF(BELL) CALL PRSELL(CBUF,IEB,BELL)
         IF(BLAT) CALL PRSLAT(CBUF,IEB,BLAT)
 10   CONTINUE
      
      CALL TOFF('RDMAD0, increase MAXINP')

 11   CONTINUE

      IF(ITYPE(LATELE(NLAT)).NE.IMR) THEN
         NELE        = NELE + 1
         IF(NELE.GT.MAXELE) CALL TOFF('RDMAD0, increase MAXELE')         
         CNAME(NELE) = 'X______EOL'
         ITYPE(NELE) = IMR
         STRX(NELE)  = 0.D+0
         STRY(NELE)  = 0.D+0
         RLEN(NELE)  = 0.D+0
         TWIS(NELE)  = 0.D+0
         ISEC(NELE)  = 0
         IFLAG(NELE) = 0

         NLAT         = NLAT + 1
         IF(NELE.GT.MAXELE) CALL TOFF('RDMAD0, increase MAXLIST')
         LATELE(NLAT) = NELE 
         PREDRF(NLAT) = 0.D+0
      ENDIF

      IF(NCAV.GT.1) THEN
         WRITE(CSDUM,'(I4)') NCAV
         CALL WARN('RDMAD0, NCAV = '
     +        //CSDUM//' might be time consuming at ramp')
      ENDIF

      WRITE(16,'(2(A,I8/),A,F16.8)')
     +     ' NUMBER OF DIFFERENT ELEMENTS : ',NELE,
     +     ' TOT. NUMBER OF ELEM. IN LAT. : ',NLAT,
     +     ' TOTAL LENGTH :         '        ,DRINGL

      RETURN 
      E N D

      S U B R O U T I N E PRSGEN(CBUF,IEB,BELL,BLAT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  PRSGEN
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     mad0 parser / general stuff
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF),
     +     CKEYW1*(MKEYW), CKEYW2*(MKEYW), CKEYW3*(MKEYW)

C%%      write(*,*)'prsgen'

      CALL MADWRD(CKEYW1,IE1,CBUF,IEB,1)
      CALL MADWRD(CKEYW2,IE2,CBUF,IEB,2)

      IF(.NOT.(BELL.OR.BLAT)) THEN
C=====so we're in mad--intro--section : ignore (but write to 16)

         IF(CKEYW1(IE1:IE1).EQ.':' .AND. IKWELL(CKEYW2).NE.0) THEN
            BELL = .TRUE.
            RETURN
         ELSE
            WRITE(16,'(A)') CBUF(1:IEB)
            RETURN
         ENDIF

      ELSEIF(BELL.AND..NOT.BLAT) THEN
C=====so we're (should be) in the 'ELement List'

         IF(CKEYW1(IE1:IE1).EQ.':') THEN

            IF(CKEYW2.EQ.'SEQUENCE') THEN
C===========switch to lattice list
               BELL = .FALSE.
               BLAT = .TRUE.
               RETURN
            ELSEIF(IKWELL(CKEYW2).NE.0) THEN
C===========still in 'ELement List'
               RETURN
            ENDIF

         ELSE
C========ignore
            WRITE(16,'(A)') CBUF(1:IEB)
            RETURN
         ENDIF

      ELSEIF(BLAT.AND..NOT.BELL) THEN
C=====so we're (should be) in the 'LATice list'
         CALL MADWRD(CKEYW3,IE3,CBUF,IEB,3)

         IF(CKEYW3(1:2).EQ.'AT') THEN
C========still in 'LATice list'
            RETURN
         ELSEIF(CKEYW1.EQ.'ENDSEQUENCE') THEN
C========finished with lattice list
            BELL = .FALSE.
            RETURN
         ELSE
C========ignore
            WRITE(16,'(A)') CBUF(1:IEB)
            RETURN
         ENDIF
         
      ELSE
C=====OOPS!
         CALL TOFF('PRSGEN, internal SPRINT-error')
      ENDIF

      RETURN
      E N D 

      S U B R O U T I N E PRSELL(CBUF,IEB,BELL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  PRSELL
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     mad-0 parser / el_ement l_ist
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====local stuff
      DIMENSION XN(MXSTRN)
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CN*(MKEYW), CTP*(MKEYW)

C%%      write(*,*)'prsell'

      SZ = 0.D+0
      K0 = 0
      K1 = 0
C-BUGFIX-VOGTM-15.02.00 : initialize IPH to 0 sinc some MAD??? routines
C-BUGFIX-VOGTM-15.02.00   don't use it
      IPH= 0
C-FUTURE-VOGTM-15.02.00 : make IPH parameter of  A L L   MAD%%% routines
C-FUTURE-VOGTM-15.02.00   and DON'T initialize it in PRSELL !!!!!!!!!!!!!

      CALL MADWRD(CN ,IE1 ,CBUF,IEB,1)
      CALL MADWRD(CTP,IDUM,CBUF,IEB,2)
      ITPM = IKWELL(CTP)

      GOTO(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     +     21,22,23,24,25,26,27,28,29,30) ITPM
      CALL WARN('PRSELL, unknown MAD-element : '//CN(1:IE1))
      BELL = .FALSE.
      RETURN
 1    CALL MADMRK(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 2    CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 3    CALL MADSBE(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,IPH,XN,*1000)
 4    CALL MADRBE(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,IPH,XN,*1000)
 5    CALL MADQUA(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,IPH,*1000)
 6    CALL MADSEX(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,IPH,*1000)
 7    CALL WARN('PRSELL, no octupoles -> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 8    CALL WARN('PRSELL, only thin quads supported yet: '//
     +     CN(1:IE1))
      CALL MADTMM(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,IPH,*1000)
 9    CALL MADSOL(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,IPH,*1000)
 10   CALL MADKIX(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,IPH,*1000)
 11   CALL MADKIY(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,IPH,*1000)
 12   CALL WARN(
     +     'PRSELL, no general kickers -> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 13   CALL MADCAV(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,IPH,*1000)
 14   CALL WARN(
     +     'PRSELL, no electrostatic septa -> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 15   CALL WARN('PRSELL, no specific monitors -> general: '//CN(1:IE1))
      CALL MADMON(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 16   CALL WARN('PRSELL, no specific monitors -> general: '//CN(1:IE1))
      CALL MADMON(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 17   CALL MADMON(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 18   CALL WARN('PRSELL, no instruments -> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 19   CALL WARN('PRSELL, no collimators -> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 20   CALL WARN('PRSELL, no collimators -> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 21   CALL WARN('PRSELL, no rotations -> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 22   CALL WARN('PRSELL, no rotations -> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 23   CALL WARN('PRSELL, no beambeam -> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 24   CALL WARN(
     +     'PRSELL, no user supp. matrices -> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 25   CALL WARN(
     +     'PRSELL, no lumps of elem.-> drift: '//CN(1:IE1))
      CALL MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*1000)
 26   CALL MADSSN(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,IPH,*1000)
 27   CALL MADASH(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,*1000)
 28   CALL MADMAT(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,*1000)
 29   CALL MADSPM(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,BN,*1000)
 30   CALL MADSRM(CBUF,IEB,ITP,SX,SY,SZ,RL,IS,TW,K0,K1,IFH,BN,IPH,*1000)

 1000 CONTINUE

      IF(IE1.LE.MAELNM+1) THEN
         NELE        = NELE + 1
         CNAME(NELE) = CN(1:INDEX(CN,':')-1)
      ELSE
         CALL TOFF('PRSELL, element name too long')
      ENDIF

      ITYPE(NELE)   = ITP
      STRX(NELE)    = SX
      STRY(NELE)    = SY
      STRZ(NELE)    = SZ
      RLEN(NELE)    = RL
      ISEC(NELE)    = IS
      TWIS(NELE)    = TW
      ISTR0(NELE)   = K0
      ISTR1(NELE)   = K1
      IFLAG(NELE)   = IFH
      IPS(NELE)     = IPH
      STRN(1,NELE)  = XN(1)  
      STRN(2,NELE)  = XN(2)  
      STRN(3,NELE)  = XN(3)  
      STRN(4,NELE)  = XN(4)  

      RETURN
      E N D


      S U B R O U T I N E PRSLAT(CBUF,IEB,BLAT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  PRSLAT
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     mad-0 parser / lat_tice list
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)


      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CPOS*(MKEYW), CDDD*(MKEYW),
     +     CN*(MKEYW),CV*(MKEYW),CDUMPR*14
      SAVE REFMD,EPOS,BIMPLP

C%%      write(*,*)'prslat'

      CALL MADWRD(CDUM,IEE,CBUF,IEB,1)

      IF(CDUM(1:IEE).EQ.'ENDSEQUENCE') THEN
         BLAT = .FALSE.
         DRINGL = DBLE(NSUPER) * EPOS
         RETURN
      ENDIF

      CALL MADWRD(CDUM,IEE,CBUF,IEB,2)

      IF(CDUM(1:IEE).EQ.'SEQUENCE') THEN
         EPOS   = 0.D+0
         BIMPLP = .TRUE. 

         DO 10 I=3,IMADWD(CBUF,IEB)
            CALL MADWRD(CPOS,IEP,CBUF,IEB,I)
            CALL PRSASS(CPOS,IEP,CN,CV,BASS)


            IF(CN.EQ.'REFER') THEN
               IF(CV.EQ.'ENTRY' ) THEN
                  REFMD =  0.D+0
               ELSEIF(CV.EQ.'CENTRE'.OR.CV.EQ.'CENTER') THEN
C==============UK & US spelling... (for Desmond: proper & wrong english )
                  REFMD = 0.5D+0
               ELSEIF(CV.EQ.'EXIT'  ) THEN
                  REFMD = 1.0D+0
               ELSE
                  CALL TOFF(
     +                 'PRSLAT, don''t know where to refer to: '//
     +                 CV//' ?')
               ENDIF
            ELSEIF(CN.EQ.'SUPER') THEN
               NSUPER = INFCS(CV)
            ELSEIF(CN.EQ.'PREDRIFTS') THEN
               IF(CV.EQ.'EXPLICIT') THEN
                  BIMPLP = .FALSE.
               ELSEIF(CV.EQ.'IMPLICIT') THEN
                  BIMPLP = .TRUE.
               ELSE
                  CALL TOFF('PRSLAT, don''t know how to create'//
     +                 ' predrifts: '//CV//' ?')
               ENDIF
            ENDIF

 10      CONTINUE

         RETURN
      ENDIF

      RPOS   = 0.D+0
      PRETST = 0.D+0

      DO 20 I=3,IMADWD(CBUF,IEB)
         CALL MADWRD(CPOS,IEP,CBUF,IEB,I)
         CALL PRSASS(CPOS,IEP,CN,CV,BASS)

         IF(CN.EQ.'AT') THEN
            RPOS   = REFCS(CV)
            CDDD   = CV
         ELSEIF(CN.EQ.'PRDRF') THEN
            PRETST = REFCS(CV)
         ELSE
            CALL TOFF(
     +        'PRSLAT, invalid lattice position descriptor')
         ENDIF
 20   CONTINUE

      LELE = NUMELE(CDUM)
      
      IF(LELE.EQ.0) THEN
         CALL WARN(
     +        'PRSLAT, unknown element : '//CDUM//' -> drift')
         RETURN
C=====maybe later create a new element ??? (like in 'RDMDTW')
C=====also : if known ele. name but with modified parameters (-> use 1-st
C===== madword as new name ) !!!!!!!!!
      ENDIF
      
      NLAT = NLAT + 1
      IF(NLAT.GT.MAXLIST) CALL TOFF(
     +     'PRSLAT, increase MAXLIST')
      LATELE(NLAT) = LELE

Ctest>
C      write(17,'(A,4I6)')'PRSLAT-1: ',nlat,lele,itype(lele),isr
Ctest<

      IF(ITYPE(LELE).EQ.ISR) THEN
         LSLIC(NLAT) = -1
         RLHLP = RLEN(LELE)
      ELSEIF(ISEC(LELE)*IGLBSL.GT.1) THEN
         LSLIC(NLAT) = ISEC(LELE)*IGLBSL
         RLHLP = RLEN(LELE) * DBLE(LSLIC(NLAT))
      ELSE
         LSLIC(NLAT) = 1
         RLHLP = RLEN(LELE)
      ENDIF

Ctest>
C      write(17,'(A,3I6)')'PRSLAT-2: ',nlat,lele,lslic(nlat)
Ctest<

      NTOTSL = NTOTSL + MAX0(1,LSLIC(NLAT))

      IF(BIMPLP) PRETST = RPOS - EPOS - REFMD*RLHLP

      IF(PRETST.LT.-1.D-10) THEN
         WRITE(CDUMPR,'(D14.6)') PRETST
         CALL WARN('PRSLAT, negative predrift'//
     +        CDUMPR//' before '//
     +        CDUM(:IEE)//' at '// CDDD)
      ENDIF

      IF(DABS(PRETST).GT.1.D-9) THEN
         PREDRF(NLAT) = PRETST
      ELSE
         PREDRF(NLAT) = 0.D+0
      ENDIF

      EPOS = EPOS + PREDRF(NLAT) + RLHLP

      IF(ITYPE(LELE).EQ.IRF) NCAV = NCAV + NSUPER

      DRINGL = DBLE(NSUPER) * EPOS

      RETURN
      E N D

      S U B R O U T I N E MADMRK(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADMRK
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type marker from mad-0
C     !!! whole thing still in developement !!!
      
      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF)

C%%      write(*,*)'madmrk'

      ITP  = IMR
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 0                  ! no global slicing possible
      IFH  = 0

C      DO 10 I = 3, IMADWD(CBUF,IEB)
C         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
C         IF(CDUM.EQ.'????') IFH = IBSET(IFH,JB???)
C 10   CONTINUE

      RETURN
      E N D

      S U B R O U T I N E MADDRF(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADDRF
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type drift from mad-0
C     !!! whole thing still in developement !!!
      
      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'maddrf'

      ITP  = IDL
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 1                  ! global slicing possible
      IFH  = 0

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'SLICE') THEN
               IS = INFCS(CV)
            ENDIF
         ENDIF

 10   CONTINUE

      IF(IS*IGLBSL.GT.1) RL = RL/DBLE(IS*IGLBSL)

      RETURN
      E N D

      S U B R O U T I N E MADMON(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADMON
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type monitor from mad-0
C     ==>> treated as drift !!!!
C     !!! whole thing still in developement !!!
      
      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'madmon'

      ITP  = IMO
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 1                  ! global slicing possible
      IFH  = 0

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'SLICE') THEN
               IS = INFCS(CV)
            ENDIF
         ENDIF

 10   CONTINUE

      IF(IS*IGLBSL.GT.1) RL = RL/DBLE(IS*IGLBSL)

      RETURN
      E N D

      S U B R O U T I N E MADSBE(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,IPH,XN, *)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADSBE
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type sector bend or
C     sector combined function mag. from mad-0
C     !!! whole thing still in developement !!!
C=====for combined function magnet :
C     JBSPC8 set -> IPS scales strx
C     JBSPC9 set -> IPS scales stry

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      DIMENSION XN(MXSTRN)
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW),
     +     CSCALE*(MKEYW)

C%%      write(*,*)'madsbe'

      BBQX = .FALSE.
      BBQY = .FALSE.
      ITP  = IBX
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 1
      IFH  = 0
      IPH  = 0
      CSCALE = ' '
      BSWOFF = .FALSE.
      CALL ZERO(MXSTRN,XN)

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN

            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'ANGLE') THEN
               SX = REFCS(CV)
            ELSEIF(CN.EQ.'K1') THEN
               IF(ITP.EQ.IBX) ITP = ICX
               IF(ITP.EQ.IBY) ITP = ICY
               SY  = REFCS(CV)
            ELSEIF(CN.EQ.'SLICE') THEN
               IS = INFCS(CV)
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ELSEIF(CN.EQ.'SCALE') THEN
               CSCALE = CV
            ELSEIF(CN.EQ.'E1') THEN
               XN(1) = REFCS(CV)
            ELSEIF(CN.EQ.'E2') THEN
               XN(2) = REFCS(CV)
            ELSEIF(CN.EQ.'HGAP') THEN
               XN(3) = REFCS(CV)
            ELSEIF(CN.EQ.'FINT') THEN
               XN(4) = REFCS(CV)
            ELSEIF(CN.EQ.'TILT') THEN
               CALL TOFF('MADSBE, no tilted bends yet')
            ENDIF
         ELSE
            IF(CDUM.EQ.'TILT') THEN
               IF(ITP.EQ.IBX) ITP = IBY
               IF(ITP.EQ.ICX) ITP = ICY
            ELSEIF(CDUM.EQ.'NOSPIN') THEN
               IFH = IBSET(IFH,JBNOSP)
            ELSEIF(CDUM.EQ.'FINT') THEN
               XN(4) = 0.5D+0
            ELSEIF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      SY = SY*RL

      IF(ITP.EQ.IBY .OR. ITP.EQ.ICY) SX = -SX

      IF(IS*IGLBSL.GT.1) THEN
         D  =1.0D+0/DBLE(IS*IGLBSL) 
         RL = RL * D
         SX = SX * D
         SY = SY * D
      ENDIF

      IF(IPH.NE.0) THEN
         IF(ITP.EQ.ICX .OR. ITP.EQ.ICY) THEN
            IF(CSCALE.EQ.' ' .OR. CSCALE.EQ.'ALL') THEN
               IFH = IBSET(IFH,JBSPC8)
               IFH = IBSET(IFH,JBSPC9)
            ELSEIF(CSCALE.EQ.'K1') THEN
               IFH = IBSET(IFH,JBSPC9)
            ELSEIF(CSCALE.EQ.'ANGLE') THEN
               IFH = IBSET(IFH,JBSPC8)
            ELSE
               CALL TOFF('MADSBE, invalid SCALE parameter: '//CSCALE)
            ENDIF
         ENDIF

         IF(IPH.LT.0) THEN
            IPH = -IPH
            IFH = IBSET(IFH,JBSNEG)
         ENDIF
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
         SY = 0.D+0
C         TW = 0.D+0
      ENDIF

      RETURN
      E N D 
                                
      S U B R O U T I N E MADRBE(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,IPH,XN, *)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADRBE
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type RECTANGULAR bend or
C     RECTANGULAR combined function mag. from mad-0
C     !!! whole thing still in developement !!!
C=====for combined function magnet :
C     JBSPC8 set -> IPS scales strx
C     JBSPC9 set -> IPS scales stry

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      DIMENSION XN(MXSTRN)
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW),
     +     CSCALE*(MKEYW)

C%%      write(*,*)'madrbe'

      BBQX = .FALSE.
      BBQY = .FALSE.
      ITP  = IRX
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 1
      IFH  = 0
      IPH  = 0
      CSCALE = ' '
      BSWOFF = .FALSE.
      CALL ZERO(MXSTRN,XN)

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'ANGLE') THEN
               SX = REFCS(CV)
            ELSEIF(CN.EQ.'K1') THEN
               IF(ITP.EQ.IRX) ITP = ICH
               IF(ITP.EQ.IRY) ITP = ICV
               SY  = REFCS(CV)
            ELSEIF(CN.EQ.'SLICE') THEN
               IS = INFCS(CV)
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ELSEIF(CN.EQ.'SCALE') THEN
               CSCALE = CV
            ELSEIF(CN.EQ.'E1') THEN
               XN(1) = REFCS(CV)
            ELSEIF(CN.EQ.'E2') THEN
               XN(2) = REFCS(CV)
            ELSEIF(CN.EQ.'HGAP') THEN
               XN(3) = REFCS(CV)
            ELSEIF(CN.EQ.'FINT') THEN
               XN(4) = REFCS(CV)
            ELSEIF(CN.EQ.'TILT') THEN
               CALL TOFF('MADRBE, no tilted bends yet')
            ENDIF
         ELSE
            IF(CDUM.EQ.'TILT') THEN
               IF(ITP.EQ.IRX) ITP = IRY
               IF(ITP.EQ.ICX) ITP = ICV
            ELSEIF(CDUM.EQ.'NOSPIN') THEN
               IFH = IBSET(IFH,JBNOSP)
            ELSEIF(CDUM.EQ.'FINT') THEN
               XN(4) = 0.5D+0
            ELSEIF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      SY = SY*RL

      IF(ITP.EQ.IRY .OR. ITP.EQ.ICY) SX = -SX

      IF(IS*IGLBSL.GT.1) THEN
         D  =1.0D+0/DBLE(IS*IGLBSL) 
         RL = RL * D
         SX = SX * D
         SY = SY * D
      ENDIF

      IF(IPH.NE.0) THEN
         IF(ITP.EQ.ICX .OR. ITP.EQ.ICY) THEN
            IF(CSCALE.EQ.' ' .OR. CSCALE.EQ.'ALL') THEN
               IFH = IBSET(IFH,JBSPC8)
               IFH = IBSET(IFH,JBSPC9)
            ELSEIF(CSCALE.EQ.'K1') THEN
               IFH = IBSET(IFH,JBSPC9)
            ELSEIF(CSCALE.EQ.'ANGLE') THEN
               IFH = IBSET(IFH,JBSPC8)
            ELSE
               CALL TOFF('MADRBE, invalid SCALE parameter: '//CSCALE)
            ENDIF
         ENDIF

         IF(IPH.LT.0) THEN
            IPH = -IPH
            IFH = IBSET(IFH,JBSNEG)
         ENDIF
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
         SY = 0.D+0
C         TW = 0.D+0
      ENDIF

      RETURN
      E N D 
                                
      S U B R O U T I N E MADQUA(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,IPH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADQUA
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type quadrupole from mad-0
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'madqua'

      ITP  = IMQ
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 1
      IFH  = 0
      IPH  = 0
      BSWOFF = .FALSE.
      BUNDY  = .FALSE.          ! yeah, the shoe salesman !

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL  = REFCS(CV)
            ELSEIF(CN.EQ.'K1') THEN
               SX  = REFCS(CV)
            ELSEIF(CN.EQ.'K1L') THEN
               BUNDY = .TRUE.
               SX  = REFCS(CV)
            ELSEIF(CN.EQ.'TILT') THEN
               ITP = IRQ
               TW  = REFCS(CV)
            ELSEIF(CN.EQ.'SLICE') THEN
               IS = INFCS(CV)
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ENDIF
         ELSE
            IF(CDUM.EQ.'TILT') THEN
               ITP = ISQ
            ELSEIF(CDUM.EQ.'NOSPIN') THEN
               IFH = IBSET(IFH,JBNOSP)
            ELSEIF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      IF(BUNDY) THEN 
         CALL WARN('MADQUA, use only MULTIPOLE with K1L --- RTFM !')
         IF(RL.NE.0.D+0) CALL TOFF(
     +        'MADQUA, MULTIPOLES are thin --- RTFM !')
         ITP = ITQ
      ELSE
         SX = SX*RL
      ENDIF

      IF(IS*IGLBSL.GT.1) THEN
         D  =1.0D+0/DBLE(IS*IGLBSL) 
         RL = RL * D
         SX = SX * D
      ENDIF

      IF(IPH.LT.0) THEN
         IPH = -IPH
         IFH = IBSET(IFH,JBSNEG)
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
C         SY = 0.D+0
         TW = 0.D+0
      ENDIF
      
      RETURN
      E N D  

      S U B R O U T I N E MADTMM(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,IPH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADTMM
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type  thin multipole 
C     from mad-0
C     ONLY THIN QUADS SUPPORTED YET !!!!
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'madtmm'
      BTQ  = .FALSE.
      BDL  = .FALSE.
      ITP  = ITQ
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 1
      IFH  = 0
      IPH  = 0
      BSWOFF = .FALSE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'K1L') THEN
               BTQ = .TRUE.
               SX  = REFCS(CV)
            ELSEIF(CN.EQ.'L') THEN
               CALL TOFF(
     +              'MADTMM, no length allowed for thin multipole')
            ELSEIF(CN.EQ.'K0L'.OR.CN.EQ.'K2L'.OR.CN.EQ.'K3L'.OR.
     +              CN.EQ.'K4L'.OR.CN.EQ.'K5L'.OR.CN.EQ.'K6L'.OR.
     +              CN.EQ.'K7L'.OR.CN.EQ.'K8L') THEN
               CALL WARN('MADTMM, no multipoles except quads yet'//
     +              '-> drift')
               ITP = IDL
               BDL = .TRUE.
            ELSEIF(CN.EQ.'T0'.OR.CN.EQ.'T1'.OR.CN.EQ.'T2'.OR.
     +              CN.EQ.'T3'.OR.CN.EQ.'T4'.OR.CN.EQ.'T5'.OR.
     +              CN.EQ.'T6'.OR.CN.EQ.'T7'.OR.CN.EQ.'T8') THEN
               CALL WARN('MADTMM no rotated thin multipoles yet'//
     +              '-> ignore or drift')
               BDL = .TRUE.
            ELSEIF(CN.EQ.'SLICE') THEN
               IS = INFCS(CV)
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ENDIF
         ELSE
            IF(CDUM.EQ.'T0'.OR.CDUM.EQ.'T1'.OR.CDUM.EQ.'T2'.OR.
     +              CDUM.EQ.'T3'.OR.CDUM.EQ.'T4'.OR.CDUM.EQ.'T5'.OR.
     +              CDUM.EQ.'T6'.OR.CDUM.EQ.'T7'.OR.CDUM.EQ.'T8') THEN
               CALL WARN('MADTMM, no rotated thin multipoles yet'//
     +              '-> ignore or drift')
               BDL = .TRUE.
            ELSEIF(CDUM.EQ.'NOSPIN') THEN
               IFH = IBSET(IFH,JBNOSP)
            ELSEIF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      IF(BDL.AND..NOT.BTQ) THEN
         ITP = IDL
         SX  = 0.D+0
      ELSEIF(IS*IGLBSL.GT.1) THEN
         D  = 1.0D+0/DBLE(IS*IGLBSL) 
         SX = SX * D
      ENDIF

      IF(IPH.LT.0) THEN
         IPH = -IPH
         IFH = IBSET(IFH,JBSNEG)
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
         SY = 0.D+0
         TW = 0.D+0
      ENDIF

      RETURN
      E N D  

      S U B R O U T I N E MADSEX(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,IPH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADSEX
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type sextupole from mad-0
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'madsex'

      ITP  = IMS
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 1
      IFH  = 0
      IPH  = 0
      BSWOFF = .FALSE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'K2') THEN
               SX = REFCS(CV)
            ELSEIF(CN.EQ.'SLICE') THEN
               IS = INFCS(CV)
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ELSEIF(CN.EQ.'TILT') THEN
               CALL TOFF('MADSEX, no rotated sextupoles yet')
            ENDIF
         ELSE
            IF(CDUM.EQ.'TILT') THEN
              CALL TOFF('MADSEX, no skew sextupoles yet')
            ELSEIF(CDUM.EQ.'NOSPIN') THEN
               IFH = IBSET(IFH,JBNOSP)
            ELSEIF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      SX = SX*RL
C     hope that'll be right once we use sextupoles......

      IF(IS*IGLBSL.GT.1) THEN
         D  =1.0D+0/DBLE(IS*IGLBSL)
         RL = RL * D
         SX = SX * D
      ENDIF

      IF(IPH.LT.0) THEN
         IPH = -IPH
         IFH = IBSET(IFH,JBSNEG)
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
C         SY = 0.D+0
C         TW = 0.D+0
      ENDIF

      RETURN
      E N D

      S U B R O U T I N E MADSOL(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,IPH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADSOL
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type solenoid from mad-0
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'madsol'

      ITP  = ISO
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 0
      IFH  = 0
      IPH  = 0
      BSWOFF = .FALSE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'KS') THEN
               SX = REFCS(CV)
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ENDIF
         ELSE
            IF(CDUM.EQ.'NOSPIN') THEN
               IFH = IBSET(IFH,JBNOSP)
            ELSEIF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      IF(IPH.LT.0) THEN
         IPH = -IPH
         IFH = IBSET(IFH,JBSNEG)
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
C         SY = 0.D+0
C         TW = 0.D+0
      ENDIF

      RETURN
      E N D

      S U B R O U T I N E MADKIX(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,IPH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADKIX
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type hor. kicker from mad-0
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'madkix'

      ITP  = IKX
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 1
      IFH  = 0
      IPH  = 0
      BSWOFF = .FALSE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'KICK') THEN
               SX = REFCS(CV)
            ELSEIF(CN.EQ.'SLICE') THEN
               IS = INFCS(CV)
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ELSEIF(CN.EQ.'TILT') THEN
               CALL TOFF('MADKIX, no rotated correction coils yet')
            ENDIF
         ELSE
            IF(CDUM.EQ.'TILT') THEN
               CALL TOFF('MADKIX, no rotated correction coils yet')
            ELSEIF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      IF(IS*IGLBSL.GT.1) THEN
         D  =1.0D+0/DBLE(IS*IGLBSL) 
         RL = RL * D
         SX = SX * D
      ENDIF

      IF(IPH.LT.0) THEN
         IPH = -IPH
         IFH = IBSET(IFH,JBSNEG)
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
C         SY = 0.D+0
C         TW = 0.D+0
      ENDIF

      RETURN
      E N D

      S U B R O U T I N E MADKIY(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,IPH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADKIY
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type vert. kicker
C     from mad-0
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'madkiy'

      ITP  = IKY
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 1
      IFH  = 0
      IPH  = 0
      BSWOFF = .FALSE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'KICK') THEN
               SX = REFCS(CV)
            ELSEIF(CN.EQ.'SLICE') THEN
               IS = INFCS(CV)
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ELSEIF(CN.EQ.'TILT') THEN
               CALL TOFF('MADKIY, no rotated correction coils yet')
            ENDIF
         ELSE
            IF(CDUM.EQ.'TILT') THEN
               CALL TOFF('MADKIY, no rotated correction coils yet')
            ELSEIF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      IF(IS*IGLBSL.GT.1) THEN
         D  =1.0D+0/DBLE(IS*IGLBSL) 
         RL = RL * D
         SX = SX * D
      ENDIF

      IF(IPH.LT.0) THEN
         IPH = -IPH
         IFH = IBSET(IFH,JBSNEG)
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
C         SY = 0.D+0
C         TW = 0.D+0
      ENDIF

      RETURN
      E N D

      S U B R O U T I N E MADCAV(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,IPH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADCAV
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type cavity from mad-0
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW),
     +     CSCALE*(MKEYW)

C%%      write(*,*)'madcav'

      ITP  = IRF
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 0
      IFH  = 0
      IPH  = 0
      BSWOFF = .FALSE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'VOLT') THEN
               SX = REFCS(CV)*1.D-3
            ELSEIF(CN.EQ.'FREQ') THEN
               SY = REFCS(CV)*1.D+6
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ELSEIF(CN.EQ.'SCALE') THEN
               CSCALE = CV
            ENDIF
         ELSE
            IF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      IF(IPH.NE.0) THEN
         IF(CSCALE.EQ.' ' .OR. CSCALE.EQ.'ALL') THEN
            IFH = IBSET(IFH,JBSPC8)
            IFH = IBSET(IFH,JBSPC9)
         ELSEIF(CSCALE.EQ.'VOLT') THEN
            IFH = IBSET(IFH,JBSPC8)
         ELSEIF(CSCALE.EQ.'FREQ') THEN
            IFH = IBSET(IFH,JBSPC9)
         ELSE
            CALL TOFF('MADCAV, invalid SCALE parameter: '//CSCALE)
         ENDIF
         
         IF(IPH.LT.0) THEN
            IPH = -IPH
            IFH = IBSET(IFH,JBSNEG)
         ENDIF
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
         SY = 0.D+0
C         TW = 0.D+0
      ENDIF

      RETURN
      E N D

      S U B R O U T I N E MADSSN(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,IPH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADSSN
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type spinrotator from mad-0
C     the existence of the 'SPINROTATOR' keyword is an extension to mad !
C     !!! whole thing still in developement !!!
C=====for spinrotator :
C     JBSPC7 set -> IPS scales phi
C     JBSPC8 set -> IPS scales theta
C     JBSPC9 set -> IPS scales psi

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW),
     +     CSCALE*(MKEYW), CWRN1*16, CWRN2*16

C%%      write(*,*)'madssn'

      ITP  = ISS
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 180.D+0 
      IS   = 0
      IFH  = 0
      IPH  = 0
      CSCALE = ' '
      BSWOFF = .FALSE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'PHI') THEN
               SX = REFCS(CV)
               IF(DABS(SX).GE.360.D+0) THEN
                  WRITE(CWRN1,'(1X,E14.7,1X)') SX
                  IANGHLP = INT(SX/360.D+0)
                  SX = SX - DBLE(IANGHLP)*360.D+0
                  WRITE(CWRN2,'(1X,E14.7,1X)') SX
                  CALL WARN('MADSSN, angle |PHI='//CWRN1//
     +                 '| >=360 -> reduce to'//CWRN2)
               ENDIF
            ELSEIF(CN.EQ.'THETA') THEN
               SY = REFCS(CV)
               IF(DABS(SY).GE.360.D+0) THEN
                  WRITE(CWRN1,'(1X,E14.7,1X)') SY
                  IANGHLP = INT(SY/360.D+0)
                  SY = SY - DBLE(IANGHLP)*360.D+0
                  WRITE(CWRN2,'(1X,E14.7,1X)') SY
                  CALL WARN('MADSSN, angle |THETA='//CWRN1//
     +                 '| >=360 -> reduce to'//CWRN2)
               ENDIF
            ELSEIF(CN.EQ.'PSI') THEN
               TW = REFCS(CV)
               IF(DABS(TW).GE.360.D+0) THEN
                  WRITE(CWRN1,'(1X,E14.7,1X)') TW
                  IANGHLP = INT(TW/360.D+0)
                  TW = TW - DBLE(IANGHLP)*360.D+0
                  WRITE(CWRN2,'(1X,E14.7,1X)') TW
                  CALL WARN('MADSSN, angle |PSI='//CWRN1//
     +                 '| >=360 -> reduce to'//CWRN2)
               ENDIF
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ELSEIF(CN.EQ.'SCALE') THEN
               CSCALE = CV
            ENDIF

         ELSE
            IF(CDUM.EQ.'FREE') THEN
                ITP = ISF
C               SX  = 0.D+0
C               SY  = 0.D+0
C               TW  = 0.D+0
            ELSEIF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      IF(IPH.NE.0 .AND. ITP.EQ.ISS) THEN
         IF(CSCALE.EQ.'ALL') THEN
            IFH = IBSET(IFH,JBSPC7)
            IFH = IBSET(IFH,JBSPC8)
            IFH = IBSET(IFH,JBSPC9)
         ELSEIF(CSCALE.EQ.'PHI') THEN
            IFH = IBSET(IFH,JBSPC7)
         ELSEIF(CSCALE.EQ.'THETA') THEN
            IFH = IBSET(IFH,JBSPC8)
         ELSEIF(CSCALE.EQ.'PSI') THEN
            IFH = IBSET(IFH,JBSPC9)
         ENDIF
         
      ENDIF
      

      IF(IPH.LT.0) THEN
         IPH = -IPH
         IFH = IBSET(IFH,JBSNEG)
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
         SY = 0.D+0
         TW = 0.D+0
      ENDIF

      RETURN
      E N D

      S U B R O U T I N E MADSPM(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADSPM
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type
C     'insantanious linear spin match' from mad-0
C     !!! whole thing still in developement !!!
      
      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'madspm'

      ITP  = ISM
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 0                  ! no global slicing possible
      IFH  = 0
      BSWOFF = .FALSE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)
         IF(BASS) THEN
            CALL TOFF('MADSPM no parameters for spin match element')
         ELSE
            IF(CDUM.EQ.'OFF') THEN
               CALL TOFF('MADSPM, can''t yet switch off spin match')
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      RETURN
      E N D

      S U B R O U T I N E MADASH(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADASH
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type
C     'shifter of orbital phases' from mad-0
C     !!! whole thing still in developement !!!
      
      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'madash'

      ITP  = IAS
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 0                  ! no global slicing possible
      IFH  = 0
      BSWOFF = .FALSE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'PHX') THEN
               SX = REFCS(CV)
               IF(DABS(SX).GE.1.D+0) CALL WARN('MADASH, angle > 1')
               SX = SX-INT(SX)
            ELSEIF(CN.EQ.'PHY') THEN
               SY = REFCS(CV)
               IF(DABS(SY).GE.1.D+0) CALL WARN('MADASH, angle > 1')
               SY = SY-INT(SY)
            ELSEIF(CN.EQ.'PHT') THEN
               TW = REFCS(CV)
               IF(DABS(TW).GE.1.D+0) CALL WARN('MADASH, angle > 1')
               TW = TW-INT(TW)
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ENDIF
            IF(CDUM.EQ.'OFF') THEN
               CALL TOFF('MADSPM, can''t yet switch off angle shifter')
               BSWOFF = .TRUE.
            ENDIF
         ENDIF
 10   CONTINUE

      RETURN
      E N D

      S U B R O U T I N E MADMAT(CBUF,IEB,ITP,SX,SY,RL,IS,TW,IFH,
     +     BSWOFF,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADMAT
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     fills internal input arrays for an element of type
C     'specified 2X2 matrix' from mad-0
C     !!! whole thing still in developement !!!
      
      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW)

C%%      write(*,*)'madmat'

      ITP  = IMA
      SX   = 0.D+0
      SY   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0
      IS   = 0                  ! no global slicing possible
      IFH  = 0
      BSWOFF = .FALSE.
      BPLANE = .TRUE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
C               RL = REFCS(CV)
               CALL TOFF('MADMAT, no length allowed for TWISSMAT')
            ELSEIF(CN.EQ.'PLANE') THEN
               IF(CV.EQ.'X') THEN
                  IFH = IBSET(IFH,JBSPC7)
                  BPLANE = .FALSE.
               ELSEIF(CV.EQ.'Y') THEN
                  IFH = IBSET(IFH,JBSPC8)
                  BPLANE = .FALSE.
               ELSEIF(CV.EQ.'Z') THEN
                  IFH = IBSET(IFH,JBSPC9)
                  BPLANE = .FALSE.
               ELSE
                  CALL TOFF('MADMAT, plane must be X,Y,Z')
               ENDIF
            ELSEIF(CN.EQ.'ALPHA') THEN
               SX = REFCS(CV)
            ELSEIF(CN.EQ.'BETA') THEN
               SY = REFCS(CV)
               IF(DABS(SY).LE.0.D+0) CALL WARN('MADMAT, beta <= 0')
            ELSEIF(CN.EQ.'MU') THEN
               TW = REFCS(CV)
               IF(DABS(TW).GE.1.D+0) CALL WARN('MADMAT, angle > 1')
               TW = TW-DBLE(INT(TW))
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ENDIF
            IF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF
 10   CONTINUE


      IF(BSWOFF) TW = 0.D+0 
C     cos(0)=1, sin(0)=0 => sets matrix to identity independent of alpha/beta
      IF(BPLANE) IFH = IBSET(IFH,JBSPC8) 
C     default plane is vertical 

      RETURN
      E N D


      S U B R O U T I N E MADSRM(CBUF,IEB,ITP,SX,SY,SZ,RL,IS,TW,K0,K1,
     +     IFH,BSWOFF,IPH,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADSRM
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     fills internal input arrays for an element of type SRMCELL from mad-0
C     the existence of the 'SRMCELL' keyword is an extension to mad !
C     !!! whole thing still in developement !!!
C=====for spinrotator :
C     JBSPC7 set -> IPS scales MUY
C     JBSPC8 set -> IPS scales ANGLE
C     JBSPC9 set -> IPS scales MUEPS

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CBUF*(MBUF), CDUM*(MKEYW), CN*(MKEYW), CV*(MKEYW),
     +     CSCALE*(MKEYW), CWRN1*16, CWRN2*16

C%%      write(*,*)'madsrm'

      ITP  = ISR
      BSRM = .TRUE.
      SX   = 0.D+0
      SY   = 0.D+0
      SZ   = 0.D+0
      RL   = 0.D+0
      TW   = 0.D+0 
      IS   = 0
      K0   = 0
      K0   = 0
      IFH  = 0
      IPH  = 0
      CSCALE = ' '
      BSWOFF = .FALSE.

      DO 10 I = 3, IMADWD(CBUF,IEB)
         CALL MADWRD(CDUM,IE1,CBUF,IEB,I)
         CALL PRSASS(CDUM,IE1,CN,CV,BASS)

         IF(BASS) THEN
            IF(CN.EQ.'L') THEN
               RL = REFCS(CV)
            ELSEIF(CN.EQ.'MUY') THEN
               SX = REFCS(CV)
               IF(DABS(SX).GE.1.D+0) CALL WARN('MADSRM, MUY > 1')
               SX = SX - DBLE(INT(SX))
            ELSEIF(CN.EQ.'ANGLE') THEN
               SY = REFCS(CV)
               IF(DABS(SY).GT.360.D+0) THEN
                  WRITE(CWRN1,'(1X,E14.7,1X)') SY
                  SY = SY - DBLE(INT(SY/360.D+0))*360.D+0
                  WRITE(CWRN2,'(1X,E14.7,1X)') SY
                  CALL WARN('MADSRM, angle |ANGLE='//CWRN1//
     +                 '| >360 -> reduce to'//CWRN2)
               ENDIF
            ELSEIF(CN.EQ.'EPS0') THEN
               SZ = REFCS(CV)
               IF(SZ.LT.0) CALL WARN('MADSRM, EPS0 < 0')
               SZ = DABS(SZ)
            ELSEIF(CN.EQ.'MUEPS') THEN
               TW = REFCS(CV)
               IF(DABS(TW).GE.1.D+0) CALL WARN('MADSRM, MUEPS > 1')
               TW = TW - DBLE(INT(TW))
            ELSEIF(CN.EQ.'OFFSET') THEN
               K0 = INFCS(CV)
            ELSEIF(CN.EQ.'ORDER') THEN
               K1 = REFCS(CV)
            ELSEIF(CN.EQ.'IFLAG') THEN
               IFH = INFCS(CV)
            ELSEIF(CN.EQ.'IPS') THEN
               IPH = INFCS(CV)
            ELSEIF(CN.EQ.'SCALE') THEN
               CSCALE = CV
            ENDIF
         ELSE
            IF(CDUM.EQ.'OFF') THEN
               BSWOFF = .TRUE.
            ENDIF
         ENDIF

 10   CONTINUE

      IF(IPH.NE.0) THEN
         
         IF(CSCALE.EQ.'ALL') THEN
            IFH = IBSET(IFH,JBSPC7)
            IFH = IBSET(IFH,JBSPC8)
            IFH = IBSET(IFH,JBSPC9)
            CALL WARN(
     +           'MADSRM, can only scale MUY,ANGLE,MUEPS -> SET')
         ELSEIF(CSCALE.EQ.'MUY') THEN
            IFH = IBSET(IFH,JBSPC7)
         ELSEIF(CSCALE.EQ.'ANGLE') THEN
            IFH = IBSET(IFH,JBSPC8)
         ELSEIF(CSCALE.EQ.'MUEPS') THEN
            IFH = IBSET(IFH,JBSPC9)
         ELSE
            CALL WARN(
     +           'MADSRM, can only scale MUY,ANGLE,MUEPS -> ignore')
         ENDIF
         
      ENDIF
      
      IF(IPH.LT.0) THEN
         IPH = -IPH
         IFH = IBSET(IFH,JBSNEG)
      ENDIF

      IF(BSWOFF) THEN
         SX = 0.D+0
         SY = 0.D+0
         SZ = 0.D+0
         TW = 0.D+0
         k0 = 0
         k1 = 0
      ENDIF

Ctest>
C      write(17,'(A,4D12.4,2I6)'  ) 'MADSRM: ',sx,sy,sz,tw,k0,k1
C      write(17,'(A,3I6,D12.4,L2)') '        ',is,iph,ifh,rl,bswoff
Ctest<

      RETURN
      E N D

      S U B R O U T I N E EXTMAD(CBUF,ITYPE,IPS,ISL,IFL,SY,BXT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  EXTMAD
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     auxillary function that writes the specific mad0--extensions
C     to CBUF  o n l y  if they are relevant (i.g. more than one slice)
C     can easily be extended (hopefully......)

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

      PARAMETER(MBUF=400)
      CHARACTER CBUF*(MBUF)

C%%      write(*,*)'extmad'

      BXT2  = .FALSE.
      BPS   = .FALSE.
      BSL   = .FALSE.
      BFL   = .FALSE.
      BK1LY = .FALSE.
      IE    = ICSEND(CBUF)

C=====checking for extensions other than default
      IF(IPS.NE.0) THEN
C=====power supply
         BXT2 = .TRUE.
         BPS  = .TRUE.
      ENDIF

      IF(ISL.GT.1) THEN
C=====slices
         BXT2 = .TRUE.
         BSL  = .TRUE.
      ENDIF

      IF(IFL.NE.0) THEN
C=====for flags
         BXT2 = .TRUE.
         BFL  = .TRUE.
      ENDIF

      IF(SY.NE.0.D+0) THEN
C=====for flags
         BXT2   = .TRUE.
         BK1LY  = .TRUE.
      ENDIF

C=====there are extensions and this line contains standard mad also
      IF(BXT.AND.BXT2) THEN
         WRITE(CBUF(IE+1:IE+3),'(A3)') ' !@'
         IE = IE + 3
      ENDIF

C=====writing out if not set to default 
      IF(BPS) THEN
         IF(BTEST(IFL,JBSNEG)) THEN
            IPS0 = -IPS
         ELSE
            IPS0 = IPS
         ENDIF

         WRITE(CBUF(IE+1:IE+9),'(A6,I3)') ', IPS=',IPS0
         IE = IE + 9            
      ENDIF

      IF(BSL) THEN
         WRITE(CBUF(IE+1:IE+11),'(A8,I3)') ', SLICE=',ISL
         IE = IE + 11            
      ENDIF
      
      IF(BFL) THEN
C         WRITE(CBUF(IE+1:IE+18),'(A8,I10)') ', IFLAG=',IFL
C         IE = IE + 18

C========type independent flags (bit 0 to ....)
         IF(BTEST(IFL,JBNOSP)) THEN
            WRITE(CBUF(IE+1:IE+8),'(A8)') ', NOSPIN'
            IE = IE + 8
         ENDIF

         IF(BTEST(IFL,JBSNEG)) THEN
            IPS = -IPS
         ENDIF

C========type dependent flags (bit ... to 30)
         IF    (ITYPE.EQ.ICX) THEN
            B8 = BTEST(IFL,JBSPC8)
            B9 = BTEST(IFL,JBSPC9)            

            IF    (B8.AND..NOT.B9) THEN
               WRITE(CBUF(IE+1:IE+13),'(A13)') ', SCALE=ANGLE'
               IE = IE + 13               
            ELSEIF(B9.AND..NOT.B8) THEN
               WRITE(CBUF(IE+1:IE+10),'(A10)') ', SCALE=K1'
               IE = IE + 10               
            ENDIF

         ELSEIF(ITYPE.EQ.ICY) THEN
            B8 = BTEST(IFL,JBSPC8)
            B9 = BTEST(IFL,JBSPC9)            

            IF    (B8.AND..NOT.B9) THEN
               WRITE(CBUF(IE+1:IE+13),'(A13)') ', SCALE=ANGLE'
               IE = IE + 13               
            ELSEIF(B9.AND..NOT.B8) THEN
               WRITE(CBUF(IE+1:IE+10),'(A10)') ', SCALE=K1'
               IE = IE + 10               
            ENDIF

         ELSEIF(ITYPE.EQ.IRF) THEN
            B8 = BTEST(IFL,JBSPC8)
            B9 = BTEST(IFL,JBSPC9)            

            IF    (B8.AND..NOT.B9) THEN
               WRITE(CBUF(IE+1:IE+12),'(A12)') ', SCALE=VOLT'
               IE = IE + 12               
            ELSEIF(B9.AND..NOT.B8) THEN
               WRITE(CBUF(IE+1:IE+12),'(A12)') ', SCALE=FREQ'
               IE = IE + 12               
            ENDIF

         ELSEIF(ITYPE.EQ.ISS) THEN

            IF    (BTEST(IFL,JBSPC7)) THEN
               WRITE(CBUF(IE+1:IE+11),'(A11)') ', SCALE=PHI'
               IE = IE + 11               
            ELSEIF(BTEST(IFL,JBSPC8)) THEN
               WRITE(CBUF(IE+1:IE+13),'(A13)') ', SCALE=THETA'
               IE = IE + 13
            ELSEIF(BTEST(IFL,JBSPC9)) THEN
               WRITE(CBUF(IE+1:IE+11),'(A11)') ', SCALE=PSI'
               IE = IE + 11 
            ENDIF

         ELSEIF(ITYPE.EQ.IMA) THEN

            IF    (BTEST(IFL,JBSPC7)) THEN
               WRITE(CBUF(IE+1:IE+9),'(A9)') ', PLANE=X'
               IE = IE + 9               
            ELSEIF(BTEST(IFL,JBSPC8)) THEN
               WRITE(CBUF(IE+1:IE+9),'(A9)') ', PLANE=Y'
               IE = IE + 9
            ELSEIF(BTEST(IFL,JBSPC9)) THEN
               WRITE(CBUF(IE+1:IE+9),'(A9)') ', PLANE=Z'
               IE = IE + 9 
            ENDIF

         ENDIF

      ENDIF

      IF(BK1LY) THEN
         IF(ITYPE.EQ.ITQ) THEN
            WRITE(CBUF(IE+1:IE+22),'(A7,E15.9)') ', K1LY=',SY 
            IE = IE + 22
         ENDIF
      ENDIF

      E N D

      INTEGER F U N C T I O N IKWELL(CKW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  IKWELL
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     interpretes a type keyword in a elelment list to return an integer
C     works the same as 'IFCMD'

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====local stuff
      PARAMETER(MKEYW=32,NKWLA=30)
      CHARACTER CKW*(MKEYW)
      CHARACTER CKWLS(NKWLA)*(MKEYW)

      DATA (CKWLS(I),I=1,NKWLA) /'MARKER','DRIFT','SBEND','RBEND',
     +     'QUADRUPOLE','SEXTUPOLE','OCTUPOLE','MULTIPOLE','SOLENOID',
     +     'HKICKER','VKICKER','KICKER','RFCAVITY','ELSEPARATOR',
     +     'HMONITOR','VMONITOR','MONITOR','INSTRUMENT',
     +     'ECOLLIMATOR','RCOLLIMATOR','YROT','SROT','BEAMBEAM',
     +     'MATRIX','LUMP','SPINROTATOR','SHIFTER','TWISSMAT',
     +     'SPINMATCH','SRMCELL'/

C%%      write(*,*)'ikwell'

      IKWELL = 0
      
      DO 10 I=1,NKWLA
         IF(CKW.EQ.CKWLS(I)) THEN
            IKWELL = I
            RETURN
         ENDIF
 10   CONTINUE

C=====unfortunately 'LINE' is an  a s s i g n m e n t  'LINE=(x,y,z,...)'
      IF(CKW(1:4).EQ.'LINE') CALL TOFF(
     +     'IKWELL, no LINEs supported yet')

      RETURN
      E N D

      CHARACTER*2 F U N C T I O N CREGMR(BALPH,NNEW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  CREGMR
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     auxillary function Create_REGion_MaRk : returns char*2
C     coresponding to the NNEWest .TRUE. item in logical array BALPH
C     => So no already existing names from the typelist can be
C        choosen for mad--sequence--element names

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MALPH=26)

      DIMENSION BALPH(MALPH,MALPH)
      CHARACTER CALPH*(MALPH)

      DATA CALPH /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      CREGMR = '!!'
      NOLD   = 0

      DO 10 I=1,MALPH
         DO 110 J=1,MALPH
            IF(BALPH(I,J)) NOLD = NOLD +1

            IF(NOLD.EQ.NNEW) THEN
               CREGMR = CALPH(I:I)//CALPH(J:J)
               RETURN
            ENDIF

 110     CONTINUE
 10   CONTINUE

      E N D

      S U B R O U T I N E WRM80X(IUMAD0,CBUF)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  WRM80X
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     formates CBUF with respect to mad & mad0 convetions, i.e. (\n=newrecord):
C     - no more than 80 cols. per line (RDMAD0 could read up to 132 !)
C     - continuation a`la : '' FistLine & \n SecondLine & \n ThirdLine.....''
C     - mad0--extensions look like comments to mad : ''madbla !@ mad0bla''
C     - if a continuation line starts with a mad0--extension then the
C       continuation character '&' has to be mad-commented :
C       '' FirstLine !@ & \n !@ SecondLine.... ''
C     - indention '(3X)' before continuation lines
 
      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MBUF=400)
      CHARACTER CBUF*(MBUF)

C%%      write(*,*)'wrm80x'

C========preparing for 80 cols. &-continuation 
C========and !@-extension[once per line!]
      IX  = INDEX(CBUF,'!@')
      IF(IX.EQ.0) IX = MBUF + 2
      IB  = 1
      IE  = ICSEND(CBUF)
      IE1 = MIN(80,MBUF)
      
C========re-entry for continuation
  10  CONTINUE
      
      IF(IE.GT.IE1) THEN
C========too long => cut !
         IE1 = IB + 72
         
         DO 110 J=IE1,IB,-1
C===========searching (backwards) for a good place to cut
            IF    (CBUF(J:J).EQ.',') THEN
C==============cut after ','
               IE1 = J
               GOTO 111
            ELSEIF(CBUF(J:J).EQ.'!') THEN
C==============cut before '!@'
               IE1 = J-1
               GOTO 111
            ENDIF
            
  110    CONTINUE
         CALL TOFF('WRMAD0, internal error')
  111    CONTINUE
         
         IF(IB.EQ.1) THEN
            IF(IX.EQ.IE1+1 .OR. IX.EQ.IE1+2) THEN
C==============cont. only for extension => std.-mad : NO cont.
               WRITE(IUMAD0,'(2A)') CBUF(IB:IE1),' !@&'
            ELSE
C==============cont. of std.-mad stuff
               WRITE(IUMAD0,'(2A)') CBUF(IB:IE1),' &'
            ENDIF
         ELSE
C===========``\setlength{\MadContIndent}{3ex}'' [speak LaTeX or die!]
            IF    (IX.EQ.IE1+1 .OR. IX.EQ.IE1+2) THEN
C==============cont. only for extension => std.-mad : NO cont.
               WRITE(IUMAD0,'(3X,2A)')      CBUF(IB:IE1),' !@&'
            ELSEIF(IX.LT.IB) THEN
C==============cont. of mad-extension => std-mad : comment
               WRITE(IUMAD0,'(3X,3A)') '!@',CBUF(IB:IE1),' &'
            ELSE
C==============cont. of std.-mad stuff
               WRITE(IUMAD0,'(3X,2A)')      CBUF(IB:IE1),' &'
            ENDIF
         ENDIF
         
         IB  = IE1 + 1
         IE1 = IE1 + 73
         
C===========continuate this one !
         GOTO 10
         
      ELSE
C========fits into 80 cols ! => no (FURTHER!!!!) cont.
         IF    (IB.EQ.1 ) THEN
C===========just one line with .LE. 80 chars
            WRITE(IUMAD0,    '(A)')      CBUF(IB:IE)
         ELSEIF(IX.LT.IB) THEN
C===========cont. of mad-extension => std-mad : comment
            WRITE(IUMAD0,'(3X,2A)') '!@',CBUF(IB:IE)
         ELSE
C==============cont. of std.-mad stuff
            WRITE(IUMAD0, '(3X,A)')      CBUF(IB:IE)
         ENDIF
         
      ENDIF
      
      E N D


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C29                                                                           C
C@ PROGRAM=SPRINT, MODULE=MAD_TWISS, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      S U B R O U T I N E RDMDTW1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  RDMDTW // o r i g i n a l taken from SPRINT
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads the mad output (twiss/survey/chrom files) as lattice input
C     ignores all the twiss parameters at the moment
C     works fine with unix & twiss/.../-files created under unix
C     trips on 'VMS' with files generated under unix and ftp-ed to the 'VAX'
C     maybe that's due to the fact that the first char of a line is used
C     quite normally under unix but might be either ignored or interpreted
C     as control char on real computers ?????
C     gotta fix that myself, i think ! don't care too much about that ! ;-)
C     who volunteers to implement a half--way consistent MAD-TWISS
C     implementation of srm-cell ?????????????  (me ???)

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====elements / internal !!!
      DIMENSION TMAP(6,6,MAXELE), VGAM(4,7,MAXELE)
      COMMON /ELEMAP/ TMAP, VGAM

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MAXINP=1000000,MLIN=80,DSKPRC=1.D-3)
      CHARACTER CLINE*(MLIN), CDUMPR*14, CDUMAT*14
      CHARACTER*8 CPRGV,CDATA,CDATM,CTIMM,CJBNM,CKW,CNM
      CHARACTER*4 CSDUM

      DSKEW(I,X) = DABS(PI/DBLE(I)-DABS(X))

C%%      write(*,*)'rdmdtw'

      NELE   = 0
      NLAT   = 0
      EPOS   = 0.D+0
      NCAV   = 0
      NTOTSL = 0

      WRITE(16,'(A/)')' READ MAD TWISS- OR SURVEY- TABLE :'
      READ(15,'(5A8,I8,L8,I8)')CPRGV,CDATA,CDATM,CTIMM,CJBNM,
     +     NSUPER,BASYM,NREC
      WRITE(16,'(5(A,A8/),A,I8/A,L8/A,I8/)')
     +     ' MAD VERSION      : ',CPRGV,
     +     ' DATA TYPE        : ',CDATA,
     +     ' CREATION DATE    : ',CDATM,
     +     ' CREATION TIME    : ',CTIMM,
     +     ' JOB NAME         : ',CJBNM,
     +     ' SUPERPERIOD      : ',NSUPER,
     +     ' ASYM. FLAG       : ',BASYM,
     +     ' NUM OF POS. RECS : ',NREC

      READ(15,'(A80)')CLINE
      WRITE(16,'(A,A80)')' TITLE : ',CLINE
      READ(15,'(2A8,A4,F12.6,3E16.9)')CKW,CNM,CT,RL,AN,DK1,DK2

      IF(CNM.EQ.'INITIAL') THEN
         READ(15,'(A80)')CLINE

         IF(CDATA.EQ.'SURVEY') THEN
            READ(15,'(4E16.9/3E16.9)',END=99)
     +           XMAD,YMAD,ZMAD,EPOS,
     +           THEMAD,PHIMAD,PSIMAD
         ELSEIF(CDATA.EQ.'TWISS') THEN
            READ(15,'(5E16.9/5E16.9/5E16.9))',END=99)
     +           ALX,BEX,DMX,DIX,DPX,
     +           ALY,BEY,DMY,DIY,DPY,
     +           XMAD,AMAD,YMAD,DBMD,EPOS
         ELSEIF(CDATA.EQ.'CHROM') THEN
            READ(15,'(5E16.9/5E16.9/5E16.9))',END=99)
     +           WX,PHX,DMX,DDX,DDPX
     +           WY,PHY,DMY,DDY,DDPY
     +           XMAD,AMAD,YMAD,DBMD,EPOS     
         ENDIF

         NREC = NREC - 1
      ELSE
         BACKSPACE(15)
      ENDIF

      DO 10 I=1,NREC
         READ(15,'(2A8,A4,F12.6,3E16.9/5E16.9)',END=99)
     +        CKW,CNM,CT,RL,AN,DK1,DK2
     +        ,TW,FR,VO,DUM1,DUM2

         IF(CDATA.EQ.'SURVEY') THEN
            READ(15,'(4E16.9/3E16.9)',END=99)
     +           XMAD,YMAD,ZMAD,AT,
     +           THEMAD,PHIMAD,PSIMAD
         ELSEIF(CDATA.EQ.'TWISS') THEN
            READ(15,'(5E16.9/5E16.9/5E16.9))',END=99)
     +           ALX,BEX,DMX,DIX,DPX,
     +           ALY,BEY,DMY,DIY,DPY,
     +           XMAD,AMAD,YMAD,DBMD,AT
         ELSEIF(CDATA.EQ.'CHROM') THEN
            READ(15,'(5E16.9/5E16.9/5E16.9))',END=99)
     +           WX,PHX,DMX,DDX,DDPX,
     +           WY,PHY,DMY,DDY,DDPY,
     +           XMAD,AMAD,YMAD,DBMD,AT     
         ENDIF

         LELE = NUMELE(CNM)

         IF(LELE.EQ.0) THEN
            NELE = NELE + 1
            IF(NELE.GT.MAXELE) CALL TOFF('RDMDTW, increase MAXELE')
            CNAME(NELE) = CNM
            ITYPE(NELE) = 0
            STRX(NELE)  = 0.D+0
            STRY(NELE)  = 0.D+0
            RLEN(NELE)  = 0.D+0
            TWIS(NELE)  = 0.D+0
            ISEC(NELE)  = 1
            IFLAG(NELE) = 0
            IPS(NELE)   = 0

            IF    (CKW.EQ.'MARKER') THEN
               ITYPE(NELE) = IMR
               ISEC (NELE) = 0
            ELSEIF(CKW.EQ.'DRIFT') THEN
               ITYPE(NELE) = IDL
               RLEN(NELE)  = RL
            ELSEIF(CKW.EQ.'RBEND') THEN
               ITDUM       = IRX
               RLEN(NELE)  = RL
               STRX(NELE)  = AN

               IF(DSKEW(2,TW).LT.DSKPRC) THEN
                  ITDUM = IRY
                  STRX(NELE) = -STRX(NELE)
               ELSEIF(TW.NE.0.D+0) THEN
                  CALL TOFF('RDMDTW, no rotated bends yet !')
               ENDIF

               IF(DK1.NE.0.D+0) THEN
                  STRY(NELE) = DK1 * RL
                  IF(ITDUM.EQ.IRY) THEN
                     ITDUM = ICV
                  ELSE
                     ITDUM = ICH
                  ENDIF
               ENDIF

               ITYPE(NELE) = ITDUM                     
            ELSEIF(CKW.EQ.'SBEND') THEN
               ITDUM       = IBX
               RLEN(NELE)  = RL
               STRX(NELE)  = AN

               IF(DSKEW(2,TW).LT.DSKPRC) THEN
                  ITDUM = IBY
                  STRX(NELE) = -STRX(NELE)
               ELSEIF(TW.NE.0.D+0) THEN
                  CALL TOFF('RDMDTW, no rotated bends yet !')
               ENDIF

               IF(DK1.NE.0.D+0) THEN
                  STRY(NELE) = DK1 * RL
                  IF(ITDUM.EQ.IBY) THEN
                     ITDUM = ICY
                  ELSE
                     ITDUM = ICX
                  ENDIF
               ENDIF

               ITYPE(NELE) = ITDUM    
            ELSEIF(CKW.EQ.'QUADRUPO') THEN
               ITDUM       = IMQ
               RLEN(NELE)  = RL
               STRX(NELE)  = DK1 * RL

               IF(DSKEW(4,TW).LT.DSKPRC) THEN
                  ITDUM = ISQ
                  IF(TW.LT.0.D+0) STRX(NELE) = -STRX(NELE)
C     right convention?
               ELSEIF(TW.NE.0.D+0) THEN
                  ITDUM = IRQ
                  TWIS(NELE) = TW
               ENDIF

               ITYPE(NELE) = ITDUM  
            ELSEIF(CKW.EQ.'SEXTUPOL') THEN

               IF(TW.NE.0.D+0) THEN
                  CALL WARN(
     +                 'RDMDTW, no rotated sextupoles yet -> drift')
                  NELE = NELE - 1
                  GOTO 10       ! skip
               ENDIF

               ITYPE(NELE) = IMS
               RLEN(NELE)  = RL
               STRX(NELE)  = DK2 * RL
            ELSEIF(CKW.EQ.'OCTUPOLE') THEN
               CALL WARN('RDMDTW, no octupoles yet -> drift')
               NELE = NELE - 1
               GOTO 10          ! skip
            ELSEIF(CKW.EQ.'MULTIPOL') THEN

               IF(DK1.NE.0.D+0) THEN
                  ITYPE(NELE) = ITQ
                  STRX(NELE)  = DK1
                  IF(AN .NE. 0.D+0) CALL TOFF(
     +                 'RMDMTW, no thin lense quad + dipole yet')
                  IF(VO .NE. 0.D+0) CALL TOFF(
     +                 'RMDMTW, no rotated thin lense quad yet')
                  IF(DABS(DK2)+DABS(FR) .NE. 0.D+0) CALL WARN(
     +                 'RMDMTW, no thin lense quad + h.o.multip. yet')
               ELSE
                  IF(AN .NE. 0.D+0) CALL TOFF(
     +                 'RMDMTW, no multipoles containing dipols yet') 
                  CALL WARN('RDMDTW, no multipoles yet -> drift')
                  NELE = NELE - 1
                  GOTO 10       ! skip
               ENDIF

            ELSEIF(CKW.EQ.'SOLENOID') THEN
               ITYPE(NELE) = ISO
               RLEN(NELE)  = RL
               STRX(NELE)  = FR
               ISEC (NELE) = 0
            ELSEIF(CKW.EQ.'RFCAVITY') THEN
               ITYPE(NELE) = IRF
               RLEN(NELE)  = RL
               STRX(NELE)  = FR
               STRY(NELE)  = VO
               ISEC (NELE) = 0
            ELSEIF(CKW.EQ.'ELSEPARA') THEN
               CALL WARN('RDMDTW, no el. sep. yet -> drift')
               NELE = NELE - 1
               GOTO 10          ! skip
            ELSEIF(CKW.EQ.'HKICK') THEN

               IF(TW.NE.0.D+0) THEN
                  CALL WARN(
     +                 'RDMDTW, no rotated kickers yet -> drift')
                  NELE = NELE - 1
                  GOTO 10       ! skip
               ENDIF

               ITYPE(NELE) = IKX
               RLEN(NELE)  = RL
               STRX(NELE)  = FR
            ELSEIF(CKW.EQ.'VKICK') THEN

               IF(TW.NE.0.D+0) THEN
                  CALL WARN(
     +                 'RDMDTW, no rotated kickers yet -> drift')
                  NELE = NELE - 1
                  GOTO 10       ! skip
               ENDIF

               ITYPE(NELE) = IKY
               RLEN(NELE)  = RL
               STRX(NELE)  = FR
            ELSEIF(CKW.EQ.'RCOLLIMA') THEN
               CALL WARN('RDMDTW, collimators ignored -> drift')
               NELE = NELE - 1
               GOTO 10          ! skip
            ELSEIF(CKW.EQ.'SROT') THEN
               CALL WARN('RDMDTW, no long. rotation yet -> drift')
               NELE = NELE - 1
               GOTO 10          ! skip
            ELSEIF(CKW.EQ.'YROT') THEN
               CALL WARN('RDMDTW, no vert. rotation yet -> drift')
               NELE = NELE - 1
               GOTO 10          ! skip
            ELSEIF(CKW.EQ.'MONITOR') THEN
               ITYPE(NELE) = IMO
               RLEN(NELE)  = RL
            ELSEIF(CKW.EQ.'HMONITOR') THEN
               CALL WARN('RDMDTW, no specific monitors -> general')

               ITYPE(NELE) = IMO
               RLEN(NELE)  = RL
            ELSEIF(CKW.EQ.'VMONITOR') THEN
               CALL WARN('RDMDTW, no specific monitors -> general')
               ITYPE(NELE) = IMO
               RLEN(NELE)  = RL
            ELSEIF(CKW.EQ.'SPINROTA') THEN
               ITYPE(NELE) = ISS
               RLEN(NELE)  = RL
               ISEC (NELE) = 0

               IF(AN.EQ.-1.D+0) THEN
                  ITYPE(NELE) = ISF
               ELSE
                  STRX(NELE)  = DK1
                  STRY(NELE)  = DK2
                  TWIS(NELE) = TW
               ENDIF

            ELSE
               CALL TOFF('RDMDTW, invalid keyword : '//CKW)
            ENDIF

            LELE = NELE

            IF(ISEC(NELE)*IGLBSL.GT.1) THEN
               D          = 1.0D+0/DBLE(IGLBSL)
               STRX(NELE) = STRX(NELE) * D
               STRY(NELE) = STRY(NELE) * D
               RLEN(NELE) = RLEN(NELE) * D
            ENDIF

         ENDIF

            NLAT = NLAT + 1
            IF(NLAT.GT.MAXLIST) CALL TOFF(
     +           'RDMDTW, increase MAXLIST')
            LATELE(NLAT) = LELE
            
            IF(ISEC(LELE)*IGLBSL.GT.1) THEN
               LSLIC(NLAT) =  ISEC(LELE)*IGLBSL
               RLHLP = RLEN(LELE) * DBLE(LSLIC(NLAT))
            ELSE
               LSLIC(NLAT) = 1
               RLHLP = RLEN(LELE)
            ENDIF
            
            NTOTSL = NTOTSL + MAX0(1,LSLIC(NLAT))
            
            REFMD        = 1.D+0 ! entry:0, center:.5, exit:1
            PRETST       = AT - EPOS - REFMD*RLHLP
            
            IF(PRETST.LT.-1.D-10) THEN
               WRITE(CDUMPR,'(D14.6)') PRETST
               WRITE(CDUMAT,'(D14.6)') AT
               CALL WARN('RDMDTW, negative predrift'//
     +              CDUMPR//' before '//
     +              CNM//' at '// CDUMAT)
            ENDIF
            
            IF(DABS(PRETST).GT.1.D-9) THEN
               PREDRF(NLAT) = PRETST
            ELSE
               PREDRF(NLAT) = 0.D+0
            ENDIF
            
            EPOS = EPOS + PREDRF(NLAT) + RLHLP
            IF(ITYPE(NELE).EQ.IRF) NCAV = NCAV + NSUPER

 10   CONTINUE

      IF(CDATA.EQ.'SURVEY') THEN
         READ(15,'(3E16.9/3E16.9)',END=99)
     +        XMAD,YMAD,ZMAD,
     +        RMIN,RMAX,CIRC
         
         WRITE(16,'(A)')       ' TRAILER RECORD :'
         WRITE(16,'(A,3D16.8)')' MACHINE CENTRE   : ',XMAD,YMAD,ZMAD
         WRITE(16,'(A,3F20.6)')' Rmin, Rmax, CIRC.: ',RMIN,RMAX,CIRC
         WRITE(16,'(A/)')      ' END OF MAD/SURVEY INPUT'
         
      ELSEIF(CDATA.EQ.'TWISS') THEN
         READ(15,'(3E16.9/5E16.9/5E16.9))',END=99)
     +        DELTAP,GAMTR,CIRC,
     +        COSMUX,QX,QXP,BXMAX,DXMAX,
     +        COSMUY,QY,QYP,BYMAX,DYMAX
         
         WRITE(16,'(A)')       ' TRAILER RECORD :'
         WRITE(16,'(A,3D16.8)')' DELTAP , GAMMAtr, CIRC : '
     +        ,DELTAP,GAMTR,CIRC
         WRITE(16,'(A)') ' COS(MU), Q, Qprime, BETAmax, DISPmax ::'
         WRITE(16,'(A,5D16.8)')'   X : ',COSMUX,QX,QXP,BXMAX,DXMAX
         WRITE(16,'(A,5D16.8)')'   Y : ',COSMUY,QY,QYP,BYMAX,DYMAX         
         WRITE(16,'(A/)')      ' END OF MAD/TWISS INPUT'
      ELSEIF(CDATA.EQ.'CHROM') THEN
         WRITE(16,'(A/)') ' END OF MAD/CHROM INPUT'
      ENDIF

      IF(BASYM) THEN
         IF(2*NLAT.GT.MAXLIST) CALL TOFF('RDMDTW, increase MAXLIST')
         LATELE(NLAT+1) = LATELE(NLAT)
         PREDRF(NLAT+1) = 0.D+0

         DO 20 I=2,NLAT
            IP         = NLAT + I
            IM         = NLAT + 1 - I
            LATELE(IP) = LATELE(IM)
            PREDRF(IP) = PREDRF(IM + 1)
 20      CONTINUE

         NLAT = IP

         IF(PREDRF(1).NE.0.D+0) THEN
            NLAT = NLAT + 1
            IF(NLAT.GT.MAXLIST) CALL TOFF('RDMDTW, increase MAXLIST')
            NELE = NELE + 1 
            IF(NELE.GT.MAXELE) CALL TOFF('RDMDTW, increase MAXELE')
            CNAME(NELE)  = 'X______EOL'
            ITYPE(NELE)  = IMR
            ISEC (NELE)  = 0
            LATELE(NLAT) = NELE
            PREDRF(NLAT) = PREDRF(1)
         ENDIF

         EPOS = 2.D+0 * EPOS
      ENDIF

      IF(ITYPE(LATELE(NLAT)).NE.IMR) THEN
         NELE        = NELE + 1
         IF(NELE.GT.MAXELE) CALL TOFF('RDMDTW, increase MAXELE')         
         CNAME(NELE) = 'X______EOL'
         ITYPE(NELE) = IMR
         STRX(NELE)  = 0.D+0
         STRY(NELE)  = 0.D+0
         RLEN(NELE)  = 0.D+0
         TWIS(NELE)  = 0.D+0
         ISEC(NELE)  = 0
         IFLAG(NELE) = 0

         NLAT         = NLAT + 1
         IF(NELE.GT.MAXELE) CALL TOFF('RDMDTW, increase MAXLIST')
         LATELE(NLAT) = NELE 
         PREDRF(NLAT) = 0.D+0
      ENDIF

      IF(NCAV.GT.1) THEN
         WRITE(CSDUM,'(I4)') NCAV
         CALL WARN('RDMDTW, NCAV = '
     +        //CSDUM//' might be time consuming at ramp')
      ENDIF

      DRINGL = DBLE(NSUPER) * EPOS
      WRITE(16,'(2(A,I8/),A,F16.8)')
     +     ' NUMBER OF DIFFERENT ELEMENTS : ',NELE,
     +     ' TOT. NUMBER OF ELEM. IN LAT. : ',NLAT,
     +     ' TOTAL LENGTH :         '        ,DRINGL

      RETURN

 99   CALL TOFF('RDMDTW, EOF read before NREC positions !')

      E N D



      S U B R O U T I N E RDMDTW
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  RDMDTW // c h a n g e d : CKW*4 CNM*16
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     reads the mad output (twiss/survey/chrom files) as lattice input
C     ignores all the twiss parameters at the moment
C     works fine with unix & twiss/.../-files created under unix
C     trips on 'VMS' with files generated under unix and ftp-ed to the 'VAX'
C     maybe that's due to the fact that the first char of a line is used
C     quite normally under unix but might be either ignored or interpreted
C     as control char on real computers ?????
C     gotta fix that myself, i think ! don't care too much about that ! ;-)
C     who volunteers to implement a half--way consistent MAD-TWISS
C     implementation of srm-cell ?????????????  (me ???)

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MAXELE=50000, MAELNM=10, MXSTRN=4)

C=====elements / input !!!
      CHARACTER CNAME*(MAELNM)
      DIMENSION ITYPE(MAXELE),CNAME(MAXELE),RLEN(MAXELE),
     +     STRX(MAXELE),STRY(MAXELE),STRZ(MAXELE),TWIS(MAXELE),
     +     STRN(MXSTRN,MAXELE),
     +     ISTR0(MAXELE),ISTR1(MAXELE),
     +     ISEC(MAXELE),IFLAG(MAXELE),ITPEFF(MAXELE),IPS(MAXELE)
      COMMON /ELEIN/ STRX,STRY,RLEN,TWIS,STRZ,STRN,
     +     NELE,ITYPE,ISEC,IFLAG,ITPEFF,IPS,ISTR0,ISTR1,CNAME

C=====elements / internal !!!
      DIMENSION TMAP(6,6,MAXELE), VGAM(4,7,MAXELE)
      COMMON /ELEMAP/ TMAP, VGAM

C=====element-types / internal !!!
      PARAMETER(IDL=1,IBX=2,IMQ=3,ISQ=4,IRF=5,IKX=6,IKY=7,IMS=8,IBY=9,
     +     ISO=10,IRQ=11,IRX=12,IRY=13,IRO=14,ICX=15,ICY=16,ISS=17,
     +     ISF=18,IMR=19,ISM=20,IMO=21,ICH=22,ICV=23,IAS=24,IMA=25,
     +     ISR=26,
     +     ITQ=97,I00=98)

C=====lattice
      PARAMETER(MAXLIST=50000)
      DIMENSION LATELE(MAXLIST),PREDRF(MAXLIST),LSLIC(MAXLIST)
      COMMON /LATTICE/ PREDRF, DRINGL, NLAT, NSUPER,
     +     LATELE, LSLIC

C=====job parameters
      PARAMETER(PI=3.1415926535897932D+0,DEG2RD=.1745329251994329555D-1,
     +     TINY=1.D-99)
      DIMENSION IFLIP(3)
      COMMON /GLBPAR/ GEVP0, GEVM0, ELECO, ASPIN, GA0,
     +     NTRK, IOUTMD, NCAV, IGLBSL, NTOTSL, NRMPSV, IUAUTO,IFLIP,
     +     BTHICK, BSPIN

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====filenames
      CHARACTER CJOB*(MAXNME), CJOBV*(MAXNME),
     +     CLATI*(MAXNME), CLATO*(MAXNME),CLATOV*(MAXNME)
      CHARACTER*(32) COPT(10)  
      COMMON /GLBNAM/ CJOB, CJOBV, CLATI, CLATO, CLATOV, COPT

C=====local stuff
      PARAMETER(MAXINP=1000000,MLIN=80,DSKPRC=1.D-3)
      CHARACTER CLINE*(MLIN), CDUMPR*14, CDUMAT*14
      CHARACTER*8 CPRGV,CDATA,CDATM,CTIMM,CJBNM
      CHARACTER*4 CKW
      CHARACTER*16 CNM
      CHARACTER*4 CSDUM

      DSKEW(I,X) = DABS(PI/DBLE(I)-DABS(X))

C%%      write(*,*)'rdmdtw'

      NELE   = 0
      NLAT   = 0
      EPOS   = 0.D+0
      NCAV   = 0
      NTOTSL = 0

      WRITE(16,'(A/)')' READ MAD TWISS- OR SURVEY- TABLE :'
      READ(15,'(5A8,I8,L8,I8)')CPRGV,CDATA,CDATM,CTIMM,CJBNM,
     +     NSUPER,BASYM,NREC
      
      I = ICSBEG(CDATA)
      CDATA(1:) = CDATA(I:)//'                   '
      

      WRITE(16,'(5(A,A8/),A,I8/A,L8/A,I8/)')
     +     ' MAD VERSION      : ',CPRGV,
     +     ' DATA TYPE        : ',CDATA,
     +     ' CREATION DATE    : ',CDATM,
     +     ' CREATION TIME    : ',CTIMM,
     +     ' JOB NAME         : ',CJBNM,
     +     ' SUPERPERIOD      : ',NSUPER,
     +     ' ASYM. FLAG       : ',BASYM,
     +     ' NUM OF POS. RECS : ',NREC

      READ(15,'(A80)')CLINE
      WRITE(16,'(A,A80)')' TITLE : ',CLINE
      READ(15,'(A4,A16,F12.6,3E16.9)')CKW,CNM,RL,AN,DK1,DK2

      IF(CNM.EQ.'INITIAL') THEN
         READ(15,'(A80)')CLINE

         IF(CDATA.EQ.'SURVEY') THEN
            READ(15,'(4E16.9/3E16.9)',END=99)
     +           XMAD,YMAD,ZMAD,EPOS,
     +           THEMAD,PHIMAD,PSIMAD
         ELSEIF(CDATA.EQ.'TWISS') THEN
            READ(15,'(5E16.9/5E16.9/5E16.9))',END=99)
     +           ALX,BEX,DMX,DIX,DPX,
     +           ALY,BEY,DMY,DIY,DPY,
     +           XMAD,AMAD,YMAD,DBMD,EPOS
         ELSEIF(CDATA.EQ.'CHROM') THEN
            READ(15,'(5E16.9/5E16.9/5E16.9))',END=99)
     +           WX,PHX,DMX,DDX,DDPX
     +           WY,PHY,DMY,DDY,DDPY
     +           XMAD,AMAD,YMAD,DBMD,EPOS     
         ENDIF

         NREC = NREC - 1
      ELSE
         BACKSPACE(15)
      ENDIF

      DO 10 I=1,NREC
         READ(15,'(A4,A16,F12.6,3E16.9/5E16.9)',END=99)
     +        CKW,CNM,RL,AN,DK1,DK2
     +        ,TW,FR,VO,DUM1,DUM2

         IF(CDATA.EQ.'SURVEY') THEN
            READ(15,'(4E16.9/3E16.9)',END=99)
     +           XMAD,YMAD,ZMAD,AT,
     +           THEMAD,PHIMAD,PSIMAD
         ELSEIF(CDATA.EQ.'TWISS') THEN
            READ(15,'(5E16.9/5E16.9/5E16.9))',END=99)
     +           ALX,BEX,DMX,DIX,DPX,
     +           ALY,BEY,DMY,DIY,DPY,
     +           XMAD,AMAD,YMAD,DBMD,AT
         ELSEIF(CDATA.EQ.'CHROM') THEN
            READ(15,'(5E16.9/5E16.9/5E16.9))',END=99)
     +           WX,PHX,DMX,DDX,DDPX,
     +           WY,PHY,DMY,DDY,DDPY,
     +           XMAD,AMAD,YMAD,DBMD,AT     
         ENDIF

         IF(CNM(1:1).EQ.'[') GOTO 10 ! ignore generated drifts

         LELE = NUMELE(CNM)

         IF(LELE.EQ.0) THEN
            NELE = NELE + 1
            IF(NELE.GT.MAXELE) CALL TOFF('RDMDTW, increase MAXELE')
            CNAME(NELE) = CNM
            ITYPE(NELE) = 0
            STRX(NELE)  = 0.D+0
            STRY(NELE)  = 0.D+0
            RLEN(NELE)  = 0.D+0
            TWIS(NELE)  = 0.D+0
            ISEC(NELE)  = 1
            IFLAG(NELE) = 0
            IPS(NELE)   = 0

            IF    (CKW.EQ.'MARK') THEN
               ITYPE(NELE) = IMR
               ISEC (NELE) = 0
            ELSEIF(CKW.EQ.'DRIF') THEN
               ITYPE(NELE) = IDL
               RLEN(NELE)  = RL
            ELSEIF(CKW.EQ.'RBEN') THEN
               ITDUM       = IRX
               RLEN(NELE)  = RL
               STRX(NELE)  = AN

               IF(DSKEW(2,TW).LT.DSKPRC) THEN
                  ITDUM = IRY
                  STRX(NELE) = -STRX(NELE)
               ELSEIF(TW.NE.0.D+0) THEN
                  CALL TOFF('RDMDTW, no rotated bends yet !')
               ENDIF

               IF(DK1.NE.0.D+0) THEN
                  STRY(NELE) = DK1 * RL
                  IF(ITDUM.EQ.IRY) THEN
                     ITDUM = ICV
                  ELSE
                     ITDUM = ICH
                  ENDIF
               ENDIF

               ITYPE(NELE) = ITDUM                     
            ELSEIF(CKW.EQ.'SBEN') THEN
               ITDUM       = IBX
               RLEN(NELE)  = RL
               STRX(NELE)  = AN

               IF(DSKEW(2,TW).LT.DSKPRC) THEN
                  ITDUM = IBY
                  STRX(NELE) = -STRX(NELE)
               ELSEIF(TW.NE.0.D+0) THEN
                  CALL TOFF('RDMDTW, no rotated bends yet !')
               ENDIF

               IF(DK1.NE.0.D+0) THEN
                  STRY(NELE) = DK1 * RL
                  IF(ITDUM.EQ.IBY) THEN
                     ITDUM = ICY
                  ELSE
                     ITDUM = ICX
                  ENDIF
               ENDIF

               ITYPE(NELE) = ITDUM    
            ELSEIF(CKW.EQ.'QUAD') THEN
               ITDUM       = IMQ
               RLEN(NELE)  = RL
               STRX(NELE)  = DK1 * RL

               IF(DSKEW(4,TW).LT.DSKPRC) THEN
                  ITDUM = ISQ
                  IF(TW.LT.0.D+0) STRX(NELE) = -STRX(NELE)
C     right convention?
               ELSEIF(TW.NE.0.D+0) THEN
                  ITDUM = IRQ
                  TWIS(NELE) = TW
               ENDIF

               ITYPE(NELE) = ITDUM  
            ELSEIF(CKW.EQ.'SEXT') THEN

               IF(TW.NE.0.D+0) THEN
                  CALL WARN(
     +                 'RDMDTW, no rotated sextupoles yet -> drift')
                  NELE = NELE - 1
                  GOTO 10       ! skip
               ENDIF

               ITYPE(NELE) = IMS
               RLEN(NELE)  = RL
               STRX(NELE)  = DK2 * RL
            ELSEIF(CKW.EQ.'OCTU') THEN
               CALL WARN('RDMDTW, no octupoles yet -> drift')
               NELE = NELE - 1
               GOTO 10          ! skip
            ELSEIF(CKW.EQ.'MULT') THEN

               IF(DK1.NE.0.D+0) THEN
                  ITYPE(NELE) = ITQ
                  STRX(NELE)  = DK1
                  IF(AN .NE. 0.D+0) CALL TOFF(
     +                 'RMDMTW, no thin lense quad + dipole yet')
                  IF(VO .NE. 0.D+0) CALL TOFF(
     +                 'RMDMTW, no rotated thin lense quad yet')
                  IF(DABS(DK2)+DABS(FR) .NE. 0.D+0) CALL WARN(
     +                 'RMDMTW, no thin lense quad + h.o.multip. yet')
               ELSE
                  IF(AN .NE. 0.D+0) CALL TOFF(
     +                 'RMDMTW, no multipoles containing dipols yet') 
                  CALL WARN('RDMDTW, no multipoles yet -> drift')
                  NELE = NELE - 1
                  GOTO 10       ! skip
               ENDIF

            ELSEIF(CKW.EQ.'SOLE') THEN
               ITYPE(NELE) = ISO
               RLEN(NELE)  = RL
               STRX(NELE)  = FR
               ISEC (NELE) = 0
            ELSEIF(CKW.EQ.'RFCA') THEN
               ITYPE(NELE) = IRF
               RLEN(NELE)  = RL
               STRX(NELE)  = FR
               STRY(NELE)  = VO
               ISEC (NELE) = 0
            ELSEIF(CKW.EQ.'ELSE') THEN
               CALL WARN('RDMDTW, no el. sep. yet -> drift')
               NELE = NELE - 1
               GOTO 10          ! skip
            ELSEIF(CKW.EQ.'HKIC') THEN

               IF(TW.NE.0.D+0) THEN
                  CALL WARN(
     +                 'RDMDTW, no rotated kickers yet -> drift')
                  NELE = NELE - 1
                  GOTO 10       ! skip
               ENDIF

               ITYPE(NELE) = IKX
               RLEN(NELE)  = RL
               STRX(NELE)  = FR
            ELSEIF(CKW.EQ.'VKIC') THEN

               IF(TW.NE.0.D+0) THEN
                  CALL WARN(
     +                 'RDMDTW, no rotated kickers yet -> drift')
                  NELE = NELE - 1
                  GOTO 10       ! skip
               ENDIF

               ITYPE(NELE) = IKY
               RLEN(NELE)  = RL
               STRX(NELE)  = FR
            ELSEIF(CKW.EQ.'RCOL') THEN
               CALL WARN('RDMDTW, collimators ignored -> drift')
               NELE = NELE - 1
               GOTO 10          ! skip
            ELSEIF(CKW.EQ.'SROT') THEN
               CALL WARN('RDMDTW, no long. rotation yet -> drift')
               NELE = NELE - 1
               GOTO 10          ! skip
            ELSEIF(CKW.EQ.'YROT') THEN
               CALL WARN('RDMDTW, no vert. rotation yet -> drift')
               NELE = NELE - 1
               GOTO 10          ! skip
            ELSEIF(CKW.EQ.'MONI') THEN
               ITYPE(NELE) = IMO
               RLEN(NELE)  = RL
            ELSEIF(CKW.EQ.'HMON') THEN
               CALL WARN('RDMDTW, no specific monitors -> general')
               ITYPE(NELE) = IMO
               RLEN(NELE)  = RL
            ELSEIF(CKW.EQ.'VMON') THEN
               CALL WARN('RDMDTW, no specific monitors -> general')
               ITYPE(NELE) = IMO
               RLEN(NELE)  = RL
            ELSEIF(CKW.EQ.'SPIN') THEN
               ITYPE(NELE) = ISS
               RLEN(NELE)  = RL
               ISEC (NELE) = 0

               IF(AN.EQ.-1.D+0) THEN
                  ITYPE(NELE) = ISF
               ELSE
                  STRX(NELE)  = DK1
                  STRY(NELE)  = DK2
                  TWIS(NELE) = TW
               ENDIF

            ELSE
               CALL TOFF('RDMDTW, invalid keyword : '//CKW)
            ENDIF

            LELE = NELE

            IF(ISEC(NELE)*IGLBSL.GT.1) THEN
               D          = 1.0D+0/DBLE(IGLBSL)
               STRX(NELE) = STRX(NELE) * D
               STRY(NELE) = STRY(NELE) * D
               RLEN(NELE) = RLEN(NELE) * D
            ENDIF

         ENDIF

            NLAT = NLAT + 1
            IF(NLAT.GT.MAXLIST) CALL TOFF(
     +           'RDMDTW, increase MAXLIST')
            LATELE(NLAT) = LELE
            
            IF(ISEC(LELE)*IGLBSL.GT.1) THEN
               LSLIC(NLAT) =  ISEC(LELE)*IGLBSL
               RLHLP = RLEN(LELE) * DBLE(LSLIC(NLAT))
            ELSE
               LSLIC(NLAT) = 1
               RLHLP = RLEN(LELE)
            ENDIF
            
            NTOTSL = NTOTSL + MAX0(1,LSLIC(NLAT))
            
            REFMD        = 1.D+0 ! entry:0, center:.5, exit:1
            PRETST       = AT - EPOS - REFMD*RLHLP
            
            IF(PRETST.LT.-1.D-10) THEN
               WRITE(CDUMPR,'(D14.6)') PRETST
               WRITE(CDUMAT,'(D14.6)') AT
               CALL WARN('RDMDTW, negative predrift'//
     +              CDUMPR//' before '//
     +              CNM//' at '// CDUMAT)
            ENDIF
            
            IF(DABS(PRETST).GT.1.D-9) THEN
               PREDRF(NLAT) = PRETST
            ELSE
               PREDRF(NLAT) = 0.D+0
            ENDIF
            
            EPOS = EPOS + PREDRF(NLAT) + RLHLP
            IF(ITYPE(NELE).EQ.IRF) NCAV = NCAV + NSUPER

 10   CONTINUE

      IF(CDATA.EQ.'SURVEY') THEN
         READ(15,'(3E16.9/3E16.9)',END=99)
     +        XMAD,YMAD,ZMAD,
     +        RMIN,RMAX,CIRC
         
         WRITE(16,'(A)')       ' TRAILER RECORD :'
         WRITE(16,'(A,3D16.8)')' MACHINE CENTRE   : ',XMAD,YMAD,ZMAD
         WRITE(16,'(A,3F20.6)')' Rmin, Rmax, CIRC.: ',RMIN,RMAX,CIRC
         WRITE(16,'(A/)')      ' END OF MAD/SURVEY INPUT'
         
      ELSEIF(CDATA.EQ.'TWISS') THEN
         READ(15,'(3E16.9/5E16.9/5E16.9))',END=99)
     +        DELTAP,GAMTR,CIRC,
     +        COSMUX,QX,QXP,BXMAX,DXMAX,
     +        COSMUY,QY,QYP,BYMAX,DYMAX
         
         WRITE(16,'(A)')       ' TRAILER RECORD :'
         WRITE(16,'(A,3D16.8)')' DELTAP , GAMMAtr, CIRC : '
     +        ,DELTAP,GAMTR,CIRC
         WRITE(16,'(A)') ' COS(MU), Q, Qprime, BETAmax, DISPmax ::'
         WRITE(16,'(A,5D16.8)')'   X : ',COSMUX,QX,QXP,BXMAX,DXMAX
         WRITE(16,'(A,5D16.8)')'   Y : ',COSMUY,QY,QYP,BYMAX,DYMAX         
         WRITE(16,'(A/)')      ' END OF MAD/TWISS INPUT'
      ELSEIF(CDATA.EQ.'CHROM') THEN
         WRITE(16,'(A/)') ' END OF MAD/CHROM INPUT'
      ENDIF

      IF(BASYM) THEN
         IF(2*NLAT.GT.MAXLIST) CALL TOFF('RDMDTW, increase MAXLIST')
         LATELE(NLAT+1) = LATELE(NLAT)
         PREDRF(NLAT+1) = 0.D+0

         DO 20 I=2,NLAT
            IP         = NLAT + I
            IM         = NLAT + 1 - I
            LATELE(IP) = LATELE(IM)
            PREDRF(IP) = PREDRF(IM + 1)
 20      CONTINUE

         NLAT = IP

         IF(PREDRF(1).NE.0.D+0) THEN
            NLAT = NLAT + 1
            IF(NLAT.GT.MAXLIST) CALL TOFF('RDMDTW, increase MAXLIST')
            NELE = NELE + 1 
            IF(NELE.GT.MAXELE) CALL TOFF('RDMDTW, increase MAXELE')
            CNAME(NELE)  = 'X______EOL'
            ITYPE(NELE)  = IMR
            ISEC (NELE)  = 0
            LATELE(NLAT) = NELE
            PREDRF(NLAT) = PREDRF(1)
         ENDIF

         EPOS = 2.D+0 * EPOS
      ENDIF

      IF(ITYPE(LATELE(NLAT)).NE.IMR) THEN
         NELE        = NELE + 1
         IF(NELE.GT.MAXELE) CALL TOFF('RDMDTW, increase MAXELE')         
         CNAME(NELE) = 'X______EOL'
         ITYPE(NELE) = IMR
         STRX(NELE)  = 0.D+0
         STRY(NELE)  = 0.D+0
         RLEN(NELE)  = 0.D+0
         TWIS(NELE)  = 0.D+0
         ISEC(NELE)  = 0
         IFLAG(NELE) = 0

         NLAT         = NLAT + 1
         IF(NELE.GT.MAXELE) CALL TOFF('RDMDTW, increase MAXLIST')
         LATELE(NLAT) = NELE 
         PREDRF(NLAT) = 0.D+0
      ELSE
         CNAME(NELE) = 'X______EOL'
      ENDIF

      IF(NCAV.GT.1) THEN
         WRITE(CSDUM,'(I4)') NCAV
         CALL WARN('RDMDTW, NCAV = '
     +        //CSDUM//' might be time consuming at ramp')
      ENDIF

      DRINGL = DBLE(NSUPER) * EPOS
      WRITE(16,'(2(A,I8/),A,F16.8)')
     +     ' NUMBER OF DIFFERENT ELEMENTS : ',NELE,
     +     ' TOT. NUMBER OF ELEM. IN LAT. : ',NLAT,
     +     ' TOTAL LENGTH :         '        ,DRINGL

      RETURN

 99   CALL TOFF('RDMDTW, EOF read before NREC positions !')

      E N D


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C30                                                                           C
C@ PROGRAM=SPRINT, MODULE=STRING, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      S U B R O U T I N E STRIP0(CLINE,IE,BCONT,*,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  STRIP0
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     strips a string from all its blanks (-> mad !!)
C     at the moment also strips quoted strings ==>>
C     change but care for continuation !!!!!!!!!!!!!
C     return   : non-blank line (default)
C     return 1 : blank line 
C     return 2 : blank line in between continuation
C     !!! whole thing still in developement !!!


      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====local stuff
      PARAMETER(MLIN=132)
      CHARACTER CLINE*(MLIN)

C%%      write(*,*)'strip0'

      IE = 0

      DO 10 I=1,MLIN
         IF(CLINE(I:I).EQ.' ') GOTO 10
         IE = IE +1
         CLINE(IE:IE) = CLINE(I:I)
 10   CONTINUE

      IF(IE.EQ.0) THEN
         IF(BCONT) RETURN
         RETURN 1
      ELSE
         CLINE(IE+1:MLIN) = ' '
      ENDIF

      RETURN
      E N D

      S U B R O U T I N E STRIPX(CLINE,IE,BCONT,*,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  STRIPX
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     strips a string from its mad0-extension markers: '!@'
C     return   : command line (default)
C     return 1 : comment line 
C     return 2 : comment line in between continuation
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====local stuff
      PARAMETER(MLIN=132)
      CHARACTER CLINE*(MLIN)

C%%      write(*,*)'stripx'

 10   CONTINUE
      IX = INDEX(CLINE(1:IE),'!@')
      IF(IX.NE.0) THEN
         CLINE(IX:IE) = CLINE(IX+2:IE)
         IE = IE - 2
         GOTO 10
      ENDIF

      IF(IE.EQ.0) THEN
         IF(BCONT) RETURN
         RETURN 1
      ELSE
         CLINE(IE+1:MLIN) = ' '
      ENDIF

      RETURN
      E N D

      S U B R O U T I N E STRIPC(CLINE,IE,BCONT,*,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  STRIPC
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     strips a string from its mad/FORTRAN comments: '! something'
C     return   : command line (default)
C     return 1 : comment line 
C     return 2 : comment line in between continuation
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====local stuff
      PARAMETER(MLIN=132)
      CHARACTER CLINE*(MLIN)

C%%      write(*,*)'stripc'

      IC = INDEX(CLINE(1:IE),'!')

      IF(IC.GT.0) THEN
         IE = IC - 1
      ENDIF

      IF(IE.EQ.0) THEN
         IF(BCONT) RETURN
         RETURN 1
      ELSE
         CLINE(IE+1:MLIN) = ' '
      ENDIF

      RETURN
      E N D

      INTEGER F U N C T I O N IMADWD(CBUF,IEB)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  IMADWD
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     counts the number of tokens in mad-0 input line

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====local stuff
      PARAMETER(MBUF=400)
      CHARACTER CDELIM*2
      CHARACTER CBUF*(MBUF)

      DATA CDELIM /':,'/

C%%      write(*,*)'imadwd'

      IMADWD = 1

      DO 10 I=1,IEB
         IF(INDEX(CDELIM,CBUF(I:I)).NE.0)  IMADWD = IMADWD + 1
 10   CONTINUE

      RETURN
      E N D
  
      S U B R O U T I N E MADWRD(CKW,IEK,CBUF,IEB,IW)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  MADWRD
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     extracts a token from a strip0-ed string
C     an assignment is (at the moment// for the zeroth order mad-0 reader !!)
C     considered one token
C     !!! whole thing still in developement !!!
 
      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)
 
C=====local stuff
      PARAMETER(MBUF=400,MKEYW=32)
      CHARACTER CDELIM*2
      CHARACTER CBUF*(MBUF), CKW*(MKEYW)

      DATA CDELIM /':,'/

C%%      write(*,*)'madwrd'
      CKW = ' '
      IEK = 1
      IL  = 1
      IR  = 0
      JW  = 0
 
      DO 10 I=1,IEB
         IND = INDEX(CDELIM,CBUF(I:I))
         IF(IND.NE.0 .OR. I.EQ.IEB) THEN
            JW = JW + 1
            IL = IR + 1
            IR = I
            IF(JW.EQ.IW) GOTO 11
         ENDIF
 10   CONTINUE
      RETURN
 11   CONTINUE
 
      IF(IND.EQ.2) IR = IR-1
      IEK = IR - IL + 1
      IF(IEK.GE.MKEYW) CALL TOFF('MADWRD, keyword too long')
      CKW = CBUF(IL:IR)

C      write(*,*) JW,' >',CKW(:IEK),'<'
 
      RETURN
      E N D

      S U B R O U T I N E PRSASS(CW,IE1,CN,CV,BASS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  PRSASS
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     parses an assigment 'CW' of the form name=value
C     to 'CN' <- name ; 'CV' <- value
C     !!! whole thing still in developement !!!

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER*(*) CW,CN,CV

C%%      write(*,*)'prsass'

      CN   = ' '
      CV   = ' '
      BASS = .TRUE.
      IEQ = INDEX(CW,'=')

      IF(IEQ.EQ.0) THEN
         BASS = .FALSE.
         RETURN
      ENDIF

      CN = CW(1    : IEQ-1)
      CV = CW(IEQ+1: IE1  )

      RETURN
      E N D

      S U B R O U T I N E STRIP(CIN,COUT,LOUT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  STRIP
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Reduces multiple blanks to single blanks except in 'strings'
C     and takes out leading blanks
C     Parameter list:
C     CIN : input character string possibley containing multiple blanks
C     COUT: output character string equals CIN except all not 'quoted'
C           multiple blanks are reduced to single blanks and leading blanks
C           are eliminated.
C     LOUT: length of COUT without trailing blanks
C     Internal parameters:
C     BUNQUO: .TRUE.  if the next character in CLINE is not in a 'string' and
C             .FALSE. if the next character in CLINE is in a 'string'.
C
C     This routine will be called several times on individual cmd-strings,
C     so internal vars must not be initialized by 'DATA'
C     and not kept by 'SAVE'

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CIN*(*), COUT*(*), C*1

C%%      write(*,*)'strip'

      BUNQUO = .TRUE.
      LIN    = LEN(CIN)
      NBL    = 1
      LOUT   = 1

      DO 10 I=1,LIN
         C = CIN(I:I)
         IF(C.EQ.'''') THEN
            BUNQUO          = .NOT.BUNQUO
            COUT(LOUT:LOUT) = C
            LOUT            = LOUT+1
         ELSEIF(C.EQ.' ') THEN

            IF(.NOT.BUNQUO .OR. NBL.EQ.0) THEN
               COUT(LOUT:LOUT) = C
               LOUT            = LOUT+1
               NBL             = NBL+1
            ENDIF

         ELSE
            COUT(LOUT:LOUT) = C
            LOUT            = LOUT+1
            NBL             = 0
         ENDIF

 10   CONTINUE

      COUT(LOUT:) = ' '
      LOUT        = LOUT-1

      RETURN
      E N D

      LOGICAL F U N C T I O N BCS(CL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  BCS
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Does CL describe a character string -> .true., .false.

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====local stuff
      CHARACTER CL*(*)

      IE =ICSEND(CL)
      BCS = CL(1:1).EQ.''''.AND.CL(IE:IE).EQ.''''

      RETURN
      E N D

      INTEGER F U N C T I O N IBFCS(CSTR)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  IBFCS
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     converts a char*31 (maybe later *32) stream of 'T' and 'F'
C     into bitarray
C     right most symbol revers to bit 0
C     if CSTR is too short it will be alligned to the right

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(MBITS=31)
C     PARAMETER(MBITS=32) ! once we trust signed integers !!!!

      CHARACTER CSTR*(*), CS*(MBITS)

C%%      write(*,*)'ibfcs'

      IB = ICSBEG(CSTR)
      IE = ICSEND(CSTR)
      LS = IE - IB + 1
      IF(LS.GT.MBITS) CALL TOFF(
     +     'IBFCS, bolean stream too large')
      IBFCS = 0
      IF(IE.EQ.0) RETURN
      CS = CSTR(IB:IE)

      DO 10 I=LS,1,-1
         IF(CS(I:I).EQ.'T') THEN
            IBFCS = IBSET(IBFCS,LS-I)
         ELSEIF(CS(I:I).NE.'F') THEN
            CALL TOFF('IBFCS, invallid symbol: '//CS(I:I))
         ENDIF
 10   CONTINUE

      RETURN
      E N D

      CHARACTER*31 F U N C T I O N CSFIB(IFL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  CSFIB
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     returns a stream of 'T' and 'F' referring to the bits in IFL
C     returns ' ' if IFL is 0 
C     rightmost symbol refers to bit 0
C     maybe switch to CHARACTER*32 once we trust signed integers !

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====bitarray codes for 'IFLAG' [be carefull, stay signed_positive !]
      PARAMETER(MBITS=31)
      PARAMETER(JBNOSP=0,JBSNEG=1,JBSPC7=28,JBSPC8=29,JBSPC9=30)

C%%      write(*,*)'csfib'

      CSFIB = ' '
      IF(IFL.EQ.0) RETURN

      DO 10 I=0,MBITS-1
         J = MBITS - I
         IF(BTEST(IFL,I)) THEN
            CSFIB(J:J) = 'T'
         ELSE
            CSFIB(J:J) = 'F'
         ENDIF
 10   CONTINUE

      E N D

      CHARACTER*(*) F U N C T I O N CFINFL(INUM,CBLANK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  CFINFL
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Charcter_From_INteger_FiLled : returns LEN(CFINFL) digits
C     of INUM, leading blanks patched with CBLANK
C     EXAMPLE:  char*5 ctest ; ctest=cfinfl(123,'x') -> ctest='xx123'

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CNUM*10, CBLANK*1

      DATA CNUM /'0123456789'/

C%%      write(*,*)'cfinfl'

      MP = LEN(CFINFL)
      NP = MP 
      IH = INUM

      DO 10 I=MP,1,-1
         IM = MOD(IH,10) + 1
         CFINFL(I:I) = CNUM(IM:IM)
         NP = NP - 1
         IH = IH/10         
         IF(IH.EQ.0) GOTO 11
         IF(NP.LT.1) CALL TOFF('CFINFL, integer too long')
 10   CONTINUE
 11   CONTINUE

      DO 20 I=NP,1,-1
         CFINFL(I:I) = CBLANK         
 20   CONTINUE

      E N D

      REAL*8 F U N C T I O N REFCS(CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  REFCS
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     turns a string into a real number

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(LDEC=308)

      CHARACTER CLINE*(*),NUM*10,DEDE*4, CDUM*32
      DIMENSION CDEC(0:LDEC)
      SAVE ICALL,CDEC
      DATA NUM / '1234567890' /
      DATA DEDE / 'DEde' /
      DATA ICALL / 0 /

C%%      write(*,*)'refcs'

      IF(ICALL.EQ.0) THEN
         ICALL = 1
         CDEC(0) = 1.D+0
         DO 5 I=1,LDEC
            CDEC(I) = CDEC(I-1)*10.D+0
 5       CONTINUE
      ENDIF

      LENGTH = LEN(CLINE)
C=====length of the incomming character
      REFCS  = 0.D+0
C=====final real number
      IFRA = 0
C=====counts numbers after '.'
      LPOI = 0
C=====becomes 1 with '.'
      IEXP = 0
C=====number of the exponent
      LEXP = 0
C=====becomes 1 for the exponent
      IPME = 1
C=====sign of exponent
      IPMM = 1
C=====sign of number
      IPOS = 0
C=====non blank position in string
      IGAP = 0
C=====counts gaps

      DO 10 I=1,LENGTH
         IF(CLINE(I:I).EQ.' ') THEN
            IF(IPOS.GT.0) IGAP = 1
            GOTO 10
         ENDIF
         IPOS = IPOS+1
         IF(IGAP.GT.0) GOTO 500
         IF(INDEX(NUM,CLINE(I:I)).NE.0) THEN
            IF(LEXP.NE.0) THEN
               IEXP = IEXP*10 + INDEX(NUM(1:9),CLINE(I:I))
            ELSEIF(LPOI.EQ.0) THEN
               REFCS  = REFCS*10.D+0 + INDEX(NUM(1:9),CLINE(I:I))
            ELSEIF(LPOI.EQ.1) THEN
               IFRA = IFRA + 1
               IF(IFRA.GT.LDEC) CALL TOFF('REFCS, increase LDEC')
               REFCS = REFCS + INDEX(NUM(1:9),CLINE(I:I))/CDEC(IFRA)
            ELSE
               CALL TOFF('@@@ internal ERROR in REFCS')
            ENDIF
         ELSEIF(CLINE(I:I).EQ.'.') THEN
            IF(LPOI.NE.0.OR.LEXP.NE.0) GOTO 500
            LPOI = 1
         ELSEIF(INDEX(DEDE,CLINE(I:I)).NE.0) THEN
            IF(LEXP.NE.0) GOTO 500
            LEXP = 1
         ELSEIF(INDEX('+-',CLINE(I:I)).NE.0) THEN
            IF(IPOS.EQ.1) THEN
               IF(CLINE(I:I).EQ.'-') IPMM = -1
            ELSEIF(INDEX(DEDE,CLINE(I-1:I-1)).NE.0) THEN
               IF(CLINE(I:I).EQ.'-') IPME = -1
            ELSE
               GOTO 500
            ENDIF
         ELSE
            GOTO 500
         ENDIF
 10   CONTINUE

      REFCS = IPMM*REFCS
      IF(LEXP.NE.0) THEN
         IF(IEXP.GT.LDEC) GOTO 500
         IF(IPME.EQ.+1) THEN
            REFCS = REFCS*CDEC(IEXP)
         ELSEIF(IPME.EQ.-1) THEN
            REFCS = REFCS/CDEC(IEXP)
         ENDIF
      ENDIF

      RETURN

 500  CONTINUE
      CDUM = CLINE
      CALL TOFF('REFCS, string is not a real number: '//CDUM)

      RETURN
      E N D

      CHARACTER*(*) F U N C T I O N CSFIN(INT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  CSFIN
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     turns an integer into a string

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CFMT*3, CDUM*10

C%%      write(*,*)'csfin'

      WRITE(CFMT,'(I3)') LEN(CSFIN)
C     next is a G77 workaround
      CDUM = '(I'
      I = ICSBEG(CFMT)
      CDUM(3:) = CFMT(I:)
      CDUM(7-I:) = ')'
      WRITE(CSFIN,CDUM) INT

      RETURN
      E N D

      INTEGER F U N C T I O N INFCS(CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  INFCS
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     turns a string into an integer

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CLINE*(*),NUM*10,CLDUM*32
      DATA NUM / '1234567890' /

C%%      write(*,*)'infcs'

      LENGTH = LEN(CLINE)
C=====length of the incomming character
      INFCS  = 0
C=====final integer
      IPMM = 1
C=====sign of number
      IPOS = 0
C=====non black position in string
      IGAP = 0
C=====counts gaps

      DO 10 I=1,LENGTH
         IF(CLINE(I:I).EQ.' ') THEN
            IF(IPOS.GT.0) IGAP = 1
            GOTO 10
         ENDIF
         IPOS = IPOS+1
         IF(IGAP.GT.0) GOTO 500

         IF(INDEX(NUM,CLINE(I:I)).NE.0) THEN
            INFCS = INFCS*10 + INDEX(NUM(1:9),CLINE(I:I))
         ELSEIF(INDEX('+-',CLINE(I:I)).NE.0) THEN
            IF(IPOS.EQ.1) THEN
               IF(CLINE(I:I).EQ.'-') IPMM = -1
            ELSE
               GOTO 500
            ENDIF
         ELSE
            GOTO 500
         ENDIF

 10   CONTINUE

      INFCS = IPMM*INFCS
      RETURN

 500  CONTINUE
C     ! next is a g77 workaround !
      CLDUM = CLINE(ICSBEG(CLINE):ICSEND(CLINE))
      CALL TOFF('INFCS, string '
     +     //CLDUM//' is not an integer')

      RETURN
      E N D

      CHARACTER*(*) F U N C T I O N COMPOS(CLINE,IP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  COMPOS
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     returns the 'IP'-th word in 'CLINE'
C     needs preprocessing by 'COMMENTS' and 'COMMANDS' (esp. 'STRIP'!)

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CLINE*(*), C*1

C%%      write(*,*)'compos'

      COMPOS = ' '
      IB   = ICSBEG(CLINE)
      IE   = ICSEND(CLINE)
      IF(IE-IB .EQ. 0) RETURN
      IC   = 1
      J    = 1
      BU   = .TRUE.

      DO 10 I=IB,IE
         C=CLINE(I:I)
         IF(C.EQ.'''') BU = .NOT.BU

         IF(((BU .AND. C.NE.' ') .OR. .NOT.BU)
     +        .AND. IC.EQ.IP) THEN
            COMPOS(J:J) = C
            J         = J+1
         ENDIF

         IF(BU .AND. C.EQ.' ') IC = IC+1
         IF(IC.GT.IP) RETURN
 10   CONTINUE

      RETURN
      E N D

      CHARACTER*(*) F U N C T I O N CIFCL(CSTR,ILINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  CIFCL
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Takes the leading and trailing quotes off the 'sting' CSTR
C     converts a 'L'anguage-'C'har_str
C     to an 'I'nternal- ('FORTRAN') -'C'har_str, 'NOT' enclosed in quotes !
C     ILINE : is the command line in which CSTR appeared.
C             It is just supplied for making debugging more easy !

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====maximum word length
      PARAMETER(MAXNME=132)

C=====local stuff
      CHARACTER CSTR*(MAXNME)

C%%      write(*,*)'cifcl'

      CIFCL = ' '
      IE    = ICSEND(CSTR)

      IF(CSTR(1:1).EQ.'''' .AND. CSTR(IE:IE).EQ.'''') THEN
         CIFCL(:IE-2) = CSTR(2:IE-1)
      ELSE
         CALL TNOFF(
     +        ILINE, 'CIFCL, '//CSTR(:IE)//
     +        ' is not enclosed in quotes')
      ENDIF

      RETURN
      E N D

      CHARACTER*(*) F U N C T I O N CPOS(CLINE,IPOS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  CPOS
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     returns the 'IPOS'-th word in 'CLINE'

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CLINE*(*)

C%%      write(*,*)'cpos'

      LENGTH = LEN(CLINE)
      IF(IPOS.LT.1) CALL TOFF('CPOS, non natural position')

      IF(CLINE(1:1).EQ.' ') THEN
         ICOUNT = 0
      ELSE
         ICOUNT = 1
      ENDIF
      IBEG = 1
      IF(ICOUNT.EQ.IPOS) GOTO 20

      DO 10 I=1,LENGTH-1
         IBEG = I+1
         IF((CLINE(I:I).EQ.' ').AND.(CLINE(IBEG:IBEG).NE.' '))
     +        ICOUNT=ICOUNT+1
         IF(ICOUNT.EQ.IPOS) GOTO 20
 10   CONTINUE

      CPOS = ' '
      RETURN

 20   CONTINUE

      DO 30 I=IBEG+1,LENGTH
         IF(CLINE(I:I).EQ.' ') THEN
            IEND = I-1
            GOTO 40
         ENDIF
 30   CONTINUE

      IEND = LENGTH
 40   CONTINUE

      IF((IEND-IBEG+1).GT.LEN(CPOS))
     +     CALL TOFF('CPOS, word too long')
      CPOS = ' '
      CPOS(1:IEND-IBEG+1) = CLINE(IBEG:IEND)

      RETURN
      E N D

      CHARACTER*(*) F U N C T I O N CASEUP(CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  CASEUP
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     returnes the upcased string of 'CLINE'

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CLINE*(*),CALPH*26,CBETA*26
      DATA CALPH / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
      DATA CBETA / 'abcdefghijklmnopqrstuvwxyz' /

C%%      write(*,*)'caseup'

      CASEUP = CLINE

      DO 10 I=1,MIN(LEN(CASEUP),LEN(CLINE))
         J = INDEX(CBETA,CLINE(I:I))
         IF(J.NE.0) CASEUP(I:I) = CALPH(J:J)
 10   CONTINUE

      RETURN
      E N D

      CHARACTER*(*) F U N C T I O N CASEDOWN(CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  CASEDOWN
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     returnes the downcased string of 'CLINE'

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CLINE*(*),CALPH*26,CBETA*26

      DATA CALPH / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
      DATA CBETA / 'abcdefghijklmnopqrstuvwxyz' /

C%%      write(*,*)'casedown'

      CASEDOWN = CLINE
      DO 10 I=1,MIN(LEN(CASEDOWN),LEN(CLINE))
         J = INDEX(CALPH,CLINE(I:I))
         IF(J.NE.0) CASEDOWN(I:I) = CBETA(J:J)
 10   CONTINUE

      RETURN
      E N D

      INTEGER F U N C T I O N ICSEND(CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  ICSEND
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CLINE  (input) : Character string
C     ICSEND (output): possition of the last non blank character in CLINE
C                      ICSEND is 0 if CLINE is blank.

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CLINE*(*)

C%%      write(*,*)'icsend'

      DO 10 I=LEN(CLINE),1,-1
         IF(CLINE(I:I).NE.' ') THEN
            ICSEND = I
            GOTO 19
         ENDIF
 10   CONTINUE
      ICSEND = 0
 19   CONTINUE

      RETURN
      E N D

      INTEGER F U N C T I O N ICSBEG(CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  ICSBEG
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     returns the position of the last non blank
C     in 'CLINE'

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CLINE*(*)

C%%      write(*,*)'icsbeg'

      DO 10 I=1,LEN(CLINE)
         IF(CLINE(I:I).NE.' ') THEN
            ICSBEG = I
            GOTO 19
         ENDIF
 10   CONTINUE
      ICSBEG = 0
 19   CONTINUE

      RETURN
      E N D

      INTEGER F U N C T I O N NUMCS(CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  NUMCS
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     returns the number of words in 'CLINE'

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER CLINE*(*)

C%%      write(*,*)'numcs'

      IBEG=ICSBEG(CLINE)
      IEND=ICSEND(CLINE)

      IF(IBEG.EQ.IEND) THEN
         NUMCS  = 0
      ELSE
         BUNQUO = .TRUE.
         NUMCS=1

         DO 10 I=IBEG+1,IEND-1
            IF(CLINE(I:I).EQ.'''') BUNQUO = .NOT.BUNQUO
            IF(CLINE(I-1:I-1).NE.' '.AND.CLINE(I:I).EQ.' '
     +           .AND.BUNQUO)
     +           NUMCS=NUMCS+1
 10      CONTINUE
      ENDIF

      RETURN
      E N D


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C33                                                                           C
C@ PROGRAM=SPRINT, MODULE=TOOLS, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      S U B R O U T I N E TOFF(CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  TOFF
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     terminates with the message 'CLINE'

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER*(*) CLINE

C%%      write(*,*)'toff'

      WRITE( *,'(A,A)') ' *** ERROR in ',CLINE
      WRITE(16,'(A,A)') ' *** ERROR in ',CLINE
      STOP

      E N D

      S U B R O U T I N E TNOFF(L, CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  TNOFF
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     terminates with the number L and the
C     message 'CLINE'

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER*(*) CLINE

C%%      write(*,*)'tnoff'

      WRITE( *,'(A,I6,2A)') ' *** ERROR at LINE:',L,' in ',CLINE
      WRITE(16,'(A,I6,2A)') ' *** ERROR at LINE:',L,' in ',CLINE
      STOP

      E N D

      S U B R O U T I N E WARN(CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  WARN
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     outputs a warning message "CLINE"

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER*(*) CLINE

C%%      write(*,*)'warn'

      WRITE( *,'(2A)') ' +++ WARNING in ',CLINE
      WRITE(16,'(2A)') ' +++ WARNING in ',CLINE

      RETURN
      E N D

      S U B R O U T I N E INFO(CLINE)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  INFO
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     writes header with 'CLINE' to unit 'IOUT' & stdout

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

C=====local stuff
      CHARACTER*(*) CLINE
      CHARACTER CINF*14

C%%      write(*,*)'info'

      CINF = ' '
      IB   = ICSBEG(CLINE)
      IE   = ICSEND(CLINE)
      L    = IE - IB + 1

      IF(L.EQ.14) THEN
         CINF = CLINE(IB:IE)
      ELSEIF(L.GT.14) THEN
         CINF = CLINE(IB:IB+13)
      ELSE
         IS             = (16-L)/2
         CINF(IS:IS+L-1) = CLINE(IB:IE)
      ENDIF

      WRITE(16,'(/3(A/))')
     +     '         ==========================',
     +     '         ======' // CINF // '======',
     +     '         =========================='
      
      WRITE(*,'(/3(A/))')
     +     '         ==========================',
     +     '         ======' // CINF // '======',
     +     '         =========================='

      RETURN
      E N D


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C34                                                                           C
C@ PROGRAM=SPRINT, MODULE=IO, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      S U B R O U T I N E BACKSTEP(IUNIT,CMARK,IMAX,IERR,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  BACKSTEP
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     winds a unit back to the last occurrence of 'CMARK' within
C     maximal 'IMAX' i/o records.
C     returns to '*' with 'IERR' .ne. 0 on errors

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      PARAMETER(LMM=32)

      CHARACTER*(*) CMARK
      CHARACTER*(LMM) CTEST

C%%      write(*,*)'backstep'

      IERR  = 0
      LMARK = MIN(LEN(CMARK),LMM)
      WRITE(IUNIT,'(A)') '#'

      DO 10 I=1,IMAX
         BACKSPACE (IUNIT, ERR=9999, IOSTAT=IOS)
         READ(IUNIT,'(A)') CTEST

         IF(CTEST(:LMARK).EQ.CMARK(:LMARK)) RETURN
         
         BACKSPACE (IUNIT, ERR=9999, IOSTAT=IOS)
 10   CONTINUE
      
      IERR = -10 000
      RETURN 1

 9999 CONTINUE
      IERR = IOS
      IF(IERR.EQ.-10 000) IERR = -IERR
      RETURN 1

      E N D

      S U B R O U T I N E NEWVNR(FIN,MAXVNR,FOUT,IVNR,LOUT,IERR,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  NEWVNR
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     finds a new version number for filen with name.extension stored in 'FIN'
C     concatenates '.IVNR' to 'FIN' until this fin.vnr is unique.
C     absolutely necessary to make a unix system behave halfway like a
C     real computer :-)

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      CHARACTER*(*) FIN, FOUT, C1*2, C2*3, C3*4

C%%      write(*,*)'newvnr'

      IERR = 0
C      FOUT = ' '

      MAXVVV = MAXVNR           ! maxvnr might be a constant

      IF(MAXVVV.GT.999) MAXVVV = 999
      IF(MAXVVV.LT.1)   MAXVVV =   1

      LIN  = ICSEND(FIN)
      LOUT = LEN(FOUT)

      IF(MAXVVV.GT.99) THEN
         IF(LIN.GT.LOUT+4) GOTO 9991
         IDEC = 3
      ELSEIF(MAXVVV.GT.9) THEN
         IF(LIN.GT.LOUT+3) GOTO 9991
         IDEC = 2
      ELSE
         IF(LIN.GT.LOUT+2) GOTO 9991
         IDEC = 1
      ENDIF

C 10  'CONTINUE'
C=====type 'file.ext.n'
      C1 = '. '

      DO 110 IVNR=1,MIN(9,MAXVVV)
         WRITE(C1(2:2),'(I1)') IVNR
         FOUT = FIN(:LIN)//C1
         INQUIRE(FILE=FOUT, EXIST=BEX)
         IF(.NOT.BEX) GOTO 111
 110  CONTINUE
      IF(IDEC.GT.1) GOTO 20
      GOTO 9992
 111  CONTINUE

      LOUT = LIN+2
      RETURN

 20   CONTINUE
C=====type 'file.ext.nn'
      C2 = '.  '

      DO 210 IVNR=10,MIN(99,MAXVVV)
         WRITE(C2(2:3),'(I2)') IVNR
         FOUT = FIN(:LIN)//C2
         INQUIRE(FILE=FOUT, EXIST=BEX)
         IF(.NOT.BEX) GOTO 211
 210  CONTINUE
      IF(IDEC.GT.2) GOTO 30
      GOTO 9992
 211  CONTINUE

      LOUT = LIN+3
      RETURN

 30   CONTINUE
C=====type 'file.ext.nnn'
      C3 = '.   '

      DO 310 IVNR=100,MIN(999,MAXVVV)
         WRITE(C3(2:4),'(I3)') IVNR
         FOUT = FIN(:LIN)//C3
         INQUIRE(FILE=FOUT, EXIST=BEX)
         IF(.NOT.BEX) GOTO 311
 310  CONTINUE
      GOTO 9992
 311  CONTINUE

      LOUT = LIN+3
      RETURN

C====='ERROR'--handling....

C====='LOUT' to small
 9991 CONTINUE
      IERR = 1
      RETURN 1

C=====more than 'MAXVNR' versions alredy exist
 9992 CONTINUE
      IERR = 2
      RETURN 1

      E N D

      S U B R O U T I N E NEWUNT(IUN)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  NEWUNT
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     finds a new unconnected unit with 10<=IUN<=90

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)
      CHARACTER*10 CDUM

C%%      write(*,*)'newunt'

      IERR = 0

      DO 10 I=10,90
         INQUIRE(UNIT=I, OPENED=BOP, IOSTAT=IERR, ERR=11)
         IF(.NOT.BOP) THEN
            IUN = I
            RETURN
         ENDIF
 10   CONTINUE

      CALL TOFF('NEWUNT, no free unit')

 11   CONTINUE
      IF(IERR.EQ.1) IERR=1000

      WRITE(CDUM,'(A10)') IERR
      CALL TOFF('NEWUNT, system i/o error number : '//CDUM//
     +        ' while seeking unit')

      RETURN
      E N D


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C32                                                                           C
C@ PROGRAM=SPRINT, MODULE=VECTORS, VERSION=1.00
C                                                                             C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL*8 F U N C T I O N DIVI0(DNUMER,DENOMI)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  DIVI0
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     returns  0.D+0 if DNUMER is  z e r o  ,
C     otherwise DNUMER/DENOMI if DENOMI  is  n o n - z e r o  ,
C     otherwise calls TOFF 
C     => 0/0 -> 0   but  x/0 (x.ne.0) -> TOFF

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)
 
C%%      write(*,*)'divi0'
 
      DIVI0 = 0.D+0
      IF(DNUMER.EQ.0.D+0) RETURN
      IF(DENOMI.EQ.0.D+0) CALL TOFF('DIVI0, division by zero')
      DIVI0 = DNUMER/DENOMI
 
      E N D



      S U B R O U T I N E ZERO(N,U)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C@ SUBPROGRAM :  ZERO
C@ AUTHORS    :  G.H.Hoffstaetter, M.Vogt
C@ ARGS/IN    :     
C@ ARGS/INOUT :  
C@ ARGS/OUT   :  
C@ EXTERNAL   :  
C@ IO_UNITS   :  
C@ DESCRIPTION:  
C@ REMARKS    :  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     initializes a vector/matrix/tensor,....
C     of overall length n*'REAL*8' to 0

      IMPLICIT REAL*8 (A,C-H,O-Z)
      IMPLICIT LOGICAL (B)

      DIMENSION U(N)

C%%      write(*,*)'zero'

      DO 10 I=1,N
         U(I) = 0.D+0
 10   CONTINUE

      RETURN
      E N D



