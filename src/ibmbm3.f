

C=====Now get the integrated small spin rotations for each section. 
C     For linear beam-beam effect, just use the linear section maps. 
C     Otherwise use the round or elliptical beam-beam maps.

      IF(IBBSECT(II).EQ.1)THEN    !If the section is a type 17 element. 
      IF(MODIBMBM.GE.2)THEN       !If we want full nonlinear b-b.
                                    !If IBMBM=0, this is skipped. 
      CALL NLBEAMBEAM9(II,XXBB(II),X2BB(II),YYBB(II),NU,ZWBB,SORBVEC,
     +                                   SPINKICKA,NSOLBB(II),MODIBMBM)
C    Or test that we get exactly the same results as when using a type 3 
C    b-b quad, by inserting the linear b-b SECTMAP that is always set up. 
C    even if it's a type 17 element. If IBMBM=0, NLBEAMBEAM is skipped.
C    That's fine since the SECTMAP = identity.
      ELSE
      GMAT9= SECTMAP(7:9,1:6,II)
      SPINKICKA(:,1:NPART3) = MATMUL(GMAT9,SORBVEC(1:6,1:NPART3))
      IF(MODES.EQ.1)THEN
      SPINKICK1(:,1:NPART3) = MATMUL(GMAT9,ORBVEC1(1:6,1:NPART3))
      SPINKICK2(:,1:NPART3) = MATMUL(GMAT9,ORBVEC2(1:6,1:NPART3))
      SPINKICK3(:,1:NPART3) = MATMUL(GMAT9,ORBVEC3(1:6,1:NPART3))
      ENDIF
C    End of test.
      ENDIF      


      ELSE
      GMAT9= SECTMAP(7:9,1:6,II)

C      IF(IDNT.GT.1000)GMAT9 = GMAT9 *0.D0
C      GMAT9(9,1:6)  = 0.D0


C    All orbit modes together 
      SPINKICKA(:,1:NPART3)  = MATMUL(GMAT9,SORBVEC(1:6,1:NPART3))
C    Single orbit modes
      IF(MODES.EQ.1)THEN
      SPINKICK1(:,1:NPART3)  = MATMUL(GMAT9,ORBVEC1(1:6,1:NPART3))
      SPINKICK2(:,1:NPART3)  = MATMUL(GMAT9,ORBVEC2(1:6,1:NPART3))
      SPINKICK3(:,1:NPART3)  = MATMUL(GMAT9,ORBVEC3(1:6,1:NPART3))
      ENDIF
      ENDIF
C     
