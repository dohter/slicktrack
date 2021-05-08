
C=====Now transport the various vectors.
C     For linear beam-beam effect, just use the linear section maps. 
C     Otherwise use the round or elliptical beam-beam maps.

      IF(IBBSECT(II).EQ.1)THEN    !If the section is a type 17 element. 
      IF(MODIBMBM.GE.2)THEN       !If we want full nonlinear b-b.
                                        !If IBMBM=0, this is skipped. 
      CALL NLBEAMBEAM(II,XXBB(II),X2BB(II),YYBB(II),NU,ZWBB,SORBVEC,
     +                                             NSOLBB(II),MODIBMBM)
C    Or test that we get exactly the same results as when using a type 3 
C    b-b quad, by inserting the linear b-b SECTMAP that is always set up. 
C    even if it's a type 17 element. If IBMBM=0, NLBEAMBEAM is skipped.
C    That's fine since the SECTMAP = identity.
      ELSE
      SORBVEC(:,1:NPART2) = MATMUL(SECTMAP(:,:,II),SORBVEC(:,1:NPART2))
      IF(MODES.EQ.1)THEN
      GMAT = SECTMAP(7:8,1:6,II)
      SPINVEC1(:,1:NPART2) = MATMUL(GMAT,ORBVEC1(:,1:NPART2))
     +                                      + SPINVEC1(:,1:NPART2)
      SPINVEC2(:,1:NPART2) = MATMUL(GMAT,ORBVEC2(:,1:NPART2))
     +                                      + SPINVEC2(:,1:NPART2)
      SPINVEC3(:,1:NPART2) = MATMUL(GMAT,ORBVEC3(:,1:NPART2))
     +                                      + SPINVEC3(:,1:NPART2)      
      ENDIF
C    End of test.
      ENDIF      


      ELSE
      SORBVEC(:,1:NPART2) = MATMUL(SECTMAP(:,:,II),SORBVEC(:,1:NPART2))
      IF(MODES.EQ.1)THEN
      GMAT = SECTMAP(7:8,1:6,II)
      SPINVEC1(:,1:NPART2) = MATMUL(GMAT,ORBVEC1(:,1:NPART2))
     +                                      + SPINVEC1(:,1:NPART2)
      SPINVEC2(:,1:NPART2) = MATMUL(GMAT,ORBVEC2(:,1:NPART2))
     +                                      + SPINVEC2(:,1:NPART2)
      SPINVEC3(:,1:NPART2) = MATMUL(GMAT,ORBVEC3(:,1:NPART2))
     +                                      + SPINVEC3(:,1:NPART2)      
      ENDIF
      ENDIF
