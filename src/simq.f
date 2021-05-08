C   15/11/79            MEMBER NAME  SIMQ     (SLIMS)       FORTRAN
       SUBROUTINE SIMQ(A,B,N,KS)
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION A(*),B(*)
       TOL=0.D0
       KS=0
       JJ=-N
       DO 65 J=1,N
C       DO 650 J=1,N                             ! 24/6/09: replaced obsolete code.
       JY=J+1                      
       JJ=JJ+N+1
       BIGA=0.D0
       IT=JJ-J
       DO 30 I=J,N
       IJ=IT+I

C       write(*,*)BIGA,A(IJ)
C      IF(DABS(BIGA)-DABS(A(IJ))) 20,30,30
       IF(DABS(BIGA)-DABS(A(IJ)).LT.0.D0)THEN   ! 24/6/09: replaced obsolete IF code. 
 20    BIGA=A(IJ)                               ! The old version is  simq.f_old_IFs 
       IMAX=I
       ENDIF
 30    CONTINUE

C      IF(DABS(BIGA)-TOL) 35,35,40
       IF(DABS(BIGA)-TOL.LE.0.D0)THEN           ! 24/6/09: replaced obsolete IF code.
 35    KS=1
       RETURN
       ENDIF

 40    I1=J+N*(J-2)
       IT=IMAX-J
       DO 50 K=J,N
       I1=I1+N
       I2=I1+IT
       SAVE=A(I1)
       A(I1)=A(I2)
       A(I2)=SAVE
 50    A(I1)=A(I1)/BIGA
       SAVE=B(IMAX)
       B(IMAX)=B(J)
       B(J)=SAVE/BIGA

C       IF(J-N) 55,70,55
       IF((J-N).EQ.0)GO TO 70                   ! 24/6/09: replaced obsolete IF code.
 55    IQS=N*(J-1)
       DO 65 IX=JY,N
       IXJ=IQS+IX
       IT=J-IX
       DO 60 JX=JY,N
       IXJX=N*(JX-1)+IX
       JJX=IXJX+IT
 60    A(IXJX)=A(IXJX)-(A(IXJ)*A(JJX))
 65    B(IX)=B(IX)-(B(J)*A(IXJ))
C 650   CONTINUE 

 70    NY=N-1
       IT=N*N
       DO 80 J=1,NY
       IA=IT-J
       IB=N-J
       IC=N
       DO 80 K=1,J
       B(IB)=B(IB)-A(IA)*B(IC)
       IA=IA-N
 80    IC=IC-1
       RETURN
       END
