C   15/11/79 311071918  MEMBER NAME  SOL6     (S)           FORTRAN
        SUBROUTINE SOL6(XX,YY,A)
        IMPLICIT REAL*8(A-H,O-Z)
        DATA ICALL/0/
        DIMENSION A(4,4)
C=======BARBER JULY 1982: THIS ROUTINE GOES CRAZY IF THE STRENGTH IS
C=======ZERO. SO IF XX=0,EXPAND THE (1-COS) & SIN AND CANCEL XX TO
C=======GET A TRANSFER MATRIX WHICH IS AS FOR A DRIFT SPACE.
        CALL UNITM(4,A)
        XK=XX/YY
        XC=DCOS(XX)
        XS=DSIN(XX)
        A(1,2)=YY
        IF(XX.NE.0.D0)A(1,2)=XS/XK
        A(1,4)=0.D0
        IF(XX.NE.0.D0)A(1,4)=(1.D0-XC)/XK
        A(2,2)=XC
        A(2,4)=XS
        A(3,2)=-A(1,4)
        A(3,4)=A(1,2)
        A(4,2)=-XS
        A(4,4)=XC
C       IF(ICALL.EQ.0)WRITE(53,9)
C   9   FORMAT(' ','SOL6 TRANSFER MATRIX')
C       IF(ICALL.EQ.0)WRITE(53,10)A
C  10   FORMAT(' ',4F10.4)
C       ICALL=1
        RETURN
        END
