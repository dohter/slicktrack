      SUBROUTINE CHROTN(E0,NTY,ICALL)                                           
C=====routine to reset the e3 rotators for use at a different energy.           
C=====parameters from j.buon's 'rotator' file .                                 
C=====see end of member 'rotator' for details. angles are read in rads.         
C     icall=0: read file                                                        
C     icall=1: change energy,find nearby data lines.                            
C     icall=2: set magnets                                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      COMMON/CROT/H4(1500),V3(1500),H3(1500),V2(1500),H2(1500),V1(1500),        
     +     H1(1500),EROT(1500),FC,NL,JL,JL1                                     
      CHARACTER  TNAME*8,A*4                                                    
C=====extra bends for dispersionless s. section files.                          
C     data dba/0.003200120d+0/,dbb/0.d+0/                                       
C     data dba/0.003200120d+0/,dbb/0.000025d+0/                                 
C     data dba/0.d+0/,dbb/0.000025d+0/                                          
      DATA DBA/0.D+0/,DBB/0.D+0/                                                
      INCLUDE "cnlist.for"                                                           
      INCLUDE "clatic.for"                                                           
      INCLUDE "csol.for"                                                             
      DATA IRANGE/0/                                                            
C23456X89|10**|****|20**|****|30**|****|40**|****|50**|****|60**|****|70|       
  
      PI=3.1415926535897932D+0                                                  
      ESET=E0                                                                   
      FH =0.5D0                           !OLD/NEW HERA FILES

C     eset=29.23d+0                                                             
**********************************************************************          
*=====read the rotator data file.                                    *          
**********************************************************************          
  
      IF(ICALL.EQ.0) THEN                                                       
         NL=1                                                                   
  
         DO 12 I=1,1000000                                                      
            READ(51,'(6E13.7)',END=13)                                          
     +           V3(NL),H3(NL),V2(NL),H2(NL),V1(NL),H1(NL)                      
            READ(51,'(6E13.7)')H4(NL),DNU,EROT(NL)                              
            NL=NL+1                                                             
  
 12      CONTINUE                                                               
  
 13      CONTINUE                                                               
**********************************************************************          
*     use the energy eset to select the correct settings for each new*          
*     rotator energy: rotator settings come in about 20mev steps.    *          
*     try interpolating to see how well it performs.                 *          
**********************************************************************          
  
      ELSEIF(ICALL.EQ.1) THEN                                                   
C     jl=(eset+0.001d+0-erot(1))/0.02d+0+1                                      
         JL1=1                                                                  
  
         DO 20 K=1,NL                                                           
  
            IF(ESET.GT.EROT(K))GO TO 20                                         
            JL1=K                                                               
  
            GOTO 21                                                             
  
 20      CONTINUE                                                               
  
 21      CONTINUE                                                               
  
         IF(JL1.EQ.1) THEN                                                      
            WRITE(53,'(a7,F6.2,2a32)') ' ENERGY',ESET,                           
     +           ' OUT OF DESIGN RANGE OF ROTATOR',                             
     +           '--USE THE INPUT SETTINGS'                                     
            IRANGE=1                                                            
  
            RETURN                                                              
  
         ENDIF                                                                  
         JL=JL1-1                                                               
         FC=(ESET-EROT(JL))/(EROT(JL1)-EROT(JL))                                


**********************************************************************          
*     fill in the dipole strengths.                                  *          
**********************************************************************          
  
      ELSEIF(ICALL.EQ.2) THEN                                                   
  
         IF (IRANGE.EQ.1)RETURN                                                 
         NN=NUNT(NTY)                                                           
         TNAME=NAME(NTY)                                                        
         A=TNAME(1:4)                                                           
  
         IF(A.EQ.'BF06')XX(NTY)= (((V1(JL1)-V1(JL))*FC+V1(JL))/NN*1.D+0          
     +        +DBB)*FH                                                              
  
         IF(A.EQ.'BF06'.AND.TNAME(8:8).EQ.'M')XX(NTY)=-XX(NTY)                  


         IF(A.EQ.'BG06')XX(NTY)= (((V3(JL1)-V3(JL))*FC+V3(JL))/NN*1.D+0          
     +        +DBB)*FH                                                              


         IF(A.EQ.'BG06'.AND.TNAME(8:8).EQ.'M')XX(NTY)=-XX(NTY)                  
  

         IF(A.EQ.'BA05') THEN                                                   
            XX(NTY)= (((H1(JL1)-H1(JL))*FC+H1(JL)+DBA)/NN)*FH                        


Cc*   xx(nty)= xx(nty)*1.004d+0                                                 
Cc*   baerr  = xx(nty)*1.004d+0                                                 
         END IF                                                                 
  
         IF(A.EQ.'BB01')XX(NTY)= (((H2(JL1)-H2(JL))*FC+H2(JL))/NN/2.D+0)
     +                                                           *FH          

  
         IF(A.EQ.'BB01'.AND.(TNAME(5:8).EQ.'OR  '                               
     +        .OR.TNAME(5:8).EQ.'OL  ')) THEN                                   
            XX(NTY)=(XX(NTY)+DBA/NN/2.D+0)*FH                                        


            WRITE(53,'(A8,F20.10)') ' BBCHECK',XX(NTY)                          
  
         ENDIF                                                                  
  
         IF(A.EQ.'BC05')XX(NTY)= (((H3(JL1)-H3(JL))*FC+H3(JL))/NN/2.D+0          
     +                                                     *FH)

         IF(A.EQ.'BD05')XX(NTY)= (((H4(JL1)-H4(JL))*FC+H4(JL))/NN/2.D+0          
     +                                                        *FH)


         ANGDEG=XX(NTY)*62.5D+0*180.D+0/PI*NN                                   
  
         IF(NTWIST(NTY).EQ.1)                                                   
     +        WRITE(6,'(2A4,2X,F10.5)') A,TNAME(5:8),ANGDEG                     
C     write(53,'(a7,i10)') 'icall  ',icall                                      
C=====correct total horizontal bend at a harmless place.                        
Cc*   if(aname(1).eq.bd05)xx(nty)= xx(nty)-baerr                                
Cc*   curv=xx(nty)/yy(nty)                                                    
C 200 format(' ','chrotn--curv:',f10.4,f20.10)                                  
  
      ENDIF                                                                     
  
      RETURN                                                                    
      END                                                                       
C%%  end_sub ::  CHROTN(E0,NTY
C23456X89|10**|****|20**|****|30**|****|40**|****|50**|****|60**|****|70|       
  
