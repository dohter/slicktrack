      SUBROUTINE DATE_AND_TIME
      CHARACTER CDATE*9,CTIME*8
C23456X89|10**|****|20**|****|30**|****|40**|****|50**|****|60**|****|70|
      CDATE = '???-??-??'
      CTIME = '??:??:??'
      CALL DATE(CDATE)
      CALL TIME(CTIME)
      WRITE(53,100) CDATE,CTIME
  100 FORMAT('0',///,45X,'S L I C K  RUN ON: ',A,
     &     ' at ',A,5x,'Version: SLICK.MAY97',////)
      RETURN                                                                    
      END                                                                       
C%%  end_sub ::  DATE         
