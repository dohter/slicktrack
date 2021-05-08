      SUBROUTINE SFDAT

      INTEGER date_time(8)
      CHARACTER*10 clock(3)
      CHARACTER*3 cmonth

      CALL date_and_time(clock(1),clock(2),clock(3),date_time)

      IF(date_time(2).eq.1)THEN
         cmonth ='Jan'
      ELSEIF(date_time(2).eq.2)THEN
         cmonth ='Feb'
      ELSEIF(date_time(2).eq.3)THEN
         cmonth ='Mar'
      ELSEIF(date_time(2).eq.4)THEN
         cmonth ='Apr'
      ELSEIF(date_time(2).eq.5)THEN
         cmonth ='May'
      ELSEIF(date_time(2).eq.6)THEN
         cmonth ='Jun'
      ELSEIF(date_time(2).eq.7)THEN
         cmonth ='Jul'
      ELSEIF(date_time(2).eq.8)THEN
         cmonth ='Aug'
      ELSEIF(date_time(2).eq.9)THEN
         cmonth ='Sep'
      ELSEIF(date_time(2).eq.10)THEN
         cmonth ='Oct'
      ELSEIF(date_time(2).eq.11)THEN
         cmonth ='Nov'
      ELSEIF(date_time(2).eq.12)THEN
         cmonth ='Dec'
      ELSE cmonth ='???'
      ENDIF

      WRITE(53,100) date_time(3),cmonth,date_time(1),
     &     date_time(5),date_time(6),date_time(7)
  100 FORMAT('0',///,45X,'S P I N F L I P run on:',i2,'-',A3,'-',i4,
     &     ' at ',i2,':',i2,':',i2////)

      RETURN
      END

