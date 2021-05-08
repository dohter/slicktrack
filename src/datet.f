      SUBROUTINE DATET

      INTEGER date_time(8),start_date_time(8)
      CHARACTER*10 clock(3)
      CHARACTER*3 cmonth
      SAVE ICALL,start_date_time
      DATA ICALL/0/

      CALL date_and_time(clock(1),clock(2),clock(3),date_time)

      IF(ICALL.EQ.0)start_date_time = date_time

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
      ELSE
         cmonth ='???'
      ENDIF

C      WRITE(53,101) date_time(3),cmonth,date_time(1),
C     &     date_time(5),date_time(6),date_time(7)

      WRITE(53,103)
      WRITE(*, 103)
  103 FORMAT(/,/,/,/,/,'  ')

      WRITE(53,100) start_date_time(3),cmonth,start_date_time(1),
     &     start_date_time(5),start_date_time(6),start_date_time(7)
      WRITE(* ,100) start_date_time(3),cmonth,start_date_time(1),
     &     start_date_time(5),start_date_time(6),start_date_time(7)
C
      IF(ICALL.EQ.1)THEN
      WRITE(53,100) date_time(3),cmonth,date_time(1),
     &     date_time(5),date_time(6),date_time(7)


      IF(date_time(5) - start_date_time(5).LT.0)THEN
      start_date_time(5) = start_date_time(5) - 24
      ENDIF
      RUNTIME = 3600 * (date_time(5) - start_date_time(5)) +
     +            60 * (date_time(6) - start_date_time(6)) +
     +                 (date_time(7) - start_date_time(7))
      RUNTIME = RUNTIME/3600.D0

      WRITE(53,'(T47,A,T105,F10.2,A)')
     +                        'Elapsed run time :',RUNTIME,' hours'
      WRITE(* ,'(T47,A,T105,F10.2,A)')
     +                        'Elapsed run time :',RUNTIME,' hours'
      ENDIF


  100 FORMAT('0',45X,'SLICKTRACK run on  '
     & , i2,'-',A3,'-',i4, ' at ',i2,':',i2,':',i2)



      WRITE(53,103)

      ICALL = 1



      RETURN
      END

