subroutine string_trim (in_string, out_string, word_len)



implicit none



character(*) in_string, out_string

character(len(in_string)) temp_str   ! In case in_str and out_str are the same actual arg

character(1) tab



integer i, j, word_len, len_in, len_out



parameter (tab = char(9))



! find number of leading blenks



len_in = len(in_string)

len_out = len(out_string)



do i = 1, len_in

  if (in_string(i:i) .ne. ' ' .and. in_string(i:i) .ne. tab) exit

enddo



! If IN_STRING is entirely blanks and/or tab characters then ...



if (i == len_in+1) then

  word_len = 0

  out_string = ' '

  return

endif



! Left shift in_string and put in out_string.



j = min(len_in, len_out + i - 1)

temp_str = in_string(i:j)

out_string = temp_str



! Count characters in first word



do i = 1, len_out

  if (out_string(i:i) .eq. ' ' .or. out_string(i:i) .eq. tab) then

    word_len = i - 1

    return

  endif

enddo



! here if no blanks or tabs



word_len = len_out



end subroutine
