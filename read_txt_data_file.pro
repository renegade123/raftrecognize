;+
; :Description:
;    Read digital number sotred in a text file, and the
;    separater of the data in each line must be a 'Space' or 'Tab'.
; :Params:
;    infilename : Input filename of the text file.
;
; :Uses:
;    data = read_txt_data_file('c:test.txt')
;-
FUNCTION read_txt_data_file, infilename
  ;Get the number of lines
  nlines = FILE_LINES(infilename)

  OPENR, lun1, infilename, /GET_LUN

  ;Used to store a line
  tmp_str = ''

  ;Get columns of the input file
  READF, lun1, tmp_str
  tmp = STRSPLIT(tmp_str, COUNT = col_count)
  POINT_LUN, lun1, 0

  ;Allocate memory
  data = FLTARR(col_count, nlines)

  row_count = 0L
  WHILE ~EOF(lun1) DO BEGIN
    READF, lun1, tmp_str
    IF ~STRCMP(tmp_str, '') THEN BEGIN
      tmp_str_split = STRSPLIT(tmp_str, /EXTRACT)
      data_line = FLOAT(tmp_str_split)
      data[*, row_count] = data_line
      row_count = row_count + 1
    ENDIF
  ENDWHILE
  FREE_LUN, lun1
  RETURN, data[*, 0 : (row_count - 1)]
END