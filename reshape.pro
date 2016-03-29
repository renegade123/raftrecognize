;+
; PURPOSE:
;  This function is a wrapper to the builtin REFORM procedure. It
;  reforms an array to match the shape of template array.
;
; INPUTS:
;  array: The array to reshape. Can be any dimension.
;  template: The array whose shape you wish to match
;
; KEYWORD PARAMETERS:
;  over: If set, array is overwritten on output. This can save
;        significant time and memory with lines like
;        'a = reshape(a, template, /over)' if a is large
;
; OUTPUTS:
;  A reformed version of array whose shape matches that of template.
;
; SIDE EFFECTS:
;  The structure of array is changed if /over is set
;
; MODIFICATION HISTORY
;  April 2010: Written by Chris Beaumont
;-
FUNCTION reshape, array, template, over = over
  COMPILE_OPT idl2
  ON_ERROR, 2

  ;- check inputs
  IF N_PARAMS() NE 2 THEN BEGIN
    PRINT, 'calling sequence'
    PRINT, 'result = reshape(array, template, [/over])'
    RETURN, !values.f_nan
  ENDIF

  IF N_ELEMENTS(array) NE N_ELEMENTS(template) THEN $
    MESSAGE, 'reshape must not change number of elements in array'

  nd = SIZE(template, /n_dim)
  sz = SIZE(template)

  CASE nd OF
    1 : RETURN, REFORM(array, sz[1], over = KEYWORD_SET(over))
    2 : RETURN, REFORM(array, sz[1], sz[2], over = KEYWORD_SET(over))
    3 : RETURN, REFORM(array, sz[1], sz[2], sz[3], over = KEYWORD_SET(over))
    4 : RETURN, REFORM(array, sz[1], sz[2], sz[3], $
      sz[4], over = KEYWORD_SET(over))
    5 : RETURN, REFORM(array, sz[1], sz[2], sz[3], $
      sz[4], sz[5], over = KEYWORD_SET(over))
    6 : RETURN, REFORM(array, sz[1], sz[2], sz[3], $
      sz[4], sz[5], sz[6], over = KEYWORD_SET(over))
    7 : RETURN, REFORM(array, sz[1], sz[2], sz[3], $
      sz[4], sz[5], sz[6], sz[7], over = KEYWORD_SET(over))
    8 : RETURN, REFORM(array, sz[1], sz[2], sz[3], sz[4], $
      sz[5], sz[6], sz[7], sz[8], over = KEYWORD_SET(over))
  ENDCASE
  ;- can't ever get here
  MESSAGE, 'bug in reshape!'
END