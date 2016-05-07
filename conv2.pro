FUNCTION CONV2,a1,a2i,a_interp,PLOT=PLOT,NORMAL=normal,FWHM=fwhm,$
  CORRELATION=correlation,METHOD=method
  ;+
  ; NAME:
  ; CONV2
  ;
  ; PURPOSE:
  ; This function calculates the mathematical convolution or
  ; correlation function of two functions a1(2,npoints1) and
  ; a2(2,npoints2).
  ; Optionally, the second function (a2) may be a gaussian, a lorentzian,
  ; or a combination between them (pseudo-voigt function).
  ; The abscissas arrays of a1 and a2 may be different.
  ;
  ; CATEGORY:
  ; Mathematics.
  ;
  ; CALLING SEQUENCE:
  ; result = conv2(a1,a2 [,ai]) or result = conv2(a1,FWHM=f)
  ;
  ; INPUTS:
  ; a1 the first function. a1 is a two-column array (fltarr(2,npts))
  ;   containing the x array in its first column and the y array
  ;   in its second column.
  ; a2 the second function. Same type of variable as a1. a2 is
  ;   not used when the FWHM keyword is set (see Optional
  ;   Outputs).
  ;
  ; OPTIONAL INPUTS:
  ;
  ; KEYWORD PARAMETERS:
  ; PLOT if set to 1, plot the intermediate interpolated data
  ; NORMAL  If set to 1, normalize the convolution with the integral
  ;    of the FIRST input function: conv = conv/integ(a1).
  ;   If set to 2, normalize the convolution with the integral
  ;    of the SECOND input function: conv = conv/integ(a2).
  ;   If set to 3, normalize the convolution with the integral
  ;    of BOTH input function: conv = conv/(integ(a1)*integ(a2)).
  ; FWHM  set this value to a scalar or to a 2-dim array to
  ;   indicate that the convolution will be done between a1 and
  ;   a pseudo voight function (gaussian, lorentzian or a combination
  ;   between them).
  ;   When defined as float scalar, the value means the FWHM
  ;   of the Gaussian to convolute with.
  ;   When defined as 2-D float array, the first value means the
  ;   FWHM of the Gaussian to convolute with and the second one
  ;   [0<=r<=1 ]the ration between the Gaussian/Lorentzian
  ;   contribution of the pseudo voigt function (1=Pure Lorentzian,
  ;   0=Pure gaussian)
  ;
  ; CORRELATION: When this keyword is set, makes a correlation
  ;   instead of a convolution (i.e. uses INFL routine instead
  ;   of CONV)
  ;
  ;        METHOD: 0 (default) calculates convolution using the FFT
  ;                1 calculates convolution performing the integrals
  ;
  ;
  ; OUTPUTS:
  ; result: a fltarr(2,npoints3)
  ;       the convoluted function defined as:
  ;   CONV(t) = integral[ a1(y) a2(t-y) dy]
  ;       the correlation function defined as:
  ;   CORR(t) = integral[ a1(y) a2(y-t) dy]
  ;
  ;
  ; OPTIONAL OUTPUTS:
  ; a2 is optiomal when FWHM keyword is set. In this case the voigt
  ; function used for the convolution is returned in a named variable.
  ;
  ; ai the interpolated input set
  ;
  ; PROCEDURE:
  ; interpolate both curves (using INTERPOL) to the biggest interval
  ; and smallest step between both abscissas and convolute them by
  ; using CONV (or INFL is the /CORRELATION keyword is set).
  ; Note that that a call conv2(a1,a2) means a coll to conv in the
  ; opposite order: conv(a2,a1).
  ;
  ; EXAMPLE:
  ;
  ;   ; create some data
  ;   x=((findgen(100)/99.)-.5)*.3
  ;   y=exp(-x^2/.01)*abs(sin(x*100))
  ;   a1 = fltarr(2,n_elements(x))
  ;   a1(0,*) = x  &  a1(1,*) = y
  ;
  ;   resul = conv2(a1,a1) ; autoconvolution
  ;   resul = conv2(a1,FWHM=0.1) ; conv with a Gaussian of FWHM=0.1
  ;   resul = conv2(a1,FWHM=[0.1,1]) ;conv with a Lorentzian FWHM=0.1
  ;
  ; MODIFICATION HISTORY:
  ;   Written by: M. Sanchez del Rio
  ; October, 1993
  ; 94-01-21 MSR corrects a bug in the normalization factor:
  ;   (I forgot to calculate the sqrt!!!)
  ; 95-08-30 MSR corrects an error in the definition of the normalization
  ;   factor: it must be 1/integset(set1), it was before
  ;   1/sqrt(integset(set1)*integset(b2))
  ; 97-10-25 srio@esrf.fr adds the convolution with pseudo-voigt
  ;   function option. Adds the CORRELATION kw (thus obsoleting
  ;   infl2) and the METHOD kwyword. Improves doc.
  ;-
  ON_ERROR,2
  ;
  ; interpolation
  ;
  ;if n_params() EQ 1 then a2=a1
  IF KEYWORD_SET(fwhm) THEN a2=a1 ELSE a2=a2i
  xmin = MIN( [MIN(a1(0,*)),MIN(a2(0,*))]  )
  xmax = MAX( [MAX(a1(0,*)),MAX(a2(0,*))]  )
  ;
  ;print,'minima are: ',[a1(0,1)-a1(0,0),a2(0,1)-a2(0,0)]
  step = ABS( MIN([a1(0,1)-a1(0,0),a2(0,1)-a2(0,0)]) )
  ;print,' step = ',step
  npoints = LONG( (xmax-xmin)/step )
  ;print,' npoints = ',npoints
  xnew = FLOAT(FINDGEN(npoints))*step + xmin
  ynew = xnew*0.+1.
  ;
  ; interpolates the new array
  ;
  b1 = FLTARR(2,npoints)
  b1(0,*) = xnew
  b1(1,*) = interpol(a1(1,*),a1(0,*),b1(0,*))

  b2 = FLTARR(2,npoints)
  b2(0,*) = xnew
  IF KEYWORD_SET(fwhm) THEN BEGIN
    IF N_ELEMENTS(fwhm) EQ 1 THEN fwhm=[fwhm,0.0]
    IF fwhm(1) EQ 0 THEN BEGIN
      MESSAGE,/info,'Using a pure Gaussian, Sigma: '+$
        STRCOMPRESS(fwhm(0)/2.35,/rem)
    ENDIF
    IF fwhm(1) EQ 1 THEN MESSAGE,/info,'Using a pure Lorentzian.'
    b2(1,*) = voigt1(b2(0,*),[1.0,0.0,fwhm(0),fwhm(1)])
  ENDIF ELSE BEGIN
    b2(1,*) = interpol(a2(1,*),a2(0,*),b2(0,*))
  ENDELSE
  ;
  ; create the optional output with the interpolated data
  ;
  npar = N_PARAMS()
  IF (npar GT 2) THEN  BEGIN
    a_interp  = FLTARR(3,npoints)
    a_interp(0,*) = xnew
    a_interp(1,*) = b1(1,*)
    a_interp(2,*) = b2(1,*)
  ENDIF
  ;
  IF KEYWORD_SET(plot) THEN BEGIN
    plotfile,b1,title='interpolated data'
    plotfile,b2,kind=2
  ENDIF
  ;
  ; convolution/correlation
  ;
;  IF KEYWORD_SET(correlation) THEN  $
;    resul = infl(b2,b1,method=method) ELSE $
;    resul = conv(b2,b1,method=method)
  ;plotfile,resul,title='convolution'
  ;
  ; normalize the result (if selected)
  ;
  IF KEYWORD_SET(normal) THEN BEGIN
    ;norm_fac = integset(b1) * integset(b2)
    ;resul(1,*) = resul(1,*)/sqrt(norm_fac)
    ;resul(1,*) = resul(1,*)/integset(b1)
    ;resul(1,*) = resul(1,*)/int_tabulated(b1(0,*),b1(1,*))
    CASE normal OF
      1:resul(1,*) = resul(1,*)/int_tabulated(b1(0,*),b1(1,*))
      2:resul(1,*) = resul(1,*)/int_tabulated(b2(0,*),b2(1,*))
      3:resul(1,*) = resul(1,*)/int_tabulated(b1(0,*),b1(1,*)) / $
        int_tabulated(b2(0,*),b2(1,*))
    ENDCASE
  ENDIF
  ;
  RETURN,resul
END