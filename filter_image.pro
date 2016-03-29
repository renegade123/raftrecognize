FUNCTION filter_image, image, SMOOTH=width_smooth, ITERATE_SMOOTH=iterate, $
  MEDIAN=width_median, ALL_PIXELS=all_pixels, $
  FWHM_GAUSSIAN=fwhm, NO_FT_CONVOL=no_ft, PSF=psf
  ;+
  ; NAME:
  ;       FILTER_IMAGE
  ;
  ; PURPOSE:
  ;       Identical to MEDIAN or SMOOTH but handle edges and allow iterations.
  ; EXPLANATION:
  ;       Computes the average and/or median of pixels in moving box,
  ;       replacing center pixel with the computed average and/or median,
  ;       (using the IDL SMOOTH() or MEDIAN() functions).
  ;       The main reason for using this function is the options to
  ;       also process the pixels at edges and corners of image, and,
  ;       to apply iterative smoothing simulating convolution with Gaussian,
  ;       and/or to convolve image with a Gaussian kernel.    Users might also
  ;       look at the function ESTIMATOR_FILTER() introduced in IDL 7.1.
  ;
  ; CALLING SEQUENCE:
  ;       Result = filter_image( image, SMOOTH=width, MEDIAN = width, /ALL_PIXELS
  ;                               /ITERATE, FWHM =,  /NO_FT_CONVOL)
  ;
  ; INPUT:
  ;       image = 2-D array (matrix)
  ;
  ; OPTIONAL INPUT KEYWORDS:
  ;       SMOOTH = scalar (odd) integer specifying the width of a square box
  ;               for moving average, in # pixels.  /SMOOTH  means use box
  ;               width = 3 pixels for smoothing.
  ;
  ;        MEDIAN = scalar (usually odd) integer specifying the width of square
  ;               moving box for median filter, in # pixels.   /MEDIAN  means use
  ;               box width = 3 pixels for median filter.
  ;
  ;       /ALL_PIXELS causes the edges of image to be filtered as well.   This
  ;               is accomplished by reflecting pixels adjacent to edges outward
  ;               (similar to the /EDGE_WRAP keyword in CONVOL).
  ;               Note that this is a different algorithm from the /EDGE_TRUNCATE
  ;               keyword to SMOOTH or CONVOL, which duplicates the nearest pixel.
  ;
  ;       /ITERATE means apply smooth(image,3) iteratively for a count of
  ;               (box_width-1)/2 times (=radius), when box_width >= 5.
  ;               This is equivalent to convolution with a Gaussian PSF
  ;               of FWHM = 2 * sqrt( radius ) as radius gets large.
  ;               Note that /ALL_PIXELS is automatically applied,
  ;               giving better results in the iteration limit.
  ;               (also, MEDIAN keyword is ignored when /ITER is specified).
  ;
  ;       FWHM_GAUSSIAN = Full-width half-max of Gaussian to convolve with image.
  ;                       FWHM can be a single number (circular beam),
  ;                       or 2 numbers giving axes of elliptical beam.
  ;
  ;       /NO_FT_CONVOL causes the convolution to be computed directly,
  ;               with intrinsic IDL CONVOL function.   The default is to use
  ;               FFT when factors of size are all LE 13.   Note that
  ;               external function convolve.pro handles both cases)
  ;
  ; OPTIONAL INPUT/OUTPUT KEYWORD:
  ;     PSF = Array containing the PSF used during the convolution.   This
  ;           keyword is only active if the FWHM_GAUSSIAN keyword is also
  ;           specified.     If PSF is undefined on input, then upon output it
  ;           contains the Gaussian convolution specified by the FWHM_GAUSSIAN
  ;           keyword.    If the PSF array is defined on input then it is used
  ;           as the convolution kernel,  the value of the  FWHM_GAUSSIAN keyword
  ;           is ignored.      Typically, on a first call set PSF to an undefined
  ;           variable, which can be reused for subsequent calls to prevent
  ;           recalculation of the Gaussian PSF.
  ; RESULT:
  ;       Function returns the smoothed, median filtered, or convolved image.
  ;       If both SMOOTH and MEDIAN are specified, median filter is applied first.
  ;       If only SMOOTH is applied, then output is of same type as input.  If
  ;       either MEDIAN or FWHM_GAUSSIAN is supplied than the output is at least
  ;       floating (double if the input image is double).
  ;
  ; EXAMPLES:
  ;       To apply 3x3 moving median filter and
  ;       then 3x3 moving average, both applied to all pixels:
  ;
  ;               Result = filter_image( image, /SMOOTH, /MEDIAN, /ALL )
  ;
  ;       To iteratively apply 3x3 moving average filter for 4 = (9-1)/2 times,
  ;       thus approximating convolution with Gaussian of FWHM = 2*sqrt(4) = 4 :
  ;
  ;               Result = filter_image( image, SMOOTH=9, /ITER )
  ;
  ;       To convolve all pixels with Gaussian of FWHM = 3.7 x 5.2 pixels:
  ;
  ;               Result = filter_image( image, FWHM=[3.7,5.2], /ALL )
  ;
  ; EXTERNAL CALLS:
  ;       function psf_gaussian
  ;       function convolve
  ;       pro factor
  ;       function prime          ;all these called only if FWHM is specified
  ;
  ; PROCEDURE:
  ;       If both /ALL_PIXELS (or /ITERATE)  keywords are set then
  ;       create a larger image by reflecting the edges outward, then call the
  ;       IDL MEDIAN() or SMOOTH() function on the larger image, and just return
  ;       the central part (the original size image).
  ;
  ;       NAN values are recognized during calls to MEDIAN() or SMOOTH(), but
  ;       not for convolution with a Gaussian (FWHM keyword supplied).
  ; HISTORY:
  ;       Written, 1991, Frank Varosi, NASA/GSFC.
  ;       FV, 1992, added /ITERATE option.
  ;       FV, 1993, added FWHM_GAUSSIAN= option.
  ;       Use /EVEN call to median, recognize NAN values in SMOOTH
  ;                  W. Landsman   June 2001
  ;       Added PSF keyword,   Bjorn Heijligers/WL, September 2001
  ;       Keep same output data type if /ALL_PIXELS supplied A. Steffl Mar 2011
  ;-
  COMPILE_OPT idl2

  IF N_PARAMS() LT 1 THEN BEGIN
    PRINT,'Syntax - Result = filter_image( image, SMOOTH=width, /ALL_PIXELS'
    PRINT,'                 MEDIAN= width, ITERATE, FWHM=,  /NO_FT_CONVOL'
    RETURN, -1
  ENDIF

  sim = SIZE( image )
  Lx = sim[1]-1
  Ly = sim[2]-1

  IF (sim[0] NE 2) || (sim[4] LE 4) THEN BEGIN
    MESSAGE,"input must be an image (a matrix)",/INFO
    RETURN,image
  ENDIF

  IF KEYWORD_SET( iterate ) THEN BEGIN
    IF N_ELEMENTS( width_smooth ) NE 1 THEN RETURN,image
    IF (width_smooth LT 1) THEN RETURN,image
    imf = image
    nit = (width_smooth>3)/2
    FOR i=1,nit DO  imf = filter_image( imf, /SMOOTH, /ALL )
    RETURN,imf
  ENDIF

  box_wid = 0
  IF KEYWORD_SET( width_smooth ) THEN box_wid = width_smooth > 3
  IF KEYWORD_SET( width_median ) THEN box_wid = (width_median > box_wid)>3

  IF KEYWORD_SET( fwhm ) THEN BEGIN
    npix = ( 3 * fwhm[ 0: ( (N_ELEMENTS( fwhm )-1) < 1 ) ] ) > 3
    npix = 2 * FIX( npix/2 ) + 1    ;make # pixels odd.
    box_wid = box_wid > MAX( [npix] )
  ENDIF

  IF (box_wid LT 3) THEN RETURN, image

  IF KEYWORD_SET(all_pixels) THEN BEGIN

    box_wid = FIX( box_wid )
    radius = (box_wid/2) > 1
    Lxr = Lx+radius
    Lyr = Ly+radius
    rr = 2*radius
    imf = MAKE_ARRAY(sim[1]+rr, sim[2]+rr, type = sim[3])
    imf[radius,radius] = image              ; reflect edges outward
    ; to make larger image.
    imf[  0,0] = ROTATE( imf[radius:rr,*], 5 )      ;Left
    imf[Lxr,0] = ROTATE( imf[Lx:Lxr,*], 5 )         ;right
    imf[0,  0] = ROTATE( imf[*,radius:rr], 7 )      ;bottom
    imf[0,Lyr] = ROTATE( imf[*,Ly:Lyr], 7 )         ;top

  ENDIF ELSE BEGIN
    radius=0
    imf = image
  ENDELSE

  IF KEYWORD_SET( width_median ) THEN $
    imf = MEDIAN(/even, imf, width_median>3 )

  IF KEYWORD_SET( width_smooth ) THEN $
    imf = SMOOTH( imf, width_smooth>3, /NAN )

  IF KEYWORD_SET( fwhm ) THEN BEGIN

    IF N_ELEMENTS( no_ft ) NE 1 THEN BEGIN
      sim = SIZE( imf )
      factor,sim[1],pfx,nfx,/quiet
      factor,sim[2],pfy,nfy,/quiet
      no_ft = MAX( [pfx,pfy] ) GT 13
    ENDIF

    IF N_ELEMENTS(PSF) EQ 0 THEN $
      psf=psf_gaussian( NP=npix,FWHM=fwhm,/NORM )

    imf = convolve( imf,  NO_FT=no_ft, psf)
  ENDIF

  IF radius GT 0 THEN $
    RETURN, imf[ radius:(Lx+radius), radius:(Ly+radius) ] $
  ELSE RETURN, imf
END