FUNCTION convolve, image, psf, FT_PSF=psf_FT, FT_IMAGE=imFT, NO_FT=noft, $
  CORRELATE=correlate, AUTO_CORRELATION=auto, $
  NO_PAD = no_pad
  ;+
  ; NAME:
  ;       CONVOLVE
  ; PURPOSE:
  ;       Convolution of an image with a Point Spread Function (PSF)
  ; EXPLANATION:
  ;       The default is to compute the convolution using a product of
  ;       Fourier transforms (for speed).
  ;
  ;       The image is padded with zeros so that a large PSF does not
  ;       overlap one edge of the image with the opposite edge of the image.
  ;
  ;       This routine is now partially obsolete due to the introduction of  the
  ;       intrinsic CONVOL_FFT() function in IDL 8.1
  ;
  ; CALLING SEQUENCE:
  ;
  ;       imconv = convolve( image1, psf, FT_PSF = psf_FT )
  ;  or:
  ;       correl = convolve( image1, image2, /CORREL )
  ;  or:
  ;       correl = convolve( image, /AUTO )
  ;
  ; INPUTS:
  ;       image = 2-D array (matrix) to be convolved with psf
  ;       psf = the Point Spread Function, (size < or = to size of image).
  ;
  ;       The PSF *must* be symmetric about the point
  ;       FLOOR((n_elements-1)/2), where n_elements is the number of
  ;       elements in each dimension.  For Gaussian PSFs, the maximum
  ;       of the PSF must occur in this pixel (otherwise the convolution
  ;       will shift everything in the image).
  ;
  ; OPTIONAL INPUT KEYWORDS:
  ;
  ;       FT_PSF = passes out/in the Fourier transform of the PSF,
  ;               (so that it can be re-used the next time function is called).
  ;       FT_IMAGE = passes out/in the Fourier transform of image.
  ;
  ;       /CORRELATE uses the conjugate of the Fourier transform of PSF,
  ;               to compute the cross-correlation of image and PSF,
  ;               (equivalent to IDL function convol() with NO rotation of PSF)
  ;
  ;       /AUTO_CORR computes the auto-correlation function of image using FFT.
  ;
  ;       /NO_FT overrides the use of FFT, using IDL function convol() instead.
  ;               (then PSF is rotated by 180 degrees to give same result)
  ;
  ;       /NO_PAD - if set, then do not pad the image to avoid edge effects.
  ;               This will improve memory and speed of the computation at the
  ;               expense of edge effects.   This was the default method prior
  ;               to October 2009
  ; METHOD:
  ;       When using FFT, PSF is centered & expanded to size of image.
  ; HISTORY:
  ;       written, Frank Varosi, NASA/GSFC 1992.
  ;       Appropriate precision type for result depending on input image
  ;                               Markus Hundertmark February 2006
  ;       Fix the bug causing the recomputation of FFT(psf) and/or FFT(image)
  ;                               Sergey Koposov     December 2006
  ;       Fix the centering bug
  ;                               Kyle Penner        October 2009
  ;       Add /No_PAD keyword for better speed and memory usage when edge effects
  ;            are not important.    W. Landsman      March 2010
  ;       Add warning when kernel type does not match integer array
  ;             W. Landsman Feb 2012
  ;       Don't force double precision output   W. Landsman July 2014
  ;-
  COMPILE_OPT idl2
  sp = SIZE( psf_FT,/str )  &  sif = SIZE( imFT, /str )
  sim = SIZE( image )


  IF (sim[0] NE 2) || KEYWORD_SET( noft ) THEN BEGIN
    IF KEYWORD_SET( auto ) THEN BEGIN
      MESSAGE,"auto-correlation only for images with FFT",/INF
      RETURN, image
    ENDIF
    dtype = SIZE(image,/type)
    IF dtype LE 3 THEN IF SIZE(psf,/type) NE dtype THEN $
      MESSAGE,/CON, $
      'WARNING - ' + SIZE(psf,/TNAME) +  $
      ' kernel converted to type ' + SIZE(image,/tname)
    IF KEYWORD_SET( correlate ) THEN $
      RETURN, CONVOL( image, psf ) $
    ELSE    RETURN, CONVOL( image, ROTATE( psf, 2 ) )
  ENDIF

  IF KEYWORD_SET(No_Pad) THEN BEGIN

    sc = sim/2  &  npix = N_ELEMENTS( image )
    IF (sif.N_dimensions NE 2) || ((sif.type NE 6) && (sif.type NE 9)) || $
      (sif.dimensions[0] NE sim[1]) || (sif.dimensions[1] NE sim[2]) THEN imFT = FFT( image,-1 )

    IF KEYWORD_SET( auto ) THEN $
      RETURN, SHIFT( npix*real_part(FFT( imFT*CONJ( imFT ),1 )), sc[1],sc[2] )

    IF (sp.N_dimensions NE 2) || ((sp.type NE 6) && (sp.type NE 9)) || $
      (sp.dimensions[0] NE sim[1]) || (sp.dimensions[1] NE sim[2]) THEN BEGIN
      sp = SIZE( psf )
      IF (sp[0] NE 2) THEN BEGIN
        MESSAGE,"must supply PSF matrix (2nd arg.)",/INFO
        RETURN, image
      ENDIF
      Loc = ( sc - sp/2 ) > 0         ;center PSF in new array,
      s = (sp/2 - sc) > 0        ;handle all cases: smaller or bigger
      L = (s + sim-1) < (sp-1)
      psf_FT = CONJ(image)*0 ;initialise with correct size+type according
      ;to logic of conj and set values to 0 (type of psf_FT is conserved)
      psf_FT[ Loc[1], Loc[2] ] = psf[ s[1]:L[1], s[2]:L[2] ]
      psf_FT = FFT( psf_FT, -1, /OVERWRITE )
    ENDIF

    IF KEYWORD_SET( correlate ) THEN $
      conv = npix * real_part(FFT( imFT * CONJ( psf_FT ), 1 ))  $
    ELSE  conv = npix * real_part(FFT( imFT * psf_FT, 1 ))

    sc = sc + (sim MOD 2)   ;shift correction for odd size images.

    RETURN, SHIFT( conv, sc[1], sc[2] )
  ENDIF ELSE BEGIN


    sc = FLOOR((sim-1)/2) & npix = N_ELEMENTS(image)*4.
    ; the spooky factor of 4 in npix is because we're going to pad the image
    ; with zeros

    IF (sif.N_dimensions NE 2) || ((sif.type NE 6) && (sif.type NE 9)) || $
      (sif.dimensions[0] NE sim[1]) || (sif.dimensions[1] NE sim[2]) THEN BEGIN

      ; here is where we make an array with twice the dimensions of image and
      ; pad with zeros -- thanks to Daniel Eisenstein for this fix

      image_big = MAKE_ARRAY(type = sim[sim[0]+1], sim[1]*2, sim[2]*2)
      image_big[0:sim[1]-1,0:sim[2]-1] = image
      imFT = FFT( image_big,-1 )
      npix = N_ELEMENTS(image_big)

    ENDIF

    IF KEYWORD_SET( auto ) THEN BEGIN
      intermed = SHIFT( npix*real_part(FFT( imFT*CONJ( imFT ),1 )), sc[1],sc[2] )
      RETURN, intermed[0:sim[1]-1,0:sim[2]-1]
    ENDIF


    IF (sp.N_dimensions NE 2) || ((sp.type NE 6) && (sp.type NE 9)) OR $
      (sp.dimensions[0] NE sim[1]) || (sp.dimensions[1] NE sim[2]) THEN BEGIN
      sp = SIZE( psf )
      IF (sp[0] NE 2) THEN BEGIN
        MESSAGE,"must supply PSF matrix (2nd arg.)",/INFO
        RETURN, image
      ENDIF
      ; this obfuscated line determines the offset between the center of the
      ; image and the center of the PSF
      Loc = ( sc - FLOOR((sp-1)/2) )  > 0

      psf_image = MAKE_ARRAY(type = sim[sim[0]+1],sim[1]*2,sim[2]*2)
      psf_image[Loc[1]:Loc[1]+sp[1]-1, Loc[2]:Loc[2]+sp[2]-1] = psf
      psf_FT = FFT(psf_image, -1)
    ENDIF

    IF KEYWORD_SET( correlate ) THEN BEGIN
      conv = npix * real_part(FFT( imFT * CONJ( psf_FT ), 1 ))
      conv = SHIFT(conv, sc[1], sc[2])
    ENDIF ELSE BEGIN
      conv = npix * real_part(FFT( imFT * psf_FT, 1 ))
      conv = SHIFT(conv, -sc[1], -sc[2])

    ENDELSE


    RETURN, conv[0:sim[1]-1,0:sim[2]-1]
  ENDELSE
END