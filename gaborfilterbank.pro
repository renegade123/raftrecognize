;TODO gabor滤波器
FUNCTION gaborFilterBank,u,v,m,n
  ; GABORFILTERBANK generates a custum Gabor filter bank.
  ; It creates a u by v array, whose elements are m by n matries;
  ; each matrix being a 2-D Gabor filter.
  ;
  ;
  ; Inputs:
  ;       u : No. of scales (usually set to 5)
  ;       v : No. of orientations (usually set to 8)
  ;       m : No. of rows in a 2-D Gabor filter (an odd integer number usually set to 39)
  ;       n : No. of columns in a 2-D Gabor filter (an odd integer number usually set to 39)
  ;
  ; Output:
  ;       gaborArray: A u by v array, element OF which are m by n
  ;                   matries; each matrix being a 2-D Gabor filter
  ;
  ;
  ; Sample use:
  ;
  ; gaborArray = gaborFilterBank(5,8,39,39);
  gaborArray = PTRARR(u,v,/ALLOCATE_HEAP);
  fmax = 0.25;
  gama = SQRT(2);
  eta = SQRT(2);

  FOR i = 0,u-1 DO BEGIN
    fu = fmax/((SQRT(2))^i);
    alpha = fu/gama;
    beta = fu/eta;
    FOR j = 0,v-1 DO BEGIN
      tetav = (DOUBLE(j)/DOUBLE(v))*!pi;
      gFilterReal = MAKE_ARRAY(m,n,VALUE=0,/double);
      gFilterImaginary = MAKE_ARRAY(m,n,VALUE=0,/double);
      gFilter = MAKE_ARRAY(m,n,VALUE=0,/double);
      FOR x = 0,m-1 DO BEGIN
        FOR y = 0,n-1 DO BEGIN
          xprime = (x+1-((m+1)/2))*COS(tetav)+(y+1-((n+1)/2))*SIN(tetav);
          yprime = -(x+1-((m+1)/2))*SIN(tetav)+(y+1-((n+1)/2))*COS(tetav);
          ;gFilterReal[x,y] = (fu^2/(!pi*gama*eta))*EXP(-((alpha^2)*(xprime^2)+(beta^2)*(yprime^2)))*EXP((i+1)*2*!pi*fu*xprime);

          gFilterReal[x,y] = ((alpha^2)*(xprime^2)+(beta^2)*(yprime^2))
          gFilterImaginary[x,y] = 2*!pi*fu*xprime
          gFilter[x,y] =(fu^2/(!pi*gama*eta))*EXP(DCOMPLEX(-gFilterReal[x,y],gFilterImaginary[x,y]))
        ENDFOR
      ENDFOR
      *(gaborArray[i,j]) = gFilter;
    ENDFOR
  ENDFOR
  RETURN,gaborArray
END