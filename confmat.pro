FUNCTION confmat,Y,T

  ysize = SIZE(Y,/dimension);
  n = ysize[0]
  c = ysize[1]
  Tsize=SIZE(T,/dimension);
  n2 = Tsize[0]
  c2 = Tsize[1]



  ;% Assume two classes with 0-1 encoding
  c = 2;
  class1 = WHERE(T gt 0.5);
  TL = MAKE_ARRAY(n, 1,value=1,/double);
  TL(class1) = 2;
  class2 = WHERE(Y gt 0.5);
  Yclass = MAKE_ARRAY(n, 1,VALUE=1,/double);
  Yclass(class2) = 2;

  ;% Compute
  correct = (Yclass EQ TL);
  sum=TOTAL(TOTAL(correct));
  rate=[sum*100/n,sum];

  Carray=MAKE_ARRAY(c,c,VALUE=0,/double);
  FOR i=0,c-1 DO BEGIN
    FOR j=0,c-1 DO BEGIN
      ci = double(i+1)
      cj = double(j+1)
      Carray[i,j] = TOTAL((Yclass EQ cj)*(TL EQ ci));
    ENDFOR
  ENDFOR
  return,Carray
END