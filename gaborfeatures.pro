
;todo gabor特征
FUNCTION gaborFeatures,img,gaborArray,d1,d2
  img = DOUBLE(img);
  ;; Filtering
  ; Filter input image by each Gabor filter
  gsize= SIZE(gaborArray);
  u = gsize[1]
  v = gsize[2]
  gaborResult = PTRARR(u,v,/ALLOCATE_HEAP);
  dimg = dcomplex(img,MAKE_ARRAY(300,300,value=0,/double))
  gaborReal = MAKE_ARRAY(300,300,value=0,/double)
  gaborIm = MAKE_ARRAY(300,300,value=0,/double)
  gabor1 = MAKE_ARRAY(300,300,value=0,/double)
  gabor2 = MAKE_ARRAY(300,300,value=0,/double)
  gabor3 = MAKE_ARRAY(300,300,value=0,/double)
  gabor4 = MAKE_ARRAY(300,300,value=0,/double)
  
  FOR i = 0,u-1 DO BEGIN
    FOR j = 0,v-1 DO BEGIN
      ;*(gaborResult[i,j]) = CONVol(dimg,*(gaborArray[i,j]),/EDGE_TRUNCATE);
      a = *(gaborArray[i,j])
      gabor1 = convol_fft(real_part(dimg),real_part(*(gaborArray[i,j])))
      gabor2 = convol_fft(real_part(dimg),IMAGINARY(*(gaborArray[i,j])))
      gabor3 = convol_fft(IMAGINARY(dimg),IMAGINARY(*(gaborArray[i,j])))
      gabor4 = convol_fft(IMAGINARY(dimg),real_part(*(gaborArray[i,j])))
      gaborReal = gabor1 - gabor3
      gaborIm = gabor2 + gabor4
      *(gaborResult[i,j]) = dcomplex(gaborReal,gaborIm)
      ;print,"conv2"+1
      ; J{u,v} = filter2(G{u,v},I);
    ENDFOR
  ENDFOR
  ;Feature Extraction
  featureMap = PTRARR(u,v,/ALLOCATE_HEAP)
  ;Extract feature vector from input image
  imgsize = SIZE(img);
  n = imgsize[1]
  m = imgsize[2]
  s = (n*m)/(d1*d2);
  l = s*u*v;
  featureVector = MAKE_ARRAY(l,1,value=0,/double);
  c = 0;
  FOR i = 0,u-1 DO BEGIN
    FOR j = 0,v-1 DO BEGIN
      c = c+1;
      gaborAbsptr = *(gaborResult[i,j])
      gaborAbs = ABS(gaborAbsptr);
      ;      gaborAbs = downsample(gaborAbs,d1);
      ;      gaborAbs = downsample(TRANSPOSE(gaborAbs),d2);
      ;gaborAbs = downsample(gaborAbs,1,1)
      gsize = SIZE(gaborAbs)
      gaborAbs = REFORM(TRANSPOSE(gaborAbs),gsize[1]*gsize[2],1);

      ; Normalized to zero mean AND unit variance. (if NOT applicable, please comment this line)
      gaborAbs = (gaborAbs-mean(gaborAbs))/stddev(gaborAbs);
      ;*(featureMap[i,j]) = REFORM(TRANSPOSE(gaborAbs),n,m);
      *(featureMap[i,j]) = REFORM(TRANSPOSE(gaborAbs),n,m);
      featureVector[((c-1)*s):(c*s-1)] = gaborAbs;
    ENDFOR
  ENDFOR
  featureMapPTR = PTR_NEW(featureMap);
  featureVectorPTR = PTR_NEW(featureVector);
  RETURN,[featureMapPTR,featureVectorPTR]
END
