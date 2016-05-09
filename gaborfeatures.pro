
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
  FOR i = 0,u-1 DO BEGIN
    FOR j = 0,v-1 DO BEGIN
      *(gaborResult[i,j]) = CONVol(dimg,*(gaborArray[i,j]),/EDGE_MIRROR);
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
      gaborAbs = REFORM(gaborAbs,gsize[1]*gsize[2],1);

      ; Normalized to zero mean AND unit variance. (if NOT applicable, please comment this line)
      gaborAbs = (gaborAbs-mean(gaborAbs))/stddev(gaborAbs);
      *(featureMap[i,j]) = REFORM(TRANSPOSE(gaborAbs),n,m);
      featureVector[((c-1)*s):(c*s-1)] = gaborAbs;
    ENDFOR
  ENDFOR
  featureMapPTR = PTR_NEW(featureMap);
  featureVectorPTR = PTR_NEW(featureVector);
  RETURN,[featureMapPTR,featureVectorPTR]
END
