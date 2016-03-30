PRO raftrecgnize
  PRINT,"hello,world!"
  originimg = read_image("C:\Users\name\IDLWorkspace83\raftrecognize\data\testme.bmp")
  HELP,originimg
  img = originimg[501:800,501:800]
  im=image(img, TITLE='Raft',/OVERPLOT)
  ;下采样窗大小
  winsize=3
  ;类别两类
  numClasses=2
  isize = SIZE(img)
  ;图片列
  icol = isize[1]
  ;图片行
  irow = isize[2]
  ;设计Gabor滤波器
  gaborArray = gaborFilterBank(5,8,39,39)
  ;提取Gabor特征
  feature = gaborFeatures(img,gaborArray,1,1)
  featureMap = *feature[0]
  featureVector = *feature[1]
  ;将cell存储特征转变为矩阵格式
  test_gabor=MAKE_ARRAY(irow,icol,40);
  t=0;
  FOR i=0,4 DO BEGIN
    FOR j=0,7 DO BEGIN
      test_gabor[*,*,t]=*(featureMap[i,j]);
      t=t+1;
    ENDFOR
  ENDFOR
END
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
      tetav = (j/v)*!pi;
      gFilterReal = make_array(m,n,VALUE=0,/double);
      gFilterImaginary = make_array(m,n,VALUE=0,/double);
      FOR x = 0,m-1 DO BEGIN
        FOR y = 0,n-1 DO BEGIN
          xprime = (x+1-((m+1)/2))*COS(tetav)+(y+1-((n+1)/2))*SIN(tetav);
          yprime = -(x+1-((m+1)/2))*SIN(tetav)+(y+1-((n+1)/2))*COS(tetav);
          gFilterReal[x,y] = (fu^2/(!pi*gama*eta))*EXP(-((alpha^2)*(xprime^2)+(beta^2)*(yprime^2)))*EXP((i+1)*2*!pi*fu*xprime);
          gFilterImaginary[x,y] = 
          gFilter = complex(gFilterReal,gFilterImaginary)
        ENDFOR
      ENDFOR
      *(gaborArray[i,j]) = gFilter;
    ENDFOR
  ENDFOR
  return,gaborArray
END

;todo gabor特征
FUNCTION gaborFeatures,img,gaborArray,d1,d2
  img = DOUBLE(img);
  ;; Filtering
  ; Filter input image by each Gabor filter
  gsize= SIZE(gaborArray);
  u = gsize[1]
  v = gsize[2]
  gaborResult = PTRARR(u,v,/ALLOCATE_HEAP);
  FOR i = 0,u-1 DO BEGIN
    FOR j = 0,v-1 DO BEGIN
      *(gaborResult[i,j]) = convol(img,*(gaborArray[i,j]));
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
  featureVector = make_array(l,1,value=0,/double);
  c = 0;
  FOR i = 0,u-1 DO BEGIN
    FOR j = 0,v-1 DO BEGIN
      c = c+1;
      gaborAbsptr = *(gaborResult[i,j])
      gaborAbs = ABS(gaborAbsptr);
;      gaborAbs = downsample(gaborAbs,d1);
;      gaborAbs = downsample(TRANSPOSE(gaborAbs),d2);
      gaborAbs = congrid(gaborAbs,300,300)
      gsize = size(gaborAbs)
      gaborAbs = REFORM(gaborAbs,gsize[1]*gsize[2],1);

      ; Normalized to zero mean AND unit variance. (if NOT applicable, please comment this line)
      gaborAbs = (gaborAbs-mean(gaborAbs))/stddev(gaborAbs);
      *(featureMap[i,j]) = REFORM(TRANSPOSE(gaborAbs),n,m);
      featureVector[((c-1)*s):(c*s-1)] = gaborAbs;
    ENDFOR
  ENDFOR
  featureMapPTR = ptr_new(featureMap);
  featureVectorPTR = ptr_new(featureVector);
  return,[featureMapPTR,featureVectorPTR]
END
;+
;NAME:
;  DOWNSAMPLE
;PURPOSE:
;  Downsample 1D, 2D, or 3D array, summing the original array over the
;  domain of the new pixels. Works like REBIN, except that the new size
;  does not have to be an integer multiple of the original size. For most
;  purposes, this procedure is preferable to interpolation when downsampling.
;  This routine will also upsample if desired, in which case the action is
;  equivalent to using CONGRID with cubic convolutional interpolation.
;CALLING SEQUENCE:
;  result = downsample(array, Nx [, Ny [, Nz]])
;RETURN VALUE:
;  The result is an array downsampled to the requested dimensions.
;ARGUMENTS:
;  Nx = Size of first dimension of the result.
;  Ny = Size of the second dimension of the result (required if array is 2-3D,
;     ignored otherwise).
;  Nz = Size of the third dimension of the result (required if array is 3D,
;     ignored otherwise).
;ALGORITHM:
;  First, the array is upsampled to a size that is an integer multiple of
;  (Nx, Ny, Nz) by CONGRID using cubic convolutional interpolation. Then
;  the image is decimated to the desired size using REBIN. Note that for 3D
;  arrays, CONGRID should revert to linear interpolation (I haven't actually
;  tested this).
;MODIFICATION HISTORY:
;  2009-Aug-11  C. Kankelborg
FUNCTION downsample, array, Nx, Ny, Nz
  ;-

  asize = SIZE(array)
  CASE asize[0] OF
    1: BEGIN
      Nx0 = asize[1]
      ratio = CEIL(FLOAT(Nx0)/FLOAT(Nx))
      Nx1 = ratio * Nx
      array1 = congrid(array, Nx1, cubic=-0.5, /center)
      result = REBIN(array1, Nx)
    END
    2: BEGIN
      Nx0 = asize[1]
      ratio = CEIL(FLOAT(Nx0)/FLOAT(Nx))
      Nx1 = ratio * Nx
      Ny0 = asize[2]
      ratio = CEIL(FLOAT(Ny0)/FLOAT(Ny))
      Ny1 = ratio * Ny
      array1 = congrid(array, Nx1, Ny1, cubic=-0.5, /center)
      result = REBIN(array1, Nx, Ny)
    END
    3: BEGIN
      Nx0 = asize[1]
      ratio = CEIL(FLOAT(Nx0)/FLOAT(Nx))
      Nx1 = ratio * Nx
      Ny0 = asize[2]
      ratio = CEIL(FLOAT(Ny0)/FLOAT(Ny))
      Ny1 = ratio * Ny
      Nz0 = asize[3]
      ratio = CEIL(FLOAT(Nz0)/FLOAT(Nz))
      Nz1 = ratio * Nz
      array1 = congrid(array, Nx1, Ny1, Nz1, cubic=-0.5, /center)
      result = REBIN(array1, Nx, Ny, Nz)
    END
    ELSE: MESSAGE,'Input array must be 1D, 2D, or 3D.'
  ENDCASE

  RETURN, result
END
;TODO getarray_mean
function getarray_mean,images,label,winsize
  sz= winsize;
  num_patches=FLOOR(image_size/sz)*FLOOR(image_size2/sz);
  totalsamples = 0;
end
