PRO raftrecgnize
  PRINT,"hello,world!"
  originimg = read_image("C:\Users\name\IDLWorkspace83\raftrecognize\data\testme.bmp")
  HELP,originimg
  img = originimg[501:800,501:800]
  ;tvscl,img
  im=image(img, TITLE='Raft',/OVERPLOT)
  groundall = read_txt_data_file('C:\Users\name\IDLWorkspace83\raftrecognize\data\groundall.txt');%导入标签
  ;groundall=groundall(1:100,1:100);
  groundall=groundall(501:800,501:800);
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
  ;[TestX TestY]
  Gabor= getarray_mean(test_gabor,groundall,winsize);%Gabor特征下采样，并将每个像素点特征转变为一列存储
  GaborX = *Gabor[0]
  GaborY = *Gabor[1]
  GLCM= getarray_mean_std(img,groundall,winsize);%GLCM特征提取并下采样，并将每个像素点特征转变为一列存储
  GlcmX = *GLCM[0]
  GlcmY = *GLCM[1]
  ;tempGaborX = transpose(GaborX)
  ;tempGlcmX =  transpose(GlcmX)
  TestX=[GaborX,GlcmX];
  TestX=transpose(TestX);%Gabor特征、GLCM特征合并到一个矩阵
  TestX=min_max_norm(0,1,TestX);%特征归一化
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
      gFilterReal = MAKE_ARRAY(m,n,VALUE=0,/double);
      gFilterImaginary = MAKE_ARRAY(m,n,VALUE=0,/double);
      FOR x = 0,m-1 DO BEGIN
        FOR y = 0,n-1 DO BEGIN
          xprime = (x+1-((m+1)/2))*COS(tetav)+(y+1-((n+1)/2))*SIN(tetav);
          yprime = -(x+1-((m+1)/2))*SIN(tetav)+(y+1-((n+1)/2))*COS(tetav);
          ;gFilterReal[x,y] = (fu^2/(!pi*gama*eta))*EXP(-((alpha^2)*(xprime^2)+(beta^2)*(yprime^2)))*EXP((i+1)*2*!pi*fu*xprime);

          gFilterReal[x,y] = ((alpha^2)*(xprime^2)+(beta^2)*(yprime^2))
          gFilterImaginary[x,y] = 2*!pi*fu*xprime
          gFilter =(fu^2/(!pi*gama*eta))*EXP(DCOMPLEX(-gFilterReal,gFilterImaginary))
        ENDFOR
      ENDFOR
      *(gaborArray[i,j]) = gFilter;
    ENDFOR
  ENDFOR
  RETURN,gaborArray
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
      *(gaborResult[i,j]) = CONVOL(img,*(gaborArray[i,j]));
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
      gaborAbs = congrid(gaborAbs,300,300)
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
FUNCTION getarray_mean,images,label,winsize
  sz= winsize;
  imgsize = SIZE(images)
  image_size = imgsize[1]
  image_size2 = imgsize[2]
  num_images = imgsize[3]
  num_patches=FLOOR(image_size/sz)*FLOOR(image_size/sz);
  totalsamples = 0;
  ;% extract subimages at random from this image to make data vector X
  ;% Step through the images
  X= MAKE_ARRAY(num_images,num_patches,value=0,/double)
  Y= MAKE_ARRAY(1,num_patches,value=0,/double)
  ;% Extract patches at random from this image to make data vector X
  FOR r=0,sz,(FLOOR(image_size/sz)-1)*sz+1 DO BEGIN
    FOR c=0,sz,(FLOOR(image_size2/sz)-1)*sz+1 DO BEGIN
      totalsamples = totalsamples + 1
      FOR i=0,num_images-1 DO BEGIN
        X(i,totalsamples)=mean(REFORM(IMAGES[r:r+sz-1,c:c+sz-1,i],sz^2,1));
      ENDFOR
      Y(*,totalsamples)=ROUND(TOTAL(TOTAL(REFORM(Label[r:r+sz-1,c:c+sz-1],sz^2,1)))/(sz^2));
    ENDFOR
  ENDFOR
  XPTR = PTR_NEW(X);
  YPTR = PTR_NEW(Y);
  RETURN,[XPTR,YPTR]
END
FUNCTION getarray_mean_std,stdimg,label,winsize
  sz= winsize;
  imgsize = SIZE(stdimg)
  image_size = imgsize[1]
  image_size2 = imgsize[2]
  num_images = 1
  num_patches=FLOOR(image_size/sz)*FLOOR(image_size/sz);
  totalsamples = 0;
  ;% extract subimages at random from this image to make data vector X
  ;% Step through the images
  X= MAKE_ARRAY(num_images*2,num_patches,value=0,/double)
  Y= MAKE_ARRAY(1,num_patches,value=0,/double)
  ;% Extract patches at random from this image to make data vector X
  FOR r=0,sz,(FLOOR(image_size/sz)-1)*sz+1 DO BEGIN
    FOR c=0,sz,(FLOOR(image_size2/sz)-1)*sz+1 DO BEGIN
      totalsamples = totalsamples + 1
      FOR i=0,num_images-1 DO BEGIN
        X[i,totalsamples]=mean(REFORM(stdimg[r:r+sz-1,c:c+sz-1,i],sz^2,1));
        X[i+num_images,totalsamples]=stddev(REFORM(stdimg[r:r+sz-1,c:c+sz-1,i],sz^2,1));
      ENDFOR
      Y[*,totalsamples]=ROUND(TOTAL(TOTAL(REFORM(Label[r:r+sz-1,c:c+sz-1],sz^2,1)))/(sz^2));
    ENDFOR
  ENDFOR
  XPTR = PTR_NEW(X);
  YPTR = PTR_NEW(Y);
  RETURN,[XPTR,YPTR]
END
FUNCTION min_max_norm,min_value,max_value,x
  ;%normalize each column OF the input matrix x using MIN MAX normalization
  ;%min_value is the lower bound after normalization AND max_value is the upper bound after normalization
  IF max_value LE min_value then begin
    ;DIALOG_MESSAGE('max value can"t be lower than min value');
  END
  size_x=SIZE(x);
  y=MAKE_ARRAY(size_x[1],size_x[2],value=0,/double)
  FOR col=0,size_x[2]-1 DO BEGIN
    max_col=MAX(x[*,col]);
    min_col=MIN(x[*,col]);
    FOR line=0,size_x[1]-1 DO BEGIN
      IF max_col eq min_col THEN BEGIN
        y[line,col]=(max_value+min_value)/2;
      ENDIF ELSE BEGIN
        y[line,col]=((x[line,col]-min_col)/(max_col-min_col))*(max_value-min_value)+min_value;
      END
    ENDFOR
  ENDFOR
END