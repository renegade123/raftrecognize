PRO raftrecgnize
  PRINT,"hello,world!"
  originimg = read_image("C:\Users\name\IDLWorkspace83\raftrecognize\data\testme.bmp")
  ;originimg = read_image('F:\IDLworkspace\raftrecognize\data\testme.bmp')
  HELP,originimg
  img = originimg[501:800,501:800]
  ;tvscl,img
  im=image(img, TITLE='Raft',/OVERPLOT)
  groundall = read_txt_data_file('C:\Users\name\IDLWorkspace83\raftrecognize\data\groundall.txt');%导入标签
  ;groundall = read_txt_data_file('F:\IDLworkspace\raftrecognize\data\groundall.txt');%导入标签
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
  test_gabor=MAKE_ARRAY(irow,icol,40,VALUE=0,/DOUBLE);
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
  TestX=TRANSPOSE(TestX);%Gabor特征、GLCM特征合并到一个矩阵
  TestX=min_max_norm(0,1,TestX);%特征归一化
  ;%随机选取训练样本
  gaborsize = SIZE(GaborY)
  PRINT,gaborsize[2]
  rand=FIX(gaborsize[2]*RANDOMU(seed,gaborsize[2]));%产生随机数
  trainingnum=CEIL(gaborsize[2]*0.3);  %取30%的点
  index=rand[0:trainingnum-1];%训练样本的对应的序号
  TrainX=TestX[index,*];%选取训练样本的特征
  TrainY=GaborY[*,index];%选取训练样本的标签
  index_one=WHERE(TrainY EQ 1);%训练样本中浮筏的标签
  index_zero=WHERE(TrainY EQ 0);%训练样本中背景的标签
  ;%稀疏表示分类器
  res=MAKE_ARRAY(1,2,VALUE=0,/DOUBLE);
  PredictY=MAKE_ARRAY(1,gaborsize[2],VALUE=0,/DOUBLE);%预测标签
  ;todo:%稀疏表示算法
  FOR i=0,gaborsize[2]-1 DO BEGIN
    print,i
    X = SimulOMP(TestX(i,*), TrainX, 1e-8, 5,1);
    print,i
    seedD_one=TrainX[index_one,*];
    seedD_zero=TrainX[index_zero,*];
    seedX_one=X[index_one,*];
    seedX_zero=X[index_zero,*];
    Res_one=TestX[i,*]-seedD_one##seedX_one;
    Res_zero=TestX[i,*]-seedD_zero##seedX_zero;
    res[0]=norm(Res_zero);
    res[1]=norm(Res_one);
    
    mres=MIN(res,location);
    index_lab = location
    ww = mres
    PredictY(*,i)=index_lab
  ENDFOR
  map=vectortoimage(irow,icol,PredictY,winsize);%预测标签向量变为矩阵，并上采样
  ;******************************************************************************
  ;%后处理：腐蚀、膨胀
  ;  fg=DOUBLE(bwareaopen(map,100,8));
  ;  SE1=strel('square',8);
  ;  SE2=strel('square',4);
  ;MORPH_CLOSE和MORPH_OPEN
  fg_dilate=DILATE(fg,SE1);%膨胀  腐蚀是erode膨胀是dilate
  fg_erode=ERODE(fg_dilate,SE2);%腐蚀
  ;  map2=fg_erode(4:row+3,4:col+3);%最终分类结果
  ground=image(groundall);
  map = image(map);
  ;****************************************************************************
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
      *(gaborResult[i,j]) = CONVOL(img,*(gaborArray[i,j]),/EDGE_ZERO);
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
  FOR r=0,(FLOOR(image_size/sz)-1)*sz,sz DO BEGIN
    FOR c=0,(FLOOR(image_size2/sz)-1)*sz,sz DO BEGIN
      totalsamples = totalsamples + 1
      FOR i=0,num_images-1 DO BEGIN
        X[i,totalsamples-1]=mean(REFORM(IMAGES[r:r+sz-1,c:c+sz-1,i],sz^2,1));
      ENDFOR
      labelreform = REFORM(Label[r:r+sz-1,c:c+sz-1],sz^2,1)
      Y[*,totalsamples-1]=ROUND(TOTAL(TOTAL(labelreform))/(sz^2));
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
  IF max_value LE min_value THEN BEGIN
    ;DIALOG_MESSAGE('max value can"t be lower than min value');
  END
  size_x=SIZE(x);
  y=MAKE_ARRAY(size_x[1],size_x[2],value=0,/double)
  FOR col=0,size_x[2]-1 DO BEGIN
    max_col=MAX(x[*,col]);
    min_col=MIN(x[*,col]);
    FOR line=0,size_x[1]-1 DO BEGIN
      IF max_col EQ min_col THEN BEGIN
        y[line,col]=(max_value+min_value)/2;
      ENDIF ELSE BEGIN
        y[line,col]=((x[line,col]-min_col)/(max_col-min_col))*(max_value-min_value)+min_value;
      END
    ENDFOR
  ENDFOR
  RETURN,y
END
FUNCTION kappa,confusion_max
  ;[a,b]=SIZE(confusion_max);
  n=TOTAL(TOTAL(confusion_max));
  u=TOTAL(confusion_max,2);
  v=TOTAL(confusion_max);
  kc=(n*TOTAL(DIAG_MATRIX(confusion_max))-v*u)/(n*n-v*u);
  RETURN,kc
END
FUNCTION vectortoimage,image_size,image_size2, label, winsize
  ;num_patches=size(label,2);
  sz= winsize;
  Y= MAKE_ARRAY(image_size,image_size2,value=0,/double)
  totalsamples = 0;
  FOR r=0,(FLOOR(image_size/sz)-1)*sz,sz DO BEGIN
    FOR c=0,(FLOOR(image_size2/sz)-1)*sz,sz DO BEGIN
      totalsamples = totalsamples + 1;
      Y[r:r+sz-1,c:c+sz-1]= label[*,totalsamples-1];
    ENDFOR
  ENDFOR
  RETURN,Y
END
FUNCTION SimulOMP,S, Phi, sigma, T, normType
  ;  % Simultaneous OMP, specify the type OF norm used FOR selecting the atoms
  ;  % based on paper (normType = 1 in Tropp's algorithm)
  ;  % "Algorithms for Simultaneous Sparse Approximation Part I: Greedy  Pursuit"
  ;  % J. Tropp, A. Gilbert, and M. Strauss
  ;  % MIN ||Coeff||_{row,0}   subject to: S = Phi * Coeff
  ;  % Input:       S -- signal to be approximated, d x K matrix
  ;  %            Phi -- DICTIONARY
  ;  %          sigma -- error tolerance
  ;  %              T -- number OF iterations
  ;  %       normType -- type OF norm, 1 FOR L1, 2 FOR L2, 3 FOR infinity norm
  ;  % Ouptut:  Coeff -- coefficients, a sparse matrix with only T nonzero rows
  ;  %         indSet -- index set FOR selected atoms (common support)
  ;  %         Approx -- approximation OF S
  ;  %            Res -- residuals
  ;
  ;  % normalizing columns OF the DICTIONARY
  ;  %Phi0 = Phi;
  ;  %Phi = zeros(SIZE(Phi0));
  phiSize = SIZE(Phi)
  sSize = SIZE(S)
  N = phiSize[1]; % the number of atoms
  d = sSize[1];
  K = sSize[2];
  Coeff = MAKE_ARRAY(N,d,VALUE=0,/DOUBLE)
  Res = S;
  indSet = MAKE_ARRAY(T,1,VALUE=0,/DOUBLE)

  iter = 0;
  norm_res = MAKE_ARRAY(T,1,VALUE=0,/DOUBLE);
  WHILE ((norm(Res[*]) GT sigma) && (iter LT T)) DO BEGIN
    ;% compute the projection
    IF normType EQ 1 THEN BEGIN
      temp = SORT(TOTAL(ABS(TRANSPOSE(Phi)##Res),1))
      dim = reverse(temp)
    ENDIF ELSE BEGIN
      IF normType EQ 2 THEN temp = reverse(SORT(TOTAL(ABS(TRANSPOSE(Phi)##Res)^2,2))) ELSE temp = reverse(SORT(MAX(ABS(TRANSPOSE(Phi)##Res), DIMENSION=2)))
    ENDELSE
    ;% update the index set
    val = temp[0]
    idx = temp
    indSet[iter] = idx[0];
    Coeff[indSet[0:iter],*] =  pinv(Phi[indSet[0:iter],*])##S;
    Approx = Phi##Coeff;
    Res = S - Approx;
    norm_res(iter) = norm(Res[*]);
    iter = iter + 1;
  END
  indSet = indSet[1:iter-1];
  RETURN,Coeff
END
;todo pinv
FUNCTION pinv,a
  apinv = invert(transpose(a)##a)##transpose(a)
  return,apinv
END
