PRO raftrecgnize
  PRINT,"hello,world!"
  originimg = double(read_bmp("C:\Users\name\IDLWorkspace83\raftrecognize\data\testme.bmp"))
  ;originimg = read_image('F:\IDLworkspace\raftrecognize\data\testme.bmp')
  HELP,originimg
  originimg = originimg/255.0
  originimg = transpose(rotate(originimg,1))
  img = originimg[501:800,501:800]
  ;tvscl,img
  ;im=image(img, TITLE='Raft',/OVERPLOT)
  groundall = read_txt_data_file('C:\Users\name\IDLWorkspace83\raftrecognize\data\groundall.txt');%导入标签
  ;groundall = read_txt_data_file('F:\IDLworkspace\raftrecognize\data\groundall.txt');%导入标签
  ;groundall=groundall(1:100,1:100);
  
  groundall=groundall[501:800,501:800];
  aimg = image(img)
  gimg = image(groundall)
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
    PredictY[*,i]=index_lab
  ENDFOR
  map=vectortoimage(irow,icol,PredictY,winsize);%预测标签向量变为矩阵，并上采样
  ;******************************************************************************
  ;%后处理：腐蚀、膨胀
  fg=DOUBLE(MORPH_OPEN(map, REPLICATE(100,8)));
  SE1=strel('square',8);
  SE2=strel('square',4);
  ;MORPH_CLOSE和MORPH_OPEN
  fg_dilate=DILATE(fg,SE1);%膨胀  腐蚀是erode膨胀是dilate
  fg_erode=ERODE(fg_dilate,SE2);%腐蚀
  ;  map2=fg_erode(4:row+3,4:col+3);%最终分类结果
  ground=image(groundall);
  aaa = image(map);
  ;****************************************************************************
END




FUNCTION min_max_norm,min_value,max_value,x
  ;%normalize each column OF the input matrix x using MIN MAX normalization
  ;%min_value is the lower bound after normalization AND max_value is the upper bound after normalization
  IF max_value LE min_value THEN BEGIN
    ;DIALOG_MESSAGE('max value can"t be lower than min value');
  END
  size_x=SIZE(x);
  y=MAKE_ARRAY(size_x[1],size_x[2],value=0,/double)
  FOR col=0,size_x[1]-1 DO BEGIN
    max_col=double(MAX(x[col,*]));
    min_col=double(MIN(x[col,*]));
    FOR line=0,size_x[2]-1 DO BEGIN
      IF max_col EQ min_col THEN BEGIN
        y[col,line]=(max_value+min_value)/2;
      ENDIF ELSE BEGIN
        y[col,line]=((x[col,line]-min_col)/(max_col-min_col))*(max_value-min_value)+min_value;
      END
    ENDFOR
  ENDFOR
  RETURN,y
END
