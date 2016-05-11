;+
; 《IDL语言程序设计》
; --数据可视化与ENVI二次开发
;
; 示例程序
;
; 作者: 董彦卿
;
; 联系方式：sdlcdyq@sina.com
;
;todo:销毁析构
PRO raftrecognize::CLEANUP
  IF  PTR_VALID(self.ORIDATA) THEN PTR_FREE,self.ORIDATA
  ;IF  PTR_VALID(self.conData) THEN PTR_FREE,self.conData
  ;IF OBJ_VALID(self.OWINDOW) THEN OBJ_DESTROY, self.OWINDOW
  IF OBJ_VALID(self.OIMAGE) THEN OBJ_DESTROY, self.OIMAGE
  IF OBJ_VALID(self.ORECT) THEN OBJ_DESTROY, self.ORECT
  ;IF OBJ_VALID(self.oContour) THEN OBJ_DESTROY, self.oContour
  ;IF OBJ_VALID(self.INITLINE) THEN OBJ_DESTROY, self.INITLINE
END

;todo:修改组件大小是
PRO raftrecognize::ChangeDrawSize,width,height
  IF N_ELEMENTS(width) THEN BEGIN
    self.OWINDOW.GETPROPERTY, graphics_tree = oView
    oView.GetProperty,ViewPlane_Rect = viewP,dimensions = dims
    oriWL = viewP[2:3]
    viewP[2:3] =viewP[2:3]*[width,height]/dims
    viewP[0:1]+=(oriWL-viewP[2:3])/2

    oView.SETPROPERTY,dimension = [width,height],viewPlane_Rect= viewP
    self.OWINDOW.Draw
  ENDIF
END
;TODO 提取海岸线
PRO raftrecognize::getline,inifile
  ;originImg=READ_IMAGE('E:\IDLWorkspace83\newraftrecognize\12.png')
  ;originImg=READ_IMAGE('E:\IDLWorkspace83\newraftrecognize\cc.JPG')
  ;originImg=READ_tiff(self.infile)
  originImg=READ_TIFF('C:\Users\name\IDLWorkspace83\raftrecognize\subset_VV.tif',R, G, B,GEOTIFF=GeoKeys,INTERLEAVE = 0)
  ;Img=TRANSPOSE(ROTATE(DOUBLE(REFORM(originImg[0,2000:2399,2000:2399])),1))
  ;Img=TRANSPOSE(ROTATE(DOUBLE(originImg),1))
  Img=DOUBLE(originImg)
  HELP,IMG
  ;plot, transpose(img,[0,2,1])
  timestep=5;   time step
  mu=0.2/timestep;coefficient of the distance regularization term R(phi)
  iter_inner=5
  ;iter_outer=150
  iter_outer=300
  lambda=5;coefficient of the weighted length term L(phi)
  ;alfa=1.5;coefficient of the weighted area term A(phi)
  alfa=5000;
  epsilon=1.5;papramater that specifies the width of the DiracDelta function
  sigma=1.5;scale parameter in Gaussian kernel
  ;  Gauss=GAUSSIAN_FUNCTION([sigma,sigma],/double,width=15);gauss核
  ;  Img_smooth=convol(Img,gauss,/CENTER, /EDGE_TRUNCATE);卷积
  Img_smooth = GAUSS_SMOOTH(Img , sigma,/EDGE_ZERO); /EDGE_ZERO gauss平滑
  iGrad=self.gradient(Img_smooth,/vector);
  Ix=iGrad[*,*,0]
  Iy=iGrad[*,*,1]
  f=Ix^2+Iy^2;
  g=1/(1+f);
  c0=2
  imgsize=SIZE(Img)
  initialLSF=c0*MAKE_ARRAY(imgsize[1],imgsize[2],/double,VALUE=1)
  ;initialLSF[0:500,0:370]=-c0
  ;initialLSF[0:289,0:279]=-c0
  ;initialLSF[0:399,40:399]=-c0
  initialLSF[2:1478,2:1559]=-c0
  phi=initialLSF


  ;显示初始level——set
  ;im = IMAGE(originImg[*,2000:2399,2000:2399], RGB_TABLE=13, TITLE='raftrecognize')
  TileData = BYTSCL(originImg)
  self.ORIDATA = PTR_NEW(TileData,/no_Copy)
  ;self.ORIDATA = PTR_NEW(TRANSPOSE(ROTATE(TileData,1)),/no_Copy)
  ;self.CONDATA = PTR_NEW(TRANSPOSE(ROTATE(phi,1)),/no_Copy)
  self.CONDATA = PTR_NEW(initialLSF,/no_Copy)
  idata = *(self.ORIDATA)
  cdata = *(self.CONDATA)
  self.OIMAGE.SETPROPERTY, data = idata
  self.oContour.SETPROPERTY, hide =0,data = cdata
  self.OWINDOW.draw
  ;im = IMAGE(originImg, RGB_TABLE=13, TITLE='raftrecognize',/OVERPLOT)
  ;c=CONTOUR(TRANSPOSE(ROTATE(phi,1)), C_LINESTYLE=0,c_label_show=0,COLOR=[0,255,0] ,c_value=[0,0] , /OVERPLOT)
  potential=2;
  IF potential EQ 1 THEN BEGIN
    potentialFunction = 'single-well';use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model
  ENDIF ELSE BEGIN
    IF potential EQ 2 THEN potentialFunction = 'double-well' ELSE potentialFunction = 'double-well'
  END
  FOR n=1,iter_outer DO BEGIN
    tic
    *(self.CONDATA)=self.drlse_edge(*(self.CONDATA), g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction)
    toc
    IF (n MOD 2) EQ 0 THEN BEGIN
      tic
      ;      c.erase
      ;      im = IMAGE(originImg[*,2000:2399,2000:2399], RGB_TABLE=13, TITLE='raftrecognize',/OVERPLOT)
      ;      ;im = IMAGE(originImg, RGB_TABLE=13, TITLE='raftrecognize',/OVERPLOT)
      ;      c = CONTOUR(TRANSPOSE(ROTATE(phi,1)), C_LINESTYLE=0,c_label_show=0,COLOR=[0,255,0] ,c_value=[0,0] ,/OVERPLOT)
      ;self.CONDATA = PTR_NEW(TRANSPOSE(ROTATE(phi,1)),/no_Copy)
      ;tempphi = phi
      ;self.CONDATA = PTR_NEW(tempphi,/no_Copy)
      cdata = *(self.CONDATA)
      self.OCONTOUR.SETPROPERTY, hide =0,data = cdata
      self.OWINDOW.draw
      toc
    ENDIF
  ENDFOR
  alfa=0;
  iter_refine = 10;
  *(self.CONDATA) = self.drlse_edge(*(self.CONDATA), g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction)
  ;  c.erase
  ;  im = IMAGE(originImg[*,2000:2399,2000:2399], RGB_TABLE=13, TITLE='raftrecognize',/OVERPLOT)
  ;  ;im = IMAGE(originImg, RGB_TABLE=13, TITLE='raftrecognize',/OVERPLOT)
  ;  c = CONTOUR(TRANSPOSE(ROTATE(phi,1)), C_LINESTYLE=0,c_label_show=0,COLOR=[0,255,0] ,c_value=[0,0] ,/OVERPLOT)
  ;self.CONDATA = PTR_NEW(TRANSPOSE(ROTATE(phi,1)),/no_Copy)
  ;self.CONDATA = PTR_NEW(phi,/no_Copy)
  cdata = *(self.CONDATA)
  self.OCONTOUR.SETPROPERTY, hide =0,data = cdata
  self.OWINDOW.draw
  B = WHERE(phi GT 0, count, COMPLEMENT=B_C, NCOMPLEMENT=count_c)
  phi[b] = -1
  phi[b_c] = 1
  WRITE_TIFF,"d:\test.tif",phi
END

;todo 云掩膜
pro raftrecognize::cloudMask
  COMPILE_OPT idl2
  ENVI, /restore_base_save_files
  ENVI_BATCH_INIT
  self.hdffile = dialog_pickfile(title='Pick .hdf File',filter = ['*.hdf'])
  self.hbasename = FILE_BASENAME(self.hdffile, '.hdf')
END

PRO raftrecognize::example_glt_usage
  COMPILE_OPT strictarr
  ENVI, /RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT, LOG_FILE = 'batch.log'
  ;  Select an image and it's associated input geometry.
  ;  convert mask to tiff
  ENVI_SELECT, title='Please select an image to georeference', $
    dims=dims, fid=fid, pos=pos
  IF fid EQ -1 THEN RETURN
  ENVI_SELECT, title='Please select the longitude band', $
    fid=x_fid, dims=dims, pos=x_pos, /band_only, /no_dims
  IF x_fid EQ -1 THEN RETURN
  ENVI_SELECT, title='Please select the latitude band', $
    fid=y_fid, dims=dims, pos=y_pos, /band_only, /no_dims
  IF y_fid EQ -1 THEN RETURN

  ;  Figure out what UTM zone we're in.
  query_dims = [-1L, LONG(dims[2]/2), LONG(dims[2]/2), LONG(dims[4]/2),$
    LONG(dims[4]/2)]
  lat_lon = [ENVI_GET_DATA(fid=y_fid, dims=query_dims, pos=y_pos), $
    ENVI_GET_DATA(fid=x_fid, dims=query_dims, pos=x_pos)]
  zone = FIX(31.0 + lat_lon[1]/6.0)
  south = (lat_lon[0] LT 0)

  ;  Make the GLT.
  out_name=self.filedir+self.basename+'_glt.dat'
  pixel_size=6000.0D
  rotation=0.0
  i_proj = ENVI_PROJ_CREATE(/geographic)
  o_proj = ENVI_PROJ_CREATE(/utm, south=south,ZONE=zone)
  ENVI_DOIT, 'ENVI_GLT_DOIT', i_proj=i_proj, $
    o_proj=o_proj, out_name=out_name, $
    pixel_size=pixel_size, r_fid=glt_fid, $
    rotation=0.0, x_fid=x_fid, y_fid=y_fid, $
    x_pos=x_pos, y_pos=y_pos
  IF glt_fid EQ -1 THEN RETURN

  ;  Georeference the image from the GLT.
  out_name=self.filedir+self.basename+'_georef.tif'
  ENVI_DOIT, 'envi_georef_from_glt_doit', fid=fid, $
    glt_fid=glt_fid, out_name=out_name, pos=pos, $
    subset=dims, r_fid=r_fid
END
;参数设置
PRO raftrecognize::GetProperty, initFlag = initFlag,hdffile = hdffile,hbasename = hbasename 
  initFlag= self.INITFLAG
  hdffile= self.hdffile
  hbasename= self.hbasename
END
;**********************************************************************
;todo主程序：还要拆分
;**********************************************************************
PRO raftrecgnize::superclassfy
  originimg = DOUBLE(read_bmp("C:\Users\name\IDLWorkspace83\raftrecognize\data\testme.bmp"))
  ;originimg = read_image('F:\IDLworkspace\raftrecognize\data\testme.bmp')
  HELP,originimg
  originimg = originimg/255.0
  originimg = TRANSPOSE(ROTATE(originimg,1))
  img = originimg[500:799,500:799]
  ;tvscl,img
  ;im=image(img, TITLE='Raft',/OVERPLOT)
  groundall = read_txt_data_file('C:\Users\name\IDLWorkspace83\raftrecognize\data\groundall.txt');%导入标签
  ;groundall = read_txt_data_file('F:\IDLworkspace\raftrecognize\data\groundall.txt');%导入标签
  ;groundall=groundall(1:100,1:100);

  groundall=groundall[501:800,501:800];

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
    PRINT,i
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
  fg=DOUBLE(bwareaopen(map));
  SE1 = REPLICATE(1, 8, 8);
  SE2 = REPLICATE(1, 4, 4);
  ;MORPH_CLOSE和MORPH_OPEN
  fg_map = MAKE_ARRAY(irow+16,irow+16,VALUE=0,/DOUBLE)
  fg_map[8:irow+7,8:irow+7] = fg
  fg_dilate=DILATE(fg_map,SE1);%膨胀  腐蚀是erode膨胀是dilate
  fg_erode=ERODE(fg_dilate,SE1);%腐蚀
  map2=DOUBLE(fg_erode[8:irow+7,8:irow+7]);%最终分类结果
  aimg = image(img,/CURRENT, LAYOUT=[3,2,1], TITLE='原图')
  gimg = image(groundall,/CURRENT, LAYOUT=[3,2,2], TITLE='groundall')
  amap = image(map,/CURRENT, LAYOUT=[3,2,3], TITLE="初始图像")
  afg = image(fg,/CURRENT, LAYOUT=[3,2,4], TITLE="斑块过滤完成")
  aaa = image(map2,/CURRENT,LAYOUT=[3,2,5],TITLE="最终结果");
  ;****************************************************************************
  ;conf_max2=confmat(reform(groundall[4:irow-4,4:icol-4],N_Elements(groundall[4:irow-4,4:icol-4]),1),reform(map2,N_Elements(map2),1));
  conf_max2=confmat(REFORM(groundall,N_ELEMENTS(groundall),1),REFORM(map2,N_ELEMENTS(map2),1));
  PRINT,conf_max2
  overallacc2=TOTAL(DIAG_MATRIX(conf_max2))/TOTAL(TOTAL(conf_max2))
  every_classacc2=DIAG_MATRIX(conf_max2)/TOTAL(conf_max2,1);
  averageacc2=mean(every_classacc2)
  kapparate2=(conf_max2)
END

;*********************************************************************
; :Description:
;    Read digital number sotred in a text file, and the
;    separater of the data in each line must be a 'Space' or 'Tab'.
; :Params:
;    infilename : Input filename of the text file.
; :Uses:
;    data = read_txt_data_file('c:test.txt')
;**********************************************************************
FUNCTION read_txt_data_file, infilename
  ;Get the number of lines
  nlines = FILE_LINES(infilename)

  OPENR, lun1, infilename, /GET_LUN

  ;Used to store a line
  tmp_str = ''

  ;Get columns of the input file
  READF, lun1, tmp_str
  tmp = STRSPLIT(tmp_str, COUNT = col_count)
  POINT_LUN, lun1, 0

  ;Allocate memory
  data = FLTARR(col_count, nlines)

  row_count = 0L
  WHILE ~EOF(lun1) DO BEGIN
    READF, lun1, tmp_str
    IF ~STRCMP(tmp_str, '') THEN BEGIN
      tmp_str_split = STRSPLIT(tmp_str, /EXTRACT)
      data_line = FLOAT(tmp_str_split)
      data[*, row_count] = data_line
      row_count = row_count + 1
    ENDIF
  ENDWHILE
  FREE_LUN, lun1
  RETURN, data[*, 0 : (row_count - 1)]
END
;**********************************************************************
;TODO gabor滤波器
;**********************************************************************
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
  fmax = 0.25D;
  gama = DOUBLE(SQRT(2));
  eta = DOUBLE(SQRT(2));

  FOR i = 0,u-1 DO BEGIN
    fu = DOUBLE(fmax/((SQRT(2))^i));
    alpha = fu/gama;
    beta = fu/eta;
    FOR j = 0,v-1 DO BEGIN
      tetav = (DOUBLE(j)/DOUBLE(v))*!const.pi;
      gFilterReal = MAKE_ARRAY(m,n,VALUE=0,/double);
      gFilterIM = MAKE_ARRAY(m,n,VALUE=0,/double);
      gFilterImaginary = MAKE_ARRAY(m,n,VALUE=0,/double);
      gFilter =DCOMPLEXARR(m,n)
      FOR x = 0,m-1 DO BEGIN
        FOR y = 0,n-1 DO BEGIN
          xprime = (x+1-((m+1)/2))*COS(tetav)+(y+1-((n+1)/2))*SIN(tetav);
          yprime = -(x+1-((m+1)/2))*SIN(tetav)+(y+1-((n+1)/2))*COS(tetav);
          ;gFilterReal[x,y] = (fu^2/(!pi*gama*eta))*EXP(-((alpha^2)*(xprime^2)+(beta^2)*(yprime^2)))*EXP((i+1)*2*!pi*fu*xprime);

          gFilterReal[x,y] = ((alpha^2)*(xprime^2)+(beta^2)*(yprime^2))
          gFilterImaginary[x,y] = 2*!const.pi*fu*xprime
          gFilter[x,y] =DOUBLE(fu^2/(!const.pi*gama*eta))*EXP(-gFilterReal[x,y])*EXP(DCOMPLEX(0.0,gFilterImaginary[x,y]))
        ENDFOR
      ENDFOR
      *(gaborArray[i,j]) = gFilter;
    ENDFOR
  ENDFOR
  RETURN,gaborArray
END
;**********************************************************************
;todo gabor特征
;**********************************************************************
FUNCTION raftrecgnize::gaborFeatures,img,gaborArray,d1,d2
  img = DOUBLE(img);
  ;; Filtering
  ; Filter input image by each Gabor filter
  gsize= SIZE(gaborArray);
  u = gsize[1]
  v = gsize[2]
  gaborResult = PTRARR(u,v,/ALLOCATE_HEAP);
  dimg = DCOMPLEX(img,MAKE_ARRAY(300,300,value=0,/double))
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
      *(gaborResult[i,j]) = DCOMPLEX(gaborReal,gaborIm)
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
;*********************************************************************
;TODO getarray_mean gabor特征采样
;*********************************************************************
FUNCTION raftrecgnize::getarray_mean,images,label,winsize
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
        X[i,totalsamples-1]=mean(REFORM(IMAGES[r:r+sz-1,c:c+sz-1,i],sz^2,1),/DOUBLE);
      ENDFOR
      labelreform = REFORM(Label[r:r+sz-1,c:c+sz-1],sz^2,1)
      Y[*,totalsamples-1]=ROUND(TOTAL(TOTAL(labelreform))/(sz^2));
    ENDFOR
  ENDFOR
  XPTR = PTR_NEW(X);
  YPTR = PTR_NEW(Y);
  RETURN,[XPTR,YPTR]
END
;*********************************************************************
;TODO getarray_mean GLCM特征采样
;*********************************************************************
FUNCTION raftrecgnize::getarray_mean_std,stdimg,label,winsize
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
  FOR r=0,(FLOOR(image_size/sz)-1)*sz,sz DO BEGIN
    FOR c=0,(FLOOR(image_size2/sz)-1)*sz,sz DO BEGIN
      totalsamples = totalsamples + 1
      FOR i=0,num_images-1 DO BEGIN
        X[i,totalsamples-1]=mean(REFORM(stdimg[r:r+sz-1,c:c+sz-1,i],sz^2,1),/DOUBLE);
        X[i+num_images,totalsamples-1]=stddev(REFORM(stdimg[r:r+sz-1,c:c+sz-1,i],sz^2,1),/DOUBLE);
      ENDFOR
      Y[*,totalsamples-1]=ROUND(TOTAL(TOTAL(REFORM(Label[r:r+sz-1,c:c+sz-1],sz^2,1)))/(sz^2));
    ENDFOR
  ENDFOR
  XPTR = PTR_NEW(X);
  YPTR = PTR_NEW(Y);
  RETURN,[XPTR,YPTR]
END
;******************************************************************************
;TODO 特征归一化
;******************************************************************************
FUNCTION raftrecgnize::min_max_norm,min_value,max_value,x
  ;%normalize each column OF the input matrix x using MIN MAX normalization
  ;%min_value is the lower bound after normalization AND max_value is the upper bound after normalization
  IF max_value LE min_value THEN BEGIN
    ;DIALOG_MESSAGE('max value can"t be lower than min value');
  END
  size_x=SIZE(x);
  y=MAKE_ARRAY(size_x[1],size_x[2],value=0,/double)
  FOR col=0,size_x[1]-1 DO BEGIN
    max_col=DOUBLE(MAX(x[col,*]));
    min_col=DOUBLE(MIN(x[col,*]));
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
;******************************************************************************
;TODO Simultaneous OMP, specify the type of norm used for selecting the atoms
;******************************************************************************
FUNCTION raftrecgnize::SimulOMP,S, Phi, sigma, T, normType
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
    idx = dim
    indSet[iter] = idx[0];
    Coeff[indSet[0:iter],*] =  PSEUDO_INVERSE(Phi[indSet[0:iter],*])##S;
    Approx = Phi##Coeff;
    Res = S - Approx;
    norm_res(iter) = norm(Res[*]);
    iter = iter + 1;
  END
  ;indSet = indSet[0:iter-1];
  RETURN,Coeff
END
;******************************************************************************
;TODO pseudo_inverse实现与matlab中pinv相同的功能
;******************************************************************************
FUNCTION raftrecgnize::pseudo_inverse,A,U,V,S,TOL=tol,SORT=sort,_extra=extra
  ;+
  ;B=pseudo_inverse(A,[U,V,S,tol,EXTRA=extra)
  ;computes the pseudo-inverse using SVDC.
  ;B=pseudo_inverse(A,U,V,S) returns the SVDC values
  ;U,V and S such that A=U##S##Transpose(V). This is
  ;especially useful in examining the sensitivity of
  ;the solution. Small values along the diagonal of
  ;S can be zeroed out in computing the pseudo-inverse
  ;to reduce numerical sensitivity.
  ;
  ;PARAMETERS
  ;Set TOL to a threshold value on the singular values
  ;to be accepted in computing the inverse. Example
  ;Ap=pseudo_inverse(A,TOL=1e-5) causes singular values
  ;smaller than 1e-5 to be excluded from the inverse
  ;calculation. Default: tol=1e-10
  ;
  ;Set SORT to sort the singular values in descending
  ;order before computing the inverse.
  ;
  ;2004-05-13 HR Added SORT and TOL keywords and replaced
  ;use of the special DIAG() function with DIAG_MATRIX()
  ;
  ;Reference:
  ;Moon and Stirling,
  ;Mathematical Methods and Algorithms
  ;Prentice-Hall, 2000
  ;Section 7.4
  ;
  ;H. Rhody  December 30, 2003
  ;-
  type=SIZE(A,/TYPE)
  dim=SIZE(A,/DIM)

  IF N_ELEMENTS(dim) EQ 1 THEN BEGIN
    ;Case of a row vector
    ;Since SVDC does not handle 1D arrays, we
    ;can coerce a solution by making a column vector
    ;and extracting a solution from it.
    SVDC,TRANSPOSE(A),S,U1,V,_EXTRA=extra
    U=V
    V=U1
  ENDIF ELSE SVDC,A,S,U,V,_extra=extra

  IF N_ELEMENTS(tol) LE 0 THEN tol=1e-10

  IF KEYWORD_SET(sort) THEN BEGIN
    I=REVERSE(SORT(S))
    S=S[I]
    U=U[I,*]
    V=V[I,*]
  END

  W=MAKE_ARRAY(SIZE=SIZE(S),VALUE=0)
  I=WHERE(ABS(S) GE tol)
  IF MIN(I) LT 0 THEN MESSAGE,'No singular values larger than TOL'

  W[I]=1.0/S[I]

  W=DIAG_MATRIX(W)

  IF type LE 5 THEN $
    RETURN,V##W##TRANSPOSE(u) ELSE $
    RETURN,V##W##TRANSPOSE(CONJ(U))
END
;*********************************************************************
;todo:vectortoimage
;向量转矩阵
;*********************************************************************
FUNCTION raftrecgnize::vectortoimage,image_size,image_size2, label, winsize
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
;*********************************************************************
;todo:实现与 Matlab 的 bwareaopen 相同的功能
;过滤斑块
;*********************************************************************
FUNCTION raftrecognize::bwareaopen,map
  COMPILE_OPT idl2

  ;显示原始二值图像
  ;i=image(BYTSCL(map), DIMENSIONS=[1200,400], LAYOUT=[3,1,1], $
  ;  TITLE='原始二值图像', WINDOW_TITLE='bwareaopen')

  ;进行分块标注
  region = LABEL_REGION(map, /ULONG, /ALL_NEIGHBORS)
  ;对标注进行直方图统计
  hist = HISTOGRAM(region, BINSIZE=1, LOCATIONS=locs, min=1)

  ;找到最大的块，并显示（这里没有考虑有多个最大块的情况）
  result = region EQ ((WHERE(hist EQ MAX(hist)))[0]+1)
  ;i=image(BYTSCL(result), /CURRENT, LAYOUT=[3,1,2], TITLE='最大的块')

  ;找到像元数大于1000的块，并显示
  label = locs[WHERE(hist GE 800)]
  sz = SIZE(map, /DIMENSIONS)
  result = BYTARR(sz[0], sz[1])
  FOR i=0,N_ELEMENTS(label)-1 DO BEGIN
    result = result OR (region EQ label[i])
  ENDFOR
  ;i=image(BYTSCL(result), /CURRENT, LAYOUT=[3,1,3], TITLE='大于1000的块')
  RETURN,result
END
;******************************************************************
;todo:实现与 Matlab 的 confusionmat 相同的功能
;混淆矩阵
;******************************************************************
FUNCTION raftrecognize::confmat,Y,T

  ysize = SIZE(Y,/dimension);
  n = ysize[0]
  c = ysize[1]
  Tsize=SIZE(T,/dimension);
  n2 = Tsize[0]
  c2 = Tsize[1]



  ;% Assume two classes with 0-1 encoding
  c = 2;
  class1 = WHERE(T GT 0.5);
  TL = MAKE_ARRAY(n, 1,value=1,/double);
  TL(class1) = 2;
  class2 = WHERE(Y GT 0.5);
  Yclass = MAKE_ARRAY(n, 1,VALUE=1,/double);
  Yclass(class2) = 2;

  ;% Compute
  correct = (Yclass EQ TL);
  sum=TOTAL(TOTAL(correct));
  rate=[sum*100/n,sum];

  Carray=MAKE_ARRAY(c,c,VALUE=0,/double);
  FOR i=0,c-1 DO BEGIN
    FOR j=0,c-1 DO BEGIN
      ci = DOUBLE(i+1)
      cj = DOUBLE(j+1)
      Carray[i,j] = TOTAL((Yclass EQ cj)*(TL EQ ci));
    ENDFOR
  ENDFOR
  RETURN,Carray
END
;参数设置
PRO raftrecognize::SetProperty, mouseType = mouseType,tmpfile = tmpfile,filedir = filedir
  self.MOUSETYPE= mouseType
  self.tmpfile = tmpfile
  self.filedir = filedir
END

;TODO 鼠标滚轮时的事件
PRO raftrecognize::WheelEvents,wType,xPos,yPos
  COMPILE_OPT idl2

  self.OWINDOW.GETPROPERTY, dimensions = winDims,graphics_tree = oView
  oView.GETPROPERTY, viewPlane_Rect = viewRect

  IF wType GT 0 THEN rate = 0.8 ELSE rate = 1.125


  oriDis =[xPos,yPos]*viewRect[2:3]/winDims
  viewRect[0:1]+=(1-rate)*oriDis
  viewRect[2:3]= viewRect[2:3]*rate
  ;
  oView.SETPROPERTY, viewPlane_Rect = viewRect
  self.OWINDOW.Draw
END

;todo:根据两个点返回矩形的四个点坐标
FUNCTION raftrecognize::CalRectPoints,ulPos,drPos

  self.OWINDOW.GETPROPERTY, dimensions = winDims,graphics_tree = oView
  oView.GETPROPERTY, viewPlane_Rect = viewRect
  lLoc = viewRect[0:1]+[ulPos[0],ulPos[1]]*viewRect[2:3]/winDims
  rLoc = viewRect[0:1]+[drPos[0],drPos[1]]*viewRect[2:3]/winDims
  IF ABS(rLoc[0]-lLoc[0]) EQ 0 THEN rLoc[0]= lLoc[0]+1
  IF ABS(rLoc[1]-lLoc[1]) EQ 0 THEN rLoc[1]= lLoc[1]+1
  RETURN,[[lLoc],[lLoc[0],rLoc[1]],$
    [rLoc],[rLoc[0],lLoc[1]]]
END
;鼠标点击时的事件
PRO raftrecognize::MousePress,xpos,ypos
  COMPILE_OPT idl2
  self.MOUSELOC[0:1] = [xPos,yPos]
  ;*(self.initlinedata) = [*(self.initlinedata),self.MOUSELOC[0:1]]
  CASE self.MOUSETYPE OF
    ;放大
    2: BEGIN
      data = self.CALRECTPOINTS(self.MOUSELOC[0:1],self.MOUSELOC[0:1])
      self.ORECT.SETPROPERTY, hide =0,data = data
    END
    ;缩小
    3: BEGIN
      data = self.CALRECTPOINTS(self.MOUSELOC[0:1],self.MOUSELOC[0:1])
      ;void = dialog_Message(string(data),/infor)
      self.ORECT.SETPROPERTY, hide =0,data = data
    END
    ;画折线
    4: BEGIN
      data = self.initlinedata
      self.INITLINE.SETPROPERTY, data = data
      self.OWINDOW.Draw
    END
    ELSE:
  ENDCASE
END


;鼠标弹起时的操作
PRO raftrecognize::MouseRelease,xpos,ypos
  COMPILE_OPT idl2

  self.ORECT.SETPROPERTY, hide =1
  curLoc = [xPos,yPos]
  self.OWINDOW.GETPROPERTY, dimensions = winDims,graphics_tree = oView
  oView.GETPROPERTY, viewPlane_Rect = viewRect

  CASE self.MOUSETYPE OF
    ;放大
    2: BEGIN
      data = self.CALRECTPOINTS(self.MOUSELOC,curLoc)
      maxV = MAX(data[0,*],min= minV)
      xRange = [minV,maxV]

      maxV = MAX(data[1,*],min= minV)
      yRange = [minV,maxV]
      ;
      viewRate= viewRect[3]/viewRect[2]
      rectRate = (yRange[1]-yRange[0])/(xRange[1]-xRange[0])
      ;
      IF viewRate GT rectRate THEN BEGIN
        width = xRange[1]-xRange[0]
        height = (xRange[1]-xRange[0])*winDims[1]/winDims[0]
        viewStartLoc= [TOTAL(xRange),TOTAL(yRange)]/2-[width,height]/2
      ENDIF ELSE BEGIN
        width = (yRange[1]-yRange[0])*winDims[0]/winDims[1]
        height = yRange[1]-yRange[0]
        viewStartLoc= [TOTAL(xRange),TOTAL(yRange)]/2-[width,height]/2
      ENDELSE
      newVP =[viewStartLoc,width,height ]
      oView.SETPROPERTY, viewPlane_Rect = newVP
    END
    ;缩小
    3: BEGIN
      data = self.CALRECTPOINTS(self.MOUSELOC,curLoc)
      maxV = MAX(data[0,*],min= minV)
      xRange = [minV,maxV]

      maxV = MAX(data[1,*],min= minV)
      yRange = [minV,maxV]
      ;
      viewRate= viewRect[3]/viewRect[2]
      rectRate = (yRange[1]-yRange[0])/(xRange[1]-xRange[0])
      ;
      IF viewRate GT rectRate THEN BEGIN
        viewRect[2:3] = viewRect[2:3]*(viewRect[2])/(2*(xRange[1]-xRange[0]))
        ViewRect[0:1]= [TOTAL(xRange),TOTAL(yRange)]/2 - viewRect[2:3]/2

      ENDIF ELSE BEGIN
        viewRect[2:3] = viewRect[2:3]*(viewRect[3])/(2*(yRange[1]-yRange[0]))
        ViewRect[0:1]= [TOTAL(xRange),TOTAL(yRange)]/2 - viewRect[2:3]/2

      ENDELSE
      oView.SETPROPERTY, viewPlane_Rect = ViewRect

    END
    ;画折线
    4: BEGIN
      data = self.CALRECTPOINTS(self.MOUSELOC,curLoc)
      self.INITLINE.SETPROPERTY, data = data
      self.OWINDOW.Draw
    END
    ELSE:
  ENDCASE

  self.OWINDOW.Draw
END


;鼠标双击时的操作
PRO raftrecognize::DbClick,drawId,button,xpos,ypos
  self.ORIGINALSHOW
END

PRO raftrecognize::MouseMotion,xpos,ypos
  ;
  curLoc = [xPos,yPos]
  ;
  self.OWINDOW.GETPROPERTY, dimensions = winDims,graphics_tree = oView
  oView.GETPROPERTY, viewPlane_Rect = viewRect

  CASE self.MOUSETYPE OF
    ;平移
    1: BEGIN
      ;屏幕偏移量
      offset = curLoc- self.MOUSELOC
      ;对应偏移量
      viewRect[0:1]-=offset*viewRect[2:3]/WinDims
      oView.SETPROPERTY, viewPlane_Rect = viewRect
      self.OWINDOW.Draw
      ;
      self.MOUSELOC = curLoc
      self.OWINDOW.SETCURRENTCURSOR, 'Move'
    END
    ;放大
    2: BEGIN
      data = self.CALRECTPOINTS(self.MOUSELOC,curLoc)
      self.ORECT.SETPROPERTY, data = data
      self.OWINDOW.Draw
    END
    ;缩小
    3: BEGIN
      data = self.CALRECTPOINTS(self.MOUSELOC,curLoc)
      self.ORECT.SETPROPERTY, data = data
      self.OWINDOW.Draw
    END
    ELSE:
  ENDCASE
END
;todo初始化图像显示，注意XY方向同比例变换
PRO raftrecognize::originalShow
  ;
  self.OWINDOW.GETPROPERTY, dimensions = windowDims,graphics_tree = oView
  imageDims = self.IMAGEDIMS
  ;
  imgRate = FLOAT(imageDims[0])/imageDims[1]
  viewRate = FLOAT(windowDims[0])/windowDims[1]
  ;
  IF imgRate GT viewRate THEN BEGIN
    viewWidth = imageDims[0]
    viewHeight = imageDims[0]/viewRate
    viewPlant_rect = [0, -(viewHeight-imageDims[1])/2,viewWidth,viewHeight]

  ENDIF ELSE BEGIN
    viewHeight = imageDims[1]
    viewwidth = imageDims[1]*ViewRate
    viewPlant_rect = [-(viewwidth-imageDims[0])/2,0,viewWidth,viewHeight]

  ENDELSE
  oView.SETPROPERTY, viewPlane_Rect = viewPlant_rect,dimensions = WindowDims
  self.OWINDOW.draw
END
;todo构建显示图像体系
PRO raftrecognize::CreateDrawImage
  oView = OBJ_NEW('IDLgrView',color = [255,255,255])
  self.OWINDOW.SETPROPERTY, graphics_tree = oView

  queryStatus = QUERY_TIFF(self.INFILE, imageInfo)
  IF queryStatus EQ 0 THEN BEGIN
    self.INITFLAG= 0
    RETURN
  ENDIF

  self.IMAGEDIMS = imageInfo.DIMENSIONS
  data =READ_TIFF(self.INFILE,GEOTIFF=GeoKeys)
  ;data =READ_TIFF(self.INFILE,GEOTIFF=GeoKeys)
  TileData = BYTSCL(data)
  self.ORIDATA = PTR_NEW(TileData,/no_Copy)
  ;initialLSF=2*MAKE_ARRAY(290,280,/double,VALUE=1)
  ;initialLSF[0:289,0:279]=-2
  ;self.CONDATA = PTR_NEW(TRANSPOSE(initialLSF),/no_Copy)

  ;
  IF imageInfo.CHANNELS EQ 1 THEN BEGIN
    ;
    self.RGBTYPE =0
    self.OIMAGE = OBJ_NEW('IDLgrImage',*(self.ORIDATA) )

  ENDIF ELSE BEGIN
    self.RGBTYPE =1
    self.OIMAGE = OBJ_NEW('IDLgrImage',*(self.ORIDATA) ,INTERLEAVE =0)
  ENDELSE
  ;辅助红色矩形，初始化为隐藏
  self.ORECT = OBJ_NEW('IDLgrPolygon', $
    style =1,$
    thick=1,$
    color = [230,0,0])
  
  ;oTopModel = OBJ_NEW('IDLgrModel')
  oModel = OBJ_NEW('IDLgrModel')
  ;lModel = OBJ_NEW('IDLgrModel')
  oModel.ADD, [self.OIMAGE,self.ORECT]
  ;lModel.add,self.OPRELINE
  ;  lModel->Rotate,[1,0,0],180
  ;  lModel->Rotate,[0,0,1],180
  ;oTopModel.ADD,oModel
  oView.Add,oModel
  self.ORIGINALSHOW
  self.INITFLAG= 1
END
FUNCTION raftrecognize::INIT,infile,drawID
  ;
  self.INFILE = infile
   ;传入的drawID
  WIDGET_CONTROL, drawID,GET_VALUE = oWindow
  ;Widget_Info(drawID,GET_ = oWindow)
  self.OWINDOW = oWindow
  ;self.OWINDOW = OBJ_NEW("IDLgrWindow")
  ;调用CreateImage方法创建显示图像
  self.CREATEDRAWIMAGE
  RETURN, self.INITFLAG
END

;对象类定义
PRO raftrecognize__define
  struct = {raftrecognize, $
    initFlag  : 0b, $
    mouseType : 0B, $;鼠标状态，0-默认,1-平移,2-放大,3-缩小。
    mouseLoc : FLTARR(2), $ ;

    infile: '' , $
    hdffile: '' , $
    tmpfile: '' , $
    filedir: '' , $
    hbasename: '' , $
    rgbType : 0, $
    imageDims : LONARR(2), $
    oriData  : PTR_NEW(), $;图像数据
    ;conData  : PTR_NEW(), $;生成的海岸线数据
    oWindow : OBJ_NEW(), $;显示窗口
    oImage  : OBJ_NEW(), $;显示图像的对象
    oRect   : OBJ_NEW(), $;放大缩小矩形对象
    ;oContour   : OBJ_NEW(), $;生成海岸线对象
    DrawID: 0L $
  }
END

