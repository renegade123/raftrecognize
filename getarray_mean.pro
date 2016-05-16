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
  FOR c=0,(FLOOR(image_size/sz)-1)*sz,sz DO BEGIN
    FOR r=0,(FLOOR(image_size2/sz)-1)*sz,sz DO BEGIN
      totalsamples = totalsamples + 1
      FOR i=0,num_images-1 DO BEGIN
        X[i,totalsamples-1]=mean(REFORM(IMAGES[c:c+sz-1,r:r+sz-1,i],sz^2,1),/DOUBLE);
      ENDFOR
      labelreform = REFORM(Label[r:r+sz-1,c:c+sz-1],sz^2,1)
      Y[*,totalsamples-1]=ROUND(TOTAL(TOTAL(labelreform))/(sz^2));
    ENDFOR
  ENDFOR
  XPTR = PTR_NEW(X);
  YPTR = PTR_NEW(Y);
  RETURN,[XPTR,YPTR]
END