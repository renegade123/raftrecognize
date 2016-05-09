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