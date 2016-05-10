;+
; :Description:
;    实现与 Matlab 的 bwareaopen 相同的功能
;    将 bindata.sav 与 bwareaopen.pro 放在同一路径下即可运行
;
; :Author: duhj@esrichina.com.cn
;-
function bwareaopen,map
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
  return,result
END