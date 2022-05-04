 ; Copyright (C) 2016 by The HDF Group.
;   All rights reserved.
;
;   This example code illustrates how to access and visualize GES DISC AIRS v6
; Swath HDF-EOS2 file in NCL. 
;
;  If you have any questions, suggestions, comments on this example, please 
; use the HDF-EOS Forum (http://hdfeos.org/forums). 
;
;  If you would like to see an  example of any other NASA HDF/HDF-EOS data 
; product that is not listed in the HDF-EOS Comprehensive Examples page 
; (http://hdfeos.org/zoo), feel free to contact us at eoshelp@hdfgroup.org or 
; post it at the HDF-EOS Forum (http://hdfeos.org/forums).
; Usage:save this script and run  
; $ncl AIRS.2002.08.30.227.L2.RetStd_H.v6.0.12.0.G14101125810.hdf.ncl
; Tested under: NCL 6.3.0
; Last updated: 2016-11-22
load "$NCARG_ROOT/lib/ncarg/nclex/gsun/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; This is an example of a 2-D data field.
; We assume that users know how to obtain data field information from HDFView.
; For information about HDFView, visit [1].
begin
      ;configPath="$/home/sw/Download/NCLlearn/work/sat"
     ; configName="awx.config"
  file_name = "/home/sw/Download/NCLlearn/1work/sat/dat/1-130_41km_ 649.61-681.99cm-1/AIRS.2013.07.12.172.L1B.SUBIBRAD.v5.0.22.0.G13311152621.hdf"
; To read HDF-EOS2 files, .he2 is appended to the argument.
; For more information, consult section 4.3.2 of [2].
  eos_file=addfile(file_name+".he2", "r") ; Read file.
  ra=eos_file->radiances_L1B_AIRS_Science(:,:,74)

;ra=eos_file->radiances_L1B_AIRS_Science(:,:,channel15)
; NCL. For more information, visit section 4.3.2 of [2].
; Associate longitude and latitude.
  lat2d=eos_file->Latitude_L1B_AIRS_Science(:,:)
  lon2d=eos_file->Longitude_L1B_AIRS_Science(:,:)
  ang=eos_file->scanang_L1B_AIRS_Science(:,:)
; See AIRS L1 user guide [3] for unit specification.
; all weights are set to one 
  c1 = 1.191066e-5
  c2 = 1.438833
  la1 = 667.782
  ;la1 = (10^4)/la
; la2 = 8.1;
; la2 = 1e4/la2;
  tb_c= (c2*la1)/(log(1+(c1*la1^3/ra)))
  n=5
  c = lspoly(ang,tb_c, 1, n)    ; all weights are set to on
  dim=dimsizes(tb_c);weishu1
  mean1=new((/dim(0),dim(1)/),"float")
  do i=0,dim(0)-1
	do j=0,dim(1)-1
	mean1(i,j)=c(i,0) + c(i,1)*ang(i,j) + c(i,2)*ang(i,j)^2+ c(i,3)*ang(i,j)^3+c(i,4)*ang(i,j)^4
	;print(mean1)
	disturb=tb_c-mean1
	end do
  end do  
   
;创建工作台
wks_type = "png"
;wks_type@wkWidth                                = 1000          ;工作台宽度
;wks_type@wkHeight                               = 1000        ;工作台高度
xwks = gsn_open_wks(wks_type,"disturbchannel15")       


;***********************************设置数据属性****************************
        lat2d@units               = "degrees_north"
        lon2d@units               = "degrees_east"
        disturb@lat2d                 = lat2d
        disturb@lon2d                 = lon2d
        disturb@units                 = "K"
        ;disturb@long_name             = ""
        disturb@coordinates           = "lat2d lon2d"
        ;minvalue=min(disturb)
        ;maxvalue=max(disturb)
        ;disturb@_FillValue            = -999. 
;************************************设置配置文件****************************
        ;configFile= configPath+"/"+configName
        ;等经纬投影经纬度配置
        ;truelat1 = stringtofloat(readParamInfo(configFile,"TRUELAT1"))
        ;truelat2 = stringtofloat(readParamInfo(configFile,"TRUELAT2"))
        ;stdlon   = stringtofloat(readParamInfo(configFile,"STDLON"))
        ;print(truelat1)
; *************************************设置工作台****************************
        res                                         = True
        ;res@gsnDraw                                 = False ;不再单独创建每个res的frame,写成false没有图
        res@gsnFrame                                = False ;不再单独对一个res　draw（与上句一起常用在overlay）
        res@gsnMaximize                             = True
        res@gsnPaperOrientation                     = "portrait"
        res@gsnAddCyclic                            = False                                 
        res@pmTickMarkDisplayMode                   = "Always"
;*************************************网格设置****************************
  res@cnFillOn=True;打开存在的等值线填充颜色
  res@cnLinesOn=False  ;关掉等值线填充
  res@cnLineLabelsOn=False  ;关掉等值线说明标签
  res@gsnSpreadColors=False
  ;设置地图经纬大小
  res@mpGridAndLimbOn                     = False
;通过corners方法正确的限定地图的范围  
  res@mpLimitMode                         = "LatLon"
  res@mpMinLatF=0.
  res@mpMaxLatF=50.
  res@mpMinLonF=100.
  res@mpMaxLonF=150.

  res@trGridType                              = "TriangularMesh"
  res@cnFillOn                                = True;打开存在的等值线填充颜色
  res@cnFillMode                              = "RasterFill"
  res@cnLinesOn                               = False;关掉等值线填充
  res@cnLineLabelsOn                          = False ;关掉等值线说明标签 

  ; 坐标轴属性
  ;res@tmXBOn                                     = True   ;关掉坐标刻度隐藏
  ;res@tmXTOn                                     = False ;上部刻度隐藏
  ;res@tmYLOn                                     = True
  ;res@tmYROn                                     = False ;右部分刻度隐藏
  ;res@tmXBBorderOn                            = True    ;四周边框
  ;res@tmXTBorderOn                            = True
  ;res@tmYLBorderOn                            = True
  ;res@tmYRBorderOn                            = True

;色标显示
  res@lbLabelBarOn                            = True
  res@lbLabelAutoStride                       = True
; 垂直色标条
  res@lbOrientation                           = "Vertical";lb的方向
;自定义等值线
; res@cnLevelSelectionMode                    = "ExplicitLevels" 
 
;*************************************绘制设置****************************
        ;图片输出目录
        ;outputPath="fy2e"
        ;checkOrCreateDir(outputPath)
        ;产品文件名

        ;gsn_define_colormap(xwks,"amwg")
; *************************************色标设置****************************  
 ;---Using the color map to define colors for contours or vectors
 
res@cnFillPalette = "BlueDarkRed18"

;---Reading the colormap into an N x 4 (RGBA) array
  cmap = read_colormap_file("BlueDarkRed18")
;控制色标开始值与结束值
res@gsnSpreadColorStart =min(disturb)
res@gsnSpreadColorEnd   =max(disturb)

;---Using the color map to define a color map for the workstation
  gsn_define_colormap(xwks,"BlueDarkRed18") 
        res@tiMainString="2013.07.13 41Km"
       ;ra@long_name=data@hdfeos_name      
;叠加shp底图
        mapres=True
       ; addEastChinaShpMap(mapres)
        ;chinamap = add_china_map(wks,contour,mapres)
        ;draw(contour)
        contour= gsn_csm_contour_map(xwks,disturb,res)
        frame(xwks)
        ;draw(contour)
        delete([/xwks,res,contour/])  
end
; References
;
; [1] http://www.hdfgroup.org/hdf-java-html/hdfview/
; [2] http://hdfeos.org/software/ncl.php
; [3] http://disc.sci.gsfc.nasa.gov/AIRS/documentation/v6_docs/v6releasedocs-1/V6_L2_Product_User_Guide.pdf
