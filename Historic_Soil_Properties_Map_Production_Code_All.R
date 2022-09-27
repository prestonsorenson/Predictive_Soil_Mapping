library(raster)
library(ranger)
library(svMisc)
normalize <- function(x) {
  return ((x - min(na.omit(x))) / (max(na.omit(x)) - min(na.omit(x))))
}

load('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/Analysis_2022_04_27.RData')

################Production Models###################
#############Prepare Data#############
soc_all_prod=soc_all[,carbon_features]
soc_all_prod[is.na(soc_all_prod)] <- 0
soc_all_prod=do.call(data.frame, lapply(soc_all_prod, function(x) replace(x, is.infinite(x), 0)))

tn_all_prod=tn_all[,tn_features]
tn_all_prod=tn_all_prod[tn_all_prod$tn<0.5,]
tn_all_prod[is.na(tn_all_prod)] <- 0
tn_all_prod=do.call(data.frame, lapply(tn_all_prod, function(x) replace(x, is.infinite(x), 0)))

cec_all_prod=cec_all[,cec_features]
cec_all_prod[is.na(cec_all_prod)] <- 0
cec_all_prod=do.call(data.frame, lapply(cec_all_prod, function(x) replace(x, is.infinite(x), 0)))
  
ec_all_prod=ec_all[,ec_all_features]
ec_all_prod[is.na(ec_all_prod)] <- 0
ec_all_prod=do.call(data.frame, lapply(ec_all_prod, function(x) replace(x, is.infinite(x), 0)))

ioc_all_prod=ioc_all[,ioc_all_features]
ioc_all_prod[is.na(ioc_all_prod)] <- 0
ioc_all_prod=do.call(data.frame, lapply(ioc_all_prod, function(x) replace(x, is.infinite(x), 0)))

clay_all_prod=clay_all[,clay_all_features]
clay_all_prod[is.na(clay_all_prod)] <- 0
clay_all_prod=do.call(data.frame, lapply(clay_all_prod, function(x) replace(x, is.infinite(x), 0)))

silt_all_prod=silt_all[,silt_all_features]
silt_all_prod[is.na(silt_all_prod)] <- 0
silt_all_prod=do.call(data.frame, lapply(silt_all_prod, function(x) replace(x, is.infinite(x), 0)))

sand_all_prod=sand_all[,sand_all_features]
sand_all_prod[is.na(sand_all_prod)] <- 0
sand_all_prod=do.call(data.frame, lapply(sand_all_prod, function(x) replace(x, is.infinite(x), 0)))

hzn_all_prod=hzn_all[,hzn_all_features]
hzn_all_prod[is.na(hzn_all_prod)] <- 0
hzn_all_prod=do.call(data.frame, lapply(hzn_all_prod, function(x) replace(x, is.infinite(x), 0)))


stock_all_prod=soc_stock_all[,stock__features]
stock_all_prod[is.na(stock_all_prod)] <- 0
stock_all_prod=do.call(data.frame, lapply(stock_all_prod, function(x) replace(x, is.infinite(x), 0)))


model_cec_all_prod=ranger(cec~., data=cec_all_prod, importance='impurity', quantreg = TRUE)
model_clay_all_prod=ranger(clay~., data=clay_all_prod,importance='impurity', quantreg = TRUE)
model_ec_all_prod=ranger(ec~.,data=ec_all_prod,importance='impurity', quantreg = TRUE)
model_ioc_all_prod=ranger(ioc~., data=ioc_all_prod,importance='impurity', quantreg = TRUE)
model_hzn_all_prod=ranger(hzn_thickness~., data=hzn_all_prod,importance='impurity', quantreg = TRUE)
model_sand_all_prod=ranger(sand~., data=sand_all_prod,importance='impurity', quantreg = TRUE)
model_silt_all_prod=ranger(silt~., silt_all_prod,importance='impurity', quantreg = TRUE)
model_soc_all_prod=ranger(soc~., data=soc_all_prod,importance='impurity', quantreg = TRUE)
model_tn_all_prod=ranger(tn~., data=tn_all_prod,importance='impurity', quantreg = TRUE)
model_stock_all_prod=ranger(soc_stock~., data=stock_all_prod,importance='impurity', quantreg = TRUE)


models=list(model_cec_all_prod, model_clay_all_prod, model_ec_all_prod, model_ioc_all_prod, model_hzn_all_prod, model_sand_all_prod, model_soc_all_prod, model_tn_all_prod, model_stock_all_prod)


###################Mapping#######################
bare_soil=stack('/media/preston/My Book/Saskatchewan/sk_bare_soil/sk_ls5_bare_soil_ndvi3_ndsi0_nbr1_focal10_filt/sk_ls5_bare_soil_ndsi_ndvi3_ndsi0_nbr1_focal10.tif')
b1=bare_soil[[1]]
b2=bare_soil[[2]]
b3=bare_soil[[3]]
b4=bare_soil[[4]]
b5=bare_soil[[5]]
b6=bare_soil[[6]]
b7=bare_soil[[7]]
ari=raster('/media/preston/My Book/Saskatchewan/sk_l5_ARI_median_focal10_100m/sk_l5_ARI_median_focal10_100m.tif')
ari_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_ARI__noBareSoil_median_focal10_100m/sk_l5_ARI_noBareSoil_median_focal10_100m.tif')
CRSI=raster('/media/preston/My Book/Saskatchewan/sk_l5_CSRI_median_focal10_100m/sk_l5_CSRI_median_focal10_100m.tif')
CRSI_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_CSRI_noBareS0il_median_focal10_100m/sk_l5_CSRI_noBareS0il_median_focal10_100m.tif')
NDVI_July_Aug=raster('/media/preston/My Book/Saskatchewan/sk_l5_ndvi_median_july_aug_focal10_100m/sk_l5_ndvi_median_july_aug_focal10_100m.tif')
NDVI_July_Aug_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_ndvi_noBareSoil_median_july_aug_focal10_100m/sk_l5_ndvi_noBareSoil_median_july_aug_focal10_100m.tif')
NDVI_Sept_Oct=raster('/media/preston/My Book/Saskatchewan/sk_l5_ndvi_median_sept_oct_focal10_100m/sk_l5_ndvi_median_sept_oct_focal10_100m.tif')
NDVI_Sept_Oct_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_ndvi_noBareSoil_median_sept_oct_focal10_100m/sk_l5_ndvi_noBareSoil_median_sept_oct_focal10_100m.tif')
NDVI_SD=raster('/media/preston/My Book/Saskatchewan/sk_l5_ndvi_sd_focal10_100m/sk_l5_ndvi_sd_focal10_100m.tif')
SAVI_Jul_Aug=raster('/media/preston/My Book/Saskatchewan/sk_l5_SAVI_median_july_aug_focal10_100m/sk_l5_SAVI_median_july_aug_focal10_100m.tif')
SAVI_Jul_Aug_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_SAVI_noBareSoil_median_july_aug_focal10_100m/sk_l5_SAVI_noBareSoil_median_july_aug_focal10_100m.tif')
SAVI_Sept_Oct=raster('/media/preston/My Book/Saskatchewan/sk_l5_SAVI_median_sept_oct_focal10_100m/sk_l5_SAVI_median_sept_oct_focal10_100m.tif')
SAVI_Sept_Oct_noBareSoil=raster('/media/preston/My Book/Saskatchewan/sk_l5_SAVI_noBareSoil_median_sept_oct_focal10_100m/sk_l5_SAVI_noBareSoil_median_sept_oct_focal10_100m.tif')
dem_3x3_sd3x3=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_3x3_sd3x3.tif')
dem_3x3_sd5x5=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_3x3_sd5x5.tif')
dem_3x3_sd9x9=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_3x3_sd9x9.tif')
dem_3x3_sd21x21=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_3x3_sd21x21.tif')
dem_9x9_sd21x21=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_9x9_sd21x21.tif')
dem_9x9_sd101x101=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_9x9_sd101x101.tif')
dem_3x3_tri=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_3x3_tri.tif') 
dem_5x5_tri=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_5x5_tri.tif')
dem_9x9_tri=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_9x9_tri.tif')
dem_9x9_tri_20=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/100m/sk_elevation_alos_9x9_tri_20.tif')
precip=raster('/media/preston/My Book/Saskatchewan/Saskatchewan_Climate/sk_era5_tp_median_100m_bicubic.tif')
temperature=raster('/media/preston/My Book/Saskatchewan/Saskatchewan_Climate/sk_era5_2mt_median_100m_bicubic.tif') 
MidSlope_Pos_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/MidSlope_Pos_100m.tif')
Norm_Height_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/Norm_Height_100m.tif') 
Slope_Height_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/Slope_Height_100m.tif')
Stand_Height_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/Stand_Height_100m.tif') 
SWI_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/SWI_100m.tif')
Valley_Depth_100m=raster('/media/preston/My Book/Saskatchewan/sk_elevation_alos/Terrain_Derivatives/Valley_Depth_100m.tif')

NDWI=readAll(raster('/media/preston/My Book/Saskatchewan/sk_sen2_ndwi_med_july/sk_sen2_ndwi_med_july_100m.tif'))
NDWI@data@values[NDWI@data@values>0] <- 1

raster_stack=stack(b1, b2, b3, b4, b5, b6, b7, ari, ari_noBareSoil, CRSI, CRSI_noBareSoil, NDVI_July_Aug, NDVI_July_Aug_noBareSoil, NDVI_Sept_Oct, NDVI_Sept_Oct_noBareSoil, NDVI_SD, SAVI_Jul_Aug, SAVI_Jul_Aug_noBareSoil, SAVI_Sept_Oct, SAVI_Sept_Oct_noBareSoil, dem_3x3_sd3x3, dem_3x3_sd5x5, dem_3x3_sd9x9, dem_3x3_sd21x21, dem_9x9_sd21x21, dem_9x9_sd101x101, dem_3x3_tri, dem_5x5_tri, dem_9x9_tri, dem_9x9_tri_20, precip, temperature, MidSlope_Pos_100m, Norm_Height_100m, Slope_Height_100m, Stand_Height_100m, SWI_100m, Valley_Depth_100m)

#############predictions############################
setwd('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/Maps/3D')
soil_parameters=c('cec', 'clay', 'ec','ioc', 'hzn_thickness', 'sand', 'soc', 'tn_a', 'soc_stock')

tiles=shapefile("/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/analysis_tiles.shp")

pred_results=vector('list')
pred_results[[1]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), horizon=character(), x=integer(), y=integer())
pred_results[[2]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), horizon=character(), x=integer(), y=integer())
pred_results[[3]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), horizon=character(), x=integer(), y=integer())
pred_results[[4]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), horizon=character(), x=integer(), y=integer())
pred_results[[5]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), horizon=character(), x=integer(), y=integer())
pred_results[[6]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), horizon=character(), x=integer(), y=integer())
pred_results[[7]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), horizon=character(), x=integer(), y=integer())
pred_results[[8]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), horizon=character(), x=integer(), y=integer())
pred_results[[9]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), horizon=character(), x=integer(), y=integer())




for (i in 1:length(tiles)){
progress(i, max.value=length(tiles))
try({
tile_sub=tiles[tiles$id==i,]
NDWI_sub=crop(NDWI, tile_sub)
raster_sub=crop(raster_stack, tile_sub)
raster_sub=rasterToPoints(raster_sub)
xy=raster_sub[,1:2]
raster_sub=raster_sub[,-c(1:2)]
raster_sub=data.frame(raster_sub)
colnames(raster_sub)=c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'ari', 'ari_noBareSoil', 'CRSI', 'CRSI_noBareSoil', 'NDVI_July_Aug', 'NDVI_July_Aug_noBareSoil', 'NDVI_Sept_Oct', 'NDVI_Sept_Oct_noBareSoil', 'NDVI_SD', 'SAVI_Jul_Aug', 'SAVI_Jul_Aug_noBareSoil', 'SAVI_Sept_Oct', 'SAVI_Sept_Oct_noBareSoil', 'dem_3x3_sd3x3', 'dem_3x3_sd5x5', 'dem_3x3_sd9x9', 'dem_3x3_sd21x21', 'dem_9x9_sd21x21', 'dem_9x9_sd101x101', 'dem_3x3_tri', 'dem_5x5_tri', 'dem_9x9_tri','dem_9x9_tri_20', 'precip', 'temperature', 'MidSlope_Pos_100m', 'Norm_Height_100m', 'Slope_Height_100m', 'Stand_Height_100m', 'SWI_100m', 'Valley_Depth_100m')
raster_sub[is.na(raster_sub)] <- 0
raster_sub=do.call(data.frame, lapply(raster_sub, function(x) replace(x, is.infinite(x), 0)))
raster_sub$hzn='A'
raster_sub_b=raster_sub
raster_sub_b$hzn='B'
raster_sub_c=raster_sub
raster_sub_c$hzn='C'

raster_sub=rbind(raster_sub, raster_sub_b, raster_sub_c)

#soc stock
model_sub=models[[9]]
pred_sub=predict(model_sub, data=raster_sub[raster_sub$hzn=='A',], type='quantiles', quantiles=c(0.25, 0.5, 0.75))
pred_sub=pred_sub$predictions
pred_sub=data.frame(pred_sub)
colnames(pred_sub)=c("pred_25", "pred_50", "pred_75")
pred_sub$iqr=pred_sub$pred_75-pred_sub$pred_25
pred_sub=data.frame(pred_sub, xy)
colnames(pred_sub)=c("pred_25", "pred_50", "pred_75", 'iqr', 'x', 'y')

wd=paste('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/Maps/3D/soc_stock')
setwd(wd)
pred_results=pred_sub
pred_results_25=pred_results[,c(1,5:6)]
pred_results_50=pred_results[,c(2,5:6)]
pred_results_75=pred_results[,c(3,5:6)]
pred_results_iqr=pred_results[,c(4,5:6)]
pred_results_iqr_ratio=pred_results[,c(4,5:6)]
pred_results_iqr_ratio$iqr=normalize(pred_results_iqr_ratio$iqr)

coordinates(pred_results_25)=~x+y
coordinates(pred_results_50)=~x+y
coordinates(pred_results_75)=~x+y
coordinates(pred_results_iqr)=~x+y
coordinates(pred_results_iqr_ratio)=~x+y

pred_results_25=rasterFromXYZ(pred_results_25)
pred_results_50=rasterFromXYZ(pred_results_50)
pred_results_75=rasterFromXYZ(pred_results_75)
pred_results_iqr=rasterFromXYZ(pred_results_iqr)
pred_results_iqr_ratio=rasterFromXYZ(pred_results_iqr_ratio)

pred_results_25=mask(pred_results_25, NDWI_sub, maskvalue=1)
pred_results_50=mask(pred_results_50, NDWI_sub, maskvalue=1)
pred_results_75=mask(pred_results_75, NDWI_sub, maskvalue=1)
pred_results_iqr=mask(pred_results_iqr, NDWI_sub, maskvalue=1)
pred_results_iqr_ratio=mask(pred_results_iqr_ratio, NDWI_sub, maskvalue=1)

pred_stack=stack(pred_results_25, pred_results_50, pred_results_75, pred_results_iqr, pred_results_iqr_ratio)

names(pred_stack)=c("pred_results_25","pred_results_50","pred_results_75","pred_results_iqr","pred_results_iqr_ratio")

figure_name=paste("sk_predicted",soil_parameters[9], i, sep='_')
crs(pred_stack)=crs(ari)
writeRaster(pred_stack, figure_name, format="GTiff", overwrite=TRUE)


for (j in 1:(length(models)-1))
{
model_sub=models[[j]]
pred_sub=predict(model_sub, data=raster_sub, type='quantiles', quantiles=c(0.25, 0.5, 0.75))
pred_sub=pred_sub$predictions
pred_sub=data.frame(pred_sub)
colnames(pred_sub)=c("pred_25", "pred_50", "pred_75")
pred_sub$iqr=pred_sub$pred_75-pred_sub$pred_25
pred_sub=data.frame(pred_sub, raster_sub$hzn, xy)
colnames(pred_sub)=c("pred_25", "pred_50", "pred_75", 'iqr', 'horizon', 'x', 'y')

for (k in c('A', 'B', 'C')){
wd=paste('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/Maps/3D', soil_parameters[j], k, sep='/')
setwd(wd)
pred_results=pred_sub[pred_sub$horizon==k,]

pred_results_25=pred_results[,c(1,6:7)]
pred_results_50=pred_results[,c(2,6:7)]
pred_results_75=pred_results[,c(3,6:7)]
pred_results_iqr=pred_results[,c(4,6:7)]
pred_results_iqr_ratio=pred_results[,c(4,6:7)]
pred_results_iqr_ratio$iqr=normalize(pred_results_iqr_ratio$iqr)

coordinates(pred_results_25)=~x+y
coordinates(pred_results_50)=~x+y
coordinates(pred_results_75)=~x+y
coordinates(pred_results_iqr)=~x+y
coordinates(pred_results_iqr_ratio)=~x+y

pred_results_25=rasterFromXYZ(pred_results_25)
pred_results_50=rasterFromXYZ(pred_results_50)
pred_results_75=rasterFromXYZ(pred_results_75)
pred_results_iqr=rasterFromXYZ(pred_results_iqr)
pred_results_iqr_ratio=rasterFromXYZ(pred_results_iqr_ratio)

pred_results_25=mask(pred_results_25, NDWI_sub, maskvalue=1)
pred_results_50=mask(pred_results_50, NDWI_sub, maskvalue=1)
pred_results_75=mask(pred_results_75, NDWI_sub, maskvalue=1)
pred_results_iqr=mask(pred_results_iqr, NDWI_sub, maskvalue=1)
pred_results_iqr_ratio=mask(pred_results_iqr_ratio, NDWI_sub, maskvalue=1)

pred_stack=stack(pred_results_25, pred_results_50, pred_results_75, pred_results_iqr, pred_results_iqr_ratio)

names(pred_stack)=c("pred_results_25","pred_results_50","pred_results_75","pred_results_iqr","pred_results_iqr_ratio")

figure_name=paste("sk_predicted",soil_parameters[j],k, i, sep='_')
crs(pred_stack)=crs(ari)
writeRaster(pred_stack, figure_name, format="GTiff", overwrite=TRUE)
}
}
})
}


