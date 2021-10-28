library(raster)
library(ranger)
library(svMisc)
normalize <- function(x) {
  return ((x - min(na.omit(x))) / (max(na.omit(x)) - min(na.omit(x))))
}

load('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/Analysis_2021_10_04.RData')

################Production Models###################
#############Prepare Data#############
soc_a_prod=soc_a[,carbon_features]
soc_a_prod[is.na(soc_a_prod)] <- 0
soc_a_prod=do.call(data.frame, lapply(soc_a_prod, function(x) replace(x, is.infinite(x), 0)))

tn_a_prod=tn_a[,tn_features]
tn_a_prod=tn_a_prod[tn_a_prod$tn<0.5,]
tn_a_prod[is.na(tn_a_prod)] <- 0
tn_a_prod=do.call(data.frame, lapply(tn_a_prod, function(x) replace(x, is.infinite(x), 0)))

soc_stock_prod=soc_stock[,stock_subon_features]
soc_stock_prod[is.na(soc_stock_prod)] <- 0
soc_stock_prod=do.call(data.frame, lapply(soc_stock_prod, function(x) replace(x, is.infinite(x), 0)))
soc_stock_prod$soc_stock=soc_stock_prod$soc_stock*10

cec_a_prod=cec_a[,cec_features]
cec_a_prod[is.na(cec_a_prod)] <- 0
cec_a_prod=do.call(data.frame, lapply(cec_a_prod, function(x) replace(x, is.infinite(x), 0)))
  
ec_a_prod=ec_a[,ec_a_features]
ec_a_prod[is.na(ec_a_prod)] <- 0
ec_a_prod=do.call(data.frame, lapply(ec_a_prod, function(x) replace(x, is.infinite(x), 0)))

ec_solum_prod=ec_solum[,ec_solum_features]
ec_solum_prod[is.na(ec_solum_prod)] <- 0
ec_solum_prod=do.call(data.frame, lapply(ec_solum_prod, function(x) replace(x, is.infinite(x), 0)))

ec_profile_prod=ec_profile[,ec_profile_features]
ec_profile_prod[is.na(ec_profile_prod)] <- 0
ec_profile_prod=do.call(data.frame, lapply(ec_profile_prod, function(x) replace(x, is.infinite(x), 0)))

ioc_a_prod=ioc_a[,ioc_a_features]
ioc_a_prod[is.na(ioc_a_prod)] <- 0
ioc_a_prod=do.call(data.frame, lapply(ioc_a_prod, function(x) replace(x, is.infinite(x), 0)))

ioc_all_prod=ioc_all[,ioc_all_features]
ioc_all_prod[is.na(ioc_all_prod)] <- 0
ioc_all_prod=do.call(data.frame, lapply(ioc_all_prod, function(x) replace(x, is.infinite(x), 0)))

clay_a_prod=clay_a[,clay_a_features]
clay_a_prod[is.na(clay_a_prod)] <- 0
clay_a_prod=do.call(data.frame, lapply(clay_a_prod, function(x) replace(x, is.infinite(x), 0)))

silt_a_prod=silt_a[,silt_a_features]
silt_a_prod[is.na(silt_a_prod)] <- 0
silt_a_prod=do.call(data.frame, lapply(silt_a_prod, function(x) replace(x, is.infinite(x), 0)))

sand_a_prod=sand_a[,sand_a_features]
sand_a_prod[is.na(sand_a_prod)] <- 0
sand_a_prod=do.call(data.frame, lapply(sand_a_prod, function(x) replace(x, is.infinite(x), 0)))


clay_solum_prod=clay_solum[,clay_solum_features]
clay_solum_prod[is.na(clay_solum_prod)] <- 0
clay_solum_prod=do.call(data.frame, lapply(clay_solum_prod, function(x) replace(x, is.infinite(x), 0)))

silt_solum_prod=silt_solum[,silt_solum_features]
silt_solum_prod[is.na(silt_solum_prod)] <- 0
silt_solum_prod=do.call(data.frame, lapply(silt_solum_prod, function(x) replace(x, is.infinite(x), 0)))

sand_solum_prod=sand_solum[,sand_solum_features]
sand_solum_prod[is.na(sand_solum_prod)] <- 0
sand_solum_prod=do.call(data.frame, lapply(sand_solum_prod, function(x) replace(x, is.infinite(x), 0)))


clay_all_prod=clay_all[,clay_all_features]
clay_all_prod[is.na(clay_all_prod)] <- 0
clay_all_prod=do.call(data.frame, lapply(clay_all_prod, function(x) replace(x, is.infinite(x), 0)))

silt_all_prod=silt_all[,silt_all_features]
silt_all_prod[is.na(silt_all_prod)] <- 0
silt_all_prod=do.call(data.frame, lapply(silt_all_prod, function(x) replace(x, is.infinite(x), 0)))

sand_all_prod=sand_all[,sand_all_features]
sand_all_prod[is.na(sand_all_prod)] <- 0
sand_all_prod=do.call(data.frame, lapply(sand_all_prod, function(x) replace(x, is.infinite(x), 0)))


model_cec_a_prod=ranger(cec~., data=cec_a_prod, importance='impurity', quantreg = TRUE)
model_clay_a_prod=ranger(clay~., data=clay_a_prod,importance='impurity', quantreg = TRUE)
model_clay_all_prod=ranger(clay~.,data=clay_all_prod, importance='impurity', quantreg = TRUE)
model_clay_solum_prod=ranger(clay~., data=clay_solum_prod,importance='impurity', quantreg = TRUE)
model_ec_a_prod=ranger(ec~.,data=ec_a_prod,importance='impurity', quantreg = TRUE)
model_ec_profile_prod=ranger(ec~.,data=ec_profile_prod,importance='impurity', quantreg = TRUE)
model_ec_solum_prod=ranger(ec~., data=ec_solum_prod,importance='impurity', quantreg = TRUE)
model_ioc_a_prod=ranger(ioc~., data=ioc_a_prod,importance='impurity', quantreg = TRUE)
model_ioc_all_prod=ranger(ioc~., data=ioc_all_prod,importance='impurity', quantreg = TRUE)
model_sand_a_prod=ranger(sand~., data=sand_a_prod,importance='impurity', quantreg = TRUE)
model_sand_all_prod=ranger(sand~.,data=sand_all_prod,importance='impurity', quantreg = TRUE)
model_sand_solum_prod=ranger(sand~.,data=sand_solum_prod, importance='impurity', quantreg = TRUE)
model_silt_a_prod=ranger(silt~., silt_a_prod,importance='impurity', quantreg = TRUE)
model_silt_all_prod=ranger(silt~., silt_all_prod,importance='impurity', quantreg = TRUE)
model_silt_solum_prod=ranger(silt~., silt_solum_prod,importance='impurity', quantreg = TRUE)
model_soc_a_prod=ranger(soc~., data=soc_a_prod,importance='impurity', quantreg = TRUE)
model_soc_stock_prod=ranger(soc_stock~.,soc_stock_prod,importance='impurity', quantreg = TRUE)
model_tn_a_prod=ranger(tn~., data=tn_a_prod,importance='impurity', quantreg = TRUE)

models=list(model_cec_a_prod, model_clay_a_prod, model_clay_all_prod, model_clay_solum_prod, model_ec_profile_prod, model_ec_a_prod, model_ec_solum_prod,
            model_ioc_a_prod, model_ioc_all_prod, model_sand_a_prod, model_sand_all_prod, model_sand_solum_prod, model_silt_a_prod, model_silt_all_prod, model_silt_solum_prod,
            model_soc_a_prod, model_soc_stock_prod)

models=list(model_cec_a_prod, model_clay_a_prod, model_ec_profile_prod, model_ioc_all_prod, model_sand_a_prod, model_silt_a_prod, model_soc_a_prod, model_soc_stock_prod, model_tn_a_prod)


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
setwd('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/Maps')
soil_parameters=c('cec_a', 'clay_a', 'clay_all', 'clay_solum', 'ec_profile', 'ec_a','ec_solum', 'ioc_a', 'ioc_all', 'sand_a', 'sand_all',
                  'sand_solum', 'silt_a', 'silt_all', 'silt_solum', 'soc_a', 'soc_stock', 'tn_a')

soil_parameters=c('cec_a', 'clay_a', 'ec_profile','ioc_all', 'sand_a', 'silt_a', 'soc_a', 'soc_stock', 'tn_a')

models=list(model_cec_a_prod, model_clay_a_prod, model_ec_profile_prod, model_ioc_all_prod, model_sand_a_prod, model_silt_a_prod, model_soc_a_prod, model_soc_stock_prod, model_tn_a_prod)

tiles=shapefile("/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/analysis_tiles.shp")

pred_results=vector('list')
pred_results[[1]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), x=integer(), y=integer())
pred_results[[2]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), x=integer(), y=integer())
pred_results[[3]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), x=integer(), y=integer())
pred_results[[4]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), x=integer(), y=integer())
pred_results[[5]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), x=integer(), y=integer())
pred_results[[6]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), x=integer(), y=integer())
pred_results[[7]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), x=integer(), y=integer())
pred_results[[8]]=data.frame(pred_25=double(), pred_50=double(), pred_75=double(), iqr=double(), x=integer(), y=integer())

for (i in 1:length(tiles)){
progress(i, max.value=length(tiles))
tile_sub=tiles[tiles$id==i,]
raster_sub=crop(raster_stack, tile_sub)
raster_sub=rasterToPoints(raster_sub)
xy=raster_sub[,1:2]
raster_sub=raster_sub[,-c(1:2)]
raster_sub=data.frame(raster_sub)
colnames(raster_sub)=c('b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'ari', 'ari_noBareSoil', 'CRSI', 'CRSI_noBareSoil', 'NDVI_July_Aug', 'NDVI_July_Aug_noBareSoil', 'NDVI_Sept_Oct', 'NDVI_Sept_Oct_noBareSoil', 'NDVI_SD', 'SAVI_Jul_Aug', 'SAVI_Jul_Aug_noBareSoil', 'SAVI_Sept_Oct', 'SAVI_Sept_Oct_noBareSoil', 'dem_3x3_sd3x3', 'dem_3x3_sd5x5', 'dem_3x3_sd9x9', 'dem_3x3_sd21x21', 'dem_9x9_sd21x21', 'dem_9x9_sd101x101', 'dem_3x3_tri', 'dem_5x5_tri', 'dem_9x9_tri','dem_9x9_tri_20', 'precip', 'temperature', 'MidSlope_Pos_100m', 'Norm_Height_100m', 'Slope_Height_100m', 'Stand_Height_100m', 'SWI_100m', 'Valley_Depth_100m')
raster_sub[is.na(raster_sub)] <- 0
raster_sub=do.call(data.frame, lapply(raster_sub, function(x) replace(x, is.infinite(x), 0)))
for (j in 1:length(models))
{
model_sub=models[[j]]
pred_sub=predict(model_sub, data=raster_sub, type='quantiles', quantiles=c(0.25, 0.5, 0.75))
pred_sub=pred_sub$predictions
pred_sub=data.frame(pred_sub)
colnames(pred_sub)=c("pred_25", "pred_50", "pred_75")
pred_sub$iqr=pred_sub$pred_75-pred_sub$pred_25
pred_sub=data.frame(pred_sub, xy)
pred_results[[j]]=rbind(pred_results[[j]], pred_sub)
}
}

for (j in 1:length(pred_results)){
pred_results_25=pred_results[[j]][,c(1,5:6)]
pred_results_50=pred_results[[j]][,c(2,5:6)]
pred_results_75=pred_results[[j]][,c(3,5:6)]
pred_results_iqr=pred_results[[j]][,c(4,5:6)]
pred_results_iqr_ratio=pred_results[[j]][,c(4,5:6)]
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

pred_results_25=mask(pred_results_25, NDWI, maskvalue=1)
pred_results_50=mask(pred_results_50, NDWI, maskvalue=1)
pred_results_75=mask(pred_results_75, NDWI, maskvalue=1)
pred_results_iqr=mask(pred_results_iqr, NDWI, maskvalue=1)
pred_results_iqr_ratio=mask(pred_results_iqr_ratio, NDWI, maskvalue=1)

pred_stack=stack(pred_results_25, pred_results_50, pred_results_75, pred_results_iqr, pred_results_iqr_ratio)

names(pred_stack)=c("pred_results_25","pred_results_50","pred_results_75","pred_results_iqr","pred_results_iqr_ratio")

figure_name=paste("sk_predicted_",soil_parameters[j], sep='')
crs(pred_stack)=crs(ari)
writeRaster(pred_stack, figure_name, format="GTiff", overwrite=TRUE)
}



