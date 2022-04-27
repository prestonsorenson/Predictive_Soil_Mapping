library(raster)
library(clhs)
library(ranger)
library(Metrics)
library(ggplot2)
library(DescTools)
library(prospectr)
library(caret)
library(Cubist)
library(brnn)

setwd('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Figures')
#feature selection script
psm_ranger_feature_selection <- function (x, y){
  #feature selection
  #select bare soil
  data_temp=y[,c(x,'b1', 'b2', 'b3', 'b4','b5','b7', 'precip', 'temperature')]
  model_temp_bare_soil=ranger(dependent.variable.name = x, data=data_temp, importance='impurity')
  model_temp_bare_soil
  sort(importance(model_temp_bare_soil), decreasing=TRUE)
  
  #forward feature selection
  bare_soil_var=names(sort(importance(model_temp_bare_soil), decreasing=TRUE))
  bare_soil_var=c(x, bare_soil_var)
  
  val=vector('list')
  bs_feature=vector('list')
  q=0
  for (i in 2:length(bare_soil_var)){
    q=q+1
    temp=data_temp[,bare_soil_var]
    model=ranger(dependent.variable.name = x, data=temp, importance='impurity')
    val=c(val, model$prediction.error)
    bs_feature[[q]]=bare_soil_var
    bare_soil_var=c(x, names(sort(importance(model), decreasing=TRUE))[-length(sort(importance(model), decreasing=TRUE))])
    rm(model)
  }
  val_bare_soil=unlist(val)
  which.min(val_bare_soil)
  
  bare_soil_var=unlist(bs_feature[which.min(val)])[-1]
  
  
  
  #band ratios
  data_temp=y[,c(x,'ari','ari_noBareSoil','CRSI','CRSI_noBareSoil','NDVI_July_Aug','NDVI_July_Aug_noBareSoil','NDVI_Sept_Oct','NDVI_Sept_Oct_noBareSoil','NDVI_SD','SAVI_Jul_Aug','SAVI_Jul_Aug_noBareSoil','precip','temperature')]
  model_temp_band_ratios=ranger(dependent.variable.name = x, data=data_temp, importance='impurity')
  model_temp_band_ratios
  sort(importance(model_temp_band_ratios), decreasing=TRUE)
  
  #forward feature selection
  band_ratios_var=names(sort(importance(model_temp_band_ratios), decreasing=TRUE))
  band_ratios_var=c(x, band_ratios_var)
  
  
  val=vector('list')
  band_ratios_feature=vector('list')
  q=0
  for (i in 2:length(band_ratios_var)){
    q=q+1
    temp=data_temp[,band_ratios_var]
    model=ranger(dependent.variable.name = x, data=temp, importance='impurity')
    val=c(val, model$prediction.error)
    band_ratios_feature[[q]]=band_ratios_var
    band_ratios_var=c(x, names(sort(importance(model), decreasing=TRUE))[-length(sort(importance(model), decreasing=TRUE))])
    rm(model)
  }
  val_band_ratios=unlist(val)
  which.min(val_band_ratios)
  
  band_ratios_var=unlist(band_ratios_feature[which.min(val)])[-1]
  
  
  #terrain attributes
  data_temp=y[,c(x,'precip', 'temperature','dem_3x3_sd3x3','dem_3x3_sd5x5','dem_3x3_sd9x9','dem_3x3_sd21x21','dem_9x9_sd21x21','dem_9x9_sd101x101','dem_3x3_tri','dem_5x5_tri','dem_9x9_tri','dem_9x9_tri_20', 'MidSlope_Pos_100m', 'Norm_Height_100m', 'Slope_Height_100m', 'Stand_Height_100m', 'SWI_100m', 'Valley_Depth_100m')]
  model_temp_terrain=ranger(dependent.variable.name = x, data=data_temp, importance='impurity')
  model_temp_terrain
  sort(importance(model_temp_terrain), decreasing=TRUE)
  
  #forward feature selection
  terrain_var=names(sort(importance(model_temp_terrain), decreasing=TRUE))
  terrain_var=c(x, terrain_var)
  
  val=vector('list')
  terrain_feature=vector('list')
  q=0
  for (i in 2:length(terrain_var)){
    q=q+1
    temp=data_temp[,terrain_var]
    model=ranger(dependent.variable.name = x, data=temp, importance='impurity')
    val=c(val, model$prediction.error)
    terrain_feature[[q]]=terrain_var
    terrain_var=c(x, names(sort(importance(model), decreasing=TRUE))[-length(sort(importance(model), decreasing=TRUE))])
    rm(model)
  }
  val_terrain=unlist(val)
  which.min(val_terrain)
  
  terrain_var=unlist(terrain_feature[which.min(val)])[-1]
  
  #final features
  temp_init_var=c(x, bare_soil_var, band_ratios_var, terrain_var)
  temp_init_var=temp_init_var[!duplicated(temp_init_var)]
  
  var_train_a_init=y[,temp_init_var]
  
  #remove correlated features
  var_cor=findCorrelation(cor(var_train_a_init[,-1]), cutoff=0.9)
  var_cor
  var_cor=var_cor+1
  
  var_train_a_init=var_train_a_init[,-var_cor]
  
  
  
  model_temp_final=ranger(dependent.variable.name = x, data=var_train_a_init, importance='impurity')
  model_temp_final
  sort(importance(model_temp_final), decreasing=TRUE)
  
  #forward feature selection
  final_var=names(sort(importance(model_temp_final), decreasing=TRUE))
  final_var=c(x, final_var)
  val=vector('list')
  final_feature=vector('list')
  q=0
  for (i in 2:length(final_var)){
    q=q+1
    temp=y[,final_var]
    model=ranger(dependent.variable.name = x, data=temp, importance='impurity')
    val=c(val, model$prediction.error)
    final_feature[[q]]=final_var
    final_var=c(x, names(sort(importance(model), decreasing=TRUE))[-length(sort(importance(model), decreasing=TRUE))])
    rm(model)
  }
  val_final=unlist(val)
  which.min(val_final)
  
  final_var=unlist(final_feature[which.min(val)])[-1]
  final_var=c(x, final_var)
  
  return(final_var)
}

load('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/Analysis_2021_10_03.RData')

#####################Bulk Density Pedotransfer Function################################
npdb=read.csv('/home/preston/OneDrive/Graduate Studies/Post_Doc/Shared/Data/Sask_Soil_Data/npdb/NPDB_Combined_Data.csv')
npdb_sub=c('BULK_DEN', 'CARB_ORG', 'T_CLAY', 'T_SILT', "T_SAND")
npdb_bd=npdb[,npdb_sub]
npdb_bd=na.omit(npdb_bd)
npdb_bd=npdb_bd[npdb_bd$BULK_DEN<2.6,]

npdb_train=kenStone(npdb_bd[,-1], k=1431)
npdb_test=npdb_bd[npdb_train$test,]
npdb_train=npdb_bd[npdb_train$model,]

#train test builds
model_bd=ranger(BULK_DEN~., data=npdb_train, importance='impurity')
#model_bd=cubist(x=npdb_train[,-1], y=npdb_train[,1], committees=10, neighbors=5)
write.csv(sort(importance(model_bd), decreasing=TRUE), "bd_feature.csv")


pred_bd=predict(model_bd, data=npdb_test)
pred_bd=pred_bd$predictions

#model performance
cor.test(pred_bd, npdb_test$BULK_DEN)
sd(npdb_test$BULK_DEN)/rmse(npdb_test$BULK_DEN, pred_bd)
summary(lm(pred_bd~npdb_test$BULK_DEN))$r.squared
rmse(npdb_test$BULK_DEN, pred_bd)
CCC(npdb_test$BULK_DEN, pred_bd)$rho.c
bias(npdb_test$BULK_DEN, pred_bd)

#create plot
plot_data_bd=data.frame(npdb_test$BULK_DEN, pred_bd)
colnames(plot_data_bd)=c("actual", 'predicted')
bd_plot=ggplot(plot_data_bd, (aes(x=actual, y=predicted))) + geom_point(size=2) + xlim(0, 2.65) + ylim(0, 2.65) + geom_abline() + xlab(expression(paste('Observed Bulk Densityy (g cm'^"-3",")"))) + ylab(expression(paste('Predicted Cation Exchange Capacity (g cm'^"-3",")")))
bd_plot = bd_plot + annotate("text", x=0.1, y = 2.65, label=expression(paste("R"^"2 ", "= 0.27"))) + annotate("text", x=0.1, y=2.55, label=expression(paste("RMSE = 0.22 g cm"^"-3"))) + annotate("text", x=0.1, y=2.45, label=expression(paste(rho['c'], "= 0.47"))) + annotate("text", x=0.1, y=2.35, label='Bias = 0.01')
bd_plot

ggsave('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Figures/bulk_denisty.png', plot=last_plot(), width=11, height=8.5)


pred_null=rep(mean(npdb_train$BULK_DEN), length(npdb_test$BULK_DEN))
rmse(npdb_test$BULK_DEN, pred_null)



###############Data Preparation################
soil_points=shapefile('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Data/Sask_NPDB_Chemical_Physical_Data.shp')
raster_values=read.csv('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Data/Historical_Soil_Properties_Training_Data_2021_10_01.csv')

#predict bulk density values
bd_input=soil_points[,c(29,41,42,47)]
bd_input=apply(bd_input, 2, function(x) as.numeric(x))
bd_input=data.frame(bd_input)
colnames(bd_input)=colnames(npdb_bd)[-1]

bd_input[is.na(bd_input)] <- 999

bd=predict(model_bd, data=bd_input)
bd=bd$predictions

bd[bd_input$CARB_ORG==999] <- NA
bd[bd_input$T_CLAY==999] <- NA
bd[bd_input$T_SILT==999] <- NA
bd[bd_input$T_SAND==999] <- NA

soil_points$bd=bd

#Calculate Weighted Averages for Soil Properties
soil_points$thickness=abs(as.numeric(soil_points$Lower_Dept)-as.numeric(soil_points$Upper_Dept))
soil_points=data.frame(soil_points)
soil_points_A=soil_points[soil_points$Master_Hor=='A',]
soil_points_solum=soil_points[soil_points$Master_Hor %in% c('A', 'B'),]

soc_a=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_A$PEDON_ID)){
  temp=soil_points_A[soil_points_A$PEDON_ID==i,]
  temp=temp[temp$Organic_Ca!="NA",]
  soc=round(sum(as.numeric(temp$Organic_Ca)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, soc))
  temp=t(data.frame(temp))
  soc_a=rbind(soc_a, temp)
  }
colnames(soc_a)=c("PEDON_ID", 'soc')


tn_a=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_A$PEDON_ID)){
  temp=soil_points_A[soil_points_A$PEDON_ID==i,]
  temp=temp[temp$Total_Nitr!="NA",]
  tn=round(sum(as.numeric(temp$Total_Nitr)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, tn))
  temp=t(data.frame(temp))
  tn_a=rbind(tn_a, temp)
}
colnames(tn_a)=c("PEDON_ID", 'tn')


soc_stock=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points$PEDON_ID)){
  temp=soil_points[soil_points$PEDON_ID==i,]
  temp=temp[temp$Organic_Ca!="NA",]
  temp=temp[temp$bd!="NA",]
  soc=round(sum((as.numeric(temp$Organic_Ca)/100)*(temp$thickness)*(temp$bd)),2)*10
  temp=as.numeric(c(i, soc))
  temp=t(data.frame(temp))
  soc_stock=rbind(soc_stock, temp)
}
colnames(soc_stock)=c("PEDON_ID", 'soc_stock')

cec_a=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_A$PEDON_ID)){
  temp=soil_points_A[soil_points_A$PEDON_ID==i,]
  temp=temp[temp$CEC!="NA",]
  cec=round(sum(as.numeric(temp$CEC)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, cec))
  temp=t(data.frame(temp))
  cec_a=rbind(cec_a, temp)
}
colnames(cec_a)=c("PEDON_ID", 'cec')

ec_a=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_A$PEDON_ID)){
  temp=soil_points_A[soil_points_A$PEDON_ID==i,]
  temp=temp[temp$EC!="NA",]
  ec=round(sum(as.numeric(temp$EC)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, ec))
  temp=t(data.frame(temp))
  ec_a=rbind(ec_a, temp)
}
colnames(ec_a)=c("PEDON_ID", 'ec')

ec_profile=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points$PEDON_ID)){
  temp=soil_points[soil_points$PEDON_ID==i,]
  temp=temp[temp$EC!="NA",]
  ec=round(sum(as.numeric(temp$EC)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, ec))
  temp=t(data.frame(temp))
  ec_profile=rbind(ec_profile, temp)
}
colnames(ec_profile)=c("PEDON_ID", 'ec')



ec_solum=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_solum$PEDON_ID)){
  temp=soil_points_solum[soil_points_solum$PEDON_ID==i,]
  temp=temp[temp$EC!="NA",]
  ec=round(sum(as.numeric(temp$EC)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, ec))
  temp=t(data.frame(temp))
  ec_solum=rbind(ec_solum, temp)
}
colnames(ec_solum)=c("PEDON_ID", 'ec')


ioc_a=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_A$PEDON_ID)){
  temp=soil_points_A[soil_points_A$PEDON_ID==i,]
  temp=temp[temp$Inorganic_!="NA",]
  ioc=round(sum(as.numeric(temp$Inorganic_)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, ioc))
  temp=t(data.frame(temp))
  ioc_a=rbind(ioc_a, temp)
}
colnames(ioc_a)=c("PEDON_ID", 'ioc')

ioc_all=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points$PEDON_ID)){
  temp=soil_points[soil_points$PEDON_ID==i,]
  temp=temp[temp$Inorganic_!="NA",]
  ioc=round(sum(as.numeric(temp$Inorganic_)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, ioc))
  temp=t(data.frame(temp))
  ioc_all=rbind(ioc_all, temp)
}
colnames(ioc_all)=c("PEDON_ID", 'ioc')

clay_a=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_A$PEDON_ID)){
  temp=soil_points_A[soil_points_A$PEDON_ID==i,]
  temp=temp[temp$clay_0.2!="NA",]
  clay=round(sum(as.numeric(temp$clay_0.2)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, clay))
  temp=t(data.frame(temp))
  clay_a=rbind(clay_a, temp)
}
colnames(clay_a)=c("PEDON_ID", 'clay')

silt_a=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_A$PEDON_ID)){
  temp=soil_points_A[soil_points_A$PEDON_ID==i,]
  temp=temp[temp$silt!="NA",]
  silt=round(sum(as.numeric(temp$silt)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, silt))
  temp=t(data.frame(temp))
  silt_a=rbind(silt_a, temp)
}
colnames(silt_a)=c("PEDON_ID", 'silt')

sand_a=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_A$PEDON_ID)){
  temp=soil_points_A[soil_points_A$PEDON_ID==i,]
  temp=temp[temp$total_sand!="NA",]
  sand=round(sum(as.numeric(temp$total_sand)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, sand))
  temp=t(data.frame(temp))
  sand_a=rbind(sand_a, temp)
}
colnames(sand_a)=c("PEDON_ID", 'sand')


clay_all=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points$PEDON_ID)){
  temp=soil_points[soil_points$PEDON_ID==i,]
  temp=temp[temp$clay_0.2!="NA",]
  clay=round(sum(as.numeric(temp$clay_0.2)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, clay))
  temp=t(data.frame(temp))
  clay_all=rbind(clay_all, temp)
}
colnames(clay_all)=c("PEDON_ID", 'clay')

silt_all=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points$PEDON_ID)){
  temp=soil_points[soil_points$PEDON_ID==i,]
  temp=temp[temp$silt!="NA",]
  silt=round(sum(as.numeric(temp$silt)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, silt))
  temp=t(data.frame(temp))
  silt_all=rbind(silt_all, temp)
}
colnames(silt_all)=c("PEDON_ID", 'silt')

sand_all=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points$PEDON_ID)){
  temp=soil_points[soil_points$PEDON_ID==i,]
  temp=temp[temp$total_sand!="NA",]
  sand=round(sum(as.numeric(temp$total_sand)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, sand))
  temp=t(data.frame(temp))
  sand_all=rbind(sand_all, temp)
}
colnames(sand_all)=c("PEDON_ID", 'sand')

clay_solum=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_solum$PEDON_ID)){
  temp=soil_points_solum[soil_points_solum$PEDON_ID==i,]
  temp=temp[temp$clay_0.2!="NA",]
  clay=round(sum(as.numeric(temp$clay_0.2)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, clay))
  temp=t(data.frame(temp))
  clay_solum=rbind(clay_solum, temp)
}
colnames(clay_solum)=c("PEDON_ID", 'clay')

silt_solum=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_solum$PEDON_ID)){
  temp=soil_points_solum[soil_points_solum$PEDON_ID==i,]
  temp=temp[temp$silt!="NA",]
  silt=round(sum(as.numeric(temp$silt)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, silt))
  temp=t(data.frame(temp))
  silt_solum=rbind(silt_solum, temp)
}
colnames(silt_solum)=c("PEDON_ID", 'silt')

sand_solum=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_solum$PEDON_ID)){
  temp=soil_points_solum[soil_points_solum$PEDON_ID==i,]
  temp=temp[temp$total_sand!="NA",]
  sand=round(sum(as.numeric(temp$total_sand)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, sand))
  temp=t(data.frame(temp))
  sand_solum=rbind(sand_solum, temp)
}
colnames(sand_solum)=c("PEDON_ID", 'sand')


#create train-test split
raster_pca=raster_values[,54:84]
raster_pca[is.na(raster_pca)] <- 0


env_pca=princomp(raster_pca)
env_pca=env_pca$scores[,1:2]

nrow(env_pca)

#clhs
#test=clhs(data.frame(raster_pca),size=144)
#test_pedons=raster_values$PEDON_ID[test]

#create training files
raster_values=raster_values[,c(3,54:91)]

soc_a=merge(soc_a,raster_values, by='PEDON_ID')
tn_a=merge(tn_a,raster_values, by='PEDON_ID')
soc_stock=merge(soc_stock,raster_values, by='PEDON_ID')
cec_a=merge(cec_a,raster_values, by='PEDON_ID')
ec_a=merge(ec_a,raster_values, by='PEDON_ID')
ec_solum=merge(ec_solum,raster_values, by='PEDON_ID')
ec_profile=merge(ec_profile,raster_values, by='PEDON_ID')
ioc_a=merge(ioc_a,raster_values, by='PEDON_ID')
ioc_all=merge(ioc_all,raster_values, by='PEDON_ID')
clay_a=merge(clay_a,raster_values, by='PEDON_ID')
silt_a=merge(silt_a,raster_values, by='PEDON_ID')
sand_a=merge(sand_a,raster_values, by='PEDON_ID')
clay_solum=merge(clay_solum,raster_values, by='PEDON_ID')
silt_solum=merge(silt_solum,raster_values, by='PEDON_ID')
sand_solum=merge(sand_solum,raster_values, by='PEDON_ID')
clay_all=merge(clay_all,raster_values, by='PEDON_ID')
silt_all=merge(silt_all,raster_values, by='PEDON_ID')
sand_all=merge(sand_all,raster_values, by='PEDON_ID')

############Soil Organic Carbon Model#################
setwd('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Figures/')
'%!in%' <- function(x,y)!('%in%'(x,y))
#remove SOC values above 5, 95th percentile is 4.95, and extreme values won't be reliably predictable
#soc_a=soc_a[soc_a$soc<5,]

soc_train_a=soc_a[soc_a$PEDON_ID %!in% test_pedons,]
soc_test_a=soc_a[soc_a$PEDON_ID %in% test_pedons,]

soc_train_a=soc_train_a[,-1]

soc_train_a[is.na(soc_train_a)] <- 0
soc_test_a[is.na(soc_test_a)] <- 0

soc_train_a=do.call(data.frame, lapply(soc_train_a, function(x) replace(x, is.infinite(x), 0)))
soc_test_a=do.call(data.frame, lapply(soc_test_a, function(x) replace(x, is.infinite(x), 0)))

#feature selection
carbon_features=psm_ranger_feature_selection(x="soc", soc_train_a)
carbon_features

#final model
soc_train_a_final=soc_train_a[,carbon_features]
model_carbon=ranger(soc~., data=soc_train_a_final, importance='impurity')
write.csv(sort(importance(model_carbon), decreasing=TRUE), "soc_a_feature.csv")

pred=predict(model_carbon, data=soc_test_a)
pred=pred$predictions

#model performance
cor.test(pred, soc_test_a$soc)
sd(soc_test_a$soc)/rmse(soc_test_a$soc, pred)
summary(lm(pred~soc_test_a$soc))$r.squared
rmse(soc_test_a$soc, pred)
CCC(soc_test_a$soc, pred)$rho.c
bias(soc_test_a$soc, pred)

#create plot
plot_data_soc=data.frame(soc_test_a$soc, pred)
plot_data_soc$bare_soil=ifelse(soc_test_a$b4==0, 'N', "Y")
colnames(plot_data_soc)=c("actual", 'predicted', 'Bare_Soil_Data')
carbon_plot=ggplot(plot_data_soc, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + xlim(0, 6.5) + ylim(0, 6.5) + geom_abline() + xlab("Observed Soil Organic Carbon (%)") + ylab("Predicted Soil Organic Carbon (%)")
carbon_plot = carbon_plot + annotate("text", x=0.5, y = 6.5, label=expression(paste("R"^"2 ", "= 0.50"))) + annotate("text", x=0.5, y=6.3, label="RMSE = 0.87%") + annotate("text", x=0.5, y=6.1, label=expression(paste(rho['c'], "= 0.65"))) + annotate("text", x=0.5, y=5.9, label='Bias = -0.18')
carbon_plot

ggsave('soc_a.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(soc_train_a$soc), length(soc_test_a$soc))
rmse(soc_test_a$soc, pred_null)


############Total Nitrogen Model#################
setwd('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Figures/')
'%!in%' <- function(x,y)!('%in%'(x,y))
#remove tn values above 0.5, 95th percentile is 0.43, and extreme values won't be reliably predictable
tn_a=tn_a[tn_a$tn<0.5,]

tn_train_a=tn_a[tn_a$PEDON_ID %!in% test_pedons,]
tn_test_a=tn_a[tn_a$PEDON_ID %in% test_pedons,]

tn_train_a=tn_train_a[,-1]

tn_train_a[is.na(tn_train_a)] <- 0
tn_test_a[is.na(tn_test_a)] <- 0

tn_train_a=do.call(data.frame, lapply(tn_train_a, function(x) replace(x, is.infinite(x), 0)))
tn_test_a=do.call(data.frame, lapply(tn_test_a, function(x) replace(x, is.infinite(x), 0)))

#feature selection
tn_features=psm_ranger_feature_selection(x="tn", tn_train_a)
tn_features

#final model
tn_train_a_final=tn_train_a[,tn_features]
model_tn=ranger(tn~., data=tn_train_a_final, importance='impurity')
write.csv(sort(importance(model_tn), decreasing=TRUE), "tn_a_feature.csv")

pred=predict(model_tn, data=tn_test_a)
pred=pred$predictions

#model performance
cor.test(pred, tn_test_a$tn)
sd(tn_test_a$tn)/rmse(tn_test_a$tn, pred)
summary(lm(pred~tn_test_a$tn))$r.squared
rmse(tn_test_a$tn, pred)
CCC(tn_test_a$tn, pred)$rho.c
bias(tn_test_a$tn, pred)

#create plot
plot_data_tn=data.frame(tn_test_a$tn, pred)
plot_data_tn$bare_soil=ifelse(tn_test_a$b4==0, 'N', "Y")
colnames(plot_data_tn)=c("actual", 'predicted', 'Bare_Soil_Data')
tn_plot=ggplot(plot_data_tn, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + xlim(0, 0.5) + ylim(0, 0.5) + geom_abline() + xlab("Observed Total Nitrogen (%)") + ylab("Predicted Total Nitrogen (%)")
tn_plot = tn_plot + annotate("text", x=0.04, y = 0.5, label=expression(paste("R"^"2 ", "= 0.44"))) + annotate("text", x=0.04, y=0.48, label="RMSE = 0.08%") + annotate("text", x=0.04, y=0.46, label=expression(paste(rho['c'], "= 0.58"))) + annotate("text", x=0.04, y=0.44, label='Bias = 0.01')
tn_plot

ggsave('tn_a.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(tn_train_a$tn), length(tn_test_a$tn))
rmse(tn_test_a$tn, pred_null)

############Soil Organic Carbon Stock Model#################
soc_stock=soc_stock[!is.na(soc_stock$soc_stock),]

soc_train_stock=soc_stock[soc_stock$PEDON_ID %!in% test_pedons,]
soc_test_stock=soc_stock[soc_stock$PEDON_ID %in% test_pedons,]

soc_train_stock=soc_train_stock[,-1]

soc_train_stock[is.na(soc_train_stock)] <- 0
soc_test_stock[is.na(soc_test_stock)] <- 0

soc_train_stock=do.call(data.frame, lapply(soc_train_stock, function(x) replace(x, is.infinite(x), 0)))
soc_test_stock=do.call(data.frame, lapply(soc_test_stock, function(x) replace(x, is.infinite(x), 0)))

#feature selection
stock_carbon_features=psm_ranger_feature_selection(x="soc_stock", soc_train_stock)
stock_carbon_features

#final model
soc_train_stock_final=soc_train_stock[,stock_carbon_features]
model_carbon_stock=ranger(soc_stock~., data=soc_train_stock_final, importance='impurity')
write.csv(sort(importance(model_carbon_stock), decreasing=TRUE), "soc_stock_feature.csv")

pred=predict(model_carbon_stock, data=soc_test_stock)
pred=pred$predictions

#model performance
cor.test(pred, soc_test_stock$soc_stock)
sd(soc_test_stock$soc_stock)/rmse(soc_test_stock$soc_stock, pred)
summary(lm(pred~soc_test_stock$soc_stock))$r.squared
rmse(soc_test_stock$soc_stock, pred)
CCC(soc_test_stock$soc_stock, pred)$rho.c
bias(soc_test_stock$soc_stock, pred)

#create plot
plot_data_soc_stock=data.frame(soc_test_stock$soc_stock, pred)
plot_data_soc_stock$bare_soil=ifelse(soc_test_stock$b4==0, 'N', "Y")
colnames(plot_data_soc_stock)=c("actual", 'predicted', 'Bare_Soil_Data')
carbon_stock_plot=ggplot(plot_data_soc_stock, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + xlim(0,30) + ylim(0, 30) + geom_abline() +  xlab(expression(paste('Observed Soil Organic Carbon Stock (kg m'^"-2",")"))) + ylab(expression(paste('Predicted Soil Organic Carbon Stock (kg m'^"-2",")")))
carbon_stock_plot = carbon_stock_plot + annotate("text", x=2, y = 30, label=expression(paste("R"^"2 ", "= 0.26"))) + annotate("text", x=2, y=29, label=expression(paste("RMSE = 8.0 kg m"^"-2"))) + annotate("text", x=2, y=28, label=expression(paste(rho['c'], "= 0.41"))) + annotate("text", x=2, y=27, label='Bias = 1.12')
carbon_stock_plot

ggsave('soc_stock.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(soc_train_stock$soc_stock), length(soc_test_stock$soc_stock))
rmse(soc_test_stock$soc_stock, pred_null)


############IOC a#################
ioc_train=ioc_a[ioc_a$PEDON_ID %!in% test_pedons,]
ioc_test=ioc_a[ioc_a$PEDON_ID %in% test_pedons,]

ioc_train=ioc_train[,-1]

ioc_train[is.na(ioc_train)] <- 0
ioc_test[is.na(ioc_test)] <- 0

ioc_train=do.call(data.frame, lapply(ioc_train, function(x) replace(x, is.infinite(x), 0)))
ioc_test=do.call(data.frame, lapply(ioc_test, function(x) replace(x, is.infinite(x), 0)))

#feature selection
ioc_a_features=psm_ranger_feature_selection(x="ioc", ioc_train)
ioc_a_features

#final model
ioc_train_final=ioc_train[,ioc_a_features]
model_ioc_a=ranger(ioc~., data=ioc_train_final, importance='impurity')
write.csv(sort(importance(model_ioc_a), decreasing=TRUE), "ioc_a_feature.csv")

pred=predict(model_ioc_a, data=ioc_test)
pred=pred$predictions

#model performance
cor.test(pred, ioc_test$ioc)
sd(ioc_test$ioc)/rmse(ioc_test$ioc, pred)
summary(lm(pred~ioc_test$ioc))$r.squared
rmse(ioc_test$ioc, pred)
CCC(ioc_test$ioc, pred)$rho.c
bias(ioc_test$ioc, pred)

#create plot
plot_data_ioc_a=data.frame(ioc_test$ioc, pred)
plot_data_ioc_a$bare_soil=ifelse(ioc_test$b4==0, 'N', "Y")
colnames(plot_data_ioc_a)=c("actual", 'predicted', 'Bare_Soil_Data')
ioc_a_plot=ggplot(plot_data_ioc_a, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + xlim(0,30) + ylim(0, 30) + geom_abline() +  xlab('Observed A Horizon Inorganic Carbon (%)') + ylab('Predicted A Horizon Inorganic Carbon (%)')
ioc_a_plot = ioc_a_plot + annotate("text", x=2, y = 30, label=expression(paste("R"^"2 ", "= 0.11"))) + annotate("text", x=2, y=29, label=expression(paste("RMSE = 2.9%"))) + annotate("text", x=2, y=28, label=expression(paste(rho['c'], "= 0.27"))) + annotate("text", x=2, y=27, label='Bias = -0.43')
ioc_a_plot

ggsave('ioc_a.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(ioc_train$ioc), length(ioc_test$ioc))
rmse(ioc_test$ioc, pred_null)

############IOC All#################
ioc_train_all=ioc_all[ioc_all$PEDON_ID %!in% test_pedons,]
ioc_test_all=ioc_all[ioc_all$PEDON_ID %in% test_pedons,]

ioc_train_all=ioc_train_all[,-1]

ioc_train_all[is.na(ioc_train_all)] <- 0
ioc_test_all[is.na(ioc_test_all)] <- 0

ioc_train_all=do.call(data.frame, lapply(ioc_train_all, function(x) replace(x, is.infinite(x), 0)))
ioc_test_all=do.call(data.frame, lapply(ioc_test_all, function(x) replace(x, is.infinite(x), 0)))

#feature selection
ioc_all_features=psm_ranger_feature_selection(x="ioc", ioc_train_all)
ioc_all_features

#final model
ioc_train_all_final=ioc_train_all[,ioc_all_features]
model_ioc_all=ranger(ioc~., data=ioc_train_all_final, importance='impurity')
write.csv(sort(importance(model_ioc_all), decreasing=TRUE), "ioc_all_feature.csv")

pred=predict(model_ioc_all, data=ioc_test_all)
pred=pred$predictions

#model performance
cor.test(pred, ioc_test_all$ioc)
sd(ioc_test_all$ioc)/rmse(ioc_test_all$ioc, pred)
summary(lm(pred~ioc_test_all$ioc))$r.squared
rmse(ioc_test_all$ioc, pred)
CCC(ioc_test_all$ioc, pred)$rho.c
bias(ioc_test_all$ioc, pred)

#create plot
plot_data_ioc_all=data.frame(ioc_test_all$ioc, pred)
plot_data_ioc_all$bare_soil=ifelse(ioc_test_all$b4==0, 'N', "Y")
colnames(plot_data_ioc_all)=c("actual", 'predicted', 'Bare_Soil_Data')
ioc_all_plot=ggplot(plot_data_ioc_all, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + xlim(0,60) + ylim(0, 60) + geom_abline() + xlab('Observed Profile Inorganic Carbon (%)') + ylab('Predicted Profile Inorganic Carbon (%)')
ioc_all_plot = ioc_all_plot + annotate("text", x=4, y = 60, label=expression(paste("R"^"2 ", "= 0.59"))) + annotate("text", x=4, y=58, label=expression(paste("RMSE = 4.9%"))) + annotate("text", x=4, y=56, label=expression(paste(rho['c'], "= 0.73"))) + annotate("text", x=4, y=54, label='Bias = 0.27')
ioc_all_plot

ggsave('ioc_all.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(ioc_train_all$ioc), length(ioc_test_all$ioc))
rmse(ioc_test_all$ioc, pred_null)

##############CEC#################
cec_train_a=cec_a[cec_a$PEDON_ID %!in% test_pedons,]
cec_test_a=cec_a[cec_a$PEDON_ID %in% test_pedons,]

cec_train_a=cec_train_a[,-1]

cec_train_a[is.na(cec_train_a)] <- 0
cec_test_a[is.na(cec_test_a)] <- 0

#feature selection with ranger
cec_features=psm_ranger_feature_selection(x="cec", cec_train_a)
cec_features

#final model
cec_train_a_final=cec_train_a[,cec_features]
model_cec=ranger(cec~., data=cec_train_a_final, importance='impurity')
sort(importance(model_cec), decreasing=TRUE)
write.csv(sort(importance(model_cec), decreasing=TRUE), "cec_a_feature.csv")

pred=predict(model_cec, data=cec_test_a)
pred=pred$predictions

#model performance
cor.test(pred, cec_test_a$cec)
sd(cec_test_a$cec)/rmse(cec_test_a$cec, pred)
summary(lm(pred~cec_test_a$cec))$r.squared
rmse(cec_test_a$cec, pred)
CCC(cec_test_a$cec, pred)$rho.c
bias(cec_test_a$cec, pred)

#create plot
plot_data_cec=data.frame(cec_test_a$cec, pred)
plot_data_cec$bare_soil=ifelse(cec_test_a$b4==0, 'N', "Y")
colnames(plot_data_cec)=c("actual", 'predicted', 'Bare_Soil_Data')
cec_plot=ggplot(plot_data_cec, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + xlim(0, 60) + ylim(0, 60) + geom_abline() + xlab(expression(paste('Observed Cation Exchange Capacity (meq 100g'^"-1",")"))) + ylab(expression(paste('Predicted Cation Exchange Capacity (meq 100g'^"-1",")")))
cec_plot = cec_plot + annotate("text", x=5, y = 60, label=expression(paste("R"^"2 ", "= 0.43"))) + annotate("text", x=5, y=57, label=expression(paste("RMSE = 9.2 meq 100g"^"-1"))) + annotate("text", x=5, y=54, label=expression(paste(rho['c'], "= 0.60"))) + annotate("text", x=5, y=51, label='Bias = -0.90')
cec_plot

ggsave('cec_a.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(cec_train_a$cec), length(cec_test_a$cec))
rmse(cec_test_a$cec, pred_null)

#############EC a###################
ec_train_a=ec_a[ec_a$PEDON_ID %!in% test_pedons,]
ec_test_a=ec_a[ec_a$PEDON_ID %in% test_pedons,]

ec_train_a=ec_train_a[,-1]

ec_train_a[is.na(ec_train_a)] <- 0
ec_test_a[is.na(ec_test_a)] <- 0

#feature selection with ranger
ec_a_features=psm_ranger_feature_selection(x="ec", ec_train_a)
ec_a_features

#final model
ec_train_a_final=ec_train_a[,ec_a_features]
model_ec_a=ranger(ec~., data=ec_train_a_final, importance='impurity')
sort(importance(model_ec_a), decreasing=TRUE)
write.csv(sort(importance(model_ec_a), decreasing=TRUE), "ec_a_feature.csv")

pred=predict(model_ec_a, data=ec_test_a)
pred=pred$predictions

#model performance
cor.test(pred, ec_test_a$ec)
sd(ec_test_a$ec)/rmse(ec_test_a$ec, pred)
summary(lm(pred~ec_test_a$ec))$r.squared
rmse(ec_test_a$ec, pred)
CCC(ec_test_a$ec, pred)$rho.c
bias(ec_test_a$ec, pred)

#create plot
plot_data_ec_a=data.frame(ec_test_a$ec, pred)
plot_data_ec_a$bare_soil=ifelse(ec_test_a$b4==0, 'N', "Y")
colnames(plot_data_ec_a)=c("actual", 'predicted', 'Bare_Soil_Data')
ec_plot_a=ggplot(plot_data_ec_a, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + xlim(0, 18) + ylim(0, 18) + geom_abline() + xlab(expression(paste('Observed A Horizon Electrical Conductivity (dS m'^"-1",")"))) + ylab(expression(paste('Predicted A Horizon Electrical Conductivity (dS m'^"-1",")")))
ec_plot_a = ec_plot_a + annotate("text", x=2, y = 18, label=expression(paste("R"^"2 ", "= 0.09"))) + annotate("text", x=2, y=17, label=expression(paste("RMSE = 1.6 dS m"^"-1"))) + annotate("text", x=2, y=16, label=expression(paste(rho['c'], "= 0.13"))) + annotate("text", x=2, y=15, label='Bias = 0.01')
ec_plot_a

ggsave('ec_a.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(ec_train_a$ec), length(ec_test_a$ec))
rmse(ec_test_a$ec, pred_null)


###########EC solum#############
ec_train_solum=ec_solum[ec_solum$PEDON_ID %!in% test_pedons,]
ec_test_solum=ec_solum[ec_solum$PEDON_ID %in% test_pedons,]

ec_train_solum=ec_train_solum[,-1]

ec_train_solum[is.na(ec_train_solum)] <- 0
ec_test_solum[is.na(ec_test_solum)] <- 0

#feature selection with ranger
ec_solum_features=psm_ranger_feature_selection(x="ec", ec_train_solum)
ec_solum_features

#final model
ec_train_solum_final=ec_train_solum[,ec_solum_features]
model_ec_solum=ranger(ec~., data=ec_train_solum_final, importance='impurity')
sort(importance(model_ec_solum), decreasing=TRUE)
write.csv(sort(importance(model_ec_solum), decreasing=TRUE), "ec_solum_feature.csv")

pred=predict(model_ec_solum, data=ec_test_solum)
pred=pred$predictions

#model performance
cor.test(pred, ec_test_solum$ec)
sd(ec_test_solum$ec)/rmse(ec_test_solum$ec, pred)
summary(lm(pred~ec_test_solum$ec))$r.squared
rmse(ec_test_solum$ec, pred)
CCC(ec_test_solum$ec, pred)$rho.c
bias(ec_test_solum$ec, pred)

#create plot
plot_data_ec_solum=data.frame(ec_test_solum$ec, pred)
plot_data_ec_solum$bare_soil=ifelse(ec_test_solum$b4==0, 'N', "Y")
colnames(plot_data_ec_solum)=c("actual", 'predicted', 'Bare_Soil_Data')
ec_plot_solum= ggplot(plot_data_ec_solum, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + xlim(0, 18) + ylim(0, 18) + geom_abline() + xlab(expression(paste('Observed Solum Electrical Conductivity (dS m'^"-1",")"))) + ylab(expression(paste('Predicted Solum Electrical Conductivity (dS m'^"-1",")")))
ec_plot_solum = ec_plot_solum + annotate("text", x=2, y = 18, label=expression(paste("R"^"2 ", "= 0.14"))) + annotate("text", x=2, y=17, label=expression(paste("RMSE = 1.8 dS m"^"-1"))) + annotate("text", x=2, y=16, label=expression(paste(rho['c'], "= 0.28"))) + annotate("text", x=2, y=15, label='Bias = -0.20')
ec_plot_solum

ggsave('ec_solum.png', plot=last_plot(), width=11, height=8.5)


#null model
pred_null=rep(mean(ec_train_solum$ec), length(ec_test_solum$ec))
rmse(ec_test_solum$ec, pred_null)

###########EC Profile#############
ec_train_profile=ec_profile[ec_profile$PEDON_ID %!in% test_pedons,]
ec_test_profile=ec_profile[ec_profile$PEDON_ID %in% test_pedons,]

ec_train_profile=ec_train_profile[,-1]

ec_train_profile[is.na(ec_train_profile)] <- 0
ec_test_profile[is.na(ec_test_profile)] <- 0

#feature selection with ranger
ec_profile_features=psm_ranger_feature_selection(x="ec", ec_train_profile)
ec_profile_features

#final model
ec_train_profile_final=ec_train_profile[,ec_profile_features]
model_ec_profile=ranger(ec~., data=ec_train_profile_final, importance='impurity')
sort(importance(model_ec_profile), decreasing=TRUE)
write.csv(sort(importance(model_ec_profile), decreasing=TRUE), "ec_profile_feature.csv")

pred=predict(model_ec_profile, data=ec_test_profile)
pred=pred$predictions

#model performance
cor.test(pred, ec_test_profile$ec)
sd(ec_test_profile$ec)/rmse(ec_test_profile$ec, pred)
summary(lm(pred~ec_test_profile$ec))$r.squared
rmse(ec_test_profile$ec, pred)
CCC(ec_test_profile$ec, pred)$rho.c
bias(ec_test_profile$ec, pred)

#create plot
plot_data_ec_profile=data.frame(ec_test_profile$ec, pred)
plot_data_ec_profile$bare_soil=ifelse(ec_test_profile$b4==0, 'N', "Y")
colnames(plot_data_ec_profile)=c("actual", 'predicted', 'Bare_Soil_Data')
ec_plot_profile= ggplot(plot_data_ec_profile, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + xlim(0, 18) + ylim(0, 18) + geom_abline() + xlab(expression(paste('Observed Profile Electrical Conductivity (dS m'^"-1",")"))) + ylab(expression(paste('Predicted Profile Electrical Conductivity (dS m'^"-1",")")))
ec_plot_profile = ec_plot_profile + annotate("text", x=2, y = 18, label=expression(paste("R"^"2 ", "= 0.37"))) + annotate("text", x=2, y=17, label=expression(paste("RMSE = 2.3 dS m"^"-1"))) + annotate("text", x=2, y=16, label=expression(paste(rho['c'], "= 0.57"))) + annotate("text", x=2, y=15, label='Bias = -0.78')
ec_plot_profile

ggsave('ec_profile.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(ec_train_profile$ec), length(ec_test_profile$ec))
rmse(ec_test_profile$ec, pred_null)

##############Clay A#################
clay_train_a=clay_a[clay_a$PEDON_ID %!in% test_pedons,]
clay_test_a=clay_a[clay_a$PEDON_ID %in% test_pedons,]

clay_train_a=clay_train_a[,-1]

clay_train_a[is.na(clay_train_a)] <- 0
clay_test_a[is.na(clay_test_a)] <- 0

#feature selection with ranger
clay_a_features=psm_ranger_feature_selection(x="clay", clay_train_a)
clay_a_features

#final model
clay_train_a_final=clay_train_a[,clay_a_features]

#final model
model_clay=ranger(clay~., data=clay_train_a_final, importance='impurity')
sort(importance(model_clay), decreasing=TRUE)
write.csv(sort(importance(model_clay), decreasing=TRUE), "clay_a_feature.csv")

pred=predict(model_clay, data=clay_test_a)
pred=pred$predictions

#model performance
cor.test(pred, clay_test_a$clay)
sd(clay_test_a$clay)/rmse(clay_test_a$clay, pred)
summary(lm(pred~clay_test_a$clay))$r.squared
rmse(clay_test_a$clay, pred)
CCC(clay_test_a$clay, pred)$rho.c
bias(clay_test_a$clay, pred)

#create plot
plot_data_clay_a=data.frame(clay_test_a$clay, pred)
plot_data_clay_a$bare_soil=ifelse(clay_test_a$b4==0, 'N', "Y")
colnames(plot_data_clay_a)=c("actual", 'predicted', 'Bare_Soil_Data')
clay_plot_a= ggplot(plot_data_clay_a, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Clay Content (%)'))) + ylab(expression(paste('Predicted Clay Content (%)')))
clay_plot_a = clay_plot_a + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.66"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 8.2%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.75"))) + annotate("text", x=5, y=70, label='Bias = 0.33')
clay_plot_a

ggsave('clay_a.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(clay_train_a$clay), length(clay_test_a$clay))
rmse(clay_test_a$clay, pred_null)

##############Silt A#################
silt_train_a=silt_a[silt_a$PEDON_ID %!in% test_pedons,]
silt_test_a=silt_a[silt_a$PEDON_ID %in% test_pedons,]

silt_train_a=silt_train_a[,-1]

silt_train_a[is.na(silt_train_a)] <- 0
silt_test_a[is.na(silt_test_a)] <- 0

#feature selection with ranger
silt_a_features=psm_ranger_feature_selection(x="silt", silt_train_a)
silt_a_features

#final model
silt_train_a_final=silt_train_a[,silt_a_features]

#final model
model_silt=ranger(silt~., data=silt_train_a_final, importance='impurity')
sort(importance(model_silt), decreasing=TRUE)
write.csv(sort(importance(model_silt), decreasing=TRUE), "silt_a_feature.csv")

pred=predict(model_silt, data=silt_test_a)
pred=pred$predictions

#model performance
cor.test(pred, silt_test_a$silt)
sd(silt_test_a$silt)/rmse(silt_test_a$silt, pred)
summary(lm(pred~silt_test_a$silt))$r.squared
rmse(silt_test_a$silt, pred)
CCC(silt_test_a$silt, pred)$rho.c
bias(silt_test_a$silt, pred)

#create plot
plot_data_silt_a=data.frame(silt_test_a$silt, pred)
plot_data_silt_a$bare_soil=ifelse(silt_test_a$b4==0, 'N', "Y")
colnames(plot_data_silt_a)=c("actual", 'predicted', 'Bare_Soil_Data')
silt_plot_a= ggplot(plot_data_silt_a, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Silt Content (%)'))) + ylab(expression(paste('Predicted Silt Content (%)')))
silt_plot_a = silt_plot_a + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.31"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 11.9"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.50"))) + annotate("text", x=5, y=70, label='Bias = 1.19')
silt_plot_a

ggsave('silt_a.png', plot=last_plot(), width=11, height=8.5)


#null model
pred_null=rep(mean(silt_train_a$silt), length(silt_test_a$silt))
rmse(silt_test_a$silt, pred_null)

##############Sand A#################
sand_train_a=sand_a[sand_a$PEDON_ID %!in% test_pedons,]
sand_test_a=sand_a[sand_a$PEDON_ID %in% test_pedons,]

sand_train_a=sand_train_a[,-1]

sand_train_a[is.na(sand_train_a)] <- 0
sand_test_a[is.na(sand_test_a)] <- 0

#feature selection with ranger
sand_a_features=psm_ranger_feature_selection(x="sand", sand_train_a)
sand_a_features

#final model
sand_train_a_final=sand_train_a[,sand_a_features]

#final model
model_sand=ranger(sand~., data=sand_train_a_final, importance='impurity')
sort(importance(model_sand), decreasing=TRUE)
write.csv(sort(importance(model_sand), decreasing=TRUE), "sand_a_feature.csv")

pred=predict(model_sand, data=sand_test_a)
pred=pred$predictions

#model performance
cor.test(pred, sand_test_a$sand)
sd(sand_test_a$sand)/rmse(sand_test_a$sand, pred)
summary(lm(pred~sand_test_a$sand))$r.squared
rmse(sand_test_a$sand, pred)
CCC(sand_test_a$sand, pred)$rho.c
bias(sand_test_a$sand, pred)

#create plot
plot_data_sand_a=data.frame(sand_test_a$sand, pred)
plot_data_sand_a$bare_soil=ifelse(sand_test_a$b4==0, 'N', "Y")
colnames(plot_data_sand_a)=c("actual", 'predicted', 'Bare_Soil_Data')
sand_plot_a= ggplot(plot_data_sand_a, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Sand Content (%)'))) + ylab(expression(paste('Predicted Sand Content (%)')))
sand_plot_a = sand_plot_a + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.50"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 15.3%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.67"))) + annotate("text", x=5, y=70, label='Bias = -1.63')
sand_plot_a

ggsave('sand_a.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(sand_train_a$sand), length(sand_test_a$sand))
rmse(sand_test_a$sand, pred_null)

##############Clay Solum#################
clay_train_solum=clay_solum[clay_solum$PEDON_ID %!in% test_pedons,]
clay_test_solum=clay_solum[clay_solum$PEDON_ID %in% test_pedons,]

clay_train_solum=clay_train_solum[,-1]

clay_train_solum[is.na(clay_train_solum)] <- 0
clay_test_solum[is.na(clay_test_solum)] <- 0

#feature selection with ranger
clay_solum_features=psm_ranger_feature_selection(x="clay", clay_train_solum)
clay_solum_features

#final model
clay_train_solum_final=clay_train_solum[,clay_solum_features]

#final model
model_clay_solum=ranger(clay~., data=clay_train_solum_final, importance='impurity')
sort(importance(model_clay_solum), decreasing=TRUE)
write.csv(sort(importance(model_clay_solum), decreasing=TRUE), "clay_solum_feature.csv")

pred=predict(model_clay_solum, data=clay_test_solum)
pred=pred$predictions

#model performance
cor.test(pred, clay_test_solum$clay)
sd(clay_test_solum$clay)/rmse(clay_test_solum$clay, pred)
summary(lm(pred~clay_test_solum$clay))$r.squared
rmse(clay_test_solum$clay, pred)
CCC(clay_test_solum$clay, pred)$rho.c
bias(clay_test_solum$clay, pred)

#create plot
plot_data_clay_solum=data.frame(clay_test_solum$clay, pred)
plot_data_clay_solum$bare_soil=ifelse(clay_test_solum$b4==0, 'N', "Y")
colnames(plot_data_clay_solum)=c("actual", 'predicted', 'Bare_Soil_Data')
clay_plot_solum= ggplot(plot_data_clay_solum, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Clay Content (%)'))) + ylab(expression(paste('Predicted Clay Content (%)')))
clay_plot_solum = clay_plot_solum + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.58"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 9.6%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.71"))) + annotate("text", x=5, y=70, label='Bias = 0.75')
clay_plot_solum

ggsave('clay_solum.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(clay_train_solum$clay), length(clay_test_solum$clay))
rmse(clay_test_solum$clay, pred_null)

##############Silt Solum#################
silt_train_solum=silt_solum[silt_solum$PEDON_ID %!in% test_pedons,]
silt_test_solum=silt_solum[silt_solum$PEDON_ID %in% test_pedons,]

silt_train_solum=silt_train_solum[,-1]

silt_train_solum[is.na(silt_train_solum)] <- 0
silt_test_solum[is.na(silt_test_solum)] <- 0

#feature selection with ranger
silt_solum_features=psm_ranger_feature_selection(x="silt", silt_train_solum)
silt_solum_features

#final model
silt_train_solum_final=silt_train_solum[,silt_solum_features]

#final model
model_silt_solum=ranger(silt~., data=silt_train_solum_final, importance='impurity')
sort(importance(model_silt_solum), decreasing=TRUE)
write.csv(sort(importance(model_silt_solum), decreasing=TRUE), "silt_solum_feature.csv")

pred=predict(model_silt_solum, data=silt_test_solum)
pred=pred$predictions

#model performance
cor.test(pred, silt_test_solum$silt)
sd(silt_test_solum$silt)/rmse(silt_test_solum$silt, pred)
summary(lm(pred~silt_test_solum$silt))$r.squared
rmse(silt_test_solum$silt, pred)
CCC(silt_test_solum$silt, pred)$rho.c
bias(silt_test_solum$silt, pred)

#create plot
plot_data_silt_solum=data.frame(silt_test_solum$silt, pred)
plot_data_silt_solum$bare_soil=ifelse(silt_test_solum$b4==0, 'N', "Y")
colnames(plot_data_silt_solum)=c("actual", 'predicted', 'Bare_Soil_Data')
silt_plot_solum= ggplot(plot_data_silt_solum, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Silt Content (%)'))) + ylab(expression(paste('Predicted Silt Content (%)')))
silt_plot_solum = silt_plot_solum + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.34"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 11.5%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.54"))) + annotate("text", x=5, y=70, label='Bias = 1.46')
silt_plot_solum

ggsave('silt_solum.png', plot=last_plot(), width=11, height=8.5)


#null model
pred_null=rep(mean(silt_train_solum$silt), length(silt_test_solum$silt))
rmse(silt_test_solum$silt, pred_null)


##############Sand Solum#################
sand_train_solum=sand_solum[sand_solum$PEDON_ID %!in% test_pedons,]
sand_test_solum=sand_solum[sand_solum$PEDON_ID %in% test_pedons,]

sand_train_solum=sand_train_solum[,-1]

sand_train_solum[is.na(sand_train_solum)] <- 0
sand_test_solum[is.na(sand_test_solum)] <- 0

#feature selection with ranger
sand_solum_features=psm_ranger_feature_selection(x="sand", sand_train_solum)
sand_solum_features

#final model
sand_train_solum_final=sand_train_solum[,sand_solum_features]

#final model
model_sand_solum=ranger(sand~., data=sand_train_solum_final, importance='impurity')
sort(importance(model_sand_solum), decreasing=TRUE)
write.csv(sort(importance(model_sand_solum), decreasing=TRUE), "sand_solum_feature.csv")

pred=predict(model_sand_solum, data=sand_test_solum)
pred=pred$predictions

#model performance
cor.test(pred, sand_test_solum$sand)
sd(sand_test_solum$sand)/rmse(sand_test_solum$sand, pred)
summary(lm(pred~sand_test_solum$sand))$r.squared
rmse(sand_test_solum$sand, pred)
CCC(sand_test_solum$sand, pred)$rho.c
bias(sand_test_solum$sand, pred)

#create plot
plot_data_sand_solum=data.frame(sand_test_solum$sand, pred)
plot_data_sand_solum$bare_soil=ifelse(sand_test_solum$b4==0, 'N', "Y")
colnames(plot_data_sand_solum)=c("actual", 'predicted', 'Bare_Soil_Data')
sand_plot_solum= ggplot(plot_data_sand_solum, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Sand Content (%)'))) + ylab(expression(paste('Predicted Sand Content (%)')))
sand_plot_solum = sand_plot_solum + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.48"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 16.0%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.66"))) + annotate("text", x=5, y=70, label='Bias = -2.46')
sand_plot_solum

ggsave('sand_solum.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(sand_train_solum$sand), length(sand_test_solum$sand))
rmse(sand_test_solum$sand, pred_null)


##############Clay all#################
clay_train_all=clay_all[clay_all$PEDON_ID %!in% test_pedons,]
clay_test_all=clay_all[clay_all$PEDON_ID %in% test_pedons,]

clay_train_all=clay_train_all[,-1]

clay_train_all[is.na(clay_train_all)] <- 0
clay_test_all[is.na(clay_test_all)] <- 0

#feature selection with ranger
clay_all_features=psm_ranger_feature_selection(x="clay", clay_train_all)
clay_all_features

#final model
clay_train_all_final=clay_train_all[,clay_all_features]

#final model
model_clay_all=ranger(clay~., data=clay_train_all_final, importance='impurity')
sort(importance(model_clay_all), decreasing=TRUE)
write.csv(sort(importance(model_clay_all), decreasing=TRUE), "clay_all_feature.csv")

pred=predict(model_clay_all, data=clay_test_all)
pred=pred$predictions

#model performance
cor.test(pred, clay_test_all$clay)
sd(clay_test_all$clay)/rmse(clay_test_all$clay, pred)
summary(lm(pred~clay_test_all$clay))$r.squared
rmse(clay_test_all$clay, pred)
CCC(clay_test_all$clay, pred)$rho.c
bias(clay_test_all$clay, pred)

#create plot
plot_data_clay_all=data.frame(clay_test_all$clay, pred)
plot_data_clay_all$bare_soil=ifelse(clay_test_all$b4==0, 'N', "Y")
colnames(plot_data_clay_all)=c("actual", 'predicted', 'Bare_Soil_Data')
clay_plot_all= ggplot(plot_data_clay_all, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Clay Content (%)'))) + ylab(expression(paste('Predicted Clay Content (%)')))
clay_plot_all = clay_plot_all + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.34"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 14.5%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.51"))) + annotate("text", x=5, y=70, label='Bias = 1.77')
clay_plot_all

ggsave('clay_all.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(clay_train_all$clay), length(clay_test_all$clay))
rmse(clay_test_all$clay, pred_null)


##############Silt all#################
silt_train_all=silt_all[silt_all$PEDON_ID %!in% test_pedons,]
silt_test_all=silt_all[silt_all$PEDON_ID %in% test_pedons,]

silt_train_all=silt_train_all[,-1]

silt_train_all[is.na(silt_train_all)] <- 0
silt_test_all[is.na(silt_test_all)] <- 0

#feature selection with ranger
silt_all_features=psm_ranger_feature_selection(x="silt", silt_train_all)
silt_all_features

#final model
silt_train_all_final=silt_train_all[,silt_all_features]

#final model
model_silt_all=ranger(silt~., data=silt_train_all_final, importance='impurity')
sort(importance(model_silt_all), decreasing=TRUE)
write.csv(sort(importance(model_silt_all), decreasing=TRUE), "silt_all_feature.csv")

pred=predict(model_silt_all, data=silt_test_all)
pred=pred$predictions

#model performance
cor.test(pred, silt_test_all$silt)
sd(silt_test_all$silt)/rmse(silt_test_all$silt, pred)
summary(lm(pred~silt_test_all$silt))$r.squared
rmse(silt_test_all$silt, pred)
CCC(silt_test_all$silt, pred)$rho.c
bias(silt_test_all$silt, pred)

#create plot
plot_data_silt_all=data.frame(silt_test_all$silt, pred)
plot_data_silt_all$bare_soil=ifelse(silt_test_all$b4==0, 'N', "Y")
colnames(plot_data_silt_all)=c("actual", 'predicted', 'Bare_Soil_Data')
silt_plot_all= ggplot(plot_data_silt_all, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Silt Content (%)'))) + ylab(expression(paste('Predicted Silt Content (%)')))
silt_plot_all = silt_plot_all + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.34"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 14.1%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.52"))) + annotate("text", x=5, y=70, label='Bias = 2.11')
silt_plot_all

ggsave('silt_all.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(silt_train_all$silt), length(silt_test_all$silt))
rmse(silt_test_all$silt, pred_null)


##############Sand all#################
sand_train_all=sand_all[sand_all$PEDON_ID %!in% test_pedons,]
sand_test_all=sand_all[sand_all$PEDON_ID %in% test_pedons,]

sand_train_all=sand_train_all[,-1]

sand_train_all[is.na(sand_train_all)] <- 0
sand_test_all[is.na(sand_test_all)] <- 0

#feature selection with ranger
sand_all_features=psm_ranger_feature_selection(x="sand", sand_train_all)
sand_all_features

#final model
sand_train_all_final=sand_train_all[,sand_all_features]

#final model
model_sand_all=ranger(sand~., data=sand_train_all_final, importance='impurity')
sort(importance(model_sand_all), decreasing=TRUE)
write.csv(sort(importance(model_sand_all), decreasing=TRUE), "sand_all_feature.csv")

pred=predict(model_sand_all, data=sand_test_all)
pred=pred$predictions

#model performance
cor.test(pred, sand_test_all$sand)
sd(sand_test_all$sand)/rmse(sand_test_all$sand, pred)
summary(lm(pred~sand_test_all$sand))$r.squared
rmse(sand_test_all$sand, pred)
CCC(sand_test_all$sand, pred)$rho.c
bias(sand_test_all$sand, pred)

#create plot
plot_data_sand_all=data.frame(sand_test_all$sand, pred)
plot_data_sand_all$bare_soil=ifelse(sand_test_all$b4==0, 'N', "Y")
colnames(plot_data_sand_all)=c("actual", 'predicted', 'Bare_Soil_Data')
sand_plot_all= ggplot(plot_data_sand_all, (aes(x=actual, y=predicted, shape=Bare_Soil_Data))) + geom_point(size=2) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Sand Content (%)'))) + ylab(expression(paste('Predicted Sand Content (%)')))
sand_plot_all = sand_plot_all + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.42"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 18.4%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.61"))) + annotate("text", x=5, y=70, label='Bias = -1.70')
sand_plot_all

ggsave('sand_all.png', plot=last_plot(), width=11, height=8.5)

#null model
pred_null=rep(mean(sand_train_all$sand), length(sand_test_all$sand))
rmse(sand_test_all$sand, pred_null)
