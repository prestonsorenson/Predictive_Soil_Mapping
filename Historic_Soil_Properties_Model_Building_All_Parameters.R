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
library(gridExtra)

setwd('C:\\Users\\prest\\OneDrive\\Papers\\Historic_Soil_Properties\\Figures')
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
  
  if(ncol(var_train_a_init)==0){
    var_train_a_init=y[,temp_init_var]
  }
  
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

load('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Analysis/Analysis_2022_04_27.RData')

#####################Bulk Density Pedotransfer Function################################
info=read.csv('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Data/npdb_v3_csv/Info.csv')
phys=read.csv("/home/preston/OneDrive/Papers/Historic_Soil_Properties/Data/npdb_v3_csv/Physical.csv")
chem=read.csv("/home/preston/OneDrive/Papers/Historic_Soil_Properties/Data/npdb_v3_csv/Chemical.csv")


npdb=merge(info, phys, by='PEDON_ID')
npdb=merge(npdb, chem, by='PEDON_ID')
npdb=npdb[!is.na(npdb$BULK_DEN),]
npdb=npdb[!duplicated(npdb),]

npdb_sub=c('BULK_DEN', 'CARB_ORG', 'T_CLAY', 'T_SILT', "T_SAND")
npdb_bd=npdb[,npdb_sub]
npdb_bd=na.omit(npdb_bd)
npdb_bd=npdb_bd[npdb_bd$BULK_DEN<2.6,]
npdb_bd=npdb_bd[npdb_bd$CARB_ORG<10,]

npdb_train=kenStone(npdb_bd[,2:5], k=round((10257*0.75), 0))
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
bd_plot=ggplot(plot_data_bd, (aes(x=actual, y=predicted))) + geom_point(size=4) + xlim(0, 2.65) + ylim(0, 2.65) + geom_abline() + xlab(expression(paste('Observed Bulk Density (g cm'^"-3",")"))) + ylab(expression(paste('Predicted Bulk Density (g cm'^"-3",")"))) + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))+ggtitle("Bulk Density")
#bd_plot = bd_plot + annotate("text", x=0.1, y = 2.65, label=expression(paste("R"^"2 ", "= 0.52"))) + annotate("text", x=0.1, y=2.55, label=expression(paste("RMSE = 0.22 g cm"^"-3"))) + annotate("text", x=0.1, y=2.45, label=expression(paste(rho['c'], "= 0.66"))) + annotate("text", x=0.1, y=2.35, label='Bias = 0.002')
bd_plot

ggsave('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Figures/3D/bulk_denisty.png', plot=last_plot(), width=11, height=8.5)


pred_null=rep(mean(npdb_train$BULK_DEN), length(npdb_test$BULK_DEN))
rmse(npdb_test$BULK_DEN, pred_null)


###############Data Preparation################
soil_points=shapefile('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Data/Sask_NPDB_Chemical_Physical_Data.shp')
raster_values=read.csv('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Data/Historical_Soil_Properties_Training_Data_2021_10_01.csv')

#predict bulk density values
bd_input=soil_points[,c(29,41,42,47)]
bd_input[,1:4]=apply(bd_input[,1:4], 2, function(x) as.numeric(x))
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
soil_points_B=soil_points[soil_points$Master_Hor=='B',]
soil_points_C=soil_points[soil_points$Master_Hor=='C',]

soil_points_solum=soil_points[soil_points$Master_Hor %in% c('A', 'B'),]

#soc
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

soc_b=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_B$PEDON_ID)){
  temp=soil_points_B[soil_points_B$PEDON_ID==i,]
  temp=temp[temp$Organic_Ca!="NA",]
  soc=round(sum(as.numeric(temp$Organic_Ca)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, soc))
  temp=t(data.frame(temp))
  soc_b=rbind(soc_b, temp)
}
colnames(soc_b)=c("PEDON_ID", 'soc')

soc_c=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_C$PEDON_ID)){
  temp=soil_points_C[soil_points_C$PEDON_ID==i,]
  temp=temp[temp$Organic_Ca!="NA",]
  soc=round(sum(as.numeric(temp$Organic_Ca)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, soc))
  temp=t(data.frame(temp))
  soc_c=rbind(soc_c, temp)
}
colnames(soc_c)=c("PEDON_ID", 'soc')

#tn
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

tn_b=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_B$PEDON_ID)){
  temp=soil_points_B[soil_points_B$PEDON_ID==i,]
  temp=temp[temp$Total_Nitr!="NA",]
  tn=round(sum(as.numeric(temp$Total_Nitr)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, tn))
  temp=t(data.frame(temp))
  tn_b=rbind(tn_b, temp)
}
colnames(tn_b)=c("PEDON_ID", 'tn')

tn_c=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_C$PEDON_ID)){
  temp=soil_points_C[soil_points_C$PEDON_ID==i,]
  temp=temp[temp$Total_Nitr!="NA",]
  tn=round(sum(as.numeric(temp$Total_Nitr)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, tn))
  temp=t(data.frame(temp))
  tn_c=rbind(tn_c, temp)
}
colnames(tn_c)=c("PEDON_ID", 'tn')


#cec
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

cec_b=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_B$PEDON_ID)){
  temp=soil_points_B[soil_points_B$PEDON_ID==i,]
  temp=temp[temp$CEC!="NA",]
  cec=round(sum(as.numeric(temp$CEC)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, cec))
  temp=t(data.frame(temp))
  cec_b=rbind(cec_b, temp)
}
colnames(cec_b)=c("PEDON_ID", 'cec')

cec_c=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_C$PEDON_ID)){
  temp=soil_points_C[soil_points_C$PEDON_ID==i,]
  temp=temp[temp$CEC!="NA",]
  cec=round(sum(as.numeric(temp$CEC)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, cec))
  temp=t(data.frame(temp))
  cec_c=rbind(cec_c, temp)
}
colnames(cec_c)=c("PEDON_ID", 'cec')


#EC
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

ec_b=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_B$PEDON_ID)){
  temp=soil_points_B[soil_points_B$PEDON_ID==i,]
  temp=temp[temp$EC!="NA",]
  ec=round(sum(as.numeric(temp$EC)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, ec))
  temp=t(data.frame(temp))
  ec_b=rbind(ec_b, temp)
}
colnames(ec_b)=c("PEDON_ID", 'ec')

ec_c=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_C$PEDON_ID)){
  temp=soil_points_C[soil_points_C$PEDON_ID==i,]
  temp=temp[temp$EC!="NA",]
  ec=round(sum(as.numeric(temp$EC)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, ec))
  temp=t(data.frame(temp))
  ec_c=rbind(ec_c, temp)
}
colnames(ec_c)=c("PEDON_ID", 'ec')

#ioc
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

ioc_b=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_B$PEDON_ID)){
  temp=soil_points_B[soil_points_B$PEDON_ID==i,]
  temp=temp[temp$Inorganic_!="NA",]
  ioc=round(sum(as.numeric(temp$Inorganic_)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, ioc))
  temp=t(data.frame(temp))
  ioc_b=rbind(ioc_b, temp)
}
colnames(ioc_b)=c("PEDON_ID", 'ioc')

ioc_c=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_C$PEDON_ID)){
  temp=soil_points_C[soil_points_C$PEDON_ID==i,]
  temp=temp[temp$Inorganic_!="NA",]
  ioc=round(sum(as.numeric(temp$Inorganic_)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, ioc))
  temp=t(data.frame(temp))
  ioc_c=rbind(ioc_c, temp)
}
colnames(ioc_c)=c("PEDON_ID", 'ioc')

#clay
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

clay_b=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_B$PEDON_ID)){
  temp=soil_points_B[soil_points_B$PEDON_ID==i,]
  temp=temp[temp$clay_0.2!="NA",]
  clay=round(sum(as.numeric(temp$clay_0.2)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, clay))
  temp=t(data.frame(temp))
  clay_b=rbind(clay_b, temp)
}
colnames(clay_b)=c("PEDON_ID", 'clay')

clay_c=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_C$PEDON_ID)){
  temp=soil_points_C[soil_points_C$PEDON_ID==i,]
  temp=temp[temp$clay_0.2!="NA",]
  clay=round(sum(as.numeric(temp$clay_0.2)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, clay))
  temp=t(data.frame(temp))
  clay_c=rbind(clay_c, temp)
}
colnames(clay_c)=c("PEDON_ID", 'clay')

#silt
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

silt_b=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_B$PEDON_ID)){
  temp=soil_points_B[soil_points_B$PEDON_ID==i,]
  temp=temp[temp$silt!="NA",]
  silt=round(sum(as.numeric(temp$silt)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, silt))
  temp=t(data.frame(temp))
  silt_b=rbind(silt_b, temp)
}
colnames(silt_b)=c("PEDON_ID", 'silt')

silt_c=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_C$PEDON_ID)){
  temp=soil_points_C[soil_points_C$PEDON_ID==i,]
  temp=temp[temp$silt!="NA",]
  silt=round(sum(as.numeric(temp$silt)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, silt))
  temp=t(data.frame(temp))
  silt_c=rbind(silt_c, temp)
}
colnames(silt_c)=c("PEDON_ID", 'silt')

#sand
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

sand_b=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_B$PEDON_ID)){
  temp=soil_points_B[soil_points_B$PEDON_ID==i,]
  temp=temp[temp$total_sand!="NA",]
  sand=round(sum(as.numeric(temp$total_sand)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, sand))
  temp=t(data.frame(temp))
  sand_b=rbind(sand_b, temp)
}
colnames(sand_b)=c("PEDON_ID", 'sand')

sand_c=data.frame(matrix(nrow=0, ncol=2))
for (i in unique(soil_points_C$PEDON_ID)){
  temp=soil_points_C[soil_points_C$PEDON_ID==i,]
  temp=temp[temp$total_sand!="NA",]
  sand=round(sum(as.numeric(temp$total_sand)*(temp$thickness/sum(temp$thickness))),2)
  temp=as.numeric(c(i, sand))
  temp=t(data.frame(temp))
  sand_c=rbind(sand_c, temp)
}
colnames(sand_c)=c("PEDON_ID", 'sand')

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


#horizon thickness
a_thick=aggregate(soil_points_A$thickness, by=list(soil_points_A$PEDON_ID), FUN = function(x) sum(x)/100)
b_thick=aggregate(soil_points_B$thickness, by=list(soil_points_B$PEDON_ID), FUN = function(x) sum(x)/100)
c_thick=aggregate(soil_points_C$thickness, by=list(soil_points_C$PEDON_ID), FUN = function(x) sum(x)/100)

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
raster_values=apply(raster_values, 2, FUN = function(x) as.numeric(x))
raster_values=data.frame(raster_values)

#add horizon labels
soc_a$hzn='A'
tn_a$hzn='A'
cec_a$hzn='A'
ec_a$hzn='A'
ioc_a$hzn='A'
clay_a$hzn='A'
silt_a$hzn='A'
sand_a$hzn='A'
a_thick$hzn='A'
  
soc_b$hzn='B'
tn_b$hzn='B'
cec_b$hzn='B'
ec_b$hzn='B'
ioc_b$hzn='B'
clay_b$hzn='B'
silt_b$hzn='B'
sand_b$hzn='B'
b_thick$hzn='B'

soc_c$hzn='C'
tn_c$hzn='C'
cec_c$hzn='C'
ec_c$hzn='C'
ioc_c$hzn='C'
clay_c$hzn='C'
silt_c$hzn='C'
sand_c$hzn='C'
c_thick$hzn='C'

soc_all=rbind(soc_a, soc_b, soc_c)
tn_all=rbind(tn_a, tn_b, tn_c)
cec_all=rbind(cec_a, cec_b, cec_c)
ec_all=rbind(ec_a, ec_b, ec_c)
ioc_all=rbind(ioc_a, ioc_b, ioc_c)
clay_all=rbind(clay_a, clay_b, clay_c)
silt_all=rbind(silt_a, silt_b, silt_c)
sand_all=rbind(sand_a, sand_b, sand_c)
hzn_all=rbind(a_thick, b_thick, c_thick)
colnames(hzn_all)=c("PEDON_ID", 'hzn_thickness', 'hzn')

soc_all=merge(soc_all,raster_values, by='PEDON_ID')
tn_all=merge(tn_all,raster_values, by='PEDON_ID')
cec_all=merge(cec_all,raster_values, by='PEDON_ID')
ec_all=merge(ec_all,raster_values, by='PEDON_ID')
ioc_all=merge(ioc_all,raster_values, by='PEDON_ID')
clay_all=merge(clay_all,raster_values, by='PEDON_ID')
silt_all=merge(silt_all,raster_values, by='PEDON_ID')
sand_all=merge(sand_all,raster_values, by='PEDON_ID')
hzn_all=merge(hzn_all,raster_values, by='PEDON_ID')

soc_stock_all=merge(soc_stock, raster_values, by="PEDON_ID")


############Soil Organic Carbon Model#################
setwd('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Figures/3D')
'%!in%' <- function(x,y)!('%in%'(x,y))
soc_train_all=soc_all[soc_all$PEDON_ID %!in% test_pedons,]
soc_test_all=soc_all[soc_all$PEDON_ID %in% test_pedons,]

soc_train_all=soc_train_all[,-1]

soc_train_all[is.na(soc_train_all)] <- 0
soc_test_all[is.na(soc_test_all)] <- 0

soc_train_all=do.call(data.frame, lapply(soc_train_all, function(x) replace(x, is.infinite(x), 0)))
soc_test_all=do.call(data.frame, lapply(soc_test_all, function(x) replace(x, is.infinite(x), 0)))

soc_train_all=soc_train_all[soc_train_all$soc>0,]
soc_test_all=soc_test_all[soc_test_all$soc>0,]


#feature selection
carbon_features=psm_ranger_feature_selection(x="soc", soc_train_all[soc_train_all$hzn=='A',])
carbon_features

carbon_features=c(carbon_features, 'hzn')

#final model
soc_train_all_final=soc_train_all[,carbon_features]
model_carbon=ranger(soc~., data=soc_train_all_final, importance='impurity')
write.csv(sort(importance(model_carbon), decreasing=TRUE)/sum(importance(model_carbon)), "soc_all_feature.csv")

pred=predict(model_carbon, data=soc_test_all)
pred=pred$predictions

#model performance
cor.test(pred, soc_test_all$soc)
sd(soc_test_all$soc)/rmse(soc_test_all$soc, pred)
pred_null=rep(mean(soc_train_all$soc), length(soc_test_all$soc))


soc_performance=data.frame(

  t(data.frame(summary(lm(pred~soc_test_all$soc))$r.squared,
              rmse(soc_test_all$soc, pred),
              CCC(soc_test_all$soc, pred)$rho.c,
              bias(soc_test_all$soc, pred))),
  #performance by horizon
  
  t(data.frame(summary(lm(pred[soc_test_all$hzn=='A']~soc_test_all$soc[soc_test_all$hzn=='A']))$r.squared,
              rmse(soc_test_all$soc[soc_test_all$hzn=='A'], pred[soc_test_all$hzn=='A']),
              CCC(soc_test_all$soc[soc_test_all$hzn=='A'], pred[soc_test_all$hzn=='A'])$rho.c,
              bias(soc_test_all$soc[soc_test_all$hzn=='A'], pred[soc_test_all$hzn=='A']))),
  
  t(data.frame(summary(lm(pred[soc_test_all$hzn=='B']~soc_test_all$soc[soc_test_all$hzn=='B']))$r.squared,
              rmse(soc_test_all$soc[soc_test_all$hzn=='B'], pred[soc_test_all$hzn=='B']),
              CCC(soc_test_all$soc[soc_test_all$hzn=='B'], pred[soc_test_all$hzn=='B'])$rho.c,
              bias(soc_test_all$soc[soc_test_all$hzn=='B'], pred[soc_test_all$hzn=='B']))),
  
  t(data.frame(summary(lm(pred[soc_test_all$hzn=='C']~soc_test_all$soc[soc_test_all$hzn=='C']))$r.squared,
              rmse(soc_test_all$soc[soc_test_all$hzn=='C'], pred[soc_test_all$hzn=='C']),
              CCC(soc_test_all$soc[soc_test_all$hzn=='C'], pred[soc_test_all$hzn=='C'])$rho.c,
              bias(soc_test_all$soc[soc_test_all$hzn=='C'], pred[soc_test_all$hzn=='C']))),
  
  t(data.frame(rmse(soc_test_all$soc, pred_null),
               rmse(soc_test_all$soc[soc_test_all$hzn=='A'], pred_null[soc_test_all$hzn=='A']),
               rmse(soc_test_all$soc[soc_test_all$hzn=='B'], pred_null[soc_test_all$hzn=='B']),
               rmse(soc_test_all$soc[soc_test_all$hzn=='C'], pred_null[soc_test_all$hzn=='C'])))
  )


row.names(soc_performance) = c('r2', 'rmse', 'rhoc', 'rho.lwr', 'rho.high', 'bias')
colnames(soc_performance)=c('overall', 'A', 'B', 'C', 'null')

write.csv(soc_performance, 'soc_model_performance_metrics.csv')

#create plot
plot_data_soc=data.frame(soc_test_all$soc, pred, soc_test_all$hzn)
plot_data_soc$bare_soil=ifelse(soc_test_all$b4==0, 'N', "Y")
colnames(plot_data_soc)=c("actual", 'predicted','Horizon','Bare_Soil_Data')
soc_plot=ggplot(plot_data_soc, (aes(x=actual, y=predicted, shape=Horizon, color=Bare_Soil_Data))) + geom_point(size=4) + xlim(0, 6.5) + ylim(0, 6.5) + geom_abline() + xlab("Observed Soil Organic Carbon (%)") + ylab("Predicted Soil Organic Carbon (%)") + ggtitle("Soil Organic Carbon") + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))
#soc_plot = soc_plot + annotate("text", x=0.5, y = 6.5, label=expression(paste("R"^"2 ", "= 0.73"))) + annotate("text", x=0.5, y=6.2, label="RMSE = 0.60%") + annotate("text", x=0.5, y=5.9, label=expression(paste(rho['c'], "= 0.85"))) + annotate("text", x=0.5, y=5.6, label='Bias = -0.05') + update_geom_defaults("text", list(size = 6))
soc_plot

ggsave('soc_all.png', plot=last_plot(), width=11, height=8.5)


############Soil Organic Carbon Stock Model#################
setwd('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Figures/3D/')
'%!in%' <- function(x,y)!('%in%'(x,y))
soc_stock_all=soc_stock_all[soc_stock_all$soc_stock>0]
soc_stock_all=soc_stock_all[!is.na(soc_stock_all$soc_stock),]

stock_train_all=soc_stock_all[soc_stock_all$PEDON_ID %!in% test_pedons,]
stock_test_all=soc_stock_all[soc_stock_all$PEDON_ID %in% test_pedons,]

stock_train_all=stock_train_all[,-1]

stock_train_all[is.na(stock_train_all)] <- 0
stock_test_all[is.na(stock_test_all)] <- 0

stock_train_all=do.call(data.frame, lapply(stock_train_all, function(x) replace(x, is.infinite(x), 0)))
stock_test_all=do.call(data.frame, lapply(stock_test_all, function(x) replace(x, is.infinite(x), 0)))

#feature selection
stock__features=psm_ranger_feature_selection(x="soc_stock", stock_train_all)
stock__features

#final model
stock_train_all_final=stock_train_all[,stock__features]
model_carbon=ranger(soc_stock~., data=stock_train_all_final, importance='impurity')
write.csv(sort(importance(model_carbon), decreasing=TRUE)/sum(importance(model_carbon)), "stock_all_feature.csv")

pred=predict(model_carbon, data=stock_test_all)
pred=pred$predictions

summary(lm(pred~stock_test_all$soc_stock))$r.squared

#model performance
pred_null=rep(mean(stock_train_all$soc_stock), length(stock_test_all$soc_stock))


stock_performance=data.frame(
  
  t(data.frame(summary(lm(pred~stock_test_all$soc_stock))$r.squared,
               rmse(stock_test_all$soc_stock, pred),
               CCC(stock_test_all$soc_stock, pred)$rho.c,
               bias(stock_test_all$soc_stock, pred))),
              rmse(stock_test_all$soc_stock, pred_null))

row.names(stock_performance) = c('r2', 'rmse', 'rhoc', 'rho.lwr', 'rho.high', 'bias')
colnames(stock_performance)=c('overall',  'null')

write.csv(stock_performance, 'stock_model_performance_metrics.csv')


#create plot
plot_data_soc_stock=data.frame(stock_test_all$soc_stock, pred)
plot_data_soc_stock$bare_soil=ifelse(stock_test_all$b4==0, 'N', "Y")
colnames(plot_data_soc_stock)=c("actual", 'predicted', 'Bare_Soil_Data')
carbon_stock_plot=ggplot(plot_data_soc_stock, (aes(x=actual, y=predicted, shape=Bare_Soil_Data, color=Bare_Soil_Data))) + geom_point(size=4) + xlim(0,30) + ylim(0, 30) + geom_abline() +  xlab(expression(paste('Observed Soil Organic Carbon Stock (kg m'^"-2",")"))) + ylab(expression(paste('Predicted Soil Organic Carbon Stock (kg m'^"-2",")"))) +ggtitle("Soil Organic Carbon Stock") + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))
#carbon_stock_plot = carbon_stock_plot + annotate("text", x=3, y = 30, label=expression(paste("R"^"2 ", "= 0.27"))) + annotate("text", x=3, y=28, label=expression(paste("RMSE = 4.8 kg m"^"-2"))) + annotate("text", x=3, y=26, label=expression(paste(rho['c'], "= 0.47"))) + annotate("text", x=3, y=24, label='Bias = 0.31')+ update_geom_defaults("text", list(size = 6))
carbon_stock_plot

ggsave('soc_stock.png', plot=last_plot(), width=11, height=8.5)


############Total Nitrogen Model#################
setwd('/home/preston/OneDrive/Papers/Historic_Soil_Properties/Figures/3D')
'%!in%' <- function(x,y)!('%in%'(x,y))
#remove tn values above 0.5, 95th percentile is 0.43, and extreme values won't be reliably predictable
tn_all=tn_all[tn_all$tn<0.5,]

tn_train_all=tn_all[tn_all$PEDON_ID %!in% test_pedons,]
tn_test_all=tn_all[tn_all$PEDON_ID %in% test_pedons,]

tn_train_all=tn_train_all[,-1]

tn_train_all[is.na(tn_train_all)] <- 0
tn_test_all[is.na(tn_test_all)] <- 0

tn_train_all=do.call(data.frame, lapply(tn_train_all, function(x) replace(x, is.infinite(x), 0)))
tn_test_all=do.call(data.frame, lapply(tn_test_all, function(x) replace(x, is.infinite(x), 0)))


tn_train_all=tn_train_all[tn_train_all$tn>0,]
tn_test_all=tn_test_all[tn_test_all$tn>0,]

#feature selection
tn_features=psm_ranger_feature_selection(x="tn", tn_train_all[tn_train_all$hzn=='A',])

tn_features

tn_features=c(tn_features, 'hzn')

#final model
tn_train_all_final=tn_train_all[,tn_features]
model_tn=ranger(tn~., data=tn_train_all_final, importance='impurity')
write.csv(sort(importance(model_tn), decreasing=TRUE)/sum(importance(model_tn)), "tn_all_feature.csv")

pred=predict(model_tn, data=tn_test_all)
pred=pred$predictions

summary(lm(pred[tn_test_all$hzn=='A']~tn_test_all$tn[tn_test_all$hzn=='A']))$r.squared

#model performance
pred_null=rep(mean(tn_train_all$tn), length(tn_test_all$tn))


tn_performance=data.frame(
  
  t(data.frame(summary(lm(pred~tn_test_all$tn))$r.squared,
               rmse(tn_test_all$tn, pred),
               CCC(tn_test_all$tn, pred)$rho.c,
               bias(tn_test_all$tn, pred))),
  #performance by horizon
  
  t(data.frame(summary(lm(pred[tn_test_all$hzn=='A']~tn_test_all$tn[tn_test_all$hzn=='A']))$r.squared,
               rmse(tn_test_all$tn[tn_test_all$hzn=='A'], pred[tn_test_all$hzn=='A']),
               CCC(tn_test_all$tn[tn_test_all$hzn=='A'], pred[tn_test_all$hzn=='A'])$rho.c,
               bias(tn_test_all$tn[tn_test_all$hzn=='A'], pred[tn_test_all$hzn=='A']))),
  
  t(data.frame(summary(lm(pred[tn_test_all$hzn=='B']~tn_test_all$tn[tn_test_all$hzn=='B']))$r.squared,
               rmse(tn_test_all$tn[tn_test_all$hzn=='B'], pred[tn_test_all$hzn=='B']),
               CCC(tn_test_all$tn[tn_test_all$hzn=='B'], pred[tn_test_all$hzn=='B'])$rho.c,
               bias(tn_test_all$tn[tn_test_all$hzn=='B'], pred[tn_test_all$hzn=='B']))),
  
  t(data.frame(summary(lm(pred[tn_test_all$hzn=='C']~tn_test_all$tn[tn_test_all$hzn=='C']))$r.squared,
               rmse(tn_test_all$tn[tn_test_all$hzn=='C'], pred[tn_test_all$hzn=='C']),
               CCC(tn_test_all$tn[tn_test_all$hzn=='C'], pred[tn_test_all$hzn=='C'])$rho.c,
               bias(tn_test_all$tn[tn_test_all$hzn=='C'], pred[tn_test_all$hzn=='C']))),
  
  t(data.frame(summary(lm(pred_null~tn_test_all$tn))$r.squared,
               rmse(tn_test_all$tn, pred_null),
               CCC(tn_test_all$tn, pred_null)$rho.c,
               bias(tn_test_all$tn, pred_null)))
)


row.names(tn_performance) = c('r2', 'rmse', 'rhoc', 'rho.lwr', 'rho.high', 'bias')
colnames(tn_performance)=c('overall', 'A', 'B', 'C', 'null')

write.csv(tn_performance, 'tn_model_performance_metrics.csv')

#create plot
plot_data_tn=data.frame(tn_test_all$tn, pred, tn_test_all$hzn)
plot_data_tn$bare_soil=ifelse(tn_test_all$b4==0, 'N', "Y")
colnames(plot_data_tn)=c("actual", 'predicted', 'Horizon', 'Bare_Soil_Data')
tn_plot=ggplot(plot_data_tn, (aes(x=actual, y=predicted, shape=Horizon, color=Bare_Soil_Data))) + geom_point(size=4) + xlim(0, 0.5) + ylim(0, 0.5) + geom_abline() + xlab("Observed Total Nitrogen (%)") + ylab("Predicted Total Nitrogen (%)")+ggtitle("Total Nitrogen") + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))
#tn_plot = tn_plot + annotate("text", x=0.04, y = 0.5, label=expression(paste("R"^"2 ", "= 0.44"))) + annotate("text", x=0.04, y=0.47, label="RMSE = 0.08%") + annotate("text", x=0.04, y=0.44, label=expression(paste(rho['c'], "= 0.58"))) + annotate("text", x=0.04, y=0.41, label='Bias = 0.01')+ update_geom_defaults("text", list(size = 6))
tn_plot

ggsave('tn_all.png', plot=last_plot(), width=11, height=8.5)

############IOC#################
ioc_train=ioc_all[ioc_all$PEDON_ID %!in% test_pedons,]
ioc_test=ioc_all[ioc_all$PEDON_ID %in% test_pedons,]

ioc_train=ioc_train[!is.na(ioc_train$ioc),]
ioc_test=ioc_test[!is.na(ioc_test$ioc),]

ioc_train=ioc_train[,-1]

ioc_train[is.na(ioc_train)] <- 0
ioc_test[is.na(ioc_test)] <- 0

ioc_train=do.call(data.frame, lapply(ioc_train, function(x) replace(x, is.infinite(x), 0)))
ioc_test=do.call(data.frame, lapply(ioc_test, function(x) replace(x, is.infinite(x), 0)))

#feature selection
ioc_all_features=psm_ranger_feature_selection(x="ioc", ioc_train[ioc_train$hzn=='C',])
ioc_all_features

ioc_all_features=c(ioc_all_features, 'hzn')

#final model
ioc_train_final=ioc_train[,ioc_all_features]
model_ioc_all=ranger(ioc~., data=ioc_train_final, importance='impurity')
write.csv(sort(importance(model_ioc_all), decreasing=TRUE)/sum(importance(model_ioc_all)), "ioc_all_feature.csv")

pred=predict(model_ioc_all, data=ioc_test)
pred=pred$predictions

#model performance
pred_null=rep(mean(ioc_train$ioc), length(ioc_test$ioc))


ioc_performance=data.frame(
  
  t(data.frame(summary(lm(pred~ioc_test$ioc))$r.squared,
               rmse(ioc_test$ioc, pred),
               CCC(ioc_test$ioc, pred)$rho.c,
               bias(ioc_test$ioc, pred))),
  #performance by horizon
  
  t(data.frame(summary(lm(pred[ioc_test$hzn=='A']~ioc_test$ioc[ioc_test$hzn=='A']))$r.squared,
               rmse(ioc_test$ioc[ioc_test$hzn=='A'], pred[ioc_test$hzn=='A']),
               CCC(ioc_test$ioc[ioc_test$hzn=='A'], pred[ioc_test$hzn=='A'])$rho.c,
               bias(ioc_test$ioc[ioc_test$hzn=='A'], pred[ioc_test$hzn=='A']))),
  
  t(data.frame(summary(lm(pred[ioc_test$hzn=='B']~ioc_test$ioc[ioc_test$hzn=='B']))$r.squared,
               rmse(ioc_test$ioc[ioc_test$hzn=='B'], pred[ioc_test$hzn=='B']),
               CCC(ioc_test$ioc[ioc_test$hzn=='B'], pred[ioc_test$hzn=='B'])$rho.c,
               bias(ioc_test$ioc[ioc_test$hzn=='B'], pred[ioc_test$hzn=='B']))),
  
  t(data.frame(summary(lm(pred[ioc_test$hzn=='C']~ioc_test$ioc[ioc_test$hzn=='C']))$r.squared,
               rmse(ioc_test$ioc[ioc_test$hzn=='C'], pred[ioc_test$hzn=='C']),
               CCC(ioc_test$ioc[ioc_test$hzn=='C'], pred[ioc_test$hzn=='C'])$rho.c,
               bias(ioc_test$ioc[ioc_test$hzn=='C'], pred[ioc_test$hzn=='C']))),
  
  t(data.frame(summary(lm(pred_null~ioc_test$ioc))$r.squared,
               rmse(ioc_test$ioc, pred_null),
               CCC(ioc_test$ioc, pred_null)$rho.c,
               bias(ioc_test$ioc, pred_null)))
)


row.names(ioc_performance) = c('r2', 'rmse', 'rhoc', 'rho.lwr', 'rho.high', 'bias')
colnames(ioc_performance)=c('overall', 'A', 'B', 'C', 'null')

write.csv(ioc_performance, 'ioc_model_performance_metrics.csv')


#create plot
plot_data_ioc_all=data.frame(ioc_test$ioc, pred, ioc_test$hzn)
plot_data_ioc_all$bare_soil=ifelse(ioc_test$b4==0, 'N', "Y")
colnames(plot_data_ioc_all)=c("actual", 'predicted', 'Horizon', 'Bare_Soil_Data')
ioc_all_plot=ggplot(plot_data_ioc_all, (aes(x=actual, y=predicted, shape=Horizon, color=Bare_Soil_Data))) + geom_point(size=4) + xlim(0,30) + ylim(0, 30) + geom_abline() +  xlab('Observed A Horizon Inorganic Carbon (%)') + ylab('Predicted A Horizon Inorganic Carbon (%)') + ggtitle("Inorganic Carbon") + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))
#ioc_all_plot = ioc_all_plot + annotate("text", x=2, y = 30, label=expression(paste("R"^"2 ", "= 0.07"))) + annotate("text", x=2, y=28, label=expression(paste("RMSE = 3.1%"))) + annotate("text", x=2, y=26, label=expression(paste(rho['c'], "= 0.22"))) + annotate("text", x=2, y=24, label='Bias = -0.44')+ update_geom_defaults("text", list(size = 6))
ioc_all_plot

ggsave('ioc_all.png', plot=last_plot(), width=11, height=8.5)

##############CEC#################
cec_train_all=cec_all[cec_all$PEDON_ID %!in% test_pedons,]
cec_test_all=cec_all[cec_all$PEDON_ID %in% test_pedons,]

cec_train_all=cec_train_all[!is.na(cec_train_all$cec),]
cec_test_all=cec_test_all[!is.na(cec_test_all$cec),]

cec_train_all=cec_train_all[,-1]

cec_train_all[is.na(cec_train_all)] <- 0
cec_test_all[is.na(cec_test_all)] <- 0

#feature selection with ranger
cec_features=psm_ranger_feature_selection(x="cec", cec_train_all[cec_train_all$hzn=='A',])
cec_features

cec_features=c(cec_features, 'hzn')

#final model
cec_train_all_final=cec_train_all[,cec_features]

model_cec=ranger(cec~., data=cec_train_all_final, importance='impurity')

sort(importance(model_cec), decreasing=TRUE)
write.csv(sort(importance(model_cec), decreasing=TRUE)/sum(importance(model_cec)), "cec_all_feature.csv")

pred=predict(model_cec, data=cec_test_all)
pred=pred$predictions


#model performance
pred_null=rep(mean(cec_train_all$cec), length(cec_test_all$cec))


cec_performance=data.frame(
  
  t(data.frame(summary(lm(pred~cec_test_all$cec))$r.squared,
               rmse(cec_test_all$cec, pred),
               CCC(cec_test_all$cec, pred)$rho.c,
               bias(cec_test_all$cec, pred))),
  #performance by horizon
  
  t(data.frame(summary(lm(pred[cec_test_all$hzn=='A']~cec_test_all$cec[cec_test_all$hzn=='A']))$r.squared,
               rmse(cec_test_all$cec[cec_test_all$hzn=='A'], pred[cec_test_all$hzn=='A']),
               CCC(cec_test_all$cec[cec_test_all$hzn=='A'], pred[cec_test_all$hzn=='A'])$rho.c,
               bias(cec_test_all$cec[cec_test_all$hzn=='A'], pred[cec_test_all$hzn=='A']))),
  
  t(data.frame(summary(lm(pred[cec_test_all$hzn=='B']~cec_test_all$cec[cec_test_all$hzn=='B']))$r.squared,
               rmse(cec_test_all$cec[cec_test_all$hzn=='B'], pred[cec_test_all$hzn=='B']),
               CCC(cec_test_all$cec[cec_test_all$hzn=='B'], pred[cec_test_all$hzn=='B'])$rho.c,
               bias(cec_test_all$cec[cec_test_all$hzn=='B'], pred[cec_test_all$hzn=='B']))),
  
  t(data.frame(summary(lm(pred[cec_test_all$hzn=='C']~cec_test_all$cec[cec_test_all$hzn=='C']))$r.squared,
               rmse(cec_test_all$cec[cec_test_all$hzn=='C'], pred[cec_test_all$hzn=='C']),
               CCC(cec_test_all$cec[cec_test_all$hzn=='C'], pred[cec_test_all$hzn=='C'])$rho.c,
               bias(cec_test_all$cec[cec_test_all$hzn=='C'], pred[cec_test_all$hzn=='C']))),
  
  t(data.frame(summary(lm(pred_null~cec_test_all$cec))$r.squared,
               rmse(cec_test_all$cec, pred_null),
               CCC(cec_test_all$cec, pred_null)$rho.c,
               bias(cec_test_all$cec, pred_null)))
)


row.names(cec_performance) = c('r2', 'rmse', 'rhoc', 'rho.lwr', 'rho.high', 'bias')
colnames(cec_performance)=c('overall', 'A', 'B', 'C', 'null')

write.csv(cec_performance, 'cec_model_performance_metrics.csv')

#create plot
plot_data_cec=data.frame(cec_test_all$cec, pred, cec_test_all$hzn)
plot_data_cec$bare_soil=ifelse(cec_test_all$b4==0, 'N', "Y")
colnames(plot_data_cec)=c("actual", 'predicted', 'Horizon', 'Bare_Soil_Data')
cec_plot=ggplot(plot_data_cec, (aes(x=actual, y=predicted, shape=Horizon, color=Bare_Soil_Data))) + geom_point(size=4) + xlim(0, 60) + ylim(0, 60) + geom_abline() + xlab(expression(paste('Observed Cation Exchange Capacity (meq 100g'^"-1",")"))) + ylab(expression(paste('Predicted Cation Exchange Capacity (meq 100g'^"-1",")")))+ggtitle("Cation Exchange Capacity") + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))
#cec_plot = cec_plot + annotate("text", x=5, y = 60, label=expression(paste("R"^"2 ", "= 0.48"))) + annotate("text", x=5, y=57, label=expression(paste("RMSE = 8.8 meq 100g"^"-1"))) + annotate("text", x=5, y=54, label=expression(paste(rho['c'], "= 0.62"))) + annotate("text", x=5, y=51, label='Bias = -0.71')+ update_geom_defaults("text", list(size = 6))
cec_plot

ggsave('cec_all.png', plot=last_plot(), width=11, height=8.5)

#############EC###################
ec_train_all=ec_all[ec_all$PEDON_ID %!in% test_pedons,]
ec_test_all=ec_all[ec_all$PEDON_ID %in% test_pedons,]


ec_train_all=ec_train_all[!is.na(ec_train_all$ec),]
ec_test_all=ec_test_all[!is.na(ec_test_all$ec),]

ec_train_all=ec_train_all[,-1]

ec_train_all[is.na(ec_train_all)] <- 0
ec_test_all[is.na(ec_test_all)] <- 0

#feature selection with ranger
ec_all_features=psm_ranger_feature_selection(x="ec", ec_train_all[ec_train_all$hzn=='C',])
ec_all_features

ec_all_features=c(ec_all_features, 'hzn')

#final model
ec_train_all_final=ec_train_all[,ec_all_features]
model_ec_all=ranger(ec~., data=ec_train_all_final, importance='impurity')
sort(importance(model_ec_all), decreasing=TRUE)
write.csv(sort(importance(model_ec_all), decreasing=TRUE)/sum(importance(model_ec_all)), "ec_all_feature.csv")

pred=predict(model_ec_all, data=ec_test_all)
pred=pred$predictions

#model performance
pred_null=rep(mean(ec_train_all$ec), length(ec_test_all$ec))


ec_performance=data.frame(
  
  t(data.frame(summary(lm(pred~ec_test_all$ec))$r.squared,
               rmse(ec_test_all$ec, pred),
               CCC(ec_test_all$ec, pred)$rho.c,
               bias(ec_test_all$ec, pred))),
  #performance by horizon
  
  t(data.frame(summary(lm(pred[ec_test_all$hzn=='A']~ec_test_all$ec[ec_test_all$hzn=='A']))$r.squared,
               rmse(ec_test_all$ec[ec_test_all$hzn=='A'], pred[ec_test_all$hzn=='A']),
               CCC(ec_test_all$ec[ec_test_all$hzn=='A'], pred[ec_test_all$hzn=='A'])$rho.c,
               bias(ec_test_all$ec[ec_test_all$hzn=='A'], pred[ec_test_all$hzn=='A']))),
  
  t(data.frame(summary(lm(pred[ec_test_all$hzn=='B']~ec_test_all$ec[ec_test_all$hzn=='B']))$r.squared,
               rmse(ec_test_all$ec[ec_test_all$hzn=='B'], pred[ec_test_all$hzn=='B']),
               CCC(ec_test_all$ec[ec_test_all$hzn=='B'], pred[ec_test_all$hzn=='B'])$rho.c,
               bias(ec_test_all$ec[ec_test_all$hzn=='B'], pred[ec_test_all$hzn=='B']))),
  
  t(data.frame(summary(lm(pred[ec_test_all$hzn=='C']~ec_test_all$ec[ec_test_all$hzn=='C']))$r.squared,
               rmse(ec_test_all$ec[ec_test_all$hzn=='C'], pred[ec_test_all$hzn=='C']),
               CCC(ec_test_all$ec[ec_test_all$hzn=='C'], pred[ec_test_all$hzn=='C'])$rho.c,
               bias(ec_test_all$ec[ec_test_all$hzn=='C'], pred[ec_test_all$hzn=='C']))),
  
  t(data.frame(summary(lm(pred_null~ec_test_all$ec))$r.squared,
               rmse(ec_test_all$ec, pred_null),
               CCC(ec_test_all$ec, pred_null)$rho.c,
               bias(ec_test_all$ec, pred_null)))
)


row.names(ec_performance) = c('r2', 'rmse', 'rhoc', 'rho.lwr', 'rho.high', 'bias')
colnames(ec_performance)=c('overall', 'A', 'B', 'C', 'null')

write.csv(ec_performance, 'ec_model_performance_metrics.csv')

#create plot
plot_data_ec_all=data.frame(ec_test_all$ec, pred, ec_test_all$hzn)
plot_data_ec_all$bare_soil=ifelse(ec_test_all$b4==0, 'N', "Y")
colnames(plot_data_ec_all)=c("actual", 'predicted', 'Horizon', 'Bare_Soil_Data')
ec_plot_all=ggplot(plot_data_ec_all, (aes(x=actual, y=predicted, shape=Horizon, color=Bare_Soil_Data))) + geom_point(size=4) + xlim(0, 18) + ylim(0, 18) + geom_abline() + xlab(expression(paste('Observed A Horizon Electrical Conductivity (dS m'^"-1",")"))) + ylab(expression(paste('Predicted A Horizon Electrical Conductivity (dS m'^"-1",")")))+ggtitle("Electrical Conductivity") + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))
#ec_plot_all = ec_plot_all + annotate("text", x=1.5, y = 18, label=expression(paste("R"^"2 ", "= 0.09"))) + annotate("text", x=1.5, y=17, label=expression(paste("RMSE = 1.6 dS m"^"-1"))) + annotate("text", x=1.5, y=16, label=expression(paste(rho['c'], "= 0.14"))) + annotate("text", x=1.5, y=15, label='Bias = 0.01')+ update_geom_defaults("text", list(size = 6))
ec_plot_all

ggsave('ec_all.png', plot=last_plot(), width=11, height=8.5)

##############Clay#################
clay_train_all=clay_all[clay_all$PEDON_ID %!in% test_pedons,]
clay_test_all=clay_all[clay_all$PEDON_ID %in% test_pedons,]

clay_train_all=clay_train_all[!is.na(clay_train_all$clay),]
clay_test_all=clay_test_all[!is.na(clay_test_all$clay),]


clay_train_all=clay_train_all[,-1]

clay_train_all[is.na(clay_train_all)] <- 0
clay_test_all[is.na(clay_test_all)] <- 0

#feature selection with ranger
clay_all_features=psm_ranger_feature_selection(x="clay", clay_train_all[clay_train_all$hzn=='A',])
clay_all_features

clay_all_features=c(clay_all_features, 'hzn')

#final model
clay_train_all_final=clay_train_all[,clay_all_features]

#final model
model_clay=ranger(clay~., data=clay_train_all_final, importance='impurity')
sort(importance(model_clay), decreasing=TRUE)
write.csv(sort(importance(model_clay), decreasing=TRUE)/sum(importance(model_clay)), "clay_all_feature.csv")

pred=predict(model_clay, data=clay_test_all)
pred=pred$predictions

summary(lm(pred[clay_test_all$hzn=='A']~clay_test_all$clay[clay_test_all$hzn=='A']))$r.squared

#model performance
pred_null=rep(mean(clay_train_all$clay), length(clay_test_all$clay))


clay_performance=data.frame(
  
  t(data.frame(summary(lm(pred~clay_test_all$clay))$r.squared,
               rmse(clay_test_all$clay, pred),
               CCC(clay_test_all$clay, pred)$rho.c,
               bias(clay_test_all$clay, pred))),
  #performance by horizon
  
  t(data.frame(summary(lm(pred[clay_test_all$hzn=='A']~clay_test_all$clay[clay_test_all$hzn=='A']))$r.squared,
               rmse(clay_test_all$clay[clay_test_all$hzn=='A'], pred[clay_test_all$hzn=='A']),
               CCC(clay_test_all$clay[clay_test_all$hzn=='A'], pred[clay_test_all$hzn=='A'])$rho.c,
               bias(clay_test_all$clay[clay_test_all$hzn=='A'], pred[clay_test_all$hzn=='A']))),
  
  t(data.frame(summary(lm(pred[clay_test_all$hzn=='B']~clay_test_all$clay[clay_test_all$hzn=='B']))$r.squared,
               rmse(clay_test_all$clay[clay_test_all$hzn=='B'], pred[clay_test_all$hzn=='B']),
               CCC(clay_test_all$clay[clay_test_all$hzn=='B'], pred[clay_test_all$hzn=='B'])$rho.c,
               bias(clay_test_all$clay[clay_test_all$hzn=='B'], pred[clay_test_all$hzn=='B']))),
  
  t(data.frame(summary(lm(pred[clay_test_all$hzn=='C']~clay_test_all$clay[clay_test_all$hzn=='C']))$r.squared,
               rmse(clay_test_all$clay[clay_test_all$hzn=='C'], pred[clay_test_all$hzn=='C']),
               CCC(clay_test_all$clay[clay_test_all$hzn=='C'], pred[clay_test_all$hzn=='C'])$rho.c,
               bias(clay_test_all$clay[clay_test_all$hzn=='C'], pred[clay_test_all$hzn=='C']))),
  
  t(data.frame(summary(lm(pred_null~clay_test_all$clay))$r.squared,
               rmse(clay_test_all$clay, pred_null),
               CCC(clay_test_all$clay, pred_null)$rho.c,
               bias(clay_test_all$clay, pred_null)))
)


row.names(clay_performance) = c('r2', 'rmse', 'rhoc', 'rho.lwr', 'rho.high', 'bias')
colnames(clay_performance)=c('overall', 'A', 'B', 'C', 'null')

write.csv(clay_performance, 'clay_model_performance_metrics.csv')

#create plot
plot_data_clay_all=data.frame(clay_test_all$clay, pred, clay_test_all$hzn)
plot_data_clay_all$bare_soil=ifelse(clay_test_all$b4==0, 'N', "Y")
colnames(plot_data_clay_all)=c("actual", 'predicted', 'Horizon', 'Bare_Soil_Data')
clay_plot_all= ggplot(plot_data_clay_all, (aes(x=actual, y=predicted, shape=Horizon, color=Bare_Soil_Data))) + geom_point(size=4) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Clay Content (%)'))) + ylab(expression(paste('Predicted Clay Content (%)')))+ggtitle("Clay Content") + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))
#clay_plot_all = clay_plot_all + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.68"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 8.0%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.76"))) + annotate("text", x=5, y=70, label='Bias = 0.04')
clay_plot_all

ggsave('clay_all.png', plot=last_plot(), width=11, height=8.5)


##############Silt#################
silt_train_all=silt_all[silt_all$PEDON_ID %!in% test_pedons,]
silt_test_all=silt_all[silt_all$PEDON_ID %in% test_pedons,]

silt_train_all=silt_train_all[!is.na(silt_train_all$silt),]
silt_test_all=silt_test_all[!is.na(silt_test_all$silt),]


silt_train_all=silt_train_all[,-1]

silt_train_all[is.na(silt_train_all)] <- 0
silt_test_all[is.na(silt_test_all)] <- 0

#feature selection with ranger
silt_all_features=psm_ranger_feature_selection(x="silt", silt_train_all[silt_train_all$hzn=='A',])
silt_all_features

silt_all_features=c(silt_all_features, 'hzn')

#final model
silt_train_all_final=silt_train_all[,silt_all_features]

#final model
model_silt=ranger(silt~., data=silt_train_all_final, importance='impurity')
sort(importance(model_silt), decreasing=TRUE)
write.csv(sort(importance(model_silt), decreasing=TRUE)/sum(importance(model_silt)), "silt_all_feature.csv")

pred=predict(model_silt, data=silt_test_all)
pred=pred$predictions

summary(lm(pred[silt_test_all$hzn=='A']~silt_test_all$silt[silt_test_all$hzn=='A']))$r.squared

#model performance
pred_null=rep(mean(silt_train_all$silt), length(silt_test_all$silt))


silt_performance=data.frame(
  
  t(data.frame(summary(lm(pred~silt_test_all$silt))$r.squared,
               rmse(silt_test_all$silt, pred),
               CCC(silt_test_all$silt, pred)$rho.c,
               bias(silt_test_all$silt, pred))),
  #performance by horizon
  
  t(data.frame(summary(lm(pred[silt_test_all$hzn=='A']~silt_test_all$silt[silt_test_all$hzn=='A']))$r.squared,
               rmse(silt_test_all$silt[silt_test_all$hzn=='A'], pred[silt_test_all$hzn=='A']),
               CCC(silt_test_all$silt[silt_test_all$hzn=='A'], pred[silt_test_all$hzn=='A'])$rho.c,
               bias(silt_test_all$silt[silt_test_all$hzn=='A'], pred[silt_test_all$hzn=='A']))),
  
  t(data.frame(summary(lm(pred[silt_test_all$hzn=='B']~silt_test_all$silt[silt_test_all$hzn=='B']))$r.squared,
               rmse(silt_test_all$silt[silt_test_all$hzn=='B'], pred[silt_test_all$hzn=='B']),
               CCC(silt_test_all$silt[silt_test_all$hzn=='B'], pred[silt_test_all$hzn=='B'])$rho.c,
               bias(silt_test_all$silt[silt_test_all$hzn=='B'], pred[silt_test_all$hzn=='B']))),
  
  t(data.frame(summary(lm(pred[silt_test_all$hzn=='C']~silt_test_all$silt[silt_test_all$hzn=='C']))$r.squared,
               rmse(silt_test_all$silt[silt_test_all$hzn=='C'], pred[silt_test_all$hzn=='C']),
               CCC(silt_test_all$silt[silt_test_all$hzn=='C'], pred[silt_test_all$hzn=='C'])$rho.c,
               bias(silt_test_all$silt[silt_test_all$hzn=='C'], pred[silt_test_all$hzn=='C']))),
  
  t(data.frame(summary(lm(pred_null~silt_test_all$silt))$r.squared,
               rmse(silt_test_all$silt, pred_null),
               CCC(silt_test_all$silt, pred_null)$rho.c,
               bias(silt_test_all$silt, pred_null)))
)


row.names(silt_performance) = c('r2', 'rmse', 'rhoc', 'rho.lwr', 'rho.high', 'bias')
colnames(silt_performance)=c('overall', 'A', 'B', 'C', 'null')

write.csv(silt_performance, 'silt_model_performance_metrics.csv')

#create plot
plot_data_silt_all=data.frame(silt_test_all$silt, pred, silt_test_all$hzn)
plot_data_silt_all$bare_soil=ifelse(silt_test_all$b4==0, 'N', "Y")
colnames(plot_data_silt_all)=c("actual", 'predicted', 'Horizon', 'Bare_Soil_Data')
silt_plot_all= ggplot(plot_data_silt_all, (aes(x=actual, y=predicted, shape=Horizon, color=Bare_Soil_Data))) + geom_point(size=4) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Silt Content (%)'))) + ylab(expression(paste('Predicted Silt Content (%)')))+ggtitle("Silt Content") + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))
#silt_plot_all = silt_plot_all + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.68"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 8.0%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.76"))) + annotate("text", x=5, y=70, label='Bias = 0.04')
silt_plot_all

ggsave('silt_all.png', plot=last_plot(), width=11, height=8.5)


##############Sand#################
sand_train_all=sand_all[sand_all$PEDON_ID %!in% test_pedons,]
sand_test_all=sand_all[sand_all$PEDON_ID %in% test_pedons,]

sand_train_all=sand_train_all[!is.na(sand_train_all$sand),]
sand_test_all=sand_test_all[!is.na(sand_test_all$sand),]


sand_train_all=sand_train_all[,-1]

sand_train_all[is.na(sand_train_all)] <- 0
sand_test_all[is.na(sand_test_all)] <- 0

#feature selection with ranger
sand_all_features=psm_ranger_feature_selection(x="sand", sand_train_all[sand_train_all$hzn=='A',])
sand_all_features

sand_all_features=c(sand_all_features, 'hzn')

#final model
sand_train_all_final=sand_train_all[,sand_all_features]

#final model
model_sand=ranger(sand~., data=sand_train_all_final, importance='impurity')
sort(importance(model_sand), decreasing=TRUE)
write.csv(sort(importance(model_sand), decreasing=TRUE)/sum(importance(model_sand)), "sand_all_feature.csv")

pred=predict(model_sand, data=sand_test_all)
pred=pred$predictions

summary(lm(pred[sand_test_all$hzn=='A']~sand_test_all$sand[sand_test_all$hzn=='A']))$r.squared

#model performance
pred_null=rep(mean(sand_train_all$sand), length(sand_test_all$sand))


sand_performance=data.frame(
  
  t(data.frame(summary(lm(pred~sand_test_all$sand))$r.squared,
               rmse(sand_test_all$sand, pred),
               CCC(sand_test_all$sand, pred)$rho.c,
               bias(sand_test_all$sand, pred))),
  #performance by horizon
  
  t(data.frame(summary(lm(pred[sand_test_all$hzn=='A']~sand_test_all$sand[sand_test_all$hzn=='A']))$r.squared,
               rmse(sand_test_all$sand[sand_test_all$hzn=='A'], pred[sand_test_all$hzn=='A']),
               CCC(sand_test_all$sand[sand_test_all$hzn=='A'], pred[sand_test_all$hzn=='A'])$rho.c,
               bias(sand_test_all$sand[sand_test_all$hzn=='A'], pred[sand_test_all$hzn=='A']))),
  
  t(data.frame(summary(lm(pred[sand_test_all$hzn=='B']~sand_test_all$sand[sand_test_all$hzn=='B']))$r.squared,
               rmse(sand_test_all$sand[sand_test_all$hzn=='B'], pred[sand_test_all$hzn=='B']),
               CCC(sand_test_all$sand[sand_test_all$hzn=='B'], pred[sand_test_all$hzn=='B'])$rho.c,
               bias(sand_test_all$sand[sand_test_all$hzn=='B'], pred[sand_test_all$hzn=='B']))),
  
  t(data.frame(summary(lm(pred[sand_test_all$hzn=='C']~sand_test_all$sand[sand_test_all$hzn=='C']))$r.squared,
               rmse(sand_test_all$sand[sand_test_all$hzn=='C'], pred[sand_test_all$hzn=='C']),
               CCC(sand_test_all$sand[sand_test_all$hzn=='C'], pred[sand_test_all$hzn=='C'])$rho.c,
               bias(sand_test_all$sand[sand_test_all$hzn=='C'], pred[sand_test_all$hzn=='C']))),
  
  t(data.frame(summary(lm(pred_null~sand_test_all$sand))$r.squared,
               rmse(sand_test_all$sand, pred_null),
               CCC(sand_test_all$sand, pred_null)$rho.c,
               bias(sand_test_all$sand, pred_null)))
)


row.names(sand_performance) = c('r2', 'rmse', 'rhoc', 'rho.lwr', 'rho.high', 'bias')
colnames(sand_performance)=c('overall', 'A', 'B', 'C', 'null')

write.csv(sand_performance, 'sand_model_performance_metrics.csv')

#create plot
plot_data_sand_all=data.frame(sand_test_all$sand, pred, sand_test_all$hzn)
plot_data_sand_all$bare_soil=ifelse(sand_test_all$b4==0, 'N', "Y")
colnames(plot_data_sand_all)=c("actual", 'predicted', 'Horizon', 'Bare_Soil_Data')
sand_plot_all= ggplot(plot_data_sand_all, (aes(x=actual, y=predicted, shape=Horizon, color=Bare_Soil_Data))) + geom_point(size=4) + geom_point() + xlim(0, 85) + ylim(0, 85) + geom_abline() + xlab(expression(paste('Observed Sand Content (%)'))) + ylab(expression(paste('Predicted Sand Content (%)')))+ggtitle("Sand Content") + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))
#sand_plot_all = sand_plot_all + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.68"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 8.0%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.76"))) + annotate("text", x=5, y=70, label='Bias = 0.04')
sand_plot_all

ggsave('sand_all.png', plot=last_plot(), width=11, height=8.5)

######################Horizon Thickness###################
hzn_train_all=hzn_all[hzn_all$PEDON_ID %!in% test_pedons,]
hzn_test_all=hzn_all[hzn_all$PEDON_ID %in% test_pedons,]

hzn_train_all=hzn_train_all[!is.na(hzn_train_all$hzn),]
hzn_test_all=hzn_test_all[!is.na(hzn_test_all$hzn),]


hzn_train_all=hzn_train_all[,-1]

hzn_train_all[is.na(hzn_train_all)] <- 0
hzn_test_all[is.na(hzn_test_all)] <- 0

#feature selection with ranger
hzn_all_features=psm_ranger_feature_selection(x="hzn_thickness", hzn_train_all[hzn_train_all$hzn=='C',])
hzn_all_features

hzn_all_features=c(hzn_all_features, 'hzn')

#final model
hzn_train_all_final=hzn_train_all[,hzn_all_features]

#final model
model_hzn=ranger(hzn_thickness~., data=hzn_train_all_final, importance='impurity')
sort(importance(model_hzn), decreasing=TRUE)
write.csv(sort(importance(model_hzn), decreasing=TRUE)/sum(importance(model_hzn)), "hzn_all_feature.csv")

pred=predict(model_hzn, data=hzn_test_all)
pred=pred$predictions

summary(lm(pred[hzn_test_all$hzn=='A']~hzn_test_all$hzn_thickness[hzn_test_all$hzn=='A']))$r.squared

#model performance
pred_null=rep(mean(hzn_train_all$hzn_thickness), length(hzn_test_all$hzn))


hzn_performance=data.frame(
  
  t(data.frame(summary(lm(pred~hzn_test_all$hzn_thickness))$r.squared,
               rmse(hzn_test_all$hzn_thickness, pred),
               CCC(hzn_test_all$hzn_thickness, pred)$rho.c,
               bias(hzn_test_all$hzn_thickness, pred))),
  #performance by horizon
  
  t(data.frame(summary(lm(pred[hzn_test_all$hzn=='A']~hzn_test_all$hzn_thickness[hzn_test_all$hzn=='A']))$r.squared,
               rmse(hzn_test_all$hzn_thickness[hzn_test_all$hzn=='A'], pred[hzn_test_all$hzn=='A']),
               CCC(hzn_test_all$hzn_thickness[hzn_test_all$hzn=='A'], pred[hzn_test_all$hzn=='A'])$rho.c,
               bias(hzn_test_all$hzn_thickness[hzn_test_all$hzn=='A'], pred[hzn_test_all$hzn=='A']))),
  
  t(data.frame(summary(lm(pred[hzn_test_all$hzn=='B']~hzn_test_all$hzn_thickness[hzn_test_all$hzn=='B']))$r.squared,
               rmse(hzn_test_all$hzn_thickness[hzn_test_all$hzn=='B'], pred[hzn_test_all$hzn=='B']),
               CCC(hzn_test_all$hzn_thickness[hzn_test_all$hzn=='B'], pred[hzn_test_all$hzn=='B'])$rho.c,
               bias(hzn_test_all$hzn_thickness[hzn_test_all$hzn=='B'], pred[hzn_test_all$hzn=='B']))),
  
  t(data.frame(summary(lm(pred[hzn_test_all$hzn=="C"]~hzn_test_all$hzn_thickness[hzn_test_all$hzn=="C"]))$r.squared,
               rmse(hzn_test_all$hzn_thickness[hzn_test_all$hzn=="C"], pred[hzn_test_all$hzn=="C"]),
               CCC(hzn_test_all$hzn_thickness[hzn_test_all$hzn=="C"], pred[hzn_test_all$hzn=="C"])$rho.c,
               bias(hzn_test_all$hzn_thickness[hzn_test_all$hzn=="C"], pred[hzn_test_all$hzn=="C"]))),
  

  t(data.frame(summary(lm(pred_null~hzn_test_all$hzn_thickness))$r.squared,
               rmse(hzn_test_all$hzn_thickness, pred_null),
               CCC(hzn_test_all$hzn_thickness, pred_null)$rho.c,
               bias(hzn_test_all$hzn_thickness, pred_null)))
)


row.names(hzn_performance) = c('r2', 'rmse', 'rhoc', 'rho.lwr', 'rho.high', 'bias')
colnames(hzn_performance)=c('overall', 'A', 'B', 'C', 'null')

write.csv(hzn_performance, 'hzn_model_performance_metrics.csv')

#create plot
plot_data_hzn_all=data.frame(hzn_test_all$hzn_thickness, pred, hzn_test_all$hzn)
plot_data_hzn_all$bare_soil=ifelse(hzn_test_all$b4==0, 'N', "Y")
colnames(plot_data_hzn_all)=c("actual", 'predicted', 'Horizon', 'Bare_Soil_Data')
plot_data_hzn_all[,1:2]=plot_data_hzn_all[,1:2]*100
hzn_plot_all= ggplot(plot_data_hzn_all, (aes(x=actual, y=predicted, shape=Horizon, color=Bare_Soil_Data))) + geom_point(size=4) + geom_point() + xlim(0, 100) + ylim(0, 100) + geom_abline() + xlab(expression(paste('Observed Horizon Thickness (%)'))) + ylab(expression(paste('Predicted Horizon Thickness (%)')))+ggtitle("Horizon Thickness") + theme(axis.text=element_text(size=20), axis.title=element_text(size=20), title=element_text(size=20),legend.key.size = unit(1, 'cm'), legend.text=element_text(size=20))
#hzn_plot_all = hzn_plot_all + annotate("text", x=5, y = 85, label=expression(paste("R"^"2 ", "= 0.68"))) + annotate("text", x=5, y=80, label=expression(paste("RMSE = 8.0%"))) + annotate("text", x=5, y=75, label=expression(paste(rho['c'], "= 0.76"))) + annotate("text", x=5, y=70, label='Bias = 0.04')
hzn_plot_all

ggsave('hzn_all.png', plot=last_plot(), width=11, height=8.5)

####################SOC Stock###################
soc_stock_test_all=hzn_test_all

#bulk densities
pred_soc=predict(model_carbon, data=hzn_test_all)$predictions
pred_clay=predict(model_clay, data=hzn_test_all)$predictions
pred_silt=predict(model_silt, data=hzn_test_all)$predictions
pred_sand=predict(model_sand, data=hzn_test_all)$predictions

pred_bd_test=data.frame(pred_soc,pred_clay, pred_silt, pred_sand)
colnames(pred_bd_test)=colnames(bd_input)

pred_bd_stock=predict(model_bd, pred_bd_test)$predictions


#thickness
pred_thick=predict(model_hzn, data=hzn_test_all)$predictions
hzn_test_all

pred_soc_stock=data.frame()
for (i in unique(hzn_test_all$PEDON_ID)){
  temp_1=temp[temp$PEDON_ID==i,]
  temp_1$pred_thick=temp_1$pred_thick/sum(temp_1$pred_thick)*100
  stock_temp=sum(round(temp_1$pred_soc/100*(temp_1$pred_thick)*temp_1$pred_bd_stock, 2)*10)
  temp_1=data.frame(i, stock_temp)
  pred_soc_stock=rbind(pred_soc_stock, temp_1)
}


colnames(pred_soc_stock)=c("PEDON_ID", 'pred_soc_stock')

pred_soc_stock=merge(soc_stock, pred_soc_stock, by='PEDON_ID')

summary(lm(pred_soc_stock$pred_soc_stock~pred_soc_stock$soc_stock))


#no bd transfer in predictions
library(reshape2)
stock_ptf=data.frame(soc_all$PEDON_ID, soc_all$soc, clay_all$clay, sand_all$sand, hzn_all$hzn_thickness, soc_all$hzn)
colnames(stock_ptf)=c("PEDON_ID", 'soc', 'clay', 'sand', 'hzn_thickness', 'hzn')
stock_ptf=reshape(stock_ptf, idvar='PEDON_ID', timevar='hzn', direction='wide')
stock_ptf=na.omit(stock_ptf)

stock_ptf=merge(soc_stock, stock_ptf, by='PEDON_ID')


stock_ptf_train=stock_ptf[stock_ptf$PEDON_ID %!in% test_pedons,]
stock_ptf_test=stock_ptf[stock_ptf$PEDON_ID %in% test_pedons,]

stock_ptf_train=na.omit(stock_ptf_train)

model_stock_ptf=ranger(soc_stock~., data=stock_ptf_train, importance='impurity')

#create prediction data
pred_soc=predict(model_carbon, data=hzn_test_all)$predictions
pred_clay=predict(model_clay, data=hzn_test_all)$predictions
pred_silt=predict(model_silt, data=hzn_test_all)$predictions
pred_thick=predict(model_hzn, data=hzn_test_all)$predictions

stock_ptf_test_input=data.frame(hzn_test_all$PEDON_ID, pred_soc, pred_clay, pred_silt, pred_thick, hzn_test_all$hzn)

colnames(stock_ptf_test_input)=c("PEDON_ID", 'soc', 'clay', 'sand', 'hzn_thickness', 'hzn')

stock_ptf_test_input=reshape(stock_ptf_test_input, idvar='PEDON_ID', timevar='hzn', direction='wide')
stock_ptf_test_input=na.omit(stock_ptf_test_input)

stock_ptf_test_input=stock_ptf_test_input[,colnames(stock_ptf_train[,-1])]

pred_stock_ptf=predict(model_stock_ptf, data=stock_ptf_test_input)$predictions

pred_stock_ptf=data.frame(stock_ptf_test_input$PEDON_ID, pred_stock_ptf)
colnames(pred_stock_ptf)=c("PEDON_ID", 'pred_stock')

pred_stock_ptf=merge(pred_stock_ptf, stock_ptf_test[,1:2], by='PEDON_ID')





###############Multipanel Plots##############
#2600 by 3200
grid.arrange(carbon_plot, tn_plot, ioc_all_plot, cec_plot, ec_plot_all, clay_plot_all, sand_plot_all,hzn_plot_all, ncol=2)
#1600 by 2400
grid.arrange(bd_plot, carbon_stock_plot, ncol=1)

#########Combined Performance#################
soc_performance$model='soc'
tn_performance$model='tn'
ioc_performance$model='ioc'
ec_performance$model='ec'
cec_performance$model='cec'
clay_performance$model='clay'
silt_performance$model='silt'
sand_performance$model='sand'
hzn_performance$model='hzn'

performance_table=rbind(soc_performance, tn_performance, ioc_performance, ec_performance, cec_performance, clay_performance, silt_performance, sand_performance, hzn_performance)
write.csv(performance_table, 'all_models_performance.csv')



