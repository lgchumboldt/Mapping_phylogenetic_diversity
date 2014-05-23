library(maptools)
library(raster)
library(rgdal)

#Load study area masks
    #Paramos
ruta_paramos<-"C:\\Users\\GIC 40\\Documents\\Andrea\\Informe_Estado_Tendencias"
paramos<-readOGR(dsn=ruta_paramos,layer="paramos_colombia_project1")
    #Colombia
colombia_shape<-"C:\\Users\\GIC 66\\Documents\\Andrea\\SIG\\Colombia\\adm"
colombia<-readOGR(dsn=colombia_shape,layer="COL_adm0")

    #Load pre-existing phylogenetic diversity (PD) and species richness (SR) maps
   
    #Rasters created using IUCN distribution maps
Birds_PD<-raster("PD_aves_IUCN_083.asc")
Birds_SR<-raster("C:\\Users\\GIC 66\\Dropbox\\TD_aves_IUCN_083.asc")


    #Create a dataframe including both PD and SR
  
dataframe_Birds_PD<-as.data.frame(Birds_PD)
dataframe_Birds_SR<-as.data.frame(Birds_SR)
Birds_dataframe<-data.frame(dataframe_Birds_PD$PD_aves_IUCN_083,dataframe_Birds_SR$TD_aves_IUCN_083)


    #Generate regression models
    #LINEAR
    linear_model<-lm(Birds_dataframe)
    #LOESS
    loess_model<-loess(Birds_dataframe)

    #Compute and save residuals from both regression models
    #LINEAR
    linear_residuals<-resid(linear_model)
    #LOESS
    loess_residuals<-resid(loess_model)
  
    #Add residuals to dataframe
    #LINEAR
    Birds_dataframe$PD_linear_residuals<-rep(NA,length(Birds_dataframe[,1]))
    Birds_dataframe$PD_linear_residuals[as.numeric(names(linear_residuals))]<-linear_residuals[names(linear_residuals)]
    #LOESS
    Birds_dataframe$PD_loess_residuals<-rep(NA,length(Birds_dataframe[,1]))
    Birds_dataframe$PD_loess_residuals[as.numeric(names(loess_residuals))]<-loess_residuals[names(loess_residuals)]

    #Generate rasters with residual values
    PD_birds<-raster("C:\\Users\\GIC 40\\Dropbox\\PD_aves_IUCN_083.asc")
    
    #LINEAR
    linear_residuals_rasterresiduos_lineal_raster<-PD_birds
    values(linear_residuals_raster)<-NA
    values(linear_residuals_raster)<- Birds_dataframe$PD_linear_residuals
     writeRaster(linear_residuals_raster,"C:\\Users\\GIC 40\\Dropbox\\residuals_lineal_PD.asc","ascii")
    #LOESS
    loess_residuals_raster<-PD_birds
    values(loess_residuals_raster)<-NA
    values(loess_residuals_raster)<- Birds_dataframe$PD_loess_residuals
    writeRaster(loess_residuals_raster,"C:\\Users\\GIC 40\\Dropbox\\residuals_loess_PD.asc","ascii")




    #Graphs and maps
        #use this part of the code if you wish to have in the same image regression models, residuals and residual maps
    par(mfrow=c(3,2),mar=c(4,4,3,3))
    #Regressions
  

    plot(Birds_dataframe$dataframe_Birds_SR.TD_aves_083,Birds_dataframe$dataframe_Birds_PD.PD_aves_IUCN_083,xlab="Species richness",ylab="Phylogenetic diversity (PD)", main="Linear regression")
    abline(modelo_lineal,col="red")
     plot(Birds_dataframe$dataframe_Birds_SR.TD_aves_083,Birds_dataframe$dataframe_Birds_PD.PD_aves_IUCN_083,xlab="Species richness",ylab="Phylogenetic diversity (PD)", main="Loess regression")
     SR_vals<-seq(0,600,1)
     loess_vals<-predict(loess_model,SR_vals)
     lines(loess_vals ~ SR_vals,col="red")

    #Residuals graph

    plot(Birds_dataframe$dataframe_Birds_SR.TD_aves_083,Birds_dataframe$PD_linear_residuals,xlab="Species richness",ylab="Phylogenetic diversity (PD) residuals")
    abline(0,0)
    plot(Birds_dataframe$dataframe_Birds_SR.TD_aves_083,Birds_dataframe$residuos_Pd_loess,xlab="Species richness",ylab="Phylogenetic diversity (PD) residuals")
    abline(0,0)

    #Maps
    plot(linear_residuals_raster,box=F,axes=F)
    plot(loess_residuals_raster,box=F,axes=F)
    plot(colombia)
    extent(PD_birds)
