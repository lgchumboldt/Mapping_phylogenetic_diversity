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
    Birds_dataframe$residuos_Pd_lineal<-rep(NA,length(Birds_dataframe[,1]))
    Birds_dataframe$residuos_Pd_lineal[as.numeric(names(linear_residuals))]<-linear_residuals[names(linear_residuals)]
    #LOESS
    Birds_dataframe$residuos_Pd_loess<-rep(NA,length(Birds_dataframe[,1]))
    Birds_dataframe$residuos_Pd_loess[as.numeric(names(loess_residuals))]<-loess_residuals[names(loess_residuals)]

    #Generate rasters with residual values
    Pd_aves<-raster("C:\\Users\\GIC 40\\Dropbox\\PD_aves_IUCN_083.asc")
    residuos_lineal_raster<-Pd_aves
    values(residuos_lineal_raster)<-NA
    values(residuos_lineal_raster)<-marco_aves$residuos_Pd_lineal

    Pd_aves<-raster("C:\\Users\\GIC 40\\Dropbox\\PD_aves_IUCN_083.asc")
    residuos_loess_raster<-Pd_aves
    values(residuos_loess_raster)<-NA
    values(residuos_loess_raster)<-marco_aves$residuos_Pd_loess
    writeRaster(residuos_lineal_raster,"C:\\Users\\GIC 40\\Dropbox\\residuals_lineal_PD.asc","ascii")
    writeRaster(residuos_loess_raster,"C:\\Users\\GIC 40\\Dropbox\\residuals_loess_PD.asc","ascii")




    #Graphs and maps
        #use this part of the code if you wish to have in the same image regression models, residuals and residual maps
    par(mfrow=c(3,2),mar=c(4,4,3,3))
    #Regressions
  

    plot(Birds_dataframe$dataframe_Birds_SR.TD_aves_083,Birds_dataframe$dataframe_Birds_PD.PD_aves_IUCN_083,xlab="N�mero de especies",ylab="Diversidad filogen�tica (PD)", main="Regresi�n lineal")
    abline(modelo_lineal,col="red")
     plot(Birds_dataframe$dataframe_Birds_SR.TD_aves_083,Birds_dataframe$dataframe_Birds_PD.PD_aves_IUCN_083,xlab="N�mero de especies",ylab="Diversidad filogen�tica (PD)", main="Regresi�n loess")
     Td_vals<-seq(0,600,1)
     loess_vals<-predict(modelo_loess,Td_vals)
     lines(loess_vals ~ Td_vals,col="red")

    #Residuals graph

    plot(Birds_dataframe$dataframe_Birds_SR.TD_aves_083,Birds_dataframe$residuos_Pd_lineal,xlab="N�mero de especies",ylab="Residuos diversidad filogen�tica (PD)")
    abline(0,0)
    plot(Birds_dataframe$dataframe_Birds_SR.TD_aves_083,Birds_dataframe$residuos_Pd_loess,xlab="N�mero de especies",ylab="Residuos diversidad filogen�tica (PD)")
    abline(0,0)

    #Maps
    plot(residuos_lineal_raster,box=F,axes=F)
    plot(residuos_loess_raster,box=F,axes=F)
    plot(colombia)
    extent(Pd_aves)
