library(maptools)
library(raster)
library(rgdal)

#Cargar paramos
ruta_paramos<-"C:\\Users\\GIC 40\\Documents\\Andrea\\Informe_Estado_Tendencias"
paramos<-readOGR(dsn=ruta_paramos,layer="paramos_colombia_project1")
#Cargar Colombia
colombia_shape<-"C:\\Users\\GIC 66\\Documents\\Andrea\\SIG\\Colombia\\adm"
colombia<-readOGR(dsn=colombia_shape,layer="COL_adm0")
#Cargar Raster PD, TD 
    #AVES LBAB
Pd_aves<-raster("C:\\Users\\GIC 66\\Dropbox\\PD_aves_LBAB_agg8.asc")
Td_aves<-raster("C:\\Users\\GIC 66\\Dropbox\\TD_aves_LBAB_agg8.asc")
#AVES IUCN
Pd_aves<-raster("PD_aves_IUCN_083.asc")
Td_aves<-raster("C:\\Users\\GIC 66\\Dropbox\\TD_aves_IUCN_083.asc")


#Crear un marco de datos con valores de Pd y Td
  #AVES
marco_Pd_aves<-as.data.frame(Pd_aves)
marco_Td_aves<-as.data.frame(Td_aves)
marco_aves<-data.frame(marco_Pd_aves$PD_aves_IUCN_083,marco_Td_aves$TD_aves_IUCN_083)


#Generar modelos de regresi�n 
  #LINEAL
  modelo_lineal<-lm(marco_aves)
  #LOESS
  modelo_loess<-loess(marco_aves)

#Calcular y guardar los residuales de las regresiones
  #LINEAL
  residuos_lineal<-resid(modelo_lineal)
  #LOESS
  residuos_loess<-resid(modelo_loess)
#Agregar los resiudales al marco de datos
marco_aves$residuos_Pd_lineal<-rep(NA,length(marco_aves[,1]))
marco_aves$residuos_Pd_lineal[as.numeric(names(residuos_lineal))]<-residuos_lineal[names(residuos_lineal)]


marco_aves$residuos_Pd_loess<-rep(NA,length(marco_aves[,1]))
marco_aves$residuos_Pd_loess[as.numeric(names(residuos_loess))]<-residuos_loess[names(residuos_loess)]

#Generar rasters con los datos de residuos
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




#Graficas y mapas
  #Regresiones
  par(mfrow=c(3,2),mar=c(4,4,3,3))

  plot(marco_aves$marco_Td_aves.TD_aves_083,marco_aves$marco_Pd_aves.PD_aves_IUCN_083,xlab="N�mero de especies",ylab="Diversidad filogen�tica (PD)", main="Regresi�n lineal")
  abline(modelo_lineal,col="red")
  plot(marco_aves$marco_Td_aves.TD_aves_083,marco_aves$marco_Pd_aves.PD_aves_IUCN_083,xlab="N�mero de especies",ylab="Diversidad filogen�tica (PD)", main="Regresi�n loess")
  Td_vals<-seq(0,600,1)
  loess_vals<-predict(modelo_loess,Td_vals)
  lines(loess_vals ~ Td_vals,col="red")

#Residuos

plot(marco_aves$marco_Td_aves.TD_aves_083,marco_aves$residuos_Pd_lineal,xlab="N�mero de especies",ylab="Residuos diversidad filogen�tica (PD)")
abline(0,0)
plot(marco_aves$marco_Td_aves.TD_aves_083,marco_aves$residuos_Pd_loess,xlab="N�mero de especies",ylab="Residuos diversidad filogen�tica (PD)")
abline(0,0)

#Mapas
plot(residuos_lineal_raster,box=F,axes=F)
plot(residuos_loess_raster,box=F,axes=F)
plot(colombia)
extent(Pd_aves)
