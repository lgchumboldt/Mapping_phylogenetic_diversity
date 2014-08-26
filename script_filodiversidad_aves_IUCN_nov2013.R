######Phylogenetic diversity from distribution maps iin ESRI shapefile format and phylogenies in nexus format######
####In case you have distribution maps in raster format use the other script####
#############Written by Andrea Paz and David Urbina#################

###For using this script you need a folder containing distrbution files in ESRI shapefile format. If you have distribution maps in raster formats use the SDM script
#Clean workspace
rm(list=ls())

#Load required packages
library(maptools)
library(raster)
library(rgdal)
library(picante)

########Function: rasterize_species############
####This function converts IUCN distribution files (.shp) to raster files while cutting them to the selected area######
####It elminates species for which the distribution does not include the area of the selected mask (values null and 0)#############

rasterize_species= function (x,mask=selected_mask) {
  r<-raster(ncol=1462,nrow=624) #This is based on Colombia extent if you have a bigger area you should change it
  res(r)<-resolution 
  r<-crop(r, extent(selected_mask))
  values(r)<-0
  map<-readOGR(dsn=maps_folder,layer=x)
  r<-rasterize(map,r,1,update=T,background=0)
  r<-mask(r,selected_mask)
  valor<-unique(getValues(r))
  
  if(length(valor)==1&&is.na(valor)==TRUE){
    
    
  }
  else if (length(valor)==2&&valor[2]==0){
    
  }
  else {
    writeRaster(r,paste(x,".asc"))
    return (raster(paste(x,".asc")))
  }
}

#######END OF FUNCTION######

########################################USER DEFINED VARIABLES#########################################
######Prompt user to determine location of files and selected geographic areas, resolution etc.########
#######################################################################################################


  {print("Remember phylogeny file must be stored in working directory and named phylogeny.nex")

  # set working directory
  print("Choose working directory:")
  working_directory<-choose.dir()
  setwd(working_directory)

  #Select folder containing distribution maps
  print("Select folder containing distribution maps")
  maps_folder<-choose.dir()
  
  #File containing shape of the area to be analyzed (can be a country, province, ecosystem etc.). It will be treated as a mask
  working_mask<-file.choose()

  #Determine resolution
  resolution1<-readline("PD raster resolution: ")
   resolution<-as.numeric(resolution1)
  #Determine names for outputs
  PD_map<-readline("Name of phylogenetic diversity raster: Do not use spaces and must include .asc extention (Ex: PD_country.asc): ")
  species_richness<-readline("Name of the species richness map: Do not use spaces and must include .asc extention  (Ex: Species_richness_country.asc): ")
  
  #Read distribution files and obtain species names for phylogeny

  distribution_files<-list.files(path=maps_folder, pattern= "*.shp$")
  species_names<-sub(".shp","",distribution_files)
  table1<-as.data.frame(species_names)
  colnames(table1)<-"Grid"
  write.table(table1,"species_list.txt",quote=F,row.names=F) 

  #Read Shape of the geographic Mask 
  selected_mask<-readShapePoly(working_mask)

#########################
####START ANALYSES#######
#########################
   
#Create empty raster of the desired area and resolution to assign pixel numbers
r<-raster(ncol=1462,nrow=624)
res(r)<-resolution #resolution is user-defined
r<-crop(r,extent(selected_mask))
#Generate a raster that has the pixel number as the pixel value
grid=r
names(grid)="grid"
grid[1:ncell(grid)]<-1:ncell(grid) 
grid_list<-list(grid)
names(grid_list)<-"Grid"

#Rasterize distribution files and keep only those of interest
layers<-lapply(species_names,rasterize_species)

#assign species names to layers with species distributions
names(layers)<-as.vector(table1$Grid)
layers[sapply(layers,is.null)]<-NULL

#Combine species distribution rasters with pixel number raster
complete_list<-c(grid_list,layers)

#Stack all distribution files and pixel numbers for PD computation
Stack<-stack(complete_list)

#Convert to data frame for storage and posterior analysis

community_data<-as.data.frame(Stack)

#Store 
write.table(community_data, file = "communities_iucn_0_083.txt", append = FALSE,row.names=F,quote=F,sep="\t") 

#Generate the list of species present in the mask area and write to file to obtain corresponding phylogeny

distribution_files1<-list.files(path=maps_folder, pattern= "*.asc$")
species_names1<-sub(" .asc","",distribution_files1)
table2<-as.data.frame(species_names1)
colnames(table2)<-"Grid"
write.table(table2,"final_species_list.txt",quote=F,row.names=F) 

#Generate a raster of species richness (SR) for comparison
all_species<-stack(layers)
SR<-calc(all_species,sum)
writeRaster(SR,species_richness)


#Phylogenetic diversity computation

#Read phylogeny file in nexus format
  ##To read a phylogeny in newick format use read.tree instead of read.nexus
      #In working directory always a file with the same name: phylogeny.nex

  user_phylogeny=read.nexus("phylogeny.nex")


#In the community data frame NA must be eliminated
community_data=na.omit(community_data)
head(community_data)

#III-Phylogenetic diversity computation 
  #computes only Faith´s PD others may be added

pd.result <-pd(community_data[,2:ncol(community_data)],user_phylogeny,include.root = TRUE) 


#Add the pixel PD value to data frame
community_data$pd<-pd.result[[1]]

#Write the new matrix to a file to avoid rerunning the script for potential further analyses
write.table(community_data, file = "communities_and_pd.txt", append = FALSE,row.names=F,quote=F,sep="\t")

#IV-Generate a raster containing PD information per pixel

#1-First generate an empty raster using a base model for resolution and area

r<-raster(ncol=1462,nrow=624)
res(r)<-resolution
r<-crop(r,extent(selected_mask))
values(r)<-0
map<-readOGR(dsn=maps_folder,layer=species_names[[1]])
r<-rasterize(map,r,map$PRESENCE,update=T)
r<-mask(r,selected_mask)
pd_ras<-r
values(pd_ras)<-NA #Eliminate every value they will be replaced by PD values further down


#2- Assign PD values to raster
pd_ras[community_data$Grid]<-community_data$pd

#3- Save raster to file 
  
writeRaster(pd_ras,PD_map)

#4-Optional plotting map in R (write to pdf?)

plot(pd_ras)
}
