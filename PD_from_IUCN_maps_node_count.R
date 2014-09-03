######Phylogenetic diversity as node count from distribution maps iin ESRI shapefile format and phylogenies in nexus format######
#############Script written by Andrea Paz in August 2014#################

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
  r<-raster(ncol=1462,nrow=624)
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

########Function to obtain node track###########
##http://rstudio-pubs-static.s3.amazonaws.com/15766_a6d320bf9a754d4ea43c3d8cdf4b6b04.html

####################### FastXtreePhylo #######################
#
# Based on an algorithm by Jens Schumacher developed in 2003. The
# orginal version can be found at:
#      http://www.thetrophiclink.org/resources/calculating-fd
# - last accessed 25 March 2014
#
# Adapted to work with objects of class "phylo" defined in the
# package ape. From the old FastXtree algorithm we borrow the principle
# that when we descend from tip (==terminal node) level towards the root,
# row entries in the currently focussed edge (ie col) of H1 are equal
# to the sum of all H1 columns representing edges above the focal edge
# in the tree. Using the rowSums() function vectorises this step making
# it very fast.
#
# PDW 10 March 2014

FastXtreePhylo <- function(thisTree)
{
  # In edge matrix in a phylo object, OTUs are indexed 1:Ntip, the root node has
  # index Ntip + 1, and internal nodes are indexed from root.node to 
  # total.nodes. This is clever because it lets us very efficiently
  # cross-reference between edge indicies and the left and right terminal node
  # indices of each edge. So, compute the key values:
  root.node <- Ntip(thisTree) + 1
  total.nodes <- max(thisTree$edge)
  
  # Matrix H1 of Schumacher's Xtree method, which is exactly the same as
  # David Nipperess's phylomatrix:
  H1 <- matrix(0,Ntip(thisTree),Nedge(thisTree),dimnames=list(thisTree$tip.label,NULL))
  
  # A short-cut: we can instantly populate H1 with terminal edges:
  rc.ind <- cbind("row"=1:Ntip(thisTree),"col"=which(thisTree$edge[,2] < root.node))
  H1[rc.ind] <- 1
  
  # Make a vector of internal node indices in descending order so we can traverse
  # tree from first nodes below tip nodes down to the root:
  internal.nodes <- seq(total.nodes,root.node,-1)
  
  # Now, visit each internal node accumulating subtended child edge records:
  for (thisNode in internal.nodes)
  {
    # Find the edge which has thisNode on its left:
    next.edge <- which(thisTree$edge[,2]==thisNode)
    
    # Which edges are children of the node to the right:
    child.edges <- which(thisTree$edge[,1]==thisNode)
    
    # Do the magic rowSums() trick to accumulate edges subtended by the current edge:
    H1[,next.edge] <- rowSums(H1[,child.edges])
  }
  
  # All done...package results as a named list:
  return(list("h2.prime"=thisTree$edge.length,"H1"=H1))  
}

####End of function#####

#######Funcion: count nodes involved in a community######

compute_node_based_length <- function( community,matrix_nodes) {
  
  community_matrix <- matrix_nodes[community,]
  if (length(community) == 1) {
    result <- community_matrix
  }
  else if (dim(community_matrix)[1] > 3 ) {
    a <- community_matrix[1,]
    b <- community_matrix[2,]
    result <- bitwOr(a,b)
    community_matrix <- community_matrix[3:dim(community_matrix)[1],]
    community_matrix <- rbind(community_matrix, result)
    compute_node_based_length(community_matrix,matrix_nodes)
  }
  else {
    a <- community_matrix[1,]
    b <- community_matrix[2,]
    result <- bitwOr(a,b)
  }
  node_length <- sum(result)
  return(node_length)
}

####End of function####

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
  print("Select shapefile mask to use for anlyses")
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
  write.table(table1,"species_list.txt",quote=F,row.names=F) #seria bueno que agregara el nombre de la mascara para saber si es colombia o paramo

  #Read Shape of the geographic Mask 
  selected_mask<-readShapePoly(working_mask)

#########################
####START ANALYSES#######
#########################
   
#Create empty raster of the desired area and resolution to assign pixel numbers
r<-raster(ncol=1462,nrow=624)
res(r)<-resolution #resolution
r<-crop(r,extent(selected_mask))
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

#Combine distributions with pixel number raster
complete_list<-c(grid_list,layers)

#Stack all distribution files and pixel numbers for PD computation
Stack<-stack(complete_list)

#Convert to data frame for storage and posterior analysis

community_data<-as.data.frame(Stack)

#Store 
write.table(community_data, file = "communities_iucn_0_083.txt", append = FALSE,row.names=F,quote=F,sep="\t") #tratar de agregar nombre de mascara

#Generate the list of species present in the mask area and write to file to obtain corresponding phylogeny

distribution_files1<-list.files(path=maps_folder, pattern= "*.asc$")
species_names1<-sub(" .asc","",distribution_files1)
table2<-as.data.frame(species_names1)
colnames(table2)<-"Grid"
write.table(table2,"final_species_list.txt",quote=F,row.names=F) #seria bueno incluir en el nombre del archivo el nombre de la mascara


#Phylogenetic diversity computation

#Read phylogeny file in nexus format
  ##To read a phylogeny in newick format use read.tree instead of read.nexus
      #In working directory always a file with the same name: phylogeny.nex

  user_phylogeny=read.nexus("phylogeny.nex")
  anfibios<-collapse.singles(user_phylogeny)

#In the community data frame NA must be eliminated
community_data=na.omit(community_data)
head(community_data)


#community_data
lista<-match.phylo.comm(anfibios,community_data)
community_data1<-as.data.frame(lista[[2]])
community_data1$Grid<-community_data$Grid
taxa_elim<-lista$taxa.not.tree
taxa_elim<-gsub("[.]"," ",taxa_elim)
layers[taxa_elim]<-NULL
#Generate a raster of taxonomic diversity (TD) for comparison
all_species<-stack(layers)
TD1<-calc(all_species,sum)
writeRaster(TD1,"SR_anfibios_solo_filo_0_083.asc") 
writeRaster(TD1,species_richness) 

##fix phylogeny 
anfibios2<-lista[[1]]


#III-Phylogenetic diversity computation 
####Obtain node tracks

matriz_nodos<-FastXtreePhylo(anfibios2)
matriz_nodos<-matriz_nodos[[2]]

###Obtain a list with species present in each community####
comm<-list()
for (i in 1:dim(community_data1)[1]){
  comunidad<-community_data1[i,]
  especies<-colnames(comunidad[which(comunidad==1)])
  pixel<-comunidad$Grid
  if (length(especies)==0){
    #print(i)
  }
  else if(length(especies==1)){
    comm[[i]]<-especies
    names(comm)[i]<-pixel
  }
}
comm[sapply(comm,is.null)]<-NULL

####Compute node based pd value for every community 
valores_nodos<-lapply(comm,compute_node_based_length,matrix_nodes=matriz_nodos)


##Assign node based pd values to each community in the matrix
community_data$pd_node<-rep(NA,dim(community_data)[1]) 
for(i in 1:dim(community_data)[1]){
  posicion<-which(as.numeric(names(valores_nodos))==community_data$Grid[i])
  if(length(posicion==1)){
    community_data$pd_node[i]<-valores_nodos[[posicion]]
  }
  else {
    #print(posicion)
  }
}
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
pd_ras[community_data$Grid]<-community_data$pd_node

#Write the new matrix to a file to avoid rerunning the script for potential further analyses
write.table(community_data, file = "communities_and_pd.txt", append = FALSE,row.names=F,quote=F,sep="\t")
#3- Save raster to file 

#writeRaster(pd_ras,PD_map)
writeRaster(pd_ras,PD_map)

#4-Optional plotting map in R (write to pdf?)
plot(pd_ras)

}
