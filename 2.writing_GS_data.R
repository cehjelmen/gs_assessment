library(reshape)
#read in data for all things
#data can be collected from
#genomesize.com
#http://www.zbi.ee/fungal-genomesize/
#https://cvalues.science.kew.org/
#and supplemental data from
#Sylvester T, Hjelmen CE, Hanrahan SJ, Lenhart PA, Johnston JS, Blackmon H. Lineage-specific patterns of chromosome evolution are the rule not the exception in Polyneoptera insects. Proceedings B DOI: 10.1098/rspb.2020.1388
#Hjelmen CE, Blackmon H, Holmes VR#, Burrus CG#, Johnston JS. Genome Size evolution differs between Drosophila subgenera with striking differences in male and female genome size in Sophophora. G3
animal.gs.dat<-read.csv("gs_estimates/genome_size_data_Feb_10_24.csv")
plant.gs.dat<-read.csv("gs_estimates/kew_GS_plants_feb_20_2024.csv")
fungus.gs.dat<-read.csv("gs_estimates/fungus_GS_feb_20.csv")
polyneop<-read.csv("gs_estimates/polyneoptera.csv")
dros<-read.csv("gs_estimates/dros.data.csv")

#next is to get gs data
####cleaning animal genome size data####
#remove extra columns
#such as referneces, standard, cell type, common name, referneces
animal.gs.dat[9:15]<-NULL
animal.gs.dat[7]<-NULL
animal.gs.dat$Sub.Phylum<-NULL
colnames(animal.gs.dat)

#removing values which provide ranges
char_to_remove<-"-"
animal.gs.dat$C.value<-ifelse(grepl(char_to_remove, animal.gs.dat$C.value), "", animal.gs.dat$C.value)
#make numeric
animal.gs.dat$C.value<-as.numeric(animal.gs.dat$C.value)
animal.gs.dat <- animal.gs.dat %>%
  mutate(across(c(Phylum, Class, Order, Family), as.factor))
#remove.nas
animal.gs.dat<-na.omit(animal.gs.dat)
str(animal.gs.dat)
#we're to 8398 records
#across 6525 unique species
length(unique(animal.gs.dat$Species))

#convert C.value into bases
basepair<-matrix(,nrow=length(animal.gs.dat$Phylum), ncol=2)
colnames(basepair)<-c("GS_Species", "Bp")

for(i in 1: length(animal.gs.dat$Phylum)){
  basepair[i,1]<-animal.gs.dat$Species[i]
  #c-value is in pg, so multiply by 978 to get Mbp
  #then multiply by 1 million to get bases
  basepair[i,2]<-animal.gs.dat$C.value[i]*978*1000000
}
basepair<-as.data.frame(basepair)
basepair$Bp<-as.numeric(basepair$Bp)

animal.gs.dat<-data.frame(animal.gs.dat, basepair)
#remove extra columns that other data doesn't have
colnames(animal.gs.dat)
animal.gs.dat[1:4]<-NULL
#str(animal.gs.dat.bp)

####polyneoptera data####
#make factors
polyneop <- polyneop %>%
  mutate(across(c(Phylum, Class, Order, Family, Method), as.factor))
#make matrix for names
polynames<-matrix(,nrow=length(polyneop$Phylum), ncol=1)
colnames(polynames)<-"Species"
for(i in 1:length(polyneop$Phylum)){
  polynames[i,1]<-paste(polyneop$Genus[i],polyneop$Species[i])
}
polynames<-as.data.frame(polynames)
polyneop$Species<-polynames$Species
#converte Mbp into bp
basepair<-matrix(,nrow=length(polyneop$Phylum), ncol=2)
colnames(basepair)<-c("GS_Species", "Bp")

for(i in 1: length(polyneop$Phylum)){
  basepair[i,1]<-polyneop$Species[i]
  #c-value is in Mbp so multiply by 1 million to get bases
  basepair[i,2]<-polyneop$C.value[i]*1000000
}
basepair<-as.data.frame(basepair)
basepair$Bp<-as.numeric(basepair$Bp)

poly.gs.dat<-data.frame(polyneop, basepair)
colnames(poly.gs.dat)
#remove extra columns
poly.gs.dat[1:5]<-NULL
poly.gs.dat$Method<-NULL

####Drosophila data#####
#make new object where you'll remove underreplication values
dros.gs<-dros
dros.gs$Female_UR<-NULL
dros.gs$Male_UR<-NULL
#melt into one dataframe
dros.gs<-melt(dros.gs, id=c("Phylum", "Class", "Order", "Family", "Species", "Method"))
#rename gs value
colnames(dros.gs)[8]<-"C.value"
dros.gs$variable<-NULL
dros.gs <- dros.gs %>%
  mutate(across(c(Phylum, Class, Order, Family, Method), as.factor))
#get basepair values
basepair<-matrix(,nrow=length(dros.gs$Phylum), ncol=2)
colnames(basepair)<-c("GS_Species", "Bp")

for(i in 1: length(dros.gs$Phylum)){
  basepair[i,1]<-dros.gs$Species[i]
  #c-value is in Mbp so multiply by 1 million to get bases
  basepair[i,2]<-dros.gs$C.value[i]*1000000
}
basepair<-as.data.frame(basepair)
basepair$Bp<-as.numeric(basepair$Bp)
dros.gs.dat<-data.frame(dros.gs, basepair)
length(unique(dros.gs.dat$Species))
#remove extra columns
colnames(dros.gs.dat)
dros.gs.dat[1:4]<-NULL
dros.gs.dat$Method<-NULL
#genome size data for plants was collected from https://cvalues.science.kew.org/ 
#on feb 20, 2024
#all data was copy and pasted into csv and read in object called "plant.gs.dat"
####Plant GS Data####
colnames(plant.gs.dat)
#remove columns we won't use
plant.gs.dat[8:9]<-NULL
plant.gs.dat[5:6]<-NULL
plant.gs.dat$Subspecies<-NULL

#make combined genus/species name
plant.names<-matrix(,nrow=length(plant.gs.dat$Family), ncol=1)
colnames(plant.names)<-"Species"
for(i in 1:length(plant.gs.dat$Family)){
  plant.names[i,1]<-paste(plant.gs.dat$Genus[i],plant.gs.dat$Species[i])
}
plant.names<-as.data.frame(plant.names)
plant.gs.dat$Species<-plant.names$Species

#str(plant.gs.dat)
plant.gs.dat$Mbp<-as.numeric(plant.gs.dat$Mbp)

#plant.gs.dat$Species[plant.gs.dat$Mbp=="NA"]
#converte Mbp into bp
basepair<-matrix(,nrow=length(plant.gs.dat$Family), ncol=2)
colnames(basepair)<-c("GS_Species", "Bp")
colnames(plant.gs.dat)[4]<-"C.value"
for(i in 1: length(plant.gs.dat$Family)){
  basepair[i,1]<-plant.gs.dat$Species[i]
  #c-value is in Mbp so multiply by 1 million to get bases
  basepair[i,2]<-plant.gs.dat$C.value[i]*1000000
}
basepair<-as.data.frame(basepair)
basepair$Bp<-as.numeric(basepair$Bp)

plant.gs.dat<-data.frame(plant.gs.dat, basepair)
plant.gs.dat$Genus<-NULL
#remove extra
colnames(plant.gs.dat)
plant.gs.dat$Family<-NULL
length(unique(plant.gs.dat$Species))


####Fungus GS data####
#fungus data was copy and pasted into csv on Feb 20, 2024
#from http://www.zbi.ee/fungal-genomesize/
colnames(fungus.gs.dat)
fungus.gs.dat[8:10]<-NULL
fungus.gs.dat[5:6]<-NULL
colnames(fungus.gs.dat)<-c("Phylum", "Order", "Genus", "Species", "C.value")

#make combined genus/species name
fungus.names<-matrix(,nrow=length(fungus.gs.dat$Phylum), ncol=1)
colnames(fungus.names)<-"Species"
for(i in 1:length(fungus.gs.dat$Phylum)){
  fungus.names[i,1]<-paste(fungus.gs.dat$Genus[i],fungus.gs.dat$Species[i])
}
fungus.names<-as.data.frame(fungus.names)
fungus.gs.dat$Species<-fungus.names$Species
#str(fungus.gs.dat)
#converte Mbp into bp
basepair<-matrix(,nrow=length(fungus.gs.dat$Phylum), ncol=2)
colnames(basepair)<-c("GS_Species", "Bp")

for(i in 1: length(fungus.gs.dat$Phylum)){
  basepair[i,1]<-fungus.gs.dat$Species[i]
  #c-value is in Mbp so multiply by 1 million to get bases
  basepair[i,2]<-fungus.gs.dat$C.value[i]*1000000
}
basepair<-as.data.frame(basepair)
basepair$Bp<-as.numeric(basepair$Bp)

fungus.gs.dat<-data.frame(fungus.gs.dat, basepair)
str(fungus.gs.dat)
#remove extra columns
colnames(fungus.gs.dat)
fungus.gs.dat[1:3]<-NULL

#leads to 8762 animal records for 6619 species
gs.an.dros.combined<-merge(animal.gs.dat, dros.gs.dat,all=TRUE)
gs.an.dros.poly<-merge(gs.an.dros.combined, poly.gs.dat, all=TRUE)
#length(unique(gs.an.dros.poly$Species))

#plant data is 12,273 records for 11634 species
#length(unique(plant.gs.dat$Species))
#fungus data is 2412 records for 1324 species
#length(unique(fungus.gs.dat$Species))
gs.animal.fungus<-merge(gs.an.dros.poly, fungus.gs.dat, all=TRUE)

#including plant data
gs.all.data<-merge(gs.animal.fungus, plant.gs.dat, all=TRUE)

#totals to 23,447 records across 19577 species
length(unique(gs.all.data$Species))

write.csv(gs.all.data, "all.gs.data.bp.csv")
