library(reshape)
library(dplyr)
#Dros Data analysis
#read in file
genome.seq.info.df<-read.csv("all_genome_seq_info_Feb_27.csv")
dros.genome<-subset(genome.seq.info.df, subset=(Family=="Drosophilidae"))
#633 genomes from Drosophilidae
#read in drosophila genome and ur data
dros<-read.csv("gs_estimates/dros.data.csv")
#make new object where you'll remove underreplication values
dros.gs<-dros
dros.gs$Female_UR<-NULL
dros.gs$Male_UR<-NULL
#melt into one dataframe
dros.data.test<-melt(dros.gs, id=c("Phylum", "Class", "Order", "Family", "Species", "Method"))
#rename gs value
colnames(dros.data.test)[8]<-"C.value"
dros.data.test$variable<-NULL
dros.data.test <- dros.data.test %>%
  mutate(across(c(Phylum, Class, Order, Family, Method), as.factor))

basepair.dros<-matrix(,nrow=length(dros.data.test$Phylum), ncol=2)
colnames(basepair.dros)<-c("GS_Species", "Bp")

for(i in 1: length(dros.data.test$Phylum)){
  basepair.dros[i,1]<-dros.data.test$Species[i]
  #c-value is in Mbp so multiply by 1 million to get bases
  basepair.dros[i,2]<-dros.data.test$C.value[i]*1000000
}
basepair.dros.df<-as.data.frame(basepair.dros)
basepair.dros.df$Bp<-as.numeric(basepair.dros.df$Bp)
dros.gs.dat.bp<-data.frame(dros.data.test, basepair.dros.df)

#make new object where you'll remove gs values
dros.ur<-dros
dros.ur$MbDNA_Female<-NULL
dros.ur$MbDNA_Male<-NULL
#melt into one dataframe
dros.data.ur<-melt(dros.ur, id=c("Phylum", "Class", "Order", "Family", "Species", "Method"))
#rename gs value
colnames(dros.data.ur)[8]<-"UR"
dros.data.test$variable<-NULL
dros.data.test <- dros.data.test %>%
  mutate(across(c(Phylum, Class, Order, Family, Method), as.factor))

dros.gs.dat.ur<-data.frame(dros.gs.dat.bp, dros.data.ur$UR)
colnames(dros.gs.dat.ur)[10]<-"UR"

write.csv(dros.gs.dat.ur, "dros_UR_GS_data.csv")

