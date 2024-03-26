library(ggplot2)
library(dplyr)
library(ggpubr)
library(forcats)
library(viridis)
#read in drosophila data
dros.ur<-read.csv("dros_UR_GS_data.csv")
genome.seq.info.df<-read.csv("all_GS_assembly.cleaned.csv")

#get out just drosophila info
dros.genome<-subset(genome.seq.info.df, subset=(genome.seq.info.df$Family=="Drosophilidae"))

#make table to have all relevant info in it for UR analysis
dattest<-matrix(,nrow=length(dros.genome$Species), ncol=5)
colnames(dattest)<-c("Species", "GS", "AssblySz", "Seq_Method", "UR")
for(i in 1:length(dros.genome$Species)){
  dattest[i,1]<-dros.genome$Species[i]
  dattest[i,3]<-dros.genome$Assembly_Length_Total[i]
  dattest[i,2]<-mean(dros.ur$Bp[dros.ur$Species==dros.genome$Species[i]])
  dattest[i,4]<-dros.genome$Sequencing_method_1[i]
  dattest[i,5]<-mean(dros.ur$UR[dros.ur$Species==dros.genome$Species[i]])
}

#prepping dataframe
dattest<-as.data.frame(dattest)
dattest$GS<-as.numeric(dattest$GS)
dattest$UR<-as.numeric(dattest$UR)
dattest[sapply(dattest, is.nan)] <- NA
dattest$AssblySz<-as.numeric(dattest$AssblySz)

#figuring out Difference in GS and assembly
Diff<-(dattest$GS-dattest$AssblySz)
#proportional differnece
DiffProp<-Diff/dattest$GS
#combine info
dattest<-data.frame(dattest, Diff, DiffProp)
#combine all data
alldat<-data.frame(dros.genome, dattest)

#get rid of NAs
dat.clean<-subset(alldat, !is.na(UR))
#write.csv(dat.clean, "dros_GS_assembly_ur.cleaned.csv")
dat.clean<-read.csv("dros_GS_assembly_ur.cleaned.csv")

alldat<-dat.clean

alldat <- alldat %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Organism, Assembly_Status, Sequencing_method_1, Sequencing_method_2), as.factor))

sequencers<-c("NA", "ILLUMINA", "OXFORD_NANOPORE", "PACBIO_SMRT")
levels(alldat$Sequencing_method_1)<-c("NA", "Illumina", "Oxford Nanopore", "PacBio SMRT")

#getting heterochromatin percentage
ur.neg<-(1-alldat$UR)
urmod<-lm(alldat$Diff~ur.neg)
#urmod2<-lm(alldat$DiffProp~ur.neg)
#urmod.out<-summary(urmod)
summary(urmod)

urmod.out$adj.r.squared
urmod.out$coefficients

seq.colors<-c("#DF7C18", "#2D87A9", "#AD004F", "grey")
reg.plot<-ggplot(data=alldat, aes(x=(1-UR), y=Diff))+
  geom_point(aes(fill=Sequencing_method_1, shape=Sequencing_method_1),size=5, alpha=0.9, stroke=1)+
  scale_shape_manual(values=c(25, 21, 22, 24))+
  scale_fill_manual(values=c(seq.colors[4],
                             seq.colors[1],
                             seq.colors[2],
                             seq.colors[3]))+
  theme_minimal()+
  #one outlier way up there
  ylim(-1e+08, 2e+08)+
  ylab("GS - Assembly Size (bp)")+
  xlab("Proportion late-replicating heterochromatin")+
  stat_smooth(method=lm,color="black")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold", color="black"),
        axis.text=element_text(size=10, color="black"))+
  annotate("text", x=0.07, y=150000000,
            label=expression("Adj."~R^2~"= 0.333"), size=4, parse=TRUE)+
  annotate("text",x=0.07, y=140000000, label=paste("p < 0.0001"), size=4)
#differnce plot
diff.plot<-ggplot(data=alldat, aes(x=Species, y=Diff, fill=DiffProp>0))+
  geom_bar(stat="identity", width=0.8,position=position_dodge2(), linewidth=0.2)+
  coord_flip()+
  ylim(-1e+08,2e+08)+
  theme_bw()+
  ylab("GS - Assembly size (bp)")+
  scale_fill_viridis_d(end=0.8, direction=-1)+
  geom_hline(aes(yintercept=30000000),
             color="red", linewidth=0.9,
             linetype="dashed")+
  geom_hline(aes(yintercept=0), color="black", linewidth=1.1)+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size=12, face="bold"),
        plot.title = element_text(size=10),
        axis.text.y=element_text(size=5,colour="black", face="italic"),
        axis.text.x=element_text(size=10),
        legend.position = "none",
        panel.grid.major = element_blank())

#final plot for paper
#figure4
ggarrange(diff.plot,theme_void(),reg.plot,
          nrow=1,
          widths=c(0.4,0.05,1),
          labels=c("a", "", "b"),
          font.label = list(size=20, face="bold"))

