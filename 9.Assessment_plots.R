####using chromosome information in the genome assessment####
library(ggplot2)
library(forcats)
library(dplyr)
library(ggpubr)
library(viridis)
library(reshape)
library(EnvStats)
library(cowplot)
library(ggExtra)
####Read in data####
dat.clean<-read.csv("all_GS_assembly.cleaned.csv")
alldat<-dat.clean
chrom<-read.csv("chromosome.info.csv")
str(alldat)
alldat <- alldat %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Organism, Assembly_Status, Sequencing_method_1, Sequencing_method_2), as.factor))

####Get haploid numbers to match####
dattest<-matrix(,nrow=length(dat.clean$Species), ncol=2)
colnames(dattest)<-c("Species", "HaploidNum")
for(i in 1:length(dat.clean$Species)){
  dattest[i,1]<-dat.clean$Species[i]
  dattest[i,2]<-mean(chrom$HaploidNum[chrom$Species==dat.clean$Species[i]])
}
dattest2<-as.data.frame(dattest)
dattest2$HaploidNum<-as.numeric(dattest2$HaploidNum)
dattest2[sapply(dattest2, is.nan)] <- NA
str(dattest2)
alldat<-data.frame(alldat, dattest2)
length(na.omit(alldat$HaploidNum))
#2761 with chormosome number for metazoa

####Subset by sequencing method####
seq.info<-subset(alldat, subset=(Sequencing_method_1=="ILLUMINA" | Sequencing_method_1=="PACBIO_SMRT" | Sequencing_method_1=="OXFORD_NANOPORE"))
seq.info$Sequencing_method_1<-as.character(seq.info$Sequencing_method_1)
seq.info$Sequencing_method_1<-as.factor(seq.info$Sequencing_method_1)
sequencers<-c("ILLUMINA", "OXFORD_NANOPORE", "PACBIO_SMRT")
levels(seq.info$Sequencing_method_1)<-c("Illumina", "Oxford Nanopore", "PacBio SMRT")
names(sequencers)<-c("Illumina", "Oxford Nanopore", "PacBio SMRT")
viridis_colors <- viridis_pal(end=0.8, option="D")(3)
color_mapping<-setNames(viridis_colors, c("Smaller", "Within", "Larger"))
seq.colors<-c("#DF7C18", "#2D87A9", "#AD004F")

####new assessment info####
#make new columns with assessment info
seq.info$ContigN50GSHap<-(seq.info$ContigN50/(seq.info$GS/seq.info$HaploidNum))
seq.info$ScaffoldN50GSHap<-(seq.info$ScaffoldN50/(seq.info$GS/seq.info$HaploidNum))
seq.info$BuscoGSHap<-(seq.info$BUSCO_Complete/(seq.info$GS/seq.info$HaploidNum))

####regression analysis for diffprop####
#contig n50 has significant negative relationship with p=0.0394
contig.lm.gs<-lm(abs(DiffProp)~ContigN50, data=seq.info)
summary(contig.lm.gs)
#contig n50 with gs correction is more strongly significant and negative
contig.lm2.gs<-lm(abs(DiffProp)~ContigN50GSHap, data=seq.info)
summary(contig.lm2.gs)

#scaffold has a slight negative relationship, but not significant at all
scaff.lm<-lm(abs(DiffProp)~ScaffoldN50, data=seq.info)
summary(scaff.lm)
#with GS info, it becomes quite significant and negative
scaff.lm2.gs<-lm(abs(DiffProp)~ScaffoldN50GSHap, data=seq.info)
summary(scaff.lm2.gs)

#with busco, there is a relationship.  
busco.lm<-lm(abs(DiffProp)~BUSCO_Complete, data=seq.info)
summary(busco.lm)


####Making scatterplots for assessment vs DiffProp####
contig<-ggplot(seq.info, aes(x=ContigN50, y=abs(DiffProp)))+
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=0, ymax=0.1, alpha=.2, fill='black')+
  geom_point(aes(shape=Sequencing_method_1,
                 fill=Sequencing_method_1), alpha=0.8, size=2.5)+
  scale_shape_manual(values=c(21,22,24))+
  ylim(0,1.2)+
  #xlim(0,1.6e+08)+
  geom_smooth(method="lm", se=F, col="red",  linewidth=1.5)+
  #geom_hline(yintercept=0, size=2, linetype="dashed")+
  scale_fill_manual(values=seq.colors)+
  theme_minimal()+
  ylab("Prop. diff.")+
  xlab("Contig N50 (bp)")+
  annotate("text", x=2e+08, y=1.18,
           label=paste("n = 10,005"), size=4)+
  annotate("text", x=2e+08, y=1.1,
           label=expression("Adj."~R^2~"= 0.0003242"), size=4, parse=TRUE)+
  annotate("text", x=2e+08, y=1.03,
           label=paste("p = 0.0394"), size=4)+
  
  theme( axis.title = element_text(size=12, face="bold", color="black"),
         axis.text = element_text(size=10, color="black"),
         legend.title = element_blank(),
         legend.text=element_text(size=10, face="bold"),
         legend.position = "none")#+geom_smooth(method="lm")


contig.marg<-ggMarginal(contig, type="density", size=10, fill="lightgrey")

contig.chrom<-ggplot(seq.info, aes(x=ContigN50GSHap, y=abs(DiffProp)))+
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=0, ymax=0.1, alpha=.2, fill='black')+
  geom_point(aes(shape=Sequencing_method_1,
                 fill=Sequencing_method_1), alpha=0.8, size=2.5)+
  scale_shape_manual(values=c(21,22,24))+
  ylim(0,1.2)+
  #xlim(0,1.6e+08)+
  #geom_hline(yintercept=0, size=2, linetype="dashed")+
  geom_smooth(method="lm", se=F, col="red",linewidth=1.5)+
  scale_fill_manual(values=seq.colors)+
  theme_minimal()+
  ylab("Prop. diff.")+
  xlab("Contig PN50")+
  annotate("text", x=0.75, y=1.18,
           label=paste("n = 2,092"), size=4)+
  annotate("text", x=0.75, y=1.1,
           label=expression("Adj."~R^2~"= 0.031"), size=4, parse=TRUE)+
  annotate("text", x=0.75, y=1.03,
           label=paste("p < 0.00001"), size=4)+
  theme( axis.title = element_text(size=12, face="bold", color="black"),
         axis.text = element_text(size=10, color="black"),
         legend.title = element_blank(),
         legend.text=element_text(size=10, face="bold"),
         legend.position = "none")

contig.chrom.marg<-ggMarginal(contig.chrom, type="density", size=10, fill="lightgrey")

#just contig info next to eachother

scaff<-ggplot(seq.info, aes(x=ScaffoldN50,y=DiffProp))+
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=0, ymax=0.1, alpha=.2, fill='black')+
  geom_point(aes(shape=Sequencing_method_1,
                 fill=Sequencing_method_1), alpha=0.8, size=2.5)+
  scale_shape_manual(values=c(21,22,24))+
  #facet_wrap(~Kingdom, scale="free_x")+
  ylim(0,1.2)+
  #xlim(0,1.5e+09)+
  #geom_hline(yintercept=0, size=2, linetype="dashed")+
  geom_smooth(method="lm", se=F, color="red", linewidth=1.5)+
  scale_fill_manual(values=seq.colors)+theme_minimal()+
  ylab("")+
  xlab("Scaffold N50 (bp)")+
  annotate("text", x=1.5e+09, y=1.18,
           label=paste("n = 10,005"), size=4)+
  annotate("text", x=1.5e+09, y=1.1,
           label=expression("Adj."~R^2~"= 0.000"), size=4, parse=TRUE)+
  annotate("text", x=1.5e+09, y=1.03,
           label=paste("p = 0.9882"), size=4)+
  theme(
    axis.title = element_text(size=12, face="bold", color="black"),
    axis.text = element_text(size=10, color="black"),
    legend.title = element_blank(),
    legend.text=element_text(size=10, face="bold"),
    legend.position = "none")

scaff.marg<-ggMarginal(scaff, type="density", size=10, fill="lightgrey")


scaff.chrom<-ggplot(seq.info, aes(x=ScaffoldN50GSHap,y=DiffProp))+
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=0, ymax=0.1, alpha=.2, fill='black')+
  geom_point(aes(shape=Sequencing_method_1,
                 fill=Sequencing_method_1), alpha=0.8, size=2.5)+
  #facet_wrap(~Kingdom, scale="free_x")+
  scale_shape_manual(values=c(21,22,24))+
  ylim(0,1.2)+
  #xlim(0,2)+
 # geom_hline(yintercept=0, size=2, linetype="dashed")+
  geom_smooth(method="lm", se=F, color="red", linewidth=1.5)+
  scale_fill_manual(values=seq.colors)+theme_minimal()+
  ylab("")+
  annotate("text", x=1.2, y=1.18,
           label=paste("n = 2,092"), size=4)+
  annotate("text", x=1.2, y=1.1,
           label=expression("Adj."~R^2~"= 0.083"), size=4, parse=TRUE)+
  annotate("text", x=1.2, y=1.03,
           label=paste("p < 0.00001"), size=4)+
  xlab("Scaffold PN50")+
  theme(
    axis.title = element_text(size=12, face="bold", color="black"),
    axis.text = element_text(size=10, color="black"),
    legend.title = element_blank(),
    legend.text=element_text(size=10, face="bold"),
    legend.position = "none")

scaff.chrom.marg<-ggMarginal(scaff.chrom, type="density", size=10, fill="lightgrey")

busco<-ggplot(seq.info, aes(x=BUSCO_Complete,y=DiffProp))+
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=0, ymax=0.1, alpha=.2, fill='black')+
  geom_point(aes(shape=Sequencing_method_1,
                 fill=Sequencing_method_1), alpha=0.8, size=2.5)+
  scale_shape_manual(values=c(21,22,24))+
  ylim(0,1.2)+
  #geom_hline(yintercept=0, size=2, linetype="dashed")+
  geom_smooth(method="lm", se=F, color="red", linewidth=1.5)+
  scale_fill_manual(values=seq.colors)+theme_minimal()+
  ylab(" ")+
  xlab("Complete BUSCO %")+
  annotate("text", x=0.95, y=1.18,
           label=paste("n = 446"), size=4)+
  annotate("text", x=0.95, y=1.1,
           label=expression("Adj."~R^2~"= 0.016"), size=4, parse=TRUE)+
  annotate("text", x=0.95, y=1.03,
           label=paste("p = 0.00423"), size=4)+
  theme(
    axis.title = element_text(size=12, face="bold", color="black"),
    axis.text = element_text(size=10, color="black"),
    legend.title = element_blank(),
    legend.text=element_text(size=10, face="bold"),
    legend.position = "none")

busco.marg<-ggMarginal(busco, type="density", size=10, fill="lightgrey")


#busco chrom doesn't make sense to do. I'm doing it so I have the legend info to plot
busco.chrom<-ggplot(seq.info, aes(x=(BUSCO_Complete/(GS/HaploidNum)),y=DiffProp))+
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=-0.1, ymax=0.1, alpha=.2, fill='black')+
  geom_point(aes(shape=Sequencing_method_1,
                 fill=Sequencing_method_1), alpha=0.8, size=4)+
  scale_shape_manual(values=c(21,22,24))+
  ylim(-1.2,1.2)+
  geom_hline(yintercept=0, size=2, linetype="dashed")+
  scale_fill_manual(values=seq.colors)+theme_minimal()+
  ylab(" ")+
  xlab("Complete BUSCO %/(GS/Haploid Number)")+
  theme(
    axis.title = element_text(size=12, face="bold", color="black"),
    axis.text = element_text(size=10, color="black"),
    legend.title = element_blank(),
    legend.text=element_text(size=16, face="bold"),
    legend.position = "right")+
  guides(shape=guide_legend(override.aes = list(size=8)))

#this gets the legend for the combo plot
lgd<-ggpubr::get_legend(busco.chrom)
lgd.plot<-as_ggplot(lgd)

#this is the plot I want to use for the paper
ggarrange(contig, scaff, busco, contig.chrom, scaff.chrom,lgd.plot,
          ncol=3, nrow=2, labels=c("a", "b","c","d", "e",""), common.legend = FALSE,
          font.label = list(size=24, face="bold"))

ggarrange(contig.marg, scaff.marg, busco.marg, contig.chrom.marg, scaff.chrom.marg,lgd.plot,
          ncol=3, nrow=2, labels=c("a", "b","c","d", "e",""), common.legend = FALSE,
          font.label = list(size=24, face="bold"))

#info for table
length(unique(na.omit(seq.info$Species[seq.info$ScaffoldN50GSHap>0.90])))
(na.omit(seq.info$Sequencing_method_1[seq.info$ScaffoldN50GSHap>0.98]))
(na.omit(seq.info$Sequencing_method_1[seq.info$ScaffoldN50GSHap>0.98]))
na.omit(seq.info$Species[seq.info$ContigN50GSHap>0.85])

pn50tab<-subset(seq.info, subset=(ScaffoldN50GSHap>0.85 & abs(DiffProp) <0.1))
pn50tab$X.1<-NULL
pn50tab$X<-NULL

length((na.omit(seq.info$Species[seq.info$Species=="Zea mays"])))
zea<-subset(seq.info, subset=(Species=="Zea mays"))
ory<-subset(seq.info, subset=(Species=="Oryza sativa"))


brach<-subset(seq.info, subset=(Species=="Brachypodium distachyon"))
tig<-subset(pn50tab, subset=Species=="Panthera tigris")
ara<-subset(seq.info, subset=Species=="Arabidopsis thaliana")
mel<-subset(seq.info, subset=(Species=="Drosophila melanogaster"))
