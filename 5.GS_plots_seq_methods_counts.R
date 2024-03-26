library(ggplot2)
library(forcats)
library(dplyr)
library(ggpubr)
library(viridis)
library(reshape)
library(EnvStats)
library(rphylopic)
library(dunn.test)

####Making file that contains genome size and assembly information####
#read in data for assemblies
genome.seq.info.df<-read.csv("all_genome_seq_info_Feb_27.csv")
genome.seq.info.df <- genome.seq.info.df %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Organism, Assembly_Status, Sequencing_method_1, Sequencing_method_2), as.factor))




#this makes a file we already made, so I'm annotating it out
# #all data that are illumina pacbio or nanopore are seq.info
# #gs data is from genomesize.com, all animals
# #first few lines of excel file are removed, and then pg is converted, etc./
# gs<-read.csv("all.gs.data.bp.csv")
# # 
# dattest<-matrix(,nrow=length(genome.seq.info.df$Species), ncol=4)
# colnames(dattest)<-c("Species", "GS", "AssblySz", "Seq_Method")
# for(i in 1:length(genome.seq.info.df$Species)){
#   dattest[i,1]<-genome.seq.info.df$Species[i]
#   dattest[i,3]<-genome.seq.info.df$Assembly_Length_Total[i]
#   dattest[i,2]<-mean(gs$Bp[gs$Species==genome.seq.info.df$Species[i]])
#   dattest[i,4]<-genome.seq.info.df$Sequencing_method_1[i]
# }
# 
# dattest2<-as.data.frame(dattest)
# dattest2$GS<-as.numeric(dattest2$GS)
# dattest2[sapply(dattest2, is.nan)] <- NA
# dattest2$AssblySz<-as.numeric(dattest2$AssblySz)
# Diff<-(dattest2$GS-dattest2$AssblySz)
# DiffProp<-Diff/dattest2$GS
# dattest2<-data.frame(dattest2, Diff, DiffProp)
# #combine all data
# all.dat<-data.frame(genome.seq.info.df, dattest2)
# 
# 
# #make cleaned file
# dat.clean<-subset(all.dat, !is.na(GS))
# write.csv(dat.clean, "all_GS_assembly.cleaned.csv")


#####read in cleandata and format it####
dat.clean<-read.csv("all_GS_assembly.cleaned.csv")

alldat<-dat.clean
#naming all blank kingdoms and phylums as "Other"
alldat$Kingdom <- replace(alldat$Kingdom, is.na(alldat$Kingdom), "Other")
alldat$Phylum <- ifelse(alldat$Phylum == "", "Other", alldat$Phylum)
alldat <- alldat %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Organism, Assembly_Status, Sequencing_method_1, Sequencing_method_2), as.factor))
alldat$Kingdom<-factor(alldat$Kingdom, levels=c("Metazoa", "Viridiplantae", "Fungi", "Other"))

alldat_reordered <- alldat %>%
  arrange(Kingdom, Phylum)

#giving us "class" info for coloring plots
alldat_reordered$class <- cut(alldat_reordered$DiffProp, 
                              breaks = c(-Inf, -0.1, 0.1, Inf), 
                              labels = c("Below -0.1", "Between -0.1 and 0.1", "Above 0.1"))
viridis_colors <- viridis_pal(end=0.8, option="D")(3)
color_mapping<-setNames(viridis_colors, c("Smaller", "Within", "Larger"))
seq.colors<-c("#DF7C18", "#2D87A9", "#AD004F")
#how many within 10%
#levels(alldat_reordered$class)
#7310 within 10%
#length(alldat_reordered$Species[alldat_reordered$class=="Between -0.1 and 0.1"])
#6538 lower than GS estimate
#length(alldat_reordered$Species[alldat_reordered$class=="Above 0.1"])/15133


#how many species for each kingdom we have
#100 class
length(unique(alldat_reordered$Class))
#1913 that are yeast
length(alldat_reordered$Species[alldat_reordered$Species=="Saccharomyces cerevisiae"])
#1169 species of Metazoa, mean DiffProp 0.1609353
length(unique(alldat_reordered$Species[alldat_reordered$Kingdom=="Metazoa"]))
#man diffprop of 0.16 for metazoa
mean(alldat_reordered$DiffProp[alldat_reordered$Kingdom=="Metazoa"])
#mean gs of 2215.029 Mbp for metazoa
mean(alldat_reordered$GS[alldat_reordered$Kingdom=="Metazoa"])/1000000
#660 species of viridiplantae 0.3321464
length(unique(alldat_reordered$Species[alldat_reordered$Kingdom=="Viridiplantae"]))
#mean diffprop of 0.3321464 for Plants
mean(alldat_reordered$DiffProp[alldat_reordered$Kingdom=="Viridiplantae"])
#mean gs of 2587.655Mb forplants
mean(alldat_reordered$GS[alldat_reordered$Kingdom=="Viridiplantae"])/1000000
#480 species of Fungi, -0.0006581992
length(unique(alldat_reordered$Species[alldat_reordered$Kingdom=="Fungi"]))
#-0.0006582 diffprop for fungi
mean(alldat_reordered$DiffProp[alldat_reordered$Kingdom=="Fungi"])
#31.63158 mb mean gs for fungi
mean(alldat_reordered$GS[alldat_reordered$Kingdom=="Fungi"])/1000000

#dunn test is significant, all different from each other
dunn.test(x=alldat_reordered$GS, g=alldat_reordered$Kingdom, method="bonferroni")



#7770 illumina genomes, mean diffprop of 0.1241223
length(alldat_reordered$Species[alldat_reordered$Sequencing_method_1=="ILLUMINA"])
mean(alldat_reordered$DiffProp[alldat_reordered$Sequencing_method_1=="ILLUMINA"])
#609 ONT, mean diffproip of 0.1423506
length(alldat_reordered$Species[alldat_reordered$Sequencing_method_1=="OXFORD_NANOPORE"])
mean(alldat_reordered$DiffProp[alldat_reordered$Sequencing_method_1=="OXFORD_NANOPORE"])
#1627 PB, mean diffprop of 0.138995
length(alldat_reordered$Species[alldat_reordered$Sequencing_method_1=="PACBIO_SMRT"])
mean(alldat_reordered$DiffProp[alldat_reordered$Sequencing_method_1=="PACBIO_SMRT"])


####statistical differences between groups####
dunn.test(x=alldat_reordered$DiffProp, g=alldat_reordered$Kingdom, method="bonferroni")


#336 orders
length(levels(alldat$Order))
#28 phylums
length(levels(alldat$Phylum))


####Diff plot faceted by kingdoms with sample size. adding in within% barplot####
#get counts by kingdom
king.counts <- c(
  "Metazoa" = paste("Metazoa\nn =", length(alldat_reordered$Species[alldat_reordered$Kingdom=="Metazoa"])),
  "Viridiplantae" = paste("Viridiplantae\nn =",length(alldat_reordered$Species[alldat_reordered$Kingdom=="Viridiplantae"])),
  "Fungi" = paste("Fungi\nn =",length(alldat_reordered$Species[alldat_reordered$Kingdom=="Fungi"])),
  "Other" = paste("Other\nn =", length(alldat_reordered$Species[alldat_reordered$Kingdom=="Other"]))
)


king.diff.prop<-ggplot(data=alldat_reordered, aes(x=Species, y=DiffProp))+
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=-0.1, ymax=0.1, alpha=.2, fill='black')+
  geom_point(aes(color=class),size=1, alpha=0.6)+
  coord_flip()+
  theme_bw()+
  #ggtitle(paste("Proportional difference between est. GS and assembly size"))+
  scale_color_viridis_d(end=0.8, direction=-1)+theme(legend.position = "none")+
  geom_hline(aes(yintercept=0), color="red", linewidth=1.5, linetype="dashed")+
  ylim(c(-1.2,1))+
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.title=element_blank(),
        plot.title = element_text(size=14),
        axis.text.x=element_text(size=12, colour="black"),
        panel.grid.major = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face="bold", size=14, hjust=0))+
  facet_wrap(~Kingdom, scale="free_y", ncol=1, labeller=labeller(Kingdom=king.counts))

#getting count names for each kingdom for above and below
prop.kingdom<-matrix(, nrow=4, ncol=4)
colnames(prop.kingdom)<-c("Kingdom", "Within", "Above", "Below")

for(i in 1:length(levels(alldat_reordered$Kingdom))){
  testing<-subset(alldat_reordered, subset=(Kingdom==levels(alldat_reordered$Kingdom)[i]))
  prop.kingdom[i,1]<-levels(alldat_reordered$Kingdom)[i]
  prop.kingdom[i,2]<-length(testing$Kingdom[testing$class=="Between -0.1 and 0.1"])/length(testing$Kingdom)
  prop.kingdom[i,3]<-length(testing$Kingdom[testing$class=="Below -0.1"])/length(testing$Kingdom)
  prop.kingdom[i,4]<-length(testing$Kingdom[testing$class=="Above 0.1"])/length(testing$Kingdom)

}
prop.kingdom<-as.data.frame(prop.kingdom)
prop.kingdom <- prop.kingdom %>%
  mutate(across(c(Within, Above, Below), as.numeric))
prop.kingdom$Kingdom<-as.factor(prop.kingdom$Kingdom)
prop.kingdom$Kingdom<-factor(prop.kingdom$Kingdom, levels=c("Metazoa", "Viridiplantae", "Fungi", "Other"))
prop.kingdom<-melt(prop.kingdom)

viridis_colors <- viridis_pal(end=0.8, option="D")(nlevels(prop.kingdom$variable))
color_mapping<-setNames(viridis_colors, c("Below", "Within", "Above"))

king.prop.info <- c(
  "Metazoa" = paste(round(prop.kingdom$value[1]*100, digits=2),"%\n"),
  "Viridiplantae" = paste(round(prop.kingdom$value[2]*100, digits=2),"%\n"),
  "Fungi" = paste(round(prop.kingdom$value[3]*100, digits=2),"%\n"),
  "Other" = paste(round(prop.kingdom$value[4]*100, digits=2),"%\n"))

#making a barplot of how many samples are within the 10% range
king.prop.map<-ggplot(prop.kingdom, aes(x=Kingdom, y=value, fill=variable))+
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~Kingdom, scale="free_x", ncol=1,labeller=labeller(Kingdom=king.prop.info))+
  scale_fill_manual(values = color_mapping)+
  theme_bw()+
  xlab("")+
  theme(legend.position="right",
        legend.title = element_blank(),
                  axis.text=element_blank(), 
                  axis.title.y =element_blank(),
        legend.text = element_text(size=10, face="bold"),
        axis.title.x = element_text(size=12),
                  #axis.text.x=element_text(size=12, color="black"),
                  panel.grid.major = element_blank(),
                  axis.ticks = element_blank(),
                  strip.text = element_text(face="bold", size=14, hjust=0))

#figure  s4
ggarrange(king.diff.prop, king.prop.map, nrow=1,ncol=2, widths = c(1,0.4),
          align="hv") 


####Violin plot of gneome size by kingdom####
gs.violin<-ggplot(alldat_reordered, aes(x=Kingdom, y=log(GS), fill=Kingdom))+
  geom_violin()+
  coord_flip()+
  scale_fill_viridis_d(option="G",end=0.8, direction=-1)+
  theme_bw()+
  ylab("log(GS (bp))")+
  xlab("")+
  theme(legend.position="none",
        axis.text=element_text(size=12,color="black"),
        axis.title.x=element_text(size=14, face="bold", color="black"),
        axis.text.y=element_blank())

####make seq info object for when we want just illumina, pacbio, or nanopore####
seq.info<-subset(alldat, subset=(Sequencing_method_1=="ILLUMINA" | Sequencing_method_1=="PACBIO_SMRT" | Sequencing_method_1=="OXFORD_NANOPORE"))
seq.info$Sequencing_method_1<-as.character(seq.info$Sequencing_method_1)
seq.info$Sequencing_method_1<-as.factor(seq.info$Sequencing_method_1)
#change format of sequencing platform
levels(seq.info$Sequencing_method_1)
sequencers<-c("ILLUMINA", "OXFORD_NANOPORE", "PACBIO_SMRT")
levels(seq.info$Sequencing_method_1)<-c("Illumina", "Oxford Nanopore", "PacBio SMRT")
names(sequencers)<-c("Illumina", "Oxford Nanopore", "PacBio SMRT")
levels(seq.info$Sequencing_method_1)

####point plot of log(GS) vs Diffprop by kingdom
alldat_no_sac<-subset(alldat_reordered, subset=(Species!="Saccharomyces cerevisiae"))



#all data log DiffProp ~ log GS model
king.mod<-lm(abs(DiffProp)~log(GS), data=alldat_reordered)
summary(king.mod)
#all data no yeast
king.mod<-lm(abs(DiffProp)~log(GS), data=alldat_no_sac)
summary(king.mod)


#regression plot for figure 2
gs.reg<-ggplot(alldat_reordered, aes(x=log(GS), y=abs(DiffProp)))+
  geom_point(aes(fill=Kingdom),shape=21,size=2,alpha=0.6)+
  ylim(0,1)+
  theme_bw()+
  ylab("Prop. diff. from GS")+
  xlab("log(GS (bp))")+
  scale_fill_viridis_d(option="G", end=0.8, direction=-1)+
  theme(axis.text=element_text(color="black", size=12),
        axis.title=element_text(size=14, color="black", face="bold"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text=element_text(size=14, face="bold"))+
  geom_smooth(method='lm', color="red", linetype="dashed")+
  guides(fill = guide_legend(override.aes = list(shape = 22, size = 5, alpha=1)))+
  annotate("text", x=15.3, y=0.91,
           label=expression("Adj."~R^2~"= 0.08248"), size=4, parse=TRUE)+
  annotate("text", x=15.3, y=0.85,
           label=paste("p < 0.0001"), size=4)
#plot for paper figure 2
ggarrange(gs.reg, gs.violin, nrow=1, widths=c(1, 0.3), labels="auto",
          common.legend = TRUE, legend="bottom")

####Generating count data for plots####
count.data<-matrix(, nrow=length(levels(alldat_reordered$Kingdom)), ncol=2)
colnames(count.data)<-c("Kingdom", "Count")
for(i in 1:length(levels(alldat_reordered$Kingdom))){
  count.data[i,1]<-levels(alldat_reordered$Kingdom)[i]
  count.data[i,2]<-length(alldat_reordered$DiffProp[alldat_reordered$Kingdom==levels(alldat_reordered$Kingdom)[i]])
}
count.data.df<-as.data.frame(count.data)
count.data.df$Count<-as.integer(count.data.df$Count)
count.data.df$Kingdom<-as.factor(count.data.df$Kingdom)

genome.seq.info.df$Phylum<-as.factor(genome.seq.info.df$Phylum)
count.data.all<-matrix(, nrow=length(levels(genome.seq.info.df$Phylum)), ncol=2)
colnames(count.data.all)<-c("Phylum", "Count")
for(i in 1:length(levels(genome.seq.info.df$Phylum))){
  count.data.all[i,1]<-levels(genome.seq.info.df$Phylum)[i]
  count.data.all[i,2]<-length(genome.seq.info.df$Assembly_Length_Total[genome.seq.info.df$Phylum==levels(genome.seq.info.df$Phylum)[i]])
}
count.data.all.df<-as.data.frame(count.data.all)
count.data.all.df$Count<-as.integer(count.data.all.df$Count)

#how many assemblies--41358
sum(count.data.all.df$Count)
#how many assemblies have genome size estimates, 15133
sum(count.data.df$Count)

####Barplot on numbers of assemblies vs. estimates#####
#getting information for barplot on genome size estimations
count.all<-matrix(,nrow=2,ncol=3)
colnames(count.all)<-c("Data", "Status", "Count")
count.all[,1]<-"data"
count.all[1,2]<-"Without Estimate"
count.all[1,3]<-(sum(count.data.all.df$Count)-sum(count.data.df$Count))
count.all[2,2]<-"With Estimate"
count.all[2,3]<-sum(count.data.df$Count)
count.all.df<-as.data.frame(count.all)
str(count.all.df)
count.all.df$Data<-as.factor(count.all.df$Data)
count.all.df$Status<-as.factor(count.all.df$Status)
count.all.df$Count<-as.integer(count.all.df$Count)


col.values<-c("#374E8E", "#4FBBAE", "#E3B13E", "#CE4632")
#sup 1A figure
prop.count<-ggplot(data=count.all.df, aes(x=Data, y=Count, fill=Status))+
  geom_bar(position="stack", stat="identity", color="black")+
  geom_text(data=count.all.df, 
            aes(label=paste(Count),
                y=c(20000, 35000 )),size=5, color=c("black", "white"))+
  ylab("Number of Assemblies")+
  #ylim(0,17000)+
  theme_minimal()+
  #scale_fill_viridis_d(begin=0.5, end=0.9, direction=-1)+
  #scale_fill_manual(values=natparks.pals("Arches", 10))+
  scale_fill_manual(values=c(col.values[1], "grey"))+
  theme(legend.title = element_blank(),
        legend.text=element_text(face="bold", size=12),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14, face="bold"),
        legend.position = "bottom")

col.pal <- viridis_pal(end=0.8, option="G")(nlevels(count.data.df$Kingdom))
color_mapping.2<-setNames(col.pal, c("Other", "Fungi", "Viridiplantae", "Metazoa"))

#assemblies by kingdom
gs.assem<-ggplot(data=count.data.df, aes(x=fct_reorder(Kingdom, Count), y=Count))+
  geom_bar(aes(fill=Kingdom), 
           #fill=col.values, 
           color="black",
           stat="identity")+
  coord_flip()+
  theme_minimal()+
  scale_fill_manual(values=color_mapping.2)+
  ylab("Number of Assemblies with available GS estimate")+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size=12, color="black", face="bold"),
        axis.text=element_text(size=12, color="black", face="bold"),
        plot.title = element_text(size=18, face="bold"),
        legend.position = "none")

#possible figure for paper. may be redundant with table  
ggarrange(prop.count, gs.assem, 
          ncol=2, labels="auto",
          common.legend = FALSE,
          font.label=list(size=22, face="bold"),
          widths =c(0.3,1),
          align="v")


####Point plot for proportional differences####
#change format of sequencing platform
seq.method.tally <- seq.info %>% group_by(Sequencing_method_1) %>% tally()
seq.method.tally$n[1]
seq.counts <- c(
  "Illumina" = paste("Illumina\nn =", seq.method.tally$n[1]),
  "Oxford Nanopore" = paste("Oxford Nanopore\nn =",seq.method.tally$n[2]),
  "PacBio SMRT" = paste("PacBio SMRT\nn =",seq.method.tally$n[3])
  )


seq.info$Phylum<-as.character(seq.info$Phylum)
seq.info$Phylum<-as.factor(seq.info$Phylum)

seq.info$class <- cut(seq.info$DiffProp, 
                  breaks = c(-Inf, -0.1, 0.1, Inf), 
                  labels = c("Below -0.1", "Between -0.1 and 0.1", "Above 0.1"))
viridis_colors <- viridis_pal(end=0.8, option="D")(3)
color_mapping<-setNames(viridis_colors, c("Smaller", "Within", "Larger"))
seq.colors<-c("#DF7C18", "#2D87A9", "#AD004F")


seq.info.2<-subset(seq.info, subset=(seq.info$Kingdom!="Other"))

#gridded plot for seq method and kingdom.  Maybe useful for supplemental?
#figure s3
Prop.data<-ggplot(data=seq.info.2, aes(x=Species, y=DiffProp))+
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=-0.1, ymax=0.1, alpha=.2, fill='black')+
  geom_point(aes(fill=class, shape=Assembly_Status),size=2, alpha=0.6)+
  geom_hline(aes(yintercept=0), color="red", linewidth=1.5, linetype="dashed")+
  scale_shape_manual(values=c(21, 22, 24, 25), name="Assembly Status")+
  coord_flip()+
  ylab("Proportional Diff. from GS")+
  ylim(c(-1.2,1))+
  facet_grid(Sequencing_method_1~Kingdom, scale="free")+
  theme_bw()+
  scale_fill_viridis_d(end=0.8, direction=-1)+
  guides(fill="none",shape=guide_legend(override.aes=list(size=5, fill="grey")))+
  theme(axis.text.y=element_blank(), 
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=10, color="black", face="bold"),
        plot.title = element_text(size=12),
        axis.text.x=element_text(size=10, colour="black"),
        panel.grid.major = element_blank(), 
        strip.text = element_text(face="bold", size=10, hjust=0),
        legend.text=element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"),
        plot.margin=margin(0.75, 0.5, 0.5,1, "cm"),
        legend.position = "bottom",
        axis.ticks.y=element_blank())


#getting count names for each seq.method for above and below
prop.seq<-matrix(, nrow=3, ncol=4)
colnames(prop.seq)<-c("Method", "Within", "Above", "Below")

for(i in 1:length(levels(seq.info$Sequencing_method_1))){
  testing<-subset(seq.info, subset=(Sequencing_method_1==levels(seq.info$Sequencing_method_1)[i]))
  prop.seq[i,1]<-levels(seq.info$Sequencing_method_1)[i]
  prop.seq[i,2]<-length(testing$Sequencing_method_1[testing$class=="Between -0.1 and 0.1"])/length(testing$Sequencing_method_1)
  prop.seq[i,3]<-length(testing$Sequencing_method_1[testing$class=="Below -0.1"])/length(testing$Sequencing_method_1)
  prop.seq[i,4]<-length(testing$Sequencing_method_1[testing$class=="Above 0.1"])/length(testing$Sequencing_method_1)
  
}
prop.seq<-as.data.frame(prop.seq)
prop.seq <- prop.seq %>%
  mutate(across(c(Within, Above, Below), as.numeric))
prop.seq$Method<-as.factor(prop.seq$Method)
prop.seq<-melt(prop.seq)
viridis_colors <- viridis_pal(end=0.8, option="D")(nlevels(prop.seq$variable))
color_mapping<-setNames(viridis_colors, c("Below", "Within", "Above"))

seq.prop.info <- c(
  "Illumina" = paste("Illumina\n",round(prop.seq$value[1]*100, digits=2),"%"),
  "Oxford Nanopore" = paste("Oxford Nanopore\n",round(prop.seq$value[2]*100, digits=2),"%"),
  "PacBio SMRT" = paste("PacBio SMRT\n", round(prop.seq$value[3]*100, digits=2),"%"))

####what if we remove fungi####
#getting count names for each seq.method for above and below
seq.info.nofun<-subset(seq.info, subset=(Kingdom!="Fungi"))

#3428 illumina genomes with no fungus, mean diffprop of 0.2743664
length(seq.info.nofun$Species[seq.info.nofun$Sequencing_method_1=="Illumina"])
mean(seq.info.nofun$DiffProp[seq.info.nofun$Sequencing_method_1=="Illumina"])
#451 ONT, mean diffproip of 0.2134537
length(seq.info.nofun$Species[seq.info.nofun$Sequencing_method_1=="Oxford Nanopore"])
mean(seq.info.nofun$DiffProp[seq.info.nofun$Sequencing_method_1=="Oxford Nanopore"])
#1367 PB, mean diffprop of 0.1722424
length(seq.info.nofun$Species[seq.info.nofun$Sequencing_method_1=="PacBio SMRT"])
mean(seq.info.nofun$DiffProp[seq.info.nofun$Sequencing_method_1=="PacBio SMRT"])

#significant KW test differences by method excluding fungi. Illumina is different, but ONT and Nanopore arent' different
#positive value for comparison with illumina. Illumina has "greater" value (distribtuion higher?)
#p= 0 for illumina comparisons, p=0.9399 for ONT vs Pacbio
#kw p =0
dunn.test(x=seq.info.nofun$DiffProp, g=seq.info.nofun$Sequencing_method_1, kw=TRUE, method="bonferroni")

#significant KW test differences by method including fungi. Illumina is different, ONT and pb are differnet
#- value for comparison, illumina has lower values? closer to zero
#p= 0 for illumina comparisons, p=0.0.0245 for ONT vs Pacbio
#kw p =0
dunn.test(x=seq.info$DiffProp, g=seq.info$Sequencing_method_1, kw=TRUE, method="bonferroni")

seq.colors<-c("#DF7C18", "#2D87A9", "#AD004F")


#getting proportion data without fungi  
prop.seq2<-matrix(, nrow=3, ncol=4)
colnames(prop.seq2)<-c("Method", "Within", "Above", "Below")
for(i in 1:length(levels(seq.info.nofun$Sequencing_method_1))){
  testing<-subset(seq.info.nofun, subset=(Sequencing_method_1==levels(seq.info.nofun$Sequencing_method_1)[i]))
  prop.seq2[i,1]<-levels(seq.info.nofun$Sequencing_method_1)[i]
  prop.seq2[i,2]<-length(testing$Sequencing_method_1[testing$class=="Between -0.1 and 0.1"])/length(testing$Sequencing_method_1)
  prop.seq2[i,3]<-length(testing$Sequencing_method_1[testing$class=="Below -0.1"])/length(testing$Sequencing_method_1)
  prop.seq2[i,4]<-length(testing$Sequencing_method_1[testing$class=="Above 0.1"])/length(testing$Sequencing_method_1)
  
}
#formating data
prop.seq2<-as.data.frame(prop.seq2)
prop.seq2 <- prop.seq2 %>%
  mutate(across(c(Within, Above, Below), as.numeric))
prop.seq2$Method<-as.factor(prop.seq2$Method)
prop.seq2<-melt(prop.seq2)
viridis_colors <- viridis_pal(end=0.8, option="D")(nlevels(prop.seq2$variable))
color_mapping<-setNames(viridis_colors, c("Below", "Within", "Above"))

#proportion within 10% excluding fungi
seq.prop.info2 <- c(
  "Illumina" = paste("Illumina\n",round(prop.seq2$value[1]*100, digits=2),"%"),
  "Oxford Nanopore" = paste("Oxford Nanopore\n",round(prop.seq2$value[2]*100, digits=2),"%"),
  "PacBio SMRT" = paste("PacBio SMRT\n", round(prop.seq2$value[3]*100, digits=2),"%"))

#density that has overlap
density.nofun.2<-ggplot(seq.info.nofun, aes(x=DiffProp, fill=Sequencing_method_1))+
  annotate('rect', xmin=-0.1, xmax=0.1, ymin=-Inf, ymax=Inf, alpha=.2, fill='black')+
  geom_density(alpha=0.7)+
  xlab("Proportional Difference from GS")+
  geom_vline(xintercept=0, col="red", linetype="dashed")+
  scale_fill_manual(values=c(seq.colors))+
  xlim(-1.2,1)+
  theme_bw()+
  ggtitle("b.  Excluding Fungi")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.text = element_text(size=12, color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=12, face="bold"),
        legend.position="bottom",
        strip.text=element_text(hjust=0, face="bold", color="black", size=15),
        panel.grid.major = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(size=20, face="bold"))+
  guides(fill = guide_legend(nrow = 1))


density.all2<-ggplot(seq.info, aes(x=DiffProp, fill=Sequencing_method_1))+
  annotate('rect', xmin=-0.1, xmax=0.1, ymin=-Inf, ymax=Inf, alpha=.2, fill='black')+
  geom_density(alpha=0.7)+
  xlab("Proportional Difference from GS")+
  geom_vline(xintercept=0, col="red", linetype="dashed")+
  scale_fill_manual(values=c(seq.colors))+
  xlim(-1.2,1)+
  theme_bw()+
  ggtitle("a.  Including Fungi")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.text = element_text(size=12, color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=12, face="bold"),
        legend.position="bottom",
        strip.text=element_text(hjust=0, face="bold", color="black", size=15),
        panel.grid.major = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        plot.title = element_text(size=20, face="bold"))+
  guides(fill = guide_legend(nrow = 1))


seq.prop.map.all<-ggplot(prop.seq, aes(x=Method, y=value, fill=variable))+
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~Method, scale="free_x", ncol=3,labeller=labeller(Method=seq.prop.info))+
  xlab("Prop. of assemblies")+
  scale_fill_manual(values = color_mapping)+
  scale_x_discrete(labels=c("","",""))+
  ggtitle("")+
  theme_bw()+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.text.y=element_blank(),
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y =element_blank(),
        legend.text = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        axis.ticks.y  = element_blank(),
        strip.text = element_text(face="bold", size=14, hjust=0),
        #plot.title = element_text(size=18, face = "bold"),
        plot.title = element_text(size=20, face="bold"))
  
seq.prop.map.2.nofun<-ggplot(prop.seq2, aes(x=Method, y=value, fill=variable))+
  geom_bar(position="stack", stat="identity")+
  facet_wrap(~Method, scale="free_x", ncol=3,labeller=labeller(Method=seq.prop.info2))+
  scale_fill_manual(values = color_mapping)+
  theme_bw()+
  scale_x_discrete(labels=c("","",""))+
  xlab("Prop. of assemblies")+
  ggtitle("")+
  theme(legend.position="bottom",
        legend.title = element_blank(),
        axis.text.y=element_blank(), 
        axis.title.x = element_text(size=14, face="bold", color="black"),
        axis.title.y =element_blank(),
        legend.text = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(face="bold", size=14, hjust=0),
        #plot.title = element_text(size=18, face = "bold"),
        plot.title = element_text(size=20, face="bold"))
  
#figure 3
ggarrange(density.all2, seq.prop.map.all, density.nofun.2, seq.prop.map.2.nofun,
          ncol=2, nrow=2, widths=c(1,0.6))
