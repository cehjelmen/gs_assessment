library(ggplot2)
library(forcats)
library(dplyr)
library(ggpubr)
library(viridis)
library(reshape)
library(EnvStats)
library(cowplot)
#read in data
dat.clean<-read.csv("all_GS_assembly.cleaned.csv")

alldat<-dat.clean
#prepping info for the analysis
alldat <- alldat %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Organism, Assembly_Status, Sequencing_method_1, Sequencing_method_2), as.factor))
#removing extranenous info
alldat$X.1<-NULL
alldat$X<-NULL
alldat[8:32]<-NULL
alldat[9:12]<-NULL

#make column in order to get GS info for each species. one value per species (mean across records)
unique.gs<-matrix(,nrow=length(unique(alldat$Species)), ncol=ncol(alldat))
colnames(unique.gs)<-colnames(alldat)
for(i in 1:length(unique(alldat$Species))){
  unique.gs[i,1]<-as.character(unique(alldat$Kingdom[alldat$Species==unique(alldat$Species)[i]]))
  unique.gs[i,2]<-as.character(unique(alldat$Phylum[alldat$Species==unique(alldat$Species)[i]]))
  unique.gs[i,3]<-as.character(unique(alldat$Class[alldat$Species==unique(alldat$Species)[i]]))
  unique.gs[i,4]<-as.character(unique(alldat$Order[alldat$Species==unique(alldat$Species)[i]]))
  unique.gs[i,5]<-as.character(unique(alldat$Family[alldat$Species==unique(alldat$Species)[i]]))
  unique.gs[i,6]<-as.character(unique(alldat$Genus[alldat$Species==unique(alldat$Species)[i]]))
  unique.gs[i,7]<-unique(alldat$Species)[i]
  unique.gs[i,8]<-as.numeric(mean(alldat$GS[alldat$Species==unique(alldat$Species)[i]]))
}
unique.gs.df<-as.data.frame(unique.gs)
unique.gs.df <- unique.gs.df %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus), as.factor))
unique.gs.df$GS<-as.numeric(unique.gs.df$GS)


#here we are getting the GS/mean estimate for that taxonomic rank
genus.testing<-matrix(,nrow=length(unique.gs.df$Kingdom), ncol=18)
colnames(genus.testing)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", 
                           "Species", "GS", "Genus.Mean","Genus.Count", 
                           "Family.Mean","Family.Count", "Order.Mean","Order.Count", 
                           "Class.Mean","Class.Count", "Phylum.Mean", "Phylum.Count")
for(i in 1:length(unique.gs.df$Kingdom)){
  genus.testing[i,1]<-as.character(unique.gs.df$Kingdom[i])
  genus.testing[i,2]<-as.character(unique.gs.df$Phylum[i])
  genus.testing[i,3]<-as.character(unique.gs.df$Class[i])
  genus.testing[i,4]<-as.character(unique.gs.df$Order[i])
  genus.testing[i,5]<-as.character(unique.gs.df$Family[i])
  genus.testing[i,6]<-as.character(unique.gs.df$Genus[i])
  genus.testing[i,7]<-as.character(unique.gs.df$Species[i])
  genus.testing[i,8]<-unique.gs.df$GS[i]
  #genus mean
  genus.testing[i,9]<-as.numeric(genus.testing[i,8])/(mean(unique.gs.df$GS[unique.gs.df$Genus==genus.testing[i,6]]))
  #genus.count
  genus.testing[i,10]<-length((unique.gs.df$GS[unique.gs.df$Genus==genus.testing[i,6]]))
  #family mean
  genus.testing[i,11]<-as.numeric(genus.testing[i,8])/(mean(unique.gs.df$GS[unique.gs.df$Family==genus.testing[i,5]]))
  #family count
  genus.testing[i,12]<-length((unique.gs.df$GS[unique.gs.df$Family==genus.testing[i,5]]))
  #order mean
  genus.testing[i,13]<-as.numeric(genus.testing[i,8])/(mean(unique.gs.df$GS[unique.gs.df$Order==genus.testing[i,4]]))
  #order count
  genus.testing[i,14]<-length((unique.gs.df$GS[unique.gs.df$Order==genus.testing[i,4]]))
  #Class Mean
  genus.testing[i,15]<-as.numeric(genus.testing[i,8])/(mean(unique.gs.df$GS[unique.gs.df$Class==genus.testing[i,3]]))
  #class count
  genus.testing[i,16]<-length((unique.gs.df$GS[unique.gs.df$Class==genus.testing[i,3]]))
  #Phylum Mean
  genus.testing[i,17]<-as.numeric(genus.testing[i,8])/(mean(unique.gs.df$GS[unique.gs.df$Phylum==genus.testing[i,2]]))
  #phylum count
  genus.testing[i,18]<-length((unique.gs.df$GS[unique.gs.df$Phylum==genus.testing[i,2]]))
  
}
genus.testing.df<-as.data.frame(genus.testing)

genus.testing.df <- genus.testing.df %>%
  mutate(across(c(GS, Genus.Mean, Family.Mean, Order.Mean, Class.Mean, Phylum.Mean), as.numeric))
genus.testing.df <- genus.testing.df %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus), as.factor))
genus.testing.df <- genus.testing.df %>%
  mutate(across(c(Phylum.Count, Class.Count, Order.Count, Family.Count, Genus.Count), as.integer))

#making a genus.dat file, which makes sure we only look at taxonomic ranks which have at least 5
genus.dat<-subset(genus.testing.df, subset=(Genus.Count>4))
genus.dat$Phylum<-as.character(genus.dat$Phylum)
genus.dat$Phylum<-as.factor(genus.dat$Phylum)
genus.dat$Class<-as.character(genus.dat$Class)
genus.dat$Class<-as.factor(genus.dat$Class)
genus.dat$Order<-as.character(genus.dat$Order)
genus.dat$Order<-as.factor(genus.dat$Order)
genus.dat$Family<-as.character(genus.dat$Family)
genus.dat$Family<-as.factor(genus.dat$Family)
genus.dat$Genus<-as.character(genus.dat$Genus)
genus.dat$Genus<-as.factor(genus.dat$Genus)

#str(genus.dat)
#makign these for later
meta.unique.gs<-subset(genus.dat, subset=(Kingdom=="Metazoa"))
virid.unique.gs<-subset(genus.dat, subset=(Kingdom=="Viridiplantae"))
fungi.unique.gs<-subset(genus.dat, subset=(Kingdom=="Fungi"))

#genus.dat$Phylum<-NULL
genus.dat$Class<-NULL
genus.dat$Order<-NULL
genus.dat$Family<-NULL
genus.dat$Genus<-NULL
genus.dat$Species<-NULL
genus.dat$Phylum.Count<-NULL
genus.dat$Order.Count<-NULL
genus.dat$Class.Count<-NULL
genus.dat$Family.Count<-NULL
genus.dat$Genus.Count<-NULL
genus.dat$GS<-NULL

#melt the data wiht reshape package for use in later plots
var.gs.dat<-melt(genus.dat, id.vars = c("Kingdom", "Phylum"))
str(var.gs.dat)


#renaming levels
levels(var.gs.dat$variable)<-c("Genus", "Family", "Order", "Class", "Phylum")
#ordering factors for plot
var.gs.dat$variable<-factor(var.gs.dat$variable, levels =c("Phylum", "Class", "Order", "Family", "Genus"))
#we have 2535 records
length(var.gs.dat$Kingdom)

#density plots
col.pal <- viridis_pal(end=0.8, option="G")(nlevels(var.gs.dat$Kingdom))
color_mapping<-setNames(col.pal, c("Fungi", "Viridiplantae", "Metazoa"))
#levels(var.gs.dat$Phylum)
levels(var.gs.dat$Kingdom)


####subsetting data by kingdom####
metazoa<-subset(var.gs.dat, subset=(Kingdom=="Metazoa"))
viridiplantae<-subset(var.gs.dat, subset=(Kingdom=="Viridiplantae"))
fungi<-subset(var.gs.dat, subset=(Kingdom=="Fungi"))

####faceted density plots for kingdoms####
#colors and plot for metazoans
meta.col<-c("#60CEAC", "#1D5E49", "#D0F1E6", "white")
meta.tax<-ggplot(metazoa, aes(x=value, fill=Phylum))+
  ggtitle("Metazoa")+
  geom_density(alpha=0.7)+
  facet_wrap(~variable, ncol=1)+
  xlim(0,4)+
  xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=meta.col)+
  theme(legend.position="bottom",
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=10, face="bold"),
        axis.text=element_text(size=8, color="black"),
        strip.text = element_text(hjust=0, size=12, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(ncol = 2))
  
#plants only had one phylum, so color is built in
plant.tax<-ggplot(viridiplantae, aes(x=value, fill=Phylum))+
  ggtitle("Viridiplantae")+
  geom_density(alpha=0.7)+
  facet_wrap(~variable, ncol=1)+
  xlim(0,4)+
  xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values="#395D9CFF")+
  theme(legend.position="bottom",
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=10, face="bold"),
        axis.text=element_text(size=8, color="black"),
        strip.text = element_text(hjust=0, size=12, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(ncol = 2))


#fungus plot and colors
fun.col<-c("#41346F", "#7866B7", "#DAD5EB")                      
fun.tax<-ggplot(fungi, aes(x=value, fill=Phylum))+
  ggtitle("Fungi")+
  geom_density(alpha=0.7)+
  facet_wrap(~variable, ncol=1)+
  xlim(0,4)+
  xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=fun.col)+
  theme(legend.position="bottom",
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=10, face="bold"),
        axis.text=element_text(size=8, color="black"),
        strip.text = element_text(hjust=0, size=12, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(ncol = 2))
#ggarrange(meta.tax, plant.tax, fun.tax,ncol=3, align="hv")


####Larger plots with rnages and 95% confidence intervals

####metazoa plot####
meta.unique.gs$Phylum<-as.character(meta.unique.gs$Phylum)
meta.unique.gs$Phylum<-as.factor(meta.unique.gs$Phylum)

meta.line<-matrix(,nrow=length(levels(meta.unique.gs$Phylum)), ncol=27)
colnames(meta.line)<-c("Kingdom", "Phylum",
                       "Genus.Mean", "Genus.Lower", "Genus.Upper", "Genus.95L", "Genus.95U", 
                       "Family.Mean", "Family.Lower", "Family.Upper", "Family.95L", "Family.95U", 
                       "Order.Mean", "Order.Lower", "Order.Upper", "Order.95L", "Order.95U", 
                       "Class.Mean", "Class.Lower", "Class.Upper", "Class.95L", "Class.95U", 
                       "Phylum.Mean", "Phylum.Lower", "Phylum.Upper", "Phylum.95L", "Phylum.95U")
for(i in 1:length(levels(meta.unique.gs$Phylum))){
  meta.line[i,1]<-"Metazoa"
  meta.line[i,2]<-levels(meta.unique.gs$Phylum)[i]
  #genus info
  meta.line[i,3]<-mean(meta.unique.gs$Genus.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])
  meta.line[i,4]<-range(meta.unique.gs$Genus.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])[1]
  meta.line[i,5]<-range(meta.unique.gs$Genus.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])[2]
  output<-t.test(meta.unique.gs$Genus.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  meta.line[i,6]<-output$conf.int[1]
  meta.line[i,7]<-output$conf.int[2]
  #family info
  meta.line[i,8]<-mean(meta.unique.gs$Family.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])
  meta.line[i,9]<-range(meta.unique.gs$Family.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])[1]
  meta.line[i,10]<-range(meta.unique.gs$Family.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])[2]
  output<-t.test(meta.unique.gs$Family.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  meta.line[i,11]<-output$conf.int[1]
  meta.line[i,12]<-output$conf.int[2]
  #Order info
  meta.line[i,13]<-mean(meta.unique.gs$Order.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])
  meta.line[i,14]<-range(meta.unique.gs$Order.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])[1]
  meta.line[i,15]<-range(meta.unique.gs$Order.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])[2]
  output<-t.test(meta.unique.gs$Order.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  meta.line[i,16]<-output$conf.int[1]
  meta.line[i,17]<-output$conf.int[2]
  #class info
  meta.line[i,18]<-mean(meta.unique.gs$Class.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])
  meta.line[i,19]<-range(meta.unique.gs$Class.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])[1]
  meta.line[i,20]<-range(meta.unique.gs$Class.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])[2]
  output<-t.test(meta.unique.gs$Class.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  meta.line[i,21]<-output$conf.int[1]
  meta.line[i,22]<-output$conf.int[2]
  #phylum info
  meta.line[i,23]<-mean(meta.unique.gs$Phylum.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])
  meta.line[i,24]<-range(meta.unique.gs$Phylum.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])[1]
  meta.line[i,25]<-range(meta.unique.gs$Phylum.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]])[2]
  output<-t.test(meta.unique.gs$Phylum.Mean[meta.unique.gs$Phylum==levels(meta.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  meta.line[i,26]<-output$conf.int[1]
  meta.line[i,27]<-output$conf.int[2]
}
meta.line<-as.data.frame(meta.line)
meta.line <- meta.line %>%
  mutate(across(c("Kingdom", "Phylum"), as.factor))


meta.line <- meta.line %>%
  mutate(across(c("Genus.Mean", "Genus.Lower", "Genus.Upper", "Genus.95L", "Genus.95U", 
                       "Family.Mean", "Family.Lower", "Family.Upper", "Family.95L", "Family.95U", 
                       "Order.Mean", "Order.Lower", "Order.Upper", "Order.95L", "Order.95U", 
                       "Class.Mean", "Class.Lower", "Class.Upper", "Class.95L", "Class.95U", 
                       "Phylum.Mean", "Phylum.Lower", "Phylum.Upper", "Phylum.95L", "Phylum.95U"), as.numeric))

#genus info
meta.col<-c("#60CEAC", "#1D5E49", "#D0F1E6", "white")
meta.gen<-ggplot(meta.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Genus.Lower-0.0125, ymax=Genus.Upper+0.0125, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Genus.Lower, ymax=Genus.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Genus.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Genus.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Genus.Mean), size=1.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,6.2)+
  scale_color_manual(values=meta.col)+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=12, face="bold"))
meta.col<-c("#60CEAC", "#1D5E49", "#D0F1E6", "white")
meta.gen.dat<-subset(metazoa, subset=(variable=="Genus"))
meta.gen.dat$title<-"Genus"

meta.gen.dens<-ggplot(meta.gen.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Genus")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,6.2)+
  ylim(0,8)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=meta.col)+
  theme(legend.position="bottom",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

meta.gen.tot<-ggarrange(meta.gen.dens, meta.gen, nrow=2, ncol=1, common.legend = TRUE, legend="bottom", heights=c(1,0.4), align="v")

#family info
meta.fam.dat<-subset(metazoa, subset=(variable=="Family"))
meta.fam.dat$title<-"Family"

meta.fam.dens<-ggplot(meta.fam.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Family")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,6.2)+
  ylim(0,8)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=meta.col)+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

meta.fam<-ggplot(meta.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Family.Lower-0.0125, ymax=Family.Upper+0.0125, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Family.Lower, ymax=Family.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Family.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Family.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Family.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,6.2)+
  scale_color_manual(values=meta.col)+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )  

meta.fam.tot<-ggarrange(meta.fam.dens, meta.fam, nrow=2, ncol=1, common.legend = FALSE, heights=c(1,0.4), align="v")

#Order info
meta.ord.dat<-subset(metazoa, subset=(variable=="Order"))
meta.ord.dat$title<-"Order"

meta.ord.dens<-ggplot(meta.ord.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Order")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,6.2)+
  ylim(0,8)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=meta.col)+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))


meta.ord<-ggplot(meta.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Order.Lower-0.0125, ymax=Order.Upper+0.0125, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Order.Lower, ymax=Order.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Order.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Order.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Order.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,6.2)+
  scale_color_manual(values=meta.col)+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )   

meta.ord.tot<-ggarrange(meta.ord.dens, meta.ord, nrow=2, ncol=1, common.legend = FALSE, heights=c(1,0.4), align="v")

#Class info
meta.class.dat<-subset(metazoa, subset=(variable=="Class"))
meta.class.dat$title<-"Class"

meta.class.dens<-ggplot(meta.class.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Class")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,6.2)+
  ylim(0,8)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=meta.col)+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

meta.class<-ggplot(meta.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Class.Lower-0.0125, ymax=Class.Upper+0.0125, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Class.Lower, ymax=Class.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Class.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Class.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Class.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,6.2)+
  scale_color_manual(values=meta.col)+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )   

meta.class.tot<-ggarrange(meta.class.dens, meta.class, nrow=2, ncol=1, common.legend =FALSE, heights=c(1,0.4), align="v")


#phylum info
meta.phy.dat<-subset(metazoa, subset=(variable=="Phylum"))
meta.phy.dat$title<-"Phylum"

meta.phy.dens<-ggplot(meta.phy.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  ggtitle("Metazoa")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,6.2)+
  ylim(0,8)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=meta.col)+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=18, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

meta.phy<-ggplot(meta.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Phylum.Lower-0.0125, ymax=Phylum.Upper+0.0125, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Phylum.Lower, ymax=Phylum.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Phylum.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Phylum.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Phylum.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,6.2)+
  scale_color_manual(values=meta.col)+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
        )  

meta.phy.tot<-ggarrange(meta.phy.dens, meta.phy, nrow=2, ncol=1, common.legend = FALSE, heights=c(1,0.4), align="v")

meta.dist<-ggarrange(meta.phy.tot, meta.class.tot, meta.ord.tot, meta.fam.tot, meta.gen.tot,
          ncol=1,heights=c(0.85,0.75,0.75,0.75,1),
          common.legend = TRUE)



####Fungi time####
#fungi.unique.gs

fungi.unique.gs$Phylum<-as.character(fungi.unique.gs$Phylum)
fungi.unique.gs$Phylum<-as.factor(fungi.unique.gs$Phylum)

fungi.line<-matrix(,nrow=length(levels(fungi.unique.gs$Phylum)), ncol=27)
colnames(fungi.line)<-c("Kingdom", "Phylum",
                       "Genus.Mean", "Genus.Lower", "Genus.Upper", "Genus.95L", "Genus.95U", 
                       "Family.Mean", "Family.Lower", "Family.Upper", "Family.95L", "Family.95U", 
                       "Order.Mean", "Order.Lower", "Order.Upper", "Order.95L", "Order.95U", 
                       "Class.Mean", "Class.Lower", "Class.Upper", "Class.95L", "Class.95U", 
                       "Phylum.Mean", "Phylum.Lower", "Phylum.Upper", "Phylum.95L", "Phylum.95U")
for(i in 1:length(levels(fungi.unique.gs$Phylum))){
  fungi.line[i,1]<-"Fungi"
  fungi.line[i,2]<-levels(fungi.unique.gs$Phylum)[i]
  #genus info
  fungi.line[i,3]<-mean(fungi.unique.gs$Genus.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])
  fungi.line[i,4]<-range(fungi.unique.gs$Genus.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])[1]
  fungi.line[i,5]<-range(fungi.unique.gs$Genus.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])[2]
  output<-t.test(fungi.unique.gs$Genus.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  fungi.line[i,6]<-output$conf.int[1]
  fungi.line[i,7]<-output$conf.int[2]
  #family info
  fungi.line[i,8]<-mean(fungi.unique.gs$Family.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])
  fungi.line[i,9]<-range(fungi.unique.gs$Family.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])[1]
  fungi.line[i,10]<-range(fungi.unique.gs$Family.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])[2]
  output<-t.test(fungi.unique.gs$Family.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  fungi.line[i,11]<-output$conf.int[1]
  fungi.line[i,12]<-output$conf.int[2]
  #Order info
  fungi.line[i,13]<-mean(fungi.unique.gs$Order.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])
  fungi.line[i,14]<-range(fungi.unique.gs$Order.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])[1]
  fungi.line[i,15]<-range(fungi.unique.gs$Order.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])[2]
  output<-t.test(fungi.unique.gs$Order.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  fungi.line[i,16]<-output$conf.int[1]
  fungi.line[i,17]<-output$conf.int[2]
  #class info
  fungi.line[i,18]<-mean(fungi.unique.gs$Class.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])
  fungi.line[i,19]<-range(fungi.unique.gs$Class.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])[1]
  fungi.line[i,20]<-range(fungi.unique.gs$Class.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])[2]
  output<-t.test(fungi.unique.gs$Class.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  fungi.line[i,21]<-output$conf.int[1]
  fungi.line[i,22]<-output$conf.int[2]
  #phylum info
  fungi.line[i,23]<-mean(fungi.unique.gs$Phylum.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])
  fungi.line[i,24]<-range(fungi.unique.gs$Phylum.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])[1]
  fungi.line[i,25]<-range(fungi.unique.gs$Phylum.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]])[2]
  output<-t.test(fungi.unique.gs$Phylum.Mean[fungi.unique.gs$Phylum==levels(fungi.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  fungi.line[i,26]<-output$conf.int[1]
  fungi.line[i,27]<-output$conf.int[2]
}
fungi.line<-as.data.frame(fungi.line)
fungi.line <- fungi.line %>%
  mutate(across(c("Kingdom", "Phylum"), as.factor))


fungi.line <- fungi.line %>%
  mutate(across(c("Genus.Mean", "Genus.Lower", "Genus.Upper", "Genus.95L", "Genus.95U", 
                  "Family.Mean", "Family.Lower", "Family.Upper", "Family.95L", "Family.95U", 
                  "Order.Mean", "Order.Lower", "Order.Upper", "Order.95L", "Order.95U", 
                  "Class.Mean", "Class.Lower", "Class.Upper", "Class.95L", "Class.95U", 
                  "Phylum.Mean", "Phylum.Lower", "Phylum.Upper", "Phylum.95L", "Phylum.95U"), as.numeric))


fun.col<-c("#41346F", "#7866B7", "#DAD5EB")                      
fungi.gen<-ggplot(fungi.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Genus.Lower-0.01, ymax=Genus.Upper+0.01, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Genus.Lower, ymax=Genus.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Genus.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Genus.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Genus.Mean), size=1.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,3.7)+
  scale_color_manual(values=fun.col)+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=12, face="bold"))
fungi.gen.dat<-subset(fungi, subset=(variable=="Genus"))
fungi.gen.dat$title<-"Genus"
fungi.gen.dens<-ggplot(fungi.gen.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Genus")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,3.7)+
  ylim(0,6)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=fun.col)+
  theme(legend.position="bottom",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

fungi.gen.tot<-ggarrange(fungi.gen.dens, fungi.gen, nrow=2, ncol=1, common.legend = TRUE, legend="bottom", heights=c(1,0.4), align="v")

#family info
fungi.fam.dat<-subset(fungi, subset=(variable=="Family"))
fungi.fam.dat$title<-"Family"
fungi.fam.dens<-ggplot(fungi.fam.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Family")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,3.7)+
  ylim(0,6)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=fun.col)+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

fungi.fam<-ggplot(fungi.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Family.Lower-0.01, ymax=Family.Upper+0.01, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Family.Lower, ymax=Family.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Family.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Family.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Family.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,3.7)+
  scale_color_manual(values=fun.col)+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )  

fungi.fam.tot<-ggarrange(fungi.fam.dens, fungi.fam, nrow=2, ncol=1, common.legend = FALSE, heights=c(1,0.4), align="v")

#Order info
fungi.ord.dat<-subset(fungi, subset=(variable=="Order"))
fungi.ord.dat$title<-"Order"
fungi.ord.dens<-ggplot(fungi.ord.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Order")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,3.7)+
  ylim(0,6)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=fun.col)+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))


fungi.ord<-ggplot(fungi.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Order.Lower-0.01, ymax=Order.Upper+0.01, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Order.Lower, ymax=Order.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Order.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Order.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Order.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,3.7)+
  scale_color_manual(values=fun.col)+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )   

fungi.ord.tot<-ggarrange(fungi.ord.dens, fungi.ord, nrow=2, ncol=1, common.legend = FALSE, heights=c(1,0.4), align="v")

#Class info
fungi.class.dat<-subset(fungi, subset=(variable=="Class"))
fungi.class.dat$title<-"Class"
fungi.class.dens<-ggplot(fungi.class.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Class")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,3.7)+
  ylim(0,6)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=fun.col)+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

fungi.class<-ggplot(fungi.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Class.Lower-0.01, ymax=Class.Upper+0.01, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Class.Lower, ymax=Class.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Class.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Class.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Class.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,3.7)+
  scale_color_manual(values=fun.col)+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )   

fungi.class.tot<-ggarrange(fungi.class.dens, fungi.class, nrow=2, ncol=1, common.legend =FALSE, heights=c(1,0.4), align="v")


#phylum info
fungi.phy.dat<-subset(fungi, subset=(variable=="Phylum"))
fungi.phy.dat$title<-"Phylum"
fungi.phy.dens<-ggplot(fungi.phy.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  ggtitle("Fungi")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,3.7)+
  ylim(0,6)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values=fun.col)+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=18, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

fungi.phy<-ggplot(fungi.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Phylum.Lower-0.01, ymax=Phylum.Upper+0.01, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Phylum.Lower, ymax=Phylum.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Phylum.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Phylum.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Phylum.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,3.7)+
  scale_color_manual(values=fun.col)+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )  

fungi.phy.tot<-ggarrange(fungi.phy.dens, fungi.phy, nrow=2, ncol=1, common.legend = FALSE, heights=c(1,0.4), align="v")

fungi.dist<-ggarrange(fungi.phy.tot, fungi.class.tot, fungi.ord.tot, fungi.fam.tot, fungi.gen.tot,
          ncol=1,heights=c(0.85,0.75,0.75,0.75,1),
          common.legend = TRUE)





####Plant time####

#virid.unique.gs

virid.unique.gs$Phylum<-as.character(virid.unique.gs$Phylum)
virid.unique.gs$Phylum<-as.factor(virid.unique.gs$Phylum)

virid.line<-matrix(,nrow=length(levels(virid.unique.gs$Phylum)), ncol=27)
colnames(virid.line)<-c("Kingdom", "Phylum",
                        "Genus.Mean", "Genus.Lower", "Genus.Upper", "Genus.95L", "Genus.95U", 
                        "Family.Mean", "Family.Lower", "Family.Upper", "Family.95L", "Family.95U", 
                        "Order.Mean", "Order.Lower", "Order.Upper", "Order.95L", "Order.95U", 
                        "Class.Mean", "Class.Lower", "Class.Upper", "Class.95L", "Class.95U", 
                        "Phylum.Mean", "Phylum.Lower", "Phylum.Upper", "Phylum.95L", "Phylum.95U")
for(i in 1:length(levels(virid.unique.gs$Phylum))){
  virid.line[i,1]<-"Viridiplantae"
  virid.line[i,2]<-levels(virid.unique.gs$Phylum)[i]
  #genus info
  virid.line[i,3]<-mean(virid.unique.gs$Genus.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])
  virid.line[i,4]<-range(virid.unique.gs$Genus.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])[1]
  virid.line[i,5]<-range(virid.unique.gs$Genus.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])[2]
  output<-t.test(virid.unique.gs$Genus.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  virid.line[i,6]<-output$conf.int[1]
  virid.line[i,7]<-output$conf.int[2]
  #family info
  virid.line[i,8]<-mean(virid.unique.gs$Family.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])
  virid.line[i,9]<-range(virid.unique.gs$Family.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])[1]
  virid.line[i,10]<-range(virid.unique.gs$Family.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])[2]
  output<-t.test(virid.unique.gs$Family.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  virid.line[i,11]<-output$conf.int[1]
  virid.line[i,12]<-output$conf.int[2]
  #Order info
  virid.line[i,13]<-mean(virid.unique.gs$Order.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])
  virid.line[i,14]<-range(virid.unique.gs$Order.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])[1]
  virid.line[i,15]<-range(virid.unique.gs$Order.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])[2]
  output<-t.test(virid.unique.gs$Order.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  virid.line[i,16]<-output$conf.int[1]
  virid.line[i,17]<-output$conf.int[2]
  #class info
  virid.line[i,18]<-mean(virid.unique.gs$Class.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])
  virid.line[i,19]<-range(virid.unique.gs$Class.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])[1]
  virid.line[i,20]<-range(virid.unique.gs$Class.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])[2]
  output<-t.test(virid.unique.gs$Class.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  virid.line[i,21]<-output$conf.int[1]
  virid.line[i,22]<-output$conf.int[2]
  #phylum info
  virid.line[i,23]<-mean(virid.unique.gs$Phylum.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])
  virid.line[i,24]<-range(virid.unique.gs$Phylum.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])[1]
  virid.line[i,25]<-range(virid.unique.gs$Phylum.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]])[2]
  output<-t.test(virid.unique.gs$Phylum.Mean[virid.unique.gs$Phylum==levels(virid.unique.gs$Phylum)[i]], mu=1,alternative = "two.sided", var.equal = FALSE)
  virid.line[i,26]<-output$conf.int[1]
  virid.line[i,27]<-output$conf.int[2]
}
virid.line<-as.data.frame(virid.line)
virid.line <- virid.line %>%
  mutate(across(c("Kingdom", "Phylum"), as.factor))


virid.line <- virid.line %>%
  mutate(across(c("Genus.Mean", "Genus.Lower", "Genus.Upper", "Genus.95L", "Genus.95U", 
                  "Family.Mean", "Family.Lower", "Family.Upper", "Family.95L", "Family.95U", 
                  "Order.Mean", "Order.Lower", "Order.Upper", "Order.95L", "Order.95U", 
                  "Class.Mean", "Class.Lower", "Class.Upper", "Class.95L", "Class.95U", 
                  "Phylum.Mean", "Phylum.Lower", "Phylum.Upper", "Phylum.95L", "Phylum.95U"), as.numeric))


virid.gen<-ggplot(virid.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Genus.Lower-0.01, ymax=Genus.Upper+0.01, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Genus.Lower, ymax=Genus.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Genus.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Genus.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Genus.Mean), size=1.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,8.2)+
  scale_color_manual(values="#395D9CFF")+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=12, face="bold"))

virid.gen.dat<-subset(viridiplantae, subset=(variable=="Genus"))
virid.gen.dat$title<-"Genus"
virid.gen.dens<-ggplot(virid.gen.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Genus")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,8.2)+
  ylim(0,2)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values="#395D9CFF")+
  theme(legend.position="bottom",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=12, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),nrow = 2))

virid.gen.tot<-ggarrange(virid.gen.dens, virid.gen, nrow=2, ncol=1, common.legend = TRUE, legend="bottom", heights=c(1,0.4), align="v")

#family info
virid.fam.dat<-subset(viridiplantae, subset=(variable=="Family"))
virid.fam.dat$title<-"Family"
virid.fam.dens<-ggplot(virid.fam.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Family")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,8.2)+
  ylim(0,2)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values="#395D9CFF")+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

virid.fam<-ggplot(virid.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Family.Lower-0.01, ymax=Family.Upper+0.01, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Family.Lower, ymax=Family.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Family.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Family.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Family.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,8.2)+
  scale_color_manual(values="#395D9CFF")+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )  

virid.fam.tot<-ggarrange(virid.fam.dens, virid.fam, nrow=2, ncol=1, common.legend = FALSE, heights=c(1,0.4), align="v")

#Order info
virid.ord.dat<-subset(viridiplantae, subset=(variable=="Order"))
virid.ord.dat$title<-"Order"
virid.ord.dens<-ggplot(virid.ord.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Order")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,8.2)+
  ylim(0,2)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values="#395D9CFF")+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))


virid.ord<-ggplot(virid.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Order.Lower-0.01, ymax=Order.Upper+0.01, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Order.Lower, ymax=Order.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Order.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Order.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Order.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,8.2)+
  scale_color_manual(values="#395D9CFF")+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )   

virid.ord.tot<-ggarrange(virid.ord.dens, virid.ord, nrow=2, ncol=1, common.legend = FALSE, heights=c(1,0.4), align="v")

#Class info
virid.class.dat<-subset(viridiplantae, subset=(variable=="Class"))
virid.class.dat$title<-"Class"

virid.class.dens<-ggplot(virid.class.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  #ggtitle("Class")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,8.2)+
  ylim(0,2)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values="#395D9CFF")+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=16, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

virid.class<-ggplot(virid.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Class.Lower-0.01, ymax=Class.Upper+0.01, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Class.Lower, ymax=Class.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Class.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Class.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Class.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,8.2)+
  scale_color_manual(values="#395D9CFF")+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )   

virid.class.tot<-ggarrange(virid.class.dens, virid.class, nrow=2, ncol=1, common.legend =FALSE, heights=c(1,0.4), align="v")


#phylum info
virid.phy.dat<-subset(viridiplantae, subset=(variable=="Phylum"))
virid.phy.dat$title<-"Phylum"
virid.phy.dens<-ggplot(virid.phy.dat, aes(x=value, fill=Phylum))+
  geom_vline(xintercept=1, color="red", size=1, linetype="dashed")+
  ggtitle("Viridiplantae")+
  geom_density(alpha=0.7)+
  facet_wrap(~title, ncol=1)+
  xlim(0,8.2)+
  ylim(0,2)+
  #xlab("GS / Average GS for Taxonomic Rank")+
  theme_bw()+
  scale_fill_manual(values="#395D9CFF")+
  theme(legend.position="none",
        axis.title=element_blank(),
        #axis.title.x=element_text(size=8, face="bold", color="black"),
        legend.title=element_blank(),
        legend.text=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x=element_blank(),
        strip.text = element_text(hjust=0, size=16, face="bold"),
        plot.title = element_text(size=18, face="bold"))+
  guides(fill = guide_legend(override.aes = list(alpha=1),ncol = 2))

virid.phy<-ggplot(virid.line, aes(x=Phylum))+
  geom_linerange(aes(ymin=Phylum.Lower-0.01, ymax=Phylum.Upper+0.01, x=Phylum), color="black", size=2.75)+
  geom_linerange(aes(ymin=Phylum.Lower, ymax=Phylum.Upper, x=Phylum, color=Phylum), size=2)+
  geom_point(aes( y=Phylum.95L), size=2, shape=22, fill="grey")+
  geom_point(aes( y=Phylum.95U), size=2, shape=22, fill="grey")+
  geom_point(aes(y=Phylum.Mean), size=2.5, shape=24, fill="maroon")+
  coord_flip()+
  ylim(0,8.2)+
  scale_color_manual(values="#395D9CFF")+
  ylab("GS/Average GS for Taxonomic Rank")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12, color="black"),
        axis.text.y=element_blank(),
        axis.title=element_blank()#,
        #axis.title.x=element_text(size=12, face="bold")
  )  

virid.phy.tot<-ggarrange(virid.phy.dens, virid.phy, nrow=2, ncol=1, common.legend = FALSE, heights=c(1,0.4), align="v")

plant.dist<-ggarrange(virid.phy.tot, virid.class.tot, virid.ord.tot, virid.fam.tot, virid.gen.tot,
          ncol=1,heights=c(0.85,0.75,0.75,0.75,1),
          common.legend = TRUE)

#figure s5
ggarrange(meta.dist, plant.dist, fungi.dist, ncol=3)
