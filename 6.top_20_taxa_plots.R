library(forcats)
library(ggplot2)
library(viridis)
library(ggpubr)
library(dunn.test)
dat.clean<-read.csv("all_GS_assembly.cleaned.csv")
alldat<-dat.clean

#Count occurrences of each unique species
string_counts <- table(alldat$Species)
#Sort the sort by occurrence
sorted_counts <- sort(string_counts, decreasing = TRUE)
#find the top 20 occurences
N <- 20  
top_strings <- names(sorted_counts)[1:N]
top_frequencies <- sorted_counts[1:N]

# Display the top 20 most common speices and their frequencies
# top_strings
# top_frequencies

top20<-matrix(,nrow=length(top_strings), ncol=7)
colnames(top20)<-c("Kingdom", "Phylum", "Genus", "Species", "GS", "Count", "DiffProp")
i<-1
for(i in 1:length(top_strings)){
  top20[i,1]<-unique(alldat$Kingdom[alldat$Species==top_strings[i]])
  top20[i,2]<-unique(alldat$Phylum[alldat$Species==top_strings[i]])
  top20[i,3]<-unique(alldat$Genus[alldat$Species==top_strings[i]])
  top20[i,4]<-top_strings[i]
  top20[i,5]<-mean(alldat$GS[alldat$Species==top_strings[i]])
  top20[i,6]<-top_frequencies[[i]]
  top20[i,7]<-mean(alldat$DiffProp[alldat$Species==top_strings[i]])
}

top20<-as.data.frame(top20)
top20 <- top20 %>%
  mutate(across(c("Kingdom", "Phylum", "Genus", "Species"), as.factor))
top20$GS<-as.numeric(top20$GS)
top20$DiffProp<-as.numeric(top20$DiffProp)
top20$Count<-as.integer(top20$Count)
#13 of top 20 are fungi
#totalling about 4500 records
#levels(top20$Kingdom)
top20$Kingdom<-factor(top20$Kingdom,levels=c("Fungi", "Viridiplantae", "Metazoa"))
#barplot of top 20 with occurence
top20.count<-ggplot(top20, aes(x=fct_reorder(Species,Count), y=Count, fill=Kingdom))+
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  scale_fill_viridis_d(option="G", begin=0.2, end=0.8)+
  theme_bw()+
  ylim(0,2050)+
  theme(axis.title=element_blank(),
        axis.text.y=element_text(face="italic", size=10, color="black"),
        axis.text.x=element_text(face="bold", size=10, color="black"),
        legend.title = element_blank(),
        legend.text=element_text(size=12, face="bold", color="black"),
        legend.position = "bottom")
#log gs for top 20
top20.gs<-ggplot(top20, aes(x=fct_reorder(Species, Count), y=log(GS), fill=Kingdom))+
  geom_point(shape=21,size=3)+
  coord_flip()+
  ylim(16,23)+
  scale_fill_viridis_d(option="G", begin=0.2, end=0.8)+
  ylab("log(GS (bp))")+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(face="bold", size=10, color="black"),
        #axis.text.y=element_text(face="italic", size=12, color="black"),
        axis.text.y=element_blank(),
        axis.text.x=element_text(face="bold", size=10),
        legend.title = element_blank(),
        legend.text=element_text(size=10, face="bold", color="black"),
        legend.position = "bottom")
#diff prop for top 20
diffprop.plot<-ggplot(top20, aes(x=fct_reorder(Species, Count), y=DiffProp, fill=Kingdom))+
  geom_point(shape=21,size=3)+
  coord_flip()+
  ylim(NA,0.8)+
  scale_fill_viridis_d(option="G", begin=0.2, end=0.8)+
  ylab("Mean Prop. Diff.")+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(face="bold", size=10, color="black"),
        #axis.text.y=element_text(face="italic", size=12, color="black"),
        axis.text.y=element_blank(),
        axis.text.x=element_text(face="bold", size=10),
        legend.title = element_blank(),
        legend.text=element_text(size=10, face="bold", color="black"),
        legend.position = "bottom")
#make one 3 panelled fig
count.info.plot<-ggarrange(top20.count, top20.gs,diffprop.plot, nrow=1, widths=c(1,0.3,0.3), common.legend = TRUE,
          legend = "none", align="h")
#regression of gs vs prop diff
gs.diff<-ggplot(top20, aes(x=log(GS), y=DiffProp))+
  geom_point(shape=21,size=4, aes(fill=Kingdom))+
  #oord_flip()+
  scale_fill_viridis_d(option="G", begin=0.2, end=0.8)+
  ylab("Proportional Difference")+
  theme_bw()+
  theme(axis.title.y=element_text(size=10, color="black", face="bold"),
        axis.title.x=element_text(face="bold", size=10, color="black"),
        #axis.text.y=element_text(face="italic", size=12, color="black"),
        axis.text.y=element_text(face="bold", size=10, color="black"),
        axis.text.x=element_text(face="bold", size=10),
        legend.title = element_blank(),
        legend.text=element_text(size=12, face="bold", color="black"),
        legend.position = "bottom")+
  geom_smooth(method="lm", se=FALSE, color="red", linetype="dashed")

#figure s2
ggarrange(top20.count, top20.gs,diffprop.plot,gs.diff,
          nrow=2,ncol=3, widths=c(1,0.3,0.3),align="h", 
          heights=c(1,0.5), labels=c("a","","","b"),
          common.legend = TRUE, legend="bottom", hjust=-0.1,
          vjust=1.1
          )

