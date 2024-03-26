library(ggplot2)
library(forcats)
library(dplyr)
library(ggpubr)
library(viridis)
library(reshape)
library(EnvStats)
library(rphylopic)

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
alldat_reordered$class <- cut(alldat_reordered$DiffProp, 
                              breaks = c(-Inf, -0.1, 0.1, Inf), 
                              labels = c("Larger", "Within 10%", "Smaller"))
viridis_colors <- viridis_pal(end=0.8, option="D")(3)
color_mapping<-setNames(viridis_colors, c("Smaller", "Within", "Larger"))
seq.colors<-c("#DF7C18", "#2D87A9", "#AD004F")

#getting images from phylopic
uuid <- get_uuid(name = "Saccharomyces cerevisiae")
homo.sapiens <- get_phylopic(uuid = "acf1cbec-f6ef-4a82-8ab5-be2f963b93f5")
yeast<-get_phylopic(uuid=uuid)
uuid <- get_uuid(name = "Arabidopsis thaliana")
arabidopsis<-get_phylopic(uuid=uuid)
uuid <- get_uuid(name = "Drosophila melanogaster")
drosophila<-get_phylopic(uuid=uuid)
uuid <- get_uuid(name = "Mus musculus")
mus<-get_phylopic(uuid=uuid)
uuid <- get_uuid(name = "Caenorhabditis elegans")
celegans<-get_phylopic(uuid=uuid)
uuid <- get_uuid(name = "Danio rerio")
zebrafish<-get_phylopic(uuid=uuid)
uuid <- get_uuid(name = "Zea mays")
corn<-get_phylopic(uuid=uuid)
uuid<-get_uuid(name="Pinus sylvestris")
pine<-get_phylopic(uuid=uuid)
uuid<-get_uuid(name="Amanita muscaria")
mushroom<-get_phylopic(uuid=uuid)
uuid<-get_uuid(name="Hydra vulgaris")
cnidaria<-get_phylopic(uuid=uuid)

ggplot(data=alldat_reordered, aes(x=fct_inorder(Species), y=DiffProp))+
  annotate('rect', xmin=-Inf, xmax=Inf, ymin=-0.1, ymax=0.1, alpha=.2, fill='black')+
  geom_point(aes(fill=class), shape=21,size=2, alpha=0.6)+
  coord_flip()+
  theme_bw()+
  #ggtitle(paste("Proportional difference between est. GS and assembly size"))+
  scale_fill_viridis_d(end=0.8, direction=-1)+theme(legend.position = "none")+
  geom_hline(aes(yintercept=0), color="red", linewidth=1.5, linetype="dashed")+
  geom_vline(aes(xintercept=1169), linewidth=0.75, color="black", linetype="dashed")+
  geom_vline(aes(xintercept=(1169+660)), linewidth=0.75, color="black", linetype="dashed")+
  geom_vline(aes(xintercept=(1169+660+480)), linewidth=0.75, color="black", linetype="dashed")+
  # geom_vline(aes(xintercept=(493.8+262.5)), linewidth=1, color="black")+
  # geom_vline(aes(xintercept=(493.8+262.5+753.0)), linewidth=1, color="black")+
  #geom_vline(aes(xintercept="Homo sapiens"), linewidth=0.75, color="red", linetype="dashed")+
  #geom_vline(aes(xintercept="Saccharomyces cerevisiae"), linewidth=0.75, color="red", linetype="dashed")+
  annotate("text", x=1100, y =-1.25, label="Metazoa", 
           fontface="bold", size=6, alpha=0.6,
           hjust=0)+
  annotate("text", x=1750, y =-1.25, label="Viridiplantae", 
           fontface="bold", size=6, alpha=0.6,
           hjust=0)+
  annotate("text", x=2250, y =-1.25, label="Fungi", 
           fontface="bold", size=6, alpha=0.6,
           hjust=0)+
  ylim(c(-1.25,1))+
  ylab("Proportional Difference from GS")+
  theme(legend.position="none",
        axis.text.y=element_blank(), 
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=13, colour="black", face="bold"),
        axis.text.x=element_text(size=12, colour="black"),
        panel.grid.major = element_blank(),
        axis.ticks.y = element_blank())+
  add_phylopic(img=homo.sapiens, x="Homo sapiens", y=-1.05, ysize=175, alpha=0.6)+
  add_phylopic(img=yeast, x="Saccharomyces cerevisiae", y=-1.2, ysize=55, alpha=0.6)+
  add_phylopic(img=arabidopsis, x="Arabidopsis thaliana", y=-0.45, ysize=125, alpha=0.6)+
  add_phylopic(img=drosophila, x="Drosophila melanogaster", y=-0.6, ysize=70, alpha=0.6)+
  add_phylopic(img=mus, x="Mus musculus", y=-0.7, ysize=100, alpha=0.6)+
  add_phylopic(img=celegans, x="Caenorhabditis elegans", y=-0.85, ysize=150, alpha=0.6)+
  add_phylopic(img=zebrafish, x="Danio rerio", y=-0.5, ysize=65, alpha=0.6)+
  add_phylopic(img=corn, x="Zea mays", y=-1.1, ysize=80, alpha=0.6)+
  add_phylopic(img=pine, x="Pinus sylvestris", y=-0.55, ysize=150, alpha=0.6)+
  add_phylopic(img=mushroom, x="Amanita muscaria", y=-.55, ysize=100, alpha=0.6)

