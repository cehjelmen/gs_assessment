#run Feb 20-27, 2024 for eukaryotes
####load packages####
library(rentrez)
library(XML)
library(dplyr)
library(httr)

options(timeout = max(3000000, getOption("timeout")))
#load your api key for ncbi
#set_entrez_key("YOUR_KEY")
Sys.getenv("ENTREZ_KEY")
###assembly numbers###
####Functions for the loop and getting accessions####
####function to get assembly accessions
get_assembly_accessions <- function(taxonomy_id) {
  # Construct the Entrez query for the specified taxonomic group
  query <- sprintf("txid%s[Organism:exp]", taxonomy_id)
  
  # Use rentrez to search the Assembly database
  search_result <- entrez_search(db = "assembly", term = query, retmax = 50000, config = httr::config(timeout(30000)))  # Adjust retmax as needed
  
  # Extract assembly accessions from the search results
  assembly_accessions <- search_result$ids
  
  return(assembly_accessions)
}
###function to get taxonomy info
get_taxonomy_info <- function(taxid) {
  # Use rentrez to retrieve taxonomy information
  taxonomy_data <- entrez_fetch(db = "taxonomy", id = taxid, retmode = "xml", rettype ="XML", config = httr::config(timeout(30000)) )
  # Use XML package to parse the XML data
  #seeing if a pause will help the han
  #Sys.sleep(1)
  taxonomy_parsed <- xmlTreeParse(taxonomy_data, useInternalNodes = TRUE)
  
  # Extract relevant information
  taxonomy_info <- list(
    taxid = xpathSApply(taxonomy_parsed, "//TaxId", xmlValue),
    scientific_name = xpathSApply(taxonomy_parsed, "//ScientificName", xmlValue),
    rank = xpathSApply(taxonomy_parsed, "//Rank", xmlValue),
    lineage = xpathSApply(taxonomy_parsed, "//Lineage", xmlValue)
  )
  
  return(taxonomy_info)
}

#insert your taxomic ID for NCBI
taxonomic_group_id<-2759 #eukaryotes taxonomic group ID

####getting accessions and making large matrix####
assembly_accessions <- get_assembly_accessions(taxonomic_group_id)
all.genome.dat<-matrix(,nrow=length(assembly_accessions), ncol=28)
colnames(all.genome.dat)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
                        "Species", "Organism", "TaxID", "Assembly_Accession", "Assembly_Status",
                        "Accession_Date", "Assembly_Length_Total", "Assembly_Length_Ungapped",
                        "BUSCO_RefseqAnnotationRelease", "BUSCO_Lineage", "BUSCO_Cover",
                        "BUSCO_Complete","BUSCO_SingleCopy", "BUSCO_Duplicated", "BUSCO_Fragmented",
                        "BUSCO_Missing", "BUSCO_TotalCount", "ContigN50", "ScaffoldN50", "Sequencing_Method",
                        "Sequencing_Attributes", "SRA_Species")

####The loop doing all the work####

for(i in 41311:length(assembly_accessions)){
  #let's you know what you're working on 
  print(paste("Working on assembly:", i))  
  acc<-assembly_accessions[i]
  #gets summary info for your accession #
  esums2<-entrez_summary(db="assembly", id=acc, config = httr::config(timeout(30000)))
  #let's you know what species
  print(paste("Species:", esums2$speciesname))
  #get taxonomic information
  taxid <- esums2$taxid  # taxID
  taxonomy_info <- get_taxonomy_info(taxid)
  #getting assembly sizes
  assembly.meta<-esums2$meta
  #remove excess info
  assembly.meta2<-sub(" </Stats>.*$", "</Stats>", assembly.meta)
  #parse xml
  xml_doc <- xmlParse(assembly.meta2, useInternalNodes = TRUE)
  # Convert the XML tree to a list
  xml_list <- xmlToList(xml_doc)

  #getting sequencing method
  #look up SRA by project number from assembly
  assm.link<-entrez_link(dbfrom = "assembly", db="all", id=esums2$uid, config = httr::config(timeout(30000)))
  #assm.link$links
  if(length(assm.link$links)==0){
    xml_list2<-"Empty"
  }else if(any(names(assm.link$links)=="assembly_biosample")){
    testing<-entrez_summary(db="biosample", id=assm.link$links$assembly_biosample, config = httr::config(timeout(30000)))
    #get links from biosample
    biosamp.link<-entrez_link(dbfrom="biosample", db="all", id=assm.link$links$assembly_biosample, config = httr::config(timeout(30000)))
  }  else{
  xml_list2<-"Empty"
  }
  
  if(any(names(biosamp.link$links)=="biosample_sra")){
  #search SRA using link
    sra.search<-entrez_search(db="sra", term =biosamp.link$links$biosample_sra[1], config = httr::config(timeout(30000)))
    platform<-c()
    seq.info<-c()
    #length(sra.search$ids)
    for(k in 1:length(sra.search$ids)){
      sra.res<-entrez_summary(db="sra", id=sra.search$ids[k], config = httr::config(timeout(30000)))
      #sra.res$expxml
      sra.xml<-sra.res$expxml
      sra.xml.2<-sub("</Summary>.*$", "</Summary>", sra.xml)
      xml_doc.sra <- xmlParse(sra.xml.2, useInternalNodes = TRUE)
      # Convert the XML tree to a list
      xml_list.sra <- xmlToList(sra.xml.2)
      platform<-c(xml_list.sra$Platform$text, platform)
      seq.info<-c(xml_list.sra$Platform$.attrs[[1]], seq.info)
    }
    platform.info<-paste(platform, collapse = ",")
    seq.info.info<-paste(seq.info, collapse=",")
    xml_list2<-c("not empty","not empty", "not empty")
  }else if(testing$organism==esums2$speciesname){
    sra.search<-entrez_search(db="sra", term =testing$accession, config = httr::config(timeout(30000)))
    if(length(sra.search$ids)==0){
      platform.info<-""
      seq.info.info<-""
      xml_list2<-"Empty"
    }else {
    platform<-c()
    seq.info<-c()
    for(l in 1:length(sra.search$ids)){
      sra.res<-entrez_summary(db="sra", id=sra.search$ids[l], config = httr::config(timeout(30000)))
      #sra.res$expxml
      sra.xml<-sra.res$expxml
      sra.xml.2<-sub("</Summary>.*$", "</Summary>", sra.xml)
      xml_doc.sra <- xmlParse(sra.xml.2, useInternalNodes = TRUE)
      # Convert the XML tree to a list
      xml_list.sra <- xmlToList(sra.xml.2)
      platform<-c(xml_list.sra$Platform$text, platform)
      seq.info<-c(xml_list.sra$Platform$.attrs[[1]], seq.info)
    }
    platform.info<-paste(platform, collapse = ",")
    seq.info.info<-paste(seq.info, collapse=",")
    xml_list2<-c("not empty","not empty", "not empty")
  }
    }else{
    platform.info<-""
    seq.info.info<-""
    xml_list2<-"Empty"
  }
  
  
  #Kingdom
  if(length(taxonomy_info$scientific_name[which(taxonomy_info$rank=="kingdom")])>0){
  all.genome.dat[i, 1]<- taxonomy_info$scientific_name[which(taxonomy_info$rank=="kingdom")]
  }else{
    all.genome.dat[i,1]<-""
  }
  #Phylum
  if(length(taxonomy_info$scientific_name[which(taxonomy_info$rank=="phylum")])>0){
  all.genome.dat[i, 2]<-taxonomy_info$scientific_name[which(taxonomy_info$rank=="phylum")]
  }else{
    all.genome.dat[i,2]<-""
  }
  #Class
  if(length(taxonomy_info$scientific_name[which(taxonomy_info$rank=="class")])>0){
  all.genome.dat[i, 3]<-taxonomy_info$scientific_name[which(taxonomy_info$rank=="class")]
  }else{
    all.genome.dat[i,3]<-""
  }
  #Order
  if(length(taxonomy_info$scientific_name[which(taxonomy_info$rank=="order")])>0){
  all.genome.dat[i, 4]<-taxonomy_info$scientific_name[which(taxonomy_info$rank=="order")]
  }else{
    all.genome.dat[i,4]<-""
  }
  #Family
  if(length(taxonomy_info$scientific_name[which(taxonomy_info$rank=="family")])>0){
  all.genome.dat[i, 5]<-taxonomy_info$scientific_name[which(taxonomy_info$rank=="family")]
  }else{
    all.genome.dat[i,5]<-""
  }
  #Genus
  if(length(taxonomy_info$scientific_name[which(taxonomy_info$rank=="genus")])>0){
  all.genome.dat[i, 6]<-taxonomy_info$scientific_name[which(taxonomy_info$rank=="genus")]
  }else{
    all.genome.dat[i,6]<-""
  }
  #Species
  all.genome.dat[i, 7]<-esums2$speciesname
  #Organism
  all.genome.dat[i, 8]<-esums2$organism
  #TaxID
  all.genome.dat[i, 9]<-esums2$taxid
  #Assembly accession
  all.genome.dat[i, 10]<-esums2$assemblyaccession
  #Assembly Status
  all.genome.dat[i, 11]<-esums2$assemblystatus
  #Date of publication
  all.genome.dat[i, 12]<-testing$publicationdate
  #Assembly.length Total
  all.genome.dat[i, 13]<- xml_list[14]$Stat$text
  #Assembly Length ungapped
  all.genome.dat[i, 14]<-xml_list[15]$Stat$text
  #busco scores (all of them)
  ##refseqannotationrelease 15
  all.genome.dat[i,15]<-esums2$busco$refseqannotationrelease
  ##buscolineage 16
  all.genome.dat[i,16]<-esums2$busco$buscolineage
  ##buscover 17
  all.genome.dat[i,17]<-esums2$busco$buscover
  ##complete 18
  all.genome.dat[i,18]<-esums2$busco$complete
  ##singlecopy 19
  all.genome.dat[i,19]<-esums2$busco$singlecopy
  ##duplicated 20
  all.genome.dat[i,20]<-esums2$busco$duplicated
  ##fragmented 21
  all.genome.dat[i,21]<-esums2$busco$fragmented
  ##missing 22
  all.genome.dat[i,22]<-esums2$busco$missing
  ##totalcount 23
  all.genome.dat[i,23]<-esums2$busco$totalcount
  #contig n50 24
  all.genome.dat[i,24]<-esums2$contign50
  #scaffold n50 25
  all.genome.dat[i,25]<-esums2$scaffoldn50
  #sequencing method 26 
  if(xml_list2[1]=="Empty"){
    all.genome.dat[i,26]<-""
    all.genome.dat[i,27]<-""
  }else if(xml_list2[3]=="withdrawn"){
    all.genome.dat[i,26]<-""
    all.genome.dat[i,27]<-""
    } else if(xml_list2[3]=="suppressed"){
      all.genome.dat[i,26]<-""
      all.genome.dat[i,27]<-""
    }else if(xml_list2[3]=="new"){
      all.genome.dat[i,26]<-""
      all.genome.dat[i,27]<-""
    }else if(xml_list2[3]=="wait"){
      all.genome.dat[i,26]<-""
      all.genome.dat[i,27]<-""
    } else if(length(xml_list2)==2){
      all.genome.dat[i,26]<-""
      all.genome.dat[i,27]<-""
    }  else{
  all.genome.dat[i,26]<-platform.info
  #sequencing attributes
  all.genome.dat[i,27]<-seq.info.info
    }
  all.genome.dat[i,28]<-testing$organism
}

#write all data as csv
#write.csv(all.genome.dat, "alldata_feb_27_24.csv")

####Processing output####
#make data frame
all.genome.dat.df<-as.data.frame(all.genome.dat)
#make numeric things numeric
all.genome.dat.df <- all.genome.dat.df %>%
  mutate(across(c(Assembly_Length_Total, Assembly_Length_Ungapped, BUSCO_SingleCopy, BUSCO_Complete, BUSCO_Duplicated, BUSCO_Fragmented, BUSCO_Missing, BUSCO_TotalCount,
                  ContigN50, ScaffoldN50), as.numeric))
#make columns factors
all.genome.dat.df <- all.genome.dat.df %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Assembly_Status, SRA_Species), as.factor))
#kingdoms are "", Fungi, Metazoa, Viridiplantae
#56 phylums (including "")
#19 phylums (including "") when there is no kingdom listed

length(levels(all.genome.dat.df$SRA_Species))
#17374 species across 41358 assemblies

#separate out sequencing platform info
sequencing.info.matrix<-matrix(,nrow=nrow(all.genome.dat.df), ncol=3)
colnames(sequencing.info.matrix)<-c("Species", "Sequencing_method_1", "Sequencing_method_2")
#i<-771
for(i in 1:nrow(all.genome.dat.df)){
  sequencing.info.matrix[i,1]<-all.genome.dat.df$Species[i]
  seq.dat<-all.genome.dat.df$Sequencing_Method[i]
  seq.dat2<-unlist(strsplit(seq.dat, split=","))
  if(length(seq.dat2)==1){
    sequencing.info.matrix[i,2]<-seq.dat2
    sequencing.info.matrix[i,3]<-""
  }else if(length(unique(seq.dat2))==2){
    seq.dat2<-sort(seq.dat2, decreasing=TRUE)
    sequencing.info.matrix[i,2]<-unique(seq.dat2)[1]
    sequencing.info.matrix[i,3]<-unique(seq.dat2)[2]

  }else{
    seq.dat2<-sort(seq.dat2, decreasing=TRUE)
    sequencing.info.matrix[i,2]<-unique(seq.dat2)[1]
    sequencing.info.matrix[i,3]<-unique(seq.dat2)[2]
  }
}

#replace NAs with nothing
sequencing.info.matrix[is.na(sequencing.info.matrix)]<-""
length(unique(sequencing.info.matrix[,2]))
unique(sequencing.info.matrix[,2])
sequencing.info.matrix.df<-as.data.frame(sequencing.info.matrix)

all.genome.seq.info.df<-data.frame(all.genome.dat.df, sequencing.info.matrix.df)
#write csv with seq info
write.csv(all.genome.seq.info.df, "all_genome_seq_info_Feb_27.csv")

#41358 records for all data
seq.info<-subset(all.genome.seq.info.df, Sequencing_method_1=="ILLUMINA" | Sequencing_method_1=="PACBIO_SMRT" | Sequencing_method_1=="OXFORD_NANOPORE")

#27262 assemblies for just Illumina, ONT, or PacBio
#seq info narrowed to top three platforms
write.csv(seq.info, "all_genome_dat_PB_Ill_ONT.csv")

####Genome size data####
#genome size data was downloaded from genomesize.com on Feb 10 for all animals
#in order to format correctly, the first 5 rows were removed. It's a header about genome size.com
#GS data was downloaded for plants and fungus on Feb 20


