library(stringr)

#chromosome information was downloaded from
#https://coleoguy.github.io/karyotypes/
#on Feb. 20, 2024
#for diptera, mammals, polyneoptera, amphibians, coleoptera

#beetle data formatting
beetle<-read.csv("chromosomes/coleo-2024-02-20.csv")
beetle$Meioformula <- str_extract(beetle$Meioformula, "(?<=\\=)\\d+")
beetle$Meioformula<-as.numeric(beetle$Meioformula)
beetle<-na.omit(beetle)
(any(is.na(beetle$Meioformula)))
str(beetle)
colnames(beetle)[5]<-"HaploidNum"
colnames(beetle)[4]<-"Species"
beetle$Suborder<-NULL
beetle$Order<-NULL

#diptera data formatting
flies<-read.csv("chromosomes/diptera-2024-02-20.csv")
colnames(flies)[3]<-"HaploidNum"
flies[2:3]<-NULL
any(is.na(flies))

#amphibian data formatting
amphibians<-read.csv("chromosomes/amphib-2024-02-20.csv")
amphibians$HaploidNum<-amphibians$Diploid.number/2
amphibians[4:5]<-NULL
any(is.na(amphibians))
amphibians<-na.omit(amphibians)
amphibians$Order<-NULL

#polyneoptera data formatting
polyneoptera<-read.csv("chromosomes/polyneopteraDB-2024-02-20.csv")
polyneoptera$Genus<-NULL
polyneoptera$femalediploidnumber<-NULL
any(is.na(polyneoptera))
polyneoptera<-na.omit(polyneoptera)
colnames(polyneoptera)[4]<-"HaploidNum"
polyneoptera$Order<-NULL

#mammal data formatting
mammals<-read.csv("chromosomes/mammalia-2024-02-20.csv")
mammals$Genus<-NULL
colnames(mammals)[3]<-"Species"
mammals$HaploidNum<-mammals$Female2n/2
mammals$Female2n<-NULL
mammals<-na.omit(mammals)
any(is.na(mammals$HaploidNum))
mammals$Order<-NULL

#combine and save as csv
chromosome.dat<-rbind(amphibians, beetle, flies, polyneoptera, mammals)
write.csv(chromosome.dat, "chromosome.info.csv")
