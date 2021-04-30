## Initial treatment and filtering of COI metbarcoding data based on taxonomic classification
## and distribution to remove potential non-target OTUs, contaminant OTUs, cross-contaminant reads, OTUs from
## pelagic organisms and finally rare OTUs that may bias alpha-diveristy and dissimilarity

# Change to correct path for local directory (github clone) and uncomment:
# setwd("/home/alanzen/projects/Metamon2/")

require(vegan)

source('R/utils/heterogeneity_rarefaction_functions.R')
source('R/utilsfiltering.R')
source('R/utilsdiversity.r')
source('R/utilscorrelationTests.r')
source('R/utilstaxaplot.R')
source('R/utilsmergeOTUTable.R')

# ------ Read data --------

# Read metadata that is common for COI and 18S, limit to COI samples and sort alphabetically
md.all = read.csv(file="WP1_metadata.csv",header=T,row.names=5)
md.COI = md.all[md.all$COI,]
md.COI$COI_ID <- gsub("\\-","\\.",md.COI$COI_ID)
md.COI = md.COI[order(md.COI$COI_ID),]

# Read unfiltered OTU table from SWARM
otus.all.COI = read.delim("SWARM_WP1_20200722_COI/CREST_LULU/SWARM_table_curated.tsv",
                     row.names=1,header=T,sep="\t")

# Make dataframe with taxonomy information
tax.COI=data.frame(row.names=row.names(otus.all.COI), classification = otus.all.COI$classification)
bestTx = array(dim=dim(tax.COI)[1])
for (i in 1:length(bestTx)) bestTx[i] <- tail(unlist(strsplit(as.character(tax.COI[i,]), split=";", 
                                                              fixed=TRUE)), 1)
tax.COI$bestTx = bestTx

# Sort and name OTU table in the same way as metadata table and control
otus.t.COI = as.data.frame(t(otus.all.COI[,-dim(otus.all.COI)[2]]))
otus.t.COI = otus.t.COI[order(row.names(otus.t.COI)),]
table(row.names(otus.t.COI) == md.COI$COI_ID)
row.names(otus.t.COI) = row.names(md.COI)

# ----------- Filtering of potential contaminants ------------

dim(otus.t.COI) #we start with 109,160 OTUs
sum(otus.t.COI) #..and 24,104,787 reads

# Remove pelagic taxa and obvious contaminants

for (taxonDel in c("Calanoida", "Insecta","Mammalia","Aves","Myxini",
                   "Arachnida","Actinopterygii","Chondrichthyes","Collembola")) {


  toDel = row.names(tax.COI)[grep(taxonDel,tax.COI$classification)]
  if (length(toDel)>0){
    otus.t.COI = otus.t.COI[,!(names(otus.t.COI) %in% toDel)]
    print(paste("Removing ",taxonDel))
    print(dim(otus.t.COI)[2])
    print(sum(otus.t.COI))
  }
}

# Limit OTUs to metazoa classified at phylum
met = row.names(tax.COI)[grep("Metazoa",tax.COI$classification)]
otus.COI.clean = otus.t.COI[,names(otus.t.COI) %in% met]

phylumOrBetter = rep(NA,dim(otus.COI.clean)[2])
for (i in c(1:length(phylumOrBetter))){
  tt = strsplit(as.character(tax.COI[names(otus.COI.clean)[i],]$classification),";")
  phylumOrBetter[i] = (length(tt[[1]]) > 4)
}
sum(phylumOrBetter) #8245 
otus.COI.clean = otus.COI.clean[,phylumOrBetter]

sum(otus.COI.clean) #4822266 reads
sum(otus.COI.clean)/sum(otus.t.COI) # 22% 

## FIlter probable cross-contaminant reads

otus.COI.clean = filterCrossContaminants2(otus.COI.clean,100)
sum(otus.COI.clean) #4821601


## Filter contaminant OTUs using decontam with prevalence based filtering based 
# on negative controls, and then do manual control

require(decontam)
om = as.matrix(otus.COI.clean[md.COI$PlateCOI>1,])
mc = md.COI[md.COI$PlateCOI>1,]
predicted_contaminants = isContaminant(om, neg=(mc$Type=="C"),batch=mc$PlateCOI)
summary(predicted_contaminants)

# Limit predicted contaminants to those with p<.05
sc = (!is.na(predicted_contaminants$p) & predicted_contaminants$p<=.05) 
sigCont = predicted_contaminants[sc,]
sigCont$Tax = tax.COI$bestTx[row.names(tax.COI) %in% row.names(sigCont)]
sigCont$Tax


contNames = row.names(sigCont)[c(4,10)]

otus.COI.clean = otus.COI.clean[,!(names(otus.COI.clean) %in% contNames)]
sum(otus.COI.clean) #4821586 reads remaining

## Remove the blank and mock samples from further analysis
otusR.COI = otus.COI.clean[md.COI$Type == "S" | md.COI$Type=="R",]
sum(otusR.COI) #4606254 reads remaining
otusR.COI = otusR.COI[,colSums(otusR.COI)>0] 
dim(otusR.COI) #Retains 8242 OTUs from 291 stations

mdR.COI = droplevels(md.COI[md.COI$Type == "S" | md.COI$Type=="R",])
table(row.names(otusR.COI) == row.names(mdR.COI))


# ---- Plankton filtering -----------

otus.COI.np = otusR.COI

## Read taxon filtering list and make lists of taxa to remove completely or only unclassified
taxon_plankton_filtering = read.csv("taxon_filtering.csv", header=T,row.names=1)

## Remove further unclassified OTUs from indicated taxa
remove_unclass = row.names(taxon_plankton_filtering)[taxon_plankton_filtering$Remove_unclass]
for (r in remove_unclass){
  remove = row.names(tax.COI)[tax.COI$bestTx==r]
  toRemove = (names(otus.COI.np) %in% remove)
  otus.COI.np = otus.COI.np[,!toRemove]
  if(sum(toRemove)>0){
    print(paste("Removed",as.character(sum(toRemove)),"OTUs classified only as",r,
              " leaving",as.character(sum(otus.COI.np)),"reads."))
  }
}
sum(otus.COI.np) 

## Remove all OTUs from indicated taxa
remove_all = row.names(taxon_plankton_filtering)[taxon_plankton_filtering$Remove_all]
remove_all <- gsub(" \\(kingdom\\)","",remove_all)

for (r in remove_all){
  remove = row.names(tax.COI)[grep(r,tax.COI$classification)]
  toRemove=(names(otus.COI.np) %in% remove)
  otus.COI.np = otus.COI.np[,!toRemove]
  if(sum(toRemove)>0){
    print(paste("Removed",as.character(sum(toRemove)),"OTUs classified as",r,
                " leaving",as.character(sum(otus.COI.np)),"reads."))
  }
}
# [1] "Removed 55 OTUs classified only as Cnidaria  leaving 4604487 reads."
# [1] "Removed 1 OTUs classified as Cyclopoida  leaving 4604475 reads."

dim(otus.COI.np) #8186 OTUs (we lose ~7500)
sum(otus.COI.np) #4604475

## Remove samples with a too low read depth below 1000 reads 
row.names(otus.COI.np)[rowSums(otus.COI.np)<1000]
#"RIN37.1.2018" "VFR32.1.2019" to be removed
mdR.COI = mdR.COI[rowSums(otus.COI.np)>1000,]
otus.COI.np = otus.COI.np[rowSums(otus.COI.np)>1000,]

## Write plankton filtered OTU table
write.table(t(otus.COI.np), "SWARM_WP1_20200722_COI/SWARM_table_noPlankton.csv", 
            quote=F, sep="\t")

## Write diversity statistics
writeDivStats("Diversity_filtered_COI_WP1_wReps_np.csv", otus.COI.np)

divR.COI = read.csv("Diversity_filtered_COI_WP1_wReps_np.csv",row.names=1)
divR.COI = divR.COI[,c("Reads","Richness","Rarefied.richness","H","J",
                       "Chao1")]

