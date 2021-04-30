# Andes Lanzen 2021

## Initial treatment and filtering of 18S metbarcoding data based on taxonomic classification
## and distribution to remove potential non-target OTUs, contaminant OTUs, cross-contaminant reads, OTUs from
## pelagic organisms and finally rare OTUs that may bias alpha-diveristy and dissimilarity


# Change to correct path for local directory (github clone) and uncomment:
# setwd("/home/alanzen/projects/Metamon2/")

require(vegan)

source('R/utils/heterogeneity_rarefaction_functions.R')
source('R/utils/R/filtering.R')
source('R/utils/R/diversity.r')
source('R/utils/R/correlationTests.r')
source('R/utils/R/taxaplot.R')
source('R/utils/R/mergeOTUTable.R')

# ------ Read data --------

# Read metadata that is common for COI and 18S, limit to 18S samples and sort alphabetically
md.all = read.csv(file="WP1_metadata.csv",header=T,row.names=5)
mdR.18SsuR = md.all[md.all$X18S,]
mdR.18SsuR$SSU_ID <- gsub("\\-","\\.",mdR.18SsuR$SSU_ID)
mdR.18SsuR = mdR.18SsuR[order(mdR.18SsuR$SSU_ID),]

# Read unfiltered OTU table from SWARM
otus.all.18S = read.delim("SWARM_WP1_20200826_18S/CREST_LULU/SWARM_table_curated.tsv",
                     row.names=1,header=T,sep="\t")

# Make dataframe with taxonomy information
tax.18S=data.frame(row.names=row.names(otus.all.18S), classification = otus.all.18S$classification)
bestTx = array(dim=dim(tax.18S)[1])
for (i in 1:length(bestTx)) bestTx[i] <- tail(unlist(strsplit(as.character(tax.18S[i,]), split=";", 
                                                              fixed=TRUE)), 1)
tax.18S$bestTx = bestTx

# Sort and name OTU table in the same way as metadata table and control
otus.t.18S = as.data.frame(t(otus.all.18S[,-dim(otus.all.18S)[2]]))
otus.t.18S = otus.t.18S[order(row.names(otus.t.18S)),]
table(row.names(otus.t.18S) == mdR.18SsuR$SSU_ID)
row.names(otus.t.18S) = row.names(mdR.18SsuR)


# ----------- Filtering of potential contaminants ------------


dim(otus.t.18S) #we start with 33,994 OTUs

# Remove unclassified OTUs, pelagic taxa and obvious contaminants

for (taxonDel in c("No hits","Insecta","Mammalia","Aves","Myxini",
                   "Arachnida","Actinopterygii","Chondrichthyes","Collembola")) {
  toDel =  row.names(tax.18S)[grep(taxonDel,tax.18S$classification)]
  if (length(toDel)>0){
    otus.t.18S = otus.t.18S[,!(names(otus.t.18S) %in% toDel)]
    print(paste("Removing ",taxonDel))
    print(dim(otus.t.18S)[2])
    print(sum(otus.t.18S))
  }
}

# 33,479,561
# [1] "Removing  No hits"
# [1] 33089
# [1] 33443810
# [1] "Removing  Insecta"
# [1] 33084
# [1] 33443493
# [1] "Removing  Mammalia"
# [1] 33082
# [1] 33443361
# [1] "Removing  Arachnida"
# [1] 33068
# [1] 33392203
# [1] "Removing  Actinopterygii"
# [1] 33053
# [1] 33379493
# [1] "Removing  Chondrichthyes"
# [1] 33052
# [1] 33379489
# [1] "Removing  Collembola"
# [1] 33050 <- removes 944 OTUs, mostly No hits
# 33,379,482 <- removes ca 100k reads

## Filter probable cross-contaminant reads
otus.t.18S = filterCrossContaminants2(otus.t.18S,100)
sum(otus.t.18S) #33,374,039 (ca 5,400 reads less)

## Prevalence based filtering based on negative controls (decontam)
require(decontam)
om = as.matrix(otus.t.18S[mdR.18SsuR$Plate18S>1,])
mc = mdR.18SsuR[mdR.18SsuR$Plate18S>1,]
predicted_contaminants = isContaminant(om, neg=(mc$Type=="C"),batch=mc$Plate18S)
summary(predicted_contaminants) #->352 predicted contaminants
table(row.names(predicted_contaminants) == names(otus.t.18S)) # Yes
sc = (!is.na(predicted_contaminants$p) & predicted_contaminants$p<=.05) # 156
sigCont = predicted_contaminants[sc,]
sigCont$Tax = tax.18S$bestTx[row.names(tax.18S) %in% row.names(sigCont)]
summary(as.factor(sigCont$Tax))

#unclass. Alveoloata, Cercazoa, Enoplida and Eukaryota most common
# 

## Remove the suspected contaminants

cont = (names(otus.t.18S) %in% row.names(sigCont))
otus.t.18S = otus.t.18S[,-cont]
sum(otus.t.18S) #33,374,005 removed 5.5k reads

# Remove the blank and mock samples from further analysis
otusR.18S = otus.t.18S[mdR.18SsuR$Type == "S" | mdR.18SsuR$Type=="R",]
sum(otusR.18S) #32,915,206

otusR.18S = otusR.18S[,colSums(otusR.18S)>0] 
dim(otusR.18S) #Retains 33,032 OTUs / 297 stations
mdR.18S = droplevels(mdR.18SsuR[mdR.18SsuR$Type == "S" | mdR.18SsuR$Type=="R",])
table(row.names(otusR.18S) == row.names(mdR.18S))

# Remove OSF17.2.2019 w 3000 reads since others have >50k
mdR.18S = mdR.18S[rowSums(otusR.18S)>1E4,]
otusR.18S = otusR.18S[rowSums(otusR.18S)>1E4,]
dim(otusR.18S) #296 33,032 OTUs
sum(otusR.18S) #32911877


# ---- Plankton filtering -----------

## Read taxon filtering list and make lists of taxa to remove completely or only unclassified
otu_plankton_filtering = read.csv("SWARM_WP1_20200826_18S/OTU_Filtering.csv",
                                  header=T,row.names=1)
otu_plankton = row.names(otu_plankton_filtering)[otu_plankton_filtering$Remove] # 10
otus.18S.np = otusR.18S[,!(names(otusR.18S) %in% otu_plankton)]
sum(otus.18S.np) #25,437,233 (7.5M reads, 23% removed, from just 10 OTUs)

taxon_plankton_filtering = read.csv("SWARM_WP1_20200826_18S/taxon_filtering.csv",
                                    header=T,row.names=1)

## Remove further unclassified OTUs from indicated taxa
remove_unclass = row.names(taxon_plankton_filtering)[taxon_plankton_filtering$Remove_unclass]
for (r in remove_unclass){
  remove = row.names(tax.18S)[tax.18S$bestTx==r]
  otus.18S.np = otus.18S.np[,!(names(otus.18S.np) %in% remove)]
  print(paste("Removed",as.character(length(remove)),"OTUs classified only as",r,
              " leaving",as.character(sum(otus.18S.np)),"reads."))
}
sum(otus.18S.np) #21,386,341 (4.1M reads removed)

## Remove all OTUs from indicated taxa
remove_all = row.names(taxon_plankton_filtering)[taxon_plankton_filtering$Remove_all]
remove_all <- gsub(" \\(kingdom\\)","",remove_all)

for (r in remove_all){
  remove = row.names(tax.18S)[grep(r,tax.18S$classification)]
  otus.18S.np = otus.18S.np[,!(names(otus.18S.np) %in% remove)]
  print(paste("Removed",as.character(length(remove)),"OTUs classified as",r,
              " leaving",as.character(sum(otus.18S.np)),"reads."))
}

dim(otus.18S.np) #25,319 OTUs (we lose ~7500)
sum(otus.18S.np) #18,706,459 (2.6M more removed)

## Write plankton filtered OTU table
write.table(t(otus.18S.np), "SWARM_WP1_20200826_18S/SWARM_table_noPlankton.csv", 
            quote=F, sep="\t")

# ---- Make barchart and check abundant OTUs for contamination manually ----

otus.named = decostand(otus.18S.np,method="total")
names(otus.named) = paste(names(otus.named),tax.18S[names(otus.named),]$bestTx)

grouping_info<-data.frame(row.names=row.names(mdR.18S), mdR.18S$Station)
pdf("img/18S/SV_chart_noPlankton_18S.pdf",width=50, height = 10)
taxaplot(30,grouping_info,otus.named)
dev.off()
