# Andes Lanzen 2021

## This script that was used for global generation and analysis of de novo biotic indices:
## 1) identification of de novo indicators based on COI netabarcoding using TITAN2
## 2) generation of a a de novo biotic indices with quantile regression splines
## 3) analysis of taxonomic distribution across types of indicators (ecological groups)

# Change to correct path for local directory (github clone) and uncomment:
# setwd("/home/alanzen/projects/Metamon2/")

require(vegan)
require(TITAN2)
require(quantreg)
require(splines)
require(irr)
source("R/utils/ec_and_plot.R")
source("R/utils/splinePred.R")

trainOn = c("NSI","pi")

md.COI.p$col_plot = as.numeric(md.COI.p$Platform)

dim(otusR.COI.p.ra.ab) #2921 SVs, 97 samples
dim(asp) #780

TaxaToUse = 1000
titan_threshold = .95 ## corresponds to summary of TITAN when reporting "pure and reliable"

# ---- Analyse data using TITAN to find potential bioindicators ------

for (t in trainOn) {
  
  included_samples = !is.na(md.COI.p[,t])
  extra = "_np"
  
  md.train = md.COI.p[included_samples,]
  otus.train = otusR.COI.p.ra.ab[included_samples,]
  taxa.train = asp[included_samples,]
  
  # Much too slow to analyse all (min. 10% prevalence, then the 1000 most abundant)
  otus.train.pa = decostand(otus.train, method="pa")
  otus.train = otus.train[,colSums(otus.train.pa)>10]
  otus.train = otus.train[,order(colSums(otus.train), decreasing=T)]
  if (dim(otus.train)[2]>TaxaToUse) otus.train = otus.train[,c(1:TaxaToUse)]
  
  taxa.train.pa = decostand(taxa.train, method="pa")
  taxa.train = taxa.train[,colSums(taxa.train.pa)>10]
  
  # TITAN prediction
  otus.titan = titan(md.train[,t],otus.train,nBoot=100,numPerm=100)
  taxa.titan = titan(md.train[,t],taxa.train,nBoot=100,numPerm=100)
  write.csv(otus.titan$sppmax,paste("TITAN2_Splines/TITAN2_results_COI_np/OTUs_TITAN_all_",t,".csv", sep=""))
  write.csv(taxa.titan$sppmax,paste("TITAN2_Splines/TITAN2_results_COI_np/Taxa_TITAN_all_",t,".csv", sep=""))
}

# ----- Read TITAN results and analyse across training datasets ------

## Read output of TITAN, mark pure and reliable taxa

otuInd = matrix(ncol = length(trainOn), nrow = dim(otusR.COI.p.ra.ab)[2], 
                dimnames=list(names(otusR.COI.p.ra.ab),trainOn))

taxaInd = matrix(ncol = length(trainOn), nrow = dim(asp)[2],
                 dimnames=list(names(asp),trainOn))

otuGrp = matrix(ncol = length(trainOn), nrow = dim(otusR.COI.p.ra.ab)[2], 
                dimnames=list(names(otusR.COI.p.ra.ab),trainOn))

taxaGrp = matrix(ncol = length(trainOn), nrow = dim(asp)[2],
                 dimnames=list(names(asp),trainOn))

otuCP = matrix(ncol = length(trainOn), nrow = dim(otusR.COI.p.ra.ab)[2], 
               dimnames=list(names(otusR.COI.p.ra.ab),trainOn))

taxaCP = matrix(ncol = length(trainOn), nrow = dim(asp)[2],
                dimnames=list(names(asp),trainOn))

for (t in trainOn) {
  tind = read.csv(paste("TITAN2_Splines/TITAN2_results_COI_np/Taxa_TITAN_all_",t,".csv", sep=""),header=T,row.names=1)
  tindS = tind[tind$reliability>=titan_threshold & tind$purity>=titan_threshold & tind$obsiv.prob<.05,]
  for (ti in c(1:dim(tindS)[1])){
    taxaInd[row.names(tindS)[ti],t] = tindS$IndVal[ti]
    taxaGrp[row.names(tindS)[ti],t] = tindS$maxgrp[ti]
    taxaCP[row.names(tindS)[ti],t] = tindS$ienv.cp[ti]
  }
  
  oind = read.csv(paste("TITAN2_Splines/TITAN2_results_COI_np/OTUs_TITAN_all_",t,".csv", sep=""),header=T,row.names=1)
  oindS = oind[oind$reliability>=titan_threshold & oind$purity>=titan_threshold & oind$obsiv.prob<.05,]
  for (ti in c(1:dim(oindS)[1])){
    otuInd[row.names(oindS)[ti],t] = oindS$IndVal[ti]
    otuGrp[row.names(oindS)[ti],t] = oindS$maxgrp[ti]
    otuCP[row.names(oindS)[ti],t] = oindS$ienv.cp[ti]
  }
}


taxaInd[is.na(taxaInd)] = 0
otuInd[is.na(otuInd)] = 0

## ---- Spline prediction for taxa overview (and trivial, no cross-val performance eval) -----

## Read output of TITAN, mark pure and reliable taxa with p<0.05

# Groups for TITAN2 identified indicators, 
# reclassified as I -- V using quantile splines
otuGrp = matrix(ncol = length(trainOn), nrow = dim(otusR.COI.p.ra.ab)[2], 
                dimnames=list(names(otusR.COI.p.ra.ab),trainOn))

taxaGrp = matrix(ncol = length(trainOn), nrow = dim(asp)[2], 
                 dimnames=list(names(asp),trainOn))

for (t in trainOn) {
  print(paste("---",t,"---"))
  
  included_samples = !is.na(md.COI.p[,t])
  extra = "_np"
  
  md.train = md.COI.p[included_samples,]
  otus.train = otusR.COI.p.ra.ab[included_samples,]
  taxa.train = asp[included_samples,]
  
  ## Read TITAN2 output and select pure and reliable indicators
  oind = read.csv(paste("TITAN2_Splines/TITAN2_results_COI_np/OTUs_TITAN_all_",t,".csv", sep=""),header=T,row.names=1)
  oindS = oind[oind$reliability>=titan_threshold & oind$purity>=titan_threshold,]
  
  tind = read.csv(paste("TITAN2_Splines/TITAN2_results_COI_np/Taxa_TITAN_all_",t,".csv", sep=""),header=T,row.names=1)
  tindS = tind[tind$reliability>=titan_threshold & tind$purity>=titan_threshold,]
  
  ## Iterate over all selected OTU indicators
  otuEGs = splinePredict(oindS,t, otus.train, md.train, 
                         imageOutDir="img/Splines/OTUs_COI_modelled")
  otuGrp[,t] = otuEGs$value
  taxaEGs = splinePredict(tindS,t, taxa.train, md.train, 
                          imageOutDir="img/Splines/taxa_COI_modelled")
  taxaGrp[,t] = taxaEGs$value
  
  ## Calculate BI values based on TITAN2 picked spline models
  
  otu_indicators = otuEGs[otuEGs$value>0,,drop=F] 
  taxa_indicators = taxaEGs[taxaEGs$value>0,,drop=F]
  
  otu_w_Ind = otus.train[,row.names(otu_indicators)]
  taxa_w_Ind = taxa.train[,row.names(taxa_indicators)]
  
  otu_BI_values = rep(0,dim(otus.train)[1])
  for (i in c(1:dim(otu_indicators)[1])){
    w = getWeight(otu_indicators[i,],bi=t)
    otu_BI_values = otu_BI_values + w*otu_w_Ind[,i]
  }
  otu_BI_values = otu_BI_values/rowSums(otu_w_Ind)
  
  taxa_BI_values = rep(0,dim(taxa.train)[1])
  for (i in c(1:dim(taxa_indicators)[1])){
    w = getWeight(taxa_indicators[i,],bi=t)
    taxa_BI_values = taxa_BI_values + w*taxa_w_Ind[,i]
  }
  taxa_BI_values = taxa_BI_values/rowSums(taxa_w_Ind)
  
  
  opreds = otu_BI_values[!is.na(otu_BI_values)]
  mdo.test = md.train[!is.na(otu_BI_values),]
  
  tpreds = taxa_BI_values[!is.na(taxa_BI_values)]
  mdt.test = md.train[!is.na(taxa_BI_values),]
  
  
  pdf(paste("img/Splines/",t,"_COI_OTUs_all_no_crossval_pred",extra,".pdf",sep=""),
      width=8,height=5)
  plot_ml(data = opreds, metadata = mdo.test, 
          xIndex = t, yIndex=t,
          aggreg = c("Station", "Platform"), 
          title = "TITAN2+spline predicted BI values")
  dev.off()
  
  pdf(paste("img/Splines/",t,"_COI_taxa_all_no_crossval_pred",extra,".pdf",sep=""),
      width=8,height=5)
  plot_ml(data = tpreds, metadata = mdt.test, 
          xIndex = t, yIndex=t,
          aggreg = c("Station", "Platform"), 
          title = "TITAN2+spline predicted BI values")
  dev.off()
}


# ----------  Write summarised results to files --------
# 
# tx <- otus.all[row.names(otuInd),"classification"]
# otuInd$taxonomy = tx

## Write TITAN2 summarised results
write.csv(otuInd, "TITAN2_Splines/COI_np/OTUs_IndVal.csv")
write.csv(taxaInd, "TITAN2_Splines/COI_np/Taxa_IndVal.csv")
write.csv(otuGrp, "TITAN2_Splines/COI_np/OTUs_IndGroup.csv")
write.csv(taxaGrp, "TITAN2_Splines/COI_np/Taxa_IndGroup.csv")
write.csv(otuCP, "TITAN2_Splines/COI_np/OTUs_IndCP.csv")
write.csv(taxaCP, "TITAN2_Splines/COI_np/Taxa_IndCP.csv")


# ------- Make heatmap for taxon importance = IndVal -------

require(gplots)

table(row.names(taxaInd)==row.names(taxaGrp))
taxaInd[taxaGrp==0] <- 0
taxaIndNorm = decostand(t(taxaInd), method="total")
totalTaxaInd = colSums(taxaIndNorm)
taxaInd_chart = as.data.frame(taxaIndNorm[,order(totalTaxaInd,decreasing=TRUE)][,c(1:50)])
row.names(taxaInd_chart) = trainOn
rowSums(taxaInd_chart)

#taxaInd_chart = decostand(taxaInd_chart,method="total")

# Get ecogroup as labels for heatmap, convert to roman numerals, and, importantly set IndVal to 0
# where ecogroup is zero, i.e. not pure or reliable or disagrees with TITAN
taxaInd_labels = taxaGrp[names(taxaInd_chart),]
for (i in c(1:50)){
  taxaInd_labels[i,] = as.character(as.roman(taxaInd_labels[i,]))
  taxaInd_chart[is.na(taxaInd_labels[i,]),i] = NA
}
min(taxaInd_chart,na.rm=T) #.014
max(taxaInd_chart,na.rm=T) #.039


pal <- colorRampPalette(c("black", "blue", "purple","red","yellow"))(n = 256)

pdf(paste("img/COI/TITAN2/Taxon_importance_heatmap",extra,".pdf",sep=""),
    width=6,height=10)
# heatmap(as.matrix(taxaInd_chart),col=pal,Colv = NA,
#         Rowv = NA,revC=TRUE,scale="col",cexRow=.8,cexCol=.8,
#         margins = c(4,10))
heatmap.2(as.matrix(sqrt(t(taxaInd_chart))),col=pal,dendrogram="none",Colv = NA,
          Rowv = NA,revC=F,scale="none",cexRow=.7,cexCol=.8,
           margins = c(5,10),
          density.info="none",trace="none",
          cellnote = taxaInd_labels,
          notecol="white",key=T,notecex=.6,
          keysize=2.2,
          key.par=list(mar=c(6,1,12,1)),
          na.color = "black"
)

dev.off()


# --------- Make heatmap for OTU importance = IndVal --------

table(row.names(otuInd)==row.names(otuGrp))
otuInd[otuGrp==0] <- 0

otuIndNorm = decostand(t(otuInd), method="total")
totalotuInd = colSums(otuIndNorm)
otuInd = as.data.frame(otuInd)

otuInd_chart = as.data.frame(otuInd[order(totalotuInd,decreasing=TRUE),][c(1:50),])
otuInd_labels = otuGrp[row.names(otuInd_chart),]
for (i in c(1:50)){
  otuInd_labels[i,] = as.character(as.roman(otuInd_labels[i,]))
  otuInd_chart[i,is.na(otuInd_labels[i,])] = NA
}

titanOTUs = row.names(otuInd_chart)

bestTx <- tax.COI[row.names(otuInd_chart),"bestTx"]

row.names(otuInd_chart) = paste(gsub("SWARM_","",row.names(otuInd_chart)) ,bestTx)

pdf(paste("img/COI/TITAN2/OTU_importance_heatmap",extra,".pdf",sep=""),
    width=6.2,height=10)
heatmap.2(sqrt(as.matrix(otuInd_chart)),col=pal,dendrogram="none",Colv = NA,
          Rowv = NA,revC=F,scale="none",cexRow=.7,cexCol=.8,
          margins = c(5,10),
          density.info="none",trace="none",
          cellnote = otuInd_labels,
          notecol="white",key=T,notecex=.6,
          keysize=2.2,
          key.par=list(mar=c(6,1,12,1)),
          na.color="black"
)
dev.off()

# ----- Write output for KRONA ----

tg = as.data.frame(taxaGrp)
tg$taxonpath = ass.p.COI[row.names(tg),"Taxonpath"]

## Export KRONA text input ranging from min to max group for specific impact
writeKronaText = function(minGroup, maxGroup, impact, tg){
  
  groupImpact = tg[tg[,impact] >= minGroup & tg[,impact] <= maxGroup,]
  if (dim(groupImpact)[1]>0){
    groupImpact$print = "1"
    for (i in 1:dim(groupImpact)[1]){
      taxa = unlist(strsplit(as.character(groupImpact[i,]$taxonpath), split=";", 
                           fixed=TRUE))
      for (t in taxa){
        groupImpact[i,]$print = paste(groupImpact[i,]$print,t, sep="\t")
      }
    }
    write.table(as.data.frame(groupImpact$print), file=paste("img/Splines/KRONA/COI_group",
                                                             minGroup,"to",maxGroup,impact,".txt",
              sep=""),row.names = F,quote=F,col.names=F)
  }
}


for (impact in trainOn){
  writeKronaText(1,2,impact,tg)
  writeKronaText(3,3,impact,tg)
  writeKronaText(4,5,impact,tg)
}
