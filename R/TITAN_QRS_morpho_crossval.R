# Andes Lanzen 2021

## This script that was used for cross validation generation and evaluation of de novo biotic indices.
## Iterating over installations, we identify de novo indicators based on morphotaxonomy
## using TITAN2, then generate de novo biotic indices with quantile regression splines and
## apply them to the left out installation for evaluation

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
mdm.p$col_plot = as.numeric(mdm.p$Platform)

dim(morphOTUs.pra) #716 taxa, 93 samples
dim(mGenus.pra) #402

titan_threshold = .95 ## corresponds to summary of TITAN when reporting "pure and reliable"


for (t in trainOn) {
  
  print(paste("*** Training on",t,"***"))
  
  included_samples = !is.na(mdm.p[,t])
  
  md.trainAll = mdm.p[included_samples,]
  otus.trainAll = morphOTUs.pra[included_samples,]
  taxa.trainAll = mGenus.pra[included_samples,]
  
  platforms = unique(md.trainAll$Platform)
  preds = data.frame(row.names=row.names(md.trainAll), 
                     byOTUs=rep(NA,dim(md.trainAll)[1]),
                     byTaxa=rep(NA,dim(md.trainAll)[1]))
  
  ## Leave-one-out cross-validation - OTUs:
  ## Iterate over each platform, remove from training set and predict its value
  
  for (pf in platforms){
    
    print(pf)
    
    ## Use the other platforms to select indicator taxa w TITAN2
    md.train = md.trainAll[md.trainAll$Platform!=pf,]
    otus.train = otus.trainAll[md.trainAll$Platform!=pf,]
    otus.train.pa = decostand(otus.train, method="pa")
    otus.train = otus.train[,colSums(otus.train.pa)>5]
    otus.train = otus.train[,order(colSums(otus.train),decreasing=T)]
    
    otus.titan = titan(md.train[,t],otus.train,nBoot=100,numPerm=100) #ncpus=4,
    ot = as.data.frame(otus.titan$sppmax)
    
    ## Derive an AMBI-like index using regression splines with the picked indicators
    oindS = ot[ot$reliability>=titan_threshold & ot$purity>=titan_threshold,]
    otuEGs = splinePredict(oindS,t, otus.train, md.train)
    otu_indicators = otuEGs[otuEGs$value>0,,drop=F]
    
    ## Calculate values for this speicific index for the left out platform
    
    md.pred = md.trainAll[md.trainAll$Platform==pf,]
    otus.pred = otus.trainAll[md.trainAll$Platform==pf,]
    otu_w_Ind = otus.pred[,row.names(otu_indicators)]
    
    
    otu_BI_values = rep(0,dim(otus.pred)[1])
    for (i in c(1:dim(otu_indicators)[1])){
      w = getWeight(otu_indicators[i,],bi=t)
      otu_BI_values = otu_BI_values + w*otu_w_Ind[,i]
    }
    otu_BI_values = otu_BI_values/rowSums(otu_w_Ind)
    preds[md.trainAll$Platform==pf,"byOTUs"] = otu_BI_values
    
    # Write predictions to disk
    write.table(preds,paste("TITAN2_Splines/SplinePreds_morpho/",t,".tsv",sep=""),
                sep="\t",quote=F, col.names=NA)
  }
  
  
  
  ## Leave-one-out cross-validation - taxa:
  for (pf in platforms){
    
    print(pf)
    
    ## Use the other platforms to select indicator taxa w TITAN2
    md.train = md.trainAll[md.trainAll$Platform!=pf,]
    taxa.train = taxa.trainAll[md.trainAll$Platform!=pf,]
    taxa.train.pa = decostand(taxa.train, method="pa")
    taxa.train = taxa.train[,colSums(taxa.train.pa)>10]
    
    taxa.titan = titan(md.train[,t],taxa.train,nBoot=100,numPerm=100,memory=T) #ncpus=4
    tt = as.data.frame(taxa.titan$sppmax)
    
    ## Derive an AMBI-like index using regression splines with the picked indicators
    tindS = tt[tt$reliability>=titan_threshold & tt$purity>=titan_threshold,]
    taxaEGs = splinePredict(tindS,t, taxa.train, md.train)
    taxa_indicators = taxaEGs[taxaEGs$value>0,,drop=F]
    
    ## Calculate values for this speicific index for the left out platform
    md.pred = md.trainAll[md.trainAll$Platform==pf,]
    taxa.pred = taxa.trainAll[md.trainAll$Platform==pf,]
    taxa_w_Ind = taxa.pred[,row.names(taxa_indicators)]
    
    taxa_BI_values = rep(0,dim(taxa.pred)[1])
    for (i in c(1:dim(taxa_indicators)[1])){
      w = getWeight(taxa_indicators[i,],bi=t)
      taxa_BI_values = taxa_BI_values + w*taxa_w_Ind[,i]
    }
    taxa_BI_values = taxa_BI_values/rowSums(taxa_w_Ind)
    preds[md.trainAll$Platform==pf,"byTaxa"] = taxa_BI_values
    
  }
  
  
  # Write predictions to disk
  write.table(preds,paste("TITAN2_Splines/SplinePreds_morpho/",t,".tsv",sep=""), 
              sep="\t",quote=F, col.names=NA)
  
}


## ---- Draw from saved ----


for (t in trainOn) {
  
  preds = read.table(paste("TITAN2_Splines/SplinePreds_morpho/",t,".tsv",sep=""), 
                     sep="\t", header=F, skip=1, row.names=1)
  colnames(preds) = c("byOTUs","byTaxa")
  
  included_samples = !is.na(mdm.p[,t])
  
  md.trainAll = mdm.p[included_samples,]
  
  pdf(paste("img/Splines/",t,"_morpho_OTUs_crossval_pred",extra,".pdf",sep=""),width=8,height=5)
  plot_ml(data = preds$byOTUs[!is.na(preds$byOTUs)],
          metadata = md.trainAll[!is.na(preds$byOTUs),],
          xIndex = t, yIndex=t,
          aggreg = c("Station", "Platform"),
          title = "Spline predicted BI values")
  dev.off()
  
  pdf(paste("img/Splines/",t,"_morpho_Taxa_crossval_pred",extra,".pdf",sep=""),width=8,height=5)
  plot_ml(data = preds$byTaxa[!is.na(preds$byTaxa)], 
          metadata = md.trainAll[!is.na(preds$byTaxa),], 
          xIndex = t, yIndex=t,
          aggreg = c("Station", "Platform"),
          title = "Spline predicted BI values")
  dev.off()
}
