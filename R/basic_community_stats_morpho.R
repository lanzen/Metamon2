setwd("/home/alanzen/projects/Metamon/WP1")

require(vegan)
library(plyr)
source('R/utils/heterogeneity_rarefaction_functions.R')
source('~/kode/R/filtering.R')
source('~/kode/R/diversity.r')
source('~/kode/R/correlationTests.r')
source('~/kode/R/taxaplot.R')
source('~/kode/R/mergeOTUTable.R')
source('R/utils/calculatePIs.R')
source('R/utils/ec_and_plot.R')

# ------- Read and transform data --------

# Metadata
md.all = read.csv(file="WP1_metadata.csv",header=T,row.names=5)

# Read data from morpho
morpho.all = read.csv("Morphodata/Morphodata_BOLD.csv",row.names=1,header=T)

# Remove all undefiend and convert
morpho.all = morpho.all[,!is.na(morpho.all[1,])]

# Sort according to metadata
table(names(morpho.all)[-c(1:6)] %in% row.names(md.all)) #All
table(row.names(md.all) %in% names(morpho.all)) #All but 44
mdm = droplevels(md.all[names(morpho.all)[-c(1:6)],])
mdm$StationRep = paste(mdm$Station,mdm$Type,sep="")

# OTU / "species" df
morphOTUs = as.data.frame(t(morpho.all[,-c(1:6)]))
dim(morphOTUs) # 272 samples (of 296), 715 "morpho OTUs"

table(row.names(morphOTUs) == row.names(mdm))

# Merge metadata by stations
mdm.p = mergeMetadata(mdm, by="StationRep")
write.csv(mdm.p,"WP1_Metadata_pooled_morpho.csv")
mdm.p = read.csv("WP1_Metadata_pooled_morpho.csv",row.names=1,header=T)

# Fix parameters
mdm.p$TotalPAH = mdm.p$PAH*1000 #Conversion to ppb
mdm.p$TotalHC = mdm.p$THC
mdm.p$TotalPCB = NA
mdm.p$TotalDDT = NA
mdm.p$TotalHCH = NA
mdm.p$ADE = NA
mdm.p$OM = mdm.p$TOC
mdm.p = calculatePIs(mdm.p)

mdm.p$Northeast = mdm.p$WGS84E + mdm.p$WGS84N

summary(mdm.p)


# Merge "species" data by station
morphOTUs.p = mergeOTUTable(morphOTUs, mdm, by="StationRep")
morphOTUs.pra = decostand(morphOTUs.p, method="total")

# Check sanity
table(row.names(morphOTUs.p)==row.names(mdm.p))

# ------ Sum all rows with same taxon and return a transponded table like OTU table -------
sumTaxa = function(df){
  summed = data.frame(row.names = names(df)[-1])
  taxa = unique(df[,1])
  taxa = taxa[taxa!=""]
  for (taxon in taxa){
    summed[[taxon]] = colSums(df[df[,1]==taxon,-1])
  }
  return(summed)
}

# Sort out higher taxonomic ranks, summarising those with same name
mPhylum = sumTaxa(morpho.all[,-c(2:6)])
mClass = sumTaxa(morpho.all[,-c(1,3:6)])
mOrder = sumTaxa(morpho.all[,-c(1:2,4:6)])
mFam = sumTaxa(morpho.all[,-c(1:3,5:6)])
mGenus = sumTaxa(morpho.all[,-c(1:4,6)])

mPhylum.p = mergeOTUTable(mPhylum, mdm, by="StationRep")
mClass.p = mergeOTUTable(mClass, mdm, by="StationRep")
mOrder.p = mergeOTUTable(mOrder, mdm, by="StationRep")
mFam.p = mergeOTUTable(mFam, mdm, by="StationRep")
mGenus.p = mergeOTUTable(mGenus, mdm, by="StationRep")

table(rowSums(mPhylum.p) == rowSums(morphOTUs.p))
totalClassified = sum(mPhylum.p)
sum(mClass.p)/totalClassified # 98%
sum(mOrder.p)/totalClassified # 97%
sum(mFam.p)/totalClassified # 83%
sum(mGenus.p)/totalClassified # 75%

mPhylum.pra = mPhylum.p/rowSums(morphOTUs.p)
mClass.pra = mClass.p/rowSums(morphOTUs.p)
mOrder.pra = mOrder.p/rowSums(morphOTUs.p)
mFam.pra = mFam.p/rowSums(morphOTUs.p)
mGenus.pra = mGenus.p/rowSums(morphOTUs.p)

## ------ NMDS of individual replicates - total counts ----------

crp = colorRampPalette(c("black","red","orange","darkgreen","blue","purple"))
colors = crp(13)

nmds = metaMDS(morphOTUs)
ordiplot(nmds,type="none")#,xlim=range(-4,2),ylim=range(-1,1))
points(nmds,pch=as.numeric(mdm$Platform),col=colors[as.numeric(mdm$Platform)], cex=.8)
legend("bottomleft",pch=as.numeric(sort(unique(mdm$Platform))),
       col=colors[as.numeric(sort(unique(mdm$Platform)))],cex=.6,
       legend=sort(unique(mdm$Platform)), ncol=6)
ordihull(nmds,mdm$Region,label=T,col="black", kind="sd",draw="line")

## ------ NMDS of individual replicates - relative abundance counts -------

morphOTUs.h = decostand(morphOTUs,method="hellinger")
morphOTUs.ra = decostand(morphOTUs,method="total")
nmds.h = metaMDS(morphOTUs.h)
ordiplot(nmds.h,type="none")#,xlim=range(-4,2),ylim=range(-1,1))
points(nmds.h,pch=as.numeric(mdm$Platform),col=colors[as.numeric(mdm$Platform)], cex=.8)
legend("bottomleft",pch=as.numeric(sort(unique(mdm$Platform))),
       col=colors[as.numeric(sort(unique(mdm$Platform)))],cex=.6,
       legend=sort(unique(mdm$Platform)), ncol=6)
ordihull(nmds.h,mdm$Region,label=T,col="black", kind="sd",draw="line")

##  ----------- NMDS of pooled - w envfit, similar to metabarcoding  ----------- 


#13dp:

mp1 = mdm.p[,c("NPD","TotalPAH","Acenaften","Acenaftylen" ,"Antracen","Benz.a.antracen",
                  "Benzo.a.pyren","Benzo.bjk.fluoranten","Benzo.ghi.perylen",
                  "C1.alkyldibenzotiofener","C1.alkylfenantrener.antracener",
                  "C1.alkylnaftalener","C2.alkyldibenzotiofener","C2.alkylfenantrener.antracener",
                  "C3.alkyldibenzotiofener","C3.alkylfenantrener.antracener",
                  "C3.alkylnaftalener","Chrysen","Dibenz.ah.antracen","Dibenzotiofen",
                  "Fenantren","Fluoranten","Fluoren","Indeno.1.2.3.c.d.pyren",
                  "Naftalen","Pyren","PAH_PI", "THC_PI")]
# 54 missing

mp2 = mdm.p[,c("As","Ba","Cr","Cd","Cu","Hg","Pb","Zn","piMetals","pi", "piHC","THC","TOC",
                  "TotalHC","WGS84E","WGS84N","Northeast","Depth","NSI","Sand","Grus","Pelite",
               "Kornstorrelse")] #0 missing


mp4 = mdm.p[,c("Kurtosis","Pelite","Skewness","Sorting")] #1 missing

mpdiv = mdm.p[,c("MorphoS","MorphoH","ES100","Individtetthet",
                    "ISI","NSI","NQI1")] #0 missing

mp5 = mdm.p[,c("Extract.conc","Distance")] #21 missing


## NMDS, absolute counts
nmds.p = metaMDS(morphOTUs.p)
ordiplot(nmds.p,type="none")#,xlim=range(-3,5),ylim=range(-1,1))
points(nmds.p,pch=as.numeric(mdm.p$Platform),col=colors[as.numeric(mdm.p$Platform)],
       cex=.8)
legend("bottomright",pch=as.numeric(sort(unique(mdm.p$Platform))),
       col=colors[as.numeric(sort(unique(mdm.p$Platform)))],cex=.6, 
       legend=sort(unique(mdm.p$Platform)), ncol=4)
#text(nmds.p,pos=1,cex=.3,col=colors[as.numeric(mdm.p$Platform)])
ordihull(nmds.p,mdm.p$Region,label=T,col="black", kind="sd",draw="line")

efp1=envfit(nmds.p,mp1,na.rm=T)
# NPD                             0.98793  0.15489 0.6170  0.001 ***
# TotalPAH                        0.96774 -0.25196 0.2806  0.005 ** 
#   Acenaften                       0.99664 -0.08188 0.6410  0.001 ***
#   Acenaftylen                     0.95722 -0.28936 0.6552  0.001 ***
#   C1.alkylnaftalener              0.97107 -0.23880 0.6790  0.001 ***
#   C3.alkyldibenzotiofener         0.86049  0.50947 0.5228  0.001 ***
#   C3.alkylnaftalener              0.97605 -0.21754 0.6702  0.001 ***
#   Naftalen                        0.97895  0.20412 0.6105  0.001 ***
# PAH_PI                          0.92077  0.39010 0.3687  0.001 ***
#   THC_PI                          0.89389  0.44829 0.6173  0.001 ***
# (Slightly weaker when using relative abundance)  
  

efp2 = envfit(nmds.p,mp2)
# As        0.65143 -0.75870 0.1897  0.001 ***
#   Ba        0.62593 -0.77988 0.4472  0.001 ***
#   piMetals  0.34434 -0.93885 0.2722  0.001 ***
#   pi        0.76839 -0.63998 0.3281  0.001 ***
#   piHC      0.99870 -0.05088 0.3875  0.001 ***
#   THC       0.98604 -0.16654 0.2811  0.001 ***
#   TotalHC   0.98604 -0.16654 0.2811  0.001 ***
#   WGS84E   -0.67023 -0.74215 0.6670  0.001 ***
#   WGS84N   -0.56686 -0.82382 0.7109  0.001 ***
# Northeast     -0.70826 -0.70595 0.6362  0.001 *** 
#   Depth    -0.75320 -0.65780 0.6744  0.001 ***
# NSI           -0.99811  0.06146 0.7297  0.001 ***
#   Sand           0.07357 -0.99729 0.6182  0.001 ***
#   Grus          -0.27490 -0.96147 0.2017  0.001 ***
#   Pelite        -0.32493 -0.94574 0.3030  0.001 ***
#   Kornstorrelse -0.91092  0.41258 0.1407  0.001 ***
# (Comparable for relative abundance)

efp4 = envfit(nmds.p,mp4,na.rm=T)
# Kurtosis      -0.12409  0.99227 0.1825  0.001 ***
#   Pelite        -0.36996 -0.92905 0.4712  0.001 ***
# Skewness  0.10167 -0.99482 0.0983  0.011 *  
#   Sorting  -0.99712 -0.07578 0.0555  0.075 .  
#   (slightly weaker for relative abundance)

efp5 = envfit(nmds.p,mp5,na.rm=T)
# (extract conc. is not signficiant so this is effect of metabarcoding as expected)

efp6 = envfit(nmds.p,mpdiv,na.rm=T)
# Kind of obvious. Individtetthet not significant
# ---


plot(efp1,p.max=.001,cex=.5,col="purple")
plot(efp2,p.max=.001,cex=.5,col="purple")
plot(efp4,p.max=.001,cex=.5,col="purple")
# plot(efp5,p.max=.001,cex=.5,col="purple")
plot(efp6,p.max=.001,cex=.5,col="green")


# ---------- BBI calculated BI values ---------

library(BBI)

mt = as.data.frame(t(morphOTUs.pra))
mt = data.frame(taxon=row.names(mt), mt)
bbi.m = BBI(mt)
bbi.m$BBI
table(row.names(bbi.m$BBI) == row.names(mdm.p))
plot(bbi.m$BBI[,3],mdm.p$NSI,xlab="BBI calculated NSI",ylab="Reported NSI")
abline(a=0,b=1,col="grey")
#text(bbi.m$BBI[,3],mdm.p$NSI,pos=1,cex=.7,labels=row.names(mdm.p))

plot(bbi.m$BBI[,1],mdm.p$NSI)

mdm.p$col_plot = as.numeric(mdm.p$Platform)

pdf("img/morpho/NSI_BBI_v_reported.pdf",width=8,height=5)
plot_ml(data = bbi.m$BBI[,3], metadata = mdm.p, 
        xIndex = "NSI", yIndex = "NSI", aggreg = c("Station","Platform"),
        title = "BBI calculated NSI v reported")
dev.off()


pdf("img/morpho/ISI_BBI_v_reported.pdf",width=8,height=5)
plot_ml(data = bbi.m$BBI[,2], metadata = mdm.p, 
        xIndex = "ISI", yIndex = "ISI", aggreg = c("Station","Platform"),
        title = "BBI calculated ISI v reported")
dev.off()

pdf("img/morpho/AMBI_v_piHC.pdf",width=8,height=5)
plot_ml(data = bbi.m$BBI[,1], metadata = mdm.p, 
        xIndex = "piHC", yIndex = "AMBI", aggreg = c("Station","Platform"),
        title = "COI predicted AMBI v HC PI")
dev.off()


# ---- Mantel tests to compare to 18S and COI OTUs -----

## Make comparable datasets of metabarcoding OTUs
# The missing sample is removed also from 18S to make a fair comparison

md.comp = mdm.p[row.names(mdm.p) %in% row.names(otusR.COI.p.ra.ab),] # 1 missing
morpho.coiComp = morphOTUs.p[row.names(md.comp),]
otus.COI.morpho = otusR.COI.p.ra.ab[row.names(morpho.coiComp),]
otus.18S.morpho = otus.18S.np.p.ra.ab[row.names(morpho.coiComp),]

pc.comp = md.comp[,c("Depth","piHC","piMetals","Northeast","Sand","Grus","Pelite",
                     "Kornstorrelse")]
pc.comp = decostand(pc.comp,method="standardize")
pc.dist = dist(pc.comp)

pi.comp = mdm.p[row.names(morpho.coiComp),c("piHC","piMetals")]
pi.comp = decostand(pi.comp,method="standardize")
pi.dist = dist(pi.comp)


## Calculate Bray-Curtis dissimilarities

otus.18S.morpho.bc = vegdist(otus.18S.morpho)
otus.COI.morpho.bc = vegdist(otus.COI.morpho)

morpho.abs.comp.bc = vegdist(morpho.coiComp) # Absolute sp. counts
morpho.ra.comp.bc = vegdist(decostand(morpho.coiComp,method="total")) # Relative abundances

## Mantel tests

# 18S v morpho
mantel(otus.18S.morpho.bc, morpho.abs.comp.bc) # R=.83***
mantel(otus.18S.morpho.bc, morpho.ra.comp.bc) # R=.83*** (essentially equal)

# COI v morpho
mantel(otus.COI.morpho.bc, morpho.abs.comp.bc) # R=.66***
mantel(otus.COI.morpho.bc, morpho.ra.comp.bc) # R=.67*** (essentially equal)

# 18S v COI
mantel(otus.COI.morpho.bc, otus.18S.morpho.bc) #R=.72***

# Physchem v morpho, 18S and COI
mantel(pc.dist, morpho.abs.comp.bc) #R=.69***
mantel(pc.dist, otus.18S.morpho.bc) #R=.75***
mantel(pc.dist, otus.COI.morpho.bc) #R=.57***

# Impact v morpho, 18S and COI
mantel(pi.dist, morpho.abs.comp.bc) #R=.32***
mantel(pi.dist, otus.18S.morpho.bc) #R=.26***
mantel(pi.dist, otus.COI.morpho.bc) #R=.20***

## ----- Taxonomy plots ------

gi=data.frame(row.names=row.names(mdm.p), mdm.p$Platform)

pdf("img/morpho/Assignments.pdf",height=8,width=15)
taxaplot(25,gi,morphOTUs.p)
dev.off()

pdf("img/morpho/Assignments_RA.pdf",height=8,width=15)
taxaplot(25,gi,morphOTUs.pra)
dev.off()

pdf("img/morpho/phyla.pdf",height=8,width=15)
taxaplot(10,gi,mPhylum.p)
dev.off()

pdf("img/morpho/phyla_RA.pdf",height=8,width=15)
taxaplot(10,gi,mPhylum.pra)
dev.off()

pdf("img/morpho/classes_RA.pdf",height=8,width=15)
taxaplot(14,gi,mClass.pra)
dev.off()

pdf("img/morpho/orders_RA.pdf",height=8,width=15)
taxaplot(25,gi,mOrder.pra)
dev.off()

pdf("img/morpho/fam_RA.pdf",height=8,width=15)
taxaplot(25,gi,mFam.pra)
dev.off()

pdf("img/morpho/genera_RA.pdf",height=8,width=15)
taxaplot(30,gi,mGenus.pra)
dev.off()

# ----- Compare taxonomy to 18S and to COI  -----

taxa.18S.comp = taxa.p.18S[,c("Rank","Taxonpath",row.names(morpho.coiComp))]
met.18S = taxa.18S.comp["Metazoa (Animalia)",-c(1,2)]
met.18S.comp = taxa.18S.comp[grep("Metazoa",taxa.18S.comp$Taxonpath),]
names(met.18S.comp) = paste(names(met.18S.comp),"18S",sep=".")


taxa.COI.comp = taxa.p.COI[,c("Rank","Taxonpath",row.names(morpho.coiComp))]
names(taxa.COI.comp) = paste(names(taxa.COI.comp),"COI",sep=".")


mPhylum.comp = mPhylum.pra[row.names(morpho.coiComp),]
mClass.comp = mClass.pra[row.names(morpho.coiComp),]
mOrder.comp = mOrder.pra[row.names(morpho.coiComp),]
mFam.comp = mFam.pra[row.names(morpho.coiComp),]
mGenus.comp = mGenus.pra[row.names(morpho.coiComp),]

phyla.18S = as.data.frame(t(met.18S.comp[met.18S.comp$Rank=="phylum",-c(1,2)])) / t(met.18S)
class.18S = as.data.frame(t(met.18S.comp[met.18S.comp$Rank=="class",-c(1,2)])) / t(met.18S)
order.18S = as.data.frame(t(met.18S.comp[met.18S.comp$Rank=="order",-c(1,2)])) / t(met.18S)
fam.18S = as.data.frame(t(met.18S.comp[met.18S.comp$Rank=="family",-c(1,2)])) / t(met.18S)
genus.18S = as.data.frame(t(met.18S.comp[met.18S.comp$Rank=="genus",-c(1,2)])) / t(met.18S)

phyla.COI = as.data.frame(t(taxa.COI.comp[taxa.COI.comp$Rank=="phylum",-c(1,2)]))
class.COI = as.data.frame(t(taxa.COI.comp[taxa.COI.comp$Rank=="class",-c(1,2)]))
order.COI = as.data.frame(t(taxa.COI.comp[taxa.COI.comp$Rank=="order",-c(1,2)]))
fam.COI = as.data.frame(t(taxa.COI.comp[taxa.COI.comp$Rank=="family",-c(1,2)]))
genus.COI = as.data.frame(t(taxa.COI.comp[taxa.COI.comp$Rank=="genus",-c(1,2)]))

# Make merged datasets with all types of data
phyla.all.mbc = merge(phyla.18S,phyla.COI, all=T, sort=F)
phyla.all = merge(phyla.all.mbc,mPhylum.comp, all=T, sort=F)
phyla.all[is.na(phyla.all)] = 0
dim(phyla.all) # 29 phyla (out of 26+20+15)

class.all.mbc = merge(class.18S,class.COI, all=T, sort=F)
class.all = merge(class.all.mbc,mClass.comp, all=T, sort=F)
class.all[is.na(class.all)] = 0
dim(class.all) # 82 class (out of 57+58+25 -> 58 shared)

order.all.mbc = merge(order.18S,order.COI, all=T, sort=F)
order.all = merge(order.all.mbc,mOrder.comp, all=T, sort=F)
order.all[is.na(order.all)] = 0
dim(order.all) # 227 order (out of 115+127+69 -> 84 shared)

# fam.all.mbc = merge(fam.18S,fam.COI, all=T, sort=F)
# fam.all = merge(fam.all.mbc,mFam.comp, all=T, sort=F)
# fam.all[is.na(fam.all)] = 0
# dim(fam.all) # 316 fam (out of 16+176+212) -> 88 shared)
# 
# genus.all = merge(genus.COI,mGenus.comp, all=T, sort=F)
# genus.all[is.na(genus.all)] = 0
# dim(genus.all) # 521 genus (out of 222 COI +401 morpho) -> 102 shared) 
# #nothing btw 18S (n=1) and COI breaks merge

gi.m = data.frame(row.names=c(row.names(phyla.18S),row.names(phyla.COI),
                              row.names(mPhylum.comp)),
                  groups=c(rep("18S",92),rep("COI",92),rep("Morphology",92)))

row.names(phyla.all) = row.names(class.all) = row.names(order.all) = row.names(gi.m)

# Taxonomy plots side by side
pdf("img/Taxa_comparison_phyla.pdf",height=8,width=50)
taxaplot(20,gi.m,phyla.all)
dev.off()

pdf("img/Taxa_comparison_class.pdf",height=8,width=50)
taxaplot(30,gi.m,class.all)
dev.off()

pdf("img/Taxa_comparison_order.pdf",height=8,width=50)
taxaplot(30,gi.m,order.all)
dev.off()


# Numerically..

## Return row.names by groups: 18S only, COI only, morpho only, 
## 18S_COI, 18S_morpho, COI_morpho, shared_all
shareInfo = function(taxa, groups){
  sInfo = list()
  
  ssu = colSums(taxa[groups=="18S",])
  coi = colSums(taxa[groups=="COI",])
  m = colSums(taxa[groups=="Morphology",])
  
  row.names=c("18S only", "COI only", "morpho only", 
                                 "18S_COI", "18S_morpho", "COI_morpho", "shared_all")
  
  sInfo$Unique_to_18S = names(taxa)[ssu>0 & coi==0 & m==0]
  sInfo$Unique_to_COI = names(taxa)[ssu==0 & coi>0 & m==0]
  sInfo$Unique_to_morpho = names(taxa)[ssu==0 & coi==0 & m>0]
  sInfo$Shared_18S_COI = names(taxa)[ssu>0 & coi>0 & m==0]
  sInfo$Shared_18S_morpho = names(taxa)[ssu>0 & coi==0 & m>0]
  sInfo$shared_COI_morpho = names(taxa)[ssu==0 & coi>0 & m>0]
  sInfo$shared_betwen_all = names(taxa)[ssu>0 & coi>0 & m>0]
  
  return(sInfo)
}

print(shareInfo(phyla.all,gi.m))
# $Unique_to_18S (n=7)
# [1] "Tardigrada"   "Placozoa"   "Ctenophora"  "Xenacoelomorpha"   "Kinorhyncha" "Loricifera"                        
# [7] "Rhopaluridae phylum incertae sedis"
# 
# $Unique_to_COI
# [1] "Onychophora"
# 
# $Unique_to_morpho (n=6)
# [1] "Phoronida" "Sipuncula"
# 
# $Shared_18S_COI
# [1] "Rotifera"        "Porifera"        "Bryozoa"         "Gastrotricha"    "Gnathostomulida" "Entoprocta"     
# 
# Nothing shared between specific metabarcoding and morpho
# 
# $shared_betwen_all (n=13 ca. 45%, 14 unique to mbc, 1 to morphology)
# [1] "Arthropoda"      "Nematoda"        "Chordata"        "Platyhelminthes" "Mollusca"        "Annelida"        "Cnidaria"       
# [8] "Echinodermata"   "Brachiopoda"     "Hemichordata"    "Priapulida"      "Nemertea"        "Chaetognatha"   


print(shareInfo(class.all,gi.m))
# $Unique_to_18S
# [1] "Enoplea"                                 "Catenulida"                              "Cestoda"                                
# [4] "Seriata"                                 "Macrostomida"                            "Rhabdocoela class incertae sedis"       
# [7] "Haplopharyngida class incertae sedis"    "Lecithoepitheliata class incertae sedis" "Prolecithophora class incertae sedis"   
# [10] "Protobranchia"                           "Oligochaeta order incertae sedis"        "Myxosporea"                             
# [13] "Phoroniformea"                           "Heterotardigrada"                        "Calcarea"                               
# [16] "Typhlocoela"                             "Priapulidae class incertae sedis"        "Homalorhagida class incertae sedis"     
# [19] "Cyclorhagida class incertae sedis"       "Heteronemertea"                          "Palaeonemertea"                         
# [22] "Macrodasyida class incertae sedis"       "Chaetonotida class incertae sedis"       "Bursovaginoidea class incertae sedis"   
# [25] "Coloniales class incertae sedis"         "Rhopaluridae class incertae sedis"      
# 
# $Unique_to_COI
# [1] "Diplopoda"           "Chilopoda"           "Elasmobranchii"      "Staurozoa"           "Entoprocta cla"      "Gastrotricha cla"   
# [7] "Gnathostomulida cla" "Cephalopoda"         "Secernentea"         "Adenophorea"         "Anopla"              "Onychophorida"      
# [13] "Rhabditophora"       "Hexactinellida"      "Priapulimorpha"      "Bdelloidea"         
# 
# $Unique_to_morpho
# [1] "Crustacea"          "Phoronidae (class)" "Phascolosomatidea"  "Sipunculidea"      
# 
# $Shared_18S_COI
# [1] "Maxillopoda"      "Branchiopoda"     "Chromadorea"      "Thaliacea"        "Trematoda"        "Polyplacophora"   "Monogononta"     
# [8] "Hydrozoa"         "Demospongiae"     "Gymnolaemata"     "Stenolaemata"     "Acoela"           "Nemertodermatida" "Enopla"          
# [15] "Sagittoidea"     
# 
# $Shared_18S_morpho
# character(0)
# 
# $shared_COI_morpho
# [1] "Hexanauplia"   "Pycnogonida"   "Clitellata"    "Crinoidea"     "Solenogastres"
# 
# $shared_betwen_all (n=16, 21 just not in 18S = 26%, 57 unique to mbc, 4 to morpho)
# [1] "Malacostraca"   "Ostracoda"      "Ascidiacea"     "Gastropoda"     "Bivalvia"       "Scaphopoda"     "Caudofoveata"   "Polychaeta"    
# [9] "Anthozoa"       "Scyphozoa"      "Echinoidea"     "Holothuroidea"  "Asteroidea"     "Ophiuroidea"    "Rhynchonellata" "Enteropneusta" 

print(shareInfo(order.all,gi.m))
#$Unique_to_18S
# [1] "Monstrilloida"                         "Tantulocarida"                         "Eucarida"                             
# [4] "Peracarida"                            "Eumalacostraca"                        "Mermithida"                           
# [7] "Triplonchida"                          "Ascaridida"                            "Desmodorida"                          
# [10] "Araeolaimida"                          "Desmoscolecida"                        "Salpida"                              
# [13] "Catenulidae order incertae sedis"      "Trypanorhyncha"                        "Tetraphyllidea"                       
# [16] "Proseriata"                            "Azygiida"                              "Prosorhynchoides order incertae sedis"
# [19] "Rhabdocoela"                           "Haplopharyngida"                       "Lecithoepitheliata"                   
# [22] "Prolecithophora"                       "Heterobranchia"                        "Myoida"                               
# [25] "Pectinoida"                            "Mytiloida"                             "Arcoida"                              
# [28] "Veneroida (Heteroconchia)"             "Neoloricata"                           "Nuculoida"                            
# [31] "Ploimida"                              "Flosculariacea"                        "Protodrilus order incertae sedis"     
# [34] "Capitellida"                           "Dinophilida"                           "Opheliidae order incertae sedis"      
# [37] "Orbiniidae order incertae sedis"       "Paraonidae order incertae sedis"       "Palpata Incertae Sedis"               
# [40] "Flabelligerida"                        "Echiuroinea"                           "Actiniaria (Anthozoa)"                
# [43] "Narcomedusae (Hydrozoa)"               "Bivalvulida"                           "Coronatae"                            
# [46] "Tetractinellida"                       "Suberitida"                            "Bubarida"                             
# [49] "Axinellida"                            "Biemnida"                              "Verongiida"                           
# [52] "Tethyida"                              "Leucosolenida"                         "Clathrinida"                          
# [55] "Harrimaniidae order incertae sedis"    "Ptychoderidae order incertae sedis"    "Ctenostomatida"                       
# [58] "Tubuliporida"                          "Cydippida"                             "Priapulidae order incertae sedis"     
# [61] "Homalorhagida"                         "Cyclorhagida"                          "Hoplonemertea"                        
# [64] "Macrodasyida"                          "Chaetonotida"                          "Aphragmophora"                        
# [67] "Coloniales"                            "Rhopaluridae order incertae sedis"    
# 
# $Unique_to_COI
# [1] "Poecilostomatoida"                     "Pantopoda"                             "Enchytraeida"                         
# [4] "Opheliida"                             "Canalipalpata"                         "Cyclostomata"                         
# [7] "Carcharhiniformes"                     "Narcomedusae"                          "Rhizostomeae"                         
# [10] "Comatulida"                            "Echinoida"                             "Camaradonta"                          
# [13] "Forcipulatida"                         "Valvatida"                             "Ostreoida"                            
# [16] "Bivalvia or"                           "Euheterodonta"                         "Limoida"                              
# [19] "Sorbeoconcha"                          "Pulmonata"                             "Basommatophora"                       
# [22] "Stylommatophora"                       "Architaenioglossa"                     "Gymnosomata"                          
# [25] "Anaspidea"                             "Mesogastropoda"                        "Sepiida"                              
# [28] "Octopoda"                              "Chitonida"                             "Lepidopleurida"                       
# [31] "Rhabditida"                            "Strongylida"                           "Aphelenchida"                         
# [34] "Tylenchida"                            "Trefusiida"                            "Palaeonemertea"                       
# [37] "Anopla or"                             "Anura or"                              "Heteronemertea"                       
# [40] "Euonychophora"                         "Tricladida"                            "Macrostomida"                         
# [43] "Platyhelminthes or"                    "Childiidae order incertae sedis"       "Haploposthiidae order incertae sedis" 
# [46] "Convolutidae order incertae sedis"     "Diopisthoporidae order incertae sedis" "Isodiametridae order incertae sedis"  
# [49] "Mecynostomidae order incertae sedis"   "Dakuidae order incertae sedis"         "Ascopariidae order incertae sedis"    
# [52] "Spirophorida"                          "Astrophorida"                          "Verongida"                            
# [55] "Halichondrida"                         "Hadromerida"                           "Homosclerophorida"                    
# [58] "Priapulimorphida"                      "Ploima"                               
# 
# $Unique_to_morpho
# [1] "Echiuroidea"        "Echiura"            "Calanoida"          "Thoracica"          "Leptostraca"        "Tanaidacea"        
# [7] "Spirularia"         "Camarodonta"        "Amphilepidida"      "Ophiacanthida"      "Arcida"             "Carditida"         
# [13] "Galeommatida"       "Limida"             "Lucinida"           "Nuculanida"         "Pectinida"          "Anomalodesmata "   
# [19] "Lepetellida"        "Pteropoda"          "Trochida"           "Neomeniamorpha"     "Phoronidae (order)" "Priapulomorpha"    
# [25] "Aspidosiphonida"   
# 
# $Shared_18S_COI
# [1] "Harpacticoida"                          "Diplostraca"                            "Enoplida"                              
# [4] "Monhysterida"                           "Chromadorida"                           "Enterogona"                            
# [7] "Plagiorchiida"                          "Veneroida"                              "Lucinoida"                             
# [10] "Nuculanoida"                            "Chaetodermatida"                        "Haplotaxida"                           
# [13] "Zoantharia"                             "Ceriantharia"                           "Alcyonacea"                            
# [16] "Anthoathecata"                          "Siphonophorae"                          "Limnomedusae"                          
# [19] "Leptothecata"                           "Trachymedusae"                          "Semaeostomeae"                         
# [22] "Stauromedusae"                          "Haplosclerida"                          "Poecilosclerida"                       
# [25] "Cheilostomatida"                        "Cyclostomatida"                         "Nemertodermatidae order incertae sedis"
# [28] "Monostilifera"                          "Phragmophora"                           "Bursovaginoidea"                       
# 
# $Shared_18S_morpho
# [1] "Podocopida"      "Myodocopida"     "Stolidobranchia" "Gadilida"        "Golfingiida"     "Terebratulida"  
# 
# $shared_COI_morpho
# [1] "Amphipoda"       "Decapoda"        "Mysida"          "Isopoda"         "Euphausiacea"    "Cumacea"         "Polychaeta or"  
# [8] "Amphinomida"     "Scleractinia"    "Actiniaria"      "Spatangoida"     "Clypeasteroida"  "Apodida"         "Dendrochirotida"
# [15] "Paxillosida"     "Ophiurida"       "Mytilida"        "Venerida"        "Cardiida"        "Nuculida"        "Adapedonta"     
# [22] "Anomalodesmata"  "Carditoida"      "Littorinimorpha" "Gastropoda or"   "Neogastropoda"   "Nudibranchia"    "Cephalaspidea"  
# 
# $shared_betwen_all (10)
# [1] "Sessilia"        "Caenogastropoda" "Pholadomyoida"   "Dentaliida"      "Spionida"        "Terebellida"     "Phyllodocida"   
# [8] "Sabellida"       "Eunicida"        "Pennatulacea"   




# ----- PCA based on physicochemical parameters (for Fig 4) ----- 

pcPCA = princomp(pc.comp)
summary(pcPCA) #65% explained by PCA1 (42%) and PCA2 (23%)
biplot(pcPCA, col="purple",cex=c(0.01,.8), xlab="PC1", ylab="PC2")#,xlim=range(-.1,.3),ylim=range(-.2,.1))
points(pcPCA$scores[,c(1:2)], cex=.7, col=colors[as.numeric(md.comp$Platform)],
       pch=as.numeric(md.comp$Platform))
ordihull(pcPCA,md.comp$Region,label=T,col="black", kind="sd",draw="line")

## Control with only north-east geo distance to BC
bc=unlist(vegdist(otus.18S.np.p.ra.ab))
ne = unlist(dist(md.18S.p$Northeast))
propergeodist = unlist(dist(md.18S.p[,c("WGS84E","WGS84N")]))
summary(lm(ne~bc)) #R2=.34
summary(lm(log(ne+1e-6)~bc)) #R2=.53
summary(lm(propergeodist~bc)) #R2=.34
summary(lm(log(propergeodist+1e-6)~bc)) #R2=.55
