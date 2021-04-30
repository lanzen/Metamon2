# Andes Lanzen 2021

## This script that was used for analysis of alpha diversity and multivariate analysis
## based on pairwise community dissimilarity, of 18S metabarcoding data


# Change to correct path for local directory (github clone) and uncomment:
# setwd("/home/alanzen/projects/Metamon2/")

require(vegan)
source('R/utils/heterogeneity_rarefaction_functions.R')
source('R/utils/R/filtering.R')
source('R/utils/R/diversity.r')
source('R/utils/R/correlationTests.r')
source('R/utils/R/taxaplot.R')
source('R/utils/R/mergeOTUTable.R')
source('R/utils/calculatePIs.R')

# Read data from plankton filtered OTU table in directory CREST_Filtered
otus.18S.np = as.data.frame(t(read.delim("SWARM_WP1_20200826_18S/SWARM_table_noPlankton.csv",
                                      row.names=1,header=T,sep="\t")))

# Read metadata that is common for COI and 18S, limit to 18S samples and sort alphabetically
md.all = read.csv(file="WP1_metadata.csv",header=T,row.names=5)
mdR.18S = md.all[row.names(otus.18S.np),]
dim(mdR.18S) # 296 samples

dim(otus.18S.np) #25,319

# ----------- Abundance based filtering and relative abundance calc. --------

# Relative abundance
otus.ra.18S = decostand(otus.18S.np,method="total")
summary(rowSums(otus.18S.np))

# Min read depth = 6k
# Remove taxa with abundance < .0005 in at least one sample (~3 reads in sample w lowest read depth)
otus.ra.18S.ab = dropRareByMaxAbundance(otus.ra.18S,5e-4)
dim(otus.ra.18S.ab) # 4505 OTUs (out of 25k)

## Write abundance filtered OTUs to disk with taxonomy
otus.ra.18S.ab_wTaxonomy = as.data.frame(t(otus.ra.18S.ab))
otus.ra.18S.ab_wTaxonomy$classification = tax.18S[row.names(otus.ra.18S.ab_wTaxonomy),
                                                  "classification"]

write.table(otus.ra.18S.ab_wTaxonomy, "SWARM_WP1_20200826_18S/SWARM_table_noPlankton_abundance_filtered.csv", 
            quote=F, sep="\t")

# --------- Reformat metadata in individual replicates ------------

mdR.18S$TotalPAH = mdR.18S$PAH*1000 #Conversion to ppb
mdR.18S$TotalHC = mdR.18S$THC

# Calculate Pressure Index (PI) values
mdR.18S = calculatePIs(mdR.18S)

mdR.18S$Platform = droplevels(mdR.18S$Platform)
summary(mdR.18S)

# ----- Pool samples from same station -----

# Pool and keep individual extraction pools from same station as replicates
mdR.18S$StationRep = paste(mdR.18S$Station,mdR.18S$Type,sep="")
# Pool OTUs from the dataset before filtering rare OTUs likewise
otus.18S.np.p = mergeOTUTable(otus.18S.np,mdR.18S,by="StationRep")
summary(rowSums(otus.18S.np.p)) #Min 36k reads
# Write pooled OTUs to table
write.csv(as.data.frame(t(otus.18S.np.p)),"SWARM_WP1_20200826_18S/OTUs_pooled_noPlankton.csv",
          quote=F)

# Calculate relative abundance
otus.18S.np.pra = decostand(otus.18S.np.p, method="total")

# Filter rare OTUs that do not have abundance > 0.01 % in at least one sample 
# (corresponding to >3 reads in the pool with the smallest sequencing depth)
otus.18S.np.p.ra.ab=decostand(dropRareByMaxAbundance(otus.18S.np.pra,1E-4),method="total") 
dim(otus.18S.np.p.ra.ab) #6685  OTUs remain of 25k

# Pool metadata
md.18S.p = mergeMetadata(mdR.18S, by="StationRep")
# Write to disk and read, for automatic column class detection
write.csv(md.18S.p,"WP1_Metadata_pooled_18S_np.csv")
md.18S.p=read.csv("WP1_Metadata_pooled_18S_np.csv",row.names=1,header=T)

# Reformat metadata
md.18S.p$TotalPAH = md.18S.p$PAH*1000 #Conversion to ppb
md.18S.p$TotalHC = md.18S.p$THC
md.18S.p$TotalPCB = NA
md.18S.p$TotalDDT = NA
md.18S.p$TotalHCH = NA
md.18S.p$ADE = NA
md.18S.p$OM = NA#md.18S.p$TOC
md.18S.p = calculatePIs(md.18S.p)
md.18S.p$NSIneg = -md.18S.p$NSI

# Calculate means of PI (including sub-types)
mdR.18S$pi = NA
mdR.18S$piHC = NA
mdR.18S$piMetals = NA
for (st in unique(mdR.18S$StationRep)){
  mSt = md.18S.p[md.18S.p$StationRep==st,]
  mdR.18S[mdR.18S$StationRep == st,]$pi = mSt$pi
  mdR.18S[mdR.18S$StationRep == st,]$piHC = mSt$piHC
  mdR.18S[mdR.18S$StationRep == st,]$piMetals = mSt$piMetals
}

md.18S.p$Northeast = md.18S.p$WGS84E + md.18S.p$WGS84N

summary(md.18S.p)

# Control equivalency to OTU table
table(row.names(md.18S.p) == row.names(otus.18S.np.p.ra.ab))


# Write pooled metadata for networks
write.csv(mdR.18S, "Metadata_all_reps_w_PIs_from_pooled.csv")


# ---------- Create groups of parameters with same occurence throughout for envfit and similar ------

#13dp:

mp1 = md.18S.p[,c("NPD","TotalPAH","Acenaften","Acenaftylen" ,"Antracen","Benz.a.antracen",
                  "Benzo.a.pyren","Benzo.bjk.fluoranten","Benzo.ghi.perylen",
                  "C1.alkyldibenzotiofener","C1.alkylfenantrener.antracener",
                  "C1.alkylnaftalener","C2.alkyldibenzotiofener","C2.alkylfenantrener.antracener",
                  "C3.alkyldibenzotiofener","C3.alkylfenantrener.antracener",
                  "C3.alkylnaftalener","Chrysen","Dibenz.ah.antracen","Dibenzotiofen",
                  "Fenantren","Fluoranten","Fluoren","Indeno.1.2.3.c.d.pyren",
                  "Naftalen","Pyren")]
# 54 missing

mp2 = md.18S.p[,c("As","Ba","Cr","Cd","Cu","Hg","Pb","Zn","piMetals","pi", "piHC","THC",
                  "TotalHC","Northeast","Depth","NSI")] #0 missing

mp4 = md.18S.p[,c("Kurtosis","Pelite","Sand","Skewness",
             "Sorting","TOC","Grus","Kornstorrelse")] #2 missing

mpdiv = md.18S.p[,c("MorphoS","MorphoH","ES100","Individtetthet",
                    "ISI","NSI","NQI1")] #4 missing

mp5 = md.18S.p[,c("Extract.conc","Distance")] #21 missing

## NMDS analysis of pooled samples

nmds.p = metaMDS(otus.18S.np.p.ra.ab)
ordiplot(nmds.p,type="none")#,xlim=range(-3,5),ylim=range(-1,1))
points(nmds.p,pch=as.numeric(md.18S.p$Platform),col=colors[as.numeric(md.18S.p$Platform)],cex=.8)
legend("bottomright",pch=as.numeric(sort(unique(md.18S.p$Platform))),
       col=colors[as.numeric(sort(unique(md.18S.p$Platform)))],cex=.6, 
       legend=sort(unique(md.18S.p$Platform)), ncol=4)
#text(nmds.p,pos=1,cex=.3,col=colors[as.numeric(md.18S.p$Platform)])
ordihull(nmds.p,md.18S.p$Region,label=T,col="black", kind="sd",draw="line")

## Fitting of physicochemical parameters to NMDS
efp1=envfit(nmds.p,mp1,na.rm=T)


efp2 = envfit(nmds.p,mp2)
# As         0.38201 -0.92416 0.1280  0.001 ***
#   Ba         0.38082 -0.92465 0.3899  0.001 ***
#   Cr        -0.59170 -0.80616 0.0693  0.027 *  
#   Cd         0.23482 -0.97204 0.0682  0.035 *  
#   Cu         0.41905 -0.90796 0.0885  0.010 ** 
#   Hg        -0.24673 -0.96908 0.0805  0.010 ** 
#   Pb        -0.48385 -0.87515 0.1162  0.005 ** 
#   Zn        -0.02332 -0.99973 0.0759  0.025 *  
#   piMetals   0.27301 -0.96201 0.4203  0.001 ***
#   pi         0.47278 -0.88118 0.3791  0.001 ***
#   piHC       0.74691 -0.66492 0.3262  0.001 ***
#   THC        0.78867 -0.61482 0.1188  0.005 ** 
#   TotalHC    0.78867 -0.61482 0.1188  0.005 ** 
#   Northeast -0.76884 -0.63944 0.7058  0.001 ***
#   Depth     -0.83586 -0.54894 0.6663  0.001 ***
#   NSI       -0.88038  0.47426 0.4584  0.001 ***

efp4 = envfit(nmds.p,mp4,na.rm=T)
# Kurtosis      -0.00061  1.00000 0.0938  0.011 *  
#   Pelite        -0.34409 -0.93894 0.3743  0.001 *** 
#   Sand           0.09054 -0.99589 0.5704  0.001 *** <--
#   Skewness       0.17167 -0.98515 0.1462  0.001 ***
#   Sorting       -0.82326 -0.56766 0.1083  0.004 ** 
#   TOC            0.83261 -0.55385 0.0962  0.022 *  
#   piOM           0.74144 -0.67102 0.1114  0.004 ** 
#   Grus          -0.28163 -0.95952 0.2018  0.001 ***
#   Kornstorrelse -0.99619 -0.08718 0.1623  0.002 ** 

efp5 = envfit(nmds.p,mp5,na.rm=T)
# NMDS1    NMDS2     r2 Pr(>r)   
# Extract.conc  0.97094  0.23932 0.1589  0.006 **
#   Distance     -0.33233  0.94316 0.1418  0.008 **

efp6 = envfit(nmds.p,mpdiv,na.rm=T)
# MorphoS        -0.43938  0.89830 0.2597  0.001 ***
#   MorphoH        -0.72640  0.68728 0.4444  0.001 ***
#   ES100          -0.80579  0.59220 0.4573  0.001 ***
#   Individtetthet  0.38756  0.92185 0.0562  0.074 .  
# ISI            -0.95152  0.30758 0.4209  0.001 ***
#   NSI            -0.88815  0.45955 0.4860  0.001 ***
#   NQI1           -0.68752  0.72616 0.4349  0.001 ***
# ---

## Add fitted parameters to plot
#plot(efp1,p.max=.001,cex=.5,col="purple")
plot(efp2,p.max=.001,cex=.5,col="purple")
plot(efp4,p.max=.001,cex=.5,col="purple")
plot(efp5,p.max=.001,cex=.5,col="purple")
plot(efp6,p.max=.001,cex=.5,col="green")

# ---- Taxonomy -------

## Read taxonomy data from CREST

taxa.all.18S = read.table("SWARM_WP1_20200826_18S/CREST_NoPlankton/Relative_Abundance.tsv",
                          sep="\t", header=T,row.names=3, check.names = F)

taxa.p.18S = read.table("SWARM_WP1_20200826_18S/CREST_Pooled_NoPlankton/Relative_Abundance.tsv",
                          sep="\t", header=T,row.names=3, check.names = F)

ass.all.18S = read.table("SWARM_WP1_20200826_18S/CREST_NoPlankton/All_Assignments.tsv",
                          sep="\t", header=T,row.names=3, check.names = F)
ass.p.18S = read.table("SWARM_WP1_20200826_18S/CREST_Pooled_NoPlankton/All_Assignments.tsv",
                         sep="\t", header=T,row.names=3, check.names = F)

tra = taxa.all.18S[,c("Rank","Taxonpath",row.names(otus.18S.np))]
trp = taxa.p.18S[,c("Rank","Taxonpath",row.names(md.18S.p))]

## Calculate relative abundance of assignments (mix of ranks with lowest possible)
asra = decostand(as.data.frame(t(ass.all.18S[,row.names(otus.18S.np)])),method="total")
asp_18S_np = decostand(as.data.frame(t(ass.p.18S[,row.names(md.18S.p)])),method="total")
asp_18S_np = asp_18S_np[,colSums(asp_18S_np)>0]
asp = asp_18S_np
asra = asra[,colSums(asra)>0]

## Number of unique taxa across ranks
summary(tra$Rank)
# class       domain       family        genus      kingdom 
# 181            1          179          141           21
# meta        order       phylum         root      species 
# 1          291          83            1           50
# superkingdom 
# 14 

## Prepare data frames for taxonomic plots
grouping_info<-data.frame(row.names=row.names(mdR.18S), mdR.18S$Station)

ranks = data.frame(rank=c("superkingdom","kingdom","phylum","class","order","family","genus",
                          "species"),
                   levels=c(8,11,25,25,25,25,25,25))

## Generate relative abundance barcharts, all replicates
for (i in c(1:7)){
  r=as.character(ranks$rank[i])
  levelTaxa = as.data.frame(t(tra[tra$Rank==r,-c(1:2)]))
  row.names(levelTaxa) = row.names(mdR.18S)
  print(paste(mean(rowSums(levelTaxa))*100,"% classified at rank",r))
  pdf(paste("img/18S/taxon_barplots/",r,"_all_reps_noPlankton.pdf",sep=""),height=8,width=50)
  taxaplot(ranks$levels[i],grouping_info,levelTaxa)
  dev.off()
}

## Generate relative assignments (lowest) barcharts, all replicates
pdf("img/18S/Assignments_all_reps_noPlankton.pdf",height=8,width=50)
taxaplot(30,grouping_info,asra)
dev.off()

## Prepare data frames for taxonomic plots, pooled replicates
grouping_info<-data.frame(row.names=row.names(mdR.18S), mdR.18S$Station)

ranks = data.frame(rank=c("superkingdom","kingdom","phylum","class","order","family","genus",
                          "species"),
                   levels=c(8,11,25,25,25,25,25,25))
gi=data.frame(row.names=row.names(md.18S.p), md.18S.p$Platform)

## Generate relative abundance barcharts, pooled replicates

for (i in c(1:8)){
  r=as.character(ranks$rank[i])
  levelTaxa = as.data.frame(t(trp[trp$Rank==r,-c(1:2)]))
  row.names(levelTaxa) = row.names(md.18S.p)
  print(paste(mean(rowSums(levelTaxa))*100,"% classified at rank",r))
  pdf(paste("img/18S/taxon_barplots/",r,"_noPlankton.pdf",sep=""),height=8,width=15)
  taxaplot(ranks$levels[i],gi,levelTaxa)
  dev.off()
}

## Generate relative assignments (lowest) barcharts, pooled replicates

pdf("img/18S/Assignments_noPlankton.pdf",height=8,width=15)
taxaplot(30,gi,asp_18S_np)
dev.off()

# 1] "99.8994178483063 % classified at rank kingdom"
# [1] "82.4099173637112 % classified at rank phylum"
# [1] "77.4226658042589 % classified at rank class"
# [1] "56.069826920632 % classified at rank order"
# [1] "8.57348541568563 % classified at rank family"
# [1] "4.74026677754199 % classified at rank genus"
# [1] "1.05027270163993 % classified at rank species"


# -------- Prediction of BI values using BBI package -------

require(BBI)
mt = ass.p.18S[,-1]
mt$Taxonpath = row.names(mt)
bbi.18S = BBI(mt) # Found match : 56  Not found : 907 ====

bbi.18S$BBI
table(row.names(bbi.18S$BBI) == row.names(md.18S.p))
lm.nsi.18S = lm(bbi.18S$BBI[,3]~md.18S.p$NSI)
summary(lm.nsi.18S)
# R2 near 0

pdf("img/18S/NSI_18S_v_NSI.pdf",width=8,height=5)
plot_ml(data = bbi.18S$BBI[,3], metadata = md.18S.p, 
        xIndex = "NSI", yIndex = "NSI", aggreg = c("Station","Platform"),
        title = "NSI predicted from 18S metabarcoding")
dev.off()
# 18S does not work at all, likely due to lack of specific assignments


# ---------Pooled alpha div ------

# Write diversity statistics to table
writeDivStats("QCDiversity_pooled_18S_np.csv", otus.18S.np.p)
div.p.18S = read.csv("QCDiversity_pooled_18S_np.csv", row.names=1)
div.p.18S = div.p.18S[,c("Rarefied.richness","H","J")]
printVS(div.p.18S, mp1 ,a=0.01)
printVS(div.p.18S, mp2 ,a=0.01)


# ----- Figures for manuscript ---------

## Morphotaxa diversity ~ PI
mpi = lm(MorphoH~pi, data=md.18S.p) 
summary(mpi) # R2=.42, p<E-12 F=69 (HC: R2=.53, F=109)
plot(MorphoH~pi, data=md.18S.p, xlab="PI",ylab="H' morphospecies",
     pch=md.18S.p$col_plot) #col=md.18S.p$col_plot, 
abline(mpi,col="darkgrey")

legend("topright",legend=sort(unique(md.COI.p$Platform)),
       #col=sort(unique(md.COI.p$col_plot)),
       pch=sort(unique(md.COI.p$col_plot)),cex=.9)
text(0, 0.1, "R² = .42, p < 1E-12", pos=4, col="#606060")

## 18S diversity ~ PI
mpSSURS = lm(div.p.18S$Rarefied.richness~md.18S.p$pi)  
summary(mpSSURS) # .30, 44 for pi~RS; (R2=.21, p<E-6 F=28 for H~PI
#                                     .36, 57 for piHC~RS)
plot(div.p.18S$Rarefied.richness~md.18S.p$pi, xlab="PI",ylab="Rarefied richness 18S",
     pch=md.18S.p$col_plot, ylim=range(0:2000))
abline(mpSSURS,col="darkgrey")
text(0, 40, "R² = .30, p < 1E-8", pos=4,  col="#606060")

## Morphotaxa diversity ~ 18S 

morphoSSURS = lm(div.p.18S$Rarefied.richness~md.18S.p$MorphoH)  
summary(morphoSSURS) # R=.43, p<E-12, F=70
plot(div.p.18S$Rarefied.richness~md.18S.p$MorphoH, 
     xlab="H' morphospecies",ylab="Rarefied richness 18S",
     pch=md.18S.p$col_plot, ylim=range(0:2000))
abline(morphoSSURS,col="darkgrey")
text(0, 40, "R² = .42, p < 1E-12", pos=4,  col="#606060")

