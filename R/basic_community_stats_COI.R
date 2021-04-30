# Andes Lanzen 2021

## This script that was used for analysis of alpha diversity and multivariate analysis
## based on pairwise community dissimilarity, of COI metabarcoding data

# Change to correct path for local directory (github clone) and uncomment:
# setwd("/home/alanzen/projects/Metamon2/")


require(vegan)
source('R/utils/heterogeneity_rarefaction_functions.R')
source('R/utils/filtering.R')
source('R/utils/diversity.r')
source('R/utils/correlationTests.r')
source('R/utils/taxaplot.R')
source('R/utils/mergeOTUTable.R')
source('R/utils/calculatePIs.R')


# Read data from plankton filtered OTU table in directory CREST_Filtered
otusR.COI = read.delim("SWARM_WP1_20200722_COI/SWARM_table_noPlankton.csv",
                                      row.names=1,header=T,sep="\t")
otusR.COI = as.data.frame(t(otusR.COI))

# Read metadata that is common for COI and 18S, limit to COI samples and sort alphabetically
md.all = read.csv(file="WP1_metadata.csv",header=T,row.names=5)
mdR.COI = droplevels(md.all[row.names(otusR.COI),])
dim(mdR.COI) # 289 samples

# ----------- Abundance based filtering and relative abundance calc. --------

# Relative abundance
otus.ra.COI = decostand(otusR.COI,method="total")
dim(otus.ra.COI) #8186

# Remove taxa with abundance < .003 in at least one sample (~3 reads in sample w lowest read depth)
otus.ra.COI.ab = dropRareByMaxAbundance(otus.ra.COI,3e-3) #Min read depth = 1k
otus.h.COI = decostand(otus.ra.COI.ab,method="hell")
dim(otus.ra.COI.ab) # 1857 OTUs remaining out of 8186

## Write abundance filtered OTUs to disk with taxonomy
otus.ra.COI.ab_wTaxonomy = as.data.frame(t(otus.ra.COI.ab))
otus.ra.COI.ab_wTaxonomy$classification = tax.COI[row.names(otus.ra.COI.ab_wTaxonomy),
                                                  "classification"]

write.table(otus.ra.COI.ab_wTaxonomy, "SWARM_WP1_20200722_COI/SWARM_table_noPlankton_abundance_filtered.csv", 
            quote=F, sep="\t")


# ----- Pool samples from same station -----

# Pool and keep individual extraction pools from same station as replicates
mdR.COI$StationRep = paste(mdR.COI$Station,mdR.COI$Type,sep="")

# Pool OTUs from the dataset before filtering rare OTUs likewise
otusR.COI.p = mergeOTUTable(otusR.COI,mdR.COI,by="StationRep")
summary(rowSums(otusR.COI.p)) #Min 6957

# Write pooled OTUs to table
write.csv(as.data.frame(t(otusR.COI.p)),quote=F,"SWARM_WP1_20200722_COI/OTUs_pooled.csv")

# Calculate relative abundance
otusR.COI.pra = decostand(otusR.COI.p, method="total")

# Filter rare OTUs that do not have abundance > 0.05 % in at least one sample 
# (corresponding to >3 reads in the pool with the smallest sequencing depth)
otusR.COI.p.ra.ab=decostand(dropRareByMaxAbundance(otusR.COI.pra,5E-4),method="total") 
dim(otusR.COI.p.ra.ab) #2921 OTUs

# Pool metadata
md.COI.p = mergeMetadata(mdR.COI, by="StationRep")

# Write to disk and read, for automatic column class detection
write.csv(md.COI.p,"WP1_Metadata_pooled_COI.csv")
md.COI.p=read.csv("WP1_Metadata_pooled_COI.csv",row.names=1,header=T)

# --------- Reformat metadata in pooled replicates ------------
md.COI.p$TotalPAH = md.COI.p$PAH*1000 #Conversion to ppb
md.COI.p$TotalHC = md.COI.p$THC

# Calculate Pressure Index (PI) values
md.COI.p = calculatePIs(md.COI.p)
md.COI.p$NSIneg = -md.COI.p$NSI

md.COI.p$Northeast = md.COI.p$WGS84E + md.COI.p$WGS84N

summary(md.COI.p)

# Control
table(row.names(md.COI.p) == row.names(otusR.COI.p.ra.ab))

# ---------- Create groups of parameters with same occurence throughout for envfit and similar ------

mp1 = md.COI.p[,c("NPD","TotalPAH","Acenaften","Acenaftylen" ,"Antracen","Benz.a.antracen",
                 "Benzo.a.pyren","Benzo.bjk.fluoranten","Benzo.ghi.perylen",
                 "C1.alkyldibenzotiofener","C1.alkylfenantrener.antracener",
                 "C1.alkylnaftalener","C2.alkyldibenzotiofener","C2.alkylfenantrener.antracener",
                 "C3.alkyldibenzotiofener","C3.alkylfenantrener.antracener",
                 "C3.alkylnaftalener","Chrysen","Dibenz.ah.antracen","Dibenzotiofen",
                 "Fenantren","Fluoranten","Fluoren","Indeno.1.2.3.c.d.pyren",
                 "Naftalen","Pyren")]

# 54 missing (46 d.p.)

mp2 = md.COI.p[,c("As","Ba","Cr","Cd","Cu","Hg","Pb","Zn","piMetals","pi", "piHC","THC",
                  "Northeast","Depth","THC_PI")] 
#0 missing

mp4 = md.COI.p[,c("Kurtosis","Pelite","Sand","Skewness",
                  "Sorting","TOC","Grus","Kornstorrelse")] #4 missing

mpdiv = md.COI.p[,c("MorphoS","MorphoH","ES100","Individtetthet",
                    "ISI","NSI","NQI1")] #4 missing

mp5 = md.COI.p[,c("Extract.conc","Distance")] #21 missing


nmds.p = metaMDS(otusR.COI.p.ra.ab)#, engine="isoMDS")
ordiplot(nmds.p,type="none")#,xlim=range(-2,4),ylim=range(-1,1))
points(nmds.p,pch=as.numeric(md.COI.p$Platform),col=colors[as.numeric(md.COI.p$Platform)],cex=.8)
legend("topright",pch=as.numeric(sort(unique(md.COI.p$Platform))),
       col=colors[as.numeric(sort(unique(md.COI.p$Platform)))],cex=.6, 
       legend=sort(unique(md.COI.p$Platform)), ncol=4)
#text(nmds.p,pos=1,cex=.3,col=colors[as.numeric(md.COI.p$Platform)])
ordihull(nmds.p,md.COI.p$Region,label=T,col="black", kind="sd",draw="line")

efp1=envfit(nmds.p,mp1,na.rm=T)
# NPD                             0.21843 -0.97585 0.3774  0.002 ** 
#   TotalPAH                        0.23645 -0.97164 0.3007  0.002 ** 
#   Acenaften                       0.25638 -0.96658 0.2003  0.063 .  
# Acenaftylen                     0.23210 -0.97269 0.3527  0.002 ** 
#   Antracen                        0.26913 -0.96311 0.3625  0.002 ** 
#   Benz.a.antracen                 0.39970 -0.91664 0.1061  0.105    
# Benzo.a.pyren                   0.80671 -0.59095 0.0225  0.488    
# Benzo.bjk.fluoranten            0.85786 -0.51388 0.0237  0.538    
# Benzo.ghi.perylen               0.74674 -0.66512 0.0124  0.718    
# C1.alkyldibenzotiofener         0.09121 -0.99583 0.3016  0.009 ** 
#   C1.alkylfenantrener.antracener  0.16870 -0.98567 0.3259  0.002 ** 
#   C1.alkylnaftalener              0.19352 -0.98110 0.4303  0.001 ***
#   C2.alkyldibenzotiofener         0.21054 -0.97759 0.2924  0.004 ** 
#   C2.alkylfenantrener.antracener  0.19690 -0.98042 0.3326  0.002 ** 
#   C3.alkyldibenzotiofener         0.13306 -0.99111 0.4992  0.001 *** <--
#   C3.alkylfenantrener.antracener  0.21552 -0.97650 0.3388  0.002 ** 
#   C3.alkylnaftalener              0.24354 -0.96989 0.3750  0.002 ** 
#   Chrysen                         0.22347 -0.97471 0.2625  0.005 ** 
#   Dibenz.ah.antracen              0.72360 -0.69022 0.0110  0.691    
# Dibenzotiofen                   0.19617 -0.98057 0.2881  0.003 ** 
#   Fenantren                       0.20495 -0.97877 0.2714  0.006 ** 
#   Fluoranten                      0.32622 -0.94530 0.2808  0.004 ** 
#   Fluoren                         0.20774 -0.97818 0.3179  0.001 ***
#   Indeno.1.2.3.c.d.pyren          0.90179  0.43218 0.0702  0.211    
# Naftalen                        0.18456 -0.98282 0.5088  0.001 *** <--
#   Pyren                           0.29687 -0.95492 0.3477  0.002 ** 
#PAH_PI                          0.06541 -0.99786 0.5478  0.001 *** <---
  

# Typically stronger effects compared to COI!

efp2 = envfit(nmds.p,mp2)#,na.rm=T)
# NMDS1    NMDS2     r2 Pr(>r)    
# As        0.20212 -0.97936 0.4903  0.001 ***
#   Ba       -0.21087 -0.97751 0.3259  0.001 ***
#   Cr        0.75657 -0.65391 0.4349  0.001 ***
#   Cd        0.35238 -0.93586 0.5131  0.001 *** <--
#   Cu        0.25533 -0.96686 0.4275  0.001 ***
#   Hg        0.61705 -0.78693 0.5441  0.001 *** <--
#   Pb        0.75413 -0.65673 0.3981  0.001 ***
#   Zn        0.48606 -0.87393 0.4198  0.001 ***
#   piMetals -0.17047 -0.98536 0.3090  0.001 ***
# pi         -0.14111 -0.98999 0.4952  0.001 ***
#   piHC     -0.15012 -0.98867 0.5309  0.001 ***
#   THC      -0.24257 -0.97013 0.1415  0.008 ** 
# Northeast  0.97396 -0.22670 0.7009  0.001 ***
#   Depth     0.98569 -0.16855 0.7175  0.001 *** <-- (strong also in COI)
#   THC_PI   -0.18271 -0.98317 0.4638  0.001 *** ***
  
efp4 = envfit(nmds.p,mp4,na.rm=T)
# NMDS1    NMDS2     r2 Pr(>r)    
# Kurtosis      -0.18580  0.98259 0.1429  0.003 ** 
#   Pelite         0.70506 -0.70915 0.6256  0.001 *** <- Now stronger than sand, again..
#   Sand          -0.73824 -0.67453 0.2106  0.001 ***
#   Skewness      -0.54659 -0.83740 0.0499  0.114    
# Sorting        0.85965 -0.51088 0.4146  0.001 ***
#   TOC            0.05141 -0.99868 0.2443  0.004 ** 
#   Grus           0.63936 -0.76891 0.0786  0.046 *  
#   Kornstorrelse  0.98842 -0.15176 0.4768  0.001 ***

efp5 = envfit(nmds.p,mp5,na.rm=T)
# NMDS1    NMDS2     r2 Pr(>r)    
# Extract.conc  0.74335 -0.66890 0.4634  0.001 ***
#   Distance      0.94971  0.31315 0.0953  0.037 *  

efp6 = envfit(nmds.p,mpdiv,na.rm=T)
# NMDS1    NMDS2     r2 Pr(>r)    
# MorphoS         0.15482  0.98794 0.4668  0.001 ***
#   MorphoH         0.56933  0.82211 0.3993  0.001 ***
#   ES100           0.70097  0.71319 0.3822  0.001 ***
#   Individtetthet -0.38078  0.92467 0.1883  0.001 ***
#   ISI             0.39722  0.91772 0.3514  0.001 ***
#   NSI             0.62108  0.78375 0.3959  0.001 ***
#   NQI1            0.34702  0.93786 0.4607  0.001 ***


plot(efp2,p.max=.001,cex=.5,col="purple")
plot(efp4,p.max=.001,cex=.5,col="purple")
plot(efp5,p.max=.001,cex=.5,col="purple")
plot(efp6,p.max=.001,cex=.5,col="green")

# ---- Taxonomy -------

## Read taxonomy data from CREST

taxa.all.COI = read.table("SWARM_WP1_20200722_COI/CREST_NoPlankton/Relative_Abundance.tsv",
                          sep="\t", header=T,row.names=3, check.names = F)

taxa.p.COI = read.table("SWARM_WP1_20200722_COI/CREST_Pooled_NoPlankton/Relative_Abundance.tsv",
                          sep="\t", header=T,row.names=3, check.names = F)

ass.all.COI = read.table("SWARM_WP1_20200722_COI/CREST_NoPlankton/All_Assignments.tsv",
                          sep="\t", header=T,row.names=3, check.names = F)
ass.p.COI = read.table("SWARM_WP1_20200722_COI/CREST_Pooled_NoPlankton/All_Assignments.tsv",
                         sep="\t", header=T,row.names=3, check.names = F)

tra = taxa.all.COI[,c("Rank","Taxonpath",row.names(otusR.COI))]
trp = taxa.p.COI[,c("Rank","Taxonpath",row.names(md.COI.p))]

## Calculate relative abundance of assignments (mix of ranks with lowest possible)
asra = decostand(as.data.frame(t(ass.all.COI[,row.names(otusR.COI)])),method="total")
asp = decostand(as.data.frame(t(ass.p.COI[,row.names(md.COI.p)])),method="total")
asra = asra[,colSums(asra)>0]
asp = asp[,colSums(asp)>0]

## Number of unique taxa across ranks
summary(tra$Rank)
# class       domain       family        genus      kingdom         meta        order 
# 52            1          176          222            1            1          127 
# phylum         root      species superkingdom 
# 20            1          178            1 

## Prepare data frames for taxonomic plots
grouping_info<-data.frame(row.names=row.names(mdR.COI), mdR.COI$Station)

ranks = data.frame(rank=c("superkingdom","kingdom","phylum","class","order","family","genus",
                          "species"),
                   levels=c(1,1,20,25,25,25,25,25))

## Generate relative abundance barcharts, all replicates
for (i in c(3:8)){
  r=as.character(ranks$rank[i])
  levelTaxa = as.data.frame(t(tra[tra$Rank==r,-c(1:2)]))
  row.names(levelTaxa) = row.names(mdR.COI)
  print(paste(mean(rowSums(levelTaxa))*100,"% classified at rank",r))
  pdf(paste("img/COI/taxon_barplots/",r,"_all_reps.pdf",sep=""),height=8,width=50)
  taxaplot(ranks$levels[i],grouping_info,levelTaxa)
  dev.off()
}

## Generate relative assignments (lowest) barcharts, all replicates
pdf("img/COI/Assignments_all_reps.pdf",height=8,width=50)
taxaplot(30,grouping_info,asra)
dev.off()

## Prepare data frames for taxonomic plots, pooled replicates
gi=data.frame(row.names=row.names(md.COI.p), md.COI.p$Platform)

## Generate relative abundance barcharts, pooled replicates

for (i in c(3:8)){
  r=as.character(ranks$rank[i])
  levelTaxa = as.data.frame(t(trp[trp$Rank==r,-c(1:2)]))
  row.names(levelTaxa) = row.names(md.COI.p)
  print(paste(mean(rowSums(levelTaxa))*100,"% classified at rank",r))
  pdf(paste("img/COI/taxon_barplots/",r,".pdf",sep=""),height=8,width=15)
  taxaplot(ranks$levels[i],gi,levelTaxa)
  dev.off()
}
# 1] "100.000000000005 % classified at rank phylum"
# [1] "95.4246401258708 % classified at rank class"
# [1] "65.3354089686583 % classified at rank order"
# [1] "53.2711720704321 % classified at rank family"
# [1] "52.0883206629854 % classified at rank genus"
# [1] "30.4844215278197 % classified at rank species"

## Generate relative assignments (lowest) barcharts, pooled replicates

pdf("img/COI/Assignments.pdf",height=8,width=15)
taxaplot(30,gi,asp)
dev.off()


# -------- Prediction of BI values using BBI package -------

require(BBI)

mt = ass.p.COI[,-1]
mt$Taxonpath = row.names(mt)
bbi.coi = BBI(mt) # 427 found, 353 not found
bbi.coi$BBI
table(row.names(bbi.coi$BBI) == row.names(md.COI.p))
lm.nsi.coi = lm(bbi.coi$BBI[,3]~md.COI.p$NSI)
# Adj R2=.43, p<E-12, F=73

## Plot reported NSI v NSI calculated using metabarcoding data

pdf("img/COI/NSI_COI_v_NSI.pdf",width=8,height=5)
plot_ml(data = bbi.coi$BBI[,3], metadata = md.COI.p, 
        xIndex = "NSI", yIndex = "NSI", aggreg = c("Station","Platform"),
        title = "NSI predicted from COI metabarcoding")
dev.off()

## Plot metabarcoding NSI v PI

pdf("img/COI/NSI_COI_v_PI.pdf",width=8,height=5)
plot_ml(data = -bbi.coi$BBI[,3], metadata = md.COI.p, 
        xIndex = "pi", yIndex = "NSIneg", aggreg = c("Station","Platform"),
        title = "COI predicted NSI v PI")
dev.off()

# ---------Pooled alpha div ------

# Write diversity statistics to table
writeDivStats("QCDiversity_pooled_COI.csv", otusR.COI.p)

div.p.COI = read.csv("QCDiversity_pooled_COI.csv", row.names=1)

# ----- Figures for manuscript ---------

## COI diveristy ~ PI
mpCOIRS = lm(div.p.COI$H~md.COI.p$pi) 
summary(mpCOIRS) # R2=.11, p=5E-4, F=13 (.13, F=15 for piHC, .09, 10 for H')
plot(div.p.COI$Rarefied.richness~md.COI.p$pi, xlab="PI",ylab="Rarefied richness COI",
     pch=md.COI.p$col_plot, ylim=range(0:700)) #col=md.COI.p$col_plot, 
abline(mpCOIRS,col="darkgrey")

# legend("topright",legend=sort(unique(md.COI.p$Platform)),
#        #col=sort(unique(md.COI.p$col_plot)),
#        pch=sort(unique(md.COI.p$col_plot)),cex=.8)
text(0, 10, "R² = .11, p = .0005", cex=.9, pos=4, col="#606060")

## COI diversity ~ morphotaxa
morphoCOIRS = lm(div.p.COI$Rarefied.richness~md.COI.p$MorphoH)  
summary(morphoCOIRS) # R=.25, p<E-6, F=32
plot(div.p.COI$Rarefied.richness~md.COI.p$MorphoH, 
     xlab="H' morphospecies",ylab="Rarefied richness COI",
     pch=md.18S.p$col_plot, ylim=range(0:700))
abline(morphoCOIRS,col="darkgrey")
legend("topleft",ncol=3,legend=sort(unique(md.COI.p$Platform)),
       #col=sort(unique(md.COI.p$col_plot)),
       pch=sort(unique(md.COI.p$col_plot)))#,cex=.8)
text(3.5, 10, "R² = .25, p < 1E-6", pos=4,  col="#606060")

## COI diversity ~ 18S diversity
ssuRSCOI = div.p.18S[row.names(div.p.COI),"Rarefied.richness"]
SSUCOIRS = lm(ssuRSCOI~div.p.COI$Rarefied.richness)  
summary(SSUCOIRS) # R=.54, p<E-15, F=115
plot(ssuRSCOI~div.p.COI$Rarefied.richness, 
     xlab="Rarefied richness COI",ylab="Rarefied richness 18S",
     pch=md.18S.p$col_plot, ylim=range(0:2000), xlim=range(0:700))
abline(SSUCOIRS,col="darkgrey")
# legend("topleft",ncol=3,legend=sort(unique(md.COI.p$Platform)),
#        #col=sort(unique(md.COI.p$col_plot)),
#        pch=sort(unique(md.COI.p$col_plot)))#,cex=.8)
text(400, 40, "R² = .54, p < 1E-15", pos=4,  col="#606060")

