# Andrea Bagi 2021

## This script that was used to:
##  1. generate input tables for CoNet analysis (NB: post-processing was done partly in Excel)
##  2. process node and edge tables


## The script is divided into two parts:
##  1. 18S data processing
##  2. COI data processing

# Change to correct path for local directory (github clone) and uncomment:
# setwd("/home/alanzen/projects/Metamon2/")


## Load necessary libraries
library(dplyr)
library(VennDiagram)
library(gplots)
library(ggplot2)
library(ggpubr)
library(hrbrthemes)
library(vegan)



#########################
#### PART-1 18S DATA ####
#########################



#*************************************#
### GENERATE INPUT TABLES FOR CONET ###
#*************************************#


### IMPORT 18S DATA  ###
md.18S = read.table("co-occurence_networks/Metadata_all_reps_w_PIs_from_pooled.csv",
                    row.names=1,header=T,sep=",")
otus.18S = read.table("SWARM_WP1_20200826_18S/CREST_NoPlankton/All_Assignments.tsv",
                      header=T,sep="\t")
## Should contain 963 18S taxa in 296 samples with 28 metadata variables


### EXTRACT TAXONOMY INFORMATION & CROP OTU TABLE###
taxpath.18S = otus.18S[,1:3] ## moving assignments into a new data frame

otus.18S = otus.18S[,3:length(colnames(otus.18S))]
rownames(otus.18S) = otus.18S[,1]
otus.18S = otus.18S[,2:length(colnames(otus.18S))] ## Taxa are rows

otus.18S.t = as.data.frame(t(otus.18S)) ## Samples are rows
table(row.names(md.18S) == row.names(otus.18S.t)) ## check if TRUE for all 296 samples


### SPLIT THE OTU TABLES - IMPACTED & NON-IMPACTED ###
outpath = "co-occurence_networks/Manuscript/analysis_and_scripts/SplitNetwork/18S_All_assignments/commonStations/"

## Select only region II and III samples 
otus.sel = otus.18S.t[md.18S$Region=="II"|md.18S$Region=="III",]  ## Samples are rows
md.sel = md.18S[md.18S$Region=="II"|md.18S$Region=="III",]
table(row.names(md.sel) == row.names(otus.sel))


### CREATE IMPACTED OTU table and TAXONOMY file ###
otus.18S.imp = otus.sel[md.sel$pi>2|md.sel$pi==2,]  ## Select contaminated samples based on PI (samples are rows)
md.18S.imp = md.sel[md.sel$pi>2|md.sel$pi==2,]
table(row.names(md.18S.imp) == row.names(otus.18S.imp)) ## check if TRUE for 38 samples
summary(md.18S.imp$pi) ### checking if PI is indeed > 2

## Taxon filtering - remove taxa that do not occur in any of the samples
imp.taxa.inc=otus.18S.imp ## create a temporary data frame where samples are rows
imp.taxa.inc[imp.taxa.inc>0]=1 ## convert to presence absence table 
summary(colSums(imp.taxa.inc)) ## check if there are any taxa to be removed (min = 0) and that max = 38
colSums=colSums(imp.taxa.inc)
indices.prev=which(colSums>=1) ## colect the indices of taxa to be removed
otus.f <- as.data.frame(otus.18S.imp[,indices.prev]) ## remove taxa - samples are rows
otus.f.t <- as.data.frame(t(otus.f)) ## transpose the table for CoNet - taxa are rows
write.table(otus.f.t, file=file.path(outpath,"split_18S_separateRep_impacted_abundances.txt"), quote = F, sep = "\t") 
## Need to be formatted for CoNet - two taxa from Cercozoa: 7-2.3 and 7-5.4 names need to be changed Cercozoa-07.02.2003 and Cercozoa-07.05.2004


## Use the taxa list to generate Venn diagram on input taxa
venn.18S.1 = rownames(otus.f.t) ## impacted


## matching the taxonomy table with the remaining taxa 
keep <- as.data.frame(rownames(otus.f.t))
colnames(keep) <- "Taxon"
tax = taxpath %>% right_join(keep, by = "Taxon")
rownames(tax) <- tax$Taxon
table(row.names(tax) == row.names(otus.f.t))
write.table(tax, file=file.path(outpath,"18S_impacted_taxonomy.txt"), quote = F, sep = "\t")
## Need to be formatted for CoNet - two taxa from Cercozoa: 7-2.3 and 7-5.4 names need to be changed Cercozoa-07.02.2003 and Cercozoa-07.05.2004


## Clean-up
rm(indices.prev,colSums,imp.taxa.inc, keep, otus.f, otus.f.t, tax, otus.sel, md.sel)



### CREATE NON-IMPACTED OTU table and TAXONOMY file ###
## Restrict to same platforms as contaminated and then select based on PI
otus.sel = otus.18S.t[md.18S$Platform=="OSS"|md.18S$Platform=="RIN"|md.18S$Platform=="OSF"|md.18S$Platform=="VFR",]  ## Samples are rows
md.sel = md.18S[md.18S$Platform=="OSS"|md.18S$Platform=="RIN"|md.18S$Platform=="OSF"|md.18S$Platform=="VFR",]
table(row.names(md.sel) == row.names(otus.sel)) ## check if TRUE for 90 samples


## Selecting samples based on PI
otus.18S.nonimp = otus.sel[md.sel$pi<2,]  ## Samples are rows
md.18S.nonimp = md.sel[md.sel$pi<2,]
table(row.names(md.18S.nonimp) == row.names(otus.18S.nonimp)) ## check if TRUE for 52 sample
summary(md.18S.nonimp$pi)  ### check if PI is indeed < 2


## Taxon filtering - remove taxa that do not occur in any of the samples - same procedure as above
imp.taxa.inc=otus.18S.nonimp
imp.taxa.inc[imp.taxa.inc>0]=1
summary(colSums(imp.taxa.inc))
colSums=colSums(imp.taxa.inc) ## MAX is 52 = number of sample -> OK
indices.prev=which(colSums>=1)
otus.f <- as.data.frame(otus.18S.nonimp[,indices.prev]) ## Samples are rows
otus.f.t <- as.data.frame(t(otus.f)) ## Taxa are rows
write.table(otus.f.t, file=file.path(outpath,"split_18S_separateRep_nonImpacted_abundances.txt"), 
            quote = F, sep = "\t") ## Rows are taxa -  numbers
## Need to be formatted for CoNet - two taxa from Cercozoa: 7-2.3 and 7-5.4 names need to be changed Cercozoa-07.02.2003 and Cercozoa-07.05.2004


## Use the taxa list to generate Venn diagram on input taxa
venn.18S.2 = rownames(otus.f.t) ## non-impacted taxa


### matching taxonomy table with remaining taxa list
keep <- as.data.frame(rownames(otus.f.t))
colnames(keep) <- "Taxon"
tax = taxpath.18S %>% right_join(keep, by = "Taxon")
rownames(tax) <- tax$Taxon
table(row.names(tax) == row.names(otus.f.t)) ## TRUE for all 546 taxa
write.table(tax, file=file.path(outpath,"split_18S_separateRep_nonImpacted_taxonomy.txt"), quote = F, sep = "\t")
## Need to be formatted for CoNett - two taxa from Cercozoa: 7-2.3 and 7-5.4 names need to be changed Cercozoa-07.02.2003 and Cercozoa-07.05.2004


## Clean-up
rm(indices.prev,colSums,imp.taxa.inc, otus.sel, md.sel,keep, otus.f, otus.f.t, tax, otus.18S.t)

### END OF - GENERATE INPUT TABLES FOR CONET ###
#**********************************************#




#*******************#
### VENN DIAGRAMS ###
#*******************#


## Purpose: generating NODE and EDGE categories (e.g. unique vs common, bioindicator vs non-bioindicator) 
## to be imported back to Cytoscape


##-------------------------------------##
## Starting with input taxa categories ##
##-------------------------------------##

VENN.LIST <- list(venn.18S.1, venn.18S.2) ## lists were generated through PART-1 steps
venn.plot <- venn.diagram(VENN.LIST , NULL, category.names=c("Impacted","Non-Impacted"))
grid.draw(venn.plot) 

## Collect taxa names of interest and create a table to be merged with node table later / add to Cytoscape as node properties
obj <- venn(VENN.LIST, show.plot=FALSE)
inters <- attr(obj,"intersections")
lapply(inters, head) ## checking

## split input taxa into 3 categories and add columns
impTaxa = as.data.frame(inters$A)
colnames(impTaxa) = "shared.name"
impTaxa$venn = rep("uniqueImpacted", length(rownames(impTaxa)))

nonimpTaxa = as.data.frame(inters$B)
colnames(nonimpTaxa) = "shared.name"
nonimpTaxa$venn = rep("uniqueNonImpacted", length(rownames(nonimpTaxa)))

commonTaxa = as.data.frame(inters$'A:B')
colnames(commonTaxa) = "shared.name"
commonTaxa$venn = rep("common", length(rownames(commonTaxa)))

## combine into a single data frame to be used later and correct the two taxon names
taxa.df = rbind(impTaxa, nonimpTaxa, commonTaxa)
taxa.df$shared.name <- gsub("\\ ","\\-",taxa.df$shared.name) ## format taxa names to match CoNet output
taxa.df[233, 1] = "Cercozoa-07.02.2003" ## rename as done for CoNet input tables
taxa.df[238, 1] = "Cercozoa-07.05.2004" ## rename as done for CoNet input tables



##--------------------------------------------------##
## Next work with the node tables and bioindicators ##
##--------------------------------------------------##

## Import data
bioInd = read.table("co-occurence_networks/Taxa_IndGroup_18S_PI_Cytoscape_dash and Cercozoa_corrected.txt",row.names=1,header=T,sep="\t")
bioInd = filter(bioInd, bioInd$pi > 0) ## select only those assigned EG I-V
bioInd$shared.name = rownames(bioInd) ## Bioindicators- 18S

impNode = read.table("co-occurence_networks/18S_impacted_network_node_table.csv",header=T,sep=",")
rownames(impNode)=impNode$name ## Nodes from the impacted network

nonimpNode = read.table("co-occurence_networks/18S_nonimpacted_network_node_table.csv",header=T,sep=",")
rownames(nonimpNode)=nonimpNode$name ## Nodes from the non-Impacted network

## define the lists for Venn diagram
set1 = rownames(bioInd)
set2 = rownames(impNode)
set3 = rownames(nonimpNode)

VENN.LIST <- list(set1, set2, set3)
venn.plot <- venn.diagram(VENN.LIST , NULL, 
                          filename = "venn_18S_bioindicators_nodes.tiff", height = 4000, width = 6000, 
                          resolution = 500, imagetype = "tiff", units = "px", compression = "lzw", na = "stop",
                          category.names=c("Bioindicators", "Nodes in impacted network", "Nodes in non-impacted network"),
                          fill = c("#E69F00", "#56B4E9", "#009E73"),
                          lty = 'blank', lwd = rep(2, 3),
                          cex = 1.5,
                          fontface = "italic",
                          cat.cex = 1.5,
                          cat.fontface = "bold", cat.default.pos = "outer",
                          cat.dist = c(0.075, 0.075, 0.05), cat.pos = c(335, 25, 180),
                          margin = 0.2)

## Collect taxa names of interest and create a table to be merged with node table later 
## and add to Cytoscape as node properties
obj <- venn(VENN.LIST, show.plot=FALSE)
inters <- attr(obj,"intersections")
lapply(inters, head) ## checking


### split node taxa into 4 categories based on bioindicator status & network origin ###
x = as.data.frame(inters$B)
y = as.data.frame(inters$`B:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
impNonBioInd <- rbind(x, y) ## Impacted nodes that are not bioindicators = 105 -> OK

x = as.data.frame(inters$C)
y = as.data.frame(inters$`B:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
nonimpNonBioInd <- rbind(x, y) ## Non-impacted nodes that are not bioindicators = 72 taxa -> OK

x = as.data.frame(inters$`A:B`)
y = as.data.frame(inters$`A:B:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
impBioInd <- rbind(x, y) ## Impacted nodes that are bioindicators = 63 -> OK

x = as.data.frame(inters$`A:C`)
y = as.data.frame(inters$`A:B:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
nonimpBioInd <- rbind(x, y) ## Non-impacted nodes that are bioindicators = 43 -> OK



### split node taxa into 3 categories based on which network they are present in ###
x = as.data.frame(inters$C)
y = as.data.frame(inters$`A:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
nonimpUnique <- rbind(x, y) ## Non-impacted nodes unique to this network = 22 -> OK
nonimpUnique$node.venn = rep("uniqueNonimpactedNode", length(rownames(nonimpUnique)))

x = as.data.frame(inters$B)
y = as.data.frame(inters$`A:B`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
impUnique <- rbind(x, y) ## Impacted nodes unique to this network = 75 -> OK
impUnique$node.venn = rep("uniqueImpactedNode", length(rownames(impUnique)))

x = as.data.frame(inters$`B:C`)
y = as.data.frame(inters$`A:B:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
common <- rbind(x, y) ## Nodes common to both networks = 93 -> OK
common$node.venn = rep("commonNode", length(rownames(common)))


## combine into a single data frame 
## all nodes assigned a 'node.venn' category indicating which network they are present in
node.df = rbind(nonimpUnique, impUnique, common) ## Should contain 190 taxa -> OK


## clean-up
rm(inters, x, y, VENN.LIST, venn.plot, obj, set1, set2, set3)


## Continue working with the bioindicator lists divided into impacted/non-impacted
## Add node attributes category and indicator status 
## Category indicates the network origin of the node
impNonBioInd$indicator = rep("no", length(impNonBioInd))
impNonBioInd$category = rep("impacted", length(rownames(impNonBioInd)))
colnames(impNonBioInd) = c("shared.name", "indicator", "category")

impBioInd$indicator = rep("yes", length(rownames(impBioInd)))
impBioInd$category = rep("impacted", length(rownames(impBioInd)))
colnames(impBioInd) = c("shared.name", "indicator", "category")

nonimpNonBioInd$indicator = rep("no", length(nonimpNonBioInd))
nonimpNonBioInd$category = rep("nonimpacted", length(rownames(nonimpNonBioInd)))
colnames(nonimpNonBioInd) = c("shared.name", "indicator", "category")

nonimpBioInd$indicator = rep("yes", length(nonimpBioInd))
nonimpBioInd$category = rep("nonimpacted", length(rownames(nonimpBioInd)))
colnames(nonimpBioInd) = c("shared.name", "indicator", "category")


### Compile all this into one data frame containing node properties from Cytoscape as well ###

## combined the two small tables first and than add into the node data
imp.data = rbind(impBioInd, impNonBioInd) ## 168 nodes -> OK
non.imp.data = rbind(nonimpBioInd, nonimpNonBioInd) ## 115 nodes -> OK

## add node table information
imp.data = imp.data %>% left_join(impNode, by = "shared.name")
non.imp.data = non.imp.data %>% left_join(nonimpNode, by = "shared.name")

## order the columns alphabetically
imp.data = imp.data[,order(colnames(imp.data))]
non.imp.data = non.imp.data[,order(colnames(non.imp.data))]

## select columns that match in the two tables
cols_to_keep <- intersect(colnames(imp.data),colnames(non.imp.data))
imp.data <- imp.data[,cols_to_keep, drop=FALSE]
non.imp.data <- non.imp.data[,cols_to_keep, drop=FALSE]
table(colnames(imp.data) == colnames(non.imp.data)) ## true for all 36 columns

## Create a combined table with all
df = rbind(imp.data, non.imp.data)
colnames(df)
cols_to_keep = c("shared.name", "category", "indicator", "abundance", "AverageShortestPathLength", 
                 "BetweennessCentrality", "ClosenessCentrality", "ClusteringCoefficient", "NeighborhoodConnectivity",
                 "degree", "negdegree", "posdegree",
                 "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
df <- df[,cols_to_keep, drop=FALSE]

df = df %>% left_join(taxa.df, by = "shared.name") ## add info whether taxa were present in original abundance tables
df = df %>% left_join(node.df, by = "shared.name") ## add info whether taxa were unique or shared nodes



#### COMPARATIVE NODE BOXPLOTS ####


## Some examples
## Compare the properties of nodes in the impacted vs non-impacted network
p = ggboxplot(df, x = "category", y = "NeighborhoodConnectivity",
              color = "category", palette = "jco")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Environmental status", y = "Neighborhood Connectivity")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")


## Compare biondicator vs non-bioindicator nodes
## within each network => facet by network origin of the node
p = ggboxplot(df, x = "indicator", y = "NeighborhoodConnectivity",
              color = "indicator", palette = "jco", facet.by = "category")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Biondicator", y = "Neighborhood Connectivity")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")


## Do the reverse - compare bioindicators vs non-bioindicators based on their 
## network origin => facet by indicator status 
p = ggboxplot(df, x = "category", y = "negdegree",
              color = "category", palette = "jco", facet.by = "indicator")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Environmental status", y = "Negative node degree")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")


## Compare properties of node that are shared
## keep only shared nodes
df.f = filter(df, df$node.venn == "commonNode")

p = ggboxplot(df.f, x = "category", y = "NeighborhoodConnectivity",
              color = "category", palette = "jco")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Network property from", y = "Neighborhood Connectivity")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")



## clean-up at the END of NODE Venns
rm(df, df.f, taxa.df, node.df, bioInd, common, 
   imp.data, impBioInd, impNode, impNonBioInd, impUnique,
   non.imp.data, nonimpBioInd, nonimpNode, nonimpNonBioInd, nonimpUnique,
   p, cols_to_keep)



##---------------------------##
## Work with the edge tables ##
##---------------------------##

## Import edge tables
impEdge = read.table("co-occurence_networks/18S_impacted_network_edge_table.csv",header=T,sep=",")
rownames(impEdge)=impEdge$Label ## Impacted network
vennSet1 = impEdge$Label

nonimpEdge = read.table("co-occurence_networks/18S_nonimpacted_network_edge_table.csv",header=T,sep=",")
rownames(nonimpEdge)=nonimpEdge$Label ## Non-impacted network
vennSet2 = nonimpEdge$Label

## Edge Venn
VENN.LIST <- list(vennSet1, vennSet2)
venn.plot <- venn.diagram(VENN.LIST , NULL, category.names=c("Edges in impacted network", "Edges in non-impacted network"))

venn.plot <- venn.diagram(VENN.LIST , NULL, 
                          filename = "venn_18S_edges_#2.tiff", height = 4000, width = 6000, 
                          resolution = 500, imagetype = "tiff", units = "px", compression = "lzw", na = "stop",
                          category.names=c("Impacted", "Non-impacted"),
                          fill = c("#E69F00", "#56B4E9"),
                          lty = 'blank', lwd = rep(2, 2),
                          cex = 2,
                          fontface = "italic",
                          cat.cex = 2,
                          cat.fontface = "bold", cat.default.pos = "outer",
                          cat.dist = c(0.05, 0.04), cat.pos = c(345, 10),
                          margin = 0.2)


## extract information about network origin
obj <- venn(VENN.LIST, show.plot=FALSE)
inters <- attr(obj,"intersections")
lapply(inters, head)

nonimpUnique = as.data.frame(inters$B) ## 226 -> OK
colnames(nonimpUnique) = "Label"
nonimpUnique$edge.venn = rep("uniqueNonimpactedEdge", length(rownames(nonimpUnique)))

impUnique = as.data.frame(inters$A) ## 486 -> OK
colnames(impUnique) = "Label"
impUnique$edge.venn = rep("uniqueImpactedEdge", length(rownames(impUnique)))

common = as.data.frame(inters$'A:B') ## 96 -> OK
colnames(common) = "Label"
common$edge.venn = rep("commonEdge", length(rownames(common)))

edge.df = rbind(nonimpUnique, impUnique, common)
edge.df.imp = impEdge %>% left_join(edge.df, by = "Label")
edge.df.nonimp = nonimpEdge %>% left_join(edge.df, by = "Label")

write.table(edge.df.imp, "18S_edge_table_impacted_combined.txt", sep = "\t", quote = F)
write.table(edge.df.nonimp, "18S_edge_table_non-impacted_combined.txt", sep = "\t", quote = F)


## Harmonize add category and combine
impEdge$category = rep("impacted",length(rownames(impEdge)))
nonimpEdge$category = rep("nonimpacted",length(rownames(nonimpEdge)))

impEdge = impEdge[,order(colnames(impEdge))]
nonimpEdge = nonimpEdge[,order(colnames(nonimpEdge))]

df.edge = rbind(impEdge, nonimpEdge)
colnames(df.edge)

df.edge = df.edge %>% left_join(edge.df, by = "Label") ## add info whether taxa were unique or shared nodes
## write.table(df.edge, "18S_edge_tables_combined.txt", sep = "\t", quote = F)



#### COMPARATIVE EDGE BOXPLOTS ####



#### Overall comparison between impacted and non-impacted ####
p = ggboxplot(df.edge, x = "category", y = "qval",
              color = "category", palette = "jco")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Environmental status", y = "Adjusted p-values")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")


summary(impEdge$qval)
## Min.       1st Qu.    Median      Mean       3rd Qu.      Max. 
## 0.000e+00  6.000e-08  3.306e-04   7.936e-03  1.078e-02    4.969e-02 
summary(nonimpEdge$qval)
## Min.       1st Qu.    Median      Mean       3rd Qu.      Max. 
## 0.000e+00  2.870e-06  7.817e-04   7.705e-03  8.884e-03    4.879e-02 


### Compare only common edges ### 
df.edge.f = filter(df.edge, df.edge$edge.venn == "commonEdge")
## Comparison between impacted and non-impacted among the common edges ####
p = ggboxplot(df.edge.f, x = "category", y = "method_number",
              color = "category", palette = "jco")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Environmental status", y = "Method number of common edges")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")


## Clean-up
rm(common, df.edge, df.edge.f,edge.df, edge.df.imp, edge.df.nonimp,
   impEdge, impUnique, inters, nonimpEdge, nonimpUnique,
   VENN.LIST, venn.plot, obj, vennSet1, vennSet2)


#### END of 18S ANALYSIS
#********************************************************************************#






#########################
#### PART-2 COI DATA ####
#########################



#******************************#
### GENERATE INPUT FOR CONET ###
#******************************#


outpath = "co-occurence_networks/Manuscript/analysis_and_scripts/SplitNetwork/COI_All_assignments"


### IMPORT COI data ###
md.COI = read.table("co-occurence_networks/samplesCategorized.txt",header = T, sep = "\t")
otus.COI = read.table("SWARM_WP1_20200722_COI/CREST_NoPlankton/All_Assignments.tsv",
                      header=T,sep="\t")


## clean and split OTUs
taxpath.COI = otus.COI[,1:3] ## 774 taxa in COI data
otus.COI = otus.COI[,3:length(colnames(otus.COI))]
rownames(otus.COI) = otus.COI[,1]
otus.COI = otus.COI[,2:length(colnames(otus.COI))] ## Taxa are rows
otus.COI.t = as.data.frame(t(otus.COI)) ## Samples are rows

summary(otus.COI.t$root)
summary(otus.COI.t$`Cellular organisms`)
summary(otus.COI.t$Eukaryota)
summary(otus.COI.t$Opisthokonta)
summary(otus.COI.t$Metazoa)
summary(otus.COI.t$Arthropoda)

otus.COI = otus.COI[6:length(rownames(otus.COI)),] ## Taxa are rows
otus.COI.t = as.data.frame(t(otus.COI)) ## Samples are rows
taxpath.COI = taxpath.COI[6:length(rownames(taxpath.COI)),] ## 774 taxa in COI data


## Formatting and splitting samples based on sample classification in metadata
otus.COI.t$Sample = rownames(otus.COI.t)
otus.COI.all = otus.COI.t %>% semi_join(md.COI, by = "Sample")
md.COI.all = md.COI %>% semi_join(otus.COI.all,by = "Sample")

otus.COI.imp = otus.COI.all[md.COI.all$Category=="impacted",] ## 33 samples
otus.COI.nonimp = otus.COI.all[md.COI.all$Category=="nonimpacted",] ## 49 samples

otus.COI.imp = subset(otus.COI.imp, select=-c(Sample))
otus.COI.nonimp = subset(otus.COI.nonimp, select=-c(Sample))


## Checking if there are taxa with all zero abundances
summary(colSums(otus.COI.imp))
summary(colSums(otus.COI.nonimp))




#### Create impacted COI tables ####
### Taxon filtering - IMPACTED -remove taxa that do not occur in any of the samples
imp.taxa.inc=otus.COI.imp ## Samples are rows
imp.taxa.inc[imp.taxa.inc>0]=1
summary(colSums(imp.taxa.inc)) ## MAX = 33 -> OK
colSums=colSums(imp.taxa.inc)
indices.prev=which(colSums>=1)

otus.f <- as.data.frame(otus.COI.imp[,indices.prev]) ## Samples are rows - Taxa are columns
otus.f.t <- as.data.frame(t(otus.f)) ## taxa are rows - Transposed read count table
## write.table(otus.f.t, file=file.path(outpath,"COI_impacted_abundances.txt"), quote = F, sep = "\t") 
## Rows are taxa - read numbers
## Need to be formatted for CoNet - Rememeber 


venn.COI.taxa.1 = rownames(otus.f.t) ## impacted


### FILTER taxonomy - IMPACTED
keep <- as.data.frame(colnames(otus.f))
colnames(keep) <- "Taxon"
tax = taxpath.COI %>% right_join(keep, by = "Taxon")
rownames(tax) <- tax$Taxon
table(row.names(tax) == colnames(otus.f)) ## TRUE for 169 -> OK

## write.table(tax, file=file.path(outpath,"COI_impacted_taxonomy.txt"), quote = F, sep = "\t")
rm(imp.taxa.inc, tax, indices.prev, keep, colSums)



#### Create non-impacted COI tables ####
### Taxon filtering - NON-IMPACTED - remove taxa that do not occur in any of the samples
imp.taxa.inc=otus.COI.nonimp ## Samples are rows
imp.taxa.inc[imp.taxa.inc>0]=1
summary(colSums(imp.taxa.inc)) ## MAX = 49 -> OK
colSums=colSums(imp.taxa.inc)
indices.prev=which(colSums>=1)

otus.f <- as.data.frame(otus.COI.nonimp[,indices.prev]) ## Samples are rows
otus.f.t <- as.data.frame(t(otus.f)) ## taxa are rows
## write.table(otus.f.t, file=file.path(outpath,"COI_nonimpacted_relative_abundances.txt"), quote = F, sep = "\t") 
## Rows are taxa - read numbers
## Need to be formatted for CoNet - Rememeber 


venn.COI.taxa.2 = rownames(otus.f.t) ## non-impacted


### FILTER taxonomy - impacted and nonImpacted are the same and also global ...
keep <- as.data.frame(colnames(otus.f))
colnames(keep) <- "Taxon"
tax = taxpath.COI %>% right_join(keep, by = "Taxon")
rownames(tax) <- tax$Taxon
table(row.names(tax) == colnames(otus.f)) ## TRYE for 236 -> OK

## write.table(tax, file=file.path(outpath,"COI_non-impacted_taxonomy.txt"), quote = F, sep = "\t")
rm(indices.prev,colSums,imp.taxa.inc, keep, otus.f, otus.f.ra, otus.f.ra.t, otus.f.t, tax)



#***********************#
### VENN DIAGRAMs COI ###
#***********************#



VENN.LIST <- list(venn.COI.taxa.1, venn.COI.taxa.2)
venn.plot <- venn.diagram(VENN.LIST , NULL, category.names=c("Impacted","Non-Impacted"))
grid.draw(venn.plot) 


# in order to collect gene names of interest
obj <- venn(VENN.LIST, show.plot=FALSE)
inters <- attr(obj,"intersections")
lapply(inters, head)


## Collect taxa from the input tables and categorize based on their presence in each or both
## split taxa into 3 categories and add columns
impTaxa = as.data.frame(inters$A)
colnames(impTaxa) = "shared.name"
impTaxa$venn = rep("uniqueImpacted", length(rownames(impTaxa))) ## 27 -> OK

nonimpTaxa = as.data.frame(inters$B)
colnames(nonimpTaxa) = "shared.name"
nonimpTaxa$venn = rep("uniqueNonImpacted", length(rownames(nonimpTaxa))) ## 94 -> OK

commonTaxa = as.data.frame(inters$'A:B')
colnames(commonTaxa) = "shared.name"
commonTaxa$venn = rep("common", length(rownames(commonTaxa))) ## 142 -> OK


## Combine taxa info 
taxa.df = rbind(impTaxa, nonimpTaxa, commonTaxa) ## 263 + 142 = 405 -> OK
taxa.df$shared.name <- gsub("\\ ","\\-",taxa.df$shared.name) ## Replace space with dashed - due to CoNet

## Clean-up
rm(obj,VENN.LIST,venn.plot)



#------------------------------#
# Working with the node tables #
#------------------------------#

### Working with network analysis tables and bioindicators
bioInd = read.table("co-occurence_networks/Taxa_IndGroup_COI_PI.csv",row.names=1,header=T,sep=",")
bioInd = filter(bioInd, bioInd$pi > 0)
bioInd$shared.name = rownames(bioInd)
bioInd$shared.name <- gsub("\\ ","\\-",bioInd$shared.name) ## Names need to be "dashed"
rownames(bioInd) = bioInd$shared.name

## Impacted network split
impNode = read.table("co-occurence_networks/COI_impacted_network_node_table.csv",header=T,sep=",")
rownames(impNode)=impNode$name
## write.table(impNode, file=file.path(outpath,"COI_impacted_node_table.txt"), quote = F, sep = "\t")

## Non-Impacted network split
nonimpNode = read.table("co-occurence_networks/COI_nonimpacted_network_node_table.csv",header=T,sep=",")
rownames(nonimpNode)=nonimpNode$name
## write.table(nonimpNode, file=file.path(outpath,"COI_non-impacted_node_table.txt"), quote = F, sep = "\t")


## Repeat the same process as for 18S above
## define the taxa lists present in the three tables
set1 = rownames(bioInd)
set2 = rownames(impNode)
set3 = rownames(nonimpNode)

VENN.LIST <- list(set1, set2, set3)
venn.plot <- venn.diagram(VENN.LIST , NULL, 
            filename = "venn_COI_bioindicators_nodes.tiff", height = 4000, width = 6000, 
            resolution = 500, imagetype = "tiff", units = "px", compression = "lzw", na = "stop",
            category.names=c("Bioindicators", "Nodes in impacted network", "Nodes in non-impacted network"),
            fill = c("#E69F00", "#56B4E9", "#009E73"),
            lty = 'blank', lwd = rep(2, 3),
            cex = 1.5,
            fontface = "italic",
            cat.cex = 1.5,
            cat.fontface = "bold", cat.default.pos = "outer",
            cat.dist = c(0.075, 0.075, 0.05), cat.pos = c(335, 25, 180),
            margin = 0.2
            )


# in order to collect gene names of interest
obj <- venn(VENN.LIST, show.plot=FALSE)
inters <- attr(obj,"intersections")
lapply(inters, head)


### Divide taxa into 3 categories based on which network they are present in 
x = as.data.frame(inters$C)
y = as.data.frame(inters$`A:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
nonimpUnique <- rbind(x, y) ## Non-impacted nodes unique to this network = 12 -> OK
nonimpUnique$node.venn = rep("uniqueNonimpactedNode", length(rownames(nonimpUnique)))

x = as.data.frame(inters$B)
y = as.data.frame(inters$`A:B`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
impUnique <- rbind(x, y) ## Impacted nodes unique to this network = 8 -> OK
impUnique$node.venn = rep("uniqueImpactedNode", length(rownames(impUnique)))

x = as.data.frame(inters$`B:C`)
y = as.data.frame(inters$`A:B:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
common <- rbind(x, y) ## Nodes common to both networks = 36 -> OK
common$node.venn = rep("commonNode", length(rownames(common)))

## Combine node table
node.df = rbind(nonimpUnique, impUnique, common) ## Should contain 56 taxa -> OK



### split node taxa into 4 categories based on bioindicator status & network origin ###
x = as.data.frame(inters$B)
y = as.data.frame(inters$`B:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
impNonBioInd <- rbind(x, y) ## Impacted nodes that are not bioindicators = 32 -> OK

x = as.data.frame(inters$C)
y = as.data.frame(inters$`B:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
nonimpNonBioInd <- rbind(x, y) ## Non-impacted nodes that are not bioindicators = 32 taxa -> OK

x = as.data.frame(inters$`A:B`)
y = as.data.frame(inters$`A:B:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
impBioInd <- rbind(x, y) ## Impacted nodes that are bioindicators = 12 -> OK

x = as.data.frame(inters$`A:C`)
y = as.data.frame(inters$`A:B:C`)
colnames(x) = "shared.name"
colnames(y) = "shared.name"
nonimpBioInd <- rbind(x, y) ## Non-impacted nodes that are bioindicators = 16 -> OK


## Add their attributes (category and indicator status)
impNonBioInd$indicator = rep("no", length(impNonBioInd))
impNonBioInd$category = rep("impacted", length(rownames(impNonBioInd)))
colnames(impNonBioInd) = c("shared.name", "indicator", "category")

impBioInd$indicator = rep("yes", length(rownames(impBioInd)))
impBioInd$category = rep("impacted", length(rownames(impBioInd)))
colnames(impBioInd) = c("shared.name", "indicator", "category")

nonimpNonBioInd$indicator = rep("no", length(nonimpNonBioInd))
nonimpNonBioInd$category = rep("nonimpacted", length(rownames(nonimpNonBioInd)))
colnames(nonimpNonBioInd) = c("shared.name", "indicator", "category")

nonimpBioInd$indicator = rep("yes", length(nonimpBioInd))
nonimpBioInd$category = rep("nonimpacted", length(rownames(nonimpBioInd)))
colnames(nonimpBioInd) = c("shared.name", "indicator", "category")


### Bring it altogether ###
## Combine the two small tables first and than add into the node data
imp.data = rbind(impBioInd, impNonBioInd) ## 44 nodes -> OK
non.imp.data = rbind(nonimpBioInd, nonimpNonBioInd) ## 48 nodes -> OK

## Add node table information
imp.data = imp.data %>% left_join(impNode, by = "shared.name")
non.imp.data = non.imp.data %>% left_join(nonimpNode, by = "shared.name")

## Order the columns alphabetically
imp.data = imp.data[,order(colnames(imp.data))]
non.imp.data = non.imp.data[,order(colnames(non.imp.data))]

## Select columns of interest - those that match in the two tables
cols_to_keep <- intersect(colnames(imp.data),colnames(non.imp.data))
imp.data <- imp.data[,cols_to_keep, drop=FALSE]
non.imp.data <- non.imp.data[,cols_to_keep, drop=FALSE]
table(colnames(imp.data) == colnames(non.imp.data)) ## true for all 34 columns

## Create a combined table with all
df = rbind(imp.data, non.imp.data)
colnames(df)
cols_to_keep = c("shared.name", "category", "indicator", "abundance", "AverageShortestPathLength", 
                 "BetweennessCentrality", "ClosenessCentrality", "ClusteringCoefficient", "NeighborhoodConnectivity",
                 "degree", "negdegree", "posdegree",
                 "phylum", "class", "order", "family", "genus", "species")
df <- df[,cols_to_keep, drop=FALSE]

df = df %>% left_join(taxa.df, by = "shared.name") ## add info whether taxa were present in original abundance tables
df = df %>% left_join(node.df, by = "shared.name") ## add info whether taxa were unique or shared nodes

colnames(df)

## setwd("co-occurence_networks/Manuscript/results")
## write.table(df, file = "COI_all_combined.txt", sep = "\t", quote = F)




#### COMPARATIVE NODE BOXPLOTS ####



## Overall comparison between impacted and non-impacted nodes
p = ggboxplot(df, x = "category", y = "NeighborhoodConnectivity",
              color = "category", palette = "jco")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Environmental status", y = "Neighborhood Connectivity")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")



## Indicators vs non-indicators within each network
p = ggboxplot(df, x = "indicator", y = "NeighborhoodConnectivity",
              color = "indicator", palette = "jco", facet.by = "category")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Biondicator", y = "Neighborhood Connectivity")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")


#### Indicators vs non-indicators between the two network categories ####
p = ggboxplot(df, x = "category", y = "degree",
              color = "category", palette = "jco", facet.by = "indicator")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Environmental status", y = "Node degree")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")


#### Comparing node properties  of common nodes
## keep only common nodes in the table
df.f = filter(df, df$node.venn == "commonNode")

p = ggboxplot(df.f, x = "category", y = "NeighborhoodConnectivity",
              color = "category", palette = "jco")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Network property from", y = "Neighborhood Connectivity")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")




#------------------------------#
# Working with the edge tables #
#------------------------------#

## Impacted network split
impEdge = read.table("co-occurence_networks/COI_impacted_network_edge_table.csv",header=T,sep=",")
rownames(impEdge)=impEdge$Label
vennSet1 = impEdge$Label

## Non-Impacted network split
nonimpEdge = read.table("co-occurence_networks/COI_nonimpacted_network_edge_table.csv",header=T,sep=",")
rownames(nonimpEdge)=nonimpEdge$Label
vennSet2 = nonimpEdge$Label

## MAKE A NICE PLOT
VENN.LIST <- list(vennSet1, vennSet2)
## third color: "#009E73"
venn.plot <- venn.diagram(VENN.LIST , NULL, 
                          filename = "venn_COI_edges.tiff", height = 4000, width = 6000, 
                          resolution = 500, imagetype = "tiff", units = "px", compression = "lzw", na = "stop",
                          category.names=c("Impacted", "Non-impacted"),
                          fill = c("#E69F00", "#56B4E9"),
                          lty = 'blank', lwd = rep(2, 2),
                          cex = 2,
                          fontface = "italic",
                          cat.cex = 2,
                          cat.fontface = "bold", cat.default.pos = "outer",
                          cat.dist = c(0.05, 0.04), cat.pos = c(345, 10),
                          margin = 0.2)

## MAKE A REGULAR PLOT and EXPLORE ##
venn.plot <- venn.diagram(VENN.LIST , NULL, category.names=c("Edges in impacted network", "Edges in non-impacted network"))
grid.draw(venn.plot)

# in order to collect gene names of interest
obj <- venn(VENN.LIST, show.plot=FALSE)
inters <- attr(obj,"intersections")
lapply(inters, head)


## SORT BY NETWORK
nonimpUnique = as.data.frame(inters$B) ## 226 -> OK
colnames(nonimpUnique) = "Label"
nonimpUnique$edge.venn = rep("uniqueNonimpactedEdge", length(rownames(nonimpUnique)))

impUnique = as.data.frame(inters$A) ## 486 -> OK
colnames(impUnique) = "Label"
impUnique$edge.venn = rep("uniqueImpactedEdge", length(rownames(impUnique)))

common = as.data.frame(inters$'A:B') ## 96 -> OK
colnames(common) = "Label"
common$edge.venn = rep("commonEdge", length(rownames(common)))

edge.df = rbind(nonimpUnique, impUnique, common)
edge.df.imp = impEdge %>% left_join(edge.df, by = "Label")
edge.df.nonimp = nonimpEdge %>% left_join(edge.df, by = "Label")

## write.table(edge.df.imp, "COI_edge_table_impacted_combined.txt", sep = "\t", quote = F)
## write.table(edge.df.nonimp, "COI_edge_table_non-impacted_combined.txt", sep = "\t", quote = F)


## Harmonize add category and combine
impEdge$category = rep("impacted",length(rownames(impEdge)))
nonimpEdge$category = rep("nonimpacted",length(rownames(nonimpEdge)))

impEdge = impEdge[,order(colnames(impEdge))]
nonimpEdge = nonimpEdge[,order(colnames(nonimpEdge))]

df.edge = rbind(impEdge, nonimpEdge)
colnames(df.edge)

df.edge = df.edge %>% left_join(edge.df, by = "Label") ## add info whether taxa were unique or shared nodes
## write.table(df.edge, "COI_edge_tables_combined.txt", sep = "\t", quote = F)



#### COMPARATIVE EDGE BOXPLOTS ####



#### Overall comparison between impacted and non-impacted ####
p = ggboxplot(df.edge, x = "category", y = "qval",
              color = "category", palette = "jco")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Environmental status", y = "Adjusted p-values")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")


### Compare only common edges ### 
df.edge.f = filter(df.edge, df.edge$edge.venn == "commonEdge")
## Comparison between impacted and non-impacted among the common edges ####
p = ggboxplot(df.edge.f, x = "category", y = "method_number",
              color = "category", palette = "jco")
p = p + stat_compare_means(method = "t.test", paired = FALSE, label = "p")
p + 
  labs(x = "Environmental status", y = "Method number of common edges")+
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) + theme(legend.position = "none")



