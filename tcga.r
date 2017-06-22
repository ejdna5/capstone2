#install.packages(tidyverse)
library(tidyverse)

# load data from mutations_extended.txt files for each dataset into R as dataframes
# make tables (dataframes) for each dataset containing only relevant fields
cervical_scc_data <- read.delim2("R/tcga/cervical_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 1)
select_(cervical_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> cervical_scc_table

lung_scc_data <- read.delim2("R/tcga/lung_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 1)
select_(lung_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> lung_scc_table

lung2_scc_data <- read.delim2("R/tcga/lung2_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 1)
select_(lung2_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> lung2_scc_table

esoph_scc_data <- read.delim2("R/tcga/esoph_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 2)
select_(esoph_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> esoph_scc_table

esoph2_scc_data <- read.delim2("R/tcga/esoph2_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 1)
select_(esoph2_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> esoph2_scc_table

skin_scc_data <- read.delim2("R/tcga/skin_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 1)
select_(skin_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> skin_scc_table

hn_scc_data <- read.delim2("R/tcga/hn_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 1)
select_(hn_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> hn_scc_table

hn2_scc_data <- read.delim2("R/tcga/hn2_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 1)
select_(hn2_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> hn2_scc_table

hn3_scc_data <- read.delim2("R/tcga/hn3_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 1)
select_(hn3_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> hn3_scc_table

hn4_scc_data <- read.delim2("R/tcga/hn4_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 1)
select_(hn4_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> hn4_scc_table

hn5_scc_data <- read.delim2("R/tcga/hn5_scc_data_mutations_extended.txt",header = TRUE,sep = "\t", skip = 1)
select_(hn5_scc_data, Gene="Hugo_Symbol",VariantClass="Variant_Classification",VariantType="Variant_Type",SampleBC="Tumor_Sample_Barcode",MutStat="Mutation_Status",SeqSource="Sequence_Source",SampleID="Tumor_Sample_UUID",DNAchange="HGVSc",ProteinChange="HGVSp_Short",SIFT="SIFT",PolyPhen="PolyPhen",Impact="IMPACT") -> hn5_scc_table

# filter tables for mutations predicted to have at least "MODERATE" impact
filter(cervical_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> cervical_scc_table_PotentialImpact

filter(lung_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> lung_scc_table_PotentialImpact

filter(lung2_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> lung2_scc_table_PotentialImpact

filter(esoph_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> esoph_scc_table_PotentialImpact

filter(esoph2_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> esoph2_scc_table_PotentialImpact

filter(skin_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> skin_scc_table_PotentialImpact

filter(hn_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> hn_scc_table_PotentialImpact

filter(hn2_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> hn2_scc_table_PotentialImpact

filter(hn3_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> hn3_scc_table_PotentialImpact

filter(hn4_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> hn4_scc_table_PotentialImpact

filter(hn5_scc_table, Impact=="MODERATE"|Impact=="HIGH") -> hn5_scc_table_PotentialImpact

# Build lists of unique genes for each type of scc
unique(as.character(cervical_scc_table_PotentialImpact$Gene)) -> cervical_scc_Genes

unique(as.character(skin_scc_table_PotentialImpact$Gene)) -> skin_scc_Genes

unique(c(as.character(lung_scc_table_PotentialImpact$Gene),as.character(lung2_scc_table_PotentialImpact$Gene))) -> lung_scc_Genes

unique(c(as.character(esoph_scc_table_PotentialImpact$Gene),as.character(esoph2_scc_table_PotentialImpact$Gene))) -> esoph_scc_Genes

unique(c(as.character(hn_scc_table_PotentialImpact$Gene),as.character(hn2_scc_table_PotentialImpact$Gene),as.character(hn3_scc_table_PotentialImpact$Gene),as.character(hn4_scc_table_PotentialImpact$Gene),as.character(hn5_scc_table_PotentialImpact$Gene))) -> hn_scc_Genes

# Plot Venn diagram showing overlap of lists and return list of genes in common
#install.packages("gplots")
library(gplots)
VennData <- venn(list(Cervical=cervical_scc_Genes,Esophageal=esoph_scc_Genes,HeadAndNeck=hn_scc_Genes,Lung=lung_scc_Genes,Skin=skin_scc_Genes),show.plot = TRUE)
venn(list(Cervical=cervical_scc_Genes,Esophageal=esoph_scc_Genes,Head_And_Neck=hn_scc_Genes,Lung=lung_scc_Genes,Skin=skin_scc_Genes),show.plot = TRUE)
OverlappingGenes <- attr(VennData, "intersections")$`Cervical:Esophageal:HeadAndNeck:Lung:Skin`

# Run KEGG pathway enrichment analysis on common list of squamous cell carcinoma genes
# Convert Gene Names to Entrez IDs
source("https://bioconductor.org/biocLite.R")
biocLite("mygene")
biocLite("limma")
library(mygene)
queryMany(OverlappingGenes, scopes="symbol", fields="entrezgene", species="human") -> OverlappingGenesEntrezTable
# Use limma to look for enriched KEGG pathways
library(limma)
KEGGpaths <- kegga(OverlappingGenesEntrezTable$entrezgene[!is.na(OverlappingGenesEntrezTable$entrezgene)],species = "Hs") %>% 
  arrange(P.DE)

# Map list to KEGG pathway for chosen pathway
biocLite("pathview")
library(pathview)
pathview(gene.data = OverlappingGenesEntrezTable$`_id`, cpd.data = NULL, "hsa01521")

# Generate matrix of mutated genes for each study; 0 not mutated in study; 1 mutated in study
AllGenes <- unique(c(lung_scc_Genes,skin_scc_Genes,cervical_scc_Genes,hn_scc_Genes,esoph_scc_Genes))
GeneMatrix <- matrix(0, ncol = 11, nrow = 17686, dimnames = list(AllGenes,c("cervical","lung1","lung2","esoph1","esoph2","skin","hn1","hn2","hn3","hn4","hn5")))
for (i in 1:length(cervical_scc_Genes)){
  GeneMatrix[cervical_scc_Genes[i],"cervical"]=1
}
for (i in 1:length(lung_scc_table_PotentialImpact$Gene)){
  GeneMatrix[as.character(lung_scc_table_PotentialImpact$Gene[i]),"lung1"]=1
}  
for (i in 1:length(lung2_scc_table_PotentialImpact$Gene)){
  GeneMatrix[as.character(lung2_scc_table_PotentialImpact$Gene[i]),"lung2"]=1
}
for (i in 1:length(esoph_scc_table_PotentialImpact$Gene)){
  GeneMatrix[as.character(esoph_scc_table_PotentialImpact$Gene[i]),"esoph1"]=1
}
for (i in 1:length(esoph2_scc_table_PotentialImpact$Gene)){
  GeneMatrix[as.character(esoph2_scc_table_PotentialImpact$Gene[i]),"esoph2"]=1
}
for (i in 1:length(skin_scc_Genes)){
  GeneMatrix[skin_scc_Genes[i],"skin"]=1
}
for (i in 1:length(hn_scc_table_PotentialImpact$Gene)){
  GeneMatrix[as.character(hn_scc_table_PotentialImpact$Gene[i]),"hn1"]=1
}
for (i in 1:length(hn2_scc_table_PotentialImpact$Gene)){
  GeneMatrix[as.character(hn2_scc_table_PotentialImpact$Gene[i]),"hn2"]=1
}
for (i in 1:length(hn3_scc_table_PotentialImpact$Gene)){
  GeneMatrix[as.character(hn3_scc_table_PotentialImpact$Gene[i]),"hn3"]=1
}
for (i in 1:length(hn4_scc_table_PotentialImpact$Gene)){
  GeneMatrix[as.character(hn4_scc_table_PotentialImpact$Gene[i]),"hn4"]=1
}
for (i in 1:length(hn5_scc_table_PotentialImpact$Gene)){
  GeneMatrix[as.character(hn5_scc_table_PotentialImpact$Gene[i]),"hn5"]=1
}

#generate distance matrix and cluster data (hierarchical)
TgeneMatrix <- t(GeneMatrix)
d <- dist(TgeneMatrix)
hc <- hclust(d)

#K-means clustering: k = 2,3,4, and 5 
kmeans(TgeneMatrix,2) -> K2clust
kmeans(TgeneMatrix,3) -> K3clust
kmeans(TgeneMatrix,4) -> K4clust
kmeans(TgeneMatrix,5) -> K5clust

ggplot(as.data.frame(K2clust$cluster), aes(x=factor(1), fill=factor(K2clust$cluster))) + 
  geom_bar(width=1) + 
  coord_polar(theta="y") + 
  scale_fill_discrete(name="Cluster",
                      breaks=c("1", "2"),
                      labels=c("esoph1, esoph2, skin, hn1, hn2, hn3, hn5","cervical, hn4, lung1, lung2")) +
  theme(panel.background=element_blank(),
        axis.title=element_blank(), 
        axis.ticks=element_blank(), 
        panel.grid=element_blank(), 
        axis.text.y=element_blank(),
        axis.text.x=element_blank()) 

ggplot(as.data.frame(K3clust$cluster), aes(x=factor(1), fill=factor(K3clust$cluster))) + 
  geom_bar(width=1) + 
  coord_polar(theta="y") +
  scale_fill_discrete(name="Cluster",
                      breaks=c("1", "2", "3"),
                      labels=c("cervical, hn4", "esoph1, esoph2, skin, hn1, hn2, hn3, hn5", "lung1, lung2")) +
  theme(panel.background=element_blank(),
        axis.title=element_blank(), 
        axis.ticks=element_blank(), 
        panel.grid=element_blank(), 
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggplot(as.data.frame(K4clust$cluster), aes(x=factor(1), fill=factor(K4clust$cluster))) + 
  geom_bar(width=1) + 
  coord_polar(theta="y") +
  scale_fill_discrete(name="Cluster",
                      breaks=c("1", "2", "3", "4"),
                      labels=c( "hn1","esoph2","cervical, lung1, lung2, hn4", "esoph1, skin, hn2, hn3, hn5")) +
  theme(panel.background=element_blank(),
        axis.title=element_blank(), 
        axis.ticks=element_blank(), 
        panel.grid=element_blank(), 
        axis.text.y=element_blank(),
        axis.text.x=element_blank())

ggplot(as.data.frame(K5clust$cluster), aes(x=factor(1), fill=factor(K5clust$cluster))) + 
  geom_bar(width=1) + 
  coord_polar(theta="y") +
  scale_fill_discrete(name="Cluster",
                      breaks=c("1", "2", "3", "4","5"),
                      labels=c( "hn4","esoph1, skin, hn2, hn3,n5","lung1, lung2","cervical","esoph2, hn1")) +
  theme(panel.background=element_blank(),
        axis.title=element_blank(), 
        axis.ticks=element_blank(), 
        panel.grid=element_blank(), 
        axis.text.y=element_blank(),
        axis.text.x=element_blank())