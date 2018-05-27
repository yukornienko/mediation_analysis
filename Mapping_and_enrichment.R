###Smoke6 and TCDD conc####


enhansers <- read.table("./TCDD_Smoke6_bedopsed_enhansers.bed", header = FALSE, sep = "\t")
colnames(enhansers) <- c("chr", "pos_CpG", "closest_enh_start", "closest_enh_start", "targets", "dist")

promoters_ensembl <- read.table("./TCDD_Smoke6_bedopsed_ensembl.tsv", header = FALSE, sep = "\t")
colnames(promoters_ensembl) <- c("chr", "pos_CpG", "start", "stop",  "gene_name", "dist")

promoters_ensembl_good <- promoters_ensembl[abs(promoters_ensembl$dist) < 500,]

promoters_ensembl_good$gene_name

distances <- read.table("./dist_to_genes.bed", sep = "\t", header = FALSE)
colnames(distances) <- c("chr", "pos_CpG", "pos_CpG1", "chr1", "start", "stop" , "gene_name", "strand", "dist")
hist(distances$dist, breaks = 500, main = "CpG distances to genes", xlab = "distance, bp", xlim = c(-100000,100000))


intra_gene_CpGs <- distances[which(distances$dist == 0),]
intra_gene_CpGs$dist = NULL

write.table(intra_gene_CpGs, "./intra_genetic_CpGs_Sm6TCDD.tsv", sep = "\t", row.names = FALSE)

enhancers <- read.table("./dist_to_enhansers.bed", sep = "\t", header = FALSE)
colnames(enhancers) <- c("chr", "pos_CpG", "pos_CpG1", "chr1", "start", "stop" , "gene_name", "dist")
intra_enh_CpGs <- enhancers[which(enhancers$dist == 0),]

intra_enh_CpGs$pos_CpG1 = NULL

write.table(intra_enh_CpGs, "./intra_enh_CpGs_Sm6TCDD.tsv", sep = "\t", row.names = FALSE)



promoters_ensembl_good$type <- "promoter"
intra_gene_CpGs$type <- "intragenic"
intra_enh_CpGs$type <- "enhancer"
all_genes_TCDD_Smoke6 <- promoters_ensembl_good[,c(1:5,7)]
all_genes_TCDD_Smoke6 <- rbind(all_genes_TCDD_Smoke6, intra_enh_CpGs)
all_genes_TCDD_Smoke6 <- rbind(all_genes_TCDD_Smoke6, intra_gene_CpGs[,c(1:5,7)])

write.table(all_genes_TCDD_Smoke6, "./all_genes_TCDD_Smoke6.tsv", sep = "\t", row.names = F)


#smoke_only#

promoters_smoke <- read.table("./smoke_promoters_ensembl.bed", sep = "\t", header = F)
enhancers_smoke <- read.table("./smoke_enhansers.bed", sep = "\t", header = F)
intragenic_smoke <- read.table("./smoke_intra_genic.bed", sep = "\t", header = F)

colnames(promoters_smoke) <- c("chr", "pos_CpG", "pos_CpG1", "V4", "chr1", "start", "stop", "gene_name", "dist")
colnames(enhancers_smoke) <- c("chr", "pos_CpG", "pos_CpG1", "V4", "chr1", "start", "stop", "gene_name", "dist")
colnames(intragenic_smoke) <- c("chr", "pos_CpG", "pos_CpG1", "V4", "chr1", "start", "stop", "gene_name", "strand", "dist")

promoters_smoke$pos_CpG1 = NULL
promoters_smoke$V4 = NULL
promoters_smoke$chr1 = NULL
enhancers_smoke$pos_CpG1 = NULL
enhancers_smoke$V4 = NULL
enhancers_smoke$chr1 = NULL
intragenic_smoke$pos_CpG1 = NULL
intragenic_smoke$V4 = NULL
intragenic_smoke$chr1 = NULL

promoters_smoke <- promoters_smoke[abs(promoters_smoke$dist) <500,  ]
enhancers_smoke <- enhancers_smoke[enhancers_smoke$dist  == 0,  ]
intragenic_smoke <- intragenic_smoke[intragenic_smoke$dist == 0,]
promoters_smoke$type = "promoter"
enhancers_smoke$type = "enhancer"
intragenic_smoke$type = "intragenic"
all_genes_smoke <- rbind(promoters_smoke[,c(1:5,7)], enhancers_smoke[,c(1:5,7)])
all_genes_smoke <- rbind(all_genes_smoke, intragenic_smoke[,c(1:5,8)])

write.table(all_genes_smoke, "./all_genes_Smoke6.tsv", sep = "\t", row.names = F)

all_genes_smoke <- read.table("./all_genes_Smoke6.tsv", sep = "\t", header = TRUE)

#TCDD_only#

promoters_TCDD <- read.table("./TCDD_promoters_ensembl.bed", sep = "\t", header = F)
enhancers_TCDD <- read.table("./TCDD_enhansers.bed", sep = "\t", header = F)
intragenic_TCDD <- read.table("./TCDD_intra_genic.bed", sep = "\t", header = F)

colnames(promoters_TCDD) <- c("chr", "pos_CpG", "pos_CpG1", "V4", "chr1", "start", "stop", "gene_name", "dist")
colnames(intragenic_TCDD) <- c("chr", "pos_CpG", "pos_CpG1", "V4", "chr1", "start", "stop", "gene_name", "strand", "dist")

promoters_TCDD$pos_CpG1 = NULL
promoters_TCDD$V4 = NULL
promoters_TCDD$chr1 = NULL
intragenic_TCDD$pos_CpG1 = NULL
intragenic_TCDD$V4 = NULL
intragenic_TCDD$chr1 = NULL
promoters_TCDD <- promoters_TCDD[abs(promoters_TCDD$dist) <500,  ]
intragenic_TCDD <- intragenic_TCDD[intragenic_TCDD$dist == 0,]
promoters_TCDD$type = "promoter"
intragenic_TCDD$type = "intragenic"
all_genes_TCDD <- rbind(promoters_TCDD[,c(1:5,7)], intragenic_TCDD[,c(1:5,8)])
write.table(all_genes_TCDD, "./all_genes_TCDD.tsv", sep = "\t", row.names = F)

#######ENRICHMENT!##########
setwd("./project")
library(pathview)
library(reactome.db)
library(ggplot2)
library(gplots)
library(KEGGgraph)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(InterMineR)

#Smoke6 and TCDD conc#

#convert into KEGG IDs:
all_genes_TCDD_Smoke6$gene_name <- as.character(all_genes_TCDD_Smoke6$gene_name)
sym = all_genes_TCDD_Smoke6$gene_name
all_genes_TCDD_Smoke6$KEGG_ID = mget(sym, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
all_genes_TCDD_Smoke6$KEGG_ID <- as.character(all_genes_TCDD_Smoke6$KEGG_ID)
#write.table(all_genes_TCDD_Smoke6, "./all_genes_TCDD_Smoke6.tsv", sep = "\t", row.names = F)
all_genes_TCDD_Smoke6_1 <- all_genes_TCDD_Smoke6[all_genes_TCDD_Smoke6$KEGG_ID !="NA",]

#KEGG
kk <- enrichKEGG(all_genes_TCDD_Smoke6_1$KEGG_ID, organism = 'hsa', pvalueCutoff = 0.05)
head(kk) #hsa04211 10/196 Longevity regulating pathway p.adjust: 0.02762538
cnetplot(kk, categorySize="pvalue",  vertex.label.cex=0.8, vertex.label.color='black')


#REACTOME
react <-enrichPathway(all_genes_TCDD_Smoke6_1$KEGG_ID, pvalueCutoff=0.05, readable=T)
head(as.data.frame(react)) #1 - R-HSA-8948216 COL7A1/COL11A1/COL13A1/COL4A1/COL4A2/COL23A1/COL26A1/COL5A1 - Collagen chain trimerization
#8/295, p.adjust 0.02242243
cnetplot(react, categorySize="pvalue",  vertex.label.cex=0.7, vertex.label.color='black')
barplot(react,drop = TRUE, showCategory = 9)

#GO
ego <- enrichGO(all_genes_TCDD_Smoke6_1$KEGG_ID, OrgDb = org.Hs.eg.db, ont  = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.05, readable = TRUE)
head(ego)
str(ego)
cnetplot(ego, categorySize=15, category.color = "blue", fixed = T, category.label = 1, vertex.label.cex=0.5, vertex.label.color='black', colorEdge = FALSE, circular = FALSE, node_label = TRUE)
barplot(ego,drop = TRUE, showCategory = 10)


#Smoke6 only#

#convert into KEGG IDs:
all_genes_smoke$gene_name <- as.character(all_genes_smoke$gene_name)
sym = all_genes_smoke$gene_name
all_genes_smoke$KEGG_ID = mget(sym, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
all_genes_smoke$KEGG_ID <- as.character(all_genes_smoke$KEGG_ID)
write.table(all_genes_smoke, "./all_genes_Smoke6.tsv", sep = "\t", row.names = F)
all_genes_smoke_1 <- all_genes_smoke[all_genes_smoke$KEGG_ID !="NA",]

#KEGG
kk <- enrichKEGG(all_genes_smoke_1$KEGG_ID, organism = 'hsa', pvalueCutoff = 0.05)
head(kk) #hsa04211 9/187 Longevity regulating pathway p.adjust: 0.04750666
#hsa04072 12/187 Phospholipase D signaling pathway 0.04750666
cnetplot(kk, categorySize="pvalue",  vertex.label.cex=0.8, vertex.label.color='black')


#REACTOME
react <-enrichPathway(all_genes_smoke_1$KEGG_ID, pvalueCutoff=0.05, readable=T)
head(as.data.frame(react)) #0! 
cnetplot(react, categorySize="pvalue",  vertex.label.cex=0.7, vertex.label.color='black')
barplot(react,drop = TRUE, showCategory = 9)

#GO
ego <- enrichGO(all_genes_smoke_1$KEGG_ID, OrgDb = org.Hs.eg.db, ont  = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.05, readable = TRUE)
head(ego)
barplot(ego,drop = TRUE, showCategory = 20)
cnetplot(ego, categorySize=15, category.color = "blue", fixed = T, category.label = 1, vertex.label.cex=0.5, vertex.label.color='black', colorEdge = FALSE, circular = FALSE, node_label = TRUE)

#TCDD only#

#convert into KEGG IDs:
all_genes_TCDD$gene_name <- as.character(all_genes_TCDD$gene_name)
sym = all_genes_TCDD$gene_name
all_genes_TCDD$KEGG_ID = mget(sym, revmap(org.Hs.egSYMBOL),ifnotfound=NA)
all_genes_TCDD$KEGG_ID <- as.character(all_genes_TCDD$KEGG_ID)
write.table(all_genes_TCDD, "./all_genes_TCDD.tsv", sep = "\t", row.names = F)
all_genes_TCDD_1 <- all_genes_TCDD[all_genes_TCDD$KEGG_ID !="NA",]

#KEGG
kk <- enrichKEGG(all_genes_TCDD_1$KEGG_ID, organism = 'hsa', pvalueCutoff = 0.05)
head(kk) #0
cnetplot(kk, categorySize="pvalue",  vertex.label.cex=0.8, vertex.label.color='black')


#REACTOME
react <-enrichPathway(all_genes_TCDD_1$KEGG_ID, pvalueCutoff=0.05, readable=T)
head(as.data.frame(react)) #2: R-HSA-1980143 R-HSA-1980143 Signaling by NOTCH1       2/5 0.02158359 NUMB/NCOR1     2
#R-HSA-157118   R-HSA-157118  Signaling by NOTCH       2/5 ...

cnetplot(react, categorySize="pvalue",  vertex.label.cex=0.7, vertex.label.color='black')
barplot(react,drop = TRUE, showCategory = 9)

#GO
ego <- enrichGO(all_genes_TCDD_1$KEGG_ID, OrgDb = org.Hs.eg.db, ont  = "CC", pAdjustMethod = "BH", pvalueCutoff  = 0.05, readable = TRUE)
head(ego)
barplot(ego,drop = TRUE, showCategory = 10)
cnetplot(ego, categorySize=15, category.color = "blue", fixed = T, category.label = 1, vertex.label.cex=1, vertex.label.color='black', colorEdge = FALSE, circular = FALSE, node_label = TRUE)

#
