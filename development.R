library(Seurat)

out.files = dir("~/JesusProject/CellRangerOUTS/")
out.files = out.files[!out.files== "Adult"] # Remove adult data

# Read raw expression matrices into a list
raw.data = list()
for(i in out.files){
  raw.data[[i]] = Read10X(data.dir = paste0("~/JesusProject/CellRangerOUTS/",i,"/outs/filtered_feature_bc_matrix/"))
}

# Create list of seurat objects
seurat.objects = list()
for(i in names(raw.data)){
  seurat.objects[[i]] = CreateSeuratObject(counts = raw.data[[i]], project = paste0(i), min.cells = 3, min.features = 200)
}

# Combine all objects into one seurat.
combined = merge(seurat.objects[["E12"]], y = c(seurat.objects[["E14"]], seurat.objects[["E16"]],
                                                seurat.objects[["E18"]], seurat.objects[["P2"]],
                                                seurat.objects[["P8"]], seurat.objects[["P15"]],
                                                seurat.objects[["P30"]]),
                 add.cell.ids = c("E12","E14","E16","E18","P2","P8","P15","P30"), project = "Development")
rm(raw.data); rm(seurat.objects) # Clean up environment

# Add experimental and sequencing batch to the metadata
ExperimentalBatch = c("ExpBatchA","ExpBatchA","ExpBatchB","ExpBatchB","ExpBatchC","ExpBatchD","ExpBatchD","ExpBatchE")
names(ExperimentalBatch) = c("E12","E18","E14","E16","P2","P8","P15","P30")
SequencingBatch = c("SeqBatch1","SeqBatch1",rep("SeqBatch2",5),"SeqBatch3")
names(ExperimentalBatch) = c("E12","E18","E14","E16","P2","P8","P15","P30")
combined@meta.data$exp.batch = ExperimentalBatch[combined$orig.ident]
combined@meta.data$seq.batch = SequencingBatch[combined$orig.ident]

# Add mitochondrial pecentage to the metadata
combined@meta.data$percent.mito = (Matrix::colSums(combined@assays$RNA@counts[grep(rownames(combined@assays$RNA@counts), pattern = "^mt-"),])/Matrix::colSums(combined@assays$RNA@counts))*100

# Calculate cell cycle signature to regress out 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Convert the cycling gene list to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

s.genes <- convertHumanGeneList(s.genes)
g2m.genes <- convertHumanGeneList(g2m.genes)

combined <- CellCycleScoring(combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# Order developmental stages and plot features
combined@meta.data$orig.ident = factor(x = combined@meta.data$orig.ident, levels = c("E12","E14","E16","E18","P2","P8","P15","P30"))
VlnPlot(object = combined, features = "percent.mito", split.by = "orig.ident", pt.size = 0)
VlnPlot(object = combined, features = "nFeature_RNA", split.by = "orig.ident", pt.size = 0)

# Filter cells by high/low counts and mitochondrial gene percent
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 20) # Filter out by number of genes and mitochondrial content

# Apply SCTransform to normalize the data and regress out mito and batch
combined <- SCTransform(combined, vars.to.regress = c("nCount_RNA","percent.mito","S.Score","G2M.Score"))

# Standard Seurat workflow
combined <- RunPCA(combined, npcs = 50)
DimHeatmap(combined, dims = 1:6, cells = 500, balanced = TRUE)
DimHeatmap(combined, dims = 7:12, cells = 500, balanced = TRUE)
DimHeatmap(combined, dims = 13:18, cells = 500, balanced = TRUE)
DimHeatmap(combined, dims = 19:24, cells = 500, balanced = TRUE)
DimHeatmap(combined, dims = 25:30, cells = 500, balanced = TRUE)

ElbowPlot(object = combined)

combined <- RunUMAP(combined, dims = 1:24, verbose = FALSE)

combined <- FindNeighbors(combined, dims = 1:24, verbose = FALSE)
combined <- FindClusters(combined, verbose = FALSE)

DimPlot(combined, label = TRUE) + NoLegend()
DimPlot(combined,group.by = "orig.ident", label = T) + NoLegend()
DimPlot(combined,group.by = "orig.ident", label = F) + NoLegend()
DimPlot(combined,group.by = "orig.ident", label = F, cols='Dark2' ) + NoLegend()

combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)

library(dplyr)
top_20 = combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
#FeaturePlot(object = combined, features = "Sdc1")
write.csv(x = top_20, file = "/Volumes/Samsung_T5/Single Cell Analysis/Franco new analysis/24dimensions_Top_20.csv")
dev.off
save(combined, file = "~/CH/ProcessedData/development.combined.Rdata")

genes = c("Marco","Cd74","Cd24a","Car2","Ptprc","Cdh5","Pecam1")
genes = c("Hbb-bt","Hbb-bs","Gata1","Tal1","Csf1r","Pf4","Gypa","Cd9", "Pbx1")
genes = c("Kit","Hlf","Myc")
genes = c('Marco')

library(ggplot2)
library(RColorBrewer)
library(viridis)


FeaturePlot(combined, features = genes) + 
  scale_color_gradientn(colours = (brewer.pal(n = 11, name = "Reds")))
FeaturePlot(combined, 'Fbln2') + 
  scale_color_gradientn(colours = (brewer.pal(n = 11, name = "Reds")))

FeaturePlot(combined, features = genes)
FeaturePlot(combined, features = genes) + 
  scale_color_viridis_c(option = 'B')


FeaturePlot(combined, features = genes , scale_color_manual(colours = rev(brewer.pal(n = 11, name = "RdBu"))))

p1 <- FeaturePlot(combined, features = genes, cols='red')
fix.sc <- scale_color_gradientn( colours = c('lightgrey', 'blue'),  limits = c(1, 8))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

VlnPlot(combined, features = "Hbb-bs")
VlnPlot(combined, features = "Fcgr2b")
VlnPlot(combined,group.by = "orig.ident", features = "Maf")

load(file = "~/CH/ProcessedData/development.combined.Rdata")


#basic plot of clusters by replicate
library(ggplot2)
library(RColorBrewer)


pt <- table(Idents(combined), combined$orig.ident)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank())

pt$Var1 <- factor(pt$Var1,levels = c("0","1", "2", "3", "4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"))

ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 1) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12,"Paired")) +
  theme_classic()


#Violing plots by timepoint using color Paired
VlnPlot(combined,group.by = "orig.ident", 
        features = "Maf", 
        cols = brewer.pal(12,"Paired"),
        pt.size = 0)
