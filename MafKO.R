# Process data from scRNA - Maf KO vs WT endothelial cells
library(Seurat)
library(dplyr)
library(reshape2)
library(viridis)
library(ggplot2)
library(future)
options(future.globals.maxSize = 4000 * 1024^2)

# Load data and merge to create a unique object
wt = Read10X(data.dir = "~/JesusProject/scRNA_MafKO/CellRanger_WT/outs/filtered_feature_bc_matrix/")
ko = Read10X(data.dir = "~/JesusProject/scRNA_MafKO/CellRanger_MafKO/outs/filtered_feature_bc_matrix/")
wt.obj = CreateSeuratObject(counts = wt, min.cells = 10, min.features = 100)
ko.obj = CreateSeuratObject(counts = ko, min.cells = 10, min.features = 100)
rm(wt); rm(ko) # Clean up environment

data <- merge(wt.obj, y = ko.obj, add.cell.ids = c("WT", "KO"), project = "MafKO")
rm(wt.obj); rm(ko.obj) # Clean up environment

# Add genotype information
data@meta.data$Genotype = unlist(lapply(strsplit(rownames(data@meta.data),"_"), function(x) x[1]))
Idents(data) = data@meta.data$Genotype

# Check mitochondrial reads
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Subset by features and mito percent
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)

# Sex features
bad.features <- c("Ptprb", "Xist", "Ddx3y", "Eif2s3y", "Erdr1", "Gm29650", "Kdm5d", "Uty", "Kap")
hist(data@assays$RNA@counts["Xist",], breaks = 100)

data@meta.data$Sex = "Male"
data@meta.data$Sex[data@assays$RNA@counts["Xist",] > 0] = "Female"

# Cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

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

data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

data.list <- SplitObject(data, split.by = "Sex")
for (i in 1:length(data.list)) {
  data.list[[i]] <- SCTransform(data.list[[i]], verbose = TRUE)
}

#Shared features
shared.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)

#Remove the list of excluded features
integration.features <- setdiff(shared.features, bad.features)

#Verify Pearson residuals
pancreas.list <- PrepSCTIntegration(object.list = data.list,
                                    anchor.features = integration.features,
                                    assay = "SCT", verbose = TRUE)

#Find anchors for integration
anchors <- FindIntegrationAnchors(object.list = data.list,
                                  dims = 1:40,
                                  anchor.features = integration.features)

#Intergrate datasets
integrated <- IntegrateData(anchorset = anchors,
                            normalization.method = "SCT", 
                                     verbose = FALSE)

integrated <- RunPCA(integrated,npcs = 50, verbose = TRUE)
ElbowPlot(integrated, ndims = 50)
integrated <- FindNeighbors(integrated, dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.5)
integrated <- RunUMAP(integrated, dims = 1:20)

DimPlot(integrated, reduction = "umap")
DimPlot(integrated, reduction = "umap", group.by = "Sex")

m = melt(apply(table(integrated@meta.data$Sex, integrated@meta.data$seurat_clusters),2,function(x) x/sum(x)))
ggplot(m, aes(x = as.factor(Var2), y = value, fill = Var1)) + 
  geom_bar(stat = "identity")
integrated@meta.data$Sample = paste0(integrated@meta.data$Genotype,"_", integrated@meta.data$Sex)
table(integrated@meta.data$Sample)

markers = FindAllMarkers(integrated, logfc.threshold = 0.2,min.pct = 0.25, test.use = "LR", only.pos = T)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(integrated, features = top10$gene) + NoLegend()
top1 <- markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)

VlnPlot(integrated, features = top1$gene[1:5], ncol = 1, pt.size = 0)

modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin )
  return(p)
}

extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

StackedVlnPlot(obj = integrated, features = top1$gene[1:5])

DimPlot(integrated, label = T)
DimPlot(integrated, reduction = "umap", group.by = "Genotype", split.by = "Sample")

save(integrated, file = "~/JesusProject/ProcessedData/integrated.MafKO.Rdata")
load(file = "~/JesusProject/ProcessedData/integrated.MafKO.Rdata")

FeaturePlot(integrated, features = "Cd24a")
FeaturePlot(integrated, features = "Maf", split.by = "Genotype")
FeaturePlot(integrated, features = "Rspo3")
FeaturePlot(integrated, features = "Mrc1")
FeaturePlot(integrated, features = "Ly6a")
FeaturePlot(integrated, features = "Ptprc")
FeaturePlot(integrated, features = "Cdh5")

markers.int = FindAllMarkers(integrated, logfc.threshold = 0.25, only.pos = T)
top_10 = markers.int %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

write.csv(x = top_10, file = "~/JesusProject/ProcessedData/ClusterMarkers.MafKO.csv", quote = F, row.names = T, col.names = T)

View(top_10)
View(markers.int)
cluster.name = c("S_CV","Sinusoids","Cd24_2","Cd24_1","KOS_1","Cd24_3","K","C1","Tcell_Gzm","Tcell_Cd6","Tcells_Foxp3","KOS_2","Neutrophils","KOS_3","Leukocytes_3","Leukocytes_1","CLP","KOS_4","PV","RBC","Leukocytes_2","Hepatocyte_doublets","Mast_cells")
names(cluster.name) = c(0:22)

integrated@meta.data$cluster_name = cluster.name[integrated@meta.data$seurat_clusters]

Idents(integrated) = integrated@meta.data$cluster_name

DimPlot(integrated, label = T)
FeaturePlot(integrated,"Ptprc")

DoHeatmap(integrated, features = top_10$gene) + NoLegend()
top_1 = markers.int %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)
StackedVlnPlot(obj = integrated, features = top_1$gene[20:22])

integrated@meta.data$UMAP1 = integrated@reductions$umap@cell.embeddings[,1]
integrated@meta.data$UMAP2 = integrated@reductions$umap@cell.embeddings[,2]

integrated@meta.data$Genotype = factor(integrated@meta.data$Genotype, levels = c("WT","KO"))

ggplot(integrated@meta.data, aes(x = UMAP1, y = UMAP2, color = Genotype)) + 
  geom_point(size = 0.5) +
  facet_wrap(~Genotype) +
  scale_color_manual(values = c("#255ea8","#ec7014"))+
  theme_classic()

freq = melt(apply(table(integrated@meta.data$cluster_name, integrated@meta.data$Sample),1,function(x) x/sum(x)))
ggplot(freq, aes(x = as.factor(Var2), y = value, fill = Var1, group = Var1)) +
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_manual(values = c(rep("#ec7014",2),rep("#255ea8",2)))+
  labs(x = "Cluster", y = "Fraction of cells in cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

save(integrated, file = "~/JesusProject/ProcessedData/integrated.MafKO_v1.Rdata")
load("~/JesusProject/ProcessedData/integrated.MafKO_v1.Rdata")

integrated@meta.data
cluster.name = c("S2",
                 "S",
                 "B2",
                 "B1",
                 "KOS_1",
                 "B3",
                 "K",
                 "C1",
                 "T-Gzm",
                 "T-Cd6",
                 "T-Foxp3",
                 "KOS_2",
                 "Neutrophils",
                 "C-EC",
                 "Leukocytes_1",
                 "CMP",
                 "CLP",
                 "KOS_3",
                 "PV",
                 "D",
                 "Leukocytes_2",
                 "H",
                 "Mast_cells")

names(cluster.name) = c(0:22)
integrated@meta.data$cluster_name = cluster.name[integrated@meta.data$seurat_clusters]

DimPlot(integrated, group.by = "cluster_name", label = T)
write.csv(x = integrated@meta.data, file = "~/JesusProject/GEO_MafKO/Metadata.csv", row.names = T, quote = F)

tab = table(integrated@meta.data$cluster_name, integrated@meta.data$Genotype)
df.fisher = data.frame(CLP = tab["CLP",], Out = colSums(tab[!(rownames(tab) == "CLP"),]))
fisher.test(df.fisher)
f = apply(tab,2,function(x) x/sum(x)*100)
m = melt(f)

samp = table(integrated@meta.data$cluster_name, paste0(integrated@meta.data$Genotype,"_",integrated@meta.data$Sex))
fsamp = apply(samp,2,function(x) x/sum(x)*100)
msamp = melt(fsamp["CLP",])
msamp$Genotype = unlist(lapply(strsplit(rownames(msamp),"_"), function(x) x[1]))
ggplot(m[m$Var1 == "CLP",], aes(x = Var2, y = value)) +
  geom_bar(stat = "identity", width = 0.5, color = "black", fill = "#39B9C6") +
  geom_point(data = msamp,aes(x = Genotype, y = value), size = 2) +
  labs(x = "", y = "Cells in CLP cluster (%)") +
  theme_classic()

# Percentage of Maf + cells in KO
load(file = "~/JesusProject/ProcessedData/integrated.MafKO.Rdata")

dat = integrated@meta.data
dat$MAFcounts = integrated@assays$RNA@counts["Maf",]
dat = dat[dat$cluster_name %in% c("KOS","KOS_2","KOS_3","PV","C-EC","S","S2"),]
p = apply(table(dat$MAFcounts>1, dat$Sample),2,function(x) x/sum(x)*100)
median(p[2,1:2]); sd(p[2,1:2])
median(p[2,3:4])
