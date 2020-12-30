# Trajectory Analysis of CD8 T-cells

library(dplyr)
library(slingshot)
library(SeuratDisk)
library(Seurat)
library(knitr)
library(kableExtra)
library(ggplot2)
library(RColorBrewer)

# Filtering/Normalization -------------------------------------------------

load(file="GSE141784_sub_prenorm_seurat_object.RData")

GSE141784.CD8 <- subset(experiment.aggregate,
                            subset = Cells == c("CD8"))
rm(experiment.aggregate)

GSE141784.CD8 <- NormalizeData(
  object = GSE141784.CD8,
  normalization.method = "LogNormalize",
  scale.factor = 10000)


# Variable Gene Identification --------------------------------------------

GSE141784.CD8 <- FindVariableFeatures(
  object = GSE141784.CD8,
  selection.method = "vst",
  nfeatures = 3000)

length(VariableFeatures(GSE141784.CD8))

top20_vargenes <- head(VariableFeatures(GSE141784.CD8), 20)

top20_vargenes

vargenep1 <- VariableFeaturePlot(GSE141784.CD8)
vargenep1 <- LabelPoints(plot = vargenep1, points = top20_vargenes, repel = TRUE)

png("plots/top20_vargenes_cd8subset.png", units="in", width=8, height=5, res=300)
vargenep1
dev.off()

save(GSE141784.CD8, file="data/cd8_subset/GSE141784_normalized_cd8_subset.RData")
load(file="data/cd8_subset/GSE141784_normalized_cd8_subset.RData")

# Splitting the data ------------------------------------------------------

GSE141784.list <- SplitObject(GSE141784.CD8, split.by = "orig.ident")

# for now I will leave this the way it is but I will optimize with purrr
for (i in 1:length(GSE141784.list)) {
  GSE141784.list[[i]] <- NormalizeData(GSE141784.list[[i]], 
                                       normalization.method = "LogNormalize",
                                       scale.factor = 10000)
  GSE141784.list[[i]] <- FindVariableFeatures(GSE141784.list[[i]], 
                                              selection.method = "vst", 
                                              nfeatures = 3000, 
                                              verbose = FALSE)
}


reference.list <- GSE141784.list[c("NOD_4w", "NOD_8w", "NOD_15w")]
GSE141784.anchors <- FindIntegrationAnchors(object.list = reference.list, 
                                            dims = 1:50, anchor.features = 3000)
GSE141784.int_CD8 <- IntegrateData(anchorset = GSE141784.anchors, 
                                      dims = 1:50)

DefaultAssay(GSE141784.int_CD8) <- "integrated"


save(GSE141784.int_CD8, file="data/cd8_subset/GSE141784_integrated_cd8_subset.RData")
load(file="data/cd8_subset/GSE141784_integrated_cd8_subset.RData")


# Scaling Data ------------------------------------------------------------

GSE141784.int_CD8 <- ScaleData(GSE141784.int_CD8,
                                  vars.to.regress = c(#"cell.cycle", "percent.mito",
                                                      "percent.ribo", "nCount_RNA"))

GSE141784.int_CD8 <- RunPCA(GSE141784.int_CD8)
DimHeatmap(GSE141784.int_CD8, 
           dims = 1:6, 
           cells = 500, 
           balanced = TRUE)
ElbowPlot(GSE141784.int_CD8, ndims = 40)

GSE141784.int_CD8 <- JackStraw(GSE141784.int_CD8, 
                                  dims = 50)
GSE141784.int_CD8 <- ScoreJackStraw(GSE141784.int_CD8, dims = 1:50)

png("plots/cd8subset_integrated_jackstraw.png", units="in", width=10, height=5, res=300)
JackStrawPlot(GSE141784.int_CD8, dims = 1:50)
dev.off()

use.pcs <- 1:15

GSE141784.int_CD8 <- FindNeighbors(GSE141784.int_CD8, 
                                      reduction="pca",
                                   nn.method = "annoy",
                                   k.param = 20,
                                   annoy.metric = "cosine",
                                      dims = use.pcs)

GSE141784.int_CD8 <- FindClusters(
  object = GSE141784.int_CD8, 
  resolution = seq(0.25,4,0.25),
  algorithm = 2
)

sapply(grep("res",colnames(GSE141784.int_CD8@meta.data),value = TRUE),
       function(x) length(unique(GSE141784.int_CD8@meta.data[,x])))

Idents(GSE141784.int_CD8) <- "integrated_snn_res.0.75"

GSE141784.int_CD8 <- RunTSNE(object = GSE141784.int_CD8, 
                                reduction.use = "pca", 
                                dims = use.pcs,
                             perplexity = 30,
                                do.fast = TRUE)

GSE141784.int_CD8 <- RunUMAP(GSE141784.int_CD8, 
                                reduction = "pca", 
                                dims = use.pcs)
p1 <- DimPlot(GSE141784.int_CD8, 
  reduction = "tsne") +
  scale_color_brewer(type = "qual", 
    palette = "Set1")
p1_5 <- DimPlot(GSE141784.int_CD8, reduction = "tsne", 
                label = TRUE, 
                label.size = 3, 
                #label.color = "white",
                label.box = TRUE, 
                repel = TRUE) +
  scale_color_brewer(type = "qual", 
    palette = "Set1") + 
  NoLegend()
p2 <- DimPlot(GSE141784.int_CD8, 
  reduction = "umap", 
  label = TRUE, 
  repel = TRUE) + 
  scale_color_brewer(type = "qual", 
    palette = "Set1") +
  NoLegend()

png("plots/cd8subset_integrated_foundclusters_tsne_umap.png", units="in", width=10, height=5, res=300)
cowplot::plot_grid(p1, p2)
dev.off()

save(GSE141784.int_CD8, file = "data/cd8_subset/GSE141784_integrated_clustered_cd8_subset.RData")
load(file="data/cd8_subset/GSE141784_integrated_clustered_cd8_subset.RData")

DimPlot(GSE141784.int_CD8, reduction = "tsne", split.by = "Sample") + 
  scale_color_brewer(type = "qual", palette = "Set1")

FeaturePlot(GSE141784.int_CD8, features = c('nFeature_RNA', 'nCount_RNA'), pt.size=0.5, reduction = "tsne") & scale_color_viridis_c()

DimPlot(object = GSE141784.int_CD8, pt.size=0.5, group.by = "cell.cycle", reduction = "tsne" )


GSE141784.int_CD8 <- BuildClusterTree(
  GSE141784.int_CD8, dims = use.pcs)

PlotClusterTree(GSE141784.int_CD8)


# Finding Cluster Biomarkers ----------------------------------------------

DefaultAssay(GSE141784.int_CD8) <- "RNA"

markers_all <- FindAllMarkers(
  object = GSE141784.int_CD8, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  thresh.use = 0.25,
  test.use = "bimod"
)

top5 <- markers_all %>% group_by(cluster) %>% top_n(5, avg_logFC)

top10 <- markers_all %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(
  object = GSE141784.int_CD8, 
  features = top10$gene
) 

markers <- c("Entpd1","Havcr2", "Tigit", "Klrg1","Lag3", "Pdcd1", "Batf", "Tox","Tox2", "Id2",
                 "Ccl3", "Prf1", "Fasl", "Gzmb", "Gzma", "Gzmk", "Ccl5","Hmgb2","Junb",
                 "Ly6c2", "Cxcr3", "Id3", "Slamf6", "Il2ra", "Xcl1", "Tcf7",
                 "Ccr7", "S1pr1", "Lef1","Sell", "Il7r","Cd28", "Ifng", "Cd40lg")

cytotoxic_markers <- c("Gzma", "Gzmb", "Gzmk", "Klrc1", "Klrc2", "Klre1","Klrk1")
p1_5 | DotPlot(GSE141784.int_CD8, features = markers) + RotatedAxis() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave("plots/cd8subset_integrated_markerids_dotplot.png", units="in", 
       width=10, height=5, dpi=320, scale = 2)


p1_5 | DotPlot(GSE141784.int_CD8, features = cytotoxic_markers) + RotatedAxis() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
ggsave("plots/cd8subset_integrated_cytotoxicmarkers_dotplot.png", units="in", 
       width=10, height=5, dpi=320, scale = 1)

VlnPlot(GSE141784.int_CD8, features = c(cytotoxic_markers),
        cols = rev(brewer.pal(n = 12, name = "Paired")))

FeaturePlot(GSE141784.int_CD8, reduction = "tsne",
            features = cytotoxic_markers,
            order = TRUE,
            #min.cutoff = "q10",
            split.by = "Sample",
            label = TRUE) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))


VlnPlot(GSE141784.int_CD8, features = c("Tcf7", "Tox"),
        cols = rev(brewer.pal(n = 12, name = "Paired")))

# Cell Annotation ---------------------------------------------------------

DefaultAssay(GSE141784.int_CD8) <- "integrated"

# Rename all identities
GSE141784.int_CD8 <- RenameIdents(object = GSE141784.int_CD8, 
                                  "0" = "Exhausted (Pdcd1+, Lag3+)",
                                  "1" = "Effector or memory",
                                  "2" = "Cycling Effector (Hmgb2+)",
                                  "3" = "Junb+",
                                  "4" = "Naive (Ccr7+, Lef1+)",
                                  "5" = "Memory (Il7r+, Cd28+)",
                                  "6" = "Cycling Effector (Hmgb2+)",
                                  "7" = "Early Effector (Xcl1+)")

p3 <- DimPlot(GSE141784.int_CD8, reduction = "tsne", order = TRUE) + 
  scale_color_brewer(type = "qual", palette = "Set1") +
  scale_alpha_continuous(range = 3/4, guide = F) +
  NoLegend() + DarkTheme()
p4 <- DimPlot(GSE141784.int_CD8, reduction = "umap", order = TRUE) + 
  scale_color_brewer(type = "qual", palette = "Set1") + DarkTheme()

png("plots/cd8subset_integrated_marker_ids_tsne_umap_dark.png", units="in", width=10, height=5, res=300)
cowplot::plot_grid(p3, p4)
dev.off()

save(GSE141784.int_CD8, file = "data/cd8_subset/GSE141784_integrated_markerids_cd8_subset.RData")

SaveH5Seurat(GSE141784.int_CD8,
             filename = paste0("GSE141784",
                               "_integrated_markerids_cd8_subset", ".h5Seurat"),
             overwrite = FALSE)
Convert("GSE141784_integrated_markerids_cd8_subset.h5Seurat", dest = "h5ad")
file.copy("GSE141784_integrated_markerids_cd8_subset.h5Seurat", "data/cd8_subset/h5_data/")
file.copy("GSE141784_integrated_markerids_cd8_subset.h5ad", "data/cd8_subset/h5_data/")

load(file="data/cd8_subset/GSE141784_integrated_markerids_cd8_subset.RData")

# Trajectory Analysis with **Slingshot** -----------------------------------------------------
pal <- rev(brewer.pal(n = 9, name = "Set1"))

# following the workshop from NBISWEDEN: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/slingshot/slingshot.html#basic_processing_with_seurat_pipeline

# Save the objects as separate matrices for input in slingshot
dimred <- GSE141784.int_CD8@reductions$tsne@cell.embeddings
clustering <- GSE141784.int_CD8$integrated_snn_res.0.75
counts <- as.matrix(GSE141784.int_CD8@assays$RNA@counts[GSE141784.int_CD8@assays$integrated@var.features,])

# Run default Slingshot lineage identification
set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clustering)
lineages

# Plot the lineages
png("plots/cd8subset_slingshot_default.png", units="in", width=10, height=5, res=300)
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
# legend("bottomleft", 
#        legend = levels(clustering),  
#        fill = pal[as.factor(levels(clustering))])
lines(lineages, lwd = 3, col = "black")
dev.off()


#Run Slingshot
set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        start.clus = "Naive (Ccr7+, Lef1+)") #define where to start the trajectories

lineages


#Plot the lineages
png("plots/cd8subset_slingshot_start4.png", units="in", width=10, height=5, res=300)
par(mfrow=c(1,2))
plot(dimred[,1:2], col = pal[clustering],  cex=.5,pch = 16)
for(i in levels(clustering)){ 
  text( mean(dimred[clustering==i,1]),
        mean(dimred[clustering==i,2]), labels = i,font = 2) }
plot(dimred, col = pal[clustering],  pch = 16)
lines(lineages, lwd = 3, col = 'black', show.constraints = TRUE)
dev.off()


curves <- getCurves(lineages,thresh = 0.01, stretch = 2, allow.breaks = FALSE, shrink = 0.99)
curves

png("plots/cd8subset_slingshot_curves.png", units="in", width=5, height=5, res=300)
plot(dimred, col = pal[clustering], asp = 1, pch = 16)
lines(curves, lwd = 3, col = "black")
dev.off()

# Analyzing gene expression along trajectories (DEGs) ---------------------
BiocParallel::register(BiocParallel::SerialParam())

library(tradeSeq)
set.seed(5)
png("plots/cd8subset_tradeseq_evalK.png", units="in", width=10, height=5, res=320)
icMat <- evaluateK(counts = counts, sds = curves, nGenes = 100,
                   k = 3:10, verbose = T)
dev.off()

print(icMat[1:2, ])


set.seed(7)
pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)
sce <- fitGAM(counts = counts, pseudotime = pseudotime,
              cellWeights = cellWeights,
              nknots = 8, verbose = TRUE)

table(SingleCellExperiment::rowData(sce)$tradeSeq$converged)

plotGeneCount(curves, counts, clusters = clustering, models = sce) + 
  theme_classic(base_size = 20)

assoRes <- associationTest(sce)
head(assoRes)

kbl(assoRes[1:20,], booktabs = T, caption="Association Test") %>% kable_styling(latex_options = "striped", full_width = FALSE) #%>% 



startRes <- startVsEndTest(sce)

oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[2]]
plotSmoothers(sce, counts, gene = sigGeneStart)

plotGeneCount(curves, counts, gene = sigGeneStart)


endRes <- diffEndTest(sce)


patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])


plotGeneCount(curve = curves, counts = counts,
              clusters = apply(slingClusterLabels(curves), 1, which.max),
              models = sce)

earlyDERes <- earlyDETest(sce, knots = c(1, 2))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])







