# Pre-processing ------------------------------------------------------------

# Setting up the data for merging
getwd()
setwd("~/Documents/GSE141784/GSE141784_RAW/")
list.dirs()
list.files()

dir.create("NOD_4wk")
dir.create("NOD_8wk")
dir.create("NOD_15wk")


nod4w <- list.files(pattern = "NOD_4w", 
                    full.names = FALSE, 
                    recursive = TRUE)
nod8w <- list.files(pattern = "NOD_8w", 
                    full.names = FALSE,
                    recursive = TRUE)
nod15w <- list.files(pattern = "NOD_15w", 
                     full.names = FALSE,
                     recursive = TRUE)


file.copy(nod4w[1:3], "NOD_4wk")
file.copy(nod8w[1:3], "NOD_8wk")
file.copy(nod15w, "NOD_15wk")

setwd("~/Documents/GSE141784/GSE141784_RAW/NOD_4wk")
file.rename(list.files(), 
            gsub("GSM4213196_NOD_4w_2734_", "", list.files()))
file.rename("genes.tsv.gz", gsub("genes.tsv.gz","features.tsv.gz", "genes.tsv.gz"))

setwd("~/Documents/GSE141784/GSE141784_RAW/NOD_8wk")
file.rename(list.files(),
            gsub("GSM4213197_NOD_8w_2734_", "", list.files()))
file.rename("genes.tsv.gz", gsub("genes.tsv.gz","features.tsv.gz", "genes.tsv.gz"))


setwd("~/Documents/GSE141784/GSE141784_RAW/NOD_15wk")
file.rename(list.files(), 
            gsub("GSM4213198_NOD_15w_2734_", "", list.files()))
file.rename("genes.tsv.gz", gsub("genes.tsv.gz","features.tsv.gz", "genes.tsv.gz"))


setwd("~/Documents/RStudio/GSE141784_reanalysis")


# Seurat Setup ------------------------------------------------------------

# Loading in the packages
library(Seurat)
library(knitr)
library(kableExtra)
library(ggplot2)

# Loading the CellRanger Matrix Data and creating the base Seurat object.

dataset_loc <- "~/Documents/GSE141784/GSE141784_RAW/"
ids <- c("NOD_4wk", "NOD_8wk", "NOD_15wk") # the data directory names just created above

d10x.data <- sapply(ids, function(i){
  d10x <- Read10X(file.path(dataset_loc,i))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
  d10x
})

experiment.data <- do.call("cbind", d10x.data)

experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "scRNA test",
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "\\-")

# saveRDS(experiment.aggregate,"GSE141784_experiment_aggregate.rds") # dataset is ~ 1GB

experiment.aggregate <- readRDS("GSE141784_experiment_aggregate.rds")


# Annotation --------------------------------------------------------------


# The percentage of reads that map to the mitochondrial genome
experiment.aggregate$percent.mito <- PercentageFeatureSet(experiment.aggregate, pattern = "^mt-")

# Calculate cell cycle, add to meta data
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

head(experiment.aggregate$percent.mito)

# Convert to matrix for use in cycle
mat <- as.matrix(GetAssayData(experiment.aggregate))

# Convert rownames to ENSEMBL IDs, Using biomaRt

# ensembl<- biomaRt::useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
# anno <- biomaRt::getBM(values=rownames(mat), attributes=c("mgi_symbol","ensembl_gene_id") , filters= "mgi_symbol"  ,mart=ensembl)
# saveRDS(anno, "biomart_export_Nov2020.rds")
anno <- readRDS(file="biomart_export_Nov2020.rds")

ord <- match(rownames(mat), anno$mgi_symbol) # use anno$mgi_symbol if via biomaRt
rownames(mat) <- anno$ensembl_gene_id[ord] # use anno$ensembl_gene_id if via biomaRt
drop <- which(is.na(rownames(mat)))
mat <- mat[-drop,]
cycles <- scran::cyclone(mat, pairs=mm.pairs) # since the matrix is very large this takes a while to run.
tmp <- data.frame(cell.cycle = cycles$phases)
rownames(tmp) <- colnames(mat)
experiment.aggregate <- AddMetaData(experiment.aggregate, tmp)

save(experiment.aggregate,file="GSE141784_original_seurat_object.RData")

# Plotting ----------------------------------------------------------------

load(file="GSE141784_original_seurat_object.RData")

# Generating a subdirectory for plots
ifelse(!dir.exists((paste0(getwd(),"/plots"))), 
       dir.create(paste0(getwd(),"/plots")), FALSE) 

# 5% quantiles for number of genes per cell per sample

do.call("cbind", tapply(experiment.aggregate$nFeature_RNA, Idents(experiment.aggregate),quantile,probs=seq(0,1,0.05)))

# 5% quantiles for number of UMI per cell per sample

do.call("cbind", tapply(experiment.aggregate$nCount_RNA, Idents(experiment.aggregate),quantile,probs=seq(0,1,0.05)))

# 5% quantiles for number of mitochondrial percentage per cell per sample

round(do.call("cbind", tapply(experiment.aggregate$percent.mito, Idents(experiment.aggregate),quantile,probs=seq(0,1,0.05))), digits = 3)

# Table of cell cycle (Not working)

table(experiment.aggregate@meta.data$cell.cycle) %>% kable(caption = "Number of Cells in each Cell Cycle Stage", col.names = c("Stage", "Count"), align = "c") %>% kable_styling() %>%
  save_kable(file = "plots/cell_cycle_table.pdf",
           density = 320)


# Plot the number of cells each gene is represented by

png("plots/gene_representation.png", units="in", width=5, height=5, res=300)
plot(sort(Matrix::rowSums(GetAssayData(experiment.aggregate) >= 3)) , xlab="gene rank", ylab="number of cells", main="Cells per genes (reads/gene >= 3 )")
dev.off()

png("plots/exploratory_violinplots.png", units="in", width=8, height=8, res=300)
VlnPlot(
  experiment.aggregate,
  features = c("nFeature_RNA", "nCount_RNA","percent.mito"),
  ncol = 1, pt.size = 0.3)
dev.off()

png("plots/exploratory_cellcycle_scatter.png", units="in", width=5, height=5, res=300)
FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "percent.mito")
dev.off()

png("plots/exploratory_RNA_scatter.png", units="in", width=5, height=5, res=300)
FeatureScatter(
  experiment.aggregate, "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5)
dev.off()

# Using the information from the exploratory plots to subset the dataset 

experiment.aggregate <- subset(experiment.aggregate, percent.mito <= 10)

experiment.aggregate <- subset(experiment.aggregate, nCount_RNA >= 500 & nCount_RNA <= 30000)

experiment.aggregate

save(experiment.aggregate,file="GSE141784_sub_prenorm_seurat_object.RData")

# Filtering/Normalization -------------------------------------------------

load(file="GSE141784_sub_prenorm_seurat_object.RData")









