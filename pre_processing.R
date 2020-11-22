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

saveRDS(experiment.aggregate,"GSE141784_experiment_aggregate.rds") # dataset is ~ 1GB





