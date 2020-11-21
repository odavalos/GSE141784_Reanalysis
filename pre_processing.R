# Pre-processing


# Setting up the data for merging
getwd()
setwd("~/Documents/GSE141784/GSE141784_RAW/")
list.dirs()
list.files()
nod4w <- list.files(pattern = "NOD_4w")
nod8w <- list.files(pattern = "NOD_8w")
nod15w <- list.files(pattern = "NOD_15w")

dir.create(getwd())

