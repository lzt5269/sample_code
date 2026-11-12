library(Seurat)
library(data.table)
library(stringr)
library(dplyr)
setwd("/projects/p32119/tisreg/ref/GSE128423")
##################################################
### Read in the metadata table
metadata <- read.csv("GSE128423_metadata.csv")
metadata <- metadata %>%
  filter(!grepl("std7|std8", NAME))
metadata$NAME <- paste0(metadata$NAME, "-1")

### Create Seurat objects from MEX files and merge them
subsetList <- function(myList, elementNames) {
  sapply(elementNames, FUN=function(x) myList[[x]])
}
so_ls <- list()
files <- c("std1", "std2", "std3", "std4")
for (i in 1:length(files)) {
  sample_id <- files[i]
  curr_so <- CreateSeuratObject(counts = Read10X(data.dir = paste0("./", sample_id)), project = sample_id)
  curr_so[["percent.mt"]] <- PercentageFeatureSet(curr_so, pattern = "^mt-")
  so_ls[[i]] <- curr_so
}
so <- merge(x = so_ls[[1]], y = c(subsetList(so_ls, seq(2, length(so_ls)))), add.cell.ids = files, project = "Bone")

### Filter the metadata table by keeping only the cells in the Seurat object
cells <- colnames(so)
metadata_filtered <- metadata[metadata$NAME %in% cells,] # 7137 cells only
### Check available cell type annotations
table(metadata$Cluster, useNA = "always")

metadata <- metadata %>%
  mutate(label = ifelse(Cluster == 0 | Cluster == 6 | Cluster == 11, 'BMEC',
                        ifelse(Cluster == 1, 'MSC',
                               ifelse(Cluster == 2 | Cluster == 10 | Cluster == 13 | Cluster == 17, 'Chondrocyte',
                                      ifelse(Cluster == 3 | Cluster == 5 | Cluster == 9 | Cluster == 15 | Cluster == 16, 'Fibroblast',
                                             ifelse(Cluster == 4 | Cluster == 8, 'Chondrocyte-Osteoblast-Fibroblast',
                                                    ifelse(Cluster == 7 | Cluster == 8, 'OLC',
                                                           ifelse(Cluster == 12, 'Pericyte', NA)
                                                    )))))))


### Subset the Seurat object with the cells available in the metadata table
so[["cell_id"]] <- colnames(so)
so_less <- subset(so, subset = cell_id %in% metadata_filtered$NAME)

### Check if the cell IDs match
if (sum(so_less@meta.data[["cell_id"]] %in% metadata$NAME) == nrow(metadata_filtered)) {
  cat("\nCell IDs match")
} else {
  cat("\nCell IDs do not match -- need to sort or parse the cell IDs")
}

### Check if the cell IDs match
if (sum(so_less[["cell_id"]] == metadata_filtered$NAME) == nrow(metadata_filtered)) {
  cat("\nCell IDs match")
} else {
  cat("\nCell IDs do not match -- need to sort or parse the cell IDs")
}

### Reorder labels so that the cell IDs match
meta <- so_less@meta.data
meta$i <- 1:nrow(meta)
colnames(meta)[1] <- "cell_id"
metadata_filtered$cell_id <- metadata_filtered$NAME
meta <- merge(meta, metadata_filtered, by = "cell_id")
meta <- meta[order(meta$i),]
sum(so_less[["cell_id"]] == meta$cell_id) == nrow(meta) #!!! This should be true

### Add the cell type annotation
so_less[["label"]] <- meta$label
so_less[["label_finer"]] <- meta$Subcluster # Optional

### Double check the metadata stored in the Seurat object
str(so_less@meta.data)

### Save Seurat object
saveRDS(so_less, "Seurat.RDS")

###############################################################################
### Remove doublets
so_less = NormalizeData(so_less)
so_less = FindVariableFeatures(so_less, verbose = F)
so_less = ScaleData(so_less, vars.to.regress = c("nFeature_RNA", "percent.mt"),
                      verbose = F)
so_less = RunPCA(so_less, verbose = F, npcs = 20)
so_less = RunUMAP(so_less, dims = 1:10, verbose = F)

# Can run parameter optimization with paramSweep

# sweep.res <- paramSweep_v3(data.filt) sweep.stats <-
# summarizeSweep(sweep.res, GT = FALSE) bcmvn <- find.pK(sweep.stats)
# barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# define the expected number of doublet cellscells.
nExp <- round(ncol(so_less) * 0.058)  # expect 4% doublets
so_less <- doubletFinder_v3(so_less, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(so_less@meta.data)[grepl("DF.classification", colnames(so_less@meta.data))]


cowplot::plot_grid(ncol = 2, DimPlot(so_less, group.by = "orig.ident") + NoAxes(),
                   DimPlot(so_less, group.by = DF.name) + NoAxes())

VlnPlot(so_less, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

so_less = so_less[, so_less@meta.data[, DF.name] == "Singlet"]
dim(so_less)


