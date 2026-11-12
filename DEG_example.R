# cd /share/fsmresfiles/SpatialT/bone/mirae

library(Seurat)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(dplyr)
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
rm(list = ls())

# Blank_4W vs GP_4W
compare1 <- "NB_12W"
compare2 <- "GP_12W"
directory <- paste0("./DGEanalysis/", compare1, "_", compare2, "/")

# Check if the folder exists
if (!dir.exists(directory)) {
  # Create the folder if it doesn't exist
  dir.create(directory, recursive = TRUE)
  message("Folder created: ", directory)
} else {
  message("Folder already exists: ", directory)
}

so1 <- readRDS(paste0(compare1, "/Seurat.RDS"))
so2 <- readRDS(paste0(compare2, "/Seurat.RDS"))

common.features <- intersect(rownames(so1), rownames(so2))
DefaultAssay(so1) <- "Spatial"
DefaultAssay(so2) <- "Spatial"
so1[["Spatial"]] <- subset(so1[["Spatial"]], features = common.features)
so2[["Spatial"]] <- subset(so2[["Spatial"]], features = common.features)

so1 <- SCTransform(object = so1, assay = "Spatial", return.only.var.genes = FALSE)
so2 <- SCTransform(object = so2, assay = "Spatial", return.only.var.genes = FALSE)
so1[["ROI"]] <- compare1
so2[["ROI"]] <- compare2
so <- merge(so1, y = so2, add.cell.ids = c(compare1, compare2))
table(so@meta.data$ROI)

annotations <- so[["ROI"]]$ROI
names(annotations) <- rownames(so[["ROI"]])
Idents(so) <- annotations
so <- PrepSCTFindMarkers(so, assay = "SCT", verbose = TRUE)
saveRDS(so, paste0(directory, compare1, "_", compare2, "_merged_prepSCT.RDS"))

results <- FindMarkers(so, ident.1 = compare1, ident.2 = compare2, assay = "SCT", logfc.threshold = 0.1, min.pct = 0.01, verbose = FALSE)
results$gene <- rownames(results)
results$group1 <- compare1
results$group2 <- compare2
results <- results[results$gene %in% common.features,]
fwrite(results, paste0(directory, compare1, "_", compare2, "_DGE.tsv"), sep = "\t")

results <- results[results$p_val_adj < 0.05,]
results <- results[,c("gene", "group1", "group2", "avg_log2FC", "p_val", "p_val_adj")]
results <- results[order(results$avg_log2FC),]

#add bone-genesets cluster
gset <- read.csv("bone_genesets.csv")
gset <- gset %>% mutate_all(tolower)
colnames(gset) <- c("Bone", "Osteoblasts", "Osteoclasts", "Macrophage", "BoneMarrow")
result <- results
result$gene_cluster <- NA

result <- result %>%
  rowwise() %>%
  mutate(gene_cluster = {
    gene <- tolower(gene) # Convert gene to lowercase for matching
    matches <- colnames(gset)[apply(gset, 2, function(col) gene %in% col)]
    if(length(matches) > 1) {
      paste(matches, collapse = "_")
    } else if(length(matches) == 1) {
      matches
    } else {
      NA
    }
  })

fwrite(result, paste0(directory, compare1, "_", compare2, "_DGE_sig.tsv"), sep = "\t")

# so <- ScaleData(so, features = rownames(so))
# pdf(paste0(compare1, "_", compare2, "_heatmap.pdf"))
# DoHeatmap(so, features = results$gene, slot = "data") + scale_fill_gradientn(colours = SpatialColors(n = 100))
# dev.off()

#get list of expressed genes ordered from high to low
ex_gene <- AggregateExpression(object = so, assays = "Spatial")
ex_gene <- as.data.frame(ex_gene) %>%
  arrange(desc(Spatial.NB.12W))
ex_list <- as.data.frame(ex_gene[1:200, 1])
rownames(ex_list) <- rownames(ex_gene)[1:200]
ex_list$gene <- rownames(ex_list)
fwrite(ex_list, paste0(directory, compare1, "_", compare2, "_highExpress_", compare1, ".tsv"), sep = "\t")

ex_gene <- as.data.frame(ex_gene) %>%
  arrange(desc(Spatial.Blank.12W))
ex_list <- as.data.frame(ex_gene[1:200, 2])
rownames(ex_list) <- rownames(ex_gene)[1:200]
ex_list$gene <- rownames(ex_list)
fwrite(ex_list, paste0(directory, compare1, "_", compare2, "_highExpress_", compare2, ".tsv"), sep = "\t")

#####################################################################################
# three groups comparison
compare1 <- "Blank_12W"
compare2 <- "GP_12W"
compare3 <- "NB_12W"
directory <- paste0("./DGEanalysis/", compare1, "_", compare2, "_", compare3, "/")

# Check if the folder exists
if (!dir.exists(directory)) {
  # Create the folder if it doesn't exist
  dir.create(directory, recursive = TRUE)
  message("Folder created: ", directory)
} else {
  message("Folder already exists: ", directory)
}

so1 <- readRDS(paste0(compare1, "/Seurat.RDS"))
so2 <- readRDS(paste0(compare2, "/Seurat.RDS"))
so3 <- readRDS(paste0(compare3, "/Seurat.RDS"))

common.features <- intersect(rownames(so1), rownames(so2))
common.features <- intersect(common.features, rownames(so3))

DefaultAssay(so1) <- "Spatial"
DefaultAssay(so2) <- "Spatial"
DefaultAssay(so3) <- "Spatial"
so1[["Spatial"]] <- subset(so1[["Spatial"]], features = common.features)
so2[["Spatial"]] <- subset(so2[["Spatial"]], features = common.features)
so3[["Spatial"]] <- subset(so3[["Spatial"]], features = common.features)

so1 <- SCTransform(object = so1, assay = "Spatial", return.only.var.genes = FALSE)
so2 <- SCTransform(object = so2, assay = "Spatial", return.only.var.genes = FALSE)
so3 <- SCTransform(object = so3, assay = "Spatial", return.only.var.genes = FALSE)

so1[["ROI"]] <- compare1
so2[["ROI"]] <- compare2
so3[["ROI"]] <- compare3

so <- merge(so1, y = c(so2, so3), add.cell.ids = c(compare1, compare2, compare3))
table(so@meta.data$ROI)

annotations <- so[["ROI"]]$ROI
names(annotations) <- rownames(so[["ROI"]])
Idents(so) <- annotations
so <- PrepSCTFindMarkers(so, assay = "SCT", verbose = TRUE)
saveRDS(so, paste0(directory, compare1, "_", compare2, "_", compare3, "_merged_prepSCT.RDS"))

results <- FindAllMarkers(so, assay = "SCT", logfc.threshold = 0.1, min.pct = 0.01, verbose = FALSE)
results$gene <- rownames(results)
results$group1 <- compare1
results$group2 <- compare2
results$group3 <- compare3
results <- results[results$gene %in% common.features,]
fwrite(results, paste0(directory, compare1, "_", compare2, "_", compare3, "_DGE.tsv"), sep = "\t")

results <- results[results$p_val_adj < 0.05,]
results <- results[,c("gene", "group1", "group2", "group3", "avg_log2FC", "p_val", "p_val_adj", "cluster")]
results <- results%>%
  arrange(cluster, avg_log2FC)

#add bone-genesets cluster
gset <- read.csv("bone_genesets.csv")
gset <- gset %>% mutate_all(tolower)
colnames(gset) <- c("Bone", "Osteoblasts", "Osteoclasts", "Macrophage", "BoneMarrow")
result <- results
result$gene_cluster <- NA

result <- result %>%
  rowwise() %>%
  mutate(gene_cluster = {
    gene <- tolower(gene) # Convert gene to lowercase for matching
    matches <- colnames(gset)[apply(gset, 2, function(col) gene %in% col)]
    if(length(matches) > 1) {
      paste(matches, collapse = "_")
    } else if(length(matches) == 1) {
      matches
    } else {
      NA
    }
  })

fwrite(result, paste0(directory, compare1, "_", compare2, "_", compare3, "_DGE_sig.tsv"), sep = "\t")
