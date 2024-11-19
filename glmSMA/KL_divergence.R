#KL scores
#input: brain 10x Visium data, cell prediction results, cell type
#ouput: KL scores

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


KL = function(brain_,brain_sc,Cell_Type,st_reference_atlas,spot_distribution_target){
#input

selected_cell_type = Cell_Type

#plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
#plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
#wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

allen_reference <- st_reference_atlas

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE)

RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)


cortex = brain

# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
#DimPlot(allen_reference, group.by = "subclass", label = TRUE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

DefaultAssay(cortex) <- "predictions"


spot_cell_type_prediction = cortex@assays$predictions@data
selected_index = which(rownames(spot_cell_type_prediction) == selected_cell_type)

spot_distribution = spot_cell_type_prediction[selected_index,] / sum(spot_cell_type_prediction[selected_index,])

c <- rbind(spot_distribution_target,spot_distribution)

KL = KL(c, unit='log')

return(KL)
}
