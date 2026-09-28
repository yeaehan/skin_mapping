

saveRDS(chat, "data/chat.rds")
# 1. Install ALRA (used for imputation in single-cell data)
devtools::install_github("KlugerLab/ALRA")

# 2. Install MERINGUE (used for spatial transcriptomics analysis)
devtools::install_github("JEFworks-Lab/MERINGUE")
devtools::install_github("jinworks/SpatialCellChat")
install.packages("reticulate")

library(reticulate)
install_python(version = "3.10")

library(reticulate)
reticulate::py_config()

Sys.setenv(
  RETICULATE_PYTHON = "~/.virtualenvs/r-reticulate/bin/python"
)

reticulate::py_install(
  "numpy",
  pip = TRUE,
  pip_options = c(
    "--trusted-host", "pypi.org",
    "--trusted-host", "files.pythonhosted.org"
  )
)

library(SpatialCellChat)
library(Matrix)
library(patchwork)
options(stringsAsFactors = F)

# use `SetEnvironment` to do parallel computation and get Python back-end
setEnvironment(workers = 2,future.globals.maxSize = 100000*1024^2, conda_env = NULL)
#> ℹ Future strategy in use: `parallel with 2 workers`
# setEnvironment(2,future.globals.maxSize = 100000*1024^2,conda_env = "/home/xzcheng/miniconda3/envs/Chat")


data.input <- GetAssayData(hsskin.ST, assay = "SCT", layer = "data") 

# Extract Metadata (must contain your cell type annotations)
meta <- hsskin.ST@meta.data 
# Ensure your cell types are set as a factor
meta$cell_type_major <- as.factor(meta$cell_type_major) 

# Extract Spatial Coordinates 
# Seurat V5 stores Xenium coordinates in the images/FOV slot. 
# 1. Create an empty data frame to hold all coordinates
all_spatial_locs <- data.frame()

# 2. Loop through every image/sample stored in your Seurat object
for (img_name in Images(hsskin.ST)) {
  
  # Extract coordinates for this specific FOV
  coords <- GetTissueCoordinates(hsskin.ST[[img_name]])
  
  # In Seurat V5, 'coords' will have columns: x, y, and cell
  # We set the row names to the actual cell barcodes
  rownames(coords) <- coords$cell
  
  # Append just the x and y columns to our master table
  all_spatial_locs <- rbind(all_spatial_locs, coords[, c("x", "y")])
}

# 3. Reorder and subset the master table to perfectly match the Seurat object
# This ensures 1:1 matching and drops any cells that might have been filtered out
spatial.locs <- all_spatial_locs[colnames(hsskin.ST), c("x", "y")]
# Define Spatial Factors for Xenium
# Xenium coordinates are natively in micrometers (um), so ratio = 1.
# Average human cell diameter is ~10um, so tolerance (tol) is half of that (5um).
spatial.factors <- list(ratio = 1, tol = 5)

# ----------------===================================================---------------- #
# STEP 2: CREATE SPATIALCELLCHAT OBJECT
# ----------------===================================================---------------- #
# Use the new createSpatialCellChat function (bypasses manual data/meta extraction)
# 1. Manually extract the normalized data matrix using Seurat V5 'layer' syntax
data.input <- GetAssayData(hsskin.ST, assay = "SCT", layer = "data")

# 2. Extract the metadata
meta <- hsskin.ST@meta.data

# 3. Create the SpatialCellChat object using the raw matrix and metadata
# This completely bypasses the broken Seurat V4 code inside the package
chat <- createSpatialCellChat(
  object = data.input,
  meta = meta,
  group.by = "cell_type_major",
  datatype = "spatial",
  coordinates = spatial.locs,
  spatial.factors = spatial.factors
)

# Set the Ligand-Receptor database
CellChatDB <- CellChatDB.human # Use CellChatDB.mouse if working with mouse
# Subset to relevant signaling types (removes non-protein signaling)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), non_protein = F)
chat@DB <- CellChatDB.use

# ----------------===================================================---------------- #
# STEP 3: PREPROCESSING
# ----------------===================================================---------------- #
chat <- subsetData(chat) 
chat <- preProcessing(chat) 

# 1. Force SpatialCellChat to use ALL genes available in your Xenium panel
chat@var.features <- list(features = rownames(chat@data.signaling))

# 2. Force it to test ALL possible ligand-receptor pairs in the database 
# (it will automatically ignore pairs where the genes aren't in your 500-plex panel)
chat@LR$LRsig <- chat@DB$interaction

# 3. Skip the 'identifyOverExpressedInteractions' function entirely.
# Jump straight to computing the probabilities!


# # Identify spatially variable genes using the MERINGUE method (from the tutorial)
chat2 <- identifyOverExpressedGenes(
  chat,
  selection.method = "meringue",
  do.grid = F
)
# 
# print(paste("Number of variable features:", length(chat@var.features)))
# print(paste("Genes in DB overlap:", length(intersect(rownames(chat@data.signaling), CellChatDB.use$geneInfo$Symbol))))
# graphics.off()
# closeAllConnections()
# 
# # 2. progressr 패키지가 있다면 명시적으로 꺼주기
# if ("progressr" %in% loadedNamespaces()) {
#   progressr::handlers(global = FALSE)
# }
# chat <- identifyOverExpressedInteractions(chat, variable.both = F)

# ----------------===================================================---------------- #
# STEP 4: COMPUTE SPATIAL COMMUNICATION PROBABILITIES
# ----------------===================================================---------------- #
# Cell-level inference 
# scale.distance = 1 because Xenium is already in um.
# interaction.range = 25 (max diffusion distance in um for paracrine signaling)
# contact.range = 10 (max physical distance in um for juxtacrine/contact signaling)

progressr::handlers(global = FALSE)

genes.use <- rownames(chat@data.signaling)

lr <- chat@LR$LRsig

keep <- sapply(seq_len(nrow(lr)), function(i) {
  lig <- unlist(strsplit(lr$ligand[i], "_"))
  rec <- unlist(strsplit(lr$receptor[i], "_"))
  
  all(lig %in% genes.use) && all(rec %in% genes.use)
})

chat@LR$LRsig <- lr[keep, ]

chat <- computeCommunProb(
  chat,
  distance.use = TRUE,
  scale.distance = 700,
  contact.dependent = TRUE,
  interaction.range = 25,
  contact.range = 10
)



# Filter out non-significant/low-count interactions
chat <- filterProbability(chat)
chat <- filterCommunication(chat, min.cells = 10, min.links = 10, min.cells.sr = 10)

# Group-level inference (Use computeAvgCommunProb for single-cell resolution like Xenium)
chat <- computeAvgCommunProb(chat)

# Compute pathways and aggregate the network
chat <- computeCommunProbPathway(chat)
chat <- aggregateNet(chat)

# Compute Centrality (Required for advanced plotting)
chat <- netAnalysis_computeCentrality(chat, slot.name = "net", do.group = F, degree.only = T)
chat <- netAnalysis_computeCentrality(chat, slot.name = "netP", do.group = F, degree.only = T)

# ----------------===================================================---------------- #
# STEP 5: DOWNSTREAM VISUALIZATION (Examples)
# ----------------===================================================---------------- #

# 1. Circle plot of aggregated interactions between cell types
netVisual_aggregate(chat, layout = "circle")

# 2. Extract significant pathways found in your dataset
print(chat@netP$pathways)
pathway.show <- chat@netP$pathways[1] # Change to a pathway of interest like "WNT" or "MIF"

# 3. Spatial Feature Plot of the pathway
spatialFeaturePlot(chat, signaling = pathway.show, do.group = FALSE, color.heatmap = "Spectral", point.size = 1.3)

# 4. Show incoming vs outgoing hotspot scores
gg1 <- spatialVisual_scoring(chat, signaling = pathway.show, slot.name = "netP", do.group = FALSE, measure = c("outdeg"), color.heatmap = "Blues", point.size = 1)
gg2 <- spatialVisual_scoring(chat, signaling = pathway.show, slot.name = "netP", do.group = FALSE, measure = c("indeg"), color.heatmap = "Reds", point.size = 1)
patchwork::wrap_plots(gg1, gg2, ncol = 2)


# ----------------===================================================---------------- #
# No Database, just use the ligand-receptor pairs in Xenium
# ----------------===================================================---------------- #


library(SpatialCellChat)
library(future)
plan("sequential")
options(future.globals.maxSize = 100000 * 1024^2)

# 1. Create the Object
chat <- createSpatialCellChat(
  object = data.input,
  meta = meta,
  group.by = "cell_type_major",
  datatype = "spatial",
  coordinates = spatial.locs,
  spatial.factors = spatial.factors
)

# 2. Extract the 1-to-1 pairs from the Human DB
custom_DB <- CellChatDB.human

# 2. Get your exact Xenium genes
xenium_genes <- rownames(data.input)

# 3. Find interactions where BOTH the ligand and receptor exact gene symbols are in your data
# (This brilliantly auto-drops complexes like "TGFbR1_R2" because they aren't valid gene symbols!)
valid_ligands <- custom_DB$interaction$ligand %in% xenium_genes
valid_receptors <- custom_DB$interaction$receptor %in% xenium_genes

# 4. Overwrite ONLY the interaction table in our custom DB
custom_DB$interaction <- custom_DB$interaction[valid_ligands & valid_receptors, ]

# 5. Assign it to the chat object
chat@DB <- custom_DB

# 6. Now subsetData will find all the internal tables it needs
chat <- subsetData(chat)
chat <- preProcessing(chat)

# Force CellChat to use your perfectly filtered list
chat@LR$LRsig <- chat@DB$interaction

# 7. Compute Probabilities
chat <- computeCommunProb(
  chat,
  distance.use = TRUE, 
  scale.distance = 7.5,    
  contact.dependent = TRUE,
  interaction.range = 25, 
  contact.range = 10 
)