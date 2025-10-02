# load libraries
library(Seurat)
library(tidyverse)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

# read
patient1_2hr <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient1_2hr')
patient2_2hr <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient2_2hr')
patient3_2hr <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient3_2hr')
patient4_2hr <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient4_2hr')
patient5_2hr <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient5_2hr')
patient6_2hr <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient6_2hr')
patient1_cit <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient1_cit')
patient2_cit <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient2_cit')
patient3_cit <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient3_cit')
patient4_cit <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient4_cit')
patient5_cit <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient5_cit')
patient6_cit <- Read10X(
  data.dir = 'C:/Users/tanta/Downloads/GSE220797_Processed_V2/patient6_cit')

# create seurat object
patient1_2hr <- CreateSeuratObject(counts = patient1_2hr, project = "2hr")
patient2_2hr <- CreateSeuratObject(counts = patient2_2hr, project = "2hr")
patient3_2hr <- CreateSeuratObject(counts = patient3_2hr, project = "2hr")
patient4_2hr <- CreateSeuratObject(counts = patient4_2hr, project = "2hr")
patient5_2hr <- CreateSeuratObject(counts = patient5_2hr, project = "2hr")
patient6_2hr <- CreateSeuratObject(counts = patient6_2hr, project = "2hr")
patient1_cit <- CreateSeuratObject(counts = patient1_cit, project = "cit")
patient2_cit <- CreateSeuratObject(counts = patient2_cit, project = "cit")
patient3_cit <- CreateSeuratObject(counts = patient3_cit, project = "cit")
patient4_cit <- CreateSeuratObject(counts = patient4_cit, project = "cit")
patient5_cit <- CreateSeuratObject(counts = patient5_cit, project = "cit")
patient6_cit <- CreateSeuratObject(counts = patient6_cit, project = "cit")

# create list of object names
object_Names <- c(
  "patient1_2hr", "patient2_2hr", "patient3_2hr", "patient4_2hr",
  "patient5_2hr", "patient6_2hr", "patient1_cit", "patient2_cit",
  "patient3_cit", "patient4_cit", "patient5_cit", "patient6_cit"
)

# subset half the cells from each sample
for (name in object_Names) {
  
  # retrieve object
  obj <- get(name)
  
  # randomly sample half the cells
  cell_Names <- colnames(obj)
  sampled_cells <- sample(cell_Names, length(cell_Names)/2)
  obj_subset <- subset(obj, cells = sampled_cells)
  
  # assign filtered object back to original name and delete temporary variables
  assign(name, obj_subset, envir = .GlobalEnv)
  remove(obj)
  remove(cell_Names)
  remove(sampled_cells)
  remove(obj_subset)
}

# merge datasets
merged_seurat <- merge(patient1_2hr, y = c(
  patient1_cit, patient2_2hr, patient2_cit, patient3_2hr, patient3_cit,
  patient4_2hr, patient4_cit, patient5_2hr, patient5_cit, patient6_2hr,
  patient6_cit), add.cell.ids = ls()[1:12], project = "merged dataset"
)

# subset merged object
cell_names <- colnames(merged_seurat)
sampled_cells <- sample(cell_names, length(cell_names)/2)
merged_seurat_subset <- subset(merged_seurat, cells = sampled_cells)

# manually filter out low quality cells based on distribution of data
VlnPlot(patient1_2hr, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient1_2hr <- subset(
  patient1_2hr, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & nCount_RNA < 75000)

VlnPlot(patient2_2hr, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient2_2hr <- subset(
  patient2_2hr, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 50000)

VlnPlot(patient3_2hr, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient3_2hr <- subset(
  patient3_2hr, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 50000)

VlnPlot(patient4_2hr, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient4_2hr <- subset(
  patient4_2hr, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 30000)

VlnPlot(patient5_2hr, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient5_2hr <- subset(
  patient5_2hr, subset = nFeature_RNA > 200 & nFeature_RNA < 7800 & nCount_RNA < 35000)

VlnPlot(patient6_2hr, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient6_2hr <- subset(
  patient6_2hr, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 30000)

VlnPlot(patient1_cit, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient1_cit <- subset(
  patient1_cit, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA < 60000)

VlnPlot(patient2_cit, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient2_cit <- subset(
  patient2_cit, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 40000)

VlnPlot(patient3_cit, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient3_cit <- subset(
  patient3_cit, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 50000)

VlnPlot(patient4_cit, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient4_cit <- subset(
  patient4_cit, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & nCount_RNA < 30000)

VlnPlot(patient5_cit, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient5_cit <- subset(
  patient5_cit, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA < 50000)

VlnPlot(patient6_cit, features = c(
  "nFeature_RNA", "nCount_RNA", "pct_counts_Mito"), ncol = 3)
patient6_cit <- subset(
  patient6_cit, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 28000)

# merged seurat QC
VlnPlot(merged_subset, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
merged_subset <- subset(
  merged_subset, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 40000
)
merged_subset <- qc_MT(merged_subset)


# function that filters out cells with high mitochondrial content
qc_MT <- function(obj, mito_pattern = "^MT-", mads_thresh = 4, 
                  hard_thresh = 50) {
  obj <- PercentageFeatureSet(obj, pattern = mito_pattern, 
                              col.name = "pct_counts_Mito")
  mito_thresh <- median(obj$pct_counts_Mito) + mad(obj$pct_counts_Mito) * 
    mads_thresh
  drop_mito <- obj$pct_counts_Mito > mito_thresh | obj$pct_counts_Mito > 
    hard_thresh
  obj <- obj[, !drop_mito]
  return(obj)
}

# filter out cells with high MT content
for (name in object_Names) {
  # print the sample we are on
  print(name)
  
  # retrieve object from name
  obj <- get(name)
  
  # apply function qc_MT
  obj_filtered <- qc_MT(obj)
  
  # assign filtered object back to original name and delete temporary variables
  assign(name, obj_filtered, envir = .GlobalEnv)
  remove(obj)
  remove(obj_filtered)
}

# list of seurat objects
obj_List <- list(
  patient1_2hr = patient1_2hr,
  patient2_2hr = patient2_2hr,
  patient3_2hr = patient3_2hr,
  patient4_2hr = patient4_2hr,
  patient5_2hr = patient5_2hr,
  patient6_2hr = patient6_2hr,
  patient1_cit = patient1_cit,
  patient2_cit = patient2_cit,
  patient3_cit = patient3_cit,
  patient4_cit = patient4_cit,
  patient5_cit = patient5_cit,
  patient6_cit = patient6_cit
)

# clean up workspace
rm(
  patient1_2hr, patient2_2hr, patient3_2hr,
  patient4_2hr, patient5_2hr, patient6_2hr,
  patient1_cit, patient2_cit, patient3_cit,
  patient4_cit, patient5_cit, patient6_cit
)

# standard pre-processing (prepare for integration)
for (i in 1:length(obj_List)) {
  obj_List[[i]] <- NormalizeData(object = obj_List[[i]])
  obj_List[[i]] <- FindVariableFeatures(
    object = obj_List[[i]], selection.method = "vst", nfeatures = 2000)
}

# save object list as RDS file
save_dir <- paste0(
  "C:/Users/tanta/OneDrive - University of Toronto/Desktop/Bioinformatics_Training"
)
saveRDS(integrated_data, file.path(save_dir, "integratedData_BeforeDoublet.rds"))

# load object list 
integrated_data <- readRDS(paste0(save_dir, "/integratedData.rds"))

# select features that are repeatedly variable across datasets for integration
features = SelectIntegrationFeatures(object.list = obj_List)
anchors = FindIntegrationAnchors(object.list = obj_List, anchor.features = features)

# create an 'integrated' data assay
integrated_Data <- IntegrateData(anchorset = anchors)

# change default assay from "RNA" to "integrated"
DefaultAssay(integrated_data) <- "integrated"

# run standard pre-processing workflow
integrated_data <- ScaleData(integrated_data)
integrated_data <- RunPCA(integrated_data, nfeatures.print = 10)

# use elbow plot to find optimal number of PCs
ElbowPlot(integrated_data)

# finish pre-processing
integrated_data <- RunUMAP(integrated_data, dims = 1:20)
integrated_data <- FindNeighbors(object = integrated_data, dims = 1:20)              
integrated_data <- FindClusters(object = integrated_data, resolution = 0.1)

# identify doublets
# pK identification (no ground-truth)
sweep.list <- paramSweep(integrated_data, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.list)

# bcmvn = bi-modality coefficient
bcmvn <- find.pK(sweep.stats)

# optimal pK is the pK at which the bi-modal distribution of the distances 
# between cells is most pronounced
bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
optimal.pk <- bcmvn.max$pK
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]

# homotypic doublet proportion estimate 
# homotypic doublet: two cells of the same type trapped in one droplet
annotations <- integrated_data@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)

# adjust expected doublet rate based on total cell count in each sample
if (nrow(integrated_data@meta.data) < 5000) {
  doublet_rate <- 0.025  
} else if (nrow(integrated_data@meta.data) < 10000) {
  doublet_rate <- 0.05  
} else {
  doublet_rate <- 0.075  
}
nExp.poi <- round(doublet_rate * nrow(integrated_data@meta.data))
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

# run DoubletFinder
integrated_data <- doubletFinder(
  seu = integrated_data, PCs = 1:20, pK = optimal.pk, nExp = nExp.poi.adj)

View(integrated_data@meta.data)
# need to run individual samples through DoubletFinder
# create new metadata column named "Sample" with contents patient1_2hr, etc.