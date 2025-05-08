library(Seurat)
library(pryr)

# 1. Load your full object once
obj <- readRDS("data/wt_old_early.rds")
print(object_size(obj))  # ~1.1 GB

# 2. Subset to just Salivary Gland cells
sg <- subset(obj, manual_celltypes == "Salivary Gland")

# 3. Diet it down to only the RNA assay and UMAP (and drop any counts/data you don't need)
sg_slim <- DietSeurat(
  object   = sg,
  assays    = "RNA",         # drop everything but the RNA assay
  dimreducs = "umap",        # drop PCA, neighbors, etc. if you don’t need them
  graphs    = NULL           # drop neighbor/graph slots
)

# now nuking the raw slots you’ll never call in Shiny:
DefaultAssay(sg_slim) <- "RNA"
#sg_slim[["RNA"]]@counts    <- NULL   # raw counts
#sg_slim[["RNA"]]@data      <- NULL   # log‐normalized data *if* you only need FetchData; otherwise keep this
# sg_slim[["RNA"]]@scale.data <- NULL # scale.data is even bigger—drop if you never use it

print(object_size(sg_slim))  # hopefully < 100 MB now

# 4. Save it out
saveRDS(sg_slim, "data/wt_late_salivary_slim.rds", compress = TRUE)
