# Coexpression Visualizer

https://fly-singlecell.shinyapps.io/corr_app/

This is a Shiny app to interactively explore gene co-expression patterns in single-cell RNA-seq data.  
Built using Seurat and ggplot2, the app lets users select a dataset, choose any two genes, and visualize:

- Overlayed UMAP expression plots with separate color scales for each gene
- Gene–gene scatter plots with correlation statistics
- Zoomable plots for focused exploration

## Data
The dataset used here is from Peng et al. 2024 Organogenetic transcriptomes of the Drosophila embryo at single cell resolution.\
The data is subsetted to focus on just the salivary gland due to size limitations in hosting directly on Shiny.\
You can find the original scRNAseq data here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234602\
And the code that was used to analyse and process the data is hosted here: https://github.com/CahanLab/fruitfly_organogenesis


## Sample Run 
<img width="1411" alt="Screenshot 2025-05-08 at 12 16 52 PM" src="https://github.com/user-attachments/assets/c1a9b0bc-4ec2-4a84-9988-c9753cc88a4e" />

<img width="971" alt="Screenshot 2025-05-08 at 12 17 06 PM" src="https://github.com/user-attachments/assets/b6bc6efe-633e-4c12-baf8-df0409697ab7" />

## Using other datasets
To use this app to visualize other (heavier) data, you can clone the repo and run locally. Make sure that your seurat objects are in the `data` folder and the cells are annotated.
```bash
git clone https://github.com/vidyaajay1/coexp_visualizer.git
cd coexp_visualizer
Rscript -e "shiny::runApp()"
```
