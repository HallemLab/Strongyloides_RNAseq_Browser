### Overview  
The *Strongyloides* RNA-seq Browser enables users to browse publicly
available *Strongyloides* bulk RNA-seq datasets and perform on-demand
analyses including differential expression and gene set enrichment. Data
from the following species are currently included: *S. stercoralis*, *S.ratti*,
*S. papillosus*, and *S. venezuelensis*. The app permits
browsing RNA-seq data in two modes:

1.  Browse by Gene Mode
2.  Browse by Life Stage Mode

### App Features  
Features of the app include:

-   Search for gene(s) of interest using stable geneIDs or keywords
-   Extract gene expression values for genes of interest
-   Display gene expresion across life stages as a heatmap (all
genes of interest) or a boxplot (individual genes)
- Display gene expression across life stages for individual genes and their known *Strongyloides* homologs
-   Download log2 counts per million expression for genes of
interest as .xslx
-   On demand limma-voom-based pairwise differential gene expression
analysis
-   Display results as interactive volcano plots and datatables
-   Download results as .pdf (plots) or .xlsx (datatables)
-   Gene set enrichment analysis using the clusterProfiler R package and
an Ensembl Compara protein family set established by Hunt *et al* 2016. 
(Browse By Life Stage Mode)
-   Display results as bubble plot and interactive datatable
-   Download results as .pdf (plots) or .xlsx (datatables)
-   Download raw/pre-processed data using user-friendly dropdown menu
-   Study desgin files (.csv)
-   Log2 counts per million expression for all genes and all samples
(.csv)
-   Variance-stabilized DGEList object (R object; primary data input
for the app)
-   Raw expression data for genes discarded during low-count
filterering; see Pre-processing files for more information
(.csv)
