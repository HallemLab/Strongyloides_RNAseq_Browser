# Strongyloides RNAseq Browser Shiny App

# ---- Libraries ----
suppressPackageStartupMessages({
    library(shiny)
    library(shinyjs)
    library(shinyWidgets)
    library(htmltools)
    library(shinythemes)
    library(DT)
    library(Cairo)
    library(tidyverse)
    library(ggforce)
    library(limma)
    library(edgeR)
    library(gt) 
    library(plotly)
    library(magrittr)
    library(ggthemes)
    library(gplots)
    library(svglite)
    library(Cairo)
    library(heatmaply)
    library(RColorBrewer)
    library(rcartocolor)
    library(openxlsx)
    library(egg)
    library(dendextend)
    library(vctrs)
    library(clusterProfiler)
    library(markdown)
    source("Server/ggheatmap_local.R")
    source("Server/calc_DEG_tbl.R")
    source("Server/theme_Publication.R")
    source('Server/limma_ranking.R')
    source("Server/generate_excel_report.R")
})

# ---- Background ----
source('Server/preprocess_data.R', local = T)

# ---- UI ----

ui <- fluidPage(
    
    tags$head(
        HTML('<base target="_blank">')
    ),
    
    source('UI/Str_RNAseqBrowser_ui.R', local = T)$value,
    
    source('UI/custom_css.R', local = T)$value
    
)


# ---- Server  ----
server <- function(input, output, session) {
    
    vals<-reactiveValues()
    
    ## Initialize Species ----
    ## GW: Load species data in Gene-wise tab
    ## LS: Load species data in Life-stage tab
    ## GW: Generate UI for downloading study information for loaded species
    ## LS: Generate UI for downloading study information for loaded species
    ## About Tab: Generate UI for downloading study information for all species
    source('Server/load_Species.R', local = T)$value
    
    ## GW: Gene Expression Across Life Stage ----
    
    ## GW: Generate/Reset Gene File Upload
    ## GW: Clear Genes
    ## GW/LS: Reset Elements on Tab Change 
    source('Server/reset_state.R', local = T)
    
    ## GW: Generate Main Gene Input Panel
    ## GW: Parse Gene Inputs
    source('Server/parse_geneIDs_gw.R', local = T)
    
    ## GW: Generate Panel for Gene-wise Plots
    ## GW: Generate Responsive Selection for Gene to Display
    ## GW: Generate Legend Explaining the Life Stages
    source('Server/generate_geneplot_UI_gw.R', local = T)
    
    
    ## GW: Plots of Gene Expression
    ## GW: Switch what type of GW Plot is produced
    ## GW: Save Gene Plots
    ## GW: Save Excel Table with Gene Expression Data
    source('Server/plot_save_gene_expression_gw.R', local = T)
    
    
    ## GW: Pairwise Comparisons Across Life Stage ----

    ## GW: Clear Comparison Selections
    ## GW: Generate Comparison Selection Boxes
    ## GW: Parse the inputs
    ## GW: Generate Responsive Selection for Life Stage to Display
    source('Server/initialize_pairwiseComparisons_gw.R', local = T)
    
    
    ## GW: Set Contrast Matrix, Fit the Linear Model, Extract the Differentially Expressed Genes
    source('Server/set_linear_model_gw.R', local = T)
    
    
    ## GW: Assemble the Volcano Plot 
    ## GW: Volcano Plot, Generate UI
    ## GW: Save Volcano Plot
    ## GW: Volcano Hover Info
    ## GW: Data Table of Differentially Expressed Genes from User Subset
    ## GW: Save Excel Tables with DEG Tables
    source('Server/plot_save_DEGs_gw.R', local = T)
    
    ## LS: Pairwise Comparisons Across Life Stage ----
    ## LS: Clear Comparison Selections
    ## LS: Pairwise comparisons Across Life Stages
    ## LS: Generate Responsive Selection for Life Stage to Display
    ## LS: Generate Legend Explaining the Life Stages ----
    source('Server/initialize_pairwiseComparisons_ls.R', local = T)
    
    ## LS: Set Contrast Matrix, Fit the Linear Model, Extract the Differentially Expressed Genes
    source('Server/set_linear_model_ls.R', local = T)
    
    ## LS: Assemble the Volcano Plot 
    ## LS: Volcano Plot, Generate UI
    ## LS: Save Volcano Plot
    ## LS: Volcano Hover Info
    ## LS: Data Table of Differentially Expressed Genes
    ## LS: Save Excel Tables with DEG Tables
    source('Server/plot_save_DEGs_ls.R', local = T)
    
    ## LS: Functional Enrichment Analysis ----
    ## LS: Run GSEA
    ## LS: GSEA Plot
    ## LS: Save GSEA Plot
    ## LS: GSEA Data Table
    ## LS: Save GSEA Datatable
    if (!config::is_active("shinyapps")){
        source('Server/functional_enrichment_ls.R', local = T)
    } else {
        
    }
    ## Output Options ----
    outputOptions(output, "downloadbuttonsGenes", suspendWhenHidden = FALSE)
    outputOptions(output, "volcano_GW", suspendWhenHidden = FALSE)
    outputOptions(output, "pairwiseSelector_GW", suspendWhenHidden = FALSE)
    outputOptions(output, "pairwiseSelector_LS", suspendWhenHidden = FALSE)
    
    session$onSessionEnded(stopApp)
}

# Run the application 
shinyApp(ui = ui, server = server)
