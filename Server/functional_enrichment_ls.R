## LS: Run GSEA ----
perform_GSEA_LS <-reactive ({
    # Shiny version of functional enrichment analysis. 
    # Performs GSEA using clusterProfiler
    # Ability to do this depends on the availability of gene sets. 
    # 
    # In Hunt et al 2016, there is an Ensembl Compara protein family set
    # Note that this uses specific transcript information, which I throw out. 
    # (e.g. SSTP_0001137400.2 is recoded as SSTP_0001137400)
    # I preprocessed this list offline, removing genes that are not found in
    # the preprocessed list of genes for which we have RNASeq data.
    req(vals$displayedComparison_LS)
    setProgress(0.1)
    # Generate rank ordered list of genes
    mydata.df.sub <- dplyr::select(vals$list.highlight.tbl_LS[[vals$displayedComparison_LS]], geneID, logFC)
    mydata.gsea <- mydata.df.sub$logFC
    names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
    mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)
    
    # run GSEA using the 'GSEA' function from clusterProfiler
    # Given a priori defined set of gene S (e.g., genes sharing the same DO category), the goal of GSEA is to determine whether the members of S are randomly distributed throughout the ranked gene list (L) or primarily found at the top or bottom.
    # There are three key elements of the GSEA method:
    # **Calculation of an Enrichment Score.**
    # The enrichment score (ES) represent the degree to which a set S is over-represented at the top or bottom of the ranked list L. The score is calculated by walking down the list L, increasing a running-sum statistic when we encounter a gene in S and decreasing when it is not. The magnitude of the increment depends on the gene statistics (e.g., correlation of the gene with phenotype). The ES is the maximum deviation from zero encountered in the random walk; it corresponds to a weighted Kolmogorov-Smirnov-like statistic (Subramanian et al. 2005).
    # **Esimation of Significance Level of ES.**
    # The p-value of the ES is calculated using permutation test. Specifically, we permute the gene labels of the gene list L and recompute the ES of the gene set for the permutated data, which generate a null distribution for the ES. The p-value of the observed ES is then calculated relative to this null distribution.
    # **Adjustment for Multiple Hypothesis Testing.**
    # When the entire gene sets were evaluated, DOSE adjust the estimated significance level to account for multiple hypothesis testing and also q-values were calculated for FDR control.
    setProgress(0.2)
    myGSEA.res <- suppressWarnings(GSEA(mydata.gsea, TERM2GENE=ensComp, verbose=FALSE))
    myGSEA.df <- as_tibble(myGSEA.res@result)
    setProgress(0.6)
    
    myGSEA.df <- myGSEA.df %>%
        dplyr::mutate(life_stage = case_when(
            NES > 0 ~ str_split(vals$comparison_LS[vals$displayedComparison_LS],'-',simplify = T)[1,1],
            NES < 0 ~ str_split(vals$comparison_LS[vals$displayedComparison_LS],'-',simplify = T)[1,2]))
    
    myGSEA.df$ID <- myGSEA.df$ID %>%
        word(sep = ',') %>%
        word(sep = ' and')
    
    vals$myGSEA.df <- myGSEA.df
    setProgress(0.8)
    ggplot(myGSEA.df, aes(x=life_stage, y=ID)) + 
        geom_point(aes(size=setSize, color = NES, alpha=p.adjust)) +
        scale_color_gradient(low="blue", high="red") +
        scale_alpha(name = "Adjusted p-value",
                    trans = trans_reverser('log10')) +
        guides(alpha = guide_legend(override.aes = list(size=3,
                                                        color = "red"))) +
        scale_size(name = "Set Size",
                   breaks = c(0,10,50,100,200),
                   labels = c("0", "10", "50", "100", "200"))+
        labs(title = paste0('Gene Families Enriched in ', 
                            gsub('-',' vs ',
                                 vals$comparison_LS[vals$displayedComparison_LS])),
             subtitle = 'NES = Normalized Enrichment Score; Gene family assignments 
             from Ensembl Compara dataset defined in Hunt et al 2016',
             x = "Life Stage",
             y = "Family ID") +
        theme_bw() +
        theme(plot.title.position = "plot",
              plot.caption.position = "plot",
              plot.title = element_text(face = "bold",
                                        size = 13, hjust = 0),
              axis.title = element_text(face = "bold",size = 10.4),
              legend.title = element_text(face="bold",size = 10.4),
              aspect.ratio = 3/1)
    
})

## LS: GSEA Plot ----
output$GSEAPlot_LS <- renderPlot({
    withProgress({perform_GSEA_LS()}, message = "Running GSEA...")
})

## LS: Save GSEA Plot ----
output$downloadGSEAPlot_LS <- downloadHandler(
    filename = function(){
        paste('GSEAPlot_',vals$comparison_LS[vals$displayedComparison_LS], '_',Sys.Date(),'.pdf', sep='')
    },
    content = function(file){
        gseaplot<- perform_GSEA_LS()
        ggsave(file,
               plot = gseaplot, 
               width = 11, 
               height = 8, 
               units = "in", 
               device = cairo_pdf)
    }
)

## LS: GSEA Data Table ----
output$GSEATbl_LS <- renderDT({
    req(vals$myGSEA.df)
    
    myGSEA.tbl<-vals$myGSEA.df %>%
        dplyr::select(!c(Description, pvalue, enrichmentScore,life_stage))
    
    # view results as an interactive table
    enrichment.DT <- datatable(myGSEA.tbl, 
                               rownames = TRUE,
                               caption =  htmltools::tags$caption(
                                   style = 'caption-side: top; text-align: left; color: black',
                                   htmltools::tags$b('Gene Families Enriched in ', 
                                                     gsub('-',' vs ',
                                                          vals$comparison_LS[vals$displayedComparison_LS])),
                                   htmltools::tags$br(),
                                   'NES = Normalized Enrichment Score',
                                   htmltools::tags$br(),
                                   'Gene family assignments from Ensembl Compara dataset defined in Hunt',
                                   htmltools::tags$em('et al'), '2016.'),
                               options = list(
                                   autoWidth = TRUE,
                                   scrollX = TRUE,
                                   scrollY = '300px',
                                   scrollCollapse = TRUE,
                                   searchHighlight = TRUE, 
                                   order = list(3, 'desc'),
                                   pageLength = 25, 
                                   lengthMenu = c("5",
                                                  "10",
                                                  "25",
                                                  "50",
                                                  "100"),
                                   initComplete = htmlwidgets::JS(
                                       "function(settings, json) {",
                                       paste0("$(this.api().table().container()).css({'font-size': '", "10pt", "'});"),
                                       "}"),
                                   columnDefs = list(
                                       list(targets = "_all",
                                            class="dt-right"),
                                       list(
                                           targets = c(1,8),
                                           render = JS(
                                               "function(data, type, row, meta) {",
                                               "return type === 'display' && data.length > 40 ?",
                                               "'<span title=\"' + data + '\">' + data.substr(0, 40) + '...</span>' : data;",
                                               "}")
                                       )))) %>%
        formatRound(columns=c(3,5:6), digits=2) %>%
        formatRound(columns=c(4), digits=4)
    enrichment.DT
})

## LS: Save GSEA Datatable ----
output$downloadGSEAtbl_LS <- renderUI({
    req(input$goLS)
    req(vals$myGSEA.df)
    myGSEA.tbl<-vals$myGSEA.df %>%
        dplyr::select(!c(Description, pvalue, enrichmentScore,life_stage)) %>%
        dplyr::arrange(desc(NES)) %>%
        list("GSEA" = . )
    
    output$generate_GSEA_report_LS <- generate_excel_report(vals$comparison_LS[vals$displayedComparison_LS], 
                                                            myGSEA.tbl,
                                                            name = "Gene Set Enrichment Analysis",
                                                            filename_prefix = paste('GSEATable_',
                                                                                    vals$comparison_LS[vals$displayedComparison_LS],
                                                                                    '_', 
                                                                                    sep=''),
                                                            
                                                            n_header_rows = 5)
    downloadButton("generate_GSEA_report_LS",
                   "Download GSEA Table",
                   class = "btn-primary")
    
})