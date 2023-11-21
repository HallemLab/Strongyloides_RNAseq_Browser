## GW: Assemble the Volcano Plot ----
pull_DEGs_GW <- reactive({
    req(vals$comparison_GW)
    req(vals$list.myTopHits.df_GW)
    req(vals$v.DEGList.filtered.norm)
    setProgress(0.5)
    if (isTruthy(input$displayedComparison_GW)){
        vals$displayedComparison_GW <- match(input$displayedComparison_GW,
                                             vals$comparison_GW, nomatch = 1)
    } else {vals$displayedComparison_GW <- 1}
    
    if(isTruthy(input$bkgndVolcanoYN_GW)){
        vals$includeBkgnd <- TRUE
    } else {vals$includeBkgnd <- FALSE}
    
    if(isTruthy(input$scaleVolcanoYN_GW)){
        vals$scaletogenome <- FALSE
    } else {vals$scaletogenome <- TRUE}
  
    #### Volcano Plots
    point_labels <- if(nrow(vals$list.highlight.tbl_GW[[vals$displayedComparison_GW]])<
                       20){guide_legend(override.aes = list(size = 4))} else {FALSE}
    setProgress(0.6)
    vplot <- ggplot(vals$list.myTopHits.df_GW[[vals$displayedComparison_GW]]) +
        aes(y=BH.adj.P.Val, x=logFC) +
        scale_y_continuous(trans = trans_reverser('log10')) + { 
            if (vals$includeBkgnd == TRUE) geom_point(size=2, color="grey",
                   na.rm = T) } +
        geom_point(data = vals$list.highlight.tbl_GW[[vals$displayedComparison_GW]], 
                   mapping = aes(y=(BH.adj.P.Val), 
                                 x=logFC, 
                                 color = geneID),
                   size = 3,
                   na.rm = T) +
        geom_hline(yintercept = -log10(adj.P.thresh), 
                   linetype="longdash", 
                   colour="black", 
                   size=1) + 
        geom_vline(xintercept = lfc.thresh, 
                   linetype="longdash", 
                   colour="#BE684D", 
                   size=1) +
        geom_vline(xintercept = -lfc.thresh, 
                   linetype="longdash", 
                   colour="#2C467A", 
                   size=1) + {
                       if (vals$scaletogenome == TRUE) expand_limits(x = c(min(vals$list.myTopHits.df_GW[[vals$displayedComparison_GW]]$logFC), max(vals$list.myTopHits.df_GW[[vals$displayedComparison_GW]]$logFC)), y = c(min(vals$list.myTopHits.df_GW[[vals$displayedComparison_GW]]$BH.adj.P.Val), max(vals$list.myTopHits.df_GW[[vals$displayedComparison_GW]]$BH.adj.P.Val))) } +
        guides(size = "none",
               colour = point_labels)+
        labs(title = paste0('Pairwise Comparison: ',
                            gsub('-',
                                 ' vs ',
                                 vals$comparison_GW[vals$displayedComparison_GW])),
             subtitle = paste0("black dashed line: p = ",
                               adj.P.thresh, "; colored lines: log-fold change = +/-",  lfc.thresh),
             color = "GeneIDs",
             y = "BH-adjusted p-value",
             x = "log2FC") +
        theme_Publication() +
        theme(aspect.ratio=1/3)
    setProgress(0.9)
    vplot
    
})


## GW: Volcano Plot, Generate UI  ----

    
output$volcano_GW <- renderUI({
    parse_contrasts_GW()
    req(vals$genelist,vals$comparison_GW)
    
    output$volcano_UI_GW <- renderPlot({
        withProgress({
        set_linear_model_GW()
        pull_DEGs_GW()}, message = "Calculating DGE...")
    })
    
    panel(
        heading = tagList(h5(shiny::icon("fas fa-mountain"),
                             "Pairwise Differential Gene Expression: Volcano Plot")),
        status = "primary",
        
        tagList(plotOutput('volcano_UI_GW',
                           hover = hoverOpts("plot_hover", 
                                             delay = 100, 
                                             delayType = "debounce")),
                uiOutput("hover_info"),
                uiOutput("bkgndVolcanoGW"),
                uiOutput('scaleVolcanoGW'),
                uiOutput("downloadVolcanoGW")
        )
    )
})

output$bkgndVolcanoGW <- renderUI({
    req(input$goLifeStage_GW, vals$comparison_GW,vals$genelist)
    checkboxInput("bkgndVolcanoYN_GW",
                  p("Show the DEG values for the entire genome as black dots in the plot background."), value = TRUE)
})

output$scaleVolcanoGW <- renderUI({
    req(input$goLifeStage_GW, vals$comparison_GW,vals$genelist)
    checkboxInput("scaleVolcanoYN_GW",
                  p("If values for entire genome are not shown, rescale plot axis to values of seached genes only?"), value = FALSE)
})

## GW: Save Volcano Plot ----
output$downloadVolcanoGW <- renderUI({
    req(input$goLifeStage_GW, vals$comparison_GW,vals$genelist)
    
    output$downloadVolcano_GW <- downloadHandler(
        filename = function(){
            paste('VolcanoPlot_',vals$comparison_GW[vals$displayedComparison_GW], '_',Sys.Date(),'.pdf', sep='')
        },
        content = function(file){
            withProgress({
                vplot<-pull_DEGs_GW()
                ggsave(file, 
                       plot = vplot, 
                       width = 11, 
                       height = 8, 
                       units = "in", 
                       device = cairo_pdf)
            },
            message = "Saving Plot")
        }
    )
    
    downloadButton("downloadVolcano_GW",
                   "Download Plot as PDF",
                   class = "btn-primary")
    
})

## GW: Volcano Hover Info ----
output$hover_info <- renderUI({
    req(vals$displayedComparison_GW,vals$genelist,vals$comparison_GW)
    
    pointer.df <- vals$list.highlight.tbl_GW[[vals$displayedComparison_GW]] %>%
        dplyr::mutate(log10.adj.P.Val = -log10(BH.adj.P.Val))
    
    hover <- input$plot_hover
    point <- nearPoints(pointer.df, hover, 
                        xvar = "logFC",
                        yvar = "log10.adj.P.Val",
                        threshold = 5, maxpoints = 1, addDist = TRUE)
    if (nrow(point) == 0) return(NULL)
    
    # calculate point position INSIDE the image as percent of total dimensions
    # from left (horizontal) and from top (vertical)
    left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
    top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
    
    # calculate distance from left and bottom side of the picture in pixels
    left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
    top_px <- hover$coords_img$y - top_pct * (hover$range$bottom - hover$range$top)
    
    left_px <- (hover$coords_img$x + left_pct)/ hover$img_css_ratio$x 
    #top_px <- (hover$coords_img$y - top_pct) / hover$imge_css_ratio$y
    
    #left_px <- hover$coords_img$x + hover$x
    #top_px <- hover$coords_img$y + hover$y
    
    # create style property fot tooltip
    # background color is set so tooltip is a bit transparent
    # z-index is set so we are sure are tooltip will be on top
    style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                    #"right:", hover$x + 10, "px; top:", hover$y + 10, "px;")
                    #"left:", hover$range$left, "px; bottom:", hover$range$bottom , "px;")
                    "left:", left_px + 50, "px; bottom:", top_px + 5, "px;")
    
    # actual tooltip created as wellPanel
    wellPanel(
        style = style,
        p(HTML(paste0("<b> GeneID: </b>", 
                      point$geneID, 
                      "<br/>", 
                      "<b> Log FC: </b>",
                      round(point$logFC,digits = 2),
                      "<br/>",
                      "<b> p-value: </b>",
                      format(point$BH.adj.P.Val, digits = 3, scientific = TRUE))))
    )
})


## GW: Data Table of Differentially Expressed Genes from User Subset ----
assemble_DEGs_GW <- reactive({
    
    req(vals$comparison_GW,vals$displayedComparison_GW,vals$genelist)
    
    tS<- vals$targetStage_GW[vals$displayedComparison_GW,
                             ][vals$targetStage_GW[vals$displayedComparison_GW,
                                                   ]!=""]
    cS<- vals$contrastStage_GW[vals$displayedComparison_GW,
                               ][vals$contrastStage_GW[vals$displayedComparison_GW,
                                                       ]!=""]
    
    
    # Add back on genes that were submitted by the user but don't appear in the list of genes for which there is available data.
    excluded.genes <- dplyr::anti_join(vals$submitted.genelist, 
                                       vals$genelist,
                                       by = "geneID") %>%
        left_join(vals$annotations, by = "geneID") # Add gene annotations
    
    sample.num.tS <- sapply(tS, function(x) {colSums(vals$v.DEGList.filtered.norm$design)[[x]]}) %>% sum()
    sample.num.cS <- sapply(cS, function(x) {colSums(vals$v.DEGList.filtered.norm$design)[[x]]}) %>% sum()
   
    n_num_cols <- sample.num.tS + sample.num.cS + 5
    index_homologs <- length(colnames(vals$list.highlight.tbl_GW[[vals$displayedComparison_GW]])) - 6
    
    highlight.datatable <- vals$list.highlight.tbl_GW[[vals$displayedComparison_GW]] %>%
        {suppressMessages(dplyr::full_join(.,excluded.genes))} %>%
        dplyr::mutate(In.subclade_geneID = paste0("<a href='https://parasite.wormbase.org/Multi/Search/Results?species=all;idx=;q=", In.subclade_geneID,"' target = '_blank'>", In.subclade_geneID,"</a>"))%>%
        dplyr::mutate(Out.subclade_geneID = paste0("<a href='https://parasite.wormbase.org/Multi/Search/Results?species=all;idx=;q=", Out.subclade_geneID,"' target = '_blank'>", Out.subclade_geneID,"</a>"))%>%
        dplyr::mutate(Out2.subclade_geneID = paste0("<a href='https://parasite.wormbase.org/Multi/Search/Results?species=all;idx=;q=", Out2.subclade_geneID,"' target = '_blank'>", Out2.subclade_geneID,"</a>"))%>%
        dplyr::mutate(Ce_geneID = paste0("<a href='https://parasite.wormbase.org/Caenorhabditis_elegans_prjna13758/Gene/Summary?g=", Ce_geneID,"' target = '_blank'>", Ce_geneID,"</a>"))%>%
        dplyr::mutate(WBPSLink = paste0("<a href='https://parasite.wormbase.org/Multi/Search/Results?species=all;idx=;q=", geneID,"' target = '_blank'>", geneID,"</a>")) %>%
        dplyr::relocate(UniProtKB, Description, InterPro, GO_term,
                        In.subclade_geneID, In.subclade_percent_homology,
                        Out.subclade_geneID, Out.subclade_percent_homology,
                        Out2.subclade_geneID, Out2.subclade_percent_homology,
                        Ce_geneID, Ce_percent_homology, .after = last_col())  %>%
        DT::datatable(rownames = FALSE,
                      escape = FALSE,
                      caption = htmltools::tags$caption(
                          style = 'caption-side: top; text-align: left; color: black',
                          htmltools::tags$b('Differentially Expressed Genes in', 
                                            htmltools::tags$em(input$selectSpecies_GW), 
                                            gsub('-',' vs ',vals$comparison_GW[vals$displayedComparison_GW])),
                          htmltools::tags$br(),
                          "Threshold: p < ",
                          adj.P.thresh, "; log-fold change > ",
                          lfc.thresh,
                          htmltools::tags$br(),
                          'Values = log2 counts per million'),
                      options = list(autoWidth = TRUE,
                                     scrollX = TRUE,
                                     scrollY = '300px',
                                     scrollCollapse = TRUE,
                                     order = list(n_num_cols-1, 
                                                  'desc'),
                                     searchHighlight = TRUE, 
                                     pageLength = 10, 
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
                                         list(
                                             targets = ((n_num_cols + 
                                                             5):(n_num_cols + 
                                                                     6)),
                                             render = JS(
                                                 "function(data, type, row, meta) {",
                                                 "return type === 'display' && data.length > 20 ?",
                                                 "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                                                 "}")
                                         ),
                                         list(targets = "_all",
                                              class="dt-right")
                                     ),
                                     rowCallback = JS(c(
                                         "function(row, data){",
                                         "  for(var i=0; i<data.length; i++){",
                                         "    if(data[i] === null){",
                                         "      $('td:eq('+i+')', row).html('NA')",
                                         "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
                                         "    }",
                                         "  }",
                                         "}"  
                                     ))
                                     
                      )) 

    highlight.datatable <- highlight.datatable %>%
        DT::formatRound(columns=c(3:(n_num_cols-3),
                                  n_num_cols), 
                        digits=3)
    
    highlight.datatable <- highlight.datatable %>%
        DT::formatRound(columns=c(n_num_cols+2, 
                                  index_homologs+1,
                                  index_homologs+3,
                                  index_homologs+5), 
                        digits=2)
    
    highlight.datatable <- highlight.datatable %>%
        DT::formatSignif(columns=c(n_num_cols+1), 
                         digits=3)
    highlight.datatable
})

output$highlight.df <- renderDT ({
    req(input$goGW,vals$list.highlight.tbl_GW,vals$genelist)
    DEG.datatable_GW<-assemble_DEGs_GW()
    DEG.datatable_GW
})

##GW: Generate responsive pulldown list for filtering/saving DEG Tables
output$downloadSelectionMenu_GW <- renderUI({
    req(input$goGW,vals$list.highlight.tbl_GW,vals$comparison_GW)
    req(str_detect(names(vals$list.highlight.tbl_GW),paste0(gsub("\\+","\\\\+",vals$comparison_GW) %>%
                                                                gsub("\\-","\\\\-",.)  ,
                                                            collapse = "|")))
    
    selectInput("download_DGE_Selection_GW",
                label = h6("Select Contrast(s) to Download"),
                choices = c('Choose one or more' = ''
                            ,'Everything' = 'Everything'
                            ,as.list(vals$comparison_GW)),
                selected = "Everything",
                selectize = TRUE,
                multiple = TRUE)
})

##GW: Filter DEG Tables Before Saving ----
filter_DEG_tbl_GW <- reactive({
    req(input$goLifeStage_GW,vals$displayedComparison_GW,vals$list.highlight.tbl_GW,vals$comparison_GW,input$download_DGE_Selection_GW)
    req(str_detect(names(vals$list.highlight.tbl_GW),paste0(gsub("\\+","\\\\+",vals$comparison_GW) %>%
                                                                gsub("\\-","\\\\-",.)  ,
                                                            collapse = "|")))
    ## Figure out which datasets to download
    if (any(str_detect(input$download_DGE_Selection_GW, "Everything"))){
        download_DT <- vals$list.highlight.tbl_GW
    } else {
        subsetContrasts <- str_detect(names(vals$list.highlight.tbl_GW),paste0(input$download_DGE_Selection_GW, collapse = "|"))
        download_DT<- vals$list.highlight.tbl_GW[subsetContrasts]
    }
    
    if (input$download_DGEdt_across_GW == TRUE) {
        subsetIDs <-lapply(names(download_DT), function(y){
            dplyr::filter(download_DT[[y]],str_detect(DEG_Desc,paste0(input$download_DGEdt_direction_GW, collapse = "|")))
        }) %>%
            purrr::reduce(inner_join, by = c("geneID", "DEG_Desc")) %>%
            dplyr::select(geneID)
        
        filtered.list.highlight.tbl_GW <- lapply(names(download_DT), function(y){
            dplyr::filter(download_DT[[y]], geneID %in% subsetIDs$geneID)%>%
                dplyr::arrange(desc(logFC))
        })
        names(filtered.list.highlight.tbl_GW) <- names(download_DT)
    } else {
        filtered.list.highlight.tbl_GW <-lapply(names(download_DT), function(y){
            dplyr::filter(download_DT[[y]],str_detect(DEG_Desc,paste0(input$download_DGEdt_direction_GW, collapse = "|"))) %>%
                dplyr::arrange(desc(logFC))
        })
        names(filtered.list.highlight.tbl_GW)<- names(download_DT)
    }
    
    # Pass only a specific proportion of genes. Remember, the datatable is grouped by the DEG_Desc value, and ordered by descending logFC value. 
    filtered.list.highlight.tbl_GW<-lapply(names(filtered.list.highlight.tbl_GW), function(y){
        
        filtered.list.highlight.tbl_GW[[y]] %>%
            dplyr::group_map(~ {
                if (str_detect(.y, "Up")) {
                    slice_max(.x, order_by = logFC, prop = as.numeric(input$percentDGE_GW)/100)
                } else if (str_detect(.y, "NotSig")) {
                    slice_max(.x, order_by = logFC, prop = as.numeric(input$percentDGE_GW)/100)
                } else if (str_detect(.y, "Down")) {
                    slice_min(.x, order_by = logFC, prop = as.numeric(input$percentDGE_GW)/100)
                }
            }, .keep = TRUE) %>%
            bind_rows() %>%
            dplyr::arrange(desc(logFC))
    })
    names(filtered.list.highlight.tbl_GW)<- names(download_DT)
    filtered.list.highlight.tbl_GW
})

## GW: Save Excel Tables with DEG Tables ----
output$downloadbuttonGW <- renderUI({
    req(input$goLifeStage_GW,vals$comparison_GW,vals$genelist)
    
    filtered.tbl_GW <- filter_DEG_tbl_GW()
    
    if (input$download_missing_genes_GW == T) {
    # Add back on genes that were submitted by the user but don't appear in the list of genes for which there is available data.
    excluded.genes <- dplyr::anti_join(vals$submitted.genelist, 
                                       vals$genelist,
                                       by = "geneID") %>%
        left_join(vals$annotations, by = "geneID") # Add gene annotations
    
    downloadablel.tbl_GW <-lapply(filtered.tbl_GW, function (x) {
        suppressMessages(dplyr::full_join(x,excluded.genes))
    })
    } else {downloadablel.tbl_GW <- filtered.tbl_GW}
   
    expressionnotes <- "Columns labeled with <life stage - sample ID> report log2 counts per million (CPM) expression. Columns labeled avg_<life stage> are median log2CPM with IQR. The column labeled logFC reports log2 fold change."
    ### Generate some text that describes the analysis conditions
    ### 1. If p-values are adjusted for multiple pairwise comparisons, which comparisons are included in the adjustment parameters? This should be the list of selected contrasts
    if (vals$multipleCorrection_GW == TRUE){
        multiplecorrection <- paste0(vals$comparison_GW, collapse = "; ")
        multiplecorrection <- paste0("P-values corrected across the following multiple pairwise comparisons: ", multiplecorrection)
    } else {
        multiplecorrection <- "P-values *not* corrected across multiple pairwise comparisons"
    }
    ### 2. If downloading results are being filtered to show only genes that display consistent differential expression across all comparisons targeted for download, list the pairwise comparisons being used.
    if (input$download_DGEdt_across_GW == TRUE){
        filteredacross <- paste0(names(filtered.tbl_GW), collapse = "; ")
        filteredacross <- paste0("Lists only include genes with matching differential expression descriptions across the following pairwise comparisons: ", filteredacross)
    } else filteredacross <- ""
        
    ### 3. Specify which types of differential expression pattern
    DEGpattern <- paste0(input$download_DGEdt_direction_GW, collapse = "; ")
    DEGpattern <- paste0("Lists include genes that are differentially regulated in the following directions (relative to the target life stage): ", DEGpattern)
    
    ### 4. Specify the proportion of genes for each DEG Type that are being saved (e.g. top 10% of upregulated, and top10% of downregulated genes)
    proportionexport <- paste0("Percentage of genes for each differential expression pattern: ", input$percentDGE_GW)
    
    output$generate_excel_report_GW <- generate_excel_report(names(downloadablel.tbl_GW), 
                                                             downloadablel.tbl_GW,
                                                             name = paste(input$selectSpecies_GW,
                                                                          "RNA-seq Differential Gene Expression"),
                                                             filename_prefix = 'DGETable_',
                                                             expressionnotes = expressionnotes,
                                                             multiplecorrection = multiplecorrection,
                                                             filteredacross = filteredacross,
                                                             DEGpattern = DEGpattern,
                                                             proportionexport = proportionexport)
    
    downloadButton("generate_excel_report_GW",
                   "Download DGE Tables as Excel",
                   class = "btn-primary")
})
