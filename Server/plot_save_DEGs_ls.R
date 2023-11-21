## LS: Assemble the Volcano Plot ----
pull_DEGs_LS <- reactive({
  req(vals$comparison_LS)
  req(vals$list.highlight.tbl_LS)
  
  setProgress(0.5)
  if (isTruthy(input$displayedComparison_LS)){
    vals$displayedComparison_LS <- match(input$displayedComparison_LS,
                                         vals$comparison_LS, nomatch = 1)
  } else {vals$displayedComparison_LS <- 1}
  setProgress(0.6)
  #### Volcano Plots
  vplot <- ggplot(vals$list.highlight.tbl_LS[[vals$displayedComparison_LS]]) +
    aes(y=BH.adj.P.Val, x=logFC) +
    scale_y_continuous(trans = trans_reverser('log10')) +
    geom_point(size=2, color = "grey") +
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
               size=1) +
    labs(title = paste0('Pairwise Comparison: ',
                        gsub('-',
                             ' vs ',
                             vals$comparison_LS[vals$displayedComparison_LS])),
         subtitle = paste0("black dashed line: p = ",
                           adj.P.thresh, "; colored lines: log-fold change = +/-", lfc.thresh),
         color = "GeneIDs",
         y = "BH-adjusted p-value",
         x = "log2FC") +
    theme_Publication() +
    theme(aspect.ratio=1/3)
  setProgress(0.9)
  vplot
})

## LS: Volcano Plot, Generate UI ----
output$volcano_LS <- renderUI({
  req(input$goLS)
  parse_contrasts_LS()
  req(vals$comparison_LS)
  
  output$volcano_UI_LS <- renderPlot({
    withProgress({
      set_linear_model_LS()
      pull_DEGs_LS()}, 
      message = "Calculating DGE...")
  })
  
  panel(
    heading = tagList(h5(shiny::icon("fas fa-mountain"),
                         "Pairwise Differential Gene Expression: Volcano Plot")),
    status = "primary",
    plotOutput('volcano_UI_LS',
               hover = hoverOpts("plot_hover_LS",
                                 delay = 100,
                                 delayType = "debounce")),
    
    uiOutput("hover_info_LS"),
    
    uiOutput("downloadVolcanoLS")
  )
  
})

## LS: Save Volcano Plot ----
output$downloadVolcanoLS <- renderUI({
  req(input$goLS,vals$comparison_LS)
  
  output$downloadVolcano_LS <- downloadHandler(
    filename = function(){
      paste('VolcanoPlot_',vals$comparison_LS[vals$displayedComparison_LS], '_',Sys.Date(),'.pdf', sep='')
    },
    content = function(file){
      withProgress({
        setProgress(.25)
        vplot<-pull_DEGs_LS()
        setProgress(.75)
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
  
  downloadButton("downloadVolcano_LS",
                 "Download Plot as PDF",
                 class = "btn-primary")
  
})


## LS: Volcano Hover Info ----
output$hover_info_LS <- renderUI({
  req(vals$comparison_LS,vals$displayedComparison_LS)
  
  pointer.df <- vals$list.highlight.tbl_LS[[vals$displayedComparison_LS]] %>%
    dplyr::mutate(log10.adj.P.Val = -log10(BH.adj.P.Val))
  
  hover <- input$plot_hover_LS
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

## LS: Data Table of Differentially Expressed Genes ----
assemble_DEGs_LS <- reactive({
  req(vals$comparison_LS,vals$displayedComparison_LS)
  
  tS<- vals$targetStage_LS[vals$displayedComparison_LS,
                           ][vals$targetStage_LS[vals$displayedComparison_LS,
                                                 ]!=""]
  cS<- vals$contrastStage_LS[vals$displayedComparison_LS,
                             ][vals$contrastStage_LS[vals$displayedComparison_LS,
                                                     ]!=""]
  sample.num.tS <- sapply(tS, function(x) {colSums(vals$v.DEGList.filtered.norm$design)[[x]]}) %>% sum()
  sample.num.cS <- sapply(cS, function(x) {colSums(vals$v.DEGList.filtered.norm$design)[[x]]}) %>% sum()
  
  n_num_cols <- sample.num.tS + sample.num.cS + 5
  index_homologs <- length(colnames(vals$list.highlight.tbl_LS[[vals$displayedComparison_LS]])) - 6

  LS.datatable <- vals$list.highlight.tbl_LS[[vals$displayedComparison_LS]] %>%
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
                                      htmltools::tags$em(input$selectSpecies_LS), 
                                      gsub('-',' vs ',vals$comparison_LS[vals$displayedComparison_LS])),
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
                                                   1)),
                                     render = JS(
                                       "function(data, row) {",
                                       "data.toExponential(1);",
                                       "}")
                                   ),
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
                                 )
                                 
                  )) 
  LS.datatable <- LS.datatable %>%
    DT::formatRound(columns=c(3:(n_num_cols-3),
                                 n_num_cols), 
                    digits=3)
  
  LS.datatable <- LS.datatable %>%
    DT::formatRound(columns=c(n_num_cols+2, 
                              index_homologs+1,
                              index_homologs+3,
                              index_homologs+5), 
                    digits=2)
  
  LS.datatable <- LS.datatable %>%
    DT::formatSignif(columns=c(n_num_cols+1), 
                     digits=3)
  
  LS.datatable
  
})

output$tbl_LS <- renderDT ({
  req(input$goLS,vals$list.highlight.tbl_LS)
  
  DEG.datatable_LS<-assemble_DEGs_LS()
  DEG.datatable_LS
})

##LS: Generate responsive pulldown list for filtering/saving DEG Tables
output$downloadSelectionMenu_LS <- renderUI({
  req(input$goLS,vals$list.highlight.tbl_LS,vals$comparison_LS)
  req(str_detect(names(vals$list.highlight.tbl_LS),paste0(gsub("\\+","\\\\+",vals$comparison_LS) %>%
                                                            gsub("\\-","\\\\-",.)  ,
                                                          collapse = "|")))
  
  selectInput("download_DGE_Selection_LS",
              label = h6("Select Contrast(s) to Download"),
              choices = c('Choose one or more' = ''
                          ,'Everything' = 'Everything'
                          ,as.list(vals$comparison_LS)),
              selected = "Everything",
              selectize = TRUE,
              multiple = TRUE)
})

##LS: Filter DEG Tables Before Saving ----
filter_DEG_tbl_LS <- reactive({
  req(input$goLS,vals$displayedComparison_LS,vals$list.highlight.tbl_LS,vals$comparison_LS,input$download_DGE_Selection_LS)
  req(str_detect(names(vals$list.highlight.tbl_LS),paste0(gsub("\\+","\\\\+",vals$comparison_LS) %>%
                                                            gsub("\\-","\\\\-",.)  ,
                                                          collapse = "|")))
  ## Figure out which datasets to download
  if (any(str_detect(input$download_DGE_Selection_LS, "Everything"))){
    download_DT <- vals$list.highlight.tbl_LS
  } else {
    subsetContrasts <- str_detect(names(vals$list.highlight.tbl_LS),paste0(input$download_DGE_Selection_LS, collapse = "|"))
    download_DT<- vals$list.highlight.tbl_LS[subsetContrasts]
  }
  
  if (input$download_DGEdt_across_LS == TRUE) {
    subsetIDs <-lapply(names(download_DT), function(y){
      dplyr::filter(download_DT[[y]],str_detect(DEG_Desc,paste0(input$download_DGEdt_direction_LS, collapse = "|")))
    }) %>%
      purrr::reduce(inner_join, by = c("geneID", "DEG_Desc")) %>%
      dplyr::select(geneID)
    
    filtered.list.highlight.tbl_LS <- lapply(names(download_DT), function(y){
      dplyr::filter(download_DT[[y]], geneID %in% subsetIDs$geneID)%>%
        dplyr::arrange(desc(logFC))
    })
    names(filtered.list.highlight.tbl_LS) <- names(download_DT)
  } else {
    filtered.list.highlight.tbl_LS <-lapply(names(download_DT), function(y){
      dplyr::filter(download_DT[[y]],str_detect(DEG_Desc,paste0(input$download_DGEdt_direction_LS, collapse = "|"))) %>%
        dplyr::arrange(desc(logFC))
    })
    names(filtered.list.highlight.tbl_LS)<- names(download_DT)
  }
  
  # Pass only a specific proportion of genes. Remember, the datatable is grouped by the DEG_Desc value, and ordered by descending logFC value. 
  filtered.list.highlight.tbl_LS<-lapply(names(filtered.list.highlight.tbl_LS), function(y){
    
    filtered.list.highlight.tbl_LS[[y]] %>%
      dplyr::group_map(~ {
        if (str_detect(.y, "Up")) {
          slice_max(.x, order_by = logFC, prop = as.numeric(input$percentDGE_LS)/100)
        } else if (str_detect(.y, "NotSig")) {
          slice_max(.x, order_by = logFC, prop = as.numeric(input$percentDGE_LS)/100)
        } else if (str_detect(.y, "Down")) {
          slice_min(.x, order_by = logFC, prop = as.numeric(input$percentDGE_LS)/100)
        }
      }, .keep = TRUE) %>%
      bind_rows() %>%
      dplyr::arrange(desc(logFC))
  })
  names(filtered.list.highlight.tbl_LS)<- names(download_DT)
  filtered.list.highlight.tbl_LS
})

## LS: Save Excel Tables with DEG Tables ----
output$downloadbuttonLS <- renderUI({
  req(input$goLS,vals$comparison_LS)
  filtered.tbl_LS <- filter_DEG_tbl_LS()
  
  expressionnotes <- "Columns labeled with <life stage - sample ID> report log2 counts per million (CPM) expression. Columns labeled avg_<life stage> are median log2CPM with IQR. The column labeled logFC reports log2 fold change."
  ### Generate some text that describes the analysis conditions
  ### 1. If p-values are adjusted for multiple pairwise comparisons, which comparisons are included in the adjustment parameters? This should be the list of selected contrasts
  if (vals$multipleCorrection_LS == TRUE){
    multiplecorrection <- paste0(vals$comparison_LS, collapse = "; ")
    multiplecorrection <- paste0("P-values corrected across the following multiple pairwise comparisons: ", multiplecorrection)
  } else {
    multiplecorrection <- "P-values *not* corrected across multiple pairwise comparisons"
  }
  ### 2. If downloading results are being filtered to show only genes that display consistent differential expression across all comparisons targeted for download, list the pairwise comparisons being used.
  if (input$download_DGEdt_across_LS == TRUE){
    filteredacross <- paste0(names(filtered.tbl_LS), collapse = "; ")
    filteredacross <- paste0("Lists only include genes with matching differential expression descriptions across the following pairwise comparisons: ", filteredacross)
  } else filteredacross <- ""
  
  ### 3. Specify which types of differential expression pattern
  DEGpattern <- paste0(input$download_DGEdt_direction_LS, collapse = "; ")
  DEGpattern <- paste0("Lists include genes that are differentially regulated in the following directions (relative to the target life stage): ", DEGpattern)
  ### 4. Specify the proportion of genes for each DEG Type that are being saved (e.g. top 10% of upregulated, and top10% of downregulated genes)
  proportionexport <- paste0("Percentage of genes for each differential expression pattern: ", input$percentDGE_LS)
  output$generate_excel_report_LS <- generate_excel_report(names(filtered.tbl_LS), 
                                                           filtered.tbl_LS,
                                                           name = paste(input$selectSpecies_LS,
                                                                        "RNA-seq Differential Gene Expression"),
                                                           filename_prefix = 'DGETable_',
                                                           expressionnotes = expressionnotes,
                                                           multiplecorrection = multiplecorrection,
                                                           filteredacross = filteredacross,
                                                           DEGpattern = DEGpattern,
                                                           proportionexport = proportionexport)
  
  downloadButton("generate_excel_report_LS",
                 "Download DGE Tables as Excel",
                 class = "btn-primary")
  
})