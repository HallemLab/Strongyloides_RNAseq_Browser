## GW: Generate Main Gene Input Panel ----
output$genePanelinputs <- renderUI({
    tagList(
        panel(
            id = "GeneInputBox",
            heading = tagList(h5(shiny::icon("fas fa-dna"), "Step 1: Input Genes / Keywords")),
            status = "primary",
            ### GeneID (text box)
            h5('Pick Genes', class = 'text-danger', style = "margin: 0px 0px 5px 0px"),
            p(tags$em('Users may type gene stable IDs or keywords that will be matched against WormBase ParaSite Gene Descriptions, known C. elegans homologs, homologs from closely-related Strongyloides species, InterPro terms, and an Ensembl Compara database of gene families. Search terms may be separated using commas, semicolons, or new lines. Users may also upload a .csv file containing search terms.', style = "color: #7b8a8b")),
            p(tags$em("Type 'everything' or 'all genes' to display all genes in the genome. Warning: this will take a long time to process.", style = "color: #7b8a8b")),
            p(tags$em(tags$b('Note: Please hit the Clear button if switching between typing and uploading inputs.', style = "color: #F39C12"))),
            textAreaInput('idtext',
                          h6('Gene IDs or Keywords'),
                          rows = 5, 
                          resize = "vertical"),
            
            ### Upload list of GeneIDs
            uiOutput('genefile_upload'),
            
            ### Action Button
            actionButton('goGW',
                         'Submit',
                         icon = icon("fas fa-share"),
                         class = "btn-primary"),
            
            actionButton('resetGenes', 'Clear',
                         icon = icon("far fa-trash-alt"))
        )
    )
})

## GW: Parse Gene Inputs ----
parse_ids <- eventReactive(input$goGW,{
    vals$genelist <- NULL
    vals$HeatmapRowOrder <- NULL
    validate(
        need(isTruthy(vals$v.DEGList.filtered.norm), 
             "Please re-select a species for analysis")
    )
    
    validate(
        need({isTruthy(input$loadfile) | isTruthy(input$idtext)}, "Please input gene ids")
    )
    isolate({
        withProgress({
        if (isTruthy(input$idtext)){
            terms <- input$idtext %>%
                gsub("\\n",",",.) %>% #replace any new lines with commas
                trimWhiteSpace %>% #remove leading and trailing white psace from string
                str_split(pattern = ",|;") %>%
                unlist()
        } else if (isTruthy(input$loadfile)){
            file <- input$loadfile
            ext <- tools::file_ext(file$datapath)
            validate(need(ext == "csv", "Please upload a csv file"))
            suppressWarnings(
                terms <- read.csv(file$datapath, 
                                  header = FALSE, 
                                  colClasses = "character", 
                                  strip.white = T,
                                  na.strings = "") %>%
                    as_tibble() %>%
                    
                    pivot_longer(cols = everything(), 
                                 values_to = "geneID",
                                 values_drop_na = T) 
            )
            terms <- terms$geneID
            
        } 
        setProgress(0.1)
        
            genelist <- vals$annotations %>%
            dplyr::select(geneID)
        
        if (any(grepl('everything|all genes', terms, ignore.case = TRUE))) {
            # Text input matches the strings 'everything' or 'all genes'
            genelist <- genelist
        } else {
            inc <- 0.1/nrow(terms)
            # Search for gene IDs
            terms.cleaned <- gsub("\\.[0-9]$","",terms) #strip any transcript values from the inputed list
            
            geneindex.geneID<-sapply(terms.cleaned, function(y) {
                incProgress(amount = inc)
                grepl(gsub("^\\s+|\\s+$", "", y), #remove any number of whitespace from start or end
                      vals$annotations$geneID,
                      ignore.case = TRUE)
            }) %>%
                rowSums() %>%
                as.logical()
            
            
            # Search WormBase Parasite Gene Description Terms
            geneindex.description<-sapply(terms, function(y) {
                incProgress(amount = inc)
                grepl(gsub("^\\s+|\\s+$", "", y), #remove any number of whitespace from start or end
                      vals$annotations$Description,
                      ignore.case = TRUE)
            }) %>%
                rowSums() %>%
                as.logical()
            
            # Search Parasite Ensembl Compara Gene Sets
            ensComp<- ensComp %>%
                left_join(genelist, ., by = "geneID") %>%
                dplyr::relocate(gs_name, geneID)
            geneindex.ensembl<-sapply(terms, function(y) {
                incProgress(amount = inc)
                gsub("^\\s+|\\s+$", "", y) %>%
                paste0("\\<",.,"\\>") %>%
                    grepl(., 
                          ensComp$gs_name,
                          ignore.case = TRUE)
            }) %>%
                rowSums() %>%
                as.logical()
            
            # Search C. elegans homolog IDs
            geneindex.Cehomologs<-sapply(terms, function(y) {
                incProgress(amount = inc)
                gsub("^\\s+|\\s+$", "", y) %>%
                paste0("\\<",.,"\\>") %>%
                    grepl(., 
                          vals$annotations$Ce_geneID,
                          ignore.case = TRUE)
            }) %>%
                rowSums() %>%
                as.logical()
            
            # Search In-group homolog IDs
            geneindex.InGroup<-sapply(terms, function(y) {
                incProgress(amount = inc)
                gsub("^\\s+|\\s+$", "", y) %>%
                    paste0("\\<",.,"\\>") %>%
                    grepl(., 
                          vals$annotations$In.subclade_geneID,
                          ignore.case = TRUE)
            }) %>%
                rowSums() %>%
                as.logical()
            
            # Search Out-group homolog IDs
            geneindex.OutGroup<-sapply(terms, function(y) {
                incProgress(amount = inc)
                gsub("^\\s+|\\s+$", "", y) %>%
                    paste0("\\<",.,"\\>") %>%
                    grepl(., 
                          vals$annotations$Out.subclade_geneID,
                          ignore.case = TRUE)
            }) %>%
                rowSums() %>%
                as.logical()
            
            # Search InterPro Terms
            geneindex.InterPro<-sapply(terms, function(y) {
                incProgress(amount = inc)
                gsub("^\\s+|\\s+$", "", y) %>%
                    paste0("\\<",.,"\\>") %>%
                    grepl(., 
                          vals$annotations$InterPro,
                          ignore.case = TRUE)
            }) %>%
                rowSums() %>%
                as.logical()
            
            geneindex <- geneindex.geneID | geneindex.description | geneindex.ensembl | geneindex.Cehomologs | geneindex.InGroup |geneindex.OutGroup | geneindex.InterPro 
            genelist <- dplyr::filter(genelist,geneindex) 
        }
       
        if (nrow(genelist) == 0){
            disable("goLifeStage_GW")
        } else {enable("goLifeStage_GW")}
        
        # Produces error message if genelist is empty
        validate(
            need(nrow(genelist) != 0, "RNA-seq data unavailable for submitted genes (submitted names may be invalid). Please try a new search.")
        )
       
        # Save record of original genelist before filtering, removing rows that contain the word 'gene'
        vals$submitted.genelist <- genelist %>%
            dplyr::filter(!grepl("gene", geneID, ignore.case = T))
        
        # Remove genes from the list that aren't part of vals$Log2CPM
        genelist <- genelist %>%
            dplyr::filter(geneID %in% vals$Log2CPM$geneID)
        
        setProgress(0.95)
        
        if (nrow(genelist) == 0){
            disable("goLifeStage_GW")
        } else {enable("goLifeStage_GW")}
        
        # Produces error message if genelist is empty after removing genes not included in the RNAseq dataset
        validate(
            need(nrow(genelist) != 0, "RNA-seq data unavailable for submitted genes. Please try a new search.")
        )
        
        vals$genelist <- genelist
        vals$genelist.Log2CPM <- vals$Log2CPM %>%
            dplyr::filter(geneID %in% genelist$geneID)
        
        },message = "Parsing gene IDs...")
    })
})