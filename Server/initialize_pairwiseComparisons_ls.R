## LS: Clear Comparison Selections ----
observeEvent({input$resetLS
    input$tab},{
        updateSelectInput(session, "selectContrast_LS", selected = "")
        updateSelectInput(session, "selectTarget_LS", selected = "")
        updateTextAreaInput(session,"multiContrasts_LS",value = "")
    })

## LS: Generate Comparison Selection Boxes
output$pairwiseSelector_LS<- renderUI({
    req(vals$v.DEGList.filtered.norm)
    
    list(
        panel(
            heading = tagList(h5(shiny::icon("fas fa-sliders-h"), 
                                 "Step 1: Pick Life Stage Comparisons")),
            status = "primary",
            p(tags$em('Users may set up a single pairwise comparison for differential gene analysis using the selection menus below. Alternatively, users may also input multiple pairwise comparisons by typing a comma-separated list of contrasts into the textbox.', style = "color: #7b8a8b")),
            h5('A: Single Comparison', class = 'text-danger', style = "margin: 0px 0px 10.5px 0px"),
            selectInput("selectTarget_LS",
                        h6("Select Target"),
                        choices = c('Choose one or more' = ''
                                    ,as.list(levels(vals$v.DEGList.filtered.norm$targets$group))),
                        selectize = TRUE,
                        multiple = TRUE),
            # Select Contrast Life Stage
            selectInput("selectContrast_LS",
                        h6("Select Contrast"),
                        choices = c('Choose one or more' = ''
                                    ,'All Others (separately)' = 'everythingElse'
                                    ,'All Others (one group)' = 'remainingGroup'
                                    ,as.list(levels(vals$v.DEGList.filtered.norm$targets$group))),
                        selectize = TRUE,
                        multiple = TRUE),
            tags$hr(style="border-color: black;"),
            p(tags$em('If using the textbox to type comma-separated contrasts, please use the format: (Target)-(Contrast). For example typing "iL3-PF, (iL3+iL3a)-(PF+FLF)" will run two pairwise comparisons: iL3 vs PF and iL3+iL3a vs PF+FLF. To correct for multiple, closely related pairwise comparisons, use the toggle switch below.', style = "color: #7b8a8b")),
            h5('B: Multiple Comparisons', class = 'text-danger'),
            # Text Input for Multiple Contrasts
            textAreaInput('multiContrasts_LS',
                          NULL,
                          rows = 5, 
                          resize = "vertical"),
            
            h6("Correct for multiple comparisons?"),
            
            checkboxInput("multipleContrastsYN_LS",
                          p("Yes, correct p-values for multiple pairwise comparisons")),
            tags$hr(style="border-color: #2C3E50;"),
            
            ### Action Button
            actionButton('goLS',
                         'Submit',
                         icon = icon("fas fa-share"),
                         class = "btn-primary"),
            
            actionButton('resetLS', 'Clear',
                         icon = icon("far fa-trash-alt"))
            
        )
    )
})

## LS: Update checkbox input for multiple contrasts depending on contrast selection
### specifically, if users select the All Others (separately) option, the multiple contrasts
### checkbox should be updated to TRUE. Users may still de-select this option.
observe({
    req(input$selectContrast_LS)
    multiY <- any(str_detect(input$selectContrast_LS, 'everythingElse'))
    updateCheckboxInput(session,
                        "multipleContrastsYN_LS",
                        value = multiY)
})

## LS: Pairwise comparisons Across Life Stages ----
parse_contrasts_LS<-eventReactive(input$goLS,{
    # Make sure there are sufficient inputs
    shiny::validate(
        shiny::need((isTruthy(input$multiContrasts_LS) || (isTruthy(input$selectTarget_LS) && isTruthy(input$selectContrast_LS))),
                    "Not enough inputs were provided to generate contrasts. Please re-select inputs.")
    )
    
    if (isTruthy(input$multiContrasts_LS)) {
        comparison <- input$multiContrasts_LS %>%
            gsub(" ", "", ., fixed = TRUE) %>%
            str_split(pattern = "[,;]") %>%
            unlist()
        
        # Validation check - does the matrix of comparisons contain info for targets and contrasts?
        # This will fail only if a single incomplete pairwise comparison is provided 
        # using the multiContrasts input box
        shiny::validate(
            shiny::need({
                (comparison %>%
                     str_split(pattern="-", simplify = T) %>%
                     ncol()) == 2
                
            },
            message = "Not enough inputs were provided to generate contrasts. Please re-type inputs.")) 
        
        targetStage<- comparison %>%
            str_split(pattern="-", simplify = T) %>%
            .[,1] %>%
            gsub("(", "", ., fixed = TRUE) %>%
            gsub(")", "", ., fixed = TRUE) %>%
            str_split(pattern = "\\+", simplify = T)
        contrastStage<-comparison %>%
            str_split(pattern="-", simplify = T) %>%
            .[,2] %>%
            gsub("(", "", ., fixed = TRUE) %>%
            gsub(")", "", ., fixed = TRUE) %>%
            str_split(pattern = "\\+", simplify = T)
        comparison<- sapply(seq_along(comparison),function(x){
            tS <- as.vector(targetStage[x,]) %>%
                .[. != ""] 
            cS <- as.vector(contrastStage[x,]) %>%
                .[. != ""] 
            paste(paste0(tS, 
                         collapse = "+") %>%
                      paste0("(",.,")/",length(tS)),
                  paste0(cS, 
                         collapse = "+") %>%
                      paste0("(",.,")/",length(cS)),
                  sep = "-")
            
        })
        if (input$multipleContrastsYN_LS == T) {
            vals$multipleCorrection_LS <- T
        } else vals$multipleCorrection_LS <- F
    } else if (str_detect(input$selectContrast_LS[[1]], 'everythingElse')){
        targetStage <- rbind(input$selectTarget_LS)
        contrastStage <- setdiff(levels(vals$v.DEGList.filtered.norm$targets$group),targetStage) %>%
            cbind()
        targetStage <- rep_len(targetStage,length(contrastStage)) %>%
            cbind()
        comparison <- sapply(seq_along(contrastStage), function(x){
            paste(paste0(input$selectTarget_LS, 
                         collapse = "+") %>%
                      paste0("(",.,")/",length(input$selectTarget_LS)),
                  contrastStage[[x]] %>%
                      paste0("(",.,")/1"),
                  sep = "-"
            )
        })
        if (input$multipleContrastsYN_LS == T) {
            vals$multipleCorrection_LS <- T
        } else vals$multipleCorrection_LS <- F
    } else if (str_detect(input$selectContrast_LS[[1]], 'remainingGroup')){
        targetStage <- rbind(input$selectTarget_LS)
        contrastStage <- setdiff(levels(vals$v.DEGList.filtered.norm$targets$group),targetStage) %>%
            rbind()
        comparison <- paste(paste0(targetStage, 
                                   collapse = "+") %>%
                                paste0("(",.,")/",length(targetStage)),
                            paste0(contrastStage, 
                                   collapse = "+") %>%
                                paste0("(",.,")/",length(contrastStage)),
                            sep = "-"
        )
        
        if (input$multipleContrastsYN_LS == T) {
            vals$multipleCorrection_LS <- T
        } else vals$multipleCorrection_LS <- F
        
    } else {
        targetStage <- rbind(input$selectTarget_LS)
        contrastStage <- rbind(input$selectContrast_LS)
        comparison <- paste(paste0(input$selectTarget_LS, 
                                   collapse = "+") %>%
                                paste0("(",.,")/",length(input$selectTarget_LS)), 
                            paste0(input$selectContrast_LS, 
                                   collapse = "+") %>%
                                paste0("(",.,")/",length(input$selectContrast_LS)), 
                            sep = "-")
        if (input$multipleContrastsYN_LS == T) {
            vals$multipleCorrection_LS <- T
        } else vals$multipleCorrection_LS <- F
    }
    
    # Validation checks for contrast inputs:
    # 1. Are contrasts complete and formatted correctly (are there life stages on each side of the '-' side)
    ## Check for missing "target" elements 
    missing.targets <- sapply(seq(1,nrow(targetStage)), function(x){
        row <- x
        targets.clean <- targetStage[row, targetStage[row,] != ""] %>%
            is_empty()
    })
    
    shiny::validate(shiny::need(any(!missing.targets), 
                                message = "At least one pairwise comparison is missing a target element. Please check inputs.")) 
    
    ## Check for missing "contrast" elements
    missing.contrasts <-sapply(seq(1,nrow(targetStage)), function(x){
        row <- x
        contrasts.clean <- contrastStage[row, contrastStage[row,] != ""] %>%
            is_empty()
    }) 
    shiny::validate(shiny::need(any(!missing.contrasts), 
                                message = "At least one pairwise comparison is missing a contrast element. Please check inputs.")) 
    
    # 2. Do contrasts include recognized life stages (corrects for spelling mistakes); compare relative to abbreviated names in lifestage_legend
    ## Do all target elements match a life stage in this dataset? 
    error.targets.validNames <- targetStage[targetStage != ""] %in% 
        levels(vals$v.DEGList.filtered.norm$targets$group) 
    shiny::validate(shiny::need(all(error.targets.validNames), 
                                message = "At least one target name doesn't match available life stages. Please check inputs for spelling mistakes or incorrect capitalization.")) 
    
    ## Do all contrast elements match a life stage in this dataset?
    error.contrast.validNames <- contrastStage[contrastStage != ""] %in% 
        levels(vals$v.DEGList.filtered.norm$targets$group) 
    shiny::validate(shiny::need(all(error.contrast.validNames), 
                                message = "At least one contrast name doesn't match available life stages. Please check inputs for spelling mistakes or incorrect capitalization.")) 
    
    # 3. Are there repeated elements within a contrast.
    ## Produce error message if there are repeated elements in any of the target and contrast pairings. Using *apply and vectorized string detection to search rows of target/contrast stage arrays for matching elements located in the columns.
    error.matches <- sapply(seq(1,nrow(targetStage)), function(x){
        tar.row <- x
        targets.clean <- targetStage[tar.row, targetStage[tar.row,] != ""] #remove any blank values
        sapply(seq(1,length(targets.clean)), function(x, tar.row, targets.clean){
            tar.col <- x
            contrasts.clean <- contrastStage[tar.row, (contrastStage[tar.row,] != "")] #remove any blank values
            targetStage[tar.row, tar.col]
            str_detect(contrasts.clean, paste0('\\b',targets.clean[tar.col],'\\b'))
        }, tar.row, targets.clean) %>%
            any()
    })
    
    shiny::validate(shiny::need(!any(error.matches), 
                                message = "There are repeated elements in at least one target/contrast pairing. Please input new pairings.")) 
    
    
    vals$targetStage_LS <- targetStage 
    vals$contrastStage_LS <- contrastStage
    vals$limmacontrast_LS <- comparison
    vals$comparison_LS <- gsub("/[0-9]*","", comparison) %>%
        gsub("\\(|\\)","",.)
    
})

## LS: Generate Responsive Selection for Life Stage to Display ----
output$contrastDisplaySelection_LS <- renderUI({
    req(vals$v.DEGList.filtered.norm)
    input$resetLS
    input$speciesLS
    
    comparison <- parse_contrasts_LS()
    isolate({
        panel(
            heading = tagList(h4(shiny::icon("fas fa-filter"),
                                 "Pick Contrast to Display")),
            status = "danger",
            pickerInput("displayedComparison_LS",
                        NULL, 
                        comparison,
                        options = list(style = 'btn btn-danger'))
        )
    })
})



## LS: Generate Legend Explaining the Life Stages ----
output$lifeStageLegend_LS <- renderDT({
    req(vals$v.DEGList.filtered.norm)
    lifestage_legend.df <- lifestage_legend %>%
        dplyr::select(any_of(unique(vals$v.DEGList.filtered.norm$targets$group))) %>%
        as.data.frame()
    rownames(lifestage_legend.df)<- c("<b>Life Stage</b>")
    lifestage_legend.DT <- lifestage_legend.df %>%
        DT::datatable(rownames = TRUE,
                      escape = FALSE,
                      class = "table-bordered",
                      colnames = c("Abbr." = 1),
                      options = list(scrollX = TRUE,
                                     ordering = FALSE,
                                     dom = 'tS'
                      )
        )
    lifestage_legend.DT
})

output$Legend_LS <- renderUI({
    req(vals$v.DEGList.filtered.norm)
    panel(
        heading = tagList(h5(shiny::icon("fas fa-book-open"),
                             "Life Stage Legend")),
        status = "default",
        DTOutput("lifeStageLegend_LS")
    )
})