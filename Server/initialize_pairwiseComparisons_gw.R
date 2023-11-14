## GW: Clear Comparison Selections ----
observeEvent({input$resetGW
    input$tab},{
        updateSelectInput(session, "selectContrast_GW", selected = "")
        updateSelectInput(session, "selectTarget_GW", selected = "")
        updateTextAreaInput(session,"multiContrasts_GW",value = "")
        
    })

## GW: Generate Comparison Selection Boxes ----
output$pairwiseSelector_GW<- renderUI({
    req(vals$v.DEGList.filtered.norm)
    
    list(
        panel(
            # Select Life Stage Comparisons 
            heading = tagList(h5(shiny::icon("fas fa-sliders-h"), 
                                 "Step 2: Pick Life Stage Comparisons")),
            status = "primary",
            p(tags$em('Users may set up a single pairwise comparison for differential gene analysis using the selection menus below. Alternatively, users may also input multiple pairwise comparisons by typing a comma-separated list of contrasts into the textbox.', style = "color: #7b8a8b")),
            # Select Target Life Stage
            h5('A: Single Comparison', class = 'text-danger'),
            
            selectInput("selectTarget_GW",
                        h6("Select Target"),
                        choices = c('Choose one or more' = ''
                                    ,as.list(levels(vals$target.contrast.options))),
                        selectize = TRUE,
                        multiple = TRUE),
            # Select Contrast Life Stage
            selectInput("selectContrast_GW",
                        h6("Select Contrast"),
                        choices = c('Choose one or more' = ''
                                    ,'All Others (separately)' = 'everythingElse'
                                    ,'All Others (one group)' = 'remainingGroup'
                                    ,as.list(levels(vals$target.contrast.options))),
                        selectize = TRUE,
                        multiple = TRUE),
            
            
            tags$hr(style="border-color: #2C3E50;"),
            h5('B: Multiple Comparisons', class = 'text-danger'),
            p(tags$em('If using the textbox to type comma-separated contrasts, please use the format: (Target)-(Contrast). For example typing "iL3-PF, (iL3+iL3a)-(PF+FLF)" will run two pairwise comparisons: iL3 vs PF and iL3+iL3a vs PF+FLF. To correct for multiple, closely related pairwise comparisons, use the toggle switch below.', style = "color: #7b8a8b")),
            # Text Input for Multiple Contrasts
            textAreaInput('multiContrasts_GW',
                          NULL,
                          rows = 5, 
                          resize = "vertical"),
            
            h6("Correct for multiple comparisons?"),
            checkboxInput("multipleContrastsYN_GW",
                          p("Yes, correct p-values for multiple pairwise comparisons")),
            
            tags$hr(style="border-color: #2C3E50;"),
            
            ### Action Button
            actionButton('goLifeStage_GW',
                         'Process',
                         icon = icon("fas fa-share"),
                         class = "btn-primary"),
            
            actionButton('resetGW', 'Clear',
                         icon = icon("far fa-trash-alt"))
        )
        
    )
})

## GW: Update checkbox input for multiple contrasts depending on contrast selection
### specifically, if users select the All Others (separately) option, the multiple contrasts
### checkbox should be updated to TRUE. Users may still de-select this option.
observe({
    req(input$selectContrast_GW)
    multiY <- any(str_detect(input$selectContrast_GW, 'everythingElse'))
    updateCheckboxInput(session,
                        "multipleContrastsYN_GW",
                        value = multiY)
})

## GW: Parse the inputs ----
parse_contrasts_GW <- eventReactive(input$goLifeStage_GW,{
    # Make sure there are sufficient inputs
    validate(
        need((isTruthy(input$multiContrasts_GW) || (isTruthy(input$selectTarget_GW) && isTruthy(input$selectContrast_GW))),
             "Not enough inputs were provided to generate contrasts. Please re-select.")
    )
    if (isTruthy(input$multiContrasts_GW)) {
        comparison <- input$multiContrasts_GW %>%
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
            gsub(")", "", ., fixed = TRUE)  %>%
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
        if (input$multipleContrastsYN_GW == T) {
            vals$multipleCorrection_GW <- T
        } else vals$multipleCorrection_GW <- F
    } else if (str_detect(input$selectContrast_GW[[1]], 'everythingElse')){
        targetStage <- rbind(input$selectTarget_GW)
        contrastStage <- setdiff(levels(vals$v.DEGList.filtered.norm$targets$group),targetStage) %>%
            cbind()
        targetStage <- rep_len(targetStage,length(contrastStage)) %>%
            cbind()
        comparison <- sapply(seq_along(contrastStage), function(x){
            paste(paste0(input$selectTarget_GW, 
                         collapse = "+") %>%
                      paste0("(",.,")/",length(input$selectTarget_GW)),
                  contrastStage[[x]] %>%
                      paste0("(",.,")/1"),
                  sep = "-"
            )
        })
        if (input$multipleContrastsYN_GW == T) {
            vals$multipleCorrection_GW <- T
        } else vals$multipleCorrection_GW <- F
    } else if (str_detect(input$selectContrast_GW[[1]], 'remainingGroup')){
        targetStage <- rbind(input$selectTarget_GW)
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
        
        if (input$multipleContrastsYN_GW == T) {
            vals$multipleCorrection_GW <- T
        } else vals$multipleCorrection_GW <- F
    } else {
        targetStage <- rbind(input$selectTarget_GW)
        contrastStage <- rbind(input$selectContrast_GW)
        comparison <- paste(paste0(input$selectTarget_GW, 
                                   collapse = "+") %>%
                                paste0("(",.,")/",length(input$selectTarget_GW)), 
                            paste0(input$selectContrast_GW, 
                                   collapse = "+") %>%
                                paste0("(",.,")/",length(input$selectContrast_GW)), 
                            sep = "-")
        if (input$multipleContrastsYN_GW == T) {
            vals$multipleCorrection_GW <- T
        } else vals$multipleCorrection_GW <- F
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
    
    vals$targetStage_GW <- targetStage 
    vals$contrastStage_GW <- contrastStage
    vals$limmacontrast_GW <- comparison
    vals$comparison_GW <- gsub("/[0-9]*","", comparison) %>%
        gsub("\\(|\\)","",.)
    
})

## GW: Generate Responsive Selection for Life Stage to Display ----
output$contrastDisplaySelection_GW <- renderUI({
    comparison <- parse_contrasts_GW()
    req(vals$genelist)
    isolate({
        panel(
            heading = tagList(h5(shiny::icon("fas fa-filter"),
                                 "Pick Contrast to Display")),
            status = "danger",
            pickerInput("displayedComparison_GW",
                        NULL, 
                        comparison,
                        selected = comparison[[1]],
                        options = list(style = 'btn btn-danger'))
        )
    })
})