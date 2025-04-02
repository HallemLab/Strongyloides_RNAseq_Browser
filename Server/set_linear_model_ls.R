## LS: Set Contrast Matrix, Fit the Linear Model, Extract the Differentially Expressed Genes ----
## This reactive function is called by app.R
set_linear_model_LS <- eventReactive(input$goLS,{
    req(vals$limmacontrast_LS)
    limma_ranking(vals$limmacontrast_LS, 
                  vals$targetStage_LS, 
                  vals$contrastStage_LS, 
                  vals$multipleCorrection_LS, 
                  NA, 
                  vals, vals$fit, 
                  vals$v.DGEList.filtered.norm, 
                  input$adj.P.thresh_LS,
                  input$lfc.thresh_LS, 
                  vals$diffGenes.df)
    vals$list.highlight.tbl_LS <- vals$list.highlight.tbl
})