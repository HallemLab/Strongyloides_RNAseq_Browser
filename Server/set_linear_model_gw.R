## GW: Set Contrast Matrix, Fit the Linear Model, Extract the Differentially Expressed Genes ----
## This reactive function is called by app.R
set_linear_model_GW <- eventReactive(input$goLifeStage_GW,{
    req(vals$limmacontrast_GW)
    limma_ranking(vals$limmacontrast_GW, 
                  vals$targetStage_GW, 
                  vals$contrastStage_GW, 
                  vals$multipleCorrection_GW, 
                  vals$genelist, 
                  vals, vals$fit, 
                  vals$v.DGEList.filtered.norm, 
                  input$adj.P.thresh_GW,
                  input$lfc.thresh_GW,
                  vals$diffGenes.df)
    
    vals$list.myTopHits.df_GW <- vals$list.myTopHits.df
    vals$list.highlight.tbl_GW <- vals$list.highlight.tbl
    
})