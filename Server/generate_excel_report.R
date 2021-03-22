generate_excel_report <- function(comparison, tbl,
                                  name = "RNA-seq Differential Gene Expression",
                                  filename_prefix = "Stercoralis_RNAseq_Data_",
                                  subtitle_prefix = "Contrast:",
                                  expressionnotes = "",
                                  multiplecorrection = "",
                                  filteredacross = "",
                                  DEGpattern = "",
                                  proportionexport = "",
                                  n_header_rows = 10){
  
  temp <- downloadHandler(
    
    filename = function(){
      paste(filename_prefix,Sys.Date(),".xlsx",sep = "")
    },
    
    content = function(file){
      
      withProgress({
        removeModal() 
        
        # Generate Workbook ----
        to_download <<- createWorkbook()
        
        setProgress(.25)
        
        # Write Data to Workbook ----
        
        sapply(seq_along(comparison), function(y){
          
          # Create Tab
          addWorksheet(wb = to_download, sheetName = paste0("Sheet_",y))
          
          ## Write Header to Sheet
          writeData(
            to_download,
            sheet = y,
            x = c(
              name,
              paste(subtitle_prefix, comparison[[y]]),
              paste0("Report generated on ", format(Sys.Date(), "%B %d, %Y")),
              expressionnotes,
              multiplecorrection,
              DEGpattern,
              proportionexport,
              filteredacross
            )
          )
        
          ## Write Data to Sheet
          writeData(
            to_download,
            sheet = y,
            x = tbl[[y]],
            startRow = n_header_rows,
            startCol = 1,
            headerStyle = createStyle(
              textDecoration = "Bold",
              halign = "center",
              border = "bottom"
            ),
            keepNA = TRUE,
            na.string = "NA"
          )
          
          # Implement Styling ----
          ## Styling the title row
          addStyle(
            to_download,
            sheet = y,
            rows = 1,
            cols = 1:10,
            style = createStyle(
              fontSize = "14",
              textDecoration = "bold"
            )
          )
          
        })
        setProgress(.75)
        saveWorkbook(to_download, file)
        setProgress(1)
      },
      message = "Generating Excel Report"
      )
    }
  )
}