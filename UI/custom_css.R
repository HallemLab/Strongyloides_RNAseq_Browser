## Custom CSS ----
tags$head(
    tags$style(HTML("
    h3 {
    font-size: 16px;
    margin: 10.5px 0px;
    }
    
    .navbar-brand {
    height: 60px;
    padding: 10px 15px;
    }
    
    .navbar-nav>li>a{
    height: 60px;
    padding: 10px 15px;
    }
    
    .navbar-text {
    float: right;
    }
    
    h4 {
    font-size: 15px;
    }
    
    h5 {
    font-size: 15px;
    margin: 5px 0px;
    }
    
    h6 {
    font-size: 13px;
    margin: 0px;
    font-weight: 550;
    line-height: 1.4;
    }
    
    strong {
    font-size: 13px;
    font-weight: bold;
    }
    
    p{
    font-size: 13px;
    font-weight: normal;
    }
    
    ol{
    font-size: 13px
    }
    
    li{
    font-size: 13px
    }
    
    .selectize-input {
    word-wrap: break-word;
    font-size: 13px;
    overflow-x: auto;
    }
    
    .selectize-dropdown {
    word-wrap: break-word;
    font-size: 13px;
    }
    
    .form-control {
    font-size: 13px;
    height: 40px;
    }
    
    .btn {
    font-size: 13px;
    height: 40px;
    }
    
    #CPMPlotlydiv {
    text-align: center;
    color: black;
    }
    
    .shiny-output-error-validation {
    font-size: 15px;
    color: #E74C3C
    }
    
    #geneSelection_conditionalPanel .shiny-output-error-validation {
    color: white;
    }
    
    #contrastDisplaySelectionPanel_LS .shiny-output-error-validation {
    color: white;
    }
    
    #contrastSelectionPanel_GW .shiny-output-error-validation {
    color: white;
    }
    
    #lifeStageLegend_GW .shiny-html-output {
    font-size: 13px}
    
    #lifeStageLegend_LS .shiny-html-output {
    font-size: 13px}
    
    #speciesPanelID_GW .selectize-input {
    font-style: italic;
    }
    
    #speciesPanelID_GW .selectize-dropdown-content {
    font-style: italic;
    }
    
    #speciesPanelID_LS .selectize-input {
    font-style: italic;
    }
    
    #speciesPanelID_LS .selectize-dropdown-content {
    font-style: italic;
    }
    
    #About_Overview h3 {
    font-size: 15px;
    text-decoration: underline;
    font-weight: bold;
    }
    
     #About_Updates h3 {
    font-size: 15px;
    font-weight: bold;
    }
    
    #About_Preprocess h3 {
    font-size: 15px;
    text-decoration: underline;
    font-weight: bold;
    }
    
    #About_Data h3 {
    font-size: 15px;
    text-decoration: underline;
    font-weight: bold;
    }
    
    #About_Analysis_Methods h3 {
    font-size: 15px;
    text-decoration: underline;
    font-weight: bold;
    }
    
    .nav-pills li {
    font-size: 14px;
    font-weight: bold;
    }
    
    .nav-pills {
    border-bottom: 1px solid #2C3E50;
    }
    
    #download_DGEdt_direction_LS .shiny-options-group{
    font-size: 13px;
    }
    
    #download_DGEdt_direction_GW .shiny-options-group{
    font-size: 13px;
    }

                    "))
    
)