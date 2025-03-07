# DE Analysis Server


################################ Helper Functions ############################

# javascript code to collapse box
jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"

################################## Observer ##################################

observeEvent(input$deMultTest, {
  if(input$deMultTest){
    shinyjs::show(id = "deMultTestMethod")
  } else {
    shinyjs::hide(id = "deMultTestMethod")
  }
})

observeEvent(input$deLogFC, {
  if(input$deLogFC){
    shinyjs::show(id = "deLogFCThr")
  } else {
    shinyjs::hide(id = "deLogFCThr")
  }
})

observeEvent(input$deMethodInput, {
  if(input$deMethodInput == "ROTS"){
    shinyjs::show(id = "deROTSK")
    shinyjs::show(id = "deROTSB")
  } else {
    shinyjs::hide(id = "deROTSK")
    shinyjs::hide(id = "deROTSB")
  }
  if(input$deMethodInput == "DEqMS"){
    shinyjs::show(id = "deDEqMSColumn")
  } else {
    shinyjs::hide(id = "deDEqMSColumn")
  }
  if(input$deMethodInput == "limma"){
    shinyjs::show(id = "deLimmaTrend")
    shinyjs::show(id = "deLimmaRobust")
    
  } else {
    shinyjs::hide(id = "deLimmaTrend")
    shinyjs::hide(id = "deLimmaRobust")
  }
})

#observeEvent(input$deBiomarkerColumn, {
#  column <- input$deBiomarkerColumn
#  rd <- as.data.table(rowData(reactiveVals$se))
#  ids <- rd[[column]]
#  ids <- ids[ids != ""]
#  ids <- unlist(strsplit(ids, ";"))
#  ids <- unique(ids)
#  updatePickerInput(session = session, inputId =   "deBiomarkerSelection", choices = ids , selected = ids[1:10])
#})

observeEvent(input$performDEAnalysis, {
  waiter_show(id = "app",html = tagList(spinner$logo,
                                        HTML("<br>DE Analysis in Progress...<br>Please be patient"),
  ), color=spinner$color)

  # extract inputs
  ain <- input$deNormInput
  method <- input$deMethodInput
  comparisons <- input$deComparison
  condition <- input$deColComparison
  DEqMS_column <- input$deDEqMSColumn
  K <- input$deROTSK
  B <- input$deROTSB
  logfc_b <- input$deLogFC
  logfc <- input$deLogFCThr
  padj_b <- input$deMultTest
  p <- input$dePThr
  trend <- input$deLimmaTrend
  robust <- input$deLimmaRobust

  # run DE analysis
  tryCatch({
    de_results <- run_DE(reactiveVals$se,
                         comparisons = comparisons,
                         ain = ain,
                         condition = condition,
                         DE_method = method,
                         logFC = logfc_b,
                         logFC_up = logfc,
                         logFC_down = -logfc,
                         p_adj = padj_b,
                         alpha = p,
                         K = K,
                         B = B,
                         DEqMS_PSMs_column = DEqMS_column,
                         trend = as.logical(trend),
                         robust = as.logical(robust))
    
    # add to steps
    de_analysis_steps <- paste0("DE Analysis: Performed DE analysis using the following parameters: ",
                                                                                           "input data = ", paste0(ain, collapse = "; "),
                                                                                           ", comparisons = ", paste0(comparisons, collapse = "; "),
                                                                                           ", condition = ", condition,
                                                                                           ", method = ", method)
    
    if(logfc_b){
      de_analysis_steps <- paste0(de_analysis_steps, ", logFC threshold = ", logfc)
    }
    if(padj_b){
      de_analysis_steps <- paste0(de_analysis_steps, ", adjusted p-value threshold = ", p)
    } else {
      de_analysis_steps <- paste0(de_analysis_steps, ", p-value threshold = ", p)
    }
    if(method == "ROTS"){
      de_analysis_steps <- paste0(de_analysis_steps, ", K = ", K, ", B = ", B)
    } else if(method == "DEqMS"){
      de_analysis_steps <- paste0(de_analysis_steps, ", DEqMS column = ", DEqMS_column)
    } else if(method == "limma"){
      de_analysis_steps <- paste0(de_analysis_steps, ", trend = ", trend, ", robust = ", robust)
    }
    
    metadata(reactiveVals$se)$steps[[length(metadata(reactiveVals$se)$steps)+1]] <- paste0(de_analysis_steps, ".")
    
    browser()
    
                                          
    # save current settings
    de_current <- data.table("Method" = c(method),
                             "Input Data" = c(ain),
                             "Comparisons" = c(comparisons)
    )
    if(method == "ROTS"){
      de_current$K <- c(K)
      de_current$B <- c(B)
    }
    if(method == "DEqMS"){
      de_current$DEqMS_column <- c(DEqMS_column)
    }
    if(logfc_b){
      de_current[ , ":=" ("LogFC Thr" = logfc)]
    }
    if(padj_b){
      de_current[ , ":=" ("P.Adj Thr" = p)]
    } else {
      de_current[ , ":=" ("P.Value Thr" = p)]
    }
  
    reactiveVals$de_current <- de_current
    reactiveVals$de_results <- de_results
    reactiveVals$de_condition <- condition
  
    updatePickerInput(session = session, inputId = "deVisComparison", choices = comparisons, selected = comparisons)
    updatePickerInput(session = session, inputId = "deVisAin", choices = ain, selected = ain)
  
    # check if gene names column empty
    #rd <- as.data.table(rowData(reactiveVals$se))
    #if(sum(is.na(rd$Gene.Names)) == nrow(rd)){
    #  biomarker_col <- c("Protein.IDs")
    #} else {
    #  biomarker_col <- c("Protein.IDs", "Gene.Names")
    #}
    #updatePickerInput(session = session, inputId =   "deBiomarkerColumn", choices = biomarker_col, selected = biomarker_col[[1]])
    shinyjs::runjs("openBox('deEvaluation')")
    shinyjs::runjs("openBox('deOverviewVis')")
  }, error = function(e){
    shinyalert(
      inputId = "confirmDEAnalysisFailure",
      title = "DE Analysis Failed With the Following Message:",
      text = HTML(sprintf("<b>%s</b>", e$message)),
      size = "s",
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = TRUE,
      type = "error",
      showConfirmButton = FALSE,
      showCancelButton = TRUE,
      cancelButtonText = "Close",
      timer = 5000,
      imageUrl = "",
      animation = TRUE,
      immediate = FALSE
    )
  })
  waiter_hide(id ="app")
})

observeEvent(input$deColComparison, {
  if(!is.null(input$deColComparison)){
    se <- reactiveVals$se
    metadata(se)$condition <- input$deColComparison
    if(!is.null(S4Vectors::metadata(se)$refs)){
      se <- remove_reference_samples(se)
    }
    coldata <- as.data.table(colData(se))
    selected_columns <- input$deColComparison
    # Generate Possible Comparisons
    condition <- unique(coldata[, get(input$deColComparison)])
    comp <- as.data.table(t(combn(condition,2)))
    colnames(comp) <- c("Sample1", "Sample2")
    comp$Comparison <- paste0(comp$Sample1, "-", comp$Sample2)
    updatePickerInput(session = session, inputId = "deComparison", choices = unique(comp$Comparison))
  }
}, ignoreNULL=FALSE)

observeEvent({
  input$deComparison
  input$deNormInput
},{
  if(is.null(input$deComparison) || is.null(input$deNormInput)){
    updateButton(session = session, inputId = "performDEAnalysis", disabled = TRUE)
  } else {
    updateButton(session = session, inputId = "performDEAnalysis", disabled = FALSE)
  }
}, ignoreNULL = FALSE)

observeEvent({
  input$deVisAin
},{
  if(is.null(input$deVisAin)){
    shinyjs::disable("deOverviewBarPlotOptions")
  } else if(length(input$deVisAin) == 1){
    shinyjs::disable("deOverviewBarPlotOptions")
  } else {
    shinyjs::enable("deOverviewBarPlotOptions")

  }
}, ignoreNULL = FALSE)


################################## UI Output ##################################

# Intersection Output
output$de_intersection_plot_tab <- renderUI({
  if(is.null(input$deVisAin)){
    HTML("First Select Parameters and Click on The Button Perform DE Analysis")
  } else if(length(input$deVisAin)==1){
    HTML("Intersection Analysis Not Possible For Only 1 Method!")
  } else {
    fluidRow(
      column(
        1,
        div(
          dropdownButton(
            tags$h3(""),
            numericInput(
              "deIntersectionMinDegree",
              "Minimum Degree of an Intersection: ",
              min = 2,
              max = length(input$deVisAin),
              step = 1,
              value = 2
            ),
            inputId = "deOverviewBarPlotOptions",
            circle = TRUE,
            status = "custom",
            icon = tags$i(class = "fa fa-gear", style = "color: white"),
            width = "400px",
            tooltip = tooltipOptions(title = "Click to see plot options")
          ),
          style = "position:relative; height: 500px;"
        )
      ),
      column(
        11,
        shinycssloaders::withSpinner(plotOutput("de_intersections_plot",  width = "100%", height = "500px")),
        uiOutput("de_intersection_plot_download")
      )
    )
  }
})

output$de_intersection_table_tab <- renderUI({
  if(is.null(input$deVisAin)){
    HTML("First Select Parameters and Click on The Button Perform DE Analysis")
  } else if(length(input$deVisAin)==1){
    HTML("Intersection Analysis Not Possible For Only 1 Method")
  }else {
    fluidRow(
      style = "padding-left:20px;padding-right: 20px",
      shinycssloaders::withSpinner(DT::dataTableOutput("de_intersections_table")),
      uiOutput("de_intersection_table_download")
    )
  }
})


output$de_intersection_consensus_tab <- renderUI({
  if(is.null(input$deVisAin)){
    HTML("First Select Parameters and Click on The Button Perform DE Analysis")
  } else if(length(input$deVisAin)==1){
    HTML("Intersection Analysis Not Possible For Only 1 Method")
  }else {
    fluidRow(
      column(
        3,
        div(
          radioGroupButtons(
            inputId = "deConsensusPerComparison",
            label = "How should the consensus DEPs be reported?",
            choices = c("Per Comparison" = TRUE, "By Ignoring the Comparisons" = FALSE),
            selected = TRUE,
            status = "s"
          ),
          sliderInput(
            inputId = "deConsensusThr",
            label = "Consensus Threshold (Percentage of Normalization Methods that Must Agree on a DEP)",
            min = 10,
            max = 95,
            step = 5,
            value = 70,
            post = "%"
          ),
          style = "position:relative; height: 500px;"
        )
      ),
      column(
        9,
        shinycssloaders::withSpinner(DT::dataTableOutput("de_intersections_consensus_table")),
      )   
    )
  }
})

output$de_intersection_jaccard_tab <- renderUI({
  if(is.null(input$deVisAin)){
    HTML("First Select Parameters and Click on The Button Perform DE Analysis")
  } else if (length(input$deVisAin)==1){
    HTML("Intersection Analysis Not Possible For Only 1 Method")
  } else {
    fluidRow(
      column(
        1,
        div(
          dropdownButton(
            tags$h3("Plot Options"),
            selectizeInput(
              "deJaccardPlotType",
              "Plot Type: ",
              choices = c("Facet by Comparison" = "facet_comp", "Neglecting Comparisons" = "all"),
              multiple = FALSE,
              selected = "facet_comp"
            ),
            inputId = "deOverviewBarPlotOptions",
            circle = TRUE,
            status = "custom",
            icon = tags$i(class = "fa fa-gear", style = "color: white"),
            width = "400px",
            tooltip = tooltipOptions(title = "Click to see plot options")
          ),
          style = "position:relative; height: 500px;"
        )
      ),
      column(
        11,
        shinycssloaders::withSpinner(plotOutput("de_intersections_jaccard_plot",  width = "100%", height = "500px")),
        uiOutput("de_intersection_jaccard_download")
      )
    )
  }
})

output$de_intersections_plot <- renderPlot({
  req(reactiveVals$de_results)
  if(length(input$deVisAin) == 1){
    res <- plot_upset_DE(reactiveVals$de_results, ain = input$deVisAin, comparison = input$deVisComparison, plot_type = "single")[[1]]
  } else {
    nlevels <- length(input$deVisComparison)
    custom_colors <- reactiveVals$selected_palette
    if(nlevels > length(custom_colors)){
      custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
    }
    res <- plot_upset_DE(reactiveVals$de_results, ain = input$deVisAin, comparison = input$deVisComparison, plot_type = "stacked")
    res$upset[[2]] <- res$upset[[2]] + ggplot2::scale_fill_manual(name = "Comparison", values = custom_colors)
  }
  reactiveVals$de_intersection_plot <- res$upset
  reactiveVals$de_intersection_plot
})

output$de_intersections_table <- DT::renderDataTable({
  req(reactiveVals$de_results)
  if(length(input$deVisAin) == 1){
    res <- plot_upset_DE(reactiveVals$de_results, ain = input$deVisAin, comparison = input$deVisComparison, min_degree = input$deIntersectionMinDegree, plot_type = "single")[[1]]
  } else {
    res <- plot_upset_DE(reactiveVals$de_results, ain = input$deVisAin, comparison = input$deVisComparison, min_degree = input$deIntersectionMinDegree, plot_type = "stacked")
  }
  reactiveVals$de_intersection_table <- res$table
  DT::datatable(reactiveVals$de_intersection_table,
                rownames = FALSE,
                options = list(pageLength = 10, searching = TRUE, scrollX = TRUE))
})

output$de_intersections_consensus_table <- DT::renderDataTable({
  req(reactiveVals$de_results)
  reactiveVals$de_consensus_table <- extract_consensus_DE_candidates(reactiveVals$de_results, ain = input$deVisAin, comparisons = input$deVisComparison, per_comparison = input$deConsensusPerComparison, norm_thr = input$deConsensusThr / 100)
  DT::datatable(reactiveVals$de_consensus_table,
                rownames = FALSE,
                options = list(pageLength = 10, searching = TRUE, scrollX = TRUE))
})

output$de_intersections_jaccard_plot <- renderPlot({
  req(reactiveVals$de_results)
  if(length(input$deVisAin) == 1){
    reactiveVals$de_jaccard_plot <- plot_jaccard_heatmap(reactiveVals$de_results, ain = input$deVisAin, comparison = input$deVisComparison, plot_type = "single")[[1]]
  } else {
    reactiveVals$de_jaccard_plot <- plot_jaccard_heatmap(reactiveVals$de_results, ain = input$deVisAin, comparison = input$deVisComparison, plot_type = input$deJaccardPlotType)
  }
  reactiveVals$de_jaccard_plot
})

output$de_intersection_plot_download <- renderUI({
  req(reactiveVals$de_intersection_plot)
  div(
    downloadButton(
      outputId = "DE_Intersection_Plot_Download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download")
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$de_intersection_table_download <- renderUI({
  req(reactiveVals$de_intersection_table)
  div(
    downloadButton(
      outputId = "DE_Intersection_Table_Download",
      label = "Download Table",
      class = "download-butt",
      icon = icon("download")
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$de_intersection_consensus_download <- renderUI({
  req(reactiveVals$de_consensus_table)
  div(
    downloadButton(
      outputId = "DE_Consensus_Table_Download",
      label = "Download Table",
      class = "download-butt",
      icon = icon("download")
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$de_intersection_jaccard_plot_download <- renderUI({
  req(reactiveVals$de_jaccard_plot)
  div(
    downloadButton(
      outputId = "DE_Intersection_Jaccard_Plot_Download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download")
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$DE_Intersection_Table_Download <- downloadHandler(
  filename = function(){
    paste0("DE_Intersection_Table", ".csv")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    write.csv(reactiveVals$de_intersection_table, file, row.names = FALSE)
    waiter_hide(id="app")
  }
)

output$DE_Consensus_Table_Download <- downloadHandler(
  filename = function(){
    paste0("DE_Consensus_Table", ".csv")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    write.csv(reactiveVals$de_consensus_table, file, row.names = FALSE)
    waiter_hide(id="app")
  }
)

output$DE_Intersection_Plot_Download <- downloadHandler(
  filename = function(){
    paste0("DE_Intersection_Plot", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$de_intersection_plot, width=12, height=6)
    waiter_hide(id="app")
  }
)

output$DE_Intersection_Jaccard_Plot_Download <- downloadHandler(
  filename = function(){
    paste0("DE_Intersection_Jaccard_Plot", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$de_jaccard_plot, width=12, height=6)
    waiter_hide(id="app")
  }
)

# DE Analysis Overview
output$de_overview_bar_tab <- renderUI({
  de_results <- reactiveVals$de_results
  de_res <- de_results[(de_results$Assay %in% input$deVisAin & de_results$Comparison %in% input$deVisComparison),]
  de_res <- de_res[de_res$Change != "No Change",]
  if(is.null(input$deVisAin) || is.null(input$deVisComparison)){
    HTML("First Select Parameters and Click on The Button Perform DE Analysis")
  } else if(nrow(de_res) == 0){
    HTML("No DE Results for the Selected Parameters")
  } else {
    fluidRow(
      column(
        1,
        div(
          dropdownButton(
            tags$h3("Plot Options"),
            selectizeInput(
              "deOverviewBarPlotType",
              "Plot Type: ",
              choices = c("Facet by Comparison" = "facet_comp", "Stack by Comparison" = "stacked", "Stack by Comparison with Facet by Regulation" = "facet_regulation"),
              multiple = FALSE,
              selected = "stacked"
            ),
            inputId = "deOverviewBarPlotOptions",
            circle = TRUE,
            status = "custom",
            icon = tags$i(class = "fa fa-gear", style = "color: white"),
            width = "400px",
            tooltip = tooltipOptions(title = "Click to see plot options")
          ),
          style = "position:relative; height: 500px;"
        )
      ),
      column(
        11,
        shinycssloaders::withSpinner(plotOutput("deOverviewBarPlot", width = "100%", height = "500px"))
      )
    )
  }
})

output$deOverviewBarPlot <- renderPlot({
  req(reactiveVals$de_results)
  ain <- input$deVisAin
  comp <- input$deVisComparison
  type <- input$deOverviewBarPlotType
  if(length(ain) == 1){
    reactiveVals$de_overview_bar_plot <- plot_overview_DE_bar(reactiveVals$de_results, ain = ain, comparison = comp, plot_type = "single")[[1]]
  } else {
    de_overview_bar_plot <- plot_overview_DE_bar(reactiveVals$de_results, ain = ain, comparison = comp, plot_type = type)
    if(input$deOverviewBarPlotType != "facet_comp"){
      nlevels <- length(comp)
      custom_colors <- reactiveVals$selected_palette
      if(nlevels > length(custom_colors)){
        custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
      }
      reactiveVals$de_overview_bar_plot <- de_overview_bar_plot + scale_fill_manual(name = "Comparison", values = custom_colors)
    } else {
      reactiveVals$de_overview_bar_plot <- de_overview_bar_plot
    }
  }
  reactiveVals$de_overview_bar_plot
})

output$de_overview_bar_download <- renderUI({
  req(reactiveVals$de_overview_bar_plot)
  div(
    downloadButton(
      outputId = "DE_Overview_Bar_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$DE_Overview_Bar_download <- downloadHandler(
  filename = function(){
    paste0("DE_Overview_Bar", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$de_overview_bar_plot, width=12, height=6)
    waiter_hide(id="app")
  }
)

output$de_overview_tile_tab <- renderUI({
  de_results <- reactiveVals$de_results
  if(is.null(input$deVisAin) || is.null(input$deVisComparison)){
    HTML("First Select Parameters and Click on The Button Perform DE Analysis")
  } else {
    de_res <- de_results[(de_results$Assay %in% input$deVisAin & de_results$Comparison %in% input$deVisComparison),]
    de_res <- de_res[de_res$Change != "No Change",]
    if(nrow(de_res)==0){
      HTML("No DE Results for the Selected Parameters")
    } else {
    shinycssloaders::withSpinner(plotOutput("deOverviewTilePlot", width = "100%", height = "500px"))
    }
  }
})

output$deOverviewTilePlot <- renderPlot({
  req(reactiveVals$de_results)
  reactiveVals$de_overview_tile_plot <- plot_overview_DE_tile(reactiveVals$de_results, ain = input$deVisAin, comparisons = input$deVisComparison)
  reactiveVals$de_overview_tile_plot
})

output$de_overview_tile_download <- renderUI({
  req(reactiveVals$de_overview_tile_plot)
  div(
    downloadButton(
      outputId = "DE_Overview_Tile_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$DE_Overview_Tile_download <- downloadHandler(
  filename = function(){
    paste0("DE_Overview_Tile", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$de_overview_tile_plot, width=12, height=6)
    waiter_hide(id="app")
  }
)

# DE Table Results

output$de_table_results <- renderUI({
  if(is.null(input$deVisAin) || is.null(input$deVisComparison)){
    HTML("First Select Parameters and Click on The Button Perform DE Analysis")
  } else {
    fluidRow(
      style = "padding-left: 20px; padding-right: 20px;",
      radioGroupButtons(
        inputId = "de_table_signif",
        label = "Complete DE Results or Only Significant Changes?",
        choices = c("Complete DE Results", "Only Significant Changes"),
        selected = "Only Significant Changes",
        status = "s"
      ),
      DT::dataTableOutput("deTableResults")
    )
  }
})

output$de_table_results_download <- renderUI({
  req(reactiveVals$de_table)
  div(
    downloadButton(
      outputId = "DE_Results_Table_Download",
      label = "Download Table",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$DE_Results_Table_Download <- downloadHandler(
  filename = function(){
    paste0("DE_Results_Table", ".csv")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    write.csv(reactiveVals$de_table, file, row.names = FALSE)
    waiter_hide(id="app")
  }
)

output$deTableResults <- DT::renderDataTable({
  req(reactiveVals$de_results)
  type <- input$de_table_signif
  ain <- input$deVisAin
  comparison <- input$deVisComparison

  de_results <- reactiveVals$de_results
  de_results <- de_results[de_results$Assay %in% c(ain),]
  de_results <- de_results[de_results$Comparison %in% c(comparison),]

  if(type != "Complete DE Results"){
    de_results <- de_results[de_results$Change != "No Change",]
  }
  reactiveVals$de_table <- de_results
  DT::datatable(as.data.table(de_results),
                rownames = FALSE,
                options = list(pageLength = 10, searching = TRUE, scrollX = TRUE))
})


# Volcano Plots
output$de_volcano_plot_download <- renderUI({
  req(reactiveVals$de_volcano)
  div(
    downloadButton(
      outputId = "DE_Volcano_Plot_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$DE_Volcano_Plot_download <- downloadHandler(
  filename = function(){
    paste0("DE_Volcano_Plot", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$de_volcano, width=12, height=6)
    waiter_hide(id="app")
  }
)

output$de_volcano_tab <- renderUI({
  if(is.null(input$deVisAin) || is.null(input$deVisComparison)){
    HTML("First Select Parameters and Click on The Button Perform DE Analysis")
  } else if(length(input$deVisComparison) == 1 | length(input$deVisAin) == 1){
    shinycssloaders::withSpinner(plotOutput("de_volcano_plot",  width = "100%", height = "700px"))
  } else {
    HTML("This plot is only available when a single comparison or a single normalization method is selected.")
  }
})

output$de_volcano_plot <- renderPlot({
  req(reactiveVals$de_results)
  de_results <- reactiveVals$de_results
  ain <- input$deVisAin
  comparison <- input$deVisComparison
  if(!is.null(ain)){
    if((length(ain)==1) && (length(comparison) != 1)){
      reactiveVals$de_volcano <- plot_volcano_DE(de_results, comparison = comparison, facet_norm = FALSE, facet_comparison = TRUE)[[1]] + ggtitle("")
    } else if(length(comparison) == 1) {
      reactiveVals$de_volcano <- plot_volcano_DE(de_results, comparison = comparison, facet_norm = TRUE)[[1]] + ggtitle("")
    }
    reactiveVals$de_volcano
  } else {
    reactiveVals$de_volcano <- NULL
    return(NULL)
  }
})


# Biomarker Coverage Plot

# output$de_biomarker_coverage_plot_download <- renderUI({
#   req(reactiveVals$deBiomarkerCoverage)
#   div(
#     downloadButton(
#       outputId = "DE_Biomarker_Coverage_Plot_download",
#       label = "Download Plot",
#       class = "download-butt",
#       icon = icon("download"),
#     ),
#     style = "float: right; padding-top: 20px;"
#   )
# })
#
# output$DE_Biomarker_Coverage_Plot_download <- downloadHandler(
#   filename = function(){
#     paste0("DE_Biomarker_Coverage_Plot", ".pdf")
#   },
#   content = function(file){
#     waiter_show(id = "app",html = tagList(spinner$logo,
#                                           HTML("<br>Downloading...")),
#                 color=spinner$color)
#     ggsave(file, plot = reactiveVals$deBiomarkerCoverage, width=12, height=6)
#     waiter_hide(id="app")
#   }
# )
#
# output$de_biomarker_coverage_tab <- renderUI({
#   if(is.null(input$deBiomarkerSelection)){
#     HTML("No Biomarkers Selected!")
#   } else {
#     shinycssloaders::withSpinner(plotOutput("de_biomarker_coverage_plot", width = "100%", height = "500px"))
#   }
# })
#
# output$de_biomarker_coverage_plot <- renderPlot({
#   req(reactiveVals$de_results)
#   de_results <- reactiveVals$de_results
#   ain <- input$deVisAin
#   comparison <- input$deVisComparison
#   id_column <- input$deBiomarkerColumn
#   biomarkers <- input$deBiomarkerSelection
#   if(!is.null(input$deBiomarkerSelection)){
#     de_results <- de_results[de_results$Method %in% c(ain),]
#     de_results <- de_results[de_results$Comparison %in% c(comparison),]
#     reactiveVals$deBiomarkerCoverage <- plot_coverage_biomarkers_DE(reactiveVals$se, de_results, biomarkers, id_column = id_column) + ggtitle("")
#     reactiveVals$deBiomarkerCoverage
#   }
# })
#
#
# # Biomarker Boxplot Plot
# output$de_biomarker_boxplot_download <- renderUI({
#   req(reactiveVals$deBiomarkerBoxplot)
#   div(
#     downloadButton(
#       outputId = "DE_Biomarker_Boxplot_Plot_download",
#       label = "Download Plot",
#       class = "download-butt",
#       icon = icon("download"),
#     ),
#     style = "float: right; padding-top: 20px;"
#   )
# })
#
# output$DE_Biomarker_Boxplot_Plot_download <- downloadHandler(
#   filename = function(){
#     paste0("DE_Biomarker_Boxplot_Plot", ".pdf")
#   },
#   content = function(file){
#     waiter_show(id = "app",html = tagList(spinner$logo,
#                                           HTML("<br>Downloading...")),
#                 color=spinner$color)
#     ggsave(file, plot = reactiveVals$deBiomarkerBoxplot, width=12, height=6)
#     waiter_hide(id="app")
#   }
# )
#
# output$de_biomarker_boxplot_tab <- renderUI({
#   if(is.null(input$deBiomarkerSelection)){
#     HTML("No Biomarkers Selected!")
#   } else {
#     req(reactiveVals$de_results)
#     markers <- input$deBiomarkerSelection
#     fluidRow(
#       column(
#         1,
#         div(
#           dropdownButton(
#             tags$h3("Plot Options"),
#             selectizeInput(
#               "deBiomarkerBoxplotSelection",
#               "Biomarker: ",
#               choices = input$deBiomarkerSelection,
#               multiple = FALSE,
#               selected = input$deBiomarkerSelection[1]
#             ),
#             circle = TRUE,
#             status = "custom",
#             icon = tags$i(class = "fa fa-gear", style = "color: white"),
#             width = "400px",
#             tooltip = tooltipOptions(title = "Click to see plot options")
#           ),
#           style = "position:relative; height: 500px;"
#         )
#       ),
#       column(
#         11,
#         shinycssloaders::withSpinner(plotOutput("de_biomarker_boxplot_plot",  width = "100%", height = "500px")),
#       )
#     )
#   }
# })
#
# output$de_biomarker_boxplot_plot <- renderPlot({
#   req(reactiveVals$de_results)
#   de_results <- reactiveVals$de_results
#   ain <- input$deVisAin
#   comparison <- input$deVisComparison
#   id_column <- input$deBiomarkerColumn
#   biomarker <- input$deBiomarkerBoxplotSelection
#   condition <- reactiveVals$de_condition
#   if(!is.null(input$deBiomarkerSelection)){
#     samplegroups <- unlist(strsplit(comparison, "-"))
#     se <- reactiveVals$se[,reactiveVals$se[[condition]] %in% c(samplegroups)]
#     reactiveVals$deBiomarkerBoxplot <- plot_selection_markers(se, biomarker, id_column = id_column, ain = ain, group_color_by = condition, shape_by = NULL)[[1]] + ggtitle("")
#     reactiveVals$deBiomarkerBoxplot
#   }
# })


# Heatmap

output$de_heatmap_tab <- renderUI({
  de_results <- reactiveVals$de_results
  if(is.null(input$deVisAin) || is.null(input$deVisComparison)){
    HTML("First Select Parameters and Click on The Button Perform DE Analysis")
  } else if(length(input$deVisComparison) == 1 && length(input$deVisAin) == 1){
    de_res <- de_results[(de_results$Assay %in% input$deVisAin & de_results$Comparison %in% input$deVisComparison),]
    de_res <- de_res[de_res$Change != "No Change",]
    if(nrow(de_res)==0){
      HTML("No DE Results for the Selected Parameters")
    } else {
      shinycssloaders::withSpinner(plotOutput("de_heatmap_plot",  width = "100%", height = "750px"))
    }
  } else {
    HTML("This plot is only available when a single comparison and a single normalization method is selected.")
  }
})

output$de_heatmap_plot_download <- renderUI({
  req(reactiveVals$de_heatmap)
  div(
    downloadButton(
      outputId = "DE_Heatmap_Plot_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$DE_Heatmap_Plot_download <- downloadHandler(
  filename = function(){
    paste0("DE_Heatmap_Plot", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    pdf(file, width = 12, height = 8)
    ComplexHeatmap::draw(reactiveVals$de_heatmap)
    dev.off()
    waiter_hide(id="app")
  }
)

output$de_heatmap_plot <- renderPlot({
  req(reactiveVals$de_results)
  de_results <- reactiveVals$de_results
  ain <- input$deVisAin
  comparison <- input$deVisComparison
  condition <- reactiveVals$de_condition
  if(!is.null(ain)){
    if((length(ain)==1) && (length(comparison) == 1)){
      custom_colors <- reactiveVals$selected_palette
      reactiveVals$de_heatmap <- plot_heatmap_DE(se = reactiveVals$se, de_res = de_results, ain = ain, condition = condition, comparison = comparison, col_vector = custom_colors)[[ain]]
      reactiveVals$de_heatmap
    } else {
      reactiveVals$de_heatmap <- NULL
      return(NULL)
    }
  }
})
