# Normalization Server


################################ Helper Functions ############################





################################## Observer ##################################

observeEvent(input$normMethods, {
  if(is.null(input$normMethods)){
    updateButton(session = session, inputId = "performNormalization", disabled = TRUE)
  } else {
    updateButton(session = session, inputId = "performNormalization", disabled = FALSE)
  }
}, ignoreNULL = FALSE)

observeEvent(input$performNormalization,{
  waiter_show(id = "app",html = tagList(spinner$logo,
                                        HTML("<br>Normalization in Progress...<br>Please be patient"),
  ), color=spinner$color)
  ain <- input$normComb
  norm_methods <- input$normMethods
  se <- reactiveVals$se
  metadata(se)$condition <- input$normMethodsCondition # for EigenMS normalization

  if(input$normRawLog == "default"){
    on_raw <- NULL
  } else {
    on_raw <- as.logical(input$normRawLog)
  }

  if(!is.null(ain)){
    # generate combination of methods
    norm_methods_all <- as.vector(sapply(ain, function(x) paste0(norm_methods, "_on_", x)))
  } else {
    norm_methods_all <- norm_methods
  }

  # normalize data
  tryCatch({
    # add combination methods
    reactiveVals$se <- normalize_se(se = se,
                                    methods = norm_methods_all,
                                    combination_pattern = "_on_",
                                    on_raw = on_raw,
                                    gamma.0 = input$normMethodsGamma,
                                    reduce_correlation_by = input$normMethodsNormics_reduce_corr,
                                    top_x = input$normMethodsNormics_top_x,
                                    NormicsVSN_quantile = input$normMethodsNormics_quantile,
                                    VSN_quantile = input$normMethodsVSN_quantile
                                    )
    
    # save this as steps
    normalization_step <- paste0("Normalization: ", paste(norm_methods_all, collapse = "; "), " with the following parameters: data input: ", input$normRawLog)
    if("RobNorm" %in% norm_methods_all){
      normalization_step <- paste0(normalization_step, ", RobNorm gamma = ", input$normMethodsGamma)
    }
    if("EigenMS" %in% norm_methods_all){
      normalization_step <- paste0(normalization_step, ", EigenMS condition = ", input$normMethodsCondition)
    }
    if("VSN" %in% norm_methods_all){
      normalization_step <- paste0(normalization_step, ", lts.quantile for VSN = ", input$normMethodsVSN_quantile)
    }
    if("NormicsVSN" %in% norm_methods_all){
      normalization_step <- paste0(normalization_step, 
                                   ", lts.quantile for Normics = ", 
                                   input$normMethodsNormics_quantile, 
                                   ", data reduction for Normics = ", 
                                   input$normMethodsNormics_reduce_corr,
                                   ", top x proteins for Normics = ", 
                                   input$normMethodsNormics_top_x)
    }
    
    metadata(reactiveVals$se)$steps[[length(metadata(reactiveVals$se)$steps)+1]] <- paste0(normalization_step, ".")
    
    # update picker input (normInput) --> for combination of methods
    assays <- names(assays(reactiveVals$se))
    # update picker input (normMethods) --> those that have already been performed --> disable
    all_norm_methods <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                          "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "NormicsVSN", "NormicsMedian", "limBE")

    show_methods <- all_norm_methods[!all_norm_methods %in% names(assays(reactiveVals$se))]
    updatePickerInput(session = session, inputId = "normComb", choices = assays[!assays %in% c("raw", "log2")], selected = NULL)
    updatePickerInput(session = session, inputId = "normMethods", choices = show_methods, selected = NULL)
    updatePickerInput(session = session, inputId = "normVisAin", choices = assays[assays != "raw"], selected = assays[assays != "raw"])
    updatePickerInput(session = session, inputId = "deNormInput", choices = assays[assays != "raw"])
  },
  error = function(e){
    shinyalert(
      inputId = "confirmNormalizationFailure",
      title = "Normalization Failed With the Following Message:",
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
  waiter_hide(id = "app")
})

observeEvent(input$normComb,{
  ain <- input$normComb
  # update picker input (normMethods) --> those that have already been performed --> disable
  norm_methods <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                    "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "NormicsMedian", "NormicsVSN", "limBE")
  choices <- norm_methods[!norm_methods %in% ain]
  batch_choices <- choices[choices %in% c("limBE", "IRS")]
  norm_choices <- choices[!choices %in% c("limBE", "IRS")]
  if(length(batch_choices) == 0){
    new_choices <- list("Normalization" = c(norm_choices))
  } else if (length(norm_choices) == 0){
    new_choices <- list("Batch Effect Correction" = c(batch_choices))
  } else {
    new_choices <- list("Normalization" = c(norm_choices), "Batch Effect Correction" = c(batch_choices))
  }

  updatePickerInput(session = session, inputId = "normMethods", choices = new_choices, selected = NULL)
})


################################## Popover ##################################

normalizationRobNormGammaPopover <- bsPopover( id = "normalizationRobNormGammaInfo",
                                         title = "Numeric representing the exponent of the weighted density..",
                                         content = "When the sample size is small, the fitted population of some proteins could be locally trapped such that the variance of those proteins was very small under a large gamma. To avoid this, a small gamma is recommended. When sample size smaller than 40, then set gamma to 0.5 or 0.1.")

normalizationEigenMSConditionPopover <- bsPopover( id = "normalizationEigenMSConditionInfo",
                                            title = "Condition for EigenMS normalization.",
                                            content = "The column of the metadata specifying the condition for EigenMS normalization.")

normalizationVSNQuantilePopover <- bsPopover( id = "normalizationVSNQuantileInfo",
                                       title = "lts.quantile for VSN normalization.",
                                       content = "Numeric of length 1. The quantile that is used for the resistant least trimmed sum of squares regression (see vsn2 lts.quantile).")

normalizationNormicsVSNQuantilePopover <- bsPopover( id = "normalizationNormicsVSNQuantileInfo",
                                          title = "lts.quantile for NormicsVSN normalization.",
                                          content = "Numeric of length 1. The quantile that is used for the resistant least trimmed sum of squares regression (see vsn2 lts.quantile).")

normalizationNormicsReductionCorrPopover <- bsPopover( id = "normalizationNormicsReductionCorrInfo",
                                            title = "Reduce correlation calculation by only using a subset of the data.",
                                            content = "If the data is too big for the computation of the params, increase this parameter by 2,3,4.... The whole data will still be normalized, but the params are calculated on every second row etc.")

normalizationNormicsTopxPopover <- bsPopover( id = "normalizationNormicsTopxInfo",
                                      title = "Top x proteins to use for actual normalization.",
                                      content = "Number of reference proteins extracted for the calculation of parameters.")

################################## UI Output ##################################


output$normMethodsCombination <- renderUI({
  req(reactiveVals$se)
  available_input <- names(assays(reactiveVals$se))
  available_input <- available_input[! available_input %in% c("raw", "log2")]
  if(length(available_input)==0){
    return(NULL)
} else {
    pickerInput(
      "normComb",
      label = span("Select the input data that you want to perform normalization on.", tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"), id = "normalizationCombInfo"),
      multiple = TRUE,
      choices = NULL
    )
  }
})

# Render Norm Methods Parameter UI
output$normMethodsParameter <- renderUI({
  req(reactiveVals$se)
  cd <- colData(reactiveVals$se)
  cols_cs <- colnames(cd)[colnames(cd) != "Column"]
  selected_cs <- metadata(reactiveVals$se)$batch
  parameters_ui <- list()
  if(!is.null(input$normMethods)){
    if("RobNorm" %in% input$normMethods){
      pGamma <- column(
        width = 4,
        numericInput(
          "normMethodsGamma",
          label = span("Select gamma for RobNorm.", tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"), id = "normalizationRobNormGammaInfo"),
          value = 0.1,
          min = 0.1,
          max = 1,
          step = 0.1
        )
      )
      parameters_ui[["gamma"]] <- pGamma
      parameters_ui[["gamma_popover"]] <- normalizationRobNormGammaPopover
    }
    if("EigenMS" %in% input$normMethods){
      pCondition <- column(
        width = 4,
        pickerInput(
          "normMethodsCondition",
          label = span("Select condition required for EigenMS.", tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"), id = "normalizationEigenMSConditionInfo"),
          choices = cols_cs,
          multiple = FALSE
        )
      )
      parameters_ui[["condition"]] <- pCondition
      parameters_ui[["condition_popover"]] <- normalizationEigenMSConditionPopover
    }
    if("VSN" %in% input$normMethods){
      pVSN_quantile <- column(
        width = 4,
        numericInput(
          "normMethodsVSN_quantile",
          label = span("Select lts.quantile for VSN.", tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"), id = "normalizationVSNQuantileInfo"),
          value = 0.9,
          min = 0.5,
          max = 1.0,
          step = 0.1
        )
      )
      parameters_ui[["VSN_quantile"]] <- pVSN_quantile
      parameters_ui[["VSN_quantile_popover"]] <- normalizationVSNQuantilePopover
    }
    if("NormicsMedian" %in% input$normMethods | "NormicsVSN" %in% input$normMethods){
      pNormics_coefficient <- column(
        width = 4,
        numericInput(
          "normMethodsNormics_reduce_corr",
          label = span("Select data reduction parameter for Normics.", tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"), id = "normalizationNormicsReductionCorrInfo"),
          value = 1,
          min = 1,
          max = 20, # TODO: change -> dynamic based on data set
          step = 1,
        )
      )
      parameters_ui[["Normics_reduction_corr"]] <- pNormics_coefficient
      parameters_ui[["Normics_reduction_corr_popover"]] <- normalizationNormicsReductionCorrPopover
      pNormics_top_x <- column(
        width = 4,
        numericInput(
          "normMethodsNormics_top_x",
          label = span("Select top x proteins for Normics.", tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"), id = "normalizationNormicsTopxInfo"),
          value = 50,
          min = 10,
          max = 200, # TODO: change -> dynamic based on data set
          step = 10,
        )
      )
      parameters_ui[["Normics_top_x"]] <- pNormics_top_x
      parameters_ui[["Normics_top_x_popover"]] <- normalizationNormicsTopxPopover
    }
    if("NormicsVSN" %in% input$normMethods){
      pNormics_quantile <- column(
        width = 4,
        numericInput(
          "normMethodsNormics_quantile",
          label = span("Select lts.quantile for NormicsVSN.", tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"), id = "normalizationNormicsVSNQuantileInfo"),
          value = 0.9,
          min = 0.5,
          max = 1.0,
          step = 0.1
        )
      )
      parameters_ui[["NormicsVSN_quantile"]] <- pNormics_quantile
      parameters_ui[["NormicsVSN_quantile_popover"]] <- normalizationNormicsVSNQuantilePopover

    }
  }
  if(length(parameters_ui) != 0){
    return(parameters_ui)
  }
})

# Boxplot
output$norm_boxplots_tab <- renderUI({
  req(reactiveVals$se)
  ain <- names(assays(reactiveVals$se))
  ain <- ain[ain != c("raw")]
  coldata <- as.data.table(colData(reactiveVals$se))
  color_by <- colnames(coldata)
  # TODO:
  color_by <- color_by[color_by != "Column"]
  selected_color_by <- input$boxplots_color
  fluidRow(
    column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "boxplots_color",
            "Color By: ",
            choices = color_by,
            multiple = FALSE,
            selected = selected_color_by
          ),
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
      shinycssloaders::withSpinner(plotOutput("norm_boxplots_plot", width = "100%", height = "900px")),
      uiOutput("norm_boxplot_download")
    )
  )
})

output$norm_boxplots_plot <- renderPlot({
  req(reactiveVals$se)
  color_by <- input$boxplots_color
  ain <- input$normVisAin
  nlevels <- length(unique(colData(reactiveVals$se)[[color_by]]))
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  reactiveVals$norm_boxplots <- plot_boxplots(reactiveVals$se, ain = ain, color_by = color_by) + scale_fill_manual(name = color_by, values = custom_colors)
  reactiveVals$norm_boxplots
})

output$norm_boxplot_download <- renderUI({
  req(reactiveVals$norm_boxplots)
  div(
    downloadButton(
      outputId = "Norm_Boxplot_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$Norm_Boxplot_download <- downloadHandler(
  filename = function(){
    paste0("Boxplots", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$norm_boxplots, width=12, height=6)
    waiter_hide(id="app")
  }
)

# Densities
output$norm_densities_tab <- renderUI({
  req(reactiveVals$se)
  ain <- names(assays(reactiveVals$se))
  ain <- ain[ain != "raw"]
  coldata <- as.data.table(colData(reactiveVals$se))
  color_by <- colnames(coldata)
  # TODO:
  color_by <- color_by[color_by != "Column"]
  selected_color_by <- input$densities_color
  fluidRow(
    column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "densities_color",
            "Color By: ",
            choices = color_by,
            multiple = FALSE,
            selected = selected_color_by
          ),
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
      shinycssloaders::withSpinner(plotOutput("norm_densities_plot", width = "100%", height = "900px")),
      uiOutput("norm_densities_download")
    )
  )
})

output$norm_densities_plot <- renderPlot({
  req(reactiveVals$se)
  color_by <- input$densities_color
  ain <- input$normVisAin
  nlevels <- length(unique(colData(reactiveVals$se)[[color_by]]))
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  reactiveVals$norm_densities <- plot_densities(reactiveVals$se, ain = ain, color_by = color_by) + scale_color_manual(name = color_by, values = custom_colors)
  reactiveVals$norm_densities
})

output$norm_densities_download <- renderUI({
  req(reactiveVals$norm_densities)
  div(
    downloadButton(
      outputId = "Norm_Densities_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$Norm_Densities_download <- downloadHandler(
  filename = function(){
    paste0("Densities", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$norm_densities, width=12, height=6)
    waiter_hide(id="app")
  }
)

# PCA
output$norm_PCA_tab <- renderUI({
  req(reactiveVals$se)
  # check if only one single normalization method selected -> then no faceting according to method but rather to column
  ain <- input$normVisAin
  ain <- ain[ain != "raw"]
  coldata <- as.data.table(colData(reactiveVals$se))
  color_by <- colnames(coldata)
  # TODO:
  color_by <- color_by[color_by != "Column"]
  selected_color_by <- input$pca_color
  selected_label_by <- input$pca_label
  selected_facet_by <- input$pca_facet
  selected_shape_by <- input$pca_shape

  if(length(ain) == 1){
    selected_facet_by <- input$pca_facet
    global_facet_input <- selectizeInput(
      "pca_facet",
      "Facet By: ",
      choices = c("No Faceting", color_by),
      multiple = FALSE,
      selected = selected_facet_by
    )
  } else {
    global_facet_input <- NULL
  }
  fluidRow(
    column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "pca_color",
            "Color By: ",
            choices = color_by,
            multiple = FALSE,
            selected = selected_color_by
          ),
          selectizeInput(
            "pca_shape",
            "Shape By: ",
            choices = c("No Shaping", color_by),
            multiple = FALSE,
            selected = selected_shape_by
          ),
          selectizeInput(
            "pca_label",
            "Label By: ",
            choices = c("No Labeling",color_by),
            multiple = FALSE,
            selected = selected_label_by
          ),
          global_facet_input,
          selectizeInput(
            "pca_ellipses",
            "Show Ellipses: ",
            choices = c("Yes", "No"),
            multiple = FALSE,
            selected = "No"
          ),
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
      shinycssloaders::withSpinner(plotOutput("norm_PCA_plot", width = "100%", height = "900px")),
      uiOutput("norm_PCA_download")
    )
  )
})

output$norm_PCA_plot <- renderPlot({
  req(reactiveVals$se)
  color_by <- input$pca_color
  facet_by <- input$pca_facet
  label_by <- input$pca_label
  shape_by <- input$pca_shape
  ellipses <- input$pca_ellipses
  ain <- input$normVisAin
  if(ellipses == "Yes"){
    ellipses <- TRUE
  } else {
    ellipses <- FALSE
  }
  if(label_by == "No Labeling"){
    label_by <- NULL
  }
  if(shape_by == "No Shaping"){
    shape_by <- NULL
  }
  if(facet_by == "No Faceting"){
    facet_by <- NULL
  }
  nlevels <- length(unique(colData(reactiveVals$se)[[color_by]]))
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  if(length(ain) != 1){
    reactiveVals$norm_pca <- plot_PCA(reactiveVals$se, ain = ain, label_by = label_by, shape_by = shape_by, color_by = color_by, ellipse = ellipses) + scale_color_manual(name = color_by, values = custom_colors)
  } else {
    reactiveVals$norm_pca <- plot_PCA(reactiveVals$se, ain = ain, label_by = label_by, shape_by = shape_by, color_by = color_by, facet_by = facet_by, ellipse = ellipses, facet_norm = FALSE)[[1]] + scale_color_manual(name = color_by, values = custom_colors)
  }
  reactiveVals$norm_pca
})

output$norm_PCA_download <- renderUI({
  req(reactiveVals$norm_pca)
  div(
    downloadButton(
      outputId = "Norm_PCA_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$Norm_PCA_download <- downloadHandler(
  filename = function(){
    paste0("PCA", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$norm_pca, width=12, height=6)
    waiter_hide(id="app")
  }
)


# PEV
output$norm_pev_tab <- renderUI({
  req(reactiveVals$se)
  fluidRow(
    column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "PEV_diff",
            "Visualize Reduction of PEV Relative to log2:",
            choices = c("Yes" = TRUE, "No" = FALSE),
            multiple = FALSE,
            selected = FALSE
          ),
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
      shinycssloaders::withSpinner(plotOutput("norm_pev_plot", width = "100%", height = "500px")),
      uiOutput("norm_pev_download")
    )
  )
})

output$norm_pev_download <- renderUI({
  req(reactiveVals$PEV_plot_intra)
  div(
    downloadButton(
      outputId = "Norm_PEV_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$norm_pev_plot <- renderPlot({
  req(reactiveVals$se)
  nlevels <- length(input$normVisAin)
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  reactiveVals$PEV_plot_intra <- plot_intragroup_PEV(se = reactiveVals$se, ain = input$normVisAin, condition = input$normVisCondition, diff = as.logical(input$PEV_diff)) + scale_fill_manual(name = "Normalization", values = custom_colors)
  reactiveVals$PEV_plot_intra
})

output$Norm_PEV_download <- downloadHandler(
  filename = function(){
    paste0("PEV", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$PEV_plot_intra, width=12, height=6)
    waiter_hide(id="app")
  }
)

# PCV
output$norm_pcv_tab <- renderUI({
  req(reactiveVals$se)
  fluidRow(
    column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "PCV_diff",
            "Visualize Reduction of PCV Relative to log2:",
            choices = c("Yes" = TRUE, "No" = FALSE),
            multiple = FALSE,
            selected = FALSE
          ),
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
      shinycssloaders::withSpinner(plotOutput("norm_pcv_plot", width = "100%", height = "500px")),
      uiOutput("norm_pcv_download")
    )
  )
})

output$norm_pcv_download <- renderUI({
  req(reactiveVals$PCV_plot_intra)
  div(
    downloadButton(
      outputId = "Norm_PCV_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$norm_pcv_plot <- renderPlot({
  req(reactiveVals$se)
  nlevels <- length(input$normVisAin)
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  reactiveVals$PCV_plot_intra <- plot_intragroup_PCV(reactiveVals$se, ain = input$normVisAin, condition = input$normVisCondition, diff = as.logical(input$PCV_diff)) + scale_fill_manual(name = "Normalization", values = custom_colors)
  reactiveVals$PCV_plot_intra
})

output$Norm_PCV_download <- downloadHandler(
  filename = function(){
    paste0("PCV", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$PCV_plot_intra, width=12, height=6)
    waiter_hide(id="app")
  }
)

# PMAD
output$norm_pmad_tab <- renderUI({
  req(reactiveVals$se)
  fluidRow(
    column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "PMAD_diff",
            "Visualize Reduction of PMAD Relative to log2:",
            choices = c("Yes" = TRUE, "No" = FALSE),
            multiple = FALSE,
            selected = FALSE
          ),
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
      shinycssloaders::withSpinner(plotOutput("norm_pmad_plot", width = "100%", height = "500px")),
      uiOutput("norm_pmad_download")
    )
  )
})

output$norm_pmad_download <- renderUI({
  req(reactiveVals$PMAD_plot_intra)
  div(
    downloadButton(
      outputId = "Norm_PMAD_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$norm_pmad_plot <- renderPlot({
  req(reactiveVals$se)
  nlevels <- length(input$normVisAin)
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  reactiveVals$PMAD_plot_intra <- plot_intragroup_PMAD(se = reactiveVals$se, ain = input$normVisAin, condition = input$normVisCondition, diff = as.logical(input$PMAD_diff)) + scale_fill_manual(name = "Normalization", values = custom_colors)
  reactiveVals$PMAD_plot_intra
})

output$Norm_PMAD_download <- downloadHandler(
  filename = function(){
    paste0("PMAD", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$PMAD_plot_intra, width=12, height=6)
    waiter_hide(id="app")
  }
)

# Intragroup Correlation
output$norm_intr_corr_tab <- renderUI({
  req(reactiveVals$se)
  selected_method <- input$intr_corr_method
  if(is.null(selected_method)){
    selected_method <- "pearson"
  }
  fluidRow(
    column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "intr_corr_method",
            "Correlation Method:",
            choices = c("pearson", "spearman"),
            multiple = FALSE,
            selected = selected_method,
          ),
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
      shinycssloaders::withSpinner(plotOutput("norm_intr_corr_plot", width = "100%", height = "500px")),
      uiOutput("norm_intr_corr_download")
    )
  )
})

output$norm_intr_corr_download <- renderUI({
  req(reactiveVals$intr_corr_plot_intra)
  div(
    downloadButton(
      outputId = "Norm_Intr_Corr_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$norm_intr_corr_plot <- renderPlot({
  req(reactiveVals$se)
  nlevels <- length(input$normVisAin)
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  reactiveVals$intr_corr_plot_intra <- plot_intragroup_correlation(reactiveVals$se, ain = input$normVisAin, condition = input$normVisCondition, method = input$intr_corr_method) + scale_fill_manual(name = "Normalization", values = custom_colors)
  reactiveVals$intr_corr_plot_intra
})

output$Norm_Intr_Corr_download <- downloadHandler(
  filename = function(){
    paste0("Intragroup_Correlation", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    ggsave(file, plot = reactiveVals$intr_corr_plot_intra, width=12, height=6)
    waiter_hide(id="app")
  }
)
