# Outlier Detection Server


################################ Helper Functions ############################


################################## Observer ##################################

observeEvent(input$removeSamplesColumn, {
  updatePickerInput(session = session, inputId = "removeSamples", choices = colData(reactiveVals$se)[[input$removeSamplesColumn]])
})

observeEvent(input$removeSamples,{
  if(length(input$removeSamples) > 0 | !is.null(input$removeSamples)){
    updateButton(session, "removeSamplesManuallyButton", disabled = FALSE)
  } else {
    updateButton(session, "removeSamplesManuallyButton", disabled = TRUE)
  }
}, ignoreNULL = FALSE)

observeEvent(input$performOutlierDetectionButton,{
  waiter_show(id = "app",html = tagList(spinner$logo,
                                        HTML("<br>Outlier Detection in Progress...<br>Please be patient"),
  ), color=spinner$color)
  condition <- input$pomaGroup
  method <- input$pomaMethod
  type <- input$pomaType
  metadata(reactiveVals$se)$condition <- condition
  coeff <- input$pomaCoefficient
  poma_res <- detect_outliers_POMA(reactiveVals$se, ain="log2", method = method, type = type, group = TRUE, coeff = coeff)
  polygon_plot <- poma_res$polygon_plot
  distance_plot <- poma_res$distance_boxplot
  outliers_dt <- as.data.table(poma_res$outliers)
  reactiveVals$poma_type <- type
  reactiveVals$poma_coeff <- coeff
  reactiveVals$poma_method <- method
  reactiveVals$poma_condition <- condition
  reactiveVals$poma_polygon <- polygon_plot
  reactiveVals$poma_boxplot <- distance_plot
  reactiveVals$outliers_dt <- outliers_dt
  # update tab view
  updateTabsetPanel(inputId = "outlier_plots", selected = "outlier_boxplot")
  waiter_hide(id="app")
})

observeEvent(input$removeSamplesManuallyButton, {
  samples_to_remove <- input$removeSamples
  shinyalert(
    inputId = "confirmSamplesManuallyRemoval",
    title = "Confirm the removal of samples",
    text = paste0("You are removing ", length(samples_to_remove), " samples from your data. This step can not be reverted."),
    size = "s", 
    closeOnEsc = TRUE,
    closeOnClickOutside = FALSE,
    html = FALSE,
    type = "warning",
    showConfirmButton = TRUE,
    showCancelButton = TRUE,
    confirmButtonText = "Remove Samples Manually",
    confirmButtonCol = "#53A93F",
    cancelButtonText = "Cancel",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  )
})

observeEvent(input$confirmSamplesManuallyRemoval, {
  if(input$confirmSamplesManuallyRemoval){
    samples_to_remove <- input$removeSamples
    column_rm <- input$removeSamplesColumn
    
    reactiveVals$se <- reactiveVals$se[, ! reactiveVals$se[[column_rm]] %in% samples_to_remove]
    # add removal of sample to steps
    metadata(reactiveVals$se)$steps[[length(metadata(reactiveVals$se)$steps)+1]] <- paste0("Outlier Detection: Manually removed sample(s): ", paste0(samples_to_remove, collapse=", "), ".")
    
    reactiveVals$se_filtered <- reactiveVals$se
    # run outlier detection without samples --> only if outlier detection have already been run
    if(!is.null(reactiveVals$polygon_plot)){
      condition <- input$pomaGroup
      method <- input$pomaMethod
      type <- input$pomaType
      coeff <- input$pomaCoefficient
      metadata(reactiveVals$se)$condition <- condition
      poma_res <- detect_outliers_POMA(reactiveVals$se, ain="log2", method = method, type = type, group = TRUE, coeff = coeff)
      reactiveVals$polygon_plot <- poma_res$polygon_plot
      reactiveVals$distance_plot <- poma_res$distance_boxplot
      reactiveVals$outliers_dt <- as.data.table(poma_res$outliers)
      reactiveVals$poma_type <- type
      reactiveVals$poma_coeff <- coeff
      reactiveVals$poma_method <- method
      reactiveVals$poma_condition <- condition
    }
    # update
    updatePickerInput(session = session, inputId = "removeSamples", choices =   colData(reactiveVals$se)[[input$removeSamplesColumn]])

  }
})

observeEvent(input$removeOutliersButton,{
  dt <- reactiveVals$outliers_dt
  shinyalert(
    inputId = "confirmOutlierRemoval",
    title = "Confirm the removal of outliers",
    text = paste0("You are removing ", nrow(dt), " samples from your data. This step can not be reverted."),
    size = "s", 
    closeOnEsc = TRUE,
    closeOnClickOutside = FALSE,
    html = FALSE,
    type = "warning",
    showConfirmButton = TRUE,
    showCancelButton = TRUE,
    confirmButtonText = "Remove Outlier Samples",
    confirmButtonCol = "#53A93F",
    cancelButtonText = "Cancel",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  )
})

observeEvent(input$confirmOutlierRemoval, {
  if(input$confirmOutlierRemoval){
    reactiveVals$se <- remove_POMA_outliers(reactiveVals$se, reactiveVals$outliers_dt)
    metadata(reactiveVals$se)$steps[[length(metadata(reactiveVals$se)$steps)+1]] <- paste0("Outlier Detection: Removed sample(s) ", 
                                                                                           paste0(reactiveVals$outliers_dt$sample, collapse=", "), 
                                                                                           " with POMA using these parameters: distance method = ", 
                                                                                           reactiveVals$poma_method, ", type = ", 
                                                                                           reactiveVals$poma_type, ", condition = ", 
                                                                                           reactiveVals$poma_condition, ", and coeff = ", 
                                                                                           reactiveVals$poma_coeff, ".")
    reactiveVals$se_filtered <- reactiveVals$se
    # run outlier detection without outliers
    condition <- input$pomaGroup
    method <- input$pomaMethod
    type <- input$pomaType
    metadata(reactiveVals$se)$condition <- condition
    coeff <- input$pomaCoefficient
    poma_res <- detect_outliers_POMA(reactiveVals$se, ain="log2", method = method, type = type, group = TRUE, coeff = coeff)
    reactiveVals$polygon_plot <- poma_res$polygon_plot
    reactiveVals$distance_plot <- poma_res$distance_boxplot
    reactiveVals$outliers_dt <- as.data.table(poma_res$outliers)
    reactiveVals$poma_type <- type
    reactiveVals$poma_coeff <- coeff
    reactiveVals$poma_method <- method
    reactiveVals$poma_condition <- condition
    # update
    updatePickerInput(session = session, inputId = "removeSamples", choices = colData(reactiveVals$se)$Column)
  }
})


################################## UI Output ##################################

# Outlier Samples Table
output$pomaOutlierSamples <- renderUI({
  if(is.null(reactiveVals$outliers_dt)){
    updateButton(session, "removeOutliersButton", disabled = TRUE)
    output_box <- fluidRow(
      column(
        HTML("<b>First Select Parameters and Click on The Button Perform Outlier Detection</b>"),
        width = 12
      )
    )
  } else {
    dt <- reactiveVals$outliers_dt
    if(nrow(dt) == 0){
      updateButton(session, "removeOutliersButton", disabled = TRUE)
      output_box <- fluidRow(
        column(
          HTML("<b>No outliers detected !</b>"),
          width = 12
        )
      )
    } else {
      updateButton(session, "removeOutliersButton", disabled = FALSE)
      output_box <- fluidRow(
        column(
          renderTable(
            dt,
            caption = "Detected Outlier Samples",
            caption.placement = "top"
          ),
          width = 12,
        )
      )
    }
  }
  output_box
})

# POMA Polygon Plot
output$outlier_polygon_plot <- renderPlot({
  condition <- reactiveVals$poma_condition
  nlevels <- length(unique(colData(reactiveVals$se)[[condition]]))
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  reactiveVals$poma_polygon + scale_fill_manual(name = condition, values = custom_colors)
})

# POMA Boxplot Plot
output$outlier_boxplot_plot <- renderPlot({
  condition <- reactiveVals$poma_condition
  nlevels <- length(unique(colData(reactiveVals$se)[[condition]]))
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  reactiveVals$poma_boxplot + scale_fill_manual(name = condition, values = custom_colors)
})

# Total Intensity 
output$outlier_tot_int_tab <- renderUI({
  req(reactiveVals$se)
  coldata <- as.data.table(colData(reactiveVals$se))
  cols <- sapply(coldata, function(x) if(length(unique(x)) == length(x)) return(TRUE) else return(FALSE))
  label_options <- colnames(coldata)[cols]
  # TODO:
  color_by <- colnames(coldata)
  color_by <- color_by[color_by != "Column"]
  selected_color_by <- input$outlier_tot_int_colors
  selected_label_by <- input$outlier_tot_int_labels
  fluidRow(
    column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "outlier_tot_int_colors",
            "Color By: ",
            choices = color_by,
            multiple = FALSE,
            selected = selected_color_by
          ),
          selectizeInput(
            "outlier_tot_int_labels",
            "Label By: ",
            choices = label_options,
            multiple = FALSE,
            selected = selected_label_by
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
      shinycssloaders::withSpinner(plotOutput("outlier_tot_int_plot", width = "100%", height = "750px"))
    )
  )
})

output$outlier_tot_int_plot <- renderPlot({
  req(reactiveVals$se)
  color_by <- input$outlier_tot_int_colors
  label_by <- input$outlier_tot_int_labels
  nlevels <- length(unique(colData(reactiveVals$se)[[color_by]]))
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  reactiveVals$outlier_tot_int_plot <- plot_tot_int_samples(reactiveVals$se, ain = "raw", color_by = color_by, label_by = label_by) + scale_fill_manual(name = color_by, values = custom_colors)
  reactiveVals$outlier_tot_int_plot
})

# Number of Proteins
output$outlier_nr_prot_tab <- renderUI({
  req(reactiveVals$se)
  coldata <- as.data.table(colData(reactiveVals$se))
  cols <- sapply(coldata, function(x) if(length(unique(x)) == length(x)) return(TRUE) else return(FALSE))
  label_options <- colnames(coldata)[cols]
  # TODO:
  color_by <- colnames(coldata)
  color_by <- color_by[color_by != "Column"]
  selected_color_by <- input$outlier_nr_prot_colors
  selected_label_by <- input$outlier_nr_prot_labels
  fluidRow(
    column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "outlier_nr_prot_colors",
            "Color By: ",
            choices = color_by,
            multiple = FALSE,
            selected = selected_color_by
          ),
          selectizeInput(
            "outlier_nr_prot_labels",
            "Label By: ",
            choices = label_options,
            multiple = FALSE,
            selected = selected_label_by
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
      shinycssloaders::withSpinner(plotOutput("outlier_nr_prot_plot", width = "100%", height = "750px"))
    )
  )
})

output$outlier_nr_prot_plot <- renderPlot({
  req(reactiveVals$se)
  color_by <- input$outlier_nr_prot_colors
  label_by <- input$outlier_nr_tot_labels
  nlevels <- length(unique(colData(reactiveVals$se)[[color_by]]))
  custom_colors <- reactiveVals$selected_palette
  if(nlevels > length(custom_colors)){
    custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
  }
  reactiveVals$outlier_nr_prot_plot <- plot_nr_prot_samples(reactiveVals$se, ain = "raw", color_by = color_by, label_by = label_by) + scale_fill_manual(name = color_by, values = custom_colors)
  reactiveVals$outlier_nr_prot_plot
})

output$outlier_poma_polygon_download <- renderUI({
  req(reactiveVals$poma_polygon)
  div(
    downloadButton(
      outputId = "Outlier_POMA_polygon_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$outlier_poma_boxplot_download <- renderUI({
  req(reactiveVals$poma_boxplot)
  div(
    downloadButton(
      outputId = "Outlier_POMA_boxplot_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})



output$outlier_nr_prot_download <- renderUI({
  req(reactiveVals$outlier_nr_prot_plot)
  div(
    downloadButton(
      outputId = "Outlier_nr_prot_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$outlier_tot_int_download <- renderUI({
  req(reactiveVals$outlier_tot_int_plot)
  div(
    downloadButton(
      outputId = "Outlier_tot_int_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})


# Visualization Tabs
output$outlierPlots <- renderUI({
  if(is.null(reactiveVals$poma_boxplot)){
    poma_polygon <- HTML("<b>First Select Parameters and Click on The Button</b>")
    poma_boxplot <- HTML("<b>First Select Parameters and Click on The Button</b>")
  } else {
    poma_polygon <- list(plotOutput("outlier_polygon_plot", width="100%", height="500px"),
                         uiOutput("outlier_poma_polygon_download"))
    poma_boxplot <- list(plotOutput("outlier_boxplot_plot", width="100%", height="500px"),
                         uiOutput("outlier_poma_boxplot_download"))
  }

  tabs <- fluidRow(
    column(
      width = 12,
      tabBox(
        tabPanel(
          uiOutput("outlier_tot_int_tab"),
          div(
            uiOutput("outlier_tot_int_download"),
            style = "float: right; padding-top: 20px;"
          ),
          value = "outlier_tot_int",
          title = "Total Protein Intensity"
        ),
        tabPanel(
          uiOutput("outlier_nr_prot_tab"),
          div(
            uiOutput("outlier_nr_prot_download"),
            style = "float: right; padding-top: 20px;"
          ),
          value = "outlier_nr_prot",
          title = "Number of Proteins Per Sample"
        ),
        tabPanel(
          poma_boxplot,
          value = "outlier_boxplot",
          title = "POMA:Distance to Centroid"
        ),
        tabPanel(
          poma_polygon,
          value = "outlier_polygon",
          title = "POMA:Polygon Plot"
        ),
        id = "outlier_plots",
        title = "Outlier Detection Plots",
        width = 12
      )
    )
    )
  
  tabs
})


output$Outlier_POMA_polygon_download <- downloadHandler(
  filename = function(){
    paste0("Outlier_POMA_Polygon", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(file, plot = reactiveVals$poma_polygon, width=12, height=6)
    waiter_hide(id="app")
  }
)

output$Outlier_POMA_boxplot_download <- downloadHandler(
  filename = function(){
    paste0("Outlier_POMA_Boxplot", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(file, plot = reactiveVals$poma_boxplot, width=12, height=6)
    waiter_hide(id="app")
  }
)

output$Outlier_tot_int_download <- downloadHandler(
  filename = function(){
    paste0("Outlier_Total_Intensity", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(file, plot = reactiveVals$outlier_tot_int_plot, width=12, height=6)
    waiter_hide(id="app")
  }
)

output$Outlier_nr_prot_download <- downloadHandler(
  filename = function(){
    paste0("Outlier_Number_Proteins", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(file, plot = reactiveVals$outlier_nr_prot_plot, width=12, height=6)
    waiter_hide(id="app")
  }
)



