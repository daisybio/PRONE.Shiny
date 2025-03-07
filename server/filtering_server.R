# Filtering Server


################################ Helper Functions ############################

remove_proteins_in_column <- function(se, column, value){
  rowdata <- as.data.table(rowData(se))
  rowdata_subset <- rowdata[! rowdata[, get(column)] %in% c(value),]
  se <- se[rowdata_subset$ID, ]
  return(se)
}

################################## Observer ##################################

observeEvent(input$NAheatmap_cluster_samples, {
  if(as.logical(input$NAheatmap_cluster_samples)){
    shinyjs::show("NAheatmap_col_dend")
  } else {
    shinyjs::hide("NAheatmap_col_dend")
  }
})

observeEvent(input$NAheatmap_cluster_proteins, {
  if(as.logical(input$NAheatmap_cluster_proteins)){
    shinyjs::show("NAheatmap_row_dend")
  } else {
    shinyjs::hide("NAheatmap_row_dend")
  }
})
  
  

observeEvent(input$filterProteinsColumn,{
  req(reactiveVals$se)
  col <- input$filterProteinsColumn
  rowdata <- as.data.table(rowData(reactiveVals$se))
  values <- unique(rowdata[[col]])
  values <- values[(!is.na(values)) & (values != "")]
  if(length(values)!=0){
    updatePickerInput(session, "filterProteinsValue", choices = values, clearOptions = TRUE)
  } else {
    updatePickerInput(session, "filterProteinsValue", choices = c("No valid values"), choicesOpt = list(
      disabled = c(TRUE)
    ))
  }
})

observeEvent(input$filterProteinsValue,{
  if (is.null(input$filterProteinsValue)){
    updateButton(session, "filterProteinsColumnButton", disabled = TRUE)
  } else {
    updateButton(session, "filterProteinsColumnButton", disabled = FALSE)
  }
}, ignoreNULL = FALSE)

observeEvent(input$filterProteinsNAThr, {
  req(reactiveVals$curr_NA_thr)
  thr <- input$filterProteinsNAThr
  if(thr != reactiveVals$curr_NA_thr){
    updateButton(session, "filterProteinsNAButton", label = " Filter Proteins", disabled = TRUE)
  } else {
    updateButton(session, "filterProteinsNAButton", label = " Filter Proteins", disabled = FALSE)
  }
})

observeEvent(input$filterProteinsColumnButton,{
  nr_original <- nrow(reactiveVals$se)
  se_removed <- remove_proteins_in_column(reactiveVals$se, input$filterProteinsColumn, input$filterProteinsValue)
  nr_removed <- nrow(se_removed)
  shinyalert(
    inputId = "confirmColumnRemoval",
    title = "Confirm the removal of proteins",
    text = paste0("You are removing ", nr_original - nr_removed, " proteins from your data. This step can not be reverted."),
    size = "s", 
    closeOnEsc = TRUE,
    closeOnClickOutside = FALSE,
    html = FALSE,
    type = "warning",
    showConfirmButton = TRUE,
    showCancelButton = TRUE,
    confirmButtonText = "Remove Proteins",
    confirmButtonCol = "#53A93F",
    cancelButtonText = "Cancel",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  )
})

observeEvent(input$confirmColumnRemoval, {
  if(input$confirmColumnRemoval){
    req(reactiveVals$se)
    metadata(reactiveVals$se)$steps[[length(metadata(reactiveVals$se)$steps)+1]] <- paste0("Filter Proteins: Removed proteins with ", input$filterProteinsValue, " value in column ", input$filterProteinsColumn, ".")
    reactiveVals$se <- filter_out_proteins_by_value(reactiveVals$se, input$filterProteinsColumn, input$filterProteinsValue)
    reactiveVals$se_filtered <- reactiveVals$se
  }
})

observeEvent(input$filterProteinsNAUpdateButton, {
  thr <- input$filterProteinsNAThr * 0.01
  reactiveVals$se_filtered <- filter_out_NA_proteins_by_threshold(reactiveVals$se, thr)
  reactiveVals$curr_NA_thr <- input$filterProteinsNAThr
  updateButton(session, "filterProteinsNAButton", label = " Filter Proteins", disabled = FALSE)
})

observeEvent(input$filterProteinsNAButton, {
  thr <- input$filterProteinsNAThr * 0.01
  se_removed <- filter_out_NA_proteins_by_threshold(reactiveVals$se, thr)
  nr_original <- nrow(reactiveVals$se)
  nr_removed <- nrow(se_removed)
  shinyalert(
    inputId = "confirmNARemoval",
    title = "Confirm the removal of proteins",
    text = paste0("You are removing ", nr_original - nr_removed," proteins from your data. This step can not be reverted."),
    size = "s", 
    closeOnEsc = TRUE,
    closeOnClickOutside = FALSE,
    html = FALSE,
    type = "warning",
    showConfirmButton = TRUE,
    showCancelButton = TRUE,
    confirmButtonText = "Remove Proteins",
    confirmButtonCol = "#53A93F",
    cancelButtonText = "Cancel",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  )
})

observeEvent(input$confirmNARemoval, {
  if(input$confirmNARemoval){
    thr <- input$filterProteinsNAThr * 0.01
    reactiveVals$se <- filter_out_NA_proteins_by_threshold(reactiveVals$se, thr)
    reactiveVals$se_filtered <-  reactiveVals$se
    metadata(reactiveVals$se)$steps[[length(metadata(reactiveVals$se)$steps)+1]] <- paste0("Filter Proteins: Filtered proteins with NA threshold of ", thr, "%.")
  }
})

# TODO: shinyalert if only NA proteins


################################## UI Output ##################################

output$NA_heatmap_plot <- renderPlot({
  req(reactiveVals$se_filtered)
  color_by <- input$NAheatmap_color
  label_by <- input$NAheatmap_label
  cluster_samples <- as.logical(input$NAheatmap_cluster_samples)
  cluster_proteins <- as.logical(input$NAheatmap_cluster_proteins)
  show_row_dend <- as.logical(input$NAheatmap_row_dend)
  show_col_dend <- as.logical(input$NAheatmap_col_dend)
  if(color_by != ""){
    nlevels <- length(unique(colData(reactiveVals$se)[[color_by]]))
    custom_colors <- reactiveVals$selected_palette
    if(nlevels > length(custom_colors)){
      custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
    }
    reactiveVals$NA_heatmap <- plot_NA_heatmap(se = reactiveVals$se_filtered,
                                              color_by = color_by,
                                              label_by = label_by,
                                              cluster_samples = cluster_samples,
                                             cluster_proteins = cluster_proteins,
                                              show_row_dend = show_row_dend,
                                             show_column_dend = show_col_dend,
                                             col_vector = custom_colors)
  }
  reactiveVals$NA_heatmap
})

output$NA_density_plot <- renderPlot({
  req(reactiveVals$se_filtered)
  reactiveVals$NA_density <- plot_NA_density(reactiveVals$se_filtered)
  reactiveVals$NA_density
})

output$NA_frequency_plot <- renderPlot({
  req(reactiveVals$se_filtered)
  reactiveVals$NA_frequency <- plot_NA_frequency(reactiveVals$se_filtered)
  reactiveVals$NA_frequency
})

output$NA_heatmap_download <- renderUI({
  req(reactiveVals$NA_heatmap)
  div(
    downloadButton(
      outputId = "NA_heatmap_handler_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$NA_heatmap_handler_download <- downloadHandler(
  filename = function(){
    paste0("NA_Heatmap", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    pdf(file, width = 12, height = 8)
    ComplexHeatmap::draw(reactiveVals$NA_heatmap)
    dev.off()
    waiter_hide(id="app")
  }
)

output$NA_density_download <- renderUI({
  req(reactiveVals$NA_density)
  div(
    downloadButton(
      outputId = "NA_density_handler_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$NA_density_handler_download <- downloadHandler(
  filename = function(){
    paste0("NA_Density", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(file, plot = reactiveVals$NA_density, width=12, height=6)
    waiter_hide(id="app")
  }
)


output$NA_frequency_download <- renderUI({
  req(reactiveVals$NA_frequency)
  div(
    downloadButton(
      outputId = "NA_frequency_handler_download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download"),
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$NA_frequency_handler_download <- downloadHandler(
  filename = function(){
    paste0("NA_Frequency", ".pdf")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo, 
                                          HTML("<br>Downloading...")), 
                color=spinner$color)
    ggsave(file, plot = reactiveVals$NA_frequency, width=12, height=6)
    waiter_hide(id="app")
  }
)





