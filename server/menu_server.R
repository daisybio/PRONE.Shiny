tab_ids <-
  c("welcome",
    "upload",
    "filterProteins",
    "outlierDetection",
    "normalization",
    "de_analysis",
    "goodbye"
    )

reactiveVals$current_tab <- 1
reactiveVals$max_tab <- 1
reactiveVals$continue <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)

reactiveVals$selected_palette <- RColorBrewer::brewer.pal(8, "Set2")

tabs <- list(
  menuItem("Welcome",
           tabName = tab_ids[1],
           icon = icon("house")),
  menuItem("Upload",
           tabName = tab_ids[2],
           icon = icon("circle-play")),
  menuItem(
    "Filter Proteins",
    tabName = tab_ids[3],
    icon = icon("filter")
  ),
  menuItem(
    "Outlier Detection",
    tabName = tab_ids[4],
    icon = icon("magnifying-glass")
  ),
  menuItem(
    "Normalization",
    tabName = tab_ids[5],
    icon = icon("screwdriver-wrench")
  ),
  menuItem(
    "DE Analysis",
    tabName = tab_ids[6],
    icon = icon("object-group")
  ),
  menuItem(
    "Goodbye",
    tabName = tab_ids[7],
    icon = icon("hands")
  )
)

#observeEvent(reactiveVals$max_tab,{
#  # check NA values
#  if((reactiveVals$max_tab == 3) && (!is.null(reactiveVals$se))){
#    tmp <- remove_complete_NA_proteins(reactiveVals$se)
#    if(nrow(tmp[[2]]) > 0){
#      shinyalert(
#        inputId = "confirmCompleteNAs",
#        title = "Protein Groups with NA values for all samples have been identified",
#        text = paste0(nrow(tmp[[2]])," proteins have only NA values. These proteins should be removed for the analysis."),
#        size = "s",
#        closeOnEsc = TRUE,
#        closeOnClickOutside = TRUE,
#        html = FALSE,
#        type = "warning",
#        showConfirmButton = TRUE,
#        showCancelButton = TRUE,
#        confirmButtonText = "Remove Proteins",
#        confirmButtonCol = "#53A93F",
#        cancelButtonText = "Cancel",
#        timer = 0,
#        imageUrl = "",
#        animation = TRUE,
#        callbackJS =
#
#      )
#    }
#  }
#})

observeEvent({
  reactiveVals$current_tab
  reactiveVals$continue
},
{
  if (reactiveVals$current_tab == 1)
    shinyjs::hide("previousTab")
  else if (reactiveVals$current_tab > 1)
    shinyjs::show("previousTab")
  else
    stop("what is the current id?")
  if (reactiveVals$current_tab < length(tab_ids))
    shinyjs::show("nextTab")
  else if (reactiveVals$current_tab == length(tab_ids))
    shinyjs::hide("nextTab")
  #shinyBS::updateButton(session, inputId = "nextTab", icon = icon("times-circle"), style = "warning", disabled = FALSE)
  else
    stop("what is the current id?")
  if (!reactiveVals$continue[reactiveVals$current_tab] &&
      reactiveVals$current_tab == reactiveVals$max_tab)
    shinyjs::disable("nextTab")
  else if (reactiveVals$continue[reactiveVals$current_tab] || reactiveVals$current_tab != reactiveVals$max_tab) {
    shinyjs::enable("nextTab")
  }
})

prevNames <- c(" Previous", " Welcome", " Upload Data", " Filter Proteins", " Outlier Detection", " Normalization", " DE Analysis")
nextNames <- c(" Upload Data", " Filter Proteins", " Outlier Detection", " Normalization"," DE Analysis", " Goodbye")

observeEvent(input$tabs, {
  reactiveVals$current_tab <- match(input$tabs, tab_ids)
  shinyBS::updateButton(session, inputId = "previousTab", label = prevNames[reactiveVals$current_tab])
  shinyBS::updateButton(session, inputId = "nextTab", label = nextNames[reactiveVals$current_tab])
})

observeEvent(input$previousTab, {
  shinyjs::enable("nextTab")
  reactiveVals$current_tab <- reactiveVals$current_tab - 1
  updateTabItems(session, "tabs", tab_ids[reactiveVals$current_tab]) # here we update in case the new current tab is not a new max tab
  shinyjs::runjs("window.scrollTo(0, 0)")
})

observeEvent(input$nextTab, {
  shinyjs::enable("previousTab")
  reactiveVals$current_tab <- reactiveVals$current_tab + 1
  if (reactiveVals$current_tab > reactiveVals$max_tab){
    reactiveVals$max_tab <- reactiveVals$current_tab
  } else
    updateTabItems(session, "tabs", tab_ids[reactiveVals$current_tab]) # here we update in case the new current tab is not a new max tab
  shinyjs::runjs("window.scrollTo(0, 0)")
})

observeEvent(input$confirmCompleteNAs, {
  if(input$confirmCompleteNAs){
    req(reactiveVals$se)
    tmp <- remove_complete_NA_proteins(reactiveVals$se)
    reactiveVals$se <- tmp[[1]]
    reactiveVals$se_filtered <- tmp[[1]]
  }
})


# grow menu depending on current tab
output$sidebar <- renderUI({
  curr_menu <- sidebarMenu(id = "tabs",
                           tabs[1:reactiveVals$max_tab])

  # here we update in case the new current tab is also a new max tab
  updateTabItems(session, "tabs", tab_ids[reactiveVals$max_tab])

  return(curr_menu)
})

observeEvent(input$download_data, {
  if(is.null(reactiveVals$se)){
    output <- HTML("No Data Uploaded")
  } else if(length(metadata(reactiveVals$se)$steps) == 0){
    methods <- names(assays(reactiveVals$se))
    output <- fluidRow(
      style = "padding-left: 10px; padding-right: 10px",
      awesomeCheckboxGroup(
        inputId = "download_assays",
        label = "Select Data You Want To Download:",
        choices = methods,
        selected = methods
      ),
      radioGroupButtons(
        inputId = "download_assays_type",
        label = "Select Type of Data To Download:",
        choices = c("As SummarizedExperiment", "As CSV"),
        selected = ("As SummarizedExperiment"),
        justified = TRUE
      ),
      downloadButton(
        "downloadData",
        " Download Data ",
        icon = icon("download"),
        style = "success",
        disabled = FALSE,
      ),
      HTML("<br><br>"),
      HTML("<br><b>No Steps Executed with the Current Data !</b>")
    )
  } else {
    methods <- names(assays(reactiveVals$se))
    output <- fluidRow(
        style = "padding-left: 10px; padding-right: 10px",
        awesomeCheckboxGroup(
          inputId = "download_assays",
          label = "Select Data You Want To Download:",
          choices = methods,
          selected = methods
        ),
        radioGroupButtons(
          inputId = "download_assays_type",
          label = "Select Type of Data To Download:",
          choices = c("As SummarizedExperiment", "As CSV"),
          selected = ("As SummarizedExperiment"),
          justified = TRUE
        ),
        downloadButton(
          "downloadData",
          " Download Data ",
          icon = icon("download"),
          style = "success",
          disabled = FALSE,
        ),
        
        # add some space here
        HTML("<br><br>"),
        HTML("<br><b>Executed Data Steps:</b>"),
        tags$ol(
          lapply(metadata(reactiveVals$se)$steps, function(x) {
            tags$li(x)
          })
        ),
        downloadButton(
          "downloadSteps",
          " Download Steps ",
          icon = icon("download"),
          style = "success",
          disabled = FALSE,
        ),
    )
  }
  showModal(
    modalDialog(title = "Download Data",
                output,
                easyClose = TRUE,
                footer = modalButton("Close"))
  )

})

observeEvent(input$download_assays, {
  if(is.null(input$download_assays)){
    # disable button
    disable("downloadData")
    disable("downloadSteps")
  } else if(length(input$download_assays) == 0){
    # disable
    disable("downloadData")
    disable("downloadSteps")
  } else {
    # enable
    enable("downloadData")
    enable("downloadSteps")
    
  }
}, ignoreNULL = FALSE)

output$downloadData <- downloadHandler(
  filename = function(){
    ifelse(input$download_assays_type == "As CSV", "Data_CSV.zip", "Data_SE.rds")
  },
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    if(input$download_assays_type == "As CSV"){
      # Cip Directory
      temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(temp_directory)
      rowdata <- as.data.table(rowData(reactiveVals$se))
      write.csv(as.data.table(colData(reactiveVals$se)), file.path(temp_directory, "experimental_design.csv"))
      for(ain in input$download_assays){
        dt <- as.data.table(assays(reactiveVals$se)[[ain]])
        dt <- cbind(rowdata[, c("Protein.IDs"), with = FALSE], dt)
        file_name <- paste0(ain, "_data.csv")
        write.csv(dt, file.path(temp_directory, file_name))
      }
      zip::zip(
        zipfile = file,
        files = dir(temp_directory),
        root = temp_directory
      )
    } else {
      # SummarizedExperiment
      se <- PRONE::subset_SE_by_norm(reactiveVals$se, input$download_assays)
      saveRDS(se, file)
    }
    waiter_hide(id="app")
  }
)

output$downloadSteps <- downloadHandler(
  filename = "Data_Steps.txt",
  content = function(file){
    waiter_show(id = "app",html = tagList(spinner$logo,
                                          HTML("<br>Downloading...")),
                color=spinner$color)
    steps <- as.character(metadata(reactiveVals$se)$steps)  # Convert to character
    writeLines(steps, file)
    waiter_hide(id="app")
  }
)

output$licenseTable <- DT::renderDataTable({
  deps_table <- sapply(unique(renv::dependencies()$Package), packageDescription, fields="License")
  DT::datatable(data.frame(Package = names(deps_table), License = deps_table),
                rownames = FALSE,
                options = list(pageLength = 10, searching = TRUE))
})

observeEvent(input$dependencies_modal, {
  showModal(
    modalDialog(title = "Dependencies",
                DT::dataTableOutput("licenseTable"),
                easyClose = TRUE,
                footer = modalButton("Close"))
  )
})


observeEvent(input$change_colors, {
  if(is.null(reactiveVals$se)){
    output <- HTML("No Data Uploaded")
  } else {
    methods <- names(assays(reactiveVals$se))
    output <- fluidRow(
      style = "padding-left: 10px; padding-right: 10px",
      HTML("Select the color palette you want to use for plotting. Some of the plots have hard coded colors because they do not depend on the data and colors represent predefined categories."),
      column(
        width = 4,
        radioGroupButtons(
          inputId = "chosenColorPalette",
          label = "Select the Color Palette You Want to Use:",
          choices = c("Colorblind1", "Colorblind2: Wong", "Colorblind3: Tol bright", "Colorblind4: Tol vibrant", "Colorblind5: Tol muted", "RColorBrewer: Set1", "RColorBrewer: Set2", "RColorBrewer: Set3", "RColorBrewer: Pastel1", "RColorBrewer: Pastel2", "RColorBrewer: Paired", "RColorBrewer: Dark2", "RColorBrewer: Accent"),
          selected = "RColorBrewer: Set2",
          justified = FALSE,
          direction = "vertical"
        )
      ),
      column(
        width = 8,
        plotOutput("showColorPalette", width = "100%", height = "500px")
      ),
      column(
        width = 10,
        bsButton(
          "confirmPaletteSelection",
          " Confirm Palette Selection ",
          icon = icon("check"),
          style = "success",
        )
      )
    )
  }
  showModal(
    modalDialog(title = "Select Colors",
                output,
                easyClose = TRUE,
                footer = modalButton("Close"))
  )
})

output$showColorPalette <- renderPlot({
  if(input$chosenColorPalette == "Colorblind1"){
    reactiveVals$chosen_palette <- c("#ff6db6", "#004949", "#db6d00",  "#B2DF8A", 
                                       "#FDB462", "#490092", "#009999", "#8f4e00", 
                                       "#ffdf4d", "#171723","#b66dff")
  }else if(input$chosenColorPalette == "Colorblind2: Wong"){
    reactiveVals$chosen_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  }else if(input$chosenColorPalette == "Colorblind3: Tol bright"){
    reactiveVals$chosen_palette <- c('#4477AA', '#EE6677', '#228833', '#CCBB44', 
                                       '#66CCEE', '#AA3377', '#BBBBBB')
  }else if(input$chosenColorPalette == "Colorblind4: Tol vibrant"){
    reactiveVals$chosen_palette <- c('#EE7733', '#0077BB', '#33BBEE', '#EE3377', 
                                       '#CC3311', '#009988', '#BBBBBB')
  }else if(input$chosenColorPalette == "Colorblind5: Tol muted"){
    reactiveVals$chosen_palette <- c('#CC6677', '#332288', '#DDCC77', '#117733', 
                                       '#88CCEE', '#882255', '#44AA99', '#999933', 
                                       '#AA4499', '#DDDDDD')
  }else if(input$chosenColorPalette  == "RColorBrewer: Set1"){
    reactiveVals$chosen_palette <- RColorBrewer::brewer.pal(9, "Set1")
  }else if(input$chosenColorPalette  == "RColorBrewer: Set2"){
    reactiveVals$chosen_palette <- RColorBrewer::brewer.pal(8, "Set2")
  }else if(input$chosenColorPalette  == "RColorBrewer: Set3"){
    reactiveVals$chosen_palette <- RColorBrewer::brewer.pal(12, "Set3")
  }else if(input$chosenColorPalette  == "RColorBrewer: Pastel1"){
    reactiveVals$chosen_palette <- RColorBrewer::brewer.pal(9, "Pastel1")
  }else if(input$chosenColorPalette  == "RColorBrewer: Pastel2"){
    reactiveVals$chosen_palette <- RColorBrewer::brewer.pal(8, "Pastel2")
  }else if(input$chosenColorPalette  == "RColorBrewer: Paired"){
    reactiveVals$chosen_palette <- RColorBrewer::brewer.pal(12, "Paired")
  }else if(input$chosenColorPalette  == "RColorBrewer: Dark2"){
    reactiveVals$chosen_palette <- RColorBrewer::brewer.pal(8, "Dark2")
  }else if(input$chosenColorPalette  == "RColorBrewer: Accent"){
    reactiveVals$chosen_palette <- RColorBrewer::brewer.pal(8, "Accent")
  }
  return(scales::show_col(reactiveVals$chosen_palette))
})

observeEvent(input$confirmPaletteSelection,{
  reactiveVals$selected_palette <- reactiveVals$chosen_palette
  removeModal()
})
