################################## Example Data Info ##################################


exampleDataVector <- c( "Li et al. on hPDLCs (TMT)" = "data/li_et_al.rds",
                        "Vehmas et al. on AROM + mice (LFQ)" = "data/mouse_liver_cytochrome_P450_se.rds"
)

exampleDataDescription <- c("Li et al. on hPDLCs (TMT)" = "Li et al. Dynamic proteomic profiling of human periodontal ligament stem cells during osteogenic differentiation. Stem Cell Research & Therapy 12.1 (Feb. 2021), p. 98. doi: <a href='https://doi.org/10.1186/s13287-020-02123-6'>10.1186/s13287-020-02123-6</a>.
                            </br>
                            </br>
                            Li et al. conducted a TMT experiment to assess the temporal protein expression changes during osteogenic differentiation of human periodontal ligament stem cells (hPDLSCs) at four different time points. The four time points (D0, D3, D7, and D14) with triplicate biological replication were analyzed on three 6-plex TMT sets. The same pooled reference sample was labeled with Label 130 and included in each TMT plex. LC-MS/MS analysis was performed on a Q Exactive mass spectrometer that was coupled to an Easy-nanoLC. MS/MS spectra were searched using Mascot, embedded into Proteome Discoverer.
                            ",
                            "Vehmas et al. on AROM + mice (LFQ)" = "Vehmas et al. Liver lipid metabolism is altered by increased circulating estrogen to androgen ratio in male mouse. J. Proteomics. Bd. 133. S. 66–75. Feb. 2016. doi: <a href='https://doi.org/10.1016/j.jprot.2015.12.009'>10.1016/j.jprot.2015.12.009</a>.
                            </br> 
                            </br> 
                            Vehmas et al. investigated the effect of an increased estrogen-to-androgen ratio on liver lipid metabolism in males. AROM+ transgenic mice were utilized as a model as these mice possess an increased estrogen-to-androgen ratio due to the overexpression of the human P450 aromatose enzyme, which converts androgens to estrogens. The authors studied the impact of this hormonal imbalance on the liver proteome and transcriptome, and plasma phospholipid profile. Liver samples of five AROM+ and seven wild-type male mice were analyzed on an LTQ Obritrap Velos Pro mass spectrometer coupled to an EASY-nanonLC, and the protein data was searched using Proteome Discoverer and the UniProtKB/Swiss-Prot mouse database. ")


################################## Helper Functions ##################################

checkNullTable <- function(toCheck) {
  if (is.null(toCheck))
    return(data.frame("Nothing" = ""))
  else
    return(toCheck)
}

loadData <- function(data, md, metadata_column, protein_column, gene_column, ref_samples_bool, ref_samples, batch_bool, batch_column){
  # check if all samples of metadata are in data column
  # check names
  samples <- make.names(md[,get(metadata_column),])
  md$Column <- samples
  # check if in colnames of data V1
  if("V1" %in% colnames(data)){
    data$V1 <- NULL
  }
  if(!ref_samples_bool){
    ref_samples <- NULL
  }
  if(!batch_bool){
    batch_column <- NULL
  }
  nlevels <- sort(sapply(md, function(x) length(unique(x))))
  md <- md[, names(nlevels), with = FALSE]
  # make columns characters
  for(column in names(md)){
    md[[column]] <- as.character(md[[column]])
  }
  se <- PRONE::load_data(data, md, protein_column = protein_column, gene_column = gene_column, ref_samples = ref_samples, batch_column = batch_column, condition_column = NULL, label_column = NULL)
  return(se)
}

################################## Observer ##################################

observeEvent(input$seFile, {
  if(is.null(input$seFile)){
    updateButton(session, "loadDataSE", label = " Load Data", disabled = TRUE)
  } else if(endsWith(input$seFile$datapath, ".rds")){
    reactiveVals$data$se <- readRDS(input$seFile$datapath)
    updateButton(session, "loadDataSE", label = " Load Data", disabled = FALSE)
  } else {
    updateButton(session, "loadDataSE", label = " Load Data", disabled = TRUE)
  }
}, ignoreNULL = FALSE)

observeEvent(input$metaFile, {
  if(endsWith(input$metaFile$datapath, ".csv") | endsWith(input$metaFile$datapath, ".txt")){
    reactiveVals$data$upload$md <- data.table::fread(input$metaFile$datapath, check.names = TRUE)
    updatePickerInput(session, "columnName", choices = colnames(reactiveVals$data$upload$md), selected = "")
    updatePickerInput(session, "batchColumn", choices = colnames(reactiveVals$data$upload$md), selected = "")
  }
})

observeEvent(input$columnName, {
  md <- reactiveVals$data$upload$md
  column <- input$columnName
  samples <- md[,..column]
  updatePickerInput(session, "refSamples", choices = samples, selected = "")
})

observeEvent(input$dataFile, {
  if(endsWith(input$dataFile$datapath, ".csv") | endsWith(input$dataFile$datapath, "txt")){
    reactiveVals$data$upload$data <- data.table::fread(input$dataFile$datapath, check.names = TRUE)
    updatePickerInput(session, "protIDsColumn", selected = "", choices = colnames(reactiveVals$data$upload$data))
    updatePickerInput(session, "geneNamesColumn", selected = "", choices = c("No Gene Name Column", colnames(reactiveVals$data$upload$data)))
  }
})

observeEvent(input$batchBool,{
  if(input$batchBool){
    shinyjs::show("refSamplesBool")
    shinyjs::show("refSamples")
    shinyjs::show("batchColumn")
  } else {
    updatePickerInput(session, "refSamples", selected = "")
    shinyjs::hide("refSamplesBool")
    shinyjs::hide("refSamples")
    shinyjs::hide("batchColumn")
  }
})

observeEvent(input$refSamplesBool,{
  if(input$refSamplesBool){
    shinyjs::show("refSamples")
  } else {
    shinyjs::hide("refSamples")
    updatePickerInput(session, "refSamples", selected = "")
  }
})

observe({
  md <- !is.null(reactiveVals$data$upload$md)
  data <- !is.null(reactiveVals$data$upload$data)
  column <- !is.null(input$columnName)
  protids <- !is.null(input$protIDsColumn)
  geneNames <- !is.null(input$geneNamesColumn)
  refs <- input$refSamplesBool
  batches <- input$batchBool
  if(batches){
    batches <- !is.null(input$batchColumn)
    if(refs){
      refs <- !is.null(input$refSamples)
    } else {
      refs <- TRUE
    }
  } else {
    batches <- TRUE
    refs <- TRUE
  }
  if(md & data & column & protids & geneNames & refs & batches){
    updateButton(session, "loadDataOwn", disabled = FALSE)
  } else {
    updateButton(session, "loadDataOwn", disabled = TRUE)
  }
})

observeEvent(input$exampleData, {
  if (input$exampleData != ""){
    updateButton(session, "loadDataExample", label = " Load Data", disabled = FALSE)
  } else {
    updateButton(session, "loadDataExample", label = " Load Data", disabled = TRUE)
  }
})

observeEvent(input$loadDataExample, {
  # reset reactiveVals
  ids <- reactiveVal()
  ids(NULL)
  updateButton(session, "loadDataExample", label = " Loading...", disabled = TRUE)
  status <- "warning"
  if(input$exampleData != ""){
    se <- readRDS(input$exampleData)
  }
  reactiveVals$se <- se
  metadata(reactiveVals$se)$exampleData <- input$exampleData
  updatePickerInput(session = session, inputId = "filterProteinsColumn", choices = colnames(rowData(se))) # only show those columns that are not the intensities
  updatePickerInput(session = session, inputId = "pomaGroup", choices = colnames(colData(se))[colnames(colData(se)) != "Column"])
  updatePickerInput(session = session, inputId = "removeSamplesColumn", choices = colnames(colData(se)), selected = colnames(colData(se))[1])
  updatePickerInput(session = session, inputId = "removeSamples", choices = colData(se)[[colnames(colData(se))[1]]])
  updatePickerInput(session = session, inputId = "normInput", choices = names(assays(se))[names(assays(se))!= "raw"])
  updatePickerInput(session = session, inputId = "normVisAin", choices = names(assays(se))[names(assays(se)) != "raw"], selected = names(assays(se))[names(assays(se)) != "raw"])
  updatePickerInput(session = session, inputId = "deNormInput", choices = names(assays(se))[names(assays(se)) != "raw"])
  updatePickerInput(session = session, inputId = "deColComparison", choices = colnames(colData(se))[colnames(colData(se)) != "Column"])
  updateSelectizeInput(session = session, inputId = "NAheatmap_color", choices = colnames(colData(se)), selected = colnames(colData(se))[2])
  all_norm_methods <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                        "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "NormicsVSN", "NormicsMedian", "limBE")
  updatePickerInput(session = session, inputId = "normMethods", choices = all_norm_methods, selected = NULL)
  reactiveVals$de_enrichment_intersection_plot <- NULL
  reactiveVals$de_enrichment_heatmap_plot <- NULL
  reactiveVals$de_enrichment_table <- NULL
  

  # check unique columns (for labeling rows of heatmap)
  cols <- sapply(as.data.table(colData(se)), function(x) if(length(unique(x)) == length(x)) return(TRUE) else return(FALSE))
  label_options <- colnames(colData(se))[cols]
  updateSelectizeInput(session = session, inputId = "NAheatmap_label", choices = label_options, selected = label_options[1])
  # only those columns with at least one sample group with more than 1 value
  choices <- colnames(colData(se))[sapply(as.data.table(colData(se)), function(x) {
    if(length(unique(x)) < length(x)) return(TRUE) else FALSE
  })]
  updatePickerInput(session = session, inputId = "normVisCondition", choices = choices)
  reactiveVals$se_filtered <- se
  reactiveVals$curr_NA_thr <- 0
  reactiveVals$poma_polygon <- NULL
  reactiveVals$poma_boxplot <- NULL
  reactiveVals$outliers_dt <- NULL
  status <- "success"
  if(status == "success"){
    upload_tab <- which(tab_ids == "upload")
    if (isolate(reactiveVals$max_tab) > upload_tab)
      reactiveVals$max_tab <- upload_tab
    reactiveVals$continue[upload_tab] <- TRUE
  } else if (status=="error"){
    shinyalert::shinyalert("The selected column does not match with the column names of the uploaded data!", type ="error")
  }
  updateButton(session, "loadDataExample", label = " Load Data", disabled = FALSE)

  # check if empty gene names --> if then show ProHarMeD alert
  rowdata <- as.data.table(rowData(se))
  empty_gn <- sum(rowdata$Gene.Names == "")
  na_gn <- sum(is.na(rowdata$Gene.Names))
  if((empty_gn != 0) | (na_gn != 0)){
    proharmed_alert <-  shinyalert(
      inputId = "confirmProharmedExample",
      title = "Empty Gene Names Found !",
      text = div(
        HTML("To remap gene names and filter protein IDs, go check "),
        HTML('<img src="ProHarMeD_Logo.png" width="250px" />'),
      ),
      size = "s",
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = TRUE,
      type = "warning",
      showConfirmButton = TRUE,
      showCancelButton = TRUE,
      confirmButtonText = "Visit ProHarMeD",
      confirmButtonCol = "#53A93F",
      cancelButtonText = "Cancel",
      timer = 10000/2,
      imageUrl = "",
      animation = TRUE
    )
  }

})

observeEvent(input$loadDataSE, {
  # reset reactiveVals
  ids <- reactiveVal()
  ids(NULL)
  updateButton(session, "loadDataSE", label = " Loading...", disabled = TRUE)
  status <- "warning"
  if(!is.null(reactiveVals$data$se)){
    se <- reactiveVals$data$se
  }
  
  reactiveVals$se <- se
  reactiveVals$poma_polygon <- NULL
  reactiveVals$poma_boxplot <- NULL
  reactiveVals$outliers_dt <- NULL
  updatePickerInput(session = session, inputId = "filterProteinsColumn", choices = colnames(rowData(se))) # only show those columns that are not the intensities
  updatePickerInput(session = session, inputId = "pomaGroup", choices = colnames(colData(se))[colnames(colData(se)) != "Column"])
  updatePickerInput(session = session, inputId = "removeSamplesColumn", choices = colnames(colData(se)), selected = colnames(colData(se))[1])
  updatePickerInput(session = session, inputId = "removeSamples", choices = colData(se)[[colnames(colData(se))[1]]])
  updatePickerInput(session = session, inputId = "normInput", choices = names(assays(se))[names(assays(se))!= "raw"])
  updatePickerInput(session = session, inputId = "normVisAin", choices = names(assays(se))[names(assays(se)) != "raw"], selected = names(assays(se))[names(assays(se)) != "raw"])
  updatePickerInput(session = session, inputId = "deNormInput", choices = names(assays(se))[names(assays(se)) != "raw"])
  updatePickerInput(session = session, inputId = "deColComparison", choices = colnames(colData(se))[colnames(colData(se)) != "Column"])
  updateSelectizeInput(session = session, inputId = "NAheatmap_color", choices = colnames(colData(se)), selected = colnames(colData(se))[2])
  updatePickerInput(session = session, inputId = "deDEqMSColumn", choices = colnames(rowData(se)))
  all_norm_methods <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                        "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "NormicsVSN", "NormicsMedian", "limBE")
  updatePickerInput(session = session, inputId = "normMethods", choices = all_norm_methods, selected = NULL)
  reactiveVals$de_enrichment_intersection_plot <- NULL
  reactiveVals$de_enrichment_heatmap_plot <- NULL
  reactiveVals$de_enrichment_table <- NULL
  

  # check unique columns (for labeling rows of heatmap)
  cols <- sapply(as.data.table(colData(se)), function(x) if(length(unique(x)) == length(x)) return(TRUE) else return(FALSE))
  label_options <- colnames(colData(se))[cols]
  updateSelectizeInput(session = session, inputId = "NAheatmap_label", choices = label_options, selected = label_options[1])
  updatePickerInput(session = session, inputId = "normVisCondition", choices = colnames(colData(se))[colnames(colData(se)) != "Column"])

  reactiveVals$se_filtered <- se
  reactiveVals$curr_NA_thr <- 0
  status <- "success"
  if(status == "success"){
    upload_tab <- which(tab_ids == "upload")
    if (isolate(reactiveVals$max_tab) > upload_tab)
      reactiveVals$max_tab <- upload_tab
    reactiveVals$continue[upload_tab] <- TRUE
  } else if (status=="error"){
    shinyalert::shinyalert("The selected column does not match with the column names of the uploaded data!", type ="error")
  }
  # check if empty gene names --> if then show ProHarMeD alert
  rowdata <- as.data.table(rowData(se))
  empty_gn <- sum(rowdata$Gene.Names == "")
  na_gn <- sum(is.na(rowdata$Gene.Names))
  if((empty_gn != 0) | (na_gn != 0)){
    shinyalert(
      inputId = "confirmProharmedSE",
      title = "Empty Gene Names Found !",
      text = div(
        HTML("To remap gene names and filter protein IDs, go check "),
        HTML('<img src="ProHarMeD_Logo.png" width="250px" />'),
      ),
      size = "s",
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = TRUE,
      type = "warning",
      showConfirmButton = TRUE,
      showCancelButton = TRUE,
      confirmButtonText = "Visit ProHarMeD",
      confirmButtonCol = "#53A93F",
      cancelButtonText = "Cancel",
      timer = 10000/2,
      imageUrl = "",
      animation = FALSE
    )
  }
  updateButton(session, "loadDataSE", label = " Load Data", disabled = FALSE)
})

observeEvent(input$loadDataOwn, {
  # reset reactiveVals
  ids <- reactiveVal()
  ids(NULL)
  updateButton(session, "loadDataOwn", label = " Loading...", disabled = TRUE)
  status <- "warning"

  tryCatch(
    {
      se <- loadData(reactiveVals$data$upload$data, reactiveVals$data$upload$md, input$columnName, input$protIDsColumn, input$geneNamesColumn, input$refSamplesBool, input$refSamples, input$batchBool, input$batchColumn)
      reactiveVals$se <- se
      reactiveVals$poma_polygon <- NULL
      reactiveVals$poma_boxplot <- NULL
      reactiveVals$outliers_dt <- NULL
      updatePickerInput(session = session, inputId = "filterProteinsColumn", choices = colnames(rowData(se))) # only show those columns that are not the intensities
      updatePickerInput(session = session, inputId = "pomaGroup", choices = colnames(colData(se))[colnames(colData(se)) != "Column"])
      updatePickerInput(session = session, inputId = "removeSamplesColumn", choices = colnames(colData(se)), selected = colnames(colData(se))[1])
      updatePickerInput(session = session, inputId = "removeSamples", choices = colData(se)[[colnames(colData(se))[1]]])
      updatePickerInput(session = session, inputId = "normInput", choices = names(assays(se))[names(assays(se))!= "raw"])
      updatePickerInput(session = session, inputId = "normVisAin", choices = names(assays(se))[names(assays(se)) != "raw"], selected = names(assays(se))[names(assays(se)) != "raw"])
      updatePickerInput(session = session, inputId = "deNormInput", choices = names(assays(se))[names(assays(se)) != "raw"])
      updatePickerInput(session = session, inputId = "deColComparison", choices = colnames(colData(se))[colnames(colData(se)) != "Column"])
      updateSelectizeInput(session = session, inputId = "NAheatmap_color", choices = colnames(colData(se)), selected = colnames(colData(se))[2])
      all_norm_methods <- c("GlobalMean","GlobalMedian", "Median", "Mean", "IRS", "Quantile", "VSN",
                            "LoessF", "LoessCyc", "RLR", "RlrMA", "RlrMACyc", "EigenMS", "MAD", "RobNorm", "TMM", "NormicsVSN", "NormicsMedian", "limBE")
      updatePickerInput(session = session, inputId = "normMethods", choices = all_norm_methods, selected = NULL)
      reactiveVals$de_enrichment_intersection_plot <- NULL
      reactiveVals$de_enrichment_heatmap_plot <- NULL
      reactiveVals$de_enrichment_table <- NULL

      # check unique columns (for labeling rows of heatmap)
      cols <- sapply(as.data.table(colData(se)), function(x) if(length(unique(x)) == length(x)) return(TRUE) else return(FALSE))
      label_options <- colnames(colData(se))[cols]
      updateSelectizeInput(session = session, inputId = "NAheatmap_label", choices = label_options, selected = label_options[1])
      updatePickerInput(session = session, inputId = "normVisCondition", choices = colnames(colData(se))[colnames(colData(se)) != "Column"])

      reactiveVals$se_filtered <- se
      reactiveVals$curr_NA_thr <- 0
      status <- "success"
    },
    error = function(e){
      status <- "error"
      sendSweetAlert(
        session = session,
        title = "Data Loading Failed !",
        text = e,
        type = "error"
      )
    }
  )

  if(status == "success"){
    upload_tab <- which(tab_ids == "upload")
    if (isolate(reactiveVals$max_tab) > upload_tab)
      reactiveVals$max_tab <- upload_tab
    reactiveVals$continue[upload_tab] <- TRUE
  }
  updateButton(session, "loadDataOwn", label = " Load Data", disabled = FALSE)

  # check if empty gene names --> if then show ProHarMeD alert
  rowdata <- as.data.table(rowData(se))
  empty_gn <- sum(rowdata$Gene.Names == "")
  na_gn <- sum(is.na(rowdata$Gene.Names))
  if((empty_gn != 0) | (na_gn != 0)){
    proharmed_alert <-  shinyalert(
      inputId = "confirmProharmedOwn",
      title = "Empty Gene Names Found !",
      text = div(
        HTML("To remap gene names and filter protein IDs, go check "),
        HTML('<img src="ProHarMeD_Logo.png" width="250px" />'),
      ),
      size = "s",
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = TRUE,
      type = "warning",
      showConfirmButton = TRUE,
      showCancelButton = TRUE,
      confirmButtonText = "Visit ProHarMeD",
      confirmButtonCol = "#53A93F",
      cancelButtonText = "Cancel",
      timer = 10000/2,
      imageUrl = "",
      animation = TRUE
    )
  }
})

observeEvent(input$confirmProharmedOwn, {
  if(input$confirmProharmedOwn){
    browseURL("https://proharmed.zbh.uni-hamburg.de")
  }
})

observeEvent(input$confirmProharmedExample, {
  if(input$confirmProharmedExample){
    browseURL("https://proharmed.zbh.uni-hamburg.de")
  }
})

observeEvent(input$confirmProharmedSE, {
  if(input$confirmProharmedSE){
    browseURL("https://proharmed.zbh.uni-hamburg.de")
  }
})


################################## UI Output ##################################


output$metadataTable <- DT::renderDataTable({
  req(reactiveVals$se)
  dt <- DT::datatable(as.data.table(colData(reactiveVals$se)),
                rownames = FALSE,
                options = list(pageLength = 8, searching = TRUE, scrollX = TRUE))
  return(dt)
})

# Output Selected Data Box
output$currentData <- renderInfoBox({
  refSamples <- NULL
  info <- NULL
  se <- NULL
  batchColumn <- NULL
  if(!is.null(reactiveVals$se)){
    se <- reactiveVals$se
    mq <- assays(se)["raw"]
    md <- as.data.table(colData(se))
    refSamples <- metadata(se)$refs
    batchColumn <- metadata(se)$batch
    if(is.null(refSamples)){
      info <- data.frame("Number of Samples" = dim(se)[2],
                         "Number of Proteins" = dim(se)[1]
      )
    } else {
      info <- data.frame("Number of Samples" = dim(se)[2],
                         "Number of Proteins" = dim(se)[1],
                         "Reference Samples" = "Yes"
      )
    }
  }

  status <- "warning"
  # Store Info Table
  overview_current <- NULL
  if(!is.null(info)){
    info_dt <- renderTable(
      checkNullTable(info),
      caption = "Data Overview",
      caption.placement = "top")
    overview_current <- c(list(info_dt))
  }
  # Store NA Overview Table
  if(!is.null(se)){
    NA_overview <- get_NA_overview(se)
    NA_overviewText <- fluidRow(column(
      width = 12,
      renderTable(
        checkNullTable(NA_overview),
        caption = "Overview of Entries",
        caption.placement = "top"
      )
    ))
    overview_current <- c(overview_current, list(NA_overviewText))
  }
  # Store Reference Samples as Table if Reference Samples in Data
  if(!is.null(refSamples)){
    refSamples <- data.frame(Samples = refSamples)
    refSamplesText <- fluidRow(column(
      width = 12,
      renderTable(
        checkNullTable(refSamples),
        caption = "Reference Samples",
        caption.placement = "top"
      )
    ))
    overview_current <- c(overview_current, list(refSamplesText))
  }
  if(!is.null(batchColumn)){
    overview_current <- c(overview_current, list(h5(HTML("<p style=\"color:grey;\">Batch Column</p>")), batchColumn))
  }
  metadatalist <- NULL
  if(is.null(overview_current)){
    status <- "warning"
    collapsed <- TRUE
  } else {
    status <- "success"
    collapsed <- FALSE
    metadatalist <- list(h5(HTML("<p style=\"color:grey;\">Metadata</p>")),
                         DT::dataTableOutput("metadataTable"))
  }
  if((!is.null(reactiveVals$se)) &&  (!is.null(metadata(reactiveVals$se)$exampleData))){
    shinydashboard::box(
      collapsible = TRUE,
      collapsed = collapsed,
      fluidRow(
        column(12, 
               h5(HTML("<p style=\"color:grey;\">Data Description</p>")),
               HTML(exampleDataDescription[names(exampleDataVector[exampleDataVector == metadata(reactiveVals$se)$exampleData])])
        ),
        column(4, list(overview_current)),
        column(7, metadatalist),
        style='padding-top:10px;'
      ),
      title = "Selected Data", status = status)
  } else {
    shinydashboard::box(
      collapsible = TRUE,
      collapsed = collapsed,
      fluidRow(
        column(4, list(overview_current)),
        column(8, metadatalist),
        style='padding-top:10px;'
      ),
      title = "Selected Data", status = status)
  }

})


