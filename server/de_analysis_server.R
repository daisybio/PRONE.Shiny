# DE Analysis Server



DEEnrichmentOrganismPopover <- bsPopover(id = "DEEnrichmentOrganismInfo",
                                         title = "Organism for functional enrichment analysis using gProfiler.",
                                         content = "If you want to execute your enrichment on. If the organism is not yet available, please open a GitHub issue.")

DEEnrichmentIntersectionPopover <-  bsPopover(id = "DEEnrichmentIntersectionPlotInfo",
                                              title = "Jaccard index of the intersection of enrichment terms",
                                              content = "The Jaccard index is a measure of similarity between two sets. It is defined as the size of the intersection divided by the size of the union of the sets.")

DEEnrichmentHeatmapPopover <-  bsPopover(id = "DEEnrichmentHeatmapInfo",
                                         title = "Heatmap storing the terms and their overlap over normalization methods.",
                                         content = "This plot shows the normalization methods on the x-axis with the enrichment terms on the y-axis and is colored according to whether the specific term was significantly enriched for a specific normalization method.")



################################ Helper Functions ############################

# javascript code to collapse box
jscode <- "
shinyjs.collapse = function(boxid) {
$('#' + boxid).closest('.box').find('[data-widget=collapse]').click();
}
"

plot_intersection_enrichment_shiny <- function(se,
                                               de_res,
                                               comparison,
                                               ain,
                                               id_column,
                                               organism,
                                               source,
                                               signif_thr) {
  # Check parameters
  tmp <- check_plot_DE_parameters(de_res, ain, c(comparison))
  de_res <- tmp[[1]]
  ain <- tmp[[2]]
  comparison <- tmp[[3]]
  
  # Check id_column
  if (!id_column %in% colnames(de_res)) {
    # if id_column not in rowData
    stop(paste0(id_column, " not in DE results!"))
  }
  
  # Prepare DE results
  de_res <- de_res[de_res$Change != "No Change", ]
  de_res <- de_res[de_res$Assay %in% ain, ]
  
  de_res <- de_res[de_res$Comparison %in% c(comparison), ]
  
  if (nrow(de_res) == 0) {
    stop("No significant DEPs for selected parameters!")
  }
  
  # check Gene.names
  if(sum(is.na(de_res$Gene.Names))==nrow(de_res)){
    stop("No valid Gene.Names in the dataset!")
  }
  dt <- de_res
  
  assay_ordering <- levels(de_res$Assay)
  if (is.null(assay_ordering)) {
    assay_ordering <- unique(de_res$Assay)
  }
  
  # Prepare queries for gProfiler
  queries <- lapply(ain, function(method) {
    dt <- dt[dt$Assay == method, ]
    query <- dt[[id_column]]
    query <- unique(query[query != ""])
    # split by ;
    query <- unlist(strsplit(query, ";"))
    query <- query[!is.na(query)]
    query
  })
  names(queries) <- ain
  
  # Run gProfiler
  gres <- gprofiler2::gost(
    queries,
    organism = organism,
    sources = c(source),
    user_threshold = signif_thr
  )
  enrichment_table <- gres$result
  
  if (nrow(de_res) == 0) {
    stop("No significantly enriched terms for selected parameters!")
  }
  
  colnames(enrichment_table)[1] <- "Assay"
  
  gres <- gres$result
  gres <- gres[, c("query", "term_name", "source")]
  gres$present <- 1
  dt <- data.table::dcast(data.table::as.data.table(gres), ... ~ query, value.var = "present")
  dt[is.na(dt)] <- 0
  
  
  # Prepare results for clustering and Jaccard
  cluster_dt <- dt
  terms <- cluster_dt$term_name
  cluster_dt$term_name <- NULL
  cluster_dt$source <- NULL
  cluster_dt <- as.data.frame(cluster_dt)
  rownames(cluster_dt) <- terms
  
  # Clustering for heatmap
  tryCatch({
    dist.mat <- vegan::vegdist(cluster_dt, method = "jaccard")
    clust.res <- stats::hclust(dist.mat)
    ordering_terms <- rownames(cluster_dt)[clust.res$order]
    
    dist.mat <- vegan::vegdist(t(cluster_dt), method = "jaccard")
    clust.res <- stats::hclust(dist.mat)
    ordering_assays <- rownames(t(cluster_dt))[clust.res$order]
  }, error = function(e) {
    message("Error in clustering: ", e)
    message("Rows and columns in enrichment heatmap not clustered!")
    ordering_terms <- rownames(cluster_dt)
    ordering_assays <- colnames(cluster_dt)
  })
  
  # Jaccard distance
  jaccard <- function(x, y) {
    M.11 <- sum(x == 1 & y == 1)
    M.10 <- sum(x == 1 & y == 0)
    M.01 <- sum(x == 0 & y == 1)
    return(M.11 / (M.11 + M.10 + M.01))
  }
  
  assays_without_terms <- assay_ordering[!assay_ordering %in% colnames(dt)]
  # add tables with 0 for assays without terms
  for (assay in assays_without_terms) {
    cluster_dt[, assay] <- 0
  }
  
  cluster_dt <- cluster_dt[, assay_ordering]
  
  m <- matrix(
    data = NA,
    nrow = ncol(cluster_dt),
    ncol = ncol(cluster_dt)
  )
  for (i in seq_len(ncol(cluster_dt))) {
    for (j in seq_len(ncol(cluster_dt))) {
      col1 <- colnames(cluster_dt)[i]
      col2 <- colnames(cluster_dt)[j]
      if (col1 == col2) {
        m[i, j] <- 1
      } else if (i > j) {
        m[i, j] <- jaccard(cluster_dt[, col1], cluster_dt[, col2])
      }
    }
  }
  colnames(m) <- colnames(cluster_dt)
  rownames(m) <- colnames(cluster_dt)
  
  melted_m <- reshape2::melt(m, measure.vars = colnames(m), na.rm = TRUE)
  
  # Jaccard heatmap
  melted_m$Var1 <- factor(melted_m$Var1, levels = assay_ordering)
  melted_m$Var2 <- factor(melted_m$Var2, levels = assay_ordering)
  
  jaccard_heatmap <- ggplot2::ggplot(melted_m, ggplot2::aes(x = Var1, y = Var2, fill = value)) +
    ggplot2::geom_tile(color = "black") +
    ggplot2::geom_text(ggplot2::aes(label = round(value, digits = 2), color = value > 0.5)) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#0072B2",
      limits = c(0, 1)
    ) +
    ggplot2::scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"),
                                guide = "none") +
    ggplot2::labs(x = "Normalization Method", y = "Normalization Method", fill = "Jaccard Similarity") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust =
                                                         0.5))
  
  # Prepare data for intersection heatmap
  dt[dt == 1] <- "Yes"
  dt[dt == 0] <- "No"
  
  # Number of unique terms
  nr_unique_terms <- nrow(dt)
  
  melted_dt <- data.table::melt(
    data.table::as.data.table(dt),
    value.name = "Present",
    variable.name = "Assay",
    measure.vars = colnames(dt)[!colnames(dt) %in%
                                  c("term_name", "source")]
  )
  
  melted_dt$term_name <- factor(melted_dt$term_name, levels = ordering_terms)
  melted_dt$Assay <- factor(melted_dt$Assay, levels = ordering_assays)
  enrich_heatmap <- ggplot2::ggplot(melted_dt, ggplot2::aes(
    x = get("Assay"),
    y = get("term_name"),
    fill = get("Present")
  )) + ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(name = "Significant",
                               values = c(No = "grey80", Yes = "#D55E00")) +
    ggplot2::labs(x = "Normalization Method", y = "Terms") +  ggplot2::theme_bw() +
    
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
  
  return(
    list(
      "enrichment_term_plot" = enrich_heatmap,
      "jaccard_intersection_plot" = jaccard_heatmap,
      "enrichment_table" = enrichment_table,
      "nr_unique_terms" = nr_unique_terms
    )
  )
}

################################## Observer ##################################

observeEvent(input$deMultTest, {
  if (input$deMultTest) {
    shinyjs::show(id = "deMultTestMethod")
  } else {
    shinyjs::hide(id = "deMultTestMethod")
  }
})

observeEvent(input$deLogFC, {
  if (input$deLogFC) {
    shinyjs::show(id = "deLogFCThr")
  } else {
    shinyjs::hide(id = "deLogFCThr")
  }
})

observeEvent(input$deMethodInput, {
  if (input$deMethodInput == "ROTS") {
    shinyjs::show(id = "deROTSK")
    shinyjs::show(id = "deROTSB")
  } else {
    shinyjs::hide(id = "deROTSK")
    shinyjs::hide(id = "deROTSB")
  }
  if (input$deMethodInput == "DEqMS") {
    shinyjs::show(id = "deDEqMSColumn")
  } else {
    shinyjs::hide(id = "deDEqMSColumn")
  }
  if (input$deMethodInput == "limma") {
    shinyjs::show(id = "deLimmaTrend")
    shinyjs::show(id = "deLimmaRobust")
    
  } else {
    shinyjs::hide(id = "deLimmaTrend")
    shinyjs::hide(id = "deLimmaRobust")
  }
})

observeEvent(input$performDEAnalysis, {
  waiter_show(
    id = "app",
    html = tagList(
      spinner$logo,
      HTML("<br>DE Analysis in Progress...<br>Please be patient"),
    ),
    color = spinner$color
  )
  
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
    de_results <- run_DE(
      reactiveVals$se,
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
      robust = as.logical(robust)
    )
    
    # add to steps
    de_analysis_steps <- paste0(
      "DE Analysis: Performed DE analysis using the following parameters: ",
      "input data = ",
      paste0(ain, collapse = "; "),
      ", comparisons = ",
      paste0(comparisons, collapse = "; "),
      ", condition = ",
      condition,
      ", method = ",
      method
    )
    
    if (logfc_b) {
      de_analysis_steps <- paste0(de_analysis_steps, ", logFC threshold = ", logfc)
    }
    if (padj_b) {
      de_analysis_steps <- paste0(de_analysis_steps, ", adjusted p-value threshold = ", p)
    } else {
      de_analysis_steps <- paste0(de_analysis_steps, ", p-value threshold = ", p)
    }
    if (method == "ROTS") {
      de_analysis_steps <- paste0(de_analysis_steps, ", K = ", K, ", B = ", B)
    } else if (method == "DEqMS") {
      de_analysis_steps <- paste0(de_analysis_steps, ", DEqMS column = ", DEqMS_column)
    } else if (method == "limma") {
      de_analysis_steps <- paste0(de_analysis_steps, ", trend = ", trend, ", robust = ", robust)
    }
    
    metadata(reactiveVals$se)$steps[[length(metadata(reactiveVals$se)$steps) +
                                       1]] <- paste0(de_analysis_steps, ".")
    
    
    # save current settings
    de_current <- data.table(
      "Method" = c(method),
      "Input Data" = c(ain),
      "Comparisons" = c(comparisons)
    )
    if (method == "ROTS") {
      de_current$K <- c(K)
      de_current$B <- c(B)
    }
    if (method == "DEqMS") {
      de_current$DEqMS_column <- c(DEqMS_column)
    }
    if (logfc_b) {
      de_current[, ":=" ("LogFC Thr" = logfc)]
    }
    if (padj_b) {
      de_current[, ":=" ("P.Adj Thr" = p)]
    } else {
      de_current[, ":=" ("P.Value Thr" = p)]
    }
    
    reactiveVals$de_current <- de_current
    reactiveVals$de_results <- de_results
    reactiveVals$de_condition <- condition
    
    updatePickerInput(
      session = session,
      inputId = "deVisComparison",
      choices = comparisons,
      selected = comparisons
    )
    updatePickerInput(
      session = session,
      inputId = "deVisAin",
      choices = ain,
      selected = ain
    )
    
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
  }, error = function(e) {
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
  waiter_hide(id = "app")
})

observeEvent(input$deColComparison, {
  if (!is.null(input$deColComparison)) {
    se <- reactiveVals$se
    metadata(se)$condition <- input$deColComparison
    if (!is.null(S4Vectors::metadata(se)$refs)) {
      se <- remove_reference_samples(se)
    }
    coldata <- as.data.table(colData(se))
    selected_columns <- input$deColComparison
    # Generate Possible Comparisons
    condition <- unique(coldata[, get(input$deColComparison)])
    comp <- as.data.table(t(combn(condition, 2)))
    colnames(comp) <- c("Sample1", "Sample2")
    comp$Comparison <- paste0(comp$Sample1, "-", comp$Sample2)
    updatePickerInput(
      session = session,
      inputId = "deComparison",
      choices = unique(comp$Comparison)
    )
  }
}, ignoreNULL = FALSE)

observeEvent({
  input$deComparison
  input$deNormInput
}, {
  if (is.null(input$deComparison) || is.null(input$deNormInput)) {
    updateButton(session = session,
                 inputId = "performDEAnalysis",
                 disabled = TRUE)
  } else {
    updateButton(session = session,
                 inputId = "performDEAnalysis",
                 disabled = FALSE)
  }
}, ignoreNULL = FALSE)

observeEvent({
  input$deVisAin
}, {
  if (is.null(input$deVisAin)) {
    shinyjs::disable("deOverviewBarPlotOptions")
  } else if (length(input$deVisAin) == 1) {
    shinyjs::disable("deOverviewBarPlotOptions")
  } else {
    shinyjs::enable("deOverviewBarPlotOptions")
    
  }
}, ignoreNULL = FALSE)

observeEvent({
  input$deEnrichmentComparison
}, {
  req(reactiveVals$de_results)
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% c(input$deEnrichmentComparison)
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  # check if comparison, organism, or source specified
  if (is.null(input$deEnrichmentComparison)){
    updateButton(session = session,
                 inputId = "performEnrichmentAnalysis",
                 disabled = TRUE)
  } else if(nrow(de_res) == 0) {
      updateButton(session = session,
                   inputId = "performEnrichmentAnalysis",
                   disabled = TRUE)
  } else {
    updateButton(session = session,
                 inputId = "performEnrichmentAnalysis",
                 disabled = FALSE)
  }
}, ignoreNULL = FALSE)

observeEvent(input$performEnrichmentAnalysis, {
  waiter_show(
    id = "app",
    html = tagList(
      spinner$logo,
      HTML("<br>Enrichment Analysis in Progress...<br>Please be patient"),
    ),
    color = spinner$color
  )
  
  reactiveVals$de_enrichment_intersection_plot <- NULL
  reactiveVals$de_enrichment_heatmap_plot <- NULL
  reactiveVals$de_enrichment_table <- NULL
  reactiveVals$de_enrichment_nr_unique_terms <- NULL
  
  
  # extract inputs
  ain <- input$deNormInput
  comparison <- input$deEnrichmentComparison
  source <- input$deEnrichmentSource
  id_column <- "Gene.Names"
  organism <- input$deEnrichmentOrganism
  sign_thr <- input$deEnrichmentThr
  de_results <- reactiveVals$de_results
  se <- reactiveVals$se
  # run DE analysis
  tryCatch({
    # run enrichment analysis
    enrich_results <- plot_intersection_enrichment_shiny(
      se = se,
      de_res = de_results,
      comparison = comparison,
      ain = ain,
      id_column = id_column,
      organism = organism,
      source = source,
      signif_thr = sign_thr
    )
    
    reactiveVals$de_enrichment_intersection_plot <- enrich_results$jaccard_intersection_plot
    reactiveVals$de_enrichment_heatmap_plot <- enrich_results$enrichment_term_plot
    reactiveVals$de_enrichment_table <- enrich_results$enrichment_table[, 1:13]
    reactiveVals$de_enrichment_nr_unique_terms <- enrich_results$nr_unique_terms
    
    # add to steps
    de_enrichment_step <- paste0(
      "Enrichment Analysis: Performed functional enrichment using the following parameters: ",
      "input data = ",
      paste0(ain, collapse = "; "),
      ", comparison = ",
      comparison,
      ", source = ",
      source,
      ", id_column = ",
      id_column,
      ", organism = ",
      organism,
      ", significance threshold = ",
      sign_thr
    )
    
    metadata(reactiveVals$se)$steps[[length(metadata(reactiveVals$se)$steps) +
                                       1]] <- paste0(de_enrichment_step, ".")
    
  }, error = function(e) {
    shinyalert(
      inputId = "confirmEnrichmentFailure",
      title = "Enrichment Failed With the Following Message (check your parameters):",
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

################################## UI Output ##################################


# Enrichment Output

output$de_enrichment_tab <- renderUI({
  req(reactiveVals$se)
  div(
    style = "position:relative; height: 800px;",
    fluidRow(
      column(
        width = 3,
        pickerInput(
          "deEnrichmentComparison",
          choices = input$deVisComparison,
          multiple = FALSE,
          options = list("actions-box" = FALSE),
          label = "Select the comparison you want to evaluate."
        ),
      ),
      column(
        width = 3,
        pickerInput(
          "deEnrichmentOrganism",
          choices = c(
            "hsapiens",
            "mmusculus",
            "rnorvegicus",
            "dmelanogaster",
            "scerevisiae",
            "athaliana"
          ),
          multiple = FALSE,
          label = span(
            "Select the organism.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEEnrichmentOrganismInfo"
          )
        ),
        DEEnrichmentOrganismPopover,
      ),
      column(
        width = 3,
        pickerInput(
          "deEnrichmentSource",
          choices = c("GO:BP", "GO:MF", "GO:CC", "KEGG"),
          multiple = FALSE,
          options = list("actions-box" = FALSE),
          label = "Select the database source."
        ),
      ),
      column(
        width = 3,
        numericInput(
          "deEnrichmentThr",
          value = 0.05,
          min = 0.01,
          max = 1,
          step = 0.01,
          label = "Select the significance threshold.",
          
        ),
      ),
    ),
    div(
      bsButton(
        "performEnrichmentAnalysis",
        " Perform Enrichment",
        icon = icon("object-group"),
        style = "success",
        disabled = FALSE,
      ),
      style = "margin-bottom: 30px;"
    ),
    fluidRow(tabBox(
      width = 12,
      tabPanel(
        uiOutput("de_enrichment_intersection_tab"),
        DEEnrichmentIntersectionPopover,
        uiOutput("de_enrichment_intersection_download"),
        # todo
        value = "de_enrichment_1",
        title = span(
          "Intersection Plot",
          tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
          id = "DEEnrichmentIntersectionInfo"
        )
      ),
      tabPanel(
        uiOutput("de_enrichment_heatmap_tab"),
        DEEnrichmentHeatmapPopover,
        uiOutput("de_enrichment_heatmap_download"),
        # todo
        value = "de_intersection_2",
        title = span(
          "Enrichment Heatmap",
          tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
          id = "DEEnrichmentHeatmapInfo"
        )
      ),
      tabPanel(
        uiOutput("de_enrichment_table_tab"),
        uiOutput("de_enrichment_table_download"),
        # todo
        value = "de_intersection_3",
        title = "Enrichment Table"
      ),
    ))
  )
})

output$de_enrichment_intersection_tab <- renderUI({
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% c(input$deEnrichmentComparison)
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if (nrow(de_res) == 0) {
    HTML("<b>No DE Results for the Selected Parameters</b>")
  } else if (is.null(reactiveVals$de_enrichment_intersection_plot)) {
    HTML("<b>First Perform Enrichment Analysis</b>")
  } else {
    fluidRow(column(12, shinycssloaders::withSpinner(
      plotOutput(
        "de_enrichment_intersection_plot",
        width = "100%",
        height = "500px"
      )
    )))
  }
})

output$de_enrichment_intersection_plot <- renderPlot({
  req(reactiveVals$de_enrichment_intersection_plot)
  reactiveVals$de_enrichment_intersection_plot
})

output$de_enrichment_heatmap_tab <- renderUI({
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% c(input$deEnrichmentComparison)
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if (nrow(de_res) == 0) {
    HTML("<b>No DE Results for the Selected Parameters</b>")
  } else if (is.null(reactiveVals$de_enrichment_intersection_plot)) {
    HTML("<b>First Perform Enrichment Analysis</b>")
  } else if (reactiveVals$de_enrichment_nr_unique_terms > 30) {
    HTML(
      "<b>Only Feasible For Less Than 30 Terms. Please Consider Checking the Table</b>"
    )
  } else {
    fluidRow(column(12, shinycssloaders::withSpinner(
      plotOutput(
        "de_enrichment_heatmap_plot",
        width = "100%",
        height = "500px"
      )
    )))
  }
})

output$de_enrichment_heatmap_plot <- renderPlot({
  req(reactiveVals$de_enrichment_heatmap_plot)
  reactiveVals$de_enrichment_heatmap_plot
})

output$de_enrichment_table_tab <- renderUI({
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% c(input$deEnrichmentComparison)
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if (nrow(de_res) == 0) {
    HTML("<b>No DE Results for the Selected Parameters</b>")
  } else if (is.null(reactiveVals$de_enrichment_intersection_plot)) {
    HTML("<b>First Perform Enrichment Analysis</b>")
  }  else {
    fluidRow(column(12, shinycssloaders::withSpinner(
      DT::dataTableOutput("de_enrichment_table")
    ), ))
  }
})

output$de_enrichment_table <- DT::renderDataTable({
  req(reactiveVals$de_enrichment_table)
  DT::datatable(
    reactiveVals$de_enrichment_table,
    rownames = FALSE,
    options = list(
      pageLength = 10,
      searching = TRUE,
      scrollX = TRUE
    )
  )
})

output$de_enrichment_intersection_download <- renderUI({
  req(reactiveVals$de_enrichment_intersection_plot)
  div(
    downloadButton(
      outputId = "DE_Enrichment_Intersection_Download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download")
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$de_enrichment_heatmap_download <- renderUI({
  req(reactiveVals$de_enrichment_heatmap_plot)
  div(
    downloadButton(
      outputId = "DE_Enrichment_Heatmap_Download",
      label = "Download Plot",
      class = "download-butt",
      icon = icon("download")
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$de_enrichment_table_download <- renderUI({
  req(reactiveVals$de_enrichment_table)
  div(
    downloadButton(
      outputId = "DE_Enrichment_Table_Download",
      label = "Download Table",
      class = "download-butt",
      icon = icon("download")
    ),
    style = "float: right; padding-top: 20px;"
  )
})

output$DE_Enrichment_Table_Download <- downloadHandler(
  filename = function() {
    paste0("DE_Enrichment_Table", ".csv")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    write.csv(reactiveVals$de_enrichment_table, file, row.names = FALSE)
    waiter_hide(id = "app")
  }
)

output$DE_Enrichment_Intersection_Download <- downloadHandler(
  filename = function() {
    paste0("DE_Enrichment_Intersection_Plot", ".pdf")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    ggsave(
      file,
      plot = reactiveVals$de_enrichment_intersection_plot,
      width = 12,
      height = 6
    )
    waiter_hide(id = "app")
  }
)

output$DE_Enrichment_Heatmap_Download <- downloadHandler(
  filename = function() {
    paste0("DE_Enrichment_Heatmap_Plot", ".pdf")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    ggsave(
      file,
      plot = reactiveVals$de_enrichment_heatmap_plot,
      width = 12,
      height = 6
    )
    waiter_hide(id = "app")
  }
)


# Intersection Output

output$de_intersection_plot_tab <- renderUI({
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% input$deVisComparison
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if (is.null(input$deVisAin) || is.null(input$deVisComparison)) {
    HTML("<b>First Select Parameters and Click on The Button Perform DE Analysis</b>")
  } else if (length(input$deVisAin) == 1) {
    HTML("<b>Intersection Analysis Not Possible For Only 1 Method!</b>")
  } else if (nrow(de_res) == 0) {
    HTML("<b>No DE Results for the Selected Parameters</b>")
  } else {
    fluidRow(column(
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
      shinycssloaders::withSpinner(
        plotOutput(
          "de_intersections_plot",
          width = "100%",
          height = "500px"
        )
      ),
      uiOutput("de_intersection_plot_download")
    ))
  }
})

output$de_intersection_table_tab <- renderUI({
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% input$deVisComparison
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if (is.null(input$deVisAin) || is.null(input$deVisComparison)) {
    HTML("<b>First Select Parameters and Click on The Button Perform DE Analysis</b>")
  } else if (length(input$deVisAin) == 1) {
    HTML("<b>Intersection Analysis Not Possible For Only 1 Method!</b>")
  } else if (nrow(de_res) == 0) {
    HTML("<b>No DE Results for the Selected Parameters</b>")
  } else {
    fluidRow(
      style = "padding-left:20px;padding-right: 20px",
      shinycssloaders::withSpinner(DT::dataTableOutput("de_intersections_table")),
      uiOutput("de_intersection_table_download")
    )
  }
})

output$de_intersection_consensus_tab <- renderUI({
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% input$deVisComparison
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if (is.null(input$deVisAin) || is.null(input$deVisComparison)) {
    HTML("<b>First Select Parameters and Click on The Button Perform DE Analysis</b>")
  } else if (length(input$deVisAin) == 1) {
    HTML("<b>Intersection Analysis Not Possible For Only 1 Method!</b>")
  } else if (nrow(de_res) == 0) {
    HTML("<b>No DE Results for the Selected Parameters</b>")
  } else {
    fluidRow(column(
      3,
      div(
        radioGroupButtons(
          inputId = "deConsensusPerComparison",
          label = "How should the consensus DEPs be reported?",
          choices = c(
            "Per Comparison" = TRUE,
            "By Ignoring the Comparisons" = FALSE
          ),
          selectfed = TRUE,
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
    column(9, shinycssloaders::withSpinner(
      DT::dataTableOutput("de_intersections_consensus_table")
    ), ))
  }
})

output$de_intersection_jaccard_tab <- renderUI({
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% input$deVisComparison
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if (is.null(input$deVisAin) || is.null(input$deVisComparison)) {
    HTML("</b>First Select Parameters and Click on The Button Perform DE Analysis</b>")
  } else if (length(input$deVisAin) == 1) {
    HTML("</b>Intersection Analysis Not Possible For Only 1 Method!</b>")
  } else if (nrow(de_res) == 0) {
    HTML("</b>No DE Results for the Selected Parameters</b>")
  } else {
    fluidRow(column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "deJaccardPlotType",
            "Plot Type: ",
            choices = c(
              "Facet by Comparison" = "facet_comp",
              "Neglecting Comparisons" = "all"
            ),
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
      shinycssloaders::withSpinner(
        plotOutput(
          "de_intersections_jaccard_plot",
          width = "100%",
          height = "500px"
        )
      ),
      uiOutput("de_intersection_jaccard_download")
    ))
  }
})

output$de_intersections_plot <- renderPlot({
  req(reactiveVals$de_results)
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% input$deVisComparison
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if(nrow(de_res) == 0) {
    return(NULL)
  }
  if (length(input$deComparison) == 1) {
    res <- plot_upset_DE(
      reactiveVals$de_results,
      ain = input$deVisAin,
      comparison = input$deVisComparison,
      plot_type = "single"
    )[[1]]
  } else {
    nlevels <- length(input$deVisComparison)
    custom_colors <- reactiveVals$selected_palette
    if (nlevels > length(custom_colors)) {
      custom_colors <- grDevices::colorRampPalette(colors = reactiveVals$selected_palette)(nlevels)
    }
    res <- plot_upset_DE(
      reactiveVals$de_results,
      ain = input$deVisAin,
      comparison = input$deVisComparison,
      plot_type = "stacked"
    )
    res$upset[[2]] <- res$upset[[2]] + ggplot2::scale_fill_manual(name = "Comparison", values = custom_colors)
  }
  reactiveVals$de_intersection_plot <- res$upset
  reactiveVals$de_intersection_plot
})

output$de_intersections_table <- DT::renderDataTable({
  req(reactiveVals$de_results)
  if (length(input$deVisAin) == 1) {
    res <- plot_upset_DE(
      reactiveVals$de_results,
      ain = input$deVisAin,
      comparison = input$deVisComparison,
      min_degree = input$deIntersectionMinDegree,
      plot_type = "single"
    )[[1]]
  } else {
    res <- plot_upset_DE(
      reactiveVals$de_results,
      ain = input$deVisAin,
      comparison = input$deVisComparison,
      min_degree = input$deIntersectionMinDegree,
      plot_type = "stacked"
    )
  }
  reactiveVals$de_intersection_table <- res$table
  DT::datatable(
    reactiveVals$de_intersection_table,
    rownames = FALSE,
    options = list(
      pageLength = 10,
      searching = TRUE,
      scrollX = TRUE
    )
  )
})

output$de_intersections_consensus_table <- DT::renderDataTable({
  req(reactiveVals$de_results)
  reactiveVals$de_consensus_table <- extract_consensus_DE_candidates(
    reactiveVals$de_results,
    ain = input$deVisAin,
    comparisons = input$deVisComparison,
    per_comparison = input$deConsensusPerComparison,
    norm_thr = input$deConsensusThr / 100
  )
  DT::datatable(
    reactiveVals$de_consensus_table,
    rownames = FALSE,
    options = list(
      pageLength = 10,
      searching = TRUE,
      scrollX = TRUE
    )
  )
})

output$de_intersections_jaccard_plot <- renderPlot({
  req(reactiveVals$de_results)
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% input$deVisComparison
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if(nrow(de_res) == 0) {
    return(NULL)
  }
  if (length(input$deVisAin) == 1) {
    reactiveVals$de_jaccard_plot <- plot_jaccard_heatmap(
      reactiveVals$de_results,
      ain = input$deVisAin,
      comparison = input$deVisComparison,
      plot_type = "single"
    )[[1]]
  } else {
    reactiveVals$de_jaccard_plot <- plot_jaccard_heatmap(
      reactiveVals$de_results,
      ain = input$deVisAin,
      comparison = input$deVisComparison,
      plot_type = input$deJaccardPlotType
    )
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
  filename = function() {
    paste0("DE_Intersection_Table", ".csv")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    write.csv(reactiveVals$de_intersection_table, file, row.names = FALSE)
    waiter_hide(id = "app")
  }
)

output$DE_Consensus_Table_Download <- downloadHandler(
  filename = function() {
    paste0("DE_Consensus_Table", ".csv")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    write.csv(reactiveVals$de_consensus_table, file, row.names = FALSE)
    waiter_hide(id = "app")
  }
)

output$DE_Intersection_Plot_Download <- downloadHandler(
  filename = function() {
    paste0("DE_Intersection_Plot", ".pdf")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    ggsave(
      file,
      plot = reactiveVals$de_intersection_plot,
      width = 12,
      height = 6
    )
    waiter_hide(id = "app")
  }
)

output$DE_Intersection_Jaccard_Plot_Download <- downloadHandler(
  filename = function() {
    paste0("DE_Intersection_Jaccard_Plot", ".pdf")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    ggsave(
      file,
      plot = reactiveVals$de_jaccard_plot,
      width = 12,
      height = 6
    )
    waiter_hide(id = "app")
  }
)

# DE Analysis Overview
output$de_overview_bar_tab <- renderUI({
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% input$deVisComparison
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if (is.null(input$deVisAin) ||
      is.null(input$deVisComparison)) {
    HTML("<b>First Select Parameters and Click on The Button Perform DE Analysis</b>")
  } else if (nrow(de_res) == 0) {
    HTML("<b>No DE Results for the Selected Parameters</b>")
  } else {
    fluidRow(column(
      1,
      div(
        dropdownButton(
          tags$h3("Plot Options"),
          selectizeInput(
            "deOverviewBarPlotType",
            "Plot Type: ",
            choices = c(
              "Facet by Comparison" = "facet_comp",
              "Stack by Comparison" = "stacked",
              "Stack by Comparison with Facet by Regulation" = "facet_regulation"
            ),
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
    column(11, shinycssloaders::withSpinner(
      plotOutput("deOverviewBarPlot", width = "100%", height = "500px")
    )))
  }
})

output$deOverviewBarPlot <- renderPlot({
  req(reactiveVals$de_results)
  ain <- input$deVisAin
  comp <- input$deVisComparison
  type <- input$deOverviewBarPlotType
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% input$deVisComparison
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if(nrow(de_res) == 0) {
    return(NULL)
  }
  if (length(ain) == 1) {
    reactiveVals$de_overview_bar_plot <- plot_overview_DE_bar(
      reactiveVals$de_results,
      ain = ain,
      comparison = comp,
      plot_type = "single"
    )[[1]]
  } else {
    de_overview_bar_plot <- plot_overview_DE_bar(
      reactiveVals$de_results,
      ain = ain,
      comparison = comp,
      plot_type = type
    )
    if (input$deOverviewBarPlotType != "facet_comp") {
      nlevels <- length(comp)
      custom_colors <- reactiveVals$selected_palette
      if (nlevels > length(custom_colors)) {
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
  filename = function() {
    paste0("DE_Overview_Bar", ".pdf")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    ggsave(
      file,
      plot = reactiveVals$de_overview_bar_plot,
      width = 12,
      height = 6
    )
    waiter_hide(id = "app")
  }
)

output$de_overview_tile_tab <- renderUI({
  de_results <- reactiveVals$de_results
  if (is.null(input$deVisAin) ||
      is.null(input$deVisComparison)) {
    HTML("<b>First Select Parameters and Click on The Button Perform DE Analysis</b>")
  } else {
    de_res <- de_results[(
      de_results$Assay %in% input$deVisAin &
        de_results$Comparison %in% input$deVisComparison
    ), ]
    de_res <- de_res[de_res$Change != "No Change", ]
    if (nrow(de_res) == 0) {
      HTML("<b>No DE Results for the Selected Parameters</b>")
    } else {
      shinycssloaders::withSpinner(plotOutput(
        "deOverviewTilePlot",
        width = "100%",
        height = "500px"
      ))
    }
  }
})

output$deOverviewTilePlot <- renderPlot({
  req(reactiveVals$de_results)
  de_results <- reactiveVals$de_results
  de_res <- de_results[(
    de_results$Assay %in% input$deVisAin &
      de_results$Comparison %in% input$deVisComparison
  ), ]
  de_res <- de_res[de_res$Change != "No Change", ]
  if(nrow(de_res) == 0) {
    return(NULL)
  }
  reactiveVals$de_overview_tile_plot <- plot_overview_DE_tile(
    reactiveVals$de_results,
    ain = input$deVisAin,
    comparisons = input$deVisComparison
  )
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
  filename = function() {
    paste0("DE_Overview_Tile", ".pdf")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    ggsave(
      file,
      plot = reactiveVals$de_overview_tile_plot,
      width = 12,
      height = 6
    )
    waiter_hide(id = "app")
  }
)

# DE Table Results

output$de_table_results <- renderUI({
  if (is.null(input$deVisAin) || is.null(input$deVisComparison)) {
    HTML("<b>First Select Parameters and Click on The Button Perform DE Analysis</b>")
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
  filename = function() {
    paste0("DE_Results_Table", ".csv")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    write.csv(reactiveVals$de_table, file, row.names = FALSE)
    waiter_hide(id = "app")
  }
)

output$deTableResults <- DT::renderDataTable({
  req(reactiveVals$de_results)
  type <- input$de_table_signif
  ain <- input$deVisAin
  comparison <- input$deVisComparison
  
  de_results <- reactiveVals$de_results
  de_results <- de_results[de_results$Assay %in% c(ain), ]
  de_results <- de_results[de_results$Comparison %in% c(comparison), ]
  
  if (type != "Complete DE Results") {
    de_results <- de_results[de_results$Change != "No Change", ]
  }
  reactiveVals$de_table <- de_results
  DT::datatable(
    as.data.table(de_results),
    rownames = FALSE,
    options = list(
      pageLength = 10,
      searching = TRUE,
      scrollX = TRUE
    )
  )
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
  filename = function() {
    paste0("DE_Volcano_Plot", ".pdf")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    ggsave(file,
           plot = reactiveVals$de_volcano,
           width = 12,
           height = 6)
    waiter_hide(id = "app")
  }
)

output$de_volcano_tab <- renderUI({
  if (is.null(input$deVisAin) || is.null(input$deVisComparison)) {
    HTML("<b>First Select Parameters and Click on The Button Perform DE Analysis</b>")
  } else if (length(input$deVisComparison) == 1 |
             length(input$deVisAin) == 1) {
    shinycssloaders::withSpinner(plotOutput("de_volcano_plot", width = "100%", height = "700px"))
  } else {
    HTML("<b>Only Available For Single Comparison or Single Normalization Method</b>")
  }
})

output$de_volcano_plot <- renderPlot({
  req(reactiveVals$de_results)
  de_results <- reactiveVals$de_results
  ain <- input$deVisAin
  comparison <- input$deVisComparison
  if (!is.null(ain)) {
    if ((length(ain) == 1) && (length(comparison) != 1)) {
      reactiveVals$de_volcano <- plot_volcano_DE(
        de_results,
        comparison = comparison,
        facet_norm = FALSE,
        facet_comparison = TRUE
      )[[1]] + ggtitle("")
    } else if (length(comparison) == 1) {
      reactiveVals$de_volcano <- plot_volcano_DE(de_results,
                                                 comparison = comparison,
                                                 facet_norm = TRUE)[[1]] + ggtitle("")
    }
    reactiveVals$de_volcano
  } else {
    reactiveVals$de_volcano <- NULL
    return(NULL)
  }
})


# Heatmap

output$de_heatmap_tab <- renderUI({
  de_results <- reactiveVals$de_results
  if (is.null(input$deVisAin) ||
      is.null(input$deVisComparison)) {
    HTML("<b>First Select Parameters and Click on The Button Perform DE Analysis<b>")
  } else if (length(input$deVisComparison) == 1 &&
             length(input$deVisAin) == 1) {
    de_res <- de_results[(
      de_results$Assay %in% input$deVisAin &
        de_results$Comparison %in% input$deVisComparison
    ), ]
    de_res <- de_res[de_res$Change != "No Change", ]
    if (nrow(de_res) == 0) {
      HTML("<b>No DE Results for the Selected Parameters<b>")
    } else {
      shinycssloaders::withSpinner(plotOutput("de_heatmap_plot", width = "100%", height = "750px"))
    }
  } else {
    HTML(
      "<b>This plot is only available when a single comparison and a single normalization method is selected.<b>"
    )
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
  filename = function() {
    paste0("DE_Heatmap_Plot", ".pdf")
  },
  content = function(file) {
    waiter_show(
      id = "app",
      html = tagList(spinner$logo, HTML("<br>Downloading...")),
      color = spinner$color
    )
    pdf(file, width = 12, height = 8)
    ComplexHeatmap::draw(reactiveVals$de_heatmap)
    dev.off()
    waiter_hide(id = "app")
  }
)

output$de_heatmap_plot <- renderPlot({
  req(reactiveVals$de_results)
  de_results <- reactiveVals$de_results
  ain <- input$deVisAin
  comparison <- input$deVisComparison
  condition <- reactiveVals$de_condition
  if (!is.null(ain)) {
    if ((length(ain) == 1) && (length(comparison) == 1)) {
      custom_colors <- reactiveVals$selected_palette
      reactiveVals$de_heatmap <- plot_heatmap_DE(
        se = reactiveVals$se,
        de_res = de_results,
        ain = ain,
        condition = condition,
        comparison = comparison,
        col_vector = custom_colors
      )[[ain]]
      reactiveVals$de_heatmap
    } else {
      reactiveVals$de_heatmap <- NULL
      return(NULL)
    }
  }
})
