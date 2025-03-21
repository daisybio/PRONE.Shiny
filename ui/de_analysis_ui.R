deAnalysisBody <- function() {
  box_height <- "24em"
  
  DEInputPopover <- bsPopover(id = "DEInputInfo",
                              title = "Input data to use for performing DE analyis",
                              content = "DE analayis can be performed simulateously on different normalized intensities. Select the normalizations to compare their ability of detecting DE proteins.")
  
  DEMethodPopover <- bsPopover(id = "DEMethodInfo",
                               title = "Method for performing DE analysis",
                               content = "In limma, a separate linear model is fitted for each gene and an additional parametric empirical Bayes approach is applied. ROTS optimizes the parameters of a family of modified t-statistics by maximizing the reproducibility of the top-ranked features in group-preserving bootstrap datasets.")
  
  DEColumnPopover <-  bsPopover(id = "DEColumnInfo",
                                title = "Column of metadata specifiying the condition to analyze",
                                content = "You need to specify the column of the metadata for which you want to compare sample groups.")
  
  DEComparisonPopover <-  bsPopover(id = "DEComparisonInfo",
                                    title = "Comparisons that should be executed.",
                                    content = "DE analysis is typically performed pairwise. Select here all pairwise comparisons that should be evaluated. To ensure clarity, a comparison in the form of Condition-Reference, where the logFC is positive, signifies that the protein intensity is greater in the Condition compared to the Reference.")
  
  DEROTSBPopover <-  bsPopover(id = "DEROTSBInfo",
                               title = "Number of bootstrapping.",
                               content = "An integer specifying the number of bootstrap and permutation resamplings. Default = 1000")
  
  DEROTSKPopover <-  bsPopover(id = "DEROTSKInfo",
                               title = "Number indicating the largest top list size considered.",
                               content = "Lowering this size can markedly reduce computation time. In large data matrices with thousands of proteins, use a size of several thousands. Default = 1/4 of the features are used.")
  
  DEEvalMethodPopover <- bsPopover(id = "DEEvalMethodInfo",
                                   title = "Select normalization methods you want to include in evaluation.",
                                   content = "If you performed DE analysis on several normalization methods, you can here easily deselect some for visualizations.")
  
  DEEvalComparisonPopover <- bsPopover(id = "DEEvalComparisonInfo",
                                       title = "Select the pairwise comparison you want to analyze in more detail.",
                                       content = "The performance of the normalization methods is available separately for each pairwise comparison that was executed in DE analysis.")
  
  DEEvalHeatmapPopover <- bsPopover(id = "DEEvalHeatmapInfo",
                                    title = "Heatmap of normalized intensities with p-values or p-adjusted as row annotation.",
                                    content = "Since the heatmap can only include normalized intensities of one normalization approach, you need to switch between methods using the settings button.")
  
  DEBiomarkerColumnPopover <- bsPopover(id = "DEBiomarkerColumnInfo",
                                        title = "Select the column of which you want to select the biomarkers.",
                                        content = "Protein.IDs or Gene.Names can be selected here, and the biomarker selection menu automatically updates.")
  
  DEBiomarkerSelectionPopover <- bsPopover(id = "DEBiomarkerSelectionInfo",
                                           title = "Select the biomarkers you want to evaluate the normalization methods on.",
                                           content = "The performance of normalization method can be assessed by selecting known disease-related biomarkers and checking how many of these were detected as DEP by each normalization method.")
  
  DEBiomarkerCoveragePopover <- bsPopover(id = "DEBiomarkerCoverageInfo",
                                          title = "The coverage of significant biomarkers is reported for each selected normalization approach.",
                                          content =  "The number of biomarkers found as significantly DE by each normalization method is divided by the total number of selected biomarkers.")
  
  DEBiomarkerBoxplotPopover <- bsPopover(id = "DEBiomarkerBoxplotInfo",
                                         title = "The normalized intensities of each biomarker are plotted, colored and grouped by the sample groups of the selected pairwise comparison.",
                                         content =  "A biomarker needs to be selected separately here since the plots are only rendered for one biomarker at the time.")
  
  DEIntersectionPlotPopover <- bsPopover(id = "DEIntersectionPlotInfo",
                                         title = "All DEPs of the specified pairwise comparison of all normalization methods are intersected.",
                                         content =  "The minimal degree of an intersection for it to be included is 2. All results can be extracted from the table at the next tab.")
  
  DEIntersectionTablePopover <- bsPopover(id = "DEIntersectionTableInfo",
                                          title = "All DEPs of the specified pairwise comparison of all normalization methods are intersected.",
                                          content =  "For each protein, the number of intersected normalization methods and the methods themselves are provided.")
  
  DEDEqMSColumnPopover <- bsPopover(id = "DEDEqMSColumnInfo",
                                    title = "Column for accounting for variance in DEqMS",
                                    content = "DEqMS is built on top of limma and accounts for variance dependence on the number of quantified peptides or PSMs for statistical testing in differential protein expression. Select the column of the data that should be used to correct bias of variance in DEqMS.")
  
  DEIntersectionJaccardPopover <- bsPopover(id = "DEIntersectionJaccardInfo",
                                            title = "Jaccard index of the intersection of DEPs.",
                                            content = "The Jaccard index is a measure of similarity between two sets. It is defined as the size of the intersection divided by the size of the union of the sets.")
  
  DEVolcanoPopover <- bsPopover(id = "DEVolcanoInfo",
                                title = "Volcano plot of DE results.",
                                content = "This plot can only be generated for either a single normalization method or a single comparison. If one normalization method is selected, the plot is facet by comparison and vice versa when selecting a single comparison.")
  
  DELimmaTrendPopover <- bsPopover(id = "DELimmaTrendInfo",
                                   title = "Limma trend parameter.",
                                   content = "Should an intensity-dependent trend be allowed for the prior variance? If FALSE then the prior variance is constant.")
  
  DELimmaRobustPopover <- bsPopover(id = "DELimmaRobustInfo",
                                    title = "Limma robust parameter.",
                                    content = "Should the estimation of df.prior and var.prior be robustified against outlier sample variances?")
  
  DEConsensusDEPPopover <- bsPopover(id = "DEConsensusDEPInfo",
                                     title = "Use a consensus set of DEPs.",
                                     content = "If normalization methods agree on a high amount of DEPs, a consensus set of DE proteins can be used for further analysis.")

  deSelectionBox <- shinydashboard::box(
    width = 12,
    height = "40em",
    title = "Run DE Analysis",
    fluidRow(
      column(
        width = 6,
        pickerInput(
          "deNormInput",
          label = span(
            "Select the input data for DE analysis.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEInputInfo"
          ),
          multiple = TRUE,
          choices = NULL,
          options = list("actions-box" = TRUE)
        ),
        pickerInput(
          "deMethodInput",
          label = span(
            "Select the DE analysis method.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEMethodInfo"
          ),
          choices = c("ROTS", "limma", "DEqMS"),
          multiple = FALSE,
          selected = "limma"
        ),
        numericInput(
          "deROTSB",
          label = span(
            "Select the number of bootstrap and permutation resamplings.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEROTSBInfo"
          ),
          value = 1000,
          step = 100,
          min = 100,
          max = 3000
        ),
        pickerInput(
          "deDEqMSColumn",
          label = span(
            "Select the column for accounting for variance in DEqMS.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEDEqMSColumnInfo"
          ),
          choices = NULL,
          multiple = FALSE,
        )
      ),
      column(
        width = 6,
        pickerInput(
          "deColComparison",
          label = span(
            "Select the column for comparison.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEColumnInfo"
          ),
          choices = NULL,
          multiple = FALSE,
        ),
        pickerInput(
          "deComparison",
          "Select Comparisons: ",
          label = span(
            "Select comparisons.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEComparisonInfo"
          ),
          choices = NULL,
          multiple = TRUE,
          options = list("actions-box" = TRUE)
        ),
        numericInput(
          "deROTSK",
          label = span(
            "Specify the largest top list size considered.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEROTSKInfo"
          ),
          value = 500,
          step = 100,
          min = 100,
          max = 3000
        ),
        radioGroupButtons(
          inputId = "deLimmaTrend",
          label = span(
            "Limma trend parameter.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DELimmaTrendInfo"
          ),
          choices = c("TRUE" = TRUE, "FALSE" = FALSE),
          selected = TRUE,
          status = "s"
        ),
        radioGroupButtons(
          inputId = "deLimmaRobust",
          label = span(
            "Limma robust parameter.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DELimmaRobustInfo"
          ),
          choices = c("TRUE" = TRUE, "FALSE" = FALSE),
          selected = TRUE,
          status = "s"
        ),
        DEInputPopover,
        DEMethodPopover,
        DEColumnPopover,
        DEComparisonPopover,
        DEROTSBPopover,
        DEROTSKPopover,
        DEDEqMSColumnPopover,
        DELimmaTrendPopover,
        DELimmaRobustPopover,
        DEConsensusDEPPopover
      ),
    ),
    fluidRow(
      style = "padding-top:30px;",
      column(
        width = 6,
        materialSwitch(
          "deLogFC",
          value = TRUE,
          label = HTML("<b>Do you want to set a threshold for logFC?</b>"),
          status = "primary",
        ),
        numericInput(
          "deLogFCThr",
          label = "Select threshold for logFC",
          value = 1,
          step = 0.5,
          min = 0,
          max = 10,
        ),
      ),
      column(
        width = 6,
        materialSwitch(
          "deMultTest",
          value = TRUE,
          label = HTML("<b>Do you want to adjust for multiple testing?</b>"),
          status = "primary"
        ),
        column(
          width = 6,
          style = "padding-right: 5px; padding-left: 0px",
          numericInput(
            "dePThr",
            label = "Select threshold for p-value or p.adjust respectively",
            value = 0.05,
            step = 0.05,
            min = 0,
            max = 0.1
          ),
        )
      ),
      
      
    ),
    
    div(
      bsButton(
        "performDEAnalysis",
        " Perform DE Analysis ",
        icon = icon("object-group"),
        style = "success",
        disabled = TRUE,
      ),
      style = "float: right; position: absolute; bottom: 15px;"
    )
  )
  
  
  deTableBox <- shinydashboard::box(
    width = 12,
    id = "deTable",
    collapsible = TRUE,
    collapsed = TRUE,
    title = "Table of DE Results",
    fluidRow(
      style = "padding-left: 20px; padding-right: 20px;",
      uiOutput("de_table_results"),
      uiOutput("de_table_results_download")
    )
  )
  
  deOverviewBox <- shinydashboard::box(
    width = 12,
    id = "deOverviewVis",
    title = "Overview of DE Results",
    collapsible = TRUE,
    collapsed = TRUE,
    fluidRow(
      tabBox(
        tabPanel(
          uiOutput("de_overview_bar_tab"),
          uiOutput("de_overview_bar_download"),
          value = "de_overview_1",
          title = "Bar Plot"
        ),
        tabPanel(
          uiOutput("de_overview_tile_tab"),
          uiOutput("de_overview_tile_download"),
          value = "de_overview_2",
          title = "Tile Plot",
        ),
        id = "deOverviewTab",
        title = "Overview DE Plots",
        width = 12
      )
    )
  )
  
  deResCompVisBox <- shinydashboard::box(
    width = 12,
    id = "deCompVis",
    title = "Comparison- or Assay-Specific Visualizations",
    collapsible = TRUE,
    collapsed = TRUE,
    fluidRow(
      tabBox(
        tabPanel(
          uiOutput("de_volcano_tab"),
          uiOutput("de_volcano_plot_download"),
          DEVolcanoPopover,
          value = "de_volcano_1",
          title = span(
            "Volcano Plot",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEVolcanoInfo"
          )
        ),
        tabPanel(
          uiOutput("de_heatmap_tab"),
          uiOutput("de_heatmap_plot_download"),
          DEEvalHeatmapPopover,
          value = "de_heatmap",
          title = span(
            "Heatmap",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEEvalHeatmapInfo"
          ),
        ),
        id = "quan_norm_plots",
        title = "Comparison- and Assay-Specific DE Plots",
        width = 12
      )
    )
  )
  
  deIntersectionBox <- shinydashboard::box(
    width = 12,
    id = "deIntersection",
    title = "Intersection Analysis of DE Results",
    collapsible = TRUE,
    collapsed = TRUE,
    fluidRow(
      tabBox(
        tabPanel(
          uiOutput("de_intersection_plot_tab"),
          DEIntersectionPlotPopover,
          uiOutput("de_intersection_plot_download"),
          value = "de_intersection_1",
          title = span(
            "Intersection Plot",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEIntersectionPlotInfo"
          )
        ),
        tabPanel(
          uiOutput("de_intersection_table_tab"),
          DEIntersectionTablePopover,
          uiOutput("de_intersection_table_download"),
          value = "de_intersection_2",
          title = span(
            "Intersection Table",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEIntersectionTableInfo"
          )
        ),
        tabPanel(
          uiOutput("de_intersection_jaccard_tab"),
          DEIntersectionJaccardPopover,
          uiOutput("de_intersection_jaccard_plot_download"),
          value = "de_intersection_3",
          title = span(
            "Jaccard Similarity Coefficient",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEIntersectionJaccardInfo"
          )
        ),
        tabPanel(
          uiOutput("de_intersection_consensus_tab"),
          uiOutput("de_intersection_consensus_download"),
          value = "de_intersection_4",
          title = span(
            "Consensus DEPs",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEConsensusDEPInfo"
          )
        ),
        id = "DE_biomarker_plots",
        title = "Intersection Plots",
        width = 12
      )
    )
    
  )
  
  deEnrichmentBox <- shinydashboard::box(
    width = 12,
    id = "deEnrichment",
    title = "Functional Enrichment of DE Results",
    collapsible = TRUE,
    collapsed = TRUE,
    uiOutput("de_enrichment_tab")
  )
  
  deEvaluationBox <- shinydashboard::box(
    width = 12,
    id = "deEvaluation",
    title = "Evaluation of DE Results",
    collapsible = TRUE,
    collapsed = TRUE,
    fluidRow(
      column(
        width = 6,
        pickerInput(
          "deVisAin",
          choices = NULL,
          multiple = TRUE,
          options = list("actions-box" = TRUE),
          label = span(
            "Select the normalization methods you want to evaluate.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEEvalMethodInfo"
          )
        ),
        DEEvalMethodPopover,
      ),
      column(
        width = 6,
        pickerInput(
          "deVisComparison",
          choices = NULL,
          multiple = TRUE,
          label = span(
            "Select the comparisons you want to visualize.",
            tags$i(class = "fa fa-circle-info", style = "color: rgb(0,166,90)"),
            id = "DEEvalComparisonInfo"
          )
        ),
        DEEvalComparisonPopover,
      )
    ),
    fluidRow(deOverviewBox),
    fluidRow(deResCompVisBox),
    fluidRow(deIntersectionBox),
    fluidRow(deTableBox),
    fluidRow(deEnrichmentBox),
  )
  
  
  
  deanalysis <- tabItem(
    tabName = "de_analysis",
    fluidRow(shinydashboard::box(
      div(
        HTML(
          "Among the various computational methods utilized in the examination of biological data, DE analysis is a highly prevalent approach that leverages statistical inference to
        identify genes or proteins exhibiting significant changes in expression across different sample groups.
        <br />
        This study adopts two widely used algorithms, limma and ROTS, to compare the effectiveness of normalization methods in identifying DE features. limma is a popular R/Bioconductor package developed primarily for microarray data but lately extended to other high-throughput omics data. A separate linear model is fitted for each gene and an additional parametric empirical Bayes approach is particularly useful for studies with small sample sizes. The R/Bioconductor package ROTS addresses the challenge of selecting an appropriate test statistic and optimizes the parameters of a family of modified t-statistics by maximizing the reproducibility of the top-ranked features in group-preserving bootstrap datasets.
        <br />
        An overview barplot is rendered for easy comparison of DE results of the different normalization methods.
        <br />
        The visualization tab provides volcano plots and heatmaps for a more detailed inspection of specific comparisons. Additionally, the intersection of DEPs of multiple normalized intensities are provided for comparison of different normalization methods. Finally, functional enrichment analysis can be performed on the DE results of the different normalization methods to biologically interpret the DE results of the different normalization methods."
        )
      ),
      title = h2("DE Analysis"),
      width = 12
    )),
    fluidRow(deSelectionBox),
    fluidRow(deEvaluationBox)
  )
  return(deanalysis)
  
}
