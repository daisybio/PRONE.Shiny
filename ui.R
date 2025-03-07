# Define UI for application that draws a histogram

library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyWidgets)
library(shinyBS)
library(shinyalert)
library(waiter)
library(shinycssloaders)

jscode <- "shinyjs.closewindow = function() { window.close(); }"

# read all ui files
sapply(list.files("ui", full.names = TRUE), source, environment())

# init spinner
options(spinner.color = "#95AAD3")
spinner <- list(logo = list(spin_loaders(id = 5, color = "#95AAD3")), color="rgb(153, 199, 191, .5)")



header <-
  dashboardHeader(
    title = span(img(src="PRONE_Text_Logo.png", height = "60px")),
    tags$li(
      actionLink(
        inputId = "change_colors",
        label = div(icon("palette", style = "color:white;"), " Change Colors"),
      ),
      class = "dropdown"
    ),
    tags$li(
      actionLink(
        inputId = "download_data",
        label = div(icon("download", style = "color:white;"), " Download Data"),
      ),
      class = "dropdown"
    ),
    tags$li(
      a(
        div(icon("circle-info", style = "color:white;"), " Vignette"),
        href="https://daisybio.github.io/PRONE/",
        title="",
        icon = icon("circle-info"),
        target = "_blank"
      ),
      class = "dropdown"
    ),
    tags$li(
      a(
        div(icon("question", style = "color:white;"), " Issues"),
        href="https://github.com/lisiarend/PRONE.Shiny/issues",
        title="",
        icon = icon("question"),
        target = "_blank"
      ),
      class = "dropdown"
    ),
    tags$li(
      a(
        div(icon("github", style = "color:white;"), " Github"),
        href="https://github.com/lisiarend/PRONE.Shiny",
        title="",
        target = "_blank"
      ),
      class = "dropdown"
    ),
    tags$li(
      actionLink(
        inputId = "dependencies_modal",
        label = div(icon("shield", style = "color:white;"), " Dependencies"),
      ),
      class = "dropdown"
    )

  )


sidebar <- dashboardSidebar(useShinyjs(),
                            extendShinyjs(text = jscode, functions = c("closewindow")),
                            tags$head(
                              tags$style(".inactiveLink {
                            pointer-events: none;
                           cursor: default;
                           }"),
                            ),
                            uiOutput("sidebar"),
                            collapsed = FALSE)

body <-
  dashboardBody(
    shinyalert::useShinyalert(force = TRUE),
    shinyjs::useShinyjs(),
    useWaiter(), # include dependencies
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
      tags$script(src = "custom.js")
    ),
    tags$style( # manually make waiter overly for the whole page
      ".waiter-overlay-content{
      color: white;
      font-weight: bold;
      height: 100vh;
      position: absolute;
      top: 50vh; /*30 pixels from the top*/
      right: 50%; /*48% from the right*/
      }

      .btn-custom{
        background-color: #95AAD3 !important;
      }

      .download-butt{
      background: #95AAD3 !important;
      color: white;
      }


      "
    ),
    waiter_show_on_load(html = tagList(spinner$logo,
                                       HTML("<br>Loading App ...")),
                        color=spinner$color),
    div(
      id = "loading",
      fluidRow(
        shinydashboard::box(fluidRow(
          column(width = 6,
                 div(
                   HTML(
                     "
                                            Latest, high-throughput technologies, such as DNA microarrays or mass spectrometry, have made substantial advancements in several ways, including instrument detection accuracy and data generation speed.
                                            These developments result in massive amounts of information-rich transcriptomics, proteomics, and metabolomics data.
                                            However, high-throughput OMICs data frequently comprise systematic biases introduced throughout various steps of a clinical study, from biological sample collection to quantification.
                                            Neglecting these biases could result in erroneous conclusions drawn from quantitative analysis.
                                            </br>
                                            </br>
                                            Data pre-processing techniques, in particular, normalization of data post-acquisition, aim to account for these biases and improve sample comparability.
                                            There are several approaches for processing and normalizing OMICs data generally and mass spectrometry (MS)-based proteomics data specifically.
                                            However, since the origin of these biases is usually unknown, selecting an appropriate normalization technique for a given dataset is challenging.
                                            </br>
                                            </br>
                                            Here, we present PRONE, a user-friendly R package that comes with a Shiny app that employs state-of-the-art normalization methods and enables simple evaluation of
                                            normalization methods through both quantitative and qualitative evaluation metrics and DE analysis.
                                            </br>
                                            </br>
                                            A detailed description of the PRONE package that is also useful for the navigation through the Shiny app is available <a href='https://lisiarend.github.io/PRONE/'>here</a>.
                                            </br>
                                            </br>
                                            If you are using either the R package or the Shiny app, please cite the following paper:
                                            </br>
                                            </br>
                                            <p style='margin-left: 10px;'>Systematic Evaluation of Normalization Approaches in Tandem Mass Tag and Label-Free Protein Quantification Data Using PRONE
                                            </p>
                                            <p style='margin-left: 10px;'>Lis Arend, Klaudia Adamowicz, Johannes R. Schmidt, Yuliya Burankova, Olga Zolotareva, Olga Tsoy, Josch K. Pauling, Stefan Kalkhof, Jan Baumbach, Markus List, Tanja Laske 
                                            </p>
                                            <p style='margin-left: 10px;'>bioRxiv 2025.01.27.634993; doi: <a href='https://doi.org/10.1101/2025.01.27.634993'>https://doi.org/10.1101/2025.01.27.634993</a>
                                            </p>
                     "
                   ),
                   style = "font-size: large;"
                 )
          ),
          column(
            width = 6,
            img(
              src = "PRONE_Workflow.png",
              height = "auto",
              width = "100%"
            ),
            style = "vertical-align: middle; text-align:center;"
          )
        ),
        title = h2("Welcome to PRONE, the PROteomics Normalization Evaluator."),
        width = 12
        ),
        width = 12)
  ),
    tabItems(
      welcomeBody,
      uploadBody(),
      filterProteinsBody(),
      outlierDetectionBody(),
      normalizationBody(),
      deAnalysisBody(),
      goodbyeBody
    ),
    fluidRow(column(
      6,
      shinyjs::hidden(bsButton(
        inputId = "previousTab",
        label = "",
        icon = icon("spinner"),
        style = "warning",
        block = TRUE,
        disabled = TRUE
      ))
    ),
    column(
      6,
      bsButton(
        inputId = "nextTab",
        label = "Get Started",
        icon = icon("spinner"),
        style = "warning",
        block = TRUE,
        disabled = TRUE
      )
    ))
  )


ui <- tags$div(id = "app", dashboardPage(header, sidebar, body, title = "PRONE", skin = "black"))
