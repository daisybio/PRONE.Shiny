goodbyeBody <-  tabItem(tabName = "goodbye", fluidRow(shinydashboard::box(
  fluidRow(
    column(
      width = 12,
      div(HTML("<b>Curious about the magic behind your data? We believe in reproducibility!</b></br>"), style = "font-size: large;font-weight: bold;"),
      uiOutput("reproducibilitySteps"),
      div(HTML("</br></br>"), style = "font-size: large;"),
      div(HTML("<b>Thank you for using PRONE!</b></br>"), style = "font-size: large;font-weight: bold;"),
      div(
        HTML(
          "We sincerely appreciate your trust in our tool to assist you in achieving reliable and accurate proteomics data normalization and performing proteomics data analysis."
        ),
        HTML(
          "We are committed to continuously improving PRONE and eagerly anticipate your valuable feedback, which is crucial for driving future enhancements. Therefore, please feel free to write an issue to the GitHub repository."
        ),
        HTML("</br></br>"),
        style = "font-size: large;"
      ),
      div(HTML("<b>How to cite PRONE?</b>"), style = "font-size: large;font-weight: bold;"),
      div(
        HTML(
          "
                                            <p style='margin-left: 20px;'>Systematic Evaluation of Normalization Approaches in Tandem Mass Tag and Label-Free Protein Quantification Data Using PRONE
                                            </p>
                                            <p style='margin-left: 20px;'>Lis Arend, Klaudia Adamowicz, Johannes R. Schmidt, Yuliya Burankova, Olga Zolotareva, Olga Tsoy, Josch K. Pauling, Stefan Kalkhof, Jan Baumbach, Markus List, Tanja Laske
                                            </p>
                                            <p style='margin-left: 20px;'>bioRxiv 2025.01.27.634993; doi: <a href='https://doi.org/10.1101/2025.01.27.634993'>https://doi.org/10.1101/2025.01.27.634993</a>
                                            </p>"
        ),
        style = "font-size: large;"
      )
    ),
    column(
      width = 12,
      img(
        src = "PRONE_Text_Logo.png",
        height = "170px",
        width = "auto"
      ),
      style = "vertical-align: middle; text-align:center;"
    )
  ),
  title = h2(
    "Thank you for using PRONE, the PROteomics Normalization Evaluator."
  ),
  width = 12
), width = 12))
