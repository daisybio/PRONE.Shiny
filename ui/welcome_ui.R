welcomeBody <-  tabItem(tabName = "welcome",
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
                                            A detailed description of the PRONE package that is also useful for the navigation through the Shiny app is available <a href='https://daisybio.github.io/PRONE/'>here</a>.
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
                                  
)