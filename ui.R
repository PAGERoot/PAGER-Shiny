# Guillaume Lobet - University of Liege


library(shiny)
# library(shinyFiles)

shinyUI(fluidPage(
  
  # Application title
  titlePanel(h1("--| PAGE-Root |--")),
  
  fluidRow(
    column(3, wellPanel(
      helpText("PAGE-Root (Pattern Analysis of Gene Expression in Root) is a pipeline designed to analyse and visualize gene and reporter expression patterns in root."),

      h3("1. Load your data"),
      tags$hr(),      
      
      
      #textInput('dirPath', "Select folder containing the reporter xlsx files", placeholder = "Select folder"),
      #shinyDirButton('directory', 'Select folder', 'Please select a folder'),
      
      fileInput('rep_file', 'Choose reporter expression data File', accept=c('text/tab-separated-values', 'csv')),

      tags$hr(), 
      
      fileInput('gene_file', 'Choose gene expression data File', accept=c('text/tab-separated-values', '.txt')),

      checkboxInput('use_example', "Use example files", value = FALSE, width = NULL),

      actionButton(inputId = "load_data", label="Load and analyse data"),
      
      tags$hr(),      

      
      
      h3("2. Plot your data"),
      tags$hr(),      

      selectInput("ref_reps", label = "Reference genotype", choices = c("Please load datafile")), # updated with the datafile
      selectInput("to_plot", label = "Genotype to compare", choices = c("Please load datafile"), selected = 2), # updated with the datafile
      
      tags$hr(),      
      
      selectInput("ref_genes", label = "Promoter to analyse", choices = c("Please load datafile")), # updated with the datafile
      
      # tags$hr(),      
      helpText("Wait for the data loading and analysis to be done before hitting the plot button. 
               It will close the application if you do so. 
               The analysis is done when you can see the different line name in the dropdown menu above."),
      
      actionButton(inputId = "runROOTEXP", label="Plot data")
      
      )),  
    
    
    
    # Show a plot of the generated distribution
    column(9,
           tabsetPanel( 
             tabPanel("Reporter - Reporter",
                      
                      helpText("Expression pattern in the root"),
                      tags$hr(),                                  
                      fluidRow(
                        column(2,
                          h5(textOutput("nameWt")),
                          tags$hr(),
                          imageOutput("root_1", width = "100px")),
                        column(2,
                          h5(textOutput("nameMt")),
                          tags$hr(),
                          imageOutput("root_2", width = "100px")),
                        column(1,
                               h5("scale"),
                               tags$hr(),
                               plotOutput("scalePlot")),
                        column(2,
                          h5(textOutput("nameDiff")),
                          tags$hr(),
                          imageOutput("root_diff", width = "100px")),

                        column(5,
                               tags$hr(),                       
                               textOutput("line_comp_text"),
                               textOutput("line_comp_pval"),
                               tags$hr(),                       
                               h5("Different cell types:"),
                               textOutput("cell_comp"),
                               tags$hr(), 
                               tabsetPanel( 
                                 tabPanel("LDA plot",
                                          
                                    
                                          plotOutput("ldaPlot"),
                                          tags$hr(),                       
                                          helpText("Graphical representation of the two first dimentions on the 
                                                  Linear Discriminant Analysis done on the reporter lines dataset.
                                                  The plotted lines are shown in color, while the other lines are displayed in grey"),
                                          value=1
                                 ),
                                 tabPanel("Barplot",
                                          
                                          
                                          plotOutput("barplot_comp"),
                                          tags$hr(),                       
                                          helpText("Normalized level of fluorescence fopr the different cell layers"),
                                          value=1
                                          ),
                                 tabPanel("MAOV results",
                                          
                                          tags$hr(), 
                                          downloadButton('download_moav', 'Download'),
                                          tags$hr(), 
                                          tableOutput('maov_results'),  
                                          value=2
                                  ),  
                                 tabPanel("AOV results",
                                          
                                          tags$hr(), 
                                          downloadButton('download_oav', 'Download'),
                                          tags$hr(), 
                                          tableOutput('aov_results'),  
                                          value=2 
                                 )  
                                )
                        ),
                        width="100%",
                        height="100%"
                      ),
                      value=1
             ),              

            tabPanel("Reporter - Gene",
                     helpText("Expression pattern in the root"),
                     tags$hr(),    
                     fluidRow(
                       column(2,
                              h5(textOutput("name_gene_1")),
                              tags$hr(),
                              imageOutput("root_gene_1", width = "100px")),
                       column(2,
                              h5(textOutput("name_pred_1")),
                              tags$hr(),
                              imageOutput("root_pred_1", width = "100px")),
                       column(1,
                              h5("  "),
                              tags$hr(),
                              plotOutput("scalePlot2")),
                       column(5,
                              tags$hr(),                       
                              tags$hr(),                       

                              tags$hr(), 
                              tabsetPanel( 

                                tabPanel("Barplot",
                                         plotOutput("barplot_comp_1"),
                                         tags$hr(),                       
                                         helpText("Graphical representation of the two first dimentions on the 
                                                  Linear Discriminant Analysis done on the reporter lines dataset.
                                                  The plotted lines are shown in color, while the other lines are displayed in grey"),
                                         value=1
                                         ),
                                 tabPanel("Comparison results",
                                         
                                         tags$hr(), 
                                         downloadButton('download_comp', 'Download'),
                                         tags$hr(), 
                                         tableOutput('comp_results'),  
                                         value=2 
                                  )  
                                )
                              ),
                       width="100%",
                       height="100%"
                       ),
                     value=2
            ), 

            
            tabPanel("Promoter - Reporter",
                     helpText("Expression pattern in the root"),
                     tags$hr(),                                  
                     flowLayout(
                       verticalLayout(
                         h4(textOutput("nameProm")),
                         tags$hr(),
                         imageOutput("root_p_2", width = "120px")),
                       verticalLayout(
                         h4(textOutput("name_p_Wt")),
                         tags$hr(),
                         imageOutput("root_p_1", width = "120px")),
                       width="100%"
                     ),
                     value=3
            ), 

              tabPanel("Reporter lines comparisons",
                       fluidRow(
                         column(5,
                           plotOutput("heatmap", width = "100%"),
                           tags$hr(),                       
                           helpText("Heatmap of the MANOVA results between the different lines. 
                                    For each line combinaison, a MANOVA analysis was performed in order to determine if there
                                    expression patterns were different."),
                           downloadButton('downloadPlotHeat', 'Download Heatmap')                 
                         ),
                         column(5,
                           plotOutput("histoPlot", width = "100%"),
                           tags$hr(), 
                           helpText("Distribution of the expression data. Data from each replicate was normalized internally."),
                           downloadButton('downloadPlot', 'Download Plot')
                         )
                       ),                       
                       value=4
              ), 
#              tabPanel("Tissue Expression (detail)",
#                       helpText("Difference in expression for each tissue"),
#                       downloadButton('downloadPlot1', 'Download Plot1'),                 
#                       tags$hr(),                       
#                       plotOutput("diffPlot", width = "100%"),
#                       value=1
#              ),             
#              
#              tabPanel("Distribution of expression levels",
#                       helpText("Distribution of the expression data. Normalized in the right panel"),
#                       downloadButton('downloadPlot', 'Download Plot'),                 
#                       tags$hr(),                       
#                       plotOutput("histoPlot", width = "100%"),
#                       value=5
#              ),
             
   
             
             id="tabs1"
           )
    )
    )
))