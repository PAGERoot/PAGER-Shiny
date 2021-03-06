#
# Copyright © 2017, Université catholique de Louvain
# All rights reserved.
# 
# Copyright © 2017 Forschungszentrum Jülich GmbH
# All rights reserved.
# 
# Developers: Guillaume Lobet
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted under the GNU General Public License v3 and provided that the following conditions are met:
#   
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# 
# Disclaimer
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# You should have received the GNU GENERAL PUBLIC LICENSE v3 with this file in license.txt but can also be found at http://www.gnu.org/licenses/gpl-3.0.en.html
# 
# NOTE: The GPL.v3 license requires that all derivative work is distributed under the same license. That means that if you use this source code in any other program, you can only distribute that program with the full source code included and licensed under a GPL license.



library(shiny)
# library(shinyFiles)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("PAGE-Root"),
  
  fluidRow(
    column(3, 
      h3("1. Load your data"),
      wellPanel(
      tabsetPanel(
        tabPanel("Your data",
             fileInput('rep_file', 'Choose reporter expression data File', accept=c('text', '.rsml')),
             fileInput('gene_file', 'Choose gene expression data File', accept=c('text/comma-separated-values', '.csv')),
             h3("2. Choose the options"),
             checkboxInput('use_absolute', "Use absolute values (no scaling)", value = FALSE, width = NULL),
             checkboxInput('log2', "Log2 transform your data", value = TRUE, width = NULL),
             selectInput("type_to_analyse", label="Cell types to use in analysis", choices = c("Load datafile"), selected = NULL, multiple = TRUE),
             selectInput("method", label = "How the average the data", choices = c("Mean", "Median", "Min", "Max"))

             # selectInput("method", label = "Method used for the analysis", choices = c("Mean", "Median", "Min", "Max")), # updated with the datafile
             # helpText("Define which method to use to aggregate the data at the line x root x cell type level")
        ),
        tabPanel("Sample data",
              checkboxInput('use_example', "Use example data", value = F, width = NULL),
              selectInput("reporters", label = "Select reporter dataset", choices = c("Load datafile")), # updated with the datafile
              htmlOutput("littTitle"),
              htmlOutput("littAuth"),
              htmlOutput("littRef"),
              htmlOutput("doi"),
              checkboxInput('log2bis', "Use log2 transform data", value = TRUE, width = NULL),
              tags$hr(),
              selectInput("microarrays", label = "Select microarray dataset", choices = c("Load datafile")),
              htmlOutput("littTitle1"),
              htmlOutput("littAuth1"),
              htmlOutput("littRef1"),
              htmlOutput("doi1")#,
        )# updated with the datafile
        )
      ),
      tags$hr(),
      actionButton(inputId = "load_data", label=" Launch PAGE-Root", icon("paper-plane"), 
                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4")

    ),  
    
    
#------------------------------------------------------------------
#------------------------------------------------------------------
    
    # Show a plot of the generated distribution
    column(9,
           tabsetPanel( 
             tabPanel("Compare reporters expressions",
                      tags$hr(),
                      fluidRow(
                        
                        column(6,
                               fluidRow(
                                 selectInput("type_to_plot", label="Cell types to plot", choices = c("Load datafile"), 
                                             selected = NULL, multiple = TRUE, width="100%")
                               ),
                               fluidRow(
                                 column(5,
                                        selectInput("ref_reps", label = NULL, choices = c("Load datafile"))
                                        ),
                                 column(5,
                                        selectInput("to_plot", label = NULL, choices = c("Load datafile"), selected = 2)
                                 ),
                                 column(2,
                                        checkboxInput('show_diff', "", value = FALSE, width = NULL)
                                 )
                               ),
                               fluidRow(
                                plotOutput("plotRoot", width="100%")#,
                               )
                               # sliderInput("display_range", "Display range:",min = 0, max = 1, value = c(0,1))
                               
                        ),
                        column(6,
                               h4("Overall difference between the selected reporter lines"),
                               helpText("Results from the MANOVA analysis, revealing if there is a global differences between both lines"),
                               textOutput("line_comp_text"),
                               textOutput("line_comp_pval"),
                               textOutput("line_comp_fit"),
                               textOutput("line_comp_pears"),
                               textOutput("line_comp_spear"),
                               tags$hr(),
                               h4("Single differences between the root types of the selected reporter lines"),
                               helpText("It should be noted, that even if two lines are not globally different, they might have local, single cell-type differences."),
                               textOutput("cell_comp"),
                               tags$hr(),
                               tabsetPanel(
                                 tabPanel("Boxplot",
                                          plotOutput("barplot_comp"),
                                          tags$hr(),
                                          helpText("Normalized level of fluorescence for the different cell layers"),
                                          value=1
                                          ),
                                 tabPanel("Regression plot",
                                          plotOutput("fitPlot"),
                                          tags$hr(),
                                          helpText("Comparison of the fluorescence elvels for both reporter lines. The dotted represents the 1:1 diagonal."),
                                          value=1
                                          ),                                 
                                 tabPanel("PCA plot",
                                          plotlyOutput("ldaPlot"),
                                          tags$hr(),
                                          helpText("Graphical representation of the two first dimentions on the
                                                   Linear Discriminant Analysis done on the reporter lines dataset.
                                                   The plotted lines are shown in color, while the other lines are displayed in grey"),
                                          value=1
                                          ),

                                 tabPanel("AOV results",
                                          tags$hr(),
                                          helpText("This table contains all the comparison between the different reporter lines, for each cell type, using 2-way ANOVA analyses"),
                                          downloadButton('download_oav', 'Download table'),
                                          tags$hr(),
                                          tableOutput('aov_results'),
                                          value=2
                                 )
                              )
                        )
                      ),
                      value=1
             ),              

#------------------------------------------------------------------
#------------------------------------------------------------------


                                 

          tabPanel("Correlations heatmaps",
                   tabsetPanel(
                     tabPanel("MANOVA",
                       plotOutput("heatmap", width = "100%", height="800px"),
                       tags$hr(),
                       helpText("Heatmap of the MANOVA results between the different lines.
                                                  For each line combinaison, a MANOVA analysis was performed in order to determine if there
                                                  expression patterns were different.")         
                     ),
                    tabPanel("Linear regression",
                      plotOutput("heatmap_fit", width = "100%", height="800px"),
                      tags$hr(),
                      helpText("Heatmap of the r-squared results between the different lines.
                               For each line combinaison, a linear regression analysis was performed in order to determine if there
                               expression patterns were different.")                                     
                    ),
                    tabPanel("Pearson correlations",
                             plotOutput("heatmap_pearson", width = "100%", height="800px"),
                             tags$hr(),
                             helpText("Heatmap of the pearson correlation results between the different lines.")                                     
                    ),
                    tabPanel("Spearman correlations",
                             plotOutput("heatmap_spearman", width = "100%", height="800px"),
                             tags$hr(),
                             helpText("Heatmap of the spearman correlation results between the different lines.")                                     
                             )
                   )
          ),
#------------------------------------------------------------------
#------------------------------------------------------------------



             
            tabPanel("Compare Reporter to Gene",
                     tags$hr(),
                     fluidRow(
                       column(6,
                          fluidRow(
                            column(5,
                                   selectInput("ref_reps_2", label = "Reporter", choices = c("Load datafile"))
                            ),
                            column(7,
                                   selectInput("type_to_plot_2", label="Cell types to plot", choices = c("Load datafile"), 
                                               selected = NULL, multiple = TRUE, width="100%")
                            )
                          ),
                          fluidRow(
                            plotOutput("plotRootGene", width="100%")
                          )
                          # sliderInput("display_range_1", "Display range:",min = 0, max = 1, value = c(0,1))
                          
                       ),
                       column(6,
                              h4("Overall difference between the selected reporter lines"),
                              textOutput("gene_comp_fit"),
                              textOutput("gene_comp_pears"),
                              textOutput("gene_comp_spear"),                              
                              h4("Single differences between the root types of the selected reporter lines"),
                              textOutput("gene_comp"),
                              tags$hr(),
                              tabsetPanel(
                                tabPanel("Boxplot",
                                         plotOutput("barplot_comp_1"),
                                         tags$hr(),
                                         helpText("Graphical representation of the two first dimentions on the
                                                  Linear Discriminant Analysis done on the reporter lines dataset.
                                                  The plotted lines are shown in color, while the other lines are displayed in grey"),
                                         value=1
                                         ),
                                tabPanel("Regression plot",
                                         plotOutput("fitPlot_1"),
                                         tags$hr(),
                                         helpText("Comparison of the fluorescence and experession levels for the selected gene.The dotted represents the 1:1 diagonal."),
                                         value=1
                                )                        
                                )
                              ),
                       width="100%",
                       height="100%"
                       ),
                     value=2
            ),

#------------------------------------------------------------------
#------------------------------------------------------------------

              tabPanel("Download processed data",
                  tabsetPanel(
                     
                     tabPanel("MAOV results",
                              helpText("This table contains all the comparison between the different reporter lines, using MANOVA analyses"),
                              downloadButton('download_moav_all', 'Download full table'),
                              tags$hr(),
                              tableOutput('maov_results_all'),
                              value=2
                     ),
                     
                     tabPanel("Correlation results",
                              helpText("This table contains all the comparison between the different reporter lines, using a linear regression, Pearson and Spearman correlations"),
                              downloadButton('download_fit_all', 'Download full table'),
                              tags$hr(),
                              tableOutput('fit_results_all'),
                              value=2
                     ), 
                     
                     # tabPanel("Pearson correlation results",
                     #          helpText("This table contains all the comparison between the different reporter lines, using Pearson correlation"),
                     #          downloadButton('download_pears_all', 'Download full table'),
                     #          tags$hr(),
                     #          tableOutput('pears_results_all'),
                     #          value=2
                     # ), 
                     # 
                     # tabPanel("Spearman correlation results",
                     #          helpText("This table contains all the comparison between the different reporter lines, using Spearman correlation"),
                     #          downloadButton('download_spear_all', 'Download full table'),
                     #          tags$hr(),
                     #          tableOutput('spear_results_all'),
                     #          value=2
                     # ),                      
                     
                     tabPanel("AOV results",
                              helpText("This table contains all the comparison between the different reporter lines, for each cell type, using 2-way ANOVA analyses"),
                              downloadButton('download_oav_all', 'Download full table'),
                              tags$hr(),
                              tableOutput('aov_results_all'),
                              value=2
                     )  ,
                     
                     tabPanel("Gene AOV results",
                              helpText("This table contains all the comparison between the different reporter lines and their corresponding gene, for each cell type, using 2-way ANOVA analyses"),
                              downloadButton('download_oav_gene_all', 'Download full table'),
                              tags$hr(),
                              tableOutput('aov_gene_results_all'),
                              value=2
                     ),    
                     
                     tabPanel("Gene correlation results",
                              helpText("This table contains all the comparison between the different reporter lines and their corresponding gene, using linear regression, Pearson correlation and Spearman correlation"),
                              downloadButton('download_fit_gene_all', 'Download full table'),
                              tags$hr(),
                              tableOutput('fit_gene_results_all'),
                              value=2
                     )                         
                   ),                       
                   value=3
              ),


#------------------------------------------------------------------
#------------------------------------------------------------------

              tabPanel("About",
                  tags$hr(),
                   fluidRow(
                     column(4,
                            h3("Cell types legend"),
                            plotOutput("plotRootKey", height="800px", width="100%")
                     ),
                     column(4,
                        h3("What is PAGE-Root?"),
                        helpText("PAGE-Root (Pattern Analysis of Gene Expression in Root) is a pipeline designed to analyse and visualize gene and reporter expression patterns in root.")
                     ),
                     column(4,
                        h3("DISCLAIMER"),
                        helpText("THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.")
                     )
                   ),
                   value=4
              ),
             id="tabs1"
           )
    )
    )
))
