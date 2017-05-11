#User Interface for PALO FILTER Script!
#options(shiny.deprecation.messages = FALSE)
#require(shiny)
library(shinydashboard)
#library(rmarkdown)

function(request){

header <- dashboardHeader(title = "Biostatistical Analysis Reporting Tool (BART)",
                          dropdownMenuOutput("messageMenu"),
                          titleWidth = 475)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("About BART", icon = icon("info-circle"), tabName = "intro"),
    menuItem("Upload", icon = icon("upload"), tabName = "upload"),
    menuItemOutput("SumStat"),
    menuItemOutput("Unsupervised"),
    menuItemOutput("diffge"),
    menuItemOutput("qusage"),
    menuItemOutput("roast"),
    menuItemOutput("flow.data"),
    menuItemOutput("metab.data"),
    menuItemOutput("correlations")
  )
)

body <-  dashboardBody(
  shinyjs::useShinyjs(),
  tabItems(
    uiOutput("css_colors"),
    tabItem(tabName = "intro",
            htmltools::includeMarkdown("bart-vignette.md")),
    tabItem(tabName = "upload",
            sidebarLayout(         
              sidebarPanel(
                verbatimTextOutput('version'),
                hr(),
                fileInput('file1', 'Upload BART Files and any QC reports for the Project:', multiple = TRUE,
                          accept=c(".Rdata",".png",".html",".csv")),
                bookmarkButton(),
                tags$script('
                            Shiny.addCustomMessageHandler("resetFileInputHandler", function(x) {      
                            var id = "#" + x + "_progress";
                            var idBar = id + " .bar";
                            $(id).css("visibility", "hidden");
                            $(idBar).css("width", "0%");
                            });
                            '),
                tags$hr(),
                helpPopup("Refresh the page", "Please refresh this page before upload/load another project to have a better performance.",
                          placement = "right", trigger = "click")
              ),
              mainPanel(
                helpText('Name of Uploaded Project: '),
                textOutput('path'),
                helpText(" "),
                #uiOutput("video"),
                helpText("The tabs on the left contain reports of various aspects of the project. Be aware that
                         the app utilizes 'lazy' loads, meaning that it only uses the necessary pieces of information to produce the image at hand.
                         So when you arrive at each tab for the first time there will be some small red error messages that are residuals of the lazy loading.
                         Give the app a few seconds (maybe a minute for the heatmaps) and the errors will go away and the reports generated.  Thank you.")
                ))
    ),
    
    tabItem(tabName = "design",
            tags$head(
              tags$style(type="text/css", "tfoot {display: table-header-group}") # Move filters to top of datatable
              #tags$style("table {background-color: grey !important;}",type="text/css")
            ),
            downloadButton('downloadDesign', "Download Table"),
            div(style = "height: 85vh; overflow: auto", dataTableOutput('designDataTable'))
    ),
    
    tabItem(tabName = "summary",
            sidebarLayout(
              sidebarPanel(
                uiOutput('summaryName'),
                uiOutput('respVar'),   
                selectInput("digits", 
                            label = "Number of decimal places:",
                            choices = c(0,1,2,3,4),
                            selected = 1)
              ),
              
              mainPanel(
                div(style = "display:inline-block", uiOutput("summary0Text")),
                div(style = "display:inline-block", helpPopup("Summary Table 1", "Given the summary variable selected by the user, statistics will be provided for 
                                                              each category of the users BY variable selection.  In longitudinal settings, this table will display 
                                                              a warning if the summary variable is changing over time.  Refer to Table 2.", placement = "right",
                                                              trigger = "click")),
                br(),
                downloadButton('downloadSummary0', 'Download Table'),
                tableOutput('summary0'),
                div(style = "display:inline-block", uiOutput("summary1Text")),
                div(style = "display:inline-block", helpPopup("Summary Table 2", "Given the summary variable selected by the user, statistics will be provided for 
                                                              each category of the users BY variable selection and for each time point in longitudinal settings.", placement = "right",
                                                              trigger = "click")),
                br(),
                downloadButton('downloadSummary1', 'Download Table'),
                tableOutput('summary1'),
                div(style = "display:inline-block", uiOutput("summary2Text")),
                div(style = "display:inline-block", helpPopup("Summary Table 3", "For longitudinal data, this table provides the counts of subjects categorized by 
                                                              their observed longitudinal profiles and by the specified BY variable.", placement = "right",
                                                              trigger = "click")),
                br(),
                downloadButton('downloadSummary2', 'Download Table'),
                tableOutput('summary2'),
                tableOutput('summary_tab')
                #textOutput('summary3Text')
              )
            )
          ),
    
    tabItem(tabName = "pvca",
            helpText("Data Filtering and Batch Assesment Summary:"),
            textOutput("PVCAtext"), #reactive text with the new content
            imageOutput('PVCA1'),  
            imageOutput('PVCA2')
    ),
    
    tabItem(tabName = "fastqc",
            sidebarLayout(
              sidebarPanel(
                uiOutput('qc_select'),
                uiOutput('qc_select2'),
                width = 3
              ),
              mainPanel(style="height: 94vh",
                        div(style="overflow: auto", uiOutput("fqc")),
                        #uiOutput("fqc"),
                        plotOutput("ex"),
                        plotOutput("ex2"),
                        plotOutput("ex3"),
                        width = 9
              )
            )
    ),
    
    tabItem(tabName = "qcupload",
            sidebarLayout(
              sidebarPanel(
                uiOutput("dropdown"),
                uiOutput("dropdowncolumns"),
                uiOutput("coordflip"),
                br(),
                div(style = "display:inline-block", checkboxInput("graphicsqc", strong("Graphing options", style = "color:blue"), TRUE)),
                div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                conditionalPanel(condition = "input.graphicsqc",
                                 sliderInput('qcplotsizew', "Grid width", min = 500, max = 2000, value = 800, step = 25),
                                 sliderInput('qcplotsizeh', "Grid height", min = 500, max = 2000, value = 600, step = 25),
                                 numericInput("axis_text_sizeqc", "Axis text size:", min = 10, value = 10, step = 1),
                                 numericInput("axis_label_sizeqc", "Axis label size:", min = 12, value = 12, step = 1),
                                 numericInput('PlotResolutionqc', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))),
              mainPanel(
                downloadButton('downloadQC', 'Download Figure'),
                uiOutput("QCbox")
              )
            )
    ),
   
    tabItem(tabName = "probeheatmap",
                     sidebarLayout(         
                       sidebarPanel(
                         conditionalPanel(condition="input.start2=='Probe Level Heat Maps'",
                                          div(style = "display:inline-block", actionButton("go", "Plot")),
                                          div(style = "display:inline-block", infoPopup("Plot", 'The heatmap will not update until the "Plot" button is clicked.
                                                                                        This allows the user to make multiple adjustments at once, without having to wait for 
                                                                                        each individual adjustment to update on the heatmap.', placement = "right", trigger = "click")),
                                          uiOutput('test'),
                                          br(),
                                          div(style = "display:inline-block", checkboxInput("uploadprobes", strong("Upload probes:", style = "color:blue"), FALSE)),
                                          div(style = "display:inline-block", infoPopup("Upload Probes", "Allows the user to provide their own list of genes (CSV) to plot. The CSV file should 
                                                                                        contain a single column named 'PROBE_ID' or 'SYMBOL', depending on whether the list provided is the PROBE 
                                                                                        ID's or gene symbols.", placement = "right", trigger = "click")),
                                          conditionalPanel(condition="input.uploadprobes",
                                                           fileInput('probe_select', '', multiple = FALSE,
                                                                     accept=c(".csv")),
                                                           div(style = "display:inline-block",checkboxInput("row_cluster", strong("Row cluster"), FALSE)),
                                                           div(style = "display:inline-block", infoPopup("Row cluster", 'If unchecked, the order of the probes will be exactly the same as the order
                                                                                                                         of the probes you upload. If checked, the probes will be row clustered'))),
                                          numericInput("setcutoff","Max value on color key:",2),
                                          br(),
                                          div(style = "display:inline-block", checkboxInput("subsetProbe",strong("Subset column options", style = "color:blue"),FALSE)),
                                          div(style = "display:inline-block", helpPopup("Subset column options", "These options allow the user to plot a subset of the samples within the expression data set.")),
                                          conditionalPanel(condition="input.subsetProbe",
                                                           uiOutput('subsetProbeVariable'),
                                                           uiOutput('subsetProbeValue'),
                                                           uiOutput('subsetProbeVariable2'),
                                                           uiOutput('subsetProbeValue2')),
                                          br(),
                                          div(style = "display:inline-block", checkboxInput('ordercolumns1', strong("Ordering column options", style = "color:blue"), FALSE)),
                                          div(style = "display:inline-block", helpPopup("Ordering column options", "These options allow the user to specify column ordering by column variable.")),
                                          conditionalPanel(condition = "input.ordercolumns1",
                                                           uiOutput('TopTierProbe'),
                                                           uiOutput('MidTierProbe'),
                                                           uiOutput('LowTierProbe'),
                                                           uiOutput('group_label_probe')),
                                          br(),
                                          div(style = "display:inline-block", checkboxInput('clusteroptions1', strong("Clustering column options", style = "color:blue"), FALSE)),
                                          div(style = "display:inline-block", helpPopup("Cluster column options", "These options allow the user to cluster the samples and label the groups they are clustered in.")),
                                          conditionalPanel(condition = "input.clusteroptions1",
                                                           checkboxInput("ColClustProbe","Cluster samples (columns)",FALSE),
                                                           #uiOutput('ColumnClusterProbe'),
                                                           checkboxInput('ClusterChoice4', "Show cluster groups", FALSE),
                                                           conditionalPanel(condition = "input.ClusterChoice4",
                                                                            div(style = "display:inline-block", uiOutput('ClusterCuts4')),
                                                                            div(style = "display:inline-block", infoPopup("Cluster Cuts", "Optimal number calculated using Dunn's Index.",
                                                                                                                        placement = "right", trigger = "click")))),
                                          br(),
                                          div(style = "display:inline-block", checkboxInput("graphics1", strong("Graphing options", style = "color:blue"), FALSE)),
                                          div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                          conditionalPanel(condition = "input.graphics1",
                                                           sliderInput('HeatMapSize1', "Grid width", min = 500, max = 10000, value = 750, step = 50),
                                                           sliderInput('HeatMapSize2', "Grid height", min = 500, max = 10000, value = 650, step = 50),
                                                           numericInput('FontSize', "Font size:", min = 10, value = 10, step = 2),
                                                           numericInput('LegendSize', "Legend size:", min = 1, max = 5, value = 1, step = .50),
                                                           numericInput('TreeHeight', "Tree height:", min = 0, value = 50, step = 5),
                                                           numericInput('PlotResolution2', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))),
                         
                         conditionalPanel(condition = "input.start2 == 'Cluster Number Diagnostics'"),
                         
                         conditionalPanel(condition = "input.start2=='Cluster Association Analysis'",
                                          uiOutput('ResponderStatus4'),
                                          uiOutput("clusternumber")
                         )
                                          )
                       ,
                       mainPanel(style = "overflow: auto; height: 92vh",
                         tabsetPanel(id="start2",
                                     tabPanel('Probe Level Heat Maps',
                                              div(style = "display:inline-block", strong("Help Message")),
                                              div(style = "display:inline-block", infoPopup("Probe Level Heat Maps", 'The heat map below is constructed on individual samples for a number of scenarios, 
                                                                                            specifically for baseline samples only (cross sectional) or all samples at all time points.  For baseline samples, 
                                                                                            heatmaps can be generated based on normalizing the expression data to the median, or to a control group if applicable.  
                                                                                            When all samples at all time points are to be plotted, heatmaps can be generated by normalizing to the median, 
                                                                                            a control group, or each subjects own baseline value. Only samples that have a corresponding baseline sample are 
                                                                                            included in the map. The initial graph of the heat map may not be very appealing depending on the number of samples. 
                                                                                            The inputs on the left have a wide variety of user options that range from the type of normalized data, clustering rows 
                                                                                            and columns, subsetting samples and probes, etc. These options are consistent across all the unsupervised analysis plots.  
                                                                                            One addition, unique to the probe level heat map is the "Max value on color key" option that specifies the coloring legend 
                                                                                            for the heatmap index. The default is set to +/- 2. If the user would like the index to have the most extreme red and blue 
                                                                                            to be set to a value of +/- 4 then user simply needs to enter 4. The user can also enter 0 which will make the limits based 
                                                                                            on the max and min of the entire expression file.',
                                                                                            placement = "bottom", trigger = "click")),
                                              helpText(""),
                                              textOutput('OptimalNumber'),
                                              downloadButton('downloadHeatmap', 'Download Data'),
                                              downloadButton("downloadModPlot4", "Download Figure"),
                                              helpText(""),
                                              div(style = "display:inline-block", checkboxInput("modulemeans", strong("Download module scores:", style = "color:blue"), FALSE)),
                                              div(style = "display:inline-block", infoPopup("Download module scores", "The downloaded data is calculated by averaging all the probes within each module.",
                                                                                            placement = "right", trigger = "click")),
                                              #div(style="height: 94vh", plotOutput('heatmap')),
                                              plotOutput('heatmap'),
                                              plotOutput('heatmap_select')),
                                     
                                     tabPanel('Cluster Number Diagnostics',
                                              downloadButton('downloadClusterPlot3', 'Download Figure'),
                                              plotOutput('clusterplot3')),
                                     
                                     tabPanel('Cluster Association Analysis',
                                              tableOutput('cluster_output4'),
                                              helpText("The Chi-square test below is used to test for association between the cluster groups and user specified variable. 
                                                       The Chi-square test is ideal when expected cell counts are large (expected values greater than 5). 
                                                       Because the test only makes sense for 2x2 contingency tables and above, test statistics are not shown 
                                                       when the data is subsetted on only one value."),
                                              helpText(""),
                                              helpText("Fisher's exact test is also shown. This test is ideal when cell counts are small. Like the Chi-square test, 
                                                       when the data is subsetted on only one value, test statistics are not shown."),
                                              helpText(""),
                                              textOutput('explanation3'),
                                              tableOutput('cluster_tab3'),
                                              textOutput('chisquare_test4'),
                                              textOutput('fisher4'))
                                              )
                         )
                       )
                         ),
    
    tabItem(tabName = "modulemap1",
            
            sidebarLayout(         
              sidebarPanel(
                conditionalPanel(condition="input.start1=='Individual Baseline Module Map'",
                                 uiOutput('modselTF3'),
                                 div(style = "display:inline-block",checkboxInput(inputId = "rowselect3", label = strong("Upload modules:", style = "color:blue"), value = FALSE)),
                                 div(style = "display:inline-block",helpPopup("Upload Modules", 'Allows the user to provide their own list of Modules (CSV) to plot. The CSV file should contain a single column named "Module". The modules in the list should begin with a capital "M" and a decimal point should separate the module number (e.g. M1.1).',
                                                                              placement = "right", trigger = "click")),
                                 br(),
                                 conditionalPanel(condition = "input.rowselect3 == true",
                                                  fileInput('modsel3', '', accept = ".csv")),
                                 div(style = "display:inline-block", checkboxInput("MMorder_3","Order Modules",FALSE)),
                                 div(style = "display:inline-block", infoPopup("Order Modules", 'The modules will be ordered by module number unless you select to upload a list of modules to plot from below (the order of modules will be exactly the same as the order of the modules you upload).',
                                                                               placement = "right", trigger = "click")),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("subsetMod",strong("Subset column options", style = "color:blue"),FALSE)),
                                 div(style = "display:inline-block", helpPopup("Subset column options", "These options allow the user to plot a subset of the samples within the expression data set.")),
                                 conditionalPanel(condition="input.subsetMod",
                                                  uiOutput('subsetModVariable'),
                                                  uiOutput('subsetModValue'),
                                                  uiOutput('subsetModVariable2'),
                                                  uiOutput('subsetModValue2')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('ordercolumns2', strong("Ordering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Ordering column options", "These options allow the user to specify column ordering by column variable.")),
                                 conditionalPanel(condition = "input.ordercolumns2",
                                                  uiOutput('TopTier'),
                                                  uiOutput('MidTier'),
                                                  uiOutput('LowTier'),
                                                  uiOutput('group_label_mod2')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('clusteroptions2', strong("Clustering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Cluster column options", "These options allow the user to cluster the samples and label the groups they are clustered in.")),
                                 conditionalPanel(condition = "input.clusteroptions2",
                                                  checkboxInput("ColClust","Cluster samples (columns)",FALSE),
                                                  checkboxInput('ClusterChoice2', "Show cluster groups", FALSE),
                                                  conditionalPanel(condition = "input.ClusterChoice2",
                                                                   div(style = "display:inline-block", uiOutput('ClusterCuts2')),
                                                                   div(style = "display:inline-block", infoPopup("Cluster Cuts", "Optimal number calculated using Dunn's Index.",
                                                                                                                 placement = "right", trigger = "click"))))
                                ,
                                br(),
                                div(style = "display:inline-block", checkboxInput("graphics2", strong("Graphing options", style = "color:blue"), FALSE)),
                                div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                conditionalPanel(condition = "input.graphics2",
                                                  sliderInput('modmap2size', "Grid width", min=500, max=10000,value=900,step=50),
                                                  sliderInput('modmap2size1', "Grid height", min = 500, max = 10000, value = 850, step = 50),
                                                  numericInput('modmap2radius',"Circle size (.25 to 4):",min=.25,max=4,value=1,step=.25),
                                                  #numericInput('rowfontsize', "Row Font size:", min = 1, max = 3, value = 1, step = .25),
                                                  #numericInput('colfontsize', "Column Font size:", min = 1, max = 3, value = 1, step = .25)
                                                  numericInput('Fontsize2', "Font size:", min = 10, value = 10, step = 2),
                                                  numericInput('PlotResolution', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                                 
                ),
                
                conditionalPanel(condition = "input.start1 == 'Cluster Number Diagnostics'"),
                
                conditionalPanel(condition = "input.start1=='Cluster Association Analysis'",
                                 uiOutput('ResponderStatus'),
                                 uiOutput("clusternumber1")
                                 
                )),
                
                mainPanel(
                  tabsetPanel(id="start1",
                              
                              
                              tabPanel('Individual Baseline Module Map',
                                       div(style = "display:inline-block", strong("Help Message")),
                                       div(style = "display:inline-block", infoPopup("Individual Baseline Module Map", "The module map below is constructed on individual samples at the baseline time point ONLY. Healthy Control samples are used to determine an upper and lower threshold (mean healthy controls +/- 2 sd). The module proportion for each sample is then calculated based on the percentage of probes within a module that are above or below this threshold.  If there are no healthy controls in the study this graph will be blank.
                                                 
                                                 The initial graph of the modules may not be very appealing depending on the number of samples and modules. The inputs on the left are designed to let the user modify the map to their liking. The first 3 inputs allow the user to order the samples how they would like.
                                                 The default setting is to order the samples by responder status and then by donor id. The ordering variables are labeled across the top and included in the legend.
                                                 The next box allows for additional labeling of the sample across the top but do not contribute to the ordering of the samples. 
                                                 The remaining options are either graphical parameters or subsetting paramaters. The user can specify to look at all the Baylor modules at once, only the first 6 rounds, or only annotated ones.
                                                 If there are too many samples the user can subset the data (Subset the module map checkbox) by selecting up to two variables to subset on. This allows the user to only look at certain subgroups like just the male, high responders.
                                                 Use the graphing inputs to make the heatmap pretty. It is suggested to start with the size and aspect ratio and only edit the circle sizes if needed after.",
                                                 placement = "bottom", trigger = "click")),
                                       helpText(""),
                                       textOutput('OptimalNumber1'),
                                       downloadButton('downloadModMap2', 'Download Data'),
                                       downloadButton("downloadModPlot2", "Download Figure"),
                                       plotOutput('modmap2'),
                                       plotOutput('modmap2plus')),
    
                              tabPanel('Cluster Number Diagnostics',
                                       downloadButton('downloadClusterPlot', 'Download Figure'),
                                       plotOutput('clusterplot')),
                              
                              tabPanel('Cluster Association Analysis',
                                       tableOutput('cluster_output'),
                                       helpText("The Chi-square test below is used to test for association between the cluster groups and user specified variable. 
                                                The Chi-square test is ideal when expected cell counts are large (expected values greater than 5). 
                                                Because the test only makes sense for 2x2 contingency tables and above, test statistics are not shown 
                                                when the data is subsetted on only one value."),
                                       helpText(""),
                                       helpText("Fisher's exact test is also shown. This test is ideal when cell counts are small. Like the Chi-square test, 
                                                when the data is subsetted on only one value, test statistics are not shown."),
                                       helpText(""),
                                       textOutput('explanation'),
                                       tableOutput('cluster_tab'),
                                       textOutput('chisquare_test'),
                                       textOutput('fisher')))))),
    
    tabItem(tabName = "modulemap2",
            sidebarLayout(
              sidebarPanel(
                conditionalPanel(condition="input.start11=='Individual Longitudinal Module Map'",
                                 uiOutput("BaseOrHealthy1"),
                                 uiOutput('modselTF2'),
                                 div(style = "display:inline-block", checkboxInput(inputId = "rowselect2", label = strong("Upload modules:", style = "color:blue"), value = FALSE)),
                                 div(style = "display:inline-block", helpPopup("Upload Modules", 'Allows the user to provide their own list of Modules (CSV) to plot. The CSV file should contain a single column named "Module". The modules in the list should begin with a capital "M" and a decimal point should separate the module number (e.g. M1.1).',
                                                                               placement = "right", trigger = "click")),
                                 br(),
                                 conditionalPanel(condition = "input.rowselect2 == true",
                                                  fileInput('modsel2', '', accept = ".csv")),
                                 div(style = "display:inline-block", checkboxInput("MMorder_2","Order Modules",FALSE)),
                                 div(style = "display:inline-block", infoPopup("Order Modules", 'The modules will be ordered by module number unless you select to upload a list of modules to plot from below (the order of modules will be exactly the same as the order of the modules you upload).',
                                                                               placement = "right", trigger = "click")),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("subsetMod3",strong("Subset column options", style = "color:blue"),FALSE)),
                                 div(style = "display:inline-block", helpPopup("Subset column options", "These options allow the user to plot a subset of the samples within the expression data set.")),
                                 conditionalPanel(condition="input.subsetMod3",
                                                  uiOutput('subsetMod3Variable'),
                                                  uiOutput('subsetMod3Value'),
                                                  uiOutput('subsetMod3Variable2'),
                                                  uiOutput('subsetMod3Value2')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('ordercolumns3', strong("Ordering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Ordering column options", "These options allow the user to specify column ordering by column variable.")),
                                 conditionalPanel(condition = "input.ordercolumns3",
                                                  uiOutput('TopTier3'),
                                                  uiOutput('MidTier3'),
                                                  uiOutput('LowTier3'),
                                                  uiOutput('group_label_mod3')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('clusteroptions3', strong("Clustering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Cluster column options", "These options allow the user to cluster the samples and label the groups they are clustered in.")),
                                 conditionalPanel(condition = "input.clusteroptions3",
                                                  checkboxInput("ColClust3","Cluster samples (columns)",FALSE),
                                                  #uiOutput('ColumnClusterProbe'),
                                                  checkboxInput('ClusterChoice3', "Show cluster groups", FALSE),
                                                  conditionalPanel(condition = "input.ClusterChoice3",
                                                                   div(style = "display:inline-block", uiOutput('ClusterCuts3')),
                                                                   div(style = "display:inline-block", infoPopup("Cluster Cuts", "Optimal number calculated using Dunn's Index.",
                                                                                                                 placement = "right", trigger = "click")))),
                                 
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("graphics3", strong("Graphing options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphics3",
                                                  sliderInput('modmap3size', "Grid width", min=500, max=10000,value=900,step=50),
                                                  sliderInput('modmap3size1', "Grid height", min = 500, max = 10000, value = 850, step = 50),
                                                  numericInput('modmap3radius',"Circle size (.25 to 4):",min=.25,max=4,value=1,step=.25),
                                                  numericInput('Fontsize', "Font size:", min = 10, value = 10, step = 2),
                                                  numericInput('PlotResolution1', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))),
                
                conditionalPanel(condition = "input.start11 == 'Cluster Number Diagnostics'"),
                
                conditionalPanel(condition = "input.start11=='Cluster Association Analysis'",
                                 uiOutput('ResponderStatus3'),
                                 uiOutput("clusternumber2")
                                 
                )),
                
                mainPanel(
                  tabsetPanel(id="start11",
                              
                              tabPanel('Individual Longitudinal Module Map',
                                       div(style = "display:inline-block", strong("Help Message")),
                                       div(style = "display:inline-block", infoPopup("Individual Longitudinal Module Map", "The module map below is constructed on individual samples for all time points. Healthy Control samples are used to determine an upper and lower threshold (mean healthy controls +/- 2 sd). The module proportion for each sample is then calculated based on the percentage of probes within a module that are above or below this threshold.  If there are no healthy controls in the study, the baseline samples will be treated as a healthy control and will serve as the threshold.
                                           
                                           The initial graph of the modules may not be very appealing depending on the number of samples and modules. The inputs on the left are designed to let the user modify the map to their liking. The first 3 inputs allow the user to order the samples how they would like.
                                           The default setting is to order the samples by responder status and then by donor id. The ordering variables are labeled across the top and included in the legend.
                                           The next box allows for additional labeling of the sample across the top but do not contribute to the ordering of the samples.
                                           The remaining options are either graphical parameters or subsetting paramaters. The user can specify to look at all the Baylor modules at once, only the first 6 rounds, or only annotated ones.
                                           If there are too many samples the user can subset the data (Subset the module map checkbox) by selecting up to two variables two subset on. This allows the user to only look at certain subgroups like just the male, high responders.
                                           Use the graphing inputs to make the heatmap pretty. It is suggested to start with the size and aspect ratio and only edit the circle sizes if needed after.",
                                                 placement = "bottom", trigger = "click")),
                                       helpText(""),
                                       textOutput('OptimalNumber2'),
                                       downloadButton('downloadModMap3', 'Download Data'),
                                       downloadButton("downloadModPlot3", "Download Figure"),
                                       plotOutput('modmap3'),
                                       plotOutput('modmap3plus')),
                              
                              tabPanel('Cluster Number Diagnostics',
                                       downloadButton('downloadClusterPlot2', 'Download Figure'),
                                       plotOutput('clusterplot2')),
                              
                              tabPanel('Cluster Association Analysis',
                                       tableOutput('cluster_output3'),
                                       helpText("The Chi-square test below is used to test for association between the cluster groups and user specified variable. 
                                          The Chi-square test is ideal when expected cell counts are large (expected values greater than 5). 
                                          Because the test only makes sense for 2x2 contingency tables and above, test statistics are not shown 
                                          when the data is subsetted on only one value."),
                                       helpText(""),
                                       helpText("Fisher's exact test is also shown. This test is ideal when cell counts are small. Like the Chi-square test, 
                                          when the data is subsetted on only one value, test statistics are not shown."),
                                       helpText(""),
                                       textOutput('explanation2'),
                                       tableOutput('cluster_tab2'),
                                       textOutput('chisquare_test3'),
                                       textOutput('fisher3')))))),
    
    tabItem(tabName = "othermodmap1",
            
            sidebarLayout(         
              sidebarPanel(
                conditionalPanel(condition="input.start1Other=='Individual Baseline Module Map'",
                                 #uiOutput('modselTF3'),
                                 div(style = "display:inline-block",checkboxInput(inputId = "rowselect3Other", label = strong("Upload modules:", style = "color:blue"), value = FALSE)),
                                 div(style = "display:inline-block",helpPopup("Upload Modules", 'Allows the user to provide their own list of Modules (CSV) to plot. The CSV file should contain a single column named "Module". The modules in the list should begin with a capital "M" and a decimal point should separate the module number (e.g. M1.1).',
                                                                              placement = "right", trigger = "click")),
                                 br(),
                                 conditionalPanel(condition = "input.rowselect3Other == true",
                                                  fileInput('modsel3Other', '', accept = ".csv")),
                                 div(style = "display:inline-block", checkboxInput("MMorder_3Other","Order Modules",FALSE)),
                                 div(style = "display:inline-block", infoPopup("Order Modules", 'The modules will be ordered by module number unless you select to upload a list of modules to plot from below (the order of modules will be exactly the same as the order of the modules you upload).',
                                                                               placement = "right", trigger = "click")),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("subsetModOther",strong("Subset column options", style = "color:blue"),FALSE)),
                                 div(style = "display:inline-block", helpPopup("Subset column options", "These options allow the user to plot a subset of the samples within the expression data set.")),
                                 conditionalPanel(condition="input.subsetModOther",
                                                  uiOutput('subsetModVariableOther'),
                                                  uiOutput('subsetModValueOther'),
                                                  uiOutput('subsetModVariable2Other'),
                                                  uiOutput('subsetModValue2Other')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('ordercolumns2Other', strong("Ordering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Ordering column options", "These options allow the user to specify column ordering by column variable.")),
                                 conditionalPanel(condition = "input.ordercolumns2Other",
                                                  uiOutput('TopTierOther'),
                                                  uiOutput('MidTierOther'),
                                                  uiOutput('LowTierOther'),
                                                  uiOutput('group_label_mod2Other')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('clusteroptions2Other', strong("Clustering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Cluster column options", "These options allow the user to cluster the samples and label the groups they are clustered in.")),
                                 conditionalPanel(condition = "input.clusteroptions2Other",
                                                  checkboxInput("ColClustOther","Cluster samples (columns)",FALSE),
                                                  #uiOutput('ColumnClusterProbe'),
                                                  checkboxInput('ClusterChoice2Other', "Show cluster groups", FALSE),
                                                  conditionalPanel(condition = "input.ClusterChoice2Other",
                                                                   div(style = "display:inline-block", uiOutput('ClusterCuts2Other')),
                                                                   div(style = "display:inline-block", infoPopup("Cluster Cuts", "Optimal number calculated using Dunn's Index.",
                                                                                                                 placement = "right", trigger = "click"))))
                                 ,
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("graphics2Other", strong("Graphing options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphics2Other",
                                                  sliderInput('modmap2sizeOther', "Grid width", min=500, max=10000,value=900,step=50),
                                                  sliderInput('modmap2size1Other', "Grid height", min = 500, max = 10000, value = 850, step = 50),
                                                  numericInput('modmap2radiusOther',"Circle size (.25 to 4):",min=.25,max=4,value=1,step=.25),
                                                  #numericInput('rowfontsize', "Row Font size:", min = 1, max = 3, value = 1, step = .25),
                                                  #numericInput('colfontsize', "Column Font size:", min = 1, max = 3, value = 1, step = .25)
                                                  numericInput('Fontsize2Other', "Font size:", min = 10, value = 10, step = 2),
                                                  numericInput('PlotResolutionOther', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                                 
                ),
                
                conditionalPanel(condition = "input.start1Other == 'Cluster Number Diagnostics'"),
                
                conditionalPanel(condition = "input.start1Other=='Cluster Association Analysis'",
                                 uiOutput('ResponderStatusOther'),
                                 uiOutput("clusternumber1Other")
                                 
                )),
              
              mainPanel(
                tabsetPanel(id="start1Other",
                            
                            
                            tabPanel('Individual Baseline Module Map',
                                     div(style = "display:inline-block", strong("Help Message")),
                                     div(style = "display:inline-block", infoPopup("Individual Baseline Module Map", "The module map below is constructed on individual samples at the baseline time point ONLY. Healthy Control samples are used to determine an upper and lower threshold (mean healthy controls +/- 2 sd). The module proportion for each sample is then calculated based on the percentage of probes within a module that are above or below this threshold.  If there are no healthy controls in the study this graph will be blank.
                                                                                   
                                                                                   The initial graph of the modules may not be very appealing depending on the number of samples and modules. The inputs on the left are designed to let the user modify the map to their liking. The first 3 inputs allow the user to order the samples how they would like.
                                                                                   The default setting is to order the samples by responder status and then by donor id. The ordering variables are labeled across the top and included in the legend.
                                                                                   The next box allows for additional labeling of the sample across the top but do not contribute to the ordering of the samples. 
                                                                                   The remaining options are either graphical parameters or subsetting paramaters. The user can specify to look at all the Baylor modules at once, only the first 6 rounds, or only annotated ones.
                                                                                   If there are too many samples the user can subset the data (Subset the module map checkbox) by selecting up to two variables to subset on. This allows the user to only look at certain subgroups like just the male, high responders.
                                                                                   Use the graphing inputs to make the heatmap pretty. It is suggested to start with the size and aspect ratio and only edit the circle sizes if needed after.",
                                                                                   placement = "bottom", trigger = "click")),
                                     helpText(""),
                                     textOutput('OptimalNumber1Other'),
                                     downloadButton('downloadModMap2Other', 'Download Data'),
                                     downloadButton("downloadModPlot2Other", "Download Figure"),
                                     plotOutput('modmap2Other'),
                                     plotOutput('modmap2plusOther')),
                            #verbatimTextOutput("info")),
                            
                            tabPanel('Cluster Number Diagnostics',
                                     downloadButton('downloadClusterPlotOther', 'Download Figure'),
                                     plotOutput('clusterplotOther')),
                            
                            tabPanel('Cluster Association Analysis',
                                     tableOutput('cluster_outputOther'),
                                     helpText("The Chi-square test below is used to test for association between the cluster groups and user specified variable. 
                                              The Chi-square test is ideal when expected cell counts are large (expected values greater than 5). 
                                              Because the test only makes sense for 2x2 contingency tables and above, test statistics are not shown 
                                              when the data is subsetted on only one value."),
                                     helpText(""),
                                     helpText("Fisher's exact test is also shown. This test is ideal when cell counts are small. Like the Chi-square test, 
                                              when the data is subsetted on only one value, test statistics are not shown."),
                                     helpText(""),
                                     textOutput('explanationOther'),
                                     tableOutput('cluster_tabOther'),
                                     textOutput('chisquare_testOther'),
                                     textOutput('fisherOther')))))),
    
    tabItem(tabName = "othermodmap2",
            sidebarLayout(
              sidebarPanel(
                conditionalPanel(condition="input.start11Other=='Individual Longitudinal Module Map'",
                                 uiOutput("BaseOrHealthy1Other"),
                                 #uiOutput('modselTF2'),
                                 div(style = "display:inline-block", checkboxInput(inputId = "rowselect2Other", label = strong("Upload modules:", style = "color:blue"), value = FALSE)),
                                 div(style = "display:inline-block", helpPopup("Upload Modules", 'Allows the user to provide their own list of Modules (CSV) to plot. The CSV file should contain a single column named "Module". The modules in the list should begin with a capital "M" and a decimal point should separate the module number (e.g. M1.1).',
                                                                               placement = "right", trigger = "click")),
                                 br(),
                                 conditionalPanel(condition = "input.rowselect2Other == true",
                                                  fileInput('modsel2Other', '', accept = ".csv")),
                                 div(style = "display:inline-block", checkboxInput("MMorder_2Other","Order Modules",FALSE)),
                                 div(style = "display:inline-block", infoPopup("Order Modules", 'The modules will be ordered by module number unless you select to upload a list of modules to plot from below (the order of modules will be exactly the same as the order of the modules you upload).',
                                                                               placement = "right", trigger = "click")),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("subsetMod3Other",strong("Subset column options", style = "color:blue"),FALSE)),
                                 div(style = "display:inline-block", helpPopup("Subset column options", "These options allow the user to plot a subset of the samples within the expression data set.")),
                                 conditionalPanel(condition="input.subsetMod3Other",
                                                  uiOutput('subsetMod3VariableOther'),
                                                  uiOutput('subsetMod3ValueOther'),
                                                  uiOutput('subsetMod3Variable2Other'),
                                                  uiOutput('subsetMod3Value2Other')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('ordercolumns3Other', strong("Ordering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Ordering column options", "These options allow the user to specify column ordering by column variable.")),
                                 conditionalPanel(condition = "input.ordercolumns3Other",
                                                  uiOutput('TopTier3Other'),
                                                  uiOutput('MidTier3Other'),
                                                  uiOutput('LowTier3Other'),
                                                  uiOutput('group_label_mod3Other')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('clusteroptions3Other', strong("Clustering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Cluster column options", "These options allow the user to cluster the samples and label the groups they are clustered in.")),
                                 conditionalPanel(condition = "input.clusteroptions3Other",
                                                  checkboxInput("ColClust3Other","Cluster samples (columns)",FALSE),
                                                  #uiOutput('ColumnClusterProbe'),
                                                  checkboxInput('ClusterChoice3Other', "Show cluster groups", FALSE),
                                                  conditionalPanel(condition = "input.ClusterChoice3Other",
                                                                   div(style = "display:inline-block", uiOutput('ClusterCuts3Other')),
                                                                   div(style = "display:inline-block", infoPopup("Cluster Cuts", "Optimal number calculated using Dunn's Index.",
                                                                                                                 placement = "right", trigger = "click")))),
                                 
                                 #checkboxInput(inputId = "rowselect2", label = strong("Upload modules:", style = "color:blue"), value = FALSE),
                                 #conditionalPanel(condition = "input.rowselect2 == true",
                                 #fileInput('modsel2', '', accept = ".csv")
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("graphics3Other", strong("Graphing options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphics3Other",
                                                  sliderInput('modmap3sizeOther', "Grid width", min=500, max=10000,value=900,step=50),
                                                  sliderInput('modmap3size1Other', "Grid height", min = 500, max = 10000, value = 850, step = 50),
                                                  numericInput('modmap3radiusOther',"Circle size (.25 to 4):",min=.25,max=4,value=1,step=.25),
                                                  #numericInput('rowfontsize', "Row Font size::", min = 1, max = 3, value = 1, step = .25),
                                                  #numericInput('colfontsize', "Column Font size::", min = 1, max = 3, value = 1, step = .25)
                                                  numericInput('FontsizeOther', "Font size:", min = 10, value = 10, step = 2),
                                                  numericInput('PlotResolution1Other', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))),
                
                conditionalPanel(condition = "input.start11Other == 'Cluster Number Diagnostics'"),
                
                conditionalPanel(condition = "input.start11Other=='Cluster Association Analysis'",
                                 uiOutput('ResponderStatus3Other'),
                                 uiOutput("clusternumber2Other")
                                 
                )),
              
              mainPanel(
                tabsetPanel(id="start11Other",
                            
                            tabPanel('Individual Longitudinal Module Map',
                                     div(style = "display:inline-block", strong("Help Message")),
                                     div(style = "display:inline-block", infoPopup("Individual Longitudinal Module Map", "The module map below is constructed on individual samples for all time points. Healthy Control samples are used to determine an upper and lower threshold (mean healthy controls +/- 2 sd). The module proportion for each sample is then calculated based on the percentage of probes within a module that are above or below this threshold.  If there are no healthy controls in the study, the baseline samples will be treated as a healthy control and will serve as the threshold.
                                                                                   
                                                                                   The initial graph of the modules may not be very appealing depending on the number of samples and modules. The inputs on the left are designed to let the user modify the map to their liking. The first 3 inputs allow the user to order the samples how they would like.
                                                                                   The default setting is to order the samples by responder status and then by donor id. The ordering variables are labeled across the top and included in the legend.
                                                                                   The next box allows for additional labeling of the sample across the top but do not contribute to the ordering of the samples.
                                                                                   The remaining options are either graphical parameters or subsetting paramaters. The user can specify to look at all the Baylor modules at once, only the first 6 rounds, or only annotated ones.
                                                                                   If there are too many samples the user can subset the data (Subset the module map checkbox) by selecting up to two variables two subset on. This allows the user to only look at certain subgroups like just the male, high responders.
                                                                                   Use the graphing inputs to make the heatmap pretty. It is suggested to start with the size and aspect ratio and only edit the circle sizes if needed after.",
                                                                                   placement = "bottom", trigger = "click")),
                                     helpText(""),
                                     textOutput('OptimalNumber2Other'),
                                     downloadButton('downloadModMap3Other', 'Download Data'),
                                     downloadButton("downloadModPlot3Other", "Download Figure"),
                                     plotOutput('modmap3Other'),
                                     plotOutput('modmap3plusOther')),
                            
                            tabPanel('Cluster Number Diagnostics',
                                     downloadButton('downloadClusterPlot2Other', 'Download Figure'),
                                     plotOutput('clusterplot2Other')),
                            
                            tabPanel('Cluster Association Analysis',
                                     tableOutput('cluster_output3Other'),
                                     helpText("The Chi-square test below is used to test for association between the cluster groups and user specified variable. 
                                              The Chi-square test is ideal when expected cell counts are large (expected values greater than 5). 
                                              Because the test only makes sense for 2x2 contingency tables and above, test statistics are not shown 
                                              when the data is subsetted on only one value."),
                                     helpText(""),
                                     helpText("Fisher's exact test is also shown. This test is ideal when cell counts are small. Like the Chi-square test, 
                                              when the data is subsetted on only one value, test statistics are not shown."),
                                     helpText(""),
                                     textOutput('explanation2Other'),
                                     tableOutput('cluster_tab2Other'),
                                     textOutput('chisquare_test3Other'),
                                     textOutput('fisher3Other')))))),
    
    tabItem(tabName = "generalinfo",
            textOutput("readme"), helpText(" "), textOutput("MixedModelText")
    ),
    
    tabItem(tabName = "overview",
            sidebarLayout(         
              sidebarPanel(
                numericInput('alphalevel2',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025),                                                            
                selectInput("sigsign", "Fold change sign:", c("All (include 0)", "+", "-"), selected = "all (include 0)"),
                tags$hr() # add a horizontal line
              ),
              mainPanel(
                helpText("Cell value represents the number of significant probes under the selected significance level."),
                downloadButton('downloadSC', 'Download Table'),
                dataTableOutput("sigcomptable")
              ))
    ),
    
    tabItem(tabName = "genelistmaker",
            sidebarLayout(
              sidebarPanel(
                
                conditionalPanel(condition="input.start=='Significant Probe Level Heat Map'",
                                 div(style = "display:inline-block", actionButton("go2", "Plot")),
                                 div(style = "display:inline-block", infoPopup("Plot", 'The heatmap will not update until the "Plot" button is clicked.
                                                                               This allows the user to make multiple adjustments at once, without having to wait for 
                                                                               each individual adjustment to update on the heatmap.', placement = "right", trigger = "click"))),
                
                conditionalPanel(condition = "(input.start == 'DGE: Gene Lists')||(input.start == 'Significant Probe Level Heat Map')",
                                 uiOutput("Comparison1")
                ),
                
                conditionalPanel(condition = "input.start != 'DGE: Diagnostics'",
                                 div(style = "display:inline-block", checkboxInput("TSoptions", strong("Testing and fold change subsetting options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Testing and fold change subsetting options", "These options allow the user to filter on signifance threshold and fold change.")),
                                 conditionalPanel(condition = "input.TSoptions",
                                                  selectInput("correction_method1","Multiple testing correction:",
                                                              c("FDR","Bonferroni","Raw"),"Raw"),
                                                  numericInput('alphalevel1',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025),
                                                  checkboxInput(inputId = "showfc",label = strong("Filter on fold change"),value = FALSE),
                                                  conditionalPanel(condition = "input.showfc == true",
                                                                   numericInput(inputId = "fcval", label = "FC cut off:",
                                                                                min = 0, max = 10, value = 0, step = 0.25),
                                                                   selectInput("sign","Fold change sign:", c("+", "-", "Both"), c("Both")))
                                 )
                ),
                
                conditionalPanel(condition = "input.start == 'DGE: Gene Lists'",
                                 checkboxInput("merge", "Merge module information", FALSE),
 				 numericInput('modmaplmmres', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1),
                                 #downloadButton('downloadData', 'Download Table')
                                 #downloadButton('downloadAllComp', 'Download genelists for all comparisons')
 				              checkboxInput("DownloadOptions", strong("Download multiple genelists at once:", style = "color:blue"), FALSE),
 				              conditionalPanel(condition = "input.DownloadOptions",
 				                               uiOutput("selectComps"),
 				                               textInput("ziptext", label = "Filename Input", value = "All_Comparisons_Raw.05"),
 				                               downloadButton('downloadSelComp', 'Download selected genelists'))
                ),
                
                conditionalPanel(condition = "input.start == 'Cluster Number Diagnostics'"
                ),
                
                #uiOutput("DGE_Module_Analysis_Overview1"),
                conditionalPanel(condition = "input.start=='DGE: Module Analysis Overview'",
                                 uiOutput("modselTF"),
                                 #numericInput("LMMxcex","X Label Font size::",1),
                                 #numericInput("LMMycex","Y Label Fond Size:",1),
                                 div(style = "display:inline-block", checkboxInput("LMMmodorder","Order Modules",FALSE)),
                                 div(style = "display:inline-block", infoPopup("Order Modules", 'The modules will be ordered by module number unless you select to upload a list of modules to plot from below (the order of modules will be exactly the same as the order of the modules you upload).',
                                                                               placement = "right", trigger = "click")),
				 br(),
 				 div(style = "display:inline-block", checkboxInput("LMMdeleterows", "Delete all zero rows", FALSE)),
                                 div(style = "display:inline-block", infoPopup("Delete all zero rows", 'Rows with all zeros (all white rows) will be deleted from the heat map.',
                                                                               placement = "right", trigger = "click")),
                                 numericInput("percentage","Percentage cut off value for Heat index:",100),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput(inputId = "comparselect", label = strong("Comparison(s) selection:", style = "color:blue"), value = FALSE)),
                                 div(style = "display:inline-block", helpPopup("Comparison(s) selection", "This option allows the user to select specific comparisons of interest")),
                                 conditionalPanel(condition = "input.comparselect == true",
                                                  uiOutput("Comparison2")),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput(inputId = "rowselect", label = strong("Upload modules:", style = "color:blue"), value = FALSE)),
                                 div(style = "display:inline-block", helpPopup("Upload list of modules", 'Allows the user to provide their own list of Modules (CSV) to plot. The CSV file should contain a single column named "Module". The modules in the list should begin with a capital "M" and a decimal point should separate the module number (e.g. M1.1).')),
                                 conditionalPanel(condition = "input.rowselect == true",
                                                  fileInput('modsel', '', accept = ".csv")
                                                  ),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("graphics5", strong("Graphing options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphics5",
                                                  sliderInput('LMMmap_1', "Grid width", min = 500, max = 10000,value = 1000, step = 50),
                                                  sliderInput('LMMmap_2', "Grid height", min = 500, max = 10000, value = 850, step = 50),
                                                  numericInput('LMMradius',"Circle size (.25 to 4)", min = .25, max = 4,value = 1,step= .25),
                                                  numericInput('LMMFont', "Font size:", min = 10, value = 10, step = 2),
                                                  numericInput('DGE_PlotRes1', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                ),
                
                conditionalPanel(condition = "input.start == 'Venn Diagram'",
                                 uiOutput("vennComparison"),
                                 selectInput("UorI", "Intersection or union:", c("Intersection" = 1, "Union" = 2), selected = 1),
                                 uiOutput("include"),
                                 uiOutput("exclude")),
                conditionalPanel(condition="input.start=='Significant Probe Level Heat Map'",
                                 uiOutput('test1'),
                                 numericInput("setcutoff1","Select the maximum values on color key:",2),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("hm1rclu", strong("Row Cluster"), value = FALSE)),
                                 div(style = "display:inline-block", infoPopup("Row Cluster", "Based on the PC performance, it might not be recommended to cluster the rows if the number of rows exceed 7000",
                                                                               placement = "right", trigger = "click")),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("subsetProbe1",strong("Subset column options", style = "color:blue"),FALSE)),
                                 div(style = "display:inline-block", helpPopup("Subset column options", "These options allow the user to plot a subset of the samples within the expression data set.")),
                                 conditionalPanel(condition="input.subsetProbe1",
                                                  uiOutput('subsetProbeVariable1'),
                                                  uiOutput('subsetProbeValue1'),
                                                  uiOutput('subsetProbeVariable3'),
                                                  uiOutput('subsetProbeValue3')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('ordercolumns4', strong("Ordering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Ordering column options", "These options allow the user to specify column ordering by column variable.")),
                                 conditionalPanel(condition = "input.ordercolumns4",
                                                  uiOutput('TopTierProbe1'),
                                                  uiOutput('MidTierProbe1'),
                                                  uiOutput('LowTierProbe1'),
                                                  uiOutput('group_label_probe1')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('clusteroptions4', strong("Clustering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Cluster column options", "These options allow the user to cluster the samples and label the groups they are clustered in.")),
                                 conditionalPanel(condition = "input.clusteroptions4",
                                                  checkboxInput("ColClustProbe1","Cluster samples (columns)",FALSE),
                                                  #uiOutput('ColumnClusterProbe'),
                                                  checkboxInput('ClusterChoice5', "Show cluster groups", FALSE),
                                                  conditionalPanel(condition = "input.ClusterChoice5",
                                                                   div(style = "display:inline-block", uiOutput('ClusterCuts5')),
                                                                   div(style = "display:inline-block", infoPopup("Cluster Cuts", "Optimal number calculated using Dunn's Index.",
                                                                                                                 placement = "right", trigger = "click")))),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("graphics4", strong("Graphing options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphics4",
                                                  sliderInput('DGE_HeatMapSize1', "Grid width", min = 500, max = 10000, value = 750, step = 50),
                                                  sliderInput('DGE_HeatMapSize2', "Grid height", min = 500, max = 10000, value = 650, step = 50),
                                                  numericInput('DGE_FontSize', "Font size:", min = 10, value = 10, step = 2),
                                                  numericInput('DGE_LegendSize', "Legend size:", min = 1, max = 5, value = 1, step = .50),
                                                  numericInput('DGE_TreeHeight', "Tree height:", min = 0, value = 50, step = 5),
                                                  numericInput('DGE_PlotRes', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))),
                
                conditionalPanel(condition = "input.start=='Cluster Association Analysis'",
                                 uiOutput('ResponderStatus5'),
                                 uiOutput("clusternumber3")
                )
                
                                 ),
              
              mainPanel(
                uiOutput("mytabs")
                
              ))
        ),
    
      tabItem(tabName = "genesearch",
               sidebarLayout(
                 sidebarPanel(
                   uiOutput("specgene1"),
                   uiOutput("probeid1"),
                   uiOutput("genTable")
                 ),
                 mainPanel("DGE: Comparisons for specific Genes",
                           helpText("The processing speed depends on the number of gene symbols and your computer's performance. It may take up to 30 seconds to proceed."),
                           downloadButton('downloadCG','Download Table'),
                           dataTableOutput("specgenetable")
                 ))
        ),
    
    tabItem(tabName = "qusage",
            sidebarLayout(
              sidebarPanel(
                conditionalPanel(condition = "input.start3 == 'Significant Comparisons Overview'",
                                 #uiOutput("genesets"),
                                 numericInput("SigLevel", "Significance threshold:", min = 0, max = 1, step = .025, value = .05),
                                 selectInput("foldchange.q", strong("Fold change sign:"), choices = c("All (include 0)", "+","-"), selected = "All (include 0)")
                ),
                
                conditionalPanel(condition = "input.start3 == 'Fold Change Plot and Venn Diagram'",
                                 infoPopup("Comparisons", "The user may select up to five comparisons.", placement = "right",
                                           trigger = "click"),
                                 uiOutput("Comparisons1"),
                                 #div(style = "display:inline-block", uiOutput("Comparisons1")),
                                 #div(style = "display:inline-block", infoPopup("Comparisons", "The user may select up to five comparisons.", placement = "right",
                                                                               #trigger = "click")),
                                 uiOutput("PaloOrFirst1"),
                                 checkboxInput("only_annotated", strong("Plot only annotated modules"), FALSE),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("test_subset", strong("Testing and fold change subsetting options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Testing and fold change subsetting options", "These options allow the user to filter on signifance threshold and fold change.", placement = "bottom")),
                                 conditionalPanel(condition = "input.test_subset",
                                                  selectInput("Multi_testing1", "Multiple testing correction:", choices = c("FDR","Bonferroni","Raw"), selected = "Raw"),
                                                  #uiOutput("MultipleTesting1"),
                                                  numericInput("FilterOnPValues1", "Significance threshold:", min = 0, max = 1, step = .01, value = .05),
                                                  div(style = "display:inline-block",numericInput("FilterOnFoldchange1", "cut off:", min = 0, max = 1, step = .01, value = 0)),
                                                  div(style = "display:inline-block",helpPopup("FilterOnFoldchange1", "Help Message: Throw out modules with a fold change between +- the cutt off value",
                                                            placement = "right", trigger = "click"))),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("graphics6", strong("Graphing options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphics6",
                                                  div(style = "display:inline-block",uiOutput("ColorChoice")), 
                                                  div(style = "display:inline-block",helpPopup("Help Message", "The order of the colors corresponds to the order in the legend, not the order selected.",
                                                         placement = "right", trigger = "click")),
                                                  sliderInput("PlotWidth1", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                                  numericInput("PlotRes1", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                                 
                ),
                
                conditionalPanel(condition = "input.start3 == 'Individual Gene Set Plot and Table'",
                                 uiOutput("Module_Select"),
                                 uiOutput("Comparisons2"),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("test_subset2", strong("Testing and fold change subsetting options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Testing and fold change subsetting options", "These options allow the user to filter on signifance threshold and fold change.")),
                                 conditionalPanel(condition = "input.test_subset2",
                                                  selectInput("Multi_testing2", "Multiple testing correction:", choices = c("FDR","Bonferroni","Raw"), selected = "Raw"),
                                                  numericInput("FilterOnPValues2", "Significance threshold:", min = 0, max = 1, step = .01, value = 1),
                                                  div(style = "display:inline-block",numericInput("FilterOnFoldchange2", "cut off:", min = 0, max = 1, step = .01, value = 0)),
                                                  div(style = "display:inline-block",helpPopup("FilterOnFoldchange2","Help Message: Throw out modules with a fold change between +- the cut off value",
                                                            placement = "right", trigger = "click"))),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("graphics7", strong("Graphing options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphics7",
                                                  sliderInput("PlotWidth2", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                                  numericInput("PlotRes2", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                                 #sliderInput("PlotWidth2", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                 #numericInput("PlotRes2", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1)
                )
                
                
                
              ),
              
              mainPanel(
                
                tabsetPanel(id = "start3",
                            
                            tabPanel("Significant Comparisons Overview",
                                     downloadButton("downloadSigComps", "Download Table"),
                                     dataTableOutput('CompOverview')),
                            
                            tabPanel(style="height: 85vh; overflow: auto","Fold Change Plot and Venn Diagram",
                                     downloadButton("downloadPlot2", "Download Figure"),
                                     #plotOutput('foldchange1', dblclick = "plot_dbleclick", brush = brushOpts(
                                     #id = "plot_brush", 
                                     #resetOnNew = TRUE)),
                                     plotOutput("foldchange1"),
                                     uiOutput("venn.download"),
                                     #verbatimTextOutput("info"),
                                     uiOutput('venn'),
                                     downloadButton("downloadTable2", "Download Table"),
                                     dataTableOutput('MultipleCompTab')),
                            
                            tabPanel(style="height: 85vh; overflow: auto","Individual Gene Set Plot and Table",
                                     downloadButton("downloadPlot3", "Download Figure"),
                                     plotOutput('GeneSetPlot'),#, dblclick = "plot1_dbleclick", brush = brushOpts(
                                     #id = "plot1_brush",
                                     #resetOnNew = TRUE)),
                                     #verbatimTextOutput("info"),
                                     downloadButton("downloadTable3", "Download Table"),
                                     dataTableOutput('GeneSetTab'))        
                )
              )
            )
    ),
    
      tabItem(tabName = "roast",
              sidebarLayout(
                sidebarPanel(
                  conditionalPanel(condition = "input.startR == 'Significant Comparisons Overview'",
                                   uiOutput("setStat"),
                                   numericInput("SigLevelR", "Significance threshold:", min = 0, max = 1, step = .025, value = .05)
                                   #selectInput("foldchange.r", strong("Fold change sign:"), choices = c("All (include 0)", "+","-"), selected = "All (include 0)")
                  )
                ),
                mainPanel(
                  tabsetPanel(id = "startR",
                              tabPanel(style="height: 85vh; overflow: auto","Significant Comparisons Overview",
                                       downloadButton("downloadSigCompsR", "Download Table"),
                                       dataTableOutput('CompOverviewR')))
              )
      )),
		
      tabItem(tabName = "flow",
              sidebarLayout(
                sidebarPanel(
                  conditionalPanel(condition = "(input.start4 == 'Flow Analysis Lists')",
                                   uiOutput("CompF"),
                                   selectInput("methodF", "Multiple testing correction:",
                                               c("FDR", "Bonferroni", "Raw"), "Raw"),
                                   numericInput('alphaFlow_2',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025)
                  ),
                  conditionalPanel(condition="input.start4 == 'Significant Comparisons Overview'",
                                   uiOutput("FlowVariableSets"),
                                   numericInput('alphaFlow_1',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025)
                  ),
                  
                  conditionalPanel(condition = "input.start4 == 'Flow Analysis Lists'"#,
                                   
                  ),          
                  
                  conditionalPanel(condition="input.start4 == 'Flow Analysis Figures'",
                                   uiOutput("FlowDataNames"),
                                   uiOutput("flowsub"),               
                                   checkboxInput("FlowTransform","lg2 Transform",FALSE),
                                   helpText("Options for 1st Plot:"),
                                   uiOutput("flowmax"),
                                   uiOutput("flowmin"),
                                   checkboxInput("FlowSamples","Individual curves",FALSE),
                                   helpText("Options for 2nd Plot:"),
                                   checkboxInput("flowbox","Box Plot View",FALSE),
                                   downloadButton('downloadFlowData', 'Download Data')
                  )),
                 
                mainPanel(
                  
                  tabsetPanel(id = "start4",
                              tabPanel("Significant Comparisons Overview",
                                       helpText("Cell value represents the number of significant probes."),
                                       downloadButton('downloadFC', 'Download Table'),
                                       dataTableOutput("FlowOverview")
                              ),
                              
                              tabPanel(style="height: 85vh; overflow: auto","Flow Analysis Lists",
                                       helpText("Right click on hyperlinks to open in new window"),
                                       downloadButton('downloadFL', 'Download Table'),
                                       dataTableOutput("flowlisttable")
                              ),
                              tabPanel(style="height: 85vh; overflow: auto","Flow Analysis Figures",
                                       downloadButton('downloadFlowPlot', 'Download Figure'),
                                       plotOutput("FlowPlot"),
                                       downloadButton('downloadFlowPlot2', 'Download Figure'),
                                       plotOutput("FlowPlot2"),
                                       downloadButton('downloadFlowSummaries','Download Table'),
                                       dataTableOutput("FlowPlotSummary")
                              )
                  )))),
    
    tabItem(tabName = "pcaM",
            sidebarLayout(
              sidebarPanel(
                uiOutput("MetabVariableSets1"),#),   
                uiOutput("pcaAnnotM"),
                uiOutput("pcaAnnotValsM"),
                uiOutput("pcaBlockingM"),
                selectInput("PCSnumM", "Select number of PCs to plot:", choices = c(2,3), selected = 2),
                div(style = "display:inline-block", checkboxInput("graphicsPCAM", strong("Graphing options", style = "color:blue"), FALSE)),
                div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                conditionalPanel(condition = "input.graphicsPCAM",
                                 sliderInput("PlotWidthPCAM", "Plot Width", min = 500, max = 2000, value = 800, step = 50),
                                 sliderInput("PlotHeightPCAM", "Plot Height", min = 500, max = 2000, value = 500, step = 50),
                                 numericInput("CircleSizePCAM","Circle size",min=3,max=15,value=7,step=1),
                                 numericInput("PlotResPCAM", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
              ),
              mainPanel(style="height: 90vh; overflow: auto",
                        actionButton("showM", "View/hide screeplot", style = "color: white; background-color: green"),
                        uiOutput("screePlotM"),
                        br(),
                        plotOutput("PCAplotM"),
                        br(),
                        uiOutput("PCAplot3dM")
              )
            )
    ),
    
    
    tabItem(tabName = "metab",
            sidebarLayout(
              sidebarPanel(
                conditionalPanel(condition = "(input.startM == 'Metabolomics Analysis Lists')||(input.startM == 'DGE Heatmap')",
                                 uiOutput("CompM"),
                                 selectInput("methodM", "Multiple testing correction:",
                                             c("FDR", "Bonferroni", "Raw"), "Raw"),
                                 numericInput('alphaMetab_2',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025)
                ),
                conditionalPanel(condition="input.startM == 'Significant Comparisons Overview'",
                                 uiOutput("MetabVariableSets"),
                                 numericInput('alphaMetab_1',"Significance threshold:", min = 0, max = 1, value = 0.05, step = 0.025)
                ),
                
                conditionalPanel(condition = "input.startM == 'Metabolomics Analysis Lists'"#,
                                 
                ),          
                
                conditionalPanel(condition="input.startM == 'Metabolomics Analysis Figures'",
                                 uiOutput("MetabDataNames"),
                                 uiOutput("metabsub"),               
                                 checkboxInput("MetabTransform","lg2 Transform",FALSE),
                                 helpText("Options for 1st Plot:"),
                                 uiOutput("metabmax"),
                                 uiOutput("metabmin"),
                                 checkboxInput("MetabSamples","Individual curves",FALSE),
                                 helpText("Options for 2nd Plot:"),
                                 checkboxInput("metabbox","Box Plot View",FALSE),
                                 downloadButton('downloadMetabData', 'Download Data')
                ),
                conditionalPanel(condition="input.startM == 'DGE Heatmap'",
                                 uiOutput("normMetab"),
                                 numericInput("setcutoff2","Select the maximum values on color key:",2),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("rowClustM", strong("Row Cluster"), value = FALSE)),
                                 div(style = "display:inline-block", infoPopup("Row Cluster", "Based on the PC performance, it might not be recommended to cluster the rows if the number of rows exceed 7000",
                                                                               placement = "right", trigger = "click")),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("subsetMetab",strong("Subset column options", style = "color:blue"),FALSE)),
                                 div(style = "display:inline-block", helpPopup("Subset column options", "These options allow the user to plot a subset of the samples within the expression data set.")),
                                 conditionalPanel(condition="input.subsetMetab",
                                                  uiOutput('subsetMetabVariable'),
                                                  uiOutput('subsetMetabValue'),
                                                  uiOutput('subsetMetabVariable2'),
                                                  uiOutput('subsetMetabValue2')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('ordercolumnsM', strong("Ordering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Ordering column options", "These options allow the user to specify column ordering by column variable.")),
                                 conditionalPanel(condition = "input.ordercolumnsM",
                                                  uiOutput('TopTierMetab'),
                                                  uiOutput('MidTierMetab'),
                                                  uiOutput('LowTierMetab'),
                                                  uiOutput('group_label_metab')),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput('clusteroptionsM', strong("Clustering column options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Cluster column options", "These options allow the user to cluster the samples and label the groups they are clustered in.")),
                                 conditionalPanel(condition = "input.clusteroptionsM",
                                                  checkboxInput("ColClustMetab","Cluster samples (columns)",FALSE),
                                                  #uiOutput('ColumnClusterProbe'),
                                                  checkboxInput('ClusterChoiceM', "Show cluster groups", FALSE),
                                                  conditionalPanel(condition = "input.ClusterChoiceM",
                                                                   div(style = "display:inline-block", uiOutput('ClusterCutsM')),
                                                                   div(style = "display:inline-block", infoPopup("Cluster Cuts", "Optimal number calculated using Dunn's Index.",
                                                                                                                 placement = "right", trigger = "click")))),
                                 br(),
                                 div(style = "display:inline-block", checkboxInput("graphicsM", strong("Graphing options", style = "color:blue"), FALSE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphicsM",
                                                  sliderInput('Metab_HeatMapSizew', "Grid width", min = 500, max = 10000, value = 750, step = 50),
                                                  sliderInput('Metab_HeatMapSizeh', "Grid height", min = 500, max = 10000, value = 650, step = 50),
                                                  numericInput('Metab_FontSize', "Font size:", min = 10, value = 10, step = 2),
                                                  numericInput('Metab_LegendSize', "Legend size:", min = 1, max = 5, value = 1, step = .50),
                                                  numericInput('Metab_TreeHeight', "Tree height:", min = 0, value = 50, step = 5),
                                                  numericInput('Metab_PlotRes', "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                )
                
              ),
              mainPanel(
                
                tabsetPanel(id = "startM",
                            tabPanel("Significant Comparisons Overview",
                                     helpText("Cell value represents the number of significant probes."),
                                     downloadButton('downloadFC2', 'Download Table'),
                                     dataTableOutput("MetabOverview")
                            ),
                            
                            tabPanel(style="height: 85vh; overflow: auto","Metabolomics Analysis Lists",
                                     helpText("Right click on hyperlinks to open in new window"),
                                     downloadButton('downloadME', 'Download Table'),
                                     dataTableOutput("metablisttable")
                            ),
                            tabPanel(style="height: 85vh; overflow: auto","Metabolomics Analysis Figures",
                                     downloadButton('downloadMetabPlot', 'Download Figure'),
                                     plotOutput("MetabPlot"),
                                     downloadButton('downloadMetabPlot2', 'Download Figure'),
                                     plotOutput("MetabPlot2"),
                                     downloadButton('downloadMetabSummaries','Download Table'),
                                     dataTableOutput("MetabPlotSummary")
                            ),
                            tabPanel(style="height: 100vh; overflow: auto","DGE Heatmap",
                                     downloadButton("downloadHeatmapm"),
                                     plotOutput("heatmapm"))
                )))),
    
    tabItem(tabName = "corr",
            sidebarLayout(
              sidebarPanel(
                conditionalPanel(condition = "input.start5 == 'Correlation Overview Heatmap'",
                                 uiOutput("TypeVariable"),
                                 checkboxInput("var_switch", strong("Swap 'Variable' and 'With' columns"), value = FALSE),
                                 uiOutput("uploadmod1"),
                                 uiOutput("fileupload1"),
                                 uiOutput("WithVariable1"),
                                 checkboxInput("Base_subtract", strong("Baseline subtracted"), value = FALSE),
                                 div(style = "display:inline-block", checkboxInput("TSoptions2", strong("Testing and correlation value subsetting options", style = "color:blue"), FALSE)),
                                 conditionalPanel(condition = "input.TSoptions2",
                                                  selectInput("correction_method.corr","Multiple testing correction:",c("FDR","Bonferroni","Raw"),"Raw"),
                                                  numericInput("Alpha1", "Significance threshold:", min = 0, max = 1, step = .01, value = 1),
                                                  checkboxInput("corrval2",label = strong("Filter on correlation value:"),value = FALSE),
                                                  conditionalPanel(condition = "input.corrval2 == true",
                                                                   numericInput("corrval3", label = "Correlation value cut off:",
                                                                                min = 0, max = 1, value = 0, step = 0.1),
                                                                   selectInput("corrsign1","Correlation value sign:", c("+", "-", "Both"), c("Both")))
                                 ),
                                 checkboxInput("subsetModcorr",strong("Subset column options:", style = "color:blue"),FALSE),
                                 conditionalPanel(condition="input.subsetModcorr",
                                                  uiOutput("subsetcorr")),
                                 checkboxInput("ordercorr", strong("Order columns by time"), FALSE),
                                 div(style = "display:inline-block", checkboxInput("rowclustcorr", strong("Row cluster"),FALSE)),
                                 div(style = "display:inline-block",infoPopup("Row Cluster", "Based on the PC performance, it might not be recommended to cluster the rows if the number of rows exceed 7000",
                                                                               placement = "right", trigger = "click")),
                                 checkboxInput("colclustcorr", strong("Column cluster"), FALSE),
                                 div(style = "display:inline-block", checkboxInput("graphics9", strong("Graphing options", style = "color:blue"), TRUE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options:", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphics9",
                                                  sliderInput("PlotWidth4", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                                  sliderInput("PlotHeight4", "Plot Height", min = 500, max = 5000, value = 600, step = 100),
                                                  numericInput("PlotRes4", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                ),
                
                conditionalPanel(condition = "input.start5 == 'Correlation Analysis Table'",
                                 uiOutput("TypeVariable2"),
                                 checkboxInput("var_switch2", strong("Swap 'Variable' and 'With' columns"), value = FALSE),
                                 #selectInput("TypeVariable2", "Choose type:", choices = c("Micro vs Serology", "Flow vs Serology", "Micro vs Flow", "Module vs Serology", "Module vs Flow"), selected = "Micro vs Serology"),
                                 uiOutput("WithVariable2"),
                                 checkboxInput("Base_subtract2", "Baseline subtracted", value = FALSE),
                                 #selectInput("WithVariable2", "Choose 'with' variable:", choices = as.character(unique(correlations$With)), selected = unique(correlations$With)[1]),
                                 div(style = "display:inline-block", checkboxInput("TSoptions3", strong("Testing and correlation value subsetting options", style = "color:blue"), FALSE)),
                                 conditionalPanel(condition = "input.TSoptions3",
                                                  selectInput("correction_method.corr2","Multiple testing correction:",c("FDR","Bonferroni","Raw"),"Raw"),
                                                  numericInput("Alpha", "Significance threshold:", min = 0, max = 1, step = .01, value = .05),
                                                  checkboxInput("corrval",label = strong("Filter on correlation value:"),value = FALSE),
                                                  conditionalPanel(condition = "input.corrval == true",
                                                                   numericInput("corrval1", label = "Correlation value cut off:",
                                                                                min = 0, max = 1, value = 0, step = 0.1),
                                                                   selectInput("corrsign","Correlation value sign:", c("+", "-", "Both"), c("Both")))
                                 ),
                                 uiOutput("Visit")
                ),
                
                conditionalPanel(condition = "input.start5 == 'Correlation Analysis Plot'",
                                 uiOutput("TypeVariable4"),
                                 checkboxInput("var_switch3", strong("Swap 'Variable' and 'With' columns"), value = FALSE),
                                 uiOutput("Visit2"),
                                 uiOutput("uploadmod3"),
                                 uiOutput("fileupload2"),
                                 checkboxInput("Base_subtract3", "Baseline subtracted", value = FALSE),
                                 div(style = "display:inline-block", checkboxInput("TSoptions4", strong("Testing and correlation value subsetting options", style = "color:blue"), FALSE)),
                                 conditionalPanel(condition = "input.TSoptions4",
                                                  selectInput("correction_method.corr3","Multiple testing correction:",c("FDR","Bonferroni","Raw"),"Raw"),
                                                  numericInput("Alpha2", "Significance threshold:", min = 0, max = 1, step = .01, value = .05),
                                                  checkboxInput("corrval4",label = strong("Filter on correlation value:"),value = FALSE),
                                                  conditionalPanel(condition = "input.corrval4 == true",
                                                                   numericInput("corrval5", label = "Correlation value cut off:",
                                                                                min = 0, max = 1, value = 0, step = 0.1),
                                                                   selectInput("corrsign2","Correlation value sign:", c("+", "-", "Both"), c("Both")))
                                 ),
                                 checkboxInput("subsetModcorr2",strong("Subset column options:", style = "color:blue"),FALSE),
                                 conditionalPanel(condition="input.subsetModcorr2",
                                                  uiOutput("subsetcorr2")),
                                 div(style = "display:inline-block", checkboxInput("graphics8", strong("Graphing options", style = "color:blue"), TRUE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options:", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphics8",
                                                  sliderInput("PlotWidth3", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                                  sliderInput("PlotHeight3", "Plot Height", min = 500, max = 5000, value = 600, step = 100),
                                                  numericInput("PlotRes3", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                ),
                
                conditionalPanel(condition = "input.start5 == 'Pairwise Scatterplots'",
                                 uiOutput("TypeVariable6"),
                                 uiOutput("Visit3"),
                                 uiOutput("corr_Variable"),
                                 uiOutput("WithVariable4"),
                                 checkboxInput("log_scale", "lg2 Transform", value = FALSE),
                                 checkboxInput("plot_reg", "Plot Regression Line", value = FALSE),
                                 checkboxInput("plot_loess", "Fit loess curve", value = FALSE),
                                 conditionalPanel(condition = "input.plot_loess",
                                                  numericInput("span", "Span (0 to 1):", min = 0, max = 1, value = 1, step = .05)),
                                 div(style = "display:inline-block", checkboxInput("graphics10", strong("Graphing options", style = "color:blue"), TRUE)),
                                 div(style = "display:inline-block", helpPopup("Graphing options:", "These options allow the user to make adjustments to the plot (i.e. width and height).")),
                                 conditionalPanel(condition = "input.graphics10",
                                                  numericInput("Point_size", "Circle size (1 to 5):",min=1,max=5,value=2,step=.5),
                                                  numericInput("axis_text_size", "Axis text size:", min = 10, value = 10, step = 1),
                                                  numericInput("axis_label_size", "Axis label size:", min = 12, value = 12, step = 1),
                                                  sliderInput("PlotWidth5", "Plot Width", min = 500, max = 5000, value = 800, step = 100),
                                                  sliderInput("PlotHeight5", "Plot Height", min = 500, max = 5000, value = 600, step = 100),
                                                  numericInput("PlotRes5", "Plot resolution (for downloaded plot):", min = 72, max = 400, value = 72, step = 1))
                )
              ),
              
              mainPanel(
                tabsetPanel(id = "start5",
                            tabPanel(style="height: 100vh; overflow: auto","Correlation Overview Heatmap",
                                     downloadButton('download_heatmap_data', "Download Data"),
                                     downloadButton('download_corr_heatmap', "Download Figure"),
                                     plotOutput("correlations_plotOverview"),
                                     dataTableOutput("overviewTable")
                            ),
                            
                            tabPanel(style="height: 85vh; overflow: auto","Correlation Analysis Table",
                                     downloadButton('download_data', "Download Table"),
                                     dataTableOutput("correlation_table")#,
                                     #plotOutput("correlations_plot")
                            ),
                            
                            tabPanel(style="height: 100vh; overflow: auto","Correlation Analysis Plot",
                                     downloadButton('download_plot_data', "Download Data"),
                                     downloadButton('download_corr_plot', "Download Figure"),
                                     plotOutput("correlations_plot"),
                                     tableOutput("correlations_table")
                            ),
                            
                            tabPanel(style="height: 100vh; overflow: auto","Pairwise Scatterplots",
                                     downloadButton('download_scatter_data', "Download Data"),
                                     downloadButton('download_scatter', "Download Figure"),
                                     plotOutput("correlations_scatter_plot")
                            )
                            
              ))))))


shinyUI(dashboardPage(skin="blue",
  header,
  sidebar,
  body
))

}


