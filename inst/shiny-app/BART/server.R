#require(shiny)
#library(ggplot2)
#library(RColorBrewer)
#library(fastcluster)
#library(NMF)
#library(grid)
#library(clValid)
#library(qusage)
#library(VennDiagram)
#library(shinydashboard)
#library(gtools)
#library(scales)
#library(reshape2)
#library(data.table)
#library(pca3d)
#library(shinyjs)
#library(stringr)
library(ggplot2)

#Increase maximum file size input; currently set to 1GB -- Not sure if we need more.
options(shiny.maxRequestSize=1000*1024^2)
shinyServer(function(input,output, session){

  setBookmarkExclude(c("file1"))
  output$test<-renderUI({

    if(is.null(values$h3_rowdendro)){
      if(hc ==TRUE){
        try<-list("Baseline Median Normalized" = 1,
                  "Baseline Healthy Normalized" = 2)
      }

      if(hc == FALSE){
        try<-list("Baseline Median Normalized" = 1)
      }
    }

    if(!is.null(values$h3_rowdendro)){
      if(values$hc == TRUE){
        try<-list("Baseline Median Normalized" = 1,
                  "Baseline Healthy Normalized" = 2,
                  "All Samples Median Normalized" = 3,
                  "All Samples Healthy Normalized"= 4,
                  "All Samples Baseline Normalized"=5)
      }

      if(values$hc == FALSE){
        try<-list("Baseline Median Normalized" = 1,
                  "All Samples Median Normalized" = 3,
                  "All Samples Baseline Normalized"=5)
      }
    }

    selectInput("set", "Select heatmap:",as.list(try))
  })

  dir.projects <- reactive({
    "C:/Users/e89628/Documents/Intern/BartDownloadFilesForTesting/BART"
  })

  list.projects <- reactive({
    list.files(dir.projects())
  })

  get.wd <- eventReactive(input$uploadserver,{return(getwd())})


  fileindex<-reactive({
    index1<-c()
    index1[1]<-which(input$file1$name == "pvca1.png")
    index1[2]<-which(input$file1$name == "pvca2.png")
    index1[3]<-which(input$file1$name == "Unsupervised.Rdata")
    index1
  })

  values <- reactiveValues()

  updateData <- function(path, route){
    vars <- load(file = path, envir = .GlobalEnv)
    values$final_expression <- NULL
    values$qusage_results <- NULL
    values$flow_data <- NULL
    values$results_file <- NULL
    values$results_file2 <- NULL
    values$mod1 <- NULL
    values$mod2 <- NULL
    values$illumina <- NULL
    values$moduleinfo2 <- NULL
    values$correlations <- NULL
    values$ModulesTF <- NULL
    values$BaylorTF <- NULL
    values$base_mod <- NULL
    values$long_mod <- NULL
    values$long_mod2 <- NULL
    values$roast_results <- NULL
    values$numbGenes <- NULL
    values$PaloValue <- NULL
    values$quantileValue <- NULL
    values$SDcutoff <- NULL
    values$Approx_num_probes <- NULL
    values$Filter_text <- NULL
    values$BatchCorrectionTF <- NULL
    values$Batch_Variable <- NULL
    values$Batch_text <- NULL
    values$MixedModelDescription <- NULL
    values$results_file2 <- NULL
    values$FlockFlowMMR <- NULL
    values$SpadeFlowMMR <- NULL
    values$FlockFlowData <- NULL
    values$SpadeFlowData <- NULL
    setwd(route)

    for (var in vars){
      values[[var]] <- get(var, .GlobalEnv)
    }
    return(values)
  }

  observe({
    if(is.null(input$file1)) return(NULL)
    mypath<-input$file1[[fileindex()[3], 'datapath']]
    updateData(mypath, tempdir())
  })

  output$versionbox <- renderValueBox({
    valueBox(
      paste0("Software Version: BETA \nRelease Date: TBD"), version$version.string, icon = icon("th-list"),
      color = "purple"
    )
  })

  output$version <- renderPrint({
    writeLines("Software Version: BETA \nRelease Date: TBD")
    writeLines(version$version.string)
  })

  output$path <- renderText({
    path <- values$project_name
    path
  })


  ############ Design File and Summary Statistics ##############################

  output$SumStat <- renderMenu({
    if(is.null(values$design)){
      return(strong(""))
    }
    if(is.null(values$design) == FALSE){
      if(length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE))) == 0){
        return(menuItem("Summary Stats and QC", icon = icon("th-list"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) == 0){
        return(menuItem("Summary Stats and QC", icon = icon("th-list"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("PVCA", tabName = "pvca")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) == 0){
        return(menuItem("Summary Stats and QC", icon = icon("th-list"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("PVCA", tabName = "pvca"),
                        menuSubItem("Fast QC Report", tabName = "fastqc")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) > 0){
        return(menuItem("Summary Stats and QC", icon = icon("th-list"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("PVCA", tabName = "pvca"),
                        menuSubItem("QC files upload", tabName = "qcupload")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) == 0){
        return(menuItem("Summary Stats and QC", icon = icon("th-list"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("Fast QC Report", tabName = "fastqc")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) > 0){
        return(menuItem("Summary Stats and QC", icon = icon("th-list"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("QC files upload", tabName = "qcupload")))
      }
      if(length(which(input$file1$name %in% grep(".png", input$file1$name, value = TRUE))) == 0 &
         length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))) > 0 &
         length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) > 0){
        return(menuItem("Summary Stats and QC", icon = icon("th-list"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("Fast QC Report", tabName = "fastqc"),
                        menuSubItem("QC files upload", tabName = "qcupload")))
      }
      else{
        return(menuItem("Summary Stats and QC", icon = icon("th-list"), tabName = "summarystats",
                        menuSubItem("Design file", tabName = "design"),
                        menuSubItem("Summary Statistics", tabName = "summary"),
                        menuSubItem("PVCA", tabName = "pvca"),
                        menuSubItem("Fast QC Report", tabName = "fastqc"),
                        menuSubItem("QC files upload", tabName = "qcupload")))
      }
    }
  })


  output$designDataTable<- renderDataTable({
    values$design
  })

  output$downloadDesign <- downloadHandler(

    filename = function() {paste(values$project_name,'_Design','.csv', sep='')  },
    content = function(file) {
      write.csv(values$design, file, row.names = FALSE)
    }
  )

  output$downloadSummary0 <- downloadHandler(

    filename = function() {'SummaryStats_Table1.csv'},
    content = function(file) {
      write.csv(table_0(), file, row.names = FALSE)
    }
  )

  output$downloadSummary1 <- downloadHandler(

    filename = function() {'SummaryStats_Table2.csv'},
    content = function(file){
      write.csv(table_1(), file,row.names = FALSE)
    }
  )

  output$downloadSummary2 <- downloadHandler(

    filename = function() {'SummaryStats_Table3.csv'},
    content = function(file) {
      write.csv(table_2(), file,row.names = FALSE)
    }
  )

  output$summaryName<-renderUI({
    selectInput("summary_var", "Variable for summary statistics:", names(values$design),selected=values$summary_var)
  })

  output$timeVar<-renderUI({
    selectInput("time_var", "Time variable:*", values$time_var,selected=values$time_var)
  })

  output$respVar<-renderUI({
    selectInput("responder_var", "By variable (i.e. responder status):", names(values$design),selected=values$responder_var)
  })

  output$patientVar<-renderUI({
    selectInput("patient_id", "Patient ID variable:*", values$patient_id,selected=values$patient_id)
  })

  output$summary0Text<-renderUI({
    HTML(paste(strong("Table 1:"), "Summary statistics by", input$responder_var, "for",input$summary_var,sep=" "))
  })

  output$summary1Text<-renderUI({
    HTML(paste(strong("Table 2: "), "Summary statistics by Time and", input$responder_var, "for",input$summary_var,sep=" "))
  })

  output$summary2Text <- renderUI({
    HTML(paste(strong("Table 3: "), "Frequency of Longitudinal Profiles by", input$responder_var, sep = " "))
  })

  data_unique <- reactive({
    values$design[!duplicated(values$design[, values$patient_id]), ]
  })

  table0_cont <- reactive({
    sum0 <- aggregate(as.formula(paste(input$summary_var, "~", values$patient_id)), data = values$design, sd)
    sum0 <- as.data.frame(sum0)
    if(sum(sum0[[2]], na.rm = T) == 0){
      mysummary<-function(x){y<-c(length(x),mean(x),median(x),sd(x),min(x),max(x))
      names(y)<-c("N","Mean","Median","Sd","Min.","Max.")
      return(y)}

      sum1<-aggregate(as.formula(paste(input$summary_var,"~",input$responder_var)), data= data_unique(), mysummary)
      summaries<- as.data.frame(sum1[,2])
      result<-as.data.frame(cbind(sum1[,c(1)],summaries))
      names(result)<-c("Responder Status",colnames(summaries))
      result
    }
    else{
      result <- data.frame(Warning = "Variable is time dependent. Look at table below.")
    }
  })

  table0_cat <- reactive({
    sum0 <- aggregate(as.formula(paste(input$summary_var, "~", values$patient_id)), data = values$design, summary)
    sum01 <- sum0[,2]
    sum01[sum01 > 0] <- 1
    sums <- apply(sum01, 1, sum)
    check <- F
    if(all(sums == 1)) {check <- T}

    if(check == T){

      mysummary<-function(x){y<-summary(x)
      y<-c(length(x),y)
      names(y)<-c("N", names(summary(x)))
      return(y)}

      sum1<-aggregate(as.formula(paste(input$summary_var,"~",input$responder_var)), data = data_unique(), mysummary)
      summaries<- as.data.frame(sum1[,2])
      summaries2 = list()
      for(i in 1:nrow(summaries[,-1])){
        summaries2[[i]] = paste(summaries[i,-1], " ", "(", round(100*(summaries[i,-1]/summaries[i,1]), 2), "%",")", sep = "")
      }

      summaries2 = do.call("rbind", summaries2)

      result<- cbind(as.data.frame(sum1[,1]),summaries[,1], summaries2)
      names(result)<-c("Responder Status", "N", names(summaries[,-1]))
      result
    }

    else{
      result <- data.frame(Warning = "Variable is time dependent. Look at table below.")
    }
  })

  table_0 <- reactive({
    if(is.numeric(values$design[[input$summary_var]])){
      table0_cont()
    }

    else{
      table0_cat()
    }
  })

  table_1 <- reactive({
    if(is.numeric(values$design[[input$summary_var]])){
      mysummary<-function(x){y <- c(length(x), mean(x), median(x), sd(x), min(x), max(x))
      names(y)<-c("N", "Mean", "Median", "Sd", "Min.","Max.")
      return(y)}
      sum1<-aggregate(as.formula(paste(input$summary_var,"~",values$time_var,"+",input$responder_var)), data = values$design, mysummary)
      summaries<- as.data.frame(sum1[,3])

      result<- data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("Time","Responder Status",colnames(summaries))
      result
    }

    else{
      mysummary<-function(x){y<-summary(x)
      y<-c(length(x),y)
      names(y)<-c("N",names(summary(x)))
      return(y)}
      sum1<-aggregate(as.formula(paste(input$summary_var,"~",values$time_var,"+",input$responder_var)), data= values$design, mysummary)
      summaries<- as.data.frame(sum1[,3])
      result<- data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("Time","Responder Status",colnames(summaries))
      result
    }
  })

  table_2 <- reactive({
    x<-values$design
    names(x)[which(names(x)==values$time_var)]<-"Time"
    names(x)[which(names(x)==input$responder_var)]<-"Responder"
    names(x)[which(names(x)==values$patient_id)]<-"PATIENT_ID"

    timelevel<-unique(x$Time[order(x$Time)])
    x$Time2<-factor(x$Time,levels=timelevel)
    freqs<-aggregate(Time2~PATIENT_ID,data=x,table)

    y<-unlist(apply(as.matrix(freqs[,-1]),1,function(x){paste(timelevel[which(x>0)],collapse=",")}))

    newx<-data.frame(freqs[,1],y)
    names(newx)<-c("PATIENT_ID","course")

    dum<-c()
    for(i in 1:dim(newx)[1]){
      dum[i]<-which(x$PATIENT_ID==newx$PATIENT_ID[i])[1]}
    n1<-names(newx)
    n2<-"Responder"
    newx<-data.frame(newx,x$Responder[dum])
    names(newx)<-c(n1,n2)
    z<-aggregate(Responder~course,data=newx,table)

    if(length(z[,2]) == 1){
      z2<- as.data.frame(z[,2,drop = F])
    }

    if(length(z[,2]) > 1){
      z2 <- as.data.frame(z[,2])
    }

    z[,1] <- as.character(z[,1])
    sum_cols <- colSums(z2)
    result <- rbind(data.frame(z[,1], stringsAsFactors = F), "Total")

    if(length(z[,2] == 1)){
      result2 <- rbind(data.frame(z[,2,drop = F]), sum_cols)
    }

    if(length(z[,2]) > 1){
      result2 <- rbind(data.frame(z[,2]), sum_cols)
    }

    result3<- data.frame(cbind(result,result2))
    names(result3)<-c("Time Course",colnames(z2))
    result3$Total <- rowSums(result3[,-1])
    result3
  })

  dig <- function(){
    as.numeric(input$digits)
  }

  table_select <- reactive({
    output$summary0 <- renderTable({
      table_0()
    }, digits = dig(), include.rownames = FALSE)

    output$summary1 <- renderTable({
      table_1()
    }, digits = dig(), include.rownames = FALSE)

    output$summary2 <- renderTable({
      table_2()
    }, digits = dig(), include.rownames = FALSE)

  })

  output$summary_tab <- renderTable({
    table_select()
  })


  ################# PVCA ##################################

  output$PVCAtext <- renderText({
    if(BatchCorrectionTF == TRUE){
      return(paste0("The data set was processed as follows. Expression data is background subtracted and scaled,
                    raw values less than 10 are set to 10, and then log2 transformed. PALO filtering (Present in At
                    Least One sample) using a p-value detection threshold of 0.01 was conducted. Standard deviation
                    filtering info: ",  SDcutoff, ". This resulted in a final expression set containing ", Approx_num_probes,
                    " The expression file was then assessed for batch effects that contribute to technical variablility
                    in the data such as chip to chip variability or hybridization chamber. The results of that assessment
                    are presented in a PVCA analysis below. Batch correction is implicated if the figures provide a before
                    and after batch correction assessment."))
    }
    return(paste0("The data set was processed as follows. Expression data is background subtracted and scaled, raw
                  values less than 10 are set to 10, and then log2 transformed. PALO filtering (Present in At Least One sample)
                  using a p-value detection threshold of 0.01 was conducted. Additional filtering on standard deviation was not
                  conducted.  This resulted in a final expression set containing ", Approx_num_probes, " No batch effect correction
                  was applied to this data set. The results of the batch effect assessment are presented in a PVCA analysis below
                  if available."))
  })


  output$PVCA1 <- renderImage({
    if (is.null(input$file1)){return(NULL)}
    if(is.na(fileindex()[1]) == TRUE){
      filename <- 'datapath'
    }
    else{
      filename<-input$file1[[fileindex()[1], 'datapath']]
    }

    # Return a list containing the filename and alt text
    list(src = filename,
         alt = "Principle Component 3D PLOT",
         height=300,
         width=700)

  }, deleteFile = FALSE)



  output$PVCA2 <- renderImage({
    if (is.null(input$file1)){return(NULL)}
    if(is.na(fileindex()[2]) == TRUE){
      filename <- 'datapath'
    }
    else{
      filename<-input$file1[[fileindex()[2], 'datapath']]
    }

    # Return a list containing the filename and alt text
    list(src = filename,
         alt = "PVCA Summary")

  }, deleteFile = FALSE)


  ######################### QC Metrics ################################
  output$qc_select <- renderUI(
    if(length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))) < 1){
      return(NULL)
    }
    else{
      qcnames <- grep(".html", input$file1$name, value = TRUE)
      if(length(qcnames) > 1){
        names <- c()
        for(i in 1:length(qcnames)){
          names[i] <- sub("_S.*", "", qcnames[i])
        }
        selectInput("fqc_samples", "Select sample:", names)
      }
    }
  )

  output$qc_select2 <- renderUI(
    if(length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))) == 0){return(NULL)}
    else{
      qcnames <- grep(".html", input$file1$name, value = TRUE)
      if(length(qcnames) > 1){
        selectInput("fqc_read", "Select read:", c("R1", "R2"))
      }
    }
  )

  output$fqc <- renderUI(
    if(length(which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))) == 0){return(NULL)}
    else{
      qcnames <- grep(".html", input$file1$name, value = TRUE)
      index <- which(input$file1$name %in% grep(".html", input$file1$name, value = TRUE))
      if(length(qcnames) == 1){
        return(includeHTML(input$file1$datapath[index]))
      }
      if(length(qcnames) > 1){
        names <- c()
        check <- c()
        for(i in 1:length(qcnames)){
          names[i] <- sub("_S.*", "", qcnames[i])
        }
        for(i in 1:length(qcnames)){
          check[i] <- grepl(input$fqc_read, qcnames[i]) & grepl(input$fqc_samples, qcnames[i])
        }
        return(includeHTML(input$file1$datapath[index][check]))
      }
    }
  )

  output$css_colors <- renderUI(
    if(is.null(input$file1)){return(NULL)}
    else{
      tags$head(
        tags$style(type="text/css", "table th{background-color:#ecf0f5;color:black;}"),
        tags$style(type = "text/css", "table td{background-color:#ecf0f5;")
      )
    }
  )


  output$dropdown <- renderUI(
    if(length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) == 0){return(NULL)}
    else{
      qcnames <- grep(".csv", input$file1$name, value = TRUE)
      index <- which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))
      selectInput("qcTab", "Select QC table:", qcnames)
    }
  )

  output$dropdowncolumns <- renderUI(
    if(length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) == 0){return(NULL)}
    else{
      qctab <- read.csv(input$file1$datapath[which(input$file1$name == input$qcTab)], header = TRUE)
      numeric <- sapply(qctab, is.numeric)
      qctab <- qctab[,numeric]
      selectInput("qcCols", "Select metrics to plot:", colnames(qctab))
    }
  )

  output$QCbox <- renderUI(
    if(length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) == 0){return(NULL)}
    else{
      return(plotOutput("qcbox"))
    }
  )

  output$coordflip <- renderUI(
    if(length(which(input$file1$name %in% grep(".csv", input$file1$name, value = TRUE))) == 0){return(NULL)}
    else{
      checkboxInput("coord_flip", "Flip coordinates", FALSE)
    }
  )

  qc.width <- function(){
    input$qcplotsizew
  }

  qc.height <- function(){
    input$qcplotsizeh
  }

  axis_label_size_qc <- reactive({
    input$axis_label_sizeqc
  })

  axis_text_size_qc <- reactive({
    input$axis_text_sizeqc
  })



  output$qcbox <- renderPlot({
    qctab <- read.csv(input$file1$datapath[which(input$file1$name == input$qcTab)], header = TRUE)
    qcplot <- ggplot(data = qctab[,c(1,which(colnames(qctab) %in% input$qcCols))], aes(x = X)) + geom_bar(aes_string(weight = input$qcCols), fill="blue4")
    qcplot <- qcplot + theme(axis.text.x = element_text(hjust = 1, angle = 90), axis.text = element_text(size = axis_text_size_qc()),
                             axis.title = element_text(size = axis_label_size_qc())) + xlab("Samples") + ylab(input$qcCols)
    if(input$coord_flip){
      qcplot <- qcplot + coord_flip()
      qcplot
    }
    else{
      qcplot
    }
  }, height = qc.height, width = qc.width)

  output$downloadQC <- downloadHandler(

    filename = function() {paste("QC_Plot_",input$qcCols, ".png", sep = "")},
    content = function(file) {
      png(file, width = input$qcplotsizew, height = input$qcplotsizeh)
      qctab <- read.csv(input$file1$datapath[which(input$file1$name == input$qcTab)], header = TRUE)
      qcplot <- ggplot(data = qctab[,c(1,which(colnames(qctab) %in% input$qcCols))], aes(x = X)) + geom_bar(aes_string(weight = input$qcCols), fill="blue4")
      qcplot <- qcplot + theme(axis.text.x = element_text(hjust = 1, angle = 90), axis.text = element_text(size = axis_text_size_qc()),
                               axis.title = element_text(size = axis_label_size_qc())) + xlab("Samples") + ylab(input$qcCols)
      if(input$coord_flip){
        qcplot <- qcplot + coord_flip()
        print(qcplot)
      }
      else{
        print(qcplot)
      }
      dev.off()
    }
  )

  ######################### Create Module Percentage Matrices ############################

  module_matrices <- reactive({
    if(is.null(values$final_expression)){
      mod1 <- mod2 <- mod3 <- NULL
    }
    else{
    if(values$data_type == "rna"){
      moduleinfo <- moduleinfo_rna
    }
    data_type <- values$data_type
    design <- values$design
    final_expression <- values$final_expression
    baseline_var <- values$baseline_var
    baseline_val <- values$baseline_val
    control_val <- values$control_val
    control_var <- values$control_var
    sample_id <- values$sample_id
    patient_id <- values$patient_id
    time_var <- values$time_var
    hc <- values$hc
    if(is.null(values$illumina)){illumina <- TRUE}
    else{illumina <- values$illumina}
    LongitudinalTF = values$long
    mod3 <- NULL

    if(is.null(control_var)){ mod1<-NULL}

    if(!is.null(control_var)){
      des<-design

      if(!illumina){
        if(data_type == "rna"){
          modexp<-merge(final_expression,moduleinfo,by=c("PROBE_ID"))
          index<-which(modexp$Module=="")
          if(length(index) > 0){
            modexp<-modexp[-index,]
          }
        }
        else{
          modexp<-merge(final_expression,moduleinfo2,by=c("SYMBOL"))
          index<-which(modexp$Module=="")
          if(length(index)>0){
            modexp<-modexp[-index,]}
        }
      }

      if(illumina){
        modexp<-merge(final_expression,moduleinfo,by=c("PROBE_ID"))
        index<-which(modexp$Module=="")
        if(length(index) > 0){
          modexp<-modexp[-index,]
        }
      }

      #Creating an ordered Module List and the total number of modules
      modnames<-unique(moduleinfo$Module)[-1]
      modnum<-gsub("M","",modnames)
      modvec<-as.numeric(unlist(strsplit(modnum,".",fixed=TRUE)))
      modmat<-matrix(modvec,ncol=2,byrow=T)
      modordered<-factor(modnames[order(modmat[,1],modmat[,2])],levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE)
      modordnum<-table(factor(moduleinfo$Module,levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE))
      #modordered is a list of the modules in order
      #modordnum is a list of the number of probes in each module "latest chip version"

      #inputs for function include des,reduced exp of just module annotated probes
      #and the following


      subjectid=patient_id

      #conforming to Xuans code
      base_sample_name=des$columnname[des[,baseline_var]==baseline_val] ### get columnname of baseline samples
      base_sample_grp=des[,control_var][des[,baseline_var]==baseline_val] ### the group info of all baseline samples in design file
      exp_base_sam=modexp[,colnames(modexp)%in%base_sample_name] ### get expression data of baseline samples

      #### connect expression file with design file by key information--group, donor, timepoint
      grp=base.donorid=base.sample=c()
      for (i in 1:ncol(exp_base_sam))
      {
        grp[i]=as.character(base_sample_grp[base_sample_name==colnames(exp_base_sam)[i]])
        ### baseline sample group info in expression file
        base.donorid[i]=as.character((des[,subjectid][des[,baseline_var]==baseline_val])[base_sample_name==colnames(exp_base_sam)[i]])
        ### baseline sample donorid info in expression file
        base.sample[i]=as.character((des[,sample_id][des[,baseline_var]==baseline_val])[base_sample_name==colnames(exp_base_sam)[i]])
        ### baseline sample timepoint info in expression file
      }
      unique.module<-modordered
      n.module<-length(unique.module)
      n.case=sum(grp!=control_val)
      exp_base_case=exp_base_sam[,grp!=control_val]
      n.base.case=ncol(exp_base_case)

      bhc_mean=apply(exp_base_sam[,grp==control_val],1,mean) ### baylor healthy control mean
      bhc_sd=apply(exp_base_sam[,grp==control_val],1,sd) ### baylor healthy control std


      module.l=c()
      bhc.lower.prop=bhc.upper.prop=bhc.sign=matrix(NA,n.module,n.base.case)
      for (i in 1:n.module)
      {
        index.i=which(modexp$Module%in%unique.module[i]) ### for ith module
        #module.l[i]=sum(index.i)
        module.l[i]=modordnum[i]
        for (j in 1:n.base.case) ### for jth case baseline
        {
          bhc.lower.prop[i,j]=sum(exp_base_case[index.i,j]<bhc_mean[index.i]-2*bhc_sd[index.i])/module.l[i]
          bhc.upper.prop[i,j]=sum(exp_base_case[index.i,j]>bhc_mean[index.i]+2*bhc_sd[index.i])/module.l[i]
          if (bhc.lower.prop[i,j]>0.1 & bhc.upper.prop[i,j]<0.1) {bhc.sign[i,j]=1-bhc.lower.prop[i,j]}
          else if (bhc.lower.prop[i,j]<0.1 & bhc.upper.prop[i,j]>0.1) {bhc.sign[i,j]=1+bhc.upper.prop[i,j]}
          else {bhc.sign[i,j]=1+bhc.upper.prop[i,j]-bhc.lower.prop[i,j]}
        }
      }

      base.case.sample=base.sample[grp!=control_val]
      #base.case.sample_s=base.case.sample[order(base.case.sample)]
      #bhc.sign_s=bhc.sign[,order(base.case.sample)]
      #colnames(bhc.sign_s)=base.case.sample_s
      colnames(bhc.sign)=base.case.sample
      #rownames(bhc.sign_s)=unique.module
      rownames(bhc.sign)=unique.module
      mod1<-bhc.sign
    }

    if(!LongitudinalTF){
      mod2<-NULL}


    if(LongitudinalTF){
      if(is.null(control_var)){
        des<-design
        if(!illumina){
          if(data_type == "rna"){
            modexp<-merge(final_expression,moduleinfo,by=c("PROBE_ID"))
            index<-which(modexp$Module=="")
            if(length(index) > 0){
              modexp<-modexp[-index,]
            }
          }
          else{
            modexp<-merge(final_expression,moduleinfo2,by=c("SYMBOL"))
            index<-which(modexp$Module=="")
            if(length(index)>0){
              modexp<-modexp[-index,]}
          }
        }
        if(illumina){
          modexp<-merge(final_expression,moduleinfo,by=c("PROBE_ID"))
          index<-which(modexp$Module=="")
          if(length(index) > 0){
            modexp<-modexp[-index,]
          }
        }

        #Creating an ordered Module List and the total number of modules
        modnames<-unique(moduleinfo$Module)[-1]
        modnum<-gsub("M","",modnames)
        modvec<-as.numeric(unlist(strsplit(modnum,".",fixed=TRUE)))
        modmat<-matrix(modvec,ncol=2,byrow=T)
        modordered<-factor(modnames[order(modmat[,1],modmat[,2])],levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE)
        modordnum<-table(factor(moduleinfo$Module,levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE))
        subjectid=patient_id

        #Order expression by timepoint and subjectid
        sort_des<-des[order(des[,time_var],des[,subjectid]),]
        data_ordered = modexp[, match(sort_des[,"columnname" ], names(modexp), nomatch=0)]
        names(data_ordered)<-sort_des[,sample_id]
        baseline_data<-data_ordered[,which(sort_des[,baseline_var]==baseline_val)]
        baseline_mean<-apply(baseline_data,1,mean)
        baseline_sd<-apply(baseline_data,1,sd)

        final_data<-as.matrix(data_ordered[,-which(names(data_ordered) %in% names(baseline_data))])
        rownames(final_data)<-modexp$Module

        SignMatrix<-matrix(rep(0,dim(final_data)[1]*dim(final_data)[2]),dim(final_data)[1],dim(final_data)[2])
        SignMatrix[final_data<(baseline_mean-2*baseline_sd)]<- -1
        SignMatrix[final_data>(baseline_mean+2*baseline_sd)]<- 1
        rownames(SignMatrix)<-rownames(final_data)
        colnames(SignMatrix)<-colnames(final_data)

        CountMatrix<-c()
        for (i in 1:length(modordered)){
          subset<-SignMatrix[which(rownames(SignMatrix) == modordered[i]),]

          if(is.vector(subset)){
            CountMatrix<-rbind(CountMatrix,subset)}

          if(is.matrix(subset)){
            if(dim(subset)[1]==0){
              CountMatrix<-rbind(CountMatrix,rep(0,dim(SignMatrix)[2]))}


            if(dim(subset)[1]>1){
              CountMatrix<-rbind(CountMatrix,apply(subset,2,sum))}
          }
        }
        rownames(CountMatrix)<-modordered
        PercentMatrix<-CountMatrix/as.vector(modordnum)
        mod2<-PercentMatrix+1
        #mod2 <- values$mod2
      }

      if(!is.null(control_var)){
        des<-design
        if(!illumina){
          if(data_type == "rna"){
            modexp<-merge(final_expression,moduleinfo,by=c("PROBE_ID"))
            index<-which(modexp$Module=="")
            if(length(index) > 0){
              modexp<-modexp[-index,]
            }
          }
          else{
            modexp<-merge(final_expression,moduleinfo2,by=c("SYMBOL"))
            index<-which(modexp$Module=="")
            if(length(index)>0){
              modexp<-modexp[-index,]}
          }
        }
        if(illumina){
          modexp<-merge(final_expression,moduleinfo,by=c("PROBE_ID"))
          index<-which(modexp$Module=="")
          if(length(index) > 0){
            modexp<-modexp[-index,]
          }
        }

        #Creating an ordered Module List and the total number of modules
        modnames<-unique(moduleinfo$Module)[-1]
        modnum<-gsub("M","",modnames)
        modvec<-as.numeric(unlist(strsplit(modnum,".",fixed=TRUE)))
        modmat<-matrix(modvec,ncol=2,byrow=T)
        modordered<-factor(modnames[order(modmat[,1],modmat[,2])],levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE)
        modordnum<-table(factor(moduleinfo$Module,levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE))
        subjectid=patient_id

        #conforming to Xuans code
        base_sample_name=des$columnname[des[,baseline_var]==baseline_val] ### get columnname of baseline samples
        base_sample_grp=des[,control_var][des[,baseline_var]==baseline_val] ### the group info of all baseline samples in design file
        exp_base_sam=modexp[,colnames(modexp)%in%base_sample_name] ### get expression data of baseline samples

        #### connect expression file with design file by key information--group, donor, timepoint
        grp=base.donorid=base.sample=c()
        for (i in 1:ncol(exp_base_sam))
        {
          grp[i]=as.character(base_sample_grp[base_sample_name==colnames(exp_base_sam)[i]])
          ### baseline sample group info in expression file
          base.donorid[i]=as.character((des[,subjectid][des[,baseline_var]==baseline_val])[base_sample_name==colnames(exp_base_sam)[i]])
          ### baseline sample donorid info in expression file
          base.sample[i]=as.character((des[,sample_id][des[,baseline_var]==baseline_val])[base_sample_name==colnames(exp_base_sam)[i]])
          ### baseline sample timepoint info in expression file
        }
        unique.module<-modordered
        n.module<-length(unique.module)
        n.case=sum(grp!=control_val)
        exp_base_case=exp_base_sam[,grp!=control_val]
        n.base.case=ncol(exp_base_case)
        bhc_mean=apply(exp_base_sam[,grp==control_val],1,mean) ### baylor healthy control mean
        bhc_sd=apply(exp_base_sam[,grp==control_val],1,sd) ### baylor healthy control std
        exp_sam=modexp[,colnames(modexp)%in%des$columnname] #only expression, not probeid etc
        column_name=colnames(exp_sam)
        grp=tp=donorid=sample_name=c()
        for (i in 1:ncol(exp_sam))
        {
          grp[i]=as.character(des[,control_var][des$columnname==column_name[i]])
          tp[i]=as.character(des[,baseline_var][des$columnname==column_name[i]])
          donorid[i]=as.character(des[,subjectid][des$columnname==column_name[i]])
          sample_name[i]=as.character(des[,sample_id][des$columnname==column_name[i]])
        }
        exp_case_sample=exp_sam[,grp!=control_val]
        n.case=ncol(exp_case_sample)
        module.l=c()
        bhc.lower.prop=bhc.upper.prop=bhc.sign=matrix(NA,n.module,n.case)
        for (i in 1:n.module)
        {
          index.i=(modexp$Module%in%unique.module[i])
          module.l[i]=modordnum[i]
          for (j in 1:n.case)
          {
            bhc.lower.prop[i,j]=sum(exp_case_sample[index.i,j]<bhc_mean[index.i]-2*bhc_sd[index.i])/module.l[i]
            bhc.upper.prop[i,j]=sum(exp_case_sample[index.i,j]>bhc_mean[index.i]+2*bhc_sd[index.i])/module.l[i]
          }
        }
        case.grp=grp[grp!=control_val]
        case.column=column_name[grp!=control_val]
        case.donor=donorid[grp!=control_val]
        case.tp=tp[grp!=control_val]
        case.wk=as.numeric(case.tp)
        case.donor.unique=unique(case.donor)
        case.sample=sample_name[grp!=control_val]
        bhc.sign=matrix(NA,n.module,n.case)

        for (i in 1:n.module)
        {
          for (j in 1:length(case.donor))
          {
            if (bhc.lower.prop[i,j]>0.1 & bhc.upper.prop[i,j]<0.1) {bhc.sign[i,j]=1-bhc.lower.prop[i,j]}
            else if (bhc.lower.prop[i,j]<0.1 & bhc.upper.prop[i,j]>0.1) {bhc.sign[i,j]=1+bhc.upper.prop[i,j]}
            else {bhc.sign[i,j]=1-bhc.lower.prop[i,j]+bhc.upper.prop[i,j]}
          }
        }
        case.sample_s=case.sample[order(case.grp,case.donor,case.wk)]
        bhc.lower.prop_s=bhc.lower.prop[,order(case.grp,case.donor,case.wk)]
        bhc.upper.prop_s=bhc.upper.prop[,order(case.grp,case.donor,case.wk)]
        bhc.sign_s=bhc.sign[,order(case.grp,case.donor,case.wk)]
        colnames(bhc.lower.prop_s)=colnames(bhc.upper.prop_s)=colnames(bhc.sign_s)=case.sample_s
        rownames(bhc.sign_s)=unique.module
        mod2<-bhc.sign_s
      }
##### Create separate with respect to baseline module map ########
    if(!is.null(control_var)){
      final_expression <- final_expression[,-which(colnames(final_expression) %in% design$columnname[which(design[,control_var] == control_val)])]
      des<-design[-which(design[,control_var] == control_val),]

      if(!illumina){
        if(data_type == "rna"){
          modexp<-merge(final_expression,moduleinfo,by=c("PROBE_ID"))
          index<-which(modexp$Module=="")
          if(length(index) > 0){
            modexp<-modexp[-index,]
          }
        }
        else{
          modexp<-merge(final_expression,moduleinfo2,by=c("SYMBOL"))
          index<-which(modexp$Module=="")
          if(length(index)>0){
            modexp<-modexp[-index,]}
        }
      }
      if(illumina){
        modexp<-merge(final_expression,moduleinfo,by=c("PROBE_ID"))
        index<-which(modexp$Module=="")
        if(length(index) > 0){
          modexp<-modexp[-index,]
        }
      }
      #Creating an ordered Module List and the total number of modules
      modnames<-unique(moduleinfo$Module)[-1]
      modnum<-gsub("M","",modnames)
      modvec<-as.numeric(unlist(strsplit(modnum,".",fixed=TRUE)))
      modmat<-matrix(modvec,ncol=2,byrow=T)
      modordered<-factor(modnames[order(modmat[,1],modmat[,2])],levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE)
      modordnum<-table(factor(moduleinfo$Module,levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE))
      subjectid=patient_id

      #Order expression by timepoint and subjectid
      sort_des<-des[order(des[,time_var],des[,subjectid]),]
      data_ordered = modexp[, match(sort_des[,"columnname" ], names(modexp), nomatch=0)]
      names(data_ordered)<-sort_des[,sample_id]
      baseline_data<-data_ordered[,which(sort_des[,baseline_var]==baseline_val)]
      baseline_mean<-apply(baseline_data,1,mean)
      baseline_sd<-apply(baseline_data,1,sd)
      final_data<-as.matrix(data_ordered[,-which(names(data_ordered) %in% names(baseline_data))])
      rownames(final_data)<-modexp$Module

      SignMatrix<-matrix(rep(0,dim(final_data)[1]*dim(final_data)[2]),dim(final_data)[1],dim(final_data)[2])
      SignMatrix[final_data<(baseline_mean-2*baseline_sd)]<- -1
      SignMatrix[final_data>(baseline_mean+2*baseline_sd)]<- 1
      rownames(SignMatrix)<-rownames(final_data)
      colnames(SignMatrix)<-colnames(final_data)

      CountMatrix<-c()
      for (i in 1:length(modordered)){
        subset<-SignMatrix[which(rownames(SignMatrix) == modordered[i]),]

        if(is.vector(subset)){
          CountMatrix<-rbind(CountMatrix,subset)}

        if(is.matrix(subset)){
          if(dim(subset)[1]==0){
            CountMatrix<-rbind(CountMatrix,rep(0,dim(SignMatrix)[2]))}


          if(dim(subset)[1]>1){
            CountMatrix<-rbind(CountMatrix,apply(subset,2,sum))}
        }
      }

      rownames(CountMatrix)<-modordered

      PercentMatrix<-CountMatrix/as.vector(modordnum)

      mod3<-PercentMatrix+1
####### module map 3 (with respect to baseline) ##########
    }
  }
}
    z = list(mod1 = mod1, mod2 = mod2, mod3 = mod3)
    return(z)
})

  ######################### Unsupervised Side Menu ###############################

output$Unsupervised <- renderMenu({
    if(is.null(values$final_expression)){
      return(strong(""))
    }

    if(is.null(values$base_mod) & is.null(values$long_mod) & is.null(values$long_mod2)){
      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) & is.null(values$final_expression) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1")))
        ))
      }

      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) & is.null(values$final_expression) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap")))
      }

      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression)){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) & is.null(values$final_expression)){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression)){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")))
        ))
      }
    }

    else{
      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) == FALSE & is.null(values$long_mod) &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1")))
        ))
      }

      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) == FALSE & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod)){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) == FALSE & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) & is.null(values$final_expression) == FALSE & is.null(values$long_mod) &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) & is.null(values$final_expression) == FALSE & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod)){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) & is.null(values$final_expression) == FALSE & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) & is.null(values$final_expression) == FALSE & is.null(values$long_mod) &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1")))
        ))
      }

      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) & is.null(values$final_expression) == FALSE & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod)){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }


      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) & is.null(values$final_expression) == FALSE & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) & is.null(values$long_mod) &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1")))
        ))
      }

      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod)){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }


      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) & is.null(values$final_expression) & is.null(values$long_mod) &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1")))
        ))
      }


      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) & is.null(values$final_expression) & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod)){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) & is.null(values$final_expression) & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) & is.null(values$long_mod) &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod)){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) == FALSE & is.null(values$long_mod) &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) == FALSE & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod)){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }

      if(is.null(module_matrices()$mod1) == FALSE & is.null(module_matrices()$mod2) == FALSE & is.null(values$final_expression) == FALSE & is.null(values$long_mod) == FALSE &
         is.null(values$base_mod) == FALSE){
        return(menuItem("Unupervised Analysis", icon = icon("th-list"), tabName = "unsupervised",
                        menuSubItem("Probe Level Heat Maps", tabName = "probeheatmap"),
                        menuItem("Geneset Maps", icon = icon("angle-double-right"), tabName = "genesetmap",
                                 menuItem("Baylor Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "modulemap1"),
                                          menuSubItem("Longitudinal", tabName = "modulemap2")),
                                 menuItem("Other Module Maps", icon = icon("angle-double-right"),
                                          menuSubItem("Baseline", tabName = "othermodmap1"),
                                          menuSubItem("Longitudinal", tabName = "othermodmap2")))
        ))
      }
    }
  })


  ######################### Baseline Module Map ###############################

  output$group_label_mod2<-renderUI({
    selectInput("LabelMod2","Select variables to label samples (additional to order variables):",choices=names(values$design),selected=c(input$responder_var,input$patient_id),multiple=T)
  })

  output$ResponderStatus <- renderUI({
    level_length <- lapply(values$design, function(x) length(levels(x)))
    level_length <- unname(unlist(level_length))
    selectInput("responderStatus", "Choose a Group for Association Test:", choices = names(values$design[level_length > 1 & level_length < length(values$design[,1])]), selected = input$responder_var)
  })

  output$TopTier<-renderUI({selectInput("top","Ordering columns: variable 1",
                                        c(names(values$design),"NA"),values$responder_var)
  })

  output$MidTier<-renderUI({selectInput("mid","Ordering columns: variable 2",
                                        c(names(values$design),"NA"), values$patient_id)
  })

  output$LowTier<-renderUI({selectInput("bottom","Ordering columns: variable 3",
                                        c(names(values$design),"NA"),"NA")
  })

  order_vars<-reactive({
    x<-c(values$responder_var, values$patient_id, "NA")
    y <- x[which(x!="NA")]

    if(!is.null(input$top) & !is.null(input$mid) & !is.null(input$bottom)){
      x<-c(input$top,input$mid,input$bottom)
      y <- x[which(x!="NA")]
      return(y)
    }
    return(y)
  })

  output$subsetModVariable<-renderUI({selectInput("subsetModVar","Variable 1 to subset heatmap:",c(names(values$design)),values$responder_var)})
  mod2values <- reactive({
    vals <- c(unique(as.character(values$design[,input$subsetModVar])))
    delete <- which(vals == "")
    if(identical(delete, integer(0))){
      vals <- vals
    }
    if(identical(delete, integer(0)) == FALSE){
      vals <- vals[-delete]
    }
    return(vals)
  })
  output$subsetModValue<-renderUI({selectInput("subsetModVal","Value(s) of Variable to subset heatmap:",mod2values(),mod2values()[1],multiple=TRUE)})
  output$subsetModVariable2<-renderUI({selectInput("subsetModVar2","Variable 2 to subset heatmap:",c(names(values$design)),values$time_var)})
  mod2values2 <- reactive({
    vals <- c(unique(as.character(values$design[,input$subsetModVar2])))
    delete <- which(vals == "")
    if(identical(delete, integer(0))){
      vals <- vals
    }
    if(identical(delete, integer(0)) == FALSE){
      vals <- vals[-delete]
    }
    return(vals)
  })
  output$subsetModValue2<-renderUI({selectInput("subsetModVal2","Value(s) of variable to subset heatmap:",mod2values2(),mod2values2()[1:2],multiple=TRUE)})

  modmap2 <-reactive({
    if(is.null(module_matrices()$mod1)){return(NULL)}

    if(input$rowselect3 == FALSE){
      data<-data.frame(cbind(rownames(module_matrices()$mod1),module_matrices()$mod1))
      names(data)=c("Module",colnames(module_matrices()$mod1))
      design1<-values$design[which(values$design[,values$sample_id]%in%colnames(module_matrices()$mod1)),]


      if(input$subsetMod){

        design1<-design1[which((design1[,as.character(input$subsetModVar)] %in% as.character(input$subsetModVal)) & (design1[,as.character(input$subsetModVar2)] %in% as.character(input$subsetModVal2))),]
        data<-data.frame(cbind(rownames(module_matrices()$mod1),module_matrices()$mod1[,match(design1[,values$sample_id], colnames(module_matrices()$mod1), nomatch=0)]))
        names(data)=c("Module",colnames(module_matrices()$mod1[,match(design1[,values$sample_id], colnames(module_matrices()$mod1), nomatch=0)]))
      }

      group_order<-order_vars()
      z<-unique(c(order_vars(),input$LabelMod2))
      color_groups<-z[which(z!="NA")]

      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"

      inside<-paste("design1$",group_order,sep="")
      des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
      design_ordered<-design1[des_order,]
      data_ordered = data[, match(design_ordered[,values$sample_id ], names(data), nomatch=0)]

      data_ordered = cbind(data.frame(data[, 1]), data.frame(data_ordered))
      names(data_ordered)[1] = "Module"

      if(input$FirstSix=="First Seven Rounds"){
        round<-97 #62
        data_ordered = data_ordered[1:round, ]
      }

      if(input$FirstSix=="Only Annotated"){
        anno_index<-which(data_ordered$Module %in% module_annotations$Module)
        data_ordered = data_ordered[anno_index, ]
      }

      groups=as.data.frame(design_ordered[,color_groups, drop = F])

      pc_clean = data_ordered[, -1]

      x = as.matrix(pc_clean)
      x=t(100*(apply(x,1,as.numeric)-1))
      num_col<-ncol(x)
      colnames(x) = design_ordered[, values$sample_id]
      ddm = as.dendrogram(fastcluster::hclust(dist(x)))
      x = x[order.dendrogram(ddm),]
      myval<-max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))]<- myval
        x[which(x==min(x))]<- -myval
      }

      if(input$MMorder_3 == TRUE){
        rownames(x) <- sub(".",",",rownames(x), fixed = TRUE)
        x <- x[gtools::mixedsort(rownames(x)),]
        rownames(x) <- sub(",",".",rownames(x), fixed = TRUE)
      }

      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(rownames(x)==module_annotations$Module[i])
        rownames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],rownames(x)[anno_index2],sep=" ")
      }

      if(input$ColClust == TRUE){
        colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      }
      if(input$ColClust == FALSE){
        colddm = NA
      }

      y = list(color_groups = color_groups, design_ordered = design_ordered, groups = groups, colddm = colddm, ddm = ddm, x = x,color_palette = color_palette)
      return(y)
    }

    if(input$rowselect3 == TRUE){
      data<-data.frame(cbind(rownames(module_matrices()$mod1),module_matrices()$mod1))
      names(data)=c("Module",colnames(module_matrices()$mod1))
      design1<-values$design[which(values$design[,values$sample_id]%in%colnames(module_matrices()$mod1)),]


      if(input$subsetMod){
        design1<-design1[which((design1[,as.character(input$subsetModVar)] %in% as.character(input$subsetModVal)) & (design1[,as.character(input$subsetModVar2)] %in% as.character(input$subsetModVal2))),]
        data<-data.frame(cbind(rownames(module_matrices()$mod1),module_matrices()$mod1[,match(design1[,values$sample_id], colnames(module_matrices()$mod1), nomatch=0)]))
        names(data)=c("Module",colnames(module_matrices()$mod1[,match(design1[,values$sample_id], colnames(module_matrices()$mod1), nomatch=0)]))
      }

      group_order<-order_vars()
      z<-unique(c(order_vars(),input$LabelMod2))
      color_groups<-z[which(z!="NA")]

      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"

      inside<-paste("design1$",group_order,sep="")
      des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
      design_ordered<-design1[des_order,]
      data_ordered = data[, match(design_ordered[,values$sample_id ], names(data), nomatch=0)]

      data_ordered = cbind(data.frame(data[, 1]), data.frame(data_ordered))
      names(data_ordered)[1] = "Module"

      groups=as.data.frame(design_ordered[,color_groups, drop = F])

      pc_clean = data_ordered[, -1]

      x = as.matrix(pc_clean)

      x=t(100*(apply(x,1,as.numeric)-1))
      num_col<-ncol(x)
      colnames(x) = design_ordered[, values$sample_id]

      v_mdnam <- read.csv(input$modsel3$datapath, header = TRUE)
      modnames <- as.character(v_mdnam[, 1])

      index44 <- match(modnames, rownames(x))

      x <- x[index44, ]

      ddm = as.dendrogram(fastcluster::hclust(dist(x)))

      if(input$MMorder_3 == FALSE){
        x = x[order.dendrogram(ddm),]
      }

      myval<-max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))]<- myval
        x[which(x==min(x))]<- -myval
      }

      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(rownames(x)==module_annotations$Module[i])
        rownames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],rownames(x)[anno_index2],sep=" ")
      }

      if(input$ColClust == TRUE){
        colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      }
      if(input$ColClust == FALSE){
        colddm = NA
      }

      y = list(color_groups = color_groups, design_ordered = design_ordered, groups = groups, colddm = colddm, ddm = ddm, x = x,color_palette = color_palette)
      return(y)
    }
  })

  opt_numClust <- reactive({
    if(is.null(module_matrices()$mod1)){return(NULL)}
    design_ordered = modmap2()$design_ordered
    x = modmap2()$x
    stdev_x <- apply(x,2,sd)
    stdev0 <- which(stdev_x == 0)
    if(length(stdev0) >0){
      x = modmap2()$x
    }
    else{
      x = apply(x, 2, function(y) (y - mean(y))/sd(y))
    }
    dist_x = dist(t(x))
    hcl = fastcluster::hclust(dist_x)
    colddm = as.dendrogram(hcl)
    d = sapply(2:round(nrow(design_ordered)/2), function(y) clValid::dunn(dist_x, cutree(hcl,y)))
    opt_num = which(d == max(d)) + 1

    y = list(hcl = hcl, opt_num = opt_num, d = d, colddm = colddm)
    return(y)

  })

  output$clusternumber1 <- renderUI({
    numericInput("clustnumber1", "Number of clusters:", min = 2, value = 2, step = 1)
  })

  output$ClusterCuts2 <- renderUI({
    numericInput('ClustCut2', "Number of clusters", min = 2, value = opt_numClust()$opt_num, step = 1)
  })

  clusterx <- reactive({
    if(is.null(module_matrices()$mod1)){return(NULL)}
    hcl = opt_numClust()$hcl
    design_ordered = modmap2()$design_ordered
    groups = modmap2()$groups

    if(input$ClusterChoice2 == TRUE){
      clusters = cutree(hcl, input$ClustCut2)
      design_ordered$Clusters = as.character(clusters)
      groups <- cbind(groups, Cluster = design_ordered$Clusters)
    }

    y = list(groups = groups)
    return(y)
  })

  modmap2_data <- reactive({
    if(is.null(module_matrices()$mod1)){return(NULL)}
    color_groups = modmap2()$color_groups
    groups = clusterx()$groups
    colddm = opt_numClust()$colddm
    ddm = modmap2()$ddm
    x = modmap2()$x
    color_palette = modmap2()$color_palette

    if(length(which(names(groups) %in% values$patient_id))){
      groups[,which(names(groups) %in% values$patient_id)]<-as.numeric(as.factor(groups[,which(names(groups) %in% values$patient_id)]))
    }

    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i]) & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    else{
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i])){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    for(i in 1:ncol(groups)){
      if(is.numeric(groups[,i]) == F){
        groups[,i] <- as.character((groups[,i]))
      }
    }

    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10 & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    else{
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    if(values$time_var %in% colnames(groups)){
      groups[,values$time_var] = as.factor(groups[,values$time_var])
    }

    n_groups = ncol(groups)
    palette_set = setPalettes(n_groups)
    factors = unlist(apply(groups, 2, function(x) as.data.frame(factor(x))), recursive=F)
    factor_levels = lapply(factors, getLevels)
    factor_lengths = lapply(factor_levels, function(x) length(x))

    color_gen = NULL
    add_list = list()

    for (i in 1:n_groups) {
      j = palette_set[[i]](as.numeric(factor_lengths[i]))
      color_gen = c(color_gen, j)
      add_list[[i]] = list(rect = list(col = "transparent", fill = j[factors[[i]]]))
    }

    annotation<-c()
    for(i in 1:n_groups){
      annotation<-cbind(annotation,add_list[[i]]$rect$fill)}

    annotation<-annotation[,n_groups:1]
    names(annotation)<-color_groups[n_groups:1]

    anno1<-colorRampPalette(c("navy", "yellow", "firebrick3"))(length(unique(groups[,1])))
    names(anno1)<-unique(groups[,1])
    eval(parse(text=paste(names(groups)[1],"=","anno1",sep="")))
    eval(parse(text=paste("first_color=list(",names(groups)[1],"=",names(groups)[1],")")))

    if(input$ColClust == TRUE){
      colddm = colddm
    }
    if(input$ColClust == FALSE){
      colddm = NA
    }

    y<-list(dat=x,groups = groups, colddm = colddm, ddm2=ddm,add_list2=add_list, first_color = first_color, color_palette2=color_palette, factor_levels = factor_levels, color_gen = color_gen)
    return(y)
  })

  plotsize0<-function(){input$modmap2size}
  plotsize1 <- function(){input$modmap2size1}
  plotresolution <- function(){input$PlotResolution}


  output$downloadModMap2 <- downloadHandler(

    filename = function() {paste(values$project_name,'_BaslineModuleMap','.csv', sep='')  },
    content = function(file) {
      write.csv(t(((modmap2_data()$dat/100)+1)), file)
    }
  )

  modtab2_data <- reactive({
    if(is.null(module_matrices()$mod1)){return(NULL)}
    design_ordered = modmap2()$design_ordered
    opt_num = opt_numClust()$opt_num
    hcl = opt_numClust()$hcl
    clusters = cutree(hcl, input$clustnumber1)
    design_ordered$Clusters<-clusters

    tab <- aggregate(design_ordered[[input$responderStatus]] ~ design_ordered$Clusters, data = design_ordered, FUN = table)
    tab <- as.data.frame(tab[,2])
    tab1 <- cbind("Cluster" = rownames(tab), tab)
    tab2 <- tab

    low_exp_count <- sapply(1:length(tab), function(y) if(sum(tab[,y]) == 0) y)
    low_exp_count <- unlist(low_exp_count)

    if(is.null(low_exp_count) == F){
      tab2 <- as.data.frame(tab[, -low_exp_count, drop = F])
    }

    tab4 <- cbind("Cluster" = rownames(tab2), tab2)
    if(ncol(tab2) > 1){
      chi_sqr <- chisq.test(tab2)
      fish_exact <- fisher.test(tab2, simulate.p.value = T, B = 10000)
      tab3 = list()
      for(i in 1:nrow(tab4)){
        tab3[[i]] = paste(tab4[i,-1], " ", "(", round(100*(tab4[i,-1]/sum(tab4[i,-1])), 2), "%",")", sep = "")
      }

      tab3 = do.call("rbind", tab3)
      tab3 <- cbind(tab4[,1,drop = F], tab3)
      names(tab3) <- names(tab4)

    }

    else{
      chi_sqr <- list()
      fish_exact <- list()
      tab3 = tab4
    }

    tab1[,1] <- as.character(tab1[,1])
    tab1_total <- rbind(data.frame(tab1[,1,drop = F], stringsAsFactors = FALSE), "Total")
    tab1_colsums <- rbind(tab1[,-1], colSums(tab1[,-1]))
    tab1 <- cbind(tab1_total, tab1_colsums)
    tab1$Total <- rowSums(tab1[,-1])

    y <- list(tab1 = tab1, tab2 = tab2, tab3 = tab3, chi_sqr = chi_sqr, fish_exact = fish_exact)
    return(y)

  })


  output$modmap2<-renderPlot({
    if(is.null(modmap2_data())){return(NULL)}
    if(min(modmap2_data()$dat) >= 0 & max(modmap2_data()$dat) > 0){
      aheatmap2(modmap2_data()$dat,Rowv=NA, circle_size = input$modmap2radius, fontsize = input$Fontsize2,annCol = modmap2_data()$groups, Colv = modmap2_data()$colddm , annColors = modmap2_data()$first_color , color = modmap2_data()$color_palette2[500:1000])
    }
    if(max(modmap2_data()$dat) <= 0 & min(modmap2_data()$dat) < 0){
      aheatmap2(modmap2_data()$dat,Rowv=NA, circle_size = input$modmap2radius, fontsize = input$Fontsize2,annCol = modmap2_data()$groups, Colv = modmap2_data()$colddm , annColors = modmap2_data()$first_color , color = modmap2_data()$color_palette2[1:500])
    }
    if(min(modmap2_data()$dat) < 0 & max(modmap2_data()$dat) > 0){
      aheatmap2(modmap2_data()$dat,Rowv=NA, circle_size = input$modmap2radius, fontsize = input$Fontsize2,annCol = modmap2_data()$groups, Colv = modmap2_data()$colddm , annColors = modmap2_data()$first_color , color = modmap2_data()$color_palette2)
    }
    if(all(modmap2_data()$dat == 0)){
      return(NULL)
    }

  },height=plotsize1,width=plotsize0)

  output$downloadModPlot2 <- downloadHandler(
    filename = function() {paste(values$project_name,'_','ModulePlot_Baseline','.png', sep = '')},
    content = function(file){
      png(file, width = (plotresolution()/72)*plotsize0(), height = (plotresolution()/72)*plotsize1(), res = plotresolution())
      if(min(modmap2_data()$dat) >= 0 & max(modmap2_data()$dat) > 0){
        print(aheatmap2(modmap2_data()$dat,Rowv=NA, circle_size = input$modmap2radius, fontsize = input$Fontsize2,annCol = modmap2_data()$groups, Colv = modmap2_data()$colddm , annColors = modmap2_data()$first_color , color = modmap2_data()$color_palette2[500:1000]))
      }
      if(max(modmap2_data()$dat) <= 0 & min(modmap2_data()$dat) < 0){
        print(aheatmap2(modmap2_data()$dat,Rowv=NA, circle_size = input$modmap2radius, fontsize = input$Fontsize2,annCol = modmap2_data()$groups, Colv = modmap2_data()$colddm , annColors = modmap2_data()$first_color , color = modmap2_data()$color_palette2[1:500]))
      }

      if(min(modmap2_data()$dat) < 0 & max(modmap2_data()$dat) > 0){
        print(aheatmap2(modmap2_data()$dat,Rowv=NA, circle_size = input$modmap2radius, fontsize = input$Fontsize2,annCol = modmap2_data()$groups, Colv = modmap2_data()$colddm , annColors = modmap2_data()$first_color , color = modmap2_data()$color_palette2))
      }

      if(all(modmap2_data()$dat == 0)){
        return(NULL)
      }
      dev.off()
    }
  )

  output$clusterplot <- renderPlot({
    if(is.null(modmap2_data())){return(NULL)}
    barplot(opt_numClust()$d, names.arg = 2:round(nrow(modmap2()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index")
  })

  output$downloadClusterPlot <- downloadHandler(
    filename = function() {paste(values$project_name,'_','Cluster_Plot_Baseline','.png', sep = '')},
    content = function(file){
      png(file, width = 800)
      print(barplot(opt_numClust()$d, names.arg = 2:round(nrow(modmap2()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index"))
      dev.off()
    }
  )

  output$OptimalNumber1 <- renderText({
    if(input$ClusterChoice2){
      paste("Optimal number of clusters =", opt_numClust()$opt_num)
    }
  })

  output$explanation <- renderText({
    if(is.null(module_matrices()$mod1)){return(NULL)}
    if(ncol(modtab2_data()$tab2) > 1){
      print("The table below is the table the tests are run on. It is the same as the table above, except the columns with zero counts have been deleted.")
    }
  })

  output$cluster_output <- renderTable({
    if(is.null(module_matrices()$mod1)){return(NULL)}
    modtab2_data()$tab1
  }, include.rownames = FALSE, digits = 0)

  output$cluster_tab <- renderTable({
    if(is.null(module_matrices()$mod1)){return(NULL)}
    if(ncol(modtab2_data()$tab2) > 1){
      modtab2_data()$tab3
    }
  }, include.rownames = FALSE)

  output$chisquare_test <- renderText({
    if(is.null(module_matrices()$mod1)){return(NULL)}
    if(ncol(modtab2_data()$tab2) > 1){

      if(modtab2_data()[[4]]$p.value < .001){
        paste("Chi_square statistic = ", round(modtab2_data()[[4]]$statistic, 2), ",", "p-value < .001")
      }

      else{
        paste("Chi-Square Test Statistic = ", round(modtab2_data()[[4]]$statistic, 2), ",", "p-value =", round(modtab2_data()[[4]]$p.value, 3))
      }
    }
  })

  output$fisher <- renderText({
    if(is.null(module_matrices()$mod1)){return(NULL)}
    if(ncol(modtab2_data()$tab2) > 1){

      if(modtab2_data()[[5]]$p.value < .001){
        paste("Fishers Exact Test: p-value < .001")
      }
      else{
        paste("Fishers Exact Test: p-value = ", round(modtab2_data()[[5]]$p.value, 3))
      }
    }
  })


  ########################### Longitudinal Module Maps ###################################

  output$group_label_mod3<-renderUI({
    selectInput("LabelMod3","Select variables to label samples (additional to order variables):",choices= c(names(values$design), "Clusters"),selected=c(input$responder_var,input$patient_id),multiple=T)
  })

  output$ResponderStatus3 <- renderUI({
    level_length <- lapply(values$design, function(x) length(levels(x)))
    level_length <- unname(unlist(level_length))
    selectInput("responderStatus3", "Choose a Group for Association Test", choices = names(values$design[level_length > 1 & level_length < length(values$design[,1])]), selected = input$responder_var)
  })

  res_stat1 <- reactive({
    input$responderStatus3
  })

  output$TopTier3<-renderUI({selectInput("top3","Ordering columns: variable 1",
                                         c(names(values$design),"NA"),values$responder_var)
  })

  output$MidTier3<-renderUI({selectInput("mid3","Ordering columns: variable 2",
                                         c(names(values$design), "NA"),values$patient_id)
  })

  output$LowTier3<-renderUI({selectInput("bottom3","Ordering columns: variable 3",
                                         c(names(values$design),"NA"),values$time_var)
  })

  order_vars3 <- reactive({
    x<-c(values$responder_var,values$patient_id,values$time_var)
    y <- x[which(x!="NA")]

    if(!is.null(input$top3) & !is.null(input$mid3) & !is.null(input$bottom3)){
      x<-c(input$top3,input$mid3,input$bottom3)
      y <- x[which(x!="NA")]
      return(y)
    }

    return(y)
  })

  output$subsetMod3Variable<-renderUI({selectInput("subsetMod3Var","Variable 1 to subset heatmap:",c(names(values$design)),values$responder_var)})
  mod3values <- reactive({
    vals <- c(unique(as.character(values$design[,input$subsetMod3Var])))
    delete <- which(vals == "")
    if(identical(delete, integer(0))){
      vals <- vals
    }
    if(identical(delete, integer(0)) == FALSE){
      vals <- vals[-delete]
    }
    return(vals)
  })
  output$subsetMod3Value<-renderUI({selectInput("subsetMod3Val","Value(s) of variable to subset heatmap:",mod3values(),mod3values()[1],multiple=TRUE)})
  output$subsetMod3Variable2<-renderUI({selectInput("subsetMod3Var2","Variable 2 to subset heatmap:",c(names(values$design)),values$time_var)})
  mod3values2 <- reactive({
    vals <- c(unique(as.character(values$design[,input$subsetMod3Var2])))
    delete <- which(vals == "")
    if(identical(delete, integer(0))){
      vals <- vals
    }
    if(identical(delete, integer(0)) == FALSE){
      vals <- vals[-delete]
    }
    return(vals)
  })
  output$subsetMod3Value2<-renderUI({selectInput("subsetMod3Val2","Value(s) of variable to subset heatmap:",mod3values2(),mod3values2()[2:3],multiple=TRUE)})

  radius <- reactive({input$modmap3radius})

  output$BaseOrHealthy1 <- renderUI({
    if(values$hc == TRUE){
      return(selectInput("BaseOrHealthy", "Normalized With Respect to Baseline or HCs:", choices = c("With Respect to Baseline", "With Respect to HCs", "Percentage Differences"), selected = "With Respect to Baseline"))
    }
    else{
      return(NULL)
    }
  })

  modmap3 <-reactive({

    if(is.null(module_matrices()$mod2)){return(NULL)}

    if(input$rowselect2 == FALSE){
      if(values$hc == TRUE){
        if(input$BaseOrHealthy == "Percentage Differences" & is.null(module_matrices()$mod1) == FALSE){
          mod1 = module_matrices()$mod1
          mod2 = module_matrices()$mod2

          samp_numbers_mod2 = which(values$design[,values$sample_id] %in% colnames(mod2))
          samp_names_mod2 = values$design[,values$sample_id][samp_numbers_mod2]
          special = values$design[,values$patient_id][samp_numbers_mod2]
          nbase = which(values$design[,values$baseline_var][which(values$design[,values$patient_id] %in% special)] != values$baseline_val)
          special1 = unique(special[order(match(samp_names_mod2, colnames(mod2)))])
          special = special[-nbase]
          special = special1[which(special1 %in% special)]
          samp_names_mod2 = samp_names_mod2[order(match(samp_names_mod2, colnames(mod2)))]
          samp_names_mod1 = unique(samp_names_mod2[which(samp_names_mod2 %in% colnames(mod1))])
          mod1 = mod1[,order(match(colnames(mod1),samp_names_mod1))]
          mod1 = mod1 - 1
          mod2 = mod2 - 1

          x = list()

          for(i in 1:ncol(mod1)){
            x[[i]] = which(colnames(mod2) %in% design[,sample_id][which(design[,patient_id] == special[i])])
          }

          for(i in 1:length(x)){
            mod2[,x[[i]]] = mod2[,x[[i]]] - mod1[,i]
          }

          mod2 = mod2[,unlist(x)]
        }

        if(input$BaseOrHealthy == "With Respect to HCs"){
          mod2 = module_matrices()$mod2
        }
        if(input$BaseOrHealthy == "With Respect to Baseline"){
          mod2 <- module_matrices()$mod3
        }
      }

      if(values$hc == FALSE){
        mod2 = module_matrices()$mod2
      }

      data<-data.frame(cbind(rownames(mod2),mod2))
      names(data)=c("Module",colnames(mod2))
      design<-values$design[which(values$design[,values$sample_id]%in%colnames(mod2)),]

      if(input$subsetMod3){

        design<-design[which((design[,as.character(input$subsetMod3Var)] %in% as.character(input$subsetMod3Val)) & (design[,as.character(input$subsetMod3Var2)] %in% as.character(input$subsetMod3Val2))),]
        data<-data.frame(cbind(rownames(mod2),mod2[,match(design[,values$sample_id], colnames(mod2), nomatch=0)]))
        names(data)=c("Module",colnames(mod2[,match(design[,values$sample_id], colnames(mod2), nomatch=0)]))
      }

      group_order<-order_vars3()
      z<-unique(c(order_vars3(),input$LabelMod3))
      color_groups<-z[which(z!="NA")]

      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      inside<-paste("design$",group_order,sep="")
      des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
      design_ordered<-design[des_order,]
      data_ordered = data[, match(design_ordered[,values$sample_id ], names(data), nomatch=0)]

      data_ordered = cbind(data.frame(data[, 1]), data.frame(data_ordered))
      names(data_ordered)[1] = "Module"

      if(input$FirstSix3=="First Seven Rounds"){
        round<-97#62
        data_ordered = data_ordered[1:round, ]
      }

      if(input$FirstSix3=="Only Annotated"){
        anno_index<-which(data_ordered$Module %in% module_annotations$Module)
        data_ordered = data_ordered[anno_index, ]
      }

      groups= as.data.frame(design_ordered[,color_groups, drop = F])
      pc_clean = data_ordered[, -1]

      x = as.matrix(pc_clean)
      if(values$hc == TRUE){
        if(input$BaseOrHealthy == "Percentage Differences"){
          x = t(100*apply(x,1,as.numeric))
        }
        else{
          x=t(100*(apply(x,1,as.numeric)-1))
        }
      }
      else{
        x=t(100*(apply(x,1,as.numeric)-1))
      }

      num_col<-ncol(x)
      colnames(x) = design_ordered[, values$sample_id]
      ddm = as.dendrogram(fastcluster::hclust(dist(x)))
      x = x[order.dendrogram(ddm),]
      myval<-max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))]<- myval
        x[which(x==min(x))]<- -myval
      }

      if(input$MMorder_2 == TRUE){
        rownames(x) <- sub(".",",",rownames(x), fixed = TRUE)
        x <- x[gtools::mixedsort(rownames(x)),]
        rownames(x) <- sub(",",".",rownames(x), fixed = TRUE)
      }

      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(rownames(x)==module_annotations$Module[i])
        rownames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],rownames(x)[anno_index2],sep=" ")
      }

      if(input$ColClust3 == TRUE){
        colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      }
      if(input$ColClust3 == FALSE){
        colddm = NA
      }
      y = list(color_groups = color_groups, design_ordered = design_ordered, groups = groups, colddm = colddm, ddm = ddm, x = x, color_palette = color_palette)
      return(y)
    }


    if(input$rowselect2 == TRUE){
      if(values$hc == TRUE){
        if(input$BaseOrHealthy == "Percentage Differences" & is.null(module_matrices()$mod1) == FALSE){
          mod1 = module_matrices()$mod1
          mod2 = module_matrices()$mod2

          samp_numbers_mod2 = which(values$design[,values$sample_id] %in% colnames(mod2))
          samp_names_mod2 = values$design[,values$sample_id][samp_numbers_mod2]
          special = values$design[,values$patient_id][samp_numbers_mod2]
          nbase = which(values$design[,values$baseline_var][which(values$design[,values$patient_id] %in% special)] != values$baseline_val)
          special1 = unique(special[order(match(samp_names_mod2, colnames(mod2)))])
          special = special[-nbase]
          special = special1[which(special1 %in% special)]
          samp_names_mod2 = samp_names_mod2[order(match(samp_names_mod2, colnames(mod2)))]
          samp_names_mod1 = unique(samp_names_mod2[which(samp_names_mod2 %in% colnames(mod1))])
          mod1 = mod1[,order(match(colnames(mod1),samp_names_mod1))]
          mod1 = mod1 - 1
          mod2 = mod2 - 1

          x = list()

          for(i in 1:ncol(mod1)){
            x[[i]] = which(colnames(mod2) %in% design[,sample_id][which(design[,patient_id] == special[i])])
          }

          for(i in 1:length(x)){
            mod2[,x[[i]]] = mod2[,x[[i]]] - mod1[,i]
          }
          mod2 = mod2[,unlist(x)]
        }

        if(input$BaseOrHealthy == "With Respect to HCs"){
          mod2 = module_matrices()$mod2
        }
        if(input$BaseOrHealthy == "With Respect to Baseline"){
          mod2 <- module_matrices()$mod3
        }
      }

      if(values$hc == FALSE){
        mod2 <- module_matrices()$mod2
      }

      data<-data.frame(cbind(rownames(mod2),mod2))
      names(data)=c("Module",colnames(mod2))
      design<-values$design[which(values$design[,values$sample_id]%in%colnames(mod2)),]

      if(input$subsetMod3){

        design<-design[which((design[,as.character(input$subsetMod3Var)] %in% as.character(input$subsetMod3Val)) & (design[,as.character(input$subsetMod3Var2)] %in% as.character(input$subsetMod3Val2))),]
        data<-data.frame(cbind(rownames(mod2),mod2[,match(design[,values$sample_id], colnames(mod2), nomatch=0)]))
        names(data)=c("Module",colnames(mod2[,match(design[,values$sample_id], colnames(mod2), nomatch=0)]))
      }

      group_order<-order_vars3()
      z<-unique(c(order_vars3(),input$LabelMod3))
      color_groups<-z[which(z!="NA")]

      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      inside<-paste("design$",group_order,sep="")
      des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
      design_ordered<-design[des_order,]
      data_ordered = data[, match(design_ordered[,values$sample_id ], names(data), nomatch=0)]

      data_ordered = cbind(data.frame(data[, 1]), data.frame(data_ordered))
      names(data_ordered)[1] = "Module"

      groups= as.data.frame(design_ordered[,color_groups, drop = F])
      pc_clean = data_ordered[, -1]

      x = as.matrix(pc_clean)
      if(values$hc == TRUE){
        if(input$BaseOrHealthy == "Percentage Differences"){
          x = t(100*apply(x,1,as.numeric))
        }
        else{
          x=t(100*(apply(x,1,as.numeric)-1))
        }
      }
      else{
        x=t(100*(apply(x,1,as.numeric)-1))
      }
      num_col<-ncol(x)
      colnames(x) = design_ordered[, values$sample_id]

      v_mdnam <- read.csv(input$modsel2$datapath, header = TRUE)
      modnames <- as.character(v_mdnam[, 1])

      index44 <- match(modnames, rownames(x))

      x <- x[index44, ]

      ddm = as.dendrogram(fastcluster::hclust(dist(x)))

      if(input$MMorder_2 == FALSE){
        x = x[order.dendrogram(ddm),]
      }

      myval<-max(c(abs(min(x)),max(x)))
      if(min(x)<0 & max(x)>0){
        x[which(x==max(x))]<- myval
        x[which(x==min(x))]<- -myval
      }

      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(rownames(x)==module_annotations$Module[i])
        rownames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],rownames(x)[anno_index2],sep=" ")
      }

      if(input$ColClust3 == TRUE){
        colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      }
      if(input$ColClust3 == FALSE){
        colddm = NA
      }
      y = list(color_groups = color_groups, design_ordered = design_ordered, groups = groups, colddm = colddm, ddm = ddm, x = x, color_palette = color_palette)
      return(y)
    }
  })

  opt_NumClust <- reactive({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    design_ordered = modmap3()$design_ordered
    x = modmap3()$x
    stdev_x <- apply(x,2,sd)
    stdev0 <- which(stdev_x == 0)
    if(length(stdev0) >0){
      x = modmap3()$x
    }
    else{
      x = apply(x, 2, function(y) (y - mean(y))/sd(y))
    }

    dist_x = dist(t(x))
    hcl = fastcluster::hclust(dist_x)
    colddm = as.dendrogram(hcl)
    d = sapply(2:round(nrow(design_ordered)/2), function(y) clValid::dunn(dist_x, cutree(hcl,y)))
    opt_num = which(d == max(d)) + 1
    y = list(hcl = hcl, opt_num = opt_num, d = d, colddm = colddm)
    return(y)

  })

  output$clusternumber2 <- renderUI({
    numericInput("clustnumber2", "Number of clusters:", min = 2, value = 2, step = 1)
  })

  output$ClusterCuts3 <- renderUI({
    numericInput('ClustCut3', "Number of clusters", min = 2, value = opt_NumClust()$opt_num, step = 1)
  })

  cluster_x <- reactive({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    hcl = opt_NumClust()$hcl
    design_ordered = modmap3()$design_ordered
    groups = modmap3()$groups
    if(input$ClusterChoice3 == TRUE){
      clusters = cutree(hcl, input$ClustCut3)
      design_ordered$Clusters = as.character(clusters)
      groups <- cbind(groups, Cluster = design_ordered$Clusters)
    }
    y = list(groups = groups)
    return(y)
  })

  modmap3_data <- reactive({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    color_groups = modmap3()$color_groups
    groups = cluster_x()$groups
    colddm = opt_NumClust()$colddm
    ddm = modmap3()$ddm
    x = modmap3()$x
    color_palette = modmap3()$color_palette

    if(length(which(names(groups) %in% values$patient_id))){
      groups[,which(names(groups) %in% values$patient_id)]<-as.numeric(as.factor(groups[,which(names(groups) %in% values$patient_id)]))
    }

    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i]) & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    else{
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i])){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    for(i in 1:ncol(groups)){
      if(is.numeric(groups[,i]) == F){
        groups[,i] <- as.character((groups[,i]))
      }
    }

    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10 & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    else{
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    if(values$time_var %in% colnames(groups)){
      groups[,values$time_var] = as.factor(groups[,values$time_var])
    }

    n_groups = ncol(groups)
    palette_set = setPalettes(n_groups)
    factors = unlist(apply(groups, 2, function(x) as.data.frame(factor(x))), recursive=F)
    factor_levels = lapply(factors, getLevels)
    factor_lengths = lapply(factor_levels, function(x) length(x))

    color_gen = NULL
    add_list = list()

    for (i in 1:n_groups) {
      j = palette_set[[i]](as.numeric(factor_lengths[i]))
      color_gen = c(color_gen, j)
      add_list[[i]] = list(rect = list(col = "transparent", fill = j[factors[[i]]]))
    }

    annotation<-c()
    for(i in 1:n_groups){
      annotation<-cbind(annotation,add_list[[i]]$rect$fill)}

    annotation<-annotation[,n_groups:1]
    names(annotation)<-color_groups[n_groups:1]

    anno1<-colorRampPalette(c("navy", "yellow", "firebrick3"))(length(unique(groups[,1])))
    names(anno1)<-unique(groups[,1])
    eval(parse(text=paste(names(groups)[1],"=","anno1",sep="")))
    eval(parse(text=paste("first_color=list(",names(groups)[1],"=",names(groups)[1],")")))

    if(input$ColClust3 == TRUE){
      colddm = colddm
    }
    if(input$ColClust3 == FALSE){
      colddm = NA
    }

    y<-list(dat=x, groups = groups, first_color = first_color, colddm = colddm,ddm2=ddm,add_list2=add_list,color_palette2=color_palette, factor_levels = factor_levels, color_gen = color_gen)
    return(y)
  })

  plotsize00<-function(){input$modmap3size}
  plotsize11 <- function(){input$modmap3size1}
  plotresolution1 <- function(){input$PlotResolution1}

  output$downloadModMap3 <- downloadHandler(

    filename = function() {paste(values$project_name,'_LongitudinalModuleMaps','.csv', sep='')  },
    content = function(file) {
      write.csv(t(((modmap3_data()$dat/100))), file)
    }
  )

  modtab3_data <- reactive({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    design_ordered = modmap3()$design_ordered
    opt_num = opt_NumClust()$opt_num
    hcl = opt_NumClust()$hcl
    clusters = cutree(hcl, input$clustnumber2)
    design_ordered$Clusters<-clusters

    tab <- aggregate(design_ordered[[res_stat1()]] ~ design_ordered$Clusters, data = design_ordered, FUN = table)
    tab <- as.data.frame(tab[,2])
    tab1 <- cbind("Cluster" = rownames(tab), tab)
    tab2 <- tab

    low_exp_count <- sapply(1:length(tab), function(y) if(sum(tab[,y]) == 0) y)
    low_exp_count <- unlist(low_exp_count)

    if(is.null(low_exp_count) == F){
      tab2 <- as.data.frame(tab[, -low_exp_count, drop = F])
    }

    tab4 <- cbind("Cluster" = rownames(tab2), tab2)
    if(ncol(tab2) > 1){
      chi_sqr <- chisq.test(tab2)
      fish_exact <- fisher.test(tab2, simulate.p.value = T, B = 10000)
      tab3 = list()
      for(i in 1:nrow(tab4)){
        tab3[[i]] = paste(tab4[i,-1], " ", "(", round(100*(tab4[i,-1]/sum(tab4[i,-1])), 2), "%",")", sep = "")
      }

      tab3 = do.call("rbind", tab3)
      tab3 <- cbind(tab4[,1,drop = F], tab3)
      names(tab3) <- names(tab4)

    }

    else{
      chi_sqr <- list()
      fish_exact <- list()
      tab3 = tab4
    }

    tab1[,1] <- as.character(tab1[,1])
    tab1_total <- rbind(data.frame(tab1[,1,drop = F], stringsAsFactors = FALSE), "Total")
    tab1_colsums <- rbind(tab1[,-1], colSums(tab1[,-1]))
    tab1 <- cbind(tab1_total, tab1_colsums)
    tab1$Total <- rowSums(tab1[,-1])

    y <- list(tab1 = tab1, tab2 = tab2, tab3 = tab3, chi_sqr = chi_sqr, fish_exact = fish_exact)
    return(y)

  })


  output$modmap3<-renderPlot({
    if(is.null(module_matrices()$mod2)){return(NULL)}

    if(min(modmap3_data()$dat) >= 0 & max(modmap3_data()$dat) > 0){
      aheatmap2(modmap3_data()$dat,Rowv=NA, circle_size = radius(), fontsize = input$Fontsize, annCol = modmap3_data()$groups, annColors = modmap3_data()$first_color, Colv = modmap3_data()$colddm , color = modmap3_data()$color_palette2[500:1000])
    }
    if(max(modmap3_data()$dat) <= 0 & min(modmap3_data()$dat) < 0){
      aheatmap2(modmap3_data()$dat,Rowv=NA, circle_size = radius(), fontsize = input$Fontsize, annCol = modmap3_data()$groups, annColors = modmap3_data()$first_color, Colv = modmap3_data()$colddm , color = modmap3_data()$color_palette2[1:500])
    }
    if(min(modmap3_data()$dat) < 0 & max(modmap3_data()$dat) > 0){
      aheatmap2(modmap3_data()$dat,Rowv=NA, circle_size = radius(), fontsize = input$Fontsize, annCol = modmap3_data()$groups, annColors = modmap3_data()$first_color, Colv = modmap3_data()$colddm , color = modmap3_data()$color_palette2)
    }
    if(all(modmap3_data()$dat == 0)){
      return(NULL)
    }
  }, height= plotsize11, width = plotsize00)

  output$downloadModPlot3 <- downloadHandler(
    filename = function() {paste(values$project_name,'_','ModulePlot_Longitudinal','.png', sep = '')},
    content = function(file){
      png(file, width = (plotresolution1()/72)*plotsize00(), height = (plotresolution1()/72)*plotsize11(), res = plotresolution1())
      if(min(modmap3_data()$dat) >= 0 & max(modmap3_data()$dat) > 0){
        print(aheatmap2(modmap3_data()$dat,Rowv=NA, circle_size = radius(), fontsize = input$Fontsize, annCol = modmap3_data()$groups, annColors = modmap3_data()$first_color, Colv = modmap3_data()$colddm , color = modmap3_data()$color_palette2[500:1000]))
      }
      if(max(modmap3_data()$dat) <= 0 & min(modmap3_data()$dat) < 0){
        print(aheatmap2(modmap3_data()$dat,Rowv=NA, circle_size = radius(), fontsize = input$Fontsize, annCol = modmap3_data()$groups, annColors = modmap3_data()$first_color, Colv = modmap3_data()$colddm , color = modmap3_data()$color_palette2[1:500]))
      }
      if(min(modmap3_data()$dat) < 0 & max(modmap3_data()$dat) > 0){
        print(aheatmap2(modmap3_data()$dat,Rowv=NA, circle_size = radius(), fontsize = input$Fontsize, annCol = modmap3_data()$groups, annColors = modmap3_data()$first_color, Colv = modmap3_data()$colddm , color = modmap3_data()$color_palette2))
      }
      if(all(modmap3_data()$dat == 0)){
        return(NULL)
      }
      dev.off()
    }
  )

  output$OptimalNumber2 <- renderText({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    if(input$ClusterChoice3){
      paste("Optimal number of clusters =", opt_NumClust()$opt_num)
    }
  })

  output$clusterplot2 <- renderPlot({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    barplot(opt_NumClust()$d, names.arg = 2:round(nrow(modmap3()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index")
  })

  output$downloadClusterPlot2 <- downloadHandler(
    filename = function() {paste(values$project_name,'_','Cluster_Plot_Longitudinal','.png', sep = '')},
    content = function(file){
      png(file, width = 800)
      print(barplot(opt_NumClust()$d, names.arg = 2:round(nrow(modmap3()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index"))
      dev.off()
    }
  )

  output$cluster_output3 <- renderTable({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    modtab3_data()$tab1
  }, include.rownames = FALSE, digits = 0)

  output$explanation2 <- renderText({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    if(ncol(modtab3_data()$tab2) > 1){
      print("The table below is the table the tests are run on. It is the same as the table above, except the columns with zero counts have been deleted.")
    }
  })

  output$cluster_tab2 <- renderTable({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    if(ncol(modtab3_data()$tab2) > 1){
      modtab3_data()$tab3
    }
  }, include.rownames = FALSE)

  output$chisquare_test3 <- renderText({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    if(ncol(modtab3_data()$tab2) > 1){
      if(modtab3_data()[[4]]$p.value < .001){
        paste("Chi_square statistic = ", round(modtab3_data()[[4]]$statistic, 2), ",", "p-value < .001")
      }
      else{
        paste("Chi-Square Test Statistic = ", round(modtab3_data()[[4]]$statistic, 2), ",", "p-value =", round(modtab3_data()[[4]]$p.value, 3))
      }
    }
  })

  output$fisher3 <- renderText({
    if(is.null(module_matrices()$mod2)){return(NULL)}
    if(ncol(modtab3_data()$tab2) > 1){

      if(modtab3_data()[[5]]$p.value < .001){
        paste("Fishers Exact Test: p-value < .001")
      }
      else{
        paste("Fishers Exact Test: p-value = ", round(modtab3_data()[[5]]$p.value, 3))
      }

    }

  })


#################################### Other Module Map Baseline #################################

  output$group_label_mod2Other <- renderUI({
    selectInput("LabelMod2Other","Select variables to label samples (additional to order variables):",choices=names(values$design),selected=c(input$responder_var,input$patient_id),multiple=T)
  })

  output$ResponderStatusOther <- renderUI({
    level_length <- lapply(values$design, function(x) length(levels(x)))
    level_length <- unname(unlist(level_length))

    selectInput("responderStatusOther", "Choose a Group for Association Test:", choices = names(values$design[level_length > 1 & level_length < length(values$design[,1])]), selected = input$responder_var)
  })

  output$TopTierOther <- renderUI({selectInput("topOther","Ordering columns: variable 1",
                                        c(names(values$design),"NA"),values$responder_var)
  })

  output$MidTierOther<-renderUI({selectInput("midOther","Ordering columns: variable 2",
                                        c(names(values$design),"NA"), values$patient_id)
  })

  output$LowTierOther<-renderUI({selectInput("bottomOther","Ordering columns: variable 3",
                                        c(names(values$design),"NA"),"NA")
  })

  order_varsOther<-reactive({
    x<-c(values$responder_var, values$patient_id, "NA")
    y <- x[which(x!="NA")]

    if(!is.null(input$topOther) & !is.null(input$midOther) & !is.null(input$bottomOther)){
      x<-c(input$topOther,input$midOther,input$bottomOther)
      y <- x[which(x!="NA")]
      return(y)
    }
    return(y)
  })

  output$subsetModVariableOther<-renderUI({selectInput("subsetModVarOther","Variable 1 to subset heatmap:",c(names(values$design)),values$responder_var)})
  mod2valuesOther <- reactive({
    vals <- c(unique(as.character(values$design[,input$subsetModVarOther])))
    delete <- which(vals == "")
    if(identical(delete, integer(0))){
      vals <- vals
    }
    if(identical(delete, integer(0)) == FALSE){
      vals <- vals[-delete]
    }
    return(vals)
  })

  output$subsetModValueOther<-renderUI({selectInput("subsetModValOther","Value(s) of Variable to subset heatmap:",mod2valuesOther(),mod2valuesOther()[1],multiple=TRUE)})
  output$subsetModVariable2Other<-renderUI({selectInput("subsetModVar2Other","Variable 2 to subset heatmap:",c(names(values$design)),values$time_var)})
  mod2values2Other <- reactive({
    vals <- c(unique(as.character(values$design[,input$subsetModVar2Other])))
    delete <- which(vals == "")
    if(identical(delete, integer(0))){
      vals <- vals
    }
    if(identical(delete, integer(0)) == FALSE){
      vals <- vals[-delete]
    }
    return(vals)
  })

  output$subsetModValue2Other<-renderUI({selectInput("subsetModVal2Other","Value(s) of variable to subset heatmap:",mod2values2Other(),mod2values2Other()[1:2],multiple=TRUE)})

  modmap2Other <-reactive({
    if(is.null(values$base_mod)){return(NULL)}

    if(input$rowselect3Other == FALSE){
      data<-data.frame(cbind(rownames(values$base_mod),values$base_mod))
      names(data)=c("Module",colnames(values$base_mod))
      design1<-values$design[which(values$design[,values$sample_id]%in%colnames(values$base_mod)),]


      if(input$subsetModOther){

        design1<-design1[which((design1[,as.character(input$subsetModVarOther)] %in% as.character(input$subsetModValOther)) & (design1[,as.character(input$subsetModVar2Other)] %in% as.character(input$subsetModVal2Other))),]
        data<-data.frame(cbind(rownames(values$base_mod),values$base_mod[,match(design1[,values$sample_id], colnames(values$base_mod), nomatch=0)]))
        names(data)=c("Module",colnames(values$base_mod[,match(design1[,values$sample_id], colnames(values$base_mod), nomatch=0)]))
      }

      group_order<-order_varsOther()
      z<-unique(c(order_varsOther(),input$LabelMod2Other))
      color_groups<-z[which(z!="NA")]

      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"

      inside<-paste("design1$",group_order,sep="")
      des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
      design_ordered<-design1[des_order,]
      data_ordered = data[, match(design_ordered[,values$sample_id ], names(data), nomatch=0)]

      data_ordered = cbind(data.frame(data[, 1]), data.frame(data_ordered))
      names(data_ordered)[1] = "Module"

      groups=as.data.frame(design_ordered[,color_groups, drop = F])

      pc_clean = data_ordered[, -1]

      x = as.matrix(pc_clean)
      x=t(100*(apply(x,1,as.numeric)-1))
      num_col<-ncol(x)
      colnames(x) = design_ordered[, values$sample_id]
      ddm = as.dendrogram(fastcluster::hclust(dist(x)))
      x = x[order.dendrogram(ddm),]
      myval<-max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))]<- myval
        x[which(x==min(x))]<- -myval
      }

      if(input$MMorder_3Other == TRUE){
        x <- x[gtools::mixedsort(rownames(x)),]
      }

      if(input$ColClustOther == TRUE){
        colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      }
      if(input$ColClustOther == FALSE){
        colddm = NA
      }

      y = list(color_groups = color_groups, design_ordered = design_ordered, groups = groups, colddm = colddm, ddm = ddm, x = x,color_palette = color_palette)
      return(y)
    }

    if(input$rowselect3Other == TRUE){
      data<-data.frame(cbind(rownames(values$base_mod),values$base_mod))
      names(data)=c("Module",colnames(values$base_mod))
      design1<-values$design[which(values$design[,values$sample_id]%in%colnames(values$base_mod)),]


      if(input$subsetModOther){
        design1<-design1[which((design1[,as.character(input$subsetModVarOther)] %in% as.character(input$subsetModValOther)) & (design1[,as.character(input$subsetModVar2Other)] %in% as.character(input$subsetModVal2Other))),]
        data<-data.frame(cbind(rownames(values$base_mod),values$base_mod[,match(design1[,values$sample_id], colnames(values$base_mod), nomatch=0)]))
        names(data)=c("Module",colnames(values$base_mod[,match(design1[,values$sample_id], colnames(values$base_mod), nomatch=0)]))
      }

      group_order<-order_varsOther()
      z<-unique(c(order_varsOther(),input$LabelMod2Other))
      color_groups<-z[which(z!="NA")]

      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"

      inside<-paste("design1$",group_order,sep="")
      des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
      design_ordered<-design1[des_order,]
      data_ordered = data[, match(design_ordered[,values$sample_id ], names(data), nomatch=0)]

      data_ordered = cbind(data.frame(data[, 1]), data.frame(data_ordered))
      names(data_ordered)[1] = "Module"

      groups=as.data.frame(design_ordered[,color_groups, drop = F])

      pc_clean = data_ordered[, -1]

      x = as.matrix(pc_clean)

      x=t(100*(apply(x,1,as.numeric)-1))
      num_col<-ncol(x)
      colnames(x) = design_ordered[, values$sample_id]

      v_mdnam <- read.csv(input$modsel3Other$datapath, header = TRUE)
      modnames <- as.character(v_mdnam[, 1])

      index44 <- match(modnames, rownames(x))

      x <- x[index44, ]

      ddm = as.dendrogram(fastcluster::hclust(dist(x)))

      if(input$MMorder_3Other == FALSE){
        x = x[order.dendrogram(ddm),]
      }

      myval<-max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))]<- myval
        x[which(x==min(x))]<- -myval
      }

      if(input$ColClustOther == TRUE){
        colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      }

      if(input$ColClustOther == FALSE){
        colddm = NA
      }

      y = list(color_groups = color_groups, design_ordered = design_ordered, groups = groups, colddm = colddm, ddm = ddm, x = x,color_palette = color_palette)
      return(y)
    }
  })

  opt_numClustOther <- reactive({
    if(is.null(values$base_mod)){return(NULL)}
    design_ordered = modmap2Other()$design_ordered
    x = modmap2Other()$x
    stdev_x <- apply(x,2,sd)
    stdev0 <- which(stdev_x == 0)
    if(length(stdev0) >0){
      x = modmap2Other()$x
    }
    else{
      x = apply(x, 2, function(y) (y - mean(y))/sd(y))
    }
    dist_x = dist(t(x))
    hcl = fastcluster::hclust(dist_x)
    colddm = as.dendrogram(hcl)
    d = sapply(2:round(nrow(design_ordered)/2), function(y) clValid::dunn(dist_x, cutree(hcl,y)))
    opt_num = which(d == max(d)) + 1
    y = list(hcl = hcl, opt_num = opt_num, d = d, colddm = colddm)
    return(y)

  })

  output$clusternumber1Other <- renderUI({
    numericInput("clustnumber1Other", "Number of clusters:", min = 2, value = 2, step = 1)
  })

  output$ClusterCuts2Other <- renderUI({
    numericInput('ClustCut2Other', "Number of clusters", min = 2, value = opt_numClustOther()$opt_num, step = 1)
  })

  clusterxOther <- reactive({
    if(is.null(values$base_mod)){return(NULL)}
    hcl = opt_numClustOther()$hcl
    design_ordered = modmap2Other()$design_ordered
    groups = modmap2Other()$groups

    if(input$ClusterChoice2Other == TRUE){
      clusters = cutree(hcl, input$ClustCut2Other)
      design_ordered$Clusters = as.character(clusters)
      groups <- cbind(groups, Cluster = design_ordered$Clusters)
    }

    y = list(groups = groups)
    return(y)
  })

  modmap2_dataOther <- reactive({
    if(is.null(values$base_mod)){return(NULL)}
    color_groups = modmap2Other()$color_groups
    groups = clusterxOther()$groups
    colddm = opt_numClustOther()$colddm
    ddm = modmap2Other()$ddm
    x = modmap2Other()$x
    color_palette = modmap2Other()$color_palette

    if(length(which(names(groups) %in% values$patient_id))){
      groups[,which(names(groups) %in% values$patient_id)]<-as.numeric(as.factor(groups[,which(names(groups) %in% values$patient_id)]))
    }

    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i]) & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    else{
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i])){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    for(i in 1:ncol(groups)){
      if(is.numeric(groups[,i]) == F){
        groups[,i] <- as.character((groups[,i]))
      }
    }

    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10 & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    else{
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    if(values$time_var %in% colnames(groups)){
      groups[,values$time_var] = as.factor(groups[,values$time_var])
    }

    n_groups = ncol(groups)
    palette_set = setPalettes(n_groups)
    factors = unlist(apply(groups, 2, function(x) as.data.frame(factor(x))), recursive=F)
    factor_levels = lapply(factors, getLevels)
    factor_lengths = lapply(factor_levels, function(x) length(x))

    color_gen = NULL
    add_list = list()

    for (i in 1:n_groups) {
      j = palette_set[[i]](as.numeric(factor_lengths[i]))
      color_gen = c(color_gen, j)
      add_list[[i]] = list(rect = list(col = "transparent", fill = j[factors[[i]]]))
    }

    annotation<-c()
    for(i in 1:n_groups){
      annotation<-cbind(annotation,add_list[[i]]$rect$fill)}

    annotation<-annotation[,n_groups:1]
    names(annotation)<-color_groups[n_groups:1]

    anno1<-colorRampPalette(c("navy", "yellow", "firebrick3"))(length(unique(groups[,1])))
    names(anno1)<-unique(groups[,1])
    eval(parse(text=paste(names(groups)[1],"=","anno1",sep="")))
    eval(parse(text=paste("first_color=list(",names(groups)[1],"=",names(groups)[1],")")))

    if(input$ColClustOther == TRUE){
      colddm = colddm
    }
    if(input$ColClustOther == FALSE){
      colddm = NA
    }

    y<-list(dat=x,groups = groups, colddm = colddm, ddm2=ddm,add_list2=add_list, first_color = first_color, color_palette2=color_palette, factor_levels = factor_levels, color_gen = color_gen)
    return(y)
  })

  plotsize0Other<-function(){input$modmap2sizeOther}
  plotsize1Other <- function(){input$modmap2size1Other}
  plotresolutionOther <- function(){input$PlotResolutionOther}


  output$downloadModMap2Other <- downloadHandler(

    filename = function() {paste(values$project_name,'_BaslineModuleMap','.csv', sep='')  },
    content = function(file) {
      write.csv(t(((modmap2_dataOther()$dat/100)+1)), file)
    }
  )

  modtab2_dataOther <- reactive({
    if(is.null(values$base_mod)){return(NULL)}
    design_ordered = modmap2Other()$design_ordered
    opt_num = opt_numClustOther()$opt_num
    hcl = opt_numClustOther()$hcl
    clusters = cutree(hcl, input$clustnumber1Other)
    design_ordered$Clusters<-clusters

    tab <- aggregate(design_ordered[[input$responderStatusOther]] ~ design_ordered$Clusters, data = design_ordered, FUN = table)
    tab <- as.data.frame(tab[,2])
    tab1 <- cbind("Cluster" = rownames(tab), tab)
    tab2 <- tab

    low_exp_count <- sapply(1:length(tab), function(y) if(sum(tab[,y]) == 0) y)
    low_exp_count <- unlist(low_exp_count)

    if(is.null(low_exp_count) == F){
      tab2 <- as.data.frame(tab[, -low_exp_count, drop = F])
    }

    tab4 <- cbind("Cluster" = rownames(tab2), tab2)
    if(ncol(tab2) > 1){
      chi_sqr <- chisq.test(tab2)
      fish_exact <- fisher.test(tab2, simulate.p.value = T, B = 10000)
      tab3 = list()
      for(i in 1:nrow(tab4)){
        tab3[[i]] = paste(tab4[i,-1], " ", "(", round(100*(tab4[i,-1]/sum(tab4[i,-1])), 2), "%",")", sep = "")
      }

      tab3 = do.call("rbind", tab3)
      tab3 <- cbind(tab4[,1,drop = F], tab3)
      names(tab3) <- names(tab4)

    }

    else{
      chi_sqr <- list()
      fish_exact <- list()
      tab3 = tab4
    }

    tab1[,1] <- as.character(tab1[,1])
    tab1_total <- rbind(data.frame(tab1[,1,drop = F], stringsAsFactors = FALSE), "Total")
    tab1_colsums <- rbind(tab1[,-1], colSums(tab1[,-1]))
    tab1 <- cbind(tab1_total, tab1_colsums)
    tab1$Total <- rowSums(tab1[,-1])

    y <- list(tab1 = tab1, tab2 = tab2, tab3 = tab3, chi_sqr = chi_sqr, fish_exact = fish_exact)
    return(y)

  })


  output$modmap2Other<-renderPlot({
    if(is.null(modmap2_dataOther())){return(NULL)}
    if(min(modmap2_dataOther()$dat) >= 0 & max(modmap2_dataOther()$dat) > 0){
      aheatmap2(modmap2_dataOther()$dat,Rowv=NA, circle_size = input$modmap2radiusOther, fontsize = input$Fontsize2Other,annCol = modmap2_dataOther()$groups, Colv = modmap2_dataOther()$colddm , annColors = modmap2_dataOther()$first_color , color = modmap2_dataOther()$color_palette2[500:1000])
    }
    if(max(modmap2_dataOther()$dat) <= 0 & min(modmap2_dataOther()$dat) < 0){
      aheatmap2(modmap2_dataOther()$dat,Rowv=NA, circle_size = input$modmap2radiusOther, fontsize = input$Fontsize2Other,annCol = modmap2_dataOther()$groups, Colv = modmap2_dataOther()$colddm , annColors = modmap2_dataOther()$first_color , color = modmap2_dataOther()$color_palette2[1:500])
    }
    if(min(modmap2_dataOther()$dat) < 0 & max(modmap2_dataOther()$dat) > 0){
      aheatmap2(modmap2_dataOther()$dat,Rowv=NA, circle_size = input$modmap2radiusOther, fontsize = input$Fontsize2Other,annCol = modmap2_dataOther()$groups, Colv = modmap2_dataOther()$colddm , annColors = modmap2_dataOther()$first_color , color = modmap2_dataOther()$color_palette2)
    }
    if(all(modmap2_dataOther()$dat == 0)){
      return(NULL)
    }

  },height=plotsize1Other,width=plotsize0Other)

  output$downloadModPlot2Other <- downloadHandler(
    filename = function() {paste(values$project_name,'_','ModulePlot_Baseline','.png', sep = '')},
    content = function(file){
      png(file, width = (plotresolutionOther()/72)*plotsize0Other(), height = (plotresolutionOther()/72)*plotsize1Other(), res = plotresolutionOther())
      if(min(modmap2_dataOther()$dat) >= 0 & max(modmap2_dataOther()$dat) > 0){
        print(aheatmap2(modmap2_dataOther()$dat,Rowv=NA, circle_size = input$modmap2radiusOther, fontsize = input$Fontsize2Other,annCol = modmap2_dataOther()$groups, Colv = modmap2_dataOther()$colddm , annColors = modmap2_dataOther()$first_color , color = modmap2_dataOther()$color_palette2[500:1000]))
      }
      if(max(modmap2_dataOther()$dat) <= 0 & min(modmap2_dataOther()$dat) < 0){
        print(aheatmap2(modmap2_dataOther()$dat,Rowv=NA, circle_size = input$modmap2radiusOther, fontsize = input$Fontsize2Other,annCol = modmap2_dataOther()$groups, Colv = modmap2_dataOther()$colddm , annColors = modmap2_dataOther()$first_color , color = modmap2_dataOther()$color_palette2[1:500]))
      }
      if(min(modmap2_dataOther()$dat) < 0 & max(modmap2_dataOther()$dat) > 0){
        print(aheatmap2(modmap2_dataOther()$dat,Rowv=NA, circle_size = input$modmap2radiusOther, fontsize = input$Fontsize2Other,annCol = modmap2_dataOther()$groups, Colv = modmap2_dataOther()$colddm , annColors = modmap2_dataOther()$first_color , color = modmap2_dataOther()$color_palette2))
      }
      if(all(modmap2_dataOther()$dat == 0)){
        return(NULL)
      }
      dev.off()
    }
  )

  output$clusterplot <- renderPlot({
    if(is.null(modmap2_dataOther())){return(NULL)}
    barplot(opt_numClustOther()$d, names.arg = 2:round(nrow(modmap2Other()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index")
  })

  output$downloadClusterPlot <- downloadHandler(
    filename = function() {paste(values$project_name,'_','Cluster_Plot_Baseline','.png', sep = '')},
    content = function(file){
      png(file, width = 800)
      print(barplot(opt_numClustOther()$d, names.arg = 2:round(nrow(modmap2Other()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index"))
      dev.off()
    }
  )

  output$OptimalNumber1Other <- renderText({
    if(input$ClusterChoice2Other){
      paste("Optimal number of clusters =", opt_numClustOther()$opt_num)
    }
  })

  output$explanationOther <- renderText({
    if(is.null(values$base_mod)){return(NULL)}
    if(ncol(modtab2_dataOther()$tab2) > 1){
      print("The table below is the table the tests are run on. It is the same as the table above, except the columns with zero counts have been deleted.")
    }
  })

  output$cluster_outputOther <- renderTable({
    if(is.null(values$base_mod)){return(NULL)}
    modtab2_dataOther()$tab1
  }, include.rownames = FALSE, digits = 0)

  output$cluster_tabOther <- renderTable({
    if(is.null(values$base_mod)){return(NULL)}
    if(ncol(modtab2_dataOther()$tab2) > 1){
      modtab2_dataOther()$tab3
    }
  }, include.rownames = FALSE)

  output$chisquare_testOther <- renderText({
    if(is.null(values$base_mod)){return(NULL)}
    if(ncol(modtab2_dataOther()$tab2) > 1){

      if(modtab2_dataOther()[[4]]$p.value < .001){
        paste("Chi_square statistic = ", round(modtab2_dataOther()[[4]]$statistic, 2), ",", "p-value < .001")
      }

      else{
        paste("Chi-Square Test Statistic = ", round(modtab2_dataOther()[[4]]$statistic, 2), ",", "p-value =", round(modtab2_dataOther()[[4]]$p.value, 3))
      }
    }
  })

  output$fisherOther <- renderText({
    if(is.null(values$base_mod)){return(NULL)}
    if(ncol(modtab2_dataOther()$tab2) > 1){

      if(modtab2_dataOther()[[5]]$p.value < .001){
        paste("Fishers Exact Test: p-value < .001")
      }
      else{
        paste("Fishers Exact Test: p-value = ", round(modtab2_dataOther()[[5]]$p.value, 3))
      }
    }
  })


  #################################### Other Module Maps Longitudinal ############################


  output$group_label_mod3Other<-renderUI({
    selectInput("LabelMod3Other","Select variables to label samples (additional to order variables):",choices= c(names(values$design), "Clusters"),selected=c(input$responder_var,input$patient_id),multiple=T)
  })

  output$ResponderStatus3Other <- renderUI({
    level_length <- lapply(values$design, function(x) length(levels(x)))
    level_length <- unname(unlist(level_length))
    selectInput("responderStatus3Other", "Choose a Group for Association Test", choices = names(values$design[level_length > 1 & level_length < length(values$design[,1])]), selected = input$responder_var)
  })

  res_stat1Other <- reactive({
    input$responderStatus3Other
  })

  output$TopTier3Other<-renderUI({selectInput("top3Other","Ordering columns: variable 1",
                                         c(names(values$design),"NA"),values$responder_var)
  })

  output$MidTier3Other<-renderUI({selectInput("mid3Other","Ordering columns: variable 2",
                                         c(names(values$design), "NA"),values$patient_id)
  })

  output$LowTier3Other<-renderUI({selectInput("bottom3Other","Ordering columns: variable 3",
                                         c(names(values$design),"NA"),values$time_var)
  })

  order_vars3Other <- reactive({
    x<-c(values$responder_var,values$patient_id,values$time_var)
    y <- x[which(x!="NA")]

    if(!is.null(input$top3Other) & !is.null(input$mid3Other) & !is.null(input$bottom3Other)){
      x<-c(input$top3Other,input$mid3Other,input$bottom3Other)
      y <- x[which(x!="NA")]
      return(y)
    }

    return(y)
  })

  output$subsetMod3VariableOther<-renderUI({selectInput("subsetMod3VarOther","Variable 1 to subset heatmap:",c(names(values$design)),values$responder_var)})
  mod3valuesOther <- reactive({
    vals <- c(unique(as.character(values$design[,input$subsetMod3VarOther])))
    delete <- which(vals == "")
    if(identical(delete, integer(0))){
      vals <- vals
    }
    if(identical(delete, integer(0)) == FALSE){
      vals <- vals[-delete]
    }
    return(vals)
  })
  output$subsetMod3ValueOther<-renderUI({selectInput("subsetMod3ValOther","Value(s) of variable to subset heatmap:",mod3valuesOther(),mod3valuesOther()[1],multiple=TRUE)})
  output$subsetMod3Variable2Other<-renderUI({selectInput("subsetMod3Var2Other","Variable 2 to subset heatmap:",c(names(values$design)),values$time_var)})
  mod3values2Other <- reactive({
    vals <- c(unique(as.character(values$design[,input$subsetMod3Var2Other])))
    delete <- which(vals == "")
    if(identical(delete, integer(0))){
      vals <- vals
    }
    if(identical(delete, integer(0)) == FALSE){
      vals <- vals[-delete]
    }
    return(vals)
  })
  output$subsetMod3Value2Other<-renderUI({selectInput("subsetMod3Val2Other","Value(s) of variable to subset heatmap:",mod3values2Other(),mod3values2Other()[2:3],multiple=TRUE)})

  radiusOther <- reactive({input$modmap3radiusOther})

  output$BaseOrHealthy1Other <- renderUI({
    if(values$hc == TRUE){
      return(selectInput("BaseOrHealthyOther", "Normalized With Respect to Baseline or HCs:", choices = c("With Respect to Baseline", "With Respect to HCs", "Percentage Differences"), selected = "With Respect to Baseline"))
    }
    else{
      return(NULL)
    }
  })

  modmap3Other <-reactive({

    if(is.null(values$long_mod)){return(NULL)}

    if(input$rowselect2Other == FALSE){
      if(values$hc == TRUE){
        if(input$BaseOrHealthyOther == "Percentage Differences" & is.null(values$base_mod) == FALSE){
          mod1 = values$base_mod
          mod2 = values$long_mod

          samp_numbers_mod2 = which(values$design[,values$sample_id] %in% colnames(mod2))
          samp_names_mod2 = values$design[,values$sample_id][samp_numbers_mod2]
          special = values$design[,values$patient_id][samp_numbers_mod2]
          nbase = which(values$design[,values$baseline_var][which(values$design[,values$patient_id] %in% special)] != values$baseline_val)
          special1 = unique(special[order(match(samp_names_mod2, colnames(mod2)))])
          special = special[-nbase]
          special = special1[which(special1 %in% special)]
          samp_names_mod2 = samp_names_mod2[order(match(samp_names_mod2, colnames(mod2)))]
          samp_names_mod1 = unique(samp_names_mod2[which(samp_names_mod2 %in% colnames(mod1))])
          mod1 = mod1[,order(match(colnames(mod1),samp_names_mod1))]
          mod1 = mod1 - 1
          mod2 = mod2 - 1

          x = list()

          for(i in 1:ncol(mod1)){
            x[[i]] = which(colnames(mod2) %in% design[,sample_id][which(design[,patient_id] == special[i])])
          }

          for(i in 1:length(x)){
            mod2[,x[[i]]] = mod2[,x[[i]]] - mod1[,i]
          }

          mod2 = mod2[,unlist(x)]
        }

        if(input$BaseOrHealthy == "With Respect to HCs"){
          mod2 = values$long_mod
        }
        if(input$BaseOrHealthy == "With Respect to Baseline"){
          mod2 <- values$long_mod2
        }
      }

      if(values$hc == FALSE){
        mod2 = values$long_mod
      }

      data<-data.frame(cbind(rownames(mod2),mod2))
      names(data)=c("Module",colnames(mod2))
      design<-values$design[which(values$design[,values$sample_id]%in%colnames(mod2)),]

      if(input$subsetMod3Other){

        design<-design[which((design[,as.character(input$subsetMod3VarOther)] %in% as.character(input$subsetMod3ValOther)) & (design[,as.character(input$subsetMod3Var2Other)] %in% as.character(input$subsetMod3Val2Other))),]
        data<-data.frame(cbind(rownames(mod2),mod2[,match(design[,values$sample_id], colnames(mod2), nomatch=0)]))
        names(data)=c("Module",colnames(mod2[,match(design[,values$sample_id], colnames(mod2), nomatch=0)]))
      }

      group_order<-order_vars3Other()
      z<-unique(c(order_vars3Other(),input$LabelMod3Other))
      color_groups<-z[which(z!="NA")]

      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      inside<-paste("design$",group_order,sep="")
      des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
      design_ordered<-design[des_order,]
      data_ordered = data[, match(design_ordered[,values$sample_id ], names(data), nomatch=0)]

      data_ordered = cbind(data.frame(data[, 1]), data.frame(data_ordered))
      names(data_ordered)[1] = "Module"

      groups= as.data.frame(design_ordered[,color_groups, drop = F])
      pc_clean = data_ordered[, -1]

      x = as.matrix(pc_clean)
      if(values$hc == TRUE){
        if(input$BaseOrHealthyOther == "Percentage Differences"){
          x = t(100*apply(x,1,as.numeric))
        }
        else{
          x=t(100*(apply(x,1,as.numeric)-1))
        }
      }
      else{
        x=t(100*(apply(x,1,as.numeric)-1))
      }
      num_col<-ncol(x)
      colnames(x) = design_ordered[, values$sample_id]
      ddm = as.dendrogram(fastcluster::hclust(dist(x)))
      x = x[order.dendrogram(ddm),]
      myval<-max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))]<- myval
        x[which(x==min(x))]<- -myval
      }

      if(input$MMorder_2Other == TRUE){
        x <- x[gtools::mixedsort(rownames(x)),]
      }


      if(input$ColClust3Other == TRUE){
        colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      }
      if(input$ColClust3Other == FALSE){
        colddm = NA
      }
      y = list(color_groups = color_groups, design_ordered = design_ordered, groups = groups, colddm = colddm, ddm = ddm, x = x, color_palette = color_palette)
      return(y)
    }

    if(input$rowselect2Other == TRUE){
      if(values$hc == TRUE){
        if(input$BaseOrHealthyOther == "Percentage Differences" & is.null(values$base_mod) == FALSE){
          mod1 = values$base_mod
          mod2 = values$long_mod

          samp_numbers_mod2 = which(values$design[,values$sample_id] %in% colnames(mod2))
          samp_names_mod2 = values$design[,values$sample_id][samp_numbers_mod2]
          special = values$design[,values$patient_id][samp_numbers_mod2]
          nbase = which(values$design[,values$baseline_var][which(values$design[,values$patient_id] %in% special)] != values$baseline_val)
          special1 = unique(special[order(match(samp_names_mod2, colnames(mod2)))])
          special = special[-nbase]
          special = special1[which(special1 %in% special)]
          samp_names_mod2 = samp_names_mod2[order(match(samp_names_mod2, colnames(mod2)))]
          samp_names_mod1 = unique(samp_names_mod2[which(samp_names_mod2 %in% colnames(mod1))])
          mod1 = mod1[,order(match(colnames(mod1),samp_names_mod1))]
          mod1 = mod1 - 1
          mod2 = mod2 - 1

          x = list()

          for(i in 1:ncol(mod1)){
            x[[i]] = which(colnames(mod2) %in% design[,sample_id][which(design[,patient_id] == special[i])])
          }

          for(i in 1:length(x)){
            mod2[,x[[i]]] = mod2[,x[[i]]] - mod1[,i]
          }
          mod2 = mod2[,unlist(x)]
        }

        if(input$BaseOrHealthyOther == "With Respect to HCs"){
          mod2 = values$long_mod
        }
        if(input$BaseOrHealthy == "With Respect to Baseline"){
          mod2 <- values$long_mod2
        }
      }

      if(values$hc == FALSE){
        mod2 <- values$long_mod
      }

      data<-data.frame(cbind(rownames(mod2),mod2))
      names(data)=c("Module",colnames(mod2))
      design<-values$design[which(values$design[,values$sample_id]%in%colnames(mod2)),]

      if(input$subsetMod3Other){

        design<-design[which((design[,as.character(input$subsetMod3VarOther)] %in% as.character(input$subsetMod3ValOther)) & (design[,as.character(input$subsetMod3Var2Other)] %in% as.character(input$subsetMod3Val2Other))),]
        data<-data.frame(cbind(rownames(mod2),mod2[,match(design[,values$sample_id], colnames(mod2), nomatch=0)]))
        names(data)=c("Module",colnames(mod2[,match(design[,values$sample_id], colnames(mod2), nomatch=0)]))
      }

      group_order<-order_vars3Other()
      z<-unique(c(order_vars3Other(),input$LabelMod3Other))
      color_groups<-z[which(z!="NA")]

      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      inside<-paste("design$",group_order,sep="")
      des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
      design_ordered<-design[des_order,]
      data_ordered = data[, match(design_ordered[,values$sample_id ], names(data), nomatch=0)]

      data_ordered = cbind(data.frame(data[, 1]), data.frame(data_ordered))
      names(data_ordered)[1] = "Module"

      groups= as.data.frame(design_ordered[,color_groups, drop = F])
      pc_clean = data_ordered[, -1]

      x = as.matrix(pc_clean)
      if(values$hc == TRUE){
        if(input$BaseOrHealthyOther == "Percentage Differences"){
          x = t(100*apply(x,1,as.numeric))
        }
        else{
          x=t(100*(apply(x,1,as.numeric)-1))
        }
      }
      else{
        x=t(100*(apply(x,1,as.numeric)-1))
      }
      num_col<-ncol(x)
      colnames(x) = design_ordered[, values$sample_id]

      v_mdnam <- read.csv(input$modsel2Other$datapath, header = TRUE)
      modnames <- as.character(v_mdnam[, 1])

      index44 <- match(modnames, rownames(x))

      x <- x[index44, ]

      ddm = as.dendrogram(fastcluster::hclust(dist(x)))

      if(input$MMorder_2Other == FALSE){
        x = x[order.dendrogram(ddm),]
      }

      myval<-max(c(abs(min(x)),max(x)))
      if(min(x)<0 & max(x)>0){
        x[which(x==max(x))]<- myval
        x[which(x==min(x))]<- -myval
      }

      if(input$ColClust3Other == TRUE){
        colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      }

      if(input$ColClust3Other == FALSE){
        colddm = NA
      }
      y = list(color_groups = color_groups, design_ordered = design_ordered, groups = groups, colddm = colddm, ddm = ddm, x = x, color_palette = color_palette)
      return(y)
    }
  })

  opt_NumClustOther <- reactive({
    if(is.null(values$long_mod)){return(NULL)}
    design_ordered = modmap3Other()$design_ordered
    x = modmap3Other()$x
    stdev_x <- apply(x,2,sd)
    stdev0 <- which(stdev_x == 0)
    if(length(stdev0) >0){
      x = modmap3Other()$x
    }
    else{
      x = apply(x, 2, function(y) (y - mean(y))/sd(y))
    }
    dist_x = dist(t(x))
    hcl = fastcluster::hclust(dist_x)
    colddm = as.dendrogram(hcl)
    d = sapply(2:round(nrow(design_ordered)/2), function(y) clValid::dunn(dist_x, cutree(hcl,y)))
    opt_num = which(d == max(d)) + 1
    y = list(hcl = hcl, opt_num = opt_num, d = d, colddm = colddm)
    return(y)

  })

  output$clusternumber2Other <- renderUI({
    numericInput("clustnumber2Other", "Number of clusters:", min = 2, value = 2, step = 1)
  })

  output$ClusterCuts3Other <- renderUI({
    numericInput('ClustCut3Other', "Number of clusters", min = 2, value = opt_NumClustOther()$opt_num, step = 1)
  })

  cluster_xOther<- reactive({
    if(is.null(values$long_mod)){return(NULL)}
    hcl = opt_NumClustOther()$hcl
    design_ordered = modmap3Other()$design_ordered
    groups = modmap3Other()$groups

    if(input$ClusterChoice3Other == TRUE){
      clusters = cutree(hcl, input$ClustCut3Other)
      design_ordered$Clusters = as.character(clusters)
      groups <- cbind(groups, Cluster = design_ordered$Clusters)
    }
    y = list(groups = groups)
    return(y)
  })

  modmap3_dataOther <- reactive({
    if(is.null(values$long_mod)){return(NULL)}
    color_groups = modmap3Other()$color_groups
    groups = cluster_xOther()$groups
    colddm = opt_NumClustOther()$colddm
    ddm = modmap3Other()$ddm
    x = modmap3Other()$x
    color_palette = modmap3Other()$color_palette

    if(length(which(names(groups) %in% values$patient_id))){
      groups[,which(names(groups) %in% values$patient_id)]<-as.numeric(as.factor(groups[,which(names(groups) %in% values$patient_id)]))
    }

    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i]) & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    else{
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i])){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    for(i in 1:ncol(groups)){
      if(is.numeric(groups[,i]) == F){
        groups[,i] <- as.character((groups[,i]))
      }
    }

    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10 & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    else{
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    if(values$time_var %in% colnames(groups)){
      groups[,values$time_var] = as.factor(groups[,values$time_var])
    }

    n_groups = ncol(groups)
    palette_set = setPalettes(n_groups)
    factors = unlist(apply(groups, 2, function(x) as.data.frame(factor(x))), recursive=F)
    factor_levels = lapply(factors, getLevels)
    factor_lengths = lapply(factor_levels, function(x) length(x))

    color_gen = NULL
    add_list = list()

    for (i in 1:n_groups) {
      j = palette_set[[i]](as.numeric(factor_lengths[i]))
      color_gen = c(color_gen, j)
      add_list[[i]] = list(rect = list(col = "transparent", fill = j[factors[[i]]]))
    }

    annotation<-c()
    for(i in 1:n_groups){
      annotation<-cbind(annotation,add_list[[i]]$rect$fill)}

    annotation<-annotation[,n_groups:1]
    names(annotation)<-color_groups[n_groups:1]

    anno1<-colorRampPalette(c("navy", "yellow", "firebrick3"))(length(unique(groups[,1])))
    names(anno1)<-unique(groups[,1])
    eval(parse(text=paste(names(groups)[1],"=","anno1",sep="")))
    eval(parse(text=paste("first_color=list(",names(groups)[1],"=",names(groups)[1],")")))

    if(input$ColClust3Other == TRUE){
      colddm = colddm
    }
    if(input$ColClust3Other == FALSE){
      colddm = NA
    }

    y<-list(dat=x, groups = groups, first_color = first_color, colddm = colddm,ddm2=ddm,add_list2=add_list,color_palette2=color_palette, factor_levels = factor_levels, color_gen = color_gen)
    return(y)
  })

  plotsize00Other<-function(){input$modmap3sizeOther}
  plotsize11Other <- function(){input$modmap3size1Other}
  plotresolution1Other <- function(){input$PlotResolution1Other}

  output$downloadModMap3Other <- downloadHandler(

    filename = function() {paste(values$project_name,'_LongitudinalModuleMaps','.csv', sep='')  },
    content = function(file) {
      write.csv(t(((modmap3_dataOther()$dat/100))), file)
    }
  )

  modtab3_dataOther <- reactive({
    if(is.null(values$long_mod)){return(NULL)}
    design_ordered = modmap3Other()$design_ordered
    opt_num = opt_NumClustOther()$opt_num
    hcl = opt_NumClustOther()$hcl
    clusters = cutree(hcl, input$clustnumber2Other)
    design_ordered$Clusters<-clusters
    tab <- aggregate(design_ordered[[res_stat1Other()]] ~ design_ordered$Clusters, data = design_ordered, FUN = table)
    tab <- as.data.frame(tab[,2])
    tab1 <- cbind("Cluster" = rownames(tab), tab)
    tab2 <- tab
    low_exp_count <- sapply(1:length(tab), function(y) if(sum(tab[,y]) == 0) y)
    low_exp_count <- unlist(low_exp_count)

    if(is.null(low_exp_count) == F){
      tab2 <- as.data.frame(tab[, -low_exp_count, drop = F])
    }

    tab4 <- cbind("Cluster" = rownames(tab2), tab2)
    if(ncol(tab2) > 1){
      chi_sqr <- chisq.test(tab2)
      fish_exact <- fisher.test(tab2, simulate.p.value = T, B = 10000)
      tab3 = list()
      for(i in 1:nrow(tab4)){
        tab3[[i]] = paste(tab4[i,-1], " ", "(", round(100*(tab4[i,-1]/sum(tab4[i,-1])), 2), "%",")", sep = "")
      }

      tab3 = do.call("rbind", tab3)
      tab3 <- cbind(tab4[,1,drop = F], tab3)
      names(tab3) <- names(tab4)

    }

    else{
      chi_sqr <- list()
      fish_exact <- list()
      tab3 = tab4
    }

    tab1[,1] <- as.character(tab1[,1])
    tab1_total <- rbind(data.frame(tab1[,1,drop = F], stringsAsFactors = FALSE), "Total")
    tab1_colsums <- rbind(tab1[,-1], colSums(tab1[,-1]))
    tab1 <- cbind(tab1_total, tab1_colsums)
    tab1$Total <- rowSums(tab1[,-1])

    y <- list(tab1 = tab1, tab2 = tab2, tab3 = tab3, chi_sqr = chi_sqr, fish_exact = fish_exact)
    return(y)

  })


  output$modmap3Other<-renderPlot({
    if(is.null(values$long_mod)){return(NULL)}
    if(min(modmap3_dataOther()$dat) >= 0 & max(modmap3_dataOther()$dat) > 0){
      aheatmap2(modmap3_dataOther()$dat,Rowv=NA, circle_size = radiusOther(), fontsize = input$FontsizeOther, annCol = modmap3_dataOther()$groups, annColors = modmap3_dataOther()$first_color, Colv = modmap3_dataOther()$colddm , color = modmap3_dataOther()$color_palette2[500:1000])
    }
    if(max(modmap3_dataOther()$dat) <= 0 & min(modmap3_dataOther()$dat) < 0){
      aheatmap2(modmap3_dataOther()$dat,Rowv=NA, circle_size = radiusOther(), fontsize = input$FontsizeOther, annCol = modmap3_dataOther()$groups, annColors = modmap3_dataOther()$first_color, Colv = modmap3_dataOther()$colddm , color = modmap3_dataOther()$color_palette2[1:500])
    }
    if(min(modmap3_dataOther()$dat) < 0 & max(modmap3_dataOther()$dat) > 0){
      aheatmap2(modmap3_dataOther()$dat,Rowv=NA, circle_size = radiusOther(), fontsize = input$FontsizeOther, annCol = modmap3_dataOther()$groups, annColors = modmap3_dataOther()$first_color, Colv = modmap3_dataOther()$colddm , color = modmap3_dataOther()$color_palette2)
    }
    if(all(modmap3_dataOther()$dat == 0)){
      return(NULL)
    }
  }, height= plotsize11Other, width = plotsize00Other)

  output$downloadModPlot3Other <- downloadHandler(
    filename = function() {paste(values$project_name,'_','ModulePlot_Longitudinal','.png', sep = '')},
    content = function(file){
      png(file, width = (plotresolution1Other()/72)*plotsize00Other(), height = (plotresolution1Other()/72)*plotsize11Other(), res = plotresolution1Other())
      if(min(modmap3_data()$dat) >= 0 & max(modmap3_data()$dat) > 0){
        print(aheatmap2(modmap3_dataOther()$dat,Rowv=NA, circle_size = radiusOther(), fontsize = input$FontsizeOther, annCol = modmap3_dataOther()$groups, annColors = modmap3_dataOther()$first_color, Colv = modmap3_dataOther()$colddm , color = modmap3_dataOther()$color_palette2[500:1000]))
      }
      if(max(modmap3_dataOther()$dat) <= 0 & min(modmap3_dataOther()$dat) < 0){
        print(aheatmap2(modmap3_dataOther()$dat,Rowv=NA, circle_size = radiusOther(), fontsize = input$FontsizeOther, annCol = modmap3_dataOther()$groups, annColors = modmap3_dataOther()$first_color, Colv = modmap3_dataOther()$colddm , color = modmap3_dataOther()$color_palette2[1:500]))
      }
      if(min(modmap3_dataOther()$dat) < 0 & max(modmap3_dataOther()$dat) > 0){
        print(aheatmap2(modmap3_dataOther()$dat,Rowv=NA, circle_size = radiusOther(), fontsize = input$FontsizeOther, annCol = modmap3_dataOther()$groups, annColors = modmap3_dataOther()$first_color, Colv = modmap3_dataOther()$colddm , color = modmap3_dataOther()$color_palette2))
      }
      if(all(modmap3_dataOther()$dat == 0)){
        return(NULL)
      }
      dev.off()
    }
  )

  output$OptimalNumber2Other <- renderText({
    if(is.null(values$long_mod)){return(NULL)}
    if(input$ClusterChoice3Other){
      paste("Optimal number of clusters =", opt_NumClustOther()$opt_num)
    }
  })

  output$clusterplot2Other <- renderPlot({
    if(is.null(values$long_mod)){return(NULL)}
    barplot(opt_NumClustOther()$d, names.arg = 2:round(nrow(modmap3Other()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index")
  })

  output$downloadClusterPlot2Other <- downloadHandler(
    filename = function() {paste(values$project_name,'_','Cluster_Plot_Longitudinal','.png', sep = '')},
    content = function(file){
      png(file, width = 800)
      print(barplot(opt_NumClustOther()$d, names.arg = 2:round(nrow(modmap3Other()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index"))
      dev.off()
    }
  )

  output$cluster_output3Other <- renderTable({
    if(is.null(values$long_mod)){return(NULL)}
    modtab3_dataOther()$tab1
  }, include.rownames = FALSE, digits = 0)

  output$explanation2Other <- renderText({
    if(is.null(values$long_mod)){return(NULL)}
    if(ncol(modtab3_dataOther()$tab2) > 1){
      print("The table below is the table the tests are run on. It is the same as the table above, except the columns with zero counts have been deleted.")
    }
  })

  output$cluster_tab2Other <- renderTable({
    if(is.null(values$long_mod)){return(NULL)}
    if(ncol(modtab3_dataOther()$tab2) > 1){
      modtab3_dataOther()$tab3
    }
  }, include.rownames = FALSE)

  output$chisquare_test3Other <- renderText({
    if(is.null(values$long_mod)){return(NULL)}
    if(ncol(modtab3_dataOther()$tab2) > 1){
      if(modtab3_dataOther()[[4]]$p.value < .001){
        paste("Chi_square statistic = ", round(modtab3_dataOther()[[4]]$statistic, 2), ",", "p-value < .001")
      }
      else{
        paste("Chi-Square Test Statistic = ", round(modtab3_dataOther()[[4]]$statistic, 2), ",", "p-value =", round(modtab3_dataOther()[[4]]$p.value, 3))
      }
    }
  })

  output$fisher3Other <- renderText({
    if(is.null(values$long_mod)){return(NULL)}
    if(ncol(modtab3_dataOther()$tab2) > 1){

      if(modtab3_dataOther()[[5]]$p.value < .001){
        paste("Fishers Exact Test: p-value < .001")
      }
      else{
        paste("Fishers Exact Test: p-value = ", round(modtab3_dataOther()[[5]]$p.value, 3))
      }

    }

  })



  ##################################### Probe Level Heat Map #####################################

  output$ResponderStatus4 <- renderUI({
    level_length <- lapply(values$design, function(x) length(levels(x)))
    level_length <- unname(unlist(level_length))
    selectInput("responderStatus4", "Choose a Group for Association Test", choices = names(values$design[level_length > 1 & level_length < length(values$design[,1])]), selected = input$responder_var)
  })

  res_stat <- reactive({
    input$responderStatus4
  })

  output$group_label_probe<-renderUI({
    selectInput("Labelprobe","Select variables to label samples (additional to order variables):",choices=names(values$design),selected=c(input$responder_var,input$patient_id),multiple=T)
  })

  output$TopTierProbe<-renderUI({selectInput("topprobe","Ordering columns: variable 1",
                                             c(names(values$design),"NA"),values$responder_var)
  })

  output$MidTierProbe<-renderUI({selectInput("midprobe","Ordering columns: variable 2",
                                             c(names(values$design),"NA"),values$patient_id)
  })

  output$LowTierProbe<-renderUI({selectInput("bottomprobe","Ordering columns: variable 3",
                                             c(names(values$design),"NA"),"NA")
  })

  order_varsProbe<- eventReactive(input$go,{
    x<-c(values$responder_var, values$patient_id, "NA")
    y <- x[which(x!="NA")]

    if(!is.null(input$topprobe) & !is.null(input$midprobe) & !is.null(input$bottomprobe)){
      x<-c(input$topprobe,input$midprobe,input$bottomprobe)
      y <- x[which(x!="NA")]
      return(y)
    }
    return(y)
  })

  output$subsetProbeVariable<-renderUI({selectInput("subsetProbeVar","Variable 1 to subset heatmap:",c(names(values$design)), values$responder_var)})
  output$subsetProbeValue<-renderUI({selectInput("subsetProbeVal","Value(s) of Variable to subset heatmap:",c(unique(as.character(values$design[,input$subsetProbeVar]))),c(unique(as.character(values$design[,input$subsetProbeVar])))[1],multiple=TRUE)})
  output$subsetProbeVariable2<-renderUI({selectInput("subsetProbeVar2","Variable 2 to subset heatmap:",c(names(values$design)), values$time_var)})
  output$subsetProbeValue2<-renderUI({selectInput("subsetProbeVal2","Value(s) of Variable to subset heatmap:",c(unique(as.character(values$design[,input$subsetProbeVar2]))),c(unique(as.character(values$design[,input$subsetProbeVar2])))[1],multiple=TRUE)})

  heatmapdata<-reactive({
    if (input$set==1){#heatmapbase1
      base_sample_name=values$design$columnname[values$design[,baseline_var]==values$baseline_val]
      ind1<-which(colnames(values$final_expression)%in%c("PROBE_ID","SYMBOL"))
      ind2<-which(colnames(values$final_expression)%in%base_sample_name)
      if(input$uploadprobes == FALSE){
        exp_base_sam=values$final_expression[,c(ind1,ind2)]
        des_base_sam=values$design[which(values$design$columnname%in%colnames(exp_base_sam)),]
        y<-data.manipulate(exp=exp_base_sam,des=des_base_sam,values$baseline_val,longitudinal=FALSE,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=TRUE)
        ddm<-values$h1b_rowdendro
        colddm<-values$h1b_coldendro
      }

      if(input$uploadprobes == TRUE){
        probes <- read.csv(input$probe_select$datapath, header = TRUE)
        probes.vector <- as.character(probes[,1])
        exp_base_sam=values$final_expression[which(values$final_expression[,names(probes)] %in% probes.vector), c(ind1,ind2)]
        des_base_sam=values$design[which(values$design$columnname%in%colnames(exp_base_sam)),]
        y<-data.manipulate(exp=exp_base_sam,des=des_base_sam,values$baseline_val,longitudinal=FALSE,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=TRUE)
        ddm<-values$h1b_rowdendro
        colddm<-values$h1b_coldendro
      }
    }
    if (input$set==2){#heatmapbase2
      base_sample_name=values$design$columnname[values$design[,baseline_var]==values$baseline_val]
      ind1<-which(colnames(values$final_expression)%in%c("PROBE_ID","SYMBOL"))
      ind2<-which(colnames(values$final_expression)%in%base_sample_name)

      if(input$uploadprobes == FALSE){
        exp_base_sam=values$final_expression[,c(ind1,ind2)]
        des_base_sam=values$design[which(values$design$columnname%in%colnames(exp_base_sam)),]
        y<-data.manipulate(exp=exp_base_sam,des=des_base_sam,values$control_var,values$control_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=FALSE)
      }

      if(input$uploadprobes == TRUE){
        probes <- read.csv(input$probe_select$datapath, header = TRUE)
        probes.vector <- as.character(probes[,1])
        exp_base_sam=values$final_expression[which(values$final_expression[,names(probes)] %in% probes.vector), c(ind1,ind2)]
        des_base_sam=values$design[which(values$design$columnname%in%colnames(exp_base_sam)),]
        y<-data.manipulate(exp=exp_base_sam,des=des_base_sam,values$control_var,values$control_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=FALSE)
      }

      ddm<-values$h2b_rowdendro
      colddm<-values$h2b_coldendro

    }
    if (input$set==3){#heatmap1

      if(input$uploadprobes == FALSE){
        y<-data.manipulate(exp=values$final_expression,des=values$design[which(values$design$columnname %in% colnames(values$final_expression)),],values$baseline_var,values$baseline_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=TRUE)
        ddm<-values$h1_rowdendro
        colddm<-values$h1_coldendro
      }

      if(input$uploadprobes == TRUE){
        probes <- read.csv(input$probe_select$datapath, header = TRUE)
        probes.vector <- as.character(probes[,1])
        y<-data.manipulate(exp=values$final_expression[which(values$final_expression[,names(probes)] %in% probes.vector),],des=values$design[which(values$design$columnname %in% colnames(values$final_expression)),],values$baseline_var,values$baseline_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=TRUE)
        ddm<-values$h1_rowdendro
        colddm<-values$h1_coldendro
      }

    }
    if (input$set==4){#heatmap2

      if(input$uploadprobes == FALSE){
        y<-data.manipulate(exp=values$final_expression,des=values$design[which(values$design$columnname %in% colnames(values$final_expression)),],values$control_var,values$control_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=FALSE)
        ddm<-values$h2_rowdendro
        colddm<-values$h2_coldendro
      }

      if(input$uploadprobes == TRUE){
        probes <- read.csv(input$probe_select$datapath, header = TRUE)
        probes.vector <- as.character(probes[,1])
        y<-data.manipulate(exp=values$final_expression[which(values$final_expression[,names(probes)] %in% probes.vector), ],des=values$design[which(values$design$columnname %in% colnames(values$final_expression)),],values$control_var,values$control_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=FALSE)
        ddm<-values$h2_rowdendro
        colddm<-values$h2_coldendro
      }
    }
    if(input$set==5){#heatmap3

      if(input$uploadprobes == FALSE){
        if(values$hc==TRUE){
          des_w_controls<-values$design[which(values$design$columnname %in% colnames(values$final_expression)),]
          des_wo_controls<-values$design[-which(values$design[,values$control_var]==values$control_val),]
          h5index<-c(1,2,which(colnames(values$final_expression) %in% des_wo_controls$columnname))
          y<-data.manipulate(exp=values$final_expression[,h5index],des=des_wo_controls,values$baseline_var,values$baseline_val,longitudinal=TRUE,subjects=values$patient_id,lg2=FALSE,keepbase=FALSE,format="Probes",allsamples=FALSE)
        }

        if(values$hc==FALSE){
          y<-data.manipulate(exp=values$final_expression,des=values$design[which(values$design$columnname %in% colnames(values$final_expression)),],values$baseline_var,values$baseline_val,longitudinal=TRUE,subjects=values$patient_id,lg2=FALSE,keepbase=FALSE,format="Probes",allsamples=FALSE)
        }
        ddm<-values$h3_rowdendro
        colddm<-values$h3_coldendro
      }

      if(input$uploadprobes == TRUE){
        probes <- read.csv(input$probe_select$datapath, header = TRUE)
        probes.vector <- as.character(probes[,1])

        if(values$hc==TRUE){
          des_w_controls<-values$design[which(values$design$columnname %in% colnames(values$final_expression)),]
          des_wo_controls<-values$design[-which(values$design[,values$control_var]==values$control_val),]
          h5index<-c(1,2,which(colnames(values$final_expression) %in% des_wo_controls$columnname))
          y<-data.manipulate(exp=values$final_expression[which(values$final_expression[,names(probes)] %in% probes.vector),h5index],des=des_wo_controls,values$baseline_var,values$baseline_val,longitudinal=TRUE,subjects=values$patient_id,lg2=FALSE,keepbase=FALSE,format="Probes",allsamples=FALSE)
        }

        if(values$hc==FALSE){
          y<-data.manipulate(exp=values$final_expression[which(values$final_expression[,names(probes)] %in% probes.vector),],des=values$design[which(values$design$columnname %in% colnames(values$final_expression)),],values$baseline_var,values$baseline_val,longitudinal=TRUE,subjects=values$patient_id,lg2=FALSE,keepbase=FALSE,format="Probes",allsamples=FALSE)
        }

        ddm<-values$h3_rowdendro
        colddm<-values$h3_coldendro
      }
    }
    z<-list(y=y,ddm=ddm,colddm=colddm)
    return(z)})

  heatmapname<-reactive({
    if(input$set==1) heattxt<-"Baseline Median Normalized"
    if(input$set==2) heattxt<-"Baseline Healthy Normalized"
    if(input$set==3) heattxt<-"All Samples Median Normalized"
    if(input$set==4) heattxt<-"All Samples Healthy Normalized"
    if(input$set==5) heattxt<-"All Samples Normalized to each Subjects Baseline"
    heattxt
  })

  Allprobes <- eventReactive(input$go,{input$uploadprobes})

  Width <- eventReactive(input$go,{input$HeatMapSize1})
  Height <- eventReactive(input$go, {input$HeatMapSize2})
  plot_width <- function(){Width()}
  plot_height <- function(){Height()}
  Legend_Size <- eventReactive(input$go, {input$LegendSize})
  plotresolution2 <- function(){input$PlotResolution2}
  tree_height <- eventReactive(input$go, {input$TreeHeight})
  font_size <- eventReactive(input$go,{input$FontSize})
  row.clust <- eventReactive(input$go,{
    if(input$row_cluster == TRUE){
      rowclust <- TRUE
    }
    if(input$row_cluster == FALSE){
      rowclust <- NA
    }
    return(rowclust)
  })

  heatmap_order <- reactive({
    y<-heatmapdata()$y

    if(input$subsetProbe){

      y$heatdes<-y$heatdes[which((y$heatdes[,as.character(input$subsetProbeVar)] %in% as.character(input$subsetProbeVal)) & (y$heatdes[,as.character(input$subsetProbeVar2)] %in% as.character(input$subsetProbeVal2))),]
      y$heatexp<-y$heatexp[,match(y$heatdes[,"columnname"], colnames(y$heatexp), nomatch=0)]

    }

    group_order<-order_varsProbe()
    z<-unique(c(order_varsProbe(),input$Labelprobe))
    color_groups<-z[which(z!="NA")]
    inside<-paste("y$heatdes$",group_order,sep="")
    des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
    design_ordered<-y$heatdes[des_order,]
    groups= as.data.frame(design_ordered[,color_groups, drop = F])
    z = list(y = y, groups = groups, color_groups = color_groups, design_ordered = design_ordered)
    return(z)
  })

  expression_matrix <- eventReactive(input$go,{
    probeids <- values$final_expression[,1,drop = FALSE]
    symb <- values$final_expression[,2,drop = FALSE]
    y = heatmap_order()$y
    design_ordered = heatmap_order()$design_ordered
    ddm = heatmapdata()$ddm

    x = y$heatexp[, match(design_ordered[,"columnname" ], colnames(y$heatexp), nomatch=0)]
    y$heatexp <- y$heatexp[, match(design_ordered[,"columnname" ], colnames(y$heatexp), nomatch=0)]
    colnames(y$heatexp) <- design_ordered[,values$sample_id]
    colnames(x)<-design_ordered[,values$sample_id]

    if(input$setcutoff!=0){
      cut1<-as.numeric(input$setcutoff)
      x[x>cut1]<-cut1
      x[x<(-cut1)]<--cut1}

    if(input$uploadprobes == TRUE){
      if(input$row_cluster == TRUE){
        ddm = as.dendrogram(fastcluster::hclust(dist(x)))
        x <- x[order.dendrogram(ddm),]
        y$heatexp <- y$heatexp[order.dendrogram(ddm),]
        probeids <- probeids[order.dendrogram(ddm),,drop = FALSE]
        symb <- symb[order.dendrogram(ddm),,drop = FALSE]
      }
      if(input$row_cluster == FALSE){
        ddm <- NA
        probeids <- probeids[,,drop = FALSE]
        symb <- symb[,,drop = FALSE]
      }
    }

    if(input$uploadprobes == FALSE){
      x=x[order.dendrogram(ddm),]
      y$heatexp <- y$heatexp[order.dendrogram(ddm),]
      probeids <- probeids[order.dendrogram(ddm),,drop = FALSE]
      symb <- symb[order.dendrogram(ddm),,drop = FALSE]
    }
    z = list(x = x, probeids = probeids, y = y, symb = symb, ddm = ddm)
    return(z)
  })


  opt_numclust2 <- eventReactive(input$go, {
    x = expression_matrix()$x
    design_ordered = heatmap_order()$design_ordered
    colddm = NA
    stdev_x <- apply(x,2,sd)
    stdev0 <- which(stdev_x == 0)
    if(length(stdev0) >0){
      x = modmap3()$x
    }
    else{
      x = apply(x, 2, function(y) (y - mean(y))/sd(y))
    }

    dist_x = dist(t(x))
    hcl = fastcluster::hclust(dist_x)
    if(input$ColClustProbe == TRUE){
      colddm <- as.dendrogram(hcl)
    }
    d = sapply(2:round(nrow(design_ordered)/2), function(y) clValid::dunn(dist_x, cutree(hcl,y)))
    opt_num = which(d == max(d)) + 1
    z = list(opt_num = opt_num, colddm = colddm, hcl = hcl, d = d)
    return(z)
  })

  output$clusternumber <- renderUI({
    numericInput("clustnumber", "Number of clusters:", min = 2, value = 2, step = 1)
  })

  opt_NumClust2 <- reactive({
    opt_num = opt_numclust2()$opt_num
    colddm = opt_numclust2()$colddm
    hcl = opt_numclust2()$hcl
    d = opt_numclust2()$d
    z = list(opt_num = opt_num, colddm = colddm, hcl = hcl, d = d)
    return(z)
  })


  output$ClusterCuts4 <- renderUI({
    numericInput('ClustCut4', "Number of clusters", min = 2, value = opt_NumClust2()$opt_num, step = 1)
  })

  cluster_x2 <- eventReactive(input$go, {
    groups = heatmap_order()$groups
    design_ordered = heatmap_order()$design_ordered
    hcl = opt_NumClust2()$hcl

    if(input$ClusterChoice4 == TRUE){
      clusters = cutree(hcl, input$ClustCut4)
      design_ordered$Clusters = as.character(clusters)
      groups <- cbind(groups, Cluster = design_ordered$Clusters)
      z = list(groups = groups, design_ordered = design_ordered, clusters = clusters)
    }
    else{
      z = list(groups = groups, design_ordered = design_ordered)
    }
    return(z)
  })

  heatmap_colors <- eventReactive(input$go,{
    groups = cluster_x2()$groups
    color_groups = heatmap_order()$color_groups

    if(length(which(names(groups) %in% values$patient_id))){
      groups[,which(names(groups) %in% values$patient_id)]<-as.numeric(as.factor(groups[,which(names(groups) %in% values$patient_id)]))
    }

    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i]) & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    else{
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i])){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    for(i in 1:ncol(groups)){
      if(is.numeric(groups[,i]) == F){
        groups[,i] <- as.character((groups[,i]))
      }
    }

    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10 & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    else{
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }

    if(values$time_var %in% colnames(groups)){
      groups[,values$time_var] = as.factor(groups[,values$time_var])
    }

    n_groups = ncol(groups)
    palette_set = setPalettes(n_groups)
    factors = unlist(apply(groups, 2, function(x) as.data.frame(factor(x))), recursive=F)
    factor_levels = lapply(factors, getLevels)
    factor_lengths = lapply(factor_levels, function(x) length(x))

    color_gen = NULL
    add_list = list()

    for (i in 1:n_groups) {
      j = palette_set[[i]](as.numeric(factor_lengths[i]))
      color_gen = c(color_gen, j)
      add_list[[i]] = list(rect = list(col = "transparent", fill = j[factors[[i]]]))
    }

    annotation<-c()
    for(i in 1:n_groups){
      annotation<-cbind(annotation,add_list[[i]]$rect$fill)
    }

    annotation<-annotation[,n_groups:1]
    names(annotation)<-color_groups[n_groups:1]
    anno1<-colorRampPalette(c("navy", "yellow", "firebrick3"))(length(unique(groups[,1])))
    names(anno1)<-unique(groups[,1])
    eval(parse(text=paste(names(groups)[1],"=","anno1",sep="")))
    eval(parse(text=paste("first_color=list(",names(groups)[1],"=",names(groups)[1],")")))
    z <- list(groups = groups, first_color = first_color)
    return(z)
  })


  gen_clustTab <- reactive({
    if(is.null(opt_NumClust2())){return(NULL)}
    design_ordered = heatmap_order()$design_ordered
    opt_num = opt_NumClust2()$opt_num
    hcl = opt_NumClust2()$hcl
    clusters = cutree(hcl, input$clustnumber)
    design_ordered$Clusters = clusters
    tab <- aggregate(design_ordered[[res_stat()]] ~ design_ordered$Clusters, data = design_ordered, FUN = table)
    tab <- as.data.frame(tab[,2])
    tab1 <- cbind("Cluster" = rownames(tab), tab)
    tab2 <- tab
    low_exp_count <- sapply(1:length(tab), function(y) if(sum(tab[,y]) == 0) y)
    low_exp_count <- unlist(low_exp_count)

    if(is.null(low_exp_count) == F){
      tab2 <- as.data.frame(tab[, -low_exp_count, drop = F])
    }

    tab4 <- cbind("Cluster" = rownames(tab2), tab2)
    if(ncol(tab2) > 1){
      chi_sqr <- chisq.test(tab2)
      fish_exact <- fisher.test(tab2, simulate.p.value = T, B = 10000)
      tab3 = list()
      for(i in 1:nrow(tab4)){
        tab3[[i]] = paste(tab4[i,-1], " ", "(", round(100*(tab4[i,-1]/sum(tab4[i,-1])), 2), "%",")", sep = "")
      }
      tab3 = do.call("rbind", tab3)
      tab3 <- cbind(tab4[,1,drop = F], tab3)
      names(tab3) <- names(tab4)
    }

    else{
      chi_sqr <- list()
      fish_exact <- list()
      tab3 = tab4
    }

    tab1[,1] <- as.character(tab1[,1])
    tab1_total <- rbind(data.frame(tab1[,1,drop = F], stringsAsFactors = FALSE), "Total")
    tab1_colsums <- rbind(tab1[,-1], colSums(tab1[,-1]))
    tab1 <- cbind(tab1_total, tab1_colsums)
    tab1$Total <- rowSums(tab1[,-1])
    y <- list(tab1 = tab1, tab2 = tab2, tab3 = tab3, chi_sqr = chi_sqr, fish_exact = fish_exact)
    return(y)

  })

  output$heatmap <- renderPlot({
    if(Allprobes() == FALSE){
      withProgress(message = 'Making plot',
                   detail = 'This may take a while...', value = 1,{
                     aheatmap2(expression_matrix()$x,Rowv=NA, border_color = "grey60", Colv=opt_NumClust2()$colddm, treeheight = tree_height(), fontsize = font_size(), cexRow=1.2, annheight = Legend_Size(),color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),annCol = heatmap_colors()$groups,annColors= heatmap_colors()$first_color,labRow = NA, breaks=0)
                   })
    }

    if(Allprobes() == TRUE){
      withProgress(message = 'Making plot',
                   detail = 'This may take a while...', value = 1,{
                     aheatmap2(expression_matrix()$x,Rowv=row.clust(), border_color = "grey60", Colv=opt_NumClust2()$colddm, treeheight = tree_height(), fontsize = font_size(), cexRow=1.2, annheight = Legend_Size(),color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),annCol = heatmap_colors()$groups,annColors= heatmap_colors()$first_color,breaks=0)
                   })
    }

  }, height = plot_height, width = plot_width)

  heatmap_download <- reactive({

    ddm <- expression_matrix()$ddm

    if(input$uploadprobes == FALSE){
      heatmapdata <- cbind(expression_matrix()$symb, expression_matrix()$y$heatexp)
      heatmapdata <- cbind(expression_matrix()$probeids, heatmapdata)
    }

    if(input$uploadprobes == TRUE){
      probes <- read.csv(input$probe_select$datapath, header = TRUE)
      probes.vector <- as.character(probes[,1])
      exp.symb <- values$final_expression[which(values$final_expression[,names(probes)] %in% probes.vector),2,drop = FALSE]
      exp.probes <- values$final_expression[which(values$final_expression[,names(probes)] %in% probes.vector),1,drop = FALSE]
      if(input$row_cluster == TRUE){
        exp.symb <- exp.symb[order.dendrogram(ddm),,drop = FALSE]
        exp.probes <- exp.probes[order.dendrogram(ddm),,drop = FALSE]
      }
      heatmapdata <- cbind(exp.symb, expression_matrix()$y$heatexp)
      heatmapdata <- cbind(exp.probes, heatmapdata)
    }

    if(input$modulemeans == TRUE){
      if(is.null(values$illumina) == FALSE){
        if(values$illumina == TRUE){
          heatmapdata <- heatmapdata[which(heatmapdata$PROBE_ID %in% moduleinfo1$PROBE_ID),]
          mod.info <- moduleinfo1[which(moduleinfo1$PROBE_ID %in% heatmapdata$PROBE_ID),]
          heatmapdata <- heatmapdata[match(mod.info$PROBE_ID, heatmapdata$PROBE_ID, nomatch = 0),]
          heatmapdata$Module <- mod.info$Module
        }
        else{
          if(is.null(values$moduleinfo2)){
            heatmapdata <- heatmapdata[which(heatmapdata$PROBE_ID %in% moduleinfo2$SYMBOL),]
            mod.info <- moduleinfo2[which(moduleinfo2$SYMBOL %in% heatmapdata$PROBE_ID),]
            heatmapdata <- heatmapdata[match(mod.info$SYMBOL, heatmapdata$PROBE_ID, nomatch = 0),]
            heatmapdata$Module <- mod.info$Module
          }
          else{
            heatmapdata <- heatmapdata[which(heatmapdata$PROBE_ID %in% values$moduleinfo2$SYMBOL),]
            mod.info <- values$moduleinfo2[which(values$moduleinfo2$SYMBOL %in% heatmapdata$PROBE_ID),]
            heatmapdata <- heatmapdata[match(mod.info$SYMBOL, heatmapdata$PROBE_ID, nomatch = 0),]
            heatmapdata$Module <- mod.info$Module
          }
        }
      }

      if(is.null(values$illumina)){
        heatmapdata <- heatmapdata[which(heatmapdata$PROBE_ID %in% moduleinfo1$PROBE_ID),]
        mod.info <- moduleinfo1[which(moduleinfo1$PROBE_ID %in% heatmapdata$PROBE_ID),]
        heatmapdata <- heatmapdata[match(mod.info$PROBE_ID, heatmapdata$PROBE_ID, nomatch = 0),]
        heatmapdata$Module <- mod.info$Module
      }

      num_col <- ncol(heatmapdata)-1
      new.expression <- list()

      for(i in 3:num_col){
        new.expression[[i-2]] <- aggregate(heatmapdata[,i] ~ Module, data = heatmapdata, mean)
      }

      n <- length(new.expression)
      for(i in 2:n){
        new.expression[[i]] <- data.frame(new.expression[[i]][,-1])
      }

      new.expression <- do.call(cbind, new.expression)
      colnames(new.expression)[-1] <- colnames(heatmapdata)[3:num_col]
      heatmapdata <- new.expression
    }

    return(heatmapdata)

  })

  output$downloadHeatmap <- downloadHandler(

    filename = function() {paste(values$project_name,'_',heatmapname(),'.csv', sep='')  },
    content = function(file) {
      write.csv(heatmap_download(), file, row.names = FALSE)
    }
  )

  output$downloadModPlot4 <- downloadHandler(
    filename = function() {paste('ProbeLevel_Heatmap','.png', sep = '')},
    content = function(file){
      png(file, width = (plotresolution2()/72)*plot_width(), height = (plotresolution2()/72)*plot_height(), res = plotresolution2())
      print(aheatmap2(expression_matrix()$x,Rowv=NA,Colv=opt_NumClust2()$colddm, treeheight = input$TreeHeight, fontsize = input$FontSize, cexRow=1.2, annheight = Legend_Size(),color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),annCol = heatmap_colors()$groups,annColors= heatmap_colors()$first_color,labRow=NA,breaks=0))
      dev.off()
    }
  )

  output$clusterplot3 <- renderPlot({
    barplot(opt_NumClust2()$d, names.arg = 2:round(nrow(heatmap_order()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index")
  })

  output$downloadClusterPlot3 <- downloadHandler(
    filename = function() {paste('Cluster_Plot_ProbeLevel','.png', sep = '')},
    content = function(file){
      png(file, width = 800)
      print(barplot(opt_NumClust2()$d, names.arg = 2:round(nrow(heatmap_order()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index"))
      dev.off()
    }
  )

  output$OptimalNumber <- renderText({
    if(input$ClusterChoice4 == TRUE){
      paste("Optimal number of clusters =", opt_NumClust2()$opt_num )
    }
  })

  output$cluster_output4 <- renderTable({
    if(is.null(gen_clustTab())){return(NULL)}
    gen_clustTab()$tab1
  }, include.rownames = FALSE, digits = 0)

  output$explanation3 <- renderText({
    if(is.null(gen_clustTab())){return(NULL)}
    if(ncol(gen_clustTab()$tab2) > 1){
      print("The table below is the table the tests are run on. It is the same as the table above, except the columns with zero counts have been deleted.")
    }
  })

  output$cluster_tab3 <- renderTable({
    if(is.null(gen_clustTab())){return(NULL)}
    if(ncol(gen_clustTab()$tab2) > 1){
      gen_clustTab()$tab3
    }
  }, include.rownames = FALSE)

  output$chisquare_test4 <- renderText({
    if(is.null(gen_clustTab())){return(NULL)}
    if(ncol(gen_clustTab()$tab2) > 1){
      if(gen_clustTab()[[4]]$p.value < .001){
        paste("Chi_square statistic = ", round(gen_clustTab()[[4]]$statistic, 2), ",", "p-value < .001")
      }
      else{
        paste("Chi-Square Test Statistic = ", round(gen_clustTab()[[4]]$statistic, 2), ",", "p-value =", round(gen_clustTab()[[4]]$p.value, 3))
      }
    }
  })

  output$fisher4 <- renderText({
    if(is.null(gen_clustTab())){return(NULL)}
    if(ncol(gen_clustTab()$tab2) > 1){
      if(gen_clustTab()[[5]]$p.value < .001){
        paste("Fishers Exact Test: p-value < .001")
      }
      else{
        paste("Fishers Exact Test: p-value = ", round(gen_clustTab()[[5]]$p.value, 3))
      }
    }
  })


  ################################# DGE ###########################################

  output$diffge <- renderMenu({
    if(is.null(values$results_file)){
      return(strong(""))
    }
    else{
      return(menuItem("DGE", icon = icon("th-list"), tabName = "dge",
                      menuSubItem("General Info", tabName = "generalinfo"),
                      menuSubItem("Overview", tabName = "overview"),
                      menuSubItem("Gene List Maker", tabName = "genelistmaker"),
                      menuSubItem("Gene Search", tabName = "genesearch")))
    }
  })

  output$mytabs <- renderUI({
    if(is.null(values$ModulesTF) == FALSE){
      if(values$ModulesTF == TRUE){
        return(tabsetPanel(id="start",

                    tabPanel(style="height: 80vh; overflow: auto","DGE: Gene Lists", downloadButton("download_modmap_comparison", "Download Figure"), uiOutput('modmap_comparison2'),#plotOutput("modmap_comparison"),
                             helpText("Right click on hyperlinks to open in new window"),
                             downloadButton('downloadData', 'Download Table'),
                             dataTableOutput("genelisttable")),
                    tabPanel(style="height: 80vh; overflow: auto","DGE: Diagnostics", plotOutput("distplot"), plotOutput("genelistgraph"), textOutput("intro"),
                             dataTableOutput("numtable")),
                    tabPanel("Significant Probe Level Heat Map",style = "overflow: auto; height: 80vh",
                             tags$style(type="text/css", ".tab-content { overflow: visible !important; }"),
                             infoPopup("Help Message", 'The heat map below is constructed on individual samples for a number of scenarios,
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
                                       placement = "bottom", trigger = "click"),
                             helpText(""),
                             downloadButton('downloadHeatmap1', 'Download Data'),
                             downloadButton('downloadHeatmap2', 'Download Figure'),
                             plotOutput("heatmap1"),
                             plotOutput('heatmap_select2')),
                    tabPanel("Cluster Number Diagnostics",
                             downloadButton('downloadClusterPlot4', 'Download Figure'),
                             plotOutput('clusterplot4')),
                    tabPanel(style="height: 80vh; overflow: auto","Cluster Association Analysis",
                             tableOutput('cluster_output5'),
                             helpText("The Chi-square test below is used to test for association between the cluster groups and user specified variable.
                                      The Chi-square test is ideal when expected cell counts are large (expected values greater than 5).
                                      Because the test only makes sense for 2x2 contingency tables and above, test statistics are not shown
                                      when the data is subsetted on only one value."),
                             helpText(""),
                             helpText("Fisher's exact test is also shown. This test is ideal when cell counts are small. Like the Chi-square test,
                                      when the data is subsetted on only one value, test statistics are not shown."),
                             helpText(""),
                             textOutput('explanation4'),
                             tableOutput('cluster_tab4'),
                             textOutput('chisquare_test5'),
                             textOutput('fisher5')),
                    tabPanel(style="height: 80vh; overflow: auto","DGE: Module Analysis Overview",
                              downloadButton('downloadLMMModMap', 'Download Data'),
                              downloadButton('downloadLMMModMap2', 'Download Figure'),
                              uiOutput("LMMmodule_text"),
                              plotOutput("LMMmodule"),
                              plotOutput("LMMmoduleplus")),
                    tabPanel(style="height: 80vh; overflow: auto","Venn Diagram",
                             downloadButton('downloadVennPic', 'Download Figure'),
                             plotOutput("vennDiagram"),
                             downloadButton('downloadVennData', 'Download Data'),
                             dataTableOutput("venn.intersection"))
                             )
               )
      }
      else{
        return(tabsetPanel(id="start",
                    tabPanel(style="height: 80vh; overflow: auto","DGE: Gene Lists", downloadButton("download_modmap_comparison", "Download Figure"), uiOutput('modmap_comparison2'),#plotOutput("modmap_comparison"),
                             helpText("Right click on hyperlinks to open in new window"),
                             downloadButton('downloadData', 'Download Table'),
                             dataTableOutput("genelisttable")),
                    tabPanel(style="height: 80vh; overflow: auto","DGE: Diagnostics", plotOutput("distplot"), plotOutput("genelistgraph"), textOutput("intro"),
                             dataTableOutput("numtable")),
                    tabPanel("Significant Probe Level Heat Map",style = "overflow: auto; height: 80vh",
                             tags$style(type="text/css", ".tab-content { overflow: visible !important; }"),
                             infoPopup("Help Message", 'The heat map below is constructed on individual samples for a number of scenarios,
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
                                       placement = "bottom", trigger = "click"),
                             helpText(""),
                             downloadButton('downloadHeatmap1', 'Download Data'),
                             downloadButton('downloadHeatmap2', 'Download Figure'),
                             plotOutput("heatmap1"),
                             plotOutput('heatmap_select2')),
                    tabPanel("Cluster Number Diagnostics",
                             downloadButton('downloadClusterPlot4', 'Download Figure'),
                             plotOutput('clusterplot4')),
                    tabPanel(style="height: 80vh; overflow: auto","Cluster Association Analysis",
                             tableOutput('cluster_output5'),
                             helpText("The Chi-square test below is used to test for association between the cluster groups and user specified variable.
                                      The Chi-square test is ideal when expected cell counts are large (expected values greater than 5).
                                      Because the test only makes sense for 2x2 contingency tables and above, test statistics are not shown
                                      when the data is subsetted on only one value."),
                             helpText(""),
                             helpText("Fisher's exact test is also shown. This test is ideal when cell counts are small. Like the Chi-square test,
                                      when the data is subsetted on only one value, test statistics are not shown."),
                             helpText(""),
                             textOutput('explanation4'),
                             tableOutput('cluster_tab4'),
                             textOutput('chisquare_test5'),
                             textOutput('fisher5')),
                    tabPanel(style="height: 80vh; overflow: auto","Venn Diagram",
                             downloadButton('downloadVennPic', 'Download Figure'),
                             plotOutput("vennDiagram"),
                             downloadButton('downloadVennData', 'Download Data'),
                             dataTableOutput("venn.intersection"))
                             )
               )
      }
    }
    else{
      return(
        tabsetPanel(id="start",
                          tabPanel(style="height: 80vh; overflow: auto","DGE: Gene Lists", downloadButton("download_modmap_comparison", "Download Figure"), uiOutput('modmap_comparison2'),#plotOutput("modmap_comparison"),
                                    helpText("Right click on hyperlinks to open in new window"),
                                    downloadButton('downloadData', 'Download Table'),
                                    dataTableOutput("genelisttable")),
                           tabPanel(style="height: 80vh; overflow: auto","DGE: Diagnostics", plotOutput("distplot"), plotOutput("genelistgraph"), textOutput("intro"),
                                    dataTableOutput("numtable")),
                           tabPanel(style="height: 80vh; overflow: auto","Significant Probe Level Heat Map",
                                    tags$style(type="text/css", ".tab-content { overflow: visible !important; }"),
                                    infoPopup("Help Message", 'The heat map below is constructed on individual samples for a number of scenarios,
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
                                       placement = "bottom", trigger = "click"),
                             helpText(""),
                             downloadButton('downloadHeatmap1', 'Download Data'),
                             downloadButton('downloadHeatmap2', 'Download Figure'),
                             div(style = "height: 70vh; overflow: auto", plotOutput("heatmap1")),
                             plotOutput('heatmap_select2')),
                  tabPanel("Cluster Number Diagnostics",
                           downloadButton('downloadClusterPlot4', 'Download Figure'),
                           plotOutput('clusterplot4')),
                  tabPanel(style="height: 80vh; overflow: auto","Cluster Association Analysis",
                           tableOutput('cluster_output5'),
                           helpText("The Chi-square test below is used to test for association between the cluster groups and user specified variable.
                                    The Chi-square test is ideal when expected cell counts are large (expected values greater than 5).
                                    Because the test only makes sense for 2x2 contingency tables and above, test statistics are not shown
                                    when the data is subsetted on only one value."),
                           helpText(""),
                           helpText("Fisher's exact test is also shown. This test is ideal when cell counts are small. Like the Chi-square test,
                                    when the data is subsetted on only one value, test statistics are not shown."),
                           helpText(""),
                           textOutput('explanation4'),
                           tableOutput('cluster_tab4'),
                           textOutput('chisquare_test5'),
                           textOutput('fisher5')),
                  tabPanel(style="height: 80vh; overflow: auto","DGE: Module Analysis Overview",
                            downloadButton('downloadLMMModMap', 'Download Data'),
                            downloadButton('downloadLMMModMap2', 'Download Figure'),
                            uiOutput("LMMmodule_text"),
                            plotOutput("LMMmodule"),
                            plotOutput("LMMmoduleplus")),
                  tabPanel(style="height: 80vh; overflow: auto","Venn Diagram",
                           downloadButton('downloadVennPic', 'Download Figure'),
                           plotOutput("vennDiagram"),
                           downloadButton('downloadVennData', 'Download Data'),
                           dataTableOutput("venn.intersection"))
                           )
        )
    }
  })

  observeEvent(input$graphics4,{
    if(input$graphics4 == TRUE){
      removeUI(selector =  "div:has(> #heatmap1)")
    }
  })

  observeEvent(input$graphics4,{
    if(input$graphics4 == TRUE){
      insertUI(selector = "#downloadHeatmap2", where = "afterEnd", ui = div(style = "height: 130vh; overflow: auto", plotOutput("heatmap1")))
    }
  })

  observeEvent(input$graphics4,{
    if(input$graphics4 == FALSE){
      removeUI(selector = "div:has(> #heatmap1)")
    }
  })

  observeEvent(input$graphics4,{
    if(input$graphics4 == FALSE){
      insertUI(selector = "#downloadHeatmap2", where = "afterEnd", ui = div(style = "height: 70vh; overflow: auto", plotOutput("heatmap1")))
    }
  })

  dataset <- reactive({
    ind_fct <- grep("factor", sapply(values$results_file[, -c(1, 2)], class))
    if(sum(ind_fct) <= 1) return(values$results_file)
    testdat <- values$results_file
    testdat[, ind_fct + 2] <- apply(testdat[, ind_fct + 2], 2, function(x) as.numeric(gsub("<.0001", 0, x)))
    testdat
  })

  siglist <- reactive({
    results_file <- dataset()
    sigcomp0 <- results_file[grep(("PROBE_ID|Estimate of|Test.statistic|P.Value"), names(results_file))]
    sigcomp1 <- sigcomp0[grep("Estimate of", names(sigcomp0))]
    sigcomp2 <- sigcomp0[grep("P.Value", names(sigcomp0))]
    sigcomp3 <- data.frame(apply(sigcomp2, 2, p.adjust, method="fdr"))
    colnames(sigcomp3) <- gsub("P.Value","FDR.P.Value",colnames(sigcomp3))
    sigcomp4 <- data.frame(apply(sigcomp2, 2, p.adjust, method="bonferroni"))
    colnames(sigcomp4) <- gsub("P.Value","Bonf.P.Value",colnames(sigcomp4))
    if(input$showfc == TRUE){
      if(input$sign == "Both"){
        sigcomp5 <- data.frame(sigcomp1 >= input$fcval | sigcomp1 <= -input$fcval)
      }
      if(input$sign == "+"){
        sigcomp5 <- data.frame(sigcomp1 >= input$fcval)
      }
      if(input$sign == "-"){
        sigcomp5 <- data.frame(sigcomp1 <= -input$fcval)
      }
    }
    TF_Raw <- data.frame(sigcomp2 <= input$alphalevel2)
    TF_FDR <- data.frame(sigcomp3 <= input$alphalevel2)
    TF_Bonf <- data.frame(sigcomp4 <= input$alphalevel2)
    TF_Raw1 <- data.frame(sigcomp2 <= input$alphalevel1)
    TF_FDR1 <- data.frame(sigcomp3 <= input$alphalevel1)
    TF_Bonf1 <- data.frame(sigcomp4 <= input$alphalevel1)
    TF_p <- data.frame(sigcomp1 > 0)
    TF_pRaw <- TF_Raw + TF_p
    TF_pFDR <- TF_FDR + TF_p
    TF_pBonf <- TF_Bonf + TF_p
    TF_n <- data.frame(sigcomp1 < 0)
    TF_nRaw <- TF_Raw + TF_n
    TF_nFDR <- TF_FDR + TF_n
    TF_nBonf <- TF_Bonf + TF_n
    if(input$sigsign == "All (include 0)"){
      Raw <- apply(sigcomp2, 2, function(x) sum(x[!is.na(x)] < input$alphalevel2))
      names(Raw) <- gsub("P.Value for ", "", names(Raw))
      FDR <- apply(sigcomp3, 2, function(x) sum(x[!is.na(x)] < input$alphalevel2))
      names(FDR) <- gsub("FDR.P.Value.for.", "", names(FDR))
      Bonf <- apply(sigcomp4, 2, function(x) sum(x[!is.na(x)] < input$alphalevel2))
      names(Bonf) <- gsub("Bonf.P.Value.for.", "", names(Bonf))
      sigtab <- rbind(Raw, FDR, Bonf)
    }

    if(input$sigsign == "+"){
      Raw <- apply(TF_pRaw, 2, function(x) sum(x[!is.na(x)] == 2))
      names(Raw) <- gsub("P.Value for ", "", names(Raw))
      FDR <- apply(TF_pFDR, 2, function(x) sum(x[!is.na(x)] == 2))
      names(FDR) <- gsub("FDR.P.Value.for.", "", names(FDR))
      Bonf <- apply(TF_pBonf, 2, function(x) sum(x[!is.na(x)] == 2))
      names(Bonf) <- gsub("Bonf.P.Value.for.", "", names(Bonf))
      sigtab <- rbind(Raw, FDR, Bonf)
    }

    if(input$sigsign == "-"){
      Raw <- apply(TF_nRaw, 2, function(x) sum(x[!is.na(x)] == 2))
      names(Raw) <- gsub("P.Value for ", "", names(Raw))
      FDR <- apply(TF_nFDR, 2, function(x) sum(x[!is.na(x)] == 2))
      names(FDR) <- gsub("FDR.P.Value.for.", "", names(FDR))
      Bonf <- apply(TF_nBonf, 2, function(x) sum(x[!is.na(x)] == 2))
      names(Bonf) <- gsub("Bonf.P.Value.for.", "", names(Bonf))
      sigtab <- rbind(Raw, FDR, Bonf)
    }
    x <- t(sigtab)
    y <- cbind(Comparison = rownames(x), x)
    rownames(y) <- NULL
    z <- data.frame(y)
    for(i in 2:4){
      z[, i] <- as.numeric(levels(z[, i])[z[, i]])
    }
    if(input$showfc == TRUE){
      y <- list(z = z, TF_Raw = TF_Raw, TF_FDR = TF_FDR, TF_Bonf = TF_Bonf, TF_Raw1 = TF_Raw1, TF_FDR1 = TF_FDR1, TF_Bonf1 = TF_Bonf1, sigcomp5 = sigcomp5)
    }
    else{
      y <- list(z = z, TF_Raw = TF_Raw, TF_FDR = TF_FDR, TF_Bonf = TF_Bonf, TF_Raw1 = TF_Raw1, TF_FDR1 = TF_FDR1, TF_Bonf1 = TF_Bonf1)
    }
    return(y)
  })

  output$sigcomptable <- renderDataTable({
    withProgress(message = 'Making the table',
                 detail = 'This may take a while...', value = 1,{
                   siglist()$z[order(siglist()$z$Raw,decreasing=TRUE),]
                 })
  })

  output$MixedModelText <- renderText({
    paste(values$MixedModelDescription)
  })

  index <- reactive({
    which(names(dataset())==pcomp())
  })

  plottitle1 <- reactive({
    paste("Distribution of Raw p-values for", pcomp())
  })

  output$distplot<-renderPlot({
    y<-max(hist(dataset()[, index()])$density)
    hist(dataset()[,index()], freq=FALSE, xlim=c(0,1), ylim=c(0,y),
         main=paste("Distribution of Raw p-values for", substring(pcomp(),21)),
         xlab="Raw p-value's", ylab="Density")
    lines(c(0,1),c(1,1),lwd=2,lty=2)
  })

  plotdata <-reactive({
    mymatrix<-c()
    index2<-c(0.001,0.01,1:19/20)
    for (i in 1:21){mymatrix<-rbind(mymatrix,c(mycorrection(dataset()[,index()],index2[i],"RAW"),mycorrection(dataset()[,index()],index2[i],"FDR"),mycorrection(dataset()[,index()],index2[i],"BONF")))}
    mydata<-data.frame(cbind(index2,mymatrix))
    names(mydata)<-c("alpha","raw","fdr","bonf")
    data.frame(mydata)
  })

  output$numtable <-renderDataTable({
    numtab<-cbind(plotdata()[,1],plotdata()[,2:4]*dim(dataset())[1])
    names(numtab)<-c("Alpha", "Raw", "FDR", "Bonf")
    numtab[, 2:dim(numtab)[2]] = round(numtab[, 1:dim(numtab)[2]], digits = 0)
    numtab
  })

  pcomp <- reactive({
    x <- sub("^", "P.Value for ", input$comparison)
    x
  })

  genelist<-reactive({
    if(is.null(values$illumina) == FALSE){
      if(values$illumina == FALSE){
        if(!is.null(values$moduleinfo2)){
          moduleinfo <- values$moduleinfo2
        }
        ann <- list()
        moduleinfo$Modulev2_Annotation <- c("")
        for(i in 1:length(module_annotations$Module)){
          ann[[i]] <- which(moduleinfo$Module %in% module_annotations$Module[i])
        }
        for(i in 1:length(ann)){
          moduleinfo$Modulev2_Annotation[ann[[i]]] <- as.character(module_annotations$Modulev2_Annotation[i])
        }
      }
    }

    x <- dataset()
    x <- x[,c(1,2,grep(input$comparison, colnames(x), fixed = TRUE))]
    FDR <- p.adjust(x[,5], method = "fdr")
    Bonf <- p.adjust(x[,5], method = "bonferroni")
    x <- do.call("cbind", list(x, FDR = FDR, Bonf = Bonf))
    colnames(x) <- c("PROBE_ID", "SYMBOL", "Log2FC", "Test.Statistic", "P.Value", "FDR", "Bonferroni")
    x <- as.data.frame(x)
    if(input$correction_method1 == "Raw"){
      x <- x[which(x$P.Value <= input$alphalevel1),]
      x <- x[order(x$P.Value, decreasing = FALSE),]
    }
    if(input$correction_method1 == "FDR"){
      x <- x[which(x$FDR <= input$alphalevel1),]
      x <- x[order(x$FDR, decreasing = FALSE),]
    }
    if(input$correction_method1 == "Bonferroni"){
      x <- x[which(x$Bonferroni <= input$alphalevel1),]
      x <- x[order(x$Bonferroni, decreasing = FALSE),]
    }
    if(input$merge==TRUE){
      PROBE_ID<-merge(x,moduleinfo,by="PROBE_ID",all.x = TRUE)
      index4<- which( !(names(PROBE_ID)%in%c("PROBE_ID","SYMBOL","Module","Module_V3","Modulev2_Annotation","Modulev3_Annotation")))
      name2<-names(PROBE_ID)[index4]
      PROBE_ID<-cbind(PROBE_ID[,c("PROBE_ID","SYMBOL","Module","Modulev2_Annotation")],PROBE_ID[,index4])
      names(PROBE_ID)<-c("PROBE_ID","SYMBOL","Module","Modulev2_Annotation",name2)
      x<-PROBE_ID
    }
    dummy <- data.frame("No Genes Present")
    names(dummy) <- "PROBE_ID"
    if(identical(x,dummy) == FALSE){
      if(input$showfc){
        if(input$sign=="Both"){x<-x[which(abs(x$Log2FC)>input$fcval),]}
        if(input$sign=="+"){x<-x[which(x$Log2FC>input$fcval),]}
        if(input$sign=="-"){x<-x[which(x$Log2FC<(-input$fcval)),]}
      }
    }
    if(nrow(x) == 0){
      x <- data.frame("No Genes Present")
      names(x) <- "PROBE_ID"
    }
    return(x)
  })
  
  sel_genelists <- reactive({
    if(length(grep("All", input$comparisons_download)) > 0){
      nam <- names(dataset())
      index <- grep("P.Value",nam,fixed=T)
      p.names <- nam[index]
      comps <- gsub("P.Value for ", "", p.names)
    }
    else{
      comps <- input$comparisons_download
    }
    x.all <- list()
    for(i in 1:length(comps)){
      x <- dataset()
      x <- x[,c(1,2,grep(comps[i], colnames(x), fixed = TRUE))]
      FDR <- p.adjust(x[,5], method = "fdr")
      Bonf <- p.adjust(x[,5], method = "bonferroni")
      x <- do.call("cbind", list(x, FDR = FDR, Bonf = Bonf))
      colnames(x) <- c("PROBE_ID", "SYMBOL", "Log2FC", "Test.Statistic", "P.Value", "FDR", "Bonferroni")
      x <- as.data.frame(x)
      if(input$correction_method1 == "Raw"){
        x <- x[which(x$P.Value <= input$alphalevel1),]
        x <- x[order(x$P.Value, decreasing = FALSE),]
      }
      if(input$correction_method1 == "FDR"){
        x <- x[which(x$FDR <= input$alphalevel1),]
        x <- x[order(x$FDR, decreasing = FALSE),]
      }
      if(input$correction_method1 == "Bonferroni"){
        x <- x[which(x$Bonferroni <= input$alphalevel1),]
        x <- x[order(x$Bonferroni, decreasing = FALSE),]
      }
      if(input$merge==TRUE){
        PROBE_ID<-merge(x,moduleinfo,by="PROBE_ID",all.x = TRUE)
        index4<- which( !(names(PROBE_ID)%in%c("PROBE_ID","SYMBOL","Module","Module_V3","Modulev2_Annotation","Modulev3_Annotation")))
        name2<-names(PROBE_ID)[index4]
        PROBE_ID<-cbind(PROBE_ID[,c("PROBE_ID","SYMBOL","Module","Modulev2_Annotation")],PROBE_ID[,index4])
        names(PROBE_ID)<-c("PROBE_ID","SYMBOL","Module","Modulev2_Annotation",name2)
        x<-PROBE_ID
      }
      dummy <- data.frame("No Genes Present")
      names(dummy) <- "PROBE_ID"
      if(identical(x,dummy) == FALSE){
        if(input$showfc){
          if(input$sign=="Both"){x<-x[which(abs(x$Log2FC)>input$fcval),]}
          if(input$sign=="+"){x<-x[which(x$Log2FC>input$fcval),]}
          if(input$sign=="-"){x<-x[which(x$Log2FC<(-input$fcval)),]}
        }
      }
      if(nrow(x) == 0){
        x <- data.frame("No Genes Present")
        names(x) <- "PROBE_ID"
      }
      x.all[[i]] <- x
    }
    names(x.all) <- comps
    return(x.all)
  })



  sig_ind<-reactive({
    if(is.null(values$illumina) == FALSE){
      if(values$illumina == FALSE){
        if(is.null(values$moduleinfo2)){
          modnames<-unique(moduleinfo2$Module)
          modnum<-gsub("M","",modnames)
          modvec<-as.numeric(unlist(strsplit(modnum,".",fixed=TRUE)))
          modmat<-matrix(modvec,ncol=2,byrow=T)
          modordered<-factor(modnames[order(modmat[,1],modmat[,2])],levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE)
          modordnum<-table(factor(moduleinfo2$Module,levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE))
          mod_for_merge<-data.frame(Module=names(modordnum),Size=as.vector(modordnum))
        }
        else{
          modnames<-unique(values$moduleinfo2$Module)
          modnum<-gsub("M","",modnames)
          modvec<-as.numeric(unlist(strsplit(modnum,".",fixed=TRUE)))
          modmat<-matrix(modvec,ncol=2,byrow=T)
          modordered<-factor(modnames[order(modmat[,1],modmat[,2])],levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE)
          modordnum<-table(factor(values$moduleinfo2$Module,levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE))
          mod_for_merge<-data.frame(Module=names(modordnum),Size=as.vector(modordnum))
        }
      }
      else{
        modnames<-unique(moduleinfo$Module)[-1]
        modnum<-gsub("M","",modnames)
        modvec<-as.numeric(unlist(strsplit(modnum,".",fixed=TRUE)))
        modmat<-matrix(modvec,ncol=2,byrow=T)
        modordered<-factor(modnames[order(modmat[,1],modmat[,2])],levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE)
        modordnum<-table(factor(moduleinfo$Module,levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE))
        mod_for_merge<-data.frame(Module=names(modordnum),Size=as.vector(modordnum))
      }
    }
    if(is.null(values$illumina)){
      modnames<-unique(moduleinfo$Module)[-1]
      modnum<-gsub("M","",modnames)
      modvec<-as.numeric(unlist(strsplit(modnum,".",fixed=TRUE)))
      modmat<-matrix(modvec,ncol=2,byrow=T)
      modordered<-factor(modnames[order(modmat[,1],modmat[,2])],levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE)
      modordnum<-table(factor(moduleinfo$Module,levels=modnames[order(modmat[,1],modmat[,2])],ordered=TRUE))
      mod_for_merge<-data.frame(Module=names(modordnum),Size=as.vector(modordnum))
    }
    #Creating significant indexes for all comparisons and caluclating up down counts
    dataset <- dataset()
    colnames(dataset) <- gsub(" ", ".", colnames(dataset))
    nam<-names(dataset)
    index<-grep("P.Value",nam,fixed=T)
    p.names<-nam[index]
    y<-as.matrix(dataset[,p.names])
    if(input$correction_method1=="Raw"){dat<-y
      dat[dat==1]<-.999
      dat[dat<=input$alphalevel1]<-1
      dat[dat<1]<-0
    }
    if(input$correction_method1=="FDR"){dat<-apply(y,2,p.adjust,method="fdr")
      dat[dat==1]<-.999
      dat[dat<=input$alphalevel1]<-1
      dat[dat<1]<-0
    }
    if(input$correction_method1=="Bonferroni"){dat<-apply(y,2,p.adjust,method="bonferroni")
      dat[dat==1]<-.999
      dat[dat<=input$alphalevel1]<-1
      dat[dat<1]<-0
    }
    fcindex<-grep("Estimate",colnames(dataset),fixed=T)
    dataset.fcindex <- dataset[,fcindex]
    if(input$showfc == TRUE){
      dataset.fcindex[dataset.fcindex < input$fcval & dataset.fcindex > -input$fcval] <- 0
    }
    up_down<-sign(dataset.fcindex)
    dat<-up_down*dat
    if(input$showfc == TRUE){
      if(input$sign == "+"){
        dat[dat == -1] <- 0
      }
      if(input$sign == "-"){
        dat[dat == 1] <- 0
      }
    }
    dat<-data.frame(cbind(as.character(dataset$PROBE_ID),dat))
    if(is.null(values$illumina) == FALSE){
      if(values$illumina == FALSE){
        if(is.null(values$moduleinfo2)){
          names(dat)<-c("SYMBOL",p.names)
          dat_mod <- merge(dat,moduleinfo2, by = "SYMBOL")
          dat_mod<-dat_mod[,c("SYMBOL","Module",p.names)]
        }
        else{
          names(dat)<-c("SYMBOL",p.names)
          dat_mod <- merge(dat,values$moduleinfo2, by = "SYMBOL")
          dat_mod<-dat_mod[,c("SYMBOL","Module",p.names)]
        }
      }
      else{
        names(dat)<-c("PROBE_ID",p.names)
        dat_mod<-merge(dat,moduleinfo,by="PROBE_ID")
        dat_mod<-dat_mod[,c("PROBE_ID","Module",p.names)]
      }
    }

    if(is.null(values$illumina)){
      names(dat)<-c("PROBE_ID",p.names)
      dat_mod<-merge(dat,moduleinfo,by="PROBE_ID")
      dat_mod<-dat_mod[,c("PROBE_ID","Module",p.names)]
    }

    if(length(p.names)==1){
      dat_mod[,-c(1,2)]<-as.numeric(as.character(dat_mod[,-c(1,2)]))
    }
    count_matrix<-aggregate(as.formula(paste(names(dat_mod)[3],"~","Module",sep="")),data=dat_mod,sum)
    if(length(p.names)>1){
      for(i in 2:(length(p.names))){
        count_matrix<-cbind(count_matrix,aggregate(as.formula(paste(names(dat_mod)[i+2],"~","Module",sep="")),data=dat_mod,sum)[,2])
      }
    }
    #code for missing modules if present
    missing_mods<-setdiff(modordered,count_matrix$Module)
    miss_matrix<-cbind(missing_mods,matrix(rep(0,length(missing_mods)*(dim(count_matrix)[2]-1)),length(missing_mods),(dim(count_matrix)[2]-1)))
    count_matrix<-rbind(as.matrix(count_matrix),as.matrix(miss_matrix))
    count_frame<-data.frame(count_matrix)
    names(count_frame)<-c("Module",p.names)
    mod_blanks <- which(count_frame$Module == "")
    if(identical(mod_blanks,integer(0)) == TRUE){
      count_frame <- count_frame
    }
    if(identical(mod_blanks,integer(0)) == FALSE){
      count_frame<-count_frame[-which(count_frame$Module==""),]
    }
    count_frame_merge<-merge(mod_for_merge,count_frame,by="Module",sort=FALSE)
    if(length(p.names)>1){
      prop_matrix<-(1/count_frame_merge$Size)*t(as.matrix(apply(as.matrix(count_frame_merge[,p.names]),1,as.numeric)))
    }
    if(length(p.names)==1){
      prop_matrix<-(1/count_frame_merge$Size)*as.matrix(apply(as.matrix(count_frame_merge[,p.names]),1,as.numeric))
    }
    rownames(prop_matrix)<-count_frame_merge$Module
    colnames(prop_matrix)<-p.names
    prop_matrix2 <- prop_matrix
    prop_matrix[(prop_matrix<.1 & prop_matrix> -.1)]<-0
    prop_matrix
    z <- list(prop_matrix = prop_matrix, prop_matrix2 = prop_matrix2)
    return(z)
  })


  output$genelisttable <-renderDataTable({
    y<-genelist()
    y$SYMBOL <- paste("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",y$SYMBOL," target = '_blank'",'>',y$SYMBOL,"</a>",sep='')
    escape = FALSE
    y
  },escape=FALSE)

  sig_exp_file <- reactive({
    g <- genelist()
    g <- values$final_expression[which(values$final_expression$PROBE_ID %in% g$PROBE_ID),]
    z = list(g = g)
    return(z)
  })

  output$test1<-renderUI({
    if(is.null(values$h3_rowdendro)){
      if(values$hc == TRUE){
        try<-list("Baseline Median Normalized" = 1,
                  "Baseline Healthy Normalized" = 2)
      }
      if(values$hc == FALSE){
        try<-list("Baseline Median Normalized" = 1)
      }
    }
    if(!is.null(values$h3_rowdendro)){
      if(values$hc == TRUE){
        try<-list("Baseline Median Normalized" = 1,
                  "Baseline Healthy Normalized" = 2,
                  "All Samples Median Normalized" = 3,
                  "All Samples Healthy Normalized"=4,
                  "All Samples Baseline Normalized"=5)
      }
      if(values$hc == FALSE){
        try<-list("Baseline Median Normalized" = 1,
                  "All Samples Median Normalized" = 3,
                  "All Samples Baseline Normalized"=5)
      }
    }
    selectInput("set1", "Select heatmap:",as.list(try))
  })

  output$ResponderStatus5 <- renderUI({
    level_length <- lapply(values$design, function(x) length(levels(x)))
    level_length <- unname(unlist(level_length))
    selectInput("responderStatus5", "Choose a Group for Association Test", choices = names(values$design[level_length > 1 & level_length < length(values$design[,1])]), selected = input$responder_var)
  })

  res_stat2 <- reactive({
    input$responderStatus5
  })

  output$group_label_probe1<-renderUI({
    selectInput("Labelprobe1","Select variables to label samples (additional to order variables):",choices=names(values$design),selected=c(input$responder_var,input$patient_id),multiple=T)
  })

  output$TopTierProbe1<-renderUI({
    selectInput("topprobe1","Ordering columns: variable 1",c(names(values$design),"NA"),values$responder_var)
  })

  output$MidTierProbe1<-renderUI({
    selectInput("midprobe1","Ordering columns: variable 2",c(names(values$design),"NA"),values$patient_id)
  })

  output$LowTierProbe1<-renderUI({
    selectInput("bottomprobe1","Ordering columns: variable 3",c(names(values$design),"NA"),"NA")
  })

  order_varsProbe1 <- reactive({
    x<-c(values$responder_var, values$patient_id, "NA")
    y <- x[which(x!="NA")]
    if(!is.null(input$topprobe1) & !is.null(input$midprobe1) & !is.null(input$bottomprobe1)){
      x<-c(input$topprobe1,input$midprobe1,input$bottomprobe1)
      y <- x[which(x!="NA")]
      return(y)
    }
    return(y)
  })

  output$subsetProbeVariable1<-renderUI({selectInput("subsetProbeVar1","Variable 1 to subset heatmap:",c(names(values$design)), values$responder_var)})
  output$subsetProbeValue1<-renderUI({selectInput("subsetProbeVal1","Value(s) of Variable to subset heatmap:",c(unique(as.character(values$design[,input$subsetProbeVar1]))),c(unique(as.character(values$design[,input$subsetProbeVar1])))[1],multiple=TRUE)})
  output$subsetProbeVariable3<-renderUI({selectInput("subsetProbeVar3","Variable 2 to subset heatmap:",c(names(values$design)), values$time_var)})
  output$subsetProbeValue3<-renderUI({selectInput("subsetProbeVal3","Value(s) of Variable to subset heatmap:",c(unique(as.character(values$design[,input$subsetProbeVar3]))),c(unique(as.character(values$design[,input$subsetProbeVar3])))[1],multiple=TRUE)})

  heatmapdata2 <-reactive({
    final_expression = sig_exp_file()$g
    if (input$set1==1){#heatmapbase1
      base_sample_name=values$design$columnname[values$design[,baseline_var]==values$baseline_val]
      ind1<-which(colnames(final_expression)%in%c("PROBE_ID","SYMBOL"))
      ind2<-which(colnames(final_expression)%in%base_sample_name)
      exp_base_sam=final_expression[,c(ind1,ind2)]
      des_base_sam=values$design[which(values$design$columnname%in%colnames(exp_base_sam)),]
      y<-data.manipulate(exp=exp_base_sam,des=des_base_sam,values$baseline_val,longitudinal=FALSE,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=TRUE)
    }
    if (input$set1==2){
      base_sample_name=values$design$columnname[values$design[,baseline_var]==values$baseline_val]
      ind1<-which(colnames(final_expression)%in%c("PROBE_ID","SYMBOL"))
      ind2<-which(colnames(final_expression)%in%base_sample_name)
      exp_base_sam=final_expression[,c(ind1,ind2)]
      des_base_sam=values$design[which(values$design$columnname%in%colnames(exp_base_sam)),]
      y<-data.manipulate(exp=exp_base_sam,des=des_base_sam,values$control_var,values$control_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=FALSE)
    }
    if (input$set1==3){
      y<-data.manipulate(exp=final_expression,des=values$design[which(values$design$columnname %in% colnames(final_expression)),],values$baseline_var,values$baseline_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=TRUE)
    }
    if (input$set1==4){
      y<-data.manipulate(exp=final_expression,des=values$design[which(values$design$columnname %in% colnames(final_expression)),],values$control_var,values$control_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=FALSE)
    }
    if (input$set1==5){#heatmap3
      if(values$hc==TRUE){
        des_w_controls<-values$design[which(values$design$columnname %in% colnames(final_expression)),]
        des_wo_controls<-values$design[-which(values$design[,values$control_var]==values$control_val),]
        h5index<-c(1,2,which(colnames(final_expression) %in% des_wo_controls$columnname))
        y<-data.manipulate(exp=final_expression[,h5index],des=des_wo_controls,values$baseline_var,values$baseline_val,longitudinal=TRUE,subjects=values$patient_id,lg2=FALSE,keepbase=FALSE,format="Probes",allsamples=FALSE)
      }
      if(values$hc==FALSE){
        y<-data.manipulate(exp=values$final_expression,des=values$design[which(values$design$columnname %in% colnames(final_expression)),],values$baseline_var,values$baseline_val,longitudinal=TRUE,subjects=values$patient_id,lg2=FALSE,keepbase=FALSE,format="Probes",allsamples=FALSE)
      }
    }
    z<-list(y=y)
    return(z)
  })

  heatmapname1<-reactive({
    if(input$set1==1) heattxt<-"Baseline Median Normalized"
    if(input$set1==2) heattxt<-"Baseline Healthy Normalized"
    if(input$set1==3) heattxt<-"All Samples Median Normalized"
    if(input$set1==4) heattxt<-"All Samples Healthy Normalized"
    if(input$set1==5) heattxt<-"All Samples Normalized to each Subjects Baseline"
    heattxt
  })

  heatmap_order2 <- eventReactive(input$go2,{
    y<-heatmapdata2()$y
    if(input$subsetProbe1){
      y$heatdes<-y$heatdes[which((y$heatdes[,as.character(input$subsetProbeVar1)] %in% as.character(input$subsetProbeVal1)) & (y$heatdes[,as.character(input$subsetProbeVar3)] %in% as.character(input$subsetProbeVal3))),]
      y$heatexp<-y$heatexp[,match(y$heatdes[,"columnname"], colnames(y$heatexp), nomatch=0)]

    }
    group_order<-order_varsProbe1()
    z<-unique(c(order_varsProbe1(),input$Labelprobe1))
    color_groups<-z[which(z!="NA")]
    inside<-paste("y$heatdes$",group_order,sep="")
    des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
    design_ordered<-y$heatdes[des_order,]
    groups= as.data.frame(design_ordered[,color_groups, drop = F])
    z = list(y = y, groups = groups, color_groups = color_groups, design_ordered = design_ordered)
    return(z)
  })

  expression_matrix2 <- eventReactive(input$go2,{
    final_expression = sig_exp_file()$g
    probeids <- final_expression[,1,drop = FALSE]
    symb <- final_expression[,2,drop = FALSE]
    y = heatmap_order2()$y
    design_ordered = heatmap_order2()$design_ordered
    x = y$heatexp[, match(design_ordered[,"columnname" ], colnames(y$heatexp), nomatch=0)]
    y$heatexp <- y$heatexp[, match(design_ordered[,"columnname" ], colnames(y$heatexp), nomatch=0)]
    colnames(y$heatexp) <- design_ordered[,values$sample_id]
    colnames(x)<-design_ordered[,values$sample_id]
    if(input$setcutoff1!=0){
      cut1<-as.numeric(input$setcutoff1)
      x[x>cut1]<-cut1
      x[x<(-cut1)]<--cut1
    }
    z = list(x = x, probeids = probeids, y = y, symb = symb)
    return(z)
  })

  opt_numClust3 <- eventReactive(input$go2,{
    x = expression_matrix2()$x
    design_ordered = heatmap_order2()$design_ordered
    colddm = NA
    dist_x = dist(t(x))
    hcl = fastcluster::hclust(dist_x)
    if(input$ColClustProbe1==TRUE){
      colddm <- as.dendrogram(hcl)
    }
    d = sapply(2:round(nrow(design_ordered)/2), function(y) clValid::dunn(dist_x, cutree(hcl,y)))
    opt_num = which(d == max(d)) + 1
    z = list(opt_num = opt_num, colddm = colddm, hcl = hcl, d = d)
    return(z)
  })

  output$ClusterCuts5 <- renderUI({
    numericInput('ClustCut5', "Number of clusters", min = 2, value = opt_numClust3()$opt_num, step = 1)
  })

  output$clusternumber3 <- renderUI({
    numericInput("clustnumber3", "Number of clusters:", min = 2, value = 2, step = 1)
  })

  clusterx3 <- eventReactive(input$go2,{
    groups = heatmap_order2()$groups
    design_ordered = heatmap_order2()$design_ordered
    hcl = opt_numClust3()$hcl
    if(input$ClusterChoice5){
      clusters = cutree(hcl, input$ClustCut5)
      design_ordered$Clusters = as.character(clusters)
      groups <- cbind(groups, Cluster = design_ordered$Clusters)
    }
    z = list(groups = groups, design_ordered = design_ordered)
    return(z)
  })

  heatmap_colors2 <- eventReactive(input$go2,{
    groups = clusterx3()$groups
    color_groups = heatmap_order2()$color_groups
    if(length(which(names(groups) %in% values$patient_id))){
      groups[,which(names(groups) %in% values$patient_id)]<-as.numeric(as.factor(groups[,which(names(groups) %in% values$patient_id)]))
    }
    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i]) & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }
    else{
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i])){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }
    for(i in 1:ncol(groups)){
      if(is.numeric(groups[,i]) == F){
        groups[,i] <- as.character((groups[,i]))
      }
    }
    if(values$time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10 & groups[,i] != groups[,values$time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }
    else{
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }
    if(values$time_var %in% colnames(groups)){
      groups[,values$time_var] = as.factor(groups[,values$time_var])
    }
    n_groups = ncol(groups)
    palette_set = setPalettes(n_groups)
    factors = unlist(apply(groups, 2, function(x) as.data.frame(factor(x))), recursive=F)
    factor_levels = lapply(factors, getLevels)
    factor_lengths = lapply(factor_levels, function(x) length(x))
    color_gen = NULL
    add_list = list()
    for (i in 1:n_groups) {
      j = palette_set[[i]](as.numeric(factor_lengths[i]))
      color_gen = c(color_gen, j)
      add_list[[i]] = list(rect = list(col = "transparent", fill = j[factors[[i]]]))
    }
    annotation<-c()
    for(i in 1:n_groups){
      annotation<-cbind(annotation,add_list[[i]]$rect$fill)
    }
    annotation<-annotation[,n_groups:1]
    names(annotation)<-color_groups[n_groups:1]
    anno1<-colorRampPalette(c("navy", "yellow", "firebrick3"))(length(unique(groups[,1])))
    names(anno1)<-unique(groups[,1])
    eval(parse(text=paste(names(groups)[1],"=","anno1",sep="")))
    eval(parse(text=paste("first_color=list(",names(groups)[1],"=",names(groups)[1],")")))
    z <- list(groups = groups, first_color = first_color)
    return(z)
  })

  gen_clustTab2 <- reactive({
    if(is.null(opt_numClust3())){return(NULL)}
    design_ordered = heatmap_order2()$design_ordered
    opt_num = opt_numClust3()$opt_num
    hcl = opt_numClust3()$hcl
    clusters = cutree(hcl, input$clustnumber3)
    design_ordered$Clusters = clusters
    tab <- aggregate(design_ordered[[res_stat2()]] ~ design_ordered$Clusters, data = design_ordered, FUN = table)
    tab <- as.data.frame(tab[,2])
    tab1 <- cbind("Cluster" = rownames(tab), tab)
    tab2 <- tab
    low_exp_count <- sapply(1:length(tab), function(y) if(sum(tab[,y]) == 0) y)
    low_exp_count <- unlist(low_exp_count)
    if(is.null(low_exp_count) == F){
      tab2 <- as.data.frame(tab[, -low_exp_count, drop = F])
    }
    tab4 <- cbind("Cluster" = rownames(tab2), tab2)
    if(ncol(tab2) > 1){
      chi_sqr <- chisq.test(tab2)
      fish_exact <- fisher.test(tab2, simulate.p.value = T, B = 10000)
      tab3 = list()
      for(i in 1:nrow(tab4)){
        tab3[[i]] = paste(tab4[i,-1], " ", "(", round(100*(tab4[i,-1]/sum(tab4[i,-1])), 2), "%",")", sep = "")
      }
      tab3 = do.call("rbind", tab3)
      tab3 <- cbind(tab4[,1,drop = F], tab3)
      names(tab3) <- names(tab4)
    }
    else{
      chi_sqr <- list()
      fish_exact <- list()
      tab3 = tab4
    }
    tab1[,1] <- as.character(tab1[,1])
    tab1_total <- rbind(data.frame(tab1[,1,drop = F], stringsAsFactors = FALSE), "Total")
    tab1_colsums <- rbind(tab1[,-1], colSums(tab1[,-1]))
    tab1 <- cbind(tab1_total, tab1_colsums)
    tab1$Total <- rowSums(tab1[,-1])
    y <- list(tab1 = tab1, tab2 = tab2, tab3 = tab3, chi_sqr = chi_sqr, fish_exact = fish_exact)
    return(y)
  })

  plotw <- eventReactive(input$go2,{input$DGE_HeatMapSize1})
  plotH <- eventReactive(input$go2, {input$DGE_HeatMapSize2})
  plotsize_W <- function(){plotw()}
  plotsize_H <- function(){plotH()}
  Treeheight <- eventReactive(input$go2, {input$DGE_TreeHeight})
  fontSize <- eventReactive(input$go2, {input$DGE_FontSize})
  leg_size <- eventReactive(input$go2,{input$DGE_LegendSize})
  resolution <- function(){input$DGE_PlotRes}

  row_cluster <- eventReactive(input$go2,{
    if(input$hm1rclu == TRUE){
      rowclust <- TRUE
    }
    if(input$hm1rclu == FALSE){
      rowclust <- NA
    }
    return(rowclust)
  })

  output$heatmap1 <- renderPlot({
    if(all(expression_matrix2()$x == 0)){
      withProgress(message = '',
                   detail = 'Generating the Options...', value = 1,{
                     aheatmap2(expression_matrix2()$x,Rowv=row_cluster(),Colv = opt_numClust3()$colddm,cexRow=1.2, treeheight = Treeheight(), fontsize = fontSize(), annheight = leg_size(), color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),annCol = heatmap_colors2()$groups,annColors= heatmap_colors2()$first_color,labRow=NA,breaks=0,legend = FALSE)
                   })
    }
    else{
      withProgress(message = '',
                   detail = 'Generating the Options...', value = 1,{
                     aheatmap2(expression_matrix2()$x,Rowv=row_cluster(),Colv = opt_numClust3()$colddm, cexRow=1.2, treeheight = Treeheight(), fontsize = fontSize(), annheight = leg_size(), color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),annCol = heatmap_colors2()$groups,annColors= heatmap_colors2()$first_color,labRow=NA,breaks=0)
                   })
    }
  }, height = plotsize_H, width = plotsize_W)

  output$downloadHeatmap2 <- downloadHandler(
    filename = function() {paste('SignificantProbes_Heatmap','.png', sep = '')},
    content = function(file){
      png(file, width = plotsize_W(), height = plotsize_H(), res = resolution())
      print(aheatmap2(expression_matrix2()$x,Rowv=TRUE,Colv=opt_numClust3()$colddm,cexRow=1.2, treeheight = input$DGE_TreeHeight, fontsize = input$DGE_FontSize, annheight = input$DGE_LegendSize, color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),annCol = heatmap_colors2()$groups,annColors= heatmap_colors2()$first_color,labRow=NA,breaks=0))
      dev.off()
    }
  )

  heatmap_download2 <- reactive({
    heatmapdata <- cbind(expression_matrix2()$symb, expression_matrix2()$y$heatexp)
    heatmapdata <- cbind(expression_matrix2()$probeids, heatmapdata)
    return(heatmapdata)
  })

  output$downloadHeatmap1 <- downloadHandler(

    filename = function() {paste(values$project_name,'_',heatmapname1(),'.csv', sep='')  },
    content = function(file) {
      write.csv(heatmap_download2(), file, row.names = FALSE)
    }
  )

  output$clusterplot4 <- renderPlot({
    barplot(opt_numClust3()$d, names.arg = 2:round(nrow(heatmap_order2()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index")
  })

  output$downloadClusterPlot4 <- downloadHandler(
    filename = function() {paste('Cluster_Plot_SignificantProbes','.png', sep = '')},
    content = function(file){
      png(file, width = 800)
      print(barplot(opt_numClust3()$d, names.arg = 2:round(nrow(heatmap_order2()$design_ordered)/2), xlab = "Number of Clusters", ylab = "Dunn's Index"))
      dev.off()
    }
  )

  output$cluster_output5 <- renderTable({
    if(is.null(gen_clustTab2())){return(NULL)}
    gen_clustTab2()$tab1
  }, include.rownames = FALSE, digits = 0)

  output$explanation4 <- renderText({
    if(is.null(gen_clustTab2())){return(NULL)}
    if(ncol(gen_clustTab2()$tab2) > 1){
      print("The table below is the table the tests are run on. It is the same as the table above, except the columns with zero counts have been deleted.")
    }
  })

  output$cluster_tab4 <- renderTable({
    if(is.null(gen_clustTab2())){return(NULL)}
    if(ncol(gen_clustTab2()$tab2) > 1){
      gen_clustTab2()$tab3
    }
  }, include.rownames = FALSE)

  output$chisquare_test5 <- renderText({
    if(is.null(gen_clustTab2())){return(NULL)}
    if(ncol(gen_clustTab2()$tab2) > 1){
      if(gen_clustTab2()[[4]]$p.value < .001){
        paste("Chi_square statistic = ", round(gen_clustTab2()[[4]]$statistic, 2), ",", "p-value < .001")
      }
      else{
        paste("Chi-Square Test Statistic = ", round(gen_clustTab2()[[4]]$statistic, 2), ",", "p-value =", round(gen_clustTab2()[[4]]$p.value, 3))
      }
    }
  })

  output$fisher5 <- renderText({
    if(is.null(gen_clustTab2())){return(NULL)}
    if(ncol(gen_clustTab2()$tab2) > 1){
      if(gen_clustTab2()[[5]]$p.value < .001){
        paste("Fishers Exact Test: p-value < .001")
      }
      else{
        paste("Fishers Exact Test: p-value = ", round(gen_clustTab2()[[5]]$p.value, 3))
      }
    }

  })

  output$Comparison1 <- renderUI({
    nam <- names(dataset())
    index <- grep("P.Value",nam,fixed=T)
    p.names <- nam[index]
    p.names <- gsub("P.Value for ", "", p.names)
    selectInput("comparison", "Comparison:", p.names, p.names[1])
  })

  output$Comparison2 <- renderUI({
    nam <- names(dataset())
    index <- grep("P.Value",nam,fixed=T)
    p.names <- nam[index]
    p.names <- gsub("P.Value for ", "", p.names)
    selectInput("comparison_LMM", "Comparison:", p.names, p.names[1], multiple = TRUE)
  })
  
  output$selectComps <- renderUI({
    nam <- names(dataset())
    index <- grep("P.Value",nam,fixed=T)
    p.names <- nam[index]
    p.names <- gsub("P.Value for ", "", p.names)
    p.names <- c("All", p.names)
    selectInput("comparisons_download", "Comparison:", p.names, p.names[1], multiple = TRUE)
  })

  output$modselTF <- renderUI({
    if(input$rowselect == TRUE){
      return(selectInput("FirstSixLMM","Modules to include:", c("Customized"), "Customized"))
    }
    selectInput("FirstSixLMM","Modules to include:",c("All","First Seven Rounds","Only Annotated"), "First Seven Rounds")
  })

  output$modselTF2 <- renderUI({
    if(input$rowselect2 == TRUE){
      return(selectInput("FirstSix3","Modules to include:", c("Customized"), "Customized"))
    }
    selectInput("FirstSix3","Modules to include:",c("All","First Seven Rounds","Only Annotated"),"First Seven Rounds")
  })

  output$modselTF3 <- renderUI({
    if(input$rowselect3 == TRUE){
      return(selectInput("FirstSix","Modules to include:", c("Customized"), "Customized"))
    }
    selectInput("FirstSix","Modules to include:",c("All","First Seven Rounds","Only Annotated"), "First Seven Rounds")
  })

  gensymb.search <- reactive({
    dat <- dataset()
    genlist <- unique(dat$SYMBOL)
    genlist
  })

  probe.search <- reactive({
    dat <- dataset()
    selectedprobes <- dat$PROBE_ID[which(dat$SYMBOL %in% input$specgene)]
    selectedprobes
  })

  output$specgene1 <- renderUI({
    genlist <- gensymb.search()
    withProgress(message = '',
                 detail = 'Generating the Options...', value = 1,{
                   selectizeInput("specgene", "Gene symbol selection:",
                                  choices = genlist,
                                  selected = genlist[1],options = list(maxOptions = 30))
                 })
  })

  output$probeid1 <- renderUI({
    if(is.null(input$specgene)){return(NULL)}
    else{
      selectedprobes <- probe.search()
      selectizeInput("probeid", "Probe ID selection:", multiple = TRUE, selectedprobes,selected=selectedprobes[1], options = list(maxOptions = 30))
    }
  })

  output$genTable <- renderUI({
    actionButton("go3", "Generate table")
  })

  sgl_flat <- reactive({
    results_file <- dataset()
    sgl <- results_file[grep("P.Value", names(results_file))]
    sgl_1 <- apply(sgl, 2, p.adjust, method = "fdr")
    colnames(sgl_1) <- gsub("P.Value","FDR.P.Value",colnames(sgl_1))
    sgl_2 <- apply(sgl, 2, p.adjust, method = "bonferroni")
    colnames(sgl_2) <- gsub("P.Value","Bonf.P.Value",colnames(sgl_2))
    sgl_flat <- data.frame(results_file, sgl_1, sgl_2)
    sgl_flat
  })

  specgenelist <- eventReactive(input$go3,{
    m_newsgl <- function(sss){
      newsgl_SYMBOL <- sgl_flat()[which(sgl_flat()$SYMBOL %in% sss),][,-2]
      comparisons <- gsub("Estimate.of.", "", colnames(newsgl_SYMBOL)[grep("Estimate", colnames(newsgl_SYMBOL))])
      dat <- reshape2::melt(newsgl_SYMBOL, id.vars = "PROBE_ID")
      dat$variable <- gsub(".for.",".of.",dat$variable, fixed = TRUE)
      splt <- as.data.frame(stringr::str_split_fixed(dat$variable, ".of.",2))
      colnames(splt) <- c("variable", "Comparison")
      dat <- do.call("cbind", list(PROBE_ID = as.character(dat[,1]), splt, value = dat$value))
      newdat <- reshape2::dcast(dat, PROBE_ID + Comparison ~ variable)
      newdat <- newdat[,c(1,2,4,7,6,5,3)]
      return(newdat)
    }
    xx_0 <- m_newsgl(input$specgene)
    xx <- xx_0[which(xx_0$PROBE_ID %in% input$probeid),]
    xx
  })

  output$specgenetable <- renderDataTable({
    specgenelist()
  })

  output$genelistgraph<-renderPlot({
    plot(plotdata()$alpha,plotdata()$raw,type="l",main="Multiple Testing Comparison Plot",col="black",xlim=c(0,1),ylim=c(0,1),xlab=expression(alpha),ylab="% of Total Probes in Gene List")
    lines(plotdata()$alpha,plotdata()$fdr,col="red",lwd=2)
    lines(plotdata()$alpha,plotdata()$bonf,col="green")
    legend("bottomright",legend=c("Raw P-value","FDR","Bonf."), lty=c(1,1,1),col=c("black","red","green") )
    axis(4,at=1:10/10,labels=round(1:10/10*dim(dataset())[1],0))
    lines(c(0,1),c(0,1),lty=2)
  })

  output$downloadCG <- downloadHandler(
    filename = function(){paste0("Gene_Search_List","_", input$specgene,"_with_selected_probes", ".csv")},
    content = function(file){
      write.csv(specgenelist(), file,row.names = FALSE)
    }
  )

  output$downloadSC <- downloadHandler(
    filename = function() {paste("Significance_Comparison_Overview", '_', input$alphalevel2, '.csv', sep = '') },
    content = function(file){
      write.csv(siglist()$z, file,row.names = FALSE)
    }
  )

  output$downloadData <- downloadHandler(
    filename = function() {paste(substring(pcomp(),12),'_',input$correction_method1,input$alphalevel1,'.csv', sep='')  },
    content = function(file) {
      write.csv(genelist(), file,row.names = FALSE)
    }
  )
  
  output$downloadSelComp <- downloadHandler(
    filename = function() {paste(input$ziptext,".zip", sep = "")},
    content = function(file) {
      file.names <- c()
      for(i in 1:length(sel_genelists())){
        file.names[i] <- paste(names(sel_genelists())[i], "_", input$correction_method1, input$alphalevel1, ".csv", sep = "")
        write.csv(sel_genelists()[[i]], file = paste(names(sel_genelists())[i], "_", input$correction_method1, input$alphalevel1, ".csv", sep = ""), row.names = FALSE)
      }
      zip(zipfile = file, files = file.names)
    },
    contentType = "application/zip"
  )

  output$readme <-renderText({paste("The following reports are designed to investigate the results of Differential Gene Expression Analysis.
                                    Below is a summary by a biostatistics team member of the statistical model that was conducted for this project.
                                    The following tabs provide the user with the ability to create their own gene lists for any comparison the analysis conducted.
                                    In the DGE: Gene Lists tab, the user can create the list using raw p-value thresholds or multiple testing corrections by tweaking the inputs on the left hand side of interface.
                                    In addition to the list, a module map of the first Seven rounds of the Baylor Modules is produced.  The percentages are calculated by taking the significant probes in the current gene list
                                    and calculating the percentages of significant probes that are significantly up or down regulated.") })
  output$intro <-renderText({paste("Number of probes in Gene list by method and alpha parameter.")})

  LMMmap1 <- function(){input$LMMmap_1}
  LMMmap2 <- function(){input$LMMmap_2}
  resolution1 <- function(){input$DGE_PlotRes1}

  lmm_module <- reactive({
    rowv = NA
    if(input$comparselect == FALSE & input$rowselect == FALSE){
      data<-data.frame(cbind(rownames(sig_ind()$prop_matrix),sig_ind()$prop_matrix))
      names(data)=c("Module",colnames(sig_ind()$prop_matrix))
      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      if(input$FirstSixLMM=="First Seven Rounds"){
        round<-97#62
        data = data[1:round, ]
      }
      if(input$FirstSixLMM=="Only Annotated"){
        anno_index<-which(data$Module %in% module_annotations$Module)
        data = data[anno_index, ]
      }
      x = as.matrix(t(data[,-1]))
      if(dim(x)[1]==1){
        x<-100*t(sapply(data.frame(x), function(x) as.numeric(as.character(x))) )
        num_col<-ncol(x)
        x=matrix(x[,1:num_col],1,num_col)
      }
      if(dim(x)[1]>1){
        x=t(100*(apply(x,1,as.numeric)))
        num_col<-ncol(x)
        x=x[,1:num_col]
      }
      colnames(x) = data[, 1][1:num_col]
      rownames(x) = names(data)[-1]
      rownames(x) = gsub("P.Value.for.","",rownames(x))
      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(colnames(x)==module_annotations$Module[i])
        colnames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],colnames(x)[anno_index2],sep=" ")
      }
      colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      Col.means <- colMeans(x, na.rm = TRUE)
      colddm <- reorder(colddm, Col.means)
      colddm_ord = rev(order.dendrogram(colddm))
      myval <- max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))] <- myval
        x[which(x==min(x))] <- -myval
      }
      if(input$LMMmodorder==FALSE){
        if(dim(x)[1]>1){
          dat = t(x[,colddm_ord])
        }
        if(dim(x)[1]==1){
          newx<-t(as.matrix(x[,colddm_ord]))
          rownames(newx) = names(data)[-1]
          rownames(newx) = gsub("P.Value.for.","",rownames(newx))
          dat = t(newx)
        }
      }
      if(input$LMMmodorder==TRUE){
        rowv = NA
        dat = t(x)
      }
    }

    if(input$comparselect == TRUE & input$rowselect == FALSE){
      namm <- names(dataset())
      index11 <- grep("P.Value",namm,fixed=T)
      nammm <- namm[index11]
      nammm <- gsub("P.Value for ", "", nammm)
      index33 <- match(input$comparison_LMM, nammm)
      data_0 <- data.frame(cbind(rownames(sig_ind()$prop_matrix),sig_ind()$prop_matrix))
      names(data_0)=c("Module",colnames(sig_ind()$prop_matrix))
      data <- data_0[, c(1, index33+1)]
      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      if(input$FirstSixLMM=="First Seven Rounds"){
        round<-97
        data = data[1:round, ]
      }
      if(input$FirstSixLMM=="Only Annotated"){
        anno_index<-which(data$Module %in% module_annotations$Module)
        data = data[anno_index, ]
      }
      x = as.matrix(t(data[,-1]))
      if(dim(x)[1]==1){
        x<-100*t(sapply(data.frame(x), function(x) as.numeric(as.character(x))) )
        num_col<-ncol(x)
        x=matrix(x[,1:num_col],1,num_col)
      }
      if(dim(x)[1]>1){
        x=t(100*(apply(x,1,as.numeric)))
        num_col<-ncol(x)
        x=x[,1:num_col]
      }
      colnames(x) = data[, 1][1:num_col]
      rownames(x) = names(data)[-1]
      rownames(x) = gsub("P.Value.for.", "", rownames(x), fixed = TRUE)
      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(colnames(x)==module_annotations$Module[i])
        colnames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],colnames(x)[anno_index2],sep=" ")
      }
      colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      Col.means <- colMeans(x, na.rm = TRUE)
      colddm <- reorder(colddm, Col.means)
      colddm_ord = rev(order.dendrogram(colddm))
      myval <- max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))] <- myval
        x[which(x==min(x))] <- -myval
      }
      if(input$LMMmodorder==FALSE){
        if(dim(x)[1]>1){
          dat = t(x[,colddm_ord])
        }
        if(dim(x)[1]==1){
          newx<-t(as.matrix(x[,colddm_ord]))
          rownames(newx) = names(data)[-1]
          rownames(newx) = gsub("P.Value.for.","",rownames(newx))
          dat = t(newx)
        }

      }
      if(input$LMMmodorder==TRUE){
        rowv = NA
        dat = t(x)
      }
    }

    if(input$comparselect == TRUE & input$rowselect == TRUE){
      namm <- names(dataset())
      index11 <- grep("P.Value",namm,fixed=T)
      nammm <- namm[index11]
      nammm <- gsub("P.Value for ", "", nammm)
      v_mdnam <- read.csv(input$modsel$datapath, header = TRUE)
      modnames <- as.character(v_mdnam[, 1])
      index33 <- match(input$comparison_LMM, nammm)
      data_0 <- data.frame(cbind(rownames(sig_ind()$prop_matrix),sig_ind()$prop_matrix))
      names(data_0)=c("Module",colnames(sig_ind()$prop_matrix))
      index44 <- match(modnames, data_0$Module)
      data <- data_0[index44, c(1, index33+1)]
      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      x = as.matrix(t(data[,-1]))
      if(dim(x)[1]==1){
        x<-100*t(sapply(data.frame(x), function(x) as.numeric(as.character(x))) )
        num_col<-ncol(x)
        x=matrix(x[,1:num_col],1,num_col)
      }
      if(dim(x)[1]>1){
        x=t(100*(apply(x,1,as.numeric)))
        num_col<-ncol(x)
        x=x[,1:num_col]
      }
      colnames(x) = data[, 1][1:num_col]
      rownames(x) = names(data)[-1]
      rownames(x) = gsub("P.Value.for.", "", rownames(x), fixed = TRUE)
      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(colnames(x)==module_annotations$Module[i])
        colnames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],colnames(x)[anno_index2],sep=" ")
      }
      colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      Col.means <- colMeans(x, na.rm = TRUE)
      colddm <- reorder(colddm, Col.means)
      colddm_ord = rev(order.dendrogram(colddm))
      myval <- max(c(abs(min(x)),max(x)))
      if(min(x) <0 & max(x) > 0){
        x[which(x==max(x))] <- myval
        x[which(x==min(x))] <- -myval
      }
      if(input$LMMmodorder==FALSE){
        if(dim(x)[1]>1){
          dat = t(x[,colddm_ord])
        }
        if(dim(x)[1]==1){
          newx<-t(as.matrix(x[,colddm_ord]))
          rownames(newx) = names(data)[-1]
          rownames(newx) = gsub("P.Value.for.","",rownames(newx))
          dat = t(newx)
        }
      }
      if(input$LMMmodorder==TRUE){
        rowv = NA
        dat = t(x)
      }
    }

    if(input$comparselect == FALSE & input$rowselect == TRUE){
      v_mdnam <- read.csv(input$modsel$datapath, header = TRUE)
      modnames <- as.character(v_mdnam[, 1])
      data_0 <- data.frame(cbind(rownames(sig_ind()$prop_matrix),sig_ind()$prop_matrix))
      names(data_0)=c("Module",colnames(sig_ind()$prop_matrix))
      index44 <- match(modnames, data_0$Module)
      data <- data_0[index44, ]
      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      x = as.matrix(t(data[,-1]))
      if(dim(x)[1]==1){
        x<-100*t(sapply(data.frame(x), function(x) as.numeric(as.character(x))) )
        num_col<-ncol(x)
        x=matrix(x[,1:num_col],1,num_col)
      }
      if(dim(x)[1]>1){
        x=t(100*(apply(x,1,as.numeric)))
        num_col<-ncol(x)
        x=x[,1:num_col]
      }
      colnames(x) = data[, 1][1:num_col]
      rownames(x) = names(data)[-1]
      rownames(x) = gsub("P.Value.for.", "", rownames(x), fixed = TRUE)
      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(colnames(x)==module_annotations$Module[i])
        colnames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],colnames(x)[anno_index2],sep=" ")
      }
      colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      Col.means <- colMeans(x, na.rm = TRUE)
      colddm <- reorder(colddm, Col.means)
      colddm_ord = rev(order.dendrogram(colddm))
      myval <- max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))] <- myval
        x[which(x==min(x))] <- -myval
      }
      if(input$LMMmodorder==FALSE){
        if(dim(x)[1]>1){
          dat = t(x[,colddm_ord])
        }
        if(dim(x)[1]==1){
          newx<-t(as.matrix(x[,colddm_ord]))
          rownames(newx) = names(data)[-1]
          rownames(newx) = gsub("P.Value.for.","",rownames(newx))
          dat = t(newx)
        }
      }
      if(input$LMMmodorder==TRUE){
        rowv = NA
        dat = t(x)
      }
    }
    if(input$LMMdeleterows == TRUE){
      dat <- dat[-which(rowSums(dat) == 0),,drop = FALSE]
    }
    z = list(dat = dat, rowv = rowv, color_palette = color_palette)
    return(z)
  })

  lmm_module_download <- reactive({
    rowv = NA
    if(input$comparselect == FALSE & input$rowselect == FALSE){
      data<-data.frame(cbind(rownames(sig_ind()$prop_matrix2),sig_ind()$prop_matrix2))
      names(data)=c("Module",colnames(sig_ind()$prop_matrix2))
      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      if(input$FirstSixLMM=="First Seven Rounds"){
        round<-97#62
        data = data[1:round, ]
      }
      if(input$FirstSixLMM=="Only Annotated"){
        anno_index<-which(data$Module %in% module_annotations$Module)
        data = data[anno_index, ]
      }
      x = as.matrix(t(data[,-1]))
      if(dim(x)[1]==1){
        x<-100*t(sapply(data.frame(x), function(x) as.numeric(as.character(x))) )
        num_col<-ncol(x)
        x=matrix(x[,1:num_col],1,num_col)
      }
      if(dim(x)[1]>1){
        x=t(100*(apply(x,1,as.numeric)))
        num_col<-ncol(x)
        x=x[,1:num_col]
      }
      colnames(x) = data[, 1][1:num_col]
      rownames(x) = names(data)[-1]
      rownames(x) = gsub("P.Value.for.", "", rownames(x), fixed = TRUE)
      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(colnames(x)==module_annotations$Module[i])
        colnames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],colnames(x)[anno_index2],sep=" ")
      }
      colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      Col.means <- colMeans(x, na.rm = TRUE)
      colddm <- reorder(colddm, Col.means)
      colddm_ord = rev(order.dendrogram(colddm))
      myval <- max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))] <- myval
        x[which(x==min(x))] <- -myval
      }
      if(input$LMMmodorder==FALSE){
        if(dim(x)[1]>1){
          dat = t(x[,colddm_ord])
        }
        if(dim(x)[1]==1){
          newx<-t(as.matrix(x[,colddm_ord]))
          rownames(newx) = names(data)[-1]
          rownames(newx) = gsub("P.Value.for.","",rownames(newx))
          dat = t(newx)
        }
      }
      if(input$LMMmodorder==TRUE){
        rowv = NA
        dat = t(x)
      }
    }

    if(input$comparselect == TRUE & input$rowselect == FALSE){
      namm <- names(dataset())
      index11 <- grep("P.Value",namm,fixed=T)
      nammm <- namm[index11]
      nammm <- gsub("P.Value for ", "", nammm)
      index33 <- match(input$comparison_LMM, nammm)
      data_0 <- data.frame(cbind(rownames(sig_ind()$prop_matrix2),sig_ind()$prop_matrix2))
      names(data_0)=c("Module",colnames(sig_ind()$prop_matrix2))
      data <- data_0[, c(1, index33+1)]
      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      if(input$FirstSixLMM=="First Seven Rounds"){
        round<-97
        data = data[1:round, ]
      }
      if(input$FirstSixLMM=="Only Annotated"){
        anno_index<-which(data$Module %in% module_annotations$Module)
        data = data[anno_index, ]
      }
      x = as.matrix(t(data[,-1]))
      if(dim(x)[1]==1){
        x<-100*t(sapply(data.frame(x), function(x) as.numeric(as.character(x))) )
        num_col<-ncol(x)
        x=matrix(x[,1:num_col],1,num_col)
      }
      if(dim(x)[1]>1){
        x=t(100*(apply(x,1,as.numeric)))
        num_col<-ncol(x)
        x=x[,1:num_col]
      }
      colnames(x) = data[, 1][1:num_col]
      rownames(x) = names(data)[-1]
      rownames(x) = gsub("P.Value.for.", "", rownames(x), fixed = TRUE)
      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(colnames(x)==module_annotations$Module[i])
        colnames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],colnames(x)[anno_index2],sep=" ")
      }
      colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      Col.means <- colMeans(x, na.rm = TRUE)
      colddm <- reorder(colddm, Col.means)
      colddm_ord = rev(order.dendrogram(colddm))
      myval <- max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))] <- myval
        x[which(x==min(x))] <- -myval
      }
      if(input$LMMmodorder==FALSE){
        if(dim(x)[1]>1){
          dat = t(x[,colddm_ord])
        }
        if(dim(x)[1]==1){
          newx<-t(as.matrix(x[,colddm_ord]))
          rownames(newx) = names(data)[-1]
          rownames(newx) = gsub("P.Value.for.","",rownames(newx))
          dat = t(newx)
        }
      }
      if(input$LMMmodorder==TRUE){
        rowv = NA
        dat = t(x)
      }
    }

    if(input$comparselect == TRUE & input$rowselect == TRUE){
      namm <- names(dataset())
      index11 <- grep("P.Value",namm,fixed=T)
      nammm <- namm[index11]
      nammm <- gsub("P.Value for ", "", nammm)
      v_mdnam <- read.csv(input$modsel$datapath, header = TRUE)
      modnames <- as.character(v_mdnam[, 1])
      index33 <- match(input$comparison_LMM, nammm)
      data_0 <- data.frame(cbind(rownames(sig_ind()$prop_matrix2),sig_ind()$prop_matrix2))
      names(data_0)=c("Module",colnames(sig_ind()$prop_matrix2))
      index44 <- match(modnames, data_0$Module)
      data <- data_0[index44, c(1, index33+1)]
      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      x = as.matrix(t(data[,-1]))
      if(dim(x)[1]==1){
        x<-100*t(sapply(data.frame(x), function(x) as.numeric(as.character(x))) )
        num_col<-ncol(x)
        x=matrix(x[,1:num_col],1,num_col)
      }
      if(dim(x)[1]>1){
        x=t(100*(apply(x,1,as.numeric)))
        num_col<-ncol(x)
        x=x[,1:num_col]
      }
      colnames(x) = data[, 1][1:num_col]
      rownames(x) = names(data)[-1]
      rownames(x) = gsub("P.Value.for.", "", rownames(x), fixed = TRUE)
      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(colnames(x)==module_annotations$Module[i])
        colnames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],colnames(x)[anno_index2],sep=" ")
      }
      colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      Col.means <- colMeans(x, na.rm = TRUE)
      colddm <- reorder(colddm, Col.means)
      colddm_ord = rev(order.dendrogram(colddm))
      myval <- max(c(abs(min(x)),max(x)))
      if(min(x) <0 & max(x) > 0){
        x[which(x==max(x))] <- myval
        x[which(x==min(x))] <- -myval
      }
      if(input$LMMmodorder==FALSE){
        if(dim(x)[1]>1){
          dat = t(x[,colddm_ord])
        }
        if(dim(x)[1]==1){
          newx<-t(as.matrix(x[,colddm_ord]))
          rownames(newx) = names(data)[-1]
          rownames(newx) = gsub("P.Value.for.","",rownames(newx))
          dat = t(newx)
        }
      }
      if(input$LMMmodorder==TRUE){
        rowv = NA
        dat = t(x)
      }
    }

    if(input$comparselect == FALSE & input$rowselect == TRUE){
      v_mdnam <- read.csv(input$modsel$datapath, header = TRUE)
      modnames <- as.character(v_mdnam[, 1])
      data_0 <- data.frame(cbind(rownames(sig_ind()$prop_matrix2),sig_ind()$prop_matrix2))
      names(data_0)=c("Module",colnames(sig_ind()$prop_matrix2))
      index44 <- match(modnames, data_0$Module)
      data <- data_0[index44, ]
      color_palette = colorRampPalette(c("blue", "white", "red"))(1000)
      color_palette[c(495:505)] = "#FFFFFF"
      x = as.matrix(t(data[,-1]))
      if(dim(x)[1]==1){
        x<-100*t(sapply(data.frame(x), function(x) as.numeric(as.character(x))) )
        num_col<-ncol(x)
        x=matrix(x[,1:num_col],1,num_col)
      }
      if(dim(x)[1]>1){
        x=t(100*(apply(x,1,as.numeric)))
        num_col<-ncol(x)
        x=x[,1:num_col]
      }
      colnames(x) = data[, 1][1:num_col]
      rownames(x) = names(data)[-1]
      rownames(x) = gsub("P.Value.for.", "", rownames(x), fixed = TRUE)
      for(i in 1:(dim(module_annotations)[1])){
        anno_index2<-which(colnames(x)==module_annotations$Module[i])
        colnames(x)[anno_index2]<-paste(module_annotations$Modulev2_Annotation[i],colnames(x)[anno_index2],sep=" ")
      }
      colddm = as.dendrogram(fastcluster::hclust(dist(t(x))))
      Col.means <- colMeans(x, na.rm = TRUE)
      colddm <- reorder(colddm, Col.means)
      colddm_ord = rev(order.dendrogram(colddm))
      myval <- max(c(abs(min(x)),max(x)))
      if(min(x) < 0 & max(x) > 0){
        x[which(x==max(x))] <- myval
        x[which(x==min(x))] <- -myval
      }
      if(input$LMMmodorder==FALSE){
        if(dim(x)[1]>1){
          dat = t(x[,colddm_ord])
        }
        if(dim(x)[1]==1){
          #rowv = FALSE
          newx<-t(as.matrix(x[,colddm_ord]))
          rownames(newx) = names(data)[-1]
          rownames(newx) = gsub("P.Value.for.","",rownames(newx))
          dat = t(newx)
        }
      }
      if(input$LMMmodorder==TRUE){
        rowv = NA
        dat = t(x)
      }
    }
    if(input$LMMdeleterows == TRUE){
      dat <- dat[-which(rowSums(dat) == 0),,drop = FALSE]
    }
    z = list(dat = dat, rowv = rowv, color_palette = color_palette)
    return(z)
  })

  output$downloadLMMModMap <- downloadHandler(

    filename = function() {paste(values$project_name,"_LMM_Module_Overview",'_',input$correction_method1,input$alphalevel1,'.csv', sep='')},
    content = function(file) {
      write.csv((lmm_module_download()$dat/100)+1, file)
    }
  )

  output$vennComparison <- renderUI({
    nam <- names(dataset())
    index <- grep("P.Value",nam,fixed=T)
    p.names <- nam[index]
    p.names <- gsub("P.Value for ", "", p.names)
    selectInput("Vcomparison", "Comparison:", p.names, p.names[1], multiple = TRUE)
  })

  output$include <- renderUI({
    n = length(input$Vcomparison)
    Vcomparison = input$Vcomparison
    selectInput("Include", "Include:", Vcomparison, Vcomparison[1:n], multiple = TRUE)
  })

  output$exclude <- renderUI({
    comparisons <- setdiff(input$Vcomparison, input$Include)
    selectInput("Exclude", "Exclude:", comparisons, comparisons[1], multiple = TRUE)
  })

  venndata <- reactive({
    results_file <- dataset()
    TF_Raw = siglist()$TF_Raw1
    names(TF_Raw) <- gsub("P.Value.for.","",names(TF_Raw))
    TF_FDR = siglist()$TF_FDR1
    names(TF_FDR) <- gsub("FDR.P.Value.for.","",names(TF_FDR))
    TF_Bonf = siglist()$TF_Bonf1
    names(TF_Bonf) <- gsub("Bonf.P.Value.for.","",names(TF_Bonf))
    if(input$showfc == TRUE){
      sigcomp <- siglist()$sigcomp5
      names(sigcomp) <- gsub("Estimate.of.", "", names(sigcomp))
    }

    venn.Data <- list()
    for(i in 1:length(input$Vcomparison)){
      if(input$correction_method1 == "Raw"){
        if(input$showfc == TRUE){
          venn.Data[[i]] = results_file[,1][which(TF_Raw[[input$Vcomparison[i]]] == TRUE & sigcomp[[input$Vcomparison[i]]] == TRUE)]
        }
        if(input$showfc == FALSE){
          venn.Data[[i]] = results_file[,1][which(TF_Raw[[input$Vcomparison[i]]] == TRUE)]
        }
      }
      if(input$correction_method1 == "FDR"){
        if(input$showfc == TRUE){
          venn.Data[[i]] = results_file[,1][which(TF_FDR[[input$Vcomparison[i]]] == TRUE & sigcomp[[input$Vcomparison[i]]] == TRUE)]
        }
        if(input$showfc == FALSE){
          venn.Data[[i]] = results_file[,1][which(TF_FDR[[input$Vcomparison[i]]] == TRUE)]
        }
      }
      if(input$correction_method1 == "Bonferroni"){
        if(input$showfc == TRUE){
          venn.Data[[i]] = results_file[,1][which(TF_Bonf[[input$Vcomparison[i]]] == TRUE & sigcomp[[input$Vcomparison[i]]] == TRUE)]
        }
        if(input$showfc == FALSE){
          venn.Data[[i]] = results_file[,1][which(TF_Bonf[[input$Vcomparison[i]]] == TRUE)]
        }
      }
    }
    names(venn.Data) = input$Vcomparison
    if(length(venn.Data) > 5){
      venn.Data = venn.Data[1:5]
    }
    return(venn.Data)
  })

  Venn.intersection <- reactive({
    if(is.null(venndata())){return(NULL)}
    venndata1 <- venndata()[which(names(venndata()) %in% input$Include)]
    intersections <- Reduce(intersect, venndata1)
    n = length(input$Exclude)
    if(n > 0){
      venndata2 <- venndata()[which(names(venndata()) %in% input$Exclude)]
      venndata2 = unlist(venndata2)
      excl <- list(intersections, venndata2)
      intersections <- Reduce(setdiff, excl)
    }
    return(intersections)
  })

  Venn.union <- reactive({
    if(is.null(venndata())){return(NULL)}
    venndata1 <- venndata()[which(names(venndata()) %in% input$Include)]
    unions <- Reduce(union, venndata1)
    n = length(input$Exclude)
    if(n > 0){
      venndata2 <- venndata()[which(names(venndata()) %in% input$Exclude)]
      venndata2 <- unlist(venndata2)
      excl <- list(unions, venndata2)
      unions <- Reduce(setdiff, excl)
    }
    return(unions)
  })

  genelist2 <- reactive({
    n = length(input$Include)
    comp <- input$Include
    fcomp <- c()
    for(i in 1:n){
      fcomp[i] <- paste0("Estimate of ", input$Include[i],
                         "|Test.statistic for ", input$Include[i],
                         "|P.Value for ", input$Include[i])
    }
    if(input$UorI == 1){
      y <- which(dataset()[,1] %in% Venn.intersection())
      cols <- list()
      for(i in 1:n){
        cols[[i]] <- grep(fcomp[i], colnames(dataset()))
      }
    }
    if(input$UorI == 2){
      y <- which(dataset()[,1] %in% Venn.union())
      cols <- list()
      for(i in 1:n){
        cols[[i]] <- grep(fcomp[i], colnames(dataset()))
      }
    }
    x = list()
    pcols = list()
    pvals_fdr = list()
    pvals_bonf = list()
    pvals <- list()
    estimates <- list()
    tstats <- list()
    for(i in 1:n){
      if(i == 1){
        x[[i]] <- dataset()[y,c(1,2,cols[[i]])]
      }
      if(i > 1){
        x[[i]] <- dataset()[y,cols[[i]]]
      }
      pvals[[i]] <- grep("^P.Value", colnames(x[[i]]))
      estimates[[i]] <- grep("^Estimate", colnames(x[[i]]))
      tstats[[i]] <- grep("^Test.statistic", colnames(x[[i]]))
      if(i == 1){
        x[[i]] <- x[[i]][,c(1,2,estimates[[i]],pvals[[i]])]
      }
      if(i > 1){
        x[[i]] <- x[[i]][,c(estimates[[i]],pvals[[i]])]
      }
    }
    for(i in 1:n){
      if(i == 1){
        colnames(x[[i]]) = c("PROBE_ID", "SYMBOL", "Log2FC","P.Value")
      }
      if(i > 1){
        colnames(x[[i]]) = c("Log2FC", "P.Value")
      }
    }
    if(n > 1){
      for(i in 1:n){
        if(i == 1){
          colnames(x[[i]]) = c("PROBE_ID", "SYMBOL", paste("Log2FC for", comp[i], sep = " "), paste("P.Value for", comp[i], sep =" "))
        }
        if(i > 1){
          colnames(x[[i]]) = c(paste("Log2FC for", comp[i], sep = " "), paste("P.Value for", comp[i], sep = " "))
        }
      }
    }
    x = do.call("cbind", x)
    if(nrow(x) == 0){
      x <- data.frame("No Genes Present")
      names(x) <- "PROBE_ID"
    }
    return(x)
  })

  output$venn.intersection <- renderDataTable({
    y <- genelist2()
    y$SYMBOL <- paste("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",y$SYMBOL," target = '_blank'",'>',y$SYMBOL,"</a>",sep='')
    escape = FALSE
    y
  }, escape = FALSE)

  output$vennDiagram <- renderPlot({
    dummy <- data.frame("No Genes Present")
    names(dummy) <- "PROBE_ID"
    if(identical(genelist2(),dummy)){return(NULL)}
    color.choices = c("blue", "green", "red", "orange","purple")
    if(length(venndata()) == 1){
      color.choices = color.choices[1]
    }
    if(length(venndata()) == 2){
      color.choices = color.choices[1:2]
    }
    if(length(venndata()) == 3){
      color.choices = color.choices[1:3]
    }
    if(length(venndata()) == 4){
      color.choices = color.choices[1:4]
    }
    if(length(venndata()) == 5){
      color.choices = color.choices
    }
    grid.draw(VennDiagram::venn.diagram(venndata(), filename = NULL, lwd = 1, col = color.choices, cat.cex = .9, fil = color.choices, margin = .12, ext.text = FALSE,
                           euler.d = TRUE))
  })

  output$downloadVennData <- downloadHandler(
    filename = function() {paste0("VennDiagram_", input$UorI, '.csv')  },
    content = function(file) {
      dummy <- data.frame("No Genes Present")
      names(dummy) <- "PROBE_ID"
      if(identical(genelist2(),dummy)){return(NULL)}
      write.csv(genelist2(), file,row.names = FALSE)
    }
  )

  output$downloadVennPic <- downloadHandler(
    filename = function() {paste0("VennDiagram_", input$UorI, '.png')  },
    content = function(file) {
      dummy <- data.frame("No Genes Present")
      names(dummy) <- "PROBE_ID"
      if(identical(genelist2(),dummy)){return(NULL)}
      color.choices = c("blue", "green", "red", "orange","purple")
      if(length(venndata()) == 1){
        color.choices = color.choices[1]
      }
      if(length(venndata()) == 2){
        color.choices = color.choices[1:2]
      }
      if(length(venndata()) == 3){
        color.choices = color.choices[1:3]
      }
      if(length(venndata()) == 4){
        color.choices = color.choices[1:4]
      }
      if(length(venndata()) == 5){
        color.choices = color.choices
      }
      png(file)
      grid.draw(VennDiagram::venn.diagram(venndata(), filename = NULL, lwd = 1, col = color.choices, cat.cex = .9, fil = color.choices, margin = .12, ext.text = FALSE,euler.d = TRUE))
      dev.off()
    }
  )

  output$LMMmodule_text <- renderUI({
    if(all(lmm_module()$dat == 0)){
      return(helpText(strong("All module percentages are zero, so no plot is given.")))
    }
    else{
      return(NULL)
    }
  })

  output$LMMmodule<-renderPlot({
    if(min(lmm_module()$dat) >= 0 & max(lmm_module()$dat) > 0){
      aheatmap2(lmm_module()$dat, Rowv = lmm_module()$rowv, circle_size = input$LMMradius, fontsize = input$LMMFont, Colv = NA , color = lmm_module()$color_palette[500:1000])
    }
    if(max(lmm_module()$dat) <= 0 & min(lmm_module()$dat) < 0){
      aheatmap2(lmm_module()$dat, Rowv = lmm_module()$rowv, circle_size = input$LMMradius, fontsize = input$LMMFont, Colv = NA , color = lmm_module()$color_palette[1:500])
    }
    if(min(lmm_module()$dat) < 0 & max(lmm_module()$dat) >0){
      aheatmap2(lmm_module()$dat, Rowv = lmm_module()$rowv, circle_size = input$LMMradius, fontsize = input$LMMFont, Colv = NA , color = lmm_module()$color_palette)
    }
    if(all(lmm_module()$dat == 0)){
      return(NULL)
    }
  },height = LMMmap2, width = LMMmap1)

  output$downloadLMMModMap2 <- downloadHandler(
    filename = function() {paste('LMMModMap','.png', sep = '')},
    content = function(file){
      png(file, width = (resolution1()/72)*LMMmap1(), height = (resolution1()/72)*LMMmap2(), res = resolution1())
      if(min(lmm_module()$dat) >= 0 & max(lmm_module()$dat) > 0){
        print(aheatmap2(lmm_module()$dat, Rowv = lmm_module()$rowv, circle_size = input$LMMradius, fontsize = input$LMMFont, Colv = NA , color = lmm_module()$color_palette[500:1000]))
      }
      if(max(lmm_module()$dat) <= 0 & min(lmm_module()$dat) < 0){
        print(aheatmap2(lmm_module()$dat, Rowv = lmm_module()$rowv, circle_size = input$LMMradius, fontsize = input$LMMFont, Colv = NA , color = lmm_module()$color_palette[1:500]))
      }
      if(min(lmm_module()$dat) < 0 & max(lmm_module()$dat) >0){
        print(aheatmap2(lmm_module()$dat, Rowv = lmm_module()$rowv, circle_size = input$LMMradius, fontsize = input$LMMFont, Colv = NA , color = lmm_module()$color_palette))
      }
      if(all(lmm_module()$dat == 0)){
        return(NULL)
      }
      dev.off()
    }
  )

  output$modmap_comparison2 <- renderUI({
    if(is.null(values$ModulesTF) == FALSE){
      if(values$ModulesTF == FALSE){return(NULL)}
      if(values$ModulesTF == TRUE){return(plotOutput("modmap_comparison"))}
    }
    else{
      return(plotOutput("modmap_comparison"))
    }
  })


ModMap_Comparison <- reactive({
    numgen <- function(x){
      d<-c()
      odd<-2*(1:20)-1
      if(length(x)==1){
        return(odd[1:x])}
      if(length(x)>1){
        for(i in 1:length(x)){d<-c(d,odd[1:x[i]])}
        return(d)}
    }
    numgen2 <- function(x){
      d<-c()
      odd<-2*(1:8)-1
      for (i in 1:length(x)){d<-c(d,rep(odd[9-i],x[i]))}
      return(d)
    }
    pcomp <- pcomp()
    pcomp <- gsub(" ", ".", pcomp)
    colorvar = sig_ind()$prop_matrix[1:97,c(pcomp)]+1
    colgrp <- findInterval(colorvar,seq(0,2,length.out=10))
    colfunc <- colorRampPalette(c("blue","white", "red"))
    collist <- colfunc(length(unique(colgrp)))
    mycolors <- collist[colgrp]
    mycolors = colorvar
    df2 <- data.frame(
      x = numgen(c(2,3,6,16,15,20,20,15)),
      y = c(numgen2(c(2,3,6,16,15,20,20,15))),
      z = mycolors
    )

    test<-ggplot(df2, aes(xmin=x-1,xmax=x+1,ymin=y-1,ymax=y+1),environment=environment())+
      theme(axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),panel.grid.major = element_blank()
            ,panel.grid.minor = element_blank()
            ,panel.border = element_blank())+geom_rect(fill="white", colour="black")+
      scale_y_continuous(breaks=c(15,13,11,9,7,5,3,1), labels=c("M1","M2","M3","M4","M5","M6","M7","M7 (21-35)"),limits=c(0,40))+
      scale_x_continuous(breaks=(2*(1:20)-1),labels=1:20,limits=c(0,40))

    test2<-test+coord_cartesian(xlim = c(0, 40), ylim=c(0, 16))
    test3<-test2+geom_point(aes(x=x,y=y,size=500,colour=z))+
      scale_colour_gradient2(low="blue", high="red", guide="colorbar",midpoint=1)+
      theme(axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position = "none") + scale_size_continuous(limits=c(1,500))

    if(is.null(values$ModulesTF) == FALSE){
      if(values$ModulesTF == FALSE){return(NULL)}
      if(values$ModulesTF == TRUE){return(test3)}
    }
    else{
      return(test3)
    }
  })
  mod_map_lmm_res <- reactive({
    input$modmaplmmres
  })


  output$modmap_comparison<-renderPlot({
    if(is.null(ModMap_Comparison())){return(NULL)}

    withProgress(message = 'Making the table',
                 detail = 'This may take a while...', value = 1,{
                   ModMap_Comparison()
                 })

  }, height = 190, width = 700)

  output$download_modmap_comparison <- downloadHandler(
    filename = function() {paste('ModMap_Comparison','.png', sep = '')},
    content = function(file){
      png(file, width = (mod_map_lmm_res()/72)*700, height = (mod_map_lmm_res()/72)*190, res = mod_map_lmm_res())
      print(ModMap_Comparison())
      dev.off()
    }
  )


############## QUSAGE #####################

  output$qusage <- renderMenu({
    if(is.null(values$qusage_results)){
      return(strong(""))
    }
    if(is.null(values$qusage_results) == FALSE){
      menuItem("Q-Gen", icon = icon("th-list"), tabName = "qusage")
    }
  })

  output$genesets <- renderUI({
    if(is.null(values$qusage_results)){
      return(NULL)
    }
    else{
     if(is.list(values$qusage_results) & !is.data.frame(values$qusage_results)){
       if(is.null(names(values$qusage_results))){
         return(selectInput("genesets1", "Select geneset definition:", choices = c(1:length(values$qusage_results))))
       }
       if(!is.null(names(values$qusage_results))){
         return(selectInput("genesets1", "Select geneset definition:", choices = c(names(values$qusage_results))))
       }
     }
      else{
        return(NULL)
      }
    }
  })

  master <- reactive({
    if(is.null(values$qusage_results)){return(NULL)}
    if(is.null(input$genesets1)){
      qusage_results <- values$qusage_results
    }
    if(!is.null(input$genesets1)){
      if(!is.null(names(values$qusage_results))){
        qusage_results <- values$qusage_results[[input$genesets1]]
      }
      else{
        qusage_results <- values$qusage_results[[as.numeric(input$genesets1)]]
      }
    }
    master = qusage_results
    colnames(master)[which(colnames(master) == "log.fold.change")] <- "Log2FC"
    unique_comp <- unique(master$Comparison)
    sub_uniquecomp <- list()
    n = length(unique_comp)
    for(i in 1:n){
      sub_uniquecomp[[i]] <- which(master$Comparison == unique_comp[i])
    }
    bonf <- list()
    for(i in 1:n){
      bonf[[i]] <- p.adjust(master$p.Value[sub_uniquecomp[[i]]], method = "bonferroni")
    }
    bonf <- unlist(bonf)
    master$Bonf <- bonf
    if(!is.null(values$illumina)){
      if(values$illumina == FALSE){
        if(is.null(values$moduleinfo2)){
          master$pathway.name <- gsub("_", ".", master$pathway.name)
        }
      }
    }
    n2 <- length(module_annotations[,1])
    mod.ann <- list()
    master$Modulev2_Annotation <- ""
    for(i in 1:n2){
      mod.ann[[i]] <- which(master$pathway.name %in% module_annotations[,1][i])
    }
    for(i in 1:n2){
      master$Modulev2_Annotation[mod.ann[[i]]] <- as.character(module_annotations[,2][i])
    }
    if(is.null(values$BaylorTF)){
      mod.annot <- module_annotations
      mod.annot <- mod.annot[match(master$pathway.name[which(master$pathway.name %in% mod.annot$Module)], mod.annot$Module, nomatch = 0),]
      master$pathway.name <- as.character(master$pathway.name)
      master$pathway.name[which(master$pathway.name %in% mod.annot$Module)] <- paste(mod.annot$Modulev2_Annotation,mod.annot$Module,sep=" ")
      master$pathway.name <- as.factor(master$pathway.name)
    }
    return(master)
  })

  mod_list <- reactive({
    if(is.null(values$illumina) == FALSE){
      if(values$illumina == FALSE){
        if(is.null(values$moduleinfo2)){
          modnames<-unique(moduleinfo2$Module)
          y<-unlist(strsplit(as.character(modnames),".",fixed=TRUE))
          y <- unlist(modnames)
          y<-matrix(y,byrow=T,nrow=260,ncol=2)
          sortmodnames<-as.character(modnames)[order(y[,1],y[,2])]
          Modulelist<-list(as.character(moduleinfo2$SYMBOL[which(moduleinfo2$Module==sortmodnames[1])]),as.character(moduleinfo2$SYMBOL[which(moduleinfo2$Module==sortmodnames[2])]))
          for(i in 3:length(sortmodnames)){Modulelist<-c(Modulelist,list(as.character(moduleinfo2$SYMBOL[which(moduleinfo2$Module==sortmodnames[i])])))}
          names(Modulelist)<-sortmodnames
        }
        else{
          modnames<-unique(values$moduleinfo2$Module)
          y<-unlist(strsplit(as.character(modnames),".",fixed=TRUE))
          y<-matrix(y,byrow=T,nrow=260,ncol=2)
          sortmodnames<-as.character(modnames)[order(y[,1],y[,2])]
          Modulelist<-list(as.character(values$moduleinfo2$SYMBOL[which(values$moduleinfo2$Module==sortmodnames[1])]),as.character(values$moduleinfo2$SYMBOL[which(values$moduleinfo2$Module==sortmodnames[2])]))
          for(i in 3:length(sortmodnames)){Modulelist<-c(Modulelist,list(as.character(values$moduleinfo2$SYMBOL[which(values$moduleinfo2$Module==sortmodnames[i])])))}
          names(Modulelist)<-sortmodnames
        }
      }
      else{
        Modulelist <- Modulelist
        mod.anno <- names(Modulelist)[which(names(Modulelist) %in% module_annotations[,1])]
        module_annotations2 <- module_annotations[match(mod.anno, module_annotations[,1], nomatch=0),]
        names(Modulelist)[which(names(Modulelist) %in% module_annotations2[,1])] <- paste(module_annotations2[,2], module_annotations2[,1], sep=" ")
      }
    }
    if(is.null(values$illumina)){
      Modulelist <- Modulelist
      mod.anno <- names(Modulelist)[which(names(Modulelist) %in% module_annotations[,1])]
      module_annotations2 <- module_annotations[match(mod.anno, module_annotations[,1], nomatch=0),]
      names(Modulelist)[which(names(Modulelist) %in% module_annotations2[,1])] <- paste(module_annotations2[,2], module_annotations2[,1], sep=" ")
    }
    return(Modulelist)
  })


  output$Comparisons <- renderUI({
    selectInput("set_comp", "Comparison selection:", choices = levels(master()$Comparison), selected = levels(master()$Comparison)[1])
  })

  output$Comparisons1 <- renderUI({
    selectInput("set_comp1", "Comparison(s) selection:", choices = levels(master()$Comparison), selected = levels(master()$Comparison)[1], multiple = T)
  })

  output$Comparisons2 <- renderUI({
    selectInput("set_comp2", "Comparison selection:", choices = levels(master()$Comparison), selected = levels(master()$Comparison)[1])
  })

  output$Comparisons3 <- renderUI({
    selectInput("set_comp3", "Comparison selection:", choices = levels(master()$Comparison), selected = levels(master()$Comparison)[1])
  })

  output$Module_Select1 <- renderUI({
    selectInput("Mod_Select1", "Module selection:", choices = names(mod_list()), selected = names(mod_list()[1]))
  })

  output$PaloOrFirst1  <- renderUI({
    if(length(input$set_comp1) <= 1){return(NULL)}
    else{
      selectInput("PaloOrFirst", "Significant gene sets:", choices = c("In at least one comparison chosen" = 1, "In first comparison chosen" = 2), selected = 1)
    }
  })

  output$ColorChoice <- renderUI({
    if(length(input$set_comp1) <= 1){return(NULL)}
    else{
      cols = c("blue", "red","green","lightskyblue","purple","hotpink","brown","gold")
      selectInput("ColorChoice1", "Choose line colors:", choices = cols, selected = cols[1:length(input$set_comp1)], multiple = T)
    }
  })

  multi_testing <- reactive({
    Multi_testing1 = input$Multi_testing1
    Multi_testing1
  })

  multi_testing2 <- reactive({
    Multi_testing2 = input$Multi_testing2
    Multi_testing2
  })

  comp_subset1 <- reactive({
    if(is.null(master())){return(NULL)}
    if(length(input$set_comp1) < 1){return(NULL)}
    master = master()[which(master()$Comparison %in% input$set_comp1),]
    first_comp = master[which(master$Comparison %in% input$set_comp1[1]),]
    first_cut = unique(first_comp[which(first_comp$Log2FC <= -input$FilterOnFoldchange1 | first_comp$Log2FC >= input$FilterOnFoldchange1), ]$pathway.name)
    palo_cut = unique(master[which(master$Log2FC <= -input$FilterOnFoldchange1 | master$Log2FC >= input$FilterOnFoldchange1), ]$pathway.name)
    if(length(input$set_comp1) == 1 || input$PaloOrFirst == 1){
      master = master[which(master$pathway.name %in% palo_cut), ]
    }
    else{
      master = master[which(master$pathway.name %in% first_cut), ]
    }
    if(multi_testing() == "Raw"){
      if(length(input$set_comp1) == 1 || input$PaloOrFirst == 1){
        palo = unique(master[which(master$p.Value <= input$FilterOnPValues1),]$pathway.name)
        if(length(palo) == 0){return(NULL)}
        master1 = master[which(master$Comparison %in% input$set_comp1[1]),]
        master1 = master1[which(master1$pathway.name %in% palo), ]
        index <- order(master1$Log2FC, decreasing = T)
        newindex <- rep(index, length(input$set_comp1))
        if(length(input$set_comp1) == 1){
          shift <- rep(0, each = length(palo))
        }
        if(length(input$set_comp1) == 2){
          shift <- rep(c(0, length(palo)), each = length(palo))
        }
        if(length(input$set_comp1) == 3){
          shift <- rep(c(0, length(palo), 2*length(palo)), each = length(palo))
        }
        if(length(input$set_comp1) == 4){
          shift <- rep(c(0, length(palo), 2*length(palo), 3*length(palo)), each = length(palo))
        }
        if(length(input$set_comp1) == 5){
          shift <- rep(c(0, length(palo), 2*length(palo), 3*length(palo), 4*length(palo)), each = length(palo))
        }
        newindex <- newindex + shift
        master = master[which(master$pathway.name %in% palo),]
        master = master[newindex,]
        master = cbind(1:length(newindex), master)
        names(master) = c("Index", names(master[,-1]))
        master[,2] = factor(master[,2], levels = unique(master[,2]))
        master$Index = rep(1:length(palo), length(input$set_comp1))
        sig = which(master$p.Value <= input$FilterOnPValues1)
        master$SIG = "Not Significant"
        master$SIG[sig] = "Significant"
        return(master)
      }

      if(input$PaloOrFirst == 2){
        master1 = master[which(master$Comparison %in% input$set_comp1[1]),]
        first = unique(master1[which(master1$p.Value <= input$FilterOnPValues1),]$pathway.name)
        if(length(first) == 0){return(NULL)}
        master1 = master1[which(master1$pathway.name %in% first), ]
        index <- order(master1$Log2FC, decreasing = T)
        newindex <- rep(index, length(input$set_comp1))
        if(length(input$set_comp1) == 1){
          shift <- rep(0, each = length(first))
        }
        if(length(input$set_comp1) == 2){
          shift <- rep(c(0, length(first)), each = length(first))
        }
        if(length(input$set_comp1) == 3){
          shift <- rep(c(0, length(first), 2*length(first)), each = length(first))
        }
        if(length(input$set_comp1) == 4){
          shift <- rep(c(0, length(first), 2*length(first), 3*length(first)), each = length(first))
        }
        if(length(input$set_comp1) == 5){
          shift <- rep(c(0, length(first), 2*length(first), 3*length(first), 4*length(first)), each = length(first))
        }
        newindex <- newindex + shift
        master = master[which(master$pathway.name %in% first),]
        master = master[newindex,]
        master = cbind(1:length(newindex), master)
        names(master) = c("Index", names(master[,-1]))
        master[,2] = factor(master[,2], levels = unique(master[,2]))
        master$Index = rep(1:length(first), length(input$set_comp1))
        sig = which(master$p.Value <= input$FilterOnPValues1)
        master$SIG = "Not Significant"
        master$SIG[sig] = "Significant"
        return(master)
      }
    }

    if(multi_testing() == "FDR"){
      if(length(input$set_comp1) == 1 || input$PaloOrFirst == 1){
        palo = unique(master[which(master$FDR <= input$FilterOnPValues1),]$pathway.name)
        if(length(palo) == 0){return(NULL)}
        master1 = master[which(master$Comparison %in% input$set_comp1[1]), ]
        master1 = master1[which(master1$pathway.name %in% palo), ]
        index <- order(master1$Log2FC, decreasing = T)
        newindex <- rep(index, length(input$set_comp1))
        if(length(input$set_comp1) == 1){
          shift <- rep(0, each = length(palo))
        }
        if(length(input$set_comp1) == 2){
          shift <- rep(c(0, length(palo)), each = length(palo))
        }
        if(length(input$set_comp1) == 3){
          shift <- rep(c(0, length(palo), 2*length(palo)), each = length(palo))
        }
        if(length(input$set_comp1) == 4){
          shift <- rep(c(0, length(palo), 2*length(palo), 3*length(palo)), each = length(palo))
        }
        if(length(input$set_comp1) == 5){
          shift <- rep(c(0, length(palo), 2*length(palo), 3*length(palo), 4*length(palo)), each = length(palo))
        }
        newindex <- newindex + shift
        master = master[which(master$pathway.name %in% palo),]
        master = master[newindex,]
        master <- cbind(1:length(newindex), master)
        names(master) = c("Index", names(master[,-1]))
        master[,2] = factor(master[,2], levels = unique(master[,2]))
        master$Index = rep(1:length(palo), length(input$set_comp1))
        sig = which(master$FDR <= input$FilterOnPValues1)
        master$SIG = "Not Significant"
        master$SIG[sig] = "Significant"
        return(master)
      }

      if(input$PaloOrFirst == 2){
        master1 = master[which(master$Comparison %in% input$set_comp1[1]),]
        first = unique(master1[which(master1$FDR <= input$FilterOnPValues1),]$pathway.name)
        if(length(first) == 0){return(NULL)}
        master1 = master1[which(master1$pathway.name %in% first), ]
        index <- order(master1$Log2FC, decreasing = T)
        newindex <- rep(index, length(input$set_comp1))
        if(length(input$set_comp1) == 1){
          shift <- rep(0, each = length(first))
        }
        if(length(input$set_comp1) == 2){
          shift <- rep(c(0, length(first)), each = length(first))
        }
        if(length(input$set_comp1) == 3){
          shift <- rep(c(0, length(first), 2*length(first)), each = length(first))
        }
        if(length(input$set_comp1) == 4){
          shift <- rep(c(0, length(first), 2*length(first), 3*length(first)), each = length(first))
        }
        if(length(input$set_comp1) == 5){
          shift <- rep(c(0, length(first), 2*length(first), 3*length(first), 4*length(first)), each = length(first))
        }
        newindex <- newindex + shift
        master = master[which(master$pathway.name %in% first),]
        master = master[newindex,]
        master = cbind(1:length(newindex), master)
        names(master) = c("Index", names(master[,-1]))
        master[,2] = factor(master[,2], levels = unique(master[,2]))
        master$Index = rep(1:length(first), length(input$set_comp1))
        sig = which(master$FDR <= input$FilterOnPValues1)
        master$SIG = "Not Significant"
        master$SIG[sig] = "Significant"
        return(master)
      }
    }

    if(multi_testing() == "Bonferroni"){
      if(length(input$set_comp1) == 1 || input$PaloOrFirst == 1){
        palo = unique(master[which(master$Bonf <= input$FilterOnPValues1),]$pathway.name)
        if(length(palo) == 0){return(NULL)}
        master1 = master[which(master$Comparison %in% input$set_comp1[1]), ]
        master1 = master1[which(master1$pathway.name %in% palo), ]
        index <- order(master1$Log2FC, decreasing = T)
        newindex <- rep(index, length(input$set_comp1))
        if(length(input$set_comp1) == 1){
          shift <- rep(0, each = length(palo))
        }
        if(length(input$set_comp1) == 2){
          shift <- rep(c(0, length(palo)), each = length(palo))
        }
        if(length(input$set_comp1) == 3){
          shift <- rep(c(0, length(palo), 2*length(palo)), each = length(palo))
        }
        if(length(input$set_comp1) == 4){
          shift <- rep(c(0, length(palo), 2*length(palo), 3*length(palo)), each = length(palo))
        }
        if(length(input$set_comp1) == 5){
          shift <- rep(c(0, length(palo), 2*length(palo), 3*length(palo), 4*length(palo)), each = length(palo))
        }
        newindex <- newindex + shift
        master = master[which(master$pathway.name %in% palo),]
        master = master[newindex,]
        master <- cbind(1:length(newindex), master)
        names(master) = c("Index", names(master[,-1]))
        master[,2] = factor(master[,2], levels = unique(master[,2]))
        master$Index = rep(1:length(palo), length(input$set_comp1))
        sig = which(master$Bonf <= input$FilterOnPValues1)
        master$SIG = "Not Significant"
        master$SIG[sig] = "Significant"
        return(master)
      }

      if(input$PaloOrFirst == 2){
        master1 = master[which(master$Comparison %in% input$set_comp1[1]),]
        first = unique(master1[which(master1$Bonf <= input$FilterOnPValues1),]$pathway.name)
        if(length(first) == 0){return(NULL)}
        master1 = master1[which(master1$pathway.name %in% first), ]
        index <- order(master1$Log2FC, decreasing = T)
        newindex <- rep(index, length(input$set_comp1))
        if(length(input$set_comp1) == 1){
          shift <- rep(0, each = length(first))
        }
        if(length(input$set_comp1) == 2){
          shift <- rep(c(0, length(first)), each = length(first))
        }
        if(length(input$set_comp1) == 3){
          shift <- rep(c(0, length(first), 2*length(first)), each = length(first))
        }
        if(length(input$set_comp1) == 4){
          shift <- rep(c(0, length(first), 2*length(first), 3*length(first)), each = length(first))
        }
        if(length(input$set_comp1) == 5){
          shift <- rep(c(0, length(first), 2*length(first), 3*length(first), 4*length(first)), each = length(first))
        }
        newindex <- newindex + shift
        master = master[which(master$pathway.name %in% first),]
        master = master[newindex,]
        master = cbind(1:length(newindex), master)
        names(master) = c("Index", names(master[,-1]))
        master[,2] = factor(master[,2], levels = unique(master[,2]))
        master$Index = rep(1:length(first), length(input$set_comp1))
        sig = which(master$Bonf <= input$FilterOnPValues1)
        master$SIG = "Not Significant"
        master$SIG[sig] = "Significant"
        return(master)
      }
    }
  })

output$Module_Select <- renderUI({
    master <- master()[which(master()$Comparison %in% input$set_comp2),]
    pathways <- master$pathway.name
    selectInput("Mod_Select", "Module selection:", choices = gtools::mixedsort(as.character(pathways)), selected = gtools::mixedsort(as.character(pathways))[1])
  })

  venn.data <- reactive({
    if(is.null(master())){return(NULL)}
    master <- comp_subset1()
    unique.comparison <- unique(master$Comparison)
    n = length(unique.comparison)
    master1 <- list()
    for(i in 1:n){
      master1[[i]] <- subset(master, Comparison == unique.comparison[i])
    }
    names(master1) <- as.character(unique.comparison)
    for(i in 1:n){
      master1[[i]] <- master1[[i]]$pathway.name[which(master1[[i]]$SIG == "Significant")]
    }
    master1
  })

  comp_overview <- reactive({
    if(is.null(master())){return(NULL)}
    master1 = master()[which(master()$p.Value <= input$SigLevel),]
    master2 = master()[which(master()$FDR <= input$SigLevel),]
    master3 = master()[which(master()$Bonf <= input$SigLevel),]
    if(nrow(master1) == 0 & nrow(master2) == 0 & nrow(master3) == 0){return(NULL)}
    if(input$foldchange.q == "+"){
      master1 = master1[which(master1$Log2FC > 0),]
      master2 = master2[which(master2$Log2FC > 0),]
      master3 = master3[which(master3$Log2FC > 0),]
    }
    if(input$foldchange.q == "-"){
      master1 = master1[which(master1$Log2FC < 0),]
      master2 = master2[which(master2$Log2FC < 0),]
      master3 = master3[which(master3$Log2FC < 0),]
    }
    if(nrow(master1) > 0){
      master_sig1 = aggregate(Comparison ~ Comparison, data = master1, FUN = table)
      master_sig1 = as.matrix(master_sig1)
      master_sig1 = as.data.frame(master_sig1)
      colnames(master_sig1) = sub("Comparison.", "",colnames(master_sig1))
      master_sig1 = t(master_sig1)
      master_sig1 = as.data.frame(master_sig1)
      master_sig1 = cbind(rownames(master_sig1), master_sig1[,1])
      master_sig1 = as.data.frame(master_sig1)
    }
    if(nrow(master2) > 0){
      master_sig2 = aggregate(Comparison ~ Comparison, data = master2, FUN = table)
      master_sig2 = as.matrix(master_sig2)
      master_sig2 = as.data.frame(master_sig2)
      colnames(master_sig2) = sub("Comparison.", "",colnames(master_sig2))
      master_sig2 = t(master_sig2)
      master_sig2 = as.data.frame(master_sig2)
      master_sig2 = cbind(rownames(master_sig2), master_sig2[,1])
      master_sig2 = as.data.frame(master_sig2)
    }
    if(nrow(master3) > 0){
      master_sig3 = aggregate(Comparison ~ Comparison, data = master3, FUN = table)
      master_sig3 = as.matrix(master_sig3)
      master_sig3 = as.data.frame(master_sig3)
      colnames(master_sig3) = sub("Comparison.", "",colnames(master_sig3))
      master_sig3 = t(master_sig3)
      master_sig3 = as.data.frame(master_sig3)
      master_sig3 = cbind(rownames(master_sig3), master_sig3[,1])
      master_sig3 = as.data.frame(master_sig3)
    }
    if(nrow(master1) > 0 & nrow(master2) > 0 & nrow(master3) > 0){
      master_sig = cbind(master_sig1, master_sig2[,-1], master_sig3[,-1])
      master_sig[,2] = as.numeric(as.character(master_sig[,2]))
      master_sig[,3] = as.numeric(as.character(master_sig[,3]))
      master_sig[,4] = as.numeric(as.character(master_sig[,4]))
      colnames(master_sig) = c("Comparison", "Raw", "FDR", "Bonf")
      master_sig = master_sig[order(master_sig$Raw, decreasing = T),]
    }
    if(nrow(master1) > 0 & nrow(master2) > 0 & nrow(master3) == 0){
      master_sig = cbind(master_sig1, master_sig2[,-1], 0)
      master_sig[,2] = as.numeric(as.character(master_sig[,2]))
      master_sig[,3] = as.numeric(as.character(master_sig[,3]))
      master_sig[,4] = as.numeric(as.character(master_sig[,4]))
      colnames(master_sig) = c("Comparison", "Raw", "FDR", "Bonf")
      master_sig = master_sig[order(master_sig$Raw, decreasing = T),]
    }
    if(nrow(master1) > 0 & nrow(master2) == 0 & nrow(master3) > 0){
      master_sig = cbind(master_sig1, 0, master_sig3[,-1])
      master_sig[,2] = as.numeric(as.character(master_sig[,2]))
      master_sig[,3] = as.numeric(as.character(master_sig[,3]))
      master_sig[,4] = as.numeric(as.character(master_sig[,4]))
      colnames(master_sig) = c("Comparison", "Raw", "FDR", "Bonf")
      master_sig = master_sig[order(master_sig$Raw, decreasing = T),]
    }
    if(nrow(master1) > 0 & nrow(master2) == 0 & nrow(master3) == 0){
      master_sig = cbind(master_sig1, 0, 0)
      master_sig[,2] = as.numeric(as.character(master_sig[,2]))
      master_sig[,3] = as.numeric(as.character(master_sig[,3]))
      master_sig[,4] = as.numeric(as.character(master_sig[,4]))
      colnames(master_sig) = c("Comparison", "Raw", "FDR", "Bonf")
      master_sig = master_sig[order(master_sig$Raw, decreasing = T),]
    }
    return(master_sig)
  })

  individual_modselect <- reactive({
    if(is.null(master())){return(NULL)}
    lowerCI = values$lowerCI
    upperCI = values$upperCI
    if(is.null(values$results_file2) == FALSE){
      results_file = values$results_file2
    }
    else{
      results_file = dataset()
    }
    mmr_comp <- grep("^Estimate", colnames(results_file))
    mmr_comparisons <- results_file[,c(1,2,mmr_comp)]
    NAs <- apply(mmr_comparisons[,-c(1,2)], 2, function(x) is.na(x))
    NAs <- grep("TRUE", NAs)
    if(identical(NAs, integer(0))){
      mmr_comparisons = mmr_comparisons
      mmr_pvals <- grep("^P.Value", colnames(results_file))
      pvals <- results_file[,c(1,2,mmr_pvals)]
      colnames(pvals) <- sub("P.Value for ", "", colnames(pvals))
      rownames(pvals) <- rownames(mmr_comparisons)
      t_cols <- grep("^Test.statistic", colnames(results_file))
      tstats <- results_file[,c(1,2,t_cols)]
      std_error <- abs(mmr_comparisons[,-c(1,2)]/tstats[,-c(1,2)])
      colnames(std_error) <- paste("Std_error.", colnames(std_error))
    }
    else {
      mmr_comparisons = mmr_comparisons[-NAs,]
      mmr_pvals <- grep("^P.Value", colnames(results_file))
      pvals <- results_file[-NAs,c(1,2,mmr_pvals)]
      colnames(pvals) <- sub("P.Value for ", "", colnames(pvals))
      rownames(pvals) <- rownames(mmr_comparisons)
      t_cols <- grep("^Test.statistic", colnames(results_file))
      tstats <- results_file[-NAs,c(1,2,t_cols)]
      std_error <- abs(mmr_comparisons[,-c(1,2)]/tstats[,-c(1,2)])
      colnames(std_error) <- paste("Std_error.", colnames(std_error))
    }
    master = master()
    master = master[which(master$Comparison == input$set_comp2 & master$pathway.name == input$Mod_Select), ,drop = F]
    mmr_modselection <- mmr_comparisons
    mmr_modselection <- mmr_modselection[which(mmr_modselection[,1] %in% mod_list()[[input$Mod_Select]]),]
    mmr_names <- as.character(mmr_modselection[,2])
    mmr_probes <- as.character(mmr_modselection[,1])
    mod_pvals <- which(pvals[,1] %in% mod_list()[[input$Mod_Select]])
    lowerCI <- lowerCI[which(mmr_comparisons[,1] %in% mod_list()[[input$Mod_Select]]),]
    upperCI <- upperCI[which(mmr_comparisons[,1] %in% mod_list()[[input$Mod_Select]]),]
    colnames(mmr_modselection) <- sub("Estimate of ", "", colnames(mmr_modselection))
    mmr_modselection <- as.data.frame(mmr_modselection[,which(colnames(mmr_modselection) %in% input$set_comp2), drop = F])
    pvals2 <- as.data.frame(pvals[,which(colnames(pvals) %in% input$set_comp2)])
    fdr_pvals <- apply(pvals2, 2, p.adjust, method = "fdr")
    bonf_pvals <- apply(pvals2,2,p.adjust,method = "bonferroni")
    pvals <- pvals2[mod_pvals,]
    fdr_pvals <- fdr_pvals[mod_pvals,]
    bonf_pvals <- bonf_pvals[mod_pvals,]
    lowerCI <- as.data.frame(lowerCI[,which(colnames(lowerCI) %in% input$set_comp2)])
    upperCI <- as.data.frame(upperCI[,which(colnames(upperCI) %in% input$set_comp2)])
    missing_symb <- which(sapply(mmr_names, function(x) identical("",x)))
    if(length(missing_symb) > 0){
      mmr_names[missing_symb] = mmr_probes[missing_symb]
    }
    duplicates = which(duplicated(mmr_names))
    if(sum(duplicates) > 0){
      mmr_names[duplicates] = mmr_probes[duplicates]
    }
    mmr_modselection <- cbind(mmr_probes,mmr_names,mmr_modselection, pvals, fdr_pvals, bonf_pvals,lowerCI, upperCI)
    for(i in 3:ncol(mmr_modselection)){
      mmr_modselection[,i] <- as.numeric(mmr_modselection[,i])
    }
    colnames(mmr_modselection) <- c("PROBE_ID","SYMBOL","Log2FC", "P.Value", "FDR", "Bonf","low", "up")
    if(multi_testing2() == "Raw"){
      mmr_modselection <- mmr_modselection[which(mmr_modselection$P.Value <= input$FilterOnPValues2), ]
    }
    if(multi_testing2() == "FDR"){
      mmr_modselection <- mmr_modselection[which(mmr_modselection$FDR <= input$FilterOnPValues2), ]
    }
    if(multi_testing2() == "Bonferroni"){
      mmr_modselection <- mmr_modselection[which(mmr_modselection$Bonf <= input$FilterOnPValues2), ]
    }
    mmr_modselection <- mmr_modselection[which(mmr_modselection$Log2FC <= -input$FilterOnFoldchange2 | mmr_modselection$Log2FC >= input$FilterOnFoldchange2), ]
    if(nrow(mmr_modselection) == 0){return(NULL)}
    Index <- order(mmr_modselection$Log2FC, decreasing = T)
    mmr_modselection <- mmr_modselection[Index,]
    mmr_modselection <- cbind(1:length(Index), mmr_modselection)
    colnames(mmr_modselection) <- c("Index", colnames(mmr_modselection[,-1]))
    mmr_modselection <- mmr_modselection
    y = list(mmr_modselection = mmr_modselection, master = master)
    return(y)
  })

  plotHeight = function(){input$PlotHeight}
  plotWidth = function(){input$PlotWidth}
  plotHeight1 = function(){input$PlotHeight1}
  plotWidth1 = function(){input$PlotWidth1}
  plotHeight2 = function(){input$PlotHeight2}
  plotWidth2 = function(){input$PlotWidth2}
  plotres = function(){input$PlotRes}
  plotres1 = function(){input$PlotRes1}
  plotres2 = function(){input$PlotRes2}

  ranges <- reactiveValues(x = NULL, y = NULL)
  ranges2 <- reactiveValues(x = NULL, y = NULL)
  comp_plot1 <- reactive({
    if(is.null(master())){return(NULL)}
    if(is.null(comp_subset1())){return(NULL)}

    if(input$only_annotated){
      comp_subset1 <- comp_subset1()[which(comp_subset1()$pathway.name %in% paste(module_annotations$Modulev2_Annotation,module_annotations$Module,sep=" ")),]
    }
    else{
      comp_subset1 <- comp_subset1()
    }

    ylow = c()
    yhigh = c()

    if(min(comp_subset1$low) < 0 & max(comp_subset1$up) > 0){
      if(abs(min(comp_subset1$low)) >= max(comp_subset1$up)){
        ylow = min(comp_subset1$low)
        yhigh = -min(comp_subset1$low)
      }
      else{
        ylow = -max(comp_subset1$up)
        yhigh = max(comp_subset1$up)
      }
    }

    if(min(comp_subset1$low) > 0){
      ylow = -max(comp_subset1$up)
      yhigh = max(comp_subset1$up)
    }

    if(max(comp_subset1$up) < 0){
      ylow = min(comp_subset1$low)
      yhigh = -min(comp_subset1$low)
    }

    if(length(input$set_comp1) > 1){
      if(min(comp_subset1$Log2FC) < 0 & max(comp_subset1$Log2FC) > 0){
        if(abs(min(comp_subset1$Log2FC)) >= max(comp_subset1$Log2FC)){
          ylow = min(comp_subset1$Log2FC)
          yhigh = -min(comp_subset1$Log2FC)
        }
        else{
          ylow = -max(comp_subset1$Log2FC)
          yhigh = max(comp_subset1$Log2FC)
        }
      }
      if(min(comp_subset1$Log2FC) > 0){
        ylow = -max(comp_subset1$Log2FC)
        yhigh = max(comp_subset1$Log2FC)
      }
      if(max(comp_subset1$Log2FC) < 0){
        ylow = min(comp_subset1$Log2FC)
        yhigh = -min(comp_subset1$Log2FC)
      }
    }

    comp_subset1$ylow = round(ylow,1)
    comp_subset1$yhigh = round(yhigh,1)

    if(length(input$set_comp1) > 1){
      n = length(input$set_comp1)
      if(input$graphics6 == FALSE){
        line_colors = c("blue", "red","green","lightskyblue","purple","hotpink","brown","gold")
      }
      else{
        line_colors = input$ColorChoice1
      }
      qusage_plot = ggplot(data = comp_subset1, aes(x = Index, y = Log2FC, ymin = ylow, ymax = yhigh))
      qusage_plot = qusage_plot + scale_x_discrete(breaks = 1:length(comp_subset1$Index), labels = as.character.factor(comp_subset1$pathway.name[comp_subset1$Index]))
      qusage_plot = qusage_plot + xlab("Modules") + ylab("Pathway Activity")
      qusage_plot = qusage_plot + scale_colour_manual(values = line_colors) + geom_line(aes(colour = Comparison), size = 1) + geom_hline(yintercept = 0, colour = "black", size = .5)
      qusage_plot = qusage_plot + geom_point(aes(shape = SIG), size = 2) + scale_shape_manual(values = c(16,8))
    }

    else{
      qusage_plot = qplot(x = pathway.name, y = Log2FC, data = comp_subset1, ymin = ylow, ymax = yhigh, ylab = "Pathway Activity", xlab = "Modules")
      qusage_plot = qusage_plot + geom_point() + geom_errorbar(data = comp_subset1, aes(x = pathway.name, y = Log2FC, ymin = low, ymax = up))
      qusage_plot = qusage_plot + geom_hline(yintercept = 0, colour = "red", size = 1)
    }

    qusage_plot = qusage_plot + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ coord_cartesian(xlim = ranges$x, ylim = ranges$y)
    return(qusage_plot)

  })


  individual_modplot <- reactive({
    if(is.null(master())){return(NULL)}
    if(is.null(individual_modselect())){return(NULL)}
    mmr_modselection <- individual_modselect()$mmr_modselection
    master <- individual_modselect()$master
    ylow = c()
    yhigh = c()
    if(min(mmr_modselection$low) < 0 & max(mmr_modselection$up) > 0){
      if(abs(min(mmr_modselection$low)) >= max(mmr_modselection$up)){
        ylow = min(mmr_modselection$low)
        yhigh = -min(mmr_modselection$low)
      }
      else{
        ylow = -max(mmr_modselection$up)
        yhigh = max(mmr_modselection$up)
      }
    }
    if(min(mmr_modselection$low) > 0){
      ylow = -max(mmr_modselection$up)
      yhigh = max(mmr_modselection$up)
    }
    if(max(mmr_modselection$up) < 0){
      ylow = min(mmr_modselection$low)
      yhigh = -min(mmr_modselection$low)
    }
    mmr_modselection$ylow = round(ylow, 1)
    mmr_modselection$yhigh = round(yhigh, 1)
    mmr_modselection$Pathway.Activity = master$Log2FC
    mod_plot = ggplot(data = mmr_modselection, aes(x = Index, y = Log2FC, ymin = ylow, ymax = yhigh)) + xlim(0, length(mmr_modselection$Index))
    mod_plot = mod_plot + annotate("rect", xmin = -Inf, xmax = Inf, ymin = master$low, ymax = master$up, alpha = .3, fill = "lightblue")
    mod_plot = mod_plot + scale_x_discrete(limits = as.character(mmr_modselection$SYMBOL[mmr_modselection$Index]))
    mod_plot = mod_plot + geom_point()
    mod_plot = mod_plot + geom_errorbar(data = mmr_modselection, aes(x = Index, y = Log2FC, ymin = low, ymax = up))
    mod_plot = mod_plot + xlab("Symbol ID") + ylab("Probe Level Fold Change") + theme_bw() +
               theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_hline(yintercept = 0, colour = "red", size = 1)
    mod_plot = mod_plot + geom_hline(aes(yintercept = Pathway.Activity), linetype = "dashed",colour = "black") + theme(plot.margin = unit(c(1,10,0,0), "lines"))
    mod_plot = mod_plot + annotation_custom(textGrob(label = "Pathway Activity", hjust = 0), xmin = length(mmr_modselection$Index) + 1.08,
                                            xmax = length(mmr_modselection$Index) + 1.08, ymin = master$Log2FC, ymax = master$Log2FC)
    gt <- ggplot_gtable(ggplot_build(mod_plot))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    gt
  })

  output$GeneSetTab <- renderDataTable({
    if(is.null(master())){return(NULL)}
    if(is.null(individual_modselect())){return(NULL)}
    y <- individual_modselect()$mmr_modselection
    y$Symbol <- paste("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=",y$Symbol," target = '_blank'",'>',y$Symbol,"</a>",sep='')
    y <- y[,-1]
    escape = FALSE
    y
  }, escape = FALSE)

  output$GeneSetPlot <- renderPlot({
    if(is.null(master())){return(NULL)}
    if(is.null(individual_modplot())){return(NULL)}
    grid.draw(individual_modplot())
  },width = plotWidth2)

  output$CompOverview <- renderDataTable({
    if(is.null(master())){return(NULL)}
    comp_overview()
  })

  output$foldchange1 <- renderPlot({
    if(is.null(master())){return(NULL)}
    comp_plot1()
  }, width = plotWidth1)

  output$venn <- renderUI({
    if(is.null(master())){return(NULL)}
    if(is.null(comp_subset1())){return(NULL)}
    if(length(venn.data()) == 1){return(NULL)}
    else{
      plotOutput("venndiagram")
    }
  })

  output$venn.download <- renderUI({
    if(is.null(master())){return(NULL)}
    if(is.null(comp_subset1())){return(NULL)}
    if(length(venn.data()) == 1){return(NULL)}
    else{
      downloadButton("downloadVenn2", "Download Figure")
    }
  })

  venn.d <- reactive({
    if(is.null(master())){return(NULL)}
    color.choices = c("blue", "red","green","lightskyblue","purple","hotpink","brown","gold")
    if(length(venn.data()) == 2){
      color.choices = color.choices[1:2]
    }
    if(length(venn.data()) == 3){
      color.choices = color.choices[1:3]
    }
    if(length(venn.data()) == 4){
      color.choices = color.choices[1:4]
    }
    if(length(venn.data()) == 5){
      color.choices = color.choices[1:5]
    }
    if(input$graphics6 == FALSE){
      return(VennDiagram::venn.diagram(venn.data(), filename = NULL, lwd = 1, col = color.choices, cat.cex = .9, fil = color.choices, margin = .12, ext.text = FALSE,
                          euler.d = TRUE))
    }
    else{
      return(VennDiagram::venn.diagram(venn.data(), filename = NULL, lwd = 1, col = input$ColorChoice1, cat.cex = .9, fil = input$ColorChoice1, margin = .12, ext.text = FALSE,
                          euler.d = TRUE))
    }
  })

  output$venndiagram <- renderPlot({
    if(is.null(comp_subset1())){return(NULL)}
    grid.draw(venn.d())
  }, height = 350, width = 750)

  output$MultipleCompTab <- renderDataTable({
    if(is.null(master())){return(NULL)}
    if(input$only_annotated){
      if(is.null(values$BaylorTF)){
        comp_subset1 <- comp_subset1()[which(comp_subset1()$pathway.name %in% paste(module_annotations$Modulev2_Annotation,module_annotations$Module,sep=" ")),]
        comp_subset1$pathway.name <- as.character(comp_subset1$pathway.name)
        annot <- which(paste(module_annotations$Modulev2_Annotation,module_annotations$Module,sep=" ") %in% comp_subset1()$pathway.name)
      }
      else{
        if(values$BaylorTF == TRUE){
          comp_subset1 <- comp_subset1()[which(comp_subset1()$pathway.name %in% paste(module_annotations$Modulev2_Annotation,module_annotations$Module,sep=" ")),]
        }
        else{
          comp_subset1 <- comp_subset1()[which(comp_subset1()$pathway.name %in% module_annotations$Module),]
        }
      }
    }
    else{
      comp_subset1 <- comp_subset1()
    }
    comp_subset1 <- comp_subset1[,-c(1,11)]
    colnames(comp_subset1)[which(colnames(comp_subset1) == "p.Value")] <- "P.Value"
    comp_subset1[,c(1,9,2,3,4,8,5,6,7)]
  })

  output$downloadSigComps <- downloadHandler(

    filename = function() {paste(values$project_name,'_','Qusage_Sig_Comps_Table','.csv', sep='')  },
    content = function(file) {
      write.csv(comp_overview(), file,row.names = FALSE)
    }
  )

  output$downloadPlot2 <- downloadHandler(
    filename = function() {paste(values$project_name,'_','Multi_Comparisons_Plot','.png', sep = '')},
    content = function(file){
      png(file, width = (plotres1()/72)*plotWidth1(), height = (plotres1()/72)*480, res = plotres1())
      print(comp_plot1())
      dev.off()
    }
  )

  output$downloadVenn2 <- downloadHandler(
    filename = function() {paste(values$project_name,'_','Multi_Comparisons_Venn', '.png', sep = '')},
    content = function(file){
      png(file)
      grid.draw(venn.d())
      dev.off()
    }
  )

  output$downloadTable2 <- downloadHandler(
    filename = function() {paste(values$project_name,'_','Qusage_Data_Table_Multi_Comp','.csv', sep='')  },
    content = function(file) {
      comp_subset1 <- comp_subset1()[,-c(1,11)]
      comp_subset1 <- comp_subset1[,c(1,9,2,3,4,8,5,6,7)]
      write.csv(comp_subset1, file,row.names = FALSE)
    }
  )

  output$downloadPlot3 <- downloadHandler(
    filename = function() {paste(values$project_name,'_', 'Individual_GeneSet_Plot','.png', sep = '')},
    content = function(file){
      png(file, width = (plotres2()/72)*plotWidth2(), height = (plotres2()/72)*480, res = plotres2())
      grid.draw(individual_modplot())
      dev.off()
    }
  )

  output$downloadTable3 <- downloadHandler(
    filename = function() {paste(values$project_name,'_','Qusage_Individual_GeneSet_Data_Table','.csv', sep='')  },
    content = function(file) {
      write.csv(individual_modselect()$mmr_modselection[,-1], file,row.names = FALSE)
    }
  )

  ################### ROAST ############################

  output$roast <- renderMenu({
    if(is.null(values$roast_results)){
      return(strong(""))
    }
    if(is.null(values$roast_results) == FALSE){
      menuItem("Roast", icon = icon("th-list"), tabName = "roast")
    }
  })

  output$setStat <- renderUI({
    selectInput("setStat1", "Select gene set statistic:", choices = names(values$roast_results), selected = names(values$roast_results)[1])
  })

  roast.overview <- reactive({
    results <- values$roast_results[[which(names(values$roast_results) %in% input$setStat1)]]
    Comparison <- names(results)
    Raw <- FDR <- Bonf <- c()
    for(i in 1:length(results)){
      Bonferroni <- p.adjust(results[[i]]$PValue, method = "bonferroni")
      Raw[i] <- length(which(results[[i]]$PValue <= input$SigLevelR))
      FDR[i] <- length(which(results[[i]]$FDR <= input$SigLevelR))
      Bonf[i] <- length(which(Bonferroni <= input$SigLevelR))
    }
    overview <- data.frame(Comparison = Comparison, Raw = Raw, FDR = FDR, Bonf = Bonf)
    overview
  })

  roast.results <- reactive({
    results <- values$roast_results[[which(names(values$roast_results) %in% input$setStat1)]]
    for(i in 2:length(results)){
      results[[i]] <- results[[i]][match(rownames(results[[1]]), rownames(results[[i]]), nomatch = 0),]
    }
    Gene.set <- rownames(results[[1]])
    results <- do.call("cbind", results)
    results <- cbind(Gene.set, results)
    results
  })

  output$CompOverviewR <- renderDataTable({
    roast.overview()
  })

  ###################FLOW PART##########################

  output$flow.data <- renderMenu({
    if(is.null(values$flow_data)){
      #return(menuItem(strong("Flow", style = "color:red"), tabName = ""))
      return(strong(""))
    }
    else{
      return(menuItem("Flow", icon = icon("th-list"), tabName = "flow"))
    }
  })

  FlowMMR<-reactive({
    if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    if(input$FlowVarSet=="Manual") {x<-values$flow_results}
    if(input$FlowVarSet=="Flock") {x<-values$FlockFlowMMR}
    if(input$FlowVarSet=="Spade") {x<-values$SpadeFlowMMR}
    x
  })

  FlowMMR1<-reactive({
    if(input$FlowVarSet1=="No Flow Data"){return(NULL)}
    if(input$FlowVarSet1=="Manual") {x<-values$flow_results}
    if(input$FlowVarSet1=="Flock") {x<-values$FlockFlowMMR}
    if(input$FlowVarSet1=="Spade") {x<-values$SpadeFlowMMR}
    x
  })

  FlowVariableSetNames<-reactive({ind <- which(c("flow_results","FlockFlowMMR","SpadeFlowMMR") %in% names(values))
  if(length(ind)<3){if(length(ind)==0){x<-c("No Flow Data")}
    if(length(ind)>0){x<-c("Manual","Flock","Spade")[ind]}}
  if(length(ind)==3){check1<-is.null(values$flow_results)
  check2<-is.null(values$FlockFlowMMR)
  check3<-is.null(values$SpadeFlowMMR)
  x<-c("Manual","Flock","Spade")[!c(check1,check2,check3)]}
  x
  })

  output$FlowVariableSets<-renderUI({
    selectInput("FlowVarSet", "Select flow variable set (manual, flock, or spade):", FlowVariableSetNames())
  })

  output$FlowVariableSets1<-renderUI({
    selectInput("FlowVarSet1", "Select flow variable set (manual, flock, or spade):", FlowVariableSetNames())
  })

  FClist <- reactive({
    if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    fc_1 <- FlowMMR()[grep("P.Value", names(FlowMMR()))]
    fc_2 <- apply(fc_1, 2, p.adjust, method = "fdr")
    fc_3 <- apply(fc_1, 2, p.adjust, method = "bonferroni")
    sumfc <- function(alpha){
      Raw <- apply(fc_1, 2, function(x) sum(x[!is.na(x)] <= alpha))
      FDR <- apply(fc_2, 2, function(x) sum(x[!is.na(x)] <= alpha))
      Bonf <- apply(fc_3, 2, function(x) sum(x[!is.na(x)] <= alpha))
      sumfc <- rbind(Raw, FDR, Bonf)
      return(sumfc)
    }
    x <- t(sumfc(input$alphaFlow_1))
    y <- data.frame(Comparison = substring(rownames(x), 13), x)
    rownames(y) <- NULL
    y
  })

  output$FlowOverview <- renderDataTable({
    FClist()
  })

  output$downloadFC <- downloadHandler(
    filename = function() {paste0('Flow_Analysis_Overview-Significant_Probes_under_', input$alphaFlow_1, '.csv')  },
    content = function(file) {
      write.csv(FClist(), file,row.names = FALSE)
    }
  )

  output$CompF <- renderUI({
    nam <- names(FlowMMR())
    p.names <- substring(nam[grep("Estimate of ", nam, fixed = TRUE)], 13)
    selectInput("compflow", "Comparison Selection:", p.names, p.names[1])
  })

  flowlist <- reactive({
    if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    dat <- FlowMMR()
    x <- flowlistmaker(dat, input$compflow, input$alphaFlow_2, input$methodF)
    dummy <- data.frame("No Variables Present")
    names(dummy) <- "Flow_Variable"
    if(identical(x,dummy)){return(x)}
    else{
      ind1 <- grep("Estimate", names(x), fixed = TRUE)
      ind2 <- grep("Test.statistic", names(x), fixed = TRUE)
      ind3 <- grep("P.Value", names(x), fixed = TRUE)
      ind4 <- grep("FDR", names(x), fixed = TRUE)
      ind5 <- grep("BONF", names(x), fixed = TRUE)
      x <- x[,c(1,ind1, ind2, ind3, ind4, ind5)]
      names(x) <- c("Flow_Variable","Estimate", "Test.statistic", "P.Value", "FDR", "BONF")
      return(x)
    }
  })

  output$flowlisttable <-renderDataTable({
    if(is.null(flowlist())){return(NULL)}
    flowlist()
  })

  output$downloadFL <- downloadHandler(
    filename = function() {paste0('Flow_Analysis_Lists-', input$compflow, '_under_',  input$alphaFlow_2, '_', input$methodF, '.csv')  },
    content = function(file) {
      write.csv(flowlist(), file, row.names = FALSE)
    }
  )

  FlowData<-reactive({  if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    if (input$FlowVarSet=="Manual") {x<-values$flow_data}
    if (input$FlowVarSet=="Flock") {x<-values$FlockFlowData}
    if (input$FlowVarSet=="Spade") {x<-values$SpadeFlowData}
    subsetindex<-which(x[,which(names(x) %in% values$f_responder_var)] %in% input$FlowSub)
    x[subsetindex,]
  })

  FlowDataForDownload<-reactive({  if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    if (input$FlowVarSet=="Manual") {x<-values$flow_data}
    if (input$FlowVarSet=="Flock") {x<-values$FlockFlowData}
    if (input$FlowVarSet=="Spade") {x<-values$SpadeFlowData}
    subsetindex<-which(x[,which(names(x) %in% values$f_responder_var)] %in% input$FlowSub)
    x[subsetindex,]
    x
  })

  FlowResult<-reactive({
    if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    FlowMMR()
  })

  output$FlowResultTable<-renderDataTable({
    if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    head(FlowResult())
  })

  output$FlowDataNames<-renderUI({
    if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    selectInput("FlowPlotVars","Variable(s):",as.character(FlowMMR()$Flow.variable), multiple = TRUE, selected = as.character(FlowMMR()$Flow.variable)[1])
  })

  output$FlowSummary<-renderTable({
    if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    if (input$FlowVarSet=="Manual") {x<-values$flow_results}
    if (input$FlowVarSet=="Flock") {x<-values$FlockFlowMMR}
    if (input$FlowVarSet=="Spade") {x<-values$SpadeFlowMMR}
    x[is.na(x)]<-1
    head(x)
  }, include.rownames = FALSE)

  output$flowmax<-renderUI({ numericInput("FlowMax","Max value to be plotted for time:",max(FlowData()[,which(names(FlowData()) %in% values$f_time_var)]))  })
  output$flowmin<-renderUI({ numericInput("FlowMin","Min value to be plotted for time:",min(FlowData()[,which(names(FlowData()) %in% values$f_time_var)]))  })
  output$flowsub<-renderUI({ if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    if (input$FlowVarSet=="Manual") {x<-values$flow_data}
    if (input$FlowVarSet=="Flock") {x<-values$FlockFlowData}
    if (input$FlowVarSet=="Spade") {x<-values$SpadeFlowData}
    y<-unique(as.character(x[,which(names(x) %in% values$f_responder_var)]))
    selectInput("FlowSub","Responder levels:",y,y[1:length(y)],multiple=TRUE)})

  flowplot <- reactive({
    if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    if(length(input$FlowPlotVars) == 0){return(NULL)}
    if(length(input$FlowSub) == 0){return(NULL)}
    if(is.numeric(input$FlowMin) == FALSE & is.numeric(input$FlowMax) == FALSE){return(NULL)}

    if(length(input$FlowPlotVars) == 1){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])

      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
    }

    if(length(input$FlowPlotVars) == 2){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      mydata2 <- data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])

      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
        mydata2 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

    }

    if(length(input$FlowPlotVars) == 3){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      mydata2 <- data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      mydata3 <- data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])

      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
        mydata2 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
        mydata3 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
    }

    if(length(input$FlowPlotVars) >= 4){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      mydata2 <- data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      mydata3 <- data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      mydata4 <- data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[4])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])

      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
        mydata2 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
        mydata3 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
        mydata4 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[4])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

    }

    mysummary<-function(x){
      y<-c(mean(x),sd(x),length(x))
      names(y)<-c("Mean","Sd","N")
      return(y)
    }

    if(length(input$FlowPlotVars) == 1){
      sum1<-aggregate(FlowVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))
      myplot<-ggplot(data=result,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +
        scale_size(range=c(.1, 2),guide=FALSE)+ xlim(input$FlowMin, input$FlowMax)+ggtitle(input$FlowPlotVars[1])
      if(input$FlowSamples){
        myplot<-myplot+geom_point(data=mydata,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=.25))
      }
    }
    if(length(input$FlowPlotVars) == 2){
      sum1<-aggregate(FlowVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))
      myplot<-ggplot(data=result,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
              geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) + scale_size(range=c(.1, 2),guide=FALSE)+
              xlim(input$FlowMin, input$FlowMax)+ggtitle(input$FlowPlotVars[1])
      sum2 <- aggregate(FlowVar~responder+time, data=mydata2, mysummary)
      summaries2 <-sum2[,3]
      result2 <-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))
      myplot2 <-ggplot(data=result2,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
                geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
                xlim(input$FlowMin, input$FlowMax)+ggtitle(input$FlowPlotVars[2])
      if(input$FlowSamples){
        myplot<-myplot+geom_point(data=mydata,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=1))+
                geom_line(data=mydata,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=.25))

        myplot2 <-myplot2+geom_point(data=mydata2,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=1))+
                  geom_line(data=mydata2,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=.25))
      }
    }

    if(length(input$FlowPlotVars) == 3){
      sum1<-aggregate(FlowVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))
      myplot<-ggplot(data=result,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
              geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
              xlim(input$FlowMin, input$FlowMax)+ggtitle(input$FlowPlotVars[1])
      sum2 <- aggregate(FlowVar~responder+time, data=mydata2, mysummary)
      summaries2 <-sum2[,3]
      result2 <-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))
      myplot2 <-ggplot(data=result2,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
                geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
                xlim(input$FlowMin, input$FlowMax)+ggtitle(input$FlowPlotVars[2])
      sum3 <- aggregate(FlowVar~responder+time, data=mydata3, mysummary)
      summaries3 <-sum3[,3]
      result3 <-as.data.frame(cbind(sum3[,c(1,2)],summaries3))
      names(result3)<-c("ResponderStatus","Time",colnames(summaries3))
      myplot3 <-ggplot(data=result3,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +
        scale_size(range=c(.1, 2),guide=FALSE)+ xlim(input$FlowMin, input$FlowMax)+ggtitle(input$FlowPlotVars[3])
      if(input$FlowSamples){
        myplot<-myplot+geom_point(data=mydata,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=1))+
                geom_line(data=mydata,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=.25))
        myplot2 <-myplot2+geom_point(data=mydata2,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=1))+
                  geom_line(data=mydata2,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=.25))
        myplot3 <-myplot3+geom_point(data=mydata3,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=1))+
                  geom_line(data=mydata3,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=.25))
      }
    }

    if(length(input$FlowPlotVars) >= 4){
      sum1<-aggregate(FlowVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))
      myplot<-ggplot(data=result,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
              geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7))+scale_size(range=c(.1, 2),guide=FALSE)+
              xlim(input$FlowMin, input$FlowMax)+ggtitle(input$FlowPlotVars[1])
      sum2 <- aggregate(FlowVar~responder+time, data=mydata2, mysummary)
      summaries2 <-sum2[,3]
      result2 <-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))
      myplot2 <-ggplot(data=result2,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
                geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
                xlim(input$FlowMin, input$FlowMax)+ggtitle(input$FlowPlotVars[2])
      sum3 <- aggregate(FlowVar~responder+time, data=mydata3, mysummary)
      summaries3 <-sum3[,3]
      result3 <-as.data.frame(cbind(sum3[,c(1,2)],summaries3))
      names(result3)<-c("ResponderStatus","Time",colnames(summaries3))
      myplot3 <-ggplot(data=result3,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
                geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
                xlim(input$FlowMin, input$FlowMax)+ggtitle(input$FlowPlotVars[3])
      sum4 <- aggregate(FlowVar~responder+time, data=mydata4, mysummary)
      summaries4 <-sum4[,3]
      result4 <-as.data.frame(cbind(sum4[,c(1,2)],summaries4))
      names(result4)<-c("ResponderStatus","Time",colnames(summaries4))
      myplot4 <-ggplot(data=result4,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
                geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
                xlim(input$FlowMin, input$FlowMax)+ggtitle(input$FlowPlotVars[4])

      if(input$FlowSamples){
        myplot <- myplot+geom_point(data=mydata,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=1))+
                  geom_line(data=mydata,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=.25))
        myplot2 <- myplot2+geom_point(data=mydata2,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=1))+
                   geom_line(data=mydata2,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=.25))
        myplot3 <- myplot3+geom_point(data=mydata3,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=1))+
                   geom_line(data=mydata3,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=.25))
        myplot4 <- myplot4+geom_point(data=mydata4,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=1))+
                   geom_line(data=mydata4,aes(x=time,y=FlowVar,group=subjectid,colour=responder,size=.25))
      }
    }

    if(length(input$FlowPlotVars) == 1){
      z = myplot
    }
    if(length(input$FlowPlotVars) == 2){
      z = multiplot(myplot, myplot2, cols = 2)
    }
    if(length(input$FlowPlotVars) == 3){
      z = multiplot(myplot, myplot2, myplot3, cols = 2)
    }
    if(length(input$FlowPlotVars) >= 4){
      z = multiplot(myplot, myplot2, myplot3, myplot4, cols = 2)
    }
    z
  })

  flowplot2 <- reactive({
    if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    if(length(input$FlowPlotVars) == 0){return(NULL)}
    if(length(input$FlowPlotVars) == 1){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      mydata<-mydata[order(mydata$time,mydata$responder),]
      mydata$time<-factor(mydata$time,ordered=T)
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),ordered=T)
      dummy<-data.frame(do.call(rbind,strsplit(levels(mydata$Time_Responder),"_")))
      index<-order(as.numeric(as.character(dummy$X1)))
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),levels(mydata$Time_Responder)[index])
      if(input$flowbox){
        myplot<-ggplot(data=mydata,aes(x=time,y=FlowVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$FlowPlotVars[1])
      }
      if(!input$flowbox){
        myplot<-ggplot(data=mydata,aes(x=Time_Responder,y=FlowVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$FlowPlotVars[1])
      }
    }
    if(length(input$FlowPlotVars) == 2){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      mydata<-mydata[order(mydata$time,mydata$responder),]
      mydata$time<-factor(mydata$time,ordered=T)
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),ordered=T)
      dummy<-data.frame(do.call(rbind,strsplit(levels(mydata$Time_Responder),"_")))
      index<-order(as.numeric(as.character(dummy$X1)))
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),levels(mydata$Time_Responder)[index])
      mydata2 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata2 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      mydata2 <-mydata2[order(mydata2$time,mydata2$responder),]
      mydata2$time<-factor(mydata2$time,ordered=T)
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),ordered=T)
      dummy2 <-data.frame(do.call(rbind,strsplit(levels(mydata2$Time_Responder),"_")))
      index2 <-order(as.numeric(as.character(dummy2$X1)))
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),levels(mydata2$Time_Responder)[index2])
      if(input$flowbox){
        myplot<-ggplot(data=mydata,aes(x=time,y=FlowVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$FlowPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=time,y=FlowVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$FlowPlotVars[2])
      }
      if(!input$flowbox){
        myplot<-ggplot(data=mydata,aes(x=Time_Responder,y=FlowVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$FlowPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=Time_Responder,y=FlowVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$FlowPlotVars[2])
      }
    }
    if(length(input$FlowPlotVars) == 3){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      mydata<-mydata[order(mydata$time,mydata$responder),]
      mydata$time<-factor(mydata$time,ordered=T)
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),ordered=T)
      dummy<-data.frame(do.call(rbind,strsplit(levels(mydata$Time_Responder),"_")))
      index<-order(as.numeric(as.character(dummy$X1)))
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),levels(mydata$Time_Responder)[index])
      mydata2 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata2 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      mydata2 <-mydata2[order(mydata2$time,mydata2$responder),]
      mydata2$time<-factor(mydata2$time,ordered=T)
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),ordered=T)
      dummy2 <-data.frame(do.call(rbind,strsplit(levels(mydata2$Time_Responder),"_")))
      index2 <-order(as.numeric(as.character(dummy2$X1)))
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),levels(mydata2$Time_Responder)[index2])
      mydata3 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata3 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      mydata3 <-mydata3[order(mydata3$time,mydata3$responder),]
      mydata3$time<-factor(mydata3$time,ordered=T)
      mydata3$Time_Responder<-factor(paste(mydata3$time,mydata3$responder,sep="_"),ordered=T)
      dummy3 <-data.frame(do.call(rbind,strsplit(levels(mydata3$Time_Responder),"_")))
      index3 <-order(as.numeric(as.character(dummy3$X1)))
      mydata3$Time_Responder<-factor(paste(mydata3$time,mydata3$responder,sep="_"),levels(mydata3$Time_Responder)[index3])
      if(input$flowbox){
        myplot<-ggplot(data=mydata,aes(x=time,y=FlowVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$FlowPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=time,y=FlowVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$FlowPlotVars[2])
        myplot3 <-ggplot(data=mydata3,aes(x=time,y=FlowVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$FlowPlotVars[3])
      }
      if(!input$flowbox){
        myplot<-ggplot(data=mydata,aes(x=Time_Responder,y=FlowVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$FlowPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=Time_Responder,y=FlowVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$FlowPlotVars[2])
        myplot3 <-ggplot(data=mydata3,aes(x=Time_Responder,y=FlowVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$FlowPlotVars[3])
      }
    }
    if(length(input$FlowPlotVars) >= 4){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      mydata<-mydata[order(mydata$time,mydata$responder),]
      mydata$time<-factor(mydata$time,ordered=T)
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),ordered=T)
      dummy<-data.frame(do.call(rbind,strsplit(levels(mydata$Time_Responder),"_")))
      index<-order(as.numeric(as.character(dummy$X1)))
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),levels(mydata$Time_Responder)[index])

      mydata2 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata2 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      mydata2 <-mydata2[order(mydata2$time,mydata2$responder),]
      mydata2$time<-factor(mydata2$time,ordered=T)
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),ordered=T)
      dummy2 <-data.frame(do.call(rbind,strsplit(levels(mydata2$Time_Responder),"_")))
      index2 <-order(as.numeric(as.character(dummy2$X1)))
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),levels(mydata2$Time_Responder)[index2])

      mydata3 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata3 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      mydata3 <-mydata3[order(mydata3$time,mydata3$responder),]
      mydata3$time<-factor(mydata3$time,ordered=T)
      mydata3$Time_Responder<-factor(paste(mydata3$time,mydata3$responder,sep="_"),ordered=T)
      dummy3 <-data.frame(do.call(rbind,strsplit(levels(mydata3$Time_Responder),"_")))
      index3 <-order(as.numeric(as.character(dummy3$X1)))
      mydata3$Time_Responder<-factor(paste(mydata3$time,mydata3$responder,sep="_"),levels(mydata3$Time_Responder)[index3])

      mydata4 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[4])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata4 <-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[4])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      mydata4 <-mydata4[order(mydata4$time,mydata4$responder),]
      mydata4$time<-factor(mydata4$time,ordered=T)
      mydata4$Time_Responder<-factor(paste(mydata4$time,mydata4$responder,sep="_"),ordered=T)
      dummy4 <-data.frame(do.call(rbind,strsplit(levels(mydata4$Time_Responder),"_")))
      index4 <-order(as.numeric(as.character(dummy4$X1)))
      mydata4$Time_Responder<-factor(paste(mydata4$time,mydata4$responder,sep="_"),levels(mydata4$Time_Responder)[index4])

      if(input$flowbox){
        myplot<-ggplot(data=mydata,aes(x=time,y=FlowVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$FlowPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=time,y=FlowVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$FlowPlotVars[2])
        myplot3 <-ggplot(data=mydata3,aes(x=time,y=FlowVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$FlowPlotVars[3])
        myplot4 <-ggplot(data=mydata4,aes(x=time,y=FlowVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$FlowPlotVars[4])
      }

      if(!input$flowbox){
        myplot<-ggplot(data=mydata,aes(x=Time_Responder,y=FlowVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$FlowPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=Time_Responder,y=FlowVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$FlowPlotVars[2])
        myplot3 <-ggplot(data=mydata3,aes(x=Time_Responder,y=FlowVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$FlowPlotVars[3])
        myplot4 <-ggplot(data=mydata4,aes(x=Time_Responder,y=FlowVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$FlowPlotVars[4])
      }
    }
    if(length(input$FlowPlotVars) == 1){
      z = myplot
    }
    if(length(input$FlowPlotVars) == 2){
      z = multiplot(myplot, myplot2, cols = 2)
    }
    if(length(input$FlowPlotVars) == 3){
      z = multiplot(myplot, myplot2, myplot3, cols = 2)
    }
    if(length(input$FlowPlotVars) >= 4){
      z = multiplot(myplot, myplot2, myplot3, myplot4, cols = 2)
    }
    z
  })

  output$FlowPlot<-renderPlot({
    print(flowplot())
  })

  output$FlowPlot2<-renderPlot({
    print(flowplot2())
  })

  FlowPlotSummary1<-reactive({
    if(input$FlowVarSet=="No Flow Data"){return(NULL)}
    if(length(input$FlowPlotVars) == 0){return(NULL)}
    if(length(input$FlowSub) == 0){return(NULL)}
    mysummary<-function(x){
      y<-c(mean(x),sd(x),length(x))
      names(y)<-c("Mean","Sd","N")
      return(y)
    }
    if(length(input$FlowPlotVars) == 1){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }
      sum1<-aggregate(FlowVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))
    }

    if(length(input$FlowPlotVars) == 2){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

      sum1<-aggregate(FlowVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))

      mydata2 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata2<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

      sum2<-aggregate(FlowVar~responder+time, data=mydata2, mysummary)
      summaries2<-sum2[,3]
      result2<-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))
    }

    if(length(input$FlowPlotVars) == 3){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

      sum1<-aggregate(FlowVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))

      mydata2 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata2<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

      sum2<-aggregate(FlowVar~responder+time, data=mydata2, mysummary)
      summaries2<-sum2[,3]
      result2<-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))

      mydata3 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata3<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

      sum3<-aggregate(FlowVar~responder+time, data=mydata3, mysummary)
      summaries3<-sum3[,3]
      result3<-as.data.frame(cbind(sum3[,c(1,2)],summaries3))
      names(result3)<-c("ResponderStatus","Time",colnames(summaries3))
    }

    if(length(input$FlowPlotVars) >= 4){
      mydata<-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[1])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

      sum1<-aggregate(FlowVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))

      mydata2 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata2<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[2])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

      sum2<-aggregate(FlowVar~responder+time, data=mydata2, mysummary)
      summaries2<-sum2[,3]
      result2<-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))

      mydata3 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata3<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[3])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

      sum3<-aggregate(FlowVar~responder+time, data=mydata3, mysummary)
      summaries3<-sum3[,3]
      result3<-as.data.frame(cbind(sum3[,c(1,2)],summaries3))
      names(result3)<-c("ResponderStatus","Time",colnames(summaries3))

      mydata4 <-data.frame(FlowVar=FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[4])],time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      if(input$FlowTransform){
        mydata4<-data.frame(FlowVar=log2(FlowData()[,which(names(FlowData()) %in% input$FlowPlotVars[4])]),time=FlowData()[,which(names(FlowData()) %in% values$f_time_var)],subjectid=FlowData()[,which(names(FlowData()) %in% values$f_patient_id)],responder=FlowData()[,which(names(FlowData()) %in% values$f_responder_var)])
      }

      sum4<-aggregate(FlowVar~responder+time, data=mydata4, mysummary)
      summaries4<-sum4[,3]
      result4<-as.data.frame(cbind(sum4[,c(1,2)],summaries4))
      names(result4)<-c("ResponderStatus","Time",colnames(summaries4))
    }

    if(length(input$FlowPlotVars) == 1){
      z <- result
    }

    if(length(input$FlowPlotVars) == 2){
      FlowVariable <- c(rep(input$FlowPlotVars[1], length(result[,1])), rep(input$FlowPlotVars[2], length(result2[,1])))
      z <- rbind(result, result2)
      z <- cbind(FlowVariable, z)
    }

    if(length(input$FlowPlotVars) == 3){
      FlowVariable <- c(rep(input$FlowPlotVars[1], length(result[,1])), rep(input$FlowPlotVars[2], length(result2[,1])), rep(input$FlowPlotVars[3], length(result3[,1])))
      z <- rbind(result, result2, result3)
      z <- cbind(FlowVariable, z)
    }

    if(length(input$FlowPlotVars) >= 4){
      FlowVariable <- c(rep(input$FlowPlotVars[1], length(result[,1])), rep(input$FlowPlotVars[2], length(result2[,1])), rep(input$FlowPlotVars[3], length(result3[,1])), rep(input$FlowPlotVars[4], length(result4[,1])))
      z <- rbind(result, result2, result3, result4)
      z <- cbind(FlowVariable, z)
    }

    z

  })

  output$FlowPlotSummary<-renderDataTable({
    FlowPlotSummary1()
  })

  output$downloadFlowResults <- downloadHandler(
    filename = function() {paste(values$project_name,'_',input$comparisonflow,'.csv', sep='')  },
    content = function(file) {
      write.csv(FlowResult(), file,row.names = FALSE)
    }
  )

  output$downloadFlowData <- downloadHandler(
    filename = function() {paste(values$project_name,'_FlowData','.csv', sep='')  },
    content = function(file) {
      write.csv(FlowDataForDownload(), file, row.names = FALSE)
    }
  )

  output$downloadFlowSummaries <- downloadHandler(
    filename = function() {paste(values$project_name,'_',"FlowSummary_",input$FlowPlotVars[1],'.csv', sep='')  },
    content = function(file) {
      write.csv(FlowPlotSummary1(), file, row.names = FALSE)
    }
  )

  output$downloadFlowPlot <- downloadHandler(
    filename = function() {paste(values$project_name,'_','FlowPlot1','.png', sep = '')},
    content = function(file){
      png(file, width = 900)
      print(flowplot())
      dev.off()
    }
  )

  output$downloadFlowPlot2 <- downloadHandler(
    filename = function() {paste(values$project_name,'_','FlowPlot2','.png', sep = '')},
    content = function(file){
      png(file, width = 900)
      print(flowplot2())
      dev.off()
    }
  )


  ###############END OF FLOW PART####################

  ############### METABOLOMICS PART #################

  output$metab.data <- renderMenu({
    if(is.null(values$metab_data)){
      #return(menuItem(strong("Metab", style = "color:red"), tabName = ""))
      return(strong(""))
    }
    else{
      return(menuItem("Metab", icon = icon("th-list"), tabName = "metab2",
                      menuSubItem("PCA", tabName = "pcaM"),
                      menuSubItem("Differential Analysis", tabName = "metab")))
    }
  })

  MetabMMR<-reactive({
    if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    if(input$MetabVarSet=="Manual") {x<-values$metab_results}
    if(input$MetabVarSet=="Flock") {x<-values$FlockMetabMMR}
    if(input$MetabVarSet=="Spade") {x<-values$SpadeMetabMMR}
    x
  })

  MetabMMR1<-reactive({
    if(input$MetabVarSet1=="No Metab Data"){return(NULL)}
    if(input$MetabVarSet1=="Manual") {x<-values$metab_results}
    if(input$MetabVarSet1=="Flock") {x<-values$FlockMetabMMR}
    if(input$MetabVarSet1=="Spade") {x<-values$SpadeMetabMMR}
    x
  })

  MetabVariableSetNames<-reactive({ind <- which(c("metab_results","FlockMetabMMR","SpadeMetabMMR") %in% names(values))
  if(length(ind)<3){if(length(ind)==0){x<-c("No Metab Data")}
    if(length(ind)>0){x<-c("Manual","Flock","Spade")[ind]}}
  if(length(ind)==3){check1<-is.null(values$metab_results)
  check2<-is.null(values$FlockMetabMMR)
  check3<-is.null(values$SpadeMetabMMR)
  x<-c("Manual","Flock","Spade")[!c(check1,check2,check3)]}
  x
  })

  output$MetabVariableSets<-renderUI({
    selectInput("MetabVarSet", "Select metab variable set (manual, flock, or spade):", MetabVariableSetNames())
  })

  output$MetabVariableSets1<-renderUI({
    selectInput("MetabVarSet1", "Select metab variable set (manual, flock, or spade):", MetabVariableSetNames())
  })

  Mlist <- reactive({
    if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    fc_1 <- MetabMMR()[grep("P.Value", names(MetabMMR()))]
    fc_2 <- apply(fc_1, 2, p.adjust, method = "fdr")
    fc_3 <- apply(fc_1, 2, p.adjust, method = "bonferroni")
    sumfc <- function(alpha){
      Raw <- apply(fc_1, 2, function(x) sum(x[!is.na(x)] <= alpha))
      FDR <- apply(fc_2, 2, function(x) sum(x[!is.na(x)] <= alpha))
      Bonf <- apply(fc_3, 2, function(x) sum(x[!is.na(x)] <= alpha))
      sumfc <- rbind(Raw, FDR, Bonf)
      return(sumfc)
    }
    x <- t(sumfc(input$alphaMetab_1))
    y <- data.frame(Comparison = substring(rownames(x), 13), x)
    rownames(y) <- NULL
    y
  })

  output$MetabOverview <- renderDataTable({
    Mlist()
  })

  output$downloadFC2 <- downloadHandler(
    filename = function() {paste0('Metab_Analysis_Overview-Significant_Probes_under_', input$alphaMetab_1, '.csv')  },
    content = function(file) {
      write.csv(Mlist(), file,row.names = FALSE)
    }
  )

  output$CompM <- renderUI({
    nam <- names(MetabMMR())
    p.names <- substring(nam[grep("Estimate of ", nam, fixed = TRUE)], 13)
    selectInput("compmetab", "Comparison Selection:", p.names, p.names[1])
  })

  metablist <- reactive({
    if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    dat <- MetabMMR()
    x <- flowlistmaker(dat, input$compmetab, input$alphaMetab_2, input$methodM)
    dummy <- data.frame("No Variables Present")
    names(dummy) <- "Metab_Variable"
    if(identical(x,dummy)){return(x)}
    else{
      ind1 <- grep("Estimate", names(x), fixed = TRUE)
      ind2 <- grep("Test.statistic", names(x), fixed = TRUE)
      ind3 <- grep("P.Value", names(x), fixed = TRUE)
      ind4 <- grep("FDR", names(x), fixed = TRUE)
      ind5 <- grep("BONF", names(x), fixed = TRUE)
      x <- x[,c(1,ind1, ind2, ind3, ind4, ind5)]
      names(x) <- c("Metab_Variable","Estimate", "Test.statistic", "P.Value", "FDR", "BONF")
      return(x)
    }
  })

  output$metablisttable <-renderDataTable({
    if(is.null(metablist())){return(NULL)}
    metablist()
  })

  output$downloadME <- downloadHandler(
    filename = function() {paste0('Metab_Analysis_Lists-', input$compmetab, '_under_',  input$alphaMetab_2, '_', input$methodM, '.csv')  },
    content = function(file) {
      write.csv(metablist(), file, row.names = FALSE)
    }
  )

  MetabData<-reactive({  if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    if (input$MetabVarSet=="Manual") {x<-values$metab_data}
    if (input$MetabVarSet=="Flock") {x<-values$FlockMetabData}
    if (input$MetabVarSet=="Spade") {x<-values$SpadeMetabData}
    subsetindex<-which(x[,which(names(x) %in% values$m_responder_var)] %in% input$MetabSub)
    x[subsetindex,]
  })

  MetabDataForDownload<-reactive({  if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    if (input$MetabVarSet=="Manual") {x<-values$metab_data}
    if (input$MetabVarSet=="Flock") {x<-values$FlockMetabData}
    if (input$MetabVarSet=="Spade") {x<-values$SpadeMetabData}
    subsetindex<-which(x[,which(names(x) %in% values$m_responder_var)] %in% input$MetabSub)
    x[subsetindex,]
    x
  })

  MetabResult<-reactive({
    if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    MetabMMR()
  })

  output$MetabResultTable<-renderDataTable({
    if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    head(MetabResult())
  })

  output$MetabDataNames<-renderUI({
    if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    selectInput("MetabPlotVars","Variable(s):",as.character(MetabMMR()$Metab.variable), multiple = TRUE, selected = as.character(MetabMMR()$Metab.variable)[1])
  })

  output$pcaAnnotM <- renderUI({
    if(input$MetabVarSet1=="No Metab Data"){return(NULL)}
    dat <- values$metab_data
    num <- sapply(dat, is.numeric)
    groups <- dat[,which(num == FALSE),drop = FALSE]
    selectInput("pcaAnnot2M", "Select grouping factor:", colnames(groups), selected = colnames(groups)[which(colnames(groups) %in% values$m_responder_var)])
  })

  output$pcaAnnotValsM <- renderUI({
    if(input$MetabVarSet1=="No Metab Data"){return(NULL)}
    dat <- values$metab_data
    num <- sapply(dat, is.numeric)
    groups <- dat[,which(num == FALSE),drop = FALSE]
    values <- unique(groups[,input$pcaAnnot2M])
    selectInput("pcaAnnValsM", "Select grouping factor values:", values, selected = values, multiple = TRUE)
  })

  output$pcaBlockingM <- renderUI({
    if(input$MetabVarSet1=="No Metab Data" || is.null(values$m_time_var)){return(NULL)}
    dat <- values$metab_data
    selectInput("pcaBlockM", "PCA by time:", c("All timepoints", unique(dat[,values$m_time_var])), selected = "All timepoints")
  })

  PCAM <- reactive({
    dat <- values$metab_data
    if(!is.null(values$m_time_var)){
      if(input$pcaBlockM == "All timepoints"){
        metabDat <- dat[which(dat[,input$pcaAnnot2M] %in% input$pcaAnnValsM),which(colnames(dat) %in% as.character(MetabMMR1()$Metab.variable))]
        dat <- dat[which(dat[,input$pcaAnnot2M] %in% input$pcaAnnValsM),]
      }
      else{
        metabDat <- dat[which(dat[,input$pcaAnnot2M] %in% input$pcaAnnValsM & dat[,values$m_time_var] %in% input$pcaBlockM),which(colnames(dat) %in% as.character(MetabMMR1()$Metab.variable))]
        dat <- dat[which(dat[,input$pcaAnnot2M] %in% input$pcaAnnValsM & dat[,values$m_time_var] %in% input$pcaBlockM),]
      }
    }
    else{
      metabDat <- dat[which(dat[,input$pcaAnnot2M] %in% input$pcaAnnValsM),which(colnames(dat) %in% as.character(MetabMMR1()$Metab.variable))]
      dat <- dat[which(dat[,input$pcaAnnot2M] %in% input$pcaAnnValsM),]
    }
    sds <- sapply(metabDat, sd)
    if(length(which(sds == 0)) > 0){
      metabDat <- metabDat[,-which(sds == 0)]
    }
    pca <- prcomp(metabDat, scale = TRUE)
    z <- list(dat = dat, pca = pca)
  })

  PCAdatM <- reactive({
    pca <- PCAM()$pca
    dat <- PCAM()$dat
    num <- sapply(dat, is.numeric)
    groups = dat[,which(num == FALSE),drop = FALSE]
    annotations <- groups[,input$pcaAnnot2M]
    pcs <- data.frame(pca$x, annotations)
    colnames(pcs)[ncol(pcs)] <- "Group"
    pcs
  })

  plotw.pcaM <- function(){
    input$PlotWidthPCAM
  }

  ploth.pcaM <- function(){
    input$PlotHeightPCAM
  }

  observeEvent(input$showM, {
    shinyjs::toggle("screePlot2M")
  })

  output$screePlotM <- renderUI({
    plotOutput("screePlot2M")
  })

  output$PCAplot3dM <- renderUI({
    if(input$PCSnumM != 3){
      return(NULL)
    }
    else{
      return(plotOutput("pcaPlot3dM"))
    }
  })

  output$screePlot2M <- renderPlot({
    dat <- as.data.frame(summary(PCAM()$pca)[[length(summary(PCAM()$pca))]])
    dat <- as.data.frame(t(dat))
    dat <- cbind(PC = rownames(dat), dat)
    dat[,2] <- NULL
    dat <- reshape2::melt(dat, id = "PC")
    dat$PC <- factor(dat$PC, levels = dat$PC)
    p <- ggplot(data = dat[which(dat$PC %in% levels(dat$PC)[1:10]),], aes(x = PC, y = value, group = variable)) + geom_line(size = 1.5, aes(colour = variable))
    p <- p + geom_point(size = 3, shape = 21, fill = "white") + ylim(0,1) + ylab("Proportion of Variance") + theme_minimal()
    p
  })


  output$PCAplotM <- renderPlot({
    qplot(PC1, PC2, data = PCAdatM(), colour = Group, size = I(input$CircleSizePCAM))
  }, width = plotw.pcaM, height = ploth.pcaM)

  output$pcaPlot3dM <- renderPlot({
    pca3d::pca3d(PCAM()$pca, group = PCAdatM()$Group, bg = "black", legend = "topright")
  })

  output$MetabSummary<-renderTable({
    if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    if (input$MetabVarSet=="Manual") {x<-values$metab_results}
    if (input$MetabVarSet=="Flock") {x<-values$FlockMetabMMR}
    if (input$MetabVarSet=="Spade") {x<-values$SpadeMetabMMR}
    x[is.na(x)]<-1
    head(x)
  }, include.rownames = FALSE)

  output$metabmax<-renderUI({ numericInput("MetabMax","Max value to be plotted for time:",max(MetabData()[,which(names(MetabData()) %in% values$m_time_var)]))  })
  output$metabmin<-renderUI({ numericInput("MetabMin","Min value to be plotted for time:",min(MetabData()[,which(names(MetabData()) %in% values$m_time_var)]))  })
  output$metabsub<-renderUI({ if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    if (input$MetabVarSet=="Manual") {x<-values$metab_data}
    if (input$MetabVarSet=="Flock") {x<-values$FlockMetabData}
    if (input$MetabVarSet=="Spade") {x<-values$SpadeMetabData}
    y<-unique(as.character(x[,which(names(x) %in% values$m_responder_var)]))
    selectInput("MetabSub","Responder levels:",y,y[1:length(y)],multiple=TRUE)})

  metabplot <- reactive({
    if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    if(length(input$MetabPlotVars) == 0){return(NULL)}
    if(length(input$MetabSub) == 0){return(NULL)}
    if(is.numeric(input$MetabMin) == FALSE & is.numeric(input$MetabMax) == FALSE){return(NULL)}

    if(length(input$MetabPlotVars) == 1){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])

      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
    }

    if(length(input$MetabPlotVars) == 2){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      mydata2 <- data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])

      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
        mydata2 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

    }

    if(length(input$MetabPlotVars) == 3){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      mydata2 <- data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      mydata3 <- data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])

      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
        mydata2 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
        mydata3 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
    }

    if(length(input$MetabPlotVars) >= 4){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      mydata2 <- data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      mydata3 <- data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      mydata4 <- data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[4])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])

      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
        mydata2 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
        mydata3 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
        mydata4 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[4])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

    }

    mysummary<-function(x){
      y<-c(mean(x),sd(x),length(x))
      names(y)<-c("Mean","Sd","N")
      return(y)
    }

    if(length(input$MetabPlotVars) == 1){
      sum1<-aggregate(MetabVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))
      myplot<-ggplot(data=result,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +
        scale_size(range=c(.1, 2),guide=FALSE)+ xlim(input$MetabMin, input$MetabMax)+ggtitle(input$MetabPlotVars[1])
      if(input$MetabSamples){
        myplot<-myplot+geom_point(data=mydata,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=.25))
      }
    }
    if(length(input$MetabPlotVars) == 2){
      sum1<-aggregate(MetabVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))
      myplot<-ggplot(data=result,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) + scale_size(range=c(.1, 2),guide=FALSE)+
        xlim(input$MetabMin, input$MetabMax)+ggtitle(input$MetabPlotVars[1])
      sum2 <- aggregate(MetabVar~responder+time, data=mydata2, mysummary)
      summaries2 <-sum2[,3]
      result2 <-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))
      myplot2 <-ggplot(data=result2,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
        xlim(input$MetabMin, input$MetabMax)+ggtitle(input$MetabPlotVars[2])
      if(input$MetabSamples){
        myplot<-myplot+geom_point(data=mydata,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=.25))

        myplot2 <-myplot2+geom_point(data=mydata2,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata2,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=.25))
      }
    }

    if(length(input$MetabPlotVars) == 3){
      sum1<-aggregate(MetabVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))
      myplot<-ggplot(data=result,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
        xlim(input$MetabMin, input$MetabMax)+ggtitle(input$MetabPlotVars[1])
      sum2 <- aggregate(MetabVar~responder+time, data=mydata2, mysummary)
      summaries2 <-sum2[,3]
      result2 <-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))
      myplot2 <-ggplot(data=result2,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
        xlim(input$MetabMin, input$MetabMax)+ggtitle(input$MetabPlotVars[2])
      sum3 <- aggregate(MetabVar~responder+time, data=mydata3, mysummary)
      summaries3 <-sum3[,3]
      result3 <-as.data.frame(cbind(sum3[,c(1,2)],summaries3))
      names(result3)<-c("ResponderStatus","Time",colnames(summaries3))
      myplot3 <-ggplot(data=result3,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +
        scale_size(range=c(.1, 2),guide=FALSE)+ xlim(input$MetabMin, input$MetabMax)+ggtitle(input$MetabPlotVars[3])
      if(input$MetabSamples){
        myplot<-myplot+geom_point(data=mydata,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=.25))
        myplot2 <-myplot2+geom_point(data=mydata2,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata2,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=.25))
        myplot3 <-myplot3+geom_point(data=mydata3,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata3,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=.25))
      }
    }

    if(length(input$MetabPlotVars) >= 4){
      sum1<-aggregate(MetabVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))
      myplot<-ggplot(data=result,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7))+scale_size(range=c(.1, 2),guide=FALSE)+
        xlim(input$MetabMin, input$MetabMax)+ggtitle(input$MetabPlotVars[1])
      sum2 <- aggregate(MetabVar~responder+time, data=mydata2, mysummary)
      summaries2 <-sum2[,3]
      result2 <-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))
      myplot2 <-ggplot(data=result2,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
        xlim(input$MetabMin, input$MetabMax)+ggtitle(input$MetabPlotVars[2])
      sum3 <- aggregate(MetabVar~responder+time, data=mydata3, mysummary)
      summaries3 <-sum3[,3]
      result3 <-as.data.frame(cbind(sum3[,c(1,2)],summaries3))
      names(result3)<-c("ResponderStatus","Time",colnames(summaries3))
      myplot3 <-ggplot(data=result3,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
        xlim(input$MetabMin, input$MetabMax)+ggtitle(input$MetabPlotVars[3])
      sum4 <- aggregate(MetabVar~responder+time, data=mydata4, mysummary)
      summaries4 <-sum4[,3]
      result4 <-as.data.frame(cbind(sum4[,c(1,2)],summaries4))
      names(result4)<-c("ResponderStatus","Time",colnames(summaries4))
      myplot4 <-ggplot(data=result4,aes(x=Time, y=Mean,group=ResponderStatus,colour=ResponderStatus))+geom_point()+geom_line(aes(size=.7))+
        geom_errorbar(aes(ymin = Mean - 2*Sd, ymax = Mean + 2*Sd,size=.7)) +scale_size(range=c(.1, 2),guide=FALSE)+
        xlim(input$MetabMin, input$MetabMax)+ggtitle(input$MetabPlotVars[4])

      if(input$MetabSamples){
        myplot <- myplot+geom_point(data=mydata,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=.25))
        myplot2 <- myplot2+geom_point(data=mydata2,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata2,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=.25))
        myplot3 <- myplot3+geom_point(data=mydata3,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata3,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=.25))
        myplot4 <- myplot4+geom_point(data=mydata4,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=1))+
          geom_line(data=mydata4,aes(x=time,y=MetabVar,group=subjectid,colour=responder,size=.25))
      }
    }

    if(length(input$MetabPlotVars) == 1){
      z = myplot
    }
    if(length(input$MetabPlotVars) == 2){
      z = multiplot(myplot, myplot2, cols = 2)
    }
    if(length(input$MetabPlotVars) == 3){
      z = multiplot(myplot, myplot2, myplot3, cols = 2)
    }
    if(length(input$MetabPlotVars) >= 4){
      z = multiplot(myplot, myplot2, myplot3, myplot4, cols = 2)
    }
    z
  })

  metabplot2 <- reactive({
    if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    if(length(input$MetabPlotVars) == 0){return(NULL)}
    if(length(input$MetabPlotVars) == 1){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      mydata<-mydata[order(mydata$time,mydata$responder),]
      mydata$time<-factor(mydata$time,ordered=T)
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),ordered=T)
      dummy<-data.frame(do.call(rbind,strsplit(levels(mydata$Time_Responder),"_")))
      index<-order(as.numeric(as.character(dummy$X1)))
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),levels(mydata$Time_Responder)[index])
      if(input$metabbox){
        myplot<-ggplot(data=mydata,aes(x=time,y=MetabVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$MetabPlotVars[1])
      }
      if(!input$metabbox){
        myplot<-ggplot(data=mydata,aes(x=Time_Responder,y=MetabVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$MetabPlotVars[1])
      }
    }
    if(length(input$MetabPlotVars) == 2){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      mydata<-mydata[order(mydata$time,mydata$responder),]
      mydata$time<-factor(mydata$time,ordered=T)
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),ordered=T)
      dummy<-data.frame(do.call(rbind,strsplit(levels(mydata$Time_Responder),"_")))
      index<-order(as.numeric(as.character(dummy$X1)))
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),levels(mydata$Time_Responder)[index])
      mydata2 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata2 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      mydata2 <-mydata2[order(mydata2$time,mydata2$responder),]
      mydata2$time<-factor(mydata2$time,ordered=T)
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),ordered=T)
      dummy2 <-data.frame(do.call(rbind,strsplit(levels(mydata2$Time_Responder),"_")))
      index2 <-order(as.numeric(as.character(dummy2$X1)))
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),levels(mydata2$Time_Responder)[index2])
      if(input$metabbox){
        myplot<-ggplot(data=mydata,aes(x=time,y=MetabVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$MetabPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=time,y=MetabVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$MetabPlotVars[2])
      }
      if(!input$metabbox){
        myplot<-ggplot(data=mydata,aes(x=Time_Responder,y=MetabVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$MetabPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=Time_Responder,y=MetabVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$MetabPlotVars[2])
      }
    }
    if(length(input$MetabPlotVars) == 3){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      mydata<-mydata[order(mydata$time,mydata$responder),]
      mydata$time<-factor(mydata$time,ordered=T)
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),ordered=T)
      dummy<-data.frame(do.call(rbind,strsplit(levels(mydata$Time_Responder),"_")))
      index<-order(as.numeric(as.character(dummy$X1)))
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),levels(mydata$Time_Responder)[index])
      mydata2 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata2 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      mydata2 <-mydata2[order(mydata2$time,mydata2$responder),]
      mydata2$time<-factor(mydata2$time,ordered=T)
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),ordered=T)
      dummy2 <-data.frame(do.call(rbind,strsplit(levels(mydata2$Time_Responder),"_")))
      index2 <-order(as.numeric(as.character(dummy2$X1)))
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),levels(mydata2$Time_Responder)[index2])
      mydata3 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata3 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      mydata3 <-mydata3[order(mydata3$time,mydata3$responder),]
      mydata3$time<-factor(mydata3$time,ordered=T)
      mydata3$Time_Responder<-factor(paste(mydata3$time,mydata3$responder,sep="_"),ordered=T)
      dummy3 <-data.frame(do.call(rbind,strsplit(levels(mydata3$Time_Responder),"_")))
      index3 <-order(as.numeric(as.character(dummy3$X1)))
      mydata3$Time_Responder<-factor(paste(mydata3$time,mydata3$responder,sep="_"),levels(mydata3$Time_Responder)[index3])
      if(input$metabbox){
        myplot<-ggplot(data=mydata,aes(x=time,y=MetabVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$MetabPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=time,y=MetabVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$MetabPlotVars[2])
        myplot3 <-ggplot(data=mydata3,aes(x=time,y=MetabVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$MetabPlotVars[3])
      }
      if(!input$metabbox){
        myplot<-ggplot(data=mydata,aes(x=Time_Responder,y=MetabVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$MetabPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=Time_Responder,y=MetabVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$MetabPlotVars[2])
        myplot3 <-ggplot(data=mydata3,aes(x=Time_Responder,y=MetabVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$MetabPlotVars[3])
      }
    }
    if(length(input$MetabPlotVars) >= 4){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      mydata<-mydata[order(mydata$time,mydata$responder),]
      mydata$time<-factor(mydata$time,ordered=T)
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),ordered=T)
      dummy<-data.frame(do.call(rbind,strsplit(levels(mydata$Time_Responder),"_")))
      index<-order(as.numeric(as.character(dummy$X1)))
      mydata$Time_Responder<-factor(paste(mydata$time,mydata$responder,sep="_"),levels(mydata$Time_Responder)[index])

      mydata2 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata2 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      mydata2 <-mydata2[order(mydata2$time,mydata2$responder),]
      mydata2$time<-factor(mydata2$time,ordered=T)
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),ordered=T)
      dummy2 <-data.frame(do.call(rbind,strsplit(levels(mydata2$Time_Responder),"_")))
      index2 <-order(as.numeric(as.character(dummy2$X1)))
      mydata2$Time_Responder<-factor(paste(mydata2$time,mydata2$responder,sep="_"),levels(mydata2$Time_Responder)[index2])

      mydata3 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata3 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      mydata3 <-mydata3[order(mydata3$time,mydata3$responder),]
      mydata3$time<-factor(mydata3$time,ordered=T)
      mydata3$Time_Responder<-factor(paste(mydata3$time,mydata3$responder,sep="_"),ordered=T)
      dummy3 <-data.frame(do.call(rbind,strsplit(levels(mydata3$Time_Responder),"_")))
      index3 <-order(as.numeric(as.character(dummy3$X1)))
      mydata3$Time_Responder<-factor(paste(mydata3$time,mydata3$responder,sep="_"),levels(mydata3$Time_Responder)[index3])

      mydata4 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[4])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata4 <-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[4])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      mydata4 <-mydata4[order(mydata4$time,mydata4$responder),]
      mydata4$time<-factor(mydata4$time,ordered=T)
      mydata4$Time_Responder<-factor(paste(mydata4$time,mydata4$responder,sep="_"),ordered=T)
      dummy4 <-data.frame(do.call(rbind,strsplit(levels(mydata4$Time_Responder),"_")))
      index4 <-order(as.numeric(as.character(dummy4$X1)))
      mydata4$Time_Responder<-factor(paste(mydata4$time,mydata4$responder,sep="_"),levels(mydata4$Time_Responder)[index4])

      if(input$metabbox){
        myplot<-ggplot(data=mydata,aes(x=time,y=MetabVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$MetabPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=time,y=MetabVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$MetabPlotVars[2])
        myplot3 <-ggplot(data=mydata3,aes(x=time,y=MetabVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$MetabPlotVars[3])
        myplot4 <-ggplot(data=mydata4,aes(x=time,y=MetabVar,fill=responder,colour=responder))+geom_boxplot(position = position_dodge(width = .9))+ggtitle(input$MetabPlotVars[4])
      }

      if(!input$metabbox){
        myplot<-ggplot(data=mydata,aes(x=Time_Responder,y=MetabVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$MetabPlotVars[1])
        myplot2 <-ggplot(data=mydata2,aes(x=Time_Responder,y=MetabVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$MetabPlotVars[2])
        myplot3 <-ggplot(data=mydata3,aes(x=Time_Responder,y=MetabVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$MetabPlotVars[3])
        myplot4 <-ggplot(data=mydata4,aes(x=Time_Responder,y=MetabVar,fill=responder,colour=responder))+geom_point()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle(input$MetabPlotVars[4])
      }
    }
    if(length(input$MetabPlotVars) == 1){
      z = myplot
    }
    if(length(input$MetabPlotVars) == 2){
      z = multiplot(myplot, myplot2, cols = 2)
    }
    if(length(input$MetabPlotVars) == 3){
      z = multiplot(myplot, myplot2, myplot3, cols = 2)
    }
    if(length(input$MetabPlotVars) >= 4){
      z = multiplot(myplot, myplot2, myplot3, myplot4, cols = 2)
    }
    z
  })

  output$MetabPlot<-renderPlot({
    print(metabplot())
  })

  output$MetabPlot2<-renderPlot({
    print(metabplot2())
  })

  MetabPlotSummary1<-reactive({
    if(input$MetabVarSet=="No Metab Data"){return(NULL)}
    if(length(input$MetabPlotVars) == 0){return(NULL)}
    if(length(input$MetabSub) == 0){return(NULL)}
    mysummary<-function(x){
      y<-c(mean(x),sd(x),length(x))
      names(y)<-c("Mean","Sd","N")
      return(y)
    }
    if(length(input$MetabPlotVars) == 1){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }
      sum1<-aggregate(MetabVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))
    }

    if(length(input$MetabPlotVars) == 2){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

      sum1<-aggregate(MetabVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))

      mydata2 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata2<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

      sum2<-aggregate(MetabVar~responder+time, data=mydata2, mysummary)
      summaries2<-sum2[,3]
      result2<-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))
    }

    if(length(input$MetabPlotVars) == 3){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

      sum1<-aggregate(MetabVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))

      mydata2 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata2<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

      sum2<-aggregate(MetabVar~responder+time, data=mydata2, mysummary)
      summaries2<-sum2[,3]
      result2<-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))

      mydata3 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata3<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

      sum3<-aggregate(MetabVar~responder+time, data=mydata3, mysummary)
      summaries3<-sum3[,3]
      result3<-as.data.frame(cbind(sum3[,c(1,2)],summaries3))
      names(result3)<-c("ResponderStatus","Time",colnames(summaries3))
    }

    if(length(input$MetabPlotVars) >= 4){
      mydata<-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[1])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

      sum1<-aggregate(MetabVar~responder+time, data=mydata, mysummary)
      summaries<-sum1[,3]
      result<-as.data.frame(cbind(sum1[,c(1,2)],summaries))
      names(result)<-c("ResponderStatus","Time",colnames(summaries))

      mydata2 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata2<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[2])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

      sum2<-aggregate(MetabVar~responder+time, data=mydata2, mysummary)
      summaries2<-sum2[,3]
      result2<-as.data.frame(cbind(sum2[,c(1,2)],summaries2))
      names(result2)<-c("ResponderStatus","Time",colnames(summaries2))

      mydata3 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata3<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[3])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

      sum3<-aggregate(MetabVar~responder+time, data=mydata3, mysummary)
      summaries3<-sum3[,3]
      result3<-as.data.frame(cbind(sum3[,c(1,2)],summaries3))
      names(result3)<-c("ResponderStatus","Time",colnames(summaries3))

      mydata4 <-data.frame(MetabVar=MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[4])],time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      if(input$MetabTransform){
        mydata4<-data.frame(MetabVar=log2(MetabData()[,which(names(MetabData()) %in% input$MetabPlotVars[4])]),time=MetabData()[,which(names(MetabData()) %in% values$m_time_var)],subjectid=MetabData()[,which(names(MetabData()) %in% values$m_patient_id)],responder=MetabData()[,which(names(MetabData()) %in% values$m_responder_var)])
      }

      sum4<-aggregate(MetabVar~responder+time, data=mydata4, mysummary)
      summaries4<-sum4[,3]
      result4<-as.data.frame(cbind(sum4[,c(1,2)],summaries4))
      names(result4)<-c("ResponderStatus","Time",colnames(summaries4))
    }

    if(length(input$MetabPlotVars) == 1){
      z <- result
    }

    if(length(input$MetabPlotVars) == 2){
      MetabVariable <- c(rep(input$MetabPlotVars[1], length(result[,1])), rep(input$MetabPlotVars[2], length(result2[,1])))
      z <- rbind(result, result2)
      z <- cbind(MetabVariable, z)
    }

    if(length(input$MetabPlotVars) == 3){
      MetabVariable <- c(rep(input$MetabPlotVars[1], length(result[,1])), rep(input$MetabPlotVars[2], length(result2[,1])), rep(input$MetabPlotVars[3], length(result3[,1])))
      z <- rbind(result, result2, result3)
      z <- cbind(MetabVariable, z)
    }

    if(length(input$MetabPlotVars) >= 4){
      MetabVariable <- c(rep(input$MetabPlotVars[1], length(result[,1])), rep(input$MetabPlotVars[2], length(result2[,1])), rep(input$MetabPlotVars[3], length(result3[,1])), rep(input$MetabPlotVars[4], length(result4[,1])))
      z <- rbind(result, result2, result3, result4)
      z <- cbind(MetabVariable, z)
    }

    z

  })

  output$normMetab<-renderUI({

    if(is.null(values$h3_rowdendro)){
      if(values$m_hc == TRUE){
        try<-list("Baseline Median Normalized" = 1,
                  "Baseline Healthy Normalized" = 2)
      }

      if(values$m_hc == FALSE){
        try<-list("Baseline Median Normalized" = 1)
      }
    }

    if(!is.null(values$h3_rowdendro)){
      if(values$m_hc == TRUE){
        try<-list("Baseline Median Normalized" = 1,
                  "Baseline Healthy Normalized" = 2,
                  "All Samples Median Normalized" = 3,
                  "All Samples Healthy Normalized"=4,
                  "All Samples Baseline Normalized"=5)
      }

      if(values$m_hc == FALSE){
        try<-list("Baseline Median Normalized" = 1,
                  "All Samples Median Normalized" = 3,
                  "All Samples Baseline Normalized"=5)
      }
    }

    selectInput("set3", "Select heatmap:",as.list(try))
  })

  design <- reactive({
    values$m_design
  })

  output$group_label_metab<-renderUI({
    selectInput("LabelMetab","Select variables to label samples (additional to order variables):",choices=names(values$m_design),selected=c(input$responder_var,input$patient_id),multiple=T)
  })

  output$TopTierMetab<-renderUI({selectInput("topmetab","Ordering columns: variable 1",
                                             c(names(values$m_design),"NA"),values$m_responder_var)
  })

  output$MidTierMetab<-renderUI({selectInput("midmetab","Ordering columns: variable 2",
                                             c(names(values$m_design),"NA"),values$m_patient_id)
  })

  output$LowTierMetab<-renderUI({selectInput("bottommetab","Ordering columns: variable 3",
                                             c(names(values$m_design),"NA"),"NA")
  })

  order_varsMetab <- reactive({
    x<-c(values$m_responder_var, values$m_patient_id, "NA")
    y <- x[which(x!="NA")]

    if(!is.null(input$topmetab) & !is.null(input$midmetab) & !is.null(input$bottommetab)){
      x<-c(input$topmetab,input$midmetab,input$bottommetab)
      y <- x[which(x!="NA")]
      return(y)
    }
    return(y)
  })

  output$subsetMetabVariable<-renderUI({selectInput("subsetMetabVar","Variable 1 to subset heatmap:",c(names(values$m_design)), values$m_responder_var)})
  output$subsetMetabValue<-renderUI({selectInput("subsetMetabVal","Value(s) of Variable to subset heatmap:",c(unique(as.character(values$m_design[,input$subsetMetabVar]))),c(unique(as.character(values$m_design[,input$subsetMetabVar])))[1],multiple=TRUE)})
  output$subsetMetabVariable2<-renderUI({selectInput("subsetMetabVar2","Variable 2 to subset heatmap:",c(names(values$m_design)), values$m_time_var)})
  output$subsetMetabValue2<-renderUI({selectInput("subsetMetabVal2","Value(s) of Variable to subset heatmap:",c(unique(as.character(values$m_design[,input$subsetMetabVar2]))),c(unique(as.character(values$m_design[,input$subsetMetabVar2])))[1],multiple=TRUE)})

  heatmapdatam <-reactive({
    design <- values$m_design
    metabdat <- values$metab_data
    metabvars <- as.character(values$metab_results$Metab.variable)
    final_expression = metabdat[,which(colnames(metabdat) %in% metabvars)]
    final_expression <- final_expression[,match(metabvars, colnames(final_expression), nomatch = 0)]
    rownames(final_expression) <- as.character(design$columnname)
    final_expression <- as.data.frame(t(final_expression))
    metabvars <- metabvars[match(rownames(final_expression), metabvars,nomatch = 0)]
    PROBE_ID <- as.character(metabvars)
    SYMBOL <- as.character(metabvars)
    final_expression <- log2(final_expression)
    final_expression <- do.call("cbind", list(PROBE_ID = PROBE_ID,SYMBOL=SYMBOL, final_expression))
    final_expression <- final_expression[which(final_expression$PROBE_ID %in% metablist()[,1]),]

    if (input$set3==1){#heatmapbase1
      base_sample_name=design$columnname[design[,m_baseline_var]==values$m_baseline_val]
      ind1<-which(colnames(final_expression)%in%c("PROBE_ID","SYMBOL"))
      ind2<-which(colnames(final_expression)%in%base_sample_name)
      exp_base_sam=final_expression[,c(ind1,ind2)]
      des_base_sam=design[which(design$columnname%in%colnames(exp_base_sam)),]
      y<-data.manipulate(exp=exp_base_sam,des=des_base_sam,values$m_baseline_val,longitudinal=FALSE,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=TRUE)
    }

    if (input$set3==2){
      base_sample_name=design$columnname[design[,m_baseline_var]==values$m_baseline_val]
      ind1<-which(colnames(final_expression)%in%c("PROBE_ID","SYMBOL"))
      ind2<-which(colnames(final_expression)%in%base_sample_name)
      exp_base_sam=final_expression[,c(ind1,ind2)]
      des_base_sam=design[which(design$columnname%in%colnames(exp_base_sam)),]
      y<-data.manipulate(exp=exp_base_sam,des=des_base_sam,values$m_control_var,values$m_control_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=FALSE)
    }
    if (input$set3==3){
      y<-data.manipulate(exp=final_expression,des=design[which(design$columnname %in% colnames(final_expression)),],values$m_baseline_var,values$m_baseline_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=TRUE)
    }
    if (input$set3==4){
      y<-data.manipulate(exp=final_expression,des=design[which(design$columnname %in% colnames(final_expression)),],values$m_control_var,values$m_control_val,longitudinal=FALSE,subjects,lg2=FALSE,keepbase=TRUE,format="Probes",allsamples=FALSE)
    }
    if (input$set3==5){#heatmap3
      if(values$m_hc==TRUE){
        des_w_controls<-design[which(design$columnname %in% colnames(final_expression)),]
        des_wo_controls<-design[-which(design[,values$m_control_var]==values$m_control_val),]
        h5index<-c(1,2,which(colnames(final_expression) %in% des_wo_controls$columnname))
        y<-data.manipulate(exp=final_expression[,h5index],des=des_wo_controls,values$m_baseline_var,values$m_baseline_val,longitudinal=TRUE,subjects=values$m_patient_id,lg2=FALSE,keepbase=FALSE,format="Probes",allsamples=FALSE)
      }
      if(values$m_hc==FALSE){
        y<-data.manipulate(exp=final_expression,des=design[which(design$columnname %in% colnames(final_expression)),],values$m_baseline_var,values$m_baseline_val,longitudinal=TRUE,subjects=values$m_patient_id,lg2=FALSE,keepbase=FALSE,format="Probes",allsamples=FALSE)
      }
    }
    z<-list(y=y)
    return(z)
  })

  heatmapnamem<-reactive({
    if(input$set3==1) heattxt<-"Baseline Median Normalized"
    if(input$set3==2) heattxt<-"Baseline Healthy Normalized"
    if(input$set3==3) heattxt<-"All Samples Median Normalized"
    if(input$set3==4) heattxt<-"All Samples Healthy Normalized"
    if(input$set3==5) heattxt<-"All Samples Normalized to each Subjects Baseline"
    heattxt
  })

  heatmap_orderm <- reactive({
    y<-heatmapdatam()$y
    if(input$subsetMetab){
      y$heatdes<-y$heatdes[which((y$heatdes[,as.character(input$subsetMetabVar)] %in% as.character(input$subsetMetabVal)) & (y$heatdes[,as.character(input$subsetMetabVar2)] %in% as.character(input$subsetMetabVal2))),]
      y$heatexp<-y$heatexp[,match(y$heatdes[,"columnname"], colnames(y$heatexp), nomatch=0)]
    }

    group_order<-order_varsMetab()
    z<-unique(c(order_varsMetab(),input$LabelMetab))
    color_groups<-z[which(z!="NA")]
    groups=y$heatdes[,color_groups]
    inside<-paste("y$heatdes$",group_order,sep="")
    des_order<-eval(parse(text=paste("order(",paste(inside,collapse=","),")",sep="")))
    design_ordered<-y$heatdes[des_order,]
    groups= as.data.frame(design_ordered[,color_groups, drop = F])
    z = list(y = y, groups = groups, color_groups = color_groups, design_ordered = design_ordered)
    return(z)
  })

  expression_matrixm <- reactive({
    y = heatmap_orderm()$y
    design_ordered = heatmap_orderm()$design_ordered
    x = y$heatexp[, match(design_ordered[,"columnname" ], colnames(y$heatexp), nomatch=0)]
    y$heatexp <- y$heatexp[, match(design_ordered[,"columnname" ], colnames(y$heatexp), nomatch=0)]
    colnames(y$heatexp) <- design_ordered[,values$m_sample_id]
    colnames(x)<-design_ordered[,values$m_sample_id]
    if(input$setcutoff2!=0){
      cut1<-as.numeric(input$setcutoff2)
      x[x>cut1]<-cut1
      x[x<(-cut1)]<--cut1}
    z = list(x = x, y = y)
    return(z)
  })

  opt_numClustm <- reactive({
    x = expression_matrixm()$x
    design_ordered = heatmap_orderm()$design_ordered
    colddm = NA
    dist_x = dist(t(x))
    hcl = fastcluster::hclust(dist_x)
    if(input$ColClustMetab==TRUE){
      colddm <- as.dendrogram(hcl)
    }
    d = sapply(2:round(nrow(design_ordered)/2), function(y) clValid::dunn(dist_x, cutree(hcl,y)))
    opt_num = which(d == max(d)) + 1
    z = list(opt_num = opt_num, colddm = colddm, hcl = hcl, d = d)
    return(z)
  })

  output$ClusterCutsM <- renderUI({
    numericInput('ClustCutMetab', "Number of clusters", min = 2, value = opt_numClustf()$opt_num, step = 1)
  })

  clusterxm <- reactive({
    groups = heatmap_orderm()$groups
    design_ordered = heatmap_orderm()$design_ordered
    hcl = opt_numClustm()$hcl
    if(input$ClusterChoiceM){
      clusters = cutree(hcl, input$ClustCutMetab)
      design_ordered$Clusters = as.character(clusters)
      groups <- cbind(groups, Cluster = design_ordered$Clusters)
    }
    z = list(groups = groups, design_ordered = design_ordered)
    return(z)
  })

  heatmap_colorsm <- reactive({
    groups = clusterxm()$groups
    color_groups = heatmap_orderm()$color_groups
    if(length(which(names(groups) %in% values$m_patient_id))){
      groups[,which(names(groups) %in% values$m_patient_id)]<-as.numeric(as.factor(groups[,which(names(groups) %in% values$m_patient_id)]))
    }
    if(values$m_time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i]) & groups[,i] != groups[,values$m_time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }
    else{
      for(i in 1:ncol(groups)){
        if(is.numeric(groups[,i])){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }
    for(i in 1:ncol(groups)){
      if(is.numeric(groups[,i]) == F){
        groups[,i] <- as.character((groups[,i]))
      }
    }
    if(values$m_time_var %in% colnames(groups)){
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10 & groups[,i] != groups[,values$m_time_var]){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }
    else{
      for(i in 1:ncol(groups)){
        if(length(unique(groups[,i])) > 10){
          groups[,i] <- as.numeric(as.factor(groups[,i]))
        }
      }
    }
    if(values$m_time_var %in% colnames(groups)){
      groups[,values$m_time_var] = as.factor(groups[,values$m_time_var])
    }
    n_groups = ncol(groups)
    palette_set = setPalettes(n_groups)
    factors = unlist(apply(groups, 2, function(x) as.data.frame(factor(x))), recursive=F)
    factor_levels = lapply(factors, getLevels)
    factor_lengths = lapply(factor_levels, function(x) length(x))
    color_gen = NULL
    add_list = list()
    for (i in 1:n_groups) {
      j = palette_set[[i]](as.numeric(factor_lengths[i]))
      color_gen = c(color_gen, j)
      add_list[[i]] = list(rect = list(col = "transparent", fill = j[factors[[i]]]))
    }
    annotation<-c()
    for(i in 1:n_groups){
      annotation<-cbind(annotation,add_list[[i]]$rect$fill)
    }

    annotation<-annotation[,n_groups:1]
    names(annotation)<-color_groups[n_groups:1]
    anno1<-colorRampPalette(c("navy", "yellow", "firebrick3"))(length(unique(groups[,1])))
    names(anno1)<-unique(groups[,1])
    eval(parse(text=paste(names(groups)[1],"=","anno1",sep="")))
    eval(parse(text=paste("first_color=list(",names(groups)[1],"=",names(groups)[1],")")))
    z <- list(groups = groups, first_color = first_color)
    return(z)
  })


  plotmw <- reactive({input$Metab_HeatMapSizew})
  plotmH <- reactive({input$Metab_HeatMapSizeh})
  plotsize_mW <- function(){plotmw()}
  plotsize_mH <- function(){plotmH()}
  Treeheightm <- reactive({input$Metab_TreeHeight})
  fontSizem <- reactive({input$Metab_FontSize})
  leg_sizem <- reactive({input$Metab_LegendSize})
  resolutionm <- function(){input$Metab_PlotRes}
  row_clusterm <- reactive({
    if(input$rowClustM == TRUE){
      rowclust <- TRUE
    }
    if(input$rowClustM == FALSE){
      rowclust <- NA
    }
    return(rowclust)
  })

  output$heatmapm <- renderPlot({
    if(all(expression_matrixm()$x == 0)){
      withProgress(message = '',
                   detail = 'Generating the Options...', value = 1,{
                     aheatmap2(expression_matrixm()$x,Rowv=row_clusterm(),Colv = opt_numClustm()$colddm,cexRow=1.2, treeheight = Treeheightm(), fontsize = fontSizem(), annheight = leg_sizem(), color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),annCol = heatmap_colorsm()$groups,annColors= heatmap_colorsm()$first_color,labRow=NA,breaks=0,legend = FALSE)
                   })
    }
    else{
      withProgress(message = '',
                   detail = 'Generating the Options...', value = 1,{
                     aheatmap2(expression_matrixm()$x,Rowv=row_clusterm(),Colv = opt_numClustm()$colddm, cexRow=1.2, treeheight = Treeheightm(), fontsize = fontSizem(), annheight = leg_sizem(), color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),annCol = heatmap_colorsm()$groups,annColors= heatmap_colorsm()$first_color,labRow=NA,breaks=0)
                   })
    }
  }, height = plotsize_mH, width = plotsize_mW)

  output$downloadHeatmapm <- downloadHandler(
    filename = function() {paste('SignificantMetab_Heatmap','.png', sep = '')},
    content = function(file){
      png(file, width = plotsize_mW(), height = plotsize_mH(), res = resolutionm())
      print(aheatmap2(expression_matrixm()$x,Rowv=TRUE,Colv=opt_numClustm()$colddm,cexRow=1.2, treeheight = input$Metab_TreeHeight, fontsize = input$Metab_FontSize, annheight = input$Metab_LegendSize, color = colorRampPalette(c("navy", "yellow", "firebrick3"))(50),annCol = heatmap_colorsm()$groups,annColors= heatmap_colorsm()$first_color,labRow=NA,breaks=0))
      dev.off()
    }
  )

  heatmap_download2M <- reactive({
    heatmapdata <- cbind(expression_matrix2()$symb, expression_matrix2()$y$heatexp)
    heatmapdata <- cbind(expression_matrix2()$probeids, heatmapdata)
    return(heatmapdata)
  })

  output$downloadHeatmap1 <- downloadHandler(
    filename = function() {paste(values$project_name,'_',heatmapname1(),'.csv', sep='')  },
    content = function(file) {
      write.csv(heatmap_download2M(), file, row.names = FALSE)
    }
  )

  output$MetabPlotSummary<-renderDataTable({
    MetabPlotSummary1()
  })

  output$downloadMetabResults <- downloadHandler(
    filename = function() {paste(values$project_name,'_',input$comparisonmetab,'.csv', sep='')  },
    content = function(file) {
      write.csv(MetabResult(), file,row.names = FALSE)
    }
  )

  output$downloadMetabData <- downloadHandler(
    filename = function() {paste(values$project_name,'_MetabData','.csv', sep='')  },
    content = function(file) {
      write.csv(MetabDataForDownload(), file, row.names = FALSE)
    }
  )

  output$downloadMetabSummaries <- downloadHandler(
    filename = function() {paste(values$project_name,'_',"MetabSummary_",input$MetabPlotVars[1],'.csv', sep='')  },
    content = function(file) {
      write.csv(MetabPlotSummary1(), file, row.names = FALSE)
    }
  )

  output$downloadMetabPlot <- downloadHandler(
    filename = function() {paste(values$project_name,'_','MetabPlot1','.png', sep = '')},
    content = function(file){
      png(file, width = 900)
      print(metabplot())
      dev.off()
    }
  )

  output$downloadMetabPlot2 <- downloadHandler(
    filename = function() {paste(values$project_name,'_','MetabPlot2','.png', sep = '')},
    content = function(file){
      png(file, width = 900)
      print(metabplot2())
      dev.off()
    }
  )

  ############# END OF METABOLOMICS PART ################

###############Correlations########################

  output$correlations <- renderMenu({
    if(is.null(values$correlations)){
      #return(menuItem(strong("Flow", style = "color:red"), tabName = ""))
      return(strong(""))
    }
    else{
      return(menuItem("Correlations", icon = icon("th-list"), tabName = "corr"))
    }
  })

  output$TypeVariable <- renderUI({
    selectInput("TypeVariable1", "Choose type:", choices = values$correlation_names, selected = values$correlation_names[1])
  })

  output$TypeVariable2 <- renderUI({
    selectInput("TypeVariable3", "Choose type:", choices = values$correlation_names, selected = values$correlation_names[1])
  })

  output$TypeVariable4 <- renderUI({
    selectInput("TypeVariable5", "Choose type:", choices = values$correlation_names, selected = values$correlation_names[1])
  })

  output$TypeVariable6 <- renderUI({
    selectInput("TypeVariable7", "Choose type:", choices = values$correlation_names, selected = values$correlation_names[1])
  })

  output$subsetcorr <- renderUI({
    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable1)]]
    if(input$var_switch == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    correlations <- correlations[which(correlations$Base_subtracted == input$Base_subtract),]

    visit.name = colnames(correlations)[3]
    visit <- unique(correlations[,3])
    visit <- gtools::mixedsort(as.character(visit))
    selectInput("subsetcorr1", paste("Subset by", " ",visit.name,":",sep = ""), choices = visit, selected = visit, multiple = TRUE)
  })

  output$WithVariable1 <- renderUI({
    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable1)]]
    if(input$var_switch == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    selectInput("WithVariable", "Choose 'with' variable:", choices = as.character(gtools::mixedsort(unique(correlations$With))), selected = gtools::mixedsort(unique(correlations$With))[1], multiple = TRUE)
  })

  output$WithVariable2 <- renderUI({
    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable3)]]
    if(input$var_switch2 == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    selectInput("WithVariable3", "Choose 'with' variable:", choices = as.character(gtools::mixedsort(unique(correlations$With))), selected = as.character(gtools::mixedsort(unique(correlations$With)))[1])
  })

  output$corr_Variable <- renderUI({
    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable7)]]
    selectInput("corr_Variable2", paste("Choose", values$y_var, "variable:", sep = " "), choices = as.character(gtools::mixedsort(unique(correlations$Variable))), selected = as.character(gtools::mixedsort(unique(correlations$With)))[1])
  })

  output$WithVariable4 <- renderUI({
    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable7)]]
    selectInput("WithVariable5", "Choose 'with' variable:", choices = as.character(gtools::mixedsort(unique(correlations$With))), selected = as.character(gtools::mixedsort(unique(correlations$With)))[1])
  })

  output$uploadmod1 <- renderUI({
    if(input$var_switch == FALSE){
      if(values$y_var[which(values$correlation_names == input$TypeVariable1)] == "Flow"){
        return(checkboxInput("uploadmod", strong(paste("Upload"," ",tolower(values$y_var[which(values$correlation_names == input$TypeVariable1)])," ", "variables",":",sep = ""), style = "color:blue"), FALSE))
      }
      else{
        return(checkboxInput("uploadmod", strong(paste("Upload"," ",tolower(values$y_var[which(values$correlation_names == input$TypeVariable1)]),"s",":",sep = ""), style = "color:blue"), FALSE))
      }
    }
    else{
      if(values$x_var[which(values$correlation_names == input$TypeVariable1)] == "Flow"){
        return(checkboxInput("uploadmod", strong(paste("Upload"," ",tolower(values$y_var[which(values$correlation_names == input$TypeVariable1)])," ", "variables",":",sep = ""), style = "color:blue"), FALSE))
      }
      else{
        return(checkboxInput("uploadmod", strong(paste("Upload"," ",tolower(values$x_var[which(values$correlation_names == input$TypeVariable1)]),"s",":",sep = ""), style = "color:blue"), FALSE))
      }
    }
  })

  output$fileupload1 <- renderUI({
    if(input$uploadmod == TRUE){
      return(fileInput('modselect', '', accept = ".csv"))
    }
    else{
      return(NULL)
    }
  })

  output$uploadmod3 <- renderUI({
    if(input$var_switch3 == FALSE){
      if(values$y_var[which(values$correlation_names == input$TypeVariable5)] == "Flow"){
        return(checkboxInput("uploadmod2", strong(paste("Upload"," ",tolower(values$y_var[which(values$correlation_names == input$TypeVariable5)])," ", "variables",":",sep = ""), style = "color:blue"), FALSE))
      }
      else{
        return(checkboxInput("uploadmod2", strong(paste("Upload"," ",tolower(values$y_var[which(values$correlation_names == input$TypeVariable5)]),"s",":",sep = ""), style = "color:blue"), FALSE))
      }
    }
    else{
      if(values$x_var[which(values$correlation_names == input$TypeVariable5)] == "Flow"){
        return(checkboxInput("uploadmod2", strong(paste("Upload"," ",tolower(values$x_var[which(values$correlation_names == input$TypeVariable5)])," ", "variables",":",sep = ""), style = "color:blue"), FALSE))
      }
      else{
        return(checkboxInput("uploadmod2", strong(paste("Upload"," ",tolower(values$x_var[which(values$correlation_names == input$TypeVariable5)]),"s",":",sep = ""), style = "color:blue"), FALSE))
      }
    }
  })

  output$fileupload2 <- renderUI({
    if(input$uploadmod2 == TRUE){
      return(fileInput('modselect2', '', accept = ".csv"))
    }
    else{
      return(NULL)
    }
  })

  h.dat <- reactive({
    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable1)]]
    correlation.type <- colnames(correlations)[4]

    correlations <- correlations[which(correlations$Base_subtracted == input$Base_subtract),]
    correlations <- data.table::as.data.table(correlations)
    correlations$FDR <- 1
    correlations$Bonf <- 1
    correlations[,FDR := p.adjust(Raw.P.Value, method = "fdr"), by = c(names(correlations)[c(2,3)])]
    correlations[,Bonf := p.adjust(Raw.P.Value, method = "bonferroni"), by = c(names(correlations)[c(2,3)])]

    correlations <- as.data.frame(correlations)
    z <- list(correlations = correlations)
  })

  heatmap.data <- reactive({

    if(input$uploadmod == FALSE){
      correlations <- h.dat()$correlations
      if(input$var_switch){
        correlations <- correlations[,c(2,1,3:ncol(correlations))]
        colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
      }
    }

    if(input$uploadmod == TRUE){
      correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable1)]]
      if(input$var_switch == TRUE){
        correlations <- correlations[,c(2,1,3:ncol(correlations))]
        colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
      }
      correlation.type <- colnames(correlations)[4]

      modnames.corr <- read.csv(input$modselect$datapath, header = TRUE)
      modnames.corr <- as.character(modnames.corr[,1])
      modnames.corr <- gsub(".", "_", modnames.corr,fixed = TRUE)

      correlations <- correlations[which(correlations$Base_subtracted == input$Base_subtract),]
      correlations <- correlations[which(correlations$Variable %in% modnames.corr),]

      col_anno <- list()
      col_anno2 <- list()
      dat <- list()
      dat2 <- list()
      keep_names <- list()

      names <- gtools::mixedsort(unique(as.character(correlations[,3])))

      if(input$subsetModcorr == TRUE){
        names <- input$subsetcorr1
      }

      for(j in 1:length(input$WithVariable)){
        correlations1 <- correlations[which(correlations$With %in% input$WithVariable[j]),]
        for(i in 1:length(names)){
          dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          dat[[i]] <- dat[[i]][which(dat[[i]]$Variable %in% modnames.corr),]
          dat[[i]] <- as.character(dat[[i]]$Variable)
        }

        keep_names[[j]] <- Reduce(intersect, dat)
      }

      keep_names <- unlist(keep_names)

      for(j in 1:length(input$WithVariable)){
        correlations1 <- correlations[which(correlations$With %in% input$WithVariable[j]),]
        for(i in 1:length(names)){
          dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          dat[[i]] <- dat[[i]][which(dat[[i]]$Variable %in% keep_names),]
          dat[[i]] <- dat[[i]][order(dat[[i]]$Variable),]
          dat[[i]] <- dat[[i]][,4]
          col_anno[[i]] <- input$WithVariable[j]
        }

        col_anno2[[j]] <- unlist(col_anno)
        dat2[[j]] <- do.call("cbind", dat)
        rownames(dat2[[j]]) <- sort(unique(keep_names))
        colnames(dat2[[j]]) <- names
      }

      col_anno <- data.frame(unlist(col_anno2))
      colnames(col_anno)[1] <- paste(values$x_var, "variable", sep = " ")
      dat <- do.call("cbind", dat2)
      z <- list(col_anno = col_anno, dat = dat)
    }

    else{
      correlations <- h.dat()$correlations
      if(input$var_switch){
        correlations <- correlations[,c(2,1,3:ncol(correlations))]
        colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
      }

      dat3 <- list()
      colnames3 <- list()
      rownames3 <- list()

      for(j in 1:length(input$WithVariable)){
        correlations1 <- correlations[which(correlations$With %in% input$WithVariable[j]),]

        dat <- list()
        sig <- list()
        cor_vals <- list()
        names <- gtools::mixedsort(unique(as.character(correlations1[,3])))

        if(input$subsetModcorr == TRUE){
          names <- input$subsetcorr1
        }

        if(input$corrval2 == TRUE){
          if(input$corrsign1 == "+"){
            for(i in 1:length(names)){
              dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
              cor_vals[[i]] <- dat[[i]][which(dat[[i]][,4] >= input$corrval3),]
              cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
            }
          }
          if(input$corrsign1 == "-"){
            for(i in 1:length(names)){
              dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
              cor_vals[[i]] <- dat[[i]][which(dat[[i]][,4] <= -input$corrval3),]
              cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
            }
          }
          if(input$corrsign1 == "Both"){
            for(i in 1:length(names)){
              dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
              cor_vals[[i]] <- dat[[i]][c(which(dat[[i]][,4] <= -input$corrval3),which(dat[[i]][,4] >= input$corrval3)),]
              cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
            }
          }

          cor_vals <- Reduce(union,cor_vals)
          correlations1 <- correlations1[which(correlations1$Variable %in% cor_vals),]
        }

        for(i in 1:length(names)){
          dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]

          if(input$correction_method.corr == "Raw"){
            sig[[i]] <- dat[[i]][which(dat[[i]]$Raw.P.Value <= input$Alpha1),]
          }

          if(input$correction_method.corr == "FDR"){
            sig[[i]] <- dat[[i]][which(dat[[i]]$FDR <= input$Alpha1),]
          }

          if(input$correction_method.corr == "Bonferroni"){
            sig[[i]] <- dat[[i]][which(dat[[i]]$Bonf <= input$Alpha1),]
          }

          sig[[i]] <- as.character(sig[[i]]$Variable)
        }

        sigvars <- Reduce(union,sig)

        correlations1 <- correlations1[which(correlations1$Variable %in% sigvars),]

        dat2 <- list()
        for(i in 1:length(names)){
          dat2[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          dat2[[i]] <- dat2[[i]][which(dat2[[i]]$Variable %in% sigvars),]$Variable
        }

        dat2 <- Reduce(intersect, dat2)

        for(i in 1:length(names)){
          dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          dat[[i]] <- dat[[i]][which(dat[[i]]$Variable %in% dat2),]
          dat[[i]] <- dat[[i]][order(dat[[i]]$Variable),]
          dat[[i]] <- dat[[i]][,4]
        }

        dat3[[j]] <- do.call("cbind", dat)
        colnames(dat3[[j]]) <- names
        rownames(dat3[[j]]) <- sort(dat2)
        colnames3[[j]] <- names
        rownames3[[j]] <- sort(dat2)
      }
      rownames3 <- Reduce(union, rownames3)
      dat2 <- list()
      col_anno <- list()
      col_anno2 <- list()

      for(j in 1:length(input$WithVariable)){
        correlations1 <- correlations[which(correlations$With %in% input$WithVariable[j]),]
        dat <- list()
        for(i in 1:length(names)){
          dat[[i]] <- correlations1[which(correlations1[,3] == names[i]),]
          dat[[i]] <- dat[[i]][which(dat[[i]]$Variable %in% rownames3),]
          dat[[i]] <- dat[[i]][order(dat[[i]]$Variable),]
          dat[[i]] <- dat[[i]][,4]
          col_anno[[i]] <- input$WithVariable[j]
        }

        col_anno2[[j]] <- unlist(col_anno)
        dat2[[j]] <- do.call("cbind",dat)
        rownames(dat2[[j]]) <- sort(rownames3)
        colnames(dat2[[j]]) <- names
      }
      dat <- do.call("cbind", dat2)
      col_anno <- as.data.frame(unlist(col_anno2))
      colnames(col_anno)[1] <- paste(values$x_var[which(values$correlation_names == input$TypeVariable1)], "variable", sep = " ")
      #dat <- do.call("cbind", dat2)
      if(input$ordercorr == TRUE){
        order_time <- mixedorder(colnames(dat))
        dat <- dat[,order_time]
        col_anno <- as.data.frame(col_anno[order_time,,drop = F])
      }
      z <- list(col_anno = col_anno, dat = dat)
    }

    return(z)
  })

  output$download_heatmap_data <- downloadHandler(
    filename = function() {paste(values$project_name, "_", "Correlation_heatmap_data.csv",sep = "")},
    content = function(file) {
      write.csv(heatmap.data()$dat, file, row.names = TRUE)
    }
  )

  row_clust.corr <- reactive({
    if(input$rowclustcorr == TRUE){return(TRUE)}
    else{
      return(NA)
    }
  })

  col_clust.corr <- reactive({
    if(input$colclustcorr == TRUE){return(TRUE)}
    else{
      return(NA)
    }
  })

  height1 <- function(){
    input$PlotHeight4
  }

  width1 <- function(){
    input$PlotWidth4
  }

  rows.plot <- reactive({
    if(nrow(heatmap.data()$dat) > 300){return(NA)}
    else{return(NULL)}
  })

  color.heatmap <- reactive({
    color_palette <- colorRampPalette(c("blue", "white", "red"))(1000)
    color_palette[c(495:505)] = "#FFFFFF"
    color_palette
  })

  plotresolution.corr <- reactive({
    input$PlotRes4
  })

  output$correlations_plotOverview <- renderPlot({
    if(ncol(heatmap.data()$dat) > 1){
      return(aheatmap2(heatmap.data()$dat,Colv = col_clust.corr(),annCol = heatmap.data()$col_anno,labRow = rows.plot(), color = color.heatmap(), border_color = "grey60", Rowv = row_clust.corr(), main = paste(input$TypeVariable1, "across visits", sep = " "), breaks = 0))
    }
    else{
      return(aheatmap2(heatmap.data()$dat,Colv = col_clust.corr(),labRow = rows.plot(), color = color.heatmap(), border_color = "grey60", Rowv = row_clust.corr(), main = paste(input$TypeVariable1, "across visits", sep = " "), breaks = 0))
    }
  }, height = height1, width = width1)

  output$download_corr_heatmap <- downloadHandler(
    filename = function() {paste(values$project_name,"_","Heatmap.png",sep = "")},
    content = function(file){
      png(file, width = (plotresolution.corr()/72)*width1(), height = (plotresolution.corr()/72)*height1(), res = plotresolution.corr())
      print(aheatmap2(heatmap.data()$dat, annCol = heatmap.data()$col_anno, Colv = col_clust.corr(),labRow = rows.plot(), color = color.heatmap(), border_color = "grey60", Rowv = row_clust.corr(), main = paste(input$TypeVariable1, "across visits", sep = " "), breaks = 0))
      dev.off()
    }
  )

  output$Visit <- renderUI({

    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable3)]]
    if(input$var_switch2 == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    correlations <- correlations[which(correlations$Base_subtracted == input$Base_subtract2),]

    visit.name = colnames(correlations)[3]
    visit <- unique(correlations[,3])
    visit <- gtools::mixedsort(as.character(visit))
    selectInput("visit", paste("Correlations by"," ",visit.name, ":", sep = ""), choices = visit, selected = visit[1])
  })

  sub.dat <- reactive({
    if(input$TypeVariable3 == input$TypeVariable1 & input$Base_subtract == input$Base_subtract2){
      correlations <- h.dat()$correlations
      if(input$var_switch2){
        correlations <- correlations[,c(2,1,3:ncol(correlations))]
        colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
      }
      correlations <- correlations[which(correlations$With == input$WithVariable3),]
      correlations <- correlations[which(correlations[,3] == input$visit),]
      z <- list(correlations = correlations)
    }
    else{
      correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable3)]]
      correlations <- correlations[which(correlations$Base_subtracted == input$Base_subtract2),]
      if(input$var_switch2){
        correlations <- correlations[,c(2,1,3:ncol(correlations))]
        colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
      }
      correlations <- correlations[which(correlations$With == input$WithVariable3),]
      correlations <- correlations[which(correlations[,3] == input$visit),]

      correlations$FDR <- p.adjust(correlations$Raw.P.Value, method = "fdr")
      correlations$Bonf <- p.adjust(correlations$Raw.P.Value, method = "bonferroni")
      z <- list(correlations = correlations)
    }
    return(z)
  })

  subset_correlations <- reactive({
    correlations <- sub.dat()$correlations

    if(input$corrval == TRUE){
      if(input$corrsign == "+"){
        correlations <- correlations[which(correlations[,4] >= input$corrval1),]
      }
      if(input$corrsign == "-"){
        correlations <- correlations[which(correlations[,4] <= -input$corrval1),]
      }
      if(input$corrsign == "Both"){
        correlations <- correlations[c(which(correlations[,4] <= -input$corrval1), which(correlations[,4] >= input$corrval1)),]
      }
    }

    if(input$correction_method.corr2 == "Raw"){
      correlations <- correlations[which(correlations$Raw.P.Value <= input$Alpha),]
    }

    if(input$correction_method.corr2 == "FDR"){
      correlations <- correlations[which(correlations$FDR <= input$Alpha),]
    }

    if(input$correction_method.corr2 == "Bonferroni"){
      correlations <- correlations[which(correlations$Bonf <= input$Alpha),]
    }

    correlations <- correlations[,-c(which(colnames(correlations) == "Base_subtracted"), which(colnames(correlations) == "NObs"),
                                     which(colnames(correlations) == "Sign_NegLog10_p"))]
    return(correlations)
  })

  output$download_data <- downloadHandler(
    filename = function() {paste(values$project_name, "_", "Correlations_subset.csv")},
    content = function(file) {
      write.csv(subset_correlations(), file, row.names = FALSE)
    }
  )

  output$correlation_table <- renderDataTable({
    subset_correlations()
  })

  output$Visit2 <- renderUI({

    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable5)]]
    if(input$var_switch3 == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    correlations <- correlations[which(correlations$Base_subtracted == input$Base_subtract3),]

    visit.name = colnames(correlations)[3]
    visit <- unique(correlations[,3])
    visit <- gtools::mixedsort(as.character(visit))
    selectInput("visit2", paste("Correlations by", visit.name, ":", sep = " "), choices = visit, selected = visit[1])
  })

  output$subsetcorr2 <- renderUI({

    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable5)]]
    if(input$var_switch3 == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    correlations <- correlations[which(correlations$Base_subtracted == input$Base_subtract3),]
    var <- as.character(sort(unique(correlations$With)))
    selectInput("subsetcorr3", paste("Subset by 'with' variable:",sep = ""), choices = var, selected = var[1:5], multiple = TRUE)
  })

  correlation.data <- reactive({

    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable5)]]
    if(input$var_switch3 == TRUE){
      correlations <- correlations[,c(2,1,3:ncol(correlations))]
      colnames(correlations) <- colnames(correlations)[c(2,1,3:ncol(correlations))]
    }
    correlations <- correlations[which(correlations$Base_subtracted == input$Base_subtract3),]
    correlations <- correlations[which(correlations[,3] == input$visit2),]

    if(input$subsetModcorr2 == TRUE){
      correlations <- correlations[which(correlations$With %in% input$subsetcorr3),]
    }

    correlations$FDR <- 1
    correlations$Bonf <- 1
    if(input$var_switch3 == FALSE){
      for(i in 1:length(unique(correlations$With))){
        correlations$FDR[which(correlations$With == unique(correlations$With)[i])] <- p.adjust(correlations$Raw.P.Value[which(correlations$With == unique(correlations$With)[i])], "fdr")
        correlations$Bonf[which(correlations$With == unique(correlations$With)[i])] <- p.adjust(correlations$Raw.P.Value[which(correlations$With == unique(correlations$With)[i])], "bonferroni")
      }
    }
    else{
      for(i in 1:length(unique(correlations$Variable))){
        correlations$FDR[which(correlations$Variable == unique(correlations$Variable)[i])] <- p.adjust(correlations$Raw.P.Value[which(correlations$Variable == unique(correlations$Variable)[i])], "fdr")
        correlations$Bonf[which(correlations$Variable == unique(correlations$Variable)[i])] <- p.adjust(correlations$Raw.P.Value[which(correlations$Variable == unique(correlations$Variable)[i])], "bonferroni")
      }
    }

    if(input$var_switch3 == TRUE){
      y <- values$x_var[which(values$correlation_names == input$TypeVariable5)]
      x <- values$y_var[which(values$correlation_names == input$TypeVariable5)]
    }
    else{
      y <- values$y_var[which(values$correlation_names == input$TypeVariable5)]
      x <- values$x_var[which(values$correlation_names == input$TypeVariable5)]
    }
    dat <- list()
    cor_vals <- list()
    sig <- list()
    names <- as.character(unique(correlations$With))

    if(input$corrval4 == TRUE){
      if(input$corrsign2 == "+"){
        for(i in 1:length(names)){
          dat[[i]] <- correlations[which(correlations$With == names[i]),]
          cor_vals[[i]] <- dat[[i]][which(dat[[i]][,4] >= input$corrval5),]
          cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
        }
      }
      if(input$corrsign2 == "-"){
        for(i in 1:length(names)){
          dat[[i]] <- correlations[which(correlations$With == names[i]),]
          cor_vals[[i]] <- dat[[i]][which(dat[[i]][,4] <= -input$corrval5),]
          cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
        }
      }
      if(input$corrsign2 == "Both"){
        for(i in 1:length(names)){
          dat[[i]] <- correlations[which(correlations$With == names[i]),]
          cor_vals[[i]] <- dat[[i]][c(which(dat[[i]][,4] <= -input$corrval5),which(dat[[i]][,4] >= input$corrval5)),]
          cor_vals[[i]] <- as.character(cor_vals[[i]]$Variable)
        }
      }

      cor_vals <- Reduce(union,cor_vals)
      correlations <- correlations[which(correlations$Variable %in% cor_vals),]
    }


    for(i in 1:length(names)){
      dat[[i]] <- correlations[which(correlations$With == names[i]),]
      if(input$correction_method.corr3 == "Raw"){
        sig[[i]] <- dat[[i]][which(dat[[i]]$Raw.P.Value <= input$Alpha2),]
      }
      if(input$correction_method.corr3 == "FDR"){
        sig[[i]] <- dat[[i]][which(dat[[i]]$FDR <= input$Alpha2),]
      }
      if(input$correction_method.corr3 == "Bonferroni"){
        sig[[i]] <- dat[[i]][which(dat[[i]]$Bonf <= input$Alpha2),]
      }
      sig[[i]] <- as.character(sig[[i]]$Variable)
    }

    sigvars <- Reduce(union,sig)

    correlations <- correlations[which(correlations$Variable %in% sigvars),]

    correlation.type <- colnames(correlations)[4]

    if(input$uploadmod2 == TRUE){
      modnames.corr <- read.csv(input$modselect2$datapath, header = TRUE)
      modnames.corr <- as.character(modnames.corr[,1])
      modnames.corr <- gsub(".", "_", modnames.corr,fixed = TRUE)
      correlations <- correlations[which(correlations$Variable %in% modnames.corr),]
    }

    n = length(unique(correlations$With))
    cor <- list()
    cor_sub <- list()
    names <- unique(correlations$With)

    for(i in 1:n){
      cor[[i]] <- correlations[which(correlations$With == names[i]),]

      if(i > 1){
        cor[[i]] <- cor[[i]][match(cor[[i-1]]$Variable, cor[[i]]$Variable, nomatch = 0),]
      }

      if(input$correction_method.corr3 == "Raw"){
        cor_sub[[i]] <- cor[[i]][,which(colnames(cor[[i]]) %in% c(correlation.type, "Raw.P.Value"))]
      }
      if(input$correction_method.corr3 == "FDR"){
        cor_sub[[i]] <- cor[[i]][,which(colnames(cor[[i]]) %in% c(correlation.type, "FDR"))]
      }
      if(input$correction_method.corr3 == "Bonferroni"){
        cor_sub[[i]] <- cor[[i]][,which(colnames(cor[[i]]) %in% c(correlation.type, "Bonf"))]
      }
      colnames(cor_sub[[i]]) <- c(paste(names[i], "r.val", sep = "."), paste(names[i], "pVal", sep = "."))
    }

    mat1 <- do.call("cbind", cor_sub)
    mat1 <- as.matrix(mat1)
    rownames(mat1) <- cor[[1]]$Variable
    mat1 <- as.data.frame(mat1)

    t = mat1 #mat1 is what you created with the correlation script.
    t$modName<-rownames(t)
    t.pval=reshape2::melt(t[c(length(t),grep("pVal",colnames(t)))]) #combine all columns with "pVal" into one loooooong object.
    t.rval=reshape2::melt(t[c(length(t),grep("r.val",colnames(t)))]) #combine all columns with "r.val" into one loooooong object.
    GraphFrame=cbind(t.pval,t.rval[,3])  #concaternate them into one matrix
    if(input$correction_method.corr3 == "Raw"){
      colnames(GraphFrame)=c(y,"Correlation","Raw.pVal","rVal")  #change column names
    }
    if(input$correction_method.corr3 == "FDR"){
      colnames(GraphFrame)=c(y,"Correlation","FDR.pVal","rVal")  #change column names
    }
    if(input$correction_method.corr3 == "Bonferroni"){
      colnames(GraphFrame)=c(y,"Correlation","Bonf.pVal","rVal")  #change column names
    }

    GraphFrame[,which(colnames(GraphFrame) == y)] = factor(GraphFrame[,which(colnames(GraphFrame) == y)],levels=unique(GraphFrame[,which(colnames(GraphFrame) == y)],order=T))  #this is super important as it makes sure the order on the axis is correct!
    GraphFrame$Correlation <- gsub(".pVal", "",GraphFrame$Correlation, fixed = TRUE)
    z <- list(GraphFrame = GraphFrame, y = y, x = x)
    return(z)
  })

  output$download_plot_data <- downloadHandler(
    filename = function() {paste(values$project_name, "_", "correlations_data.csv")},
    content = function(file) {
      write.csv(correlation.data()$GraphFrame, file, row.names = FALSE)
    }
  )

  height <- function(){
    input$PlotHeight3
  }

  width <- function(){
    input$PlotWidth3
  }

  plotresolution.corr1 <- reactive({
    input$PlotRes3
  })


  correlations_makePlot <- reactive({
    y <- correlation.data()$y
    x <- correlation.data()$x
    GraphFrame <- correlation.data()$GraphFrame
    if(input$correction_method.corr3 == "Raw"){
      GraphFrame <- GraphFrame[-which(GraphFrame$Raw.pVal >= input$Alpha2),]
      p <- GraphFrame$Raw.pVal
    }
    if(input$correction_method.corr3 == "FDR"){
      GraphFrame <- GraphFrame[-which(GraphFrame$FDR.pVal >= input$Alpha2),]
      p <- GraphFrame$FDR.pVal
    }
    if(input$correction_method.corr3 == "Bonferroni"){
      GraphFrame <- GraphFrame[-which(GraphFrame$Bonf.pVal >= input$Alpha2),]
      p <- GraphFrame$Bonf.pVal
    }
    p[p < .0001] <- .0001 #set log-10 transformed p-values <4 (=0.0001) to 4
    RV <- ifelse(GraphFrame$rVal>0,GraphFrame$rVal^4,(GraphFrame$rVal^4)*-1)

    if(length(p)/length(unique(GraphFrame$Correlation)) <= 260){
      graph <- ggplot(data = GraphFrame, aes(Correlation, GraphFrame[,1],colour=RV, size = p), environment = environment()) + labs(y = y, x = x)  + scale_size_continuous("P-value",range=c(12,4),limits = c(min(p), .06), breaks = c(round(min(p),4),round((min(p)+.05)/2,4),.05)) + geom_point(colour="black",aes(size = p),shape=21,alpha=I(1)) + geom_point(alpha=I(.75))  + theme_minimal()  +
        theme(panel.background=element_rect(fill="white"),
              panel.grid = element_line(),panel.grid.major.y=element_blank(),
              panel.grid.major.x = element_blank(),
              axis.text.y=element_text(hjust=0),
              axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
        scale_x_discrete(expand=c(0,1)) +
        scale_color_gradient2(paste(values$correlation_method,"'s r",sep = ""),GraphFrame$rVal, low="navy", high="red", breaks=c(min(RV),0,max(RV)), labels = c(as.character(round(-1*(-1*min(RV))^(1/4),3)),"0",as.character(round(max(RV)^(1/4),3))))
    }

    if(length(p)/length(unique(GraphFrame$Correlation)) > 260){
      graph <- ggplot(data = GraphFrame, aes(Correlation, GraphFrame[,1],colour=RV, size = p), environment = environment()) + labs(y = y, x = x)  + scale_size_continuous("P-value",range=c(12,4),limits = c(min(p), .06), breaks = c(round(min(p),4),round((min(p)+.05)/2,4),.05)) + geom_point(colour="black",aes(size = p),shape=21,alpha=I(1)) + geom_point(alpha=I(.75))  + theme_minimal()  +
        theme(panel.background=element_rect(fill="white"),
              panel.grid = element_line(),panel.grid.major.y=element_blank(),
              panel.grid.major.x = element_blank(),
              axis.text.y=element_blank(),
              axis.text.x=element_text(angle=45,vjust=1,hjust=1)) +
        scale_x_discrete(expand=c(0,1)) + scale_y_discrete(breaks = NULL) +
        scale_color_gradient2(paste(values$correlation_method,"'s r",sep = ""),GraphFrame$rVal, low="navy", high="red", breaks=c(min(RV),0,max(RV)), labels = c(as.character(round(-1*(-1*min(RV))^(1/4),3)),"0",as.character(round(max(RV)^(1/4),3))))
    }

    return(graph)
  })

  output$correlations_plot <- renderPlot({
    correlations_makePlot()
  }, height = height, width = width)

  output$download_corr_plot <- downloadHandler(
    filename = function() {paste(values$project_name, "_","significant_correlations_plot.png")},
    content = function(file){
      png(file, width = (plotresolution.corr1()/72)*width(), height = (plotresolution.corr1()/72)*height(), res = plotresolution.corr1())
      print(correlations_makePlot())
      dev.off()
    }
  )

  output$Visit3 <- renderUI({

    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable7)]]
    visit.name = colnames(correlations)[3]
    visit <- unique(correlations[,3])
    visit <- gtools::mixedsort(as.character(visit))
    selectInput("visit3", paste("Correlations by"," ",visit.name, ":", sep = ""), choices = visit, selected = visit[1])
  })


  scatter_plot <- reactive({
    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable7)]]
    correlation_file <- values$correlation_files[[which(values$correlation_names == input$TypeVariable7)]]
    x <- correlation_file[,input$WithVariable5][which(correlation_file[,colnames(correlations)[3]] == input$visit3)]
    y <- correlation_file[,input$corr_Variable2][which(correlation_file[,colnames(correlations)[3]] == input$visit3)]
    z <- list(x = x, y = y)
    return(z)
  })

  scatter_width <- function(){
    input$PlotWidth5
  }

  scatter_height <- function(){
    input$PlotHeight5
  }

  plot_res5 <- reactive({
    input$PlotRes5
  })


  axis_text_size <- reactive({
    input$axis_text_size
  })

  axis_label_size <- reactive({
    input$axis_label_size
  })

  scatter_makePlot <- reactive({
    correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable7)]]
    if(input$log_scale == TRUE){
      dat <- data.frame(x = scatter_plot()$x, y = log2(scatter_plot()$y + 1))
    }
    else{
      dat <- data.frame(x = scatter_plot()$x, y = scatter_plot()$y)
    }
    scat.plot <-  ggplot(data = dat, aes(x = x, y = y)) + geom_point(size = input$Point_size) + theme_bw()+
      labs(x = paste(input$WithVariable5, colnames(correlations)[3], input$visit3, sep = " "), y = paste(input$corr_Variable2, colnames(correlations)[3], input$visit3, sep = " ")) +
      theme(axis.text = element_text(size = axis_text_size()), axis.title = element_text(size = axis_label_size()))

    if(input$plot_reg == TRUE){
      scat.plot <- scat.plot + geom_smooth(method = "lm", se = FALSE)
    }
    if(input$plot_loess == TRUE){
      if(input$plot_reg == TRUE){
        scat.plot <- scat.plot + geom_smooth(method = "lm", se = FALSE)
      }
      else{
        scat.plot <- scat.plot + geom_smooth(method = "loess", se = FALSE, span = input$span)
      }
    }
    z <- list(scatplot = scat.plot, dat = dat)
    return(z)
  })


  output$correlations_scatter_plot <- renderPlot({
    scatter_makePlot()$scatplot
  }, width = scatter_width, height = scatter_height)


  output$download_scatter_data <- downloadHandler(
    filename = function() {correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable7)]]
    paste(values$project_name, "_",input$corr_Variable2, "vs", input$WithVariable5, "_", colnames(correlations)[3], input$visit3,"_data.csv")},
    content = function(file) {
      dat <- scatter_makePlot()$dat
      colnames(dat) <- c(input$WithVariable5, input$corr_Variable2)
      write.csv(dat, file, row.names = FALSE)
    }
  )

  output$download_scatter <- downloadHandler(
    filename = function() {correlations <- values$correlations[[which(values$correlation_names == input$TypeVariable7)]]
    paste(values$project_name, "_",input$corr_Variable2, "vs", input$WithVariable5, "_", colnames(correlations)[3], input$visit3,".png")},
    content = function(file){
      png(file, width = (plot_res5()/72)*width(), height = (plot_res5()/72)*height(), res = plot_res5())
      print(scatter_makePlot()$scatplot)
      dev.off()
    }
  )



################End of Correlations Part###################

  })
