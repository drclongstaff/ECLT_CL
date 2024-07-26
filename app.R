
library(shiny)
#library(purrr)
#library(dplyr)

#function for truncating data
trunc_plateFun <- function(tPLATE, TR){
  shortPlate <- head(tPLATE, n=nrow(tPLATE)-TR)
}
#function for adjusting baseline
BaselineOff2 <- function(n, off) {
  n <- n-min(n)-off
}

#function for raw data plot
base_plotFun <- function(PLATE, NUMR) {
  
  Time<-PLATE[[1]] #time in the first column
  plateData<-PLATE[,-1] #the absorbance data without the time column
  absWells <- length(plateData[1,]) #the no. of columns of the absorbance data
  mint<-min(Time, na.rm = TRUE)
  maxt<-max(Time, na.rm = TRUE)
  maxy<-max(plateData, na.rm = TRUE) #min and max values for scaling the plots
  samples <- colnames(plateData) #names of columns to go on the plots
  #plot dimensions
  par(mfrow=c(NUMR,length(samples)/NUMR))
  par(mar=c(0.2,0.2,0.2,0.2)) # dimensions for figure
  #function to generate multiple plots
  lapply(seq_along(plateData), function(i) {
    
    plot(x = Time, y = plateData[, i], type = "l", col="steelblue", lwd=3,
         ylim=c(0, maxy*1.2), xaxt="n", yaxt="n", xlab = "" , ylab="")
    legend("topright", bty="n", paste0(samples[i],"=",i), cex=1.5 )
    
  })
}
#function for analysed plots
simple_plotFun <- function(PLATE, NUMR, TABRES) {
  
  Time<-PLATE[[1]] #time in the first column
  plateData<-PLATE[,-1] #the absorbance data without the time column
  absWells <- length(plateData[1,]) #the no. of columns of the absorbance data
  mint<-min(Time, na.rm = TRUE)
  maxt<-max(Time, na.rm = TRUE)
  maxy<-max(plateData, na.rm = TRUE) #min and max values for scaling the plots
  samples <- colnames(plateData) #names of columns to go on the plots
 
  par(mfrow=c(NUMR,length(samples)/NUMR))
  par(mar=c(0.2,0.2,0.2,0.2)) # dimensions for figure
  
  lapply(seq_along(plateData), function(i) {
    
    plot(x = Time, y = plateData[, i], type = "l", col="blue2", lwd=3,
    ylim=c(0, maxy*1.2), xaxt="n", yaxt="n", xlab = "" , ylab="")
    legend("topright", bty="n", paste0(samples[i],"=",i), cex=1.5 )
    endPoint <- TABRES[i,3]
    #endPoint <- which(Time==TABRES[i,2])
    #ifelse(!is.na(endPoint), endPoint <-endPoint, 2 )
    lines(Time[1:endPoint], plateData[,i][1:endPoint],col="tomato", lwd=3)
    abline("v"= TABRES[i,2], lty=2)
    abline("h"= TABRES[i,4], lty=2)
  })
}

#Function for  downcurve analysis  
downy<-function(d, Time, ini, thresh, off){ #d is the absorbance data
 
  minAbs <- min(d, na.rm = TRUE)+off #off is baseline offset
  maxAbs <- max(d, na.rm = TRUE)
  pointmax<-which.max(d)
  pcChange<-ini*(maxAbs-minAbs)+minAbs #ini is the |> % lysis, set at 50% here
  #need to deal with curves that don't go to 50% lysis
  ifelse(d[length(d)]>=pcChange, downTime <- Time, downTime<-Time[-c(1:pointmax)])
  ifelse(d[length(d)]>=pcChange, downAbs <- d,  downAbs<-d[-c(1:pointmax)] )
  #TC deals with curves that increase at the end, after lysis and should be ignored
  TC <- which.min(downAbs)
  downTime <- downTime[1:TC] #only go to min after lysis and discard later points
  downAbs <- downAbs[1:TC]
  #decaypoint is where set% lysis occurs
  ifelse(d[length(d)]>=pcChange, decayPoint <- length(d), decayPoint<-which(abs(downAbs-pcChange)==min(abs(downAbs-pcChange)))[1] )
  
  ifelse((max(d)-min(d)<thresh | min(downAbs)>=pcChange),
         decayAbs <- NA,
         decayAbs<-downAbs[decayPoint]#,using nearest point not interpolation
         # decayAbs<-round(approx(downTime, downAbs, xout = pcChange, ties = mean)$x,3)
  )
  
  #StartTime is fitted if abs > threshold, otherwise is closest point
  #This prevents crashing if there are blank wells
  ifelse((max(d)-min(d)<thresh | min(downAbs)>=pcChange),
         decayTime <- NA,
         decayTime<-downTime[decayPoint]#,
         # decayTime<-round(approx(downAbs, downTime, xout = decayAbs, ties = mean)$y,3)
  )
  #vector of outputs from the function
  #downcurve <- c(decayAbs, decayTime, decayPoint+pointmax,lastPoint, endTime, minAbs)
  downcurve <- c(decayAbs, decayTime, decayPoint+pointmax,minAbs)
}

ui <- fluidPage(
  #a js to measure the speed of the shiny app operation
  tags$script(
    src = "https://cdn.jsdelivr.net/gh/Appsilon/shiny.tictoc@v0.2.0/shiny-tic-toc.min.js"
  ),
  h3(id="Title", "Simple analysis of clot lysis curves, version 0.13"),
  helpText(
    tags$a(href = "https://github.com/drclongstaff/shiny-clots/blob/master/docs/ECLT-app-notes.pdf", 
           "help notes", target = "_blank")
          ),
  helpText(h5("Load a csv file, check the raw data and remove noisy wells")),
  fluidRow( column(4,fileInput("file", "Upload data file (CSV)")),
            column(4,textInput("remove_cols", "Remove column nos (comma-separated):", "-1")),
            column(4, helpText("Removed"),textOutput("remove_txt"))
          ),
  helpText(h5("Modify the baseline and truncate the data as necessary")),
  fluidRow(
    column(3,numericInput("thresh", "Threshold", "0.02", step=0.01)),
    column(3,numericInput("trunc", "Truncate points", "0", step=10)),
    column(3, numericInput("off", "offset zero baseline", "0", step=0.01)),
    column(3,numericInput("numr", "Plot number of rows", "8"))
          ),
  
  fluidRow(
  column(3, radioButtons(
    inputId = "plotsab",
    label = NULL,#"Analysis generates curves and table of times to 50% lysis",
    choices = c("Raw", 
                "Analysed"),
    inline = TRUE)
    #helpText(h4("Analysed generates curves and table of times to 50% lysis")),
          ),
  column(6, helpText(h5("Analysed generates curves and table of times to 50% lysis"))),
          ),
  
  fluidRow(
  column(12,plotOutput(outputId = "simpleplot"))
          ),
 
  tableOutput("data_table")
  
    )

server <- function(input, output) {
    #simple loading data function
   data <- reactive({
        req(input$file)
        inputFile <- input$file
        df <- read.csv(input$file$datapath, header = TRUE)
        #df <- openxlsx::read.xlsx(input$file$datapath, sheet = 1)
        #remove any non-numeric columns
        df_numeric <- df[, sapply(df, is.numeric)]
        return(df_numeric)
        #return(df)
          })
   
    #function to remove selected columns
    remove_cols <- reactive({
    req(input$remove_cols)
    rem <- as.numeric(strsplit(input$remove_cols, ",")[[1]])
    rem+1
          })
    #make list of removed wells to show
    col_names <- reactive({
      req(input$remove_cols, data())
          colnames(data())[remove_cols()]
          })
    #remove the columns of data selected
    filtered_data <- reactive({
      req(data(), remove_cols())
      
      filtered <- data()
      
      # Ensure we're only working with valid column indices
      valid_cols <- remove_cols()[remove_cols() > 0 & remove_cols() <= ncol(filtered)]
      
      # Set the selected columns to 0
      if (length(valid_cols) > 0) {
        filtered[, valid_cols] <- 0
      }
      
      filtered
    })
    #text of removed wells
    output$remove_txt <- renderText({
    col_names()
          })
    #truncate the data
    truncated_data <- reactive({
    trunc_plateFun(filtered_data(), input$trunc)
          })
    #the zeroing function
    zerod_data <- reactive({
      # Apply BaselineOff2 to all columns except the first one
      zerod_list <- lapply(truncated_data()[,-1], function(x) BaselineOff2(x, input$off))
      # Convert the list to a data frame
      zerod_df <- as.data.frame(zerod_list)
      # Add the Time column back
      cbind("Time" = truncated_data()[[1]], zerod_df)
    })
    #choose the appropriate plot
    output$simpleplot <- renderPlot({
    switch(input$plotsab,
             "Raw"=base_plotFun(filtered_data(), input$numr),
             "Analysed"=simple_plotFun(zerod_data(), input$numr, TabResdown())
                   )
          })
    #perform the calculation using downy function and purrr
    TabResdown <- reactive({
      args_list <- list(Time = zerod_data()[[1]], ini = 0.5, thresh = input$thresh, off = input$off)
      
      result <- lapply(zerod_data()[-1], function(x) {
        res <- do.call(downy, c(list(x), args_list))
        data.frame(
          lys.abs = res[1],
          lys.time = res[2],
          decayPoint = res[3],
          min.Abs = res[4]
        )
      })
      
      do.call(rbind, result)
    })
    #make the appropriate table
    output$data_table <- renderTable({
      switch(input$plotsab,
           "Raw"=data(),
           #"Analysed"=matrix((TabResdown() %>% select(lys.time) %>% pull()), nrow = input$numr, byrow = TRUE)
           "Analysed"=matrix(TabResdown()$lys.time, nrow = input$numr, byrow = TRUE)
           )
        })
}

shinyApp(ui, server)