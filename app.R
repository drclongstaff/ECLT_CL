
library(shiny)
library(purrr)
library(dplyr)


trunc_plateFun <- function(tPLATE, TR){
  shortPlate <- head(tPLATE, n=nrow(tPLATE)-TR)
}

BaselineOff2 <- function(n, off) {
  n <- n-min(n)-off
}

simple_plotFun <- function(PLATE, NUMR, TABRES) {
  
  Time<-PLATE[[1]] #time in the first column
  plateData<-PLATE[,-1] #the absorbance data without the time column
  absWells <- length(plateData[1,]) #the no. of columns of the absorbance data
  mint<-min(Time, na.rm = TRUE)
  maxt<-max(Time, na.rm = TRUE)
  maxy<-max(plateData, na.rm = TRUE) #min and max values for scaling the plots
  samples <- colnames(plateData) #names of columns to go on the plots
 
  #par(mfrow=c(8,12))
  par(mfrow=c(NUMR,length(samples)/NUMR))
  par(mar=c(0.2,0.2,0.2,0.2)) # dimensions for figure
  #ifelse((is.na(TABRES$lys.time) | TABRES$lys.time==0), 10, TABRES$lys.time)
  lapply(seq_along(plateData), function(i) {
    
    plot(x = Time, y = plateData[, i], type = "l", col="blue2", lwd=3,
    ylim=c(0, maxy*1.2), xaxt="n", yaxt="n", xlab = "" , ylab="")
    legend("topright", bty="n", paste0(samples[i],"=",i), cex=1.5 )
    endPoint <- TABRES[i,3]
    #endPoint <- which(Time==TABRES[i,2])
    #ifelse(!is.na(endPoint), endPoint <-endPoint, 2 )
    lines(Time[1:endPoint], plateData[,i][1:endPoint],col="tomato", lwd=3)
    abline("v"= TABRES[i,2], lty=2)
    abline("h"= TABRES[i,6], lty=2)
  })
 
}

base_plotFun <- function(PLATE, NUMR) {
  
  Time<-PLATE[[1]] #time in the first column
  plateData<-PLATE[,-1] #the absorbance data without the time column
  absWells <- length(plateData[1,]) #the no. of columns of the absorbance data
  mint<-min(Time, na.rm = TRUE)
  maxt<-max(Time, na.rm = TRUE)
  maxy<-max(plateData, na.rm = TRUE) #min and max values for scaling the plots
  samples <- colnames(plateData) #names of columns to go on the plots
  
  #par(mfrow=c(8,12))
  par(mfrow=c(NUMR,length(samples)/NUMR))
  par(mar=c(0.2,0.2,0.2,0.2)) # dimensions for figure
  
  lapply(seq_along(plateData), function(i) {
    
    plot(x = Time, y = plateData[, i], type = "l", col="steelblue", lwd=3,
         ylim=c(0, maxy*1.2), xaxt="n", yaxt="n", xlab = "" , ylab="")
    legend("topright", bty="n", paste0(samples[i],"=",i), cex=1.5 )
    
  })
  
}
#Function for  downcurve analysis  

downy<-function(d, Time, ini, thresh, off){
  
  minAbs <- min(d, na.rm = TRUE)+off
  
  #Need to define how to get the max abs
  maxAbs <- max(d, na.rm = TRUE)
  
  pointmax<-which.max(d)
  pcChange<-ini*(maxAbs-minAbs)+minAbs
  #ifelse statements are used to differentiate between clotting or clotlysis curves
  ifelse(d[length(d)]>=pcChange, downTime <- Time, downTime<-Time[-c(1:pointmax)])
  ifelse(d[length(d)]>=pcChange, downAbs <- d,  downAbs<-d[-c(1:pointmax)] )
  #NEED TO LIMIT THIS TO THE MAX -> MIN NOT MAX TO END
  #downTime<-Time[-c(1:pointmax)]
  #downAbs<-d[-c(1:pointmax)] 
  TC <- which.min(downAbs)
  #ifelse(TC==0, TC <- length(downTime), TC <- which.min(downAbs))
  downTime <- downTime[1:TC]
  downAbs <- downAbs[1:TC]
  #decaypoint is where set% lysis occurs
  if_else(d[length(d)]>=pcChange, decayPoint <- length(d), decayPoint<-which(abs(downAbs-pcChange)==min(abs(downAbs-pcChange)))[1] )
  #end point is where 100% lysis occurs
  if_else(d[length(d)]>=pcChange, endPoint <- length(d),endPoint <- which(downAbs<=minAbs)[1] )
  
  if_else(d[length(d)]>=pcChange, endTime <- Time[length(d)], endTime <- downTime[endPoint])
  
  if_else(d[length(d)]>=pcChange, lastPoint <- endPoint,lastPoint <- endPoint+pointmax )
  #will crash if lastPoint is NA
  ifelse(is.na(lastPoint), lastPoint <- length(d), lastPoint <- endPoint+pointmax)
  
  #will this avoid crashes due to endtime = NA
  ifelse(is.na(endTime), endTime <- Time[length(d)], endTime <- endTime)
  
  #In the simple clotlysis app don't use AUC it can be crashy
  #AUC<-sum(diff(Time[1:lastPoint])*(head(d[1:lastPoint],-1)+tail(d[1:lastPoint],-1)))/2
  
  ifelse((max(d)-min(d)<thresh | min(downAbs)>=pcChange),
         decayAbs <- NA,
         decayAbs<-downAbs[decayPoint]#,
         # decayAbs<-round(approx(downTime, downAbs, xout = pcChange, ties = mean)$x,3)
  )
  
  #StartTime is fitted if abs > threshold, otherwise is closest point
  #This prevents crashing if there are blank wells
  ifelse((max(d)-min(d)<thresh | min(downAbs)>=pcChange),
         decayTime <- NA,
         decayTime<-downTime[decayPoint]#,
         # decayTime<-round(approx(downAbs, downTime, xout = decayAbs, ties = mean)$y,3)
  )
  
  downcurve <- c(decayAbs, decayTime, decayPoint+pointmax,lastPoint, endTime, minAbs)
  
}

ui <- fluidPage(
 
  fluidRow( column(4,fileInput("file", "Upload data file (CSV)")),
            column(4,textInput("remove_cols", "Remove column nos (comma-separated):", "-1")),
            column(4, helpText("Removed"),textOutput("remove_txt"))
  ),
  
  fluidRow(
    column(3,numericInput("thresh", "Threshold", "0.02", step=0.01)),
    column(3,numericInput("trunc", "Truncate points", "0", step=10)),
    column(3, numericInput("off", "zero offset", "0", step=0.01)),
    column(3,numericInput("numr", "Plot number of rows", "8"))
  ),
  
  radioButtons(
    inputId = "plotsab",
    label = NULL,
    choices = c("Raw", 
                "Analysed"),
    inline = TRUE
    
  ),
  
  fluidRow(
  
  column(12,plotOutput(outputId = "simpleplot"))
  ),
 
  tableOutput("data_table")
)

# server.R

server <- function(input, output) {
 
   data <- reactive({
    
     req(input$file)
     inputFile <- input$file
     #if (is.null(inputFile)) 
     #  read.csv("./Data/clotlysistrim.csv")
     #else(
    read.csv(input$file$datapath, header = TRUE)
      #   )
  })
  
  remove_cols <- reactive({
    req(input$remove_cols)
    rem <- as.numeric(strsplit(input$remove_cols, ",")[[1]])
    rem+1
  })
  
  col_names <- reactive({
    
    colnames(data())[remove_cols()]
  })
  
  filtered_data <- reactive({
    req(data(), remove_cols())
    data() %>%
    #select(-col_names())
    #mutate_at(vars(remove_cols()), funs(ifelse(is.numeric(.), 0, .)))
    mutate_at(vars(remove_cols()), ~ ifelse(is.numeric(.), 0, .))
  })
  
  output$remove_txt <- renderText({
   
    #paste("Removed:",col_names())#[remove_cols()]
    col_names()
  })
  
  truncated_data <- reactive({
    trunc_plateFun(filtered_data(), input$trunc)
    
  })
  
  zerod_data <- reactive({
    zerod_data <- map_df(truncated_data()[,-1], ~BaselineOff2(.x, input$off) )
    zerod_data <- cbind("Time"=truncated_data()[[1]], zerod_data)
  })
  
  #splot <- reactive({simple_plotFun(zerod_data(), input$numr, TabResdown())})
  #bplot <- reactive({base_plotFun(data(), input$numr)})
  
  output$simpleplot <- renderPlot({
    
                        switch(input$plotsab,
                               "Raw"=base_plotFun(filtered_data(), input$numr),
                               "Analysed"=simple_plotFun(zerod_data(), input$numr, TabResdown())
                                )
                        })
  
  
  TabResdown <- reactive({
    
    args_list <- list(Time = zerod_data()[[1]], ini = 0.5, thresh = input$thresh, off=input$off)
    
    zerod_data()[-1] %>%
      map_df(~ {
        res <- do.call(downy, c(list(.x), args_list))
        data.frame(
          lys.abs   = res[1],
          lys.time  = res[2],
          decayPoint= res[3],
          endPoint  = res[4],
          end.time  = res[5],
          min.Abs   = res[6]#+input$off
          
        )
      })
    
  })
  
  output$data_table <- renderTable({
    
    #zerod_data()
    #TabResdown()
    #myRes <- "lys.time"
    #lysRes <- matrix((TabResdown() %>% select(all_of(myRes)) %>% pull()), nrow = input$numr, byrow = TRUE) 
    #write_clip(lysRes, allow_non_interactive = TRUE)
    #colnames(lysRes) <- 1:(nrow(TabResdown())/input$numr)
    #colnames(lysRes) <- 1:12
    #lysRes
    #TabResdown()
    #data()
    switch(input$plotsab,
           "Raw"=data(),
           "Analysed"=matrix((TabResdown() %>% select(lys.time) %>% pull()), nrow = input$numr, byrow = TRUE)
    )
  })
  
}

shinyApp(ui, server)
