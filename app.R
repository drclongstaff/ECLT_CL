library(shiny)

source("./Functions/baseDowny.R")
source("./Functions/smooth.R")
source("./Functions/loadFile.R")
ThisApp <- "Shiny App for ECLT curve analysis"
ThisVer <- "2.1a"

ui <- fluidPage(
  includeCSS("./www/styles3.css"), # make a few changes to the colours and fonts

  tags$h2(ThisApp, ThisVer, align = "center"),
  #tags$a(href = "./docs/ECLT-app-notes.pdf", h5("Help notes")),
  #tags$a(href = "https://drclongstaff.github.io/shiny-clots/docs/ECLT-app-notes.pdf", h5("Help notes")),
  tags$h4("Load a data file, set plotting and fitting parameters", align = "center"),
  fluidRow(
    column(4, fileInput("file", "Upload csv or txt")),
    column(2, radioButtons("zero", "Raw data or zeroed", choices = c("raw", "zeroed"), inline = TRUE)),
    column(2, numericInput("perc", "%clot", value = 50, min = 0, step = 5)),
    column(2, numericInput("crit", "adjust fit", value = 80, min = 20, step = 20)),
    column(2, numericInput("numrow", "plot n rows", value = 8, min = 1))
  ),
  # Set up blank plots
  tags$head(
    tags$style(HTML("

      .plot-grid {
        display: grid;
        gap: 1px;
        width: 100%;
      }

      .plot-cell {
        cursor: pointer;
        border: 1px solid #ddd;

      }
      .plot-cell:hover {
        border: 3px solid #4CAF50;

      }
    "))
  ),
  mainPanel(
    width = 12,
    splitLayout(
      cellWidths = c("50%", "50%"),

      # Left panel - plots
      div(
        tags$h4("Click a well to get a detailed view", align = "center"),
        div(
          style = "overflow-y: auto; max-height: 800px; border: 1px solid #ddd; padding: 5px;",
          uiOutput("dynamicPlotGrid")
        ),
        downloadButton("downloadFig", "Download figure", style = "margin-top: 10px;")
      ),

      # Right panel - table
      div(
        tags$h4("Results table", align = "center"),
        div(
          style = "overflow-y: auto; max-height: 800px;",
          DT::DTOutput("resultsTable")
        ),
        downloadButton("downloadData", "Download table", style = "margin-top: 10px;"),
        helpText(h3(" ")),
        #some blurb and promotional stuff
        helpText(h5("Please cite this reference in publications:")),
        helpText(h5("Longstaff C, ", 
                    tags$a(href="https://doi.org/10.1111/jth.13656","J Thromb Haemost, 15: 1044-6, 2017")
        )),
        
        #tags$i("Please contact me with issues relating to:"),
        helpText(h5("Please contact me",
                    tags$a(href="mailto: drclongstaff@gmail.com", "drclongstaff@gmail.com"), 
                    "for issues relating to:")),
        helpText(h5(ThisApp, ThisVer,
                    " last accessed", Sys.Date()),
        ),
        
        tags$a(href="https://drclongstaff.github.io/shiny-clots/", "Links to other apps and help notes"),
        tags$br(),
        tags$a(href="https://www.youtube.com/@colinlongstaff7270", "Youtube channel of help videos")
      )
    )
  )
)

server <- function(input, output, session) {
  # Load provided data or user data
  myData <- reactive({
    inputFile <- input$file
    if (is.null(input$file)) {
      mD <- read.csv("./data/Copy of Clot Lysis Data.csv")
    } else {
      mD <- load_file(input$file$name, input$file$datapath, input$sheet)
    }
    names(mD)[1] <- "time"
    return(mD)
  })

  # Get the filename
  fileName <- reactive({
    if (is.null(input$file)) {
      return("Copy of Clot Lysis Data.csv")
    } else {
      return(input$file$name)
    }
  })

  # Calculate number of data columns (excluding time column)
  numCols <- reactive({
    ncol(myData()) - 1
  })

  # Calculate number of columns in grid based on rows
  numGridCols <- reactive({
    ceiling(numCols() / input$numrow)
  })
  # Calculate the zeroed data
  dfz <- reactive({
    myD <- myData()
    absz <- data.frame(lapply(myD[, -1], function(x) fun_baseline(x, 0)))
    dfz <- cbind("time" = myD[[1]], absz)
  })

  # Generate the smoothed data from raw or zerod data
  dfs <- reactive({
    myD <- switch(input$zero,
      "raw" = myData(),
      "zeroed" = dfz()
    )
    abss <- data.frame(lapply(myD[, -1], function(x) fun_splsmooth(x, myD, input$crit)))
  })

  # Calculate the lysis times from the smoothed data
  lysTime <- reactive({
    dfzs <- dfs()
    lysis <- data.frame(lapply(dfzs[, -1], function(x) fun_Downy(x, dfzs, input$perc)))
    lysis_n <- signif(lysis[seq(from = 1, to = length(lysis), by = 2)], 4)
  })

  # Dynamically create the plot grid UI
  output$dynamicPlotGrid <- renderUI({
    n <- numCols()
    gridCols <- numGridCols()

    # Update CSS grid columns dynamically
    tags$div(
      tags$style(HTML(sprintf(".plot-grid { grid-template-columns: repeat(%d, 1fr); }", gridCols))),
      div(
        class = "plot-grid",
        lapply(1:n, function(i) {
          div(
            class = "plot-cell",
            plotOutput(paste0("plot_", i),
              height = "80px",
              click = paste0("click_", i)
            )
          )
        })
      )
    )
  })

  # Observe changes and create plots dynamically
  observeEvent(c(myData(), input$zero, input$crit, input$perc, input$numrow), {
    n <- numCols()

    # Create all plots
    lapply(1:n, function(i) {
      output[[paste0("plot_", i)]] <- renderPlot({
        req(i <= numCols()) # Only render if this plot index exists
        req(myData()) # Make sure data is loaded

        myD <- switch(input$zero,
          "raw" = myData(),
          "zeroed" = dfz()
        )
        dfs <- dfs()
        maxy <- max(myD[, -1], na.rm = TRUE)
        samples <- colnames(myD[, -1])

        par(mar = c(0.5, 0.1, 0.5, 0.1))
        plot(myD[[1]], myD[[i + 1]],
          type = "l",
          lwd = 3,
          col = "blue",
          main = paste(samples[i], i),
          xlab = "",
          ylab = "",
          ylim = c(0, maxy * 1.2),
          cex.main = 0.7,
          axes = FALSE
        )

        lines(
          x = dfs[[i * 2 - 1]], y = dfs[[i * 2]],
          type = "l",
          lty = 2,
          lwd = 3,
          col = "red",
          xlab = "",
          ylab = "",
          ylim = c(0, maxy * 1.2)
        )

        abline(v = as.numeric(lysTime()[i]), lty = 2, lwd = 2, col = "olivedrab")
        box()
      })
    })

    # Create click observers for all plots
    lapply(1:n, function(i) {
      observeEvent(input[[paste0("click_", i)]], {
        req(i <= numCols()) # Only respond if this plot index exists

        showModal(modalDialog(
          title = paste("Well", colnames(myData()[, -1])[i], "number", i, "- Detailed View"),
          plotOutput("expandedPlot", height = "500px"),
          size = "l",
          easyClose = TRUE,
          footer = modalButton("Close")
        ))

        output$expandedPlot <- renderPlot({
          req(i <= numCols()) # Only render if this plot index exists

          myD <- switch(input$zero,
            "raw" = myData(),
            "zeroed" = dfz()
          )
          par(mar = c(4, 4, 3, 2))
          plot(myD[[1]], myD[[i + 1]],
            ylab = "Reading",
            xlab = "Time",
            col = "blue",
            cex.main = 1.5,
            cex.lab = 1.2
          )

          grid()

          dfs <- dfs()
          lines(
            x = dfs[[i * 2 - 1]], y = dfs[[i * 2]],
            col = "red",
            lwd = 2,
            lty = 2
          )

          abline(v = as.numeric(lysTime()[i]), lty = 2, col = "olivedrab", lwd = 3)

          legend("topright",
            legend = c("Original Data", "Smoothed Trend"),
            col = c("blue", "red"),
            lwd = 2,
            lty = c(1, 2),
            bty = "n"
          )
        })
      })
    })
  })

  # Make table for results display in ui
  output$resultsTable <- DT::renderDT({
    matrix(as.numeric(lysTime()), nrow = input$numrow, byrow = TRUE)
  })

  # Make table for download
  resultsTable <- reactive({
    matRes <- matrix(as.numeric(lysTime()), nrow = input$numrow, byrow = TRUE)
    # Convert to data frame with row names
    dfRes <- as.data.frame(matRes)
    return(dfRes)
  })

  # Download handler for the results table with metadata
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("results-", format(Sys.time(), "%Y%m%d_%H%M"), ".csv", sep = "")
    },
    content = function(file) {
      # Create metadata header
      metadata <- data.frame(
        Parameter = c("Filename", "%lysis", "Adjust fit"),
        Value = c(fileName(), input$perc, input$crit)
      )

      # Write metadata
      write.table(metadata, file, sep = ",", row.names = FALSE, col.names = TRUE)

      # Add blank line
      write.table("", file, sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)

      # Add header for results
      write.table("Lysis Times:", file, sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)

      # Write results table
      write.table(resultsTable(), file, sep = ",", append = TRUE, row.names = FALSE, col.names = TRUE)
    }
  )

  # Download handler for the plot grid figure
  output$downloadFig <- downloadHandler(
    filename = function() {
      paste("plot-grid-", format(Sys.time(), "%Y%m%d_%H%M"), ".png", sep = "")
    },
    content = function(file) {
      # Calculate grid dimensions
      n <- numCols()
      gridCols <- numGridCols()
      gridRows <- input$numrow

      # Open PNG device with appropriate dimensions
      png(file, width = gridCols * 200, height = gridRows * 200, res = 100)

      # Set up the plot layout
      par(mfrow = c(gridRows, gridCols), mar = c(0.5, 0.1, 0.5, 0.1))

      myD <- switch(input$zero,
        "raw" = myData(),
        "zeroed" = dfz()
      )
      dfs <- dfs()
      maxy <- max(myD[, -1], na.rm = TRUE)
      samples <- colnames(myD[, -1])

      # Create all plots
      for (i in 1:n) {
        plot(myD[[1]], myD[[i + 1]],
          type = "l",
          lwd = 3,
          col = "blue",
          main = paste(samples[i], i),
          xlab = "",
          ylab = "",
          ylim = c(0, maxy * 1.2),
          cex.main = 0.7,
          axes = FALSE
        )

        lines(
          x = dfs[[i * 2 - 1]], y = dfs[[i * 2]],
          type = "l",
          lty = 2,
          lwd = 3,
          col = "red"
        )

        abline(v = as.numeric(lysTime()[i]), lty = 2, lwd = 2, col = "olivedrab")
        box()
      }

      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)
