makeSummary <- function(input,refSim,testSim,refPopDf,testPopDf,outputPathToSimulate,wsvParameterPathsList){

  ICVForAUC <- input$ICV_AUC
  ICVForCmax <- input$ICV_Cmax
  if (input$WSVformChoice == 'Parameter-wise within-subject variability'){
    ICVForAUC <- 0
    ICVForCmax <- 0
  }

  paste("*********",
        "\n*Summary*",
        "\n*********",
        "\n",
        "\nSimulation settings",
        "\n=====================",
        "\nOutput quantity path to simulate:", outputPathToSimulate,
        "\nSimulation start time time resolution:", input$simStartTime, input$simTimeUnit,
        "\nSimulation start time time resolution:", input$simEndTime, input$simTimeUnit,
        "\nSimulation time resolution:", input$simResolution,
        "\n",
        "\nBSV and WSV",
        "\n=====================",
        "\nBSV source:", input$BSVformChoice,
        "\nWSV source:", input$WSVformChoice,
        "\nICV for AUC:", ICVForAUC,
        "\nICV for Cmax:", ICVForCmax,
        "\n",
        "\nClinical trial",
        "\n=====================",
        "\nNumber of crossover trial replicates:", input$CT_numberOfReplicates,
        "\nNumber of trials per sample size:", input$CT_n_trials,
        "\nMinimum sample size:", input$CT_subj_min,
        "\nMaximum sample size:", input$CT_subj_max,
        "\nSample size increment:", input$CT_subj_step)
}

nn <- function(input){
  if(is.na(input)){
    return(NULL)
  }
  return(input)
}

runVBEGui <- function(input,refSim,testSim,refPopDf,testPopDf,outputPathToSimulate,wsvParameterPathsList){
  summary <- makeSummary(input,refSim,testSim,refPopDf,testPopDf,outputPathToSimulate,wsvParameterPathsList)
  if(input$BSVformChoice == "Upload PK-Sim population CSV file"){
    demographyRanges <- NULL
    vp <- VirtualPopulation$new(demographyRanges = demographyRanges)
    vp$addPopulationDataFrames(referencePopulationDataframe = refPopDf,testPopulationDataframe = testPopDf)
  } else {
    demographyRanges <- ospsuite::createPopulationCharacteristics(species = nn(input$species),
                                                                  population = nn(input$population),
                                                                  numberOfIndividuals = nn(input$numberOfIndividuals),
                                                                  proportionOfFemales = nn(input$propFemales),
                                                                  weightMin = nn(input$weightMin), weightMax = nn(input$weightMax),weightUnit = nn(input$weightUnit),
                                                                  heightMin = nn(input$heightMin), heightMax = nn(input$heightMax),heightUnit = nn(input$heightUnit),
                                                                  ageMin = nn(input$ageMin),ageMax = nn(input$ageMax),ageUnit = nn(input$ageUnit),
                                                                  BMIMin = nn(input$BMIMin),BMIMax = nn(input$BMIMax),BMIUnit = nn(input$BMIUnit),
                                                                  gestationalAgeMin = nn(input$gestAgeMin),gestationalAgeMax = nn(input$gestAgeMax),gestationalAgeUnit = nn(input$gestAgeUnit),
                                                                  seed = nn(input$seed))
    vp <- VirtualPopulation$new(demographyRanges = demographyRanges)
  }

  if(input$WSVformChoice ==  "Parameter-wise within-subject variability"){
    for (pth in wsvParameterPathsList){
      if(input[[paste0("wsvChoice_",pth)]] == "Multiplicative (CV)"){
        CV <- input[[paste0("CV_",pth)]]
        vp$addMultiplicativeIOV(parameterPath = pth, CV = CV)
      } else {
        CV <- input[[paste0("CV_",pth)]]
        vp$addAdditiveIOV(parameterPath = pth, SD = CV)
      }
    }
    aucICV  = 0
    cmaxICV = 0
  } else {
    aucICV  = input$ICV_AUC
    cmaxICV = input$ICV_Cmax
  }



  ## ---------------------------------------------------------------------------
  ## Virtual Population Simulation
  ## ---------------------------------------------------------------------------
  progressSensitivity <- shiny::Progress$new()
  progressSensitivity$set(message = "Simulating virtual population", value = 0)
  on.exit(progressSensitivity$close())
  refAndTestSimulationsInVirtualPopulation <- vp$simulateVirtualPopulationFromPKML(
    referenceSimulation = refSim,
    testSimulation      = testSim,
    numberOfReplicates  = input$CT_numberOfReplicates,
    outputPath          = outputPathToSimulate,
    startTime           = toBaseUnit(quantityOrDimension = ospDimensions$Time,values = input$simStartTime,unit = input$simTimeUnit),
    endTime             = toBaseUnit(quantityOrDimension = ospDimensions$Time,values = input$simEndTime,unit = input$simTimeUnit),
    resolutionPtsMin    = input$simResolution
  )
  progressSensitivity$set(message = "Simulating virtual population", value = 1)
  pkParametersData <- getAucCmaxDf(refAndTestSimulationsInVirtualPopulation)


  ## ---------------------------------------------------------------------------
  ## Clinical Trial Simulation and VBE Evaluation
  ## ---------------------------------------------------------------------------
  progressCTS <- shiny::Progress$new()
  progressCTS$set(message = "Running clinical trial simulation", value = 0)
  on.exit(progressCTS$close())
  ctsResults <- runClinicalTrialSimulation(
    virtualPopulationSimulationResults = refAndTestSimulationsInVirtualPopulation,
    subj_min   = input$CT_subj_min,
    subj_max   = input$CT_subj_max,
    subj_step  = input$CT_subj_step,
    n_trials   = input$CT_n_trials,
    aucICV     = aucICV,
    cmaxICV    = cmaxICV
  )
  progressCTS$set(message = "Running clinical trial simulation", value = 1)

  # Plot clinical trial simulation results
  bandsPlot <- plotVBEBands(ctsResults)
  show(bandsPlot)

  probabilityPlot <- plotVBEProbability(ctsResults)
  show(probabilityPlot)

  vbeResults <- list(summary = summary, pkParametersData = pkParametersData, ctsResults = ctsResults , bandsPlot = bandsPlot , probabilityPlot = probabilityPlot)

  return(vbeResults)
}

#' @export
runQuickVBE <- function(){

  options(shiny.maxRequestSize = 100*1024^2)

  # =========================
  #           UI
  # =========================
  ui <- fluidPage(

    tags$head(
      tags$title("OSP VBEToolbox Express")
    ),


    h1("OSP VBEToolbox Express", style = "color: #3379b7; font-weight: bold;"),

    tabsetPanel(

      # --- Simulations ---
      tabPanel("Simulations",
               br(),
               h3("Load simulations"),
               fluidRow(
                 column(3, fileInput("refFile", "Reference  simulation")),
                 column(3, fileInput("testFile", "Test simulation"))
               ),
               br(),
               h3("Simulation output quantity to simulate"),
               fluidRow(
                 column(6, uiOutput("selectOutputUI"))
               ),
               fluidRow(
                 column(6, verbatimTextOutput(outputId = "verbatimSelectedOutputPath",placeholder = TRUE))
               ),
               br(),
               h3("Simulation settings"),
               fluidRow(
                 column(3,numericInput(inputId = "simStartTime",label = "Simulation start time",value = 0 , min = 0)),
                 column(3,numericInput(inputId = "simEndTime",label = "Simulation end time",value = 24 , min = 0))),
               fluidRow(
                 column(3,selectInput(inputId = "simTimeUnit",label = "Time unit", selected = "h", choices = ospsuite::getUnitsForDimension(ospDimensions$Time))),
                 column(3,numericInput(inputId = "simResolution",label = "Time resolution (points per minute)", value = 0.25 , min = 0.01,step = 0.01))
               )

      ),

      # --- BSV ---
      tabPanel("Between-subject variability",
               br(),
               radioButtons("BSVformChoice", "Select BSV input method:",
                            choices = c("Upload PK-Sim population CSV file","Input PK-Sim population Characteristics"),
                            inline = TRUE),

               # --- CSV ---
               conditionalPanel(
                 condition = "input.BSVformChoice == 'Upload PK-Sim population CSV file'",
                 fileInput("csvFileRef", "Population CSV for reference simulation", accept = ".csv"),
                 fileInput("csvFileTest", "Population CSV for test simulation *", accept = ".csv"),
                 helpText(HTML("* Use the same population CSV for the reference and test simulations<br/>unless they contain different paths to corresponding parameters.")),
                 verbatimTextOutput("csvRefOut"),
                 verbatimTextOutput("csvTestOut")
               ),

               # --- Manual input ---
               conditionalPanel(
                 condition = "input.BSVformChoice == 'Input PK-Sim population Characteristics'",
                 fluidRow(
                   column(2, selectInput("species", "Species", choices = names(Species),selected = "Human")),
                   column(2, selectInput("population", "Human population", choices = names(HumanPopulation))),
                   column(2, numericInput("numberOfIndividuals", "Number of Individuals", 100, min = 1)),
                   column(2, numericInput("propFemales", "Proportion of Females (%)", 50, min = 0, max = 100))
                 ),
                 fluidRow(
                   column(4, numericInput("weightMin", "Weight Min", NA)),
                   column(4, numericInput("weightMax", "Weight Max", NA)),
                   column(4, selectInput("weightUnit", "Weight Unit", choices = getUnitsForDimension(ospDimensions$Mass)))
                 ),
                 fluidRow(
                   column(4, numericInput("heightMin", "Height Min", NA)),
                   column(4, numericInput("heightMax", "Height Max", NA)),
                   column(4, selectInput("heightUnit", "Height Unit", choices =  getUnitsForDimension(ospDimensions$Length)))
                 ),
                 fluidRow(
                   column(4, numericInput("ageMin", "Age Min", NA)),
                   column(4, numericInput("ageMax", "Age Max", NA)),
                   column(4, selectInput("ageUnit", "Age Unit", choices = getUnitsForDimension(ospDimensions$`Age in years`)))
                 ),
                 fluidRow(
                   column(4, numericInput("BMIMin", "BMI Min", NA)),
                   column(4, numericInput("BMIMax", "BMI Max", NA)),
                   column(4, selectInput("BMIUnit", "BMI Unit", choices = getUnitsForDimension(ospDimensions$BMI)))
                 ),
                 fluidRow(
                   column(4, numericInput("gestAgeMin", "Gestational Age Min", NA)),
                   column(4, numericInput("gestAgeMax", "Gestational Age Max", NA)),
                   column(4, selectInput("gestAgeUnit", "Gestational Age Unit", choices = getUnitsForDimension(ospDimensions$`Age in weeks`)))
                 ),
                 fluidRow(column(4, numericInput("seed", "Seed", NA))),
                 verbatimTextOutput("popOut")
               )


      ),

      # --- WSV ---
      tabPanel("Within-subject variability",
               br(),
               radioButtons("WSVformChoice", "Select WSV input method:",
                            choices = c("Intra-subject coefficients of variation" , "Parameter-wise within-subject variability"),
                            inline = TRUE),

               # CV option
               conditionalPanel(
                 condition = "input.WSVformChoice == 'Intra-subject coefficients of variation'",

                 fluidRow(
                   column(8,
                          numericInput(
                            inputId = paste0("ICV_AUC"),
                            value = 0.147,
                            step = 0.01,
                            min = 0,
                            width = 400,
                            label = "AUC intra-subject coefficients of variation"
                          )
                   )
                 ),

                 fluidRow(
                   column(8,
                          numericInput(
                            inputId = paste0("ICV_Cmax"),
                            value = 0.217,
                            step = 0.01,
                            min = 0,
                            width = 400,
                            label = "Cmax intra-subject coefficients of variation"
                          )
                   )
                 ),

               ),

               # Parameter-wise option
               conditionalPanel(
                 condition = "input.WSVformChoice == 'Parameter-wise within-subject variability'",
                 h4(strong("WSV parameters"), style = "color: #888888;"),
                 actionButton("addParameterButton", "Add parameter",
                              style = "background-color: #3379b7; color: white; font-weight: bold;"),
                 br(), br(),
                 fluidRow(
                   column(6, verbatimTextOutput(outputId = "verbatimSelectedWSVPaths"))
                 ),
                 uiOutput("wsvParameterMenu")
               )

      ),

      # --- Virtual clinical trial ---
      tabPanel("Virtual clinical trial",
               br(),
               column(12,
                      fluidRow(
                        column(2,
                               fluidRow(
                                 numericInput(
                                   inputId = paste0("CT_numberOfReplicates"),
                                   value = 1,
                                   step = 1,
                                   min = 1,
                                   width = "100%",
                                   label = "Number of crossover trial replicates"
                                 )
                               ),
                               hr(),
                               fluidRow(
                                 numericInput(
                                   inputId = paste0("CT_n_trials"),
                                   value = 50,
                                   step = 1,
                                   min = 1,
                                   width = "100%",
                                   label = "Number of trials per sample size"
                                 )
                               ),
                               hr(),
                               fluidRow(
                                 numericInput(
                                   inputId = paste0("CT_subj_min"),
                                   value = 5,
                                   step = 1,
                                   min = 1,
                                   width = "100%",
                                   label = "Minimum sample size"
                                 )
                               ),
                               hr(),
                               fluidRow(
                                 numericInput(
                                   inputId = paste0("CT_subj_max"),
                                   value = 50,
                                   step = 1,
                                   min = 1,
                                   width = "100%",
                                   label = "Maximum sample size"
                                 )
                               ),
                               hr(),
                               fluidRow(
                                 numericInput(
                                   inputId = paste0("CT_subj_step"),
                                   value = 5,
                                   step = 1,
                                   min = 1,
                                   width = "100%",
                                   label = "Sample size increment"
                                 )
                               ),
                               hr(),
                               fluidRow(
                                 downloadButton(
                                   inputId = "runVBE", outputId = "runVBE", label = "Run and download VBE results",
                                   style = "background-color: #3379b7; color: white; font-weight: bold;"
                                 )
                               )
                        )
                      )
               )
      )
    )
  )


  # =========================
  #         SERVER
  # =========================
  server <- function(input, output, session) {



    # -----------------------------
    # Simulations
    # -----------------------------
    refSim <- reactiveValues()
    testSim <- reactiveValues()
    outputPathToSimulate <- reactiveValues()
    dateTime <- reactiveValues()

    # Dynamic UI for the button
    output$selectOutputUI <- renderUI({
      if (is.null(input$refFile)) {
        actionButton("selectOutput", "Select output quantity", disabled = TRUE)
      } else {
        actionButton("selectOutput", "Select output quantity")
      }
    })

    observeEvent(input$refFile,{
      refSim$sim <- ospsuite::loadSimulation(filePath = input$refFile$datapath)
      refSim$outTree <- ospsuite::getSimulationTree(simulationOrFilePath = input$refFile$datapath, quantityType = c("Molecule", "Observer"))
      refSim$parTree <- ospsuite::getSimulationTree(simulationOrFilePath = input$refFile$datapath, quantityType = "Parameter")

    })

    observeEvent(input$testFile,{
      testSim$sim <- ospsuite::loadSimulation(filePath = input$testFile$datapath)
    })

    observeEvent(eventExpr = input$selectOutput,
                 handlerExpr = {
                   print(refSim$sim)
                   output$outputTree <- renderTree({  refSim$outTree })
                   showModal(
                     modalDialog(
                       shinyTree(outputId = "outputTree", checkbox = TRUE),
                       footer = tagList(
                         modalButton("Cancel"),
                         actionButton(inputId = "selectOutputPathButton", label = "Select output path")
                       ),
                       easyClose = TRUE
                     )
                   )
                 }
    )

    observeEvent(eventExpr = input$selectOutputPathButton,
                 {
                   outputPath <- shinyTree::get_selected(tree = input[["outputTree"]], format = "names")
                   outputPath <- sapply(outputPath[outputPath == "path"], function(pp) {
                     ospsuite::toPathString(attr(x = pp, "ancestry"))
                   })
                   print(outputPath)
                   outputPathToSimulate$path <- outputPath[1]
                   output$verbatimSelectedOutputPath <- renderText({
                     outputPathToSimulate$path
                   })

                 }
    )

    # -----------------------------
    # BSV
    # -----------------------------

    observeEvent(input$csvFileRef , {
      req(input$csvFileRef)
      refSim$popDf <- read.csv(input$csvFileRef$datapath, check.names = FALSE)
    })

    observeEvent(input$csvFileTest , {
      req(input$csvFileTest)
      testSim$popDf <- read.csv(input$csvFileTest$datapath, check.names = FALSE)
    })

    # -----------------------------
    # WSV
    # -----------------------------

    wsvParameters <- reactiveValues(pathsList = character(0))

    # Button that opens the modal with the parameter tree
    observeEvent(input$addParameterButton, {
      print("Opening parameter selection modal")
      output$parameterTree <- renderTree({ refSim$parTree })    # ensure the tree exists for the modal
      showModal(
        modalDialog(
          shinyTree(outputId = "parameterTree", checkbox = TRUE),
          footer = tagList(
            modalButton("Cancel"),
            actionButton(inputId = "selectWSVPathsButton", label = "Select paths of parameters with WSV")
          ),
          easyClose = TRUE
        )
      )
    })


    # Handler for the modal button - stores the selected paths and closes the modal
    observeEvent(input$selectWSVPathsButton, {
      # get_selected returns a list-like structure; debug-print it first
      sel <- shinyTree::get_selected(tree = input[["parameterTree"]], format = "names")
      print("raw selection from shinyTree:")
      print(sel)

      # Extract only nodes labeled "path" (same approach as you used before)
      if (length(sel) > 0) {
        chosen <- sel[sel == "path"]
        if (length(chosen) > 0) {
          wsvParameters$pathsList <- sapply(chosen, function(pp) {
            ospsuite::toPathString(attr(x = pp, "ancestry"))
          })
        } else {
          wsvParameters$pathsList <- character(0)
        }
      } else {
        wsvParameters$pathsList <- character(0)
      }

      print("parsed wsvParameters$pathsList:")
      print(wsvParameters$pathsList)

      # update the small text box so user sees what was selected
      output$verbatimSelectedWSVPaths <- renderText({
        if (length(wsvParameters$pathsList) > 0) {
          paste(wsvParameters$pathsList, collapse = "\n")
        } else {
          "No parameter paths selected."
        }
      })

      # close the modal so the main UI (and the new inputs) are visible
      removeModal()
    })


    # Render the dynamic UI once and let it react to wsvParameters$pathsList
    output$wsvParameterMenu <- renderUI({
      n <- length(wsvParameters$pathsList)
      if (n > 0) {
        interaction <- vector("list", n)
        for (i in seq_along(wsvParameters$pathsList)) {
          interaction[[i]] <- tagList(
            fluidRow(
              column(4,
                     tagList(
                       strong("Parameter path"),
                       div(wsvParameters$pathsList[i], style = "margin-top:5px;")
                     )
              ),
              column(4,
                     radioButtons(
                       inputId = paste0("wsvChoice_", wsvParameters$pathsList[i]),
                       label = "Perturbation type",
                       choices = c("Multiplicative (CV)" , "Additive (SD)")
                     )
              ),
              column(4,
                     numericInput(
                       inputId = paste0("CV_", wsvParameters$pathsList[i]),
                       value = 0,
                       step = 0.01,
                       min = 0,
                       width = 200,
                       label = "CV or SD"
                     )
              )
            ),
            hr()
          )
        }
        do.call(tagList, interaction)
      } else {
        tagList(
          p("No WSV parameters selected yet. Click 'Add parameter' to pick parameters from the tree.")
        )
      }
    })



    output$runVBE <- downloadHandler(
      filename = function() {
        dateTime$now <- paste0(format(Sys.Date(), "%Y%m%d"), "_", format(Sys.time(), "%H%M%S"))
        paste("vbe-gui-results-", dateTime$now, ".zip", sep = "")
      },
      content = function(file) {
        results <- runVBEGui(input = input,refSim = refSim$sim, testSim = testSim$sim, refPopDf = refSim$popDf , testPopDf = testSim$popDf , outputPathToSimulate = outputPathToSimulate$path, wsvParameterPathsList = wsvParameters$pathsList)
        fileNames <- NULL

        fname <- paste0("summary", "_" , dateTime$now, ".txt")
        write(x = results$summary, file = fname)
        fileNames <- c(fileNames, fname)

        fname <- paste0("pkParametersData", "_" , dateTime$now, ".csv")
        write.csv(x = results$pkParametersData ,file = fname ,row.names = FALSE)
        fileNames <- c(fileNames, fname)

        # Band plots
        for (pltName in names(results$bandsPlot)) {
          plt <- results$bandsPlot[[pltName]]
          fname <- paste0(pltName, "_", dateTime$now, ".png")
          ggsave(
            filename = fname,
            plot = plt, height = 6, width = 8, units = "in",
            device = "png"
          )
          fileNames <- c(fileNames, fname)
        }

        # Probability plot
        fname <- paste0("probabilityPlot", "_", dateTime$now, ".png")
        ggsave(
          filename = fname,
          plot = results$probabilityPlot, height = 6, width = 8, units = "in",
          device = "png"
        )
        fileNames <- c(fileNames, fname)

        for (trial in names(results$ctsResults$summary)){
          if(is.na(results$ctsResults$summary[trial])){
            next
          }
          fname <- paste0("summary", "_" , trial , "_", dateTime$now, ".csv")
          write.csv(x = results$ctsResults$summary[[trial]],file = fname ,row.names = FALSE)
          fileNames <- c(fileNames, fname)
        }

        zip(file, files = fileNames)
        sapply(fileNames, file.remove)
      }
    )




  }

  shinyApp(ui, server)

}
