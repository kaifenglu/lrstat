library(shiny)
library(shinyMatrix)
library(shinyFeedback)
library(shinyjs)
library(dplyr)
library(tidyr)
library(prompter)
library(plotly)
library(lrstat)

# popover texts
followupTimeText = "Follow-up time for the last enrolled subject"

rho1Text = "First parameter of Fleming-Harrington weight"

rho2Text = "Second parameter of Fleming-Harrington weight"

informationRatesText = "Information rates in terms of event fractions"

informationRatesText = gsub(
  "\\s+", " ", "Information rates in terms of variance of
 logrank score statistic under the null hypothesis")

# palette for color-blinded
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# number of time points for power curves
ntpts = 11

# Bracket interval for root finding when the calculation target is
# accrual or follow-up duration
interval = c(0.001, 120)


# reduced style fileInput
fileInputNoExtra<-function(inputId, label, multiple = FALSE, accept = NULL,
                           width = NULL, buttonLabel = "Browse...",
                           placeholder = "No file selected"){

  restoredValue <- restoreInput(id = inputId, default = NULL)
  if (!is.null(restoredValue) && !is.data.frame(restoredValue)) {
    warning("Restored value for ", inputId, " has incorrect format.")
    restoredValue <- NULL
  }
  if (!is.null(restoredValue)) {
    restoredValue <- toJSON(restoredValue, strict_atomic = FALSE)
  }
  inputTag <- tags$input(id = inputId, name = inputId, type = "file",
                         style = "display: none;",
                         `data-restore` = restoredValue)
  if (multiple)
    inputTag$attribs$multiple <- "multiple"
  if (length(accept) > 0)
    inputTag$attribs$accept <- paste(accept, collapse = ",")

  tags$label(
    class = "input-group-btn",
    type="button",
    style=if (!is.null(width))
      paste0("width: ", validateCssUnit(width),";",
             "padding-right: 5px; padding-bottom: 0px; display:inline-block;"),

    span(class = "btn btn-default btn-file",type="button",
         buttonLabel, inputTag,
         style=if (!is.null(width))
           paste0("width: ", validateCssUnit(width),";",
                  "border-radius: 4px; padding-bottom:5px;"))
  )
}


# conditional panels for interim timing
f <- function(i) {
  conditionalPanel(
    condition = paste0("input.kMax == ", i),

    shinyMatrix::matrixInput(
      paste0("xIA_",i),
      label = tags$span(
        "Timing of interim analyses",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = informationRatesText, position = "right")
      ),

      value = matrix(seq_len(i-1)/i,
                     ncol = 1,
                     dimnames = list(paste0("Look ", seq_len(i-1)),
                                     "Information rate")),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)),



    fluidRow(
      column(6, checkboxGroupInput(
        paste0("xES_",i),
        label = "Stop for efficacy",
        choices = paste0("Look ", seq_len(i-1)),
        selected = paste0("Look ", seq_len(i-1)))),

      column(6, checkboxGroupInput(
        paste0("xFS_",i),
        label = "Stop for futility",
        choices = paste0("Look ", seq_len(i-1)),
        selected = paste0("Look ", seq_len(i-1)))),
    ),

  )

}




# initial values for user-specified alpha spending
g <- function(i) {
  conditionalPanel(
    condition = paste0("input.kMax == ", i),

    shinyMatrix::matrixInput(
      paste0("xUA_",i),
      value = matrix(0.0125*(1-exp(4*seq_len(i-1)/i))/(1-exp(4)),
                     ncol = 1,
                     dimnames = list(paste0("Look ", seq_len(i-1)),
                                     c("Cumulative alpha"))),
      inputClass = "numeric",
      rows = list(names=TRUE, extend=FALSE),
      cols = list(names=TRUE, extend=FALSE)
    )
  )
}





# design panel -------------------
designPanel <- tabPanel(
  "Design",

  radioButtons(
    "target",
    "Calculation target",
    choices = c("Power" = "power",
                "Accrual duration" = "accrualDuration",
                "Follow-up duration" = "followupTime"),
    selected = "accrualDuration",
    inline = TRUE
  ),


  fluidRow(
    column(6, numericInput(
      "alpha",
      label = tags$span(
        "alpha",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = "1-sided significance level",
                     position = "right")),
      value = 0.025,
      min = 0.0001, max = 0.4999, step = 0.0001)
    ),

    column(6, numericInput(
      "power",
      "Power",
      value = 0.9,
      min = 0.001, max = 0.999, step = 0.001)
    ),
  ),



  fluidRow(
    column(6, numericInput(
      "accrualDuration",
      "Accrual duration",
      value = 11.643,
      min = 0.01, max = 120, step = 0.001)
    ),

    column(6, numericInput(
      "followupTime",
      label = tags$span(
        "Follow-up duration",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = followupTimeText, position = "right")
      ),
      value = 18,
      min = 0.01, max = 120, step = 0.001)
    ),

  ),




  fluidRow(
    column(6, numericInput(
      "allocationRatioPlanned",
      "Allocation ratio",
      value=1, min=0.1, max=10, step=0.01),
    ),

    column(6, checkboxInput(
      "fixedFollowup",
      "Fixed follow-up", value=FALSE)
    ),

  ),


  fluidRow(
    column(6, numericInput(
      "rho1",
      label = tags$span(
        "rho1",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = rho1Text, position = "right")
      ),
      value = 0,
      min = 0, max = 2, step = 0.01)
    ),

    column(6, numericInput(
      "rho2",
      label = tags$span(
        "rho2",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = rho2Text, position = "right")
      ),
      value = 0,
      min = 0, max = 2, step = 0.01)
    ),
  ),

  checkboxInput(
    "rounding",
    "Rounding up sample size and events", value=TRUE)
)



# test boundaries panel ----------------
boundariesPanel <- tabPanel(
  "Boundaries",

  lapply(2:6, f),

  selectInput(
    "asf", "Alpha spending",
    choices = c(
      "O'Brien-Fleming boundaries" = "OF",
      "Pocock boundaries" = "P",
      "Wang & Tsiatis boundaries" = "WT",
      "O'Brien-Fleming type spending" = "sfOF",
      "Pocock type spending" = "sfP",
      "Kim & DeMets spending" = "sfKD",
      "Hwang, Shi & DeCani spending" = "sfHSD",
      "User defined spending" = "user",
      "No early efficacy stopping" = "none"
    ),
    selected = "sfOF"
  ),

  conditionalPanel(
    condition = "input.asf == 'WT'",
    sliderInput(
      "deltaAlpha",
      "Delta for Wang & Tsiatis efficacy boundaries",
      min=-0.3, max=0.5, value=0.25, step=0.01)
  ),

  conditionalPanel(
    condition = "input.asf == 'sfKD'",
    sliderInput(
      "rhoAlpha",
      "rho for Kim & DeMets alpha spending",
      min=0.4, max=8, value=1, step=0.01)
  ),

  conditionalPanel(
    condition = "input.asf == 'sfHSD'",
    sliderInput(
      "gammaAlpha",
      "gamma for Hwang, Shi & DeCani alpha spending",
      min=-10, max=5, value=1, step=0.01)
  ),

  conditionalPanel(
    condition = "input.asf == 'user'",

    lapply(2:6, g),
  ),



  selectInput(
    "bsf", "Beta spending",
    choices = c(
      "O'Brien-Fleming type spending" = "sfOF",
      "Pocock type spending" = "sfP",
      "Kim & DeMets spending" = "sfKD",
      "Hwang, Shi & DeCani spending" = "sfHSD",
      "No early futility stopping" = "none"
    ),
    selected = "none"
  ),


  conditionalPanel(
    condition = "input.bsf == 'sfKD'",
    sliderInput(
      "rhoBeta",
      "rho for Kim & DeMets beta spending",
      min=0.4, max=8, value=1, step=0.01)
  ),

  conditionalPanel(
    condition = "input.bsf == 'sfHSD'",
    sliderInput(
      "gammaBeta",
      "gamma for Hwang, Shi & DeCani beta spending",
      min=-10, max=5, value=1, step=0.01)
  ),


)



# settings panel -------------------
settingsPanel <- tabPanel(
  "Settings",


  shinyMatrix::matrixInput(
    "survival",
    "Piecewise exponential survival",
    value = matrix(c(0, 0.0309, 0.0533),
                   nrow = 1,
                   dimnames = list(
                     NULL, c("Starting time",
                             "Treatment hazard rate",
                             "Control hazard rate"))
    ),
    inputClass = "numeric",
    rows = list(names=FALSE, extend=FALSE),
    cols = list(names=TRUE, extend=FALSE)
  ),

  actionButton("add_x", label=NULL, icon=icon("plus") ),
  actionButton("del_x", label=NULL, icon=icon("minus")),


  shinyMatrix::matrixInput(
    "accrual",
    "Piecewise constant accrual",
    value = matrix(c(0, 20),
                   nrow=1,
                   dimnames = list(
                     NULL, c("Starting time",
                             "Accrual intensity"))
    ),
    inputClass = "numeric",
    rows = list(names=FALSE, extend=FALSE),
    cols = list(names=TRUE, extend=FALSE)
  ),

  actionButton("add_y", label=NULL, icon=icon("plus") ),
  actionButton("del_y", label=NULL, icon=icon("minus")),


  shinyMatrix::matrixInput(
    "dropout",
    "Exponential dropout",
    value = matrix(c(0, 0),
                   nrow = 1,
                   dimnames = list(
                     NULL, c("Treatment hazard rate",
                             "Control hazard rate"))
    ),
    inputClass = "numeric",
    rows = list(names=FALSE, extend=FALSE),
    cols = list(names=TRUE, extend=FALSE)
  ),

)


# summary panel -------------------
summaryPanel <- tabPanel(
  "Summary",

  htmlOutput("design"),

  htmlOutput("text"),

  conditionalPanel(
    condition = "input.kMax > 1",
    verbatimTextOutput("table"),
  ),


  htmlOutput("text0"),

  conditionalPanel(
    condition = "input.kMax > 1",
    verbatimTextOutput("table0"),
  ),
)


# plot panel -------------------
plotPanel <- tabPanel(
  "Plot",

  conditionalPanel(
    condition = "input.kMax > 1",

    selectInput(
      "plottype",
      "Plot type",
      choices = c("Boundaries (Z)" = "boundaryZ",
                  "Boundaries (HR)" = "boundaryHR",
                  "Boundaries (p)" = "boundaryP",
                  "Error Spending" = "errorSpend",
                  "Sample Size / Events vs. Time" = "eventPred",
                  "Power vs. Follow-up Duration" = "powerVsTf",
                  "Power vs. Sample Size" = "powerVsN",
                  "Study Duration vs. Sample Size" = "TsVsN"),
      selected = "boundaryZ"
    ),


  ),


  conditionalPanel(
    condition = "input.kMax == 1",

    selectInput(
      "plottype2",
      "Plot type",
      choices = c("Sample Size / Events vs. Time" = "eventPred",
                  "Power vs. Follow-up Duration" = "powerVsTf",
                  "Power vs. Sample Size" = "powerVsN",
                  "Study Duration vs. Sample Size" = "TsVsN"),
      selected = "boundaryZ"
    ),




  ),

  plotlyOutput("plot"),
)


# simulation panel -------------------
simulationPanel <- tabPanel(
  "Simulation",

  shinyMatrix::matrixInput(
    "boundaries",
    label = tags$span(
      "Analysis timing and boundaries",
      tags$span(icon(name = "question-circle")) %>%
        add_prompt(message = informationRatesText, position = "right")
    ),
    value = matrix(),
    inputClass = "numeric",
    rows = list(names = TRUE, extend = FALSE),
    cols = list(names = TRUE, extend = FALSE)
  ),


  fluidRow(
    column(3, numericInput(
      "nIterations",
      label = "# simulations",
      value = 1000,
      min = 200, max = 10000, step = 1)
    ),

    column(3, numericInput(
      "nRawDatasets",
      label = tags$span(
        "Raw datasets",
        tags$span(icon(name = "question-circle")) %>%
          add_prompt(message = "# raw datasets to extract per stage",
                     position = "right")),
      value = 1,
      min = 0, max = 100, step = 1)
    ),

    column(3, numericInput(
      "seed",
      label = "Seed",
      value = 100,
      min = 0, max = 100000, step = 1
    )),

    column(3, style = "margin-top: 25px;", actionButton(
      "sim", "Simulate",
      style="color: #fff; background-color: #337ab7; border-color: #2e6da4"
    )),

  ),


  htmlOutput("simtext"),

  conditionalPanel(
    condition = "input.kMax > 1",

    verbatimTextOutput("simtable"),
  ),


  conditionalPanel(
    condition = "input.sim > 0",

    fluidRow(
      column(6, downloadButton(
        "downloadSumdata", "Download summary data under H1")),

      conditionalPanel(
        condition = "input.nRawDatasets > 0",

        column(6, downloadButton(
          "downloadRawdata", "Download raw data under H1")),
      ),
    ),
  ),

  htmlOutput("simtext0"),

  conditionalPanel(
    condition = "input.kMax > 1",

    verbatimTextOutput("simtable0"),
  ),

  conditionalPanel(
    condition = "input.sim > 0",

    fluidRow(
      column(6, downloadButton(
        "downloadSumdata0", "Download summary data under H0")),

      conditionalPanel(
        condition = "input.nRawDatasets > 0",

        column(6, downloadButton(
          "downloadRawdata0", "Download raw data under H0")),
      ),
    ),
  ),

)


# code panel -------------------
codePanel <- tabPanel(
  "Code",

  helpText('library(lrstat)'),

  fluidRow(
    column(6, htmlOutput("lrp")),
    column(6, htmlOutput("lrs")),
  ),
)


# user interface ----------------
ui <- fluidPage(

  shinyFeedback::useShinyFeedback(),
  shinyjs::useShinyjs(),
  prompter::use_prompt(),

  titlePanel(tagList(
    span(HTML(paste(tags$span(style="font-size:14pt",
                              paste("Power and Sample Size Calculation",
                                    "for Non-Proportional Hazards")))),
         span(downloadButton("saveInputs", "Save inputs"),
              fileInputNoExtra("loadInputs", label=NULL, accept=".RData",
                               buttonLabel=list(icon("upload"),
                                                "Load inputs"),
                               width="116px"),
              style="position:absolute;right:0.5em;",
              tags$style(type='text/css', "#saveInputs{margin-top: -5px;}")
         ))),
    windowTitle = paste("Power and Sample Size Calculation",
                        "for Non-Proportional Hazards")
  ),


  sidebarLayout(
    sidebarPanel(

      fluidRow(
        column(7, selectInput(
          "kMax", "Number of stages",
          choices = seq_len(6), selected = 2),
        ),

        column(5,  style = "margin-top: 25px;", actionButton(
          "calc", "Calculate",
          style="color: #fff; background-color: #337ab7; border-color: #2e6da4"
        ),
        ),

      ),



      tabsetPanel(
        id = "Home",
        designPanel,
        boundariesPanel,
        settingsPanel,
      ),



    ),


    mainPanel(

      tabsetPanel(
        id = "Results",
        summaryPanel,
        plotPanel,
        simulationPanel,
        codePanel,
      ),
    )

  )
)




# server function -------------
server <- function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })

  kMax <- reactive(as.numeric(input$kMax))

  observeEvent(input$kMax, {
    if (input$kMax > 1) {
      showTab(inputId = "Home", target = "Boundaries")
    } else {
      hideTab(inputId = "Home", target = "Boundaries")
    }
  })

  observeEvent(input$target, {
    shinyjs::toggleState("power",
                         input$target != "power")
    shinyjs::toggleState("accrualDuration",
                         input$target != "accrualDuration")
    shinyjs::toggleState("followupTime",
                         input$target != "followupTime")
  })


  observe({
    shinyjs::toggleState("saveInputs", input$calc > 0)

    shinyjs::toggle(paste0("xES_", kMax()),
                    condition = {input$asf != 'none'})

    shinyjs::toggle(paste0("xFS_", kMax()),
                    condition = {input$bsf != 'none'})
  })


  # edit check for alpha
  alpha <- reactive({
    req(input$alpha)
    valid <- (input$alpha >= 0.00001) && (input$alpha < 0.5)
    shinyFeedback::feedbackWarning(
      "alpha", !valid,
      "alpha is out of bounds [0.00001, 0.5)")
    req(valid)
    as.numeric(input$alpha)
  })


  # edit check for power
  beta <- reactive({
    if (input$target != "power") {
      req(input$power)
      valid <- (input$power <= 0.9999) &&
        (input$power > input$alpha)
      shinyFeedback::feedbackWarning(
        "power", !valid,
        "Power is out of bounds (alpha, 0.9999]")
      req(valid)
      1 - as.numeric(input$power)
    }
  })




  # edit check for accrual duration
  accrualDuration <- reactive({
    if (input$target != "accrualDuration") {
      req(input$accrualDuration)
      valid <- (input$accrualDuration > 0)
      shinyFeedback::feedbackWarning(
        "accrualDuration", !valid,
        "Accrual duration must be positive")
      req(valid)
      as.numeric(input$accrualDuration)
    }
  })


  # edit check for follow-up duration
  followupTime <- reactive({
    if (input$target != 'followupTime') {
      req(input$followupTime)
      if (input$fixedFollowup) {
        valid <- (input$followupTime > 0)
        shinyFeedback::feedbackWarning(
          "followupTime", !valid,
          "Follow-up duration must be positive for fixed follow-up")
      } else {
        valid <- (input$followupTime >= 0)
        shinyFeedback::feedbackWarning(
          "followupTime", !valid,
          "Follow-up duration must be nonnegative for variable follow-up")
      }
      req(valid)
      as.numeric(input$followupTime)
    }
  })


  # edit check for allocation ratio
  allocationRatioPlanned <- reactive({
    req(input$allocationRatioPlanned)
    valid <- (input$allocationRatioPlanned > 0)
    shinyFeedback::feedbackWarning(
      "allocationRatioPlanned", !valid,
      "Allocation ratio must be positive")
    req(valid)
    as.numeric(input$allocationRatioPlanned)
  })


  # edit check for rho1
  rho1 <- reactive({
    req(input$rho1)
    valid <- (input$rho1 >= 0)
    shinyFeedback::feedbackWarning(
      "rho1", !valid,
      "rho1 must be nonnegative")
    req(valid)
    as.numeric(input$rho1)
  })


  # edit check for rho2
  rho2 <- reactive({
    req(input$rho2)
    valid <- (input$rho2 >= 0)
    shinyFeedback::feedbackWarning(
      "rho2", !valid,
      "rho2 must be nonnegative")
    req(valid)
    as.numeric(input$rho2)
  })


  # edit check for interim timing
  informationRates <- eventReactive(input$calc, {
    if (kMax()>1) {
      d = input[[paste0("xIA_", kMax())]]
      t = c(0, as.vector(d[,1], "numeric"), 1)
      valid = all(diff(t) > 0)
      if (!valid) {
        showNotification(
          paste0("Information rates must be strictly ",
                 "increasing and lie between 0 and 1")
        )
      }
      req(valid)
      t[-1]
    } else {
      1
    }
  })


  efficacyStopping <- eventReactive(input$calc, {
    if (kMax()>1) {
      d = input[[paste0("xES_", kMax())]]
      checked = c(d, paste0("Look ", kMax()))
      t = rep(0, kMax())
      for (i in 1:kMax()) {
        if (paste0("Look ", i) %in% checked) t[i] = 1
      }
      t
    } else {
      1
    }
  })


  futilityStopping <- eventReactive(input$calc, {
    if (kMax()>1) {
      d = input[[paste0("xFS_", kMax())]]
      checked = c(d, paste0("Look ", kMax()))
      t = rep(0, kMax())
      for (i in 1:kMax()) {
        if (paste0("Look ", i) %in% checked) t[i] = 1
      }
      t
    } else {
      1
    }
  })


  # alpha spending
  typeAlphaSpending <- reactive({
    if (kMax()>1) {
      input$asf
    } else {
      "none"
    }
  })


  parameterAlphaSpending <- reactive({
    if (kMax()>1) {
      if (input$asf == "WT") {
        as.numeric(input$deltaAlpha)
      } else if (input$asf == "sfKD") {
        as.numeric(input$rhoAlpha)
      } else if (input$asf == "sfHSD") {
        as.numeric(input$gammaAlpha)
      } else {
        NA_real_
      }
    } else {
      NA_real_
    }
  })


  # edit check for user-specified alpha spending
  userAlphaSpending = eventReactive(input$calc, {
    if (kMax() > 1) {
      d = input[[paste0("xUA_", kMax())]]
      t = c(0, as.vector(d[,1], "numeric"), alpha())
      valid = all(diff(t) >= 0)
      if (!valid) {
        showNotification(
          "Cumulative alpha must be non-decreasing in [0, alpha]"
        )
      }
      req(valid)
      t[-1]
    } else {
      alpha()
    }
  })


  # beta-spending
  typeBetaSpending <- reactive({
    if (kMax()>1) {
      input$bsf
    } else {
      "none"
    }
  })


  parameterBetaSpending <- reactive({
    if (kMax()>1) {
      if (input$bsf == "sfKD") {
        as.numeric(input$rhoBeta)
      } else if (input$bsf == "sfHSD") {
        as.numeric(input$gammaBeta)
      } else {
        NA_real_
      }
    } else {
      NA_real_
    }
  })




  # survival information

  observeEvent(input$add_x, {
    a = matrix(as.numeric(input$survival),
               ncol=ncol(input$survival))
    b = matrix(a[nrow(a),], nrow=1)
    b[,1] = b[,1] + 1
    c = rbind(a, b)
    colnames(c) = colnames(input$survival)
    updateMatrixInput(session, "survival", c)
  })



  observeEvent(input$del_x, {
    if (nrow(input$survival) >= 2) {
      a = matrix(as.numeric(input$survival),
                 ncol=ncol(input$survival))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      colnames(b) = colnames(input$survival)
      updateMatrixInput(session, "survival", b)
    }
  })


  piecewiseSurvivalTime <- reactive({
    t = as.vector(input$survival[,1], "numeric")
    valid = all(diff(t) > 0) && (t[1]==0)
    if (!valid) {
      showNotification(
        "Starting time should be increasing and start at zero"
      )
    }
    req(valid)
    t
  })



  lambda1 <- reactive({
    lam1 = as.vector(input$survival[,2], "numeric")
    valid = all(lam1 >= 0)
    if (!valid) {
      showNotification(
        "Treatment hazard rate must be non-negative")
    }
    req(valid)
    lam1
  })



  lambda2 <- reactive({
    lam2 = as.vector(input$survival[,3], "numeric")
    valid = all(lam2 >= 0)
    if (!valid) {
      showNotification(
        "Control hazard rate must be non-negative")
    }
    req(valid)
    lam2
  })


  # edit check for dropout
  gamma1 <- reactive({
    gam1 = as.numeric(input$dropout[,1])
    valid = (gam1 >= 0)
    if (!valid) {
      showNotification(
        "Treatment hazard rate must be non-negative")
    }
    req(valid)
    gam1
  })


  gamma2 <- reactive({
    gam2 = as.numeric(input$dropout[,2])
    valid = (gam2 >= 0)
    if (!valid) {
      showNotification(
        "Control hazard rate must be non-negative")
    }
    req(valid)
    gam2
  })



  # accrual information
  observeEvent(input$add_y, {
    a = matrix(as.numeric(input$accrual),
               ncol=ncol(input$accrual))
    b = matrix(a[nrow(a),], nrow=1)
    b[,1] = b[,1] + 1
    c = rbind(a, b)
    colnames(c) = colnames(input$accrual)
    updateMatrixInput(session, "accrual", c)
  })


  observeEvent(input$del_y, {
    if (nrow(input$accrual) >= 2) {
      a = matrix(as.numeric(input$accrual),
                 ncol=ncol(input$accrual))
      b = matrix(a[-nrow(a),], ncol=ncol(a))
      colnames(b) = colnames(input$accrual)
      updateMatrixInput(session, "accrual", b)
    }
  })


  accrualTime <- reactive({
    t = as.vector(input$accrual[,1], "numeric")
    valid = all(diff(t) > 0) && (t[1]==0)
    if (!valid) {
      showNotification(
        "Starting time should be increasing and start at zero"
      )
    }
    req(valid)
    t
  })


  accrualIntensity <- reactive({
    a = as.vector(input$accrual[,2], "numeric")
    valid = all(a >= 0)
    if (!valid) {
      showNotification(
        "Accrual intensity must be non-negative")
    }
    req(valid)
    a
  })


  typeOfComputation <- reactive({
    hr = lambda1()/lambda2()
    ph = all(abs(hr - hr[1]) <= 1e-8)
    ifelse(rho1() == 0 && rho2() == 0 && ph, "schoenfeld", "direct")
  })


  lr <- eventReactive(input$calc, {
    if (input$target == "power") {
      tryCatch({
        l = lrpower(
          kMax = kMax(),
          informationRates = informationRates(),
          efficacyStopping = efficacyStopping(),
          futilityStopping = futilityStopping(),
          alpha = alpha(),
          typeAlphaSpending = typeAlphaSpending(),
          parameterAlphaSpending = parameterAlphaSpending(),
          userAlphaSpending = userAlphaSpending(),
          typeBetaSpending = typeBetaSpending(),
          parameterBetaSpending = parameterBetaSpending(),
          allocationRatioPlanned = allocationRatioPlanned(),
          accrualTime = accrualTime(),
          accrualIntensity = accrualIntensity(),
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          lambda1 = lambda1(),
          lambda2 = lambda2(),
          gamma1 = gamma1(),
          gamma2 = gamma2(),
          accrualDuration = accrualDuration(),
          followupTime = followupTime(),
          fixedFollowup = input$fixedFollowup,
          rho1 = rho1(),
          rho2 = rho2(),
          typeOfComputation = typeOfComputation())
      }, error = function(e) {
        shiny:::reactiveStop(conditionMessage(e))
      })

      updateNumericInput(
        session, "power",
        value = round(l$overallResults$overallReject, 3)
      )
    } else if (input$target == "accrualDuration") {
      tryCatch({
        l = lrsamplesize(
          beta = beta(),
          kMax = kMax(),
          informationRates = informationRates(),
          efficacyStopping = efficacyStopping(),
          futilityStopping = futilityStopping(),
          alpha = alpha(),
          typeAlphaSpending = typeAlphaSpending(),
          parameterAlphaSpending = parameterAlphaSpending(),
          userAlphaSpending = userAlphaSpending(),
          typeBetaSpending = typeBetaSpending(),
          parameterBetaSpending = parameterBetaSpending(),
          allocationRatioPlanned = allocationRatioPlanned(),
          accrualTime = accrualTime(),
          accrualIntensity = accrualIntensity(),
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          lambda1 = lambda1(),
          lambda2 = lambda2(),
          gamma1 = gamma1(),
          gamma2 = gamma2(),
          accrualDuration = NA,
          followupTime = followupTime(),
          fixedFollowup = input$fixedFollowup,
          rho1 = rho1(),
          rho2 = rho2(),
          typeOfComputation = typeOfComputation(),
          interval = interval,
          rounding = input$rounding)
      }, error = function(e) {
        shiny:::reactiveStop(conditionMessage(e))
      })

      updateNumericInput(
        session, "accrualDuration",
        value = round(l$resultsUnderH1$overallResults$accrualDuration, 3)
      )
    } else if (input$target == "followupTime") {
      tryCatch({
        l = lrsamplesize(
          beta = beta(),
          kMax = kMax(),
          informationRates = informationRates(),
          efficacyStopping = efficacyStopping(),
          futilityStopping = futilityStopping(),
          alpha = alpha(),
          typeAlphaSpending = typeAlphaSpending(),
          parameterAlphaSpending = parameterAlphaSpending(),
          userAlphaSpending = userAlphaSpending(),
          typeBetaSpending = typeBetaSpending(),
          parameterBetaSpending = parameterBetaSpending(),
          allocationRatioPlanned = allocationRatioPlanned(),
          accrualTime = accrualTime(),
          accrualIntensity = accrualIntensity(),
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          lambda1 = lambda1(),
          lambda2 = lambda2(),
          gamma1 = gamma1(),
          gamma2 = gamma2(),
          accrualDuration = accrualDuration(),
          followupTime = NA,
          fixedFollowup = input$fixedFollowup,
          rho1 = rho1(),
          rho2 = rho2(),
          typeOfComputation = typeOfComputation(),
          interval = interval,
          rounding = input$rounding)
      }, error = function(e) {
        shiny:::reactiveStop(conditionMessage(e))
      })

      updateNumericInput(
        session, "followupTime",
        value = round(l$resultsUnderH1$overallResults$followupTime, 3)
      )
    }

    l
  })


  # power vs. Tf

  tmin1 <- eventReactive(input$calc, {
    a <- lr()$resultsUnderH1$overallResults
    if (rho1() == 0 && rho2() == 0) {
      criticalValues = lr()$resultsUnderH1$byStageResults$efficacyBounds
    } else {
      criticalValues = NA
    }

    if (a$fixedFollowup) {
      tmin <- lrsamplesize(
        beta = 0.8,
        kMax = kMax(),
        informationRates = informationRates(),
        efficacyStopping = efficacyStopping(),
        futilityStopping = futilityStopping(),
        criticalValues = criticalValues,
        alpha = alpha(),
        typeAlphaSpending = typeAlphaSpending(),
        parameterAlphaSpending = parameterAlphaSpending(),
        userAlphaSpending = userAlphaSpending(),
        typeBetaSpending = typeBetaSpending(),
        parameterBetaSpending = parameterBetaSpending(),
        allocationRatioPlanned = allocationRatioPlanned(),
        accrualTime = accrualTime(),
        accrualIntensity = accrualIntensity(),
        piecewiseSurvivalTime = piecewiseSurvivalTime(),
        lambda1 = lambda1(),
        lambda2 = lambda2(),
        gamma1 = gamma1(),
        gamma2 = gamma2(),
        accrualDuration = a$accrualDuration,
        followupTime = NA,
        fixedFollowup = a$fixedFollowup,
        rho1 = a$rho1,
        rho2 = a$rho2,
        estimateHazardRatio = 0,
        typeOfComputation = typeOfComputation(),
        interval = interval,
        rounding = input$rounding
      )$resultsUnderH1$overallResults$followupTime
    } else {
      tmin = 0
    }
    tmin
  })


  powerVsTf <- eventReactive(input$calc, {
    if (kMax()>1) {
      a = lr()$resultsUnderH1$overallResults
      if (rho1() == 0 && rho2() == 0) {
        criticalValues = lr()$resultsUnderH1$byStageResults$efficacyBounds
      } else {
        criticalValues = NA
      }

      time = seq(tmin1(), a$followupTime, length.out=ntpts)
      power = rep(0, ntpts)
      for (i in 1:ntpts) {
        power[i] <- lrpower(
          kMax = kMax(),
          informationRates = informationRates(),
          efficacyStopping = efficacyStopping(),
          futilityStopping = futilityStopping(),
          criticalValues = criticalValues,
          alpha = alpha(),
          typeAlphaSpending = typeAlphaSpending(),
          parameterAlphaSpending = parameterAlphaSpending(),
          userAlphaSpending = userAlphaSpending(),
          typeBetaSpending = typeBetaSpending(),
          parameterBetaSpending = parameterBetaSpending(),
          allocationRatioPlanned = allocationRatioPlanned(),
          accrualTime = accrualTime(),
          accrualIntensity = accrualIntensity(),
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          lambda1 = lambda1(),
          lambda2 = lambda2(),
          gamma1 = gamma1(),
          gamma2 = gamma2(),
          accrualDuration = a$accrualDuration,
          followupTime = time[i],
          fixedFollowup = a$fixedFollowup,
          rho1 = a$rho1,
          rho2 = a$rho2,
          estimateHazardRatio = 0,
          typeOfComputation = typeOfComputation()
        )$overallResults$overallReject
      }

      data.frame(time, power)
    }
  })


  powerVsTf2 <- eventReactive(input$calc, {
    if (kMax()==1) {
      a = lr()$resultsUnderH1$overallResults
      time = seq(tmin1(), a$followupTime, length.out=ntpts)
      power = rep(0, ntpts)
      for (i in 1:ntpts) {
        power[i] <- lrpower(
          kMax = kMax(),
          alpha = alpha(),
          allocationRatioPlanned = allocationRatioPlanned(),
          accrualTime = accrualTime(),
          accrualIntensity = accrualIntensity(),
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          lambda1 = lambda1(),
          lambda2 = lambda2(),
          gamma1 = gamma1(),
          gamma2 = gamma2(),
          accrualDuration = a$accrualDuration,
          followupTime = time[i],
          fixedFollowup = a$fixedFollowup,
          rho1 = a$rho1,
          rho2 = a$rho2,
          estimateHazardRatio = 0,
          typeOfComputation = typeOfComputation()
        )$overallResults$overallReject
      }

      data.frame(time, power)
    }
  })


  # power vs. N

  tmin2 <- eventReactive(input$calc, {
    a <- lr()$resultsUnderH1$overallResults
    if (rho1() == 0 && rho2() == 0) {
      criticalValues = lr()$resultsUnderH1$byStageResults$efficacyBounds
    } else {
      criticalValues = NA
    }

    lrsamplesize(
      beta = 0.8,
      kMax = kMax(),
      informationRates = informationRates(),
      efficacyStopping = efficacyStopping(),
      futilityStopping = futilityStopping(),
      criticalValues = criticalValues,
      alpha = alpha(),
      typeAlphaSpending = typeAlphaSpending(),
      parameterAlphaSpending = parameterAlphaSpending(),
      userAlphaSpending = userAlphaSpending(),
      typeBetaSpending = typeBetaSpending(),
      parameterBetaSpending = parameterBetaSpending(),
      allocationRatioPlanned = allocationRatioPlanned(),
      accrualTime = accrualTime(),
      accrualIntensity = accrualIntensity(),
      piecewiseSurvivalTime = piecewiseSurvivalTime(),
      lambda1 = lambda1(),
      lambda2 = lambda2(),
      gamma1 = gamma1(),
      gamma2 = gamma2(),
      accrualDuration = NA,
      followupTime = a$followupTime,
      fixedFollowup = a$fixedFollowup,
      rho1 = a$rho1,
      rho2 = a$rho2,
      estimateHazardRatio = 0,
      typeOfComputation = typeOfComputation(),
      interval = interval,
      rounding = input$rounding
    )$resultsUnderH1$overallResults$accrualDuration
  })



  powerVsN <- eventReactive(input$calc, {
    if (kMax()>1) {
      a = lr()$resultsUnderH1$overallResults
      if (rho1() == 0 && rho2() == 0) {
        criticalValues = lr()$resultsUnderH1$byStageResults$efficacyBounds
      } else {
        criticalValues = NA
      }

      time = seq(tmin2(), a$accrualDuration, length.out=ntpts)
      N = rep(0, ntpts)
      power = rep(0, ntpts)
      for (i in 1:ntpts) {
        lrx <- lrpower(
          kMax = kMax(),
          informationRates = informationRates(),
          efficacyStopping = efficacyStopping(),
          futilityStopping = futilityStopping(),
          criticalValues = criticalValues,
          alpha = alpha(),
          typeAlphaSpending = typeAlphaSpending(),
          parameterAlphaSpending = parameterAlphaSpending(),
          userAlphaSpending = userAlphaSpending(),
          typeBetaSpending = typeBetaSpending(),
          parameterBetaSpending = parameterBetaSpending(),
          allocationRatioPlanned = allocationRatioPlanned(),
          accrualTime = accrualTime(),
          accrualIntensity = accrualIntensity(),
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          lambda1 = lambda1(),
          lambda2 = lambda2(),
          gamma1 = gamma1(),
          gamma2 = gamma2(),
          accrualDuration = time[i],
          followupTime = a$followupTime,
          fixedFollowup = a$fixedFollowup,
          rho1 = a$rho1,
          rho2 = a$rho2,
          estimateHazardRatio = 0,
          typeOfComputation = typeOfComputation())

        N[i] = lrx$overallResults$numberOfSubjects
        power[i] = lrx$overallResults$overallReject
      }

      data.frame(N, power)
    }
  })


  powerVsN2 <- eventReactive(input$calc, {
    if (kMax()==1) {
      a = lr()$resultsUnderH1$overallResults
      time = seq(tmin2(), a$accrualDuration, length.out=ntpts)
      N = rep(0, ntpts)
      power = rep(0, ntpts)
      for (i in 1:ntpts) {
        lrx <- lrpower(
          kMax = kMax(),
          alpha = alpha(),
          allocationRatioPlanned = allocationRatioPlanned(),
          accrualTime = accrualTime(),
          accrualIntensity = accrualIntensity(),
          piecewiseSurvivalTime = piecewiseSurvivalTime(),
          lambda1 = lambda1(),
          lambda2 = lambda2(),
          gamma1 = gamma1(),
          gamma2 = gamma2(),
          accrualDuration = time[i],
          followupTime = a$followupTime,
          fixedFollowup = a$fixedFollowup,
          rho1 = a$rho1,
          rho2 = a$rho2,
          estimateHazardRatio = 0,
          typeOfComputation = typeOfComputation()
        )

        N[i] = lrx$overallResults$numberOfSubjects
        power[i] = lrx$overallResults$overallReject
      }

      data.frame(N, power)
    }
  })


  # Ts vs. N

  tmin3 <- eventReactive(input$calc, {
    a <- lr()$resultsUnderH1$overallResults
    if (rho1() == 0 && rho2() == 0) {
      criticalValues = lr()$resultsUnderH1$byStageResults$efficacyBounds
    } else {
      criticalValues = NA
    }

    lrsamplesize(
      beta = 1-a$overallReject,
      kMax = kMax(),
      informationRates = informationRates(),
      efficacyStopping = efficacyStopping(),
      futilityStopping = futilityStopping(),
      criticalValues = criticalValues,
      alpha = alpha(),
      typeAlphaSpending = typeAlphaSpending(),
      parameterAlphaSpending = parameterAlphaSpending(),
      userAlphaSpending = userAlphaSpending(),
      typeBetaSpending = typeBetaSpending(),
      parameterBetaSpending = parameterBetaSpending(),
      allocationRatioPlanned = allocationRatioPlanned(),
      accrualTime = accrualTime(),
      accrualIntensity = accrualIntensity(),
      piecewiseSurvivalTime = piecewiseSurvivalTime(),
      lambda1 = lambda1(),
      lambda2 = lambda2(),
      gamma1 = gamma1(),
      gamma2 = gamma2(),
      accrualDuration = NA,
      followupTime = 60,
      fixedFollowup = a$fixedFollowup,
      rho1 = a$rho1,
      rho2 = a$rho2,
      estimateHazardRatio = 0,
      typeOfComputation = typeOfComputation(),
      interval = interval,
      rounding = input$rounding
    )$resultsUnderH1$overallResults$accrualDuration
  })


  TsVsN <- eventReactive(input$calc, {
    if (kMax()>1) {
      a = lr()$resultsUnderH1$overallResults
      time = seq(tmin3(), a$accrualDuration, length.out=ntpts)
      N = rep(0, ntpts)
      Ts = rep(0, ntpts)
      for (i in 1:ntpts) {
        if (typeOfComputation() == "schoenfeld") {
          N[i] = accrual(
            time = time[i],
            accrualTime = accrualTime(),
            accrualIntensity = accrualIntensity(),
            accrualDuration = time[i])

          Ts[i] = caltime(
            nevents = a$numberOfEvents,
            allocationRatioPlanned = allocationRatioPlanned(),
            accrualTime = accrualTime(),
            accrualIntensity = accrualIntensity(),
            piecewiseSurvivalTime = piecewiseSurvivalTime(),
            lambda1 = lambda1(),
            lambda2 = lambda2(),
            gamma1 = gamma1(),
            gamma2 = gamma2(),
            accrualDuration = time[i],
            followupTime = 1200,
            fixedFollowup = a$fixedFollowup)
        } else {
          lrx <- lrsamplesize(
            beta = 1-a$overallReject,
            kMax = kMax(),
            informationRates = informationRates(),
            efficacyStopping = efficacyStopping(),
            futilityStopping = futilityStopping(),
            alpha = alpha(),
            typeAlphaSpending = typeAlphaSpending(),
            parameterAlphaSpending = parameterAlphaSpending(),
            userAlphaSpending = userAlphaSpending(),
            typeBetaSpending = typeBetaSpending(),
            parameterBetaSpending = parameterBetaSpending(),
            allocationRatioPlanned = allocationRatioPlanned(),
            accrualTime = accrualTime(),
            accrualIntensity = accrualIntensity(),
            piecewiseSurvivalTime = piecewiseSurvivalTime(),
            lambda1 = lambda1(),
            lambda2 = lambda2(),
            gamma1 = gamma1(),
            gamma2 = gamma2(),
            accrualDuration = time[i],
            followupTime = NA,
            fixedFollowup = a$fixedFollowup,
            rho1 = a$rho1,
            rho2 = a$rho2,
            estimateHazardRatio = 0,
            typeOfComputation = typeOfComputation(),
            interval = interval,
            rounding = input$rounding)$resultsUnderH1

          N[i] = lrx$overallResults$numberOfSubjects
          Ts[i] = lrx$overallResults$studyDuration
        }
      }

      data.frame(N, Ts)
    }
  })


  TsVsN2 <- eventReactive(input$calc, {
    if (kMax()==1) {
      a = lr()$resultsUnderH1$overallResults
      time = seq(tmin3(), a$accrualDuration, length.out=ntpts)
      N = rep(0, ntpts)
      Ts = rep(0, ntpts)
      for (i in 1:ntpts) {
        if (rho1() == 0 && rho2() == 0) {
          N[i] = accrual(
            time = time[i],
            accrualTime = accrualTime(),
            accrualIntensity = accrualIntensity(),
            accrualDuration = time[i])

          Ts[i] = caltime(
            nevents = a$numberOfEvents,
            allocationRatioPlanned = allocationRatioPlanned(),
            accrualTime = accrualTime(),
            accrualIntensity = accrualIntensity(),
            piecewiseSurvivalTime = piecewiseSurvivalTime(),
            lambda1 = lambda1(),
            lambda2 = lambda2(),
            gamma1 = gamma1(),
            gamma2 = gamma2(),
            accrualDuration = time[i],
            followupTime = 1200,
            fixedFollowup = a$fixedFollowup)
        } else {
          lrx <- lrsamplesize(
            beta = 1 - a$overallReject,
            kMax = kMax(),
            alpha = alpha(),
            allocationRatioPlanned = allocationRatioPlanned(),
            accrualTime = accrualTime(),
            accrualIntensity = accrualIntensity(),
            piecewiseSurvivalTime = piecewiseSurvivalTime(),
            lambda1 = lambda1(),
            lambda2 = lambda2(),
            gamma1 = gamma1(),
            gamma2 = gamma2(),
            accrualDuration = time[i],
            followupTime = NA,
            fixedFollowup = a$fixedFollowup,
            rho1 = a$rho1,
            rho2 = a$rho2,
            estimateHazardRatio = 0,
            typeOfComputation = typeOfComputation(),
            interval = interval,
            rounding = input$rounding)$resultsUnderH1

          N[i] = lrx$overallResults$numberOfSubjects
          Ts[i] = lrx$overallResults$studyDuration
        }
      }

      data.frame(N, Ts)
    }
  })



  output$design <- renderText({
    if (kMax()>1) {
      str1 <- paste0("Group-sequential trial with ",
                     kMax(), " stages")
    } else {
      str1 <- "Fixed design"
    }

    str2 <- paste0("Fleming-Harrington weight FH(",
                   rho1(), ", ", rho2(), ")")

    if (rho1() != 0 || rho2() != 0) {
      HTML(paste(tags$h4(str1), tags$h4(str2)))
    } else {
      HTML(paste(tags$h4(str1)))
    }
  })


  output$text <- renderText({
    a <- lr()$resultsUnderH1$overallResults
    req(a$kMax == kMax())

    str1 <- paste0("Results under H1 (alternative hypothesis)")

    str2 <- paste0("Overall power: ",
                   round(a$overallReject, 3), ", ",
                   "overall significance level (1-sided): ",
                   round(a$alpha, 4))

    if (kMax()>1) {
      str3 <- paste0("Maximum # events: ",
                     round(a$numberOfEvents, 1), ", ",
                     "expected # events: ",
                     round(a$expectedNumberOfEvents, 1))

      str4 <- paste0("Maximum # dropouts: ",
                     round(a$numberOfDropouts, 1), ", ",
                     "expected # dropouts: ",
                     round(a$expectedNumberOfDropouts, 1))

      str5 <- paste0("Maximum # subjects: ",
                     round(a$numberOfSubjects, 1), ", ",
                     "expected # subjects: ",
                     round(a$expectedNumberOfSubjects, 1))

      str6 <- paste0("Total study duration: ",
                     round(a$studyDuration, 1), ", ",
                     "expected study duration: ",
                     round(a$expectedStudyDuration, 1))

    } else {
      str3 <- paste0("Number of events: ",
                     round(a$numberOfEvents, 1))

      str4 <- paste0("Number of dropouts: ",
                     round(a$numberOfDropouts, 1))

      str5 <- paste0("Number of subjects: ",
                     round(a$numberOfSubjects, 1))

      str6 <- paste0("Study duration: ",
                     round(a$studyDuration, 1))
    }

    str7 <- paste0("Accrual duration: ",
                   round(a$accrualDuration, 1), ", ",
                   "follow-up duration: ",
                   round(a$followupTime, 1), ", ",
                   "fixed follow-up: ", a$fixedFollowup)

    paste(paste("<b>", str1, "</b>", "<br>"),
          paste(str2, str3, str4, str5, str6, str7, sep='<br/>'))

  })


  output$text0 <- renderText({
    a <- lr()$resultsUnderH0$overallResults
    req(a$kMax == kMax())

    str1 <- paste0("Results under H0 (null hypothesis)")

    str2 <- paste0("Overall type I error: ",
                   round(a$overallReject, 3), ", ",
                   "overall significance level (1-sided): ",
                   round(a$alpha, 4))

    if (kMax()>1) {
      str3 <- paste0("Maximum # events: ",
                     round(a$numberOfEvents, 1), ", ",
                     "expected # events: ",
                     round(a$expectedNumberOfEvents, 1))

      str4 <- paste0("Maximum # dropouts: ",
                     round(a$numberOfDropouts, 1), ", ",
                     "expected # dropouts: ",
                     round(a$expectedNumberOfDropouts, 1))

      str5 <- paste0("Maximum # subjects: ",
                     round(a$numberOfSubjects, 1), ", ",
                     "expected # subjects: ",
                     round(a$expectedNumberOfSubjects, 1))

      str6 <- paste0("Total study duration: ",
                     round(a$studyDuration, 1), ", ",
                     "expected study duration: ",
                     round(a$expectedStudyDuration, 1))

    } else {
      str3 <- paste0("Number of events: ",
                     round(a$numberOfEvents, 1))

      str4 <- paste0("Number of dropouts: ",
                     round(a$numberOfDropouts, 1))

      str5 <- paste0("Number of subjects: ",
                     round(a$numberOfSubjects, 1))

      str6 <- paste0("Study duration: ",
                     round(a$studyDuration, 1))
    }

    str7 <- paste0("Accrual duration: ",
                   round(a$accrualDuration, 1), ", ",
                   "follow-up duration: ",
                   round(a$followupTime, 1), ", ",
                   "fixed follow-up: ", a$fixedFollowup)

    paste(paste("<b>", str1, "</b>", "<br>"),
          paste(str2, str3, str4, str5, str6, str7, sep='<br/>'))

  })



  output$table <- renderPrint({
    a <- lr()$resultsUnderH1$overallResults
    req(a$kMax == kMax())

    if (kMax()>1) {
      b <- lr()$resultsUnderH1$byStageResults %>%
        mutate(efficacyBounds = ifelse(efficacyStopping, efficacyBounds, NA),
               futilityBounds = ifelse(futilityStopping, futilityBounds, NA),
               efficacyHR = ifelse(efficacyStopping, efficacyHR, NA),
               futilityHR = ifelse(futilityStopping, futilityHR, NA),
               efficacyP = ifelse(efficacyStopping, efficacyP, NA),
               futilityP = ifelse(futilityStopping, futilityP, NA))

      # only keep necessary variables
      b <- b[, c("informationRates", "efficacyBounds", "futilityBounds",
                 "cumulativeRejection", "cumulativeFutility",
                 "cumulativeAlphaSpent",
                 "numberOfEvents", "numberOfDropouts", "numberOfSubjects",
                 "analysisTime", "efficacyHR", "futilityHR",
                 "efficacyP", "futilityP", "information", "HR")]

      # format number of digits after decimal for each column
      j1 <- c(7,8,9,10)
      j2 <- 15
      j3 <- c(1,2,3,4,5,11,12,16)
      j4 <- c(6,13,14)

      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)


      if (input$bsf != 'none') {
        df = t(b)
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z-scale)",
                         "Futility boundary (Z-scale)",
                         "Cumulative rejection",
                         "Cumulative futility",
                         "Cumulative alpha spent",
                         "Number of events",
                         "Number of dropouts",
                         "Number of subjects",
                         "Analysis time",
                         "Efficacy boundary (HR-scale)",
                         "Futility boundary (HR-scale)",
                         "Efficacy boundary (p-scale)",
                         "Futility boundary (p-scale)",
                         "Information",
                         "HR")

      } else {
        df = t(b[,c(1,2,4,6,7,8,9,10,11,13,15,16)])
        rownames(df) = c("Information rate",
                         "Efficacy boundary (Z-scale)",
                         "Cumulative rejection",
                         "Cumulative alpha spent",
                         "Number of events",
                         "Number of dropouts",
                         "Number of subjects",
                         "Analysis time",
                         "Efficacy boundary (HR-scale)",
                         "Efficacy boundary (p-scale)",
                         "Information",
                         "HR")
      }

      colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")

      print(df, quote=FALSE)

    }
  })



  output$table0 <- renderPrint({
    a <- lr()$resultsUnderH0$overallResults
    req(a$kMax == kMax())

    if (kMax()>1) {
      b <- lr()$resultsUnderH0$byStageResults %>%
        mutate(efficacyBounds = ifelse(efficacyStopping, efficacyBounds, NA),
               futilityBounds = ifelse(futilityStopping, futilityBounds, NA),
               efficacyP = ifelse(efficacyStopping, efficacyP, NA),
               futilityP = ifelse(futilityStopping, futilityP, NA))

      # only keep necessary variables
      b <- b[, c("informationRates", "efficacyBounds", "futilityBounds",
                 "cumulativeRejection", "cumulativeFutility",
                 "cumulativeAlphaSpent", "numberOfEvents",
                 "numberOfDropouts", "numberOfSubjects", "analysisTime",
                 "efficacyP", "futilityP", "information")]

      # format number of digits after decimal for each column
      j1 <- c(7,8,9,10)
      j2 <- 13
      j3 <- c(1,2,3,4,5)
      j4 <- c(6,11,12)

      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j2] <- lapply(b[j2], formatC, format = "f", digits = 2)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)
      b[j4] <- lapply(b[j4], formatC, format = "f", digits = 4)


      if (input$bsf != 'none') {
        df = t(b[,c(4,5,7,8,9,10,13)])
        rownames(df) = c("Cumulative rejection",
                         "Cumulative futility",
                         "Number of events",
                         "Number of dropouts",
                         "Number of subjects",
                         "Analysis time",
                         "Information"
        )
      } else {
        df = t(b[,c(4,7,8,9,10,13)])
        rownames(df) = c("Cumulative rejection",
                         "Number of events",
                         "Number of dropouts",
                         "Number of subjects",
                         "Analysis time",
                         "Information"
        )
      }

      colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")

      print(df, quote=FALSE)

    }
  })



  output$plot <- renderPlotly({
    a <- lr()$resultsUnderH1$overallResults
    req(a$kMax == kMax())
    kMax <- a$kMax

    if (lr()$resultsUnderH1$settings$typeAlphaSpending == "none") {
      b <- lr()$resultsUnderH1$byStageResults
      b[-nrow(b), "efficacyStopping"] = FALSE
      b <- b %>%
        mutate(efficacyBounds = ifelse(efficacyStopping, efficacyBounds, NA),
               cumulativeRejection = ifelse(efficacyStopping,
                                            cumulativeRejection, NA),
               efficacyHR = ifelse(efficacyStopping, efficacyHR, NA),
               efficacyP = ifelse(efficacyStopping, efficacyP, NA))
    } else {
      b <- lr()$resultsUnderH1$byStageResults %>%
        mutate(efficacyBounds = ifelse(efficacyStopping, efficacyBounds, NA),
               efficacyHR = ifelse(efficacyStopping, efficacyHR, NA),
               efficacyP = ifelse(efficacyStopping, efficacyP, NA))
    }

    if (lr()$resultsUnderH1$settings$typeBetaSpending == "none") {
      b <- b %>%
        mutate(futilityStopping = FALSE,
               futilityBounds = ifelse(futilityStopping, futilityBounds, NA),
               cumulativeFutility = ifelse(futilityStopping,
                                           cumulativeFutility, NA),
               futilityHR = ifelse(futilityStopping, futilityHR, NA),
               futilityP = ifelse(futilityStopping, futilityP, NA))
    } else {
      b <- b %>%
        mutate(futilityBounds = ifelse(futilityStopping, futilityBounds, NA),
               futilityHR = ifelse(futilityStopping, futilityHR, NA),
               futilityP = ifelse(futilityStopping, futilityP, NA))
    }



    if (kMax > 1) {

      if (input$plottype == "boundaryZ") {
        af <- b %>%
          pivot_longer(cols = c("efficacyBounds", "futilityBounds"),
                       names_to = "type",
                       values_to = "bounds") %>%
          mutate(type = gsub("Bounds", "", type)) %>%
          filter(!is.na(bounds))

        plotly::plot_ly(af, x=~numberOfEvents, y=~bounds,
                        type="scatter", mode="markers+lines",
                        linetype=~type) %>%
          plotly::layout(xaxis = list(title = "Events"),
                         yaxis = list(title = "Boundaries (Z)"))
        } else if (input$plottype == "boundaryHR") {
        af <- b %>%
          pivot_longer(cols = c("efficacyHR", "futilityHR"),
                       names_to = "type",
                       values_to = "bounds") %>%
          mutate(type = gsub("HR", "", type)) %>%
          filter(!is.na(bounds))

        plotly::plot_ly(af, x=~numberOfEvents, y=~bounds,
                        type="scatter", mode="markers+lines",
                        linetype=~type) %>%
          plotly::layout(xaxis = list(title = "Events"),
                         yaxis = list(title = "Boundaries (HR)"))
      } else if (input$plottype == "boundaryP") {
        af <- b %>%
          pivot_longer(cols = c("efficacyP", "futilityP"),
                       names_to = "type",
                       values_to = "bounds") %>%
          mutate(type = gsub("P", "", type)) %>%
          filter(!is.na(bounds))

        plotly::plot_ly(af, x=~numberOfEvents, y=~bounds,
                        type="scatter", mode="markers+lines",
                        linetype=~type) %>%
          plotly::layout(xaxis = list(title = "Events"),
                         yaxis = list(title = "Boundaries (p)",
                                      zeroline = FALSE))
      } else if (input$plottype == "errorSpend") {
        af <- b %>%
          pivot_longer(cols = c("cumulativeAlphaSpent", "cumulativeFutility"),
                       names_to = "type",
                       values_to = "bounds") %>%
          mutate(type = ifelse(type=="cumulativeAlphaSpent",
                               "alpha", "beta")) %>%
          filter(!is.na(bounds))

        plotly::plot_ly(af, x=~numberOfEvents, y=~bounds,
                        type="scatter", mode="markers+lines",
                        linetype=~type) %>%
          plotly::layout(xaxis = list(title = "Events"),
                         yaxis = list(title = "Cumulative error spent",
                                      zeroline = FALSE))
      } else if (input$plottype == "eventPred") {
        time <- seq(0, a$studyDuration, length.out=100)
        df <- lrstat(time,
                     allocationRatioPlanned = allocationRatioPlanned(),
                     accrualTime = accrualTime(),
                     accrualIntensity = accrualIntensity(),
                     piecewiseSurvivalTime = piecewiseSurvivalTime(),
                     lambda1 = lambda1(),
                     lambda2 = lambda2(),
                     gamma1 = gamma1(),
                     gamma2 = gamma2(),
                     accrualDuration = a$accrualDuration,
                     followupTime = a$followupTime,
                     fixedFollowup = a$fixedFollowup,
                     rho1 = a$rho1,
                     rho2 = a$rho2,
                     predictTarget = 1)

        af <- df %>%
          pivot_longer(cols = c("subjects", "nevents"),
                       names_to = "type",
                       values_to = "n") %>%
          mutate(type = ifelse(type=="subjects", "subjects", "events"))

        plotly::plot_ly(af, x=~time, y=~n, linetype=~type,
                        type="scatter", mode="lines") %>%
          plotly::layout(xaxis = list(title = "Time", zeroline = FALSE),
                         yaxis = list(zeroline=FALSE),
                         legend = list(x=0, y=1.1, orientation="h"))
      } else if (input$plottype == "powerVsTf") {
        plotly::plot_ly(powerVsTf(), x=~time, y=~power,
                        type="scatter", mode="lines") %>%
          plotly::layout(
            xaxis = list(title = "Follow-up duration", zeroline = FALSE),
            yaxis = list(title = "Power"),
            title = list(text = paste0(
              "For ", round(a$numberOfSubjects), " subjects",
              " enrolled over ", round(a$accrualDuration, 1),
              " time units"), x = 0, xref='paper')
          )
      } else if (input$plottype == "powerVsN") {
        plotly::plot_ly(powerVsN(), x=~N, y=~power,
                        type="scatter", mode="lines") %>%
          plotly::layout(
            xaxis = list(title = "Sample size", zeroline = FALSE),
            yaxis = list(title = "Power"),
            title = list(text = paste0(
              "For follow-up duration of ", round(a$followupTime, 1),
              " time units"), x = 0, xref='paper')
          )
      } else if (input$plottype == "TsVsN") {
        plotly::plot_ly(TsVsN(), x=~N, y=~Ts,
                        type="scatter", mode="lines") %>%
          plotly::layout(
            xaxis = list(title = "Sample size", zeroline = FALSE),
            yaxis = list(title = "Study duration"),
            title = list(text =  paste0(
              "For power of ", round(a$overallReject, 3)),
              x = 0, xref='paper')
          )
      }
    } else {
      if (input$plottype2 == "eventPred") {
        time <- seq(0, a$studyDuration, length.out=100)
        df <- lrstat(time,
                     allocationRatioPlanned = allocationRatioPlanned(),
                     accrualTime = accrualTime(),
                     accrualIntensity = accrualIntensity(),
                     piecewiseSurvivalTime = piecewiseSurvivalTime(),
                     lambda1 = lambda1(),
                     lambda2 = lambda2(),
                     gamma1 = gamma1(),
                     gamma2 = gamma2(),
                     accrualDuration = a$accrualDuration,
                     followupTime = a$followupTime,
                     fixedFollowup = a$fixedFollowup,
                     rho1 = a$rho1,
                     rho2 = a$rho2,
                     predictTarget = 1)
        af <- df %>%
          pivot_longer(cols = c("subjects", "nevents"),
                       names_to = "type",
                       values_to = "n") %>%
          mutate(type = ifelse(type=="subjects", "subjects", "events"))

        plotly::plot_ly(af, x=~time, y=~n, linetype=~type,
                        type="scatter", mode="lines") %>%
          plotly::layout(xaxis = list(title = "Time", zeroline = FALSE),
                         yaxis = list(zeroline=FALSE),
                         legend = list(x=0, y=1.1, orientation="h"))
      } else if (input$plottype2 == "powerVsTf") {
        plotly::plot_ly(powerVsTf2(), x=~time, y=~power,
                        type="scatter", mode="lines") %>%
          plotly::layout(
            xaxis = list(title = "Follow-up duration", zeroline = FALSE),
            yaxis = list(title = "Power"),
            title = list(text = paste0(
              "For ", round(a$numberOfSubjects), " subjects",
              " enrolled over ", round(a$accrualDuration, 1),
              " time units"), x = 0, xref='paper')
          )
      } else if (input$plottype2 == "powerVsN") {
        plotly::plot_ly(powerVsN2(), x=~N, y=~power,
                        type="scatter", mode="lines") %>%
          plotly::layout(
            xaxis = list(title = "Sample size", zeroline = FALSE),
            yaxis = list(title = "Power"),
            title = list(text = paste0(
              "For follow-up duration of ", round(a$followupTime, 1),
              " time units"), x = 0, xref='paper')
          )
      } else if (input$plottype2 == "TsVsN") {
        plotly::plot_ly(TsVsN2(), x=~N, y=~Ts,
                        type="scatter", mode="lines") %>%
          plotly::layout(
            xaxis = list(title = "Sample size", zeroline = FALSE),
            yaxis = list(title = "Study duration"),
            title = list(text =  paste0(
              "For power of ", round(a$overallReject, 3)),
              x = 0, xref='paper')
          )
      }
    }

  })



  # simulation

  observe({
    freezeReactiveValue(input, "boundaries")

    b <- lr()$resultsUnderH1$byStageResults

    b0 <- lr()$resultsUnderH0$byStageResults

    kMax <- nrow(b)

    updateMatrixInput(
      session,
      inputId = "boundaries",
      value = matrix(c(round(b0$information/b0$information[kMax], 4),
                       round(b$numberOfEvents),
                       round(b$efficacyBounds, 3),
                       round(b$futilityBounds, 3)),
                     ncol = 4,
                     dimnames = list(paste0("Look ", seq_len(kMax)),
                                     c("Information rates",
                                       "Planned events",
                                       "Efficacy boundary",
                                       "Futility boundary"))))
  })



  # edit check for number of replications
  nIterations <- reactive({
    req(input$nIterations)
    valid <- (input$nIterations >= 1)
    shinyFeedback::feedbackWarning(
      "nIterations", !valid,
      "Number of replications must be positive integers")
    req(valid)
    round(as.numeric(input$nIterations))
  })


  # edit check for number of replications
  nRawDatasets <- reactive({
    req(input$nRawDatasets, input$nIterations)
    valid <- (input$nRawDatasets >= 0) &&
      (input$nRawDatasets <= input$nIterations)
    shinyFeedback::feedbackWarning(
      "nRawDatasets", !valid,
      paste0("Number of raw datasets per stage must be nonnegative integers",
             " and less than or equal to the number of replications"))
    req(valid)
    round(as.numeric(input$nRawDatasets))
  })


  sim <- eventReactive(input$sim, {
    req(nrow(input$boundaries)==as.numeric(input$kMax))
    r <- as.numeric(input$allocationRatioPlanned)
    if (r!= round(r)) {
      allocation <- as.vector(unlist(strsplit(attr(
        fractions(r), "fracs"), "/")), "numeric")
    } else {
      allocation <- c(r, 1)
    }

    l <- lrsim(
      kMax = as.numeric(input$kMax),
      informationRates = as.vector(input$boundaries[,1], "numeric"),
      criticalValues = as.vector(input$boundaries[,3], "numeric"),
      futilityBounds = as.vector(input$boundaries[,4], "numeric"),
      allocation1 = allocation[1],
      allocation2 = allocation[2],
      accrualTime = as.vector(input$accrual[,1], "numeric"),
      accrualIntensity = as.vector(input$accrual[,2], "numeric"),
      piecewiseSurvivalTime = as.vector(input$survival[,1], "numeric"),
      lambda1 = as.vector(input$survival[,2], "numeric"),
      lambda2 = as.vector(input$survival[,3], "numeric"),
      gamma1 = as.vector(input$dropout[,1], "numeric"),
      gamma2 = as.vector(input$dropout[,2], "numeric"),
      accrualDuration = input$accrualDuration,
      followupTime = input$followupTime,
      fixedFollowup = input$fixedFollowup,
      rho1 = input$rho1,
      rho2 = input$rho2,
      plannedEvents = as.vector(input$boundaries[,2], "numeric"),
      maxNumberOfIterations = input$nIterations,
      maxNumberOfRawDatasetsPerStage = input$nRawDatasets,
      seed = input$seed
    )
    l
  })


  sim0 <- eventReactive(input$sim, {
    req(nrow(input$boundaries)==as.numeric(input$kMax))
    r <- as.numeric(input$allocationRatioPlanned)
    if (r!= round(r)) {
      allocation <- as.vector(unlist(strsplit(attr(
        fractions(r), "fracs"), "/")), "numeric")
    } else {
      allocation <- c(r, 1)
    }

    l <- lrsim(
      kMax = as.numeric(input$kMax),
      informationRates = as.vector(input$boundaries[,1], "numeric"),
      criticalValues = as.vector(input$boundaries[,3], "numeric"),
      futilityBounds = as.vector(input$boundaries[,4], "numeric"),
      allocation1 = allocation[1],
      allocation2 = allocation[2],
      accrualTime = as.vector(input$accrual[,1], "numeric"),
      accrualIntensity = as.vector(input$accrual[,2], "numeric"),
      piecewiseSurvivalTime = as.vector(input$survival[,1], "numeric"),
      lambda1 = as.vector(input$survival[,3], "numeric"),
      lambda2 = as.vector(input$survival[,3], "numeric"),
      gamma1 = as.vector(input$dropout[,1], "numeric"),
      gamma2 = as.vector(input$dropout[,2], "numeric"),
      accrualDuration = input$accrualDuration,
      followupTime = input$followupTime,
      fixedFollowup = input$fixedFollowup,
      rho1 = input$rho1,
      rho2 = input$rho2,
      plannedEvents = as.vector(input$boundaries[,2], "numeric"),
      maxNumberOfIterations = input$nIterations,
      maxNumberOfRawDatasetsPerStage = input$nRawDatasets,
      seed = input$seed
    )
    l
  })



  output$simtext <- renderText({
    req(nrow(input$boundaries)==sim()$overview$kMax)
    a <- sim()$overview
    str1 <- "Simulation results under H1"
    str2 <- paste0("Overall rejection, ",
                   round(a$overallReject,3), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents, 1), ", ",
                   "expected # dropouts: ",
                   round(a$expectedNumberOfDropouts, 1), ", ",
                   "# subjects: ",
                   round(a$expectedNumberOfSubjects, 1), ", ",
                   "study duration: ",
                   round(a$expectedStudyDuration, 1))
    paste(paste("<b>", str1, "</b>"),
          str2, sep = '<br/>')
  })



  output$simtable <- renderPrint({
    req(nrow(input$boundaries)==sim()$overview$kMax)
    if (kMax()>1) {
      a <- sim()$overview

      b <- data.frame(a$cumulativeRejection,
                      a$cumulativeFutility,
                      a$numberOfEvents,
                      a$numberOfDropouts,
                      a$numberOfSubjects,
                      a$analysisTime)

      j1 <- c(3,4,5,6)
      j3 <- c(1,2)

      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)

      df = t(b)
      rownames(df) = c("Cumulative rejection",
                       "Cumulative futility",
                       "Number of events",
                       "Number of dropouts",
                       "Number of subjects",
                       "Analysis time")

      colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")

      print(df, quote=FALSE)

    }
  })



  output$downloadSumdata <- downloadHandler(
    filename = function() {
      paste0("sim-sumdata-h1-", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(sim()$sumdata, file, row.names=FALSE)
    }
  )

  output$downloadRawdata <- downloadHandler(
    filename = function() {
      paste0("sim-rawdata-h1-", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(sim()$rawdata, file, row.names=FALSE)
    }
  )


  output$simtext0 <- renderText({
    req(nrow(input$boundaries)==sim0()$overview$kMax)
    a <- sim0()$overview
    str1 <- "Simulation results under H0"
    str2 <- paste0("Overall rejection: ",
                   round(a$overallReject,3), ", ",
                   "expected # events: ",
                   round(a$expectedNumberOfEvents), ", ",
                   "expected # dropouts: ",
                   round(a$expectedNumberOfDropouts, 1), ", ",
                   "# subjects: ",
                   round(a$expectedNumberOfSubjects), ", ",
                   "study duration: ",
                   round(a$expectedStudyDuration, 1))
    paste(paste("<b>", str1, "</b>"),
          str2, sep = '<br/>')
  })


  output$simtable0 <- renderPrint({
    req(nrow(input$boundaries)==sim0()$overview$kMax)
    if (kMax()>1) {
      a <- sim0()$overview

      b <- data.frame(a$cumulativeRejection,
                      a$cumulativeFutility,
                      a$numberOfEvents,
                      a$numberOfDropouts,
                      a$numberOfSubjects,
                      a$analysisTime)

      j1 <- c(3,4,5,6)
      j3 <- c(1,2)

      b[j1] <- lapply(b[j1], formatC, format = "f", digits = 1)
      b[j3] <- lapply(b[j3], formatC, format = "f", digits = 3)

      df = t(b)
      rownames(df) = c("Cumulative rejection",
                       "Cumulative futility",
                       "Number of events",
                       "Number of dropouts",
                       "Number of subjects",
                       "Analysis time")

      colnames(df) <- paste("Stage", seq_len(ncol(df)), sep=" ")

      print(df, quote=FALSE)

    }
  })


  output$downloadSumdata0 <- downloadHandler(
    filename = function() {
      paste0("sim-sumdata-h0-", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(sim0()$sumdata, file, row.names=FALSE)
    }
  )

  output$downloadRawdata0 <- downloadHandler(
    filename = function() {
      paste0("sim-rawdata-h0-", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(sim0()$rawdata, file, row.names=FALSE)
    }
  )



  output$lrp <- renderText({
    a <- lr()$resultsUnderH1$overallResults
    req(a$kMax == kMax())

    b <- lr()$resultsUnderH1$byStageResults
    s <- lr()$resultsUnderH1$settings

    str2 <- paste0("kMax = ", a$kMax, ",")

    if (a$kMax > 1) {
      str3 <- paste0("informationRates = c(",
                     paste(b$informationRates, collapse=", "), "),")

      if (s$typeAlphaSpending != "none") {
        str3 <- paste(str3,
                      paste0("efficacyStopping = c(",
                             paste(as.vector(b$efficacyStopping, "numeric"),
                                   collapse=", "), "),"), sep="<br/>")
      }

      if (s$typeBetaSpending != "none") {
        str3 <- paste(str3,
                      paste0("futilityStopping = c(",
                             paste(as.vector(b$futilityStopping, "numeric"),
                                   collapse=", "), "),"), sep="<br/>")
      }

      str7 <- paste0("typeAlphaSpending = '", s$typeAlphaSpending, "',")
      if (s$typeAlphaSpending %in% c("WT", "sfKD", "sfHSD")) {
        str7 <- paste(str7, paste0("parameterAlphaSpending = ",
                                   s$parameterAlphaSpending, ","), sep="<br/>")
      }
      if (s$typeAlphaSpending == "user") {
        str7 <- paste(str7, paste0("userAlphaSpending = c(",
                                   paste(s$userAlphaSpending, collapse=", "),
                                   ","), sep="<br/>")
      }

      str8 <- paste0("typeBetaSpending = '", s$typeBetaSpending, "',")
      if (s$typeBetaSpending %in% c("sfKD", "sfHSD")) {
        str8 <- paste(str8, paste0("parameterBetaSpending = ",
                                   s$parameterBetaSpending, ","), sep="<br/>")
      }

    }

    str6 <- paste0("alpha = ", round(a$alpha, 5), ",")
    str9 <- paste0("allocationRatioPlanned = ", s$allocationRatioPlanned, ",")

    if (length(s$accrualTime) == 1) {
      str10 <- paste0("accrualTime = ", s$accrualTime, ",")
      str10 <- paste(str10, paste0("accrualIntensity = ",
                                   s$accrualIntensity, ","), sep="<br/>")
    } else {
      str10 <- paste0("accrualTime = c(",
                      paste(s$accrualTime, collapse=", "), "),")
      str10 <- paste(str10, paste0("accrualIntensity = c(",
                                   paste(s$accrualIntensity, collapse=", "),
                                   "),"), sep="<br/>")
    }

    if (length(s$piecewiseSurvivalTime) == 1) {
      str11 <- paste0("piecewiseSurvivalTime = ", s$piecewiseSurvivalTime, ",")
      str11 <- paste(str11, paste0("lambda1 = ", s$lambda1, ","), sep="<br/>")
      str11 <- paste(str11, paste0("lambda2 = ", s$lambda2, ","), sep="<br/>")
    } else {
      str11 <- paste0("piecewiseSurvivalTime = c(",
                      paste(s$piecewiseSurvivalTime, collapse=", "), "),")
      str11 <- paste(str11, paste0("lambda1 = c(",
                                   paste(s$lambda1, collapse=", "),
                                   "),"), sep="<br/>")
      str11 <- paste(str11, paste0("lambda2 = c(",
                                   paste(s$lambda2, collapse=", "),
                                   "),"), sep="<br/>")
    }


    str12 <- paste0("gamma1 = ", s$gamma1, ",")
    str12 <- paste(str12, paste0("gamma2 = ", s$gamma2, ","), sep="<br/>")

    str15 <- paste0("fixedFollowup = ", a$fixedFollowup, ",")
    str16 <- paste0("rho1 = ", a$rho1, ",")
    str16 <- paste(str16, paste0("rho2 = ", a$rho2, ","), sep="<br/>")


    if (input$target == 'power') {
      str1 <- paste0("(lr <- lrpower(")

      str13 <- paste0("accrualDuration = ", a$accrualDuration, ",")
      str14 <- paste0("followupTime = ", a$followupTime, ",")
      str17 <- paste0("typeOfComputation = '", a$typeOfComputation, "'))")

      if (a$kMax == 1) {
        paste(str1, str2, str6, str9, str10,
              str11, str12, str13, str14, str15, str16, str17, sep = '<br/>')
      } else {
        paste(str1, str2, str3, str6, str7, str8, str9, str10,
              str11, str12, str13, str14, str15, str16, str17, sep = '<br/>')
      }
    } else if (input$target == 'accrualDuration') {
      str0 <- paste0("(lr <- lrsamplesize(")
      str1 <- paste0("beta = ", round(1-a$overallReject, 4), ",")

      str13 <- paste0("accrualDuration = NA,")
      str14 <- paste0("followupTime = ", a$followupTime, ",")
      str17 <- paste0("typeOfComputation = '", a$typeOfComputation, "',")
      str18 <- paste0("interval = c(", paste(interval, collapse=", "), ")))")

      if (a$kMax == 1) {
        paste(str0, str1, str2, str6, str9, str10,
              str11, str12, str13, str14, str15, str16, str17,
              str18, sep = '<br/>')
      } else {
        paste(str0, str1, str2, str3, str6, str7, str8, str9, str10,
              str11, str12, str13, str14, str15, str16, str17,
              str18, sep = '<br/>')
      }
    } else if (input$target == 'followupTime') {
      str0 <- paste0("(lr <- lrsamplesize(")
      str1 <- paste0("beta = ", round(1-a$overallReject, 4), ",")

      str13 <- paste0("accrualDuration = ", a$accrualDuration, ",")
      str14 <- paste0("followupTime = NA,")
      str17 <- paste0("typeOfComputation = '", a$typeOfComputation, "',")
      str18 <- paste0("interval = c(", paste(interval, collapse=", "), ")),")
      str19 <- paste("rounding = ", input$rounding, "))")

      if (a$kMax == 1) {
        paste(str0, str1, str2, str6, str9, str10,
              str11, str12, str13, str14, str15, str16, str17,
              str18,  str19, sep = '<br/>')
      } else {
        paste(str0, str1, str2, str3, str6, str7, str8, str9, str10,
              str11, str12, str13, str14, str15, str16, str17,
              str18, str19, sep = '<br/>')
      }
    }

  })


  output$lrs <- renderText({
    a <- lr()$resultsUnderH1$overallResults
    req(a$kMax == kMax())

    s <- lr()$resultsUnderH1$settings

    str1 <- paste0("(sim <- lrsim(")

    str2 <- paste0("kMax = ", a$kMax, ",")

    if (a$kMax > 1) {
      str3 <- paste0("informationRates = c(",
                     paste(as.vector(input$boundaries[,1], "numeric"),
                           collapse=", "), "),")
      str4 <- paste0("criticalValues = c(",
                     paste(as.vector(input$boundaries[,3], "numeric"),
                           collapse=", "), "),")
      str5 <- paste0("futilityBounds = c(",
                     paste(as.vector(input$boundaries[,4], "numeric"),
                           collapse=", "), "),")
    } else {
      str4 <- paste0("criticalValues = ", input$boundaries[,3], ",")
    }

    r <-  lr()$resultsUnderH1$settings$allocationRatioPlanned
    if (r!= round(r)) {
      allocation <- as.vector(unlist(strsplit(attr(
        fractions(r), "fracs"), "/")), "numeric")
    } else {
      allocation <- c(r, 1)
    }

    str6 <- paste0("allocation1 = ", allocation[1], ",")
    str6 <- paste(str6, paste0("allocation2 = ",
                               allocation[2], ","), sep='<br/>')

    if (length(s$accrualTime) == 1) {
      str7 <- paste0("accrualTime = ", s$accrualTime, ",")
      str7 <- paste(str7, paste0("accrualIntensity = ",
                                 s$accrualIntensity, ","), sep="<br/>")
    } else {
      str7 <- paste0("accrualTime = c(", paste(s$accrualTime,
                                               collapse=", "), "),")
      str7 <- paste(str7, paste0("accrualIntensity = c(",
                                 paste(s$accrualIntensity, collapse=", "),
                                 "),"), sep="<br/>")
    }

    if (length(s$piecewiseSurvivalTime) == 1) {
      str8 <- paste0("piecewiseSurvivalTime = ", s$piecewiseSurvivalTime, ",")
      str8 <- paste(str8, paste0("lambda1 = ", s$lambda1, ","), sep="<br/>")
      str8 <- paste(str8, paste0("lambda2 = ", s$lambda2, ","), sep="<br/>")
    } else {
      str8 <- paste0("piecewiseSurvivalTime = c(",
                     paste(s$piecewiseSurvivalTime, collapse=", "), "),")
      str8 <- paste(str8, paste0("lambda1 = c(",
                                 paste(s$lambda1, collapse=", "),
                                 "),"), sep="<br/>")
      str8 <- paste(str8, paste0("lambda2 = c(",
                                 paste(s$lambda2, collapse=", "),
                                 "),"), sep="<br/>")
    }

    str9 <- paste0("gamma1 = ", s$gamma1, ",")
    str9 <- paste(str9, paste0("gamma2 = ", s$gamma2, ","), sep="<br/>")

    str10 <- paste0("accrualDuration = ", round(a$accrualDuration, 2), ",")
    str11 <- paste0("followupTime = ", round(a$followupTime, 2), ",")

    str12 <- paste0("fixedFollowup = ", a$fixedFollowup, ",")
    str13 <- paste0("rho1 = ", a$rho1, ",")
    str13 <- paste(str13, paste0("rho2 = ", a$rho2, ","), sep="<br/>")

    if (a$kMax == 1) {
      str14 <- paste("plannedEvents = ", input$bounardies[,2], ",")
    } else {
      str14 <- paste0("plannedEvents = c(",
                      paste(as.vector(input$boundaries[,2], "numeric"),
                            collapse=", "), "),")
    }

    str15 <- paste0("maxNumberOfIterations = ", input$nIterations, ",")
    str16 <- paste0("maxNumberOfRawDatasetsPerStage = ",
                    input$nRawDatasets, ",")
    str17 <- paste0("seed = ", input$seed, "))")

    if (a$kMax == 1) {
      paste(str1, str2, str4, str6, str7, str8, str9, str10,
            str11, str12, str13, str14, str15, str16, str17, sep = '<br/>')
    } else {
      paste(str1, str2, str3, str4, str5, str6, str7, str8, str9, str10,
            str11, str12, str13, str14, str15, str16, str17, sep = '<br/>')
    }
  })



  # save inputs
  output$saveInputs <- downloadHandler(
    filename = function() {
      paste0("inputs-", Sys.Date(), ".RData")
    },

    content = function(file) {
      x <- list(
        kMax = kMax(),
        target = input$target,
        alpha = alpha(),
        power = 1-beta(),
        accrualDuration = accrualDuration(),
        followupTime = followupTime(),
        allocationRatioPlanned = allocationRatioPlanned(),
        fixedFollowup = input$fixedFollowup,
        rho1 = rho1(),
        rho2 = rho2(),
        rounding = input$rounding,
        xIA = matrix(informationRates(), ncol=1,
                     dimnames = list(paste0("Look ", seq_len(kMax())),
                                     c("Information rate"))
        ),
        xES = efficacyStopping(),
        xFS = futilityStopping(),
        xUA = matrix(userAlphaSpending(), ncol=1,
                     dimnames = list(paste0("Look ", seq_len(kMax())),
                                     c("Cumulative alpha"))
        ),
        asf = typeAlphaSpending(),
        deltaAlpha = input$deltaAlpha,
        rhoAlpha = input$rhoAlpha,
        gammaAlpha = input$gammaAlpha,
        bsf = typeBetaSpending(),
        rhoBeta = input$rhoBeta,
        gammaBeta = input$gammaBeta,
        survival = matrix(as.numeric(input$survival),
                          ncol=ncol(input$survival),
                          dimnames = list(
                            NULL, c("Starting time",
                                    "Treatment hazard rate",
                                    "Control hazard rate"))
        ),
        accrual = matrix(as.numeric(input$accrual),
                         ncol=ncol(input$accrual),
                         dimnames = list(
                           NULL, c("Starting time",
                                   "Accrual intensity"))
        ),
        dropout = matrix(as.numeric(input$dropout),
                         ncol=ncol(input$dropout),
                         dimnames = list(
                           NULL, c("Treatment hazard rate",
                                   "Control hazard rate"))
        ),
        plottype = input$plottype,
        plottype2 = input$plottype2,
        boundaries = matrix(
          as.numeric(input$boundaries),
          ncol=ncol(input$boundaries),
          dimnames = list(paste0("Look ", seq_len(kMax())),
                          c("Information rates",
                            "Planned events",
                            "Efficacy boundary",
                            "Futility boundary"))),
        nIterations = nIterations(),
        nRawDatasets = nRawDatasets(),
        seed = input$seed
      )

      save(x, file = file)
    }
  )



  # load inputs
  observeEvent(input$loadInputs, {
    file <- input$loadInputs
    ext <- tools::file_ext(file$datapath)

    req(file)

    valid <- (ext == "RData")
    if (!valid) showNotification("Please upload an RData file")
    req(valid)

    load(file=file$datapath)

    updateSelectInput(session, "kMax", selected=x$kMax)
    updateRadioButtons(session, "target", selected=x$target)
    updateNumericInput(session, "alpha", value=x$alpha)
    if (x$target != "power") {
      updateNumericInput(session, "power", value=x$power)
    }
    if (x$target != "accrualDuration") {
      updateNumericInput(session, "accrualDuration", value=x$accrualDuration)
    }
    if (x$target != "followupTime") {
      updateNumericInput(session, "followupTime", value=x$followupTime)
    }
    updateNumericInput(session, "allocationRatioPlanned",
                       value=x$allocationRatioPlanned)
    updateCheckboxInput(session, "fixedFollowup", value=x$fixedFollowup)
    updateNumericInput(session, "rho1", value=x$rho1)
    updateNumericInput(session, "rho2", value=x$rho2)
    updateCheckboxInput(session, "rounding", value=x$rounding)

    if (x$kMax > 1) {
      updateMatrixInput(
        session, paste0("xIA_", x$kMax),
        value=matrix(x$xIA[-x$kMax,], ncol = 1,
                     dimnames = list(paste0("Look ", seq_len(x$kMax-1)),
                                     c("Information rate"))))
      updateCheckboxGroupInput(
        session, paste0("xES_", x$kMax),
        selected=paste0("Look ", which(x$xES[-x$kMax]==1)))
      updateCheckboxGroupInput(
        session, paste0("xFS_", x$kMax),
        selected=paste0("Look ", which(x$xFS[-x$kMax]==1)))
    }
    if (x$kMax >1 && x$asf=="user") {
      updateMatrixInput(
        session, paste0("xUA_", x$kMax),
        value=matrix(x$xUA[-x$kMax,], ncol=1,
                     dimnames = list(paste0("Look ", seq_len(x$kMax-1)),
                                     c("Cumulative alpha"))))
    }
    updateSelectInput(session, "asf", selected=x$asf)
    updateNumericInput(session, "deltaAlpha", value=x$deltaAlpha)
    updateNumericInput(session, "gammaAlpha", value=x$gammaAlpha)
    updateSelectInput(session, "bsf", selected=x$bsf)
    updateNumericInput(session, "deltaBeta", value=x$deltaBeta)
    updateNumericInput(session, "gammaBeta", value=x$gammaBeta)
    updateMatrixInput(session, "survival", value=x$survival)
    updateMatrixInput(session, "accrual", value=x$accrual)
    updateMatrixInput(session, "dropout", value=x$dropout)
    updateSelectInput(session, "plottype", selected=x$plottype)
    updateSelectInput(session, "plottype2", selected=x$plottype2)

    updateMatrixInput(session, "boundaries", value=x$boundaries)

    updateNumericInput(session, "nIterations", value=x$nIterations)
    updateNumericInput(session, "nRawDatasets", value=x$nRawDatasets)
    updateNumericInput(session, "seed", value=x$seed)

  })



}

# Run the application
shinyApp(ui = ui, server = server)
