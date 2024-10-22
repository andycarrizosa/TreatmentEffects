library(ggplot2)
library(parallel)
library(shiny)

ui <- fluidPage(
  fluidRow(
    column(width = 4, wellPanel(
            textInput("mydat", "Data Path"),
            textInput("DV", "Dependent Variable"),
            textInput("treats", "Treatments"),
            textInput("subgroups", "Subgroups", value=""),
            textInput("model", "Model", value="linear"),
            textInput("controls", "Controls", value=""),
            numericInput("sims", "Simulations", value=1000),
            radioButtons("comb", "Perm or Comb",
                         c("perm", "comb"), selected="perm"),
            sliderInput("clust", "Clusters", min=1, max=detectCores()-2, value=1, round=TRUE, step=1),
            actionButton("go", "Go!")
            )),
    column(width = 5,
           textOutput("text"),
      plotOutput("plot1"),
      tableOutput("table")
    )
  )
)


server <- function(input, output) {

    gg<-eventReactive(input$go, {
           start<-Sys.time()
           gg<-TreatmentEffects(data=read.csv(input$mydat),
                                        treats=input$treats,
                                        DV=input$DV,
                                        model = input$model,
                                        controls = input$controls,
                                        subgroups = input$subgroups,
                                        sims = input$sims,
                                        comb = input$comb,
                                        clust = input$clust)
           end<-Sys.time()
           gg[[5]]<-paste0("Estimation time= ", round(end-start,3), " seconds")
           gg
           })

  output$plot1<-renderPlot({gg()[[4]]})
  output$table<-renderTable({gg()[[2]]})
  output$text<-renderText({gg()[[5]]})
}


shinyApp(ui, server)
