library(stringr)
library(Hmisc)
library(ggplot2)
library(clarify)
library(tidyverse)
library(RColorBrewer)
library(gtools)
library(doParallel)
library(foreach)
library(shiny)
library(knitr)
library(rmarkdown)

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
      plotOutput("plot1", height="800px", width="800px"),
      tableOutput("table"),
      tableOutput("specification")
    )
  )
)


server <- function(input, output) {

    gg<-eventReactive(input$go, {
            library(stringr)
            library(Hmisc)
            library(ggplot2)
            library(clarify)
            library(tidyverse)
            library(RColorBrewer)
            library(gtools)
            library(doParallel)
            library(foreach)
            library(shiny)
            library(knitr)
            library(rmarkdown)

           start<-Sys.time()

           # standardize control variables to allow for the inclusion of multiple controls
           if(str_detect(input$controls, "[,]")){
                   ctrl<-input$controls %>% strsplit(",") %>%  unlist() %>% trimws
           }else{
                   ctrl<-input$controls
           }

           gg<-TreatmentEffects(data=read.csv(input$mydat),
                                        treats=input$treats,
                                        DV=input$DV,
                                        model = input$model,
                                        controls = ctrl,
                                        subgroups = input$subgroups,
                                        sims = input$sims,
                                        comb = input$comb,
                                        clust = input$clust)

           end<-Sys.time()
           gg[[5]]<-paste0("Time difference of ", round(difftime(end,start, units="secs"),3), " seconds")
           gg
           })

  output$plot1<-renderPlot({gg()[[4]]})
  output$table<-renderTable({gg()[[2]]})
  output$text<-renderText({gg()[[5]]})
  output$specification<-renderTable({gg()[[3]]})
}


shiny_TE<-shinyApp(ui, server)
