
ui <- fluidPage(
  fluidRow(
    column(width = 4, wellPanel(
            textInput("dir", "Directory (Windows: Replace any \\ with /)"),
            textInput("mydat", "Dataset"),
            textInput("DV", "Dependent Variable (single value)"),
            textInput("treats", "Treatment (single value)"),
            textInput("subgroups", "Subgroup (single value)", value=""),
            radioButtons("model", "Model", c("linear", "logit"), selected="linear"),
            textInput("controls", "Controls (comma separated)", value=""),
            numericInput("sims", "Simulations", value=1000),
            radioButtons("comb", "Permutations or Combinations?",
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
           start<-Sys.time()

           # standardize control variables to allow for the inclusion of multiple controls
           if(str_detect(input$controls, "[,]")){
                   ctrl<-input$controls %>% strsplit(",") %>%  unlist() %>% trimws
           }else{
                   ctrl<-input$controls
           }
           dat<-paste0(input$dir,"/", input$mydat)
           gg<-TreatmentEffects(data=read.csv(dat),
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
