library(shiny)
library(tidyverse)
library(knitr)
library(rmarkdown)

ui <- fluidPage(
  fluidRow(
    column(width = 4, wellPanel(
            textInput("dir", "Directory"),
            textInput("mydata", "Dataset"),
            textInput("nm", "Title for IFD Report")
            )),
    column(width = 4, wellPanel(textInput("contDVs", "Continuous DVs"),
                                textInput("binDVs", "Binary DVs"),
                                textInput("covariates", "Controls"),
                                textInput("subgroups", "Subgroups"),
                                textInput("treatment", "Treatment"),
                                numericInput("sims", "Simulations", value=1000),
                                sliderInput("clust", "Clusters", min=1, max=detectCores()-2, value=1, round=TRUE, step=1),
                                actionButton("go", "Go!")
    )),
    column(width = 2,
           textOutput("text"))
  )
)


server <- function(input, output) {

        gg<-eventReactive(input$go, {
                library(shiny)
                library(tidyverse)
                library(knitr)
                library(rmarkdown)

                start<-Sys.time()
                home<-getwd()
                setwd(input$dir)
                if(!any(list.files()=="rmd_IFD.Rmd")){
                        download.file("https://github.com/andycarrizosa/TreatmentEffects/blob/main/R/rmd_IFD.Rmd", destfile="rmd_IFD.Rmd")
                        Sys.sleep(10)

                }
                rmd<-readLines("rmd_IFD.Rmd")
                mydata<-read.csv(input$mydata)
                treatment<-input$treatment
                sims<-input$sims
                clust<-input$clust
                nm<-input$nm

                covariates<-input$covariates %>% strsplit(",") %>% unlist %>% trimws
                subgroups<-input$subgroups %>% strsplit(",") %>% unlist %>% trimws
                contDVs<-input$contDVs %>% strsplit(",") %>% unlist %>% trimws
                binDVs<-input$binDVs %>% strsplit(",") %>% unlist %>% trimws

                if(length(contDVs)==0) contDVs<-""
                if(length(binDVs)==0) binDVs<-""
                if(length(covariates)==0) covariates<-""
                if(length(subgroups)==0) subgroups<-"" else subgroups<-c("", subgroups)

                newrmd<-rmd %>% gsub("PLACEHOLDER", nm,.)
                cat(newrmd, file="newtest.Rmd", sep="\n")
                render("newtest.Rmd", "html_document")
                unlink("newtest.Rmd")

                end<-Sys.time()
                setwd(home)
                print(paste0("IFD production took ", round(difftime(end,start, units="secs"),3), " seconds"))
                })

        output$text<-renderText({gg()})

}


shiny_IFD<-shinyApp(ui, server)
