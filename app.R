#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(dplyr)
library(ggplot2)

# function ----------------------------------------------------------------

zprime <- function(p, n, dispersion, center) {
    1 - 3*(dispersion(p, na.rm = TRUE) + 
               dispersion(n, na.rm = TRUE))/abs(center(p, na.rm = TRUE) - 
                                                center(n, na.rm = TRUE))
}

ssmd <- function(p, n, dispersion, center){
    (center(p, na.rm = TRUE) - center(n, na.rm = TRUE)) / sqrt(dispersion(p, na.rm = TRUE)^2 + dispersion(n, na.rm = TRUE)^2)
}

cv <- function(p, n, dispersion, center){
    pos <- (dispersion(p, na.rm = TRUE) / center(p, na.rm = TRUE))*100
    neg <- (dispersion(n, na.rm = TRUE) / center(n, na.rm = TRUE))*100
    return(list(positive = pos, negative = neg))
}

getGraph <- function(data, value, controls, title){
    data.plot <- ggplot(data, aes(x = value, color = controls)) + 
        geom_density(show.legend = TRUE) + 
        scale_color_manual(values = c("#377eb8", "#e41a1c")) +
        labs(title=title, x="", y = "")
    #scale_x_continuous(limits = c(0,xlimit)) +
    #scale_y_continuous(limits = c(0,ylimit))
    return(data.plot)
}



# ui ----------------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- bootstrapPage(
    tags$head(includeHTML("gtag.html")),
    navbarPage(
        theme = shinytheme("flatly"),
        collapsible = TRUE,
        HTML(
            '<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">Quality Control</a>'
        ),
        id = "nav",
        windowTitle = "Quality control for high-throughput screening",
        
        tabPanel("Data",
                 # Sidebar with a slider input for upload file
                 sidebarLayout(
                     sidebarPanel(
                         fileInput(
                             "file",
                             "Choose CSV File",
                             accept = c("text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")
                         ),
                         tags$hr(),
                         checkboxInput("header", "Header", TRUE)
                     ),
                     
                     # Show a file content
                     mainPanel(dataTableOutput("contents"))
                 )),
        
        tabPanel("Analyze",
                 # Sidebar with a slider input for upload file
                 sidebarLayout(
                     sidebarPanel(
                         uiOutput("columns"),
                         selectInput("inSelect", "Select positive Control",
                                     choices = NULL, multiple = FALSE),
                         selectInput("inSelect_2", "Select negative Control",
                                     choices = NULL, multiple = FALSE),
                         uiOutput("qcolumns"),
                         checkboxGroupInput("qc", "Statistical parameters:",
                                            c("Z-prime" = "zp",
                                              "SSMD" = "ssmd",
                                              "CV" = "cv"),
                                            selected = "zp"),
                         tags$hr(),
                         tags$h5("Options"),
                         checkboxInput("robust", "use robust statistics (median, mad)")
                     ),
                     
                     # Show a data content
                     mainPanel(plotOutput("plot"),
                               hr(), 
                               tableOutput("stat"))
                 ))
        
    )
)


# server ------------------------------------------------------------------

# Define server logic required to display file
server <- function(input, output, session) {

    hideTab(inputId = "nav", target = "Analyze")
    observeEvent(input$file, {
        if(!is.null(input$file)){
            showTab(inputId = "nav", target = "Analyze")
        }
    })
    
    
    output$contents <- renderDataTable({
        # input$file will be NULL initially. After the user selects
        # and uploads a file, it will be a data frame with 'name',
        # 'size', 'type', and 'datapath' columns. The 'datapath'
        # column will contain the local filenames where the data can
        # be found.
        inFile <- input$file
        
        if(is.null(inFile))
            return(NULL)

        # read the file 
        read.csv(inFile$datapath, header = input$header)
    })
    
    df <- reactive({
        inFile <- input$file
        
        if(is.null(inFile))
            return(NULL)
        
        # read the file 
        read.csv(inFile$datapath, header = input$header)
    })

    output$columns <- renderUI({
        selectInput("col", "Choose a column to select controls:",
                    choices = colnames(df())) 
                    #selected = df()[1])
    })
    
    output$qcolumns <- renderUI({
        selectInput("qcol", "Choose a quantitative variable to analyze:",
                    choices = colnames(df() %>% select(where(is.numeric)))) 
                    #selected = df()[1])
    })
    
    observe({
        x <- input$col
        # Can use character(0) to remove all choices
        if (is.null(x))
            x <- character(0)
        
        updateSelectInput(session, "inSelect",
                          label = "Select positive control:",
                          choices = unique(df()[x]))
                          #selected = tail(x, 1)
        updateSelectInput(session, "inSelect_2",
                          label = "Select negative control:",
                          choices = unique(df()[x]))
                          #selected = tail(x, 1))
    })
    
    subdata <- reactive({
        validate(
            need(input$inSelect != "", "Please select a positive control")
        )
        pos <- df() %>% 
            filter(!!as.symbol(input$col) == input$inSelect) %>% 
            select(input$qcol) %>% pull()
        validate(
            need(input$inSelect_2 != "", "Please select a negative control")
        )
        neg <- df() %>% 
            filter(!!as.symbol(input$col) == input$inSelect_2) %>% 
            select(input$qcol) %>% pull()
        return(list(positive = pos, negative = neg))
    })
    
    output$stat <- renderTable({
        x <- subdata()
        if (is.null(x))
            x <- character(0)
        validate(
            need(input$qc != "", "Please select a statistical parameter")
        )
        if(input$robust & input$qc == "zp"){
            z <- zprime(p = subdata()$positive, n = subdata()$negative,
                        dispersion = mad, center = median)
            data.frame(`Positive Control` = input$inSelect, `Negative Control` = input$inSelect_2 , "ZPrime" = z, check.names = FALSE)
        }
        else if(input$robust & input$qc == "ssmd"){
            ssmd <- ssmd(p = subdata()$positive, n = subdata()$negative,
                         dispersion = mad, center = median)
            data.frame(`Positive Control` = input$inSelect, `Negative Control` = input$inSelect_2 , "SSMD" = ssmd, check.names = FALSE)
        }
        else if(input$robust & input$qc == "cv"){
            cv <- cv(p = subdata()$positive, n = subdata()$negative,
                         dispersion = mad, center = median)
            data.frame(`Controls` = c(input$inSelect,input$inSelect_2) , "CV" = c(cv$positive,cv$negative), check.names = FALSE)
        }
        else if(input$qc == "zp"){
            z <- zprime(p = subdata()$positive, n = subdata()$negative,
                        dispersion = sd, center = mean)
            data.frame(`Positive Control` = input$inSelect, `Negative Control` = input$inSelect_2 , "ZPrime" = z, check.names = FALSE)
        }
        else if(input$qc == "ssmd"){
            ssmd <- ssmd(p = subdata()$positive, n = subdata()$negative,
                        dispersion = sd, center = mean)
            data.frame(`Positive Control` = input$inSelect, `Negative Control` = input$inSelect_2 , "SSMD" = ssmd, check.names = FALSE)
        }
        else if(input$qc == "cv"){
            cv <- cv(p = subdata()$positive, n = subdata()$negative,
                         dispersion = sd, center = mean)
            data.frame(`Controls` = c(input$inSelect,input$inSelect_2) , "CV" = c(cv$positive,cv$negative), check.names = FALSE)
        }
    })
    output$plot <- renderPlot({
        seldata <- df() %>% filter((!!as.symbol(input$col) == input$inSelect) | ( !!as.symbol(input$col) == input$inSelect_2))
        values <- as.vector(seldata %>% select(input$qcol)) %>% pull()
        ctrls <- as.vector(seldata %>% select(input$col)) %>% pull()
        getGraph(data = seldata, value = values, controls = ctrls, title = "Quality control distribution of positive and negative controls")
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
