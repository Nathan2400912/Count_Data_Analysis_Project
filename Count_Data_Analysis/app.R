#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinyjs)
library(DT)
library(tidyverse)
library(rlang)

ui <- fluidPage(
    useShinyjs(),
    # Application title
    titlePanel(
      div(h2("Count Data Analysis Application"),
          h5("Thanks for using this application! Please upload different csv files of your study to continue"))),
    # Sidebar with a slider input for number of bins 
    tabsetPanel(
      id = "tabs",
      tabPanel('Samples',
               sidebarLayout(
                 sidebarPanel(
                   fileInput('sampleinfo', 'Please upload your sample information below', accept = '.csv',
                             buttonLabel = 'Browse...', placeholder = 'Example_count_data.csv'),
                   uiOutput('samplecolumns'),
                   uiOutput('samplexcolumn')
                 ),
                 mainPanel(
                   tabsetPanel(
                     id = "sampletabs",
                     tabPanel('Summary',
                              tableOutput('summarytable')),
                     tabPanel('Table',
                              DT::dataTableOutput('sampletable')),
                     tabPanel('Plots',
                              plotOutput('densityplots'))
                   )
                 )
               )),
      tabPanel('Counts',
               sidebarLayout(
                 sidebarPanel(
                   fileInput('countdata', 'Please upload your normalized count matrix below', accept = '.csv',
                             buttonLabel = 'Browse...', placeholder = 'Example_count_data.csv'),
                   uiOutput('variance_slider'),
                   uiOutput('nonzero_slider')
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel('Summary'),
                     tabPanel('Table'),
                     tabPanel('Plots')
                   )
                 )
               )),
      tabPanel('DE',
               sidebarLayout(
                 sidebarPanel(
                   fileInput('DEmatrix', 'Please upload your sample information below', accept = 'csv',
                             buttonLabel = 'Browse...', placeholder = 'Example_count_data.csv')
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel('Summary'),
                     tabPanel('Table'),
                     tabPanel('Plots')
                   )
                 )
               )),
      tabPanel('Other',
               sidebarLayout(
                 sidebarPanel(
                   fileInput('otherinfo', 'Please upload your sample information below', accept = 'csv',
                             buttonLabel = 'Browse...', placeholder = 'Example_count_data.csv')
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel('Summary'),
                     tabPanel('Table'),
                     tabPanel('Plots')
                   )
                 )
               ))
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    # Function to load in metadata
    load_metadata <- reactive({
      req(input$sampleinfo)  # Ensure the file input is available before proceeding
      metad <- read.csv(input$sampleinfo$datapath, row.names = 1)
      metad <- as_tibble(metad)
      metad$pmi <- as.numeric(metad$pmi)
      colnames(metad)[colnames(metad) == "age.of.death"] <- "age of death"
      colnames(metad)[colnames(metad) == "mrna.seq.reads"] <- "number of mrna reads"
      return(metad)
    })
    # Function to manipulate metadata for summary view
    summarize_column <- function(col) {
      if (is.numeric(col)) {
        summary <- paste0(round(mean(col, na.rm = TRUE)), " ± ", round(sd(col, na.rm = TRUE)))
      } else {
        summary <- paste(unique(col), collapse = ", ")
      }
      return(summary)
    }
    # Create summary table 
    summary_table <- function(data) {
      data_filtered <- data %>%
        select(-sample_id, -geo_accession, -source_name_ch1, -organism_ch1, -tissue)
      summary_tibble <- tibble(
        'Column Name' = colnames(data_filtered),
        'Type' = sapply(data_filtered, class),
        'Mean (sd) or Distinct Values' = sapply(data_filtered, summarize_column)
      )
      return(summary_tibble)
    }
    output$summarytable <- renderTable({
      metadata <- load_metadata()  # Load the data reactively
      summary_table(metadata)
    })
    
    # Create radiobuttons based on column names when file is uploaded
    output$samplecolumns <- renderUI({
      req(input$sampleinfo)
      checkboxGroupInput('samplecolumns', 'Please select desired sample columns',
                         choices = colnames(load_metadata()), 
                         selected = colnames(load_metadata()))
    })
    # Sample table output according to checkboxes
    summary_choice <- function(data, cols){
      data <- data %>%
        select(cols)
      return(data)
    }
    output$sampletable <- DT::renderDataTable({
      metadata <- load_metadata() 
      selected_columns <- input$samplecolumns 
      selected_data <- summary_choice(metadata, selected_columns)
      DT::datatable(selected_data, options = list(pageLength = 5, autoWidth = TRUE))
    })
    
    #Selecting columns and groups for the plots
    output$samplexcolumn <- renderUI({
      req(input$sampleinfo)
      metadata <- load_metadata() %>%
        select(pmi, 'age of death', rin, 'number of mrna reads')
      radioButtons('samplexcolumn', 'Please select desired column for the x-axis',
                   choices = colnames(metadata))
    })
    # Only show this table if the Plots tab is selected
    observe({
      # Show checkboxes for table selection in "Table" sub-tab
      if (input$tabs == "Samples" && input$sampletabs == "Table") {
        shinyjs::show("samplecolumns")  # Show checkboxes
        shinyjs::hide("samplexcolumn")  # Hide x-axis radio buttons
      } 
      # Show radio buttons for x and y axis selection in "Plots" sub-tab
      else if (input$tabs == "Samples" && input$sampletabs == "Plots") {
        shinyjs::show("samplexcolumn")  # Show x-axis radio buttons
        shinyjs::hide("samplecolumns")  # Hide checkboxes
      } 
      else {
        # Hide all when not in the appropriate tab
        shinyjs::hide("samplecolumns") 
        shinyjs::hide("samplexcolumn")
      }
    })
    # Density plot 
    output$densityplots <- renderPlot({
      req(input$sampleinfo, input$samplexcolumn)
      metadata <- load_metadata()
      # Get the selected x-axis column
      x_col <- input$samplexcolumn
      # Ensure the selected column is numeric
      if (!is.numeric(metadata[[x_col]])) {
        metadata[[x_col]] <- as.numeric(metadata[[x_col]])
      }
      ggplot(metadata, aes(x = !!sym(x_col), fill = diagnosis)) +
        geom_density(alpha = 0.5) +  # Use alpha for transparency
        theme_minimal() +
        labs(title = paste("Density Plot of", x_col, "Grouped by Diagnosis"),
             x = x_col,
             y = "Density")
    })
    
    # Function for loading count matrix file 
    load_countdata <- reactive({
      req(input$countdata)
      counts <- read.csv(input$countdata, header=TRUE, row.names =1)
      return(counts)
    })
    # Sliders for count matrix
    output$variance_slider <- renderUI({
      req(input$countdata)
      sliderInput('variance', 'Please adjust variance to filter out genes', min = 0, max = 100, 10)
    })
    output$nonzero_slider <- renderUI({
      req(input$countdata)
      counts <- loadcountdata()
      sliderInput('nonzero', 'Please select genes to inlcude with at least X samples that are non-zero', 
                  min = 0, max = ncol(counts), round(ncol(counts)/2))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
