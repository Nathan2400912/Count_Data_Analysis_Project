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
library(pheatmap)
library(colourpicker)


options(shiny.maxRequestSize = 100 * 1024^2)
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
                   uiOutput('nonzero_slider'),
                   actionButton('apply_filters', 'Apply Filters'), 
                   uiOutput('scatterplotaxis'),
                   uiOutput("pc_selection_ui"),
                   actionButton("plot_pca", "Plot PCA")
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel('Summary',
                              tableOutput('count_summary')),
                     tabPanel('Scatter Plot',
                              plotOutput('scatter_plots')),
                     tabPanel('Heatmap',
                              plotOutput('heatmap')),
                     tabPanel('PCA',
                              plotOutput('pca_plot'))
                   )
                 )
               )),
      tabPanel('DE',
               sidebarLayout(
                 sidebarPanel(
                   fileInput('DEresults', 'Please upload your DE results below', accept = '.csv',
                             buttonLabel = 'Browse...', placeholder = 'Example_DE_results.csv'),
                   uiOutput('DE_x'),
                   uiOutput('DE_y'),
                   colourInput('basecolor', 'Base point color', value = '#22577A'),
                   colourInput('highlightcolor','Highlight point color','#FFCF56'),
                   sliderInput('magnitude','Select the magnitude of the p adjusted coloring:',
                               min = -35, max = 0, value = -17),
                   actionButton('plot_button', 'Plot', style = "color: black; background-color: #D7A0DE; border-color: black; padding: 10px 20px; font-size: 15px; width: 400px",
                                icon = icon('chart-area'))
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel('Table',
                              DT::dataTableOutput('resultstable')),
                     tabPanel('Plots',
                              plotOutput('volcano', height = '600px'))
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
      counts <- read.csv(input$countdata$datapath, header=TRUE, row.names =1)
      return(counts)
    })
    
    # Sliders for count matrix
    output$variance_slider <- renderUI({
      req(input$countdata)
      counts <- load_countdata()  # Load the count data
      sliderInput('variance', 'Please adjust variance to filter out genes',
                  min = 0, max = 100, value = 50)
    })
    # Create non-zero slider UI once count data is loaded
    output$nonzero_slider <- renderUI({
      req(input$countdata)
      
      counts <- load_countdata()  # Load the count data
      non_zero_count <- rowSums(counts > 0)  # Count how many non-zero samples per gene
      max_nonzero <- max(non_zero_count)  # Find the maximum non-zero count
      
      sliderInput('nonzero', 'Please select genes to include with at least X samples that are non-zero',
                  min = 0, max = max_nonzero, value = round(max_nonzero / 2))
    })
    
    # Y-axis choice for the scatterplot
    output$scatterplotaxis <- renderUI({
      req(input$countdata, input$variance, input$nonzero)
      radioButtons('scattery', 'Please select desired filter for the y-axis',
                   choices = c('variance', 'non-zero'),
                   selected = input$scattery %||% 'variance')
    })
    # Create a table describing samples and genes that passed filters
    observeEvent(input$apply_filters, {
      req(input$countdata, input$variance, input$nonzero)
      
      counts <- load_countdata()  # Load the count data
      
      # Calculate variance per gene
      variance <- apply(counts, 1, var, na.rm = TRUE)
      variance_percentile_threshold <- quantile(variance, probs = input$variance / 100, na.rm = TRUE)
      
      # Filter genes based on variance
      genes_passing_variance <- variance >= variance_percentile_threshold
      
      # Count how many non-zero samples per gene
      non_zero_count <- rowSums(counts > 0)
      genes_passing_nonzero <- non_zero_count >= input$nonzero
      
      # Apply both filters
      filtered_genes <- genes_passing_variance & genes_passing_nonzero
      
      # Get the filtered count data
      filtered_counts <- counts[filtered_genes, , drop = FALSE]
      
      # Create a summary table for the counts
      count_table <- tibble(
        'Number of samples' = ncol(counts),
        'Number of genes' = nrow(counts),
        'Number of genes passing filters' = sum(filtered_genes),
        'Number of genes not passing filters' = sum(!filtered_genes),
        'Percent of genes passing filters' = round(sum(filtered_genes) / nrow(counts) * 100, 2),
        'Percent of genes not passing filters' = round(sum(!filtered_genes) / nrow(counts) * 100, 2)
      )
      
      # Output the count table
      output$count_summary <- renderTable({
        count_table
      })
      # Scatter plot for count data passing the filters
      output$scatter_plots <- renderPlot({
        # Get median count per gene
        median_per_gene <- apply(counts, 1, median, na.rm = TRUE)
        
        # Create the plot data
        plot_data <- tibble(
          median_count = median_per_gene,
          variance = variance,
          non_zero_count = non_zero_count,
          genes_passing_variance = genes_passing_variance,
          genes_passing_nonzero = genes_passing_nonzero,
          genes_passing = filtered_genes
        )
        
        # Create scatter plot based on the selected y-axis choice
        if (input$scattery == 'variance') {
          ggplot(plot_data, aes(x = median_count, y = variance, color = genes_passing)) +
            geom_point() +
            scale_x_log10() +  # Apply log scale for x-axis (Median Count)
            scale_y_log10() +  # Optional: Apply log scale for y-axis (Variance)
            labs(x = "Median Count", y = "Variance", 
                 title = "Median Count vs Variance (Filtered Genes)") +
            theme_minimal() +
            scale_color_manual(values = c("lightgrey", "darkred"))  # Lighter color for non-passing genes
        } else {
          ggplot(plot_data, aes(x = median_count, y = non_zero_count, color = genes_passing)) +
            geom_point() +
            scale_x_log10() +
            scale_y_continuous(limits = c(0, NA)) +
            labs(x = "Median Count", y = "Number of Non-Zero Samples", 
                 title = "Median Count vs Number of Non-Zero Samples (Filtered Genes)") +
            theme_minimal() +
            scale_color_manual(values = c("lightblue", "darkblue"))  # Lighter color for non-passing genes
        }
      })
      
      # Heatmap code
      output$heatmap <- renderPlot({
        req(filtered_counts)
        log_counts <- log1p(filtered_counts) #Do i need to log transform?
        pheatmap(log_counts,
                 cluster_rows = TRUE,  
                 cluster_cols = TRUE,
                 scale = "row", # Should i scale or not?
                 show_rownames = FALSE, 
                 show_colnames = TRUE,
                 color = colorRampPalette(c("blue", "white", "red"))(50))
      })
    })
    # Let user select PC to plot 
    pca_result <- reactive({
      req(input$countdata)  # Ensure count data is available
      counts <- load_countdata()  # Load count data
      counts <- as.matrix(counts)  # Ensure it's a matrix
      # Perform PCA on transposed count data
      pca_res <- prcomp(t(counts), center = TRUE, scale = FALSE) # scale or Not?
      return(pca_res)
    })
    
    # Dynamically update the PC selection UI based on the number of PCs
    output$pc_selection_ui <- renderUI({
      req(pca_result())  # Ensure PCA result is available
      num_pcs <- ncol(pca_result()$x)  # Number of PCs
      pc_choices <- paste0("PC", 1:num_pcs)  # Create dynamic choices
      
      tagList(
        selectInput("pc_x", "Select X-axis PC:", choices = pc_choices, selected = "PC1"),
        selectInput("pc_y", "Select Y-axis PC:", choices = pc_choices, selected = "PC2")
      )
    })
    
    # Plot the selected PCA plot only when the button is clicked
    
    # vst/rlog or not?
    
    observeEvent(input$plot_pca, {
      req(pca_result(), input$pc_x, input$pc_y)  # Ensure PCA result and selected PCs are available
      
      pc_x <- as.numeric(gsub("PC", "", input$pc_x))  # Extract selected PC number
      pc_y <- as.numeric(gsub("PC", "", input$pc_y))  # Extract selected PC number
      
      pca_res <- pca_result()  # Get PCA result
      pca_data <- as.data.frame(pca_res$x)  # Get PCA scores for each sample
      metadata <- load_metadata()
      pca_data$Condition <- metadata$diagnosis  # Assuming 'diagnosis' column in metadata
      # Generate the PCA plot
      output$pca_plot <- renderPlot({
        ggplot(pca_data, aes_string(x = paste0("PC", pc_x), y = paste0("PC", pc_y), color = "Condition")) +
          geom_point(size = 3) +
          labs(
            title = paste("PCA: PC", pc_x, "vs PC", pc_y),
            x = paste("PC", pc_x, " - ", round(100 * summary(pca_res)$importance[2, pc_x]), "% variance"),
            y = paste("PC", pc_y, " - ", round(100 * summary(pca_res)$importance[2, pc_y]), "% variance")
          ) +
          theme_minimal() +
          scale_color_manual(values = c("red", "blue"))  # Customize color palette if needed
      })
    })
    
    
    # DE Section
    
    # Load in DE results 
    load_de <- reactive({
      req(input$DEresults)
      de <- read.csv(input$DEresults$datapath, header=TRUE, row.names =1)
      return(de)
    })
    
    # Table of DE results
    output$resultstable <- DT::renderDataTable({
      de <- load_de()
      DT::datatable(de, options = list(pageLength = 5, autoWidth = TRUE))
    })
    
    # Render UI output for the radiobuttons to select columns to plot
    output$DE_x <- renderUI({
      req(input$DEresults)
      de <- load_de() %>%
        select(-1)
      radioButtons('xbutton', 'Choose the column for the x-axis', colnames(de))
    })
    output$DE_y <- renderUI({
      req(input$DEresults)
      de <- load_de() %>%
        select(-1)
      radioButtons('ybutton', 'Choose the column for the y-axis', colnames(de))
    })
    
    # Plotting parameters
    plot_params <- eventReactive(input$plot_button, {
      list(
        x_name = input$xbutton,
        y_name = input$ybutton,
        slider = input$magnitude,
        color1 = input$basecolor,
        color2 = input$highlightcolor
      )
    })
    # DE Plot Function
    volcano_plot <-
      function(dataf, x_name, y_name, slider, color1, color2) {
        req(input$DEresults)
        volcano <- dataf %>%
          mutate(colour = !!sym(y_name) < (1*10^slider))
        volplot <- ggplot(volcano, aes(x = !!sym(x_name), y = -log10(!!sym(y_name)), color = colour)) +
          geom_point() +
          scale_color_manual(values = c(color1, color2),
                             name = paste('padj < 1 x 10^', plot_params()$slider, sep = ""),
                             labels = c('FALSE', 'TRUE')) +
          labs(x = x_name, y = paste("-log10(", y_name, ")", sep = "")) +
          theme_minimal() +
          theme(legend.position = 'bottom')
        return(volplot)
      }
    # Render the plot
    output$volcano <- renderPlot(volcano_plot(load_de(),plot_params()$x_name, plot_params()$y_name,
                                              plot_params()$slider, plot_params()$color1, plot_params()$color2))
}

# Run the application 
shinyApp(ui = ui, server = server)
