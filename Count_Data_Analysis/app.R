library(shiny)
library(shinyjs)
library(DT)
library(tidyverse)
library(rlang)
library(pheatmap)
library(colourpicker)
library(shinydashboard)

options(shiny.maxRequestSize = 100 * 1024^2)

ui <- dashboardPage(
  dashboardHeader(title = "Count Data Analysis Application"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Samples", tabName = "samples", icon = icon("table")),
      menuItem("Counts", tabName = "counts", icon = icon("chart-line")),
      menuItem("DE", tabName = "DE", icon = icon("chart-bar")),
      menuItem("GSEA", tabName = "GSEA", icon = icon("search"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "samples",
              fluidRow(
                box(
                  title = "Upload Sample Information", width = 12, status = "primary",
                  fileInput('sampleinfo', 'Please upload your sample information below', accept = '.csv',
                            buttonLabel = 'Browse...', placeholder = 'Example_count_data.csv')
                ),
                box(
                  title = "Sample Data Summary", width = 12, status = "info",
                  tabsetPanel(
                    id = "sampletabs",
                    tabPanel('Summary', tableOutput('summarytable')),
                    tabPanel('Table', DT::dataTableOutput('sampletable'),
                             conditionalPanel(
                               condition = "input.sampletabs == 'Table'",
                               uiOutput('samplecolumns')
                             )),
                    tabPanel('Plots', plotOutput('densityplots'),
                             conditionalPanel(
                               condition = "input.sampletabs == 'Plots'",
                               uiOutput('samplexcolumn')))
                  )
                )
              )
      ),
      tabItem(tabName = "counts",
              fluidRow(
                box(
                  title = "Upload Count Data", width = 12, status = "primary",
                  fileInput('countdata', 'Please upload your normalized count matrix below', 
                            accept = '.csv',
                            buttonLabel = 'Browse...', 
                            placeholder = 'Example_count_data.csv'),
                  
                  # Variance and Non-zero sliders + Apply Filters button appear only in Summary, Scatter Plot, and Heatmap tabs
                  conditionalPanel(
                    condition = "input.countstabs == 'Summary' || input.countstabs == 'Scatter Plot' || input.countstabs == 'Heatmap'",
                    uiOutput('variance_slider'),
                    uiOutput('nonzero_slider'),
                    actionButton('apply_filters', 'Apply Filters', 
                                 style = "color: white; background-color: #5FAEE3; border-color: black; padding: 10px 20px; font-size: 15px; width: 200px",
                                 icon = icon('filter'))
                  ),
                  
                  # Scatter Plot axis selection appears only in Scatter Plot tab
                  conditionalPanel(
                    condition = "input.countstabs == 'Scatter Plot'", 
                    uiOutput('scatterplotaxis')
                  ),
                  
                  # PCA selection and button appear only in PCA tab
                  conditionalPanel(
                    condition = "input.countstabs == 'PCA'", 
                    uiOutput("pc_selection_ui"),
                    actionButton("plot_pca", "Plot PCA", 
                                 style = "color: white; background-color: #5FAEE3; border-color: black; padding: 10px 20px; font-size: 15px; width: 200px",
                                 icon = icon('chart-area'))
                  )
                ),
                box(
                  title = "Counts Data Summary", width = 12, status = "info",
                  tabsetPanel(id = "countstabs",
                              tabPanel('Summary', tableOutput('count_summary')),
                              tabPanel('Scatter Plot', plotOutput('scatter_plots')),
                              tabPanel('Heatmap', plotOutput('heatmap')),
                              tabPanel('PCA', plotOutput('pca_plot'))
                  )
                )
              )
      ),
      
      tabItem(tabName = "DE",
              fluidRow(
                box(
                  title = "Upload DE Results", width = 12, status = "primary",
                  fileInput('DEresults', 'Please upload your DE results below', accept = '.csv',
                            buttonLabel = 'Browse...', placeholder = 'Example_DE_results.csv'),
                  conditionalPanel(
                    condition = "input.DEtabs == 'Plots'",
                    fluidRow(  
                      column(6, uiOutput('DE_x')),
                      column(6, uiOutput('DE_y'))  
                    ),
                    fluidRow(  
                      column(6, colourInput('basecolor', 'Base point color', value = '#22577A')),
                      column(6, colourInput('highlightcolor', 'Highlight point color','#FFCF56'))  
                    ),
                    sliderInput('magnitude','Select the magnitude of the p adjusted coloring:',
                                min = -35, max = 0, value = -17),
                    actionButton('plot_button', 'Plot', style = "color: white; background-color: #5FAEE3; border-color: black; padding: 10px 20px; font-size: 15px; width: 100px",
                                 icon = icon('chart-area'))
                  )
                ),
                box(
                  title = "DE Results", width = 12, status = "info",
                  tabsetPanel(id = "DEtabs",
                    tabPanel('Table', DT::dataTableOutput('resultstable')),
                    tabPanel('Plots', plotOutput('volcano', height = '600px'))
                  )
                )
              )
      ),
      
      tabItem(tabName = "GSEA",
              fluidRow(
                box(
                  title = "Upload GSEA Results", width = 12, status = "primary",
                  fileInput('gsearesults', 'Please upload your gsea results below', accept = '.csv',
                            buttonLabel = 'Browse...', placeholder = 'Example_gsea_results.csv'),
                  sliderInput('adj_pvalue','Filter gene sets based on adjusted p-value threshold:',
                              min = -20, max = 0, value = 0),
                  conditionalPanel(
                    condition = "input.gseasubtabs == 'Results'",
                    downloadButton('download_gsea', 'Download Results', style = "color: white; background-color: #5FAEE3; border-color: black; padding: 10px 20px; font-size: 15px; width: 200px",
                                   icon = icon('download')),
                    radioButtons(inputId = "nes_filter",label = "Select pathways:",
                      choices = list("All" = "all", "Positive NES" = "positive", "Negative NES" = "negative"),
                      selected = "all"
                    )
                  )
                ),
                box(
                  title = "GSEA Results", width = 12, status = "info",
                  tabsetPanel(
                    id = "gseasubtabs",
                    tabPanel('Top Pathways', plotOutput('enrichment_plot')),
                    tabPanel('Results', DT::dataTableOutput('gsea_table')),
                    tabPanel('Scatter Plot', plotOutput('nes_plot'))
                  )
                )
              )
      )
    )
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
    zero_count <- ncol(counts) - non_zero_count
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
        zero_count = zero_count,
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
          scale_color_manual(values = c("darkred","lightgrey"))  # Lighter color for non-passing genes
      } else {
        ggplot(plot_data, aes(x = median_count, y = zero_count, color = genes_passing)) +
          geom_point() +
          scale_x_log10() +
          scale_y_continuous(limits = c(0, NA)) +
          labs(x = "Median Count", y = "Number of Zero Samples", 
               title = "Median Count vs Number of Zero Samples (Filtered Genes)") +
          theme_minimal() +
          scale_color_manual(values = c("darkblue","lightblue"))  # Lighter color for non-passing genes
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
    numeric_columns <- sapply(de, is.numeric)
    DT::datatable(de, options = list(pageLength = 5, autoWidth = TRUE)) %>%
      DT::formatSignif(columns = which(numeric_columns), digits = 3) %>%
      DT::formatStyle(columns = which(numeric_columns), 
                      target = 'cell', 
                      props = 'white-space: nowrap;')
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
  
  
  # GSEA section
  
  # Reading the GSEA csv
  load_gsea <- reactive({
    gsea <- read.csv(input$gsearesults$datapath, header = TRUE)
    return(gsea)
  })
  
  # Enrichmnet plot
  output$enrichment_plot <- renderPlot({
    req(input$gsearesults)
    fgsea_results <- load_gsea() %>%
      arrange(desc(NES)) %>%
      filter(padj <= (1 * 10^input$adj_pvalue))
    positive_pathways <- fgsea_results %>% filter(NES > 0)
    negative_pathways <- fgsea_results %>% filter(NES < 0)
    num_top <- 10
    top_pos <- positive_pathways %>% slice_head(n = num_top)
    top_neg <- negative_pathways %>% slice_head(n = num_top)
    if (nrow(positive_pathways) < num_top) {
      top_pos <- positive_pathways
    }
    if (nrow(negative_pathways) < num_top) {
      top_neg <- negative_pathways
    }
    tog <- bind_rows(top_pos, top_neg)
    tog <- tog %>%
      mutate(condition = ifelse(NES > 0, 'positive', 'negative'))
    ggplot(tog, aes(x = reorder(pathway, NES), y = NES, fill = condition)) +  
      geom_col() +                                          
      coord_flip() +
      scale_fill_manual(values = c('positive' = '#9DEEF5', 'negative' = '#FF8282')) +
      theme(axis.text.y = element_text(size = 6),
            axis.title.x = element_text(size = 8),
            plot.title = element_text(size = 10)) +
      labs(
        y = 'Normalized Enrichment Score (NES)',
        x = '',
        title = 'fgsea results for Hallmark MSigDB gene set'
      )
  })
  
  # Results Table 
  filtered_gsea <- reactive({
    req(input$gsearesults)
    results <- load_gsea() %>%
      filter(padj <= (1 * 10^input$adj_pvalue))
    # Apply NES filtering
    if (input$nes_filter == "positive") {
      results <- results %>% filter(NES > 0)
    } else if (input$nes_filter == "negative") {
      results <- results %>% filter(NES < 0)
    }
    results
  })
  # Render the table
  output$gsea_table <- DT::renderDataTable({
    DT::datatable(filtered_gsea(), options = list(pageLength = 5, autoWidth = TRUE))
  })
  # Download the file
  output$download_gsea <- downloadHandler(
    filename = function() {
      paste('gsea_results.csv')
    },
    content = function(file) {
      write.csv(filtered_gsea(), file, row.names = FALSE)
    }
  )
  
  # Code for NES scatterplot
  output$nes_plot <- renderPlot({
    req(input$gsearesults)
    results <- load_gsea() %>%
      mutate(threshold = ifelse(padj <= (1 * 10^input$adj_pvalue), 'True', 'False'))
    ggplot(results, aes(x = NES, y= -log10(padj), color=threshold)) +
      geom_point() +
      labs(x = "NES", y = "-log10(padj)", 
           title = "NES vs padj (Filtered Gene Sets)") +
      theme_minimal() +
      scale_color_manual(values = c("True" = "darkred", "False" = "lightgrey"))
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
