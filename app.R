# =============================================================================
# LIBRARIES
# =============================================================================
library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(shinycssloaders)
library(tidyr)

# =============================================================================
# INITIAL SETUP
# =============================================================================
# Load the pre-computed results. The app will fail gracefully if the file is not found.
precomputed_results <- if (file.exists("precomputed_results.rds")) {
  readRDS("precomputed_results.rds")
} else {
  NULL
}

# =============================================================================
# USER INTERFACE (UI)
# =============================================================================
ui <- fluidPage(
  titlePanel("Single-Cell Gene Co-regulation Analysis (Server Version)"),
  
  # Conditional UI: Show error if data is not loaded, otherwise show the app
  if (is.null(precomputed_results)) {
    fluidRow(
      column(12,
             h3("Error: Pre-computed data not found."),
             p("Please run the `prepare_data.R` script on the server first to generate the necessary `precomputed_results.rds` file."),
             p("This script performs all the heavy calculations offline so the app can be fast.")
      )
    )
  } else {
    sidebarLayout(
      sidebarPanel(
        h4("Co-regulation Analysis"),
        p("Displaying pre-computed geneCOCOA results."),
        selectizeInput("gene_choice", 
                       "Select a Pre-Analyzed Gene:", 
                       choices = sort(names(precomputed_results)),
                       options = list(placeholder = 'Type to search...')),
        hr(),
        # The download button is rendered dynamically in the server
        uiOutput("download_plot_ui")
      ),
      mainPanel(
        h3("geneCOCOA Differential Analysis: AMI vs. Control"),
        # Use a spinner to indicate that the plot is being generated
        withSpinner(plotOutput("cocoa_plot", height = "800px"))
      )
    )
  }
)

# =============================================================================
# SERVER LOGIC
# =============================================================================
server <- function(input, output, session) {
  
  # Reactive value to store the plot object for download
  current_plot <- reactiveVal()
  
  output$cocoa_plot <- renderPlot({
    req(precomputed_results, input$gene_choice)
    
    gene_data <- precomputed_results[[input$gene_choice]]
    
    # Safeguard against malformed data for the selected gene
    if (is.null(gene_data) || !is.data.frame(gene_data$ami) || !is.data.frame(gene_data$control)) {
      p <- ggplot() +
        theme_void() +
        labs(title = "Error: Invalid data structure for the selected gene.") +
        theme(plot.title = element_text(hjust = 0.5, size = 16))
      current_plot(NULL)
      return(p)
    }
    
    # Combine and transform data for plotting
    # - Calculate -log10(p-value) for significance
    # - Make control values negative for the divergent bar plot
    plot_data <- bind_rows(
      mutate(gene_data$ami, neg_log_p = -log10(adj_p_value), Condition = "AMI"),
      mutate(gene_data$control, neg_log_p = -log10(adj_p_value) * -1, Condition = "Control")
    ) %>%
      select(geneset, neg_log_p, Condition) %>%
      filter(is.finite(neg_log_p)) # Remove non-finite values that can break plotting
    
    # Identify gene sets that are significant in at least one condition
    sig_genesets <- full_join(
      gene_data$ami, gene_data$control,
      by = "geneset",
      suffix = c("_ami", "_control")
    ) %>%
      filter(adj_p_value_ami < 0.05 | adj_p_value_control < 0.05) %>%
      pull(geneset)
    
    # If no gene sets are significant, display a clean, informative message
    if (length(sig_genesets) == 0) {
      p <- ggplot() +
        geom_text(aes(x = 0, y = 0, label = paste("No significant co-regulation pathways found for", input$gene_choice)), size = 5) +
        theme_void() +
        labs(title = paste("Co-regulation Analysis for:", input$gene_choice)) +
        theme(plot.title = element_text(hjust = 0.5, size = 16))
      current_plot(NULL) # Clear the plot for download
      return(p)
    }
    
    # Filter the main plot data to only include significant gene sets
    plot_data_filtered <- plot_data %>%
      filter(geneset %in% sig_genesets) %>%
      mutate(Term = gsub("HALLMARK_", "", geneset, fixed = TRUE)) # Clean up term names

    # Order terms by their total significance score for a more intuitive plot
    term_order <- plot_data_filtered %>%
      group_by(Term) %>%
      summarise(value = sum(neg_log_p)) %>%
      arrange(value) %>%
      pull(Term)
    
    plot_data_filtered$Term <- factor(plot_data_filtered$Term, levels = term_order)
    
    # Generate the divergent bar plot
    p <- ggplot(plot_data_filtered, aes(x = neg_log_p, y = Term, fill = Condition)) +
      geom_col(width = 0.8) +
      geom_vline(xintercept = 0, color = "grey20") +
      scale_fill_manual(values = c("AMI" = "#b30000", "Control" = "grey70")) +
      labs(
        x = "Significance [-log10(Adjusted P-value)]",
        y = "Hallmark Gene Set",
        title = paste("Co-regulation Analysis for:", input$gene_choice),
        subtitle = "AMI (right, red) vs. Control (left, grey)",
        fill = "Condition"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      )
    
    current_plot(p) # Save the plot for the download handler
    return(p)
  })
  
  # Dynamically render the download button only when a plot is available
  output$download_plot_ui <- renderUI({
    req(current_plot())
    downloadButton("download_plot", "Download Plot")
  })
  
  # Server logic for the download button
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("cocoa_plot_", input$gene_choice, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      ggsave(file, plot = current_plot(), width = 10, height = 12, device = "png")
    }
  )
}

# =============================================================================
# RUN THE APPLICATION
# =============================================================================
shinyApp(ui = ui, server = server)