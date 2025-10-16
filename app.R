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
        uiOutput("download_plot_ui")
      ),
      mainPanel(
        h3("geneCOCOA Differential Analysis: AMI vs. Control"),
        withSpinner(plotOutput("cocoa_plot", height = "800px"))
      )
    )
  }
)

# =============================================================================
# SERVER LOGIC
# =============================================================================
server <- function(input, output, session) {
  
  # Reactive to hold the currently generated plot
  current_plot <- reactiveVal()
  
  output$cocoa_plot <- renderPlot({
    req(precomputed_results, input$gene_choice)
    
    gene_data <- precomputed_results[[input$gene_choice]]
    
    # This should not happen if the prep script worked, but it's a good safeguard.
    if (is.null(gene_data) || !is.data.frame(gene_data$ami) || !is.data.frame(gene_data$control)) {
      return(ggplot() + theme_void() + labs(title = "Error: Invalid data structure for the selected gene."))
    }
    
    res_ami <- gene_data$ami
    res_control <- gene_data$control
    title_gene <- input$gene_choice
    
    ami_plot_data <- res_ami %>% mutate(neg_log_p = -log10(adj_p_value), Condition = "AMI") %>% select(geneset, neg_log_p, Condition)
    control_plot_data <- res_control %>% mutate(neg_log_p = -log10(adj_p_value) * -1, Condition = "Control") %>% select(geneset, neg_log_p, Condition)
    
    combined_data <- bind_rows(ami_plot_data, control_plot_data) %>% filter(is.finite(neg_log_p))
    sig_genesets <- res_ami %>% full_join(res_control, by = "geneset", suffix = c("_ami", "_control")) %>% filter(adj_p_value_ami < 0.05 | adj_p_value_control < 0.05) %>% pull(geneset)
    
    if (length(sig_genesets) == 0) {
      p <- ggplot() + theme_void() + labs(title = paste("Analysis complete: No significant co-regulation found for", title_gene))
      current_plot(NULL) # Clear download button
      return(p)
    }
    
    plot_data <- combined_data %>%
      filter(geneset %in% sig_genesets) %>%
      mutate(Term = gsub("HALLMARK_", "", geneset, fixed = TRUE))
    
    order <- plot_data %>% group_by(Term) %>% summarise(value = sum(neg_log_p)) %>% arrange(value) %>% pull(Term)
    plot_data$Term <- factor(plot_data$Term, levels = order)
    
    p <- ggplot(plot_data, aes(x = neg_log_p, y = Term, fill = Condition)) +
      geom_col(width = 0.8) +
      geom_vline(xintercept = 0, color = "grey20") +
      scale_fill_manual(values = c("AMI" = "#b30000", "Control" = "grey70")) +
      labs(x = "Significance [-log10(Adjusted P-value)]", y = "Hallmark Gene Set", title = paste("Co-regulation Analysis for:", title_gene), subtitle = "AMI (right, red) vs. Control (left, grey)", fill = "Condition") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom", panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
    
    current_plot(p) # Save the plot for download
    return(p)
  })
  
  # UI for the download button
  output$download_plot_ui <- renderUI({
    req(current_plot())
    downloadButton("download_plot", "Download Plot")
  })
  
  # Logic to download the plot
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

