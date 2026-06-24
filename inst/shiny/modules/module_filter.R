# module_filter.R - Preprocessing: Filter & Normalize module

filterUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        h3("Preprocessing Settings"),
        p("Configure filters to reduce noise and select normalization methods."),
        
        # Filtering Inputs
        div(class = "card",
          div(class = "card-header", "Filter Settings"),
          div(class = "card-body",
            sliderInput(ns("size_range"), "Nucleotide Size Range", 
                        min = 10, max = 100, value = c(20, 30), step = 1),
            numericInput(ns("count_threshold"), "Min Counts (Threshold)", value = 5, min = 0),
            sliderInput(ns("coverage_threshold"), "Min Coverage (% of samples)", 
                        min = 0, max = 100, value = 20, step = 5),
            p(style = "font-size: 0.85em; color: #64748b;", 
              "Default: Keep sequences of 20-30 nt with >= 5 counts in at least 20% of samples.")
          )
        ),
        
        # Normalization Inputs
        div(class = "card",
          div(class = "card-header", "Normalization Methods"),
          div(class = "card-body",
            checkboxGroupInput(ns("norm_methods"), "Select Normalizations to Run:",
                               choices = c("Counts per Million (CPM)" = "cpm",
                                           "Variance Stabilizing (VST)" = "vst",
                                           "Regularized Log (RLOG)" = "rlog"),
                               selected = c("cpm", "vst"))
          )
        ),
        
        actionButton(ns("apply_filter_btn"), "Apply Filter & Normalize", class = "btn-primary w-100")
      ),
      
      mainPanel(
        div(class = "card",
          div(class = "card-header", "Filter Results & Stats"),
          div(class = "card-body",
            fluidRow(
              column(6,
                h4("Sequence Counts Summary"),
                verbatimTextOutput(ns("filter_summary")),
                br(),
                uiOutput(ns("normalization_status"))
              ),
              column(6,
                h4("Sequence Length Distribution"),
                plotOutput(ns("length_dist_plot"), height = "300px")
              )
            )
          )
        )
      )
    )
  )
}

filterServer <- function(id, raw_pac_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Store processed PAC locally
    processed_pac <- reactiveVal(NULL)
    
    # Enable/disable controls based on PAC availability
    observe({
      if (is.null(raw_pac_reactive())) {
        shinyjs::disable("apply_filter_btn")
      } else {
        shinyjs::enable("apply_filter_btn")
      }
    })
    
    # Filter and Normalize action
    observeEvent(input$apply_filter_btn, {
      pac_obj <- raw_pac_reactive()
      req(pac_obj)
      
      showNotification("Filtering PAC object...", type = "message")
      
      tryCatch({
        # 1. Filter
        filtered <- PAC_filter(
          pac_obj, 
          nucleotide_range = input$size_range, 
          threshold = input$count_threshold, 
          coverage = input$coverage_threshold,
          norm = "counts"
        )
        
        # 2. Normalize
        norms <- input$norm_methods
        if (length(norms) > 0) {
          showNotification("Running normalizations...", type = "message")
          for (nm in norms) {
            filtered <- PAC_norm(filtered, norm = nm)
          }
        }
        
        # Save to reactive value
        processed_pac(filtered)
        showNotification("Filtering and normalization completed successfully!", type = "default")
        
      }, error = function(e) {
        showNotification(paste("Error during preprocessing:", e$message), type = "error")
      })
    })
    
    # Populate stats & graphs
    output$filter_summary <- renderText({
      raw_pac <- raw_pac_reactive()
      filt_pac <- processed_pac()
      
      if (is.null(raw_pac)) {
        return("Please load a PAC object first on the Load tab.")
      }
      
      out <- paste0(
        "--- Original Data ---\n",
        sprintf("Sequences: %d\n", nrow(counts(raw_pac))),
        sprintf("Total reads sum: %d\n\n", sum(counts(raw_pac)))
      )
      
      if (!is.null(filt_pac)) {
        pct_seqs <- (nrow(counts(filt_pac)) / nrow(counts(raw_pac))) * 100
        pct_reads <- (sum(counts(filt_pac)) / sum(counts(raw_pac))) * 100
        out <- paste0(
          out,
          "--- Filtered Data ---\n",
          sprintf("Sequences remaining: %d (%.2f%% of original)\n", nrow(counts(filt_pac)), pct_seqs),
          sprintf("Reads remaining: %d (%.2f%% of original)\n", sum(counts(filt_pac)), pct_reads)
        )
      } else {
        out <- paste0(out, "Filter not applied yet. Click 'Apply Filter & Normalize' to start.")
      }
      return(out)
    })
    
    output$normalization_status <- renderUI({
      filt_pac <- processed_pac()
      req(filt_pac)
      
      # Check which normalizations are in the list
      norm_slots <- names(filt_pac@norm)
      
      tagList(
        h5("Available Normalizations:"),
        if (length(norm_slots) > 0) {
          tags$ul(
            lapply(norm_slots, function(n) tags$li(strong(n)))
          )
        } else {
          p(style = "color: #b91c1c;", "No normalized tables available.")
        }
      )
    })
    
    output$length_dist_plot <- renderPlot({
      # Plot length distribution of sequences
      pac_obj <- processed_pac()
      if (is.null(pac_obj)) {
        pac_obj <- raw_pac_reactive()
      }
      req(pac_obj)
      
      # Get sequences
      seqs <- rownames(counts(pac_obj))
      lengths <- nchar(seqs)
      
      df <- data.frame(Length = lengths)
      
      ggplot(df, aes(x = Length)) +
        geom_bar(fill = "#6366f1", alpha = 0.8, color = "#4f46e5") +
        theme_minimal(base_family = "Inter") +
        labs(
          title = "Read Length Distribution",
          x = "Sequence Length (nt)",
          y = "Count of Unique Sequences"
        ) +
        theme(
          plot.title = element_text(family = "Outfit", face = "bold", size = 14),
          panel.grid.minor = element_blank()
        )
    })
    
    # Return processed PAC
    return(processed_pac)
  })
}
