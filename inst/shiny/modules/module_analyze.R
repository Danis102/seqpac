# module_analyze.R - Analysis & Visualization module

analyzeUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    tabsetPanel(
      id = ns("analysis_tabs"),
      
      # 1. PCA Tab
      tabPanel("Principal Component Analysis (PCA)",
        br(),
        sidebarLayout(
          sidebarPanel(
            h4("PCA Settings"),
            uiOutput(ns("pca_group_ui")),
            checkboxInput(ns("pca_labels"), "Show Sample Labels", value = TRUE),
            selectInput(ns("pca_style"), "PCA Target", choices = c("Samples" = "samples", "Sequences" = "anno")),
            uiOutput(ns("pca_anno_target_ui")),
            actionButton(ns("run_pca_btn"), "Run PCA", class = "btn-primary w-100")
          ),
          mainPanel(
            div(class = "card",
              div(class = "card-header", "PCA Plot"),
              div(class = "card-body",
                plotOutput(ns("pca_plot"), height = "500px")
              )
            )
          )
        )
      ),
      
      # 2. DESeq2 Tab
      tabPanel("Differential Expression (DESeq2)",
        br(),
        sidebarLayout(
          sidebarPanel(
            h4("DESeq2 Model Settings"),
            p("Define the experimental design model formula (e.g. ~ stage)."),
            uiOutput(ns("deseq_factor_ui")),
            textInput(ns("deseq_formula"), "Design Formula", value = "~ stage"),
            p(style = "font-size: 0.85em; color: #64748b;", 
              "Hint: Specify factors present in the Pheno table, prefixed by '~'."),
            actionButton(ns("run_deseq_btn"), "Run DESeq2 Analysis", class = "btn-primary w-100")
          ),
          mainPanel(
            div(class = "card",
              div(class = "card-header", "DESeq2 Results Table"),
              div(class = "card-body",
                DT::dataTableOutput(ns("deseq_results_tbl"))
              )
            )
          )
        )
      ),
      
      # 3. Size & Nucleotide Bias Tab
      tabPanel("Size & Nucleotide Bias",
        br(),
        sidebarLayout(
          sidebarPanel(
            h4("Plot Settings"),
            uiOutput(ns("bias_anno_col_ui")),
            selectInput(ns("plot_type"), "Select Plot Type", 
                        choices = c("Nucleotide Bias" = "nbias", "Size Distribution" = "sizedist")),
            conditionalPanel(
              condition = sprintf("input['%s'] == 'nbias'", ns("plot_type")),
              numericInput(ns("bias_position"), "Nucleotide Position", value = 1, min = 1, max = 50)
            ),
            actionButton(ns("run_plots_btn"), "Generate Plots", class = "btn-primary w-100")
          ),
          mainPanel(
            div(class = "card",
              div(class = "card-header", "Distribution & Bias Plots"),
              div(class = "card-body",
                plotOutput(ns("bias_sizedist_plot"), height = "550px")
              )
            )
          )
        )
      ),
      
      # 4. Composition Tab
      tabPanel("Composition (Bar & Pie)",
        br(),
        sidebarLayout(
          sidebarPanel(
            h4("Composition Settings"),
            uiOutput(ns("comp_anno_col_ui")),
            selectInput(ns("comp_style"), "Chart Style", choices = c("Stacked Bar" = "bar", "Pie Chart" = "pie")),
            actionButton(ns("run_comp_btn"), "Generate Composition Chart", class = "btn-primary w-100")
          ),
          mainPanel(
            div(class = "card",
              div(class = "card-header", "Composition Chart"),
              div(class = "card-body",
                plotOutput(ns("comp_plot"), height = "500px")
              )
            )
          )
        )
      )
    )
  )
}

analyzeServer <- function(id, pac_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ------------------ Dynamic UI Helpers ------------------
    output$pca_group_ui <- renderUI({
      pac_obj <- pac_reactive()
      if (is.null(pac_obj)) return(NULL)
      cols <- colnames(pheno(pac_obj))
      selectInput(ns("pca_group"), "Color by Group:", choices = cols, selected = cols[1])
    })
    
    output$pca_anno_target_ui <- renderUI({
      pac_obj <- pac_reactive()
      req(pac_obj)
      req(input$pca_style == "anno")
      cols <- colnames(anno(pac_obj))
      selectInput(ns("pca_anno_col"), "Annotation Column:", choices = cols)
    })
    
    output$deseq_factor_ui <- renderUI({
      pac_obj <- pac_reactive()
      if (is.null(pac_obj)) return(NULL)
      cols <- colnames(pheno(pac_obj))
      # Auto-fill design formula default
      if (length(cols) > 0) {
        updateTextInput(session, "deseq_formula", value = paste0("~ ", cols[1]))
      }
      selectInput(ns("deseq_factor"), "Primary Factor (for plot targets):", choices = cols)
    })
    
    output$bias_anno_col_ui <- renderUI({
      pac_obj <- pac_reactive()
      if (is.null(pac_obj)) return(NULL)
      cols <- colnames(anno(pac_obj))
      selectInput(ns("bias_anno_col"), "Annotation Classification Column:", choices = cols)
    })
    
    output$comp_anno_col_ui <- renderUI({
      pac_obj <- pac_reactive()
      if (is.null(pac_obj)) return(NULL)
      cols <- colnames(anno(pac_obj))
      selectInput(ns("comp_anno_col"), "Annotation Column:", choices = cols)
    })
    
    # ------------------ PCA Logic ------------------
    pca_plot_val <- reactiveVal(NULL)
    
    observeEvent(input$run_pca_btn, {
      pac_obj <- pac_reactive()
      if (is.null(pac_obj)) {
        showNotification("Please load and filter data first.", type = "error")
        return()
      }
      
      showNotification("Running PCA...", type = "message")
      tryCatch({
        lbl <- if (input$pca_labels) pheno(pac_obj)$Sample_ID else NULL
        
        # Determine targets
        pheno_tgt <- list(input$pca_group)
        
        if (input$pca_style == "anno") {
          anno_tgt <- list(input$pca_anno_col)
          pca_res <- PAC_pca(pac_obj, style = "anno", anno_target = anno_tgt, label = lbl)
        } else {
          pca_res <- PAC_pca(pac_obj, pheno_target = pheno_tgt, label = lbl)
        }
        
        # Save plots
        # FactoMineR/ggplot combination
        if (!is.null(pca_res$graphs)) {
          # Combine the list of 3 PCA ggplots into a nice grid
          p <- cowplot::plot_grid(plotlist = pca_res$graphs, ncol = 2, nrow = 2)
          pca_plot_val(p)
        } else {
          showNotification("PCA ran, but no plots were returned.", type = "warning")
        }
      }, error = function(e) {
        showNotification(paste("PCA failed:", e$message), type = "error")
      })
    })
    
    output$pca_plot <- renderPlot({
      p <- pca_plot_val()
      if (is.null(p)) {
        # Return a placeholder plot
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Click 'Run PCA' to generate plot") + 
          theme_void()
      } else {
        p
      }
    })
    
    # ------------------ DESeq2 Logic ------------------
    deseq_res_val <- reactiveVal(NULL)
    
    observeEvent(input$run_deseq_btn, {
      pac_obj <- pac_reactive()
      if (is.null(pac_obj)) {
        showNotification("Please load and filter data first.", type = "error")
        return()
      }
      
      formula_str <- input$deseq_formula
      showNotification("Running DESeq2. This may take a minute...", type = "message")
      
      tryCatch({
        # Run DESeq2 wrapper
        model_formula <- as.formula(formula_str)
        
        # Run
        de_res <- PAC_deseq(
          pac_obj, 
          model = model_formula, 
          threads = 1,
          pheno_target = list(input$deseq_factor)
        )
        
        if (!is.null(de_res$result)) {
          deseq_res_val(de_res$result)
          showNotification("DESeq2 finished successfully!", type = "default")
        } else {
          showNotification("DESeq2 completed, but no results table was generated.", type = "warning")
        }
      }, error = function(e) {
        showNotification(paste("DESeq2 run failed:", e$message), type = "error")
      })
    })
    
    output$deseq_results_tbl <- DT::renderDataTable({
      res <- deseq_res_val()
      if (is.null(res)) {
        return(DT::datatable(data.frame(Message = "Run DESeq2 to populate table.")))
      }
      DT::datatable(res, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    # ------------------ Bias & Size Distribution Logic ------------------
    bias_plot_val <- reactiveVal(NULL)
    
    observeEvent(input$run_plots_btn, {
      pac_obj <- pac_reactive()
      if (is.null(pac_obj)) {
        showNotification("Please load and filter data first.", type = "error")
        return()
      }
      
      showNotification("Generating plots...", type = "message")
      tryCatch({
        anno_col <- input$bias_anno_col
        
        if (input$plot_type == "nbias") {
          pos <- input$bias_position
          res <- PAC_nbias(pac_obj, position = pos, anno_target = list(anno_col))
          # Cowplot grid output
          p <- cowplot::plot_grid(plotlist = res$Histograms)
          bias_plot_val(p)
        } else {
          res <- PAC_sizedist(pac_obj, anno_target = list(anno_col))
          p <- cowplot::plot_grid(plotlist = res$Histograms)
          bias_plot_val(p)
        }
      }, error = function(e) {
        showNotification(paste("Plotting failed:", e$message), type = "error")
      })
    })
    
    output$bias_sizedist_plot <- renderPlot({
      p <- bias_plot_val()
      if (is.null(p)) {
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Click 'Generate Plots' to display results.") + 
          theme_void()
      } else {
        p
      }
    })
    
    # ------------------ Composition Logic ------------------
    comp_plot_val <- reactiveVal(NULL)
    
    observeEvent(input$run_comp_btn, {
      pac_obj <- pac_reactive()
      if (is.null(pac_obj)) {
        showNotification("Please load and filter data first.", type = "error")
        return()
      }
      
      showNotification("Generating composition chart...", type = "message")
      tryCatch({
        anno_col <- input$comp_anno_col
        
        if (input$comp_style == "bar") {
          p <- PAC_stackbar(pac_obj, anno_target = list(anno_col))
          comp_plot_val(p)
        } else {
          p <- PAC_pie(pac_obj, anno_target = list(anno_col))
          comp_plot_val(p)
        }
      }, error = function(e) {
        showNotification(paste("Composition plot failed:", e$message), type = "error")
      })
    })
    
    output$comp_plot <- renderPlot({
      p <- comp_plot_val()
      if (is.null(p)) {
        ggplot() + 
          annotate("text", x = 0.5, y = 0.5, label = "Click 'Generate' to display composition.") + 
          theme_void()
      } else {
        p
      }
    })
  })
}
