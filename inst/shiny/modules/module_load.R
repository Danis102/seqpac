# module_load.R - Load/Create PAC module

loadUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        h3("Data Selection"),
        p("Choose a method to import or create a PAC object:"),
        
        tabsetPanel(
          id = ns("load_tabs"),
          type = "pills",
          
          # Option 1: Built-in Example
          tabPanel("Example Data",
            br(),
            p("Load the Drosophila small RNA dataset included in the package (filtered & annotated)."),
            actionButton(ns("load_example_btn"), "Load Drosophila Dataset", class = "btn-primary w-100")
          ),
          
          # Option 2: Upload PAC RData
          tabPanel("Import PAC",
            br(),
            p("Upload an existing PAC object (.RData or .rds file)."),
            fileInput(ns("pac_file"), "Select PAC File", accept = c(".RData", ".rds", ".Rdata")),
            actionButton(ns("load_pac_btn"), "Load File", class = "btn-primary w-100")
          ),
          
          # Option 3: Create from Raw FASTQ
          tabPanel("Create from FASTQ",
            br(),
            p("Generate a new PAC object from FASTQ reads and a phenotypic CSV."),
            
            # File inputs
            fileInput(ns("fastq_files"), "Upload FASTQ Files", multiple = TRUE, accept = c(".fastq", ".fq", ".gz")),
            
            # Local Docker scan path if available
            htmlOutput(ns("docker_data_ui")),
            
            fileInput(ns("pheno_file"), "Upload Phenotype CSV (with Sample_ID)", accept = c(".csv")),
            
            # Trimming settings
            selectInput(ns("trim_type"), "Trimming Method", 
                        choices = c("Seqpac Internal" = "seqpac", "None (Already Trimmed)" = "none")),
            selectInput(ns("adapter_parse"), "Adapter / Protocol", 
                        choices = c("NEBNext Small RNA" = "default_neb", "Illumina" = "default_illumina")),
            numericInput(ns("evidence_exp"), "Evidence: Min samples", value = 2, min = 1),
            numericInput(ns("evidence_samp"), "Evidence: Min counts per sample", value = 1, min = 1),
            
            br(),
            actionButton(ns("create_pac_btn"), "Generate PAC Object", class = "btn-primary w-100")
          )
        )
      ),
      
      mainPanel(
        div(class = "card",
          div(class = "card-header", "PAC Object Overview"),
          div(class = "card-body",
            verbatimTextOutput(ns("pac_status")),
            uiOutput(ns("pac_preview_tabs"))
          )
        )
      )
    )
  )
}

loadServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive value to store the active PAC object
    current_pac <- reactiveVal(NULL)
    
    # Check if `/data` exists and contains fastq files
    output$docker_data_ui <- renderUI({
      docker_dir <- "/data"
      if (dir.exists(docker_dir)) {
        files <- list.files(path = docker_dir, pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE)
        if (length(files) > 0) {
          tagList(
            div(class = "alert alert-info",
              p(strong("Docker Volume Detected:"), sprintf("Found %d FASTQ files in %s.", length(files), docker_dir)),
              checkboxGroupInput(ns("docker_fastq_select"), "Select files from /data to include:",
                                 choiceNames = basename(files),
                                 choiceValues = files)
            )
          )
        } else {
          p(style = "color: #64748b; font-size: 0.9em;", "Note: Empty /data volume folder detected.")
        }
      } else {
        NULL
      }
    })
    
    # 1. Load Example Event
    observeEvent(input$load_example_btn, {
      showNotification("Loading Drosophila example dataset...", type = "message")
      tryCatch({
        env <- new.env()
        load(system.file("extdata", "drosophila_sRNA_pac_filt_anno.Rdata", package = "seqpac", mustWork = TRUE), envir = env)
        if (exists("pac", envir = env)) {
          current_pac(env$pac)
          showNotification("Drosophila example loaded successfully!", type = "default")
        } else {
          showNotification("Error: 'pac' object not found in the RData file.", type = "error")
        }
      }, error = function(e) {
        showNotification(paste("Error loading example data:", e$message), type = "error")
      })
    })
    
    # 2. Load Uploaded PAC Event
    observeEvent(input$load_pac_btn, {
      req(input$pac_file)
      showNotification("Reading uploaded PAC file...", type = "message")
      
      file_path <- input$pac_file$datapath
      ext <- tools::file_ext(input$pac_file$name)
      
      tryCatch({
        if (tolower(ext) == "rds") {
          pac_obj <- readRDS(file_path)
          if (inherits(pac_obj, "PAC") || is.list(pac_obj)) {
            current_pac(as.PAC(pac_obj))
            showNotification("PAC loaded successfully!", type = "default")
          } else {
            showNotification("Loaded object is not a PAC class object.", type = "error")
          }
        } else {
          env <- new.env()
          load(file_path, envir = env)
          # Find any PAC or list object that looks like a PAC
          objs <- ls(envir = env)
          found <- FALSE
          for (o in objs) {
            val <- get(o, envir = env)
            if (inherits(val, "PAC") || (is.list(val) && all(c("Pheno", "Anno", "Counts") %in% names(val)))) {
              current_pac(as.PAC(val))
              showNotification(sprintf("Loaded PAC object '%s' successfully!", o), type = "default")
              found <- TRUE
              break
            }
          }
          if (!found) {
            showNotification("Could not find a valid PAC object in the loaded RData.", type = "error")
          }
        }
      }, error = function(e) {
        showNotification(paste("Error reading PAC file:", e$message), type = "error")
      })
    })
    
    # 3. Create from FASTQ Event
    observeEvent(input$create_pac_btn, {
      # We need either uploaded fastq files or selected docker fastq files
      fastq_paths <- NULL
      if (!is.null(input$fastq_files)) {
        # Copy to clean temp files retaining their names for sample IDs
        fastq_paths <- input$fastq_files$datapath
        # Rename them to original filenames
        temp_dir <- file.path(tempdir(), "uploaded_fastq")
        dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
        new_paths <- file.path(temp_dir, input$fastq_files$name)
        file.copy(fastq_paths, new_paths, overwrite = TRUE)
        fastq_paths <- new_paths
      }
      
      # Combine with docker path selections if any
      if (!is.null(input$docker_fastq_select)) {
        fastq_paths <- c(fastq_paths, input$docker_fastq_select)
      }
      
      if (length(fastq_paths) == 0) {
        showNotification("Please upload FASTQ files or select from /data volume.", type = "error")
        return()
      }
      
      req(input$pheno_file)
      
      showModal(modalDialog(
        title = "Creating PAC Object...",
        "Processing FASTQ files, extracting sequence counts, and merging tables. This may take several minutes...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      tryCatch({
        # 1. read pheno
        pheno_df <- read.csv(input$pheno_file$datapath, stringsAsFactors = FALSE)
        if (!"Sample_ID" %in% colnames(pheno_df)) {
          removeModal()
          showNotification("Phenotype CSV must contain a 'Sample_ID' column matching FASTQ file prefixes.", type = "error")
          return()
        }
        
        # 2. make counts
        trim_val <- if (input$trim_type == "none") NULL else "seqpac"
        
        count_list <- make_counts(
          input = fastq_paths, 
          plot = FALSE, 
          trimming = trim_val, 
          parse = input$adapter_parse,
          threads = 1,
          evidence = c(experiment = input$evidence_exp, sample = input$evidence_samp)
        )
        
        # 3. make pheno
        pheno_obj <- make_pheno(
          pheno = pheno_df,
          progress_report = count_list$progress_report,
          counts = count_list$counts
        )
        
        # 4. make PAC
        pac_obj <- make_PAC(pheno = pheno_obj, counts = count_list$counts)
        
        current_pac(pac_obj)
        removeModal()
        showNotification("PAC generated successfully!", type = "default")
        
      }, error = function(e) {
        removeModal()
        showNotification(paste("Failed to generate PAC object:", e$message), type = "error")
      })
    })
    
    # Status display
    output$pac_status <- renderText({
      pac_obj <- current_pac()
      if (is.null(pac_obj)) {
        "No PAC object loaded yet. Please select an option on the left panel."
      } else {
        paste0(
          "--- PAC Object Active ---\n",
          sprintf("Class: S4 PAC\n"),
          sprintf("Number of Sequences (Rows): %d\n", nrow(counts(pac_obj))),
          sprintf("Number of Samples (Columns): %d\n", ncol(counts(pac_obj))),
          sprintf("Annotations columns: %s\n", paste(colnames(anno(pac_obj)), collapse = ", ")),
          sprintf("Phenotypic metadata factors: %s\n", paste(colnames(pheno(pac_obj)), collapse = ", "))
        )
      }
    })
    
    # Preview tables when data loaded
    output$pac_preview_tabs <- renderUI({
      req(current_pac())
      tabsetPanel(
        tabPanel("Pheno Table", DT::dataTableOutput(ns("tbl_pheno"))),
        tabPanel("Counts Preview (Top 100)", DT::dataTableOutput(ns("tbl_counts"))),
        tabPanel("Anno Preview (Top 100)", DT::dataTableOutput(ns("tbl_anno")))
      )
    })
    
    output$tbl_pheno <- DT::renderDataTable({
      req(current_pac())
      DT::datatable(pheno(current_pac()), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$tbl_counts <- DT::renderDataTable({
      req(current_pac())
      DT::datatable(head(counts(current_pac()), 100), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$tbl_anno <- DT::renderDataTable({
      req(current_pac())
      DT::datatable(head(anno(current_pac()), 100), options = list(pageLength = 10, scrollX = TRUE))
    })
    
    # Return the reactive PAC object
    return(current_pac)
  })
}
