# module_annotate.R - Annotate module

annotateUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        h3("Sequence Annotation"),
        p("Review annotations attached to small RNA sequences in the PAC object."),
        
        uiOutput(ns("anno_col_select_ui")),
        br(),
        div(class = "card",
          div(class = "card-header", "Annotation Statistics"),
          div(class = "card-body",
            tableOutput(ns("anno_stats_tbl"))
          )
        )
      ),
      
      mainPanel(
        div(class = "card",
          div(class = "card-header", "Annotation Table (Anno)"),
          div(class = "card-body",
            p("Rownames are the actual sequences. Columns show classifications and mapping properties:"),
            DT::dataTableOutput(ns("anno_table"))
          )
        )
      )
    )
  )
}

annotateServer <- function(id, pac_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Render select input for annotation columns
    output$anno_col_select_ui <- renderUI({
      pac_obj <- pac_reactive()
      if (is.null(pac_obj)) return(NULL)
      
      cols <- colnames(anno(pac_obj))
      selectInput(ns("selected_anno_col"), "Select Annotation Column to Analyze:",
                  choices = cols, selected = if ("Biotypes_mis0" %in% cols) "Biotypes_mis0" else cols[1])
    })
    
    # Calculate stats for the selected annotation column
    output$anno_stats_tbl <- renderTable({
      pac_obj <- pac_reactive()
      req(pac_obj)
      req(input$selected_anno_col)
      
      anno_col <- anno(pac_obj)[[input$selected_anno_col]]
      if (is.null(anno_col)) return(NULL)
      
      # Frequency table
      tbl <- table(anno_col, useNA = "always")
      df <- as.data.frame(tbl)
      colnames(df) <- c("Category", "Sequence Count")
      
      # Add percentage
      total <- sum(df$`Sequence Count`)
      df$Percentage <- sprintf("%.2f%%", (df$`Sequence Count` / total) * 100)
      
      df
    }, rownames = FALSE)
    
    # Main DT table
    output$anno_table <- DT::renderDataTable({
      pac_obj <- pac_reactive()
      if (is.null(pac_obj)) {
        return(DT::datatable(data.frame(Message = "Please load data first.")))
      }
      DT::datatable(anno(pac_obj), options = list(pageLength = 10, scrollX = TRUE))
    })
  })
}
