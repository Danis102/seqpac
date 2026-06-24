# module_launch.R - Launch / Welcome page module

launchUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        div(class = "hero-panel",
          h1("Seqpac sRNA Analysis Dashboard"),
          p("An interactive platform for small RNA sequence analysis using sequence-based counting. Preserve read integrity, visualize alignments, perform differential expression, and explore small RNA composition."),
          actionButton(ns("start_btn"), "Get Started", class = "btn-primary btn-lg")
        )
      )
    ),
    
    fluidRow(
      column(3,
        div(class = "card workflow-step-card",
          div(class = "workflow-icon", icon("upload")),
          h4("1. Load & Create"),
          p("Upload raw FASTQ sequencing files and a pheno table CSV, or import a pre-saved PAC object (RData) to initialize the workflow.")
        )
      ),
      column(3,
        div(class = "card workflow-step-card",
          div(class = "workflow-icon", icon("sliders-h")),
          h4("2. Filter & Normalize"),
          p("Apply nucleotide size filters and minimum count/coverage thresholds. Normalize counts using CPM, VST, or RLOG.")
        )
      ),
      column(3,
        div(class = "card workflow-step-card",
          div(class = "workflow-icon", icon("tags")),
          h4("3. Annotate"),
          p("Inspect sequence annotations, map small RNA reads against biotype reference databases, and evaluate alignment statistics.")
        )
      ),
      column(3,
        div(class = "card workflow-step-card",
          div(class = "workflow-icon", icon("chart-bar")),
          h4("4. Analyze & Plot"),
          p("Generate PCA, run DESeq2 differential expression, inspect size distribution & nucleotide bias, and view stacked composition plots.")
        )
      )
    ),
    
    hr(),
    
    fluidRow(
      column(12, class = "text-center",
        p(style = "color: #64748b;", "Seqpac R package dashboard — Created by Daniel Nätt, Signe Isacson, Lovisa Örkenby Kämpe, Alessandro Gozzo, Anna Asratian")
      )
    )
  )
}

launchServer <- function(id, mainTabsetSession) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$start_btn, {
      updateNavbarPage(mainTabsetSession, "tabs", selected = "Load / Create PAC")
    })
  })
}
