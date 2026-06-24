# server.R - Main server routing logic for the seqpac Shiny App

server <- function(input, output, session) {
  
  # 1. Launch Page Server
  launchServer("launch", session)
  
  # 2. Data Loading Module Server (returns reactive raw PAC object)
  raw_pac <- loadServer("load")
  
  # 3. Filtering & Normalization Module Server (returns reactive filtered PAC object)
  filtered_pac <- filterServer("filter", raw_pac)
  
  # Determine active PAC to pass down to downstream tabs
  # If filtered PAC exists, use it. Otherwise fall back to raw PAC.
  active_pac <- reactive({
    filt <- filtered_pac()
    if (!is.null(filt)) {
      return(filt)
    } else {
      return(raw_pac())
    }
  })
  
  # 4. Annotation Explorer Server
  annotateServer("annotate", active_pac)
  
  # 5. Post-Filtering Analysis Server
  analyzeServer("analyze", active_pac)
  
}
