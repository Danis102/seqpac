# ui.R - Main user interface shell for the seqpac Shiny App

ui <- page_navbar(
  id = "tabs",
  title = "Seqpac Dashboard",
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    primary = "#4f46e5",
    secondary = "#475569",
    base_font = font_google("Inter"),
    heading_font = font_google("Outfit")
  ),
  
  # Custom Header containing shinyjs and custom style sheet
  header = tagList(
    shinyjs::useShinyjs(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    )
  ),
  
  # Tabs
  nav_panel("Home", launchUI("launch")),
  nav_panel("Load / Create PAC", loadUI("load")),
  nav_panel("Filter & Normalize", filterUI("filter")),
  nav_panel("Annotation Explorer", annotateUI("annotate")),
  nav_panel("Post-Filtering Analysis", analyzeUI("analyze"))
)
