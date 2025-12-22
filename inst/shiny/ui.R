library(shiny)

ui <- fluidPage(

  titlePanel("CellTFusion Shiny App"),

  sidebarLayout(
    sidebarPanel(

      fileInput(
        "counts",
        "Upload raw counts (CSV)",
        accept = c(".csv")
      ),

      checkboxInput(
        "normalized",
        "Counts already normalized",
        value = TRUE
      ),

      hr(),

      h4("Step-by-step execution"),

      actionButton("run_norm", "1. Normalize counts"),
      actionButton("run_deconv", "2. Deconvolution"),
      actionButton("run_tfs", "3. TF activity"),
      actionButton("run_network", "4. TF network"),
      actionButton("run_all", "Run EVERYTHING"),

      hr(),

      verbatimTextOutput("status")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Deconvolution", tableOutput("deconv_out")),
        tabPanel("TF activity", tableOutput("tfs_out")),
        tabPanel("Network", verbatimTextOutput("network_out"))
      )
    )
  )
)
