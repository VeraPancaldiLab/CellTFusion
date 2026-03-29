library(shiny)

ui <- fluidPage(
  tags$head(
    tags$link(
      rel = "stylesheet",
      href = "https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@400;600;700&family=Source+Sans+3:wght@400;600&display=swap"
    ),
    tags$style(HTML(" 
      :root {
        --ink: #0f172a;
        --teal: #0d9488;
        --teal-2: #14b8a6;
        --sand: #f8fafc;
        --card: #ffffff;
        --muted: #475569;
        --gold: #f59e0b;
      }
      body {
        background: radial-gradient(circle at 10% 0%, #d1fae5 0%, #f8fafc 35%, #e2e8f0 100%);
        font-family: 'Source Sans 3', sans-serif;
        color: var(--ink);
      }
      h1, h2, h3, h4, h5 {
        font-family: 'Space Grotesk', sans-serif;
        letter-spacing: 0.2px;
      }
      .hero {
        border-radius: 16px;
        background: linear-gradient(120deg, #0f766e 0%, #115e59 40%, #164e63 100%);
        padding: 18px 24px;
        color: #f8fafc;
        margin: 10px 0 18px 0;
        box-shadow: 0 14px 32px rgba(15, 23, 42, 0.18);
      }
      .hero-logo {
        max-width: 95px;
        width: 100%;
        height: auto;
        display: block;
        border-radius: 10px;
        background: #ffffff;
        padding: 6px;
      }
      .hero-title {
        margin: 0;
        font-size: 30px;
        font-weight: 700;
      }
      .hero-subtitle {
        margin: 4px 0 0 0;
        opacity: 0.95;
        font-size: 16px;
      }
      .control-card, .main-card {
        background: var(--card);
        border: 1px solid #dbe4ef;
        border-radius: 14px;
        padding: 14px;
        margin-bottom: 12px;
        box-shadow: 0 8px 20px rgba(15, 23, 42, 0.06);
      }
      .section-title {
        font-size: 16px;
        font-weight: 700;
        color: #0f766e;
        margin-bottom: 8px;
      }
      .btn-run {
        background: linear-gradient(90deg, #0d9488, #0f766e);
        border: none;
        color: #fff;
        font-weight: 700;
      }
      .btn-demo {
        background: #f8fafc;
        border: 1px solid #94a3b8;
        color: #0f172a;
        font-weight: 600;
      }
      .status-note {
        color: var(--muted);
      }
      .nav-tabs {
        margin-top: 8px;
      }
      .tab-content {
        margin-top: 12px;
      }
      .citation-box {
        background: #f8fafc;
        border-left: 4px solid var(--gold);
        border-radius: 10px;
        padding: 12px;
      }
      @media (max-width: 900px) {
        .hero-title {
          font-size: 25px;
        }
      }
    "))
  ),

  fluidRow(
    column(
      width = 12,
      div(
        class = "hero",
        fluidRow(
          column(2, tags$img(src = "logo.png", class = "hero-logo", alt = "CellTFusion logo")),
          column(
            10,
            tags$h1("CellTFusion Shiny Studio", class = "hero-title"),
            tags$p(
              "Integrated immune deconvolution and TF-network profiling for tumor microenvironment cell states.",
              class = "hero-subtitle"
            )
          )
        )
      )
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 4,

      div(
        class = "control-card",
        div(class = "section-title", "Data Inputs"),
        fileInput("counts", "Upload raw counts (CSV)", accept = ".csv"),
        fileInput("coldata", "Upload metadata/coldata (CSV)", accept = ".csv"),
        fileInput("deconv_file", "Upload deconvolution matrix (optional, CSV)", accept = ".csv"),
        checkboxInput("normalized", "Normalize counts (TPM + log2)", value = TRUE),
        actionButton("load_demo", "Load tutorial data", class = "btn-demo")
      ),

      div(
        class = "control-card",
        div(class = "section-title", "Pipeline Settings"),
        selectInput("task", "Task mode", choices = c("unsupervised", "supervised"), selected = "unsupervised"),
        checkboxGroupInput(
          "deconv_methods",
          "Deconvolution methods",
          choices = c("Quantiseq", "DWLS", "Epidish", "DeconRNASeq", "CBSX"),
          selected = c("Quantiseq", "DWLS", "Epidish", "DeconRNASeq")
        ),
        textInput("cbsx_mail", "CBSX email (if CBSX selected)", value = ""),
        passwordInput("cbsx_token", "CBSX token"),
        textInput("file_name", "Output file prefix", value = "CellTFusion_Run")
      ),

      div(
        class = "control-card",
        div(class = "section-title", "Model Parameters"),
        selectInput("TF.collection", "TF collection", choices = c("CollecTRI", "Dorothea", "ARACNE"), selected = "CollecTRI"),
        numericInput("min_targets_size", "Minimum TF targets", value = 3, min = 1),
        numericInput("minMod", "Minimum module size", value = 10, min = 2),
        sliderInput("corr_mod", "Module merge correlation", min = 0, max = 1, value = 0.9, step = 0.01),
        sliderInput("corr", "Deconvolution grouping correlation", min = 0, max = 1, value = 0.7, step = 0.01),
        selectInput("corr_type", "Correlation type", choices = c("spearman", "pearson"), selected = "spearman"),
        numericInput("pval", "P-value threshold", value = 0.05, min = 0, max = 1, step = 0.01),
        sliderInput("high_corr_groups", "High-correlation cutoff", min = 0, max = 1, value = 0.8, step = 0.01),
        textInput("cells_extra", "Additional cell labels (comma separated)", value = "")
      ),

      div(
        class = "control-card",
        div(class = "section-title", "Advanced"),
        checkboxInput("batch", "Use batch-aware analysis", value = FALSE),
        textInput("batch_id", "Batch column in coldata", value = ""),
        textInput("contrast", "Supervised contrast column", value = ""),
        textInput("ref_level", "Supervised reference level", value = "")
      ),

      actionButton("run_all", "Run CellTFusion", class = "btn-run")
    ),

    mainPanel(
      width = 8,

      div(
        class = "main-card",
        tags$h4("Run Status"),
        tags$p("Use the tutorial button for a fast first run based on package example data.", class = "status-note"),
        verbatimTextOutput("status")
      ),

      tabsetPanel(
        tabPanel(
          "Overview",
          div(
            class = "main-card",
            tags$h4("CellTFusion Workflow"),
            tags$img(src = "CellTFusion_pipeline.png", style = "width:100%; border-radius:10px;"),
            tags$p(
              "This app follows the package pipeline: deconvolution, TF activity inference, WTCNA module construction, pathway scoring, cell-group computation, latent-space projection, and TME-state definition."
            )
          )
        ),
        tabPanel(
          "Results",
          div(
            class = "main-card",
            tags$h4("Interactive Deconvolution PCA"),
            tags$p("Use click, hover, and brush to inspect sample-level deconvolution patterns."),
            plotOutput(
              "deconv_pca_plot",
              height = "340px",
              click = "deconv_pca_click",
              hover = hoverOpts("deconv_pca_hover", delay = 100, delayType = "debounce"),
              brush = brushOpts("deconv_pca_brush")
            ),
            verbatimTextOutput("deconv_interaction")
          ),
          div(
            class = "main-card",
            tags$h4("Interactive TF Feature Explorer"),
            uiOutput("tf_feature_selector"),
            plotOutput(
              "tf_feature_plot",
              height = "320px",
              click = "tf_feature_click",
              hover = hoverOpts("tf_feature_hover", delay = 100, delayType = "debounce")
            ),
            tableOutput("tf_feature_stats"),
            verbatimTextOutput("tf_feature_interaction")
          ),
          div(class = "main-card", tags$h4("Deconvolution (preview)"), tableOutput("deconv_out")),
          div(class = "main-card", tags$h4("TF Activity (preview)"), tableOutput("tfs_out")),
          div(class = "main-card", tags$h4("Cell Groups (summary)"), tableOutput("cell_groups_out")),
          div(
            class = "main-card",
            tags$h4("Interactive TME States Heatmap"),
            plotOutput(
              "tme_heatmap_plot",
              height = "380px",
              click = "tme_heatmap_click",
              hover = hoverOpts("tme_heatmap_hover", delay = 100, delayType = "debounce")
            ),
            verbatimTextOutput("tme_interaction")
          ),
          div(class = "main-card", tags$h4("TME States (preview)"), tableOutput("tme_out")),
          div(class = "main-card", tags$h4("TF Network Structure"), verbatimTextOutput("network_out"))
        ),
        tabPanel(
          "Downloads",
          div(
            class = "main-card",
            tags$h4("Export Results"),
            tags$p("Download key outputs generated by CellTFusion."),
            fluidRow(
              column(6, downloadButton("download_deconv", "Download Deconvolution")),
              column(6, downloadButton("download_tfs", "Download TF Activity"))
            ),
            br(),
            fluidRow(
              column(6, downloadButton("download_groups", "Download Cell Groups")),
              column(6, downloadButton("download_tme", "Download TME States"))
            ),
            br(),
            fluidRow(
              column(6, downloadButton("download_rds", "Download Full Results (.rds)")),
              column(6, downloadButton("download_citation", "Download Citation (.bib)"))
            )
          )
        ),
        tabPanel(
          "Citation & Team",
          div(
            class = "main-card",
            tags$h4("How To Cite"),
            div(class = "citation-box", verbatimTextOutput("citation_txt"))
          ),
          div(
            class = "main-card",
            tags$h4("Authors and Maintainer"),
            htmlOutput("authors_info")
          )
        )
      )
    )
  )
)
