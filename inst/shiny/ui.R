library(shiny)

ui <- fluidPage(

  titlePanel("CellTFusion"),

  sidebarLayout(

    sidebarPanel(

      ## =========================
      ## INPUT DATA
      ## =========================
      h4("Input data"),

      fileInput(
        "counts",
        "Upload raw counts (CSV)",
        accept = ".csv"
      ),

      fileInput(
        "coldata",
        "Upload sample metadata (coldata, CSV)",
        accept = ".csv"
      ),

      checkboxInput(
        "normalized",
        "Normalize counts (TPM + log)",
        value = TRUE
      ),

      hr(),

      ## =========================
      ## DECONVOLUTION
      ## =========================
      h4("Deconvolution"),

      checkboxGroupInput(
        "deconv_methods",
        "Deconvolution methods",
        choices = c(
          "Quantiseq",
          "DWLS",
          "Epidish",
          "DeconRNASeq",
          "CibersortX"
        ),
        selected = c("Quantiseq", "DWLS", "Epidish", "DeconRNASeq", "CibersortX")
      ),

      fileInput(
        "deconv_file",
        "Upload deconvolution matrix (optional)",
        accept = ".csv"
      ),

      hr(),

      ## =========================
      ## CIBERSORTX
      ## =========================
      h4("CibersortX credentials"),

      textInput(
        "cbsx_mail",
        "CibersortX email",
        value = ""
      ),

      passwordInput(
        "cbsx_token",
        "CibersortX token"
      ),

      textInput(
        "file_name",
        "Output file name prefix",
        value = ""
      ),

      hr(),

      ## =========================
      ## TF ACTIVITY
      ## =========================
      h4("TF activity"),

      selectInput(
        "TF.collection",
        "TF collection",
        choices = c("CollecTRI", "DoRothEA"),
        selected = "CollecTRI"
      ),

      numericInput(
        "min_targets_size",
        "Minimum TF targets",
        value = 10,
        min = 1
      ),

      checkboxInput(
        "tfs.pruned",
        "Prune TFs",
        value = FALSE
      ),

      hr(),

      ## =========================
      ## NETWORK
      ## =========================
      h4("TF network"),

      numericInput(
        "minMod",
        "Minimum module size",
        value = 10,
        min = 2
      ),

      sliderInput(
        "corr_mod",
        "Module correlation threshold",
        min = 0,
        max = 1,
        value = 0.9
      ),

      hr(),

      ## =========================
      ## TRAIT & BATCH
      ## =========================
      h4("Trait & batch"),

      textInput(
        "trait",
        "Trait column name (in coldata)",
        value = ""
      ),

      textInput(
        "trait_positive",
        "Positive trait value",
        value = ""
      ),

      checkboxInput(
        "batch",
        "Use batch correction",
        value = FALSE
      ),

      textInput(
        "batch_id",
        "Batch column name (in coldata)",
        value = ""
      ),

      hr(),

      ## =========================
      ## CORRELATION & GROUPS
      ## =========================
      h4("Correlation & cell groups"),

      selectInput(
        "corr_type",
        "Correlation type",
        choices = c("spearman", "pearson"),
        selected = "spearman"
      ),

      textInput(
        "cells_extra",
        "Extra cell types (comma-separated)",
        value = ""
      ),

      numericInput(
        "pval",
        "P-value threshold",
        value = 0.05,
        min = 0,
        max = 1,
        step = 0.01
      ),

      hr(),

      ## =========================
      ## RUN
      ## =========================
      actionButton(
        "run_all",
        "Run CellTFusion",
        class = "btn-primary"
      )
    ),

    mainPanel(

      h4("Status"),
      verbatimTextOutput("status"),

      hr(),

      h4("Deconvolution (preview)"),
      tableOutput("deconv_out"),

      hr(),

      h4("TF activity (preview)"),
      tableOutput("tfs_out"),

      hr(),

      h4("TF network"),
      verbatimTextOutput("network_out")
    )
  )
)
