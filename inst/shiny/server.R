library(shiny)
library(CellTFusion)

server <- function(input, output, session) {

  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0) y else x
  }

  read_csv_matrix <- function(path) {
    as.data.frame(read.csv(path, row.names = 1, check.names = FALSE))
  }

  parse_nullable <- function(x) {
    if (is.null(x) || !nzchar(trimws(x))) {
      return(NULL)
    }
    trimws(x)
  }

  parse_cells_extra <- function(x) {
    if (is.null(x) || !nzchar(trimws(x))) {
      return(NULL)
    }
    vals <- trimws(unlist(strsplit(x, ",")))
    vals[vals != ""]
  }

  as_preview_df <- function(x, n = 8) {
    if (is.null(x)) {
      return(data.frame())
    }

    if (is.list(x) && !is.data.frame(x)) {
      if (length(x) == 0) {
        return(data.frame())
      }
      return(as_preview_df(x[[1]], n = n))
    }

    obj <- as.data.frame(x, check.names = FALSE)
    utils::head(obj, n)
  }

  first_frame <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    if (is.data.frame(x) || is.matrix(x)) {
      return(as.data.frame(x, check.names = FALSE))
    }
    if (is.list(x)) {
      if (length(x) == 0) {
        return(NULL)
      }
      return(first_frame(x[[1]]))
    }
    NULL
  }

  numeric_df <- function(df) {
    if (is.null(df) || !is.data.frame(df)) {
      return(data.frame())
    }
    num_cols <- vapply(df, is.numeric, logical(1))
    out <- df[, num_cols, drop = FALSE]
    as.data.frame(out, check.names = FALSE)
  }

  orient_samples_features <- function(df, sample_names = NULL) {
    if (is.null(df)) {
      return(data.frame())
    }

    m <- numeric_df(as.data.frame(df, check.names = FALSE))
    if (nrow(m) == 0 || ncol(m) == 0) {
      return(m)
    }

    samples <- sample_names %||% character(0)
    if (length(samples) > 0) {
      row_match <- sum(rownames(m) %in% samples)
      col_match <- sum(colnames(m) %in% samples)
      if (col_match > row_match) {
        m <- as.data.frame(t(m), check.names = FALSE)
      }
    }
    m
  }

  flatten_cell_groups <- function(cell_groups) {
    if (is.null(cell_groups) || !is.list(cell_groups)) {
      return(data.frame())
    }

    out <- list()
    idx <- 1
    for (module_name in names(cell_groups)) {
      module_obj <- cell_groups[[module_name]]
      if (!is.list(module_obj) || length(module_obj) < 2) {
        next
      }
      group_scores <- module_obj[[1]]
      group_comp <- module_obj[[2]]
      if (!is.list(group_scores) || !is.list(group_comp)) {
        next
      }
      group_names <- names(group_scores)
      if (is.null(group_names)) {
        group_names <- paste0("group_", seq_along(group_scores))
      }
      for (j in seq_along(group_scores)) {
        composition_size <- if (length(group_comp) >= j && !is.null(group_comp[[j]])) {
          length(group_comp[[j]])
        } else {
          NA_integer_
        }
        out[[idx]] <- data.frame(
          Module = module_name,
          Group = group_names[[j]],
          Features_in_group = composition_size,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }

    if (length(out) == 0) {
      return(data.frame())
    }
    do.call(rbind, out)
  }

  rv <- reactiveValues(
    counts = NULL,
    coldata = NULL,
    deconv = NULL,
    res = NULL,
    run_time = NULL,
    error = NULL
  )

  observeEvent(input$counts, {
    rv$counts <- read_csv_matrix(input$counts$datapath)
    rv$error <- NULL
    showNotification("Raw counts uploaded.", type = "message")
  })

  observeEvent(input$coldata, {
    rv$coldata <- read_csv_matrix(input$coldata$datapath)
    rv$error <- NULL
    showNotification("Metadata uploaded.", type = "message")
  })

  observeEvent(input$deconv_file, {
    rv$deconv <- read_csv_matrix(input$deconv_file$datapath)
    rv$error <- NULL
    showNotification("Deconvolution matrix uploaded.", type = "message")
  })

  observeEvent(input$load_demo, {
    data("raw.counts.tuto", package = "CellTFusion", envir = environment())
    data("traitdata.tuto", package = "CellTFusion", envir = environment())
    rv$counts <- get("raw.counts.tuto", envir = environment())
    rv$coldata <- get("traitdata.tuto", envir = environment())
    rv$deconv <- NULL
    rv$error <- NULL
    showNotification("Tutorial data loaded from CellTFusion.", type = "message")
  })

  observeEvent(input$run_all, {
    req(rv$counts)

    methods <- input$deconv_methods
    if (length(methods) == 0 && is.null(rv$deconv)) {
      showNotification("Select at least one deconvolution method or upload a deconvolution matrix.", type = "error")
      return(NULL)
    }

    file_name <- parse_nullable(input$file_name)
    coldata <- rv$coldata
    if (!is.null(coldata)) {
      coldata <- as.data.frame(coldata, check.names = FALSE)
    }

    args <- list(
      raw.counts = as.data.frame(rv$counts, check.names = FALSE),
      deconv = rv$deconv,
      normalized = isTRUE(input$normalized),
      coldata = coldata,
      batch = isTRUE(input$batch),
      batch_id = parse_nullable(input$batch_id),
      deconv_methods = methods,
      cbsx.mail = parse_nullable(input$cbsx_mail),
      cbsx.token = parse_nullable(input$cbsx_token),
      file_name = file_name,
      task = input$task,
      contrast = parse_nullable(input$contrast),
      ref_level = parse_nullable(input$ref_level),
      TF.collection = input$TF.collection,
      min_targets_size = input$min_targets_size,
      minMod = input$minMod,
      corr_mod = input$corr_mod,
      corr = input$corr,
      corr_type = input$corr_type,
      cells_extra = parse_cells_extra(input$cells_extra),
      pval = input$pval,
      high_corr_groups = input$high_corr_groups,
      return = TRUE,
      verbose = TRUE
    )

    if (isTRUE(input$batch) && (is.null(args$batch_id) || is.null(args$coldata))) {
      showNotification("Batch mode requires coldata and a valid batch column name.", type = "error")
      return(NULL)
    }

    if (identical(input$task, "supervised") && (is.null(args$contrast) || is.null(args$ref_level) || is.null(args$coldata))) {
      showNotification("Supervised mode requires coldata, contrast, and reference level.", type = "error")
      return(NULL)
    }

    rv$error <- NULL
    start_time <- Sys.time()

    withProgress(message = "Running CellTFusion pipeline", value = 0, {
      incProgress(0.15, detail = "Preparing inputs")
      run_res <- tryCatch({
        incProgress(0.65, detail = "Computing deconvolution, TFs, and modules")
        do.call(CellTFusion::CellTFusion, args)
      }, error = function(e) {
        rv$error <- conditionMessage(e)
        NULL
      })
      incProgress(0.2, detail = "Finalizing outputs")
      rv$res <- run_res
      rv$run_time <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
    })

    if (is.null(rv$res) && !is.null(rv$error)) {
      showNotification(paste("CellTFusion failed:", rv$error), type = "error", duration = 12)
    } else {
      showNotification("CellTFusion run completed.", type = "message")
    }
  })

  output$status <- renderPrint({
    list(
      counts_loaded = !is.null(rv$counts),
      coldata_loaded = !is.null(rv$coldata),
      deconv_uploaded = !is.null(rv$deconv),
      results_ready = !is.null(rv$res),
      run_time_seconds = rv$run_time,
      last_error = rv$error
    )
  })

  output$deconv_out <- renderTable({
    req(rv$res)
    as_preview_df(rv$res$Deconvolution)
  }, rownames = TRUE)

  output$tfs_out <- renderTable({
    req(rv$res)
    as_preview_df(rv$res$TFs_matrix)
  }, rownames = TRUE)

  output$cell_groups_out <- renderTable({
    req(rv$res)
    utils::head(flatten_cell_groups(rv$res$Cell_groups), 12)
  }, rownames = FALSE)

  output$tme_out <- renderTable({
    req(rv$res)
    as_preview_df(rv$res$TME_states)
  }, rownames = TRUE)

  output$network_out <- renderPrint({
    req(rv$res)
    str(rv$res$TF_network, max.level = 2)
  })

  deconv_data <- reactive({
    req(rv$res)
    first_frame(rv$res$Deconvolution)
  })

  deconv_pca <- reactive({
    df <- deconv_data()
    m <- numeric_df(df)

    validate(
      need(nrow(m) > 1, "Need at least 2 samples for PCA."),
      need(ncol(m) > 1, "Need at least 2 numeric deconvolution features for PCA.")
    )

    pca <- stats::prcomp(m, center = TRUE, scale. = TRUE)
    pcs <- as.data.frame(pca$x[, 1:2, drop = FALSE], check.names = FALSE)
    pcs$Sample <- rownames(m)

    if (!is.null(rv$coldata) && nrow(rv$coldata) > 0) {
      cdf <- as.data.frame(rv$coldata, check.names = FALSE)
      group_col <- names(cdf)[1]
      idx <- match(pcs$Sample, rownames(cdf))
      pcs$Group <- as.character(cdf[[group_col]][idx])
      pcs$Group[is.na(pcs$Group)] <- "Unknown"
    } else {
      pcs$Group <- "Samples"
    }

    var_exp <- summary(pca)$importance[2, 1:2] * 100
    list(data = pcs, var = var_exp)
  })

  output$deconv_pca_plot <- renderPlot({
    obj <- deconv_pca()
    pcs <- obj$data
    vx <- round(obj$var[[1]], 1)
    vy <- round(obj$var[[2]], 1)

    ggplot2::ggplot(pcs, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
      ggplot2::geom_point(size = 3, alpha = 0.85) +
      ggplot2::geom_text(ggplot2::aes(label = Sample), vjust = -0.7, size = 3, show.legend = FALSE) +
      ggplot2::labs(
        x = paste0("PC1 (", vx, "%)"),
        y = paste0("PC2 (", vy, "%)"),
        color = "Group"
      ) +
      ggplot2::theme_minimal(base_size = 12)
  })

  output$deconv_interaction <- renderPrint({
    obj <- deconv_pca()
    pcs <- obj$data

    clicked <- if (!is.null(input$deconv_pca_click)) {
      shiny::nearPoints(pcs, input$deconv_pca_click, xvar = "PC1", yvar = "PC2", maxpoints = 1)
    } else {
      data.frame()
    }

    hovered <- if (!is.null(input$deconv_pca_hover)) {
      shiny::nearPoints(pcs, input$deconv_pca_hover, xvar = "PC1", yvar = "PC2", maxpoints = 1)
    } else {
      data.frame()
    }

    brushed <- if (!is.null(input$deconv_pca_brush)) {
      shiny::brushedPoints(pcs, input$deconv_pca_brush, xvar = "PC1", yvar = "PC2")
    } else {
      data.frame()
    }

    list(
      clicked_sample = if (nrow(clicked) > 0) clicked$Sample[[1]] else NULL,
      hovered_sample = if (nrow(hovered) > 0) hovered$Sample[[1]] else NULL,
      brushed_samples = if (nrow(brushed) > 0) brushed$Sample else character(0)
    )
  })

  tf_feature_data <- reactive({
    req(rv$res)
    sample_names <- rownames(deconv_data()) %||% rownames(rv$counts)
    tf <- first_frame(rv$res$TFs_matrix)
    orient_samples_features(tf, sample_names = sample_names)
  })

  output$tf_feature_selector <- renderUI({
    df <- tf_feature_data()
    validate(need(ncol(df) > 0, "No numeric TF activity features available."))
    shiny::selectInput("tf_feature", "Select TF feature", choices = colnames(df), selected = colnames(df)[1])
  })

  tf_plot_data <- reactive({
    req(input$tf_feature)
    df <- tf_feature_data()
    validate(need(input$tf_feature %in% names(df), "Selected TF feature not found."))

    out <- data.frame(
      Sample = rownames(df),
      Value = as.numeric(df[[input$tf_feature]]),
      stringsAsFactors = FALSE
    )

    if (!is.null(rv$coldata) && nrow(rv$coldata) > 0) {
      cdf <- as.data.frame(rv$coldata, check.names = FALSE)
      group_col <- names(cdf)[1]
      idx <- match(out$Sample, rownames(cdf))
      out$Group <- as.character(cdf[[group_col]][idx])
      out$Group[is.na(out$Group)] <- "Unknown"
    } else {
      out$Group <- "Samples"
    }

    out
  })

  output$tf_feature_plot <- renderPlot({
    dat <- tf_plot_data()

    ggplot2::ggplot(dat, ggplot2::aes(x = Group, y = Value, color = Group)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25, width = 0.5) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.85, size = 2) +
      ggplot2::labs(x = "Group", y = input$tf_feature, color = "Group") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
  })

  output$tf_feature_stats <- renderTable({
    dat <- tf_plot_data()
    out <- stats::aggregate(Value ~ Group, data = dat, FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE)))
    data.frame(
      Group = out$Group,
      Mean = round(out$Value[, "mean"], 4),
      SD = round(out$Value[, "sd"], 4),
      stringsAsFactors = FALSE
    )
  })

  output$tf_feature_interaction <- renderPrint({
    dat <- tf_plot_data()

    clicked <- if (!is.null(input$tf_feature_click)) {
      shiny::nearPoints(dat, input$tf_feature_click, xvar = "Group", yvar = "Value", maxpoints = 1)
    } else {
      data.frame()
    }

    hovered <- if (!is.null(input$tf_feature_hover)) {
      shiny::nearPoints(dat, input$tf_feature_hover, xvar = "Group", yvar = "Value", maxpoints = 1)
    } else {
      data.frame()
    }

    list(
      selected_tf = input$tf_feature,
      clicked_sample = if (nrow(clicked) > 0) clicked$Sample[[1]] else NULL,
      hovered_sample = if (nrow(hovered) > 0) hovered$Sample[[1]] else NULL
    )
  })

  tme_long <- reactive({
    req(rv$res)
    tme <- first_frame(rv$res$TME_states)
    tme <- orient_samples_features(tme, sample_names = rownames(deconv_data()))
    validate(
      need(nrow(tme) > 0, "No TME state matrix found."),
      need(ncol(tme) > 0, "No numeric TME state features found.")
    )

    samples <- rownames(tme)
    states <- colnames(tme)

    long <- data.frame(
      Sample = rep(samples, times = length(states)),
      State = rep(states, each = length(samples)),
      Score = as.numeric(as.matrix(tme)),
      stringsAsFactors = FALSE
    )

    long$State <- factor(long$State, levels = states)
    long$Sample <- factor(long$Sample, levels = rev(samples))
    long$X <- as.numeric(long$State)
    long$Y <- as.numeric(long$Sample)
    long
  })

  output$tme_heatmap_plot <- renderPlot({
    long <- tme_long()
    ggplot2::ggplot(long, ggplot2::aes(x = State, y = Sample, fill = Score)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#2563eb", mid = "#f8fafc", high = "#dc2626", midpoint = 0) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
        panel.grid = ggplot2::element_blank()
      ) +
      ggplot2::labs(x = "TME state", y = "Sample", fill = "Score")
  })

  output$tme_interaction <- renderPrint({
    long <- tme_long()

    clicked <- if (!is.null(input$tme_heatmap_click)) {
      shiny::nearPoints(long, input$tme_heatmap_click, xvar = "X", yvar = "Y", maxpoints = 1)
    } else {
      data.frame()
    }

    hovered <- if (!is.null(input$tme_heatmap_hover)) {
      shiny::nearPoints(long, input$tme_heatmap_hover, xvar = "X", yvar = "Y", maxpoints = 1)
    } else {
      data.frame()
    }

    list(
      clicked_cell = if (nrow(clicked) > 0) {
        list(sample = as.character(clicked$Sample[[1]]), state = as.character(clicked$State[[1]]), score = round(clicked$Score[[1]], 4))
      } else {
        NULL
      },
      hovered_cell = if (nrow(hovered) > 0) {
        list(sample = as.character(hovered$Sample[[1]]), state = as.character(hovered$State[[1]]), score = round(hovered$Score[[1]], 4))
      } else {
        NULL
      }
    )
  })

  citation_bundle <- reactive({
    txt <- tryCatch(capture.output(utils::citation(package = "CellTFusion")), error = function(e) NULL)
    bib <- tryCatch(capture.output(utils::toBibtex(utils::citation(package = "CellTFusion"))), error = function(e) NULL)

    if (is.null(txt) || length(txt) == 0) {
      txt <- c(
        "CellTFusion: Integration of immune-cell deconvolution and TF networks for TME states.",
        "Please cite the package and associated publication when available."
      )
    }
    if (is.null(bib) || length(bib) == 0) {
      bib <- txt
    }

    list(text = txt, bib = bib)
  })

  output$citation_txt <- renderPrint({
    cat(paste(citation_bundle()$text, collapse = "\n"))
  })

  output$authors_info <- renderUI({
    desc <- tryCatch(utils::packageDescription("CellTFusion"), error = function(e) NULL)

    if (is.null(desc)) {
      return(HTML("<p>Package metadata not found in this R session.</p>"))
    }

    maintainer <- if (!is.null(desc$Maintainer)) desc$Maintainer else "Not available"
    authors <- if (!is.null(desc$Author)) desc$Author else "Not available"

    HTML(paste0(
      "<p><strong>Maintainer:</strong> ", maintainer, "</p>",
      "<p><strong>Authors:</strong> ", authors, "</p>",
      "<p><strong>Package:</strong> ", desc$Package, " (v", desc$Version, ")</p>",
      "<p>This Shiny interface follows the workflow described in the README and vignette, including deconvolution, TF activity, WTCNA modules, cell groups, and TME states.</p>"
    ))
  })

  download_data <- function(id, filename, extractor) {
    output[[id]] <- downloadHandler(
      filename = function() {
        paste0(filename, "_", format(Sys.Date(), "%Y%m%d"), ".csv")
      },
      content = function(file) {
        req(rv$res)
        obj <- extractor(rv$res)
        if (is.null(obj)) {
          write.csv(data.frame(), file, row.names = FALSE)
        } else {
          write.csv(as.data.frame(obj, check.names = FALSE), file, row.names = TRUE)
        }
      }
    )
  }

  download_data("download_deconv", "deconvolution", function(x) x$Deconvolution)
  download_data("download_tfs", "tf_activity", function(x) x$TFs_matrix)
  download_data("download_tme", "tme_states", function(x) x$TME_states)

  output$download_groups <- downloadHandler(
    filename = function() {
      paste0("cell_groups_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(rv$res)
      write.csv(flatten_cell_groups(rv$res$Cell_groups), file, row.names = FALSE)
    }
  )

  output$download_rds <- downloadHandler(
    filename = function() {
      paste0("CellTFusion_results_", format(Sys.Date(), "%Y%m%d"), ".rds")
    },
    content = function(file) {
      req(rv$res)
      saveRDS(rv$res, file = file)
    }
  )

  output$download_citation <- downloadHandler(
    filename = function() {
      paste0("CellTFusion_citation_", format(Sys.Date(), "%Y%m%d"), ".bib")
    },
    content = function(file) {
      writeLines(citation_bundle()$bib, con = file)
    }
  )
}
