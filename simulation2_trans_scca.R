# Simulation study for the split-sample Trans-sparse CCA estimator.
#
# Main design:
#   K = 10 is fixed;
#   K_A = |A_h| varies from 0 to 10;
#   (p,m,n) is one of
#     (200,200,200), (200,200,400), (500,500,200), (500,500,400);
#   target marginals are Toeplitz with rho = 0.3;
#   target sparse support is S = {1,3,...,19}, s = 10;
#   exact and approximately low-rank target cross-covariances are considered;
#   informative sources have joint-covariance distance 0.04 from the target;
#   # non-informative sources are generated as 0.5 * Sigma^(0) + 0.5 * Sigma_tilde^(k).
#
# Three methods:
#   1. Detected    : sources are selected by detect_trans_scca_sources();
#   2. Pooled-All  : all sources are used and transfer is forced;
#   3. Target-only : no source is used.
#
# IMPORTANT:
# For a first debugging run, use n_rep = 2 or 5.  For final experiments,
# increase n_rep to 100 or 200 after checking runtime and convergence.

rm(list = ls())
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

.find_existing_file <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Cannot find any of: ", paste(candidates, collapse = ", "))
  }
  hit[1L]
}

main_file <- .find_existing_file(c(
  "Trans-sparse_CCA-cor_CV_fullsample.R"
))
utils_file <- .find_existing_file(c(
  "trans_scca_utils.R"
))

source(main_file)
source(utils_file)

if (!requireNamespace("foreach", quietly = TRUE) ||
    !requireNamespace("doSNOW", quietly = TRUE)) {
  stop("Install packages 'foreach' and 'doSNOW' for the parallel progress bar.")
}
library(foreach)
library(doSNOW)

  simulation_config <- list(
    settings = data.frame(
      p = c(150L, 150L, 300L, 300L),
      m = c(150L, 150L, 300L, 300L),
      n = c(200L, 400L, 200L, 400L)
    ),
  
  rank_types = c("exact", "approximate"),

  r = 2L,
  s = 10L,
  K = 10L,
  K_A_values = 0:10,

  rho_target = 0.3,
  tail_decay = 2,

  h = 0.05,
  delta_informative = 0.05,
  noninformative_weight = 1,
  # Adaptive fitting: no true s_u, s_v, or h is supplied to the estimator.
  tune_initialization = FALSE,
  tune_penalties = TRUE,
  cv_grid = c(0.5, 1, 1.5, 2),
  cv_folds = 3L,

  initialization_method = "ADMM",
  
  rho_admm = 2,
  admm_tol = 1e-5,
  admm_maxit = 500L,
  admm_inner_tol = 1e-5,
  admm_inner_maxit = 50L,
  
  normalization_ridge = 1e-8,

  # Practical source detection:
  # 5-fold CV = 80% target training + 20% target validation.
  detect_repeats = 1L,
  detect_folds = 5L,
  
  # Final estimator:
  # TRUE  = full-sample practical refit;
  # FALSE = original three-way sample splitting.
  full_sample_fit = TRUE,

  n_rep = 5L,        # Change to 100/200 for the final Monte Carlo study.
  master_seed = 20260904L,

  output_dir = "simulation2_split_results"
)

.make_data_from_Z <- function(Z, chol_sigma, p) {
  W <- Z %*% chol_sigma
  list(
    X = W[, seq_len(p), drop = FALSE],
    Y = W[, p + seq_len(ncol(W) - p), drop = FALSE]
  )
}

.fit_one_split_method <- function(
    method,
    target_data,
    sources,
    informative,
    config,
    fit_seed) {

  mode <- switch(
    method,
    "Detected" = "adaptive",
    "Pooled-All" = "transfer",
    "Target-only" = "target",
    stop("Unknown method.")
  )
  
  normalization_mode <- if (method == "Pooled-All") {
    "target"
  } else {
    "truncated"
  }

  if (method == "Target-only") {
    sources_use <- list()
    informative_use <- integer(0)
  } else {
    sources_use <- sources
    informative_use <- informative
  }

  safe_trans_fit(
    trans_scca(
      target = target_data,
      sources = sources_use,
      informative = informative_use,
      r = config$r,
      s_u = NULL,
      s_v = NULL,
      full_sample = config$full_sample_fit,
      seed = fit_seed,
      initialization_method = config$initialization_method,
      rho_admm = config$rho_admm,
      admm_tol = config$admm_tol,
      admm_maxit = config$admm_maxit,
      admm_inner_tol = config$admm_inner_tol,
      admm_inner_maxit = config$admm_inner_maxit,
      tune_initialization = config$tune_initialization,
      tune_penalties = config$tune_penalties,
      cv_grid = config$cv_grid,
      cv_folds = config$cv_folds,
      normalization_ridge = config$normalization_ridge,
      aggregation_mode = mode,
      normalization_mode = normalization_mode
    )
  )
}

.result_row <- function(
    fit,
    population,
    method,
    p,
    m,
    n,
    rank_type,
    K_A,
    rep_id,
    detected_metrics = c(TPR = NA_real_, FPR = NA_real_,
                         selected_count = NA_real_)) {

  if (fit_failed(fit)) {
    return(data.frame(
      rep = rep_id,
      p = p,
      m = m,
      n = n,
      rank_type = rank_type,
      K_A = K_A,
      method = method,
      error_U = NA_real_,
      error_V = NA_real_,
      TPR = detected_metrics["TPR"],
      FPR = detected_metrics["FPR"],
      selected_count = detected_metrics["selected_count"],
      fit_failed = TRUE,
      error_message = fit$error_message,
      stringsAsFactors = FALSE
    ))
  }

  err <- evaluate_cca_fit(fit, population)

  data.frame(
    rep = rep_id,
    p = p,
    m = m,
    n = n,
    rank_type = rank_type,
    K_A = K_A,
    method = method,
    error_U = unname(err["error_U"]),
    error_V = unname(err["error_V"]),
    TPR = detected_metrics["TPR"],
    FPR = detected_metrics["FPR"],
    selected_count = detected_metrics["selected_count"],
    fit_failed = FALSE,
    error_message = "",
    stringsAsFactors = FALSE
  )
}

.run_one_rep_split <- function(
    rep_id,
    population,
    source_bank,
    p,
    m,
    n,
    rank_type,
    config) {

  K <- config$K
  d <- p + m
  rep_seed <- config$master_seed +
    100000L * rep_id +
    1000L * p +
    10L * n +
    ifelse(rank_type == "exact", 1L, 2L)

  set.seed(rep_seed)

  # Use common Gaussian innovations across K_A values so the K_A curves are
  # paired within each Monte Carlo replication.
  Z0 <- matrix(rnorm(n * d), nrow = n, ncol = d)
  Z_sources <- lapply(seq_len(K), function(k) {
    matrix(rnorm(n * d), nrow = n, ncol = d)
  })

  target_data <- .make_data_from_Z(
    Z0,
    population$chol,
    p
  )

  # Random nested ordering: the first K_A entries are informative.
  informative_order <- sample(seq_len(K), K, replace = FALSE)

  fit_seed <- rep_seed + 300000L

  # Target-only does not depend on K_A. Fit it once per replication.
  fit_target <- .fit_one_split_method(
    method = "Target-only",
    target_data = target_data,
    sources = list(),
    informative = integer(0),
    config = config,
    fit_seed = fit_seed
  )

  rows <- list()
  row_id <- 1L

  for (K_A in config$K_A_values) {
    A_true <- if (K_A == 0L) {
      integer(0)
    } else {
      sort(informative_order[seq_len(K_A)])
    }

    sources_data <- lapply(seq_len(K), function(k) {
      source_type <- if (k %in% A_true) "informative" else "noninformative"
      chol_k <- source_bank[[k]][[source_type]]$chol
      .make_data_from_Z(Z_sources[[k]], chol_k, p)
    })
    names(sources_data) <- paste0("source_", seq_len(K))

    # Detect sources:
    # 1. rank sources by K-fold held-out target canonical-correlation gain;
    # 2. select q in {0,...,K} by another K-fold held-out CC score.
    det <- tryCatch(
      detect_trans_scca_sources(
        target = target_data,
        sources = sources_data,
        r = config$r,
        repeats = config$detect_repeats,
        nfolds = config$detect_folds,
        seed = fit_seed + 5000L,
        initialization_method = config$initialization_method,
        rho = config$rho_admm,
        tol = config$admm_tol,
        maxit = config$admm_maxit,
        inner_tol = config$admm_inner_tol,
        inner_maxit = config$admm_inner_maxit
      ),
      error = function(e) {
        list(
          informative = integer(0),
          selected_q = NA_integer_,
          detection_error = conditionMessage(e)
        )
      }
    )

    if (is.null(det$detection_error)) {
      det_metrics <- source_detection_metrics(
        selected = det$informative,
        truth = A_true,
        K = K
      )

      fit_detected <- if (length(det$informative) == 0L) {
        fit_target
      } else {
        .fit_one_split_method(
          method = "Detected",
          target_data = target_data,
          sources = sources_data,
          informative = det$informative,
          config = config,
          fit_seed = fit_seed
        )
      }
    } else {
      det_metrics <- c(
        TPR = NA_real_,
        FPR = NA_real_,
        selected_count = NA_real_
      )
      fit_detected <- structure(
        list(error_message = paste0(
          "Detection failed: ", det$detection_error
        )),
        class = "trans_scca_failed"
      )
    }

    fit_pooled <- .fit_one_split_method(
      method = "Pooled-All",
      target_data = target_data,
      sources = sources_data,
      informative = seq_len(K),
      config = config,
      fit_seed = fit_seed
    )

    rows[[row_id]] <- .result_row(
      fit_target, population, "Target-only",
      p, m, n, rank_type, K_A, rep_id
    )
    row_id <- row_id + 1L

    rows[[row_id]] <- .result_row(
      fit_pooled, population, "Pooled-All",
      p, m, n, rank_type, K_A, rep_id
    )
    row_id <- row_id + 1L

    rows[[row_id]] <- .result_row(
      fit_detected, population, "Detected",
      p, m, n, rank_type, K_A, rep_id,
      detected_metrics = det_metrics
    )
    row_id <- row_id + 1L
  }
  do.call(rbind, rows)
}

run_one_scenario_split <- function(
    task,
    config,
    ncores_rep = 1L) {

  p <- task$p
  m <- task$m
  n <- task$n
  rank_type <- task$rank_type

  population_seed <- config$master_seed +
    1000L * p +
    10L * n +
    ifelse(rank_type == "exact", 1L, 2L)

  population <- make_target_population(
    p = p,
    m = m,
    r = config$r,
    s = config$s,
    rank_type = rank_type,
    rho = config$rho_target,
    tail_decay = config$tail_decay,
    seed = population_seed
  )

  source_bank <- prepare_source_population_bank(
    target_population = population,
    K = config$K,
    delta_informative =
      config$delta_informative,
    noninformative_weight =
      config$noninformative_weight,
    seed = population_seed + 10000L
  )
  
  source_diagnostics <- data.frame(
    p = p,
    m = m,
    n = n,
    rank_type = rank_type,
    source = seq_len(config$K),
    informative_distance = vapply(
      source_bank,
      function(z) z$informative$achieved_distance,
      numeric(1)
    ),
    noninformative_distance = vapply(
      source_bank,
      function(z) z$noninformative$achieved_distance,
      numeric(1)
    ),
    endpoint_distance = vapply(
      source_bank,
      function(z) z$informative$endpoint_distance,
      numeric(1)
    )
  )

  run_rep <- function(rep_id) {
    .run_one_rep_split(
      rep_id = rep_id,
      population = population,
      source_bank = source_bank,
      p = p,
      m = m,
      n = n,
      rank_type = rank_type,
      config = config
    )
  }

  if (ncores_rep <= 1L) {

    pb <- utils::txtProgressBar(
      min = 0,
      max = config$n_rep,
      style = 3
    )

    rep_results <- vector("list", config$n_rep)

    for (rep_id in seq_len(config$n_rep)) {
      rep_results[[rep_id]] <- run_rep(rep_id)
      utils::setTxtProgressBar(pb, rep_id)
    }

    close(pb)

  } else {

    ncores_rep <- min(as.integer(ncores_rep), config$n_rep)

    cl <- parallel::makeCluster(ncores_rep)
    doSNOW::registerDoSNOW(cl)

    main_abs <- normalizePath(
      main_file, winslash = "/", mustWork = TRUE
    )
    utils_abs <- normalizePath(
      utils_file, winslash = "/", mustWork = TRUE
    )

    parallel::clusterExport(
      cl,
      c("main_abs", "utils_abs"),
      envir = environment()
    )

    parallel::clusterEvalQ(
      cl,
      {
        source(main_abs)
        source(utils_abs)
        NULL
      }
    )

    parallel::clusterExport(
      cl,
      c(
        "population",
        "source_bank",
        "p",
        "m",
        "n",
        "rank_type",
        "config",
        ".run_one_rep_split",
        ".fit_one_split_method",
        ".result_row",
        ".make_data_from_Z"
      ),
      envir = environment()
    )

    pb <- utils::txtProgressBar(
      min = 0,
      max = config$n_rep,
      style = 3
    )

    progress <- function(n_done) {
      utils::setTxtProgressBar(pb, n_done)
    }

    opts <- list(progress = progress)

    rep_results <- foreach::foreach(
      rep_id = seq_len(config$n_rep),
      .combine = "c",
      .multicombine = TRUE,
      .options.snow = opts
    ) %dopar% {
      list(
        .run_one_rep_split(
          rep_id = rep_id,
          population = population,
          source_bank = source_bank,
          p = p,
          m = m,
          n = n,
          rank_type = rank_type,
          config = config
        )
      )
    }

    close(pb)
    
    try(
      parallel::stopCluster(cl),
      silent = TRUE
    )
    
    foreach::registerDoSEQ()
    
    cl <- NULL
    
    invisible(gc())
    
    # Reset foreach backend after the PSOCK cluster is stopped.
    foreach::registerDoSEQ()
    
    # Remove the stale cluster object so that RStudio does not
    # recursively inspect already-closed socket connections.
    cl <- NULL
    
    invisible(gc())
  }
  
  scenario_result_df <- do.call(rbind, rep_results)

  scenario_summary <- summarize_simulation_results(
    scenario_result_df
  )

  cat(
    "\nResults for completed setting:",
    sprintf(
      " (p,m,n)=(%d,%d,%d), rank=%s\n",
      p, m, n, rank_type
    )
  )

  print(
    scenario_summary[
      ,
      c(
        "K_A",
        "method",
        "mean_error_U",
        "mean_error_V",
        "success_rate"
      ),
      drop = FALSE
    ],
    row.names = FALSE
  )

  detected_rows <- scenario_result_df[
    scenario_result_df$method == "Detected",
    ,
    drop = FALSE
  ]
  if (nrow(detected_rows)) {
    cat("\nDetection diagnostics (raw replications):\n")
    print(
      detected_rows[
        ,
        c("rep", "K_A", "TPR", "FPR", "selected_count",
          "fit_failed", "error_message"),
        drop = FALSE
      ],
      row.names = FALSE
    )
  }

  list(
    results = scenario_result_df,
    diagnostics = source_diagnostics
  )
}

plot_summary_8panel <- function(
    summary_table,
    rank_type = c("exact", "approximate"),
    output_file = NULL,
    pointsize = 8.5
) {
  rank_type <- match.arg(rank_type)
  
  settings_order <- data.frame(
    p = c(150L, 150L, 300L, 300L),
    m = c(150L, 150L, 300L, 300L),
    n = c(200L, 400L, 200L, 400L)
  )
  
  methods <- c("Detected", "Pooled-All", "Target-only")
  method_labels <- c(
    "Trans-sCCA",
    "Pooled-sCCA",
    "Target-sCCA"
  )
  names(method_labels) <- methods
  
  method_colors <- c(
    "purple",
    "orange",
    "blue"
  )
  names(method_colors) <- methods
  
  line_types <- c(1, 1, 1)
  names(line_types) <- methods
  
  point_types <- c(15, 17, 16)
  names(point_types) <- methods
  
  metric_info <- list(
    U = list(
      mean_col = "mean_error_U",
      sd_col   = "sd_error_U",
      ylab     = expression("Prediction error of " * U^{(0)}),
      colhead  = expression(U^{(0)})
    ),
    V = list(
      mean_col = "mean_error_V",
      sd_col   = "sd_error_V",
      ylab     = expression("Prediction error of " * V^{(0)}),
      colhead  = expression(V^{(0)})
    )
  )
  
  d0 <- summary_table[
    summary_table$rank_type == rank_type &
      summary_table$method %in% methods,
    ,
    drop = FALSE
  ]
  
  if (nrow(d0) == 0L) {
    stop("No rows found for rank_type = ", rank_type)
  }
  
  if (is.null(output_file)) {
    output_file <- sprintf("main_%s_8panel.pdf", rank_type)
  }
  
  grDevices::pdf(
    output_file,
    width = 8.4,
    height = 11.0,
    pointsize = pointsize
  )
  
  layout(
    matrix(
      c(
        1, 2,
        3, 4,
        5, 6,
        7, 8,
        9, 9
      ),
      nrow = 5,
      byrow = TRUE
    ),
    heights = c(
      1, 1, 1, 1,
      0.32
    )
  )
  
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  
  par(
    oma = c(0.6, 0.6, 1.6, 0.3),
    mgp = c(2, 0.5, 0),
    tcl = -0.25
  )
  
  panel_id <- 0L
  
  for (ii in seq_len(nrow(settings_order))) {
    pp <- settings_order$p[ii]
    mm <- settings_order$m[ii]
    nn <- settings_order$n[ii]
    
    for (metric in c("U", "V")) {
      panel_id <- panel_id + 1L
      
      info <- metric_info[[metric]]
      
      d <- d0[
        d0$p == pp &
          d0$m == mm &
          d0$n == nn,
        ,
        drop = FALSE
      ]
      
      d <- d[order(d$K_A), , drop = FALSE]
      
      if (nrow(d) == 0L) {
        par(mar = c(2.2, 3.2, 1.8, 0.8))
        plot.new()
        box()
        text(0.5, 0.5, labels = sprintf("No data\n(p,m,n)=(%d,%d,%d)", pp, mm, nn), cex = 0.9)
        next
      }
      
      y_all <- c()
      
      for (method in methods) {
        dm <- d[d$method == method, , drop = FALSE]
        if (nrow(dm) == 0L) next
        
        y_all <- c(
          y_all,
          dm[[info$mean_col]]
        )
      }
      
      y_all <- y_all[
        is.finite(y_all) &
          y_all > 0
      ]
      
      if (length(y_all) == 0L) {
        y_range <- c(1e-3, 1)
      } else {
        y_range <- range(y_all)
        
        if (y_range[1] == y_range[2]) {
          y_range <- y_range * c(0.5, 2)
        } else {
          y_range <- c(
            y_range[1] / 1.3,
            y_range[2] * 1.3
          )
        }
      }
      
      par(mar = c(3.4, 3.6, if (ii == 1) 2.5 else 1.4, 0.8))
      
      plot(
        NA,
        xlim = range(d$K_A, finite = TRUE),
        ylim = y_range,
        log = "y",
        xaxt = "n",
        yaxt = "n",
        xlab = expression("|" * A[h] * "|"),
        ylab = "",
        main = sprintf("(p,m,n)=(%d,%d,%d)", pp, mm, nn),
        cex.main = 1.5,
        cex.lab = 1.8
      )
      
      ## light-gray panel background
      usr <- par("usr")
      rect(
        usr[1], 10 ^ usr[3],
        usr[2], 10 ^ usr[4],
        col = "#F2F2F2",
        border = NA
      )
      
      ## white grid lines
      abline(
        v = sort(unique(d$K_A)),
        col = "white",
        lwd = 0.8
      )
      
      abline(
        h = axTicks(2),
        col = "white",
        lwd = 0.8
      )
      
      ## axes
      axis(
        side = 1,
        at = sort(unique(d$K_A)),
        labels = sort(unique(d$K_A)),
        cex.axis = 1.3
      )
      
      axis(
        side = 2,
        cex.axis = 1
      )
      
      ## redraw panel border
      box(
        col = "#B0B0B0",
        lwd = 0.8
      )
      
      mtext(
        info$ylab,
        side = 2,
        line = 2.2,
        cex = 1.1
      )
      
      #if (ii == 1) {
      #  mtext(info$colhead, side = 3, line = 1.0, cex = 1.0, font = 2)
      #}
      
      for (method in methods) {
        dm <- d[d$method == method, , drop = FALSE]
        dm <- dm[order(dm$K_A), , drop = FALSE]
        if (nrow(dm) == 0L) next
        
        xx <- dm$K_A
        yy <- dm[[info$mean_col]]
        
        lines(
          xx, yy,
          type = "o",
          col = method_colors[method],
          lty = line_types[method],
          pch = point_types[method],
          lwd = 1.2,
          cex = 1.8
        )
      }
    }
  }
  
  
  ## shared legend
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = unname(method_labels[methods]),
    col = unname(method_colors[methods]),
    lty = unname(line_types[methods]),
    pch = unname(point_types[methods]),
    lwd = 1.2,
    horiz = TRUE,
    bty = "n",
    cex = 2,
    xpd = NA
  )
  
  ### outer title
  #mtext(
  #  if (rank_type == "exact") "Exactly low-rank" else "Approximately low-rank",
  # side = 3,
  # outer = TRUE,
  # line = 0.2,
  #  cex = 1.2,
  #  font = 2
  #)
  
  invisible(output_file)
}



run_split_simulation <- function(
    config = simulation_config,
    ncores_rep = 1L) {
  
  tasks <- do.call(
    rbind,
    lapply(config$rank_types, function(rank_type) {
      cbind(
        config$settings,
        rank_type = rank_type,
        stringsAsFactors = FALSE
      )
    })
  )
  
  task_list <- lapply(seq_len(nrow(tasks)), function(i) {
    list(
      p = as.integer(tasks$p[i]),
      m = as.integer(tasks$m[i]),
      n = as.integer(tasks$n[i]),
      rank_type = as.character(tasks$rank_type[i])
    )
  })
  
  # -------------------------------------------------------------
  # Run the 8 scenarios sequentially.
  # Within each scenario, Monte Carlo replications are parallelized.
  # -------------------------------------------------------------
  scenario_results <- vector(
    "list",
    length(task_list)
  )
  
  for (ii in seq_along(task_list)) {
    
    task <- task_list[[ii]]
    
    cat(
      "\n====================================================\n",
      sprintf(
        "Scenario %d/%d: p=%d, m=%d, n=%d, rank=%s\n",
        ii,
        length(task_list),
        task$p,
        task$m,
        task$n,
        task$rank_type
      ),
      sprintf(
        "Parallel workers for replications: %d\n",
        ncores_rep
      ),
      "====================================================\n",
      sep = ""
    )
    
    scenario_start <- Sys.time()
    
    scenario_results[[ii]] <- run_one_scenario_split(
      task = task,
      config = config,
      ncores_rep = ncores_rep
    )
    
    # -------------------------------------------------------------
    # Immediately summarize, save, and plot THIS completed setting.
    # -------------------------------------------------------------
    
    scenario_raw <- scenario_results[[ii]]$results
    
    scenario_summary <- summarize_simulation_results(
      scenario_raw
    )
    
    scenario_diagnostics <-
      scenario_results[[ii]]$diagnostics
    
    
    # Make output directory if needed
    if (!dir.exists(config$output_dir)) {
      dir.create(
        config$output_dir,
        recursive = TRUE
      )
    }
    
    
    # Unique tag for this setting
    setting_tag <- sprintf(
      "%s_p%d_m%d_n%d",
      task$rank_type,
      task$p,
      task$m,
      task$n
    )
    
    
    # -------------------------------------------------------------
    # Save this setting immediately
    # -------------------------------------------------------------
    
    saveRDS(
      scenario_raw,
      file.path(
        config$output_dir,
        paste0(
          "split_raw_",
          setting_tag,
          ".rds"
        )
      )
    )
    
    saveRDS(
      scenario_summary,
      file.path(
        config$output_dir,
        paste0(
          "split_summary_",
          setting_tag,
          ".rds"
        )
      )
    )
    
    utils::write.csv(
      scenario_raw,
      file.path(
        config$output_dir,
        paste0(
          "split_raw_",
          setting_tag,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    utils::write.csv(
      scenario_summary,
      file.path(
        config$output_dir,
        paste0(
          "split_summary_",
          setting_tag,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    utils::write.csv(
      scenario_diagnostics,
      file.path(
        config$output_dir,
        paste0(
          "split_source_distance_",
          setting_tag,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    
    # -------------------------------------------------------------
    # Draw PDFs immediately
    # -------------------------------------------------------------
    
    plot_simulation_summary(
      scenario_summary,
      output_dir = config$output_dir,
      prefix = "split"
    )
    
    cat(
      sprintf(
        "\nPDF and result files saved for setting: %s\n",
        setting_tag
      )
    )
    
    scenario_end <- Sys.time()
    
    elapsed_minutes <- as.numeric(
      difftime(
        scenario_end,
        scenario_start,
        units = "mins"
      )
    )
    
    cat(
      sprintf(
        "\nFinished scenario %d/%d in %.2f minutes.\n",
        ii,
        length(task_list),
        elapsed_minutes
      )
    )
  }
  
  
  # -------------------------------------------------------------
  # Combine results from all scenarios.
  # -------------------------------------------------------------
  results <- do.call(
    rbind,
    lapply(
      scenario_results,
      `[[`,
      "results"
    )
  )
  
  diagnostics <- do.call(
    rbind,
    lapply(
      scenario_results,
      `[[`,
      "diagnostics"
    )
  )
  
  summary_table <- summarize_simulation_results(
    results
  )
  
  
  # -------------------------------------------------------------
  # Save results.
  # -------------------------------------------------------------
  if (!dir.exists(config$output_dir)) {
    dir.create(
      config$output_dir,
      recursive = TRUE
    )
  }
  
  saveRDS(
    results,
    file.path(
      config$output_dir,
      "split_raw_results.rds"
    )
  )
  
  saveRDS(
    summary_table,
    file.path(
      config$output_dir,
      "split_summary.rds"
    )
  )
  
  utils::write.csv(
    results,
    file.path(
      config$output_dir,
      "split_raw_results.csv"
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    summary_table,
    file.path(
      config$output_dir,
      "split_summary.csv"
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    diagnostics,
    file.path(
      config$output_dir,
      "source_distance_diagnostics.csv"
    ),
    row.names = FALSE
  )
  
  
  # -------------------------------------------------------------
  # Plot the final 8-panel figures.
  # -------------------------------------------------------------
  
  if ("exact" %in% summary_table$rank_type) {
    plot_summary_8panel(
      summary_table = summary_table,
      rank_type = "exact",
      output_file = file.path(
        config$output_dir,
        "main_exact_8panel.pdf"
      )
    )
  }
  
  if ("approximate" %in% summary_table$rank_type) {
    plot_summary_8panel(
      summary_table = summary_table,
      rank_type = "approximate",
      output_file = file.path(
        config$output_dir,
        "main_approximate_8panel.pdf"
      )
    )
  }
  
  
  invisible(
    list(
      raw = results,
      summary = summary_table,
      diagnostics = diagnostics
    )
  )
}

# -------------------------------------------------------------------------
# RUN
# -------------------------------------------------------------------------
# For a quick test:

#simulation_config$settings <- data.frame(
#  p = 200L,
#  m = 200L,
#  n = 200L
#)

#simulation_config$rank_types <- "exact"
simulation_config$initialization_method <- "ADMM"
simulation_config$tune_initialization <- FALSE
simulation_config$n_rep <- 60L

start_time <- Sys.time()

ans <- run_split_simulation(
  simulation_config,
  ncores_rep = 20L
)

end_time <- Sys.time()

elapsed_minutes <- as.numeric(
  difftime(
    end_time,
    start_time,
    units = "mins"
  )
)

cat(
  "\nADMM total time =",
  as.numeric(
    difftime(end_time, start_time, units = "mins")
  ),
  "minutes\n"
)

#cat(
#  "\n========================================\n",
#  "All simulations finished.\n",
#  "Total elapsed time: ",
#  round(elapsed_minutes, 2),
#  " minutes\n",
#  "Total elapsed time: ",
#  round(elapsed_minutes / 60, 2),
# " hours\n",
#  "========================================\n",
#  sep = ""
#)
#
# For a larger Windows run, for example:
# simulation_config$n_rep <- 100L
# ans <- run_split_simulation(simulation_config, ncores = 8L)
#
# The x-axis in the generated figures is K_A = |A_h|, not total K, because
# total K is fixed at 10 in this design.
