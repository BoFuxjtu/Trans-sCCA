rm(list = ls())
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
source("Trans-sparse_CCA-cor_CV_fullsample.R")
if (!requireNamespace("foreach", quietly = TRUE) ||
    !requireNamespace("doSNOW", quietly = TRUE)) {
  stop(
    "Install packages 'foreach' and 'doSNOW' first."
  )
}

library(foreach)
library(doSNOW)
## ============================================================
## Load processed TCGA data
## ============================================================

tcga_processed <- readRDS(
  #"E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed_ESCA.rds"
  #"E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed_LUAD.rds"
  #"E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed_LUSC.rds"
  "E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed_HNSC.rds"
  #"E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed_BLCA.rds"
)


domains_screened <-
  tcga_processed$domains

target_name <-
  tcga_processed$target_name

source_names <-
  tcga_processed$source_names

selected_genes <-
  tcga_processed$selected_genes


print(
  tcga_processed$preprocessing
)

## ============================================================
## Construct target and candidate sources
## ============================================================

target <- list(
  X = domains_screened[[target_name]]$X,
  Y = domains_screened[[target_name]]$Y
)

sources <- lapply(
  source_names,
  function(cc) {
    
    list(
      X = domains_screened[[cc]]$X,
      Y = domains_screened[[cc]]$Y
    )
  }
)

names(sources) <- source_names


## ============================================================
## Final checks
## ============================================================

print(
  sapply(
    domains_screened,
    function(z) {
      c(
        n = nrow(z$X),
        p = ncol(z$X),
        m = ncol(z$Y),
        NA_X = sum(is.na(z$X)),
        NA_Y = sum(is.na(z$Y))
      )
    }
  )
)


## ============================================================
## Real-data configuration
## ============================================================

real_config <- list(
  r = 2L,
  
  n_rep = 50L,
  #train_n = 100L, # 100 for ESCA
  #train_n = 288L, # 288 for LUAD
  #train_n = 260L, # 260 for LUSC
  train_n = 170L, # 170 for HNSC
  #train_n = 272L, # 272 for HNSC
  master_seed = 20260922L,
  ncores_rep = 20L,
  
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
  
  detect_repeats = 1L,
  detect_folds = 5L,
  
  full_sample_fit = TRUE,
  
  center = TRUE
)

n_test <- nrow(target$X) -
  real_config$train_n

if (
  real_config$r >= n_test
) {
  stop(
    "r must be smaller than the held-out target test sample size."
  )
}

## ============================================================
## Fit one of the three estimators
## ============================================================

.fit_real_method <- function(
    method,
    target_train,
    sources,
    informative,
    config,
    fit_seed) {
  
  aggregation_mode <- switch(
    method,
    "Detected" = "adaptive",
    "Pooled-All" = "transfer",
    "Target-only" = "target",
    stop("Unknown method.")
  )
  
  normalization_mode <- if (
    method == "Pooled-All"
  ) {
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
  
  tryCatch(
    trans_scca(
      target = target_train,
      sources = sources_use,
      informative = informative_use,
      r = config$r,
      s_u = NULL,
      s_v = NULL,
      full_sample = config$full_sample_fit,
      center = config$center,
      seed = fit_seed,
      initialization_method =
        config$initialization_method,
      rho_admm = config$rho_admm,
      admm_tol = config$admm_tol,
      admm_maxit = config$admm_maxit,
      admm_inner_tol =
        config$admm_inner_tol,
      admm_inner_maxit =
        config$admm_inner_maxit,
      tune_initialization =
        config$tune_initialization,
      tune_penalties =
        config$tune_penalties,
      cv_grid = config$cv_grid,
      cv_folds = config$cv_folds,
      normalization_ridge =
        config$normalization_ridge,
      aggregation_mode =
        aggregation_mode,
      normalization_mode =
        normalization_mode
    ),
    error = function(e) {
      structure(
        list(
          error_message =
            conditionMessage(e)
        ),
        class = "trans_scca_failed"
      )
    }
  )
}

.fit_failed <- function(fit) {
  inherits(
    fit,
    "trans_scca_failed"
  )
}

## ============================================================
## Held-out target canonical correlation
## ============================================================

.evaluate_test_cc2 <- function(
    fit,
    target_test,
    r) {
  
  metric_names <- c(
    paste0("CC", seq_len(r)),
    "Total_CC",
    "Mean_CC"
  )
  
  if (.fit_failed(fit)) {
    out <- rep(
      NA_real_,
      length(metric_names)
    )
    names(out) <- metric_names
    return(out)
  }
  
  X_score <- target_test$X %*% fit$U
  Y_score <- target_test$Y %*% fit$V
  
  cc <- tryCatch(
    stats::cancor(
      X_score,
      Y_score
    ),
    error = function(e) NULL
  )
  
  if (is.null(cc)) {
    out <- rep(
      NA_real_,
      length(metric_names)
    )
    names(out) <- metric_names
    return(out)
  }
  
  rho <- rep(
    NA_real_,
    r
  )
  
  nr <- min(
    r,
    length(cc$cor)
  )
  
  if (nr > 0L) {
    rho[seq_len(nr)] <-
      cc$cor[seq_len(nr)]
  }
  
  out <- c(
    rho,
    Total_CC = if (
      all(is.finite(rho))
    ) {
      sum(rho)
    } else {
      NA_real_
    },
    Mean_CC = if (
      all(is.finite(rho))
    ) {
      mean(rho)
    } else {
      NA_real_
    }
  )
  
  names(out) <- metric_names
  
  out
}

.evaluate_test_cc <- function(
    fit,
    target_test,
    r,
    tol = sqrt(.Machine$double.eps)) {
  
  metric_names <- c(
    paste0("CC", seq_len(r)),
    "Total_CC",
    "Mean_CC"
  )
  
  if (.fit_failed(fit)) {
    out <- rep(
      NA_real_,
      length(metric_names)
    )
    names(out) <- metric_names
    return(out)
  }
  
  
  ## --------------------------------------------------------
  ## Estimated test canonical score matrices
  ## --------------------------------------------------------
  
  X_score <-
    target_test$X %*% fit$U
  
  Y_score <-
    target_test$Y %*% fit$V
  
  
  ## Center on the held-out test set, as cancor() does
  X_score <- sweep(
    X_score,
    2,
    colMeans(X_score),
    "-"
  )
  
  Y_score <- sweep(
    Y_score,
    2,
    colMeans(Y_score),
    "-"
  )
  
  
  ## --------------------------------------------------------
  ## SVD-based orthonormal bases
  ## --------------------------------------------------------
  
  sx <- svd(
    X_score,
    nu = min(
      nrow(X_score),
      ncol(X_score)
    ),
    nv = 0
  )
  
  sy <- svd(
    Y_score,
    nu = min(
      nrow(Y_score),
      ncol(Y_score)
    ),
    nv = 0
  )
  
  
  tol_x <- tol *
    max(dim(X_score)) *
    max(sx$d, 0)
  
  tol_y <- tol *
    max(dim(Y_score)) *
    max(sy$d, 0)
  
  
  rank_x <- sum(
    sx$d > tol_x
  )
  
  rank_y <- sum(
    sy$d > tol_y
  )
  
  
  ## --------------------------------------------------------
  ## Canonical correlations between score subspaces
  ## --------------------------------------------------------
  
  rho <- numeric(r)
  
  if (
    rank_x > 0L &&
    rank_y > 0L
  ) {
    
    Qx <- sx$u[
      ,
      seq_len(rank_x),
      drop = FALSE
    ]
    
    Qy <- sy$u[
      ,
      seq_len(rank_y),
      drop = FALSE
    ]
    
    rho_nonzero <- svd(
      crossprod(
        Qx,
        Qy
      ),
      nu = 0,
      nv = 0
    )$d
    
    rho_nonzero <- pmin(
      1,
      pmax(
        0,
        rho_nonzero
      )
    )
    
    nr <- min(
      r,
      length(rho_nonzero)
    )
    
    rho[
      seq_len(nr)
    ] <- rho_nonzero[
      seq_len(nr)
    ]
  }
  
  
  out <- c(
    rho,
    Total_CC = sum(rho),
    Mean_CC = mean(rho)
  )
  
  names(out) <- metric_names
  
  out
}

.make_real_result_row <- function(
    rep_id,
    method,
    cc_result,
    fit) {
  
  failed <- .fit_failed(fit)
  
  metric_df <- as.data.frame(
    as.list(cc_result),
    check.names = FALSE
  )
  
  cbind(
    data.frame(
      rep = rep_id,
      method = method,
      fit_failed = failed,
      error_message =
        if (failed) {
          fit$error_message
        } else {
          ""
        },
      stringsAsFactors = FALSE
    ),
    metric_df
  )
}
## ============================================================
## Repeated 100/25 target train-test splits
## ============================================================

.run_one_real_rep <- function(
    rep_id,
    target,
    sources,
    source_names,
    config) {
  
  n0 <- nrow(target$X)
  
  ## ----------------------------------------------------------
  ## 1. Random 100/25 split
  ## ----------------------------------------------------------
  
  split_seed <-
    config$master_seed +
    10000L * rep_id
  
  set.seed(split_seed)
  
  train_id <- sort(
    sample(
      seq_len(n0),
      size = config$train_n,
      replace = FALSE
    )
  )
  
  test_id <- setdiff(
    seq_len(n0),
    train_id
  )
  
  target_train <- list(
    X = target$X[
      train_id,
      ,
      drop = FALSE
    ],
    Y = target$Y[
      train_id,
      ,
      drop = FALSE
    ]
  )
  
  target_test <- list(
    X = target$X[
      test_id,
      ,
      drop = FALSE
    ],
    Y = target$Y[
      test_id,
      ,
      drop = FALSE
    ]
  )
  
  split_info <- list(
    train_id = train_id,
    test_id = test_id,
    train_patient =
      rownames(target$X)[train_id],
    test_patient =
      rownames(target$X)[test_id]
  )
  
  
  ## ----------------------------------------------------------
  ## 2. Target-sCCA
  ## ----------------------------------------------------------
  
  fit_seed <-
    config$master_seed +
    100000L * rep_id
  
  fit_target <- .fit_real_method(
    method = "Target-only",
    target_train = target_train,
    sources = list(),
    informative = integer(0),
    config = config,
    fit_seed = fit_seed
  )
  
  
  ## ----------------------------------------------------------
  ## 3. Source detection
  ## ----------------------------------------------------------
  
  det <- tryCatch(
    detect_trans_scca_sources(
      target = target_train,
      sources = sources,
      r = config$r,
      repeats =
        config$detect_repeats,
      nfolds =
        config$detect_folds,
      center =
        config$center,
      seed =
        fit_seed + 5000L,
      initialization_method =
        config$initialization_method,
      rho =
        config$rho_admm,
      tol =
        config$admm_tol,
      maxit =
        config$admm_maxit,
      inner_tol =
        config$admm_inner_tol,
      inner_maxit =
        config$admm_inner_maxit
    ),
    error = function(e) {
      list(
        informative = integer(0),
        selected_q = NA_integer_,
        detection_error =
          conditionMessage(e)
      )
    }
  )
  
  
  ## ----------------------------------------------------------
  ## 4. Detected-Trans-sCCA
  ## ----------------------------------------------------------
  
  if (is.null(det$detection_error)) {
    
    if (length(det$informative) == 0L) {
      
      fit_detected <- fit_target
      
    } else {
      
      fit_detected <-
        .fit_real_method(
          method = "Detected",
          target_train =
            target_train,
          sources = sources,
          informative =
            det$informative,
          config = config,
          fit_seed = fit_seed
        )
    }
    
  } else {
    
    fit_detected <- structure(
      list(
        error_message =
          paste0(
            "Detection failed: ",
            det$detection_error
          )
      ),
      class = "trans_scca_failed"
    )
  }
  
  
  ## ----------------------------------------------------------
  ## 5. Pooled-sCCA
  ## ----------------------------------------------------------
  
  fit_pooled <- .fit_real_method(
    method = "Pooled-All",
    target_train = target_train,
    sources = sources,
    informative =
      seq_along(sources),
    config = config,
    fit_seed = fit_seed
  )
  
  
  ## ----------------------------------------------------------
  ## 6. Held-out canonical correlations
  ## ----------------------------------------------------------
  
  cc_detected <-
    .evaluate_test_cc(
      fit_detected,
      target_test,
      config$r
    )
  
  cc_pooled <-
    .evaluate_test_cc(
      fit_pooled,
      target_test,
      config$r
    )
  
  cc_target <-
    .evaluate_test_cc(
      fit_target,
      target_test,
      config$r
    )
  
  
  ## ----------------------------------------------------------
  ## 7. Detection information
  ## ----------------------------------------------------------
  
  selected_source_names <-
    if (length(det$informative) > 0L) {
      paste(
        source_names[
          det$informative
        ],
        collapse = ","
      )
    } else {
      ""
    }
  
  detection_row <- data.frame(
    rep = rep_id,
    selected_q =
      det$selected_q,
    selected_sources =
      selected_source_names,
    detection_failed =
      !is.null(
        det$detection_error
      ),
    detection_error =
      if (
        is.null(
          det$detection_error
        )
      ) {
        ""
      } else {
        det$detection_error
      },
    stringsAsFactors = FALSE
  )
  
  
  ## ----------------------------------------------------------
  ## 8. Three result rows
  ## ----------------------------------------------------------
  
  result_rows <- list(
    
    .make_real_result_row(
      rep_id = rep_id,
      method =
        "Detected-Trans-sCCA",
      cc_result = cc_detected,
      fit = fit_detected
    ),
    
    .make_real_result_row(
      rep_id = rep_id,
      method =
        "Pooled-sCCA",
      cc_result = cc_pooled,
      fit = fit_pooled
    ),
    
    .make_real_result_row(
      rep_id = rep_id,
      method =
        "Target-sCCA",
      cc_result = cc_target,
      fit = fit_target
    )
  )
  
  
  list(
    results = do.call(
      rbind,
      result_rows
    ),
    detection = detection_row,
    split = split_info
  )
}



n0 <- nrow(target$X)

#stopifnot(
#  n0 == 125L,
#  real_config$train_n < n0
#)

ncores_rep <- min(
  as.integer(real_config$ncores_rep),
  real_config$n_rep
)

cl <- parallel::makeCluster(
  ncores_rep
)

doSNOW::registerDoSNOW(cl)


## ------------------------------------------------------------
## Load the main Trans-sCCA functions on every worker
## ------------------------------------------------------------

main_file <- normalizePath(
  "Trans-sparse_CCA-cor_CV_fullsample.R",
  winslash = "/",
  mustWork = TRUE
)

parallel::clusterExport(
  cl,
  "main_file",
  envir = environment()
)

parallel::clusterEvalQ(
  cl,
  {
    source(main_file)
    NULL
  }
)


## ------------------------------------------------------------
## Progress bar
## ------------------------------------------------------------

pb <- utils::txtProgressBar(
  min = 0,
  max = real_config$n_rep,
  style = 3
)

progress <- function(n_done) {
  utils::setTxtProgressBar(
    pb,
    n_done
  )
}

opts <- list(
  progress = progress
)


## ------------------------------------------------------------
## Parallel repetitions
## ------------------------------------------------------------

rep_results <- foreach::foreach(
  rep_id =
    seq_len(real_config$n_rep),

  .combine = "c",
  .multicombine = TRUE,

  .export = c(
    ".run_one_real_rep",
    ".fit_real_method",
    ".fit_failed",
    ".evaluate_test_cc",
    ".make_real_result_row",
    "target",
    "sources",
    "source_names",
    "real_config"
  ),

  .options.snow = opts
) %dopar% {

  list(
    .run_one_real_rep(
      rep_id = rep_id,
      target = target,
      sources = sources,
      source_names = source_names,
      config = real_config
    )
  )
}


close(pb)

parallel::stopCluster(cl)

foreach::registerDoSEQ()

rm(cl)

invisible(gc())


## ============================================================
## Collect parallel results
## ============================================================

real_results <- do.call(
  rbind,
  lapply(
    rep_results,
    `[[`,
    "results"
  )
)

detection_results <- do.call(
  rbind,
  lapply(
    rep_results,
    `[[`,
    "detection"
  )
)

split_list <- lapply(
  rep_results,
  `[[`,
  "split"
)

.summarize_real_results <- function(
    results,
    r) {
  
  methods <- unique(
    results$method
  )
  
  metrics <- c(
    paste0(
      "CC",
      seq_len(r)
    ),
    "Total_CC",
    "Mean_CC"
  )
  
  rows <- list()
  row_id <- 1L
  
  for (mm in methods) {
    
    d <- results[
      results$method == mm,
      ,
      drop = FALSE
    ]
    
    for (metric in metrics) {
      
      x <- d[[metric]]
      
      rows[[row_id]] <-
        data.frame(
          method = mm,
          metric = metric,
          mean = mean(
            x,
            na.rm = TRUE
          ),
          sd = stats::sd(
            x,
            na.rm = TRUE
          ),
          n_success =
            sum(
              is.finite(x)
            ),
          stringsAsFactors = FALSE
        )
      
      row_id <- row_id + 1L
    }
  }
  
  do.call(
    rbind,
    rows
  )
}

real_summary <-
  .summarize_real_results(
    results = real_results,
    r = real_config$r
  )

print(
  real_summary,
  row.names = FALSE
)

print(
  real_summary,
  row.names = FALSE
)

## ============================================================
## Save real-data results
## ============================================================

output_file <- sprintf(
  "TCGA_realdata_r%d_rep%d_results_HNSC.rds",
  real_config$r,
  real_config$n_rep
)

saveRDS(
  list(
    results = real_results,
    summary = real_summary,
    detection = detection_results,
    splits = split_list,
    config = real_config,
    selected_genes = selected_genes
  ),
  output_file
)

write.csv(
  real_results,
  sprintf(
    "TCGA_realdata_r%d_rep%d_raw_HNSC.csv",
    real_config$r,
    real_config$n_rep
  ),
  row.names = FALSE
)

write.csv(
  real_summary,
  sprintf(
    "TCGA_realdata_r%d_rep%d_summary_HNSC.csv",
    real_config$r,
    real_config$n_rep
  ),
  row.names = FALSE
)

write.csv(
  detection_results,
  sprintf(
    "TCGA_realdata_r%d_rep%d_detection_HNSC.csv",
    real_config$r,
    real_config$n_rep
  ),
  row.names = FALSE
)

cat(
  "\nResults saved to:",
  normalizePath(output_file),
  "\n"
)


plot_tcga_cc <- function(
    output_file = NULL,
    width = 9,
    height = 5.5) {
  
  ## ============================================================
  ## Data
  ## ============================================================
  
  targets <- c(
    "ESCA",
    "LUAD",
    "LUSC",
    "HNSC",
    "BLCA"
  )
  
  methods <- c(
    "Detected-Trans-sCCA",
    "Pooled-sCCA",
    "Target-sCCA"
  )
  
  ## rows = target datasets
  ## columns = methods
  
  CC1 <- rbind(
    ESCA = c(0.8563509, 0.7904811, 0.7154606),
    LUAD = c(0.7422207, 0.7142059, 0.6864240),
    LUSC = c(0.5928740, 0.6182290, 0.5739226),
    HNSC = c(0.6011428, 0.5017765, 0.6382757),
    BLCA = c(0.8514379, 0.8784440, 0.7822140)
  )
  
  CC2 <- rbind(
    ESCA = c(0.2463730, 0.1841894, 0.2282427),
    LUAD = c(0.5087347, 0.5146817, 0.5636632),
    LUSC = c(0.4307235, 0.4416419, 0.3168082),
    HNSC = c(0.3539746, 0.1269888, 0.3544192),
    BLCA = c(0.2646080, 0.0682418, 0.3690979)
  )
  
  colnames(CC1) <- methods
  colnames(CC2) <- methods
  
  
  ## ============================================================
  ## PDF device, if requested
  ## ============================================================
  
  if (!is.null(output_file)) {
    grDevices::pdf(
      output_file,
      width = width,
      height = height
    )
  }
  
  
  ## ============================================================
  ## Plot settings
  ## ============================================================
  
  method_col <- c(
    "purple",
    "orange",
    "blue"
  )
  
  n_target <- length(targets)
  n_method <- length(methods)
  
  group_gap <- 1.25
  bar_width <- 0.68
  inner_gap <- 0.14
  
  group_centers <- seq(
    from = 1,
    by = group_gap + n_method * (bar_width + inner_gap),
    length.out = n_target
  )
  
  offsets <- (
    seq_len(n_method) -
      mean(seq_len(n_method))
  ) * (bar_width + inner_gap)
  
  x_pos <- outer(
    group_centers,
    offsets,
    "+"
  )
  
  ymax <- 1.35
  
  
  ## ============================================================
  ## Empty plotting region
  ## ============================================================
  
  ## ============================================================
  ## Layout: upper panel for plot, lower panel for legends
  ## ============================================================
  
  layout(
    matrix(c(1, 2), nrow = 2),
    heights = c(5, 1.15)
  )
  
  old_par <- par(
    mar = c(4.8, 5.2, 1.0, 1.0)
  )
  
  on.exit({
    par(old_par)
    layout(1)
    
    if (!is.null(output_file)) {
      dev.off()
    }
  })
  
  ## ============================================================
  ## Main plot
  ## ============================================================
  
  plot(
    NA,
    xlim = c(
      min(x_pos) - 0.7,
      max(x_pos) + 0.7
    ),
    ylim = c(0, ymax),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    bty = "l"
  )
  
  axis(
    1,
    at = group_centers,
    labels = targets,
    cex.axis = 1
  )
  
  axis(
    2,
    las = 1,
    cex.axis = 0.9
  )
  
  mtext(
    "Target dataset",
    side = 1,
    line = 3.2,
    cex = 1.4
  )
  
  mtext(
    "Out-of-sample canonical correlation",
    side = 2,
    line = 3.5,
    cex = 1
  )
  
  
  ## ============================================================
  ## Bars
  ## CC1 = solid lower part
  ## CC2 = striped upper part
  ## ============================================================
  
  for (i in seq_len(n_target)) {
    
    for (j in seq_len(n_method)) {
      
      xleft <- x_pos[i, j] - bar_width / 2
      xright <- x_pos[i, j] + bar_width / 2
      
      ## CC1: solid
      rect(
        xleft = xleft,
        ybottom = 0,
        xright = xright,
        ytop = CC1[i, j],
        col = method_col[j],
        border = "black",
        lwd = 0.7
      )
      
      ## CC2: same background colour
      rect(
        xleft = xleft,
        ybottom = CC1[i, j],
        xright = xright,
        ytop = CC1[i, j] + CC2[i, j],
        col = method_col[j],
        border = "black",
        lwd = 0.7
      )
      
      ## diagonal stripes over CC2
      rect(
        xleft = xleft,
        ybottom = CC1[i, j],
        xright = xright,
        ytop = CC1[i, j] + CC2[i, j],
        col = "black",
        border = NA,
        density = 16,
        angle = 45
      )
    }
  }
  
  
  ## ============================================================
  ## Legends
  ## ============================================================
  
  ## ============================================================
  ## Legend panel
  ## ============================================================
  
  par(
    mar = c(0, 0, 0, 0)
  )
  
  plot.new()
  
  plot.window(
    xlim = c(0, 1),
    ylim = c(0, 1)
  )
  
  legend(
    x = 0.5,
    y = 0.72,
    legend = c(
      "Trans-sCCA",
      "Pooled-sCCA",
      "Target-sCCA"
    ),
    fill = method_col,
    border = "black",
    horiz = TRUE,
    bty = "n",
    cex = 1.4,
    xjust = 0.5,
    yjust = 0.5
  )
  
  legend(
    x = 0.5,
    y = 0.25,
    legend = c("CC1", "CC2"),
    fill = c("grey75", "black"),
    density = c(NA, 16),
    angle = c(NA, 45),
    border = "black",
    horiz = TRUE,
    bty = "n",
    cex = 1.4,
    xjust = 0.5,
    yjust = 0.5
  )
  
  invisible(
    list(
      CC1 = CC1,
      CC2 = CC2,
      Total_CC = CC1 + CC2
    )
  )
}

plot_tcga_cc()

#X_score <- X_test %*% U_hat
#Y_score <- Y_test %*% V_hat

#cc <- stats::cancor(X_score, Y_score)

#cc$cor          # canonical correlation
#sum(cc$cor)     # total canonical correlation