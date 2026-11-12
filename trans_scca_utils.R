# Shared utilities for the Trans-sparse CCA simulation study.
#
# This file generates the target/source population covariance matrices,
# simulates Gaussian samples, computes prediction errors, and summarizes
# Monte Carlo results.  It is shared by the split-sample and full-refit
# simulation scripts.

.sim_sym <- function(A) {
  (A + t(A)) / 2
}

.sim_inv_root <- function(A, tol = 1e-10) {
  ee <- eigen(.sim_sym(A), symmetric = TRUE)
  scale <- max(1, max(abs(ee$values)))
  if (min(ee$values) <= tol * scale) {
    stop("Matrix is not sufficiently positive definite.")
  }
  tcrossprod(
    sweep(ee$vectors, 2, 1 / sqrt(ee$values), `*`),
    ee$vectors
  )
}

.sim_root <- function(A, tol = 1e-10) {
  ee <- eigen(.sim_sym(A), symmetric = TRUE)
  scale <- max(1, max(abs(ee$values)))
  if (min(ee$values) < -tol * scale) {
    stop("Matrix is not positive semidefinite.")
  }
  d <- sqrt(pmax(ee$values, 0))
  tcrossprod(
    sweep(ee$vectors, 2, d, `*`),
    ee$vectors
  )
}

.sim_spectral_norm_sym <- function(A) {
  A <- .sim_sym(A)

  if (requireNamespace("RSpectra", quietly = TRUE) && nrow(A) >= 100L) {
    value <- RSpectra::eigs_sym(
      A,
      k = 1L,
      which = "LM"
    )$values[1L]
    return(abs(value))
  }

  max(abs(eigen(A, symmetric = TRUE, only.values = TRUE)$values))
}

toeplitz_covariance <- function(d, rho = 0.3) {
  index <- seq_len(d)
  rho ^ abs(outer(index, index, `-`))
}

.default_leading_lambdas <- function(r) {
  if (r == 1L) return(0.9)
  if (r == 2L) return(c(0.9, 0.8))
  if (r == 3L) return(c(0.9, 0.8, 0.7))

  # For r > 3, continue a mild linear decay while keeping all values below 1.
  values <- 0.9 - 0.1 * (seq_len(r) - 1L)
  if (min(values) <= 0) {
    values <- seq(0.9, 0.4, length.out = r)
  }
  values
}

.make_sparse_sigma_orthonormal <- function(
    Sigma,
    support,
    r,
    value_set = c(-2, -1, 0, 1, 2),
    max_attempts = 1000L) {

  d <- nrow(Sigma)
  s <- length(support)

  if (s < r) {
    stop("The support size must be at least r.")
  }

  for (attempt in seq_len(max_attempts)) {
    raw <- matrix(0, d, r)
    active <- matrix(
      sample(value_set, s * r, replace = TRUE),
      nrow = s,
      ncol = r
    )

    # Avoid an active row being identically zero.
    zero_rows <- which(rowSums(active^2) == 0)
    if (length(zero_rows)) next

    raw[support, ] <- active

    gram <- crossprod(raw, Sigma %*% raw)
    ee <- eigen(.sim_sym(gram), symmetric = TRUE, only.values = TRUE)$values

    if (min(ee) > 1e-8) {
      return(raw %*% .sim_inv_root(gram))
    }
  }

  stop("Could not generate a full-rank sparse Sigma-orthonormal loading matrix.")
}

.complete_sigma_orthonormal_basis <- function(
    Sigma,
    leading,
    seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  d <- nrow(Sigma)
  r <- ncol(leading)

  Sigma_half <- .sim_root(Sigma)
  Sigma_inv_half <- .sim_inv_root(Sigma)

  Q0 <- Sigma_half %*% leading

  # Numerical cleanup: force the leading block to Euclidean orthonormality
  # without changing its span.
  Q0 <- qr.Q(qr(Q0))[, seq_len(r), drop = FALSE]

  if (r == d) {
    return(Sigma_inv_half %*% Q0)
  }

  for (attempt in seq_len(50L)) {
    Z <- matrix(rnorm(d * (d - r)), d, d - r)
    Z <- Z - Q0 %*% crossprod(Q0, Z)

    Qperp <- qr.Q(qr(Z))[, seq_len(d - r), drop = FALSE]

    if (ncol(Qperp) == d - r) {
      Qfull <- cbind(Q0, Qperp)
      return(Sigma_inv_half %*% Qfull)
    }
  }

  stop("Could not complete the Sigma-orthonormal basis.")
}

.make_joint_covariance <- function(
    Sigma_x,
    Sigma_y,
    U,
    V,
    lambda) {

  Sigma_xy <- Sigma_x %*% U %*% diag(lambda, length(lambda)) %*%
    t(V) %*% Sigma_y

  rbind(
    cbind(Sigma_x, Sigma_xy),
    cbind(t(Sigma_xy), Sigma_y)
  )
}

make_target_population <- function(
    p,
    m = p,
    r = 2L,
    s = 10L,
    rank_type = c("exact", "approximate"),
    rho = 0.3,
    support_u = seq(1L, 2L * s - 1L, by = 2L),
    support_v = seq(1L, 2L * s - 1L, by = 2L),
    tail_decay = 0.5,
    seed = 1L) {

  rank_type <- match.arg(rank_type)

  if (p != m) {
    stop("The current simulation design assumes p = m.")
  }
  if (max(support_u) > p || max(support_v) > m) {
    stop("The requested support exceeds the variable dimension.")
  }
  if (length(support_u) != s || length(support_v) != s) {
    stop("support_u and support_v must each contain exactly s indices.")
  }

  set.seed(seed)

  Sigma_x <- toeplitz_covariance(p, rho)
  Sigma_y <- toeplitz_covariance(m, rho)

  U0 <- .make_sparse_sigma_orthonormal(
    Sigma_x, support_u, r
  )
  V0 <- .make_sparse_sigma_orthonormal(
    Sigma_y, support_v, r
  )

  leading_lambda <- .default_leading_lambdas(r)

  if (rank_type == "exact") {
    lambda_full <- leading_lambda
    U_full <- U0
    V_full <- V0
  } else {
    q <- min(p, m)

    if (q <= r) {
      stop("Approximate low-rank setting requires min(p,m) > r.")
    }

    tail_index <- seq_len(q - r)
    lambda_tail <- leading_lambda[r] * exp(-tail_decay * tail_index)
    lambda_full <- c(leading_lambda, lambda_tail)

    U_full <- .complete_sigma_orthonormal_basis(
      Sigma_x, U0, seed = seed + 101L
    )
    V_full <- .complete_sigma_orthonormal_basis(
      Sigma_y, V0, seed = seed + 202L
    )
  }

  Sigma_xy <- Sigma_x %*% U_full %*%
    diag(lambda_full, length(lambda_full)) %*%
    t(V_full) %*% Sigma_y

  Sigma_joint <- rbind(
    cbind(Sigma_x, Sigma_xy),
    cbind(t(Sigma_xy), Sigma_y)
  )

  list(
    p = p,
    m = m,
    r = r,
    s = s,
    rank_type = rank_type,
    rho = rho,
    support_u = support_u,
    support_v = support_v,
    Sigma_x = Sigma_x,
    Sigma_y = Sigma_y,
    Sigma_xy = Sigma_xy,
    Sigma_joint = Sigma_joint,
    U0 = U0,
    V0 = V0,
    U_full = U_full,
    V_full = V_full,
    lambda = lambda_full,
    leading_lambda = leading_lambda,
    chol = chol(.sim_sym(Sigma_joint))
  )
}

.make_source_endpoint <- function(
    target_population,
    source_id,
    rho_min = 0.1,
    rho_max = 0.6,
    rho_gap = 0.1,
    lambda1_range = c(0.85, 0.95),
    lambda2_range = c(0.65, 0.85),
    c_range = c(0.5, 2),
    seed = 1L) {

  p <- target_population$p
  m <- target_population$m
  r <- target_population$r
  s <- target_population$s

  if (p != m) {
    stop("The current source design assumes p = m.")
  }

  # S^(k) is sampled from
  # [p] \ {p-5s+1,...,p} = {1,...,p-5s},
  # with |S^(k)| = s.
  if (p - 3L * s < s) {
    stop("The source-support design requires p >= 6s.")
  }

  set.seed(seed + 1000L * source_id)

  # Source-specific marginal covariance parameters:
  # rho_x,k, rho_y,k in [0.1,0.6], separated from rho_target by at least 0.1.
  sample_rho <- function() {
    repeat {
      rho <- runif(1L, rho_min, rho_max)
      if (abs(rho - target_population$rho) >= rho_gap) {
        return(rho)
      }
    }
  }

  rho_x_alt <- sample_rho()
  rho_y_alt <- sample_rho()

  Sigma_x_alt <- toeplitz_covariance(p, rho_x_alt)
  Sigma_y_alt <- toeplitz_covariance(m, rho_y_alt)

  # The same random row support is used for U^(k) and V^(k).
  # ---------------------------------------------------------
  # Far-away sparse directions.
  # ---------------------------------------------------------
  support_alt <- sort(
    sample(
      seq.int(
        from = 3L * s + 1L,
        to = p
      ),
      size = s,
      replace = FALSE
    )
  )
  
  U_far <- .make_sparse_sigma_orthonormal(
    Sigma_x_alt,
    support_alt,
    r
  )
  
  V_far <- .make_sparse_sigma_orthonormal(
    Sigma_y_alt,
    support_alt,
    r
  )
  
  # All non-informative sources share the same wrong sparse subspaces.
  
  set.seed(seed + 900002L)
  
  U_alt <- .make_sparse_sigma_orthonormal(
    Sigma_x_alt,
    support_alt,
    r
  )
  
  set.seed(seed + 900003L)
  
  V_alt <- .make_sparse_sigma_orthonormal(
    Sigma_y_alt,
    support_alt,
    r
  )
  
  support_u_alt <- support_alt
  support_v_alt <- support_alt
  
  # Source-specific leading canonical correlations.
  if (r != 2L) {
    stop("The current source-lambda design is written for r = 2.")
  }

  set.seed(
    seed + 1000L * source_id + 500L
  )
  
  lambda_alt_leading <- c(
    runif(1L, lambda1_range[1L], lambda1_range[2L]),
    runif(1L, lambda2_range[1L], lambda2_range[2L])
  )

  # With the default intervals lambda_1 > lambda_2 almost surely.
  if (lambda_alt_leading[2L] >= lambda_alt_leading[1L]) {
    lambda_alt_leading[2L] <- runif(
      1L,
      lambda2_range[1L],
      min(lambda2_range[2L], lambda_alt_leading[1L] - 1e-3)
    )
  }

  if (target_population$rank_type == "exact") {

    U_alt_full <- U_alt
    V_alt_full <- V_alt
    lambda_alt <- lambda_alt_leading
    c_alt <- NA_real_

  } else {

    q <- min(p, m)

    U_alt_full <- .complete_sigma_orthonormal_basis(
      Sigma_x_alt,
      U_alt,
      seed = seed + 1000L * source_id + 101L
    )

    V_alt_full <- .complete_sigma_orthonormal_basis(
      Sigma_y_alt,
      V_alt,
      seed = seed + 1000L * source_id + 202L
    )

    # Each source has its own exponential tail-decay parameter.
    c_alt <- runif(1L, c_range[1L], c_range[2L])
    tail_index <- seq_len(q - r)

    lambda_tail_alt <-
      lambda_alt_leading[r] *
      exp(-c_alt * tail_index)

    lambda_alt <- c(
      lambda_alt_leading,
      lambda_tail_alt
    )
  }

  Sigma_xy_alt <- Sigma_x_alt %*%
    U_alt_full %*%
    diag(lambda_alt, length(lambda_alt)) %*%
    t(V_alt_full) %*%
    Sigma_y_alt

  Sigma_alt <- .sim_sym(
    rbind(
      cbind(Sigma_x_alt, Sigma_xy_alt),
      cbind(t(Sigma_xy_alt), Sigma_y_alt)
    )
  )

  list(
    Sigma_alt = Sigma_alt,
    Sigma_x_alt = Sigma_x_alt,
    Sigma_y_alt = Sigma_y_alt,
    U_alt = U_alt,
    V_alt = V_alt,
    U_alt_full = U_alt_full,
    V_alt_full = V_alt_full,
    support_u_alt = support_u_alt,
    support_v_alt = support_v_alt,
    rho_x_alt = rho_x_alt,
    rho_y_alt = rho_y_alt,
    lambda_alt = lambda_alt,
    leading_lambda_alt = lambda_alt_leading,
    c_alt = c_alt
  )
}

.scale_source_to_distance <- function(
    Sigma0,
    Sigma_alt,
    distance,
    tolerance = 1e-8) {

  D <- .sim_sym(Sigma_alt - Sigma0)
  raw_distance <- .sim_spectral_norm_sym(D)

  if (!is.finite(raw_distance) || raw_distance <= 0) {
    stop("The source endpoint has zero or invalid distance from the target.")
  }
  if (distance > raw_distance * (1 + tolerance)) {
    stop(
      sprintf(
        "Requested source distance %.3f exceeds endpoint distance %.3f.",
        distance, raw_distance
      )
    )
  }

  weight <- distance / raw_distance

  Sigma_source <- (1 - weight) * Sigma0 + weight * Sigma_alt
  Sigma_source <- .sim_sym(Sigma_source)

  # Because 0 <= weight <= 1, this is a convex combination of two
  # positive-definite covariance matrices.
  achieved <- .sim_spectral_norm_sym(Sigma_source - Sigma0)

  list(
    Sigma = Sigma_source,
    chol = chol(Sigma_source),
    requested_distance = distance,
    achieved_distance = achieved,
    endpoint_distance = raw_distance,
    interpolation_weight = weight
  )
}

prepare_source_population_bank <- function(
    target_population,
    K = 10L,
    delta_informative = 0.05,
    noninformative_weight = 0.5,
    seed = 100L) {
  
  bank <- vector(
    "list",
    K
  )
  
  Sigma0 <- target_population$Sigma_joint
  
  for (k in seq_len(K)) {
    
    endpoint <- .make_source_endpoint(
      target_population,
      source_id = k,
      seed = seed
    )
    
    Sigma_alt <- endpoint$Sigma_alt
    
    D <- .sim_sym(
      Sigma_alt - Sigma0
    )
    
    endpoint_distance <-
      .sim_spectral_norm_sym(D)
    
    
    # ---------------------------------------------------------
    # Informative source:
    #
    # Sigma^(k)
    # =
    # Sigma^(0)
    # +
    # delta_informative *
    # (Sigma_tilde^(k)-Sigma^(0))
    # /
    # ||Sigma_tilde^(k)-Sigma^(0)||_2.
    # ---------------------------------------------------------
    informative_weight <-
      delta_informative /
      endpoint_distance
    
    if (informative_weight > 1) {
      stop(
        sprintf(
          paste0(
            "Source %d endpoint is too close to the target: ",
            "endpoint distance = %.4f < informative distance = %.4f."
          ),
          k,
          endpoint_distance,
          delta_informative
        )
      )
    }
    
    Sigma_informative <-
      (1 - informative_weight) * Sigma0 +
      informative_weight * Sigma_alt
    
    Sigma_informative <-
      .sim_sym(
        Sigma_informative
      )
    
    
    # ---------------------------------------------------------
    # Non-informative source:
    #
    # Sigma^(k)
    # =
    # (1-w) Sigma^(0)
    # +
    # w Sigma_tilde^(k),
    #
    # with w = 0.5.
    # ---------------------------------------------------------
    Sigma_noninformative <-
      (1 - noninformative_weight) * Sigma0 +
      noninformative_weight * Sigma_alt
    
    Sigma_noninformative <-
      .sim_sym(
        Sigma_noninformative
      )
    
    
    informative_distance <-
      .sim_spectral_norm_sym(
        Sigma_informative - Sigma0
      )
    
    noninformative_distance <-
      .sim_spectral_norm_sym(
        Sigma_noninformative - Sigma0
      )
    
    
    bank[[k]] <- list(
      
      endpoint = endpoint,
      
      informative = list(
        Sigma = Sigma_informative,
        chol = chol(Sigma_informative),
        requested_distance =
          delta_informative,
        achieved_distance =
          informative_distance,
        endpoint_distance =
          endpoint_distance,
        interpolation_weight =
          informative_weight
      ),
      
      noninformative = list(
        Sigma = Sigma_noninformative,
        chol = chol(Sigma_noninformative),
        achieved_distance =
          noninformative_distance,
        endpoint_distance =
          endpoint_distance,
        interpolation_weight =
          noninformative_weight
      )
    )
  }
  
  names(bank) <- paste0(
    "source_",
    seq_len(K)
  )
  
  bank
}
simulate_from_chol <- function(n, chol_sigma, p) {
  d <- ncol(chol_sigma)
  W <- matrix(rnorm(n * d), nrow = n, ncol = d) %*% chol_sigma

  list(
    X = W[, seq_len(p), drop = FALSE],
    Y = W[, p + seq_len(d - p), drop = FALSE]
  )
}

prediction_error_cca <- function(U_hat, U_true, Sigma_x) {
  Sigma_x <- .sim_sym(Sigma_x)

  M <- crossprod(U_hat, Sigma_x %*% U_true)
  M_svd <- svd(M)
  W <- M_svd$u %*% t(M_svd$v)

  difference <- U_hat %*% W - U_true

  max(
    0,
    sum(difference * (Sigma_x %*% difference))
  )
}

evaluate_cca_fit <- function(fit, population) {
  c(
    error_U = prediction_error_cca(
      fit$U, population$U0, population$Sigma_x
    ),
    error_V = prediction_error_cca(
      fit$V, population$V0, population$Sigma_y
    )
  )
}

source_detection_metrics <- function(selected, truth, K) {
  selected <- sort(unique(as.integer(selected)))
  truth <- sort(unique(as.integer(truth)))

  n_good <- length(truth)
  n_bad <- K - n_good

  tpr <- if (n_good == 0L) NA_real_ else
    length(intersect(selected, truth)) / n_good

  fpr <- if (n_bad == 0L) NA_real_ else
    length(setdiff(selected, truth)) / n_bad

  c(
    TPR = tpr,
    FPR = fpr,
    selected_count = length(selected)
  )
}

safe_trans_fit <- function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      structure(
        list(error_message = conditionMessage(e)),
        class = "trans_scca_failed"
      )
    }
  )
}

fit_failed <- function(x) {
  inherits(x, "trans_scca_failed")
}

summarize_simulation_results <- function(results) {
  key <- interaction(
    results$p,
    results$m,
    results$n,
    results$rank_type,
    results$K_A,
    results$method,
    drop = TRUE
  )

  pieces <- split(results, key)

  out <- do.call(rbind, lapply(pieces, function(d) {
    valid_u <- is.finite(d$error_U)
    valid_v <- is.finite(d$error_V)

    data.frame(
      p = d$p[1L],
      m = d$m[1L],
      n = d$n[1L],
      rank_type = d$rank_type[1L],
      K_A = d$K_A[1L],
      method = d$method[1L],
      mean_error_U = if (any(valid_u))
        mean(d$error_U[valid_u]) else NA_real_,
      sd_error_U = if (sum(valid_u) > 1L)
        sd(d$error_U[valid_u]) else NA_real_,
      mean_error_V = if (any(valid_v))
        mean(d$error_V[valid_v]) else NA_real_,
      sd_error_V = if (sum(valid_v) > 1L)
        sd(d$error_V[valid_v]) else NA_real_,
      success_rate = mean(valid_u & valid_v),
      stringsAsFactors = FALSE
    )
  }))

  rownames(out) <- NULL
  out[order(out$p, out$n, out$rank_type, out$K_A, out$method), ]
}

plot_simulation_summary <- function(
    summary_table,
    output_dir = "simulation_plots",
    prefix = "split") {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  rank_types <- unique(summary_table$rank_type)
  settings <- unique(summary_table[c("p", "m", "n")])
  methods <- c("Target-only", "Pooled-All", "Oracle", "Detected")

  line_types <- c(1, 2, 3, 4)
  point_types <- c(16, 17, 15, 18)

  for (rank_type in rank_types) {
    for (metric in c("U", "V")) {
      filename <- file.path(
        output_dir,
        sprintf("%s_%s_error_%s.pdf", prefix, rank_type, metric)
      )

      grDevices::pdf(filename, width = 10, height = 8)
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par), add = TRUE)

      par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1.5))

      for (ii in seq_len(nrow(settings))) {
        pp <- settings$p[ii]
        mm <- settings$m[ii]
        nn <- settings$n[ii]

        d <- summary_table[
          summary_table$p == pp &
            summary_table$m == mm &
            summary_table$n == nn &
            summary_table$rank_type == rank_type,
          ,
          drop = FALSE
        ]

        mean_name <- paste0("mean_error_", metric)

        y_range <- range(d[[mean_name]], finite = TRUE)
        if (!all(is.finite(y_range))) y_range <- c(0, 1)

        first <- TRUE

        for (jj in seq_along(methods)) {
          method <- methods[jj]
          dm <- d[d$method == method, , drop = FALSE]
          dm <- dm[order(dm$K_A), , drop = FALSE]

          if (!nrow(dm)) next

          if (first) {
            plot(
              dm$K_A,
              dm[[mean_name]],
              type = "o",
              lty = line_types[jj],
              pch = point_types[jj],
              ylim = y_range,
              xlab = expression(K[A]),
              ylab = sprintf("Prediction error of %s", metric),
              main = sprintf("(p,m,n)=(%d,%d,%d)", pp, mm, nn)
            )
            first <- FALSE
          } else {
            lines(
              dm$K_A,
              dm[[mean_name]],
              type = "o",
              lty = line_types[jj],
              pch = point_types[jj]
            )
          }
        }

        legend(
          "topright",
          legend = methods,
          lty = line_types,
          pch = point_types,
          bty = "n",
          cex = 0.8
        )
      }

      grDevices::dev.off()
      on.exit(NULL, add = FALSE)
    }
  }

  invisible(TRUE)
}
