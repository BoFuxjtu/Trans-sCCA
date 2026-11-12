# Trans sparse CCA estimator
#
# This file implements Algorithm 2 in five stages:
#   1. Gao--Ma--Zhou convex initialization by ADMM or CVXR;
#   2. pooled group-Lasso transfer;
#   3. target-only group-Lasso candidate;
#   4. target-loss aggregation;
#   5. spectrally thresholded target normalization.
# Packages do the two generic optimization tasks:
#   CVXR   -> Gao--Ma--Zhou convex initialization;
#   glmnet -> pooled and target-only row-group Lasso.

#library(glmnet)
.sym <- function(A) (A + t(A)) / 2

.inv_root <- function(A, ridge = 0, tol = 1e-10) {
  ee <- eigen(.sym(A), symmetric = TRUE)
  d <- ee$values + ridge
  if (min(d) <= tol * max(1, max(abs(d)))) {
    stop("A normalization matrix is singular; tune the penalties or add ridge.")
  }
  tcrossprod(sweep(ee$vectors, 2, 1 / sqrt(d), `*`), ee$vectors)
}

.row_support_size <- function(A, tol = 1e-8) {
  row_norms <- sqrt(rowSums(A^2))
  
  threshold <- tol * max(1, max(row_norms))
  
  sum(row_norms > threshold)
}

# Optional centering is performed separately within each domain and fold.
.prepare_fold <- function(domain, rows, center) {
  X <- as.matrix(domain$X[rows, , drop = FALSE])
  Y <- as.matrix(domain$Y[rows, , drop = FALSE])
  if (center) {
    X <- sweep(X, 2, colMeans(X), `-`)
    Y <- sweep(Y, 2, colMeans(Y), `-`)
  }
  #if (center) {
  #  if (nrow(X) < 2) stop("Every centered domain-fold needs at least two rows.")
  #  X <- sweep(X, 2, colMeans(X), `-`)
   # Y <- sweep(Y, 2, colMeans(Y), `-`)
  #  correction <- sqrt(nrow(X) / (nrow(X) - 1))
   # X <- correction * X
   # Y <- correction * Y
  #}
  list(X = X, Y = Y)
}

.stack_domains <- function(domains) {
  list(
    X = do.call(rbind, lapply(domains, `[[`, "X")),
    Y = do.call(rbind, lapply(domains, `[[`, "Y"))
  )
}

.moments <- function(data) {
  n <- nrow(data$X)
  list(
    Sx = crossprod(data$X) / n,
    Sy = crossprod(data$Y) / n,
    Sxy = crossprod(data$X, data$Y) / n,
    n = n
  )
}

# Initialization step, delegated to CVXR.  eta is exactly the coefficient of ||F||_1
# in the displayed objective.
.scca_initialization <- function(
    Sx, Sy, Sxy, r, eta,
    method = c("ADMM", "CVXR"),
    rho = 2, tol = 1e-5, maxit = 500,
    inner_tol = 1e-7, inner_maxit = 200,
    solver = "SCS", solver_options = list(),
    verbose = FALSE) {
  
  method <- match.arg(method)
  sym <- function(A) (A + t(A)) / 2
  fnorm <- function(A) sqrt(sum(A^2))
  soft <- function(A, a) sign(A) * pmax(abs(A) - a, 0)
  
  root <- function(A) {
    ee <- eigen(sym(A), symmetric = TRUE)
    if (min(ee$values) < -1e-8 * max(1, max(abs(ee$values))))
      stop("Covariance matrix is not PSD.")
    d <- sqrt(pmax(ee$values, 0))
    tcrossprod(sweep(ee$vectors, 2, d, "*"), ee$vectors)
  }
  
  # Projection of singular values onto
  # {0 <= d_j <= 1, sum d_j <= r}.
  proj_sv <- function(d) {
    z <- pmin(1, pmax(d, 0))
    if (sum(z) <= r) return(z)
    
    lo <- 0
    hi <- max(d)
    while (hi - lo > 1e-10 * max(1, hi)) {
      a <- (lo + hi) / 2
      z <- pmin(1, pmax(d - a, 0))
      if (sum(z) <= r) hi <- a else lo <- a
    }
    pmin(1, pmax(d - hi, 0))
  }
  
  svcst <- function(A) {
    ss <- svd(A)
    d <- proj_sv(ss$d)
    ss$u %*% diag(d, length(d)) %*% t(ss$v)
  }
  
  Rx <- root(Sx)
  Ry <- root(Sy)
  
  ## ---------------------------------------------------------
  ## CVXR
  ## ---------------------------------------------------------
  if (method == "CVXR") {
    if (!requireNamespace("CVXR", quietly = TRUE))
      stop("Install package 'CVXR'.")
    
    p <- nrow(Sxy)
    m <- ncol(Sxy)
    F <- CVXR::Variable(c(p, m))
    G <- Rx %*% F %*% Ry
    
    problem <- CVXR::Problem(
      CVXR::Maximize(
        CVXR::sum_entries(Sxy * F) -
          eta * CVXR::sum_entries(abs(F))
      ),
      list(
        CVXR::norm_nuc(G) <= r,
        CVXR::sigma_max(G) <= 1
      )
    )
    
    args <- utils::modifyList(
      list(problem = problem, solver = solver, verbose = verbose),
      solver_options
    )
    do.call(CVXR::psolve, args)
    
    status <- CVXR::status(problem)
    if (!status %in% c("optimal", "optimal_inaccurate"))
      stop("CVXR failed: ", status)
    
    return(list(
      F = as.matrix(CVXR::value(F)),
      method = "CVXR",
      status = status
    ))
  }
  
  ## ---------------------------------------------------------
  ## ADMM
  ## ---------------------------------------------------------
  Aop  <- function(F) Rx %*% F %*% Ry
  Aadj <- function(G) Rx %*% G %*% Ry
  
  L <- rho *
    max(eigen(sym(Sx), symmetric = TRUE, only.values = TRUE)$values) *
    max(eigen(sym(Sy), symmetric = TRUE, only.values = TRUE)$values)
  
  # F-subproblem by FISTA.
  update_F <- function(F, G, H) {
    Z <- F
    t0 <- 1
    
    for (j in seq_len(inner_maxit)) {
      grad <- -Sxy + Aadj(H + rho * (Aop(Z) - G))
      Fnew <- soft(Z - grad / L, eta / L)
      
      if (fnorm(Fnew - F) / max(1, fnorm(F)) <= inner_tol)
        return(Fnew)
      
      t1 <- (1 + sqrt(1 + 4 * t0^2)) / 2
      Z <- Fnew + (t0 - 1) / t1 * (Fnew - F)
      F <- Fnew
      t0 <- t1
    }
    F
  }
  
  F <- svcst(Sxy)
  G <- matrix(0, nrow(Sxy), ncol(Sxy))
  H <- G
  
  converged <- FALSE
  
  for (iter in seq_len(maxit)) {
    Gold <- G
    
    F <- update_F(F, G, H)
    
    AF <- Aop(F)
    G <- svcst(AF + H / rho)
    
    R <- AF - G
    H <- H + rho * R
    
    primal <- fnorm(R)
    dual <- rho * fnorm(Aadj(G - Gold))
    
    if (verbose && (iter == 1 || iter %% 10 == 0))
      message(sprintf("ADMM %d: primal %.3e, dual %.3e",
                      iter, primal, dual))
    
    if (primal <= tol * max(1, fnorm(AF), fnorm(G)) &&
        dual <= tol * max(1, fnorm(Aadj(H)))) {
      converged <- TRUE
      break
    }
  }
  
  list(
    F = F,
    method = "ADMM",
    converged = converged,
    iterations = iter,
    primal = primal,
    dual = dual
  )
}

# glmnet solves 1/(2n)||R-XL||_F^2 + lambda sum_j ||L_j||_2.
# Hence lambda = eta/2 matches the paper's
# 1/n||R-XL||_F^2 + eta ||L||_{2,1} objective.
.glmnet_group_lasso <- function(X, R, eta) {
  lambda <- eta / 2
  if (lambda <= 0) stop("The glmnet penalty must be positive.")
  lambda_path <- lambda * c(4, 2, 1)
  
  if (ncol(R) == 1L) {
    fit <- glmnet::glmnet(
      x = X, y = drop(R), family = "gaussian", alpha = 1,
      lambda = lambda_path, intercept = FALSE, standardize = FALSE,
      thresh = 1e-10
    )
    L <- as.matrix(stats::coef(fit, s = lambda))[-1, , drop = FALSE]
  } else {
    fit <- glmnet::glmnet(
      x = X, y = R, family = "mgaussian", alpha = 1,
      lambda = lambda_path, intercept = FALSE, standardize = FALSE,
      standardize.response = FALSE, thresh = 1e-10
    )
    coefficients <- stats::coef(fit, s = lambda)
    L <- do.call(cbind, lapply(coefficients, function(B) {
      as.matrix(B)[-1, 1]
    }))
  }
  list(L = as.matrix(L), fit = fit)
}

.target_loss <- function(L, Sx, Sxy, V) {
  sum(L * (Sx %*% L)) - 2 * sum(L * (Sxy %*% V))
}

.normalize_side <- function(L, S_pool, S_target, epsilon, ridge = 0) {
  B <- crossprod(L, S_pool %*% L)
  Z <- L %*% .inv_root(B, ridge)
  
  H <- .sym(
    crossprod(Z, S_target %*% Z) - diag(ncol(L))
  )
  
  ee <- eigen(H, symmetric = TRUE)
  d <- ifelse(abs(ee$values) > 2 * epsilon, ee$values, 0)
  
  Ht <- tcrossprod(
    sweep(ee$vectors, 2, d, `*`),
    ee$vectors
  )
  
  list(
    final = Z %*% .inv_root(diag(ncol(L)) + Ht, ridge),
    B = B,
    Z = Z,
    H = H,
    H_thresholded = Ht
  )
}

.obs_target_loss <- function(A, B, L, R) {
  AL <- A %*% L
  BR <- B %*% R
  rowSums(AL^2) - 2 * rowSums(AL * BR)
}

.validation_cc_score <- function(
    U,
    V,
    validation,
    r,
    ridge = 1e-8) {
  
  XU <- validation$X %*% U
  YV <- validation$Y %*% V
  
  n <- nrow(validation$X)
  
  Bx <- crossprod(XU) / n
  By <- crossprod(YV) / n
  Cxy <- crossprod(XU, YV) / n
  
  C <- tryCatch(
    .inv_root(
      Bx,
      ridge = ridge
    ) %*%
      Cxy %*%
      .inv_root(
        By,
        ridge = ridge
      ),
    error = function(e) NULL
  )
  
  if (is.null(C)) {
    return(-Inf)
  }
  
  rho <- svd(
    C,
    nu = 0,
    nv = 0
  )$d
  
  sum(
    rho[
      seq_len(
        min(r, length(rho))
      )
    ]
  )
}

.cv_scca_initialization <- function(
    fold1,
    r,
    grid = c(0.5, 1, 1.5, 2),
    nfolds = 3,
    seed = 1001,
    method = "CVXR",
    rho = 2,
    tol = 1e-5,
    maxit = 500,
    inner_tol = 1e-7,
    inner_maxit = 200,
    solver = "SCS",
    solver_options = list(),
    verbose = FALSE) {
  
  target1 <- fold1[[1]]
  sources1 <- fold1[-1]
  
  n0 <- nrow(target1$X)
  p <- ncol(target1$X)
  m <- ncol(target1$Y)
  
  if (nfolds < 2 || nfolds > n0)
    stop("Invalid number of initialization CV folds.")
  
  set.seed(seed)
  cv_index <- split(
    sample.int(n0),
    rep(seq_len(nfolds), length.out = n0)
  )
  
  cv_score <- numeric(length(grid))
  
  for (g in seq_along(grid)) {
    
    b <- grid[g]
    score <- 0
    
    for (ell in seq_len(nfolds)) {
      
      val_id <- cv_index[[ell]]
      train_id <- setdiff(seq_len(n0), val_id)
      
      target_train <- list(
        X = target1$X[train_id, , drop = FALSE],
        Y = target1$Y[train_id, , drop = FALSE]
      )
      
      target_val <- list(
        X = target1$X[val_id, , drop = FALSE],
        Y = target1$Y[val_id, , drop = FALSE]
      )
      
      pooled_train <- .stack_domains(
        c(list(target_train), sources1)
      )
      
      Mtrain <- .moments(pooled_train)
      
      eta <- b * sqrt(
        log(p + m) / nrow(pooled_train$X)
      )
      
      init <- .scca_initialization(
        Mtrain$Sx,
        Mtrain$Sy,
        Mtrain$Sxy,
        r = r,
        eta = eta,
        method = method,
        rho = rho,
        tol = tol,
        maxit = maxit,
        inner_tol = inner_tol,
        inner_maxit = inner_maxit,
        solver = solver,
        solver_options = solver_options,
        verbose = FALSE
      )
      
      ss <- svd(init$F, nu = r, nv = r)
      
      if (ss$d[r] <= 1e-8 * max(ss$d)) {
        next
      }
      
      U <- ss$u[, seq_len(r), drop = FALSE]
      V <- ss$v[, seq_len(r), drop = FALSE]
      
      cc <- .validation_cc_score(
        U,
        V,
        target_val,
        r
      )
      
      score <- score +
        length(val_id) * cc
    }
    
    cv_score[g] <- score / n0
  }
  
  best <- which.max(cv_score)
  b_initial <- grid[best]
  
  N1 <- sum(vapply(
    fold1,
    function(z) nrow(z$X),
    integer(1)
  ))
  
  eta_initial <- b_initial *
    sqrt(log(p + m) / N1)
  
  list(
    b_initial = b_initial,
    eta_initial = eta_initial,
    cv_table = data.frame(
      b = grid,
      score = cv_score
    )
  )
}

.cv_trans_scca_penalties <- function(
    target,
    sources,
    r,
    c_initial = 0.55,
    grid = c(0.5, 1, 1.5, 2),
    nfolds = 5,
    seed = 1001,
    center = FALSE,
    initialization_method = "ADMM",
    rho = 2,
    tol = 1e-5,
    maxit = 500,
    inner_tol = 1e-5,
    inner_maxit = 50,
    solver = "SCS",
    solver_options = list()) {
  
  n0 <- nrow(target$X)
  p <- ncol(target$X)
  m <- ncol(target$Y)
  
  complexity_u <- r + log(p)
  complexity_v <- r + log(m)
  
  sources_full <- lapply(
    sources,
    function(z) {
      .prepare_fold(
        z,
        seq_len(nrow(z$X)),
        center
      )
    }
  )
  
  set.seed(seed)
  
  cv_index <- split(
    sample.int(n0),
    rep(seq_len(nfolds), length.out = n0)
  )
  
  cv_transfer <- numeric(length(grid))
  cv_target <- numeric(length(grid))
  
  for (ell in seq_len(nfolds)) {
    
    val_id <- cv_index[[ell]]
    train_id <- setdiff(seq_len(n0), val_id)
    
    target_train <- .prepare_fold(
      target,
      train_id,
      center
    )
    
    target_val <- .prepare_fold(
      target,
      val_id,
      center
    )
    
    pooled_train <- .stack_domains(
      c(
        list(target_train),
        sources_full
      )
    )
    
    n_transfer <- nrow(pooled_train$X)
    n_target <- nrow(target_train$X)
    
    # Fixed 0.55 initialization: NOT tuned.
    Mtrain <- .moments(pooled_train)
    
    eta_initial_cv <- c_initial *
      sqrt(log(p + m) / n_transfer)
    
    init_cv <- .scca_initialization(
      Mtrain$Sx,
      Mtrain$Sy,
      Mtrain$Sxy,
      r = r,
      eta = eta_initial_cv,
      method = initialization_method,
      rho = rho,
      tol = tol,
      maxit = maxit,
      inner_tol = inner_tol,
      inner_maxit = inner_maxit,
      solver = solver,
      solver_options = solver_options,
      verbose = FALSE
    )
    
    ss <- svd(
      init_cv$F,
      nu = r,
      nv = r
    )
    
    if (length(ss$d) < r ||
        ss$d[r] <= 1e-8 * max(ss$d)) {
      next
    }
    
    U_initial_cv <-
      ss$u[, seq_len(r), drop = FALSE]
    
    V_initial_cv <-
      ss$v[, seq_len(r), drop = FALSE]
    
    for (g in seq_along(grid)) {
      
      b <- grid[g]
      
      eta_tu <- b *
        sqrt(complexity_u / n_transfer)
      
      eta_tv <- b *
        sqrt(complexity_v / n_transfer)
      
      eta_0u <- b *
        sqrt(complexity_u / n_target)
      
      eta_0v <- b *
        sqrt(complexity_v / n_target)
      
      U_transfer <- .glmnet_group_lasso(
        pooled_train$X,
        pooled_train$Y %*% V_initial_cv,
        eta_tu
      )$L
      
      V_transfer <- .glmnet_group_lasso(
        pooled_train$Y,
        pooled_train$X %*% U_initial_cv,
        eta_tv
      )$L
      
      U_target <- .glmnet_group_lasso(
        target_train$X,
        target_train$Y %*% V_initial_cv,
        eta_0u
      )$L
      
      V_target <- .glmnet_group_lasso(
        target_train$Y,
        target_train$X %*% U_initial_cv,
        eta_0v
      )$L
      
      cc_transfer <- .validation_cc_score(
        U_transfer,
        V_transfer,
        target_val,
        r
      )
      
      cc_target <- .validation_cc_score(
        U_target,
        V_target,
        target_val,
        r
      )
      
      if (is.finite(cc_transfer)) {
        cv_transfer[g] <- cv_transfer[g] +
          length(val_id) * cc_transfer
      }
      
      if (is.finite(cc_target)) {
        cv_target[g] <- cv_target[g] +
          length(val_id) * cc_target
      }
    }
  }
  
  cv_transfer <- cv_transfer / n0
  cv_target <- cv_target / n0
  
  best_transfer <- which.max(cv_transfer)
  best_target <- which.max(cv_target)
  
  list(
    b_transfer = grid[best_transfer],
    b_target = grid[best_target],
    cv_table = data.frame(
      b = grid,
      transfer = cv_transfer,
      target = cv_target
    )
  )
}


detect_trans_scca_sources <- function(
    target, sources, r,
    repeats = 1,
    c_initial = 0.55,
    c_u = 1,
    c_v = 1,
    c_det = 1,                   
    rule = c("sum", "both"),     
    top_n = NULL,
    center = FALSE,
    seed = 1,
    initialization_method = c("ADMM", "CVXR"),
    ...) {
  
  rule <- match.arg(rule)
  initialization_method <- match.arg(initialization_method)
  
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Install package 'glmnet'.")
  
  target <- list(
    X = as.matrix(target$X),
    Y = as.matrix(target$Y)
  )
  
  sources <- lapply(sources, function(z) {
    list(
      X = as.matrix(z$X),
      Y = as.matrix(z$Y)
    )
  })
  
  n0 <- nrow(target$X)
  p <- ncol(target$X)
  m <- ncol(target$Y)
  K <- length(sources)
  
  if (K == 0L)
    return(list(
      informative = integer(0),
      source_order = integer(0),
      selected_q = 0L
    ))
  
  # ==========================================================
  # Use the full target sample for cyclic 3-fold cross-fitting.
  # ==========================================================
  
  target_train_raw <- target
  n_train <- n0
  
  source_data <- lapply(sources, function(z) {
    .prepare_fold(
      z,
      seq_len(nrow(z$X)),
      center
    )
  })
  
  
  # ==========================================================
  # Step 1: use the canonical correlations to rank sources
  # ==========================================================
  
  gap_cc_sum <- numeric(K)
  gap_cc_weight <- numeric(K)
  
  cyclic_index <- function(ell, a) {
    (ell + a - 1L) %% 3L + 1L
  }
  
  for (b in seq_len(repeats)) {
    
    set.seed(seed + b)
    
    folds <- split(
      sample.int(n_train),
      rep(1:3, length.out = n_train)
    )
    
    for (ell in 1:3) {
      
      validation_index <- folds[[ell]]
      pilot_index <- folds[[cyclic_index(ell, 1L)]]
      training_index <- folds[[cyclic_index(ell, 2L)]]
      
      validation <- .prepare_fold(
        target_train_raw,
        validation_index,
        center
      )
      
      pilot <- .prepare_fold(
        target_train_raw,
        pilot_index,
        center
      )
      
      training <- .prepare_fold(
        target_train_raw,
        training_index,
        center
      )
      
      pilot_moments <- .moments(pilot)
      
      eta_initial <- c_initial *
        sqrt(log(p + m) / nrow(pilot$X))
      
      initialization <- .scca_initialization(
        pilot_moments$Sx,
        pilot_moments$Sy,
        pilot_moments$Sxy,
        r = r,
        eta = eta_initial,
        method = initialization_method,
        ...
      )
      
      ss <- svd(
        initialization$F,
        nu = r,
        nv = r
      )
      
      if (length(ss$d) < r ||
          ss$d[r] <= 1e-8 * max(ss$d))
        stop("A target pilot has numerical rank below r.")
      
      U_pilot <- ss$u[, seq_len(r), drop = FALSE]
      V_pilot <- ss$v[, seq_len(r), drop = FALSE]
      
      eta_u0 <- c_u *
        sqrt((r + log(p)) / nrow(training$X))
      
      eta_v0 <- c_v *
        sqrt((r + log(m)) / nrow(training$Y))
      
      U0 <- .glmnet_group_lasso(
        training$X,
        training$Y %*% V_pilot,
        eta_u0
      )$L
      
      V0 <- .glmnet_group_lasso(
        training$Y,
        training$X %*% U_pilot,
        eta_v0
      )$L
      
      cc0 <- .validation_cc_score(
        U0,
        V0,
        validation,
        r
      )
      
      for (k in seq_len(K)) {
        
        pairwise <- .stack_domains(
          list(
            training,
            source_data[[k]]
          )
        )
        
        nk <- nrow(pairwise$X)
        
        eta_uk <- c_u *
          sqrt((r + log(p)) / nk)
        
        eta_vk <- c_v *
          sqrt((r + log(m)) / nk)
        
        Uk <- .glmnet_group_lasso(
          pairwise$X,
          pairwise$Y %*% V_pilot,
          eta_uk
        )$L
        
        Vk <- .glmnet_group_lasso(
          pairwise$Y,
          pairwise$X %*% U_pilot,
          eta_vk
        )$L
        
        cck <- .validation_cc_score(
          Uk,
          Vk,
          validation,
          r
        )
        
        if (is.finite(cc0) &&
            is.finite(cck)) {
          
          gap_cc_sum[k] <-
            gap_cc_sum[k] +
            length(validation_index) *
            (cc0 - cck)
          
          gap_cc_weight[k] <-
            gap_cc_weight[k] +
            length(validation_index)
        }
      }
    }
  }
  
  gap_cc <- ifelse(
    gap_cc_weight > 0,
    gap_cc_sum / gap_cc_weight,
    Inf
  )
  
  source_order <- order(
    gap_cc
  )
  
  # ==========================================================
  # Step 2: fit top-q, q = 0,...,K
  # ==========================================================
  
  q_max <- if (is.null(top_n)) {
    K
  } else {
    min(as.integer(top_n), K)
  }
  
  q_grid <- 0:q_max
  
  cc_score_sum <- numeric(length(q_grid))
  cc_score_weight <- numeric(length(q_grid))
  
  # One additional cyclic 3-fold split for selecting q.
  set.seed(seed + repeats + 1000L)
  
  folds_q <- split(
    sample.int(n0),
    rep(1:3, length.out = n0)
  )
  
  for (ell in 1:3) {
    
    validation_index <- folds_q[[ell]]
    pilot_index <- folds_q[[cyclic_index(ell, 1L)]]
    training_index <- folds_q[[cyclic_index(ell, 2L)]]
    
    validation <- .prepare_fold(
      target,
      validation_index,
      center
    )
    
    pilot <- .prepare_fold(
      target,
      pilot_index,
      center
    )
    
    training <- .prepare_fold(
      target,
      training_index,
      center
    )
    
    # --------------------------------------------------------
    # Target-only initialization on fold ell oplus 1.
    # This initializer is shared by all q.
    # --------------------------------------------------------
    
    pilot_moments <- .moments(pilot)
    
    eta_initial_q <- c_initial *
      sqrt(log(p + m) / nrow(pilot$X))
    
    init_q <- tryCatch(
      .scca_initialization(
        pilot_moments$Sx,
        pilot_moments$Sy,
        pilot_moments$Sxy,
        r = r,
        eta = eta_initial_q,
        method = initialization_method,
        ...
      ),
      error = function(e) NULL
    )
    
    if (is.null(init_q))
      next
    
    ss_q <- svd(
      init_q$F,
      nu = r,
      nv = r
    )
    
    if (length(ss_q$d) < r ||
        ss_q$d[r] <= 1e-8 * max(ss_q$d))
      next
    
    U_initial <- ss_q$u[, seq_len(r), drop = FALSE]
    V_initial <- ss_q$v[, seq_len(r), drop = FALSE]
    
    # --------------------------------------------------------
    # Fit all top-q candidates on fold ell oplus 2.
    # --------------------------------------------------------
    
    for (ii in seq_along(q_grid)) {
      
      q <- q_grid[ii]
      
      active <- if (q == 0L) {
        integer(0)
      } else {
        source_order[seq_len(q)]
      }
      
      domains_q <- c(
        list(training),
        source_data[active]
      )
      
      pooled_q <- .stack_domains(domains_q)
      Nq <- nrow(pooled_q$X)
      
      eta_uq <- c_u *
        sqrt((r + log(p)) / Nq)
      
      eta_vq <- c_v *
        sqrt((r + log(m)) / Nq)
      
      Uq <- .glmnet_group_lasso(
        pooled_q$X,
        pooled_q$Y %*% V_initial,
        eta_uq
      )$L
      
      Vq <- .glmnet_group_lasso(
        pooled_q$Y,
        pooled_q$X %*% U_initial,
        eta_vq
      )$L
      
      ccq <- .validation_cc_score(
        Uq,
        Vq,
        validation,
        r
      )
      
      if (is.finite(ccq)) {
        
        cc_score_sum[ii] <-
          cc_score_sum[ii] +
          length(validation_index) * ccq
        
        cc_score_weight[ii] <-
          cc_score_weight[ii] +
          length(validation_index)
      }
    }
  }
  
  cc_score <- ifelse(
    cc_score_weight > 0,
    cc_score_sum / cc_score_weight,
    -Inf
  )
  
  
  # ==========================================================
  # Step 3: choose q giving largest target canonical correlation
  # ==========================================================
  
  if (all(!is.finite(cc_score))) {
    selected_q <- 0L
  } else {
    selected_q <- q_grid[
      which.max(cc_score)
    ]
  }
  
  informative <- if (selected_q == 0L) {
    integer(0)
  } else {
    sort(
      source_order[
        seq_len(selected_q)
      ]
    )
  }
  
  selected <- seq_len(K) %in% informative
  
  
  scores <- data.frame(
    source = seq_len(K),
    
    cc_gap = gap_cc,
    
    cc_improvement = -gap_cc,
    
    rank = match(
      seq_len(K),
      source_order
    ),
    
    selected = selected
  )
  
  rank_path <- lapply(
    0:K,
    function(q) {
      head(source_order, q)
    }
  )
  
  names(rank_path) <-
    paste0("top_", 0:K)
  
  
  q_scores <- data.frame(
    q = q_grid,
    canonical_correlation_score = cc_score
  )
  
  
  list(
    informative = informative,
    selected_q = selected_q,
    
    source_order = source_order,
    
    ranked_scores =
      scores[
        source_order,
        ,
        drop = FALSE
      ],
    
    scores = scores,
    
    rank_path = rank_path,
    
    q_scores = q_scores
  )
}
#' Compact Trans sparse CCA estimator for U.
#'
#' The theoretical tuning coefficients can be supplied explicitly.  If they
#' are omitted, s_u and h are used to construct their rate-calibrated defaults.
trans_scca <- function(
    target,
    sources = list(),
    informative = seq_along(sources),
    r,
    s_u = NULL,
    s_v = NULL,
    h = NULL,
    eta_initial = NULL,
    eta_transfer = NULL,
    eta_target = NULL,
    aggregation_slack = NULL,
    epsilon_target = NULL,
    eta_transfer_v = NULL,
    eta_target_v = NULL,
    aggregation_slack_v = NULL,
    epsilon_target_v = NULL,
    c_initial = 0.55,
    c_transfer = 1,
    c_target = 1,
    c_aggregation = 1,
    c_normalization = 1,
    center = FALSE,
    seed = 1,
    initialization_method = c("ADMM", "CVXR"),
    rho_admm = 2,
    admm_tol = 1e-5,
    admm_maxit = 500,
    admm_inner_tol = 1e-5,
    admm_inner_maxit = 50,
    tune_initialization = FALSE,
    tune_penalties = FALSE,
    cv_grid = c(0.5, 1, 1.5, 2),
    cv_folds = 5,
    cv_seed = NULL,
    cvxr_solver = "SCS",
    cvxr_options = list(),
    normalization_ridge = 0,
    verbose = FALSE,
    aggregation_mode = c(
      "adaptive",
      "transfer",
      "target"
    ),
    normalization_mode = c(
      "truncated",
      "target"
    )) {
  aggregation_mode <- match.arg(aggregation_mode)
  normalization_mode <- match.arg(normalization_mode)
  initialization_method <- match.arg(initialization_method)
  if (is.null(cv_seed))
    cv_seed <- seed + 1000
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Install package 'glmnet'.")
  
  target$X <- as.matrix(target$X)
  target$Y <- as.matrix(target$Y)
  p <- ncol(target$X)
  m <- ncol(target$Y)
  if (r < 1 || r > min(p, m)) stop("Invalid r.")
  if (!is.null(s_u) &&
      (length(s_u) != 1L || !is.finite(s_u) || s_u < r || s_u > p)) {
    stop("s_u must be one finite number between r and p.")
  }
  if (!is.null(s_v) &&
      (length(s_v) != 1L || !is.finite(s_v) || s_v < r || s_v > m)) {
    stop("s_v must be one finite number between r and m.")
  }
  if (!is.null(h) && (length(h) != 1L || !is.finite(h) || h < 0)) {
    stop("h must be one finite nonnegative number.")
  }
  if (!is.logical(tune_penalties) ||
      length(tune_penalties) != 1L ||
      is.na(tune_penalties)) {
    stop("tune_penalties must be TRUE or FALSE.")
  }
  
  if (tune_penalties &&
      (length(cv_folds) != 1L ||
       !is.finite(cv_folds) ||
       cv_folds < 2 ||
       cv_folds != as.integer(cv_folds))) {
    stop("cv_folds must be an integer at least 2.")
  }
  
  if (tune_penalties &&
      (length(cv_grid) == 0L ||
       any(!is.finite(cv_grid)) ||
       any(cv_grid <= 0))) {
    stop("cv_grid must contain positive finite values.")
  }
  if (length(informative)) {
    if (any(informative != as.integer(informative)) ||
        any(informative < 1) ||
        any(informative > length(sources)) ||
        anyDuplicated(informative)) {
      stop("Invalid informative source indices.")
    }
    
    informative <- as.integer(informative)
  }
  selected_sources <- lapply(informative, function(k) {
    list(X = as.matrix(sources[[k]]$X), Y = as.matrix(sources[[k]]$Y))
  })
  domains <- c(list(target), selected_sources)
  if (any(vapply(domains, function(d) ncol(d$X) != p || ncol(d$Y) != m ||
                 nrow(d$X) != nrow(d$Y), logical(1)))) {
    stop("All selected domains must have compatible paired X and Y matrices.")
  }
  
  # ---------------------------------------------------------
  # Full-target CV for the second-stage penalty multipliers.
  # This is performed before the theoretical 3-way split.
  # ---------------------------------------------------------
  cv_tuning <- NULL
  b_transfer <- c_transfer
  b_target <- c_target
  
  if (tune_penalties) {
    
    cv_tuning <- .cv_trans_scca_penalties(
      target = target,
      sources = selected_sources,
      r = r,
      c_initial = c_initial,
      grid = cv_grid,
      nfolds = cv_folds,
      seed = cv_seed,
      center = center,
      initialization_method =
        initialization_method,
      rho = rho_admm,
      tol = admm_tol,
      maxit = admm_maxit,
      inner_tol = admm_inner_tol,
      inner_maxit = admm_inner_maxit,
      solver = cvxr_solver,
      solver_options = cvxr_options
    )
    
    b_transfer <-
      cv_tuning$b_transfer
    
    b_target <-
      cv_tuning$b_target
  }
  
  set.seed(seed)
  folds <- lapply(domains, function(d) {
    split(sample.int(nrow(d$X)), rep(1:3, length.out = nrow(d$X)))
  })
  get_fold <- function(j) {
    lapply(seq_along(domains), function(k) {
      .prepare_fold(domains[[k]], folds[[k]][[j]], center)
    })
  }
  fold1 <- get_fold(1)
  fold2 <- get_fold(2)
  fold3 <- get_fold(3)
  pooled1 <- .stack_domains(fold1)
  pooled2 <- .stack_domains(fold2)
  pooled3 <- .stack_domains(fold3)
  target2 <- fold2[[1]]
  target3 <- fold3[[1]]
  
  # Initialization tuning parameter.
  cv_initialization <- NULL
  
  if (tune_initialization) {
    
    cv_initialization <- .cv_scca_initialization(
      fold1 = fold1,
      r = r,
      grid = cv_grid,
      nfolds = cv_folds,
      seed = cv_seed + 1,
      method = initialization_method,
      rho = rho_admm,
      tol = admm_tol,
      maxit = admm_maxit,
      inner_tol = admm_inner_tol,
      inner_maxit = admm_inner_maxit,
      solver = cvxr_solver,
      solver_options = cvxr_options,
      verbose = verbose
    )
    
    eta_initial <- cv_initialization$eta_initial
    
  } else if (is.null(eta_initial)) {
    
    eta_initial <- c_initial *
      sqrt(log(p + m) / nrow(pooled1$X))
  }
  
  # If penalties are not tuned by CV, construct them from s_u and s_v.
  if (!tune_penalties) {
    
    # U-side
    if ((is.null(eta_transfer) || is.null(eta_target)) &&
        is.null(s_u)) {
      stop(
        "Supply s_u, or supply eta_transfer and eta_target, ",
        "or set tune_penalties = TRUE."
      )
    }
    
    if (!is.null(s_u)) {
      complexity_u <- r + log(exp(1) * p / s_u)
      
      if (is.null(eta_transfer)) {
        eta_transfer <- c_transfer *
          sqrt(complexity_u / nrow(pooled2$X))
      }
      
      if (is.null(eta_target)) {
        eta_target <- c_target *
          sqrt(complexity_u / nrow(target2$X))
      }
    }
    
    
    # V-side
    if ((is.null(eta_transfer_v) || is.null(eta_target_v)) &&
        is.null(s_v)) {
      stop(
        "Supply s_v, or supply eta_transfer_v and eta_target_v, ",
        "or set tune_penalties = TRUE."
      )
    }
    
    if (!is.null(s_v)) {
      complexity_v <- r + log(exp(1) * m / s_v)
      
      if (is.null(eta_transfer_v)) {
        eta_transfer_v <- c_transfer *
          sqrt(complexity_v / nrow(pooled2$Y))
      }
      
      if (is.null(eta_target_v)) {
        eta_target_v <- c_target *
          sqrt(complexity_v / nrow(target2$Y))
      }
    }
  }
  
  # eta_initial is already available at this point.
  if (length(eta_initial) != 1L ||
      !is.finite(eta_initial) ||
      eta_initial < 0) {
    stop("eta_initial must be finite and nonnegative.")
  }
  
  # Step 1: pooled convex initialization.
  M1 <- .moments(pooled1)
  initialization <- .scca_initialization(
    M1$Sx, M1$Sy, M1$Sxy,
    r = r,
    eta = eta_initial,
    method = initialization_method,
    rho = rho_admm,
    tol = admm_tol,
    maxit = admm_maxit,
    inner_tol = admm_inner_tol,
    inner_maxit = admm_inner_maxit,
    solver = cvxr_solver,
    solver_options = cvxr_options,
    verbose = verbose
  )
  if (identical(initialization$method, "ADMM") &&
      !isTRUE(initialization$converged)) {
    warning("ADMM initialization did not converge.")
  }
  initial_svd <- svd(initialization$F, nu = r, nv = r)
  if (sum(initial_svd$d > 1e-8 * max(initial_svd$d)) < r) {
    stop("The initializer has numerical rank below r; reduce eta_initial.")
  }
  U_initial <- initial_svd$u[, seq_len(r), drop = FALSE]
  V_initial <- initial_svd$v[, seq_len(r), drop = FALSE]

  # ---------------------------------------------------------
  # Convert the CV-selected multipliers into the penalties
  # appropriate for the actual fold-2 sample sizes.
  # ---------------------------------------------------------
  if (tune_penalties) {
    
    complexity_u <- r + log(p)
    complexity_v <- r + log(m)
    
    eta_transfer <-
      b_transfer *
      sqrt(
        complexity_u /
          nrow(pooled2$X)
      )
    
    eta_transfer_v <-
      b_transfer *
      sqrt(
        complexity_v /
          nrow(pooled2$Y)
      )
    
    eta_target <-
      b_target *
      sqrt(
        complexity_u /
          nrow(target2$X)
      )
    
    eta_target_v <-
      b_target *
      sqrt(
        complexity_v /
          nrow(target2$Y)
      )
  }
  penalty_values <- c(
    eta_transfer,
    eta_target,
    eta_transfer_v,
    eta_target_v
  )
  
  if (any(!is.finite(penalty_values)) ||
      any(penalty_values <= 0)) {
    stop("All transfer and target penalties must be positive and finite.")
  }
  
  
  # Steps 2--3: pooled and target-only grouped multi-response regression.
  # U-side
  transfer_u <- .glmnet_group_lasso(
    pooled2$X,
    pooled2$Y %*% V_initial,
    eta_transfer
  )
  
  target_u <- .glmnet_group_lasso(
    target2$X,
    target2$Y %*% V_initial,
    eta_target
  )
  
  U_transfer <- transfer_u$L
  U_target <- target_u$L
  
  # V-side
  transfer_v <- .glmnet_group_lasso(
    pooled2$Y,
    pooled2$X %*% U_initial,
    eta_transfer_v
  )
  
  target_v <- .glmnet_group_lasso(
    target2$Y,
    target2$X %*% U_initial,
    eta_target_v
  )
  
  V_transfer <- transfer_v$L
  V_target <- target_v$L
  
  # Data-driven working sparsity levels.
  s_u_used <- if (is.null(s_u)) {
    max(r, .row_support_size(U_target))
  } else {
    s_u
  }
  
  s_v_used <- if (is.null(s_v)) {
    max(r, .row_support_size(V_target))
  } else {
    s_v
  }
  
  complexity_aggregation_u <-
    r + log(exp(1) * p / s_u_used)
  
  complexity_aggregation_v <-
    r + log(exp(1) * m / s_v_used)
  
  
  if (is.null(aggregation_slack)) {
    aggregation_slack <- c_aggregation *
      s_u_used * complexity_aggregation_u /
      nrow(target3$X)
  }
  
  if (is.null(epsilon_target)) {
    epsilon_target <- c_normalization * sqrt(
      (
        s_u_used +
          log(exp(1) * p / s_u_used)
      ) /
        nrow(target3$X)
    )
  }
  
  
  if (is.null(aggregation_slack_v)) {
    aggregation_slack_v <- c_aggregation *
      s_v_used * complexity_aggregation_v /
      nrow(target3$Y)
  }
  
  if (is.null(epsilon_target_v)) {
    epsilon_target_v <- c_normalization * sqrt(
      (
        s_v_used +
          log(exp(1) * m / s_v_used)
      ) /
        nrow(target3$Y)
    )
  }
  
  downstream_tuning <- c(
    aggregation_slack,
    epsilon_target,
    aggregation_slack_v,
    epsilon_target_v
  )
  
  if (any(!is.finite(downstream_tuning)) ||
      any(downstream_tuning < 0)) {
    stop("Aggregation and normalization tuning values must be nonnegative and finite.")
  }
  
  
  # Step 4: aggregation on target fold 3.
  M03 <- .moments(target3)
  loss_transfer_u <- .target_loss(
    U_transfer,
    M03$Sx,
    M03$Sxy,
    V_initial
  )
  
  loss_target_u <- .target_loss(
    U_target,
    M03$Sx,
    M03$Sxy,
    V_initial
  )
  loss_transfer_v <- .target_loss(
    V_transfer,
    M03$Sy,
    t(M03$Sxy),
    U_initial
  )
  
  loss_target_v <- .target_loss(
    V_target,
    M03$Sy,
    t(M03$Sxy),
    U_initial
  )
  
  use_transfer_u <- switch(
    aggregation_mode,
    adaptive =
      loss_transfer_u <=
      loss_target_u + aggregation_slack,
    transfer = TRUE,
    target = FALSE
  )
  
  U_aggregated <-
    if (use_transfer_u) U_transfer else U_target
  
  use_transfer_v <- switch(
    aggregation_mode,
    adaptive =
      loss_transfer_v <=
      loss_target_v + aggregation_slack_v,
    transfer = TRUE,
    target = FALSE
  )
  
  V_aggregated <-
    if (use_transfer_v) V_transfer else V_target
  
  # Step 5: pooled normalization and thresholded target correction.
  MP3 <- .moments(pooled3)
  if (normalization_mode == "target") {
    
    # Direct normalization using the target sample covariance.
    B_u <- crossprod(
      U_aggregated,
      M03$Sx %*% U_aggregated
    )
    
    B_v <- crossprod(
      V_aggregated,
      M03$Sy %*% V_aggregated
    )
    
    U_final <- U_aggregated %*%
      .inv_root(
        B_u,
        ridge = normalization_ridge
      )
    
    V_final <- V_aggregated %*%
      .inv_root(
        B_v,
        ridge = normalization_ridge
      )
    
    norm_u <- list(
      final = U_final,
      B = B_u
    )
    
    norm_v <- list(
      final = V_final,
      B = B_v
    )
    
  } else {
    
    # Proposed truncated normalization.
    norm_u <- .normalize_side(
      U_aggregated,
      MP3$Sx,
      M03$Sx,
      epsilon_target,
      normalization_ridge
    )
    
    norm_v <- .normalize_side(
      V_aggregated,
      MP3$Sy,
      M03$Sy,
      epsilon_target_v,
      normalization_ridge
    )
    
    U_final <- norm_u$final
    V_final <- norm_v$final
  }
  
  list(
    U = U_final,
    V = V_final,
    
    working_sparsity = c(
      U = s_u_used,
      V = s_v_used
    ),
    
    U_aggregated = U_aggregated,
    V_aggregated = V_aggregated,
    
    U_transfer = U_transfer,
    V_transfer = V_transfer,
    
    U_target = U_target,
    V_target = V_target,
    
    U_initial = U_initial,
    V_initial = V_initial,
    
    F_initial = initialization$F,
    
    selected = c(
      U = if (use_transfer_u) "transfer" else "target",
      V = if (use_transfer_v) "transfer" else "target"
    ),
    
    target_losses = list(
      U = c(
        transfer = loss_transfer_u,
        target = loss_target_u
      ),
      V = c(
        transfer = loss_transfer_v,
        target = loss_target_v
      )
    ),
    tuning = list(
      eta_initial = eta_initial,
      
      eta_transfer_u = eta_transfer,
      eta_target_u = eta_target,
      aggregation_slack_u = aggregation_slack,
      epsilon_target_u = epsilon_target,
      
      eta_transfer_v = eta_transfer_v,
      eta_target_v = eta_target_v,
      aggregation_slack_v = aggregation_slack_v,
      epsilon_target_v = epsilon_target_v
    ),
    cv_initialization = cv_initialization,
    cv_tuning = cv_tuning,
    initialization = initialization,
    
    glmnet_fits = list(
      transfer_u = transfer_u$fit,
      target_u = target_u$fit,
      transfer_v = transfer_v$fit,
      target_v = target_v$fit
    ),
    
    folds = folds
  )
}
