#install.packages("BiocManager")
# BiocManager::install("curatedTCGAData")
#BiocManager::install("TCGAutils")

library(curatedTCGAData)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(TCGAutils)


dir.create(
  "E:/fb/R-code/Trans-sparse CCA/TCGA/ExperimentHub",
  recursive = TRUE,
  showWarnings = FALSE
)

ExperimentHub::setExperimentHubOption(
  "CACHE",
  "E:/fb/R-code/Trans-sparse CCA/TCGA/ExperimentHub"
)

ExperimentHub::setExperimentHubOption(
  "ASK",
  FALSE
)


tcga <- curatedTCGAData::curatedTCGAData(
  diseaseCode = c(
    "LUAD",
    "LUSC",
    "HNSC",
    "BLCA",
    "ESCA"
  ),
  assays = c(
    "RNASeq2GeneNorm",
    "RPPAArray"
  ),
  version = "2.1.1",
  dry.run = FALSE
)

tcga <- readRDS(
  "E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_5cancers_RNA_RPPA.rds"
)

#domains_raw <- readRDS(
#  "E:/fb/R-code/Trans-sparse CCA/TCGA/domains_raw.rds"
#)

exps <- MultiAssayExperiment::experiments(tcga)

sm <- as.data.frame(
  MultiAssayExperiment::sampleMap(tcga)
)

sm$assay <- as.character(sm$assay)
sm$primary <- as.character(sm$primary)
sm$colname <- as.character(sm$colname)


cancers <- c(
  "ESCA",
  "LUAD",
  "LUSC",
  "HNSC",
  "BLCA"
)

rna_common <- Reduce(
  intersect,
  lapply(cancers, function(cc) {
    rownames(
      exps[[
        paste0(
          cc,
          "_RNASeq2GeneNorm-20160128"
        )
      ]]
    )
  })
)

rppa_common <- Reduce(
  intersect,
  lapply(cancers, function(cc) {
    rownames(
      exps[[
        paste0(
          cc,
          "_RPPAArray-20160128"
        )
      ]]
    )
  })
)

length(rna_common)
length(rppa_common)



extract_paired_domain <- function(
    cancer,
    rna_features = rna_common,
    rppa_features = rppa_common) {
  
  rna_name <- paste0(
    cancer,
    "_RNASeq2GeneNorm-20160128"
  )
  
  rppa_name <- paste0(
    cancer,
    "_RPPAArray-20160128"
  )
  
  rna_se <- exps[[rna_name]]
  rppa_se <- exps[[rppa_name]]
  
  rna_map <- sm[
    sm$assay == rna_name &
      substr(sm$colname, 14, 15) == "01",
    ,
    drop = FALSE
  ]
  
  rppa_map <- sm[
    sm$assay == rppa_name &
      substr(sm$colname, 14, 15) == "01",
    ,
    drop = FALSE
  ]
  
  patients <- sort(
    intersect(
      rna_map$primary,
      rppa_map$primary
    )
  )
  
  rna_cols <- rna_map$colname[
    match(
      patients,
      rna_map$primary
    )
  ]
  
  rppa_cols <- rppa_map$colname[
    match(
      patients,
      rppa_map$primary
    )
  ]
  
  stopifnot(
    all(rna_cols %in% colnames(rna_se)),
    all(rppa_cols %in% colnames(rppa_se))
  )
  
  X <- t(
    SummarizedExperiment::assay(rna_se)[
      rna_features,
      rna_cols,
      drop = FALSE
    ]
  )
  
  Y <- t(
    SummarizedExperiment::assay(rppa_se)[
      rppa_features,
      rppa_cols,
      drop = FALSE
    ]
  )
  
  rownames(X) <- patients
  rownames(Y) <- patients
  
  stopifnot(
    identical(rownames(X), rownames(Y))
  )
  
  list(
    X = X,
    Y = Y,
    patient = patients
  )
}


domains_raw <- setNames(
  lapply(
    cancers,
    extract_paired_domain
  ),
  cancers
)


na_table <- t(
  sapply(
    domains_raw,
    function(z) {
      c(
        RNA_NA = sum(is.na(z$X)),
        RNA_NaN = sum(is.nan(z$X)),
        RPPA_NA = sum(is.na(z$Y)),
        RPPA_NaN = sum(is.nan(z$Y))
      )
    }
  )
)

na_table
na_feature_table <- t(
  sapply(
    domains_raw,
    function(z) {
      c(
        RNA_genes_with_NA =
          sum(colSums(is.na(z$X)) > 0),
        
        RPPA_proteins_with_NA =
          sum(colSums(is.na(z$Y)) > 0)
      )
    }
  )
)
na_feature_table

saveRDS(
  domains_raw,
  "E:/fb/R-code/Trans-sparse CCA/TCGA/domains_raw.rds"
)

file.exists(
  "E:/fb/R-code/Trans-sparse CCA/TCGA/domains_raw.rds"
)


## ============================================================
## Preprocessing configuration
## ============================================================

#target_name <- "ESCA"
#target_name <- "LUAD"
#target_name <- "LUSC"
#target_name <- "HNSC"
target_name <- "BLCA"

source_names <- setdiff(
  cancers,
  target_name
)

impute_k <- 10L
screen_fraction <- 0.02


## ============================================================
## KNN imputation within each cancer
## ============================================================

if (!requireNamespace("impute", quietly = TRUE)) {
  stop(
    "Package 'impute' is required. ",
    "Install it once using BiocManager::install('impute')."
  )
}

impute_assay <- function(
    Z,
    k = 10L) {
  
  if (!anyNA(Z)) {
    return(Z)
  }
  
  Z_imp <- impute::impute.knn(
    t(Z),
    k = k
  )$data
  
  Z_imp <- t(Z_imp)
  
  dimnames(Z_imp) <- dimnames(Z)
  
  Z_imp
}

domains_imp <- lapply(
  domains_raw,
  function(z) {
    
    z$X <- impute_assay(
      z$X,
      k = 10L
    )
    
    z$Y <- impute_assay(
      z$Y,
      k = 10L
    )
    
    z
  }
)

stopifnot(
  all(
    sapply(
      domains_imp,
      function(z) {
        !anyNA(z$X) &&
          !anyNA(z$Y)
      }
    )
  )
)


## ============================================================
## Target-only RNA screening
## Keep top 2% most variable ESCA genes
## ============================================================

target_var <- apply(
  domains_imp[[target_name]]$X,
  2,
  stats::var
)

stopifnot(
  !anyNA(target_var),
  all(is.finite(target_var))
)


p_screen <- ceiling(
  screen_fraction *
    length(target_var)
)


selected_genes <- names(
  sort(
    target_var,
    decreasing = TRUE
  )
)[seq_len(p_screen)]


cat(
  "Number of common RNA genes before screening:",
  ncol(domains_imp[[target_name]]$X),
  "\n"
)

cat(
  "Number of selected RNA genes:",
  length(selected_genes),
  "\n"
)

## ============================================================
## Apply ESCA-selected genes to all domains
## ============================================================

domains_screened <- lapply(
  domains_imp,
  function(z) {
    
    list(
      X = z$X[
        ,
        selected_genes,
        drop = FALSE
      ],
      
      Y = z$Y,
      
      patient = z$patient
    )
  }
)

processed_dimensions <- sapply(
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

print(
  processed_dimensions
)

## ============================================================
## Save processed TCGA data
## ============================================================

tcga_processed <- list(
  
  domains = domains_screened,
  
  target_name = target_name,
  
  source_names = source_names,
  
  selected_genes = selected_genes,
  
  preprocessing = list(
    impute_method = "KNN",
    impute_k = impute_k,
    screen_method =
      "target variance",
    screen_fraction =
      screen_fraction,
    p_before =
      ncol(
        domains_imp[[target_name]]$X
      ),
    p_after =
      length(selected_genes),
    m =
      ncol(
        domains_imp[[target_name]]$Y
      )
  )
)


saveRDS(
  tcga_processed,
  #"E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed_ESCA.rds"
  #"E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed_LUAD.rds"
  #"E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed_LUSC.rds"
  #"E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed_HNSC.rds"
  "E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed_BLCA.rds"
)

success_by_rep <- reshape(
  real_results[
    ,
    c("rep", "method", "fit_failed")
  ],
  idvar = "rep",
  timevar = "method",
  direction = "wide"
)

common_reps <- success_by_rep$rep[
  rowSums(
    success_by_rep[
      ,
      grep("^fit_failed\\.", names(success_by_rep))
    ]
  ) == 0
]

common_results <- real_results[
  real_results$rep %in% common_reps,
  ,
  drop = FALSE
]

.summarize_real_results(
  common_results,
  r = 2L
)

#file.exists(
#  "E:/fb/R-code/Trans-sparse CCA/TCGA/tcga_processed.rds"
#)