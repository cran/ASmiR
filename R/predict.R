#' Locate bundled ASmiR model and feature-ranking files
#' @param stress One of cold, drought, heat, salt, or all.
#' @return Named character vector of file paths.
#' @export
asmir_model_files <- function(stress = c("cold", "drought", "heat", "salt", "all")) {
  stress <- match.arg(stress)
  ss <- if (stress == "all") c("cold", "drought", "heat", "salt") else stress
  out <- unlist(lapply(ss, function(s) c(
    model = system.file("extdata", "models", paste0("mirbase_mla_", s, ".rds"), package = "ASmiR"),
    feature_ranking = system.file("extdata", "models", paste0("mirbase_", s, "_svmrfe_feature.txt"), package = "ASmiR")
  )), use.names = FALSE)
  names(out) <- unlist(lapply(ss, function(s) paste(s, c("model", "feature_ranking"), sep = ".")))
  if (any(!nzchar(out)) || any(!file.exists(out))) stop("One or more bundled model files are missing. Reinstall ASmiR.", call. = FALSE)
  out
}

#' Predict stress association for plant miRNA sequences
#' @param x Input supported by asmir_read_fasta.
#' @param stress Stress model: cold, drought, heat, salt, or all.
#' @param threshold Probability threshold for stress-associated class Y.
#' @param model_file Optional override for a single-stress model RDS file.
#' @param feature_file Optional override for a single-stress ranking file.
#' @return Data frame with sequence name, stress, raw Y probability, confidence-like probability and label.
#' @export
asmir_predict <- function(x, stress = c("cold", "drought", "heat", "salt", "all"), threshold = 0.5,
                          model_file = NULL, feature_file = NULL) {
  stress <- match.arg(stress)
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) || threshold < 0 || threshold > 1) stop("threshold must be a number between 0 and 1.", call. = FALSE)
  if (!requireNamespace("e1071", quietly = TRUE)) stop("Package 'e1071' is required for SVM prediction.", call. = FALSE)
  d <- asmir_read_fasta(x)
  val <- asmir_validate_sequences(d)
  if (!all(val$valid)) stop(paste("Invalid RNA sequence(s):", paste(val$sequence_id[!val$valid], collapse = ", ")), call. = FALSE)
  raw <- asmir_features(d)
  stresses <- if (stress == "all") c("cold", "drought", "heat", "salt") else stress
  nfeat <- c(cold = 380L, drought = 272L, heat = 174L, salt = 340L)
  labels <- c(cold = "cold-stress", drought = "drought-stress", heat = "heat-stress", salt = "salt-stress")
  negative <- c(cold = "no-cold-stress", drought = "no-drought-stress", heat = "no-heat-stress", salt = "no-salt-stress")
  out <- vector("list", length(stresses))
  for (i in seq_along(stresses)) {
    s <- stresses[i]
    paths <- asmir_model_files(s)
    mf <- if (!is.null(model_file) && length(stresses) == 1L) model_file else unname(paths[paste0(s, ".model")])
    ff <- if (!is.null(feature_file) && length(stresses) == 1L) feature_file else unname(paths[paste0(s, ".feature_ranking")])
    idx <- base::scan(ff, what = integer(), quiet = TRUE)
    if (length(idx) < nfeat[[s]]) stop("Feature ranking for ", s, " contains fewer than ", nfeat[[s]], " indices.", call. = FALSE)
    idx <- idx[seq_len(nfeat[[s]])]
    if (anyNA(idx) || any(idx < 1L | idx > ncol(raw))) stop("Invalid feature indices in ranking file for ", s, ".", call. = FALSE)
    selected <- raw[, idx, drop = FALSE]
    colnames(selected) <- paste0("Species_", seq.int(2L, ncol(selected) + 1L))
    model <- readRDS(mf)
    # Use the exported predict generic; R dispatches to e1071's predict.svm method.
    pred <- tryCatch(
      stats::predict(
        model,
        newdata = as.data.frame(selected, check.names = FALSE),
        probability = TRUE
      ),
      error = function(e) stop(
        "SVM prediction failed for ", s, ": ", conditionMessage(e),
        call. = FALSE
      )
    )
    probs <- attr(pred, "probabilities")
    if (is.null(probs) || !"Y" %in% colnames(probs)) stop("SVM model for ", s, " did not return a probability column named Y.", call. = FALSE)
    p <- as.numeric(probs[, "Y"])
    label <- ifelse(p >= threshold, labels[[s]], negative[[s]])
    out[[i]] <- data.frame(Name = d$sequence_id, Stress = s, p_stress = p,
                           Probability = pmax(p, 1 - p), Label = label,
                           stringsAsFactors = FALSE)
  }
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}
