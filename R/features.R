#' Check Python feature-generation runtime
#' @return Invisible logical TRUE if Python can be invoked.
#' @export
asmir_check_python <- function() {
  py <- Sys.getenv("ASMiR_PYTHON", unset = "python")
  out <- tryCatch(system2(py, "--version", stdout = TRUE, stderr = TRUE), error = function(e) e)
  if (inherits(out, "error") || !length(out) || any(grepl("not found|not recognized|No such file", out, ignore.case = TRUE))) {
    stop("Python 3 was not found. Install Python 3 and/or set ASMiR_PYTHON to its executable path.", call. = FALSE)
  }
  message(paste(out, collapse = " "))
  invisible(TRUE)
}

#' Generate original ASmiR PseKNC features
#' @param x Input supported by asmir_read_fasta.
#' @return Numeric matrix with 1,379 columns, in original k=1..5 block order.
#' @export
asmir_features <- function(x) {
  d <- asmir_read_fasta(x)
  val <- asmir_validate_sequences(d)
  if (!all(val$valid)) stop(paste("Invalid RNA sequence(s):", paste(val$sequence_id[!val$valid], collapse = ", ")), call. = FALSE)
  py <- Sys.getenv("ASMiR_PYTHON", unset = "python")
  asmir_check_python()
  pse <- system.file("python", "pseknc", "pse.py", package = "ASmiR")
  if (!nzchar(pse) || !file.exists(pse)) stop("Bundled PseKNC script is missing from the installed package.", call. = FALSE)
  td <- tempfile("asmir_"); dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  in_file <- file.path(td, "input_dna.fasta")
  seq_dna <- chartr("U", "T", d$sequence)
  lines <- unlist(Map(function(id, sq) c(paste0(">", gsub("[\r\n]", "_", id)), sq), d$sequence_id, seq_dna), use.names = FALSE)
  writeLines(lines, in_file, useBytes = TRUE)
  blocks <- vector("list", 5L)
  for (k in 1:5) {
    out_file <- file.path(td, paste0("pseknc_k", k, ".txt"))
    args <- c(shQuote(pse), shQuote(in_file), shQuote(out_file), "DNA", "PseKNC", "-lamada", "3", "-w", "0.2", "-k", as.character(k), "-f", "tab")
    out <- tryCatch(system2(py, args = args, stdout = TRUE, stderr = TRUE), error = function(e) e)
    if (inherits(out, "error") || !file.exists(out_file)) {
      stop("PseKNC feature generation failed for k=", k, ". Python output: ", paste(out, collapse = "\n"), call. = FALSE)
    }
    b <- tryCatch(as.matrix(utils::read.delim(out_file, header = FALSE, check.names = FALSE)), error = function(e) e)
    if (inherits(b, "error")) stop("Could not read PseKNC output for k=", k, ": ", conditionMessage(b), call. = FALSE)
    storage.mode(b) <- "numeric"
    if (nrow(b) != nrow(d)) stop("PseKNC returned ", nrow(b), " rows for k=", k, "; expected ", nrow(d), ".", call. = FALSE)
    blocks[[k]] <- b
  }
  feat <- do.call(cbind, blocks)
  if (ncol(feat) != 1379L) stop("Unexpected PseKNC feature count: ", ncol(feat), " (expected 1379).", call. = FALSE)
  rownames(feat) <- d$sequence_id
  colnames(feat) <- as.character(seq_len(ncol(feat)))
  feat
}
