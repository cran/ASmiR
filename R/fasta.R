#' Read ASmiR FASTA input
#' @param x FASTA path, data.frame with sequence_id and sequence columns, or named/unnamed character vector.
#' @return A data.frame with sequence_id and sequence.
#' @export
asmir_read_fasta <- function(x) {
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    ss <- Biostrings::readRNAStringSet(x, format = "fasta")
    return(data.frame(sequence_id = names(ss), sequence = toupper(as.character(ss)),
                      stringsAsFactors = FALSE, check.names = FALSE))
  }
  if (is.data.frame(x)) {
    nms <- names(x)
    id_col <- intersect(c("sequence_id", "Name", "name", "ID", "id"), nms)[1]
    seq_col <- intersect(c("sequence", "Sequence", "seq", "RNA"), nms)[1]
    if (is.na(id_col) || is.na(seq_col)) stop("Data frame must contain sequence_id/name and sequence columns.", call. = FALSE)
    return(data.frame(sequence_id = as.character(x[[id_col]]), sequence = toupper(gsub("\\s+", "", as.character(x[[seq_col]]))), stringsAsFactors = FALSE))
  }
  if (is.character(x)) {
    ids <- names(x)
    if (is.null(ids)) ids <- paste0("sequence_", seq_along(x))
    return(data.frame(sequence_id = ids, sequence = toupper(gsub("\\s+", "", unname(x))), stringsAsFactors = FALSE))
  }
  stop("x must be a FASTA path, a compatible data.frame, or a character vector.", call. = FALSE)
}

#' Validate ASmiR RNA sequences
#' @param x Input supported by asmir_read_fasta.
#' @return Validation data.frame with sequence ID, length, valid flag and message.
#' @export
asmir_validate_sequences <- function(x) {
  d <- asmir_read_fasta(x)
  valid <- !is.na(d$sequence) & nzchar(d$sequence) & grepl("^[AUGC]+$", d$sequence)
  msg <- ifelse(is.na(d$sequence) | !nzchar(d$sequence), "Empty sequence",
                ifelse(valid, "OK", "Non-standard residues: only A, U, G and C are accepted"))
  data.frame(sequence_id = d$sequence_id, length = nchar(d$sequence), valid = valid,
             message = msg, stringsAsFactors = FALSE)
}
