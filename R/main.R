#' Encode DNA sequences by nucleotide position
#'
#' Split each input sequence into individual bases and store each
#' nucleotide position in a separate column.
#'
#' @param dna_strings A character vector of DNA sequences of equal length,
#'   using uppercase A, T, C, and G.
#' @return A data frame with one row per sequence. Columns are named
#'   nt_pos1, nt_pos2, and so on. Each column is a factor with levels
#'   A, T, C, and G.
#' @keywords internal
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}
#' Predict m6A status for multiple candidate sites
#'
#' Use a fitted random forest model to predict the probability and
#' classification of each candidate site. Each input row represents one site.
#'
#' @param ml_fit A fitted random forest classifier with a class named Positive.
#' @param feature_df A data frame containing gc_content, RNA_type, RNA_region,
#'   exon_length, distance_to_junction, evolutionary_conservation, and DNA_5mer.
#'   Use the same numeric scales as the training data. RNA_type must be one
#'   of "mRNA", "lincRNA", "lncRNA", or "pseudogene". RNA_region must be one
#'   of "CDS", "intron", "3'UTR", or "5'UTR". Each DNA_5mer must contain
#'   five uppercase bases drawn from A, T, C, and G.
#' @param positive_threshold Probability cutoff, default 0.5. A site is
#'   classified as Positive only when its probability exceeds the cutoff.
#'   A probability equal to the cutoff is classified as Negative.
#' @return The input data frame with two additional columns:
#'   predicted_m6A_prob (numeric probability of the Positive class) and
#'   predicted_m6A_status (Positive or Negative).
#' @seealso [prediction_single()] for predicting one site.
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' example_df <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv",
#'   package = "m6APrediction", mustWork = TRUE
#' ), stringsAsFactors = FALSE)
#' batch_result <- prediction_multiple(
#'   ml_fit = ml_fit,
#'   feature_df = example_df[1:3, , drop = FALSE],
#'   positive_threshold = 0.6
#' )
#' batch_result[, c("predicted_m6A_prob", "predicted_m6A_status")]
#' @importFrom stats predict
#' @importFrom randomForest randomForest
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}

#' Predict m6A status for a single candidate site
#'
#' Use a fitted random forest model to estimate the probability that one
#' candidate site is m6A-positive and classify it using a probability cutoff.
#'
#' @param ml_fit A fitted random forest classifier with a class named Positive.
#' @param gc_content Numeric GC content on the same scale as the training data.
#' @param RNA_type One of "mRNA", "lincRNA", "lncRNA", or "pseudogene".
#' @param RNA_region One of "CDS", "intron", "3'UTR", or "5'UTR".
#' @param exon_length Numeric exon length in the units used for training.
#' @param distance_to_junction Numeric distance to the junction, measured
#'   in the same units as the training data.
#' @param evolutionary_conservation Numeric conservation score for the site.
#' @param DNA_5mer Five-base DNA sequence in uppercase A, T, C, and G.
#' @param positive_threshold Probability cutoff, default 0.5. A probability
#'   strictly above the cutoff is Positive. A tie is Negative.
#' @return A named character vector with elements `predicted_m6A_prob`
#'   and `predicted_m6A_status`. Use `as.numeric()` to extract the
#'   probability as a number.
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' prediction_single(
#'   ml_fit = ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
