#!/usr/bin/env Rscript
# DADA2 amplicon denoising pipeline.
# Local variables use descriptive snake_case names. Scientific operations,
# parameter values, output schemas and filenames preserve the archived script.
# Arguments: input_directory output_directory maxEE_f maxEE_r truncLen_f
#            truncLen_r minOverlap threads
# Runtime recorded in the archive: R 4.4.3, dada2 1.34.0, digest 0.6.37.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2) {
  stop("Usage: Rscript dada2_denoise.R <trimmed_dir> <output_dir> [maxEE_f] [maxEE_r] [truncLen_f] [truncLen_r] [minOverlap] [threads]")
}

input_directory <- arguments[1]
output_directory <- arguments[2]
max_expected_errors_forward <- ifelse(length(arguments) >= 3, as.numeric(arguments[3]), 2.0)
max_expected_errors_reverse <- ifelse(length(arguments) >= 4, as.numeric(arguments[4]), 2.0)
truncation_length_forward <- ifelse(length(arguments) >= 5, as.integer(arguments[5]), 0)
truncation_length_reverse <- ifelse(length(arguments) >= 6, as.integer(arguments[6]), 0)
minimum_overlap <- ifelse(length(arguments) >= 7, as.integer(arguments[7]), 12)
threads <- ifelse(length(arguments) >= 8, as.integer(arguments[8]), 4)

library(dada2)

dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

cat("=== DADA2 Denoise ===\n")
cat(sprintf("Input: %s\n", input_directory))
cat(sprintf("Output: %s\n", output_directory))

forward_read_files <- sort(list.files(
  input_directory,
  pattern = "_R1\\.fastq\\.gz$|_R1\\.fq\\.gz$|_1\\.fastq\\.gz$|_1\\.fq\\.gz$",
  full.names = TRUE
))
reverse_read_files <- sort(list.files(
  input_directory,
  pattern = "_R2\\.fastq\\.gz$|_R2\\.fq\\.gz$|_2\\.fastq\\.gz$|_2\\.fq\\.gz$",
  full.names = TRUE
))

if (length(forward_read_files) == 0) {
  forward_read_files <- sort(list.files(input_directory, pattern = "trimmed_R1", full.names = TRUE))
  reverse_read_files <- sort(list.files(input_directory, pattern = "trimmed_R2", full.names = TRUE))
}

if (length(forward_read_files) == 0 || length(reverse_read_files) == 0) {
  stop("No paired fastq files found in trimmed_dir")
}

cat(sprintf("Found %d paired samples\n", length(forward_read_files)))

sample_ids <- sub("_trimmed_R[12]\\.fastq\\.gz$", "", basename(forward_read_files))
sample_ids <- sub("_R[12]\\.fastq\\.gz$", "", sample_ids)
sample_ids <- sub("_[12]\\.fastq\\.gz$", "", sample_ids)
sample_ids <- sub("_[12]\\.fq\\.gz$", "", sample_ids)

filtered_forward_files <- file.path(
  output_directory, "filtered", paste0(sample_ids, "_F_filt.fastq.gz")
)
filtered_reverse_files <- file.path(
  output_directory, "filtered", paste0(sample_ids, "_R_filt.fastq.gz")
)
names(filtered_forward_files) <- sample_ids
names(filtered_reverse_files) <- sample_ids

cat("\n--- filterAndTrim ---\n")
forward_truncation <- if (truncation_length_forward > 0) truncation_length_forward else 0
reverse_truncation <- if (truncation_length_reverse > 0) truncation_length_reverse else 0

filter_statistics <- filterAndTrim(
  forward_read_files, filtered_forward_files, reverse_read_files, filtered_reverse_files,
  maxEE = c(max_expected_errors_forward, max_expected_errors_reverse),
  truncLen = c(forward_truncation, reverse_truncation),
  truncQ = 2,
  maxN = 0,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = threads
)
cat("Filter stats:\n")
print(filter_statistics)

write.csv(filter_statistics, file.path(output_directory, "filter_stats.csv"))

forward_files_exist <- file.exists(filtered_forward_files)
reverse_files_exist <- file.exists(filtered_reverse_files)
retained_sample_mask <- forward_files_exist & reverse_files_exist
if (sum(retained_sample_mask) == 0) {
  stop("No reads passed filtering")
}

filtered_forward_files <- filtered_forward_files[retained_sample_mask]
filtered_reverse_files <- filtered_reverse_files[retained_sample_mask]
sample_ids <- sample_ids[retained_sample_mask]

cat("\n--- learnErrors ---\n")
forward_error_model <- learnErrors(filtered_forward_files, multithread = threads)
reverse_error_model <- learnErrors(filtered_reverse_files, multithread = threads)

cat("\n--- dada ---\n")
denoised_forward_reads <- dada(
  filtered_forward_files, err = forward_error_model, multithread = threads
)
denoised_reverse_reads <- dada(
  filtered_reverse_files, err = reverse_error_model, multithread = threads
)

cat("\n--- mergePairs ---\n")
merged_read_pairs <- mergePairs(
  denoised_forward_reads, filtered_forward_files,
  denoised_reverse_reads, filtered_reverse_files,
  minOverlap = minimum_overlap
)

cat("\n--- makeSequenceTable ---\n")
sequence_table <- makeSequenceTable(merged_read_pairs)
cat(sprintf("ASV table: %d samples x %d ASVs\n", nrow(sequence_table), ncol(sequence_table)))

cat("\n--- removeBimeraDenovo ---\n")
nonchimeric_sequence_table <- removeBimeraDenovo(
  sequence_table, method = "consensus", multithread = threads
)
cat(sprintf("After chimera removal: %d ASVs\n", ncol(nonchimeric_sequence_table)))
cat(sprintf("Retained %.1f%% of reads\n", sum(nonchimeric_sequence_table) / sum(sequence_table) * 100))

sequences <- colnames(nonchimeric_sequence_table)
asv_ids <- sapply(
  sequences, function(sequence) digest::digest(sequence, algo = "md5", serialize = FALSE)
)
colnames(nonchimeric_sequence_table) <- asv_ids

asv_table <- as.data.frame(nonchimeric_sequence_table)
asv_table$SampleID <- rownames(asv_table)
asv_table <- asv_table[, c(ncol(asv_table), 1:(ncol(asv_table) - 1))]
write.table(
  asv_table, file.path(output_directory, "asv_table.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

sequence_metadata <- data.frame(ASV_ID = asv_ids, Sequence = sequences, stringsAsFactors = FALSE)
fasta_file <- file.path(output_directory, "rep_seqs.fasta")
fasta_connection <- file(fasta_file, "w")
for (sequence_index in 1:nrow(sequence_metadata)) {
  writeLines(paste0(">", sequence_metadata$ASV_ID[sequence_index]), fasta_connection)
  writeLines(sequence_metadata$Sequence[sequence_index], fasta_connection)
}
close(fasta_connection)

write.table(
  sequence_metadata, file.path(output_directory, "asv_sequences.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

if (is(denoised_forward_reads, "dada")) denoised_forward_reads <- list(denoised_forward_reads)
if (is(denoised_reverse_reads, "dada")) denoised_reverse_reads <- list(denoised_reverse_reads)
if (is.data.frame(merged_read_pairs)) merged_read_pairs <- list(merged_read_pairs)

count_unique_reads <- function(read_result) sum(getUniques(read_result))

read_tracking_table <- data.frame(
  sample = sample_ids,
  input = filter_statistics[retained_sample_mask, 1],
  filtered = filter_statistics[retained_sample_mask, 2],
  denoisedF = sapply(denoised_forward_reads, count_unique_reads),
  denoisedR = sapply(denoised_reverse_reads, count_unique_reads),
  merged = sapply(merged_read_pairs, function(read_result) sum(read_result$abundance)),
  nonchim = rowSums(nonchimeric_sequence_table)
)
write.table(
  read_tracking_table, file.path(output_directory, "track_stats.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

cat("\n=== DADA2 Complete ===\n")
cat(sprintf("Output files in: %s\n", output_directory))
cat(sprintf("  asv_table.tsv: %d samples x %d ASVs\n", nrow(nonchimeric_sequence_table), ncol(nonchimeric_sequence_table)))
