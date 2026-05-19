## ===== Simple BLAST clustering (LibreOffice-friendly) =====
## Place files "R1.tabular" and "R2.tabular" in project folder (or edit FILES).

options(stringsAsFactors = FALSE)

# ---- Config ----
FILES <- c(R1 = "1a R1 vs Thermus thermophillus HB8.tabular", R2 = "1a R2 vs Thermus thermophillus HB8.tabular")  # change names if needed
PID_MIN   <- 98
EVAL_MAX  <- 1e-5
LEN_MIN   <- 50
BIN       <- 200          # cluster window size (bp)
SUPPORT_N <- 1            # keep clusters with > SUPPORT_N reads (i.e., ≥2)
DROP_TIES <- TRUE         # drop reads whose top hit is tied (ambiguous)

pp <- function(...) { cat(sprintf(...), "\n") }

read_blast <- function(path) {
  if (!file.exists(path)) stop("File not found: ", normalizePath(path))
  pp("Reading: %s", path)
  dat <- read.delim(path, header = FALSE, sep = "\t", quote = "", check.names = FALSE)
  if (ncol(dat) < 12L) stop("Expected ≥12 columns (BLAST outfmt 6). Got: ", ncol(dat))
  # keep needed cols: 1,2,3,4,7,8,9,10,11,12
  pick <- c(1,2,3,4,7,8,9,10,11,12)
  dat  <- dat[, pick, drop = FALSE]
  names(dat) <- c("query_id","sseqid","pident","alignment_length",
                  "qstart","qend","sstart","send","evalue","bitscore")
  # numeric coercion
  numcols <- c("pident","alignment_length","qstart","qend","sstart","send","evalue","bitscore")
  for (nm in numcols) dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
  dat <- dat[complete.cases(dat[, numcols]), ]
  pp("  Rows after coercion: %d", nrow(dat))
  dat
}

best_hit <- function(x, drop_ties = TRUE) {
  # order by: query_id, -bitscore, evalue, -alignment_length, -pident
  ord <- with(x, order(query_id, -bitscore, evalue, -alignment_length, -pident))
  x <- x[ord, ]
  first_idx <- which(!duplicated(x$query_id))
  if (!drop_ties) return(x[first_idx, , drop = FALSE])
  # Drop ties for top (bitscore,evalue)
  top_key  <- paste(x$query_id, x$bitscore, x$evalue, sep = "|")
  counts   <- table(top_key)
  firstkey <- top_key[first_idx]
  keep     <- counts[firstkey] == 1L
  keepset  <- x$query_id[first_idx[keep]]
  x[x$query_id %in% keepset & !duplicated(x$query_id), , drop = FALSE]
}

cluster_summ <- function(best, bin = 200, support_n = 1) {
  if (!nrow(best)) {
    return(data.frame(sseqid=character(), cluster_start=integer(), cluster_end_bin=integer(),
                      read_count=integer(), unique_reads=integer(), start_min=integer(),
                      end_max=integer(), mean_pident=numeric(), median_len=numeric()))
  }
  min_pos <- pmin(best$sstart, best$send)
  max_pos <- pmax(best$sstart, best$send)
  cstart  <- (min_pos %/% bin) * bin
  cend    <- cstart + (bin - 1)
  
  df <- data.frame(
    sseqid = best$sseqid,
    cluster_start = cstart,
    cluster_end_bin = cend,
    min_pos = min_pos,
    max_pos = max_pos,
    pident = best$pident,
    alignment_length = best$alignment_length,
    query_id = best$query_id,
    stringsAsFactors = FALSE
  )
  
  # aggregate per contig & bin
  key <- interaction(df$sseqid, df$cluster_start, drop = TRUE)
  idx_list <- split(seq_len(nrow(df)), key)
  
  out <- lapply(idx_list, function(ix) {
    g <- df[ix, , drop = FALSE]
    data.frame(
      sseqid        = g$sseqid[1],
      cluster_start = g$cluster_start[1],
      cluster_end_bin = g$cluster_end_bin[1],
      read_count    = nrow(g),
      unique_reads  = length(unique(g$query_id)),
      start_min     = min(g$min_pos),
      end_max       = max(g$max_pos),
      mean_pident   = mean(g$pident),
      median_len    = stats::median(g$alignment_length),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[out$read_count > support_n, , drop = FALSE]
  out[order(out$sseqid, out$cluster_start), , drop = FALSE]
}

process_one <- function(file, label) {
  pp("\n=== %s ===", label)
  x <- read_blast(file)
  
  # filters
  x <- x[x$pident >= PID_MIN & x$evalue <= EVAL_MAX & x$alignment_length >= LEN_MIN, ]
  pp("  After filters: %d rows", nrow(x))
  if (!nrow(x)) {
    out_empty <- paste0("clusters_", label, "_EMPTY.tsv")
    write.table(x, out_empty, sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
    pp("  Wrote empty TSV: %s", normalizePath(out_empty))
    return(invisible(NULL))
  }
  
  # best hit per read
  b <- best_hit(x, DROP_TIES)
  pp("  After best-hit%s: %d reads",
     if (DROP_TIES) " (drop ties)" else "", nrow(b))
  
  # cluster
  s <- cluster_summ(b, BIN, SUPPORT_N)
  pp("  Clusters kept (> %d reads): %d", SUPPORT_N, nrow(s))
  
  # ---- Write outputs (LibreOffice friendly) ----
  base <- paste0("clusters_", label)
  
  # 1) TSV (recommended)
  tsv_path <- paste0(base, ".tsv")
  write.table(s, tsv_path, sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  
  # 2) Semicolon CSV (EU locale-friendly)
  sc_path <- paste0(base, "_sc.csv")
  write.csv2(s, sc_path, row.names = FALSE, fileEncoding = "UTF-8")
  
  # 3) XLSX (optional)
  xlsx_path <- paste0(base, ".xlsx")
  if (requireNamespace("writexl", quietly = TRUE)) {
    sheets <- setNames(list(s), label)   # <-- fixed: dynamic sheet name without :=
    writexl::write_xlsx(sheets, xlsx_path)
    pp("  ✅ Wrote: %s, %s, %s",
       normalizePath(tsv_path), normalizePath(sc_path), normalizePath(xlsx_path))
  } else {
    pp("  ✅ Wrote: %s, %s  (xlsx skipped; install.packages('writexl') to enable)",
       normalizePath(tsv_path), normalizePath(sc_path))
  }
}

# ---- Run for R1 and R2 ----
pp("Working directory: %s", getwd())
pp("Looking for files: %s", paste(FILES, collapse = ", "))
for (nm in names(FILES)) process_one(FILES[[nm]], nm)
pp("\nDone.")

