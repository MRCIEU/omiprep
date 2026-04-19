## ── QC Pipeline Timing Benchmark ─────────────────────────────────────────────
## Runs quality_control() across sample sizes × platform-typical feature counts,
## records wall-clock time, prints summary tables, and writes a line plot PNG.
##
## Usage: Rscript scripts/qc_timing_benchmark.R
##        or source() from an interactive session with metaboprep loaded.
## ─────────────────────────────────────────────────────────────────────────────
library(devtools)
load_all()
#library(metaboprep)
library(ggplot2)

## ── 1. Benchmark grid ────────────────────────────────────────────────────────

n_samples_vec <- c(100, 1000, 10000)

platforms <- data.frame(
  platform   = factor(c("Nightingale", "Metabolon (subset)", "Olink", "Metabolon (full)", "SomaLogic"), 
                      levels = c("Nightingale", "Metabolon (subset)", "Olink", "Metabolon (full)",  "SomaLogic")),
  n_features = c(250,       1250,    3000,  5000,     11000),
  stringsAsFactors = FALSE
)

## ── 2. Helper: build a Metaboprep object of given dimensions ─────────────────

make_metaboprep <- function(n_samples, n_features, seed = 42) {
  set.seed(seed)
  
  sample_ids  <- paste0("s_", seq_len(n_samples))
  feature_ids <- paste0("f_", seq_len(n_features))
  
  # ── Correlated feature structure ──────────────────────────────────────────
  # ~30 % of features belong to correlated clusters (sizes 3–8), mimicking
  # pathway or protein-complex co-regulation. A shared latent factor drives
  # each cluster (loading 0.6–0.9); remaining features are independent noise.
  n_clustered  <- floor(n_features * 0.30)
  cluster_sizes <- sample(3:8, ceiling(n_clustered / 5), replace = TRUE)
  cluster_sizes <- cluster_sizes[cumsum(cluster_sizes) <= n_clustered]
  n_factors     <- length(cluster_sizes)
  
  # latent factors (one per cluster)
  latent <- matrix(rnorm(n_samples * n_factors),
                   nrow = n_samples, ncol = n_factors)
  
  mat <- matrix(rnorm(n_samples * n_features),
                nrow = n_samples, ncol = n_features)
  
  col <- 1L
  for (k in seq_len(n_factors)) {
    loading <- runif(cluster_sizes[k], 0.6, 0.9)
    for (j in seq_len(cluster_sizes[k])) {
      mat[, col] <- loading[j] * latent[, k] + sqrt(1 - loading[j]^2) * mat[, col]
      col <- col + 1L
    }
  }
  # remaining features are already independent standard-normal noise
  
  dimnames(mat) <- list(sample_ids, feature_ids)
  
  # ~15 % missingness applied after generating the signal
  n_miss <- floor(0.15 * n_samples * n_features)
  mat[sample(length(mat), n_miss)] <- NA
  
  samples <- data.frame(
    sample_id = sample_ids,
    age       = sample(18:70, n_samples, replace = TRUE),
    sex       = sample(c("male", "female"), n_samples, replace = TRUE),
    pos       = sample(c("batch1", "batch2"), n_samples, replace = TRUE),
    neg       = sample(c("batch1", "batch2"), n_samples, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  features <- data.frame(
    feature_id         = feature_ids,
    platform           = sample(c("pos", "neg"), n_features, replace = TRUE,
                                prob = c(0.20, 0.80)),
    pathway            = NA_character_,
    derived_feature    = sample(c(TRUE, FALSE), n_features, replace = TRUE,
                                prob = c(0.05, 0.95)),
    xenobiotic_feature = sample(c(TRUE, FALSE), n_features, replace = TRUE,
                                prob = c(0.10, 0.90)),
    stringsAsFactors = FALSE
  )
  
  Metaboprep(data = mat, samples = samples, features = features)
}

## ── 3. Run benchmark ─────────────────────────────────────────────────────────

grid <- merge(
  data.frame(n_samples = n_samples_vec),
  platforms,
  by = NULL  # cross join
)
grid <- grid[order(grid$platform, grid$n_samples), ]

results <- vector("list", nrow(grid))

for (i in seq_len(nrow(grid))) {
  ns   <- grid$n_samples[i]
  nf   <- grid$n_features[i]
  plat <- grid$platform[i]
  
  cat(sprintf("[%d/%d]  %-12s  samples = %6s  |  features = %s ",
              i, nrow(grid), plat,
              format(ns, big.mark = ","),
              format(nf, big.mark = ",")))
  
  obj <- make_metaboprep(ns, nf)
  
  result  <- NULL
  elapsed <- NA_real_
  result_mb <- NA_real_
  
  tryCatch({
    t <- system.time(
      suppressMessages(
        result <- quality_control(
          obj,
          source_layer        = "input",
          sample_missingness  = 0.2,
          feature_missingness = 0.2,
          total_peak_area_sd  = 5,
          outlier_udist       = 5,
          outlier_treatment   = "leave_be",
          winsorize_quantile  = 1.0,
          tree_cut_height     = 0.5,
          pc_outlier_sd       = 5,
          sample_ids          = NULL,
          feature_ids         = NULL, 
          cores               = 10, 
          fast                = TRUE
        )
      )
    )
    elapsed   <- unname(t["elapsed"])
    result_mb <- as.numeric(object.size(result)) / 1024^2
  }, error = function(e) {
    cat(sprintf("  ! ERROR: %s\n", conditionMessage(e)))
  })
  
  cat(sprintf("  → %.1f s  |  %.1f MB\n",
              ifelse(is.na(elapsed), 0, elapsed),
              ifelse(is.na(result_mb), 0, result_mb)))
  
  results[[i]] <- data.frame(
    platform   = plat,
    n_samples  = ns,
    n_features = nf,
    seconds    = elapsed,
    minutes    = elapsed / 60,
    result_mb  = result_mb
  )
}

timing <- do.call(rbind, results)

## ── 4. Print tables ───────────────────────────────────────────────────────────

cat("\n── Timing results (minutes) ──────────────────────────────────────────────\n")
timing_wide <- reshape(
  timing[, c("platform", "n_samples", "minutes")],
  idvar     = "n_samples",
  timevar   = "platform",
  direction = "wide"
)
colnames(timing_wide) <- sub("minutes\\.", "", colnames(timing_wide))
print(timing_wide, row.names = FALSE)

cat("\n── Output object size (MB) ───────────────────────────────────────────────\n")
mem_wide <- reshape(
  timing[, c("platform", "n_samples", "result_mb")],
  idvar     = "n_samples",
  timevar   = "platform",
  direction = "wide"
)
colnames(mem_wide) <- sub("result_mb\\.", "", colnames(mem_wide))
print(mem_wide, row.names = FALSE)

out_csv <- file.path("scripts", "qc_timing_results.csv")
write.csv(timing, out_csv, row.names = FALSE)
cat(sprintf("\nTable written to: %s\n", out_csv))

## ── 5. Line plot ──────────────────────────────────────────────────────────────

platform_colours <- c(
  "Nightingale"        = "#2196F3",
  "Metabolon (subset)" = "#4CAF50",
  "Metabolon (full)"   = "darkgreen",
  "Olink"              = "#FF9800",
  "SomaLogic"          = "#E91E63"
)


timing$platform <- factor(timing$platform, levels = platforms$platform)

# label for legend: platform + feature count
timing$platform_label <- factor(
  sprintf("%s (~%s features)", timing$platform,
          format(timing$n_features, big.mark = ",")),
  levels = sprintf("%s (~%s features)", platforms$platform,
                   format(platforms$n_features, big.mark = ","))
)

p <- ggplot(timing[!is.na(timing$minutes), ],
            aes(x = n_samples, y = minutes,
                colour = platform_label, group = platform_label)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5) +
  scale_x_log10(
    breaks = n_samples_vec,
    labels = function(x) format(x, big.mark = ",", scientific = FALSE)
  ) +
  scale_y_log10(
    breaks = c(0, 1, 10, 100, 1000)
  ) +
  scale_colour_manual(
    values = setNames(
      platform_colours[platforms$platform],
      levels(timing$platform_label)
    ),
    name = "Platform"
  ) +
  labs(
    title    = "metaboprep QC pipeline — runtime by platform",
    subtitle = "Wall-clock time for quality_control() at platform-typical feature counts",
    x        = "Number of samples (log scale)",
    y        = "Time (minutes, log scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(colour = "grey40", size = 10,
                                    margin = margin(b = 10)),
    axis.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.title     = element_text(face = "bold"),
    legend.position  = "right"
  )
p

out_png <- file.path("scripts", "qc_timing_by_platform.png")
ggsave(out_png, plot = p, width = 8, height = 5.5, dpi = 150)
cat(sprintf("Plot written to: %s\n", out_png))

## ── 6. metaboprep1 benchmark ──────────────────────────────────────────────────
## Runs the legacy metaboprep1 pipeline for both Nightingale (~250 features)
## and Metabolon (~1,000 features) at sample sizes n_samples_vec[-last].
## The GitHub repo is cloned once and reused for every run.
##
## Notes:
##   • 100,000 samples is omitted; the legacy pipeline would be impractically
##     slow and writing a 100k-sample Excel file is itself a bottleneck.
##   • Each run spawns a fresh Rscript subprocess, so timings include ~5–10 s
##     of R startup overhead.
##   • Metabolon: synthetic v1.1 Excel layout (single-sheet, feature rows).
##   • Nightingale: synthetic single-sheet "Worksheet" layout with sampleid /
##     success % anchors.

mp1_n_samples <- c(100, 1000)

mp1_platforms <- data.frame(
  platform              = c("Nightingale", "Olink"),
  n_features            = c(250,           3000),
  Nightingale_OR_Metabolon = c("Nightingale", "Metabolon"),
  feat_anno_run_mode_col   = c(NA_character_, "platform"),
  stringsAsFactors = FALSE
)

mp1_grid <- merge(data.frame(n_samples = mp1_n_samples), mp1_platforms, by = NULL)
mp1_grid <- mp1_grid[order(mp1_grid$platform, mp1_grid$n_samples), ]

## ── 6a. Helpers ───────────────────────────────────────────────────────────────

# Write synthetic data as a Metabolon v1.1-compatible Excel workbook.
# Layout mirrors what read_metabolon() expects:
#   rows 1:(data_header_row-1)  — sample metadata (col 6 = variable name)
#   row   data_header_row       — feature metadata column headers
#   rows (data_header_row+1):nrow — feature rows; cols > batch_header_col = data
write_metabolon_v1.1 <- function(n_samples, n_features, filepath, seed = 42) {
  set.seed(seed)
  
  sample_ids  <- paste0("S", seq_len(n_samples))
  feature_ids <- paste0("COMP_", seq_len(n_features))
  
  # numeric data (features × samples), log-normal-ish positive values
  mat <- matrix(
    round(abs(rnorm(n_features * n_samples, mean = 5000, sd = 2000)), 2),
    nrow = n_features, ncol = n_samples
  )
  n_miss <- floor(0.15 * n_features * n_samples)
  mat[sample(length(mat), n_miss)] <- NA
  
  # full worksheet matrix (character): (7 + n_features) × (6 + n_samples)
  # rows 1-6 = sample metadata, row 7 = feature header, rows 8+ = feature data
  n_rows <- 7 + n_features
  n_cols <- 6 + n_samples
  ws     <- matrix(NA_character_, nrow = n_rows, ncol = n_cols)
  
  # ── sample metadata (rows 1-6) ──────────────────────────────────────────────
  ws[1, 6] <- "SAMPLE NAME"
  ws[2, 6] <- "NEG"
  ws[3, 6] <- "POS"
  ws[4, 6] <- "age"
  ws[5, 6] <- "sex"
  ws[6, 6] <- "BOX_ID"
  
  ws[1, 7:(6 + n_samples)] <- sample_ids
  ws[2, 7:(6 + n_samples)] <- sample(c("batch1", "batch2"), n_samples, replace = TRUE)
  ws[3, 7:(6 + n_samples)] <- sample(c("batch1", "batch2"), n_samples, replace = TRUE)
  ws[4, 7:(6 + n_samples)] <- as.character(sample(18:70, n_samples, replace = TRUE))
  ws[5, 7:(6 + n_samples)] <- sample(c("M", "F"), n_samples, replace = TRUE)
  ws[6, 7:(6 + n_samples)] <- "box1"
  
  # ── feature metadata header (row 7) ─────────────────────────────────────────
  ws[7, 1:6] <- c("metabolite_id", "COMP_ID", "platform", "pathway", "kegg", "Group HMDB")
  
  # ── feature metadata values (rows 8+) ───────────────────────────────────────
  feat_rows <- 8:(7 + n_features)
  ws[feat_rows, 1] <- paste0("met", seq_len(n_features))
  ws[feat_rows, 2] <- feature_ids
  ws[feat_rows, 3] <- sample(c("pos", "neg"), n_features, replace = TRUE, prob = c(0.3, 0.7))
  
  # ── numeric data (rows 8+, cols 7+) ─────────────────────────────────────────
  ws[feat_rows, 7:(6 + n_samples)] <- ifelse(is.na(mat), NA_character_, as.character(mat))
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Explanation")
  openxlsx::addWorksheet(wb, "OrigScale")
  openxlsx::addWorksheet(wb, "ScaledImp")
  df <- as.data.frame(ws, stringsAsFactors = FALSE)
  openxlsx::writeData(wb, "OrigScale", df, colNames = FALSE, rowNames = FALSE)
  openxlsx::writeData(wb, "ScaledImp", df, colNames = FALSE, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, filepath, overwrite = TRUE)
  invisible(filepath)
}


# Write synthetic data as a Nightingale single-sheet ("Worksheet") Excel workbook.
# Layout mirrors what read_nightingale() expects for the single-sheet format:
#   row (head_row-1)  — sample metadata column names in cols 1:(sampleid_col-1)
#   row  head_row     — "sampleid" at sampleid_col; feature IDs in cols to the right
#   row  head_row+1   — units row (ignored by reader)
#   row  head_row+2   — "success %" at sampleid_col; success values to the right
#   rows head_row+3+  — sample data: flag cols, sample_id, numeric data
write_nightingale_worksheet <- function(n_samples, n_features, filepath, seed = 42) {
  set.seed(seed)
  
  sample_ids  <- paste0("S", seq_len(n_samples))
  feature_ids <- paste0("NMR_", seq_len(n_features))
  
  mat <- matrix(
    round(abs(rnorm(n_features * n_samples, mean = 1, sd = 0.3)), 4),
    nrow = n_samples, ncol = n_features
  )
  n_miss <- floor(0.15 * n_features * n_samples)
  mat[sample(length(mat), n_miss)] <- NA
  
  # worksheet dimensions: 12 header rows + n_samples data rows; 3 + n_features cols
  # col layout: col1="QC1", col2="QC2", col3="sampleid", cols4+= features
  n_header <- 12L
  n_rows   <- n_header + n_samples
  n_cols   <- 3L + n_features
  samp_col <- 3L
  feat_col <- 4L
  
  ws <- matrix(NA_character_, nrow = n_rows, ncol = n_cols)
  
  # row (head_row - 1) = row 9: sample-flag column names
  ws[9L, 1L] <- "QC1"
  ws[9L, 2L] <- "QC2"
  
  # row head_row = 10: sampleid anchor + feature IDs
  ws[10L, samp_col]             <- "sampleid"
  ws[10L, feat_col:(3+n_features)] <- feature_ids
  
  # row 11: units (ignored)
  ws[11L, feat_col:(3+n_features)] <- "mmol/L"
  
  # row 12: success % anchor + all-pass values
  ws[12L, samp_col]             <- "success %"
  ws[12L, feat_col:(3+n_features)] <- "1"
  
  # data rows 13+: cols 1-2 NA (flags=0), col 3 = sample_id, cols 4+ = data
  data_rows <- (n_header + 1L):n_rows
  ws[data_rows, samp_col] <- sample_ids
  ws[data_rows, feat_col:(3+n_features)] <-
    ifelse(is.na(mat), NA_character_, as.character(mat))
  
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Worksheet")
  df <- as.data.frame(ws, stringsAsFactors = FALSE)
  openxlsx::writeData(wb, "Worksheet", df, colNames = FALSE, rowNames = FALSE)
  openxlsx::saveWorkbook(wb, filepath, overwrite = TRUE)
  invisible(filepath)
}


# Parse key=value pairs from a metaboprep1 parameter_file.sh template.
parse_params <- function(lines) {
  params <- list()
  for (line in lines) {
    line_trim <- trimws(line)
    if (nchar(line_trim) == 0 || startsWith(line_trim, "#")) next
    if (grepl("=", line_trim, fixed = TRUE)) {
      parts         <- strsplit(line_trim, "=", fixed = TRUE)[[1]]
      params[[trimws(parts[1])]] <- trimws(parts[2])
    }
  }
  params
}

# Update parameter lines and write to outfile.
write_params <- function(lines, params, outfile) {
  out_lines <- lines
  for (i in seq_along(lines)) {
    line_trim <- trimws(lines[i])
    if (nchar(line_trim) == 0 || startsWith(line_trim, "#")) next
    if (grepl("=", line_trim, fixed = TRUE)) {
      key <- trimws(strsplit(line_trim, "=", fixed = TRUE)[[1]][1])
      if (key %in% names(params)) {
        out_lines[i] <- paste0(key, "=", params[[key]])
      }
    }
  }
  writeLines(out_lines, outfile)
}

## ── 6b. Clone metaboprep1 once ───────────────────────────────────────────────

mp1_tmp  <- tempfile("mp1_lib_")
mp1_ok   <- FALSE

cat("\nCloning metaboprep1 repository (once) ...\n")
clone_status <- system2("git",
                        c("clone", "https://github.com/MRCIEU/metaboprep.git", mp1_tmp),
                        stdout = FALSE, stderr = FALSE)

if (clone_status == 0) {
  system2("git", c("-C", mp1_tmp, "checkout", "bbe1f85"),
          stdout = FALSE, stderr = FALSE)
  
  # Inject load_all() + truncate at report step
  pipeline_script <- file.path(mp1_tmp, "run_metaboprep_pipeline.R")
  script_lines    <- readLines(pipeline_script)
  cut_index       <- grep("\\(X\\) Make Report", script_lines, fixed = FALSE)[1]
  script_lines    <- script_lines[seq_len(cut_index)]
  inject_lines    <- c("library(devtools)", sprintf('load_all("%s")', mp1_tmp))
  writeLines(c(inject_lines, "", script_lines), pipeline_script)
  
  # Read parameter template
  param_template <- readLines(file.path(mp1_tmp, "parameter_file.sh"))
  
  mp1_ok <- TRUE
  cat("Clone complete.\n\n")
} else {
  cat("! Clone failed — skipping metaboprep1 benchmark.\n\n")
}

## ── 6c. Run benchmark ────────────────────────────────────────────────────────

mp1_results <- vector("list", nrow(mp1_grid))

if (mp1_ok) {
  for (i in seq_len(nrow(mp1_grid))) {
    ns   <- mp1_grid$n_samples[i]
    nf   <- mp1_grid$n_features[i]
    plat <- mp1_grid$platform[i]
    nom  <- mp1_grid$Nightingale_OR_Metabolon[i]
    ramc <- mp1_grid$feat_anno_run_mode_col[i]
    
    cat(sprintf("[%d/%d]  metaboprep1  %-12s  samples = %6s  |  features = %s ",
                i, nrow(mp1_grid), plat,
                format(ns, big.mark = ","),
                format(nf, big.mark = ",")))
    
    elapsed <- NA_real_
    
    tryCatch({
      data_dir  <- tempfile("mp1data_")
      dir.create(data_dir)
      
      if (nom == "Metabolon") {
        xlsx_path <- file.path(data_dir, "synthetic_metabolon.xlsx")
        write_metabolon_v1.1(ns, nf, xlsx_path)
      } else {
        xlsx_path <- file.path(data_dir, "synthetic_nightingale.xlsx")
        write_nightingale_worksheet(ns, nf, xlsx_path)
      }
      
      params <- parse_params(param_template)
      params[["projectname"]]              <- "BENCH"
      params[["datadirectory"]]            <- data_dir
      params[["metabolite_data_file"]]     <- basename(xlsx_path)
      params[["Nightingale_OR_Metabolon"]] <- nom
      params[["feature_missingness"]]      <- "0.2"
      params[["sample_missingness"]]       <- "0.2"
      params[["total_peak_area_SD"]]       <- "5"
      params[["outlier_udist"]]            <- "5"
      params[["outlier_treatment"]]        <- "leave_be"
      params[["tree_cut_height"]]          <- "0.5"
      params[["PC_outlier_SD"]]            <- "5"
      if (!is.na(ramc)) params[["feat_anno_run_mode_col"]] <- ramc
      
      param_out <- file.path(data_dir, "parameter_file.sh")
      write_params(param_template, params, param_out)
      
      t <- system.time(
        system2(file.path(R.home("bin"), "Rscript"),
                c(file.path(mp1_tmp, "run_metaboprep_pipeline.R"), param_out),
                stdout = TRUE, stderr = TRUE)
      )
      elapsed <- unname(t["elapsed"])
      unlink(data_dir, recursive = TRUE)
      
    }, error = function(e) {
      cat(sprintf("  ! ERROR: %s\n", conditionMessage(e)))
    })
    
    cat(sprintf("  → %.1f s\n", ifelse(is.na(elapsed), 0, elapsed)))
    
    mp1_results[[i]] <- data.frame(
      pipeline   = "metaboprep1",
      platform   = plat,
      n_samples  = ns,
      n_features = nf,
      seconds    = elapsed,
      minutes    = elapsed / 60
    )
  }
}

# cleanup
unlink(mp1_tmp, recursive = TRUE, force = TRUE)

mp1_timing <- do.call(rbind, mp1_results)

## ── 6d. Print metaboprep1 table ──────────────────────────────────────────────

if (nrow(mp1_timing) > 0) {
  cat("\n── metaboprep1 timing results ────────────────────────────────────────────\n")
  print(mp1_timing[, c("platform", "n_samples", "n_features", "seconds", "minutes")],
        row.names = FALSE)
  
  mp1_csv <- file.path("scripts", "mp1_timing_results.csv")
  write.csv(mp1_timing, mp1_csv, row.names = FALSE)
  cat(sprintf("\nTable written to: %s\n", mp1_csv))
}

## ── 7. metaboprep1 line plot ──────────────────────────────────────────────────

if (nrow(mp1_timing) > 0 && any(!is.na(mp1_timing$minutes))) {
  
  mp1_colours <- c("Nightingale" = "#2196F3", "Metabolon (subset)" = "#4CAF50")
  mp1_timing$platform_label <- factor(
    sprintf("%s (~%s features)", mp1_timing$platform,
            format(mp1_timing$n_features, big.mark = ",")),
    levels = sprintf("%s (~%s features)", mp1_platforms$platform,
                     format(mp1_platforms$n_features, big.mark = ","))
  )
  
  p_mp1 <- ggplot(mp1_timing[!is.na(mp1_timing$minutes), ],
                  aes(x = n_samples, y = minutes,
                      colour = platform_label, group = platform_label)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_x_log10(
      breaks = mp1_n_samples,
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    scale_colour_manual(
      values = setNames(mp1_colours[mp1_platforms$platform],
                        levels(mp1_timing$platform_label)),
      name = "Platform"
    ) +
    labs(
      title    = "metaboprep1 pipeline — runtime by platform",
      subtitle = "Wall-clock time for legacy pipeline (incl. R startup ~5–10 s per run)",
      x        = "Number of samples (log scale)",
      y        = "Time (minutes)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      plot.subtitle    = element_text(colour = "grey40", size = 10,
                                      margin = margin(b = 10)),
      axis.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.title     = element_text(face = "bold"),
      legend.position  = "right"
    )
  
  out_png_mp1 <- file.path("scripts", "mp1_timing_by_platform.png")
  ggsave(out_png_mp1, plot = p_mp1, width = 8, height = 5.5, dpi = 150)
  cat(sprintf("Plot written to: %s\n", out_png_mp1))
}

## ── 8. Combined comparison bar chart (metaboprep2 vs metaboprep1, per platform) ──
## Bar chart: x = sample size (factor), y = minutes, faceted by platform,
## bars dodged by pipeline version.  coord_trans(y = "log10") is used so bars
## render correctly on a log scale without the zero-floor problem of
## scale_y_log10().  A small floor (1e-3 min) is applied so near-zero values
## don't break the transform.


## ── 2. Platform reference (needed for feat_lookup) ───────────────────────────
n_samples_vec <- c(100, 1000, 10000)
mp1_n_samples <- c(100, 1000)

mp1_platforms <- data.frame(
  platform   = c("Nightingale", "Olink"),
  n_features = c(250,            3000),
  stringsAsFactors = FALSE
)

## ── 3. Build combined data frame ─────────────────────────────────────────────
shared_platforms <- intersect(as.character(timing$platform),
                              as.character(mp1_timing$platform))

cmp_list <- lapply(shared_platforms, function(plat) {
  mp2_rows <- timing[as.character(timing$platform) == plat, ]
  mp1_rows <- mp1_timing[as.character(mp1_timing$platform) == plat, ]
  if (nrow(mp2_rows) == 0 || nrow(mp1_rows) == 0) return(NULL)
  mp2_rows$pipeline <- "metaboprep2"
  rbind(mp2_rows[, c("pipeline", "platform", "n_samples", "minutes")],
        mp1_rows[,  c("pipeline", "platform", "n_samples", "minutes")])
})
combined <- do.call(rbind, Filter(Negate(is.null), cmp_list))

# add NA placeholder rows for mp1 at sample sizes it wasn't run at
mp2_only_sizes <- setdiff(n_samples_vec, mp1_n_samples)
na_rows <- do.call(rbind, lapply(shared_platforms, function(plat) {
  data.frame(pipeline  = "metaboprep1",
             platform  = plat,
             n_samples = mp2_only_sizes,
             minutes   = NA_real_)
}))
combined <- rbind(combined, na_rows)

combined$pipeline  <- factor(combined$pipeline,
                             levels = c("metaboprep1", "metaboprep2"))
combined$platform  <- factor(combined$platform,
                             levels = shared_platforms)
combined$n_samples <- factor(combined$n_samples,
                             levels = sort(unique(combined$n_samples)),
                             labels = format(sort(unique(combined$n_samples)),
                                             big.mark = ",", scientific = FALSE))

# apply floor so coord_trans(log10) doesn't hit -Inf; NA rows kept for labels
combined$minutes_floored <- ifelse(is.na(combined$minutes), 1e-3,
                                   pmax(combined$minutes, 1e-3))

feat_lookup <- setNames(mp1_platforms$n_features, mp1_platforms$platform)
combined$platform_label <- factor(
  sprintf("%s\n(~%s features)", combined$platform,
          format(feat_lookup[as.character(combined$platform)], big.mark = ",")),
  levels = sprintf("%s\n(~%s features)", shared_platforms,
                   format(feat_lookup[shared_platforms], big.mark = ","))
)

## ── 8. Combined comparison bar chart (metaboprep2 vs metaboprep1, per platform) ──
## Bar chart: x = sample size (factor), y = minutes, faceted by platform,
## bars dodged by pipeline version.  coord_trans(y = "log10") is used so bars
## render correctly on a log scale without the zero-floor problem of
## scale_y_log10().  A small floor (1e-3 min) is applied so near-zero values
## don't break the transform.

# shared platforms (Nightingale, Metabolon) at the overlapping sample sizes
shared_platforms <- intersect(as.character(timing$platform),
                              as.character(mp1_timing$platform))

cmp_list <- lapply(shared_platforms, function(plat) {
  mp2_rows <- timing[as.character(timing$platform) == plat, ]
  mp1_rows <- mp1_timing[as.character(mp1_timing$platform) == plat, ]
  if (nrow(mp2_rows) == 0 || nrow(mp1_rows) == 0) return(NULL)
  mp2_rows$pipeline <- "metaboprep2"
  rbind(mp2_rows[, c("pipeline", "platform", "n_samples", "minutes")],
        mp1_rows[,  c("pipeline", "platform", "n_samples", "minutes")])
})
combined <- do.call(rbind, Filter(Negate(is.null), cmp_list))

# add NA placeholder rows for mp1 at sample sizes it wasn't run at
mp2_only_sizes <- setdiff(n_samples_vec, mp1_n_samples)
na_rows <- do.call(rbind, lapply(shared_platforms, function(plat) {
  data.frame(pipeline  = "metaboprep1",
             platform  = plat,
             n_samples = mp2_only_sizes,
             minutes   = NA_real_)
}))
combined <- rbind(combined, na_rows)

combined$pipeline  <- factor(combined$pipeline,
                             levels = c("metaboprep1", "metaboprep2"))
combined$platform  <- factor(combined$platform,
                             levels = shared_platforms)
combined$n_samples <- factor(combined$n_samples,
                             levels = sort(unique(combined$n_samples)),
                             labels = format(sort(unique(combined$n_samples)),
                                             big.mark = ",", scientific = FALSE))

# apply floor so log10 transform doesn't hit -Inf; NA rows kept for labels
combined$minutes_floored <- ifelse(is.na(combined$minutes), 1e-3,
                                   pmax(combined$minutes, 1e-3))
# pre-compute log10 so geom_col bars run from log10(floor) to log10(value),
# avoiding the coord_trans(log10) + ymin=0 → -Inf rendering problem
combined$log_minutes <- log10(combined$minutes_floored)

feat_lookup <- setNames(mp1_platforms$n_features, mp1_platforms$platform)
combined$platform_label <- factor(                                                                                                                                                                       
  sprintf("%s features", format(feat_lookup[as.character(combined$platform)], big.mark = ",")),                                                                                                                         
  levels = sprintf("%s features", format(feat_lookup[shared_platforms], big.mark = ","))                                                                                                                                
)                                   

pipeline_colours <- c("metaboprep2" = "#1565C0", "metaboprep1" = "#9E9E9E")

# Shift log10 values so that every bar is non-negative (geom_col always draws
# from 0 upward, so we set 0 = y_floor_log and express all heights relative to
# that baseline).  The axis labels are back-translated to real time units.
y_floor_log  <- log10(1e-3)                          # baseline ≈ 0.06 s
combined$bar_height <- combined$log_minutes - y_floor_log  # always >= 0

y_breaks_min <- c(1/60, 10/60, 1, 10, 1000)           # 1 s, 10 s, 1 min, 10 min, 100 min
y_breaks_adj <- log10(y_breaks_min) - y_floor_log     # shifted to match bar_height scale
y_labels     <- ifelse(y_breaks_min >= 1,
                       paste0(round(y_breaks_min), " min"),
                       paste0(round(y_breaks_min * 60), " sec"))

p_cmp <- ggplot(combined[combined$platform=="Olink",],
                aes(x = n_samples, y = bar_height, fill = pipeline)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(
    aes(label = ifelse(is.na(minutes), "not\nrun", ""),
        y     = 0),
    position = position_dodge(width = 0.75),
    size = 2.8, colour = "grey50", vjust = -0.5, lineheight = 0.85
  ) +
  facet_wrap(~ platform_label, ncol = 2, scales = "free_y") +
  scale_y_continuous(
    breaks = y_breaks_adj,
    labels = y_labels,
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_manual(
    values = pipeline_colours,
    name   = "Pipeline",
    labels = c("metaboprep2" = "metaboprep 2", "metaboprep1" = "metaboprep 1")
  ) +
  labs(
    title    = "metaboprep 1 vs 2 — QC pipeline runtime",
    x        = "Number of samples",
    y        = "Time (log scale)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(colour = "grey40", size = 10,
                                    margin = margin(b = 10)),
    axis.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text       = element_text(face = "bold", size = 11),
    legend.title     = element_text(face = "bold"),
    legend.position  = "right"
  )
p_cmp

out_png_cmp <- file.path("scripts", "timing_metaboprep1_vs_2.png")
ggsave(out_png_cmp, plot = p_cmp, width = 9, height = 9, dpi = 150)
cat(sprintf("Plot written to: %s\n", out_png_cmp))




