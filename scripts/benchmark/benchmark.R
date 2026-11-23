# Benchmark of omiprep quality_control()
#
# Part 1: runtime and peak memory by platform, exact and fast mode (Fig 3)
# Part 2: speed-up over legacy metaboprep v1.0.1, one change at a time (Fig 4)
# Part 3: the PCA step on its own, full vs truncated SVD (S7 Fig)
#
# All data are simulated with omiprep::simulate_omics() using fixed seeds, at
# feature counts typical of each platform (Nightingale 250, Metabolon 1,250,
# Olink 3,000, SomaLogic 11,000). 10% of values per feature are missing on
# average (half below a detection limit, half at random) and 10% of features
# are complete.
#
# Run from this folder on an otherwise idle machine:
#   Rscript benchmark.R
# Each QC run is timed in a fresh R session under /usr/bin/time -l, which also
# gives the peak memory (macOS; on Linux use /usr/bin/time -v and change the
# line that reads "maximum resident set size").
#
# Options (environment variables):
#   BENCH_PART     parts to run, e.g. "1,3" (default "1,2,3")
#   BENCH_CORES    cores for the multi-core runs (default 10)
#   BENCH_TIMEOUT  time limit per run in seconds (default 3600)
#   BENCH_MISS     average missingness per feature (default 0.10)
#   BENCH_ONLY     comma-separated run keys (without the seed) to run; others are skipped
#   FIG_ONLY=1     just redraw the figures from the saved results
#
# Results are added to results/benchmark_raw.csv after every run, and runs
# already in that file are not repeated, so the script can be stopped and
# restarted. Runs over the time limit are recorded as timeouts.
#
# Figure 3 and 4 runs were done with BENCH_CORES=10. The longest runs were
# rerun with longer limits (8 h for legacy at 1,000 x 11,000; 24 h for exact
# mode at 10,000 x 11,000); the limit given to each run is kept in its status.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---- paths and settings ----
args <- commandArgs(trailingOnly = FALSE)
SCRIPT <- normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1]))
HERE <- dirname(SCRIPT)
REPO <- normalizePath(file.path(HERE, "..", ".."))
RES_DIR <- file.path(HERE, "results")
FIG_DIR <- file.path(HERE, "figures")
for (d in c(RES_DIR, FIG_DIR)) dir.create(d, showWarnings = FALSE)

suppressMessages(pkgload::load_all(REPO, quiet = TRUE, export_all = TRUE, helpers = FALSE))
options(cli.default_handler = function(...) invisible(NULL))

PARTS   <- as.integer(strsplit(Sys.getenv("BENCH_PART", "1,2,3"), ",")[[1]])
CORES   <- as.integer(Sys.getenv("BENCH_CORES", 10))
TIMEOUT <- as.integer(Sys.getenv("BENCH_TIMEOUT", 3600))
MISS    <- as.numeric(Sys.getenv("BENCH_MISS", 0.10))
ONLY    <- strsplit(Sys.getenv("BENCH_ONLY"), ",")[[1]]
MP1_REF <- "v1.0.1"
MP1     <- Sys.getenv("MP1_DIR", file.path(tempdir(), paste0("metaboprep_", MP1_REF)))
RAW     <- file.path(RES_DIR, "benchmark_raw.csv")
SVD_RAW <- file.path(RES_DIR, "svd_timing.csv")

PLATFORMS <- data.frame(platform = c("Nightingale", "Metabolon", "Olink", "SomaLogic"),
                        p = c(250, 1250, 3000, 11000), stringsAsFactors = FALSE)
plat_lab <- function(p) sprintf("%s (%s features)", PLATFORMS$platform[match(p, PLATFORMS$p)],
                                formatC(p, format = "d", big.mark = ","))

# QC thresholds, as in legacy metaboprep
QC <- list(sample_missingness = 0.2, feature_missingness = 0.2, total_sum_abundance_sd = 5,
           outlier_udist = 5, tree_cut_height = 0.5, pc_outlier_sd = 5)

# Part 2 steps, each adding one change to the one before
STEPS <- data.frame(step = c("legacy metaboprep", "1 core, exact, full SVD", "+ 10 cores", "+ fast = TRUE",
                             "+ truncated SVD (irlba)"),
                    engine = c("legacy", rep("omiprep", 4)), cores = c(1, 1, CORES, CORES, CORES),
                    fast = c(FALSE, FALSE, FALSE, TRUE, TRUE), pca = c("full", "full", "full", "full", "irlba"),
                    stringsAsFactors = FALSE)
step_id <- function(d) paste(d$engine, d$cores, as.logical(d$fast), d$pca)

simulate <- function(n, p, seed) {
  simulate_omics(n, p, miss_rate = MISS, miss_mech = "mixed", complete_frac = 0.1, seed = seed)$data
}

theme_sub <- function() theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95"),
        legend.position = "bottom")
write_res <- function(df, name) utils::write.csv(df, file.path(RES_DIR, paste0(name, ".csv")), row.names = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

environment_record <- function() {
  sh <- function(cmd) tryCatch(trimws(system(cmd, intern = TRUE)[1]), error = function(e) NA_character_)
  si <- sessionInfo()
  data.frame(
    item = c("cpu", "logical_cores", "ram_gb", "os", "r_version", "blas", "lapack",
             "omiprep_version", "omiprep_commit", "irlba", "nFactors", "date"),
    value = c(
      sh("sysctl -n machdep.cpu.brand_string 2>/dev/null || lscpu | grep 'Model name' | sed 's/.*: *//'"),
      parallel::detectCores(),
      round(as.numeric(sh("sysctl -n hw.memsize 2>/dev/null || echo 0")) / 1024^3, 1),
      paste(si$running), R.version.string, si$BLAS, si$LAPACK,
      as.character(packageVersion("omiprep", lib.loc = NULL)),
      sh(sprintf("git -C '%s' rev-parse --short HEAD", REPO)),
      as.character(packageVersion("irlba")), as.character(packageVersion("nFactors")),
      format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
    stringsAsFactors = FALSE)
}

# ---- one timed run (called as: Rscript benchmark.R worker key=value ...) ----
# engine=omiprep|legacy, n, p, cores, fast, pca=irlba|full, settings=default|matched, seed, out=<csv>
# "matched" uses the settings that reproduce legacy metaboprep (least-missingness, feature selection, all informative PCs).
worker <- function(a) {
  arg <- function(k, d) if (!is.null(a[[k]])) a[[k]] else d
  engine <- arg("engine", "omiprep"); n <- as.integer(arg("n", 1000)); p <- as.integer(arg("p", 1000))
  cores <- as.integer(arg("cores", 10)); fast <- as.logical(arg("fast", "TRUE")); pca <- arg("pca", "irlba")
  settings <- arg("settings", "default"); seed <- as.integer(arg("seed", 1)); out <- arg("out", "bench_row.csv")

  mat <- simulate(n, p, seed); invisible(gc())
  input_mb <- as.numeric(object.size(mat)) / 1024^2
  irlba_used <- NA; n_indep <- NA; steps <- NA_character_; n_layers <- NA

  if (engine == "omiprep") {
    if (pca == "full") {
      # swap the truncated SVD in pc_and_outliers() for a full one
      f <- S7::method(pc_and_outliers, Omiprep)
      b <- deparse(body(f), width.cutoff = 500L)
      target <- "irlba::prcomp_irlba(pcadata, n = k, center = FALSE, scale = FALSE)"
      stopifnot("could not find the irlba call to replace" = any(grepl(target, b, fixed = TRUE)))
      b <- sub(target, "stats::prcomp(pcadata, center = FALSE, scale = FALSE)", b, fixed = TRUE)
      body(f) <- parse(text = paste(b, collapse = "\n"))[[1]]
      S7::method(pc_and_outliers, Omiprep) <- f
    }
    obj <- Omiprep(data = mat, samples = data.frame(sample_id = rownames(mat)),
                   features = data.frame(feature_id = colnames(mat)))
    extra <- if (settings == "matched") list(feature_selection = "least_missingness", max_num_pcs = NULL) else list()
    invisible(gc(reset = TRUE))
    res <- NULL
    t <- system.time(txt <- utils::capture.output(suppressMessages(suppressWarnings(
      res <- do.call(quality_control, c(list(obj, cores = cores, fast = fast), QC, extra))))))[["elapsed"]]
    g <- gc()
    # step timings printed by quality_control()
    st <- grep("^\\s*[a-z_]+\\s+[0-9.]+\\s+[0-9.]+\\s*$", txt, value = TRUE)
    st <- st[!grepl("^\\s*step", st)]
    if (length(st)) steps <- paste(trimws(gsub("\\s+", "=", trimws(sub("\\s+[0-9.]+\\s*$", "", st)))), collapse = ";")
    fs <- res@feature_summary[, , "qc"]
    n_indep <- sum(fs["independent_features", ] == 1, na.rm = TRUE)
    n_kept  <- sum(!res@samples$excluded)
    # irlba is only used when 50 PCs is no more than 10% of the smaller dimension
    irlba_used <- min(50, n_kept - 1, n_indep) / min(n_kept, n_indep) <= 0.1 && pca == "irlba"
    result_mb <- as.numeric(object.size(res)) / 1024^2
    n_layers  <- dim(res@data)[3]
  } else {
    # legacy metaboprep v1.0.1, core QC function only
    mp <- suppressMessages(pkgload::load_all(MP1, quiet = TRUE, export_all = TRUE, attach = FALSE))$env
    invisible(gc(reset = TRUE))
    t <- system.time(invisible(utils::capture.output(suppressMessages(suppressWarnings(
      res <- mp$perform.metabolite.qc(wdata = mat, fmis = QC$feature_missingness, smis = QC$sample_missingness,
                                      tpa_out_SD = QC$total_sum_abundance_sd, outlier_udist = QC$outlier_udist,
                                      tree_cut_height = QC$tree_cut_height, PC_out_SD = QC$pc_outlier_sd))))))[["elapsed"]]
    g <- gc()
    result_mb <- as.numeric(object.size(res)) / 1024^2
  }
  row <- data.frame(engine = engine, n = n, p = p, cores = cores, fast = fast, pca = pca, settings = settings,
                    missing = MISS, seed = seed, seconds = t, input_mb = input_mb, result_object_mb = result_mb,
                    r_heap_peak_mb = sum(g[, ncol(g)]), n_layers = n_layers, independent_features = n_indep,
                    irlba_used = irlba_used, step_seconds = steps)
  utils::write.csv(row, out, row.names = FALSE)
}

cmd <- commandArgs(trailingOnly = TRUE)
if (length(cmd) && cmd[1] == "worker") {
  kv <- cmd[-1]
  worker(as.list(setNames(sub("^[^=]*=", "", kv), sub("=.*$", "", kv))))
  quit(save = "no")
}

# ---- running and recording ----
done <- if (file.exists(RAW)) utils::read.csv(RAW, stringsAsFactors = FALSE) else NULL
key_of <- function(a) paste(a$part, a$engine, a$n, a$p, a$cores, a$fast, a$pca, a$settings, a$seed, sep = "|")

run_one <- function(a) {
  tmp <- tempfile(fileext = ".csv"); err <- tempfile()
  wargs <- a[c("engine", "n", "p", "cores", "fast", "pca", "settings", "seed")]
  argv <- c(file.path(R.home("bin"), "Rscript"), SCRIPT, "worker",
            paste0(names(wargs), "=", unlist(wargs)), paste0("out=", tmp))
  status <- system2("/usr/bin/time", c("-l", argv), stdout = FALSE, stderr = err, timeout = TIMEOUT,
                    env = paste0("MP1_DIR=", MP1))
  e <- readLines(err, warn = FALSE)
  rss <- suppressWarnings(as.numeric(sub("^\\s*([0-9]+)\\s+maximum resident set size.*$", "\\1",
                                         grep("maximum resident set size", e, value = TRUE))))
  r <- if (file.exists(tmp)) {
    x <- utils::read.csv(tmp, stringsAsFactors = FALSE); x$status <- "ok"; x
  } else data.frame(wargs, seconds = NA, status = if (identical(status, 124L)) sprintf("timeout > %d s", TIMEOUT)
                    else paste("error:", paste(head(c(grep("Error", e, value = TRUE), tail(e, 1)), 2), collapse = " ")))
  r$peak_rss_mb <- if (length(rss)) rss[1] / 1024^2 else NA   # macOS reports bytes
  r$part <- a$part; r$step <- a$step %||% NA; r$key <- key_of(a)
  r
}

record <- function(r) {
  if (!is.null(done)) {
    for (cn in setdiff(names(done), names(r))) r[[cn]] <- NA
    for (cn in setdiff(names(r), names(done))) done[[cn]] <<- NA
    r <- r[, names(done)]
  }
  done <<- rbind(done, r)
  utils::write.csv(done, RAW, row.names = FALSE)
}

# up to 3 replicates; only 1 if a run takes over 10 minutes
run_reps <- function(a, max_reps = 3) {
  out <- list()
  for (s in seq_len(max_reps)) {
    a$seed <- s; k <- key_of(a)
    r <- if (!is.null(done) && k %in% done$key) done[done$key == k, ][1, ]
         else if (length(ONLY) && !sub("\\|[0-9]+$", "", k) %in% ONLY) data.frame(status = "skipped", seconds = NA)
         else { x <- run_one(a); record(x); x }
    out[[s]] <- r
    cat(sprintf("[part %d] %-28s n=%6d p=%5d cores=%2d fast=%-5s rep %d | %s\n", a$part, a$step %||% a$engine, a$n, a$p,
                a$cores, a$fast, s, if (r$status == "ok") sprintf("%.1f s, %.0f MB", r$seconds, r$peak_rss_mb) else r$status))
    if (r$status != "ok" || r$seconds > 600) break
  }
  do.call(rbind, lapply(out, function(x) x[, c("status", "seconds")]))
}

if (!nzchar(Sys.getenv("FIG_ONLY"))) {
  write_res(environment_record(), "benchmark_environment")

  # Part 1: if a size times out, the larger sizes for that platform are not run
  if (1 %in% PARTS) for (fast in c(TRUE, FALSE)) for (p in PLATFORMS$p) {
    for (n in c(100, 1000, 10000)) {
      r <- run_reps(list(part = 1, engine = "omiprep", n = n, p = p, cores = CORES, fast = fast, pca = "irlba",
                         settings = "default"))
      if (r$status[1] == "skipped") next
      if (r$status[1] != "ok") {
        for (n2 in c(100, 1000, 10000)[c(100, 1000, 10000) > n]) record(data.frame(part = 1, engine = "omiprep", n = n2, p = p,
          cores = CORES, fast = fast, pca = "irlba", settings = "default", seed = 1, seconds = NA, peak_rss_mb = NA, step = NA,
          status = paste0("not run (smaller n: ", r$status[1], ")"), key = paste("skip", n2, p, fast)))
        break
      }
    }
  }

  # Part 2: legacy metaboprep is fetched from GitHub at the pinned version
  if (2 %in% PARTS) {
    if (!dir.exists(MP1)) system2("git", c("clone", "-q", "https://github.com/MRCIEU/metaboprep.git", MP1))
    system2("git", c("-C", MP1, "fetch", "-q", "--tags", "--force"))
    stopifnot(system2("git", c("-C", MP1, "checkout", "-q", MP1_REF)) == 0)
    for (p in PLATFORMS$p) for (j in seq_len(nrow(STEPS)))
      run_reps(list(part = 2, step = STEPS$step[j], engine = STEPS$engine[j], n = 1000, p = p, cores = STEPS$cores[j],
                    fast = STEPS$fast[j], pca = STEPS$pca[j], settings = "matched"))
  }

  # Part 3: the PCA step alone, with the same preparation as pc_and_outliers()
  # (median imputation, scaled features). All p features go into the PCA, so
  # these are upper bounds; in QC only the independent features do.
  # 1 core for R, BLAS as installed. Seeds set the irlba start vector.
  if (3 %in% PARTS) {
    K <- 50
    sizes <- expand.grid(p = PLATFORMS$p, n = c(1000, 10000))[, c("n", "p")]
    time_it <- function(f) {
      invisible(gc(reset = TRUE)); t <- system.time(r <- f())[["elapsed"]]
      g <- gc(); list(r = r, seconds = t, heap_mb = sum(g[, which(colnames(g) == "max used") + 1]))
    }
    res <- list()
    for (i in seq_len(nrow(sizes))) {
      n <- sizes$n[i]; p <- sizes$p[i]
      x <- median_impute(data = simulate(n, p, seed = 1))
      x <- apply(x, 2, function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
      for (method in c("full SVD (prcomp)", "truncated SVD (irlba, 50 PCs)")) {
        f <- if (startsWith(method, "full")) function() stats::prcomp(x, center = FALSE, scale. = FALSE)$sdev[1:K]
             else function() irlba::prcomp_irlba(x, n = K, center = FALSE, scale. = FALSE)$sdev
        for (rep in 1:3) {
          set.seed(rep)
          o <- time_it(f)
          res[[length(res) + 1]] <- data.frame(n = n, p = p, method = method, rep = rep, seconds = o$seconds,
                                               heap_peak_mb = o$heap_mb, sdev1 = o$r[1], sdev50 = o$r[K],
                                               max_rel_diff_top50_sdev = NA)
          if (startsWith(method, "full")) full_sdev <- o$r
          else res[[length(res)]]$max_rel_diff_top50_sdev <- max(abs(o$r - full_sdev) / full_sdev)
          cat(sprintf("[part 3] n=%5d p=%5d %-30s rep %d | %8.2f s\n", n, p, method, rep, o$seconds))
          if (o$seconds > 600) break
        }
      }
      utils::write.csv(do.call(rbind, res), SVD_RAW, row.names = FALSE)
    }
  }
}

# ---- summaries ----
fmt_limit <- function(s) ifelse(s < 7200, paste(s / 60, "min"), paste(s / 3600, "h"))
fmt_dur <- function(s) ifelse(s < 60, paste(signif(s, 3), "s"),
                       ifelse(s < 3600, paste(signif(s / 60, 3), "min"), paste(signif(s / 3600, 3), "h")))
TIME_BREAKS <- c(1, 10, 60, 600, 3600, 8 * 3600, 24 * 3600)
TIME_LABELS <- c("1 s", "10 s", "1 min", "10 min", "1 h", "8 h", "24 h")

if (file.exists(RAW)) {
  raw <- utils::read.csv(RAW, stringsAsFactors = FALSE)
  raw$step[is.na(raw$step)] <- ""
  # label part 2 runs by their settings; runs from other configurations are kept but not plotted
  p2r <- raw$part == 2
  raw$step[p2r] <- STEPS$step[match(step_id(raw[p2r, ]), step_id(STEPS))]
  raw$step[p2r & is.na(raw$step)] <- with(raw[p2r & is.na(raw$step), ],
    sprintf("not in figure: %s cores, %s, %s SVD", cores, ifelse(as.logical(fast), "fast", "exact"),
            ifelse(pca == "irlba", "truncated", "full")))
  ok <- subset(raw, status == "ok")
  summ <- do.call(rbind, lapply(split(ok, list(ok$part, ok$engine, ok$n, ok$p, ok$cores, ok$fast, ok$pca, ok$settings, ok$step),
                                      drop = TRUE), function(d) data.frame(
    d[1, c("part", "step", "engine", "n", "p", "cores", "fast", "pca", "settings")], reps = nrow(d),
    median_s = median(d$seconds), q25_s = unname(quantile(d$seconds, 0.25)), q75_s = unname(quantile(d$seconds, 0.75)),
    min_s = min(d$seconds), max_s = max(d$seconds),
    median_peak_rss_mb = median(d$peak_rss_mb), median_r_heap_peak_mb = median(d$r_heap_peak_mb),
    input_mb = d$input_mb[1], result_object_mb = d$result_object_mb[1], n_layers = d$n_layers[1],
    irlba_used = d$irlba_used[1])))
  rownames(summ) <- NULL
  write_res(summ, "benchmark_summary")

  # runs that did not finish, with the time limit each was given
  not_run <- subset(raw, status != "ok")
  not_run$limit_s <- suppressWarnings(as.numeric(sub("^.*timeout > ([0-9]+) s.*$", "\\1", not_run$status)))
  not_run$limit_s[is.na(not_run$limit_s)] <- TIMEOUT

  env <- utils::read.csv(file.path(RES_DIR, "benchmark_environment.csv"), stringsAsFactors = FALSE)
  env_txt <- with(env, sprintf("%s, %s logical cores, %s GB RAM; %s; %s", value[item == "cpu"], value[item == "logical_cores"],
                               value[item == "ram_gb"], value[item == "os"], value[item == "r_version"]))

  # ---- Fig 3: runtime and memory by platform ----
  s1 <- subset(summ, part == 1)
  if (nrow(s1)) {
    s1$mode <- factor(ifelse(s1$fast, "fast = TRUE", "exact (fast = FALSE)"), c("exact (fast = FALSE)", "fast = TRUE"))
    s1$platform <- factor(plat_lab(s1$p), plat_lab(PLATFORMS$p))
    nr1 <- subset(not_run, part == 1 & grepl("^timeout", status))
    nr1$mode <- factor(ifelse(as.logical(nr1$fast), "fast = TRUE", "exact (fast = FALSE)"), levels(s1$mode))
    nr1$platform <- factor(plat_lab(nr1$p), levels(s1$platform))
    xs <- scale_x_log10(breaks = c(100, 1000, 10000), labels = scales::comma)
    pa <- ggplot(s1, aes(n, median_s, colour = platform, linetype = mode, shape = mode)) +
      geom_line() + geom_point(size = 2) + geom_errorbar(aes(ymin = min_s, ymax = max_s), width = 0.06) +
      (if (nrow(nr1)) geom_point(data = nr1, aes(as.numeric(n), limit_s), shape = 4, size = 3, show.legend = FALSE)) +
      xs + scale_y_log10(breaks = TIME_BREAKS, labels = TIME_LABELS) +
      labs(x = "samples", y = "QC runtime (log scale; median, bars = min-max)", colour = NULL, linetype = NULL, shape = NULL,
           title = "A. Runtime") + theme_sub()
    pb <- ggplot(s1, aes(n, median_peak_rss_mb, colour = platform, linetype = mode, shape = mode)) +
      geom_line() + geom_point(size = 2) + xs + scale_y_log10(labels = scales::comma) +
      labs(x = "samples", y = "peak memory of the R process, MB", colour = NULL, linetype = NULL, shape = NULL,
           title = "B. Peak memory") + theme_sub() +
      guides(colour = "none", linetype = "none", shape = "none")
    p1 <- (pa + pb + plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
      plot_annotation(title = sprintf("omiprep quality_control() runtime and memory by platform (%d cores)", CORES),
                      caption = paste0("Synthetic data, ", 100 * MISS,
                                       "% missing per feature on average (half below a detection limit, half at random). Default QC settings. ",
                                       "Each run in a fresh R process; replicates: 3 (1 if over 10 min).\n",
                                       if (nrow(nr1)) paste0("x: exact run did not finish within ", paste(unique(fmt_limit(nr1$limit_s)), collapse = " / "),
                                                             ", plotted at that limit; larger n not run after a timeout. "),
                                       "Peak memory = maximum resident set size of the R process (/usr/bin/time -l).\n", env_txt),
                      theme = theme(plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30")))
    ggsave(file.path(FIG_DIR, "benchmark_runtime_memory.png"), p1, width = 11, height = 5.5, dpi = 200, bg = "white")
  }

  # ---- Fig 4: where the speed-up comes from ----
  s2 <- subset(summ, part == 2 & step %in% STEPS$step)
  if (nrow(s2)) {
    s2$step <- factor(s2$step, rev(STEPS$step))
    s2$platform <- factor(plat_lab(s2$p), plat_lab(PLATFORMS$p))
    base <- setNames(s2$median_s[s2$step == "legacy metaboprep"], as.character(s2$platform[s2$step == "legacy metaboprep"]))
    s2$label <- sprintf("  %s%s", fmt_dur(s2$median_s),
                        ifelse(is.na(base[as.character(s2$platform)]) | s2$step == "legacy metaboprep", "",
                               sprintf("  (%.0fx)", base[as.character(s2$platform)] / s2$median_s)))
    nr2 <- subset(not_run, part == 2 & step %in% STEPS$step); nr2$step <- factor(nr2$step, levels(s2$step))
    nr2$platform <- factor(plat_lab(nr2$p), levels(s2$platform))
    p2 <- ggplot(s2, aes(median_s, step)) + geom_col(fill = "grey65") +
      geom_errorbarh(aes(xmin = min_s, xmax = max_s), height = 0.3) +
      geom_text(aes(label = label), hjust = -0.1, size = 2.8) +
      geom_text(data = nr2, aes(1, step, label = paste("did not finish within", fmt_limit(limit_s))),
                hjust = 0, size = 2.8, colour = "#D95F02") +
      facet_wrap(~ platform, ncol = 2) +
      scale_x_log10(breaks = TIME_BREAKS, labels = TIME_LABELS, expand = expansion(mult = c(0.02, 0.6))) +
      labs(x = "QC runtime (log scale; median, bars = min-max)", y = NULL,
           title = "QC pipeline processing speed breakdown (1,000 samples; settings matched to legacy metaboprep)",
           caption = paste0("Each row adds one change to the row above. Speed-up relative to legacy metaboprep v1.0.1 in brackets. ",
                            "Truncated SVD is only used when the PCA is large enough\n(50 PCs <= 10% of the smaller of samples and independent features); ",
                            "otherwise rows 4 and 5 run identical code. ", env_txt)) +
      theme_sub() + theme(plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30"))
    ggsave(file.path(FIG_DIR, "benchmark_speedup.png"), p2, width = 11, height = 6, dpi = 200, bg = "white")
  }
}

# ---- S7 Fig: PCA step, full vs truncated SVD ----
if (file.exists(SVD_RAW)) {
  out <- utils::read.csv(SVD_RAW, stringsAsFactors = FALSE)
  s <- do.call(rbind, lapply(split(out, list(out$n, out$p, out$method), drop = TRUE), function(d)
    data.frame(d[1, c("n", "p", "method")], median_s = median(d$seconds), min_s = min(d$seconds), max_s = max(d$seconds),
               heap_mb = median(d$heap_peak_mb), max_rel_diff = suppressWarnings(max(d$max_rel_diff_top50_sdev, na.rm = TRUE)))))
  s$samples <- factor(paste(formatC(s$n, format = "d", big.mark = ","), "samples"),
                      paste(formatC(sort(unique(s$n)), format = "d", big.mark = ","), "samples"))
  write_res(s, "svd_timing_summary")

  TB <- c(0.01, 0.1, 1, 10, 60, 600, 3600); TL <- c("0.01 s", "0.1 s", "1 s", "10 s", "1 min", "10 min", "1 h")
  pa <- ggplot(s, aes(p, median_s, colour = method)) + geom_line() + geom_point(size = 2) +
    geom_errorbar(aes(ymin = min_s, ymax = max_s), width = 0.04, show.legend = FALSE) + facet_wrap(~ samples) +
    scale_x_log10(breaks = PLATFORMS$p, labels = scales::comma) +
    scale_y_log10(breaks = TB, labels = TL) +
    labs(x = "features in the PCA", y = "PCA runtime (log scale; median, bars = min-max)", colour = NULL,
         title = "A. Runtime of the PCA step") + theme_sub()
  pb <- ggplot(s, aes(p, heap_mb, colour = method)) + geom_line() + geom_point(size = 2) + facet_wrap(~ samples) +
    scale_x_log10(breaks = PLATFORMS$p, labels = scales::comma) + scale_y_log10(labels = scales::comma) +
    labs(x = "features in the PCA", y = "peak R memory during the PCA, MB", colour = NULL,
         title = "B. Memory of the PCA step") + theme_sub()
  agree <- s$max_rel_diff[is.finite(s$max_rel_diff)]
  p3 <- (pa + pb + plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
    plot_annotation(caption = paste0(
      "Same preparation and calls as pc_and_outliers(): median imputation, z-scored features, stats::prcomp() vs irlba::prcomp_irlba(n = 50). ",
      "Synthetic data, 10% missing; all features enter the PCA\n(in QC only independent features do, so these are upper bounds). ",
      "Single R process, 1 core for R; BLAS as installed. Largest relative difference in the top 50 PC standard deviations: ",
      if (length(agree)) formatC(max(agree), format = "e", digits = 1) else "NA", ".\n",
      "In QC, irlba is used only when 50 PCs <= 10% of the smaller dimension; below that prcomp is used."),
      theme = theme(plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey30")))
  ggsave(file.path(FIG_DIR, "svd_timing.png"), p3, width = 11, height = 5.5, dpi = 200, bg = "white")
}
cat("done\n")
