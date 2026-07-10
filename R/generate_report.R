#' @title Generate Output Report
#' @description
#' This function writes an output report
#' @param omiprep an object of class Omiprep
#' @param output_dir character, the directory to save to
#' @param output_filename character, default NULL i.e. create from input object
#' @param project character, name for the current project
#' @param format character, write either 'html' or 'pdf' report
#' @param template character, type of report to output only current option is "qc_report"
#' @include class_omiprep.R
#' @importFrom rmarkdown render
#' @export
generate_report <- new_generic("generate_report", c("omiprep"), function(omiprep, output_dir, output_filename=NULL, project = "Project", format="pdf", template="qc_report") { S7_dispatch() })
#' @name generate_report
method(generate_report, Omiprep) <- function(omiprep, output_dir, output_filename=NULL, project = "Project", format="pdf", template="qc_report") {

  # testing
  if (FALSE) {
    output_dir="/Users/xx20081/git/omiprep/inst/rmarkdown/templates/qc_report/skeleton"
    output_filename=NULL
    project = "Project"
    format="pdf"
    template="qc_report"
  }

  # checks
  format   <- match.arg(format, choices = c("pdf","html"))
  template <- match.arg(template, choices = available_report_templates())
  stopifnot("\n'qc' data layer not found, have you run the quality_control() function on your Omiprep object? \n Run `dimnames(omiprep@data)[[3]]` to see current data layers" = "qc" %in% dimnames(omiprep@data)[[3]])

  # the qc_report template requires these packages to be installed; fail early with
  # an actionable message rather than deep inside rmarkdown::render()
  required_pkgs <- c("kableExtra", "dendextend", "glue", "scales")
  missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "\nThe following packages are required to generate a report but are not installed: ", paste(missing_pkgs, collapse = ", "),
      "\n Install them with:\n   install.packages(c(", paste(sprintf("'%s'", missing_pkgs), collapse = ", "), "))",
      call. = FALSE
    )
  }
  if (!rmarkdown::pandoc_available()) {
    stop(
      "\nGenerating a report requires pandoc, which was not found on this system.",
      "\n Install RStudio (which bundles pandoc) or pandoc itself: https://pandoc.org/installing.html",
      call. = FALSE
    )
  }
  if (format == "pdf" && !(requireNamespace("tinytex", quietly = TRUE) && tinytex::is_tinytex())) {
    stop(
      "\nGenerating a PDF report requires a working LaTeX installation via tinytex, which was not found.",
      "\n Install it with:\n   install.packages('tinytex')\n   tinytex::install_tinytex()",
      "\n Alternatively, use format = 'html' to generate a report without LaTeX.",
      call. = FALSE
    )
  }

  # ensure dir exists and normalise
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # name the report
  if (is.null(output_filename)) {
    outpath <- file.path(output_dir, clean_names(paste0(project, "_omiprep_", template)))
  } else {
    output_filename <- basename(output_filename)
    outpath <- file.path(output_dir, output_filename)
  }

  # correct file extension
  if (!grepl("(?i)\\.(html|pdf)$", outpath)) {
    outpath <- paste(outpath, format, sep=".")
  }

  # ensure dir exists
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)

  # get the template
  template_path <- system.file("rmarkdown", "templates", template, "skeleton", "skeleton.Rmd", package="omiprep")

  # render the report; clean = TRUE (the default) removes all rendering
  # intermediates on success, including sibling template files (child .Rmds,
  # latex_styles.sty, styles.css, omiprep_workflow.png) that rmarkdown copies
  # into intermediates_dir so relative references resolve. On a LaTeX failure,
  # tinytex already prints the relevant .log excerpt into the error itself, so
  # nothing further needs to be preserved on disk.
  rmarkdown::render(
    input             = template_path,
    output_file       = outpath,
    output_dir        = dirname(outpath),
    knit_root_dir     = dirname(outpath),
    intermediates_dir = dirname(outpath),
    params            = list(project = project, omiprep = omiprep),
    output_format     = paste(format, "document", sep="_"),
    envir             = new.env()  # Use a new environment to avoid conflicts
  )

  invisible(omiprep)
}
