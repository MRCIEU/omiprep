library(omiprep)
library(devtools)
load_all()


fp  <- file.path(Sys.getenv("BBS_METABOLON_RAW_DIR"), "2025_04_04_bbs_raw_metabolon_data.xlsx")
out <- file.path(Sys.getenv("TMP_DIR"), "omiprep_testing")

omi <- read_metabolon(fp, 
                      sheet = 5,
                      feature_sheet  = "Chemical Annotation", 
                      feature_id_col = "CHEM_ID",
                      sample_sheet   = "Sample Meta Data", 
                      sample_id_col  = "PARENT_SAMPLE_NAME")

omi <- quality_control(omi, 
                       cores = 10, 
                       fast  = TRUE)

generate_report(omi, 
                output_dir = out, 
                output_filename="omiprep_example", 
                project = "OmiPrep Example", 
                format="pdf", 
                template="qc_report")
