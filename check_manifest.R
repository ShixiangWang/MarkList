#!/usr/bin/env Rscript
# Check if all files recorded in GDC manifest file have been downloaded

message("Usage: check_manifest.R <manifest_file_path> <file_dir_path>")
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  message("At least 2 parameters required, see usage above")
  quit(status = -1)
}
path = args[1:2]
#path = c("~/Downloads/gdc_manifest_20210723_055904.txt", "~/Downloads")
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("cli"))

df = read_tsv(path[1], col_names = TRUE, col_types = cols())
dir_path = path.expand(path[2])

# Check files in batch
cli_alert_info("Total {nrow(df)} records")
cli_alert_info("Checking bam file stats")
bam_files = file.path(dir_path, df$id, df$filename)
exists = file.exists(bam_files)
cli_alert_info("Non-exist bam file number: {nrow(df) - sum(exists)}")
notdone = file.info(bam_files[exists])$size != df$size[exists]
cli_alert_info("Unfinished bam file number: {sum(notdone)}")
cli_alert_success("Totol bam files to be downloaded: {nrow(df) - sum(exists) + sum(notdone)}")
message()
cli_alert_info("Checking bai file stats")
bai_files = sub("bam$", "bai", bam_files)
exists = file.exists(bai_files)
cli_alert_success("Totol bai files to be downloaded: {nrow(df) - sum(exists)}")

