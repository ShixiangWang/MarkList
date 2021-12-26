#!/usr/bin/env Rscript
# Check GDC manifest file, remove downloaded IDs from the file

message("Usage: update_manifest.R <manifest_file_path> <file_dir_path> [--update]")
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

# Check file one by one
cli_alert_info("Checking file stats")
notOK = vector(length = nrow(df))
cli_progress_bar("Walking files: ", total = nrow(df), type = "tasks")

for (i in seq_len(nrow(df))) {
  f = file.path(dir_path, df$id[i], df$filename[i])
  full_size = df$size[i]
  if (file.exists(f)) {
    current_size = file.info(f)$size
    if (full_size != current_size) {
      msg = glue::glue("File with idx {i} has been downlowned {round(current_size/full_size, 2)*100}%")
      notOK[i] = TRUE
    } else {
      msg = "File is downloaded"
    }
  } else {
    msg = glue::glue("File with idx {i} does not exist yet")
    notOK[i] = TRUE
  }
  Sys.sleep(0.1)
  cli_progress_update(status = msg)
}
cli_progress_done()

todo = sum(notOK)
done = length(notOK) - todo
cli_alert_info("{done} file(s) done, {todo} file(s) to download")

if (length(args) > 2 && args[3] == "--update") {
  if (todo == 0) {
    cli_alert_success("All files recorded in manifest have been downloaded, no need to go.")
    quit()
  }
  
  if (done == 0) {
    cli_alert_success("All files recorded in manifest have not been downloaded, please stick to current manifest.")
    quit()
  }
  
  cli_alert_info("Generating new manifest file for unfinished tasks")
  bk_path = path.expand(file.path("~", ".gdc_manifest_bk", paste0(basename(path[1]), "_bk")))
  cli_alert_info("Backup current manifest file to {bk_path}")
  if (!dir.exists(dirname(bk_path))) dir.create(dirname(bk_path), recursive = TRUE)
  write_tsv(df, bk_path)
  cli_alert_info("Update current manifest file {path[1]}")
  write_tsv(df[notOK, ], path.expand(path[1]))
  
  cli_alert_success("Done")
}
