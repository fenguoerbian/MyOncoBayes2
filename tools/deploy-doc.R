#! /usr/bin/env Rscript

args <-  commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (job id with artifacts)\n")
}

job_id <- args[1]

cat("Must be run on DaVinci for appropiate rsconnect user setup.\nYou have been warned!\n")

cat("Downloading artifacts from job", job_id, "...\n")

target_dir <- tempfile("docs-")
dir.create(target_dir)

artifacts <- file.path(target_dir, "artifacts.zip")

h <- curl::new_handle()
source("~/.ssh/gitlab_token")
curl::handle_setheaders(h, "PRIVATE-TOKEN" = gitlab_token)
curl::curl_download(paste0("https://gitlabce.apps.dit-prdocp.novartis.net/api/v4/projects/6/jobs/", job_id, "/artifacts"), destfile=artifacts, handle=h)

cwd <- getwd()
setwd(target_dir)

cat("Unpacking artifacts", job_id, "...\n")

contents <- unzip("artifacts.zip", list=TRUE)
doc_files <- grep("^docs", contents$Name, value=TRUE)
unzip("artifacts.zip", files=doc_files)

library(here)

cat("Uploading new web-site...\n")

rsconnect::deployApp(
  appDir = "./docs",
  # appFileManifest = app_files,
  appPrimaryDoc = "index.html",
  appName = "OncoBayes2",
  appTitle = "OncoBayes2",
  contentCategory = "site",
  account = "weberse2",
  server = "rsconnect-prod.dit.eu.novartis.net",
  upload = TRUE,
  forceUpdate = TRUE,
  ##appId = "825",
  logLevel = "verbose"
)

cat("Cleaning up...\n")
setwd(cwd)
unlink(target_dir, recursive=TRUE, force=TRUE)

cat("Done.\n")
