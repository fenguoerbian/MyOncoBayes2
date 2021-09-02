library(devtools)
library(testthat)

load_all()

Sys.setenv(VDIFFR_LOG_PATH="vdiffr.log")

snaps <- test_path("_snaps")

cat("INFO: Removing reference plots in", snaps, "..\n")

snap_dirs <- dir(snaps, "plot.*", full.names=TRUE)

for(s in snap_dirs) {
    file.remove(dir(s, ".*svg", full.names=TRUE))
}

Sys.setenv(TEST_VDIFFR="true")
Sys.setenv(NOT_CRAN="true")

cat("INFO: Running plot tests to recreate reference plots...\n")
test(".", filter="plot", reporter = "summary")

cat("INFO: Running plot tests and check against created reference...\n")
test(".", filter="plot", reporter = "summary")

cat("Done.\n")

