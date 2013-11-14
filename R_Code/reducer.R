#! /usr/bin/env Rscript

con <- file("stdin", open = "r")
while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
  if(grepl("library", line, fixed=TRUE)) next
  cat(line, "\n")
}
close(con)
