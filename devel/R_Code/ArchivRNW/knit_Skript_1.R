library(knitr)
knit("../../dissRepo/abstracts/Chapter/6_rfh_vs_fh.Rnw", 
     "../../dissRepo/abstracts/Chapter/6_rfh_vs_fh.tex", 
     encoding = "UTF-8", envir = new.env())


pdfFiles <- list.files("figure/", full.names = TRUE)

for (fileName in pdfFiles) {
  file.copy(fileName, 
            paste("../../dissRepo/abstracts/Issue 3 - Robust FH Models/", fileName, sep = ""), 
            overwrite = TRUE)
}
