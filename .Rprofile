.First <- function() {
  if (grepl("Windows", Sys.getenv("OS"))) {
    .libPaths(paste(getwd(), "libWin", sep = "/")) 
    Sys.setenv(USERPROFILE=Sys.getenv("HOME"))
  } else {
    .libPaths(paste(getwd(), "libLinux", sep = "/"))
  }
}
.First()