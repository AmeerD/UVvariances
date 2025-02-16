library(dplyr)
library(tidyr)
library(purrr)

resdir <- "./res/"

resfiles <- list.files(resdir)

loadres <- function(fpath, cols) {
  res <- read.delim(fpath, sep=" ", header=F)
  colnames(res) <- cols
  res <- res %>% mutate(sim = strsplit(gsub(".txt", "", gsub("./res/", "", fpath)), "_")[[1]][2])
  return(res)
}

estdf <- map(resfiles[grepl("^est", resfiles)], 
               ~loadres(paste(resdir, .x, sep=""),
                        c("idx", "sample.size", "ss.obs", "design", "outcome", "value"))) %>%
  list_rbind() %>%
  filter(!is.na(idx))

vardf <- map(resfiles[grepl("^var", resfiles)], 
             ~loadres(paste(resdir, .x, sep=""),
                      c("idx", "sample.size", "method", "design", "outcome", "runtime", "value"))) %>%
  list_rbind() %>%
  filter(!is.na(idx))


save(estdf, vardf, file="ACSresults.rda")
