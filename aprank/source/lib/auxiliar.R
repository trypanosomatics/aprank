############################-
#### Required libraries ####
############################-
if (!require(data.table, quietly = TRUE)) {
    writeLines("Installing library 'data.table' for R")
    install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = T)
    library(data.table)
}

###################-
#### Functions ####
###################-
getSizeOfLongestSegmentAboveThreshold <- function(data, threshold) {
    x <- as.numeric(data >= threshold)
    
    max <- 0
    count <- 0
    if (length(x) > 0) {
        for (i in 1:length(x)) {
            if (x[i] == 1) {
                count <- count + 1
            } else {
                if (count > max) {
                    max <- count
                }
                count <- 0
            }
        }
    }
    
    max
}
