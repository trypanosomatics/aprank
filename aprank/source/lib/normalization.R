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
normalizeLinear <- function(x) {
    #Normalize data as a linear function between its max and min
    M = max(x)
    m = min(x)
    
    if (M == m) {
        output <- rep(0, length(x))
    } else {
        output <- (x - m) / (M - m)
    }
    
    output
}
normalizeLinearFixed <- function(x, m, M) {
    #Normalize data as a linear function between a given max and min
    if (M == m) {
        output <- rep(0, length(x))
    } else {
        output <- (x - m) / (M - m)
    }

    output[output > 1] <- 1
    output[output < 0] <- 0
    
    output
}
normalizeLinearMinRange <- function(x, min_m, min_M) {
    #Give a minimum range to normalize data with for when its own range is smaller
    M = max(x)
    m = min(x)
    
    if ((M - m) >= (min_M - min_m)) {
        output <- normalizeLinear(x)
    } else {
        output <- normalizeLinearFixed(x, min_m, min_M)
    }
    
    output
}

normalizeSigmoid09 <- function(x, b) {
    #b value is the one who will be 0.9 after the normalization
    output <- -1 + 2/(1 + 20^(-x/b))
    
    output
} 
normalizeSigmoid05 <- function(x, b) {
    #b value is the one who will be 0.5 after the normalization
    output <- -1 + 2/(1 + 3^( -x/b))
    
    output
} 