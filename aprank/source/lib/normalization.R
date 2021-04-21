# BSD 2-Clause License
# 
# Copyright (c) 2021, Alejandro Ricci (aricci@iib.unsam.edu.ar), Fernán Agüero (fernan@iib.unsam.edu.ar)
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   
#   1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#          SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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