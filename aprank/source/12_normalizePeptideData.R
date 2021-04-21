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

#######################-
#### **LIBRARIES** ####
#######################-
if (!require(data.table, quietly = TRUE)) {
    writeLines("Installing library 'data.table' for R")
    install.packages("data.table", repos = "http://cran.rstudio.com/", dependencies = T)
    library(data.table)
}

###########################-
#### **AUX FUNCTIONS** ####
###########################-
getOffsettedColumnsData <- function(output_per_aminoacid,
                                    column_to_offset,
                                    peptide_length,
                                    peptide_offset) {
    aminoacids_aux <- output_per_aminoacid[, .(start, end, protein_length)]
    for (i in 1:peptide_length) {
        new_column_aux <- paste(column_to_offset, "_", i, sep = "")
        aminoacids_aux[[new_column_aux]] <- c(output_per_aminoacid[[column_to_offset]][i:output_per_aminoacid[, .N]], rep(0, i - 1))
    }
    #Remove the incomplete rows (if peptide_lenght is 15 this would be the last 14 rows)
    aminoacids_aux <- aminoacids_aux[end <= protein_length]
    
    #Keep only the peptides I want based on the peptide_overlap
    aminoacids_aux <- aminoacids_aux[((start - 1) %% peptide_offset) == 0] 
    
    aminoacids_aux
}

##################-
#### **MAIN** ####
##################-
normalizePeptideData <- function(output_per_aminoacid, output_per_protein, output_per_kmer, output_per_peptide, output_per_repeat,
                                 output_method, normalized_output_per_peptide_file = "",
                                 peptide_length = 15, peptide_overlap = 14,
                                 Paircoil2_threshold = 0.5,
                                 Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty = 0, Coendemicity_peptide_start_penalty_proportion = 0, Coendemicity_peptide_max_penalty_proportion = 2/15,
                                 use_BepiPred = 1, use_Paircoil2 = 1, use_PredGPI = 1, use_SignalP = 1,
                                 use_TMHMM = 1, use_NetSurfp = 1,
                                 use_Iupred = 1, use_NetOglyc = 1, use_Xstream = 1, use_NetMHCIIpan = 1,
                                 use_IsoelectricPoint = 1, use_MolecularWeight = 1,
                                 use_SelfSimilarity = 1, use_CrossReactivity = 1, use_Coendemicity = 1) {
    ####################-
    #### **CONFIG** ####
    ####################-
    Iupred_threshold <- 0.5
    
    round_decimals <- 4
    
    ###################################################-
    #### **SET UP DATA FOR PEPTIDE NORMALIZATION** ####
    ###################################################-
    writeLines("Setting up data for peptide normalization...")
    
    output_per_aminoacid_aux <- output_per_aminoacid
    output_per_aminoacid_aux[1] <- output_per_aminoacid_aux[1] #unlink
    
    #Add the start, end and protein_length variables
    setnames(output_per_aminoacid_aux, "pos", "start")
    output_per_aminoacid_aux[, end := start + peptide_length - 1]
    protein_length_data <- output_per_aminoacid_aux[, .(protein_length = .N), by = .(id)]
    output_per_aminoacid_aux <- merge(output_per_aminoacid_aux,
                                      protein_length_data,
                                      by = c("id"))
    
    #Fetch the data for the different predictors from protein level
    output_per_aminoacid_aux <- merge(output_per_aminoacid_aux,
                                      output_per_protein[, c("id", "PredGPI_start", "SignalP_end")],
                                      by = c("id"))
    
    ################-
    #### Iupred ####
    ################-
    if (use_Iupred) {
        #Transform the Iupred data to later on calculate the proportion of values above a threshold
        output_per_aminoacid_aux[, Iupred := as.integer(Iupred >= Iupred_threshold)]
    }
    
    #####################-
    #### NetMHCIIpan ####
    #####################-
    if (use_NetMHCIIpan) {
        NetMHCIIpan_peptide_length <- nchar(output_per_peptide[1]$peptide)
        NetMHCIIpan_peptides <- output_per_peptide[, .(id, start, NetMHCIIpan_mean_rank = mean_rank)]
        output_per_aminoacid_aux <- merge(output_per_aminoacid_aux,
                                          NetMHCIIpan_peptides,
                                          by = c("id", "start"),
                                          all.x = T)
        output_per_aminoacid_aux[is.na(NetMHCIIpan_mean_rank), NetMHCIIpan_mean_rank := 50] #50 is the worst rank possible
        output_per_aminoacid_aux[, NetMHCIIpan := 1 - normalizeLinearFixed(NetMHCIIpan_mean_rank / 100, 0.05, 0.5)]
        output_per_aminoacid_aux[, NetMHCIIpan := round(NetMHCIIpan, round_decimals)]
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c("NetMHCIIpan_mean_rank")]
    }
    
    ###################-
    #### Paircoil2 ####
    ###################-
    if (use_Paircoil2) {
        #Transform the Paircoil2 data to later on calculate the proportion of values above a threshold
        output_per_aminoacid_aux[, Paircoil2 := as.integer(Paircoil2 >= Paircoil2_threshold)]
    }
    
    #################-
    #### Xstream ####
    #################-
    if (use_Xstream) {
        #Assign the Xstream data to each aa (each aa has the highest copy number repeat he's involved in)
        if (output_per_repeat[, .N] > 0) {
            aminoacids_with_repeats <- data.table(id = character(), start = numeric(), copy_number = numeric())
            
            for (repeat_i in 1:output_per_repeat[, .N]) {
                # repeat_for <- output_per_repeat[1008]
                # writeLines(sprintf("%s", repeat_i))
                repeat_for <- output_per_repeat[repeat_i]

                aminoacids_with_repeats <- rbindlist(list(aminoacids_with_repeats, 
                                                          data.table(id = repeat_for$id,
                                                                     start = c(repeat_for$start:repeat_for$end),
                                                                     copy_number = repeat_for$copy_number)))
            }
            
            #Keep only the greatest copy_number for each aa
            aminoacids_with_repeats <- aminoacids_with_repeats[, .(Xstream = max(copy_number)), by = .(id, start)]
            
            output_per_aminoacid_aux <- merge(output_per_aminoacid_aux,
                                              aminoacids_with_repeats,
                                              by = c("id", "start"),
                                              all.x = T)
            output_per_aminoacid_aux[is.na(Xstream), Xstream := 0]
        } else {
            output_per_aminoacid_aux$Xstream <- 0    
        }
    }
    
    ################################################-
    #### **CALCULATE & NORMALIZE PEPTIDE DATA** ####
    ################################################-
    writeLines("Normalizing peptide data...")
    
    ###Since looping for each possible peptide and subsetting its data was very slow, what I'm doing now is
    ###copying several times the columns for each predictor but with a shift so each row ends up having the data
    ###of all the aa from a given peptide. It's not pretty, but it works better.
    
    ###While copying all the columns at the beginning was cleaner, it also had problems with really big datasets,
    ###so for peptides I'm doing it in pieces
    
    #################-
    #### General ####
    #################-
    peptide_offset <- peptide_length - peptide_overlap
    
    #Create the output
    normalized_output_per_peptide <- output_per_aminoacid_aux[, .(id, start, end, protein_length, PredGPI_start, SignalP_end)]
    
    #Remove the incomplete rows (if peptide_lenght is 15 this would be the last 14 rows)
    normalized_output_per_peptide <- normalized_output_per_peptide[end <= protein_length]
    
    #Keep only the peptides I want based on the peptide_overlap
    normalized_output_per_peptide <- normalized_output_per_peptide[((start - 1) %% peptide_offset) == 0]
    
    #Calculate the sequence of each peptide
    column_aux <- "aa"
    aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                              column_to_offset = column_aux,
                                              peptide_length = peptide_length,
                                              peptide_offset = peptide_offset)
    #Calculate the sequence
    cols <- paste(column_aux, "_", c(1:peptide_length), sep = "")
    normalized_output_per_peptide$peptide <- apply(aminoacids_aux[, cols, with=FALSE], 1, paste, collapse = "")
    
    #Remove the extra columns to make space in memory
    output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c(column_aux), with=FALSE]
    
    ##################-
    #### BepiPred ####
    ##################-
    column_aux <- "BepiPred"
    if (use_BepiPred) {
        ## Calculated as: mean of values for each aa
        ## Normalization: the distribution of scores for the aa range between -2 and 2, so I'll normalize it between -1.5 and 1.5
        aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                                  column_to_offset = column_aux,
                                                  peptide_length = peptide_length,
                                                  peptide_offset = peptide_offset)
        #Normalize the data
        cols <- paste(column_aux, "_", c(1:peptide_length), sep ="")
        normalized_output_per_peptide$BepiPred <- rowMeans(aminoacids_aux[, cols, with=FALSE])
        normalized_output_per_peptide$BepiPred <- normalizeLinearFixed(normalized_output_per_peptide$BepiPred, -1.5, 1.5)
        normalized_output_per_peptide[, BepiPred := round(BepiPred, round_decimals)]
        
        #Remove the extra columns to make space in memory
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c(column_aux), with=FALSE]
    } else {
        normalized_output_per_peptide$BepiPred <- 0
    }
    
    ###################-
    #### Paircoil2 ####
    ###################-
    column_aux <- "Paircoil2"
    if (use_Paircoil2) {
        ## Calculated as: since it was transformed above, this calculates the proportion of values greater or equal than a threshold
        aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                                  column_to_offset = column_aux,
                                                  peptide_length = peptide_length,
                                                  peptide_offset = peptide_offset)
        #Normalize the data
        cols <- paste(column_aux, "_", c(1:peptide_length), sep ="")
        normalized_output_per_peptide$Paircoil2 <- rowMeans(aminoacids_aux[, cols, with=FALSE])
        normalized_output_per_peptide[, Paircoil2 := round(Paircoil2, round_decimals)]
        
        #Remove the extra columns to make space in memory
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c(column_aux), with=FALSE]
    } else {
        normalized_output_per_peptide$Paircoil2 <- 0
    }
    
    #################-
    #### PredGPI ####
    #################-
    if (use_PredGPI) {
        ##Calcualted as: Checks if the protein has a GPI and if this fragment is inside of it at least in part
        normalized_output_per_peptide[, PredGPI := as.integer((PredGPI_start > 0) & (end > PredGPI_start))]
    } else {
        normalized_output_per_peptide$PredGPI <- 0
    }
    normalized_output_per_peptide[, PredGPI_start := NULL]
    
    #################-
    #### SignalP ####
    #################-
    if (use_SignalP) {
        ##Calculated as: Checks if the protein has a SignalP and if this fragment is inside of it at least in part
        normalized_output_per_peptide[, SignalP := as.numeric((SignalP_end > 0) & (start < SignalP_end))]
    } else {
        normalized_output_per_peptide$SignalP <- 0
    }
    normalized_output_per_peptide[, SignalP_end := NULL]
    
    ###############-
    #### TMHMM ####
    ###############-
    column_aux <- "TMHMM"
    if (use_TMHMM) {
        aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                                  column_to_offset = column_aux,
                                                  peptide_length = peptide_length,
                                                  peptide_offset = peptide_offset)
        #Checks if the peptide is some part inside a membrane helix
        cols <- paste(column_aux, "_", c(1:peptide_length), sep = "")
        normalized_output_per_peptide$TMHMM <- rowSums(aminoacids_aux[, cols, with=FALSE])
        normalized_output_per_peptide[, TMHMM := as.integer(TMHMM > 0)]
        
        #Remove the extra columns to make space in memory
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c(column_aux), with=FALSE]
    } else {
        normalized_output_per_peptide$TMHMM <- 0
    }
    
    ##################-
    #### NetSurfp ####
    ##################-
    if (use_NetSurfp) {
        ## NetSurfp_RSA
        column_aux <- "NetSurfp_RSA"
        aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                                  column_to_offset = column_aux,
                                                  peptide_length = peptide_length,
                                                  peptide_offset = peptide_offset)
        #Normalize the data
        cols <- paste(column_aux, "_", c(1:peptide_length), sep = "")
        normalized_output_per_peptide$NetSurfp_RSA <- rowMeans(aminoacids_aux[, cols, with=FALSE])
        normalized_output_per_peptide[, NetSurfp_RSA := round(NetSurfp_RSA, round_decimals)]
        
        ## NetSurfp_AH
        column_aux <- "NetSurfp_AH"
        aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                                  column_to_offset = column_aux,
                                                  peptide_length = peptide_length,
                                                  peptide_offset = peptide_offset)
        #Normalize the data
        cols <- paste(column_aux, "_", c(1:peptide_length), sep = "")
        normalized_output_per_peptide$NetSurfp_AH <- rowMeans(aminoacids_aux[, cols, with=FALSE])
        normalized_output_per_peptide[, NetSurfp_AH := round(NetSurfp_AH, round_decimals)]
        
        ## NetSurfp_BS
        column_aux <- "NetSurfp_BS"
        aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                                  column_to_offset = column_aux,
                                                  peptide_length = peptide_length,
                                                  peptide_offset = peptide_offset)
        #Normalize the data
        cols <- paste(column_aux, "_", c(1:peptide_length), sep = "")
        normalized_output_per_peptide$NetSurfp_BS <- rowMeans(aminoacids_aux[, cols, with=FALSE])
        normalized_output_per_peptide[, NetSurfp_BS := round(NetSurfp_BS, round_decimals)]
        
        #Remove the extra columns to make space in memory
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c("NetSurfp_RSA"), with=FALSE]
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c("NetSurfp_AH"), with=FALSE]
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c("NetSurfp_BS"), with=FALSE]
    } else {
        normalized_output_per_peptide$NetSurfp_RSA <- 0
        normalized_output_per_peptide$NetSurfp_AH <- 0
        normalized_output_per_peptide$NetSurfp_BS <- 0
    }
    
    ################-
    #### Iupred ####
    ################-
    column_aux <- "Iupred"
    if (use_Iupred) {
        ## Calcualted as: since it was transformed above, this calculates the proportion of values greater or equal than 0.5
        aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                                  column_to_offset = column_aux,
                                                  peptide_length = peptide_length,
                                                  peptide_offset = peptide_offset)
        #Normalize the data
        cols <- paste(column_aux, "_", c(1:peptide_length), sep = "")
        normalized_output_per_peptide$Iupred <- rowMeans(aminoacids_aux[, cols, with=FALSE])
        normalized_output_per_peptide[, Iupred := round(Iupred, round_decimals)]
        
        #Remove the extra columns to make space in memory
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c(column_aux), with=FALSE]
    } else {
        normalized_output_per_peptide$Iupred <- 0
    }
    
    ##################-
    #### NetOglyc ####
    ##################-
    column_aux <- "NetOglyc"
    if (use_NetOglyc) {
        aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                                  column_to_offset = column_aux,
                                                  peptide_length = peptide_length,
                                                  peptide_offset = peptide_offset)
        #Normalize the data
        cols <- paste(column_aux, "_", c(1:peptide_length), sep = "")
        normalized_output_per_peptide$NetOglyc <- as.integer(rowSums(aminoacids_aux[, cols, with=FALSE]) > 0) #finds out if at least 1 column is 1
        
        #Remove the extra columns to make space in memory
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c(column_aux), with=FALSE]
    } else {
        normalized_output_per_peptide$NetOglyc <- 0
    }
    
    #################-
    #### XStream ####
    #################-
    column_aux <- "Xstream"
    if (use_Xstream) {
        aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                                  column_to_offset = column_aux,
                                                  peptide_length = peptide_length,
                                                  peptide_offset = peptide_offset)
        #Normalize the data
        cols <- paste(column_aux, "_", c(1:peptide_length), sep = "")
        normalized_output_per_peptide$Xstream <- rowMeans(aminoacids_aux[, cols, with=FALSE])
        normalized_output_per_peptide[, Xstream := normalizeSigmoid09(Xstream, 5)] #5 repeats = 0.9
        normalized_output_per_peptide[, Xstream := round(Xstream, round_decimals)]
        
        #Remove the extra columns to make space in memory
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c(column_aux), with=FALSE]        
    } else {
        normalized_output_per_peptide$Xstream <- 0
    }
    
    #####################-
    #### NetMHCIIpan ####
    #####################-
    column_aux <- "NetMHCIIpan"
    if (use_NetMHCIIpan) {
        aminoacids_aux <- getOffsettedColumnsData(output_per_aminoacid = output_per_aminoacid_aux,
                                                  column_to_offset = column_aux,
                                                  peptide_length = peptide_length,
                                                  peptide_offset = peptide_offset)
        #Normalize the data
        #Here I limit the columns because I only want the MHCII fragments that are inside the peptide
        ##By subtracting NetMHCIIpan_peptide_length I end up with the data from the shorter peptides that are inside the longer one (of peptide_length length)
        cols <- paste(column_aux, "_", c(1:(peptide_length - NetMHCIIpan_peptide_length + 1)), sep = "")
        normalized_output_per_peptide$NetMHCIIpan <- rowMeans(aminoacids_aux[, cols, with=FALSE])
        normalized_output_per_peptide[, NetMHCIIpan := round(NetMHCIIpan, round_decimals)]
        
        #Remove the extra columns to make space in memory
        output_per_aminoacid_aux <- output_per_aminoacid_aux[, -c(column_aux), with=FALSE]
    } else {
        normalized_output_per_peptide$NetMHCIIpan <- 0
    }
    
    ##############################################-
    #### Isoelectric point & Molecular weight ####
    ##############################################-
    #These aren't calculated at peptide level
    
    ##########################################################-
    #### Self Similarity, Cross Reactivity & Coendemicity ####
    ##########################################################-
    if (use_SelfSimilarity | use_CrossReactivity | use_Coendemicity) {
        kmer_length <- nchar(output_per_kmer[1]$kmer)
        kmer_amount <- peptide_length - kmer_length + 1
        #Split all the peptides in their possible kmers, create a new table with peptide start kmer
        unique_peptides <- unique(normalized_output_per_peptide$peptide)
        output_per_kmer_per_peptide <- data.table(peptide = rep(unique_peptides, each = kmer_amount),
                                                  kmer_start = rep(c(1:kmer_amount), length(unique_peptides)))
        output_per_kmer_per_peptide[, kmer_end := kmer_start + kmer_length - 1]
        output_per_kmer_per_peptide[, kmer := substr(peptide, kmer_start, kmer_end)]
        #Add kmer data
        output_per_kmer_per_peptide <- merge(output_per_kmer_per_peptide,
                                             unique(output_per_kmer[, .(kmer, quantity_in_proteome, quantity_in_host_proteome, quantity_in_coendemic_proteome)]),
                                             by = "kmer")

        #Save the data I'll need to calculate Coendemicity
        output_per_kmer_per_peptide_aux <- output_per_kmer_per_peptide[, .(kmer, peptide, quantity_in_coendemic_proteome)]
        
        #Combine data per peptide
        #Here, instead of calculating the amount of each kmer in the peptide (like I did in the protein) I'll just add them all up
        output_per_kmer_per_peptide <- output_per_kmer_per_peptide[, .(SelfSimilarity = sum(normalizeSigmoid05(quantity_in_proteome - 1, 1)) / .N,
                                                                       CrossReactivity = sum(normalizeSigmoid05(quantity_in_host_proteome, 1)) / .N), 
                                                                   by = .(peptide)]
        
        #Calculate the Coendemicity Penalty
        output_per_kmer_per_peptide_aux <- output_per_kmer_per_peptide_aux[, .(quantity = .N),
                                                                           by = .(kmer, peptide, quantity_in_coendemic_proteome)]
        output_per_kmer_per_peptide_aux$coendemic_quantity_aux <- 0L
        output_per_kmer_per_peptide_aux[(quantity_in_coendemic_proteome >= Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty) &
                                (quantity_in_coendemic_proteome > 0), coendemic_quantity_aux := quantity]
        output_per_kmer_per_peptide_aux <- output_per_kmer_per_peptide_aux[, .(kmer_amount = sum(quantity),
                                                       coendemic_kmer_amount = sum(coendemic_quantity_aux)),
                                                   by = .(peptide)]
        output_per_kmer_per_peptide_aux[, coendemic_proportion := coendemic_kmer_amount / kmer_amount]
        output_per_kmer_per_peptide_aux[, Coendemicity_Penalty := normalizeLinearFixed(coendemic_proportion,
                                                                           Coendemicity_peptide_start_penalty_proportion,
                                                                           Coendemicity_peptide_max_penalty_proportion)]
        output_per_kmer_per_peptide <- merge(output_per_kmer_per_peptide,
                                             output_per_kmer_per_peptide_aux[, .(peptide, Coendemicity_Penalty)],
                                             by = "peptide")
        
        #Round the data
        output_per_kmer_per_peptide[, SelfSimilarity := round(SelfSimilarity, round_decimals)]
        output_per_kmer_per_peptide[, CrossReactivity := round(CrossReactivity, round_decimals)]
        output_per_kmer_per_peptide[, Coendemicity_Penalty := round(Coendemicity_Penalty, round_decimals)]
        
        #Make sure that the ones that aren't used are 0
        if (!use_SelfSimilarity) {
            output_per_kmer_per_peptide$SelfSimilarity <- 0
        }
        if (!use_CrossReactivity) {
            output_per_kmer_per_peptide$CrossReactivity <- 0
        }
        if (!use_Coendemicity) {
            output_per_kmer_per_peptide$Coendemicity_Penalty <- 0
        }
        
        #Add to output
        normalized_output_per_peptide <- merge(normalized_output_per_peptide,
                                               output_per_kmer_per_peptide,
                                               by = "peptide")
    } else {
        normalized_output_per_peptide$SelfSimilarity <- 0
        normalized_output_per_peptide$CrossReactivity <- 0
        normalized_output_per_peptide$Coendemicity_Penalty <- 0
    }
    
    #########################-
    #### **OUTPUT DATA** ####
    #########################-
    writeLines("Readying Output...")
    
    #Sort data
    normalized_output_per_peptide <- normalized_output_per_peptide[order(id, start)]
    setcolorder(normalized_output_per_peptide, c("id", "start", "end", "peptide"))
    
    if (output_method %in% c("write", "both")) {
        write.table(normalized_output_per_peptide, file = normalized_output_per_peptide_file, col.names = T, row.names = F, sep = "\t", quote = T)
    }
    
    list_output <- list()
    if (output_method %in% c("list", "both")) {
        #Create the output
        list_output[["normalized_output_per_peptide"]] <- normalized_output_per_peptide
    }
    
    writeLines("Done!")
    
    list_output
}
normalizePeptideData_fromFile <- function(output_per_aminoacid_file, output_per_protein_file, output_per_kmer_file, output_per_peptide_file, output_per_repeat_file,
                                          output_method, normalized_output_per_peptide_file = "",
                                          peptide_length = 15, peptide_overlap = 14,
                                          Paircoil2_threshold = 0.5,
                                          Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty = 0, Coendemicity_peptide_start_penalty_proportion = 0, Coendemicity_peptide_max_penalty_proportion = 2/15,
                                          use_BepiPred = 1, use_Paircoil2 = 1, use_PredGPI = 1, use_SignalP = 1,
                                          use_TMHMM = 1, use_NetSurfp = 1,
                                          use_Iupred = 1, use_NetOglyc = 1, use_Xstream = 1, use_NetMHCIIpan = 1,
                                          use_IsoelectricPoint = 1, use_MolecularWeight = 1,
                                          use_SelfSimilarity = 1, use_CrossReactivity = 1, use_Coendemicity = 1) {
    output_per_aminoacid <- fread(output_per_aminoacid_file, header = T, sep = "\t", na.strings = NULL)
    output_per_protein <- fread(output_per_protein_file, header = T, sep = "\t", na.strings = NULL)
    output_per_kmer <- fread(output_per_kmer_file, header = T, sep = "\t", na.strings = NULL)
    output_per_peptide <- fread(output_per_peptide_file, header = T, sep = "\t", na.strings = NULL)
    output_per_repeat <- fread(output_per_repeat_file, header = T, sep = "\t", na.strings = NULL)
    
    normalizePeptideData(output_per_aminoacid = output_per_aminoacid, output_per_protein = output_per_protein, output_per_kmer = output_per_kmer, output_per_peptide = output_per_peptide, output_per_repeat = output_per_repeat,
                         output_method = output_method, normalized_output_per_peptide_file = normalized_output_per_peptide_file,
                         peptide_length = peptide_length, peptide_overlap = peptide_overlap,
                         Paircoil2_threshold = Paircoil2_threshold,
                         Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty = Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_peptide_start_penalty_proportion = Coendemicity_peptide_start_penalty_proportion, Coendemicity_peptide_max_penalty_proportion = Coendemicity_peptide_max_penalty_proportion,
                         use_BepiPred = use_BepiPred, use_Paircoil2 = use_Paircoil2, use_PredGPI = use_PredGPI, use_SignalP = use_SignalP,
                         use_TMHMM = use_TMHMM, use_NetSurfp = use_NetSurfp,
                         use_Iupred = use_Iupred, use_NetOglyc = use_NetOglyc, use_Xstream = use_Xstream, use_NetMHCIIpan = use_NetMHCIIpan,
                         use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight = use_MolecularWeight,
                         use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
}

# ##################-
# #### **CALL** ####
# ##################-
# ####################################-
# #### Config - Select Predictors ####
# ####################################-
# ## General
# use_BepiPred <- 1
# use_IsoelectricPoint <- 1
# use_Iupred <- 1
# use_MolecularWeight <- 1
# use_NetMHCIIpan <- 1
# use_NetOglyc <- 1
# use_NetSurfp <- 1
# use_Paircoil2 <- 1
# use_PredGPI <- 1
# use_SignalP <- 1
# use_TMHMM <- 1
# use_Xstream <- 1
# 
# use_SelfSimilarity <- 1
# use_CrossReactivity <- 1
# use_Coendemicity <- 1
# 
# #################################-
# #### Config - Normalize data ####
# #################################-
# ## General
# setwd("/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK")
# 
# input_data_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK"
# temp_data_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK/tmp"
# output_data_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK"
# scripts_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK/source2"
# 
# source(sprintf("%s/lib/normalization.R", scripts_folder))
# source(sprintf("%s/lib/aux.R", scripts_folder))
# 
# #The input files are the outputs from the predictors
# output_per_aminoacid_file <- sprintf("%s/output_per_aminoacid.tsv", output_data_folder)
# output_per_protein_file <- sprintf("%s/output_per_protein.tsv", output_data_folder)
# output_per_kmer_file <- sprintf("%s/output_per_kmer.tsv", output_data_folder)
# output_per_peptide_file <- sprintf("%s/output_per_peptide.tsv", output_data_folder)
# output_per_repeat_file <- sprintf("%s/output_per_repeat.tsv", output_data_folder)
# 
# #The output method can be "write" (makes files), "list" (returns a list), or "both"
# output_method <- "write"
# normalized_output_per_peptide_file <- sprintf("%s/normalized_output_per_peptide.tsv", output_data_folder)
# 
# peptide_length <- 15
# peptide_overlap <- 14
# 
# ## Paircoil2
# #This is to set to 0 or 1 the value of each aa and then calculate an average for the peptide
# Paircoil2_threshold <- 0.5
# 
# ## Coendemicity
# # Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty is the min amount of times a given kmer has to appear in the
# # coendemic proteome to apply the penalty to peptides containing it
# Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty <- 0
# # Coendemicity_peptide_start_penalty_proportion is the starting proportion to apply a penalty to the score
# # This penalty goes up until the proportion reaches the Coendemicity_protein_max_penalty_proportion, where the penalty
# # becomes 1 and the score 0
# Coendemicity_peptide_start_penalty_proportion <- 0
# Coendemicity_peptide_max_penalty_proportion <- 2/15
# 
# normalizePeptideData_fromFile(output_per_aminoacid_file = output_per_aminoacid_file, output_per_protein_file = output_per_protein_file, output_per_kmer_file = output_per_kmer_file, output_per_peptide_file = output_per_peptide_file, output_per_repeat_file = output_per_repeat_file,
#                               output_method = output_method, normalized_output_per_peptide_file = normalized_output_per_peptide_file,
#                               peptide_length = peptide_length, peptide_overlap = peptide_overlap,
#                               Paircoil2_threshold = Paircoil2_threshold,
#                               Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty = Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_peptide_start_penalty_proportion = Coendemicity_peptide_start_penalty_proportion, Coendemicity_peptide_max_penalty_proportion = Coendemicity_peptide_max_penalty_proportion,
#                               use_BepiPred = use_BepiPred, use_Paircoil2 = use_Paircoil2, use_PredGPI = use_PredGPI, use_SignalP = use_SignalP,
#                               use_TMHMM = use_TMHMM, use_NetSurfp = use_NetSurfp,
#                               use_Iupred = use_Iupred, use_NetOglyc = use_NetOglyc, use_Xstream = use_Xstream, use_NetMHCIIpan = use_NetMHCIIpan,
#                               use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight = use_MolecularWeight,
#                               use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
