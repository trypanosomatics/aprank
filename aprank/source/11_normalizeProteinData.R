# BSD 2-Clause License
# 
# Copyright (c) 2021, Alejandro Ricci (aricci@iib.unsam.edu.ar), Fernán Agüero (fernan@iib.unsam.edu.ar)
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#   
# 1. Redistributions of source code must retain the above copyright notice, this
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
##################-
#### **MAIN** ####
##################-
normalizeProteinData <- function(output_per_aminoacid, output_per_protein, output_per_kmer, output_per_peptide, output_per_repeat,
                                 output_method, normalized_output_per_protein_file = "",
                                 Paircoil2_fragment_length = 50, Paircoil2_threshold = 0.5,
                                 Xstream_min_period = 1, Xstream_min_copy_number = 1, Xstream_max_consensus_error = 1,
                                 Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty = 0, Coendemicity_protein_start_penalty_proportion = 0, Coendemicity_protein_max_penalty_proportion = 1/3,
                                 use_BepiPred = 1, use_Paircoil2 = 1, use_PredGPI = 1, use_SignalP = 1,
                                 use_TMHMM = 1, use_NetSurfp = 1,
                                 use_Iupred = 1, use_NetOglyc = 1, use_Xstream = 1, use_NetMHCIIpan = 1,
                                 use_IsoelectricPoint = 1, use_MolecularWeight = 1,
                                 use_SelfSimilarity = 1, use_CrossReactivity = 1, use_Coendemicity = 1) {
    ####################-
    #### **CONFIG** ####
    ####################-
    round_decimals <- 4
    
    ################################################-
    #### **CALCULATE & NORMALIZE PROTEIN DATA** ####
    ################################################-
    writeLines("Normalizing protein data...")
    
    normalized_output_per_protein <- output_per_protein[, .(id, original_id)]

    ##################-
    #### BepiPred ####
    ##################-
    if (use_BepiPred) {
        ## Calculated as: mean of the aa scores
        ## Normalization: the distribution of scores for the aa range between -2 and 2, so I'll normalize it between -1.5 and 1.5
        BepiPred_aux <- output_per_aminoacid[, mean(BepiPred), by = .(id)]
        normalized_output_per_protein$BepiPred <- normalizeLinearFixed(BepiPred_aux$V1, -1.5, 1.5)
        normalized_output_per_protein[, BepiPred := round(BepiPred, round_decimals)]
    } else {
        normalized_output_per_protein$BepiPred <- 0
    }
    
    ##############################################-
    #### Isoelectric point & Molecular weight ####
    ##############################################-
    #Isoelectronic Point
    if (use_IsoelectricPoint) {
        ##Calculated as: PI 0 = 0, PI 14 = 1
        normalized_output_per_protein$IsoelectricPoint <- output_per_protein$IsoelectricPoint / 14
        normalized_output_per_protein[, IsoelectricPoint := round(IsoelectricPoint, round_decimals)]
    } else {
        normalized_output_per_protein$IsoelectricPoint <- 0
    }
    
    #Molecular Weight
    if (use_MolecularWeight) {
        #### Option 1 (the median is 0.5)
        # normalized_proteins$MolecularWeight <- normalizationSigmoid05(proteins$MolecularWeight, median(proteins$MolecularWeight)) #median = 0.5
        
        #### Option 2
        ## Calculated as: 30.000 is 0.5 (which is the median if I look at a lot of data)
        normalized_output_per_protein$MolecularWeight <- normalizeSigmoid05(output_per_protein$MolecularWeight, 30000) #30000 = 0.5
        normalized_output_per_protein[, MolecularWeight := round(MolecularWeight, round_decimals)]
    } else {
        normalized_output_per_protein$MolecularWeight <- 0
    }
    
    ################-
    #### Iupred ####
    ################-
    ##Iupred
    if (use_Iupred) {
        #### Option 1 (mean) (better AUC)
        # Iupred_aux <- aminoacids[, mean(Iupred), by = .(protein)]
        # normalized_proteins$Iupred <- Iupred_aux$V1
        
        #### Option 2
        ## Calculated as: proportion of values greater or equal than 0.5
        Iupred_aux <- output_per_aminoacid[, sum(Iupred >= 0.5) / .N, by = .(id)]
        normalized_output_per_protein$Iupred <- Iupred_aux$V1    
        normalized_output_per_protein[, Iupred := round(Iupred, round_decimals)]
    } else {
        normalized_output_per_protein$Iupred <- 0
    }
    
    #####################-
    #### NetMHCIIpan ####
    #####################-
    if (use_NetMHCIIpan) {
        #### Option 1 (proportion of fragments below threshold)
        # NetMHCIIpan_peptides$rank_score <- as.numeric(NetMHCIIpan_peptides$rank <= NetMHCIIpan_rank_threshold)
        # NetMHCIIpan_peptides <- NetMHCIIpan_peptides[, .(rank_score = mean(rank_score)), by = .(protein)]
        # normalized_proteins$NetMHCIIpan <- NetMHCIIpan_peptides$rank_score
        
        #### Option 2 (amount of fragments below threshold)
        # NetMHCIIpan_peptides$rank_score <- as.numeric(NetMHCIIpan_peptides$rank <= NetMHCIIpan_rank_threshold)
        # NetMHCIIpan_peptides <- NetMHCIIpan_peptides[, .(rank_score = sum(rank_score)), by = .(protein)]
        # normalized_proteins$NetMHCIIpan <- normalizeLinear(NetMHCIIpan_peptides$rank_score, 0, 150) #check range
        
        #### Option 3
        ## Calculated as: 1 - rank mean (best AUC)
        NetMHCIIpan_aux <- output_per_peptide[, .(NetMHCIIpan = mean(mean_rank)), by = .(id)]
        normalized_output_per_protein <- merge(normalized_output_per_protein,
                                               NetMHCIIpan_aux,
                                               by = c("id"),
                                               all.x = T)
        normalized_output_per_protein[, NetMHCIIpan := 1 - normalizeLinearFixed((NetMHCIIpan / 100), 0.05, 0.5)]
        normalized_output_per_protein[, NetMHCIIpan := round(NetMHCIIpan, round_decimals)]
        
        #### Option 4 (1 - best rank mean between the two alleles)
        # NetMHCIIpan_peptides_aux <- NetMHCIIpan_peptides[, .(rank_score = mean(rank)), by = .(protein, MHCII_allele)]
        # NetMHCIIpan_peptides_aux <- NetMHCIIpan_peptides_aux[, .(rank_score = min(rank_score)), by = protein]
        # normalized_proteins$NetMHCIIpan <- 1 - (NetMHCIIpan_peptides_aux$rank_score / 100)        
    } else {
        normalized_output_per_protein$NetMHCIIpan <- 0
    }
    
    ##################-
    #### NetOglyc ####
    ##################-
    if (use_NetOglyc) {
        #### Option 1
        ## Calculated as: proportion of glycosilated residues
        NetOglyc_aux <- output_per_aminoacid[, mean(NetOglyc), by = .(id)]
        normalized_output_per_protein$NetOglyc <- normalizeLinearFixed(NetOglyc_aux$V1, 0, 0.05)
        normalized_output_per_protein[, NetOglyc := round(NetOglyc, round_decimals)]
        
        #### Option 2 (has at least 1 glycosilated residue)
        # NetOglyc_aux <- aminoacids[, as.numeric(any(NetOglyc == 1)), by = .(protein)]
        # normalized_proteins$NetOglyc <- NetOglyc_aux$V1        
    } else {
        normalized_output_per_protein$NetOglyc <- 0
    }
    
    ##################-
    #### NetSurfp ####
    ##################-
    if (use_NetSurfp) {
        ##NetSurfp RSA
        NetSurfp_RSA_aux <- output_per_aminoacid[, mean(NetSurfp_RSA), by = .(id)]
        normalized_output_per_protein$NetSurfp_RSA <- NetSurfp_RSA_aux$V1
        normalized_output_per_protein[, NetSurfp_RSA := round(NetSurfp_RSA, round_decimals)]
        
        ##NetSurfp AH and BS
        #### Option 1
        NetSurfp_AH_aux <- output_per_aminoacid[, mean(NetSurfp_AH), by = .(id)]
        NetSurfp_BS_aux <- output_per_aminoacid[, mean(NetSurfp_BS), by = .(id)]
        normalized_output_per_protein$NetSurfp_AH <- NetSurfp_AH_aux$V1
        normalized_output_per_protein$NetSurfp_BS <- NetSurfp_BS_aux$V1
        normalized_output_per_protein[, NetSurfp_AH := round(NetSurfp_AH, round_decimals)]
        normalized_output_per_protein[, NetSurfp_BS := round(NetSurfp_BS, round_decimals)]
        
        #### Option 2 (for some reason AH * BS gives a better AUC than either of those two alone)
        # NetSurfp_AH_aux <- aminoacids[, mean(NetSurfp_AH), by = .(id)]
        # NetSurfp_BS_aux <- aminoacids[, mean(NetSurfp_BS), by = .(id)]
        # normalized_proteins$NetSurfp_AH_BS <- NetSurfp_AH_aux$V1 * NetSurfp_BS_aux$V1
    } else {
        normalized_output_per_protein$NetSurfp_RSA <- 0
        normalized_output_per_protein$NetSurfp_AH <- 0
        normalized_output_per_protein$NetSurfp_BS <- 0
    }
    
    ###################-
    #### Paircoil2 ####
    ###################-
    if (use_Paircoil2) {
        #### Option 1
        # ## Calculated as: mean of scores (has a bit better AUC than proportion above threshold, but less meaning)
        # Paircoil2_aux <- aminoacids[, mean(Paircoil2), by = .(protein)]
        # normalized_proteins$Paircoil2 <- Paircoil2_aux$V1
        
        #### Option 2
        ## Calculated as: check if there is a fragment of a given length that has all the values above a given threshold
        Paircoil2_aux <- output_per_aminoacid[, as.integer(getSizeOfLongestSegmentAboveThreshold(Paircoil2, Paircoil2_threshold) >= Paircoil2_fragment_length), by = .(id)]
        normalized_output_per_protein$Paircoil2 <- Paircoil2_aux$V1
        
        #### Option 3
        # ## Calculated as: proportion of fragments above threshold
        # Paircoil2_aux <- aminoacids[, sum(Paircoil2 > Paircoil2_threshold) / .N, by = .(protein)]
        # normalized_proteins$Paircoil2 <- Paircoil2_aux$V1
    } else {
        normalized_output_per_protein$Paircoil2 <- 0
    }
    
    #################-
    #### PredGPI ####
    #################-
    if (use_PredGPI) {
        #### Option 1 (the score obtained)
        # PredGPI_aux <- proteins$PredGPI_score
        # normalized_proteins$PredGPI <- PredGPI_aux
        
        #### Option 2 (if it has a GPI)
        normalized_output_per_protein$PredGPI <- as.integer(output_per_protein$PredGPI_start > 0)
    } else {
        normalized_output_per_protein$PredGPI <- 0
    }
    
    #################-
    #### SignalP ####
    #################-
    if (use_SignalP) {
        #### Option 1 (if it has a SignalP)
        # SignalP_aux <- as.integer(proteins$SignalP_end > 0)
        # normalized_proteins$SignalP <- SignalP_aux
        
        #### Option 2
        ## Calculated as: mean of the scores (most of the score comes from the S part, but the C part has good score alone too)
        ### (I remove the 0 because those were the aa outside the SignalP)
        ## Normalization: I'll normalize it between 0 and 0.05 (after multipling the result is a small number)
        SignalP_S_aux <- output_per_aminoacid[SignalP_S > 0, mean(SignalP_S), by = .(id)]
        SignalP_C_aux <- output_per_aminoacid[SignalP_C > 0, mean(SignalP_C), by = .(id)]
        normalized_output_per_protein$SignalP <- SignalP_S_aux$V1 * SignalP_C_aux$V1
        normalized_output_per_protein[, SignalP := normalizeLinearFixed(SignalP, 0, 0.05)]
        normalized_output_per_protein[, SignalP := round(SignalP, round_decimals)]
    } else {
        normalized_output_per_protein$SignalP <- 0
    }
    
    ###############-
    #### TMHMM ####
    ###############-
    ##TMHMM (AUC near 0.5 at protein level, but I can't do much about it)
    if (use_TMHMM) {
        normalized_output_per_protein$TMHMM <- output_per_protein$TMHMM
    } else {
        normalized_output_per_protein$TMHMM <- 0
    }
    
    #################-
    #### XStream ####
    #################-
    if (use_Xstream) {
        ## Calculated as: filter the repeats, find the biggest copy number for each protein, and apply a sigmoid normalization to it
        if (output_per_repeat[, .N] > 0) {
            Xstream_aux <- output_per_repeat[(period >= Xstream_min_period) & (copy_number >= Xstream_min_copy_number) & (consensus_error <= Xstream_max_consensus_error)]
            Xstream_aux <- Xstream_aux[, .(Xstream = max(copy_number)), by = .(id)]
            normalized_output_per_protein <- merge(normalized_output_per_protein,
                                                   Xstream_aux,
                                                   by = c("id"),
                                                   all.x = T)
            normalized_output_per_protein[is.na(Xstream), Xstream := 0]
            normalized_output_per_protein$Xstream <- normalizeSigmoid09(normalized_output_per_protein$Xstream, 5) #5 repeats = 0.9
            normalized_output_per_protein[, Xstream := round(Xstream, round_decimals)]
        } else {
            normalized_output_per_protein$Xstream <- 0
        }
    } else {
        normalized_proteins$Xstream <- 0
    }
    
    #########################-
    #### Self Similarity ####
    #########################-
    if (use_SelfSimilarity) {
        ## Calculated as: each kmer inside a protein it's assigned a score between 0 or 1 relative to how many other times it appears in the proteome
        ### no other time it's 0
        ### 1 other time it's 0.5
        ### many other times it's close to 1
        ### The scores of all this kmers are added up for every protein and then it's divided by the total amount of kmers
        ### inside that protein (which it's the maximum value that sum can achieve)
        SelfSimilarity_aux <- output_per_kmer[, sum(quantity * normalizeSigmoid05(quantity_in_proteome - 1, 1)) / sum(quantity), by = .(id)]
        normalized_output_per_protein$SelfSimilarity <- SelfSimilarity_aux$V1            
        normalized_output_per_protein[, SelfSimilarity := round(SelfSimilarity, round_decimals)]
    } else {
        normalized_output_per_protein$SelfSimilarity <- 0
    }
    
    #########################-
    #### CrossReactivity ####
    #########################-
    if (use_CrossReactivity) {
        ## Calculated as: each kmer inside a protein it's assigned a score between 0 or 1 relative to how many times it appears in the host proteome
        ### no other time it's 0
        ### 1 other time it's 0.5
        ### many other times it's close to 1
        ### The scores of all this kmers are added up for every protein and then it's divided by the total amount of kmers
        ### inside that protein (which it's the maximum value that sum can achieve)
        CrossReactivity_aux <- output_per_kmer[, sum(quantity * normalizeSigmoid05(quantity_in_host_proteome, 1)) / sum(quantity), by = .(id)]
        normalized_output_per_protein$CrossReactivity <- CrossReactivity_aux$V1    
        normalized_output_per_protein[, CrossReactivity := round(CrossReactivity, round_decimals)]
    } else {
        normalized_output_per_protein$CrossReactivity <- 0    
    }
    
    ######################-
    #### Coendemicity ####
    ######################-
    if (use_Coendemicity) {
        ## Coendemicity is a penalty to apply to the future scores
        ## Calculated as: Find the kmers in my proteome that appear at least a given amount of times in the coendemic proteome
        ### For those kmers, calculate their amount in each of the proteins, and then calculate that as a proportion based on the total of kmers in that protein
        ## Normalized as: Said proportion is then normalized linearly between two values of min and max penalty
        output_per_kmer_aux <- output_per_kmer[, .(kmer, id, quantity, quantity_in_coendemic_proteome)]
        output_per_kmer_aux$coendemic_quantity_aux <- 0L
        output_per_kmer_aux[(quantity_in_coendemic_proteome >= Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty) &
                                (quantity_in_coendemic_proteome > 0), coendemic_quantity_aux := quantity]
        output_per_kmer_aux <- output_per_kmer_aux[, .(kmer_amount = sum(quantity),
                                                       coendemic_kmer_amount = sum(coendemic_quantity_aux)),
                                                   by = .(id)]
        output_per_kmer_aux[, coendemic_proportion := coendemic_kmer_amount / kmer_amount]
        output_per_kmer_aux[, Coendemicity_Penalty := normalizeLinearFixed(coendemic_proportion,
                                                                           Coendemicity_protein_start_penalty_proportion,
                                                                           Coendemicity_protein_max_penalty_proportion)]
        
        normalized_output_per_protein <- merge(normalized_output_per_protein,
                                               output_per_kmer_aux[, .(id, Coendemicity_Penalty)],
                                               by = "id")
        normalized_output_per_protein[, Coendemicity_Penalty := round(Coendemicity_Penalty, round_decimals)]
    } else {
        normalized_output_per_protein$Coendemicity_Penalty <- 0
    }
    
    #########################-
    #### **OUTPUT DATA** ####
    #########################-
    writeLines("Readying Output...")
    
    if (output_method %in% c("write", "both")) {
        write.table(normalized_output_per_protein, file = normalized_output_per_protein_file, col.names = T, row.names = F, sep = "\t", quote = T)
    }
    
    list_output <- list()
    if (output_method %in% c("list", "both")) {
        #Create the output
        list_output[["normalized_output_per_protein"]] <- normalized_output_per_protein
    }
    
    writeLines("Done!")
    
    list_output
}
normalizeProteinData_fromFile <- function(output_per_aminoacid_file, output_per_protein_file, output_per_kmer_file, output_per_peptide_file, output_per_repeat_file,
                                          output_method, normalized_output_per_protein_file = "",
                                          Paircoil2_fragment_length = 50, Paircoil2_threshold = 0.5,
                                          Xstream_min_period = 1, Xstream_min_copy_number = 1, Xstream_max_consensus_error = 1,
                                          Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty = 0, Coendemicity_protein_start_penalty_proportion = 0, Coendemicity_protein_max_penalty_proportion = 1/3,
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
    
    normalizeProteinData(output_per_aminoacid = output_per_aminoacid, output_per_protein = output_per_protein, output_per_kmer = output_per_kmer, output_per_peptide = output_per_peptide, output_per_repeat = output_per_repeat,
                         output_method = output_method, normalized_output_per_protein_file = normalized_output_per_protein_file,
                         Paircoil2_fragment_length = Paircoil2_fragment_length, Paircoil2_threshold = Paircoil2_threshold,
                         Xstream_min_period = Xstream_min_period, Xstream_min_copy_number = Xstream_min_copy_number, Xstream_max_consensus_error = Xstream_max_consensus_error,
                         Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty = Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_protein_start_penalty_proportion = Coendemicity_protein_start_penalty_proportion, Coendemicity_protein_max_penalty_proportion = Coendemicity_protein_max_penalty_proportion,
                         use_BepiPred = use_BepiPred, use_Paircoil2 = use_Paircoil2, use_PredGPI = use_PredGPI, use_SignalP = use_SignalP,
                         use_TMHMM = use_TMHMM, use_NetSurfp = use_NetSurfp,
                         use_Iupred = use_Iupred, use_NetOglyc = use_NetOglyc, use_Xstream = use_Xstream, use_NetMHCIIpan = use_NetMHCIIpan,
                         use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight = use_MolecularWeight,
                         use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
}
