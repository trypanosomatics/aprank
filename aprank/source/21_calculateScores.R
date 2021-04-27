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

##################-
#### **MAIN** ####
##################-
calculateProteinScores <- function(normalized_output_per_protein,
                                   protein_model,
                                   output_method, score_output_per_protein_file = "",
                                   use_BepiPred = 1, use_Paircoil2 = 1, use_PredGPI = 1, use_SignalP = 1,
                                   use_TMHMM = 1, use_NetSurfp = 1,
                                   use_Iupred = 1, use_NetOglyc = 1, use_Xstream = 1, use_NetMHCIIpan = 1,
                                   use_IsoelectricPoint = 1, use_MolecularWeight = 1,
                                   use_SelfSimilarity = 1, use_CrossReactivity = 1, use_Coendemicity = 1) {
    ####################-
    #### **CONFIG** ####
    ####################-
    round_decimals <- 4
    
    ##############################-
    #### **CALCULATE SCORES** ####
    ##############################-
    writeLines("Calculating protein scores...")
    
    score_output_per_protein <- normalized_output_per_protein[, .(id, original_id, Coendemicity_Penalty)]
    
    #######################################-
    #### Modify the model if necessary ####
    #######################################-
    #Check if I have to modify the model
    if ((use_BepiPred == 0) | (use_Paircoil2 == 0) | (use_PredGPI == 0) | (use_SignalP == 0) |
        (use_TMHMM == 0) | (use_NetSurfp == 0) |
        (use_Iupred == 0) | (use_NetOglyc == 0) | (use_Xstream == 0) | (use_NetMHCIIpan == 0) |
        (use_IsoelectricPoint == 0) | (use_MolecularWeight == 0) |
        (use_SelfSimilarity == 0) | (use_CrossReactivity == 0)) {
        
        #Parse the data used to make that model
        model_data <- as.data.table(protein_model$model)
        
        #Select which columns to use in the new prediction
        columns_used <- c("in_antigenic_cluster")
        if (use_BepiPred == 1) {columns_used <- c(columns_used, "BepiPred")}
        if (use_Paircoil2 == 1) {columns_used <- c(columns_used, "Paircoil2")}
        if (use_PredGPI == 1) {columns_used <- c(columns_used, "PredGPI")}
        if (use_SignalP == 1) {columns_used <- c(columns_used, "SignalP")}
        if (use_TMHMM == 1) {columns_used <- c(columns_used, "TMHMM")}
        if (use_NetSurfp == 1) {columns_used <- c(columns_used, "NetSurfp_RSA", "NetSurfp_AH", "NetSurfp_BS")}
        if (use_Iupred == 1) {columns_used <- c(columns_used, "Iupred")}
        if (use_NetOglyc == 1) {columns_used <- c(columns_used, "NetOglyc")}
        if (use_Xstream == 1) {columns_used <- c(columns_used, "Xstream")}
        if (use_NetMHCIIpan == 1) {columns_used <- c(columns_used, "NetMHCIIpan")}
        if (use_IsoelectricPoint == 1) {columns_used <- c(columns_used, "IsoelectricPoint")}
        if (use_MolecularWeight == 1) {columns_used <- c(columns_used, "MolecularWeight")}
        if (use_SelfSimilarity == 1) {columns_used <- c(columns_used, "SelfSimilarity")}
        if (use_CrossReactivity == 1) {columns_used <- c(columns_used, "CrossReactivity")}
        
        #Keep only the columns I want
        model_data <- model_data[, columns_used, with=FALSE]
        
        #Recalculate the model that uses only this new columns
        protein_model <- glm(in_antigenic_cluster ~ ., family = binomial(link = 'logit'), data = model_data)
    }
    
    #########################-
    #### Apply the model ####
    #########################-
    score_output_per_protein$protein_score <- predict(protein_model, newdata = normalized_output_per_protein, type = "response")
    score_output_per_protein[, protein_score := round(protein_score, round_decimals)]
    score_output_per_protein<- score_output_per_protein[order(-protein_score, original_id)]
    
    ########################################-
    #### Apply the Coendemicity Penalty ####
    ########################################-
    if (use_Coendemicity) {
        score_output_per_protein[, protein_score_wCoendemicPenalty := protein_score * (1 - Coendemicity_Penalty)]
        score_output_per_protein[, protein_score_wCoendemicPenalty := round(protein_score_wCoendemicPenalty, round_decimals)]
        
        score_output_per_protein <- score_output_per_protein[, .(id, original_id, protein_score, protein_score_wCoendemicPenalty)]
    } else {
        score_output_per_protein <- score_output_per_protein[, .(id, original_id, protein_score)]
    }
    
    #########################-
    #### **OUTPUT DATA** ####
    #########################-
    writeLines("Readying Output...")
    
    if (output_method %in% c("write", "both")) {
        write.table(score_output_per_protein, file = score_output_per_protein_file, col.names = T, row.names = F, sep = "\t", quote = T)
    }
    
    list_output <- list()
    if (output_method %in% c("list", "both")) {
        #Create the output
        list_output[["score_output_per_protein"]] <- score_output_per_protein
    }
    
    writeLines("Done!")
    
    list_output
}
calculateProteinScores_fromFile <- function(normalized_output_per_protein_file,
                                            protein_model_file,
                                            output_method, score_output_per_protein_file = "",
                                            use_BepiPred = 1, use_Paircoil2 = 1, use_PredGPI = 1, use_SignalP = 1,
                                            use_TMHMM = 1, use_NetSurfp = 1,
                                            use_Iupred = 1, use_NetOglyc = 1, use_Xstream = 1, use_NetMHCIIpan = 1,
                                            use_IsoelectricPoint = 1, use_MolecularWeight = 1,
                                            use_SelfSimilarity = 1, use_CrossReactivity = 1, use_Coendemicity = 1) {
    normalized_output_per_protein <- fread(normalized_output_per_protein_file, header = T, sep = "\t", na.strings = NULL)
    
    load(protein_model_file) # > protein_balanced_generic_model
    
    calculateProteinScores(normalized_output_per_protein = normalized_output_per_protein,
                           protein_model = protein_balanced_generic_model,
                           output_method = output_method, score_output_per_protein_file =  score_output_per_protein_file,
                           use_BepiPred = use_BepiPred, use_Paircoil2 =  use_Paircoil2, use_PredGPI =  use_PredGPI, use_SignalP =  use_SignalP,
                           use_TMHMM = use_TMHMM, use_NetSurfp =  use_NetSurfp,
                           use_Iupred = use_Iupred, use_NetOglyc =  use_NetOglyc, use_Xstream =  use_Xstream, use_NetMHCIIpan =  use_NetMHCIIpan,
                           use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight =  use_MolecularWeight,
                           use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
}

calculatePeptideScores <- function(normalized_output_per_peptide,
                                   peptide_model,
                                   score_output_per_protein,
                                   output_method, score_output_per_peptide_file = "",
                                   use_BepiPred = 1, use_Paircoil2 = 1, use_PredGPI = 1, use_SignalP = 1,
                                   use_TMHMM = 1, use_NetSurfp = 1,
                                   use_Iupred = 1, use_NetOglyc = 1, use_Xstream = 1, use_NetMHCIIpan = 1,
                                   use_SelfSimilarity = 1, use_CrossReactivity = 1, use_Coendemicity = 1) {
    ####################-
    #### **CONFIG** ####
    ####################-
    round_decimals <- 4
    
    ##############################-
    #### **CALCULATE SCORES** ####
    ##############################-
    writeLines("Calculating peptide scores...")
    
    score_output_per_peptide <- normalized_output_per_peptide[, .(id, start, end, peptide, Coendemicity_Penalty)]
    
    #######################################-
    #### Modify the model if necessary ####
    #######################################-
    #Check if I have to modify the model
    if ((use_BepiPred == 0) | (use_Paircoil2 == 0) | (use_PredGPI == 0) | (use_SignalP == 0) |
        (use_TMHMM == 0) | (use_NetSurfp == 0) |
        (use_Iupred == 0) | (use_NetOglyc == 0) | (use_Xstream == 0) | (use_NetMHCIIpan == 0) |
        (use_SelfSimilarity == 0) | (use_CrossReactivity == 0)) {
        
        #Parse the data used to make that model
        model_data <- as.data.table(peptide_model$model)
        
        #Select which columns to use in the new prediction
        columns_used <- c("is_antigenic")
        if (use_BepiPred == 1) {columns_used <- c(columns_used, "BepiPred")}
        if (use_Paircoil2 == 1) {columns_used <- c(columns_used, "Paircoil2")}
        if (use_PredGPI == 1) {columns_used <- c(columns_used, "PredGPI")}
        if (use_SignalP == 1) {columns_used <- c(columns_used, "SignalP")}
        if (use_TMHMM == 1) {columns_used <- c(columns_used, "TMHMM")}
        if (use_NetSurfp == 1) {columns_used <- c(columns_used, "NetSurfp_RSA", "NetSurfp_AH", "NetSurfp_BS")}
        if (use_Iupred == 1) {columns_used <- c(columns_used, "Iupred")}
        if (use_NetOglyc == 1) {columns_used <- c(columns_used, "NetOglyc")}
        if (use_Xstream == 1) {columns_used <- c(columns_used, "Xstream")}
        if (use_NetMHCIIpan == 1) {columns_used <- c(columns_used, "NetMHCIIpan")}
        # if (use_IsoelectricPoint == 1) {columns_used <- c(columns_used, "IsoelectricPoint")}
        # if (use_MolecularWeight == 1) {columns_used <- c(columns_used, "MolecularWeight")}
        if (use_SelfSimilarity == 1) {columns_used <- c(columns_used, "SelfSimilarity")}
        if (use_CrossReactivity == 1) {columns_used <- c(columns_used, "CrossReactivity")}
        
        #Keep only the columns I want
        model_data <- model_data[, columns_used, with=FALSE]
        
        #Recalculate the model that uses only this new columns
        peptide_model <- glm(is_antigenic ~ ., family = binomial(link = 'logit'), data = model_data)
    }
    
    #########################-
    #### Apply the model ####
    #########################-
    score_output_per_peptide$peptide_score <- predict(peptide_model, newdata = normalized_output_per_peptide, type = "response")
    score_output_per_peptide[, peptide_score := round(peptide_score, round_decimals)]
    
    #Add the protein data
    score_output_per_peptide <- merge(score_output_per_peptide,
                                      score_output_per_protein,
                                      by = "id")
    score_output_per_peptide[, combined_score := (protein_score + peptide_score) / 2]
    score_output_per_peptide[, combined_score := round(combined_score, round_decimals)]
    score_output_per_peptide <- score_output_per_peptide[order(-combined_score, original_id)]
    
    ########################################-
    #### Apply the Coendemicity Penalty ####
    ########################################-
    if (use_Coendemicity) {
        score_output_per_peptide[, peptide_score_wCoendemicPenalty := peptide_score * (1 - Coendemicity_Penalty)]
        score_output_per_peptide[, peptide_score_wCoendemicPenalty := round(peptide_score_wCoendemicPenalty, round_decimals)]
        
        score_output_per_peptide[, combined_score_wCoendemicPenalty := (protein_score_wCoendemicPenalty + peptide_score_wCoendemicPenalty) / 2]
        score_output_per_peptide[, combined_score_wCoendemicPenalty := round(combined_score_wCoendemicPenalty, round_decimals)]
        
        score_output_per_peptide <- score_output_per_peptide[, .(id, original_id, start, end, peptide,
                                                                 peptide_score, protein_score, combined_score,
                                                                 peptide_score_wCoendemicPenalty, protein_score_wCoendemicPenalty, combined_score_wCoendemicPenalty)]
    } else {
        score_output_per_peptide <- score_output_per_peptide[, .(id, original_id, start, end, peptide,
                                                                 peptide_score, protein_score, combined_score)]
    }
    
    #########################-
    #### **OUTPUT DATA** ####
    #########################-
    writeLines("Readying Output...")
    
    if (output_method %in% c("write", "both")) {
        write.table(score_output_per_peptide, file = score_output_per_peptide_file, col.names = T, row.names = F, sep = "\t", quote = T)
    }
    
    list_output <- list()
    if (output_method %in% c("list", "both")) {
        #Create the output
        list_output[["score_output_per_peptide"]] <- score_output_per_peptide
    }
    
    writeLines("Done!")
    
    list_output
}
calculatePeptideScores_fromFile <- function(normalized_output_per_peptide_file,
                                            peptide_model_file,
                                            score_output_per_protein_file,
                                            output_method, score_output_per_peptide_file = "",
                                            use_BepiPred = 1, use_Paircoil2 = 1, use_PredGPI = 1, use_SignalP = 1,
                                            use_TMHMM = 1, use_NetSurfp = 1,
                                            use_Iupred = 1, use_NetOglyc = 1, use_Xstream = 1, use_NetMHCIIpan = 1,
                                            use_SelfSimilarity = 1, use_CrossReactivity = 1, use_Coendemicity = 1) {
    normalized_output_per_peptide <- fread(normalized_output_per_peptide_file, header = T, sep = "\t", na.strings = NULL)
    
    load(peptide_model_file) # > peptide_balanced_generic_model
    
    score_output_per_protein <- fread(score_output_per_protein_file, header = T, sep = "\t", na.strings = NULL)
    
    calculatePeptideScores(normalized_output_per_peptide = normalized_output_per_peptide,
                           peptide_model = peptide_balanced_generic_model,
                           score_output_per_protein = score_output_per_protein,
                           output_method = output_method, score_output_per_peptide_file =  score_output_per_peptide_file,
                           use_BepiPred = use_BepiPred, use_Paircoil2 =  use_Paircoil2, use_PredGPI =  use_PredGPI, use_SignalP =  use_SignalP,
                           use_TMHMM = use_TMHMM, use_NetSurfp =  use_NetSurfp,
                           use_Iupred = use_Iupred, use_NetOglyc =  use_NetOglyc, use_Xstream =  use_Xstream, use_NetMHCIIpan =  use_NetMHCIIpan,
                           use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
}
