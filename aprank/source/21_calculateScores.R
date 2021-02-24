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
        columns_used <- c("is_antigen")
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
        protein_model <- glm(is_antigen ~ ., family = binomial(link = 'logit'), data = model_data)
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
    
    load(protein_model_file) # > model
    
    calculateProteinScores(normalized_output_per_protein = normalized_output_per_protein,
                           protein_model = model,
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
        columns_used <- c("is_antigen")
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
        peptide_model <- glm(is_antigen ~ ., family = binomial(link = 'logit'), data = model_data)
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
    
    load(peptide_model_file) # > model
    
    score_output_per_protein <- fread(score_output_per_protein_file, header = T, sep = "\t", na.strings = NULL)
    
    calculatePeptideScores(normalized_output_per_peptide = normalized_output_per_peptide,
                           peptide_model = model,
                           score_output_per_protein = score_output_per_protein,
                           output_method = output_method, score_output_per_peptide_file =  score_output_per_peptide_file,
                           use_BepiPred = use_BepiPred, use_Paircoil2 =  use_Paircoil2, use_PredGPI =  use_PredGPI, use_SignalP =  use_SignalP,
                           use_TMHMM = use_TMHMM, use_NetSurfp =  use_NetSurfp,
                           use_Iupred = use_Iupred, use_NetOglyc =  use_NetOglyc, use_Xstream =  use_Xstream, use_NetMHCIIpan =  use_NetMHCIIpan,
                           use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
}

##################-
#### **CALL** ####
##################-
####################################-
#### Config - Select Predictors ####
####################################-
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
# #### Config - Protein scores ####
# #################################-
# ## General
# setwd("/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK")
# 
# model_data_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK/models"
# output_data_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK"
# scripts_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK/source2"
# 
# #The input files are the outputs from the normalizations
# normalized_output_per_protein_file <- sprintf("%s/normalized_output_per_protein.tsv", output_data_folder)
# 
# #The linear model file
# protein_model_file <- sprintf("%s/generic_protein_model.rda", model_data_folder)
# 
# #The output method can be "write" (makes files), "list" (returns a list), or "both"
# protein_scores_output_method <- "write"
# score_output_per_protein_file <- sprintf("%s/score_output_per_protein.tsv", output_data_folder)
# 
# #################################-
# #### Config - Peptide scores ####
# #################################-
# #The input files are the outputs from the normalizations
# normalized_output_per_peptide_file <- sprintf("%s/normalized_output_per_peptide.tsv", output_data_folder)
# 
# #The linear model file
# peptide_model_file <- sprintf("%s/generic_peptide_model.rda", model_data_folder)
# 
# #The output method can be "write" (makes files), "list" (returns a list), or "both"
# peptide_scores_output_method <- "write"
# score_output_per_peptide_file <- sprintf("%s/score_output_per_peptide.tsv", output_data_folder)
# 
# calculateProteinScores_fromFile(normalized_output_per_protein_file = normalized_output_per_protein_file,
#                                 protein_model_file = protein_model_file,
#                                 output_method = protein_scores_output_method, score_output_per_protein_file =  score_output_per_protein_file,
#                                 use_BepiPred = use_BepiPred, use_Paircoil2 =  use_Paircoil2, use_PredGPI =  use_PredGPI, use_SignalP =  use_SignalP,
#                                 use_TMHMM = use_TMHMM, use_NetSurfp =  use_NetSurfp,
#                                 use_Iupred = use_Iupred, use_NetOglyc =  use_NetOglyc, use_Xstream =  use_Xstream, use_NetMHCIIpan =  use_NetMHCIIpan,
#                                 use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight =  use_MolecularWeight,
#                                 use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
# 
# calculatePeptideScores_fromFile(normalized_output_per_peptide_file = normalized_output_per_peptide_file,
#                                 peptide_model_file = peptide_model_file,
#                                 score_output_per_protein_file = score_output_per_protein_file,
#                                 output_method = peptide_scores_output_method, score_output_per_peptide_file =  score_output_per_peptide_file,
#                                 use_BepiPred = use_BepiPred, use_Paircoil2 =  use_Paircoil2, use_PredGPI =  use_PredGPI, use_SignalP =  use_SignalP,
#                                 use_TMHMM = use_TMHMM, use_NetSurfp =  use_NetSurfp,
#                                 use_Iupred = use_Iupred, use_NetOglyc =  use_NetOglyc, use_Xstream =  use_Xstream, use_NetMHCIIpan =  use_NetMHCIIpan,
#                                 use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
