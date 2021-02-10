####################-
#### **CONFIG** ####
####################-
aprank_development_folder <- "/home/user/Desktop/APRANK/ModelDevelopment"
aprank_folder <- "/home/user/Desktop/APRANK"

setwd(aprank_folder)

library(data.table)
library(ROSE)
library(pROC)

source("source/lib/normalization.R")

###########################################################-
#### Config - Select Organisms to include in the model ####
###########################################################-
organisms <- c("Borrelia_burgdorferi",
               "Brucella_melitensis",
               "Coxiella_burnetii",
               "Escherichia_coli_H10407",
               "Francisella_tularensis",
               "Leptospira_interrogans",
               "Porphyromonas_gingivalis_837",
               "Mycobacterium_leprae",
               "Mycobacterium_tuberculosis",
               "Staphylococcus_aureus",
               "Streptococcus_pyogenes_serotype_M1",
               "Leishmania_brasilensis",
               "Plasmodium_falciparum",
               "Toxoplasma_gondi",
               "Trypanosoma_cruzi")

##########################-
#### Config - General ####
##########################-
input_data_folder <- sprintf("%s/21_outputs", aprank_development_folder)
output_data_folder <- sprintf("%s/21_outputs", aprank_development_folder)
scripts_folder <- sprintf("%s/source", aprank_folder)

protein_model_output_file <- "protein_balanced_generic_model.rda"
peptide_model_output_file <- "peptide_balanced_generic_model.rda"

protein_scores_file_name <- "scores_per_protein.tsv"
peptide_scores_file_name <- "scores_per_peptide.tsv"

main_seed <- 1234

###############################-
#### Config - Balance Data ####
###############################-
ROSE_protein_amount_per_organism <- 3000
ROSE_peptide_amount_per_organism <- 100000
ROSE_seed <- 42

##################-
#### **MAIN** ####
##################-
set.seed(main_seed)

###############################-
#### **CREATE THE MODELS** ####
###############################-
#For RAM space reasons I'm going to be reading each file once to make the models and then a second time
#(one at a time) to calculate the scores
peptide_file_initializated <- 0
for (organism_for in organisms) {
    # organism_for <- organisms[1]
    
    ###########################-
    #### Read protein data ####
    ###########################-
    #Protein data
    normalized_protein_data_file <- sprintf("%s/%s/normalized_output_per_protein.tsv", input_data_folder, organism_for)
    normalized_protein_data_for <- fread(normalized_protein_data_file, header = T, sep = "\t", na.strings = NULL)
    
    normalized_protein_data_for$organism <- organism_for
    
    ##Add the antigenic tag data
    protein_antigenic_tag_data_file <- sprintf("%s/%s/protein_antigenic_tag.tsv", input_data_folder, organism_for)
    protein_antigenic_tag_data <- fread(protein_antigenic_tag_data_file, header = T, sep = "\t", na.strings = NULL)
    
    normalized_protein_data_for <- merge(normalized_protein_data_for,
                                         protein_antigenic_tag_data,
                                         by = "id")

    ##############################-
    #### Balance protein data ####
    ##############################-
    #Extract the important columns for training
    training_set <- normalized_protein_data_for[, c("BepiPred", "IsoelectricPoint", "MolecularWeight", "Iupred", "NetMHCIIpan", "NetOglyc", 
                                                    "NetSurfp_RSA", "NetSurfp_AH", "NetSurfp_BS", "Paircoil2", "PredGPI", 
                                                    "SignalP", "TMHMM", "Xstream", "SelfSimilarity", "CrossReactivity",
                                                    "in_antigenic_cluster")]
    
    #Balance the data
    balanced_training_set <- ROSE(in_antigenic_cluster ~ ., data = training_set,
                                  N = ROSE_protein_amount_per_organism,
                                  seed = ROSE_seed)
    balanced_training_set <- as.data.table(balanced_training_set$data)
    
    #Add the organism data
    balanced_training_set$organism <- organism_for
    
    #Save the data
    if (organism_for == organisms[1]) {
        full_balanced_training_set <- balanced_training_set
    } else {
        full_balanced_training_set <- rbindlist(list(full_balanced_training_set, balanced_training_set))
    }
    
    ###########################-
    #### Read peptide data ####
    ###########################-
    #Peptide data
    peptide_antigenic_tag_data_file <- sprintf("%s/%s/peptide_antigenic_tag.tsv", input_data_folder, organism_for)
    if (file.exists(peptide_antigenic_tag_data_file)) {
        normalized_peptide_data_file <- sprintf("%s/%s/normalized_output_per_peptide.tsv", input_data_folder, organism_for)
        normalized_peptide_data_for <- fread(normalized_peptide_data_file, header = T, sep = "\t", na.strings = NULL)
        
        normalized_peptide_data_for <- merge(normalized_peptide_data_for,
                                             normalized_protein_data_for[, .(id, original_id)],
                                             by = "id",
                                             all.x = T)
        
        normalized_peptide_data_for$organism <- organism_for
        
        ##Add the antigenic tag data
        peptide_antigenic_tag_data <- fread(peptide_antigenic_tag_data_file, header = T, sep = "\t", na.strings = NULL)
        
        normalized_peptide_data_for <- merge(normalized_peptide_data_for,
                                             peptide_antigenic_tag_data,
                                             by = c("id", "start"))
        
        ##############################-
        #### Balance peptide data ####
        ##############################-
        #Extract the important columns for training
        peptide_training_set <- normalized_peptide_data_for[, c("BepiPred", "Paircoil2", "PredGPI", "SignalP",
                                                                "TMHMM", "NetSurfp_RSA", "NetSurfp_AH", "NetSurfp_BS", 
                                                                "Iupred", "NetOglyc", "Xstream", "NetMHCIIpan", 
                                                                "SelfSimilarity", "CrossReactivity",
                                                                "is_antigenic")]
        
        #Balance the data
        balanced_peptide_training_set <- ROSE(is_antigenic ~ ., data = peptide_training_set,
                                              N = ROSE_peptide_amount_per_organism,
                                              seed = ROSE_seed)
        balanced_peptide_training_set <- as.data.table(balanced_peptide_training_set$data)
        
        #Add the organism data
        balanced_peptide_training_set$organism <- organism_for
        
        #Save the data
        if (peptide_file_initializated == 0) {
            full_balanced_peptide_training_set <- balanced_peptide_training_set
            
            peptide_file_initializated <- 1
        } else {
            full_balanced_peptide_training_set <- rbindlist(list(full_balanced_peptide_training_set, balanced_peptide_training_set))
        }
    }
}
rm(normalized_protein_data_for)
rm(normalized_peptide_data_for)
rm(protein_antigenic_tag_data)
rm(peptide_antigenic_tag_data)
gc()

#Create the generic models
protein_balanced_generic_model <- glm(in_antigenic_cluster ~ ., family = binomial(link='logit'), data = full_balanced_training_set[, -c("organism")])
peptide_balanced_generic_model <- glm(is_antigenic ~ ., family = binomial(link='logit'), data = full_balanced_peptide_training_set[, -c("organism")])

#Save the generic models
save(protein_balanced_generic_model, file = protein_model_output_file)
save(peptide_balanced_generic_model, file = peptide_model_output_file)

#Create the generic noBepipred models
protein_balanced_generic_noBepipred_model <- glm(in_antigenic_cluster ~ ., family = binomial(link='logit'), data = full_balanced_training_set[, -c("organism", "BepiPred")])
peptide_balanced_generic_noBepipred_model <- glm(is_antigenic ~ ., family = binomial(link='logit'), data = full_balanced_peptide_training_set[, -c("organism", "BepiPred")])

##################################-
#### **CALCULATE THE SCORES** ####
##################################-
for (organism_for in organisms) {
    # organism_for <- organisms_to_calculate_scores[1]
    ######################-
    #### **PROTEINS** ####
    ######################-
    #############################################-
    #### Read specific organism protein data ####
    #############################################-
    #Protein data
    normalized_protein_data_file <- sprintf("%s/%s/normalized_output_per_protein.tsv", input_data_folder, organism_for)
    normalized_protein_data_for <- fread(normalized_protein_data_file, header = T, sep = "\t", na.strings = NULL)
    
    normalized_protein_data_for$organism <- organism_for
    
    #Add Bepipred data
    Bepipred_protein_data_file <- sprintf("%s/%s/Bepipred_raw_output_per_protein.tsv", input_data_folder, organism_for)
    Bepipred_protein_data <- fread(Bepipred_protein_data_file, header = T, sep = "\t", na.strings = NULL)
    
    normalized_protein_data_for <- merge(normalized_protein_data_for,
                                         Bepipred_protein_data[, .(id, protein_Bepipred_raw_score = BepiPred)],
                                         by = "id")
    
    ##Add the antigenic tag data
    protein_antigenic_tag_data_file <- sprintf("%s/%s/protein_antigenic_tag.tsv", input_data_folder, organism_for)
    protein_antigenic_tag_data <- fread(protein_antigenic_tag_data_file, header = T, sep = "\t", na.strings = NULL)
    
    normalized_protein_data_for <- merge(normalized_protein_data_for,
                                         protein_antigenic_tag_data,
                                         by = "id")
    
    #########################################-
    #### Create the leave one out models ####
    #########################################-
    leaveOneOut_full_balanced_training_set <- full_balanced_training_set[organism != organism_for]
    
    protein_leaveOneOut_model <- glm(in_antigenic_cluster ~ ., family = binomial(link='logit'), data = leaveOneOut_full_balanced_training_set[, -c("organism")])
    
    protein_leaveOneOut_noBepipred_model <- glm(in_antigenic_cluster ~ ., family = binomial(link='logit'), data = leaveOneOut_full_balanced_training_set[, -c("organism", "BepiPred")])
    
    #############################################-
    #### Predict the scores for the organism ####
    #############################################-
    normalized_protein_data_for$protein_generic_score <- predict(protein_balanced_generic_model, newdata = normalized_protein_data_for, type="response")
    normalized_protein_data_for[, protein_generic_score := round(protein_generic_score, 4)]
    normalized_protein_data_for[, protein_normalized_generic_score := round(normalizeLinear(protein_generic_score), 4)]
    
    normalized_protein_data_for$protein_leaveOneOut_score <- predict(protein_leaveOneOut_model, newdata = normalized_protein_data_for, type="response")
    normalized_protein_data_for[, protein_leaveOneOut_score := round(protein_leaveOneOut_score, 4)]
    normalized_protein_data_for[, protein_normalized_leaveOneOut_score := round(normalizeLinear(protein_leaveOneOut_score), 4)]    
    
    normalized_protein_data_for[, protein_Bepipred_raw_score := round(protein_Bepipred_raw_score, 4)]
    normalized_protein_data_for[, protein_normalized_Bepipred_raw_score := round(normalizeLinear(protein_Bepipred_raw_score), 4)]
    
    normalized_protein_data_for$protein_generic_noBepipred_score <- predict(protein_balanced_generic_noBepipred_model, newdata = normalized_protein_data_for, type="response")
    normalized_protein_data_for[, protein_generic_noBepipred_score := round(protein_generic_noBepipred_score, 4)]
    normalized_protein_data_for[, protein_normalized_generic_noBepipred_score := round(normalizeLinear(protein_generic_noBepipred_score), 4)]
    
    normalized_protein_data_for$protein_leaveOneOut_noBepipred_score <- predict(protein_leaveOneOut_noBepipred_model, newdata = normalized_protein_data_for, type="response")
    normalized_protein_data_for[, protein_leaveOneOut_noBepipred_score := round(protein_leaveOneOut_noBepipred_score, 4)]
    normalized_protein_data_for[, protein_normalized_leaveOneOut_noBepipred_score := round(normalizeLinear(protein_leaveOneOut_noBepipred_score), 4)]
    
    #########################-
    #### Save the scores ####
    #########################-
    output_organism_folder_aux <- sprintf("%s/%s", output_data_folder, organism_for)
    if (!(dir.exists(output_organism_folder_aux))) {
        dir.create(output_organism_folder_aux)
    }
    
    output_aux <- normalized_protein_data_for[, .(id,
                                                  protein_generic_score, protein_leaveOneOut_score, protein_Bepipred_raw_score,
                                                  protein_normalized_generic_score, protein_normalized_leaveOneOut_score, protein_normalized_Bepipred_raw_score,
                                                  protein_generic_noBepipred_score, protein_leaveOneOut_noBepipred_score,
                                                  protein_normalized_generic_noBepipred_score, protein_normalized_leaveOneOut_noBepipred_score)]
    protein_scores_file <- sprintf("%s/%s", output_organism_folder_aux, protein_scores_file_name)
    write.table(output_aux, file = protein_scores_file, col.names = T, row.names = F, sep = "\t", quote = T)
    
    ######################-
    #### **PEPTIDES** ####
    ######################-
    #############################################-
    #### Read specific organism peptide data ####
    #############################################-
    peptide_antigenic_tag_data_file <- sprintf("%s/%s/peptide_antigenic_tag.tsv", input_data_folder, organism_for)
    if (file.exists(peptide_antigenic_tag_data_file)) {
        #Peptide data
        normalized_peptide_data_file <- sprintf("%s/%s/normalized_output_per_peptide.tsv", input_data_folder, organism_for)
        normalized_peptide_data_for <- fread(normalized_peptide_data_file, header = T, sep = "\t", na.strings = NULL)
        
        normalized_peptide_data_for$organism <- organism_for
        
        #Add Bepipred data
        Bepipred_peptide_data_file <- sprintf("%s/%s/Bepipred_raw_output_per_peptide.tsv", input_data_folder, organism_for)
        Bepipred_peptide_data <- fread(Bepipred_peptide_data_file, header = T, sep = "\t", na.strings = NULL)
        
        normalized_peptide_data_for <- merge(normalized_peptide_data_for,
                                             Bepipred_peptide_data[, .(id, start, peptide_Bepipred_raw_score = BepiPred)],
                                             by = c("id", "start"))
        
        ##Add the antigenic tag data
        peptide_antigenic_tag_data_file <- sprintf("%s/%s/peptide_antigenic_tag.tsv", input_data_folder, organism_for)
        peptide_antigenic_tag_data <- fread(peptide_antigenic_tag_data_file, header = T, sep = "\t", na.strings = NULL)
        
        normalized_peptide_data_for <- merge(normalized_peptide_data_for,
                                             peptide_antigenic_tag_data,
                                             by = c("id", "start"))
        
        #########################################-
        #### Create the leave one out models ####
        #########################################-
        leaveOneOut_full_balanced_peptide_training_set <- full_balanced_peptide_training_set[organism != organism_for]
        
        peptide_leaveOneOut_model <- glm(is_antigenic ~ ., family = binomial(link='logit'), data = leaveOneOut_full_balanced_peptide_training_set[, -c("organism")])
        
        peptide_leaveOneOut_noBepipred_model <- glm(is_antigenic ~ ., family = binomial(link='logit'), data = leaveOneOut_full_balanced_peptide_training_set[, -c("organism", "BepiPred")])
        
        #############################################-
        #### Predict the scores for the organism ####
        #############################################-
        normalized_peptide_data_for$peptide_generic_score <- predict(peptide_balanced_generic_model, newdata = normalized_peptide_data_for, type="response")
        normalized_peptide_data_for[, peptide_generic_score := round(peptide_generic_score, 4)]    
        normalized_peptide_data_for[, peptide_normalized_generic_score := round(normalizeLinear(peptide_generic_score), 4)]
        
        normalized_peptide_data_for$peptide_leaveOneOut_score <- predict(peptide_leaveOneOut_model, newdata = normalized_peptide_data_for, type="response")
        normalized_peptide_data_for[, peptide_leaveOneOut_score := round(peptide_leaveOneOut_score, 4)]    
        normalized_peptide_data_for[, peptide_normalized_leaveOneOut_score := round(normalizeLinear(peptide_leaveOneOut_score), 4)]    
        
        normalized_peptide_data_for[, peptide_Bepipred_raw_score := round(peptide_Bepipred_raw_score, 4)]
        normalized_peptide_data_for[, peptide_normalized_Bepipred_raw_score := round(normalizeLinear(peptide_Bepipred_raw_score), 4)]
        
        normalized_peptide_data_for$peptide_generic_noBepipred_score <- predict(peptide_balanced_generic_noBepipred_model, newdata = normalized_peptide_data_for, type="response")
        normalized_peptide_data_for[, peptide_generic_noBepipred_score := round(peptide_generic_noBepipred_score, 4)]    
        normalized_peptide_data_for[, peptide_normalized_generic_noBepipred_score := round(normalizeLinear(peptide_generic_noBepipred_score), 4)]
        
        normalized_peptide_data_for$peptide_leaveOneOut_noBepipred_score <- predict(peptide_leaveOneOut_noBepipred_model, newdata = normalized_peptide_data_for, type="response")
        normalized_peptide_data_for[, peptide_leaveOneOut_noBepipred_score := round(peptide_leaveOneOut_noBepipred_score, 4)]    
        normalized_peptide_data_for[, peptide_normalized_leaveOneOut_noBepipred_score := round(normalizeLinear(peptide_leaveOneOut_noBepipred_score), 4)]
        
        ############################################-
        #### Add Protein Scores to Peptide Data ####
        ############################################-
        normalized_peptide_data_for <- merge(normalized_peptide_data_for,
                                             normalized_protein_data_for[, .(id,
                                                                             protein_normalized_generic_score = protein_normalized_generic_score,
                                                                             protein_normalized_leaveOneOut_score = protein_normalized_leaveOneOut_score,
                                                                             protein_normalized_Bepipred_raw_score = protein_normalized_Bepipred_raw_score,
                                                                             protein_normalized_generic_noBepipred_score = protein_normalized_generic_noBepipred_score,
                                                                             protein_normalized_leaveOneOut_noBepipred_score = protein_normalized_leaveOneOut_noBepipred_score)],
                                             by = "id")
        
        normalized_peptide_data_for[, peptide_normalized_combined_generic_score := protein_normalized_generic_score + peptide_normalized_generic_score]
        normalized_peptide_data_for[, peptide_normalized_combined_generic_score := round(normalizeLinear(peptide_normalized_combined_generic_score), 4)]
        
        normalized_peptide_data_for[, peptide_normalized_combined_leaveOneOut_score := protein_normalized_leaveOneOut_score + peptide_leaveOneOut_score]
        normalized_peptide_data_for[, peptide_normalized_combined_leaveOneOut_score := round(normalizeLinear(peptide_normalized_combined_leaveOneOut_score), 4)]
        
        normalized_peptide_data_for[, peptide_normalized_combined_Bepipred_raw_score := protein_normalized_Bepipred_raw_score + peptide_normalized_Bepipred_raw_score]
        normalized_peptide_data_for[, peptide_normalized_combined_Bepipred_raw_score := round(normalizeLinear(peptide_normalized_combined_Bepipred_raw_score), 4)]
        
        normalized_peptide_data_for[, peptide_normalized_combined_generic_noBepipred_score := protein_normalized_generic_noBepipred_score + peptide_normalized_generic_noBepipred_score]
        normalized_peptide_data_for[, peptide_normalized_combined_generic_noBepipred_score := round(normalizeLinear(peptide_normalized_combined_generic_noBepipred_score), 4)]
        
        normalized_peptide_data_for[, peptide_normalized_combined_leaveOneOut_noBepipred_score := protein_normalized_leaveOneOut_noBepipred_score + peptide_normalized_leaveOneOut_noBepipred_score]
        normalized_peptide_data_for[, peptide_normalized_combined_leaveOneOut_noBepipred_score := round(normalizeLinear(peptide_normalized_combined_leaveOneOut_noBepipred_score), 4)]

        #########################-
        #### Save the scores ####
        #########################-
        output_aux <- normalized_peptide_data_for[, .(id, start,
                                                      peptide_generic_score, peptide_leaveOneOut_score, peptide_Bepipred_raw_score,
                                                      peptide_normalized_generic_score, peptide_normalized_leaveOneOut_score, peptide_normalized_Bepipred_raw_score,
                                                      peptide_normalized_combined_generic_score, peptide_normalized_combined_leaveOneOut_score, peptide_normalized_combined_Bepipred_raw_score,
                                                      peptide_generic_noBepipred_score, peptide_leaveOneOut_noBepipred_score,
                                                      peptide_normalized_generic_noBepipred_score, peptide_normalized_leaveOneOut_noBepipred_score,
                                                      peptide_normalized_combined_generic_noBepipred_score, peptide_normalized_combined_leaveOneOut_noBepipred_score)]
        peptide_scores_file <- sprintf("%s/%s/%s", output_data_folder, organism_for, peptide_scores_file_name)
        write.table(output_aux, file = peptide_scores_file, col.names = T, row.names = F, sep = "\t", quote = T)
    }
}
