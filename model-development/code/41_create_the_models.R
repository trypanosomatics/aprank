####################-
#### **CONFIG** ####
####################-
aprank_development_folder <- "/home/user/Desktop/APRANK/ModelDevelopment"
aprank_folder <- "/home/user/Desktop/APRANK"

setwd(aprank_folder)

library(data.table)
library(ROSE)

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
