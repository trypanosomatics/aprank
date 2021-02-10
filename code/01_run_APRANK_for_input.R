####################-
#### **CONFIG** ####
####################-
aprank_development_folder <- "/home/user/Desktop/APRANK/ModelDevelopment"
aprank_folder <- "/home/user/Desktop/APRANK"
    
library(data.table)

############################################-
#### Config - Select Organisms to parse ####
############################################-
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

####################################-
#### Config - Select Predictors ####
####################################-
## General
use_BepiPred <- 1
use_IsoelectricPoint <- 1
use_Iupred <- 1
use_MolecularWeight <- 1
use_NetMHCIIpan <- 1
use_NetOglyc <- 1
use_NetSurfp <- 1
use_Paircoil2 <- 1
use_PredGPI <- 1
use_SignalP <- 1
use_TMHMM <- 1
use_Xstream <- 1

use_SelfSimilarity <- 1
use_CrossReactivity <- 1
use_Coendemicity <- 0 #no coendemicity here

##########################-
#### Config - General ####
##########################-
setwd(aprank_folder)

input_data_folder <- sprintf("%s/01_inputs", aprank_development_folder)
temp_data_folder <- sprintf("%s/88_tmp", aprank_development_folder)
output_data_folder <- sprintf("%s/21_outputs", aprank_development_folder)
scripts_folder <- sprintf("%s/source", aprank_folder)

peptide_length <- 15
peptide_overlap <- 14

max_protein_length <- 9999
replace_nonAA_chars_by <- "X"

#################################-
#### Config - Run Predictors ####
#################################-
## General
source(sprintf("%s/01_runPredictors.R", scripts_folder))

#This has to be 1 or more (1 meaning just 1 process without parallelization)
number_of_parallel_processes <- 4
#NetSurfp uses multiple cores per process. I would give 4 cores to each process
#This has to be 1 or more (1 meaning just 1 process without parallelization)
number_of_parallel_processes_for_NetSurfp <- max(1, floor(number_of_parallel_processes / 4))

## SignalP
#The organism group of the genome to analyze. It can be: "euk", "gram+" or "gram-".
SignalP_organism_group_list <- list()
SignalP_organism_group_list[["Borrelia_burgdorferi"]] <- "gram-"
SignalP_organism_group_list[["Brucella_melitensis"]] <- "gram-"
SignalP_organism_group_list[["Coxiella_burnetii"]] <- "gram-"
SignalP_organism_group_list[["Escherichia_coli_H10407"]] <- "gram-"
SignalP_organism_group_list[["Francisella_tularensis"]] <- "gram-"
SignalP_organism_group_list[["Leptospira_interrogans"]] <- "gram-"
SignalP_organism_group_list[["Porphyromonas_gingivalis_837"]] <- "gram-"
SignalP_organism_group_list[["Mycobacterium_leprae"]] <- "gram+"
SignalP_organism_group_list[["Mycobacterium_tuberculosis"]] <- "gram+"
SignalP_organism_group_list[["Staphylococcus_aureus"]] <- "gram+"
SignalP_organism_group_list[["Streptococcus_pyogenes_serotype_M1"]] <- "gram+"
SignalP_organism_group_list[["Leishmania_brasilensis"]] <- "euk"
SignalP_organism_group_list[["Plasmodium_falciparum"]] <- "euk"
SignalP_organism_group_list[["Trypanosoma_cruzi"]] <- "euk"
SignalP_organism_group_list[["Toxoplasma_gondi"]] <- "euk"

## Xstream
Xstream_path <- "/usr/local/xstream/xstream.jar"

## NetMHCIIPan
# The length of the fragment used when analyzing with NetMHCIIpan (integer between 9 and 50).
NetMHCIIpan_binding_peptide_length <- 9
# The names of the MHC class II alleles to consider when analyzing
NetMHCIIPan_alleles <- c("DRB1_0101",
                         "DRB3_0101",
                         "DRB4_0101",
                         "DRB5_0101")

## Self Similarity, CrossReactivity & Coendemicity
KmerSimilarity_kmer_length <- 6

## CrossReactivity (Pathogen vs Host)
CrossReactivity_fasta_file <- sprintf("%s/host_fasta/human.PRJNA178030.prot.fasta", input_data_folder)

#########################################-
#### Config - Normalize protein data ####
#########################################-
## General
source(sprintf("%s/11_normalizeProteinData.R", scripts_folder))

source(sprintf("%s/lib/normalization.R", scripts_folder))
source(sprintf("%s/lib/aux.R", scripts_folder))

## Paircoil2
#Here I'm asking if there is a fragment of length Paircoil2_fragment_length of aa above Paircoil2_threshold
Paircoil2_fragment_length <- 50
Paircoil2_threshold <- 0.5

## Xstream
Xstream_min_period <- 1
Xstream_min_copy_number <- 1
Xstream_max_consensus_error <- 1

## Coendemicity
# Coendemicity_protein_min_coendemic_amount_for_penalty is the min amount of times a given kmer has to appear in the
# coendemic proteome to apply the penalty to proteins containing it
Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty <- 0
# Coendemicity_protein_start_penalty_proportion is the starting proportion to apply a penalty to the score
# This penalty goes up until the proportion reaches the Coendemicity_protein_max_penalty_proportion, where the penalty
# becomes 1 and the score 0
Coendemicity_protein_start_penalty_proportion <- 0
Coendemicity_protein_max_penalty_proportion <- 1/3

#########################################-
#### Config - Normalize peptide data ####
#########################################-
## General
source(sprintf("%s/12_normalizePeptideData.R", scripts_folder))

source(sprintf("%s/lib/normalization.R", scripts_folder))
source(sprintf("%s/lib/aux.R", scripts_folder))

## Paircoil2
#Here I'm asking if there is a fragment of length Paircoil2_fragment_length of aa above Paircoil2_threshold
# Paircoil2_threshold <- 0.5

## Coendemicity
# Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty is the min amount of times a given kmer has to appear in the
# coendemic proteome to apply the penalty to peptides containing it
Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty <- 0
# Coendemicity_peptide_start_penalty_proportion is the starting proportion to apply a penalty to the score
# This penalty goes up until the proportion reaches the Coendemicity_protein_max_penalty_proportion, where the penalty
# becomes 1 and the score 0
Coendemicity_peptide_start_penalty_proportion <- 0
Coendemicity_peptide_max_penalty_proportion <- 2/15

##################-
#### **MAIN** ####
##################-
for (organism_for in organisms) {
    #organism_for <- organisms[1]
    
    ###############################-
    #### Create Needed Folders ####
    ###############################-
    organism_temp_data_folder <- sprintf("%s/%s", temp_data_folder, organism_for)
    organism_output_data_folder <- sprintf("%s/%s", output_data_folder, organism_for)
    if (!dir.exists(organism_temp_data_folder)) {
        dir.create(organism_temp_data_folder)
    }
    if (!dir.exists(organism_output_data_folder)) {
        dir.create(organism_output_data_folder)
    }
    
    ##################################-
    #### Specific Organism Config ####
    ##################################-
    input_fasta_file <- sprintf("%s/%s.fasta", input_data_folder, organism_for)
    
    temp_fasta_file <- sprintf("%s/temp.fasta", organism_temp_data_folder)
    temp_tabbed_fasta_file <- sprintf("%s/temp_fasta_tabbed.tsv", organism_temp_data_folder)
    temp_splitted_fasta_file <- sprintf("%s/temp_fasta_splitted.tsv", organism_temp_data_folder)
    temp_aminoacid_fasta_file <- sprintf("%s/temp_fasta_aminoacids.tsv", organism_temp_data_folder)
    
    #The output method can be "write" (makes files), "list" (returns a list), or "both"
    predictors_output_method <- "both"
    output_per_aminoacid_file <- sprintf("%s/output_per_aminoacid.tsv", organism_output_data_folder)
    output_per_protein_file <- sprintf("%s/output_per_protein.tsv", organism_output_data_folder)
    output_per_kmer_file <- sprintf("%s/output_per_kmer.tsv", organism_output_data_folder)
    output_per_peptide_file <- sprintf("%s/output_per_peptide.tsv", organism_output_data_folder)
    output_per_repeat_file <- sprintf("%s/output_per_repeat.tsv", organism_output_data_folder)
    
    #SignalP organism
    SignalP_organism_group <- SignalP_organism_group_list[[organism_for]]
    
    #The output method can be "write" (makes files), "list" (returns a list), or "both"
    protein_normalization_output_method <- "write"
    normalized_output_per_protein_file <- sprintf("%s/normalized_output_per_protein.tsv", organism_output_data_folder)
    
    #The output method can be "write" (makes files), "list" (returns a list), or "both"
    peptide_normalization_output_method <- "write"
    normalized_output_per_peptide_file <- sprintf("%s/normalized_output_per_peptide.tsv", organism_output_data_folder)
    
    ########################-
    #### Run Predictors ####
    ########################-
    list_output <- 
        runPredictors(input_fasta_file = input_fasta_file,
                      input_data_folder = input_data_folder, temp_data_folder = organism_temp_data_folder, output_data_folder = organism_output_data_folder, scripts_folder = scripts_folder,
                      temp_fasta_file = temp_fasta_file, temp_tabbed_fasta_file = temp_tabbed_fasta_file, temp_splitted_fasta_file = temp_splitted_fasta_file, temp_aminoacid_fasta_file = temp_aminoacid_fasta_file,
                      output_method = predictors_output_method,
                      output_per_aminoacid_file = output_per_aminoacid_file, output_per_protein_file = output_per_protein_file, output_per_kmer_file = output_per_kmer_file, output_per_peptide_file = output_per_peptide_file, output_per_repeat_file = output_per_repeat_file,
                      number_of_parallel_processes = number_of_parallel_processes, number_of_parallel_processes_for_NetSurfp = number_of_parallel_processes_for_NetSurfp,
                      peptide_length = peptide_length, peptide_overlap = peptide_overlap, max_protein_length = max_protein_length, replace_nonAA_chars_by = replace_nonAA_chars_by,
                      SignalP_organism_group = SignalP_organism_group,
                      Xstream_path = Xstream_path,
                      NetMHCIIpan_binding_peptide_length = NetMHCIIpan_binding_peptide_length, NetMHCIIPan_alleles = NetMHCIIPan_alleles,
                      KmerSimilarity_kmer_length = KmerSimilarity_kmer_length, CrossReactivity_fasta_file = CrossReactivity_fasta_file, Coendemicity_fasta_file = Coendemicity_fasta_file,
                      use_BepiPred = use_BepiPred, use_Paircoil2 = use_Paircoil2, use_PredGPI = use_PredGPI, use_SignalP = use_SignalP,
                      use_TMHMM = use_TMHMM, use_NetSurfp = use_NetSurfp,
                      use_Iupred = use_Iupred, use_NetOglyc = use_NetOglyc, use_Xstream = use_Xstream, use_NetMHCIIpan = use_NetMHCIIpan,
                      use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight = use_MolecularWeight,
                      use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
    
    output_per_aminoacid <- list_output[["output_per_aminoacid"]]
    output_per_protein <- list_output[["output_per_protein"]]
    output_per_kmer <- list_output[["output_per_kmer"]]
    output_per_peptide <- list_output[["output_per_peptide"]]
    output_per_repeat <- list_output[["output_per_repeat"]]
    rm(list_output)
    gc()
    
    ################################-
    #### Normalize Protein Data ####
    ################################-
    normalizeProteinData(output_per_aminoacid = output_per_aminoacid, output_per_protein = output_per_protein, output_per_kmer = output_per_kmer, output_per_peptide = output_per_peptide, output_per_repeat = output_per_repeat,
                         output_method = protein_normalization_output_method,
                         normalized_output_per_protein_file = normalized_output_per_protein_file,
                         Paircoil2_fragment_length = Paircoil2_fragment_length, Paircoil2_threshold = Paircoil2_threshold,
                         Xstream_min_period = Xstream_min_period, Xstream_min_copy_number = Xstream_min_copy_number, Xstream_max_consensus_error = Xstream_max_consensus_error,
                         Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty = Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_protein_start_penalty_proportion = Coendemicity_protein_start_penalty_proportion, Coendemicity_protein_max_penalty_proportion = Coendemicity_protein_max_penalty_proportion,
                         use_BepiPred = use_BepiPred, use_Paircoil2 = use_Paircoil2, use_PredGPI = use_PredGPI, use_SignalP = use_SignalP,
                         use_TMHMM = use_TMHMM, use_NetSurfp = use_NetSurfp,
                         use_Iupred = use_Iupred, use_NetOglyc = use_NetOglyc, use_Xstream = use_Xstream, use_NetMHCIIpan = use_NetMHCIIpan,
                         use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight = use_MolecularWeight,
                         use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
    
    ################################-
    #### Normalize Peptide Data ####
    ################################-
    normalizePeptideData(output_per_aminoacid = output_per_aminoacid, output_per_protein = output_per_protein, output_per_kmer = output_per_kmer, output_per_peptide = output_per_peptide, output_per_repeat = output_per_repeat,
                         output_method = peptide_normalization_output_method,
                         normalized_output_per_peptide_file = normalized_output_per_peptide_file,
                         peptide_length = peptide_length, peptide_overlap = peptide_overlap,
                         Paircoil2_threshold = Paircoil2_threshold,
                         Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty = Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_peptide_start_penalty_proportion = Coendemicity_peptide_start_penalty_proportion, Coendemicity_peptide_max_penalty_proportion = Coendemicity_peptide_max_penalty_proportion,
                         use_BepiPred = use_BepiPred, use_Paircoil2 = use_Paircoil2, use_PredGPI = use_PredGPI, use_SignalP = use_SignalP,
                         use_TMHMM = use_TMHMM, use_NetSurfp = use_NetSurfp,
                         use_Iupred = use_Iupred, use_NetOglyc = use_NetOglyc, use_Xstream = use_Xstream, use_NetMHCIIpan = use_NetMHCIIpan,
                         use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight = use_MolecularWeight,
                         use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity)
}
