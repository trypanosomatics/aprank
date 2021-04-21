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

####################-
#### **CONFIG** ####
####################-
aprank_development_folder <- "/home/user/Desktop/APRANK/ModelDevelopment"
aprank_folder <- "/home/user/Desktop/APRANK"

library(data.table)
library(ROSE)

###########################################################-
#### Config - Select Organisms to include in the model ####
###########################################################-
#Here I don't use all organisms because for some of them we don't have the antigenic information at peptide level
organisms <- c("Borrelia_burgdorferi",
               "Leptospira_interrogans",
               "Mycobacterium_leprae",
               "Mycobacterium_tuberculosis",
               "Plasmodium_falciparum",
               "Toxoplasma_gondi",
               "Trypanosoma_cruzi",
               "Leishmania_brasilensis",
               "Staphylococcus_aureus",
               "Streptococcus_pyogenes_serotype_M1",
               "Escherichia_coli_H10407",
               "Porphyromonas_gingivalis_837")

##########################-
#### Config - General ####
##########################-
output_data_folder <- sprintf("%s/21_outputs", aprank_development_folder)
scripts_folder <- sprintf("%s/source", aprank_folder)

peptide_antigenic_tag_file_name <- "peptide_antigenic_tag.tsv"

#################################-
#### Config - Antigenic Tags ####
#################################-
antigenic_tags_folder <- sprintf("%s/11_antigens", aprank_development_folder)
antigenic_kmer_length <- 8

##################-
#### **MAIN** ####
##################-
###################-
#### Read data ####
###################-
for (organism_for in organisms) {
    # organism_for <- organisms[1]
    #Protein data (I need it for the protein name)
    normalized_protein_data_file <- sprintf("%s/%s/normalized_output_per_protein.tsv", output_data_folder, organism_for)
    normalized_protein_data_for <- fread(normalized_protein_data_file, header = T, sep = "\t", na.strings = NULL)
    
    normalized_protein_data_for <- normalized_protein_data_for[, .(id, original_id)]
    
    #Peptide data
    normalized_peptide_data_file <- sprintf("%s/%s/normalized_output_per_peptide.tsv", output_data_folder, organism_for)
    normalized_peptide_data_for <- fread(normalized_peptide_data_file, header = T, sep = "\t", na.strings = NULL)
    
    normalized_peptide_data_for <- merge(normalized_peptide_data_for,
                                         normalized_protein_data_for,
                                         by = "id",
                                         all.x = T)
    
    normalized_peptide_data_for$organism <- organism_for
    
    if (organism_for == organisms[1]) {
        normalized_peptide_data <- normalized_peptide_data_for
    } else {
        normalized_peptide_data <- rbindlist(list(normalized_peptide_data, normalized_peptide_data_for))
    }
    
    #Antigen data
    antigenic_tag_peptide_data_file <- sprintf("%s/%s_antigens_peptides", antigenic_tags_folder, organism_for)
    antigenic_tag_peptide_data_for <- fread(antigenic_tag_peptide_data_file, header = T, sep = "\t", na.strings = NULL)
    
    antigenic_tag_peptide_data_for$organism <- organism_for
    
    if (organism_for == organisms[1]) {
        antigenic_tag_peptide_data <- antigenic_tag_peptide_data_for
    } else {
        antigenic_tag_peptide_data <- rbindlist(list(antigenic_tag_peptide_data, antigenic_tag_peptide_data_for))
    }
}
rm(normalized_protein_data_for)
rm(normalized_peptide_data_for)
rm(antigenic_tag_peptide_data_for)
gc()

###########################################-
#### Assign positive and negative tags ####
###########################################-
#Here I'm going to do a "kmer expansion" of the antigens, where each antigen is divided in kmers of a given length
#and then I assign as antigenic all peptides inside that antigen's protein that includes said kmers

#Keep only antigens at least as long as the kmers
antigenic_tag_peptide_data[, sequence_length := nchar(sequence)]
antigenic_tag_peptide_data <- antigenic_tag_peptide_data[sequence_length >= antigenic_kmer_length]

unique_proteins_with_antigenic_peptides <- unique(antigenic_tag_peptide_data$protein)
for (protein_for in unique_proteins_with_antigenic_peptides) {
    # protein_for <- unique_proteins_with_antigenic_peptides[1]
    antigenic_peptides <- antigenic_tag_peptide_data[protein == protein_for]
    posible_antigenic_peptides <- normalized_peptide_data[original_id == protein_for, .(original_id, start, peptide)]
    
    #Calculate the kmers for the known antigens for the protein
    antigenic_peptides[, kmer_amount := sequence_length - antigenic_kmer_length + 1]
    antigenic_kmers <- antigenic_peptides[rep(seq(.N), kmer_amount), -c("kmer_amount")]
    antigenic_kmers[, kmer_start := .SD[, .I], by = .(protein, start)]
    antigenic_kmers[, kmer := substring(sequence, kmer_start, kmer_start + antigenic_kmer_length - 1)]
    antigenic_kmers <- unique(antigenic_kmers$kmer)
    
    #See what peptides in the protein contains the antigenic kmers
    posible_antigenic_peptides[, peptide_length := nchar(peptide)]
    posible_antigenic_peptides[, kmer_amount := peptide_length - antigenic_kmer_length + 1]
    posible_antigenic_peptides <- posible_antigenic_peptides[rep(seq(.N), kmer_amount), -c("kmer_amount")]
    posible_antigenic_peptides[, kmer_start := .SD[, .I], by = .(original_id, start)]
    posible_antigenic_peptides[, kmer := substring(peptide, kmer_start, kmer_start + antigenic_kmer_length - 1)]
    new_antigenic_peptides <- posible_antigenic_peptides[kmer %in% antigenic_kmers]
    new_antigenic_peptides <- unique(new_antigenic_peptides[, .(original_id, start)])
    
    if (protein_for == unique_proteins_with_antigenic_peptides[1]) {
        full_new_antigenic_peptides <- new_antigenic_peptides
    } else {
        full_new_antigenic_peptides <- rbindlist(list(full_new_antigenic_peptides, new_antigenic_peptides))
    }
}
full_new_antigenic_peptides[, is_antigenic := 1]

normalized_peptide_data <- merge(normalized_peptide_data,
                                 full_new_antigenic_peptides,
                                 by = c("original_id", "start"),
                                 all.x = T)
normalized_peptide_data[is.na(is_antigenic), is_antigenic := 0]

rm(antigenic_tag_peptide_data)
rm(full_new_antigenic_peptides)
gc()

#######################-
#### Save the tags ####
#######################-
for (organism_for in organisms) {
    # organism_for <- organisms[1]
    sub_normalized_peptide_data <- normalized_peptide_data[organism == organism_for]
    
    #Calculate the output file
    organism_output_data_folder <- sprintf("%s/%s", output_data_folder, organism_for)
    peptide_antigenic_tag_file <- sprintf("%s/%s", organism_output_data_folder, peptide_antigenic_tag_file_name)
    
    #Keep the important columns
    sub_normalized_peptide_data <- sub_normalized_peptide_data[, .(id, start, is_antigenic)]
    
    #Save
    write.table(sub_normalized_peptide_data, file = peptide_antigenic_tag_file, col.names = T, row.names = F, sep = "\t", quote = T)
}
