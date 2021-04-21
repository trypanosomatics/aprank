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
library(foreach) #parallel processing
library(parallel) #parallel processing
library(doParallel) #parallel processing

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

setwd(aprank_folder)

input_data_folder <- sprintf("%s/01_inputs", aprank_development_folder)
temp_data_folder <- sprintf("%s/88_tmp", aprank_development_folder)
output_data_folder <- sprintf("%s/01_inputs/blast/outputs", aprank_development_folder)
scripts_folder <- sprintf("%s/source", aprank_folder)

all_organisms_fasta_file <- sprintf("%s/blast/all_organisms.fasta", input_data_folder)
all_organisms_tabbed_fasta_file <- sprintf("%s/blast/all_organisms_tabbed.tsv", input_data_folder)

number_of_parallel_processes <- 4

##################-
#### **MAIN** ####
##################-
all_organisms_tabbed_fasta <- fread(all_organisms_tabbed_fasta_file, header = T, sep = "\t", na.strings = NULL)

for (organism_for in organisms) {
    # organism_for <- organisms[1]
    ##################################-
    #### Specific Organism Config ####
    ##################################-
    input_fasta_file <- sprintf("%s/%s.fasta", input_data_folder, organism_for)
    output_BLAST_file <- sprintf("%s/%s_BLAST.tsv", output_data_folder, organism_for)
    
    tabbed_fasta <- all_organisms_tabbed_fasta[organism == organism_for, .(protein = id, sequence)]
    
    ########################################-
    #### Set up the Parallel Processing ####
    ########################################-
    #Create the subFASTAs for parallel processing
    sub_temp_fastas_files <- c()
    tabbed_fasta$parallel_process_i <- 1
    if (number_of_parallel_processes > 1) {
        tabbed_fasta[, output_aux := sprintf(">%s\n%s", protein, sequence)]
        
        unique_protein_amount <- tabbed_fasta[, .N]
        
        division_aux <- unique_protein_amount %/% number_of_parallel_processes
        reminder_aux <- unique_protein_amount %% number_of_parallel_processes
        protein_per_process_amount <-c(0, rep(division_aux, number_of_parallel_processes))
        protein_per_process_amount <- protein_per_process_amount + c(0, rep(1, reminder_aux), rep(0, number_of_parallel_processes - reminder_aux))
        for (i in 1:number_of_parallel_processes) {
            # i <- 1
            protein_subset_min_for_parallel <- 1 + sum(protein_per_process_amount[0:i])
            protein_subset_max_for_parallel <- protein_subset_min_for_parallel + protein_per_process_amount[i + 1] - 1
            tabbed_fasta[protein_subset_min_for_parallel:protein_subset_max_for_parallel]$parallel_process_i <- i
            
            sub_fasta_output <- paste(tabbed_fasta[protein_subset_min_for_parallel:protein_subset_max_for_parallel]$output_aux, collapse = "\n")
            sub_temp_fasta_file <- sprintf("%s/sub_temp_%s.fasta", temp_data_folder, i)
            
            write(sub_fasta_output, file = sub_temp_fasta_file)
            
            sub_temp_fastas_files <- c(sub_temp_fastas_files, sub_temp_fasta_file)
        }
        
        tabbed_fasta <- tabbed_fasta[, -c("output_aux")]
    } else {
        sub_temp_fastas_files <- c(temp_fasta_file)
    }
    
    ########################-
    #### Run the BLASTs ####
    ########################-
    numCores <- detectCores() #set the cores
    
    registerDoParallel(numCores) #use multicore, set to the number of our cores
    
    full_BLAST_data <- foreach (i = 1:number_of_parallel_processes, .combine = rbind) %dopar% {
        BLAST_output_file <- sprintf("%s/blast_output_temp_%s", temp_data_folder, i)
        # command_aux <- sprintf("blastp -db '%s' -query '%s' -outfmt '6 qseqid sseqid pident length mismatch gapopen evalue bitscore qseq sseq sstart send' -word_size 2 -max_target_seqs 10000 -matrix BLOSUM80",
        #                        all_organisms_fasta_file,
        #                        sub_temp_fastas_files[i])
        # BLAST_output <- system(command_aux, intern = T)
        command_aux <- sprintf("blastp -db '%s' -query '%s' -outfmt '6 qseqid sseqid pident length mismatch gapopen evalue bitscore qseq sseq sstart send' -word_size 2 -max_target_seqs 10000 -matrix BLOSUM80 > '%s'",
                               all_organisms_fasta_file,
                               sub_temp_fastas_files[i],
                               BLAST_output_file)
        system(command_aux)
        
        # BLAST_parsed_data <- read.delim(text = BLAST_output, header = F, sep = "", na.strings = NULL, comment.char = "")
        BLAST_parsed_data <- read.delim(file = BLAST_output_file, header = F, sep = "", na.strings = NULL, comment.char = "")
        BLAST_col_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "evalue", "bitscore", "qseq", "sseq", "sstart", "send")
        colnames(BLAST_parsed_data) <- BLAST_col_names
        BLAST_parsed_data <- as.data.table(BLAST_parsed_data)
        
        BLAST_parsed_data
    }
    #End the parallel processing
    stopImplicitCluster()
    
    write.table(full_BLAST_data, file = output_BLAST_file, col.names = T, row.names = F, sep = "\t", quote = T)
}
