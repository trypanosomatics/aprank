####################-
#### **CONFIG** ####
####################-
aprank_development_folder <- "/home/user/Desktop/APRANK/ModelDevelopment"
aprank_folder <- "/home/user/Desktop/APRANK"

library(data.table)

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
output_data_folder <- sprintf("%s/21_outputs", aprank_development_folder)
scripts_folder <- sprintf("%s/source", aprank_folder)

peptide_length <- 15
peptide_overlap <- 14

max_protein_length <- 9999
replace_nonAA_chars_by <- "X"

output_fasta_file <- sprintf("%s/blast/all_organisms.fasta", input_data_folder)
output_tabbed_file <- sprintf("%s/blast/all_organisms_tabbed.tsv", input_data_folder)

##################-
#### **MAIN** ####
##################-
for (organism_for in organisms) {
    # organism_for <- organisms[1]
    ###############################-
    #### Create Needed Folders ####
    ###############################-
    organism_temp_data_folder <- sprintf("%s/%s-%s-%s", temp_data_folder, organism_for,
                                         peptide_length, peptide_overlap)
    if (!dir.exists(organism_temp_data_folder)) {
        dir.create(organism_temp_data_folder)
    }
    
    ##################################-
    #### Specific Organism Config ####
    ##################################-
    input_fasta_file <- sprintf("%s/%s.fasta", input_data_folder, organism_for)
    
    temp_fasta_file <- sprintf("%s/temp.fasta", organism_temp_data_folder)
    temp_tabbed_fasta_file <- sprintf("%s/temp_fasta_tabbed.tsv", organism_temp_data_folder)
    temp_splitted_fasta_file <- sprintf("%s/temp_fasta_splitted.tsv", organism_temp_data_folder)
    temp_aminoacid_fasta_file <- sprintf("%s/temp_fasta_aminoacids.tsv", organism_temp_data_folder)
    
    ########################################-
    #### Get the FASTA for the Organism ####
    ########################################-
    #I honestly just want the tabbed fasta, but reusing code is easier
    command_aux <- sprintf("perl '%s/lib/tab_and_split_proteome.pl' '%s' %s %s %s '%s' '%s' '%s' '%s' '%s'",
                           scripts_folder,
                           input_fasta_file,
                           peptide_length,
                           peptide_overlap,
                           max_protein_length,
                           replace_nonAA_chars_by,
                           temp_fasta_file,
                           temp_tabbed_fasta_file,
                           temp_splitted_fasta_file,
                           temp_aminoacid_fasta_file)
    system(command_aux)
    
    tabbed_fasta <- fread(temp_tabbed_fasta_file, header = T, sep = "\t", na.strings = NULL)
    tabbed_fasta$organism <- organism_for
    
    if (organism_for == organisms[1]) {
        full_tabbed_fasta <- tabbed_fasta   
    } else {
        full_tabbed_fasta <- rbindlist(list(full_tabbed_fasta, tabbed_fasta))
    }
}

#Assign an unique ID for all proteins across all proteomes
setnames(full_tabbed_fasta, "original_protein", "original_id")
setnames(full_tabbed_fasta, "protein", "inProteome_id")
full_tabbed_fasta[, id := sprintf("p_%s", .I)]

#Save the mixed FASTA
full_tabbed_fasta[, output_aux:= sprintf(">%s\n%s", id, sequence)]
output <- paste(full_tabbed_fasta$output_aux, collapse = "\n")
full_tabbed_fasta <- full_tabbed_fasta[, -c("output_aux")]
write(output, file = output_fasta_file)

#Save the data of which protein is which
full_tabbed_fasta <- full_tabbed_fasta[, .(organism, id, inProteome_id, original_id, sequence)]
write.table(full_tabbed_fasta, file = output_tabbed_file, col.names = T, row.names = F, sep = "\t", quote = T)
