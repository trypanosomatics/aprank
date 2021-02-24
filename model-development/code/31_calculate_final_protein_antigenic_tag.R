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

##########################-
#### Config - General ####
##########################-
setwd(aprank_folder)

output_data_folder <- sprintf("%s/21_outputs", aprank_development_folder)
scripts_folder <- sprintf("%s/source", aprank_folder)

protein_antigenic_tag_file_name <- "protein_antigenic_tag.tsv"

######################################-
#### Config - Similarity Clusters ####
######################################-
BLAST_output_data_folder <- sprintf("%s/01_inputs/blast/outputs", aprank_development_folder)
BLAST_protein_information_file <- sprintf("%s/01_inputs/blast/all_organisms_tabbed.tsv", aprank_development_folder)

min_pident_to_keep_match <- 0.75 #between 0 and 1
max_evalue_to_keep_match <- 1.0e-12 #length close to 30
#The min lenght fraction means that the match has to be a % of at least one of both sequences to be kept
min_match_length_proportion_to_keep_match <- 0.50 #between 0 and 1

#################################-
#### Config - Antigenic Tags ####
#################################-
antigenic_tags_folder <- sprintf("%s/11_antigens", aprank_development_folder)

###########################-
#### **AUX FUNCTIONS** ####
###########################-
BLAST_to_similarityClusters <- function(BLAST_data,
                                        normalized_protein_data,
                                        clustering_method = "single",
                                        min_pident_to_keep_match) {
    ### BLAST_data has to have the columns qseqid, sseqid and pident
    
    #Since all the matches should be mirrored, check this and calculate a single pident for both mirrors
    #completing some cases where there is only one of both cases
    mirror_BLAST_data <- BLAST_data[qseqid > sseqid]
    BLAST_data <- BLAST_data[qseqid < sseqid]
    
    mirror_BLAST_data[, aux := qseqid]
    mirror_BLAST_data[, qseqid := sseqid]
    mirror_BLAST_data[, sseqid := aux]
    setnames(mirror_BLAST_data, "pident", "mirror_pident")
    mirror_BLAST_data <- mirror_BLAST_data[, c("qseqid", "sseqid", "mirror_pident")]
    
    #Keep only the best match between proteins (for long proteins it's possible more than 1 good match)
    BLAST_data <- BLAST_data[, .(pident = max(pident)), by = .(qseqid, sseqid)]
    mirror_BLAST_data <- mirror_BLAST_data[, .(mirror_pident = max(mirror_pident)), by = .(qseqid, sseqid)]
    
    #Add the mirror data to the original and calculate the avarage (if both exists)
    BLAST_data <- merge(BLAST_data,
                        mirror_BLAST_data,
                        by = c("qseqid", "sseqid"),
                        all = TRUE)
    
    BLAST_data[is.na(pident), pident := mirror_pident]
    BLAST_data[is.na(mirror_pident), mirror_pident := pident]
    BLAST_data[, pident := round((pident + mirror_pident) / 2, 4)]
    
    #Add the matches against itself just for the proteins that have matches against other things
    proteins_with_hits <- unique(c(BLAST_data$qseqid, BLAST_data$sseqid))
    proteins_without_hits <- setdiff(normalized_protein_data$BLAST_id, proteins_with_hits)
    
    BLAST_data <- rbindlist(list(BLAST_data[, c("qseqid", "sseqid", "pident")],
                                 data.table(qseqid = proteins_with_hits,
                                            sseqid = proteins_with_hits,
                                            pident = 1)))
    
    #Create the matrix and make clusters of similarity
    BLAST_data$sseqid <- as.factor(BLAST_data$sseqid)
    BLAST_data$qseqid <- as.factor(BLAST_data$qseqid)
    
    similarity_matrix <- with(BLAST_data, {
        output <- matrix(nrow = nlevels(sseqid),
                         ncol = nlevels(qseqid),
                         dimnames = list(levels(sseqid),
                                         levels(qseqid)))
        output[cbind(sseqid, qseqid)] <- pident
        output
    })
    
    #Fill the blanks
    similarity_matrix[is.na(similarity_matrix)] <- 0
    
    #Mirror it
    similarity_matrix <- similarity_matrix + t(similarity_matrix)
    similarity_matrix[similarity_matrix > 1] <- 1 #I think this was mostly to fix the diagonal
    
    #Make it distances
    similarity_matrix <- 1 - similarity_matrix
    distance_matrix <- similarity_matrix #this seems weird, but it allows me to change the name without copying the full matrix
    rm(similarity_matrix)
    gc()
    
    distance_matrix <- round(distance_matrix, 3)
    distance_matrix <- as.dist(distance_matrix)
    
    #Create the clusters
    cj <- hclust(distance_matrix, method = clustering_method)
    
    cutGroups <- cutree(cj, h = 1 - min_pident_to_keep_match) #1 - similarity
    similarity_clusters_data <- data.table(protein = names(cutGroups), group = cutGroups)
    
    #Add back the proteins_without_hits in their own new groups
    max_group <- max(similarity_clusters_data$group)
    similarity_clusters_data <- rbindlist(list(similarity_clusters_data,
                                               data.table(protein = proteins_without_hits,
                                                          group = max_group + c(1:length(proteins_without_hits)))))
    
    similarity_clusters_data
}

##################-
#### **MAIN** ####
##################-
###################-
#### Read data ####
###################-
BLAST_protein_information <- fread(BLAST_protein_information_file, header = T, sep = "\t", na.strings = NULL)

for (organism_for in organisms) {
    #Protein data
    normalized_protein_data_file <- sprintf("%s/%s/normalized_output_per_protein.tsv", output_data_folder, organism_for)
    normalized_protein_data_for <- fread(normalized_protein_data_file, header = T, sep = "\t", na.strings = NULL)
    
    normalized_protein_data_for$organism <- organism_for
    
    if (organism_for == organisms[1]) {
        normalized_protein_data <- normalized_protein_data_for
    } else {
        normalized_protein_data <- rbindlist(list(normalized_protein_data, normalized_protein_data_for))
    }
    
    #BLAST data
    BLAST_data_file <- sprintf("%s/%s_BLAST.tsv", BLAST_output_data_folder, organism_for)
    BLAST_data_for <- fread(BLAST_data_file, header = T, sep = "\t", na.strings = NULL)
    
    if (organism_for == organisms[1]) {
        BLAST_data <- BLAST_data_for
    } else {
        BLAST_data <- rbindlist(list(BLAST_data, BLAST_data_for))
    }
    
    #Antigen data
    antigenic_tag_protein_data_file <- sprintf("%s/%s_antigens_proteins", antigenic_tags_folder, organism_for)
    antigenic_tag_protein_data_for <- fread(antigenic_tag_protein_data_file, header = T, sep = "\t", na.strings = NULL)
    
    antigenic_tag_protein_data_for$organism <- organism_for
    
    if (organism_for == organisms[1]) {
        antigenic_tag_protein_data <- antigenic_tag_protein_data_for
    } else {
        antigenic_tag_protein_data <- rbindlist(list(antigenic_tag_protein_data, antigenic_tag_protein_data_for))
    }
}
rm(normalized_protein_data_for)
rm(BLAST_data_for)
rm(antigenic_tag_protein_data_for)
gc()

####################################-
#### Create similarity clusters ####
####################################-
#Add the BLAST id to the protein data (everything should match, the all.x it's there just in case)
normalized_protein_data <- merge(normalized_protein_data,
                                 BLAST_protein_information[, .(organism, id = inProteome_id, BLAST_id = id)],
                                 by = c("organism", "id"),
                                 all.x = T)

#Keep only the BLASTS hits inside the organisms in the model
proteins_in_model <- unique(normalized_protein_data$BLAST_id)
BLAST_data <- BLAST_data[(qseqid %in% proteins_in_model) & (sseqid %in% proteins_in_model)]

#Filter out the bad matchs (removes some matches vs itself because the evalue)
BLAST_data$pident <- BLAST_data$pident / 100 #so it's between 0 and 1
BLAST_data <- BLAST_data[(pident >= min_pident_to_keep_match) & (evalue <= max_evalue_to_keep_match)]

#Remove the matches against the same protein (because in the end I want to expand antigenic tags between proteins)
BLAST_data <- BLAST_data[qseqid != sseqid]

#Add the sequence length to the BLAST data
BLAST_data <- merge(BLAST_data,
                    BLAST_protein_information[, .(qseqid = id, qseqlength = nchar(sequence))],
                    by = "qseqid")
BLAST_data <- merge(BLAST_data,
                    BLAST_protein_information[, .(sseqid = id, sseqlength = nchar(sequence))],
                    by = "sseqid")

#Keep only matches where the match length is a given proportion of at least one of the sequences' length
BLAST_data <- BLAST_data[((length / qseqlength) >= min_match_length_proportion_to_keep_match) | ((length / sseqlength) >= min_match_length_proportion_to_keep_match)]

#Apparently some matches are found more than once in the same direction, remove the repeats
BLAST_data <- unique(BLAST_data)

#Create the similarity clusters
similarity_clusters_data <- BLAST_to_similarityClusters(BLAST_data = BLAST_data,
                                                        normalized_protein_data = normalized_protein_data,
                                                        clustering_method = "single",
                                                        min_pident_to_keep_match = min_pident_to_keep_match)

#Add the groups to the proteins and peptides (again, this should match all, the all.x just in case)
normalized_protein_data <- merge(normalized_protein_data,
                                 similarity_clusters_data[, .(BLAST_id = protein, similarity_cluster = group)],
                                 by = "BLAST_id",
                                 all.x = TRUE)

###########################################-
#### Assign positive and negative tags ####
###########################################-
#Add the antigenic tag to the protein data
normalized_protein_data <- merge(normalized_protein_data,
                                 antigenic_tag_protein_data[, .(organism, original_id = protein, is_antigenic = 1)],
                                 by = c("organism", "original_id"),
                                 all.x = T)
normalized_protein_data[is.na(is_antigenic), is_antigenic := 0]
rm(antigenic_tag_protein_data)
gc()

#Expand the antigenic tags according to similarity groups
positive_clusters <- unique(normalized_protein_data[is_antigenic == 1]$similarity_cluster)
negative_clusters <- setdiff(unique(normalized_protein_data$similarity_cluster), positive_clusters)

normalized_protein_data$in_antigenic_cluster <- 0
normalized_protein_data[similarity_cluster %in% positive_clusters, in_antigenic_cluster := 1]

#######################-
#### Save the tags ####
#######################-
for (organism_for in organisms) {
    # organism_for <- organisms[1]
    sub_normalized_protein_data <- normalized_protein_data[organism == organism_for]
    
    #Calculate the output file
    organism_output_data_folder <- sprintf("%s/%s", output_data_folder, organism_for)
    protein_antigenic_tag_file <- sprintf("%s/%s", organism_output_data_folder, protein_antigenic_tag_file_name)
    
    #Keep the important columns
    sub_normalized_protein_data <- sub_normalized_protein_data[, .(id, is_antigenic, similarity_cluster, in_antigenic_cluster)]
    
    #Save
    write.table(sub_normalized_protein_data, file = protein_antigenic_tag_file, col.names = T, row.names = F, sep = "\t", quote = T)
}
