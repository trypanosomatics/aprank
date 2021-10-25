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
if (!require(foreach, quietly = TRUE)) {
    writeLines("Installing library 'foreach' for R")
    install.packages("foreach", repos = "http://cran.rstudio.com/", dependencies = T)
    library(foreach) #parallel processing
}
if (!require(parallel, quietly = TRUE)) {
    writeLines("Installing library 'parallel' for R")
    install.packages("parallel", repos = "http://cran.rstudio.com/", dependencies = T)
    library(parallel) #parallel processing
}
if (!require(doParallel, quietly = TRUE)) {
    writeLines("Installing library 'doParallel' for R")
    install.packages("doParallel", repos = "http://cran.rstudio.com/", dependencies = T)
    library(doParallel) #parallel processing
}

###########################-
#### **AUX FUNCTIONS** ####
###########################-
runAndParse_BepiPred <- function(input_fasta_files,
                                 number_of_parallel_processes = 1) {
    numCores <- detectCores() #set the cores
    
    registerDoParallel(numCores)  #use multicore, set to the number of our cores
    
    full_BepiPred_parsed_data <- foreach (i = 1:number_of_parallel_processes, .combine = rbind) %dopar% {
        command_aux <- sprintf("bepipred '%s' -k",
                               input_fasta_files[i])
        BepiPred_output <- system(command_aux, intern = T)
        
        BepiPred_parsed_data <- read.delim(text = BepiPred_output, header = F, sep = "", na.strings = NULL, comment.char = "#")
        colnames(BepiPred_parsed_data) <- c("seqname", "source", "feature", "start", "end", "score", "filler1", "filler2", "is_epitope")
        BepiPred_parsed_data <- as.data.table(BepiPred_parsed_data)
        BepiPred_parsed_data <- BepiPred_parsed_data[, .(id = seqname, pos = start, score)]
        
        BepiPred_parsed_data
    }
    #End the parallel processing
    stopImplicitCluster()
    
    full_BepiPred_parsed_data
}

runAndParse_IP_MW <- function(input_fasta_files,
                              temp_data_folder,
                              number_of_parallel_processes = 1) {
    numCores <- detectCores() #set the cores
    
    registerDoParallel(numCores)  #use multicore, set to the number of our cores
    
    full_PW_MW_parsed_data <- foreach (i = 1:number_of_parallel_processes, .combine = rbind) %dopar% {
        # i <- 1
        IP_MW_output_file <- sprintf("%s/IP_MW_output_file_%s", temp_data_folder, i)
        
        command_aux <- sprintf("pepstats -sequence '%s' -sprotein1 -aadata Eamino.dat -mwdata Emolwt.dat -termini -nomono -auto -outfile '%s'",
                               input_fasta_files[i],
                               IP_MW_output_file)
        system(command_aux)
        
        #Extract important data
        proteins <- c()
        molecular_weights <- c()
        isoelectric_points <- c()
        con <- file(IP_MW_output_file, "r")
        while(TRUE) {
            line = readLines(con, 1)

            if (length(line) == 0) {
                break #to exit the loop
            } else if (grepl("PEPSTATS of", line)) {
                # PEPSTATS of protein_1 from 1 to 609
                protein <- gsub("^PEPSTATS of (.+) from .*$", "\\1", line)
                proteins <- c(proteins, protein)
            } else if (grepl("Molecular weight", line)) {
                molecular_weight <- gsub("^.*Molecular weight\\s+=\\s+(\\S+)\\s+Residues.*$", "\\1", line)
                molecular_weight <- as.numeric(molecular_weight)
                molecular_weights <- c(molecular_weights, molecular_weight)
            } else if (grepl("Average Residue Weight.*Charge", line)) {
                charge <- gsub("^.*Average Residue Weight.*Charge\\s+=\\s+(\\S+)\\s*$", "\\1", line)
                charge <- as.numeric(charge)
            } else if (grepl("Isoelectric Point", line)) {
                isoelectric_point <- gsub("^.*Isoelectric Point\\s+=\\s+(\\S+)\\s*$", "\\1", line)
                #In some strange cases IP is "None", set as 0 or as 14 based on the charge
                if (isoelectric_point != "None") {
                    isoelectric_point <- as.numeric(isoelectric_point)    
                } else {
                    if (charge > 0) {
                        isoelectric_point <- 14
                    } else {
                        isoelectric_point <- 0
                    }
                }
                isoelectric_points <- c(isoelectric_points, isoelectric_point)
            }
        }
        close(con)

        IP_MW_parsed_data <- data.table(id = proteins,
                                        molecular_weight = as.numeric(molecular_weights),
                                        isoelectric_point = as.numeric(isoelectric_points))
        
        IP_MW_parsed_data
    }
    #End the parallel processing
    stopImplicitCluster()
    
    full_PW_MW_parsed_data
}

runAndParse_NetSurfp <- function(input_singleProteinFasta_file) {
    #NetSurfp already has parallelization, so it's convenient to use less parallel processes, such as
    #leaving 4 cores free per process
    command_aux <- sprintf("netsurfp -i '%s' -a",
                           input_singleProteinFasta_file)
    NetSurfp_output <- system(command_aux, intern = T)
    
    #There is a problem with strange AA such as X where it doesn't give a class_asignment, meaning
    #there is one less column, let's fix those cases
    for (i in 1:length(NetSurfp_output)) {
        # i <- 1
        line <- NetSurfp_output[i]
        
        if ((line != "") & (!grepl("^#.*$", line))) {
            if (grepl("^\\s+[^ACDEFGHIKLMNPQRSTVWY]\\s+.*$", line)) {
                line <- sprintf("X %s", gsub("^\\s+(\\S.*)$", "\\1", line))
                NetSurfp_output[i] <- line
            }
        }
    }
    
    NetSurfp_parsed_data <- read.delim(text = NetSurfp_output, header = F, sep = "", na.strings = NULL, comment.char = "#")
    colnames(NetSurfp_parsed_data) <- c("class_asignment", "aa", "id", "pos", "RSA", "ASA", "Z_fit_score", "prob_of_alpha_helix", "prob_of_beta_strand", "prob_of_coil")
    NetSurfp_parsed_data <- as.data.table(NetSurfp_parsed_data)
    NetSurfp_parsed_data <- NetSurfp_parsed_data[, .(id, pos, relative_surface_accessibility = RSA,
                                                     probability_for_alpha_helix = prob_of_alpha_helix,
                                                     probability_for_beta_strand = prob_of_beta_strand)]
    
    NetSurfp_parsed_data
}

runAndParse_Paircoil2 <- function(input_fasta_files,
                                  temp_data_folder,
                                  aminoacid_fasta,
                                  tabbed_fasta,
                                  number_of_parallel_processes = 1) {
    numCores <- detectCores() #set the cores
    
    registerDoParallel(numCores)  #use multicore, set to the number of our cores
    
    full_Paircoil_parsed_data <- foreach (i = 1:number_of_parallel_processes, .combine = rbind) %dopar% {
        # i <- 1
        Paircoil_output_file <- sprintf("%s/Paircoil_output_%s_file", temp_data_folder, i)
        Paircoil_error1_file <- sprintf("%s/Paircoil_error1_%s_file", temp_data_folder, i)
        Paircoil_error2_file <- sprintf("%s/Paircoil_error2_%s_file", temp_data_folder, i)
        
        command_aux <- sprintf("paircoil2 '%s' '%s' >>'%s' 2>>'%s'",
                               input_fasta_files[i],
                               Paircoil_output_file,
                               Paircoil_error1_file,
                               Paircoil_error2_file)
        system(command_aux)
        
        Paircoil_parsed_data <- read.delim(file = Paircoil_output_file, header = F, sep = "", na.strings = NULL, comment.char = "#")
        colnames(Paircoil_parsed_data) <- c("position", "residue", "register", "P-score", "parenthesis", "score")
        Paircoil_parsed_data <- as.data.table(Paircoil_parsed_data)
        
        #Paircoil doesn't process all proteins (they have to be at least 28aa), so I'll look for the proteins that it actually processed
        con <- file(Paircoil_output_file, "r")
        Paircoil_proteins <- c()
        while(TRUE) {
            line = readLines(con, 1)
            
            if (length(line) == 0) {
                break #to exit the loop
            } else if (grepl("^# Sequence Code:", line)) {
                Paircoil_proteins <- c(Paircoil_proteins,
                                       gsub("^# Sequence Code: (.+)$", "\\1", line))
            }
        }
        close(con)
        
        #Now, I'll add these protein ids to the data and keep what I care about
        Paircoil_parsed_data$id <- rep(Paircoil_proteins, tabbed_fasta[protein %in% Paircoil_proteins]$sequence_length)
        Paircoil_parsed_data <- Paircoil_parsed_data[, .(id, pos = position, p_score = `P-score`)]
        
        Paircoil_parsed_data
    }
    #End the parallel processing
    stopImplicitCluster()
    
    #And I'll add any protein that wasn't processed
    Paircoil_other_proteins <- aminoacid_fasta[!(protein %in% unique(full_Paircoil_parsed_data$id))]
    if (Paircoil_other_proteins[, .N] > 0) {
        Paircoil_other_proteins$p_score <- 0
        
        full_Paircoil_parsed_data <- rbindlist(list(full_Paircoil_parsed_data,
                                               Paircoil_other_proteins[, .(id = protein, pos, p_score)]))
    }
    
    full_Paircoil_parsed_data
}

runAndParse_PredGPI <- function(temp_data_folder,
                                tabbed_fasta,
                                PredGPI_min_length = 41,
                                number_of_parallel_processes = 1) {
    numCores <- detectCores() #set the cores
    
    registerDoParallel(numCores)  #use multicore, set to the number of our cores
    
    full_PredGPI_parsed_data <- foreach (i = 1:number_of_parallel_processes, .combine = rbind) %dopar% {
        # i <- 4
        # PredGPI_min_length <- 41
        PredGPI_fasta_file <- sprintf("%s/PredGPI_temp_%s.fasta", temp_data_folder, i)
        
        #PredGPI needs a clean FASTA with long sequences to work, so let's do that
        PredGPI_tabbed_fasta <- tabbed_fasta[(parallel_process_i == i) & (sequence_length >= PredGPI_min_length) & (grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", sequence))]
        PredGPI_tabbed_fasta[, output_aux := sprintf(">%s\n%s", protein, sequence)]
        PredGPI_output_aux <- paste(PredGPI_tabbed_fasta$output_aux, collapse = "\n")
        write(PredGPI_output_aux, file = PredGPI_fasta_file)
        
        command_aux <- sprintf("PredGPI.py '%s'",
                               PredGPI_fasta_file)
        PredGPI_output <- system(command_aux, intern = T)
        
        #Each proteins gives two rows, the id and the data so I have to parse that
        PredGPI_parsed_data <- read.delim(text = PredGPI_output, header = F, sep = "\t", na.strings = NULL, comment.char = "")
        PredGPI_parsed_data_aux <- PredGPI_parsed_data[c(F, T),]
        PredGPI_parsed_data_aux$id <- PredGPI_parsed_data[c(T, F),]$V1
        PredGPI_parsed_data <- PredGPI_parsed_data_aux
        colnames(PredGPI_parsed_data) <- c("has_gpi_detailed", "score", "gpi_start_string", "id")
        PredGPI_parsed_data <- as.data.table(PredGPI_parsed_data)
        
        #The start is the first part of the string "586-S" (for example) so I'll extract that
        PredGPI_parsed_data[, gpi_start := as.numeric(gsub("^(.+)-.+", "\\1", gpi_start_string))]
        
        #has_dpi_detailed can be "GPI Highly probable", "GPI Lowly probable" and "Non-GPI" (among others); I'll transform that to numbers
        PredGPI_parsed_data$has_gpi <- 1
        PredGPI_parsed_data[has_gpi_detailed == "Non-GPI", has_gpi := 0]
        PredGPI_parsed_data[has_gpi_detailed == "GPI Lowly probable", has_gpi := 0]
        
        PredGPI_parsed_data <- PredGPI_parsed_data[, .(id, has_gpi, score, gpi_start)]
        
        PredGPI_parsed_data
    }
    #End the parallel processing
    stopImplicitCluster()
    
    #Add any protein that wasn't analyzed
    PredGPI_proteins <- full_PredGPI_parsed_data$id
    PredGPI_other_proteins <- tabbed_fasta[!(protein %in% PredGPI_proteins)]
    if (PredGPI_other_proteins[, .N] > 0) {
        PredGPI_other_proteins$has_gpi <- 0
        PredGPI_other_proteins$score <- 0
        PredGPI_other_proteins$gpi_start <- 0
        
        full_PredGPI_parsed_data <- rbindlist(list(full_PredGPI_parsed_data,
                                                   PredGPI_other_proteins[, .(id = protein, has_gpi, score, gpi_start)]))
    }
    
    #If it doesn't have gpi then set the start and the score as 0
    full_PredGPI_parsed_data[has_gpi == 0, score := 0]
    full_PredGPI_parsed_data[has_gpi == 0, gpi_start := 0]
    
    full_PredGPI_parsed_data
}

runAndParse_SignalP <- function(input_fasta_files,
                                SignalP_organism_group,
                                aminoacid_fasta,
                                number_of_parallel_processes = 1) {
    numCores <- detectCores() #set the cores
    
    registerDoParallel(numCores)  #use multicore, set to the number of our cores
    
    full_list_output <- foreach (i = 1:number_of_parallel_processes, .combine = c) %dopar% {
        command_aux <- sprintf("signalp -f long -t %s '%s'",
                               SignalP_organism_group,
                               input_fasta_files[i])
        SignalP_output <- system(command_aux, intern = T)
        
        #SignalP output is a mess with 3 different formats, extract the important data
        SignalP_protein_parsed_data <- data.table(id = character(),
                                                  has_signal_peptide = integer(),
                                                  cleavage_site = integer())
        SignalP_aminoacid_parsed_data <- data.table(id = character(),
                                                    pos = integer(),
                                                    C_score = numeric(),
                                                    S_score = numeric())
        SignalP_protein_aminoacids_data <- c()
        for (line in SignalP_output) {
            if (grepl("^\\s*\\d", line)) {
                #If it starts with a number (aa data)
                SignalP_protein_aminoacids_data <- c(SignalP_protein_aminoacids_data,
                                                     line)
            } else if (grepl("^\\s*Name=", line)) {
                #These are the lines to extract protein data from
                id_aux <- gsub("^\\s*Name=(\\S+)\\t.*$", "\\1", line)
                
                SignalP_protein_row_aux <- data.table(id = id_aux,
                                                      has_signal_peptide = 0,
                                                      cleavage_site = 0)
                if (grepl("SP=.YES", line)) {
                    cleavage_site_aux <- gsub("^.*between pos\\. (\\d+) and.*$", "\\1", line)
                    
                    SignalP_protein_row_aux$has_signal_peptide <- 1
                    SignalP_protein_row_aux$cleavage_site <- as.numeric(cleavage_site_aux)
                }
                SignalP_protein_parsed_data <- rbindlist(list(SignalP_protein_parsed_data,
                                                              SignalP_protein_row_aux))
                
                #This line appears at the end of each protein, so save the aminoacid data for that protein
                SignalP_protein_aminoacids_parsed_data <- read.delim(text = SignalP_protein_aminoacids_data, header = F, sep = "", na.strings = NULL, comment.char = "#")
                SignalP_protein_aminoacids_parsed_data <- as.data.table(SignalP_protein_aminoacids_parsed_data)
                colnames(SignalP_protein_aminoacids_parsed_data) <- c("pos", "aa", "C_score", "S_score", "Y_score")
                
                SignalP_protein_aminoacids_parsed_data$id <- id_aux
                SignalP_protein_aminoacids_parsed_data <- SignalP_protein_aminoacids_parsed_data[, .(id, pos, C_score, S_score)]
                
                SignalP_aminoacid_parsed_data <- rbindlist(list(SignalP_aminoacid_parsed_data,
                                                                SignalP_protein_aminoacids_parsed_data))
                
                SignalP_protein_aminoacids_data <- c()
            }
        }
        
        #Create the output
        list_output <- list()
        list_output[[sprintf("SignalP_protein_parsed_data_%s", i)]] <- SignalP_protein_parsed_data
        list_output[[sprintf("SignalP_aminoacid_parsed_data_%s", i)]] <- SignalP_aminoacid_parsed_data
        
        list_output
    }
    #End the parallel processing
    stopImplicitCluster()
    
    #Paste the outputs together
    for (i in 1:number_of_parallel_processes) {
        #i <- 2
        if (i == 1) {
            full_SignalP_protein_parsed_data <- full_list_output[[sprintf("SignalP_protein_parsed_data_%s", i)]]
            full_SignalP_aminoacid_parsed_data <- full_list_output[[sprintf("SignalP_aminoacid_parsed_data_%s", i)]]
        } else {
            full_SignalP_protein_parsed_data <- rbindlist(list(full_SignalP_protein_parsed_data,
                                                               full_list_output[[sprintf("SignalP_protein_parsed_data_%s", i)]]))
            full_SignalP_aminoacid_parsed_data <- rbindlist(list(full_SignalP_aminoacid_parsed_data,
                                                                 full_list_output[[sprintf("SignalP_aminoacid_parsed_data_%s", i)]]))
        }
    }

    #Add the missing positions with score 0
    full_SignalP_aminoacid_parsed_data <- merge(aminoacid_fasta[, .(id = protein, pos)],
                                                full_SignalP_aminoacid_parsed_data,
                                                by = c("id", "pos"),
                                                all.x = T)
    full_SignalP_aminoacid_parsed_data[is.na(C_score), C_score := 0]
    full_SignalP_aminoacid_parsed_data[is.na(S_score), S_score := 0]
    

    
    #Create the output
    list_output <- list()
    list_output[["SignalP_protein_parsed_data"]] <- full_SignalP_protein_parsed_data
    list_output[["SignalP_aminoacid_parsed_data"]] <- full_SignalP_aminoacid_parsed_data
    
    list_output
}

runAndParse_TMHMM <- function(input_fasta_files,
                              aminoacid_fasta,
                              number_of_parallel_processes = 1) {
    numCores <- detectCores() #set the cores
    
    registerDoParallel(numCores)  #use multicore, set to the number of our cores
    
    full_TMHMM_protein_parsed_data <- foreach (i = 1:number_of_parallel_processes, .combine = rbind) %dopar% {
        command_aux <- sprintf("tmhmm '%s'",
                               input_fasta_files[i])
        TMHMM_output <- system(command_aux, intern = T)
        
        #Extract important data
        TMHMM_data <- c()
        for (line in TMHMM_output) {
            if (grepl("^.*Topology=.*$", line)) {
                TMHMM_data <- c(TMHMM_data,
                                line)
            }
        }
        TMHMM_protein_parsed_data <- read.delim(text = TMHMM_data, header = F, sep = "\t", na.strings = NULL, comment.char = "")
        TMHMM_protein_parsed_data <- as.data.table(TMHMM_protein_parsed_data)
        colnames(TMHMM_protein_parsed_data) <- c("id", "len", "expAA", "first60", "predHel", "topology")
        
        #If the whole sequence is labeled as inside or outside, the prediction is that it contains no membrane helices
        TMHMM_protein_parsed_data[, topology := gsub("Topology=", "", topology)]
        TMHMM_protein_parsed_data$has_transmembrane_helix <- 1
        TMHMM_protein_parsed_data[topology == "i", has_transmembrane_helix := 0]
        TMHMM_protein_parsed_data[topology == "o", has_transmembrane_helix := 0]
        
        TMHMM_protein_parsed_data <- TMHMM_protein_parsed_data[, .(id, has_transmembrane_helix, topology)]
        
        TMHMM_protein_parsed_data
    }
    #End the parallel processing
    stopImplicitCluster()
    
    #The topology indicates which aa are inside transmembrane helixes, calculate the aa data
    full_TMHMM_aminoacid_parsed_data <- aminoacid_fasta[, .(id = protein, pos)]
    full_TMHMM_aminoacid_parsed_data$in_transmembrane_helix <- 0
    
    proteins_with_helix <- full_TMHMM_protein_parsed_data[has_transmembrane_helix == 1]
    if (proteins_with_helix[, .N] > 0) {
        for (i in 1:proteins_with_helix[, .N]) {
            #i <- 1
            protein_for <- proteins_with_helix[i]
            topology_for <- protein_for$topology
            
            if (topology_for == "i") {
                full_TMHMM_aminoacid_parsed_data[id == protein_for$id, in_transmembrane_helix := 1]
            } else {
                #Topology example: Topology=i7-29o44-66i87-109o (a set of ranges with i and o outside of them)
                # topology_for <- "i7-29o44-66i87-109o"
                ranges_aux <- unlist(strsplit(topology_for, "[i|o]"))
                
                positions_in_helix <- c()
                for (range_for in ranges_aux) {
                    # range_for <- ranges_aux[2]
                    if (range_for != "") {
                        ranges_borders_aux <- unlist(strsplit(range_for, "-"))
                        min_border_aux <- as.numeric(ranges_borders_aux[1])
                        max_border_aux <- as.numeric(ranges_borders_aux[2])
                        
                        positions_in_helix <- c(positions_in_helix,
                                                c(min_border_aux:max_border_aux))
                    }
                }
                
                full_TMHMM_aminoacid_parsed_data[(id == protein_for$id) & (pos %in% positions_in_helix),
                                            in_transmembrane_helix := 1]
            }
        }
    }
    
    #Create the output
    list_output <- list()
    list_output[["TMHMM_protein_parsed_data"]] <- full_TMHMM_protein_parsed_data
    list_output[["TMHMM_aminoacid_parsed_data"]] <- full_TMHMM_aminoacid_parsed_data
    
    list_output
}

runAndParse_Iupred <- function(input_singleProteinFasta_file,
                               protein_par) {
    command_aux <- sprintf("iupred '%s' short",
                           input_singleProteinFasta_file)
    Iupred_output <- system(command_aux, intern = T)
    
    Iupred_parsed_data <- read.delim(text = Iupred_output, header = F, sep = "", na.strings = NULL, comment.char = "#")
    colnames(Iupred_parsed_data) <- c("pos", "aa", "disorder")
    Iupred_parsed_data <- as.data.table(Iupred_parsed_data)
    Iupred_parsed_data$id <- protein_par
    Iupred_parsed_data <- Iupred_parsed_data[, .(id, pos, disorder)]
    
    Iupred_parsed_data
}

runAndParse_NetMHCIIpan <- function(input_singleProteinFasta_file,
                                    NetMHCIIpan_alleles,
                                    NetMHCIIpan_binding_peptide_length,
                                    list_output_suffix = "") {
    #Run the program for each allele and merge the data in one file
    for (allele_for in NetMHCIIpan_alleles) {
        # allele_for <- NetMHCIIpan_alleles[1]
        command_aux <- sprintf("netMHCIIpan -a '%s' -f '%s' -l %s",
                               allele_for,
                               input_singleProteinFasta_file,
                               NetMHCIIpan_binding_peptide_length)
        NetMHCIIpan_output_aux <- system(command_aux, intern = T)
        
        if (allele_for == NetMHCIIpan_alleles[1]) {
            NetMHCIIpan_output <- NetMHCIIpan_output_aux
        } else {
            NetMHCIIpan_output <- c(NetMHCIIpan_output, NetMHCIIpan_output_aux)
        }
    }
    
    #Read the important lines of the file
    NetMHCIIpan_allele_data <- c()
    for (line in NetMHCIIpan_output) {
        if (grepl("^\\s*\\d", line)) {
            #Some lines have a <= at the end, remove it
            line_aux <- gsub("<= .*$", "", line)
            #Record the line
            NetMHCIIpan_allele_data <- c(NetMHCIIpan_allele_data,
                                         line_aux)
        }
    }
    
    NetMHCIIpan_parsed_allele_data <- read.delim(text = NetMHCIIpan_allele_data, header = F, sep = "", na.strings = NULL, comment.char = "")
    NetMHCIIpan_parsed_allele_data <- as.data.table(NetMHCIIpan_parsed_allele_data)
    colnames(NetMHCIIpan_parsed_allele_data) <- c("pos", "HLA", "peptide", "id", "other_pos", "core", "aff", "affinity_nM", "rank")
    
    NetMHCIIpan_parsed_allele_data[, start := pos + 1]
    
    #Get data by peptide
    NetMHCIIpan_parsed_peptide_data <- NetMHCIIpan_parsed_allele_data[, .(min_rank = min(rank),
                                                                          max_rank = max(rank),
                                                                          mean_rank = mean(rank)),
                                                                      by = .(id, peptide, start)]
    
    #Create the output
    list_output <- list()
    list_output[[sprintf("NetMHCIIpan_parsed_allele_data%s", list_output_suffix)]] <- NetMHCIIpan_parsed_allele_data
    list_output[[sprintf("NetMHCIIpan_parsed_peptide_data%s", list_output_suffix)]] <- NetMHCIIpan_parsed_peptide_data
    
    list_output
}

runAndParse_NetOglyc <- function(input_singleProteinFasta_file,
                                 protein_par,
                                 aminoacid_fasta) {
    ### OJO, el input ID tiene que tener menos de 12 caracteres...
    command_aux <- sprintf("netOglyc '%s'",
                           input_singleProteinFasta_file)
    
    #This gives an "seq2seq: empty input line, not allowed", but the output it's correct, so I'll ignore it
    NetOglyc_output <- system(command_aux, intern = T, ignore.stderr = T)
    
    #Extract important data
    NetOglyc_data <- c()
    in_data_aux <- 0
    for (line in NetOglyc_output) {
        if (length(line) == 0) {
            break #to exit the loop
        } else if (grepl("Name\\s+S/T\\s+", line)) {
            in_data_aux <- 1
        } else if ((in_data_aux == 1) & grepl("\\s+[ST]\\s+", line)) {
            NetOglyc_data <- c(NetOglyc_data, line)
        }
    }
    
    #Sometimes (for short sequences?) the data comes back empty
    if (length(NetOglyc_data) > 0) {
        NetOglyc_parsed_data <- read.delim(text = NetOglyc_data, header = F, sep = "", na.strings = NULL, comment.char = "")
        NetOglyc_parsed_data <- as.data.table(NetOglyc_parsed_data)
        colnames(NetOglyc_parsed_data) <- c("id", "aa", "pos", "G_score", "I_score", "glycosylation_site_aa", "comments")
        
        NetOglyc_parsed_data$is_glycosylation_site <- 0
        NetOglyc_parsed_data[glycosylation_site_aa %in% c("S", "T"), is_glycosylation_site := 1]
    } else {
        NetOglyc_parsed_data <- data.table(id = character(),
                                           pos = integer(),
                                           is_glycosylation_site = integer())
    }
    
    #Add the missing positions with score 0
    NetOglyc_parsed_data <- merge(aminoacid_fasta[protein == protein_par, .(id = protein, pos)],
                                  NetOglyc_parsed_data,
                                  by = c("id", "pos"),
                                  all.x = T)
    NetOglyc_parsed_data[is.na(is_glycosylation_site), is_glycosylation_site := 0]
    
    NetOglyc_parsed_data <- NetOglyc_parsed_data[, .(id, pos, is_glycosylation_site)]
    
    NetOglyc_parsed_data
}

runAndParse_Xstream <- function(input_singleProteinFasta_file,
                                protein,
                                Xstream_temp_folder,
                                Xstream_path) {
    #I can't change the name of the output, but I can change it's location
    Xstream_output_file <- sprintf("%s/\\XSTREAM__i0.7_g3_m5_e2.0_out_2.html", Xstream_temp_folder)
    
    #I don't care about the log, but otherwise it's shown on console
    command_aux <- sprintf("java -jar '%s' '%s' -d'%s/' >'%s/xstream.log'",
                           Xstream_path,
                           input_singleProteinFasta_file,
                           Xstream_temp_folder,
                           Xstream_temp_folder)
    system(command_aux)
    
    #Extract important data
    #There is a warning about the final line not ending properly, but that's Xstream's problem and it doesn't affect this
    con <- file(Xstream_output_file, "r")
    Xstream_data <- c()
    while(TRUE) {
        line = suppressWarnings(readLines(con, 1)) #To hide the final line warning
        
        if (length(line) == 0) {
            break #to exit the loop
        } else if (grepl("No Repeats Found!", line)) {
            break #this means there are no repeats in this protein
        } else if (grepl("Positions", line)) {
            Xstream_data <- c(Xstream_data,
                              line)
        }
    }
    close(con)
    
    if (length(Xstream_data) > 0) {
        for (repeat_for in Xstream_data) {
            # repeat_for <- Xstream_data[1]
            general_pattern <- "^.*<FONT COLOR=\"00BFFF\">(\\d+)-(\\d+)</FONT></CENTER></TD><TD WIDTH=80><CENTER><FONT COLOR=\\\"00BFFF\\\">(\\d+)</FONT></CENTER></TD><TD WIDTH=80><CENTER><FONT COLOR=\\\"00BFFF\\\">([^<]+)</FONT></CENTER></TD><TD WIDTH=80><CENTER><FONT COLOR=\\\"00BFFF\\\">([^<]+)</FONT>.*$"
            start_aux <- gsub(general_pattern, "\\1", repeat_for)
            end_aux <- gsub(general_pattern, "\\2", repeat_for)
            period_aux <- gsub(general_pattern, "\\3", repeat_for)
            copy_number_aux <- gsub(general_pattern, "\\4", repeat_for)
            consensus_error_aux <- gsub(general_pattern, "\\5", repeat_for)
            
            Xstream_row_data <- data.table(id = protein,
                                           start = as.integer(start_aux),
                                           end = as.integer(end_aux),
                                           period = as.numeric(period_aux),
                                           copy_number = as.numeric(copy_number_aux),
                                           consensus_error = as.numeric(consensus_error_aux))
            
            if (repeat_for == Xstream_data[1]) {
                Xstream_parsed_data <- Xstream_row_data
            } else {
                Xstream_parsed_data <- rbindlist(list(Xstream_parsed_data, Xstream_row_data))
            }
        }
    } else {
        #There are no repeats, return empty datatable
        # Xstream_parsed_data <- data.table(id = character(),
        #                                   start = integer(),
        #                                   end = integer(),
        #                                   period = numeric(),
        #                                   copy_number = numeric(),
        #                                   consensus_error = numeric())
        Xstream_parsed_data <- data.table(id = protein,
                                          start = 0,
                                          end = 0,
                                          period = 0,
                                          copy_number = 0,
                                          consensus_error = 0)
    }
    
    Xstream_parsed_data
}

##################-
#### **MAIN** ####
##################-
runPredictors <- function(input_fasta_file,
                          temp_data_folder, output_data_folder, scripts_folder,
                          temp_fasta_file, temp_tabbed_fasta_file, temp_splitted_fasta_file, temp_aminoacid_fasta_file,
                          output_method, output_per_aminoacid_file = "", output_per_protein_file = "", output_per_kmer_file = "", output_per_peptide_file = "", output_per_repeat_file = "",
                          number_of_parallel_processes = 1, number_of_parallel_processes_for_NetSurfp = 1,
                          peptide_length = 15, peptide_overlap = 14, max_protein_length = 9999, replace_nonAA_chars_by = "X",
                          SignalP_organism_group = "euk",
                          Xstream_path = "/usr/local/xstream/xstream.jar",
                          NetMHCIIpan_binding_peptide_length = 9, NetMHCIIpan_alleles = c("DRB1_0101", "DRB3_0101", "DRB4_0101", "DRB5_0101"),
                          KmerSimilarity_kmer_length = 6, CrossReactivity_fasta_file = "", Coendemicity_fasta_file = "",
                          use_BepiPred = 1, use_Paircoil2 = 1, use_PredGPI = 1, use_SignalP = 1,
                          use_TMHMM = 1, use_NetSurfp = 1,
                          use_Iupred = 1, use_NetOglyc = 1, use_Xstream = 1, use_NetMHCIIpan = 1,
                          use_IsoelectricPoint = 1, use_MolecularWeight = 1,
                          use_SelfSimilarity = 1, use_CrossReactivity = 1, use_Coendemicity = 1) {
    #####################-
    #### Parse FASTA ####
    #####################-
    #I limit the length of each protein to a number (9999 by default) and replace non aminoacid characters
    #to X to make it easier for the predictors (there are some that simply cant handle "strange" things or long
    #proteins)
    if (file.exists(temp_tabbed_fasta_file) == F) {
        min_protein_length <- max(c(10, peptide_length))
        
        command_aux <- sprintf("perl '%s/lib/tab_and_split_proteome.pl' '%s' %s %s %s %s '%s' '%s' '%s' '%s' '%s'",
                               scripts_folder,
                               input_fasta_file,
                               peptide_length,
                               peptide_overlap,
                               min_protein_length,
                               max_protein_length,
                               replace_nonAA_chars_by,
                               temp_fasta_file,
                               temp_tabbed_fasta_file,
                               temp_splitted_fasta_file,
                               temp_aminoacid_fasta_file)
        system(command_aux)
    }
    tabbed_fasta <- fread(temp_tabbed_fasta_file, header = T, sep = "\t", na.strings = NULL)
    tabbed_fasta[, sequence_length := nchar(sequence)]
    aminoacid_fasta <- fread(temp_aminoacid_fasta_file, header = T, sep = "\t", na.strings = NULL)
    
    ########################################-
    #### Set up the Parallel Processing ####
    ########################################-
    #Check if I have more processors than proteins
    if (number_of_parallel_processes > tabbed_fasta[, .N]) {
        number_of_parallel_processes <- tabbed_fasta[, .N]
    }
    
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
            # i <- 4
            protein_subset_min_for_parallel <- 1 + sum(protein_per_process_amount[0:i])
            protein_subset_max_for_parallel <- protein_subset_min_for_parallel + protein_per_process_amount[i + 1] - 1
            tabbed_fasta[protein_subset_min_for_parallel:protein_subset_max_for_parallel]$parallel_process_i <- i
            
            sub_temp_fasta_file <- sprintf("%s/sub_temp_%s.fasta", temp_data_folder, i)
            if (file.exists(sub_temp_fasta_file) == F) {
                sub_fasta_output <- paste(tabbed_fasta[protein_subset_min_for_parallel:protein_subset_max_for_parallel]$output_aux, collapse = "\n")
                write(sub_fasta_output, file = sub_temp_fasta_file)    
            }
            
            sub_temp_fastas_files <- c(sub_temp_fastas_files, sub_temp_fasta_file)
        }
        
        tabbed_fasta <- tabbed_fasta[, -c("output_aux")]
    } else {
        sub_temp_fastas_files <- c(temp_fasta_file)
    }
    
    #Assign the parallel process id for NetSurfp
    tabbed_fasta$parallel_process_NetSurfp_i <- 1
    if ((use_NetSurfp) & (number_of_parallel_processes_for_NetSurfp > 1)) {
        unique_protein_amount <- tabbed_fasta[, .N]
        
        division_aux <- unique_protein_amount %/% number_of_parallel_processes_for_NetSurfp
        reminder_aux <- unique_protein_amount %% number_of_parallel_processes_for_NetSurfp
        protein_per_process_amount <-c(0, rep(division_aux, number_of_parallel_processes_for_NetSurfp))
        protein_per_process_amount <- protein_per_process_amount + c(0, rep(1, reminder_aux), rep(0, number_of_parallel_processes_for_NetSurfp - reminder_aux))
        for (i in 1:number_of_parallel_processes_for_NetSurfp) {
            # i <- 4
            protein_subset_min_for_parallel <- 1 + sum(protein_per_process_amount[0:i])
            protein_subset_max_for_parallel <- protein_subset_min_for_parallel + protein_per_process_amount[i + 1] - 1
            tabbed_fasta[protein_subset_min_for_parallel:protein_subset_max_for_parallel]$parallel_process_NetSurfp_i <- i
        }
    }
    
    ###################################-
    #### **RUN ONCE PER PROTEOME** ####
    ###################################-
    ##################-
    #### BepiPred ####
    ##################-
    if (use_BepiPred) {
        BepiPred_temp_file <- sprintf("%s/temp_output_BepiPred.tsv", temp_data_folder)
        if (file.exists(BepiPred_temp_file) == F) {
            writeLines("Running BepiPred...")
            
            BepiPred_parsed_data <- runAndParse_BepiPred(input_fasta_files = sub_temp_fastas_files,
                                                         number_of_parallel_processes = number_of_parallel_processes)    
            
            write.table(BepiPred_parsed_data, file = BepiPred_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
        } else {
            writeLines("Recovering BepiPred Data...")
            
            BepiPred_parsed_data <- fread(file = BepiPred_temp_file, header = T, sep = "\t", na.strings = NULL)
        }
    }
    
    ##############################################-
    #### Isoelectric point & Molecular weight ####
    ##############################################-
    if (use_IsoelectricPoint | use_MolecularWeight) {
        IP_MW_temp_file <- sprintf("%s/temp_output_IP_MW.tsv", temp_data_folder)
        if (file.exists(IP_MW_temp_file) == F) {
            writeLines("Running EMBOSS...")
            
            IP_MW_parsed_data <- runAndParse_IP_MW(input_fasta_files = sub_temp_fastas_files,
                                                   temp_data_folder = temp_data_folder,
                                                   number_of_parallel_processes = number_of_parallel_processes)
            
            write.table(IP_MW_parsed_data, file = IP_MW_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
        } else {
            writeLines("Recovering EMBOSS Data...")
            IP_MW_parsed_data <- fread(file = IP_MW_temp_file, header = T, sep = "\t", na.strings = NULL)
        }
    }
    
    ###################-
    #### Paircoil2 ####
    ###################-
    if (use_Paircoil2) {
        Paircoil_temp_file <- sprintf("%s/temp_output_Paircoil.tsv", temp_data_folder)
        if (file.exists(Paircoil_temp_file) == F) {
            writeLines("Running Paircoil2...")
            
            Paircoil_parsed_data <- runAndParse_Paircoil2(input_fasta_files = sub_temp_fastas_files,
                                                          temp_data_folder = temp_data_folder,
                                                          aminoacid_fasta = aminoacid_fasta,
                                                          tabbed_fasta = tabbed_fasta,
                                                          number_of_parallel_processes = number_of_parallel_processes)
            
            write.table(Paircoil_parsed_data, file = Paircoil_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
        } else {
            writeLines("Recovering Paircoil2 Data...")
            Paircoil_parsed_data <- fread(file = Paircoil_temp_file, header = T, sep = "\t", na.strings = NULL)
        }
    }
    
    #################-
    #### PredGPI ####
    #################-
    if (use_PredGPI) {
        PredGPI_temp_file <- sprintf("%s/temp_output_PredGPI.tsv", temp_data_folder)
        if (file.exists(PredGPI_temp_file) == F) {
            writeLines("Running PredGPI...")
            
            PredGPI_parsed_data <- runAndParse_PredGPI(temp_data_folder = temp_data_folder,
                                                       tabbed_fasta = tabbed_fasta,
                                                       PredGPI_min_length = 41,
                                                       number_of_parallel_processes = number_of_parallel_processes)  
            
            write.table(PredGPI_parsed_data, file = PredGPI_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
        } else {
            writeLines("Recovering PredGPI Data...")
            PredGPI_parsed_data <- fread(file = PredGPI_temp_file, header = T, sep = "\t", na.strings = NULL)
        }
    }
    
    #################-
    #### SignalP ####
    #################-
    if (use_SignalP) {
        SignalP_proteins_temp_file <- sprintf("%s/temp_output_SignalP_proteins.tsv", temp_data_folder)
        SignalP_aminoacids_temp_file <- sprintf("%s/temp_output_SignalP_aminoacids.tsv", temp_data_folder)
        if (file.exists(SignalP_proteins_temp_file) == F) {
            writeLines("Running SignalP...")
            
            list_output <- runAndParse_SignalP(input_fasta_files = sub_temp_fastas_files,
                                               SignalP_organism_group = SignalP_organism_group,
                                               aminoacid_fasta = aminoacid_fasta,
                                               number_of_parallel_processes = number_of_parallel_processes)
            SignalP_protein_parsed_data <- list_output[["SignalP_protein_parsed_data"]]
            SignalP_aminoacid_parsed_data <- list_output[["SignalP_aminoacid_parsed_data"]] 
            
            write.table(SignalP_protein_parsed_data, file = SignalP_proteins_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
            write.table(SignalP_aminoacid_parsed_data, file = SignalP_aminoacids_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
        } else {
            writeLines("Recovering SignalP Data...")
            
            SignalP_protein_parsed_data <- fread(file = SignalP_proteins_temp_file, header = T, sep = "\t", na.strings = NULL)
            SignalP_aminoacid_parsed_data <- fread(file = SignalP_aminoacids_temp_file, header = T, sep = "\t", na.strings = NULL)
        }
    }
    
    ###############-
    #### TMHMM ####
    ###############-
    if (use_TMHMM) {
        TMHMM_proteins_temp_file <- sprintf("%s/temp_output_TMHMM_proteins.tsv", temp_data_folder)
        TMHMM_aminoacids_temp_file <- sprintf("%s/temp_output_TMHMM_aminoacids.tsv", temp_data_folder)
        if (file.exists(TMHMM_proteins_temp_file) == F) {
            writeLines("Running TMHMM...")
            
            list_output <- runAndParse_TMHMM(input_fasta_files = sub_temp_fastas_files,
                                             aminoacid_fasta = aminoacid_fasta,
                                             number_of_parallel_processes = number_of_parallel_processes)
            TMHMM_protein_parsed_data <- list_output[["TMHMM_protein_parsed_data"]]
            TMHMM_aminoacid_parsed_data <- list_output[["TMHMM_aminoacid_parsed_data"]]
            
            write.table(TMHMM_protein_parsed_data, file = TMHMM_proteins_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
            write.table(TMHMM_aminoacid_parsed_data, file = TMHMM_aminoacids_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
        } else {
            writeLines("Recovering TMHMM Data...")
            
            TMHMM_protein_parsed_data <- fread(file = TMHMM_proteins_temp_file, header = T, sep = "\t", na.strings = NULL)
            TMHMM_aminoacid_parsed_data <- fread(file = TMHMM_aminoacids_temp_file, header = T, sep = "\t", na.strings = NULL)
        }
    }
    
    ##################################-
    #### **RUN ONCE PER PROTEIN** ####
    ##################################-
    #Declare variable to save data between iterations (its a good idea for it to be a power of 2)
    write_data_every_X_proteins <- 128

    if (use_Iupred | use_NetMHCIIpan | use_NetOglyc | use_Xstream | use_NetSurfp) {
        #Create the temp saving files
        Iupred_temp_file <- sprintf("%s/temp_output_Iupred.tsv", temp_data_folder)
        NetMHCIIpan_alleles_temp_file <- sprintf("%s/temp_output_NetMHCIIpan_alleles.tsv", temp_data_folder)
        NetMHCIIpan_peptides_temp_file <- sprintf("%s/temp_output_NetMHCIIpan_peptides.tsv", temp_data_folder)
        NetOglyc_temp_file <- sprintf("%s/temp_output_NetOglyc.tsv", temp_data_folder)
        Xstream_temp_file <- sprintf("%s/temp_output_Xstream.tsv", temp_data_folder)
        NetSurfp_temp_file <- sprintf("%s/temp_output_NetSurfp.tsv", temp_data_folder)
        
        #Create the base variables
        full_Iupred_parsed_data <- data.table()
        full_NetMHCIIpan_parsed_allele_data <- data.table()
        full_NetMHCIIpan_parsed_peptide_data <- data.table()
        full_NetOglyc_parsed_data <- data.table()
        full_Xstream_parsed_data <- data.table()
        full_NetSurfp_parsed_data <- data.table()
        
        already_parsed_proteins_Iupred <- c()
        already_parsed_proteins_NetMHCIIpan <- c()
        already_parsed_proteins_NetOglyc <- c()
        already_parsed_proteins_Xstream <- c()
        already_parsed_proteins_NetSurfp <- c()
        
        #If they exist, load temp data
        if (use_Iupred & file.exists(Iupred_temp_file)) {
            writeLines("Recovering Iupred Data...")
            full_Iupred_parsed_data <- fread(file = Iupred_temp_file, header = T, sep = "\t", na.strings = NULL)
            already_parsed_proteins_Iupred <- unique(full_Iupred_parsed_data$id)
        }
        if (use_NetMHCIIpan & file.exists(NetMHCIIpan_alleles_temp_file)) {
            writeLines("Recovering NetMHCIIpan Data...")
            full_NetMHCIIpan_parsed_allele_data <- fread(file = NetMHCIIpan_alleles_temp_file, header = T, sep = "\t", na.strings = NULL)
            full_NetMHCIIpan_parsed_peptide_data <- fread(file = NetMHCIIpan_peptides_temp_file, header = T, sep = "\t", na.strings = NULL)
            already_parsed_proteins_NetMHCIIpan <- unique(full_NetMHCIIpan_parsed_peptide_data$id)
        }
        if (use_NetOglyc & file.exists(NetOglyc_temp_file)) {
            writeLines("Recovering NetOglyc Data...")
            full_NetOglyc_parsed_data <- fread(file = NetOglyc_temp_file, header = T, sep = "\t", na.strings = NULL)
            already_parsed_proteins_NetOglyc <- unique(full_NetOglyc_parsed_data$id)
        }
        if (use_Xstream & file.exists(Xstream_temp_file)) {
            writeLines("Recovering Xstream Data...")
            full_Xstream_parsed_data <- fread(file = Xstream_temp_file, header = T, sep = "\t", na.strings = NULL)
            already_parsed_proteins_Xstream <- unique(full_Xstream_parsed_data$id)
        }
        if (use_NetSurfp & file.exists(NetSurfp_temp_file)) {
            writeLines("Recovering NetSurfp Data...")
            full_NetSurfp_parsed_data <- fread(file = NetSurfp_temp_file, header = T, sep = "\t", na.strings = NULL)
            already_parsed_proteins_NetSurfp <- unique(full_NetSurfp_parsed_data$id)
        }
        
        #Calculate the proteins that were already parsed everywhere
        already_parsed_proteins_all <- already_parsed_proteins_Iupred
        already_parsed_proteins_all <- intersect(already_parsed_proteins_all, already_parsed_proteins_NetMHCIIpan)
        already_parsed_proteins_all <- intersect(already_parsed_proteins_all, already_parsed_proteins_NetOglyc)
        already_parsed_proteins_all <- intersect(already_parsed_proteins_all, already_parsed_proteins_Xstream)
        already_parsed_proteins_all <- intersect(already_parsed_proteins_all, already_parsed_proteins_NetSurfp)
        
        #Filter the proteins that were already parsed for all out of the tabbed_fasta
        if (length(already_parsed_proteins_all) > 0) {
            filtered_tabbed_fasta <- tabbed_fasta[!(protein %in% already_parsed_proteins_all)]
        } else {
            filtered_tabbed_fasta <- tabbed_fasta
        }
        
        if (filtered_tabbed_fasta[, .N] > 0) {
            #Loop for each group of proteins to be saved
            loop_amount <- ceiling(filtered_tabbed_fasta[, .N] / write_data_every_X_proteins)
            for (i in 1:loop_amount) {
                # i <- 2
                min_prot_index <- 1 + write_data_every_X_proteins * (i - 1)
                max_prot_index <- min_prot_index + write_data_every_X_proteins - 1
                if (max_prot_index > filtered_tabbed_fasta[, .N]) {
                    max_prot_index <- filtered_tabbed_fasta[, .N]
                }
                
                #Extract the sub_fasta
                sub_tabbed_fasta <- filtered_tabbed_fasta[min_prot_index:max_prot_index]
                
                #Create the mini fastas
                for (protein_index_for in 1:sub_tabbed_fasta[, .N]) {
                    # prot_index_for <- 2
                    fasta_for <- sub_tabbed_fasta[protein_index_for]
                    
                    mini_fasta <- sprintf(">%s\n%s\n", fasta_for$protein, fasta_for$sequence)
                    mini_fasta_file <- sprintf("%s/temp_mini_%s.fasta", temp_data_folder, protein_index_for)
                    write(mini_fasta, file = mini_fasta_file)
                }
                
                ################-
                #### Iupred ####
                ################-
                if (use_Iupred) {
                    writeLines(sprintf("Running Iupred (%s/%s)...", i, loop_amount))
                    
                    #Initialize the parallel processing
                    numCores <- detectCores() #set the cores
                    registerDoParallel(numCores)  #use multicore, set to the number of our cores
                    
                    loop_Iupred_parsed_data <- foreach (protein_index_for = 1:sub_tabbed_fasta[, .N], .combine = rbind) %dopar% {
                        dt_protein_for <- sub_tabbed_fasta[protein_index_for]
                        
                        if (!(dt_protein_for$protein %in% already_parsed_proteins_Iupred)) {
                            mini_fasta_file_for <- sprintf("%s/temp_mini_%s.fasta", temp_data_folder, protein_index_for)
                            
                            Iupred_parsed_data <- runAndParse_Iupred(input_singleProteinFasta_file = mini_fasta_file_for,
                                                                     protein_par = dt_protein_for$protein)
                            
                            Iupred_parsed_data
                        }
                    }
                    
                    #End the parallel processing
                    stopImplicitCluster()
                    
                    #Save or append the data
                    if (full_Iupred_parsed_data[, .N] == 0) {
                        full_Iupred_parsed_data <- loop_Iupred_parsed_data
                        write.table(loop_Iupred_parsed_data, file = Iupred_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
                    } else {
                        full_Iupred_parsed_data <- rbindlist(list(full_Iupred_parsed_data, loop_Iupred_parsed_data))
                        write.table(loop_Iupred_parsed_data, file = Iupred_temp_file, col.names = F, row.names = F, sep = "\t", quote = T, append = T)
                    }   
                }
                
                #####################-
                #### NetMHCIIpan ####
                #####################-
                if (use_NetMHCIIpan) {
                    writeLines(sprintf("Running NetMHCIIpan (%s/%s)...", i, loop_amount))
                    
                    #Initialize the parallel processing
                    numCores <- detectCores() #set the cores
                    registerDoParallel(numCores)  #use multicore, set to the number of our cores
                    
                    loop_NetMHCIIpan_list_output <- foreach (protein_index_for = 1:sub_tabbed_fasta[, .N], .combine = c) %dopar% {
                        # protein_index_for <- 1
                        dt_protein_for <- sub_tabbed_fasta[protein_index_for]
                        
                        if (!(dt_protein_for$protein %in% already_parsed_proteins_NetMHCIIpan)) {
                            mini_fasta_file_for <- sprintf("%s/temp_mini_%s.fasta", temp_data_folder, protein_index_for)
                            
                            list_output <- runAndParse_NetMHCIIpan(input_singleProteinFasta_file = mini_fasta_file_for,
                                                                   NetMHCIIpan_alleles = NetMHCIIpan_alleles,
                                                                   NetMHCIIpan_binding_peptide_length = NetMHCIIpan_binding_peptide_length,
                                                                   list_output_suffix = sprintf("_%s", protein_index_for))
                            
                            list_output
                        }
                    }
                    
                    #End the parallel processing
                    stopImplicitCluster()
                    
                    #Combine the data in data.tables
                    for (protein_index_for in 1:sub_tabbed_fasta[, .N]) {
                        # protein_index_for <- 1
                        NetMHCIIpan_parsed_allele_data_aux <- loop_NetMHCIIpan_list_output[[sprintf("NetMHCIIpan_parsed_allele_data_%s", protein_index_for)]]
                        NetMHCIIpan_parsed_peptide_data_aux <- loop_NetMHCIIpan_list_output[[sprintf("NetMHCIIpan_parsed_peptide_data_%s", protein_index_for)]]
                        
                        if (protein_index_for == 1) {
                            loop_NetMHCIIpan_parsed_allele_data <- NetMHCIIpan_parsed_allele_data_aux
                            loop_NetMHCIIpan_parsed_peptide_data <- NetMHCIIpan_parsed_peptide_data_aux
                        } else {
                            loop_NetMHCIIpan_parsed_allele_data <- rbindlist(list(loop_NetMHCIIpan_parsed_allele_data, NetMHCIIpan_parsed_allele_data_aux))
                            loop_NetMHCIIpan_parsed_peptide_data <- rbindlist(list(loop_NetMHCIIpan_parsed_peptide_data, NetMHCIIpan_parsed_peptide_data_aux))
                        }
                    }
                    rm(loop_NetMHCIIpan_list_output)
                    gc()
                    
                    #Save or append the data
                    if (full_NetMHCIIpan_parsed_allele_data[, .N] == 0) {
                        full_NetMHCIIpan_parsed_allele_data <- loop_NetMHCIIpan_parsed_allele_data
                        full_NetMHCIIpan_parsed_peptide_data <- loop_NetMHCIIpan_parsed_peptide_data
                        write.table(loop_NetMHCIIpan_parsed_allele_data, file = NetMHCIIpan_alleles_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
                        write.table(loop_NetMHCIIpan_parsed_peptide_data, file = NetMHCIIpan_peptides_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
                    } else {
                        full_NetMHCIIpan_parsed_allele_data <- rbindlist(list(full_NetMHCIIpan_parsed_allele_data, loop_NetMHCIIpan_parsed_allele_data))
                        full_NetMHCIIpan_parsed_peptide_data <- rbindlist(list(full_NetMHCIIpan_parsed_peptide_data, loop_NetMHCIIpan_parsed_peptide_data))
                        write.table(loop_NetMHCIIpan_parsed_allele_data, file = NetMHCIIpan_alleles_temp_file, col.names = F, row.names = F, sep = "\t", quote = T, append = T)
                        write.table(loop_NetMHCIIpan_parsed_peptide_data, file = NetMHCIIpan_peptides_temp_file, col.names = F, row.names = F, sep = "\t", quote = T, append = T)
                    }   
                }
                
                ##################-
                #### NetOglyc ####
                ##################-
                if (use_NetOglyc) {
                    writeLines(sprintf("Running NetOglyc (%s/%s)...", i, loop_amount))
                    
                    #Initialize the parallel processing
                    numCores <- detectCores() #set the cores
                    registerDoParallel(numCores)  #use multicore, set to the number of our cores
                    
                    loop_NetOglyc_parsed_data <- foreach (protein_index_for = 1:sub_tabbed_fasta[, .N], .combine = rbind) %dopar% {
                        dt_protein_for <- sub_tabbed_fasta[protein_index_for]
                        
                        if (!(dt_protein_for$protein %in% already_parsed_proteins_NetOglyc)) {
                            mini_fasta_file_for <- sprintf("%s/temp_mini_%s.fasta", temp_data_folder, protein_index_for)
                            
                            NetOglyc_parsed_data <- runAndParse_NetOglyc(input_singleProteinFasta_file = mini_fasta_file_for,
                                                                         protein_par = dt_protein_for$protein,
                                                                         aminoacid_fasta = aminoacid_fasta)
                            
                            NetOglyc_parsed_data
                        }
                    }
                    
                    #End the parallel processing
                    stopImplicitCluster()
                    
                    #Save or append the data
                    if (full_NetOglyc_parsed_data[, .N] == 0) {
                        full_NetOglyc_parsed_data <- loop_NetOglyc_parsed_data
                        write.table(loop_NetOglyc_parsed_data, file = NetOglyc_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
                    } else {
                        full_NetOglyc_parsed_data <- rbindlist(list(full_NetOglyc_parsed_data, loop_NetOglyc_parsed_data))
                        write.table(loop_NetOglyc_parsed_data, file = NetOglyc_temp_file, col.names = F, row.names = F, sep = "\t", quote = T, append = T)
                    }   
                }
                
                #################-
                #### Xstream ####
                #################-
                if (use_Xstream) {
                    writeLines(sprintf("Running Xstream (%s/%s)...", i, loop_amount))
                    
                    #Initialize the parallel processing
                    numCores <- detectCores() #set the cores
                    registerDoParallel(numCores)  #use multicore, set to the number of our cores
                    
                    loop_Xstream_parsed_data <- foreach (protein_index_for = 1:sub_tabbed_fasta[, .N], .combine = rbind) %dopar% {
                        # protein_index_for <- 2
                        dt_protein_for <- sub_tabbed_fasta[protein_index_for]
                        
                        if (!(dt_protein_for$protein %in% already_parsed_proteins_Xstream)) {
                            mini_fasta_file_for <- sprintf("%s/temp_mini_%s.fasta", temp_data_folder, protein_index_for)
                            
                            #Create Xstream folder
                            Xstream_temp_folder <- sprintf("%s/Xstream_temp_%s", temp_data_folder, protein_index_for)
                            if (!dir.exists(Xstream_temp_folder)) {
                                dir.create(Xstream_temp_folder)
                            }
                            
                            Xstream_parsed_data <- runAndParse_Xstream(input_singleProteinFasta_file = mini_fasta_file_for,
                                                                       protein = dt_protein_for$protein,
                                                                       Xstream_temp_folder = Xstream_temp_folder,
                                                                       Xstream_path = Xstream_path)
                            
                            #Delete Xstream folder
                            unlink(Xstream_temp_folder, recursive = T)
                            
                            Xstream_parsed_data
                        }
                    }
                    
                    #End the parallel processing
                    stopImplicitCluster()
                    
                    #Save or append the data
                    if (full_Xstream_parsed_data[, .N] == 0) {
                        full_Xstream_parsed_data <- loop_Xstream_parsed_data
                        write.table(loop_Xstream_parsed_data, file = Xstream_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
                    } else {
                        full_Xstream_parsed_data <- rbindlist(list(full_Xstream_parsed_data, loop_Xstream_parsed_data))
                        write.table(loop_Xstream_parsed_data, file = Xstream_temp_file, col.names = F, row.names = F, sep = "\t", quote = T, append = T)
                    }   
                }
                
                ##################-
                #### NetSurfp ####
                ##################-
                if (use_NetSurfp) {
                    writeLines(sprintf("Running NetSurfp (%s/%s)...", i, loop_amount))
                    
                    #Initialize the parallel processing
                    numCores <- detectCores() #set the cores
                    registerDoParallel(numCores)  #use multicore, set to the number of our cores
                    
                    loop_NetSurfp_parsed_data <- foreach (protein_index_for = 1:sub_tabbed_fasta[, .N], .combine = rbind) %dopar% {
                        dt_protein_for <- sub_tabbed_fasta[protein_index_for]
                        
                        if (!(dt_protein_for$protein %in% already_parsed_proteins_NetSurfp)) {
                            mini_fasta_file_for <- sprintf("%s/temp_mini_%s.fasta", temp_data_folder, protein_index_for)
                            
                            NetSurfp_parsed_data <- runAndParse_NetSurfp(input_singleProteinFasta_file = mini_fasta_file_for)
                            
                            NetSurfp_parsed_data
                        }
                    }
                    
                    #End the parallel processing
                    stopImplicitCluster()
                    
                    #Save or append the data
                    if (full_NetSurfp_parsed_data[, .N] == 0) {
                        full_NetSurfp_parsed_data <- loop_NetSurfp_parsed_data
                        write.table(loop_NetSurfp_parsed_data, file = NetSurfp_temp_file, col.names = T, row.names = F, sep = "\t", quote = T)
                    } else {
                        full_NetSurfp_parsed_data <- rbindlist(list(full_NetSurfp_parsed_data, loop_NetSurfp_parsed_data))
                        write.table(loop_NetSurfp_parsed_data, file = NetSurfp_temp_file, col.names = F, row.names = F, sep = "\t", quote = T, append = T)
                    }   
                }
            }
        }
    }
    
    ########################-
    #### **MERGE DATA** ####
    ########################-
    round_decimals <- 4
    
    ################-
    #### Per aa ####
    ################-
    output_per_aminoacid <- aminoacid_fasta[, .(id = protein, pos, aa)]
    
    # ### DEBUG
    # writeLines(sprintf("%s tiene %s IDs", "BepiPred_parsed_data", length(unique(BepiPred_parsed_data$id))))
    # writeLines(sprintf("%s tiene %s IDs", "Paircoil_parsed_data", length(unique(Paircoil_parsed_data$id))))
    # writeLines(sprintf("%s tiene %s IDs", "SignalP_aminoacid_parsed_data", length(unique(SignalP_aminoacid_parsed_data$id))))
    # writeLines(sprintf("%s tiene %s IDs", "TMHMM_aminoacid_parsed_data", length(unique(TMHMM_aminoacid_parsed_data$id))))
    # writeLines(sprintf("%s tiene %s IDs", "full_Iupred_parsed_data", length(unique(full_Iupred_parsed_data$id))))
    # writeLines(sprintf("%s tiene %s IDs", "full_NetOglyc_parsed_data", length(unique(full_NetOglyc_parsed_data$id))))
    # ###
    
    #BepiPred
    if (use_BepiPred) {
        output_per_aminoacid <- merge(output_per_aminoacid,
                                      BepiPred_parsed_data[, .(id, pos, BepiPred = round(score, round_decimals))],
                                      by = c("id", "pos"))    
    }
    
    #Paircoil2
    if (use_Paircoil2) {
        output_per_aminoacid <- merge(output_per_aminoacid,
                                      Paircoil_parsed_data[, .(id, pos, Paircoil2 = round(p_score, round_decimals))],
                                      by = c("id", "pos"))    
    }
    
    #SignalP
    if (use_SignalP) {
        output_per_aminoacid <- merge(output_per_aminoacid,
                                      SignalP_aminoacid_parsed_data[, .(id, pos,
                                                                        SignalP_C = round(C_score, round_decimals),
                                                                        SignalP_S = round(S_score, round_decimals))],
                                      by = c("id", "pos"))
    }
    
    #TMHMM
    if (use_TMHMM) {
        output_per_aminoacid <- merge(output_per_aminoacid,
                                      TMHMM_aminoacid_parsed_data[, .(id, pos, TMHMM = in_transmembrane_helix)],
                                      by = c("id", "pos"))    
    }
    
    #NetSurfp
    if (use_NetSurfp) {
        output_per_aminoacid <- merge(output_per_aminoacid,
                                      full_NetSurfp_parsed_data[, .(id, pos, 
                                                               NetSurfp_RSA = round(relative_surface_accessibility, round_decimals),
                                                               NetSurfp_AH = round(probability_for_alpha_helix, round_decimals),
                                                               NetSurfp_BS = round(probability_for_beta_strand, round_decimals))],
                                      by = c("id", "pos"))
    }
    
    #Iupred
    if (use_Iupred) {
        output_per_aminoacid <- merge(output_per_aminoacid,
                                      full_Iupred_parsed_data[, .(id, pos, Iupred = round(disorder, round_decimals))],
                                      by = c("id", "pos"))    
    }
    
    #NetOglyc
    if (use_NetOglyc) {
        output_per_aminoacid <- merge(output_per_aminoacid,
                                      full_NetOglyc_parsed_data[, .(id, pos, NetOglyc = is_glycosylation_site)],
                                      by = c("id", "pos"))    
    }

    #####################-
    #### Per peptide ####
    #####################-
    output_per_peptide <- full_NetMHCIIpan_parsed_peptide_data
    
    ####################-
    #### Per repeat ####
    ####################-
    full_Xstream_parsed_data <- full_Xstream_parsed_data[(start != 0) | (end != 0)] #remove non repeats
    output_per_repeat <- full_Xstream_parsed_data
    
    #####################-
    #### Per protein ####
    #####################-
    output_per_protein <- tabbed_fasta[, .(id = protein, original_id = original_protein)]

    #PredGPI
    if (use_PredGPI) {
        output_per_protein <- merge(output_per_protein,
                                    PredGPI_parsed_data[, .(id,
                                                            PredGPI_start = gpi_start,
                                                            PredGPI_score = round(score, round_decimals))],
                                    by = c("id"))    
    }
    
    #SignalP
    if (use_SignalP) {
        output_per_protein <- merge(output_per_protein,
                                    SignalP_protein_parsed_data[, .(id, SignalP_end = cleavage_site)],
                                    by = c("id"))
    }
    
    #TMHMM
    if (use_TMHMM) {
        output_per_protein <- merge(output_per_protein,
                                    TMHMM_protein_parsed_data[, .(id, TMHMM = has_transmembrane_helix)],
                                    by = c("id"))
    }
    
    #Isoelectric point & Molecular weight
    if (use_IsoelectricPoint | use_MolecularWeight) {
        output_per_protein <- merge(output_per_protein,
                                    IP_MW_parsed_data[, .(id,
                                                               MolecularWeight = round(molecular_weight, round_decimals),
                                                               IsoelectricPoint = round(isoelectric_point, round_decimals))],
                                    by = c("id"))
    }
    
    ######################################-
    #### **KMER SIMILARITY ANALYSIS** ####
    ######################################-
    #########################-
    #### Self Similarity ####
    #########################-
    #Pathogen vs itself
    if (use_SelfSimilarity | use_CrossReactivity | use_Coendemicity) {
        if (use_SelfSimilarity) {
            writeLines("Calculating Self Similarity...")
        } else {
            writeLines("Calculating kmers...")
        }
        
        SelfSimilarity_output_file <- sprintf("%s/SelfSimilarity_output_file", temp_data_folder)
        command_aux <- sprintf("perl '%s/lib/split_proteins_into_kmers.pl' '%s' %s %s '%s' '%s'",
                               scripts_folder,
                               temp_fasta_file,
                               KmerSimilarity_kmer_length,
                               max_protein_length,
                               replace_nonAA_chars_by,
                               SelfSimilarity_output_file)
        system(command_aux)
        
        SelfSimilarity_parsed_data <- read.delim(file = SelfSimilarity_output_file, header = T, sep = "", na.strings = NULL, comment.char = "")
        SelfSimilarity_parsed_data <- as.data.table(SelfSimilarity_parsed_data)
        setnames(SelfSimilarity_parsed_data, "protein", "id")
        
        output_per_kmer <- SelfSimilarity_parsed_data
        
        rm(SelfSimilarity_parsed_data)
        gc()
    } else {
        output_per_kmer <- data.table(id = character(),
                                      kmer = character(),
                                      quantity = integer(),
                                      quantity_in_proteome = integer())
    }
    
    #########################-
    #### CrossReactivity ####
    #########################-
    #Pathogen vs host
    if (use_CrossReactivity) {
        writeLines("Calculating Cross Reactivity...")
        
        CrossReactivity_output_file <- sprintf("%s/CrossReactivity_output_file", temp_data_folder)
        command_aux <- sprintf("perl '%s/lib/split_proteome_into_kmers.pl' '%s' %s %s '%s' '%s'",
                               scripts_folder,
                               CrossReactivity_fasta_file,
                               KmerSimilarity_kmer_length,
                               max_protein_length,
                               replace_nonAA_chars_by,
                               CrossReactivity_output_file)
        system(command_aux)
        
        CrossReactivity_parsed_data <- read.delim(file = CrossReactivity_output_file, header = T, sep = "", na.strings = NULL, comment.char = "")
        CrossReactivity_parsed_data <- as.data.table(CrossReactivity_parsed_data)
        
        output_per_kmer <- merge(output_per_kmer,
                                 CrossReactivity_parsed_data[, .(kmer, quantity_in_host_proteome = quantity)],
                                 by = "kmer",
                                 all.x = T)
        output_per_kmer[is.na(quantity_in_host_proteome), quantity_in_host_proteome := 0]
        
        rm(CrossReactivity_parsed_data)
        gc()
    } else {
        output_per_kmer$quantity_in_host_proteome <- 0
    }
    
    ######################-
    #### Coendemicity ####
    ######################-
    #Pathogen vs other pathogen
    if (use_Coendemicity) {
        writeLines("Calculating Coendemicity...")
        
        Coendemicity_output_file <- sprintf("%s/Coendemicity_output_file", temp_data_folder)
        command_aux <- sprintf("perl '%s/lib/split_proteome_into_kmers.pl' '%s' %s %s '%s' '%s'",
                               scripts_folder,
                               Coendemicity_fasta_file,
                               KmerSimilarity_kmer_length,
                               max_protein_length,
                               replace_nonAA_chars_by,
                               Coendemicity_output_file)
        system(command_aux)
        
        Coendemicity_parsed_data <- read.delim(file = Coendemicity_output_file, header = T, sep = "", na.strings = NULL, comment.char = "")
        Coendemicity_parsed_data <- as.data.table(Coendemicity_parsed_data)
        
        output_per_kmer <- merge(output_per_kmer,
                                 Coendemicity_parsed_data[, .(kmer, quantity_in_coendemic_proteome = quantity)],
                                 by = "kmer",
                                 all.x = T)
        output_per_kmer[is.na(quantity_in_coendemic_proteome), quantity_in_coendemic_proteome := 0]
        
        rm(Coendemicity_parsed_data)
        gc()
    } else {
        output_per_kmer$quantity_in_coendemic_proteome <- 0
    }
    
    ############################################-
    #### **ADD THE ORIGINAL PROTEIN NAMES** ####
    ############################################-
    
    #########################-
    #### **OUTPUT DATA** ####
    #########################-
    writeLines("Readying Output...")
    
    #Fix some data types
    output_per_aminoacid$aa <- as.character(output_per_aminoacid$aa)
    if (output_per_kmer[, .N] > 0) {
        output_per_kmer$kmer <- as.character(output_per_kmer$kmer)    
    }
    if (output_per_peptide[, .N] > 0) {
        output_per_peptide$peptide <- as.character(output_per_peptide$peptide)
    }
    
    if (output_method %in% c("write", "both")) {
        write.table(output_per_aminoacid, file = output_per_aminoacid_file, col.names = T, row.names = F, sep = "\t", quote = T)
        write.table(output_per_protein, file = output_per_protein_file, col.names = T, row.names = F, sep = "\t", quote = T)
        if (output_per_kmer_file != "") {
            write.table(output_per_kmer, file = output_per_kmer_file, col.names = T, row.names = F, sep = "\t", quote = T)
        }
        if (output_per_peptide_file != "") {
            write.table(output_per_peptide, file = output_per_peptide_file, col.names = T, row.names = F, sep = "\t", quote = T)
        }
        if (output_per_repeat_file != "") {
            write.table(output_per_repeat, file = output_per_repeat_file, col.names = T, row.names = F, sep = "\t", quote = T)
        }
    }
    
    list_output <- list()
    if (output_method %in% c("list", "both")) {
        #Create the output
        list_output[["output_per_aminoacid"]] <- output_per_aminoacid
        list_output[["output_per_protein"]] <- output_per_protein
        list_output[["output_per_kmer"]] <- output_per_kmer
        list_output[["output_per_peptide"]] <- output_per_peptide
        list_output[["output_per_repeat"]] <- output_per_repeat
    }
    
    writeLines("Done!")
    
    list_output
}
