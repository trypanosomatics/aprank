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
            } else if (grepl("Charge", line)) {
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

runAndParse_NetSurfp_Parallel <- function(input_fasta_files,
                                 number_of_parallel_processes) {
    #NetSurfp already has parallelization, so it's convenient to use less parallel processes, such as
    #leaving 4 cores free per process
    numCores <- detectCores() #set the cores
    
    registerDoParallel(numCores)  #use multicore, set to the number of our cores
    
    full_NetSurfp_parsed_data <- foreach (i = 1:number_of_parallel_processes, .combine = rbind) %dopar% {
        #This can give core dumped errors, but I think NetSurfp just do those again because the outputs are full
        command_aux <- sprintf("netsurfp -i '%s' -a",
                               input_fasta_files[i])
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
    #End the parallel processing
    stopImplicitCluster()
    
    full_NetSurfp_parsed_data
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

runAndParse_NetMHCIIPan <- function(input_singleProteinFasta_file,
                                    NetMHCIIPan_alleles,
                                    NetMHCIIpan_binding_peptide_length) {
    #Run the program for each allele and merge the data in one file
    for (allele_for in NetMHCIIPan_alleles) {
        # allele_for <- NetMHCIIPan_alleles[1]
        command_aux <- sprintf("netMHCIIpan -a '%s' -f '%s' -l %s",
                               allele_for,
                               input_singleProteinFasta_file,
                               NetMHCIIpan_binding_peptide_length)
        NetMHCIIpan_output_aux <- system(command_aux, intern = T)
        
        if (allele_for == NetMHCIIPan_alleles[1]) {
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
    list_output[["NetMHCIIpan_parsed_allele_data"]] <- NetMHCIIpan_parsed_allele_data
    list_output[["NetMHCIIpan_parsed_peptide_data"]] <- NetMHCIIpan_parsed_peptide_data
    
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
                                temp_data_folder,
                                Xstream_path) {
    #I can't change the name of the output, but I can change it's location
    Xstream_output_file <- sprintf("%s/\\XSTREAM__i0.7_g3_m5_e2.0_out_2.html", temp_data_folder)
    
    #I don't care about the log, but otherwise it's shown on console
    command_aux <- sprintf("java -jar '%s' '%s' -d'%s/' >'%s/xstream.log'",
                           Xstream_path,
                           input_singleProteinFasta_file,
                           temp_data_folder,
                           temp_data_folder)
    system(command_aux)
    
    #Extract important data
    #There is a warning about the final line not ending properly, but that's Xstream's problem and it doesn't affect this
    con <- file(Xstream_output_file, "r")
    Xstream_data <- c()
    while(TRUE) {
        line = readLines(con, 1)
        
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
        Xstream_parsed_data <- data.table(id = character(),
                                          start = integer(),
                                          end = integer(),
                                          period = numeric(),
                                          copy_number = numeric(),
                                          consensus_error = numeric())
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
                          NetMHCIIpan_binding_peptide_length = 9, NetMHCIIPan_alleles = c("DRB1_0101", "DRB3_0101", "DRB4_0101", "DRB5_0101"),
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
    tabbed_fasta[, sequence_length := nchar(sequence)]
    aminoacid_fasta <- fread(temp_aminoacid_fasta_file, header = T, sep = "\t", na.strings = NULL)
    
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
            # i <- 4
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
    
    #Create the subFASTAs for parallel processing NetSurfp
    sub_temp_fastas_files_for_NetSurfp <- c()
    tabbed_fasta$parallel_process_NetSurfp_i <- 1
    if ((use_NetSurfp) & (number_of_parallel_processes_for_NetSurfp > 1)) {
        # tabbed_fasta[, output_NetSurfp_aux := sprintf(">%s\n%s", protein, sequence)]
        
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
            
            # sub_fasta_output <- paste(tabbed_fasta[protein_subset_min_for_parallel:protein_subset_max_for_parallel]$output_NetSurfp_aux, collapse = "\n")
            # sub_temp_fasta_file <- sprintf("%s/sub_temp_NetSurfp_%s.fasta", temp_data_folder, i)
            # 
            # write(sub_fasta_output, file = sub_temp_fasta_file)
            # 
            # sub_temp_fastas_files_for_NetSurfp <- c(sub_temp_fastas_files_for_NetSurfp, sub_temp_fasta_file)
        }
        
        # tabbed_fasta <- tabbed_fasta[, -c("output_NetSurfp_aux")]
    }
    
    #Create the folders for each Xstream output (can't change the name)
    sub_Xstream_temp_folders <- c()
    if (use_Xstream) {
        for (i in 1:number_of_parallel_processes) {
            sub_Xstream_temp_folder <- sprintf("%s/Xstream_temp_%s", temp_data_folder, i)
            if (!dir.exists(sub_Xstream_temp_folder)) {
                dir.create(sub_Xstream_temp_folder, recursive = T)
            }
            
            sub_Xstream_temp_folders <- c(sub_Xstream_temp_folders, sub_Xstream_temp_folder)
        }
    }
    
    ###################################-
    #### **RUN ONCE PER PROTEOME** ####
    ###################################-
    ##################-
    #### BepiPred ####
    ##################-
    if (use_BepiPred) {
        writeLines("Running BepiPred...")
        BepiPred_parsed_data <- runAndParse_BepiPred(input_fasta_files = sub_temp_fastas_files,
                                                     number_of_parallel_processes = number_of_parallel_processes)    
    }
    
    ##############################################-
    #### Isoelectric point & Molecular weight ####
    ##############################################-
    if (use_IsoelectricPoint | use_MolecularWeight) {
        writeLines("Running EMBOSS...")
        IP_MW_parsed_data <- runAndParse_IP_MW(input_fasta_files = sub_temp_fastas_files,
                                               temp_data_folder = temp_data_folder,
                                               number_of_parallel_processes = number_of_parallel_processes)
    }
    
    ###################-
    #### Paircoil2 ####
    ###################-
    if (use_Paircoil2) {
        writeLines("Running Paircoil2...")
        Paircoil_parsed_data <- runAndParse_Paircoil2(input_fasta_files = sub_temp_fastas_files,
                                                      temp_data_folder = temp_data_folder,
                                                      aminoacid_fasta = aminoacid_fasta,
                                                      tabbed_fasta = tabbed_fasta,
                                                      number_of_parallel_processes = number_of_parallel_processes)    
    }
    
    #################-
    #### PredGPI ####
    #################-
    if (use_PredGPI) {
        writeLines("Running PredGPI...")
        PredGPI_parsed_data <- runAndParse_PredGPI(temp_data_folder = temp_data_folder,
                                                   tabbed_fasta = tabbed_fasta,
                                                   PredGPI_min_length = 41,
                                                   number_of_parallel_processes = number_of_parallel_processes)    
    }
    
    #################-
    #### SignalP ####
    #################-
    if (use_SignalP) {
        writeLines("Running SignalP...")
        list_output <- runAndParse_SignalP(input_fasta_files = sub_temp_fastas_files,
                                           SignalP_organism_group = SignalP_organism_group,
                                           aminoacid_fasta = aminoacid_fasta,
                                           number_of_parallel_processes = number_of_parallel_processes)
        SignalP_protein_parsed_data <- list_output[["SignalP_protein_parsed_data"]]
        SignalP_aminoacid_parsed_data <- list_output[["SignalP_aminoacid_parsed_data"]]
        rm(list_output)
        gc()    
    }
    
    ###############-
    #### TMHMM ####
    ###############-
    if (use_TMHMM) {
        writeLines("Running TMHMM...")
        list_output <- runAndParse_TMHMM(input_fasta_files = sub_temp_fastas_files,
                                         aminoacid_fasta = aminoacid_fasta,
                                         number_of_parallel_processes = number_of_parallel_processes)
        TMHMM_protein_parsed_data <- list_output[["TMHMM_protein_parsed_data"]]
        TMHMM_aminoacid_parsed_data <- list_output[["TMHMM_aminoacid_parsed_data"]]
        rm(list_output)
        gc()    
    }
    
    ##################################-
    #### **RUN ONCE PER PROTEIN** ####
    ##################################-
    numCores <- detectCores() #set the cores
    
    registerDoParallel(numCores)  #use multicore, set to the number of our cores
    
    full_list_output <- foreach (i = 1:number_of_parallel_processes, .combine = c) %dopar% {
        full_Iupred_parsed_data <- data.table()
        full_NetMHCIIPan_parsed_allele_data <- data.table()
        full_NetMHCIIPan_parsed_peptide_data <- data.table()
        full_NetOglyc_parsed_data <- data.table()
        full_Xstream_parsed_data <- data.table()
        
        sub_tabbed_fasta <- tabbed_fasta[parallel_process_i == i]
        sub_protein_amount <- sub_tabbed_fasta[, .N]
        if (use_Iupred | use_NetMHCIIpan | use_NetOglyc | use_Xstream) {
            for (j in 1:sub_protein_amount) {
                #j <- 1
                #fasta_for <- sub_tabbed_fasta[original_protein == "XP_001348321.1"]
                fasta_for <- sub_tabbed_fasta[j]
                
                writeLines(sprintf("Processing data for %s (%s/%s of core %s/%s)...", fasta_for$original_protein, j, sub_protein_amount, i, number_of_parallel_processes))
                
                mini_fasta <- sprintf(">%s\n%s\n", fasta_for$protein, fasta_for$sequence)
                mini_fasta_file <- sprintf("%s/temp_mini_%s.fasta", temp_data_folder, i)
                write(mini_fasta, file = mini_fasta_file)
                
                ################-
                #### Iupred ####
                ################-
                if (use_Iupred) {
                    writeLines(sprintf("Running Iupred for %s...", fasta_for$original_protein))
                    
                    Iupred_parsed_data <- runAndParse_Iupred(input_singleProteinFasta_file = mini_fasta_file,
                                                             protein_par = fasta_for$protein)
                    
                    if (full_Iupred_parsed_data[, .N] == 0) {
                        full_Iupred_parsed_data <- Iupred_parsed_data
                    } else {
                        full_Iupred_parsed_data <- rbindlist(list(full_Iupred_parsed_data, Iupred_parsed_data))
                    }    
                }
                
                #####################-
                #### NetMHCIIpan ####
                #####################-
                if (use_NetMHCIIpan) {
                    writeLines(sprintf("Running NetMHCIIpan for %s...", fasta_for$original_protein))
                    
                    list_output <- runAndParse_NetMHCIIPan(input_singleProteinFasta_file = mini_fasta_file,
                                                           NetMHCIIPan_alleles = NetMHCIIPan_alleles,
                                                           NetMHCIIpan_binding_peptide_length = NetMHCIIpan_binding_peptide_length)
                    NetMHCIIpan_parsed_allele_data <- list_output[["NetMHCIIpan_parsed_allele_data"]]
                    NetMHCIIpan_parsed_peptide_data <- list_output[["NetMHCIIpan_parsed_peptide_data"]]
                    rm(list_output)
                    gc()
                    
                    if (full_NetMHCIIPan_parsed_allele_data[, .N] == 0) {
                        full_NetMHCIIPan_parsed_allele_data <- NetMHCIIpan_parsed_allele_data
                    } else {
                        full_NetMHCIIPan_parsed_allele_data <- rbindlist(list(full_NetMHCIIPan_parsed_allele_data, NetMHCIIpan_parsed_allele_data))
                    }
                    if (full_NetMHCIIPan_parsed_peptide_data[, .N] == 0) {
                        full_NetMHCIIPan_parsed_peptide_data <- NetMHCIIpan_parsed_peptide_data
                    } else {
                        full_NetMHCIIPan_parsed_peptide_data <- rbindlist(list(full_NetMHCIIPan_parsed_peptide_data, NetMHCIIpan_parsed_peptide_data))
                    }
                }
                
                ##################-
                #### NetOglyc ####
                ##################-
                if (use_NetOglyc) {
                    writeLines(sprintf("Running NetOglyc for %s...", fasta_for$original_protein))
                    
                    NetOglyc_parsed_data <- runAndParse_NetOglyc(input_singleProteinFasta_file = mini_fasta_file,
                                                                 protein_par = fasta_for$protein,
                                                                 aminoacid_fasta = aminoacid_fasta)
                    
                    if (full_NetOglyc_parsed_data[, .N] == 0) {
                        full_NetOglyc_parsed_data <- NetOglyc_parsed_data
                    } else {
                        full_NetOglyc_parsed_data <- rbindlist(list(full_NetOglyc_parsed_data, NetOglyc_parsed_data))
                    }    
                }
                
                #################-
                #### XStream ####
                #################-
                if (use_Xstream) {
                    writeLines(sprintf("Running XStream for %s...", fasta_for$original_protein))
                    
                    Xstream_parsed_data <- runAndParse_Xstream(input_singleProteinFasta_file = mini_fasta_file,
                                                               protein = fasta_for$protein,
                                                               temp_data_folder = sub_Xstream_temp_folders[i],
                                                               Xstream_path = Xstream_path)
                    
                    if (full_Xstream_parsed_data[, .N] == 0) {
                        full_Xstream_parsed_data <- Xstream_parsed_data
                    } else {
                        full_Xstream_parsed_data <- rbindlist(list(full_Xstream_parsed_data, Xstream_parsed_data))
                    }
                }
            }
        }
        
        #Create the list output
        list_output <- list()
        list_output[[sprintf("full_Iupred_parsed_data_%s", i)]] <- full_Iupred_parsed_data
        list_output[[sprintf("full_NetMHCIIPan_parsed_allele_data_%s", i)]] <- full_NetMHCIIPan_parsed_allele_data
        list_output[[sprintf("full_NetMHCIIPan_parsed_peptide_data_%s", i)]] <- full_NetMHCIIPan_parsed_peptide_data
        list_output[[sprintf("full_NetOglyc_parsed_data_%s", i)]] <- full_NetOglyc_parsed_data
        list_output[[sprintf("full_Xstream_parsed_data_%s", i)]] <- full_Xstream_parsed_data
        
        list_output
    }
    
    #End the parallel processing
    stopImplicitCluster()
    
    #Paste the outputs together
    for (i in 1:number_of_parallel_processes) {
        #i <- 2
        if (i == 1) {
            full_Iupred_parsed_data <- full_list_output[[sprintf("full_Iupred_parsed_data_%s", i)]]
            full_NetMHCIIPan_parsed_allele_data <- full_list_output[[sprintf("full_NetMHCIIPan_parsed_allele_data_%s", i)]]
            full_NetMHCIIPan_parsed_peptide_data <- full_list_output[[sprintf("full_NetMHCIIPan_parsed_peptide_data_%s", i)]]
            full_NetOglyc_parsed_data <- full_list_output[[sprintf("full_NetOglyc_parsed_data_%s", i)]]
            full_Xstream_parsed_data <- full_list_output[[sprintf("full_Xstream_parsed_data_%s", i)]]
        } else {
            full_Iupred_parsed_data <- rbindlist(list(full_Iupred_parsed_data,
                                                      full_list_output[[sprintf("full_Iupred_parsed_data_%s", i)]]))
            full_NetMHCIIPan_parsed_allele_data <- rbindlist(list(full_NetMHCIIPan_parsed_allele_data,
                                                                  full_list_output[[sprintf("full_NetMHCIIPan_parsed_allele_data_%s", i)]]))
            full_NetMHCIIPan_parsed_peptide_data <- rbindlist(list(full_NetMHCIIPan_parsed_peptide_data,
                                                                   full_list_output[[sprintf("full_NetMHCIIPan_parsed_peptide_data_%s", i)]]))
            full_NetOglyc_parsed_data <- rbindlist(list(full_NetOglyc_parsed_data,
                                                        full_list_output[[sprintf("full_NetOglyc_parsed_data_%s", i)]]))
            full_Xstream_parsed_data <- rbindlist(list(full_Xstream_parsed_data,
                                                       full_list_output[[sprintf("full_Xstream_parsed_data_%s", i)]]))
        }
    }
    rm(full_list_output)
    gc()
    
    ###########################################-
    #### **RUN NETSURFP ONCE PER PROTEIN** ####
    ###########################################-
    ### NetSurfp can theoretically run using large FASTAs, however it gave a lot of errors (core dumps) which
    ### made it have to start over, and if it had many proteins it had to do that with all of them. This way
    ### even if it has to do something all over again it will be a lot shorter
    ### The bad news about this is that I'm creating 1 file per protein twice
    NetSurfp_write_data_every_X_proteins <- 100
    
    if (use_NetSurfp) {
        numCores <- detectCores() #set the cores
        
        registerDoParallel(numCores)  #use multicore, set to the number of our cores
        
        #I tried using combine as I did for the other predictors, but NetSurfp uses too much memory and that
        #resulted in errors, so now I'll write the data into files as it is being processed
        foreach (i = 1:number_of_parallel_processes_for_NetSurfp) %dopar% {
            # i <- 1
            full_NetSurfp_parsed_data <- data.table()
            
            sub_NetSurfp_parsed_data_file <- sprintf("%s/NetSurfp_sub_parsed_data_temp_%s.tsv", temp_data_folder, i)
            
            sub_tabbed_fasta <- tabbed_fasta[parallel_process_NetSurfp_i == i]
            sub_protein_amount <- sub_tabbed_fasta[, .N]
            for (j in 1:sub_protein_amount) {
                # j <- 1
                fasta_for <- sub_tabbed_fasta[j]
                
                writeLines(sprintf("Running NetSurfp for %s (%s/%s of core %s/%s)...", fasta_for$original_protein, j, sub_protein_amount, i, number_of_parallel_processes_for_NetSurfp))
                
                mini_fasta <- sprintf(">%s\n%s\n", fasta_for$protein, fasta_for$sequence)
                mini_fasta_file <- sprintf("%s/temp_mini_%s.fasta", temp_data_folder, i)
                write(mini_fasta, file = mini_fasta_file)
                
                ##################-
                #### NetSurfp ####
                ##################-
                NetSurfp_parsed_data <- runAndParse_NetSurfp(input_singleProteinFasta_file = mini_fasta_file)
                
                if (j == 1) {
                    #The first time just create the file and the empty structure for the data
                    sub_full_NetSurfp_parsed_data <- NetSurfp_parsed_data
                    write.table(sub_full_NetSurfp_parsed_data,file = sub_NetSurfp_parsed_data_file,
                                col.names = T, row.names = F, sep = "\t", quote = T)
                    sub_full_NetSurfp_parsed_data <- NetSurfp_parsed_data[0,]
                } else if ((j == sub_protein_amount) | ((j %% NetSurfp_write_data_every_X_proteins) == 0)) {
                    #If it's the last time, or time to write the data, then write it and empty it
                    sub_full_NetSurfp_parsed_data <- rbindlist(list(sub_full_NetSurfp_parsed_data, NetSurfp_parsed_data))
                    write.table(sub_full_NetSurfp_parsed_data, file = sub_NetSurfp_parsed_data_file,
                                col.names = F, row.names = F, sep = "\t", quote = T, append = T)
                    sub_full_NetSurfp_parsed_data <- NetSurfp_parsed_data[0,]
                } else {
                    #Otherwise, just save the data in the data table to be saved later on
                    sub_full_NetSurfp_parsed_data <- rbindlist(list(sub_full_NetSurfp_parsed_data, NetSurfp_parsed_data))
                }
            }
        }
        
        #End the parallel processing
        stopImplicitCluster()
        
        #Paste the outputs together
        for (i in 1:number_of_parallel_processes_for_NetSurfp) {
            sub_NetSurfp_parsed_data_file <- sprintf("%s/NetSurfp_sub_parsed_data_temp_%s.tsv", temp_data_folder, i)
            sub_NetSurfp_parsed_data <- fread(sub_NetSurfp_parsed_data_file, header = T, sep = "\t", na.strings = NULL)
            
            if (i == 1) {
                full_NetSurfp_parsed_data <- sub_NetSurfp_parsed_data
            } else {
                full_NetSurfp_parsed_data <- rbindlist(list(full_NetSurfp_parsed_data,
                                                            sub_NetSurfp_parsed_data))
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
    output_per_peptide <- full_NetMHCIIPan_parsed_peptide_data
    
    ####################-
    #### Per repeat ####
    ####################-
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

##################-
#### **CALL** ####
##################-
##################################-
#### Call - Select Predictors ####
##################################-
# ## General
# use_BepiPred <- 1
# use_IsoelectricPoint <- 1
# use_MolecularWeight <- 1
# use_Paircoil2 <- 1
# use_PredGPI <- 1
# use_SignalP <- 1
# use_TMHMM <- 1
# use_NetSurfp <- 1
# 
# use_Iupred <- 1
# use_NetOglyc <- 1
# use_Xstream <- 1
# use_NetMHCIIpan <- 1
# 
# use_SelfSimilarity <- 1
# use_CrossReactivity <- 1
# use_Coendemicity <- 0
# 
# ###############################-
# #### Call - Run Predictors ####
# ###############################-
# ## General
# setwd("/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK")
# 
# input_data_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK"
# temp_data_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK/tmp"
# output_data_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK"
# scripts_folder <- "/home/alejandror/Workspace/Dropbox_Link/Programas/20170002_Pepranker/APRANK/source2"
# 
# # input_fasta_file <- sprintf("%s/88_tests_and_stuff/test_2.fasta", input_data_folder)
# input_fasta_file <- sprintf("%s/88_tests_and_stuff/cox/fasta to check/mix.fasta", input_data_folder)
# 
# temp_fasta_file <- sprintf("%s/temp.fasta", temp_data_folder)
# temp_tabbed_fasta_file <- sprintf("%s/temp_fasta_tabbed.tsv", temp_data_folder)
# temp_splitted_fasta_file <- sprintf("%s/temp_fasta_splitted.tsv", temp_data_folder)
# temp_aminoacid_fasta_file <- sprintf("%s/temp_fasta_aminoacids.tsv", temp_data_folder)
# 
# #The output method can be "write" (makes files), "list" (returns a list), or "both"
# output_method <- "write"
# output_per_aminoacid_file <- sprintf("%s/output_per_aminoacid.tsv", output_data_folder)
# output_per_protein_file <- sprintf("%s/output_per_protein.tsv", output_data_folder)
# output_per_kmer_file <- sprintf("%s/output_per_kmer.tsv", output_data_folder)
# output_per_peptide_file <- sprintf("%s/output_per_peptide.tsv", output_data_folder)
# output_per_repeat_file <- sprintf("%s/output_per_repeat.tsv", output_data_folder)
# 
# #This has to be 1 or more (1 meaning just 1 process without parallelization)
# number_of_parallel_processes <- 1
# 
# peptide_length <- 15
# peptide_overlap <- 14
# 
# ## SignalP
# #The organism group of the genome to analyze. It can be: "euk", "gram+" or "gram-".
# SignalP_organism_group <- "gram-"
# 
# ## Xstream
# Xstream_path <- "/usr/local/xstream/xstream.jar"
# 
# ## NetMHCIIPan
# # The length of the fragment used when analyzing with NetMHCIIpan (integer between 9 and 50).
# NetMHCIIpan_binding_peptide_length <- 9
# # The names of the MHC class II alleles to consider when analyzing
# NetMHCIIPan_alleles <- c("DRB1_0101",
#                          "DRB3_0101",
#                          "DRB4_0101",
#                          "DRB5_0101")
# 
# ## Self Similarity, CrossReactivity & Coendemicity
# KmerSimilarity_kmer_length <- 6
# 
# ## CrossReactivity (Pathogen vs Host)
# CrossReactivity_fasta_file <- sprintf("%s/88_tests_and_stuff/test_crossreactivity.fasta", input_data_folder)
# 
# ## Coendemicity (Pathogen vs Other pathogen)
# Coendemicity_fasta_file <- sprintf("%s/88_tests_and_stuff/test_coendemicity.fasta", input_data_folder)
# 
# runPredictors(input_fasta_file = input_fasta_file,
#               temp_data_folder = temp_data_folder, output_data_folder = output_data_folder, scripts_folder = scripts_folder,
#               temp_fasta_file = temp_fasta_file, temp_tabbed_fasta_file = temp_tabbed_fasta_file, temp_splitted_fasta_file = temp_splitted_fasta_file, temp_aminoacid_fasta_file = temp_aminoacid_fasta_file,
#               output_method = output_method,
#               output_per_aminoacid_file = output_per_aminoacid_file, output_per_protein_file = output_per_protein_file, output_per_kmer_file = output_per_kmer_file, output_per_peptide_file = output_per_peptide_file, output_per_repeat_file = output_per_repeat_file,
#               number_of_parallel_processes = number_of_parallel_processes,
#               peptide_length = peptide_length, peptide_overlap = peptide_overlap,
#               SignalP_organism_group = SignalP_organism_group,
#               Xstream_path = Xstream_path,
#               NetMHCIIpan_binding_peptide_length = NetMHCIIpan_binding_peptide_length, NetMHCIIPan_alleles = NetMHCIIPan_alleles,
#               KmerSimilarity_kmer_length = KmerSimilarity_kmer_length, CrossReactivity_fasta_file = CrossReactivity_fasta_file, Coendemicity_fasta_file = Coendemicity_fasta_file,
#               use_BepiPred = use_BepiPred, use_Paircoil2 = use_Paircoil2, use_PredGPI = use_PredGPI, use_SignalP = use_SignalP,
#               use_TMHMM = use_TMHMM, use_NetSurfp = use_NetSurfp,
#               use_Iupred = use_Iupred, use_NetOglyc = use_NetOglyc, use_Xstream = use_Xstream, use_NetMHCIIpan = use_NetMHCIIpan,
#               use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight = use_MolecularWeight,
#               use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
