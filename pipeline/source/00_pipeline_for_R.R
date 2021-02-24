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
runPipeline <- 
    function(aprank_folder,
             use_BepiPred, use_IsoelectricPoint, use_Iupred, use_MolecularWeight,
             use_NetMHCIIpan, use_NetOglyc, use_NetSurfp, use_Paircoil2,
             use_PredGPI, use_SignalP, use_TMHMM, use_Xstream,
             use_SelfSimilarity, use_CrossReactivity, use_Coendemicity,
             output_data_folder, input_fasta_file, CrossReactivity_fasta_file, Coendemicity_fasta_file,
             number_of_parallel_processes,
             peptide_length, peptide_overlap, KmerSimilarity_kmer_length,
             SignalP_organism_group,
             NetMHCIIpan_binding_peptide_length, NetMHCIIPan_alleles,
             Paircoil2_fragment_length, Paircoil2_threshold,
             Xstream_min_period, Xstream_min_copy_number, Xstream_max_consensus_error,
             Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_protein_start_penalty_proportion, Coendemicity_protein_max_penalty_proportion,
             Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_peptide_start_penalty_proportion, Coendemicity_peptide_max_penalty_proportion) {
        #############################################-
        #### **INTERNAL CONFIG (BE CAREFUL WHEN EDITING THESE)** ####
        #############################################-
        ##########################-
        #### Config - General ####
        ##########################-
        setwd(aprank_folder)
        
        #Set the base folders
        temp_data_folder <- sprintf("%s/tmp", aprank_folder)
        scripts_folder <-  sprintf("%s/source", aprank_folder)
        model_data_folder <- sprintf("%s/models", aprank_folder)
        
        if ((output_data_folder == "") | (dir.exists(output_data_folder) == FALSE)) {
            output_data_folder <- sprintf("%s/output", aprank_folder)
        }
        
        #Set the subfolders
        aux_name <- basename(input_fasta_file)
        if (grepl("\\.", aux_name)) {
            aux_name <- gsub("(^.+)\\.[^\\.]*$", "\\1", aux_name)
        }
        aux_name <- gsub("[^a-zA-z0-9 \\.]", "", aux_name)
        aux_name <- sprintf("%s_%s", aux_name, format(Sys.time(), "%Y%m%d-%H%M%S"))
        
        temp_data_folder <- sprintf("%s/%s", temp_data_folder, aux_name)
        dir.create(temp_data_folder)
        output_data_folder <- sprintf("%s/%s", output_data_folder, aux_name)
        dir.create(output_data_folder)
        
        #Set the rest
        max_protein_length <- 9999
        replace_nonAA_chars_by <- "X"
        remove_temp_folder <- 1
        
        #################################-
        #### Config - Run Predictors ####
        #################################-
        ## General
        source(sprintf("%s/01_runPredictors.R", scripts_folder))
        
        temp_fasta_file <- sprintf("%s/temp.fasta", temp_data_folder)
        temp_tabbed_fasta_file <- sprintf("%s/temp_fasta_tabbed.tsv", temp_data_folder)
        temp_splitted_fasta_file <- sprintf("%s/temp_fasta_splitted.tsv", temp_data_folder)
        temp_aminoacid_fasta_file <- sprintf("%s/temp_fasta_aminoacids.tsv", temp_data_folder)
        
        #The output method can be "write" (makes files), "list" (returns a list), or "both"
        predictors_output_method <- "both"
        output_per_aminoacid_file <- sprintf("%s/output_per_aminoacid.tsv", output_data_folder)
        output_per_protein_file <- sprintf("%s/output_per_protein.tsv", output_data_folder)
        output_per_kmer_file <- sprintf("%s/output_per_kmer.tsv", output_data_folder)
        output_per_peptide_file <- sprintf("%s/output_per_peptide.tsv", output_data_folder)
        output_per_repeat_file <- sprintf("%s/output_per_repeat.tsv", output_data_folder)
        
        #NetSurfp uses multiple cores per process. I would give 4 cores to each process
        #This has to be 1 or more (1 meaning just 1 process without parallelization)
        number_of_parallel_processes_for_NetSurfp <- max(1, floor(number_of_parallel_processes / 4))
        
        ## Xstream
        Xstream_path <- "/usr/local/xstream/xstream.jar"
        
        #########################################-
        #### Config - Normalize protein data ####
        #########################################-
        ## General
        source(sprintf("%s/11_normalizeProteinData.R", scripts_folder))
        
        source(sprintf("%s/lib/normalization.R", scripts_folder))
        source(sprintf("%s/lib/auxiliar.R", scripts_folder))
        
        #The output method can be "write" (makes files), "list" (returns a list), or "both"
        protein_normalization_output_method <- "both"
        normalized_output_per_protein_file <- sprintf("%s/normalized_output_per_protein.tsv", output_data_folder)
        
        #########################################-
        #### Config - Normalize peptide data ####
        #########################################-
        ## General
        source(sprintf("%s/12_normalizePeptideData.R", scripts_folder))
        
        source(sprintf("%s/lib/normalization.R", scripts_folder))
        source(sprintf("%s/lib/auxiliar.R", scripts_folder))
        
        #The output method can be "write" (makes files), "list" (returns a list), or "both"
        peptide_normalization_output_method <- "both"
        normalized_output_per_peptide_file <- sprintf("%s/normalized_output_per_peptide.tsv", output_data_folder)
        
        ## Paircoil2
        #Here I'm asking if there is a fragment of length Paircoil2_fragment_length of aa above Paircoil2_threshold
        # Paircoil2_threshold <- 0.5
        
        ###########################################-
        #### Config - Calculate Protein Scores ####
        ###########################################-
        ## General
        source(sprintf("%s/21_calculateScores.R", scripts_folder))
        
        #The linear model file
        protein_model_file <- sprintf("%s/protein_balanced_generic_model.rda", model_data_folder)
        
        #The output method can be "write" (makes files), "list" (returns a list), or "both"
        protein_scores_output_method <- "both"
        score_output_per_protein_file <- sprintf("%s/score_output_per_protein.tsv", output_data_folder)
        
        ###########################################-
        #### Config - Calculate Peptide Scores ####
        ###########################################-
        #The linear model file
        peptide_model_file <- sprintf("%s/peptide_balanced_generic_model.rda", model_data_folder)
        
        #The output method can be "write" (makes files), "list" (returns a list), or "both"
        peptide_scores_output_method <- "both"
        score_output_per_peptide_file <- sprintf("%s/score_output_per_peptide.tsv", output_data_folder)
        
        ##############################-
        #### **CHECK PARAMETERS** ####
        ##############################-
        validateArgument <- function(arg_name, arg_val, comparison, ...) {
            more_args <- list(...)
            
            output <- ""
            if (comparison == ">") {
                value <- more_args[[1]]
                if (arg_val <= value) {
                    output <- sprintf("%s should be greater than %s", arg_name, value)
                }
            }
            else if (comparison == ">=") {
                value <- more_args[[1]]
                if (arg_val < value) {
                    output <- sprintf("%s should be greater or equal than %s", arg_name, value)
                }
            }
            else if (comparison == "<") {
                value <- more_args[[1]]
                if (arg_val >= value) {
                    output <- sprintf("%s should be less than %s", arg_name, value)
                }
            }
            else if (comparison == "<=") {
                value <- more_args[[1]]
                if (arg_val > value) {
                    output <- sprintf("%s should be less or equal than %s", arg_name, value)
                }
            }
            else if (comparison == "<>") {
                min_value <- more_args[[1]]
                max_value <- more_args[[2]]
                if ((arg_val <= min_value) | (arg_val >= max_value)) {
                    output <- sprintf("%s should be between %s and %s (not including them)", arg_name, min_value, max_value)
                }
            }
            else if (comparison == "<=>") {
                min_value <- more_args[[1]]
                max_value <- more_args[[2]]
                if ((arg_val < min_value) | (arg_val > max_value)) {
                    output <- sprintf("%s should be between %s and %s", arg_name, min_value, max_value)
                }
            }
            else if (comparison == "inS") {
                words <- unlist(more_args)
                
                if ((arg_val %in% words) == FALSE) {
                    output <- sprintf("%s has to have one of the following values: %s", arg_name, paste(words, collapse = ", "))
                }
            }
            else if (comparison == ">V") {
                other_arg_name <- more_args[[1]]
                value <- more_args[[2]]
                if (arg_val <= value) {
                    output <- sprintf("%s should be greater than %s", arg_name, other_arg_name)
                }
            }
            else if (comparison == ">=V") {
                other_arg_name <- more_args[[1]]
                value <- more_args[[2]]
                if (arg_val < value) {
                    output <- sprintf("%s should be greater or equal than %s", arg_name, other_arg_name)
                }
            }
            else if (comparison == "<V") {
                other_arg_name <- more_args[[1]]
                value <- more_args[[2]]
                if (arg_val >= value) {
                    output <- sprintf("%s should be less than %s", arg_name, other_arg_name)
                }
            }
            else if (comparison == "<=V") {
                other_arg_name <- more_args[[1]]
                value <- more_args[[2]]
                if (arg_val > value) {
                    output <- sprintf("%s should be less or equal than %s", arg_name, other_arg_name)
                }
            }
            else if (comparison == "<>V") {
                min_other_arg_name <- more_args[[1]]
                min_value <- more_args[[2]]
                max_other_arg_name <- more_args[[3]]
                max_value <- more_args[[4]]
                if ((arg_val <= min_value) | (arg_val >= max_value)) {
                    output <- sprintf("%s should be between %s and %s (not including them)", arg_name, min_other_arg_name, max_other_arg_name)
                }
            }
            else if (comparison == "<=>V") {
                min_other_arg_name <- more_args[[1]]
                min_value <- more_args[[2]]
                max_other_arg_name <- more_args[[3]]
                max_value <- more_args[[4]]
                if ((arg_val < min_value) | (arg_val > max_value)) {
                    output <- sprintf("%s should be between %s and %s", arg_name, min_other_arg_name, max_other_arg_name)
                }
            }
            
            output
        }
        
        #Set use_X based on file path being ""
        if (CrossReactivity_fasta_file == "") {
            use_CrossReactivity <- 0
        }
        if (Coendemicity_fasta_file == "") {
            use_Coendemicity <- 0
        }
        
        #Check that the files exist
        validation <- c()
        if (file.exists(input_fasta_file) == FALSE) {
            validation <- c(validation, sprintf("Couldn't find the file: ", input_fasta_file))
        }
        if ((use_CrossReactivity) & (file.exists(CrossReactivity_fasta_file) == FALSE)) {
            validation <- c(validation, sprintf("Couldn't find the file: ", CrossReactivity_fasta_file))
        }
        if ((use_Coendemicity) & (file.exists(Coendemicity_fasta_file) == FALSE)) {
            validation <- c(validation, sprintf("Couldn't find the file: ", Coendemicity_fasta_file))
        }
        
        #Check the parameters
        validation <- c(validation, validateArgument("Peptide Length", peptide_length, "<=>", 4, 30))
        validation <- c(validation, validateArgument("Peptide Overlap", peptide_overlap, ">", 0))
        validation <- c(validation, validateArgument("Peptide Overlap", peptide_overlap, "<V", "Peptide Length", peptide_length))
        validation <- c(validation, validateArgument("Kmer Length", KmerSimilarity_kmer_length, "<=>", 4, 30))
        validation <- c(validation, validateArgument("NetMHCIIpan - Peptide Length", NetMHCIIpan_binding_peptide_length, "<=>", 9, 30))
        validation <- c(validation, validateArgument("Paircoil2 - Fragment Length", Paircoil2_fragment_length, ">=", 28))
        validation <- c(validation, validateArgument("Paircoil2 - Threshold", Paircoil2_threshold, "<=>", 0, 1))
        validation <- c(validation, validateArgument("SignalP - Organism Group", SignalP_organism_group, "inS", "euk", "gram+", "gram-"))
        validation <- c(validation, validateArgument("Xstream - Min Period", Xstream_min_period, ">", 0))
        validation <- c(validation, validateArgument("Xstream - Min Copy Number", Xstream_min_copy_number, ">", 0))
        validation <- c(validation, validateArgument("Xstream - Max Consensus Error", Xstream_max_consensus_error, "<=>", 0, 1))
        validation <- c(validation, validateArgument("Coendemicity - Protein - Min Amount in Coendemic Proteome for Penalty", Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty, ">=", 0))
        validation <- c(validation, validateArgument("Coendemicity - Protein - Start Penalty Proportion", Coendemicity_protein_start_penalty_proportion, "<=>", 0, 1))
        validation <- c(validation, validateArgument("Coendemicity - Protein - Max Penalty Proportion", Coendemicity_protein_max_penalty_proportion, "<=>", 0, 1))
        validation <- c(validation, validateArgument("Coendemicity - Peptide - Min Amount in Coendemic Proteome for Penalty", Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty, ">=", 0))
        validation <- c(validation, validateArgument("Coendemicity - Peptide - Start Penalty Proportion", Coendemicity_peptide_start_penalty_proportion, "<=>", 0, 1))
        validation <- c(validation, validateArgument("Coendemicity - Peptide - Max Penalty Proportion", Coendemicity_peptide_max_penalty_proportion, "<=>", 0, 1))
        
        validation <- validation[validation != ""]
        
        if (length(validation) > 0) {
            writeLines(validation)
            stop("Parameters validation failed")
        }
        
        ##################-
        #### **MAIN** ####
        ##################-
        ########################-
        #### Run Predictors ####
        ########################-
        output_method <- predictors_output_method
        list_output <-
            runPredictors(input_fasta_file = input_fasta_file,
                          temp_data_folder = temp_data_folder, output_data_folder = output_data_folder, scripts_folder = scripts_folder,
                          temp_fasta_file = temp_fasta_file, temp_tabbed_fasta_file = temp_tabbed_fasta_file, temp_splitted_fasta_file = temp_splitted_fasta_file, temp_aminoacid_fasta_file = temp_aminoacid_fasta_file,
                          output_method = output_method,
                          output_per_aminoacid_file = output_per_aminoacid_file, output_per_protein_file = output_per_protein_file, output_per_kmer_file = output_per_kmer_file, output_per_peptide_file = output_per_peptide_file, output_per_repeat_file = output_per_repeat_file,
                          number_of_parallel_processes = number_of_parallel_processes, number_of_parallel_processes_for_NetSurfp = number_of_parallel_processes_for_NetSurfp,
                          peptide_length = peptide_length, peptide_overlap = peptide_overlap,
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
        output_method <- protein_normalization_output_method
        list_output <-
            normalizeProteinData(output_per_aminoacid = output_per_aminoacid, output_per_protein = output_per_protein, output_per_kmer = output_per_kmer, output_per_peptide = output_per_peptide, output_per_repeat = output_per_repeat,
                                 output_method = output_method,
                                 normalized_output_per_protein_file = normalized_output_per_protein_file,
                                 Paircoil2_fragment_length = Paircoil2_fragment_length, Paircoil2_threshold = Paircoil2_threshold,
                                 Xstream_min_period = Xstream_min_period, Xstream_min_copy_number = Xstream_min_copy_number, Xstream_max_consensus_error = Xstream_max_consensus_error,
                                 Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty = Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_protein_start_penalty_proportion = Coendemicity_protein_start_penalty_proportion, Coendemicity_protein_max_penalty_proportion = Coendemicity_protein_max_penalty_proportion,
                                 use_BepiPred = use_BepiPred, use_Paircoil2 = use_Paircoil2, use_PredGPI = use_PredGPI, use_SignalP = use_SignalP,
                                 use_TMHMM = use_TMHMM, use_NetSurfp = use_NetSurfp,
                                 use_Iupred = use_Iupred, use_NetOglyc = use_NetOglyc, use_Xstream = use_Xstream, use_NetMHCIIpan = use_NetMHCIIpan,
                                 use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight = use_MolecularWeight,
                                 use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
        
        normalized_output_per_protein <- list_output[["normalized_output_per_protein"]]
        rm(list_output)
        gc()
        
        ################################-
        #### Normalize Peptide Data ####
        ################################-
        output_method <- peptide_normalization_output_method
        list_output <-
            normalizePeptideData(output_per_aminoacid = output_per_aminoacid, output_per_protein = output_per_protein, output_per_kmer = output_per_kmer, output_per_peptide = output_per_peptide, output_per_repeat = output_per_repeat,
                                 output_method = output_method,
                                 normalized_output_per_peptide_file = normalized_output_per_peptide_file,
                                 peptide_length = peptide_length, peptide_overlap = peptide_overlap,
                                 Paircoil2_threshold = Paircoil2_threshold,
                                 Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty = Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_peptide_start_penalty_proportion = Coendemicity_peptide_start_penalty_proportion, Coendemicity_peptide_max_penalty_proportion = Coendemicity_peptide_max_penalty_proportion,
                                 use_BepiPred = use_BepiPred, use_Paircoil2 = use_Paircoil2, use_PredGPI = use_PredGPI, use_SignalP = use_SignalP,
                                 use_TMHMM = use_TMHMM, use_NetSurfp = use_NetSurfp,
                                 use_Iupred = use_Iupred, use_NetOglyc = use_NetOglyc, use_Xstream = use_Xstream, use_NetMHCIIpan = use_NetMHCIIpan,
                                 use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight = use_MolecularWeight,
                                 use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity)
        
        normalized_output_per_peptide <- list_output[["normalized_output_per_peptide"]]
        rm(list_output)
        gc()
        
        ##################################-
        #### Calculate Protein Scores ####
        ##################################-
        load(protein_model_file) # > protein_balanced_generic_model
        
        output_method <- protein_scores_output_method
        list_output <-
            calculateProteinScores(normalized_output_per_protein = normalized_output_per_protein,
                                   protein_model = protein_balanced_generic_model,
                                   output_method = output_method, score_output_per_protein_file =  score_output_per_protein_file,
                                   use_BepiPred = use_BepiPred, use_Paircoil2 =  use_Paircoil2, use_PredGPI =  use_PredGPI, use_SignalP =  use_SignalP,
                                   use_TMHMM = use_TMHMM, use_NetSurfp =  use_NetSurfp,
                                   use_Iupred = use_Iupred, use_NetOglyc =  use_NetOglyc, use_Xstream =  use_Xstream, use_NetMHCIIpan =  use_NetMHCIIpan,
                                   use_IsoelectricPoint = use_IsoelectricPoint, use_MolecularWeight =  use_MolecularWeight,
                                   use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity)
        
        score_output_per_protein <- list_output[["score_output_per_protein"]]
        rm(list_output)
        gc()
        
        ##################################-
        #### Calculate Peptide Scores ####
        ##################################-
        load(peptide_model_file) # > peptide_balanced_generic_model
        
        output_method <- peptide_scores_output_method
        list_output <-
            calculatePeptideScores(normalized_output_per_peptide = normalized_output_per_peptide,
                                   peptide_model = peptide_balanced_generic_model,
                                   score_output_per_protein = score_output_per_protein,
                                   output_method = output_method, score_output_per_peptide_file =  score_output_per_peptide_file,
                                   use_BepiPred = use_BepiPred, use_Paircoil2 =  use_Paircoil2, use_PredGPI =  use_PredGPI, use_SignalP =  use_SignalP,
                                   use_TMHMM = use_TMHMM, use_NetSurfp =  use_NetSurfp,
                                   use_Iupred = use_Iupred, use_NetOglyc =  use_NetOglyc, use_Xstream =  use_Xstream, use_NetMHCIIpan =  use_NetMHCIIpan,
                                   use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity)
        
        score_output_per_peptide <- list_output[["score_output_per_peptide"]]
        rm(list_output)
        gc()
        
        #####################-
        #### Remove Temp ####
        #####################-
        if (remove_temp_folder) {
            unlink(temp_data_folder, recursive = T)
        }
    }