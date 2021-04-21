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

############################-
#### **READ ARGUMENTS** ####
############################-
args = commandArgs(trailingOnly = TRUE)

args_names <- c("aprank_folder",
                "use_BepiPred", "use_IsoelectricPoint", "use_Iupred", "use_MolecularWeight",
                "use_NetMHCIIpan", "use_NetOglyc", "use_NetSurfp", "use_Paircoil2",
                "use_PredGPI", "use_SignalP", "use_TMHMM", "use_Xstream",
                "use_SelfSimilarity", "use_CrossReactivity", "use_Coendemicity",
                "output_data_folder", "input_fasta_file", "CrossReactivity_fasta_file", "Coendemicity_fasta_file",
                "number_of_parallel_processes",
                "peptide_length", "peptide_overlap", "KmerSimilarity_kmer_length",
                "SignalP_organism_group",
                "NetMHCIIpan_binding_peptide_length", "NetMHCIIPan_alleles",
                "Paircoil2_fragment_length", "Paircoil2_threshold",
                "Xstream_min_period", "Xstream_min_copy_number", "Xstream_max_consensus_error",
                "Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty", "Coendemicity_protein_start_penalty_proportion", "Coendemicity_protein_max_penalty_proportion",
                "Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty", "Coendemicity_peptide_start_penalty_proportion", "Coendemicity_peptide_max_penalty_proportion")

integer_args_names <- c("use_BepiPred", "use_IsoelectricPoint", "use_Iupred", "use_MolecularWeight",
                        "use_NetMHCIIpan", "use_NetOglyc", "use_NetSurfp", "use_Paircoil2",
                        "use_PredGPI", "use_SignalP", "use_TMHMM", "use_Xstream",
                        "use_SelfSimilarity", "use_CrossReactivity", "use_Coendemicity",
                        "number_of_parallel_processes",
                        "peptide_length", "peptide_overlap", "KmerSimilarity_kmer_length",
                        "NetMHCIIpan_binding_peptide_length",
                        "Paircoil2_fragment_length",
                        "Xstream_min_period", "Xstream_min_copy_number",
                        "Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty",
                        "Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty")

numeric_args_names <- c("Paircoil2_threshold",
                        "Xstream_max_consensus_error",
                        "Coendemicity_protein_start_penalty_proportion", "Coendemicity_protein_max_penalty_proportion",
                        "Coendemicity_peptide_start_penalty_proportion", "Coendemicity_peptide_max_penalty_proportion")

if (length(args) != length(args_names)) { stop("Wrong number of arguments") }

args_values <- list()
for (i in 1:length(args_names)) {
    arg_for <- args_names[i]
  
    if (arg_for %in% integer_args_names) {
        args_values[[arg_for]] <- as.numeric(args[i])
    } else if (arg_for %in% numeric_args_names) {
        args_values[[arg_for]] <- as.numeric(args[i])
    } else {
        args_values[[arg_for]] <- args[i]
    }
    # print(sprintf("%s | %s", arg_for, args[i]))
}

###############################-
#### **PROCESS ARGUMENTS** ####
###############################-
args_values[["NetMHCIIPan_alleles"]] <- unlist(strsplit(args_values[["NetMHCIIPan_alleles"]], " "))

##########################-
#### **RUN PIPELINE** ####
##########################-
source(sprintf("%s/source/00_pipeline_for_R.R", args_values[["aprank_folder"]]))

runPipeline(aprank_folder = args_values[["aprank_folder"]],
            use_BepiPred = args_values[["use_BepiPred"]], use_IsoelectricPoint = args_values[["use_IsoelectricPoint"]], use_Iupred = args_values[["use_Iupred"]], use_MolecularWeight = args_values[["use_MolecularWeight"]],
            use_NetMHCIIpan = args_values[["use_NetMHCIIpan"]], use_NetOglyc = args_values[["use_NetOglyc"]], use_NetSurfp = args_values[["use_NetSurfp"]], use_Paircoil2 = args_values[["use_Paircoil2"]],
            use_PredGPI = args_values[["use_PredGPI"]], use_SignalP = args_values[["use_SignalP"]], use_TMHMM = args_values[["use_TMHMM"]], use_Xstream = args_values[["use_Xstream"]],
            use_SelfSimilarity = args_values[["use_SelfSimilarity"]], use_CrossReactivity = args_values[["use_CrossReactivity"]], use_Coendemicity = args_values[["use_Coendemicity"]],
            output_data_folder = args_values[["output_data_folder"]], input_fasta_file = args_values[["input_fasta_file"]], CrossReactivity_fasta_file = args_values[["CrossReactivity_fasta_file"]], Coendemicity_fasta_file = args_values[["Coendemicity_fasta_file"]],
            number_of_parallel_processes = args_values[["number_of_parallel_processes"]],
            peptide_length = args_values[["peptide_length"]], peptide_overlap = args_values[["peptide_overlap"]], KmerSimilarity_kmer_length = args_values[["KmerSimilarity_kmer_length"]],
            SignalP_organism_group = args_values[["SignalP_organism_group"]],
            NetMHCIIpan_binding_peptide_length = args_values[["NetMHCIIpan_binding_peptide_length"]], NetMHCIIPan_alleles = args_values[["NetMHCIIPan_alleles"]],
            Paircoil2_fragment_length = args_values[["Paircoil2_fragment_length"]], Paircoil2_threshold = args_values[["Paircoil2_threshold"]],
            Xstream_min_period = args_values[["Xstream_min_period"]], Xstream_min_copy_number = args_values[["Xstream_min_copy_number"]], Xstream_max_consensus_error = args_values[["Xstream_max_consensus_error"]],
            Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty = args_values[["Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty"]], Coendemicity_protein_start_penalty_proportion = args_values[["Coendemicity_protein_start_penalty_proportion"]], Coendemicity_protein_max_penalty_proportion = args_values[["Coendemicity_protein_max_penalty_proportion"]],
            Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty = args_values[["Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty"]], Coendemicity_peptide_start_penalty_proportion = args_values[["Coendemicity_peptide_start_penalty_proportion"]], Coendemicity_peptide_max_penalty_proportion = args_values[["Coendemicity_peptide_max_penalty_proportion"]])