# APRANK v1.0 by Alejandro Ricci
# Based on Pepranker v4.0 by Mauricio Brunner, Diego Ramoa, Santiago Carmona and Fernán Agüero

# Needs Perl (Config::General, Bio::SeqIO::fasta),
# R (data.table, foreach, parallel, doParallel),
# BepiPred 1.0, Iupred 1.0, netMHCIIpan 2.0,
# NetOGlyc 3.1d, NetSurfp 1.0, Paircoil2,
# PredGPI 1.4.3, SignalP 4.0, TMHMM 2.0c and Xstream 1.71

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
aprank_folder <- "/home/user/Desktop/APRANK"

##############################-
#### Config - Input Files ####
##############################-
## FASTA File
input_fasta_file <- sprintf("%s/test/CHO.fsa", aprank_folder)

## CrossReactivity (Optional, for Pathogen vs Host)
CrossReactivity_fasta_file <- sprintf("%s/test/test_crossreactivity.fasta", aprank_folder)

## Coendemicity (Optional, for Pathogen vs Other pathogen)
Coendemicity_fasta_file <- sprintf("%s/test/test_coendemicity.fasta", aprank_folder)

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
use_Coendemicity <- 1

##########################-
#### Config - General ####
##########################-
output_data_folder <- sprintf("%s/output", aprank_folder)

peptide_length <- 15
peptide_overlap <- 14

#################################-
#### Config - Run Predictors ####
#################################-
## SignalP
#The organism group of the genome to analyze. It can be: "euk", "gram+" or "gram-".
SignalP_organism_group <- "euk"

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

## General
#This has to be 1 or more (1 meaning just 1 process without parallelization)
number_of_parallel_processes <- 1

#########################################-
#### Config - Normalize protein data ####
#########################################-
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
## Coendemicity
# Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty is the min amount of times a given kmer has to appear in the
# coendemic proteome to apply the penalty to peptides containing it
Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty <- 0
# Coendemicity_peptide_start_penalty_proportion is the starting proportion to apply a penalty to the score
# This penalty goes up until the proportion reaches the Coendemicity_protein_max_penalty_proportion, where the penalty
# becomes 1 and the score 0
Coendemicity_peptide_start_penalty_proportion <- 0
Coendemicity_peptide_max_penalty_proportion <- 2/15

##########################-
#### **RUN PIPELINE** ####
##########################-
source(sprintf("%s/source/00_pipeline_for_R.R", aprank_folder))

runPipeline(aprank_folder = aprank_folder,
            use_BepiPred = use_BepiPred, use_IsoelectricPoint = use_IsoelectricPoint, use_Iupred = use_Iupred, use_MolecularWeight = use_MolecularWeight,
            use_NetMHCIIpan = use_NetMHCIIpan, use_NetOglyc = use_NetOglyc, use_NetSurfp = use_NetSurfp, use_Paircoil2 = use_Paircoil2,
            use_PredGPI = use_PredGPI, use_SignalP = use_SignalP, use_TMHMM = use_TMHMM, use_Xstream = use_Xstream,
            use_SelfSimilarity = use_SelfSimilarity, use_CrossReactivity = use_CrossReactivity, use_Coendemicity = use_Coendemicity,
            output_data_folder = output_data_folder, input_fasta_file = input_fasta_file, CrossReactivity_fasta_file = CrossReactivity_fasta_file, Coendemicity_fasta_file = Coendemicity_fasta_file,
            number_of_parallel_processes = number_of_parallel_processes,
            peptide_length = peptide_length, peptide_overlap = peptide_overlap, KmerSimilarity_kmer_length = KmerSimilarity_kmer_length,
            SignalP_organism_group = SignalP_organism_group,
            NetMHCIIpan_binding_peptide_length = NetMHCIIpan_binding_peptide_length, NetMHCIIPan_alleles = NetMHCIIPan_alleles,
            Paircoil2_fragment_length = Paircoil2_fragment_length, Paircoil2_threshold = Paircoil2_threshold,
            Xstream_min_period = Xstream_min_period, Xstream_min_copy_number = Xstream_min_copy_number, Xstream_max_consensus_error = Xstream_max_consensus_error,
            Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty = Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_protein_start_penalty_proportion = Coendemicity_protein_start_penalty_proportion, Coendemicity_protein_max_penalty_proportion = Coendemicity_protein_max_penalty_proportion,
            Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty = Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty, Coendemicity_peptide_start_penalty_proportion = Coendemicity_peptide_start_penalty_proportion, Coendemicity_peptide_max_penalty_proportion)
