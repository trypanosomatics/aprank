# APRANK v1.0

## DESCRIPTION

This program takes a Fasta file as input and analyzes several aspects for each protein which are then translated into a score that is related with 
the probability of that protein being antigenic

## AUTHORS

APRANK v1.0 by Alejandro Ricci

Based on Pepranker v4.0 by Mauricio Brunner, Diego Ramoa, Santiago Carmona and Fernán Agüero

## DEPENDENCES

Needs Perl (Config::General, Bio::SeqIO::fasta, Parallel::ForkManager), R (data.table, pROC, ggplot2), BepiPred 1.0, Iupred 1.0, netMHCIIpan 2.0, NetOGlyc 3.1d, NetSurfp 1.0, Paircoil2, PredGPI 1.4.3, SignalP 4.0, TMHMM 2.0c and Xstream 1.71

## GETTING THE MODELS

Due to their size, the models aren't present in this repository. You will need to download them from [here](https://mega.nz/folder/3U00RLBJ#8H8hGMuLVuEn_SaqXG1bRw) and place them inside the ***/models*** folder.

## INSTALLING THE PREDICTORS

APRANK requires many third-party predictors to work. See the **INSTALL** file for more detailed guide on how to get, install, configure and test each of the predictors needed.

## RUNNING APRANK VIA R

Edit "R_pipeline.R" using RStudio or other text editor. Adjust the config variables as necessary.

Either run the code in RStudio or run the following command in the APRANK folder
```bash
Rscript --vanilla R_pipeline.R
```

## RUNNING APRANK VIA PERL

Edit the "Perl Pipeline Config.txt" file using a text editor. Adjust the config variables as necessary.

Run the "pipeline.pl" from console, passing the input file
```bash
perl pipeline.pl test/CHO.fsa
```

Remember that if you are using the Perl pipeline you can also change settings directly from console, such as
```bash
perl pipeline.pl test/CHO.fsa -l 16 -o 12
```

### PERL COMMAND LINE PARAMETERS
```bash
pipeline.pl fasta_file [-hf host_fasta_file] [-cf coendemic_fasta_file] [-h|help] [-l peptide_length] [-o peptide_overlap]
    [-kl kmer_length] [-mpr protein_model_file] [-mpe peptide_model_file] [-of output_folder] [-bp use_BepiPred]
    [-ip use_IsoelectricPoint] [-ir use_Iupred] [-mw use_MolecularWeight] [-mhc use_NetMHCIIpan] [-no use_NetOglyc] [-ns use_NetSurfp]
    [-p2 use_Paircoil2] [-gpi use_PredGPI] [-ss use_SelfSimilarity] [-sp use_SignalP] [-tm use_TMHMM] [-xs use_Xstream]
    [-cem Coendemicity_max_matches_allowed] [-cep Coendemicity_penalty_per_match] [-mhcl NetMHCIIpan_binding_peptide_length]
    [-mhca NetMHCIIpan_alleles] [-p2l Paircoil2_fragment_length] [-p2t Paircoil2_threshold]
    [-spg SignalP_organism_group] [-xsmp Xstream_min_period] [-xsmc Xstream_min_copy_number] [-xsme Xstream_max_consensus_error]
```

The meaning of each parameter can be seen in the **README** file. It's recomended to change the parameters values from **Perl Pipeline Config.txt**.
