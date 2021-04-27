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

Due to their size, the models aren't present in this repository. You will need to download them from [Dryad](https://doi.org/10.5061/dryad.zcrjdfnb1) and place them inside the ***/models*** folder.

You should end up with the following files:
- /models/peptide_balanced_generic_model.rda
- /models/protein_balanced_generic_model.rda

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

The meaning of each parameter can be seen in the **USAGE** file. It's recomended to change the parameters values from **Perl Pipeline Config.txt**.

## READING THE OUTPUTS

By default, APRANK will save the outputs in the ***/outputs*** folder, inside a new folder with the *FastaName_Date* format.

APRANK will produce two set of scores, one for the proteins and one for the peptides. All scores go from 0 (predicted to be non-antigenic) to 1 (predicted to be antigenic).

The data for the protein predictions can be found inside the **score_output_per_protein.tsv** file, and its columns are:

- **id:** an internal id assigned ID APRANK to deal with possible repeated proteins IDs
- **original_id:** the original ID found in the FASTA for this protein
- **protein_score:** the antigenicity score predicted by APRANK

The data for the protein predictions can be found inside the **score_output_per_peptide.tsv** file, and its columns are:

- **id:** the corresponding protein ID
- **original_id:** the corresponding protein original ID
- **start:** the position in the protein of the peptide's first amino acid
- **end:** the position in the protein of the peptide's last amino acid
- **peptide:** the peptide's sequence
- **peptide_score:** the antigenicity score predicted by APRANK
- **protein_score:** the corresponding protein antigenicity score
- **combined_score:** the combined antigenicity score, including info of both protein and peptide scores

Also, if calculating the Coendemicity, an extra set of scores will appear with a *_wCoendemicPenalty* suffix in the name.

APRANK will also save the processed data from the predictors inside the ***/predictors_outputs*** folder just in case you need it.
