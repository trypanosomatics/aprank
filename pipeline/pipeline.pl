#!/usr/bin/perl
use warnings; #remove once finished?
use strict;

use Getopt::Long qw(GetOptions); #to get arguments from console
Getopt::Long::Configure("pass_through"); #to get arguments from console
use File::Basename qw(basename dirname); #to get the file name without PATH_CLUSTER and to get the dirname
use File::Spec qw(rel2abs); #to get the absolute path of a file

=pod

=head1 NAME
    APRANK v1.0
 
=head1
USAGE
    pipeline.pl fasta_file [-h|help] [-hf host_fasta_file]
        [-cf coendemic_fasta_file] [-l peptide_length] [-o peptide_overlap]
        [-kl kmer_length] [-of output_folder]
        [-pp number_of_parallel_processes] [-bp use_BepiPred]
        [-ip use_IsoelectricPoint] [-ir use_Iupred] [-mw use_MolecularWeight]
        [-mhc use_NetMHCIIpan] [-no use_NetOglyc] [-ns use_NetSurfp]
        [-p2 use_Paircoil2] [-gpi use_PredGPI] [-sp use_SignalP]
        [-tm use_TMHMM] [-xs use_Xstream] [-ss use_SelfSimilarity]
        [-cr use_CrossReactivity] [-ce use_Coendemicity]
        [-mhcl NetMHCIIpan_binding_peptide_length] [-mhca NetMHCIIpan_alleles]
        [-p2l Paircoil2_fragment_length] [-p2t Paircoil2_threshold]
        [-spg SignalP_organism_group] [-xsp Xstream_min_period]
        [-xscn Xstream_min_copy_number] [-xse Xstream_max_consensus_error]
        [-cea Coendemicity_protein_min_kmer_amount_for_penalty]
        [-cesp Coendemicity_protein_start_penalty_proportion]
        [-cemp Coendemicity_protein_max_penalty_proportion]
        [-cepepa Coendemicity_peptide_min_kmer_amount_for_penalty]
        [-cepepsp Coendemicity_peptide_start_penalty_proportion]
        [-cepepmp Coendemicity_peptide_max_penalty_proportion]


=head1
DEPENDENCES
    Needs Perl (Config::General, Bio::SeqIO::fasta), R (data.table, foreach, 
    parallel, doParallel), BepiPred 1.0, Iupred 1.0, netMHCIIpan 2.0, 
    NetOGlyc 3.1d, NetSurfp 1.0, Paircoil2, PredGPI 1.4.3, SignalP 4.0, TMHMM 
    2.0c and Xstream 1.71


=head1
DESCRIPTION
    This program takes a Fasta file and analyze several aspects for each 
    protein which are then translated to a score that is related with the 
    probability of that protein being antigenic


=head1
AUTHORS
    APRANK v1.0 by Alejandro Ricci
    Based on Pepranker v4.0 by Mauricio Brunner, Diego Ramoa, Santiago Carmona and Fernán Agüero

=cut

################################################
### Internal Config
################################################
my $main_dir_aux;
my %internal_config;
#The BEGIN section is because some of this variables are needed on load (the lib paths)
BEGIN {
    $main_dir_aux = dirname(File::Spec->rel2abs($0));

    $internal_config{"MAIN_DIR"}                           = $main_dir_aux;

    $internal_config{"LIB_LOCATION"}                       = $main_dir_aux . "/source/perl_pipeline";

    $internal_config{"CONFIG_FILE"}                        = $main_dir_aux . "/Perl Pipeline Config.txt";
    $internal_config{"HELP_TEXT_FILE"}                     = $main_dir_aux . "/source/perl_pipeline/help_text.txt";
}

#Set the quote character
my $q = "'";

################################################
### Modules
################################################
my $eval_config_general = eval {
    require Config::General; #to read default config values from file
    Config::General->import();
    1;
};
if (!$eval_config_general) { die "Perl's module Config::General is needed for this program to work.\n"; }

#Custom Modules
use lib $internal_config{"LIB_LOCATION"};

use My::Console qw(readOptionalArgument readArgumentOrConfig);

################################################
### Read and validate the arguments
################################################
#Get the default values from the config file
my $conf = Config::General->new($internal_config{"CONFIG_FILE"});
my %config = $conf->getall();

my $fasta_file;
my $help                                             = readOptionalArgument(\@ARGV, "Help"                                                                 , "h|help" , ""                                                                );
#Additional Input Files
my $host_fasta_file                                  = readArgumentOrConfig(\@ARGV, "Host FASTA File"                                                      , "hf"     , "s" , \%config, "host_fasta_file"                                 );
my $coendemic_fasta_file                             = readArgumentOrConfig(\@ARGV, "Coendemic Organism FASTA File"                                        , "cf"     , "s" , \%config, "coendemic_fasta_file"                            );
#General
my $peptide_length                                   = readArgumentOrConfig(\@ARGV, "Peptide Length"                                                       , "l"      , "i" , \%config, "peptide_length"                                  );
my $peptide_overlap                                  = readArgumentOrConfig(\@ARGV, "Peptide Overlap"                                                      , "o"      , "i" , \%config, "peptide_overlap"                                 );
my $kmer_length                                      = readArgumentOrConfig(\@ARGV, "Kmer Length"                                                          , "kl"     , "i" , \%config, "kmer_length"                                     );
my $output_folder                                    = readArgumentOrConfig(\@ARGV, "Output Folder Name"                                                   , "of"     , "s" , \%config, "output_folder"                                   );
my $number_of_parallel_processes                     = readArgumentOrConfig(\@ARGV, "Number of Parallel Processes"                                         , "pp"     , "i" , \%config, "number_of_parallel_processes"                    );
#Predictors
my $use_BepiPred                                     = readArgumentOrConfig(\@ARGV, "Use BepiPred"                                                         , "bp"     , "bi", \%config, "use_BepiPred"                                    );
my $use_IsoelectricPoint                             = readArgumentOrConfig(\@ARGV, "Use IsoelectricPoint"                                                 , "ip"     , "bi", \%config, "use_IsoelectricPoint"                            );
my $use_Iupred                                       = readArgumentOrConfig(\@ARGV, "Use Iupred"                                                           , "ir"     , "bi", \%config, "use_Iupred"                                      );
my $use_MolecularWeight                              = readArgumentOrConfig(\@ARGV, "Use MolecularWeight"                                                  , "mw"     , "bi", \%config, "use_MolecularWeight"                             );
my $use_NetMHCIIpan                                  = readArgumentOrConfig(\@ARGV, "Use NetMHCIIpan"                                                      , "mhc"    , "bi", \%config, "use_NetMHCIIpan"                                 );
my $use_NetOglyc                                     = readArgumentOrConfig(\@ARGV, "Use NetOglyc"                                                         , "no"     , "bi", \%config, "use_NetOglyc"                                    );
my $use_NetSurfp                                     = readArgumentOrConfig(\@ARGV, "Use NetSurfp"                                                         , "ns"     , "bi", \%config, "use_NetSurfp"                                    );
my $use_Paircoil2                                    = readArgumentOrConfig(\@ARGV, "Use Paircoil2"                                                        , "p2"     , "bi", \%config, "use_Paircoil2"                                   );
my $use_PredGPI                                      = readArgumentOrConfig(\@ARGV, "Use PredGPI"                                                          , "gpi"    , "bi", \%config, "use_PredGPI"                                     );
my $use_SignalP                                      = readArgumentOrConfig(\@ARGV, "Use SignalP"                                                          , "sp"     , "bi", \%config, "use_SignalP"                                     );
my $use_TMHMM                                        = readArgumentOrConfig(\@ARGV, "Use TMHMM"                                                            , "tm"     , "bi", \%config, "use_TMHMM"                                       );
my $use_Xstream                                      = readArgumentOrConfig(\@ARGV, "Use Xstream"                                                          , "xs"     , "bi", \%config, "use_Xstream"                                     );
#Custom Predictors
my $use_SelfSimilarity                               = readArgumentOrConfig(\@ARGV, "Use SelfSimilarity"                                                   , "ss"     , "bi", \%config, "use_SelfSimilarity"                              );
my $use_CrossReactivity                              = readArgumentOrConfig(\@ARGV, "Use CrossReactivity"                                                  , "cr"     , "bi", \%config, "use_CrossReactivity"                             );
my $use_Coendemicity                                 = readArgumentOrConfig(\@ARGV, "Use Coendemicity"                                                     , "ce"     , "bi", \%config, "use_Coendemicity"                                );
#NetMHCIIpan
my $NetMHCIIpan_binding_peptide_length               = readArgumentOrConfig(\@ARGV, "NetMHCIIpan Peptide Length"                                           , "mhcl"   , "i" , \%config, "NetMHCIIpan_binding_peptide_length"              );
my $NetMHCIIpan_alleles                              = readArgumentOrConfig(\@ARGV, "NetMHCIIpan Alleles"                                                  , "mhca"   , "s" , \%config, "NetMHCIIpan_alleles"                             );
#Paircoil2
my $Paircoil2_fragment_length                        = readArgumentOrConfig(\@ARGV, "Paircoil2 Fragment Length"                                            , "p2l"    , "i" , \%config, "Paircoil2_fragment_length"                       );
my $Paircoil2_threshold                              = readArgumentOrConfig(\@ARGV, "Paircoil2 Threshold"                                                  , "p2t"    , "f" , \%config, "Paircoil2_threshold"                             );
#SignalP
my $SignalP_organism_group                           = readArgumentOrConfig(\@ARGV, "SignalP Organism Group"                                               , "spg"    , "s" , \%config, "SignalP_organism_group"                          );
#Xstream
my $Xstream_min_period                               = readArgumentOrConfig(\@ARGV, "Xstream Min Period"                                                   , "xsp"    , "i" , \%config, "Xstream_min_period"                              );
my $Xstream_min_copy_number                          = readArgumentOrConfig(\@ARGV, "Xstream Min Copy Number"                                              , "xscn"   , "i" , \%config, "Xstream_min_copy_number"                         );
my $Xstream_max_consensus_error                      = readArgumentOrConfig(\@ARGV, "Xstream Max Consensus Error"                                          , "xse"    , "f" , \%config, "Xstream_max_consensus_error"                     );
#Coendemicity
my $Coendemicity_protein_min_kmer_amount_for_penalty = readArgumentOrConfig(\@ARGV, "Coendemicity - Protein - Min Amount in Coendemic Proteome for Penalty", "cea"    , "i" , \%config, "Coendemicity_protein_min_kmer_amount_for_penalty");
my $Coendemicity_protein_start_penalty_proportion    = readArgumentOrConfig(\@ARGV, "Coendemicity - Protein - Start Penalty Proportion"                    , "cesp"   , "f" , \%config, "Coendemicity_protein_start_penalty_proportion"   );
my $Coendemicity_protein_max_penalty_proportion      = readArgumentOrConfig(\@ARGV, "Coendemicity - Protein - Max Penalty Proportion"                      , "cemp"   , "f" , \%config, "Coendemicity_protein_max_penalty_proportion"     );
my $Coendemicity_peptide_min_kmer_amount_for_penalty = readArgumentOrConfig(\@ARGV, "Coendemicity - Peptide - Min Amount in Coendemic Proteome for Penalty", "cepepa" , "i" , \%config, "Coendemicity_peptide_min_kmer_amount_for_penalty");
my $Coendemicity_peptide_start_penalty_proportion    = readArgumentOrConfig(\@ARGV, "Coendemicity - Peptide - Start Penalty Proportion"                    , "cepepsp", "f" , \%config, "Coendemicity_peptide_start_penalty_proportion"   );
my $Coendemicity_peptide_max_penalty_proportion      = readArgumentOrConfig(\@ARGV, "Coendemicity - Peptide - Max Penalty Proportion"                      , "cepepmp", "f" , \%config, "Coendemicity_peptide_max_penalty_proportion"     );

#Flags
if ($help) {
    open(FILE, $internal_config{"HELP_TEXT_FILE"});
    print(<FILE>);
    print("\n\n");
    close(FILE);    
    exit;
}

#Every argument should had been removed from the array by now except $fasta_file
if (scalar(@ARGV) > 1) { die "Too many arguments\n"; }

($fasta_file) = @ARGV;

################################################
### Run APRANK
################################################
system("Rscript --vanilla " . $internal_config{"LIB_LOCATION"} . "/00_pipeline_for_perl.R" . 
	" $q" . $internal_config{"MAIN_DIR"} . "$q" . 
    " " . $use_BepiPred . " " . $use_IsoelectricPoint . " " . $use_Iupred . " " . $use_MolecularWeight .
    " " . $use_NetMHCIIpan . " " . $use_NetOglyc . " " . $use_NetSurfp . " " . $use_Paircoil2 .
    " " . $use_PredGPI . " " . $use_SignalP . " " . $use_TMHMM . " " . $use_Xstream .
    " " . $use_SelfSimilarity . " " . $use_CrossReactivity . " " . $use_Coendemicity .
    " $q" . $output_folder . "$q $q" . $fasta_file . "$q $q" . $host_fasta_file . "$q $q" . $coendemic_fasta_file . "$q" .
    " " . $number_of_parallel_processes .
	" " . $peptide_length . " " . $peptide_overlap . " " . $kmer_length .
    " $q" . $SignalP_organism_group . "$q" .
    " " . $NetMHCIIpan_binding_peptide_length . " $q" . $NetMHCIIpan_alleles . "$q" .
    " " . $Paircoil2_fragment_length . " " . $Paircoil2_threshold .
    " " . $Xstream_min_period . " " . $Xstream_min_copy_number . " " . $Xstream_max_consensus_error .
    " " . $Coendemicity_protein_min_kmer_amount_for_penalty . " " . $Coendemicity_protein_start_penalty_proportion . " " . $Coendemicity_protein_max_penalty_proportion .
    " " . $Coendemicity_peptide_min_kmer_amount_for_penalty . " " . $Coendemicity_peptide_start_penalty_proportion . " " . $Coendemicity_peptide_max_penalty_proportion);

#"aprank_folder",
#"use_BepiPred", "use_IsoelectricPoint", "use_Iupred", "use_MolecularWeight",
#"use_NetMHCIIpan", "use_NetOglyc", "use_NetSurfp", "use_Paircoil2",
#"use_PredGPI", "use_SignalP", "use_TMHMM", "use_Xstream",
#"use_SelfSimilarity", "use_CrossReactivity", "use_Coendemicity",
#"output_data_folder", "input_fasta_file", "CrossReactivity_fasta_file", "Coendemicity_fasta_file",
#"number_of_parallel_processes",
#"peptide_length", "peptide_overlap", "KmerSimilarity_kmer_length",
#"SignalP_organism_group",
#"NetMHCIIpan_binding_peptide_length", "NetMHCIIpan_alleles",
#"Paircoil2_fragment_length", "Paircoil2_threshold",
#"Xstream_min_period", "Xstream_min_copy_number", "Xstream_max_consensus_error",
#"Coendemicity_protein_min_amount_in_coendemic_proteome_for_penalty", "Coendemicity_protein_start_penalty_proportion", "Coendemicity_protein_max_penalty_proportion",
#"Coendemicity_peptide_min_amount_in_coendemic_proteome_for_penalty", "Coendemicity_peptide_start_penalty_proportion", "Coendemicity_peptide_max_penalty_proportion"