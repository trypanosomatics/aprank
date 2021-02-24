#!/usr/bin/perl
use warnings;
use strict;

use Bio::SeqIO::fasta;

#Read the arguments
my $num_args = $#ARGV + 1;
if ($num_args != 9) {
    print "\nUsage: tab_and_split_proteome.pl input.fasta sequence_length sequence_overlap max_protein_length replace_nonAA_chars_by output_file tabbed_output_file splitted_output_file aminoacid_output_file\n";
    exit;
}

my $input_file = $ARGV[0];
my $sequence_length = $ARGV[1];
my $sequence_overlap = $ARGV[2];
my $max_protein_length = $ARGV[3]; #0 means disabled
my $replace_nonAA_chars_by = $ARGV[4]; #"no" means disabled
my $output_file = $ARGV[5];
my $tabbed_output_file = $ARGV[6];
my $splitted_output_file = $ARGV[7];
my $aminoacid_output_file = $ARGV[8];

#Read the new file as a proteome using the Bio package
my $proteome = new Bio::SeqIO::fasta(-file => "$input_file", '-format' => 'multifasta');

#Loop for each protein and process them
my %used_IDs;
my $output = "";
my $output_tab = qq("protein"\t"original_protein"\t"sequence"\n);
my $output_splitted = qq("protein"\t"sequence"\t"start"\t"end"\n);
my $output_aminoacid = qq("protein"\t"pos"\t"aa"\n);
my $new_protein_ID_index = 1;
while (my $seq = $proteome->next_seq()) {
    my $ID = $seq -> id();
    my $sequence = $seq -> seq();    

    $sequence = uc $sequence;
    if (($max_protein_length > 0) && (length($sequence) > $max_protein_length)) {
        $sequence = substr($sequence, 0, $max_protein_length);
    }
    if ($replace_nonAA_chars_by ne "no") {
        #Remove any strange characters from the sequence
        $sequence =~ s/[^ACDEFGHIKLMNPQRSTVWY]/$replace_nonAA_chars_by/g;    
    }
    
    my $new_protein_ID = "p_$new_protein_ID_index";
    $new_protein_ID_index++;

    $output .= ">$new_protein_ID\n$sequence\n";
    $output_tab .= "\"$new_protein_ID\"\t\"$ID\"\t\"$sequence\"\n";
    $output_splitted .= splitProtein($new_protein_ID, $sequence, $sequence_length, $sequence_overlap);
    $output_aminoacid .= splitProteinByAminoacid($new_protein_ID, $sequence);
}

open(OUTPUT, ">" . "$output_file") || die "Could not open output file\n";
print OUTPUT $output;
close OUTPUT;

open(OUTPUT, ">" . "$tabbed_output_file") || die "Could not open output file\n";
print OUTPUT $output_tab;
close OUTPUT;

open(OUTPUT, ">" . "$splitted_output_file") || die "Could not open output file\n";
print OUTPUT $output_splitted;
close OUTPUT;

open(OUTPUT, ">" . "$aminoacid_output_file") || die "Could not open output file\n";
print OUTPUT $output_aminoacid;
close OUTPUT;

########################
### Auxiliar Functions
########################
sub splitProtein {
    #Arguments => ID sequence sequence_length sequence_overlap
    my ($ID, $sequence, $sequence_length, $sequence_overlap) = @_;

    $sequence = uc $sequence;

    my $sequence_aux = "";
    my $sequence_length_aux = 0;
    my $start = 1;
    my $end;
    my $output = "";
    foreach my $char (split //, $sequence) {
        if ($char ne "") {
            $sequence_aux .= $char;
            $sequence_length_aux++;

            if ($sequence_length_aux == $sequence_length) {
                $end = $start + $sequence_length - 1;
                $output .= "\"" . $ID . "\"\t\"" . $sequence_aux . "\"\t" . $start . "\t" . $end . "\n";

                $sequence_aux = substr($sequence_aux, $sequence_length - $sequence_overlap);
                $sequence_length_aux -= $sequence_length - $sequence_overlap;
                $start += $sequence_length - $sequence_overlap;
            }
        }        
    }
    
    return $output;
}
sub splitProteinByAminoacid {
    #Arguments => ID sequence
    my ($ID, $sequence) = @_;

    $sequence = uc $sequence;

    my $pos = 1;
    my $output = "";
    foreach my $char (split //, $sequence) {
        if ($char ne "") {
            $output .= "\"" . $ID . "\"\t" . $pos . "\t\"" . $char . "\"\n";
            $pos++;
        }        
    }
    
    return $output;
}