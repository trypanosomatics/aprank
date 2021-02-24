#!/usr/bin/perl
use warnings;
use strict;

use Bio::SeqIO::fasta;

#Read the arguments
my $num_args = $#ARGV + 1;
if ($num_args != 5) {
    print "\nUsage: split_proteome_into_kmers.pl input.fasta kmer_length max_protein_length replace_nonAA_chars_by output_file\n";
    exit;
}

my $input_file = $ARGV[0];
my $kmer_length = $ARGV[1];
my $max_protein_length = $ARGV[2]; #0 means disabled
my $replace_nonAA_chars_by = $ARGV[3]; #"no" means disabled
my $output_file = $ARGV[4];

#Fetch the kmers found in the genome
my %protein_kmers;
my %total_kmers;    
my $host_proteome = new Bio::SeqIO::fasta(-file => $input_file, '-format' => 'multifasta');
while (my $seq = $host_proteome->next_seq()) {
    my $id = $seq->id();
    my $sequence = $seq->seq();
    my $sequence_length = length($sequence);
    my $position = 0;

    $sequence = uc $sequence;
    if (($max_protein_length > 0) && ($sequence_length > $max_protein_length)) {
        $sequence = substr($sequence, 0, $max_protein_length);
        $sequence_length = length($sequence);
    }
    if ($replace_nonAA_chars_by ne "no") {
        #Remove any strange characters from the sequence
        $sequence =~ s/[^ACDEFGHIKLMNPQRSTVWY]/$replace_nonAA_chars_by/g;    
    }

    #Count the number of times each kmer appears in the sequence
    while ($kmer_length  <= ($sequence_length - $position)) {
        my $fragment = substr($sequence, $position, $kmer_length);
        $fragment = uc($fragment); 

        #Check for the protein
        if (exists $protein_kmers{$id}{$fragment}) {
            $protein_kmers{$id}{$fragment} += 1;
        } else {
            $protein_kmers{$id}{$fragment} = 1;
        }
        #Check for the total kmers
        if (exists $total_kmers{$fragment}) {
            $total_kmers{$fragment} += 1;
        } else {
            $total_kmers{$fragment} = 1;
        }

        $position++;
    }
}

#Write the output at protein level
my $output_aux = "protein\tkmer\tquantity\tquantity_in_proteome\n";
foreach my $protein (sort {lc $a cmp lc $b} keys %protein_kmers) {
    foreach my $fragment (sort {lc $a cmp lc $b} keys %{$protein_kmers{$protein}}) {
        $output_aux .= $protein . "\t" . $fragment . "\t" . $protein_kmers{$protein}{$fragment} . "\t" . $total_kmers{$fragment} . "\n";
    }
}

open(OUTPUT, ">" . "$output_file") || die "Could not open output file\n";
print OUTPUT $output_aux;
close(OUTPUT);
