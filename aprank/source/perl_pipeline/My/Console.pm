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

package My::Console;

use strict;
use warnings;

use Exporter qw(import);

=pod

=head1 NAME
    Console Module

=head1
DEPENDENCES
    Needs Perl (Modules: My::Format)

=head1
DESCRIPTION
    Module with Console related functions for Perl

=head1
AUTHORS
    Console.pm by Alejandro Ricci

=cut

our @EXPORT_OK = qw(readArgument
                    readOptionalArgument
                    readArgumentOrConfig
                    validateArgument
                    );

use My::Format qw(checkType checkTypeOutput);

sub readArgument {
    #Arguments => \@ARGV arg_name console_code type
    my ($ARGV_ref, $arg_name, $console_code, $type) = @_;

    return readArgumentAux($ARGV_ref, $arg_name, $console_code, $type, 0);
}
sub readOptionalArgument {
    #Arguments => \@ARGV arg_name console_code type
    my ($ARGV_ref, $arg_name, $console_code, $type) = @_;

    return readArgumentAux($ARGV_ref, $arg_name, $console_code, $type, 1);
}
sub readArgumentAux {
    #Arguments => \@ARGV arg_name console_code type
    my ($ARGV_ref, $arg_name, $console_code, $type, $optional) = @_;

    my %error_message = checkTypeOutput();

    my $output;
    my $i = 0;
    while ($i < scalar(@$ARGV_ref)) {
        my $argFor = @$ARGV_ref[$i];

        if (substr($argFor, 0, 1) eq "-") {
            #Check if is a console command
            my $argv_code = substr($argFor, 1);

            my @codes = split(/\|/, $console_code); #to allow multiple codes
            foreach my $codeFor (@codes) {
                $codeFor = "^" . $codeFor . "\$";
            }
            my $match = join("|", @codes);

            if ($argv_code =~ m/$match/g) {
                #I've found the code I was looking for
                if ($type ne "") {
                    #Check if there is no next index
                    if (($i + 1) >= scalar(@$ARGV_ref)) { die "Value expected for $arg_name\n"; }

                    my $valueFor = @$ARGV_ref[$i + 1];

                    #Give an error if the value type isn't the expected one
                    #Types are i (int), f (float), bi (boolean integer), c (char), s (string)
                    unless(checkType($valueFor, $type)) { die $error_message{$type} . " expected for $arg_name\n"; }

                    $output = $valueFor;

                    splice(@$ARGV_ref, $i, 2); #2 to delete both the code and the value
                } else {
                    #If type is "" then that means is a flag argument, so just return 1
                    $output = 1;

                    splice(@$ARGV_ref, $i, 1);
                }

                last;
            }
        }

        $i += 1;
    }

    if (($optional == 0) && (not defined $output)) {
        die "Missing argument: $arg_name\n";
    }

    return $output;
}
sub readArgumentOrConfig {
    #Arguments => \@ARGV arg_name console_code type \%config config_name
    my ($ARGV_ref, $arg_name, $console_code, $type, $config_ref, $config_name) = @_;

    my %error_message = checkTypeOutput();

    my $output;

    #This check allows arguments that only have Config
    if ($console_code ne "") {
        #Here I pass optional always so it doesn't give the error without first checking the config file
        $output = readArgumentAux($ARGV_ref, $arg_name, $console_code, $type, 1);
    }

    if ((not defined $output) && (exists $$config_ref{$config_name})) {
        $output = $$config_ref{$config_name};

        #Give an error if the value type isn't the expected one
        #Types are i (int), f (float), bi (boolean integer), c (char), s (string)
        unless(checkType($output, $type)) { die $error_message{$type} . " expected for $arg_name (value gotten from Config)\n"; }
    }

    if (not defined $output) {
        die "Missing argument: $arg_name\n";
    }    

    return $output;
}
sub validateArgument {
    #Arguments => arg_to_val comparison @more_args
    my ($arg_name, $arg_val, $comparison, @more_args) = @_;

    if ($comparison eq ">") {
        my ($value) = @more_args;
        if ($arg_val <= $value) {
            die "$arg_name should be greater than $value\n";
        }
    }
    elsif ($comparison eq ">=") {
        my ($value) = @more_args;
        if ($arg_val < $value) {
            die "$arg_name should be greater or equal than $value\n";
        }
    }       
    elsif ($comparison eq "<") {
        my ($value) = @more_args;
        if ($arg_val >= $value) {
            die "$arg_name should be less than $value\n";
        }
    }   
    elsif ($comparison eq "<=") {
        my ($value) = @more_args;
        if ($arg_val > $value) {
            die "$arg_name should be less or equal than $value\n";
        }
    }
    elsif ($comparison eq "<>") {
        my ($min_value, $max_value) = @more_args;
        if (($arg_val <= $min_value) || ($arg_val >= $max_value)) {
            die "$arg_name should be between $min_value and $max_value (not including them)\n";
        }
    }          
    elsif ($comparison eq "<=>") {
        my ($min_value, $max_value) = @more_args;
        if (($arg_val < $min_value) || ($arg_val > $max_value)) {
            die "$arg_name should be between $min_value and $max_value\n";
        }
    }      
    elsif ($comparison eq "inS") {
        #S for string
        my @words = @more_args;

        my $found = 0;
        foreach my $word_for (@words) {
            if ($arg_val eq $word_for) {
                $found = 1;
                last;
            }
        }
        if ($found == 0) {
            die "$arg_name has to have one of the following values: " . join(", ", @words) . "\n";
        }
    }    
    elsif ($comparison eq ">V") {
        my ($other_arg_name, $value) = @more_args;
        if ($arg_val <= $value) {
            die "$arg_name should be greater than $other_arg_name\n";
        }
    }        
    elsif ($comparison eq ">=V") {
        my ($other_arg_name, $value) = @more_args;
        if ($arg_val < $value) {
            die "$arg_name should be greater or equal than $other_arg_name\n";
        }
    }       
    elsif ($comparison eq "<V") {
        my ($other_arg_name, $value) = @more_args;
        if ($arg_val >= $value) {
            die "$arg_name should be less than $other_arg_name\n";
        }
    }   
    elsif ($comparison eq "<=V") {
        my ($other_arg_name, $value) = @more_args;
        if ($arg_val > $value) {
            die "$arg_name should be less or equal than $other_arg_name\n";
        }
    }
    elsif ($comparison eq "<>V") {
        my ($other_arg_name_1, $min_value, $other_arg_name_2, $max_value) = @more_args;
        if (($arg_val <= $min_value) || ($arg_val >= $max_value)) {
            die "$arg_name should be between $other_arg_name_1 and $other_arg_name_2 (not including them)\n";
        }
    }          
    elsif ($comparison eq "<=>V") {
        my ($other_arg_name_1, $min_value, $other_arg_name_2, $max_value) = @more_args;
        if (($arg_val < $min_value) || ($arg_val > $max_value)) {
            die "$arg_name should be between $other_arg_name_1 and $other_arg_name_2\n";
        }
    }
}

1;