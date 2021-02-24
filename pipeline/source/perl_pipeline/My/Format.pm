package My::Format;

use strict;
use warnings;

use Exporter qw(import);

=pod

=head1 NAME
    Format Module

=head1
DEPENDENCES
    Needs Perl

=head1
DESCRIPTION
    Module with Format related functions for Perl

=head1
AUTHORS
    Format.pm by Alejandro Ricci

=cut

our @EXPORT_OK = qw(cleanLineBorders
                    escapeText
                    prefixText
                    indentText
                    wrapTextBase
                    wrapText
                    checkType
                    checkTypeOutput
                    alignText
                    alignTextByMany
                    );

sub cleanLineBorders {
    #Arguments => line
    my ($line) = @_;

    $line =~ s/^\s+|\s+$//g; #remove trailing spaces
    $line =~ s/^\t+|\t+$//g; #remove trailing tabs
    $line =~ s/^\n+|\n+$//g; #remove trailing new lines

    return $line;
}

sub escapeText {
    #Arguments => text
    my ($text) = @_;

    $text =~ s/\\/\\\\/g;
    $text =~ s/\$/\\\$/g;
    $text =~ s/\"/\\\"/g;

    return $text;
}

sub prefixText {
    #Arguments => text prefix
    my ($text, $prefix) = @_;

    my $ends_with_newline = ($text =~ m/\n$/g);

    my @lines = split(/\n/, $text);

    foreach my $lineFor (@lines) {
        $lineFor = $prefix . $lineFor;
    }

    my $output = join("\n", @lines);
    $output .= "\n" if $ends_with_newline;

    return $output;
}
sub indentText {
    #Arguments => text number_of_spaces
    my ($text, $number_of_spaces) = @_;

    my $indent = " " x $number_of_spaces;

    return prefixText($text, $indent);
}

sub wrapTextBase {
    #Arguments => text width first_line_prefix other_lines_prefix character_to_jump_after
    my ($text, $width, $first_line_prefix, $other_lines_prefix, $character_to_jump_after) = @_;

    my $first_line_wrap_width = $width - length($first_line_prefix);
    my $other_lines_wrap_width = $width - length($other_lines_prefix);

    my @lines = split(/\n/, $text);

    my $output = "";
    for my $lineFor (@lines) {       
        if (length($lineFor) > $first_line_wrap_width) {
            #First line
            my $wrapped_line = substr($lineFor, 0, $first_line_wrap_width);
            my $last_jump_char_index = rindex($wrapped_line, $character_to_jump_after);

            if ($last_jump_char_index != -1) {
                $last_jump_char_index += 1; #to let the char stay in the line that had it
            } else {
                $last_jump_char_index = $first_line_wrap_width; #there was no such char, just cut the words
            }

            $output .= $first_line_prefix . substr($lineFor, 0, $last_jump_char_index) . "\n";
            $lineFor = substr($lineFor, $last_jump_char_index); #remove the part I already added to the output
            $lineFor =~ s/^\s+//g; #remove spaces at the beginning

            #Rest of the lines
            while (length($lineFor) > $other_lines_wrap_width) {
                $wrapped_line = substr($lineFor, 0, $other_lines_wrap_width);
                $last_jump_char_index = rindex($wrapped_line, $character_to_jump_after);

                if ($last_jump_char_index != -1) {
                    $last_jump_char_index += 1; #to let the space stay in the line that had it
                } else {
                    $last_jump_char_index = $other_lines_wrap_width; #there was no space in those chars, just cut the words
                }

                $output .= $other_lines_prefix . substr($lineFor, 0, $last_jump_char_index) . "\n";
                $lineFor = substr($lineFor, $last_jump_char_index); #remove the part I already added to the output
                $lineFor =~ s/^\s+//g; #remove spaces at the beginning
            }            

            $output .= $other_lines_prefix . $lineFor . "\n";
        } else {
            $output .= $first_line_prefix . $lineFor . "\n";
        }        
    }

    return $output;
}
sub wrapText {
    #Arguments => text width prefix
    my ($text, $width, $prefix) =  @_;

    return wrapTextBase($text, $width, $prefix, $prefix, " ");
}

sub checkType {
    #Arguments => text type
    my ($text, $type) =  @_;

    #Types are i (int), f (float), bi (boolean integer), c (char), s (string)
    my $is_valid = 0;
    if ($type eq "i") {
        $is_valid = 1 if ($text =~ m/^-?[0-9]+$/g); #Apparently I don't have to escape the $ if it's the last char
    } elsif ($type eq "f") {
        $is_valid = 1 if ($text =~ m/^-?[0-9.,]+$/g);
    } elsif ($type eq "bi") {
        $is_valid = 1 if ($text =~ m/^[01]$/g);
    } elsif ($type eq "c") {
        $is_valid = 1 if ($text =~ m/^.$/g);
        $is_valid = 1 if ($text eq ""); #also true for empty strings
    } elsif ($type eq "s") {
        $is_valid = 1; #there is no real control to do if it's a string
    }

    return $is_valid;
}
sub checkTypeOutput {
    #It's not exactly Format, but it's related with the function above
    my %error_message = (
        "i"  => "Integer",
        "f"  => "Float",
        "bi" => "Boolean (0 or 1)",
        "c"  => "Char",
        "s"  => "String",
    );

    return %error_message;
}

sub alignText {
    #Arguments => text string_to_align ignore_string_inside_quotes
    my ($text, $string_to_align, $ignore_string_inside_quotes) = @_;

    my $ends_with_newline = ($text =~ m/\n$/g);

    my @lines = split(/\n/, $text);
    
    my @ended;
    my @current_index;    
    for (my $i = 0; $i < scalar(@lines); $i++) {
        $ended[$i] = 0;
        $current_index[$i] = 0;
    }

    my $finished = 0;
    while (!$finished) {
        my $amount_to_align = 0;
        my $max_index = 0;
        
        for (my $i = 0; $i < scalar(@lines); $i++) {
            unless ($ended[$i]) {
                my $index_aux;

                if ($ignore_string_inside_quotes == 1) {
                    my $current_index_aux = $current_index[$i];
                    my $found = 0;
                    do {
                        $index_aux = index($lines[$i], $string_to_align, $current_index_aux);

                        if ($index_aux != -1) { 
                            my $left_side = substr($lines[$i], 0, $index_aux);
                            my $quotes_left = () = $left_side =~ /"/gi;
                            my $right_side = substr($lines[$i], $index_aux);
                            my $quotes_right = () = $right_side =~ /"/gi;

                            if (($quotes_right % 2) == 0) {
                                $found = 1;
                            } else {
                                $current_index_aux = $index_aux + 1;
                            }                            
                        }                       
                    } while ((!$found) && ($index_aux != -1));
                } else {
                    $index_aux = index($lines[$i], $string_to_align, $current_index[$i]);
                }                

                if ($index_aux != -1) {
                    $current_index[$i] = $index_aux;
                    $amount_to_align += 1;
                    $max_index = $current_index[$i] if $current_index[$i] > $max_index;
                } else {
                    $ended[$i] = 1;
                }
            }
        }

        if ($amount_to_align >= 2) {
            for (my $i = 0; $i < scalar(@lines); $i++) {
                unless ($ended[$i]) {
                    my $index_difference = $max_index - $current_index[$i];
                    my $spaces_to_add = " " x ($max_index - $current_index[$i]);
                    $lines[$i] = substr($lines[$i], 0, $current_index[$i]) . (" " x $index_difference) . substr($lines[$i], $current_index[$i]);
                    $current_index[$i] = $max_index + 1;
                }
            }
        } else {
            $finished = 1;
        }
    }

    my $output = join("\n", @lines);
    $output .= "\n" if $ends_with_newline;

    return $output;
}
sub alignTextByMany {
    #Arguments => text ignore_string_inside_quotes string_to_align string_to_align string_to_align string_to_align ...
    my ($text, $ignore_string_inside_quotes, @list_of_string_to_align) = @_;

    foreach my $string_to_align (@list_of_string_to_align) {
        $text = alignText($text, $string_to_align, $ignore_string_inside_quotes);
    }

    return $text;
}

1;