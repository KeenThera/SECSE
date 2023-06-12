#!/usr/bin/perl
# -*- coding:utf-8 _*-
# @author: Lu Chong
# @file: filter_sdf_by_titles.pl
# @time: 2023/06/15/17:58

use strict;
use warnings;

# Check command line arguments
my ($input, $title_file, $output) = @ARGV;
die "Usage: perl filter_sdf_by_titles.pl <input> <title_file> <output>\n" unless defined $input;

# Open output file
open(my $out, ">$output") or die "Cannot open $output for writing: $!\n";

# Read titles
my %titles;
if ($title_file) {
    open(my $fh, "<", $title_file) or die "Cannot open $title_file for reading: $!\n";
    while (my $title = <$fh>) {
        chomp($title);
        $title =~ s/^\s+|\s+$//g;
        $titles{$title} = 1;
    }
    close($fh);
}

open(my $in, "<", $input) or die "Cannot open $input for reading: $!\n";

# Process input file
local $/ = '$$$$'; # Set the input record separator

my $numstructs = 0;
my @title_indices;
my $buffer;
my $index = 0;
my $first_structure = 1;

while ($buffer = <$in>) {
    $numstructs++;

    # Get the title of the current structure
    my $title = get_title($buffer);

    if (exists $titles{$title}) {
        $index++;
        if (!$title_indices[$index]) {
            $title_indices[$index] = $numstructs;
        } else {
            $title_indices[$index] .= ",$numstructs";
        }
        # if the first structure starts with a newline, then strip the newline
        if ($first_structure && $buffer =~ /^\n/) {
            $buffer =~ s/^\n//;
            $first_structure = 0;
        }
        print $out $buffer;
    }
}

# Add newline at the end of the output file
print $out "\n";

close($in);
close($out);

# Extract the first line as the CT title
sub get_title {
    my ($ct) = @_;
    $ct =~ s/^\s+//;
    my ($title) = $ct =~ /^(.+)$/m;
    return $title || '';
}
