#!/usr/bin/env perl

# Script to round log2 fold-changes from differential expression analysis to
# the nearest 0.1.

use strict;
use warnings;

use Math::Round;

my $analyticsFile = shift;

# Read in file.
open my $in, "<", $analyticsFile or die "Cannot open $analyticsFile for reading: $!\n";
my @lines = <$in>;
close $in;

# Write out, replacing log2 fold-changes with rounded values.
open my $out, ">", "$analyticsFile.rounded" or die "Cannot open $analyticsFile.rounded for writing: $!\n";

# Flag for first line.
my $firstline = 1;
# Array to store indices of log2 fold-change columns.
my @indices = ();
my $reordered_indexes;

# Go through each line.
foreach my $line ( @lines ) {

	# Remove newline.
	chomp $line;

	# Split on tabs.
	my @splitLine = split "\t", $line;

	# If this is the first line, don't need to change it so write it back out,
	# but also get the indices of log2 fold-change columns.
	if( $firstline ) {
		@indices = ( grep { $splitLine[$_] =~ /log2foldchange/ } 0..$#splitLine );
		# Check if we need to reverse logfold and p-value
		$reordered_indexes = &index_order_for_pvalue_logfoldchange($line);
		my $r_line = join "\t", @splitLine[@$reordered_indexes];

		print $out "$r_line\n";
		# Unset flag.
		$firstline = 0;
	}
	# For all other lines, round all the log2 fold-changes to the nearest 0.1,
	# then write them back out.
	else {

		foreach my $index ( @indices ) {
			$splitLine[ $index ] = nearest( .1, $splitLine[ $index ] );
			if( $splitLine[ $index ] == -0 ) {
				$splitLine[ $index ] = 0;
			}
		}

		$line = join "\t", @splitLine[@$reordered_indexes];
		print $out "$line\n";
	}
}

close $out;

sub index_order_for_pvalue_logfoldchange {
	# The web application requires that p-value fields are always
	# before their matching log2foldchange field (for the same group)
	# we have seen proteomics differential experiments not respecting this
	my ( $firstline ) = @_;

	my @splitLine = split "\t", $firstline;
	my $starting_data_index = 2; # provided this is a decorated file unless...
	# the second round doesn't contain gene symbols (undecorated):
	$starting_data_index = 1 if( $splitLine[1] =~ /g\d+_g\d+\./ );

	# first part doesn't get reordered
	my @reordered_indexes = 0..($starting_data_index-1);
	# now for every pair (i,i+1) starting from the first data index
	# we need to check if the p-value column is first or not, and swap if needed.
	for (my $i = $starting_data_index; $i < scalar @splitLine; $i+=2) {
		# in proteomics, where we see the issue, column says <contrast>.foldChange
		# instead of log2foldchange
		if( $splitLine[$i] =~ /foldchange/i ) {
			# log2foldchange is first, we need to swap
			my $contrast_lf = $splitLine[$i] =~ s/\..*//r;
			my $contrast_pv = $splitLine[$i+1] =~ s/\..*//r;
			if($contrast_lf eq $contrast_pv) {
				# we have confirmed they are the same contrast
				push @reordered_indexes, ($i+1, $i);
				next;
			}
		}
		# if the first one is not log2foldchange or contrast check failed
		push @reordered_indexes, ($i, $i+1);
	}

	return \@reordered_indexes;
}