#!/usr/bin/env perl

=pod

=head1 NAME

Summarize expression, either into median per biological replicate or into quartile per assay group.

=head1 SYNOPSIS

Assays -> assay groups:
cat my-genes.tsv.undecorated | gxa_summarize_expression.pl --configuration <path to configuration file> --aggregate-quartiles

Assays -> technical replicates:
cat my-transcripts.tsv.undecorated | gxa_summarize_expression.pl  --configuration <path to configuration file>

Pass input file through stdin. Output will be sent to stdout, and error messages - if any - to stderr.

=cut

use 5.10.0;

use Pod::Usage;
use Getopt::Long;
use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
use Atlas::FPKMmatrix;

my %args;

GetOptions(
    "c|configuration=s" => \$args{ "xml_config" },
    "a|aggregate-quartiles" => \$args{ "do_quartiles" }
) or pod2usage(-exitval => 255, -verbose => 1);

pod2usage(
    -message => "You must specify the path to configuration xml .\n",
    -exitval => 255,
    -verbose => 1
) unless $args{ "xml_config" } ;

Atlas::FPKMmatrix::run(
	parseAtlasConfig($args{"xml_config"}),
	$args{"do_quartiles"},
	STDIN,
	STDOUT
);

1;