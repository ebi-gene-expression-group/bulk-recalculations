#!/usr/bin/env perl
#
=pod

=head1 NAME

get_experiment_info.pl - retrieve certain information about an experiment.

=head1 SYNOPSIS

get_experiment_info.pl --experiment E-MTAB-1066 --xmlfile exp/E-MTAB-1066-configuration.xml --organism
get_experiment_info.pl --experiment E-MTAB-1066 --xmlfile exp/E-MTAB-1066-configuration.xml --arraydesign
get_experiment_info.pl --experiment E-MTAB-1066 --xmlfile exp/E-MTAB-1066-configuration.xml --rawdatafiles

=head1 DESCRIPTION

This script accepts an ArrayExpress experiment accession as an argument as well
as a list of options, to return the organism, array design accession(s), and/or
microarray raw data file names.

=head1 OPTIONS

=over 2

=item -e --experiment

Required. ArrayExpress accession of experiment.

=item -x --xmlfile

Required. Full path to Atlas XML configuration file.

=item -o --organism

Optional. Organism used in experiment. Only one is allowed.

=item -a --arraydesign

Optional. ArrayExpress accession of array designs used in a microarray experiment.

=item -r --rawdatafiles

Optional. Filenames of raw data files for a microarray experiment.

=item -m --magetabfiles

Optional. Filenames of IDF and SDRF files, comma separated.

=item -h --help

Optional. Display this help.

=back

=head1 AUTHOR

Expression Atlas team <arrayexpress-atlas@ebi.ac.uk>

=cut

use strict;
use warnings;
use 5.10.0;

use Atlas::Magetab4Atlas;
use EBI::FGPT::Config qw( $CONFIG );
use Atlas::Common qw(
    create_atlas_site_config
    make_ae_idf_path
    get_idfFile_path
);
use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
use File::Spec;
use Getopt::Long;
use Pod::Usage;
use Log::Log4perl;
use Data::Dumper;


# Logger config. This one only logs fatal errors, so it doesn't interfere with
# calling processes.
my $logger_config = q(
	log4perl.rootlogger			       = FATAL, SCREEN
    log4perl.appender.SCREEN             = Log::Log4perl::Appender::Screen
	log4perl.appender.SCREEN.stderr      = 0
	log4perl.appender.SCREEN.layout      = Log::Log4perl::Layout::PatternLayout
	log4perl.appender.SCREEN.layout.ConversionPattern = %-5p - %m%n
);

# Initialise logger.
Log::Log4perl::init(\$logger_config);
my $logger = Log::Log4perl::get_logger;

my $args = parse_args();

my $expAcc = $args->{ "experiment_accession" };

my $idfFile = get_idfFile_path( $expAcc );

my $magetab4atlas = Atlas::Magetab4Atlas->new( idf_filename => $idfFile, strict => !$args->{ "not_strict" } );

my $experimentType = $magetab4atlas->get_experiment_type;

# Get the assays.
my $assays = $magetab4atlas->get_assays;

# Remove assays that aren't in the XML config from the array.
$assays = remove_non_xml_assays( $assays, $args->{ "xml_filename" } );

if( $args->{ "organism" } ) { say_experiment_organism( $assays, $expAcc ); }
elsif( $args->{ "array_design" } ) {
    unless( $experimentType =~ /array/ ) {
        $logger->logdie( "$expAcc is not a microarray experiment, cannot find array designs." );
    }
    say_array_designs( $assays, $expAcc );
}
elsif( $args->{ "raw_data_files" } ) {
    unless( $experimentType =~ /array/ ) {
        $logger->logdie( "$expAcc is not a microarray experiment, cannot find raw microarray data files." );
    }
    say_raw_data_files( $assays, $expAcc );
}
elsif( $args->{ "magetabfiles" } ) {
  (my $sdrfFile = $idfFile) =~ s/\.idf\./\.sdrf\./g;
  say $idfFile.",".$sdrfFile;
}


sub parse_args {
    my %args;
    my $want_help;
    GetOptions(
        "h|help"            => \$want_help,
        "e|experiment=s"    => \$args{ "experiment_accession" },
        "o|organism"        => \$args{ "organism" },
        "a|arraydesign"     => \$args{ "array_design" },
        "r|rawdatafiles"    => \$args{ "raw_data_files" },
        "m|magetabfiles"    => \$args{ "magetabfiles" },
	"n|not_strict"      => \$not_strict, # added as negation to keep the default true as it was.
        "x|xmlfile=s"         => \$args{ "xml_filename" }
    );
    if( $want_help ) {
        pod2usage(
            -exitval    => 255,
            -output     => \*STDOUT,
            -verbose    => 1
        );
    }
    unless( $args{ "experiment_accession" } ) {
        pod2usage(
            -message    => "You must specify an experiment accession.\n",
            -exitval    => 255,
            -output     => \*STDOUT,
            -verbose    => 1
        );
    }
    unless( $args{ "xml_filename" } ) {
        pod2usage(
            -message    => "You must specify the full path to the Atlas XML configuration file.\n",
            -exitval    => 255,
            -output     => \*STDOUT,
            -verbose    => 1
        );
    }
    unless( $args{ "experiment_accession" } =~ /^E-\w{4}-\d+$/ ) {
        pod2usage(
            -message    => "\"" . $args{ "experiment_accession" } . "\" does not look like an ArrayExpress experiment accession.\n",
            -exitval    => 255,
            -output     => \*STDOUT,
            -verbose    => 1
        );
    }
    # Only do one attribute per script run.
    my $definedArgs;
    foreach my $arg ( keys %args ) {
        if( defined( $args{ $arg } ) ) { $definedArgs++; }
    }
    if( $definedArgs > 3 ) {
        pod2usage(
            -message    => "Please specify only one attribute type to retrieve.\n",
            -exitval    => 255,
            -output     => \*STDOUT,
            -verbose    => 1
        );
    }
    unless( $args{ "organism" } || $args{ "array_design" } || $args{ "raw_data_files" } || $args{ "magetabfiles" } ) {
        pod2usage(
            -message    => "You must specify an attribute to retrieve.\n",
            -exitval    => 255,
            -output     => \*STDOUT,
            -verbose    => 1
        );
    }

    return \%args;
}


sub say_experiment_organism {
    my ( $assays, $expAcc ) = @_;
    my $organisms = {};
    foreach my $assay ( @{ $assays } ) {
        my $characteristics = $assay->get_characteristics;
        my @organisms = keys %{ $characteristics->{ "organism" } };
        foreach my $organism ( @organisms ) {
            $organisms->{ $organism } = 1;
        }
    }

    if( ( keys %{ $organisms } ) > 1 ) {
        $logger->logdie( "ERROR - $expAcc has more than one organism " . Dumper($organisms));
    }
    else {
        my $organism = ( keys %{ $organisms } )[ 0 ];
        say $organism;
    }
}

sub say_array_designs {
    my ( $assays, $expAcc ) = @_;
    my $arrayDesigns = {};
    foreach my $assay ( @{ $assays } ) {
        my $arrayDesign = $assay->get_array_design;
        $arrayDesigns->{ $arrayDesign } = 1;
    }
    unless( keys %{ $arrayDesigns } ) {
        $logger->logdie( "ERROR - No array designs found for $expAcc." );
    }
    foreach my $arrayDesign ( keys %{ $arrayDesigns } ) {
        say $arrayDesign;
    }
}


sub say_raw_data_files {
    my ( $assays, $expAcc ) = @_;
    my $rawDataFiles = {};
    foreach my $assay ( @{ $assays } ) {
        if( $assay->has_array_data_file ) {
            my $rawDataFile = $assay->get_array_data_file;
            $rawDataFiles->{ $rawDataFile } = 1;
        }
    }

    unless( keys %{ $rawDataFiles } ) {
        $logger->logdie( "ERROR - No raw microarray data files found for $expAcc." );
    }

    foreach my $rawDataFile ( keys %{ $rawDataFiles } ) {
        say $rawDataFile;
    }
}


sub remove_non_xml_assays {
    my ( $assays, $xmlFilename ) = @_;
    my $experimentConfig = parseAtlasConfig( $xmlFilename );
    my $allAnalytics = $experimentConfig->get_atlas_analytics;
    my $xmlAssayNames = {};
    foreach my $analytics ( @{ $allAnalytics } ) {
        my $analyticsAssays = $analytics->get_assays;
        foreach my $analyticsAssay ( @{ $analyticsAssays } ) {
            my $assayName = $analyticsAssay->get_name;
            $xmlAssayNames->{ $assayName } = 1;
        }
    }

    my $assaysToKeep = [];
    foreach my $assay ( @{ $assays } ) {
        if( $xmlAssayNames->{ $assay->get_name } ) {

            push @{ $assaysToKeep }, $assay;
        }
    }

    return $assaysToKeep;
}
