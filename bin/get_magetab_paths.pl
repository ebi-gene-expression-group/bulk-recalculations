#!/usr/bin/env perl

# Script to get IDF path and Array Express MirBase load directory
# -i flag will retrieve path to idf filename
# -d flag will retrieve path to Array Express load directory
# -m flag will retrieve path to mirbase directory

# usage : get_magetab_paths.pl -e $expTargetDir -i

use strict;
use warnings;
use 5.10.0;

use File::Spec;
use File::Basename;
use File::Copy;
use File::Path;
use Getopt::Long;
use Pod::Usage;
use Log::Log4perl;
use Atlas::Common qw( create_atlas_site_config );

my $logger_config = q(
  log4perl.rootlogger           = INFO, SCREEN
  log4perl.appender.SCREEN        = Log::Log4perl::Appender::Screen
  log4perl.appender.SCREEN.stderr     = 0
  log4j.PatternLayout.cspec.Q       = sub { return "[QC]" }
  log4perl.appender.SCREEN.layout     = Log::Log4perl::Layout::PatternLayout
  log4perl.appender.SCREEN.layout.ConversionPattern = %-5p - %Q %m%n
);

# Initialise logger.
Log::Log4perl::init( \$logger_config );
my $logger = Log::Log4perl::get_logger;

my $abs_path = dirname(File::Spec->rel2abs(__FILE__));

my $atlasSiteConfig = create_atlas_site_config;

sub parse_args
{
  my ( %args);
  GetOptions(
    "e|exp=s" => \$args{ "experiment_directory" },
    "i|idf" => \$args{ "idf_path" },
    "d|ae-directory" => \$args{ "arrayexpress_dir" },
    "m|mirbase-directory" => \$args{ "mirbase_dir" },
  );

  unless ( $args{experiment_directory} )
  {
  pod2usage(
    -message => 'You must provide an experiment directory.',
    -exitval => 255,
    -output  => \*STDOUT,
    -verbose => 0,
  );
  }

  unless ($args{idf_path} or $args{arrayexpress_dir} or $args{mirbase_dir})
  {
  pod2usage(
    -message => 'Provide flags for IDF path (-i) flag or Array Express load dir (-d) or mirbase dir (-m)',
    -exitval => 255,
    -output  => \*STDOUT,
    -verbose => 0,
  );
  }
 return ( \%args );
}

# Get our arguments
my $args = parse_args();

my $atlasExperimentDir = $args->{experiment_directory};

my $exptAccession = (split '\/', $atlasExperimentDir)[-1];

# Path to directory with ATLAS_PROD, ArrayExpress/Atlas load directories.
my $atlasProdDir = $ENV{ "ATLAS_PROD" };
my $aeLoadDir = $ENV{"AE2_BASE_DIR"};

unless( $atlasProdDir ) {
  $logger->logdie( "ATLAS_PROD environment variable is not defined, cannot continue." );
}
unless( $aeLoadDir ) {
  $logger->logdie( "AE2_BASE_DIR environment variable is not defined, cannot continue." );
}

# Get the pipeline (e.g. MEXP, MTAB, GEOD, ...) for this experiment.
(my $pipeline = $exptAccession) =~ s/E-(\w{4})-\d+/$1/;

## if any GEOD or ENAD studies, retrieve raw fles from AtlasProd and curated magetab files
## or else look in ArrayExpress
my $rawFilesDir;
if ( $pipeline eq "GEOD" ) {
  $rawFilesDir = File::Spec->catdir( $atlasProdDir, "GEO_import" );
}
elsif ( $pipeline eq "ENAD" ) {
  $rawFilesDir = File::Spec->catdir( $atlasProdDir, "ENA_import" );
}
else {
  $rawFilesDir = $aeLoadDir;
}
# Experiment load directory and IDF filename and mirBase directory.
my $loadDir = File::Spec->catdir( $rawFilesDir, $pipeline, $exptAccession );
my $idfFilename = File::Spec->catfile( $loadDir, "$exptAccession.idf.txt" );
my $miRBaseDirectory = File::Spec->catdir( $atlasProdDir, $atlasSiteConfig->get_mirbase_mappings_directory );

# Check for GEO datasets to see if they are new (idf exists in the current set)
# directory or not (old GEO, needs to be retrieved from AE).
if(! -e $idfFilename && $pipeline eq "GEOD") {
  $loadDir = File::Spec->catdir( $aeLoadDir, $pipeline, $exptAccession );
  $idfFilename = File::Spec->catfile( $loadDir, "$exptAccession.idf.txt" );
}

if ( $args->{idf_path} ) {
    say $idfFilename;
}

if ($args->{arrayexpress_dir}){
     say $loadDir;
}

if ($args->{mirbase_dir}){
     say $miRBaseDirectory;
}