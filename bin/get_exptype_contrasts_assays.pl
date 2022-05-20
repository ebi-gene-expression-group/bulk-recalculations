#!/usr/bin/env perl

use strict;
use warnings;
use 5.10.0;

use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );
use YAML qw(DumpFile);

my $accession = shift;
my $xmlFilename = shift;
my $output_yaml = shift;

my $experimentConfig = parseAtlasConfig( $xmlFilename );
my $experimentType = $experimentConfig->get_atlas_experiment_type;

# print "experiment_type\t".$experimentType."\n";
my %data = ( 'experiment_type' => $experimentType );

if ($experimentType =~ /differential/) {
  my %contrasts;
  foreach my $analytics ( @{ $experimentConfig->get_atlas_analytics } ) {
      foreach my $contrast ( @{ $analytics->get_atlas_contrasts } ) {
          $contrasts{$contrast->get_contrast_id} = $contrast->get_contrast_name;
      }
  }
  $data{'contrasts'} = \%contrasts;
  # keys and values have the same order as per
  # http://perldoc.perl.org/functions/keys.html
  # my $ids = join "::", keys %contrasts;
  # my $values = join "&&", values %contrasts;

  # print "contrasts_ids\t".$ids."\n";
  # print "contrasts_labels\t".$values."\n";
}

my %assays;
foreach my $analytics ( @{ $experimentConfig->get_atlas_analytics } ) {
  foreach my $assay ( @{ $analytics->get_atlas_assay_groups  } ) {
    $assays{$assay->get_assay_group_id} = $assay->get_label;
  }
}

$data{'assays'} = \%assays;

# my $ids = join "::", keys %assays;
# my $values = join "&&", values %assays;

# print "assays_ids\t".$ids."\n";
# print "assays_labels\t".$values."\n";

print DumpFile($output_yaml, \%data);
