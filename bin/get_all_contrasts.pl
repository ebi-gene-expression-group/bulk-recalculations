#!/usr/bin/env perl

use strict;
use warnings;
use 5.10.0;

use Atlas::AtlasConfig::Reader qw( parseAtlasConfig );

my $xmlFilename = shift;
my $experimentConfig = parseAtlasConfig( $xmlFilename );
my @contrastIDs = ();

my %contrasts = ();
foreach my $analytics ( @{ $experimentConfig->get_atlas_analytics } ) {
    foreach my $contrast ( @{ $analytics->get_atlas_contrasts } ) {
        $contrasts{$contrast->get_contrast_id} = $contrast->get_contrast_name;
    }
}
print join " ", @contrastIDs;
