#!/usr/bin/env perl

use strict;
use warnings;
use YAML qw(LoadFile);

my ($input_sp, $YAML_file) = @ARGV;

# get total arg passed to this script
my $total = $#ARGV + 1;
print "Total args passed to the script: $total\n";

if ($total != 2) {
  die "The total arg passed is not equal to 2 (<organism> <yaml_file>)\n";
}

print "Searching '$input_sp' in YAML file '$YAML_file'\n";

#Read YAML
my $data = LoadFile($YAML_file);

#Search
my $found = 0;
my $output_sp = ' ';

for my $k1 (keys %{$data}) {
    for my $k2 (keys %{$data->{$k1}}) {
        # transform arrays into strings
        my $string1 = join('-' ,$k2);
        my $string2 = join( '-', @{$data->{$k1}{$k2}} );
        print "Key:  '$string1' -> '$string2' \n";

        if ($string1 eq $input_sp){
          print "MATCH ! \n";
          $found +=1;
          $output_sp = $string2 ;
        }
    }
}

print "Number of occurences: $found\n";

if ($found eq 0){
  print "$input_sp\n";
} else {
  print "$output_sp\n";
}


