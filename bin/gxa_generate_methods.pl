#!/usr/bin/env perl

use strict;
use warnings;

use File::Spec;
use Atlas::Common qw(
    create_atlas_site_config
);

$_ = shift for my ( $versionReferenceFile, $expAcc, $species, $template, $baselineMapper, $baselineQuantMethod, $deMapper, $deQuantMethod, $deDEMethod );

my $atlasProd = $ENV{'ATLAS_PROD'};
my $irapSingleLib = $ENV{'IRAP_SINGLE_LIB'};
my $irapProd = $ENV{'IRAP_DIR'};

if (!$irapSingleLib) { die "Set environment variable \$IRAP_SINGLE_LIB.\n" }
if (!$irapProd) { die "Set environment variable \$IRAP_DIR.\n" }
if (!$atlasProd) { die "Set environment variable \$ATLAS_PROD.\n" }
if (!$baselineMapper) { die "Set variable \$baselineMapper.\n" }
if (!$baselineQuantMethod) { die "Set variable \$baselineQuantMethod.\n" }
if (!$deMapper) { die "Set variable \$deMapper.\n" }
if (!$deQuantMethod) { die "Set variable \$deQuantMethod.\n" }
if (!$deDEMethod) { die "Set variable \$deDEMethod.\n" }


my $atlasSiteConfig = create_atlas_site_config;

# gxa.references.conf
my $gxaReferencesConf = File::Spec->catfile( $atlasProd, $atlasSiteConfig->get_genome_references_config );

# differential or baseline?
my $diff = 0;
if($template =~ /differential/) { $diff++; }

# A hash for the info that needs to go into the output.
# This will have a key for each tag in the template and the corresponding
# value it finds.
# Tags are:
# 	<EXP_ACC>: experiment accession
# 	<PIPELINE_VERSION>: `irap | grep '* IRAP' | awk '{print $3}'` (stderr TBC) irap version
#       <ENS_SECTION>: Ensembl Genomes site (if available; e.g. Plants, Fungi, ...)
# 	<REF_REL>: Ensembl release version
# 	<MAPPER>: mapper name e.g. tophat2
# 	<MAPPER_VERSION>: mapper version number
# 	<QUANT_METHOD>: quantification tool e.g. cufflinks1
# 	<QUANT_METHOD_VERSION>: quantification tool version number
# 	<GSEA_VERSION>: version of the piano package

# Optionally for Differential Atlas, also:
# 	<DE_PKG>: package used for differential expression
# 	<DE_PKG_VERSION>>: differential expression package version number
my $tagsHash = {};

$tagsHash->{"EXP_ACC"} = $expAcc;

# Get QUANT_METHOD and DE_PKG from environmental variables (set prior to calling this script by isl/gxa_preirap.conf
if($diff) {
  $tagsHash->{"QUANT_METHOD"} = $deQuantMethod;
  $tagsHash->{"DE_PKG"} = $deDEMethod;
} else {
  $tagsHash->{"QUANT_METHOD"} = $baselineQuantMethod;
}

# Get REF_NAME (Ensembl or WBPS) based on whether annotation source file lives in WBPS dir or not.
if( -e File::Spec->catfile(
        $atlasProd,
        "sw",
        "atlasinstall_prod",
        "atlasprod",
        "bioentity_annotations",
        "fetch",
        "annsrcs",
        "wbps",
        $species
    ) ) {
    $tagsHash->{ "REF_NAME" } = "Wormbase ParaSite";
}
elsif( -e File::Spec->catfile(
        $atlasProd,
        "sw",
        "atlasinstall_prod",
        "atlasprod",
        "bioentity_annotations",
        "fetch",
        "annsrcs",
        "ensembl",  
        $species
    ) ) {
    $tagsHash->{ "REF_NAME" } = "Ensembl";
}

unless( $tagsHash->{ "REF_NAME" } ) { die "No reference genome resource name found (e.g. Ensembl or Wormbase ParaSite)."; }

# Get REF_REL from $irapOrganismConfigFile
my $irapOrganismConfigFile="$irapSingleLib/configs/$species.conf";
open(CONF, "<", $irapOrganismConfigFile) or die("Can't open $irapOrganismConfigFile: $!\n");
while(defined(my $line = <CONF>)) {
	chomp $line;
	# Get Ensembl release version
    if( $tagsHash->{ "REF_NAME" } eq "Ensembl" ) {
        if($line =~ /^gtf_file=.*\.(\d+)\.gtf/) {
            $tagsHash->{"REF_REL"} = $1;
        }
    }
    elsif( $tagsHash->{ "REF_NAME" } eq "Wormbase ParaSite" ) {
        if( $line =~ /^gtf_file=.*WBPS(\d+).canonical_geneset\.gtf\.gz/ ) {
            $tagsHash->{ "REF_REL" } = $1;
        }
    }
}
close(CONF);

unless( $tagsHash->{ "REF_REL" } ) { die "No reference genome resource release number found."; }

# Check we have all required info from the experiment config
checkTags("EXP_ACC", "REF_REL", "QUANT_METHOD");
if($diff) { checkTags("DE_PKG"); }

############################################################################
# Check if the reference is an Ensembl Genomes species and get the relevant
# Ensembl Genomes site name if so (i.e. Plants, Fungi, Bacteria, Protists and
# Metazoa).
############################################################################
open(GXAREFS, $gxaReferencesConf) or die("Can't open $gxaReferencesConf: $!\n");
while(defined(my $line = <GXAREFS>)) {

	chomp $line;

	if($line =~ /^$species.*(ensembl.*)/) {

		my $ensURL = $1;

		# See if this is an Ensembl Genomes species
		if($ensURL =~ /ensemblgenomes/) {

			# Append genomes
			my $ensSect ='genomes';

			# Capitalise the first letter :)
			$ensSect =~ s/^(\w)/\U$1/;

            # Replace RELNO with REF_REL as a fix to try produce valid links
            my $releaseNumber = $tagsHash->{ "REF_REL" };
            $ensSect =~ s/RELNO/$releaseNumber/g;

			# Add it to %tagsHash with a space at the end.
			$tagsHash->{"ENS_SECTION"} = "$ensSect ";
		}
		# Otherwise this isn't an Ensembl Genomes species, so just add an empty
		# string to %tagsHash
		else {
			$tagsHash->{"ENS_SECTION"} = "";
		}
	}
}
close(GXAREFS);
# Check we got something for ENS_SECTION tag.
if( $tagsHash->{ "REF_NAME" } eq "Ensembl" ) {
    checkTags("ENS_SECTION");
}

my @irapOut = `cat $versionReferenceFile`;
my $gseaVersion = `grep '^piano_version=' $irapProd/aux/mk/irap_versions.mk | awk -F"=" '{print \$NF}'`;

# Get differential expression package version from Atlas installation.
my $deMethodVersion = _get_atlas_rnaseq_de_version( $deDEMethod, $atlasSiteConfig );

# Set what mappers to possibly expect
# Robert: we used to use star, but since tophat2 was fixed to deal with large genomes, we use tophat2 for all our species
# Old experiments could nevertheless have been mapped with the big genomes mapper
# big genomes mapper not supported
my @mappers;
if($diff) {
    @mappers = ($deMapper);
} else {
    @mappers = ($baselineMapper);
  }

# only want the first IRAP line as that's the one with the version number
my $gotIrap = 0;
foreach my $line (@irapOut) {
	if($line =~ /IRAP\t([\d\w.]+)/ && !$gotIrap) {
		$tagsHash->{"PIPELINE_VERSION"} = $1;
		$gotIrap++;
	}
	elsif($line =~ /^[\w\s,]+\t(\w+)\t(.*)\t/) {
		my ($pkg, $vers) = ($1, $2);
        foreach my $mapper (@mappers) {
            if($pkg =~ /$mapper/i) {
                $tagsHash->{"MAPPER_VERSION"} = $vers;
                $tagsHash->{"MAPPER"} = $mapper;
            }
        }
		if($pkg =~ /$tagsHash->{"QUANT_METHOD"}/i) { $tagsHash->{"QUANT_METHOD_VERSION"} = $vers; }
		# special part for htseq/htseq2
		else {
			(my $qm_noNum = $tagsHash->{"QUANT_METHOD"}) =~ s/\d+$//;
			if($pkg =~ /$qm_noNum/i) { $tagsHash->{"QUANT_METHOD_VERSION"} = $vers; }
		}
	}
}
$tagsHash->{"GSEA_VERSION"} = $gseaVersion;
$tagsHash->{"DE_PKG_VERSION"} = $deMethodVersion;
checkTags("PIPELINE_VERSION","MAPPER", "MAPPER_VERSION", "QUANT_METHOD_VERSION");
if($diff) { checkTags("DE_PKG_VERSION", "GSEA_VERSION"); }

############################################################################
# Read template and write output
############################################################################
open(TEMPLATE, "<", $template) or die("Can't open $template: $!\n");
while(defined(my $line = <TEMPLATE>)) {

	# If we come across a potential tag (enclosed in <>)
	if($line =~ /<[A-Z_]+>/) {

		# Split on "<"
		my @splitLine = split "<", $line;
		$line = "";
		foreach my $part (@splitLine) {
			# find the potential tag, at the start of the line
			if($part =~ /^([A-Z_]+)>/) {
				my $tag = $1;

				# see if there is a key which matches the potential tag
				if( grep $_ eq $tag, keys %{ $tagsHash } ) {

					# if so, replace it with the value from the hash
					$part =~ s/$tag>/$tagsHash->{$tag}/;
				}
                elsif( $tagsHash->{ "REF_NAME" } eq "Wormbase ParaSite" && $tag eq "ENS_SECTION" ) {

                    $part =~ s/$tag>//;
                }

			} elsif (length $line > 0) {
			    $line .= "<";
			}
			$line .= $part;
		}

	}
        print $line;
}
close(TEMPLATE);





############################################################################
# Subroutines
############################################################################

# checkTags
# 	- Takes a tag e.g. "MAPPER_VERSION" (or array of these) as argument.
# 	- If it doesn't exist in the hash, dies with message.
sub checkTags {
	my @tags = @_;
	foreach my $tag (@tags) {
		unless(exists($tagsHash->{$tag})) {
			die("Error: couldn't find value for tag \"$tag\"\n");
		}
	}
}


sub _get_atlas_rnaseq_de_version {

    my ( $dePkgName, $atlasSiteConfig ) = @_;

    # Get the path to the Atlas production directory.
    my $atlasProdDir = $ENV{ "ATLAS_PROD" };

    # Get the path to the Atlas R installation
    my $Rinstallation = $atlasSiteConfig->get_atlas_r_installation;

    # Create the path to the DE package's DESCRIPTION file, which contains the
    # version.
    my $dePackageDescFile = File::Spec->catfile(
        $atlasProdDir,
        $Rinstallation,
        "lib64", "R", "library",
        $dePkgName,
        "DESCRIPTION"
    );

   # Make sure this file exists and is readable.
    unless( -r $dePackageDescFile ) {
        die( "Error: cannot read file $dePackageDescFile. Please check it exists and is readable.\n" );
    }

    # Grep for the line containing the version.
    my $versionLine = `grep -i version $dePackageDescFile`;
    # Remove newline.
    chomp $versionLine;

    # Check the grepping worked.
    if( $? ) {
        die( "Error: problem grepping in $dePackageDescFile for $dePkgName version.\n" );
    }

    unless( $versionLine ) {
        die( "Error: did not find version in $dePackageDescFile\n" );
    }

    # If we're still here, we got a version line. Parse out the version.
    ( my $dePkgVersion = $versionLine ) =~ s/^\D*//;

    return $dePkgVersion;
}
