#!/usr/bin/perl
use strict;
use warnings;
use Locus;
use SVG;
use Math::Round;
use Getopt::Long;
use List::Util 'sum';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $width     = 800;
my $height    = 600;
my $collapse  = 1000;
my $margin    = 10000;
my $verbosity = WARN;
my ( $locus, $gff3, $source, $nexus, $reference, $accessions );
my @ingroup;
GetOptions(
	'width=i'      => \$width,
	'height=i'     => \$height,
	'collapse=i'   => \$collapse,
	'locus=s'      => \$locus,
	'gff3=s'       => \$gff3,
	'source=s'     => \$source,
	'verbose+'     => \$verbosity,
	'nexus=s'      => \$nexus,
	'reference=s'  => \$reference,
	'margin=i'     => \$margin,
	'ingroup=s'    => \@ingroup,
	'accessions=s' => \$accessions,
);

# instantiate helper objects
my $svg = SVG->new( 
	'width'       => $width, 
	'height'      => $height 
);
my $log = Bio::Phylo::Util::Logger->new(
	'-level'      => $verbosity
);
my $locusDrawer = Locus->new(
	'svg'         => $svg,
	'log'         => $log,
	'bp_collapse' => $collapse,
	'locus'       => $locus,
	'source'      => $source,
	'reference'   => $reference,
	'margin'      => $margin,
	'ingroup'     => \@ingroup,
	'accessions'  => $accessions,
);

# make the drawing
$locusDrawer->read_features($gff3);
$locusDrawer->read_breakpoints($nexus);
$locusDrawer->read_crossovers($nexus,@ingroup);
$locusDrawer->draw;

# wide characters otherwise
binmode(STDOUT, ":utf8");
print $svg->xmlify;