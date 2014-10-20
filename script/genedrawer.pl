#!/usr/bin/perl
package Feature;
use Moose;

has [ qw[x1 y1 x2 y2 collapse] ] => ( is => 'rw', isa => 'Int' );
has [ qw[start end]  ] => ( is => 'ro', isa => 'Int' );
has [ qw[tag strand] ] => ( is => 'ro', isa => 'Str' );
has 'svg'    => ( is => 'ro', isa => 'SVG' );
has 'log'    => ( is => 'ro', isa => 'Str', default => 'Bio::Phylo::Util::Logger' );
has 'stroke' => ( is => 'rw', isa => 'Str', default => 'black' ); 
has 'fill'   => ( is => 'rw', isa => 'Str', default => 'white' ); 
has 'height' => ( is => 'rw', isa => 'Int', default => 50 );

sub draw {
	my $self = shift;
	$self->svg->tag( 
		$self->tag, 
		$self->attributes,
		'class' => ref($self),
	);
}

##########################################################################################
package exon;
use Moose;
use Math::Round;
extends 'Feature';

has '+tag'  => ( is => 'ro', isa => 'Str', lazy => 1, default => 'rect' );

sub attributes {
	my $self = shift;	
	return (
		'x'      => $self->x1,
		'y'      => $self->y1,
		'width'  => ( $self->x2 - $self->x1 ),
		'height' => $self->height,
		'style'  => {
			'stroke' => $self->stroke,
			'fill'   => $self->fill,		
		}
	);
}

##########################################################################################
package intron;
use Moose;
use Math::Round;
extends 'Feature';

has '+height' => ( is => 'rw', isa => 'Int', lazy => 1, default => 0 );
has '+tag' => ( is => 'ro', isa => 'Str', lazy => 1, default => 'line' );

sub attributes {
	my $self = shift;
	return (
		'x1' => $self->x1,
		'y1' => $self->y1,
		'x2' => $self->x2,
		'y2' => $self->y2,
		'style' => {
			'stroke' => $self->stroke,
			'fill'   => $self->fill,
		}
	);
}

##########################################################################################
package UTR;
use Moose;
extends 'exon';
has '+fill' => ( is => 'rw', isa => 'Str', lazy => 1, default => 'gray' ); 
##########################################################################################
package five_prime_UTR;
use Moose;
extends 'UTR';
##########################################################################################
package three_prime_UTR;
use Moose;
extends 'UTR';
##########################################################################################
package breakpoint;
use Moose;
extends 'Feature';
##########################################################################################
package arrow;
use Moose;
extends 'Feature';

has '+height' => ( is => 'rw', isa => 'Int', lazy => 1, default => 20 );
has '+tag'  => ( is => 'ro', isa => 'Str', lazy => 1, default => 'polygon' );

sub attributes {
	my $self = shift;
	my ( $x1, $x2, $y1, $y2, $h ) = ( $self->x1, $self->x2, $self->y1, $self->y2, $self->height );
	my ( @x, @y );
	
	# left to right ==>
	if ( $self->strand eq '+' ) {
		push @x, $x1, $x1, $x2-$h, $x2-$h, $x2, $x2-$h, $x2-$h, $x1;
		push @y, $y1+$h/4, $y1-$h/4, $y2-$h/4, $y2-$h/2, $y2, $y2+$h/2, $y2+$h/4, $y1+$h/4;
	}
	
	# right to left
	else {
		push @x, $x1+$h, $x1+$h, $x1, $x1+$h, $x1+$h, $x2, $x2, $x1+$h;
		push @y, $y1+$h/4, $y1+$h/2, $y1, $y1-$h/2, $y1-$h/4, $y2-$h/4, $y2+$h/4, $y1+$h/4;
	}
	
	# make the path
	my $path = $self->svg->get_path(
		'-type' => $self->tag,
		'x'     => \@x,
		'y'     => \@y,
	);
	
	return ( 
		%$path,
		'style'  => {
			'stroke' => $self->stroke,
			'fill'   => $self->fill,
		}
	);
}

##########################################################################################
package main;
use strict;
use warnings;
use SVG;
use Math::Round;
use Getopt::Long;
use List::Util 'sum';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $width     = 800;
my $height    = 600;
my $collapse  = 1000;
my $verbosity = WARN;
my ( $locus, $gff3, $source, $nexus, $reference );
GetOptions(
	'width=i'     => \$width,
	'height=i'    => \$height,
	'collapse=i'  => \$collapse,
	'locus=s'     => \$locus,
	'gff3=s'      => \$gff3,
	'source=s'    => \$source,
	'verbose+'    => \$verbosity,
	'nexus=s'     => \$nexus,
	'reference=s' => \$reference,
);

# instantiate helper objects
my $svg = SVG->new( 
	'width'  => $width, 
	'height' => $height 
);
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity
);

# read features from GFF3
$log->info("going to read (sub)features for locus $locus from $gff3");
my @features;
my $strandedness;
{
	open my $fh, '<', $gff3 or die $!;
	LINE: while(<$fh>) {
		chomp;
		
		# skip over header lines and separators
		next LINE if /^#/;
		
		# parse tabular lines
		my @line = split /\t/, $_;
		
		# skip over untrusted lines
		next LINE if $line[1] ne $source;
		
		# skip over lines that aren't part of the locus, terminate when at the next locus
		next LINE if $line[8] !~ /ID=(?:exon|intron|three_prime_UTR|five_prime_UTR):$locus\./;
		last LINE if $line[8] !~ /ID=(?:exon|intron|three_prime_UTR|five_prime_UTR):$locus\./ and @features;
		
		# parse salient details out of line
		my ( $type, $start, $end, $strand ) = @line[2,3,4,6];
		my $length = $end - $start;
		my $ccount = int( $length / $collapse );	
		$strandedness = $strand;	
		
		# create object
		push @features, $type->new(
			'start'    => $start,
			'end'      => $end,
			'strand'   => $strand,
			'svg'      => $svg,
			'log'      => $log,
			'collapse' => $ccount,
		);
		$log->info("read $type at $start..$end on $strand strand, $length with $ccount collapses");
	}
}

# calculate how many pixels a base pos corresponds with
my $start   = [ sort { $a <=> $b } map { $_->start  } @features ]->[0]; # overall start
my $end     = [ sort { $b <=> $a } map { $_->end    } @features ]->[0]; # overall end
my $fheight = [ sort { $b <=> $a } map { $_->height } @features ]->[0]; # tallest feature
my $ccount  = sum( map { $_->collapse } @features ); # count of all collapsed regions
my $bp_width = ( $end - $start - ( $ccount * $collapse ) ) / $width;
$log->info("locus $start..$end has $ccount collapsed regions, 1 pixel = ${bp_width} bp");

# draw features
my $offset = 0; # running offset due to collapsed features
for my $f ( sort { $a->start <=> $b->start } @features ) {

	# compute starting coordinates	
	$f->x1( round( ( $f->start - ( $start + $offset ) ) / $bp_width ) );
	$f->y1( round( $height / 2 - $f->height / 2 ) );
	
	# increment offset
	$offset += $f->collapse * $collapse;
	
	# compute ending coordinates	
	$f->x2( round( ( $f->end - ( $start + $offset ) ) / $bp_width ) );	
	$f->y2( round( $height / 2 - $f->height / 2 ) );	
	
	# draw feature
	$f->draw;
}

# draw arrow
my $arrow = arrow->new(
	'start'  => $start,
	'end'    => $end,
	'x1'     => 0,
	'x2'     => $width,
	'strand' => $strandedness,
	'svg'    => $svg,
);
$arrow->y1( round( $height / 2 + $fheight / 2 + $arrow->height ) );
$arrow->y2( round( $height / 2 + $fheight / 2 + $arrow->height ) );
$arrow->draw;

# done!
print $svg->xmlify;