package Locus;
use Moose;
use Feature;
use Math::Round;
use Data::Dumper;
use List::Util 'sum';
use List::MoreUtils 'firstidx';
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::Util::Logger ':levels';

# constructor args
has [ qw[locus source reference accessions] ] => ( is => 'ro', isa => 'Str' );
has [ qw[bp_collapse margin] ] => ( is => 'ro', isa => 'Int' );
has 'log' => ( is => 'ro', isa => 'Str', default => 'Bio::Phylo::Util::Logger' );
has 'svg' => ( is => 'ro', isa => 'SVG' );

# read from GFF3
has [ qw[strand] ] => ( is => 'rw', isa => 'Str' );
has [ qw[features] ] => ( is => 'rw', isa => 'ArrayRef' );

# _recalculated
has [ qw[start end collapse height locus_start locus_end locus_height] ] => ( is => 'rw', isa => 'Int' );
has 'bp_width' => ( is => 'rw', isa => 'Num' );

sub read_features {
	my ( $self, $gff3 ) = @_;
	my $locus  = $self->locus;
	my $source = $self->source;
	my @features;
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
			my $ccount = int( $length / $self->bp_collapse );	
			$self->strand($strand);
		
			# create object
			push @features, $type->new(
				'start'    => $start,
				'end'      => $end,
				'strand'   => $strand,
				'svg'      => $self->svg,
				'log'      => $self->log,
				'collapse' => $ccount,
			);
			$self->log->info("read $type at $start..$end on $strand strand, $length with $ccount collapses");
		}
	}
	
	# need to adjust exons that flank UTRs so they don't overlap
	for my $i ( 0 .. $#features ) {
		if ( $features[$i]->isa('UTR') ) {
			if ( $features[$i-1]->start == $features[$i]->start ) {
				$features[$i-1]->start($features[$i]->end);
			}
			if ( $features[$i-1]->end == $features[$i]->end ) {
				$features[$i-1]->end($features[$i]->start);
			}
		}
	}
	
	$self->features(\@features);
	return @features;
}

sub read_breakpoints {
	my ( $self, $nexus ) = @_;
	my @feat = @{ $self->features };
	my $ungapped = [ sort { $a <=> $b } map { $_->start } @feat ]->[0] - $self->margin;
	my @offsets;     # translates local aln coordinates to genomic coordinates
	my $ref_idx;     # index of reference taxon
	my $refseq;      # aligned sequence of reference taxon
	my $seq_counter = 0; # counts seen sequences	
	my ( $in_taxa, $in_matrix ) = ( 0, 0 ); # flags for state machine	
	
	# we're not using the Bio::Phylo NEXUS parser because 
	# HYPHY uses the obscure 'nolabels' directive
	$self->log->info("going to read breakpoints from $nexus");
	open my $fh, '<', $nexus or die $!;
	LINE: while(<$fh>) {
		chomp;
		
		# switch flags of state machine
		$in_taxa++   if /TAXLABELS/; # start of taxa
		$in_matrix++ if /MATRIX/;    # start of matrix
		$in_taxa--   if $in_taxa   and /END;/; # end of taxa
		$in_matrix-- if $in_matrix and /END;/; # end of matrix
		
		# find index of reference taxon
		if ( $in_taxa ) {
			next LINE if /TAXLABELS/;			
			my $ref = $self->reference;
			$ref_idx = firstidx { /^'$ref'$/ } grep { /^\S+$/ } split /\s+/, $_;	
			$self->log->info("reference taxon is at index $ref_idx");			
		}
		
		# find reference sequence, compute coordinate translation
		if ( $in_matrix ) {
			next LINE if /MATRIX/;
			if ( $seq_counter == $ref_idx ) {
				s/[ ;]//g; # strip whitespace and (possibly) trailing semicolon
				$refseq = $_;
				$self->log->info("refseq is ".length($refseq)." nt long");
				$self->log->info("genomic coordinates start at $ungapped");				
				
				# translate coordinates from aligned to ungapped genomic
				for my $i ( 0 .. ( length($refseq) - 1 ) ) {
					$ungapped += ( substr($refseq,$i,1) ne '-' );
					$offsets[$i] = $ungapped;
				}				
			}
			$seq_counter++;
		}
		
		# each charset defines a segment between breakpoints
		if ( /CHARSET span_\d+ = (\d+)-(\d+);/ ) {
			my ( $start, $end ) = ( $1, $2 );
			push @feat, breakpoint->new(
				'start' => $offsets[$start-1], # from aln pos to genomic coordinate
				'end'   => $offsets[$end-1],
				'svg'   => $self->svg,
				'log'   => $self->log,
			);
			$self->log->info("breakpoint at site $end (genomic: $offsets[$end])");
		}			
	}
	
	# XXX we add breakpoints as features
	$self->features( \@feat );
}

sub _read_trees {
	my ( $self, $nexus ) = @_;
	
	# reduce verbosity from Bio::Phylo
	$self->log->VERBOSE(
		'-level' => WARN,
		'-class' => [
			'Bio::Phylo',
			'Bio::Phylo::Util::OptionalInterface',
			'Bio::Phylo::Forest::TreeRole',
		],
	);
	
	# read newick strings from nexus file
	my @newick;
	{
		$self->log->info("going to read newick strings from $nexus");
		open my $fh, '<', $nexus or die $!;
		while(<$fh>) {
			chomp;
			if ( /^\s*TREE part_\d+ = (.+)$/ ) {
				my $tree = $1;
				push @newick, $tree;
			}		
		}
	}
	
	# parse to tree objects
	$self->log->info("going to parse ".scalar(@newick)." newick string");
	return parse(
		'-format' => 'newick',
		'-string' => join("\n",@newick),
	);	
}

sub read_crossovers {
	my ( $self, $nexus, @ingroup ) = @_;
	my $forest = $self->_read_trees($nexus);	
	my @feat = @{ $self->features };
	my ($start) = sort { $a <=> $b } map { $_->start } @feat;
	my ($end)   = sort { $b <=> $a } map { $_->end } @feat;
	my @bp = grep { $_->isa('breakpoint') } @feat;
	$self->log->info("going to superimpose crossovers on ".scalar(@bp)." breakpoints");
	
	# iterate over trees
	my $index = 0;
	my ( @segments, %outgroup );
	for my $tree ( @{ $forest->get_entities } ) {
		my %in = map { $_ => 1 } @ingroup;
		my @tips = @{ $tree->get_terminals };
		my %nearest;
		
		# iterate over tips in the ingroup
		SOURCE: for my $i ( 0 .. $#tips ) {
			my $source = $tips[$i]->get_name;
			next SOURCE unless $in{$source};
			my %dist; # distances to outgroup tips
			
			# iterate over tips in the outgroup
			TARGET: for my $j ( 0 .. $#tips ) {
				my $target = $tips[$j]->get_name;
				next TARGET if $in{$target};
				
				# calculate distance
				$dist{$target} = $tips[$i]->calc_patristic_distance($tips[$j]);				
			}
			($nearest{$source}) = sort { $dist{$a} <=> $dist{$b} } keys %dist;
			$outgroup{$nearest{$source}}++;
		}
		
		# create segment
		push @segments, segment->new( 
			'nearest' => \%nearest,
			'svg'     => $self->svg,
			'start'   => $bp[$index]->start,
			'end'     => $bp[$index]->end,
		);
		$index++;
		$self->log->debug("segment $index: ".Dumper(\%nearest));
	}	
	
	# instantiate the crossover object
	push @feat, crossover->new( 
		'segments' => \@segments,
		'ingroup'  => \@ingroup,
		'outgroup' => [ keys %outgroup ],
		'svg'      => $self->svg,	
		'start'    => $start,
		'end'      => $end,	
	);
	$self->features(\@feat);
}

sub _make_mapping {
	my $self = shift;
	my $acc = $self->accessions;
	my %map;
	open my $fh, '<', $acc or die $!;
	while(<$fh>) {
		chomp;
		if ( /^(.+?)=(.+)$/ ) {
			my ( $key, $value ) = ( $1, $2 );
			$map{$key} = $value;
		}	
	}	
	return \%map;
}

sub draw { 
	my $self = shift;
	$self->_recalculate;
	$self->_add_scale;	
	$self->_add_arrow;
	$self->_add_crossover;
	my $map = $self->_make_mapping;
	my @feat = @{ $self->features };
	$_->mapping($map) for @feat;
	
	# breakpoints are drawn first with other things on top of them
	map { $_->draw } grep { $_->isa('breakpoint') } @feat;
	map { $_->draw } grep { ! $_->isa('UTR') } @feat;
	map { $_->draw } grep { $_->isa('UTR') } @feat;
}

sub _recalculate {
	my $self   = shift;
	my @feat   = @{ $self->features };
	my $width  = $self->svg->{'-document'}->{'width'};
	my $height = $self->svg->{'-document'}->{'height'};

	# overall start in bp	
	$self->start( [ sort { $a <=> $b } map { $_->start } @feat ]->[0] ); 
	
	# overall end in bp
	$self->end( [ sort { $b <=> $a } map { $_->end } @feat ]->[0] ); 
	
	# locus start in bp
	$self->locus_start( [ sort { $a <=> $b } map { $_->start } grep { $_->isa('LocusFeature') } @feat ]->[0] ); 
	
	# locus end in bp
	$self->locus_end( [ sort { $b <=> $a } map { $_->end } grep { $_->isa('LocusFeature') } @feat ]->[0] ); 
	
	# tallest feature in pixels
	$self->height( [ sort { $b <=> $a } map { $_->height } @feat ]->[0] ); 
	
	# tallest locus feature in pixels
	$self->locus_height( [ sort { $b <=> $a } map { $_->height } grep { $_->isa('LocusFeature') } @feat ]->[0] ); 	
	
	# count of all collapsed regions
	$self->collapse( sum( map { $_->collapse } @feat ) ); 
	
	# width of pixel in bp
	$self->bp_width( ( $self->end - $self->start - ( $self->collapse * $self->bp_collapse ) ) / $width );
	
	# reporting
	{
		my $start    = $self->start;
		my $end      = $self->end;
		my $ccount   = $self->collapse;
		my $bp_width = $self->bp_width;
		$self->log->info("locus $start..$end has $ccount collapsed regions, 1 pixel = ${bp_width} bp");
	}

	# apply coordinates
	my $offset = 0; # running offset due to collapsed features
	for my $f ( sort { $a->start <=> $b->start } @feat ) {

		# compute starting coordinates	
		$f->x1( round( ( $f->start - ( $self->start + $offset ) ) / $self->bp_width ) );
		$f->y1( round( $height / 2 - $f->height / 2 ) );
	
		# increment offset
		$offset += $f->collapse * $self->bp_collapse;
	
		# compute ending coordinates	
		$f->x2( round( ( $f->end - ( $self->start + $offset ) ) / $self->bp_width ) );	
		$f->y2( round( $height / 2 - $f->height / 2 ) );		
	}
}

sub _add_crossover {
	my $self = shift;
	my $width       = $self->svg->{'-document'}->{'width'};	
	my $height      = $self->svg->{'-document'}->{'height'};	
	my $fheight     = $self->locus_height;
	my $bp_width    = $self->bp_width;
	my $bp_collapse = $self->bp_collapse;
	my ($crossover) = grep { $_->isa('crossover') } @{ $self->features };
	
	# compute and assign height range
	my ($arrow) = grep { $_->isa('arrow') } @{ $self->features };	
	my $y1 = round( $height / 2 + $fheight / 2 + $arrow->height * 2 + $crossover->spacing );
	my $y2 = round( $height - $crossover->spacing - 120 );
	$crossover->y1( $y1 );
	$crossover->y2( $y2 );
	
	# apply coordinates
	my $start = $self->start;	
	my $offset = 0; # running offset due to collapsed features
	for my $f ( sort { $a->start <=> $b->start } @{ $crossover->segments } ) {

		# compute starting coordinates	
		$f->x1( round( ( $f->start - ( $start + $offset ) ) / $bp_width ) );
	
		# increment offset
		$offset += $f->collapse * $bp_collapse;
	
		# compute ending coordinates	
		$f->x2( round( ( $f->end - ( $start + $offset ) ) / $bp_width ) );	
	}
}

sub _add_arrow {
	my $self    = shift;
	my $width   = $self->svg->{'-document'}->{'width'};	
	my $height  = $self->svg->{'-document'}->{'height'};	
	my $fheight = $self->locus_height;
	my @feat = @{ $self->features };

	# create arrow
	my $arrow = arrow->new(
		'start'  => $self->locus_start,
		'end'    => $self->locus_end,
		'x1'     => round( ( $self->locus_start - $self->start ) / $self->bp_width ),
		'x2'     => round( ( $self->locus_end - $self->start ) / $self->bp_width ),
		'strand' => $self->strand,
		'svg'    => $self->svg,
	);
	$arrow->y1( round( $height / 2 + $fheight / 2 + $arrow->height ) );
	$arrow->y2( round( $height / 2 + $fheight / 2 + $arrow->height ) );
	push @feat, $arrow;
	$self->log->info("added arrow glyph");
	$self->features(\@feat);
}

sub _add_scale {
	my $self    = shift;
	my $width   = $self->svg->{'-document'}->{'width'};	
	my $height  = $self->svg->{'-document'}->{'height'};
	my @feat = @{ $self->features };
	
	# create scale bar
	my $scale = scale->new(
		'start'    => $self->start,
		'end'      => $self->end,
		'x1'       => 0,
		'x2'       => $width,
		'svg'      => $self->svg,		
		'bp_width' => $self->bp_width,
	);
	$scale->y1( round( $height / 2 - $scale->height / 2 ) );
	$scale->y2( round( $height / 2 + $scale->height / 2 ) );
	unshift @feat, $scale;
	$self->log->info("added scale bar glyph");
	$self->features(\@feat);	
}

1;