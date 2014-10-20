package Locus;
use Moose;
use Feature;
use Math::Round;
use List::Util 'sum';

# constructor args
has [ qw[locus source] ] => ( is => 'ro', isa => 'Str' );
has 'log' => ( is => 'ro', isa => 'Str', default => 'Bio::Phylo::Util::Logger' );
has 'svg' => ( is => 'ro', isa => 'SVG' );
has 'bp_collapse' => ( is => 'ro', isa => 'Int' );

# read from GFF3
has [ qw[strand] ] => ( is => 'rw', isa => 'Str' );
has [ qw[features] ] => ( is => 'rw', isa => 'ArrayRef' );

# _recalculated
has [ qw[start end collapse height] ] => ( is => 'rw', isa => 'Int' );
has bp_width => ( is => 'rw', isa => 'Num' );

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
	$self->features(\@features);
	$self->_recalculate;
	$self->_add_arrow;
	return @features;
}

sub draw { map { $_->draw } @{ shift->features } }

sub _recalculate {
	my $self   = shift;
	my @feat   = @{ $self->features };
	my $width  = $self->svg->{'-document'}->{'width'};
	my $height = $self->svg->{'-document'}->{'height'};

	# overall start in bp	
	$self->start( [ sort { $a <=> $b } map { $_->start } @feat ]->[0] ); 
	
	# overall end in bp
	$self->end( [ sort { $b <=> $a } map { $_->end } @feat ]->[0] ); 
	
	# tallest feature in pixels
	$self->height( [ sort { $b <=> $a } map { $_->height } @feat ]->[0] ); 
	
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

	# draw features
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

sub _add_arrow {
	my $self    = shift;
	my $width   = $self->svg->{'-document'}->{'width'};	
	my $height  = $self->svg->{'-document'}->{'height'};	
	my $fheight = $self->height;
	my @feat = @{ $self->features };

	# draw arrow
	my $arrow = arrow->new(
		'start'  => $self->start,
		'end'    => $self->end,
		'x1'     => 0,
		'x2'     => $width,
		'strand' => $self->strand,
		'svg'    => $self->svg,
	);
	$arrow->y1( round( $height / 2 + $fheight / 2 + $arrow->height ) );
	$arrow->y2( round( $height / 2 + $fheight / 2 + $arrow->height ) );
	push @feat, $arrow;
	$self->log->info("added arrow glyph");
	$self->features(\@feat);
}

1;