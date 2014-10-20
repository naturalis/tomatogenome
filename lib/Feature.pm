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

1;