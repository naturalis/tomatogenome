package Feature;
use Moose;
use Math::Round;

has [ qw[x1 y1 x2 y2 start end] ] => ( is => 'rw', isa => 'Int' );
has [ qw[tag strand]  ] => ( is => 'ro', isa => 'Str' );
has 'svg'      => ( is => 'ro', isa => 'SVG' );
has 'log'      => ( is => 'ro', isa => 'Str', default => 'Bio::Phylo::Util::Logger' );
has 'stroke'   => ( is => 'rw', isa => 'Str', default => 'black' ); 
has 'fill'     => ( is => 'rw', isa => 'Str', default => 'white' ); 
has 'height'   => ( is => 'rw', isa => 'Int', default => 50 );
has 'collapse' => ( is => 'rw', isa => 'Int', default => 0 );
has 'rotate'   => ( is => 'rw', isa => 'Int', default => 0 );
has 'text'     => ( is => 'rw', isa => 'Str' );
has 'mapping'  => ( is => 'rw', isa => 'HashRef' );

sub draw {
	my $self = shift;
	$self->svg->tag( 
		$self->tag, 
		$self->attributes,
		'class' => ref($self),
	);
	if ( my $text = $self->text ) {
		my $x = 
		my %args = (
			'x'      => $self->x1 + 4,
			'y'      => round( ( $self->y1 + $self->y2 ) / 2 - 10 ),
			'-cdata' => $text,
			'class'  => ref($self),
			'style'  => { 'font' => '10px Verdana' },				
		);
		if ( my $angle = $self->rotate ) {
			$args{'transform'} = "rotate($angle,$args{x},$args{y})";
		}
		$self->svg->text(%args);
	}
}

##########################################################################################
package LocusFeature;
use Moose;
extends 'Feature';

##########################################################################################
package exon;
use Moose;
use Math::Round;
extends 'LocusFeature';

has '+tag'    => ( is => 'ro', isa => 'Str', lazy => 1, default => 'rect' );
has '+rotate' => ( is => 'rw', isa => 'Int', lazy => 1, default => -90 );

sub attributes {
	my $self = shift;	
	return (
		'x'      => $self->x1,
		'y'      => $self->y1,
		'width'  => ( $self->x2 - $self->x1 ) || 1,
		'height' => $self->height,
		'style'  => {
			'stroke' => $self->stroke,
			'fill'   => $self->fill,		
		}
	);
}

sub text {
	my $self  = shift;
	my $start = $self->start;
	my $end   = $self->end;
	$start =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
	$end   =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
	return $start . ' .. ' . $end;
}

##########################################################################################
package intron;
use Moose;
use Math::Round;
extends 'LocusFeature';

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
has '+fill' => ( is => 'rw', isa => 'Str', lazy => 1, default => 'silver' ); 

sub text {
	my $self  = shift;
	my $start = $self->start;
	my $end   = $self->end;
	$start =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
	$end   =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
	return $start . ' .. ' . $end;
}

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

has '+rotate' => ( is => 'rw', isa => 'Int', lazy => 1, default => -90 );
has '+tag'    => ( is => 'ro', isa => 'Str', lazy => 1, default => 'line' );
has '+height' => ( is => 'rw', isa => 'Int', lazy => 1, default => 400 );

sub text {
	my $start = shift->start;
	$start =~ s/(\d{1,3}?)(?=(\d{3})+$)/$1,/g;
	return $start;
}

sub attributes {
	my $self = shift;
	my $y2 = $self->y2 + $self->height / 2;
	return (
		'x1' => $self->x2,
		'y1' => $self->y1,
		'x2' => $self->x2,
		'y2' => $y2,
		'style' => {
			'stroke' => $self->stroke,
			'fill'   => $self->fill,
			'stroke-width' => 0.1,			
		}
	);
}

##########################################################################################
package arrow;
use Moose;
extends 'Feature';

has '+height' => ( is => 'rw', isa => 'Int', lazy => 1, default => 10 );
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
package scale;
use Moose;
use Math::Round;
extends 'Feature';

has '+height'  => ( is => 'rw', isa => 'Int', lazy => 1, default => 12 );
has '+tag'     => ( is => 'ro', isa => 'Str', lazy => 1, default => 'polyline' );
has 'interval' => ( is => 'ro', isa => 'Int', default => 100 );
has 'bp_width' => ( is => 'rw', isa => 'Num' );

sub attributes {
	my $self   = shift;
	my $inc    = $self->interval;
	my $start  = $self->start;
	my $end    = $self->end;
	my $bp     = $self->bp_width;
	my $y1     = $self->y1;
	my $y2     = $self->y2;
	my $ym     = round(($y1+$y2)/2);
	my $height = $self->height;	
	
	# populate points
	my ( @x, @y );
	my $major = 0;
	for ( my $i = $start; $i <= $end; $i += $inc ) {
		my $x1 = round(($i-$start)/$bp);
		my $x2 = ($i+$inc) < $end ? round(($i-$start+$inc)/$bp) : round(($end-$start)/$bp);
		push @x, $x1, $x1, $x1, $x2;
		
		# major cross hatch in increments of 10 x $inc
		if ( not $major % 10 ) {
			push @y, $y1, $y2, $ym, $ym;
		}
		else {
			push @y, $y1 + $height / 3, $y2 - $height / 3, $ym, $ym;
		}
		$major++;
	}
	
	# make the path
	my $path = $self->svg->get_path(
		'-type'   => $self->tag,
		'x'       => \@x,
		'y'       => \@y,
		'-closed' => 0,
	);	
	
	return ( 
		%$path,
		'style' => {
			'stroke'       => $self->stroke,
			'fill'         => 'none',
			'stroke-width' => 0.1,
		}
	);	
}

##########################################################################################
package segment;
use Moose;
extends 'Feature';

has 'nearest' => ( is => 'rw', isa => 'HashRef' );

##########################################################################################
package crossover;
use Moose;
use POSIX;
use Math::Round;
extends 'Feature';

has 'segments' => ( is => 'rw', isa => 'ArrayRef[segment]' );
has 'ingroup'  => ( is => 'rw', isa => 'ArrayRef' );
has 'outgroup' => ( is => 'rw', isa => 'ArrayRef' );
has 'spacing'  => ( is => 'rw', isa => 'Int', default => 0 );

sub draw {
	my $self = shift;
	my $svg = $self->svg;
	my @ingroup = @{ $self->ingroup };
	my @segment = @{ $self->segments };
	my $space = $self->spacing;
	my $map = $self->mapping;
	
	# calculate y coordinate lookup for the outgroup, draw line
	my %height;
	{
		my $height = $self->y2 - $self->y1;
		my @outgroup = @{ $self->outgroup };
		my $inc = round( $height / scalar(@outgroup) );
		for my $i ( 0 .. $#outgroup ) {
			$height{$outgroup[$i]} = $i * $inc + $self->y1;
			
			# axis marker
			$svg->line(
				'x1' => 0,
				'x2' => $svg->{'-document'}->{'width'},
				'y1' => $height{$outgroup[$i]},
				'y2' => $height{$outgroup[$i]},
				'style' => {
					'stroke' => 'gray',
					'fill'   => 'none',
					'stroke-width' => 0.1,
				}
			);
			
			# axis label
			$svg->text(
				'x'      => 4,
				'y'      => $height{$outgroup[$i]} - 6,
				'-cdata' => $map->{$outgroup[$i]},
				'class'  => ref($self) . ' ' . $outgroup[$i],
				'style'  => { 'font' => '10px Verdana' },
			);
		}	
	}
	
	# calculate colour table and offsets for the ingroup
	my %colour;
	my %translate;
	{
		my $inc = round( 360 / scalar(@ingroup) );
		my $w = 100;
		for my $i ( 0 .. $#ingroup ) {
			my ( $r, $g, $b ) = map { round(255 * $_) } _hsv_to_rgb( $inc * $i, 1, 1 );
			my $code = "rgb($r,$g,$b)";
			$self->log->info("$ingroup[$i] = $code");
			$colour{$ingroup[$i]} = $code;
			$translate{$ingroup[$i]} = "translate($i,$i)";
			
			# legend
			$svg->line(
				'x1' => 5,
				'x2' => 10,
				'y1' => $self->y2 + 15 * $i,
				'y2' => $self->y2 + 15 * $i,
				'style' => {
					'stroke' => $code,
					'fill'   => 'none',
				}
			);
			$svg->text(
				'x' => 15,
				'y' => $self->y2 + 15 * $i + 5,
				'-cdata' => $map->{$ingroup[$i]},
				'class'  => ref($self),
				'style'  => { 'font' => '10px Verdana' },				
			);
		}
	}
	
	# iterate over segments
	for my $i ( 0 .. $#segment ) {
		my $nearest = $segment[$i]->nearest;
		for my $in ( @ingroup ) {
			my ( $x1, $x2 ) = ( $segment[$i]->x1, $segment[$i]->x2 );		
			my $out  = $nearest->{$in};
			my $y    = $height{$out};
			my $rgb  = $colour{$in};
			my %args = (
				'class'          => 'segment ' . $in,
				'transform'	     => $translate{$in},
				'stroke-linecap' => 'round',
				'style' => {
					'stroke' => $colour{$in},
					'fill'   => 'none',
				},
			);			
			
			# need to link to next segment, decrease x2
			$x2 -= $space if $i != $#segment;
			
			# make link with preceding
			if ( $i != 0 ) {
				
				# compute coordinates
				$x1 += $space;				
				my $xp1  = $segment[$i-1]->x2 - $space;
				my $xp2  = $x1;
				my $xpm1 = $xp1 + ( ($xp2-$xp1)/3 );
				my $xpm2 = $xpm1 + ( ($xp2-$xp1)/3 );
				my $yp1  = $height{ $segment[$i-1]->nearest->{$in} };
				my $yp2  = $height{ $out };
				
				# draw curve
				$svg->path( 'd' => qq{M$xp1,$yp1 C$xpm1,$yp1 $xpm2,$yp2 $xp2,$yp2}, %args );
			}			
						
			# draw line
			$svg->line(
				'x1' => $x1,
				'x2' => $x2,
				'y1' => $height{$out},
				'y2' => $height{$out},
				%args,
			);
		}	
	}	
}


sub _hsv_to_rgb ($$$) {
	my ( $h, $s, $v ) = @_;
	my ( $r, $g, $b, $i, $f, $p, $q, $t );
    if( $s == 0 ) {
    	## achromatic (grey)
        return ($v,$v,$v);
    }

    $h /= 60;                       ## sector 0 to 5
    $i = POSIX::floor( $h );
    $f = $h - $i;                   ## factorial part of h
    $p = $v * ( 1 - $s );
    $q = $v * ( 1 - $s * $f );
    $t = $v * ( 1 - $s * ( 1 - $f ) );

	if ( $i < 1 ) {
		$r = $v;
        $g = $t;
        $b = $p;
	} 
	elsif ( $i < 2 ) {
		$r = $q;
        $g = $v;
        $b = $p;
	} 
	elsif ( $i < 3 ) {
		$r = $p;
        $g = $v;
        $b = $t;
	} 
	elsif ( $i < 4 ) {
		$r = $p;
        $g = $q;
        $b = $v;
	} 
	elsif ( $i < 5 ) {
		$r = $t;
        $g = $p;
        $b = $v;
	} 
	else {
		$r = $v;
        $g = $p;
        $b = $q;
	}
	return ($r,$g,$b);
}

1;