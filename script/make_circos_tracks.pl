#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util 'sum';
use Bio::Phylo::Util::Logger ':levels';

# global variables
my $Log; # logger

# scans a GFF3 file for one or more locus identifiers 
# (i.e. whatever value is associated with the ID=... key
# in the description column). returns a list of hashes
# with the gene ID, and chr, start, end coordinates
sub read_gff3 {
	my ( $gff3, @genes ) = @_;

	# column numbers in gff3
	my $chr   = 0;
	my $type  = 2;
	my $start = 3;
	my $end   = 4;
	my $desc  = 8;

	# we return the ordered list but use
	# the hash for lookups
	my @result;
	my %lookup = map { $_ => 1 } @genes;

	# start reading the file
	$Log->info("going to read gene locations from $gff3");
	$Log->info("genes: \n" . join("\n",@genes));
	open my $fh, '<', $gff3 or die $!;
	while(<$fh>) {
		chomp;
		my @record = split;

		# capture the locus id from the description column
		if ( $record[$desc] and $record[$desc] =~ m/ID=([^;]+)/ ) {
			my $id = $1;

			# found a locus of interest, store the location
			if ( $lookup{$id} ) {
				push @result, {
					'gene'  => $id,
					'start' => $record[$start],
					'end'   => $record[$end],
					'chr'   => $record[$chr],
				};
				$Log->info(sprintf('found %s at %s:%i-%i',$id,@record[$chr,$start,$end]));
			}
		}
	}

	# return result
	return @result;
}

# validates a command line argument
sub check_arg {
	my %args = @_; # name, value, file, list, dir
	use Data::Dumper;
	warn Dumper(\%args);

	# is there a value?
	if ( not $args{'value'} ) {
		warn "Need -$args{'name'} <value> argument\n";
	}

	# split if need be
	my @value = $args{'list'} ? grep { /\S/ } split /,/, $args{'value'} : ( $args{'value'} );

	# check files and directories
	for my $v ( @value ) {
		if ( $args{'file'} and not -f $v ) {
			warn "Need -$args{name} <file1,file2> argument\n";
		}
		if ( $args{'dir'} and not -d $v ) {
			warn "Need -$args{dir} <dir1,dir2> argument\n";
		}
	}
	return @value;

}

# process and validate all command line arguments
sub check_args {

	# process command line arguments
	my $gff3;             # gff3 file
	my $workdir;          # working directory
	my $refseq;           # reference genome in FASTA format
	my $verbosity = WARN; # default logging level
	my $mapdamage;        # output tables from mapDamage
	my $genes;            # gene identifiers in GFF3
	my $bams;             # BAM files to compute coverage
	my $size = 100000;    # bin size for cover
	my $readwindow = 10;  # window size for misincorporation
	GetOptions(
		'verbose+'     => \$verbosity,
		'gff3=s'       => \$gff3,
		'mapdamage=s'  => \$mapdamage,
		'workdir=s'    => \$workdir,
		'refseq=s'     => \$refseq,
		'genes=s'      => \$genes,
		'bams=s'       => \$bams,
		'size=i'       => \$size,
		'readwindow=i' => \$readwindow,
	);

	# instantiate logger
	$Log = Bio::Phylo::Util::Logger->new(
		'-level' => $verbosity,
		'-class' => 'main',
	);

=begin pod

	# check individual arguments
	check_arg( 'name' => 'gff3', 'value' => $gff3 );
        check_arg( 'name' => 'refseq',  'value' => $refseq, 'file' => 1 );
	check_arg( 'name' => 'workdir', 'value' => $workdir, 'dir' => 1 );

	# check and expand lists
	my @genes = check_arg( 'name' => 'genes', 'value' => $genes, 'list' => 1 );
	my @mapdamage = check_arg( 
		'name'  => 'mapdamage', 
		'value' => $mapdamage, 
		'dir'   => 1, 
		'list'  => 1 
	);
	my @bams = check_args(
		'name'  => 'bams',
		'value' => $bams,
		'file'  => 1,
		'list'  => 1,
	);
=cut

my @genes = split /,/, $genes;
my @mapdamage = split /,/, $mapdamage;
my @bams = split /,/, $bams;

	# done
	return 
		'gff3'       => $gff3,
		'workdir'    => $workdir,
		'mapdamage'  => \@mapdamage,
		'refseq'     => $refseq,
		'genes'      => \@genes,
		'bams'       => \@bams,
		'size'       => $size,
		'readwindow' => $readwindow;
}

# write a circos karyotype track to a provided file name. 
# second argument is a list of array refs where each array
# ref's first elt is the chromosome label, the second the
# length of the chromosome ( [ chr1, 213127 ], [ chr2, 72346 ] );
sub write_karyotype {
	my ( $karyotype, @chr ) = @_;

	# open file handle
	$Log->info("going to write karyotype to $karyotype");
	open my $fh, '>', $karyotype or die $!;

	# iterate over chromosomes
	for my $i ( 0 .. $#chr ) {
		my ( $label, $end ) = @{ $chr[$i] };
		printf $fh "chr - %s %s 0 %i chr%i\n", $label, $label, $end, $i + 1;
	}
}

# reads a fasta file, returns a list as required by the input
# for write_karyotype
sub read_fasta {
	my $file = shift;
	$Log->info("going to read chromosomes from FASTA $file");
	open my $fh, '<', $file or die $!;
	my ( $current, %length, @result );
	while(<$fh>) {
		chomp;
		if ( /^>(\S+)/ ) {
			my $id = $1;
			$current = $id;
			$length{$current} = 0;
			push @result, $current;
			$Log->debug("going to read chromosome $current");
		}
		else {
			my $seq = $_;
			$seq =~ s/\s//g;
			$length{$current} += length $seq;
		}
	}
	return map { [ $_ => $length{$_} ] } @result;
}

# writes a label track to the provided file name using
# the output as formatted by read_gff3
sub write_labels {
	my ( $outfile, @locations ) = @_;
	$Log->info("going to write label track to $outfile");
	open my $fh, '>', $outfile or die $!;
	for my $l ( @locations ) {
		my %location = %{ $l }; 
		printf $fh "%s %i %i %s\n", @location{qw(chr start end gene)};
	}
}

# writes a continuous-valued 2d track to the provided file name. the
# data to be written is a hash reference whose keys are chromosome names
# and values are array references with the raw data. the data is written
# out as binned over as many elements as there are in the arrays.
sub write_continuous_track {
	my ( $outfile, $data, $size, @chr ) = @_;
	$Log->info("going to write continuous data track to $outfile");
	open my $fh, '>', $outfile or die $!;
	for my $c ( @chr ) {
		my $chromo = $c->[0];
		for my $i ( 0 .. $#{ $data->{$chromo} } ) {
			my $start = $i * $size;
			my $end   = ( $start + $size ) < $c->[1] ? $start + $size : $c->[1];
			my $value = $data->{$chromo}->[$i];
			printf $fh "%s %i %i %f\n", $chromo, $start, $end, $value; 
		}
	}
}

# reads the coverage in a BAM file, averages this in bins of $size for each chromosome
sub read_coverage {
	my ( $bam, $refseq, $size ) = @_;
	$Log->info("going to read sliding window (size: $size) coverage from BAM file $bam");
	my %result;

	# this is going to create a lot of data, so we will pass by reference
	open my $fh, "samtools mpileup -BQ0 -d10000000 -f $refseq $bam |" or die "Can't run samtools: $!";
	my @window;
	my $current;
	while(<$fh>) {
		chomp;
		my @record = split /\s+/, $_;
		my ( $chr, $cover ) = @record[0,3];

		# starting a new chromosome
		if ( not $result{$chr} ) {

			# finish the previous one
			push @{ $result{$current} }, sum(@window)/scalar(@window) if $current;
			@window = ();

			# initialize the new one
			$result{$chr} = [];
			$current = $chr;
			$Log->info("going to compute sliding window coverage for chromosome $chr");
		}
		push @window, $cover;

		# compute average over window
		if ( $size == scalar @window ) {
			push @{ $result{$current} }, sum(@window)/$size;
			@window = ();
			$Log->info("done window ".scalar(@{ $result{$current} }));
		}
		
	}
	return \%result;
}

# reads mapDamages fragmentation and misincorporation results, returns these
# per chromosome
sub read_mapdamage {
	my ( $dir, $readwindow ) = @_;
	my ( %fragmentation, %misincorporation );

	$Log->info("going to read mapDamage results in dir '$dir'");
	opendir my $dh, $dir or die "Can't open dir handle: $!";
	while( my $entry = readdir $dh ) {
		if ( $entry =~ m/dnaFrag_(.+?)_(\d+?)_(\d+?)\.txt/ ) {
			my ( $chromo, $around, $length ) = ( $1, $2, $3 );
			$fragmentation{$chromo} = [ read_fragmentation( "$dir/$entry", $around, $length ) ];
		}
		elsif ( $entry =~ m/nuclComp_(.+?)\.txt/ ) {
			my $chromo = $1;
			$misincorporation{$chromo} = [ read_misincorporation( "$dir/$entry", $readwindow ) ];
		}
	}
	return \%fragmentation, \%misincorporation;
}

# reads the fraction of C>T substitutions in the first $size read positions
sub read_misincorporation {
	my ( $file, $readwindow ) = @_;
	$Log->info("going to read fraction of C>T substitions in first 10 read positions from $file");
	my @header;
	my %transition;
	open my $fh, '<', $file or die $!;
	my $counter = 0;
	RECORD: while(<$fh>) {
		chomp;
		if ( not @header ) {
			@header = split;
		}
		else {
			my @record = split;
			for my $i ( 0 .. $#header ) {
				if ( $header[$i] =~ /^[ACGT]>[ACGT]$/ ) {
					$transition{$header[$i]} += $record[$i];		
				}
			}
			last RECORD if $counter++ == $readwindow;
		}
	}
	return $transition{"C>T"} / sum(values(%transition));
}

# computes 1 base upstream excess of purines
sub read_fragmentation {
	my ( $file, $around, $length ) = @_;
	$Log->info("going to read fragmentation data from $file");
	my %result;
	my @header;
	open my $fh, '<', $file or die $!;
	while(<$fh>) {
		chomp;
		if ( not @header ) {
			@header = split;
		}
		else {
			my @record = split;
			my $strand = 0;
			my $pos    = 1;
			my $tot    = 8;
			my $A      = 2;
			my $G      = 4;
			if ( $record[$pos] == -1 ) {
				$Log->info("computing 1 base upstream excess of purines on $record[$strand] strand"); 
				$result{ $record[$strand] } = [ $record[$A] + $record[$G], $record[$tot] ];
			}
		}
	}
	my $totPurines = $result{ "0"  }->[0] + $result{ "16" }->[0];
	my $totTot     = $result{ "16" }->[1] + $result{ "16" }->[1]; 
	return $totPurines / $totTot;
}

sub main {
	my %args = check_args();

	# read karyotype from FASTA, write to karyotype file
	my @chromosomes = read_fasta( $args{'refseq'} );
	my $karyotype = $args{'workdir'} . '/karyotype.txt';
	write_karyotype( $karyotype, @chromosomes );

        # write circos config file
        my $circos_conf = $args{'workdir'} . '/circos.conf';
        open my $fh, '>', $circos_conf or die $!;

	# write header
        print $fh <<"HEADER";
<colors>
	<<include etc/colors.conf>>
</colors>
<fonts>
	<<include etc/fonts.conf>>
</fonts>
<<include ideogram.conf>>
<<include ticks.conf>>
karyotype = $karyotype
<image>
	dir = .
	file = circos.png

	# radius of inscribed circle in image
	radius         = 1500p
	background     = white

	# by default angle=0 is at 3 o'clock position
	angle_offset   = -90
</image>
<plots>
HEADER

        # read gene locations, write to labels file
        my @locations = read_gff3( $args{'gff3'}, @{ $args{'genes'} } );
	my $labels = $args{'workdir'} . '/labels.txt';
	write_labels( $labels, @locations );
	print $fh <<"LABELS";

<plot>
	type   = text
	color  = black
	file   = $labels

	# on tick scale
	r0 = 1r
	r1 = 1r+200p

	show_links     = yes
	link_dims      = 0p,0p,50p,0p,10p
	link_thickness = 2p
	link_color     = red
	label_size   = 24p
	label_font   = condensed
	padding  = 0p
	rpadding = 0p
</plot>
LABELS

	# read coverage, write 2d tracks
	my $radius = 0.9;

=begin pod

	for my $i ( 0 .. $#{ $args{'bams'} } ) {
		my $bam = $args{'bams'}->[$i];
		my $coverage = read_coverage( $bam, @args{qw[refseq size]} ); # returns hash ref
		my $coverfile = $args{'workdir'} . "/coverage${i}.txt";
		write_continuous_track( $coverfile, $coverage, $args{'size'}, @chromosomes );
		$radius -= 0.1;
		print $fh <<"COVERAGE";
<plot>
	type    = heatmap
	file    = $coverfile
	color   = spectral-11-div
	r0      = ${radius}r
	r1      = ${radius}r+25p
	stroke_thickness = 1
	stroke_color     = black
</plot>
COVERAGE

	}

=cut

	# read mapdamage tables, write 2d tracks
	my %radius = (
		'frag' => $radius,
		'mis'  => $radius - 0.1 * scalar(@{ $args{'mapdamage'} }),
	);
	for my $i ( 0 .. $#{ $args{'mapdamage'} } ) {
		my $md = $args{'mapdamage'}->[$i];
		my ( $fragmentation, $misincorporation ) = read_mapdamage($md, $args{'readwindow'});
		$radius{'frag'} -= 0.1;
		$radius{'mis'} -= 0.1;		

		# write fragmentation track
		my $frag_file = $args{'workdir'} . "/fragmentation${i}.txt";
		write_continuous_track( $frag_file, $fragmentation, $args{'size'}, @chromosomes );
		print $fh <<"FRAGMENTATION";
<plot>
	type      = histogram
	file      = $frag_file
	r1        = $radius{'frag'}r+25p
	r0        = $radius{'frag'}r
	stroke_type = outline
	thickness   = 1
	color       = black
	extend_bin  = yes
	<axes>
		<axis>
			spacing   = 5p
			color     = lgrey
			thickness = 1
		</axis>
	</axes>
</plot>
FRAGMENTATION

		# write misincorporation track
		my $mis_file = $args{'workdir'} . "/misincorporation${i}.txt";
		write_continuous_track( $mis_file, $misincorporation, $args{'size'}, @chromosomes );
		print $fh <<"MISINCORPORATION";
<plot>
        type      = histogram
        file      = $mis_file
        r1        = $radius{'mis'}r+25p
        r0        = $radius{'mis'}r
        stroke_type = outline
        thickness   = 1
        color       = black
        extend_bin  = yes
        <axes>
                <axis>
                        spacing   = 5p
                        color     = lgrey
                        thickness = 1
                </axis>
        </axes>
</plot>
MISINCORPORATION

	}

	# print footer
	print $fh "</plots>\n";
}
main();
