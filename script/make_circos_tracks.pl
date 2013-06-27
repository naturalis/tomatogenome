#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use List::Util qw'sum shuffle';
use Bio::DB::Sam;
use Bio::Phylo::Util::Logger ':levels';

# global variables
my $Log; # logger

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
	my $size = 1_000_000; # bin size for cover
	my $readwindow = 10;  # window size for misincorporation
	my $increment = 0.07; # radius increments for circos
	my $compute;          # compute data tracks
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
		'increment=f'  => \$increment,
		'compute'      => \$compute,
	);

	# instantiate logger
	$Log = Bio::Phylo::Util::Logger->new(
		'-level' => $verbosity,
		'-class' => 'main',
	);

	# split arguments that are comma-separated lists, filter
	# out all null values
	my @genes     = grep { /\S/ } split /,/, $genes;
	my @mapdamage = grep { /\S/ } split /,/, $mapdamage;
	my @bams      = grep { /\S/ } split /,/, $bams;

	# done
	return 
		'gff3'       => $gff3,
		'workdir'    => $workdir,
		'mapdamage'  => \@mapdamage,
		'refseq'     => $refseq,
		'genes'      => \@genes,
		'bams'       => \@bams,
		'size'       => $size,
		'readwindow' => $readwindow,
		'increment'  => $increment,
		'compute'    => $compute;
}

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

# reads a fasta file, returns a list as required by the input
# for write_karyotype, i.e. a list of tuples where each first
# element is a chromosome ID (actually, the first word on the
# FASTA definition line) and each second element is the size
# of that chromosome.
sub read_fasta {
	my ( $file, $size ) = @_;
	
	# open file handle
	$Log->info("going to read chromosomes from FASTA $file");
	open my $fh, '<', $file or die $!;
	
	# iterate over lines
	my ( $this, $seq, %length, %gc, %n, %ag, @result );
	while(<$fh>) {
		chomp;
		
		# line is a FASTA definition line, capture the first word
		if ( /^>(\S+)/ ) {
			my $id = $1;
			
			# process previous record
			if ( $this ) {
				$length{$this} = length $seq;
				
				# compute windows of GC content, missing data and purine content
				($gc{$this},$n{$this},$ag{$this}) = process_seq_segments( $seq, $size );
				$seq = '';
			}
			
			# move on to the current record
			$this = $id;
			push @result, $this;
			$Log->info("going to read chromosome $this");
		}
		
		# line is a sequence line
		else {
			s/\s//g;
			$seq .= $_;
		}
	}
	
	# process last record
	$length{$this} = length $seq;
	( $gc{$this}, $n{$this}, $ag{$this} ) = process_seq_segments( $seq, $size );
	$seq = '';	
	
	# create and return ordered list of lists
	return map { [ 
		$_,            # chromosome ID
		$length{ $_ }, # integer: length of this chromosome
		$gc{     $_ }, # list of floats: GC content in windows of $size
		$n{      $_ }, # list of floats: missing data fraction in windows of $size
		$ag{     $_ }, # list of floats: purine content in windows of $size
	
	] } @result;
}

# calculate percentage of GC content, missing bases (nN) and purines (aAgG) for
# sequence $seq over bins of $size. returns array refs of percentages, i.e. values
# between 0..100
sub process_seq_segments {
	my ( $seq, $size ) = @_;
	my ( @gc, @n, @purines );

	# iterate over segments		
	for my $window ( unpack("(A$size)*", $seq) ) {
	
		# compute GC content
		my $g = $window =~ tr/gG/gG/;
		my $c = $window =~ tr/cC/cC/;
		push @gc, ( ( $g + $c ) / length $window ) * 100;
		
		# compute Nn content
		my $n = $window =~ tr/nN/nN/;
		push @n, ( $n / length $window ) * 100;
		
		# compute purine content
		my $a = $window =~ tr/aA/aA/;
		push @purines, ( ( $a + $g ) / length $window ) * 100;
		
	}
	return \@gc, \@n, \@purines;
}

# reads the coverage in a BAM file, averages this in bins of $size for each chromosome
sub read_coverage {
	my ( $bam, $refseq, $size, @chromosomes ) = @_;
	$Log->info("going to read sliding window (size: $size) coverage from BAM file $bam");

    # instantiate helper objects
    my $db = Bio::DB::Sam->new(
        '-bam'   => $bam,
        '-fasta' => $refseq,
    );
	
	# iterate over chromosome tuples
	my ( %cover, %qual, %bases );
	for my $tuple ( @chromosomes ) {
		my ( $chromo, $chr_size ) = @{ $tuple };
		$cover{ $chromo } = [];
		$qual{  $chromo } = [];
		$bases{ $chromo } = [];
		
		# read reference chromosome
		$Log->info("going to read chromosome $chromo from reference $refseq");
		my $seq;
		open my $fh, '<', $refseq or die $!;
		while(<$fh>) {
			chomp;
			if ( /^>${chromo}\s*$/ ) {
				$seq = '';
				next;
			}
			$seq .= $_ if defined $seq;
		}		
		
		# iterate over bins
		for ( my $start = 1; $start <= $chr_size; $start += $size ) {
		
			# calculate interval end
			my $end = ( $start + $size ) < $chr_size ? $start + $size : $chr_size;
			
			# scan the BAM file
			my ( $qual, $cover, $locations ) = compute_coverage_and_quality(  
				'chromo' => $chromo,
				'start'  => $start,
				'end'    => $end,
				'db'     => $db,
			);

			# store results
			push @{ $cover{$chromo} }, $cover;
			push @{ $qual{$chromo} }, $qual;
			
			# get bases from reference chromosome
			push @{ $bases{ $chromo } }, compute_upstream_bias(
				'seq'       => \$seq, # this is a whole chromosome, pass by reference
				'locations' => $locations,
				'samples'   => 100,
				'chromo'    => $chromo,
			);			

		}	
	}
	return \%cover, \%qual, \%bases;
}

# for a sample of locations, fetches the base of the refseq
# in order to compute upstream AG bias due to fragmentation
sub compute_upstream_bias {
    my %args = @_;

	# we will use this to complement bases on the - strand
    my %complement = ( 
    	'A' => 'T', 
    	'C' => 'G', 
    	'G' => 'C', 
    	'T' => 'A',
    );
    
    my $sample; # tracks number of samples
    my %locations = %{ $args{'locations'} }; # keys are locations, values are strand
    my %bases; # just for logging

	# this iterates over raw locations (single integers)
    LOCATION: for my $loc ( shuffle(keys %locations) ) {
    
    	# locations are 1-based so decrement for substr
        my $raw_base = uc substr ${ $args{'seq'} }, $loc - 1, 1;
        next LOCATION if $raw_base !~ /[ACGT]/; # skip missing
        
        # bases on the +1 strand need to be complemented because
        # we are now looking at the refseq and we want what was
        # on its complement in the fragmented template
        my $base = $locations{$loc} < 1 ? $complement{$raw_base} : $raw_base;
        $bases{$base}++;
        
        last LOCATION if $args{'samples'} == ++$sample;
    }
    
    # return frequencies of bases
	my %frequencies = map { $_ => $bases{$_} / $sample } keys %bases;
	use Data::Dumper;
	$Log->info(Dumper(\%frequencies));
    return \%frequencies;
}

# iterates over reads inside interval, computes 
# coverage and mapping quality and calculates
# the locations of 1 base upstream from start of read
sub compute_coverage_and_quality {
    my %args = @_;
    
    # make interval for get_features_by_location
    my %loc = (
        '-seq_id' => $args{'chromo'},
        '-start'  => $args{'start'},
        '-end'    => $args{'end'},      
    );  

    # we compute coverage by counting the number of aligned reads with $cover
    # and computing their average length by summing the lengths
    my $reads  = 0;
    my $length = 0;

    # we record all distinct starting positions of aligned reads with %starts. 
    # subsequently we will calculate the base composition one base upstream from 
    # that on the reference. if there is fragmentation this should show an AG 
    # bias on the reference.
    my %starts;

    # we record all qualities in @quals, then average that
    my @quals;

    # iterate over alignments within interval
    ALN: for my $aln ( $args{'db'}->get_features_by_location(%loc) ) {

        # start and end coordinates of the alignment on the reference
        my ( $aln_start, $aln_end ) = ( $aln->start, $aln->end );
        my $strand = $aln->query->strand; # strandedness of the read
        
        # if the read is on the +1 strand the upstream base is aln start - 1,
        # if it is on the -1 strand the upstream base is aln end + 1 
        # and we need to complement
        my $upstream_base = $strand < 1 ? $aln_end + 1 : $aln_start - 1;
        $starts{ $upstream_base } = $strand;
    
        # running tally of qualities
        push @quals, $aln->qual;
    
        # running tally of alignment lengths
        $length += ( $aln_end - $aln_start );
    
        # number of reads
        $reads++;
        $Log->info("Read: $reads") unless $reads % 100000;
    }

    # compute coverage
    my $cover = ( $reads * ( $length / $reads ) ) / ( $args{'end'} - $args{'start'} );

    # compute average quality
    my $qual = sum(@quals) / scalar(@quals);
    
    # communicate results
    $Log->info( "Total reads: $reads" );
    $Log->info( "Average cover: $cover" );
    $Log->info( "Average mapping quality: $qual" );
    $Log->info( "Found ".scalar(keys(%starts))." distinct upstream bases" );

    return $qual, $cover, \%starts;
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
		printf $fh "chr - %s %s 1 %i chr%i\n", $label, $label, $end, $i + 1;
	}
}

# writes a label track to the provided file name using
# the output data structure returned by read_gff3
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
	
	# open file handle
	$Log->info("going to write continuous data track (bin size: $size) to $outfile");
	open my $fh, '>', $outfile or die $!;
	
	# iterate over chromosome data structure as returned by read_fasta
	for my $c ( @chr ) {
		my $chromo  = $c->[0]; # i.e., chromosome name, such as "SL2.40ch00"
		my $chrsize = $c->[1]; # chromosome size (in bp, from ref genome)
		
		# iterate over data bins for this chromosome
		my $last_index = $#{ $data->{$chromo} };
		for my $i ( 0 .. $last_index ) {
			my $start = $i * $size;

			# calculate the end coordinate:			
			# there is only one data element, assume it applies to the whole chromosome
			# OR the start + increment exceeds the chromosome size
			my $end;
			if ( $last_index == 0 || ( $start + $size ) > $chrsize ) {
				$end = $chrsize;
			}
			else {
				$end = $start + $size;
			}

			# print the datum
			printf $fh "%s %i %i %f\n", $chromo, $start, $end, $data->{$chromo}->[$i]; 
		}
	}
}

# similar to write_continuous_track, but each datum is expected to consist
# of a hash whose keys are upper case bases (ACGT), the values of which
# are listed as comma-separated, for stacked histograms
sub write_stacked_track {
	my ( $outfile, $data, $size, @chr ) = @_;
	
	# open file handle
	$Log->info("going to write continuous data track (bin size: $size) to $outfile");
	open my $fh, '>', $outfile or die $!;
	
	# iterate over chromosome data structure as returned by read_fasta
	for my $c ( @chr ) {
		my $chromo  = $c->[0]; # i.e., chromosome name, such as "SL2.40ch00"
		my $chrsize = $c->[1]; # chromosome size (in bp, from ref genome)
		
		# iterate over data bins for this chromosome
		my $last_index = $#{ $data->{$chromo} };
		for my $i ( 0 .. $last_index ) {
			my $start = $i * $size;

			# calculate the end coordinate:			
			# there is only one data element, assume it applies to the whole chromosome
			# OR the start + increment exceeds the chromosome size
			my $end;
			if ( $last_index == 0 || ( $start + $size ) > $chrsize ) {
				$end = $chrsize;
			}
			else {
				$end = $start + $size;
			}
			
			# the data hash is keyed on uppercase nucleotide symbols, values are freqs
			my %bases = %{ $data->{$chromo}->[$i] };
			my $stacked = join ',', map { $_ || 0.0 } @bases{qw[A C G T]};

			# print the datum
			printf $fh "%s %i %i %s\n", $chromo, $start, $end, $stacked; 
		}
	}
}


# the write_circos_* subroutines populate a circos.conf file with
# karyotype, label, heatmap and histogram tracks. they contain 
# hardcoded circos config language, which should be made as 
# tweakable as need be.
sub write_circos_header {
	my ( $circos_conf, $karyotype ) = @_;
	open my $fh, '>', $circos_conf or die $!;

	# write header
	print $fh <<"HEADER";
<<include etc/housekeeping.conf>>
<colors>
	<<include etc/colors.conf>>
</colors>
<fonts>
	<<include etc/fonts.conf>>
</fonts>
<ideogram>		
	<spacing>
		default = 30u
		break   = 2u
		axis_break_style = 2		
		<pairwise SL2.40ch12,SL2.40ch00>
			spacing = 10u
		</pairwise>		
		<break_style 2>
			thickness        = 50p
			stroke_color     = black
			stroke_thickness = 5p
		</break>		
	</spacing>	
	show           = yes
	thickness      = 25p
	fill           = yes
	fill_color     = black
	radius         = 0.75r
	show_label     = yes
	label_font     = condensed
	label_radius   = dims(ideogram,radius_outer) + 20p
	label_size     = 24p
	label_parallel = yes
	show_bands            = yes
	fill_bands            = yes
	band_stroke_thickness = 0
	band_stroke_color     = black
	band_transparency     = 4
</ideogram>
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

	return $fh;
}

sub write_circos_label_config {
	my ( $fh, $labels ) = @_;
	print $fh <<"LABELS";
<plot>
	type   = text
	color  = black
	file   = $labels
	r0 = 1r
	r1 = 1r+500p
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
}

sub write_circos_histogram_config {
	my ( $fh, $datafile, $radius, $min, $max ) = @_;
	my $spacing = int( ( $max - $min ) / 10 );
	print $fh <<"HISTOGRAM";
<plot>
	type      = histogram
	file      = $datafile
	r1        = ${radius}r+70p
	r0        = ${radius}r
	stroke_type = outline
	thickness   = 1
	color       = black
	extend_bin  = yes
	min         = $min
	max         = $max
	<axes>
		<axis>
			spacing   = $spacing
			color     = lgrey
			thickness = 1
		</axis>
	</axes>
	<rules>
		<rule>
			condition  = 1
			fill_color = eval(sprintf("spectral-9-div-%d",remap_int(var(value),$min,$max,1,9)))
		</rule>
	</rules>
</plot>
HISTOGRAM
}

sub write_circos_stacked_histogram_config {
	my ( $fh, $datafile, $radius, $min, $max ) = @_;
	print $fh <<"STACKED";
<plot>
	type       = histogram
	file       = $datafile
	r1         = ${radius}r+70p
	r0         = ${radius}r
	color      = white
	fill_color = red,green,black,blue
	thickness  = 0
	extend_bin = no
	axis       = no
</plot>
STACKED
}

sub main {
	my %args = check_args();

	# read karyotype from FASTA, write to karyotype file
	my @chromosomes = read_fasta( $args{'refseq'}, $args{'size'} );
	my $karyotype = $args{'workdir'} . '/karyotype.txt';
	write_karyotype( $karyotype, @chromosomes );

	# initialize circos config file
	my $circos_conf = $args{'workdir'} . '/circos.conf';
	my $fh = write_circos_header( $circos_conf, $karyotype );

    # handle gene locations
	my $labelfile = $args{'workdir'} . '/candidates.txt';
	if ( $args{'compute'} ) {    
    	my @locations = read_gff3( $args{'gff3'}, @{ $args{'genes'} } );
		write_labels( $labelfile, @locations );
	}
	write_circos_label_config( $fh, $labelfile );
	
	# these govern the radius of the outermost data track and the
	# decrements going from outermost to innermost
	my $radius = 1 - $args{'increment'};
	my $inc    = $args{'increment'};	
	
	# handle GC content
	my %gc = map { $_->[0] => $_->[2] } @chromosomes;
	my $gcfile = $args{'workdir'} . '/gccontent.txt';
	write_continuous_track( $gcfile, \%gc, $args{'size'}, @chromosomes );
	# decrement data track radius
	$radius -= $inc;		
	write_circos_histogram_config( $fh, $gcfile, $radius, 25, 40 );	

	# read coverage, write histograms tracks
	for my $i ( 0 .. $#{ $args{'bams'} } ) {
			
		# data files
		my $coverfile = $args{'workdir'} . "/coverage${i}.txt";
		my $qualfile  = $args{'workdir'} . "/quality${i}.txt";
		my $basesfile = $args{'workdir'} . "/bases${i}.txt";
		
		# this step takes long
		if ( $args{'compute'} ) {
	
			# read coverage for BAM file $i, returns hash ref
			my $bam = $args{'bams'}->[$i];
			my ( $cover, $qual, $bases ) = read_coverage( $bam, @args{qw[refseq size]}, @chromosomes ); 
						
			# write data tracks
			write_continuous_track( $coverfile, $cover, $args{'size'}, @chromosomes );
			write_continuous_track( $qualfile, $qual, $args{'size'}, @chromosomes );
			write_stacked_track( $basesfile, $bases, $args{'size'}, @chromosomes );
		}
		
		# decrement data track radius
		$radius -= $inc * 1.6;		
		write_circos_histogram_config( $fh, $coverfile, $radius, 0, 10 );
		# decrement data track radius
		$radius -= $inc;
		write_circos_histogram_config( $fh, $qualfile, $radius, 30, 60 );
		# decrement data track radius
		$radius -= $inc;
		write_circos_stacked_histogram_config( $fh, $basesfile, $radius, 0, 1 );			
	}

	# print footer
	print $fh "</plots>\n";
}
main();
