#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_matrix';
use Bio::Phylo::Util::Logger ':simple';
use Bio::Tools::Run::StandAloneBlastPlus;
use Data::Dumper;

# process command line arguments
my $verbosity = WARN;
my $type = 'gene';
my ( $indir, $outdir, $keep, $genes, $db, $gff3 );
GetOptions(
	'indir=s'  => \$indir,
	'outdir=s' => \$outdir,
	'verbose+' => \$verbosity,
	'keep=s'   => \$keep,
	'db=s'     => \$db,
	'gff=s'    => \$gff3,
	'type=s'   => \$type,
);

# report results
VERBOSE( '-level' => $verbosity, '-class' => 'main' );
INFO "input directory: $indir";
INFO "output directory: $outdir";
INFO "keeping accessions: $keep";
INFO "using local BLAST db: $db";
INFO "using GFF3: $gff3";

# instantiate blast wrapper
my $blast = Bio::Tools::Run::StandAloneBlastPlus->new(
	'-db_data' => $db,
#	'-create'  => 1,
);

# make a lookup table to match file names against
my %keep;
for ( split /,/, $keep ) {
	$keep{$_} = 1 if $_;
}

# start reading from the directory
opendir my $dh, $indir or die $!;
while( my $entry = readdir $dh ) {

	# these should be fasta files where the numerical
	# part of the file name is the accession number
	# as in http://www.tomatogenome.net/accessions.html
	if ( $entry =~ /\.(\d+)\./ ) {
		my $accession = $1;
		
		# found an accession of interest
		if ( $keep{$accession} ) {
			INFO "going to process accession $accession";
			parse_matrix(
				'-format' => 'fasta',
				'-type'   => 'dna',
				'-file'   => "$indir/$entry",
			)->visit(sub{ write_seq( shift, $accession ) });
		}
	}
}

sub write_seq {
	my ( $seq, $accession ) = @_;
	
	# do a blast and a scan of the GFF3 to find the gene id :-/
	my ( $gene_id, $strand ) = @{ find_gene($seq) };
	
	# maybe do a reverse complement
	my $raw = $seq->get_char;
	if ( $strand eq '-' ) {
		$raw =~ tr/ACGT/TGCA/;
		$raw = reverse $raw;
	}
	
	# create outfile name
	my $file = $outdir . '/' . $gene_id . '.fa';
	DEBUG "going to append to file $file";
	
	# write output
	open my $fh, '>>', $file or die $!;
	print $fh '>', $accession, "\n", $raw, "\n";
}

{
	my %cache;
	sub find_gene {
		my $seq = shift;
		my $index = $seq->get_name;
		
		# we've done this already
		if ( $cache{$index} ) {
			DEBUG "already blasted an instance of index $index";
			return $cache{$index};
		}
		
		# need to do a blast, then scan the GFF file for the region
		INFO "going to blast $index";
		my $result = $blast->blastn( '-query' => $seq );
		
		# iterate over hits
		while( my $hit = $result->next_hit ) {
			my $hitname = $hit->name;
			$hitname =~ s/^lcl\|//; # strip prefix
			INFO "found hit '$hitname'"; # should be chromosome

			# iterate over high scoring pairs
			while( my $hsp = $hit->next_hsp ) {
				my $start = $hsp->start("subject");
				my $end   = $hsp->end("subject");
				INFO "HSP range is $start..$end (".($end-$start)." bp)";

				# having found the coordinates of the focal
				# hit into the reference we will now
				# query the BAM file for the same coordinates
				($cache{$index}) = scan_gff(
					'start' => $start,
					'end'   => $end,
					'chr'   => $hitname, 
				);
				return $cache{$index} if $cache{$index};
			}
		}
	}
}

sub scan_gff {
	my %args = @_;

	# 0-based column numbers in gff3
	my $chr_idx    = 0;
	my $type_idx   = 2;
	my $start_idx  = 3;
	my $end_idx    = 4;
	my $strand_idx = 6;
	my $meta_idx   = 8;
	# SL2.40ch00      ITAG_eugene     gene    16437   18189   .       +       .       Alias=Solyc00g005000;ID=gene:Solyc00g005000.2;Name=Solyc00g005000.2;from_BOGAS=1;length=1753

	# start reading
	INFO "going to scan gff3 $gff3 for a $type on chromosome $args{chr} near $args{start}..$args{end}";
	my @result;	
	open my $fh, '<', $gff3 or die $!;
	while(<$fh>) {
		chomp;
		my @line = split /\t/, $_;
=pod
		
		# we're on the right chromosome
		if ( $line[$chr_idx] eq $args{'chr'} ) {
			
			# we've got the right type
			if ( $line[$type_idx] eq $type ) {
			
				# we're within range
				if ( $line[$start_idx] <= $args{'start'} && $line[$end_idx] >= $args{'end'} ) {
				
					# the metadata column has the right format
					if ( $line[$meta_idx] =~ m/;ID=$type:([^;]+)/ ) {
						my $id = $1;
						my $strand = $line[$strand_idx];
						INFO "found $id on $args{'chr'} at $line[$start_idx]..$line[$end_idx]";
						return [$id, $strand];
					}					
				}
			}
		}

=cut

		if ( $line[$start_idx] && $line[$start_idx] <= $args{start} && $line[$end_idx] >= $args{end} && $line[$type_idx] eq $type ) {
			# the metadata column has the right format
			if ( $line[$meta_idx] =~ m/;ID=$type:([^;]+)/ ) {
				my $id = $1;
				my $strand = $line[$strand_idx];
				INFO "found $id on $args{'chr'} at $line[$start_idx]..$line[$end_idx]";
				push @result, [$id, $strand];
			}
		}
	}
	if ( scalar(@result) > 1 ) {
		WARN "more than one exact coordinate match, using first one!\n".Dumper(\@result);
	}
	return @result;
}