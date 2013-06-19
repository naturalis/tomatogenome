#!/usr/bin/perl
use strict;
use warnings;
use File::Temp 'tempfile';
use Getopt::Long;
use Bio::DB::GenBank;
use Bio::Phylo::Factory;
use Bio::Tools::CodonTable;
use Bio::Tools::Run::RemoteBlast;
use Bio::Phylo::Matrices::Matrix;
use Bio::Phylo::IO 'parse_matrix';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Tools::Run::Alignment::Muscle;

# workaround for buggy bioperl wrapper
my @switches = @Bio::Tools::Run::Alignment::Muscle::MUSCLE_SWITCHES;
@Bio::Tools::Run::Alignment::Muscle::MUSCLE_SWITCHES = grep { $_ ne 'profile' } @switches;

# global variables
my $Log; # logger

# process command line arguments
sub check_args {
	my $verbosity = WARN; # log level
	my $divergers; # spreadsheet with pinvar per file
	my $cutoff; # max pinvar to consider
	my $reference; # reference taxon to use for BLAST searches
	
	GetOptions(
		'verbose+'    => \$verbosity,
		'divergers=s' => \$divergers,
		'cutoff=f'    => \$cutoff,
		'reference=s' => \$reference,
	);

	# instantiate logger
	$Log = Bio::Phylo::Util::Logger->new(
		'-level' => $verbosity,
		'-class' => 'main',
	);

	# communicate settings
	$Log->info("divergence spreadsheet is $divergers");
	$Log->info("cutoff is $cutoff");
	$Log->info("reference is $reference");
	
	return $divergers, $reference, $cutoff;
}

# reads tab-separated file, first column is a file name, second column is the
# proportion of invariant sites. returns list of files whose proportion of 
# invariant sites is below a specified cutoff.
sub read_spreadsheet {
	my ($divergers,$cutoff) = @_;
	my @files;

	# start reading from file
	open my $fh, '<', $divergers or die $!;
	$Log->info("going to read from $divergers, picking alignments where pinvar <= $cutoff");
	while(<$fh>) {
		chomp;
		my ( $infile, $pinvar ) = split;
	
		# found a file of interest
		if ( $pinvar <= $cutoff ) {
			$Log->info("$infile has enough variation");
			push @files, $infile;
		}
		else {
			$Log->info("not enough variation in $infile: $pinvar > $cutoff");
		}
	}
	return @files;
}

# realigns a phylip file
sub handle_file {
	my ($infile,$reference) = @_;

	# parse the file
	$Log->info("going to read $infile");
	my $matrix = parse_matrix(
		'-format'     => 'phylip',
		'-type'       => 'dna',
		'-as_project' => 1,
		'-file'       => $infile,
	);		
	
	# fetch the refseq, do the protaln if possible
	if ( my $refseq = $matrix->get_by_name($reference) ) {
		my $results = blast_refseq($refseq);
		my ($acc,$cds,$translation) = evaluate_results($results);
		my $realigned  = realign(
			'file'   => $infile,
			'matrix' => $matrix,
			'cds'    => $cds,
			'acc'    => $acc,
		);	
		protalign($realigned,$acc,$translation);
		write_outfile($infile,$realigned);	
	}
	else {
		$Log->warn("no refseq in $infile");
	}	
}

sub write_outfile {
	my ($infile,$m) = @_;

	# make new outfile name
	my $outfile = $infile;
	$outfile =~ s/[a-z]+$/realign.phy/;
	$Log->info("writing result to $outfile");
	open my $outfh, '>', $outfile or die $!;
	print $outfh $m->get_ntax, ' ', $m->get_nchar, "\n";
	$m->visit(sub{
		my $row  = shift;
		my $char = $row->get_char;
		my $name = $row->get_name;
		print $outfh $name, ' ', $char, "\n";
	});
}

# runs BLAST on reference sequence, returns results set
sub blast_refseq {
	my $refseq = shift;
	$Log->info("going to run BLAST on $refseq");
	$Log->debug(my $char = $refseq->get_char);

	# instantiate blast client			
	my $blast = Bio::Tools::Run::RemoteBlast->new(
		'-data'       => 'nr',
		'-prog'       => 'blastn',
		'-expect'     => '1e-10',
		'-readmethod' => 'SearchIO',
	);			
	my $result = $blast->submit_blast($refseq);
	
	# iterate over result ID sets
	while( my @rids = $blast->each_rid ) {
	
		# iterate over result IDs
		for my $rid ( @rids ) {
		
			# fetch results
			my $results = $blast->retrieve_blast($rid);
			
			# not an object
			if ( not ref $results ) {
				if ( $results < 0 ) {
					$blast->remove_rid($rid);
				}
				$Log->info("waiting...");
				sleep 5;
			}
			
			# results is an object
			else {
				return $results;
			}
		}
	}
}

# iterates over BLAST results, returns best hit that has a CDS that intersects
# with our query
sub evaluate_results {
	my $results = shift;
	$Log->info("going to evaluate BLAST results $results");
	while( my $result = $results->next_result ) {
	
		# this is a genbank record
		while( my $hit = $result->next_hit ) {
		
			# fetch the record
			my $acc = $hit->accession;
			$Log->info("going to fetch $acc");
			my $seq = Bio::DB::GenBank->new->get_Seq_by_acc($acc);
			
			# now we have to do the following: 
			# 1. fetch the CDS parts of the genbank record, if any
			# 2. see if they intersect with any of the HSPs
			# if there are CDS parts, one (or more?) of them must
			# intersect, as the refseq is ostensibly from a CDS
			# as well
			# 3. fetch the protein translation 
			# 4. align the protein with our local refseq, (using
			# exonerate?)
			
			for my $feat ( $seq->get_SeqFeatures ) {
				if ( $feat->primary_tag eq 'CDS' ) {
					my $start  = $feat->start;
					my $end    = $feat->end;
					my $strand = $feat->strand;
					
					# now see if the HSPs fall within this range
					while( my $hsp = $hit->next_hsp ) {
						my $hstart  = $hsp->start('hit');
						my $hend    = $hsp->end('hit');
						my $hstrand = $hsp->strand('hit');
						
						# have to be on the same strand for this
						# to make any sense
						if ( $strand == $hstrand ) {
						
							# found the intersection
							if ( $hstart >= $start && $hend <= $end ) {
								my $cds = $feat->spliced_seq->seq;
								my ($translation) = $feat->get_tag_values('translation');
								$Log->debug("spliced seq: $cds");
								$Log->debug("translation: $translation");								
								return $acc, $cds, $translation;
							}
							else {
								$Log->debug("feature $feat doesn't intersect with $hsp");
							}
						}
						else {
							$Log->debug("feature $feat not on same strand as $hsp");
						}
					}
				}								
			}
		}						
	}	
}

# runs muscle on matrix + raw sequence, returns character state matrix
sub realign {
	my %args = @_;
		
	# align sequences
	my $aln = Bio::Tools::Run::Alignment::Muscle->new->align([
		Bio::Phylo::Factory->new->create_datum(
			'-type' => 'dna',
			'-name' => $args{'acc'},
			'-char' => $args{'cds'},
		),
		@{ $args{'matrix'}->get_entities },		
	]);
	
	# convert result to Bio::Phylo character state matrix
	return Bio::Phylo::Matrices::Matrix->new_from_bioperl($aln);
}

# reconciles matrix with protein translation of BLAST hit, uses exonerate
sub protalign {

=begin pod

	my ($matrix,$reference,$translation) = @_;
	$Log->debug("input matrix: \n".$matrix->to_nexus);
	my $datum = $matrix->get_by_name($reference);
	my $char  = $datum->get_char;
	my $unaligned = $char;
	$unaligned =~ s/[\-\?]//g;
	
	# make temp files
	my ( $dnaFH, $dnaFile ) = tempfile();
	print $dnaFH ">$reference\n$unaligned\n";
	my ( $protFH, $protFile ) = tempfile();
	print $protFH ">$reference\n$translation\n";

	# run exonerate, capture result
	my $result = `exonerate -t $dnaFile -T dna -q $protFile -Q protein -m protein2dna`;
	unlink( $dnaFile, $protFile );
	$Log->debug("\n$result");
	
	# parse result
	my @protaln;
	my ( $start, $stop );	
	{
		my ( $query, $target );
		for my $line ( split /\n/, $result ) {
			if ( $line =~ /Target range: (\d+) -> (\d+)/ ) {
				( $start, $stop ) = ( $1, $2 );
			}
			elsif ( $line =~ /^\s+\d+\s+:\s+(\S+)/ and not $query ) {
				$query = $1;
			}
			elsif ( $line =~ /^\s+\d+\s+:\s+(\S+)/ and not $target ) {
				$target = $1;
				push @protaln, grep { /\S/ } split //, $target;
				undef $query;
				undef $target;
			}	
		}
	}

	# first calculate which parts to keep. here we need to translate between the
	# site numbers in the unaligned sequence (which exonerate reports) and column
	# numbers in the original alignment
	my @aligned = split //, $char;
	my $unaligned_base_idx;
	my @keep;
	for my $i ( 0 .. $#aligned ) {
	
		# count which base we've seen
		if ( $aligned[$i] ne '-' && $aligned[$i] ne '?' ) {
			$unaligned_base_idx++;		
		}
		
		# we are within the bases of interest
		if ( $unaligned_base_idx >= $start && $unaligned_base_idx <= $stop ) {
			push @keep, $i;
		}
	}
	
	# make a trimmed version of the original alignment
	$Log->debug("columns to keep:\n@keep");
	my %raw;
	$matrix->visit(sub{
		my $row = shift;
		my $name = $row->get_name;
		my @char = $row->get_char;
		my @sub = @char[@keep];
		$raw{$name} = \@sub;
	});
	
	# initialize a raw new alignment
	my %reconciled = map { $_ => [] } keys %raw;
	
	# iterate over the protein alignment
	my $j = 0;
	for my $i ( 0 .. $#protaln ) {
	
		# the protein alignment has inserted a gap, we must insert one also
		if ( $protaln[$i] eq '-' ) {
			push @{ $reconciled{$_} }, '-' for keys %raw;
		}
		
		# no gap in protein alignment, need to copy over what we have in the alignment.
		# if there are gaps in our alignment, keep copying
		else {
			do {
				push @{ $reconciled{$_} }, $raw{$_}->[$j] for keys %raw;
				$j++;
			
			} while exists $raw{$reference}->[$j] and $raw{$reference}->[$j] eq '-';
		}
	}
	
	# transform to AoA
	my @raw = map { [ $_, join('', grep { defined } @{ $reconciled{$_} }) ] } keys %reconciled;
	$matrix->set_raw(\@raw);
	$Log->debug("output matrix: \n".$matrix->to_nexus);	

=cut

}
	
sub main {
	my ($divergers,$reference,$cutoff) = check_args();
	for my $infile ( read_spreadsheet($divergers,$cutoff) ) {
		handle_file($infile,$reference);
	}
}
main();
