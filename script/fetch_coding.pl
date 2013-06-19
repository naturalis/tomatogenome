#!/usr/bin/perl
use strict;
use warnings;
use File::Temp 'tempfile';
use Getopt::Long;
use Bio::DB::GenBank;
use Bio::Phylo::Factory;
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
	my $reftaxon; # reference taxon to use for BLAST searches
	my $ignorestrand; # ignore strandedness
	
	GetOptions(
		'verbose+'     => \$verbosity,
		'reftaxon=s'   => \$reftaxon,
		'ignorestrand' => \$ignorestrand,
	);

	# instantiate logger
	$Log = Bio::Phylo::Util::Logger->new(
		'-level' => $verbosity,
		'-class' => 'main',
	);

	# communicate settings
	$Log->info("reference taxon is $reftaxon");
	
	return $reftaxon, $ignorestrand;
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
	my ($results,$ignorestrand) = @_;
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
						if ( $strand == $hstrand or $ignorestrand ) {
						
							# found the intersection
							if ( $hstart >= $start && $hend <= $end ) {
								my $cds = $feat->spliced_seq->seq;
								$Log->debug("spliced seq: $cds");
								return $acc, $cds;
							}
							else {
								$Log->info("feature $feat doesn't intersect with $hsp");
							}
						}
						else {
							$Log->info("feature $feat not on same strand as $hsp");
						}
					}					
				}
				else {
					$Log->info("primary tag: ".$feat->primary_tag);
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

sub main {
	my ($reftaxon,$ignorestrand) = check_args();
	
	# parse the input data
	my $sequences = parse_matrix(
		'-format'     => 'fasta',
		'-type'       => 'dna',
		'-as_project' => 1,
		'-handle'     => \*STDIN,
	);		
	
	# fetch the refseq, do the protaln if possible
	if ( my $refseq = $sequences->get_by_name($reftaxon) ) {
		my $results = blast_refseq($refseq);
		my ($acc,$cds) = evaluate_results($results,$ignorestrand);
		my $realigned = realign(
			'acc'    => $acc,
			'cds'    => $cds,
			'matrix' => $sequences,
		);
		$realigned->visit(sub{
			my $row = shift;
			my $name = $row->get_name;
			my $char = $row->get_char;
			print ">$name\n$char\n";
		});
	}
	else {
		$Log->error("no reference taxon $reftaxon found");
	}
}
main();
