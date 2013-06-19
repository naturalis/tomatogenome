#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::DB::Sam;
use Bio::Tools::IUPAC;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my $name = ''; # optional, a name prefix 
my ( $queryfile, $dbfile, $bamfile, $help, $printquery );
GetOptions(
	'queryfile=s' => \$queryfile,
	'dbfile=s'    => \$dbfile,
	'bamfile=s'   => \$bamfile,
	'verbose+'    => \$verbosity,
	'help'        => \$help,
	'name=s'      => \$name,
	'printquery'  => \$printquery,
);

# too complicated to use, need help message:
if ( $help ) {
	print <<"HELP";
Usage: $0 -query <fasta file> -db <local blast file> -bam <bam file> [-name <defline prefix>] [--printquery]

The query file contains one or more fasta sequences for interesting loci.
These are then blast searched against the reference genome so that we
can reconstruct the coordinates of the loci in the BAM file. Out of the
bam file we then extract a consensus sequence, which is printed to STDOUT.
HELP
exit 0;
}

# when --printquery is set, the query sequence is printed to STDOUT. Since
# we do the same with the consensus sequences we should hereby be able to
# pipe the STDOUT of this script to muscle
if ( $printquery ) {
	open my $fh, '<', $queryfile or die $!;
	while(<$fh>) {
		print $_;
	}
	close $fh;
}

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);
$log->info("instantiating SAM file handler for file $bamfile with reference $dbfile");
my $sam = Bio::DB::Sam->new(
	'-bam'   => $bamfile,
	'-fasta' => $dbfile,
);
$log->info("instantiating standalone BLAST+ with database $dbfile");
my $blast = Bio::Tools::Run::StandAloneBlastPlus->new(
	'-db_data' => $dbfile,
#	'-create'  => 1,
);
$log->info("instantiating file reader for $queryfile in FASTA format");
my $io = Bio::SeqIO->new(
	'-file'   => $queryfile,
	'-format' => 'fasta',
);
$log->info("instantiating IUPAC ambiguity code mapper");
my $iupac = Bio::Tools::IUPAC->new;
my %rev  = $iupac->iupac_rev_iub;
my %comp = ( 'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A' );

# iterate over sequences in query file. the query file contains loci
# from genbank whose coordinates into the reference genome we are 
# fetching at this step
while( my $seq = $io->next_seq ) {
	my $id = $seq->id;
	$log->info("going to BLAST sequence '$id'");

	# run blast
	my $result = $blast->blastn( '-query' => $seq );
	my $resultname = $result->query_name; # should be same as $id
	$log->info("found result for '$resultname'");

	# iterate over hits
	while( my $hit = $result->next_hit ) {
		my $hitname = $hit->name;
		$hitname =~ s/^lcl\|//; # strip prefix
		$log->info("found hit '$hitname'"); # should be chromosome

		# iterate over high scoring pairs
		while( my $hsp = $hit->next_hsp ) {
			my $start = $hsp->start("subject");
			my $end   = $hsp->end("subject");
			$log->info("HSP range is $start..$end (".($end-$start)." bp)");

			# having found the coordinates of the focal
			# hit into the reference we will now
			# query the BAM file for the same coordinates
			getreads(
				'start' => $start,
				'end'   => $end,
				'seq'   => $hitname, 
				'id'    => $id,
				'query' => $hsp->seq_str("query"),
			);
		}
	}
}

sub getreads {
	my %args = @_;
	return unless $args{'seq'} eq 'SL2.40ch07';
	$log->info("going to fetch bases from $args{seq} between $args{start} and $args{end}");

	# get the reads that intersect with this range
	my @reads = $sam->get_features_by_location(
		'-seq_id' => $args{"seq"},
		'-start'  => $args{"start"},
		'-end'    => $args{"end"},
	);

	# print the sites one by one
	my ( @seq, @refseq );
	for my $i ( $args{"start"} .. $args{"end"} ) {
		my ( %read_sites, %ref_sites );
		for my $read ( @reads ) {
			my $start = $read->start;
			my $end = $read->end;
			if ( $i >= $start and $i <= $end ) {
				my $ref  = $read->dna;
				my $dna  = $read->query->dna;				
				my $site = substr $dna, $i - $start, 1;
				my $ref_site = substr $ref, $i - $start, 1;
				$read_sites{$site}++;
				$ref_sites{$ref_site}++;
			}
		}
		push @refseq, keys %ref_sites;
		if ( %read_sites ) {
			push @seq, $rev{ join '', sort { $a cmp $b } keys %read_sites };
		}
	}

	# reverse complement
	my $seq = reverse join '', @seq;
	my $refseq = reverse join '', @refseq;
	$seq =~ tr/ACGT/TGCA/;
	$refseq =~ tr/ACGT/TGCA/;

	# print 
	print '>', $name, '-', $args{"seq"}, '-', $args{"start"}, '..', $args{"end"}, "\n";
	print $seq, "\n";

	$log->info("QUERY: $args{query}");
	$log->info("REF:   $refseq");
	$log->info("READS: $seq");
}
