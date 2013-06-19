#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::Sam;
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my $window = 1000;
my $readlength = 101;
my ( $bam, $fasta );
GetOptions(
	'window=i'     => \$window,
	'verbose+'     => \$verbosity,
	'bam=s'        => \$bam,
	'fasta=s'      => \$fasta,
	'readlength=i' => \$readlength,
);

unless ( -e $bam && -e $fasta ) {
	die "Usage: $0 -bam <aln.bam> -fasta <ref.fa> [-window 1000] [-verbose] > <outfile>\n";
}

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);
my $sam = Bio::DB::Sam->new(
	'-bam'   => $bam,
	'-fasta' => $fasta,
);

# get the identifiers of the records
# in the FASTA reference file. These
# should be the chromosomes.
my @references = $sam->seq_ids;

# iterate over the reference sequences
my $totalreads = 0;
for my $ref ( @references ) {
	my $length = $sam->length($ref);
        $log->info("analyzing reference $ref ($length bp)");

	# now jump in windows, locations appear to be 1-based
	my $chrreads = 0;
	for ( my $i = 1; $i <= $length; $i += $window ) {
		$log->debug("window: $i");

		# don't know if there can be weird "index out of bounds"
		# errors, so let's be tidy and don't exceed length
		my $end = $i + $window > $length ? $length : $i + $window;
		my @aln = grep { $_->start >= $i } $sam->get_features_by_location(
			'-seq_id' => $ref,
			'-start'  => $i,
			'-end'    => $end,
		);
		my $cover = ( $end - $i ) / ( $readlength * ( scalar(@aln) || 1 ) );
		printf "%s %i %i %f\n", $ref, $i, $end, $cover;
		$chrreads += scalar @aln;
	}
	$log->info("seen $chrreads reads mapped on $ref");
	$totalreads += $chrreads;
}
$log->info("seen $totalreads mapped reads in $bam");
