#!/usr/bin/perl
use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Data::Dumper;

=head1 TITLE

gene_extractor.pl - extracts genes from alignments

=head1 DESCRIPTION

This script extracts loci from aligned BAM/SAM files. The loci are provided on the command 
line as identifiers that occur in a provided GFF3 annotation file, from which this script 
extracts the coordinates (and strandedness). Subsequently, for each identified region this 
script will create a consensus sequence for each provided BAM/SAM file. These consensus 
sequences are then written to FASTA files, one file per locus.

=head1 SYNOPSIS

 gene_extractor.pl -id <gene id> -gff3 <gff3 file> -bam <bam file> -ref <refseq>

Note that the chromosome identifiers in the GFF3 file, the reference genome and the
BAM files must match exactly.

=head1 ARGUMENTS

All arguments can both be provided in long form (--argument=value) and
in the shortest unambiguous form (-a v).

=over

=item B<--id=gene id>

A gene identifier as it occurs in the GFF3 annotation file. This argument
can be used multiple times. Each occurrence will result in a separate FASTA
files with the consensus sequences from all BAM files for that locus.

=item B<--gff3=gff3 annotation>

An annotation file in GFF3 format.

=item B<--bam=bam file>

BAM/SAM files from which the loci are to be extracted. This argument can be
used multiple times.

=item B<--refseq=FASTA file>

The location of the reference genome in FASTA format.

=back

=head1 OPTIONS

The following arguments are optional.

=over

=item B<--samtools=location of samtools>

Location of the samtools executable. This is optional, by default the program
is assumed to be on the path.

=item B<--bcftools=location of bcftools>

Location of the bcftools executable. This is optional, by default the program
is assumed to be on the path.

=item B<--vcfutils=location of vcfutils.pl>

Location of the bcftools vcfutils.pl. This is optional, by default the program
is assumed to be on the path.

=item B<--workdir=working directory>

Directory to write BED and FASTA files to. The default is the current working directory.

=item B<--help|?>

Prints this help message and quits.

=item B<--verbose>

Increments the verbosity level.

=back

=head1 DEPENDENCIES

This script depends on the following programs.

=over

=item C<samtools>

Required is a version that has the mpileup command (not pileup) whose output can be
piped into bcftools (i.e. a version of C<samtools> >r865). Samtools includes the programs
C<bcftools> and C<vcfutils.pl>. By default this script assumes these three executables 
are on the path, otherwise their locations need to be provided on the command line.

For more about the tools that this script integrates, visit this web page:
L<http://samtools.sourceforge.net/mpileup.shtml>

=back

=head1 OUTPUT FILES

=over

=item BED files

For each locus identifier supplied on the command line, a simple BED file will be created
in the working directory. The BED file will have the name of the identifier with the 
'.bed' suffix.

=item FASTA files

For each locus identifier supplied on the command line, a FASTA file will be created that
contains the consensus sequences for that locus as extracted from all provided BAM files.
The definition line of each FASTA record will consist of the locus ID, the chromosome and
the start and stop coordinates all separated by spaces.

=back

=head1 ABOUT

This script was first written and released by Rutger Vos for the Naturalis Biodiversity 
Center's bio-informatics programme and our consortium partners. It is open source under 
a CC0 license.

=cut

# command line arguments, some with defaults
my @ids;                      # array of gene IDs of interest
my $gff3;                     # GFF3 annotation file
my @bams;                     # array of BAM files from whence to extract
my $samtools = 'samtools';    # default assumption is that it's on PATH
my $bcftools = 'bcftools';    # default assumption is that it's on PATH
my $vcfutils = 'vcfutils.pl'; # default assumption is that it's on PATH
my $help;                     # if true, print help message and quit
my $verbosity = 1;            # logging level
my $refseq;                   # the reference genome in FASTA format
my $workdir = '.';            # working directory to write output files to

# basic logging functionality
sub LOG ($$) {
	my ($msg,$method) = @_;
	my ( $package, $file1up, $line1up, $sub ) = caller( 2 );
	my ( $pack0up, $file, $line, $sub0up )    = caller( 1 );
	my $log = sprintf( "%s %s [%s, %s] - %s\n", uc $method, $sub || '', $file, $line, $msg );
	print STDERR $log;
}
sub DEBUG ($) { LOG shift, 'DEBUG' if $verbosity >= 3 }
sub INFO ($)  { LOG shift, 'INFO'  if $verbosity >= 2 }
sub WARN ($)  { LOG shift, 'WARN'  if $verbosity >= 1 }

# process command line arguments (void)
sub check_args {
	GetOptions(
		'id=s'       => \@ids,
		'gff3=s'     => \$gff3,
		'bam=s'      => \@bams,
		'samtools=s' => \$samtools,
		'bcftools=s' => \$bcftools,
		'vcfutils=s' => \$vcfutils,
		'verbose+'   => \$verbosity,
		'help|?'     => \$help,
		'refseq=s'   => \$refseq,
	);

	# print help message and quit, if requested
	$help && pod2usage({ '-verbose' => $help });

	# check sanity
	if ( not @ids or not $gff3 or not @bams ) {
		pod2usage({ '-verbose' => 1 });
	}

	# report provided arguments
	INFO "provided IDs: @ids";
	INFO "GFF3 file: $gff3";
	INFO "BAM files: @bams";
	INFO "samtools: $samtools";
}

# read annotation file, return extracted coordinates
sub read_gff3 {

	# create mapping of requested IDs
	my %ids = map { $_ => {} } @ids;
	INFO "going to read GFF3 annotation file";
	open my $fh, '<', $gff3 or die $!;
	while(<$fh>) {
		chomp;
		next if /^#/; # skip comments and headers
		
		# split the entire record
		my ( $chromo, $auth, $tag, $start, $stop, $dot1, $strand, $dot2, $desc ) = split;
		
		# split the description column
		my %meta = map { split /=/, $_ } split /;/, $desc;
		
		# we want this
		if ( $ids{$meta{'ID'}} ) {
			my $id = $meta{'ID'};
			$ids{$id}->{'start'}  = $start;
			$ids{$id}->{'stop'}   = $stop;
			$ids{$id}->{'strand'} = $strand;
			$ids{$id}->{'chromo'} = $chromo;
			INFO "found $id ($start..$stop) on the $strand strand";
		}
		else {
			DEBUG "skipping feature $meta{ID}";
		}
	}
	return %ids;
}

# write consensus sequences to FASTA files, receives extracted coordinates
sub write_consensus {
	my %args = @_;
	INFO "going to write consensus sequences to FASTA files";
	
	for my $id ( keys %args ) {
		my %region = %{ $args{$id} };
	
		# we will now create a one-line BED file
		my $bedline = join ' ', @region{qw[chromo start stop]};
		my $bedfile = "${workdir}/${id}.bed";
		open my $bedFH, '>', $bedfile or die $!;
		print $bedFH $bedline, "\n";
		INFO "wrote $bedfile containing '$bedline'";
		
		# instantiate the new FASTA file with consensus sequences
		my $fasta = "${workdir}/${id}.fas";
		open my $fastaFH, '>', $fasta or die $!;
		INFO "instantiated FASTA file $fasta";
		
		# iterate over the BAM files
		for my $bam ( @bams ) {
		
			# here we will open a handle to a pipe from which we read the FASTQ result
			my $pipe    = "| $bcftools view -cg - | $vcfutils vcf2fq";
			my $command = "$samtools mpileup -uf $refseq -l $bedfile $bam $pipe";
			my $output  = `$command`;
			DEBUG "ran command '$command'";
			
			# now extract the subsequence from the FASTQ and write to FASTA
			my $subseq = read_fastq( 'lines' => [ split /\n/, $output ], %region );
			print $fastaFH ">$bam $id $bedline\n$subseq\n";
			DEBUG "printed subseq to FASTA file: $subseq";
		}
		INFO "populated FASTA file $fasta";
	}
}

# reads from FASTQ handle, returns seq from provided coordinates
sub read_fastq {
	my %args = @_; # e.g. fh => foo, chromo => bar, start => baz, stop => quux
	INFO "going to read FASTQ data";
	
	my ( $seqlength, $phredlength ) = ( 0, 0 );
	my ( $id, $seq, $plus );
	RECORD: for my $line ( @{ $args{'lines'} } ) {
		chomp $line;
		
		# find the FASTQ id line
		if ( not $id and $line =~ /^\@(.+)$/ ) {
			$id = $1;
			INFO "found record ID $id";
		}
		
		# concatenate the sequences, output is multiline!
		elsif ( $id and not $plus and $line !~ m/^\+/ ) {
			$seq .= $line;
			$seqlength++;
		}
		
		# check for plus line
		elsif ( $id and $seq and not $plus and $line =~  m/^\+/ ) {
			$plus++;
			INFO "found + line";
		}
		
		# check for end of phred lines, output is multiline!
		elsif ( $id and $seq and $plus ) {
			if ( $phredlength < ( $seqlength - 1 ) ) {
				$phredlength++;
			}
			else {
			
				# return wanted result or continue to next chromosome
				if ( $id eq $args{'chromo'} ) {
					INFO "found chromosome of interest $id";
					
					# sequence is from 5' -->
					if ( $args{'strand'} eq '+' ) {
						INFO "locus is on + strand";
						return substr $seq, $args{'start'}, $args{'stop'} - $args{'start'};
					}
					
					# sequence is <-- 3'
					else {
						INFO "locus is on - strand, will reverse complement";
						my $sub = substr reverse($seq), $args{'start'}, $args{'stop'} - $args{'start'};
						$sub =~ tr/ACGTacgt/TGCAtgca/;
						return $sub;
					}
				}
				else {
					INFO "ignoring: chromosome $id != $args{chromo}";
					undef $id;
					undef $seq;
					undef $plus;
					( $seqlength, $phredlength ) = ( 0, 0 );
					next RECORD;
				}
			}
		}
	}
}

# to aid comprehension, this script is organized in a main that calls the other subs. 
# start here and drill down if you're trying to make sense of this.
sub main {
	check_args();
	my %coordinates = read_gff3();
	write_consensus(%coordinates);
}
main();
