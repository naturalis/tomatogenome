#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':simple';

# process command line arguments
my $verbosity = WARN;
my ( $fasta, $gff3 );
my $type = 'gene';
GetOptions(
	'fasta=s'  => \$fasta,
	'gff3=s'   => \$gff3,
	'verbose+' => \$verbosity,
	'type=s'   => \$type,
);

# communicate settings
VERBOSE( '-level' => $verbosity, '-class' => 'main' );

# open GFF3 handle
INFO "reading annotations from GFF3 file $gff3";
open my $gff3FH, '<', $gff3 or die $!;

# start reading from fasta
INFO "reading genome from FASTA file $fasta";
open my $fastaFH, '<', $fasta or die $!;
my ( $id, $seq );
while(<$fastaFH>) {
	chomp;
	if ( /^>(\S+)$/ ) {
		my $defline = $1;
		if ( $id and $seq ) {
			split_seq( $id, \$seq );
		}
		$id = $defline;
		$seq = '';		
	}
	else {
		$seq .= $_;
	}
}
split_seq( $id, \$seq );

sub split_seq {
	my ( $id, $seqref ) = @_;
	INFO "going to scan for annotations on $id";

	# 0-based column numbers in gff3
	my $chr_idx    = 0;
	my $type_idx   = 2;
	my $start_idx  = 3;
	my $end_idx    = 4;
	my $strand_idx = 6;
	my $meta_idx   = 8;
	
	# iterate over lines
	my $previous_line = 0;
	LINE: while(<$gff3FH>) {
		chomp;
		next LINE if /^#/;
		my @line = split /\t/, $_;
		
		# check to see if we need to backtrack
		if ( $line[$chr_idx] ne $id ) {
			seek $gff3FH, $previous_line, 0;
			INFO "no longer seeing $id, backtracking to byte $previous_line";
			last LINE;
		}
		
		# write the focal annotation
		if ( $line[$type_idx] eq $type ) {
			if ( $line[$meta_idx] =~ /ID=$type:([^;]+)/ ) {
				my $gene   = $1;
				my $offset = $line[$start_idx] - 1;
				my $length = $line[$end_idx] - $line[$start_idx];
				my $seq = substr $$seqref, $offset, $length;
				print '>', $gene, "\n", $seq, "\n";
			}		
		}
		
		# record where we are now, for the next iteration
		$previous_line = tell $gff3FH;
	}
}