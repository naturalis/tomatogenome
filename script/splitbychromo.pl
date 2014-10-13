#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# $indir = fasta files per individu, $outdir = fasta files per chromosoom
my ( $indir, $outdir );
GetOptions(
	'indir=s'  => \$indir,
	'outdir=s' => \$outdir,
);

opendir my $dh, $indir or die $!;
while( my $entry = readdir $dh ) {
	if ( $entry =~ /\.fas$/ ) {
		open my $fh, '<', "${indir}/${entry}" or die $!;
		my $out;
		LINE: while(<$fh>) {
			if ( /^>(.+)/ ) {
				my $chromo = $1;
				open $out, '>>', "${outdir}/${chromo}.fas" or die $!;
				print $out ">$entry\n";
				next LINE;
			}
			print $out $_;
		}
	}
}