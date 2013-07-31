#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::GenBank;
use Bio::Tools::CodonTable;
use Bio::Phylo::Util::Logger ':simple';

sub check_args {
	my $verbosity = WARN;
	my ( $infile, $refseq );
	GetOptions(
		'verbose+' => \$verbosity,
		'infile=s' => \$infile,
		'refseq=s' => \$refseq,
	);
	VERBOSE( '-level' => $verbosity, '-class' => 'main' );
	return 'infile' => $infile, 'refseq' => $refseq;
}

sub read_fasta {
	my ( $infile, $refseq ) = @_;
	INFO "going to scan fasta file $infile for seq $refseq";
	my $seq;
	open my $fh, '<', $infile or die $!;
	LINE: while(<$fh>) {
		chomp;
		if ( /^>/ ) {
			if ( /$refseq/ ) {
				INFO "found $refseq, going to start reading";
				$seq = '';				
				next LINE;
			}
			if ( $seq ) {
				INFO "encountered next sequence, stop reading";
				last LINE;
			}
		}
		$seq .= $_ if defined $seq;
	}
	return $seq;
}

sub fetch_protein {
	my $refseq = shift;
	INFO "going to fetch protein from GenBank";
	my $seq = Bio::DB::GenBank->new->get_Seq_by_acc($refseq);
	for my $feat ( $seq->get_SeqFeatures ) {
		if ( $feat->primary_tag eq 'CDS' and $feat->has_tag('translation') ) {
			INFO "found translation feature";
			my ($translation) = $feat->get_tag_values('translation');
			return $translation;
		}								
	}
}

sub translate {
	my $seq = shift;
	my $ct  = Bio::Tools::CodonTable->new;
	INFO "translating codons";
	join '', map { $ct->translate($_) } unpack "(A3)*", $seq;
}

sub main {
	my %args = check_args();
	my $seq  = read_fasta(@args{qw[infile refseq]});
	my $prot = fetch_protein($args{'refseq'});
	my $translation = translate($seq);
	$translation =~ s/\-//g;
	print "PROT $prot\n";
	print "TRAN $translation\n";
}
main();