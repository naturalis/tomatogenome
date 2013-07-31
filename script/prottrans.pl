#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::DB::GenBank;
use File::Temp 'tempfile';
use Bio::Phylo::Util::Logger ':simple';

my %aasym = (
	'Ala' => 'A', # alanine
	'Arg' => 'R', # arginine
	'Asn' => 'N', # asparagine
	'Asp' => 'D', # aspartic acid
	'Asx' => 'B', # asparagine or aspartic acid
	'Cys' => 'C', # cysteine
	'Glu' => 'E', # glutamic acid
	'Gln' => 'Q', # glutamine
	'Glx' => 'Z', # glutamine or glutamic acid
	'Gly' => 'G', # glycine
	'His' => 'H', # histidine
	'Ile' => 'I', # isoleucine
	'Leu' => 'L', # leucine
	'Lys' => 'K', # lysine
	'Met' => 'M', # methionine
	'Phe' => 'F', # phenylalanine
	'Pro' => 'P', # proline
	'Ser' => 'S', # serine
	'Thr' => 'T', # threonine
	'Trp' => 'W', # tryptophan
	'Tyr' => 'Y', # tyrosine
	'Val' => 'V', # valine
	'Unk' => 'X', # unknown
);

sub check_args {
	 my $verbosity = WARN;
	 my ( $infile, $refseq );
	 GetOptions(
	 	'verbose+' => \$verbosity,
	 	'infile=s'  => \$infile,
	 	'refseq=s' => \$refseq,
	 );
	 VERBOSE( '-level' => $verbosity, '-class' => 'main' );
	 return $infile, $refseq;
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

sub read_fasta {
	my $fh = shift;
	INFO "going to read FASTA data";
	my %result;
	my $current;
	while(<$fh>) {
		chomp;
		if ( />(\S+)/ ) {
			$current = $1;
			INFO "going to read sequence data for $current";
		}
		else {
			$result{$current} .= $_;
		}
	}	
	return \%result;
}

sub run_exonerate {
	my ( $id, $unaligned, $translation ) = @_;
	
	# make temp files
	my ( $dnaFH, $dnaFile ) = tempfile();
	write_fasta( { $id => $unaligned }, $dnaFH );
	my ( $protFH, $protFile ) = tempfile();
	write_fasta( { $id => $translation }, $protFH );
	
	# run the analysis
	INFO "going to run exonerate for $id on DNA $dnaFile Protein $protFile";
	my $result = `exonerate -t $dnaFile -T dna -q $protFile -Q protein -m protein2dna -E`;
	DEBUG $result;
	unlink( $dnaFile, $protFile );

	# parse result
	my $aa;
	my $nuc;
	my $seen_aa;
	my $line_count = 0;
	LINE: for ( split /\n/, $result ) {
		if ( /Target range/ .. /vulgar/ ) {
			chomp;
			if ( /^\s*\d+\s+:\s+\S+\s+:\s+\d+\s*$/ and not $seen_aa ) {
				$line_count = 0;
				$seen_aa = 1;
				next LINE;
			}
			$line_count++ if $seen_aa;
			if ( $line_count == 2 ) {
				s/ //g;
				$aa .= $_;
				DEBUG "AA  $_";
			}
			if ( $line_count == 3 and /\S/ and not /vulgar/ and not /Target/ ) {
				s/[0-9: ]//g;
				$nuc .= $_;
				$seen_aa = 0;
				DEBUG "DNA $_";
			}			
		}
	}
	DEBUG $aa;
	DEBUG $nuc;
	
	# translate three-letter exonerate amino acids to single character IUPAC
	my $cleaned_aa;
	my $cleaned_nuc;
	while( $aa and $nuc ) {
		if ( $aa =~ /^---/ ) {
			$cleaned_aa  .= '-';
			$cleaned_nuc .= '---';
			$aa  = substr $aa, 3;
			$nuc = substr $nuc, 3;
		}
		elsif ( $aa =~ /^\*\*\*/ ) {
			$cleaned_aa  .= '*';
			$cleaned_nuc .= substr $nuc, 0, 3;
			$aa  = substr $aa, 3;
			$nuc = substr $nuc, 3;			
		}
		elsif ( $aa =~ /^#/ ) { # XXX!
			$aa  = substr $aa, 1;
			$nuc = substr $nuc, 1;
		}
		elsif ( my $sym = $aasym{ substr($aa, 0, 3) } ) {
			$cleaned_aa  .= $sym;
			$cleaned_nuc .= substr $nuc, 0, 3;
			$aa = substr $aa, 3;
			$nuc = substr $nuc, 3;			
		}
		else {
			die $aa;
		}
	}

	# done
	INFO $translation;
	INFO $cleaned_aa;
	INFO $cleaned_nuc;

	return $cleaned_aa, $cleaned_nuc;
}

sub write_fasta {
	my ( $data, $fh ) = @_;
	for my $key ( keys %{ $data } ) {
		print $fh '>', $key, "\n", $data->{$key}, "\n";
	}
}

sub run_muscle {
	my $data = shift;
	my ( $protFH, $protFile ) = tempfile();
	write_fasta( $data, $protFH );
	my $result = `muscle -quiet -seqtype protein -in $protFile`;
	unlink $protFile;
	return $result;	
}

sub reconcile_nucprot {
	my ( $nuc, $prot ) = @_;
	my %result;
	my $i = 0;
	no warnings;
	AA: while(1) {
		for my $key ( keys %{ $prot } ) {
			last AA if $i >= length $prot->{$key};
			my $aa = substr $prot->{$key}, $i, 1;
			if ( '-' eq $aa ) {
				$result{$key} .= '---';
			}
			else {
				my $codon = substr $nuc->{$key}, $i * 3, 3;
				$result{$key} .= ( $codon || '---' );
			}
		}
		$i++;
	}
	use warnings;
	return \%result;
}

sub main {

	# process command line arguments
	my ( $infile, $refseq ) = check_args;
	INFO "input file is $infile, reference sequence is $refseq";
	
	# download the reference protein from genbank
	my $reftrans = fetch_protein($refseq);
	
	# read the sequences
	my $nucdata;
	{
		open my $fh, '<', $infile or die $!;
		$nucdata = read_fasta( $fh );
	}
	
	# iterate over sequences
	my ( %protein, %dna );
	for my $key ( keys %{ $nucdata } ) {
		my $seq = $nucdata->{$key};
		$seq =~ s/[N\-\?]//g;
		my ( $aa, $nuc ) = run_exonerate( $key, $seq, $reftrans );
		$protein{$key} = $aa;
		$dna{$key} = $nuc; # this is now aligned relative to its protein
	}

	# align proteins
	my $aaaln = run_muscle( \%protein );
	
	# read the protein alignment
	my $protdata;
	{
		open my $fh, '<', \$aaaln or die $!;
		$protdata = read_fasta( $fh );
	}
	
	# reconcile the nucleotides (which are aligned against their respective protein
	# translations) with the amino acid alignment
	my $reconciled = reconcile_nucprot( \%dna, $protdata );
	
	# done
	write_fasta( $reconciled, \*STDOUT );
}
main();