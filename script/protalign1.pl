#!/usr/bin/perl
use strict;
use warnings;
use File::Temp 'tempfile';
use Getopt::Long;
use Bio::DB::GenBank;
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
	my $infile = shift;
	INFO "going to read FASTA data from $infile";
	open my $fh, '<', $infile or die $!;
	my ( $current, %data );
	while(<$fh>) {
		chomp;
		if ( />(\S+)/ ) {
			$current = $1;
			$data{$current} = [];
			INFO "reading sequence for $current";
		}
		else {
			my @seq = split //, $_;
			push @{ $data{$current} }, @seq;
		}
	}
	my ($length) = sort { $b <=> $a } map { scalar @{ $_ } } values %data;
	return \%data, $length;
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

sub run_exonerate {
	my ( $reference, $seq, $translation ) = @_;
	
	# unalign the nucleotide sequence
	my $unaligned = $seq;
	$unaligned =~ s/[\-\?N]//g;
	
	# make temp files
	my ( $dnaFH, $dnaFile ) = tempfile();
	print $dnaFH ">$reference\n$unaligned\n";
	my ( $protFH, $protFile ) = tempfile();
	print $protFH ">$reference\n$translation\n";
	
	# run the analysis
	INFO "going to run exonerate on DNA $dnaFile Protein $protFile";
	my $result = `exonerate -t $dnaFile -T dna -q $protFile -Q protein -m protein2dna`;
	unlink( $dnaFile, $protFile );

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
	
	# return results
	INFO "start: $start";
	DEBUG $result;
	return $start, $stop, \@protaln;	
}

sub reconcile_alignments {
	my ( $data, $start, $stop, $refseq, $protaln, $length ) = @_;
	INFO "reconciling MSA and protein/$refseq alignments between $start and $stop";
	my @ids    = keys %{ $data };
	my %result = map { $_ => [] } @ids;
	my $pgap   = 0;
	my $i = 0;
	my $uidx = 0;	
	BASE: while(1) {
		last BASE if not exists $data->{$refseq}->[$i];
			
		# track the unaligned index
		$uidx++ if $data->{$refseq}->[$i] =~ /[ACGT]/i;
#		INFO $uidx;
		if ( $uidx > $start and ( $uidx + $pgap ) <= $stop ) {
#			INFO "$uidx >= $start and ( $uidx + $pgap ) <= $stop";

			# three possibilities: 
			# 1. both protaln and nucaln are bases, copy one over to other
			# 2. there is a gap in the protaln, copy over
			# 3. there is a gap in the nucaln, pad to make it a multiple of three
			
			# there is a gap in the exonerate alignment, we need to add that to all the sequences
			if ( $protaln->[$uidx-$start] !~ /[ACGT]/i ) {
				$pgap++;
				INFO "inserting protein alignment gap at $uidx-$start";
				push @{ $result{$_} }, '-', for @ids;
				next BASE;
			}
			
			# there is a gap in the nucleotide alignment
			if ( $data->{$refseq}->[$i] !~ /[ACGT]/i ) {
				INFO "inserting MSA gap at $i";
				
				# copy the gap over
				my $j = $i;
				GAP: for $j ( $i .. $length - 1 ) {
					last GAP if $data->{$refseq}->[$i] =~ /[ACGT]/i;
					push @{ $result{$_} }, $data->{$_}->[$j] for @ids;
				}
				my $diff = $j - $i + 1;
				INFO "gap was $diff base pairs";
				
				# check for frame shifts
				if ( my $shift = $diff % 3 ) {
					INFO "gap induces a $shift bp shift, padding";
					for ( 1 .. $shift ) {
						push @{ $result{$_} }, '-', for @ids;
					}
				}
				$i = $j + 1;
				next BASE;
			}
			
			# no gap in either, copy verbatim
			push @{ $result{$_} }, $data->{$_}->[$i], for @ids;
		}
		$i++;
#		INFO $i++;		
	}	
	return \%result;
}

sub write_fasta {
	my $data = shift;
	INFO "going to write result to STDOUT";
	for my $key ( keys %{ $data } ) {
		my $seq = join '', @{ $data->{$key} };
		print '>', $key, "\n", $seq, "\n";
	}
}

sub main {
	my %args = check_args();
	
	# read the fasta file
	my $file = $args{'infile'};
	my ( $data, $length ) = read_fasta($file);
	
	# fetch the protein translation
	my $refseq = $args{'refseq'};
	my $translation = fetch_protein($refseq);
	
	# run exonerate on the reference sequence
	my $seq = join '', @{ $data->{$refseq} };
	my ( $start, $stop, $protaln ) = run_exonerate( $refseq, $seq, $translation );
	
	# reconcile the alignments
	my $reconciled = reconcile_alignments($data,$start,$stop,$refseq,$protaln,$length);

	# write the output
	write_fasta($reconciled);
}

main();
