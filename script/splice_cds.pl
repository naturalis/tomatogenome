#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':simple';

sub check_args {
	my $verbosity = WARN;
	my ( $fasta, $gff3, $id, $refseq );
	GetOptions(
		'verbose+' => \$verbosity,
		'fasta=s'  => \$fasta,
		'gff3=s'   => \$gff3,
		'id=s'     => \$id,
		'refseq=s' => \$refseq,	
	);
	VERBOSE( '-level' => $verbosity, '-class' => 'main' );
	INFO "FASTA file is $fasta";
	INFO "Gene ID is $id";
	INFO "GFF3 is $gff3";
	return 
		'fasta'  => $fasta,
		'gff3'   => $gff3,
		'id'     => $id,
		'refseq' => $refseq;
}

sub read_fasta {
	my $fasta = shift;
	
	# start reading
	INFO "going to read FASTA file $fasta";
	my %data;
	my $current;
	open my $fh, '<', $fasta or die $!;
	LINE: while(<$fh>) {
		chomp;
		if ( /^>(\S+)/ ) {
			$current = $1;
			$data{$current} = [];
			INFO "going to read sequence for $current";
			next LINE;
		}
		my @seq = split //, $_;
		push @{ $data{$current} }, @seq;	
	}
	return \%data;
}

sub write_fasta {
	my $data = shift;
	INFO "going to write FASTA data to STDOUT";
	for my $id ( keys %{ $data } ) {
		print ">$id\n", @{ $data->{$id} }, "\n";
	}
}

sub splice_alignment {
	my ( $data, $refseq, @cds ) = @_;
	my @ids = keys %{ $data };
	my %spliced = map { $_ => [] } @ids;
	my $i = 0;
	CDS: for my $cds ( @cds ) {
		INFO "splicing CDS " . $cds->{cds};
	
		# insert gaps as needed to stay in frame
		for ( 1 .. $cds->{offset} ) {
			push @{ $spliced{$_} }, '-' for @ids;
		}
		
		# copy $cds->{length} columns that aren't gaps in $refseq
		my $length = 0;
		while( $length < $cds->{length} ) {
			if ( exists $data->{$refseq}->[$i] and $data->{$refseq}->[$i] ne '-' ) {
				$length++;
				push @{ $spliced{$_} }, $data->{$_}->[$i] for @ids;
			}
			last CDS if not exists $data->{$refseq}->[$i];			
			$i++;
		}
	}	
	return \%spliced;
}

sub fetch_cds {
	my ( $gff3, $id ) = @_;

	# gff3 column numbers
	my $type_idx   = 2;
	my $start_idx  = 3;
	my $end_idx    = 4;
	my $strand_idx = 6;
	my $offset_idx = 7;
	my $meta_idx   = 8;
	
	# start reading
	INFO "going to scan $gff3 for CDSs belonging to gene $id";
	my $mRNA;
	my $strand;
	my @result;
	open my $fh, '<', $gff3 or die $!;
	LINE: while(<$fh>) {
		chomp;
		next LINE if /^#/; # skip header
		my @line = split /\t/, $_;
		my $meta = $line[$meta_idx];
		
		# check if we've found our gene as parent of focal mRNA
		if ( $meta =~ /Parent=gene:$id/ and $meta =~ /ID=mRNA:([^;]+)/ ) {
			$mRNA = $1;
			INFO "found mRNA $mRNA for gene $id";
			next LINE;
		}
		
		# either we're reading CDSs for the focal mRNA or we've passed it,
		# in which case we stop scanning
		if ( $mRNA ) {
			if ( $line[$type_idx] eq 'CDS' ) {
				if ( $meta =~ /Parent=mRNA:$mRNA/ and $meta =~ /ID=CDS:([^;]+)/ ) {
					my $cds = $1;					
					push @result, {
						'start'  => $line[$start_idx],
						'end'    => $line[$end_idx],
						'strand' => $line[$strand_idx],
						'offset' => $line[$offset_idx],
						'length' => $line[$end_idx] - $line[$start_idx],
						'cds'    => $cds,
					};					
					if ( $strand and $strand ne $line[$strand_idx] ) {
						die "strand is flipped halfway in?";
					}
					$strand = $line[$strand_idx];
					INFO "found CDS $cds on strand $strand";
				}
				else {
					INFO "encountered next mRNA, stop scanning";
					last LINE;
				}
			}		
		}
	}
	return $strand eq '-' ? reverse @result : @result;
}

sub main {
	my %args = check_args();
	my @cds  = fetch_cds( @args{'gff3','id'} );
	my $data = read_fasta( $args{'fasta'} );
	my $spliced = splice_alignment( $data, $args{'refseq'}, @cds );
	write_fasta( $spliced );
}
main();