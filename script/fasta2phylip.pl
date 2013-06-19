#!/usr/bin/perl
use strict;
use warnings;
use Bio::Tools::CodonTable;
my %fasta;
my $current;

while(<>){
	chomp;
	if ( />(\S+)/ ) {
		my $name = $1;
		$current = $name;
	}
	else {
		s/!/-/g;
		$fasta{$current} = [] unless $fasta{$current};
		push @{ $fasta{$current} }, split //;
	}
}

my ($length) = map { scalar @{ $_ } } values %fasta;
my %stripped = map { $_ => '' } keys %fasta;
my $table = Bio::Tools::CodonTable->new;

for ( my $i = 0; $i < $length; $i += 3 ) {
	my $has_ter;
	for my $row ( keys %fasta ) {
		my @row = @{ $fasta{$row} };
		my $codon = join '', @row[$i,$i+1,$i+2];
		$has_ter++ if $table->is_ter_codon($codon);
	}
	if ( not $has_ter ) {
		for my $row ( keys %fasta ) {
			my @row = @{ $fasta{$row} };
			my $codon = join '', @row[$i,$i+1,$i+2];
			$stripped{$row} .= $codon;
		}		
	}
}

my ($nchar) = map { length($_) } values %stripped;
print scalar(keys(%stripped)), ' ', $nchar, "\n";
print $_, ' ', $stripped{$_}, "\n" for keys %stripped;