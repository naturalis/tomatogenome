#!/usr/bin/perl
use strict;
use warnings;

my ( $defline, @data );
while(<>) {
	chomp;
	if ( not $defline ) {
		$defline = $_;
		next;
	}
	my @line = split //, $_;
	push @data, grep { /\S/ } @line;
}

print $defline, "\n";
while(@data) {
	my $counter = 0;
	while( $counter < 80 ) {
		print shift @data;
		$counter++;
		last if not @data;
	}
	print "\n";
}
