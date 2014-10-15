#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my $wrap = 80;
my ( $infile, $start, $end, $reference );
GetOptions(
	'infile=s'    => \$infile,
	'start=i'     => \$start,
	'end=i'       => \$end,
	'wrap=i'      => \$wrap,
	'verbose+'    => \$verbosity,
	'reference=s' => \$reference,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# read the reference, compute offset and length
my ( $offset, $length );
{
	$log->info("going to compute offset and length for $reference in $infile");
	open my $fh, '<', $infile or die $!;
	my $seq;
	my $inside;
	LINE: while(<$fh>) {
		chomp;
		if ( /^>(.+)/ ) {
			my $header = $1;
			
			# start reading reference
			if ( $header eq $reference ) {
				$log->info("found $reference, start reading");
				$inside = 1;
			}
			else {
				if ( $seq ) {
					$log->info("found end of $reference, stop reading");
					last LINE;
				}
			}
		}
		
		# this will only extend the sequence if we are inside the reference
		$seq .= $_ if $inside;
	}
	$log->info("sequence length: ".length($seq));
	
	# compute where unaligned start is from beginning
	my $unaligned = 0;
	HEAD: for my $i ( 0 .. length($seq) - 1 ) {
		$unaligned++ if substr($seq,$i,1) ne '-';
		
		# have seen the first $start unaligned bases, need to aligned column index
		if ( $unaligned == $start ) {
			$offset = $i;
			$unaligned = 0;
			$log->info("offset for unaligned base $start: $offset");
			last HEAD;
		}
	}
	
	# compute where unaligned stop is from end
	TAIL: for ( my $i = length($seq) - 1; $i >= 0; $i-- ) {
		$unaligned++ if substr($seq,$i,1) ne '-';
		
		# have seen the last $end unaligned bases, need the aligned column
		# index minus the offset from the start
		if ( $unaligned == $end ) {
			$length = $i - length($seq); # is negative
			$log->info("remaining length when trimming $end bases from end: $length");
			last TAIL;
		}
	}
}

# makes template for unpack() to wrap at $wrap (default: 80) char lines
my $template = sub {
	my $length = shift;
	return "a$wrap" x ($length/$wrap) . "a*";
};

my $clean = sub {
	my $defline = shift;
	
	# remove suffixes from WUR sequences
	if ( $defline =~ /^(>RF_\d+)/ ) {
		my $clean = $1;
		return $clean;
	}
	
	# remove prefixes and suffixes from NBC sequences
	elsif ( $defline =~ /^>genes\.var\.([^\.]+)\.flt\.vcf\.fas/ ) {
		my $clean = $1;
		return ">$clean";	
	}
	
	# returned standard name for refseq
	elsif ( $defline =~ /^>ITAG2_3.+/ ) {
		return '>ITAG2_3';
	}
	
	else {
		die "Don't know what to do with $defline";
	}
};

# write out the topped and tailed alignment
{
	open my $fh, '<', $infile or die $!;
	my $seq;
	while(<$fh>) {
		chomp;
		if ( /^>/ ) {
		
			# print topped and tailed seq, if any
			if ( $seq ) {
				my $subseq = substr $seq, $offset, $length;
				print join "\n", unpack $template->(length $subseq), $subseq;
				print "\n";
			}
			
			# print cleaned definition line
			print $clean->($_), "\n";
			$seq = ''; # reset
		}
		else {
			$seq .= $_; # extend
		}
	}
	
	# print the final seq
	my $subseq = substr $seq, $offset, $length;
	print join "\n", unpack $template->(length $subseq), $subseq;
	print "\n";	
}