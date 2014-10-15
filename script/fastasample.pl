#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my ( $infile, $begin, $end );
GetOptions(
	'infile=s' => \$infile,
	'begin=i'  => \$begin,
	'end=i'    => \$end,
	'verbose+' => \$verbosity,
);

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);
$log->info("going to extract $begin..$end from $infile");

# start reading from file
open my $fh, '<', $infile or die $!;
my $pos;
my $line = 0;
LINE: while(<$fh>) {
	chomp;
	if ( /^>/ ) {
		$pos = 0;
		print "\n" if $line > 0;
		print $_, "\n";
		$log->debug("defline: $_");
	}
	else {
		my $seq = $_;
		
		# won't go into window, or already beyond it
		if ( ( $pos + length $seq ) < $begin || ( $pos > $end ) ) {
			#$log->debug("$pos outside window") unless $line % 10000;
		}
		
		# inside of window
		elsif ( ( $pos > $begin ) && ( ( $pos + length $seq ) < $end ) ) {
			print $seq;			
			#$log->debug("$pos inside window") unless $line % 100;
		}
		
		# window start
		elsif ( ( $pos <= $begin ) && ( ( $pos + length $seq ) > $begin ) ) {
			my $offset = $begin - $pos;
			print substr $seq, $offset;		
			$log->debug("$pos intering into window");
		}
		
		# window end
		elsif ( ( $pos < $end ) && ( ( $pos + length $seq ) >= $end ) ) {
			my $length = $end - $pos;
			print substr $seq, 0, $length;
			$log->debug("$pos leaving window");
		}
				
		$pos += length $seq;
	}
	$line++;
}