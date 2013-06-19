#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# global variables
my ( %FastaData, %QualData ); # concatenated data
my ( @FastaFiles, @QualFiles ); # concatenated, sorted lists of files
my $Verbosity = 1; # only emit warnings

# basic logging functionality
sub LOG ($$) {
	my ($msg,$method) = @_;
	my ( $package, $file1up, $line1up, $sub ) = caller( 2 );
	my ( $pack0up, $file, $line, $sub0up )    = caller( 1 );
	my $log = sprintf( "%s %s [%s %s] - %s\n", uc $method, $sub || '', $0, $line, $msg );
	print STDERR $log;
}
sub DEBUG ($) { LOG shift, 'DEBUG' if $Verbosity >= 3 }
sub INFO ($)  { LOG shift, 'INFO'  if $Verbosity >= 2 }
sub WARN ($)  { LOG shift, 'WARN'  if $Verbosity >= 1 }

# check command line arguments
sub check_args {

	# process command line arguments
	my ( $fastalist, $quallist );
	GetOptions(
		'fastalist=s' => \$fastalist,
		'quallist=s'  => \$quallist,
		'verbose+'    => \$Verbosity,
	);

	# make sorted lists of files
	@FastaFiles = sort { cds_index($a) <=> cds_index($b) } 
	              grep { /\S/ } 
	              split /,/, $fastalist;
	INFO "sorted FASTA files are @FastaFiles";
	
	@QualFiles  = sort { cds_index($a) <=> cds_index($b) } 
	              grep { /\S/ } 
	              split /,/, $quallist;
	INFO "sorted qual files are @QualFiles";	
}

# extracts the index number as per the ITAG conventions from a file name
sub cds_index {
	my $file = shift;
	$file =~ s/^.+\.(\d)\.[a-z]+$/$1/;
	return $file;
}

# read FASTA formatted list of files, concatenate data structures
sub read_files {

	# read the input files
	for my $i ( 0 .. $#FastaFiles ) {
		my $fasta = $FastaFiles[$i];
		my $qual  = $QualFiles[$i];
		INFO "$i - going to read FASTA file $fasta and qual file $qual";

		my %fasta;
		read_fasta($fasta,0,\%fasta);

		my %quality;
		read_fasta($qual,1,\%quality);

		# here we assemble a growing data structure	
		if ( $i > 0 ) {
			concat_data( \%FastaData, %fasta );
			concat_data( \%QualData, %quality );
			INFO "growing concatenated data structure";
		}
		else {
			%FastaData = %fasta;
			%QualData = %quality;
			INFO "initializing concatenated data structure";
		}
	}

}

# write data such that for each site we pick the base with the highest phred score
sub write_data {

	# iterate over the ids, which have now lumped duplicated records
	for my $id ( keys %FastaData ) {

		# print FASTA header
		print ">$id\n";
	
		# start iterating until no duplicate has data anymore
		my $i = 0;
		CELL: while(1) {
	
			# this will be keyed on bases, values are phred scores
			my %basequal;
			my $j = 0; # row counter, to get matching qual record
			my $data = 0; # flag to quit when no data in any row
		
			# iterate over duplicated rows
			for my $j ( 0 .. $#{ $FastaData{$id} } ) {
				if ( my $base = $FastaData{$id}->[$j]->[$i] ) {
					$basequal{$base} = $QualData{$id}->[$j]->[$i];
					$data++;
				}
			}
		
			# no data in any row, quit
			last CELL unless $data;
		
			# sort bases in phred score in descending order
			my ($base) = sort { $basequal{$b} <=> $basequal{$a} } keys %basequal;
			print $base;
			$i++; # next site
		}
		print "\n"; # end of record
		INFO "wrote data for $id";
	}
}

# concatenate hash of array of array
sub concat_data {
	INFO "going to concatenate data";
	my ( $total, %data ) = @_;	
	
	# iterate over data in current hash
	for my $id ( keys %data ) {
		for my $i ( 0 .. $#{ $data{$id} } ) {
			DEBUG "growing row $i for $id";
			push @{ $total->{$id}->[$i] }, @{ $data{$id}->[$i] };
		}
	}
}

# reads a FASTA formatted sequence or quality file
sub read_fasta {
	my ($file,$isphred,$data) = @_;
	open my $fh, '<', $file or die $!;
	INFO "going to read $file containing " . ( $isphred ? "qual scores" : "sequences" );
	my $current;
	
	# start reading
	while(<$fh>) {
		chomp;
		
		# found FASTA header line
		if ( />(\S+)/ ) {
		
			# id is first "word"
			my $id = $1;
			DEBUG "first word in FASTA definition line is $id";
			
			# strip hyphenated prefix that identifies duplicate samples
			$id =~ s/.+-//;
			DEBUG "stripped hyphenated prefix, keeping $id";
			
			# instantiate hash of arrays record
			$data->{$id} = [] unless $data->{$id};
			
			# add additional array to make 2D
			push @{ $data->{$id} }, [];
			
			# cache the ID
			$current = $id;
		}
		
		# processing sequential data
		else {
		
			# split the current line in cells, either on whitespace (if phred) 
			# or on each character (otherwise, i.e. if DNA data)
			my @cells = $isphred ? split /\s+/, $_ : split //, $_;
			
			# append the focal cells to the current record
			push @{ $data->{$current}->[-1] }, grep { /\S/ } @cells;
			DEBUG "data added to $current: @cells";
		}
	}
}

# this is the top level sub, drill down from here
sub main {
	INFO "checking command line arguments";
	check_args();
	INFO "reading lists of files";
	read_files();
	INFO "writing data";
	write_data();
}
main();