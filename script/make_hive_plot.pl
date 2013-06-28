#!/usr/bin/perl
use strict;
use warnings;
use File::Find;
use Data::Dumper;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my $format    = 'newick';
my $p         = 0.01;
my $o3        = 1.00;
my $extension = '.phy_phyml_tree.txt';
my ( @dir, @id, $workdir, %axis );
GetOptions(
	'dir=s'       => \@dir,
	'id=s'        => \@id,
	'extension=s' => \$extension,
	'verbose+'    => \$verbosity,
	'format=s'    => \$format,
	'p=f'         => \$p,
	'o3=i'        => \$o3,
	'workdir=s'   => \$workdir,
	'axis=s'      => \%axis,
);

# instantiate logger
my $Log = Bio::Phylo::Util::Logger->new(
	'-class' => 'main',
	'-level' => $verbosity,
);

# color code by sweetness / disease resistance
my %colour = (
	'Solyc02g085500.2'   => 'black',   # vorm van vrucht (rond of peervormig)
	'Solyc07g053630.2.1' => 'lgreen',  # golden 1-like TF
	'Solyc10g008160.2.1' => 'green',   # golden 2-like TF
	'Solyc09g010080.2'   => 'dgreen',  # suikergehalte van vrucht (hoog of laag)
	'Solyc10g083290'     => 'lblue',   # carbohydrate transports
	'Solyc09g010090'     => 'dblue',   # carbohydrate transports
	'Solyc10g083300'     => 'purple',  # carbohydrate transports
	'Solyc06g008300'     => 'red',     # disease resistance
	'Solyc01g009690.1.1' => 'lred',    # disease resistance
	'Solyc01g006550.2.1' => 'dred',    # disease resistance
	'Solyc03g082780.1.1' => 'orange',  # disease resistance
	'Solyc05g013300'     => 'lorange', # disease resistance
	'Solyc09g018220.1.1' => 'dorange', # resistentiegen tegen tomato mosaic virus
);

# open file handle for links file
open my $fh, '>', "${workdir}/links.txt" or die $!;

# start tracking max values
my %max;

# traverse input directories
find( \&wanted, @dir );

# write segments file
{
	open my $segFH, '>', "${workdir}/segments.txt" or die $!;
	for my $id ( @id ) {
		printf $segFH "%s 0 %f %s %s\n", $id, $max{$id}, $id, $axis{$id};
	}
}

# this is the code reference that gets passed
# to find()
sub wanted {

	# the file matches the requested extension
	if ( /\Q$extension\E$/ ) {
		my $treefile = $_;
		my $gene = $_;
		$gene =~ s/\Q$extension\E$//;
	
		# the output from datamonkey
		my $csv = 'slacreport.csv';
		$Log->info("going to read CSV $csv");
		my %slacreport = read_csv($csv);
	
		# parse the tree
		$Log->info("going to read $format tree file $treefile");
		my $tree = parse_tree(
			'-format'     => $format,
			'-file'       => $treefile,
			'-as_project' => 1,
		);
		
		# fetch the tips of interest
		my @tips;
		push @tips, $tree->get_by_name($_) for @id;
		
		# prune all other tips
		$tree->keep_tips(\@tips);
		
		# do all pairwise comparisons
		for my $i ( 0 .. $#tips - 1 ) {
			PAIR: for my $j ( $i + 1 .. $#tips ) {
				my $source_node = $tips[$i];
				my $source_name = $source_node->get_name;
				my $target_node = $tips[$j];
				my $target_name = $target_node->get_name;
				
				# compute patristic distance
				my $distance = $source_node->calc_patristic_distance($target_node);
				$Log->info("distance between $source_name and $target_name is $distance");
				next PAIR if $distance <= 0.0000000001; # i.e. identical sequences

				# these will be numbers between 0 and 1
				my $source_y = $source_node->get_branch_length / $distance * 1000;
				my $target_y = $target_node->get_branch_length / $distance * 1000;
				
				# initialize widths
				my $source_width = 0;
				my $target_width = 0;
				
				# find the source and target in the slacreport
				for my $k ( 0 .. $#{ $slacreport{'Branch'} } ) {
					my $name = $slacreport{'Branch'}->[$k];
					if ( $name eq $source_name ) {
						$source_width = compute_width( $k, %slacreport ) / 2 * 1000;
						$Log->info("width at $name is $source_width");
					}
					elsif ( $name eq $target_name ) {
						$target_width = compute_width( $k, %slacreport ) / 2 * 1000;
						$Log->info("width at $name is $target_width");
					}
				}
				
				# print the edge
				printf $fh "%s %f %f %s %f %f color=%s # %s\n",
					$source_name, $source_y + $source_width, $source_y - $source_width,
					$target_name, $target_y + $target_width, $target_y - $target_width,
					$colour{$gene} || 'black', $gene;
				
				# increment maxes
				if ( $source_y + $source_width > $max{$source_name} ) {
					$max{$source_name} = $source_y + $source_width;					
				}
				# increment maxes
				if ( $target_y + $target_width > $max{$target_name} ) {
					$max{$target_name} = $target_y + $target_width;					
				}				
			}
		}
	}
}

# width of the ribbon is proportional to the p-value of omega class 3,
# if the overall p-value is significant and omega class 3 indicates selection
sub compute_width {
	my ( $i, %slacreport ) = @_;
	my $p_value  = $slacreport{'Corrected p-value'}->[$i];
	my $o3_value = $slacreport{'omega3'}->[$i];
	
	# the result is significant
	if ( $p_value <= $p ) {
	
		# omega 3 implies directional selection
		if ( $o3_value >= $o3 ) {
			return $slacreport{'p3'}->[$i];
		}
		else {
			$Log->info("Omega class 3 shows no directional selection $o3_value < $o3");
			return 0;
		}		
	}
	else {
		$Log->info("Corrected p-value not significant: $p_value > $p");
		return 0;
	}
}

# reads the output from data monkey
sub read_csv {
	my $csv = shift;
	
	# keys will be column headers, values will be the records in read order
	my %result;
		
	# open file handle
	open my $fh, '<', $csv or die "Can't open $csv: $!";
	my @header;
	while(<$fh>) {
		chomp;
		my @line = split /,/, $_;
		
		# on the first line, create the header
		if ( not @header ) {
			@header = @line;
			%result = map { $_ => [] } @header;
		}
		
		# on subsequent lines, record as a record
		else {
			for my $i ( 0 .. $#header ) {
				push @{ $result{$header[$i]} }, $line[$i];
			}
		}
	}
	return %result;
}