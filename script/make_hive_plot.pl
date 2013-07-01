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

# traverse input directories
find( \&wanted, @dir );

# this is the code reference that gets passed
# to find()
sub wanted {

	# the file matches the requested extension
	if ( /\Q$extension\E$/ ) {
		my $treefile = $_;
		my $gene = $_;
		$gene =~ s/\Q$extension\E$//;
				
		# open file handle for links file
		my $links_file = "${workdir}/links-${gene}.txt";
		open my $fh, '>', $links_file or die "Can't open $links_file: $!";
	
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
		
		# calculate the sum of branch lengths of the remainder
		my $length = $tree->calc_tree_length;
		$Log->info("total tree length is $length");
		
		# start tracking max values
		my %max = map { $_ => 1000 } @id;		
		# do all pairwise comparisons
		for my $i ( 0 .. $#tips - 1 ) {
			PAIR: for my $j ( $i + 1 .. $#tips ) {
				my $source_node = $tips[$i];
				my $source_name = $source_node->get_name;
				my $target_node = $tips[$j];
				my $target_name = $target_node->get_name;

				# these will be numbers between 0 and 1
				my $source_y = $source_node->get_branch_length / $length * 1000;
				my $target_y = $target_node->get_branch_length / $length * 1000;
				$Log->info("$source_name y: $source_y");
				$Log->info("$target_name y: $target_y");
				
				# the widths are the proportions of the alignments that are under
				# direction selection. being proportions these are numbers between
				# 0 and 1, but we divide them by 2 and multiply by 1000 so we can
				# more concisely add them to and subtract them from the coordinates
				# of where the ribbon touches a hive axis
				my ($source_width,$target_width) = compute_ribbon_widths($source_name,$target_name,%slacreport);
				$Log->info("$source_name width: $source_width");
				$Log->info("$target_name width: $target_width");
				
				# print the edge
				printf $fh "%s %f %f %s %f %f color=%s # %s\n",
					$source_name, $source_y + $source_width, $source_y - $source_width,
					$target_name, $target_y + $target_width, $target_y - $target_width,
					'black', $gene;
				
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
		
		# write segments file
		open my $segFH, '>', "${workdir}/segments-${gene}.txt" or die $!;
		for my $id ( @id ) {
			printf $segFH "%s 0 %f %s %s\n", $id, $max{$id}, $id, $axis{$id};
		}
		
		# write config file
		print_conf($gene);	
		$Log->info("*** done with $gene ***");	
	}
}

# prints the config file for hiveplot
sub print_conf {
my $gene = shift;
my $template = '<colors>
<<include etc/colors.conf>>
<<include etc/brewer.conf>>
</colors>

<image>
	size = 1000
	dir  = .
	file = hiveplot-%s.png
	png  = yes
	background        = white
	auto_alpha_steps  = 10
</image>

<segments>
	file    = doc/hiveplot/segments-%s.txt
	width   = 5
	spacing = __$CONF{segments}{width}*3__
	radius  = __$CONF{image}{size}/20__
</segments>

<links>
	<link>
		file         = doc/hiveplot/links-%s.txt
		ribbon       = yes
		offset_start = __$CONF{segments}{width}__
		offset_end   = __$CONF{segments}{width}__
		crest        = __$CONF{image}{size}/20__
		bezier_radius_power = 1
	</link>
</links>

<scale>
	pixsize = 2.5
</scale>

<axes>

	<axis a>
		angle      = 0
		scale      = 1
		reverse    = no
		segments   = U0015717
		scale_norm = 1000
	</axis>

	<axis b>
		scale      = 1
		angle      = 120
		reverse    = no
		segments   = WAG0463703
		scale_norm = 1000
	</axis>

	<axis c>
		scale      = 1
		angle      = 240
		reverse    = no
		segments   = U0015716
		scale_norm = 1000
	</axis>
</axes>';
open my $fh, '>', "${workdir}/linnet-${gene}.conf" or die $!;
printf $fh $template, $gene, $gene, $gene;
}

# dispatches to compute_width for both the source and target taxon
sub compute_ribbon_widths {
	my ( $source_name, $target_name, %slacreport ) = @_;

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
	return $source_width, $target_width;
}

# width of the ribbon is proportional to the p-value of omega class 3,
# if the overall p-value is significant and omega class 3 indicates selection.
# returns a number between 0 and 1
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