#!/usr/bin/perl
use strict;
use warnings;
use Net::FTP;
use Getopt::Long;
use Bio::Phylo::Util::Logger ':simple';

# process command line arguments
my $verbosity = WARN "start";
my $host = 'ftp.ab.wur.nl';
my $wd; # server side working directory
my ( $user, $pass ); # login credentials
my $extension = '.fq.gz';
my $outdir; # dir to write to
GetOptions(
	'verbose+' => \$verbosity,
	'host=s'   => \$host,
	'wd=s'     => \$wd,
	'user=s'   => \$user,
	'pass=s'   => \$pass,
	'outdir=s' => \$outdir,
	'extension=s' => \$extension,
);

# update verbosity
Bio::Phylo::Util::Logger::VERBOSE( '-level' => $verbosity, '-class' => 'main' );

# instantiate ftp object
my $ftp = Net::FTP->new($host) or die "Can't connect to $host: $@";
INFO "FTP connected to host $host";

# log in
$ftp->login( $user => $pass ) or die "Can't login with u=$user p=$pass: $@";
INFO "Login successful as $user";

# change to directory
$ftp->cwd( $wd ) or die "Can't cwd to $wd: $@";
INFO "Changed working directory to $wd";

# search for matches
my @matches = grep { /\Q$extension\E$/ } $ftp->ls();
INFO "Found ".scalar(@matches)." files with extension $extension";

# start downloading
for my $file ( @matches ) {
	INFO "Going to download $file into directory $outdir";
	$ftp->get( $file, "$outdir/$file" );
}
INFO "Done.";