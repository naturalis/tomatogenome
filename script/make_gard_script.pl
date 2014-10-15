#!/usr/bin/perl
use strict;
use warnings;
use Template;
use Getopt::Long;

# process command line arguments
my ( $template, $file );
GetOptions(
	'template=s' => \$template,
	'file=s'     => \$file,
);

my $tt = Template->new;
$tt->process($template,{ 'file' => $file });