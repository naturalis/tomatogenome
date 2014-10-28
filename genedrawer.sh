#!/bin/bash

# fixed arguments
WIDTH=800
HEIGHT=600
COLLAPSE=10000
GFF3=data/reference/ITAG2.3_release/ITAG2.3_gene_models.gff3
ACCESSIONS=data/genes_with_window/map.ini
SOURCE=ITAG_eugene
VERBOSITY=-v
REFERENCE=ITAG2_3
MARGIN=10000

# RF_074 = S. pennelli
# RF_072 = S. habrochaites
# RF_058 = S. arcanum
# RF_051 = S. chiemliewskii
# RF_104 = S. galapagense
# RF_047 = S. pimpinellifolium
# ITAG2_3 = S. lycopersicum Heinz 1706
INGROUP="-i RF_015 -i RF_003 -i RF_007 -i U0015716 -i U0015717 -i RF_001 -i RF_012"

# shortened for invocation
ARGS="-a $ACCESSIONS -w $WIDTH -h $HEIGHT -c $COLLAPSE -g $GFF3 -s $SOURCE -r $REFERENCE -m $MARGIN $INGROUP $VERBOSITY"

# data dir
DATA=data/genes_with_window

# list of genes
GENES='Solyc06g008300 Solyc01g009690 Solyc01g006550 Solyc11g071430 Solyc03g082780 
Solyc06g008450 Solyc09g098130 Solyc02g062560 Solyc09g018220 Solyc05g013300 
Solyc09g005090 Solyc09g005080 Solyc09g010080 Solyc06g051550 Solyc02g090730 
Solyc07g066250 Solyc09g010090 Solyc10g083300 Solyc02g085500 Solyc03g031860 
Solyc05g012020 Solyc10g008160 Solyc10g083290'

for LOCUS in $GENES; do

	# make file stem
	STEM=$DATA/$LOCUS/$LOCUS
	
	# result from hyphy GARD
	NEXUS=$STEM.trimmed.fas.out_finalout
	
	# output file
	OUT=doc/gard/$LOCUS.svg

	# invocation
	perl -Ilib script/genedrawer.pl $ARGS -l $LOCUS -n $NEXUS > $OUT
done