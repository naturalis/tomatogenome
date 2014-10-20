#!/bin/bash

# fixed arguments
WIDTH=800
HEIGHT=600
COLLAPSE=10000
GFF3=data/reference/ITAG2.3_release/ITAG2.3_gene_models.gff3
SOURCE=ITAG_eugene
VERBOSITY=-v
REFERENCE=ITAG2_3
ARGS="-w $WIDTH -h $HEIGHT -c $COLLAPSE -g $GFF3 -s $SOURCE -r $REFERENCE $VERBOSITY"

# variable arguments, will loop these
LOCUS=Solyc10g083290
NEXUS=data/genes_with_window/$LOCUS/$LOCUS.trimmed.fas.out_finalout

# invocation
perl -Ilib script/genedrawer.pl $ARGS -l $LOCUS -n $NEXUS

