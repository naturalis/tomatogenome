#!/bin/bash
DATA=data
INDIR=$DATA/candidates_wur
OUTDIR=$DATA/wur_reformatted
DB=$DATA/reference/ITAG2.3_release/ITAG2_3_genomic.fasta
GFF=$DATA/reference/ITAG2.3_release/ITAG2.3_gene_models.gff3
KEEP="046,051,054,073"
REFORMAT="perl script/reformat_wur_candidates.pl"
VERBOSITY=
GFF3=

# the WUR data is such that for each accession there is a file
# with the consensus for all candidate genes. this needs to be
# the other way around: a file for each candidate gene, with the
# accessions of interest
# $REFORMAT -indir $INDIR -outdir $OUTDIR -keep $KEEP -db $DB -gff $GFF $VERBOSITY