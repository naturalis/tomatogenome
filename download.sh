#!/bin/bash

# 46. S. pimpinellifolium (LA1584/PY8040) TGRC White fly resistance
# 51.S. chiemliewskii (TR00007) TGRC
# 54.  cheesmaniae (G1.1615) CGN F oxysporum resistance
# 73. S. pennellii (LYC 1831) IPK Gatersleben

HOST=ftp.ab.wur.nl
USER=bgrave
PASS=******
PREFIX=/reseq/raw
SUFFIX=illumina/pairedend_500
EXTENSION=.fq.gz
ACCESSIONS="046 051 054 073"
DOWNLOAD="perl script/downloader.pl --verbose"
OUTDIR=data/reads
for ACC in $ACCESSIONS; do
	WD="${PREFIX}/${ACC}/${SUFFIX}"
	mkdir -p $OUTDIR/$ACC
	$DOWNLOAD -host $HOST -user $USER -pass $PASS -wd $WD -outdir $OUTDIR/$ACC
done



