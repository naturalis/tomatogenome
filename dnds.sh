#!/bin/bash

# this shell script creates NEXUS files containing a matrix and a tree ready
# for upload to datamonkey to do a BranchSiteREL analysis

# variables
# IDS="Solyc02g085500.2 Solyc09g010080.2 Solyc09g018220.1.1 Solyc10g083290 Solyc09g010090 Solyc10g083300 Solyc06g008300 Solyc01g009690.1.1 Solyc01g006550.2.1 Solyc03g082780.1.1 Solyc05g013300"
# IDS="Solyc07g053630.2.1" 
# IDS="Solyc10g008160.2.1"
# INTERVALS="SL2.40ch10:2292142-2295825"
# Solyc07g053630.2.1 = golden 1-like TF
# Solyc10g008160.2.1 = golden 2-like TF
# Solyc02g085500.2   = vorm van vrucht (rond of peervormig)
# Solyc09g010080.2   = suikergehalte van vrucht (hoog of laag)
# Solyc09g018220.1.1 = resistentiegen tegen tomato mosaic virus
# Solyc10g083290     = carbohydrate transports
# Solyc09g010090     = carbohydrate transports
# Solyc10g083300     = carbohydrate transports
# solyc06g008300     = disease resistance
# Solyc01g009690.1.1 = disease resistance
# Solyc01g006550.2.1 = disease resistance
# Solyc03g082780.1.1 = disease resistance
# Solyc05g013300     = disease resistance

# extra!!!
# Solyc11g071430 = I2 Fusarium oxysporum resistance
# Solyc09g010090 = cell wall invertase
# Solyc06g008450 = Meloidogyne + Macrosiphum resistance
# Solyc09g098130 = Thrips resistance
# Solyc02g062560 = ToMV resistance
# Solyc09g005080 = Verticillium resistance
IDS=Solyc05g013300.1\ Solyc07g053630.2\ Solyc10g008160.2


REFTAXON=ITAG2.3
OUTGROUP=WAG0463703
CUTOFF=1.00

# data files and directories
DATA=`pwd`/data
GFF3=$DATA/reference/ITAG2.3_release/ITAG2.3_gene_models.gff3
REF=$DATA/reference/ITAG2.3_release/ITAG2_3_genomic.fasta
BAM=$DATA/reads/bam
BAMLIST=`ls $BAM/*.bam`
BAMS=`echo $BAMLIST | sed -e 's/ /,/'`
CANDIDATES=$DATA/newgenes
PHYMLTREE=_phyml_tree.txt

# executables
MUSCLE=muscle
PHYML="phyml -q -p -m GTR -f m -a e -s BEST --quiet"
PERL=perl
MACSE="java -jar macse_v0.9b1.jar -g -7 -x -1 -f -30 -d 1 -s -100 -o_dna"
BIN=`pwd`/bin

# scripts
SCRIPT=`pwd`/"script"
MERGE="$PERL $SCRIPT/mergedata.pl"
CONSENSE="$PERL $SCRIPT/consense.pl -v"
DIVERGER="$PERL $SCRIPT/diverger.pl"
PROTALIGN="$PERL $SCRIPT/protalign.pl"
EXTRACT="$PERL $PROJECTS/fastq-simple-tools/script/gene_extractor.pl --verbose"
FASTA2PHYLIP="$PERL $SCRIPT/fasta2phylip.pl"
FETCHCODING="$PERL $SCRIPT/fetch_coding.pl --verbose --ignorestrand"

# build comma-separated list of BAM files
BAMLIST=""
for BAM in $BAMS; do
	BAMLIST="$BAMLIST,$BAM"
done

# fetch CDS identifiers
for ID in $IDS; do

	# create directory for all intermediate files for focal ID
	WORKDIR=$CANDIDATES/$ID
	if [ ! -d $WORKDIR ]; then
		mkdir -p $WORKDIR
	fi

	# grep the CDS identifiers from the GFF3
	CDSS=`grep $ID $GFF3 | cut -f 9 | grep ID=CDS | cut -f 1 -d ';' | sed -e 's/ID=//'`
	
	# extract CDSSs
	FASTALIST=""
	QUALLIST=""	
	for CDS in $CDSS; do
		SAFE=`echo $CDS | sed -e 's/:/_/'`
		FASTA=$WORKDIR/$SAFE.fa
		FASTALIST=$FASTALIST,$FASTA
		QUAL=$WORKDIR/$SAFE.qual
		QUALLIST=$QUALLIST,$QUAL
		if [ ! -f $FASTA ]; then
			$EXTRACT -id $CDS -bam $BAMLIST -ref $REF -gff3 $GFF3 -o $FASTA -q $QUAL
		fi
	done
	
# concatenate, consense, align, convert to phylip
	FASTA=$WORKDIR/$ID.fa
	ALN=$WORKDIR/$ID.aln
	PHYLIP=$WORKDIR/$ID.phy
	if [ ! -f $ALN ]; then
	
		# build a consensus by concatenating all the FASTA data and picking
		# the best supported base at each site, then fetch the nearest annotated
		# CDS from GenBank
#		$CONSENSE -f $FASTALIST -q $QUALLIST | $FETCHCODING -r $REFTAXON > $FASTA
		$CONSENSE -f $FASTALIST -q $QUALLIST > $FASTA

		# build a codon alignment
#		cd $BIN && $MACSE -i $FASTA -o $ALN && cd -
		
		# convert to PHYLIP, strip stop codons
#		cat $ALN | $FASTA2PHYLIP > $PHYLIP
	fi
	
# run phyml, convert to nexus
# 	NEXUS=$WORKDIR/$ID.nex
# 	if [ ! -f $NEXUS ]; then
# 		$PHYML -i $PHYLIP
# 		$MERGE -t $PHYLIP$PHYMLTREE -a $PHYLIP > $NEXUS
# 	fi	
done

# process intervals
for INTERVAL in $INTERVALS; do

	# make safe file path
	SAFE=`echo $INTERVAL | sed -e 's/:/_/'`
	
	# create directory for all intermediate files for focal interval
	WORKDIR=$CANDIDATES/$SAFE
	if [ ! -d $WORKDIR ]; then
		mkdir $WORKDIR
	fi	
	
	# extract consensus seq and quality scores for interval
	FASTA=$WORKDIR/$SAFE.fa
	QUAL=$WORKDIR/$SAFE.qual
	if [ ! -f $FASTA ]; then
		$EXTRACT -interval $INTERVAL -bam $BAMLIST -ref $REF -o $FASTA -q $QUAL
	fi
	
	# concatenate, consense, align, convert to phylip
	GBFASTA=$WORKDIR/$SAFE.gb.fa
	ALN=$WORKDIR/$SAFE.aln
	PHYLIP=$WORKDIR/$SAFE.phy
	if [ ! -f $ALN ]; then
	
		# build a consensus by concatenating all the FASTA data and picking
		# the best supported base at each site, then fetch the nearest annotated
		# CDS from GenBank
		$CONSENSE -f $FASTA -q $QUAL | $FETCHCODING -r $REFTAXON > $GBFASTA
		
		# build a codon alignment
		cd $BIN && $MACSE -i $GBFASTA -o $ALN && cd -
		
		# convert to PHYLIP, strip stop codons
		cat $ALN | $FASTA2PHYLIP > $PHYLIP
	fi	

	# run phyml, convert to nexus
	NEXUS=$WORKDIR/$SAFE.nex
	if [ ! -f $NEXUS ]; then
		$PHYML -i $PHYLIP
		$MERGE -t $PHYLIP$PHYMLTREE -a $PHYLIP > $NEXUS
	fi	
done