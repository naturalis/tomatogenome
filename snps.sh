#!/bin/bash

# variables
IDS="Solyc02g085500.2 Solyc09g010080.2 Solyc09g018220.1.1 Solyc10g083290 Solyc09g010090 Solyc10g083300 Solyc06g008300 Solyc01g009690.1.1 Solyc01g006550.2.1 Solyc03g082780.1.1 Solyc05g013300 Solyc07g053630.2.1 Solyc10g008160.2.1"
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

# data files and directories
DATA=`pwd`/data
GFF3=$DATA/reference/ITAG2.3_release/ITAG2.3_gene_models.gff3
REF=$DATA/reference/ITAG2.3_release/ITAG2.3_genomic.fasta
READS=$DATA/reads
CANDIDATES=$DATA/candidates

# suffix for BAM files
BAM=L005-ITAG2.3_genomic.fasta.sam.sorted.bam

# sample identifiers
U1=1-U0015717
WAG=5-WAG0463703
U2=9-U0015716

# concatenated paths
BAMS="$READS/Sample_${U1}/${U1}_GTGGCC_${BAM} $READS/Sample_${U2}/${U2}_ACTGAT_${BAM} $READS/Sample_${WAG}/${WAG}_CGTACG_${BAM}"

# iterate over top-level identifiers
for ID in $IDS; do

	# grep the CDS identifiers from the GFF3
	CDSS=`grep $ID $GFF3 | cut -f 9 | grep ID=CDS | cut -f 1 -d ';' | sed -e 's/ID=//'`
	
	# iterate over CDSs
	for CDS in $CDSS; do
		START=`grep $CDS $GFF3 | cut -f 4 | head -1`
		STOP=`grep $CDS $GFF3 | cut -f 5 | head -1`
		CHROMO=`grep $CDS $GFF3 | cut -f 1 | head -1`
		REGION="${CHROMO}:${START}-${STOP}"
		echo "*** $ID $REGION ***"	
		SAFE=`echo $CDS | sed -e 's/CDS:/CDS_/'`
		
		# BCF raw output file
		RAW=$CANDIDATES/$ID/$SAFE.raw.bcf
		
		# VCF filtered output file
		FILTERED=$CANDIDATES/$ID/$SAFE.flt.vcf
		
		# run samtools and bcftools
		samtools mpileup -r $REGION -u -f $REF $BAMS 2> /dev/null | bcftools view -bvcg - > $RAW 2> /dev/null
		bcftools view $RAW | vcfutils.pl varFilter -D100 > $FILTERED 2> /dev/null
		
		# check for indels
		grep "INDEL;" $FILTERED
	done
done