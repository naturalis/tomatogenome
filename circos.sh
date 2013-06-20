#!/bin/bash

# variables
GENES=Solyc02g085500.2,Solyc09g010080.2,Solyc09g018220.1.1,Solyc10g083290,Solyc09g010090,Solyc10g083300,Solyc06g008300,Solyc01g009690.1.1\
Solyc01g006550.2.1,Solyc03g082780.1.1,Solyc05g013300,Solyc07g053630.2.1,Solyc10g008160.2.1
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

# extensions for BAM files and mapDamage directories
BAMEXT=L005-ITAG2.3_genomic.fasta.sam.sorted.bam
MAPEXT=L005-ITAG2.3_genomic.fasta.sam.mapDamage

# sample identifiers
U1=1-U0015717
WAG=5-WAG0463703
U2=9-U0015716

# concatenated paths
BAMS=$DATA/reads/Sample_${U1}/${U1}_GTGGCC_${BAMEXT},$DATA/reads/Sample_${WAG}/${WAG}_CGTACG_${BAMEXT},$DATA/reads/Sample_${U2}/${U2}_ACTGAT_${BAMEXT}
MAPS=$DATA/reads/Sample_${U1}/${U1}_GTGGCC_${MAPEXT},$DATA/reads/Sample_${WAG}/${WAG}_CGTACG_${MAPEXT},$DATA/reads/Sample_${U2}/${U2}_ACTGAT_${MAPEXT}

# output goes here
WORKDIR=doc/circos

perl script/make_circos_tracks.pl --verbose \
	--gff3=$GFF3 \
	--mapdamage=$MAPS \
	--workdir=$WORKDIR \
	--refseq=$REF \
	--genes=$GENES \
	--bams=$BAMS

