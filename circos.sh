#!/bin/bash

# variables
GENES=gene:Solyc02g085500.2,gene:Solyc09g010080.2,gene:Solyc09g018220.1,gene:Solyc10g083290.1,\
gene:Solyc09g010090.2,gene:Solyc10g083300.1,gene:Solyc06g008300.2,gene:Solyc01g009690.1,\
gene:Solyc01g006550.2,gene:Solyc03g082780.1,gene:Solyc05g013300.1,gene:Solyc07g053630.2,\
gene:Solyc10g008160.2
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

# extensions for BAM files and mapDamage directories
BAM=L005-ITAG2.3_genomic.fasta.sam.sorted.bam
MAP=L005-ITAG2.3_genomic.fasta.sam.mapDamage

# sample identifiers
U1=1-U0015717
WAG=5-WAG0463703
U2=9-U0015716

# concatenated paths
BAMS=$READS/Sample_${U1}/${U1}_GTGGCC_${BAM},$READS/Sample_${U2}/${U2}_ACTGAT_${BAM},$READS/Sample_${WAG}/${WAG}_CGTACG_${BAM}
MAPS=$READS/Sample_${U1}/${U1}_GTGGCC_${MAP},$READS/Sample_${U2}/${U2}_ACTGAT_${MAP},$READS/Sample_${WAG}/${WAG}_CGTACG_${MAP}

# output goes here
WORKDIR=doc/circos

# add --compute to recalculate binned data files. this takes a long time!
ARGS="-gff3 $GFF3 -m $MAPS -w $WORKDIR -ref $REF -genes $GENES -b $BAMS --verbose --compute"
perl script/make_circos_tracks.pl $ARGS
circos -conf $WORKDIR/circos.conf -debug_group text
open circos.png