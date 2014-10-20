#!/bin/bash

# the general idea is that we want to have alignments that span our genes of interest
# with an additional window of $SIZE - $ENDS upstream as well as downstream. we start
# by extracting a larger window and aligning that in order to allow for indels. we 
# subsequently splice our region of interest out of that larger alignment.

# ANCESTORS=S._habrochaites,S._pennellii,S._peruvianum,S._Arcanum,S._galapagense,S._pimpinellifolium,S._chiemliewskii
# DESCENDANTS=S._lycopersicum_old1,S._lycopersicum_old2,Moneymaker,Garderners_Delight,Katinka_Cherry,Sonato,Momotaro,Heinz

TEMPLATE=script/GARD.bf.tmpl
GARD_HELPER=script/GARD_GA_CHC.ibf
GFF3=data/reference/ITAG2.3_release/ITAG2.3_gene_models.gff3
CONDIR=data/reference/deleteme/out
OUTDIR=data/genes_with_window
SIZE=11000
ENDS=1000
GENES='Solyc06g008300 Solyc01g009690 Solyc01g006550 Solyc11g071430 Solyc03g082780 
Solyc06g008450 Solyc09g098130 Solyc02g062560 Solyc09g018220 Solyc05g013300 
Solyc09g005090 Solyc09g005080 Solyc09g010080 Solyc06g051550 Solyc02g090730 
Solyc07g066250 Solyc09g010090 Solyc10g083300 Solyc02g085500 Solyc03g031860 
Solyc05g012020 Solyc10g008160 Solyc10g083290'

# iterate over genes
for GENE in $GENES; do

	# parse GFF3
	LINE=`grep "ID=gene:$GENE" $GFF3 | cut -f 1,4,5`
	CHR=`echo $LINE | cut -f 1 -d ' ' | sed -e 's/SL2.40ch//'`
	
	# numbers are off by one for consensus sequences (1=miscellaneous sequences, so
	# chromosome 1 = 2.fas)
	FASTA=`echo "$CHR+1" | bc`
	
	# begin and end coordinates of gene
	GBEGIN=`echo $LINE | cut -f 2 -d ' '`
	GEND=`echo $LINE | cut -f 3 -d ' '`
	
	# begin and end coordinates of window
	BEGIN=`echo "$GBEGIN-$SIZE" | bc`
	END=`echo "$GEND+$SIZE" | bc`
	
	# make working directory
	if [ ! -e $OUTDIR/$GENE	]; then
		mkdir $OUTDIR/$GENE	
	fi
	
	# extract and align window
	if [ ! -e $OUTDIR/$GENE/$GENE.fas ]; then
		echo "aligning $OUTDIR/$GENE/$GENE.fas"
		perl script/fastasample.pl -i $CONDIR/$FASTA.fas -b $BEGIN -e $END | muscle -diags -maxiters 1 > $OUTDIR/$GENE/$GENE.fas
	fi
	
	# top 'n' tail the alignment
	if [ ! -e $OUTDIR/$GENE/$GENE.trimmed.fas ]; then
		echo "trimming $OUTDIR/$GENE/$GENE.trimmed.fas"
		perl script/topntail.pl -i $OUTDIR/$GENE/$GENE.fas -s $ENDS -e $ENDS -wrap 80 -reference ITAG2_3_genomic.fas > $OUTDIR/$GENE/$GENE.trimmed.fas
	fi
	
	# make GARD batch file
	if [ ! -e $OUTDIR/$GENE/GARD.bf ]; then
		echo "making GARD batch file $OUTDIR/$GENE/GARD.bf"
		perl script/make_gard_script.pl -t $TEMPLATE -f `pwd`/$OUTDIR/$GENE/$GENE.trimmed.fas > $OUTDIR/$GENE/GARD.bf
		cp $GARD_HELPER $OUTDIR/$GENE/
	fi
	
	# run hyphy
	if [ ! -e $OUTDIR/$GENE/$GENE.trimmed.fas.out ]; then
		echo "launching HYPHY for $GENE"
		cd $OUTDIR/$GENE/
		mpirun -np 4 HYPHYMPI GARD.bf
		cd -
	fi		
done