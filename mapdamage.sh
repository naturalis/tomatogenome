#!/bin/bash

# run me as $ nohup time ./runall.sh &

# now using the reference from ITAG, so that
# we might be able to crossreference regions
# interesting results in our data (e.g. strange
# read depths, areas under selection) with
# known features
INDICES=data/reference/ITAG2.3_release
FASTAS=`ls $INDICES/*.fasta`
READS=data/reads

# all samples
SAMPLES=`ls $READS | grep 'naturalis_'`

# recreate BWA index if not exists
for FASTA in $FASTAS; do
	if [ ! -e $FASTA.bwt ]; then
		echo "going to index $FASTA"

		# Warning: "-a bwtsw" does not work for short genomes, 
		# while "-a is" and "-a div" do not work not for long 
		# genomes. Please choose "-a" according to the length 
		# of the genome.
		time bwa index -a bwtsw $FASTA
	else 
		echo "$FASTA already indexed"
	fi
done

# BWA align reads
for SAMPLE in $SAMPLES; do
	echo "going to process sample $SAMPLE"

	# unzip the archives
	ZIPS=`ls $READS/$SAMPLE/*.gz`
	gunzip $ZIPS

	# align the reads
	FASTQS=`ls $READS/$SAMPLE/*.fastq`
	for FASTA in $FASTAS; do
		echo "going to analyse refseq $FASTA"

		# lists of produced files
		SAIS=""
		SAM=""
		for FASTQ in $FASTQS; do

			# create new name
			LOCALFASTA=`echo $FASTA | sed -e 's/.*\///'`
			LOCALFASTQ=`echo $FASTQ | sed -e 's/.*\///'`
			OUTFILE=$READS/$SAMPLE/$LOCALFASTQ-$LOCALFASTA.sai
			SAIS="$SAIS $OUTFILE"
			SAM=`echo $OUTFILE | sed -e "s/_R.*/-$LOCALFASTA.sam/"`

			# note: we don't do basic QC here, because that might mean
			# that the mate pairs in the FASTQ files go out of order,
			# which will result in the bwa sampe step taking an inordinate
			# amount of time

			# do bwa aln if needed
			if [ ! -e $OUTFILE ]; then
				echo "going to align $FASTQ against $FASTA"

				# use 6 threads
				time bwa aln -t 6 $FASTA $FASTQ > $OUTFILE
			else
				echo "alignment $OUTFILE already created"
			fi
		done
		echo "*** DONE BUILDING ALIGNMENTS FOR $FASTQCS ***"

		# do bwa sampe if needed
		if [ ! -e $SAM ]; then

			# create paired-end SAM file
			echo "going to run bwa sampe $FASTA $SAIS $FASTQS > $SAM"
			time bwa sampe $FASTA $SAIS $FASTQS > $SAM

			# -bS   = input is SAM, output is BAM
			# -F 4  = remove unmapped reads
			# -q 25 = remove reads with mapping qual < 25
			echo "going to run samtools view -bS -F 4 -q 25 $SAM > $SAM.filtered"
			time samtools view -bS -F 4 -q 25 $SAM > $SAM.filtered

			# sorting is needed for indexing
			echo "going to run samtools sort $SAM.filtered $SAM.sorted"
			time samtools sort $SAM.filtered $SAM.sorted

			# this should result in faster processing
			echo "going to run samtools index $SAM.sorted.bam"
			time samtools index $SAM.sorted.bam
		else
			echo "sam file $SAM already created"
		fi
		echo "*** DONE BUILDING SAM FILE ***"

		# do mapDamage if needed
		MAPDAMAGEDIR=$SAM.mapDamage
		if [ ! -e $MAPDAMAGEDIR ]; then
			echo "going to run mapDamage map -c -i $SAM.sorted.bam -d $MAPDAMAGEDIR -r $FASTA -l 101 -u -k"

			# -c = complete analysis
			# -i = input file
			# -d = working directory
			# -r = reference file
			# -l = read length (should be 100, maybe off-by-one error in mapDamage?)
			## XXX probably not needed: -u = skip duplicates based on X1/XT bwa headers
			# -n = random sample of reads to analyse
			# -k = keep bed file (maybe useful for circos)
			time mapDamage map -c -i $SAM.sorted.bam -d $MAPDAMAGEDIR -r $FASTA -l 101 -n 10000 -k
		else
			echo "mapDamage already run on dir $MAPDAMAGEDIR"
		fi
		echo "*** DONE RUNNING MAPDAMAGE ***"		
	done
done
