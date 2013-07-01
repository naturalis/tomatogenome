#!/bin/bash
DATA=`pwd`/data/candidates
WORK=`pwd`/doc/hiveplot
perl script/make_hive_plot.pl -v -d $DATA -w $WORK \
	-i U0015716 -i WAG0463703 -i U0015717 \
	-a U0015717=black -a U0015716=blue -a WAG0463703=yellow
CONFS=`ls $WORK/*.conf`
for CONF in $CONFS; do
	linnet -conf $CONF
done