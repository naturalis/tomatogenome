#!/bin/bash
TREES=`ls data/candidates/*phyml_tree.txt`
FROM=ITAG2
TO=051,054,073,046,U0015717,U0015716,WAG0463703
OUTFILE=data/distances.tsv
for TREE in $TREES; do
	perl script/distances.pl -infile $TREE -from $FROM -to $TO >> $OUTFILE
done