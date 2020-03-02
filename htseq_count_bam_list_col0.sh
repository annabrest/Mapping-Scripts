#!/bin/bash


fileOfFiles=$1
gff=$2

cat $fileOfFiles | while read f; do
	echo $f;
	

	htseq-count -r name -s no -f bam -t exon -i ID $f $gff > $f.ct

done
