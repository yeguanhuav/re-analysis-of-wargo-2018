#!/usr/bin/bash
for sam in `less fecal_16S_sample_id.txt`
do
	echo $sam
	less /home/yeguanhua/Wargo/PRJEB22894/ERR2162225/fecal_16S.fastq | grep $sam -A 3 > $sam.fastq
done
