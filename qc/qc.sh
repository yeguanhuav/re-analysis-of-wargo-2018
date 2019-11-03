#!/usr/bin/bash
run_path=/home/yeguanhua/Wargo/PRJEB22894/ERR2162225/demux/
cd $run_path
ls $run_path
for file in $(ls $run_path)
do
	docker run -u $UID:$UID --rm -v $PWD:/data/ quay.io/biocontainers/fastqc:0.11.7--4 sh -c "mkdir -p /data/fastqc && fastqc --threads 4 --outdir /data/fastqc --noextract /data/$file"
done
docker run -u $UID:$UID --rm -v $PWD:/data/ quay.io/biocontainers/multiqc:1.6--py27h24bf2e0_0 sh -c "multiqc   -o /data/multiqc /data/fastqc"
