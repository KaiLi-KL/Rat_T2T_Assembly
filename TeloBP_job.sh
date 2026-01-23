#!/bin/bash

head -n 44 teloBP.bed > telo_chr_BP.bed
while read chr start end
do
    samtools faidx SHRSP_T2T_submit.fa "${chr}:${start}-${end}" >>"non_telo.fa"
done <trimed_BP.bed

minimap2 -ax map-ont -t 128 non_telo.fa ONT_all.fastq.gz |samtools sort -@ 128 -o ont-asm.sorted.bam -
samtools index ont-asm.sorted.bam

mkdir -p parm && cd parm
while read chr stat end
do
    outname="${chr}_parm.fq"
    samtools view -h ../ont-asm.sorted.bam "${chr}:${stat}-${end}" -o "${chr}.bam"
done <parm.bed

cd parm
for file in *.bam
do
    outname="${file%.bam}.parm.seq.list"
    samtools view -F 256 $file | cut -f 1 | sort | uniq > $outname
done


cd parm
for file in *.parm.seq.list
do
    outname="${file%.parm.seq.list}.parm.fq"
    seqkit grep -f $file ONT_all.fastq.gz -o $outname
done


while read chr stat end
do
    samtools view -h ont-asm.sorted.bam "${chr}:${stat}-${end}" -o "${chr}.bam"
done <qarm.bed

for file in *.bam
do
    outname="${file%.bam}.qarm.seq.list"
    samtools view -F 256 $file | cut -f 1 | sort | uniq > $outname
done


for file in *.qarm.seq.list
do
    outname="${file%.qarm.seq.list}.qarm.fq"
    seqkit grep -f $file ONT_all.fastq.gz -o $outname
done

file=$1
outname="${file%.seq.list}.fq"
seqkit grep -f $file ONT_all.fastq.gz -o $outname


for file in *.parm.fastq
do
    python3 teloBPCmd.py $file . --fileMode --teloNP -v 
done

for file in *parm.csv
do
    echo -n "$file: " >> summary.txt
    awk -F',' 'NR>1 {sum+=$3; count++; if($3>max) max=$3} END {if (count>0) print "Average: " sum/count "\tMax: " max; else print "No data"}' $file >> summary.txt
done

