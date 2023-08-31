#!/bin/bash

set -ex

proj="spruce"
mail="nicolas.delhomme@umu.se"
in=/mnt/picea/projects/spruce/laser-capture-microdissection-olga/raw
out=/mnt/picea/projects/spruce/laser-capture-microdissection-olga/
genome=/mnt/picea/storage/reference/Picea-abies/v1.1/indices/STAR/Pabies01-genome
gff3=/mnt/picea/storage/reference/Picea-abies/v1.0/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene.gff3
start=2 #start at FastQC
end=8 #end at the HTSeq step 
mem=fat

module load bioinfo-tools FastQC sortmerna Trimmomatic samtools star htseq

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*_[1,2].fastq.gz"`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do
## non strand spec data
bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -m $mem -g $genome -G $gff3 -H $gff3 -k $proj $mail ${line}_1.fastq.gz ${line}_2.fastq.gz $out
done
