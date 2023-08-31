#!/bin/bash

proj=u2015025
mail=nicolas.delhomme@umu.se
#CPU=64
CPU=20
MEM=300G

#in=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/Trinity.fasta
in=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/TransDecoder/Trinity.fasta.transdecoder.cds
out=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/GMAP
inxDir=/mnt/picea/storage/reference/Picea-abies/v1.0/indices/GMAP
inxName=Pabies1.0

if [ ! -d $out ]; then
  mkdir -p $out
fi

module load bioinfo-tools gmap-gsnap

sbatch -A $proj -n $CPU --mem=$MEM --mail-user=$mail -o $out/gmapl.out -e $out/gmapl.err $UPSCb/pipeline/runGmapl.sh $in $inxDir $inxName $out
