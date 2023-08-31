#!/bin/bash -l

## be verbose and print
set -ex

proj=u2015025
mail=nicolas.delhomme@umu.se

## define a function
usage () {
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
}

## source functions
source $UPSCb/src/bash/functions.sh

## process the argument
in=/mnt/picea/projects/spruce/facility/diginorm/LCM
ref=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/TransDecoder/Trinity.fasta.transdecoder.cds.salmon_quasi.idx
out=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/Salmon

## check vars
if [ -z $UPSCb ]; then
    abort "The UPSCb var needs to be set."
fi

## load the tool - we actually use singularity
#module load bioinfo-tools salmon

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## for every file
cd $out
for f in $(find $in -name "*_1.fq.gz"); do
  fnam=$(basename ${f/_1.fq.gz/})
  
  ## execute
  sbatch --mail-user=$mail \
  -e $out/$fnam.err -o $out/$fnam.out -J salmon.$fnam \
  $UPSCb/pipeline/runSalmon.sh $ref $f $in/${fnam}_2.fq.gz $out

done
