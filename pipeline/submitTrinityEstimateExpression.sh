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
ref=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/TransDecoder/Trinity.fasta.transdecoder.cds
out=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/Salmon

## check vars
if [ -z $UPSCb ]; then
    echo "The UPSCb var needs to be set."
    usage
fi

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
  -e $out/$fnam.err -o $out/$fnam.out \
  -d afterok:437200 -J salmon.$fnam\
  $UPSCb/pipeline/runTrinityEstimateExpression.sh -f $f -r $in/${fnam}_2.fq.gz $ref $out

done
