#!/bin/bash

## error and verbose
set -ex

# load functions
source $UPSCb/src/bash/functions.sh

## check for the UPSCb env. var.                                                                                                  
if [ -z $UPSCb ]; then
    abort "You need to set the UPSCb environment variable"
fi

## default args
in=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/GMAP
out=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/bedtools
ref=/mnt/picea/storage/reference/Picea-abies/v1.0/GBrowse/Pabies1.0/Gene_Prediction_Transcript_assemblies/Eugene-gene-only.gff3
mail=nicolas.delhomme@umu.se
prefix="Pabies1.0-Trinity.2"

# the data is prep'ed as follows (first time only):
cd $in
if [ ! -f $prefix.uniq.gene-only ]; then
  grep "gene" $prefix.uniq > $prefix.uniq.gene-only
fi

if [ ! -f $prefix.mult.gene-only ]; then
  grep "gene" $prefix.mult > $prefix.mult.gene-only
fi

if [ ! -f $prefix.transloc.gene-only ]; then
  grep "gene" $prefix.transloc > $prefix.transloc.gene-only
fi

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## module
module load bioinfo-tools BEDTools

## prepare
for f in $(find $in -name "$prefix*.gene-only"); do
  fnam=$(basename $f)
  sbatch --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out $UPSCb/pipeline/runBedToolsIntersect.sh $f $ref $out -s -f 0.3 -wo
done
