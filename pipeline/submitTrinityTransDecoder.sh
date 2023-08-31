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
in=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/Trinity.fasta
out=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity/TransDecoder

## check vars
if [ -z $UPSCb ]; then
    echo "The UPSCb var needs to be set."
    usage
fi

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## execute
sbatch --mail-user=$mail \
-e $out/transdecoder.err -o $out/transdecoder.out \
$UPSCb/pipeline/runTrinityTransDecoder.sh $in $out



