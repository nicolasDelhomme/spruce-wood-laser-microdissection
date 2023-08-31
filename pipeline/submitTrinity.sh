#!/bin/bash -l

## be verbose and print
set -ex

proj=u2015025
mail=nicolas.delhomme@umu.se
CPU=60
MEM=300G

## define a function
usage () {
    echo "The UPSCb env. var. needs to be set to your Git UPSCb checkout directory."
}

## source functions
source $UPSCb/src/bash/functions.sh

## process the argument
in=/mnt/picea/projects/spruce/facility/diginorm/LCM
out=/mnt/picea/projects/spruce/facility/trinity/LMD-trinity

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
sbatch -p nolimit -n $CPU --mem=$MEM  --mail-user=$mail \
-e $out/trinity.err -o $out/trinity.out -p core -n $CPU \
-t unlimited \
$UPSCb/pipeline/runTrinity.sh -p $CPU -m $MEM $out \
$(find $in -name "*_1.fq.gz" | sort | paste -s --delimiters=,) \
$(find $in -name "*_2.fq.gz" | sort | paste -s --delimiters=,)



