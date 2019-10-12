#!/bin/sh

WORK=$PWD
for sim in `find $1 -name \*.sim`
do
    simdir=`dirname $sim`
    cd $simdir
    for plot in `grep filename mccode.sim | cut -f2- -d\:`
    do
	echo Working on $plot
	mcplot-svg --nobrowse --autosize --libpath='.' $plot
    done
    cd $WORK
done
