#!/bin/bash

export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc472
eval `scramv1 runtime -sh`

logfile=log.txt

touch $logfile

for x in 1.01 1.05 1.25 1.5 2.0 2.5 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 16.0 26.0 51.0 76.0 101.0
do

    cat example.dat | sed "s/replaceme/$x/g" > datacard_$x.dat
    result=`combine -M Asymptotic datacard_$x.dat 2>&1 | grep "Observed Limit:" | cut -b 21-`

    echo $x $result > $logfile

    rm datacard_$x.dat
done

