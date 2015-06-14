#!/bin/bash

version=$1

outfile=table/stop-bino_$version.table

[ -e $outfile ] && rm $outfile
touch $outfile

fileName=../output/stop-bino_*$version.dat.result.txt

for file in `dir -d $fileName`
do
    mStop=`grep "stop = " $file | cut -d = -f 2`
    [ "x$mStop" == "x" ] && continue

    mBino=`grep "bino = " $file | cut -d = -f 2`
    [ "x$mBino" == "x" ] && continue

    xsec=`grep "xsec = " $file | cut -d = -f 2`
    xsecError=`grep "xsec uncertainty = " $file | cut -d = -f 2 | sed 's/%//'`

    obsLimit=`grep "CLs observed = " $file | cut -d = -f 2`
    expLimit=`grep "CLs expected = " $file | cut -d = -f 2`
    exp_m1s=`grep "CLs expected m1sigma = " $file | cut -d = -f 2`
    exp_m2s=`grep "CLs expected m2sigma = " $file | cut -d = -f 2`
    exp_p1s=`grep "CLs expected p1sigma = " $file | cut -d = -f 2`
    exp_p2s=`grep "CLs expected p2sigma = " $file | cut -d = -f 2`

    #if CLs failed for whatever reason, just grab the less accurate versions?
    [ "x$obsLimit" == "x " ] && obsLimit=`grep "CLs observed asymptotic = " $file | cut -d = -f 2`
    [ "x$expLimit" == "x " ] && expLimit=`grep "CLs expected profile likelihood = " $file | cut -d = -f 2`
    [ "x$exp_m1s" == "x " ] && exp_m1s=`grep "CLs expected m1sigma profile likelihood = " $file | cut -d = -f 2`
    [ "x$exp_m2s" == "x " ] && exp_m2s=`grep "CLs expected m2sigma profile likelihood = " $file | cut -d = -f 2`
    [ "x$exp_p1s" == "x " ] && exp_p1s=`grep "CLs expected p1sigma profile likelihood = " $file | cut -d = -f 2`
    [ "x$exp_p2s" == "x " ] && exp_p2s=`grep "CLs expected p2sigma profile likelihood = " $file | cut -d = -f 2`

    [ "x$obsLimit" == "x " ] && continue
    [ "x$expLimit" == "x " ] && continue
    [ "x$exp_m1s" == "x " ] && continue
    [ "x$exp_m2s" == "x " ] && continue
    [ "x$exp_p1s" == "x " ] && continue
    [ "x$exp_p2s" == "x " ] && continue

    #echo "mStop mBino acc xsec xsecError obsLimit expLimit exp_m1s exp_m2s exp_p1s exp_p2s"

    echo -e "$mStop\t$mBino\t$xsec\t$xsecError\t$obsLimit\t$expLimit\t$exp_m1s\t$exp_m2s\t$exp_p1s\t$exp_p2s" >> $outfile

done

sed -i 's/nan/10000/' $outfile
