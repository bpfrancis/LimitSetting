#!/bin/bash

bino=$1
jet=$2
channel=$3

dir=../output/${channel}

outfile=${bino}_${jet}.table

[ -e $outfile ] && rm $outfile
touch $outfile

prefix=$bino
[ $bino == "bino" ] && prefix="bino_mS"

fileName=$dir/
if [ $jet == 'jet' ]; then
    fileName=$dir/${prefix}*_mN375_1jet.dat.result.txt
fi
if [ $jet == 'nojet' ]; then
    fileName=$dir/${prefix}*_nojet.dat.result.txt
fi

for file in `dir -d $fileName`
do
    mS=`grep "mS :" $file | cut -d : -f 2`

    check=x$mS
    [ "$check" == "x" ] && continue

    mG=`grep "mG :" $file | cut -d : -f 2`
    mN=`grep "mN :" $file | cut -d : -f 2`
    acc=`grep "acc :" $file | cut -d : -f 2`
    xsec=`grep "xsecValue :" $file | cut -d : -f 2`
    xsecPDFError=`grep "xsecPDFError :" $file | cut -d : -f 2`
    xsecRSErrorNeg=`grep "xsecRSErrorNeg :" $file | cut -d : -f 2`
    xsecRSErrorPos=`grep "xsecRSErrorPos :" $file | cut -d : -f 2`

    obsLimit=`grep "CLs observed =" $file | cut -d = -f 2`
    expLimit=`grep "CLs expected =" $file | cut -d = -f 2`
    exp_m1s=`grep "CLs expected m1sigma =" $file | cut -d = -f 2`
    exp_m2s=`grep "CLs expected m2sigma =" $file | cut -d = -f 2`
    exp_p1s=`grep "CLs expected p1sigma =" $file | cut -d = -f 2`
    exp_p2s=`grep "CLs expected p2sigma =" $file | cut -d = -f 2`

    check=x$acc
    [ "$check" == "x 0" ] && continue

    check=x$exp_p2s
    #[ "$check" == "x " ] && continue

#if CLs failed for whatever reason, just grab the less accurate versions?
check=x$obsLimit
[ "$check" == "x " ] && obsLimit=`grep "CLs observed asymptotic =" $file | cut -d = -f 2`
check=x$expLimit
[ "$check" == "x " ] && expLimit=`grep "CLs expected profile likelihood =" $file | cut -d = -f 2`
check=x$exp_m1s
[ "$check" == "x " ] && exp_m1s=`grep "CLs expected m1sigma profile likelihood =" $file | cut -d = -f 2`
check=x$exp_m2s
[ "$check" == "x " ] && exp_m2s=`grep "CLs expected m2sigma profile likelihood =" $file | cut -d = -f 2`
check=x$exp_p1s
[ "$check" == "x " ] && exp_p1s=`grep "CLs expected p1sigma profile likelihood =" $file | cut -d = -f 2`
check=x$exp_p2s
[ "$check" == "x " ] && exp_p2s=`grep "CLs expected p2sigma profile likelihood =" $file | cut -d = -f 2`

    #echo "mS mG mN acc xsec xsecPDFError xsecRSErrorNeg xsecRSErrorPos obsLimit expLimit exp_m1s exp_m2s exp_p1s exp_p2s"

    echo -e "$mS\t$mG\t$mN\t$acc\t$xsec\t$xsecPDFError\t$xsecRSErrorNeg\t$xsecRSErrorPos\t$obsLimit\t$expLimit\t$exp_m1s\t$exp_m2s\t$exp_p1s\t$exp_p2s" >> $outfile

done

sed -i 's/nan/10000/' $outfile
