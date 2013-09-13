#!/bin/bash

let RUN_NUM=$1+1
BINO=$2
MASS_N=$3
REQ_JETS=$4

WORK_DIR=`pwd`

export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 project CMSSW CMSSW_5_3_3

mv higgsEnvironment.tgz CMSSW_5_3_3/src
cd CMSSW_5_3_3/src
tar -xzf higgsEnvironment.tgz
eval `scramv1 runtime -sh`
scramv1 b

mkdir $WORK_DIR/datacards
mv $WORK_DIR/datacards.tgz $WORK_DIR/datacards
cd $WORK_DIR/datacards
tar -xzf datacards.tgz

NOJET_CARD=`ls -1k $BINO\_*mN$MASS_N\_nojet.dat | sed -n $RUN_NUM\p`
JET_CARD=`ls -1k $BINO\_*mN$MASS_N\_1jet.dat | sed -n $RUN_NUM\p`

TEMP_DIR=$WORK_DIR/CMSSW_5_3_3/src/HiggsAnalysis/CombinedLimit/test/work
mkdir $TEMP_DIR
if ! $REQ_JETS ; then
        cp $NOJET_CARD $TEMP_DIR
fi
if $REQ_JETS ; then
        cp $JET_CARD $TEMP_DIR
fi
cd $WORK_DIR
rm -r datacards

mv $WORK_DIR/limit_V2 $TEMP_DIR
mv $WORK_DIR/FindNextR $TEMP_DIR

logfile=$WORK_DIR/$BINO\_log.$RUN_NUM\_$REQ_JETS

cd $TEMP_DIR

if ! $REQ_JETS ; then
        touch $logfile
        echo                                                     >> $logfile
        echo "-------------------------------------------------" >> $logfile
        echo "nojet start..."                                      >> $logfile
        echo "-------------------------------------------------" >> $logfile
        echo                                                     >> $logfile
        ./limit_V2 $NOJET_CARD                                        >> $logfile
	mv *.result.txt $WORK_DIR
fi

if $REQ_JETS ; then
        touch $logfile
        echo                                                     >> $logfile
        echo                                                     >> $logfile
        echo "-------------------------------------------------" >> $logfile
        echo "jet start..."                                     >> $logfile
        echo "-------------------------------------------------" >> $logfile
        echo                                                     >> $logfile
        ./limit_V2 ${JET_CARD}                                     >> $logfile
	mv *result.txt $WORK_DIR
fi

echo                                                     >> $logfile
echo                                                     >> $logfile
echo "-------------------------------------------------" >> $logfile
echo "moving result files and cleaning..."               >> $logfile
echo "-------------------------------------------------" >> $logfile

cd $WORK_DIR

rm -rf CMSSW_5_3_3

ls
