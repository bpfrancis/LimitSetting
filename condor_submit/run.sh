#!/bin/bash

let RUN_NUM=$1+1
GRID=$2
CATEGORY=$3

WORK_DIR=`pwd`

source /cvmfs/cms.cern.ch/cmsset_default.csh
export SCRAM_ARCH=slc6_amd64_gcc481
scramv1 project CMSSW CMSSW_7_1_5

mv higgsEnvironment.tgz CMSSW_7_1_5/
cd CMSSW_7_1_5/
tar -xzf higgsEnvironment.tgz
cd src/
eval `scramv1 runtime -sh`

mkdir $WORK_DIR/datacards
mv $WORK_DIR/datacards.tgz $WORK_DIR/datacards
cd $WORK_DIR/datacards
tar -xzf datacards.tgz

CARD=`ls -1k stop-bino_*$CATEGORY.dat | sed -n $RUN_NUM\p`

TEMP_DIR=$WORK_DIR/CMSSW_7_1_5/src/HiggsAnalysis/CombinedLimit/test/work
mkdir $TEMP_DIR
cp $CARD $TEMP_DIR
cd $WORK_DIR
rm -r datacards

mv $WORK_DIR/limit $TEMP_DIR

logfile=$WORK_DIR/$GRID\_log.$RUN_NUM\_$CATEGORY

cd $TEMP_DIR

touch $logfile
echo                                                     >> $logfile
echo "-------------------------------------------------" >> $logfile
echo "$CATEGORY start..."                                      >> $logfile
echo "-------------------------------------------------" >> $logfile
echo                                                     >> $logfile
./limit $CARD                                        >> $logfile
mv *.result.txt $WORK_DIR

cd $WORK_DIR

rm -rf CMSSW_7_1_5

