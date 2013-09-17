#!/bin/bash

let RUN_NUM=$1+1
GRID=$2
CATEGORY=$3

WORK_DIR=`pwd`

export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc472
scramv1 project CMSSW CMSSW_6_1_1

mv higgsEnvironmet.tgz CMSSW_6_1_1/src
cd CMSSW_6_1_1/src
tar -xzf higgsEnvironment.tgz
eval `scramv1 runtime -sh`
scramv1 b

mkdir $WORK_DIR/datacards
mv $WORK_DIR/datacards.tgz $WORK_DIR/datacards
cd $WORK_DIR/datacards
tar -xzf datacards.tgz
cd datacards/

CARD=`ls -1k stop-bino_*$CATEGORY.dat | sed -n $RUN_NUM\p`

TEMP_DIR=$WORK_DIR/CMSSW_6_1_1/src/HiggsAnalysis/CombinedLimit/test/work
mkdir $TEMP_DIR
cp $CARD $TEMP_DIR
cd $WORK_DIR
rm -r datacards

mv $WORK_DIR/limit_V2 $TEMP_DIR
mv $WORK_DIR/FindNextR $TEMP_DIR

logfile=$WORK_DIR/$GRID\_log.$RUN_NUM\_$CATEGORY

cd $TEMP_DIR

touch $logfile
echo                                                     >> $logfile
echo "-------------------------------------------------" >> $logfile
echo "$CATEGORY start..."                                      >> $logfile
echo "-------------------------------------------------" >> $logfile
echo                                                     >> $logfile
./limit_V2 $NOJET_CARD                                        >> $logfile
mv *.result.txt $WORK_DIR

cd $WORK_DIR

rm -rf CMSSW_6_1_1

