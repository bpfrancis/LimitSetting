export VO_CMS_SW_DIR=/uscmst1/prod/sw/cms
source $VO_CMS_SW_DIR/cmsset_default.sh
eval `scramv1 runtime -sh`

dest_dir=plots

#root -b -q -l 'drawContour.C+("sms_gg_1jet",true)'
#root -b -q -l 'drawContour.C+("sms_gg_nojet",true)'
root -b -q -l 'drawContour.C("bino","_jet",true)'
#root -b -q -l 'drawContour.C("bino","_nojet",true)'
#root -b -q -l 'drawContour.C("wino","_jet",true)'
#root -b -q -l 'drawContour.C("wino","_nojet",true)'
#root -q -l 'drawContour.C("bino_mNScan","jet",true)'
#root -q -l 'drawContour.C("bino_mNScan","_btag0",true)'

mkdir -p $dest_dir
mv *.gif *.pdf $dest_dir
