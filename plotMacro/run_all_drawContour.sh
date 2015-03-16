source /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scramv1 runtime -sh`

dest_dir=plots

root -b -q -l 'drawContour.C+'

mkdir -p $dest_dir
mv *.gif *.pdf $dest_dir

cd $dest_dir
tar -czf stuff.tgz *.*
scp stuff.tgz $HEP
rm stuff.tgz
