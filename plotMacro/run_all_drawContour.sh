source /cvmfs/cms.cern.ch/cmsset_default.csh
eval `scramv1 runtime -sh`

version=$1

cp table/stop-bino_$version.table table/stop-bino.table

dest_dir=plots

root -b -q -l 'drawContour.C+'

mkdir -p $dest_dir
mv *.gif *.pdf $dest_dir

cd $dest_dir
tar -czf stuff.tgz *.*
scp stuff.tgz $HEP
rm stuff.tgz
cd ..

rm table/stop-bino.table
