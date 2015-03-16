root -b -q -l 'makeTemplate.C+'

cd datacards/
tar -czf datacards.tgz *.dat
mv datacards.tgz ../condor_submit/
cd ../
