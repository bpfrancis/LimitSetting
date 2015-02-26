root -b -q -l 'makeTemplate_v3.C+'

cd datacards/
tar -czf datacards.tgz *.dat
mv datacards.tgz ../condor_submit/
cd ../
