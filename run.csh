root -b -q -l 'makeTemplate_v2.C+'

cd datacards/
tar -czf datacards.tgz *.dat
mv datacards.tgz ../condor_submit/
cd ../
