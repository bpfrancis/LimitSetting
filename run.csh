root -b -q -l 'makeTemplate.C+\(\"allChannels\"\)'
root -b -q -l 'makeTemplate.C+\(\"SR1\"\)'
root -b -q -l 'makeTemplate.C+\(\"SR2\"\)'
root -b -q -l 'makeTemplate.C+\(\"ele\"\)'
root -b -q -l 'makeTemplate.C+\(\"muon\"\)'

cd datacards/
tar -czf datacards.tgz *.dat
mv datacards.tgz ../condor_submit/
cd ../
