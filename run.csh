#root -b -q -l 'makeTemplate.C+("bino","nojet")'
root -b -q -l 'makeTemplate.C+("bino","1jet")'

#root -b -q -l 'makeTemplate.C+("wino","nojet")'
#root -b -q -l 'makeTemplate.C+("wino","1jet")'

#root -b -q -l 'makeTemplate.C+("SMSScan", "nojet")'
#root -b -q -l 'makeTemplate.C+("SMSScan", "1jet")'

#root -b -q -l 'makeTemplate.C+("bino_mNScan","nojet")'
#root -b -q -l 'makeTemplate.C+("bino_mNScan","1jet")'

cd datacards/multiChannel/
tar -czf datacards.tgz *.dat
mv datacards.tgz ../../condor_submit/
cd ../../

