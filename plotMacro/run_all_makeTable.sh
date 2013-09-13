
dest_dir="table/multiChannel"
mkdir -p $dest_dir

rm -f *.table

./makeTable.sh bino jet multiChannel
#./makeTable.sh bino nojet multiChannel
#./makeTable.sh wino jet multiChannel
#./makeTable.sh wino nojet multiChannel
#./makeTable.sh bino_mNScan 1jet multiChannel
#./makeTable.sh bino_mNScan nojet multiChannel
mv *.table $dest_dir
