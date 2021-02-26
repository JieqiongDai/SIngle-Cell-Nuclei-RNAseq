#!/bin/bash

## need to change directory,projects,flowcells and sample directory for each run

touch config.yaml
echo "directory: 'diem' ">> config.yaml

echo "projects:">> config.yaml

echo " - diem" >> config.yaml

echo "samples:">> config.yaml

for i in $(ls /data/CGR_10X/raw_data/snRNA/200929_A00423_0094_BHLM2GDMXX/200929_A00423_0094_BHLM2GDMXX-simple-1.0.0.csv)
do
cat ${i} | sed "1 d" | awk -F "," '{print " -",$2}'|sort |uniq >>config.yaml
done

echo "flowcells: 'diem' ">> config.yaml
