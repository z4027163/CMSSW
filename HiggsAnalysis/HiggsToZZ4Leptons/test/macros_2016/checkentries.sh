#!/bin/bash


# for file in `ls /lustre/cms/store/user/dburns/MonoHiggs/Spring16_merged/*Baryon*`; do
for file in `ls /lustre/cms/store/user/gminiell/MonoHiggs/Spring16_merged/*2HDM*`; do
# for file in `ls /lustre/cms/store/user/gminiell/MonoHiggs/Data2016_MonoHiggs_13TeV_merged/*`; do

echo $file
cat checkentries.C | sed "s?filename?${file}?g" > tmp.C
g++ -I $ROOTSYS/include tmp.C `root-config --glibs` `root-config --libs` `root-config --cflags` -o checkentries
./checkentries

done

