#! /bin/bash



netMHCpanOUT=$1
tempFrame=$2
CoVpath=$3


python3 $CoVpath/scripts/GetCore/V4_WithDictionariesWithMultiprocessingIncludeCoreSeqs_MorePostProcessingAnalysis_CommonCore_patterns.py $netMHCpanOUT $tempFrame


withCore="_withCore"
withCoreFile="_withCore.txt"
folder=${netMHCpanOUT%.csv*}$withCore
file=${netMHCpanOUT%.csv*}$withCoreFile

mkdir $folder

rm -f ./$folder/$file

cd ./Results_withCore

cp *withCore.csv ../$folder/

ls *withCore.csv >> ../$folder/$file

cd ..


