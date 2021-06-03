#! /bin/bash




netMHCpanOUT=$1
mutations=$2
CoVpath=$3

pathToData=$CoVpath/scripts/PostNetMHCPan/data



CurrentDIR=$(pwd)



python3 $CoVpath/scripts/PostNetMHCPan/V8_nMp_netMHCpan41_DataAnalysis_Faster_Cleaned_up_Final_VersionWith_MULTIPROCESSINGandSystem.py $netMHCpanOUT $mutations y y $pathToData

DataOutput="_DataOutput"
tempFrame="_tempFrame.txt"
folder=${netMHCpanOUT%.csv*}$DataOutput
file=${netMHCpanOUT%.csv*}$tempFrame

mkdir $folder

rm -f ./$folder/$file

cd ./${netMHCpanOUT%.csv*}

cp *temp_Frame* ../$folder/

ls *temp_Frame* >> ../$folder/$file

cd ..


