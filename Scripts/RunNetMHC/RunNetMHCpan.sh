#!/bin/bash


echo "$#"
muts=$1
echo "$muts"
HLA=$2
echo "$HLA"
OutputPATH=$3
echo "$OutputPATH"
netMHC=$4
echo "$netMHC"
CoVpath=$5
echo "This is CoVescape directory: ${CoVpath}"



CurrentDIR=$(pwd)
echo "This is current dir $CurrentDIR"

cd $CoVpath/scripts/RunNetMHC/

python3 V4Acquire_Juliesmutations_OUTPUTCOMPUTECANADA_NoDuplicates.py $CurrentDIR/$muts $HLA $OutputPATH $netMHC 

cd $OutputPATH/
chmod -R 755 $OutputPATH/
$CoVpath/scripts/RunNetMHC/RunAll_Jobs.sh

$CoVpath/scripts/RunNetMHC/Select_xlsfiles_Convert_to_csv.sh $OutputPATH $CoVpath

