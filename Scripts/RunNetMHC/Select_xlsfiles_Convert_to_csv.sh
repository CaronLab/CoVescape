#! /bin/bash


OutputPATH=$1
CoVpath=$2

cd $OutputPATH

rm -f Out.txt

ls *.xls >> Out.txt

#cd $CoVpath/scripts/RunNetMHC/
python3 $CoVpath/scripts/RunNetMHC/SaveastextRead_EXCELFile_Pandas.py Out.txt $OutputPATH
