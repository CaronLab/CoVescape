#! /bin/bash



MutationList=$1
CoVpath=$2

pathToData=$CoVpath/scripts/MutationalPatterns/data


python3 $CoVpath/scripts/MutationalPatterns/V6_Compare_Country_mutations_withReplicates_2sampleTtest.py $MutationList $CoVpath/scripts/MutationalPatterns/data/NeutralEvolution_simulations.txt $pathToData


#withCore="_withCore"
#withCoreFile="_withCore.txt"
#folder=${netMHCpanOUT%.csv*}$withCore
#file=${netMHCpanOUT%.csv*}$withCoreFile

#mkdir $folder

#rm -f ./$folder/$file

#cd ./Results_withCore

#cp *withCore.csv ../$folder/

#ls *withCore.csv >> ../$folder/$file

#cd ..


