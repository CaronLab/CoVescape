#! /bin/bash

#Enter code; raw netMHCpan 4.1 file; protein fasta sequences; mutation file; list of validated epitopes; HLA-specific gain-loss files.

MutationList=$1
coreOut=$2
CoVpath=$3

pathToData=$CoVpath/scripts/EpiMutDist/data


python3 $CoVpath/scripts/EpiMutDist/V7e_findMatchedEpitopes_epitope_Hotspots_WithMcKay_Cleanedup_Final.py $pathToData/AllProteins_netMHCpan4_1_noORF9.csv $pathToData/Reference_Seqs.txt $MutationList $pathToData/McKay_Validated_Epitopes.csv $coreOut