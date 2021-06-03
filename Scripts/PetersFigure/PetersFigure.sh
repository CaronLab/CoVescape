#! /bin/bash

#Enter code; raw netMHCpan 4.1 file; protein fasta sequences; mutation file; list of validated epitopes; HLA-specific gain-loss files.


CoreOut=$1
CoVpath=$2

python3 $CoVpath/scripts/PetersFigure/V4_Peters_ttest_DIfferent_coloring.py $CoreOut