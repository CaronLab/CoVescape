
from Bio import SeqIO #have to install Biopython before using this.
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys


def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq)]: #get frame and nuc sequence from 5'-> 3' strand only (specific to SARS)
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table)) #translate using Biopython tool translate()
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start) # '*' is STOP CODON. We are constantly looking for the next stop codon.
                if aa_end == -1: # If cant find next stop codon, we've reached end of genome.
                    aa_end = trans_len
                temp = trans.find("M", aa_start) #Once have stop codon, look for corresponding 1st start codon.
                if temp != -1:
                    aa_start = temp
                if aa_end - aa_start >= min_protein_length: # Determine start and end of ORF in genetic code.
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    answer.append((start, end, strand, trans[aa_start:aa_end])) #Store ORF info in a list of ORFs
                aa_start = aa_end + 1
    answer.sort()
    return answer #Return list of ORFs.


def generate_new_seq(mutSeq, nuc):
    if nuc < 14: # if mutation is closer than 14 a.a. to 5'end (won't have 14 a.a. n that side)
        return mutSeq[:(nuc+14)]
    elif (len(mutSeq)-1)-nuc < 14: #same thing on 3' end.
        return mutSeq[(nuc-14):(len(mutSeq)-1)]
    else: #If mutation is in the middle of the protein, no worries!
        return mutSeq[(nuc-14):(nuc+15)] 

def Compare_Seqs(mutDict, refDict):
    newShortSeqs = {}
    for mutID, mutSeq in mutDict.items(): #iterate through ORFs of 'mutated' genome
        newSeqList = []
        for refID, refSeq in refDict.items(): #for each ORF of mutated genome, compare to all reference proteins (ORFs)
            tracker = 0
            for nuc in range(len(refSeq)-1): #Comparing potentially mutated ORF to reference ORF nuc by nuc.
                if mutSeq[:4] != refSeq[:4]: #Compare the first few residues to see if its the same protein.
                    break
                elif mutSeq[:4] == refSeq[:4] and len(mutSeq) != len(refSeq): #If it's the same protein but different lengths (deletion mutation or sequencing error), skip. Will have to change that eventually...
                    break
                elif mutSeq[nuc].upper() == 'X':
                    pass
                elif mutSeq[nuc] == refSeq[nuc]: #Its the same protein but no mutation.
                    pass
                else:
                    print('%s%s%s'% (refSeq[nuc], nuc+1, mutSeq[nuc]))
                    NewSeqMut = generate_new_seq(mutSeq, nuc) #calling 'generate_new_seq' function on Mut.
                    NewSeqRef = generate_new_seq(refSeq, nuc) #calling 'generate_new_seq' function on Ref.
                    mut_info = '%s_%s%s%s'% (refID, refSeq[nuc], nuc+1, mutSeq[nuc]) #generate handle of fasta.
                    newShortSeqs[mut_info] = (NewSeqMut, NewSeqRef)
    print(newShortSeqs)


def make_Bash_Script(fileName, folder):
    job_content = f'''#!/bin/bash
#SBATCH --time=0-12:00
#SBATCH --account=def-carone
#SBATCH --job-name={fileName}
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --mail-user=david.hamelin.1@umontreal.ca
#SBATCH --mail-type=FAIL

cd {folder}

{sys.argv[4]} -a {sys.argv[2]} -f {folder}/{fileName} -xls -xlsfile {fileName[:-4]}.xls -BA'''





    with open(f'{folder}/{fileName[:-4]}_job.sh', 'w') as f: #./
        f.write(job_content)


def Make_ComputeCanada_File(sortedMuts, folder):

    index = 0
    file_number = 1
    AllFIles_Dict = {}
    File_Specific_Dict = {}
    for mut_info, seqs in sortedMuts: #Split mutations (mutated/reference mutation pairs) in groups of 400. They will all be submitted as separate jobs
        if index < 400:
            File_Specific_Dict[mut_info] = seqs
            index += 1
        
        else:
            AllFIles_Dict[file_number] = File_Specific_Dict
            File_Specific_Dict = {}
            index = 0
            file_number += 1

    if index < 400:
        AllFIles_Dict[file_number] = File_Specific_Dict #If there are less than 400 unique mutations overall, take as one netMHCpan input file.



    for fileNum, FileSpecific_Dict in AllFIles_Dict.items():
        if fileNum < 10:
            fileNum = '0' + str(fileNum)
        fileName = '%s_Compute_Allresults_NUC.txt'% (fileNum)
        with open('%s/%s'% (folder, fileName), 'w') as File: #./
            for mut_info, seqs in FileSpecific_Dict.items():
                if int(seqs[0]) == 1:
                    pass
                else:
                    File.write('>%s_%s\n%s\n\n'% (seqs[4], seqs[3], seqs[1]))
                    File.write('>%s%s%s\n%s\n\n'% (seqs[3], 'Ref', seqs[4], seqs[2]))
        make_Bash_Script(fileName, folder)  
        #For each groups of mutations (containing 400 reference/mutated pairs), write out a txt file containing the
        #info in the following format:

        #>mutant1
            #mutated sequence
        #>reference1
            #reference sequence


##########################################


    #Aso print out a file with all mutations/sequences (provided to netMHCpan), for record-keeping.
    with open('%s/%s'% (folder, 'Compute_Allresults_NUC.txt'), 'w') as File: #./
        for mut_info, seqs in sortedMuts: #.items()
            if int(seqs[0]) == 1:
                File.write('>SINGLETON%s_%s\n%s\n\n'% (seqs[4], seqs[3], seqs[1]))
                File.write('>SINGLETON%s%s%s\n%s\n\n'% (seqs[3], 'Ref', seqs[4], seqs[2]))
            else:
                File.write('>%s_%s\n%s\n\n'% (seqs[4], seqs[3], seqs[1]))
                File.write('>%s%s%s\n%s\n\n'% (seqs[3], 'Ref', seqs[4], seqs[2]))



def Create_Fasta_File(AllMutInfo, aa_MutList):

    folder = sys.argv[3]  #Create folder name 
    os.makedirs(os.path.dirname('%s/'% (folder)), exist_ok=True) #create folder

    sortedMuts = sorted(AllMutInfo.items(), key = lambda x: x[1], reverse = True) #Sort mutation dictionary by mutation frequency
    print('This is sorted: %s'% (sortedMuts))

    Make_ComputeCanada_File(sortedMuts, folder) #Create netMHCpan input file(s)
     
    with open('%s/%s'% (folder, '_Allresults_NUC.txt'), 'w') as File:  #Create txt file with all results
        for mut_info, seqs in sortedMuts: #.items()
            if int(seqs[0]) == 1:
                File.write('>SINGLETON%s_%s\n%s\n\n'% (seqs[4], seqs[3], seqs[1]))
                File.write('>SINGLETON%s%s%s\n%s\n\n'% (seqs[3], 'Ref', seqs[4], seqs[2]))
            else:
                File.write('>%s_%s\n%s\n\n'% (seqs[4], seqs[3], seqs[1]))
                File.write('>%s%s%s\n%s\n\n'% (seqs[3], 'Ref', seqs[4], seqs[2]))


    

    MutRateGenomeID_Dataframe = pd.DataFrame.from_dict(aa_MutList, orient = 'index')
    MutRateGenomeID_Dataframe = MutRateGenomeID_Dataframe.rename(columns = {0 : 'MutID', 1 : 'Total_mutations', 2: 'Mutation_Rates'})

    MutRateGenomeID_Dataframe.sort_values(by = 'Total_mutations', ascending = False, inplace = True)

    MutRateGenomeID_Dataframe.to_csv('%s/%s'% (folder, '_MutFrequency_NUC.csv')) #GEnerate mutation file, which will be used in downstream analyses. This file contains all unique mutations, sorted by mutation frequency.



while True:
    RefProteins_aa = {}

    

    try:
        with open("./data/Reference_Seqs_noORF9.txt", "rU") as File: #open file containing reference protein sequences, store in  dictionary (key = protein name, item = protein sequence)
            for record in SeqIO.parse(File, "fasta"):
                RefProteins_aa[record.id] = record.seq


        with open('./data/Refsequence_Nuc.fasta', "rU") as File: #Open file containing all potentially mutated genomes.   
            for record in SeqIO.parse(File, "fasta"): #Parse through all genomes in the fasta file.
                RefProteins_nuc = {}
                sequence = record.seq
                
                table = 1 #Translation table. 1 = standard translation.
                min_pro_len = 37 #minimum protein length is set to 100 residues.
                
                if '-' in record.seq:
                    sequence = record.seq.tomutable()
                    while(sequence.count('-')):
                        sequence.remove('-')
                    sequence = sequence.toseq()
                else:
                    sequence = record.seq

                
                
                orf_list = find_orfs_with_trans(sequence, table, min_pro_len) #Submit genome fasta, codon table and min protein length to find_orfs_with_trans function to find all ORFs.
                for start, end, strand, pro in orf_list: 
                    
                    RefProteins_nuc['%s:%s'% (start, end)] = [pro, record.seq[start:end]] # Store each ORF in a dictionary, where key = genome location of ORF, item = ORF sequence.
                    print('This is RefProteins_nuc : %s\n\n\n\n'% (RefProteins_nuc))

            Not_Proteins = []
            for key1, item1 in RefProteins_nuc.items(): #iterate through ORFs (dict contains both nuc and aa sequences)
                for key2, item2 in RefProteins_aa.items(): #Iterate through reference sequences
                    if item1[0][:10] == item2[:10]: #compare first 10 residues of translated ORF vs each ref sequences
                        item1.append(key2) #If there is a match, we have identified the ORF (append protein name to file)
                        print('This is coordinates and name for this protein: %s, %s'% (key1, key2))
                if len(item1) != 3: #See if we found a protein match (If we did, item1 should have length 3)
                    Not_Proteins.append(key1)
                    print(key1)
            
            for key in Not_Proteins:
                RefProteins_nuc.pop(key, None)

            print('This is RefProteins_nuc_Final : %s\n\n\n\n'% (RefProteins_nuc))
    
    except FileNotFoundError:
        print('Wrong file or file path')
    else:
        break




while True:
    #######REF = input('\n\nPlease, enter Reference Genome\nEx. MyRefGenomes.txt\nwhen done, press Enter: ')


  

    try: 

        with open(sys.argv[1], "rU") as File:  #Import mutation list (tab delimited txt file)
            lines = File.readlines() #Read file
        count = 0
        listOfList = []
        for line in lines: #Here, we convert tab delimited file into a csv for handling by pandas and numpy (replace \t by commas, split strings into lists, convert to dataframe)
            print(len(line))
            newline = line.replace('\t', ',')    
            #newline = list(line.split(','))
                    
            print(newline)
            if count == 0:
                Columns = newline.split(',')
                count += 1
            else:
                listOfList.append(newline.split(','))
        array = np.array(listOfList)
        print('This is ListOfList: %s'% (listOfList))
        print('This is array: %s'% (array))

        AllMutations_Frame = pd.DataFrame(data = listOfList, columns = ['Isolate', 'Position', 'N_Alleles', 'N_Genomes', 'Allele1:Freq', 'Allele2:Freq', 'Allele3:Freq', 'Allele4:Freq']) #['Isolate', 'Position', 'N_Alleles', 'N_Genomes', 'Allele:Freq'] , columns = Columns

        print(AllMutations_Frame)
        AllMutations_Frame.to_csv('JCMutations.csv')


        ####Manipulate file to calculate number of individuals carrying mutation, and reformat for further handling
        AllMutations_Frame.set_index(['Isolate', 'Position', 'N_Alleles', 'N_Genomes', 'Allele1:Freq'], inplace = True)
        AllMutations_Frame = AllMutations_Frame.stack().to_frame().reset_index().rename(columns = {'level_5': 'Allele', 0: 'Mutation'})
        
        print(AllMutations_Frame)

        AllMutations_Frame.insert(loc = 7, column = 'Formatted_Mutation', value = AllMutations_Frame.apply(lambda x: '%s%s%s'% (x['Allele1:Freq'][0], x['Position'], x['Mutation'][0]), axis = 1))

        AllMutations_Frame.insert(loc = 8, column = 'N_Mutants', value = AllMutations_Frame.apply(lambda x: round( float(x['N_Genomes']) * float(x['Mutation'][ x['Mutation'].find(':')+1 : ]) ), axis = 1))

        AllMutations_Frame.insert(loc = 9, column = 'Mutant_Freq', value = AllMutations_Frame['Mutation'].apply(lambda x: float(x[ x.find(':')+1 : ]) ))
        

        ###################
        AllMutations_Frame = AllMutations_Frame[AllMutations_Frame['N_Mutants'] > 5000]
        ###################

        #AllMutations_Frame.to_csv('JCMutations_Stacked.csv')

    except FileNotFoundError:
        print('Wrong file or file path')
    else:
        break


################
#Here, for each mutation, we will 'synthesize the mutant 'in silico'.
################

AllMutInfo = {}
aa_Mut_count = 0
AllMut_Count = 0
aa_MutList = {}
for loc, seqs in RefProteins_nuc.items(): #For each protein, acquire the start and the end of the protein. This will help select mutations that are specific to a protein
    start = int( loc[:loc.find(':')] ) + 1
    end = int( loc[loc.find(':')+1 :] )

    Protein_specific_Mutations_Frame = AllMutations_Frame[ ( AllMutations_Frame['Position'].astype(int) >= ( start ) ) & ( AllMutations_Frame['Position'].astype(int) <= ( end ) ) ] ##Do plus one b/c here, start and end are calculated in 0-9 style. So need to ad one to have 1-10 style.
    print('This is start, end and Protein_specific_Mutations: \n%s\n%s\n%s'% (start, end, Protein_specific_Mutations_Frame))


    ######Exctract numpy arrays for various mutation info
    Protein_specific_Mutations_numpy = Protein_specific_Mutations_Frame['Formatted_Mutation'].values
    Protein_specific_genomeCount_numpy = Protein_specific_Mutations_Frame['N_Mutants'].values
    Protein_specific_MutFreq_numpy = Protein_specific_Mutations_Frame['Mutant_Freq'].values

    Protein_Specific_count = 0

    for mutation in Protein_specific_Mutations_numpy:
        protein_name = seqs[2] #acquire protein
        position = int(mutation[1:-1]) - start #acquire mutation position within the protein
        print('seq length: %s, Position: %s'% (len(seqs[1]), position))
        mutated_sequence = seqs[1].tomutable()
        mutated_sequence[position] = mutation[-1] #Create a mutated sequence. Remove 1 from position b/c position is in 1-10 format, need 0-9
        mutated_sequence = mutated_sequence.toseq()
        Reference_sequence = seqs[0]
        
        print('Ref_nuc --> Mut_nuc: %s --> %s, mutation = %s'% (seqs[1][position], mutated_sequence[position], mutation))

        orf_list = find_orfs_with_trans(mutated_sequence, table, min_pro_len) #Translate mutated nuc sequence
        for start2, end2, strand, pro in orf_list:
            if len(pro) ==  len(Reference_sequence):       
                
                count = 0
                mut_aa_pos = 0
                for nuc in pro:
                    if pro[count] != Reference_sequence[count]: #Determine the amino acid mutation resulting from the nucleotide mutation event
                        mut_aa_pos = count
                        ref_aa = Reference_sequence[count]
                        mut_aa = pro[count]
                        print('Ref_aa --> Mut_aa: %s --> %s'% (ref_aa, mut_aa))
                    count += 1

                if mut_aa_pos != 0:
                    mutated_short_seq = generate_new_seq(pro, mut_aa_pos) #Generate mutated 'short sequence windows' (netMHCpan mutated input)
                    reference_short_seq = generate_new_seq(Reference_sequence, mut_aa_pos) #Generate reference 'short sequence windows' (netMHCpan reference input)
                    #mut_info = [protein_name, ref_aa, mut_aa_pos, mut_aa, aa_Mut_count]
                    mut_info = '%s_%s%s%s'% (protein_name, ref_aa, mut_aa_pos+1, mut_aa) #Generate mutation info
                    AAA = 0
                    if mut_info not in AllMutInfo.keys():
                        AAA += 1
                        AllMutInfo[mut_info] = [Protein_specific_genomeCount_numpy[Protein_Specific_count], mutated_short_seq, reference_short_seq, aa_Mut_count, mut_info]
                    if AAA == 0:
                        print('WE HAVE A DUPLICATE!!!: %s'% (mut_info))
                    aa_MutList['%s_%s'% (protein_name, mutation)] = [ '%s_%s%s%s'% (protein_name, ref_aa, mut_aa_pos+1, mut_aa), Protein_specific_genomeCount_numpy[Protein_Specific_count], Protein_specific_MutFreq_numpy[Protein_Specific_count]  ]

                    aa_Mut_count += 1
        Protein_Specific_count += 1
        
        AllMut_Count += 1




Create_Fasta_File(AllMutInfo, aa_MutList) #Generate the mutation files as well as appropriate netMHCpan input file for downstream analysis




