import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import os
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
import math
from matplotlib_venn import venn2
from multiprocessing import Process
import time
import sys


aligner = Align.PairwiseAligner()
sns.set_style(style = 'white')

def LetsAlign(x, y):
    if len(x) == len(y):
        aligned = aligner.align(x, y)
        return (aligned.score/len(y))

def Alignment(x, newFrame):
    
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -5
    
    newFrame.insert(loc = 0, column = 'alignScore', value = newFrame['Epitope.2'].apply(lambda y: LetsAlign(x, y)))
    print(newFrame)
    maxScore = newFrame['alignScore'].max()
    if math.isnan(maxScore) == False :
        corresPeptide = newFrame.iloc[newFrame['alignScore'].idxmax(), 1]
        binding_nM = newFrame.iloc[newFrame['alignScore'].idxmax(), 3]
        
        binding_Qualitative = ''
        if type(newFrame.iloc[newFrame['alignScore'].idxmax(), 2]) == str:
            if 'positive' in newFrame.iloc[newFrame['alignScore'].idxmax(), 2].lower():
                binding_Qualitative = 'pos'
            elif 'negative' in newFrame.iloc[newFrame['alignScore'].idxmax(), 2].lower():
                binding_Qualitative = 'neg'

        
        Tcell_Qualitative = ''
        if type(newFrame.iloc[newFrame['alignScore'].idxmax(), 5]) == str:
            if 'positive' in newFrame.iloc[newFrame['alignScore'].idxmax(), 5].lower():
                Tcell_Qualitative = 'pos'
            elif 'negative' in newFrame.iloc[newFrame['alignScore'].idxmax(), 5].lower():
                Tcell_Qualitative = 'neg'
        

        newFrame.drop(['alignScore'], axis = 1, inplace = True)
        print('results are %s, %s, %s, %s and %s'% (maxScore, corresPeptide, binding_nM, binding_Qualitative, Tcell_Qualitative))
        if maxScore >= 3.5:
            
            return (maxScore, corresPeptide, binding_nM, binding_Qualitative, Tcell_Qualitative)
    else:
        newFrame.drop(['alignScore'], axis = 1, inplace = True)



def Import_bindingIEDB():
    while True:
        prompt = f'{sys.argv[5]}/SARSCOV_posneg_IEDB.csv'
        try:
            bindingIEDB = pd.read_csv(prompt)
           
            
            return bindingIEDB

        except FileNotFoundError:
            print('file not found. Make sure to copy in in same path as this script.')
        else:
            break
def Import_tcellIEDB():
    while True:
        prompt = f'{sys.argv[5]}/TCELL_posneg_IEDB.csv'
        try:
            tcellIEDB = pd.read_csv(prompt)
           
            
            return tcellIEDB

        except FileNotFoundError:
            print('file not found. Make sure to copy in in same path as this script.')
        else:
            break


def PairwiseAlignment_Tcell_MHCbinding(bindingIEDB, tcellIEDB, temp_Frame, FileName):
    bindingIEDBFrame = bindingIEDB[['Epitope.2', 'Assay.4', 'Assay.6', 'MHC']].copy()
    tcellIEDBFrame = tcellIEDB[['Epitope.2', 'Assay.4', 'MHC']].copy()
    tcellIEDBFrame.rename(columns = {'Assay.4': 'Assay_Tcell'}, inplace = True)
    IEDBFrame = pd.merge(bindingIEDBFrame, tcellIEDBFrame, on = ['Epitope.2', 'MHC'],  how = 'outer', sort = False)

    print(IEDBFrame)
    IEDBFrame.to_csv('Binding_Tcell_combined_frames.csv')

    IEDBFrame = IEDBFrame.apply(lambda x: x.str.replace('*', ''))
    IEDBFrame = IEDBFrame.loc[IEDBFrame['MHC'] == temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')]]


    if IEDBFrame.empty == False:
        print('\n\nwere here%s\n\n'% (temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')]))
        IEDBFrame = IEDBFrame.drop(IEDBFrame.index[0])
        IEDBFrame.reset_index(drop = True, inplace = True)
        print('%s, %s'% (IEDBFrame, type(IEDBFrame)))
                    
        temp_Frame.insert(loc = 3, column = 'SARS_CoV_Match', value = temp_Frame['Peptide_Ref'].apply(lambda x: Alignment(x, IEDBFrame)))
        temp = pd.DataFrame()
        print(temp_Frame['SARS_CoV_Match'].dropna().tolist())
        if temp_Frame['SARS_CoV_Match'].dropna().tolist():
            temp_Frame.sort_values(by = 'SARS_CoV_Match', ascending = False, inplace = True)
            print(temp_Frame)
            temp[['AlignScore', 'SARS_CoV_Pep', 'binding_nM', 'binding_Qualitative', 'Tcell_Qualitative']] = pd.DataFrame(temp_Frame['SARS_CoV_Match'].tolist(), index = temp_Frame.index)
            

            temp_Frame.drop('SARS_CoV_Match', axis = 1, inplace = True)
            temp_Frame.insert(loc = 5, column = 'AlignScore', value = temp['AlignScore'])
            temp_Frame.insert(loc = 6, column = 'SARS_CoV_Pep', value = temp['SARS_CoV_Pep'])
            temp_Frame.insert(loc = 7, column = 'binding_nM', value = temp['binding_nM'])
            temp_Frame.insert(loc = 8, column = 'binding_Qualitative', value = temp['binding_Qualitative'])
            temp_Frame.insert(loc = 9, column = 'Tcell_Qualitative', value = temp['Tcell_Qualitative'])
            temp_Frame = temp_Frame.round({'AlignScore': 2})
            temp_Frame.sort_index(axis = 0, ascending = True, inplace = True)
            print('GoodToGo!!!')

            
            AlignmentFrame = temp_Frame[['Peptide_Ref' ,'Peptide_Mut','SARS_CoV_Pep', 'nM_Ref', 'nM_Mut', 'binding_nM', 'binding_Qualitative', 'Tcell_Qualitative']].copy()
            AlignmentFrame.rename(columns = {'nM_Ref': 'Predicted Ref nM', 'nM_Mut': 'Predicted Mut nM', 'binding_nM': 'Measured nM'}, inplace = True)

            AlignmentFrame.dropna(subset = ['SARS_CoV_Pep'],inplace = True)
            AlignmentFrame.reset_index(drop = True)
            print(AlignmentFrame)
            
            print(AlignmentFrame)
            AlignmentFrame.to_html('alignment_table.html')

            BarGraphFrame = AlignmentFrame[['Peptide_Ref', 'Predicted Ref nM', 'Measured nM', 'Predicted Mut nM']].copy().set_index('Peptide_Ref')

            BoxplotFrame = BarGraphFrame

            BarGraphFrame = BarGraphFrame.stack().to_frame().reset_index().rename(columns = {'level_1': 'data_source', 0: 'binding_affinity'})

            print(BarGraphFrame)
            print(type(BarGraphFrame))
            BarGraphFrame['binding_affinity'] = BarGraphFrame['binding_affinity'].apply(lambda x: np.log2(float(x)) if math.isnan(float(x)) == False else x)
            print('BarGraphFrame: \n%s'% (BarGraphFrame))


            BoxplotFrame['measured vs Reference'] = np.log2(BoxplotFrame['Predicted Ref nM'].astype(float) / BoxplotFrame['Measured nM'].astype(float)) #Calculate log2 fold change
            BoxplotFrame['Reference vs mutated'] = np.log2(BoxplotFrame['Predicted Mut nM'].astype(float) / BoxplotFrame['Predicted Ref nM'].astype(float))
            BoxplotFrame['measured vs mutated'] = np.log2(BoxplotFrame['Predicted Mut nM'].astype(float) / BoxplotFrame['Measured nM'].astype(float))
            BoxplotFrame.drop(['Predicted Ref nM', 'Measured nM', 'Predicted Mut nM'], axis = 1, inplace = True)
            BoxplotFrame = BoxplotFrame.stack().to_frame().reset_index().rename(columns = {'level_1': 'transition', 0: 'fold change'})

            print('BoxplotFrame: \n%s'%(BoxplotFrame))

            with sns.axes_style("ticks"):
                fig, ax = plt.subplots(figsize=(6, 3))
                sns.set(font_scale=0.6)
                sns.despine()
                sns.boxplot(x = 'transition', y = 'fold change', data = BoxplotFrame, ax = ax)
                sns.swarmplot(x = 'transition', y = 'fold change', data = BoxplotFrame, color = '.2', ax = ax)
                ax.set_ylabel('Fold Change')
                ax.set_title('comparison of predicted binding affinity (nM) of SARS-CoV-2 reference and mutated peptides to %s\n to corresponding SARS-CoV measured binding affinity'% (temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')]))
                plt.tight_layout()
                plt.savefig('%s/%s_FoldChange_SARS_IEDB_%s'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))

            with sns.axes_style("ticks"):
                plt.subplots(figsize=(50,5))
                sns.set(font_scale=0.6)
                sns.despine()
                fig = sns.catplot(x = 'Peptide_Ref', y = 'binding_affinity', hue = 'data_source', data = BarGraphFrame, kind = 'bar', height = 5, aspect = 2.5, palette = 'muted', legend = False)
                fig.set_axis_labels('SARS-CoV-2 peptides','log2 binding affinity (nM)')
                fig.set_xticklabels(rotation = 45, horizontalalignment = 'right')
                plt.title('comparison of predicted binding affinity (nM) of SARS-CoV-2 reference and mutated peptides to %s\n to corresponding SARS-CoV measured binding affinity'% (temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')]))
                fig.despine(left=True)
                plt.legend(loc = 'upper right')
                plt.grid()
                plt.tight_layout()
                plt.savefig('%s/%s_predicted_REF_MUT_nM_vs_MEASURED_%s'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))


        else:
            temp_Frame.drop('SARS_CoV_Match', axis = 1, inplace = True)
            temp_Frame.insert(loc = 5, column = 'AlignScore', value = np.nan)
            temp_Frame.insert(loc = 6, column = 'SARS_CoV_Pep', value = np.nan)
            temp_Frame.insert(loc = 7, column = 'binding_nM', value = np.nan)
            temp_Frame.insert(loc = 8, column = 'binding_Qualitative', value = np.nan)
            temp_Frame.insert(loc = 9, column = 'Tcell_Qualitative', value = np.nan)
            print('GoodToGo!!!')
    else:
        temp_Frame.insert(loc = 5, column = 'AlignScore', value = np.nan)
        temp_Frame.insert(loc = 6, column = 'SARS_CoV_Pep', value = np.nan)
        temp_Frame.insert(loc = 7, column = 'binding_nM', value = np.nan)
        temp_Frame.insert(loc = 8, column = 'binding_Qualitative', value = np.nan)
        temp_Frame.insert(loc = 9, column = 'Tcell_Qualitative', value = np.nan)
        print('No SARS-CoV peptides were validated for this HLA type!!')

    return (temp_Frame, IEDBFrame)


'''
11. inputs: two peptide sequences, one mutated, one not mutation (reference)
outputs: index of mutation in peptide.
'''

def compare_seqs(Mut, Ref):
    print('THis is Mut and Ref: %s, %s'% (Mut, Ref))
    for i in range(len(Mut)):
        if Mut[i] != Ref[i]:
            if i <= 3:
                return '%s'% (i+1)
            elif i >= (len(Mut) - 4):
                return '%s'% (i - len(Mut))
            else:
                return '(...)'

def Effect_on_Binding_NonConservative(Mut, Ref):
    if (float(Mut) <= 2.0) and (float(Ref) >= 2.0):
        return 'gain'
    elif (float(Mut) >= 2.0) and (float(Ref) <= 2.0):
        return 'loss'
    else:
        return 'no effect'

def Effect_on_Binding_SemiConservative(Mut, Ref):
    if (float(Mut) <= 1.0) and (float(Ref) >= 2.0):
        return 'gain'
    elif (float(Mut) >= 2.0) and (float(Ref) <= 1.0):
        return 'loss'
    else:
        return 'no effect'

def Effect_on_Binding(Mut, Ref):
    if (float(Mut) <= 0.5) and (float(Ref) >= 2.0):
        return 'gain'
    elif (float(Mut) >= 2.0) and (float(Ref) <= 0.5):
        return 'loss'
    else:
        return 'no effect'


def compare_seqs_ReturnMut(Mut, Ref):
    for i in range(len(Mut)):
        if Mut[i] != Ref[i]:
            return '%s-->%s'% (Ref[i], Mut[i])

def mutations_per_peptides(Mut, Ref):
    num = 0
    for i in range(len(Mut)):
        if Mut[i] != Ref[i]:
            num += 1
    return num

'''
10. Function that reads .txt file ouput from first program (contains list of unique mutations, as well as their respective mutation rates)
Outputs dictionary with key = unique mutations, values = corresponding mutation rates.

'''

def MutationRate():
    while True:
        prompt = sys.argv[2]
        Nucleotide_positions = sys.argv[3]

        try:
            MutRate = pd.read_csv(prompt, index_col = 'MutID')
            rateArray = MutRate['Mutation_Rates'].to_numpy()
            MutRate.drop('Mutation_Rates', axis = 1, inplace = True)
            MutRate.insert(loc = 1, column = 'Mutation_Rates', value = np.round(rateArray, 2))

            if Nucleotide_positions == 'y':
                MutRate.drop('Unnamed: 0', axis = 1, inplace = True)
                MutRate.reset_index(inplace = True)

                MutRate_NoDuplicates = MutRate.groupby(['MutID']).agg({'Mutation_Rates': 'sum', 'Total_mutations': 'sum'})
                MutRate_NoDuplicates.sort_values(by = 'Total_mutations', ascending = False, inplace = True)

                print(MutRate_NoDuplicates)

                return MutRate_NoDuplicates

            else:
                return MutRate
        except FileNotFoundError:
            print('file not found. Make sure to copy in in same path as this script.')
        else:
            break



'''
9. Function to print out Figures using Seaborn package.
inputs: numerous parameters allowing to create different figures fitting each figure type.
Outputs: Figures saved in folder named by string saved in FileName variable.
'''

def printFigures(X, Y_FC, Y_diff, Data, Xlabel, Ylabel, Title, FigureName, SubPlot, FileName, Hue, Kind):
    Order = ['1', '2', '3', '4', '(...)', '-4', '-3', '-2', '-1']
    MutProteinList = Data.Mut_Protein.unique().tolist()
    MutProteinList = MutProteinList.sort()
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.set(font_scale=1)
    if len(Hue) != 0 and len(Kind) != 0: #If both hue (categorical coloring) and figure type were specified.
        if X == 'Mut_Position':
            sns.swarmplot(x = X, y = Y_FC, hue = Hue, data = Data, order = Order, ax = ax) #swarm plot for fold change (on bottom or subplot)
        elif X == 'Mut_Protein':
            sns.swarmplot(x = X, y = Y_FC, hue = Hue, data = Data, order = MutProteinList, ax = ax)

        else:
            sns.swarmplot(x = X, y = Y_FC, hue = Hue, data = Data, ax = ax)
        ax.set_xlabel(Xlabel)
        ax.set_ylabel('Fold Change')
        ax.axhline(y = 0, color = 'black', linestyle = '--') #add line at y = 0, to clarify which changes are negative vs positive.
        
    elif len(Hue) == 0 and len(Kind) != 0:
        if SubPlot == 'y':
            sns.swarmplot(x = X, y = Y_FC, data = Data, ax = ax)
            ax.set_xlabel(Xlabel)
            ax.set_ylabel('Fold Change')
            ax.axhline(y = 0, color = 'black', linestyle = '--')
        else:
            fig = sns.catplot(x = X, y = Y_FC, kind = Kind, data = Data, color = 'green', legend = False)

    ax.set_title('%s'% (Title[0:Title.find('_')]), fontsize = 10)
    plt.tight_layout()
    plt.savefig('%s/%s'% (FileName, FigureName))
    plt.clf()

def Make_Mut_Position_subplot(effect_on_binding_nonConservative, effect_on_binding, Order, temp_Frame, TM, FileName):
    
    with sns.axes_style("ticks"):
        fig, axes = plt.subplots(3, 1, sharex = True, figsize=(14, 12)) # , gridspec_kw={'width_ratios':[1,0.02]}, figsize=(14, 12)
        sns.set(font_scale=1.2)
        sns.despine()
        
        sns.swarmplot(x = 'Mut_Position', y = temp_Frame.columns[-3], hue = 'Reoccuring', data = temp_Frame.sort_values(by = 'PepLength', ascending = True), order = Order, ax = axes[0])
        sns.barplot(x = 'Mut_Position', y = 'count', hue = 'effect on binding', data = effect_on_binding_nonConservative, order = Order, hue_order = ['gain', 'loss', 'no effect'], palette = 'muted', ax = axes[1])
        sns.barplot(x = 'Mut_Position', y = 'count', hue = 'effect on binding', data = effect_on_binding, order = Order, hue_order = ['gain', 'loss', 'no effect'], palette = 'muted', ax = axes[2])
        axes[1].set_xlabel('Peptide Position')
        axes[0].set_xlabel('')
        axes[1].set_xlabel('')
        axes[0].set_ylabel('difference in netMHCpan Rank\n-(Mutated %Rank - Reference %Rank')
        axes[1].set_ylabel('Number of unique mutations')
        axes[2].set_ylabel('Number of unique mutations')
        axes[0].set_title(('Position-specific analysis of all unique HLA-associated polymorphisms \nwithin netMHCpan-predicted %s epitopes')% (temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')]), fontsize = 15)
        axes[0].axhline(y = 0, color = 'black', linestyle = '--')
        plt.grid()
        plt.tight_layout()
        plt.savefig('%s/%s_Mut_Position_DOUBLESubplot_%s'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))
        plt.clf()

    

def Normalize_By_Proteinlength(count, MutProtein, protein_list):
    return count / protein_list[MutProtein]

def Protein_Analysis(temp_Frame, MutRate, FileName):
    
    hue_palette = {'0-15%': 'blue', '15-30%': 'green', '30-50%': 'orange', '50-70%': 'red', '70% +': 'purple'}

    
    Protein_lengths = {'ORF1b': 2595, 'ORF1a': 4405, 'ORF3a': 275, 'ORF7a': 121, 'ORF8': 121, 'M': 222, 'N': 419, 'S': 1273, 'ORF6': 61, 'ORF10': 38, 'ORF9': 419, 'E': 75}

    All_mutations_frame = MutRate.reset_index()
    All_mutations_frame.insert(loc = 1, column = 'Mut_Protein', value = All_mutations_frame['MutID'].apply(lambda x: x[0:(x.find('_'))]))
    
    
    effect_on_binding_byProtein_nonConservative = temp_Frame.groupby(['Mut_Protein', 'effect on binding_nonConservative'])['level_0'].nunique().reset_index(name = 'count')
    effect_on_binding_byProtein_nonConservative.rename(columns = {'effect on binding_nonConservative': 'effect on binding'}, inplace = True)
    effect_on_binding_byProtein = temp_Frame.drop(temp_Frame[temp_Frame['effect on binding'] == 'no effect'].index)
    effect_on_binding_byProtein_non_normalized = effect_on_binding_byProtein.groupby(['Mut_Protein', 'effect on binding'])['level_0'].nunique().reset_index(name = 'count')

    effect_on_binding_byProtein = effect_on_binding_byProtein.groupby(['Mut_Protein', 'effect on binding'])['level_0'].nunique().reset_index(name = 'count')

    effect_on_binding_byProtein.insert(loc = 2, column = 'Normalized count', value = effect_on_binding_byProtein.apply(lambda x: Normalize_By_Proteinlength(count = x['count'], MutProtein = x['Mut_Protein'], protein_list = Protein_lengths), axis = 1))
    

    print('This is effect_on_binding_byProtein \n\n\n\n\n\n\n\n\n%s\n\n\n\n\n\n\n\n'% (effect_on_binding_byProtein))

    effect_on_binding_byProtein.drop(['count'], axis = 1, inplace = True)



    MutProteinList = All_mutations_frame.Mut_Protein.unique().tolist()
    MutProteinList.sort()
    print('This is MutProteinList!!! \n\n\n\n%s\n\n\n'% (MutProteinList))

    with sns.axes_style("ticks"):
        fig, axes = plt.subplots(3, 1, sharex = True, figsize=(14, 12)) 
        sns.set(font_scale=1.5)
        sns.despine()
        sns.swarmplot(x = 'Mut_Protein', y = temp_Frame.columns[-3], hue = 'Reoccuring', data = temp_Frame, order = MutProteinList, palette = hue_palette, ax = axes[0])
        sns.barplot(x = 'Mut_Protein', y = 'count', hue = 'effect on binding', data = effect_on_binding_byProtein_nonConservative, order = MutProteinList, hue_order = ['gain', 'loss', 'no effect'], palette = 'muted', ax = axes[1])
        sns.barplot(x = 'Mut_Protein', y = 'Normalized count', hue = 'effect on binding', data = effect_on_binding_byProtein, order = MutProteinList, hue_order = ['gain', 'loss', 'no effect'], palette = 'muted', ax = axes[2])
        axes[0].set_xlabel('')
        axes[1].set_xlabel('')
        axes[0].set_ylabel('difference in netMHCpan Rank\n-(Mutated %Rank - Reference %Rank')
        axes[1].set_ylabel('Number of mutations')
        axes[2].set_ylabel('Number of mutations')
        axes[0].axhline(y = 0, color = 'black', linestyle = '--')
        plt.grid()
        plt.tight_layout()
        plt.savefig('%s/%s_PROTEIN_SPECIFIC_Doublesubplot_%s'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))
        plt.clf()

    with sns.axes_style("ticks"):
        fig, axes = plt.subplots(figsize=(10, 6)) 
        sns.set(font_scale=1.5)
        sns.despine()
        
        sns.barplot(x = 'Mut_Protein', y = 'count', hue = 'effect on binding', data = effect_on_binding_byProtein_non_normalized, order = MutProteinList, hue_order = ['gain', 'loss'], palette = 'muted', ax = axes)
        axes.set_xlabel('')
        axes.set_ylabel('Number of mutations')
        plt.grid()
        plt.tight_layout()
        plt.savefig('%s/%s_PROTEIN_SPECIFIC_Only_Conservative_NONNORMALIZED%s'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))
        plt.clf() 

    with sns.axes_style("ticks"):
        fig, axes = plt.subplots(figsize=(10, 6)) 
        sns.set(font_scale=1.5)
        sns.despine()
        
        sns.barplot(x = 'Mut_Protein', y = 'Normalized count', hue = 'effect on binding', data = effect_on_binding_byProtein, order = MutProteinList, hue_order = ['gain', 'loss'], palette = 'muted', ax = axes)
        axes.set_xlabel('')
        axes.set_ylabel('Number of mutations')
        plt.grid()
        plt.tight_layout()
        plt.savefig('%s/%s_PROTEIN_SPECIFIC_Only_Conservative_%s'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))
        plt.clf() 



def Print_mutation_analysis_figure(RefResidue_frame, Ref_resid_pos_count, RefResidue_count, MutResidue_frame, Mut_resid_pos_count, MutResidue_count, FileName, HLA_type):
    
    fig, axes = plt.subplots(2, 3, sharex = 'col', gridspec_kw={'width_ratios':[1,1,0.33]}, figsize = (10, 12)) #, sharey = 'row',
    sns.set(font_scale=1.5)
    sns.heatmap(RefResidue_frame, cmap = "YlGnBu", yticklabels = 1, ax = axes[0,0])
    sns.heatmap(Ref_resid_pos_count, cmap = "YlGnBu", yticklabels = 1, ax = axes[0,1])
    sns.barplot(x = 'reference residue rate', y = 'reference residue', data = RefResidue_count, color = 'green', ax = axes[0,2])
    sns.heatmap(MutResidue_frame, cmap = "YlGnBu", yticklabels = 1, ax = axes[1,0])
    sns.heatmap(Mut_resid_pos_count, cmap = "YlGnBu", yticklabels = 1, ax = axes[1,1])

    sns.barplot(x = 'mutated residue rate', y = 'mutated residue', data = MutResidue_count, color = 'green', ax = axes[1,2])
    axes[0,1].set_xlabel('')
    axes[1,1].set_xlabel('')
    axes[0,1].set_ylabel('')
    axes[1,1].set_ylabel('')
    axes[0,0].set_xlabel('')
    axes[1,0].set_xlabel('')
    plt.setp( axes[0,1].yaxis.get_majorticklabels(), rotation=0 )
    plt.setp( axes[1,1].yaxis.get_majorticklabels(), rotation=0 )
    plt.setp( axes[1,0].xaxis.get_majorticklabels(), rotation=0 )
    plt.setp( axes[1,1].xaxis.get_majorticklabels(), rotation=0 )

    plt.grid()
    plt.tight_layout()
    plt.savefig('%s/%s_MUTATION_ANALYSIS_MUTandREFResidues_QUADPLOT_%s'% (FileName, HLA_type, FileName))
    plt.clf()


'''
8. Function that sorts mutations into reoccurence categories.
inputs: mmutation rates from MutRate dictionary.
Outputs: corresponding categories, returned as strings.
'''

def MutationRate_Categories(x):
    if x < 15:
        return '0-15%'
    elif x >= 15 and x < 30:
        return '15-30%'
    elif x >= 30 and x < 50:
        return '30-50%'
    elif x >= 50 and x < 70:
        return '50-70%'
    else:
        return '70% +'


'''
7. Function that ouputs all figures.

inputs: FullTable (large dataframe with multiindex where level 0 and 1 =  mutation ID, Mutated peptide sequence respectively), FileName, RefPeptides (dataframe with same multiindex as FullTable, where column = reference peptide sequences.)
Outputs: figures, saved in folder named by string saved in FileName variable, .csv file with all data generated.
'''

def Print_HLAspecific_Plots(FullTable, RefPeptides, Additional_Data): #, FileName


    Processed_HLAs = {} 

    Concatenated_Full_Frames = pd.DataFrame()   

    #use separate function to read .txt output from first program (contains list of unique mutations and their mutation rates). Function returns Dictionary with Mutation as keys, mutation rate as item.
    
    MutRate = Additional_Data[0]
    Perform_IEDB_SARSCOV_alignment = Additional_Data[1]
    bindingIEDB = Additional_Data[2]
    tcellIEDB = Additional_Data[3]
    FileName = Additional_Data[4]

    
    RefPeptides.reset_index(level = 0, inplace = True) #set mutation ID and mutated peptide sequences as columns, leaving numerical index.
    RefPeptides.reset_index(level = 0, inplace = True)
    positions_List = [0, 1, 2, 3, '(...)', -4, -3, -2, -1]
    

    #Now, Full table has only one column (fold difference) per HLA for all mutations, and has multi-index where level 0 = mutation ID and level 1 = mutated peptide sequences.
    Rkdiff_Count = 0

    #for column in list(FullTable): #Iterate through each HLA (fold change). Here, each column has fold change for a single HLA, for peptides corresponding to all mutations.
        
    ColNames = FullTable.columns.values
    HLA_FC = ''
    for name in ColNames:
        if 'Rk_FC' in name:
            HLA_FC = name

    FC = FullTable[HLA_FC].copy()
    temp_Frame = FullTable.drop([HLA_FC], axis = 1)

    temp_Frame.insert(loc = 4, column = HLA_FC, value = FC)
    Rkdiff_Count += 1
    print('temp_frame_pre index reset\n%s'% (temp_Frame))
    temp_Frame.reset_index(level = 0, inplace = True) #set Mutation ID as a column (mutation ID = protein_X350Y). Column is now named 'level_0'
    temp_Frame.reset_index(level = 0, inplace = True) #Set Mutated peptide sequences (now 'index at level 0') as a column, named 'Peptide_Mut'. the Temp_Frame dataframe and Refpeptide dataframe now have the same index, and can be merged on 'level_0' (mutation ID), and Peptide_Mut.
    temp_Frame = pd.merge(temp_Frame, RefPeptides, on = ['Peptide_Mut', 'level_0']) #merging temp_frame and Refpeptides to add reference peptides to temp_Frame dataframe. There is only one set of RefPeptides corresponding to one set of mutated peptides. HLAs that do not present a particlar peptide pair (mut/Ref) will display N/A instead of fold change/rank diff.
    temp_Frame.insert(loc = 0, column = 'PepLength', value = temp_Frame['Peptide_Mut'].apply(len)) #insert column named 'PepLength', where we use .apply() to calculate length of mutated peptide sequences.
    temp_Frame.insert(loc = 0, column = 'Mut_Position', value = temp_Frame.apply(lambda x: compare_seqs(x['Peptide_Mut'], x['Peptide_Ref']), axis = 1))   
    temp_Frame.insert(loc = 0, column = 'Mut_per_peptides', value = temp_Frame.apply(lambda x: mutations_per_peptides(x['Peptide_Mut'], x['Peptide_Ref']), axis = 1))
    temp_Frame.insert(loc = 0, column = 'Mut_Protein', value = temp_Frame['level_0'].apply(lambda x: x[0:(x.find('_'))])) #insert Mut_Protein' column, where we use .apply() function to extract protein name from Mutation ID
    temp_Frame['level_0'] = temp_Frame['level_0'].apply(lambda x: x[0:x.find('_', x.find('_')+1)])
    print(f'This is tempframe columns after level_0 addition: {temp_Frame.columns.values}') 
    temp_Frame.insert(loc = 0, column = 'mutationRate', value = temp_Frame['level_0'].apply(lambda x: MutRate.loc[x, 'Mutation_Rates'])) 
    temp_Frame.insert(loc = 0, column = 'Total_mutations', value = temp_Frame['level_0'].apply(lambda x: MutRate.loc[x, 'Total_mutations']))
    temp_Frame.insert(loc = 1, column = 'Reoccuring', value = temp_Frame['mutationRate'].apply(MutationRate_Categories))
    temp_Frame.insert(loc = 0, column = 'mutation_type', value = temp_Frame.apply(lambda x: compare_seqs_ReturnMut(x['Peptide_Mut'], x['Peptide_Ref']), axis = 1))
    temp_Frame.insert(loc = 1, column = 'reference residue', value = temp_Frame['mutation_type'].str[0])
    temp_Frame.insert(loc = 2, column = 'mutated residue', value = temp_Frame['mutation_type'].str[-1])
    temp_Frame.insert(loc = 12, column = 'Mut rank qual', value = temp_Frame['Rank_Mut'].apply(lambda x: 1 if float(x) < 2.0 else 0))
    temp_Frame.insert(loc = 14, column = 'Ref rank qual', value = temp_Frame['Rank_Ref'].apply(lambda x: 1 if float(x) < 2.0 else 0))  
    temp_Frame.insert(loc = 15, column = 'effect on binding_nonConservative', value = temp_Frame.apply(lambda x: Effect_on_Binding_NonConservative(x['Rank_Mut'], x['Rank_Ref']), axis = 1))
    temp_Frame.insert(loc = 15, column = 'effect on binding_SemiConservative', value = temp_Frame.apply(lambda x: Effect_on_Binding_SemiConservative(x['Rank_Mut'], x['Rank_Ref']), axis = 1))
    temp_Frame.insert(loc = 15, column = 'effect on binding', value = temp_Frame.apply(lambda x: Effect_on_Binding(x['Rank_Mut'], x['Rank_Ref']), axis = 1))

    ######################
    #Drop Singletons and all mutations with less than 4 reps
    #temp_Frame.drop(temp_Frame[temp_Frame['Total_Mutations'] < 4].index, inplace = True)

    ######################

    temp_Frame = temp_Frame.drop(temp_Frame[temp_Frame['Peptide_Mut'].str.contains('X')].index)


    AllHLA_APs_per_protein = temp_Frame.groupby(['Mut_Protein'])['level_0'].nunique().reset_index(name = 'count')
    
    print('This is Temp_Frame before figures!!!! \n\n\n\n\n%s\n\n\n\n\n\n'% (temp_Frame))
    print(f'This is HLA_FC column {temp_Frame.columns[-3]}')

    Protein_Analysis(temp_Frame, MutRate, FileName)

    dataName = temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')].replace(':', '')

    print(f'TIS IS DATANAME: {dataName} ')

    temp_Frame.to_csv('%s/%s_temp_Frame_%s.csv'% (FileName, dataName, FileName))

    
    All_mutations = len(MutRate.index)
    allHLAs_only_HLA_APs_below2 = temp_Frame['level_0'].nunique()
    
    if Concatenated_Full_Frames.empty:
        Concatenated_Full_Frames = temp_Frame
        print('\n\n\n\n\n\nOMG YESSSSSSSSSSSS\n\n\n\n\n%s\n\n\n\n\n'% (Concatenated_Full_Frames))

    else:
        Concatenated_Full_Frames = pd.concat([Concatenated_Full_Frames, temp_Frame])
        print('\n\n\n\n\n\nOMGGGGGGGGGGGG\n\n\n\n\n%s\n\n\n\n\n'% (Concatenated_Full_Frames))

    

    only_currentHLA_only_HLA_APs = temp_Frame.dropna(subset = [temp_Frame.columns[-3]])['level_0'].nunique()
    
    effect_on_binding_bymutation = temp_Frame.dropna(subset = [temp_Frame.columns[-3]]).groupby(['level_0', 'effect on binding'])['mutation_type'].count().reset_index(name = 'count')
    effect_on_binding_bymutation.to_csv('%s/%s_effect_on_binding_bymutation_%s.csv'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))
    
    gain_frame = effect_on_binding_bymutation[effect_on_binding_bymutation['effect on binding'] == 'gain']
    loss_frame = effect_on_binding_bymutation[effect_on_binding_bymutation['effect on binding'] == 'loss']
    no_effect_frame = effect_on_binding_bymutation[effect_on_binding_bymutation['effect on binding'] == 'no effect']

    effect_on_binding_bymutation_non_conservative = temp_Frame.dropna(subset = [temp_Frame.columns[-3]]).groupby(['level_0', 'effect on binding_nonConservative'])['mutation_type'].count().reset_index(name = 'count')

    gain_frame_nonCons = effect_on_binding_bymutation_non_conservative[effect_on_binding_bymutation_non_conservative['effect on binding_nonConservative'] == 'gain']
    loss_frame_nonCons = effect_on_binding_bymutation_non_conservative[effect_on_binding_bymutation_non_conservative['effect on binding_nonConservative'] == 'loss']
    no_effect_frame_nonCons = effect_on_binding_bymutation_non_conservative[effect_on_binding_bymutation_non_conservative['effect on binding_nonConservative'] == 'no effect']

    effect_on_binding_bymutation_Semi_conservative = temp_Frame.dropna(subset = [temp_Frame.columns[-3]]).groupby(['level_0', 'effect on binding_SemiConservative'])['mutation_type'].count().reset_index(name = 'count')

    gain_frame_SemiCons = effect_on_binding_bymutation_Semi_conservative[effect_on_binding_bymutation_Semi_conservative['effect on binding_SemiConservative'] == 'gain']
    loss_frame_SemiCons = effect_on_binding_bymutation_Semi_conservative[effect_on_binding_bymutation_Semi_conservative['effect on binding_SemiConservative'] == 'loss']
    no_effect_frame_SemiCons = effect_on_binding_bymutation_Semi_conservative[effect_on_binding_bymutation_Semi_conservative['effect on binding_SemiConservative'] == 'no effect']

    

    gain_only_list = []
    loss_only_list = []
    noeffect_only_list = []
    gain_and_loss_list = []
    for items in effect_on_binding_bymutation.level_0.unique():
        print(items)
        print(len(effect_on_binding_bymutation.level_0.unique()))
        if items in set(gain_frame['level_0']) and items in set(loss_frame['level_0']):
            gain_and_loss_list.append(items)
        elif items in set(gain_frame['level_0']) and items not in set(loss_frame['level_0']):
            gain_only_list.append(items)
        elif items not in set(gain_frame['level_0']) and items in set(loss_frame['level_0']):
            loss_only_list.append(items)
        elif items not in set(gain_frame['level_0']) and items not in set(loss_frame['level_0']) and items in set(no_effect_frame['level_0']):
            noeffect_only_list.append(items)

        



    
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.set(font_scale=1.5)
    venn2([set(gain_frame['level_0']), set(loss_frame['level_0'])], ('Gain', 'Loss'))
    
    plt.tight_layout()
    plt.savefig('%s/%s_VENNDIAGRAM_Conservative_%s'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))
    plt.clf()

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.set(font_scale=1.5)
    venn2([set(gain_frame_nonCons['level_0']), set(loss_frame_nonCons['level_0'])], ('Gain', 'Loss'))
   
    plt.tight_layout()
    plt.savefig('%s/%s_VENNDIAGRAM_nonConservative_%s'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))
    plt.clf()

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.set(font_scale=1.5)
    venn2([set(gain_frame_SemiCons['level_0']), set(loss_frame_SemiCons['level_0'])], ('Gain', 'Loss'))
   
    plt.tight_layout()
    plt.savefig('%s/%s_VENNDIAGRAM_SemiConservative_%s'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))
    plt.clf()

    data_dict = {'categories': ['Gain and Loss', 'Gain only *', 'Loss only **', 'No Effect only ***'],
                'values': [(len(gain_and_loss_list)/only_currentHLA_only_HLA_APs)*100, (len(gain_only_list)/only_currentHLA_only_HLA_APs)*100, (len(loss_only_list)/only_currentHLA_only_HLA_APs)*100, (len(noeffect_only_list)/only_currentHLA_only_HLA_APs)*100]
                }
    mutation_analysis = pd.DataFrame(data_dict, columns = ['categories', 'values'])

    with sns.axes_style("ticks"):
        fig, ax = plt.subplots(figsize=(8, 6)) 
        sns.set(font_scale=1)
        sns.despine()
        sns.barplot(x = 'categories', y = 'values', data = mutation_analysis)
        ax.set_xlabel('')
        ax.set_ylabel('rate of %s-associated polymorphisms'% ((temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')])))
        
        plt.xticks(rotation=45, horizontalalignment = 'right')
        plt.ylim(0, 100)
        plt.grid()
        plt.tight_layout()
        plt.savefig('%s/%s_GAIN_LOSS_barGraph_%s'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))
        plt.clf()


    Order = ['1', '2', '3', '4', '-4', '-3', '-2', '-1']
    effect_on_binding_byposition_nonConservative = temp_Frame.dropna(subset = [temp_Frame.columns[-3]]).groupby(['Mut_Position', 'effect on binding_nonConservative'])['level_0'].nunique().reset_index(name = 'count')
    effect_on_binding_byposition_nonConservative.rename(columns = {'effect on binding_nonConservative': 'effect on binding'}, inplace = True)
    
    effect_on_binding_byposition = temp_Frame.dropna(subset = [temp_Frame.columns[-3]])
    effect_on_binding_byposition.drop(effect_on_binding_byposition[effect_on_binding_byposition['effect on binding'] == 'no effect'].index, inplace = True)
    effect_on_binding_byposition = effect_on_binding_byposition.groupby(['Mut_Position', 'effect on binding'])['level_0'].nunique().reset_index(name = 'count')

    print(effect_on_binding_byposition)

    effect_on_binding_byposition.to_csv('%s/%s_effect_on_binding_MutPosition_%s.csv'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))

    
    TM_perc = temp_Frame.dropna(subset = [temp_Frame.columns[-3]]).groupby(['Mut_Position'])['mutation_type'].count().reset_index(name = 'count')
    Total_Mutations = len(temp_Frame.dropna().index)
    TM_perc.to_csv('%s/%s_TM_perc_%s.csv'% (FileName, temp_Frame.columns[-3][:temp_Frame.columns[-3].find('_')], FileName))
    TM_perc['count'] = (TM_perc['count'] / Total_Mutations) * 100

    
    TM = temp_Frame.groupby(['Mut_Position', 'PepLength'], as_index = False, sort = False)['Total_mutations'].sum()
    
    Make_Mut_Position_subplot(effect_on_binding_nonConservative = effect_on_binding_byposition_nonConservative, effect_on_binding = effect_on_binding_byposition, Order = Order, temp_Frame = temp_Frame, TM = TM_perc, FileName = FileName)
    

    temp_Frame.dropna(axis = 0, inplace = True)
    temp_Frame.reset_index(drop = True, inplace = True)

    
    
    if Perform_IEDB_SARSCOV_alignment.lower() == 'y':
        temp_Frame_IEDBFrame = PairwiseAlignment_Tcell_MHCbinding(bindingIEDB, tcellIEDB, temp_Frame, FileName)
        temp_Frame = temp_Frame_IEDBFrame[0]
        IEDBFrame = temp_Frame_IEDBFrame[1]
    
    
    temp_Frame.to_csv('%s/%s_ALL_DATA.csv'% (FileName, temp_Frame.columns[-2][:temp_Frame.columns[-2].find('_')]))
        
        



'''
6. This function will bring everything together (each mutation-specific dataframe containing fold change and rank difference for all HLAs) to create
one large dataframe containing fold change and rank difference for all mutations and all HLAs, in an organized fashion.

Here, we will also extract reference peptide sequences and store them in different dataframe (seems weird, but it makes the next step a bit easier...)

Inputs: Dataframe from function 5, FIleName.
Output: two dataframes: on with fold change for all mutations/all HLAs (with mutated and ref sequences as multiindex), and one with netMHCpan rank difference.

'''

def Create_FullTable(ConcFrames, FileName, Additional_Data):
    
    frames = []
    mutations = []
    Separating_lines = []
    row_count = 0

    for mut, frame in ConcFrames.items(): #iterate through multi-HLA dataframes of each mutation.
        frames.append(frame) #create a list of all the dataframes
        mutations.append(mut) #create a list of all the corresponding mutations
        row_count += len(frame.index)
        Separating_lines.append(row_count)
     
        
    FullTable = pd.concat(frames, keys = mutations) #concatenate all mutation-specific dataframes created in previous loop, using mutation list as keys. Now have 3-leveled multi index where level 0 = mutation ID, level 1 = mutated peptide sequence, level 2 = reference sequence.
    FullTable.reset_index(level = 2, inplace = True) #Set index level 2 (reference peptide sequences) as a column.
    
    
    RefPeptides = FullTable[['Peptide_Ref']].copy() #extract reference peptides, store in another dataframe.
    FullTable.drop('Peptide_Ref', axis = 1, inplace = True)


    Print_HLAspecific_Plots(FullTable, RefPeptides, Additional_Data) #Folder
    
    
    
'''
5. This function will further process data. Output is a new dictionary where each mutation has a large dataframe with fold change and difference in netMHCpan rank for all HLAs.
We will also shed information that is no longer used, such as netMHCpan mutated and reference rank, peptide position, and nM.
Mutated and reference sequences are both set as multiindex, since they are shared by all HLAs.

Inputs: dicionary from function 4, FileName.
Outputs: dataframe where each mutation (keys) has dataframe with fold change and netMHCpan rank difference for all HLAs
'''


def Create_Mutation_Dict(ProcessedMutDict, FileName, Additional_Data):
    ConcFrames = {}
    for mut, hla in ProcessedMutDict.items(): #Iterate through dictionary of mutations/HLAs/binding affinities to further manipulate data.
        first = 0
        TempFrame = pd.DataFrame() #create temp frame

        hla[1].drop(['Pos'], axis = 1, inplace = True) #Remove all unnecessary info. from now on, we only need sequences for mutated and non mutated peptides, fold change, and rank difference.
        TempFrame = hla[1].rename(columns={'Rk_FC':'%s_Rk_FC'% (hla[0])}).set_index(['Peptide_Mut', 'Peptide_Ref']) #Set sequences (mutated and reference) as multiindex, as these are shared by peptides in all HLAs. From now on, each HLA will have two columns: fold change and rank difference.
        print('TempFrame: \n%s\n'% (TempFrame))
        
        TempFrame.index.name = 'Mut_%s'% (mut)

        ConcFrames[mut] = TempFrame

    print(ConcFrames)
    Create_FullTable(ConcFrames, FileName, Additional_Data)

            
            #print(TempFrame)

    #ProcessedMutDict = ProcessedMutDict.drop(['Pos', 'nM_Mut', 'Rank', 'nM_Ref'])

'''
4. Function that takes in two dataframes, one with mutated peptides and one with their reference peptides (each pair is HLA and mutation specific).
Dataframes containg rank and nM binding affinities for mutated and reference peptides. 
Function will return rank difference as well as fold change for each peptide in dataframe and will shortlist peptide pairs in which either the mutated or reference peptide is binding..

inputs: two dataframes (one mutated and one reference) with relevent peptide info corresponding to a specific mutation and HLA.
Outputs: rank difference as well as fold change for each peptide relevant to a specific mutation/HLA combination.
'''


def Mut_vs_Ref(MutFrame, RefFrame):
    
    
    CombinedFrames = MutFrame.join(RefFrame['Peptide'], lsuffix = '_Mut', rsuffix = '_Ref') #join Peptide and rank of ref peptides to mutated peptides.
    CombinedFrames = CombinedFrames.join(RefFrame['Rank'], lsuffix = '_Mut', rsuffix = '_Ref')
    CombinedFrames = CombinedFrames.join(RefFrame['nM'], lsuffix = '_Mut', rsuffix = '_Ref')

    print('CombinedFrames: %s'% (CombinedFrames))

    Rank_Mut = CombinedFrames['Rank_Mut'].to_numpy()
    Rank_Mut = Rank_Mut.astype(np.float)
    Rank_Ref = CombinedFrames['Rank_Ref'].to_numpy()
    Rank_Ref = Rank_Ref.astype(np.float)

    FC = Rank_Mut - Rank_Ref
    Diff = Rank_Mut - Rank_Ref

    FC = FC * -1
    Diff = Diff * -1

    CombinedFrames['Rk_FC'] = FC
    CombinedFrames['Rk_diff'] = Diff

    location = np.where((Rank_Mut != Rank_Ref) & ((Rank_Mut <= 2.0) | (Rank_Ref <= 2.0))) #Shortlist only mutated/referece peptide pairs where either the mutated OR the reference peptide is predicted to bind.
    newFrame = CombinedFrames[CombinedFrames.index.isin(location[0])]

    return newFrame

'''
3. Function that splits each HLA dataframe into each unique mutation.
As a result, each HLA will have its own dictionary of mutations, where key = mutation ID and item = dataframe of corresponding mutated peptides and relevant info.

input: Dictionary with HLA-specific dataframes (from function 2)
outputs: Dictionary where key = mutations, value = list of lists, where each inner list = HLA type at pos 0, and dataframe of corresponding mutated peptides and relevant info at pos 1.
'''


def Split_Each_HLA(HLA, Data):
    MutationDict = {}
    
    #for HLA, Data in HLA_Data_Dict.items(): #iterate through HLAs in HLA dictionary from previous function.
        #set ID as new row names, then save them into new dataframe/remove them one by one.
    print('Data is %s'% (Data))
    MutationIDs = Data['ID'].to_numpy() #Turn mutationIDs into numpy array
    print('MutationIDs is %s'% (MutationIDs))
    indexes = np.unique(MutationIDs, return_index=True)[1] #get index of unique mutationIDs from array
    uniqueMutationIDs = [MutationIDs[index] for index in sorted(indexes)]

        
    print('uniqueMutationIDs is %s'% (uniqueMutationIDs))
    Column_names = Data.columns.values
    Data_array = Data.to_numpy() #Convert HLA-specific dataframes to numpy
    count = 1
        

    for IDs in uniqueMutationIDs: #Iterate through unique mutationIDs. Here, we consecutively acquire all data lines corresponding to mutated peptide, than all data corresponding to reference peptide.
        location = np.where(MutationIDs == IDs) #Acquire all rows of numpy
        if count == 1:
            MutName = IDs
            MutFrame = Data[Data.index.isin(location[0])]
            print('This is MutFrame: %s'% (MutFrame))
            print('This is count: %s'% (count))
            MutFrame = MutFrame.reset_index(drop = True)
            MutFrame.drop('ID', inplace = True, axis = 1)

            count += 1
        elif count == 2:
            RefFrame = Data[Data.index.isin(location[0])]
            print('This is RefFrame: %s'% (RefFrame))
            print('This is count: %s'% (count))
            RefFrame = RefFrame.reset_index(drop = True)
            RefFrame.drop('ID', inplace = True, axis = 1)

            CombinedFrames = Mut_vs_Ref(MutFrame, RefFrame)

            MutationDict[MutName] = [HLA, CombinedFrames]
            
            count = 1

    for keys, values in MutationDict.items():
        print('\n%s : %s\n'% (keys, values))
    return MutationDict

    

'''
2. Function to split data into all separate HLA types.
For each HLA type, create a new dataframe with relevant information (peptide position, peptide sequence, mutation ID (includes info about protein, etc.), predicted binding affinities (rank/nM))

Inputs: dataframe with all raw data, list of all HLAs in data (in right order)
Output: Dictionary where key = HLA type, value = dataframe with all important data for that HLA (peptide position, peptide sequence, mutation ID (includes info about protein, etc.), predicted binding affinities (rank/nM))
'''
def Split_into_HLAs(RawData, HLAs):
    HLA_Data_List = {}
    CommonInfo = RawData[['Pos','Peptide','ID']].copy() #gather info common to all HLA types into separate dataframe (peptide position, sequence, and mutation info)
    #print('Hello, %s'% (type(CommonInfo)))
    RawData.drop(['Pos','Peptide','ID', 'core', 'icore', 'EL-score', 'BA-score', RawData.columns[-2], RawData.columns[-1]], axis = 1, inplace = True) #Remove all info that will not be used, to simplyfy dataframe.
    ColNum = len(list(RawData.columns.values)) #after having removed unused columns, each HLA should have two columns: rank and nM. see how many HLAs were tested.
    if ColNum  == 2: #Only one HLA found.
        newTemp = pd.concat([CommonInfo, RawData], axis = 1) #Connect HLA info (rank/nM) with rest of info (pos, peptide sequence, ID)
        HLA_Data_List[HLAs[0]] = newTemp #Store in dataframe.
    else: #Else, more than on HLA
        RawData.columns = range(len(list(RawData.columns.values)))
        HLAnum = 0
        while len(list(RawData.columns.values)) > 0: #Iterate through HLAs with while loop
            TempFrame = RawData[[RawData.columns[0], RawData.columns[1]]].copy() #For each HLA, store its two columns (nM and rank) into temporary dataframe
            RawData.drop([RawData.columns[0], RawData.columns[1]], axis = 1, inplace = True) #Remove those columns from original dataframe.
            
            TempFrame.rename(columns = {'EL_Rank': 'Rank', 'BA_Rank' : 'nM'})

            TempFrame.columns = ['Rank', 'nM'] #remane columns in temp dataframe containing rank and nM for that particular HLA
            newTemp = pd.concat([CommonInfo, TempFrame], axis = 1) #Combine rank/nM dataframe of HLA with 'commoninfo' (pos, peptide and ID). As a result, each HLA will have its rank and nM binding affinities, as well as other important info (peptide pos, sequence, and ID)
            
            HLA_Data_List[HLAs[HLAnum]] = newTemp #.T #transpos HLA dataframe to facilitate iteration through mutations, and store in dictionary (key = HLA type, item = corresponding dataframe)
            HLAnum += 1 #Iterate through HLA type names in list created in previous function.
        
    return HLA_Data_List


def Call_Functions(HLA, Data, FileName, Additional_Data):
    
    ProcessedMutationsDict = Split_Each_HLA(HLA, Data) #SPlit each HLA

    Create_Mutation_Dict(ProcessedMutationsDict, FileName, Additional_Data)


'''
1. Function to read netMHCpan output (in .csv format), and create a list with all HLA types from file.
'''

def ReadFile(FileName, Additional_Data):
    Raw_data = pd.read_csv(FileName) #Read csv file
    preHLAs = list(Raw_data.columns.values) #create list of HLA types from file
    HLAs = []
    for item in preHLAs: #Acquire names of HLAs
        if item[:3].upper() == 'HLA':
            HLAs.append(item)
    Raw_data.columns = list(Raw_data.iloc[0])
    Raw_data.drop([0], inplace = True)
    Raw_data.reset_index(drop = True, inplace = True)
    

    HLAspecific_Data_Dict = Split_into_HLAs(Raw_data, HLAs) #Split raw HLA file into HLAs

    #####Iterate through HLAs here for multuprocessing.
    procs = []
    for HLA, Data in HLAspecific_Data_Dict.items():
        
        proc = Process(target = Call_Functions, args=(HLA, Data, FileName, Additional_Data))

        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()

        







if __name__ == '__main__':
    while True:
        
        print(f'This is sys.argv: {sys.argv}')

        
        FileName = sys.argv[1]

        MutRate = MutationRate() #Import mutation info file
        Genome_IDsonly = MutRate.drop(['Total_mutations', 'Mutation_Rates'], axis = 1)

        frames = []

        
        Perform_IEDB_SARSCOV_alignment = sys.argv[4]
        if Perform_IEDB_SARSCOV_alignment.lower() == 'y':
            bindingIEDB = Import_bindingIEDB()
            tcellIEDB = Import_tcellIEDB()

        
        global Folder
        Folder = FileName[:-4] #Create new folder where we will store all the figures. It is named by the name of the netMHCpan .csv ouput file.
        
        os.makedirs(os.path.dirname('./%s/'% (Folder)), exist_ok=True)
        

        Additional_Data = [MutRate, Perform_IEDB_SARSCOV_alignment, bindingIEDB, tcellIEDB, Folder]

        try:

            starttime = time.time()

            ReadFile(FileName, Additional_Data)

            print('Time taken = {} seconds'.format(time.time() - starttime))

        except FileNotFoundError:
            print('This file does not exist')
        else:
            break


