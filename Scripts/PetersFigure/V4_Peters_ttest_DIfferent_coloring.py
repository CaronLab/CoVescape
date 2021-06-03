import pandas as pd
from pandas import DataFrame
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats
from scipy.stats import ttest_ind
import math
import sys




def Import_results_csv(file):
    gain_loss = pd.read_csv(file)
    print(gain_loss)
    return gain_loss.reset_index()

def print_heatmaps(merged, direction, folder, datatype):
    
    merged = merged.sort_index()

    with sns.axes_style("ticks"):
        if direction == 'transition':
            fig, axes = plt.subplots(figsize = (10, 20))
        else:
            fig, axes = plt.subplots(figsize = (10, 8)) #, gridspec_kw={'width_ratios':[1,0.02]}, , figsize=(14, 12)
        sns.heatmap(merged, cmap = "YlGnBu", annot = True, yticklabels = 1) #, cbar_ax = axes[0,1]
        axes.set_xlabel('mutation position')
        axes.set_ylabel('%s_residues'% (direction))
 
        plt.grid()
        plt.tight_layout()
        plt.yticks(rotation=0, horizontalalignment = 'right')
        plt.savefig('%s/%s_%s_heatmap'% (folder, direction, datatype))
        plt.clf()

def letsMerge(Dict):
    merged = pd.DataFrame()
    for position, frame in Dict.items():
        print(frame)
        if merged.empty:
            merged = frame
        else:
            merged = pd.merge(merged, frame, how = 'outer', on = ['residue_list'])
    return merged

def Combine_positions_makeFigure(pValue_Dict, FC_Dict, direction, folder):
    print('This is with pvalue: %s\n%s'% (2, pValue_Dict[2.0]))
    print('This is with no pvalue: %s\n%s'% (2, FC_Dict[2.0]))
    merged_FC = pd.DataFrame()
    merged_pValue = pd.DataFrame()

    merged_FC = letsMerge(FC_Dict)
    merged_pValue = letsMerge(pValue_Dict)

    merged_FC.set_index('residue_list', inplace = True)
    merged_pValue.set_index('residue_list', inplace = True)
    
    merged_FC.to_csv('%s/%s_fold_changes.csv'% (folder, direction))
    merged_pValue.to_csv('%s/%s_p_values.csv'% (folder, direction))

    merged_FC = merged_FC.astype(float)
    merged_pValue = merged_pValue.astype(float)

    print(merged_FC)
    print(merged_FC.dtypes)
    print(merged_pValue)
    print(merged_pValue.dtypes)

    print_heatmaps(merged_FC, direction, folder, 'FC')
    print_heatmaps(merged_pValue, direction, folder, 'pvalue')
    

def make_boxplot_figure(boxplot_frame, Transition_order, folder, bargraph_frame): #
    with sns.axes_style("ticks"):
        
        fig, axes = plt.subplots(2, 1, sharex = 'col', gridspec_kw={'height_ratios': [0.35, 1]}, figsize=(30, 10)) 
        sns.set(font_scale=1)
        #sns.set_style('white')
        sns.despine()

        position_colors = {1.0: 'black', 2.0: 'red', 3.0: 'black', 4.0: 'black', 5.0: 'black', 6.0: 'black', 7.0: 'black', 8.0: 'black', 9.0: 'limegreen'}
        
        
        sns.barplot(x = 'transition', y = 'number of mutations', data = bargraph_frame, palette = ('blue',), ax = axes[0])
        
        sns.swarmplot(x = 'transition', y = 0, hue = 'position', palette = position_colors, data = boxplot_frame, edgecolor = 'black', linewidth = .8, size = 5, ax = axes[1]) 
        axes[0].set_ylabel('number of\n mutations', fontsize = 35)
        axes[0].tick_params(axis='y', labelsize=40)
        axes[0].set_xlabel('')
        axes[1].set_ylabel('-log2 (Fold Change)', fontsize = 35)
        axes[1].tick_params(axis='y', labelsize=40)
        axes[1].set_xlabel('Number of mutations', fontsize = 40)
        axes[1].tick_params(axis='x', labelsize=25)
        
        plt.legend([], [], frameon = False)
        plt.xticks(rotation=90, horizontalalignment = 'center')
        plt.grid()
        #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig('%s/boxplot'% (folder))
        plt.savefig('%s/boxplot.pdf'% (folder), dpi = 300)
        plt.clf()
        

def Get_P_values(Combined_HLAs, Residues, Mutation_positions, direction, folder):
    dict_pvalue = {}
    dict_FC = {}
    boxplot_frame = pd.DataFrame()

    for position in Mutation_positions:
        if int(position) == 100:
            continue
        else:
            Combined_HLAs_Background = Combined_HLAs[Combined_HLAs['MutPos_in_Core'] == position]

            Mut_Rank_Background = Combined_HLAs_Background['Rank_Mut'].to_numpy()
            Ref_Rank_Background = Combined_HLAs_Background['Rank_Ref'].to_numpy()
            Background_FC = Mut_Rank_Background/Ref_Rank_Background
            print(Background_FC)
            Background_FC_log2 = np.log2(Background_FC)
            print(Background_FC_log2)

            residue_list = []
            residue_list_wNumbers = []
            FC_list = []
            FC_list_average = []
            p_value_list = []
            print('\n\n\n\nThis is RESIDUES: %s'% (Residues))
            for Residue in Residues:
                if direction == 'transition':
                    Combined_HLAs_subPopulation = Combined_HLAs[Combined_HLAs['mutation_type'] == Residue]
                else:
                    Combined_HLAs_subPopulation = Combined_HLAs[Combined_HLAs['%s residue'% (direction)] == Residue]
                Combined_HLAs_subPopulation = Combined_HLAs_subPopulation[Combined_HLAs_subPopulation['MutPos_in_Core'] == position]

                if Combined_HLAs_subPopulation.empty or len(Combined_HLAs_subPopulation.index) == 1:
                    continue
                else:
                    print('This is residue \n\n%s\n\n\nThis is subpopulation\n\n%s\n\n'% (Residue, Combined_HLAs_subPopulation))


                    Mut_Rank_SubPop = Combined_HLAs_subPopulation['Rank_Mut'].to_numpy()
                    Ref_Rank_SubPop = Combined_HLAs_subPopulation['Rank_Ref'].to_numpy()

                    Subpop_FC = Mut_Rank_SubPop/Ref_Rank_SubPop
                    Subpop_FC_log2 = np.log2(Subpop_FC)
                    #print(Subpop_FC)
                    #print(Subpop_FC_log2)
                    t, p = stats.ttest_ind(Subpop_FC_log2,Background_FC_log2, equal_var = False)

                    print('\n\n\n\n\n\nTHIS IS subpop_FC, subpop_log2, background_FC, background_log2, and p for %s, %s: %s\n%s\n%s\n%s\n%s\n\n\n\n'% (Residue, position, Subpop_FC, Subpop_FC_log2, Background_FC, Background_FC_log2, p))

                    if direction == 'transition' and Residue == 'L-->P' and position == 2.0:
                        subpop_fc_PtoL = pd.DataFrame(Subpop_FC_log2)
                        subpop_fc_PtoL.to_csv('Subpop_FC_log2_L-->P.csv')
                        
                        background = pd.DataFrame(Background_FC_log2)
                        background.to_csv('Background_FC_log2.csv')

                    residue_list.append(Residue)
                    residue_list_wNumbers.append('%s     %s'% (str(int(position)), Residue))
                    FC_list.append(Subpop_FC_log2*-1)
                    FC_list_average.append((np.average(Subpop_FC_log2))*-1)
                    p_value_list.append(p)

                print('This is residue list, FC list and pvalue: \n\n%s\n\n%s\n\n%s'% (residue_list, FC_list_average, p_value_list))
            array_pvalue = np.array([p_value_list, residue_list])
            frame_pValue = pd.DataFrame(data = array_pvalue, index = [position, 'residue_list'])
            dict_pvalue[position] = frame_pValue.T

            array_FC = np.array([FC_list_average, residue_list])
            frame_FC = pd.DataFrame(data = array_FC, index = [position, 'residue_list'])
            dict_FC[position] = frame_FC.T


            if direction == 'transition':
                array_boxplot = pd.DataFrame(data = FC_list)
                array_boxplot['average log2'] = FC_list_average
                array_boxplot['p-values'] = p_value_list
                array_boxplot['transition'] = residue_list_wNumbers
                array_boxplot['position'] = position
                array_boxplot.set_index(['transition', 'position'], inplace = True)
                print('\n\nThis is array_box plot for %s BEFORE selection: %s\n\n'% (position, array_boxplot))
                array_boxplot = array_boxplot[(array_boxplot['average log2'] < -2.5) | (array_boxplot['average log2'] > 2.5)]
                array_boxplot = array_boxplot[(array_boxplot['p-values'] < 0.001)]
                print('\n\nThis is array_box plot for %s BEFORE stack: %s\n\n'% (position, array_boxplot))
                array_boxplot.drop(['p-values'], inplace = True, axis = 1)

                if boxplot_frame.empty:
                    boxplot_frame = array_boxplot
                else:
                    boxplot_frame = boxplot_frame.append(array_boxplot)
    if direction == 'transition':
        print('\n\n\nTHIS IS boxplot_frame: %s\n\n\n'% (boxplot_frame))
        boxplot_frame.sort_values(by = ['average log2'], inplace = True)
        boxplot_frame.reset_index(inplace = True)
        Transition_order = boxplot_frame['transition'].values
        boxplot_frame.set_index(['transition', 'position'], inplace = True)
        boxplot_frame['number of mutations'] = boxplot_frame.count(axis = 1)
        bargraph_frame = boxplot_frame.reset_index()
        bargraph_frame = bargraph_frame[['transition', 'number of mutations']].copy()
        mutnum = bargraph_frame['number of mutations'].to_numpy()
        
        bargraph_frame['number of mutations'] = mutnum

        boxplot_frame.drop(['average log2', 'number of mutations'], inplace = True, axis = 1) #, 'number of mutations'

        boxplot_frame = boxplot_frame.stack().to_frame().reset_index()
        boxplot_frame.to_csv('boxplot_frame.csv')

        if boxplot_frame.empty:
            print('There are no transitions that pass the fold change OR p-value cutoff!!')
        else:
            make_boxplot_figure(boxplot_frame, Transition_order, folder, bargraph_frame) #, bargraph_frame
    

    Combine_positions_makeFigure(dict_pvalue, dict_FC, direction, folder)




def get_core_positions(Combined_HLAs):
    refpos_inCore = Combined_HLAs['RefPos_in_Core'].to_numpy() 
    mutpos_inCore = Combined_HLAs['MutPos_in_Core'].to_numpy()
    shared_pos = np.where(refpos_inCore == mutpos_inCore)
    Core_positions = np.take(mutpos_inCore, shared_pos)
    return np.unique(Core_positions)

#Combine HLA-specific results dataframes, get log2(FD), perform t-test on verctor 1 (P-->X at position 2) and vector 2 (background; everything else at position 2)

def reference_Residues_TTest(Combined_HLAs, folder):
    Ref_Residues = Combined_HLAs['reference residue'].unique()

    Mut_Residues = Combined_HLAs['mutated residue'].unique()

    transition = Combined_HLAs['mutation_type'].unique()
    #mut_Residues = Combined_HLAs['mutated residue'].unique()
    Mutation_positions = get_core_positions(Combined_HLAs)

    Get_P_values(Combined_HLAs, Ref_Residues, Mutation_positions, 'reference', folder)

    Get_P_values(Combined_HLAs, Mut_Residues, Mutation_positions, 'mutated', folder)

    Get_P_values(Combined_HLAs, transition, Mutation_positions, 'transition', folder)

def Effect_on_Binding(Mut, Ref):
    if (float(Mut) <= 0.5) and (float(Ref) >= 2.0):
        return 'gain'
    elif (float(Mut) >= 2.0) and (float(Ref) <= 0.5):
        return 'loss'
    else:
        return 'no effect'

def Make_Mut_Transitions_Heamap(Combined_HLAs, folder, effect):
    #select only strong loss
    reduced_frame = Combined_HLAs[Combined_HLAs['effect on binding'] == effect]

    peptide_seqs = reduced_frame[['Peptide_Mut', 'Peptide_Ref']].copy()
    peptide_seqs['length'] = peptide_seqs['Peptide_Mut'].apply(lambda x: len(x))
    peptide_seqs = peptide_seqs[peptide_seqs['length'] == 9]
    
    Peptide_Mut = peptide_seqs[['Peptide_Mut']].copy()
    Peptide_Ref = peptide_seqs[['Peptide_Ref']].copy()

    Peptide_Mut.to_csv('%s/Peptide_Mut_%s.csv'% (folder, effect))
    Peptide_Ref.to_csv('%s/Peptide_Ref_%s.csv'% (folder, effect))

    #Groupby mutation type
    reduced_frame = reduced_frame.groupby('mutation_type')['Mutations'].nunique().reset_index(name = 'count')
    #Pivot to have reference as y, mutated as x
    reduced_frame.insert(loc = 1, column = 'reference residue', value = reduced_frame['mutation_type'].str[0])
    reduced_frame.insert(loc = 2, column = 'mutated residue', value = reduced_frame['mutation_type'].str[-1])
    reduced_frame.drop(['mutation_type'], axis = 1, inplace = True)
    reduced_frame_heatmap = reduced_frame.pivot(index = 'reference residue', columns = 'mutated residue', values = 'count')
    #make heatmap
    print('This is reduced frame heatmap: %s'% (reduced_frame_heatmap))
    reduced_frame_heatmap.to_csv('%s/%s_heatmap.csv'% (folder, effect))
    if not reduced_frame_heatmap.empty:
        with sns.axes_style("white"):
            fig, axes = plt.subplots(figsize = (15, 13)) 
            sns.set(font_scale=4)
            sns.heatmap(reduced_frame_heatmap, cmap = "YlGnBu", yticklabels = 1, xticklabels = 1)  #   , vmin = 0, vmax = 25   , cbar_ax = axes[0,1]   , annot = True
            axes.set_xlabel('Mutated Residue')

            
            axes.set_ylabel('Reference residue', fontsize = 30)
            axes.tick_params(axis='y', labelsize=30)
            axes.set_xlabel('Mutated Residue', fontsize = 30)
            axes.tick_params(axis='x', labelsize=30)

            plt.rcParams['font.size'] = 20
            plt.grid()
            plt.tight_layout()
            plt.yticks(rotation=0, horizontalalignment = 'right')
            plt.savefig('%s/%s_heatmap.pdf'% (folder, effect), dpi = 300)
            plt.clf()


def T_test(Results_file):

    folder = 'results_ALL_Stringent_Only'
    os.makedirs(os.path.dirname('./%s/'% (folder)), exist_ok=True)

    Combined_HLAs = pd.DataFrame()

    for HLA, value in Results_file.items(): #Iterate through gain/loss results files (in other wors, iterate through HLAs queried).
        reform_HLA = HLA[:-2] + ':' + HLA[-2] + HLA[-1] #reformat HLA name.
        print(reform_HLA)
        temp_Frame = Import_results_csv(value[:-1]) #For each loss/gain file name, import corresponding csv file, which is saved in same directory.
        #print('This is temp_Frame: %s'% (Combined_HLAs))
        temp_Frame.drop(['index', 'Unnamed: 0'], axis = 1, inplace = True) #reformat columns (drop unnecessary columns)
        reduced_Frame = temp_Frame[['Mutations', 'mutation_type', 'reference residue', 'mutated residue', 'Total_mutations', 'Rank_Mut', 'Rank_Ref', 'Peptide_Mut', 'Peptide_Ref', 'RefPos_in_Core', 'MutPos_in_Core', 'Refpos_sameAs_MutPos']].copy()
        Combined_HLAs = Combined_HLAs.append(reduced_Frame)
        
        print('This is Combined_HLAs: %s'% (Combined_HLAs))
    Combined_HLAs['effect on binding'] = Combined_HLAs.apply(lambda x: Effect_on_Binding(x['Rank_Mut'], x['Rank_Ref']), axis = 1)
    print('This is Combined_HLAs with dupli: %s'% (Combined_HLAs))
    Combined_HLAs.drop_duplicates(inplace= True)

    Combined_HLAs = Combined_HLAs[Combined_HLAs['effect on binding'] != 'no effect']

    print('This is Combined_HLAs  dupli: %s'% (Combined_HLAs))
    print(Combined_HLAs.dtypes)
    Make_Mut_Transitions_Heamap(Combined_HLAs, folder, 'gain')
    Make_Mut_Transitions_Heamap(Combined_HLAs, folder, 'loss')

    Combined_HLAs.to_csv('%s/Combined_HLAs.csv'% (folder))

    reference_Residues_TTest(Combined_HLAs, folder)



while True:
    prompt = sys.argv[1]
    try:

        Results_Files = {}
        with open(prompt, 'r') as File:   #open and save to Results_Files dictionary HLA-specific processed netMHCpan files.
            for line in File:
                if len(line) > 1:
                    Results_Files[line[:line.find('_')]] = line
            Results_file = Results_Files

    except FileNotFoundError:
        print('This file does not exist')
    else:
        break

T_test(Results_file)  #Call main method to initiate processing.








