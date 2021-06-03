import pandas as pd
from pandas import DataFrame
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats
from scipy.stats import ttest_ind
import datetime
import math
from matplotlib_venn import venn2
from Bio import pairwise2
import sys
from multiprocessing import Process, Pipe

def find_core(x, Peptide_array, core_array):
    Peptide_index = np.where(Peptide_array == x)
    print('Peptide_index %s'% (Peptide_index))
    print('Peptide_index[0][0] %s'% (Peptide_index[0][0]))
    core_sequence = np.take(core_array, Peptide_index[0][0])
    print('core_sequence %s'% (core_sequence))
    if core_sequence != '#NAME?':
        print('core_sequence %s'% (core_sequence))
        #core = np.delete(core_sequence, np.where(core_sequence == '-'))
        #core = np.char.strip(core_sequence, chars = '-')
        print('core_sequence %s'% (core_sequence))
        return core_sequence

def Import_results_csv(file):
    gain_loss = pd.read_csv(file)
    return gain_loss.reset_index()

def get_value(array):
    if len(array) == 0:
        return np.nan
    else:
        return array[0]

def find_Mut_in_Core(ref_core, mut_core, ref_peptide, mut_peptide):
    
    if ref_core and mut_core:
        mut_position = np.where( np.array(list(ref_peptide)) != np.array(list(mut_peptide))  )
        if len(mut_position[0]) > 1:
            MutPos = mut_position[0][0]
        else:
            MutPos = mut_position[0]
        ref_core_dash_pos = []
        mut_core_dash_pos = []
        if '-' in ref_core:
            ref_core_dash_pos = np.where( np.array(list(ref_core)) == '-' )
            ref_core = ref_core[:ref_core.find('-')] + ref_core[ (ref_core.find('-') + 1) :]
        if '-' in mut_core:
            mut_core_dash_pos = np.where(np.array(list(mut_core)) == '-' )
            mut_core = mut_core[:mut_core.find('-')] + mut_core[ (mut_core.find('-') + 1) :]

        alignments_ref = pairwise2.align.globalxx(ref_peptide, ref_core)
        aligned_ref = np.array(list(alignments_ref[0][1]))

        alignments_mut = pairwise2.align.globalxx(mut_peptide, mut_core)
        aligned_mut = np.array(list(alignments_mut[0][1]))

        ref_dashes = np.where(aligned_ref == '-')
        print(ref_dashes)

        mut_dashes = np.where(aligned_mut == '-')
        print(mut_dashes)


        Mutation_in_core_REF = np.zeros((len(ref_peptide),), dtype=int)
        Mutation_in_core_REF[MutPos] = 1
        Mutation_in_core_REF = np.delete(Mutation_in_core_REF, ref_dashes[0])
        ref_mut_pos = np.where(Mutation_in_core_REF == 1)[0]

        Mutation_in_core_MUT = np.zeros((len(ref_peptide),), dtype=int)
        Mutation_in_core_MUT[MutPos] = 1
        Mutation_in_core_MUT = np.delete(Mutation_in_core_MUT, mut_dashes[0])
        mut_mut_pos = np.where(Mutation_in_core_MUT == 1)[0]

        if mut_mut_pos.size == 0 and ref_mut_pos.size == 0:
            return (int(np.array([100])), int(np.array([100])), False, 'NO', 'NO')
        if mut_mut_pos.size == 0: 
            return (int(ref_mut_pos), int(np.array([100])), False, 'NO', 'NO')

        if ref_mut_pos.size == 0:
            return (int(np.array([100])), int(mut_mut_pos), False, 'NO', 'NO')

        elif np.array_equal(ref_dashes, mut_dashes):
            if len(ref_core_dash_pos) > 0:
                if ref_core_dash_pos[0] < ref_mut_pos:
                    ref_mut_pos += 1
                    
            if len(mut_core_dash_pos) > 0:
                if mut_core_dash_pos[0] < mut_mut_pos:
                    mut_mut_pos += 1
                    
            equal = (ref_mut_pos == mut_mut_pos)
            return (int(ref_mut_pos)+1, int(mut_mut_pos)+1, equal, 'YES', 'YES')
        else:
            ref_core_mut_position = np.zeros((len(ref_peptide),), dtype=int)
            ref_core_mut_position[MutPos] = 1
            ref_core_mut_position = np.delete(ref_core_mut_position, ref_dashes[0])
            print(ref_core_mut_position)
            ref_mutation_in_core = np.nonzero(ref_core_mut_position)
            ref_mut_pos = ref_mutation_in_core[0]
            if len(ref_core_dash_pos) > 0:
                if ref_core_dash_pos[0] < ref_mutation_in_core[0]:
                    ref_mut_pos += 1

            mut_core_mut_position = np.zeros((len(ref_peptide),), dtype=int)
            mut_core_mut_position[MutPos] = 1
            mut_core_mut_position = np.delete(mut_core_mut_position, mut_dashes[0])
            print(mut_core_mut_position)
            mut_mutation_in_core = np.nonzero(mut_core_mut_position)
            mut_mut_pos = mut_mutation_in_core[0]
            if len(mut_core_dash_pos) > 0:
                if mut_core_dash_pos[0] < mut_mutation_in_core[0]:
                    mut_mut_pos += 1
                    

            if np.array_equal(ref_core_mut_position, mut_core_mut_position):
                print('%s_%s  YES'% (ref_mut_pos, mut_mut_pos))
                equal = (ref_mut_pos == mut_mut_pos)
                return (int(ref_mut_pos)+1, int(mut_mut_pos)+1, equal, 'YES', 'NO')
            else:
                print('%s_%s  YES'% (ref_mut_pos, mut_mut_pos))
                equal = (ref_mut_pos == mut_mut_pos)
                return (int(ref_mut_pos)+1, int(mut_mut_pos)+1, equal, 'YES', 'NO')


def Make_transition_heatmap(heatmap, HLA, folder, binding_effect, direction):
    with sns.axes_style("white"):
        fig, axes = plt.subplots(figsize = (10, 8)) 
        sns.set(font_scale=4)
        sns.heatmap(heatmap, cmap = "YlGnBu", yticklabels = 1, vmin = 0, vmax = 25) 
        
        axes.set_xlabel('Mutated Residue')

        if direction == 'allHLAs':
            axes.set_ylabel('Reference residue', fontsize = 40)
            axes.tick_params(axis='y', labelsize=40)
            axes.set_xlabel('Mutated Residue', fontsize = 40)
            axes.tick_params(axis='x', labelsize=40)

                
        else:
            plt.rcParams['font.size'] = 20
        plt.rcParams['font.size'] = 20
        plt.grid()
        plt.tight_layout()
        plt.yticks(rotation=0, horizontalalignment = 'right')
        plt.savefig('%s/%s_%s'% (folder, HLA, binding_effect))
        plt.clf()


'''
make_mutation_transitions_figure function:
This function generates heatmaps for the preferred types of substitutions (x-->z) associated with either loss or gain of epitopes.
This is done for each individual HLA types (direction = 'single') or all pooled HLA-related data (direction = allHLAs)
'''

def make_mutation_transitions_figure(temp_Frame, reform_HLA, Results_withCore, direction):
    
    
    ##################MutAnalysis_Fig = temp_Frame[temp_Frame['Refpos_sameAs_MutPos']== True] #only pick rows where reference residue and mutated residue are in same position of the binding groove.
    ##################MutAnalysis_Fig = MutAnalysis_Fig[MutAnalysis_Fig['effect on binding'] != 'no effect'] #only pick rows corresponding to to mutations leading to strong gain or loss (non-binder --> strong binder, vice versa)
    MutAnalysis_Fig = temp_Frame[temp_Frame['effect on binding'] != 'no effect'] #only pick rows corresponding to to mutations leading to strong gain or loss (non-binder --> strong binder, vice versa)

    if direction != 'single':
        reform_HLA = 'allHLAS'
    gain_fig = MutAnalysis_Fig[MutAnalysis_Fig['effect on binding'] == 'gain']
    loss_fig = MutAnalysis_Fig[MutAnalysis_Fig['effect on binding'] == 'loss']


    #save the following as pandas series (for individual HLAs if direction = 'single' or all pooled data if direction = 'allHLAs'): 

    mut_peptides_gain = gain_fig[['Peptide_Mut']].copy() #mutated peptides associated with gain
    ref_peptides_gain = gain_fig[['Peptide_Ref']].copy() #Reference peptides associated with gain
    mut_peptides_loss = loss_fig[['Peptide_Mut']].copy() #mutated peptides associated with loss
    ref_peptides_loss = loss_fig[['Peptide_Ref']].copy() #Reference peptides associated with loss

    print('This is loss fig IN_TRANSITION_FIGURE_FUNCT!!!!!%s'% (loss_fig))


    #Determine number of unique mutations correspinding to each type of substitution (x-->z) detected, for gain-associated mutations.
    Trantision_Count_gain = gain_fig.groupby('mutation_type')['level_0'].nunique().reset_index(name = 'count')
    Trantision_Count_gain.insert(loc = 1, column = 'reference residue', value = Trantision_Count_gain['mutation_type'].str[0])
    Trantision_Count_gain.insert(loc = 2, column = 'mutated residue', value = Trantision_Count_gain['mutation_type'].str[-1])
    Trantision_Count_gain.drop(['mutation_type'], axis = 1, inplace = True)

    #Pivot data to heatmap-compatible format, with reference residues along y axis, mutated residues along x axis, and number of unique mutations corresponding to each type of mutation as the values.
    print('This is Trantision_Count_gain %s\n\n\n\n\n'% (Trantision_Count_gain))
    heatmap_Table_gain = Trantision_Count_gain.pivot(index = 'reference residue', columns = 'mutated residue', values = 'count')
    print(('This is heatmap_Table %s\n\n\n')% (heatmap_Table_gain))
    #heatmap_Table_gain.to_excel('%s/%s_gain.xlsx'% ('Results_withCore', reform_HLA))
    heatmap_Table_gain.to_csv('%s/%s_gain.csv'% ('Results_withCore', reform_HLA))
    Make_transition_heatmap(heatmap_Table_gain, reform_HLA, 'Results_withCore', 'gain', direction) #Make heatmap


    #Determine number of unique mutations correspinding to each type of substitution (x-->z) detected, for gain-associated mutations.
    Trantision_Count_loss = loss_fig.groupby('mutation_type')['level_0'].nunique().reset_index(name = 'count')
    Trantision_Count_loss.insert(loc = 1, column = 'reference residue', value = Trantision_Count_loss['mutation_type'].str[0])
    Trantision_Count_loss.insert(loc = 2, column = 'mutated residue', value = Trantision_Count_loss['mutation_type'].str[-1])
    Trantision_Count_loss.drop(['mutation_type'], axis = 1, inplace = True)

    #Pivot data to heatmap-compatible format, with reference residues along y axis, mutated residues along x axis, and number of unique mutations corresponding to each type of mutation as the values.
    print('This is Trantision_Count_loss %s\n\n\n\n\n'% (Trantision_Count_loss))
    heatmap_Table_loss = Trantision_Count_loss.pivot(index = 'reference residue', columns = 'mutated residue', values = 'count')
    print(('This is heatmap_Table %s\n\n\n')% (heatmap_Table_loss))
    #heatmap_Table_loss.to_excel('%s/%s_loss.xlsx'% ('Results_withCore', reform_HLA))
    heatmap_Table_loss.to_csv('%s/%s_loss.csv'% ('Results_withCore', reform_HLA))
    Make_transition_heatmap(heatmap_Table_loss, reform_HLA, 'Results_withCore', 'loss', direction) #Make heatmap

    return [gain_fig, loss_fig, mut_peptides_gain, ref_peptides_gain, mut_peptides_loss, ref_peptides_loss]

    '''
    [0]: all HLA-specific data for mutations leading to gain (non-binder to strong binder)
	[1]: all HLA-specific data for mutations leading to loss (strong binder to non-binder to )
	[2]: mutated peptides for ref/mut peptide pairs associated with gain (for current HLA type)
	[3]: reference peptides for ref/mut peptide pairs associated with gain (for current HLA type)
	[4]: mutated peptides for ref/mut peptide pairs associated with loss (for current HLA type)
	[5]: reference peptides for ref/mut peptide pairs associated with gain (for current HLA type)
	'''


'''
make_mut_pos_figure function:
This function generates figure of association bewteen binding groove positions and gain or loss mutations.
This is done either for every single HLA individually (direction = 'single'), or for pooled data from all HLAs tested (direction = allHLAs)
'''

def make_mut_pos_figure(temp_Frame, HLA, folder, direction):
    
    MutAnalysis_Fig = temp_Frame[temp_Frame['Refpos_sameAs_MutPos']== True]  #only pick rows where reference residue and mutated residue are in same position of the binding groove.
    MutAnalysis_Fig = MutAnalysis_Fig[MutAnalysis_Fig['effect on binding'] != 'no effect'] #only pick rows corresponding to to mutations leading to strong gain or loss (non-binder --> strong binder, vice versa)
    if direction != 'single':
        HLA = 'All_HLAs'
    print('This is MutAnalysis_Fig: %s'% (MutAnalysis_Fig))
    MutAnalysis_Fig = MutAnalysis_Fig.groupby(['effect on binding', 'MutPos_in_Core'])['level_0'].nunique().reset_index(name = 'count') #_SemiConservative
    #MutAnalysis_Fig.rename(columns = {'effect on binding_SemiConservative': 'effect on binding'}, inplace = True)
    print(MutAnalysis_Fig)

    MutAnalysis_Fig.to_csv('%s/%s_Mut_Position.csv'% (folder, HLA))

    #Order = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    Order = [1, 2, 3, 4, 5, 6, 7, 8, 9] #Set x axis of bar plot.

    with sns.axes_style("ticks"):
        fig, axes = plt.subplots(figsize=(15, 10)) # , gridspec_kw={'width_ratios':[1,0.02]}, figsize=(14, 12)
        sns.set(font_scale=4)
        sns.despine()
    
       
        sns.barplot(x = 'MutPos_in_Core', y = 'count', hue = 'effect on binding', data = MutAnalysis_Fig, order = Order, hue_order = ['gain', 'loss'], palette = 'muted', ax = axes)
       
        axes.set_xlabel('Binding Groove Position')
        axes.set_ylabel('Number of unique mutations')
         
        if direction == 'allHLAs':
            axes.set_ylabel('Number of unique\n mutations', fontsize = 50)
            axes.tick_params(axis='y', labelsize=50)
            axes.set_xlabel('Binding Groove Position', fontsize = 50)
            axes.tick_params(axis='x', labelsize=50)

            
        else:
            axes.set_ylabel('Number of unique\n mutations', fontsize = 30)
            axes.tick_params(axis='y', labelsize=30)
            axes.set_xlabel('Binding Groove Position', fontsize = 30)
            axes.tick_params(axis='x', labelsize=30)
        plt.grid()
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.tight_layout()
        plt.savefig('%s/%s_Mut_Position'% (folder, HLA))
        plt.clf()


'''
Combine_Data function:
Print pooled (for all HLAs tested) mutated/reference peptides associated with gain or loss. This is to generate seq2logo from outputted files. 
This function also generates boxplot of average number of mutations associated with gain or loss for all HLA types tested (output is boxplot)
'''


def Combine_Data(Gain_quants, Loss_quants, mut_peptides_gain, ref_peptides_gain, mut_peptides_loss, ref_peptides_loss, folder):
    print('gain list: %s'% (Gain_quants))
    print('loss list: %s'% (Loss_quants))

    BoxPlot_Frame = pd.DataFrame({'gain': Gain_quants, 'loss':Loss_quants})
    BoxPlot_Frame.to_csv('%s/BoxPlot_Frame.csv'% (folder))
    BoxPlot_Frame = BoxPlot_Frame.stack().to_frame().reset_index()
    BoxPlot_Frame.rename(columns = {'level_1': 'impact', 0: 'num'}, inplace = True)

    print(BoxPlot_Frame)

    mut_peptides_gain = mut_peptides_gain.drop_duplicates()
    mut_peptides_gain.to_csv('%s/mut_peptides_gain.csv'% (folder))

    ref_peptides_gain = ref_peptides_gain.drop_duplicates()
    ref_peptides_gain.to_csv('%s/ref_peptides_gain.csv'% (folder))

    mut_peptides_loss = mut_peptides_loss.drop_duplicates()
    mut_peptides_loss.to_csv('%s/mut_peptides_loss.csv'% (folder))

    ref_peptides_loss = ref_peptides_loss.drop_duplicates()
    ref_peptides_loss.to_csv('%s/ref_peptides_loss.csv'% (folder))

    t, p = ttest_ind(Gain_quants, Loss_quants, equal_var=False)
    print('\n\n\n\n\nTHIS IS TTEST!!!!!!%s\n\n\n\n\n'% (p))

    with sns.axes_style("ticks"):
        fig, axes = plt.subplots(figsize=(12, 9)) 
        sns.set(font_scale=1)
        sns.despine()
        sns.boxplot(x = 'impact', y = 'num', data = BoxPlot_Frame, palette = 'muted', ax = axes)
        sns.swarmplot(x = 'impact', y = 'num', data = BoxPlot_Frame, color = '.2', size = 20, ax = axes)
       
        axes.set_xlabel('mutation impact on pMHC interaction')
        axes.set_ylabel('Number of mutations')

        axes.set_ylabel('mutation impact on pMHC interaction', fontsize = 30)
        axes.tick_params(axis='y', labelsize=30)
        axes.set_xlabel('Number of mutations', fontsize = 30)
        axes.set_title('P-Value: %s'% (p), fontsize = 30)
        plt.grid()
        plt.tight_layout()
        plt.savefig('%s/gain_loss_boxPlot'% (folder))
        plt.clf()

'''
Incorporte_Core function:
Incorporate from raw netMHCpan output csv file the 'binding cores' of each mut/ref peptide pairs in the HLA-specific processed netMHCpan files.

function input:
- raw netMHCpan output file, processed and split into individual HLAs. Each HLA-specific raw output file is stored in HLA_Specific_DataDict dictionary.
- gain/loss results file names, stored in Results_file dictionary.
'''

def Incorporte_Core(HLA, value, Results_file, HLA_Specific_DataDict, conn):
  
    #for HLA, value in Results_file.items(): #Iterate through gain/loss results files (in other wors, iterate through HLAs queried).
    print('\n\n\n\n\nTHIS IS HLA: %s\n\n\n\n\n'% (HLA))
    reform_HLA = HLA[:-2] + ':' + HLA[-2] + HLA[-1] #reformat HLA name.
    print(reform_HLA)
    temp_Frame = Import_results_csv(value[:-1]) #For each loss/gain file name, import corresponding csv file, which is saved in same directory.
    temp_Frame.drop(['index', 'Unnamed: 0'], axis = 1, inplace = True) #reformat columns (drop unnecessary columns)
    temp_Frame.drop([temp_Frame.columns[-2], temp_Frame.columns[-2]], axis = 1, inplace = True)

    Peptide_array = HLA_Specific_DataDict[reform_HLA]['Peptide'].values #extract (as numpy array) all peptide sequences from raw netMHCpan file corresponding to the current HLA.
    print('Peptide_array: %s'% (Peptide_array))
    core_array = HLA_Specific_DataDict[reform_HLA]['core'].values #extract (as numpy array) all core binding sequences from raw netMHCpan file corresponding to the current HLA.
    print('core_array: %s'% (core_array))

    temp_Frame['Ref_Core'] = temp_Frame['Peptide_Ref'].apply(lambda x: find_core(x, Peptide_array, core_array)) #Call find_core function to index reference (unmutated) peptide sequences of HLA-specific loss/gain peptide pairs (from loss/gain csv file) to peptide array from raw netMHCpan file in order to identify correspondiong binding core of unmutated peptide in context of corresponding HLA type. 
    temp_Frame['Mut_Core'] = temp_Frame['Peptide_Mut'].apply(lambda x: find_core(x, Peptide_array, core_array)) #Call find_core function to index mutated peptide sequences of HLA-specific loss/gain peptide pairs (from loss/gain csv file) to peptide array from raw netMHCpan file in order to identify correspondiong binding core of mutated peptide in context of corresponding HLA type.
    temp_Frame['core_Mut_position'] = temp_Frame.apply(lambda x: find_Mut_in_Core(x['Ref_Core'], x['Mut_Core'], x['Peptide_Ref'], x['Peptide_Mut']), axis = 1) #call find_Mut_in_Core function to identify index of mutation in the context of the binding groove (binding core).


    print(temp_Frame)
    test = pd.DataFrame(temp_Frame['core_Mut_position'].values.tolist(), index=temp_Frame.index) #Create separate dataframe from column of lists.
    print(test)


    temp_Frame[['RefPos_in_Core', 'MutPos_in_Core', 'Refpos_sameAs_MutPos', 'Mutation_in_Groove', 'peptide_flat_inGroove']] = pd.DataFrame(temp_Frame['core_Mut_position'].values.tolist(), index=temp_Frame.index) #Unpack dataframe of lists, insert in temp_Frame dataframe.

    temp_Frame.drop(['core_Mut_position'], inplace = True, axis = 1) #Drop column of lists.

    #Reformat temp_Frame dataframe to have more interpretable csv output.
    Frame_to_csv = temp_Frame.set_index('level_0')
    Frame_to_csv.reset_index(inplace = True)
    print(Frame_to_csv)
    mutations = Frame_to_csv['level_0'].to_numpy()
    Frame_to_csv.insert(loc = 0, column = 'Mutations', value = mutations)
    Frame_to_csv.drop(['level_0'], inplace = True, axis = 1)
    Frame_to_csv.drop(['Reoccuring', 'Mut rank qual', 'Ref rank qual'], inplace = True, axis = 1)
    Mutated_peptides = Frame_to_csv['Peptide_Mut'].to_numpy()
    Frame_to_csv.drop(['Peptide_Mut'], inplace = True, axis = 1)
    Frame_to_csv.insert(loc = 24, column = 'Peptide_Mut', value = Mutated_peptides)

    Frame_to_csv.to_csv('%s/%s_withCore.csv'% ('Results_withCore', HLA)) #Output to csv.
    print(Frame_to_csv)

    core_gain = temp_Frame[temp_Frame['effect on binding'] == 'gain'] #Create new dataframe with data for 'gain' mutations for current HLA (non-binder to strong binder).
    mut_core_gain = core_gain[['Mut_Core']].copy() #Extrapolate mutated cores for 'gain' peptide pairs, save in pandas series.
    ref_core_gain = core_gain[['Ref_Core']].copy() #Extrapolate reference. cores for 'gain' peptide pairs, save in pandas series.
    core_loss = temp_Frame[temp_Frame['effect on binding'] == 'loss'] #Create new dataframe with data for 'loss' mutations for current HLA (non-binder to strong binder).
    mut_core_loss = core_loss[['Mut_Core']].copy() #Extrapolate mutated cores for 'loss' peptide pairs, save in pandas series.
    ref_core_loss = core_loss[['Ref_Core']].copy() #Extrapolate reference cores for 'loss' peptide pairs, save in pandas series.


    
    temp_ALLHLAs = temp_Frame[temp_Frame['Refpos_sameAs_MutPos']== True] #Create new dataframe where reference residue and mutated residue are in the same position with regards to the binding motif (no loop is formed or removed as result of mutation).
    temp_ALLHLAs = temp_ALLHLAs[temp_ALLHLAs['effect on binding'] != 'no effect'] #remove all rows where mutation does not cause transition from strong binder or vice versa (no effect)
    temp_ALLHLAs['HLA'] = reform_HLA #Add new column with HLA name repeated for every row, for later indexing.

    make_mut_pos_figure(temp_Frame, reform_HLA, 'Results_withCore', 'single') #Call functino to generate figure for preferred position of loss or gain mutations for each individual HLA type ('single').
    gain_loss_fig = make_mutation_transitions_figure(temp_Frame, reform_HLA, 'Results_withCore', 'single') #Call function to make figure for substitution types (x --> z) associated with either loss or gain, for each single HLA type ('single')
    
    

    conn.send([temp_ALLHLAs, mut_core_gain, ref_core_gain, mut_core_loss, ref_core_loss, gain_loss_fig, reform_HLA])
    conn.close()



def Split_into_HLAs(RawData, HLAs):
    HLA_Data_List = {}
    CommonInfo = RawData[['Pos','Peptide', 'ID']].copy() #gather info common to all HLA types into separate dataframe (peptide position, sequence, and mutation info)
    RawData.drop(['Pos','Peptide','ID', 'icore', 'EL-score', 'BA-score', RawData.columns[-2], RawData.columns[-1]], axis = 1, inplace = True) #Remove all info that will not be used, to simplyfy dataframe. '1-log50k', '1-log50k'
    ColNum = len(list(RawData.columns.values)) #after having removed unused columns, each HLA should have two columns: rank and nM. see how many HLAs were tested.
    if ColNum  == 3: #Only one HLA found.
        newTemp = pd.concat([CommonInfo, RawData], axis = 1) #Connect HLA info (rank/nM) with rest of info (pos, peptide sequence, ID)
        HLA_Data_List[HLAs[0]] = newTemp #Store in dataframe.
    else: #Else, more than on HLA
        RawData.columns = range(len(list(RawData.columns.values)))
        HLAnum = 0
        while len(list(RawData.columns.values)) > 0: #Iterate through HLAs with while loop
            TempFrame = RawData[[RawData.columns[0], RawData.columns[1], RawData.columns[2]]].copy() #For each HLA, store its two columns (nM and rank) into temporary dataframe
            RawData.drop([RawData.columns[0], RawData.columns[1], RawData.columns[2]], axis = 1, inplace = True) #Remove those columns from original dataframe.
            
            TempFrame.rename(columns = {'EL_Rank': 'Rank', 'BA_Rank' : 'nM'})

            TempFrame.columns = ['core', 'nM', 'Rank'] #remane columns in temp dataframe containing rank and nM for that particular HLA
            ####
            TempFrame['HLA'] = HLAs[HLAnum]
            newTemp = pd.concat([CommonInfo, TempFrame], axis = 1) #Combine rank/nM dataframe of HLA with 'commoninfo' (pos, peptide and ID). As a result, each HLA will have its rank and nM binding affinities, as well as other important info (peptide pos, sequence, and ID)
            newTemp.set_index('ID', inplace = True) #Set ID as index
            HLA_Data_List[HLAs[HLAnum]] = newTemp  #.T #transpos HLA dataframe to facilitate iteration through mutations, and store in dictionary (key = HLA type, item = corresponding dataframe)
            HLAnum += 1 #Iterate through HLA type names in list created in previous function.
    return HLA_Data_List

'''
1. Function to read netMHCpan output (in .csv format), and create a list with all HLA types from file.
'''

def ReadFile(FileName, prompt):
    Raw_data = pd.read_csv(FileName) #Read csv file
    preHLAs = list(Raw_data.columns.values) #create list of HLA types from file
    HLAs = []
    for item in preHLAs:
        if item[:3].upper() == 'HLA':
            HLAs.append(item)
    Raw_data.columns = list(Raw_data.iloc[0])
    Raw_data.drop([0], inplace = True)
    Raw_data.reset_index(drop = True, inplace = True)
    
    HLAspecific_Data_Dict = Split_into_HLAs(Raw_data, HLAs)
    #return HLAspecific_Data_Dict
    Results_Files = {}
    with open(prompt, 'r') as File:   #open and save to Results_Files dictionary HLA-specific processed netMHCpan files.
        for line in File:
            if len(line) > 1:
                Results_Files[line[:line.find('_')]] = line
        Results_file = Results_Files
    
    loss_Tables = {}
    Gain_quants = []
    Loss_quants = []
    mut_peptides_gain = pd.DataFrame({'mut_peptides_gain': []})
    ref_peptides_gain = pd.DataFrame({'ref_peptides_gain': []})
    mut_peptides_loss = pd.DataFrame({'mut_peptides_loss': []})
    ref_peptides_loss = pd.DataFrame({'ref_peptides_loss': []})
    AllHLAs_DataFrame = pd.DataFrame()
    AllHLAs_loss_DataFrame = pd.DataFrame()
    AllHLAs_gain_DataFrame = pd.DataFrame()


    os.makedirs(os.path.dirname('./%s/'% ('Results_withCore')), exist_ok=True) #Open new folder named 'Results_withCore'

    #parent_conn, child_conn = Pipe()
    parent_num = 0
    Conns = {}
    procs = []
    for HLA, Data in Results_file.items():
        
        #parent_conn, child_conn = Pipe()
        Conns[parent_num] = Pipe()

        proc = Process(target = Incorporte_Core, args=(HLA, Data, Results_file, HLAspecific_Data_Dict, Conns[parent_num][1]))

        procs.append(proc)
        proc.start()

        parent_num += 1

    print(f'This is parent_num!!: {parent_num}')

    for num, Conn in Conns.items():

        items = Conn[0].recv()
        temp_ALLHLAs, mut_core_gain, ref_core_gain, mut_core_loss, ref_core_loss, gain_loss_fig, reform_HLA = items[0], items[1], items[2], items[3], items[4], items[5], items[6]


        
    #Create necessary variables before entering loop(s).
        

        AllHLAs_DataFrame = AllHLAs_DataFrame.append(temp_ALLHLAs) #As we loop through Results_file (iterate through gain/loss csv files for each HLA tested), we merge together data from all HLAs one by one to be able to perform one analysis on all HLA types at the same time (for example, get mutation positions in binding groove for all HLA types with P in motif).
        print('\n\n\n\n\nTHIS IS AllHLAs_DataFrame!!!!!%s\n\n\n\n'% (AllHLAs_DataFrame))

        loss_Tables[HLA] = gain_loss_fig[1] #Create dictionary for loss_assocated data for all HLA types.

        Gain_quants.append(gain_loss_fig[0]['level_0'].nunique()) #get all unique mutations associated with epitope gain for current HLA.
        
        Loss_quants.append(gain_loss_fig[1]['level_0'].nunique()) #get all unique mutations associated with epitope loss for current HLA.
        
        
        #As we loop through various gain/loss files (for each HLA), we append together unmutated peptides associated with gain, mutated peptides associated with gain, unmutated peptides associated with loss, and mutated peptides assocaited with loss in their respective HLA types, to be able to perform one single analysis on all HLAs at the same time.
        mut_peptides_gain = mut_peptides_gain.append(mut_core_gain) 
        ref_peptides_gain = ref_peptides_gain.append(ref_core_gain)
        
        mut_peptides_loss = mut_peptides_loss.append(mut_core_loss)
        ref_peptides_loss = ref_peptides_loss.append(ref_core_loss)



    for proc in procs:
        proc.join()


    combined_loss_table = pd.DataFrame()
    for HLA, table in loss_Tables.items(): #Combine 'loss data' from each HLA type into one.
        table[HLA] = HLA
        print('This is loss fig!!!!!%s'% (table))
        newTable = table[['level_0', table.columns[-1]]].copy()
        print('This is loss fig just level_0 and HLA!!!!!%s'% (newTable))
        newTable.drop_duplicates(inplace = True)
        print('This is loss fig just level_0 and HLA NO DUPLICATES!!!!!%s'% (newTable))
        if combined_loss_table.empty:
            combined_loss_table = newTable
        else:
            combined_loss_table = pd.merge(combined_loss_table, newTable, how = 'outer', on = ['level_0'])

    print(f'This is Gain_quants: {Gain_quants}')
    print(f'This is Loss_quants: {Loss_quants}')
    Combine_Data(Gain_quants, Loss_quants, mut_peptides_gain, ref_peptides_gain, mut_peptides_loss, ref_peptides_loss, 'Results_withCore') #Print pooled (for all HLAs tested) mutated/reference peptides associated with gain or loss. This is to generate seq2logo from outputted files. This function also generates boxplot of average number of mutations associated with gain or loss for all HLA types tested (output is boxplot)

    storage = make_mutation_transitions_figure(AllHLAs_DataFrame, reform_HLA, 'Results_withCore', 'allHLAs') #Geneate figure loss-, or gain-specific preferred substitution types (x --> z) for association with gain/loss mutations for ALL HLAs pooled together ('allHLAs').
    make_mut_pos_figure(AllHLAs_DataFrame, reform_HLA, 'Results_withCore', 'allHLAs') #Geneate figure for binding motif position-specific assssociation with gain/loss mutations for ALL HLAs pooled together ('allHLAs').


'''
def Results_file(): This function imports all necessary files:
 - netMHCpan raw output csv files for all HLAs queried here
 - processed netMHCpan files (mined for gain, loss) in csv format.
'''
def Results_file():
    while True:
        prompt = input('\n\nPlease, enter the file name for country 1, press Enter: ')
        try:
            Results_Files = {}
            with open(prompt, 'r') as File:
                for line in File:
                    Results_Files[line[:line.find('_')]] = line
                return [Results_Files, prompt]
        except FileNotFoundError:
            print('file not found. Make sure to copy in in same path as this script.')
        else:
            break

'''
Prompt user for input files:
- netMHCpan raw output csv files for all HLAs queried here
- gain/loss results file (txt) with file names for gain/loss results csv file.
'''

if __name__ == '__main__':
    while True:
        print(f'Inputs are as follows: {sys.argv}')
        FileName = sys.argv[1]
        prompt = sys.argv[2]
        try:
            ReadFile(FileName, prompt)
            
        except FileNotFoundError:
            print('This file does not exist')
        else:
            break

    