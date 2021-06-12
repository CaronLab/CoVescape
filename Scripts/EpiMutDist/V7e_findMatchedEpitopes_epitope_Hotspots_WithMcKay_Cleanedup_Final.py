import pandas as pd
from pandas import DataFrame
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats
import datetime
from Bio import Align
from Bio.Align import substitution_matrices
import math
from matplotlib_venn import venn2
from Bio import SeqIO
import sys

            
            #print(TempFrame)

    #ProcessedMutDict = ProcessedMutDict.drop(['Pos', 'nM_Mut', 'Rank', 'nM_Ref'])
def Get_Cterm(pos, peptide):
    return int(pos) + len(peptide[:peptide.find('_')])

def epitope_density(pos, peptide, num):
    #print('%s, %s, %s'% (type(int(pos)), type(len(peptide)), type(num)))
    if int(pos) <= num <= (int(pos) + len(peptide[:peptide.find('_')])):
        return 1

def Epitope_Association(pos, EpitopeFrame):
    EpitopeFrame['HasMut'] = EpitopeFrame.apply(lambda x: epitope_density(x['Pos'], x['Peptide'], int(pos)), axis = 1)
    Count = EpitopeFrame['HasMut'].sum(axis = 0)
    EpitopeFrame.drop(['HasMut'], inplace = True, axis = 1)
    if Count > 0:
        return 1

def Mutation_density(pos, Total_mutations, num):
    if int(pos) == num:
        return Total_mutations
    #else:
        #return 0


def Mut_HLA_Association(HLA, Mutation):
    if HLA > 0 and Mutation > 0:
        return 1
    else:
        return 0

def Mut_Getquant(pos, x):
    if pos == x:
        return 1

def Make_Figure(AllProteinsData, Protein_Name, lengths, width, s1, s2, fig_size, linewidth, fig_num):

    with sns.axes_style("white"):
        fig, axes = plt.subplots(3, len(width), sharex = 'col', sharey = 'row', gridspec_kw={'height_ratios': [0.5, 1, 0.7], 'width_ratios':width}, figsize=fig_size) # , gridspec_kw={'width_ratios':[1,0.02]}, figsize=(14, 12) 'height_ratios':[0.35, 1, 1], 
        sns.set(font_scale=1.5) #, rc = {lines.linewidth = .5}
        #sns.set_style('white')
        #sns.set_context('notebook', font_scale = 10)
        sns.despine(top = True, right = True, left = True, bottom = True)
        POS = 0
        for length in lengths:
            #sns.lineplot(x = 'position', y = 'Rank', hue = 'Peptide', data = AllProteinsData[length][0], legend = False, color = 'blue', alpha = 0.5, ax = axes[0, POS])
            if length == lengths[-1]:
                #sns.lineplot(x = 'position', y = 0, data = AllProteinsData[length][1], color = 'blue', ax = axes[0, POS], legend = False)
                sns.lineplot(x = 'position', y = 'HLADensity', hue = 'HLA', data = AllProteinsData[length][2], linewidth = linewidth, ax = axes[1, POS], legend = False) #size = 'HLA', sizes = {'HLA-B40:01':1, 'HLA-A02:01':1, 'HLA-B44:02':1, 'HLA-B44:03':1, 'HLA-B08:01':1, 'HLA-B07:02':1, 'HLA-B35:01':1, 'HLA-A24:02':1, 'HLA-A23:01':1, 'HLA-A11:01':1, 'HLA-A03:01':1, 'HLA-A01:01':1, 'sum': 10}, 
                sns.set_style('white')
                sns.scatterplot(x = 'position', y = 'Mutation Density', data = AllProteinsData[length][3], s = s1, ax = axes[2, POS]) #alpha = .35,
                sns.scatterplot(x = 'position', y = 'Mutation Density', data = AllProteinsData[length][5], color = 'black', s = s2, clip_on = False, ax = axes[0, POS])
                sns.scatterplot(x = 'position', y = 'Mutation Density', data = AllProteinsData[length][6], color = 'red', s = s2, clip_on = False, ax = axes[0, POS])
            else:
                #sns.lineplot(x = 'position', y = 0, data = AllProteinsData[length][1], color = 'blue', ax = axes[0, POS], legend = False)
                sns.lineplot(x = 'position', y = 'HLADensity', hue = 'HLA', data = AllProteinsData[length][2], linewidth = linewidth, ax = axes[1, POS], legend = False)# size = 'HLA', sizes = {'HLA-B40:01':1, 'HLA-A02:01':1, 'HLA-B44:02':1, 'HLA-B44:03':1, 'HLA-B08:01':1, 'HLA-B07:02':1, 'HLA-B35:01':1, 'HLA-A24:02':1, 'HLA-A23:01':1, 'HLA-A11:01':1, 'HLA-A03:01':1, 'HLA-A01:01':1, 'sum': 10}
                sns.set_style('white')
                sns.scatterplot(x = 'position', y = 'Mutation Density', data = AllProteinsData[length][3], s = s1, ax = axes[2, POS]) #alpha = .35,
                sns.scatterplot(x = 'position', y = 'Mutation Density', data = AllProteinsData[length][5], color = 'black', s = s2, clip_on = False, ax = axes[0, POS])
                sns.scatterplot(x = 'position', y = 'Mutation Density', data = AllProteinsData[length][6], color = 'red', s = s2, clip_on = False, ax = axes[0, POS])
            #sns.barplot(x = 'position', y = 0, data = AllProteinsData[length][3], color = 'blue', ax = axes[2, POS])

            #axes[0, POS].set_ylim(-.3, 1)
            #axes[1, POS].set_ylim(-20000, 60000)
            axes[0, POS].set_title('%s'% (AllProteinsData[length][4]))
            axes[1, POS].set_ylabel('Epitope Density', fontsize = 20)
            axes[1, POS].tick_params(axis='y', labelsize=20)
            axes[2, POS].set_ylabel('Mutation frequency', fontsize = 20)
            axes[2, POS].tick_params(axis='y', labelsize=20)
            axes[2, POS].set_xlabel('residue position', fontsize = 20)
            axes[2, POS].tick_params(axis='x', labelsize=20)
            #axes[0, POS].yaxis.get_label().set_fontsize(20)
            '''
            for ind, label in enumerate(axes[1, POS].get_xticklabels()):
                if ind % 100 == 0:  # every 100th label is kept
                    label.set_visible(True)
                else:
                    label.set_visible(False)
            '''
            print('HEEEEEEEEEEELLLLLLLLLLLLLLLLLLOOOOOOOOOOOOOOOOOO!!!!!')
            #plt.setp(axes[0, POS].lines,linewidth=2, color = 'blue') #fig.
            plt.setp(axes[1, POS].lines,linewidth=5) #fig.
            plt.setp( axes[2, POS].xaxis.get_majorticklabels(), rotation=45 )

            print('BBBBBBBBBBBBOOOOOOOOOONNNNNNNNNNJJJJJJJJJJOOOOOOOOUUUUUUUURRRRR!!!!!')

            POS += 1
        

        plt.grid(b = False)
        #plt.rcParams.update({'font.size': 50})
        #axes[0, 0].set_xticklabels(axes[0, 0].get_xmajorticklabels(), fontsize = 200)
        #plt.legend(bbox_to_anchor=(1.05, 1), borderaxespad=0) #, loc=2
        plt.tight_layout()
        plt.savefig('%s_Epitope_HotPsots_line_bar%s.pdf'% (Protein_Name, fig_num), dpi = 300)
        plt.clf()


def Create_Mutation_Dict(ProcessedHLAFrame, Split_HLAs, RefProteins, MutData, Validated_Data, Validated_Data_Matched):
    
    

    MutData.insert(loc = 0, column = 'protein', value = MutData['MutID'].apply(lambda x: x[: x.find('_')]))
    MutData.insert(loc = 0, column = 'position', value = MutData['MutID'].apply(lambda x: x[(x.find('_') + 2): -1]))
    
    #Validated_Data_AllInfo = pd.merge(Validated_Data, MutData, how = 'left', on = ['MutID'])
    


    
    Mut_Epitope_Association_Rate = {}
    AllProteinsData = {}
    AllPro_EpitopeSum = pd.DataFrame()
    ProcessedHLAFrame.reset_index(inplace = True)
    #Proteins = ProcessedHLAFrame['ID'].unique()
    Proteins = list(RefProteins.keys())
    

    HLA_MutAssociation = {}
    folder = 'Protein_Specific_MutVsEpitopes'
    os.makedirs(os.path.dirname('./%s/'% (folder)), exist_ok=True)

    for Protein in Proteins:
        
        ProteinSpecific_HLAFrame = ProcessedHLAFrame[ProcessedHLAFrame['ID'] == Protein].copy()
        
        Protein_Name = Protein
        Protein_length = len(RefProteins[Protein])

        
        
        MutData_ProteinSpecific = MutData[MutData['protein'] == Protein].copy()
        MutData_ProteinSpecific = MutData_ProteinSpecific[MutData_ProteinSpecific['Total_mutations'] >= 4].copy()

        Validated_Data_AllInfo_Prot_Spec = Validated_Data[Validated_Data['Protein'] == Protein].copy()
        Validated_Data_Matched_Prot_Spec = Validated_Data_Matched[Validated_Data_Matched['Mut_Protein'] == Protein].copy()
        

        ###############################Validated_Data_AllInfo_Prot_Spec = MutData_ProteinSpecific[MutData_ProteinSpecific['Total_mutations'] >= 5].copy()


        MutData_ProteinSpecific['position'] = MutData_ProteinSpecific['position'].astype(float)
        #MutData_Spike['Total_mutations'] = MutData_Spike['Total_mutations'].apply(lambda x: np.log2(float(x)))
        


        Individual_HLA_CoverageMap = {}

        
        Epitope_Density = []
        Mutation_Density = []
        Max_TCellMagn_Val = []
        Max_TCellMagn_Matched_Val = []

        HLA_EpitopeDensity = pd.DataFrame()
        for num in range(1, Protein_length+1):
            ProteinSpecific_HLAFrame.insert(loc = 0, column = 'tempCol', value = ProteinSpecific_HLAFrame.apply(lambda x: epitope_density(x['Pos'], x['Peptide'], num), axis = 1))
            Epitope_Density.append(ProteinSpecific_HLAFrame['tempCol'].sum(axis = 0))
            HLAGroupby = ProteinSpecific_HLAFrame.groupby(['HLA'])['tempCol'].sum()
            if HLA_EpitopeDensity.empty:
                HLA_EpitopeDensity = HLAGroupby
            else:
                HLA_EpitopeDensity = pd.concat([HLA_EpitopeDensity, HLAGroupby], axis = 1)
            ProteinSpecific_HLAFrame.drop(['tempCol'], inplace = True, axis = 1)
            
        
        Mut_Positions = MutData_ProteinSpecific['position'].values
        Mut_TotalMutations = MutData_ProteinSpecific['Total_mutations'].values
        for num in range(1, Protein_length+1):
            tempMut_Index = np.where(Mut_Positions == num)
            print('This is tempMut_Index %s'% (tempMut_Index))
            Total_Mutations_atPosition = np.take(Mut_TotalMutations, tempMut_Index)
            print('This is Total_Mutations_atPosition %s'% (Total_Mutations_atPosition))
            print('len is %s, size is %s'% (len(Total_Mutations_atPosition), Total_Mutations_atPosition.size))
            if Total_Mutations_atPosition.size > 0:
                
                Max_TotalMutation_atPosition = np.amax(Total_Mutations_atPosition)
            else:
                
                Max_TotalMutation_atPosition = np.nan

            Mutation_Density.append(np.log2(Max_TotalMutation_atPosition))
            
            print(Mutation_Density)


        TCell_seq_Val = Validated_Data_AllInfo_Prot_Spec['Start'].values
        TCellMagnitude_Val = Validated_Data_AllInfo_Prot_Spec['Response freq. (overall)'].values

        


        for num in range(1, Protein_length+1):
            tempMut_Index = np.where(TCell_seq_Val == num)
            print('This is tempMut_Index %s'% (tempMut_Index))
            TCellMagn_atPosition = np.take(TCellMagnitude_Val, tempMut_Index)
            print('This is Total_Mutations_atPosition %s'% (TCellMagn_atPosition))
            print('len is %s, size is %s'% (len(TCellMagn_atPosition), TCellMagn_atPosition.size))
            if TCellMagn_atPosition.size > 0:
                print('hello, above 0')
                Max_TCellMagn_atPosition = np.amax(TCellMagn_atPosition)
            else:
                print('its zero!')
                Max_TCellMagn_atPosition = np.nan

            Max_TCellMagn_Val.append(Max_TCellMagn_atPosition) #####np.log2()
            #MutData_ProteinSpecific.insert(loc = 0, column = 'tempCol', value = MutData_ProteinSpecific.apply(lambda x: Mutation_density(x['position'], x['Total_mutations'], num), axis = 1))
            #Mutation_Density.append(MutData_ProteinSpecific['tempCol'].max(axis = 0))
            #MutData_ProteinSpecific.drop(['tempCol'], inplace = True, axis = 1)
            print(Max_TCellMagn_Val)




        TCell_seq_Val_Matched = Validated_Data_Matched_Prot_Spec['Start'].values
        TCellMagnitude_Val_Matched = Validated_Data_Matched_Prot_Spec['Response freq. (overall)'].values

        for num in range(1, Protein_length+1):
            tempMut_Index = np.where(TCell_seq_Val_Matched == num)
            print('This is tempMut_Index %s'% (tempMut_Index))
            TCellMagn_Matched_atPosition = np.take(TCellMagnitude_Val_Matched, tempMut_Index)
            print('This is Total_Mutations_atPosition %s'% (TCellMagn_Matched_atPosition))
            print('len is %s, size is %s'% (len(TCellMagn_Matched_atPosition), TCellMagn_Matched_atPosition.size))
            if TCellMagn_Matched_atPosition.size > 0:
                print('hello, above 0')
                Max_TCellMagn_Matched_atPosition = np.amax(TCellMagn_Matched_atPosition)
            else:
                print('its zero!')
                Max_TCellMagn_Matched_atPosition = np.nan

            Max_TCellMagn_Matched_Val.append(Max_TCellMagn_Matched_atPosition) #np.log2()
            #MutData_ProteinSpecific.insert(loc = 0, column = 'tempCol', value = MutData_ProteinSpecific.apply(lambda x: Mutation_density(x['position'], x['Total_mutations'], num), axis = 1))
            #Mutation_Density.append(MutData_ProteinSpecific['tempCol'].max(axis = 0))
            #MutData_ProteinSpecific.drop(['tempCol'], inplace = True, axis = 1)
            print(Max_TCellMagn_Matched_Val)


        if len(Mut_Positions) == 0:
            print('Mut_Positions empty')
            MutData_ProteinSpecific.insert(loc = 0, column = 'Epitope-association rate', value = 0)

        else:
            print('Mut_Positions NOT empty')
            MutData_ProteinSpecific.insert(loc = 0, column = 'Epitope-association rate', value = MutData_ProteinSpecific.apply(lambda x: Epitope_Association(x['position'], ProteinSpecific_HLAFrame), axis = 1))
        Epitope_Association_Rate = MutData_ProteinSpecific['Epitope-association rate'].sum(axis = 0) / len(MutData_ProteinSpecific.index)
        
        Mut_Epitope_Association_Rate[Protein_Name] = Epitope_Association_Rate

        
      
        Epitope_Density_df = DataFrame(Epitope_Density)
        Mutation_Density_df = DataFrame(Mutation_Density)
        Max_TCellMagn_Val_df = DataFrame(Max_TCellMagn_Val)
        Max_TCellMagn_Matched_Val_df = DataFrame(Max_TCellMagn_Matched_Val)
        print('THIS IS Max_TCellMagn_Val_df: %s'% (Max_TCellMagn_Val_df))

        allData = [Epitope_Density_df, Mutation_Density_df, Max_TCellMagn_Val_df, Max_TCellMagn_Matched_Val_df]
        for frame in allData:
            frame['position'] = frame.index
            frame['position'] = frame['position'].apply(lambda x: x + 1)
   
        HLA_EpitopeDensity = HLA_EpitopeDensity.T
        HLA_EpitopeDensity['position'] = Epitope_Density_df.index
        HLA_EpitopeDensity['position'] = HLA_EpitopeDensity['position'].apply(lambda x: x + 1)
        

        Mutation_Density_df.rename(columns = {0: 'Mutation Density'}, inplace = True)
        Max_TCellMagn_Val_df.rename(columns = {0: 'Mutation Density'}, inplace = True)
        Max_TCellMagn_Matched_Val_df.rename(columns = {0: 'Mutation Density'}, inplace = True)


        HLA_EpitopeDensity.set_index('position', inplace = True)

        HLA_EpitopeDensity['sum'] = HLA_EpitopeDensity.sum(axis = 1)
        
        

        HLA_EpitopeDensity = HLA_EpitopeDensity.stack().to_frame()
        HLA_EpitopeDensity.reset_index(inplace = True)
        HLA_EpitopeDensity.rename(columns = {0: 'HLADensity'}, inplace = True)
        

        
        if Protein_length not in AllProteinsData.keys():
            AllProteinsData[Protein_length] = [0, Epitope_Density_df, HLA_EpitopeDensity, Mutation_Density_df, Protein_Name, Max_TCellMagn_Val_df, Max_TCellMagn_Matched_Val_df] #CoverageMap
        else:
            AllProteinsData[Protein_length + 0.1] = [0, Epitope_Density_df, HLA_EpitopeDensity, Mutation_Density_df, Protein_Name, Max_TCellMagn_Val_df, Max_TCellMagn_Matched_Val_df] #CoverageMap
        #####Somethig like this: AllProteinsData[proteinName] = [CoverageMap, Epitope_Density_df, Mutation_Density_df]
        ##############################
    
    #AllPro_EpitopeSum = AllPro_EpitopeSum.T


    

    lengths = list(AllProteinsData.keys())
    print('\n\n\n\n\n\n\n\nTHIS IS PROTEIN LIST!!!!!!!!!!!: %s\n\n\n\n\n\n\n\n\n\n'% (lengths))
    for key, item in AllProteinsData.items():
        print(item[4])
    for item in lengths:
        print(AllProteinsData[item][4])
    lengths.sort(reverse = True)
    lengths1 = lengths[:3]
    lengths2 = lengths[3:7]
    lengths3 = lengths[7:]
    widths1 = []
    widths2 = []
    widths3 = []
    Largest1 = lengths1[0]
    Largest2 = lengths2[0]
    Largest3 = lengths3[0]

    for length in lengths1:
        widths1.append(length/Largest1)
    for length in lengths2:
        widths2.append(length/Largest2)
    for length in lengths3:
        widths3.append(length/Largest3)
    


    Make_Figure(AllProteinsData = AllProteinsData, Protein_Name = Protein_Name, lengths = lengths1, width = widths1, s1 = 80, s2 = 350, fig_size = (150, 8), linewidth = 50, fig_num = 1)
    Make_Figure(AllProteinsData = AllProteinsData, Protein_Name = Protein_Name, lengths = lengths2, width = widths2, s1 = 100, s2 = 400, fig_size = (35, 9.5), linewidth = 50, fig_num = 2)
    Make_Figure(AllProteinsData = AllProteinsData, Protein_Name = Protein_Name, lengths = lengths3, width = widths3, s1 = 100, s2 = 200, fig_size = (20, 8), linewidth = 50, fig_num = 3)   



    


'''
4. Function that takes in two dataframes, one with mutated peptides and one with their reference peptides (each pair is HLA and mutation specific).
Dataframes containg rank and nM binding affinities for mutated and reference peptides. 
Function will return rank difference as well as fold change for each peptide in dataframe.

inputs: two dataframes (one mutated and one reference) with relevent peptide info corresponding to a specific mutation and HLA.
Outputs: rank difference as well as fold change for each peptide relevant to a specific mutation/HLA combination.
'''


def Join_HLAs(Frames):
    combined_Frames = pd.DataFrame()
    for key, item in Frames.items():
        item = item[(item['Rank'].astype(float) <= 1.0)] #Further shortlist the peptides to keep only the mutated/reference pairs where one of the two is binding (rank below 2, according to netMHCpan)
        item['Peptide'] = item['Peptide'].apply(lambda x: '%s_%s'% (x, key))
        if combined_Frames.empty:
            combined_Frames = item
        else:
            combined_Frames = pd.concat([combined_Frames, item])



    return combined_Frames

'''
3. Function that splits each HLA dataframe into each unique mutation.
As a result, each HLA will have its own dictionary of mutations, where key = mutation ID and item = dataframe of corresponding mutated peptides and relevant info.

input: Dictionary with HLA-specific dataframes (from function 2)
outputs: Dictionary where key = mutations, value = list of lists, where each inner list = HLA type at pos 0, and dataframe of corresponding mutated peptides and relevant info at pos 1.
'''


def Split_into_HLAs(RawData, HLAs):
    HLA_Data_List = {}
    CommonInfo = RawData[['Pos','Peptide','ID']].copy() #gather info common to all HLA types into separate dataframe (peptide position, sequence, and mutation info)
    #print('Hello, %s'% (type(CommonInfo)))
    RawData.drop(['Pos','Peptide','ID', 'core', 'icore', 'EL-score', 'BA-score', RawData.columns[-2], RawData.columns[-1]], axis = 1, inplace = True) #Remove all info that will not be used, to simplyfy dataframe. '1-log50k',
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

            TempFrame.columns = ['nM', 'Rank'] #remane columns in temp dataframe containing rank and nM for that particular HLA
            ####
            TempFrame['HLA'] = HLAs[HLAnum]
            newTemp = pd.concat([CommonInfo, TempFrame], axis = 1) #Combine rank/nM dataframe of HLA with 'commoninfo' (pos, peptide and ID). As a result, each HLA will have its rank and nM binding affinities, as well as other important info (peptide pos, sequence, and ID)
            newTemp.set_index('ID', inplace = True) #Set ID as index
            HLA_Data_List[HLAs[HLAnum]] = newTemp  #.T #transpos HLA dataframe to facilitate iteration through mutations, and store in dictionary (key = HLA type, item = corresponding dataframe)
            HLAnum += 1 #Iterate through HLA type names in list created in previous function.
            #print('This is newTemp: \n%s'% (newTemp))
        #print(HLA_Data_List)
        #print(HLA_Data_List['HLA-A02:02'].loc[:, ['S_D623G_1']])
    return HLA_Data_List

def Import_results_csv(file):
    gain_loss = pd.read_csv(file)
    gain_loss = gain_loss.drop([gain_loss.columns[-3], gain_loss.columns[-2]], axis = 1)
    print(gain_loss)
    return gain_loss.reset_index()


def T_test(Results_file, Strong_SARSCOV2_epitopes, folder, stringency, effect):

    #folder = 'results'
    os.makedirs(os.path.dirname('./%s/'% (folder)), exist_ok=True)

    Combined_files = pd.DataFrame()


    for HLA, value in Results_file.items(): #Iterate through gain/loss results files (in other wors, iterate through HLAs queried).
        #reform_HLA = HLA[:-2] + ':' + HLA[-2] + HLA[-1] #reformat HLA name.
        #print(reform_HLA)
        print(HLA)
        temp_Frame = Import_results_csv(value[:-1]) #For each loss/gain file name, import corresponding csv file, which is saved in same directory.
        #print('\n\n\nThis is temp_Frame: %s\n\n\n'% (temp_Frame))
        temp_Frame.drop(['index', 'Unnamed: 0'], axis = 1, inplace = True) #reformat columns (drop unnecessary columns)
        #reduced_Frame = temp_Frame[['level_0', 'Total_mutations', 'Rank_Mut', 'Rank_Ref', 'Peptide_Mut', 'Peptide_Ref']].copy()
        
        if stringency == 'strong':
            loss_frame = temp_Frame[temp_Frame['effect on binding'] == effect] #'loss'
        if stringency == 'intermediate':
            loss_frame = temp_Frame[temp_Frame['effect on binding_SemiConservative'] == effect] #'loss'
        if stringency == 'no effect':
            loss_frame = temp_Frame[temp_Frame['effect on binding_nonConservative'] == 'no effect']
        #Intermediate_loss_frame = temp_Frame[temp_Frame['effect on binding_SemiConservative'] == 'loss']

        #temp_SARSCOV2_epitopes = SARSCOV2_epitopes[SARSCOV2_epitopes['SARS-CoV-2 based HLA alleles'] == HLA]
        
        if Combined_files.empty:
            Combined_files = loss_frame
        else:
            Combined_files = Combined_files.append(loss_frame)
    Matched_epitopes = pd.merge(Combined_files, Strong_SARSCOV2_epitopes, how = 'right', on = ['Peptide_Ref'])
    Matched_epitopes.dropna(subset = ['Peptide_Mut'], inplace = True)

    count = Matched_epitopes['Peptide_Mut'].nunique()


    Matched_epitopes.to_csv('%s/%s_%s_All_Mutations_detail.csv'% (folder, effect, stringency))

    #Strong_SARSCOV2_epitopes.to_csv('%s/%s_%s_SARSCOV2_epitopes.csv'% (folder, effect, stringency))
    print('This is count: %s' %(count))

    return Matched_epitopes

def Import_results_Validated_csv(file):
    gain_loss = pd.read_csv(file)
    print(gain_loss)
    return gain_loss.reset_index()





def reform_HLAcolumn(x):
    #print(pd.isnull(x))
    if not pd.isnull(x):
        x = x.replace('*', '')
        x = x.replace(':', '')
    else:
        x = 'unknown'
    return x

'''
1. Function to read netMHCpan output (in .csv format), and create a list with all HLA types from file.
'''

def ReadFile(FileName):
    Raw_data = pd.read_csv(FileName) #Read csv file
    preHLAs = list(Raw_data.columns.values) #create list of HLA types from file
    HLAs = []
    for item in preHLAs:
        if item[:3].upper() == 'HLA':
            HLAs.append(item)
    Raw_data.columns = list(Raw_data.iloc[0])
    Raw_data.drop([0], inplace = True)
    Raw_data.reset_index(drop = True, inplace = True)
    #print(list(Raw_data.iloc[0]))
    #print(Raw_data)

    HLAspecific_Data_Dict = Split_into_HLAs(Raw_data, HLAs)
    ProcessedHLAFrame = Join_HLAs(HLAspecific_Data_Dict)

    while True:
        Protein_Seq = sys.argv[2]
        #Protein_Seq = input('\n\nPlease, enter txt file (Protein sequence) in .csv format\nEx. ExampleFile.csv\nWhen done, press Enter: ')
        MutFrequency = sys.argv[3]
        #MutFrequency = input('\n\nPlease, enter MutFrequency file in .csv format\nEx. ExampleFile.csv\nWhen done, press Enter: ')
        Validated_Epitopes = sys.argv[4]
        #Validated_Epitopes = input('\n\nPlease, enter nValidated_Epitopes .csv format\nEx. ExampleFile.csv\nWhen done, press Enter: ')
        #Validated_peptides_Matched = input('Please input validated Matched peptide list: ')
        prompt = sys.argv[5]
        #prompt = input('\n\nPlease, enter the text filename contaning the anme of all HLA-specific gain/loss csv files (processed netMHCpan files), press Enter: ')

        try:
            RefProteins = {}

            with open(Protein_Seq, "rU") as File: #open file containing reference protein sequences, store in  dictionary (key = protein name, item = protein sequence)
                for record in SeqIO.parse(File, "fasta"):
                    RefProteins[record.id] = record.seq

            MutData = pd.read_csv(MutFrequency)
            #Validated_Data = pd.read_csv(Validated_Epitopes)
            #Validated_Data_Matched = pd.read_csv(Validated_peptides_Matched)
            #print('This is Validated_Data_Matched: %s:' %(Validated_Data_Matched['Mut_Protein']))
            
            val_Ep_File = Import_results_Validated_csv(Validated_Epitopes) #For each loss/gain file name, import corresponding csv file, which is saved in same directory.
            val_Ep_File.drop(['index'], axis = 1, inplace = True) #reformat columns (drop unnecessary columns)
                        
            #####val_Ep_File['SARS-CoV-2 based HLA alleles'] = val_Ep_File.apply(lambda x: reform_HLAcolumn(x['SARS-CoV-2 based HLA alleles']), axis = 1)  #.str
            #####val_Ep_File.rename(columns = {'SARS-CoV-2 based HLA alleles': Validated_Epitopes}, inplace = True)
            val_Ep_File.rename(columns = {'SARS-CoV-2 derived epitope': 'Peptide_Ref'}, inplace = True)

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

    Validated_Data_Matched = T_test(Results_file, val_Ep_File, 'HLA_MATCHED_VALIDATED_EPITOPES', 'strong', 'loss')

    Create_Mutation_Dict(ProcessedHLAFrame, HLAspecific_Data_Dict, RefProteins, MutData, val_Ep_File, Validated_Data_Matched) #






#FileName = 'Germany_GENOMES.csv'

print(f'This is sys.argv: {sys.argv}')

FileName = sys.argv[1]

ReadFile(FileName)

'''
while True:
    
    FileName = sys.argv[1]
    #FileName = input('\n\nPlease, enter netMHCpan output file (protein-specific epitopes) in .csv format\nEx. ExampleFile.csv\nWhen done, press Enter: ')
    try:
        ReadFile(FileName)
    except FileNotFoundError:
        print('This file does not exist')
    else:
        break
'''


