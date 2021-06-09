import pandas as pd
import numpy as np
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy.stats import ttest_1samp
from scipy import stats
from scipy.stats import ttest_ind
import sys
import os

sns.set_style(style = 'ticks')

plt.style.use('seaborn-white')

def MutFrequency_file():
    while True:
        prompt = sys.argv[2]
        try:
            MutFrequency_Files = {}
            with open(prompt, 'r') as File:
                for line in File:
                    MutFrequency_Files[line[:line.find('_')]] = line
                return (MutFrequency_Files, prompt)
        except FileNotFoundError:
            print('file not found. Make sure to copy in in same path as this script.')
        else:
            break

def Import_MutFrequency(file, type):
    if type == 'neutral':
        Country_1 = pd.read_csv(f'{sys.argv[3]}/{file}')
    else:
        Country_1 = pd.read_csv(file)
    return Country_1.reset_index()

'''

def Get_Frame(MutFrequency_Arrays):
    Compare_Countries_Frame = pd.DataFrame()
    for key, value in MutFrequency_Arrays.items():
        MutID_Dict = {}
        for key2, value2 in MutFrequency_Arrays.items():
            if key == key2:
                MutID_Dict[key2] = 0
            else:
                a = len(np.intersect1d(value, value2))
                b = len(np.setdiff1d(value, value2))
                MutID_Dict[key2] = (a/(a + b))*100

        Frame_Column = pd.DataFrame.from_dict(MutID_Dict, orient = 'index')
        Frame_Column.rename(columns = {0: key}, inplace = True)

        print(Frame_Column)
        
        if len(Compare_Countries_Frame.columns) == 0:
            Compare_Countries_Frame = Frame_Column
        else:
            Compare_Countries_Frame = pd.concat([Compare_Countries_Frame, Frame_Column], axis = 1)

    print(Compare_Countries_Frame)
    return Compare_Countries_Frame
'''


def printbarGraphs(Ref_Res, Mut_Res, Country):
    fig, axes = plt.subplots(2, 1, figsize = (14, 12))
    sns.set(font_scale = 1.5)
    sns.barplot(x = 'reference residue rate', y = 'reference residue', data = Ref_Res, color = 'green', ax = axes[0])
    sns.barplot(x = 'mutated residue rate', y = 'mutated residue', data = Mut_Res, color = 'green', ax = axes[1])
    axes[0].set_xlabel('')
    axes[0].set_xlabel('')

    plt.grid()
    plt.tight_layout()
    #####plt.savefig('%s_Total_Mutation_Analysis.png'% (Country))

def print_Lineplot(All_Samples, folder):
    with sns.axes_style("ticks"):
        fig, ax = plt.subplots(figsize = (60, 25))  #17
        sns.set(font_scale = 2)
        sns.lineplot(x = 'residue', y = 'Mutation Difference', hue = 'Sample', data = All_Samples, linewidth = 3)

        ax.set_ylabel('Overall Residue Output', fontsize = 50)
        ax.tick_params(axis='y', labelsize=60)
        ax.set_xlabel('Residue', fontsize = 50)
        ax.tick_params(axis='x', labelsize=60)

        leg = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        for legobj in leg.legendHandles:
            legobj.set_linewidth(4.0)

        plt.grid()
        plt.tight_layout()
        plt.savefig('./%s/MutationalPatterns.png'% (folder))

def print_Lineplot_allreplicates(All_Samples, folder):

    samples = All_Samples['Sample'].unique()
    colors = ['blue', 'orange', 'green', 'red', 'brown', 'purple', 'yellow']
    line_dash = ['dotted', 'dashed', 'dashdot', 'solid']

    with sns.axes_style("ticks"):
        fig, ax = plt.subplots(figsize = (16, 7))  #17
        sns.set(font_scale = 2)
        colNum = 0
        for sample in samples:
            
            temp_Frame = All_Samples[All_Samples['Sample'] == sample]
            
            
            
            if 'Mutations' in sample:
                print('HEEEEELLLLLLLLOOOOOOOOOOOOOOOO')
                sns.lineplot(x = 'residue', y = 'Mutation Difference', hue = 'Countries', style = 'Countries', data = temp_Frame, palette = (colors[colNum],), linewidth = 2) # , dashes = ['--']  , legend = False
                colNum += 1
            else:
                sns.lineplot(x = 'residue', y = 'Mutation Difference', hue = 'Countries', data = temp_Frame, palette = ('black',), linewidth = 2, legend = False) #style = 'Countries', dashes = (line_dash[colNum],), 
            

        ax.set_ylabel('Overall Residue Output', fontsize = 30)
        ax.tick_params(axis='y', labelsize=25)
        ax.set_xlabel('Residue', fontsize = 30)
        ax.tick_params(axis='x', labelsize=25)

        leg = plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)

        plt.grid()
        plt.tight_layout()
        plt.savefig('./%s/Patterns_ALLREPLICATES.pdf'% (folder), dpi = 300)
        plt.savefig('./%s/Patterns_ALLREPLICATES.png'% (folder))


def OneSample_Ttest(One_Sample_Ttest_Dict, folder):

    print(f'This is One_Sample_Ttest_Dict: {One_Sample_Ttest_Dict}')

    samples = {}
    observed = {}

    for sampleName, sample in One_Sample_Ttest_Dict.items():
        if 'Neutral' in sampleName:
            Neutral = sample.T
        elif len(sample.columns.values) == 1:
            observed[sampleName] = sample.T
            #print('THIS IS OBSERVED SAMPLE: %s'% (observed))
        else:
            samples[sampleName] = sample.T

    print(f'This is Neutral: {Neutral}')
    print(f'This is observed: {observed}')
    print(f'This is samples: {samples}')

    T_test_aminoacids = {}
    for sampleName, sample in samples.items():
        aminoAcids = sample.columns.values
        for amino_acid in aminoAcids:
            temp_sample = sample[amino_acid].to_numpy()
            temp_Neutral = Neutral[amino_acid].to_numpy()
            t, p = stats.ttest_ind(temp_sample,temp_Neutral, equal_var = False)
            print(sampleName)
            print(amino_acid)
            T_test_aminoacids['%s_%s'% (amino_acid, sampleName[:-1])] = [t,p,temp_sample[0],np.average(temp_Neutral),np.std(temp_Neutral),amino_acid]

    for sampleName, sample in observed.items():
        aminoAcids = sample.columns.values
        for amino_acid in aminoAcids:
            temp_Neutral = Neutral[amino_acid].to_numpy()
            temp_observed = sample[amino_acid].to_numpy()
            t, p = ttest_1samp(temp_Neutral, popmean = temp_observed[0])
            T_test_aminoacids['%s_%s'% (amino_acid, sampleName[:-1])] = [t,p,temp_observed[0],np.average(temp_Neutral),np.std(temp_Neutral),amino_acid]


    

    T_Test_results = pd.DataFrame.from_dict(T_test_aminoacids, orient = 'index')
    T_Test_results.rename(columns = {0: 'TTest', 1: 'p-val', 2: 'GRSO', 3: 'Neutral average', 4: 'Neutral stdv', 5: 'Amino Acid'}, inplace = True)
    T_Test_results.to_csv('./%s/T_TEST_RESULTS.csv'% (folder))

    for sample, observed in T_test_aminoacids.items():
        print(sample)
        print(observed[0])
        print(observed[1])
    

MutFreq_return = MutFrequency_file()
MutFrequency_Files = MutFreq_return[0]
folder = 'RESULTS' #MutFreq_return[1]
os.makedirs(os.path.dirname('./%s/'% (folder)), exist_ok=True)
print(MutFrequency_Files)
MutFrequency_Arrays_ALL = {}
MutFrequency_Arrays_Common = {}
MutFrequency_Arrays_Rare = {}
MutationType_Frame = {}
Samples_Frame = {}
sample = ''

for key, value in MutFrequency_Files.items():
    if key != '':
        if 'csv' not in value:
            sample = value
            MutationType_Frame = {}
        else:
            file = Import_MutFrequency(value[:-1], 'neutral')
            MutationType_Frame[key] = file
            MutFrequency_Arrays_ALL[key] = file['MutID'].to_numpy()
            Samples_Frame[sample] = MutationType_Frame


##Import mutation list csv as datafrma, segragate into frequency segments, then save each sub-dataframe as MutationType_Frame[frequency group] = dataframe, then Samples_Frame[frequency group] = MutationType_Frame
Mutation_List = Import_MutFrequency(sys.argv[1], 'samples')
Mutation_List_1000UP = Mutation_List[ Mutation_List['Total_mutations'] > 1000  ]
Mutation_List_100to1000UP = Mutation_List[ (Mutation_List['Total_mutations'] > 100) &  (Mutation_List['Total_mutations'] < 1000) ]
Mutation_List_2to100UP = Mutation_List[ (Mutation_List['Total_mutations'] > 2) &  (Mutation_List['Total_mutations'] < 100) ]
Mutation_List_Singletons = Mutation_List[ Mutation_List['Total_mutations'] < 2  ]

print(f'This is 1000UP: \n{Mutation_List_1000UP}\n\n100to1000UP: {Mutation_List_100to1000UP}\n\n2to100UP: {Mutation_List_2to100UP}\n\nsingletons: {Mutation_List_Singletons}')

MutationType_Frame = {}



if not Mutation_List_Singletons.empty:
    MutationType_Frame['Singletons'] = Mutation_List_Singletons
    Samples_Frame['Mutations in Singletons'] = MutationType_Frame
    MutationType_Frame = {}

if not Mutation_List_100to1000UP.empty:
    MutationType_Frame['100to1000'] = Mutation_List_100to1000UP
    Samples_Frame['Mutations in 100to1000'] = MutationType_Frame
    MutationType_Frame = {}

if not Mutation_List_1000UP.empty:
    MutationType_Frame['1000UP'] = Mutation_List_1000UP
    Samples_Frame['Mutations in 1000UP'] = MutationType_Frame
    MutationType_Frame = {}

if not Mutation_List_2to100UP.empty:
    MutationType_Frame['2to100'] = Mutation_List_2to100UP
    Samples_Frame['Mutations in 2to100'] = MutationType_Frame
    MutationType_Frame = {}




########


All_Samples = pd.DataFrame()
Ref_count_mean_allSamples = pd.DataFrame()
Mut_count_mean_allSamples = pd.DataFrame()
One_Sample_Ttest_Dict = {}


for sample, replicates in Samples_Frame.items():
    print(f'this is (Samples_Frame) !!!!!!!!!!\n\n\n{key}\n\n\n\n\n{value}')
    
    All_countries_frame = pd.DataFrame()
    Ref_count_Allreps = pd.DataFrame()
    Mut_count_Allreps = pd.DataFrame()
    All_countries_frame_Ref = pd.DataFrame()
    All_countries_frame_Mut = pd.DataFrame()
    first = 1

    for key, value in replicates.items():
        print(f'this is (replicates) !!!!!!!!!!\n\n\n{key}\n\n\n\n\n{value}')
        value.insert(loc = 1, column = 'reference residue', value = value['MutID'].apply(lambda x: x[(x.find('_') + 1)]))
        value.insert(loc = 2, column = 'mutated residue', value = value['MutID'].str[-1])
        print(value)

        total_mutations = len(value.index)

        Ref_resid_count = value.groupby(['reference residue'])['Mutation_Rates'].count().reset_index(name = 'count') #['Mutation_Rates'] name = 'count'
        Ref_resid_count.insert(loc = 1, column = 'reference residue rate', value = Ref_resid_count['count'].apply(lambda x: (x/total_mutations)*100))
        Ref_count_getAverage = Ref_resid_count.reset_index()[['reference residue', 'count']].copy()
        Ref_count_getAverage.rename(columns = {'count': key}, inplace = True)
        Ref_resid_count.drop(['count'], axis = 1, inplace = True)
        
        if Ref_count_Allreps.empty:
            Ref_count_Allreps = Ref_count_getAverage
        else:
            Ref_count_Allreps = pd.merge(Ref_count_Allreps, Ref_count_getAverage, on = ['reference residue'])


        Mut_resid_count = value.groupby(['mutated residue'])['Mutation_Rates'].count().reset_index(name = 'count') #['Mutation_Rates'] name = 'count'
        Mut_resid_count.insert(loc = 1, column = 'mutated residue rate', value = Mut_resid_count['count'].apply(lambda x: (x/total_mutations)*100))
        Mut_count_getAverage = Mut_resid_count.reset_index()[['mutated residue', 'count']].copy()
        Mut_count_getAverage.rename(columns = {'count': key}, inplace = True)
        Mut_resid_count.drop(['count'], axis = 1, inplace = True)

        if Mut_count_Allreps.empty:
            Mut_count_Allreps = Mut_count_getAverage
        else:
            Mut_count_Allreps = pd.merge(Mut_count_Allreps, Mut_count_getAverage, on = ['mutated residue'])
        

        printbarGraphs(Ref_Res = Ref_resid_count, Mut_Res = Mut_resid_count, Country = key)

        Ref_resid_count.rename(columns = {'reference residue' : 'residue'}, inplace = True)
        Mut_resid_count.rename(columns = {'mutated residue' : 'residue'}, inplace = True)

        print('%s, \n%s'% (key, Ref_resid_count))
        print('%s, \n%s'% (key, Mut_resid_count))

        Combined_frame = pd.merge(Ref_resid_count, Mut_resid_count, how = 'left', on = ['residue'])
        
        Combined_frame['%s'% (key)] = Combined_frame['mutated residue rate'] - Combined_frame['reference residue rate']

        Combined_frame.drop(['reference residue rate', 'mutated residue rate'], axis = 1, inplace = True)

        Ref_resid_count.rename(columns = {'reference residue rate' : '%s'% (key)}, inplace = True)
        Mut_resid_count.rename(columns = {'mutated residue rate' : '%s'% (key)}, inplace = True)

        if first == 1:
            All_countries_frame = Combined_frame
            All_countries_frame_Ref = Ref_resid_count
            All_countries_frame_Mut = Mut_resid_count
            first += 1
        else:
            All_countries_frame = pd.merge(All_countries_frame, Combined_frame, how = 'left', on = ['residue'])
            All_countries_frame_Ref = pd.merge(All_countries_frame_Ref, Ref_resid_count, how = 'left', on = ['residue'])
            All_countries_frame_Mut = pd.merge(All_countries_frame_Mut, Mut_resid_count, how = 'left', on = ['residue'])

    print('\nTHIS IS Ref_count_Allreps and Mut_count_Allreps: %s\n\n%s\n\n\n'% (Ref_count_Allreps, Mut_count_Allreps))

    Ref_count_Allreps['mean'] = Ref_count_Allreps.mean(axis = 1, numeric_only = True)
    Ref_count_Allreps.rename(columns = {'mean': sample[:-1]}, inplace = True)
    Mut_count_Allreps['mean'] = Mut_count_Allreps.mean(axis = 1, numeric_only = True)
    Mut_count_Allreps.rename(columns = {'mean': sample[:-1]}, inplace = True)

    print('\nTHIS IS Ref_count_Allreps and Mut_count_Allreps: %s\n\n%s\n\n\n'% (Ref_count_Allreps, Mut_count_Allreps))
    
    Ref_res_count_mean = Ref_count_Allreps[['reference residue', sample[:-1]]].copy()
    print(Ref_res_count_mean)
    Mut_res_count_mean = Mut_count_Allreps[['mutated residue', sample[:-1]]].copy()

    

    if Ref_count_mean_allSamples.empty:
        Ref_count_mean_allSamples = Ref_res_count_mean
        Mut_count_mean_allSamples = Mut_res_count_mean
    else:
        Ref_count_mean_allSamples = pd.merge(Ref_count_mean_allSamples, Ref_res_count_mean, on = ['reference residue'])
        Mut_count_mean_allSamples = pd.merge(Mut_count_mean_allSamples, Mut_res_count_mean, on = ['mutated residue'])

    All_countries_frame.set_index('residue', inplace = True)

    All_countries_frame_Ref.set_index('residue', inplace = True)
    All_countries_frame_Ref = All_countries_frame_Ref.stack().reset_index(name = 'reference residue rate')
    All_countries_frame_Ref.rename(columns = {'level_1': 'Countries'}, inplace = True)

    All_countries_frame_Mut.set_index('residue', inplace = True)
    All_countries_frame_Mut = All_countries_frame_Mut.stack().reset_index(name = 'Mutated residue rate')
    All_countries_frame_Mut.rename(columns = {'level_1': 'Countries'}, inplace = True)
    

    All_countries_frameMut_Ref = pd.merge(All_countries_frame_Ref, All_countries_frame_Mut, how = 'left', on = ['residue', 'Countries'])
    
    print(All_countries_frame_Ref)
    print(All_countries_frame_Mut)
    print(All_countries_frameMut_Ref)

    One_Sample_Ttest_Dict[sample] = All_countries_frame
    #All_countries_frame.to_csv('%s_All_countries_frame.csv'% (sample))
    All_countries_frame = All_countries_frame.stack().reset_index(name = 'Mutation Difference')
    All_countries_frame.rename(columns = {'level_1': 'Countries'}, inplace = True)
    All_countries_frame['Sample'] = sample

    if All_Samples.empty:
        All_Samples = All_countries_frame
    else:
        All_Samples = pd.concat([All_Samples, All_countries_frame])

All_Samples['Mutation Difference'] = pd.to_numeric(All_Samples['Mutation Difference'], errors='coerce')
print(All_Samples.dtypes)
print(Samples_Frame.keys())


print_Lineplot(All_Samples, folder)
print_Lineplot_allreplicates(All_Samples, folder)
OneSample_Ttest(One_Sample_Ttest_Dict, folder)













