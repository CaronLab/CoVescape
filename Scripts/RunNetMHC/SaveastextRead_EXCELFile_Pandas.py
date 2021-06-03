import pandas as pd
import numpy as np
import sys
import time


######This code takes in all the netMHCpan output files in .xls format (each containing as many as 400 mutations), converts them to .csv and aggregates them all into one large .csv file

print(f'HELLOOOOOOO, this is arguments: {sys.argv}')

starttime = time.time()

Combined_Files = pd.DataFrame()

xls_files = []
with open(sys.argv[1], "rU") as File:
    lines = File.readlines()
    
for line in lines:

    line = line[:-1]
    

    if line[-4:] == '.xls':
        xls_files.append(line)
        

filenum = 0

all_files = []
ColNames = []

for filename in xls_files:
    print(f'This is filenum: {filenum}')
    if filenum > 29: #<
        filenum += 1
        continue
    elif filenum > 51:
        break

    with open(filename, "rU") as File:
        lines = File.readlines()
    print(f'Opened file {filename}')
    count = 0
    listOfList = []
    linenum = 0
    for line in lines:
        
        if filenum > 0: #29
            if linenum == 0:
                linenum += 1
                continue
            if linenum == 1:
                linenum += 1
                continue


        newline = line.replace('\t', ',')
        newline = newline[:-1]  
        
        newline = newline.split(',')
        if filenum == 0: #29
            if linenum == 0:
                #newline.append('')
                ColNames = newline
                #continue

        
        
        listOfList.append(newline) #.split(',')
        linenum += 1
    array = np.array(listOfList)
    print(f'this is array: {array}')
    
    print(f'This is filenum and filename: {filenum}, {filename}')
    

    filenum += 1

    all_files.append(array)

    

if len(all_files) > 1:
    Allarrays = np.vstack(all_files)
    Combined_Files = pd.DataFrame(Allarrays)
else:
    Combined_Files = pd.DataFrame(array)
print('We made it past the loop!')


Combined_Files.set_index(0, inplace = True)
new_colnames = Combined_Files.iloc[0]
Combined_Files = Combined_Files[1:]
Combined_Files.columns = new_colnames

Combined_Files.index.name = None

print(Combined_Files)



Combined_Files.to_csv('netMHCpanOUTPUT.csv')




