#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:11:49 2020

@author: Guillaume Petit (guillaume.petit@griffithuni.edu.au)
"""


'''
This Script opens formatted fasta files, count the number of cysteines in this sequence perform operations
on the sequences to generate a histogram of the distribution of cysteines in the extra-cytoplasmic and cytoplasmic 
proteome that is given

'''


### import of libraries used in the script
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 




##opening fast file containing all the proteins of Bps mature core genome (protein with signal sequences and all core genome)
           
path_matureSequence='/Path/to/Mature/sequences.fasta' ## mature sequences 
path_allSequences='/Path/to/all/sequences.fasta'  ## complete core genome

"""
 insert you path to data for the extra cytoplasmic protein and for the cytoplasmic ones\
reads all the sequences concatenated in Fasta file starting with '>' before each protein and ending with a stop codon "*"
For the extra cytoplasmic proteins, the files that come from SignalP5.0 concatenated together work perfeclty
(Submit sequences, download the results as fasta and use the "cat" commande  to assemble them all in one file )

"""

## putting all the mature sequences in a dataframes
## some tidying up required to have the sequence, gene name and identifier
## in different cells
prot_listMature=pd.read_table(path_matureSequence,header=None, names=['identifier'])## read the information 
prot_listMature['sequences']=prot_listMature['identifier'][1::2] ## every second cell contains the sequence and every other the descriptor
prot_listMature['sequences'][0::2]=np.nan
prot_listMature['identifier'][1::2]=np.nan
prot_listMature=prot_listMature.apply(lambda x: pd.Series(x.dropna().values))##remove empty cells



##importing fasta block of text into another dataframe
## again some tidying up required to have the right information in the right cell

text=open(path_allSequences,'r') ## start by opening the file
text_list=text.read() ## convert it to a list easier to handle 
text.close()
##separate the different gene blocks 
prot_listAll=pd.DataFrame(text_list.split('*\n'), columns=['col1'])## dataframe with the whole proteome

## populate different column with gene identifier/ name and sequence 
prot_listAll=pd.DataFrame(prot_listAll.col1.str.split(r'\s',2).tolist(), columns = ['identifier','name','sequences'])
prot_listAll.drop(prot_listAll.tail(1).index,inplace=True)## last row is empty





## counting the cysteines and putting the final number in a new column call number cys \
## this should work regardless of the indexing 

## for mature proteins 
for i in prot_listMature.index: 
   prot_listMature.loc[i,'Cys']=prot_listMature.loc[i,'sequences'].count("C")


## for all the proteins    
for i in prot_listAll.index:
   prot_listAll.loc[i,'Cys']=prot_listAll.loc[i,'sequences'].count("C")   

## comparing the two dataframe keeping only the difference, 
## effectively removing the mature sequences from the whole proteome to get cytoplasmic proteins
prot_listCyto=prot_listAll.merge(prot_listMature, on='identifier', how='outer', indicator=True).loc[lambda x : x['_merge']=='left_only']


##trimming the list a little, removing unused columns
prot_listCyto['Cys']=prot_listCyto['Cys_x']
prot_listCyto=prot_listCyto.drop(['Cys_x','Cys_y','sequences_y','_merge'], axis=1)


## generating a dataframe with all the proteins (mature + cytoplasmic)
## but using the sequences without the SP for the non cytoplasmic proteins
## Since the Signal peptide has not been removed from the whole proteome list 
## we actually need to replace these proteines with correpsonding mature sequences 
## in case there is a cysteines in the signal peptide

prot_listAll_noSP=pd.concat([prot_listMature, prot_listCyto]).drop(['name','sequences','sequences_x'],axis=1)

## Creating N bins to count each proteins with i cysteines (i between 0 and N-1) in the whole proteome
## no protein has more than 73 cysteines 
normal_base=pd.DataFrame({'index':np.arange(0,74),'Cys':np.zeros(74)})
normal_base['Cys']=prot_listAll_noSP['Cys'].value_counts()


## same operation for extra cytoplasmic proteins 
extra_cys=pd.DataFrame({'index':np.arange(0,74),'Cys':np.zeros(74)})
extra_cys['Cys']=prot_listMature['Cys'].value_counts()


## calculating the normalised frequency 
normal_extra=extra_cys['Cys']/normal_base['Cys']
normal_extra=normal_extra.fillna(value=0)
#print(normal_extra)


## simple plot, can be customised in the I python consol
## plot the normalised distribution
plt.figure(0)
plt.plot(normal_extra[0:15])
 ## Limits ofthe graph (0:15) can be changed but remember that there are very few proteins with more than 
## 10-12 cysteines meaning that at some point each single protein will generate a peak  



##saving the dataframes in Excel format with the commande to_excel('/path/to/excel/file.xlsx')
#prot_listCyto.to_excel('/path/to/file.xlsx') ## forlist of cytoplasmic prot
#prot_listMature.to_excel('/path/to/file.xlsx')## for mature proteins 


## plotting the number of cysteines vs frequency vs subcellular localisation 
plt.figure(1)
plt.hist((prot_listMature['Cys']),bins=np.arange(0,25),histtype='bar',rwidth=0.9,align='left')
plt.hist(prot_listCyto['Cys'],bins=np.arange(0,25),histtype='step', rwidth=0.9, align='left')
plt.xlabel('Number of cysteines')
plt.ylabel('Count')
plt.xticks(np.arange(0,25,1)) ## Only the 25 first bins are plotted for visibility 
plt.legend()
plt.show()




