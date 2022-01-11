# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""



''' Data sources:
    https://health-infobase.canada.ca/covid-19/epidemiological-summary-covid-19-cases.html
    https://health-infobase.canada.ca/covid-19/epidemiological-summary-covid-19-cases.html#VOC
    pulled as of Jan 10, 2022'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from datetime import date
from datetime import timedelta


lst_VOC = ['Alpha', 'Beta', 'Delta', 'Gamma', 'Omicron']


# working with Canada-wide data
df_cases = pd.read_csv('covid19-download.csv')
df_cases = df_cases.loc[df_cases['prname']=='Canada',:]
for i in df_cases.index:
    df_cases.loc[i,'date'] = date.fromisoformat(df_cases.loc[i,'date'])
df_variants = pd.read_csv('covid19-epiSummary-variants.csv')
# shift variant data from Sun to Sat to align with other data
for r in range(len(df_variants)):
    yr = int(df_variants.loc[r,'Collection (week)'].split('/')[2])
    m = int(df_variants.loc[r,'Collection (week)'].split('/')[0])
    d = int(df_variants.loc[r,'Collection (week)'].split('/')[1])
    df_variants.loc[r,'Collection (week)'] = date(yr,m,d)
df_vaccination = pd.read_csv('vaccination-coverage-map.csv')
df_vaccination = df_vaccination.loc[df_vaccination['prename']=='Canada',:]



# calculate cases by variant by week

df_variants = df_variants.groupby(['_Identifier','Collection (week)'])['%CT Count of Sample #'].sum()
df_variants = df_variants.unstack(level=0).fillna(0)

df_cases['numcases'] = df_cases['numconf'] + df_cases['numprob'].fillna(0)
df_cases['numpast'] = df_cases['numdeaths'] + df_cases['numrecover'].fillna(0)
df_varCases = pd.merge(df_cases[['date', 'numcases', 'numpast']], df_variants, 
                       left_on='date', right_on='Collection (week)')
df_varCases['newcases'] = df_varCases['numcases'].diff()
df_varCases.loc[0,'newcases'] = df_varCases.loc[0,'numcases']
df_varCases['recentpast'] = df_varCases['numpast'].diff()
df_varCases.loc[0,'recentpast'] = df_varCases.loc[0,'numpast']

lst_vars = df_variants.columns.values
for col in lst_vars:
    df_varCases[col] = df_varCases['newcases'] * df_varCases[col]
    
df_varActive = pd.DataFrame(index=df_varCases.index, columns=lst_vars)
r, r2clear = 0, 0 #FIFO
sr_varsActive = df_varCases.loc[r,lst_vars] * 0
sr_vars2clear = df_varCases.loc[r2clear,lst_vars]
while r < len(df_varCases):
    sr_varsActive = sr_varsActive + df_varCases.loc[r,lst_vars]
    num2clear = df_varCases.loc[r,'recentpast']
    numOldest = sr_vars2clear.sum()
    while numOldest <= num2clear:
        sr_varsActive = sr_varsActive - sr_vars2clear
        num2clear -= numOldest
        r2clear += 1
        sr_vars2clear = df_varCases.loc[r2clear,lst_vars]
        numOldest = sr_vars2clear.sum()
    sr_varsActive = sr_varsActive - sr_vars2clear * num2clear / numOldest
    sr_vars2clear = sr_vars2clear * (numOldest - num2clear) / numOldest
    df_varActive.iloc[r] = sr_varsActive
    for v in lst_vars:
        df_varActive.loc[df_varActive.index[r],v] = round(df_varActive.loc[df_varActive.index[r],v])
    r += 1
df_varActive.index = df_varCases['date']



# plot
'''activecases_variants = {}
for var in lst_vars:
    activecases_variants.update({var: df_varActive[var].values})
#print(df_varActive.index.values)'''
fig, ax = plt.subplots()
plt.rcParams['figure.figsize'] = [30,20]
plt.rcParams.update({'font.size': 35})
colours = []
cmap_blues = plt.get_cmap('GnBu')
cmap_purples = plt.get_cmap('PuRd')
for c in range(len(lst_vars)):
    if lst_vars[c] in lst_VOC:
        colours.append(cmap_blues(1 - c/len(lst_vars)))
    else:
        colours.append(cmap_purples(1 - c/len(lst_vars)))
        
ax.stackplot(df_varActive.index,
             df_varActive['Alpha'], df_varActive['B.1.1.318'], df_varActive['B.1.617.3'], 
             df_varActive['Beta'], df_varActive['Delta'], df_varActive['Eta'], 
             df_varActive['Gamma'], df_varActive['Iota'], df_varActive['Lambda'], 
             df_varActive['Mu'], df_varActive['Omicron'],
             df_varActive['Other'], df_varActive['Theta'],
             labels=lst_vars, colors=colours)
ax.legend(loc='upper left')
ax.set_title('Coronavirus Waves by Variant [Canada]')
ax.set_xlabel('Date', fontsize=35)
ax.set_ylabel('Active Cases', fontsize=35)

#plt.show()
plt.savefig('VariantActivity.png')
























