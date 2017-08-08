# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 09:12:11 2016

@author: pyao1
"""

import pandas as pd
import numpy as np
import glob
import os

# os.chdir(r'C:/Users/pyao1/python/n003-00500/uplc_csv/')

combined = pd.DataFrame()

for filename in glob.glob('*.txt'):
    data = pd.read_csv(filename,delimiter='\t')
    data['Result'] = filename
    combined = combined.append(data)

combined.sort_values(['Retention Time','Result'],inplace=True,kind='mergesort')

combined.dropna(inplace=True)

combined.index = np.arange(len(combined))

n_1 = combined['Retention Time'][1:]
n = combined['Retention Time'][:-1]
n_1.index = n.index

combined['Diff'] = n_1 - n
combined.set_value(combined.index[-1],'Diff',combined['Diff'].values[-2])

thresh = 0.005

low_edge = combined['Retention Time'].values[-0] - thresh
high_edge = combined['Retention Time'].values[-1] + thresh

cuts_RT = np.insert(combined[combined['Diff'] >= thresh]['Retention Time'].values, 0, low_edge)
cuts_RT = np.append(cuts_RT, high_edge)

labels = np.arange(1,len(cuts_RT))

combined['Peak'] = pd.cut(combined['Retention Time'],cuts_RT,labels=labels,include_lowest=True)

pivoted = combined.pivot(index='Result',columns='Peak',values='Area').fillna(0)