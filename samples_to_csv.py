'''
The MIT License (MIT)
Copyright 2020 John Chen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import os
import pandas as pd

# search through all NGS read files and parse their labels
path = 'input'

all_files = []
for _, dirs, _ in os.walk(path):
    for dir in dirs:
        for _, _, files in os.walk(path+'/'+dir):
            files = [ path+'/'+dir+'/'+x for x in files ]
            all_files.extend(files)
    break

print(all_files)

data = {}
i = 0
for file in all_files:
    # parse files based on naming structure
    sep = file.split('/')[-1].split('_')
    first = sep[0]
    if first in ['AMP','CTX','MEM']: # variant selected
        label = first.lower()
        conc = sep[1]
        # further distinguish concentration, specific to this experiment
        if conc=='HIGH':
            if first=='AMP':
                conc='128'
            elif first=='CTX':
                conc='4'
        elif conc=='MED':
            if first=='AMP':
                conc='16'
            elif first=='CTX':
                conc='0.5'
            elif first=='MEM':
                conc='0.031'
        elif conc=='LOW':
            conc='2' # AMP by default
        temp = sep[2]
        group = sep[3][-1]
        rep = sep[4][-1]
        read = sep[6][-1]
    elif first=='NO': # variant non-selected
        label = 'nosel'
        conc = 'nosel'
        temp = sep[2]
        group = sep[3][-1]
        rep = sep[4][-1]
        read = sep[6][-1]
    elif first=='WT': # wt 
        label = 'wt'
        conc = 'na'
        temp = 'na'
        group = sep[1][-1]
        rep = sep[2][-1]
        read = sep[4][-1]
    data[i] = {'path':file,'label':label,'conc':conc,'temp':temp,'group':group,'rep':rep,'read':read}
    i+=1
    
d = pd.DataFrame.from_dict(data, orient='index')
d.to_csv('parsed_samples.txt',sep='\t',index=False)