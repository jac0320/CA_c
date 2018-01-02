import sys
import pandas as pd
import numpy as np
import json
from heights import busEle
from operator import add, div
from seedlib import seeds


sizes = [1, 2, 5, 10, 20, 50, 100, 200, 300, 400, 500, 1000, 2000]

seeds = seeds[0:10]

df_slr = pd.read_csv('datums_slrQuantiles.csv')
df_datum = pd.read_csv('listDatums.csv')
df_hurricane = pd.read_csv('surgeTotal.csv')

total_slr = df_slr.shape[0]
total_hurricane = 100

datum = list(df_datum['x'])
slr2datum = {round(datum[i],7):int(i+1) for i in range(len(datum))}
slrdict = df_slr.to_dict()
df_hurricane = df_hurricane.fillna(0.0)

def sample_scenario(idx, dfs, dfh, s2d, fixed_slr_idx=-1, fixed_ss_idx=-1):
    one_scenario = {}
    one_scenario['meta'] = {}
    one_scenario['SS'] = []
    one_scenario['SL'] = []
    if fixed_slr_idx >= 0:
        one_slr_idx = fixed_slr_idx
    else:
        one_slr_idx = np.random.randint(0, 99)
    for i in range(1,6):
        one_slr = dfs.loc[one_slr_idx][i]
        slr_idx = int(one_slr)
        if fixed_ss_idx >= 0:
            one_hurricane_idx = fixed_ss_idx
        else:
        	one_hurricane_idx = np.random.randint(1,100)
        one_datum_idx = s2d[round(one_slr,7)]
        while len(dfh.loc[(dfh['datum'] == one_datum_idx) & (dfh['hurricane'] == one_hurricane_idx)]) == 0:
            one_datum_idx = s2d[round(one_slr,7)]     # Resample if missing
            one_hurricane_idx = np.random.randint(1,100)
        one_row = list(dfh.loc[(dfh['datum'] == one_datum_idx) & (dfh['hurricane'] == one_hurricane_idx)].values[0])
        one_scenario['meta'][i] = {}
        one_scenario['meta'][i]['idx'] = idx
        one_scenario['meta'][i]['row_idx'] = one_row[0]
        one_scenario['meta'][i]['hurricane_idx'] = one_hurricane_idx
        one_scenario['meta'][i]['datum'] = one_datum_idx
        one_scenario['meta'][i]['slr_idx'] = slr_idx
        one_scenario['meta'][i]['t'] = i
        one_scenario['SS'].append(map(div, one_row[3:], [0.3048] * len(one_row[3:])))
        one_scenario['SL'].append(one_slr * 0.3084)

    return one_scenario

def write_to_json(S, howmany=0, filename=''):
    if howmany is 0:
        howmany = len(S)
    tasks = {}
    tasks['SS'] = {}
    tasks['SL'] = {}
    tasks['META'] = {}
    for s in range(howmany):
        tasks['SS'][str(s+1)] = S[s]['SS']
        tasks['SL'][str(s+1)] = S[s]['SL']
        tasks['META'][str(s+1)] = S[s]['meta']
    if filename == '':
        return json.dumps(tasks)
    else:
        f = open(filename, 'w')
        f.write(json.dumps(tasks))
        f.close()
    pass

# Scrips that generate inputs
fssi = -1
fslri = -1
for seed in range(len(seeds)):
    print "sampling seed is ", seeds[seed]
    np.random.seed(seeds[seed])
    scenarios= []
    for scen in range(2000):
        if len(sys.argv) == 0:
            fssi = -1
            fslri = -1
        elif sys.argv[1] == "slr":
            if sys.argv[2] == "random":
                fslri = np.random.randint(0, 99)
            else:
                fslri = int(sys.argv[2])
        elif sys.argv[1] == "ss":
            if sys.argv[2] == "random":
                fssi = np.random.randint(0, 99)
            else:
                fssi = int(sys.argv[2])
        scenarios.append(sample_scenario(scen+1, df_slr, df_hurricane, slr2datum, fixed_slr_idx=fslri, fixed_ss_idx=fssi))
    print "finished sampling..."
    for size in sizes:
        print "writing file ", str('paper2_'+str(sys.argv[1])+str(fssi)+'_s'+str(seed)+'_'+str(size)+".json")
        write_to_json(scenarios, howmany=size, filename=str('paper2_'+str(sys.argv[1])+str(fssi)+'_s'+str(seed)+'_'+str(size)+".json"))
