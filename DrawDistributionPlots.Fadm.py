import os,sys
import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
from scipy import stats

def grubbs_test_p_value(data, test_value):
    data = np.array(data)
    n = len(data)
    mean = np.mean(data)
    std = np.std(data, ddof=1)

    # Grubbs test statistic
    G = abs(test_value - mean) / std

    # Critical value of Grubbs' test (two-sided)
    t_dist = stats.t.ppf(1 - 0.05 / (2 * n), n - 2)
    G_crit = ((n - 1) / np.sqrt(n)) * np.sqrt(t_dist**2 / (n - 2 + t_dist**2))

    # P-value for the observed G (two-sided)
    t_val = ((G * np.sqrt(n - 1)) / np.sqrt(n - G**2))
    p_value = 2 * (1 - stats.t.cdf(t_val, df=n - 2))

    return G, p_value





parser = argparse.ArgumentParser(prog='Draw distributions from LAI file.')
groupReq = parser.add_argument_group('Required parameters')
groupReq.add_argument('-i', '--input_fileName',required=True,type=str,help="LAI file name")
groupReq.add_argument('-val', '--fractionToHighlight',required=False,type=float,help="fraction to highlight")
args = parser.parse_args()



lai_data = pd.read_csv(args.input_fileName,sep="\t",low_memory=False,decimal=".",header=0,index_col=0)

lai_data.apply(pd.to_numeric)



currPlot = sns.kdeplot(lai_data["Fadm"],   color='red', alpha=0.6, fill = True)
currPlot.set_title(f'Distribution of Fadm')
currPlot.set_xlabel(f'Fadm')
currPlot.set_ylabel('Frequency')

if not args.fractionToHighlight is None:
    currFraction = float(args.fractionToHighlight)
    
    quantileVal = 1 - np.count_nonzero(lai_data["Fadm"]<currFraction)/ lai_data["Fadm"].size

    print(quantileVal)
    G_stat, p_val = grubbs_test_p_value(lai_data["Fadm"], currFraction)
    print(p_val)
    print(G_stat)
    currPlot.axvline(x=currFraction, linewidth=2, color='black', ls=':')

plt.tight_layout()
plt.savefig('distributions.png', dpi=300)


