import os,sys
import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd



parser = argparse.ArgumentParser(prog='CalculateSDStatsAndDrawPlotsForChrs.percentile.py')
groupReq = parser.add_argument_group('Required parameters')
groupReq.add_argument('-i', '--laiFileName',required=True,type=str,help="prefix for LAI file name")
groupReq.add_argument('-chr', '--chromosome',required=True,type=int,help="chromosome index")
groupReq.add_argument('-sn', '--sampleName',required=True,type=str,help="sample name")
groupReq.add_argument('-fai', '--fastaIndexFileName',required=True,type=str,help="fasta index file name")
groupReq.add_argument('-o', '--outFile',required=True,type=str,help='output file')
args = parser.parse_args()


chrLensDict = dict()



lai_data = pd.read_csv(args.laiFileName,sep="\t",low_memory=False,decimal=".",header=0)

lai_data.apply(pd.to_numeric)

sampleName = args.sampleName

arrSampleName = sampleName.split('_')
ancestryNam = arrSampleName[1]

# get ancestry name for y axis
if ancestryNam == "Bantu":
    ancestryNam = "Bantu"
else:
    if ancestryNam == "Khoisan":
        ancestryNam = "South Africa"
    else:
        ancestryNam = "East Africa"
    


pos_lai_data = lai_data['physical_position']
curr_lai_data = lai_data[args.sampleName]
# Calculate percentiles and sd
p99 = curr_lai_data.quantile(0.99)
p01 = curr_lai_data.quantile(0.01)

averageColumn = curr_lai_data.mean()
stdColumn = curr_lai_data.std()
std3MinusColumn = averageColumn - 3*curr_lai_data.std()
std3PlusColumn = averageColumn + 3*curr_lai_data.std()

print(averageColumn)

lai_data = lai_data[lai_data.chromosome==int(args.chromosome)]
pos_lai_data = lai_data['physical_position']
curr_lai_data = lai_data[args.sampleName]




fastaIndexFileHandle = open(args.fastaIndexFileName,"r")

for line in fastaIndexFileHandle:
    line = line[:-1]    
    arrLine = line.split()
    chrNam = arrLine[0]
    chrSize = arrLine[1]
    chrLensDict[chrNam] = int(chrSize)
fastaIndexFileHandle.close()

maxChrSize = chrLensDict["chr" + str(args.chromosome)]

sns.set_theme(rc={'figure.figsize':(10,4)})
sns.set_style("white")
fig, ax = plt.subplots()
#Start plotting
currPlot = sns.lineplot(x=pos_lai_data, y=curr_lai_data, color='red', alpha=0.6, ax=ax)
for spine in ax.spines.values():
    spine.set_visible(False)



currPlot.set_title(f'')
currPlot.set_ylabel(ancestryNam + ' related ancestry')
currPlot.set_xlabel('Position in chromosome 1, 100 Mbp')
#currPlot.set_ylim(0.3, 0.8) #you can use it to limit y axis
currPlot.axvline(x=159174123, linewidth=1.0, color='black', ls=':')
currPlot.axhline(y=p99, linewidth=1.0, color='grey', ls=':') #99 percentile line
currPlot.axhline(y=p01, linewidth=1.0, color='grey', ls=':') #1 percentile line 
currPlot.axhline(y=averageColumn, linewidth=1.0, color='black') # average line



plt.tight_layout()
plt.savefig(args.outFile, dpi=300)

