import os,sys
import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd



parser = argparse.ArgumentParser(prog='FilterCentromersTelomers.py')
groupReq = parser.add_argument_group('Required parameters')
groupReq.add_argument('-fb', '--fbFileName',required=True,type=str,help="fb file name")
groupReq.add_argument('-chr', '--chrNam',required=True,type=int,help="chr name")
groupReq.add_argument('-o', '--outFileName',required=False,type=str,help="out file name")
args = parser.parse_args()

outFileHandle = open(args.outFileName,"w")

#coordinates of centromeres by chrs
centromeres ={ 
             1 : (121535434,124535434),
             2 : (92326171,95326171),
             3 : (90504854,93504854),
             4 : (49660117,52660117),
             5 : (46405641,49405641),
             6 : (58830166,61830166),
             7:  (58054331,61054331),
             8 : (43838887,46838887),
             9 : (47367679,50367679),
             10 : (39254935,42254935),
             11 : (51644205,54644205),
             12 : (34856694,37856694),
             13 : (16000000,19000000),
             14 : (16000000,19000000),
             15 : (17000000,20000000),
             16 : (35335801,38335801),
             17 : (22263006,25263006),
             18 : (15460898,18460898),
             19 : (24681782,27681782),
             20 : (26369569,29369569),
             21 : (11288129,14288129),
             22 : (13000000,16000000)
}

#coordinates of telomers by chrs
telomers ={
        1 : 249240621,
        2 : 243189373,
        3 : 198012430,
        4 : 191144276,
        5 : 180905260,
        6 : 171105067,
        7 : 159128663,
        8 : 146354022,
        9 : 141203431,
        10 : 135524747,
        11 : 134996516,
        12 : 133841895,
        13 : 115159878,
        14 : 107339540,
        15 : 102521392,
        16 : 90344753,
        17 : 81006629,
        18 : 78067248,
        19 : 59118983,
        20 : 63015520,
        21 : 48119895,
        22 : 51294566,
}





fbFileHandle = open(args.fbFileName, "r")

isFirst = 2
for line in fbFileHandle:
    line = line[:-1]
    if isFirst > 0:
        outFileHandle.write(line+"\n")
        isFirst = isFirst - 1
    else:
        #filter borders 
        filter_value = 2100000
        first = filter_value
        last = telomers[int(args.chrNam)]
        last = last - filter_value
        left_cent =centromeres[int(args.chrNam)][0] - filter_value
        right_cent =centromeres[int(args.chrNam)][1] + filter_value
    
        arrLine = line.split("\t")
        
        chrNam = int(arrLine[0])
        pos = int(arrLine[1])
        
        
        if chrNam == args.chrNam:
            
            if pos < first:
                continue
            if pos > last:
                continue
            
            if pos > left_cent and pos < right_cent:
                continue
            
       
        outFileHandle.write(line+"\n")
    
    



fbFileHandle.close()
outFileHandle.close()