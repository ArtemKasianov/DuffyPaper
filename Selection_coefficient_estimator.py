import numpy as np
import matplotlib.pyplot as plt
import argparse
import random


def main():
    parser = argparse.ArgumentParser(description='Wright-Fisher Selection coeff Calculator')
    parser.add_argument('--p0', type=float, required=True,
                        help='Initial frequency of allele A (between 0 and 1)')
    parser.add_argument('--generations', type=int, required=True,
                        help='Generations intermediate')
    parser.add_argument('--seed', type=int, default=None,
                        help='Base random seed (optional)')
    parser.add_argument('--outputFile', type=str, default=None,
                        help='Output file')
    parser.add_argument('--freq_to_check', type=float, default=None,
                        help='Frequency to check')

    args = parser.parse_args()
    fileHandle = open(args.outputFile, "w")
    checkFreqs = float(args.freq_to_check)
    
    p_initial = float(args.p0)
    gen_inter = int(args.generations) - 1
    print(gen_inter)
    for h_index in range(0,3):
        h = h_index * 0.5
        isPrevTouched = 0
        isFound = 0
        print(h)
        for s_index in range(0,400):
            s = s_index * 0.001
            currMin = 100000
            minGen = -1
            isTouched = 0
            p = p_initial
            for gen in range(0,71):
                #print(gen)
                if gen == gen_inter:
                    p = p*0.42 + 0.58
                    print(p)
                p = p + p*(1-p)*(p*h*s + (1-p)*s*(1-h))/(1-2*p*(1-p)*s*h-(1-p)*(1-p)*s)
                if abs(p - checkFreqs) <= 0.01:
                    diff = abs(p - checkFreqs)
                    if diff < currMin:
                        currMin = diff
                        minGen = gen
                    if gen == 59:
                        isTouched = 1
            isPrevTouched = isTouched
            #print(str(s) + "\t" + str(minGen))
            if isTouched == 1:
                
                if minGen == 59:
                    fileHandle.write(str(h) + "\t" + str(s) + "\n")
                    isFound = 1
                    break

            if minGen < 59:
                if isTouched == 1:
                    fileHandle.write(str(h) + "\t" + str(s) + "\n")
                    isFound = 1
                    break
                else:
                    if isPrevTouched == 1:
                        s = s - 0.001
                        fileHandle.write(str(h) + "\t" + str(s) + "\n")
                        isFound = 1
                        break                
        
        if isFound == 0:
            fileHandle.write(str(h) + "\t" + str("not_found") + "\n")
    
    fileHandle.close()
    
    
if __name__ == '__main__':
    main()