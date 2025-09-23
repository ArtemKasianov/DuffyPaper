import os,sys
import argparse




parser = argparse.ArgumentParser(prog='Calculate Fadm values')
groupReq = parser.add_argument_group('Required parameters')
groupReq.add_argument('-i', '--inputFileName',required=True,type=str,help="input file name from plink")
groupReq.add_argument('-aNam', '--admixedPopNam',required=True,type=str,help="admixed population name")
groupReq.add_argument('-p1Nam', '--parentPopNam1',required=True,type=str,help="first parent population name")
groupReq.add_argument('-p2Nam', '--parentPopNam2',required=True,type=str,help="second parent population name")
groupReq.add_argument('-p3Nam', '--parentPopNam3',required=True,type=str,help="third parent population name")
groupReq.add_argument('-p1Fract', '--parentPopFract1',required=True,type=float,help="first parent population fraction")
groupReq.add_argument('-p2Fract', '--parentPopFract2',required=True,type=float,help="second parent population fraction")
groupReq.add_argument('-p3Fract', '--parentPopFract3',required=True,type=float,help="third parent population fraction")
groupReq.add_argument('-o', '--outFileName',required=False,type=str,help="out file name")
args = parser.parse_args()

outFileHandle = open(args.outFileName,"w")




admixedPopNam = args.admixedPopNam
parentPopNam1 = args.parentPopNam1
parentPopNam2 = args.parentPopNam2
parentPopNam3 = args.parentPopNam3

parentPopFract1 = float(args.parentPopFract1)
parentPopFract2 = float(args.parentPopFract2)
parentPopFract3 = float(args.parentPopFract3)



admixPopFreq = dict()
parent1PopFreq = dict()
parent2PopFreq = dict()
parent3PopFreq = dict()


listVariant = list()
allele1ByVarIdent = dict()
allele2ByVarIdent = dict()

# read plink file with frequencies
inputFileHandle = open(args.inputFileName,"r")
title = inputFileHandle.readline()
for line in inputFileHandle:
    line = line[:-1]
    arrLine = line.split(" ")
    
    index = 0
    
    for currVal in arrLine:        
        if currVal == "":
            continue
        if index == 0:
            chrNam = currVal
            index = index + 1
        else:
            if index == 1:
                varIdent = currVal
                index = index + 1
            else:
                if index == 2:
                    clusterIdent = currVal
                    index = index + 1
                else:
                    if index == 3:
                        allele1Val = currVal
                        index = index + 1
                    else:
                        if index == 4:
                            allele2Val = currVal
                            index = index + 1
                        else:
                            if index == 5:
                                allele1Freq = float(currVal)
                                allele2Freq = 1-float(currVal)
                                index = index + 1
    
    
    
    allele1VariantPos = varIdent+"\t"+allele1Val
    allele2VariantPos = varIdent+"\t"+allele2Val
    
    
    if not varIdent in allele1ByVarIdent.keys():
        listVariant.append(varIdent)
        allele1ByVarIdent[varIdent] = allele1Val
        allele2ByVarIdent[varIdent] = allele2Val
        
    if clusterIdent == parentPopNam1:
        parent1PopFreq[allele1VariantPos] = allele1Freq
        parent1PopFreq[allele2VariantPos] = allele2Freq
    if clusterIdent == parentPopNam2:
        parent2PopFreq[allele1VariantPos] = allele1Freq
        parent2PopFreq[allele2VariantPos] = allele2Freq
        
    if clusterIdent == parentPopNam3:
        parent3PopFreq[allele1VariantPos] = allele1Freq
        parent3PopFreq[allele2VariantPos] = allele2Freq
    
    if clusterIdent == admixedPopNam:
        admixPopFreq[allele1VariantPos] = allele1Freq
        admixPopFreq[allele2VariantPos] = allele2Freq
    
inputFileHandle.close()

outFileHandle.write("Chromosome\tPosition\tFadm\tExpFreqAllele1\tObsFreqAllele1\tExpFreqAllele2\tObsFreqAllele2\tcurrParent1PopFreqAllele1\tcurrParent2PopFreqAllele1\tcurrParent3PopFreqAllele1\tcurrParent1PopFreqAllele2\tcurrParent2PopFreqAllele2\tcurrParent3PopFreqAllele2\n")
for currVarIdent in listVariant:
    currAllele1 = allele1ByVarIdent[currVarIdent]
    currAllele2 = allele2ByVarIdent[currVarIdent]
    
    allele1VariantPos = currVarIdent+"\t"+currAllele1
    allele2VariantPos = currVarIdent+"\t"+currAllele2
    
    if allele1VariantPos in parent1PopFreq.keys():
        if allele1VariantPos in parent2PopFreq.keys():
            if allele1VariantPos in parent3PopFreq.keys():
                if allele1VariantPos in admixPopFreq.keys():
                    arrCurrVarIdent = currVarIdent.split('_')
                    chrNam = arrCurrVarIdent[0] # chromosome name
                    chrPos = arrCurrVarIdent[1] # position inside chromosome
                    
                    currParent1PopFreqAllele1 = parent1PopFreq[allele1VariantPos] # frequency of first allele in first parental population for 3 way admixture model
                    currParent2PopFreqAllele1 = parent2PopFreq[allele1VariantPos] # frequency of first allele in second parental population for 3 way admixture model
                    currParent3PopFreqAllele1 = parent3PopFreq[allele1VariantPos] # frequency of first allele in third parental population for 3 way admixture model
                    currAdmixPopFreqAllele1 = admixPopFreq[allele1VariantPos] # frequency of first allele in admixed population for 3 way admixture model
                    
                    currParent1PopFreqAllele2 = parent1PopFreq[allele2VariantPos] # frequency of second allele in first parental population for 3 way admixture model
                    currParent2PopFreqAllele2 = parent2PopFreq[allele2VariantPos] # frequency of second allele in second parental population for 3 way admixture model
                    currParent3PopFreqAllele2 = parent3PopFreq[allele2VariantPos] # frequency of second allele in third parental population for 3 way admixture model
                    currAdmixPopFreqAllele2 = admixPopFreq[allele2VariantPos] # frequency of second allele in admixed population for 3 way admixture model
                    
                    yAllele1 = parentPopFract1*currParent1PopFreqAllele1 + parentPopFract2*currParent2PopFreqAllele1 + parentPopFract3*currParent3PopFreqAllele1 # expected first allele frequency in admixed population
                    yAllele2 = parentPopFract1*currParent1PopFreqAllele2 + parentPopFract2*currParent2PopFreqAllele2 + parentPopFract3*currParent3PopFreqAllele2 # expected second allele frequency in admixed population
                    
                    #Fadm values calculation
                    FadmNom = (currAdmixPopFreqAllele1 - yAllele1)*(currAdmixPopFreqAllele1 - yAllele1)+(currAdmixPopFreqAllele2 - yAllele2)*(currAdmixPopFreqAllele2 - yAllele2)
                    FadmDeNom = 2*(1-(yAllele1*yAllele1+yAllele2*yAllele2))
                    Fadm = 0.0
                    if FadmDeNom == 0.0:
                        Fadm = 0.0
                    else:
                        Fadm = FadmNom/FadmDeNom
                    
                    if Fadm > 1000:
                        Fadm = 0.0
                    
                    outFileHandle.write(chrNam+"\t"+chrPos+"\t" + str(Fadm)+"\t" + str(yAllele1)+"\t" + str(currAdmixPopFreqAllele1)+"\t" + str(yAllele2)+"\t" + str(currAdmixPopFreqAllele2)+"\t" + str(currParent1PopFreqAllele1)+"\t" + str(currParent2PopFreqAllele1)+"\t" + str(currParent3PopFreqAllele1)+"\t" + str(currParent1PopFreqAllele2)+"\t" + str(currParent2PopFreqAllele2)+"\t" + str(currParent3PopFreqAllele2)+"\n")
                
            
outFileHandle.close()

