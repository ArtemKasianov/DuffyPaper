import os,sys
import argparse


parser = argparse.ArgumentParser(prog='combine ancestry info for individuals')
groupReq = parser.add_argument_group('Required parameters')
groupReq.add_argument('-fb', '--fbPrefixInputFileName',required=True,type=str,help="prefix for fb input file names")
groupReq.add_argument('-mn', '--minChrIndex',required=True,type=int,help='min_chr_index')
groupReq.add_argument('-mx', '--maxChrIndex',required=True,type=int,help='max_chr_index')
groupReq.add_argument('-t', '--treshold',required=True,type=float,help='treshold on snps prob')
groupReq.add_argument('-o', '--outFile',required=True,type=str,help='output file')
args = parser.parse_args()

samplesToPopDict = dict()
samplesSumByGroupsAndByPopNam = dict()
samplesCountByGroupsAndByPopNam = dict()
indexToPopDict = dict()
indexToGroupDict = dict()


minIndex = int(args.minChrIndex)
maxIndex = int(args.maxChrIndex)

listPopName = []
listGroupName = []
popNameAdded = dict()
groupNameAdded = dict()


outFileHandle = open(args.outFile,"w")

title_printed = 0

for chrIndex in range(minIndex,maxIndex):
    currFileHandle = open(args.fbPrefixInputFileName+".chr"+str(chrIndex)+".fb.tsv","r")
    title = ""
    for line in currFileHandle:
        line = line[:-1]
        if line[:1] == "#":
            continue
        if title == "":
            title = line
            if title_printed == 1:
                continue
            title_printed = 1
            arrTitle = title.split()
            chrNam_title = arrTitle[0]
            phys_pos_title = arrTitle[1]
            genetic_pos_title = arrTitle[2]
            genetic_marker_title = arrTitle[3]
            outFileHandle.write(chrNam_title+"\t"+phys_pos_title+"\t"+genetic_pos_title+"\t"+genetic_marker_title)
            for index in range(4,len(arrTitle)):
                currColumnName = arrTitle[index]
                arrColumnName = currColumnName.split(":::")
                sampleName = arrColumnName[0]
                groupName = arrColumnName[2]
                arrSampleName = sampleName.split('_')
                popName = arrSampleName[0]
                indexToPopDict[index] = popName
                indexToGroupDict[index] = groupName
                
                if not groupName in groupNameAdded.keys():
                    groupNameAdded[groupName] = 1
                    listGroupName.append(groupName)
                if not popName in popNameAdded.keys():
                    popNameAdded[popName] = 1
                    listPopName.append(popName)
                
                
                if popName in samplesSumByGroupsAndByPopNam.keys():
                    samplesSumByGroupsAndByPopNam[popName][groupName] = 0
                    samplesCountByGroupsAndByPopNam[popName][groupName] = 0
                else:
                    groupNameSumDict = dict()
                    groupNameCountDict = dict()
                    groupNameSumDict[groupName] = 0
                    groupNameCountDict[groupName] = 0
                    samplesSumByGroupsAndByPopNam[popName] = groupNameSumDict
                    samplesCountByGroupsAndByPopNam[popName] = groupNameCountDict
            
            
            for popNamTitle in listPopName:
                for groupNamTitle in listGroupName:
                    outFileHandle.write("\t"+popNamTitle + "_" + groupNamTitle)
            for groupNamTitle in listGroupName:
                    outFileHandle.write("\tBugaOKAKhweDEMI_" + groupNamTitle)
            for groupNamTitle in listGroupName:
                    outFileHandle.write("\tAniOKABugaOKAKhweDEMI_" + groupNamTitle)
            outFileHandle.write("\n")
                
            
            
                    
        else:
            for currPopNam in listPopName:
                for currGroupNam in listGroupName:
                    samplesSumByGroupsAndByPopNam[currPopNam][currGroupNam] = 0
                    samplesCountByGroupsAndByPopNam[currPopNam][currGroupNam] = 0
                
            
            arrLine = line.split()
            chrNam = arrLine[0]
            phys_pos = arrLine[1]
            genetic_pos = arrLine[2]
            genetic_marker = arrLine[3]
            outFileHandle.write(chrNam + "\t" + phys_pos + "\t" + genetic_pos + "\t" + genetic_marker)
            countGroups = len(listGroupName)
            
            for index in range(4,len(arrTitle),countGroups):
                
                columnValArr = []
                popNamArr = []
                groupNamArr = []
                
                tempPopNam = ""
                tempGroupNam = ""
                
                
                for i in range(0,countGroups):
                    currColumnVal = arrLine[index+i]
                    currPopName = indexToPopDict[index+i]
                    currGroupName = indexToGroupDict[index+i]
                    
                    columnValArr.append(currColumnVal)
                    popNamArr.append(currPopName)
                    groupNamArr.append(currGroupName)
                    
                    if tempPopNam == "":
                        tempPopNam = currPopName
                        tempGroupNam = currGroupName
                    else:
                        if currPopName != tempPopNam:
                            print("ErrorPopNam\n")
                        if tempGroupNam == currGroupName:
                            print("currGroupName\n")
                
                
                maxVal = float(arrLine[index])
                maxPopName = indexToPopDict[index]
                maxGroupName = indexToGroupDict[index]
                for i in range(1,countGroups):
                    currColumnVal = float(arrLine[index+i])
                    currPopName = indexToPopDict[index+i]
                    currGroupName = indexToGroupDict[index+i]
                    if currColumnVal > maxVal:
                        maxVal = float(arrLine[index+i])
                        maxPopName = indexToPopDict[index+i]
                        maxGroupName = indexToGroupDict[index+i]
                
                
                if maxVal >= args.treshold:
                    for groupNameCurr in listGroupName:
                        if groupNameCurr != maxGroupName:
                            samplesSumByGroupsAndByPopNam[maxPopName][groupNameCurr] = samplesSumByGroupsAndByPopNam[maxPopName][groupNameCurr] + 0
                            samplesCountByGroupsAndByPopNam[maxPopName][groupNameCurr] = samplesCountByGroupsAndByPopNam[maxPopName][groupNameCurr] + 1
                        else:
                            samplesSumByGroupsAndByPopNam[maxPopName][groupNameCurr] = samplesSumByGroupsAndByPopNam[maxPopName][groupNameCurr] + 1
                            samplesCountByGroupsAndByPopNam[maxPopName][groupNameCurr] = samplesCountByGroupsAndByPopNam[maxPopName][groupNameCurr] + 1
                
                
            for popNamTitle in listPopName:
                for groupNamTitle in listGroupName:
                    currSum = samplesSumByGroupsAndByPopNam[popNamTitle][groupNamTitle]
                    currCount = samplesCountByGroupsAndByPopNam[popNamTitle][groupNamTitle]
                    if currCount == 0:
                        currFraction = -1.0
                    else:
                        currFraction = currSum/currCount
                    
                    #currFractionStr = "%.2f" % (currFraction)
                    currFractionStr = str(currFraction)
                    outFileHandle.write("\t"+currFractionStr)
            
            for groupNamTitle in listGroupName:
                bothPopSum = 0
                bothPopCount = 0
                for popNamTitle in ("BugaOKA","KhweDEMI"):
                    currSum = samplesSumByGroupsAndByPopNam[popNamTitle][groupNamTitle]
                    currCount = samplesCountByGroupsAndByPopNam[popNamTitle][groupNamTitle]
                    bothPopSum = bothPopSum + currSum
                    bothPopCount = bothPopCount + currCount
                if currCount == 0:
                    currFraction = -1.0
                else:
                    currFraction = bothPopSum/bothPopCount
                currFractionStr = str(currFraction)
                outFileHandle.write("\t"+currFractionStr)
            for groupNamTitle in listGroupName:
                triplePopSum = 0
                triplePopCount = 0
                for popNamTitle in ("AniOKA","BugaOKA","KhweDEMI"):
                    currSum = samplesSumByGroupsAndByPopNam[popNamTitle][groupNamTitle]
                    currCount = samplesCountByGroupsAndByPopNam[popNamTitle][groupNamTitle]
                    triplePopSum = triplePopSum + currSum
                    triplePopCount = triplePopCount + currCount
                if currCount == 0:
                    currFraction = -1.0
                else:
                    currFraction = triplePopSum/triplePopCount
                currFractionStr = str(currFraction)
                outFileHandle.write("\t"+currFractionStr)
            
            outFileHandle.write("\n")
            
            
    currFileHandle.close()
    
outFileHandle.close()

