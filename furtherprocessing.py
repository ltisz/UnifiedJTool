import datetime
import numpy as np
import math
import os,sys
import linecache

try:
    os.chdir(os.path.dirname(sys.argv[0]))  #Change working folder to folder that script is in
except:
    pass

day = "20220622"
filename = "output{}.txt".format(day)
filenameSD = "sizedistouput{}.txt".format(day)

ozoneConcs = np.genfromtxt(filename, skip_header = 1, delimiter=",", usecols = 2)

i = 0
whereChunks = []
tempChunk = []
while i<len(ozoneConcs):
    if ozoneConcs[i] == 0:
        if len(tempChunk) > 0:
            whereChunks.append(tuple(tempChunk))
            tempChunk = []
    else:
        tempChunk.append(i)
    i += 1
if len(tempChunk) > 0:
    whereChunks.append(tuple(tempChunk))

print(whereChunks)
allData = np.genfromtxt(filename, skip_header = 1, delimiter=",")

numcols = len(linecache.getline(filenameSD,1).split(","))
sizeDistData = np.genfromtxt(filenameSD, delimiter = ",", usecols=range(1,numcols))
bins = linecache.getline(filenameSD,1).split(",")[1:]

avgO3s = []
avgSizes = []
avgTotals = []
avgJs = []
avgJ5s = []
avgTO = []
avgpHOM = []
avgpHOMap = []
avgpHOMiso = []
avgR = []
avgMass = []
avgCS = []
avgCoag = []
avgWL = []
avgSizeDists = []
avgSizeDist = []

for chunk in whereChunks:
    if len(chunk) < 4:
        pass
    else:
        avgO3s.append(np.average(allData[:,2][chunk[2]:chunk[-2]]))
        avgSizes.append(np.average(allData[:,3][chunk[2]:chunk[-2]]))
        avgTotals.append(np.average(allData[:,4][chunk[2]:chunk[-2]]))
        avgJs.append(np.average(allData[:,5][chunk[2]:chunk[-2]]))
        avgJ5s.append(np.average(allData[:,6][chunk[2]:chunk[-2]]))
        avgTO.append(np.average(allData[:,9][chunk[2]:chunk[-2]]))
        avgpHOM.append(np.average(allData[:,10][chunk[2]:chunk[-2]]))
        avgpHOMap.append(np.average(allData[:,11][chunk[2]:chunk[-2]]))
        avgpHOMiso.append(np.average(allData[:,12][chunk[2]:chunk[-2]]))
        avgR.append(np.average(allData[:,16][chunk[2]:chunk[-2]]))
        avgMass.append(np.average(allData[:,17][chunk[2]:chunk[-2]]))
        avgCS.append(np.average(allData[:,18][chunk[2]:chunk[-2]]))
        avgCoag.append(np.average(allData[:,19][chunk[2]:chunk[-2]]))
        avgWL.append(np.average(allData[:,20][chunk[2]:chunk[-2]]))

for chunk in whereChunks:
    if len(chunk) < 4:
        pass
    else:
        i = 0
        while i < len(sizeDistData[0]):
            avgSizeDist.append(np.average(sizeDistData[:,i][chunk[2]:chunk[-2]]))
            i = i+1
        avgSizeDists.append(avgSizeDist)
        avgSizeDist = []

with open("averagedOutputs{}.txt".format(day), "w") as f:
    f.write("{0}-O3, {0}-Sizes, {0}-Totals, {0}-J, {0}-J5, {0}-turnover, {0}-pHOM, {0}-pHOM(AP), {0}-pHOM(iso), {0}-R, {0}-Mass, {0}-CS, {0}-Coag, {0}-WL\n".format(day))
    i=0
    while i < len(avgO3s):
        f.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
            avgO3s[i],
            avgSizes[i],
            avgTotals[i],
            avgJs[i],
            avgJ5s[i],
            avgTO[i],
            avgpHOM[i],
            avgpHOMap[i],
            avgpHOMiso[i],
            avgR[i],
            avgMass[i],
            avgCS[i],
            avgCoag[i],
            avgWL[i]))
        i=i+1

with open("averagedSizeDists{}.txt".format(day), "w") as f:
    i=0
    while i<len(bins):
        j=0
        f.write(str(bins[i]) + ",")
        while j<len(avgSizeDists):
            f.write(str(avgSizeDists[j][i]) + ",")
            j=j+1
        f.write("\n")
        i=i+1