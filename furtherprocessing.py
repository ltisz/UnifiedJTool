import datetime
import numpy as np
import math
import os,sys
import linecache

try:
    os.chdir(os.path.dirname(sys.argv[0]))  #Change working folder to folder that script is in
except:
    pass

day = "20220603"
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

print(whereChunks)
allData = np.genfromtxt(filename, skip_header = 1, delimiter=",")

numcols = len(linecache.getline(filenameSD,1).split(","))
sizeDistData = np.genfromtxt(filenameSD, delimiter = ",", usecols=range(1,numcols))
bins = linecache.getline(filenameSD,1).split(",")[1:]

avgO3s = []
avgSizes = []
avgTotals = []
avgJs = []
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
    avgO3s.append(np.average(allData[:,2][chunk[2]:chunk[-2]]))
    avgSizes.append(np.average(allData[:,3][chunk[2]:chunk[-2]]))
    avgTotals.append(np.average(allData[:,4][chunk[2]:chunk[-2]]))
    avgJs.append(np.average(allData[:,5][chunk[2]:chunk[-2]]))
    avgpHOM.append(np.average(allData[:,8][chunk[2]:chunk[-2]]))
    avgpHOMap.append(np.average(allData[:,9][chunk[2]:chunk[-2]]))
    avgpHOMiso.append(np.average(allData[:,10][chunk[2]:chunk[-2]]))
    avgR.append(np.average(allData[:,14][chunk[2]:chunk[-2]]))
    avgMass.append(np.average(allData[:,15][chunk[2]:chunk[-2]]))
    avgCS.append(np.average(allData[:,16][chunk[2]:chunk[-2]]))
    avgCoag.append(np.average(allData[:,17][chunk[2]:chunk[-2]]))
    avgWL.append(np.average(allData[:,18][chunk[2]:chunk[-2]]))

for chunk in whereChunks:
    i = 0
    while i < len(sizeDistData[0]):
        avgSizeDist.append(np.average(sizeDistData[:,i][chunk[2]:chunk[-2]]))
        i = i+1
    avgSizeDists.append(avgSizeDist)
    avgSizeDist = []

with open("averagedOutputs{}.txt".format(day), "w") as f:
    f.write("{}-O3, {}-Sizes, {}-Totals, {}-J, {}-pHOM, {}-pHOM(AP), {}-pHOM(iso), {}-R, {}-Mass, {}-CS, {}-Coag, WL\n")
    i=0
    while i < len(avgO3s):
        f.write("{},{},{},{},{},{},{},{},{},{},{},{}\n".format(
            avgO3s[i],
            avgSizes[i],
            avgTotals[i],
            avgJs[i],
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