## -*- coding: utf-8 -*-
## UNIFIEDJTOOL.py --  Takes input data from O3 monitor, GC, SMPS and PSM
##					   Calculates J, size distributions, average sizes, HOMs, mass, condensation sink
##					   Saves processed data to output-HHMM.txt, plots combined size distribution from PSM/SMPS
##					 
## Categories:
##		Aerosol nucleation
## Author:
##		Lee Tiszenkel, UAH
##
##IMPORTANT - all input files need to be named and formatted by a certain convention:
##Need the following files with appropriate names:
##
##SMPS:
##   Export dN/dLogDp data by ROW, COMMA as:
##      MMDDYY-dN.txt
##   Export concentration data by ROW, COMMA as:
##      MMDDYY-conc.txt
##
##PSM:
##   Same as output files from inversion.py:
##      MMDDYY-PSMconc.txt
##      MMDDYY-PSMconcdN.txt
##
##O3:
##   Name O3 concentration file exported from 49i as MMDDYY-O3.dat
##      You may need to trim this file as needed, as the 49i exports in bulk.
##
##Organics:
##   Create a file MMDDYY-Organics.txt with each line formatted as follows:
##      MM/DD/YY HH:MM [MT peak area] [Isoprene peak area]
##      See example file for exact format
##   If there is GC data measured for the end of the tube, creat another called MMDDYY-OrganicsEnd.txt
##      with the same formatting
##
##   Set "day" variable to MMDDYY, same as files
##
##   Set "PSM" and "GC" to TRUE or FALSE to match experiment. Program will work with SMPS only and static monoterpene
##   IMPORTANT: If GC is set to FALSE, MT and ISO must be set MANUALLY.

######################### USER SETTINGS HERE -- CHANGE SETTINGS BASED ON YOUR DATA ###############################

#Set the DAY for your data

day = "05202022"

###### WHAT DATA DO YOU HAVE ######
#This program can work in several ways
#SMPS only (No PSM data)
#Static monoterpene (No GC data, MT constant and only calculated)
#
#If PSM data is included, change to PSM = True
#If GC data is included (see format above) set GC = True
#If GC data is present for both beginning and end set GCyield = True

PSM = True
GC = False
GCyield = False

#Set MT, ISO to correct values ONLY IF GC not used for experiment
if GC == False:
    MT = 243.0
    ISO = 0.0

#RESIDENCE TIMES - Based on flow Excel sheet. Input in SECONDS
#Inlet
resTime1 = 1.45

#FT-1
resTime2 = 4.12

#FT-2
resTime3 = 114.49

#TEMPERATURE
T = 295.0

######################### DO NOT CHANGE ANYTHING AFTER THIS LINE UNLESS YOU KNOW WHAT YOU'RE DOING ######################

#Imported libraries
import datetime
import numpy as np
import math
import linecache
from functools import reduce
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LogNorm
import os,sys

os.chdir(os.path.dirname(sys.argv[0]))  #Change working folder to folder that script is in

def Reverse(lst): 
    return [ele for ele in reversed(lst)] 

def isFloat(num):
    try:
        float(num)
        return True
    except:
        return False

def throwSomeDs(dee):    #input diameter in nm
    deeUm = dee*1e-3    #convert to um
    deeCm = dee*1e-7    #convert to cm
    deeM = dee*1e-9     #convert to m
    Cc = 1. + (((2.*lamb)/(deeUm))*(1.257+(0.4*(exp**((-1.1*deeUm)/(2*lamb))))))       #Slip correction
    bigDee = ((k*T*Cc)*(100**2)*1000)/(3*pi*mew*(deeCm))
    mass = ((((4/3)*pi*((deeM/2)**3))*rho)*1000)
    Ca = (((8*k*T)/(pi*mass))**0.5)*100
    MFP = (8*bigDee)/(pi*Ca)
    gfirstTerm = ((2**0.5)/(3*deeM*MFP))
    gsecondTerm = ((deeCm + MFP)**3) - (((deeCm**2)+(MFP**2))**(3/2))
    g = (gfirstTerm*gsecondTerm) - deeCm
    return bigDee, Ca, g

def KoneTwo(D1, D2, Dp1, Dp2, g1, g2, Ca1, Ca2):
    K12a = (2*pi*(D1+D2)*(Dp1+Dp2))
    K12baTop = (Dp1 + Dp2)
    K12baBottom = (Dp1 + Dp2 + (2*(((g1*g1)+(g2*g2))**0.5)))
    K12bbTop = (8*(D1+D2))
    K12bbBottom = (((Ca1*Ca1)+(Ca2*Ca2))**0.5)*(Dp1+Dp2)
    K12b = 1/((K12baTop/K12baBottom) + (K12bbTop/K12bbBottom))
    K12 = K12a*K12b
    return K12
    
#COAG SINK FUNCTION
def CoagS(di, dj, nj):
    D1, Ca1, g1 = throwSomeDs(di)
    D2, Ca2, g2 = throwSomeDs(dj)
    Kij = KoneTwo(D1, D2, di*1e-7, dj*1e-7, g1, g2, Ca1, Ca2)
    return (Kij * nj)

def WallLoss(ni, kwall):
    WL = []
    for k in kwall:
        WL.append(ni * k)
    totalLoss = reduce((lambda x, y: x * y), WL)
    return abs(totalLoss)
    
#def GetData(filenames)
#Constants
k = (1.381*(10.**-23.))          #Boltzmann constant
pi = 3.1415927                   #Pi
exp = 2.718281828459045          #e
lamb = 0.0651                    #Mean free path in air (micrometers)
mu = (1.80*(10**-5))             #Viscosity of air (Pa * s)
mew = (1.72e-4)
#Wall Loss Calc Values
D = 0.0716                       #Diffusion coefficient of alpha-pinene
alpha = 0.1                      #Accommodation coefficient
tubeRads = np.array([3.2,4.85,6.35])
resTimes = np.array([resTime1, resTime2, resTime3])
wallLossKs = 0.1*(3.65*(D/(tubeRads**2)))
dilution = 1                     #dilution factor for ft-2

#HOM Calcs
y1 = 0.029                       #AP-O3 yield (CLOUD, Trostl 2016)
k1 = 8.06e-17                    #AP-O3 reaction rate (CLOUD, Trostl 2016)
y2 = 0.01                        #Isoprene-O3 yield (Jokinen 2015, Bianchi 2019)
k2 = 9.6e-18                     #Isoprene-O3 reaction rate (Karl 2004)
ppbconv = 2.445e10                 #ppb to molecules cm^-3

#Variables
rho = 1.5                        #Density (g cm^-3)
h = 1.                           #Correction factor in gamma calculation, taken to be 1 with particle measurements
restime = sum(resTimes)          #residence time

filenameSMPS = "{}-conc.txt".format(day)
filenameSMPSdN = "{}-dN.txt".format(day)
filenameO3 = "{}-O3.dat".format(day)

print("Collecting data...")
#GETTING SMPS BINS
with open(filenameSMPS, 'r') as f:
    lines = f.readlines()
binLine = lines[17].split(',')                                      #Line 17 always contains size bins
resultBins = filter(lambda x: isFloat(x[1]), enumerate(binLine))    #Filter out all non-size bin columns
rawBins = list(resultBins)              #Rawbins is list of tuples (list index, bin)
indexes = []
SMPSbins = []
for i in rawBins:
    indexes.append(i[0])                #Create list of bin column indexes
    SMPSbins.append(i[1])               #Create list of SMPS bins
indexesTup = tuple(indexes)             #Make indexes in to tuple to retrieve data later

if GC == True:
    filenameOrg = "{}-organics.txt".format(day)

#IF WE HAVE PSM DATA:
if PSM == True:
    print("Combining PSM and SMPS size bins...")
    filenamePSM = "{}-PSMconc.txt".format(day)
    filenamePSMdN = "{}-PSMconcDN.txt".format(day)
    PSMbins = linecache.getline(filenamePSMdN,1).split(",")     #Get PSM bins from frist line of PSM dN
    PSMbins = Reverse(PSMbins)                                  #PSM dN file is in reverse size order, flip them
    PSMnewbins = []
    for h in PSMbins:
        PSMnewbins.append(round(float(h),2))                    #Get 2 decimal places
    PSMdNbins = []
    i = 0
    while i < len(PSMnewbins)-1:                                #Get center of each size bin
        PSMdNbins.append((PSMnewbins[i]+PSMnewbins[i+1])/2)
        i += 1
    bins = tuple(PSMnewbins + SMPSbins)                         #Combine PSM and SMPS bins
    binsdN = tuple(PSMdNbins + SMPSbins)
else:
    bins = tuple(SMPSbins)
    binsdN = tuple(SMPSbins)
    
#Getting the line from the SMPS output with the size bins
floatBins = list(float(x) for x in bins)                        #Make size bins floats, not strings
floatBinsdN = list(float(x) for x in binsdN)

print("Bins:")
print(floatBins)
print("Done.")

#SMPS data
#Retrieving times for scans as well as size distributions (both raw concentration and dN/dLogDp)
timelistSMPS = np.genfromtxt(filenameSMPS, skip_header=18, delimiter=',', usecols=(1,2), dtype="str")
sizedistSMPSraw = np.genfromtxt(filenameSMPS, skip_header=18, delimiter=',', usecols=indexesTup)
sizedistSMPSdNraw = np.genfromtxt(filenameSMPSdN, skip_header=18, delimiter=',', usecols=indexesTup)

#Ozone data
O3list = np.genfromtxt(filenameO3, skip_header=6,usecols=(0,1,3),dtype="str")

#PSM data
#Retrieving times for scans as well as size distributions (both raw concentration and dN/dLogDp)
if PSM == True:
    timelistPSM = np.genfromtxt(filenamePSM, delimiter=',', usecols=(0), dtype='str')
    sizedistPSMraw = np.genfromtxt(filenamePSM, delimiter=',', usecols=tuple(range(1,8)))
    sizedistPSMraw = np.flip(sizedistPSMraw, 1)
    sizedistPSMdNraw = np.genfromtxt(filenamePSMdN, skip_header=1, delimiter=',', usecols=tuple(range(1,7)))

#Organics data
if GC == True:
    timelistOrganics = np.genfromtxt(filenameOrg, usecols=(0,1), dtype="str")
    orgConcs = np.genfromtxt(filenameOrg, usecols=(2,3))
if GCyield == True:
    timelistOrganics = np.genfromtxt(filenameOrgEnd, usecols=(0,1), dtype="str")
    orgConcsEnd = np.genfromtxt(filenameOrgEnd, usecols=(2,3))
    
#Initializing lists for full size distributions
#Full size distributions are in the format [[Datetime, [Size distribution]]]
#For example, fullSizeDistSMPS[2][0] will give you the time of the second scan from the SMPS;
#             fullSizeDistSMPS[2][1][3] will give you the concentration in the 4th size bin of the second scan from the SMPS
fullSizeDistSMPS = []
fullSizeDistSMPSdN = []
datetimelistSMPS = []
timelist = []
i = 0
while i<len(timelistSMPS):
    lineTime = datetime.datetime.strptime("{} {}".format(timelistSMPS[i][0],timelistSMPS[i][1]), "%m/%d/%y %H:%M:%S")
    fullSizeDistSMPS.append([lineTime, sizedistSMPSraw[i]])
    fullSizeDistSMPSdN.append([lineTime, sizedistSMPSdNraw[i]])
    datetimelistSMPS.append(lineTime)
    if PSM == False:
        timelist.append(lineTime)       #timelist is a separate list of times used for graphing later.
    i = i+1                             #When PSM is used, it's built when SMPS and PSM scans are matched. If no PSM, it needs to be built here.

#Retrieving size distributions from PSM data
if PSM == True:
    fullSizeDistPSM = []
    fullSizeDistPSMdN = []
    i = 0
    while i<len(timelistPSM):
        lineTime = datetime.datetime.strptime(timelistPSM[i], "%m/%d/%Y %H:%M:%S")
        fullSizeDistPSM.append([lineTime, sizedistPSMraw[i]])
        fullSizeDistPSMdN.append([lineTime, sizedistPSMdNraw[i]])
        i = i+1

#Retrieving organics data
if GC == True:
    datetimelistOrg = []
    for i in timelistOrganics:
        datetimelistOrg.append(datetime.datetime.strptime("{} {}".format(i[0],i[1]), "%m/%d/%y %H:%M"))

#Getting the lower and upper time limits of the ozone data. This will be used later to match up the ozone concentrations with the time of each scan.
O3listMin = datetime.datetime.strptime("{} {}".format(O3list[0][0], O3list[0][1]), "%H:%M %m-%d-%y")
O3listMax = datetime.datetime.strptime("{} {}".format(O3list[len(O3list)-1][0], O3list[len(O3list)-1][1]), "%H:%M %m-%d-%y")

#Initializing lists for full size distributions (combined PSM and SMPS)
fullSizeDist = []
fullSizeDistdN = []
organicsList = []
organicsListEnd = []
O3ListTrim = []
rList = []

#Combining SMPS and PSM size distributions based on PSM time
if PSM == True:
    i = 0
    while i<len(fullSizeDistPSM):                   #Iterate through each PSM line
        PSMtime = fullSizeDistPSM[i][0]
        j=0
        while j<len(fullSizeDistSMPS):              #For each PSM line, iterate through each SMPS line
            SMPStime = fullSizeDistSMPS[j][0]
            threshold = datetime.timedelta(seconds=67)
            if SMPStime>(PSMtime-threshold) and SMPStime<(PSMtime+threshold):                                               #If SMPS time is within 67 seconds of PSM time
                fullSizeDist.append([PSMtime, Reverse(fullSizeDistPSM[i][1].tolist()) + fullSizeDistSMPS[j][1].tolist()])   #Combine PSM and SMPS size distribution
                fullSizeDistdN.append([PSMtime,Reverse(fullSizeDistPSMdN[i][1].tolist()) + fullSizeDistSMPSdN[j][1].tolist()])
                timelist.append(PSMtime)
            j = j+1
        i = i+1
else:
    fullSizeDist = fullSizeDistSMPS             #If no PSM, full size distribution is just the SMPS full size distribution.
    fullSizeDistdN = fullSizeDistSMPSdN

print("Done.")

#### ORGANICS AND O3 DATA ORGANIZING ####
print("Calculating organics, O3...")
apConcs = 0
isoConcs = 0
apConcsEnd = 0
isoConcsEnd = 0
p = 0
q = 0
for i in fullSizeDist:
    if GC == True:                                              #If we have GC data, scan through the organics data to match GC scans with size distributions
        if i[0] > datetimelistOrg[p]:                           #If the datetime object for a scan is greater than the time of an organics scan, retrieve it
            if orgConcs[p][0] == 0:                             #If a-p = 0, set to zero
                apConcs = 0
            else:                                               #Else, put the GC peak area in to the calibration curve equation
                #apConcs = (float(orgConcs[p][0])*37.446)+1.19755
                apConcs = orgConcs[p][0]
            if orgConcs[p][1] == 0:                             #Same for isoprene
                isoConcs = 0
            else:
                #isoConcs = (float(orgConcs[p][1])*34.646)-2.0489
                isoConcs = orgConcs[p][1]
            concs = [apConcs, isoConcs]                         #Make 2-element list of concentrations
            if p == len(datetimelistOrg)-1:                     #p is a check to see when to stop looking at the organics data (last scan)
                p = p
            else:
                p+=1
        else:
            #rList.append(0)
            concs = [apConcs, isoConcs]                         #If the datetime object for a scan is LESS than the time of an organics scan, set concs to same concs as previous.
        if GCyield == True:
            if orgConcsEnd[q][0] == 0:
                apConcsEnd = 0
            else:
                apConcsEnd = (float(orgConcsEnd[q][0])*37.446)+1.19755
            if orgConcsEnd[q][0] == 0:
                isoConcsEnd = 0
            else:
                isoConcsEnd = (float(orgConcsEnd[q][1])*34.646)-2.0489
            concsEnd = [apConcsEnd, isoConcsEnd]
            if q == len(datetimelistOrgEnd)-1:
                q=q
            else:
                p+=1
            organicsListEnd.append([i[0],concsEnd)
        organicsList.append([i[0],concs])
    else:
        organicsList.append([i[0],[MT, ISO]])                   #If no GC, set organics to static values set at top of program
        organicsListEnd.append([i[0], [0,0]])
    if (O3listMin > i[0]) or (i[0] > O3listMax):                #If the earliest ozone scan is later than the current size distribution OR if the current size dist is later than the last O3 scan
        O3ListTrim.append([i[0],0])                             #Set ozone to zero (Basically set ozone to zero if there's no ozone data for that scan)
    else:
        for m in O3list:
            if "{} {}".format(m[0],m[1]) == datetime.datetime.strftime(i[0], "%H:%M %m-%d-%y"):     #Matching ozone times to scan times; ozone monitor returns 1 reading per minute, so there should
                if float(m[2]) < 0:                                                                 #always be a match. Negative ozone values are considered zero.
                    O3ListTrim.append([i[0],0])
                else:
                    O3ListTrim.append([i[0],m[2]])
for x in organicsList:                              #Assembling R list, Isoprene Carbon divided by a-p carbon.
    if (x[1][1] == 0) or (x[1][0] == 0):            #If either are zero, R is zero (avoids divide by zero errors)
        rList.append(0)
    else:
        rList.append((x[1][1]/2)/x[1][0])

#Initializing lists and loop parameters for next calculations
b=0
Jlist = []          #Nucleation rate
CSlistty = []       #Condensation sink
Coaglisty = []      #Coagulation sink
WLListy = []        #Wall loss
avgSizeList = []    #Average size
HOMlist = []        #Total HOM
HOMlistAP = []      #AP HOM
HOMlistISO = []     #Isoprene HOM
pHOMlist = []       #Total pHOM
pHOMlistAP = []     #AP pHOM
pHOMlistISO = []    #Isoprene pHOM
Jlistprint = []     #Nucleation rate as strings for output
totallist = []      #Total concentration
totalMassList = []  #Total mass
print("Done.")

#### [HOM], LOSS, J CALCULATIONS ####
print("Calculating [HOM], CS, WL, J...")
while b<(len(fullSizeDist)-1):  #Iterate through every size distribution
    WLList = []                 #Wall loss for this individual size distribution
    CoagList = []               #Coagulation sink for this individual size distribution
    i = len(bins)-1
    j = i-1
    total = sum(fullSizeDist[b][1])     #Total is simply the sum of all size bins
    totallist.append(total)             #Append total
    if total == 0:                      #Accounting for divide by zero errors
        theAvg = 0
        avgSizeList.append(theAvg)
    else:
        avgList = []
        r = 0
        while r<(len(bins)-1):
            proportion = float(fullSizeDist[b][1][r])/float(total)
            avg = proportion*float(bins[r])
            avgList.append(avg)
            r+=1
        theAvg = sum(np.array(avgList))
        avgSizeList.append(theAvg)
    #while i >= 0:
    while j >= 0:
        CoagSink = CoagS(float(bins[0]), float(bins[j]), fullSizeDist[b][1][j])
        CoagList.append(CoagSink)
        j-=1
    #    i-=1
    #    j=i-1
    avgVol = ((4/3)*pi*((theAvg/2)**3))*(1e-21)         #average volume (cm^-3)
    avgMass = avgVol*rho
    totalMass = (avgMass*total)*(1e6)                   #total mass (ug cm^-3)
    totalMassList.append(totalMass)
    CoagSinky = sum(CoagList)
    i = len(bins)-1
    CSlist = []
    while i >= 0:
        Kni = (2.*lamb)/(float(bins[i])/1000.)
        betaM = (1+Kni)/(1+(0.337*Kni)+((4/3)*Kni)+((4/3)*Kni*Kni))
        CSlist.append(betaM*(float(bins[i])*(1e-7))*fullSizeDist[b][1][i])
        WLList.append(WallLoss((fullSizeDist[b][1][i-1]-fullSizeDist[b][1][i]), wallLossKs))
        i-=1
    CS = 2*pi*D*(sum(CSlist))
    wallLosses = sum(WLList)
    pHOMconcAP = ((y1*k1*(float(O3ListTrim[b][1])*ppbconv))*((float(organicsList[b][1][0])*ppbconv)))
    pHOMconcISO = ((y2*k2*(float(O3ListTrim[b][1])*ppbconv))*((float(organicsList[b][1][1])*ppbconv)))
    pHOMlistAP.append(pHOMconcAP)
    pHOMlistISO.append(pHOMconcISO)
    pHOMlist.append(pHOMconcAP + pHOMconcISO)
    if CS == 0:
        HOMlist.append(0)
    else:
        HOMlistAP.append(pHOMconcAP/CS)
        HOMlistISO.append(pHOMconcISO/CS)
        HOMlist.append((pHOMconcAP/CS) + (pHOMconcISO/CS))
    i+=1
    sizeDistCorr = fullSizeDist[b][1]
    CSlistty.append(CS)
    Coaglisty.append(CoagSinky)
    WLListy.append(wallLosses)

    Jcorr = (sum(sizeDistCorr[3:])/restime) + CoagSinky + wallLosses
    Jlist.append([fullSizeDist[b][0], Jcorr])
    Jlistprint.append(Jcorr)
    b+=1

#### CREATING OUTPUT.TXT ####
print("Assembling output.txt...")

timeStrs = []
for x in timelist:
    timeStrs.append(datetime.datetime.strftime(x, "%m/%d/%y, %H:%M"))
print(timeStrs)
o3Only = []
for x in O3ListTrim:
    o3Only.append(x[1])
apList = []
isoList = []
for x in organicsList:
    apList.append(x[1][0])
    isoList.append(x[1][1])
zippy = zip(timeStrs, o3Only, avgSizeList, totallist, Jlistprint, apList, isoList, pHOMlist, pHOMlistAP, pHOMlistISO, HOMlist, HOMlistAP, HOMlistISO, rList, totalMassList, CSlistty, Coaglisty, WLListy)
mylist = []
for i in zippy:
    mylist.append(i)

outputLabel = day + datetime.datetime.strftime(datetime.datetime.now(), "%H%M")
with open('output{}.txt'.format(outputLabel), 'w') as fp:
    fp.write('{0}_Date, {0}_Time, {0}_O3, {0}_AvgSize, {0}_TotalConc, {0}_J, {0}_AP, {0}_ISO, {0}_pHOM, {0}_pHOM(AP), {0}_pHOM(ISO), {0}_HOM, {0}_HOM(AP), {0}_HOM(ISO), {0}_R, {0}_Mass, {0}_CS, {0}_CoagS, {0}_WL\n'.format(day))
    for x in mylist:
        fp.write(str(x).replace("(","").replace(")","").replace("'","") + "\n")
print("Done.")

#### GRAPHING ####
print("Graphing...")
fullSizeDist3 = []
timelistContour = []
for a in fullSizeDistdN:
    timelistContour.append(a[0])
    fullSizeDist3.append(a[1])

fig = plt.figure(figsize=(30,6))
ax = fig.add_subplot(1,1,1)
ax.clear()

print(len(bins))
print(np.shape(np.array(fullSizeDist3).T))
print(len(timelist))

bananaPlot = ax.contourf(timelistContour,floatBins,np.array(fullSizeDist3).T+0.1, np.arange(0.1, 100000, 100), cmap='jet', extend='max')
ax.set_yscale('log')
cb = fig.colorbar(bananaPlot, ax = ax, extendrect=True)
cb.set_label('dN/dlogDp')
ax.set_xlabel('Time (LT)')
ax.set_ylabel('Dp (nm)')

myFmt = mdates.DateFormatter("%m/%d/%y") #Format date
myFmt2 = mdates.DateFormatter("%H:%M:%S") #Format time
myFmt3 = mdates.DateFormatter("%m/%d") #SO2 Format date


ax.xaxis.set_major_formatter(myFmt2)
ax.xaxis.set_major_locator(mdates.SecondLocator(interval=int(((max(datetimelistSMPS)-min(datetimelistSMPS)).total_seconds())/5))) #6 marks on x axis
ax.xaxis.set_minor_formatter(myFmt)
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=1))
ax.xaxis.set_tick_params(which='minor', pad=15) #Keeps major/minor axis from overlapping

ax2 = ax.twinx()
ax2.plot(timelistContour, [float(x) for x in o3Only], color='black')
ax2.set_ylabel('[O3] (ppb)')
ax2.set_xlabel('Time (LT)')
plt.savefig("sizedist{}_{}.png".format(day, outputLabel), bbox_inches='tight', pad_inches=0) #Saves