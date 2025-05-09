# Receive date time data with miliseconds in format %Y.%m.%d-%H.%M.%S.%f

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime

import re
import csv
import os
import fnmatch
from scipy import optimize
from scipy import interpolate
from scipy import signal
from scipy.interpolate import interp1d

import random
import pwlf

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       #raise Exception('lines do not intersect')
       
       return 0, 0

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y


def resample(time, xdata, ydata, newFrameRate):
    
    #re sampling time 
    duration = time[-1] - time[0]
    nbSamples = (duration * newFrameRate) + 1
    newTime = np.linspace(time[0], time[-1], len(time))
    
    x = xdata[:]-xdata[0]
    y = ydata[:]-ydata[0]
    
    #resample x data accordingly
    xInter = interp1d(time, x, kind='linear')
    new_xdata = xInter(newTime)
    
    #resample y data accordingly
    yInter = interp1d(time, y, kind='linear')
    new_ydata = yInter(newTime)
    
    #butterworth lowpass filter
    Wn = 2.0*fc/newFrameRate
    sos = signal.butter(1, Wn, 'low', output='sos')
    
    xdataFilter = signal.sosfilt(sos, new_xdata) +xdata[0]
    ydataFilter = signal.sosfilt(sos, new_ydata) +ydata[0]
    
    return (newTime, xdataFilter, ydataFilter)


def resample3D_scale(time, xdata, ydata, zdata, newFrameRate):
    
    #re sampling time 
    
    duration = time[-1] - time[0]
    nbSamples = (duration * newFrameRate) + 1
    #newTime = np.linspace(time[0], time[-1], len(time))
    newTime = np.linspace(time[0], time[-1], int(nbSamples))
    
    x = (xdata[:]-xdata[0])/max(xdata)
    y = (ydata[:]-ydata[0])/max(ydata)
    z = zdata[:]-zdata[0]
    
    #resample x data accordingly
    xInter = interp1d(time, x, kind='linear', assume_sorted = True, fill_value="extrapolate")
    new_xdata = xInter(newTime)
    #new_xdata = xInter(time)
    
    #resample y data accordingly
    yInter = interp1d(time, y, kind='linear', assume_sorted = True, fill_value="extrapolate")
    new_ydata = yInter(newTime)
    #new_ydata = yInter(time)
    
    #resample z data accordingly
    zInter = interp1d(time, z, kind='linear', assume_sorted = True,fill_value="extrapolate" )
    new_zdata = zInter(newTime)
    #new_zdata = zInter(time)
    
    #butterworth lowpass filter
    Wn = 2.0*fc/nbSamples
    sos = signal.butter(1, Wn, 'low', output='sos')
    #sos = signal.butter(1, Wn, 'low', output='ba',  fs=None)
    
    xdataFilter = max(xdata)*signal.sosfilt(sos, new_xdata) +xdata[0]
    ydataFilter = max(ydata)*signal.sosfilt(sos, new_ydata) +ydata[0]
    zdataFilter = signal.sosfilt(sos, new_zdata) +zdata[0]
    
    nx = []
    ny = []
    nz = []
    
    l = len(xdataFilter)
    
    for i in range(0, l):
        nx.append(xdataFilter[i] +(i/l) * (xdata[-1] - xdataFilter[-1])) 
        ny.append(ydataFilter[i] +(i/l) * (ydata[-1] -ydataFilter[-1])) 
        nz.append(zdataFilter[i] +(i/l) * (zdata[-1] -zdataFilter[-1])) 
    
    return (newTime, nx, ny, nz)





def resample3D(time, xdata, ydata, zdata, newFrameRate):
    
    #re s
    duration = time[-1] - time[0]
    nbSamples = (duration * newFrameRate) + 1
    #newTime = np.linspace(time[0], time[-1], len(time))
    newTime = np.linspace(time[0], time[-1], int(nbSamples))
    
    x = (xdata[:]-xdata[0])
    y = (ydata[:]-ydata[0])
    z = zdata[:]-zdata[0]
    
    #resample x data accordingly
    xInter = interp1d(time, x, kind='linear')
    new_xdata = xInter(newTime)
    #new_xdata = xInter(time)
    
    #resample y data accordingly
    yInter = interp1d(time, y, kind='linear')
    new_ydata = yInter(newTime)
    #new_ydata = yInter(time)
    
    #resample z data accordingly
    zInter = interp1d(time, z, kind='linear')
    new_zdata = zInter(newTime)
    #new_zdata = zInter(time)
    
    #butterworth lowpass filter
    Wn = 2.0*fc/nbSamples
    sos = signal.butter(1, Wn, 'low', output='sos')
    #sos = signal.butter(1, Wn, 'low', output='ba',  fs=None)
    
    xdataFilter = signal.sosfilt(sos, new_xdata) +xdata[0]
    ydataFilter = signal.sosfilt(sos, new_ydata) +ydata[0]
    zdataFilter = signal.sosfilt(sos, new_zdata) +zdata[0]

    return (newTime, xdataFilter, ydataFilter, zdataFilter)


def resample_spatial(xdata, ydata, newFrameRate):
    
    #re sampling on distance
    distArr = []
    
    #distArr.append(np.sqrt(xdata[0]*xdata[0] + ydata[0]*ydata[0])  )
    distArr.append(0)
    
    for i, x in enumerate(xdata[:-1]):
        dx = xdata[i+1] - xdata[i]
        dy = ydata[i+1] - ydata[i]
        
        dist = distArr[i] + np.sqrt(dx*dx + dy*dy)
        distArr.append(dist)
        
    newDist = np.linspace(distArr[0], distArr[-1], int(len(distArr) ))
    
    x = xdata[:]-xdata[0]
    y = ydata[:]-ydata[0]
    newDistArr = newDist[:]-newDist[0]
    #resample x data accordingly
    xInter = interp1d(distArr, x, kind='linear')
    new_xdata = xInter(newDistArr)
    
    #resample y data accordingly
    yInter = interp1d(distArr, y, kind='linear')
    new_ydata = yInter(newDistArr)
    
    #butterworth lowpass filter
    Wn = 2.0*fc/newFrameRate
    sos = signal.butter(1, Wn, 'low', output='sos')
    
    xdataFilter = signal.sosfilt(sos, new_xdata) +xdata[0]
    ydataFilter = signal.sosfilt(sos, new_ydata) +ydata[0]
    
    return (newDist, xdataFilter, ydataFilter)



def resample_spatial_all(tdata, xdata, ydata, zdata, yadata, rdata, pdata):
    
    #re sampling on distance
    distArr = []
    
    #distArr.append(np.sqrt(xdata[0]*xdata[0] + ydata[0]*ydata[0])  )
    distArr.append(0)
    
    for i, x in enumerate(xdata[1:]):
        dx = xdata[i+1] - xdata[i]
        dy = ydata[i+1] - ydata[i]
        
        dist = distArr[i] + np.sqrt(dx*dx + dy*dy)
        #dist = distArr[i] + dx + dy
        distArr.append(dist)
        
    newDist = np.linspace(distArr[0], distArr[-1], int(len(distArr) ))
    
    x = xdata[:]-xdata[0]
    y = ydata[:]-ydata[0]
    z = zdata[:]-zdata[0]
    ya = yadata[:]-yadata[0]
    r = rdata[:]-rdata[0]
    p = pdata[:]-pdata[0]
    t = tdata[:]-tdata[0]
    newDistArr = newDist[:]-newDist[0]
    
    #resample x data accordingly
    xInter = interp1d(distArr, x, kind='linear')
    new_xdata = xInter(newDistArr) +xdata[0]
    
    #resample y data accordingly
    yInter = interp1d(distArr, y, kind='linear')
    new_ydata = yInter(newDistArr) +ydata[0]
    
    #resample z data accordingly
    zInter = interp1d(distArr, z, kind='linear')
    new_zdata = zInter(newDistArr) +zdata[0]
    
    yaInter = interp1d(distArr, ya, kind='linear')
    new_yadata = yaInter(newDistArr) +yadata[0]
    
    rInter = interp1d(distArr, r, kind='linear')
    new_rdata = rInter(newDistArr) +rdata[0]
    
    pInter = interp1d(distArr, p, kind='linear')
    new_pdata = pInter(newDistArr) +pdata[0]
    
    tInter = interp1d(distArr, t, kind='linear')
    new_tdata = tInter(newDistArr) +tdata[0]
    
    
    return (new_tdata, new_xdata, new_ydata, new_zdata, new_yadata, new_rdata, new_pdata)

    

def visStat(statDF, fout):
        
    plt.clf()
    plt.figure(figsize=(10,10), dpi= 10)
    levelList = list(statDF.level.unique()) # uniquw=e levels list
    luserIdList= list(statDF.ID.unique()) 
    
    turnPointDF = pd.DataFrame() 
    titles = ["ID", "x", "y", "rep"]
    
    for l in levelList:
        titles.append(l)
    
    for i, userID in enumerate(luserIdList):
        xCoordArr = np.zeros(4 + len(levelList))
        yCoordArr = np.zeros(4 + len(levelList))
        xCoordArr[0] = userID
        xCoordArr[1] = 1 # x is used
        yCoordArr[0] = userID
        yCoordArr[2] = 1 # y coords
        
        repList= list(statDF.rep.unique()) 
        for rep in repList:
            
            xCoordArr[3] = rep # y coords
            yCoordArr[3] = rep # y coords
        
            for index, row in statDF[statDF["ID"] == userID ][statDF["rep"] == rep ].iterrows():
                lvlName = row['level']
                if lvlName in levelList:
                    lvlID = levelList.index(lvlName)
                 
                    xCoordArr[4+lvlID] = row["x"]
                    yCoordArr[4+lvlID] = row["y"]
                    
            ls = [xCoordArr,yCoordArr]
            df = pd.DataFrame(ls)
            turnPointDF =turnPointDF.append(df)
            
            userColor = i/len(luserIdList)
            
            tx = np.array(range(0, len(xCoordArr[4:])))
            
            plt.plot( xCoordArr[4:], color=[rep ,userColor,1 - userColor], linewidth= (2.2 + 2*rep), label=userID)
    
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12});    
    
    turnPointDF.columns = titles
    foutStat = os.path.join(fout+".csv")   
    turnPointDF.to_csv(foutStat)
    foutStat = os.path.join(fout+".xlsx") 
    turnPointDF.to_excel(foutStat,  index = False, header=True) 
    
    foutImg = os.path.join(fout + ".png")    
    plt.savefig(f'{foutImg}', dpi=my_dpi)
                    
    
    return



def visStatRepSeperate(statDF, fout):
        
    plt.clf()
    plt.figure(figsize=(10,10), dpi= 10)
    levelList = list(statDF.level.unique()) # uniquw=e levels list
    luserIdList= list(statDF.ID.unique()) 
    devRate = [0,0 ] # 1st component : if 3>5 -> true or if 6>4 -> true
    
    
    fig, axs = plt.subplots(1, 2, figsize=(20, 10), dpi=my_dpi)
    
    turnPointDF = pd.DataFrame() 
    titles = ["ID", "x", "y", "rep"]
    
    for l in levelList:
        titles.append(l)
    
    for i, userID in enumerate(luserIdList):
        xCoordArr = np.zeros(4 + len(levelList))
        yCoordArr = np.zeros(4 + len(levelList))
        xCoordArr[0] = userID
        xCoordArr[1] = 1 # x is used
        yCoordArr[0] = userID
        yCoordArr[2] = 1 # y coords
        
        repList= list(statDF.rep.unique()) 
        for rep in repList:
            
            xCoordArr[3] = rep # y coords
            yCoordArr[3] = rep # y coords
        
            for index, row in statDF[statDF["ID"] == userID ][statDF["rep"] == rep ].iterrows():
                lvlName = row['level']
                if lvlName in levelList:
                    lvlID = levelList.index(lvlName)
                 
                    xCoordArr[4+lvlID] = row["x"]
                    yCoordArr[4+lvlID] = row["y"]
                    
            ls = [xCoordArr,yCoordArr]
            df = pd.DataFrame(ls)
            turnPointDF =turnPointDF.append(df)
            
            #  if 3>5 -> true or if 6>4 -> true
            if xCoordArr[-1]>xCoordArr[-3]:
                devRate[0] += 1
            else:
                devRate[1] +=1
                
            if xCoordArr[-4]>xCoordArr[-2]:
                devRate[0] +=1
            else:
                devRate[1] +=1
            
            userColor = i/len(luserIdList)
            
            tx = np.array(range(0, len(xCoordArr[4:])))
            
            #axs[0,rep].plot( xCoordArr[4:], color=[userColor ,userColor,1 - userColor], linewidth= (2.2 + 2*rep), label=userID)
            x = xCoordArr[4:]
            axs[int(rep)].plot( x, color=[1 - userColor ,userColor,1 - userColor], linewidth= (2.2 ), label=userID)
                    
            #plt.plot( xCoordArr[4:], color=[userColor ,userColor,1 - userColor], linewidth= (2.2 + 2*rep), label=userID)
    
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12});    
    
    turnPointDF.columns = titles
    foutStat = os.path.join(fout+".csv")   
    turnPointDF.to_csv(foutStat)
    foutStat = os.path.join(fout+".xlsx") 
    turnPointDF.to_excel(foutStat,  index = False, header=True) 
    
    plt.title('Turning point for different levels')
    plt.xlabel('Level')
    plt.ylabel('X value')
 
# Create names on the x axis
    x_ticks = np.arange(len(levelList))

    plt.xticks(x_ticks, levelList)
    
    foutImg = os.path.join(fout + ".png")    
    plt.savefig(f'{foutImg}', dpi=my_dpi)
                    
    
    return devRate




def areaCalculation(pathDF, fout):
        
    borderDist = 30
    YborderMin = YEs+borderDist
    YborderMax = YNe-borderDist
    XborderMin = XEs
    
    levelList = list(pathDF.level.unique()) # uniquw=e levels list
    luserIdList= list(pathDF.ID.unique()) 
    
    areaAvoidDF = pd.DataFrame() 
    pointsAvoidDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in levelList:
        titles.append(l)
    
    for i, userID in enumerate(luserIdList):
        areaAvoid = np.zeros(startInd + len(levelList))
        areaAvoid[0] =  str(userID)
        
        
        pointsAvoid = np.zeros(startInd + len(levelList))
        pointsAvoid[0] = str(userID)
       
        repList= list(pathDF.rep.unique()) 
        
            
        plt.clf()
        plt.figure(figsize=(10,10), dpi= 10)
        fig, axs = plt.subplots(1, 2, figsize=(20, 10), dpi=my_dpi)
        for rep in repList:
            
            areaAvoid[1] = rep 
            pointsAvoid[1] = rep 
            
            for l, level in enumerate(levelList):
        
                
                levelColor = l/len(levelList)
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)]
                pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df[ (df["X"]>XNe) & ( df["Y"]> YborderMin) & (df["Y"]< YborderMax) ].copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if level in levelList:
                    lvlID = levelList.index(level)
                 
                    xSum = areaData['X'].sum() / len(areaData['X']) - XborderMin
                    
                    ySum = areaData['Y'].sum() / len(areaData['Y']) 
                    area = np.trapz( np.array(areaData['X']  - XborderMin), np.array(areaData['Y'])  )
                    area2 = np.trapz(  np.array(areaData['Y'] - YborderMax ),  np.array(areaData['X']) )  
                    
                    pointsAvoid[ startInd +lvlID ] = xSum
                    areaAvoid[startInd+lvlID ] = area
                    
                    #print(userID , level, ySum, area) #, area2)
                    
                    if lvlID > 2:
                        
                        axs[int(rep)].plot( np.array(pathPltDF['X']), np.array(pathPltDF['Y']) , color=[ levelColor ,1 -levelColor, 1 -levelColor], linewidth= (2.2 ), label=userID)
                        
                        axs[int(rep)].plot( np.array(areaData['X']), np.array(areaData['Y']) , color=[1 - levelColor ,levelColor, levelColor], linewidth= (4.2 ), label=userID)
                        
            ls = [pointsAvoid]
            df = pd.DataFrame(ls)
            pointsAvoidDF =pointsAvoidDF.append(df)
                    
            ls = [areaAvoid]
            df = pd.DataFrame(ls)
            areaAvoidDF =areaAvoidDF.append(df)
    
        lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12});    
        
        plt.title('Area point for different levels')
        plt.xlabel('Level')
        plt.ylabel('X value')
     
    # Create names on the x axis
    
        foutImg = os.path.join(fout + str(userID) + "_Area.png")    
        plt.savefig(f'{foutImg}', dpi=my_dpi)
        
    
    areaAvoidDF.columns = titles
    foutStat = os.path.join(fout+"_Area.csv")   
    areaAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_Area.xlsx") 
    areaAvoidDF.to_excel(foutStat,  index = False, header=True) 
          
    pointsAvoidDF.columns = titles
    foutStat = os.path.join(fout+"_pointsAvoid.csv")   
    pointsAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_pointsAvoid.xlsx") 
    pointsAvoidDF.to_excel(foutStat,  index = False, header=True) 

    

    # vis area
    
    plt.clf()
    plt.figure(figsize=(10,10), dpi= 10)
    fig, axs = plt.subplots(1, 2, figsize=(20, 10), dpi=my_dpi)
    
    areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("int")
    
    areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("string")
    
    luserIdList= list(areaAvoidDF.ID.unique()) 
    
    for i, userID in enumerate(luserIdList):
        userColor = i/len(luserIdList)
        for rep in repList:
            
            for index, row in areaAvoidDF[ (areaAvoidDF["ID"] == userID) & (areaAvoidDF["rep"] == rep) ].iterrows():
                
                xarr = row.to_numpy()
                x = xarr[startInd:]
                
                axs[int(rep)].plot( x, color=[1 - userColor ,userColor,1 - userColor], linewidth= (2.2 ), label=userID)
                
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12});    
    
    
    plt.title('Avoidence Area for different levels')
    plt.xlabel('Level')
    plt.ylabel('X value')
 
# Create names on the x axis
    x_ticks = np.arange(len(levelList))

    plt.xticks(x_ticks, levelList)
    
    foutImg = os.path.join(fout + "Area per level distribution.png")    
    plt.savefig(f'{foutImg}', dpi=my_dpi)
          
    
    
    # points Avoidence avarege
    
    plt.clf()
    plt.figure(figsize=(10,10), dpi= 10)
    fig, axs = plt.subplots(1, 2, figsize=(20, 10), dpi=my_dpi)
    
    pointsAvoidDF['ID'] = pointsAvoidDF['ID'].astype("int")
    
    pointsAvoidDF['ID'] = pointsAvoidDF['ID'].astype("string")
    
    luserIdList= list(pointsAvoidDF.ID.unique()) 
    
    for i, userID in enumerate(luserIdList):
        userColor = i/len(luserIdList)
        for rep in repList:
            
            for index, row in pointsAvoidDF[ (pointsAvoidDF["ID"] == userID) & (pointsAvoidDF["rep"] == rep) ].iterrows():
                
                xarr = row.to_numpy()
                x = xarr[startInd:]
                
                axs[int(rep)].plot( x, color=[1 - userColor ,userColor,1 - userColor], linewidth= (2.2 ), label=userID)
                
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12});    
    
    
    plt.title('Avoidence Area for different levels')
    plt.xlabel('Level')
    plt.ylabel('X value')
 
# Create names on the x axis
    x_ticks = np.arange(len(levelList))

    plt.xticks(x_ticks, levelList)
    
    foutImg = os.path.join(fout + "Avoidence points avarege per level distribution.png")    
    plt.savefig(f'{foutImg}', dpi=my_dpi)
          
    
    return 






plt.style.use('seaborn-whitegrid')

path_of_the_directory= os.path.join("C:\\Users\ypatotsk\Documents\Dev\pythonDev", "trackingData")
raw_data_csv_directory = os.path.join(path_of_the_directory, "rawData")
#path_of_the_directory= os.path.join("C:\\Dev\\pythonDev", "trackingData")

modified_csv_directory= os.path.join( path_of_the_directory, "csvCleanedData")

modified_csv_filt_dir = os.path.join( path_of_the_directory, "csvFilteredData")

img_directory= os.path.join( path_of_the_directory, "csvVisData")

stat_directory= os.path.join( path_of_the_directory, "statistics")
    
outputDir = os.path.join(path_of_the_directory, "img")

if not (os.path.exists(outputDir)) :
    os.mkdir(outputDir)
    
if not (os.path.exists(modified_csv_directory)) :
    os.mkdir(modified_csv_directory)

if not (os.path.exists(modified_csv_filt_dir)) :
    os.mkdir(modified_csv_filt_dir)

if not (os.path.exists(img_directory)) :
    os.mkdir(img_directory)
    
if not (os.path.exists(stat_directory)) :
    os.mkdir(stat_directory)
                    
                    
outputDir = os.path.join(path_of_the_directory, "imgSorted")
if not (os.path.exists(outputDir)) :
    os.mkdir(outputDir)

tXarr = []
tYarr = []
tZarr = []
            
tTarr = []
tYAarr = []
tParr = []
tRarr = []

iserIdList = []
extractId = []

for filename in os.listdir(raw_data_csv_directory): 
    ids = re.findall(r"[-+]?\d*\.*\d+", filename)
    id = ids[0]
    extractId.append(id)
    
iserIdList = set(extractId)
FPSTarget = 40
my_dpi = 100
nSlopes = 2
fc = 1.0
format_data = "%Y.%m.%d-%H.%M.%S.%f"
#format_data = "%Y.%m.%d-%H.%M.%S"
YlimMax = 150 #110
YlimMin = -170 #-130
XlimMax = -30
XlimMin = -250


XEs = -118 # ES character
YEs = -110 
XNe = -118 # Neurotic character
YNe = 90

#stat_data
statDF = pd.DataFrame(columns=['ID', 'level', 'x', 'y'])
resTSstatDF = pd.DataFrame(columns=['ID', 'level', 'x', 'y'])
resTSXYstatDF = pd.DataFrame(columns=['ID', 'level', 'x', 'y'])

#merged DF per experiment
mergedDF = pd.DataFrame()
mergedTSDF = pd.DataFrame()
mergedTSXYDF = pd.DataFrame()


for id in iserIdList:
    
    plt.figure(figsize=(20,25), dpi= my_dpi)
    levelColor = 0
    titles = []
    
    #fig, axs = plt.subplots(2, 3, figsize=(40, 40), dpi = my_dpi)
    imgCount = 0
    for filename in os.listdir(raw_data_csv_directory): 
        
        name, ext = os.path.splitext(filename)
         
        if (( id in filename ) and (ext == '.csv') and ("Coord" in filename) and not( "init" in filename ) and filename[-11:-10] != '0' )  :
        
            # destinguish repetitions
            rep = 404
            if ( "rep0" in filename ):
                rep = 0
            elif( "rep1" in filename ):
                rep = 1
            
            
            Xarr = []
            Yarr = []
            Zarr = []
            
            Tarr = []
            DTarr = []
            MSarr = []
                            
            YAarr = []
            Parr = []
            Rarr = []
            
            Xarr_DP1 = []
            Yarr_DP1 = []
            Xarr_DP2 = []
            Yarr_DP2 = []
            
            tableAngle = pd.DataFrame()
            tableTS = pd.DataFrame()
            
            lineAngle = ""    
            f = os.path.join(raw_data_csv_directory, filename)
            startIndex = 0 
            
            if os.path.isfile(f):
                
                # >>> Merge data in one DF
                fileCoord = f
                csvfileCoord = open(fileCoord) 
                    
                readerCoord = csv.reader(csvfileCoord)
                
                listCoord = ( list(readerCoord))
                
                tableCoord = pd.DataFrame(listCoord)
                for x in tableCoord.itertuples():
                    if x[1].count("=0.000") == 3:
                        startIndex = x[0]
                        
                
                fileStrAngle = filename[0:-9] + "Angle*.csv"
                
                findRes = find(fileStrAngle, raw_data_csv_directory)
                if len(findRes) > 0:
                    fileAngle = findRes[0]
                    
                    csvfileAngle = open(fileAngle) 
                        
                    readerAngle = csv.reader(csvfileAngle)
                    
                    listAngle = ( list(readerAngle))
                    
                    tableAngle = pd.DataFrame(listAngle)
                        

                fileStrTS = filename[0:-9] + "TS*.csv"
                findRes = find(fileStrTS, raw_data_csv_directory)
                
                if len(findRes) > 0:
                    fileTS = findRes[0]
                    csvfileTS = open(fileTS) 
                        
                    readerTS = csv.reader(csvfileTS)
                    
                    listTS = ( list(readerTS))
                    
                    tableTS = pd.DataFrame(listTS)
    
                    
                tempDF =tableCoord[startIndex:]
                
                
                for i, row in tempDF.iterrows():
                    
                    s = row[0]
                    coords = re.findall(r"[-+]?\d*\.*\d+", s)
                    X=float(coords[0])
                    Y=float(coords[1])
                    Z=float(coords[2])
                    
                    
                    if tableAngle.size>0:
                        s = tableAngle.iloc[i][0]
                        coords = re.findall(r"[-+]?\d*\.*\d+", s)
                        P=float(coords[0])
                        YA=float(coords[1])
                        R=float(coords[2])
                    
                    if tableTS.size>0:
                        s = tableTS.iloc[i][0]
                        dateTime = s
                        T= s.split("-",1)[1]
                           
                    if not((X< XlimMin) or (X> XlimMax) or (Y> YlimMax)  or (Y< YlimMin)):
                        Xarr.append(X)
                        Yarr.append(Y)
                        Zarr.append(Z)
                        
                        YAarr.append(YA)
                        Parr.append(P)
                        Rarr.append(R)
                        
                        Tarr.append(T)
                        DTarr.append(dateTime)
                    elif ((X> XlimMin) and (X!=0)):
                        print("!!!!!!")
                        print(i, X, Y, Z, T)
                        
            
            trackDF = pd.DataFrame({'X': Xarr,
                                    'Y': Yarr,
                                    'Z': Zarr,
                                    'YA': YAarr,
                                    'P': Parr,
                                    'R': Rarr,
                                    'T': Tarr,
                                    'DT': DTarr,
                                    })
            
            # <<< Merge data in one DF
            
            # >>> Restore path of HDMI
            indMed = int(trackDF[ (trackDF["X"]<- 50) ][ (trackDF["X"]>- 170) ][(trackDF["Y"]<50)][(trackDF["Y"]>-70)].index[-1])
            indStrat = 0
            indEnd = len(trackDF)
            
            if len(trackDF[ (trackDF["X"]<- 50) ][ (trackDF["X"]>- 170) ][(trackDF["Y"]<50)][ (trackDF["Y"]>-70) ] ):
                for i, r in trackDF.iloc[indMed:].iterrows():
                    X = r.X
                    Y = r.Y
                    if ((X< XlimMin) or (X> XlimMax) or (Y> YlimMax)  or (Y< YlimMin)):
                        indEnd = i
                        print("--end index allert--", indStrat, indEnd)
                        break
                    
                for i,r in trackDF[::-1].iterrows():
                    X = r.X
                    Y = r.Y
                    if ((X< XlimMin) or (X> XlimMax) or (Y> YlimMax)  or (Y< YlimMin)):
                        indStrat = i
                        print("--begin index allert--", indStrat, indEnd)
                        break
            
            trackDF = trackDF[indStrat:indEnd]
            #fout = os.path.join(modified_csv_directory, filename)
            #trackDF.to_csv(fout)
            # <<< Restore path of HDMI
                
            
            TDarr = np.array(trackDF['DT'])
            TimeSarr = []
            
            for i, time_data in enumerate(TDarr):
            
                date = datetime.strptime(time_data, format_data)
                timeVars = re.findall(r'\d+', time_data)
                #ts = date.timestamp() 
                ts = date.timestamp() - date.microsecond/1000000 + int(timeVars[-1]) / 1000
                
                TimeSarr.append(ts)
            
            trackDF['timestamp'] = TimeSarr
            trackDF = trackDF.sort_values('timestamp')
            
            fout = os.path.join(modified_csv_directory, filename)
            trackDF.to_csv(fout)    
            
            df = trackDF.copy()
            df["ID"] = id
            df['level']= filename[-18:-10]
            df['rep' ]=  rep
            mergedDF = mergedDF.append(df)
            
            Xarr = np.array(trackDF["X"])
            Yarr = np.array(trackDF["Y"])
            
            # <<< TS resample          
            # Linear Approximation of original signal
            lbl = "init " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                #model = np.polyfit(Xarr, Yarr, 1)
                                
                pwlf_model = pwlf.PiecewiseLinFit(Xarr, Yarr ,degree=1)
                breaks = pwlf_model.fit(nSlopes)
                
                slopes = pwlf_model.calc_slopes()
                
                randColor = random.randint(0,10)
                    
                statDF = statDF.append({
                                        'ID': id, 
                                        'level': filename[-18:-10],
                                        'rep' : rep,
                                        'x': breaks[1] , 
                                        'y': pwlf_model.predict(breaks[1])[0]
                                        
                                        }, ignore_index=True)
                
                plt.plot(Xarr, Yarr, '.', color=[0.1*randColor,1 -0.2*levelColor, 0.2*levelColor], linewidth=0.02 , label=filename[-18:-10])
                                #,figure=plt.figure());            
                randColor = random.randint(0,10)
                for xb in breaks:
                    pltInit = plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor,1 -0.2*levelColor, 0.2*levelColor], marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                    
                
                plt.plot(XEs, YEs, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
             markersize=20, markevery=12)
                plt.plot(XNe, YNe, color='red', marker='h', markerfacecolor='lightcoral', markeredgewidth=12,
             markersize=20, markevery=12)
                
                foutImg = os.path.join(img_directory,name+"_resample.png")    
                #plt.savefig(f'{foutImg}', dpi=my_dpi)
                
            
            
            # >>> Time + Space resample to a new FPS
            # Linear Approximation for filtered data
            lbl = "TXsample " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                
                (t ,x ,y, z, ya, p, r ) = resample_spatial_all(
                    np.array(trackDF["timestamp"]), 
                    np.array(trackDF["X"]), 
                    np.array(trackDF["Y"]),
                    np.array(trackDF["Z"]),
                    np.array(trackDF["YA"]), 
                    np.array(trackDF["P"]),
                    np.array(trackDF["R"])
                    )
                    
                    
                resTSXYDF = pd.DataFrame({'X': x,
                                        'Y': y,
                                        'Z': z,
                                        'YA': ya,
                                        'P': p,
                                        'R': r,
                                        'timestamp': t,
                                        })
                
                fout = os.path.join(modified_csv_filt_dir, filename+"TSXYresample.csv")
                resTSXYDF.to_csv(fout)
                
                
                df = resTSXYDF.copy()
                df["ID"] = id
                df['level']= filename[-18:-10]
                df['rep' ]=  rep
                mergedTSXYDF = mergedTSXYDF.append(df)

                pwlf_model = pwlf.PiecewiseLinFit(x, y ,degree=1)
                breaks = pwlf_model.fit(nSlopes)
                
                slopes = pwlf_model.calc_slopes()
                
                resTSXYstatDF = resTSXYstatDF.append({
                                        'ID': id, 
                                        'level': filename[-18:-10],
                                        'rep' : rep,
                                        'x': breaks[1] , 
                                        'y': pwlf_model.predict(breaks[1])[0]
                                        
                                        }, ignore_index=True)
                
                randColor = random.randint(0,10)
                for xb in breaks:
                    
                    #pltResTSXT = plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                    plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                
                
                #pltResTSXT = plt.plot(x, pwlf_model.predict(x), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= lbl)
                pltResTSXT = plt.plot(x,y, '.', color=[1,0.2*levelColor,1 - 0.2*levelColor], linewidth=0.18 , label="filtered "+filename[-18:-10])       
                plt.plot(x, pwlf_model.predict(x), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor])
                #plt.plot(x,y, '.', color=[0.524,0.2*levelColor,1 - 0.2*levelColor], linewidth=0.18 , label= lbl)       
                titles.append(filename[-18:-10])
                # <<< Time resample to a new FPS
                
            
            # >>> Time resample to a new FPS
            # Linear Approximation for filtered data
            lbl = "Tsample " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                
                (TResArr ,x ,y, z ) = resample3D(
                    np.array(resTSXYDF["timestamp"]), 
                    np.array(resTSXYDF["X"]), 
                    np.array(resTSXYDF["Y"]),
                    np.array(resTSXYDF["Z"]), FPSTarget)
                
                (TResArr ,ya ,p, r ) = resample3D(
                    np.array(resTSXYDF["timestamp"]), 
                    np.array(resTSXYDF["YA"]), 
                    np.array(resTSXYDF["P"]),
                    np.array(resTSXYDF["R"]), FPSTarget)
                
                    
                resTSDF = pd.DataFrame({'X': x,
                                        'Y': y,
                                        'Z': z,
                                        'YA': ya,
                                        'P': p,
                                        'R': r,
                                        'timestamp': TResArr,
                                        })
                fout = os.path.join(modified_csv_filt_dir, filename+"TSresample.csv")
                resTSDF.to_csv(fout)
                
                df = resTSDF.copy()
                df["ID"] = id
                df['level']= filename[-18:-10]
                df['rep' ]=  rep
                mergedTSDF = mergedTSDF.append(df)
                
                
                pwlf_model = pwlf.PiecewiseLinFit(x, y ,degree=1)
                breaks = pwlf_model.fit(nSlopes)
                
                slopes = pwlf_model.calc_slopes()
                
                resTSstatDF = resTSstatDF.append({
                                        'ID': id, 
                                        'level': filename[-18:-10],
                                        'rep' : rep,
                                        'x': breaks[1] , 
                                        'y': pwlf_model.predict(breaks[1])[0]
                                        
                                        }, ignore_index=True)
                
                randColor = random.randint(0,10)
                #for xb in breaks:
                    
                    #pltResTS = plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                
                #plt.plot(x, pwlf_model.predict(x), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= lbl)
                plt.plot(x,y, '.', color=[0.524,0.2*levelColor,1 - 0.2*levelColor], linewidth=0.18 , label="filtered "+filename[-18:-10])       
                #titles.append(filename[-18:-10])
                # <<< Time resample to a new FPS
                
            
            levelColor =+ 0.05
            
    

    
    foutImg = os.path.join(img_directory,id+"_lin"+str(nSlopes)+"_lowPass_" +str(fc) + ".png")      
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
    
    if len(lgnd.legendHandles) >0:
        lgnd.legendHandles[0]._legmarker.set_markersize(40)
        
        for l in lgnd.legendHandles:
            l._legmarker.set_markersize(40)
        #lgnd.legendHandles[1]._legmarker.set_markersize(40)
    
    plt.savefig(f'{foutImg}', dpi=my_dpi)
            
    
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+".csv")   
statDF.to_csv(foutStat)
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+".xlsx") 
statDF.to_excel(foutStat,  index = False, header=True) 


foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resampleTS"+str(FPSTarget)+".csv")   
resTSstatDF.to_csv(foutStat)
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resampleTS"+str(FPSTarget)+".xlsx") 
resTSstatDF.to_excel(foutStat,  index = False, header=True) 


foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resampleTSXY"+str(FPSTarget)+".csv")   
resTSXYstatDF.to_csv(foutStat)
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resampleTSXY"+str(FPSTarget)+".xlsx") 
resTSXYstatDF.to_excel(foutStat,  index = False, header=True) 


fout = os.path.join(stat_directory,"devStat"+"_lin"+str(nSlopes)+"_raw")   
correctRate = visStatRepSeperate(statDF, fout)
print(correctRate )

fout = os.path.join(stat_directory,"devStat"+"_lin"+str(nSlopes)+"_resTS"+str(FPSTarget))   
correctRate = visStatRepSeperate(resTSstatDF, fout)
print(correctRate )

fout = os.path.join(stat_directory,"devStat"+"_lin"+str(nSlopes)+"_resTSXY"+str(FPSTarget))   
correctRate = visStatRepSeperate(resTSXYstatDF, fout)    
print(correctRate )

fout = os.path.join(stat_directory,"_raw")   
areaCalculation(mergedDF, fout)

fout = os.path.join(stat_directory,"_resTSXY"+str(FPSTarget))   
areaCalculation(mergedTSXYDF, fout)
               
    
plt.clf()
plt.figure(figsize=(10,10), dpi= 10)
