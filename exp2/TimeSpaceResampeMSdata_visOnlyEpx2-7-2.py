# Main Data Analysis file for Exp 13-15 June 2023
# Version 2-4. Code with no fix to obst data resampling and "level in level list"
# 
# WORKING version dec 2023
#
# Included processing the case where there were no avoidance 
#
# Refer to TimeSpaceResampeMSdata.py for similar data analysis approaches.
# Working processing of raw data
# Analysis provided for
# - raw data 
# - time resampling + ascilations filtering
# - time resampling + ascilations filtering + space resampling
# Receive date time data with miliseconds in format %Y.%m.%d-%H.%M.%S.%f
# Metrics provided:
# - Turning point
# - Area
# - Avarage distance
# - Time to complete the task (in seconds)
#
# 03.01.2024
# Levels are sorted on speed inctease 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import seaborn as sns
import matplotlib.cm as cm

import re
import csv
import os
import fnmatch
from scipy import optimize
from scipy import interpolate
from scipy import signal
from scipy.interpolate import interp1d

from scipy.ndimage import gaussian_filter1d

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
    
    x = (xdata[:]-xdata[0])
    y = (ydata[:]-ydata[0])
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
    
    xdataFilter = signal.sosfilt(sos, new_xdata) +xdata[0]
    ydataFilter = signal.sosfilt(sos, new_ydata) +ydata[0]
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
    
    x = xdata[:]-xdata[0]
    y = ydata[:]-ydata[0]
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

def resample3D_time(time, xdata, ydata, zdata, newFrameRate):
    
    #re s
    duration = time[-1] - time[0]
    nbSamples = (duration * newFrameRate) + 1
    #newTime = np.linspace(time[0], time[-1], len(time))
    newTime = np.linspace(time[0], time[-1], int(nbSamples))
    
    x = xdata[:]-xdata[0]
    y = ydata[:]-ydata[0]
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
    
    xdataFilter = new_xdata +xdata[0]
    ydataFilter = new_ydata +ydata[0]
    zdataFilter = new_zdata +zdata[0]

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


def resample_spatial_all(tdata, xdata, ydata, zdata, yadata, rdata, pdata, oydata):
    
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
    
    x = xdata[:]#-xdata[0]
    y = ydata[:]#-ydata[0]
    z = zdata[:]#-zdata[0]
    ya = yadata[:]#-yadata[0]
    r = rdata[:]#-rdata[0]
    p = pdata[:]#-pdata[0]
    t = tdata[:]#-tdata[0]
    oy = oydata[:]
    newDistArr = newDist[:]#-newDist[0]
    
    #resample x data accordingly
    xInter = interp1d(distArr, x, kind='linear')
    new_xdata = xInter(newDistArr)# +xdata[0]
    
    #resample y data 
    yInter = interp1d(distArr, y, kind='linear')
    new_ydata = yInter(newDistArr)# +ydata[0]
    
    #resample z data 
    zInter = interp1d(distArr, z, kind='linear')
    new_zdata = zInter(newDistArr)# +zdata[0]
    
    yaInter = interp1d(distArr, ya, kind='linear')
    new_yadata = yaInter(newDistArr)# +yadata[0]
    
    rInter = interp1d(distArr, r, kind='linear')
    new_rdata = rInter(newDistArr)# +rdata[0]
    
    pInter = interp1d(distArr, p, kind='linear')
    new_pdata = pInter(newDistArr)# +pdata[0]
    
    tInter = interp1d(distArr, t, kind='linear')
    new_tdata = tInter(newDistArr)# +tdata[0]
    
    oyInter = interp1d(distArr, oy, kind='linear')
    new_oydata =oyInter(newDistArr)# +ydata[0]
    
    """        
    plt.clf()
    plt.figure(figsize=(40,40), dpi= 200)                         
    
    plt.plot(newDistArr, new_ydata , color=[0.1*randColor ,0.2*levelColor,0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', alpha = 0.5, markeredgewidth=2, markersize=4, markevery=2, label="XY")
    
    plt.plot(distArr, ydata-10 , color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', alpha = 0.5, markeredgewidth=2, markersize=4, markevery=2, label="TS")
    
    
    plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12});    
    
    foutImg = os.path.join(fout + str(rdata[0]) + "_ResExample_v2.pdf")     
    plt.savefig(foutImg, dpi=200)
    """
    
    return (new_tdata, new_xdata, new_ydata, new_zdata, new_yadata, new_rdata, new_pdata, new_oydata)
  
def visStat(fout):
    
    doorOpenDF
    
    
    for filename in os.listdir(raw_data_csv_directory): 
        
        name, ext = os.path.splitext(filename)
        
        rep = repNum(filename)
            
        #axs[rep,0] = plt.subplot2grid(shape=(1, 1), loc=(rep, 0), colspan=1, rowspan=1 , title = "raw")
        #axs[rep,1] = plt.subplot2grid(shape=(1, 1), loc=(rep, 1), colspan=1, rowspan=1, title = "Butterworth filter")
        
         
        if (( id in filename ) and (ext == '.csv') and ("Coord" in filename) ) :
            
            lvlID = constLevelList.index(filename[-18:-10])
            lbl = filename[-18:-10]
    
    return doorOpenDF 

def visStat(statDF, fout):
        
    plt.clf()
    plt.figure(figsize=(10,10), dpi= 10)
    levelList = list(statDF.level.unique()) # unique levels list
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
            #turnPointDF =turnPointDF.append(df)
            turnPointDF = pd.concat([turnPointDF, df ], ignore_index=True)
            
            
            userColor = i/len(luserIdList)
            
            tx = np.array(range(0, len(xCoordArr[4:])))
            
            plt.plot( xCoordArr[4:], color=[rep ,userColor,1 - userColor], linewidth= (2.2 + 2*rep), label=userID)
    
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12});    
    
    turnPointDF.columns = titles
    foutStat = os.path.join(fout+".csv")   
    turnPointDF.to_csv(foutStat)
    foutStat = os.path.join(fout+".xlsx") 
    turnPointDF.to_excel(foutStat,  index = False, header=True) 
    
    foutImg = os.path.join(fout + ".pdf")    
    plt.savefig(f'{foutImg}', dpi=4*my_dpi)
                    
    
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
    
    foutImg = os.path.join(fout + ".pdf")    
    plt.savefig(f'{foutImg}', dpi=cdpi)
           
    return devRate


# Calculate task completion Time
def completionTime(pathDF, fout):
        
    XborderMin = -280
    XborderMax = 280
    
    Xobst = 0
    Yobst = 0
    
    areaRate = [0,0]
    medRate = [0,0]
    
    luserIdList= list(pathDF.ID.unique()) 
    
    complTimeDF = pd.DataFrame()  
    trajLenDF = pd.DataFrame()  
    avgVelDF = pd.DataFrame()  
    
    pointsAvoidSumDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 
        
    
    for i, userID in enumerate(luserIdList):
        complTime = np.zeros(startInd + len(constLevelList))
        complTime[0] =  str(userID)
        
        
        trajLenArr = np.zeros(startInd + len(constLevelList))
        trajLenArr[0] =  str(userID)
        
        
        avgVelArr = np.zeros(startInd + len(constLevelList))
        avgVelArr[0] =  str(userID)
        
        
       
          
        for rep in repList:
            
            complTime[1] = rep  
            trajLenArr[1] = rep  
            avgVelArr[1] = rep  
            
            for l, level in enumerate(constLevelList):
        
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if not areaData.empty:
                    lvlID = constLevelList.index(level)
                 
                    dt = areaData['timestamp'].max() -areaData['timestamp'].min() 
                    
                    
                    lx = np.array(areaData['X']) 
                    ly = np.array(areaData['Y'])
                    
                         
                    trajLen = 0
                    avgVel = 0
                    
                    if len(ly)>0:
                        for ind in range(1, len(ly)) :
                            trajLen =trajLen+ np.sqrt( (ly[ind-1] - ly[ind])*(ly[ind-1] - ly[ind]) +  (lx[ind-1] - lx[ind])*(lx[ind-1] - lx[ind]) )
                             
                        avgVel = trajLen / len(ly)
                    
                    
                    complTime[ startInd +lvlID ] = dt
                    trajLenArr[ startInd +lvlID ] = trajLen
                    avgVelArr[ startInd +lvlID ] = avgVel
                    
                        
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'dt': [dt],
                                            'trajLen': [trajLen],
                                            'avgVel': [avgVel],
                                            
                                            })
                    
                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)
                    
                      
            ls = [complTime]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            complTimeDF = pd.concat([complTimeDF, df ], ignore_index=True)
            

    
            ls = [trajLenArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            trajLenDF = pd.concat([trajLenDF, df ], ignore_index=True)
            

            ls = [avgVelArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            avgVelDF = pd.concat([avgVelDF, df ], ignore_index=True)
            

    foutStat = os.path.join(fout+"_Time.xlsx") 
    complTimeDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"_trajLenDF.xlsx") 
    trajLenDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"_avgVelDF.xlsx") 
    avgVelDF.to_excel(foutStat,  index = False, header=True) 
        
    
    complTimeAllrepDF = pd.DataFrame(complTimeDF.groupby(['ID'], as_index = False).mean()) 
    
    trajLenAllrepDF = pd.DataFrame(trajLenDF.groupby(['ID'], as_index = False).mean()) 
    avgVelAllrepDF = pd.DataFrame(avgVelDF.groupby(['ID'], as_index = False).mean()) 
    
    
    foutStat = os.path.join(fout+"_complTimeAllrepDF.xlsx") 
    complTimeAllrepDF.to_excel(foutStat,  index = False, header=True) 
    foutStat = os.path.join(fout+"_trajLenAllrepDF.xlsx") 
    trajLenAllrepDF.to_excel(foutStat,  index = False, header=True) 
    foutStat = os.path.join(fout+"_avgVelAllrepDF.xlsx") 
    avgVelAllrepDF.to_excel(foutStat,  index = False, header=True) 
    
    
    
    
    
        
    # >>> Visualisation 
    
    fig, axs = plt.subplots(ncols=2, nrows=len(repList), figsize=(32, 32), dpi=cdpi)
                    
    
    # Create names on the x axis
    
    x_ticks = np.arange(len(constLevelList))
    
    # ------------------
    # Box with mustaches 
    plt.clf()
    
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    print(pointsAvoidSumDF)
    
    
    # >>> Viz all repetitions
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        
    ax = sns.boxplot(data=pointsAvoidSumDF ,  x="level", y = "dt", order=constLevelList, palette= colorUserList )
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot Task Completion Time per level_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    plt.clf()
    
    # <<< Viz all repetitions
    
    for rep in repList:
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "dt", order=constLevelList , palette= colorUserList )
            
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Task Completion Time per level_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "dt",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "dt", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot Task Completion Time per level_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
  
    



# Calculate task completion Time
def velocityStat(pathDF, fout):

    import matplotlib.cm as cm
    XborderMin = -284
    XborderMax = 290

    XAvoidence = -254

    Xobst = 0
    Yobst = 0

    XdecisionConst = -254

    areaRate = [0,0]
    medRate = [0,0]

    luserIdList= list(pathDF.ID.unique()) 

    complTimeDF = pd.DataFrame()  

    pointsAvoidSumDF = pd.DataFrame() 


    expTDF = pd.DataFrame() 
    expXDF = pd.DataFrame() 
    expYDF = pd.DataFrame() 

    pointsAvoidersDF = pd.DataFrame()
    
    expYAvoidersDF = pd.DataFrame() 
    expXAvoidersDF = pd.DataFrame() 
    expTAvoidersDF = pd.DataFrame() 
    
    expXExpectersDF = pd.DataFrame() 
    expYExpectersDF = pd.DataFrame() 
    expTExpectersDF = pd.DataFrame() 
    
    
    
    

    nSteps = 10

    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)

    for l in constLevelList:
        titles.append(l)



    repList= list(pathDF.rep.unique()) 

    for i, userID in enumerate(luserIdList):

        expTimeArr = np.zeros(startInd + len(constLevelList))
        expXCoordArr = np.zeros(startInd + len(constLevelList))
        expYCoordArr = np.zeros(startInd + len(constLevelList))

        fig, axs = plt.subplots(ncols=3, nrows=len(repList), figsize=(64, 48), dpi=2*cdpi)

        for rep in repList:

            expTimeArr[0] =  str(userID)
            expTimeArr[1] = rep 
            expXCoordArr[0] =  str(userID)
            expXCoordArr[1] = rep 
            expYCoordArr[0] =  str(userID)
            expYCoordArr[1] = rep 

            for l, level in enumerate(constLevelList):

                lvlID = constLevelList.index(level)

                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()

                print(userID , rep, level, df.size)
                if df.size <= 0 or ( '100' in level ):
                    continue

                areaData1 = df.copy()

                areaData = areaData1.sort_values(by=['timestamp'])

                t = np.array(areaData['timestamp'])

                tmin = t[0]

                tnorm = t - tmin

                x = np.array(areaData['X'])

                y = np.array(areaData['Y'])

                v = []


                for i, it in enumerate(t):

                    i0 = 0
                    i1 = 0

                    if ( (i> nSteps) and (i < len(t)-nSteps) ):
                        i0 = i-nSteps
                        i1 = i+nSteps
                        v.append( np.sqrt( (x[i1]-x[i0])*(x[i1]-x[i0]) + (y[i1]-y[i0])*(y[i1]-y[i0]) ) / np.abs(t[i1] - t[i0])  )
                    else:
                        v.append(0)


                varr = np.array(v)
                if len(varr) >0:

                    vnorm = 1/ np.max(v)

                    vmax = np.max(v)
                    vmaxInd = np.argmax(v)

                    v25perc = 0.25*vmax

                    v25percInd = 0

                    for ind, vval in enumerate(varr[vmaxInd:0:-1]):
                        if vval < v25perc:
                            v25percInd = vmaxInd - ind

                            break

                                        
                    #nonZeroInd= np.where(v > 0)
                    nonZeroInd = np.nonzero(varr > 0.001)
                    

                    expTimeArr[ startInd +lvlID ] = tnorm[v25percInd]
                    expXCoordArr[ startInd +lvlID ] = x[v25percInd]
                    expYCoordArr[ startInd +lvlID ] = y[v25percInd]

                    colors = cm.rainbow(vnorm * varr)

                    #axs[rep,0].scatter(x,y, c=colors, linewidth= (6.2 ) )   
                    axs[rep,0].set_title("Raw, rep = " + str(rep))

                    axs[rep,0].plot( x, y, color= colorUserList[l], alpha = 0.8, linewidth= 4.2, label= constLevelTitles[l] )
                    #        color= colorLevelListRed[l], alpha = 0.5, linewidth= 0.2 ,marker='h', markerfacecolor=colorUserList[l], markeredgewidth=12,
                    #        markersize=18, markevery=12, label= constLevelTitles[l] )


                    axs[rep,0].scatter(x[v25percInd],y[v25percInd], c="black", linewidth= (12.2 ) )  
                    axs[rep,0].scatter(x[vmaxInd],y[vmaxInd], c="magenta", linewidth= (12.2 ) )   

                    #axs[rep,0].legend( fontsize=12, loc='upper left' , ncol = 2, frameon=True, prop={'size': 10});
                    #axs[rep,0].legend(  loc='upper left' , ncol = 2 );


                    axs[rep,1].set_title("Velocity in time")


                    axs[rep,1].plot( tnorm, v, color= colorUserList[l], alpha = 0.8, linewidth= 4.2, label= constLevelTitles[l] )
                    axs[rep,1].legend(  loc='upper right' , ncol = 2 );

                    axs[rep,1].scatter(tnorm[v25percInd],v[v25percInd], color= "black", linewidth= (12.2 ) )   

                    axs[rep,1].scatter(tnorm[vmaxInd],v[vmaxInd], color= "magenta", linewidth= (12.2 ) ) 


                    axs[rep,2].set_title("Velocity in space")


                    axs[rep,2].plot( x[nonZeroInd], varr[nonZeroInd], color= colorUserList[l], alpha = 0.8, linewidth= 4.2, label= constLevelTitles[l] )
                    #axs[rep,2].legend(  loc='upper right' , ncol = 2 );


                    axs[rep,2].scatter(x[v25percInd],v[v25percInd], c="black", linewidth= (12.2 ) )  
                    axs[rep,2].scatter(x[vmaxInd],v[vmaxInd], c="magenta", linewidth= (12.2 ) )   



                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'expTime': [tnorm[v25percInd]],
                                            'expXCoord': [x[v25percInd]],
                                            'expYCoord': [y[v25percInd]],

                                            })


                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)

                    #if x[v25percInd] < XdecisionConst:

                    #pointsAvoidersDF =  pd.concat([pointsAvoidersDF, df ], ignore_index=True)



            ls = [expYCoordArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            expYDF = pd.concat([expYDF, df ], ignore_index=True)

            ls = [expXCoordArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            expXDF = pd.concat([expXDF, df ], ignore_index=True)


            ls = [expTimeArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            expTDF = pd.concat([expTDF, df ], ignore_index=True)


            #pointsAvoidSumDF
            
            
            dfAvoiders =  pointsAvoidSumDF[(pointsAvoidSumDF["ID"] == userID) & (pointsAvoidSumDF["rep"] == rep) 
                                & (pointsAvoidSumDF["expXCoord"] > XAvoidence) ].copy()
                    
            
            
            print(dfAvoiders )
            
            print(len(dfAvoiders))
            # exclude datasets with at least 3 examples of non-avoiders 
            if len(dfAvoiders) < 3:
                
                    
                ls = [expYCoordArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expYAvoidersDF = pd.concat([expYAvoidersDF, df ], ignore_index=True)
    
                ls = [expXCoordArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expXAvoidersDF = pd.concat([expXAvoidersDF, df ], ignore_index=True)
    
    
                ls = [expTimeArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expTAvoidersDF = pd.concat([expTAvoidersDF, df ], ignore_index=True)

            else:

    
                ls = [expYCoordArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expYExpectersDF = pd.concat([expYExpectersDF, df ], ignore_index=True)
    
                ls = [expXCoordArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expXExpectersDF = pd.concat([expXExpectersDF, df ], ignore_index=True)
    
    
                ls = [expTimeArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expTExpectersDF = pd.concat([expTExpectersDF, df ], ignore_index=True)
    
    



        foutImg = os.path.join(fout + userID+"nSteps_" + str(nSteps) + "_velocityProfile.pdf")      


        fig.savefig(f'{foutImg}', dpi=2*cdpi)



    foutStat = os.path.join(fout+"_expXExpectersDF.xlsx") 
    expXExpectersDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expYExpectersDF.xlsx") 
    expYExpectersDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expTExpectersDF.xlsx") 
    expTExpectersDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expXAvoidersDF.xlsx") 
    expXAvoidersDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expYAvoidersDF.xlsx") 
    expYAvoidersDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expTAvoidersDF.xlsx") 
    expTAvoidersDF.to_excel(foutStat,  index = False, header=True) 



    foutStat = os.path.join(fout+"_pointsAvoidSumDF.csv")   
    pointsAvoidSumDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_pointsAvoidSumDF.xlsx") 
    pointsAvoidSumDF.to_excel(foutStat,  index = False, header=True) 



    foutStat = os.path.join(fout+"_expYDF.csv")   
    expYDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_expYDF.xlsx") 
    expYDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expXDF.csv")   
    expXDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_expXDF.xlsx") 
    expXDF.to_excel(foutStat,  index = False, header=True) 


    foutStat = os.path.join(fout+"_expTDF.csv")   
    expTDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_expTDF.xlsx") 
    expTDF.to_excel(foutStat,  index = False, header=True) 


    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "expTime", order=constLevelList , palette= colorUserList )

        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 0.8*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Wait Time.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF,  x="level", y = "expTime",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF,  x="level", y = "expTime", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 0.8*ymax )

    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait Time.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))




    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "expXCoord", order=constLevelList , palette= colorUserList )

        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 0.8*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Wait X.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,40))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF,  x="level", y = "expXCoord",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF,  x="level", y = "expXCoord", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )


    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait X.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))





    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "expTime", order=constLevelList , palette= colorUserList )

        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    #ymin, ymax = plt.ylim()
    #plt.ylim(ymin - 0.5*ymin, 0.8*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Wait Time.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,40))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF, y = "expXCoord",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF, y = "expXCoord", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )


    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait X All.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))



    plt.clf()
    plt.figure(figsize=(40,40))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF, y = "expTime",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF, y = "expTime", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )


    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait Time All.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))

    return pointsAvoidSumDF


# >>> velocity plot





# Calculate task completion Time
def velocityStatPath(pathDF, fout):

    XborderMin = -290
    XborderMax = 290

    XAvoidence = -254

    Xobst = 0
    Yobst = 0

    XdecisionConst = -254

    areaRate = [0,0]
    medRate = [0,0]

    luserIdList= list(pathDF.ID.unique()) 

    complTimeDF = pd.DataFrame()  

    pointsAvoidSumDF = pd.DataFrame() 


    expTDF = pd.DataFrame() 
    expXDF = pd.DataFrame() 
    expYDF = pd.DataFrame() 

    pointsAvoidersDF = pd.DataFrame()
    
    expYAvoidersDF = pd.DataFrame() 
    expXAvoidersDF = pd.DataFrame() 
    expTAvoidersDF = pd.DataFrame() 
    
    expXExpectersDF = pd.DataFrame() 
    expYExpectersDF = pd.DataFrame() 
    expTExpectersDF = pd.DataFrame() 
    
    
    
    

    nSteps = 10

    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)

    for l in constLevelList:
        titles.append(l)



    repList= list(pathDF.rep.unique()) 

    for i, userID in enumerate(luserIdList):

        expTimeArr = np.zeros(startInd + len(constLevelList))
        expXCoordArr = np.zeros(startInd + len(constLevelList))
        expYCoordArr = np.zeros(startInd + len(constLevelList))

        fig, axs = plt.subplots(ncols=3, nrows=len(repList), figsize=(64, 48), dpi=2*cdpi)

        for rep in repList:

            expTimeArr[0] =  str(userID)
            expTimeArr[1] = rep 
            expXCoordArr[0] =  str(userID)
            expXCoordArr[1] = rep 
            expYCoordArr[0] =  str(userID)
            expYCoordArr[1] = rep 

            for l, level in enumerate(constLevelList):

                lvlID = constLevelList.index(level)

                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()

                print(userID , rep, level, df.size)
                if df.size <= 0 or ( '100' in level ):
                    continue

                areaData1 = df.copy()

                areaData = areaData1.sort_values(by=['timestamp'])

                t = np.array(areaData['timestamp'])

                tmin = t[0]

                tnorm = t - tmin

                x = np.array(areaData['X'])

                y = np.array(areaData['Y'])

                v = []


                for i, it in enumerate(t):

                    i0 = 0
                    i1 = 0

                    if ( (i> nSteps) and (i < len(t)-nSteps) ):
                        i0 = i-nSteps
                        i1 = i+nSteps
                        v.append( np.sqrt( (x[i1]-x[i0])*(x[i1]-x[i0]) + (y[i1]-y[i0])*(y[i1]-y[i0]) ) / np.abs(t[i1] - t[i0])  )
                    else:
                        v.append(0)


                varr = np.array(v)
                if len(varr) >0:

                    vnorm = 1/ np.max(v)

                    vmax = np.max(v)
                    vmaxInd = np.argmax(v)

                    v25perc = 0.25*vmax

                    v25percInd = 0

                    for ind, vval in enumerate(varr[vmaxInd:0:-1]):
                        if vval < v25perc:
                            v25percInd = vmaxInd - ind

                            break



                    expTimeArr[ startInd +lvlID ] = tnorm[v25percInd]
                    expXCoordArr[ startInd +lvlID ] = x[v25percInd]
                    expYCoordArr[ startInd +lvlID ] = y[v25percInd]

                    colors = cm.rainbow(vnorm * varr)

                    #axs[rep,0].scatter(x,y, c=colors, linewidth= (6.2 ) )   
                    axs[rep,0].set_title("Raw, rep = " + str(rep))

                    #axs[rep,0].plot( x, y, color= colorUserList[l], linewidth= 0.2 ,marker='h', label= constLevelTitles[l] )
                    axs[rep,0].plot( x, y, color= colorUserList[l], alpha = 0.8, linewidth= 4.2, label= constLevelTitles[l] )

                    axs[rep,0].scatter(x[v25percInd],y[v25percInd], c="black", linewidth= (10.2 ) )  
                    axs[rep,0].scatter(x[vmaxInd],y[vmaxInd], c="black", linewidth= (10.2 ) )   

                    axs[rep,0].legend( fontsize=12, loc='upper left' , ncol = 2, frameon=True, prop={'size': 10});
                    #axs[rep,0].legend(  loc='upper left' , ncol = 2 );


                    axs[rep,1].set_title("Velocity in time")


                    axs[rep,1].plot( tnorm, v, color= colorUserList[l], alpha = 0.8, linewidth= 4.2, label= constLevelTitles[l] )
                    axs[rep,1].legend(  loc='upper right' , ncol = 2 );

                    axs[rep,1].scatter(tnorm[v25percInd],v[v25percInd], color= "red", linewidth= (10.2 ) )   

                    axs[rep,1].scatter(tnorm[vmaxInd],v[vmaxInd], color= "red", linewidth= (10.2 ) ) 


                    axs[rep,2].set_title("Velocity in space")


                    axs[rep,2].plot( x, v, color= colorUserList[l], alpha = 0.8, linewidth= 4.2, label= constLevelTitles[l] )
                    axs[rep,2].legend(  loc='upper right' , ncol = 2 );


                    axs[rep,2].scatter(x[v25percInd],v[v25percInd], c="red", linewidth= (10.2 ) )  
                    axs[rep,2].scatter(x[vmaxInd],v[vmaxInd], c="red", linewidth= (10.2 ) )   



                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'expTime': [tnorm[v25percInd]],
                                            'expXCoord': [x[v25percInd]],
                                            'expYCoord': [y[v25percInd]],

                                            })


                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)

                    #if x[v25percInd] < XdecisionConst:

                    #pointsAvoidersDF =  pd.concat([pointsAvoidersDF, df ], ignore_index=True)



            ls = [expYCoordArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            expYDF = pd.concat([expYDF, df ], ignore_index=True)

            ls = [expXCoordArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            expXDF = pd.concat([expXDF, df ], ignore_index=True)


            ls = [expTimeArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            expTDF = pd.concat([expTDF, df ], ignore_index=True)


            #pointsAvoidSumDF
            
            
            dfAvoiders =  pointsAvoidSumDF[(pointsAvoidSumDF["ID"] == userID) & (pointsAvoidSumDF["rep"] == rep) 
                                & (pointsAvoidSumDF["expXCoord"] > XAvoidence) ].copy()
                    
            
            
            print(dfAvoiders )
            
            print(len(dfAvoiders))
            # exclude datasets with at least 3 examples of non-avoiders 
            if len(dfAvoiders) < 3:
                
                    
                ls = [expYCoordArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expYAvoidersDF = pd.concat([expYAvoidersDF, df ], ignore_index=True)
    
                ls = [expXCoordArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expXAvoidersDF = pd.concat([expXAvoidersDF, df ], ignore_index=True)
    
    
                ls = [expTimeArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expTAvoidersDF = pd.concat([expTAvoidersDF, df ], ignore_index=True)

            else:

    
                ls = [expYCoordArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expYExpectersDF = pd.concat([expYExpectersDF, df ], ignore_index=True)
    
                ls = [expXCoordArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expXExpectersDF = pd.concat([expXExpectersDF, df ], ignore_index=True)
    
    
                ls = [expTimeArr]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                expTExpectersDF = pd.concat([expTExpectersDF, df ], ignore_index=True)
    
    



        foutImg = os.path.join(fout + userID+"nSteps_" + str(nSteps) + "_velocityProfile.pdf")      


        fig.savefig(f'{foutImg}', dpi=2*cdpi)



    foutStat = os.path.join(fout+"_expXExpectersDF.xlsx") 
    expXExpectersDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expYExpectersDF.xlsx") 
    expYExpectersDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expTExpectersDF.xlsx") 
    expTExpectersDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expXAvoidersDF.xlsx") 
    expXAvoidersDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expYAvoidersDF.xlsx") 
    expYAvoidersDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expTAvoidersDF.xlsx") 
    expTAvoidersDF.to_excel(foutStat,  index = False, header=True) 



    foutStat = os.path.join(fout+"_pointsAvoidSumDF.csv")   
    pointsAvoidSumDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_pointsAvoidSumDF.xlsx") 
    pointsAvoidSumDF.to_excel(foutStat,  index = False, header=True) 



    foutStat = os.path.join(fout+"_expYDF.csv")   
    expYDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_expYDF.xlsx") 
    expYDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expXDF.csv")   
    expXDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_expXDF.xlsx") 
    expXDF.to_excel(foutStat,  index = False, header=True) 


    foutStat = os.path.join(fout+"_expTDF.csv")   
    expTDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_expTDF.xlsx") 
    expTDF.to_excel(foutStat,  index = False, header=True) 


    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "expTime", order=constLevelList , palette= colorUserList )

        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 0.8*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Wait Time.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF,  x="level", y = "expTime",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF,  x="level", y = "expTime", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 0.8*ymax )

    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait Time.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))




    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "expXCoord", order=constLevelList , palette= colorUserList )

        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 0.8*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Wait X.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,40))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF,  x="level", y = "expXCoord",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF,  x="level", y = "expXCoord", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )


    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait X.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))





    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "expTime", order=constLevelList , palette= colorUserList )

        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    #ymin, ymax = plt.ylim()
    #plt.ylim(ymin - 0.5*ymin, 0.8*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Wait Time.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,40))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF, y = "expXCoord",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF, y = "expXCoord", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )


    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait X All.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))



    plt.clf()
    plt.figure(figsize=(40,40))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF, y = "expTime",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF, y = "expTime", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )


    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait Time All.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))

    return pointsAvoidSumDF




# <<< velocity plot

# Calculate task completion Time
def velocityPlot(pathDF, fout):
        
    import matplotlib.cm as cm
    XborderMin = -290
    XborderMax = 290
    
    Xobst = 0
    Yobst = 0
    
    areaRate = [0,0]
    medRate = [0,0]
    
    luserIdList= list(pathDF.ID.unique()) 
    
    complTimeDF = pd.DataFrame()  
    
    pointsAvoidSumDF = pd.DataFrame() 
    
    nSteps = 10
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 

    for i, userID in enumerate(luserIdList):
        
        fig, axs = plt.subplots(ncols=1, nrows=len(repList), figsize=(32, 64), dpi=2*cdpi)
        
        for rep in repList:
            
            for l, level in enumerate(constLevelList):
      
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                t = np.array(areaData['timestamp'])
                
                x = np.array(areaData['X'])
                
                y = np.array(areaData['Y'])
                
                v = []
                
                
                for i, it in enumerate(t):
                    
                    i0 = 0
                    i1 = 0
                    
                    if ( (i> nSteps) and (i < len(t)-nSteps) ):
                        i0 = i-nSteps
                        i1 = i+nSteps
                        v.append( np.sqrt( (x[i1]-x[i0])*(x[i1]-x[i0]) + (y[i1]-y[i0])*(y[i1]-y[i0]) ) / np.abs(t[i1] - t[i0])  )
                    else:
                        v.append(0)
                
                print(userID , rep, level)
                    
                    
                
                varr = np.array(v)
                if len(varr) >0:
                    
                    vnorm = 1/ np.max(v)
                    
                    colors = cm.rainbow(vnorm * varr)
                    
                    axs[rep].scatter(x,y, c=colors, linewidth= (6.2 ) )   
                    
                    
                    
                    axs[rep].set_title("Raw, rep = " + str(rep))
                    
                    axs[rep].plot( x, y, color= colorLevelListRed[l], alpha = 0.5, linewidth= 0.2 ,marker='h', markerfacecolor=colorUserList[l], markeredgewidth=12,
                             markersize=18, markevery=12, label= titles[l] )
                    
                    #axs[rep].legend( fontsize=40, loc='lower left' , ncol = 2, frameon=True, prop={'size': 10});
                    axs[rep].legend(  loc='lower left' , ncol = 2 );
                    
                    
                
        foutImg = os.path.join(fout + userID+ "_velocity.pdf")      
        
        # not needed to bcp for now
        #fig.savefig(f'{foutImg}', dpi=2*cdpi)
        
        print(foutImg)
 




# Calculate task completion Time
def trajPlot(pathDF, fout):
        
    import matplotlib.cm as cm
    XborderMin = -290
    XborderMax = 290
    
    Xobst = 0
    Yobst = 0
    
    areaRate = [0,0]
    medRate = [0,0]
    
    luserIdList= list(pathDF.ID.unique()) 
    
    complTimeDF = pd.DataFrame()  
    
    pointsAvoidSumDF = pd.DataFrame() 
    
    nSteps = 10
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 

    for i, userID in enumerate(luserIdList):
        
        fig, axs = plt.subplots(ncols=1, nrows=len(repList), figsize=(32, 64), dpi=2*cdpi)
        
        for rep in repList:
            
            for l, level in enumerate(constLevelList):
      
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                t = np.array(areaData['timestamp'])
                
                x = np.array(areaData['X'])
                
                y = np.array(areaData['Y'])
                
                v = []
                
                
                for i, it in enumerate(t):
                    
                    i0 = 0
                    i1 = 0
                    
                    if ( (i> nSteps) and (i < len(t)-nSteps) ):
                        i0 = i-nSteps
                        i1 = i+nSteps
                        v.append( np.sqrt( (x[i1]-x[i0])*(x[i1]-x[i0]) + (y[i1]-y[i0])*(y[i1]-y[i0]) ) / np.abs(t[i1] - t[i0])  )
                    else:
                        v.append(0)
                
                print(userID , rep, level)
                    
                    
                
                varr = np.array(v)
                if len(varr) >0:
                    
                    vnorm = 1/ np.max(v)
                    
                    colors = cm.rainbow(vnorm * varr)
                    
                    
                    
                    axs[rep].set_title("Raw, rep = " + str(rep))
                    
                    axs[rep].plot( x, y, color= colorUserList[l], label= titles[l] )
                    
                    #axs[rep].legend( fontsize=40, loc='lower left' , ncol = 2, frameon=True, prop={'size': 10});
                    axs[rep].legend(  loc='lower left' , ncol = 2 );
                    
                    
                
        foutImg = os.path.join(fout + userID+ "_path.pdf")      
        
        
        fig.savefig(f'{foutImg}', dpi=2*cdpi)
        
        print(foutImg)

               
                
""" 
# Calculate task completion Time
def statisticsMinMax(pathDF, fout):
    
    XborderMin = -100
    XborderMax = 100
    
    Xobst = 0
    Yobst = 0
    
    areaRate = [0,0]
    medRate = [0,0]
    
    luserIdList= list(pathDF.ID.unique()) 
    
    minYDF = pd.DataFrame()  
    maxYDF = pd.DataFrame()  
    
    pointsAvoidSumDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 
        
    
    for i, userID in enumerate(luserIdList):
        minYarr = np.zeros(startInd + len(constLevelList))
        minYarr[0] =  str(userID)
        
        maxYarr = np.zeros(startInd + len(constLevelList))
        maxYarr[0] =  str(userID)
        
        
        for rep in repList:
            
            minYarr[1] = rep  
            
            maxYarr[1] = rep  
            
            for l, level in enumerate(constLevelList):
        
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if not areaData.empty:
                    lvlID = constLevelList.index(level)
                 
                    maxY = areaData['Y'].abs().max() 
                    minY = areaData['Y'].abs().min() 
                    
                    minYarr[ startInd +lvlID ] = minY
                          
                    maxYarr[ startInd +lvlID ] = maxY
                          
                    
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'maxY': [maxY],
                                            'minY': [minY],
                                            
                                            })
                    
                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)
                    
                
                         
            ls = [maxYarr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            maxYDF = pd.concat([maxYDF, df ], ignore_index=True)
            
            ls = [minYarr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            minYDF = pd.concat([minYDF, df ], ignore_index=True)
            
    
    foutStat = os.path.join(fout+"_Max.csv")   
    maxYDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_Max.xlsx") 
    maxYDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"_Min.csv")   
    minYDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_Min.xlsx") 
    minYDF.to_excel(foutStat,  index = False, header=True) 
    
    
    
        
    # >>> Visualisation    
    fig, axs = plt.subplots(ncols=2, nrows=len(repList), figsize=(32, 32), dpi=cdpi)
                    
    
    # Create names on the x axis
    x_ticks = np.arange(len(constLevelList))
    
    # ------------------
    # Box with mustaches 
    plt.clf()
    
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    
    for rep in repList:
        
        # >>> Max Lateral deviation vis
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "maxY", order=constLevelList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Max Lateral deviation per level_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "maxY",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "maxY",
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot Max Lateral deviation per level_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
        # >>> Min Lateral deviation vis
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "minY", 
                         order=constLevelList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Min Lateral deviation per level_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "minY",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "minY", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot Min Lateral deviation per level_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
          
    return
"""

# Calculate task completion Time
def statisticsMinMaxDist(pathDF, fout):
    
    XborderMin = -100
    XborderMax = 100
    
    Xobst = 0
    Yobst = 0
    
    areaRate = [0,0]
    medRate = [0,0]
    
    luserIdList= list(pathDF.ID.unique()) 
    
    minYDF = pd.DataFrame()  
    maxYDF = pd.DataFrame()  
    
    pointsAvoidSumDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 
        
    
    for i, userID in enumerate(luserIdList):
        minYarr = np.zeros(startInd + len(constLevelList))
        minYarr[0] =  str(userID)
        
        maxYarr = np.zeros(startInd + len(constLevelList))
        maxYarr[0] =  str(userID)
        
        
        for rep in repList:
            
            minYarr[1] = rep  
            
            maxYarr[1] = rep  
            
            for l, level in enumerate(constLevelList):
        
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if not areaData.empty:
                    lvlID = constLevelList.index(level)
                    
                    distArr = []
                    
                    x1 = areaData['X'].values
                    y1 = areaData['Y'].values
                    
                    d2 = np.square( x1 )  + np.square( y1 ) 
                    distArr = np.sqrt( d2 )
                    yAbsArr = np.abs( y1 )

                    
                 
                    maxD = yAbsArr.max() #max lateral
                    minD = distArr.min() #min eucleadian
                    
                    minYarr[ startInd +lvlID ] = minD
                          
                    maxYarr[ startInd +lvlID ] = maxD
                          
                    
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'maxY': [maxD],
                                            'minY': [minD],
                                            
                                            })
                    
                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)
                    
                
                         
            ls = [maxYarr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            maxYDF = pd.concat([maxYDF, df ], ignore_index=True)
            
            ls = [minYarr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            minYDF = pd.concat([minYDF, df ], ignore_index=True)
            
    
    foutStat = os.path.join(fout+"_MaxLateralDist.xlsx") 
    maxYDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"_MinDist.xlsx") 
    minYDF.to_excel(foutStat,  index = False, header=True) 
    
    
    maxAvgDF = pd.DataFrame(maxYDF.groupby(['ID'], as_index = False).mean()) 
    minAvgDF = pd.DataFrame(minYDF.groupby(['ID'], as_index = False).mean())  
    
    #pointsAvoidAllrepDF = pd.DataFrame(pointsAvoidSumDF.groupby(['ID'], as_index = False).mean()) 
    
    
    foutStat = os.path.join(fout+"__MaxLateralDistAllrepDF.xlsx") 
    maxAvgDF.to_excel(foutStat,  index = False, header=True) 
    foutStat = os.path.join(fout+"_MinDistAllrepDF.xlsx") 
    minAvgDF.to_excel(foutStat,  index = False, header=True) 
    
    
    
        
    # >>> Visualisation    
    fig, axs = plt.subplots(ncols=2, nrows=len(repList), figsize=(32, 32), dpi=cdpi)
                    
    
    # Create names on the x axis
    x_ticks = np.arange(len(constLevelList))
    
    # ------------------
    # Box with mustaches 
    plt.clf()
    
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    
    """
    foutStat = os.path.join(fout+"__MaxLateralDistAllrepDF.xlsx") 
    maxAvgDF.to_excel(foutStat,  index = False, header=True) 
    foutStat = os.path.join(fout+"_MinDistAllrepDF.xlsx") 
    minAvgDF.to_excel(foutStat,  index = False, header=True) 
    
    """
    
    # >>> Max Lateral deviation vis
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "maxY", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot__MaxLateralDistAllrepDF.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "minY", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot__MinDistAllrepDF.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    
    for rep in repList:
        
        # >>> Max Lateral deviation vis
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "maxY", order=constLevelList , palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Max deviation per level_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "maxY",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "maxY",
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot Max deviation per level_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
        # >>> Min Lateral deviation vis
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "minY", 
                         order=constLevelList , palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Min deviation per level_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "minY",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "minY", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot Min deviation per level_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
    
          
    return





# Calculate adea of deviation
def areaDist(pathDF, fout):

    XborderMin = -120
    XborderMax = 120

    Xobst = 0
    Yobst = 0

    areaRate = [0,0]
    medRate = [0,0]

    nSteps = 10

    luserIdList= list(pathDF.ID.unique()) 

    areaAvoidDF = pd.DataFrame() 
    lateralAvoidDF = pd.DataFrame() 
    pointsAvoidDF = pd.DataFrame() 
    zeroXDeviationDF = pd.DataFrame() 
    yOnZeroFromObstDF = pd.DataFrame() 


    areaAvoidAvgDF = pd.DataFrame() 
    lateralAvoidAvgDF = pd.DataFrame() 
    pointsAvoidAvgDF = pd.DataFrame() 
    zeroXDeviationAvgDF = pd.DataFrame() 
    yOnZeroFromObstAvgDF = pd.DataFrame() 


    pointsAvoidSumDF = pd.DataFrame() 


    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)

    for l in constLevelList:
        titles.append(l)


    repList= list(pathDF.rep.unique()) 


    for i, userID in enumerate(luserIdList):
        areaAvoid = np.zeros(startInd + len(constLevelList))
        areaAvoid[0] =  str(userID)


        pointsAvoid = np.zeros(startInd + len(constLevelList))
        pointsAvoid[0] = str(userID)

        lateralAvoid = np.zeros(startInd + len(constLevelList))
        lateralAvoid[0] = str(userID)

        zeroXDeviation = np.zeros(startInd + len(constLevelList))
        zeroXDeviation[0] =  str(userID)
        
        
        yOnZeroFromObstArr = np.zeros(startInd + len(constLevelList))
        yOnZeroFromObstArr[0] =  str(userID)
        
        

        for rep in repList:

            areaAvoid[1] = rep 
            pointsAvoid[1] = rep 
            lateralAvoid[1] = rep
            zeroXDeviation[1] = rep 
            yOnZeroFromObstArr[1] = rep 

            for l, level in enumerate(constLevelList):


                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()

                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()

                areaData1 = df.copy()

                areaData = areaData1.sort_values(by=['timestamp']).reset_index()

                if not areaData.empty:
                    lvlID = constLevelList.index(level)

                    xSum = areaData['X'].sum() / len(areaData['X']) 

                    ySum = areaData['Y'].sum() / len(areaData['Y']) 

                    lx = np.array(areaData['X']) 
                    ly = np.array(areaData['Y'])
                    
                    
                    lyObst = gaussian_filter1d( np.array(areaData["OY"]), 1)

                    
                    # >>> dist to Obst
                    x1 = areaData['X'].values
                    y1 = areaData['Y'].values
                    
                    y1ToObst = y1 - lyObst
                    d2 = np.square( x1 )  + np.square( y1ToObst ) 
                    distArr = np.sqrt( d2 )
                    avgDistToObst = np.average(distArr)
                    
                    d2 = np.square( x1 )  + np.square( y1 ) 
                    distArr = np.sqrt( d2 )
                    yAbsArr = np.abs( y1 )
                 
                    maxD = yAbsArr.max() #max lateral
                    minD = distArr.min() #min eucleadian
                    
                    # <<< dist to Obst


                    yOnZero = 0
                    yOnZeroFromObst = 0
                    

                    xzero = areaData['X'].abs().min(skipna = True, numeric_only = True)

                    zeroIndArr = areaData.index[areaData['X'].abs() ==xzero]

                    zeroInd = zeroIndArr[0]
                    if ((zeroInd> nSteps) and (zeroInd < len(ly)-nSteps )):
                        yOnZero = abs(float(ly[zeroInd-nSteps:zeroInd+nSteps].sum()) / len(ly[zeroInd-nSteps:zeroInd+nSteps]))
                        yOnZeroFromObst = abs( float(ly[zeroInd-nSteps:zeroInd+nSteps].sum())-float(lyObst[zeroInd-nSteps:zeroInd+nSteps].sum()) ) / len(ly[zeroInd-nSteps:zeroInd+nSteps])




                    avgDist = 0
                    avgLateralDist = 0
                    larea = 0
                    for ind in range(0, len(lx)) :
                        avgDist =avgDist+ np.sqrt( lx[ind]*lx[ind] +  ly[ind]*ly[ind] )
                        avgLateralDist =avgLateralDist+ np.sqrt( ly[ind]* ly[ind] )
                        if ind > 0:
                            larea = larea + abs(lx[ind]-lx[ind-1])*(abs(ly[ind-1])+abs(ly[ind]))/2

                    # Change integration
                    #area = np.trapz( np.array(areaData['X']) , np.array(areaData['Y'])  )
                    area = larea #np.trapz( np.array(areaData['Y']), np.array(areaData['X']) )

                    if len(lx)>0:
                        pointsAvoid[ startInd +lvlID ] = avgDist / len(lx)
                        lateralAvoid[ startInd +lvlID ] = avgLateralDist / len(lx)
                        zeroXDeviation[ startInd +lvlID ] = yOnZero
                        avgDist = avgDist / len(lx)
                        avgLateralDist = avgLateralDist / len(lx)
                        yOnZeroFromObstArr[ startInd +lvlID ] = yOnZeroFromObst
                    else: 
                        pointsAvoid[ startInd +lvlID ] = 0
                        lateralAvoid[ startInd +lvlID ] = 0
                        zeroXDeviation[ startInd +lvlID ] = 0
                        yOnZeroFromObstArr[ startInd +lvlID ] = 0
                        avgDist = 0
                        avgLateralDist = 0

                    area = np.abs(area)
                    areaAvoid[startInd+lvlID ] = area

                    # for summary


                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'area': [area/10000],
                                            'avgDist': [avgDist],
                                            'avgLateralDist': [avgLateralDist],
                                            'yOnZero': [yOnZero],
                                            'maxLateral': maxD,
                                            'minDist': minD,
                                            'yOnZeroFromObst':[yOnZeroFromObst],
                                            'avgDistToObst': avgDistToObst

                                            })

                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)



            ls = [pointsAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            pointsAvoidDF = pd.concat([pointsAvoidDF, df ], ignore_index=True)


            ls = [areaAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles
            df['ID'] = userID
            df['rep'] = rep
            #areaAvoidDF =areaAvoidDF.append(df)
            areaAvoidDF = pd.concat([areaAvoidDF, df ], ignore_index=True)


            ls = [lateralAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            lateralAvoidDF = pd.concat([lateralAvoidDF, df ], ignore_index=True)

            ls = [zeroXDeviation]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            zeroXDeviationDF = pd.concat([zeroXDeviationDF, df ], ignore_index=True)
             
            ls = [yOnZeroFromObstArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            yOnZeroFromObstDF = pd.concat([zeroXDeviationDF, df ], ignore_index=True)
            


    areaAvoidAvgDF = pd.DataFrame(areaAvoidDF.groupby(['ID'], as_index = False).mean())    
    lateralAvoidAvgDF = pd.DataFrame(lateralAvoidDF.groupby(['ID'], as_index = False).mean())    
    pointsAvoidAvgDF = pd.DataFrame(pointsAvoidDF.groupby(['ID'], as_index = False).mean())    
    zeroXDeviationAvgDF = pd.DataFrame(zeroXDeviationDF.groupby(['ID'], as_index = False).mean())    
    yOnZeroFromObstAvgDF = pd.DataFrame(yOnZeroFromObstDF.groupby(['ID'], as_index = False).mean())    


    """
    for i, userID in enumerate(luserIdList):

        df = areaAvoidDF[(areaAvoidDF["ID"] == userID)].copy()
        meanDF = df.mean()
        areaAvoidAvgDF = pd.concat([areaAvoidAvgDF, meanDF ], ignore_index=True)

        df = lateralAvoidDF[(lateralAvoidDF["ID"] == userID)].copy()
        meanDF = df.mean()
        lateralAvoidAvgDF = pd.concat([lateralAvoidAvgDF, meanDF ], ignore_index=True)

        df = pointsAvoidDF[(pointsAvoidDF["ID"] == userID)].copy()
        meanDF = df.mean()
        pointsAvoidAvgDF = pd.concat([pointsAvoidAvgDF, meanDF ], ignore_index=True)

        df = zeroXDeviationDF[(zeroXDeviationDF["ID"] == userID)].copy()
        meanDF = df.mean()
        zeroXDeviationAvgDF = pd.concat([zeroXDeviationAvgDF, meanDF ], ignore_index=True)
    """    

    foutStat = os.path.join(fout+"_areaAvgDF.csv")   
    areaAvoidAvgDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_areaAvgDF.xlsx") 
    areaAvoidAvgDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_lateralAvgDF.csv")   
    lateralAvoidAvgDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_lateralAvgDF.xlsx") 
    lateralAvoidAvgDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_pointsAvgDF.csv")   
    pointsAvoidAvgDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_pointsAvgDF.xlsx") 
    pointsAvoidAvgDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_zeroXDeviationAvgDF.csv")   
    zeroXDeviationAvgDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_zeroXDeviationAvgDF.xlsx") 
    zeroXDeviationAvgDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"yOnZeroFromObstAvgDF.xlsx") 
    yOnZeroFromObstAvgDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"yOnZeroFromObstDF.xlsx") 
    yOnZeroFromObstDF.to_excel(foutStat,  index = False, header=True) 
    


    foutStat = os.path.join(fout+"_Area.csv")   
    areaAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_Area.xlsx") 
    areaAvoidDF.to_excel(foutStat,  index = False, header=True) 

    #pointsAvoidDF['ID'] = pointsAvoidDF['ID'].astype("int")
    #pointsAvoidDF['ID'] = pointsAvoidDF['ID'].astype("string")
    foutStat = os.path.join(fout+"_pointsAvarageAvoid.csv")   
    pointsAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_pointsAvarageAvoid.xlsx") 
    pointsAvoidDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_lateralAvarageAvoid.csv")   
    lateralAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_lateralAvarageAvoid.xlsx") 
    lateralAvoidDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_lateralZeroXDeviationDF.csv")   
    zeroXDeviationDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_lateralZeroXDeviationDF.xlsx") 
    zeroXDeviationDF.to_excel(foutStat,  index = False, header=True) 

    # >>> Visualisation 

    fig, axs = plt.subplots(ncols=2, nrows=len(repList), figsize=(32, 32), dpi=cdpi)


    # Create names on the x axis

    x_ticks = np.arange(len(constLevelList))

    for i, userID in enumerate(luserIdList):
        userColor = i/len(luserIdList)
        for rep in repList:

            for index, row in areaAvoidDF[ (areaAvoidDF["ID"] == userID) & (areaAvoidDF["rep"] == rep) ].iterrows():

                xarr = row.to_numpy()
                x = xarr[startInd:]
                print("------areaAvoidDF----->", x)

                axs[rep, 0].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)

                axs[rep, 0].set(xticks=x_ticks, xticklabels=constLevelList)          


            for index, row in pointsAvoidDF[ (pointsAvoidDF["ID"] == userID) & (pointsAvoidDF["rep"] == rep) ].iterrows():

                xarr = row.to_numpy()
                x = xarr[startInd:]
                print("------pointsAvoidDF----->", x)
                axs[rep, 1].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)
                # Create names on the x axis

                axs[rep, 1].set(xticks=x_ticks, xticklabels=constLevelList) 
                axs[rep, 1].legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12}); 





    axs[0,0].set_title('Avoidence Area for different levels')
    axs[0,1].set_title('Average distance to the center of the obstacle origine ')

    foutImg = os.path.join(fout + "Area per level distribution.pdf")    
    plt.savefig(f'{foutImg}', dpi=cdpi)


    # ------------------
    # Box with mustaches 
    plt.clf()

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    # >>> ALL


    pointsAvoidSumDF_plot =  pd.DataFrame(pointsAvoidSumDF.groupby(['ID','level'], as_index = False).mean())    

    foutStat = os.path.join(fout+"pointsAvoidSumDF_plot.xlsx") 
    pointsAvoidSumDF_plot.to_excel(foutStat,  index = False, header=True)  
    foutStat = os.path.join(fout+"pointsAvoidSumDF_plot.csv")   
    pointsAvoidSumDF_plot.to_csv(foutStat)

    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF_plot ,  x="level", y = "area", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Area_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    sns.boxplot(data=pointsAvoidSumDF_plot, x="level", y = "avgDist", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot avgDist_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))

    #avgLateralDist      
    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    sns.boxplot(data=pointsAvoidSumDF_plot, x="level", y = "avgLateralDist", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot avgLateralDist_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))  

    # yOnZero
    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF_plot ,  x="level", y = "yOnZero", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Lateral_yOnZero_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)


    # yOnZero
    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF_plot ,  x="level", y = "yOnZeroFromObst", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Lateral_yOnZeroFromObst_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)


    # <<< All

    for rep in repList:


        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", order=constLevelList , palette= colorUserList )

        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                           jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.75 , order=constLevelList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "BoxPlot Area per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)

        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist", order=constLevelList , palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "BoxPlot avgDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))



        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()

        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )

        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )

        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "ViolinePlot avgDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))


        # >>> avgLateralDist
        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()

        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgLateralDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )

        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )

        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "ViolinePlot avgLateralDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))

        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist"
                    , order=constLevelList , palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "BoxPlot avgLateralDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))

        # <<< avgLateralDist

        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "area",
                            inner='box', scale='width', dodge=True, linewidth=1, order=constLevelList )

        sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                      jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )

        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "ViolinePlot area per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))


    # <<< Visualisation 


    return pointsAvoidSumDF






# calculate statistics on avoiders only 
# pointsAvoidSumDF
def areaDistAvoiders(pathDF, fout, lAvoidersDF):

    XborderMin = -120
    XborderMax = 120

    Xobst = 0
    Yobst = 0

    areaRate = [0,0]
    medRate = [0,0]

    nSteps = 5

    luserIdList= list(pathDF.ID.unique()) 

    areaAvoidDF = pd.DataFrame() 
    lateralAvoidDF = pd.DataFrame() 
    pointsAvoidDF = pd.DataFrame() 
    zeroXDeviationDF = pd.DataFrame() 

    pointsAvoidSumDF = pd.DataFrame() 


    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)

    for l in constLevelList:
        titles.append(l)


    repList= list(pathDF.rep.unique()) 


    for i, userID in enumerate(luserIdList):
        areaAvoid = np.zeros(startInd + len(constLevelList))
        areaAvoid[0] =  str(userID)


        pointsAvoid = np.zeros(startInd + len(constLevelList))
        pointsAvoid[0] = str(userID)

        lateralAvoid = np.zeros(startInd + len(constLevelList))
        lateralAvoid[0] = str(userID)

        zeroXDeviation = np.zeros(startInd + len(constLevelList))
        zeroXDeviation[0] =  str(userID)

        for rep in repList:

            areaAvoid[1] = rep 
            pointsAvoid[1] = rep 
            lateralAvoid[1] = rep
            zeroXDeviation[1] = rep 


            avoidersDataDF = lAvoidersDF[(lAvoidersDF["ID"] == userID) & 
                                         (lAvoidersDF["rep"] == rep)].copy()
            #>>> exclude some avoiders
            if len(avoidersDataDF) <= 9:
                continue

            for l, level in enumerate(constLevelList):


                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()

                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()

                areaData1 = df.copy()

                areaData = areaData1.sort_values(by=['timestamp']).reset_index()

                avoidersDataDF = lAvoidersDF[(lAvoidersDF["ID"] == userID) & (lAvoidersDF["rep"] == rep) & 
                            (lAvoidersDF["level"]  ==level)].copy()


                if not (areaData.empty or avoidersDataDF.empty) :
                    lvlID = constLevelList.index(level)

                    xSum = areaData['X'].sum() / len(areaData['X']) 

                    ySum = areaData['Y'].sum() / len(areaData['Y']) 

                    lx = np.array(areaData['X']) 
                    ly = np.array(areaData['Y'])


                    yOnZero = 0

                    xzero = areaData['X'].abs().min(skipna = True, numeric_only = True)

                    zeroIndArr = areaData.index[areaData['X'].abs() ==xzero]

                    zeroInd = zeroIndArr[0]
                    if ((zeroInd> nSteps) and (zeroInd < len(ly)-nSteps )):
                        yOnZero = abs(float(ly[zeroInd-nSteps:zeroInd+nSteps].sum()) / len(ly[zeroInd-nSteps:zeroInd+nSteps]))



                    avgDist = 0
                    avgLateralDist = 0
                    larea = 0
                    for ind in range(0, len(lx)) :
                        avgDist =avgDist+ np.sqrt( lx[ind]*lx[ind] +  ly[ind]*ly[ind] )
                        avgLateralDist =avgLateralDist+ np.sqrt( ly[ind]* ly[ind] )
                        if ind > 0:
                            larea = larea + abs(lx[ind]-lx[ind-1])*(abs(ly[ind-1])+abs(ly[ind]))/2

                    # Change integration
                    #area = np.trapz( np.array(areaData['X']) , np.array(areaData['Y'])  )
                    area = larea #np.trapz( np.array(areaData['Y']), np.array(areaData['X']) )

                    if len(lx)>0:
                        pointsAvoid[ startInd +lvlID ] = avgDist / len(lx)
                        lateralAvoid[ startInd +lvlID ] = avgLateralDist / len(lx)
                        zeroXDeviation[ startInd +lvlID ] = yOnZero
                        avgDist = avgDist / len(lx)
                        avgLateralDist = avgLateralDist / len(lx)
                    else: 
                        pointsAvoid[ startInd +lvlID ] = 0
                        lateralAvoid[ startInd +lvlID ] = 0
                        zeroXDeviation[ startInd +lvlID ] = 0
                        avgDist = 0
                        avgLateralDist = 0

                    area = np.abs(area)
                    areaAvoid[startInd+lvlID ] = area

                    # for summary


                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'area': [area/10000],
                                            'avgDist': [avgDist],
                                            'avgLateralDist': [avgLateralDist],
                                            'yOnZero': [yOnZero],

                                            })

                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)





            ls = [pointsAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            pointsAvoidDF = pd.concat([pointsAvoidDF, df ], ignore_index=True)


            ls = [areaAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles
            df['ID'] = userID
            df['rep'] = rep
            #areaAvoidDF =areaAvoidDF.append(df)
            areaAvoidDF = pd.concat([areaAvoidDF, df ], ignore_index=True)


            ls = [lateralAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            lateralAvoidDF = pd.concat([lateralAvoidDF, df ], ignore_index=True)

            ls = [zeroXDeviation]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            zeroXDeviationDF = pd.concat([zeroXDeviationDF, df ], ignore_index=True)



    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("int")  
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("string")





    foutStat = os.path.join(fout+"_pointsAvoidSumDF.csv")   
    pointsAvoidSumDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_pointsAvoidSumDF.xlsx") 
    pointsAvoidSumDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_Area.csv")   
    areaAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_Area.xlsx") 
    areaAvoidDF.to_excel(foutStat,  index = False, header=True) 


    #pointsAvoidDF['ID'] = pointsAvoidDF['ID'].astype("int")
    #pointsAvoidDF['ID'] = pointsAvoidDF['ID'].astype("string")
    foutStat = os.path.join(fout+"_pointsAvarageAvoid.csv")   
    pointsAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_pointsAvarageAvoid.xlsx") 
    pointsAvoidDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_lateralAvarageAvoid.csv")   
    lateralAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_lateralAvarageAvoid.xlsx") 
    lateralAvoidDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_lateralZeroXDeviationDF.csv")   
    zeroXDeviationDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_lateralZeroXDeviationDF.xlsx") 
    zeroXDeviationDF.to_excel(foutStat,  index = False, header=True) 

    # >>> Visualisation 

    fig, axs = plt.subplots(ncols=2, nrows=len(repList), figsize=(32, 32), dpi=cdpi)


    # Create names on the x axis

    x_ticks = np.arange(len(constLevelList))

    for i, userID in enumerate(luserIdList):
        userColor = i/len(luserIdList)
        for rep in repList:

            for index, row in areaAvoidDF[ (areaAvoidDF["ID"] == userID) & (areaAvoidDF["rep"] == rep) ].iterrows():

                xarr = row.to_numpy()
                x = xarr[startInd:]
                print("------areaAvoidDF----->", x)

                axs[rep, 0].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)

                axs[rep, 0].set(xticks=x_ticks, xticklabels=constLevelList)          


            for index, row in pointsAvoidDF[ (pointsAvoidDF["ID"] == userID) & (pointsAvoidDF["rep"] == rep) ].iterrows():

                xarr = row.to_numpy()
                x = xarr[startInd:]
                print("------pointsAvoidDF----->", x)
                axs[rep, 1].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)
                # Create names on the x axis

                axs[rep, 1].set(xticks=x_ticks, xticklabels=constLevelList) 
                axs[rep, 1].legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12}); 





    axs[0,0].set_title('Avoidence Area for different levels')
    axs[0,1].set_title('Average distance to the center of the obstacle origine ')

    foutImg = os.path.join(fout + "Area per level distribution.pdf")    
    plt.savefig(f'{foutImg}', dpi=cdpi)


    # ------------------
    # Box with mustaches 
    plt.clf()

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    # >>> ALL


    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["area"] > 0) ] ,  x="level", y = "area", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Area_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["avgDist"] > 0) ], x="level", y = "avgDist", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot avgDist_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))

    #avgLateralDist      
    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["avgLateralDist"] > 0) ], x="level", y = "avgLateralDist", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot avgLateralDist_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))  

    # yOnZero
    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["yOnZero"] > 0) ] ,  x="level", y = "yOnZero", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Lateral_yOnZero_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)


    # <<< All

    for rep in repList:


        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", order=constLevelList , palette= colorUserList )

        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                           jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.75 , order=constLevelList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "BoxPlot Area per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)

        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist", order=constLevelList , palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "BoxPlot avgDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))



        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()

        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )

        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )

        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "ViolinePlot avgDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))


        # >>> avgLateralDist
        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()

        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgLateralDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )

        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )

        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "ViolinePlot avgLateralDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))

        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist"
                    , order=constLevelList , palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "BoxPlot avgLateralDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))

        # <<< avgLateralDist

        plt.clf()
        plt.figure(figsize=(40,15))

        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "area",
                            inner='box', scale='width', dodge=True, linewidth=1, order=constLevelList )

        sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                      jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )

        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)

        foutImg = os.path.join(fout + "ViolinePlot area per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))


    # <<< Visualisation 


    return 




# Calculate adea of deviation , separate Avoiders from Expectors 
def areaDistForAvoidersOnly(pathDF, fout, lVelocityDF):
        
    XAvoidence = -254 #avoidence criteria 
    
    XborderMin = -120
    XborderMax = 120
    
    Xobst = 0
    Yobst = 0
    
    areaRate = [0,0]
    medRate = [0,0]
    
    nSteps = 5
    
    luserIdList= list(pathDF.ID.unique()) 
    
    areaAvoidDF = pd.DataFrame() 
    lateralAvoidDF = pd.DataFrame() 
    pointsAvoidDF = pd.DataFrame() 
    zeroXDeviationDF = pd.DataFrame() 
    maxlateralAvoidDF = pd.DataFrame()  
    distToObstAvoidDF = pd.DataFrame()     
    zeroXToObstAvoidDF = pd.DataFrame()     

    areaExpectDF = pd.DataFrame() 
    lateralExpectDF = pd.DataFrame() 
    pointsExpectDF = pd.DataFrame() 
    zeroXExpectDF = pd.DataFrame() 
    maxlateralExpectDF = pd.DataFrame()
    distToObstExpectDF = pd.DataFrame()     
    zeroXToObstExpectDF = pd.DataFrame()     
    
    
    pointsAvoidSumDF = pd.DataFrame() 
    pointsExpectSumDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 
        
    
    for i, userID in enumerate(luserIdList):
        areaAvoid = np.zeros(startInd + len(constLevelList))
        areaAvoid[0] =  str(userID)
        
        
        pointsAvoid = np.zeros(startInd + len(constLevelList))
        pointsAvoid[0] = str(userID)
        
        lateralAvoid = np.zeros(startInd + len(constLevelList))
        lateralAvoid[0] = str(userID)
        
        maxlateralAvoid = np.zeros(startInd + len(constLevelList))
        maxlateralAvoid[0] = str(userID)
       
        zeroXDeviation = np.zeros(startInd + len(constLevelList))
        zeroXDeviation[0] =  str(userID)
        
        distToObst = np.zeros(startInd + len(constLevelList))
        distToObst[0] = str(userID)
        
        zeroXToObst = np.zeros(startInd + len(constLevelList))
        zeroXToObst[0] = str(userID)
          
        for rep in repList:
            
            areaAvoid[1] = rep 
            pointsAvoid[1] = rep 
            lateralAvoid[1] = rep
            maxlateralAvoid[1] = rep
            zeroXDeviation[1] = rep 
            distToObst[1] = rep
            zeroXToObst[1] = rep 
            
               
            
            for l, level in enumerate(constLevelList):
        
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp']).reset_index()
                
                if not areaData.empty:
                    lvlID = constLevelList.index(level)
                 
                    xSum = areaData['X'].sum() / len(areaData['X']) 
                    
                    ySum = areaData['Y'].sum() / len(areaData['Y']) 
                    
                    lx = np.array(areaData['X']) 
                    ly = np.array(areaData['Y'])
                    
                    
                    yOnZero = 0
                    
                    xzero = areaData['X'].abs().min(skipna = True, numeric_only = True)
                        
                    zeroIndArr = areaData.index[areaData['X'].abs() ==xzero]
                        

                    lyObst = gaussian_filter1d( np.array(areaData["OY"]), 1)
                         
                    avgDistToObst = 0
                    yOnZeroFromObst = 0
                        
                    zeroInd = zeroIndArr[0]
                    if ((zeroInd> nSteps) and (zeroInd < len(ly)-nSteps )):
                        yOnZero = abs(float(ly[zeroInd-nSteps:zeroInd+nSteps].sum()) / len(ly[zeroInd-nSteps:zeroInd+nSteps]))
                        
                        yOnZeroFromObst = abs( float(ly[zeroInd-nSteps:zeroInd+nSteps].sum())-float(lyObst[zeroInd-nSteps:zeroInd+nSteps].sum()) ) / len(ly[zeroInd-nSteps:zeroInd+nSteps])

                            
                        
                         
                    avgDist = 0
                    avgLateralDist = 0
                    maxLateralVal = 0
                    larea = 0
                    for ind in range(0, len(lx)) :
                        
                        
                        avgDistToObst =avgDistToObst+ np.sqrt( (lyObst[ind] - ly[ind])*(lyObst[ind] - ly[ind]) +  lx[ind]*lx[ind] )
                        
                        
                        avgDist =avgDist+ np.sqrt( lx[ind]*lx[ind] +  ly[ind]*ly[ind] )
                        avgLateralDist =avgLateralDist+ np.sqrt( ly[ind]* ly[ind] )
                        if maxLateralVal < np.sqrt( ly[ind]* ly[ind] ):
                            maxLateralVal = np.sqrt( ly[ind]* ly[ind] )
                            
                        if ind > 0:
                            larea = larea + abs(lx[ind]-lx[ind-1])*(abs(ly[ind-1])+abs(ly[ind]))/2
  
                    # Change integration
                    #area = np.trapz( np.array(areaData['X']) , np.array(areaData['Y'])  )
                    area = larea #np.trapz( np.array(areaData['Y']), np.array(areaData['X']) )
                    
                    velocityDataDf = pd.DataFrame() 
                    
                    if not(level == 'level100'):
                        velocityDataDf = lVelocityDF[(lVelocityDF["ID"] == userID) 
                                                & (lVelocityDF["rep"] == rep) & (lVelocityDF["level"]  ==level) ].copy()
                        
                    print("-----------velocityDataDf------")
                    print(userID, rep, level)
                    
                    print(velocityDataDf)
                    
                    

                    if not velocityDataDf.empty:
                        
                        pointsAvoid[ startInd +lvlID ] = avgDist / len(lx)
                        lateralAvoid[ startInd +lvlID ] = avgLateralDist / len(lx)
                        zeroXDeviation[ startInd +lvlID ] = yOnZero
                        maxlateralAvoid[ startInd +lvlID ] = maxLateralVal
                        avgDist = avgDist / len(lx)
                        avgLateralDist = avgLateralDist / len(lx)
                        avgDistToObst = avgDistToObst/ len(lx)
                        distToObst[ startInd +lvlID ] = avgDistToObst
                        zeroXToObst[ startInd +lvlID ] = yOnZeroFromObst
                        
 
                    else: 
                        pointsAvoid[ startInd +lvlID ] = 0
                        lateralAvoid[ startInd +lvlID ] = 0
                        zeroXDeviation[ startInd +lvlID ] = 0
                        maxlateralAvoid[ startInd +lvlID ] = 0
                        distToObst[ startInd +lvlID ] = 0
                        zeroXToObst[ startInd +lvlID ] = 0
                        avgDist = 0
                        avgLateralDist = 0
                        avgDistToObst = 0
                        yOnZeroFromObst  = 0
                        
                        
                    area = np.abs(area)
                    areaAvoid[startInd+lvlID ] = area
                    
                    
                    lexpTime = 0   
                    lexpXCoord = 0  
                    lexpYCoord = 0   
                    
                    if not velocityDataDf.empty:
                    
                        lexpTime = np.array(velocityDataDf["expTime"])[0]   
                        lexpXCoord = np.array(velocityDataDf["expXCoord"])[0]   
                        lexpYCoord = np.array(velocityDataDf["expYCoord"])[0]  
                    
                    # for summary
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'area': [area/10000],
                                            'avgDist': [avgDist],
                                            'avgLateralDist': [avgLateralDist],
                                            
                                            'maxLateralDist': [maxLateralVal],
                                            'yOnZero': [yOnZero],
                                            'expTime':[lexpTime],
                                            'expXCoord':[lexpXCoord],
                                            'expYCoord':[lexpYCoord],
                                            
                                            'avgDistToObst' : [avgDistToObst],
                                            'yOnZeroFromObst'  : [yOnZeroFromObst],
                                            
                                            })
                    
                    
                    dfAvoiders =  lVelocityDF[(lVelocityDF["ID"] == userID) & (lVelocityDF["rep"] == rep) 
                                & (lVelocityDF["expXCoord"] > XAvoidence) ].copy()
                    
                    # exclude datasets with at least 3 examples of non-avoiders 
                    if len(dfAvoiders) < 3 :
                    
                        pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)
                    else:
                        pointsExpectSumDF = pd.concat([pointsExpectSumDF, df ], ignore_index=True)    
                    
                    
            dfAvoiders =  lVelocityDF[(lVelocityDF["ID"] == userID) & (lVelocityDF["rep"] == rep) 
                        & (lVelocityDF["expXCoord"] > XAvoidence) ].copy()
            
            
            # exclude datasets with at least 3 examples of non-avoiders 
            if len(dfAvoiders) < 3:
                  
                             
                ls = [pointsAvoid]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                pointsAvoidDF = pd.concat([pointsAvoidDF, df ], ignore_index=True)
                
                
                ls = [areaAvoid]
                df = pd.DataFrame(ls)
                df.columns = titles
                df['ID'] = userID
                df['rep'] = rep
                #areaAvoidDF =areaAvoidDF.append(df)
                areaAvoidDF = pd.concat([areaAvoidDF, df ], ignore_index=True)
                
                             
                ls = [lateralAvoid]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                lateralAvoidDF = pd.concat([lateralAvoidDF, df ], ignore_index=True)
                
                
                ls = [maxlateralAvoid]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                maxlateralAvoidDF = pd.concat([maxlateralAvoidDF, df ], ignore_index=True)
                
                ls = [zeroXDeviation]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                zeroXDeviationDF = pd.concat([zeroXDeviationDF, df ], ignore_index=True)
                
                ls = [distToObst]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                distToObstAvoidDF = pd.concat([distToObstAvoidDF, df ], ignore_index=True)
    
                ls = [zeroXToObst]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                zeroXToObstAvoidDF = pd.concat([zeroXToObstAvoidDF, df ], ignore_index=True)
                
            else:
                
                        
                ls = [pointsAvoid]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                pointsExpectDF = pd.concat([pointsExpectDF, df ], ignore_index=True)
                
                
                ls = [areaAvoid]
                df = pd.DataFrame(ls)
                df.columns = titles
                df['ID'] = userID
                df['rep'] = rep
                #areaAvoidDF =areaAvoidDF.append(df)
                areaExpectDF = pd.concat([areaExpectDF, df ], ignore_index=True)
                
                             
                ls = [lateralAvoid]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                lateralExpectDF = pd.concat([lateralExpectDF, df ], ignore_index=True)
                
                
                  
                ls = [maxlateralAvoid]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                maxlateralExpectDF = pd.concat([maxlateralExpectDF, df ], ignore_index=True)
                
                ls = [zeroXDeviation]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                zeroXExpectDF = pd.concat([zeroXExpectDF, df ], ignore_index=True)
    
    
                ls = [distToObst]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                distToObstExpectDF = pd.concat([distToObstExpectDF, df ], ignore_index=True)
    
                ls = [zeroXToObst]
                df = pd.DataFrame(ls)
                df.columns = titles         
                df['ID'] = userID
                df['rep'] = rep
                #pointsAvoidDF =pointsAvoidDF.append(df)
                zeroXToObstExpectDF = pd.concat([zeroXToObstExpectDF, df ], ignore_index=True)
    
    
    areaAvoidAvoidersAvgDF = pd.DataFrame(areaAvoidDF.groupby(['ID'], as_index = False).mean())    
    lateralAvoidAvoidersAvgDF = pd.DataFrame(lateralAvoidDF.groupby(['ID'], as_index = False).mean())   
    maxlateralAvoidAvoidersAvgDF = pd.DataFrame(maxlateralAvoidDF.groupby(['ID'], as_index = False).mean())    
    pointsAvoidAvoidersAvgDF = pd.DataFrame(pointsAvoidDF.groupby(['ID'], as_index = False).mean())    
    zeroXDeviationAvoidersAvgDF = pd.DataFrame(zeroXDeviationDF.groupby(['ID'], as_index = False).mean()) 
    distToObstAvoidersAvgDF = pd.DataFrame(distToObstAvoidDF.groupby(['ID'], as_index = False).mean())    
    zeroXToObstAvoidersAvgDF = pd.DataFrame(zeroXToObstAvoidDF.groupby(['ID'], as_index = False).mean())       

    
    pointsExpectAvgDF = pd.DataFrame(pointsExpectDF.groupby(['ID'], as_index = False).mean())    
    areaExpectAvgDF = pd.DataFrame(areaExpectDF.groupby(['ID'], as_index = False).mean())    
    lateralExpecAvgDF = pd.DataFrame(lateralExpectDF.groupby(['ID'], as_index = False).mean())  
    maxlateralExpecAvgDF = pd.DataFrame(maxlateralExpectDF.groupby(['ID'], as_index = False).mean())    
    zeroXExpectAvgDF = pd.DataFrame(zeroXExpectDF.groupby(['ID'], as_index = False).mean())   
    distToObstExpectAvgDF = pd.DataFrame(distToObstExpectDF.groupby(['ID'], as_index = False).mean())    
    zeroXToObstExpectAvgDF = pd.DataFrame(zeroXToObstExpectDF.groupby(['ID'], as_index = False).mean())    
    
    
    
    foutStat = os.path.join(fout+"__pointsAvoidSumDF.xlsx") 
    pointsAvoidSumDF.to_excel(foutStat,  index = False, header=True)     
    foutStat = os.path.join(fout+"__pointsExpectSumDF.xlsx") 
    pointsExpectSumDF.to_excel(foutStat,  index = False, header=True)    

    
    foutStat = os.path.join(fout+"_areaAvoidAvoidersAvgDF.xlsx") 
    areaAvoidAvoidersAvgDF.to_excel(foutStat,  index = False, header=True)     
    foutStat = os.path.join(fout+"_lateralAvoidAvoidersAvgDF.xlsx") 
    lateralAvoidAvoidersAvgDF.to_excel(foutStat,  index = False, header=True)    
    
    foutStat = os.path.join(fout+"_maxlateralAvoidAvoidersAvgDF.xlsx") 
    maxlateralAvoidAvoidersAvgDF.to_excel(foutStat,  index = False, header=True)    
    
    foutStat = os.path.join(fout+"_pointsAvoidAvoidersAvgDF.xlsx") 
    pointsAvoidAvoidersAvgDF.to_excel(foutStat,  index = False, header=True)       
    foutStat = os.path.join(fout+"_zeroXDeviationAvoidersAvgDF.xlsx") 
    zeroXDeviationAvoidersAvgDF.to_excel(foutStat,  index = False, header=True)     
  
    foutStat = os.path.join(fout+"_pointsExpectAvgDF.xlsx") 
    pointsExpectAvgDF.to_excel(foutStat,  index = False, header=True)    
    foutStat = os.path.join(fout+"_areaExpectAvgDF.xlsx") 
    areaExpectAvgDF.to_excel(foutStat,  index = False, header=True)    
    foutStat = os.path.join(fout+"_lateralExpecAvgDF.xlsx") 
    lateralExpecAvgDF.to_excel(foutStat,  index = False, header=True)  
    foutStat = os.path.join(fout+"_maxlateralExpecAvgDF.xlsx") 
    maxlateralExpecAvgDF.to_excel(foutStat,  index = False, header=True)  
    foutStat = os.path.join(fout+"_zeroXExpectAvgDF.xlsx") 
    zeroXExpectAvgDF.to_excel(foutStat,  index = False, header=True)  


    
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("int")  
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("string")
    
    foutStat = os.path.join(fout+"_AreaAvoiders.csv")   
    areaAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_AreaAvoiders.xlsx") 
    areaAvoidDF.to_excel(foutStat,  index = False, header=True) 
          

    #pointsAvoidDF['ID'] = pointsAvoidDF['ID'].astype("int")
    #pointsAvoidDF['ID'] = pointsAvoidDF['ID'].astype("string")

    foutStat = os.path.join(fout+"_pointsAvarageAvoidAvoiders.xlsx") 
    pointsAvoidDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"_lateralAvarageAvoidAvoiders.xlsx") 
    lateralAvoidDF.to_excel(foutStat,  index = False, header=True) 
    
    
    foutStat = os.path.join(fout+"_maxlateralAvarageAvoidAvoiders.xlsx") 
    maxlateralAvoidDF.to_excel(foutStat,  index = False, header=True) 
            
    foutStat = os.path.join(fout+"_lateralZeroXDeviationDFAvoiders.csv")   
    zeroXDeviationDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_lateralZeroXDeviationDFAvoiders.xlsx") 
    zeroXDeviationDF.to_excel(foutStat,  index = False, header=True) 
    
    
    foutStat = os.path.join(fout+"_AreaExpect.xlsx") 
    areaExpectDF.to_excel(foutStat,  index = False, header=True) 
          

    foutStat = os.path.join(fout+"_pointsAvarageExpect.xlsx") 
    pointsExpectDF.to_excel(foutStat,  index = False, header=True) 
    

    foutStat = os.path.join(fout+"_lateralAvarageExpect.xlsx") 
    lateralExpectDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_maxlateralAvarageExpect.xlsx") 
    maxlateralExpectDF.to_excel(foutStat,  index = False, header=True) 
            

    foutStat = os.path.join(fout+"_lateralZeroXExpectDF.xlsx") 
    zeroXExpectDF.to_excel(foutStat,  index = False, header=True) 
    
    
    
    foutStat = os.path.join(fout+"_pointsAvoidOnlySumDF.xlsx") 
    pointsAvoidSumDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"_pointsExpectOnlySumDF.xlsx") 
    pointsExpectSumDF.to_excel(foutStat,  index = False, header=True) 
    
    
    
    foutStat = os.path.join(fout+"_distToObstAvoidersAvgDF.xlsx") 
    distToObstAvoidersAvgDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"_zeroXToObstAvoidersAvgDF.xlsx") 
    zeroXToObstAvoidersAvgDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"_distToObstExpectAvgDF.xlsx") 
    distToObstExpectAvgDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"_zeroXToObstExpectAvgDF.xlsx") 
    zeroXToObstExpectAvgDF.to_excel(foutStat,  index = False, header=True) 
    
    
    
    # >>> Visualisation 
    
    print(" - - pointsExpectSumDF - - ")
    print(pointsExpectSumDF)
    
    
    # >>> vis expTime vs avgDist
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
        
    x = np.array(pointsAvoidSumDF["expTime"])
    y = np.array(pointsAvoidSumDF["avgDist"])
    print("-----------x-----y-------")
    print(x,y)
    ax.scatter(x,y, color= 'blue' )
    
    x = np.array(pointsExpectSumDF["expTime"])
    y = np.array(pointsExpectSumDF["avgDist"])
    print("-----------x-----y-------")
    print(x,y)
    ax.scatter(x,y, color= 'red' )
    
     
    foutImg = os.path.join(fout + "expTime vs avgDist_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    plt.clf()
    
    # >>> vis expTime vs avgDist
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
        
    x = np.array(pointsAvoidSumDF["expTime"])
    y = np.array(pointsAvoidSumDF["yOnZero"])
    print("-----------x-----y-------")
    print(x,y)
    ax.scatter(x,y, color= 'blue' )
    
    x = np.array(pointsExpectSumDF["expTime"])
    y = np.array(pointsExpectSumDF["yOnZero"])
    print("-----------x-----y-------")
    print(x,y)
    ax.scatter(x,y, color= 'red' )
    
     
    foutImg = os.path.join(fout + "expTime vs yOnZero_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    plt.clf()
    
    
    # <<< vis expTime vs avgDist
    
    
    
    
    
    
    
    
    fig, axs = plt.subplots(ncols=2, nrows=len(repList), figsize=(32, 32), dpi=cdpi)
                    
    
    # Create names on the x axis
    
    x_ticks = np.arange(len(constLevelList))
    
    for i, userID in enumerate(luserIdList):
        userColor = i/len(luserIdList)
        for rep in repList:
            
            for index, row in areaAvoidDF[ (areaAvoidDF["ID"] == userID) & (areaAvoidDF["rep"] == rep) ].iterrows():
                
                xarr = row.to_numpy()
                x = xarr[startInd:]
                print("------areaAvoidDF----->", x)
                
                axs[rep, 0].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)
                
                axs[rep, 0].set(xticks=x_ticks, xticklabels=constLevelList)          
                
                
            for index, row in pointsAvoidDF[ (pointsAvoidDF["ID"] == userID) & (pointsAvoidDF["rep"] == rep) ].iterrows():
                
                xarr = row.to_numpy()
                x = xarr[startInd:]
                print("------pointsAvoidDF----->", x)
                axs[rep, 1].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)
                # Create names on the x axis
                
                axs[rep, 1].set(xticks=x_ticks, xticklabels=constLevelList) 
                axs[rep, 1].legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12}); 
                
  

    
    
    axs[0,0].set_title('Avoidence Area for different levels')
    axs[0,1].set_title('Average distance to the center of the obstacle origine ')
    
    foutImg = os.path.join(fout + "Area per level distribution.pdf")    
    plt.savefig(f'{foutImg}', dpi=cdpi)
        
    
    # ------------------
    # Box with mustaches 
    plt.clf()
    
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    # >>> ALL
    
    
    pointsAvoidSumDF_plot = pd.DataFrame(pointsAvoidSumDF.groupby(['ID','level'], as_index = False).mean())    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF_plot ,  x="level", y = "avgDist", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot pointsAvoidSumAvgDF-avgDist.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF_plot ,  x="level", y = "maxLateralDist", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot pointsAvoidSumAvgDF-maxLateralDist.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
        
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF_plot ,  x="level", y = "yOnZero", order=constLevelList ,
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot pointsAvoidSumAvgDF-yOnZero.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    
    
    pointsExpectSumDF_plot = pd.DataFrame(pointsExpectSumDF.groupby(['ID','level'], as_index = False).mean())    
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsExpectSumDF_plot ,  x="level", y = "avgDist", order=constLevelList ,
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot pointsExpectSumDF-avgDist.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsExpectSumDF_plot ,  x="level", y = "maxLateralDist", order=constLevelList  , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot pointsExpectSumDF-maxLateralDist.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
        
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsExpectSumDF_plot ,  x="level", y = "yOnZero", order=constLevelList  , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot pointsExpectSumDF-yOnZero.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF_plot ,  x="level", y = "yOnZeroFromObst", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot pointsAvoidSumAvgDF-yOnZeroFromObst.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    #avgDistToObst
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF_plot ,  x="level", y = "avgDistToObst", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot pointsAvoidSumAvgDF-avgDistToObst.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data= pointsExpectSumDF_plot,  x="level", y = "yOnZeroFromObst", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot pointsExpectSumDF-yOnZeroFromObst.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    #avgDistToObst
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsExpectSumDF_plot ,  x="level", y = "avgDistToObst", order=constLevelList , palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot pointsExpectSumDF-avgDistToObst.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    # >>> Violin
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    
    ax = sns.violinplot(data=pointsAvoidSumDF_plot ,  x="level", y = "avgDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot(data=pointsAvoidSumDF_plot ,  x="level", y = "avgDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "vio pointsAvoidSumAvgDF-avgDist.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    
    ax = sns.violinplot(data=pointsAvoidSumDF_plot ,  x="level", y = "maxLateralDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot(data=pointsAvoidSumDF_plot ,  x="level", y = "maxLateralDist",  
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        

    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "vio pointsAvoidSumAvgDF-maxLateralDist.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
        
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    
    ax = sns.violinplot(data=pointsAvoidSumDF ,  x="level", y = "yOnZero",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot(data=pointsAvoidSumDF ,  x="level", y = "yOnZero",  
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
    
    
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "voi pointsAvoidSumAvgDF-yOnZero.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    
    
    ax = sns.violinplot(data=pointsExpectSumDF ,  x="level", y = "avgDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot(data=pointsExpectSumDF ,  x="level", y = "avgDist",  
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
    
    
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "vio pointsExpectSumDF-avgDist.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
            
    ax = sns.violinplot(data=pointsExpectSumDF ,  x="level", y = "maxLateralDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot(data=pointsExpectSumDF ,  x="level", y = "maxLateralDist",  
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
    
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "vio pointsExpectSumDF-maxLateralDist.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
        
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    ax = sns.violinplot(data=pointsExpectSumDF ,  x="level", y = "yOnZero", 
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        
    ax = sns.stripplot(data=pointsExpectSumDF ,  x="level", y = "yOnZero",  
                               jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
            
        
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "vio pointsExpectSumDF-yOnZero.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
    
    
    #<<< vio
    
    
    
    
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF ,  x="level", y = "area", order=constLevelList  , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot Area_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist", order=constLevelList  , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot avgDist_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
    #avgLateralDist      
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist", order=constLevelList  , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot avgLateralDist.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))  
        
    # yOnZero
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF ,  x="level", y = "yOnZero", order=constLevelList , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot Lateral_yOnZero_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    
    # <<< All
    
    for rep in repList:
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", order=constLevelList , 
                    palette= colorUserList )
        
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                           jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.75 , order=constLevelList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Area per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist", order=constLevelList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot avgDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
        
        # >>> avgLateralDist
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgLateralDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgLateralDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist"
                    , order=constLevelList , 
                                palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot avgLateralDist per level distribution_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
        # <<< avgLateralDist

    # <<< Visualisation 
    
    
    return 





# Calculate area of deviation and other metrics ONLY for the case of avoidence
def areaDistForAvoidenceSide(pathDF, fout):
        
    XborderMin = -120
    XborderMax = 120
    
    Xobst = 0
    Yobst = 0
    
    areaRate = [0,0]
    medRate = [0,0]
    
    nSteps = 5
    
    luserIdList= list(pathDF.ID.unique()) 
    
    areaAvoidDF = pd.DataFrame() 
    lateralAvoidDF = pd.DataFrame() 
    pointsAvoidDF = pd.DataFrame() 
    zeroXDeviationDF = pd.DataFrame() 
    
    pointsAvoidSumDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 
        
    
    for i, userID in enumerate(luserIdList):
        areaAvoid = np.zeros(startInd + len(constLevelList))
        areaAvoid[0] =  str(userID)
        
        
        pointsAvoid = np.zeros(startInd + len(constLevelList))
        pointsAvoid[0] = str(userID)
        
        lateralAvoid = np.zeros(startInd + len(constLevelList))
        lateralAvoid[0] = str(userID)
       
        zeroXDeviation = np.zeros(startInd + len(constLevelList))
        zeroXDeviation[0] =  str(userID)
          
        for rep in repList:
            
            areaAvoid[1] = rep 
            pointsAvoid[1] = rep 
            lateralAvoid[1] = rep
            zeroXDeviation[1] = rep 
            
            for l, level in enumerate(constLevelList):
        
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp']).reset_index()
                
                if not areaData.empty:
                    lvlID = constLevelList.index(level) 
                    
                    lyObst = np.array(areaData["OY"])   
                 
                    xSum = areaData['X'].sum() / len(areaData['X']) 
                    
                    ySum = areaData['Y'].sum() / len(areaData['Y']) 
                    
                    lx = np.array(areaData['X']) 
                    ly = np.array(areaData['Y'])
                    
                    
                    yOnZero = 0
                    
                    xzero = areaData['X'].abs().min(skipna = True, numeric_only = True)
                        
                    zeroIndArr = areaData.index[areaData['X'].abs() ==xzero]
                        
                    zeroInd = zeroIndArr[0]
                    if ((zeroInd> nSteps) and (zeroInd < len(ly)-nSteps )):
                        # check if the obstacle was on the same side of the scene
                        yOnZero = ly[zeroInd-nSteps:zeroInd+nSteps].sum() / len(ly[zeroInd-nSteps:zeroInd+nSteps])
                        yOnZeroObst = lyObst[zeroInd-nSteps:zeroInd+nSteps].sum() / len(lyObst[zeroInd-nSteps:zeroInd+nSteps])
                        
                        if yOnZero*yOnZeroObst < 0:
                            #NO obstacle avoidence
                            continue
                        
                         
                    avgDist = 0
                    avgLateralDist = 0
                    larea = 0
                    for ind in range(0, len(lx)) :
                        avgDist =avgDist+ np.sqrt( lx[ind]*lx[ind] +  ly[ind]*ly[ind] )
                        avgLateralDist =avgLateralDist+ np.sqrt( ly[ind]* ly[ind] )
                        if ind > 0:
                            larea = larea + abs(lx[ind]-lx[ind-1])*(abs(ly[ind-1])+abs(ly[ind]))/2
  
                    # Change integration
                    #area = np.trapz( np.array(areaData['X']) , np.array(areaData['Y'])  )
                    area = larea #np.trapz( np.array(areaData['Y']), np.array(areaData['X']) )
                    
                    if len(lx)>0:
                        pointsAvoid[ startInd +lvlID ] = avgDist / len(lx)
                        lateralAvoid[ startInd +lvlID ] = avgLateralDist / len(lx)
                        zeroXDeviation[ startInd +lvlID ] = yOnZero
                        avgDist = avgDist / len(lx)
                        avgLateralDist = avgLateralDist / len(lx)
                    else: 
                        pointsAvoid[ startInd +lvlID ] = 0
                        lateralAvoid[ startInd +lvlID ] = 0
                        zeroXDeviation[ startInd +lvlID ] = 0
                        avgDist = 0
                        avgLateralDist = 0
                        
                    area = np.abs(area)
                    areaAvoid[startInd+lvlID ] = area
                    
                    # for summary
                                 
                    
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'area': [area/10000],
                                            'avgDist': [avgDist],
                                            'avgLateralDist': [avgLateralDist],
                                            'yOnZero': [yOnZero],
                                            
                                            })
                    
                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)
                    
                    
                    
                
                         
            ls = [pointsAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            pointsAvoidDF = pd.concat([pointsAvoidDF, df ], ignore_index=True)
            
            
            ls = [areaAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles
            df['ID'] = userID
            df['rep'] = rep
            #areaAvoidDF =areaAvoidDF.append(df)
            areaAvoidDF = pd.concat([areaAvoidDF, df ], ignore_index=True)
            
                         
            ls = [lateralAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            lateralAvoidDF = pd.concat([lateralAvoidDF, df ], ignore_index=True)
            
            ls = [zeroXDeviation]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            zeroXDeviationDF = pd.concat([zeroXDeviationDF, df ], ignore_index=True)
            
        
    
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("int")  
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("string")
    
    foutStat = os.path.join(fout+"_Area_ForAvoidence.csv")   
    areaAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_Area_ForAvoidence.xlsx") 
    areaAvoidDF.to_excel(foutStat,  index = False, header=True) 
          

    #pointsAvoidDF['ID'] = pointsAvoidDF['ID'].astype("int")
    #pointsAvoidDF['ID'] = pointsAvoidDF['ID'].astype("string")
    foutStat = os.path.join(fout+"_pointsAvarageAvoid_ForAvoidence.csv")   
    pointsAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_pointsAvarageAvoid_ForAvoidence.xlsx") 
    pointsAvoidDF.to_excel(foutStat,  index = False, header=True) 
    
    foutStat = os.path.join(fout+"_lateralAvarageAvoid_ForAvoidence.csv")   
    lateralAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_lateralAvarageAvoid_ForAvoidence.xlsx") 
    lateralAvoidDF.to_excel(foutStat,  index = False, header=True) 
            
    foutStat = os.path.join(fout+"_lateralZeroXDeviationDF_ForAvoidence.csv")   
    zeroXDeviationDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_lateralZeroXDeviationDF_ForAvoidence.xlsx") 
    zeroXDeviationDF.to_excel(foutStat,  index = False, header=True) 
    
    # >>> Visualisation 
    
    fig, axs = plt.subplots(ncols=2, nrows=len(repList), figsize=(32, 32), dpi=cdpi)
                    
    
    # Create names on the x axis
    
    x_ticks = np.arange(len(constLevelList))
    
    for i, userID in enumerate(luserIdList):
        userColor = i/len(luserIdList)
        for rep in repList:
            
            for index, row in areaAvoidDF[ (areaAvoidDF["ID"] == userID) & (areaAvoidDF["rep"] == rep) ].iterrows():
                
                xarr = row.to_numpy()
                x = xarr[startInd:]
                print("------areaAvoidDF----->", x)
                
                axs[rep, 0].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)
                
                axs[rep, 0].set(xticks=x_ticks, xticklabels=constLevelList)          
                
                
            for index, row in pointsAvoidDF[ (pointsAvoidDF["ID"] == userID) & (pointsAvoidDF["rep"] == rep) ].iterrows():
                
                xarr = row.to_numpy()
                x = xarr[startInd:]
                print("------pointsAvoidDF----->", x)
                axs[rep, 1].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)
                # Create names on the x axis
                
                axs[rep, 1].set(xticks=x_ticks, xticklabels=constLevelList) 
                axs[rep, 1].legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12}); 
                
  

    
    
    axs[0,0].set_title('Avoidence Area for different levels_ForAvoidence')
    axs[0,1].set_title('Average distance to the center of the obstacle origine_ForAvoidence ')
    
    foutImg = os.path.join(fout + "Area per level distribution_ForAvoidence.pdf")    
    plt.savefig(f'{foutImg}', dpi=cdpi)
        
    
    # ------------------
    # Box with mustaches 
    plt.clf()
    
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    # >>> ALL
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF ,  x="level", y = "area", order=constLevelList , 
                palette= colorUserList )
    #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot Area_ForAvoidence_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    sns.boxplot(data=pointsAvoidSumDF, x="level", y = "avgDist", order=constLevelList , 
                palette= colorUserList )
    #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot avgDist_ForAvoidence_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
    #avgLateralDist      
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    sns.boxplot(data=pointsAvoidSumDF, x="level", y = "avgLateralDist", order=constLevelList  , 
                palette= colorUserList )
    #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot avgLateralDist_ForAvoidence.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))  
        
    # yOnZero
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ax = sns.boxplot(data=pointsAvoidSumDF ,  x="level", y = "yOnZero", order=constLevelList  , 
                palette= colorUserList )
    #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot Lateral_yOnZero_ForAvoidence_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    
    # <<< All
    
    for rep in repList:
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", order=constLevelList , 
                    palette= colorUserList )
        
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                           jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.75 , order=constLevelList )
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Area per level distribution_ForAvoidence_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist", order=constLevelList , 
                    palette= colorUserList )
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot avgDist per level distribution_ForAvoidence_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist",
                            inner='box', scale='width', dodge=False, linewidth=1)
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgDist per level distribution_ForAvoidence_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
        
        # >>> avgLateralDist
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgLateralDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgLateralDist per level distribution_ForAvoidence_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist", order=constLevelList , 
                    palette= colorUserList )
        
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot avgLateralDist per level distribution_ForAvoidence_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
        # <<< avgLateralDist
    
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "area",
                            inner='box', scale='width', dodge=True, linewidth=1, order=constLevelList )
        
        sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                      jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75, order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot area per level distribution_ForAvoidence_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
    
    # <<< Visualisation 
    
    
    return 





# Calculate MPD
def planDistFromObst(pathDF, fout):
        
    XborderMin = -120
    XborderMax = 120
    
    Xobst = 0
    Yobst = 0
    
    step = 20
    deltastep = 10 #int(step/2)
    
    luserIdList= list(pathDF.ID.unique()) 
    
    planDist = pd.DataFrame() 
     
    planDistDF = pd.DataFrame() 
    
    planDistToObstDF = pd.DataFrame() 
     
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 
    print(pathDF[ (pathDF["level"] == 'level108') ]["OY"])
    
    for i, userID in enumerate(luserIdList):
       
          
        for rep in repList:
            
            
            fig, axs = plt.subplots(ncols=1, nrows=len(constLevelList), figsize=(30, 30), dpi=200)
                        
            for l, level in enumerate(constLevelList):

                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                lvlID = constLevelList.index(level)
                  
                lx = np.array(areaData['X']) 
                ly = np.array(areaData['Y'])
                lt = np.array(areaData['timestamp'])
                    
                #lyObst = gaussian_filter1d( np.array(areaData["OY"]), 1)
                lyObst = np.array(areaData["OY"])
                    
                    
                planDistArr = np.zeros(len(ly))                    
                planDistToObstArr = np.zeros(len(ly))

                lateralDistToObst = 0
                lateralDist = 0
                    
                for ind in range(deltastep, len(ly)-deltastep-1 ) :
                        
                    y1=ly[ind-deltastep]
                    y2=ly[ind+deltastep]
                        
                    x1=lx[ind-deltastep]
                    x2=lx[ind+deltastep]
                        
                    planDistToObstArr[ind] =  y1 + (Xobst - x1)*(y2-y1)/(x2-x1) -Yobst
                        
                    planDistArr[ind] = y1 + (Xobst - x1)*(y2-y1)/(x2-x1)
                        

                
                axs[l].plot( lx, planDistArr , color= 'lightblue' , linewidth= (4.2 ), label=userID)
                
                axs[l].plot( lx, planDistToObstArr , color= 'blue' , linewidth= (4.2 ), label=userID)
                
                axs[l].plot( lx, lyObst , color= 'red' , linewidth= (4.2 ), label=userID)
                
                axs[l].plot( lx, ly , color= 'green' , linewidth= (4.2 ), label=userID)
                axs[l].set_title(level)
                
                
        
            foutImg = os.path.join(fout + "Planned path dist rep = "+ str(rep)+ " ID " + userID + ".pdf")    
            plt.savefig(f'{foutImg}', dpi=200)
                
            
            # ------------------
            # Box with mustaches 
            plt.clf()
                
                
                
        
    


# Calculate adea of deviation
def absement(pathDF, fout):
        
    XborderMin = -120
    XborderMax = 120
    
    Yobst = 0
    Yobst = 0
    
    
    luserIdList= list(pathDF.ID.unique()) 
    
    absementDF = pd.DataFrame() 
     
    absementToCenterDF = pd.DataFrame() 
     
    pointsAvoidSumDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 
    print(pathDF[ (pathDF["level"] == 'level108') ]["OY"])
    
    for i, userID in enumerate(luserIdList):
        absementArr = np.zeros(startInd + len(constLevelList))
        absementArr[0] =  str(userID)
        
        absementToCenterArr = np.zeros(startInd + len(constLevelList))
        absementToCenterArr[0] =  str(userID)
        
       
          
        for rep in repList:
            
            absementArr[1] = rep 
            absementToCenterArr[1] = rep 
            
            for l, level in enumerate(constLevelList):
        
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if not areaData.empty:
                    lvlID = constLevelList.index(level)
                  
                    lx = np.array(areaData['X']) 
                    ly = np.array(areaData['Y'])
                    lt = np.array(areaData['timestamp'])
                    
                    lyObst = gaussian_filter1d( np.array(areaData["OY"]), 1)
                         
                    lateralDispToObst = 0
                    lateralDisp = 0
                    for ind in range(0, len(ly)-1) :
                        dy1 = np.abs(lyObst[ind] - ly[ind])
                        
                        dy2 = np.abs(lyObst[ind+1] - ly[ind+1])
                        
                        lateralDispToObst = lateralDispToObst + (lt[ind+1]-lt[ind])*(np.abs(ly[ind]) + np.abs(ly[ind+1]) )/2
                        
                        lateralDisp = lateralDisp + (lt[ind+1]-lt[ind])*(dy1+dy2)/2
                        
                    
                    if len(ly)<=0:
                        lateralDispToObst = 0
                        lateralDisp = 0
                        
                        
                    absementArr[ startInd +lvlID ] = lateralDispToObst
                    absementToCenterArr[ startInd +lvlID ] = lateralDisp
                        
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'lateralDispToObst': lateralDispToObst,
                                            'lateralDispToCenter': lateralDisp,
                                            })
                    
                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)
                    
                    
                    
            # TODO : LOG pointsAvoidSumDF        
                    
            ls = [absementArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            absementDF = pd.concat([absementDF, df ], ignore_index=True)
            
                    
        
                    
            ls = [absementToCenterArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            absementToCenterDF = pd.concat([absementToCenterDF, df ], ignore_index=True)
            
            
    
    foutStat = os.path.join(fout+"_absement to obst.csv")   
    absementDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_absement to obst.xlsx") 
    absementDF.to_excel(foutStat,  index = False, header=True) 
                
    
    foutStat = os.path.join(fout+"_absement to center.csv")   
    absementDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_absement to center.xlsx") 
    absementDF.to_excel(foutStat,  index = False, header=True) 
        
    
    
    # ------------------
    # Box with mustaches 
    plt.clf()
    
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    
    
    # >>> Viz all repetitions
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        
    ax = sns.boxplot(data=pointsAvoidSumDF ,  x="level", y = "lateralDispToObst", order=constLevelList  , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot lateralDispToObst per level_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        
    ax = sns.boxplot(data=pointsAvoidSumDF ,  x="level", y = "lateralDispToCenter", order=constLevelList , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot lateralDispToCenter per level_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    plt.clf()
    
    # <<< Viz all repetitions
    
    
    
    for rep in repList:
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", 
                         y = "lateralDispToObst", order=constLevelList  , 
                                     palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "Abscement to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "lateralDispToObst",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  
                           x="level", y = "lateralDispToObst", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot Abscement to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
            
        # >>> to center
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  
                         x="level", y = "lateralDispToCenter", order=constLevelList , 
                                     palette= colorUserList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "Abscement to center_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "lateralDispToCenter",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "lateralDispToCenter", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot Abscement to center_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
    
    return 






# Calculate adea of deviation
def distToObst(pathDF, fout):
        
    XborderMin = -120
    XborderMax = 120
    
    Yobst = 0
    Yobst = 0
    
    areaRate = [0,0]
    medRate = [0,0]
    
    luserIdList= list(pathDF.ID.unique()) 
    
    areaAvoidDF = pd.DataFrame() 
    pointsAvoidDF = pd.DataFrame()
    lateralAvoidDF = pd.DataFrame() 
     
    pointsAvoidSumDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 
    print(pathDF[ (pathDF["level"] == 'level108') ]["OY"])
    
    for i, userID in enumerate(luserIdList):
        areaAvoid = np.zeros(startInd + len(constLevelList))
        areaAvoid[0] =  str(userID)
        
        
        pointsAvoid = np.zeros(startInd + len(constLevelList))
        pointsAvoid[0] = str(userID)
        
        lateralAvoid = np.zeros(startInd + len(constLevelList))
        lateralAvoid[0] = str(userID)
       
          
        for rep in repList:
            
            areaAvoid[1] = rep 
            pointsAvoid[1] = rep 
            lateralAvoid[1] = rep 
            
            for l, level in enumerate(constLevelList):
        
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if not areaData.empty:
                    lvlID = constLevelList.index(level)
                  
                    lx = np.array(areaData['X']) 
                    ly = np.array(areaData['Y'])
                    
                    #lyObst = np.array(areaData["OY"])
                    
                    lyObst = gaussian_filter1d( np.array(areaData["OY"]), 1)
                         
                    avgDist = 0
                    avgLateralDist = 0
                    for ind in range(0, len(ly)) :
                        avgDist =avgDist+ np.sqrt( (lyObst[ind] - ly[ind])*(lyObst[ind] - ly[ind]) +  lx[ind]*lx[ind] )
                        avgLateralDist =avgLateralDist+ np.sqrt( (lyObst[ind] - ly[ind])*(lyObst[ind] - ly[ind]) )
                        
                    # Change integration
                    #area = np.trapz( np.array(areaData['X']) , np.array(areaData['Y'])  )
                    area = np.trapz( np.array(areaData['Y'] - areaData['OY']), np.array(areaData['X']) )
                    
                    if len(ly)>0:
                        avgDist = avgDist / len(ly)
                        avgLateralDist = avgLateralDist/ len(ly)
                    else:
                        avgDist = 0
                        avgLateralDist = 0
                        
                    area = np.abs(area)
                    areaAvoid[startInd+lvlID ] = area
                    pointsAvoid[ startInd +lvlID ] = avgDist
                    lateralAvoid[ startInd +lvlID ] = avgLateralDist
                    
                    
                    # for summary
                    
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'area': [area/10000],
                                            'avgDist': [avgDist],
                                            'avgLateralDist': [avgLateralDist],
                                            })
                    
                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)
                    
                    
                    
                    
                    
            ls = [pointsAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            pointsAvoidDF = pd.concat([pointsAvoidDF, df ], ignore_index=True)
            
            
            ls = [areaAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles
            df['ID'] = userID
            df['rep'] = rep
            #areaAvoidDF =areaAvoidDF.append(df)
            areaAvoidDF = pd.concat([areaAvoidDF, df ], ignore_index=True)
        
                
            ls = [lateralAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            lateralAvoidDF = pd.concat([lateralAvoidDF, df ], ignore_index=True)
            
    
    

    areaAvoidAvgDF = pd.DataFrame(areaAvoidDF.groupby(['ID'], as_index = False).mean())    
    lateralAvoidAvgDF = pd.DataFrame(lateralAvoidDF.groupby(['ID'], as_index = False).mean())    
    pointsAvoidAvgDF = pd.DataFrame(pointsAvoidDF.groupby(['ID'], as_index = False).mean())  

    print(lateralAvoidAvgDF['rep'])

    foutStat = os.path.join(fout+"_obst_areaAvgDF.xlsx") 
    areaAvoidAvgDF.to_excel(foutStat,  index = False, header=True) 


    foutStat = os.path.join(fout+"_obst_lateralAvgDF.xlsx") 
    lateralAvoidAvgDF.to_excel(foutStat,  index = False, header=True) 


    foutStat = os.path.join(fout+"_obst_pointsAvgDF.xlsx") 
    pointsAvoidAvgDF.to_excel(foutStat,  index = False, header=True) 

    
    
    
    
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("int")  
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("string")
    
    foutStat = os.path.join(fout+"_obst_Area.csv")   
    areaAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_obst_Area.xlsx") 
    areaAvoidDF.to_excel(foutStat,  index = False, header=True) 
          

    foutStat = os.path.join(fout+"_obst_pointsAvarageAvoid.csv")   
    pointsAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_obst_pointsAvarageAvoid.xlsx") 
    pointsAvoidDF.to_excel(foutStat,  index = False, header=True) 
      

    foutStat = os.path.join(fout+"_obst_lateralAvarageAvoid.csv")   
    lateralAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_obst_lateralAvarageAvoid.xlsx") 
    lateralAvoidDF.to_excel(foutStat,  index = False, header=True) 
    
    
    print("Medium dist deviation rate for ", fout)
    print(medRate)
    print("Area deviation rate for ", fout)
    print(areaRate)
            
    # vis area
    
    
    
    fig, axs = plt.subplots(ncols=2, nrows=len(repList), figsize=(32, 32), dpi=cdpi)
                    
    
    # Create names on the x axis
    
    x_ticks = np.arange(len(constLevelList))
    
    for i, userID in enumerate(luserIdList):
        userColor = i/len(luserIdList)
        for rep in repList:
            
            for index, row in areaAvoidDF[ (areaAvoidDF["ID"] == userID) & (areaAvoidDF["rep"] == rep) ].iterrows():
                
                xarr = row.to_numpy()
                x = xarr[startInd:]
                
                axs[rep, 0].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)
                
                axs[rep, 0].set(xticks=x_ticks, xticklabels=constLevelList)          
                
                
            for index, row in pointsAvoidDF[ (pointsAvoidDF["ID"] == userID) & (pointsAvoidDF["rep"] == rep) ].iterrows():
                
                xarr = row.to_numpy()
                x = xarr[startInd:]
                axs[rep, 1].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)
                # Create names on the x axis
                
                axs[rep, 1].set(xticks=x_ticks, xticklabels=constLevelList) 
                axs[rep, 1].legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12}); 
                
  

    
    
    axs[0,0].set_title('Avoidence Area for different levels')
    axs[0,1].set_title('Average distance to the center of the obstacle origine ')
    
    foutImg = os.path.join(fout + "Area per level distribution.pdf")    
    plt.savefig(f'{foutImg}', dpi=cdpi)
        
    
    # ------------------
    # Box with mustaches 
    plt.clf()
    
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    # >>> all
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
        
    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "area", order=constLevelList , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot Area to obst_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF, x="level", y = "avgDist", order=constLevelList , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot avgDist to obst_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))
         
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist", order=constLevelList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot lateral to obst_All.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
    plt.clf()
    
    # <<< all
    
    print(pointsAvoidSumDF)
    
    for rep in repList:
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", order=constLevelList )
        
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                           jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.75 , order=constLevelList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Area to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist", order=constLevelList , 
                    palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot avgDist to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist"
                    , order=constLevelList  , palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot lateral to obstn_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  
                           x="level", y = "avgDist", jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgDist to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
    
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "area",
                            inner='box', scale='width', dodge=True, linewidth=1, order=constLevelList )
        
        sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                      jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot area to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
           
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgLateralDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgLateralDist to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
    
    return 





# Calculate adea of deviation
def distToObstAvoidersOnly(pathDF, fout):
        
    XborderMin = -120
    XborderMax = 120
    
    Yobst = 0
    Yobst = 0
    
    areaRate = [0,0]
    medRate = [0,0]
    
    luserIdList= list(pathDF.ID.unique()) 
    
    areaAvoidDF = pd.DataFrame() 
    pointsAvoidDF = pd.DataFrame()
    lateralAvoidDF = pd.DataFrame() 
     
    pointsAvoidSumDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 
    print(pathDF[ (pathDF["level"] == 'level108') ]["OY"])
    
    for i, userID in enumerate(luserIdList):
        areaAvoid = np.zeros(startInd + len(constLevelList))
        areaAvoid[0] =  str(userID)
        
        
        pointsAvoid = np.zeros(startInd + len(constLevelList))
        pointsAvoid[0] = str(userID)
        
        lateralAvoid = np.zeros(startInd + len(constLevelList))
        lateralAvoid[0] = str(userID)
       
          
        for rep in repList:
            
            areaAvoid[1] = rep 
            pointsAvoid[1] = rep 
            lateralAvoid[1] = rep 
            
            for l, level in enumerate(constLevelList):
        
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if not areaData.empty:
                    lvlID = constLevelList.index(level)
                  
                    lx = np.array(areaData['X']) 
                    ly = np.array(areaData['Y'])
                    
                    #lyObst = np.array(areaData["OY"])
                    
                    lyObst = gaussian_filter1d( np.array(areaData["OY"]), 1)
                         
                    avgDist = 0
                    avgLateralDist = 0
                    for ind in range(0, len(ly)) :
                        avgDist =avgDist+ np.sqrt( (lyObst[ind] - ly[ind])*(lyObst[ind] - ly[ind]) +  lx[ind]*lx[ind] )
                        avgLateralDist =avgLateralDist+ np.sqrt( (lyObst[ind] - ly[ind])*(lyObst[ind] - ly[ind]) )
                        
                    # Change integration
                    #area = np.trapz( np.array(areaData['X']) , np.array(areaData['Y'])  )
                    area = np.trapz( np.array(areaData['Y'] - areaData['OY']), np.array(areaData['X']) )
                    
                    if len(ly)>0:
                        avgDist = avgDist / len(ly)
                        avgLateralDist = avgLateralDist/ len(ly)
                    else:
                        avgDist = 0
                        avgLateralDist = 0
                        
                    area = np.abs(area)
                    areaAvoid[startInd+lvlID ] = area
                    pointsAvoid[ startInd +lvlID ] = avgDist
                    lateralAvoid[ startInd +lvlID ] = avgLateralDist
                    
                    
                    # for summary
                    
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'area': [area/10000],
                                            'avgDist': [avgDist],
                                            'avgLateralDist': [avgLateralDist],
                                            })
                    
                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)
                    
                    
                    
                    
                    
            ls = [pointsAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            pointsAvoidDF = pd.concat([pointsAvoidDF, df ], ignore_index=True)
            
            
            ls = [areaAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles
            df['ID'] = userID
            df['rep'] = rep
            #areaAvoidDF =areaAvoidDF.append(df)
            areaAvoidDF = pd.concat([areaAvoidDF, df ], ignore_index=True)
        
                
            ls = [lateralAvoid]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            lateralAvoidDF = pd.concat([lateralAvoidDF, df ], ignore_index=True)
            
    
    

    areaAvoidAvgDF = pd.DataFrame(areaAvoidDF.groupby(['ID'], as_index = False).mean())    
    lateralAvoidAvgDF = pd.DataFrame(lateralAvoidDF.groupby(['ID'], as_index = False).mean())    
    pointsAvoidAvgDF = pd.DataFrame(pointsAvoidDF.groupby(['ID'], as_index = False).mean())  

    print(lateralAvoidAvgDF['rep'])

    foutStat = os.path.join(fout+"_obst_areaAvgDF.xlsx") 
    areaAvoidAvgDF.to_excel(foutStat,  index = False, header=True) 


    foutStat = os.path.join(fout+"_obst_lateralAvgDF.xlsx") 
    lateralAvoidAvgDF.to_excel(foutStat,  index = False, header=True) 


    foutStat = os.path.join(fout+"_obst_pointsAvgDF.xlsx") 
    pointsAvoidAvgDF.to_excel(foutStat,  index = False, header=True) 

    
    
    
    
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("int")  
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("string")
    
    foutStat = os.path.join(fout+"_obst_Area.csv")   
    areaAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_obst_Area.xlsx") 
    areaAvoidDF.to_excel(foutStat,  index = False, header=True) 
          

    foutStat = os.path.join(fout+"_obst_pointsAvarageAvoid.csv")   
    pointsAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_obst_pointsAvarageAvoid.xlsx") 
    pointsAvoidDF.to_excel(foutStat,  index = False, header=True) 
      

    foutStat = os.path.join(fout+"_obst_lateralAvarageAvoid.csv")   
    lateralAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_obst_lateralAvarageAvoid.xlsx") 
    lateralAvoidDF.to_excel(foutStat,  index = False, header=True) 
    
    
    print("Medium dist deviation rate for ", fout)
    print(medRate)
    print("Area deviation rate for ", fout)
    print(areaRate)
            
    # vis area
    
    
    
    fig, axs = plt.subplots(ncols=2, nrows=len(repList), figsize=(32, 32), dpi=cdpi)
                    
    
    # Create names on the x axis
    
    x_ticks = np.arange(len(constLevelList))
    
    for i, userID in enumerate(luserIdList):
        userColor = i/len(luserIdList)
        for rep in repList:
            
            for index, row in areaAvoidDF[ (areaAvoidDF["ID"] == userID) & (areaAvoidDF["rep"] == rep) ].iterrows():
                
                xarr = row.to_numpy()
                x = xarr[startInd:]
                
                axs[rep, 0].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)
                
                axs[rep, 0].set(xticks=x_ticks, xticklabels=constLevelList)          
                
                
            for index, row in pointsAvoidDF[ (pointsAvoidDF["ID"] == userID) & (pointsAvoidDF["rep"] == rep) ].iterrows():
                
                xarr = row.to_numpy()
                x = xarr[startInd:]
                axs[rep, 1].plot( x, color=colorUserList[i], linewidth= (4.2 ), label=userID)
                # Create names on the x axis
                
                axs[rep, 1].set(xticks=x_ticks, xticklabels=constLevelList) 
                axs[rep, 1].legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12}); 
                
  

    
    
    axs[0,0].set_title('Avoidence Area for different levels')
    axs[0,1].set_title('Average distance to the center of the obstacle origine ')
    
    foutImg = os.path.join(fout + "Area per level distribution.pdf")    
    plt.savefig(f'{foutImg}', dpi=cdpi)
        
    
    # ------------------
    # Box with mustaches 
    plt.clf()
    
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
    
    # >>> all
    
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
        
    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "area", order=constLevelList  , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot Area to obst_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
    
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF, x="level", y = "avgDist", order=constLevelList  , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot avgDist to obst_ALL.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))
         
    plt.clf()
    plt.figure(figsize=(40,15))
        
    fig, ax = plt.subplots()
        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist", order=constLevelList  , 
                palette= colorUserList )
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
    plt.xticks(fontsize=7, rotation=45)
         
    foutImg = os.path.join(fout + "BoxPlot lateral to obst_All.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
    plt.clf()
    
    # <<< all
    
    print(pointsAvoidSumDF)
    
    for rep in repList:
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", order=constLevelList  , 
                    palette= colorUserList )
        
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                           jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.75 , order=constLevelList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Area to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=4*cdpi)
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist", order=constLevelList , 
                    palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot avgDist to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist"
                    , order=constLevelList , palette= colorUserList )
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)
        
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot lateral to obstn_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
        
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  
                           x="level", y = "avgDist", jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgDist to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
    
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "area",
                            inner='box', scale='width', dodge=True, linewidth=1, order=constLevelList )
        
        sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", 
                      jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot area to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
           
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralDist",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgLateralDist", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )
        
        ymin, ymax = plt.ylim()
        plt.ylim(ymin - 0.5*ymin, 1.4*ymax )
        plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
        plt.xticks(fontsize=7, rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgLateralDist to obst_rep" + str(rep) + ".pdf")    
        plt.savefig(f'{foutImg}', dpi=(4*cdpi))
    
    
    return 






# Calculate min distance 
def distCalculation(pathDF, fout):
        
    YborderMin = YEs+borderDist
    YborderMax = YNe-borderDist
    XborderMin = XlimMin
    
    timeRate = [0,0]
    
    levelList = list(pathDF.level.unique()) # uniquw=e levels list
    luserIdList= list(pathDF.ID.unique()) 
    
    distESDF = pd.DataFrame() 
    distNEDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in levelList:
        titles.append(l)
    
    for i, userID in enumerate(luserIdList):
        minDistlArrES = np.zeros(startInd + len(levelList))
        minDistlArrES[0] =  str(userID)
        
        minDistlArrNE = np.zeros(startInd + len(levelList))
        minDistlArrNE[0] =  str(userID)
       
        repList= list(pathDF.rep.unique()) 
        
            
        for rep in repList:
            
            minDistlArrES[1] = rep 
            minDistlArrNE[1] = rep 
            
            for l, level in enumerate(levelList):
        
                
                levelColor = l/len(levelList)
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)]
                pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df[ (df["X"]>(XborderMin)) & ( df["Y"]> YborderMin) & (df["Y"]< YborderMax) ].copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if level in levelList:
                    lvlID = levelList.index(level)
                    xArr =  np.array(areaData["X"])
                    yArr =  np.array(areaData["Y"])
                    
                    dNE = 0
                    dES = 0
                    if len(xArr)>0:
                        for i in range(0,len(xArr)):
                            ldNE = np.sqrt((XNe-xArr[i])*(XNe-xArr[i]) +(YNe-yArr[i])*(YNe-yArr[i]))
                            ldES = np.sqrt((XEs-xArr[i])*(XEs-xArr[i]) +(YEs-yArr[i])*(YEs-yArr[i]))
                            if (dNE == 0) or (dNE>ldNE):
                                dNE = ldNE
                            if (dES == 0) or (dES>ldES):
                                dES = ldES
                            
                            
                    minDistlArrES[ startInd +lvlID ] = dES
                    minDistlArrNE[ startInd +lvlID ] = dNE
                    
                    
            ls = [minDistlArrES]
            df = pd.DataFrame(ls)
            #distESDF = distESDF.append(df)
            
            distESDF = pd.concat([distESDF, df ], ignore_index=True) 
            
            ls = [minDistlArrNE]
            df = pd.DataFrame(ls)
            #distNEDF = distNEDF.append(df)
            
            distESDF = pd.concat([distESDF, df ], ignore_index=True) 
            
            
            #  if 3>5 -> true or if 6>4 -> true
            if abs(minDistlArrNE[-1])>abs(minDistlArrES[-3]) and (minDistlArrNE[-1]*minDistlArrES[-3] != 0) : #lvl 6 vs 4
                timeRate[0] += 1
            else:
                timeRate[1] +=1
                
            if abs(minDistlArrNE[-4])>abs(minDistlArrES[-2]) and (minDistlArrNE[-4]*minDistlArrES[-2] != 0):
                timeRate[0] +=1
            else:
                timeRate[1] +=1
                
    print("min dist to VR characters:", fout, timeRate)
    
    distESDF.columns = titles
    foutStat = os.path.join(fout+"_minDist_ES.csv")   
    distESDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_minDist_ES.xlsx") 
    distESDF.to_excel(foutStat,  index = False, header=True) 
          
    distNEDF.columns = titles
    foutStat = os.path.join(fout+"_minDist_NE.csv")   
    distNEDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_minDist_NE.xlsx") 
    distNEDF.to_excel(foutStat,  index = False, header=True) 
          
    return 

def repNum(filename):
    # destinguish repetitions
    rep = 404
    if ( "rep0" in filename ):
        rep = 0
    elif( "rep1" in filename ):
        rep = 1
    elif( "rep2" in filename ):
        rep = 2
    elif( "rep3" in filename ):
        rep = 3
    elif( "rep4" in filename ):
        rep = 4
    elif( "rep5" in filename ):
        rep = 5
    elif( "rep6" in filename ):
        rep = 6
    elif( "rep7" in filename ):
        rep = 7
    elif( "rep8" in filename ):
        rep = 8
    elif( "rep9" in filename ):
        rep = 9
        
    return rep


# Lateral deviation
def lateralDeviation(pathDF, fout):
        
    XborderMin = -200
    XborderMax = 200
    
    Xobst = 0
    Yobst = 0

    windowSize = 20
    luserIdList= list(pathDF.ID.unique()) 
    
    complTimeDF = pd.DataFrame()  
    
    pointsAvoidSumDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in constLevelList:
        titles.append(l)
    
    
    repList= list(pathDF.rep.unique()) 
        
    
    for i, userID in enumerate(luserIdList):
        complTime = np.zeros(startInd + len(constLevelList))
        complTime[0] =  str(userID)
        
        
       
          
        for rep in repList:
            
            complTime[1] = rep  
            
            for l, level in enumerate(constLevelList):
        
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if not areaData.empty:
                    lvlID = constLevelList.index(level)
                 
                    dt = areaData['timestamp'].max() -areaData['timestamp'].min() 
                    
                    
                    complTime[ startInd +lvlID ] = dt
                          
                    
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'dt': [dt],
                                            
                                            })
                    
                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)
                    
                    



plt.style.use('seaborn-whitegrid')


dir = os.getcwd()

path_of_the_directory= os.path.join(dir, "trackingData")

raw_data_csv_directory = os.path.join(path_of_the_directory, "rawData")

obst_data_csv_directory = os.path.join(path_of_the_directory, "rawDataObstCoord")

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
repIDs = []
repList = []

for filename in os.listdir(raw_data_csv_directory): 
    ids = re.findall(r"[-+]?\d*\.*\d+", filename)
    id = ids[0]
    extractId.append(id)
    
for filename in os.listdir(raw_data_csv_directory): 
    ids = re.findall(r"[-+]?\d*\.*\d+", filename)
    rep = ids[1]
    repIDs.append(rep)
    
iserIdList = set(extractId)
repList = set(repIDs)

FPSTarget = 40
cdpi = 40
nSlopes = 2
fc = 1.0
format_data = "%Y.%m.%d-%H.%M.%S.%f"
#format_data = "%Y.%m.%d-%H.%M.%S"
#YlimMax = 90 #110
#YlimMin = -110 #-130
#XlimMax = -30
#XlimMin = -250



Yshift = 10
Xshift = 350
    

borderDist = 0 #distance to compute area

XEs = 0 # ES character
YEs = 0 
XNe = 0 # Neurotic character
YNe = 0

YlimMax = 400 #110
YlimMin = -400 #-130
XlimMax = 290
XlimMin = -290

#colorLevelList = ['aquamarine2','brown1', 'burlywood2', 'cadetblue2', 'darkgoldenrod1', 'darkolivegreen2', 'darkorange','darkorchid1','darkturquoise', 'hotpink','indianred1', 'lightpink1', 'lightsalmon1','magenta2', 'maroon1','mediumorchid1','violetred1' ]

#colorLevelListRed = ['darkorchid1','darkorchid1','darkorchid1','darkturquoise', 'hotpink','indianred1', 'lightpink1', 'lightsalmon1','magenta2', 'maroon1','mediumorchid1','violetred1','palevioletred1', 'orchid1' ]
#colorLevelListBlue = ['lightblue','lightblue','lightblue','lightblue', 'darkturquoise','magenta2', 'orchid1','mediumorchid1','violetred1','palevioletred1' ]
#redList = ['orange1', 'orangered1', 'orchid1']


#colorLevelListRed = ['darkorchid','darkorchid','darkorchid','darkturquoise', 'hotpink','indianred', 'lightpink', 'lightsalmon','magenta', 'maroon','mediumorchid','darkorchid','deeppink', 'orchid' ]
colorLevelListRed = ['grey', 'darkgrey', 'royalblue','aquamarine','darkorchid','darkturquoise', 'gold','mediumvioletred', 'lightpink', 'deeppink','magenta', 'mediumorchid','darkorchid','deeppink', 'violet' ]
"""
colorUserList = ['grey', 'darkgrey', 'royalblue','aquamarine','darkorchid','darkturquoise', 'gold',
                  'mediumvioletred', 'lightpink', 'deeppink','magenta', 'mediumorchid',
                  'darkorchid','deeppink', 'violet' ,'brown', 
                  'burlywood', 'cadetblue', 'darkgoldenrod', 'darkolivegreen', 
                  'darkorange','darkturquoise', 'hotpink','indianred', 
                  'lightsalmon', 'maroon','violetred', 'salmon', 'salmon', 'salmon', 'salmon', 'salmon']
"""

colorUserList = ['grey', 'darkgrey', 'royalblue',
                 'darkturquoise','green', 'olive', 'gold',
                 'bisque', 'coral', 'deeppink', 
                 'crimson', 'brown',  'plum' ,
                 'magenta',
                  'burlywood', 'cadetblue', 'darkgoldenrod', 'darkolivegreen', 
                  'darkorange','darkturquoise', 'hotpink','indianred', 
                  'lightsalmon', 'maroon','violetred', 'salmon', 'salmon', 'salmon', 'salmon', 'salmon']


colorLevelListBlue = ['grey', 'darkgrey', 'lightblue', 
                      'royalblue', 'aquamarine','aqua', 
                      'cornflowerblue' ,'darkturquoise',
                      'khaki', 'gold', 'lightsalmon', 'coral',
                      'palevioletred', 'maroon', 'sepia', 
                      'grey', 'grey', 'grey' ]
redList = ['orange', 'orangered', 'orchid']

"""
constLevelList = ['level100', 'level110',  'vel107_2', 'vel108_2', 'level107', 'level108', 'level101',
                  'level102', 'level103', 'level104', 'level105', 'level106']
"""

constLevelList = ['level100', 'level110',  'vel107_2',
                  'level101', 'level104', 'level102', 
                  'level103','vel108_2', 'level107',
                  'level105', 'level106', 'level108' ]


constLevelTitles = ['Empty', 'Static', 'Per, Low Freq', 
                    'Slow, μ=1.34,σ=-0.23',   'Slow, μ=1.34,σ=0.39',       'Med, μ=0.90,σ=-0.60', 
                    'Med, μ=1.18,σ=0.55',     'Per, Low Freq, with noise', 'Per,  High Freq',    
                    'Fast, μ=1.18,σ=-1.34', 'Fast, μ=-0.71,σ=-1.5',      'Per, High Freq, with noise']

#constLevelList = ['level0_1', 'level0_2', 'level0_3', 'o_level1', 'o_level2', 'o_level3', 'o_level4', 'o_level5', 'o_level6']

#stat_data
statDF = pd.DataFrame(columns=['ID', 'level', 'rep', 'x', 'y'])
resTSstatDF = pd.DataFrame(columns=['ID', 'level', 'rep', 'x', 'y'])
resXYstatDF = pd.DataFrame(columns=['ID', 'level', 'rep', 'x', 'y'])

#merged DF per experiment
mergedDF = pd.DataFrame()
mergedTSDF = pd.DataFrame()
mergedXYDF = pd.DataFrame()
mergedRotDF = pd.DataFrame()
mergedRotTSDF = pd.DataFrame()

trackingDataLog = pd.DataFrame()

velocityDF  = pd.DataFrame()

plt.figure(figsize=(20,25), dpi= cdpi)

for id in iserIdList:
    
    
    
    fig, axs = plt.subplots(ncols=2, nrows=len(repList), figsize=(32, 64), dpi=cdpi)
    
    #fig = plt.figure(figsize=(80, 80), dpi=cdpi)
    
    fig.set_figheight(40)
    fig.set_figwidth(40)
    
    levelColor = 0
    titles = []
    
    #fig, axs = plt.subplots(2, 3, figsize=(40, 40), dpi = my_dpi)
    imgCount = 0
    for filename in os.listdir(raw_data_csv_directory): 
        
        name, ext = os.path.splitext(filename)
        
        rep = repNum(filename)
            
        #axs[rep,0] = plt.subplot2grid(shape=(1, 1), loc=(rep, 0), colspan=1, rowspan=1 , title = "raw")
        #axs[rep,1] = plt.subplot2grid(shape=(1, 1), loc=(rep, 1), colspan=1, rowspan=1, title = "Butterworth filter")
        
         
        if (( id in filename ) and (ext == '.csv') and ("Coord" in filename) ) :
            
            lvlID = constLevelList.index(filename[-18:-10])
            lbl = filename[-18:-10]
            
            Xarr = []
            Yarr = []
            Zarr = []
            
            OXarr = []
            OYarr = []
            OZarr = []
            
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
            tableObst = pd.DataFrame()
            
            lineAngle = ""    
            f = os.path.join(raw_data_csv_directory, filename)
            #fn = name+ "_ObstCoord.csv"
            #fobst = os.path.join(obst_data_csv_directory, fn)
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
                        startIndex = x[0]+1 #skip zeros
                        
                
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
                
                # tableObst        
                # Obstacle coordinates
                fileStrObst = filename[0:-9] + "ObstCoord.csv"
                findRes = find(fileStrObst, obst_data_csv_directory)
                
                if len(findRes) > 0:
                    fileObst = findRes[0]
                    csvfileObst = open(fileObst) 
                        
                    readerObst = csv.reader(csvfileObst)
                    
                    listObst = ( list(readerObst))
                    
                    tableObst = pd.DataFrame(listObst)
                    
    
    
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
                           
                    
                    OX = 0
                    OY = 0
                    OZ = 0
                    if tableObst.size>0:
                        s = tableObst.iloc[i][0]
                        coords = re.findall(r"[-+]?\d*\.*\d+", s)
                        OY=float(coords[0]) # Fix of coordinate backup in UE
                        OX=float(coords[1])
                        OZ=float(coords[2])
                           
                    #if not((X< XlimMin) or (X> XlimMax) or (Y> YlimMax)  or (Y< YlimMin)):
                    Xarr.append(X)
                    Yarr.append(Y)
                    Zarr.append(Z)
                    
                    YAarr.append(YA)
                    Parr.append(P)
                    Rarr.append(R)
                        
                    Tarr.append(T)
                    DTarr.append(dateTime)
                    
                    # obstacle
                    OYarr.append(OY)
            
            
            trackDF = pd.DataFrame({'X': Xarr,
                                    'Y': Yarr,
                                    'Z': Zarr,
                                    'YA': YAarr,
                                    'P': Parr,
                                    'R': Rarr,
                                    'OY': OYarr,
                                    'T': Tarr,
                                    'DT': DTarr,
                                    })
            
            # <<< Merge data in one DF
            
            
            
            
            
            trackBackup = trackDF.copy()
            
            axs[rep, 0].plot(np.array(trackBackup["X"]),np.array(trackBackup["Y"]), color=colorLevelListBlue[lvlID], linewidth= (6.2 ), label= lbl ) 
                
            
            
            # >>> Restore path of HDMI
            
            if (trackDF[ (trackDF["X"]<100) ][ (trackDF["X"]>10) ].empty):
                
                df = pd.DataFrame({'filename': [filename],
                                   'indStart':[0],
                                   'indMed': [0],
                                   'indEnd':[0],
                                   'len': [0],
                                   'lenBorder': [len(trackDF[ (trackDF["X"]<XlimMax) ][ (trackDF["X"]>XlimMin) ]) ]
                                         })
                
                print("--level empty allert--", filename)
                continue
            
            if (len(trackDF[ (trackDF["X"]<100) ][ (trackDF["X"]>10) ]) < 10 ):
                
                df = pd.DataFrame({'filename': [filename],
                                   'indStart':[0],
                                   'indMed': [0],
                                   'indEnd':[0],
                                   'len': [0],
                                   'lenBorder': [len(trackDF[ (trackDF["X"]<XlimMax) ][ (trackDF["X"]>XlimMin) ]) ]
                                         })
                print("--level empty allert--", filename)
                trackingDataLog =  pd.concat([trackingDataLog, df ])
                continue
            
            indMed = int(trackDF[ (trackDF["X"]<100) ][ (trackDF["X"]>10) ].index[-1])
            indStart = 0
            indEnd = len(trackDF)
            
            if len(trackDF[ (trackDF["X"]<100) ][ (trackDF["X"]>10) ] ):
                for i, r in trackDF.iloc[indMed:].iterrows():
                    X = r.X
                    Y = r.Y
                    if ((X< XlimMin) or (X> XlimMax) or (Y> YlimMax)  or (Y< YlimMin)):
                        indEnd = i
                        print("--end index allert--", indStart, indEnd)
                        break
                    
                for i,r in trackDF.iloc[:indMed][::-1].iterrows():
                    X = r.X
                    Y = r.Y
                    if ((X< XlimMin) or (X> XlimMax) or (Y> YlimMax)  or (Y< YlimMin)):
                        indStart = i
                        print("--begin index allert--", indStart, indEnd)
                        break
            
            trackDF = trackDF[indStart:indEnd]
            #fout = os.path.join(modified_csv_directory, filename)
            #trackDF.to_csv(fout)
            # <<< Restore path of HDMI
            
            
            
            
            trackRotDF = trackDF.copy()
            
            
            if (( np.array(trackDF['X'])[-1] - np.array(trackDF['X'])[0] ) <0) :
                
                # Rotate
                                
                Yobst = []
                theta = np.radians(180)
                c, s = np.cos(theta), np.sin(theta)
                rot_matrix = np.array(((c, s), (-s, c)))
                Xtemp = np.array(trackRotDF['X'])
                Ytemp = np.array(trackRotDF['Y'])
                OYtemp = np.array(trackRotDF['OY'])
                XY = np.array([Xtemp, Ytemp]).T @ rot_matrix
                if sum(OYtemp) != 0.0:
                    XYobst = np.array([Xtemp, OYtemp]).T @ rot_matrix
                    Yobst = XYobst[:, 1]
                    
                else:
                    Yobst = np.zeros( len(Xtemp) )
                
                trackRotDF['X'] = XY[:, 0]
                trackRotDF['Y'] = XY[:, 1]
                trackRotDF['OY'] = Yobst
               
            
            # >>> rotate
            """
            if ( np.array(trackDF['X'])[-1] - np.array(trackDF['X'])[0] ) :
                trackRotDF = trackRotDF[indStart:indEnd]
                trackRotDF['X'] = np.array(trackRotDF['X'])[indEnd:indStart]
                trackRotDF['Y'] = np.array(trackRotDF['Y'])[indEnd:indStart]
                trackRotDF['OY'] = np.array(trackRotDF['OY'])[indEnd:indStart]
            
            
                
                Yobst = []
                theta = np.radians(180)
                c, s = np.cos(theta), np.sin(theta)
                rot_matrix = np.array(((c, s), (-s, c)))
                XY = np.array([Xarr, Yarr]).T @ rot_matrix
                if sum(OYarr) != 0.0:
                    XYobst = np.array([Xarr, OYarr]).T @ rot_matrix
                    Yobst = XYobst[:, 1]
                    
                else:
                    Yobst = OYarr
                
                
                
                trackRotDF = pd.DataFrame({'X': XY[:, 0],
                                        'Y': XY[:, 1],
                                        'Z': Zarr,
                                        'YA': YAarr,
                                        'P': Parr,
                                        'R': Rarr,
                                        'OY': Yobst,
                                        'T': Tarr,
                                        'DT': DTarr,
                                        })
                
            
            """
            # <<< rotate    
            
            
            
            df = pd.DataFrame({'filename': [filename],
                               'indStart':[indStart],
                               'indMed': [indMed],
                               'indEnd':[indEnd],
                               'len': [len(trackDF[indStart:indEnd])],
                               'lenBorder': [len(trackDF[ (trackDF["X"]<XlimMax) ][ (trackDF["X"]>XlimMin) ]) ]
                                     
                                     })
            
            trackingDataLog =  pd.concat([trackingDataLog, df ])
                
            
            TDarr = np.array(trackDF['DT'])
            TimeSarr = []
            
            for i, time_data in enumerate(TDarr):
            
                date = datetime.strptime(time_data, format_data)
                timeVars = re.findall(r'\d+', time_data)
                #ts = date.timestamp() 
                ts = date.timestamp() - date.microsecond/1000000 + int(timeVars[-1]) / 1000
                
                TimeSarr.append(ts)
            
            trackDF['timestamp'] = TimeSarr
            trackDF = trackDF.sort_values(by=['timestamp'])
            
            fout = os.path.join(modified_csv_directory, filename)
            trackDF.to_csv(fout)    
            
            df = trackDF.copy()
            df["ID"] = id
            df['level']= filename[-18:-10]
            df['rep' ]=  rep
            #mergedDF = mergedDF.append(df)
            
            mergedDF = pd.concat([mergedDF, df ], ignore_index=True)
            
            
            rotdf = df.copy()
            rotdf['X']= trackRotDF['X'] 
            rotdf['Y']= trackRotDF['Y'] 
            rotdf['OY']= trackRotDF['OY'] 
            mergedRotDF = pd.concat([mergedRotDF, rotdf ], ignore_index=True)
            
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
                  
                
                statDF = pd.concat([statDF, pd.DataFrame({
                                        'ID': [id], 
                                        'level': [filename[-18:-10]],
                                        'rep' : [rep],
                                        'x': [breaks[1]] , 
                                        'y': [pwlf_model.predict(breaks[1])[0]]
                                        
                                        })])
                """
                statDF = pd.concat([statDF, pd.DataFrame({
                                        'ID': [ id ], 
                                        'level': [filename[-18:-10]],
                                        'x': [xInter], 
                                        'y': [yInter],
                                        
                                        })], ignore_index=True)
                """
                
                #plt.plot(Xarr, Yarr, '.', color=[0.1*randColor,1 -0.2*levelColor, 0.2*levelColor], linewidth=0.02 , label=filename[-18:-10])
                                #,figure=plt.figure());            
                randColor = random.randint(0,10)
                #for xb in breaks:
                #    pltInit = plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor,1 -0.2*levelColor, 0.2*levelColor], marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                    
                
                #plt.plot(XEs, YEs, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
             #markersize=20, markevery=12)
                
                foutImg = os.path.join(img_directory,name+"_resample.pdf")    
                #plt.savefig(f'{foutImg}', dpi=my_dpi)
                
            
            # >>> Time resample to a new FPS
            # Linear Approximation for filtered data
            lbl = "Tsample " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                
                (TResArr ,x ,y, z ) = resample3D_time(
                    np.array(trackDF["timestamp"]), 
                    np.array(trackDF["X"]), 
                    np.array(trackDF["Y"]),
                    np.array(trackDF["Z"]), FPSTarget)
                
                (TResArr ,ya ,p, r ) = resample3D_time(
                    np.array(trackDF["timestamp"]), 
                    np.array(trackDF["YA"]), 
                    np.array(trackDF["P"]),
                    np.array(trackDF["R"]), FPSTarget)
                
                (TResArr ,oy ,p, r ) = resample3D_time(
                        np.array(trackDF["timestamp"]), 
                        np.array(trackDF["OY"]), 
                        np.array(trackDF["P"]),
                        np.array(trackDF["R"]), FPSTarget)
                
                    
                resTSDF = pd.DataFrame({'X': x,
                                        'Y': y,
                                        'Z': z,
                                        'YA': ya,
                                        'P': p,
                                        'R': r,
                                        'timestamp': TResArr,
                                        'OY': oy
                                        })
                fout = os.path.join(modified_csv_filt_dir, filename+"TSresample.csv")
                resTSDF.to_csv(fout)
                
                df = resTSDF.copy()
                df["ID"] = id
                df['level']= filename[-18:-10]
                df['rep' ]=  rep
                mergedTSDF = pd.concat([mergedTSDF, df])
                
                (TResRotArr ,xRot ,yRot, zRot ) = resample3D_time(
                    np.array(rotdf["timestamp"]), 
                    np.array(rotdf["X"]), 
                    np.array(rotdf["Y"]),
                    np.array(rotdf["Z"]), FPSTarget)
                
                
                rotdf = df.copy()
                
                rotdf["timestamp"] = TResRotArr
                rotdf["X"] = xRot
                rotdf["Y"] = yRot
                
                mergedRotTSDF = pd.concat([mergedRotTSDF, rotdf])
                
                
                pwlf_model = pwlf.PiecewiseLinFit(x, y ,degree=1)
                breaks = pwlf_model.fit(nSlopes)
                
                slopes = pwlf_model.calc_slopes()
                
                resTSstatDF = pd.concat([resTSstatDF, pd.DataFrame({
                                        'ID': [id], 
                                        'level': [filename[-18:-10]],
                                        'rep' : [rep],
                                        'x': [breaks[1]] , 
                                        'y': [pwlf_model.predict(breaks[1])[0]]
                                        
                                        })])
                
                """
                resTSstatDF = resTSstatDF.append({
                                        'ID': id, 
                                        'level': filename[-18:-10],
                                        'rep' : rep,
                                        'x': breaks[1] , 
                                        'y': pwlf_model.predict(breaks[1])[0]
                                        
                                        }, ignore_index=True)
                """                
                randColor = random.randint(0,10)
                #for xb in breaks:
                    
                    #pltResTS = plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                
                #plt.plot(x, pwlf_model.predict(x), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= lbl)
                
                # TO visualize: 
                    
                #if rep == 0:
                #axs[rep, 0].plot(np.array(trackDF["X"]),np.array(trackDF["Y"]), color=colorLevelListBlue[lvlID], linewidth= (6.2 ), label= lbl ) 
                #else:
                #    axs[rep, 0].plot(np.array(trackDF["X"]),np.array(trackDF["Y"]), color=colorLevelListBlue[lvlID], linewidth= (6.2 )) 
                    
                x = np.array(trackRotDF['X'])
                
                y = np.array(trackRotDF['Y'])
                
                axs[rep, 1].plot(x,y, color=colorLevelListBlue[lvlID], linewidth= (6.2 ),label= lbl )       
                
                axs[rep, 0].set_title("Raw, rep = " + str(rep))
                axs[rep, 1].set_title("Rotated, rep = "+ str(rep))
                
                axs[rep, 0].plot( XEs, (YEs), color='limegreen',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
                         markersize=18, markevery=12)
                axs[rep, 1].plot( XEs, (YEs), color='limegreen',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
                         markersize=18, markevery=12)
                axs[rep, 1].legend( loc='lower left' , ncol = 2, frameon=True, prop={'size': 10});
                
                axs[rep, 0].legend( loc='lower left' , ncol = 2, frameon=True, prop={'size': 10});
                
                #titles.append(filename[-18:-10])
                # <<< Time resample to a new FPS
                
            
            # >>> Space resample 
            # Linear Approximation for filtered data
            lbl = "TXsample " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                
                (t ,x ,y, z, ya, p, r, oy ) = resample_spatial_all(
                    np.array(trackDF["timestamp"]), 
                    np.array(trackDF["X"]), 
                    np.array(trackDF["Y"]),
                    np.array(trackDF["Z"]),
                    np.array(trackDF["YA"]), 
                    np.array(trackDF["P"]),
                    np.array(trackDF["R"]), 
                    np.array(trackDF["OY"])
                    )
                
                    
                resXYDF = pd.DataFrame({'X': x,
                                        'Y': y,
                                        'Z': z,
                                        'YA': ya,
                                        'P': p,
                                        'R': r,
                                        'timestamp': t,
                                        'OY': oy
                                        })
                
                fout = os.path.join(modified_csv_filt_dir, filename+"XYresample.csv")
                resXYDF.to_csv(fout)
                
                
                df = resXYDF.copy()
                df["ID"] = id
                df['level']= filename[-18:-10]
                df['rep' ]=  rep
                #mergedXYDF = mergedXYDF.append(df)
                mergedXYDF =  pd.concat([mergedXYDF, df])

                pwlf_model = pwlf.PiecewiseLinFit(x, y ,degree=1)
                breaks = pwlf_model.fit(nSlopes)
                
                slopes = pwlf_model.calc_slopes() 
                
                resXYstatDF = pd.concat([resXYstatDF, pd.DataFrame({
                                                        'ID': [id], 
                                                        'level': [filename[-18:-10]],
                                                        'rep' : [rep],
                                                        'x': [breaks[1]] , 
                                                        'y': [pwlf_model.predict(breaks[1])[0]]
                                                        
                                                        })])
                """
                resXYstatDF = resXYstatDF.append({
                                        'ID': id, 
                                        'level': filename[-18:-10],
                                        'rep' : rep,
                                        'x': breaks[1] , 
                                        'y': pwlf_model.predict(breaks[1])[0]
                                        
                                        }, ignore_index=True)
                """
                randColor = random.randint(0,10)
                #for xb in breaks:
                    
                    #pltResTSXT = plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                    #plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                
                
                #pltResTSXT = plt.plot(x, pwlf_model.predict(x), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= lbl)
                #pltResTSXT = plt.plot(x,y, '.', color=[0.524,0.2*levelColor,1 - 0.2*levelColor], linewidth=0.18 , label="filtered "+filename[-18:-10])       
                #plt.plot(x, pwlf_model.predict(x), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor])
                #plt.plot(x,y, '.', color=[0.524,0.2*levelColor,1 - 0.2*levelColor], linewidth=0.18 , label= lbl)       
                #titles.append(filename[-18:-10])
                # <<< Time resample to a new FPS
                
            
            
            levelColor =+ 0.05
            
          
    
    foutImg = os.path.join(img_directory,id+ ".pdf")      
    
    
    #fig.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
    fig.savefig(f'{foutImg}', dpi=4*cdpi)


    
foutImg = os.path.join(img_directory,"TOTAL"+"_lin"+str(nSlopes)+"_lowPass_" +str(fc) + ".pdf")      
    
fig.savefig(f'{foutImg}', dpi=4*cdpi)

foutStat = os.path.join(stat_directory,"trackingDataLog"+".xlsx")   
trackingDataLog.to_excel(foutStat)
    
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+".csv")   
statDF.to_csv(foutStat)
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+".xlsx") 
statDF.to_excel(foutStat,  index = False, header=True) 


"""
# for resampled
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resampleTS"+str(FPSTarget)+".csv")   
resTSstatDF.to_csv(foutStat)
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resampleTS"+str(FPSTarget)+".xlsx") 
resTSstatDF.to_excel(foutStat,  index = False, header=True) 


foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resampleXY"+str(FPSTarget)+".csv")   
resXYstatDF.to_csv(foutStat)
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resampleXY"+str(FPSTarget)+".xlsx") 
resXYstatDF.to_excel(foutStat,  index = False, header=True) 
"""


"""
# statistics  vis
fout = os.path.join(stat_directory,"deviationStat"+"_lin"+str(nSlopes)+"_raw")   
correctRate = visStatRepSeperate(statDF, fout)
print( "Raw: ", correctRate )

fout = os.path.join(stat_directory,"deviationStat"+"_lin"+str(nSlopes)+"_resTS"+str(FPSTarget))   
correctRate = visStatRepSeperate(resTSstatDF, fout)
print("TS: ", correctRate )

fout = os.path.join(stat_directory,"deviationStat"+"_lin"+str(nSlopes)+"_resXY"+str(FPSTarget))   
correctRate = visStatRepSeperate(resXYstatDF, fout)    
print("XY: ", correctRate )
"""

#area



fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
areaDist(mergedTSDF, fout)

"""
fout = os.path.join(stat_directory,"_resXY"+str(FPSTarget))  
areaDist(mergedXYDF, fout)  
"""
flagStop = 1

#if (flagStop == 1) :
#    quit() 

fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
trajPlot(mergedRotTSDF, fout)


fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
velocityPlot(mergedTSDF, fout)
      

fout = os.path.join(stat_directory,"_Rot_resTS"+str(FPSTarget))   
velocityDF = velocityStat(mergedRotTSDF, fout)

"""
fout = os.path.join(stat_directory,"_AvoidersSide_resTS"+str(FPSTarget))   
areaDistForAvoidenceSide(mergedTSDF, fout)
"""

fout = os.path.join(stat_directory,"_resTS_AvoidersOnly"+str(FPSTarget))   
areaDistForAvoidersOnly(mergedTSDF, fout, velocityDF)


  

"""
# for resampled
fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
areaCalculation(mergedTSDF, fout)

fout = os.path.join(stat_directory,"_resXY"+str(FPSTarget))   
areaCalculation(mergedXYDF, fout)
"""

#distance to obstacle


fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
distToObst(mergedTSDF, fout)


#time
"""
fout = os.path.join(stat_directory,"_raw_")   
completionTime(mergedDF, fout)
"""
fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
completionTime(mergedTSDF, fout)
"""
fout = os.path.join(stat_directory,"_resXY"+str(FPSTarget))  
completionTime(mergedXYDF, fout)  

#min max
fout = os.path.join(stat_directory,"_raw_")   
statisticsMinMax(mergedDF, fout)
"""
fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
statisticsMinMaxDist(mergedTSDF, fout)
"""
fout = os.path.join(stat_directory,"_resXY"+str(FPSTarget))  
statisticsMinMax(mergedXYDF, fout) 

# absement
fout = os.path.join(stat_directory,"_raw_")   
absement(mergedDF, fout)
"""
fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
absement(mergedTSDF, fout)
"""
fout = os.path.join(stat_directory,"_resXY"+str(FPSTarget))  
absement(mergedXYDF, fout) 

# Calculate dist of deviation - MPD
fout = os.path.join(stat_directory,"_raw_")   
planDistFromObst(mergedDF, fout)
"""
fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
planDistFromObst(mergedTSDF, fout)
"""
fout = os.path.join(stat_directory,"_resXY"+str(FPSTarget))  
planDistFromObst(mergedXYDF, fout) 
"""





fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
velocityPlot(mergedTSDF, fout)
      

fout = os.path.join(stat_directory,"_Rot_resTS"+str(FPSTarget))   
velocityDF = velocityStat(mergedRotTSDF, fout)


# >>> viz trajectory

idDF = mergedDF[mergedDF["ID"] == '6601']
repDF = idDF[idDF["rep"] == 0]
             
lvlDFraw = repDF[repDF["level"] == 'level107']
                          

idDF = mergedTSDF[mergedTSDF["ID"] == '6601']
repDF = idDF[idDF["rep"] == 0]
             
lvlDFTS = repDF[repDF["level"] == 'level107']
                          

idDF = mergedXYDF[mergedXYDF["ID"] == '6601']
repDF = idDF[idDF["rep"] == 0]
             
lvlDFXY = repDF[repDF["level"] == 'level107']
                          
    
plt.clf()
plt.figure(figsize=(40,40), dpi= 200)
#plt.scatter(lvlDFXY, s=20, c="green", alpha=0.5)

t1 = np.array(lvlDFXY['timestamp'])
x1 = np.array(lvlDFXY['X']) 
y1 = np.array(lvlDFXY['Y'])
x2 = np.array(lvlDFTS['X']) 
y2 = np.array(lvlDFTS['Y'])
x3 = np.array(lvlDFraw['X']) 
y3 = np.array(lvlDFraw['Y'])

                         

plt.plot(x1, y1 , color=[0.1*randColor ,0.2*levelColor,0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', alpha = 0.5, markeredgewidth=2, markersize=4, markevery=2)
plt.plot(x2, y2 , color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', alpha = 0.5, markeredgewidth=2, markersize=4, markevery=2)
plt.plot(x3, y3 , color=[0.1*randColor ,1 - 0.1*levelColor,1 - 0.1*levelColor] ,marker='h', markerfacecolor='lightgreen', alpha = 0.5, markeredgewidth=2, markersize=4, markevery=2)


foutImg = os.path.join(fout + "_ResExample1.pdf")     
plt.savefig(foutImg, dpi=200)

                      
    
plt.clf()
plt.figure(figsize=(40,40), dpi= 200)                         

plt.plot(x1, y1 , color=[0.1*randColor ,0.2*levelColor,0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', alpha = 0.5, markeredgewidth=2, markersize=4, markevery=2, label="XY")

plt.plot(x2, y2-10 , color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', alpha = 0.5, markeredgewidth=2, markersize=4, markevery=2, label="TS")

plt.plot(x3, y3-20 , color=[0.1*randColor ,1 - 0.1*levelColor,1 - 0.1*levelColor] ,marker='h', markerfacecolor='lightgreen', alpha = 0.5, markeredgewidth=2, markersize=4, markevery=2, label="Raw")


plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12});    

foutImg = os.path.join(fout + "_ResExample2.pdf")     
plt.savefig(foutImg, dpi=200)


"""



# Calculate task completion Time
def velocityStat(pathDF, fout):

    import matplotlib.cm as cm
    XborderMin = -300
    XborderMax = 300 

    Xobst = 0
    Yobst = 0

    XdecisionConst = -254

    areaRate = [0,0]
    medRate = [0,0]

    luserIdList= list(pathDF.ID.unique()) 

    complTimeDF = pd.DataFrame()  

    pointsAvoidSumDF = pd.DataFrame() 


    expTDF = pd.DataFrame() 
    expXDF = pd.DataFrame() 
    expYDF = pd.DataFrame() 

    pointsAvoidersDF = pd.DataFrame()

    nSteps = 20

    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)

    for l in constLevelList:
        titles.append(l)



    repList= list(pathDF.rep.unique()) 

    for i, userID in enumerate(luserIdList):

        expTimeArr = np.zeros(startInd + len(constLevelList))
        expXCoordArr = np.zeros(startInd + len(constLevelList))
        expYCoordArr = np.zeros(startInd + len(constLevelList))

        fig, axs = plt.subplots(ncols=3, nrows=len(repList), figsize=(64, 48), dpi=2*cdpi)

        for rep in repList:

            expTimeArr[0] =  str(userID)
            expTimeArr[1] = rep 
            expXCoordArr[0] =  str(userID)
            expXCoordArr[1] = rep 
            expYCoordArr[0] =  str(userID)
            expYCoordArr[1] = rep 

            for l, level in enumerate(constLevelList):

                lvlID = constLevelList.index(level)

                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()

                print(userID , rep, level, df.size)
                if df.size <= 0 or ( '100' in level ):
                    continue

                areaData1 = df.copy()

                areaData = areaData1.sort_values(by=['timestamp'])

                t = np.array(areaData['timestamp'])

                tmin = t[0]

                tnorm = t - tmin

                x = np.array(areaData['X'])

                y = np.array(areaData['Y'])

                v = []


                for i, it in enumerate(t):

                    i0 = 0
                    i1 = 0

                    if ( (i> nSteps) and (i < len(t)-nSteps) ):
                        i0 = i-nSteps
                        i1 = i+nSteps
                        v.append( np.sqrt( (x[i1]-x[i0])*(x[i1]-x[i0]) + (y[i1]-y[i0])*(y[i1]-y[i0]) ) / np.abs(t[i1] - t[i0])  )
                    else:
                        v.append(0)


                varr = np.array(v)
                if len(varr) >0:

                    vnorm = 1/ np.max(v)

                    vmax = np.max(v)
                    vmaxInd = np.argmax(v)

                    v25perc = 0.25*vmax

                    v25percInd = 0

                    for ind, vval in enumerate(varr[vmaxInd:0:-1]):
                        if vval < v25perc:
                            v25percInd = vmaxInd - ind

                            break



                    expTimeArr[ startInd +lvlID ] = tnorm[v25percInd]
                    expXCoordArr[ startInd +lvlID ] = x[v25percInd]
                    expYCoordArr[ startInd +lvlID ] = y[v25percInd]

                    colors = cm.rainbow(vnorm * varr)

                    axs[rep,0].scatter(x,y, c=colors, linewidth= (6.2 ) )   
                    axs[rep,0].set_title("Raw, rep = " + str(rep))

                    axs[rep,0].plot( x, y, color= colorLevelListRed[l], alpha = 0.5, linewidth= 0.2 ,marker='h', markerfacecolor=colorLevelListRed[l], markeredgewidth=12,
                             markersize=18, markevery=12, label= constLevelTitles[l] )


                    axs[rep,0].scatter(x[v25percInd],y[v25percInd], c="black", linewidth= (10.2 ) )  
                    axs[rep,0].scatter(x[vmaxInd],y[vmaxInd], c="black", linewidth= (10.2 ) )   

                    #axs[rep].legend( fontsize=40, loc='lower left' , ncol = 2, frameon=True, prop={'size': 10});
                    axs[rep,0].legend(  loc='lower left' , ncol = 2 );


                    axs[rep,1].set_title("Velocity in time")


                    print(tnorm)
                    print(v)

                    axs[rep,1].plot( tnorm, v, color= colorLevelListRed[l], alpha = 0.8, linewidth= 4.2, label= constLevelTitles[l] )
                    axs[rep,1].legend(  loc='lower left' , ncol = 2 );

                    axs[rep,1].scatter(tnorm[v25percInd],v[v25percInd], color= "red", linewidth= (10.2 ) )   

                    axs[rep,1].scatter(tnorm[vmaxInd],v[vmaxInd], color= "red", linewidth= (10.2 ) ) 


                    axs[rep,2].set_title("Velocity in space")


                    axs[rep,2].plot( x, v, color= colorLevelListRed[l], alpha = 0.8, linewidth= 4.2, label= constLevelTitles[l] )
                    axs[rep,2].legend(  loc='lower left' , ncol = 2 );


                    axs[rep,2].scatter(x[v25percInd],v[v25percInd], c="red", linewidth= (10.2 ) )  
                    axs[rep,2].scatter(x[vmaxInd],v[vmaxInd], c="red", linewidth= (10.2 ) )   



                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'expTime': [tnorm[v25percInd]],
                                            'expXCoord': [x[v25percInd]],
                                            'expYCoord': [y[v25percInd]],

                                            })


                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)

                    if x[v25percInd] < XdecisionConst:

                        pointsAvoidersDF =  pd.concat([pointsAvoidersDF, df ], ignore_index=True)



            print(titles) 

            print(expYCoordArr) 

            print(expXCoordArr) 

            print(expTimeArr) 

            ls = [expYCoordArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            expYDF = pd.concat([expYDF, df ], ignore_index=True)

            ls = [expXCoordArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            expXDF = pd.concat([expXDF, df ], ignore_index=True)


            ls = [expTimeArr]
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            expTDF = pd.concat([expTDF, df ], ignore_index=True)



        foutImg = os.path.join(fout + userID+"nSteps_" + str(nSteps) + "_velocityProfile.pdf")      


        fig.savefig(f'{foutImg}', dpi=2*cdpi)



    foutStat = os.path.join(fout+"_avoidersDF.csv")   
    pointsAvoidersDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_avoidersDF.xlsx") 
    pointsAvoidersDF.to_excel(foutStat,  index = False, header=True) 



    foutStat = os.path.join(fout+"_expYDF.csv")   
    expYDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_expYDF.xlsx") 
    expYDF.to_excel(foutStat,  index = False, header=True) 

    foutStat = os.path.join(fout+"_expXDF.csv")   
    expXDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_expXDF.xlsx") 
    expXDF.to_excel(foutStat,  index = False, header=True) 


    foutStat = os.path.join(fout+"_expTDF.csv")   
    expTDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_expTDF.xlsx") 
    expTDF.to_excel(foutStat,  index = False, header=True) 


    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "expTime", order=constLevelList )

        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 0.8*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Wait Time.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF,  x="level", y = "expTime",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF,  x="level", y = "expTime", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 0.8*ymax )

    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait Time.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))




    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "expXCoord", order=constLevelList )

        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ymin, ymax = plt.ylim()
    plt.ylim(ymin - 0.5*ymin, 0.8*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Wait X.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,40))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF,  x="level", y = "expXCoord",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF,  x="level", y = "expXCoord", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )


    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait X.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))





    plt.clf()
    plt.figure(figsize=(40,15))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.boxplot(data=pointsAvoidSumDF,  x="level", y = "expTime", order=constLevelList )

        # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    #ymin, ymax = plt.ylim()
    #plt.ylim(ymin - 0.5*ymin, 0.8*ymax )
    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)

    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "BoxPlot Wait Time.pdf")    
    plt.savefig(f'{foutImg}', dpi=4*cdpi)

    plt.clf()
    plt.figure(figsize=(40,40))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF, y = "expXCoord",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF, y = "expXCoord", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )


    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait X All.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))



    plt.clf()
    plt.figure(figsize=(40,40))

    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelTitles)), labels=constLevelTitles)

    ax = sns.violinplot(data=pointsAvoidSumDF, y = "expTime",
                            inner='box', scale='width', dodge=False, linewidth=1, order=constLevelList )
    ax = sns.stripplot( data=pointsAvoidSumDF, y = "expTime", 
                           jitter=True, linewidth=0.4, hue = 'ID', alpha = 0.75 , order=constLevelList )


    plt.legend(title='id', fontsize='6', title_fontsize='6', loc='upper left', ncol=3)
    plt.xticks(fontsize=7, rotation=45)

    foutImg = os.path.join(fout + "ViolinePlot Wait Time All.pdf")    
    plt.savefig(f'{foutImg}', dpi=(4*cdpi))

    return pointsAvoidersDF
"""

