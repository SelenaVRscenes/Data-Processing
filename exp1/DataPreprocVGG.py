# Same as TimeSpaceResampeMSdata.py but with visualisation optimisation
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
# Resample data to 70 points in trajectory

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

import sklearn.datasets
import sklearn.metrics
from sklearn.model_selection import train_test_split
import xgboost as xgb


from ray import tune

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



def resampleTime_3D(time, xdata, ydata, zdata, newFrameRate):
    
    #re s
    duration = time[-1] - time[0]
    nbSamples = (duration * newFrameRate) + 1
    #newTime = np.linspace(time[0], time[-1], len(time))
    newTime = np.linspace(time[0], time[-1], newFrameRate)
    
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
  


def resample_spatial_all_fixFPS(tdata, xdata, ydata, zdata, yadata, rdata, pdata):
    
    #re sampling on distance
    distArr = []
    cFPS = 70
    
    #distArr.append(np.sqrt(xdata[0]*xdata[0] + ydata[0]*ydata[0])  )
    distArr.append(0)
    
    for i, x in enumerate(xdata[1:]):
        dx = xdata[i+1] - xdata[i]
        dy = ydata[i+1] - ydata[i]
        
        dist = distArr[i] + np.sqrt(dx*dx + dy*dy)
        #dist = distArr[i] + dx + dy
        distArr.append(dist)
        
    newDist = np.linspace(distArr[0], distArr[-1],cFPS)
    
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


# Calculate adea of deviation
def areaCalculation(pathDF, fout):
        
    YborderMin = YEs+borderDist
    YborderMax = YNe-borderDist
    XborderMin = XEs
    
    areaRate = [0,0]
    medRate = [0,0]
    
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
        plt.figure(figsize=(12,10), dpi= 10)
        #fig, axs = plt.subplots(1, 1, figsize=(10, 10), dpi=my_dpi)
        for rep in repList:
            
            areaAvoid[1] = rep 
            pointsAvoid[1] = rep 
            
            for l, level in enumerate(levelList):
        
                
                levelColor = l/len(levelList)
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)]
                pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df[ (df["X"]>(XborderMin)) & ( df["Y"]> YborderMin) & (df["Y"]< YborderMax) ].copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if level in levelList:
                    lvlID = levelList.index(level)
                 
                    xSum = areaData['X'].sum() / len(areaData['X']) - XborderMin
                    
                    ySum = areaData['Y'].sum() / len(areaData['Y']) 
                    
                    lx = np.array(areaData['X']) / len(areaData['X']) - XEs 
                    
                    if ySum<0:
                        ly = np.array(areaData['Y']) / len(areaData['Y']) - YEs 
                         
                    else:
                        ly = np.array(areaData['Y']) / len(areaData['Y']) - YNe
  
                    avgCount = 0
                    for ind in range(0, len(lx)) :
                        avgCount =avgCount+ np.sqrt( lx[ind]*lx[ind] +  ly[ind]*ly[ind] )
  
                    area = np.trapz( np.array(areaData['X']  - XborderMin), np.array(areaData['Y'])  )
                    area2 = np.trapz(  np.array(areaData['Y'] - YborderMax ),  np.array(areaData['X']) )  
                    
                    pointsAvoid[ startInd +lvlID ] = avgCount
                    areaAvoid[startInd+lvlID ] = area
                    
                    #print(userID , level, ySum, area) #, area2)
                    
                    if (lvlID > 0) and ('level5' in level or 'level6' in level) and int(rep)==0:
                        
                        disp = pathPltDF[ (pathPltDF["X"]>-200) ]
                        
                        if 'level5' in level:
                        
                                
                            plt.plot( (Xshift+ np.array(disp['X'])), -(Yshift+ np.array(disp['Y'])) , color= 'skyblue', linewidth= (6.2 ), label="Paths")
                            
                            plt.plot( (Xshift+ np.array(areaData['X'])), -(Yshift+ np.array(areaData['Y'])) , color= 'mediumvioletred', linewidth= (6.8  ), label="Prox. Areas")
                           
                        else:
                                    
                            plt.plot( (Xshift+ np.array(disp['X'])), -(Yshift+ np.array(disp['Y']))  , color= 'skyblue', linewidth= (6.2 ))
                            
                            plt.plot( (Xshift+ np.array(areaData['X'])), -(Yshift+ np.array(areaData['Y'])) , color= 'mediumvioletred', linewidth= (6.8 ))
                           
                         
            ls = [pointsAvoid]
            df = pd.DataFrame(ls)
            pointsAvoidDF =pointsAvoidDF.append(df)
                    
            ls = [areaAvoid]
            df = pd.DataFrame(ls)
            areaAvoidDF =areaAvoidDF.append(df)
            
            
            #  if 3>5 -> true or if 6>4 -> true
            if abs(areaAvoid[-1])>abs(areaAvoid[-3]) and (areaAvoid[-1]*areaAvoid[-3]!=0):
                areaRate[0] += 1
            else:
                areaRate[1] +=1
                
            if abs(areaAvoid[-4])>abs(areaAvoid[-2])and (areaAvoid[-2]*areaAvoid[-4]!=0):
                areaRate[0] +=1
            else:
                areaRate[1] +=1
                
            
            #  if 3>5 -> true or if 6>4 -> true
            if abs(pointsAvoid[-1])>abs(pointsAvoid[-3]):
                medRate[0] += 1
            else:
                medRate[1] +=1
                
            if abs(pointsAvoid[-4])>abs(pointsAvoid[-2]):
                medRate[0] +=1
            else:
                medRate[1] +=1
            

                
        plt.plot( Xshift+XEs, -(Yshift+YEs), color='limegreen',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
             markersize=18, markevery=12, label="ES character")
        plt.plot( Xshift+XNe, -(Yshift+YNe), color='coral', marker='h', markerfacecolor='mistyrose', markeredgewidth=12,
             markersize=18, markevery=12, label="NE character")
          
        #lgnd  = plt.legend( loc='lower left' , ncol = 1, frameon=True, prop={'size': 16});    
        lines = plt.gca().get_lines()
        include = [0,1]
        legend1 = plt.legend([lines[i] for i in include],[lines[i].get_label() for i in include], loc=2, frameon=True, prop={'size': 16});   
        legend2 = plt.legend([lines[i] for i in [-2,-1]],['ES character','NE character'],loc=3, frameon=True, prop={'size': 20});   
        
        
        plt.gca().add_artist(legend1)
        
        plt.title('Proximity Area examples for conditions 8 and 9').set_size(20)
        plt.xlabel('X-position (cm)').set_size(20)
        plt.ylabel('Y-position (cm)').set_size(20)
        
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        
        
     
    # Create names on the x axis
    
        foutImg = os.path.join(fout + str(userID) + "_Area.png")    
        plt.savefig(f'{foutImg}', dpi=my_dpi)
        
    
    areaAvoidDF.columns = titles
    foutStat = os.path.join(fout+"_Area.csv")   
    areaAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_Area.xlsx") 
    areaAvoidDF.to_excel(foutStat,  index = False, header=True) 
          
    pointsAvoidDF.columns = titles
    foutStat = os.path.join(fout+"_pointsAvarageAvoid.csv")   
    pointsAvoidDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_pointsAvarageAvoid.xlsx") 
    pointsAvoidDF.to_excel(foutStat,  index = False, header=True) 
    
    
    print("Medium dist deviation rate for ", fout)
    print(medRate)
    print("Area deviation rate for ", fout)
    print(areaRate)

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


# Calculate time of task completion 
def timeComplCalculation(pathDF, fout):
        
    YborderMin = YEs+borderDist
    YborderMax = YNe-borderDist
    XborderMin = XlimMin
    
    timeRate = [0,0]
    
    levelList = list(pathDF.level.unique()) # uniquw=e levels list
    luserIdList= list(pathDF.ID.unique()) 
    
    timeComplDF = pd.DataFrame() 
    
    
    titles = ["ID", "rep" ]
              #"Turn to top", "YEs on top"]
    startInd =  len(titles)
    
    for l in levelList:
        titles.append(l)
    
    for i, userID in enumerate(luserIdList):
        timeComplArr = np.zeros(startInd + len(levelList))
        timeComplArr[0] =  str(userID)
        
       
        repList= list(pathDF.rep.unique()) 
        
            
        for rep in repList:
            
            timeComplArr[1] = rep 
            
            for l, level in enumerate(levelList):
        
                
                levelColor = l/len(levelList)
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)]
                pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df[ (df["X"]>(XborderMin)) & ( df["Y"]> YborderMin) & (df["Y"]< YborderMax) ].copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if level in levelList:
                    lvlID = levelList.index(level)
                    tArr =  np.array(areaData["timestamp"])
                    if len(tArr)>0:
                        timeCompl = max(tArr) - min(tArr)
                    else:
                        timeCompl=0
                    
                    timeComplArr[ startInd +lvlID ] = timeCompl
                    
            ls = [timeComplArr]
            df = pd.DataFrame(ls)
            timeComplDF =timeComplDF.append(df)
                    
            
            #  if 3>5 -> true or if 6>4 -> true
            if abs(timeComplArr[-1])>abs(timeComplArr[-3]) and (timeComplArr[-1]*timeComplArr[-3]!=0):
                timeRate[0] += 1
            else:
                timeRate[1] +=1
                
            if abs(timeComplArr[-4])>abs(timeComplArr[-2]) and (timeComplArr[-4]*timeComplArr[-2]!=0):
                timeRate[0] +=1
            else:
                timeRate[1] +=1
                
    print("Time to complete", fout, timeRate)
    
    timeComplDF.columns = titles
    foutStat = os.path.join(fout+"_timeCompl.csv")   
    timeComplDF.to_csv(foutStat)
    foutStat = os.path.join(fout+"_timeCompl.xlsx") 
    timeComplDF.to_excel(foutStat,  index = False, header=True) 
          
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
            distESDF = distESDF.append(df)
                    
            ls = [minDistlArrNE]
            df = pd.DataFrame(ls)
            distNEDF = distNEDF.append(df)
                    
            
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



plt.style.use('seaborn-whitegrid')

path_of_the_directory= os.path.join("C:\\Users\ypatotsk\Documents\Dev\pythonDev", "trackingData")
raw_data_csv_directory = os.path.join(path_of_the_directory, "rawData")
#path_of_the_directory= os.path.join("C:\\Dev\\pythonDev", "trackingData")

modified_csv_directory= os.path.join( path_of_the_directory, "csvCleanedData_exp2")

modified_csv_filt_dir = os.path.join( path_of_the_directory, "csvFilteredData_exp2")

img_directory= os.path.join( path_of_the_directory, "csvVisData_exp2")

stat_directory= os.path.join( path_of_the_directory, "statistics_exp2")
    
outputDir = os.path.join(path_of_the_directory, "img_exp2")

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
                    
                    
outputDir = os.path.join(path_of_the_directory, "imgSorted_exp2")
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
#YlimMax = 90 #110
#YlimMin = -110 #-130
#XlimMax = -30
#XlimMin = -250



Yshift = 10
Xshift = 350
    

borderDist = 25 #distance to compute area

XEs = -118 # ES character
YEs = -110 
XNe = -118 # Neurotic character
YNe = 90

YlimMax = YNe #110
YlimMin = YEs #-130
XlimMax = -30
XlimMin = -250

#colorLevelList = ['aquamarine2','brown1', 'burlywood2', 'cadetblue2', 'darkgoldenrod1', 'darkolivegreen2', 'darkorange','darkorchid1','darkturquoise', 'hotpink','indianred1', 'lightpink1', 'lightsalmon1','magenta2', 'maroon1','mediumorchid1','violetred1' ]

#colorLevelListRed = ['darkorchid1','darkorchid1','darkorchid1','darkturquoise', 'hotpink','indianred1', 'lightpink1', 'lightsalmon1','magenta2', 'maroon1','mediumorchid1','violetred1','palevioletred1', 'orchid1' ]
#colorLevelListBlue = ['lightblue','lightblue','lightblue','lightblue', 'darkturquoise','magenta2', 'orchid1','mediumorchid1','violetred1','palevioletred1' ]
#redList = ['orange1', 'orangered1', 'orchid1']


#colorLevelListRed = ['darkorchid','darkorchid','darkorchid','darkturquoise', 'hotpink','indianred', 'lightpink', 'lightsalmon','magenta', 'maroon','mediumorchid','darkorchid','deeppink', 'orchid' ]
colorLevelListRed = ['royalblue','aquamarine','darkorchid','darkturquoise', 'gold','mediumvioletred', 'lightpink', 'deeppink','magenta', 'mediumorchid','darkorchid','deeppink', 'violet' ]

colorLevelListBlue = ['lightblue','lightblue','lightblue','lightblue', 'darkturquoise','mediumslateblue', 'royalblue','aquamarine','blueviolet' ]
redList = ['orange', 'orangered', 'orchid']


constLevelList = ['level0_1', 'level0_2', 'level0_3', 'o_level1', 'o_level2', 'o_level3', 'o_level4', 'o_level5', 'o_level6']

#stat_data
statDF = pd.DataFrame(columns=['ID', 'level', 'x', 'y'])
resTSstatDF = pd.DataFrame(columns=['ID', 'level', 'x', 'y'])
resTSXYstatDF = pd.DataFrame(columns=['ID', 'level', 'x', 'y'])

#merged DF per experiment
mergedDF = pd.DataFrame()
mergedTSDF = pd.DataFrame()
mergedTSXYDF = pd.DataFrame()



plt.figure(figsize=(20,25), dpi= my_dpi)


dataSet = np.array([])
labelSet = np.array([])


for id in iserIdList:
    
    levelColor = 0
    titles = []
    
    #fig, axs = plt.subplots(2, 3, figsize=(40, 40), dpi = my_dpi)
    imgCount = 0
    for filename in os.listdir(raw_data_csv_directory): 
        
        name, ext = os.path.splitext(filename)
         
        if (( id in filename ) and (ext == '.csv') and ("Coord" in filename) and not( "init" in filename ) and filename[-11:-10] != '0' ) :
        
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
            
            
            # >>> Time resample to a new FPS
            # Linear Approximation for filtered data
            lbl = "Tsample " + filename[-18:-10]
            FPSTarget = 70
            if (len(Xarr)>0 and len(Yarr)>0):
                
                
                (TResArr ,x ,y, z ) = resampleTime_3D(
                    np.array(trackDF["timestamp"]), 
                    np.array(trackDF["X"]), 
                    np.array(trackDF["Y"]),
                    np.array(trackDF["Z"]), FPSTarget)
                
                (TResArr ,ya ,p, r ) = resampleTime_3D(
                    np.array(trackDF["timestamp"]), 
                    np.array(trackDF["YA"]), 
                    np.array(trackDF["P"]),
                    np.array(trackDF["R"]), FPSTarget)
                
                print("--time--")
                
                print(max(np.array(trackDF["timestamp"])), 
                    max(TResArr))
                
                    
                resTSDF = pd.DataFrame({'X': x,
                                        'Y': y,
                                        'Z': z,
                                        'YA': ya,
                                        'P': p,
                                        'R': r,
                                        'timestamp': TResArr-TResArr[0],
                                        })
                
                fout = os.path.join(modified_csv_filt_dir, filename+"TSresample.csv")
                resTSDF.to_csv(fout)
                
                df = resTSDF.copy()
                df["ID"] = id
                df['level']= filename[-18:-10]
                df['rep' ]=  rep
                mergedTSDF = mergedTSDF.append(df)
                
                lvlID = constLevelList.index(filename[-18:-10])
                
                
                
                randColor = random.randint(0,10)
                
                plt.plot(x,y, color=colorLevelListBlue[lvlID], linewidth= (6.2 ) )       
                
                arr = resTSDF.to_numpy()
                B = np.reshape(arr, arr.shape[0] * arr.shape[1] )
                
                if (filename[-11:-10] in ('1','2') ):
                    print(filename[-11:-10])
                    print(filename[-11:-10] in ('1','2','3','4','5','6'))
                
                    if dataSet.size ==0  : 
                        dataSet = np.array([B])
                    else:
                        dataSet = np.append(dataSet, [B], axis = 0)
                      
                    clab = 0
                    if (filename[-11:-10] in ('1') ):
                        clab= 1
                        
                        
                    if labelSet.size ==0 : 
                        labelSet = np.array(clab)
                    else:
                        labelSet = np.append(labelSet, clab)
                        
                
                
                
            
            levelColor =+ 0.05
            
    
                
    plt.plot( XEs, (YEs), color='limegreen',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
             markersize=18, markevery=12, label="ES character")
    plt.plot( XNe, (YNe), color='coral', marker='h', markerfacecolor='mistyrose', markeredgewidth=12,
             markersize=18, markevery=12, label="NE character")
          
    
    foutImg = os.path.join(img_directory,id+"_lin"+str(nSlopes)+"_lowPass_" +str(fc) + ".png")      
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
    
    if len(lgnd.legendHandles) >0:
        lgnd.legendHandles[0]._legmarker.set_markersize(40)
        
        for l in lgnd.legendHandles:
            l._legmarker.set_markersize(40)
        #lgnd.legendHandles[1]._legmarker.set_markersize(40)
    
    plt.savefig(f'{foutImg}', dpi=my_dpi)
 

# Split data
train_x, test_x, train_y, test_y = train_test_split(
        dataSet, labelSet, test_size=0.25)
# Build input matrices for XGBoost
train_set = xgb.DMatrix(train_x, label=train_y)
test_set = xgb.DMatrix(test_x, label=test_y)



# Train the classifier
results = {}




from keras.applications.vgg16 import VGG16
model = VGG16()

print(model.summary())
    






foutImg = os.path.join(img_directory,"TOTAL"+"_lin"+str(nSlopes)+"_lowPass_" +str(fc) + ".png")      
    
plt.savefig(f'{foutImg}', dpi=my_dpi)

    


foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resampleTS"+str(FPSTarget)+".csv")   
resTSstatDF.to_csv(foutStat)
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resampleTS"+str(FPSTarget)+".xlsx") 
resTSstatDF.to_excel(foutStat,  index = False, header=True) 



    












