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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime
import seaborn as sns

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
    plt.savefig(f'{foutImg}', dpi=cdpi)
           
    return devRate


# Calculate adea of deviation
def areaCalculation(pathDF, fout):
        
    XborderMin = -120
    XborderMax = 120
    
    Xobst = 0
    Yobst = 0
    
    areaRate = [0,0]
    medRate = [0,0]
    
    luserIdList= list(pathDF.ID.unique()) 
    
    areaAvoidDF = pd.DataFrame() 
    pointsAvoidDF = pd.DataFrame() 
    
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
       
          
        for rep in repList:
            
            areaAvoid[1] = rep 
            pointsAvoid[1] = rep 
            
            for l, level in enumerate(constLevelList):
        
                
                #curData = pathDF[ ([pathDF["ID"] == userID ][pathDF["rep"] == rep ][ pathDF["level"]  ==level])]
                df = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level) 
                            & (pathDF["X"] > XborderMin) & (pathDF["X"] < XborderMax) ].copy()
                
                #pathPltDF = pathDF[(pathDF["ID"] == userID) & (pathDF["rep"] == rep) & (pathDF["level"]  ==level)].copy()
                
                areaData1 = df.copy()
                
                areaData = areaData1.sort_values(by=['timestamp'])
                
                if level in constLevelList:
                    lvlID = constLevelList.index(level)
                 
                    xSum = areaData['X'].sum() / len(areaData['X']) 
                    
                    ySum = areaData['Y'].sum() / len(areaData['Y']) 
                    
                    lx = np.array(areaData['X']) 
                    ly = np.array(areaData['Y'])
                         
                    avgCount = 0
                    for ind in range(0, len(lx)) :
                        avgCount =avgCount+ np.sqrt( lx[ind]*lx[ind] +  ly[ind]*ly[ind] )
  
                    # Change integration
                    #area = np.trapz( np.array(areaData['X']) , np.array(areaData['Y'])  )
                    area = np.trapz( np.array(areaData['Y']), np.array(areaData['X']) )
                    
                    if len(lx)>0:
                        pointsAvoid[ startInd +lvlID ] = avgCount / len(lx)
                        avgCount = avgCount / len(lx)
                    else: 
                        pointsAvoid[ startInd +lvlID ] = 0
                        avgCount = 0
                        
                    area = np.abs(area)
                    areaAvoid[startInd+lvlID ] = area
                    
                    print(userID , level, ySum, area, avgCount) #, area2)
                    
                    # for summary
                                 
                    
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'area': [area/10000],
                                            'avgCount': [avgCount],
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
        
        
    
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("int")  
    #areaAvoidDF['ID'] = areaAvoidDF['ID'].astype("string")
    
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
    
    foutImg = os.path.join(fout + "Area per level distribution.png")    
    plt.savefig(f'{foutImg}', dpi=cdpi)
        
    
    # ------------------
    # Box with mustaches 
    plt.clf()
    
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)
    
    print(pointsAvoidSumDF)
    
    for rep in repList:
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area")
        
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.4 )
        #sns.stripplot(data=df, x='day', y='tip', color='black', alpha=0.5)
        
        plt.xticks(rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Area per level distribution_rep" + str(rep) + ".png")    
        plt.savefig(f'{foutImg}', dpi=10*cdpi)
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgCount")
        
        plt.xticks(rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot avgCount per level distribution_rep" + str(rep) + ".png")    
        plt.savefig(f'{foutImg}', dpi=(10*cdpi))
        
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgCount",
                            inner='box', scale='width', dodge=False, linewidth=1)
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgCount", jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.4 )
        
        plt.xticks(rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgCount per level distribution_rep" + str(rep) + ".png")    
        plt.savefig(f'{foutImg}', dpi=(10*cdpi))
    
    
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)

        sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "area",
                            inner='box', scale='width', dodge=True, linewidth=1)
        
        sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.4 )
        
        plt.xticks(rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot area per level distribution_rep" + str(rep) + ".png")    
        plt.savefig(f'{foutImg}', dpi=(10*cdpi))
    
    
    
    return 






# Calculate adea of deviation
def distCalculationToObst(pathDF, fout):
        
    XborderMin = -120
    XborderMax = 120
    
    Xobst = 0
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
    print(pathDF[ (pathDF["level"] == 'level108') ]["OX"])
    
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
                
                if level in constLevelList:
                    lvlID = constLevelList.index(level)
                  
                    lx = np.array(areaData['X']) 
                    ly = np.array(areaData['Y'])
                    
                    lxObst = np.array(areaData["OX"])
                         
                    avgCount = 0
                    avgLateralCount = 0
                    for ind in range(0, len(lx)) :
                        avgCount =avgCount+ np.sqrt( (lxObst[ind] - lx[ind])*(lxObst[ind] - lx[ind]) +  ly[ind]*ly[ind] )
                        avgLateralCount =avgLateralCount+ np.sqrt( (lxObst[ind] - lx[ind])*(lxObst[ind] - lx[ind]) )
  
                    # Change integration
                    #area = np.trapz( np.array(areaData['X']) , np.array(areaData['Y'])  )
                    area = np.trapz( np.array(areaData['Y']), np.array(areaData['X'] - areaData['OX']) )
                    
                    if len(lx)>0:
                        avgCount = avgCount / len(lx)
                        avgLateralCount = avgLateralCount/ len(lx)
                    else:
                        avgCount = 0
                        avgLateralCount = 0
                        
                    print( "pointsAvoid", pointsAvoid )
                    print( "avgCount", avgCount )
                    
                    area = np.abs(area)
                    areaAvoid[startInd+lvlID ] = area
                    pointsAvoid[ startInd +lvlID ] = avgCount
                    lateralAvoid[ startInd +lvlID ] = avgLateralCount
                    
                    print(userID , level, avgLateralCount, area, avgCount) #, area2)
                    
                    # for summary
                                 
                    
                    df = pd.DataFrame({'ID': [userID],
                                            'rep': [rep],
                                            'level': [level],
                                            'area': [area/10000],
                                            'avgCount': [avgCount],
                                            'avgLateralCount': [avgLateralCount],
                                            })
                    
                    pointsAvoidSumDF = pd.concat([pointsAvoidSumDF, df ], ignore_index=True)
                    
                    
                    
            print("lateralAvoid" , lateralAvoid)
            print("pointsAvoid" , pointsAvoid)
            print("areaAvoid" , areaAvoid)
                         
            ls = [pointsAvoid]
            print("pointsAvoid -> ", ls)
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            #pointsAvoidDF =pointsAvoidDF.append(df)
            pointsAvoidDF = pd.concat([pointsAvoidDF, df ], ignore_index=True)
            
            
            ls = [areaAvoid]
            print("areaAvoid -> ", ls)
            df = pd.DataFrame(ls)
            df.columns = titles
            df['ID'] = userID
            df['rep'] = rep
            #areaAvoidDF =areaAvoidDF.append(df)
            areaAvoidDF = pd.concat([areaAvoidDF, df ], ignore_index=True)
        
                
            ls = [lateralAvoid]
            print("lateralAvoid -> ", ls)
            df = pd.DataFrame(ls)
            df.columns = titles         
            df['ID'] = userID
            df['rep'] = rep
            lateralAvoidDF = pd.concat([lateralAvoidDF, df ], ignore_index=True)
            
    
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
    
    foutImg = os.path.join(fout + "Area per level distribution.png")    
    plt.savefig(f'{foutImg}', dpi=cdpi)
        
    
    # ------------------
    # Box with mustaches 
    plt.clf()
    
    fig, ax = plt.subplots()
    # Create names on the x axis
    ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)
    
    print(pointsAvoidSumDF)
    
    for rep in repList:
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)
        
        ax = sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area")
        
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.4 )
        #sns.stripplot(data=df, x='day', y='tip', color='black', alpha=0.5)
        
        plt.xticks(rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot Area to obst per level distribution_rep" + str(rep) + ".png")    
        plt.savefig(f'{foutImg}', dpi=10*cdpi)
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgCount")
        
        plt.xticks(rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot avgCount to obst per level distribution_rep" + str(rep) + ".png")    
        plt.savefig(f'{foutImg}', dpi=(10*cdpi))
        
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)

        sns.boxplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralCount")
        
        plt.xticks(rotation=45)
         
        foutImg = os.path.join(fout + "BoxPlot lateral to obst per level distribution_rep" + str(rep) + ".png")    
        plt.savefig(f'{foutImg}', dpi=(10*cdpi))
        
        
        
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgCount",
                            inner='box', scale='width', dodge=False, linewidth=1)
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgCount", jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.4 )
        
        plt.xticks(rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgCount to obst_rep" + str(rep) + ".png")    
        plt.savefig(f'{foutImg}', dpi=(10*cdpi))
    
    
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        # Create names on the x axis
        ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)

        sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "area",
                            inner='box', scale='width', dodge=True, linewidth=1)
        
        sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "area", jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.4 )
        
        plt.xticks(rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot area to obst_rep" + str(rep) + ".png")    
        plt.savefig(f'{foutImg}', dpi=(10*cdpi))
    
           
        plt.clf()
        plt.figure(figsize=(40,15))
        
        fig, ax = plt.subplots()
        
        # Create names on the x axis
        #ax.set_xticks(np.arange(len(constLevelList)), labels=constLevelList)

        ax = sns.violinplot(data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ], x="level", y = "avgLateralCount",
                            inner='box', scale='width', dodge=False, linewidth=1)
        ax = sns.stripplot( data=pointsAvoidSumDF[ (pointsAvoidSumDF["rep"] == rep) ],  x="level", y = "avgLateralCount", jitter=True, linewidth=0.1, color = 'darkblue', alpha = 0.4 )
        
        plt.xticks(rotation=45)
         
        foutImg = os.path.join(fout + "ViolinePlot avgLateralCount to obst_rep" + str(rep) + ".png")    
        plt.savefig(f'{foutImg}', dpi=(10*cdpi))
    
    
    return 



# Calculate time of task completion 
def timeComplCalculation(pathDF, fout):
        
    
    XborderMin = -120
    XborderMax = 120
    
    Xobst = 0
    Yobst = 0
    
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
            #timeComplDF =timeComplDF.append(df)
            timeComplDF  =  pd.concat([timeComplDF, df ], ignore_index=True)
            
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
cdpi = 120
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
XlimMax = 250
XlimMin = -250

#colorLevelList = ['aquamarine2','brown1', 'burlywood2', 'cadetblue2', 'darkgoldenrod1', 'darkolivegreen2', 'darkorange','darkorchid1','darkturquoise', 'hotpink','indianred1', 'lightpink1', 'lightsalmon1','magenta2', 'maroon1','mediumorchid1','violetred1' ]

#colorLevelListRed = ['darkorchid1','darkorchid1','darkorchid1','darkturquoise', 'hotpink','indianred1', 'lightpink1', 'lightsalmon1','magenta2', 'maroon1','mediumorchid1','violetred1','palevioletred1', 'orchid1' ]
#colorLevelListBlue = ['lightblue','lightblue','lightblue','lightblue', 'darkturquoise','magenta2', 'orchid1','mediumorchid1','violetred1','palevioletred1' ]
#redList = ['orange1', 'orangered1', 'orchid1']


#colorLevelListRed = ['darkorchid','darkorchid','darkorchid','darkturquoise', 'hotpink','indianred', 'lightpink', 'lightsalmon','magenta', 'maroon','mediumorchid','darkorchid','deeppink', 'orchid' ]
colorLevelListRed = ['grey', 'darkgrey', 'royalblue','aquamarine','darkorchid','darkturquoise', 'gold','mediumvioletred', 'lightpink', 'deeppink','magenta', 'mediumorchid','darkorchid','deeppink', 'violet' ]


colorUserList = ['grey', 'darkgrey', 'royalblue',
                 'darkturquoise','green', 'olive', 'gold',
                 'bisque', 'coral', 'deeppink', 
                 'crimson', 'brown',  'plum' ,
                 'magenta',
                  'burlywood', 'cadetblue', 'darkgoldenrod', 'darkolivegreen', 
                  'darkorange','darkturquoise', 'hotpink','indianred', 
                  'lightsalmon', 'maroon','violetred', 'salmon', 'salmon', 'salmon', 'salmon', 'salmon']

colorLevelListBlue = ['lightblue','lightblue','lightblue','lightblue', 'darkturquoise',
                      'mediumslateblue', 'royalblue','aquamarine','blueviolet', 'violet', 
                      'cadetblue', 'palevioletred', 'blue' ]
redList = ['orange', 'orangered', 'orchid']


constLevelList = ['level100', 'level110', 'level108', 'level107', 'vel108_2', 'vel107_2', 'level101',
                  'level102', 'level103', 'level104', 'level105', 'level106']

#constLevelList = ['level0_1', 'level0_2', 'level0_3', 'o_level1', 'o_level2', 'o_level3', 'o_level4', 'o_level5', 'o_level6']

#stat_data
statDF = pd.DataFrame(columns=['ID', 'level', 'rep', 'x', 'y'])
resTSstatDF = pd.DataFrame(columns=['ID', 'level', 'rep', 'x', 'y'])
resTSXYstatDF = pd.DataFrame(columns=['ID', 'level', 'rep', 'x', 'y'])

#merged DF per experiment
mergedDF = pd.DataFrame()
mergedTSDF = pd.DataFrame()
mergedTSXYDF = pd.DataFrame()

trackingDataLog = pd.DataFrame()



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
                        OX=float(coords[0])
                        OY=float(coords[1])
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
                    OXarr.append(OX)
                    
                        
            
            trackDF = pd.DataFrame({'X': Xarr,
                                    'Y': Yarr,
                                    'Z': Zarr,
                                    'YA': YAarr,
                                    'P': Parr,
                                    'R': Rarr,
                                    'T': Tarr,
                                    'DT': DTarr,
                                    'OX': OXarr,
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
            trackDF = trackDF.sort_values('timestamp')
            
            fout = os.path.join(modified_csv_directory, filename)
            trackDF.to_csv(fout)    
            
            df = trackDF.copy()
            df["ID"] = id
            df['level']= filename[-18:-10]
            df['rep' ]=  rep
            #mergedDF = mergedDF.append(df)
            
            mergedDF = pd.concat([mergedDF, df ], ignore_index=True)
            
            
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
                
                foutImg = os.path.join(img_directory,name+"_resample.png")    
                #plt.savefig(f'{foutImg}', dpi=my_dpi)
                
            
            # >>> Time resample to a new FPS
            # Linear Approximation for filtered data
            lbl = "Tsample " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                
                (TResArr ,x ,y, z ) = resample3D_scale(
                    np.array(trackDF["timestamp"]), 
                    np.array(trackDF["X"]), 
                    np.array(trackDF["Y"]),
                    np.array(trackDF["Z"]), FPSTarget)
                
                (TResArr ,ya ,p, r ) = resample3D_scale(
                    np.array(trackDF["timestamp"]), 
                    np.array(trackDF["YA"]), 
                    np.array(trackDF["P"]),
                    np.array(trackDF["R"]), FPSTarget)
                
                    
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
                mergedTSDF = pd.concat([mergedTSDF, df])
                
                
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
                    
                
                axs[rep, 1].plot(x,y, color=colorLevelListBlue[lvlID], linewidth= (6.2 ),label= lbl )       
                
                axs[rep, 0].set_title("Raw, rep = " + str(rep))
                axs[rep, 1].set_title("Filtered, rep = "+ str(rep))
                
                axs[rep, 0].plot( XEs, (YEs), color='limegreen',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
                         markersize=18, markevery=12)
                axs[rep, 1].plot( XEs, (YEs), color='limegreen',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
                         markersize=18, markevery=12)
                axs[rep, 1].legend( loc='lower left' , ncol = 2, frameon=True, prop={'size': 10});
                
                axs[rep, 0].legend( loc='lower left' , ncol = 2, frameon=True, prop={'size': 10});
                
                #titles.append(filename[-18:-10])
                # <<< Time resample to a new FPS
                
            
            # >>> Time + Space resample to a new FPS
            # Linear Approximation for filtered data
            lbl = "TXsample " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                
                (t ,x ,y, z, ya, p, r ) = resample_spatial_all(
                    np.array(resTSDF["timestamp"]), 
                    np.array(resTSDF["X"]), 
                    np.array(resTSDF["Y"]),
                    np.array(resTSDF["Z"]),
                    np.array(resTSDF["YA"]), 
                    np.array(resTSDF["P"]),
                    np.array(resTSDF["R"])
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
                #mergedTSXYDF = mergedTSXYDF.append(df)
                mergedTSXYDF =  pd.concat([mergedTSXYDF, df])

                pwlf_model = pwlf.PiecewiseLinFit(x, y ,degree=1)
                breaks = pwlf_model.fit(nSlopes)
                
                slopes = pwlf_model.calc_slopes() 
                
                resTSXYstatDF = pd.concat([resTSXYstatDF, pd.DataFrame({
                                                        'ID': [id], 
                                                        'level': [filename[-18:-10]],
                                                        'rep' : [rep],
                                                        'x': [breaks[1]] , 
                                                        'y': [pwlf_model.predict(breaks[1])[0]]
                                                        
                                                        })])
                """
                resTSXYstatDF = resTSXYstatDF.append({
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
            
          
    
    foutImg = os.path.join(img_directory,id+ ".png")      
    
    
    #fig.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
    fig.savefig(f'{foutImg}', dpi=cdpi)
 

    
foutImg = os.path.join(img_directory,"TOTAL"+"_lin"+str(nSlopes)+"_lowPass_" +str(fc) + ".png")      
    
fig.savefig(f'{foutImg}', dpi=cdpi)

foutStat = os.path.join(stat_directory,"trackingDataLog"+".xlsx")   
trackingDataLog.to_excel(foutStat)
    
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


"""
# statistics  vis
fout = os.path.join(stat_directory,"deviationStat"+"_lin"+str(nSlopes)+"_raw")   
correctRate = visStatRepSeperate(statDF, fout)
print( "Raw: ", correctRate )

fout = os.path.join(stat_directory,"deviationStat"+"_lin"+str(nSlopes)+"_resTS"+str(FPSTarget))   
correctRate = visStatRepSeperate(resTSstatDF, fout)
print("TS: ", correctRate )

fout = os.path.join(stat_directory,"deviationStat"+"_lin"+str(nSlopes)+"_resTSXY"+str(FPSTarget))   
correctRate = visStatRepSeperate(resTSXYstatDF, fout)    
print("TSXY: ", correctRate )
"""

#area
fout = os.path.join(stat_directory,"_raw_")   
areaCalculation(mergedDF, fout)

fout = os.path.join(stat_directory,"_resTS"+str(FPSTarget))   
areaCalculation(mergedTSDF, fout)

fout = os.path.join(stat_directory,"_resTSXY"+str(FPSTarget))   
areaCalculation(mergedTSXYDF, fout)

#area
fout = os.path.join(stat_directory,"_raw_")   
distCalculationToObst(mergedDF, fout)


    






