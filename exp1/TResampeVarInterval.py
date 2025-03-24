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


def resample3D(time, xdata, ydata, zdata, newFrameRate):
    
    #re sampling time 
    duration = time[-1] - time[0]
    nbSamples = (duration * newFrameRate) + 1
    newTime = np.linspace(time[0], time[-1], len(time))
    
    x = xdata[:]-xdata[0]
    y = ydata[:]-ydata[0]
    z = zdata[:]-zdata[0]
    
    #resample x data accordingly
    xInter = interp1d(time, x, kind='linear')
    new_xdata = xInter(newTime)
    
    #resample y data accordingly
    yInter = interp1d(time, y, kind='linear')
    new_ydata = yInter(newTime)
    
    #resample z data accordingly
    zInter = interp1d(time, z, kind='linear')
    new_zdata = zInter(newTime)
    
    #butterworth lowpass filter
    Wn = 2.0*fc/newFrameRate
    sos = signal.butter(1, Wn, 'low', output='sos')
    
    xdataFilter = signal.sosfilt(sos, new_xdata) +xdata[0]
    ydataFilter = signal.sosfilt(sos, new_ydata) +ydata[0]
    zdataFilter = signal.sosfilt(sos, new_zdata) +zdata[0]
    
    return (newTime, xdataFilter, ydataFilter, zdataFilter)


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
YlimMax = 110
YlimMin = -130
XlimMax = -30
XlimMin = -250


XEs = -118 # ES character
YEs = -110 
XNe = -118 # Neurotic character
YNe = 90

#stat_directory
statDF = pd.DataFrame(columns=['ID', 'level', 'x', 'y'])

statDFfilt = pd.DataFrame(columns=['ID', 'level', 'x', 'y'])


for id in iserIdList:
    
    plt.figure(figsize=(20,25), dpi= my_dpi)
    levelColor = 0
    titles = []
    
    #fig, axs = plt.subplots(2, 3, figsize=(40, 40), dpi = my_dpi)
    imgCount = 0
    for filename in os.listdir(raw_data_csv_directory): 
        
        name, ext = os.path.splitext(filename)
         
        
        
        if (( id in filename ) and (ext == '.csv') and ("Coord" in filename) and not( "init" in filename ) and filename[-11:-10] != '0' ) :
            
            
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
            
            if len(trackDF[ (trackDF["X"]<- 50) ][ (trackDF["X"]>- 170) ][(trackDF["Y"]<50)][(trackDF["Y"]>-50)]):
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
            fout = os.path.join(modified_csv_directory, filename)
            trackDF.to_csv(fout)
            # <<< Restore path of HDMI
                
            # >>> Count frames, FPS
            df = trackDF.copy()
            
            df['FNum'] = df.groupby('T')['T'].transform('size')
            
            df['FPS'] = df['FNum']
            
            FPSarr = trackDF['T'].value_counts(sort = False)
            
            Tbeg = df.iloc[0].loc["T"]   
            Tend = df.iloc[-1].loc["T"]
            
            df.loc[df["T"] == Tbeg, "FPS"] = FPSarr[1]
            df.loc[df["T"] == Tend, "FPS"] = FPSarr[-2]
            
            TIntervals = list(set(trackDF['T']))
            TIntervals.sort(reverse = False)
            # <<< Count frames, FPS
            
            # >>> Compute ms 
            tdf = df.loc[df["T"] == Tbeg].copy()
            FNum =tdf.iloc[0]['FNum']
            tMSarr = np.arange(0.00, 1-float(0.5/FNum), float(1/FNum))
            MSarr = np.append(MSarr, tMSarr)        
            # <<< Compute ms 
            
            for timestamp in TIntervals[1:-1]: 
                
                tdf = df.loc[df["T"] == timestamp].copy()
                
                FNum =tdf.iloc[0]['FNum']
                FPS =tdf.iloc[0]['FPS']
                    
                tMSarr = np.arange(0.00, 1-float(0.5/FNum), float(1/FNum))
                MSarr = np.append(MSarr, tMSarr)
                
                if FNum < 2:
                    cPolyType = 'constant'
                else:
                    cPolyType = 'line'
                
            # >>> Compute ms 
            tdf = df.loc[df["T"] == Tend].copy()
            FNum =tdf.iloc[0]['FNum']
            tMSarr = np.arange(0.00, 1-float(0.5/FNum), float(1/FNum))
            MSarr = np.append(MSarr, tMSarr)      
            # <<< Compute ms 
                
            
            #fout = os.path.join(modified_csv_directory, "resampled" + filename)
            #resampleDF.to_csv(fout)
            
            
            # >>> add MS and Timestamp with MS
            trackDF['FNum'] = df['FNum']
            trackDF['MS'] = MSarr
            
            TDarr = np.array(trackDF['DT'])
            TimeSarr = []
            
            for i, time_data in enumerate(TDarr):
            
                date = datetime.strptime(time_data, format_data)
                ts = date.timestamp() 
                ts+= MSarr[i]
                TimeSarr.append(ts)
            
            trackDF['timestamp'] = TimeSarr
            
            fout = os.path.join(modified_csv_directory, filename)
            trackDF.to_csv(fout)     
            # <<< add MS and Timestamp with MS
             
            Xarr = np.array(trackDF["X"])
            Yarr = np.array(trackDF["Y"])
                      
            #Linear Approximation of original signal
            lbl = "approx. for " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                #model = np.polyfit(Xarr, Yarr, 1)
                                
                pwlf_model = pwlf.PiecewiseLinFit(Xarr, Yarr ,degree=1)
                breaks = pwlf_model.fit(nSlopes)
                
                slopes = pwlf_model.calc_slopes()
                
                randColor = random.randint(0,10)
                    
                statDF = statDF.append({'A': i, 
                                        'ID': id, 
                                        'level': filename[-18:-10],
                                        'x': breaks[1] , 
                                        'y': pwlf_model.predict(breaks[1])
                                        
                                        }, ignore_index=True)
                
                #plt.plot(Xarr, Yarr, '.', color=[0.1*randColor,0.2*levelColor,1 - 0.2*levelColor], linewidth=0.02 , label=filename[-18:-10])
                                #,figure=plt.figure());            
                #randColor = random.randint(0,10)
                #for xb in breaks:
                #    plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                    
                titles.append(filename[-18:-10])
                
                plt.plot(XEs, YEs, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
             markersize=20, markevery=12)
                plt.plot(XEs, YEs, color='red', marker='h', markerfacecolor='lightcoral', markeredgewidth=12,
             markersize=20, markevery=12)
                
                foutImg = os.path.join(img_directory,name+"_resample.png")    
                #plt.savefig(f'{foutImg}', dpi=my_dpi)
            

            #Linear Approximation for filtered data
            lbl = "approx. for " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                
                (TResArr ,x ,y, z ) = resample3D(
                    np.array(trackDF["timestamp"]), 
                    np.array(trackDF["X"]), 
                    np.array(trackDF["Y"]),
                    np.array(trackDF["Z"]), FPSTarget)
                
                (TResArr ,ya ,p, r ) = resample3D(
                    np.array(trackDF["timestamp"]), 
                    np.array(trackDF["YA"]), 
                    np.array(trackDF["P"]),
                    np.array(trackDF["R"]), FPSTarget)
                
                    
                filtDF = pd.DataFrame({'X': x,
                                        'Y': y,
                                        'Z': z,
                                        'YA': ya,
                                        'P': p,
                                        'R': r,
                                        'timestamp': TResArr,
                                        })

                pwlf_model = pwlf.PiecewiseLinFit(x, y ,degree=1)
                breaks = pwlf_model.fit(nSlopes)
                
                slopes = pwlf_model.calc_slopes()
                
                statDFfilt = statDFfilt.append({
                                        'ID': id, 
                                        'level': filename[-18:-10],
                                        'x': breaks[1] , 
                                        'y': pwlf_model.predict(breaks[1])
                                        
                                        }, ignore_index=True)
                
                randColor = random.randint(0,10)
                for xb in breaks:
                    
                    plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                
                plt.plot(x, pwlf_model.predict(x), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= lbl)
                plt.plot(x,y, '.', color=[0.524,0.2*levelColor,1 - 0.2*levelColor], linewidth=0.18 , label="filtered "+filename[-18:-10])       
                titles.append(filename[-18:-10])
                
                fout = os.path.join(modified_csv_filt_dir, str(id) + "_lowPass_" +str(fc) +"_"+ filename)
                filtDF.to_csv(fout)
                
            
            levelColor =+ 0.05
    
    foutImg = os.path.join(img_directory,id+"_lin"+str(nSlopes)+"_lowPass_" +str(fc) + ".png")      
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
    
    lgnd.legendHandles[0]._legmarker.set_markersize(40)
    
    for l in lgnd.legendHandles:
        l._legmarker.set_markersize(40)
    #lgnd.legendHandles[1]._legmarker.set_markersize(40)
    lbl
    
    plt.savefig(f'{foutImg}', dpi=my_dpi)
            
    foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resample.csv")   
    statDF.to_csv(foutStat)
    
    
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resample.csv")   
statDF.to_csv(foutStat)
foutStat = os.path.join(stat_directory,"deviation"+"_lin"+str(nSlopes)+"_resample.xlsx") 
statDF.to_excel(foutStat,  index = False, header=True) 
    
plt.clf()
plt.figure(figsize=(10,10), dpi= 10)

correctRate = [0,0] # [0] - correct ; [1] - not

for i, id in enumerate (iserIdList):
    proxArr = [0,0,0,0,0] # metrics for 5 experiments 
    
    idDF = statDF[statDF['ID'] == id]
    levelColor = int(id)/ 100000
    
    for index, row in idDF.iterrows():
        curLevel = int(row['level'][-1:])
        print(curLevel)
        if (row['level'][-2:-1] != '_') and curLevel !=0:
            proxArr[curLevel] = row['x']
        elif (row['level'][-2:-1] != '_'):
            proxArr[0] = float(row['x'])
         
    plt.plot(proxArr, '.', color=[0.2, levelColor ,1 - levelColor], linewidth=1.4 , label=id )
    plt.plot(proxArr, color=[0.2, levelColor ,1 - levelColor], linewidth=0.42  )
    
    if i == len(iserIdList) -1: 
                
        #lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
        lgnd  = plt.legend( loc='lower center' , ncol = 4, frameon=True, fontsize = 'x-small');
        
        foutImg = os.path.join(img_directory,"deviation_part1_lin"+str(nSlopes)+"_lowPass_" +str(fc) + ".png")    
        plt.savefig(f'{foutImg}', dpi=my_dpi)
                    
        plt.clf()
        plt.figure(figsize=(10,10), dpi= 10)
        
    if proxArr[3] < proxArr[4]:
        correctRate[0] =correctRate[0] +1
    else:
        correctRate[1] =correctRate[1] +1
        

print("correctRate ", correctRate)   

correctRate = [0,0] # [0] - correct ; [1] - not
for i, id in enumerate (iserIdList):
    proxArr = [0,0,0,0,0] # metrics for 5 experiments 
    
    idDF = statDFfilt[statDFfilt['ID'] == id]
    levelColor = int(id)/ 100000
    
    for index, row in idDF.iterrows():
        curLevel = int(row['level'][-1:])
        print(curLevel)
        if (row['level'][-2:-1] != '_') and curLevel !=0:
            proxArr[curLevel] = row['x']
        elif (row['level'][-2:-1] != '_'):
            proxArr[0] = float(row['x'])
         
    plt.plot(proxArr, '.', color=[0.2, levelColor ,1 - levelColor], linewidth=1.4 , label=id )
    plt.plot(proxArr, color=[0.2, levelColor ,1 - levelColor], linewidth=0.42  )
    
                
        #lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
lgnd  = plt.legend( loc='lower center' , ncol = 4, frameon=True, fontsize = 'x-small');
        
foutImg = os.path.join(img_directory,"deviation_part1_lin"+str(nSlopes)+"_lowPass_" +str(fc) + ".png")       
plt.savefig(f'{foutImg}', dpi=my_dpi)
                    
plt.clf()
plt.figure(figsize=(10,10), dpi= 10)
        
if proxArr[3] < proxArr[4]:
    correctRate[0] =correctRate[0] +1
else:
    correctRate[1] =correctRate[1] +1
        

print("correctRate ", correctRate)   


        
        
        
        