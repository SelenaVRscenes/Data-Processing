import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import re
import csv
import os
import fnmatch
from scipy import optimize
from scipy import interpolate
from scipy import signal

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


plt.style.use('seaborn-whitegrid')

path_of_the_directory= os.path.join("C:\\Users\ypatotsk\Documents\Dev\pythonDev", "trackingData")
raw_data_csv_directory = os.path.join(path_of_the_directory, "rawData")
#path_of_the_directory= os.path.join("C:\\Dev\\pythonDev", "trackingData")

modified_csv_directory= os.path.join( path_of_the_directory, "csvCleanedData")

img_directory= os.path.join( path_of_the_directory, "csvVisData")

stat_directory= os.path.join( path_of_the_directory, "statistics")
    
outputDir = os.path.join(path_of_the_directory, "img")

if not (os.path.exists(outputDir)) :
    os.mkdir(outputDir)
    
if not (os.path.exists(modified_csv_directory)) :
    os.mkdir(modified_csv_directory)

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
my_dpi = 300
nSlopes = 2

#stat_directory
statDF = pd.DataFrame(columns=['ID', 'level', 'x', 'y'])

for id in iserIdList:
    
    plt.figure(figsize=(25,25), dpi= my_dpi)
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
                           
                    if not((X< -250) or (X> -30) or (Y> 100)  or (Y< -100)):
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
            
            #name, ext = os.path.splitext(filename)
            
            # Restore path
            indMed = int(trackDF[ (trackDF["X"]<- 50) ][ (trackDF["X"]>- 170) ][(trackDF["Y"]<50)][(trackDF["Y"]>-50)].index[-1])
            indStrat = 0
            indEnd = len(trackDF)
            
            trackDFclean = pd.DataFrame()
            
            if len(trackDF[ (trackDF["X"]<- 50) ][ (trackDF["X"]>- 170) ][(trackDF["Y"]<50)][(trackDF["Y"]>-50)]):
                for i, r in trackDF.iloc[indMed:].iterrows():
                    X = r.X
                    Y = r.Y
                    if ((X< -210) or (X> -30) or (Y> 100)  or (Y< -100)):
                        indEnd = i
                        print("--end index allert--", indStrat, indEnd)
                        break
                    
                for i,r in trackDF[::-1].iterrows():
                    X = r.X
                    Y = r.Y
                    if ((X< -210) or (X> -30) or (Y> 100)  or (Y< -100)):
                        indStrat = i
                        print("--begin index allert--", indStrat, indEnd)
                        break
            
            trackDF = trackDF[indStrat:indEnd]
            fout = os.path.join(modified_csv_directory, filename)
            trackDF.to_csv(fout)
                
            # Timewise Resampling

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
            
            resampleDF = df.loc[df["T"] == TIntervals[0]].copy()
            
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
                #FNumTarget = int(FPSTarget*FNum/FPS)
                #print(FNumTarget)
                
                #tXarr = signal.resample(tdf["X"], FPSTarget)
                #tYarr = signal.resample(tdf["Y"], FPSTarget)
                #tZarr = signal.resample(tdf["Z"], FPSTarget)
                    
                #tYAarr = signal.resample(tdf["YA"], FPSTarget)
                #tParr = signal.resample(tdf["P"], FPSTarget)
                #tRarr = signal.resample(tdf["R"], FPSTarget)
                #tTarr = np.full(FPSTarget, timestamp)
                #tDTarr = np.zeros( FPSTarget)
                
                #scipy.signal.resample_poly(y, 20, 256)
                if FNum < 2:
                    cPolyType = 'constant'
                else:
                    cPolyType = 'line'
                
                tXarr = signal.resample_poly(tdf["X"],FPSTarget, FNum, padtype= cPolyType)
                tYarr = signal.resample_poly(tdf["Y"],FPSTarget, FNum, padtype=cPolyType)
                tZarr = signal.resample_poly(tdf["Z"], FPSTarget, FNum, padtype= cPolyType)
                    
                tYAarr = signal.resample_poly(tdf["YA"], FPSTarget, FNum, padtype= cPolyType)
                tParr = signal.resample_poly(tdf["P"], FPSTarget, FNum, padtype= cPolyType)
                tRarr = signal.resample_poly(tdf["R"], FPSTarget, FNum, padtype= cPolyType)
                tTarr = np.full(FPSTarget, timestamp)
                tDTarr = np.zeros( FPSTarget)
                
                resampleDF = resampleDF.append(  pd.DataFrame({'X': tXarr,
                                    'Y': tYarr,
                                    'Z': tZarr,
                                    'YA': tYAarr,
                                    'P': tParr,
                                    'R': tRarr,
                                    'T': tTarr,
                                    'DT': tDTarr,
                                    }), ignore_index=True)
            
            
            resampleDF = resampleDF.append( df.loc[df["T"] == TIntervals[-1]], ignore_index=True)    
            
            # >>> Compute ms 
            tdf = df.loc[df["T"] == Tend].copy()
            FNum =tdf.iloc[0]['FNum']
            tMSarr = np.arange(0.00, 1-float(0.5/FNum), float(1/FNum))
            MSarr = np.append(MSarr, tMSarr)      
            # <<< Compute ms 
                
            
            fout = os.path.join(modified_csv_directory, "resampled" + filename)
            
            resampleDF.to_csv(fout)
            
            
            # >>> add ms
            trackDF['FNum'] = df['FNum']
            trackDF['MS'] = MSarr
            fout = os.path.join(modified_csv_directory, filename)
            trackDF.to_csv(fout)
            
            # <<< add ms
                
                        
            #Linear Approximation
            lbl = "approx. for " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                #model = np.polyfit(Xarr, Yarr, 1)
                
                x = np.array(Xarr)      
                y = np.array(Yarr)
                
                #condlist, funclist = optimize.curve_fit(piecewise_linear, x, y)
                
                pwlf_model = pwlf.PiecewiseLinFit(x, y ,degree=1)
                breaks = pwlf_model.fit(nSlopes)
                
                slopes = pwlf_model.calc_slopes()
                
                randColor = random.randint(0,10)
                
                #for xb in breaks:
                    
                    #plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                
                #plt.plot(x, pwlf_model.predict(x), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= lbl)
                
                    
                statDF = statDF.append({'A': i, 
                                        'ID': id, 
                                        'level': filename[-18:-10],
                                        'x': breaks[1] , 
                                        'y': pwlf_model.predict(breaks[1])
                                        
                                        }, ignore_index=True)
                
        
                #plt.plot(Xarr, Yarr, '.', color=[0.1*randColor,0.2*levelColor,1 - 0.2*levelColor], linewidth=0.02 , label=filename[-18:-10])
                                #,figure=plt.figure());
                                
                    
                titles.append(filename[-18:-10])
                
                plt.plot(-110, -121, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
             markersize=20, markevery=12)
                plt.plot(-116, 112, color='red', marker='h', markerfacecolor='lightcoral', markeredgewidth=12,
             markersize=20, markevery=12)
                
                foutImg = os.path.join(img_directory,name+"_resample.png")    
                #plt.savefig(f'{foutImg}', dpi=my_dpi)
            

            Xarr = resampleDF['X'].to_numpy()
            Yarr = resampleDF['Y'].to_numpy()
            
            
            
            #Linear Approximation
            lbl = "approx. for " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                x = np.array(Xarr)      
                y = np.array(Yarr)
                
                #condlist, funclist = optimize.curve_fit(piecewise_linear, x, y)
                
                pwlf_model = pwlf.PiecewiseLinFit(x, y ,degree=1)
                breaks = pwlf_model.fit(nSlopes)
                
                slopes = pwlf_model.calc_slopes()
                
                import random
                randColor = random.randint(0,10)
                
                for xb in breaks:
                    
                    plt.plot(xb, pwlf_model.predict(xb), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor] ,marker='h', markerfacecolor='lightgreen', markeredgewidth=12, markersize=20, markevery=10)
                
                plt.plot(x, pwlf_model.predict(x), color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= lbl)
                
        
                plt.plot(Xarr, Yarr, '.', color=[0.1*randColor,0.2*levelColor,1 - 0.2*levelColor], linewidth=0.12 , label="res "+filename[-18:-10])
                                #,figure=plt.figure());
                                
                    
                titles.append(filename[-18:-10])
                
            
            levelColor = levelColor+1
    
    foutImg = os.path.join(img_directory,id+"_lin"+str(nSlopes)+"_resample.png")   
    #plt.legend(bbox_to_anchor =(1.75,2.15), ncol = 2) 
    
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
    
    if i == int (len(iserIdList) / 2): 
                
        #lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
        lgnd  = plt.legend( loc='lower center' , ncol = 4, frameon=True, fontsize = 'x-small');
        
        foutImg = os.path.join(img_directory,"deviation_part1_lin"+str(nSlopes)+"_resample.png")   
        plt.savefig(f'{foutImg}', dpi=my_dpi)
                    
        plt.clf()
        plt.figure(figsize=(10,10), dpi= 10)
        
    if proxArr[3] < proxArr[4]:
        correctRate[0] =correctRate[0] +1
    else:
        correctRate[1] =correctRate[1] +1
        

print("correctRate ", correctRate)   

#lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
lgnd  = plt.legend( loc='lower center' , ncol = 4, frameon=True, fontsize = 'x-small');

foutImg = os.path.join(img_directory,"deviation_part2_lin"+str(nSlopes)+"_resample.png")   
plt.savefig(f'{foutImg}', dpi=my_dpi)


from scipy.fftpack import rfft, irfft, fftfreq, fft

dT = float(1/FPSTarget)
freq = 0.90
Y = resampleDF["Y"].to_numpy()

X = resampleDF["X"].to_numpy()

yf = fft(Y)
 
f_signal = rfft(Y)

W = fftfreq(y.size, dT)

plt.figure(figsize=(12,12), dpi= my_dpi)

for freq  in np.arange(0.75, 2.5, 0.1):
    
    cut_f_signal = f_signal.copy()
    cut_f_signal[(W>freq)] = 0  # filter all frequencies above 0.6
    randColor = random.randint(0,10)
                
    cut_signal = irfft(cut_f_signal)
    lbl= "freq " + str(freq)
    plt.plot(Xarr, cut_signal, '.', color=[freq/5.5,0.2*levelColor,1-0.1*randColor], linewidth=0.08 , label=lbl)
 
                                   
plt.plot(Xarr, Yarr, '-', color=[1-0.1*randColor,0.2*levelColor,0.2*levelColor], linewidth=0.32 , label="initial signal")
      
#lgnd  = plt.legend( loc='top' , ncol = 2, frameon=True);
lgnd  = plt.legend()       

foutImg = os.path.join(img_directory,"filter freq "+str(freq)+"_resample.png")   
plt.savefig(f'{foutImg}', dpi=my_dpi)
     
        
        
        
        