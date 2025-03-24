import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import re
import csv
import os
import fnmatch
from scipy import optimize
from scipy import interpolate
import pwlf
import random

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


iserIdList = []
extractId = []

for filename in os.listdir(raw_data_csv_directory): 
    ids = re.findall(r"[-+]?\d*\.*\d+", filename)
    id = ids[0]
    extractId.append(id)
    
iserIdList = set(extractId)

savePerUser = 0

#stat_directory
statDF = pd.DataFrame(columns=['ID', 'level', 'x', 'y'])

for id in iserIdList:
    
    my_dpi = 100
    nSlopes = 2
    
    plt.figure(figsize=(25,25), dpi= my_dpi)
    levelColor = 0
    titles = []
    
    #fig, axs = plt.subplots(2, 3, figsize=(40, 40), dpi = my_dpi)
    imgCount = 0
    for filename in os.listdir(raw_data_csv_directory): 
        
        name, ext = os.path.splitext(filename)
         
        if (( id in filename ) and (ext == '.csv') and ("Coord" in filename) and not( "init" in filename ) and filename[-11:-10] != '0'):
            
            
            beginIndex =0
            endIndex =0
            
            indexMed =0
        
            Xarr = []
            Yarr = []
            Zarr = []
            
            Tarr = []
            DTarr = []
                            
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
                
                finIndex = 0
                startIndex = 0
                
                nonzeroElem = 0
                for x in tableCoord.iloc[::-1].itertuples():
                    print(x)
                    if x[1].count("=0.000") == 3:
                        #finIndex = x[0]
                        nonzeroElem = 0
                        print("delete zeros", x[1])
                    else:
                        nonzeroElem = nonzeroElem+1
                        if nonzeroElem > 10:
                            finIndex = x[0]
                            break
                    
                
                for x in tableCoord[:finIndex].itertuples():
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
    
                    
                tempDF =tableCoord[startIndex:finIndex]
                
                print("startIndex , finIndex")
                print(startIndex, finIndex)
                print("len(tableCoord[startIndex:])", len(tableCoord[startIndex:]))
                
                
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
                           
                    #if not((X< -250) or (X> -30) or (Y> 100)  or (Y< -100)):
                    Xarr.append(X)
                    Yarr.append(Y)
                    Zarr.append(Z)
                        
                    YAarr.append(YA)
                    Parr.append(P)
                    Rarr.append(R)
                        
                    Tarr.append(T)
                    DTarr.append(dateTime)
                    
                        
            endIndex = len(Xarr) -1   
            
            beginIndex = 0  
            
            
            for ind, x in enumerate(Xarr):
                
                if (x<=-100 and x>=-120 and Yarr[ind]<=100 and Yarr[ind] >=-100):
                    print(ind)
                    indexMed = ind
                    break
            
            if indexMed >0:
                for ind, y in enumerate(Yarr[indexMed:]):  #restore path to the end
                    if (y < - 120 or y > 120) and (Xarr[ind]>-100):
                        endIndex = ind
                        break
                    #if (ind >2) and (-Xarr[ind] + Xarr[ind-2])<0 :
                    #    print("-----------")
                    #    print(ind)
                    #    print("Xarr[ind] - Xarr[ind-2]", Xarr[ind] - Xarr[ind-2])
                    #    endIndex = ind
                    #    break
                          
                for ind, x in reversed( list(enumerate( Xarr[:indexMed] ) ) ):
                    if x < - 230 :
                         beginIndex = ind
                         break
                     
            
            print("-------Xarr--------")
            print(Xarr)
            
            
            Xarr_b = Xarr
            Yarr_b = Yarr
            
            Xarr = Xarr[beginIndex:endIndex]
            Yarr = Yarr[beginIndex:endIndex]
            Zarr = Zarr[beginIndex:endIndex]
            YAarr = YAarr[beginIndex:endIndex]
            Parr = Parr[beginIndex:endIndex]
            Rarr = Rarr[beginIndex:endIndex]
            Tarr = Tarr[beginIndex:endIndex]
            DTarr = DTarr[beginIndex:endIndex]
            
            print(" ---[beginIndex:endIndex] " , beginIndex,endIndex)
            
            trackDF = pd.DataFrame({'X': Xarr[beginIndex:endIndex],
                                    'Y': Yarr[beginIndex:endIndex],
                                    'Z': Zarr[beginIndex:endIndex],
                                    'YA': YAarr[beginIndex:endIndex],
                                    'P': Parr[beginIndex:endIndex],
                                    'R': Rarr[beginIndex:endIndex],
                                    'T': Tarr[beginIndex:endIndex],
                                    'DT': DTarr[beginIndex:endIndex],
                                    })
            
            #name, ext = os.path.splitext(filename)
            
            fout = os.path.join(modified_csv_directory, filename)
                
            
            trackDF.to_csv(fout)
            
            #Linear Approximation
            
            lbl = " trajectory for " + filename[-18:-10]
            if (len(Xarr)>0 and len(Yarr)>0):
                
                #model = np.polyfit(Xarr, Yarr, 1)
                
                x = np.array(Xarr)      
                y = np.array(Yarr)
                ya = np.array(YAarr)
                p = np.array(Parr)
                r = np.array(Rarr)
                
                randColor = random.randint(0,10)
                    
                if savePerUser > 0:
                    plt.clf()
                
                plt.plot(x, y, color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= "Sparce " + lbl, linestyle = '-', linewidth=4)
                
                plt.plot(x, ya, color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= " Y angle " + lbl,  linestyle = '--', linewidth=8)
                
                plt.plot(x, p, color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= " P angle " + lbl,  linestyle = ':', linewidth=8)
                
                plt.plot(x, r, color=[0.1*randColor ,0.2*levelColor,1 - 0.2*levelColor], label= " R angle " + lbl,  linestyle = '-.', linewidth=8)

                                
                    
                titles.append(filename[-18:-10])
                
                plt.plot(-110, -121, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
             markersize=20, markevery=12)
                plt.plot(-116, 112, color='red', marker='h', markerfacecolor='lightcoral', markeredgewidth=12,
             markersize=20, markevery=12)
                                
                
                if savePerUser > 0:    
                    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
                    
                    lgnd.legendHandles[0]._legmarker.set_markersize(40)
                    
                    for l in lgnd.legendHandles:
                        l._legmarker.set_markersize(40)
                        
                    foutImg = os.path.join(img_directory,id+"_angles_"+lbl+".png") 
                    plt.savefig(f'{foutImg}', dpi=my_dpi) 

                
            
            levelColor = levelColor+1
    
    foutImg = os.path.join(img_directory,id+"_angles"+".png")   
    #plt.legend(bbox_to_anchor =(1.75,2.15), ncol = 2) 
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 20});
    
    lgnd.legendHandles[0]._legmarker.set_markersize(40)
    
    for l in lgnd.legendHandles:
        l._legmarker.set_markersize(40)
    #lgnd.legendHandles[1]._legmarker.set_markersize(40)
    
    plt.savefig(f'{foutImg}', dpi=my_dpi)
    
    
    
plt.clf()
plt.figure(figsize=(10,10), dpi= 10)

