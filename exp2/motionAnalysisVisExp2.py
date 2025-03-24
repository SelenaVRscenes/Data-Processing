import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import re
import csv
import os
import fnmatch

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


#plt.style.use('seaborn-whitegrid')



dir = os.getcwd()

path_of_the_directory= os.path.join(dir, "trackingData")
#path_of_the_directory= os.path.join("C:\\Users\ypatotsk\Documents\Dev\pythonDev", "trackingData")
raw_data_csv_directory = os.path.join(path_of_the_directory, "rawData")
#path_of_the_directory= os.path.join("C:\\Dev\\pythonDev", "trackingData")

modified_csv_directory= os.path.join( path_of_the_directory, "csvMod")

outputDir = os.path.join(path_of_the_directory, "img")
if not (os.path.exists(outputDir)) :
    os.mkdir(outputDir)
    
if not (os.path.exists(modified_csv_directory)) :
    os.mkdir(modified_csv_directory)


                    
                    
outputDir = os.path.join(path_of_the_directory, "imgSorted")
if not (os.path.exists(outputDir)) :
    os.mkdir(outputDir)


listDir = os.listdir(modified_csv_directory)

iserIdList = []
extractId = []

for filename in os.listdir(modified_csv_directory): 
    ids = re.findall(r"[-+]?\d*\.*\d+", filename)
    id = ids[0]
    extractId.append(id)
    
iserIdList = set(extractId)

for id in iserIdList:
    #plt.clf() 
    
    my_dpi = 100
    
    fig, axs = plt.subplots(2, 3, figsize=(40, 40), dpi=my_dpi)
    imgCount = 0
    for filename in os.listdir(modified_csv_directory): 
        
        name, ext = os.path.splitext(filename)
        
        if (( id in filename ) and (ext == '.csv') and ("Coord" in filename) and not( "init" in filename )):
            
        
            f = os.path.join(modified_csv_directory, filename)
            if os.path.isfile(f):
                
                with open(f) as csvfile:
                    
                    print(f)
                    
                    #plt.figure(figsize=[25,25])
                    
                    reader = csv.reader(csvfile)
                    
                    Xarr = []
                    Yarr = []
                    
                    Parr = []
                    Rarr = []
                    
                    dataset = []
                    
                    for row in reversed( list(reader)):
                        str = row[0]
                           
                        coords = re.findall(r"[-+]?\d*\.*\d+", str)
                        X=float(coords[0])
                        Y=float(coords[1])
                        Z=float(coords[2])
                        Xarr.append(X)
                        Yarr.append(-Y)
                        
                        
                    # group images
                    (ax1, ax2) = divmod(imgCount, 3)
                    print(ax1, ax2)
                    
                    
                    axs[ax1, ax2].plot(0, 0, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
         markersize=40, markevery=12)
                    
                    axs[ax1, ax2].plot(Xarr, Yarr, 'g-o', color='blue', linewidth=2)
                   
                        
                        
                    axs[ax1, ax2].set_title(filename[-20:-14])
                    imgCount = imgCount+1;
                
                for ax in axs.flat:
                    ax.set(xlabel='x', ylabel='y')
                
                # Hide x labels and tick labels for top plots and y ticks for right plots.
                for ax in axs.flat:
                    ax.label_outer()

                fout = os.path.join(outputDir,id+".png")    
                fig.savefig(fout)
                
    imgCount = 0
        
    
    
    
    
    
for id in iserIdList:
    #plt.clf() 
    
    my_dpi = 100
    
    fig, axs = plt.subplots(2, 3, figsize=(40, 40), dpi=my_dpi)
    imgCount = 0
    for filename in os.listdir(modified_csv_directory): 
        
        name, ext = os.path.splitext(filename)
        
        if (( id in filename ) and (ext == '.csv') and ("Coord" in filename) and not( "init" in filename )):
            
        
            lineAngle = ""    
            f = os.path.join(modified_csv_directory, filename)
            if os.path.isfile(f):
                
                fileStrAngle = filename[0:-13] + "Angle*.csv"
                fileAngle = find(fileStrAngle, modified_csv_directory)[0]
                
                csvfileAngle = open(fileAngle) 
                    
                readerAngle = csv.reader(csvfileAngle)
                
                listAngle = ( list(readerAngle))
                
                tableAngle = pd.DataFrame(listAngle)
                    
                (ax1, ax2) = divmod(imgCount, 3)
                print(ax1, ax2)
                

                        
                        
                # group images
                if  ("level2" in filename):
                    axs[ax1, ax2].plot(-100, 104, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
         markersize=40, markevery=12)
                    axs[ax1, ax2].plot(-96, -90, color='red', marker='h', markerfacecolor='lightcoral', markeredgewidth=12,
         markersize=40, markevery=12)
                elif  ("level0" in filename):
                    print()
                        
                else:
                    axs[ax1, ax2].plot(-110, -121, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
         markersize=40, markevery=12)
                    axs[ax1, ax2].plot(-116, 112, color='red', marker='h', markerfacecolor='lightcoral', markeredgewidth=12,
         markersize=40, markevery=12)
                    
                
                with open(f) as csvfile:
                    
                    print(f)
                    
                    print(fileAngle)
                    
                    #plt.figure(figsize=[25,25])
                    
                    reader = csv.reader(csvfile)
                    
                    Xarr = []
                    Yarr = []
                    
                    XANGarr = []
                    YANGarr = []
                    
                    Parr = []
                    Rarr = []
                    YAarr = []
                    
                    dataset = []
                    
                    
                    for index, row in enumerate(( list(reader))):
                        str = row[0]
                           
                        coords = re.findall(r"[-+]?\d*\.*\d+", str)
                        X=float(coords[0])
                        Y=float(coords[1])
                        Z=float(coords[2])
                        Xarr.append(X)
                        Yarr.append(-Y)
                        
                        #read line for angle  
                        if index < tableAngle.size:
                            str = tableAngle.iloc[index][0]
                            
                            coords = re.findall(r"[-+]?\d*\.*\d+", str)
                                                
                            #X= Roll
                            #Y= Pitch
                            #Z= Yaw
        
                            P=float(coords[0])
                            YA=float(coords[1])
                            R=float(coords[2])
                            Parr.append(P)
                            YAarr.append(Y)
                            Rarr.append(R)
                            
                            XANGarr.append( X+np.cos(np.deg2rad(YA))) 
                            YANGarr.append( Y+np.sin(np.deg2rad(YA))) 
                            
                            print(YA,"* in Rad = ",np.deg2rad(YA))
                            
                            axs[ax1, ax2].arrow(X, -Y, 3*np.cos(np.deg2rad(YA)), 3*np.sin(np.deg2rad(YA)) , color='green', linewidth=2.8)
                            
                    axs[ax1, ax2].plot(Xarr, Yarr, 'g-o', color='blue', linewidth=2)
                    
                        
                        
                    axs[ax1, ax2].set_title(filename[-20:-14])
                    imgCount = imgCount+1;
                
                for ax in axs.flat:
                    ax.set(xlabel='x', ylabel='y')
                
                # Hide x labels and tick labels for top plots and y ticks for right plots.
                for ax in axs.flat:
                    ax.label_outer()

                fout = os.path.join(outputDir,id+"_Angle_YA.png")    
                fig.savefig(fout)
                
    imgCount = 0
    
    
    
    
   
    
    
for id in iserIdList:
    #plt.clf() 
    
    my_dpi = 100
    
    fig, axs = plt.subplots(2, 3, figsize=(40, 40), dpi=my_dpi)
    imgCount = 0
    for filename in os.listdir(modified_csv_directory): 
        
        name, ext = os.path.splitext(filename)
        
        if (( id in filename ) and (ext == '.csv') and ("Coord" in filename) and not( "init" in filename )):
            
        
            lineAngle = ""    
            f = os.path.join(modified_csv_directory, filename)
            if os.path.isfile(f):
                
                fileStrAngle = filename[0:-13] + "Angle*.csv"
                fileAngle = find(fileStrAngle, modified_csv_directory)[0]
                
                csvfileAngle = open(fileAngle) 
                    
                readerAngle = csv.reader(csvfileAngle)
                
                listAngle = ( list(readerAngle))
                
                tableAngle = pd.DataFrame(listAngle)
                    
                (ax1, ax2) = divmod(imgCount, 3)
                print(ax1, ax2)
                

                        
                        
                # group images
                if  ("level2" in filename):
                    axs[ax1, ax2].plot(-100, 104, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
         markersize=40, markevery=12)
                    axs[ax1, ax2].plot(-96, -90, color='red', marker='h', markerfacecolor='lightcoral', markeredgewidth=12,
         markersize=40, markevery=12)
                elif  ("level0" in filename):
                    print()
                        
                else:
                    axs[ax1, ax2].plot(-110, -121, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
         markersize=40, markevery=12)
                    axs[ax1, ax2].plot(-116, 112, color='red', marker='h', markerfacecolor='lightcoral', markeredgewidth=12,
         markersize=40, markevery=12)
                    
                
                with open(f) as csvfile:
                    
                    print(f)
                    
                    print(fileAngle)
                    
                    #plt.figure(figsize=[25,25])
                    
                    reader = csv.reader(csvfile)
                    
                    Xarr = []
                    Yarr = []
                    
                    XANGarr = []
                    YANGarr = []
                    
                    Parr = []
                    Rarr = []
                    YAarr = []
                    
                    dataset = []
                    
                    
                    for index, row in enumerate(( list(reader))):
                        str = row[0]
                           
                        coords = re.findall(r"[-+]?\d*\.*\d+", str)
                        X=float(coords[0])
                        Y=float(coords[1])
                        Z=float(coords[2])
                        Xarr.append(X)
                        Yarr.append(-Y)
                        
                        #read line for angle     
                        
                        if index < tableAngle.size:
                            str = tableAngle.iloc[index][0]
                            
                            coords = re.findall(r"[-+]?\d*\.*\d+", str)
                                                
                            #X= Roll
                            #Y= Pitch
                            #Z= Yaw
        
                            P=float(coords[0])
                            YA=float(coords[1])
                            R=float(coords[2])
                            Parr.append(P)
                            YAarr.append(Y)
                            Rarr.append(R)
                            
                            XANGarr.append( X+np.cos(np.deg2rad(YA))) 
                            YANGarr.append( Y+np.sin(np.deg2rad(YA))) 
                            
                            print(YA,"* in Rad = ",np.deg2rad(YA))
                            
                            axs[ax1, ax2].arrow(X, -Y, 3*np.cos(np.deg2rad(R)), 3*np.sin(np.deg2rad(R)) , color='green', linewidth=2.8)
                            
                    axs[ax1, ax2].plot(Xarr, Yarr, 'g-o', color='blue', linewidth=2)
                    
                        
                        
                    axs[ax1, ax2].set_title(filename[-20:-14])
                    imgCount = imgCount+1;
                
                for ax in axs.flat:
                    ax.set(xlabel='x', ylabel='y')
                
                # Hide x labels and tick labels for top plots and y ticks for right plots.
                for ax in axs.flat:
                    ax.label_outer()

                fout = os.path.join(outputDir,id+"_Angle_R.png")    
                fig.savefig(fout)
                
    imgCount = 0
    
    
    
   
    
    
for id in iserIdList:
    #plt.clf() 
    
    my_dpi = 100
    
    fig, axs = plt.subplots(2, 3, figsize=(40, 40), dpi=my_dpi)
    imgCount = 0
    for filename in os.listdir(modified_csv_directory): 
        
        name, ext = os.path.splitext(filename)
        
        if (( id in filename ) and (ext == '.csv') and ("Coord" in filename) and not( "init" in filename )):
            
            lineAngle = ""    
            f = os.path.join(modified_csv_directory, filename)
            if os.path.isfile(f):
                
                fileStrAngle = filename[0:-13] + "Angle*.csv"
                fileAngle = find(fileStrAngle, modified_csv_directory)[0]
                
                csvfileAngle = open(fileAngle) 
                    
                readerAngle = csv.reader(csvfileAngle)
                
                listAngle = ( list(readerAngle))
                
                tableAngle = pd.DataFrame(listAngle)
                    
                (ax1, ax2) = divmod(imgCount, 3)
                print(ax1, ax2)
                

                        
                        
                # group images
                if  ("level2" in filename):
                    axs[ax1, ax2].plot(-100, 104, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
         markersize=40, markevery=12)
                    axs[ax1, ax2].plot(-96, -90, color='red', marker='h', markerfacecolor='lightcoral', markeredgewidth=12,
         markersize=40, markevery=12)
                elif  ("level0" in filename):
                    print()
                        
                else:
                    axs[ax1, ax2].plot(-110, -121, color='green',marker='h', markerfacecolor='lightgreen', markeredgewidth=12,
         markersize=40, markevery=12)
                    axs[ax1, ax2].plot(-116, 112, color='red', marker='h', markerfacecolor='lightcoral', markeredgewidth=12,
         markersize=40, markevery=12)
                    
                
                with open(f) as csvfile:
                    
                    print(f)
                    
                    print(fileAngle)
                    
                    #plt.figure(figsize=[25,25])
                    
                    reader = csv.reader(csvfile)
                    
                    Xarr = []
                    Yarr = []
                    
                    XANGarr = []
                    YANGarr = []
                    
                    Parr = []
                    Rarr = []
                    YAarr = []
                    
                    dataset = []
                    
                    
                    for index, row in enumerate(( list(reader))):
                        str = row[0]
                           
                        coords = re.findall(r"[-+]?\d*\.*\d+", str)
                        X=float(coords[0])
                        Y=float(coords[1])
                        Z=float(coords[2])
                        Xarr.append(X)
                        Yarr.append(-Y)
                        
                        #read line for angle     
                        
                        if index < tableAngle.size:
                            str = tableAngle.iloc[index][0]
                            
                            coords = re.findall(r"[-+]?\d*\.*\d+", str)
                                                
                            #X= Roll
                            #Y= Pitch
                            #Z= Yaw
        
                            P=float(coords[0])
                            YA=float(coords[1])
                            R=float(coords[2])
                            Parr.append(P)
                            YAarr.append(Y)
                            Rarr.append(R)
                            
                            XANGarr.append( X+np.cos(np.deg2rad(P))) 
                            YANGarr.append( Y+np.sin(np.deg2rad(P))) 
                            
                            print(YA,"* in Rad = ",np.deg2rad(P))
                            
                            axs[ax1, ax2].arrow(X, -Y, 3*np.cos(np.deg2rad(P)), 3*np.sin(np.deg2rad(P)) , color='green', linewidth=2.8)
                            
                    axs[ax1, ax2].plot(Xarr, Yarr, 'g-o', color='blue', linewidth=2)
                    
                        
                        
                    axs[ax1, ax2].set_title(filename[-20:-14])
                    imgCount = imgCount+1;
                
                for ax in axs.flat:
                    ax.set(xlabel='x', ylabel='y')
                
                # Hide x labels and tick labels for top plots and y ticks for right plots.
                for ax in axs.flat:
                    ax.label_outer()

                fout = os.path.join(outputDir,id+"_Angle_P.png")    
                fig.savefig(fout)
                
    imgCount = 0