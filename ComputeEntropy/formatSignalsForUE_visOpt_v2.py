
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML
import numpy as np
import os
import random
import pandas as pd
import scipy


""""""
import numpy as np
import matplotlib.pyplot as plt

import statistics

from math import isnan

from itertools import filterfalse



statOnSignals = pd.DataFrame() 

constLevelList = ['level100', 'level110',  'vel107_2',
                      'level101', 'level104', 'level102', 
                      'level103','vel108_2', 'level107',
                      'level105', 'level106', 'level108' ]

constLevelTitles = ['Empty', 'Static', 'Per, Low Freq', 
                        'Slow, μ=1.34,σ=-0.23',   'Slow, μ=1.34,σ=0.39',       'Med, μ=0.90,σ=-0.60', 
                        'Med, μ=1.18,σ=0.55',     'Per, Low Freq, \n with noise', 'Per,  High Freq',    
                        'Fast, μ=1.18,σ=-1.34', 'Fast, μ=-0.71,σ=-1.5',      'Per, High Freq, \n with noise']

        

"""
colorLevelListBlue = ['aquamarine','cornflowerblue' ,'aqua', 'green',
                      'khaki', 'gold', 'lightsalmon', 'coral', 'palevioletred', 'maroon', 'sepia', 
                      'grey', 'grey', 'grey' ]

"""


colorUserList = ['grey', 'darkgrey', 'royalblue',
                 'darkturquoise','green', 'olive', 'gold',
                 'bisque', 'coral', 'deeppink', 
                 'crimson', 'brown',  'plum' ,
                 'magenta',
                  'burlywood', 'cadetblue', 'darkgoldenrod', 'darkolivegreen', 
                  'darkorange','darkturquoise', 'hotpink','indianred', 
                  'lightsalmon', 'maroon','violetred', 'salmon', 'salmon', 'salmon', 'salmon', 'salmon']

def legendInit(data):
    
    """
    101 - mu_1_34_sigma_m_0_23.csv
    102 - mu_0_89_sigma_m0_60.csv
    103 - mu_1_18_sigma_0_55.csv
    104 - mu_1_34_sigma_0_39.csv
    105 - res_mu_1_18_sigma_m1_34.csv
    106 - 71_sigma_m1_5_.csv
    107 - sinMotion_periods_10.csv
    108 - sinMotion_wNoise_periods_10.csv
    107_2 - sinMotion_periods_3.csv
    108_2 - sinMotion_wNoise_periods_3.csv
    """
    
    """
    constLevelList = ['level100', 'level110',  'vel107_2',
                      'level101', 'level104', 'level102', 
                      'level103','vel108_2', 'level107',
                      'level105', 'level106', 'level108' ]

    """
    
    
    
    
    
    lvlNumber = ''
    
    lvlText = '' 
    
    if 'sinMotion_wNoise_N=400 periods = 3' in filename:
        lvlNumber = '108_2'
        lvlID = 'vel108_2'
        lvlText = 'Per, Low Freq, with noise'
            
    elif 'sinMotion_wNoise_N=400 periods = 10' in filename:
        lvlNumber = '108'
        lvlID = 'level108'
        lvlText = 'Per, High Freq, with noise'
            
    elif 'sinMotion_N=400 periods = 3' in filename:
        lvlNumber = '107_2'    
        lvlID = 'vel107_2'  
        lvlText = 'Per, Low Freq'

    elif 'sinMotion_N=400 periods = 10' in filename:
        lvlNumber = '107'
        lvlID = 'level107'
        lvlText = 'Per, High Freq'
            
    elif 'mu = -0.7105263157894737, sigma = -1.5' in filename:
        lvlNumber = '106'
        lvlID = 'level106'
        lvlText = 'Fast, μ=-0.71,σ=-1.5'
                        
    elif 'mu = 1.1842105263157894, sigma = -1.3421052631578947' in filename:
        lvlNumber = '105'
        lvlID = 'level105'
        lvlText = 'Fast, μ=1.18,σ=-1.34'
    
    elif 'mu = 1.3421052631578947, sigma = 0.39473684210526305' in filename:
        lvlNumber = '104'
        lvlID = 'level104'
        lvlText = 'Slow, μ=1.34,σ=0.39'
            
    elif 'mu = 1.3421052631578947, sigma = -0.23684210526315796' in filename:
        lvlNumber = '101'
        lvlID = 'level101'
        lvlText = 'Slow, μ=1.34,σ=-0.23'
            
    elif 'mu = 1.1842105263157894, sigma = 0.5526315789473681' in filename:
        lvlNumber = '103'
        lvlID = 'level103'
        lvlText = 'Med, μ=1.18,σ=0.55'
        
    elif 'mu = 0.8999999999999999, sigma = -0.6000000000' in filename:
        lvlNumber = '102'
        lvlID = 'level102'
        lvlText = 'Med, μ=0.90,σ=-0.60'
        
        
        """
        if 'sinMotion_wNoise_N=400 periods = 3' in filename:
            legend.append('108_2')
            legendText.append('sin+Noise, slow')
            
        elif 'sinMotion_wNoise_N=400 periods = 10' in filename:
            legend.append('108')
            legendText.append('sin+Noise, fast')
            
        elif 'sinMotion_N=400 periods = 3' in filename:
            legend.append('107_2')      
            legendText.append('sin, slow')

        elif 'sinMotion_N=400 periods = 10' in filename:
            legend.append('107')
            legendText.append('sin, fast')
            
        elif 'mu = -0.7105263157894737, sigma = -1.5' in filename:
            legend.append('106')
            legendText.append('entr, fast')
                        
        elif 'mu = 1.1842105263157894, sigma = -1.3421052631578947' in filename:
            legend.append('105')
            legendText.append('entr, fast')
            
        elif 'mu = 1.3421052631578947, sigma = 0.39473684210526305' in filename:
            legend.append('104')
            legendText.append('entr, med')
            
        elif 'mu = 1.3421052631578947, sigma = -0.23684210526315796' in filename:
            legend.append('101')
            legendText.append('entr, slow')
            
        elif 'mu = 1.1842105263157894, sigma = 0.5526315789473681' in filename:
            legend.append('103')
            legendText.append('entr, med')
           
        elif 'mu = 0.8999999999999999, sigma = -0.6000000000' in filename:
            legend.append('102')
            legendText.append('entr, slow')
            
        """
    return  lvlNumber,lvlID, lvlText
        


cperiodNum = 10  # number of periods for periodic signal
Nsamles = 400  # Number of data points
cdpi = 120
legendArr = []
optSpeed = 0.05


dir = os.path.dirname(__file__)
raw_data_directory = os.path.join(dir, "selectedSets")

if not (os.path.exists(raw_data_directory)):
    os.mkdir(raw_data_directory)

raw_data_csv_directory = os.path.join(raw_data_directory, "csv")
raw_data_img_directory = os.path.join(raw_data_directory, "img")
if not (os.path.exists(raw_data_csv_directory)):
    os.mkdir(raw_data_csv_directory)
if not (os.path.exists(raw_data_img_directory)):
    os.mkdir(raw_data_img_directory)

plt.clf()
plt.figure(figsize=(20, 20), dpi=cdpi)


nSubPlots = len(os.listdir(raw_data_csv_directory))
iFile = 0


fig, axs = plt.subplots(5, 5, figsize=(32, 32), dpi=cdpi)

#fig = plt.figure(figsize=(80, 80), dpi=cdpi)

fig.set_figheight(40)
fig.set_figwidth(40)

axs[0,0] = plt.subplot2grid(shape=(5, 5), loc=(0, 0), colspan=2, rowspan=2 , title = "Path_Exp")

# axs[0,1] = plt.subplot2grid(shape=(5, 5), loc=(0, 2), colspan=2, rowspan=2, title = "Offset_Exp")
#ax3 = plt.subplot2grid(shape=(4, 4), loc=(2, 0), colspan=1, rowspan=1)
#ax4 = plt.subplot2grid((4, 4), (2, 1), colspan=1, rowspan=1)
#ax5 = plt.subplot2grid((4, 4), (2, 2), colspan=1)

for isub in range(nSubPlots):
    i = int(isub / 5)
    j = isub % 5
    print (i,j)
    axs[i,j] = plt.subplot2grid(shape=(5, 5), loc=(i, j), colspan=1, rowspan=1)


dir = os.path.dirname(__file__)
img_directory = os.path.join(dir, "plots")

if not (os.path.exists(img_directory)):
    os.mkdir(img_directory)

ivis=0
for filename in os.listdir(raw_data_csv_directory):

    lvis = len(os.listdir(raw_data_csv_directory))
    series_data = []
    resData = []

    file_path = os.path.join(raw_data_csv_directory, filename)
    name, ext = os.path.splitext(filename)
    
    if ext == ".txt":

        print(name)

        legendArr.append(name[:70])
        
        pathArr = []
        offsetArr = []
        offsetArrM = []

        with open(file_path, 'r') as file:
            for row in file:
                series_data.append(float(row.strip(' \t\n\r')))
                
        signalLen = len(series_data)
                
        sum_dx = 0
        for i, x in enumerate(series_data[0:len(series_data)-1]):
            sum_dx = sum_dx + abs(series_data[i] - series_data[i+1])

            pathArr.append(sum_dx)

            offsetArr.append(abs(series_data[i] - series_data[i+1]))
            offsetArrM.append(80*abs(series_data[i] - series_data[i+1]))
            
            
            avgSpeed = sum_dx / (len(series_data)-1)

        print(sum_dx, "avgSpeed = ", avgSpeed)
        print("offsetArr.max = ", max(offsetArr))
        print("offsetArr.min = ", min(offsetArr))
        
        print("offsetArr.mean: ", np.mean(offsetArr))
        print("offsetArr.var = ", np.var(offsetArr))
        print("offsetArr.var2 = ")
        x = abs(offsetArr - np.mean(offsetArr))**2.
        
        print(x)
        
    
        lvlID, lvlTemp ,title = legendInit(name)
                    
        df = pd.DataFrame({'name': [name],
                           'lvlID': [lvlID],
                           'title': [title],
                           'avgSpeed': [avgSpeed],
                           'mean': np.mean(offsetArr),
                           'max': [max(offsetArr)],
                           'min': [min(offsetArr)],
                           'var': np.var(offsetArr),
                           'std' : np.std(offsetArr),
                           
                           
                           'avgSpeed_m': [avgSpeed * 80],
                           'mean_m': np.mean(offsetArrM),
                           'max_m': [max(offsetArrM)],
                           'min_m': [min(offsetArrM)],
                           'var_m': np.var(offsetArrM),
                           'std_m' : np.std(offsetArrM),
                                            
                           })
                    
        statOnSignals = pd.concat([statOnSignals, df ], ignore_index=True)
                    
                    
        
        mirrorSig = 1
        if mirrorSig == 0:
            
            fullSignal = series_data
            
        else:
            rSignal = series_data.copy()
            
            if "entr" in name:
                rSignal.reverse()
            
            fullSignal = series_data + rSignal
            
        #series_data = fullSignal
        
        
        df = pd.DataFrame(fullSignal, columns=['x'])
        

        df['y'] = 0
        df['z'] = 0

        raw_data_directory_format = os.path.join(
            raw_data_directory, "formatedSets")

        if not (os.path.exists(raw_data_directory_format)):
            os.mkdir(raw_data_directory_format)

        fout = os.path.join(raw_data_directory_format, name+ '_mir_'+ str(mirrorSig) +'.csv')
        df.to_csv(fout)

        color = np.random.rand(3)
        color[0] = 1-ivis/lvis
        color[2] = ivis/lvis
        
        l=(2.2)
        if ("sin" in name) and ('10' in name) and ('wNo' in name):
            color = 'r'
            l = (4.2)
        if ("sin" in name) and ('3' in name) and not('wNo' in name):
            color = 'c'
            l = (4.2)

        #axs[0,0].plot(pathArr, c=color, linewidth=l, label='signal')

        #axs[0,1].plot(offsetArr, c=color, linewidth=l, label='signal')

        """
        i = int(iFile/4)
        j = iFile%4
        """
        
        lvlID, lvl_fromFile, lvltext = legendInit(name)
        
        
        lvli = constLevelList.index(lvl_fromFile)
        color=colorUserList[lvli]
        
        
        #axs[i,j] = plt.subplot2grid(shape=(5, 5), loc=(i, j), colspan=1, rowspan=1 , title = lvltext+ ', Level ' + str(lvlID))

        
        #axs[i,j].plot(series_data, c=color, linewidth=(4.2), label='signal')
        
        
        # >>>
        
              
        
        if 'sinMotion_wNoise_N=400 periods = 3' in filename:
            i, j = 0,1
                
        elif 'sinMotion_wNoise_N=400 periods = 10' in filename:
            i, j = 1,1
                
        elif 'sinMotion_N=400 periods = 3' in filename:
            i, j = 0,0
    
        elif 'sinMotion_N=400 periods = 10' in filename:
            i, j = 1,0
                
        elif 'mu = -0.7105263157894737, sigma = -1.5' in filename:
            i, j = 1,4
                            
        elif 'mu = 1.1842105263157894, sigma = -1.3421052631578947' in filename:
            i, j = 0,4
        
        elif 'mu = 1.3421052631578947, sigma = 0.39473684210526305' in filename:
            i, j = 1,2
                
        elif 'mu = 1.3421052631578947, sigma = -0.23684210526315796' in filename:
            i, j = 0,2
                
        elif 'mu = 1.1842105263157894, sigma = 0.5526315789473681' in filename:
            i, j = 1,3
            
        elif 'mu = 0.8999999999999999, sigma = -0.6000000000' in filename:
            i, j = 0,3
                
    
        # <<<
        
        
        #axs[i,j] = plt.subplot2grid(shape=(5, 5), loc=(i, j), colspan=1, rowspan=1 , title = lvltext+ ', Level ' + str(lvlID))
        axs[i,j] = plt.subplot2grid(shape=(5, 5), loc=(i, j), colspan=1, rowspan=1)
    
        axs[i,j].plot(series_data, c=color, linewidth=(5.2), label='signal')
        
        
        #ax.set_title('Manual ticks')
        #axs.set_yticks(np.arange(0, 1, 1/5))
        yticks = np.arange(0, 1.1, 0.25)
        #xlabels = [f'\\{x:1.2f}' for x in xticks]
        ylabels = [f'{(y*2-1):1.1f}' for y in yticks]
        axs[i,j].set_yticks(yticks, labels=ylabels,fontsize=24)
        
        xticks = np.arange(0, 440, 80)
        #xlabels = [f'\\{x:1.2f}' for x in xticks]
        xlabels = [f'{x/40:1.0f}' for x in xticks]
        axs[i,j].set_xticks(xticks, labels=xlabels,fontsize=24)
        
        
        
        
        iFile = iFile+1
    ivis = ivis+1


statOnSignals = statOnSignals.sort_values(by=['avgSpeed_m']).reset_index()


fig.legend(legendArr, ncol=1, loc='upper center')
foutImg = os.path.join(img_directory, "path_Exp_v2.pdf")

plt.savefig(f'{foutImg}', dpi=cdpi)


fout = os.path.join(raw_data_directory_format,  'speedStat.xlsx')
statOnSignals.to_excel(fout,  index = False, header=True) 







materials = np.array(statOnSignals['title']) 
x_pos = np.arange(len(materials))
CTEs = np.array(statOnSignals['avgSpeed_m']) 
error =np.array(statOnSignals['std_m']) 



# Build the plot

plt.clf()
plt.figure(figsize=(50,10))
fig, ax = plt.subplots()


ax.bar(x_pos, CTEs, yerr=error, align='center', alpha=0.5, ecolor='black', capsize=10)
#ax.set_ylabel('Coefficient of Thermal Expansion ($\degree C^{-1}$)')
ax.set_xticks(x_pos)




materials = ["Per, Low Freq", 
                        "Slow, \nμ=1.34,σ=-0.23",   "Slow, \nμ=1.34,σ=0.39",       "Med, \nμ=0.90,σ=-0.60", 
                        "Med, \nμ=1.18,σ=0.55",     "Per, Low Freq, \nwith noise", "Per,  High Freq",    
                        "Fast, \nμ=1.18,σ=-1.34", "Fast, \nμ=-0.71,σ=-1.5",      "Per, High Freq, \nwith noise"]


ax.set_xticklabels(materials)

plt.xticks(fontsize=7, rotation=45)
         

"""
def deg2rad(x_pos):
    return materials[x_pos]


def rad2deg(x_pos):
    return x_pos+1


secax = ax.secondary_xaxis('top', functions=(deg2rad, rad2deg))
"""
#plt.xticks(fontsize=7, rotation=45)

ax.set_title('Average speed of obstacle motion for 10 various conditions (in meters)')
ax.yaxis.grid(True)

# Save the figure and show
plt.tight_layout()

foutImg = os.path.join(img_directory, "bar_plot_with_error_bars.png")

plt.savefig(f'{foutImg}', dpi=cdpi)



plt.savefig('bar_plot_with_error_bars.png')

