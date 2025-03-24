from sampen import sampen2
import os

from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import EntropyHub as EH

import seaborn as sns
    

"""
colorLevelListBlue = ['lightblue', 'royalblue', 'aquamarine','cornflowerblue' ,'aqua', 'green',
                      'khaki', 'gold', 'lightsalmon', 'coral', 'palevioletred', 'maroon', 'sepia', 
                      'grey', 'grey', 'grey' ]

sortedLvlList = [ '107_2' ,'108_2' , '107' ,'108' , '101' , '102','103', '104', '105','106']
legTextSorted = ['sin, slow',  'sin+Noise, slow' ,'sin, fast' ,  'sin+Noise, fast',
                 'entr, slow','entr, slow','entr, med', 'entr, med'  , 'entr, fast' , 'entr, fast']
"""


colorUserList = [ 'royalblue',
                 'darkturquoise','green', 'olive', 'gold',
                 'bisque', 'coral', 'deeppink', 
                 'crimson', 'brown',  'plum' ,
                 'magenta',
                  'burlywood', 'cadetblue', 'darkgoldenrod', 'darkolivegreen', 
                  'darkorange','darkturquoise', 'hotpink','indianred', 
                  'lightsalmon', 'maroon','violetred', 'salmon', 'salmon', 'salmon', 'salmon', 'salmon']


sortedLvlList = ['107_2',
                      '101', '104', '102', 
                      '103','108_2', '107',
                      '105', '106', '108' ]

"""
legTextSorted = ['Per, Low Freq', 
                        'Slow, μ=1.34,σ=-0.23',   'Slow, μ=1.34,σ=0.39',       'Med, μ=0.90,σ=-0.60', 
                        'Med, μ=1.18,σ=0.55',     'Per, Low Freq, with noise', 'Per,  High Freq',    
                        'Fast, μ=1.18,σ=-1.34', 'Fast, μ=-0.71,σ=-1.5',      'Per, High Freq, with noise']

"""
legTextSorted = ['P_LF_1.2', 'Sto_1.4',   'Sto_1.9',  'Sto_2.5', 'Sto_3.1', 'PN_LF_3.8', 'P_HF_4.0', 'Sto_4.6', 'Sto_4.8', 'PN_HF_4.9']


markers = ["d", "D", "v", "s", "*", "^", "3", "4", "<", ">", "1", "2", "8", "p", "+", "P","|","_" , "o",","]


cdpi = 800
#cmaxRad = 0.22
#cnumTauSteps = 50
#lmaxM = 100

"""
cmaxRad = 0.08 #small entr
cnumTauSteps = 30
lmaxM = 60

"""


cmaxRad = 0.22
cnumTauSteps = 30
lmaxM = 150


tauStep = cmaxRad/cnumTauSteps



dir = os.path.dirname(__file__)
raw_data_csv_directory = os.path.join(dir, "entrData")
img_directory = os.path.join(dir, "plots")

legendArr = []

plt.clf()
plt.figure(figsize=(16,32), dpi= cdpi)

fig, ax = plt.subplots()



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


    constLevelTitles = ['Empty', 'Static', 'Per, Low Freq', 
                        'Slow, μ=1.34,σ=-0.23',   'Slow, μ=1.34,σ=0.39',       'Med, μ=0.90,σ=-0.60', 
                        'Med, μ=1.18,σ=0.55',     'Per, Low Freq, with noise', 'Per,  High Freq',    
                        'Fast, μ=1.18,σ=-1.34', 'Fast, μ=-0.71,σ=-1.5',      'Per, High Freq, with noise']
    """
    
    
    
    
    
    legend = []
    
    legendText = []
    for filename in data:
        
        if 'sinMotion_wNoise_N=400 periods = 3' in filename:
            legend.append('108_2')
            legendText.append('Per, Low Freq, with noise')
            
        elif 'sinMotion_wNoise_N=400 periods = 10' in filename:
            legend.append('108')
            legendText.append('Per, High Freq, with noise')
            
        elif 'sinMotion_N=400 periods = 3' in filename:
            legend.append('107_2')      
            legendText.append('Per, Low Freq')

        elif 'sinMotion_N=400 periods = 10' in filename:
            legend.append('107')
            legendText.append('Per, High Freq')
            
        elif 'mu = -0.7105263157894737, sigma = -1.5' in filename:
            legend.append('106')
            legendText.append('Fast, μ=-0.71,σ=-1.5')
                        
        elif 'mu = 1.1842105263157894, sigma = -1.3421052631578947' in filename:
            legend.append('105')
            legendText.append('Fast, μ=1.18,σ=-1.34')
            
        elif 'mu = 1.3421052631578947, sigma = 0.39473684210526305' in filename:
            legend.append('104')
            legendText.append('Slow, μ=1.34,σ=0.39')
            
        elif 'mu = 1.3421052631578947, sigma = -0.23684210526315796' in filename:
            legend.append('101')
            legendText.append('Slow, μ=1.34,σ=-0.23')
            
        elif 'mu = 1.1842105263157894, sigma = 0.5526315789473681' in filename:
            legend.append('103')
            legendText.append('Med, μ=1.18,σ=0.55')
           
        elif 'mu = 0.8999999999999999, sigma = -0.6000000000' in filename:
            legend.append('102')
            legendText.append('Med, μ=0.90,σ=-0.60')
          
            
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
    return legend, legendText


    """  
for filename in os.listdir(raw_data_csv_directory): 
    
    series_data = []
    
    file_path = os.path.join(raw_data_csv_directory, filename)
    name, ext = os.path.splitext(filename)
    
    
    entrTotal = pd.read_excel(file_path, dtype=np.float64)
    entrTotal['mm'] = entrTotal['mm'].astype('int64')
    
    legendArr = entrTotal.columns[2:].to_numpy()
    
    leg, legText = legendInit(legendArr)
    
    tauArr = np.unique(entrTotal['tau'])
    mmArr = np.unique(entrTotal['mm'])
    
    for i, level in enumerate(sortedLvlList):
        for j, tau in enumerate(tauArr):
            visEntr = entrTotal.loc[ entrTotal['tau'] == tau]
            
            ind = leg.index(level)
            lvlFile = legendArr[ind]
            
            print(level, leg[ind], lvlFile )
        
            ax.scatter(visEntr['mm'], visEntr[lvlFile], c = colorUserList[i], 
                       marker = markers[j], linewidths=0.2 , edgecolors= "white" )
   
    
    for i, level in enumerate(legendArr):
    
        ax.plot(np.NaN, np.NaN, c= colorUserList[i], label=legTextSorted[i])
    ax2 = ax.twinx()
    
    for j, tau in enumerate(tauArr):
        print(markers[j], tau)
        stau= str(tau)
        ax2.plot(np.NaN, np.NaN, marker = markers[j], label=stau[:6])
    ax2.get_yaxis().set_visible(False)
    
    ax.legend(loc='upper center', markerscale=0.75, prop={'size': 4})
    ax2.legend(loc='upper right', markerscale=0.75, prop={'size': 4})
          
    plt.xlabel("mm")
    plt.ylabel("SampleEn")
    #plt.ylim(0, 0.6)   
            
    
    foutImg = os.path.join(img_directory, filename+"total_entropy_selected_param.png")     
    
    plt.savefig(f'{foutImg}', dpi=cdpi)
            
    """
    
    

for filename in os.listdir(raw_data_csv_directory): 
    
    series_data = []
    
    file_path = os.path.join(raw_data_csv_directory, filename)
    name, ext = os.path.splitext(filename)
    
    
    entrTotal = pd.read_excel(file_path, dtype=np.float64)
    entrTotal['mm'] = entrTotal['mm'].astype('int64')
    
    legendArr = entrTotal.columns[2:].to_numpy()
    
    leg, legText = legendInit(legendArr)
    
    tauArr = np.unique(entrTotal['tau'])
    mmArr = np.unique(entrTotal['mm'])
    
    for i, level in enumerate(sortedLvlList):
        for j, tau in enumerate(tauArr):
            visEntr = entrTotal.loc[ entrTotal['tau'] == tau]
            
            ind = leg.index(level)
            lvlFile = legendArr[ind]
            
            print(level, leg[ind], lvlFile )
        
            ax.scatter(visEntr['mm'][2:30], visEntr[lvlFile][2:30], c = colorUserList[i], 
                       marker = markers[int(j/10)], linewidths=0.2 , edgecolors= "white" )
   
    
    for i, level in enumerate(legendArr):
    
        ax.plot(np.NaN, np.NaN, c= colorUserList[i], label=legTextSorted[i])
    ax2 = ax.twinx()
    
    for j in range(0,int(len(tauArr)/10)):
        stau=  str(tauArr[j*10])[:4] + " - " +str(tauArr[j*10+9])[:4]
        ax2.plot(np.NaN, np.NaN, marker = markers[j], label=stau)
    ax2.get_yaxis().set_visible(False)
    
    ax.legend(loc='upper center', markerscale=0.75, prop={'size': 4})
    ax2.legend(loc='upper right', markerscale=0.75, prop={'size': 4})

   
    
    #plt.legend(markers, tauArr, ncol=2, loc='upper center')        
    plt.xlabel("mm")
    plt.ylabel("SampleEn")
    #plt.ylim(0, 0.6)   
            
    
    foutImg = os.path.join(img_directory, filename+"__total_entropy_selected_param_02.png")     
    
    plt.savefig(f'{foutImg}', dpi=cdpi)