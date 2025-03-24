from sampen import sampen2
import os

from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import EntropyHub as EH

import seaborn as sns
    

colorLevelListBlue = ['lightblue', 'royalblue', 'aquamarine','aqua', 'cornflowerblue' ,'darkturquoise',
                      'khaki', 'gold', 'lightsalmon', 'coral', 'palevioletred', 'maroon', 'sepia', 
                      'grey', 'grey', 'grey' ]

cdpi = 100
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
lmaxM = 60


tauStep = cmaxRad/cnumTauSteps


def SampEnTotal(data):
    
    columns = []
    
    series_data = data
    lsampleSize = len(series_data)
  
    jN= int( cmaxRad / tauStep)
    
    print(lmaxM, lsampleSize)
    
    s = (lsampleSize,jN)
    resEnt = np.zeros(s )
    resDev = np.zeros(s )
    
    
    tauArr = np.arange(0, cmaxRad, tauStep)
    
    
    for j, ltau in enumerate(tauArr):
        print("j = " , j)
        sampEn = sampen2( series_data, mm =lmaxM, r = ltau, normalize=False)
        print( ltau) 
        for isamp in sampEn:
            i =isamp[0]
            resEnt[i,j] = isamp[1]
            resDev[i,j] = isamp[2]
            columns.append( ( "mm = "+ str(i) + ", tau = " + str(ltau)))
            
            #print(res)
    

    mmArr = np.arange(0, lmaxM, 1)
    
    return resEnt[:lmaxM], resDev[:lmaxM], tauArr, mmArr
    



def SampEn(Sig, m=2, tau=1, r=None, Logx=np.exp(1)):
    """
        m - Embedding Dimension
        tau - Time Delay, !Default = 1
        r - Radius Distance Threshold, !Default = 0.2
        Logx - Logarithm base, !Default ln
    """
    Sig = np.squeeze(Sig)
    N = Sig.shape[0]  
    if r is None:
        r = 0.2*np.std(Sig)
  
    assert N>10 and Sig.ndim == 1
    assert isinstance(m,int) and (m > 0)
    assert isinstance(tau,int) and (tau > 0)
    assert isinstance(r,(int,float)) and (r>=0)
    assert isinstance(Logx,(int,float)) and (Logx>0)
    
    Counter = (abs(np.expand_dims(Sig,axis=1)-np.expand_dims(Sig,axis=0))<= r)*np.triu(np.ones((N,N)),1)  
    M = np.hstack((m*np.ones(N-m*tau), np.repeat(np.arange(m-1,0,-1),tau)))
    A = np.zeros(m + 1)
    B = np.zeros(m + 1)
    A[0] = np.sum(Counter)
    B[0] = N*(N-1)/2
    
    for n in range(M.shape[0]):
        ix = np.where(Counter[n, :] == 1)[0]
        
        for k in range(1,int(M[n]+1)):              
            ix = ix[ix + (k*tau) < N]
            p1 = np.tile(Sig[n: n+1+(tau*k):tau], (ix.shape[0], 1))                       
            p2 = Sig[np.expand_dims(ix,axis=1) + np.arange(0,(k*tau)+1,tau)]
            ix = ix[np.amax(abs(p1 - p2), axis=1) <= r] 
            if ix.shape[0]:
                Counter[n, ix] += 1
            else:
                break
    
    for k in range(1, m+1):
        A[k] = np.sum(Counter > k)
        B[k] = np.sum(Counter[:,:-(k*tau)] >= k)
        print("m, k ", m, k)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        Samp = -np.log(A/B)/np.log(Logx)
        print("A = ", A)
        print("B = ", B)
        print("A/B = ", A/B)
    return Samp, A, B



dir = os.path.dirname(__file__)
raw_data_csv_directory = os.path.join(dir, "samples")
img_directory = os.path.join(dir, "plots")

entrTotal = []
legendArr = []


for filename in os.listdir(raw_data_csv_directory): 
    
    series_data = []
    
    file_path = os.path.join(raw_data_csv_directory, filename)
    name, ext = os.path.splitext(filename)
    
    print(name)
    
    with open(file_path, 'r') as file:
        for row in file:
            series_data.append(float(row.strip(' \t\n\r')))
            
    
    
    plt.clf()
    plt.figure(figsize=(10,10), dpi= cdpi)
    fig, axs = plt.subplots(1, 4, figsize=(40, 10), dpi=cdpi)
    
          
    axs[0].plot( series_data, linewidth= (2.2 ), label= 'signal')
    
    
    
    
    
            
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12});    
    
    
    plt.title(name)
    plt.xlabel('Level')
    plt.ylabel('X value')
            
        
    if not (os.path.exists(img_directory)) :
        os.mkdir(img_directory)
        
    foutImg = os.path.join(img_directory, str(name) + ".png")      
    
    plt.savefig(f'{foutImg}', dpi=cdpi)
    
    plt.clf()

    #Tau is a second soordinate
    resEnt, resDev, tauArr, mmArr = SampEnTotal(series_data)
    
    #r = resEnt[5:-5,:-1].reshape(resEnt[5:-5,:-1].size)
    r = resEnt.reshape(resEnt.size)
    
    labelsArr = []
    
    for mm in mmArr:
        for tau in tauArr:
            labelsArr.append( ("mm = "+ str(mm) + ", tau = "+ str(tau)) )
    
    if len(entrTotal) ==0:
        entrTotal = np.array(r)
    else:
        entrTotal = np.vstack([np.array(entrTotal),np.array(r)])
    
    legendArr.append(name)
    
    
    startInd = 2

    sns.set()
    plt.title("Entropy heatmap")
    plt.figure(figsize=(20,20), dpi= cdpi)
    
    
    xticklabels = tauArr[startInd:]
    yticklabels = mmArr[startInd:]
    
    ax = sns.heatmap(resEnt[startInd:,startInd:], xticklabels=np.around(xticklabels,3), annot=False)
    
      
    foutImg = os.path.join(img_directory, str(name) + "_entropy.png")      
    foutcsv = os.path.join(img_directory, str(name) + "_entropy.csv")      
    
    plt.savefig(f'{foutImg}', dpi=cdpi)
    
    dfEntr = pd.DataFrame(resEnt)
    dfEntr.columns = tauArr
    dfEntr['mm'] = mmArr
    dfEntr.to_csv(foutcsv)
    
    
    foutStat = os.path.join(img_directory, str(name) + "_entropy.xlsx") 
    dfEntr.to_excel(foutStat,  index = False, header=True) 
    
    plt.clf()
    
    sns.set()
    plt.title("Standard deviation")
    plt.figure(figsize=(20,20), dpi= cdpi)
    
    ax = sns.heatmap(resDev[startInd:,startInd:], xticklabels=np.around(xticklabels,3), annot=True)
    
      
    foutImg = os.path.join(img_directory, str(name) + "_deviation.png")      
    foutcsv = os.path.join(img_directory, str(name) + "_deviation.csv")      
    
    plt.savefig(f'{foutImg}', dpi=cdpi)
    
    pd.DataFrame(resDev).to_csv(foutcsv)
    
    plt.clf()
    
    
    np.save(os.path.join(img_directory, 'tauArr.npy'), tauArr)
    
    np.save(os.path.join(img_directory, 'mmArr.npy') , mmArr)
    
    np.save(os.path.join(img_directory, 'legendArr.npy'), legendArr)
    
    
plt.clf()

plt.figure(figsize=(60,20), dpi= cdpi)
dim = entrTotal.shape
visSize = int(dim[1]/ 4)
    

visN = 20
visSize = int(dim[1]/ visN)

mmN = mmArr.size
tauN = tauArr.size
    


for t in np.arange(visN):
        
       
    plt.clf()
    plt.figure(figsize=(60,20), dpi= cdpi)
    
    for i, entrVal in enumerate(entrTotal):
    
        sigTitle = legendArr[i]
        if 'sin' in sigTitle:
            lstyle= 'dotted'
        else:
            lstyle= 'solid'
            
        #if math.isnan(entrVal) == False:
            
        startInd = visSize*t
        endInd = visSize*(t+1)
        
        plt.plot(entrVal[startInd:endInd], c=colorLevelListBlue[i],linestyle= lstyle)
            
        plt.legend(legendArr, ncol=1, loc='upper center')
        
        x_ticks = np.arange(len(labelsArr[startInd:endInd]))
        plt.xticks(ticks = x_ticks, labels = labelsArr[startInd:endInd])
         
        plt.xticks(fontsize=6, rotation=90)
            
        foutImg = os.path.join(img_directory, "total_entropy_part"+str(t)+".png")     
          
        plt.ylim(0, 0.8)      
        plt.savefig(f'{foutImg}', dpi=cdpi)
    





plt.clf()
plt.figure(figsize=(20,20), dpi= cdpi)

       
for i, entrVal in enumerate(entrTotal):

    r = entrVal.reshape(mmN, tauN)
    
    
    for indmm, mm in enumerate(mmArr):
    
        sigTitle = legendArr[i]
        
        x= tauArr
        y = r[indmm, :]
        
        plt.scatter( x, y, s = 20, c=colorLevelListBlue[i],alpha=0.5)
            
        plt.legend(legendArr, ncol=1, loc='upper center')
        
        x_ticks = np.arange(len(tauArr))
        plt.xticks(ticks = x_ticks, labels = tauArr)
         
        plt.xticks(fontsize=6, rotation=90)
                
        plt.xlabel("tau")
        plt.ylabel("SampleEn")
            
        foutImg = os.path.join(img_directory, "total_entropy_sum_on_MM.png")     
          
        plt.ylim(0, 0.6)      
        plt.savefig(f'{foutImg}', dpi=cdpi)
    

plt.clf()
plt.figure(figsize=(30,30), dpi= cdpi)


       
for i, entrVal in enumerate(entrTotal):

    r = entrVal.reshape(mmN, tauN)
    
    
    for indtau, tau in enumerate(tauArr):
    
        sigTitle = legendArr[i]
            
        
        x = mmArr
        y = r[:, indtau]
        
        plt.scatter(x, y, s = 40, c=colorLevelListBlue[i],alpha=0.5)
            
        #plt.legend(legendArr, ncol=1, loc='upper center')
        
x_ticks = np.arange(len(mmArr))
plt.xticks(ticks = x_ticks, labels = mmArr)
         
plt.xticks(fontsize=6, rotation=90)

#ax.legend(legendArr, ncol=1, loc='upper center')
            
plt.xlabel("MM")
plt.ylabel("SampleEn")
            
foutImg = os.path.join(img_directory, "total_entropy_sum_on_tau.png")     
          
plt.ylim(0, 0.8)   


plt.savefig(f'{foutImg}', dpi=cdpi)
    
    
    

"""
for filename in os.listdir(raw_data_csv_directory): 
    
    series_data = []
    
    file_path = os.path.join(raw_data_csv_directory, filename)
    name, ext = os.path.splitext(filename)
    
    print(name)
    
    with open(file_path, 'r') as file:
        for row in file:
            series_data.append(float(row.strip(' \t\n\r')))
            
    
    
    plt.clf()
    plt.figure(figsize=(10,10), dpi= cdpi)
    fig, axs = plt.subplots(1, 4, figsize=(40, 10), dpi=cdpi)
    
          
    axs[0].plot( series_data, linewidth= (2.2 ), label= 'signal')
            
    
    lgnd  = plt.legend( loc='lower center' , ncol = 2, frameon=True, prop={'size': 12});    
    
    
    plt.title(name)
    plt.xlabel('Level')
    plt.ylabel('X value')
            
        
    if not (os.path.exists(img_directory)) :
        os.mkdir(img_directory)
        
    foutImg = os.path.join(img_directory, str(name) + ".png")      
    
    plt.savefig(f'{foutImg}', dpi=cdpi)
    
    plt.clf()
    
"""




"""
plt.clf()
plt.figure(figsize=(60,20), dpi= cdpi)

       
for i, entrVal in enumerate(entrTotal):

    r = entrVal.reshape(mmN, tauN)
     
    plt.scatter(entrVal[:, :], c=colorLevelListBlue[i],alpha=0.5)
            
    plt.legend(legendArr, ncol=1, loc='upper center')
        
    x_ticks = mmArr
    plt.xticks(ticks = x_ticks, labels = mmArr)
         
    plt.xticks(fontsize=6, rotation=90)
            
    foutImg = os.path.join(img_directory, "total_entropy_sum_on_MM.png")     
          
    plt.ylim(0, 0.8)      
    plt.savefig(f'{foutImg}', dpi=cdpi)
"""    
    
    
    

"""

plt.clf()
mu, sigma = 10, 0.1 # mean and standard deviation
s = np.random.normal(mu, sigma, 1000)
df = pd.Series(s)
p = plt.plot(df.index, df.values)

Samp = sampen2(s, mm=4, r = 0.2, normalize=False) 
print( Samp )

Samp = sampen2(s, mm=5, r = 0.2, normalize=False) 
print( Samp )


Samp = sampen2(s, mm=6, r = 0.2, normalize=False) 
print( Samp )

Samp = sampen2(s, mm=7, r = 0.2, normalize=False) 
print( Samp )
"""
