from sampen import sampen2
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import EntropyHub as EH


def SampEn(Sig, m=2, =1, r=None, Logx=np.exp(1)):
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
    assert isinstance(Logx,(int,ftauloat)) and (Logx>0)
    
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


series_data = []

dir = os.path.dirname(__file__)
file_path = os.path.join(dir, 'sampentest.txt')


with open(file_path, 'r') as file:
    for row in file:
        series_data.append(float(row.strip(' \t\n\r')))

#arr = [math.sin(i*math.pi/2) for i in range(1,21)]
#npArr = np.fromiter(arr, float)
#series_data = pd.Series(arr)

# calculate the sample entropy
sampen_of_series = sampen2(series_data, mm=4, r = 0.2, normalize=False) #, mm=5, normalize=True
print(sampen_of_series)


df = pd.Series(series_data)
p = plt.plot(df.index, df.values)

Samp = SampEn( series_data, m =4, r=0.2 )
print( Samp )


plt.clf()
mu, sigma = 10, 0.1 # mean and standard deviation
s = np.random.normal(mu, sigma, 1000)
df = pd.Series(s)
p = plt.plot(df.index, df.values)

Samp = SampEn( s, m =4, r=0.2 )
print( Samp )


plt.clf()
mu, sigma = 10, 0.1 # mean and standard deviation
s = np.random.normal(mu, sigma, 1000)
df = pd.Series(s)
p = plt.plot(df.index, df.values)

Samp = sampen2(s, mm=4, r = 0.2, normalize=False) 
print( Samp )

