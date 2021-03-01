# -*- coding: utf-8 -*-
'''
Author: 小吉心
Date: Unknown
LastEditors: Payne
LastEditTime: 2021-03-01 20:11:29
Description: 
'''


#!/usr/bin/env python
# coding: utf-8

# In[35]:


import numpy as np
import matplotlib.pyplot as plt
import math


file = "./hw2_fitting.dat"
data = np.loadtxt(file)
print(data.shape)


# In[36]:


xi =data[:,0]
yi =data[:,1]
sigmai =data[:,2]
print(xi.shape)
print(yi.shape)
print(sigmai.shape)


# In[37]:


def Loren(v,v0,L):
    return L/(math.pi*((v-v0)**2+L**2))

print(Loren(50,40,10))

plt.plot(xi,Loren(xi,46.269,7.545))
# show 只出现一次(绘制在一张图上)
# plt.show()

def pvLoren(v,v0,L):
    return 2*(v-v0)*L/(math.pi*((v-v0)**2+L**2)**2)

print(pvLoren(50,40,10))

def plLoren(v,v0,L):
    return ((v-v0)**2-L**2)/(math.pi*((v-v0)**2+L**2)**2)

print(plLoren(50,40,10))

def chisqLoren(v0,L):
    return np.sum(((yi-Loren(xi,v0,L))/sigmai)**2)

def pchisqLoren(v0,L):
    return np.matrix([np.sum(2*(Loren(xi,v0,L)-yi)*pvLoren(xi,v0,L)/sigmai**2),np.sum((Loren(xi,v0,L)-yi)*plLoren(xi,v0,L)/sigmai**2)])

print(chisqLoren(46.269,7.545))


# In[38]:


def exp(v,v0,D):
    return np.exp((-math.log(2))*(v-v0)**2/D**2)

def Gaus(v,v0,D):
    return exp(v,v0,D)*math.sqrt(math.log(2))/(math.sqrt(math.pi)*D)

print(Gaus(50,40,10))
plt.plot(xi,Gaus(xi,44.954,15.026))
# #####################################################
# plt.show()

def pvGaus(v,v0,D):
    return exp(v,v0,D)*(2*math.log(2)**(3/2)*(v-v0))/(D**3*math.sqrt(math.pi))

print(pvGaus(50,40,10))

def pdGaus(v,v0,D):
    return exp(v,v0,D)*math.sqrt(math.log(2))*(2*math.log(2)*(v-v0)**2-D**2)/(D**4*math.sqrt(math.pi))

print(pdGaus(50,40,10))

def chisqGaus(v0,D):
    return np.sum(((yi-Gaus(xi,v0,D))/sigmai)**2)

def pchisqGaus(v0,D):
    return np.matrix([np.sum(2*(Gaus(xi,v0,D)-yi)*pvGaus(xi,v0,D)/sigmai**2),np.sum((Gaus(xi,v0,D)-yi)*pdGaus(xi,v0,D)/sigmai**2)])

print(chisqGaus(44.954,15.026))

##########################################################
plt.errorbar(xi, yi, sigmai)
plt.show()
##########################################################


# In[39]:


def LMmethodLoren(v0,L,e):
    chisq_new = 0
    chisq_old = 1
    while abs(chisq_old - chisq_new)/chisq_old > 0.000001:
        chisq_old = chisqLoren(v0,L)
        A = np.zeros([2,2])
        A[0,0] = np.sum((pvLoren(xi,v0,L)/sigmai)**2)
        A[0,1] = np.sum(pvLoren(xi,v0,L)*plLoren(xi,v0,L)/sigmai**2)
        A[1,0] = np.sum(pvLoren(xi,v0,L)*plLoren(xi,v0,L)/sigmai**2)
        A[1,1] = np.sum((plLoren(xi,v0,L)/sigmai)**2)
        A1 = np.zeros(A.shape)
        for i in range(2):
            for j in range(2):
                if i == j:
                    A1[i,j] = (1+e)*A[i,j]
                else:
                    A1[i,j] = A[i,j]
        beta = np.transpose(-1*pchisqLoren(v0,L))
        da = np.dot(np.linalg.inv(A1),beta)
        v1 = v0 + da[0,0]
        L1 = L + da[1,0]
        chisq_new = chisqLoren(v1,L1)
        if chisqLoren(v1,L1) >= chisqLoren(v0,L):
            e = e*10
        else:
            e = e/10
            v0 = v1
            L = L1
        #print(k,"",v0,"",L,chisqLoren(v0,L))
    return (v0,L,chisqLoren(v0,L))
    
print(LMmethodLoren(50,15,0.001))


# In[40]:


def LMmethodGaus(v0,D,e):
    chisq_new = 0
    chisq_old = 1
    while abs(chisq_old - chisq_new)/chisq_old > 0.000001:
        chisq_old = chisqGaus(v0,D)
        A = np.zeros([2,2])
        A[0,0] = np.sum((pvGaus(xi,v0,D)/sigmai)**2)
        A[0,1] = np.sum(pvGaus(xi,v0,D)*pdGaus(xi,v0,D)/sigmai**2)
        A[1,0] = np.sum(pvGaus(xi,v0,D)*pdGaus(xi,v0,D)/sigmai**2)
        A[1,1] = np.sum((pdGaus(xi,v0,D)/sigmai)**2)
        A1 = np.zeros(A.shape)
        for i in range(2):
            for j in range(2):
                if i == j:
                    A1[i,j] = (1+e)*A[i,j]
                else:
                    A1[i,j] = A[i,j]
        beta = np.transpose(-1*pchisqGaus(v0,D))
        da = np.dot(np.linalg.inv(A1),beta)
        v1 = v0 + da[0,0]
        D1 = D + da[1,0]
        chisq_new = chisqGaus(v1,D1)
        if chisqGaus(v1,D1) >= chisqGaus(v0,D):
            e = e*10
        else:
            e = e/10
            v0 = v1
            D = D1
        #print(k,"",v0,"",D,GausLoren(v0,D))
    return (v0,D,chisqGaus(v0,D))
    
print(LMmethodGaus(50,15,0.001))


# In[ ]:




