# -*- coding: utf-8 -*-
"""
Created on Thu Nov 05 14:07:27 2015

@author: jonesad
"""

import numpy as np
a=np.random.rand(200, 50)
b=np.random.rand(200,1)
xls= np.linalg.lstsq(a,b)
xls=xls[0]
[u,s,vt]=np.linalg.svd(a)
Sd=np.zeros((vt.shape[0], u.shape[1]))
err = []
for i in range(s.size):
    for nidx in range(i + 1):
        Sd[nidx,nidx]=1/s[nidx]
    api=np.dot(np.dot(np.transpose(vt),Sd), np.transpose(u))
    xsvd=np.dot(api,b)
    err.append(np.linalg.norm(np.subtract(xls,xsvd))/np.linalg.norm(xls))
import matplotlib.pyplot as plt
plt.xkcd()
print s
plt.plot(err)
plt.show()