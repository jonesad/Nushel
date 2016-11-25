# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 11:48:14 2016

Plot the res from a directory containing a bunch of my fits for the 

@author: jonesad
"""
overdir = 'E:\initdir\Lodai_usda_srg24_fit\\'
from os import listdir
lFits = listdir(overdir)
xfit = [0]
yfit = [19]
for fit in lFits:
    spath = overdir+'\\'+fit + '\\tracking\\res.dat'
    dat = open(spath, 'r')
    nlc = fit.split('_')[-1]
    try:
        xfit.append(int(nlc))
    except ValueError:
        nlc=nlc[0]
        xfit.append(int(nlc))
    for line in dat:
        line = line.strip().split()
    yfit.append(float(line[-1]))
xfitn = sorted(xfit)
yfitn = [y for (x,y) in sorted(zip(xfit,yfit))]
from matplotlib import pyplot as plt
yold = [.6, .38, .28, .21,  .19, .170, .165, .1625, .16, .15, .13]
xold = [5,  10,  15,  20,   25,  30,   35,   40,    45,   50,  56]
plt.plot(xold, yold, marker = '.', color = 'k', label = 'usd')
plt.plot(xfitn, yfitn, marker = '+', label = 'srg24')
plt.ylim([0, .7])
plt.title('Errors as a function of fitted linear combinations')
plt.ylabel('RMSE (MeV)')
plt.xlabel('Number of linear Combinations')
plt.legend(loc = 'best')
plt.savefig('E:\\initdir\\Lodai_usda_srg24_fit\\srg24_iterations.pdf')