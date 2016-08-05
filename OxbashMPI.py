# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 14:57:36 2016
utility file for mpi execution of oxbash shell model calculations
@author: jonesad
"""

def dstributeLoad(nSize, nNames,sMethod=''):
    '''
        Deal out work to the different processor ids nased on method
    '''
    print 'distributing load'
    import numpy as np
    npaSums = np.zeros(nSize)
    llnIdxs = [[] for i in range(nSize)]
    import os
    sCWD=os.getcwd()
    bTimes = False
    if os.path.isfile(sCWD + '\\times.dat'):
        bTimes = True
        fTimes = open(sCWD + '\\times.dat', 'r')
        lfFtimes = []
        lFnames = []
        for line in fTimes:
            line=line.split()
            line[0] = int(line[0])
            line[1] = float(line[1])           
            lFnames.append(line[0])
            lfFtimes.append(line[1])
        fTimes.close()
    if (sMethod == 'bunched' or sMethod == 'spread') and bTimes ==True:
        '''
            Put largest workloads on adjacent processes
        '''
        if sMethod == 'bunching':
            print 'bunching'
        elif sMethod == 'spread':
            print 'bunching then spreading'
        lFcopy = [elem for elem in lFnames]
        for i in range(len(lFnames)):
            nIdx = np.argmin(npaSums)
            nJIdx = np.argmax(lfFtimes)
            npaSums[nIdx] += lfFtimes[nJIdx]
            llnIdxs[nIdx].append(lFcopy[nJIdx])
            del lfFtimes[nJIdx]
            del lFcopy[nJIdx]
            print len(lfFtimes), ' left to bunch.'
        if sMethod == 'spread':
            print 'spreading'
            llnIdxsnew = []
            npaScopy = np.array(npaSums)
            newSums = np.zeros(npaSums.shape)
            for i in range(npaSums.size):
                print 'i', i, i % 2
                if i % 2 == 1:
                    nIdx = np.argmin(npaScopy)
                else:
                    nIdx = np.argmax(npaScopy)
                newSums[i] = npaScopy[nIdx]
                llnIdxsnew.append(llnIdxs[nIdx])
                npaScopy = np.delete(npaScopy, nIdx)
                del llnIdxs[nIdx]
            npaSums = np.array(newSums)
            llnIdxs = [elem for elem in llnIdxsnew]
            print 'spreading  complete'               
    if sMethod =='' or bTimes == False:
        '''
            Deal out the list of names in the order they come
        '''
        print 'dealing'
        for nIdx in range(nSize):
            print 'nIdx= ', nIdx
            llnIdxs[nIdx] = list(range(nIdx, nNames, nSize))
            if bTimes:
                for nJIdx in llnIdxs[nIdx]:
                    npaSums[nIdx] += lfFtimes[nJIdx]
        print 'dealing done'
    if bTimes:
        fSF = open('SumFile.dat', 'a')
        fSF.write(str(npaSums) + '\n')
        fSF.close()
    print 'distribution complete'
    return llnIdxs


from mpi4py import MPI
oComm = MPI.COMM_WORLD
nRank = oComm.rank
oComm.Barrier()
#read input data from a file sitting in the working directory
lNames = [] 
import os
sCWD=os.getcwd()
nSize=oComm.size
llnIdxs = [[] for i in range(nSize)]
if nRank == 0:
    if not os.path.isfile('MPINames.dat'):
        print 'Error: input file MPINames.dat is not present in current directory:'
        print os.getcwd()
    else:
        fMPINames = open('MPINames.dat', 'r')
        for line in fMPINames:
            lNames.append(line.strip())
        fMPINames.close()
    llnIdxs = dstributeLoad(nSize, len(lNames), sMethod='spread')
    print 'continuing execution'
#    print lNames
llnIdxs=oComm.bcast(llnIdxs, root=0)
lNames=oComm.bcast(lNames, root=0)
mylist=[]
import time
for i in llnIdxs[nRank]:
    start = time.time()
    os.chdir(sCWD + '\\' + lNames[i])
    os.system('shell ' + lNames[i] + '.ans')
    os.system(lNames[i])
    end = time.time()
    total = end - start
    mylist.append([i, total])
mylist = oComm.gather(mylist, root=0)
if nRank == 0:
    newlist = []
    for elem in mylist:
        newlist.extend(elem)
    fTimes = open(sCWD + '\\times.dat','w')
    sform = '{:5d}{:10.2f}\n'
    for elem in newlist:
        fTimes.write(sform.format(elem[0], elem[1]))
    fTimes.close()
oComm.Barrier()
print nRank, 'made it past barrier'
MPI.Finalize()