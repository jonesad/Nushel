# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 16:47:49 2015

class for managing a single nucleus

@author: jonesad
"""
import sys
sys.path.append('C:\\PythonScripts\\OxBahshScripts')
import OxbashOptFl


class nucleus(OxbashOptFl.MEhandler):
  '''
        class for managing a single nucleus
  '''

# take the path and create a sub directory for the nucleus and run a default
# calculation for the nucleus.
  def __init__(self,  nZ, nA, nVal, fGSE, fError, sPath, sMMDR, sPar, lsShared,
               llMESpec, useGS, sForm, initialize=True, bExtrap=True,
               llStateSpec=[], sOBDir='c:\\oxbash'):
    import os
    self.sOBDir = sOBDir
    self.sPath = sPath
    self.nAZ = [nA, nZ]
    self.fGSE = fGSE
    self.fError = fError
    self.useGS = useGS
    self.bExtrap = bExtrap
    temp = self.nAZ[0]-2*self.nAZ[1]
    if temp % 2 == 0:
        self.sIsospin = str(int(temp/2))
    else:
        self.sIsospin = str(temp) + '/2'
    self.nVal = nVal
    self.sName = 'A' + str(nA) + '_Z' + str(nZ)
    self.sMS = lsShared[0]
    self.sInt = lsShared[2]
    if not os.path.exists(sPath + '\\' + self.sName):
        os.makedirs(self.sPath + '\\' + self.sName)
    if not os.path.exists(sPath + '\\' + self.sName + '\\tracking'):
        os.makedirs(self.sPath + '\\' + self.sName + '\\tracking')
    self.llMESpec = llMESpec
    self.mllspec = llStateSpec
    self.writeAns(sMMDR, sPar, lsShared)
    self.sInt = lsShared[2]
    # copy the interaction to the working directory
    import shutil
    shutil.copyfile(sOBDir + '\\sps\\' + self.sInt+'.int', self.sPath + '\\' +
                    self.sName + '\\' + self.sInt + '.int')
    if initialize:
        self.runSM()
    self.writeStatus()
# initialize the single particle matrix elements
    if len(llMESpec[0]) == 0:
        self.llMESpec[0] = range(1, self.countOBME() + 1)
    self.sForm = sForm

# make a '.ans' file for use with Nushellx
  def writeAns(self, sMMDR, sPar, lsShared):
    fAns=open(self.sPath+'\\'+self.sName+'\\'+self.sName+'.ans','w')
    sForm1='{:<15}{:<12d}{:<7d}{:<15.7f}{:<15.7f}{:<15.7f}'
    fAns.write(sForm1.format('lpe', 0, 0, 0, 0, 0)+'\n')
    fAns.write(self.sMS + '\n')
    sForm2 = '{:>12d}'
    fAns.write(sForm2.format(self.nVal)+ '\n')
    fAns.write(lsShared[1] + '\n')
    fAns.write(self.sInt +'\n')
    sForm3 = '{:>17.7f}{:>17.7f}'
    temp=sMMDR.strip().split()
    for nIdx, elem in enumerate(temp):
        temp[nIdx] = float(elem.strip(','))
    fAns.write(sForm3.format(temp[0],temp[1]) + '\n')
    fAns.write(sForm3.format(eval(self.sIsospin + '.'), 0) + '\n')
    fAns.write(sForm3.format(0, 0) + '\n')
    sForm4 = '{:>12}'
    fAns.write(sForm4.format(sPar) + '\n')
    sDenForm = '{:<15s}{:<12d}{:<7d}' + '{:<15.7f}'*3 + '\n'
    lnZero = [0]*5
    sDenLine = sDenForm.format('den', *lnZero)
    s_nJForm = '{:>12d}'*2 + '\n'
    sDecForm = '{:>15.7f}'*2 + '\n'
    for lev in self.mllspec:
        fAns.write(sDenLine)
        fAns.write('at\n')
        sEnName = self.makeEnergyName(lev[0], self.sIsospin, lev[2])
        s_nJLine = sOneForm.format(lev[1], lev[1])
        for nIdx in range(2):
            fAns.write(sEnName + '\n')
            fAns.write(s_nJLine)
        fJ = float(eval(lev[0] + '.0'))
        fTz = float(eval(self.sIsospin + '.0'))
        if lev[2] == '+1':
            fPi = 0.0
        elif lev[2] == '-1':
            fPi = 1.0
        else:
            print 'Error invalid parity: ', lev[2], 'only "+1" or "-1" allowed.'
        for nIdx in range(2):
            fAns.write(sDecForm.format(fJ, 0))            
            fAns.write(sDecForm.format(fTz, 0))            
            fAns.write(sDecForm.format(fPi, 0))
        fAns.write('y\n')
        for nIdx in range(2):
            fAns.write(sDecForm.format(0,0))
    fAns.write(sForm1.format('st', 0, 0, 0, 0, 0) + '\n\n')        
    fAns.close()

  def writeStatus(self):    
    fOut=open(self.sPath+'\\'+self.sName+'\\'+'tracking'+'\\energy.dat','a+')
    string=''
    for energy in self.getEnNu():
      string=string+str(energy)+'\t'    
    fOut.write(string+'\n')
    fOut.close()
    sLevName=self.getLevName()
    string=''
    fOut=open(self.sPath+'\\'+self.sName+'\\'+'tracking'+'\\occupation.dat','a+')
    for occ in self.getOcc(sLevName):
      string=string+str(occ)+'\t'   
    fOut.write(string+'\n')
    fOut.close()
    
  def getLevName(self):
    name = 'A' + str(self.nAZ[0]) + '_Z' + str(self.nAZ[1]) + '.lpt'
    return name

# Get Nushell energies
  def getEnNu(self,bAll=False):
    sLevName=self.getLevName()
    afETh=[]
    lsNewList=[]
    nIdx=0
    sFPath=self.sPath+'\\'+self.sName+"\\"+sLevName
    fTh=open(sFPath)
    for line in fTh:
      line=line.strip().split()
      if bAll==False:
        for lev in self.mllspec:
          if len(line)==7 and lev[0]==line[4] and lev[1]==line[1] and lev[2]==line[6]:
            if self.useGS==0:
              afETh.append(float(line[3]))
            elif self.useGS==1:
              afETh.append(float(line[2]))
            break
          #ignore first few lines and unbound states
      elif bAll==True and len(line)>=6 and nIdx>5 and float(line[2])<0.:
        lsNewList.append(line[4]+' '+line[1]+' '+line[6])
      nIdx+=1
#    change level spec to all levels and then get all levels on the list   
    if bAll==True:
      self.mllspec=list(lsNewList)
      afETh=self.getEnNu(bAll=False)
            
    if len(afETh)!=len(self.mllspec)and bAll==False:
      print "Error: # of Theory levels found does not match requested # in nAZ=", self.nAZ
      print "Requested:", len(self.mllspec)
      print "Found:", len(afETh)
    fTh.close()
    return afETh 

# run the shell model calculation        
  def runSM(self):
    import os
    os.chdir(self.sPath+'\\'+self.sName)
    print self.sName + '.ans'    
    os.system('shell '+self.sName+'.ans')
    os.system(self.sName)
    self.writeEnergies()
    self.writeOcc()

#monopole term calculation works with p-n formalism matrix element labels
  def calcMonoOcc(self):
    sLevName=self.getLevName()
    npaLabel=self.getLabel()
    npaOcc=self.getOcc(sLevName)
    
    import numpy as np      
    nMonoSize=self.getMonoME().size    
    npaMono=np.zeros([npaOcc.shape[0],nMonoSize])
    npaMonoLabel=self.getMonoLabel()
    denom=np.zeros(npaMono.shape)
    jmax=0
    for nLevIdx in range(npaMono.shape[0]):    
      nMono=0
      for nIdx in range(npaLabel.shape[0]):        
        if np.all(npaLabel[nIdx]==npaMonoLabel[nMono]):
          temp=float(2*(npaLabel[nIdx,4]+1))
#          test to see what the largest angular momentum is
          if float(npaLabel[nIdx,4])>jmax:
            jmax=float(npaLabel[nIdx,4])
          if self.sForm =='pn':
            nSPEIdx=int(npaLabel[nIdx,0])-1
          elif self.sForm =='iso':
            nSPEIdx=int(npaLabel[nIdx,0])-1+npaLabel[nIdx,-1]*self.countOBME()
          if npaLabel[nIdx,0]!=npaLabel[nIdx,1]:
            npaMono[nLevIdx,nMono]+=npaOcc[nLevIdx,nSPEIdx]*npaOcc[nLevIdx,nSPEIdx]*temp
          else:            
            npaMono[nLevIdx,nMono]+=npaOcc[nLevIdx,nSPEIdx]*(npaOcc[nLevIdx,nSPEIdx]-1.0)*temp/2.0
          denom[nLevIdx,nMono]+=temp
          nMono+=1
          if nMono>=npaMonoLabel.shape[0]:
            break
    npaMono=np.divide(npaMono,denom)
    return npaMono    

    #monopole term calculation works with p-n formalism matrix element labels
  def summedMO(self):
    sLevName=self.getLevName()
    npaLabel=self.getLabel()
    npaOcc=self.getOcc(sLevName)
    tempocc=[]
    import numpy as np
    nMonoSize=self.getMonoME().size    
    npaMono=np.zeros([npaOcc.shape[0],nMonoSize])
    npaMonoLabel=self.getMonoLabel()
    denom=np.zeros(npaMono.shape)
    for nLevIdx in range(npaMono.shape[0]):    
      nMono=0
      for nIdx in range(npaLabel.shape[0]):        
        if np.all(npaLabel[nIdx]==npaMonoLabel[nMono]):
          temp=float(2*(npaLabel[nIdx,4]+1))
          if npaLabel[nIdx,0]!=npaLabel[nIdx,1]:
            npaMono[nLevIdx,nMono]+=npaOcc[nLevIdx,int(npaLabel[nIdx,0]-1)]*npaOcc[nLevIdx,int(npaLabel[nIdx,1]-1)]*temp
          else:            
            npaMono[nLevIdx,nMono]+=npaOcc[nLevIdx,int(npaLabel[nIdx,0]-1)]*(npaOcc[nLevIdx,int(npaLabel[nIdx,1]-1)]-1.0)*temp/2.0
          denom[nLevIdx,nMono]+=temp
          nMono+=1
          if nMono>=npaMonoLabel.shape[0]:
            break
    npaMono=np.divide(npaMono,denom)
    for nIdx in self.llMESpec[0]:
      if len(tempocc) != 0:
        temp=np.array(npaOcc[:,nIdx-1])
        temp.shape=[temp.size,1]
        tempocc=np.append(tempocc, temp, axis=1)
      elif len(tempocc)==0:
        tempocc=np.array(npaOcc[:,nIdx-1])
        tempocc.shape=[tempocc.size,1]
      else:
        print 'logic err'
    if tempocc!=[]:
      npaOcc=tempocc
    npaNewLabels=[]     
    for nIdx in range(npaMonoLabel.shape[0]):
      bIsIt=False
      try:
        temp=npaNewLabels.shape[0]
      except:
        temp=len(npaNewLabels)  
      for nIdx2 in range(temp):
        if np.all(np.array([npaMonoLabel[nIdx,:4]])==npaNewLabels[nIdx2]):
          bIsIt=True
      if (not bIsIt):
        if npaNewLabels!=[]:
          npaNewLabels=np.append(npaNewLabels, [npaMonoLabel[nIdx,:4]], axis=0)          
        elif npaNewLabels==[]:
          npaNewLabels=np.array([npaMonoLabel[nIdx,:4]])
    npaNewMO=np.zeros([npaMono.shape[0], npaNewLabels.shape[0]])
    for nIdx1 in range(npaNewMO.shape[0]):
      for nIdx2 in range(npaNewLabels.shape[0]):
        for nIdx3 in range(npaMonoLabel.shape[0]):  
          if np.all(npaNewLabels[nIdx2]==npaMonoLabel[nIdx3,:4]):
            npaNewMO[nIdx1,nIdx2]+=npaMono[nIdx1][nIdx3]
    temp=np.append(npaOcc,npaNewMO,axis=1)
    return temp, npaNewLabels

#get the occupation numbers   
  def getOcc(self, sLevName):
      import numpy as np
      fIn=open(self.sPath+'\\'+self.sName+'\\'+sLevName[:-4]+'.occ', 'r')
      npaOcc=[]
      nIdx=0
      for line in fIn:
        line=line.strip().split()
        try: 
            int(line[3])
        except:
            continue
        for nlevIdx,lev in enumerate(self.mllspec):
          if (int(line[3]) == int(2 * float(eval(lev[0]+'.0'))) and
              int(line[1]) == int(lev[1]) and line[4] == lev[2]):
            temp = line[5:5 + self.countOBME()]
            temp=[temp]
            if npaOcc!=[]:              
              npaOcc=np.append(npaOcc, temp, axis=0)              
            else:
              npaOcc=temp
      fIn.close()
      
      return np.array(npaOcc,dtype=float)
#From the occupations listed return only the ones associated with SPE on the 
#llMESpec[0] list     
  def getReducedOcc(self):
    import numpy as np
    tempocc=[]
    npaOcc=self.getOcc(self.getLevName())
    for nIdx in self.llMESpec[0]:
      if len(tempocc) != 0:
        temp=np.array(npaOcc[:,nIdx-1])
        temp.shape=[temp.size,1]
        tempocc=np.append(tempocc, temp, axis=1)
      elif len(tempocc)==0:
        tempocc=np.array(npaOcc[:,nIdx-1])
        tempocc.shape=[tempocc.size,1]
    if tempocc!=[]:
      npaOcc=np.array(tempocc)
    else:
      print 'Warning getReducedOcc is returning empty!'  
    return npaOcc
      
#get errors on the levels used in the fit
  def getLevError(self, sErrorPath):
    levs=self.mllspec
    fErrors=open(sErrorPath, 'r')
    import numpy as np
    npaError=np.array([])
    for lev in levs:
      spec=[]
      spec=list(self.nAZ)
      spec.extend(lev.strip().split())
      for line in fErrors:
        line=line.strip().split()
        temp=[]
        temp.append(int(line[1]))
        temp.append(int(line[0]))
        temp.extend(line[2:-1])
        if np.array_equal(spec,temp):
          npaError=np.append(npaError,float(line[-1]))  
          break
    if len(self.mllspec)!=npaError.size:
      print "Warning: expected ",len(self.mllspec),' errors and instead found ',npaError.size, '.'
    return npaError
  
  #repeatedly calculate energies by sampling gaussians with standard deviation of 
  #matrix element errors and adding them to the orignal ME nSize times returns an
  #np array with shape [nSize,nEn]
  def accumulateEdist(self, npaErr, nSize, bAllME=False,bAllEn=True,fDefaultError=0.1):
    from numpy.random import randn
    from numpy import array
    from numpy import dot  
    from numpy import append
    if nSize>0:
      npaOrigME=array(self.getME(bAllME))
      npaEsample=array(self.getEnNu(bAllEn))
      if bAllME==True and npaErr.size<npaOrigME.size:
        npaErr=append(npaErr,[fDefaultError]*abs(npaErr.size-npaOrigME.size))
      npaEsample.shape=[1,npaEsample.size]
      for nIdx in range(nSize):
        npaGuess=array(npaOrigME+dot(randn(*(npaErr.shape)), npaErr))
        self.takeME(npaGuess)
        self.runSM()
        temp=array(self.getEnNu())
        temp.shape=[1,temp.size]
        if temp.size==npaEsample.shape[1]:
          npaEsample=append(npaEsample,temp,axis=0)
      self.takeME(npaOrigME)
      return npaEsample
    else:
      print 'nSize=', nSize, ' skip accumulation of data.'

#    normalized gaussian function
  def gaussian(self,x, mu, sig):
    import numpy as np    
    from math import pi
    from math import sqrt
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))/(sqrt(2.*pi)*sig)    
    
# calculate the error in energies by varying The ME within a given range
  def calcEThErr(self, npaErr, nSize=100, bPrev=False, bAllME=False,bAllEn=True):
    from numpy import histogram
    from math import sqrt
    from scipy.optimize import minimize
    from numpy import array
    from numpy import mean
    from numpy import std
    from numpy import dot
    from numpy import append
    if bPrev==False:
      npaSample=self.accumulateEdist(npaErr, nSize,bAllME,bAllEn)
      self.writeDat(npaSample)
    elif bPrev==True:
      npaSample=self.getPrevHist()
      if nSize>0:
        temp=self.accumulateEdist(npaErr, nSize)
        self.writeDat(temp)
        npaSample=append(npaSample,temp,axis=0)
    nSize=npaSample.shape[0]
    nBinNum=int(sqrt(float(nSize)))
    lHist=[]
    lBins=[]
    lMu=[]
    lSigma=[]
    for nEIdx in range(npaSample.shape[1]):
      npaHist,npaBinEdges=histogram(array(npaSample[:,nEIdx]),nBinNum, density=True)
      fMu0=mean(npaSample[:,nEIdx])
      fSigma0=std(npaSample[:,nEIdx])
#      save the histograms for later      
      lHist.append(npaHist)
      lBins.append(npaBinEdges)
#      fit a gaussian to the histogram
      x=[]
      for nIdx in range(npaBinEdges.size-1):
        x.append(0.5*(npaBinEdges[nIdx]+npaBinEdges[nIdx+1]))
      x=array(x)
      obj=lambda param: dot(self.gaussian(x,param[0],param[1])-array(npaHist), self.gaussian(x,param[0],param[1])-array(npaHist))
      x0=array([fMu0,fSigma0])
      from scipy.integrate import quad
      tempfun=lambda x: self.gaussian(x,x0[0],x0[1])
      print "pre",obj(x0), x0
      print quad(tempfun,-100,100)
      res=minimize(obj, x0)
      tempfun=lambda x: self.gaussian(x,x0[0],x0[1])      
      print 'post',obj(res.x), res.x
      print quad(tempfun,-100,100)
      lMu.append(res.x[0])
      lSigma.append(res.x[1])
    return lHist,lBins,lMu,lSigma, nSize
    
  def plotEthError(self, lHistm, lBins, lMu, lSigma,nSize):
    import Discrete
    import numpy as np
    
    best=Discrete.balFact(len(lHistm))
    from matplotlib import pyplot as plt
    fig, axray=plt.subplots(best[0], best[1])
    from itertools import chain
    try:
      axray= list(chain.from_iterable(axray))
    except TypeError:
      ''        
    nIdx=0
    for ax,hist,bins,mu, sigma in zip(axray, lHistm, lBins,lMu, lSigma):
      ax.bar(bins[:-1],hist, abs(bins[1]-bins[0]),label='Energy Hist N={:4d}'.format(nSize))
      xmin=max([np.amin(bins),mu-5.0*sigma])
      xmax=min([np.amax(bins),mu+5.0*sigma])
      x=np.linspace(xmin,xmax, 200)
      ax.plot(x,self.gaussian(x,mu,sigma), color='r',label='Fit: $\mu={:3.1f} \sigma={:2.1e}$'.format(mu,sigma))
      ax.legend(loc='best')
      ax.set_title('[A,Z]='+str(self.nAZ)+', ($J n_J \pi$)=('+str(self.mllspec[nIdx]).strip()+')')
      ax.set_xlabel('Energy(MeV)')
      ax.set_ylabel('$dP/dE$ (MeV$^{-1}$)')
      nIdx+=1
    plt.tight_layout(.025)
    plt.show()    
    
  def writeDat(self, npaSample):    
    fOut=open(self.sPath+'\\'+self.sName+'\\'+'tracking'+'\\energydist.dat','a+')
    sFormat='{:10.5f}\t'
    for nRIdx in range(npaSample.shape[0]):
      sNew=''
      for nCIdx in range(npaSample.shape[1]):
        sNew+=sFormat.format(npaSample[nRIdx,nCIdx])
      fOut.write(sNew+'\n')
    fOut.close()
    
  def getPrevHist(self):
    fIn=open(self.sPath+'\\'+self.sName+'\\'+'tracking'+'\\energydist.dat','r')
    import numpy as np
    npaDat=np.array([])
    for line in fIn:
      line=line.strip().split()
      line=np.array(line,dtype=float)
      line.shape=[1,line.size]
      if npaDat!=np.array([]):
        npaDat=np.append(npaDat,np.array(line),axis=0)
      elif np.all(npaDat==np.array([])):
        npaDat=np.array(line)
    return npaDat

# get rid of the interaction file and rerun the shell model calulation
  def reinitialize(self):
    import os
    os.remove(self.makeIntPath(''))
    self.runsm()

  def writeOcc(self):
      '''
          Write the single particle occupation numbers in the energy levels
          into a sing;e file.
      '''
      lsJPiList = self.makeJPiList()
      import numpy as np
      npaOutray = np.array([])
      sEnPath = self.sPath + '\\' + self.sName + '\\'
      lOldSPList=[]
      lSPList=[]
      nIter = 0
      from copy import copy
      for elem in lsJPiList:
          fIn = open(sEnPath + self.makeEnergyName(elem[0], self.sIsospin,
                                                   elem[1]) + '.lpe', 'r')
          llfOcc = []          
          nJ=0
          bIsIt = False
          lOldLine = []
          for line in fIn:
              line = line.strip().split()
              bComp1 = False
              if len(line)>=2:
                  if line[0] == 's-p' and line[1] == 'energy:':
                      lSPList = lOldLine
              if len(line) >=3:
                  bComp1 = (line[0] == 'no' and line[1] == 'energy' and
                            line[2] == 'level' and line[3] == 'average')
                  if not bIsIt and bComp1:
                      bIsIt = True
              if len(line) > 0:
                  bComp2 = bIsIt and not (len(line) == 0) and not bComp1
              if bComp2:
                  '''nj, occ1, occ2, occ3, ...'''
                  llfOcc.append(line)
              else:
                  lOldLine = line
                  continue
          fIn.close()
          if np.all(lSPList != lOldSPList) and nIter>0:
              print 'Warning: Single particle list has Changed from: '
              print lOldSPList
              print 'to:'
              print lSPList
              print 'In iteration: ', nIter, 'file: ' 
              print (self.makeEnergyName(elem[0], self.sIsospin, elem[1]) +
                     '.lpe')
          nNumSPS = len(lSPList) 
          for Occ in llfOcc:
              tOccDType = ('Occ', float, (1,nNumSPS))
              dtype = [('nJ', int), ('En', float), tOccDType, ('J', int),
                       ('Pi', 'S2')]
              values = (Occ[0], Occ[1],Occ[2:], int(2*eval(elem[0]+'.')), elem[1])
              if npaOutray.size == 0:
                  npaOutray = np.array(values, dtype=dtype)
              else:
                  temp = np.array(values, dtype=dtype)
                  npaOutray = np.append(npaOutray, temp)
          nIter += 1
          lnRmList = []
          for nIdx in range(npaOutray.size):
              for nJIdx in range(nIdx):
                  if np.all(npaOutray[nJIdx] == npaOutray[nIdx]):
                      lnRmList.append(nIdx)
      npaOutray = np.delete(npaOutray, lnRmList)
      npaOutray = np.sort(npaOutray, order=['En', 'J', 'Pi', 'nJ'])
      sOccName = self.getLevName()
      sOccName = sOccName[:-3] + 'occ'
      fOut = open(sEnPath + sOccName, 'w')
      sHFormat = '{:>5}{:>5}{:>8}{:>3}{:>3}'+'{:>8}'*nNumSPS
      sFormat = '{:>5d}{:>5d}{:>8.3f}{:>3}{:>3}'+'{:>8.2f}'*nNumSPS
      fOut.write(sHFormat.format('N', 'nJ', 'E_ex', '2J', 'p', *lSPList)+'\n')      
      for nIdx, elem in enumerate(npaOutray):
          fOut.write(sFormat.format(nIdx, elem['nJ'],
                                    elem['En'] - npaOutray[0]['En'], elem['J'],
                                    elem['Pi'],*elem['Occ'][0][:])+'\n')
      fOut.close()

  def makeJPiList(self):
      '''
          Make a  list of the spin and parities of the levels being tracked. 
      '''
      lsJPiList = []
      for elem in self.mllspec:
          bIsIt = False
          temp = [elem[0], elem[2]]
          for current in lsJPiList:
              if current[0] == temp[0] and current[1] == temp[1]:
                  bIsIt == True
                  break
          if not bIsIt:
              lsJPiList.append(temp)
      return lsJPiList
          
  def writeEnergies(self):
      '''
          Collect the energy from the files in the oxbash output and make a
          *.lpt file with all the energies.
      '''
      lsJPiList = self.makeJPiList()
      import numpy as np
      npaOutray = np.array([])
      sEnPath = self.sPath + '\\' + self.sName + '\\'
      for elem in lsJPiList:
          fIn = open(sEnPath + self.makeEnergyName(elem[0], self.sIsospin,
                                                   elem[1]) + '.lpe', 'r')
          llfEn = []          
          nJ=0
          bIsIt = False
          for line in fIn:
              line = line.strip().split()
              bComp1 = False
              if len(line) >=3:
                  bComp1 = (line[0] == 'eigenvalues' and line[1] == 'obtained' and
                            line[2] == 'in' and line[3] == 'last')
                  if not bIsIt and bComp1:
                      bIsIt = True
              if len(line) > 0:
                  bComp2 = bIsIt and not (len(line[0]) == 0) and not bComp1
              if bIsIt and len(line) == 0:
                  break
              elif bComp2:
                  for En in line:
                      if En > 0:
                          nJ += 1
                          llfEn.append([nJ, float(En)])
              else:
                  continue
          fIn.close()
          for En in llfEn:
              if npaOutray.size == 0:
                  dtype = [('nJ', int), ('En', float), ('J', 'S3'),
                           ('Pi', 'S2')]
                  values = (En[0], En[1], elem[0], elem[1])
                  npaOutray = np.array(values, dtype=dtype)
              else:
                  dtype = [('nJ', int), ('En', float), ('J', 'S3'),
                           ('Pi', 'S2')]
                  values = (En[0], En[1], elem[0], elem[1])
                  temp = np.array(values, dtype=dtype)
                  npaOutray = np.append(npaOutray, temp)
      lnRmList = []
      for nIdx in range(npaOutray.size):
          for nJIdx in range(nIdx):
              if np.all(npaOutray[nJIdx] == npaOutray[nIdx]):
                  lnRmList.append(nIdx)
      npaOutray = np.delete(npaOutray, lnRmList)
      npaOutray = np.sort(npaOutray, order=['En', 'J', 'Pi', 'nJ'])
      fOut = open(sEnPath + self.getLevName(), 'w')
      fOut.write('\n' + '-'*60 + '\n\n')
      sHFormat = '{:>5}{:>5}{:>11}{:>8}{:>5}{:>5}{:>3}'
      sFormat = '{:>5d}{:>5d}{:>11.3f}{:>8.3f}{:>5}{:>5}{:>3}'
      fOut.write(sHFormat.format('N', 'nJ', 'E(MeV)', 'E_ex', 'J', 'T_z', 'p')+'\n')      
      for nIdx, elem in enumerate(npaOutray):
          fOut.write(sFormat.format(nIdx, elem['nJ'], elem['En'],
                                    elem['En'] - npaOutray[0]['En'], elem['J'],
                                    self.sIsospin, elem['Pi'])+'\n')
      fOut.close()
  def makeEnergyName(self, sJ, sIsospin, sParity):
    '''
        Make the oxbash filename for the energy output file given the input
        from the state labels and additional information, returns a string of
        the file name without the extension.
    '''
    sName = ''
    sSpCode = self.lookupLab(self.sInt, 'MS')
    if type(sSpCode) == type(None):
        sSpCode = 'x'
    sName += sSpCode
#    if this gets used elswhere it may make sense to make it a global variable
    dCode = dict({'0': '0', '1/2': '1', '1': '2', '3/2': '3', '2': '4',
                  '5/2': '5', '3': '6', '7/2': '7', '4': '8', '9/2': '9',
                  '5': 'a', '11/2': 'b', '6': 'c', '13/2': 'd', '7': 'e',
                  '15/2': 'f', '8': 'g', '17/2': 'h', '9': 'i', '19/2': 'j',
                  '10': 'k', '21/2': 'l', '11': 'm', '23/2': 'n', '12': 'o',
                  '25/2': 'p', '12': 'q', '27/2': 'r', '14': 's', '29/2': 't',
                  '15': 'u', '31/2': 'v', '16': 'w', '33/2': 'x', '17': 'y'})
    sName += dCode[sJ]
    temp = sIsospin.split('/')
    if len(temp) > 0 and int(temp[0]) % 2 == 0:
        sIsospin = str(int(eval(sIsospin)))
    sName += dCode[sIsospin]
    dParity = dict({'+1': '0', '-1': '1'})
    sName += dParity[sParity]
    if self.nVal % 2 == 1:  
        sVal = str(self.nVal) + '/2'
    elif self.nVal % 2 == 0:
        sVal = str(self.nVal/2)
    sName += dCode[sVal]
    sIntCode = self.lookupLab(self.sInt, 'Int')
    if type(sIntCode) == type(None):
        sIntCode = 'y'        
    sName += sIntCode
    return sName

  def lookupLab(self, sLab, sType):
    '''
        Look up the label associated with the string in oxbash label.dat file.
        And return the result as a character string.
    '''
    fLabs = open(self.sOBDir + '\\sps\\label.dat','r')
    if sType == 'MS':
      nLabColIdx = 0
      nCharColIdx = 2
    elif sType == 'Int':
      nLabColIdx = 1
      nCharColIdx = 3
    else:
      print 'Error: invalid Label type. Argument 2 must be "MS" for model\
             space or "Int" for interaction.'
    for line in fLabs:
      line = line.strip().split() 
      if len(line) == 0 or line[0][0] == '!':
        continue
      elif sLab == line[nLabColIdx]:
          return line[nCharColIdx]
    print ('Error: Label "' + sLab +'" not found in ' + self.sOBDir +
           '\\sps\\label.dat')
    return None
    
#    return the two body trasition densities in an array where each row is a
#    different energy level and each column is a the transition label
    def getTBTD(self):
        npaTPTD = []
        npaLab = []
        for lev in self.mllspec:
            sLevNam = self.makeEnergyName(lev[0], self.sIsospin, lev[2])
            sDenFileName = sLevName + sLevName[1:3] +'.lbd'
            fDen = open(self.sPath + '\\' + sDenFileName, 'r')
            tempLab =[]
            tempDen = []
            for line in fDen:
                line = line.strip().split()
                try:
                    tempLab.append([int(line[nIdx]) for nIdx in range(4)])
                    tempDen.append()