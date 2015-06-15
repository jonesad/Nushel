# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 16:47:49 2015

class for managing a single nucleus

@author: jonesad
"""
import sys
sys.path.append('C:\\PythonScripts\\NushellScripts')
import ShellOptFl

class nucleus(ShellOptFl.MEhandler):
  'class for managing a single nucleus'
#take the path and create a sub directory for the nucleus and run a default calculation for the nucleus.
  def __init__(self,  nZ, nA, fGSE,sPath, sMMDR, sPar, lsShared,llMESpec,useGS,sForm,initialize=True, bExtrap=True):
    import os
    self.sPath=sPath
    self.nAZ=[nA,nZ]    
    self.sName='A'+str(nA)+'_Z'+str(nZ)
    if not os.path.exists(sPath+'\\'+self.sName):
      os.makedirs(self.sPath+'\\'+self.sName)
    if not os.path.exists(sPath+'\\'+self.sName+'\\'+'tracking'):
      os.makedirs(self.sPath+'\\'+self.sName+'\\'+'tracking')  
    self.writeAns(sMMDR, sPar,lsShared)
    self.sInt=lsShared[2]
    self.llMESpec=llMESpec
    if initialize:
      self.runSM()
    self.fGSE=fGSE
    self.useGS=useGS
    self.bExtrap=bExtrap
    #initialize the single particle matrix elements
    if len(llMESpec[0])==0:
      self.llMESpec[0]=range(1,self.countOBME()+1)
    self.sForm=sForm
    
    
#make a '.ans' file for use with Nushellx   
  def writeAns(self, sMMDR, sPar,lsShared):
    fAns=open(self.sPath+'\\'+self.sName+'\\'+self.sName+'.ans','w')
    fAns.write("--------------------------------------------------\n")
    sForm='{:21s}! '
    sFormN='{:3d}{:18s}! '
    fAns.write(sForm.format('lpe,   0')+'option (lpe or lan), neig (zero=10) \n')
    fAns.write(sForm.format(lsShared[0])+'model space (*.sp) name (a8)\n')
    fAns.write(sForm.format(lsShared[1])+'any restrictions (y/n)\n')
    fAns.write(sForm.format(lsShared[2])+'interaction (*.int) name (a8)\n')
    fAns.write(sFormN.format(self.nAZ[1],' ')+'number of protons\n')
    fAns.write(sFormN.format(self.nAZ[0],' ')+'number of nucleons\n')
    fAns.write(sForm.format(sMMDR)+'min J, max J, del J\n')
    fAns.write(sForm.format(sPar)+'parity (0 for +) (1 for -) (2 for both)\n')
    fAns.write("--------------------------------------------------\n")
    fAns.write("st                   ! option \n")
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
  
  #Store attribute LLspec for later use (separate from initialization)
  def setLevels(self, llSpec):
    self.mllspec=llSpec
    self.writeStatus()
    
  def getLevName(self):
    fIn=open(self.sPath+'\\'+self.sName+"\\list.lpt")
    sLevName=fIn.readline().strip()
    fIn.close()
    return sLevName
#  Get the energy diff for the experimental and calculated levels 
  def Ediff(self, bTrackDiff=False):
    afETh=self.getEnNu()
#    print 'Eth',afETh
    afEExp=self.getEExp()
#    print  'Eexp',afEExp
    import numpy as np
    #check if both the energies were found
    if len(afEExp)==len(afETh):
      res=np.absolute(np.array(afEExp)-np.array(afETh))
    else:
      #default value 
      res=np.array([float(abs(len(afEExp)-len(afETh))*100)])
      print afEExp
      print afETh
      print "Error: The number of Theoretical and experimental energies are not the same. Using default value of 100 MeV per missing level."
    if bTrackDiff:
      fOut=open(self.sPath+'\\'+self.sName+'\\'+'tracking'+'\\EDiff.dat','a+')
      fOut.write('Eth'+str(afETh)+'\n')
      fOut.write('EExp'+str(afEExp)+'\n')
      fOut.write('Ediff'+str(res)+'\n')      
      fOut.close()
    return res
    
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
          lev=lev.strip().split()
  #        print '\n'
  #        print self.sName, lev 
  #        print line
  #        print '\n'
          if len(line)>=6 and lev[0]==line[4] and lev[1]==line[1] and lev[2]==line[6]:
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
#    print afETh
#    raw_input("Press Enter to continue...")
    
    fTh.close()
    return afETh 

#get experimental energy     
  def getEExp(self):     
    import numpy as np
    sLevName=self.getLevName()
    afEExp=[]
#    print self.sPath+'\\'+self.sName+"\\"+sLevName[0:2]+'0'+sLevName[2:4]+'exp.lpt'
    fExp=open(self.sPath+'\\'+self.sName+"\\"+sLevName[0:2]+'0'+sLevName[2:4]+'exp.lpt')
    npaJ=np.zeros(len(self.mllspec))

    for line in fExp:
      line=line.strip().split()
      for nlevIdx,lev in enumerate(self.mllspec):
        lev=lev.strip().split()
#        print '\n\n'
#        print len(lev[1])
#        print lev
#        print '\n\n'
        for string in line:
#          print '\n\n'
#          print len(string)
#          print len(lev[0])+1
#          print '\n\n'
#          print npaJ
          if len(string)==len(lev[0])+1 and  string[:-1]==lev[0] and string[-1]==lev[2][0]:  
            npaJ[nlevIdx]+=1
#            print int(npaJ[nlevIdx]), int(lev[1])
            if  int(npaJ[nlevIdx])==int(lev[1]):
#              print 'chck 2 passed'
              afEExp.append(float(line[0]))
              break
    fExp.close()
#    print afEExp
    afEExp=np.array(afEExp)
#    print afEExp
    if self.useGS==1:
      afEExp=afEExp+self.fGSE
#    print afEExp
#      print self.fGSE
    if len(afEExp)!=len(self.mllspec):
      print "Error: # of Experimental levels found does not match requested # in nAZ=", self.nAZ
      print "Requested:", len(self.mllspec)
      print "Found:", len(afEExp)
    return afEExp
    
      
#run the shell model calculation        
  def runSM(self):
    import os
    #print os.getcwd()
    os.chdir(self.sPath+'\\'+self.sName)
    #print os.getcwd()    
    os.system('shell '+self.sName+'.ans')
    os.system(self.sName)


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
    for nLevIdx in range(npaMono.shape[0]):    
      nMono=0
      for nIdx in range(npaLabel.shape[0]):        
#        print npaLabel[nLevIdx,1], npaLabel[nLevIdx,3]
#        print npaLabel[nLevIdx,1]==npaLabel[nLevIdx,3]
#        print npaLabel[nLevIdx,0], npaLabel[nLevIdx,2]
#        print npaLabel[nLevIdx,0]==npaLabel[nLevIdx,2]
#        print npaLabel[nLevIdx,4]
#        print npaLabel[nLevIdx,:]
#        if npaLabel[nIdx,1]==npaLabel[nIdx,3] and npaLabel[nIdx,0]==npaLabel[nIdx,2] and npaLabel[nIdx,4]!=0:
        if np.all(npaLabel[nIdx]==npaMonoLabel[nMono]):
          temp=float(2*(npaLabel[nIdx,4]+1))
          if self.sForm =='pn':
            nSPEIdx=int(npaLabel[nIdx,0])-1
          elif self.sForm =='iso':
            nSPEIdx=int(npaLabel[nIdx,0])-1+npaLabel[nIdx,-1]*self.countOBME()
          if npaLabel[nIdx,0]!=npaLabel[nIdx,1]:
            npaMono[nLevIdx,nMono]+=npaOcc[nLevIdx,nSPEIdx]*npaOcc[nLevIdx,nSPEIdx]*temp
#            print npaLabel[nIdx,0],npaOcc[nLevIdx,int(npaLabel[nIdx,0]-1)]        
#            print npaLabel[nIdx,2],npaOcc[nLevIdx,int(npaLabel[nIdx,2]-1)]
#            print '\n'
          else:            
            npaMono[nLevIdx,nMono]+=npaOcc[nLevIdx,nSPEIdx]*(npaOcc[nLevIdx,nSPEIdx]-1.0)*temp/2.0
#            print npaLabel[nIdx,0],npaOcc[nLevIdx,int(npaLabel[nIdx,0]-1)]        
#            print npaLabel[nIdx,2],npaOcc[nLevIdx,int(npaLabel[nIdx,2]-1)]
#            print '\n'
          denom[nLevIdx,nMono]+=temp
#          print 'j*(j+1)=', temp
          nMono+=1
          if nMono>=npaMonoLabel.shape[0]:
            break
#    print denom
#    raw_input("Press enter to continue...")
    npaMono=np.divide(npaMono,denom)
#    print npaMono.max()
#    print npaMono
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
    
#    print npaOcc.shape[0], nMonoSize
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
#    print denom
#    print npaMono
    npaMono=np.divide(npaMono,denom)
#    print npaMono

#    print npaMono.max()
#    print npaMono
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
#          print np.all(npaNewLabels[nIdx2]==npaMonoLabel[nIdx3,:4])
          if np.all(npaNewLabels[nIdx2]==npaMonoLabel[nIdx3,:4]):
            npaNewMO[nIdx1,nIdx2]+=npaMono[nIdx1][nIdx3]
    temp=np.append(npaOcc,npaNewMO,axis=1)
#    print npaMono
#    print npaNewMO
    return temp, npaNewLabels

  
#get the occupation numbers   
  def getOcc(self, sLevName):
      import numpy as np
      fIn=open(self.sPath+'\\'+self.sName+'\\'+sLevName[:-4]+'.occ', 'r')
      npaOcc=[]
      nIdx=0
      for line in fIn:
        line=line.strip().split()
        if nIdx<2:
          nIdx+=1
          continue
        for nlevIdx,lev in enumerate(self.mllspec):
          lev=lev.strip().split()
          if int(line[3])==int(2*float(eval(lev[0]+'.0'))) and int(line[1])==int(lev[1]):
            if self.sForm=="pn":
              temp=line[5:5+self.countOBME()]
            elif self.sForm=="iso":
              temp=np.array(line[5:5+2*self.countOBME()],dtype=float)
            else:
              print 'Error: invalid Formalism specification:', self.sForm
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
#        print spec, temp
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
#    print bPrev
    if bPrev==False:
      npaSample=self.accumulateEdist(npaErr, nSize,bAllME,bAllEn)
      self.writeDat(npaSample)
    elif bPrev==True:
      npaSample=self.getPrevHist()
#      print 'init'
      if nSize>0:
        temp=self.accumulateEdist(npaErr, nSize)
        self.writeDat(temp)
        npaSample=append(npaSample,temp,axis=0)
#    print npaSample
    nSize=npaSample.shape[0]
    nBinNum=int(sqrt(float(nSize)))
#    print npaSample.shape,nBinNum
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
#      print npaHist
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