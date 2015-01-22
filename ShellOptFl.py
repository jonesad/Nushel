# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 14:09:19 2015

Self contained utility classes for Shell Opt

@author: jonesad
"""

class MEhandler:
  'Takes lists of ME Types and the and list of lists of ME to be used in methods that read and write them. The nucleus class inherits from here.'
  def __init__(self, anBody, sPath, sName, sInt, llspec=[[],[]]):
    #assign matrix element type specification    
    self.setmanBody(anBody)    
    #assign the list of ME specifications (tells shich ME to change
    self.mllspec=llspec
    #assign path associated with list for reading and writing 
    self.sPath=str(sPath)
    #the Nucleus directory name
    self.sName=str(sName)
    #interaction name string without the .int extension
    self.sInt=str(sInt)
        
    
#set manBody outside of initialization
  def setmanBody(self, anBody):
    self.manBody=anBody
    #Determine what number of each type of ME are to be read and written
    self.nMEnum=[]
    for nIdx, elem in enumerate(anBody):
      if len(self.mllspec[nIdx])!=0:
        self.nMEnum.append(len(self.mllspec[nIdx]))
      elif len(self.mllspec[nIdx])==0 and elem==1:
        self.nMEnum.append(self.countOBME())
      elif len(self.mllspec[nIdx])==0 and elem==2:
        self.nMEnum.append(self.countTBME())
            
# make the path of the interaction file   
  def makeIntPath(self, sPreExt=''):
    temp=str(self.sPath)+'\\'+str(self.sName)+'\\'+str(self.sInt)+str(sPreExt)+'.int'
    return temp
    
  def writeOBME(self, npaME):
    import shutil
    shutil.copyfile(self.makeIntPath(), self.makeIntPath('_'))
    fIntSrc=open(self.makeIntPath('_'),'r') 
    fIntOut=open(self.makeIntPath(''),'w')
    for line in fIntSrc:
      if line[0]!='!':
        line=line.strip().split()
        sStart="{0:3}{1:>11.4f}"
        sNext="{0:>10.4f}"
        if len(self.mllspec[0])==0:
          sNew=sStart.format(line[0],float(npaME[0]))
          for nElemIdx, elem in enumerate(npaME):
            if nElemIdx!=0:
              sNew=sNew+sNext.format(float(elem))
              for i in range(len(line)-len(npaME)-1):
                sNew=sNew+sNext.format(float(line[len(npaME)+i+1]))
              fIntOut.write(sNew+"\n")
            else:
              fIntOut.write(line)
        elif len(self.mllSpec[0])!=0:
          nLElIdx=0
          for nElemIdx, elem in enumerate(self.mllSpec[0]):
            if nElemIdx==0 and elem==1:
              sNew=sStart.format(line[0],float(npaME[0]))
              nLElIdx+=1
            elif nElemIdx==0 and elem!=1:
              sNew=sStart.format(line[0],line[1])
              nLElIdx+=2
            elif nElemIdx!=0 and elem==nElemIdx:
              sNew=sNew+sNext.format(float(npaME[nElemIdx]))
              nLElIdx+=1
            elif nElemIdx!=0 and elem!=nElemIdx:
              sNew=sNew+sNext.format(elem)
            else:
              print "logic error"
              break
        break            
    fIntSrc.close()
    fIntOut.close()
#    get the one body matrix elements
  def getOBME(self):
    import numpy as np
    if len(self.mllspec[0])==0:
      nOBME=self.countOBME()
    else:
      nOBME=len(self.mllspec[0])              
    fIntSrc=open(self.makeIntPath(''),'r') 
    npaME=[]    
    for line in fIntSrc:
      if line[0]!='!':
        line=line.strip().split()
        if len(self.mllspec[0])==0:
          npaME.append(line[1:1+nOBME])          
        elif len(self.mllspec[0])!=0:
          for elem in self.mllspec[0]:
            print elem
            npaME.append(line[int(elem)])
        break            
    fIntSrc.close()
    return np.array(npaME,dtype=float)
    
  def writeTBME(self,npaME):
    import shutil
    shutil.copyfile(self.makeIntPath(''),self.makeIntPath('_'))
    fIntSrc=open(self.makeIntPath('_'),'r') 
    fIntOut=open(self.makeIntPath(''),'w')
    sStart="{0:>3}"
    sNext="{0:>5}{1:>3}{2:>9.4f}"
    sNew=""          
    nElem=0
    nUnCm=0
    for line in fIntSrc:
      if line[0]!='!':             
        if  len(self.mllspec[1])==0 and nUnCm!=0:
          line=line.strip().split()
          for i in range(4):
            sNew=sNew+sStart.format(line[i])
          sNew=sNew+sNext.format(line[4],line[5], float(npaME[nElem]))
          fIntOut.write(sNew+"\n")          
          nElem=nElem+1
        elif len(self.mllspec[1])!=0 and nUnCm!=0 and self.mllspec[1][nElem]==nUnCm:
          line=line.strip().split()
          for i in range(4):
            sNew=sNew+sStart.format(line[i])
          sNew=sNew+sNext.format(line[4],line[5], float(npaME[nElem]))
          fIntOut.write(sNew+"\n")          
          nElem=nElem+1
        else:
            fIntOut.write(line)
      else:
            fIntOut.write(line)
      nUnCm+=1
      fIntSrc.close()
      fIntOut.close()

#Get TBME
  def getTBME(self):
    fIntSrc=open(self.makeIntPath(''),'r') 
    import numpy as np
    npaME=np.array([])
    nElem=0
    nUnCm=0
    for line in fIntSrc:
      if line[0]!='!':             
        if  len(self.mllspec[1])==0 and nUnCm!=0:
          line=line.strip().split()
          npaME.append(float(line[6]))          
          nElem=nElem+1
        elif len(self.mllspec[1])!=0 and nUnCm!=0 and self.mllspec[1][nElem]==nUnCm:
          line=line.strip().split()
          npaME.append(float(line[6]))         
          nElem=nElem+1
      nUnCm+=1
      fIntSrc.close()
      return npaME
      
  #return total number of OBME in the interaction file 
  def countOBME(self):
    fIntSrc=open(self.makeIntPath(''),'r')
    for line in fIntSrc:
      if line[0]!='!':
        return len(line.strip().split())-4
  #rerturn the total number of TBME in the interaction file
  def countTBME(self):
    fIntSrc=open(self.makeIntPath(''))
    nLIdx=0
    for line in fIntSrc:
      if line[0]!='!':
        nLIdx+=1
    return nLIdx-1    
          
  #Call with an npa array of matrix elements use the memebers of self to 
  #properly write the matrix elements to the supplied file path  
  def takeME(self, npaME):
    npaOBME=[]
    npaTBME=[]
    for nCt, elem in enumerate(self.nMEnum):
      if self.manBody[nCt]==1:
        npaOBME=npaME[sum(self.nMEnum[:nCt]):sum(self.nMEnum[:nCt])+1]
      elif self.manBody[nCt]==2:
        npaTBME=npaME[sum(self.nMEnum[:nCt]):sum(self.nMEnum[:nCt])+1]      
    if len(npaOBME)>0:
      self.writeOBME(npaOBME) 
    if len(npaTBME)>0:
      self.writeTBME(npaTBME) 
      
#return the specified matrix elements of the current hamiltonian.
  def getME(self):
    npaME=[]
    for elem in self.manBody:
      if elem==1:
        npaME=self.getOBME()
      elif elem==2:
        npaME=self.gettBME()      
    return npaME
    
      