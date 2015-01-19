# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 14:09:19 2015

Self contained utility classes for Shell Opt

@author: jonesad
"""

class MEhandler:
  'Takes lists of ME Types and the and list of lists of ME to be used in methods that read and write them.'
  def __init__(self, anBody, sPath, llSpec=[[],[]]):
    #assign matrix element type specification    
    self.manBody=anBody
    #assign the list of ME specifications (tells shich ME to change
    self.mllSpec=llSpec
    #assign path associated with list for reading and writing 
    self.msPath=sPath
    #Determine what number of each type of ME are to be read and written
    self.nMEnum=[]
    for nIdx, elem in enumerate(lType):
      if len(llSpec[nIdx])!=0:
        self.nMEnum.append(len(llSpec[nIdx]))
      elif len(llSpec[nIdx])==0 and elem==1:
        self.nMEnum.append(self.countOBME())
      elif len(llSpec[nIdx])==0 and elem==2:
        self.nMEnum.append(self.countTBME())  
  
  def writeOBME(self, npaME):
    import shutil
    shutil.copyfile(self.msPath, self.msPath+"_")
    fIntSrc=open(self.msPath+"_",'r') 
    fIntOut=open(self.msPath,'w')
    for line in fIntSrc:
      if line[0]!='!':
        line=line.strip().split()
        sStart="{0:3}{1:>11.4f}"
        sNext="{0:>10.4f}"
        if len(self.llspec[0])==0:
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
    
  def writeTBME(self,npaME):
    import shutil
    shutil.copyfile(self.msPath, self.msPath+"_")
    fIntSrc=open(self.msPath+"_",'r') 
    fIntOut=open(self.msPath,'w')
    sStart="{0:>3}"
    sNext="{0:>5}{1:>3}{2:>9.4f}"
    sNew=""          
    nElem=0
    nUnCm=0
    for line in fIntSrc:
      if line[0]!='!':             
        if  len(self.llspec[1])==0 and nUnCm!=0:
          line=line.strip().split()
          for i in range(4):
            sNew=sNew+sStart.format(line[i])
          sNew=sNew+sNext.format(line[4],line[5], float(npaME[nElem]))
          fIntOut.write(sNew+"\n")          
          nElem=nElem+1
        elif len(self.llspec[1])!=0 and nUnCm!=0 and self.llSpec[1][nElem]==nUnCm:
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
  #return total number of OBME in the interaction file 
  def countOBME(self):
    fIntSrc=open(self.msPath,'r')
    for line in fIntSrc:
      if line[0]!='!':
        return len(line.strip().split())-3
  #rerturn the total number of TBME in the interaction file
  def countTBME(self):
    fIntSrc=open(self.msPath)
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