# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 14:09:19 2015

Self contained utility classes for Shell Opt

@author: jonesad
"""

class MEhandler:
  'Takes lists of ME Types and the and list of lists of ME to be used in methods that read and write them. The nucleus class inherits from here.'
  def __init__(self, anBody, sPath, sName, sInt, llMESpec=[[],[]]):
    #assign matrix element type specification    
    self.setmanBody(anBody)    
    #assign the list of ME specifications (tells shich ME to change
    self.llMESpec=llMESpec
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
      if len(self.llMESpec[nIdx])!=0:
        self.nMEnum.append(len(self.llMESpec[nIdx]))
      elif len(self.llMESpec[nIdx])==0 and elem==1:
        self.nMEnum.append(self.countOBME())
      elif len(self.llMESpec[nIdx])==0 and elem==2:
        self.nMEnum.append(self.countTBME())
            
# make the path of the interaction file   
  def makeIntPath(self, sPreExt='', nExtrap=0):
    if nExtrap==0:
      temp=str(self.sPath)+'\\'+str(self.sName)+'\\'+str(self.sInt)+str(sPreExt)+'.int'
    elif nExtrap==1:
      sLevName=self.getLevName()
      temp=str(self.sPath)+'\\'+str(self.sName)+'\\'+sLevName[:-5]+'0'+str(sPreExt)+'.int'      
    return temp

#write the obme to file
  def writeOBME(self, npaME):
    import shutil
    shutil.copyfile(self.makeIntPath(), self.makeIntPath('_'))
    fIntSrc=open(self.makeIntPath('_'),'r') 
    fIntOut=open(self.makeIntPath(''),'w')
    nLine=0
    for line in fIntSrc:
      if line[0]!='!':
        sNew=''
        sStart="{0:3}{1:>11.4f}"
        sNext="{0:>10.4f}"
        templine=line.strip().split()        
        if len(self.llMESpec[0])==0 and nLine==0:
          print 'writing all ME'
          sNew=sStart.format(templine[0],float(npaME[0]))
          for nElemIdx, elem in enumerate(npaME):
            if nElemIdx!=0:
              sNew=sNew+sNext.format(float(elem))
          for i in range(len(templine)-len(npaME)-1):
              sNew=sNew+sNext.format(float(templine[len(npaME)+i+1]))
          fIntOut.write(sNew+"\n")
        elif len(self.llMESpec[0])!=0 and nLine==0:
          print "Writing only listed ME"
          for nElemIdx, elem in enumerate(self.llMESpec[0]):
            templine[elem]=npaME[nElemIdx]
          for nIdx, elem in enumerate(templine):
            if nIdx==0:
              sNew+=sStart.format(templine[0], float(templine[1]))
            elif nIdx>1:
              sNew+=sNext.format(float(elem))
          print sNew
          fIntOut.write(sNew+'\n')
        else:
          fIntOut.write(line)
        nLine+=1
      else:
        fIntOut.write(line)
        
    fIntSrc.close()
    fIntOut.close()
#    get the one body matrix elements
  def getOBME(self):
    import numpy as np
    if len(self.llMESpec[0])==0:
      nOBME=self.countOBME()
    else:
      nOBME=len(self.llMESpec[0])              
    fIntSrc=open(self.makeIntPath(''),'r') 
    npaME=[]    
    for line in fIntSrc:
      if line[0]!='!':
        line=line.strip().split()
        if len(self.llMESpec[0])==0:
          npaME.append(line[1:1+nOBME])          
        elif len(self.llMESpec[0])!=0:
          for elem in self.llMESpec[0]:
            npaME.append(line[int(elem)+1])
        break            
    fIntSrc.close()
    return np.array(npaME,dtype=float)
    
  def writeTBME(self,npaME):
    import shutil
    shutil.copyfile(self.makeIntPath(''),self.makeIntPath('_'))
    fIntSrc=open(self.makeIntPath('_'),'r') 
    fIntOut=open(self.makeIntPath(''),'w')
    sStart="{0:>4}"
    sNext="{0:>3}"
    sEnd="{0:>4}{1:>3}{2:>14.5f}"
    nElem=0
    nUnCm=0
    for line in fIntSrc:
      sNew=""          
      if line[0]!='!' and nElem<npaME.shape[0]:             
        temp1=line.strip().split()
        if nUnCm>0:
          if len(self.llMESpec[1])!=0:
            temp2=[]
            for nIdx in range(len(self.llMESpec[1][0,:])):
              if nIdx!=0:
                temp2.append(int(temp1[nIdx]))
              else:
                temp2=[int(temp1[nIdx])]
        if len(self.llMESpec[1])==0 and nUnCm!=0:
          sNew=sNew+sStart.format(temp1[0])
          for i in range(3):
            sNew=sNew+sNext.format(temp1[i])
          sNew=sNew+sEnd.format(temp1[4],temp1[5], float(npaME[nElem]))
          fIntOut.write(sNew+"\n")          
          nElem=nElem+1
        elif len(self.llMESpec[1])!=0 and nUnCm!=0 and all(x in self.llMESpec[1][nElem] for x in temp2):
          sNew=sNew+sStart.format(temp1[0])          
          for i in range(3):
            sNew=sNew+sNext.format(temp1[i])
          sNew=sNew+sEnd.format(temp1[4],temp1[5], float(npaME[nElem]))
          fIntOut.write(sNew+"\n")          
          nElem=nElem+1
        else:
          fIntOut.write(line)
        nUnCm+=1
      else:
        fIntOut.write(line)
    fIntSrc.close()
    fIntOut.close()
#Get TBME
  def getTBME(self, nExtrap=0):
    fIntSrc=open(self.makeIntPath('',nExtrap),'r') 
    import numpy as np
    npaME=[]
    nElem=0
    nUnCm=0
    for line in fIntSrc:
      if line[0][0]!='!':
        line=line.strip().split()
        if nUnCm>0:
          temp=[]
          for nIdx in range(5):
            if nIdx!=0:
              temp.append(int(line[nIdx]))
            else:
              temp=[int(line[nIdx])]
        if  len(self.llMESpec[1])==0 and nUnCm!=0:
          if npaME!=[]:
            npaME.append(float(line[6]))
          else:            
            npaME=[line[6]]
          nElem=nElem+1
        elif len(self.llMESpec[1])!=0 and nUnCm!=0 and all(x in self.llMESpec[1][nElem] for x in temp):
          npaME.append(float(line[6]))         
          nElem=nElem+1
        nUnCm+=1
        if nElem >= len(self.llMESpec[1]) and len(self.llMESpec[1])!=0:
          break
        
    fIntSrc.close()
    return np.array(npaME,dtype=float)

#set monopole matrix elements
  def writeMonoME(self, npaME):
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
        if  nUnCm!=0:
          temp=line.strip().split()
          if temp[0]==temp[2] and temp[1]==temp[3] and int(temp[4])!=0:          
            for i in range(4):
              sNew=sNew+sStart.format(temp[i])
            sNew=sNew+sNext.format(temp[4],temp[5], float(npaME[nElem]))
            fIntOut.write(sNew+"\n")          
            nElem=nElem+1
          else:
            fIntOut.write(line)
        else:
            fIntOut.write(line)
      else:
            fIntOut.write(line)
      nUnCm+=1
      fIntSrc.close()
      fIntOut.close()
#get the monopole interaction matrix elements
  def getMonoME(self):
    fIntSrc=open(self.makeIntPath(''),'r') 
    nUnCm=0
    npaMME=[]
    for line in fIntSrc:
      if line[0]!='!':
        temp=line.strip().split()        
        if nUnCm!=0 and temp[0]==temp[2] and temp[1]==temp[3] and int(temp[4])!=0:          
          npaMME.append(temp[6])
        nUnCm+=1
    fIntSrc.close()
    import numpy
    return numpy.array(npaMME,dtype=float)
     
#get label
  def getLabel(self): 
    fIntSrc=open(self.makeIntPath(''),'r') 
    import numpy as np
    npaLabel=np.array([])
    nElem=0
    nUnCm=0
    for line in fIntSrc:
      if line[0]!='!':
        if nUnCm>0: 
          line=line.strip().split()
          temp=[]
          for elem in line[0:6]:
            temp.append(int(elem))
          if npaLabel!=np.array([]):
            npaLabel=np.append(npaLabel,[temp],axis=0)
          else:
            npaLabel=[temp]
          nElem=nElem+1
        nUnCm+=1
    
    fIntSrc.close()
    return npaLabel
      
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
        npaOBME=npaME[sum(self.nMEnum[:nCt]):sum(self.nMEnum[:nCt+1])]
#        print npaOBME
      elif self.manBody[nCt]==2:
        npaTBME=npaME[sum(self.nMEnum[:nCt]):sum(self.nMEnum[:nCt+1])]      
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
        npaME=self.getTBME()      
    return npaME
    
#Check the interaction for repititions of ME
  def checkRep(self):
    fIntSrc=open(self.makeIntPath(''), 'r')
    myLines=[]    
    for line in fIntSrc:
      myLines.append(line[:21])
    if len(myLines)==len(set(myLines)):
      print "Each line in the file is unique"
    else:
      print "Repititions exist"
      print "Of", len(myLines), 'lines', len(set(myLines)), ' are unique.'