# -*- coding: utf-8 -*-
"""
Created on Thu Jan 15 14:09:19 2015

Self contained utility classes for Shell Opt

@author: jonesad
"""

class MEhandler:
  'Takes lists of ME Types and the and list of lists of ME to be used in methods that read and write them. The nucleus class inherits from here.'
  def __init__(self, anBody, sPath, sName, sInt, sForm, llMESpec=[[],[]], bExtrap=True):
    #assign matrix element type specification    
    self.setmanBody(anBody)    
    #assign the list of ME specifications (tells shich ME to change
    if len(llMESpec[0])==0:
      llMESpec[0]=list(range(7))
    self.llMESpec=llMESpec
    #assign path associated with list for reading and writing 
    self.sPath=str(sPath)
    #the Nucleus directory name
    self.sName=str(sName)
    #interaction name string without the .int extension
    self.sInt=str(sInt)
    #boolean variable that tracks whether the SPE line has entries relating to mass extrapolation of matrix elements
    self.bExtrap=bExtrap
    #initialize the single particle matrix elements
    self.llMESpec[0]=range(1,self.countOBME()+1)    
    self.sForm=sForm
    
#set manBody outside of initialization
  def setmanBody(self, anBody):
    import numpy as np
    try:
      test=[3]*len(anBody)
    except TypeError:
      test=3  
#    print anBody, test, anBody==test
    try:
      if any(anBody == test):
        self.manBody=[1,2]
        self.llMESpec[1]=np.array(self.getMonoLabel())
      else:
        self.manBody=anBody
    except TypeError:
      if anBody == test:
        self.manBody=[1,2]
        self.llMESpec[1]=np.array(self.getMonoLabel())  
      else:
        self.manBody=anBody
    self.setMEnum()
#    print self.manBody,self.nMEnum

#Determine what number of each type of ME are to be read and written
  def setMEnum(self):      
    self.nMEnum=[]
    for nIdx, elem in enumerate(self.manBody):
      if len(self.llMESpec[nIdx])!=0:
        self.nMEnum.append(len(self.llMESpec[nIdx]))
      elif len(self.llMESpec[nIdx])==0 and elem==1:
#        self.nMEnum.append(self.countOBME())
        self.nMEnum.append(0)
      elif len(self.llMESpec[nIdx])==0 and elem==2:
#        self.nMEnum.append(self.countTBME())
        self.nMEnum.append(0)
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
        if len(self.llMESpec[0])!=0 and nLine==0:
#          print "Writing only listed ME"
          for nElemIdx, elem in enumerate(self.llMESpec[0]):
            templine[elem]=npaME[nElemIdx]
#            print "templine["+str(elem)+']='+'npaME['+str(nElemIdx)+'-1]=',npaME[nElemIdx-1]
          for nIdx, elem in enumerate(templine):
            if nIdx==0:
              sNew+=sStart.format(templine[0], float(templine[1]))
            elif nIdx>1:
#              print "I am elem", elem
              sNew+=sNext.format(float(elem))
#          print sNew
          fIntOut.write(sNew+'\n')
        else:
          fIntOut.write(line)
        nLine+=1
      else:
        fIntOut.write(line)
        
    fIntSrc.close()
    fIntOut.close()
#    get the one body matrix elements
  def getOBME(self,bAll=False):
    import numpy as np
    if bAll==True:
      nOBME=self.countOBME()
    else:
      nOBME=len(self.llMESpec[0])              
#    print nOBME
    fIntSrc=open(self.makeIntPath(''),'r') 
    npaME=[]    
    for line in fIntSrc:
      if line[0]!='!':
        line=line.strip().split()
        if len(self.llMESpec[0])==0 or bAll==True:
          npaME.append(line[1:1+nOBME])
#          print "bAll is true"          
        elif len(self.llMESpec[0])!=0:
          for elem in self.llMESpec[0]:
#            print line 
#            print int(elem)
            npaME.append(line[int(elem)])
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
    nUnCm=0
    if self.llMESpec[1].shape[0]!=npaME.size:
      print 'Number of ME Labels does not match number of ME to be written: ' + str(self.llMESpec[1].shape[0])+ ' labels and ' + str(npaME.size) + ' ME.'
      print self.llMESpec[1]
      print npaME
      raw_input('press enter')
    import numpy as np
    for line in fIntSrc:
      flag=True
      sNew=""          
      if line[0]!='!':             
        temp1=line.strip().split()
        if nUnCm>0:
          if len(self.llMESpec[1])!=0:
            temp2=[]
            for nIdx in range(len(self.llMESpec[1][0,:])):
              if nIdx!=0:
                temp2.append(int(temp1[nIdx]))
              else:
                temp2=[int(temp1[nIdx])]
          for nElem, npaMELab in enumerate(self.llMESpec[1]):
            if len(self.llMESpec[1])!=0 and nUnCm!=0 and np.all(npaMELab == np.array(temp2)):
              sNew=sNew+sStart.format(temp1[0])          
              for i in range(1,4):
                sNew=sNew+sNext.format(temp1[i])
              sNew=sNew+sEnd.format(temp1[4],temp1[5], float(npaME[nElem]))
              fIntOut.write(sNew+"\n")          
              flag=False
          if flag:
            fIntOut.write(line)        
        else:
          fIntOut.write(line)
        nUnCm+=1

      else:
        fIntOut.write(line)
    fIntSrc.close()
    fIntOut.close()
#Get TBME
  def getTBME(self, bAll=False, nExtrap=0):
    fIntSrc=open(self.makeIntPath('',nExtrap),'r') 
    import numpy as np
    size= self.llMESpec[1].shape[0]
    npaME=np.zeros([size,1])
    nElem=0
    nUnCm=0
    for line in fIntSrc:
      if line[0][0]!='!':
        line=line.strip().split()
        if nUnCm>0:
          temp=[]
          for nIdx in range(6):
            if nIdx!=0:
              temp.append(int(line[nIdx]))
            else:
              temp=[int(line[nIdx])]
          temp=np.array(temp)
          temp.shape=[1,temp.size]
        if  bAll==True and nUnCm!=0:
          if npaME!=[]:
            npaME.append(float(line[6]))
          else:            
            npaME=[line[6]]
          nElem=nElem+1
        elif len(self.llMESpec[1])!=0 and nUnCm!=0:
          for nLabIdx,lab in enumerate(self.llMESpec[1]):
            if np.all(lab==temp):
              npaME[nLabIdx]=float(line[6])
              nElem=nElem+1
              break
        nUnCm+=1
    if (not bAll and self.llMESpec[1].shape[0]!=len(npaME)):
      raw_input('Warning expected: '+str(self.llMESpec[1].shape[0])+' ME. But found: '+ str(len(npaME)))        
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
#get the labels for the monopole matrix elements
  def getMonoLabel(self, npaBase=[]):
    import numpy as np
    npaLabel=self.getLabel()
    npaMonoLabel=np.array([])
    if len(npaBase)==0:
      nMono=0
      for nIdx in range(npaLabel.shape[0]):
        nIsit=0
        if npaLabel[nIdx,1]==npaLabel[nIdx,3] and npaLabel[nIdx,0]==npaLabel[nIdx,2]:
          nMono+=1
          for nJj in range(npaMonoLabel.shape[0]):
            if npaMonoLabel!=np.array([]):
              if npaMonoLabel.size==6 and np.all(npaMonoLabel==npaLabel[nIdx,:6]):
                nIsit=1
                break
              elif npaMonoLabel.size>6: 
                if np.all(npaLabel[nIdx,:6]==npaMonoLabel[nJj,:]):
                  nIsit=1
                  break
        else:
          nIsit=1        
        if nIsit==0:
          if npaMonoLabel!=np.array([]):
            npaMonoLabel=np.append(npaMonoLabel,np.array([npaLabel[nIdx,:6]]),0)
          else:
            npaMonoLabel=np.array([npaLabel[nIdx,:6]])
    else:
      for nIdx in range(npaLabel.shape[0]):
        for nJIdx in range(npaBase.shape[0]):
          if np.all(npaLabel[nIdx,:4]==npaBase[nJIdx,:]):
            if npaMonoLabel!=np.array([]):
              npaMonoLabel=np.append(npaMonoLabel,np.array([npaLabel[nIdx,:6]]),0)
            else:
              npaMonoLabel=np.array([npaLabel[nIdx,:6]])
                
    return npaMonoLabel

#get the monopole interaction matrix elements
  def getMonoME(self):
#    fIntSrc=open(self.makeIntPath(''),'r') 
#    nUnCm=0
    npaMME=[]
#    
#    for line in fIntSrc:
#      if line[0]!='!':
#        temp=line.strip().split()        
#        if nUnCm!=0 and temp[0]==temp[2] and temp[1]==temp[3] and int(temp[4])!=0:          
#          npaMME.append(temp[6])
#        nUnCm+=1
#    fIntSrc.close()
    import numpy
    self.llMESpec[1]=self.getMonoLabel()
    self.setMEnum()
    npaMME=self.getTBME()
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
      if line[0]!='!' and self.bExtrap:
#        print self.bExtrap, "so sub 4"
        return len(line.strip().split())-4
      elif line[0]!='!' and (not self.bExtrap):
#        print self.bExtrap, " so sub 1"
        return len(line.strip().split())-1

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
    print self.nMEnum
    self.setMEnum()
    print self.nMEnum
    
    for nCt, elem in enumerate(self.nMEnum):
      if self.manBody[nCt]==1:
        npaOBME=npaME[:self.nMEnum[nCt]]
#        print 'self.nMEnum[nCt]', self.nMEnum[nCt]
        print ' npaOBME', npaOBME
      elif self.manBody[nCt]==2:
        npaTBME=npaME[sum(self.nMEnum[:nCt]):sum(self.nMEnum[:nCt+1])]
#        print 'sum(self.nMEnum[:nCt])', sum(self.nMEnum[:nCt])
#        print  'sum(self.nMEnum[:nCt+1])', sum(self.nMEnum[:nCt+1])
        print ' npaTBME', npaTBME
    if len(npaOBME)>0:
      self.writeOBME(npaOBME) 
    if len(npaTBME)>0:
      self.writeTBME(npaTBME) 
      
#return the specified matrix elements of the current hamiltonian.
  def getME(self, bAll=False):
    import numpy as np    
    npaME=np.array([])
    temp=[]
    if bAll==False:
      for elem in self.manBody:
        if elem==1:
          temp=self.getOBME()
        elif elem==2 :
          temp=self.getTBME()
        npaME=np.append(npaME,temp)
    elif bAll==True:
      npaME=self.getOBME(bAll=bAll)
      npaME=np.append(npaME,self.getTBME(bAll=bAll))

    return npaME
    
#Check the interaction file for repititions of line
  def checkRep(self):
    fIntSrc=open(self.makeIntPath(''), 'r')
    myLines=[]    
    for line in fIntSrc:
      line=line.strip()
      myLines.append(line)
    if len(myLines)==len(set(myLines)):
      print "Each line in the file is unique"
    else:
      print "Repititions exist"
      print "Of", len(myLines), 'lines', len(set(myLines)), ' are unique.'

#Check the interaction file for repititions of labels
  def checkLabelRep(self):
    fIntSrc=open(self.makeIntPath(''), 'r')
    myLines=[] 
    nOffset=-2
    for line in fIntSrc:
      line=line.strip().split()
      print line 
      print line[:nOffset]
      tempstr=''
      for elem in line[:nOffset]:
        tempstr+=elem
      myLines.append(tempstr)
      
    if len(myLines)==len(set(myLines)):
      print "Each label in the file is unique"
    else:
      print "Repititions exist"
      print "Of", len(myLines), 'labels', len(set(myLines)), ' are unique.'
      
    
#   add a certain value to the monopole term for a given label
  def addDiff(self, delta, label):  
    fIntSrc=open(self.makeIntPath(''),'r') 
    import numpy as np
    npaME=[]
    npaMESpec=[]    
    nUnCm=0
    for line in fIntSrc:
      if line[0][0]!='!':
        line=line.strip().split()
        if nUnCm>0:
          temp=[]
          for nIdx in range(6):
            if nIdx!=0:
              temp.append(int(line[nIdx]))
            else:
              temp=[int(line[nIdx])]
          temp=np.array(temp)
          temp.shape=[1,temp.size]
        for nIIIdx, lab in enumerate(label):
          if nUnCm!=0 and np.all(lab==temp):
            npaME.append(float(line[6])+delta[nIIIdx])            
            for nIIdx in range(6):
              if nIIdx!=0:
                temp.append(int(line[nIIdx]))
              else:
                temp=[int(line[nIIdx])]
            npaMESpec.append(temp)
        nUnCm+=1
    fIntSrc.close()
    
    return np.array(npaME,dtype=float), np.array(npaMESpec,dtype=int)
    