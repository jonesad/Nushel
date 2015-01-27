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
  def __init__(self,  nZ, nA, fGSE,sPath, sMMDR, sPar, lsShared,llMESpec):
    import os
    self.sPath=sPath
    self.nAZ=[nA,nZ]    
    self.sName='A'+str(nA)+'_Z'+str(nZ)
    if not os.path.exists(sPath+'\\'+self.sName):
      os.makedirs(sPath+'\\'+self.sName)
    self.writeAns(sMMDR, sPar,lsShared)
    self.sInt=lsShared[2]
    self.llMESpec=llMESpec
    self.runSM()
    self.fGSE=fGSE
    
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
  
  #Store attribute LLspec for later use (separate from initialization)
  def setLevels(self, llSpec):
    self.mllspec=llSpec
  
#  Get the energy diff for the experimental and calculated levels 
  def Ediff(self):
    fIn=open(self.sPath+'\\'+self.sName+"\\list.lpt")
    sLevName=fIn.readline().strip()        
    fIn.close()
    afETh=self.getEnNu(sLevName)
    afEExp=self.getEExp(sLevName)
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
    return res
    
# Get Nushell energies
  def getEnNu(self,sLevName):
    afETh=[]
    fTh=open(self.sPath+'\\'+self.sName+"\\"+sLevName)
    for line in fTh:
      line=line.strip().split()
      for lev in self.mllspec:
        lev=lev.strip().split()
#        print '\n\n'
#        print lev
#        print line
#        print '\n\n'
        if len(line)>=6 and lev[0]==line[4] and lev[1]==line[1] and lev[2]==line[6]:
          afETh.append(float(line[3]))
          break
    fTh.close()
    return afETh 

#get experimental energy     
  def getEExp(self,sLevName):     
    import numpy as np
    afEExp=[]
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
#          print string
#          print lev
#          print '\n\n'
#          print npaJ
          if len (string)==2 and string[0]==lev[0][0] and string[1]==lev[2][0]:  
            npaJ[nlevIdx]+=1
            print int(npaJ[nlevIdx]), int(lev[1])
            if  int(npaJ[nlevIdx])==int(lev[1]):
              print 'chck 2 passed'
              afEExp.append(float(line[0]))
              break
    fExp.close()
    return afEExp
      
#run the shell model calculation        
  def runSM(self):
    import os
    #print os.getcwd()
    os.chdir(self.sPath+'\\'+self.sName)
    #print os.getcwd()    
    os.system('shell '+self.sName+'.ans')
    os.system(self.sName)    

#get the occupation numbers   
  def getOcc(self, sLevName):
      import numpy as np
      fIn=open(self.sPath+'\\'+self.sName+'\\'+sLevName[:-4]+'.occ', 'r')
      npaOcc=[]
      nIdx=0
      for line in fIn:
        line=line.strip().split()
        if nIdx<2:
          continue
        else:
          nIdx+=1
        for nlevIdx,lev in enumerate(self.mllspec):
          lev=lev.strip().split()
          print int(line[0]), int(lev[0]), int(line[2]), int(lev[1])
          if int(line[0])==int(lev[0]) and int(line[2])==int(lev[1]):
            if npaOcc!=[]:
              temp=[line[8:11]]
              npaOcc=temp
            else:
              temp=[line[8:11]]
              npaOcc=np.append(npaOcc, temp, axis=0)
      fIn.close()
      return np.array(npaOcc,dtype=float)