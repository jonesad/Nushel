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
  def __init__(self,  nZ, nA, fGSE,sPath, sMMDR, sPar, lsShared,llMESpec,useGS):
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
#    import time
#    if llMESpec==[[],[]]:
#      print 'Warning no matrix elements specified...'
#      time.sleep(2)
    self.llMESpec=llMESpec
    self.runSM()
    self.fGSE=fGSE
    self.useGS=useGS
    
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
  def Ediff(self):
    afETh=self.getEnNu()
    afEExp=self.getEExp()
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
  def getEnNu(self):
    sLevName=self.getLevName()
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
          if self.useGS==0:
            afETh.append(float(line[3]))
          elif self.useGS==1:
            afETh.append(float(line[2]))
          break
    if len(afETh)!=len(self.mllspec):
      print "Error: # of Theory levels found does not match requested # in nAZ=", self.nAZ
      print "Requested:", len(self.mllspec)
      print "Found:", len(afETh)
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
#      print afEExp
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
#get the labels for the monopole matrix elements
  def getMonoLabel(self):
    import numpy as np
    npaLabel=self.getLabel()
    npaMonoLabel=np.array([])
    nMono=0
    for nIdx in range(npaLabel.shape[0]):
      nIsit=0
      if npaLabel[nIdx,1]==npaLabel[nIdx,3] and npaLabel[nIdx,0]==npaLabel[nIdx,2] :
        nMono+=1
        for nJj in range(npaMonoLabel.shape[0]):
          if npaMonoLabel!=np.array([]):
            if npaMonoLabel.size==4 and np.all(npaMonoLabel==npaLabel[nIdx,:4]):
              nIsit=1
              break
            elif npaMonoLabel.size>4: 
              if np.all(npaLabel[nIdx,:4]==npaMonoLabel[nJj,:]):
                nIsit=1
                break
      else:
        nIsit=1        

      if nIsit==0:
        if npaMonoLabel!=np.array([]):
          npaMonoLabel=np.append(npaMonoLabel,np.array([npaLabel[nIdx,:4]]),0)
        else:
          npaMonoLabel=np.array([npaLabel[nIdx,:4]])
    return npaMonoLabel


#monopole term calculation works with p-n formalism matrix element labels
  def calcMonoOcc(self):
    sLevName=self.getLevName()
    npaLabel=self.getLabel()
    npaOcc=self.getOcc(sLevName)
    import numpy as np
    nMonoSize=self.getMonoME().size    
    npaMono=np.zeros([npaOcc.shape[0],nMonoSize])
#    print npaOcc.shape[0], nMonoSize

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
        if npaLabel[nIdx,1]==npaLabel[nIdx,3] and npaLabel[nIdx,0]==npaLabel[nIdx,2] and npaLabel[nIdx,4]!=0:
          temp=float(npaLabel[nIdx,4]*(npaLabel[nIdx,4]+1))
          if npaLabel[nIdx,0]!=npaLabel[nIdx,1]:
            npaMono[nLevIdx,nMono]+=npaOcc[nLevIdx,int(npaLabel[nIdx,0]-1)]*npaOcc[nLevIdx,int(npaLabel[nIdx,1]-1)]*temp
          else:
            npaMono[nLevIdx,nMono]+=npaOcc[nLevIdx,int(npaLabel[nIdx,0]-1)]*(npaOcc[nLevIdx,int(npaLabel[nIdx,1]-1)]-1.0)*temp/2.0
          denom[nLevIdx,nMono]+=temp
          nMono+=1
    npaMono=np.divide(npaMono,denom)
#    print npaMono.max()
#    print npaMono
    temp=np.append(npaOcc,npaMono,axis=1)
    return temp    
  
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
            temp=line[5:11]
            temp=[temp]
            if npaOcc!=[]:              
              npaOcc=np.append(npaOcc, temp, axis=0)              
            else:
              npaOcc=temp

      fIn.close()
      return np.array(npaOcc,dtype=float)
      
