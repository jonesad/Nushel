# -*- coding: utf-8 -*-
"""
Created on Mon Jan 19 16:47:49 2015

class for managing a single nucleus

@author: jonesad
"""

class nucleus:
  'class for managing a single nucleus'
#take the path and create a sub directory for the nucleus and run a default calculation for the nucleus.
  def __init__(self,  nZ, nA,sPath, sMMDR, sPar, lsShared):
    import os
    self.msPath=sPath
    self.nAZ=[nA,nZ]    
    self.msName='A'+str(nA)+'_Z'+str(nZ)
    if not os.path.exists(sPath+'\\'+self.msName):
      os.makedirs(sPath+'\\'+self.msName)
    self.writeAns(sMMDR, sPar,lsShared)
    print os.getcwd()
    os.chdir(sPath+'\\'+self.msName)
    print os.getcwd()    
    os.system('shell '+self.msName+'.ans')
    os.system(self.msName)    
#mak an .ans file for use with Nushellx   
  def writeAns(self, sMMDR, sPar,lsShared):
    fAns=open(self.msPath+'\\'+self.msName+'\\'+self.msName+'.ans','w')
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
    self.llSpec=llSpec
    