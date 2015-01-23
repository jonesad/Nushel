# -*- coding: utf-8 -*-
"""
Created on Wed Jan 14 14:03:24 2015

attempt at object orientation of nushell optimization

@author: jonesad
"""
#writes an example inputfile for the ShellOpt Class for single particle o20 optimization
def CreateInFile(sInfilePath):
  fInFile=open(sInfilePath,'w')
  # lines in me spec
  fInFile.write('1\n')
  # matrix element specification
  fInFile.write('OBME\n')
  #model space specification
  fInFile.write('sd\n')
  #restriction
  fInFile.write('n\n')
  #interaction
  fInFile.write('usdb\n')
  #number of nuclei used in optimization
  fInFile.write('1\n')
  #number of protons 
  fInFile.write('  8\n')
  #of nucleons
  fInFile.write('20\n')
  #min max delta j
  fInFile.write(' 0.0, 4.0, 1.0,\n')
  #parity
  fInFile.write('  0\n')
  #number of states used in optimization
  fInFile.write('3\n')
  #state specifications
  fInFile.write('  0  1  +1\n')
  fInFile.write('  2  1  +1\n')
  fInFile.write('  4  1  +1\n')    
  
#code to crate an input file
CreateInFile('C:/PythonScripts/NushellScripts/OptInput.in')

class ShellOpt:
   'Class for shell model hamiltonian optimization problems'
   def __init__(self, sInPath, sOutPath):
      self.sInPath = sInPath
      self.sOutPath = sOutPath
      #get info from input file and initializes the nuclei objects
      self.GetIn()
      npaGuess=self.mloNuclei[0].getME()
      from scipy.optimize import minimize
      res = minimize(self.obj, npaGuess, method='COBYLA', options={'tol':0.0001,'rhobeg':0.1,'disp': True})
      print res
      
   def GetIn(self):
     import ShellNuclei
     fIn=open(self.sInPath,'r')
     sTempL=fIn.readline()
     llMESpec=[[],[]]
     for nIdx in range(int(sTempL)):
       if nIdx==0:
         anBody=[]
         sTempL=fIn.readline()
         sTempL=sTempL.split()
         print '\n\n'
         print sTempL
         print '\n\n'
         if sTempL[0]=='OBME':
           anBody.append(1)
           if len(sTempL[1:])!=0:
             llMESpec[0].append(sTempL[1:])
         elif sTempL[0]=='TBME':
           anBody.append(2)
           if len(sTempL[1:])!=0:
             llMESpec[1].append(sTempL[1:])
         else:
           print "Error: invalid specification of matrix element class"
       else: 
         sTempL=fIn.readline().strip('\n')
         llMESpec.append(sTempL)
     self.lsShared=[]
     for nIdx in range(3):
       self.lsShared.append(fIn.readline().strip('\n'))
     self.mnNuclei=int(fIn.readline().strip('\n'))
     self.mloNuclei=[]
     for nIdx in range(self.mnNuclei):
       self.mloNuclei.append(ShellNuclei.nucleus(int(fIn.readline()),int(fIn.readline()),self.sOutPath,fIn.readline().strip('\n'),fIn.readline().strip('\n'),self.lsShared,llMESpec))
       llStateSpec=[]
       temp=int(fIn.readline())
       for iii in range(temp):
         llStateSpec.append(fIn.readline().strip('\n'))
       print llStateSpec
       self.mloNuclei[nIdx].setLevels(llStateSpec)
       self.mloNuclei[nIdx].setmanBody(anBody) 
#       print self.mloNuclei[nIdx].countOBME()
     fIn.close()
 
#The objective function takes the array of matrix elements
   def obj(self,npaME):
     import numpy as np
     res=[]
     #each nucleus
     for oNuc in self.mloNuclei:
       #write ME to file
       oNuc.takeME(npaME)
       #run Shell model calc
       oNuc.runSM()
       #get the energy difference for the releveant particles
       res.append(oNuc.Ediff())
       res=np.array(res)
     return np.sqrt(np.dot(res,res)/float(len(res)))     
import sys

sys.path.append('c:\\PythonScripts\\NushellScripts\\')     
x=ShellOpt('c:\\PythonScripts\\NushellScripts\\OptInput.in','c:\\PythonScripts\\NushellScripts\\test')