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
  fInFile.write('1')
  # matrix element specification
  fInFile.write('OBME')
  #model space specification
  fInFile.write('sd')
  #restriction
  fInFile.write('n')
  #interaction
  fInFile.write('usdb')
  #number of nuclei used in optimization
  fInFile.write('1')
  #number of protons 
  fInFile.write('  8')
  #of nucleons
  fInFile.write('20')
  #min max delta j
  fInFile.write(' 0.0, 4.0, 1.0,')
  #parity
  fInFile.write('  0')
  #number of states used in optimization
  fInFile.write('3')
  #state specifications
  fInFile.write('  0  +1')
  fInFile.write('  2  +1')
  fInFile.write('  4  +1')    
#test code to crate an input file
CreateInFile()
class ShellOpt:
   'Class for shell model hamiltonian optimization problems'

   def __init__(self, sInPath, sOutPath):
      self.sInPath = sInPath
      self.sOutPath = sOutPath
      self.GetIn
   
   def GetIn(self):
     fIn=open(self.sInPath,'r')
     sTempL=fIn.readline()
     lSpec=[]
     for nIdx in range(int(sTempL)):
       if nIdx==0:
         anBody=[]
         sTempL=fIn.readline()
         if sTempL.strip()=='OBME':
           anBody.append(1)
         elif sTempL.strip()=='TBME':
           anBody.append(2)
         else:
           print "Error: invalid specification of matrix element class"
       else: 
         sTempL=fIn.readline()
         lSpec.append(sTempL)
     self.msMspace=fIn.readline()
     self.msRest=fIn.readline()
     self.msInt=fIn.readline()
     self.mnNuclei=int(fIn.readline())
     anA=[]
     anZ=[]
     asMMDJ=[]
     asPar=[]
     llStateSpec=[]
     for nIdx in range(self.mnNuclei):
       anA.append(int(fIn.readline()))
       anZ.append(int(fIn.readline()))
       asMMDJ.append(fIn.readline())
       asPar.append(fIn.readline())
       llStateSpec.append([])
       for iii in range(int(fIn.readline())):
         llStateSpec[nIdx].append(fIn.readline())
       