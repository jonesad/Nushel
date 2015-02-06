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
  fInFile.write('4\n')

  #number of protons 
  fInFile.write('  8\n')
  #of nucleons
  fInFile.write('17\n')
  #Experimental ground state energy of the nucleus
  fInFile.write('-4.14\n')
  #min max delta j
  fInFile.write(' 0.5, 2.5, 1.0,\n')
  #parity
  fInFile.write('  0\n')
  #number of states used in optimization
  fInFile.write('2\n')
  #state specifications (J, nJ, P)
  fInFile.write('  1/2  1  +1\n')
  fInFile.write('  5/2  1  +1\n')    
  
  #number of protons 
  fInFile.write('  8\n')  
  #of nucleons
  fInFile.write('18\n')
  #Experimental ground state energy of the nucleus
  fInFile.write('-12.19\n')
  #min max delta j
  fInFile.write(' 0.0, 4.0, 1.0,\n')
  #parity
  fInFile.write('  0\n')
  #number of states used in optimization
  fInFile.write('3\n')
  #state specifications (J, nJ, P)
  fInFile.write('  0  1  +1\n')
  fInFile.write('  2  1  +1\n')
  fInFile.write('  4  1  +1\n')    
 
 #number of protons 
  fInFile.write('  8\n')
  #of nucleons
  fInFile.write('20\n')
  #Experimental ground state energy of the nucleus
  fInFile.write('-23.75\n')
  #min max delta j
  fInFile.write(' 0.0, 4.0, 1.0,\n')
  #parity
  fInFile.write('  0\n')
  #number of states used in optimization
  fInFile.write('3\n')
  #state specifications (J, nJ, P)
  fInFile.write('  0  1  +1\n')
  fInFile.write('  2  1  +1\n')
  fInFile.write('  4  1  +1\n')    
  
  #number of protons 
  fInFile.write('  8\n')
  #of nucleons
  fInFile.write('22\n')
  #Experimental ground state energy of the nucleus
  fInFile.write('-34.41\n')
  #min max delta j
  fInFile.write(' 0.0, 4.0, 1.0,\n')
  #parity
  fInFile.write('  0\n')
  #number of states used in optimization
  fInFile.write('1\n')
  #state specifications (J, nJ, P)
  fInFile.write('  0  1  +1\n')
  
#code to create an input file
CreateInFile('C:/PythonScripts/NushellScripts/OptInput.in')

class ShellOpt:
   'Class for shell model hamiltonian optimization problems'
   def __init__(self, sInPath, sOutPath):
#flag that says to use the groundstate energy in optimization
      self.useGS=1
#flag that tracks the residue and energy levels over the iterations 
      self.track=1
      self.sInPath = sInPath
      self.sOutPath = sOutPath
      #get info from input file and initializes the nuclei objects
      self.GetIn()
      import os
      if not os.path.exists(self.sOutPath+'\\'+'tracking'):
        os.makedirs(self.sOutPath+'\\'+'tracking')
   
   def testMERW(self):
     self.mloNuclei[0].takeME([1,2,3])
   
   def IterativeLSq(self, fTol=10**-5, nMaxIter=100): 
     fResLast=0
     fResNew=100
     nIter=0
     while abs(fResLast-fResNew)>fTol and nIter<nMaxIter:
       fResLast=fResNew
       npaGuess=self.singleParticleLeastSq()
       fResNew=self.obj(npaGuess)
       nIter+=1
       print fResNew
     if abs(fResLast-fResNew)<fTol:
       print 'Completed Successfully'
     else:
       print 'Iteration max reached: ', nIter       
     return fResNew
       
   def performOptimization(self, sMethod='Nelder-Mead', dOptions=None):
      npaGuess=self.mloNuclei[0].getME()
      from scipy.optimize import minimize
      res = minimize(self.obj, npaGuess, method=sMethod, options=dOptions)
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
#         print '\n\n'
#         print sTempL
#         print '\n\n'
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
       self.mloNuclei.append(ShellNuclei.nucleus(int(fIn.readline()),int(fIn.readline()),float(fIn.readline()),self.sOutPath,fIn.readline().strip('\n'),fIn.readline().strip('\n'),self.lsShared,llMESpec,self.useGS))
       llStateSpec=[]
       temp=int(fIn.readline())
       for iii in range(temp):
         llStateSpec.append(fIn.readline().strip('\n'))
#       print llStateSpec
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
#       print '[A,Z]=', oNuc.nAZ
       #write ME to file
#       print 'Setting ME:',npaME
       oNuc.takeME(npaME)
       #run Shell model calc
#       print 'running calculation...'
       oNuc.runSM()
       #get the energy difference for the releveant particles
#       print 'geting ediff'
       temp=oNuc.Ediff()
       
       for elem in temp:
         res.append(elem)

     res=np.array(res)     
#     print res
     temp=np.sqrt(np.dot(res,res)/float(len(res)))
     if self.track==1:
       self.OptStatus(temp,npaME)
     return temp   
 #
   def OptStatus(self,res, npaME):
       fOut=open(self.sOutPath+'\\'+'tracking'+'\\res.dat','a+')
       fOut.write(str(res)+'\n')
       fOut.close()
       
       fOut=open(self.sOutPath+'\\'+'tracking'+'\\ME.dat','a+')
       string=''
       for elem in npaME:
         string=string+str(elem)+'\t'
       fOut.write(string+'\n')
       fOut.close()
       
       for oNucleus in self.mloNuclei:
         oNucleus.writeStatus()
   
#returns the single particle energy lest square solution to the energy
   def singleParticleLeastSq(self):
     import numpy
     a=[]
     npaEExp=[]
     npaETh=[]
     for nucleus in self.mloNuclei:
       fIn=open(nucleus.sPath+'\\'+nucleus.sName+"\\list.lpt")
       sLevName=fIn.readline().strip()        
       fIn.close() 
       if a!=[]:
         a=numpy.append(a,nucleus.getOcc(sLevName),axis=0)
       else:
         a=numpy.array(nucleus.getOcc(sLevName))
       temp=numpy.array(nucleus.getEExp(),dtype=float)
       tempth=nucleus.getEnNu()
       if npaEExp!=[]:
         npaEExp=numpy.append(npaEExp,temp,axis=0)
         npaETh=numpy.append(npaETh,tempth,axis=0)
       else: 
         npaEExp=temp
         npaETh=tempth
     obme=self.mloNuclei[0].getOBME()
     Hexpect=npaETh-numpy.dot(a,*obme)
     target=npaEExp-Hexpect
     ans=numpy.linalg.lstsq(a, target)
     ans=ans[0]
     
     return ans#, a, target, obme

import sys

sys.path.append('c:\\PythonScripts\\NushellScripts\\')
x=ShellOpt('c:\\PythonScripts\\NushellScripts\\OptInput.in','c:\\PythonScripts\\NushellScripts\\test')

#ans,a, target, obme=x.singleParticleLeastSq()
print x.IterativeLSq()
#x.performOptimization()