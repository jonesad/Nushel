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
  #Single particle
  fInFile.write('OBME\n')
  #model space specification
  fInFile.write('sdpn\n')
  #restriction
  fInFile.write('n\n')
  #interaction
  fInFile.write('usdbpn\n')
  #formalism iso/pn
  fInFile.write('pn\n')
  
  #number of nuclei used in optimization
  fInFile.write('8\n')

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
  fInFile.write('4\n')
  #state specifications (J, nJ, P)
  fInFile.write('  0  1  +1\n')
  fInFile.write('  2  1  +1\n')
  fInFile.write('  4  1  +1\n')    
  fInFile.write('  3  1  +1\n')    
 
  #number of protons 
  fInFile.write('  8\n')
  #of nucleons
  fInFile.write('19\n')
  #Experimental ground state energy of the nucleus
  fInFile.write('-16.14\n')
  #min max delta j
  fInFile.write(' 0.5, 4.5, 1.0,\n')
  #parity
  fInFile.write('  0\n')
  #number of states used in optimization
  fInFile.write('7\n')
  #state specifications (J, nJ, P)
  fInFile.write('  5/2  1  +1\n')
  fInFile.write('  3/2  1  +1\n')    
  fInFile.write('  1/2  1  +1\n')    
  fInFile.write('  9/2  1  +1\n')    
  fInFile.write('  7/2  1  +1\n')    
  fInFile.write('  5/2  2  +1\n')
  fInFile.write('  3/2  2  +1\n')    
    
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
  fInFile.write('5\n')
  #state specifications (J, nJ, P)
  fInFile.write('  0  1  +1\n')
  fInFile.write('  2  1  +1\n')
  fInFile.write('  4  1  +1\n')    
  fInFile.write('  2  2  +1\n')
  fInFile.write('  0  2  +1\n')
#
#  #number of protons 
#  fInFile.write('  8\n')
#  #of nucleons
#  fInFile.write('21\n')
#  #Experimental ground state energy of the nucleus
#  fInFile.write('-27.56\n')
#  #min max delta j
#  fInFile.write(' 0.5, 4.5, 1.0,\n')
#  #parity
#  fInFile.write('  0\n')
#  #number of states used in optimization
#  fInFile.write('6\n')
#  #state specifications (J, nJ, P)
#  fInFile.write('  5/2  1  +1\n')
#  fInFile.write('  1/2  1  +1\n')    
#  fInFile.write('  3/2  1  +1\n')    
#  fInFile.write('  7/2  1  +1\n')    
#  fInFile.write('  5/2  2  +1\n')    
#  fInFile.write('  9/2  1  +1\n')
#  
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
  
  #number of protons 
  fInFile.write('  8\n')
  #of nucleons
  fInFile.write('23\n')
  #Experimental ground state energy of the nucleus
  fInFile.write('-37.15\n')
  #min max delta j
  fInFile.write(' 0.5, 4.5, 1.0,\n')
  #parity
  fInFile.write('  0\n')
  #number of states used in optimization
  fInFile.write('1\n')
  #state specifications (J, nJ, P)
  fInFile.write('  1/2  1  +1\n')

  #number of protons 
  fInFile.write('  8\n')
  #of nucleons
  fInFile.write('24\n')
  #Experimental ground state energy of the nucleus
  fInFile.write('-41.34\n')
  #min max delta j
  fInFile.write(' 0.0, 4.0, 1.0,\n')
  #parity
  fInFile.write('  0\n')
  #number of states used in optimization
  fInFile.write('1\n')
  #state specifications (J, nJ, P)
  fInFile.write('  0  1  +1\n')

#  #number of protons 
#  fInFile.write('  8\n')
#  #of nucleons
#  fInFile.write('25\n')
#  #Experimental ground state energy of the nucleus
#  fInFile.write('-40.46\n')
#  #min max delta j
#  fInFile.write(' 0.5, 4.5, 1.0,\n')
#  #parity
#  fInFile.write('  0\n')
#  #number of states used in optimization
#  fInFile.write('1\n')
#  #state specifications (J, nJ, P)
#  fInFile.write('  3/2  1  +1\n')

  #number of protons 
  fInFile.write('  8\n')
  #of nucleons
  fInFile.write('26\n')
  #Experimental ground state energy of the nucleus
  fInFile.write('-40.26\n')
  #min max delta j
  fInFile.write(' 0.0, 4.0, 1.0,\n')
  #parity
  fInFile.write('  0\n')
  #number of states used in optimization
  fInFile.write('1\n')
  #state specifications (J, nJ, P)
  fInFile.write('  0  1  +1\n')

#  #number of protons 
#  fInFile.write('  8\n')
#  #of nucleons
#  fInFile.write('27\n')
#  #Experimental ground state energy of the nucleus
#  fInFile.write('-39.09\n')
#  #min max delta j
#  fInFile.write(' 0.5, 2.5, 1.0,\n')
#  #parity
#  fInFile.write('  0\n')
#  #number of states used in optimization
#  fInFile.write('1\n')
#  #state specifications (J, nJ, P)
#  fInFile.write('  3/2  1  +1\n')
  
#  #number of protons 
#  fInFile.write('  8\n')
#  #of nucleons
#  fInFile.write('28\n')
#  #Experimental ground state energy of the nucleus
#  fInFile.write('-38.27\n')
#  #min max delta j
#  fInFile.write(' 0.0, 4.0, 1.0,\n')
#  #parity
#  fInFile.write('  0\n')
#  #number of states used in optimization
#  fInFile.write('1\n')
#  #state specifications (J, nJ, P)
#  fInFile.write('  0  1  +1\n')


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
   
   def IterativeLSq(self, sMethod='single', fTolin=2.5*10**-3, nMaxIter=100): 
     fResLast=0
     fResNew=100
     lRes=[100,1000, 10000, 100000]
     nIter=0
     fTol=fTolin/float(len(self.mloNuclei[0].mllspec))
     lME=[[],[],[],[]]
     if sMethod=='mono':
       npaMonoLabel=self.mloNuclei[0].getMonoLabel()
       
       for nucleus in self.mloNuclei:
         nucleus.llMESpec=[[],npaMonoLabel]
         nucleus.setmanBody([1, 2])
     while abs(lRes[0]-lRes[1])>fTol and nIter<nMaxIter:
#       fResLast=fResNew
       for nIdx in range(3):
        lRes[len(lRes)-(nIdx+1)]=lRes[len(lRes)-(nIdx+2)]
        lME[len(lRes)-(nIdx+1)]=lME[len(lRes)-(nIdx+2)]
       if sMethod=='single':
         npaGuess=self.singleParticleLeastSq()
       elif sMethod=='mono':
         npaGuess=self.monopoleLeastSq()
       import numpy
       import math
       num=numpy.random.rand()*.5 + 0.5
       if lME[1]!=[] and numpy.all(lME[1].shape==npaGuess.shape):  
         npaGuess=num*npaGuess +(1.-num)*lME[1]
         
       fResNew=self.obj(npaGuess)
       lRes[0]=fResNew
       lME[0]=npaGuess
       nIter+=1
       print lRes
       print lME
       print fResNew
       #chek for oscilating convergence and average the ME of oscilations if found 
#       if abs(lRes[0]-lRes[2])<10*fTol and abs(lRes[0]-lRes[1])>10*fTol:
#         npaGuess=(lME[0]+lME[1])*0.5
#         fResNew=self.obj(npaGuess)
#         lRes[0]=fResNew
#         lME[0]=npaGuess
#       
#       elif abs(lRes[0]-lRes[3])<10*fTol:
#         if abs(lRes[0]-lRes[1])>10*fTol:
#           npaGuess=(lME[0]+lME[1])*0.5
#           fResNew=self.obj(npaGuess)
#           lRes[0]=fResNew
#           lME[0]=npaGuess
#       
#         elif abs(lRes[0]-lRes[2])>10*fTol:
#           npaGuess=(lME[0]+lME[2])*0.5
#           fResNew=self.obj(npaGuess)
#           lRes[0]=fResNew
#           lME[0]=npaGuess      
                 
     if abs(lRes[0]-lRes[1])<fTol:
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
         elif sTempL[0]=='MONO':
           anBody.append(3)
         else:
           print "Error: invalid specification of matrix element class"
       else: 
         sTempL=fIn.readline().strip('\n')
         llMESpec.append(sTempL)
         print sTempL
     self.lsShared=[]
     for nIdx in range(3):
       self.lsShared.append(fIn.readline().strip('\n'))
     self.sForm=fIn.readline().strip('\n')
     print self.sForm
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
       sFormME='{:10.4f}\t'
       for elem in npaME:
         string=string+sFormME.format(elem)
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
         
         
#returns the single particle energy + monopole lest square solution to the energy
   def monopoleLeastSq(self):
     import numpy
     a=[]
     npaEExp=[]
     npaETh=[]
     for nucleus in self.mloNuclei:
       temp=nucleus.calcMonoOcc()
#       print temp.shape
       if a!=[]:
         a=numpy.append(a,temp,axis=0)
       else:
         a=numpy.array(temp)
       temp=numpy.array(nucleus.getEExp(),dtype=float)
       tempth=nucleus.getEnNu()
       if npaEExp!=[]:
         npaEExp=numpy.append(npaEExp,temp,axis=0)
         npaETh=numpy.append(npaETh,tempth,axis=0)
       else: 
         npaEExp=temp
         npaETh=tempth
          
#     print npaETh
#     print numpy.dot(a,npaME)
#     print a.shape
#     print npaME.shape
#     remove zero columns of the matrix and the associated matrix elements
     import sys
     sys.path.append('C:\PythonScripts\generalmath')
     import MatManip
     rmList=MatManip.getZeroCols(a)
#     print a.shape
     if rmList!=[]:
       a=MatManip.rmSlice(rmList, a, 1)
       temp1=[]
       temp2=[]
       for elem in rmList:
         if elem<=5:
           temp1.append(elem)
         if elem >5:
           temp2.append(elem-6)
#       print a.shape
#       print len(rmList)
#       print len(temp1), len(temp2)

       for nucleus in self.mloNuclei:
         nucleus.llMESpec[1]=nucleus.getMonoLabel()
         
         nucleus.llMESpec[0]=[ii for ii in range(6)  if ii not in temp1]
         nucleus.llMESpec[1]=MatManip.rmSlice(temp2, nucleus.llMESpec[1],0)
     npaME=self.mloNuclei[0].getOBME()
#     print self.mloNuclei[0].getMonoLabel().shape
#     print self.mloNuclei[0].llMESpec[1].shape
#     print len(self.mloNuclei[0].llMESpec[0])
#     print self.mloNuclei[0].getTBME()
     npaME=numpy.append(npaME,self.mloNuclei[0].getTBME())
     
     target=npaEExp-(npaETh-numpy.dot(a,npaME))
     ans=numpy.linalg.lstsq(a, target)
     ans=ans[0]
     return ans#, a, target, npaME
 
# plot the residual and the matrix elements as a funtion of iteration number
   def plotResults(self):
     fResIn=open(self.sOutPath+'\\'+'tracking'+'\\res.dat','r')       
     npaRes=[]
     for line in fResIn:
       npaRes.append(line.strip())
     fResIn.close()
     import numpy as np
     npaRes=np.array(npaRes)
     import matplotlib.pyplot as plt
     plt.figure()
     plt.plot(npaRes, label='Res')
     plt.title('RMS Energy Error by Iteration Number')
     plt.xlabel('Number of Iterations')
     plt.ylabel('RMS Energy Error (MeV)')
     plt.show()
     
#  in case the files are overwritten run the the mopole least square to rewrite the me labels
     self.monopoleLeastSq()
         
     fMEIn=open(self.sOutPath+'\\'+'tracking'+'\\ME.dat','r') 
     npaME=[]
     for line in fMEIn:
       npaME.append(line.strip().split())
     npaME=np.array(npaME)
     # construct labels
     lsLabels=[]
     for elem in self.mloNuclei[0].llMESpec[0]:
       lsLabels.append('SPE: '+str(elem+1))
     for npaMEIdx in self.mloNuclei[0].llMESpec[1]:
       tempstr='[ '
       for nIdx in npaMEIdx:
         tempstr+=str(nIdx)+' '
       lsLabels.append('TBME: '+tempstr+']')
#     plot the ME
     plt.figure()
     for nColIdx in range(npaME.shape[1]):
       plt.plot(npaME[:,nColIdx],label=lsLabels[nColIdx])
     plt.title('ME by Iteration Number')
     plt.xlabel('Number of Iterations')
     plt.ylabel('ME (MeV)')
#     plt.legend(loc='best')
     plt.show()
          
import sys
sys.path.append('c:\\PythonScripts\\NushellScripts\\')
sys.path.append('C:\PythonScripts\generalmath')

x=ShellOpt('c:\\PythonScripts\\NushellScripts\\OptInput.in','c:\\PythonScripts\\NushellScripts\\test')

#ans, a, target, npaME=x.monopoleLeastSq()
#print x.IterativeLSq(sMethod='mono')
#x.performOptimization()
x.plotResults()