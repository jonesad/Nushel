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
  #monopole
  fInFile.write('MONO\n')
#  #Single particle
#  fInFile.write('OBME\n')
#  #Two particle
#  fInFile.write('TBME\n')

  #model space specification
  fInFile.write('sdpn\n')
  #restriction
  fInFile.write('n\n')
  #interaction
  fInFile.write('usdcpn\n')
  #formalism iso/pn
  fInFile.write('pn\n')
  #does the interaction extrapolate matrix elements
  fInFile.write('False\n')
  
  #number of nuclei used in optimization
  fInFile.write('8\n')

  #number of protons 
  fInFile.write('  8\n')
  #of nucleons
  fInFile.write('17\n')
  #Experimental ground state energy of the nucleus
  fInFile.write('-4.14308\n')
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
  fInFile.write('-12.18845\n')
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
  fInFile.write('-16.14409\n')
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
  fInFile.write('-23.752104\n')
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
  fInFile.write('-34.40758\n')
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
  fInFile.write('-37.141572\n')
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
  fInFile.write('-41.333144\n')
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
  fInFile.write('-41.243138\n')
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
   def __init__(self, sInPath, sOutPath, sErrorPath='', initialize=True,conservative=False):
#flag that says to use the groundstate energy in optimization
      self.init=False
      self.useGS=1
#flag that tracks the residue and energy levels over the iterations
      self.sMethod=''
      self.track=1
      self.sInPath = sInPath
      self.sOutPath = sOutPath
      #get info from input file and initializes the nuclei objects
      self.GetIn(initialize)
      import os
      if not os.path.exists(self.sOutPath+'\\'+'tracking'):
        os.makedirs(self.sOutPath+'\\'+'tracking')
      if initialize:
        self.writeLevs(self.sOutPath+'\\'+'tracking\\')
      self.conservative=conservative
      if sErrorPath!='':
        import numpy as np
        self.npaErrors=np.array([])
        for nucleus in self.mloNuclei:
          self.npaErrors=np.append(self.npaErrors,nucleus.getLevError(sErrorPath))
      if initialize:
        self.obj(self.mloNuclei[0].getME())
      self.init=True
#        self.npaErrors.shape=[1,self.npaErrors.size]


#Write the initial state of the fit          
   def writeLevs(self, path):
     npaEExp=[]
     npaETh=[]
     lLS=[]
     lAZ=[]
     for nucleus in self.mloNuclei:
       import numpy
       temp=list(numpy.array(nucleus.getEExp(),dtype=float))
       tempth=list(nucleus.getEnNu())
       lLS.extend(nucleus.mllspec)
       lAZ.extend([nucleus.nAZ]*len(nucleus.mllspec))
       if npaEExp!=[]:
         npaEExp=numpy.append(npaEExp,temp,axis=0)
         npaETh=numpy.append(npaETh,tempth,axis=0)
       else: 
         npaEExp=list(temp)
         npaETh=list(tempth)
         
     fOut=open(path+'Levels.dat','w')
     sHFormat='{:10}{:5}{:5}{:5}{:10}{:10}'
     fOut.write(sHFormat.format('[A,Z]','J','nJ','P', 'Eexp', 'Einit'))
     sFormat='{:10}{:5}{:5}{:5}{:10.4f}{:10.4f}'
#     print lAZ, lLS, npaEExp.shape, npaETh.shape
     for AZ,LS,Exp,Th in zip(lAZ, lLS,npaEExp,npaETh):
       temp=LS.strip().split()
       fOut.write(sFormat.format(AZ,temp[0],temp[1],temp[2],Exp,Th)+'\n')
     fOut.close()

#Update the levels with their final values
   def updateLevs(self, path):
     npaETh=[]
     lLS=[]
     lAZ=[]
     for nucleus in self.mloNuclei:
       import numpy
       tempth=list(nucleus.getEnNu())
       lLS.extend(nucleus.mllspec)
       lAZ.extend([nucleus.nAZ]*len(nucleus.mllspec))
       if npaETh!=[]:
         npaETh=numpy.append(npaETh,tempth,axis=0)
       else: 
         npaETh=list(tempth)
     import os
     if not os.path.isfile(path+'Levels_.dat'):
       os.rename(path+'Levels.dat',path+'Levels_.dat')
     fIn=open(path+'Levels_.dat','r')
     fOut=open(path+'Levels.dat','w')
     sHFormat='{:10}{:5}{:5}{:5}{:10}{:10}{:10}\n'
     temp=fIn.readline().strip().split()
     temp.append('Efinal')
     fOut.write(sHFormat.format(*temp))
     sFormat='{:10}{:5}{:5}{:5}{:10}{:10}{:10}\n'
     nIdx=0
     for line in fIn:
       temp=line.strip().split()
       temp.append(npaETh[nIdx])
       fOut.write(sFormat.format(*temp))
       nIdx+=1
     fIn.close()
     fOut.close()
     os.remove(path+'Levels_.dat')
            
   def testMERW(self):
     self.mloNuclei[0].takeME([1,2,3])       
   
   def IterativeLSq(self, sMethod='single', fTolin=10**-3, nMaxIter=100, bMix=True): 
     self.sMethod=sMethod
#     fResLast=0
     # store the original matrix element specificaton so it can be restored 
     #after it is altered.
     llOriginalMESpec=list(self.mloNuclei[0].llMESpec)

     fResNew=100
     lRes=[100,1000, 10000, 100000]
     nIter=0
     fTol=fTolin/float(len(self.mloNuclei[0].mllspec))
     lME=[[],[],[],[]]
     
     if sMethod=='mono' or sMethod=='smono':
       npaMonoLabel=self.mloNuclei[0].getMonoLabel()
       for nucleus in self.mloNuclei:
         nucleus.llMESpec=[list(range(1,nucleus.countOBME()+1)),npaMonoLabel]
#         print nucleus.llMESpec
         nucleus.setmanBody([1, 2])
     while abs(lRes[0]-lRes[1])>fTol and nIter<nMaxIter:
#       fResLast=fResNew
       for nIdx in range(3):
        lRes[len(lRes)-(nIdx+1)]=lRes[len(lRes)-(nIdx+2)]
        lME[len(lRes)-(nIdx+1)]=lME[len(lRes)-(nIdx+2)]
       if sMethod=='single':
         temp=self.singleParticleLeastSq()
         npaGuess=temp[0]
       elif sMethod=='mono':
         temp=self.monopoleLeastSq()
         npaGuess=temp[0]
       elif sMethod=='smono':
         temp=self.summedMLS()
         npaGuess=temp[0]
                 
       import numpy

# test to see if the linear optimization step is too big 
       print "The guess is",npaGuess
       comp1=numpy.linalg.norm(self.mloNuclei[0].getME())/float(self.mloNuclei[0].getME().size)
       comp2=numpy.linalg.norm(npaGuess)/float(npaGuess.size)
#       print comp1, comp2
       if comp1<comp2 and self.conservative==True:
         num=0.5*comp1/comp2
#         print num
#         print npaGuess
         npaGuess=num*npaGuess + (1-num)*self.mloNuclei[0].getME()
#         print npaGuess
       if bMix:
         if lME[1]!=[]and lME[2]!=[] and numpy.all(lME[1].shape==npaGuess.shape) and numpy.all(lME[2].shape==npaGuess.shape):
           num=numpy.random.rand()*0.25 + 0.25
           npaGuess=0.5*npaGuess +num*lME[1]+(0.5-num)*lME[2]          
         elif lME[1]!=[] and numpy.all(lME[1].shape==npaGuess.shape):
           num=numpy.random.rand()*.5 + 0.5
           npaGuess=num*npaGuess +(1.-num)*lME[1]         
         
       fResNew=self.obj(npaGuess)
#       raw_input("Press Enter to continue...")
       lRes[0]=fResNew
       lME[0]=npaGuess
       nIter+=1
     
       #restore the original MESpec 
       for nucleus in self.mloNuclei:    
         nucleus.llMESpec=list(llOriginalMESpec)

#       print lRes
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
     self.updateLevs(self.sOutPath+'\\'+'tracking\\')
     return fResNew
   
   def performOptimization(self, sMethod='Nelder-Mead', dOptions=None):
      npaGuess=self.mloNuclei[0].getME()
      from scipy.optimize import minimize
      res = minimize(self.obj, npaGuess, method=sMethod, options=dOptions)
      print res
      self.updateLevs(self.sOutPath+'\\'+'tracking\\')
      
   def GetIn(self, initialize=True):
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
           print "Error: invalid specification of matrix element category"
       else: 
         sTempL=fIn.readline().strip('\n')
         llMESpec.append(sTempL)
#         print sTempL
     self.lsShared=[]
     
     for nIdx in range(3):
       self.lsShared.append(fIn.readline().strip('\n'))
     self.sForm=fIn.readline().strip('\n')
     self.bExtrap=bool( eval(fIn.readline().strip('\n')))
#     print self.bExtrap
#     print self.sForm
     self.mnNuclei=int(fIn.readline().strip('\n'))
     self.mloNuclei=[]
     for nIdx in range(self.mnNuclei):
       self.mloNuclei.append(ShellNuclei.nucleus(int(fIn.readline()),int(fIn.readline()),float(fIn.readline()),self.sOutPath,fIn.readline().strip('\n'),fIn.readline().strip('\n'),self.lsShared,llMESpec,self.useGS, self.sForm,initialize,self.bExtrap))
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
       if self.init==True:
         oNuc.takeME(npaME)
         if self.sForm=='pn':
           import os
           sLevName=oNuc.getLevName()
           temp=str(oNuc.sPath)+'\\'+str(oNuc.sName)+'\\'+sLevName[:-5]+'0'+'.int'      
           if os.path.isfile(temp+'_'):            
             os.remove(temp+'_')
           os.rename(temp, temp+'_') 
           
#           raw_input("Press Enter to continue...")
         #run Shell model calc
#       print 'running calculation...'
         oNuc.runSM()
         #get the energy difference for the releveant particles
#       print 'geting ediff'
       temp=oNuc.Ediff(bTrackDiff=True)
       for elem in temp:
         res.append(elem)

#     print res
#     raw_input("Press Enter to continue...")
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
         string=string+sFormME.format(float(elem))
       fOut.write(string+'\n')
       fOut.close()

       fOut=open(self.sOutPath+'\\'+'tracking'+'\\AllME.dat','a+')
       npaAllME=self.mloNuclei[0].getME(bAll=True)
       sFormME='{:10.4f}\t'
       string=''
       for elem in npaAllME:
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

     import sys
     sys.path.append('C:\PythonScripts\generalmath')
     import MatManip
#get the sero cols of the matrix
     rmList=MatManip.getZeroCols(a)
#add to rmList the cols for poorly determined ME
#     temp=self.getBadCol(a)
#     rmList.extend(temp)
#     rmList=sorted(list(set(rmList)))
#     print rmList
     if rmList!=[]:
       a=MatManip.rmSlice(rmList, a, 1)

       lLabSpec=self.constructLabels(rmList)
       lsOBLab=[]
       npaTBLab=[]
       print 'lab', lLabSpec
       for elem in lLabSpec:
#         print elem
         if len (elem)==1:
           lsOBLab.append(elem[0])
         else:
           if npaTBLab!=[]:
             temp=numpy.array(elem,dtype=int)
             temp.shape=[1,temp.size]
             npaTBLab=numpy.append(npaTBLab,temp, axis=0)
           else:
             npaTBLab=numpy.array(elem, dtype=int)
             npaTBLab.shape=[1, npaTBLab.size]
#       print lsOBLab
       for nucleus in self.mloNuclei:
         nucleus.llMESpec[0]=lsOBLab
         nucleus.setMEnum()
     obme=self.mloNuclei[0].getOBME()
     print 'mesp', self.mloNuclei[0].llMESpec
     obme.shape=[obme.size,1]
     npaETh.shape=[npaETh.size,1]
     Hexpect=npaETh-numpy.dot(a,obme)
#     print 'Hexpect.shape',Hexpect.shape
     npaEExp.shape=[npaEExp.size,1]
     target=npaEExp-Hexpect

     npaWeights=numpy.zeros([npaEExp.size,npaEExp.size])
     for nIdx, elem in enumerate(self.npaErrors):
#       print elem
       npaWeights[nIdx,nIdx]=1.0/(elem**2)

     a=numpy.dot(npaWeights, a)
     target=numpy.dot(npaWeights, target)

     ans=numpy.linalg.lstsq(a, target)
     ans=ans[0]
     
     return ans, a, target, obme
         
         
#returns the single particle energy + monopole lest square solution to the energy
   def monopoleLeastSq(self):
     import numpy
     npaMono=[]
     npaSPOcc=[]
     npaEExp=[]
     npaETh=[]
     for nucleus in self.mloNuclei:
       temp=nucleus.calcMonoOcc()
       if npaMono!=[]:
         npaMono=numpy.append(npaMono,temp,axis=0)
       else:
         npaMono=numpy.array(temp)
       temp=numpy.array(nucleus.getEExp(),dtype=float)
       tempth=nucleus.getEnNu()
       if npaEExp!=[]:
         npaEExp=numpy.append(npaEExp,temp,axis=0)
         npaETh=numpy.append(npaETh,tempth,axis=0)
       else: 
         npaEExp=temp
         npaETh=tempth
       temp=nucleus.getReducedOcc()
#       print temp
       if npaSPOcc!=[]:
         npaSPOcc=numpy.append(npaSPOcc,temp,axis=0)
       else:
         npaSPOcc=numpy.array(temp)
#     print npaMono

     import sys
     sys.path.append('C:\PythonScripts\generalmath')
     import MatManip
#     if self.sForm=='iso':
#       nOBME=self.mloNuclei[0].countOBME()
#       temp=numpy.array(npaSPOcc[:,:nOBME])
#       print numpy.array(npaSPOcc[:,:nOBME]).shape, numpy.array(npaSPOcc[:,nOBME:]).shape,npaSPOcc.shape
#       temp+=numpy.array(npaSPOcc[:,nOBME:])
#       npaSPOcc=numpy.array(temp)
#     
     SPRMList=MatManip.getZeroCols(npaSPOcc)
     MonoRMList=MatManip.getZeroCols(npaMono)
#     print npaMono
     if SPRMList!=[]:
       npaSPOcc=MatManip.rmSlice(SPRMList,npaSPOcc, 1)
       npaMono=MatManip.rmSlice(MonoRMList,npaMono, 1)
       for nucleus in self.mloNuclei:
#         print nucleus.llMESpec[0]
#         nucleus.llMESpec[0]=list(MatManip.rmSlice(SPRMList,numpy.array(nucleus.llMESpec[0]),0))
          nucleus.llMESpec[0]=[]
#         print nucleus.llMESpec[0]
#         print nucleus.llMESpec[1] 
          nucleus.llMESpec[1]=MatManip.rmSlice(MonoRMList,nucleus.llMESpec[1],0)
#         print nucleus.llMESpec[1] 
          nucleus.setMEnum()
     
     npaME=self.mloNuclei[0].getME()
     print npaME
#     print self.mloNuclei[0].manBody
#     print '\n',npaEExp.shape, a.shape,npaME.shape, numpy.dot(a,npaME).shape,'\n'
#     print npaSPOcc,npaMono
#     a=numpy.append(npaSPOcc,npaMono,1)
     a=npaMono

#     print a
     target=npaEExp-(npaETh-numpy.dot(a,npaME))
#     print "target is ",target
     npaWeights=numpy.zeros([npaEExp.size,npaEExp.size])
     for nIdx, elem in enumerate(self.npaErrors):
#       print elem
       npaWeights[nIdx,nIdx]=1.0/(elem**2)
#     ans=numpy.linalg.lstsq(numpy.dot(npaWeights, a), numpy.dot(npaWeights, target))
     ans=numpy.linalg.lstsq(a, target)
     ans=ans[0]

     return ans, a, target, npaME
 
#returns the single particle energy + monopole lest square solution to the energy
   def summedMLS(self):
     import numpy
     a=[]
     npaEExp=[]
     npaETh=[]
     for nucleus in self.mloNuclei:
       temp=nucleus.summedMO()
       npaNewLabels=temp[1]
#       print temp.shape
       if a!=[]:
         a=numpy.append(a,temp[0],axis=0)
       else:
         a=numpy.array(temp[0])
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
#get the sero cols of the matrix
     rmList=MatManip.getZeroCols(a)
#add to rmList the cols for poorly determined ME
#     temp=self.getBadCol(a)
#     rmList.extend(temp)
#     rmList=sorted(list(set(rmList)))
     if rmList!=[]:
       a=MatManip.rmSlice(rmList, a, 1)
       lLabSpec=self.constructML(rmList,npaNewLabels)
       lsOBLab=[]
       npaTBLab=[]
       for elem in lLabSpec:
         if len (elem)==1:
           lsOBLab.append(elem[0])
         else:
           if npaTBLab!=[]:
             temp=numpy.array(elem,dtype=int)
             temp.shape=[1,temp.size]
             npaTBLab=numpy.append(npaTBLab,temp, axis=0)
           else:
             npaTBLab=numpy.array(elem, dtype=int)
             npaTBLab.shape=[1, npaTBLab.size]
     npaME=self.mloNuclei[0].getME()
#     print self.mloNuclei[0].manBody
#     print '\n',npaEExp.shape, a.shape,npaME.shape, numpy.dot(a,npaME).shape,'\n'
     target=npaEExp-(npaETh)
#     print "target is ",target
     npaWeights=numpy.zeros([npaEExp.size,npaEExp.size])
     for nIdx, elem in enumerate(self.npaErrors):
#       print elem
       npaWeights[nIdx,nIdx]=1.0/(elem**2)
#     ans, errors=MatManip.weightedlsq(a, npaWeights, target)
#     print npaWeights.shape, a.shape, target.shape
#     ans=numpy.linalg.lstsq(numpy.dot(npaWeights, a), numpy.dot(npaWeights, target))
     diff=numpy.linalg.lstsq(a, target)
     diff=diff[0]
#     print  a.shape, target.shape, npaNewLabels.shape
     ans, temp=self.mloNuclei[0].addDiff(diff, npaTBLab)
     for nucleus in self.mloNuclei:
      nucleus.llMESpec[1]=temp
     return ans, a, target, npaME

#quickly set manbody variable 
   def initmanbody(self):
     for nucleus in self.mloNuclei:
       nucleus.setmanBody([1,2])

#Get cols of poorly determinied ME
   def getBadCol(self,npaMat):
     import numpy as np
     import MatManip
     
     rmList=MatManip.getZeroCols(npaMat)
     npaMat2=MatManip.rmSlice(rmList, npaMat, 1)
     
     npaME=self.mloNuclei[0].getME()
     rmList=MatManip.getRepeatSlices(npaMat2,1)
     npaMat3=MatManip.rmSlice(rmList,npaMat2,1)
     npaME2=MatManip.rmSlice(rmList,npaME,0)
    
     psi=np.dot(np.transpose(npaMat3), npaMat3)
     npaErr=np.diagonal(np.linalg.inv(psi))      
     
     lCond=[me/err for me, err in zip(npaME2,npaErr)]
     rmList=[nIdx for nIdx, cond in enumerate(lCond) if abs(cond)>1]
     npaME3=[npaME2[nIdx] for nIdx in rmList]
#     print lCond
     rmList=[nIdx for nIdx, me in enumerate(npaME) if me in npaME3]
     return rmList
#calculate linear fit errors for the matrix elements
   def calcError(self):
    import numpy as np
    import MatManip
    if self.sMethod=='mono':
      ans, a, target, npaME=self.monopoleLeastSq()
    elif self.sMethod=='smono':
      ans, a, target, npaME=self.monopoleLeastSq()
    elif self.sMethod=='single':
      ans, a, target, npaME=self.singleParticleLeastSq()
      
    rmList=MatManip.getRepeatSlices(a,1)
    ap=MatManip.rmSlice(rmList,a,1)
    psi=np.dot(np.transpose(ap), ap)
    print npaME
    lME=MatManip.rmSlice(rmList,npaME, 0)
    print lME
    ans=MatManip.rmSlice(rmList,ans, 0)
    return np.sqrt(np.diagonal(np.linalg.inv(psi))), lME, rmList, ans
    
 
# plot the residual and the matrix elements as a funtion of iteration number
   def plotResults(self, sMethod='single', bError=False):
     self.sMethod=sMethod
     self.initmanbody()
     lMEoriginal=list(self.mloNuclei[0].llMESpec)
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
     plt.plot()
     plt.title('RMS Energy Error by Iteration Number')
     plt.xlabel('Number of Iterations')
     plt.ylabel('RMS Energy Error (MeV)')
     plt.show()
            
     fMEIn=open(self.sOutPath+'\\'+'tracking'+'\\ME.dat','r')
     npaME=[]
     for line in fMEIn:
       temp=line.strip().split()
       if len(npaME)!=0 and len(npaME[0])== len(temp):
         npaME.append(temp)
       else:
         npaME=[temp]
     npaME=np.array(npaME)
     import MatManip
     import Discrete
     npaErr,lME, rmList, npaAns=self.calcError()
     print "Uncertaunty in ME",npaErr
     npaME=MatManip.rmSlice(rmList, npaME,1)
#     rmShortList=[int(elem)-6 for elem in rmList if elem>5]
     # construct labels
     lLabSpec=self.constructResultLabels(rmList)
     lsLabels=[]
#     re initialize the list llmespec
     for nucleus in self.mloNuclei:
       nucleus.llMESpec=list(lMEoriginal)

#     print lLabSpec
     for elem in lLabSpec: 
       print str(type(elem))
       if str(type (elem))=='<type \'numpy.int32\'>':
         lsLabels.append('SPE: '+str(elem))
       else:
         tempstr='[ '
         for nIdx in elem:
           tempstr+=str(nIdx)+' '
         lsLabels.append('TBME: '+tempstr+']')
#     plot the ME
     lBest=Discrete.balFact(npaME.shape[1])
     fig, ax=plt.subplots(nrows=lBest[0], ncols=lBest[1], sharex=True, sharey=False)
     from itertools import chain
     try:
       ax= list(chain.from_iterable(ax))
     except TypeError:
       ''
#       print ax, 'is not iterable'
     
#     print lME 
#     print npaAns
     from itertools import chain

     for nColIdx in range(npaME.shape[1]):
       ax[nColIdx].plot(npaME[:,nColIdx],label=lsLabels[nColIdx])
       if self.sMethod=='single' or self.sMethod=='mono':
         ax[nColIdx].plot([0,npaME.shape[0]-1],[npaAns[nColIdx],npaAns[nColIdx]],ls='-.',color='k', label='Next Linear Fit')
       
       sFormat="Final ME:{:3.2f}$\pm${:2.1e}"
       ax[nColIdx].plot([0,npaME.shape[0]-1],[lME[nColIdx],lME[nColIdx]],ls='--',color='r', label=sFormat.format(float(lME[nColIdx]),float(npaErr[nColIdx])))
       if bError:
         x=[0,int(npaME.shape[0])-1]
         y1=[lME[nColIdx]+npaErr[nColIdx]]*2
#         print 'y1',y1
         y1=list(chain(*y1))         
         y2=[lME[nColIdx]-npaErr[nColIdx]]*2
         
         y2=list(chain(*y2))
#         print y1, y2
         ax[nColIdx].fill_between(x, y1,y2, alpha=0.5, label='Error Band')
#     ax[1].title='ME by Iteration Number'
#       ax[nColIdx].xlabel('Number of Iterations')
#       ax[nColIdx].ylabel('ME (MeV)')
       ax[nColIdx].legend(loc='best')
     plt.show()

#consrtuct a single list of labels for use in the fitting procedure
   def constructLabels(self, rmList):
     from itertools import chain
#     print 'numob',self.mloNuclei[0].countOBME()
     for nucleus in self.mloNuclei:
       nucleus.llMESpec[0]=list(range(1,self.mloNuclei[0].countOBME()+1))
     lLabSpec=list(chain.from_iterable(self.mloNuclei[0].llMESpec))

     if self.sForm=='iso':
       for num in range(3):
         try:
           rmList.remove(num)
         except TypeError:
           ''
     temp=[]
#     print lLabSpec
     for elem in lLabSpec:
       if type(elem).__name__=='int':
#         print type(elem).__name__
         temp.append([elem])
       else:
         temp.append(elem)
#     print len(temp)
     lLabSpec=list(temp)
#     print lLabSpec
#     print 'lab',len(lLabSpec)
#     print rmList
     lLabSpec=list([elem for nIdx, elem in enumerate(lLabSpec) if nIdx not in rmList])
#     print 'lab',len(lLabSpec)
     return lLabSpec

#consrtuct a single list of labels for use in Plotting the results
#the spe are not reinitialized here
   def constructResultLabels(self, rmList):
     from itertools import chain
     lLabSpec=list(chain.from_iterable(self.mloNuclei[0].llMESpec))
     temp=[]
#     print lLabSpec
     for elem in lLabSpec:
       if type(elem).__name__=='int':
# 7        print type(elem).__name__
         temp.append([elem])
       else:
         temp.append(elem)
#     print len(temp)
     lLabSpec=list(temp)
#     print lLabSpec
#     print 'lab',len(lLabSpec)
#     print rmList
     lLabSpec=list([elem for nIdx, elem in enumerate(lLabSpec) if nIdx not in rmList])
#     print 'lab',len(lLabSpec)
     return lLabSpec


     
#consrtuct a single list of labels
   def constructML(self, rmList, labels):
     from itertools import chain
     for nucleus in self.mloNuclei:
       nucleus.llMESpec[0]=list(range(1,self.mloNuclei[0].countOBME()+1))
#     lLabSpec=list(chain.from_iterable(self.mloNuclei[0].llMESpec))
     lLabSpec=list(chain.from_iterable([range(1,self.mloNuclei[0].countOBME()+1),labels]))
     temp=[]
#     print lLabSpec
     for elem in lLabSpec:
       if type(elem).__name__=='int':
#         print type(elem).__name__
         temp.append([elem])
       else:
         temp.append(elem)
#     print len(temp)
     lLabSpec=temp
#     print lLabSpec
#     print 'lab',len(lLabSpec)
     lLabSpec=[elem for nIdx, elem in enumerate(lLabSpec) if nIdx not in rmList]
#     print 'lab',len(lLabSpec)
     return lLabSpec
     
#     add mean gaussian noise to the matrix elements with specified standard deviation
   def addMENoise(self, fStd):
    from numpy.random import randn
    npaME=self.mloNuclei[0].getME()
    npaME+=randn(*npaME.shape)*fStd
    for nucleus in self.mloNuclei: 
      nucleus.takeME(npaME)
           
import sys
sys.path.append('c:\\PythonScripts\\NushellScripts\\')
sys.path.append('C:\PythonScripts\generalmath')

x=ShellOpt('c:\\PythonScripts\\NushellScripts\\OptInput.in','c:\\PythonScripts\\NushellScripts\\test', 'c:\\PythonScripts\\NushellScripts\\errors.dat',initialize=True, conservative=False)
#ans, a, target, npaME=x.monopoleLeastSq()
#err, lME, rmList, ans=x.calcError()


#x.addMENoise(3.0)

print x.IterativeLSq(sMethod='mono',bMix=False, nMaxIter=100, fTolin=10**-5)
#x.performOptimization()
#x.sMethod='single'
#err=x.calcError()[0]
#print err
#nSpecSize=0
#en=x.mloNuclei[1].getEnNu(bAll=True)
#lHist,lBins,lMu,lSigma,nSize=x.mloNuclei[1].calcEThErr(err, nSpecSize,bPrev=True,bAllME=True)
#nStop=8
#x.mloNuclei[1].plotEthError(lHist[:nStop], lBins[:nStop], lMu[:nStop], lSigma[:nStop], nSize)

#x.plotResults(sMethod='mono', bError=True)
