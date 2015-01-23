# -*- coding: utf-8 -*-
"""
Created on Fri Jan 09 16:03:30 2015

@author: jonesad

functions for optimizing USDB
"""
#Take input matrix elements chan interaction file run the shell model code and 
#output rms energy. 
def RunSm(sBatPath, sIntPath, sEnPath, sEnPathEx, afME, nOBME):
  import NushellUtil
  from os import system
  anBody=[]
  sIntPath
  if afME.size>nOBME:
    anBody=[1,2]
  else:
    anBody=[1]
  NushellUtil.ModInt(sIntPath, anBody, afME[:nOBME],afME[nOBME:])
  system(sBatPath)
  return NushellUtil.FileRMS(sEnPath,sEnPathEx, 2)

#testcode just runs a simple test of the function RunSm 
#def testcode():
#  import NushellUtil
#  import numpy as np
#  ME1=NushellUtil.ReadInt("usdb.int", 1)
#  nOBME=len(ME1)
#  ME2=NushellUtil.ReadInt("usdb.int", 2)
#  for elem in ME2:  
#    ME1.append(elem)
#  afME=np.array(ME1,dtype=float)  
#  #print afME  
#  print RunSm("c:\\rsh-nushellx\\test\\testit.bat", "c:\\rsh-nushellx\\test\\usdb.int", "c:\\rsh-nushellx\\test\\o_20b.lpt", "c:\\rsh-nushellx\\test\\o_020exp.lpt", afME, nOBME)
#  #raw_input("Press any key...")
#testcode()

def OptUSDB(sBatPath, sIntPath, sEnPath, sEnPathEx):
  from scipy.optimize import minimize
  import numpy as np
  import NushellUtil
  import shutil  
  shutil.copy( sIntPath,sIntPath+"_start")
  ME1=NushellUtil.ReadInt(sIntPath, 1)
  nOBME=len(ME1)  
  fun=lambda x: RunSm(sBatPath, sIntPath, sEnPath, sEnPathEx, x, nOBME)  
#  ME2=NushellUtil.ReadInt(sIntPath, 2)
#  for elem in ME2:  
#    ME1.append(elem)
  afME=np.array(ME1,dtype=float)  
  res = minimize(fun, afME, method='Nelder-Mead',options={'disp': True, 'initial_trust_radius' : 0.1})
  return res
  
#test code runs simple test of OptUSDB
#OptUSDB( "c:\\rsh-nushellx\\test\\testit.bat", "c:\\rsh-nushellx\\test\\usdb.int", "c:\\rsh-nushellx\\test\\o_20b.lpt","c:\\rsh-nushellx\\test\\o_020exp.lpt")  

#iterative linear least square optimization scheme for use with shell model code.
#optimizes only the spe for first few energies.  
def ItLLSqOpt(sBatPath, sIntPath, sEnPath, sEnPathEx, sOccPath):
  fResOld=1.0
  fResNew=0.0
  fun=lambda x: RunSm(sBatPath, sIntPath, sEnPath, sEnPathEx, x, nOBME)
  import NushellUtil
  ME1=NushellUtil.ReadInt(sIntPath, 1)
  nOBME=len(ME1)
  nIter=0
  while abs(fResOld-fResNew) >= 0.00001 :  
    fResOld=fResNew
    npaOBME=NushellUtil.singleParticleLeastSq(sOccPath, sEnPathEx, 2)
    fResNew=fun(npaOBME)
    nIter+=1
    print "\n\n Number of iterations is: "+str(nIter)
    print " Res is: " +str(fResNew)
    print" Stopping Criterion is "+str(abs(fResOld-fResNew))+"\n\n"
    
  return fResNew
#test code for the Iterative llsq
import sys
sys.path.append("c:\PythonScripts\NushellScripts")
print ItLLSqOpt("c:\\rsh-nushellx\\test\\testit.bat", "c:\\rsh-nushellx\\test\\usdb.int", "c:\\rsh-nushellx\\test\\o_20b.lpt","c:\\rsh-nushellx\\test\\o_020exp.lpt","c:\\rsh-nushellx\\test\\o_20b.occ")