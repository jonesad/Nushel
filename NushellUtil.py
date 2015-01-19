# -*- coding: utf-8 -*-
"""
Created on Wed Jan 07 14:06:14 2015

@author: jonesad
contains functions generally useful for manipulating the shell model files
"""

#function for calculating the rms deviation of nNumlev lowest levels (excluding 
#gs) of nushell ouput file at sPath.
def FileRMS(sPathRes,sPathExp,nNumlev):
    import numpy as np
    import math
    fLevelsRes=open(sPathRes,'r')
    fLevelsExp=open(sPathExp,'r')
    nRownum=0    
    nIdx=0
    xERes=[]
    xEExp=[]
    for sRow in fLevelsRes:
      nRownum=nRownum+1
      if nIdx<=nNumlev:
        if nRownum>6:#skip header
          saRowray=sRow.strip().split()
          xERes.append(saRowray[3])
          nIdx=nIdx+1;
    #Reinitialize variables
    nRownum=0    
    nIdx=0
    
    for sRow in fLevelsExp:
      nRownum=nRownum+1
      if nIdx<=nNumlev:
        if nRownum>1:#skip header
          saRowray=sRow.strip().split()
          xEExp.append(saRowray[0])
          nIdx=nIdx+1;
    xERes=np.array(xERes,float)
    xEExp=np.array(xEExp,float)
    return math.sqrt(np.dot(xERes-xEExp,xERes-xEExp)/float(nNumlev))
 
#line of test code   
#print FileRMS("C:\\PythonScripts\\NushellScripts\\o_20b.lpt","C:\\PythonScripts\\NushellScripts\\o_020exp.lpt",10)
 
#Function for editing the interaction file at sPath for the shell model code. 
#Copies original file then changes all matrix element types specified by the 
#integer array anBody (1 for one body and or 2 for 2 body) to the approprate 
#value in the the float array afME. foreknowledge of the number of matrix 
#elements to be changed is necessary. The bRm variable is optional with the 
#default value 0, corresponding to leaving a copy of the original matrix 
#elements in the directory with an appended "_0". Written to work with the USDB
#interaction. Will likely need modifications for other interactions. 
 
def ModInt(sPath, anBody, afOBME=[],afTBME=[],bRm=0):
  import shutil  
  nCh=0 #number of changes to the interactrion file
  for nSpec in anBody:
    if nSpec==1 and len(afOBME)!=0:
      shutil.copyfile(sPath, sPath+"_"+str(nCh))
      fIntSrc=open(sPath+"_"+str(nCh),'r') 
      fIntOut=open(sPath,'w')
      nCh=nCh+1
      nLIdx=0
      for line in fIntSrc:
        nLIdx=nLIdx+1
        if nLIdx==7:
          line=line.strip().split()
          sStart="{0:3}{1:>11.4f}"
          sNext="{0:>10.4f}"
          sNew=sStart.format(line[0],float(afOBME[0]))
          nElemIdx=0          
          for elem in afOBME:
            if nElemIdx!=0:
                sNew=sNew+sNext.format(float(elem))
            nElemIdx=nElemIdx+1
          for i in range(len(line)-len(afOBME)-1):
            sNew=sNew+sNext.format(float(line[len(afOBME)+i+1]))
          fIntOut.write(sNew+"\n")
        else:
          fIntOut.write(line)
      fIntSrc.close()
      fIntOut.close()
    elif nSpec==2 and len(afTBME)!=0:
      shutil.copyfile(sPath, sPath+"_"+str(nCh))
      fIntSrc=open(sPath+"_"+str(nCh),'r') 
      fIntOut=open(sPath,'w')
      nCh=nCh+1
      nNumElem=len(afTBME)
      nLIdx=0
      nElem=0
      for line in fIntSrc:
        nLIdx=nLIdx+1
        if nLIdx>7 and nElem<nNumElem:
          line=line.strip().split()
          sStart="{0:>3}"
          sNext="{0:>5}{1:>3}{2:>9.4f}"
          sNew=""          
          for i in range(4):
            sNew=sNew+sStart.format(line[i])
          sNew=sNew+sNext.format(line[4],line[5], float(afTBME[nElem]))
          fIntOut.write(sNew+"\n")          
          nElem=nElem+1
        else:
          fIntOut.write(line)
        #end if
      #end for loop over fIntSrc      
      fIntSrc.close()
      fIntOut.close()
    #end if
  #end for loop over anBody
  if bRm==1:#remove file backup if requested
    from os import remove as rm
    for i in range(nCh):
      rm(sPath+"_"+str(i))

#line of test code
#ModInt("c:\\PythonScripts\\NushellScripts\\usdb.int", [1,2], [-1.111145181,2.222288212,-3.12345678],[4.123456789,5.111111111,-6.111100998,7.777788902374960])

#Function to retrieve the matrix elements from a Nushell interaction file. 
#Written to work with the usdb interaction file needs modifications for 
#other interactions. Returns the specified matrix elements only (nBody=1 for 
#One Body  and nBody=2 for Two body).
def ReadInt(sPath, nBody):
  fIntInput=open(sPath,'r')
  nLIdx=0
  out=[]
  for line in fIntInput:
    nLIdx=nLIdx+1
    line=line.strip().split()
    if nBody==1 and nLIdx==7:
      out=line[1:4]
    elif nBody==2 and nLIdx>7:
      out.append(line[6])
  return out

#2 lines of test code below
#print ReadInt("c:\PythonScripts\NushellScripts\usdb.int", 1)
#print ReadInt("c:\PythonScripts\NushellScripts\usdb.int", 2)  