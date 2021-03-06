# -*- coding: utf-8 -*-

"""
Created on Wed Jan 14 14:03:24 2015

attempt at object orientation of OxBash optimization

@author: jonesad
"""
# writes an example inputfile for the BashOpt Class for single particle o20
# optimization


def CreateInFile(sInfilePath):
    fInFile = open(sInfilePath, 'w')
# lines in me spec
    fInFile.write('1\n')
# matrix element specification
# monopole
    fInFile.write('MONO\n')
#  #Single particle
#  fInFile.write('OBME\n')
#  #Two particle
#  fInFile.write('TBME\n')

# model space specification
    fInFile.write('sd\n')
# restriction
    fInFile.write('n\n')
# interaction if two interactions present second is the background interaction
    fInFile.write('usda    sdba\n')
# formalism iso/pn
    fInFile.write('iso\n')
# does the interaction extrapolate matrix elements based on mass dependence
    fInFile.write('True\n')
# number of nuclei used in optimization
    fInFile.write('10\n')

# number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('17\n')
# number of valence nucleons
    fInFile.write('1\n')
# Experimental ground state energy of the nucleus
    fInFile.write('  -4.14308   0.00000\n')
# min max delta j
    fInFile.write(' 0.5, 2.5, 1.0,\n')
# parity
    fInFile.write('  0\n')
#   number of states used in optimization
    fInFile.write('2\n')
# state specificationss (J, nJ, P)
    fInFile.write('  5/2  1  +1   0.00000   0.00000\n')
    fInFile.write('  1/2  1  +1   0.87076   0.00010\n')

# number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('18\n')
# number of valence nucleons
    fInFile.write('2\n')
# Experimental ground state energy of the nucleus
    fInFile.write(' -12.18845   0.00000\n')
# min max delta j
    fInFile.write(' 0.0, 4.0, 1.0,\n')
# parity
    fInFile.write('  0\n')
# number of states used in optimization
    fInFile.write('4\n')
# state specifications (J, nJ, P)
    fInFile.write('  0  1  +1   0.00000   0.00000\n')
    fInFile.write('  2  1  +1   1.98207   0.00009\n')
    fInFile.write('  4  1  +1   3.55484   0.00040\n')
    fInFile.write('  3  1  +1   5.7780   0.00060\n')

# number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('19\n')
# number of valence nucleons
    fInFile.write('3\n')
# Experimental ground state energy of the nucleus
    fInFile.write(' -16.14409   0.00264\n')
# min max delta j
    fInFile.write(' 0.5, 4.5, 1.0,\n')
# parity
    fInFile.write('  0\n')
# number of states used in optimization
    fInFile.write('7\n')
# state specifications (J, nJ, P)
    fInFile.write('  5/2  1  +1   0.00000   0.00000\n')
    fInFile.write('  3/2  1  +1   0.09600   0.00050\n')
    fInFile.write('  1/2  1  +1   1.47170   0.00040\n')
    fInFile.write('  9/2  1  +1   2.37150   0.00100\n')
    fInFile.write('  7/2  1  +1   2.77900   0.00090\n')
    fInFile.write('  5/2  2  +1   3.15350   0.00170\n')
    fInFile.write('  3/2  2  +1   3.06740   0.00160\n')

# number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('20\n')
# number of valence nucleons
    fInFile.write('4\n')
# Experimental ground state energy of the nucleus
    fInFile.write('-23.752104   0.00088\n')
# min max delta j
    fInFile.write(' 0.0, 4.0, 1.0,\n')
# parity
    fInFile.write('  0\n')
# number of states used in optimization
    fInFile.write('5\n')
# state specifications (J, nJ, P)
    fInFile.write('  0  1  +1   0.00000   0.00000\n')
    fInFile.write('  2  1  +1   1.67368   0.00015\n')
    fInFile.write('  4  1  +1   3.57000   0.00700\n')
    fInFile.write('  2  2  +1   4.07200   0.00400\n')
    fInFile.write('  0  2  +1   4.45600   0.00500\n')

#  number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('21\n')
# number of valence nucleons
    fInFile.write('5\n')
# Experimental ground state energy of the nucleus
    fInFile.write(' -27.55768   0.01199\n')
# min max delta j
    fInFile.write(' 0.5, 4.5, 1.0,\n')
# parity
    fInFile.write('  0\n')
# number of states used in optimization
    fInFile.write('6\n')
#   state specifications (J, nJ, P)
    fInFile.write('  5/2  1  +1   0.00000   0.00000\n')
    fInFile.write('  1/2  1  +1   1.22000   0.00300\n')
    fInFile.write('  3/2  1  +1   2.13300   0.00500\n')
    fInFile.write('  7/2  1  +1   3.07300   0.01000\n')
    fInFile.write('  5/2  2  +1   3.02600   0.00600\n')
    fInFile.write('  9/2  1  +1   4.92700   0.01200\n')

# number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('22\n')
# number of valence nucleons
    fInFile.write('6\n')
# Experimental ground state energy of the nucleus
    fInFile.write(' -34.40758   0.05692\n')
# min max delta j
    fInFile.write(' 0.0, 4.0, 1.0,\n')
# parity
    fInFile.write('  0\n')
# number of states used in optimization
    fInFile.write('1\n')
# state specifications (J, nJ, P)
    fInFile.write('  0  1  +1   0.00000   0.00000\n')

# number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('23\n')
# number of valence nucleons
    fInFile.write('7\n')
# Experimental ground state energy of the nucleus
    fInFile.write('-37.141572   0.09000\n')
# min max delta j
    fInFile.write(' 0.5, 4.5, 1.0,\n')
# parity
    fInFile.write('  0\n')
# number of states used in optimization
    fInFile.write('1\n')
# state specifications (J, nJ, P)
    fInFile.write('  1/2  1  +1   0.00000   0.00000\n')

# number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('24\n')
# number of valence nucleons
    fInFile.write('8\n')
# Experimental ground state energy of the nucleus
    fInFile.write('-41.333144   0.10992\n')
# min max delta j
    fInFile.write(' 0.0, 4.0, 1.0,\n')
# parity
    fInFile.write('  0\n')
# number of states used in optimization
    fInFile.write('1\n')
# state specifications (J, nJ, P)
    fInFile.write('  0  1  +1   0.00000   0.00000\n')

# number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('25\n')
# number of valence nucleons
    fInFile.write('9\n')
# Experimental ground state energy of the nucleus
    fInFile.write('-40.55713   0.11093\n')
# min max delta j
    fInFile.write(' 0.5, 4.5, 1.0,\n')
# parity
    fInFile.write('  0\n')
# number of states used in optimization
    fInFile.write('1\n')
# state specifications (J, nJ, P)
    fInFile.write('  3/2  1  +1   0.00000   0.00000\n')

# number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('26\n')
# number of valence nucleons
    fInFile.write('10\n')
# Experimental ground state energy of the nucleus
    fInFile.write(' -41.24314    0.15551\n')
# min max delta j
    fInFile.write(' 0.0, 4.0, 1.0,\n')
# parity
    fInFile.write('  0\n')
# number of states used in optimization
    fInFile.write('1\n')
# state specifications (J, nJ, P)
    fInFile.write('  0  1  +1   0.00000   0.00000\n')

# code to create an input file
#CreateInFile('C:/PythonScripts/OxBashScripts/OptInput.in')


class BashOpt:
    'Class for shell model hamiltonian optimization problems'
    def __init__(self, sInPath, sOutPath, sErrorPath='', initialize=True,
                 fThError=0.10, sOBDir='c:\\oxbash', serial=False, 
                 sInitDir=''):
        self.sInitDir = sInitDir
        self.nDOF = 1
        self.sOBDir = sOBDir
        self.fLastChi = 0
        self.init = False
# flag that says to use the groundstate energy in optimization
        self.useGS = 1
# flag that tracks the residue and energy levels over the iterations
        self.sMethod = ''
        self.track = 1
        self.sInPath = sInPath
        self.sOutPath = sOutPath
# copy oxbash-dir.dat to the output folder
        import shutil
        import os
        if self.sInitDir != '':
            print 'Copying the initializaton from:"', self.sInitDir,'"...'
            initialize = False
            try:
                shutil.copytree(self.sInitDir, self.sOutPath)
                print 'Copy complete!'
            except WindowsError:
                print 'Initialization directory destination exists. Using existing directory.'
        if not os.path.isdir(self.sOutPath):
            os.makedirs(self.sOutPath)
        shutil.copyfile(self.sOBDir + '\\oxbash-dir.dat', self.sOutPath +
                        '\\oxbash-dir.dat')
        self.EExp = []
        #tbtd labels will be put here if necessary
        self.npaErrors = []
        self.serial = serial
        if not self.serial and os.path.isfile(self.sOutPath + '\\MPINames.dat'):
            os.remove(self.sOutPath + '\\MPINames.dat')
# get info from input file and initializes the nuclei objects
        self.GetIn(initialize)
        bRun=True
#       call runsm and write status if things are to be run in parallel 
        if not self.serial:
#            find last slash charachter
            i = -1
            for nIdx, char in enumerate(self.sInPath):
                if char == '\\':
                    i = nIdx
            spath = self.sInPath[:i]
            os.chdir(self.sOutPath)
            if initialize:
                from subprocess import call
                call('mpiexec -n 10 python ' + spath +'\\OxbashMPI.py')
                for nuc in self.mloNuclei:
                    nuc.runSM()
                    nuc.writeStatus()
                bRun=False
# set a value for the theory error in the determination of the energy levels
        self.fThError = fThError
        if self.useGS == 1:
            import math
            nIdx = 0
            for nucleus in self.mloNuclei:
                for state in nucleus.mllspec:
                    self.EExp[nIdx] += nucleus.fGSE
                    self.npaErrors[nIdx] = math.sqrt(self.npaErrors[nIdx]**2 +
                                                     nucleus.fError**2)
                    nIdx += 1
            self.dFitInfo = {'type': None, 'tolerance': None,
                             'iterations': None, 'iteration max': None,
                             'duration': None}
        import os
        if not os.path.exists(self.sOutPath+'\\'+'tracking'):
            os.makedirs(self.sOutPath+'\\'+'tracking')
        if initialize:
            self.writeLevs(self.sOutPath+'\\'+'tracking\\')
        self.llLastMESpec = [[], []]
        if initialize:
            self.obj(self.mloNuclei[0].getME(), bRun=bRun)
        self.init = True
#            self.EExp.shape = [self.EExp.size, 1]

# Write the initial state of the fit
    def writeLevs(self, path):
        npaETh = []
        lLS = []
        lAZ = []
        for nucleus in self.mloNuclei:
            import numpy
            tempth = list(nucleus.getEnNu())
            lLS.extend(nucleus.mllspec)
            lAZ.extend([nucleus.nAZ]*len(nucleus.mllspec))
            if npaETh != []:
                npaETh = numpy.append(npaETh, tempth, axis=0)
            else:
                npaETh = list(tempth)

        fOut = open(path + 'Levels.dat', 'w')
        sHFormat = '{:10}{:5}{:5}{:5}{:>10}{:>10}\n'
        fOut.write(sHFormat.format('[A,Z]', 'J', 'nJ', 'P', 'Eexp', 'Einit'))
        sFormat = '{:10}{:5}{:5}{:5}{:10.3f}{:10.3f}'
        sFormatAlt = '{:10}{:5}{:5}{:5}{:10.3f}{:>10}'
# print lAZ, lLS, npaEExp.shape, npaETh.shape
        for AZ, LS, Exp, Th in zip(lAZ, lLS, self.EExp, npaETh):
            try:
                fOut.write(sFormat.format(AZ, LS[0], LS[1], LS[2], float(Exp),
                                          float(Th)) + '\n')
            except:
                fOut.write(sFormatAlt.format(AZ, LS[0], LS[1], LS[2],
                                             float(Exp), Th) + '\n')
        fOut.close()

# Update the levels with their final values
    def updateLevs(self, path):
        npaETh = []
        lLS = []
        lAZ = []
        import numpy
        for nucleus in self.mloNuclei:
            tempth = list(nucleus.getEnNu())
            lLS.extend(nucleus.mllspec)
            lAZ.extend([nucleus.nAZ]*len(nucleus.mllspec))
            if npaETh != []:
                npaETh = numpy.append(npaETh, tempth, axis=0)
            else:
                npaETh = list(tempth)
        import os
        if not os.path.isfile(path + 'Levels_.dat'):
            os.rename(path + 'Levels.dat', path + 'Levels_.dat')
        fIn = open(path + 'Levels_.dat', 'r')
        fOut = open(path + 'Levels.dat', 'w')
        sHFormat = '{:10}{:5}{:5}{:5}{:>10}{:>10}{:>10}{:>10}\n'
        temp = fIn.readline().strip().split()
        temp.append('EFinal')
        temp.append('Final-Exp')
        fOut.write(sHFormat.format(*temp))
        sFormat1 = '{:5}{:5}{:5}{:5}{:5}{:>10}{:>10}'
        sFormat2 = '{:>10}'
        sFormatNew = '{:>10.3f}{:>10.3f}'
        sFormatAlt = '{:>10}{:>10}'
        nIdx = 0
        for line in fIn:
            if nIdx > -1:
                temp = line.strip().split()
                sNewLine = sFormat1.format(temp[0], temp[1], temp[2], temp[3],
                                           temp[4], temp[5], temp[6])
                for nJIdx in range(len(temp)-7):
                    sNewLine += sFormat2.format(temp[nJIdx + 7])
                try:
                    sNewLine += sFormatNew.format(float(npaETh[nIdx]),
                                                  float(npaETh[nIdx]) -
                                                  float(temp[5]))
                except:
                    sNewLine += sFormatAlt.format(npaETh[nIdx], 'N/A')
                fOut.write(sNewLine + '\n')
            nIdx += 1
        fIn.close()
        fOut.close()
        os.remove(path + 'Levels_.dat')

    def IterativeLSq(self, sMethod='single', fTolin=10**-3, nMaxIter=100,
                     bMix=False, methodArg=[]):
        '''
            start timing the optimization
        '''
        import time
        start = time.clock()
        self.sMethod = sMethod
#        fResLast=0
# store the original matrix element specificaton so it can be restored
# after it is altered.
        llOriginalMESpec = list(self.mloNuclei[0].llMESpec)
        fResNew = 100
        lRes = [100, 1000, 10000, 100000]
        nIter = 0
        fTol = fTolin #/ float(len(self.mloNuclei[0].mllspec))
        lME = [[], [], [], []]
        if sMethod == 'mono':
            npaMonoLabel = self.mloNuclei[0].getMonoLabel()
            for nucleus in self.mloNuclei:
                temp = nucleus.countOBME()
                if temp % 2 != 0:
                    print 'Warning: Odd number of OBMEs!'
                    raw_input('Press enter to continue.')
                llnOBMESpec = []
                for nIdx in range(temp/2):
                    llnOBMESpec.append([nIdx + 1, nIdx + temp / 2 + 1])
                nucleus.llMESpec = [llnOBMESpec, npaMonoLabel]
                nucleus.setmanBody([1, 2])
        while abs(lRes[0] - lRes[1]) > fTol and nIter < nMaxIter:
            for nIdx in range(3):
                lRes[len(lRes) - (nIdx+1)] = lRes[len(lRes) - (nIdx + 2)]
                lME[len(lRes) - (nIdx + 1)] = lME[len(lRes) - (nIdx + 2)]
            if sMethod == 'single':
                temp = self.singleParticleLeastSq()
                npaGuess = temp[0]
            elif sMethod == 'mono':
                temp = self.monopoleLeastSq()
                npaGuess = temp[0]
            elif sMethod == 'smono':
                temp = self.sMono(methodArg)
                npaGuess = temp[0]
                npaShortMono = temp[4]
            elif sMethod == 'csm':
                temp = self.compositeSingleMono(nMaxIter, fTolin, methodArg)
                npaGuess = temp[1]
                for nucleus in self.mloNuclei:
                    nucleus.llMESpec = temp[2]
            elif sMethod == 'TBTD':
                temp = self.TBTDLeastSq(methodArg)
                npaGuess = temp[0]
                npaGuess = npaGuess.flatten()
            else:
                print 'Error: sMethod=', sMethod, 'is not a valid argument.'
                break
            print "The guess is", npaGuess
            fResNew = self.obj(npaGuess)
#            raw_input('nIter='+str(nIter)+' press enter')
            if type(fResNew) != float:
                print ('Error: instance of ''Bashopt'' obj method ' +
                       'returning invalid type: '), type(fResNew)
            if sMethod == 'smono':
                self.sMonoIterationReport(temp, npaShortMono)
            lRes[0] = fResNew
            lME[0] = npaGuess
            nIter += 1
# restore the original MESpec
#            for nucleus in self.mloNuclei:
#                nucleus.llMESpec = list(llOriginalMESpec)
            print lME
            print fResNew
        if abs(lRes[0] - lRes[1]) < fTol:
            print 'Completed Successfully'
        else:
            print 'Iteration max reached: ', nIter
        end = time.clock()
        self.updateLevs(self.sOutPath + '\\' + 'tracking\\')
        self.dFitInfo['type'] = self.sMethod
        self.dFitInfo['tolerance'] = fTol
        self.dFitInfo['iterations'] = nIter
        self.dFitInfo['iteration max'] = nMaxIter
        self.dFitInfo['duration'] = end - start
        return fResNew, npaGuess, self.mloNuclei[0].llMESpec

    def performOptimization(self, sMethod='Nelder-Mead', llMonoBase=[],
                            dOptions=None):
        '''
            Use llMonoBase to make the list of matrix elements that will be
            optimized.
        '''
        import numpy
        if len(llMonoBase) != 0:
            lnpaShortMonoLab, lnpaMonoJLab = \
                self.constructMonoLists(llMonoBase)
# use the occupations to determine relevant degrees of freedom if the matrix
# elemnts are not specified
        npaMono = []
        npaSPOcc = []
        npaETh = []
        for nucleus in self.mloNuclei:
            temp = nucleus.calcMonoOcc()
            if npaMono != []:
                npaMono = numpy.append(npaMono, numpy.array(temp), axis=0)
            else:
                npaMono = numpy.array(temp)
        tempth = nucleus.getEnNu()
        if npaETh != []:
            npaETh = numpy.append(npaETh, tempth, axis=0)
        else:
            npaETh = tempth
        temp = nucleus.getReducedOcc()
        if npaSPOcc != []:
            npaSPOcc = numpy.append(npaSPOcc, temp, axis=0)
        else:
            npaSPOcc = numpy.array(temp)
        import sys
        sys.path.append('C:\PythonScripts\generalmath')
        import MatManip

        SPRMList = MatManip.getZeroCols(npaSPOcc)
        lNewSPList = np.array(range(self.mloNuclei[0].countOBME()))
        lNewSPList = MatManip.rmSlice(SPRMList, lNewSPList, 0)
        MonoRMList = MatManip.getZeroCols(npaMono)
        origTBME = numpy.array(self.mloNuclei[0].llMESpec[1])
        if SPRMList != []:
            npaSPOcc = MatManip.rmSlice(SPRMList, npaSPOcc, 1)
            npaMono = MatManip.rmSlice(MonoRMList, npaMono, 1)
            if len(llMonoBase) == 0:
                newTBME = numpy.array(MatManip.rmSlice(MonoRMList, origTBME,
                                                       0))
            else:
                temp = self.makeLong(lnpaShortMonoLab, lnpaMonoJLab,
                                     np.zeros(len(lnpaShortMonoLab)))
                newTBME = temp[0]
        for nucleus in self.mloNuclei:
            nucleus.llMESpec[0] = lNewSPList
            nucleus.llMESpec[1] = newTBME
            nucleus.setMEnum()
        from scipy.optimize import minimize
        if len(llMonoBase) == 0:
            npaGuess = self.mloNuclei[0].getME()
            res = minimize(self.obj, npaGuess, method=sMethod,
                           options=dOptions)
        else:
            npaGuess = self.mloNuclei[0].getOBME()
            npaGuess = np.append(npaGuess, np.zeros(len(llMonoBase)))
            newObj = lambda x: self.objSMonoNonLinear(x, lnpaShortMonoLab,
                                                      lnpaMonoJLab)
            res = minimize(newObj, npaGuess, method=sMethod,
                           options=dOptions)
        print res
        self.updateLevs(self.sOutPath+'\\'+'tracking\\')

#   non linear objective function for fitting the summed monopole matrix
#   elements. Take the npa uess composed of the single particle energies
#   followed by the monopole matrix element differences. Then use the short
#   mono lab and the jlab to make the long list of differneces. Then use that
#    to make the guess for the original objective function
    def objSMonoNonLinear(self, npaGuess, npaShortMonoLab, npaJLab):
        shortdiff = npaGuess[-len(npaShortMonoLab):]
        npaLongMonoLab, longdiff = self.makeLong(npaShortMonoLab,
                                                 npaJLab, shortdiff)
        npaOldTBME = self.mloNuclei[0].getTBME()
        longdiff.shape = npaOldTBME.shape
        npaNewTBME = npaOldTBME + longdiff
        import copy
        npaNewGuess = np.append(copy.copy(npaGuess[:-len(npaShortMonoLab)]),
                                npaNewTBME)
        return self.obj(npaNewGuess)

    def GetIn(self, initialize=True):
        import OxbashNuclei
        fIn = open(self.sInPath, 'r')
        sTempL = fIn.readline()
        llMESpec = [[], []]
        for nIdx in range(int(sTempL)):
            if nIdx == 0:
                anBody = []
                sTempL = fIn.readline()
                sTempL = sTempL.split()
                if sTempL[0] == 'OBME':
                    anBody.append(1)
                if len(sTempL[1:]) != 0:
                    llMESpec[0].append(sTempL[1:])
                elif sTempL[0] == 'TBME':
                    anBody.append(2)
                if len(sTempL[1:]) != 0:
                    llMESpec[1].append(sTempL[1:])
                elif sTempL[0] == 'MONO':
                    anBody.append(3)
                else:
                    print ("Error: invalid specification of matrix element" +
                           "category")
            else:
                sTempL = fIn.readline().strip('\n')
                llMESpec.append(sTempL)
            self.lsShared = []
            for nIdx in range(3):
                self.lsShared.append(fIn.readline().strip('\n'))
            self.sForm = fIn.readline().strip('\n')
            self.bExtrap = bool(eval(fIn.readline().strip('\n')))
# check if the a number was passed if it was read the data from this file if not
# assume it is a file and call a new script to solve the problem.
            self.mloNuclei = []
            testline = fIn.readline().strip('\n')
            try:
                self.mnNuclei = int(testline)
                for nIdx in range(self.mnNuclei):
                    '''
                        read in the variables to pass to the nucleus constructor
                        constructor signature:
                        __init__(self,  nZ, nA, fGSE, fError, sPath, sMMDR, sPar,
                                 lsShared,llMESpec,useGS,sForm,initialize=True,
                                 bExtrap=True)
                    '''
                    nZ = int(fIn.readline())
                    nA = int(fIn.readline())
                    nVal = int(fIn.readline())
                    templine = fIn.readline().strip().split()
                    fGSE = float(templine[0])
                    fError = float(templine[1])
                    sMMDR = fIn.readline().strip('\n')
                    sPar = fIn.readline().strip('\n')
                    temp = int(fIn.readline())
                    llStateSpec = []
                    for iii in range(temp):
                        line = fIn.readline().strip('\n')
                        line = line.split()
                        llStateSpec.append(line[:3])
                        self.EExp.append(float(line[3]))
                        self.npaErrors.append(float(line[4]))
                        tempnuc = OxbashNuclei.nucleus(nZ, nA, nVal, fGSE,
                                                       fError, self.sOutPath,
                                                       sMMDR, sPar,
                                                       self.lsShared, llMESpec,
                                                       self.useGS, self.sForm,
                                                       initialize,
                                                       self.bExtrap,
                                                       llStateSpec, serial=self.serial)
                    tempnuc.setmanBody(anBody)
                    self.mloNuclei.append(tempnuc)
                import numpy as np
                self.npaErrors = np.array(self.npaErrors, dtype=float)
                self.EExp = np.array(self.EExp, dtype=float)
            except ValueError:
                self.readDAI(testline, initialize)
            fIn.close()

#   reads BAB's SD shell fitting information
    def readDAI(self, sDAIpath, initialize):
        import numpy as np
        import OxbashNuclei
        from Discrete import GCDofList
        llMESpec = [[], []]
        anBody = 3
        import os.path
        if os.path.isfile(sDAIpath):
            fDAI = open(sDAIpath, 'r')
        else:
            print ('Error: In OxbashOptObj.py -> readDAI "' + sDAIpath +
                   '" is not a valid path.')
            return 'error'
        lnLastAI = []
        npa2JnJ = []
        sMMDRform = ' {:>2.1f}, {:>2.1f}, {:>2.1f}\n'
        mnNuclei = 0
        for line in fDAI:
#            print line
            line = line.strip().split()
            lnCurentAI = [int(line[0]), int(line[1])]
            if not np.all(lnLastAI == lnCurentAI):
                mnNuclei += 1
# if this is not the first nucleus write the one you just read in
                if lnLastAI != []:
                    if npa2JnJ.shape[0] > 1:
                        fMin = np.amin(npa2JnJ[:, 0])/2.
                        fMax = np.amax(npa2JnJ[:, 0])/2.
                        fDeltaJ = GCDofList(np.diff(npa2JnJ[:, 0]))/2.
                    else:
                        fMin = npa2JnJ[0, 0]/2.
                        fMax = fMin
                        fDeltaJ = 0.0
                    sMMDR = sMMDRform.format(fMin, fMax, fDeltaJ)
                    llStateSpec = []
                    sdefpar = '+1'
                    for nIdx in range(npa2JnJ.shape[0]):
                        j = int(npa2JnJ[nIdx, 0])
                        if j % 2 != 0:
                            j = str(j)+'/2'
                        else:
                            j = str(j/2)
                        nj= str(npa2JnJ[nIdx, 1])
                        llStateSpec.append([j, nj, sdefpar])
                    tempnuc = OxbashNuclei.nucleus(nZ, nA, nVal, fGSE, fError,
                                                   self.sOutPath, sMMDR, sPar,
                                                   self.lsShared, llMESpec,
                                                   self.useGS, self.sForm,
                                                   initialize, self.bExtrap,
                                                   llStateSpec, serial=self.serial)
                    tempnuc.setmanBody(anBody)
#                    print tempnuc.nAZ
                    self.mloNuclei.append(tempnuc)
# get as much info as you can about the new nucleus from first line
                nA = int(line[0])
                nVal= nA - 16
                fGSE = float(line[4]) # GSE is first energy given for a nucleus
                fError = float(line[5]) # see above
                sPar = '0'
                nZ = (nVal - int(line[1]))/2 + 8
                npa2JnJ = np.array([int(line[2]), int(line[3])])
                npa2JnJ.shape = [1, 2]
                if self.useGS == 1:
                    self.EExp.append(0.0) # first energy is the ground state
                    self.npaErrors.append(0.0)
                else:
                    self.EExp.append(fGSE) # first energy is the ground state
                    self.npaErrors.append(fError)                    
            else:
                temp = np.array([int(line[2]), int(line[3])])
                temp.shape = [1, 2]
                npa2JnJ = np.append(npa2JnJ, temp, 0)
                self.EExp.append(float(line[4])) # if it is not the first line for this nucleus this is the excitation
                self.npaErrors.append(float(line[5]))
            lnLastAI = [elem for elem in  lnCurentAI]
# do everything again to get the last nucleus
        if npa2JnJ.shape[0] > 1:
            fMin = np.amin(npa2JnJ[:, 0])/2.
            fMax = np.amax(npa2JnJ[:, 0])/2.
            fDeltaJ = GCDofList(np.diff(npa2JnJ[:, 0]))/2.
        else:
            fMin = npa2JnJ[0, 0]/2.
            fMax = fMin
            fDeltaJ = 0.0
        sMMDR = sMMDRform.format(fMin, fMax, fDeltaJ)
        llStateSpec = []
        sdefpar = '+1'
        for nIdx in range(npa2JnJ.shape[0]):
            j = int(npa2JnJ[nIdx, 0])
            if j % 2 != 0:
                j = str(j)+'/2'
            else:
                j = str(j/2)
            nj = str(npa2JnJ[nIdx, 1])
            llStateSpec.append([j, nj, sdefpar])
        tempnuc = OxbashNuclei.nucleus(nZ, nA, nVal, fGSE, fError,
                                       self.sOutPath, sMMDR, sPar,
                                       self.lsShared, llMESpec,
                                       self.useGS, self.sForm,
                                       initialize, self.bExtrap,
                                       llStateSpec, serial=self.serial)
        tempnuc.setmanBody(anBody)
        self.mloNuclei.append(tempnuc)
        self.mnNuclei = mnNuclei
        fDAI.close()


# The objective function takes the array of matrix elements
    def obj(self, npaME, bRun=True):
        res = []
        numFound = 0
        numDesired = 0
#       run the SM calculations in parallel if it is desired
        if (not self.serial) and bRun:
            import os
#            find last slash charachter
            i = -1
            for nIdx, char in enumerate(self.sInPath):
                if char == '\\':
                    i = nIdx
            spath = self.sInPath[:i]
            os.chdir(self.sOutPath)
            from subprocess import call
            call('mpiexec -n 10 python ' + spath +'\\OxbashMPI.py')
        for oNuc in self.mloNuclei:
            '''
                write ME to file
            '''
            if self.init:
#                try:
                oNuc.takeME(npaME)
#                except:
#                    print 'Error: The Matrix element list is invalid.'
#                    print npaME
#                    print 'nAZ =', oNuc.nAZ
#                    return 'Error: The Matrix element list is invalid.'
            if self.sForm == 'pn':
                import os
                sLevName = oNuc.getLevName()
                temp = (str(oNuc.sPath) + '\\' + str(oNuc.sName) + '\\' +
                        sLevName[: -5] + '0' + '.int')
                if os.path.isfile(temp + '_'):
                    os.remove(temp + '_')
                    os.rename(temp, temp + '_')
        if (not self.serial) and bRun:
            import os
#            find last slash charachter
            i = -1
            for nIdx, char in enumerate(self.sInPath):
                if char == '\\':
                    i = nIdx
            spath = self.sInPath[:i]
            os.chdir(self.sOutPath)
            from subprocess import call
            call('mpiexec -n 10 python ' + spath +'\\OxbashMPI.py')

#           run Shell model calc if it was run in parallel this will just
#           compile the tables of results.
        if bRun or not self.serial:
            for oNuc in self.mloNuclei:
                oNuc.runSM()
# get the energy difference for the relevant levels
                tempEth = oNuc.getEnNu()
                numFound += len(tempEth)
                numDesired += len(oNuc.mllspec)
                if len(tempEth) - len(oNuc.mllspec) != 0:
                    print oNuc.nAZ
                    print 'Looked for', len(oNuc.mllspec), 'levels. Foumd', len(tempEth)
                    print 'Discrepency is:', numFound - numDesired
                    raw_input('Press Enter...')
                res.extend(tempEth)
        npaETh = [e for e in res]
        fChiSq, fRes, fDelta = self.calcChiDOF(npaETh, bWeight=False)
        if self.track == 1:
            self.OptStatus(fRes, fChiSq, fDelta, npaME)
        self.fLastChi = fChiSq
        return fRes

#   calculate the chi square statistic divided by the number DOF for the fit if
#   bWeight is true return a vector of fitting weights (along with energies)
#   instead. Make sure the variable self.nDOF is sef correctly before calling.
    def calcChiDOF(self, npaETh, bWeight=False):
        import numpy as np
        lnThRm = []
        for i, elem in enumerate(npaETh):
            if elem == 'N/A':
                lnThRm.append(i)
        npaNewETh = np.array([E for i, E in enumerate(npaETh) if i not in lnThRm], dtype = float)
        npaNewETh.shape = [npaNewETh.size, 1]
        npaNewEExp = np.array([E for i, E in enumerate(self.EExp) if i not in lnThRm], dtype = float)
        npaNewEExp.shape = [npaNewEExp.size, 1]
        lfGSError = []
        for nucleus in self.mloNuclei:
            for num in range(len(nucleus.mllspec)):
                lfGSError.append(nucleus.fError**2)
        npaGSError = np.array([E**2 for i, E in enumerate(lfGSError) if i not in lnThRm], dtype = float)
        npaNewErr = np.array([E**2 for i, E in enumerate(self.npaErrors) if i not in lnThRm], dtype = float)
        if self.useGS == 1:
            npaNewErr = npaNewErr + npaGSError # add ground state error to the level error
        npaErrors = np.ones(*(npaNewErr.shape))*self.fThError**2 + npaNewErr # add the arbitrary theory error to the previous error contributions
        npaWeights = np.zeros(npaErrors.shape)
        for nIdx, elem in enumerate(npaErrors):
            npaWeights[nIdx] = 1.0 / (elem)
        npaWeights.shape = [npaWeights.size, 1]
        if bWeight:
            npaWeights.shape = [npaWeights.size]
            return npaWeights, npaNewETh, npaNewEExp
        else:
            npaRes = npaNewETh - npaNewEExp
            fChi = np.dot(np.transpose(np.power(npaRes, 2)), npaWeights)
            fChi = float(fChi[0])       
            temp = np.dot(np.transpose(npaRes), npaRes)
            temp = np.sqrt(np.sum(temp)/npaRes.size)
            term1 = 2.*((2.*10**-5/float(self.nDOF))*np.sum(np.divide(npaRes, np.power(npaErrors,2))))**2
#            term2 = ((2.*10**-5/float(self.nDOF))*np.sum(np.divide(npaRes, np.power(npaErrors,3))))**2
            ferr = np.sqrt(term1)
            return fChi / float(self.nDOF), float(temp), ferr

    def OptStatus(self, res, fChiSq, fDelta,npaME):
        sFormat = '{:14.10f}'
        fOut = open(self.sOutPath + '\\tracking\\res.dat', 'a')
        fOut.write(sFormat.format(res) + '\n')
        fOut.close()
        if self.init:
            fOut = open(self.sOutPath + '\\tracking\\chisq.dat', 'a')
            fOut.write(sFormat.format(fChiSq) + sFormat.format(fDelta) + '\n')
            fOut.close()
        fOut = open(self.sOutPath+'\\tracking\\ME.dat', 'a')
        string = ''
        sFormMELab = '{:>13}\t'
        sFormME = '{:13.4f}\t'
        import numpy as np
        bSPEequal = np.all(np.array(self.llLastMESpec[0]) ==
                           np.array(self.mloNuclei[0].llMESpec[0]))
        bTBMEequal = np.all(np.array(self.llLastMESpec[1]) ==
                            np.array(self.mloNuclei[0].llMESpec[1]))
        if not (bSPEequal and bTBMEequal):
            for elem in self.mloNuclei[0].llMESpec:
                for lab in elem:
                    string += sFormMELab.format(lab)
            fOut.write(string + '\n')
            self.llLastMESpec = [elem for elem in self.mloNuclei[0].llMESpec]
        string = ''
        for elem in npaME:
            string = string + sFormME.format(float(elem))
        fOut.write(string + '\n')
        fOut.close()
        fOut = open(self.sOutPath + '\\tracking\\AllME.dat', 'a')
        npaAllME = self.mloNuclei[0].getME(bAll=True)
        sFormME = '{:10.4f}\t'
        string = ''
        for elem in npaAllME:
            string = string + sFormME.format(elem)
        fOut.write(string + '\n')
        fOut.close()
        for oNucleus in self.mloNuclei:
            oNucleus.writeStatus()

# returns the single particle energy lest square solution to the energy
    def singleParticleLeastSq(self):
        import numpy
        a = []
        npaETh = []
        for nucleus in self.mloNuclei:
            sLevName = nucleus.getLevName()
            npaOcc = nucleus.getOcc(sLevName)
            if a != []:
                a = numpy.append(a, npaOcc, axis=0)
            else:
                a = numpy.array(npaOcc)
            tempth = nucleus.getEnNu()
            if npaETh != []:
                npaETh = numpy.append(npaETh, tempth, axis=0)
            else:
                npaETh = tempth
        import sys
        sys.path.append('C:\PythonScripts\generalmath')
        import MatManip
# get the zero cols of the matrix
        rmList = MatManip.getZeroCols(a)
        if rmList != []:
            a = MatManip.rmSlice(rmList, a, 1)
        lLabSpec = self.constructLabels(rmList)
        lsOBLab = []
        npaTBLab = []
        for elem in lLabSpec:
            if len(elem) == 1:
                lsOBLab.append(elem[0])
            else:
                if npaTBLab != []:
                    temp = numpy.array(elem, dtype=int)
                    temp.shape = [1, temp.size]
                    npaTBLab = numpy.append(npaTBLab, temp, axis=0)
                else:
                    npaTBLab = numpy.array(elem, dtype=int)
                    npaTBLab.shape = [1, npaTBLab.size]
        for nucleus in self.mloNuclei:
            nucleus.llMESpec[0] = lsOBLab
            nucleus.setMEnum()
        obme = self.mloNuclei[0].getOBME()
        obme.shape = [obme.size, 1]
        npaETh.shape = [npaETh.size, 1]
        print obme
        print self.mloNuclei[0].llMESpec
        Hexpect = npaETh - numpy.dot(a, obme)
        target = self.EExp - Hexpect
#     npaWeights=numpy.zeros([self.EExp.size,self.EExp.size])
#     for nIdx, elem in enumerate(self.npaErrors):
#       print elem
#       npaWeights[nIdx,nIdx]=1.0/(elem**2)

# adjust for weighted least squares
#     a=numpy.dot(npaWeights, a)
#     target=numpy.dot(npaWeights, target)
# end adjust for weighted least squares
        ans = numpy.linalg.lstsq(a, target)
        ans = ans[0]
        return ans, a, target, obme, self.mloNuclei[0].llMESpec[0]

# returns list of non super diagonal matrix elements
#     used to test convergence of iterative scheme when there is little data
    def getNonSuDi(self):
        ans = [row for row in self.mloNuclei[0].llMESpec[1]
               if not(row[0] == row[1] and row[0] == row[2] and
                      row[0] == row[3])]
        nonsudi = []
        for elem1 in ans:
            for iii, elem2 in zip(range(len(self.mloNuclei[0].llMESpec[1])),
                                  self.mloNuclei[0].llMESpec[1]):
                if all(elem1 == elem2):
                    nonsudi.append(iii)
        return nonsudi

# returns the single particle energy + monopole lest square solution to the
# energy
    def monopoleLeastSq(self):
        import numpy
        npaMono = []
        npaSPOcc = []
        npaETh = []
        for nucleus in self.mloNuclei:
            temp = nucleus.calcMonoOcc()
            if npaMono != []:
                npaMono = numpy.append(npaMono, numpy.array(temp), axis=0)
            else:
                npaMono = numpy.array(temp)
            tempth = nucleus.getEnNu()
            if npaETh != []:
                npaETh = numpy.append(npaETh, tempth, axis=0)
            else:
                npaETh = tempth
        temp = nucleus.getReducedOcc()
        if npaSPOcc != []:
            npaSPOcc = numpy.append(npaSPOcc, temp, axis=0)
        else:
            npaSPOcc = numpy.array(temp)
        import sys
        sys.path.append('C:\PythonScripts\generalmath')
        import MatManip

        SPRMList = MatManip.getZeroCols(npaSPOcc)
        MonoRMList = MatManip.getZeroCols(npaMono)
#     to test the iterative scheme for small data sets
        MonoRMList.extend(self.getNonSuDi())
        MonoRMList = list(set(MonoRMList))
#    cut down even more
        temp = range(max(MonoRMList))
        new = [iii for iii in temp if iii not in MonoRMList]
        MonoRMList.extend(new[-3:])
        MonoRMList = list(set(MonoRMList))
        if SPRMList != []:
            npaSPOcc = MatManip.rmSlice(SPRMList, npaSPOcc, 1)
            npaMono = MatManip.rmSlice(MonoRMList, npaMono, 1)
            origTBME = numpy.array(self.mloNuclei[0].llMESpec[1])
            newTBME = numpy.array(MatManip.rmSlice(MonoRMList, origTBME, 0))
        for nucleus in self.mloNuclei:
            nucleus.llMESpec[0] = []
            nucleus.llMESpec[1] = newTBME
            nucleus.setMEnum()
        npaME = self.mloNuclei[0].getME()
#     a=numpy.append(npaSPOcc,npaMono,1)
        a = npaMono
#     print a
        target = npaEExp - (npaETh - numpy.dot(a, npaME))
#     npaWeights=numpy.zeros([npaEExp.size,npaEExp.size])
#     for nIdx, elem in enumerate(self.npaErrors):
#       print elem
#       npaWeights[nIdx,nIdx]=1.0/(elem**2)
#     ans = numpy.linalg.lstsq(numpy.dot(npaWeights, a), numpy.dot(npaWeights,
#                                                                   target))
        ans = numpy.linalg.lstsq(a, target)
        ans = ans[0]

        return ans, a, target, npaME
    def TBTDinit(self):
        '''
            Get the two body transition densities and colect their labels for
            all nuclei
        '''
        import numpy as np
        import MatManip
        npaTBTDLabels = []
        npaTBTD = []
        for nucleus in self.mloNuclei:
    #          get and scale TBTD according to the interaction
            tempTBTD, tempLabel = nucleus.getTBTD()
            tempTBTD = np.multiply(tempTBTD, nucleus.fScale)
            if len(npaTBTDLabels) == 0 or \
            (np.all(np.equal(np.shape(npaTBTDLabels),np.shape(tempLabel))) and \
            np.all(np.equal(npaTBTDLabels, tempLabel))):
                if npaTBTDLabels == []:
                    npaTBTDLabels = np.array(tempLabel)
                if len(npaTBTD) == 0:
                    npaTBTD = tempTBTD
                else:
                    npaTBTD = np.append(npaTBTD, tempTBTD, axis=0)
            else:
                npaTBTDLabels, npaTBTD = \
                    MatManip.combinedLabeledColumns(npaTBTDLabels, npaTBTD,
                                                    tempLabel, tempTBTD)
        return np.array(npaTBTDLabels), npaTBTD
        
        
#    do a linear least square optimization of the tbme using the TBTDs
    def TBTDLeastSq(self, methodArg=[], bValid=True):
        if methodArg == []:
            nLC = 30
        else:
            nLC = methodArg[0]
        import numpy as np
        import MatManip
        npaTBTDLabels, npaTBTD = self.TBTDinit()
        npaETh = []
        for nucleus in self.mloNuclei:
            nucleus.llMESpec[0] = list(np.add(range(self.mloNuclei[0].nOBME),
                                       np.ones(self.mloNuclei[0].nOBME,
                                               dtype=int)))
            nucleus.llMESpec[1] = npaTBTDLabels
            tempth = nucleus.getEnNu()
            if len(npaETh) != 0:
                npaETh = np.append(npaETh, tempth, axis=0)
            else:
                npaETh = tempth
        npaETh = np.array(npaETh)
        npaETh.shape = [npaETh.size, 1]
        npaOcc = []
        nSought = 0
        nFound = 0
        for nucleus in self.mloNuclei:
            sLevName = nucleus.getLevName()
            npaTempOcc = nucleus.getOcc(sLevName)
            nTemp = len(nucleus.mllspec)
            nSought += nTemp
            nFound += npaTempOcc.shape[0]
            if nTemp != npaTempOcc.shape[0]:
                print 'Discrepency in sought and found Occupations in', nucleus.nAZ 
                print 'Sought:', nTemp, 'Found:', npaTempOcc.shape[0] 
            if len(npaOcc) != 0:
                npaOcc = np.append(npaOcc, npaTempOcc, axis=0)
            else:
                npaOcc = npaTempOcc
        npaME = self.mloNuclei[0].getOBME(bAll=True)
        npaME = np.append(npaME, self.mloNuclei[0].getTBME(bAll=False))
#        npaME = self.mloNuclei[0].getME(bAll=False)
        npaME.shape = [npaME.size, 1]
        a = np.append(npaOcc, npaTBTD, axis=1)
        nlRMList = MatManip.getZeroCols(a)
        newTBMEList = []
        if nlRMList != []:
            lnTBMERmList = [elem for elem in nlRMList if elem > 2]
            a= np.delete(a, lnTBMERmList, 1)
            npaME = np.delete(npaME, lnTBMERmList, 0)
            lnTBMERmList = np.subtract(lnTBMERmList, self.mloNuclei[0].nOBME)
            newTBMEList = np.delete(self.mloNuclei[0].llMESpec[1], lnTBMERmList,0)
#       check if llsq problem is valid
        if a.shape[0] <= nLC:
            print 'Warning system underdeterimined with:'
            print nLC, 'Linear combinationsto fit and only', a.shape[0], 'equations.'
            print 'Reducing nuber of linear combinations to', a.shape[0]
            nLC = a.shape[0]
        if len(newTBMEList) != 0:
            for nuc in self.mloNuclei:
                newTBMEList = np.array(newTBMEList)
                nuc.llMESpec[1] = newTBMEList
        lnThRm = []
        for i, elem in enumerate(npaETh):
            if np.all(elem == 'N/A'):
                lnThRm.append(i)
        a = np.delete(a, lnThRm, 0)
        self.nDOF = a.shape[0] - nLC - 1
        npaWeights, npaNewETh, npaNewEExp = self.calcChiDOF(npaETh,
                                                            bWeight=True)
        npaWeights = np.diag(npaWeights)
        aprime = np.dot(npaWeights, a)
        if len(methodArg) > 1 and methodArg[1] == 'fd':
            #fit differences
            target = npaNewEExp - npaNewETh
            bFitDif = True
        else:
            #fit absolute energy
            target = npaNewEExp
            bFitDif = False
        if bValid:
            num=np.linalg.norm(np.dot(a,npaME)-npaNewETh)
            fval = open(self.sOutPath+'\\val.dat', 'a')
            fval.write('{:10.4f}'.format(num))
            fval.close()
        bprime = np.dot(npaWeights, target)
        npaBGME = self.mloNuclei[0].getBGInt()
        ans = self.linCom(aprime, npaBGME, bprime, nLincoms=nLC, bFD =bFitDif)
        return ans, a, target, npaME

# lst square solution to matrix a, with initial guess xbg, target b, and an 
# epsilon threshold for determining well determined lincoms
    def linCom(self, a, xbg, b, epsilon=-0.1, nLincoms = 0, bFD=False):
        import numpy as np
        [u, s, vt] = np.linalg.svd(a)
#        use default value if no valid epsilon is given 
        if epsilon <= 0 and (nLincoms < 0 or nLincoms > s.size):
            print 'Object: oxbashopt. method: linCom.'
            print 'Using default number of linear combinations.'
            try:
                nKeep = int(.4*float(a.shape[0]))
                epsilon = s[nKeep + 1]
                print '40% of number of data points: keep', nKeep, 'variables'
            except:
                epsilon = s[-1]
            print 'cutoff is ', epsilon
        elif nLincoms == 0:
            print 'oxcbashopt.lincom.'
            print 'nLincoms = 0 returninging background interaction.'
            return xbg
        else:
            epsilon = s[nLincoms - 1]
#            epsilon= s[nLincoms]
        if a.shape[1] > a.shape[0]:
            print 'Warning: System of equations is underdetermined.'
            print a.shape[0], 'equations', a.shape[1], 'unknown quantities.'
        e = np.dot(np.transpose(u), b)
        d = np.zeros([a.shape[1], a.shape[0]])
        d[:np.amin(d.shape), :np.amin(d.shape)] = np.diag(1./s)
        yn = np.dot(d, e)
        ys = np.dot(vt, xbg)
        y = []
        numbelow = 0
        epsilon -= np.finfo(type(epsilon)).eps
        if bFD:
            y0 = self.mloNuclei[0].getME()
            y0 = np.dot(vt, y0)
        for i in range(ys.size):
            if bFD:
                if i < s.size and s[i] >= epsilon:
                    y.append(yn[i,0]+y0[i])
                    print yn[i,0], y0[i]
                else:
                    y.append(ys[i])
                    numbelow += 1
            else:
                if i < s.size and s[i] >= epsilon:
                    y.append(yn[i])                
                else:
                    y.append(ys[i])
                    numbelow += 1
        print numbelow, 'lincoms of ', s.size, 'are not well determined.'
        self.lincomreport(a.shape[0], s.size, numbelow, epsilon)
        y = np.array(y)
        y.shape = [vt.shape[0], 1]
        return np.dot(np.transpose(vt), y)

# write to file the number of linear combinations used and the cutoff that determines it
    def lincomreport(self, numdat, numvar, numignore, eps):
        import os
        sPath = self.sOutPath + '\\tracking\\lincom.dat'
        if not os.path.isfile(sPath):
            sFormHead = '{:>10}'*4+'\n'
            fOut = open(self.sOutPath + '\\tracking\\lincom.dat', 'w')
            fOut.write(sFormHead.format('Data', 'Variables', 'Ignored',
                                        'Cutoff'))
        else:
            fOut = open(self.sOutPath + '\\tracking\\lincom.dat', 'a')
        sForm = '{:>10d}{:>10d}{:>10d}{:>10.4f}\n'
        fOut.write(sForm.format(numdat, numvar, numignore, eps))
        fOut.close()
        
# quickly set manbody variable
    def initmanbody(self):
        for nucleus in self.mloNuclei:
            nucleus.setmanBody([1, 2])

# Get cols of poorly determinied ME
    def getBadCol(self, npaMat):
        import numpy as np
        import MatManip
        rmList = MatManip.getZeroCols(npaMat)
        npaMat2 = MatManip.rmSlice(rmList, npaMat, 1)
        npaME = self.mloNuclei[0].getME()
        rmList = MatManip.getRepeatSlices(npaMat2, 1)
        npaMat3 = MatManip.rmSlice(rmList, npaMat2, 1)
        npaME2 = MatManip.rmSlice(rmList, npaME, 0)
        psi = np.dot(np.transpose(npaMat3), npaMat3)
        npaErr = np.diagonal(np.linalg.inv(psi))
        lCond = [me / err for me, err in zip(npaME2, npaErr)]
        rmList = [nIdx for nIdx, cond in enumerate(lCond) if abs(cond) > 1]
        npaME3 = [npaME2[nIdx] for nIdx in rmList]
        rmList = [nIdx for nIdx, me in enumerate(npaME) if me in npaME3]
        return rmList

# calculate linear fit errors for the matrix elements
    def calcError(self):
        import numpy as np
        import MatManip
        if self.sMethod == 'mono':
            ans, a, target, npaME = self.monopoleLeastSq()
        elif self.sMethod == 'smono':
            ans, a, target, npaME, stuff1, stuff2 = self.monopoleLeastSq()
        elif self.sMethod == 'single':
            ans, a, target, npaME = self.singleParticleLeastSq()
        rmList = MatManip.getRepeatSlices(a, 1)
        ap = MatManip.rmSlice(rmList, a, 1)
        psi = np.dot(np.transpose(ap), ap)
        lME = MatManip.rmSlice(rmList, npaME, 0)
        ans = MatManip.rmSlice(rmList, ans, 0)
        return np.sqrt(np.diagonal(np.linalg.inv(psi))), lME, rmList, ans

# plot the residual and the matrix elements as a funtion of iteration number
    def plotResults(self, sMethod='single', bError=False):
        self.sMethod = sMethod
        self.initmanbody()
        lMEoriginal = list(self.mloNuclei[0].llMESpec)
        fResIn = open(self.sOutPath + '\\tracking\\res.dat', 'r')
        npaRes = []
        for line in fResIn:
            npaRes.append(line.strip())
        fResIn.close()
        import numpy as np
        npaRes = np.array(npaRes)
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(npaRes, label='Res')
        plt.plot()
        plt.title('RMS Energy Error by Iteration Number')
        plt.xlabel('Number of Iterations')
        plt.ylabel('RMS Energy Error (MeV)')
        plt.show()
        fMEIn = open(self.sOutPath + '\\tracking\\ME.dat', 'r')
        npaME = []
        for line in fMEIn:
            temp = line.strip().split()
            if len(npaME) != 0 and len(npaME[0]) == len(temp):
                npaME.append(temp)
            else:
                npaME = [temp]
                npaME = np.array(npaME)
        import MatManip
        import Discrete
        npaErr, lME, rmList, npaAns = self.calcError()
        print "Uncertaunty in ME", npaErr
        npaME = MatManip.rmSlice(rmList, npaME, 1)
#     rmShortList=[int(elem)-6 for elem in rmList if elem>5]
# construct labels
        lLabSpec = self.constructResultLabels(rmList)
        lsLabels = []
#     re initialize the list llmespec
        for nucleus in self.mloNuclei:
            nucleus.llMESpec = list(lMEoriginal)
        for elem in lLabSpec:
            if str(type(elem)) == '<type \'numpy.int32\'>':
                lsLabels.append('SPE: ' + str(elem))
            else:
                tempstr = '[ '
            for nIdx in elem:
                tempstr += str(nIdx) + ' '
                lsLabels.append('TBME: ' + tempstr + ']')
#     plot the ME
        lBest = Discrete.balFact(npaME.shape[1])
        fig, ax = plt.subplots(nrows=lBest[0], ncols=lBest[1], sharex=True,
                               sharey=False)
        from itertools import chain
        try:
            ax = list(chain.from_iterable(ax))
        except TypeError:
            print ax, 'is not iterable'
        from itertools import chain
        for nColIdx in range(npaME.shape[1]):
            ax[nColIdx].plot(npaME[:, nColIdx], label=lsLabels[nColIdx])
            if self.sMethod == 'single' or self.sMethod == 'mono':
                ax[nColIdx].plot([0, npaME.shape[0] - 1], [npaAns[nColIdx],
                                 npaAns[nColIdx]], ls='-.', color='k',
                                 label='Next Linear Fit')
            sFormat = "Final ME:{:3.2f}$\pm${:2.1e}"
            ax[nColIdx].plot([0, npaME.shape[0] - 1], [lME[nColIdx],
                             lME[nColIdx]], ls='--', color='r',
                             label=sFormat.format(float(lME[nColIdx]),
                                                  float(npaErr[nColIdx])))
            if bError:
                x = [0, int(npaME.shape[0]) - 1]
                y1 = [lME[nColIdx] + npaErr[nColIdx]]*2
            try:
                y1 = list(chain(*y1))
            except:
                ''
                print 'y1 not iterable'
            y2 = [lME[nColIdx] - npaErr[nColIdx]]*2
            try:
                y2 = list(chain(*y2))
            except:
                ''
                print 'y2 not iterable'
            ax[nColIdx].fill_between(x, y1, y2, alpha=0.5, label='Error Band')
#     ax[1].title='ME by Iteration Number'
#       ax[nColIdx].xlabel('Number of Iterations')
#       ax[nColIdx].ylabel('ME (MeV)')
            ax[nColIdx].legend(loc='best')
            plt.show()

# consrtuct a single list of labels for use in the fitting procedure
    def constructLabels(self, rmList):
        from itertools import chain
#        for nucleus in self.mloNuclei:
#            nucleus.llMESpec[0] = list(range(1, self.mloNuclei[0].countOBME() +
#                                             1))
        lLabSpec = list(chain.from_iterable(self.mloNuclei[0].llMESpec))
        temp = []
        for elem in lLabSpec:
            if type(elem).__name__ == 'int':
                temp.append([elem])
            else:
                temp.append(elem)
        lLabSpec = list(temp)
        lLabSpec = list([elem for nIdx,
                         elem in enumerate(lLabSpec) if nIdx not in rmList])
        return lLabSpec

# consrtuct a single list of labels for use in Plotting the results
# the spe are not reinitialized here
    def constructResultLabels(self, rmList):
        from itertools import chain
        lLabSpec = list(chain.from_iterable(self.mloNuclei[0].llMESpec))
        temp = []
        for elem in lLabSpec:
            if type(elem).__name__ == 'int':
                temp.append([elem])
            else:
                temp.append(elem)
        lLabSpec = list(temp)
        lLabSpec = list([elem for nIdx, elem in enumerate(lLabSpec)
                         if nIdx not in rmList])
        return lLabSpec

# consrtuct a single list of labels
    def constructML(self, rmList, labels):
        from itertools import chain
        for nucleus in self.mloNuclei:
            temp = range(1, self.mloNuclei[0].countOBME() + 1)
            nucleus.llMESpec[0] = list(temp)
        temp = range(1, self.mloNuclei[0].countOBME() + 1)
        lLabSpec = list(chain.from_iterable([temp, labels]))
        temp = []
        for elem in lLabSpec:
            if type(elem).__name__ == 'int':
                temp.append([elem])
            else:
                temp.append(elem)
        lLabSpec = temp
        lLabSpec = [elem for nIdx, elem in enumerate(lLabSpec)
                    if nIdx not in rmList]
        return lLabSpec

# add mean zero gaussian noise to the matrix elements with specified standard
# deviation
    def addMENoise(self, fStd):
        from numpy.random import randn
        npaME = self.mloNuclei[0].getME()
        npaME += randn(*npaME.shape)*fStd
        for nucleus in self.mloNuclei:
            nucleus.takeME(npaME)

#  Record how many iterations it takes to reach convergence for different
# convergence tolerances.
    def CheckConvergenceSensitivity(self, nIterMax, sMethod, display=True):
        ldInfo = []
        nPow = 1
        import shutil
        import copy
        while (self.dFitInfo['type'] == None or
                self.dFitInfo['iteration max'] != self.dFitInfo['iterations']):
            self.IterativeLSq(sMethod, 10**-nPow, nIterMax, bMix=False)
            ldInfo.append(copy.copy(self.dFitInfo))
# cleanup the directory in preparation for the next iteration (here be dragons)
# change directory so the output directory may be removed
        shutil.os.chdir('..')
        shutil.os.chdir('..')
        shutil.rmtree(self.sOutPath)
        if not shutil.os.path.exists(self.sOutPath):
            shutil.os.makedirs(self.sOutPath+'\\tracking')
# initialize the the optimization problem so that the process may be repeated
            self.GetIn(True)
            self.writeLevs(self.sOutPath + '\\tracking\\')
            nPow += 1
        if display:
            self.ConvergenceSensitivityReport(ldInfo)
        return ldInfo

#   plot the results of the convergence sensitivity test
    def ConvergenceSensitivityReport(self, ldInfo):
        import matplotlib.pyplot as plt
        lfTol = []
        lnIter = []
        lfDur = []
        for adict in ldInfo:
            lfTol.append(adict['tolerance'])
            lfDur.append(adict['duration'])
            lnIter.append(adict['iterations'])
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(lfTol, lnIter)
        ax[0].set_xlabel('Tolerance (MeV)')
        ax[0].set_ylabel('Iterations')

        ax[1].plot(lnIter, lfDur)
        ax[1].set_xlabel('Iterations')
        ax[1].set_ylabel('Duration (s)')
        plt.show()

# get full mono labels and split into monopole labels and jlabels
    def constructMonoLists(self, npaMonoList):
        npaFullMonoLab = self.mloNuclei[0].getMonoLabel(npaMonoList)
        lnpaShortMonoLab = [npaFullMonoLab[i][:4] for i in range(npaFullMonoLab.shape[0])]
        temp = []
        import numpy as np
        for i in range(len(lnpaShortMonoLab)):
            bAddIt = True
            for elem in temp:
                if np.all(lnpaShortMonoLab[i] == elem):
                    bAddIt = False
            if bAddIt:
                temp.append(lnpaShortMonoLab[i])
        import copy
        lnpaShortMonoLab = copy.copy(temp)
        lnpaMonoJLab = []

        for shortLab in lnpaShortMonoLab:
            temp = []
            for fullLab in npaFullMonoLab:
                if np.all(shortLab == fullLab[:-2]):
                    temp.append(fullLab[-2:])
                    lnpaMonoJLab.append(temp)
        return lnpaShortMonoLab, lnpaMonoJLab

#    fit the specified monopole terms
    def sMono(self, npaMonoList=[]):
        lnpaShortMonoLab, lnpaMonoJLab = self.constructMonoLists(npaMonoList)
        a = []
        import numpy as np
        for nucleus in self.mloNuclei:
            npaOcc = nucleus.getOcc(nucleus.getLevName())
            temp2 = []
            for row in npaOcc:
                temp1 = []
                for shortLab in lnpaShortMonoLab:
                    if shortLab[0] == shortLab[1]:
                        temp1.append(row[shortLab[0] - 1] *
                                     (row[shortLab[0] - 1] - 1.) / 2.)
                    elif shortLab[0] != shortLab[1]:
                        temp1.append(row[shortLab[0] - 1] *
                                     row[shortLab[1] - 1])
                if len(temp2) == 0:
                    temp2 = np.array(temp1)
                    temp2.shape = [1, len(temp2)]
                else:
                    temp1 = np.array(temp1)
                    temp1.shape = [1, len(temp1)]
                    temp2 = np.append(temp2, temp1, axis=0)
            if len(a) == 0:
                a = np.array(temp2, dtype=float)
            else:
                a = np.append(a, temp2, axis=0)
        import MatManip
        rmList = MatManip.getZeroCols(a)
# convert lnpaShortmonolab to nparray to use the rmslice utility in matmanip
        lnpaShortMonoLab = np.array(lnpaShortMonoLab)
        lnpaMonoJLab = np.array(lnpaMonoJLab)
        if rmList != []:
            a = MatManip.rmSlice(rmList, a, 1)
            lnpaShortMonoLab = MatManip.rmSlice(rmList, lnpaShortMonoLab, 0)
            lnpaMonoJLab = MatManip.rmSlice(rmList, lnpaMonoJLab, 0)
        npaETh = []
        for nucleus in self.mloNuclei:
            tempth = nucleus.getEnNu()
            if npaETh != []:
                npaETh = np.append(npaETh, tempth, axis=0)
            else:
                npaETh = tempth
        self.EExp.shape = [self.EExp.size, 1]
        npaETh.shape = [npaETh.size, 1]
        target = self.EExp - npaETh
        target.shape = [target.size, 1]
#     #start linear least squares
#     shortdiff=np.linalg.lstsq(a, target)
#     #end linear least squares
#    #start weighted least square
        npaWeights = np.zeros([self.EExp.size, self.EExp.size])
        for nIdx, elem in enumerate(self.npaErrors):
            npaWeights[nIdx, nIdx] = 1.0 / (elem**2 + self.fThError**2)
        shortdiff = np.linalg.lstsq(np.dot(npaWeights, a), np.dot(npaWeights,
                                                                  target))
#    #end weighted least squares
        shortdiff = shortdiff[0]
        npaLongMonoLab, longdiff = self.makeLong(lnpaShortMonoLab,
                                                 lnpaMonoJLab, shortdiff)
        for nucleus in self.mloNuclei:
            nucleus.llMESpec[1] = npaLongMonoLab
            nucleus.llMESpec[0] = []
            self.mloNuclei[0].setMEnum()
        npaME = self.mloNuclei[0].getME()
        npaME += longdiff
        self.nDOF = a.shape[0] - shortdiff.size - 1
        return (npaME, a, target, npaME, lnpaShortMonoLab, lnpaMonoJLab,
                shortdiff, self.mloNuclei[0].llMESpec[1])

# make the long labels and diff
    def makeLong(self, lnpaShortMonoLab, lnpaMonoJLab, shortdiff):
        '''
        Generate the full labels in order and assign the appropriate
        differences
        '''
        import numpy as np
        longdiff = []
        npaLongMonoLab = []
        for nIdx, shortLab in enumerate(lnpaShortMonoLab):
            for jLab in lnpaMonoJLab[nIdx]:
                if len(npaLongMonoLab) == 0:
                    npaLongMonoLab = np.array(np.append(shortLab, jLab))
                    npaLongMonoLab.shape = [1, npaLongMonoLab.size]
                    longdiff.append(shortdiff[nIdx])
                else:
                    temp = np.array(np.append(shortLab, jLab))
                    temp.shape = [1, temp.size]
                    npaLongMonoLab = np.append(npaLongMonoLab, temp, axis=0)
                    longdiff = np.append(longdiff, shortdiff[nIdx])
        if npaLongMonoLab.shape[0] != len(longdiff):
            print ('Error in make long: labels dont match differences: ',
                   str(npaLongMonoLab.shape[0]), 'labels and ',
                   str(len(longdiff)), ' differences.')
        return npaLongMonoLab, longdiff

# Make changes to each of the monopole matrix elements individually for
# different increments and then plot how the the expected and obtained energy
# changes depend on eachother and the increment.
    def checkMonoResponse(self, fIncLow=10.**-3, fIncHigh=1.0, nRuns=10,
                          display=True):
        import numpy as np
        npaIncs = np.linspace(fIncLow, fIncHigh, nRuns)
        import numpy as np
        npaEThOrig = []
        for nucleus in self.mloNuclei:
            tempth = nucleus.getEnNu()
            if npaEThOrig != []:
                npaEThOrig = np.append(npaEThOrig, tempth, axis=0)
            else:
                npaEThOrig = tempth
        (ans, a, target, npaOriginalME, lnpaShortMonoLab, lnpaMonoJLab,
         shortdiff) = self.sMono()
        for nIIdx, fInc in enumerate(npaIncs):
            for nJIdx, shortLab in enumerate(lnpaShortMonoLab):
                npaShortDiff = np.zeros([lnpaShortMonoLab.shape[0], 1])
# make the change for this iteration
                npaShortDiff[nJIdx] += fInc
        npaLongMonoLab, longdiff = self.makeLong(lnpaShortMonoLab,
                                                 lnpaMonoJLab, npaShortDiff)
        npaME = npaOriginalME + longdiff
        for nucleus in self.mloNuclei:
            nucleus.llMESpec[0] = []
            nucleus.llMESpec[1] = npaLongMonoLab
            nucleus.setMEnum()
# change the me
            nucleus.takeME(npaME)
            nucleus.runSM()
# get the obtained change
            npaEThNew = []
            llnAZ = []
            for nucleus in self.mloNuclei:
                tempth = nucleus.getEnNu()
                llnAZ.extend([nucleus.nAZ for i in range(len(tempth))])
                if npaEThNew != []:
                    npaEThNew = np.append(npaEThNew, tempth, axis=0)
                else:
                    npaEThNew = tempth
            npaChangeObtained = npaEThNew - npaEThOrig
            npaShortDiff.shape = [npaShortDiff.size, 1]
            npaChangeExpected = np.dot(a, npaShortDiff)
            self.writeMonoResponse(npaChangeExpected, npaChangeObtained, fInc,
                                   shortLab, llnAZ)
        self.displayMonoResponse(lnpaShortMonoLab)

# output the response for changes in the monopole term
    def writeMonoResponse(self, Eexpect, Eobtained, fInc, npaMonoLab, llnAZ):
        import os
        path = self.sOutPath + '\\tracking\\MonoResponse\\'
        sIsoPath = path + '\\ByIsotope\\'
        sMonoPath = path + '\\ByMonopole\\'
        if not os.path.isdir(path):
            os.makedirs(sIsoPath)
            os.makedirs(sMonoPath)
        sIsoForm = '{:10}{:10.4f}{:10.4f}{:10.4f}\n'
        sIsoHeadForm = '{:10}{:10}{:10}{:10}\n'
        for expect, obtain, isoLab in zip(Eexpect, Eobtained, llnAZ):
            sIsofName = 'A_' + str(isoLab[0]) + 'Z_' + str(isoLab[1]) + '.mr'
            if os.path.isfile(sIsoPath + sIsofName):
                fOut = open(sIsoPath + sIsofName, 'a')
            else:
                fOut = open(sIsoPath + sIsofName, 'w')
                fOut.write(sIsoHeadForm.format('Mono', 'Increment', 'Expected',
                                               'Obtained'))
            fOut.write(sIsoForm.format(str(npaMonoLab), float(fInc),
                                       float(expect), float(obtain)))
            fOut.close()
        sMonofName = str(npaMonoLab) + '.mr'
        sMonoForm = '{:10}{:10.4f}{:10.4f}{:10.4f}\n'
        sMonoHeadForm = '{:10}{:10}{:10}{:10}\n'
        if os.path.isfile(sMonoPath+sMonofName):
            fOut = open(sMonoPath+sMonofName, 'a')
        else:
            fOut = open(sMonoPath+sMonofName, 'w')
            fOut.write(sMonoHeadForm.format('Isotope', 'Increment', 'Expected',
                                            ' Obtained'))
        for expect, obtain, isoLab in zip(Eexpect, Eobtained, llnAZ):
            fOut.write(sMonoForm.format(str(isoLab), float(fInc),
                                        float(expect), float(obtain)))
        fOut.close()

# display the result of the monoME response test
    def displayMonoResponse(self, lnpaShortMonoLabel):
        path = self.sOutPath + '\\tracking\\MonoResponse\\'
        sMonoPath = path + '\\ByMonopole\\'
        import matplotlib.pyplot as plt
        import os
        import numpy as np
# plot by monopole term
        fig, ax = plt.subplots(2, 1)
        lfExpectedLong = []
        lfObtainedLong = []
        for npaMonoLab in lnpaShortMonoLabel:
            sMonofName = str(npaMonoLab) + '.mr'
            if os.path.isfile(sMonoPath + sMonofName):
                fOut = open(sMonoPath + sMonofName, 'r')
                lfExpected = []
                lfObtained = []
                lfInc = []
                nidx = 0
                for line in fOut:
                    temp = line.strip().split()
                    if nidx > 0:
                        lfInc.append(float(temp[2]))
                        lfExpected.append(float(temp[3]))
                        lfObtained.append(float(temp[4]))
                    nidx += 1
                    lfExpected = np.array(lfExpected)
                    lfObtained = np.array(lfObtained)
                    if len(lfExpectedLong) != 0 and len(lfObtainedLong) != 0:
                        lfExpectedLong = np.append(lfExpectedLong, lfExpected)
                        lfObtainedLong = np.append(lfObtainedLong, lfObtained)
                    else:
                        lfExpectedLong = lfExpected
                        lfObtainedLong = lfObtained
                ax[0].plot(lfExpected, lfObtained, 'o', label=npaMonoLab)
                lfRat = []
                lfIncR = []
                for fInc, fEx, fOb in zip(lfInc, lfExpected, lfObtained):
                    if fEx != 0:
                        lfRat.append(fOb / fEx)
                        lfIncR.append(fInc)
                        ax[1].plot(lfIncR, lfRat, 'o', label=str(npaMonoLab))
                    else:
                        print ('Error: displayMonoRespose: file ',
                               sMonoPath+sMonofName, ' does not exist.')
                        break
        from scipy import stats
        slope, intercept, r, p, stderr = stats.linregress(lfExpectedLong,
                                                          lfObtainedLong)
        x = ax[0].get_xlim()
        y = [slope * t + intercept for t in x]
        lform = 'Ref: ${:3.2f}x{:+3.2f}$, $r^2={:3.2f}$'

        ax[0].plot(x, y, '--', label=lform.format(slope, intercept, r**2))
        ax[0].set_title('Monpole Response by Monopole Term')
        ax[0].set_xlabel('Expected Change (MeV)')
        ax[0].set_ylabel('Change Obtained (MeV)')
        ax[0].legend()

        ax[1].set_xlabel('Increment (MeV)')
        ax[1].set_ylabel('Ratio of Obtained to Expected')
        ax[1].legend()
#     plt.tight_layout()
        plt.show()
# plot residual histogram
        temp = np.subtract(lfExpectedLong, lfObtainedLong)
        
        residual = [number for number in temp if x != 0.0]
        residual = np.absolute(residual)
        import math
        hist, bins = np.histogram(residual, bins=int(math.sqrt(residual.size)),
                                  density=False)
        print sum(hist), np.sum(hist)
        hist = hist / float(sum(hist))
        width = 0.7 * (bins[1] - bins[0])
        center = (bins[:-1] + bins[1:]) / 2
        plt.bar(center, hist, align='center', width=width)
        plt.xlabel('Residuals (MeV)')
        plt.ylabel('Relative Frequency')
        plt.title('Histogram of Residuals')
        plt.show()
# plot by Isotope
        fig, ax = plt.subplots(2, 1)
        sIsoPath = path + '\\ByIsotope\\'
        llnAZ = []
        for nucleus in self.mloNuclei:
            llnAZ.append(nucleus.nAZ)
            for isoLab in llnAZ:
                sIsofName = ('A_' + str(isoLab[0]) + 'Z_' + str(isoLab[1]) +
                             '.mr')
                if os.path.isfile(sIsoPath+sIsofName):
                    fOut = open(sIsoPath+sIsofName, 'r')
                    lfExpected = []
                    lfObtained = []
                    lfInc = []
                    nidx = 0
                    for line in fOut:
                        temp = line.strip().split()
                        if nidx > 0:
                            lfInc.append(float(temp[4]))
                            lfExpected.append(float(temp[5]))
                            lfObtained.append(float(temp[6]))
                        nidx += 1
                    lfExpected = np.array(lfExpected)
                    lfObtained = np.array(lfObtained)
                    sIsoLab = 'A_' + str(isoLab[0]) + 'Z_' + str(isoLab[1])
                    ax[0].plot(lfExpected, lfObtained, 'o', label=sIsoLab)
                lfRat = []
                lfIncR = []
                for fInc, fEx, fOb in zip(lfInc, lfExpected, lfObtained):
                    if fEx != 0:
                        lfRat.append(fOb / fEx)
                        lfIncR.append(fInc)
                        ax[1].plot(lfIncR, lfRat, 'o', label=sIsoLab)
                    else:
                        print ('Error: displayMonoRespose: file ',
                               sIsoPath+sIsofName, ' does not exist.')
                        break
                ax[0].plot(x, y, '--', label=lform.format(slope, intercept,
                                                          r**2))
                ax[0].set_title('Monpole Response by Isotope')
                ax[0].set_xlabel('Expected Change (MeV)')
                ax[0].set_ylabel('Change Obtained (MeV)')
                ax[0].legend()

                ax[1].set_xlabel('Increment (MeV)')
                ax[1].set_ylabel('Ratio of Obtained to Expected')
                ax[1].legend()
#     plt.tight_layout()
                plt.show()

# code to update a report on each iteration
    def sMonoIterationReport(self, Output, npaShortMono):
        '''
            Output=[ 0 ans, 1 a, 2 target, 3 npaME,4 lnpaShortMonoLab,
                    5 lnpaMonoJLab, 6 shortdiff]
        '''
        sMonoOutPath = self.sOutPath + '\\sMonoIteration\\'
# output the matrix
        import os
        if not os.path.isdir(sMonoOutPath):
            os.makedirs(sMonoOutPath)
# write the llsq matrix for the current iteration
        fOut = open(sMonoOutPath + 'Matrix.dat', 'w')
        sFormat = '{:10.4f}'
        for line in Output[1]:
            sLine = ''
            for elem in line:
                sLine += sFormat.format(elem)
            fOut.write(sLine + '\n')
        fOut.close()
#     Write the monodiff
        if not os.path.isfile(sMonoOutPath+'MonoDiff.dat'):
            fOut = open(sMonoOutPath + 'MonoDiff.dat', 'w')
            line = ''
            for elem in Output[4]:
                line += '{:>10}'.format(str(elem))
            fOut.write(line + '\n')
        else:
            fOut = open(sMonoOutPath + 'MonoDiff.dat', 'a')
            line = ''
            temp = Output[6]
            temp.shape = [Output[6].size]
            for elem in temp:
                line += sFormat.format(elem)
            fOut.write(line + '\n')
        fOut.close()
# write the expected energy change and the obtained energy change
        if not os.path.isfile(sMonoOutPath + 'En.dat'):
            sInPath = self.sOutPath+'\\tracking\\Levels.dat'
            fIn = open(sInPath, 'r')
            lfLastEn = []
            nLine = 0
            for line in fIn:
                line = line.strip().split()
                if nLine > 0:
                    lfLastEn.append(float(line[5]))
                nLine += 1
            fIn.close()
        else:
            fIn = open(sMonoOutPath + 'En.dat', 'r')
            lfLastEn = []
            nLine = 0
            for line in fIn:
                line = line.strip().split()
                if nLine > 0:
                    lfLastEn.append(float(line[0]))
                nLine += 1
        lfNewEn = []
        for nucleus in self.mloNuclei:
            lfNewEn.extend(nucleus.getEnNu())
        fIn.close()
        fOut = open(sMonoOutPath + 'En.dat', 'w')
        fOut.write('Previous iteration Energy \n')
        for elem in lfNewEn:
            fOut.write(sFormat.format(elem) + '\n')
        fOut.close()
        fOut = open(sMonoOutPath + 'EnDiff.dat', 'w')
        sHFormat = '{:10}{:10}{:10}\n'
        fOut.write(sHFormat.format('Expected', 'Obtained', 'diff'))
        import numpy as np
        npaDiffObtained = np.subtract(lfNewEn, lfLastEn)
        npaDiffExpect = np.dot(Output[1], Output[6])
        sFormat = '{:10.4f}{:10.4f}{:10.4f}'
        for expect, obtain in zip(npaDiffExpect, npaDiffObtained):
            line = sFormat.format(float(expect), float(obtain),
                                  float(expect) - float(obtain))
            fOut.write(line + '\n')
        fOut.close()
        fOut = open(sMonoOutPath + 'combined.dat', 'w')
        lLS = []
        lAZ = []
        for nucleus in self.mloNuclei:
            lLS.extend(nucleus.mllspec)
            lAZ.extend([nucleus.nAZ]*len(nucleus.mllspec))
        sHFormat = '{:10}{:5}{:5}{:5}{:>10}'
        sHExt = '{:>15}'
        sLine = sHFormat.format('[A,Z]', 'J', 'nJ', 'P', 'Exp-Theory')
        for label in npaShortMono:
            sLine += sHExt.format(label)
        fOut.write(sLine + '\n')
        sFormat = '{:10}{:5}{:5}{:5}{:10.4f}'
        sFormExt = '{:>15.4f}'
# print lAZ, lLS, npaEExp.shape, npaETh.shape
        for AZ, LS, Exp, Th, row in zip(lAZ, lLS, self.EExp, lfNewEn,
                                        Output[1]):
            sLine = sFormat.format(AZ, LS[0], LS[1], LS[2], float(Exp) -
                                   float(Th))
            for elem in row:
                sLine += sFormExt.format(elem)
            fOut.write(sLine + '\n')
        fOut.close()

    def multiMono(self, nMaxIter, fTolIn, lnpaBases):
        import numpy as np
        npaME = []
        npaTBMESpec = []
        for npaBase in lnpaBases:
            temp = self.IterativeLSq(sMethod='smono', bMix=False,
                                     nMaxIter=nMaxIter, fTolin=fTolIn,
                                     methodArg=npaBase)
            if len(npaME) != 0:
                npaME = np.append(npaME, temp[1], axis=0)
                npaTBMESpec = np.append(npaTBMESpec, temp[-1][1])
            else:
                npaME = temp[1]
                npaTBMESpec = temp[-1][1]
        return temp[0], npaME, npaTBMESpec

    def compositeSingleMono(self, nMaxIter, fTolIn, lnpaBases):
        import numpy as np
        temp = self.IterativeLSq(sMethod='single', bMix=False,
                                 nMaxIter=nMaxIter, fTolin=fTolIn)
        llMESpec = temp[-1]
        npaME = temp[1]
        temp = self.multiMono(nMaxIter, fTolIn, lnpaBases)
        llMESpec.append(temp[-1])
        npaME = np.append(npaME, temp[1])
        res = temp[0]
        return res, npaME, llMESpec

    def optimizeUniformError(self,  sMethod='smono', lMethodArgs=[],
                             bReset=True, sOptMethod='simple', fTolIn=10**-2):
        if sOptMethod == 'simple':
            while abs(1.0 - self.fLastChi) > fTolIn:
                self.IterativeLSq(sMethod, bMix=False, nMaxIter=60,
                                  fTolin=fTolIn, methodArg=lMethodArgs)
                if bReset:
                    for nucleus in self.mloNuclei:
                        nucleus.reinitialize
                if self.fLastChi != 0:
                    self.fThError = self.fThError * self.fLastChi
        else:
            print 'Invalid method specification.'
        sFormat = '{:+3.2e}'
        print 'The final Chi Square is: 1', sFormat.format(self.fLastChi - 1.0)
        sFormat = '{:3.2e}'
        print ('The final unform error on the matrix elements is: ',
               sFormat.format(self.fThError))

# make histogram of the errors in the levels.dat file. 
# does so for initial and final errors.
    def makeErHist(self, bReduce=False):
        fIn = open(self.sOutPath + '\\tracking\\Levels.dat', 'r')
        lFEr = []
        lIEth = []
        for nIdx, line in enumerate(fIn):
            if nIdx > 0:
                line = line.strip().split()
                lFEr.append(line[-1])
                lIEth.append(line[6])
        import matplotlib.pyplot as plt
        import math
        lGood = [nIdx for nIdx, e in enumerate(lIEth) if e != 'N/A']
        if bReduce:
            lBad = []
            for nIdx, er in enumerate(self.npaErrors):
                if er >= 1.0:
                    lBad.append(nIdx)
            lGood = [nIdx for nIdx in lGood if nIdx not in lBad]
        print max([self.npaErrors[nIdx] for nIdx in lGood])
        lFEr = [abs(float(lFEr[nIdx])) for nIdx in lGood]
        lIEr = [abs(float(lIEth[nIdx]) - float(self.EExp[nIdx])) for nIdx in lGood]
        nBins = int(math.floor(math.sqrt(len(lFEr))))
        fMin = min([min(lFEr), min(lIEr)])
        fMax = max([max(lFEr), max(lIEr)])
        fRMSI = math.sqrt(sum([er**2 for er in lIEr])/float(len(lIEr))) 
        fRMSF = math.sqrt(sum([er**2 for er in lFEr])/float(len(lFEr)))
        fakenumI = math.sqrt(sum([er**2 for er in lIEr]))/float(len(lIEr)) 
        fakenumF = math.sqrt(sum([er**2 for er in lFEr]))/float(len(lFEr))
        print 'initial', fRMSI, 'fake', fakenumI
        print 'final', fRMSF, 'fake', fakenumF
        print max(lIEr), max(lFEr)
        fInc = (fMax - fMin) / nBins
        lBins = [fMin + nIdx * fInc for nIdx in range(nBins + 1)]
        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True)
        ax1.hist(lIEr, lBins)
        ax1.set_title('Initial Error Distribution')
        ax1.plot()
        ax1.plot()
        ax2.hist(lFEr, lBins)
        ax2.set_title('Final Error Distribution')
        ylims = plt.ylim()
        ax1.plot([fRMSI, fRMSI], ylims, 'r--')
        ax1.plot([fRMSF, fRMSF], ylims, 'k--')
        ax2.plot([fRMSI, fRMSI], ylims, 'r--', label='Initial RMSE')
        ax2.plot([fRMSF, fRMSF], ylims, 'k--', label='Final RMSE')
        ax2.legend()
        plt.tight_layout()
        plt.savefig(self.sOutPath + '\\tracking\\ErrorHist')
        plt.show()

    def testScaling(self, sOutpath, nMaxNumProcesses=26):
        i = -1
        for nIdx, char in enumerate(self.sInPath):
            if char == '\\':
                i = nIdx
        spath = self.sInPath[:i]
        import os
        os.chdir(self.sOutPath)
        from subprocess import call
        import time     
        sForm = '{:5d}{:10.2e}\n'
        for nIdx in range(nMaxNumProcesses):
            fStart = time.time() 
            print 'num proc=' + str(nIdx+1)
            call('mpiexec -n '+ str(nIdx + 1) + ' python ' + spath +'\\OxbashMPI.py')
            fEnd = time.time()
            fTime = fEnd - fStart
            fOut = open(sOutpath, 'a')
            fOut.write(sForm.format(nIdx, fTime))
            fOut.close()

    def parseScalingOutput(self, sOutPath):
        '''
            Take the output files of testScaling and read in the data to the output
            variables. 
        '''
        fScaling = open(sOutPath, 'r')
        lnIdxs = []
        lfTimes = []
        for line in fScaling:
            temp = line.strip().split()
            lnIdxs.append(int(temp[0]))
            lfTimes.append(float(temp[1]))
        fScaling.close()
        fSums = open(self.sOutPath+'\\SumFile.dat', 'r')
        temp = ''
        lSums = []
        for line in fSums:
            temp = temp + line.strip('\n') + ' '
            if temp[-2] == ']':
                temp = temp.strip(' ')
                temp = temp.strip('[')
                temp = temp.strip(']').split()
                temp1 = [float(elem) for elem in temp]                
                lSums.append(temp1)
                temp = ''
        fSums.close()
        fTimes = open(self.sOutPath + '\\times.dat', 'r')
        lTimes = []
        lWUID = []
        for line in fTimes:
            line = line.strip().split()
            lWUID.append(int(line[0]))
            lTimes.append(float(line[1]))
        fTimes.close()
        return [lnIdxs, lfTimes], lSums, [lWUID, lTimes]

    def plotScalingOutput(self, plotDest='', WUdescriptor='',
                          WdistDescriptor='', scalingIntervals=[0, 12, 24]):
        '''
            Read the scaling output and produce plots to analyze the results.
        '''
        llScaling, llSums, llTimes = self.parseScalingOutput('c:\\PythonScripts\\OxBashWork\\test\\scaling.dat')
        fMaxSRT = max(llTimes[1])
        temp = [max(elem) for elem in llSums]
        fMaxSum = max(temp)
        temp = [min(elem) for elem in llSums]
        fMinSum = min(temp)
        import matplotlib.pyplot as plt
        '''
            plot scaling start by calculating regressions of the indicated intervals
        '''
        from scipy.stats import linregress as linreg
        x = [elem for elem in llScaling[0][:]]
        nMY = max(llScaling[1][1:])
        y = [nMY/elem for elem in llScaling[1][:]]
        fig0 = plt.figure()
        ax0 = fig0.add_subplot(111)
        ax0.plot(x, y, ls='', marker='o', label='Measured')
        sForm = 'R{:d}:a={:3.2f},b={:3.2f},$r^2$={:3.2f}'
        for nIdx in range(len(scalingIntervals) - 1):
            start = scalingIntervals[nIdx]
            stop = scalingIntervals[nIdx + 1]
            xlims = [x[start], x[stop]]
            temp = linreg(x[start:stop], y[start:stop])
            slope = temp[0]
            inter = temp[1]
            r = temp[2]
            ax0.plot(xlims, [xlims[0]*slope + inter, xlims[1] * slope + inter],
                     ls='--', label=sForm.format(nIdx + 1, slope, inter, r**2))
        ax0.set_title('Speed Up')
        ax0.set_ylabel('$t_{max}/{t_p}$', fontsize=14)
        ax0.set_xlabel('# of Processes')
        ax0.legend(loc='best')
        fig0.savefig(plotDest + 'Scaling.pdf')    
        '''
            Plot Sums
        '''
        fig1 = plt.figure(figsize=[12,16.5])
        ax = fig1.add_subplot(111)
        ax1 = fig1.add_subplot(211)
        ax2 = fig1.add_subplot(212)
        ax.spines['top'].set_color('none')
        ax.spines['bottom'].set_color('none')
        ax.spines['left'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
        lMarks = ['*','+']
        sForm = 'P={:d}'
        for i, series in enumerate(llSums[:12]):
            ax1.semilogy(series, marker=lMarks[i % 2], ls='', label=sForm.format(i + 1))
        lim = ax1.get_xlim()
        ax1.semilogy(lim, [fMaxSum, fMaxSum], ls='--', color='k', label='SRT')
        ax1.semilogy(lim, [fMaxSRT, fMaxSRT], ls='--', color='b', label='MaxWU')
        for i, series in enumerate(llSums[12:]):
            ax2.semilogy(series, marker=lMarks[i % 2], ls='', label=sForm.format(i + 12))
        lim = ax2.get_xlim()               
        ax2.semilogy(lim, [fMinSum, fMinSum], ls='--', color='k', label='MinWL')
        ax2.semilogy(lim, [fMaxSRT, fMaxSRT], ls='--', color='b', label='MaxWU')
        sTitle = 'WorkLoad Balancing'
        if WdistDescriptor != '':
            sTitle += 'for' + WdistDescriptor
        ax.set_title(sTitle)
        ax.set_xlabel('Process ID', fontsize=14)
        ax.set_ylabel('Serial Runtime of Workload (s)', fontsize=14)
        ax1.legend(loc='best')
        ax2.legend(loc='best')
        plt.tight_layout()
        fig1.savefig(plotDest + 'WBal.pdf')
        '''
            plot times
        '''
        from math import sqrt
        nBins = int(4 * sqrt(float(len(llTimes[1]))))
        fig2 = plt.figure()
        ax = fig2.add_subplot(111)
        ax.hist(llTimes[1], nBins)
        sTitle = 'Distribution of Work Unit Sizes'
        if WUdescriptor != '':
            sTitle += ' for ' + WUdescriptor
        ax.set_title(sTitle)
        ax.set_xlabel('Serial Runtime (s)')
        ax.set_ylabel('Frequency')
        fig2.savefig(plotDest + 'WUDist.pdf')
    
    def reInit(self):
        '''
            Reinitialize the output file using the path to a saved 
            initialization if it exists.
            Warning: Deletes the previous work in the output file so be sure to 
            save anything you want before using this method! 
        '''
        import os
        os.chdir(self.sOutPath+'\\..')
        import shutil
        shutil.rmtree(self.sOutPath)
        if self.sInitDir != '':
            print 'Copying the initializaton from:"', self.sInitDir,'"...'
            shutil.copytree(self.sInitDir, self.sOutPath)
            print 'Copy complete!'
        else:
            print 'Error: no initialization directory provided'

    def testLincomSingle(self, sOutDir, lnRng=[]):
        '''
            Test how the number of linear combinations fit in the optimization
            affects the residue of the first iteration. Take the baseline from
            the initialization, and iteratet through the range of linear 
            combinations provided. Write results to an output file in the 
            specified directory.
        '''
        import os
        if not os.path.isdir(sOutDir):
            os.makedirs(sOutDir)
        for n in lnRng:
            ans, a, target, npaME = self.TBTDLeastSq(methodArg=[n, 'fd'])
            fRes = self.obj(ans)
            fOut = open(sOutDir+'\\TLCSingle.dat', 'a')
            sForm = '{:>10d}{:>10.5f}\n'
            fOut.write(sForm.format(n, fRes))
            fOut.close()
            self.reInit()

    def bootstrappedLC(self, sOutDir, lnRng=[], bSave=True):
        '''
            Runs the iterative least square to convergence for each of a given 
            list of linear combinations. Use previous result as initialization 
            for next result and save each converged fit.  
        '''
        import os
        import shutil
        if not os.path.isdir(sOutDir):
            os.makedirs(sOutDir)
        for n in lnRng:
            self.IterativeLSq(sMethod='TBTD', bMix=False, nMaxIter=10, fTolin=10**-2,methodArg=[n, 'fd'])
            if bSave:
                shutil.copytree(self.sOutPath, sOutDir+'\\'+str(n))

    def testLincomNone(self, sOutDir, lnRng=[]):
        '''
            Test how the number of linear combinations fit in the optimization
            affects the residue of the result. Take the baseline from
            the initialization, and iteratet through the range of linear 
            combinations provided. Write results to an output file in the 
            specified directory.
        '''
        import os
        if not os.path.isdir(sOutDir):
            os.makedirs(sOutDir)
        for n in lnRng:
            ans, a, target, npaME = self.TBTDLeastSq(methodArg=[n, 'fd'])
            fRes = self.obj(ans)
            fOut = open(sOutDir+'\\TLCSingle.dat', 'a')
            sForm = '{:>10d}{:>10.5f}\n'
            fOut.write(sForm.format(n, fRes))
            fOut.close()
            
    def tlcPlot(self, sTLCDat, fRef = 0.169):
        fDat = open(sTLCDat, 'r')
        lfLC=[]
        lfDat=[]
        for line in fDat:
            line = line.strip().split()
            lfLC.append(float(line[0]))
            lfDat.append(float(line[1]))
        print len(lfLC), len(lfDat)
        print lfLC, lfDat
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(lfLC, lfDat, label = 'LinCom Data')
        ax.set_xlabel('# linear combinations fit')
        ax.set_ylabel('RMS error after 1 iteration (MeV)')
        ax.set_ylim([0.,.2])
        ax.set_title('usda w/sdba background')
        fXlim = ax.get_xlim()
        ax.plot(fXlim, [fRef,fRef], ls='--', label='sdba w/o fit')
        plt.legend(loc='best')
        plt.savefig(sTLCDat+'tlc.pdf')
    
    def twoInteractionCorrelation(self, lsIntPaths, lsIntLabs):
        if len(lsIntPaths) !=2 or len(lsIntLabs)!=2:
            print 'Error twoInteractionCorrelation: needs 2 paths and two labels.'
            return None
        llfOBME = [self.mloNuclei[0].getOBME(sBGIntPath=lsIntPaths[i], bAll=True) for i in range(len(lsIntPaths))]
        llfTBME = [self.mloNuclei[0].getTBME(sBGIntPath=lsIntPaths[i], bAll=True) for i in range(len(lsIntPaths))]
        import matplotlib.pyplot as plt
        fig,ax = plt.subplots(1,1)
        ax.plot(llfOBME[0][0][:],llfOBME[1][0][:], color= 'k', ls='', marker = 'o',label='SPE')
        ax.plot(llfTBME[0], llfTBME[1], color= 'b', ls='', marker = '*', label='TBME')
        print llfOBME
        temp=list(llfOBME[0][0][:])
        temp.extend(list(llfTBME[0]))
        print temp
        temp.sort()
        ax.plot(temp, temp, color= 'r', ls='--', marker = '',label='ideal')
        ax.set_xlabel(lsIntLabs[0])
        ax.set_ylabel(lsIntLabs[1])
        ax.legend(loc='best')
        fig.tight_layout()
        plt.show()
        return None
        
    def intSideBySide(self, sPathNew, sPathOld, sPathComp):
        lOldOBME=self.mloNuclei[0].getOBME(sBGIntPath=sPathOld,bAll=True)
        lNewOBME=self.mloNuclei[0].getOBME(sBGIntPath=sPathNew,bAll=True)
        lOldTBME=self.mloNuclei[0].getTBME(sBGIntPath=sPathOld,bAll=True)
        lNewTBME=self.mloNuclei[0].getTBME(sBGIntPath=sPathNew,bAll=True)
        
        fComp = open(sPathComp, 'w')
        sHForm = '{:22}{:>10}{:>10}{:>10}{:>10}\n'
        sDatForm = '{:22}{:10.4f}{:10.4f}{:10.4f}{:10.1f}%\n'
        fComp.write(sHForm.format('Label', 'Old', 'New', 'Diff', '%diff'))
        fComp.write('Single Particle Energies\n')
        diffs=[]
        for lLab, fOld, fNew in zip(self.mloNuclei[0].llMESpec[0], lOldOBME[0], lNewOBME[0]):
            fOld = float(fOld)
            fNew = float(fNew)
            diffs.append(fNew - fOld)
            fComp.write(sDatForm.format(str(lLab), fOld, fNew, fNew-fOld, 200.*abs(fOld-fNew)/(abs(fOld)+abs(fNew))))
        fComp.write('Two Body Matrix Elements\n')
        for lLab, fOld, fNew in zip(self.mloNuclei[0].getLabel(), lOldTBME, lNewTBME):
            fOld = float(fOld)
            fNew = float(fNew)
            diffs.append(fNew - fOld)
            fComp.write(sDatForm.format(str(lLab), fOld, fNew, fNew-fOld, 200.*abs(fOld-fNew)/(abs(fOld)+abs(fNew))))
        fComp.close()
        return diffs
        
# #########################################################
'''
Start testing code
'''
import sys
sys.path.append('c:\\PythonScripts\\OxBashScripts\\')
sys.path.append('C:\\PythonScripts\\generalmath')
#initialization = 'E:\\initdir\\Lodai_usda_sdba'
#initialization = 'E:\\initdir\\Lodai_usda_srg24_fit\\20'
#initialization = 'E:\\initdir\\shorttest'
#initialization = 'E:\\initdir\\uptoA20ex'
#initialization = 'E:\\initdir\\Lodai_imsrg_A24'
initialization =''
#
#sPathOld = 'C:\\oxbash\\sps\\usdb.int'
#sPathNew = 'C:\\PythonScripts\\OxBashWork\\Lodai_usdb\\A17_Z8\\usdb.int'
#sPathComp = 'C:\\PythonScripts\\OxBashWork\\Lodai_usdb\\tracking\\usdb_comp.dat'
x = BashOpt('c:\\PythonScripts\\OxBashScripts\\OptInput-shorttest.in',
            'c:\\PythonScripts\\OxbashWork\\shorttest',
            'c:\\PythonScripts\\OxBashScripts\\errors.dat', initialize=True, 
            serial = False, sInitDir=initialization)
#import shutil
#shutil.copytree('c:\\PythonScripts\\OxbashWork\\w_sdba_fit','e:\\initdir\\w_sdba_fit')
#lnLC=range(5,66,5)
#x.bootstrappedLC('E:\\initdir\\Lodai_w_sdba_fit\\', lnLC)
