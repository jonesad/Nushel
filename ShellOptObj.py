# -*- coding: utf-8 -*-

"""
Created on Wed Jan 14 14:03:24 2015

attempt at object orientation of nushell optimization

@author: jonesad
"""
# writes an example inputfile for the ShellOpt Class for single particle o20
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
    fInFile.write('sdpn\n')
# restriction
    fInFile.write('n\n')
# interaction
    fInFile.write('usdcpn\n')
# formalism iso/pn
    fInFile.write('pn\n')
# does the interaction extrapolate matrix elements
    fInFile.write('False\n')

# number of nuclei used in optimization
    fInFile.write('10\n')

# number of protons
    fInFile.write('  8\n')
# of nucleons
    fInFile.write('17\n')
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
CreateInFile('C:/PythonScripts/NushellScripts/OptInput.in')


class ShellOpt:
    'Class for shell model hamiltonian optimization problems'
    def __init__(self, sInPath, sOutPath, sErrorPath='', initialize=True,
                 fThError=0.07885):
        self.nDOF = 1
        self.fLastChi = 0
# flag that says to use the groundstate energy in optimization
        self.init = False
        self.useGS = 1
# flag that tracks the residue and energy levels over the iterations
        self.sMethod = ''
        self.track = 1
        self.sInPath = sInPath
        self.sOutPath = sOutPath
        self.EExp = []
        self.npaErrors = []
# get info from input file and initializes the nuclei objects
        self.GetIn(initialize)
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

            if initialize:
                self.obj(self.mloNuclei[0].getME())
            self.init = True
            self.EExp.shape = [self.EExp.size, 1]

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
        sFormat = '{:10}{:5}{:5}{:5}{:10.4f}{:10.4f}'
# print lAZ, lLS, npaEExp.shape, npaETh.shape
        for AZ, LS, Exp, Th in zip(lAZ, lLS, self.EExp, npaETh):
            fOut.write(sFormat.format(AZ, LS[0], LS[1], LS[2], Exp, Th) + '\n')
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
        sFormat1 = '{:5}{:5}{:5}{:5}{:5}'
        sFormat2 = '{:>10}'
        sFormatNew = '{:>10.4f}{:>10.4f}'
        nIdx = 0
        for line in fIn:
            if nIdx > -1:
                temp = line.strip().split()
                sNewLine = sFormat1.format(temp[0], temp[1], temp[2], temp[3],
                                           temp[4], temp[5])
                for nJIdx in range(len(temp)-6):
                    sNewLine += sFormat2.format(temp[nJIdx + 6])
                    sNewLine += sFormatNew.format(npaETh[nIdx], npaETh[nIdx] -
                                                  float(temp[6]))
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
        fTol = fTolin / float(len(self.mloNuclei[0].mllspec))
        lME = [[], [], [], []]

        if sMethod == 'mono':
            npaMonoLabel = self.mloNuclei[0].getMonoLabel()
            for nucleus in self.mloNuclei:
                nucleus.llMESpec = [list(range(1, nucleus.countOBME() + 1)),
                                    npaMonoLabel]
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
            elif sMethod == 'csm':
                temp = self.compositeSingleMono(nMaxIter, fTolin, methodArg)
                npaGuess = temp[1]
                for nucleus in self.mloNuclei:
                    nucleus.llMESpec = temp[2]
            else:
                print 'Error: sMethod=', sMethod, 'is not a valid argument.'
                break
            print "The guess is", npaGuess
            fResNew = self.obj(npaGuess)
            if type(fResNew) != float:
                print ('Error: instance of ''shellopt'' obj method ' +
                       'returning invalid type: '), type(fResNew)
            if sMethod == 'smono':
                self.sMonoIterationReport(temp)
            lRes[0] = fResNew
            lME[0] = npaGuess
            nIter += 1
# restore the original MESpec
            for nucleus in self.mloNuclei:
                nucleus.llMESpec = list(llOriginalMESpec)
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
        lNewSPList = range(self.mloNuclei[0].countOBME())
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
                                     np.zeros(lnpaShortMonoLab.shape))
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
            npaGuess = np.append(npaGuess, np.zeros(newTBME.shape))
            res = minimize(self.objSMonoNonLinear, npaGuess, method=sMethod,
                           options=dOptions)
        print res
        self.updateLevs(self.sOutPath+'\\'+'tracking\\')

#   non linear objective function for fitting the summed monopole matrix
#   elements. Take the npa uess composed of the single particle energies 
#   followed by the monopole matrix element differences. Then use the short 
#   mono lab and the jlab to make the long list of differneces. Then use that
#    to make the guess for the original objective function
    def objSMonoNonLinear(self, npaGuess, npaShortmonoLab, npaJLab):
        
    def GetIn(self, initialize=True):
        import ShellNuclei
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
            self.mnNuclei = int(fIn.readline().strip('\n'))
            self.mloNuclei = []
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
                templine = fIn.readline().strip().split()
                fGSE = float(templine[0])
                fError = float(templine[1])
                sMMDR = fIn.readline().strip('\n')
                sPar = fIn.readline().strip('\n')
                tempnuc = ShellNuclei.nucleus(nZ, nA, fGSE, fError,
                                              self.sOutPath, sMMDR, sPar,
                                              self.lsShared, llMESpec,
                                              self.useGS, self.sForm,
                                              initialize, self.bExtrap)
                self.mloNuclei.append(tempnuc)
                llStateSpec = []
                temp = int(fIn.readline())
                for iii in range(temp):
                    line = fIn.readline().strip('\n')
                    line = line.split()
                    llStateSpec.append(line[:3])
                    self.EExp.append(line[3])
                    self.npaErrors.append(line[4])
                self.mloNuclei[nIdx].setLevels(llStateSpec)
                self.mloNuclei[nIdx].setmanBody(anBody)
        import numpy as np
        self.npaErrors = np.array(self.npaErrors, dtype=float)
        self.EExp = np.array(self.EExp, dtype=float)
        fIn.close()

# The objective function takes the array of matrix elements
    def obj(self, npaME):
        import numpy as np
        res = []
        for oNuc in self.mloNuclei:
            '''
                write ME to file
            '''
            if self.init:
                try:
                    oNuc.takeME(npaME)
                except:
                    print npaME
                    print 'Error: The Matrix element list is invalid.'
                    print npaME
                    return 'Error: The Matrix element list is invalid.'
            if self.sForm == 'pn':
                import os
                sLevName = oNuc.getLevName()
                temp = (str(oNuc.sPath) + '\\' + str(oNuc.sName) + '\\' +
                        sLevName[: -5] + '0' + '.int')
                if os.path.isfile(temp + '_'):
                    os.remove(temp + '_')
                    os.rename(temp, temp + '_')
# run Shell model calc
            oNuc.runSM()
# get the energy difference for the releveant levels
            res.extend(oNuc.getEnNu())
        res = np.array(res)
        res.shape = [res.size, 1]
        self.EExp.shape = [self.EExp.size, 1]
        res = self.EExp - res
        temp = np.sqrt(np.dot(np.transpose(res), res) / float(len(res)))
        temp.shape = [temp.size]
        temp = float(temp[0])

        fChiSq = []
        for num, err in zip(res, self.npaErrors):
            fChiSq.append(num / np.sqrt(err**2 + self.fThError**2))
        fChiSq = np.array(fChiSq)
        fChiSq.shape = [fChiSq.size, 1]
        fChiSq = (np.sqrt(np.dot(np.transpose(fChiSq), fChiSq)) /
                  float(self.nDOF))
        fChiSq = float(fChiSq[0, 0])

        if self.track == 1:
            self.OptStatus(temp, fChiSq, npaME)
        self.fLastChi = fChiSq
        return temp

    def OptStatus(self, res, fChiSq, npaME):
        sFormat = '{:14.10f}'
        fOut = open(self.sOutPath + '\\tracking\\res.dat', 'a+')
        fOut.write(sFormat.format(res) + '\n')
        fOut.close()
        if self.init:
            fOut = open(self.sOutPath + '\\tracking\\chisq.dat', 'a+')
            fOut.write(sFormat.format(fChiSq) + '\n')
            fOut.close()
        fOut = open(self.sOutPath+'\\tracking\\ME.dat', 'a+')
        string = ''
        sFormME = '{:10.4f}\t'
        for elem in npaME:
            string = string + sFormME.format(float(elem))
            fOut.write(string + '\n')
        fOut.close()
        fOut = open(self.sOutPath + '\\tracking\\AllME.dat', 'a+')
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
            fIn = open(nucleus.sPath + '\\' + nucleus.sName + "\\list.lpt")
            sLevName = fIn.readline().strip()
            fIn.close()
            if a != []:
                a = numpy.append(a, nucleus.getOcc(sLevName), axis=0)
            else:
                a = numpy.array(nucleus.getOcc(sLevName))
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
        for nucleus in self.mloNuclei:
            nucleus.llMESpec[0] = list(range(1, self.mloNuclei[0].countOBME() +
                                             1))
        lLabSpec = list(chain.from_iterable(self.mloNuclei[0].llMESpec))
        if self.sForm == 'iso':
            for num in range(3):
                try:
                    rmList.remove(num)
                except TypeError:
                    ''
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
        lnpaShortMonoLab = [npaFullMonoLab[:, :4]
                            for i in range(npaFullMonoLab.shape[0])]
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
        self.nDOF = a.shape[0] - shortdiff.size
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
                fOut = open(sIsoPath + sIsofName, 'a+')
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
            fOut = open(sMonoPath+sMonofName, 'a+')
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
    def sMonoIterationReport(self, Output):
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
            fOut = open(sMonoOutPath + 'MonoDiff.dat', 'a+')
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

# #########################################################
'''
Start testing code
'''
import sys
sys.path.append('c:\\PythonScripts\\NushellScripts\\')
sys.path.append('C:\PythonScripts\generalmath')

x = ShellOpt('c:\\PythonScripts\\NushellScripts\\OptInput.in',
             'c:\\PythonScripts\\NushellScripts\\test',
             'c:\\PythonScripts\\NushellScripts\\errors.dat', initialize=True)
#x.checkMonoResponse(fIncLow=0.1, fIncHigh=0.1,nRuns=1,display=True)

#ans, a, target, npaME,lnpaShortMonoLab,lnpaMonoJLab, shortdiff=x.sMono()
#x.displayMonoResponse(lnpaShortMonoLab)

#print x.CheckConvergenceSensitivity(10, 'smono', display=True)

#x.addMENoise(3.0)

#import numpy as np
#base=np.array([[5,5,5,5],[4,4,4,4,],[6,6,6,6]])
#
#x.optimizeUniformError(sMethod='smono', lMethodArgs=base, bReset=True,
#                       sOptMethod='simple', fTolIn=10**-2)


import numpy as np
print x.IterativeLSq(sMethod='single',bMix=False, nMaxIter=60, fTolin=10**-2)
base=np.array([[5,5,5,5],[4,4,4,4,],[6,6,6,6]])
print x.IterativeLSq(sMethod='smono',bMix=False, nMaxIter=60, fTolin=10**-2,methodArg=base)
base=np.array([[5,6,5,6], [5,4,5,4],[4,6,4,6]])
print x.IterativeLSq(sMethod='smono',bMix=False, nMaxIter=60, fTolin=10**-2,methodArg=base)
print x.IterativeLSq(sMethod='single',bMix=False, nMaxIter=60, fTolin=10**-2)
base=np.array([[5,5,5,5],[4,4,4,4,],[6,6,6,6]])
print x.IterativeLSq(sMethod='smono',bMix=False, nMaxIter=60, fTolin=10**-2,methodArg=base)
base=np.array([[5,6,5,6], [5,4,5,4],[4,6,4,6]])
print x.IterativeLSq(sMethod='smono',bMix=False, nMaxIter=60, fTolin=10**-2,methodArg=base)
x.plotResults(sMethod='smono', bError=True)

#
#import numpy as np
#lnpaBases=[np.array([[5,5,5,5],[4,4,4,4,],[6,6,6,6]]),np.array([[5,6,5,6], [5,4,5,4],[4,6,4,6]])]
#print x.IterativeLSq(sMethod='csm',bMix=False, nMaxIter=60, fTolin=10**-2,methodArg=lnpaBases)


#x.performOptimization()
#x.sMethod='single'
#err=x.calcError()[0]
#print err
#nSpecSize=0
#en=x.mloNuclei[1].getEnNu(bAll=True) 
#lHist,lBins,lMu,lSigma,nSize=x.mloNuclei[1].calcEThErr(err, nSpecSize,bPrev=True,bAllME=True)
#nStop=8
#x.mloNuclei[1].plotEthError(lHist[:nStop], lBins[:nStop], lMu[:nStop], lSigma[:nStop], nSize)

#test the fortran fitting method
#sys.path.append('C:\\PythonScripts\\assorted\\')
#import os
#import fortfitinterface
#os.chdir('C:\\PythonScripts\\assorted\\drop\\')
#ans, a, target, npaME,lnpaShortMonoLab,lnpaMonoJLab, shortdiff=x.sMono()
#fakellnAZ=[ [1,1] for i in range(a.shape[0])]
#import numpy as np 
#weight= np.ones(target.shape)
#fortfitinterface.makeDatFile('C:\\PythonScripts\\assorted\\drop',fakellnAZ,target,weight,a)
#os.system('fit')
