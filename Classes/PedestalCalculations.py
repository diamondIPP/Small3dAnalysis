# from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree, TH1F
import ROOT as ro
from optparse import OptionParser
from time import time
from Settings import Settings
from RawEventReader import RawEventReader
from ADCEventReader import ADCEventReader
from HistogramSaver import HistogramSaver
# from numpy import array, zeros, extract, bitwise_and, full, isnan, reshape, repeat, append, empty, ndenumerate
import numpy as np
from collections import deque
from copy import deepcopy
import os, logging
import ipdb
from joblib import Parallel, delayed
import multiprocessing as mp

__author__ = 'DA'

class PedestalCalculations:
    def __init__(self, settings=Settings()):
        print 'Creating PedestalCalculations instance'
        self.settings = settings
        self.slide_leng = self.settings.sliding_length

    def Insert_Value_In_Buffer(self, buff=np.zeros(500), val=0):
        np.putmask(buff, np.bitwise_not(np.zeros(self.slide_leng, '?')), np.roll(buff, -1))
        buff[buff == buff.item(-1)] = val

        self.run = array(self.settings.runInfo['run'], 'I')
        self.eventNumber = array(0, 'I')
        self.timeFinalStart, self.timeFinalEnd = None, None
        self.verb = self.settings.runInfo['verbose']
        self.slidingLength = self.settings.bufferSize
        self.settings.GoToOutputPath()
        self.eventReader = ADCEventReader(self.settings)
        self.histSaver = HistogramSaver(self.settings)
        self.histSaver.SetPath(self.settings.GetAbsoluteOutputPath() + '/pedestalAnalysis/')
        self.MaxDetSigma = self.settings.siPedestalHitFactor
        self.MaxDiaSigma = self.settings.diaPedestalHitFactor
        self.maskDiaChs = array([ch in self.settings.diaChsScreened for ch in xrange(self.settings.diaDetChs)], '?') # boolean
        self.doCMC = array(self.settings.doCMC, dtype='?')
        # self.cmn_sil = zeros(self.settings.silNumDetectors, 'f4')
        self.detADCValues = zeros((self.settings.silNumDetectors, self.settings.silDetChs, self.slidingLength), 'f4')
        self.detADCValuesDeque = [[deque(zeros(self.slidingLength, 'f4'), self.slidingLength) for ch in xrange(self.settings.silDetChs)] for det in xrange(self.settings.silNumDetectors)]
        self.diaADCValues = zeros((self.settings.diaDetChs, self.slidingLength), 'f4')
        self.diaADCValuesCMC = zeros((self.settings.diaDetChs, self.slidingLength), 'f4')
        self.diaCMNValues = zeros(self.slidingLength, 'f')
        self.settings.GoToOutputPath()
        self.pedestalFilePath = self.settings.GetPedestalTreeFilePath()
        self.pedFile, self.pedTree, self.createdNewFile, self.createdNewTree = None, None, False, False
        self.meanSilValues = zeros((self.settings.silNumDetectors, self.settings.silDetChs), 'f4') # 32 bit float = Float_t
        self.meanSilValuesVect = zeros((self.settings.silNumDetectors, self.settings.silDetChs, self.slidingLength), 'f4') # 32 bit float = Float_t
        self.meanSilValuesElements = zeros((self.settings.silNumDetectors, self.settings.silDetChs), 'I')
        self.meanSqSilValues = zeros((settings.silNumDetectors, self.settings.silDetChs), 'f4')
        self.pedestalSilCalcChs = zeros((self.settings.silNumDetectors, self.settings.silDetChs, self.slidingLength), '?')
        self.pedestalSilCalcChsDeque = [[deque(zeros(self.slidingLength, '?'), self.slidingLength) for ch in xrange(self.settings.silDetChs)] for det in xrange(self.settings.silNumDetectors)]
        # self.meanSilValuesDeque = [[deque(zeros(self.slidingLength, 'f4'), self.slidingLength) for ch in xrange(self.settings.silDetChs)] for  det in xrange(self.settings.silNumDetectors)]  # 32 bit float = Float_t
        self.sigmaSilValues = zeros((self.settings.silNumDetectors, self.settings.silDetChs), 'f4') # 32 bit float = Float_t
        self.sigmaSilValuesVect = zeros((self.settings.silNumDetectors, self.settings.silDetChs, self.slidingLength), 'f4') # 32 bit float = Float_t
        # self.sigmaSilValuesDeque = [[deque(zeros(self.slidingLength, 'f4'), self.slidingLength) for ch in xrange(self.settings.silDetChs)] for det in xrange(self.settings.silNumDetectors)]  # 32 bit float = Float_t
        self.signalSilValues = zeros((self.settings.silNumDetectors, self.settings.silDetChs), 'f4') # 32 bit float = Float_t
        self.signalSilValuesVect = zeros((self.settings.silNumDetectors, self.settings.silDetChs, self.slidingLength), 'f4') # 32 bit float = Float_t
        # self.signalSilValuesDeque = [[deque(zeros(self.slidingLength, 'f'), self.slidingLength) for ch in xrange(self.settings.silDetChs)] for det in xrange(self.settings.silNumDetectors)]  # 32 bit float = Float_t
        # self.signalSilValuesCMC = zeros((self.settings.silNumDetectors, self.settings.silDetChs), 'f4') # 32 bit float = Float_t
        self.silHitChs = zeros((self.settings.silNumDetectors, self.settings.silDetChs), '?') # boolean
        self.silSeedChs = zeros((self.settings.silNumDetectors, self.settings.silDetChs), '?') # boolean
        self.diaHitChs = zeros(self.settings.diaDetChs, '?') # boolean
        self.diaSeedChs = zeros(self.settings.diaDetChs, '?') # boolean
        self.meanDiaValues = zeros(self.settings.diaDetChs, dtype='f4')
        self.meanDiaValuesVect = zeros((self.settings.diaDetChs, self.slidingLength), dtype='f4')
        # self.meanDiaValuesDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.sigmaDiaValues = zeros(self.settings.diaDetChs, dtype='f4')
        self.sigmaDiaValuesVect = zeros((self.settings.diaDetChs, self.slidingLength), dtype='f4')
        # self.sigmaDiaValuesDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.meanDiaValuesCMC = zeros(self.settings.diaDetChs, dtype='f4')
        self.meanDiaValuesCMCVect = zeros((self.settings.diaDetChs, self.slidingLength), dtype='f4')
        # self.meanDiaValuesCMCDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.sigmaDiaValuesCMC = zeros(self.settings.diaDetChs, dtype='f4')
        self.sigmaDiaValuesCMCVect = zeros((self.settings.diaDetChs, self.slidingLength), dtype='f4')
        # self.sigmaDiaValuesCMCDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.signalDiaValues = zeros(self.settings.diaDetChs, dtype='f4')
        self.signalDiaValuesVect = zeros((self.settings.diaDetChs, self.slidingLength), dtype='f4')
        # self.signalDiaValuesDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)] # 32 bits float
        self.signalDiaValuesCMC = zeros(self.settings.diaDetChs, dtype='f4')
        self.signalDiaValuesCMCVect = zeros((self.settings.diaDetChs, self.slidingLength), dtype='f4')
        # self.signalDiaValuesCMCDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)] # 32 bits float
        self.pedestalDiaCalcChs = zeros((self.settings.diaDetChs, self.slidingLength), '?')
        # self.pedestalDiaCalcChs = [deque(zeros(self.slidingLength, 'f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.pedestalDiaCalcChsCMC = zeros((self.settings.diaDetChs, self.slidingLength), '?')
        # self.pedestalDiaCalcChsCMC = [deque(zeros(self.slidingLength, 'f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.cmn_Dia = array(0, dtype='f4')
        self.CreatePedestalTree()
        self.StartHistograms()

    def CalculateSlidingPedestals(self, events=0):
        string = 'Calculate Sliding Pedestals with CMN correction' if self.doCMC else 'Calculate Sliding Pedestals without CMN correction'
        print string
        if not self.createdNewTree:
            '{f} has already the pedestal information. No need to calculate pedestals'.format(f=self.pedestalFilePath)
            self.pedFile.Close()
            self.settings.GoToPreviousDir()
            return
        t0 = time()
        self.CalculateStartingPedestal()
        self.CalculateFirstPedestals(False, 7)
        # self.settings.CreateProgressBar(self.settings.nEvents)
        # self.settings.bar.start()
        self.FillFirstEvents()
        for ev in xrange(self.slidingLength, self.settings.nEvents):
            branches = ['D0X_ADC', 'D1X_ADC', 'D2X_ADC', 'D3X_ADC', 'D0Y_ADC', 'D1Y_ADC', 'D2Y_ADC', 'D3Y_ADC', 'DiaADC']
            self.eventNumber.fill(ev)
            self.eventReader.LoadEvent(ev, branches)
            self.CalculateEventPedestal()
            # self.eventReader.LoadEvent(ev, branches)
            # self.UpdateSiliconPedestals()
            # self.DoCMNCalculation()
            # self.UpdateDiamondPedestals()
            self.pedTree.Fill()
            self.settings.bar.update(ev + 1)
        t1 = time()
        self.settings.bar.finish()
        print 'Total time calculating pedestals:', str(t1-t0), 's'
        self.pedTree.AddFriend('rawTree', self.settings.GetRawTreeFilePath())
        self.pedFile.cd()
        self.pedTree.Write()
        self.pedTree.Delete()
        self.pedFile.Close()
        self.settings.GoToPreviousDir()

    def CalculateStartingPedestal(self):
        self.eventReader.rawTree.SetBranchStatus('*', 0)
        for branch in ['D0X_ADC', 'D1X_ADC', 'D2X_ADC', 'D3X_ADC', 'D0Y_ADC', 'D1Y_ADC', 'D2Y_ADC', 'D3Y_ADC']:
            self.eventReader.rawTree.SetBranchStatus(branch, 1)
        leng = self.eventReader.rawTree.Draw('D0X_ADC:D0Y_ADC:D1X_ADC:D1Y_ADC:D2X_ADC:D2Y_ADC:D3X_ADC:D3Y_ADC', '', 'goff para', self.slidingLength)
        if leng > 1000000:
            self.eventReader.rawTree.SetEstimate(leng)
            leng = self.eventReader.rawTree.Draw('D0X_ADC:D0Y_ADC:D1X_ADC:D1Y_ADC:D2X_ADC:D2Y_ADC:D3X_ADC:D3Y_ADC', '', 'goff para', self.slidingLength)
        tempD0X = self.eventReader.rawTree.GetVal(0)
        tempD0Y = self.eventReader.rawTree.GetVal(1)
        tempD1X = self.eventReader.rawTree.GetVal(2)
        tempD1Y = self.eventReader.rawTree.GetVal(3)
        tempD2X = self.eventReader.rawTree.GetVal(4)
        tempD2Y = self.eventReader.rawTree.GetVal(5)
        tempD3X = self.eventReader.rawTree.GetVal(6)
        tempD3Y = self.eventReader.rawTree.GetVal(7)

        temp2Sil = []
        temp2Sil.insert(0, array(self.BranchToVector(tempD0X, leng, self.settings.silDetChs), 'f'))
        temp2Sil.insert(1, array(self.BranchToVector(tempD0Y, leng, self.settings.silDetChs), 'f'))
        temp2Sil.insert(2, array(self.BranchToVector(tempD1X, leng, self.settings.silDetChs), 'f'))
        temp2Sil.insert(3, array(self.BranchToVector(tempD1Y, leng, self.settings.silDetChs), 'f'))
        temp2Sil.insert(4, array(self.BranchToVector(tempD2X, leng, self.settings.silDetChs), 'f'))
        temp2Sil.insert(5, array(self.BranchToVector(tempD2Y, leng, self.settings.silDetChs), 'f'))
        temp2Sil.insert(6, array(self.BranchToVector(tempD3X, leng, self.settings.silDetChs), 'f'))
        temp2Sil.insert(7, array(self.BranchToVector(tempD3Y, leng, self.settings.silDetChs), 'f'))

        self.eventReader.rawTree.SetBranchStatus('*', 0)
        self.eventReader.rawTree.SetBranchStatus('DiaADC', 1)
        leng = self.eventReader.rawTree.Draw('DiaADC', '', 'goff', self.slidingLength)
        if leng > 1000000:
            self.eventReader.rawTree.SetEstimate(leng)
            leng = self.eventReader.rawTree.Draw('DiaADC', '', 'goff', self.slidingLength)
        tempDia = self.eventReader.rawTree.GetVal(0)
        temp2Dia = array(self.BranchToVector(tempDia, leng, self.settings.diaDetChs), 'f')

        for det, ch_value in enumerate(self.detADCValues):
            for ch, value in enumerate(ch_value):
                self.detADCValues[det, ch] = temp2Sil[det][:, ch]
        # for det in xrange(self.settings.silNumDetectors):
        #     for ch in xrange(self.settings.silDetChs):
        #         self.detADCValues[det, ch] = temp2Sil[det][:, ch]
        self.meanSilValues = self.detADCValues.mean(2)
        self.sigmaSilValues = self.detADCValues.std(2)
        self.signalSilValuesVect = self.detADCValues - reshape(repeat(self.meanSilValues, self.slidingLength, axis=1), (self.settings.silNumDetectors, self.settings.silDetChs, self.slidingLength))
        for ch, value in enumerate(self.diaADCValues):
        # for ch in xrange(self.settings.diaDetChs):
            self.diaADCValues[ch] = temp2Dia[:, ch]
        self.meanDiaValues = self.diaADCValues.mean(1)
        self.sigmaDiaValues = self.diaADCValues.std(1)
        self.signalDiaValuesVect = self.diaADCValues - reshape(repeat(self.meanDiaValues, self.slidingLength, axis=0), (self.settings.diaDetChs, self.slidingLength))

    def CalculateFirstPedestals(self, withCMC=False, iter=7):
        if not withCMC:
            print 'Calculating the pedestals for the first {s} entries. Iterating {it} times...'.format(s=self.slidingLength, it=iter)
        else:
            print 'Calculating the pedestals for the first {s} entries with CMC. Iterating {it} times...'.format(s=self.slidingLength, it=iter)
        self.settings.CreateProgressBar(iter)
        self.settings.bar.start()
        t0 = time()
        for it in xrange(iter):
            self.CalculateFirstPedestalIteration(withCMC)
            self.settings.bar.update(it + 1)
        t1 = time()
        self.settings.bar.finish()
        print 'Total time to estimate first pedestals:', str(t1 - t0), 's'
        if not withCMC:
            for det, ch_value in enumerate(self.detADCValues):
                for ch, value in enumerate(ch_value):
                    self.meanSilValuesVect[det, ch].fill(self.meanSilValues[det, ch])
                    self.sigmaSilValuesVect[det, ch].fill(self.sigmaSilValues[det, ch])
                    self.meanSqSilValues[det, ch] = self.sigmaSilValues[det, ch]**2 + self.meanSilValues[det, ch]**2
            for ch, value in enumerate(self.diaADCValues):
                self.meanDiaValuesVect[ch].fill(self.meanDiaValues[ch])
                self.sigmaDiaValuesVect[ch].fill(self.sigmaDiaValues[ch])
        else:
            for ch, value in enumerate(self.diaADCValuesCMC):
                self.meanDiaValuesCMCVect[ch].fill(self.meanDiaValuesCMC[ch])
                self.sigmaDiaValuesCMCVect[ch].fill(self.sigmaDiaValuesCMC[ch])
        # if not withCMC:
        #     self.signalSilValues = self.detADCValues - reshape(repeat(self.meanSilValues, self.slidingLength, axis=1), (self.settings.silNumDetectors, self.settings.silDetChs, self.slidingLength))
            # for det in xrange(self.settings.silNumDetectors):
                # for ch in xrange(self.settings.silDetChs):
                    # self.signalSilValuesDeque[det][ch].extend(self.detADCValues[det][ch] - self.meanSilValues[det, ch])
                    # self.meanSilValuesDeque[det][ch].extend(full(self.slidingLength, self.meanSilValues[det][ch], 'f4'))
                    # self.sigmaSilValuesDeque[det][ch].extend(full(self.slidingLength, self.sigmaSilValues[det][ch], 'f4'))
            # self.signalDiaValues = self.diaADCValues - reshape(repeat(self.meanDiaValues, self.slidingLength, axis=0), (self.settings.diaDetChs, self.slidingLength))
            # for ch in xrange(self.settings.diaDetChs):
            #     self.signalDiaValuesDeque[ch].extend(self.diaADCValues[ch] - self.meanDiaValues[ch])
            #     self.meanDiaValuesDeque[ch].extend(full(self.slidingLength, self.meanDiaValues[ch], 'f4'))
            #     self.sigmaDiaValuesDeque[ch].extend(full(self.slidingLength, self.sigmaDiaValues[ch], 'f4'))
        # else:  # only calculate for diamond, as silicon does not need cmc.
        #     self.signalDiaValuesCMC = self.diaADCValuesCMC - reshape(repeat(self.meanDiaValuesCMC, self.slidingLength, axis=0), (self.settings.diaDetChs, self.slidingLength))
            # for ch in xrange(self.settings.diaDetChs):
            #     self.signalDiaValuesCMCDeque[ch].extend(self.diaADCValuesCMC[ch] - self.meanDiaValuesCMC[ch])
            #     self.meanDiaValuesCMCDeque[ch].extend(full(self.slidingLength, self.meanDiaValuesCMC[ch], 'f4'))
            #     self.sigmaDiaValuesCMCDeque[ch].extend(full(self.slidingLength, self.sigmaDiaValuesCMC[ch], 'f4'))

    def CalculateFirstPedestalIteration(self, withCMC=False):
        if not withCMC:
            # self.signalSilValues = self.detADCValues - reshape(repeat(self.meanSilValues, self.slidingLength, axis=1), (self.settings.silNumDetectors, self.settings.silDetChs, self.slidingLength))
            for det, ch_signals in enumerate(self.signalSilValuesVect):
                for ch, signals in enumerate(ch_signals):
            # for det in xrange(self.settings.silNumDetectors):
            #     for ch in xrange(self.settings.silDetChs):
                    # self.signalSilValuesDeque[det, ch] = self.detADCValues[det, ch] - self.meanSilValues[det, ch]
                    # conditions = abs(self.signalSilValuesDeque[det, ch]) < self.MaxDetSigma * self.sigmaSilValues[det, ch]
                    conditions = abs(signals) < self.MaxDetSigma * self.sigmaSilValues[det, ch]
                    self.pedestalSilCalcChs[det, ch] = conditions
                    self.meanSilValues[det, ch] = extract(conditions, self.detADCValues[det, ch]).mean()
                    self.sigmaSilValues[det, ch] = extract(conditions, self.detADCValues[det, ch]).std()
                    self.meanSilValuesElements[det, ch] = conditions.sum()
                    # self.cmn_sil[det] = zeros(self.settings.silDetChs, 'f')
                    # TODO: INSERT HIT AND SEED FILL HERE
                    self.silHitChs = zeros((self.settings.silNumDetectors, self.settings.silDetChs), '?')
                    self.silSeedChs = zeros((self.settings.silNumDetectors, self.settings.silDetChs), '?')
                    # TODO: INSERT HIT AND SEED FILL BEFORE HERE
            self.signalSilValuesVect = self.detADCValues - reshape(repeat(self.meanSilValues, self.slidingLength, axis=1), (self.settings.silNumDetectors, self.settings.silDetChs, self.slidingLength))
            # self.signalDiaValues = self.diaADCValues - reshape(repeat(self.meanDiaValues, self.slidingLength, axis=0), (self.settings.diaDetChs, self.slidingLength))
            for ch, signals in enumerate(self.signalDiaValuesVect):
                conditions = abs(signals) < self.MaxDiaSigma * self.sigmaDiaValues[ch]
            # for ch in xrange(self.settings.diaDetChs):
                # self.signalDiaValuesDeque[ch] = self.diaADCValues[ch] - self.meanDiaValues[ch]
                # conditions = abs(self.signalDiaValuesDeque[ch]) < self.MaxDiaSigma * self.sigmaDiaValues[ch]
                self.meanDiaValues[ch] = extract(conditions, self.diaADCValues[ch]).mean()
                self.sigmaDiaValues[ch] = extract(conditions, self.diaADCValues[ch]).std()
                self.pedestalDiaCalcChs[ch] = array(conditions, '?')
                # TODO: INSERT HIT AND SEED FILL HERE
                self.diaHitChs = zeros(self.settings.diaDetChs, '?')
                self.diaSeedChs = zeros(self.settings.diaDetChs, '?')
                # TODO: INSERT HIT AND SEED FILL BEFORE HERE
                # condition_cmn = abs(self.signalDiaValues[ch]/self.sigmaDiaValues[ch]) < self.settings.cmnCut
            self.signalDiaValuesVect = self.diaADCValues - reshape(repeat(self.meanDiaValues, self.slidingLength, axis=0), (self.settings.diaDetChs, self.slidingLength))

        else: # only calculate for diamond, as silicon does not need cmc.
            for ch, signals in enumerate(self.signalDiaValuesCMCVect):
            # for ch in xrange(self.settings.diaDetChs):
            #     conditions = abs(self.signalDiaValuesCMC[ch]) < self.MaxDiaSigma * self.sigmaDiaValuesCMC[ch]
                conditions = abs(signals) < self.MaxDiaSigma * self.sigmaDiaValuesCMC[ch]
                self.meanDiaValuesCMC[ch] = extract(conditions, self.diaADCValuesCMC[ch]).mean()
                self.sigmaDiaValuesCMC[ch] = extract(conditions, self.diaADCValuesCMC[ch]).std()
                self.pedestalDiaCalcChsCMC[ch] = array(conditions, '?')
            self.signalDiaValuesCMCVect = self.diaADCValuesCMC - reshape(repeat(self.meanDiaValuesCMC, self.slidingLength, axis=0), (self.settings.diaDetChs, self.slidingLength))

    def CreatePedestalTree(self, events=0):
        temp = self.settings.OpenTree(self.pedestalFilePath, 'pedestalTree', True)
        self.pedFile, self.pedTree, self.createdNewFile, self.createdNewTree = deepcopy(temp)
        self.SetBranches()

    def BranchToVector(self, branch, leng, vecSize=1):
        if int(leng/vecSize == 1):
            vector = [branch[ch] for ch in xrange(vecSize)] if vecSize != 1 else branch[0]
        else:
            vector = [[branch[ev*vecSize + ch] for ch in xrange(vecSize)] for ev in xrange(int(leng/vecSize))] if vecSize != 1 else [branch[ev] for ev in xrange(int(leng))]
        return vector

    def FillFirstEvents(self):
        self.diaCMNValues = empty(0, 'f4')
        for ev in xrange(self.slidingLength):
            self.eventNumber.fill(ev)
            self.DoCMNCalculation()
            self.diaCMNValues = append(self.diaCMNValues, self.cmn_Dia).astype('f4')
            self.diaADCValuesCMC[:, ev] = self.diaADCValues[:, ev] - self.cmn_Dia
            # for ch in xrange(self.settings.diaDetChs):
            #     self.diaADCValuesCMC[ch].append(self.diaADCValues[ch][ev]-self.cmn_Dia)
                # THE FOLLOWING 2 LIENS ARE RIDICULOUS, AS IT WILL OVERWRITE ALWAYS THE VALUE WITHOUT STORING IT AND AT THE END, THE ONLY EVENT THAT MATTERS FOR THIS IS THE LAST ONE (the last cmnvalue)
                # self.meanDiaValuesCMC[ch] = self.meanDiaValues[ch] - self.cmn_Dia
                # self.sigmaDiaValuesCMC[ch] = self.sigmaDiaValues[ch]
        self.meanDiaValuesCMC = self.meanDiaValues - self.cmn_Dia
        self.sigmaDiaValuesCMC = self.sigmaDiaValues
        self.signalDiaValuesCMC = self.diaADCValuesCMC - reshape(repeat(self.meanDiaValuesCMC, self.slidingLength, axis=0), (self.settings.diaDetChs, self.slidingLength))
        self.CalculateFirstPedestals(True, 7)
        self.settings.CreateProgressBar(self.settings.nEvents)
        self.settings.bar.start()
        self.timeFinalStart = time()
        for ev in xrange(self.slidingLength):
            self.eventNumber.fill(ev)
            self.signalSilValues = self.signalSilValuesVect[:, :, ev]
            self.cmn_Dia = self.diaCMNValues[ev]
            self.signalDiaValues = self.signalDiaValuesVect[:, ev]
            self.signalDiaValuesCMC = self.signalDiaValuesCMCVect[:, ev]
            self.pedTree.Fill()
            self.settings.bar.update(ev + 1)
        for det, ch_adc in enumerate(self.detADCValues):
            for ch, adc in enumerate(ch_adc):
                self.pedestalSilCalcChsDeque[det][ch].extend(self.pedestalSilCalcChs[det, ch])
                self.detADCValuesDeque[det][ch].extend(self.detADCValues[det, ch])

    def DoCMNCalculation(self, diaADCs=None):
        signalVector = self.signalDiaValuesVect[:, self.eventNumber] if self.eventNumber < self.slidingLength else diaADCs - self.meanDiaValuesCMC
        sigmaVector = self.sigmaDiaValues if self.eventNumber < self.slidingLength else self.sigmaDiaValuesCMC
        snrVector = abs(signalVector / sigmaVector)
        condition1 = snrVector < self.settings.cmnCut
        condition2 = array(1 - self.maskDiaChs, '?')
        conditions = bitwise_and(condition1, condition2)
        self.cmn_Dia.fill(extract(conditions, signalVector).mean())
        if isnan(self.cmn_Dia):
            self.cmn_Dia.fill(0)
            print 'There were no entries to calculate cmn_Dia in event', self.eventNumber, '. Setting value to 0, as it should be for this case.'
        self.hCMN.Fill(self.cmn_Dia)

    def UpdateSiliconPedestals(self):
        for det, ch_value in enumerate(self.detADCValues):
            for ch, value in enumerate(ch_value):
                adc = self.eventReader.GetSilADCValue(det, ch)
                self.detADCValues[det, ch] = append(self.detADCValues[det, ch], adc)[1:]
                self.signalSilValuesVect[det, ch] = self.detADCValues[det, ch] - self.meanSilValuesVect[det, ch]
                conditions = abs(self.signalSilValuesVect[det, ch]) < self.MaxDetSigma * self.sigmaSilValuesVect[det, ch]
                self.meanSilValues[det, ch] = extract(conditions, self.detADCValues[det, ch]).mean()
                self.sigmaSilValues[det, ch] = extract(conditions, self.detADCValues[det, ch]).std()
                self.meanSilValuesVect[det, ch] = append(self.meanSilValuesVect[det, ch], self.meanSilValues[det, ch])[1:]
                self.sigmaSilValuesVect[det, ch] = append(self.sigmaSilValuesVect[det, ch], self.sigmaSilValues[det, ch])[1:]
        # for det in xrange(self.settings.silNumDetectors):
        #     for ch in xrange(self.settings.silDetChs):
        #         adc = self.eventReader.GetSilADCValue(det, ch)
        #         self.detADCValues[det, ch] = append(self.detADCValues[det, ch], adc)[1:]
        #         # self.detADCValues[det][ch].append(adc)
        #         self.signalSilValuesVect[det, ch] = self.detADCValues[det, ch] - self.meanSilValuesVect[det, ch]
        #         # self.signalSilValuesVect[det, ch] = append(self.signalSilValuesVect[det, ch], adc - self.meanSilValues[det, ch])[1:]
        #         conditions = abs(self.signalSilValuesVect[det, ch]) < self.MaxDetSigma * self.sigmaSilValuesVect[det, ch]
        #         self.meanSilValues[det, ch] = extract(conditions, array(self.detADCValues[det][ch], 'f4')).mean()
        #         self.sigmaSilValues[det, ch] = extract(conditions, array(self.detADCValues[det][ch], 'f4')).std()

    def UpdateDiamondPedestals(self):
        for ch, adc in enumerate(self.eventReader.diaADC):
            self.diaADCValues[ch] = append(self.diaADCValues[ch], adc)[1:]
            self.signalDiaValuesVect[ch] = self.diaADCValues[ch] - self.meanDiaValuesVect[ch]
            conditions = abs(self.signalDiaValuesVect[ch]) < self.MaxDiaSigma * self.sigmaDiaValuesVect[ch]
            self.meanDiaValues[ch] = extract(conditions, self.diaADCValues[ch]).mean()
            self.sigmaDiaValues[ch] = extract(conditions, self.diaADCValues[ch]).std()
            self.meanDiaValuesVect[ch] = append(self.meanDiaValuesVect[ch], self.meanDiaValues[ch])[1:]
            self.sigmaDiaValuesVect[ch] = append(self.sigmaDiaValuesVect[ch], self.sigmaDiaValues[ch])[1:]

            self.diaADCValuesCMC[ch] = append(self.diaADCValuesCMC[ch], adc - self.cmn_Dia)[1:]
            self.signalDiaValuesCMCVect[ch] = self.diaADCValuesCMC[ch] - self.meanDiaValuesCMCVect[ch]
            conditions = abs(self.signalDiaValuesCMCVect[ch]) < self.MaxDiaSigma * self.sigmaDiaValuesCMCVect[ch]
            self.meanDiaValuesCMC[ch] = extract(conditions, self.diaADCValuesCMC[ch]).mean()
            self.sigmaDiaValuesCMC[ch] = extract(conditions, self.diaADCValuesCMC[ch]).std()
            self.meanDiaValuesCMCVect[ch] = append(self.meanDiaValuesCMCVect[ch], self.meanDiaValuesCMC[ch])[1:]
            self.sigmaDiaValuesCMCVect[ch] = append(self.sigmaDiaValuesCMCVect[ch], self.sigmaDiaValuesCMC[ch])[1:]
        # for ch in xrange(self.settings.diaDetChs):
        #     adc = self.eventReader.GetDiaADCValue(ch)
        #     self.diaADCValues[ch].append(adc)
        #     self.signalDiaValuesDeque[ch].append(adc - self.meanDiaValues[ch])
        #     conditions = abs(array(self.signalDiaValuesDeque[ch], 'f4')) < self.MaxDiaSigma * array(self.sigmaDiaValuesDeque[ch], 'f4')
        #     self.meanDiaValues[ch] = extract(conditions, array(self.diaADCValues[ch], 'f4')).mean()
        #     self.sigmaDiaValues[ch] = extract(conditions, array(self.diaADCValues[ch], 'f4')).std()
        #     self.diaADCValuesCMC[ch].append(adc - self.cmn_Dia)
        #     self.signalDiaValuesCMCDeque[ch].append(adc - self.cmn_Dia - self.meanDiaValuesCMC[ch])
        #     conditions = abs(array(self.signalDiaValuesCMCDeque[ch], 'f4')) < self.MaxDiaSigma * array(self.sigmaDiaValuesCMCDeque[ch], 'f4')
        #     self.meanDiaValuesCMC[ch] = extract(conditions, array(self.diaADCValuesCMC[ch], 'f4')).mean()
        #     self.sigmaDiaValuesCMC[ch] = extract(conditions, array(self.diaADCValuesCMC[ch], 'f4')).std()

    # def CalculateBla(self, (detbla, chbla, adcbla)):
    #     self.detADCValues[detbla, chbla] = append(self.detADCValues[detbla, chbla], adcbla)[1:]
    #     self.signalSilValuesVect[detbla, chbla] = self.detADCValues[detbla, chbla] - self.meanSilValuesVect[detbla, chbla]
    #     conditions = abs(self.signalSilValuesVect[detbla, chbla]) < self.MaxDetSigma * self.sigmaSilValuesVect[detbla, chbla]
    #     self.meanSilValues[detbla, chbla] = extract(conditions, self.detADCValues[detbla, chbla]).mean()
    #     self.sigmaSilValues[detbla, chbla] = extract(conditions, self.detADCValues[detbla, chbla]).std()
    #     self.meanSilValuesVect[detbla, chbla] = append(self.meanSilValuesVect[detbla, chbla], self.meanSilValues[detbla, chbla])[1:]
    #     self.sigmaSilValuesVect[detbla, chbla] = append(self.sigmaSilValuesVect[detbla, chbla], self.sigmaSilValues[detbla, chbla])[1:]
    #     return 1

    def CalculateEventPedestal(self):
        # t0 = time()
        # self.eventReader.rawTree.SetBranchStatus('*', 0)
        # for branch in ['D0X_ADC', 'D1X_ADC', 'D2X_ADC', 'D3X_ADC', 'D0Y_ADC', 'D1Y_ADC', 'D2Y_ADC', 'D3Y_ADC']:
        #     self.eventReader.rawTree.SetBranchStatus(branch, 1)
        # leng = self.eventReader.rawTree.Draw('D0X_ADC:D0Y_ADC:D1X_ADC:D1Y_ADC:D2X_ADC:D2Y_ADC:D3X_ADC:D3Y_ADC', '', 'goff para', 1, self.eventNumber)
        # if leng > 1000000:
        #     self.eventReader.rawTree.SetEstimate(leng)
        #     leng = self.eventReader.rawTree.Draw('D0X_ADC:D0Y_ADC:D1X_ADC:D1Y_ADC:D2X_ADC:D2Y_ADC:D3X_ADC:D3Y_ADC', '', 'goff para', 1, self.eventNumber)
        # tempD0X = self.eventReader.rawTree.GetVal(0)
        # tempD0Y = self.eventReader.rawTree.GetVal(1)
        # tempD1X = self.eventReader.rawTree.GetVal(2)
        # tempD1Y = self.eventReader.rawTree.GetVal(3)
        # tempD2X = self.eventReader.rawTree.GetVal(4)
        # tempD2Y = self.eventReader.rawTree.GetVal(5)
        # tempD3X = self.eventReader.rawTree.GetVal(6)
        # tempD3Y = self.eventReader.rawTree.GetVal(7)
        #
        # temp2Sil = []
        # temp2Sil.insert(0, array(self.BranchToVector(tempD0X, leng, self.settings.silDetChs), 'f'))
        # temp2Sil.insert(1, array(self.BranchToVector(tempD0Y, leng, self.settings.silDetChs), 'f'))
        # temp2Sil.insert(2, array(self.BranchToVector(tempD1X, leng, self.settings.silDetChs), 'f'))
        # temp2Sil.insert(3, array(self.BranchToVector(tempD1Y, leng, self.settings.silDetChs), 'f'))
        # temp2Sil.insert(4, array(self.BranchToVector(tempD2X, leng, self.settings.silDetChs), 'f'))
        # temp2Sil.insert(5, array(self.BranchToVector(tempD2Y, leng, self.settings.silDetChs), 'f'))
        # temp2Sil.insert(6, array(self.BranchToVector(tempD3X, leng, self.settings.silDetChs), 'f'))
        # temp2Sil.insert(7, array(self.BranchToVector(tempD3Y, leng, self.settings.silDetChs), 'f'))
        temp2Sil = [self.eventReader.d0XADC, self.eventReader.d0YADC, self.eventReader.d1XADC, self.eventReader.d1YADC,
                    self.eventReader.d2XADC, self.eventReader.d2YADC, self.eventReader.d3XADC, self.eventReader.d3YADC]
        # t1 = time()

        # print '\nReading Silicon Tree took', str(t1-t0), 's'
        #
        # ipdb.set_trace(context=5)
        for det, ch_adc in enumerate(temp2Sil):
            for ch, adc in enumerate(ch_adc):
                # if self.eventNumber == 505 and det == 0 and ch == 176:
                #     print 'det', det, 'ch', ch
                # ipdb.set_trace(context=5)
                self.signalSilValues[det, ch] = adc - self.meanSilValues[det, ch]
                conditions = abs(self.signalSilValues[det, ch]) < self.MaxDetSigma * self.sigmaSilValues[det, ch]
                if self.pedestalSilCalcChsDeque[det][ch][0]:
                    self.meanSilValues[det, ch] = (self.meanSilValues[det, ch] * self.meanSilValuesElements[det, ch] - self.detADCValuesDeque[det][ch][0])/(self.meanSilValuesElements[det, ch] - 1)
                    self.meanSqSilValues[det, ch] = (self.meanSqSilValues[det, ch] * self.meanSilValuesElements[det, ch] - self.detADCValuesDeque[det][ch][0]**2)/(self.meanSilValuesElements[det, ch] - 1)
                    self.meanSilValuesElements[det, ch] -= 1
                self.pedestalSilCalcChsDeque[det][ch].append(conditions)
                self.detADCValuesDeque[det][ch].append(adc)
                if conditions:
                    self.meanSilValues[det, ch] = (self.meanSilValues[det, ch] * self.meanSilValuesElements[det, ch] + adc)/(self.meanSilValuesElements[det, ch] + 1)
                    self.meanSqSilValues[det, ch] = (self.meanSqSilValues[det, ch] * self.meanSilValuesElements[det, ch] + adc**2)/(self.meanSilValuesElements[det, ch] + 1)
                    self.meanSilValuesElements[det, ch] += 1
                self.sigmaSilValues[det, ch] = (self.meanSqSilValues[det, ch] - self.meanSilValues[det, ch]**2)**0.5
                # self.signalSilValuesVect[det, ch] = self.detADCValues[det, ch] - self.meanSilValuesVect[det, ch]
                # conditions = abs(self.signalSilValuesVect[det, ch]) < self.MaxDetSigma * self.sigmaSilValuesVect[det, ch]
                # self.meanSilValues[det, ch] = extract(conditions, self.detADCValues[det, ch]).mean()
                # self.sigmaSilValues[det, ch] = extract(conditions, self.detADCValues[det, ch]).std()
                # self.meanSilValuesVect[det, ch] = append(self.meanSilValuesVect[det, ch], self.meanSilValues[det, ch])[1:]
                # self.sigmaSilValuesVect[det, ch] = append(self.sigmaSilValuesVect[det, ch], self.sigmaSilValues[det, ch])[1:]


        # ipdb.set_trace(context=5)

        self.eventReader.rawTree.SetBranchStatus('*', 0)
        self.eventReader.rawTree.SetBranchStatus('DiaADC', 1)
        leng = self.eventReader.rawTree.Draw('DiaADC', '', 'goff', 1, self.eventNumber)
        if leng > 1000000:
            self.eventReader.rawTree.SetEstimate(leng)
            leng = self.eventReader.rawTree.Draw('DiaADC', '', 'goff', 1, self.eventNumber)
        tempDia = self.eventReader.rawTree.GetVal(0)
        temp2Dia = array(self.BranchToVector(tempDia, leng, self.settings.diaDetChs), 'f')

        # t3 = time()

        # print 'Reading Diamond Tree took', str(t3-t2), 's'
        self.DoCMNCalculation(temp2Dia)

        # t4 = time()

        # print 'Doing CMC calculation took', str(t4-t3), 's'

        for ch, adc in enumerate(temp2Dia):
            # for ch in xrange(self.settings.diaDetChs):
            self.diaADCValues[ch] = append(self.diaADCValues[ch], adc)[1:]
            self.signalDiaValuesVect[ch] = self.diaADCValues[ch] - self.meanDiaValuesVect[ch]
            conditions = abs(self.signalDiaValuesVect[ch]) < self.MaxDiaSigma * self.sigmaDiaValuesVect[ch]
            self.meanDiaValues[ch] = extract(conditions, self.diaADCValues[ch]).mean()
            self.sigmaDiaValues[ch] = extract(conditions, self.diaADCValues[ch]).std()
            self.meanDiaValuesVect[ch] = append(self.meanDiaValuesVect[ch], self.meanDiaValues[ch])[1:]
            self.sigmaDiaValuesVect[ch] = append(self.sigmaDiaValuesVect[ch], self.sigmaDiaValues[ch])[1:]

            self.diaADCValuesCMC[ch] = append(self.diaADCValuesCMC[ch], adc - self.cmn_Dia)[1:]
            self.signalDiaValuesCMCVect[ch] = self.diaADCValuesCMC[ch] - self.meanDiaValuesCMCVect[ch]
            conditions = abs(self.signalDiaValuesCMCVect[ch]) < self.MaxDiaSigma * self.sigmaDiaValuesCMCVect[ch]
            self.meanDiaValuesCMC[ch] = extract(conditions, self.diaADCValuesCMC[ch]).mean()
            self.sigmaDiaValuesCMC[ch] = extract(conditions, self.diaADCValuesCMC[ch]).std()
            self.meanDiaValuesCMCVect[ch] = append(self.meanDiaValuesCMCVect[ch], self.meanDiaValuesCMC[ch])[1:]
            self.sigmaDiaValuesCMCVect[ch] = append(self.sigmaDiaValuesCMCVect[ch], self.sigmaDiaValuesCMC[ch])[1:]
        # t5 = time()

        # print 'Updating Diamond pedestals took', str(t5-t4), 's'

        # print 'Total time per event:', str(t5-t0), 's\n'

    def SetBranches(self):
        self.pedTree.Branch('SilHitChs', self.silHitChs, 'SilHitChs[{ndet}][{chs}]/O'.format(ndet=self.settings.silNumDetectors, chs=self.settings.silDetChs))
        self.pedTree.Branch('SilSeedChs', self.silSeedChs, 'SilSeedChs[{ndet}][{chs}]/O'.format(ndet=self.settings.silNumDetectors, chs=self.settings.silDetChs))
        self.pedTree.Branch('PedestalMeanSil', self.meanSilValues, 'PedestalMeanSil[{ndet}][{chs}]/F'.format(ndet=self.settings.silNumDetectors, chs=self.settings.silDetChs))
        self.pedTree.Branch('PedestalSigmaSil', self.sigmaSilValues, 'PedestalSigmaSil[{ndet}][{chs}]/F'.format(ndet=self.settings.silNumDetectors, chs=self.settings.silDetChs))
        self.pedTree.Branch('SignalSil', self.signalSilValues, 'SignalSil[{ndet}][{chs}]/F'.format(ndet=self.settings.silNumDetectors, chs=self.settings.silDetChs))
        # self.pedTree.Branch('SignalSilCMC', self.signalSilValuesCMC, 'SignalSilCMC[{ndet}][{chs}]/F'.format(ndet=self.settings.silNumDetectors, chs=self.settings.silDetChs))
        # self.pedTree.Branch('CMNSil', self.cmn_sil, 'CMNSil[{ndet}]/F'.format(ndet=self.settings.silNumDetectors))
        self.pedTree.Branch('DiaHitChs', self.diaHitChs, 'DiaHitChs[{chs}]/O'.format(chs=self.settings.diaDetChs))
        self.pedTree.Branch('DiaSeedChs', self.diaSeedChs, 'DiaSeedChs[{chs}]/O'.format(chs=self.settings.diaDetChs))
        self.pedTree.Branch('PedestalMeanDia', self.meanDiaValues, 'PedestalMeanDia[{f}]/F'.format(f=self.settings.diaDetChs))
        self.pedTree.Branch('PedestalSigmaDia', self.sigmaDiaValues, 'PedestalSigmaDia[{f}]/F'.format(f=self.settings.diaDetChs))
        self.pedTree.Branch('CMNDia', self.cmn_Dia, 'CMNDia/F')
        self.pedTree.Branch('PedestalMeanDiaCMC', self.meanDiaValuesCMC, 'PedestalMeanDiaCMC[{f}]/F'.format(f=self.settings.diaDetChs))
        self.pedTree.Branch('PedestalSigmaDiaCMC', self.sigmaDiaValuesCMC, 'PedestalSigmaDiaCMC[{f}]/F'.format(f=self.settings.diaDetChs))
        self.pedTree.Branch('SignalDia', self.signalDiaValues, 'SignalDia[{chs}]/F'.format(chs=self.settings.diaDetChs))
        self.pedTree.Branch('SignalDiaCMC', self.signalDiaValuesCMC, 'SignalDiaCMC[{chs}]/F'.format(chs=self.settings.diaDetChs))
        self.pedTree.Branch('RunNumber', self.run, 'RunNumber/i')
        self.pedTree.Branch('EventNumber', self.eventNumber, 'EventNumber/i')

    def StartHistograms(self):
        self.hCMN = TH1F('CMN', 'CMN', 512, -32, 32)


def main():
    z = PedestalCalculations()


if __name__ == '__main__':
    main()
