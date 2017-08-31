from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree, TH1F
from optparse import OptionParser
from time import time
from Settings import Settings
from RawEventReader import RawEventReader
from ADCEventReader import ADCEventReader
from HistogramSaver import HistogramSaver
from numpy import array, zeros, extract, bitwise_and, full
from collections import deque
from copy import deepcopy
import os, logging

__author__ = 'DA'

class PedestalCalculation:
    def __init__(self, settings=None):
        print 'Starting Pedestal Calculation'
        self.settings = settings
        self.run = array(self.settings.runInfo['run'], 'I')
        self.eventNumber = array(0, 'I')
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
        self.detADCValues = [[deque(zeros(self.slidingLength, dtype='f'), self.slidingLength) for ch in xrange(self.settings.silDetChs)] for det in xrange(self.settings.silNumDetectors)]
        self.diaADCValues = [deque(zeros(self.slidingLength, dtype='f'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.diaADCValuesCMC = [deque(zeros(self.slidingLength, dtype='f'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.diaCMNValues = deque(zeros(self.slidingLength, 'f'), self.slidingLength)
        self.settings.GoToOutputPath()
        self.pedestalFilePath = self.settings.GetPedestalTreeFilePath()
        self.pedFile, self.pedTree, self.createdNewFile, self.createdNewTree = None, None, False, False
        self.meanSilValues = zeros((self.settings.silNumDetectors, self.settings.silDetChs), 'f4') # 32 bit float = Float_t
        self.meanSilValuesDeque = [[deque(zeros(self.slidingLength, 'f4'), self.slidingLength) for ch in xrange(self.settings.silDetChs)] for  det in xrange(self.settings.silNumDetectors)]  # 32 bit float = Float_t
        self.sigmaSilValues = zeros((self.settings.silNumDetectors, self.settings.silDetChs), 'f4') # 32 bit float = Float_t
        self.sigmaSilValuesDeque = [[deque(zeros(self.slidingLength, 'f4'), self.slidingLength) for ch in xrange(self.settings.silDetChs)] for det in xrange(self.settings.silNumDetectors)]  # 32 bit float = Float_t
        self.signalSilValues = zeros((self.settings.silNumDetectors, self.settings.silDetChs), 'f4') # 32 bit float = Float_t
        self.signalSilValuesDeque = [[deque(zeros(self.slidingLength, 'f'), self.slidingLength) for ch in xrange(self.settings.silDetChs)] for det in xrange(self.settings.silNumDetectors)]  # 32 bit float = Float_t
        # self.signalSilValuesCMC = zeros((self.settings.silNumDetectors, self.settings.silDetChs), 'f4') # 32 bit float = Float_t
        self.silHitChs = zeros((self.settings.silNumDetectors, self.settings.silDetChs), '?') # boolean
        self.silSeedChs = zeros((self.settings.silNumDetectors, self.settings.silDetChs), '?') # boolean
        self.diaHitChs = zeros(self.settings.diaDetChs, '?') # boolean
        self.diaSeedChs = zeros(self.settings.diaDetChs, '?') # boolean
        self.meanDiaValues = zeros(self.settings.diaDetChs, dtype='f4')
        self.meanDiaValuesDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.sigmaDiaValues = zeros(self.settings.diaDetChs, dtype='f4')
        self.sigmaDiaValuesDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.meanDiaValuesCMC = zeros(self.settings.diaDetChs, dtype='f4')
        self.meanDiaValuesCMCDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.sigmaDiaValuesCMC = zeros(self.settings.diaDetChs, dtype='f4')
        self.sigmaDiaValuesCMCDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.signalDiaValues = zeros(self.settings.diaDetChs, dtype='f4')
        self.signalDiaValuesDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)] # 32 bits float
        self.signalDiaValuesCMC = zeros(self.settings.diaDetChs, dtype='f4')
        self.signalDiaValuesCMCDeque = [deque(zeros(self.slidingLength, dtype='f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)] # 32 bits float
        self.pedestalCalcChs = [deque(zeros(self.slidingLength, 'f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.pedestalCalcChsCMC = [deque(zeros(self.slidingLength, 'f4'), self.slidingLength) for ch in xrange(self.settings.diaDetChs)]
        self.cmn_Dia = array(0, dtype='f4')
        self.doCMC = array(0, dtype='f4')
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
        self.settings.CreateProgressBar(self.settings.nEvents)
        self.settings.bar.start()
        self.FillFirstEvents()
        for ev in xrange(self.slidingLength, self.settings.nEvents):
            branches = ['D0X_ADC', 'D1X_ADC', 'D2X_ADC', 'D3X_ADC', 'D0Y_ADC', 'D1Y_ADC', 'D2Y_ADC', 'D3Y_ADC', 'DiaADC']
            self.eventReader.LoadEvent(ev, branches)
            self.UpdateSiliconPedestals()
            self.DoCMNCalculation()
            self.UpdateDiamondPedestals()
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

        leng = self.eventReader.rawTree.Draw('DiaADC', '', 'goff', self.slidingLength)
        if leng > 1000000:
            self.eventReader.rawTree.SetEstimate(leng)
            leng = self.eventReader.rawTree.Draw('DiaADC', '', 'goff', self.slidingLength)
        tempDia = self.eventReader.rawTree.GetVal(0)
        temp2Dia = array(self.BranchToVector(tempDia, leng, self.settings.diaDetChs), 'f')

        for det in xrange(self.settings.silNumDetectors):
            for ch in xrange(self.settings.silDetChs):
                self.detADCValues[det][ch].extend(temp2Sil[det][:, ch])
                self.meanSilValues[det, ch] = array(self.detADCValues[det][ch], 'f').mean()
                self.sigmaSilValues[det, ch] = array(self.detADCValues[det][ch], 'f').std()
        for ch in xrange(self.settings.diaDetChs):
            self.diaADCValues[ch].extend(temp2Dia[:, ch])
            self.meanDiaValues[ch] = array(self.diaADCValues[ch], 'f').mean()
            self.sigmaDiaValues[ch] = array(self.diaADCValues[ch], 'f').std()

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
            for det in xrange(self.settings.silNumDetectors):
                for ch in xrange(self.settings.silDetChs):
                    self.signalSilValuesDeque[det][ch].extend(self.detADCValues[det][ch] - self.meanSilValues[det, ch])
                    self.meanSilValuesDeque[det][ch].extend(full(self.slidingLength, self.meanSilValues[det][ch], 'f4'))
                    self.sigmaSilValuesDeque[det][ch].extend(full(self.slidingLength, self.sigmaSilValues[det][ch], 'f4'))
            for ch in xrange(self.settings.diaDetChs):
                self.signalDiaValuesDeque[ch] = self.diaADCValues[ch] - self.meanDiaValues[ch]
                self.meanDiaValuesDeque[ch].extend(full(self.slidingLength, self.meanDiaValues[ch], 'f4'))
                self.sigmaDiaValuesDeque[ch].extend(full(self.slidingLength, self.sigmaDiaValues[ch], 'f4'))
        else:  # only calculate for diamond, as silicon does not need cmc.
            for ch in xrange(self.settings.diaDetChs):
                self.signalDiaValuesCMCDeque[ch] = self.diaADCValuesCMC[ch] - self.meanDiaValuesCMC[ch]
                self.meanDiaValuesCMCDeque[ch].extend(full(self.slidingLength, self.meanDiaValuesCMC[ch], 'f4'))
                self.sigmaDiaValuesCMCDeque[ch].extend(full(self.slidingLength, self.sigmaDiaValuesCMC[ch], 'f4'))

    def CalculateFirstPedestalIteration(self, withCMC=False):
        if not withCMC:
            for det in xrange(self.settings.silNumDetectors):
                for ch in xrange(self.settings.silDetChs):
                    self.signalSilValuesDeque[det][ch].extend(self.detADCValues[det][ch] - self.meanSilValues[det, ch])
                    condition = abs(self.signalSilValuesDeque[det][ch]) < self.MaxDetSigma * self.sigmaSilValues[det, ch]
                    self.meanSilValues[det, ch] = extract(condition, array(self.detADCValues[det][ch], 'f')).mean()
                    self.sigmaSilValues[det, ch] = extract(condition, array(self.detADCValues[det][ch], 'f')).std()
                    # self.cmn_sil[det] = zeros(self.settings.silDetChs, 'f')
                    # TODO: INSERT HIT AND SEED FILL HERE
                    self.silHitChs = zeros((self.settings.silNumDetectors, self.settings.silDetChs), '?')
                    self.silSeedChs = zeros((self.settings.silNumDetectors, self.settings.silDetChs), '?')
                    # TODO: INSERT HIT AND SEED FILL BEFORE HERE
            for ch in xrange(self.settings.diaDetChs):
                self.signalDiaValuesDeque[ch] = self.diaADCValues[ch] - self.meanDiaValues[ch]
                condition = abs(self.signalDiaValuesDeque[ch]) < self.MaxDiaSigma * self.sigmaDiaValues[ch]
                self.meanDiaValues[ch] = extract(condition, array(self.diaADCValues[ch], 'f')).mean()
                self.sigmaDiaValues[ch] = extract(condition, array(self.diaADCValues[ch], 'f')).std()
                self.pedestalCalcChs[ch] = array(condition, '?')
                # TODO: INSERT HIT AND SEED FILL HERE
                self.diaHitChs = zeros(self.settings.diaDetChs, '?')
                self.diaSeedChs = zeros(self.settings.diaDetChs, '?')
                # TODO: INSERT HIT AND SEED FILL BEFORE HERE
                # condition_cmn = abs(self.signalDiaValues[ch]/self.sigmaDiaValues[ch]) < self.settings.cmnCut

        else: # only calculate for diamond, as silicon does not need cmc.
            for ch in xrange(self.settings.diaDetChs):
                self.signalDiaValuesCMCDeque[ch] = self.diaADCValuesCMC[ch] - self.meanDiaValuesCMC[ch]
                condition = abs(self.signalDiaValuesCMCDeque[ch]) < self.MaxDiaSigma * self.sigmaDiaValuesCMC[ch]
                self.meanDiaValuesCMC[ch] = extract(condition, array(self.diaADCValuesCMC[ch], 'f')).mean()
                self.sigmaDiaValuesCMC[ch] = extract(condition, array(self.diaADCValuesCMC[ch], 'f')).std()
                self.pedestalCalcChsCMC[ch] = array(condition, '?')

    def CreatePedestalTree(self, events=0):
        temp = self.settings.OpenTree(self.pedestalFilePath, 'pedestalTree', True)
        self.pedFile, self.pedTree, self.createdNewFile, self.createdNewTree = deepcopy(temp)
        self.SetBranches()

    def BranchToVector(self, branch, leng, vecSize=1):
        vector = [[branch[ev*vecSize + ch] for ch in xrange(vecSize)] for ev in xrange(int(leng/vecSize))] if vecSize != 1 else [branch[ev] for ev in xrange(int(leng))]
        return deepcopy(vector)

    def FillFirstEvents(self):
        self.diaCMNValues.clear()
        for ev in xrange(self.slidingLength):
            self.DoCMNCalculation()
            self.diaCMNValues.append(self.cmn_Dia)
            for ch in xrange(self.settings.diaDetChs):
                self.diaADCValuesCMC[ch].append(self.diaADCValues-self.cmn_Dia)
                self.meanDiaValuesCMC[ch] = self.meanDiaValues[ch] - self.cmn_Dia
                self.sigmaDiaValuesCMC[ch] = self.sigmaDiaValues[ch]
        self.CalculateFirstPedestals(True, 7)
        for ev in xrange(self.slidingLength):
            self.eventNumber.fill(ev)
            self.signalSilValues = self.signalSilValuesDeque[:, :, ev]
            self.cmn_Dia = self.diaCMNValues[ev]
            self.signalDiaValues = self.signalDiaValuesDeque[:, ev]
            self.signalDiaValuesCMC = self.signalDiaValuesCMCDeque[:, ev]
            self.pedTree.Fill()
            self.settings.bar.update(ev + 1)

    def DoCMNCalculation(self):
        signalVector = array(self.signalDiaValuesDeque, 'f4')[:, self.eventNumber] if self.eventNumber < self.slidingLength else self.eventReader.diaADC - self.meanDiaValuesCMC
        sigmaVector = self.sigmaDiaValues if self.eventNumber < self.slidingLength else self.sigmaDiaValuesCMC
        snrVector = abs(signalVector / sigmaVector)
        condition1 = snrVector < self.settings.cmnCut
        condition2 = array(1 - self.maskDiaChs, '?')
        condition = bitwise_and(condition1, condition2)
        self.cmn_Dia.fill(extract(condition, signalVector).mean())
        self.hCMN.Fill(self.cmn_Dia)

    def UpdateSiliconPedestals(self):
        for det in xrange(self.settings.silNumDetectors):
            for ch in xrange(self.settings.silDetChs):
                adc = self.eventReader.GetSilADCValue(det, ch)
                self.detADCValues[det][ch].append(adc)
                self.signalSilValuesDeque[det][ch].append(adc - self.meanSilValues[det, ch])
                condition = abs(array(self.signalSilValuesDeque[det][ch], 'f4')) < self.MaxDetSigma * array(self.sigmaSilValuesDeque[det][ch], 'f4')
                self.meanSilValues[det, ch] = extract(condition, array(self.detADCValues[det][ch], 'f4')).mean()
                self.sigmaSilValues[det, ch] = extract(condition, array(self.detADCValues[det][ch], 'f4')).std()

    def UpdateDiamondPedestals(self):
        for ch in xrange(self.settings.diaDetChs):
            adc = self.eventReader.GetDiaADCValue(ch)
            self.diaADCValues[ch].append(adc)
            self.signalDiaValuesDeque[ch].append(adc - self.meanDiaValues[ch])
            condition = abs(array(self.signalDiaValuesDeque[ch], 'f4')) < self.MaxDiaSigma * array(self.sigmaDiaValuesDeque[ch], 'f4')
            self.meanDiaValues[ch] = extract(condition, array(self.diaADCValues[ch], 'f4')).mean()
            self.sigmaDiaValues[ch] = extract(condition, array(self.diaADCValues[ch], 'f4')).std()
            self.diaADCValuesCMC[ch].append(adc - self.cmn_Dia)
            self.signalDiaValuesCMCDeque[ch].append(adc - self.cmn_Dia - self.meanDiaValuesCMC[ch])
            condition = abs(array(self.signalDiaValuesCMCDeque[ch], 'f4')) < self.MaxDiaSigma * array(self.sigmaDiaValuesCMCDeque[ch], 'f4')
            self.meanDiaValuesCMC[ch] = extract(condition, array(self.diaADCValuesCMC[ch], 'f4')).mean()
            self.sigmaDiaValuesCMC[ch] = extract(condition, array(self.diaADCValuesCMC[ch], 'f4')).std()


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

if __name__ == '__main__':
    z = PedestalCalculation()