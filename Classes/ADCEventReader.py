from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree
from numpy import array, zeros
from copy import deepcopy
import os, logging
import ipdb

__author__ = 'DA'

class ADCEventReader:
    def __init__(self, settings=None):
        print 'Starting Pedestal Calculation'
        self.settings = settings
        self.run = array(self.settings.runInfo['run'], dtype='I')
        self.eventNumber = array(int(0), dtype='I')
        self.verb = self.settings.runInfo['verbose']
        self.rawFile, self.rawTree, self.createdNewFile, self.createdNewTree = None, None, False, False
        self.OpenADCTree()
        self.d0XADC, self.d1XADC, self.d2XADC, self.d3XADC = zeros(256, dtype='B'), zeros(256, dtype='B'), zeros(256, dtype='B'), zeros(256, dtype='B')
        self.d0YADC, self.d1YADC, self.d2YADC, self.d3YADC = zeros(256, dtype='B'), zeros(256, dtype='B'), zeros(256, dtype='B'), zeros(256, dtype='B')
        self.diaADC = zeros(128, 'H')
        self.SetBranches()
        self.activeBranches = ['*']
        self.LoadEvent(0)
        self.etaIntegralFile = None
        self.hEtaIntegralSil = {i: None for i in xrange(self.settings.silNumDetectors)}
        self.hEtaIntegralDia = None
        self.LoadEtaDistributions()
        self.eventObject = None

    def OpenADCTree(self):
        temp = self.settings.OpenTree(self.settings.GetRawTreeFilePath(), 'rawTree', False)
        self.rawFile, self.rawTree, self.createdNewFile, self.createdNewTree = temp

    def SetBranches(self):
        if self.rawTree.FindBranch('D0X_ADC'):
            self.rawTree.SetBranchAddress('D0X_ADC', self.d0XADC)
        if self.rawTree.FindBranch('D1X_ADC'):
            self.rawTree.SetBranchAddress('D1X_ADC', self.d1XADC)
        if self.rawTree.FindBranch('D2X_ADC'):
            self.rawTree.SetBranchAddress('D2X_ADC', self.d2XADC)
        if self.rawTree.FindBranch('D3X_ADC'):
            self.rawTree.SetBranchAddress('D3X_ADC', self.d3XADC)
        if self.rawTree.FindBranch('D0Y_ADC'):
            self.rawTree.SetBranchAddress('D0Y_ADC', self.d0YADC)
        if self.rawTree.FindBranch('D1Y_ADC'):
            self.rawTree.SetBranchAddress('D1Y_ADC', self.d1YADC)
        if self.rawTree.FindBranch('D2Y_ADC'):
            self.rawTree.SetBranchAddress('D2Y_ADC', self.d2YADC)
        if self.rawTree.FindBranch('D3Y_ADC'):
            self.rawTree.SetBranchAddress('D3Y_ADC', self.d3YADC)
        if self.rawTree.FindBranch('DiaADC'):
            self.rawTree.SetBranchAddress('DiaADC', self.diaADC)
        if self.rawTree.FindBranch('RunNumber'):
            self.rawTree.SetBranchAddress('RunNumber', self.run)
        if self.rawTree.FindBranch('EventNumber'):
            self.rawTree.SetBranchAddress('EventNumber', self.eventNumber)

    def LoadEvent(self, event=0, branches=['*']):
        # if not self.rawTree:
        #     return False
        # if event <= self.rawTree.GetEntries():
        if self.activeBranches != branches:
            self.rawTree.SetBranchStatus('*', 0)
            self.activeBranches = branches
            for branch in branches:
                self.rawTree.SetBranchStatus(branch, 1)
        self.rawTree.GetEvent(event)
        #     return True
        # return False

    def LoadEtaDistributions(self):
        temp = self.settings.OpenFile(self.settings.GetEtaIntegralFilePath(), True)
        self.etaIntegralFile, createdFile = temp
        if createdFile:
            return
        elif self.etaIntegralFile == -1:
            del self.etaIntegralFile
            print 'File', self.settings.GetEtaIntegralFilePath(), 'is corrupted. Deleting it and creating it again'
            os.remove(self.settings.GetEtaIntegralFilePath())
            self.LoadEtaDistributions()
        self.hEtaIntegralSil = {i: self.etaIntegralFile.Get('hEtaIntegral_{d}'.format(d=i)) for i in xrange(self.settings.silNumDetectors)}
        self.hEtaIntegralDia = self.etaIntegralFile.Get('hEtaIntegral_Dia')
        return

    def GetSilADCValue(self, det, ch):
        if det == 0:
            return self.d0XADC[ch]
        elif det == 1:
            return self.d0YADC[ch]
        elif det == 2:
            return self.d1XADC[ch]
        elif det == 3:
            return self.d1YADC[ch]
        elif det == 4:
            return self.d2XADC[ch]
        elif det == 5:
            return self.d2YADC[ch]
        elif det == 6:
            return self.d3XADC[ch]
        elif det == 7:
            return self.d3YADC[ch]

    def GetDiaADCValue(self, ch):
        return self.diaADC[ch]

if __name__ == '__main__':
    z = ADCEventReader()