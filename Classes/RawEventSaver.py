# from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree
import ROOT as ro
from optparse import OptionParser
from time import time
from Settings import Settings
from RawEventReader import RawEventReader
# from numpy import array, zeros
import numpy as np
from copy import deepcopy
import os, logging

__author__ = 'DA'

class RawEventSaver:
    def __init__(self, settings):
        print 'Starting Raw event saver'
        self.settings = settings
        self.run = np.array(self.settings.run, 'I')
        self.eventNumber = np.zeros(1, 'I')
        self.det_ADC = np.zeros((, self.settings.telDetChs), 'B')  # unsigned char
        self.dia_ADC = np.zeros(self.settings.dutDetChs, 'H')  # unsigned short (int16)
        self.rawEventReader = RawEventReader(self.settings)
        


        self.settings.GoToOutputPath()
        self.rawFilePath = self.settings.GetRawTreeFilePath()
        print 'Raw file path: {p}'.format(p=self.rawFilePath)
        self.createdNewFile = False
        self.rawFile = ro.TFile(self.rawFilePath, 'NEW')
        if self.rawFile.IsOpen():
            self.createdNewFile = True
            print 'Created new Raw file'
        else:
            self.rawFile = ro.TFile(self.rawFilePath, 'READ')
            print 'Opened existing Raw file'
        self.rawTree = self.rawFile.Get('rawTree')
        self.createdNewTree = False
        if not self.rawTree:
            if not self.createdNewFile:
                print 'File had no Tree. Recreating file'
                self.rawFile = ro.TFile(self.rawFilePath, 'RECREATE')
                self.createdNewFile = True
            self.rawTree = ro.TTree('rawTree', 'Raw Data of run {r}'.format(r=self.run))
            self.createdNewTree = True
        # TODO: is this necessary? gSystem.cd('..') it will go up one folder

    def SaveEvents(self, nEvents=0):
        print 'RawEventSaver: Save {n} events'.format(n=nEvents)
        if not self.settings:
            print 'Settings not initialized. Initializing...'
            self.settings = Settings()
        if not self.createdNewTree:
            if self.rawTree.GetEntries() >= nEvents:
                print 'Tree already has enough entries for analysis'
                self.rawFile.Close()
                self.settings.GoToPreviousDir()
                return
            else:
                print 'Tree has not enough entries for analysis. Recreating'
                self.rawFile.Close()
                self.rawFile = ro.TFile(self.rawFilePath, 'RECREATE')
                self.rawTree = ro.TTree('rawTree', 'Raw Data of run {r}'.format(r=self.run))
                self.createdNewTree, self.createdNewFile = True, True
        if self.createdNewTree:
            self.rawTree.Reset()
            self.SetBranches()
            self.settings.CreateProgressBar(nEvents)
            self.settings.bar.start()
            t0 = time()
            for i in xrange(nEvents):
                succeed = self.rawEventReader.ReadRawEvent(i)
                if not succeed:
                    print 'Could not open file: break!'
                    exit()
                self.LoadEvent()
                self.eventNumber = np.array(int(i), 'I')
                self.rawTree.Fill()
                if i == nEvents - 1:
                    self.rawEventReader.current_rz_file.close()
                self.settings.bar.update(i + 1)
            self.rawTree.Write()
            t1 = time()
            self.settings.bar.finish()
            print 'Total time saving raw events:', str(t1-t0), 's'
        self.rawFile.Close()
        self.settings.GoToPreviousDir()

    def TreeExists(self, nEvents=0):
        if (not self.createdNewFile) and (not self.createdNewTree):
            return self.rawTree.GetEntries() >= nEvents
        return False

    def SetBranches(self):
        self.rawTree.Branch('D0X_ADC', self.det_ADC[0], 'D0X_ADC[256]/b')
        self.rawTree.Branch('D0Y_ADC', self.det_ADC[1], 'D0Y_ADC[256]/b')
        self.rawTree.Branch('D1X_ADC', self.det_ADC[2], 'D1X_ADC[256]/b')
        self.rawTree.Branch('D1Y_ADC', self.det_ADC[3], 'D1Y_ADC[256]/b')
        self.rawTree.Branch('D2X_ADC', self.det_ADC[4], 'D2X_ADC[256]/b')
        self.rawTree.Branch('D2Y_ADC', self.det_ADC[5], 'D2Y_ADC[256]/b')
        self.rawTree.Branch('D3X_ADC', self.det_ADC[6], 'D3X_ADC[256]/b')
        self.rawTree.Branch('D3Y_ADC', self.det_ADC[7], 'D3Y_ADC[256]/b')
        self.rawTree.Branch('DiaADC', self.dia_ADC, 'DiaADC[128]/s')
        self.rawTree.Branch('RunNumber', self.run, 'RunNumber/i')
        self.rawTree.Branch('EventNumber', self.eventNumber, 'EventNumber/i')

    def LoadEvent(self):
        diaInput = int(self.settings.diaInput)  # TODO: diaInput in settings class
        for det in xrange(8):  # TODO: 8 is number of detectors in silicon telescope
            # for ch in xrange(self.rawEventReader.silTelChs):
            self.det_ADC[det] = (self.rawEventReader.GetSilPlane(det))
        # self.det_ADC = self.rawEventReader.GetSilPlanes()
        # for ch in xrange(self.rawEventReader.diaChs):
        self.dia_ADC = self.rawEventReader.GetDiaDet(self.settings.diaInput)

if __name__ == '__main__':
    z = RawEventSaver()