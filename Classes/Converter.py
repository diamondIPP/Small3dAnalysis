# from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree
import ROOT as ro
from optparse import OptionParser
from time import time
from Settings import Settings
from RawEventSaver import RawEventSaver
# from numpy import array, zeros
import numpy as np
from copy import deepcopy
import os, logging
import subprocess as subp

__author__ = 'DA'

class Converter:
    def __init__(self, settings):
        print 'Starting Converter...'
        self.settings = settings
        self.tree_name = self.settings.tree_name
        self.file_name = self.settings.file_name
        self.ana_events = self.settings.ana_events
        self.num_parallel = self.settings.num_parallel
        self.out_dir = self.settings.output_dir
        self.ana_dir = self.settings.analysis_path
        self.do_conversion = self.Check_If_Already_Converted()
        self.num_events_per_job = 0
        self.event_saver_processes = []

    def Check_If_Already_Converted(self):
        if not os.path.isfile('{o}/{f}.root'.format(o=self.out_dir, f=self.file_name)):
            print 'The file does not exist. Conversion started.'
            return False
        else:
            tempf = ro.TFile('{o}/{f}.root'.format(o=self.out_dir, f=self.file_name))
            if tempf.IsZombie():
                print 'The existing file cannot be opened. It will be re-written. Conversion started.'
                tempf.Close()
                return False
            tempt = tempf.Get(self.tree_name)
            if not tempt:
                print 'The file does not have the tree', self.tree_name, '. It will be created. Conversion started.'
                tempf.Close()
                return False
            elif tempt.GetEntries() < self.ana_events:
                print 'The file does not have enough events. It will be re-written. Conversion started.'
                tempf.Close()
                return False
            else:
                print 'The file has enough events. Skipping conversion'
                tempf.Close()
                return True

    def Convert(self):
        self.num_events_per_job = int(self.ana_events) / int(self.num_parallel)
        for job_i in xrange(self.num_parallel):
            self.event_saver_processes.append(subp.Popen(['{ad}/RawEventSaver.py', -]))  # TODO save settings file in output dir as pickle and read it. Then use the path as input for the rawsaver. :D

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
    z = Converter()