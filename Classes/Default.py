from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree
from optparse import OptionParser
from time import time
from Settings import Settings
from RawEventReader import RawEventReader
from ADCEventReader import ADCEventReader
from numpy import array, zeros
from copy import deepcopy
import os, logging

__author__ = 'DA'

class ADCEventReader:
    def __init__(self, settings=None):
        print 'Starting Pedestal Calculation'
        self.settings = settings
        self.run = array(self.settings.runInfo['run'], 'I')
        self.eventNumber = array(int(0), 'I')
        self.verb = self.settings.runInfo['verbose']

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

if __name__ == '__main__':
    z = ADCEventReader()