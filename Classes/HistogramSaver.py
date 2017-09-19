from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree, TPaveText, gStyle
from optparse import OptionParser
from time import time
from Settings import Settings
from RawEventReader import RawEventReader
from ADCEventReader import ADCEventReader
from numpy import array, zeros
from copy import deepcopy
import os, logging
import ipdb

__author__ = 'DA'

class HistogramSaver:
    def __init__(self, settings=None):
        print 'Starting Histogram Saver'
        self.settings = settings
        self.run = array(self.settings.runInfo['run'], 'I')
        self.eventNumber = array(int(0), 'I')
        self.plotPath = '.'
        self.plotRootPath = './root/'
        self.plotPdfPath = './pdf/'
        self.paveTextOptions = {}
        self.dateTime = TDatime()
        self.paveText = TPaveText(0.07,0,0.22,0.10,'NDC')
        self.UpdatePaveText()
        self.verb = self.settings.runInfo['verbose']
        self.optStat1D = 'nemr'
        self.optStat2D = 'ne'
        self.DefaultPlotStyle()

    def DefaultPlotStyle(self):
        gStyle.SetPalette(55)  # 55 is kRainBow. 53 is kDarkBodyRadiator
        gStyle.SetOptStat(self.optStat1D)
        gStyle.SetOptFit(11111)
        gStyle.SetStatH(0.12)
        gStyle.SetStatW(0.15)
        gStyle.SetPadBottomMargin(0.15)
        gStyle.SetPadTopMargin(0.15)

    def SetPath(self, path='.'):
        self.plotPath = self.RemoveExtraBackSlashes(path, 2)
        self.plotRootPath = self.RemoveExtraBackSlashes(path + '/root/', 3)
        self.plotPdfPath = self.RemoveExtraBackSlashes(path + '/pdf/', 3)
        self.settings.CheckDirExistence(self.plotPath, True)
        self.settings.CheckDirExistence(self.plotRootPath, True)
        self.settings.CheckDirExistence(self.plotPdfPath, True)

    def RemoveExtraBackSlashes(self, string, times=1):
        if times == 0:
            return string
        string = string.replace('//', '/')
        return self.RemoveExtraBackSlashes(string, times - 1)


    def UpdatePaveText(self):
        self.paveText.Clear()
        self.paveText.SetTextSize(0.025)
        self.paveTextOptions['svn'] = 'Rev: ' + self.settings.version
        self.paveTextOptions['run'] = 'Run ' + str(self.run)
        self.paveTextOptions['nEvents'] = 'with ' + str(self.settings.runInfo['nEvents']) + 'Events'
        self.paveTextOptions['DateTime'] = self.dateTime.AsSQLString()
        self.paveText.AddText(self.paveTextOptions['svn'])
        self.paveText.AddText(self.paveTextOptions['run'])
        self.paveText.AddText(self.paveTextOptions['nEvents'])
        self.paveText.AddText(self.paveTextOptions['DateTime'])
        self.paveText.SetBorderSize(0)
        self.paveText.SetFillColor(0)


if __name__ == '__main__':
    z = HistogramSaver()