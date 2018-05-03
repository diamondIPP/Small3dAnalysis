#!/usr/bin/env python
# from ROOT import TFile, gSystem, TStopwatch, TDatime
import os, logging, sys
sys.path.append('/home/sandiego/Small3dAnalysis/Classes')  # TODO: HARDCODED!!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
import ROOT as ro
from optparse import OptionParser
from time import time
from Settings import Settings
from Converter import Converter
from PedestalCalculation import PedestalCalculation
import numpy as np
from PedestalCalculation2 import PedestalCalculation2
from copy import deepcopy
from Utils import *

__author__ = 'DA'

class SmallRD42Analysis:
    def __init__(self, settings_file='', do_parallel=False):
        print 'Running SmallRD42Analysis...:'
        self.currentDir = ro.gSystem.pwd()
        self.settings = Settings(do_parallel)
        if settings_file == '' or not os.path.isfile(settings_file):
            ExitMessage('Did not finde settings file. Exiting...')
        else:
            self.settings.ReadSettingsFile(settings_file)
        self.runWatch = ro.TStopwatch()
        self.runWatch.Start(True)
        time = ro.TDatime()
        logging.basicConfig(filename='analysis_log_{r}_{y}_{m}_{d}_{h}.{min}.{sec}.log'.format(r=self.settings.run, y=time.GetYear(), m=time.GetMonth(), d=time.GetDay(), h=time.GetHour(), min=time.GetMinute(), sec=time.GetSecond()), level=logging.DEBUG)

        self.converter = None
        # TODO Results class from settings
        # CREATE RAW ROOT TREE
        # self.rawEventSaver = RawEventSaver(self.settings)
        # self.rawEventSaver.SaveEvents(self.settings.)  # TODO: self.settings.nEvents
        # del self.rawEventSaver
        # # PEDESTAL CALCULATION
        # self.pedestalCalculation2 = PedestalCalculation2(self.settings)
        # self.pedestalCalculation2.CalculateSlidingPedestals(runInfo['nEvents'])

    def Run_Analysis(self):
        CreateDirectoryIfNecessary(self.settings.output_dir + '/' + self.settings.sub_dir)
        self.converter = Converter(self.settings)
        if self.converter.do_conversion:
            self.converter.Convert()

        self.rawEventSaver = RawEventSaver(self.settings)
        self.rawEventSaver.SaveEvents()


def main():
    st = time()
    parser = OptionParser()
    parser.add_option('-s', '--settings', dest='settings', default='', type='string', help='settings file')
    parser.add_option('--dont_run', dest='do_run', default=True, action='store_false', help='toggles option for not running the analysis.')
    parser.add_option('--parallelize', dest='parallelize', default=False, action='store_true', help='toggles option for running the analysis in parallel when possible.')

    (options, args) = parser.parse_args()
    settings = str(options.settings)
    do_run = bool(options.do_run)
    do_parallel = bool(options.parallelize)

    rd42 = SmallRD42Analysis(settings_file=settings, do_parallel=do_parallel)
    if do_run:
        rd42.Run_Analysis()


if __name__ == '__main__':
    main()
