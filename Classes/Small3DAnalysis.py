# from ROOT import TFile, gSystem, TStopwatch, TDatime
import ROOT as ro
from optparse import OptionParser
from time import time
from Settings import Settings
from RawEventSaver import RawEventSaver
from PedestalCalculation import PedestalCalculation
# from numpy import array
import numpy as np
from copy import deepcopy
import os, logging, sys

__author__ = 'DA'

class Small3DAnalysis:
    def __init__(self, run=22012, input='./', output='./', settings='./', runinfo='', doPedestal=False, doAlignment=False, doClustering=False, doFiducial=True, doTransparent=False, do3D=False):
        self.run = run
        print 'Running Small3DAnalysis with:'
        string = 'run='+str(self.run)+', '
        string += 'input Dir='+input+', '
        string += 'output Dir='+output+', '
        string += 'settings Dir='+settings+', '
        string += 'do Pedestal, ' if doPedestal else 'no Pedestal, '
        string += 'do Clustering, ' if doClustering else 'no Clustering, '
        string += 'do Fiducial' if doFiducial else 'no Fiducial, '
        string += 'do Alignment, ' if doAlignment else 'no Alignment, '
        string += 'do Transparent' if doTransparent else 'no Transparent, '
        string += 'do 3D Short' if do3D else 'no 3D Short'
        print string
        self.currentDir = ro.gSystem.pwd()
        self.runInfoFile = runinfo if runinfo != '' else 'RunList.ini'
        self.defRunInfo = {'run': int(run), 'run_des': 0, 'verbose': 0, 'nEvents': 100000, 'nStart': 0,
                        'do_pedestal_ana': doPedestal, 'do_cluster_ana': doClustering, 'do_selec_ana': doFiducial,
                        'do_align': doAlignment, 'do_align_ana': False, 'do_trans_ana': doTransparent, 'short_3D_ana': do3D, 'long_3D_ana': False}
        self.runsInfo = []
        self.ReadRunList(self.runInfoFile)
        for runInfo in self.runsInfo:
            self.settings = Settings(runInfo, self.currentDir, input, output, settings)
            self.runWatch = ro.TStopwatch()
            self.runWatch.Start(True)
            time = ro.TDatime()
            logging.basicConfig(filename='analysis_log_{r}_{y}_{m}_{d}_{h}.{min}.{sec}.log'.format(r=self.run, y=time.GetYear(), m=time.GetMonth(), d=time.GetDay(), h=time.GetHour(), min=time.GetMinute(), sec=time.GetSecond()), level=logging.DEBUG)
            # TODO Results class from settings
            # CREATE RAW ROOT TREE
            self.rawEventSaver = RawEventSaver(self.settings)
            self.rawEventSaver.SaveEvents(self.settings.nEvents)  # TODO: self.settings.nEvents
            del self.rawEventSaver
            # PEDESTAL CALCULATION
            self.pedestalCalculation = PedestalCalculation(self.settings)
            self.pedestalCalculation.CalculateSlidingPedestals(runInfo['nEvents'])
        exit()

    def ReadRunList(self, runinfo):
        file = runinfo if runinfo != '' else 'RunList.ini'
        ro.gSystem.cd(self.currentDir)
        if os.path.isfile(file):
            with open(file) as f:
                lines = f.readlines()
            lines = [line.strip() for line in lines]
            for line in list(lines):
                if self.IsAComment(line):
                    lines.remove(line)
            lines = [line.replace('\n', '') for line in lines]
            lines = [line.split('\t') for line in lines]
            for line in lines:
                if len(line) == 11:
                    line.append('0')  # for Short 3D
                    line.append('0')  # for Long 3D
            default = self.defRunInfo
            self.runsInfo = [{'run': int(line[0]),
                             'run_des': int(line[1]),
                             'verbose': int(line[2]),
                             'nEvents': int(line[3]),
                             'nStart': int(line[4]),
                             'do_pedestal_ana': bool(int(line[5])) or default['do_pedestal_ana'],
                             'do_cluster_ana': bool(int(line[6])) or default['do_cluster_ana'],
                             'do_selec_ana': bool(int(line[7])) or default['do_selec_ana'],
                             'do_align': bool(int(line[8])) or default['do_align'],
                             'do_align_ana': bool(int(line[9])) or default['do_align_ana'],
                             'do_trans_ana': bool(int(line[10])) or default['do_trans_ana'],
                             'short_3D_ana': bool(int(line[11])) or default['short_3D_ana'],
                             'long_3D_ana': bool(int(line[12])) or default['long_3D_ana']}
                            for line in lines]
        else:
            self.runsInfo = [self.defRunInfo]

    def IsAComment(self, line):
        return line.startswith('#') or line.startswith('/') or line.startswith('$') or line.startswith('%')

if __name__ == '__main__':
    st = time()
    parser = OptionParser()
    parser.add_option('-r', '--run', dest='run', default=22012, type='int', help='Run to be analysed {e.g. 22012}')
    parser.add_option('-i', '--input', dest='input', default='/User/diegoalejandro/Data/Cern/', type='string', help='Directory path to the raw input files')
    parser.add_option('-o', '--output', dest='output', default='/User/diegoalejandro/Data/Cern/Results', type='string', help='Directory path where the output will be stored')
    parser.add_option('-s', '--settings', dest='settings', default='/User/diegoalejandro/Dropbox/Small3DAnalysis/Results', type='string', help='Directory path where the settings files are')
    parser.add_option('-n', '--runinfo', dest='runinfo', default='', type='string', help='Name of file with the coarse analysis parameters (e.g. RunList.ini)')
    parser.add_option('-p', '--pedestal', action='store_true', dest='doPedestal', default=False, help='enables pedestal analysis')
    parser.add_option('-a', '--align', action='store_true', dest='doAlignment', default=False, help='enables alignment')
    parser.add_option('-c', '--cluster', action='store_true', dest='doClustering', default=False, help='enables clustering')
    parser.add_option('-f', '--fiducial', action='store_true', dest='doFiducial', default=False, help='enables fiducial')
    parser.add_option('-t', '--transparent', action='store_true', dest='doTransparent', default=False, help='enables transparent')
    parser.add_option('-d', '--threeD', action='store_true', dest='do3D', default=False, help='enables 3D short')

    (options, args) = parser.parse_args()
    run = int(options.run)
    input = str(options.input)
    output = str(options.output)
    settings = str(options.settings)
    runinfo = str(options.runinfo)
    doPedestal = bool(options.doPedestal)
    doAlignment = bool(options.doAlignment)
    doClustering = bool(options.doClustering)
    doFiducial = bool(options.doFiducial)
    doTransparent = bool(options.doTransparent)
    do3D = bool(options.do3D)

    z = Small3DAnalysis(run, input, output, settings, runinfo, doPedestal, doAlignment, doClustering, doFiducial, doTransparent, do3D)