#!/usr/bin/env python
import ROOT as ro
import os, sys, logging
sys.path.append('/home/sandiego/Small3dAnalysis/Classes')
from optparse import OptionParser
from Settings import Settings
from RawEventReader import RawEventReader
from Utils import *
import numpy as np
import pickle

__author__ = 'DA'

class RawEventSaver:
    def __init__(self, settings_bin_path='', job_num=0, show_progressbar=False):
        print 'Starting Raw event saver'
        self.show_pb = show_progressbar
        self.settings = Settings()
        self.job = job_num
        self.LoadSettingsBinary(settings_bin_path)
        self.run = self.settings.run
        self.out_dir = self.settings.output_dir
        self.sub_dir = self.settings.sub_dir
        self.file_name = self.settings.file_name
        self.tree_name = self.settings.tree_name

        self.eventNumber = np.zeros(1, 'I')
        self.tel_ADC = np.zeros((self.settings.telDetectors, self.settings.telDetChs), 'B')  # unsigned char
        self.dut_ADC = np.zeros(self.settings.dutDetChs, 'H')  # unsigned short (int16)
        self.dut_chs = np.arange(0, self.settings.dutDetChs, dtype='B')  # unsigned char

        self.rawEventReader = RawEventReader(self.settings)
        self.rawFile = ro.TFile('{o}/{s}/{r}/{f}_{j}.root'.format(o=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name, j=self.job), 'RECREATE')
        self.rawTree = ro.TTree(self.tree_name, self.tree_name)
        self.utils = Utils()

        self.num_events_to_save_firsts = int(self.settings.ana_events) / int(self.settings.num_parallel)
        self.num_events_to_save_last = self.settings.ana_events - self.job * self.num_events_to_save_firsts
        self.num_events_to_save = self.num_events_to_save_firsts if self.job != self.settings.num_parallel - 1 else self.num_events_to_save_last
        self.ev_ini = self.settings.first_event + self.job * self.num_events_to_save_firsts
        self.ev_fin = self.settings.first_event + self.job * self.num_events_to_save_firsts + self.num_events_to_save  # not inclusive, as range or xrange do not take the last one

        self.SaveEvents()

    def LoadSettingsBinary(self, settings_path):
        if os.path.isfile(settings_path):
            with open(settings_path, 'rb') as fs:
                self.settings = pickle.load(fs)
        else:
            ExitMessage('Settings file does not exist!!!!', value=os.EX_OSFILE)

    def SaveEvents(self):
        self.rawTree.Reset()
        self.SetBranches()
        if self.show_pb:
            print 'Showing last job, saving from event {evi} until event {evf}:'.format(evi=self.ev_ini, evf=self.ev_fin-1)
            self.utils.CreateProgressBar(self.num_events_to_save)
            self.utils.bar.start()
            i = 0
        for eve in xrange(self.ev_ini, self.ev_fin):
            succeed = self.rawEventReader.ReadRawEvent(eve)
            if not succeed:
                ExitMessage('Could not open the raw file to save event: {ev}!. Exiting...'.format(ev=eve), os.EX_OSFILE)
            self.LoadEvent(eve)
            self.rawTree.Fill()
            if eve == self.ev_fin - 1:
                self.rawEventReader.current_rz_file.close()
            if self.show_pb:
                self.utils.bar.update(i + 1)
                i += 1
        self.CloseAll()

    def SetBranches(self):
        self.rawTree.Branch('D0X_ADC', self.tel_ADC[0], 'D0X_ADC[{chs}]/b'.format(chs=self.settings.telDetChs))
        self.rawTree.Branch('D0Y_ADC', self.tel_ADC[1], 'D0Y_ADC[{chs}]/b'.format(chs=self.settings.telDetChs))
        self.rawTree.Branch('D1X_ADC', self.tel_ADC[2], 'D1X_ADC[{chs}]/b'.format(chs=self.settings.telDetChs))
        self.rawTree.Branch('D1Y_ADC', self.tel_ADC[3], 'D1Y_ADC[{chs}]/b'.format(chs=self.settings.telDetChs))
        self.rawTree.Branch('D2X_ADC', self.tel_ADC[4], 'D2X_ADC[{chs}]/b'.format(chs=self.settings.telDetChs))
        self.rawTree.Branch('D2Y_ADC', self.tel_ADC[5], 'D2Y_ADC[{chs}]/b'.format(chs=self.settings.telDetChs))
        self.rawTree.Branch('D3X_ADC', self.tel_ADC[6], 'D3X_ADC[{chs}]/b'.format(chs=self.settings.telDetChs))
        self.rawTree.Branch('D3Y_ADC', self.tel_ADC[7], 'D3Y_ADC[{chs}]/b'.format(chs=self.settings.telDetChs))
        self.rawTree.Branch('DiaADC', self.dut_ADC, 'DiaADC[{chs}]/s'.format(chs=self.settings.dutDetChs))
        self.rawTree.Branch('diaChannel', self.dut_chs, 'diaChannel[{chs}]/b'.format(chs=self.settings.dutDetChs))
        self.rawTree.Branch('event', self.eventNumber, 'event/i')

    def LoadEvent(self, ev):
        self.eventNumber.fill(ev)
        for det in xrange(self.settings.telDetectors):
            np.putmask(self.tel_ADC[det], np.bitwise_not(np.zeros(self.settings.telDetChs, '?')), self.rawEventReader.GetTelPlane(det))
        np.putmask(self.dut_ADC, np.bitwise_not(np.zeros(self.settings.dutDetChs, '?')), self.rawEventReader.GetDutDet(self.settings.dut_input))

    def CloseAll(self):
        outmessage = ''
        if self.show_pb:
            self.utils.bar.finish()
            outmessage = 'Finished converting from event {evi} until event {evf}'.format(evi=self.ev_ini, evf=self.ev_fin-1)
        self.rawFile.Write()
        self.rawFile.Close()
        ExitMessage(outmessage, os.EX_OK)


def main():
    parser = OptionParser()
    parser.add_option('-s', '--settings', dest='settings', type='string', help='Settings binary file e.g. run_23002_full.settings')
    parser.add_option('-j', '--job', dest='job', type='int', default=0, help='job number')
    parser.add_option('-p', '--progress', dest='progressbar', default=False, action='store_true', help='Shows progress bar of the process')

    (options, args) = parser.parse_args()
    settings_bin_path = str(options.settings)
    job_num = int(options.job)
    do_pb = bool(options.progressbar)

    z = RawEventSaver(settings_bin_path=settings_bin_path, job_num=job_num, show_progressbar=do_pb)


if __name__ == '__main__':
    main()
