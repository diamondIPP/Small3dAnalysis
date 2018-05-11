# from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree, TH1F
import ROOT as ro
from optparse import OptionParser
import os, logging, sys, shutil
sys.path.append('/home/sandiego/Small3dAnalysis/Classes')  # TODO: HARDCODED!!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
from Settings import Settings
from Utils import *
import numpy as np
import pickle
import ipdb
import multiprocessing as mp

__author__ = 'DA'

class PedestalChannelCalculations:
    def __init__(self, settings_bin_path='', channel_num=0, device='dut', show_progressbar=False):
        print 'Creating PedestalCalculations instance'
        self.show_pb = show_progressbar
        self.settings = Settings()
        self.ch = channel_num
        self.device = device
        self.LoadSettingsBinary(settings_bin_path)
        self.out_dir = self.settings.output_dir
        self.sub_dir = self.settings.sub_dir
        self.file_name = self.settings.file_name
        self.tree_name = self.settings.tree_name
        self.slide_leng = self.settings.sliding_length
        self.run = self.settings.run
        self.ped_branches = ['diaPed', 'diaPedSigma', 'cm', 'diaPedCmc', 'diaPedSigmaCmc', 'diaSignal', 'diaSignalCmc']
        self.raw_tel_branches_dic = {0: 'D0X_ADC', 1: 'D0Y_ADC', 2: 'D1X_ADC', 3: 'D1Y_ADC', 4: 'D2X_ADC', 5: 'D2Y_ADC', 6: 'D3X_ADC', 7: 'D3Y_ADC'}
        self.raw_dut_branch = 'DiaADC'
        self.rootFile = None
        self.rootTree = None


        self.read_branch = self.raw_tel_branches_dic[0] if self.device == 'telx0' else self.raw_tel_branches_dic[1] if self.device == 'tely0' else self.raw_tel_branches_dic[2] if self.device == 'telx1' else self.raw_tel_branches_dic[3] if self.device == 'tely1' else self.raw_tel_branches_dic[4] if self.device == 'telx2' else self.raw_tel_branches_dic[5] if self.device == 'tely2' else self.raw_tel_branches_dic[6] if self.device == 'telx3' else self.raw_tel_branches_dic[7] if self.device == 'tely3' else self.raw_dut_branch
        self.read_branch += '[{c}]'.format(c=self.ch)
        self.np_type = self.settings.dut_np_data_type if self.device == 'dut' else self.settings.tel_np_data_type
        self.hit_factor = self.settings.dut_hit_factor if self.device == 'dut' else self.settings.tel_hit_factor
        self.device_ADC = np.zeros(self.slide_leng, dtype=self.np_type)
        self.device_ped = np.zeros(self.slide_leng, dtype='float32')
        self.device_sigma = np.zeros(self.slide_leng, dtype='float32')
        self.device_signal = np.zeros(self.slide_leng, dtype='float32')
        self.device_ADC_mean = []
        self.device_ADC_sigma = []

        self.dut_ADC = np.zeros((self.settings.dutDetChs, self.slide_leng), dtype=self.settings.dut_np_data_type)
        self.tel_ADC = {det: np.zeros((self.settings.telDetChs, self.slide_leng), dtype=self.settings.tel_np_data_type) for det in xrange(self.settings.telDetectors)}
        self.tel_ADC_mean = {det: np.zeros(self.settings.telDetectors, dtype='float32') for det in xrange(self.settings.telDetectors)}
        self.tel_ADC_sigma = {det: np.zeros(self.settings.telDetectors, dtype='float32') for det in xrange(self.settings.telDetectors)}
        self.tel_ADC_snr = {det: np.zeros(self.settings.telDetectors, dtype='float32') for det in xrange(self.settings.telDetectors)}
        self.dut_ADC_mean = np.zeros(self.settings.dutDetChs, dtype='float32')
        self.dut_ADC_sigma = np.zeros(self.settings.dutDetChs, dtype='float32')
        self.dut_ADC_snr = np.zeros(self.settings.dutDetChs, dtype='float32')

        self.CalculatePedestals()

    def LoadSettingsBinary(self, settings_path):
        if os.path.isfile(settings_path):
            with open(settings_path, 'rb') as fs:
                self.settings = pickle.load(fs)
        else:
            ExitMessage('Settings file does not exist!!!!', value=os.EX_OSFILE)

    def CalculatePedestals(self):
        self.CalculateStartingPedestals()

    def CalculateStartingPedestals(self):
        self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), treename=self.tree_name, mode='READ')
        Draw_Branches_For_Get_Val(self.rootTree, [self.read_branch], start_ev=0, n_entries=self.slide_leng, option='goff')
        Get_Branches_Value_To_Numpy(self.rootTree, [self.read_branch], [self.device_ADC], self.slide_leng, 1)

        self.device_ADC_mean.append(self.device_ADC.mean())
        self.device_ADC_sigma.append(self.device_ADC.std())
        self.device_ped.fill(self.device_ADC.mean())
        self.device_sigma.fill(self.device_ADC.std())
        np.putmask(self.device_signal, np.ones(self.slide_leng, '?'), self.device_ADC - self.device_ped)
        for it in xrange(7):
            conditions = abs(self.device_signal) < self.hit_factor * self.device_ADC_sigma[0]
            adc_cond = np.extract(conditions, self.device_ADC)
            self.device_ADC_mean = [adc_cond.mean()]
            self.device_ADC_sigma = [adc_cond.std()]
            self.device_ped.fill(self.device_ADC_mean[0])
            np.putmask(self.device_signal, np.ones(self.slide_leng, '?'), self.device_ADC - self.device_ped)

        self.device_signal.fill(self.device_ADC_sigma[0])
        self.device_ADC_mean = [self.device_ped[0] for i in xrange(self.slide_leng)]
        self.device_ADC_sigma = [self.device_sigma[0] for i in xrange(self.slide_leng)]

def main():
    parser = OptionParser()
    parser.add_option('-s', '--settings', dest='settings', type='string', help='Settings binary file e.g. run_23002_full.settings')
    parser.add_option('-c', '--channel', dest='ch', type='int', default=0, help='channel number')
    parser.add_option('-d', '--device', dest='device', type='string', default='dut', help='Device to analyse (dut, telx0, telx1, telx2, telx3, tely0, tely1, tely2, tely3)')
    parser.add_option('-p', '--progress', dest='progressbar', default=False, action='store_true', help='Shows progress bar of the process')

    (options, args) = parser.parse_args()
    settings_bin_path = str(options.settings)
    channel = int(options.ch)
    device = str(options.device)
    if device not in ['dut', 'telx0', 'telx1', 'telx2', 'telx3', 'tely0', 'tely1', 'tely2', 'tely3']:
        ExitMessage('device must be "dut" or "telx#" or "tely#". Exiting...', os.EX_IOERR)
    do_pb = bool(options.progressbar)
    z = PedestalChannelCalculations(settings_bin_path=settings_bin_path, channel_num=channel, device=device, show_progressbar=do_pb)
    ipdb.set_trace()


if __name__ == '__main__':
    main()
