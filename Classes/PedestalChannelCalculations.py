#!/usr/bin/env python
import ROOT as ro
from optparse import OptionParser
import os, logging, sys, shutil
sys.path.append('/home/sandiego/Small3dAnalysis/Classes')  # TODO: HARDCODED!!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
from Settings import Settings
from Utils import *
import numpy as np
import pickle
import time
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
        self.ana_events = self.settings.ana_events
        self.ped_branches = ['diaPed', 'diaPedSigma', 'cm', 'diaPedCmc', 'diaPedSigmaCmc', 'diaSignal', 'diaSignalCmc']
        self.raw_tel_branches_dic = {0: 'D0X_ADC', 1: 'D0Y_ADC', 2: 'D1X_ADC', 3: 'D1Y_ADC', 4: 'D2X_ADC', 5: 'D2Y_ADC', 6: 'D3X_ADC', 7: 'D3Y_ADC'}
        self.raw_dut_branch = 'DiaADC'
        self.rootFile = None
        self.rootTree = None
        self.utils = Utils()


        self.read_branch = self.raw_tel_branches_dic[0] if self.device == 'telx0' else self.raw_tel_branches_dic[1] if self.device == 'tely0' else self.raw_tel_branches_dic[2] if self.device == 'telx1' else self.raw_tel_branches_dic[3] if self.device == 'tely1' else self.raw_tel_branches_dic[4] if self.device == 'telx2' else self.raw_tel_branches_dic[5] if self.device == 'tely2' else self.raw_tel_branches_dic[6] if self.device == 'telx3' else self.raw_tel_branches_dic[7] if self.device == 'tely3' else self.raw_dut_branch
        self.hit_factor = self.settings.clust_hit[0] if self.device == 'telx0' else self.settings.clust_hit[1] if self.device == 'tely0' else self.settings.clust_hit[2] if self.device == 'telx1' else self.settings.clust_hit[3] if self.device == 'tely1' else self.settings.clust_hit[4] if self.device == 'telx2' else self.settings.clust_hit[5] if self.device == 'tely2' else self.settings.clust_hit[6] if self.device == 'telx3' else self.settings.clust_hit[7] if self.device == 'tely3' else self.settings.clust_hit[8]
        self.seed_factor = self.settings.clust_seed[0] if self.device == 'telx0' else self.settings.clust_seed[1] if self.device == 'tely0' else self.settings.clust_seed[2] if self.device == 'telx1' else self.settings.clust_seed[3] if self.device == 'tely1' else self.settings.clust_seed[4] if self.device == 'telx2' else self.settings.clust_seed[5] if self.device == 'tely2' else self.settings.clust_seed[6] if self.device == 'telx3' else self.settings.clust_seed[7] if self.device == 'tely3' else self.settings.clust_seed[8]
        self.read_branch += '[{c}]'.format(c=self.ch)
        self.np_type = self.settings.dut_np_data_type if self.device == 'dut' else self.settings.tel_np_data_type
        # self.hit_factor = self.settings.dut_hit_factor if self.device == 'dut' else self.settings.tel_hit_factor
        self.device_ADC = np.zeros(self.slide_leng, dtype=self.np_type)
        self.device_ped = np.zeros(self.slide_leng, dtype='float32')
        self.device_sigma = np.zeros(self.slide_leng, dtype='float32')
        self.device_signal = np.zeros(self.slide_leng, dtype='float32')
        self.device_ADC_mean = np.zeros(self.ana_events, dtype='float32')
        self.device_ADC_sigma = np.zeros(self.ana_events, dtype='float32')
        self.device_ADC_is_ped = np.zeros(self.ana_events, dtype='?')
        self.device_ADC_is_hit = np.zeros(self.ana_events, dtype='?')
        self.device_ADC_is_seed = np.zeros(self.ana_events, dtype='?')

        self.adc = 0
        self.mean = 0.0
        self.sigma = 0.0
        self.mean_sq = 0.0
        self.elem = 0

        self.device_ADC_all = np.zeros(self.ana_events, dtype=self.np_type)

        self.CalculatePedestals()

    def LoadSettingsBinary(self, settings_path):
        if os.path.isfile(settings_path):
            with open(settings_path, 'rb') as fs:
                self.settings = pickle.load(fs)
        else:
            ExitMessage('Settings file does not exist!!!!', value=os.EX_OSFILE)

    # def CalculatePedestals3(self):
    #     self.CalculateStartingPedestals3()
    #     if self.show_pb:
    #         self.utils.CreateProgressBar(maxVal=self.ana_events)
    #         self.utils.bar.start()
    #         self.utils.bar.update(self.slide_leng)
    #     for ev in xrange(self.slide_leng, self.ana_events):
    #         # val = Read_Single_Value_Event(self.rootTree, ev, self.read_branch)
    #         # Insert_Value_In_Buffer_Array(self.device_ADC, val, self.slide_leng)
    #         self.device_ADC = self.device_ADC_all[ev+1-self.slide_leng:ev+1]
    #         self.device_signal = self.device_ADC - self.device_ADC_mean[ev-self.slide_leng:ev]
    #         self.device_sigma = self.device_ADC_sigma[ev-self.slide_leng:ev]
    #         condition_p = self.device_signal < self.hit_factor * self.device_sigma
    #         condition_h = np.bitwise_and(self.hit_factor * self.device_sigma <= self.device_signal, self.device_signal < self.seed_factor * self.device_sigma)
    #         condition_s = np.bitwise_and(np.bitwise_not(condition_p), np.bitwise_not(condition_h))
    #         adc_cond = np.extract(condition_p, self.device_ADC)
    #         self.device_ADC_mean.itemset(ev, adc_cond.mean())
    #         self.device_ADC_sigma.itemset(ev, adc_cond.std())
    #         self.device_ADC_is_ped.itemset(ev, condition_p.item(-1))
    #         self.device_ADC_is_hit.itemset(ev, condition_h.item(-1))
    #         self.device_ADC_is_seed.itemset(ev, condition_s.item(-1))
    #         if self.show_pb:
    #             self.utils.bar.update(ev + 1)
    #     self.utils.bar.finish()
    #
    # def CalculateStartingPedestals3(self):
    #     self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), treename=self.tree_name, mode='READ')
    #     time0 = time.time()
    #     print 'Getting all events...', ; sys.stdout.flush()
    #     Draw_Branches_For_Get_Val(self.rootTree, [self.read_branch], start_ev=0, n_entries=self.ana_events, option='goff')
    #     Get_Branches_Value_To_Numpy(self.rootTree, [self.read_branch], [self.device_ADC_all], self.ana_events, 1)
    #     print 'Done in', time.time() - time0, 'seconds'
    #     self.device_ADC = self.device_ADC_all[:self.slide_leng]
    #
    #     self.device_ADC_mean[:self.slide_leng] = self.device_ADC.mean()
    #     self.device_ADC_sigma[:self.slide_leng] = self.device_ADC.std()
    #     self.device_ped.fill(self.device_ADC_mean[0])
    #     self.device_sigma.fill(self.device_ADC_sigma[0])
    #     self.device_signal = self.device_ADC - self.device_ped
    #     for it in xrange(7):
    #         # ipdb.set_trace(context=7)
    #         condition_p = self.device_signal < self.hit_factor * self.device_ADC_sigma[-1]
    #         condition_h = np.bitwise_and(self.hit_factor * self.device_sigma <= self.device_signal, self.device_signal < self.seed_factor * self.device_sigma)
    #         condition_s = np.bitwise_and(np.bitwise_not(condition_p), np.bitwise_not(condition_h))
    #         adc_cond = np.extract(condition_p, self.device_ADC)
    #         self.device_ADC_mean.fill(adc_cond.mean())
    #         self.device_ADC_sigma.fill(adc_cond.std())
    #         self.device_ped.fill(adc_cond.mean())
    #         self.device_sigma.fill(adc_cond.std())
    #         self.device_signal = self.device_ADC - self.device_ped
    #         self.device_ADC_is_ped[:self.slide_leng] = condition_p
    #         self.device_ADC_is_hit[:self.slide_leng] = condition_h
    #         self.device_ADC_is_seed[:self.slide_leng] = condition_s

    def CalculatePedestals(self):
        self.CalculateStartingPedestals()
        if self.show_pb:
            self.utils.CreateProgressBar(maxVal=self.ana_events)
            self.utils.bar.start()
            self.utils.bar.update(self.slide_leng)
        for ev in xrange(self.slide_leng, self.ana_events):
            adc_new = self.device_ADC_all.item(ev)
            adc_old = self.device_ADC_all.item(ev - self.slide_leng)
            condition_p = ((adc_new - self.mean) < self.hit_factor * self.sigma)
            condition_h = (self.hit_factor * self.sigma <= (adc_new - self.mean)) and ((adc_new - self.mean) < self.seed_factor * self.sigma)
            condition_s = (not condition_p) and (not condition_h)

            if self.device_ADC_is_ped.item(ev - self.slide_leng):
                self.mean = (self.mean * self.elem - adc_old) / (self.elem - 1.0)
                self.mean_sq = (self.mean_sq * self.elem - adc_old ** 2) / (self.elem - 1.0)
                self.elem -= 1

            if condition_p:
                self.mean = (self.mean * self.elem + adc_new) / (self.elem + 1.0)
                self.mean_sq = (self.mean_sq * self.elem + adc_new ** 2) / (self.elem + 1.0)
                self.elem += 1

            self.sigma = (self.mean_sq - self.mean ** 2) ** 0.5

            self.device_ADC_mean.itemset(ev, self.mean)
            self.device_ADC_sigma.itemset(ev, self.sigma)
            self.device_ADC_is_ped.itemset(ev, condition_p)
            self.device_ADC_is_hit.itemset(ev, condition_h)
            self.device_ADC_is_seed.itemset(ev, condition_s)
            if self.show_pb:
                self.utils.bar.update(ev + 1)
        self.utils.bar.finish()

    def CalculateStartingPedestals(self):
        self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), treename=self.tree_name, mode='READ')
        time0 = time.time()
        print 'Getting all events...', ; sys.stdout.flush()
        Draw_Branches_For_Get_Val(self.rootTree, [self.read_branch], start_ev=0, n_entries=self.ana_events, option='goff')
        Get_Branches_Value_To_Numpy(self.rootTree, [self.read_branch], [self.device_ADC_all], self.ana_events, 1)
        print 'Done in', time.time() - time0, 'seconds'
        self.device_ADC = self.device_ADC_all[:self.slide_leng]

        self.device_ADC_mean[:self.slide_leng] = self.device_ADC.mean()
        self.device_ADC_sigma[:self.slide_leng] = self.device_ADC.std()
        self.device_ped.fill(self.device_ADC_mean[0])
        self.device_sigma.fill(self.device_ADC_sigma[0])
        self.device_signal = self.device_ADC - self.device_ped
        for it in xrange(7):
            # ipdb.set_trace(context=7)
            condition_p = self.device_signal < self.hit_factor * self.device_ADC_sigma[-1]
            condition_h = np.bitwise_and(self.hit_factor * self.device_sigma <= self.device_signal, self.device_signal < self.seed_factor * self.device_sigma)
            condition_s = np.bitwise_and(np.bitwise_not(condition_p), np.bitwise_not(condition_h))
            adc_cond = np.extract(condition_p, self.device_ADC)
            self.device_ADC_mean.fill(adc_cond.mean())
            self.device_ADC_sigma.fill(adc_cond.std())
            self.device_ped.fill(adc_cond.mean())
            self.device_sigma.fill(adc_cond.std())
            self.device_signal = self.device_ADC - self.device_ped
            self.device_ADC_is_ped[:self.slide_leng] = condition_p
            self.device_ADC_is_hit[:self.slide_leng] = condition_h
            self.device_ADC_is_seed[:self.slide_leng] = condition_s
        self.mean = self.device_ADC_mean[0]
        self.sigma = self.device_ADC_sigma[0]
        self.mean_sq = self.sigma ** 2 + self.mean ** 2
        self.elem = self.device_ADC_is_ped[:self.slide_leng].sum()


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
