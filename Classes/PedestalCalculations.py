# from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree, TH1F
import ROOT as ro
from optparse import OptionParser
import os, logging, sys, shutil
sys.path.append('/home/sandiego/Small3dAnalysis/Classes')  # TODO: HARDCODED!!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
from Settings import Settings
from Utils import *
from PedestalDeviceCalculations2 import PedestalDeviceCalculations2
import time
import numpy as np
import subprocess as subp
# from threading import Thread
from multiprocessing import Manager
import ipdb

__author__ = 'DA'

class PedestalCalculations:
    def __init__(self, settings=Settings()):
        print 'Creating PedestalCalculations instance'
        self.settings = settings
        self.out_dir = self.settings.output_dir
        self.sub_dir = self.settings.sub_dir
        self.file_name = self.settings.file_name
        self.run = self.settings.run
        self.tree_name = self.settings.tree_name
        self.num_parallel = self.settings.num_parallel
        self.slide_leng = self.settings.sliding_length
        self.ana_dir = self.settings.analysis_path
        self.ped_branches = ['diaPed', 'diaPedSigma', 'cm', 'diaPedCmc', 'diaPedSigmaCmc', 'diaSignal', 'diaSignalCmc']
        self.has_ped_branch = {branch: CheckBranchExistence('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), tree_name=self.tree_name, branch=branch) for branch in self.ped_branches}
        self.do_pedestal_calculations = not np.all(self.has_ped_branch.values())
        self.raw_tel_branches = ['D0X_ADC', 'D0Y_ADC', 'D1X_ADC', 'D1Y_ADC', 'D2X_ADC', 'D2Y_ADC', 'D3X_ADC', 'D3Y_ADC']
        self.raw_dut_branch = 'DiaADC'
        self.devices = ['telx0', 'tely0', 'telx1', 'tely1', 'telx2', 'tely2', 'telx3', 'tely3', 'dut']
        self.device_ped_process = {}
        self.rootFile = ro.TFile()
        self.rootTree = ro.TTree()
        self.utils = Utils()

        # self.tel_ADC_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype=self.settings.tel_np_data_type)
        #
        # self.tel_ADC_is_ped_buff = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='?')
        # self.tel_ADC_is_hit_buff = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='?')
        # self.tel_ADC_is_seed_buff = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='?')
        # self.tel_ADC_mean = np.zeros((self.settings.telDetectors, self.settings.telDetChs), dtype='float64')
        # self.tel_ADC_mean_sq = np.zeros((self.settings.telDetectors, self.settings.telDetChs), dtype='float64')
        # self.tel_ADC_sigma = np.zeros((self.settings.telDetectors, self.settings.telDetChs), dtype='float64')
        # self.tel_ADC_mean_elem = np.zeros((self.settings.telDetectors, self.settings.telDetChs), dtype='uint16')
        # self.tel_ADC_is_hit = np.zeros((self.settings.telDetectors, self.settings.telDetChs), dtype='?')
        # self.tel_ADC_is_seed = np.zeros((self.settings.telDetectors, self.settings.telDetChs), dtype='?')
        #
        # self.dut_ADC_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype=self.settings.dut_np_data_type)
        #
        # self.dut_ADC_is_ped_buff = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='?')
        # self.dut_ADC_is_hit_buff = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='?')
        # self.dut_ADC_is_seed_buff = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='?')
        # self.dut_ADC_mean = np.zeros(self.settings.dutDetChs, dtype='float64')
        # self.dut_ADC_mean_sq = np.zeros(self.settings.dutDetChs, dtype='float64')
        # self.dut_ADC_sigma = np.zeros(self.settings.dutDetChs, dtype='float64')
        # self.dut_ADC_mean_elem = np.zeros(self.settings.dutDetChs, dtype='uint16')
        # self.dut_ADC_is_hit = np.zeros(self.settings.dutDetChs, dtype='?')
        # self.dut_ADC_is_seed = np.zeros(self.settings.dutDetChs, dtype='?')
        #
        # tel_Hf = [np.ones((self.settings.telDetChs, self.slide_leng), dtype='uint8') * self.settings.clust_hit[i] for i in xrange(8)]
        # tel_Sf = [np.ones((self.settings.telDetChs, self.slide_leng), dtype='uint8') * self.settings.clust_seed[i] for i in xrange(8)]
        # self.tel_hf_array = np.append(np.append(np.append(np.append(np.append(np.append(np.append(tel_Hf[0], tel_Hf[1]), tel_Hf[2]), tel_Hf[3]), tel_Hf[4]), tel_Hf[5]), tel_Hf[6]), tel_Hf[7]).reshape(self.settings.telDetectors, self.settings.telDetChs, self.slide_leng)
        # self.tel_sf_array = np.append(np.append(np.append(np.append(np.append(np.append(np.append(tel_Sf[0], tel_Sf[1]), tel_Sf[2]), tel_Sf[3]), tel_Sf[4]), tel_Sf[5]), tel_Sf[6]), tel_Sf[7]).reshape(self.settings.telDetectors, self.settings.telDetChs, self.slide_leng)
        self.dic_device_to_pos = {'telx0': 0, 'tely0': 1, 'telx1': 2, 'tely1': 3, 'telx2': 4, 'tely2': 5, 'telx3': 6, 'tely3': 7}

        self.tel_ADC_mean_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='float32')
        self.tel_ADC_sigma_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='float32')
        self.tel_ADC_is_hit_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='?')
        self.tel_ADC_is_seed_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='?')

        self.dut_ADC_mean_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='float32')
        self.dut_ADC_sigma_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='float32')
        self.dut_ADC_is_hit_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='?')
        self.dut_ADC_is_seed_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='?')
        # self.dut_ADC_mean_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='float32')


    # def CalculateDevicesPedestals(self):
    #     jobs_chunks = [self.devices[i:i+self.num_parallel] for i in xrange(0, len(self.devices), self.num_parallel)]
    #     for row in jobs_chunks:
    #         self.device_ped_process = []
    #         for it, device in enumerate(row):
    #             if it == len(row) - 1:
    #                 print 'Calculating pedestals for', device, '. The progress is shown below:'
    #                 self.device_ped_process.append(subp.Popen(['{ad}/PedestalDeviceCalculations.py'.format(ad=self.ana_dir), '-s', '{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), '-d', device, '-p'], bufsize=-1, close_fds=True))
    #             else:
    #                 print 'Calculating pedestals for', device
    #                 self.device_ped_process.append(subp.Popen(['{ad}/PedestalDeviceCalculations.py'.format(ad=self.ana_dir), '-s', '{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), '-d', device], bufsize=-1, close_fds=True))
    #         for job in self.device_ped_process:
    #             while job.poll() is None:
    #                 continue
    #             CloseSubprocess(job)

    def CalculateDevicesPedestals2(self):
        jobs_chunks = [self.devices[i:i+self.num_parallel] for i in xrange(0, len(self.devices), self.num_parallel)]
        for row in jobs_chunks:
            self.device_ped_process = []
            manager = Manager()
            out_dic = manager.dict()
            for it, device in enumerate(row):
                if it == len(row) - 1:
                    print 'Calculating pedestals for', device, '. The progress is shown below:'
                    temp = PedestalDeviceCalculations2('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, True, out_dic)
                    # self.device_ped_process.append()
                    # self.device_ped_process.append(subp.Popen(['{ad}/PedestalDeviceCalculations.py'.format(ad=self.ana_dir), '-s', '{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), '-d', device, '-p'], bufsize=-1, close_fds=True))

                else:
                    print 'Calculating pedestals for', device
                    temp = PedestalDeviceCalculations2('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, False,out_dic)
                    # self.device_ped_process.append(PedestalDeviceCalculations2('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, False))
                    # self.device_ped_process.append(subp.Popen(['{ad}/PedestalDeviceCalculations.py'.format(ad=self.ana_dir), '-s', '{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), '-d', device], bufsize=-1, close_fds=True))
                    # for job in self.device_ped_process:
                    #     while job.poll() is None:
                    #         continue
                    #     CloseSubprocess(job)
                temp.start()
                self.device_ped_process.append(temp)
            for j in self.device_ped_process:
                j.join()
            del self.device_ped_process
            for key_i in out_dic.keys():
                if key_i.startswith('tel'):
                    self.tel_ADC_mean_all[self.dic_device_to_pos[key_i]] = out_dic[key_i]['mean']
                    self.tel_ADC_sigma_all[self.dic_device_to_pos[key_i]] = out_dic[key_i]['sigma']
                    self.tel_ADC_is_hit_all[self.dic_device_to_pos[key_i]] = out_dic[key_i]['is_hit']
                    self.tel_ADC_is_seed_all[self.dic_device_to_pos[key_i]] = out_dic[key_i]['is_seed']
                else:
                    self.dut_ADC_mean_all = out_dic[key_i]['mean']
                    self.dut_ADC_sigma_all = out_dic[key_i]['sigma']
                    self.dut_ADC_is_hit_all = out_dic[key_i]['is_hit']
                    self.dut_ADC_is_seed_all = out_dic[key_i]['is_seed']
            del manager, out_dic
            # while any(thread.is_alive() for thread in self.device_ped_process):
            #     continue
        ipdb.set_trace()



            # def CalculatePedestals(self):
    #     self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), self.tree_name, 'UPDATE')
    #     self.rootTree.SetEstimate(256 * self.settings.ana_events)
    #     time0 = time.time()
    #     print 'Getting all events...', ; sys.stdout.flush()
    #     Draw_Branches_For_Get_Val(self.rootTree, self.raw_tel_branches, 0, self.settings.ana_events, option='goff para')
    #     Get_Branches_Value_To_Numpy(self.rootTree, self.raw_tel_branches, self.tel_ADC_all, self.settings.ana_events, self.settings.telDetChs)
    #     self.rootFile.Close()
    #     self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), self.tree_name, 'UPDATE')
    #     Draw_Branches_For_Get_Val(self.rootTree, [self.raw_dut_branch], 0, self.settings.ana_events, option='goff')
    #     Get_Branches_Value_To_Numpy(self.rootTree, [self.raw_dut_branch], [self.dut_ADC_all], self.settings.ana_events, self.settings.dutDetChs)
    #     self.rootFile.Close()
    #     print 'Done in', time.time() - time0, 'seconds'
    #     self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), self.tree_name, 'UPDATE')
    #     self.SetBranches()
    #     self.CalculateStartingPedestals()
    #
    #     for ev in xrange(self.slide_leng, self.settings.ana_events):
    #         tel_adc_new = self.tel_ADC_all[:, :, ev]
    #         tel_signal_new = np.subtract(tel_adc_new, self.tel_ADC_mean, dtype='float64')
    #         tel_adc_old = self.tel_ADC_all[:, :, ev - self.slide_leng]
    #         tel_condition_p = tel_signal_new < np.multiply(self.tel_hf_array[:, :, -1], self.tel_ADC_sigma, dtype='float64')
    #         tel_condition_h = np.bitwise_and(np.multiply(self.tel_hf_array[:, :, -1], self.tel_ADC_sigma, dtype='float64') <= tel_signal_new, tel_signal_new < np.multiply(self.tel_sf_array[:, :, -1], self.tel_ADC_sigma, dtype='float64'))
    #         tel_condition_s = np.bitwise_and(np.bitwise_not(tel_condition_p), np.bitwise_not(tel_condition_h))
    #         tel_mean1 = self.tel_ADC_mean
    #         tel_mean2 = np.divide(np.subtract(np.multiply(self.tel_ADC_mean, self.tel_ADC_mean_elem, dtype='float64'), tel_adc_old, dtype='float64'), np.subtract(self.tel_ADC_mean_elem, 1.0, dtype='float64'), dtype='float64')
    #         tel_mean1_sq = self.tel_ADC_mean_sq
    #         tel_mean2_sq = np.divide(np.subtract(np.multiply(self.tel_ADC_mean_sq, self.tel_ADC_mean_elem, dtype='float64'), np.power(tel_adc_old, 2.0, dtype='float64'), dtype='float64'), np.subtract(self.tel_ADC_mean_elem, 1.0, dtype='float64'), dtype='float64')
    #         tel_elem1 = self.tel_ADC_mean_elem
    #         tel_elem2 = self.tel_ADC_mean_elem - 1
    #         tel_condition_old = self.tel_ADC_is_ped_buff[:, :, ev - self.slide_leng]
    #         np.putmask(self.tel_ADC_mean, np.ones((self.settings.telDetectors, self.settings.telDetChs), '?'), np.where(tel_condition_old, tel_mean2, tel_mean1))
    #         np.putmask(self.tel_ADC_mean_sq, np.ones((self.settings.telDetectors, self.settings.telDetChs), '?'), np.where(tel_condition_old, tel_mean2_sq, tel_mean1_sq))
    #         self.tel_ADC_mean_elem = np.where(tel_condition_old, tel_elem2, tel_elem1)
    #         tel_mean1 = self.tel_ADC_mean
    #         tel_mean2 = np.divide(np.add(np.multiply(self.tel_ADC_mean, self.tel_ADC_mean_elem, dtype='float64'), tel_adc_new, dtype='float64'), np.add(self.tel_ADC_mean_elem, 1.0, dtype='float64'), dtype='float64')
    #         tel_mean1_sq = self.tel_ADC_mean_sq
    #         tel_mean2_sq = np.divide(np.add(np.multiply(self.tel_ADC_mean_sq, self.tel_ADC_mean_elem, dtype='float64'), np.power(tel_adc_new, 2.0, dtype='float64'), dtype='float64'), np.add(self.tel_ADC_mean_elem, 1.0, dtype='float64'), dtype='float64')
    #         tel_elem1 = self.tel_ADC_mean_elem
    #         tel_elem2 = self.tel_ADC_mean_elem + 1
    #         np.putmask(self.tel_ADC_mean, np.ones((self.settings.telDetectors, self.settings.telDetChs), '?'), np.where(tel_condition_p, tel_mean2, tel_mean1))
    #         np.putmask(self.tel_ADC_mean_sq, np.ones((self.settings.telDetectors, self.settings.telDetChs), '?'), np.where(tel_condition_p, tel_mean2_sq, tel_mean1_sq))
    #         self.tel_ADC_mean_elem = np.where(tel_condition_p, tel_elem2, tel_elem1)
    #         np.putmask(self.tel_ADC_sigma, np.ones((self.settings.telDetectors, self.settings.telDetChs), '?'), np.sqrt(np.subtract(self.tel_ADC_mean_sq, np.power(self.tel_ADC_mean, 2.0, dtype='float64'), dtype='float64'), dtype='float64'))
    #         self.tel_ADC_is_ped_buff[:, :, ev] = tel_condition_p
    #         np.putmask(self.tel_ADC_is_hit, np.ones((self.settings.telDetectors, self.settings.telDetChs), '?'), tel_condition_h)
    #         np.putmask(self.tel_ADC_is_seed, np.ones((self.settings.telDetectors, self.settings.telDetChs), '?'), tel_condition_s)
    #
    #         dut_adc_new = self.dut_ADC_all[:, ev]
    #         dut_signal_new = np.subtract(dut_adc_new, self.dut_ADC_mean, dtype='float64')
    #         dut_adc_old = self.dut_ADC_all[:, ev - self.slide_leng]
    #         dut_condition_p = dut_signal_new < np.multiply(self.settings.clust_hit[8], self.dut_ADC_sigma, dtype='float64')
    #         dut_condition_h = np.bitwise_and(np.multiply(self.settings.clust_hit[8], self.dut_ADC_sigma, dtype='float64') <= dut_signal_new, dut_signal_new < np.multiply(self.settings.clust_seed[8], self.dut_ADC_sigma, dtype='float64'))
    #         dut_condition_s = np.bitwise_and(np.bitwise_not(dut_condition_p), np.bitwise_not(dut_condition_h))
    #         dut_mean1 = self.dut_ADC_mean
    #         dut_mean2 = np.divide(np.subtract(np.multiply(self.dut_ADC_mean, self.dut_ADC_mean_elem, dtype='float64'), dut_adc_old, dtype='float64'), np.subtract(self.dut_ADC_mean_elem, 1.0, dtype='float64'), dtype='float64')
    #         dut_mean1_sq = self.dut_ADC_mean_sq
    #         dut_mean2_sq = np.divide(np.subtract(np.multiply(self.dut_ADC_mean_sq, self.dut_ADC_mean_elem, dtype='float64'), np.power(dut_adc_old, 2.0, dtype='float64'), dtype='float64'), np.subtract(self.dut_ADC_mean_elem, 1.0, dtype='float64'), dtype='float64')
    #         dut_elem1 = self.dut_ADC_mean_elem
    #         dut_elem2 = self.dut_ADC_mean_elem - 1
    #         dut_condition_old = self.dut_ADC_is_ped_buff[:, ev - self.slide_leng]
    #         np.putmask(self.dut_ADC_mean, np.ones(self.settings.dutDetChs, '?'), np.where(dut_condition_old, dut_mean2, dut_mean1))
    #         np.putmask(self.dut_ADC_mean_sq, np.ones(self.settings.dutDetChs, '?'), np.where(dut_condition_old, dut_mean2_sq, dut_mean1_sq))
    #         self.dut_ADC_mean_elem = np.where(dut_condition_old, dut_elem2, dut_elem1)
    #         dut_mean1 = self.dut_ADC_mean
    #         dut_mean2 = np.divide(np.add(np.multiply(self.dut_ADC_mean, self.dut_ADC_mean_elem, dtype='float64'), dut_adc_new, dtype='float64'), np.add(self.dut_ADC_mean_elem, 1.0, dtype='float64'), dtype='float64')
    #         dut_mean1_sq = self.dut_ADC_mean_sq
    #         dut_mean2_sq = np.divide(np.add(np.multiply(self.dut_ADC_mean_sq, self.dut_ADC_mean_elem, dtype='float64'), np.power(dut_adc_new, 2.0, dtype='float64'), dtype='float64'), np.add(self.dut_ADC_mean_elem, 1.0, dtype='float64'), dtype='float64')
    #         dut_elem1 = self.dut_ADC_mean_elem
    #         dut_elem2 = self.dut_ADC_mean_elem + 1
    #         np.putmask(self.dut_ADC_mean, np.ones(self.settings.dutDetChs, '?'), np.where(dut_condition_p, dut_mean2, dut_mean1))
    #         np.putmask(self.dut_ADC_mean_sq, np.ones(self.settings.dutDetChs, '?'), np.where(dut_condition_p, dut_mean2_sq, dut_mean1_sq))
    #         self.dut_ADC_mean_elem = np.where(dut_condition_p, dut_elem2, dut_elem1)
    #         np.putmask(self.dut_ADC_sigma, np.ones(self.settings.dutDetChs, '?'), np.sqrt(np.subtract(self.dut_ADC_mean_sq, np.power(self.dut_ADC_mean, 2.0, dtype='float64'), dtype='float64'), dtype='float64'))
    #         self.dut_ADC_is_ped_buff[:, ev] = dut_condition_p
    #         np.putmask(self.dut_ADC_is_hit, np.ones(self.settings.dutDetChs, '?'), dut_condition_h)
    #         np.putmask(self.dut_ADC_is_seed, np.ones(self.settings.dutDetChs, '?'), dut_condition_s)
    #
    #         self.rootTree.Fill()
    #         self.utils.bar.update(ev + 1)
    #     self.rootTree.Write()
    #     self.utils.bar.finish()
    #     self.rootFile.Close()
    #
    #     ipdb.set_trace()

    def SetBranches(self):
        self.rootTree.Branch('silHitChs', self.tel_ADC_is_hit, 'silHitChs[{ndet}][{chs}]/O'.format(ndet=self.settings.telDetectors, chs=self.settings.telDetChs))
        self.rootTree.Branch('silSeedChs', self.tel_ADC_is_seed, 'silSeedChs[{ndet}][{chs}]/O'.format(ndet=self.settings.telDetectors, chs=self.settings.telDetChs))
        self.rootTree.Branch('silPedestalMeanSil', self.tel_ADC_mean, 'silPedestalMean[{ndet}][{chs}]/F'.format(ndet=self.settings.telDetectors, chs=self.settings.telDetChs))
        self.rootTree.Branch('silPedestalSigmaSil', self.tel_ADC_sigma, 'silPedestalSigma[{ndet}][{chs}]/F'.format(ndet=self.settings.telDetectors, chs=self.settings.telDetChs))
        # self.pedTree.Branch('SignalSil', self.signalSilValues, 'SignalSil[{ndet}][{chs}]/F'.format(ndet=self.settings.silNumDetectors, chs=self.settings.silDetChs))
        # self.pedTree.Branch('SignalSilCMC', self.signalSilValuesCMC, 'SignalSilCMC[{ndet}][{chs}]/F'.format(ndet=self.settings.silNumDetectors, chs=self.settings.silDetChs))
        # self.pedTree.Branch('CMNSil', self.cmn_sil, 'CMNSil[{ndet}]/F'.format(ndet=self.settings.silNumDetectors))
        self.rootTree.Branch('diaHitChs', self.dut_ADC_is_hit, 'diaHitChs[{chs}]/O'.format(chs=self.settings.dutDetChs))
        self.rootTree.Branch('diaSeedChs', self.dut_ADC_is_seed, 'diaSeedChs[{chs}]/O'.format(chs=self.settings.dutDetChs))
        self.rootTree.Branch('diaPedestalMean', self.dut_ADC_mean, 'diaPedestalMean[{f}]/F'.format(f=self.settings.dutDetChs))
        self.rootTree.Branch('diaPedestalSigma', self.dut_ADC_sigma, 'diaPedestalSigma[{f}]/F'.format(f=self.settings.dutDetChs))
        # self.pedTree.Branch('CMNDia', self.cmn_Dia, 'CMNDia/F')
        # self.pedTree.Branch('PedestalMeanDiaCMC', self.meanDiaValuesCMC, 'PedestalMeanDiaCMC[{f}]/F'.format(f=self.settings.diaDetChs))
        # self.pedTree.Branch('PedestalSigmaDiaCMC', self.sigmaDiaValuesCMC, 'PedestalSigmaDiaCMC[{f}]/F'.format(f=self.settings.diaDetChs))
        # self.pedTree.Branch('SignalDia', self.signalDiaValues, 'SignalDia[{chs}]/F'.format(chs=self.settings.diaDetChs))
        # self.pedTree.Branch('SignalDiaCMC', self.signalDiaValuesCMC, 'SignalDiaCMC[{chs}]/F'.format(chs=self.settings.diaDetChs))
        # self.rootTree.Branch('RunNumber', self.run, 'RunNumber/i')
        # self.rootTree.Branch('EventNumber', self.eventNumber, 'EventNumber/i')

    # def CalculateStartingPedestals(self):
    #     temp_tel_ADC = self.tel_ADC_all[:, :, :self.slide_leng]
    #     self.tel_ADC_mean = temp_tel_ADC.mean(axis=2, dtype='float64')
    #     self.tel_ADC_sigma = temp_tel_ADC.std(axis=2, dtype='float64')
    #     temp_tel_mean = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.slide_leng), 'float64')
    #     temp_tel_sigma = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.slide_leng), 'float64')
    #     temp_tel_signal = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.slide_leng), 'float64')
    #     for det in xrange(self.settings.telDetectors):
    #         for ch, val in enumerate(self.tel_ADC_mean[det]):
    #             temp_tel_mean[det, ch].fill(val)
    #         for ch, val in enumerate(self.tel_ADC_sigma[det]):
    #             temp_tel_sigma[det, ch].fill(val)
    #     temp_tel_signal = np.subtract(temp_tel_ADC, temp_tel_mean, dtype='float64')
    #     for it in xrange(7):
    #         condition_p = temp_tel_signal < np.multiply(self.tel_hf_array, temp_tel_sigma, dtype='float64')
    #         condition_h = np.bitwise_and(np.multiply(self.tel_hf_array, temp_tel_sigma, dtype='float64') <= temp_tel_signal, temp_tel_signal < np.multiply(self.tel_sf_array, temp_tel_sigma, dtype='float64'))
    #         condition_s = np.bitwise_and(np.bitwise_not(condition_p), np.bitwise_not(condition_h))
    #         adc_cond = [[np.extract(condition_p[det, ch], temp_tel_ADC) for ch in xrange(self.settings.telDetChs)] for det in xrange(self.settings.telDetectors)]
    #         for det, adc_chs in enumerate(adc_cond):
    #             for ch, adcs in enumerate(adc_chs):
    #                 temp_tel_mean[det, ch].fill(adcs.mean(dtype='float64'))
    #                 temp_tel_sigma[det, ch].fill(adcs.std(dtype='float64'))
    #         temp_tel_signal = np.subtract(temp_tel_ADC, temp_tel_mean, dtype='float64')
    #         self.tel_ADC_is_ped_buff[:, :, :self.slide_leng] = condition_p
    #         self.tel_ADC_is_hit_buff[:, :, :self.slide_leng] = condition_h
    #         self.tel_ADC_is_seed_buff[:, :, :self.slide_leng] = condition_s
    #
    #     self.tel_ADC_mean = temp_tel_mean[:, :, -1]
    #     self.tel_ADC_sigma = temp_tel_sigma[:, :, -1]
    #     self.tel_ADC_mean_sq = np.add(np.power(self.tel_ADC_sigma, 2.0, dtype='float64'), np.power(self.tel_ADC_mean, 2.0, dtype='float64'), dtype='float64')
    #     self.tel_ADC_mean_elem = self.tel_ADC_is_ped_buff[:, :, :self.slide_leng].sum(axis=2)
    #
    #     temp_dut_ADC = self.dut_ADC_all[:, :self.slide_leng]
    #     self.dut_ADC_mean = temp_dut_ADC.mean(axis=1, dtype='float64')
    #     self.dut_ADC_sigma = temp_dut_ADC.std(axis=1, dtype='float64')
    #     temp_dut_mean = np.zeros((self.settings.dutDetChs, self.slide_leng), dtype='float64')
    #     temp_dut_sigma = np.zeros((self.settings.dutDetChs, self.slide_leng), dtype='float64')
    #     temp_dut_signal = np.zeros((self.settings.dutDetChs, self.slide_leng), dtype='float64')
    #     for ch, val in enumerate(self.dut_ADC_mean):
    #         temp_dut_mean[ch].fill(val)
    #     for ch, val in enumerate(self.dut_ADC_sigma):
    #         temp_dut_sigma[ch].fill(val)
    #     temp_dut_signal = np.subtract(temp_dut_ADC, temp_dut_mean, dtype='float64')
    #     for it in xrange(7):
    #         condition_p = temp_dut_signal < np.multiply(self.settings.clust_hit[8], temp_dut_sigma, dtype='float64')
    #         condition_h = np.bitwise_and(np.multiply(self.settings.clust_hit[8], temp_dut_sigma, dtype='float64') <= temp_dut_signal, temp_dut_signal < np.multiply(self.settings.clust_seed[8], temp_dut_sigma, dtype='float64'))
    #         condition_s = np.bitwise_and(np.bitwise_not(condition_p), np.bitwise_not(condition_h))
    #         adc_cond = [np.extract(condition_p[ch], temp_dut_ADC) for ch in xrange(self.settings.dutDetChs)]
    #         for ch, adc_cond_i in enumerate(adc_cond):
    #             temp_dut_mean[ch].fill(adc_cond_i.mean(dtype='float64'))
    #             temp_dut_sigma[ch].fill(adc_cond_i.std(dtype='float64'))
    #         temp_dut_signal = np.subtract(temp_dut_ADC, temp_dut_mean, dtype='float64')
    #         self.dut_ADC_is_ped_buff[:, :self.slide_leng] = condition_p
    #         self.dut_ADC_is_hit_buff[:, :self.slide_leng] = condition_h
    #         self.dut_ADC_is_seed_buff[:, :self.slide_leng] = condition_s
    #
    #     self.dut_ADC_mean = temp_dut_mean[:, -1]
    #     self.dut_ADC_sigma = temp_dut_sigma[:, -1]
    #     self.dut_ADC_mean_sq = np.add(np.power(self.dut_ADC_sigma, 2, dtype='float64'), np.power(self.dut_ADC_mean, 2, dtype='float64'), dtype='float64')
    #     self.dut_ADC_mean_elem = self.dut_ADC_is_ped_buff[:, :self.slide_leng].sum(axis=1)
    #
    #     self.utils.CreateProgressBar(self.settings.ana_events)
    #     self.utils.bar.start()
    #     for ev in xrange(self.slide_leng):
    #         self.rootTree.GetEntry(ev)
    #         np.putmask(self.tel_ADC_is_hit, np.ones((self.settings.telDetectors, self.settings.telDetChs), '?'), self.tel_ADC_is_hit_buff[:, :, ev])
    #         np.putmask(self.tel_ADC_is_seed, np.ones((self.settings.telDetectors, self.settings.telDetChs), '?'), self.tel_ADC_is_seed_buff[:, :, ev])
    #         np.putmask(self.dut_ADC_is_hit, np.ones(self.settings.dutDetChs, '?'), self.dut_ADC_is_hit_buff[:, ev])
    #         np.putmask(self.dut_ADC_is_seed, np.ones(self.settings.dutDetChs, '?'), self.dut_ADC_is_seed_buff[:, ev])
    #         self.rootTree.Fill()
    #         self.utils.bar.update(ev + 1)
        

def main():
    z = PedestalCalculations()


if __name__ == '__main__':
    main()
