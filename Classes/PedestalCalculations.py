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
# from multiprocessing import Manager, sharedctypes
import multiprocessing as mp
import multiprocessing.sharedctypes as mpsc
import ctypes
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

        tel_ADC_mean_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='float32')
        tel_ADC_sigma_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='float32')
        tel_ADC_is_hit_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='uint8')
        tel_ADC_is_seed_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='uint8')

        dut_ADC_mean_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='float32')
        dut_ADC_sigma_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='float32')
        dut_ADC_is_hit_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='uint8')
        dut_ADC_is_seed_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='uint8')

        self.tel_ADC_mean_all_mp = mp.Array(np.ctypeslib.as_ctypes(tel_ADC_mean_all)._type_, np.ctypeslib.as_ctypes(tel_ADC_mean_all))
        self.tel_ADC_sigma_all_mp = mp.Array(np.ctypeslib.as_ctypes(tel_ADC_sigma_all)._type_, np.ctypeslib.as_ctypes(tel_ADC_sigma_all))
        self.tel_ADC_is_hit_all_mp = mp.Array(np.ctypeslib.as_ctypes(tel_ADC_is_hit_all)._type_, np.ctypeslib.as_ctypes(tel_ADC_is_hit_all))
        self.tel_ADC_is_seed_all_mp = mp.Array(np.ctypeslib.as_ctypes(tel_ADC_is_seed_all)._type_, np.ctypeslib.as_ctypes(tel_ADC_is_seed_all))

        self.dut_ADC_mean_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_mean_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_mean_all))
        self.dut_ADC_sigma_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_sigma_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_sigma_all))
        self.dut_ADC_is_hit_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_hit_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_hit_all))
        self.dut_ADC_is_seed_all_mp = mp.Array(np.ctypeslib.as_ctypes(dut_ADC_is_hit_all)._type_, np.ctypeslib.as_ctypes(dut_ADC_is_seed_all))
        del tel_ADC_mean_all, tel_ADC_sigma_all, tel_ADC_is_hit_all, tel_ADC_is_seed_all, dut_ADC_mean_all, dut_ADC_sigma_all, dut_ADC_is_hit_all, dut_ADC_is_seed_all

        self.tel_ADC_all_mp = self.LoadTelescopeADCs()
        self.dut_ADC_all_mp = self.LoadDutADCs()

    def LoadTelescopeADCs(self):
        self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), treename=self.tree_name, mode='READ')
        self.rootTree.SetEstimate(256 * self.rootTree.GetEntries())
        time0 = time.time()
        print 'Getting all events for telescope...', ; sys.stdout.flush()
        Draw_Branches_For_Get_Val(self.rootTree, self.raw_tel_branches, 0, self.settings.ana_events, option='goff para')
        temp = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype=self.settings.tel_np_data_type)
        Get_Branches_Value_To_Numpy(self.rootTree, self.raw_tel_branches, temp, self.settings.ana_events, self.settings.telDetChs)
        self.rootFile.Close()
        print 'Done in', time.time() - time0, 'seconds'
        return mp.Array(np.ctypeslib.as_ctypes(temp)._type_, np.ctypeslib.as_ctypes(temp))

    def LoadDutADCs(self):
        self.rootFile, self.rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), treename=self.tree_name, mode='READ')
        self.rootTree.SetEstimate(256 * self.rootTree.GetEntries())
        time0 = time.time()
        print 'Getting all events for DUT...', ;
        sys.stdout.flush()
        Draw_Branches_For_Get_Val(self.rootTree, [self.raw_dut_branch], 0, self.settings.ana_events, option='goff')
        temp = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype=self.settings.dut_np_data_type)
        Get_Branches_Value_To_Numpy(self.rootTree, [self.raw_dut_branch], [temp], self.settings.ana_events, self.settings.dutDetChs)
        self.rootFile.Close()
        print 'Done in', time.time() - time0, 'seconds'
        return mp.Array(np.ctypeslib.as_ctypes(temp)._type_, np.ctypeslib.as_ctypes(temp))

    def CalculateDevicesPedestals2(self):
        jobs_chunks = [self.devices[i:i+self.num_parallel] for i in xrange(0, len(self.devices), self.num_parallel)]
        for row in jobs_chunks:
            self.device_ped_process = []
            for it, device in enumerate(row):
                if it == len(row) - 1:
                    print 'Calculating pedestals for', device, '. The progress is shown below:'
                    if device.startswith('tel'):
                        temp = PedestalDeviceCalculations2('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, True, self.tel_ADC_all_mp, self.tel_ADC_mean_all_mp, self.tel_ADC_sigma_all_mp, self.tel_ADC_is_hit_all_mp, self.tel_ADC_is_seed_all_mp, self.dic_device_to_pos[device])
                    else:
                        temp = PedestalDeviceCalculations2('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, True, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, -1)

                else:
                    print 'Calculating pedestals for', device
                    if device.startswith('tel'):
                        temp = PedestalDeviceCalculations2('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, False, self.tel_ADC_all_mp, self.tel_ADC_mean_all_mp, self.tel_ADC_sigma_all_mp, self.tel_ADC_is_hit_all_mp, self.tel_ADC_is_seed_all_mp, self.dic_device_to_pos[device])
                    else:
                        temp = PedestalDeviceCalculations2('{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), device, False, self.dut_ADC_all_mp, self.dut_ADC_mean_all_mp, self.dut_ADC_sigma_all_mp, self.dut_ADC_is_hit_all_mp, self.dut_ADC_is_seed_all_mp, -1)
                temp.start()
                self.device_ped_process.append(temp)
            for j in self.device_ped_process:
                j.join()
        ipdb.set_trace()




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



def main():
    z = PedestalCalculations()


if __name__ == '__main__':
    main()
