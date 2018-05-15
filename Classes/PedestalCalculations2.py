# from ROOT import TFile, gSystem, TStopwatch, TDatime, TTree, TH1F
import ROOT as ro
from optparse import OptionParser
import os, logging, sys, shutil
sys.path.append('/home/sandiego/Small3dAnalysis/Classes')  # TODO: HARDCODED!!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
from Settings import Settings
from Utils import *
# from PedestalDeviceCalculations2 import PedestalDeviceCalculations2
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

class PedestalCalculations2:
    def __init__(self, settings=Settings()):
        print 'Creating PedestalCalculations instance'
        self.settings = settings
        self.out_dir = self.settings.output_dir
        self.sub_dir = self.settings.sub_dir
        self.file_name = self.settings.file_name
        self.run = self.settings.run
        self.tree_name = self.settings.tree_name
        self.num_parallel = self.settings.num_parallel - 2
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

        self.dic_device_to_pos = {'telx0': 0, 'tely0': 1, 'telx1': 2, 'tely1': 3, 'telx2': 4, 'tely2': 5, 'telx3': 6, 'tely3': 7}

        self.tel_ADC_mean_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='float32')
        # self.tel_ADC_sigma_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='float32')
        # self.tel_ADC_is_hit_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='?')
        # self.tel_ADC_is_seed_all = np.zeros((self.settings.telDetectors, self.settings.telDetChs, self.settings.ana_events), dtype='?')
        #
        self.dut_ADC_mean_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='float32')
        # self.dut_ADC_sigma_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='float32')
        # self.dut_ADC_is_hit_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='?')
        # self.dut_ADC_is_seed_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='?')
        # self.dut_ADC_mean_all = np.zeros((self.settings.dutDetChs, self.settings.ana_events), dtype='float32')

        # self.tel_ADC_mean_all_c = np.ctypeslib.as_ctypes(self.tel_ADC_mean_all)

        # self.tel_ADC_mean_all_mp = mpsc.Array(self.tel_ADC_mean_all_c._type_, self.tel_ADC_mean_all_c, lock=False)
        self.tel_ADC_mean_all_mp = mp.Array(np.ctypeslib.as_ctypes(self.tel_ADC_mean_all)._type_, np.ctypeslib.as_ctypes(self.tel_ADC_mean_all), lock=False)

        # self.dut_ADC_mean_all_c = np.ctypeslib.as_ctypes(self.dut_ADC_mean_all)
        self.dut_ADC_mean_all_mp = mp.Array(np.ctypeslib.as_ctypes(self.dut_ADC_mean_all)._type_, np.ctypeslib.as_ctypes(self.dut_ADC_mean_all), lock=False)

    def CalculateDevicesPedestals2(self):
        jobs_chunks = [self.devices[i:i+self.num_parallel] for i in xrange(0, len(self.devices), self.num_parallel)]
        for row in jobs_chunks:
            self.device_ped_process = []
            # manager = mp.Manager()
            # out_dic = manager.dict()
            for it, device in enumerate(row):
                if it == len(row) - 1:
                    print 'Calculating pedestals for', device, '. The progress is shown below:'
                    args0 = (device, True, self.dic_device_to_pos[device]) if device.startswith('tel') else (device, True, 0)
                    temp = mp.Process(target=self.PedestalDeviceCalculations, args=args0)
                else:
                    print 'Calculating pedestals for', device
                    args0 = (device, False, self.dic_device_to_pos[device]) if device.startswith('tel') else (device, False, 0)
                    temp = mp.Process(target=self.PedestalDeviceCalculations, args=args0)
                temp.start()
                self.device_ped_process.append(temp)
            for j in self.device_ped_process:
                j.join()
            ipdb.set_trace()
            
    def PedestalDeviceCalculations(self, device='dut', show_pb=False, det_index=0):
        do_cmc = True if device == 'dut' else False
        out_dir = self.settings.output_dir
        sub_dir = self.settings.sub_dir
        file_name = self.settings.file_name
        tree_name = self.settings.tree_name
        slide_leng = self.settings.sliding_length
        run_no = self.settings.run
        ana_events = self.settings.ana_events
        ped_branches = ['diaPed', 'diaPedSigma', 'cm', 'diaPedCmc', 'diaPedSigmaCmc', 'diaSignal', 'diaSignalCmc']
        raw_tel_branches_dic = {0: 'D0X_ADC', 1: 'D0Y_ADC', 2: 'D1X_ADC', 3: 'D1Y_ADC', 4: 'D2X_ADC', 5: 'D2Y_ADC', 6: 'D3X_ADC', 7: 'D3Y_ADC'}
        device_to_position = {'telx0': 0, 'tely0': 1, 'telx1': 2, 'tely1': 3, 'telx2': 4, 'tely2': 5, 'telx3': 6, 'tely3': 7}
        raw_dut_branch = 'DiaADC'
        rootFile = ro.TFile()
        rootTree = ro.TTree()
        utils = Utils()
        dest_path_stem = out_dir + '/' + sub_dir + '/' + str(run_no) + '/' + device

        read_branch = raw_tel_branches_dic[0] if device == 'telx0' else raw_tel_branches_dic[1] if device == 'tely0' else raw_tel_branches_dic[2] if device == 'telx1' else raw_tel_branches_dic[3] if device == 'tely1' else raw_tel_branches_dic[
            4] if device == 'telx2' else raw_tel_branches_dic[5] if device == 'tely2' else raw_tel_branches_dic[6] if device == 'telx3' else raw_tel_branches_dic[7] if device == 'tely3' else raw_dut_branch
        hit_factor = self.settings.clust_hit[0] if device == 'telx0' else self.settings.clust_hit[1] if device == 'tely0' else self.settings.clust_hit[2] if device == 'telx1' else self.settings.clust_hit[3] if device == 'tely1' else self.settings.clust_hit[
            4] if device == 'telx2' else self.settings.clust_hit[5] if device == 'tely2' else self.settings.clust_hit[6] if device == 'telx3' else self.settings.clust_hit[7] if device == 'tely3' else self.settings.clust_hit[8]
        seed_factor = self.settings.clust_seed[0] if device == 'telx0' else self.settings.clust_seed[1] if device == 'tely0' else self.settings.clust_seed[2] if device == 'telx1' else self.settings.clust_seed[3] if device == 'tely1' else self.settings.clust_seed[
            4] if device == 'telx2' else self.settings.clust_seed[5] if device == 'tely2' else self.settings.clust_seed[6] if device == 'telx3' else self.settings.clust_seed[7] if device == 'tely3' else self.settings.clust_seed[8]
        np_type = self.settings.dut_np_data_type if device == 'dut' else self.settings.tel_np_data_type
        chs = self.settings.dutDetChs if device == 'dut' else self.settings.telDetChs
        device_ADC = np.zeros((chs, slide_leng), dtype=np_type)
        device_ped = np.zeros((chs, slide_leng), dtype='float64')
        device_sigma = np.zeros((chs, slide_leng), dtype='float64')
        device_signal = np.zeros((chs, slide_leng), dtype='float64')
        device_ADC_mean = np.zeros((chs, ana_events), dtype='float64')
        device_ADC_sigma = np.zeros((chs, ana_events), dtype='float64')
        device_ADC_is_ped = np.zeros((chs, ana_events), dtype='?')
        device_ADC_is_hit = np.zeros((chs, ana_events), dtype='?')
        device_ADC_is_seed = np.zeros((chs, ana_events), dtype='?')

        adc = np.zeros(chs, dtype=np_type)
        mean = np.zeros(chs, dtype='float64')
        sigma = np.zeros(chs, dtype='float64')
        mean_sq = np.zeros(chs, dtype='float64')
        elem = np.zeros(chs, dtype='uint16')

        device_ADC_all = np.zeros((chs, ana_events), dtype=np_type)
        det_index = det_index
        
        rootFile, rootTree = Open_RootFile_Load_Tree('{d}/{s}/{r}/{f}.root'.format(d=out_dir, s=sub_dir, r=run_no, f=file_name), treename=tree_name, mode='READ')
        if show_pb:
            time0 = time.time()
            print 'Getting all events...', ; sys.stdout.flush()
        Draw_Branches_For_Get_Val(rootTree, [read_branch], start_ev=0, n_entries=ana_events, option='goff')
        Get_Branches_Value_To_Numpy(rootTree, [read_branch], [device_ADC_all], ana_events, chs)
        rootFile.Close()
        if show_pb:
            print 'Done in', time.time() - time0, 'seconds'
        device_ADC = device_ADC_all[:, :slide_leng]

        for ch, value in enumerate(device_ADC.mean(1, dtype='float64')):
            device_ADC_mean[ch,:slide_leng] = value
        for ch, value in enumerate(device_ADC.std(1, dtype='float64')):
            device_ADC_sigma[ch,:slide_leng] = value
        for ch in xrange(chs):
            device_ped[ch].fill(device_ADC_mean[ch, 0])
            device_sigma[ch].fill(device_ADC_sigma[ch, 0])
        device_signal = device_ADC - device_ped
        for it in xrange(7):
            condition_p = device_signal < hit_factor * device_sigma
            condition_h = np.bitwise_and(hit_factor * device_sigma <= device_signal, device_signal < seed_factor * device_sigma)
            condition_s = np.bitwise_and(np.bitwise_not(condition_p), np.bitwise_not(condition_h))
            adc_cond = [np.extract(condition_p[ch], device_ADC[ch]) for ch in xrange(chs)]
            for ch, adc_cond_i in enumerate(adc_cond):
                device_ADC_mean[ch].fill(adc_cond_i.mean(dtype='float64'))
                device_ADC_sigma[ch].fill(adc_cond_i.std(dtype='float64'))
                device_ped[ch].fill(adc_cond_i.mean(dtype='float64'))
                device_sigma[ch].fill(adc_cond_i.std(dtype='float64'))
            device_signal = device_ADC - device_ped
            device_ADC_is_ped[:, :slide_leng] = condition_p
            device_ADC_is_hit[:, :slide_leng] = condition_h
            device_ADC_is_seed[:, :slide_leng] = condition_s

        mean = device_ADC_mean[:, 0]
        sigma = device_ADC_sigma[:, 0]
        mean_sq = np.add(np.power(sigma, 2), np.power(mean, 2), dtype='float64')
        elem = device_ADC_is_ped[:, :slide_leng].sum(axis=1)

        del device_ADC, device_sigma, device_ped, device_signal
        
        if show_pb:
            utils.CreateProgressBar(maxVal=ana_events)
            utils.bar.start()
            utils.bar.update(slide_leng)
        for ev in xrange(slide_leng, ana_events):
            adc_new = device_ADC_all[:, ev]
            signal_new = adc_new - mean
            adc_old = device_ADC_all[:, ev - slide_leng]
            condition_p = signal_new < hit_factor * sigma
            condition_h = np.bitwise_and(hit_factor * sigma <= signal_new, signal_new < seed_factor * sigma)
            condition_s = np.bitwise_and(np.bitwise_not(condition_p), np.bitwise_not(condition_h))

            mean1 = mean
            # mean2 = (mean * elem - adc_old) / (elem - 1.0)
            mean2 = np.divide(np.subtract(np.multiply(mean, elem, dtype='float64'), adc_old, dtype='float64'), np.subtract(elem, 1.0, dtype='float64'), dtype='float64')

            mean1_sq = mean_sq
            # mean2_sq = (mean_sq * elem - adc_old ** 2) / (elem - 1.0)
            mean2_sq = np.divide(np.subtract(np.multiply(mean_sq, elem, dtype='float64'), np.power(adc_old, 2.0, dtype='float64'), dtype='float64'), np.subtract(elem, 1.0, dtype='float64'), dtype='float64')

            elem1 = elem
            elem2 = elem - 1

            condition_old = device_ADC_is_ped[:, ev - slide_leng]
            mean = np.where(condition_old, mean2, mean1)
            mean_sq = np.where(condition_old, mean2_sq, mean1_sq)
            elem = np.where(condition_old, elem2, elem1)

            mean1 = mean
            # mean2 = (mean * elem + adc_new) / (elem + 1.0)
            mean2 = np.divide(np.add(np.multiply(mean, elem, dtype='float64'), adc_new, dtype='float64'), np.add(elem, 1.0, dtype='float64'), dtype='float64')

            mean1_sq = mean_sq
            # mean2_sq = (mean_sq * elem + adc_new ** 2) / (elem + 1.0)
            mean2_sq = np.divide(np.add(np.multiply(mean_sq, elem, dtype='float64'), np.power(adc_new, 2.0, dtype='float64'), dtype='float64'), np.add(elem, 1.0, dtype='float64'), dtype='float64')

            elem1 = elem
            elem2 = elem + 1

            mean = np.where(condition_p, mean2, mean1)
            mean_sq = np.where(condition_p, mean2_sq, mean1_sq)
            elem = np.where(condition_p, elem2, elem1)

            # sigma = (mean_sq - mean ** 2) ** 0.5
            sigma = np.sqrt(np.subtract(mean_sq, np.power(mean, 2.0, dtype='float64'), dtype='float64'), dtype='float64')

            device_ADC_mean[:, ev] = mean
            device_ADC_sigma[:, ev] = sigma
            device_ADC_is_ped[:, ev] = condition_p
            device_ADC_is_hit[:, ev] = condition_h
            device_ADC_is_seed[:, ev] = condition_s
            if show_pb:
                utils.bar.update(ev + 1)
        if show_pb:
            utils.bar.finish()
            t0 = time.time()
            print 'Appending arrays to output...', ; sys.stdout.flush()
        # out_dic[device] = {'mean': device_ADC_mean.astype('float32'), 'sigma': device_ADC_sigma.astype('float32'), 'is_hit': device_ADC_is_hit, 'is_seed': device_ADC_is_seed}
        if device.startswith('tel'):
            self.tel_ADC_mean_all_mp[det_index] = np.ctypeslib.as_ctypes(device_ADC_mean.astype(dtype='float32'))
        else:
            self.dut_ADC_mean_all_mp = np.ctypeslib.as_ctypes(device_ADC_mean.astype(dtype='float32'))
        if show_pb:
            print 'Done in', time.time() - t0, 'seconds'



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
    z = PedestalCalculations2()


if __name__ == '__main__':
    main()
