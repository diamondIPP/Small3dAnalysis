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
# from threading import Thread
# from multiprocessing import Process, Lock
import multiprocessing as mp
# import multiprocessing.sharedctypes as mpsc

__author__ = 'DA'



class PedestalDeviceCalculations2(mp.Process):
    def __init__(self, settings_bin_path, device, show_progressbar, input_adc_array, out_array_mean, out_array_sigma, out_array_is_ped, out_array_is_hit, out_array_is_seed, out_array_chs_cm, out_array_cm, out_array_mean_cmc, out_array_sigma_cmc, out_array_is_ped_cmc, out_array_is_hit_cmc, out_array_is_seed_cmc, det_index):
        mp.Process.__init__(self)
        print 'Creating PedestalCalculations instance'
        self.show_pb = show_progressbar
        self.settings = Settings()
        self.device = device
        self.do_cmc = self.settings.do_cmc and (self.device == 'dut')
        self.LoadSettingsBinary(settings_bin_path)
        self.out_dir = self.settings.output_dir
        self.sub_dir = self.settings.sub_dir
        self.file_name = self.settings.file_name
        self.tree_name = self.settings.tree_name
        self.slide_leng = self.settings.sliding_length
        self.run_no = self.settings.run
        self.ana_events = self.settings.ana_events
        self.ped_branches = ['diaPed', 'diaPedSigma', 'cm', 'diaPedCmc', 'diaPedSigmaCmc', 'diaSignal', 'diaSignalCmc']
        self.raw_tel_branches_dic = {0: 'D0X_ADC', 1: 'D0Y_ADC', 2: 'D1X_ADC', 3: 'D1Y_ADC', 4: 'D2X_ADC', 5: 'D2Y_ADC', 6: 'D3X_ADC', 7: 'D3Y_ADC'}
        self.device_to_position = {'telx0': 0, 'tely0': 1, 'telx1': 2, 'tely1': 3, 'telx2': 4, 'tely2': 5, 'telx3': 6, 'tely3': 7}
        self.raw_dut_branch = 'DiaADC'
        self.rootFile = ro.TFile()
        self.rootTree = ro.TTree()
        self.utils = Utils()
        self.dest_path_stem = self.out_dir + '/' + self.sub_dir + '/' + str(self.run_no) + '/' + self.device

        self.read_branch = self.raw_tel_branches_dic[0] if self.device == 'telx0' else self.raw_tel_branches_dic[1] if self.device == 'tely0' else self.raw_tel_branches_dic[2] if self.device == 'telx1' else self.raw_tel_branches_dic[3] if self.device == 'tely1' else self.raw_tel_branches_dic[4] if self.device == 'telx2' else self.raw_tel_branches_dic[5] if self.device == 'tely2' else self.raw_tel_branches_dic[6] if self.device == 'telx3' else self.raw_tel_branches_dic[7] if self.device == 'tely3' else self.raw_dut_branch
        self.hit_factor = self.settings.clust_hit[0] if self.device == 'telx0' else self.settings.clust_hit[1] if self.device == 'tely0' else self.settings.clust_hit[2] if self.device == 'telx1' else self.settings.clust_hit[3] if self.device == 'tely1' else self.settings.clust_hit[4] if self.device == 'telx2' else self.settings.clust_hit[5] if self.device == 'tely2' else self.settings.clust_hit[6] if self.device == 'telx3' else self.settings.clust_hit[7] if self.device == 'tely3' else self.settings.clust_hit[8]
        self.seed_factor = self.settings.clust_seed[0] if self.device == 'telx0' else self.settings.clust_seed[1] if self.device == 'tely0' else self.settings.clust_seed[2] if self.device == 'telx1' else self.settings.clust_seed[3] if self.device == 'tely1' else self.settings.clust_seed[4] if self.device == 'telx2' else self.settings.clust_seed[5] if self.device == 'tely2' else self.settings.clust_seed[6] if self.device == 'telx3' else self.settings.clust_seed[7] if self.device == 'tely3' else self.settings.clust_seed[8]
        self.np_type = self.settings.dut_np_data_type if self.device == 'dut' else self.settings.tel_np_data_type
        self.chs = self.settings.dutDetChs if self.device == 'dut' else self.settings.telDetChs

        # Numpy arrays for the calculations
        self.device_ADC_mean = np.zeros((self.chs, self.ana_events), dtype='float32')
        self.device_ADC_sigma = np.zeros((self.chs, self.ana_events), dtype='float32')
        self.device_ADC_is_ped = np.zeros((self.chs, self.ana_events), dtype='?')
        self.device_ADC_is_hit = np.zeros((self.chs, self.ana_events), dtype='?')
        self.device_ADC_is_seed = np.zeros((self.chs, self.ana_events), dtype='?')

        if self.do_cmc:
            self.device_channels_cm = np.zeros((self.chs, self.ana_events), dtype='?')
            self.device_cm = np.zeros(self.ana_events, dtype='float32')
            self.device_ADC_mean_cmc = np.zeros((self.chs, self.ana_events), dtype='float32')
            self.device_ADC_sigma_cmc = np.zeros((self.chs, self.ana_events), dtype='float32')
            self.device_is_cm = np.zeros((self.chs, self.ana_events), dtype='?')
            self.device_is_ped_cmc = np.zeros((self.chs, self.ana_events), dtype='?')
            self.device_is_hit_cmc = np.zeros((self.chs, self.ana_events), dtype='?')
            self.device_is_seed_cmc = np.zeros((self.chs, self.ana_events), dtype='?')

        # Numpy arrays for event based calculations
        self.adc = np.zeros(self.chs, dtype=self.np_type)
        self.mean = np.zeros(self.chs, dtype='float32')
        self.sigma = np.zeros(self.chs, dtype='float32')
        self.mean_sq = np.zeros(self.chs, dtype='float32')
        self.elem = np.zeros(self.chs, dtype='uint16')
        if self.do_cmc:
            self.cm = np.zeros(1, dtype='float32')
            self.adc_cmc = np.zeros(self.chs, dtype='float32')
            self.mean_cmc = np.zeros(self.chs, dtype='float32')
            self.sigma_cmc = np.zeros(self.chs, dtype='float32')
            self.mean_sq_cmc = np.zeros(self.chs, dtype='float32')
            self.elem_cmc = np.zeros(self.chs, dtype='uint16')

        # channels that are masked for common mode calculation because they are screened or noisy

        self.is_not_masked = np.ones(self.chs, '?') if self.device.startswith('tel') else np.array([ch not in set(self.settings.noisy) | set(self.settings.screened) for ch in xrange(self.chs)], '?')

        # global variables passed as reference for shared memory ctype vectors
        global in_adc_array, out_array_m, out_array_s, out_array_is_p, out_array_is_h, out_array_is_s
        in_adc_array = input_adc_array
        out_array_m = out_array_mean
        out_array_s = out_array_sigma
        out_array_is_h = out_array_is_hit
        out_array_is_s = out_array_is_seed
        out_array_is_p = out_array_is_ped
        self.det_index = det_index

        if self.do_cmc:
            global out_array_chs_c, out_array_c, out_array_mean_cm, out_array_sigma_cm, out_array_is_ped_cm, out_array_is_hit_cm, out_array_is_seed_cm
            out_array_chs_c = out_array_chs_cm
            out_array_c = out_array_cm
            out_array_mean_cm = out_array_mean_cmc
            out_array_sigma_cm = out_array_sigma_cmc
            out_array_is_ped_cm = out_array_is_ped_cmc
            out_array_is_hit_cm = out_array_is_hit_cmc
            out_array_is_seed_cm = out_array_is_seed_cmc

        if self.device.startswith('tel'):
            self.device_ADC_all = np.ctypeslib.as_array(in_adc_array.get_obj())[self.det_index]
        else:
            self.device_ADC_all = np.ctypeslib.as_array(in_adc_array.get_obj())

    def LoadSettingsBinary(self, settings_path):
        if os.path.isfile(settings_path):
            with open(settings_path, 'rb') as fs:
                self.settings = pickle.load(fs)
        else:
            ExitMessage('Settings file does not exist!!!!', value=os.EX_OSFILE)

    def run(self):
        self.CalculatePedestals()

    def CalculatePedestals(self):
        self.CalculateStartingPedestals()
        if self.show_pb:
            self.utils.CreateProgressBar(maxVal=self.ana_events)
            self.utils.bar.start()
            self.utils.bar.update(self.slide_leng)

        for ev in xrange(self.slide_leng, self.ana_events):
            adc_new = self.device_ADC_all[:, ev]
            signal_new = np.subtract(adc_new, self.mean, dtype='float64')
            adc_old = self.device_ADC_all[:, ev - self.slide_leng]
            condition_p = np.abs(signal_new, dtype='float64') < np.multiply(self.hit_factor, self.sigma, dtype='float64')
            condition_h = np.bitwise_and(np.multiply(self.hit_factor, self.sigma, dtype='float64') <= signal_new, signal_new < np.multiply(self.seed_factor, self.sigma, dtype='float64'))
            condition_s = signal_new >= np.multiply(self.seed_factor, self.sigma, dtype='float64')
            # condition_s = np.bitwise_and(np.bitwise_not(condition_p), np.bitwise_not(condition_h))

            mean1 = self.mean
            # mean2 = (self.mean * self.elem - adc_old) / (self.elem - 1.0)
            mean2 = np.divide(np.subtract(np.multiply(self.mean, self.elem, dtype='float64'), adc_old, dtype='float64'), np.subtract(self.elem, 1.0, dtype='float64'), dtype='float64')

            mean1_sq = self.mean_sq
            # mean2_sq = (self.mean_sq * self.elem - adc_old ** 2) / (self.elem - 1.0)
            mean2_sq = np.divide(np.subtract(np.multiply(self.mean_sq, self.elem, dtype='float64'), np.power(adc_old, 2.0, dtype='float64'), dtype='float64'), np.subtract(self.elem, 1.0, dtype='float64'), dtype='float64')

            elem1 = self.elem
            elem2 = self.elem - 1

            condition_old = self.device_ADC_is_ped[:, ev - self.slide_leng]
            self.mean = np.where(condition_old, mean2, mean1)
            self.mean_sq = np.where(condition_old, mean2_sq, mean1_sq)
            self.elem = np.where(condition_old, elem2, elem1)

            mean1 = self.mean
            # mean2 = (self.mean * self.elem + adc_new) / (self.elem + 1.0)
            mean2 = np.divide(np.add(np.multiply(self.mean, self.elem, dtype='float64'), adc_new, dtype='float64'), np.add(self.elem, 1.0, dtype='float64'), dtype='float64')

            mean1_sq = self.mean_sq
            # mean2_sq = (self.mean_sq * self.elem + adc_new ** 2) / (self.elem + 1.0)
            mean2_sq = np.divide(np.add(np.multiply(self.mean_sq, self.elem, dtype='float64'), np.power(adc_new, 2.0, dtype='float64'), dtype='float64'), np.add(self.elem, 1.0, dtype='float64'), dtype='float64')

            elem1 = self.elem
            elem2 = self.elem + 1

            self.mean = np.where(condition_p, mean2, mean1)
            self.mean_sq = np.where(condition_p, mean2_sq, mean1_sq)
            self.elem = np.where(condition_p, elem2, elem1)

            self.sigma = np.sqrt(np.subtract(self.mean_sq, np.power(self.mean, 2.0, dtype='float64'), dtype='float64'), dtype='float64')

            self.device_ADC_mean[:, ev] = self.mean
            self.device_ADC_sigma[:, ev] = self.sigma
            self.device_ADC_is_ped[:, ev] = condition_p
            self.device_ADC_is_hit[:, ev] = condition_h
            self.device_ADC_is_seed[:, ev] = condition_s

            if self.do_cmc:
                adc_old_cmc = np.subtract(self.device_ADC_all[:, ev - self.slide_leng], self.device_cm[ev - self.slide_leng])
                signal_new_cmc = np.subtract(adc_new, self.mean_cmc, dtype='float64')
                condition_cm = np.bitwise_and((abs(signal_new_cmc) < np.multiply(self.settings.cm_cut, self.sigma_cmc, dtype='float64')), self.is_not_masked)
                cm = np.extract(condition_cm, signal_new_cmc).mean(dtype='float64') if condition_cm.any() else 0
                self.cm.fill(cm)
                adc_new_cmc = np.subtract(adc_new, self.cm, dtype='float64')
                signal_new_cmc = np.subtract(adc_new_cmc, self.mean_cmc, dtype='float64')
                condition_p_cmc = np.abs(signal_new_cmc, dtype='float64') < np.multiply(self.hit_factor, self.sigma_cmc, dtype='float64')
                condition_h_cmc = np.bitwise_and(np.multiply(self.hit_factor, self.sigma_cmc, dtype='float64') <= signal_new_cmc, signal_new_cmc < np.multiply(self.seed_factor, self.sigma_cmc, dtype='float64'))
                condition_s_cmc = signal_new_cmc >= np.multiply(self.seed_factor, self.sigma_cmc, dtype='float64')
                # condition_s_cmc = np.bitwise_and(np.bitwise_not(condition_p_cmc), np.bitwise_not(condition_h_cmc))

                mean1_cmc = self.mean_cmc
                mean2_cmc = np.divide(np.subtract(np.multiply(self.mean_cmc, self.elem_cmc, dtype='float64'), adc_old_cmc, dtype='float64'), np.subtract(self.elem_cmc, 1.0, dtype='float64'), dtype='float64')

                mean1_sq_cmc = self.mean_sq_cmc
                mean2_sq_cmc = np.divide(np.subtract(np.multiply(self.mean_sq_cmc, self.elem_cmc, dtype='float64'), np.power(adc_old_cmc, 2.0, dtype='float64'), dtype='float64'), np.subtract(self.elem_cmc, 1.0, dtype='float64'), dtype='float64')

                elem1_cmc = self.elem_cmc
                elem2_cmc = self.elem_cmc - 1

                condition_old_cmc = self.device_is_ped_cmc[:, ev - self.slide_leng]
                self.mean_cmc = np.where(condition_old_cmc, mean2_cmc, mean1_cmc)
                self.mean_sq_cmc = np.where(condition_old_cmc, mean2_sq_cmc, mean1_sq_cmc)
                self.elem_cmc = np.where(condition_old_cmc, elem2_cmc, elem1_cmc)

                mean1_cmc = self.mean_cmc
                mean2_cmc = np.divide(np.add(np.multiply(self.mean_cmc, self.elem_cmc, dtype='float64'), adc_new_cmc, dtype='float64'), np.add(self.elem_cmc, 1.0, dtype='float64'), dtype='float64')

                mean1_sq_cmc = self.mean_sq_cmc
                mean2_sq_cmc = np.divide(np.add(np.multiply(self.mean_sq_cmc, self.elem_cmc, dtype='float64'), np.power(adc_new_cmc, 2.0, dtype='float64'), dtype='float64'), np.add(self.elem_cmc, 1.0, dtype='float64'), dtype='float64')

                elem1_cmc = self.elem_cmc
                elem2_cmc = self.elem_cmc + 1

                self.mean_cmc = np.where(condition_p_cmc, mean2_cmc, mean1_cmc)
                self.mean_sq_cmc = np.where(condition_p_cmc, mean2_sq_cmc, mean1_sq_cmc)
                self.elem_cmc = np.where(condition_p_cmc, elem2_cmc, elem1_cmc)

                self.sigma_cmc = np.sqrt(np.subtract(self.mean_sq_cmc, np.power(self.mean_cmc, 2.0, dtype='float64'), dtype='float64'), dtype='float64')

                self.device_channels_cm[:, ev] = condition_cm
                self.device_cm[ev] = self.cm
                self.device_ADC_mean_cmc[:, ev] = self.mean_cmc
                self.device_ADC_sigma_cmc[:, ev] = self.sigma_cmc
                self.device_is_ped_cmc[:, ev] = condition_p_cmc
                self.device_is_hit_cmc[:, ev] = condition_h_cmc
                self.device_is_seed_cmc[:, ev] = condition_s_cmc

            if self.show_pb:
                self.utils.bar.update(ev + 1)
        if self.show_pb:
            self.utils.bar.finish()
            t0 = time.time()
            print 'Appending arrays to output...', ; sys.stdout.flush()

        if self.device.startswith('tel'):
            np.ctypeslib.as_array(out_array_m.get_obj())[self.det_index] = self.device_ADC_mean
            np.ctypeslib.as_array(out_array_s.get_obj())[self.det_index] = self.device_ADC_sigma
            np.ctypeslib.as_array(out_array_is_p.get_obj())[self.det_index] = self.device_ADC_is_ped.astype('uint8')
            np.ctypeslib.as_array(out_array_is_h.get_obj())[self.det_index] = self.device_ADC_is_hit.astype('uint8')
            np.ctypeslib.as_array(out_array_is_s.get_obj())[self.det_index] = self.device_ADC_is_seed.astype('uint8')
        else:
            np.ctypeslib.as_array(out_array_m.get_obj())[:] = self.device_ADC_mean
            np.ctypeslib.as_array(out_array_s.get_obj())[:] = self.device_ADC_sigma
            np.ctypeslib.as_array(out_array_is_p.get_obj())[:] = self.device_ADC_is_ped.astype('uint8')
            np.ctypeslib.as_array(out_array_is_h.get_obj())[:] = self.device_ADC_is_hit.astype('uint8')
            np.ctypeslib.as_array(out_array_is_s.get_obj())[:] = self.device_ADC_is_seed.astype('uint8')
            if self.do_cmc:
                np.ctypeslib.as_array(out_array_chs_c.get_obj())[:] = self.device_channels_cm
                np.ctypeslib.as_array(out_array_c.get_obj())[:] = self.device_cm
                np.ctypeslib.as_array(out_array_mean_cm.get_obj())[:] = self.device_ADC_mean_cmc
                np.ctypeslib.as_array(out_array_sigma_cm.get_obj())[:] = self.device_ADC_sigma_cmc
                np.ctypeslib.as_array(out_array_is_ped_cm.get_obj())[:] = self.device_is_ped_cmc
                np.ctypeslib.as_array(out_array_is_hit_cm.get_obj())[:] = self.device_is_hit_cmc
                np.ctypeslib.as_array(out_array_is_seed_cm.get_obj())[:] = self.device_is_seed_cmc
        if self.show_pb:
            print 'Done in', time.time() - t0, 'seconds'

    def CalculateStartingPedestals(self):
        device_ADC = self.device_ADC_all[:, :self.slide_leng]

        for ch, value in enumerate(device_ADC.mean(1, dtype='float64')):
            self.device_ADC_mean[ch, :self.slide_leng] = value
        for ch, value in enumerate(device_ADC.std(1, dtype='float64')):
            self.device_ADC_sigma[ch, :self.slide_leng] = value

        device_signal = device_ADC - self.device_ADC_mean[:, :self.slide_leng]

        if self.do_cmc:
            not_masked_array = np.zeros((self.chs, self.slide_leng), '?')
            for ch, val in enumerate(self.is_not_masked):
                not_masked_array[ch].fill(val)
            self.device_cm[:self.slide_leng] = (np.array([np.extract(not_masked_array[:, ev], device_signal[:, ev]) for ev in xrange(self.slide_leng)]).T).mean(axis=0, dtype='float64')
            cm_array = np.array([self.device_cm[:self.slide_leng] for ch in xrange(self.chs)], 'float32')
            device_ADC_cmc = np.subtract(device_ADC, cm_array, dtype='float64')
            for ch, value in enumerate(device_ADC_cmc.mean(axis=1, dtype='float64')):
                self.device_ADC_mean_cmc[ch, :self.slide_leng] = value
            for ch, value in enumerate(device_ADC_cmc.std(axis=1, dtype='float64')):
                self.device_ADC_sigma_cmc[ch, :self.slide_leng] = value
            device_signal_cmc = device_ADC_cmc - self.device_ADC_mean_cmc[:, :self.slide_leng]

        for it in xrange(7):
            condition_p = np.abs(device_signal, dtype='float64') < np.multiply(self.hit_factor,  self.device_ADC_sigma[:, :self.slide_leng], dtype='float64')
            condition_h = np.bitwise_and(np.multiply(self.hit_factor,  self.device_ADC_sigma[:, :self.slide_leng], dtype='float64') <= device_signal, device_signal < np.multiply(self.seed_factor,  self.device_ADC_sigma[:, :self.slide_leng], dtype='float64'))
            condition_s = device_signal >= np.multiply(self.seed_factor, self.device_ADC_sigma[:, :self.slide_leng], dtype='float64')
            # condition_s = np.bitwise_and(np.bitwise_not(condition_p), np.bitwise_not(condition_h))
            adc_cond = [np.extract(condition_p[ch], device_ADC[ch]) for ch in xrange(self.chs)]

            if self.do_cmc:
                condition_cm = np.bitwise_and((abs(device_signal_cmc) < np.multiply(self.settings.cm_cut, self.device_ADC_sigma_cmc[:, :self.slide_leng], dtype='float64')), not_masked_array)
                self.device_cm[:self.slide_leng] = np.array([np.extract(condition_cm[:, ev], device_signal_cmc[:, ev]).mean(dtype='float64') for ev in xrange(self.slide_leng)])
                cm_array = np.array([self.device_cm[:self.slide_leng] for ch in xrange(self.chs)], 'float32')
                device_ADC_cmc = np.subtract(device_ADC, cm_array, dtype='float64')
                device_signal_cmc = np.subtract(device_ADC_cmc, self.device_ADC_mean_cmc[:, :self.slide_leng], dtype='float64')
                condition_p_cmc = np.abs(device_signal_cmc, dtype='float64') < np.multiply(self.hit_factor, self.device_ADC_sigma_cmc[:, :self.slide_leng], dtype='float64')
                condition_h_cmc = np.bitwise_and(np.multiply(self.hit_factor, self.device_ADC_sigma_cmc[:, :self.slide_leng], dtype='float64') <= device_signal_cmc, device_signal_cmc < np.multiply(self.seed_factor, self.device_ADC_sigma_cmc[:, :self.slide_leng], dtype='float64'))
                condition_s_cmc = device_signal_cmc >= np.multiply(self.seed_factor, self.device_ADC_sigma_cmc[:, :self.slide_leng], dtype='float64')
                # condition_s_cmc = np.bitwise_and(np.bitwise_not(condition_p_cmc), np.bitwise_not(condition_h_cmc))
                adc_cond_cmc = [np.extract(condition_p_cmc[ch], device_ADC_cmc[ch]) for ch in xrange(self.chs)]

            for ch in xrange(self.chs):
                self.device_ADC_mean[ch, :self.slide_leng].fill(adc_cond[ch].mean(dtype='float64'))
                self.device_ADC_sigma[ch, :self.slide_leng].fill(adc_cond[ch].std(dtype='float64'))

                if self.do_cmc:
                    self.device_ADC_mean_cmc[ch, :self.slide_leng].fill(adc_cond_cmc[ch].mean(dtype='float64'))
                    self.device_ADC_sigma_cmc[ch, :self.slide_leng].fill(adc_cond_cmc[ch].std(dtype='float64'))

            device_signal = device_ADC - self.device_ADC_mean[:, :self.slide_leng]
            self.device_ADC_is_ped[:, :self.slide_leng] = condition_p
            self.device_ADC_is_hit[:, :self.slide_leng] = condition_h
            self.device_ADC_is_seed[:, :self.slide_leng] = condition_s

            if self.do_cmc:
                self.device_channels_cm[:, :self.slide_leng] = condition_cm
                device_signal_cmc = np.subtract(device_ADC, self.device_ADC_mean_cmc[:, :self.slide_leng], dtype='float64')
                self.device_is_ped_cmc[:, :self.slide_leng] = condition_p_cmc
                self.device_is_hit_cmc[:, :self.slide_leng] = condition_h_cmc
                self.device_is_seed_cmc[:, :self.slide_leng] = condition_s_cmc

        self.mean = self.device_ADC_mean[:, 0]
        self.sigma = self.device_ADC_sigma[:, 0]
        self.mean_sq = np.add(np.power(self.sigma, 2, dtype='float64'), np.power(self.mean, 2, dtype='float64'), dtype='float64')
        self.elem = self.device_ADC_is_ped[:, :self.slide_leng].sum(axis=1)

        if self.do_cmc:
            self.mean_cmc = self.device_ADC_mean_cmc[:, self.slide_leng - 1]
            self.sigma_cmc = self.device_ADC_sigma_cmc[:, self.slide_leng - 1]
            self.mean_sq_cmc = np.add(np.power(self.sigma_cmc, 2, dtype='float64'), np.power(self.mean_cmc, 2, dtype='float64'), dtype='float64')
            self.elem_cmc = self.device_is_ped_cmc[:, :self.slide_leng].sum(axis=1)
            self.cm.fill(self.device_cm[self.slide_leng - 1])


def main():
    parser = OptionParser()
    parser.add_option('-s', '--settings', dest='settings', type='string', help='Settings binary file e.g. run_23002_full.settings')
    parser.add_option('-d', '--device', dest='device', type='string', default='dut', help='Device to analyse (dut, telx0, telx1, telx2, telx3, tely0, tely1, tely2, tely3)')
    parser.add_option('-p', '--progress', dest='progressbar', default=False, action='store_true', help='Shows progress bar of the process')

    (options, args) = parser.parse_args()
    settings_bin_path = str(options.settings)
    device = str(options.device)
    if device not in ['dut', 'telx0', 'telx1', 'telx2', 'telx3', 'tely0', 'tely1', 'tely2', 'tely3']:
        ExitMessage('device must be "dut" or "telx#" or "tely#". Exiting...', os.EX_IOERR)
    do_pb = bool(options.progressbar)
    z = PedestalDeviceCalculations2(settings_bin_path=settings_bin_path, device=device, show_progressbar=do_pb)


if __name__ == '__main__':
    main()
