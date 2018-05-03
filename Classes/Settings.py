# from ROOT import gSystem, TFile, TTree
import ROOT as ro
from copy import deepcopy
import os, sys
import progressbar
from ConfigParser import ConfigParser
from collections import namedtuple
# import ipdb
from Utils import *
import multiprocessing as mp

class Settings:
    def __init__(self, do_parallel=False):
        self.do_parallel = do_parallel
        if not self.do_parallel:
            self.num_parallel = 1
        else:
            num_cores = mp.cpu_count()
            self.num_parallel = 1 if num_cores < 4 else 2 if num_cores < 8 else 4 if num_cores < 16 else 8 if num_cores < 32 else 16

        self.run = 0
        self.total_events = 10000
        self.ana_events = 10000
        self.first_event = 0
        self.repeater_card = 2
        self.voltage = 100
        self.current_begin = 0
        self.current_end = 0
        self.dut_input = 0
        self.dut_saturation = 3367
        self.tel_saturation = 255
        self.data_dir = '/data/diamondSetup/diamondTestbeam/RD42-TestBeams/2017/cern_RD42_06_2017'
        self.output_dir = '/data/diamondSetup/diamondTestbeam/RD42-TestBeams/2017/cern_RD42_06_2017/output'
        self.sub_dir = 'default'
        self.analysis_path = '/home/sandiego/Small3dAnalysis'

        self.dut_num = 1
        self.not_connected = [0, 127]
        self.screened = [0, 127]
        self.noisy = [0, 127]

        self.dut_name = {i: 'diamond' for i in xrange(self.dut_num)}
        self.dut_pos = {i: 100 for i in xrange(self.dut_num)}
        self.dut_pitch = {i: 50 for i in xrange(self.dut_num)}
        self.dut_first_ch = {i: 1 for i in xrange(self.dut_num)}
        self.dut_last_ch = {i: 126 for i in xrange(self.dut_num)}
        self.dut_skip_ch = {i: [] for i in xrange(self.dut_num)}

        self.do_pedestal = False
        self.do_cluster = False
        self.do_cluster_ana = False
        self.do_alignment = False
        self.do_alignment_ana = False
        self.do_transparent = False
        self.do_3d = False

        self.tel_hit_factor = 5
        self.dut_hit_factor = 3
        self.do_cmc = True
        self.cm_cut = 4

        self.clust_seed = {0: 14, 1: 19, 2: 23, 3: 23, 4: 14, 5: 13, 6: 11, 7: 10, 8: 5}
        self.clust_hit = {0: 12, 1: 14, 2: 17, 3: 17, 4: 8, 5: 7, 6: 7, 7: 6, 8: 3}

        self.scint_fid_cut = {'xlow': 0, 'xhigh': 255, 'ylow': 0, 'yhigh': 255}

        self.fid_region = {i: {'xlow': 0, 'xhigh': 255, 'ylow': 0, 'yhigh': 255} for i in xrange(self.dut_num)}

        self.do_align_dut = True
        self.align_dut = 0
        self.z_coordinates = {0: 0.725, 1: 0.725, 2: 1.625, 3: 1.625, 4: 18.725, 5: 18.725, 6: 19.625, 7: 19.625, 8: 10.2}
        self.align_method = 'events'
        self.align_factor = 1000
        self.no_align_dut_chs = {0, 127}
        self.align_chi2_cut = 4

        self.max_transp_clust_size = 10
        self.save_transp_clust_size = 10
        self.do_analyse_alignment = False

        self.telDetChs = 256
        self.dutDetChs = 128
        self.telPlanes = 4
        self.telDetectors = 8
        self.event_per_file = 10000
        self.tel_mem = 2048
        self.dut_mem = 256
        self.file_name = ''
        self.tree_name = ''

    def ReadSettingsFile(self, in_file=''):
        if in_file != '':
            if os.path.isfile(in_file):
                pars = ConfigParser()
                pars.read(in_file)
                print 'Loading settings from file:', in_file

                if pars.has_section('BASIC'):
                    if pars.has_option('BASIC', 'run'):
                        self.run = pars.getint('BASIC', 'run')
                    else:
                        ExitMessage('Must specify run under [BASIC]. Exiting...')
                    if pars.has_option('BASIC', 'events'):
                        self.total_events = pars.getint('BASIC', 'events')
                    if pars.has_option('BASIC', 'ana_events'):
                        self.ana_events = pars.getint('BASIC', 'ana_events')
                    if pars.has_option('BASIC', 'first_event'):
                        self.first_event = pars.getint('BASIC', 'first_event')
                    if pars.has_option('BASIC', 'repeater_card'):
                        self.repeater_card = pars.getint('BASIC', 'repeater_card')
                    if pars.has_option('BASIC', 'voltage'):
                        self.voltage = pars.getfloat('BASIC', 'voltage')
                    if pars.has_option('BASIC', 'current_begin'):
                        self.current_begin = pars.getfloat('BASIC', 'current_begin')
                    if pars.has_option('BASIC', 'current_end'):
                        self.current_end = pars.getfloat('BASIC', 'current_end')
                    if pars.has_option('BASIC', 'dut_input'):
                        self.dut_input = pars.getint('BASIC', 'dut_input')
                    if pars.has_option('BASIC', 'dut_saturation'):
                        self.dut_saturation = pars.getint('BASIC', 'dut_saturation')
                    if pars.has_option('BASIC', 'data_dir'):
                        self.data_dir = pars.get('BASIC', 'data_dir')
                        if not os.path.isdir(self.data_dir):
                            ExitMessage('The specified data directory does not exist. Exiting...')
                    if pars.has_option('BASIC', 'output_dir'):
                        self.output_dir = pars.get('BASIC', 'output_dir')
                    if pars.has_option('BASIC', 'sub_dir'):
                        self.sub_dir = pars.get('BASIC', 'sub_dir')
                    if pars.has_option('BASIC', 'analysis_path'):
                        self.analysis_path = pars.get('BASIC', 'analysis_path')
                        if not os.path.isdir(self.analysis_path):
                            ExitMessage('Must give a valid analysis_path under [BASIC]. Exiting...')
                    if pars.has_option('BASIC', 'num_parallel') and self.do_parallel:
                        self.num_parallel = pars.getint('BASIC', 'num_parallel') if pars.getint('BASIC', 'num_parallel') <= mp.cpu_count() else self.num_parallel

                if pars.has_section('DUTS'):
                    if pars.has_option('DUTS', 'num'):
                        self.dut_num = pars.getint('DUTS', 'num')
                    if pars.has_option('DUTS', 'not_connected'):
                        self.not_connected = ChannelStringToArray(pars.get('DUTS', 'not_connected'))
                    if pars.has_option('DUTS', 'screened'):
                        self.screened = ChannelStringToArray(pars.get('DUTS', 'screened'))
                    if pars.has_option('DUTS', 'noisy'):
                        self.noisy = ChannelStringToArray(pars.get('DUTS', 'noisy'))

                for dut in xrange(self.dut_num):
                    if pars.has_section('DUT{d}'.format(d=dut)):
                        if pars.has_option('DUT{d}'.format(d=dut), 'name'):
                            self.dut_name[dut] = pars.get('DUT{d}'.format(d=dut), 'name')
                        if pars.has_option('DUT{d}'.format(d=dut), 'x0'):
                            self.dut_pos[dut] = pars.getfloat('DUT{d}'.format(d=dut), 'x0')
                        if pars.has_option('DUT{d}'.format(d=dut), 'pitch'):
                            self.dut_pitch[dut] = pars.getfloat('DUT{d}'.format(d=dut), 'pitch')
                        if pars.has_option('DUT{d}'.format(d=dut), 'first'):
                            self.dut_first_ch[dut] = pars.getint('DUT{d}'.format(d=dut), 'first')
                        if pars.has_option('DUT{d}'.format(d=dut), 'skip'):
                            self.dut_skip_ch[dut] = ChannelStringToArray(pars.get('DUT{d}'.format(d=dut), 'skip'))
                        if pars.has_option('DUT{d}'.format(d=dut), 'last'):
                            self.dut_last_ch[dut] = pars.getint('DUT{d}'.format(d=dut), 'last')

                if pars.has_section('ANALYSIS'):
                    if pars.has_option('ANALYSIS', 'do_pedestal'):
                        self.do_pedestal = pars.getboolean('ANALYSIS', 'do_pedestal')
                    if pars.has_option('ANALYSIS', 'do_cluster'):
                        self.do_cluster = pars.getboolean('ANALYSIS', 'do_cluster')
                    if pars.has_option('ANALYSIS', 'do_cluster_analysis'):
                        self.do_cluster_ana = pars.getboolean('ANALYSIS', 'do_cluster_analysis')
                    if pars.has_option('ANALYSIS', 'do_alignment'):
                        self.do_alignment = pars.getboolean('ANALYSIS', 'do_alignment')
                    if pars.has_option('ANALYSIS', 'do_alignment_analysis'):
                        self.do_alignment_ana = pars.getboolean('ANALYSIS', 'do_alignment_analysis')
                    if pars.has_option('ANALYSIS', 'do_transparent'):
                        self.do_transparent = pars.getboolean('ANALYSIS', 'do_transparent')
                    if pars.has_option('ANALYSIS', 'do_3d'):
                        self.do_3d = pars.getboolean('ANALYSIS', 'do_3d')

                if pars.has_section('PEDESTAL'):
                    if pars.has_option('PEDESTAL', 'tel_ped_hit_factor'):
                        self.tel_hit_factor = pars.getfloat('PEDESTAL', 'tel_ped_hit_factor')
                    if pars.has_option('PEDESTAL', 'dut_ped_hit_factor'):
                        self.dut_hit_factor = pars.getfloat('PEDESTAL', 'dut_ped_hit_factor')
                    if pars.has_option('PEDESTAL', 'do_cmc'):
                        self.do_cmc = pars.getboolean('PEDESTAL', 'do_cmc')
                    if pars.has_option('PEDESTAL', 'cm_cut'):
                        self.cm_cut = pars.getfloat('PEDESTAL', 'cm_cut')

                if pars.has_section('CLUSTER'):
                    if pars.has_option('CLUSTER', 'clust_seed_facts'):
                        self.clust_seed = eval(pars.get('CLUSTER', 'clust_seed_facts'))
                    if pars.has_option('CLUSTER', 'clust_hit_facts'):
                        self.clust_hit = eval(pars.get('CLUSTER', 'clust_hit_facts'))

                if pars.has_section('SELECTION_SCINT'):
                    if pars.has_option('SELECTION_SCINT', 'xlow'):
                        self.scint_fid_cut['xlow'] = pars.getfloat('SELECTION_SCINT', 'xlow')
                    if pars.has_option('SELECTION_SCINT', 'xhigh'):
                        self.scint_fid_cut['xhigh'] = pars.getfloat('SELECTION_SCINT', 'xhigh')
                    if pars.has_option('SELECTION_SCINT', 'ylow'):
                        self.scint_fid_cut['ylow'] = pars.getfloat('SELECTION_SCINT', 'ylow')
                    if pars.has_option('SELECTION_SCINT', 'yhigh'):
                        self.scint_fid_cut['yhigh'] = pars.getfloat('SELECTION_SCINT', 'yhigh')

                for dut in xrange(self.dut_num):
                    if pars.has_section('SELECTION{d}'.format(d=dut)):
                        if pars.has_option('SELECTION{d}'.format(d=dut), 'xlow') and pars.has_option('SELECTION{d}'.format(d=dut), 'xhigh') and pars.has_option('SELECTION{d}'.format(d=dut), 'ylow') and pars.has_option('SELECTION{d}'.format(d=dut), 'yhigh'):
                            self.fid_region[dut] = {'xlow': pars.getfloat('SELECTION{d}'.format(d=dut), 'xlow'), 'xhigh': pars.getfloat('SELECTION{d}'.format(d=dut), 'xhigh'), 'ylow': pars.getfloat('SELECTION{d}'.format(d=dut), 'ylow'), 'yhigh': pars.getfloat('SELECTION{d}'.format(d=dut), 'yhigh')}
                        elif pars.has_option('SELECTION{d}'.format(d=dut), 'xlow') or pars.has_option('SELECTION{d}'.format(d=dut), 'xhigh') or pars.has_option('SELECTION{d}'.format(d=dut), 'ylow') or pars.has_option('SELECTION{d}'.format(d=dut), 'yhigh'):
                            print 'setting default fiducial region for nos specifying all the parameters xlow, xhigh, ylow, yhigh'
                            self.fid_region[dut] = {'xlow': 0, 'xhigh': 255, 'ylow': 0, 'yhigh': 255}
                        else:
                            self.fid_region[dut] = {'xlow': 0, 'xhigh': 255, 'ylow': 0, 'yhigh': 255}
                    else:
                        ExitMessage('You should have a SELECTION section for each dut i.e. SELECTION0 and SELECTION1 if num or DUTS is 2')
                if pars.has_section('ALIGNMENT'):
                    if pars.has_option('ALIGNMENT', 'align_dut'):
                        self.align_dut = min(0, pars.getint('ALIGNMENT', 'align_dut'))
                    if pars.has_option('ALIGNMENT', 'z_coordinates'):
                        self.z_coordinates = eval(pars.get('ALIGNMENT', 'z_coordinates'))
                    if pars.has_option('ALIGNMENT', 'alignment_method'):
                        if pars.get('ALIGNMENT', 'alignment_method') in ['events', 'percentage']:
                            self.align_method = pars.get('ALIGNMENT', 'alignment_method')
                    if pars.has_option('ALIGNMENT', 'alignment_factor'):
                        self.align_factor = pars.getint('ALIGNMENT', 'alignment_factor')
                    if pars.has_option('ALIGNMENT', 'do_align_dut'):
                        self.do_align_dut = pars.getboolean('ALIGNMENT', 'do_align_dut')
                    if pars.has_option('ALIGNMENT', 'no_align_dut_chs'):
                        self.no_align_dut_chs = ChannelStringToArray(pars.get('ALIGNMENT', 'no_align_dut_chs'))
                    if pars.has_option('ALIGNMENT', 'alignment_chi2_cut'):
                        self.align_chi2_cut = pars.getfloat('ALIGNMENT', 'alignment_chi2_cut')

                if pars.has_section('TRANSPARENT'):
                    if pars.has_option('TRANSPARENT', 'max_transp_cluster_size'):
                        self.max_transp_clust_size = pars.getint('TRANSPARENT', 'max_transp_cluster_size')
                    if pars.has_option('TRANSPARENT', 'save_transp_cluster_size'):
                        self.save_transp_clust_size = pars.getint('TRANSPARENT', 'save_transp_cluster_size')
                    if pars.has_option('TRANSPARENT', 'analyse_align'):
                        self.do_analyse_alignment = pars.getboolean('TRANSPARENT', 'analyse_align')
        self.SetFileAndTreeName()

    def SetFileAndTreeName(self):
        self.tree_name = 'run_{r}_{s}'.format(r=self.run, s=self.sub_dir)
        self.file_name = 'run_{r}_{s}'.format(r=self.run, s=self.sub_dir)

if __name__ == '__main__':
    z = Settings()