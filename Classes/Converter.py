#!/usr/bin/env python
import ROOT as ro
import os, logging, shutil, sys
sys.path.append('/home/sandiego/Small3dAnalysis/Classes')
from Settings import Settings
import subprocess as subp
from Utils import *

__author__ = 'DA'

class Converter:
    def __init__(self, settings=Settings()):
        print 'Starting Converter...'
        self.settings = settings
        self.tree_name = self.settings.tree_name
        self.file_name = self.settings.file_name
        self.ana_events = self.settings.ana_events
        self.first_ev = self.settings.first_event
        self.num_parallel = self.settings.num_parallel
        self.out_dir = self.settings.output_dir
        self.sub_dir = self.settings.sub_dir
        self.ana_dir = self.settings.analysis_path
        self.run = self.settings.run
        self.do_conversion = self.Check_If_Already_Converted()
        self.num_events_per_job = int(self.ana_events) / int(self.num_parallel)
        self.event_saver_processes = []
        self.merge_process = None

    def Check_If_Already_Converted(self):
        if not os.path.isfile('{o}/{s}/{r}/{f}.root'.format(o=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name)):
            print 'The root file does not exist. Conversion started from event {ev0} until event {evf}.'.format(ev0=self.first_ev, evf=self.first_ev + self.ana_events - 1)
            return True
        else:
            tempf = ro.TFile('{o}/{s}/{r}/{f}.root'.format(o=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name))
            if tempf.IsZombie():
                print 'The existing root file cannot be opened. It will be re-written. Conversion started from event {ev0} until event {evf}.'.format(ev0=self.first_ev, evf=self.first_ev + self.ana_events - 1)
                tempf.Close()
                return True
            tempt = tempf.Get(self.tree_name)
            if not tempt:
                print 'The root file does not have the tree', self.tree_name, '. It will be created. Conversion started from event {ev0} until event {evf}.'.format(ev0=self.first_ev, evf=self.first_ev + self.ana_events - 1)
                tempf.Close()
                return True
            elif tempt.GetEntries() < self.ana_events:
                print 'The root file does not have enough events. It will be re-written. Conversion started from event {ev0} until event {evf}.'.format(ev0=self.first_ev, evf=self.first_ev + self.ana_events - 1)
                tempf.Close()
                return True
            elif tempt.GetEntries() == self.ana_events:
                print 'The root file has enough events. Skipping conversion'
                tempf.Close()
                return False
            else:
                print 'The root file has more events than requested. The requested events will be extracted from this file, assuming that the root file contains all the events measured during the run. No conversion will be done.'
                tempf.Close()
                self.ExtractEventsFromOriginalFile()
                return False

    def ExtractEventsFromOriginalFile(self):
        if os.path.isfile('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name)):
            tempf = ro.TFile('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), 'READ')
            tempt = tempf.Get(self.tree_name)
            print 'Extracting only', self.ana_events, 'events starting from', self.first_ev, 'for analysis...', ; sys.stdout.flush()
            leng = tempt.Draw('>>evlist', 'abs(2*event-{eva}-2*{evi}+1)<=({eva}-1)'.format(evi=self.first_ev, eva=self.ana_events))
            while leng > tempt.GetEstimate():
                tempt.SetEstimate(leng)
                leng = tempt.Draw('>>evlist', 'abs(2*event-{eva}-2*{evi}+1)<=({eva}-1)'.format(evi=self.first_ev, eva=self.ana_events))
            evlist = ro.gDirectory.Get('evlist')
            tempt.SetEventList(evlist)
            tempnf = ro.TFile('{d}/{s}/{r}/{f}_new.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), 'RECREATE')
            tempnt = tempt.CopyTree('')
            tempnt.Write()
            tempnf.Close()
            tempf.Close()
            print 'Done'
            print 'Moving original file as ....root.bkp ...', ; sys.stdout.flush()
            shutil.move('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), '{d}/{s}/{r}/{f}.original.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name))
            shutil.move('{d}/{s}/{r}/{f}_new.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), '{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name))
            print 'Done'
        else:
            ExitMessage('File {f} does not exist. Quitting...'.format(f='{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name)))

    def Convert(self):
        with open(os.devnull, 'w') as FNULL:
            for job_i in xrange(self.num_parallel):
                if job_i != self.num_parallel - 1:
                    self.event_saver_processes.append(subp.Popen(['{ad}/RawEventSaver.py'.format(ad=self.ana_dir), '-s', '{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), '-j', str(job_i)], bufsize=-1, stdin=subp.PIPE, stdout=FNULL, close_fds=True))
                else:
                    self.event_saver_processes.append(subp.Popen(['{ad}/RawEventSaver.py'.format(ad=self.ana_dir), '-s', '{d}/{s}/{r}/{f}.settings'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name), '-j', str(job_i), '-p'], bufsize=-1, stdin=subp.PIPE, close_fds=True))
            for job_i in xrange(self.num_parallel):
                while self.event_saver_processes[job_i].poll() is None:
                    continue
                if job_i != self.num_parallel - 1:
                    CloseSubprocess(self.event_saver_processes[job_i], stdin=True, stdout=False)
                else:
                    CloseSubprocess(self.event_saver_processes[job_i], stdin=True, stdout=False)
                print 'Done with events:', self.first_ev + job_i * self.num_events_per_job, '-', self.first_ev + (job_i + 1) * self.num_events_per_job - 1
        self.MergeOutputFiles()

    def MergeOutputFiles(self):
        command = ['hadd', '{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name)]
        for job_i in xrange(self.num_parallel):
            if os.path.isfile('{d}/{s}/{r}/{f}_{j}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name, j=job_i)):
                command.append('{d}/{s}/{r}/{f}_{j}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name, j=job_i))
            else:
                ExitMessage('The file: {d}/{s}/{r}/{f}_{j}.root does not exist! Exiting...'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name, j=job_i), os.EX_OSFILE)
        if os.path.isfile('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name)):
            os.remove('{d}/{s}/{r}/{f}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name))
        print 'Merging files with the following command:', command, '...', ; sys.stdout.flush()
        self.merge_process = subp.Popen(command, bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)
        while self.merge_process.poll() is None:
            continue
        print 'Done'
        CloseSubprocess(self.merge_process, stdin=True, stdout=True)
        self.DeleteUnMergedFiles()

    def DeleteUnMergedFiles(self):
        for file_i in ['{d}/{s}/{r}/{f}_{j}.root'.format(d=self.out_dir, s=self.sub_dir, r=self.run, f=self.file_name, j=job_i) for job_i in xrange(self.num_parallel)]:
            print 'Deleting partial file', file_i, '...', ; sys.stdout.flush()
            if os.path.isfile(file_i):
                os.remove(file_i)
            print 'Done'


def main():
    conv = Converter()

if __name__ == '__main__':
    main()
