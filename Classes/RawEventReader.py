#!/usr/bin/env python
import os, logging, sys
sys.path.append('/home/sandiego/Small3dAnalysis/Classes')  # TODO: HARDCODED!!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
from Settings import Settings
import numpy as np
import struct

__author__ = 'DA'

class RawEventReader:
    def __init__(self, settings=Settings()):
        print 'Starting Raw event reader'
        self.settings = settings
        self.raw_data_dir = self.settings.data_dir
        self.eventsPerFile = self.settings.event_per_file
        self.telMem = self.settings.tel_mem
        self.dutMem = self.settings.dut_mem
        self.telDetChs = self.settings.telDetChs
        self.dutDetChs = self.settings.dutDetChs
        self.run = self.settings.run
        self.current_rz_filename = ''
        self.current_rz_file = None
        self.struct_fmt = self.settings.struct_fmt
        self.struct_len = struct.calcsize(self.struct_fmt)
        self.D0X, self.D0Y, self.D1X, self.D1Y = np.zeros(self.telDetChs, 'B'), np.zeros(self.telDetChs, 'B'), np.zeros(self.telDetChs, 'B'), np.zeros(self.telDetChs, 'B')
        self.D2X, self.D2Y, self.D3X, self.D3Y = np.zeros(self.telDetChs, 'B'), np.zeros(self.telDetChs, 'B'), np.zeros(self.telDetChs, 'B'), np.zeros(self.telDetChs, 'B')
        self.Dut0, self.Dut1 = np.zeros(self.dutDetChs, 'H'), np.zeros(self.dutDetChs, 'H')
        self.tel, self.dut = None, None

    def ReadRawEvent(self, eventNo):
        filename = '{aip}/RUN_{r}_{p}.rz'.format(aip=self.raw_data_dir, r=self.run, p=int(eventNo/self.eventsPerFile))
        if self.current_rz_filename != filename:
            if self.current_rz_file:
                self.current_rz_file.close()
            self.current_rz_filename = filename
            try:
                self.current_rz_file = open(filename, 'rb')  # read and binary mode
            except IOError:
                print 'Error: File does not appear to exist'
                return False
        self.current_rz_file.seek((int(eventNo) % int(self.eventsPerFile)) * self.struct_len)
        data = self.current_rz_file.read(self.struct_len)
        if not data:
            print 'File currupt or event not recorded. Try again with less events than', str(eventNo + 1)
            return False
        s = struct.Struct(self.struct_fmt).unpack_from(data)
        # evtrig, evNo, evPos, evTag, evDate = s[0], s[1], s[2], s[3], s[4]
        # evTime, trigCnt, evVmeTime, vFasCnt = s[5], s[6], s[7], array(s[8:16], 'i')
        # vFasReg, evNetTime, measNo, evInMeasNo = array(s[16:24], 'i'), s[24], s[25], s[26]
        # res, self.tel, self.dut, eor = array(s[27:29], 'i'), array(s[29:29+self.silTelMem], 'B'), array(s[29+self.silTelMem:29+self.silTelMem+self.diaMem], 'H'), s[29+self.silTelMem+self.diaMem]
        # self.tel, self.dut = array(s[29:29+self.telMem], 'B'), array(s[29+self.telMem:29+self.telMem+self.dutMem], 'H')
        self.tel, self.dut = s[29:29 + self.telMem], s[29 + self.telMem:29 + self.telMem + self.dutMem]
        self.DecodeDetectors()
        return True

    def DecodeDetectors(self):
        self.DecodeTelescope()
        self.DecodeDutDetectors()

    def DecodeTelescope(self):
        self.D0X = [self.tel[i * 4 + 3] for i in xrange(self.telDetChs)]
        self.D1X = [self.tel[i * 4 + 2] for i in xrange(self.telDetChs)]
        self.D2X = [self.tel[i * 4 + 1] for i in xrange(self.telDetChs)]
        self.D3X = [self.tel[i * 4] for i in xrange(self.telDetChs)]
        self.D0Y = [self.tel[(self.telDetChs + i) * 4 + 3] for i in xrange(self.telDetChs)]
        self.D1Y = [self.tel[(self.telDetChs + i) * 4 + 2] for i in xrange(self.telDetChs)]
        self.D2Y = [self.tel[(self.telDetChs + i) * 4 + 1] for i in xrange(self.telDetChs)]
        self.D3Y = [self.tel[(self.telDetChs + i) * 4] for i in xrange(self.telDetChs)]

    def DecodeDutDetectors(self):
        self.Dut0 = [self.dut[i * 2 + 1] for i in xrange(self.dutDetChs)]
        self.Dut1 = [self.dut[i * 2] for i in xrange(self.dutDetChs)]

    def GetTelPlane(self, det):
        return {0: self.D0X, 1: self.D0Y, 2: self.D1X, 3: self.D1Y, 4: self.D2X, 5: self.D2Y, 6: self.D3X, 7: self.D3Y}[det]

    def GetDutDet(self, det):
        return {0: self.Dut0, 1: self.Dut1}[det]


if __name__ == '__main__':
    z = RawEventReader()