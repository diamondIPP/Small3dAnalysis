from ROOT import TFile, gSystem, TStopwatch, TDatime
from optparse import OptionParser
from time import time
from Settings import Settings
from numpy import array, zeros
from copy import deepcopy
import os, logging
import struct
from numba import jit, autojit, jitclass
import ipdb

__author__ = 'DA'

class RawEventReader:
    def __init__(self, settings=None):
        print 'Starting Raw event reader'
        self.eventsPerFile = 10000
        self.silTelMem = 2048
        self.diaMem = 256
        self.silTelChs = 256
        self.diaChs = 128
        self.settings = settings
        self.run = self.settings.runInfo['run']
        self.current_rz_filename = ''
        self.current_rz_file = None
        self.struct_fmt = '>IIIiiiII8i8iihh2i2048B256HI'  # big-endian
        self.struct_len = struct.calcsize(self.struct_fmt)
        self.D0X, self.D0Y, self.D1X, self.D1Y = zeros(self.silTelChs, 'B'), zeros(self.silTelChs, 'B'), zeros(self.silTelChs, 'B'), zeros(self.silTelChs, 'B')
        self.D2X, self.D2Y, self.D3X, self.D3Y = zeros(self.silTelChs, 'B'), zeros(self.silTelChs, 'B'), zeros(self.silTelChs, 'B'), zeros(self.silTelChs, 'B')
        self.Dia0, self.Dia1 = zeros(self.diaChs, 'H'), zeros(self.diaChs, 'H')
        self.silTel, self.dia = None, None

    def ReadRawEvent(self, eventNo):
        filename = '{aip}/RUN_{r}_{p}.rz'.format(aip=self.settings.GetAbsoluteInputPath(), r=self.run, p=int(eventNo/self.eventsPerFile))
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
        # res, self.silTel, self.dia, eor = array(s[27:29], 'i'), array(s[29:29+self.silTelMem], 'B'), array(s[29+self.silTelMem:29+self.silTelMem+self.diaMem], 'H'), s[29+self.silTelMem+self.diaMem]
        self.silTel, self.dia = array(s[29:29+self.silTelMem], 'B'), array(s[29+self.silTelMem:29+self.silTelMem+self.diaMem], 'H')
        self.DecodeDetectors()
        return True

    def DecodeDetectors(self):
        self.DecodeSiliconTelescope()
        self.DecodeDiamondDetectors()

    # @jit
    def DecodeSiliconTelescope(self):
        for i in xrange(self.silTelChs):
            self.D0X[i] = self.silTel[i*4+3]
            self.D1X[i] = self.silTel[i*4+2]
            self.D2X[i] = self.silTel[i*4+1]
            self.D3X[i] = self.silTel[i*4]
            self.D0Y[i] = self.silTel[(self.silTelChs + i)*4 + 3]
            self.D1Y[i] = self.silTel[(self.silTelChs + i)*4 + 2]
            self.D2Y[i] = self.silTel[(self.silTelChs + i)*4 + 1]
            self.D3Y[i] = self.silTel[(self.silTelChs + i)*4]

    def DecodeDiamondDetectors(self):
        for i in xrange(self.diaChs):
            self.Dia0[i] = self.dia[i*2+1]
            self.Dia1[i] = self.dia[i*2]

    def GetSilPlane(self, det):
        return {0: self.D0X, 1: self.D0Y, 2: self.D1X, 3: self.D1Y, 4: self.D2X, 5: self.D2Y, 6: self.D3X, 7: self.D3Y}[det]

    def GetDiaDet(self, det):
        return {0: self.Dia0, 1: self.Dia1}[det]

if __name__ == '__main__':
    st = time()
    z = RawEventReader()