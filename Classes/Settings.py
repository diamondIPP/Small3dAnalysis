from ROOT import gSystem, TFile, TTree
from copy import deepcopy
import os
import progressbar
from ConfigParser import ConfigParser
from collections import namedtuple
import ipdb

class Settings:
    def __init__(self, runInfo=None, current='./', input='./', output='./', settings='./'):
        self.runInfo = runInfo
        self.currentDir = current
        self.inputDir = input
        self.outputDir = output
        self.settingsDir = settings
        self.CheckDirsExistence()
        self.LoadSettings()
        self.nEvents = self.runInfo['nEvents']
        if self.nEvents > self.totalEvents:
            print 'The run only has', self.totalEvents, 'events. Analysing only', self.totalEvents, 'events instead of', self.nEvents,'.'
            self.nEvents = self.totalEvents
        self.previousDir = self.currentDir
        self.bar = None
        # TODO do self.LoadSettingsFile() that calls self.LoadDefaultSettings() when it does not exist and delete following lines:
        self.diaInput = 0  # TODO: delete!!!
        self.version = '0'
        self.silDetChs = 256
        self.diaDetChs = 128

    def CheckDirsExistence(self):
        self.CheckDirExistence(self.inputDir)
        self.CheckDirExistence(self.outputDir, True)
        self.CheckDirExistence(self.settingsDir)

    def CheckDirExistence(self, directory, create=False):
        if not os.path.isdir(directory):
            print 'Dir: {d} does not exist!'.format(d=directory)
            if create:
                print 'Trying to create: {d}'.format(d=directory)
                if gSystem.mkdir(directory, True) == -1:
                    print 'Could not create it'
                    exit()
                print 'Successfully created dir: {d}'.format(d=directory)
            else:
                exit()

    def CheckFileExistence(self, file, create=False):
        if not os.path.isfile(file):
            print 'File: {f} does not exist!'.format(f=file)
            if create:
                print 'Trying to create: {f}'.format(f=file)
                temp = TFile(file, )

    def CreateProgressBar(self, maxVal=0):
        widgets = [
            'Processed: ', progressbar.Counter(),
            ' of {mv} '.format(mv=maxVal), progressbar.Percentage(),
            ' ', progressbar.Bar(marker='>'),
            ' ', progressbar.Timer(),
            ' ', progressbar.ETA()
            # ' ', progressbar.AdaptativeETA(),
            #  ' ', progressbar.AdaptativeTransferSpeed()
            ]
        self.bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)

    def GetAbsoluteOutputPath(self):
        string = self.outputDir + '/' + str(self.runInfo['run']) + '/'
        return string

    def GoToDir(self, directory):
        self.CheckDirExistence(directory, True)
        gSystem.cd(directory)
        self.currentDir = gSystem.pwd()

    def GoToOutputPath(self):
        self.previousDir = self.currentDir
        self.GoToDir(self.GetAbsoluteOutputPath())

    def GoToPreviousDir(self):
        temp = self.currentDir
        self.GoToDir(self.previousDir)
        self.previousDir = temp
        self.currentDir = gSystem.pwd()

    def GetRawTreeFilePath(self):
        path = self.GetAbsoluteOutputPath() + 'rawData.{r}.root'.format(r=self.runInfo['run'])
        return path

    def GetAbsoluteInputPath(self):
        string = self.inputDir + '/' + str(self.runInfo['run']) + '/'
        return string

    def GetPedestalTreeFilePath(self):
        path = self.GetAbsoluteOutputPath() + 'pedestalData.{r}.root'.format(r=self.runInfo['run'])
        return path

    def GetEtaIntegralFilePath(self):
        path = self.GetAbsoluteOutputPath() + 'etaCorrection.{r}.root'.format(r=self.runInfo['run'])
        return path

    def LoadSettings(self):
        self.LoadDefaults()
        if os.path.isfile(self.settingsDir + '/settings.{r}.ini'.format(r=self.runInfo['run'])):
            parser = ConfigParser()
            parser.read(self.settingsDir + '/settings.{r}.ini'.format(r=self.runInfo['run']))
            print 'Loading settings from file: {s}/settings.{r}.ini'.format(s=self.settingsDir, r=self.runInfo['run'])
            self.totalEvents = parser.getint('BASIC', 'Events') if parser.has_option('BASIC', 'Events') else self.totalEvents
            self.repeaterCard = parser.getint('BASIC', 'repeaterCardNo') if parser.has_option('BASIC', 'repeaterCardNo') else self.repeaterCard
            self.voltage = parser.getint('BASIC', 'voltage') if parser.has_option('BASIC', 'voltage') else self.voltage
            self.diamondName = parser.get('BASIC', 'diamondName') if parser.has_option('BASIC', 'diamondName') else self.diamondName
            self.currentBegin = parser.getfloat('BASIC', 'currentBegin') if parser.has_option('BASIC', 'currentBegin') else self.currentBegin
            self.currentEnd = parser.getfloat('BASIC', 'currentEnd') if parser.has_option('BASIC', 'currentEnd') else self.currentEnd
            self.diaInput = parser.getint('BASIC', 'dia_input') if parser.has_option('BASIC', 'dia_input') else self.diaInput
            self.SetDiamondPatterns(parser)
            self.SetScreenedChannels(parser)
            self.SetNoisyChannels(parser)
            self.SetNotConnectedChannels(parser)
            self.siPedestalHitFactor = parser.getint('PEDESTAL', 'SiPedestalHitFactor') if parser.has_option('PEDESTAL', 'SiPedestalHitFactor') else self.siPedestalHitFactor
            self.diaPedestalHitFactor = parser.getint('PEDESTAL', 'DiaPedestalHitFactor') if parser.has_option('PEDESTAL', 'DiaPedestalHitFactor') else self.diaPedestalHitFactor
            self.doCMC = parser.getboolean('PEDESTAL', 'doCMC') if parser.has_option('PEDESTAL', 'doCMC') else self.doCMC
            self.cmnCut = parser.getint('PEDESTAL', 'cmnCut') if parser.has_option('PEDESTAL', 'cmnCut') else self.cmnCut
            self.bufferSize = parser.getint('PEDESTAL', 'Iter_Size') if parser.has_option('PEDESTAL', 'Iter_Size') else self.bufferSize
            if parser.has_option('CLUSTER', 'clusterSeedFactors'):
                exec('self.clusterSeedFactors = {v}'.format(v=parser.get('CLUSTER', 'clusterSeedFactors')))
            if parser.has_option('CLUSTER', 'clusterHitFactors'):
                exec('self.clusterHitFactors = {v}'.format(v=parser.get('CLUSTER', 'clusterHitFactors')))
            self.nDia = parser.getint('SELECTION', 'nDia') if parser.has_option('SELECTION', 'nDia') else self.nDia
            self.alignDia = parser.getint('SELECTION', 'alignDia') if parser.has_option('SELECTION', 'alignDia') else self.alignDia
            self.threeDDia = parser.getint('SELECTION', '3DDia') if parser.has_option('SELECTION', '3DDia') else self.threeDDia
            self.phantomDia = parser.getint('SELECTION', 'phantomDia') if parser.has_option('SELECTION', 'phantomDia') else self.phantomDia
            self.SetScintillatorFidCut(parser)
            self.SetDiamondFidCut(parser)
            self.SetDiamondSelections(parser)
            if parser.has_option('ALIGNMENT', 'zCoordinates'):
                exec('self.zCoordinates = {v}'.format(v=parser.get('ALIGNMENT', 'zCoordinates')))
            self.alignMethod = parser.get('ALIGNMENT', 'alignmentMethod') if parser.has_option('BASIC', 'alignmentMethod') else self.alignMethod
            self.doAlignDia = parser.getboolean('ALIGNMENT', 'doAlignDiamond') if parser.has_option('ALIGNMENT', 'doAlignDiamond') else self.doAlignDia
            self.SetAlignIgnoreDiaChs(parser)
            self.alignmentChi2Cut = parser.getint('ALIGNMENT', 'alignmentChi2Cut') if parser.has_option('ALIGNMENT', 'alignmentChi2Cut') else self.alignmentChi2Cut
            self.alignmentFactor = parser.getfloat('ALIGNMENT', 'alignmentFactor') if parser.has_option('ALIGNMENT', 'alignmentFactor') else self.alignmentFactor
        else:
            print 'Loaded default settings'

    def SetDiamondSelections(self, parser):
        if parser.has_option('SELECTION', 'selectionFidCut'):
            dic = {}
            exec('dic = {v}'.format(v=parser.get('SELECTION', 'selectionFidCut')))
            self.selectionFidCut = {key: {'xlow': value['x'][0], 'xhigh': value['x'][-1], 'ylow': value['y'][0], 'yhigh': value['y'][-1]} for key, value in dic.iteritems()}

    def SetScintillatorFidCut(self, parser):
        if parser.has_option('SELECTION', 'ScintillatorFidCutX') and parser.has_option('SELECTION', 'ScintillatorFidCutY'):
            listX, listY = [], []
            exec('listX = {v}'.format(v=parser.get('SELECTION', 'ScintillatorFidCutX')))
            exec('listY = {v}'.format(v=parser.get('SELECTION', 'ScintillatorFidCutY')))
            self.scintFidCut = {'xlow': listX[0], 'xhigh': listX[-1], 'ylow': listY[0], 'yhigh': listY[-1]}

    def SetDiamondFidCut(self, parser):
        if parser.has_option('SELECTION', 'DiaFidCutX') and parser.has_option('SELECTION', 'DiaFidCutY'):
            listX, listY = [], []
            exec('listX = {v}'.format(v=parser.get('SELECTION', 'DiaFidCutX')))
            exec('listY = {v}'.format(v=parser.get('SELECTION', 'DiaFidCutY')))
            self.diaFidCut = {'xlow': listX[0], 'xhigh': listX[-1], 'ylow': listY[0], 'yhigh': listY[-1]}

    def SetDiamondPatterns(self, parser):
        if parser.has_option('DIAMOND', 'diamondPattern'):
            dic = {}
            exec('dic = {v}'.format(v=parser.get('DIAMOND', 'diamondPattern')))
            self.diaPattern = {key: {'x0': value[0], 'pitch': value[1], 'first': value[2], 'last': value[3]} for key, value in dic.iteritems()}

    def SetAlignIgnoreDiaChs(self, parser):
        if parser.has_option('ALIGNMENT', 'alignIgnoreDiaChs'):
            string = parser.get('ALIGNMENT', 'alignIgnoreDiaChs')
            self.alignIgnoreDiaChs = self.GetChannelsFromString(string)

    def SetScreenedChannels(self, parser):
        if parser.has_option('DIAMOND', 'diaChannelsScreened'):
            string = parser.get('DIAMOND', 'diaChannelsScreened')
            self.diaChsScreened = self.GetChannelsFromString(string)

    def SetNoisyChannels(self, parser):
        if parser.has_option('DIAMOND', 'diaNoisyChannels'):
            string = parser.get('DIAMOND', 'diaNoisyChannels')
            self.diaChsNoisy = self.GetChannelsFromString(string)

    def SetNotConnectedChannels(self, parser):
        if parser.has_option('DIAMOND', 'diaChannelsNotConnected'):
            string = parser.get('DIAMOND', 'diaChannelsNotConnected')
            self.diaChsNotConnected = self.GetChannelsFromString(string)

    def GetChannelsFromString(self, string):
        string.strip()
        string = string.replace(' ', '')
        strips = string.split(',')
        listChs = [range(int(i[0]), int(i[-1]) + 1) for i in [elem.split(':') for elem in strips]]
        return deepcopy([ch for group in listChs for ch in group])

    def OpenFile(self, filePath='', create=True):
        created = False
        if create:
            tfile = TFile(filePath, 'NEW')
            if tfile.IsOpen():
                created = True
                print 'Created new root file:', filePath
            else:
                tfile = TFile(filePath, 'READ')
                print 'Opened existing file:', filePath
        else:
            tfile = TFile(filePath, 'READ')
        fileBool = namedtuple('File_bool', ['file', 'boolFile'])
        fileBoolReturn = fileBool(tfile, created)
        return fileBoolReturn

    def OpenTree(self, filePath='', treeName='tree', create=True):
        if create:
            createdNewFile = False
            tfile = TFile(filePath, 'NEW')
            if tfile.IsOpen():
                createdNewFile = True
                print 'Created new root file:', filePath
            else:
                tfile= TFile(filePath, 'READ')
                print 'Opened existing file:', filePath
            tree = tfile.Get(treeName)
            createdNewTree = False
            if not tree:
                if not createdNewFile:
                    print 'File had no Tree. Recreating file'
                    tfile.Close()
                    tfile= TFile(filePath, 'RECREATE')
                    createdNewFile = True
                tree = TTree(treeName, '{t} data of run {r}'.format(t=treeName, r=self.runInfo['run']))
                createdNewTree = True
            else:
                if tree.GetEntries() < self.nEvents:
                    print 'Tree exists but does not have the desired number of events. Recreating'
                    tree.Delete()
                    tree = TTree(treeName, '{t} data of run {r}'.format(t=treeName, r=self.runInfo['run']))
                    createdNewTree = True
        else:
            tfile= TFile(filePath, 'READ')
            tree = tfile.Get(treeName)
            createdNewFile, createdNewTree = False, False
        fileTree = namedtuple('File_Tree_data', ['file', 'tree', 'boolFile', 'boolTree'])
        fileTreeReturn = fileTree(tfile, tree, createdNewFile, createdNewTree)
        return fileTreeReturn

    def IsMasked(self, type='dia', det=0, ch=0):
        if type == 'dia':
            return ch in self.diaChsScreened
        else:
            return False

    def IsConnected(self, ch):
        return ch not in self.diaChsNotConnected

    def LoadDefaults(self):
        self.totalEvents = 100000
        self.repeaterCard = 2
        self.voltage = 0
        self.diamondName = 'diamond'
        self.currentBegin = -1
        self.currentEnd = -1
        self.diaInput = 0
        self.diaPattern = {0: {'x0': 25, 'pitch': 50, 'first': 1, 'last': 126}}
        self.diaChsScreened = [0, 127]
        self.diaChsNoisy = []
        self.diaChsNotConnected = [0, 127]
        self.siPedestalHitFactor = 5
        self.diaPedestalHitFactor = 3
        self.doCMC = True
        self.cmnCut = 4
        self.bufferSize = 500
        self.clusterSeedFactors = {0: 14, 1: 19, 2: 23, 3: 23, 4: 14, 5: 13, 6: 11, 7: 10, 8: 5}
        self.clusterHitFactors = {0: 12, 1: 14, 2: 17, 3: 17, 4: 8, 5: 7, 6: 7, 7: 6, 8: 3}
        self.nDia = 1
        self.alignDia = 0
        self.threeDDia = -1
        self.phantomDia = -1
        self.scintFidCut = {'xlow': 1, 'xhigh': 254, 'ylow': 1, 'yhigh': 254}
        self.diaFidCut = {'xlow': 1, 'xhigh': 254, 'ylow': 1, 'yhigh': 254}
        self.selectionFidCut = {0: {'xlow': 1, 'xhigh': 254, 'ylow': 1, 'yhigh': 254}}
        self.zCoordinates = {0: 0.725, 1: 0.725, 2: 1.625, 3: 1.625, 4: 18.725, 5: 18.725, 6: 19.625, 7: 19.625, 8: 10.2}
        self.alignMethod = 'events'
        self.doAlignDia = True
        self.alignIgnoreDiaChs = [0, 127]
        self.alignmentChi2Cut = 4
        self.alignmentFactor = 10000
        self.silNumPlanes = 8
        self.diaNumPlanes = 1
        self.diaMaxADC = 4095
        self.silMaxADC = 255
        self.silNumPlanes = 4
        self.silNumDetectors = 8
        self.maxTranspCluster = 10

if __name__ == '__main__':
    z = Settings()