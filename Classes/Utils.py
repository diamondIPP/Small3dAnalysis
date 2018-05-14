#!/usr/bin/env python
import os, shutil, sys
import progressbar
import ROOT as ro
import numpy as np

class Utils:
	def __init__(self):
		self.bar = None

	def CreateProgressBar(self, maxVal=1):
		widgets = [
			'Processed: ', progressbar.Counter(),
			' out of {mv} '.format(mv=maxVal), progressbar.Percentage(),
			' ', progressbar.Bar(marker='>'),
			' ', progressbar.Timer(),
			' ', progressbar.ETA()
			# ' ', progressbar.AdaptativeETA(),
			#  ' ', progressbar.AdaptativeTransferSpeed()
		]
		self.bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)

def ChannelStringToArray(chstr):
	arr = []
	parts = chstr.split(',')
	for part in parts:
		sub_part = part.split(':')
		len_sp = len(sub_part)
		if 1 <= len_sp <= 2:
			arr.extend(range(int(sub_part[0]), int(sub_part[len_sp-1])+1))
	return sorted(list(set(arr)))

def CreateDirectoryIfNecessary(dir):
	if not os.path.isdir(dir):
		print 'Creating directory', dir, '...', ; sys.stdout.flush()
		os.makedirs(dir)
		print 'Done'

def ExitMessage(txt='', value=os.EX_SOFTWARE):
	if txt != '': print txt
	sys.exit(value)

def CloseSubprocess(p, stdin=False, stdout=False):
	pid = p.pid
	if stdin:
		p.stdin.close()
	if stdout:
		p.stdout.close()
	if p.wait() is None:
		print 'Could not terminate subprocess... forcing termination'
		p.kill()
		if p.wait() is None:
			ExitMessage('Could not kill subprocess... quitting')
	try:
		os.kill(pid, 0)
	except OSError:
		pass
	else:
		print 'The subprocess is still running. Killing it with os.kill'
		os.kill(pid, 15)
		try:
			os.kill(pid, 0)
		except OSError:
			pass
		else:
			ExitMessage('The process does not die... quitting program')
	del p, pid
	p = None

def CheckBranchExistence(file_name, tree_name, branch):
	if not os.path.isfile(file_name):
		ExitMessage('The root file {f} does not exist. Needs conversion. Exiting...'.format(f=file_name), os.EX_OSFILE)
	else:
		tempf = ro.TFile(file_name, 'READ')
		if tempf.IsZombie():
			tempf.Close()
			ExitMessage('The root file {f} cannot be opened. Exiting'.format(f=file_name), os.EX_OSFILE)
		tempt = tempf.Get(tree_name)
		if not tempt:
			ExitMessage('The root file {f} does not have the tree {t}. Exiting'.format(f=file_name, t=tree_name))
		elif not tempt.GetLeaf(branch):
			print 'The tree', tree_name, 'does not have the branch', branch
			return False
		else:
			print 'The tree', tree_name, 'has the branch', branch
			return True

def Open_RootFile_Load_Tree(path='', treename='', mode='read'):
	file_f = ro.TFile(path, mode)
	if file_f.IsZombie():
		ExitMessage('The file {f} could not be opened. Exiting...'.format(f=path), os.EX_OSFILE)
	tree_t = Get_Tree_From_TFile(file_f, treename)
	return file_f, tree_t

def Get_Tree_From_TFile(file_f, treename=''):
	tree_t = file_f.Get(treename)
	if not tree_t:
		ExitMessage('The file does not have the tree {t}. Exiting...'.format(t=treename))
	return tree_t

def Select_Branches_For_Get_Val(tree_r, list_bra=[]):
	tree_r.SetBranchStatus('*', 0)
	for bra in list_bra:
		if '[' in bra:
			tree_r.SetBranchStatus(bra.split('[')[0], 1)
		else:
			tree_r.SetBranchStatus(bra, 1)

def Draw_Branches_For_Get_Val(tree_r, list_bra=[], start_ev=0, n_entries=1, cut='', option='goff'):
	draw_opt = option + ' para' if len(list_bra) > 3 and 'para' not in option else option
	draw_opt = draw_opt + ' goff' if 'goff' not in draw_opt else draw_opt
	Select_Branches_For_Get_Val(tree_r, list_bra)
	leng = tree_r.Draw(':'.join(list_bra), cut, draw_opt, n_entries, start_ev)
	while leng > tree_r.GetEstimate():
		tree_r.SetEstimate(leng)
		leng = tree_r.Draw(':'.join(list_bra), cut, draw_opt, n_entries, start_ev)

def Get_Branches_Value_To_Numpy(tree_r, list_bra=[], list_numpy_arrays=[], n_entries=1, elements_per_ev=1):
	for val, arr in enumerate(list_numpy_arrays):
		Branch_List_To_Array(tree_r.GetVal(val), arr, n_entries, elements_per_ev)

def Branch_List_To_Array(bra_list, array_v=np.zeros((128,500)), n_entries=1, vec_size=1):
	if vec_size == 1:
		np.putmask(array_v, np.ones(n_entries, dtype='?'), [bra_list[ev] for ev in xrange(n_entries)])
	else:
		np.putmask(array_v, np.ones((vec_size, n_entries), dtype='?'), [[bra_list[ev * vec_size + ch] for ev in xrange(n_entries)] for ch in xrange(vec_size)])

def Read_Single_Value_Event(tree_t, event=0, branch='event'):
	tree_t.GetEntry(event)
	value = eval('tree_t.{b}'.format(t=tree_t, b=branch))
	if type(value) == type('a'):
		return ord(value)
	return value

# def Insert_Value_In_Buffer_Array(array_v=np.zeros(500), val=0, slide_leng=500):
# 	np.putmask(array_v, np.ones(slide_leng, '?'), np.roll(array_v, -1))
# 	array_v[-1] = val

def Insert_Value_In_Buffer_Array(array_v=np.zeros(500), val=0):
	array_v = np.roll(array_v, -1)
	array_v[-1] = val
	return array_v

# def IsInt(i):
# 	try:
# 		int(i)
# 		return True
# 	except ValueError:
# 		return False
#
# def IsFloat(f):
# 	try:
# 		float(f)
# 		return True
# 	except ValueError:
# 		return False
#
# def Set_voltage(s_file):
# 	Replace_Settings_Line(s_file, 'voltage')
#
# def Mark_channels_as_NC(s_file):
# 	Replace_Settings_Line(s_file, 'Dia_channel_not_connected')
#
# def Mark_channels_as_Screened(s_file):
# 	Replace_Settings_Line(s_file, 'Dia_channel_screen_channels')
#
# def Mark_channels_as_Noisy(s_file):
# 	Replace_Settings_Line(s_file, 'Dia_channel_noisy')
#
# def Set_Diamond_pattern(s_file):
# 	Replace_Settings_Line(s_file, 'diamondPattern')
#
# def Set_Diamond_name(s_file):
# 	Replace_Settings_Line(s_file, 'diamondName')
#
# def Mark_channels_as_no_alignment(s_file):
# 	with open(s_file, 'r') as f0:
# 		for line in f0:
# 			if line.lower().startswith('Dia_channel_not_connected'.lower()):
# 				print 'The current value(s) for not connected channels are:'
# 				print line.split('=')[1].split('\n')[0]
# 			if line.lower().startswith('Dia_channel_screen_channels'.lower()):
# 				print 'The current value(s) for screened channels are:'
# 				print line.split('=')[1].split('\n')[0]
# 			if line.lower().startswith('Dia_channel_noisy'.lower()):
# 				print 'The current value(s) for noisy channels are:'
# 				print line.split('=')[1].split('\n')[0]
# 	Replace_Settings_Line(s_file, 'AlignmentIgnoreChannelDia')
#
# def Select_fiducial_region(s_file):
# 	xlow = Replace_Settings_Line(s_file, 'si_avg_fidcut_xlow')
# 	xhigh = Replace_Settings_Line(s_file, 'si_avg_fidcut_xhigh')
# 	ylow = Replace_Settings_Line(s_file, 'si_avg_fidcut_ylow')
# 	yhigh = Replace_Settings_Line(s_file, 'si_avg_fidcut_yhigh')
# 	num_patt = 0
# 	with open(s_file, 'r') as f0:
# 		for line in f0:
# 			if line.lower().startswith('diamondPattern'.lower()):
# 				num_patt += 1
# 	if num_patt == 1:
# 		Replace_Settings_Line(s_file, 'selectionFidCut', action='value', value=('{' + xlow + '-' + xhigh + ',' + ylow + '-' + yhigh + '}').replace(' ', ''))
#
# def Replace_Settings_Line(s_file, setting_option, action='user', value=''):
# 	if not os.path.isfile(s_file):
# 		print 'The file', s_file, 'does not exist'
# 		return
# 	updated = False
# 	with open(s_file, 'r') as f0:
# 		with open('temp00.ini', 'w') as f1:
# 			for line in f0:
# 				if not line.lower().startswith(setting_option.lower()):
# 					f1.write(line)
# 				else:
# 					if action == 'user':
# 						new_values = Get_From_User_Line(setting_option, line.split('=')[1].split('\n')[0])
# 					elif action == 'even':
# 						new_values = Only_Even_Channels(line)
# 					elif action == 'odd':
# 						new_values = Only_Odd_Channels(line)
# 					elif action == 'value':
# 						new_values = value
# 					f1.write(setting_option + ' = ' + new_values + '\n')
# 					updated = True
# 	if not updated:
# 		with open('temp00.ini', 'a') as f1:
# 			if action == 'user':
# 				new_values = Get_From_User_Line(setting_option, update=False)
# 			elif action == 'even':
# 				new_values = Only_Even_Channels('Dia_channel_screen_channels = {1,127}')
# 			elif action == 'odd':
# 				new_values = Only_Odd_Channels('Dia_channel_screen_channels = {0,126}')
# 			f1.write(setting_option + ' = ' + new_values + '\n')
# 	os.remove(s_file)
# 	shutil.move('temp00.ini', s_file)
# 	if setting_option.startswith('si_avg_fidcut'):
# 		return new_values
#
# def Get_From_User_Line(option, old_value='', update=True):
# 	if update:
# 		print 'The current value(s) for:', option, 'is(are):'
# 		print old_value
# 		temp = raw_input('Enter the new value(s) for ' + option + ' in the same format as above (Press enter to leave as it is): ')
# 		if temp == '' or temp == '\n':
# 			temp = old_value
# 	else:
# 		temp = raw_input('Enter the value(s) for ' + option + ' in the correct format: ')
# 	return temp.replace(' ', '')
#
# def Only_Even_Channels(old_value):
# 	channel_str = old_value[old_value.find('{') + 1:old_value.find('}')].split(',')
# 	channel_str = ['0', '126'] if len(channel_str) == 1 and channel_str[0] == '' else channel_str
# 	return Modify_String_Even_or_Odd(channel_str, 'even')
#
# def Only_Odd_Channels(old_value):
# 	channel_str = old_value[old_value.find('{') + 1:old_value.find('}')].split(',')
# 	channel_str = ['1', '127'] if len(channel_str) == 1 and channel_str[0] == '' else channel_str
# 	return Modify_String_Even_or_Odd(channel_str, 'odd')
#
# def Modify_String_Even_or_Odd(channel_old, type='even'):
# 	if type == 'even':
# 		modulo = 1
# 		if '1' != channel_old[0] and not '1-' in channel_old[0]:
# 			if '0' != channel_old[0] and not '0-' in channel_old[0]:
# 				channel_old.insert(0, '1')
# 		if '127' != channel_old[-1] and not '-127' in channel_old[-1]:
# 			channel_old.append('127')
# 		channel_old.append('129')
# 	else:
# 		modulo = 0
# 		if '0' != channel_old[0] and not '0-' in channel_old[0]:
# 			channel_old.insert(0, '0')
# 		if '126' != channel_old[-1] and not '-126' in channel_old[-1]:
# 			if '127' != channel_old[-1] and not '-127' in channel_old[-1]:
# 				channel_old.append('126')
# 		channel_old.append('128')
# 	channel_new = ''
# 	for i in xrange(1, len(channel_old)):
# 		if IsInt(channel_old[i - 1]):
# 			prev = int(channel_old[i - 1])
# 			if IsInt(channel_old[i]):
# 				th = int(channel_old[i])
# 				if th - prev < 2:
# 					channel_new += str(prev) + ','
# 				elif prev % 2 == modulo:
# 					for ch in xrange(prev, th, 2):
# 						channel_new += str(ch) + ','
# 				else:
# 					channel_new += str(prev) + ','
# 					for ch in xrange(prev + 1, th, 2):
# 						channel_new += str(ch) + ','
# 			else:
# 				temp = channel_old[i].split('-')
# 				if IsInt(temp[0]):
# 					th = int(temp[0])
# 					if th - prev < 2:
# 						channel_new += str(prev) + ','
# 					elif prev % 2 == modulo:
# 						for ch in xrange(prev, th, 2):
# 							channel_new += str(ch) + ','
# 					else:
# 						channel_new += str(prev) + ','
# 						for ch in xrange(prev + 1, th, 2):
# 							channel_new += str(ch) + ','
# 		else:
# 			channel_new += channel_old[i - 1] + ','
# 			prev = int(channel_old[i - 1].split('-')[-1])
# 			if not IsInt(channel_old[i]):
# 				th = int(channel_old[i].split('-')[0])
# 			else:
# 				th = int(channel_old[i])
# 			if th - prev >= 2:
# 				if prev % 2 == modulo:
# 					for ch in xrange(prev + 2, th, 2):
# 						channel_new += str(ch) + ','
# 				else:
# 					for ch in xrange(prev + 1, th, 2):
# 						channel_new += str(ch) + ','
#
# 	channel_new = channel_new[:-1]
# 	return '{' + channel_new + '}'
#
# def RecreateLink(source, dest, name, doSymlink=True, doCopy=True):
# 	successful = True
# 	if os.path.isdir(source):
# 		if doSymlink:
# 			successful = RecreateSoftLink(source, dest, name, type='dir', doCopy=doCopy)
# 		elif doCopy:
# 			if os.path.islink(dest + '/' + name):
# 				os.unlink(dest + '/' + name)
# 			elif os.path.isdir(dest + '/' + name):
# 				shutil.rmtree(dest + '/' + name, True)
# 			shutil.copytree(source, dest + '/' + name)
# 			successful = True
# 		else:
# 			successful = False
# 			print 'Could not link', source, 'to', dest, 'and copy was disabled'
#
# 	elif os.path.isfile(source):
# 		if os.path.isfile(dest + '/' + name):
# 			os.unlink(dest + '/' + name)
# 		if doSymlink:
# 			successful = RecreateSoftLink(source, dest, name, type='file', doCopy=doCopy)
# 		if not doSymlink or not successful:
# 			if os.path.isfile(dest + '/' + name):
# 				os.unlink(dest + '/' + name)
# 			try:
# 				os.link(source, dest + '/' + name)
# 				successful = True
# 			except OSError:
# 				if doCopy:
# 					shutil.copy2(source, dest + '/' + name)
# 					successful = True
# 				else:
# 					successful = False
# 					print 'Could not link', source, 'to', dest, 'and copy was disabled'
# 	else:
# 		successful = False
# 	return successful
#
# def RecreateSoftLink(source, dest, name, type='dir', doCopy=True):
# 	if type == 'dir':
# 		if os.path.isdir(source):
# 			if os.path.islink(dest + '/' + name):
# 				os.unlink(dest + '/' + name)
# 			elif os.path.isdir(dest + '/' + name):
# 				shutil.rmtree(dest + '/' + name, True)
# 			try:
# 				os.symlink(source, dest + '/' + name)
# 				return True
# 			except OSError:
# 				if doCopy:
# 					shutil.copytree(source, dest + '/' + name)
# 					return True
# 	else:
# 		if os.path.isfile(source):
# 			if os.path.isfile(dest + '/' + name):
# 				os.unlink(dest + '/' + name)
# 			try:
# 				os.symlink(source, dest + '/' + name)
# 				return True
# 			except OSError:
# 				if doCopy:
# 					shutil.copy2(source, dest + '/' + name)
# 					return True
# 	return False
# 	# print 'Subprocess terminated successfully'
#
# def ReturnRGB(val, min, max):
# 		hue = (val - min) * 1.0 / (max - min) if max != min else 0
# 		huep = hue * 6.0
# 		i = int(huep)
# 		v1 = 0.
# 		v2 = 1. - huep + i
# 		v3 = huep - i
# 		if i == 0:
# 			r, g, b = 1., v3, v1
# 		elif i == 1:
# 			r, g, b = v2, 1., v1
# 		elif i == 2:
# 			r, g, b = v1, 1., v3
# 		elif i == 3:
# 			r, g, b = v1, v2, 1.
# 		elif i == 4:
# 			r, g, b = v3, v1, 1.
# 		else:
# 			r, g, b = 1., v1, v2
# 		return r, g, b
#
# def Correct_Path(path):
# 	abs_path = ''
# 	if path[0] == '~':
# 		abs_path += os.path.expanduser('~')
# 		abs_path += path[1:]
# 	elif os.path.isabs(path):
# 		abs_path += path
# 	else:
# 		abs_path += os.path.abspath(path)
# 	return abs_path
#
# def DeleteDirectoryContents(dir):
# 	if os.path.isdir(dir):
# 		print 'Deleting old analysis...', ; sys.stdout.flush()
# 		for file in os.listdir(dir):
# 			file_path = os.path.join(dir, file)
# 			try:
# 				if os.path.islink(file_path) or os.path.isfile(file_path):
# 					os.unlink(file_path)
# 				elif os.path.isdir(file_path):
# 					shutil.rmtree(file_path)
# 			except Exception as e:
# 				print(e)
# 		print 'Done'
