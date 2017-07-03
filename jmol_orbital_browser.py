#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 14:10:18 2017

@author: mzwier
"""

import subprocess, argparse, os, sys, socket, re, collections, time
import tkinter as tk

JVXLInfo = collections.namedtuple('JVXLInfo', ['filename', 'jmol_label', 'list_entry'])

class JmolGone(RuntimeError):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class MatchCapturer:
    def __init__(self, regexp):
        self.regexp = regexp
        self._match = None
        self._searched = False
        
    def search(self,*args,**kwargs):
        self._searched = True
        self._match = self.regexp.search(*args,**kwargs)
        return self._match
        
    def match(self,*args,**kwargs):
        self._searched = True
        self._match = self.regexp.match(*args,**kwargs)
        return self._match        

def get_free_port():
    s = socket.socket()
    s.bind(('', 0))
    port = s.getsockname()[1]
    s.close()
    return port

class ScrolledList(tk.Frame):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.pack(expand=tk.YES, fill=tk.BOTH)
        
        self.scroll_bar = tk.Scrollbar(self)
        self.list_box = tk.Listbox(self, relief=tk.SUNKEN)
        self.scroll_bar.config(command=self.list_box.yview)
        self.list_box.config(yscrollcommand=self.scroll_bar.set)
        self.scroll_bar.pack(side=tk.RIGHT, fill=tk.Y)
        self.list_box.pack(side=tk.LEFT, expand=tk.YES, fill=tk.BOTH)

class JmolOrbitalControl:
    default_jmol_command = 'jmol'
    default_jmol_encoding = 'UTF8'
    default_jmol_startup_timeout = 10
    default_jmol_shutdown_timeout = 2
    re_jvxl_file = re.compile('^(.*?)\.mo(\d+)([^\.]*)\.jvxl')    
    
    def __init__(self):
        self.window = None
        self.scrolled_list = None
        self.jmol_command = self.default_jmol_command
        self.jmol_encoding = self.default_jmol_encoding
        self.jmol_popen = None
        self.port = get_free_port()
        self.jmol_socket = None
        self.jmol_startup_timeout = self.default_jmol_startup_timeout
        
        # coordinate file
        self.coord_file = None        
        
        # A correctly sorted list of JVXLInfo tuples
        self.jvxl_files = None
        
        self.jvxl_loaded = set()    
        self.jvxl_displayed = set()
        
    def find_xyzfile(self):
        xyzfiles = [filename for filename in os.listdir() if filename.endswith('.xyz')]
        if not xyzfiles:
            sys.stderr.write('Cannot find an xyz file in the current directory.\n'
                            +'Specify one with the -c option.\n')
            sys.exit(1)
        elif len(xyzfiles) > 1:
            sys.stderr.write('Found multiple xyz files in the current directory.\n'
                            +'Specify one with the -c option.\n')
            sys.exit(1)        
        else:
            return xyzfiles[0]
            
    def get_jvxlfiles(self, filenames):
        jvxl_dict = dict()
        if not filenames:
            filenames = [filename for filename in os.listdir() if filename.endswith('.jvxl')]
    
        # form a dictionary with keys (basename, MO #, a vs b etc) pointing to filenames
        stems = set()
        suffixes = set()
        for filename in filenames:
            m = self.re_jvxl_file.search(filename)
            if not m:
                raise ValueError('Cannot get MO number of JVXL file {}'.format(filename))
            stems.add(m.group(1))
            suffixes.add(m.group(3))
            jvxl_dict[m.group(1), int(m.group(2)), m.group(3)] = filename
            
        # Now we have a mapping of (stem, MO number, suffix) -> filename
        # Turn it into a mapping of (stem, MO number, suffix) -> (list label, jmol label, filename)
        
        # Jmol labels are simply mo<NUMBER> if all stems and suffixes are the same
        # Otherwise if stems are the same but suffixes are different, mo<NUMBER><SUFFIX>
        # Otherwise, orbset<MOLECULE NUMBER>mo<NUMBER><SUFFIX>
    
        jvxl_files = []
        if len(stems) == 1 and len(suffixes) == 1:
            for key in sorted(jvxl_dict.keys()):
                (stem, nmo, suffix) = key
                jmol_label = 'mo{:d}'.format(nmo)
                list_entry = 'MO {:d}'.format(nmo)
                jvxl_files.append(JVXLInfo(filename = jvxl_dict[key],
                                           jmol_label = jmol_label,
                                           list_entry = list_entry))
        elif len(stems) == 1: # and len(suffixes) > 1
            for key in sorted(jvxl_dict.keys()):
                (stem, nmo, suffix) = key
                jmol_label = 'mo{:d}{}'.format(nmo, suffix)
                list_entry = 'MO {:d}{}'.format(nmo, suffix)
                jvxl_files.append(JVXLInfo(filename = jvxl_dict[key],
                                           jmol_label = jmol_label,
                                           list_entry = list_entry))
        else: # len(stems) > 1
            stem_map = {stem: istem for istem, stem in enumerate(stems)}
            for key in sorted(jvxl_dict.keys()):
                (stem, nmo, suffix) = key
                jmol_label = 'os{}mo{:d}{}'.format(stem_map[stem],nmo, suffix)
                list_entry = '{} MO {:d}{}'.format(stem, nmo, suffix)
                jvxl_files.append(JVXLInfo(filename = jvxl_dict[key],
                                           jmol_label = jmol_label,
                                           list_entry = list_entry))
                                           
        self.jvxl_files = jvxl_files
        
        
    def run_jmol(self):
        self.jmol_popen = subprocess.Popen([self.jmol_command, '-I'],
                                           stdin = subprocess.PIPE)
        self.jmol_send_pipe('SYNC -{:d} \n'.format(self.port))
        assert self.coord_file is not None
        self.jmol_send_socket('load {}'.format(self.coord_file))
        
    def jmol_send_pipe(self, script, flush=True):
        if not script.endswith(' \n'):
            script += ' \n'
        #print('sending {!r}'.format(script))
        self.jmol_popen.stdin.write(script.encode(self.jmol_encoding))
        if flush:
            self.jmol_popen.stdin.flush()
            
    def jmol_send_socket(self, command):
        if self.jmol_popen.poll() is not None:
            raise JmolGone('jmol has quit')
        if self.jmol_socket is None:
            self.jmol_socket = socket.socket()
            for n in range(self.jmol_startup_timeout):
                try:
                    self.jmol_socket.connect(('localhost', self.port))
                    break
                except ConnectionRefusedError as e:
                    time.sleep(1)
                    continue
            else:
                raise e
        
        message =  '{"type" : "login", "source" : "Jmol"}\n'
        message += '{"type" : "command", "command" : "%s" }\n' % command
        message += '{"type" : "script", "event" : "done"}\n'
        
        self.jmol_socket.send(message.encode(self.jmol_encoding))
        
    def build_window(self):
        self.window = tk.Tk()
        self.window.title('Jmol orbital control')
        self.scrolled_list = ScrolledList(self.window)
        for jvxl_file in self.jvxl_files:
            self.scrolled_list.list_box.insert(tk.END, jvxl_file.list_entry)
        self.scrolled_list.list_box.config(selectmode=tk.MULTIPLE)
        self.scrolled_list.list_box.bind('<<ListboxSelect>>',
            (lambda event: 
                self.selection_change(set(int(index) 
                for index in self.scrolled_list.list_box.curselection()))))
        tk.Button(self.window, text='Exit', command=self.window.quit).pack(side=tk.BOTTOM)        
        tk.Button(self.window, text='Deselect all', command=self.deselect_all)\
          .pack(side=tk.BOTTOM)
        
        
    
    def deselect_all(self):
        self.scrolled_list.list_box.select_clear(0,tk.END)
        self.hide_mos(set(self.jvxl_displayed))
    
    def mainloop(self):
        self.window.mainloop()
        self.shutdown()

    def shutdown(self):
        try:
            self.jmol_send_socket('exitJmol')        
        except JmolGone:
            pass
        self.jmol_popen.wait()
        self.jmol_socket.close()
        
    def load_mos(self, indices):
        for index in indices:
            jvxl_file = self.jvxl_files[index]
            cmd = 'isosurface ID {jmol_label:s} {filename:s}'\
                  .format(jmol_label=jvxl_file.jmol_label, filename=jvxl_file.filename)
            self.jmol_send_socket(cmd)
            self.jvxl_loaded.add(index)
            self.jvxl_displayed.add(index)
            
    def show_mos(self, indices):
        for index in indices:
            cmd = 'isosurface {} on'.format(self.jvxl_files[index].jmol_label)
            self.jmol_send_socket(cmd)
            self.jvxl_displayed.add(index)
            
    def hide_mos(self, indices):
        for index in indices:
            cmd = 'isosurface {} off'.format(self.jvxl_files[index].jmol_label)
            self.jmol_send_socket(cmd)
            self.jvxl_displayed.remove(index)
        
    def selection_change(self, indices):        
        selected = set(indices)
        
        # Find displayed but not selected MOs, then hide        
        to_hide = self.jvxl_displayed - selected
        self.hide_mos(to_hide)
        
        # Find selected but not loaded MOs, then display
        unloaded = selected - self.jvxl_loaded
        self.load_mos(unloaded)
        
        # Find selected but not displayed MOs, then display
        undisplayed = selected - self.jvxl_displayed
        self.show_mos(undisplayed)

    
def vprint(*args, **kwargs):
    global verbose
    if verbose:
        print(*args, **kwargs)
        
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--coordinates', 
                    help='Coordinate file to load. Defaults to any file ending in .xyz. '
                        +'If multiple .xyz files are found, the program will not run.')
parser.add_argument('jvxlfile', nargs='*',
                    help='jvxl MO files to load. Defaults to all files ending in .jvxl. '
                        +'An attempt is made to sort these files numerically.')
parser.add_argument('-v', '--verbose', action='store_true',
                    help="Be descriptive about what's going on")
                        
args = parser.parse_args()
verbose = args.verbose


joc = JmolOrbitalControl()

joc.coord_file = args.coordinates or joc.find_xyzfile()
vprint('Using {} as structure file'.format(joc.coord_file))

joc.get_jvxlfiles(args.jvxlfile)
if verbose:
    vprint('Found the following orbitals:')
    for jvxl_file in joc.jvxl_files:
        vprint('  {} : {}'.format(jvxl_file.list_entry, jvxl_file.filename))



joc.run_jmol()
joc.build_window()
joc.mainloop()
