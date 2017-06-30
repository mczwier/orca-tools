#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:51:03 2017

@author: mzwier
"""

import argparse, subprocess, re, os, sys, textwrap, select

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
        

class OrcaPlotException(RuntimeError):
    def __init__(self, *args, stdout, stderr):
        super(OrcaPlotException).__init__(*args)
        self.stdout = stdout
        self.stderr = stderr
        
class OrcaPlotInterface:
    default_orca_plot_command = 'orca_plot'
    default_orca_encoding = 'UTF8'
    
    '''Number of available orbitals : 56'''
    re_norbs = re.compile(r'Number of available orbitals\s*:\s*(\d+)')
    
    '''Output file: minimal.mo0a.cube'''
    re_outputfile = re.compile(r'Output file:\s+(.+)$')
    
    def __init__(self, gbwfile, n_points=None, logfile=None):
        self.gbwfile = gbwfile
        self.orca_encoding = self.default_orca_encoding
        self.orca_plot_command = self.default_orca_plot_command
        self.n_orbitals = None
        self.n_points = n_points
        if logfile:
            self.logfile = open(logfile, 'wb')
        else:
            self.logfile = None
        
        basename, ext = os.path.splitext(os.path.basename(self.gbwfile))
        self.xyzfile = basename + '.xyz'
        
    def run_orca_plot(self, input_text):
        if not input_text.endswith('11\n'):
            input_text += '\n11\n'
        input_bytes = input_text.encode(self.orca_encoding)
        
        if self.logfile:
            self.logfile.write('<<<\n'.encode(self.orca_encoding))
            self.logfile.write(input_bytes)
        
        args = [self.orca_plot_command, self.gbwfile, '-i']
        pop = subprocess.Popen(args,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = pop.communicate(input_bytes)
        if self.logfile:
            self.logfile.write('>>>\n'.encode(self.orca_encoding))
            self.logfile.write(stdout)
            self.logfile.write('EEE\n'.encode(self.orca_encoding))
            self.logfile.write(stderr)
        stdout_text = stdout.decode(self.orca_encoding)
        stderr_text = stderr.decode(self.orca_encoding)
        rc = pop.wait()
        if rc:
            print(stdout_text)
            print(stderr_text)
            raise subprocess.CalledProcessError(rc, ' '.join(args), 
                                                output=stdout, stderr=stderr)
        return stdout_text, stderr_text
        
    def get_n_orbitals(self):
        stdout, stderr = self.run_orca_plot('')
        for line in stdout.splitlines():
            m = self.re_norbs.search(line)
            if m:
                return int(m.group(1))
        else:
            raise OrcaPlotException('Cannot determine number of orbitals', 
                                    stdout=stdout, stderr=stderr)
            
    def save_orbital_cube(self, orbital):
        script='2\n{orbital:d}\n5\n7\n'.format(orbital=orbital)
        if self.n_points:
            script+='4\n{}\n'.format(self.n_points)
        script += '10\n'
        stdout, stderr = self.run_orca_plot(script)
        for line in stdout.splitlines():
            m = self.re_outputfile.search(line)
            if m:
                return m.group(1)
        else:
            raise OrcaPlotException('Cannot determine output file name',
                                    stdout=stdout, stderr=stderr)

class JmolCmdLineInterface:
    default_jmol_command = 'jmol'
    default_jmol_encoding = 'UTF8'
    
    def __init__(self):
        self.jmol_command = self.default_jmol_command
        self.jmol_encoding = self.default_jmol_encoding
        
    def run_jmol_script(self, script):
        script = script.encode(self.jmol_encoding)
        
        args = [self.jmol_command, '-i', '-o', '-n', '-s', '-']
        pop = subprocess.Popen(args,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = pop.communicate(script)
        stdout_text = stdout.decode(self.jmol_encoding)
        stderr_text = stderr.decode(self.jmol_encoding)
        rc = pop.wait()
        if rc:
            print(stdout_text)
            print(stderr_text)
            raise subprocess.CalledProcessError(rc, ' '.join(args), 
                                                output=stdout, stderr=stderr)
        return stdout_text, stderr_text
        
    def cube_to_jvxl(self, cubefile, jvxlfile = None):
        jvxlfile = jvxlfile or os.path.splitext(cubefile)[0] + '.jvxl'
        
        script = '''\n
            zap            
            isosurface SIGN {cubefile:s}
            write isosurface {jvxlfile:s}
            isosurface DELETE
            '''.format(cubefile=cubefile, jvxlfile=jvxlfile)
        self.run_jmol_script(script)
        return jvxlfile
         
def parse_orbital_range(n_mos, mos):
    '''Parse a list of orbital specifications. Each orbital specification
    is either an integer or a range of integers specified as "M-N" indicating
    that all orbitals M-N (inclusive) should be generated. N may be omitted to specify
    all MOs >= M.'''
    
    if not mos:
        return list(range(n_mos))
    
    nmos = set()
    
    
    for item in mos:
        item = item.strip()
        if item[-1] == '-':
            M = int(item[:-1])
            nmos.update(range(M,n_mos))
        elif '-' in item:
            M,N = (int(bound.strip()) for bound in item.split('-'))
            nmos.update(range(M,N+1))
        else:
            nmos.add(int(item))
    
    return list(sorted(nmos))
            
            
def vprint(*args, **kwargs):
    global verbose
    if verbose:
        print(*args,**kwargs)
    

parser = argparse.ArgumentParser()
parser.add_argument('gbwfile', help='Wavefunction file')
parser.add_argument('mos', nargs='*', 
                    help='Molecular orbitals to plot. Each MOS entry can be an integer, '
                        +'a range of integers M-N (inclusive), or the open range M- to '
                        +'plot all orbitals >= M. (Default: all)')
parser.add_argument('-N', '--npts', type=int, help='Number of grid points.')
parser.add_argument('-s', '--jmol-script', help='Filename of Jmol script to write.')
parser.add_argument('-S', '--no-jmol-script', action='store_true', 
                    help='Do not write Jmol script')
parser.add_argument('-L', '--orcalog', help='Write orca_plot output to ORCALOG.')                    
parser.add_argument('-v', '--verbose', action='store_true', 
                    help='Display progress and extra information.')

args = parser.parse_args()
gbwfilename = args.gbwfile
jmolfilename = args.jmol_script or os.path.splitext(os.path.basename(gbwfilename))[0] + '.jmol'
verbose = args.verbose
write_jmol = not args.no_jmol_script

orcaplot = OrcaPlotInterface(gbwfilename, n_points = args.npts, logfile=args.orcalog)
jmol = JmolCmdLineInterface()
n_mos = orcaplot.get_n_orbitals()

vprint('Orca file {} contains {:d} orbitals'.format(gbwfilename, n_mos))
vprint('Structure file is {}'.format(orcaplot.xyzfile))

mos = parse_orbital_range(n_mos, args.mos)
for nmo in mos:
    if nmo < 0 or nmo >= n_mos:
        print('Invalid MO specified: {}'.format(nmo), file=sys.stderr)
        sys.exit(1)

if verbose:
    vprint('Plotting the following orbitals:')
    if sys.stdout.isatty():
        print('\n'.join(textwrap.wrap(str(mos))))
    else:
        print(mos)


if write_jmol:
    jmolfile = open(jmolfilename, 'wt')
    jmolfile.write('load {}\n'.format(orcaplot.xyzfile))
    
for nmo in mos:
    vprint('Processing MO {}'.format(nmo))
    cubefile = orcaplot.save_orbital_cube(nmo)
    vprint('  Saved {}'.format(cubefile))
    jvxlfile = jmol.cube_to_jvxl(cubefile)
    vprint('  Converted {} to {}'.format(cubefile,jvxlfile))
    os.unlink(cubefile)
    vprint('  Deleted {}'.format(cubefile))
    
    if write_jmol:
        mo_id = 'mo{:d}'.format(nmo)
        jmolfile.write('isosurface ID {} {}\n'.format(mo_id, jvxlfile))
        jmolfile.write('isosurface {} off\n'.format(mo_id))

if write_jmol:
    vprint('Wrote Jmol script to {}'.format(jmolfilename))
    jmolfile.close()
    