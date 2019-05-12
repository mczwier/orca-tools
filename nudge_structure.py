#!/usr/bin/env python3

import sys
import numpy as np
import argparse

import scipy.constants
ANGSTROMS_PER_BOHR = scipy.constants.value('Bohr radius') * 1e10

class OrcaHessian:

    def __init__(self):
        pass

    def read(self, hessfile):
        for line in hessfile:
            if line.startswith('$vibrational_frequencies'):
                self._parse_vibrational_frequencies(hessfile)
            elif line.startswith('$normal_modes'):
                self._parse_normal_modes(hessfile)
            elif line.startswith('$atoms'):
                self._parse_atoms(hessfile)
            else:
                continue
        self._finalize()

    def _parse_vibrational_frequencies(self, hessfile):
        #nfreq = int(hessfile.readline())
        #freqs = [float(hessfile.readline().strip().split()[1]) for ifreq in range(nfreq)]
        #self.frequencies = np.array(freqs)
        self.frequencies = self.parse_block(hessfile)

    def _parse_normal_modes(self, hessfile):
        nmodearray = self.parse_block(hessfile)
        nmode = nmodearray.shape[1]
        natom = nmode // 3

        nmode_displacements = np.empty((nmodearray.shape[1], natom, 3), np.float64)
        for iatom in range(natom):
            for imode in range(nmode):
                nmode_displacements[imode, iatom, :] = nmodearray[iatom*3:(iatom+1)*3, imode]

        # These displacements are mass-weighted and normalized
        self.nmode_displacements = nmode_displacements

    def _parse_atoms(self, hessfile):
        natoms = int(next(hessfile).strip())

        self.atoms = [None]*natoms
        self.masses = np.empty((natoms,), np.float32)
        self.coords = np.empty((natoms, 3), np.float64)

        for iatom in range(natoms):
            fields = next(hessfile).strip().split()
            self.atoms[iatom] = fields[0]
            self.masses[iatom] = float(fields[1])
            self.coords[iatom] = [float(field) for field in fields[2:]]

        self.coords *= ANGSTROMS_PER_BOHR

    def _finalize(self):
        # Normal modes are mass-weighted and normalized; provide
        # unweighted displacements
        
        unweighted_displacements = np.copy(self.nmode_displacements)
        for imode in range(unweighted_displacements.shape[0]):
            unweighted_displacements[imode] = self.nmode_displacements[imode] * self.masses[:, None]**0.5
            norm = ((unweighted_displacements[imode]**2).sum())**0.5

            if norm > 1e-7:
                unweighted_displacements[imode] /= norm

        self.unweighted_nmode_displacements = unweighted_displacements

    @staticmethod
    def parse_block(infile):
        shape = [int(dim) for dim in next(infile).strip().split()]
        array = np.empty(shape, np.float64)

        if array.ndim == 1:
            for i in range(array.shape[0]):
                array[i] = float(next(infile).strip().split()[1])
        elif array.ndim == 2:
            col = 0
            while col < array.shape[1]:
                next(infile)  # Skip column labels
                for row in range(array.shape[0]):
                    # Parse row, discarding row label
                    fields = [float(field) for field in next(infile).strip().split()][1:]
                    array[row, col:col+len(fields)] = fields
                col += len(fields)
        else:
            raise NotImplementedError('only rank 1 and 2 data can be parsed')
        return array
    

parser = argparse.ArgumentParser(description='Nudge a structure along one or more modes. '
                                             'Produces an XYZ file in Angstroms. '
                                             'By default, nudges along all imaginary modes.')
parser.add_argument('hessian_file', help='ORCA hessian (.hess) file')
parser.add_argument('output_file', help='XYZ output file', nargs='?')
parser.add_argument('-l', '--list', action='store_true',
                    help='List frequencies and exit.')
parser.add_argument('-m', '--mode', '--modes', nargs='+',
                    help='Mode(s) to nudge along. Can be an integer, or N:d, where N is '
                         'the mode number and d is the displacement in Angstroms. '
                         'By default, nudges along all imaginary modes.')
parser.add_argument('-d', '--displacement', type=float, default=0.1,
                    help='Default distance for the nudge(s) in Angstroms (default: %(default)s)')
args = parser.parse_args()

hessfile = open(args.hessian_file, 'rt')
orcaparser = OrcaHessian()
orcaparser.read(hessfile)

if args.list:
    print('Frequencies:')
    for ifreq, freq in enumerate(orcaparser.frequencies):
        print('    {:6d}: {:10.1f} cm-1{}'.format(ifreq, freq, ' (*)' if freq < 0 else ''))
    sys.exit(0)


modes = []
displacements = []

if not args.mode:
    for imode in range(len(orcaparser.frequencies)):
        if orcaparser.frequencies[imode] < 0:
            modes.append(imode)
            displacements.append(args.displacement)

    if not modes:
        print('No imaginary frequencies and no modes specified. Nothing to do.')
        sys.exit(0)
else:
    for modestr in args.mode:
        imode, displacement = modestr.split(':')
        imode = int(imode)
        if displacement:
            displacement = float(displacement)
        else:
            displacement = args.displacement
        modes.append(imode)
        displacements.append(displacement)

if not args.output_file:
    print('Output file required unless only listing modes.', file=sys.stderr)
    sys.exit(1)

coords = np.copy(orcaparser.coords)
for imode, displacement in zip(modes, displacements):
    print('  Displacing mode {} by {} Angstroms'.format(imode, displacement))

    coords += displacement * orcaparser.unweighted_nmode_displacements[imode]

natoms = len(orcaparser.atoms)
with open(args.output_file, 'wt') as xyzfile:
    xyzfile.write('{:d}\n'.format(natoms))
    xyzfile.write('  Generated from {} by displacing {} modes\n'.format(args.hessian_file, len(modes)))
    for iatom in range(natoms):
        xyzfile.write('{:3s} {:f} {:f} {:f}\n'.format(orcaparser.atoms[iatom], *coords[iatom]))
                  
