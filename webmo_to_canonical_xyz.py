#!/usr/bin/env python3

import argparse, sys

import numpy as np

from mendeleev import element

def read_webmo_xyz(xyzfile):
    points = []
    energies = []
    titles = []
    all_symbols = []

    
    for line in xyzfile:
        all_coords = []
        symbols = []
        title = line
        titles.append(title)
        
        while True:
            coord_line = xyzfile.readline()
            if not coord_line.strip():
                break
                
            (atom, x, y, z) = coord_line.split()

            if not all_symbols:
                try:
                    atom = int(atom)
                except ValueError:
                    symbol = atom
                else:
                    symbol = element(atom).symbol
                symbols.append(symbol)
            
            coords = np.array((x,y,z), dtype=np.float64)                
            all_coords.append(coords)

        if symbols:
            all_symbols = symbols
        points.append(np.array(all_coords))
    
    return (np.array(points), all_symbols, titles)

def write_xyz(outfile, coords, symbols, titles):
    n_atoms = coords.shape[1]
    
    for coord_set, title in zip(coords, titles):
        outfile.write('{:d}\n'.format(n_atoms))
        outfile.write('{:s}\n'.format(title.strip()))

        for iatom in range(len(coord_set)):
            outfile.write('{:>3s} {:14.6f} {:14.6f} {:14.6f}\n'.format(symbols[iatom], *list(coord_set[iatom,:])))
    
    
parser = argparse.ArgumentParser()
parser.add_argument('webmo_xyz', help='WebMO XYZ trajectory')
parser.add_argument('-o', '--output', help='Output filename (default stdout)')
args = parser.parse_args()

xyzfile = open(args.webmo_xyz, 'rt')
points, symbols, titles = read_webmo_xyz(xyzfile)


if args.output:
    outfile = open(args.output, 'wt')
else:
    outfile = sys.stdout

write_xyz(outfile, points, symbols, titles)
