#!/usr/bin/env python3

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('xyz1')
parser.add_argument('xyz2')
parser.add_argument('--noalign', action='store_true',
                    help='Do not align structures prior to calculating RMSD')
parser.add_argument('--nocenter', action='store_true',
                    help='Do not center structures prior to calculating RMSD. Implies --noalign.')
args = parser.parse_args()

xyz1 = np.loadtxt(args.xyz1, skiprows=2, usecols=[1,2,3])
xyz2 = np.loadtxt(args.xyz2, skiprows=2, usecols=[1,2,3])

if args.nocenter:
    args.noalign = True
else:
    xyz1 -= xyz1.mean(axis=0)
    xyz2 -= xyz2.mean(axis=0)
    
if args.noalign:
    rmsd = ((xyz2 - xyz1)**2).sum(axis=1).mean(axis=0) ** 0.5
else:
    import pyximport
    pyximport.install(inplace=True, language_level=3,
                      setup_args={"include_dirs":np.get_include()})
    import qcprot 
    
    rot_matrix = np.empty((9,), dtype=np.float64)
    
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, len(xyz1), rot_matrix, None)

print('{}'.format(rmsd))

