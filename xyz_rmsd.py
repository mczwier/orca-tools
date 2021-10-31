#!/usr/bin/env python3

import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('xyz1')
parser.add_argument('xyz2')
parser.add_argument('--noalign', action='store_true',
                    help='Do not align structures prior to calculating RMSD')
parser.add_argument('--nocenter', action='store_true',
                    help='Do not center structures prior to calculating RMSD')
args = parser.parse_args()

xyz1 = np.loadtxt(args.xyz1, skiprows=2, usecols=[1,2,3])
xyz2 = np.loadtxt(args.xyz2, skiprows=2, usecols=[1,2,3])

if not args.nocenter:
    xyz1 -= xyz1.mean(axis=0)
    xyz2 -= xyz2.mean(axis=0)

if args.noalign:
    rmsd = ((xyz2 - xyz1)**2).mean() ** 0.5
else:
    from MDAnalysis.lib import qcprot
    rot_matrix = np.empty((9,), dtype=np.float64)
    rmsd = qcprot.CalcRMSDRotationalMatrix(xyz1, xyz2, len(xyz1), rot_matrix, None)

print('{}'.format(rmsd))
