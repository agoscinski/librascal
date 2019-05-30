#!/usr/bin/env/python3

from matplotlib import pylab as plt

import os, sys
from ase.io import read
sys.path.insert(0,"../build/")

import sys
import time
import rascal
import json

import ase
from ase.io import read, write
from ase.build import make_supercell
from ase.visualize import view
import numpy as np
import sys

import json

from rascal.representations import SOAP

structure = read('data/problematic_structure_numbers.json', format='json')

view(structure)
write('prob_structure.png', structure)
write('data/problematic_structure.xyz', structure, format='xyz')

frames = read('./data/problematic_structure.xyz', ':')

hypers = dict(soap_type="PowerSpectrum",
              interaction_cutoff=3.5,
              max_radial=6,
              max_angular=4,
              gaussian_sigma_constant=0.4,
              gaussian_sigma_type="Constant",
              cutoff_smooth_width=0.5,
              )
soap = SOAP(**hypers)

representation = soap.transform(frames)
