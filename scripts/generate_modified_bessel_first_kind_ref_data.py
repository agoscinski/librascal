
import sys,os
import json
import ase
import argparse
from mpmath import mp,besseli,pi,sqrt,exp

import numpy as np
import ubjson

def get_leggauss(order, a, b):
    x,w = leggauss(order)
    x = (b-a)*0.5 * x + 0.5*(a+b)
    w = (b-a)*0.5 * w
    return x,w

mp.dps = 20; mp.prec = 100;

def sbesseli(n,z):
    """e^{-x}*i_n(x)"""
    return sqrt(pi/(2*z))*besseli(n+0.5,z)*exp(-z)

#dump radial and power spectra for methane
def dump_reference_json():
    path = '../'
    sys.path.insert(0, os.path.join(path, 'build/'))
    sys.path.insert(0, os.path.join(path, 'tests/'))
    data = []
    orders = list(range(1,20))
    for x in np.linspace(10,6000,500):
        vals = []
        for order in orders:
            val = sbesseli(order, x)
            vals.append(float(val))
        data.append(dict(x=x,max_order=orders[-1],vals=vals))
    print(len(data))
    with open(path+"tests/reference_data/modified_bessel_first_kind_reference.ubjson",'wb') as f:
        ubjson.dump(data,f)

##########################################################################################
##########################################################################################

def main(json_dump):
    if json_dump == True:
        dump_reference_json()

##########################################################################################
##########################################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-json_dump', action='store_true', help='Switch for dumping json')

    args = parser.parse_args()
    main(args.json_dump)
