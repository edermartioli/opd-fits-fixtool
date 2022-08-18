# Description: Script to fix raw data
# Author: Eder Martioli
# Laboratorio Nacional de Astrofisica, Brazil
# 19 Jul 2022
#
# Example:
# python /Volumes/Samsung_T5/sparc4-pipeline/fixtool/checkkeys.py --input=*.fits

__version__ = "1.0"

__copyright__ = """
    Copyright (c) ...  All rights reserved.
    """

from optparse import OptionParser
import os, sys
import glob
import astropy.io.fits as fits
import numpy as np

from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

#from astroquery.simbad import Simbad

import fix_params


parser = OptionParser()
parser.add_option("-r", "--datadir", dest="datadir", help="data directory",type='string',default="./")
parser.add_option("-i", "--input", dest="input", help="input data pattern",type='string',default="*.fits")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="verbose",default=False)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print("Error: check usage with fixrawdata.py -h ")
    sys.exit(1)

if options.verbose:
    print('Data directory: ', options.datadir)
    print('Data input: ', options.input)

currdir = os.getcwd()
os.chdir(options.datadir)

inputdata = sorted(glob.glob(options.input))

p = fix_params.load_obsfix_parameters()

logstr = 'INDEX\tIMAGE'
for key in p['IMPORTANT_KEYS'] :
    logstr += "\t{}".format(key)
logstr += '\n'

for i in range(len(inputdata)) :
    hdr = fits.getheader(inputdata[i])
    logstr += "{}\t{}".format(i,os.path.basename(inputdata[i]))
    for key in p['IMPORTANT_KEYS'] :
        if type(hdr[key]) == str :
            logstr += "\t'{}'".format(hdr[key])
        elif type(hdr[key]) == float and key != 'BJD':
            logstr += "\t{:.2f}".format(hdr[key])
        elif type(hdr[key]) == float and key == 'BJD':
            logstr += "\t{:.6f}".format(hdr[key])
        else :
            logstr += "\t{}".format(hdr[key])

    logstr += '\n'
    
print(logstr)

os.chdir(currdir)

