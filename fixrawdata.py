# Description: Script to fix raw data
# Author: Eder Martioli
# Laboratorio Nacional de Astrofisica, Brazil
# 13 Jul 2022
#
# Examples:

"""
Using default configuration: python /Volumes/Samsung_T5/sparc4-pipeline/fixtool/fixrawdata.py --input=*.fits

ZE_22jul17_E -> python /Volumes/Samsung_T5/sparc4-pipeline/fixtool/fixrawdata.py --input=aumic_i_*.fits --telescope='ZE0.60m' --instrument='CAM4+Ixon4269' --acqsys="OPDAcquisition" --gain=3.8 --readnoise=8.2 --north='RIGHT' --east='DOWN' --timekey="DATE-OBS" --timetype="UT" --timezone=-3 -v

ZE_22jul17_W -> python /Volumes/Samsung_T5/sparc4-pipeline/fixtool/fixrawdata.py --input=aumic_i_*.fits --telescope='ZE0.60m' --instrument='CAM4+Ixon4269' --acqsys="OPDAcquisition" --gain=3.8 --readnoise=8.2 --north='LEFT' --east='UP' --timekey="DATE-OBS" --timetype="UT" --timezone=-3 -v


BC_22jul05:
Science: python /Volumes/Samsung_T5/sparc4-pipeline/fixtool/fixrawdata.py --input=*_c1_AUMIC_I.fits --telescope='BC0.60m' --instrument='CAM2+IxonUltra9915' --acqsys="GEI+ACS" --gain=3.3 --readnoise=6.57 --north='RIGHT' --east='DOWN' --timekey="DATE" --timetype="UT" --timezone=-3 -v
Bias: python /Volumes/Samsung_T5/sparc4-pipeline/fixtool/fixrawdata.py --input=*_c1_zero.fits --object="ZERO" --obstype="ZERO" --filter="None" --ra="00:00:00.00" --dec="00:00:00.00" --telescope='BC0.60m' --instrument='CAM2+IxonUltra9915' --acqsys="GEI+ACS" --gain=3.3 --readnoise=6.57 --north='RIGHT' --east='DOWN' --timekey="DATE" --timetype="UT" --timezone=-3 -v
Flat: python /Volumes/Samsung_T5/sparc4-pipeline/fixtool/fixrawdata.py --input=*_c1_flat_I.fits --object="FLAT" --obstype="FLAT" --filter="I" --ra="00:00:00.00" --dec="00:00:00.00" --telescope='BC0.60m' --instrument='CAM2+IxonUltra9915' --acqsys="GEI+ACS" --gain=3.3 --readnoise=6.57 --north='RIGHT' --east='DOWN' --timekey="DATE" --timetype="UT" --timezone=-3 -v
"""

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

from astropy.wcs import WCS


def convert_str2float_with_comma(input_string, default_value=-9999) :
    
    if type(input_string) == float :
        return input_string

    outputfloat = default_value
    
    if type(input_string) == str :
        if "," in input_string :
            try :
                s,f = input_string.split(",")
                f = "0."+f
                outputfloat = float(s) + float(f)
            except :
                print("WARNING: could not convert input string to float")
        else :
            try :
                outputfloat = float(input_string)
            except :
                print("WARNING: could not convert input string to float")

    return outputfloat


def set_siteobs_keys(p, hdr) :

    hdr.set("OBSTYPE",p['OBSTYPE'],"Observation type: OBJECT, ZERO, FLAT, DARK")
    hdr.set("OBJECT",p['OBJECT'],"Object name")
    hdr.set("OBSERVER",p['OBSERVER'],"Name of the observer")
    hdr.set("TELESCOP",p['TELESCOP'],"Telescope")
    hdr.set("FILTER",p['FILTER'],"Filter name")
    hdr.set("INSTRUME",p['INSTRUME'],"Instrument")
    hdr.set("ACQSYS",p['ACQSYS'],"Acquisition system")

    hdr.set("OBSALT",p['OBSALT'],"Observatory elevation above sea level (m)")
    hdr.set("OBSLONG",p['OBSLONG'],"Observatory East Longitude (DEG, 0 to 360)")
    hdr.set("OBSLAT",p['OBSLAT'],"Observatory North Latitude (DEG, -90 to 90)")
    hdr.set("EQUINOX",p['EQUINOX'],"Equinox")

    # calculate ambient temperature:
    ambienttemp = p['CRAZYVALUE']
    if "W-TEMP" in hdr.keys() :
        ambienttemp = convert_str2float_with_comma(hdr["W-TEMP"], default_value=p['CRAZYVALUE'])
    elif "TEMPEXT" in hdr.keys() :
        ambienttemp = convert_str2float_with_comma(hdr["TEMPEXT"], default_value=p['CRAZYVALUE'])

    hdr.set("TEMPEXT",ambienttemp,"External temperature (deg C)")

    # calculate ambient pressure:
    ambientpress = p['CRAZYVALUE']
    if "W-TEMP" in hdr.keys() :
        ambientpress = convert_str2float_with_comma(hdr["W-BAR"], default_value=p['CRAZYVALUE'])
    elif "PRESSURE" in hdr.keys() :
        ambientpress = convert_str2float_with_comma(hdr["PRESSURE"], default_value=p['CRAZYVALUE'])
    hdr.set("PRESSURE",ambientpress,"Pressure (mbar)")

    # calculate humidity:
    humidity = p['CRAZYVALUE']
    if "W-HUM" in hdr.keys() :
        humidity = convert_str2float_with_comma(hdr["W-HUM"], default_value=p['CRAZYVALUE'])
    elif "HUMIDITY" in hdr.keys() :
        humidity = convert_str2float_with_comma(hdr["HUMIDITY"], default_value=p['CRAZYVALUE'])
    hdr.set("HUMIDITY",humidity,"Humidity (per cent)")

    camfoc = p['CRAZYVALUE']
    try :
        if "FOCUSVAL" in hdr.keys() :
            camfoc = float(str(hdr["FOCUSVAL"]).replace("S",""))
        elif "CAMFOC" in hdr.keys() :
            camfoc = float(str(hdr["CAMFOC"]).replace("S",""))
    except :
        print("WARNING: could not get focus value")
    hdr.set("CAMFOC",camfoc,"Camera focus position (mm)")

    # calculate ccd temperature:
    ccdtemp = p['CRAZYVALUE']
    try :
        ccdtemp = convert_str2float_with_comma(hdr["TEMP"], default_value=p['CRAZYVALUE'])
    except :
        print("WARNING: could not get CCD temperature value")
        
    hdr.set("CCDTEMP",ccdtemp,"CCD temperature")

    hdr.set("GAIN",p['GAIN'],"Gain (e-/ADU)")
    #hdr.set("GAINERR",0.0,"Gain error (e-/ADU)")
    hdr.set("RDNOISE",p['RDNOISE'],"Read noise (e-)")
    #hdr.set("RNOISERR",0.0,"Read noise error (e-)")

    return hdr


def set_timecoords_keys(p, hdr) :

    # Set OPD geographic coordinates
    longitude = p['OBSLONG']
    latitude = p['OBSLAT']
    altitude = p['OBSALT']*u.m #hdr['OBSALT']
    opd_location = EarthLocation.from_geodetic(lat=latitude, lon=longitude, height=altitude)

    equinox="J{:.1f}".format(p['EQUINOX'])
    # set source observed
    source = SkyCoord(p['RA'], p['DEC'], unit=(u.hourangle, u.deg), frame='icrs', equinox=equinox)

    radeg, decdeg = source.ra.value, source.dec.value
    hdr.set("RA_DEG",radeg,"Requested Right Ascension (deg)")
    hdr.set("DEC_DEG",decdeg,"Requested Declination (deg)")
    hdr.set("TELRA",p['RA'],"TCS right ascension: HH:MM:SS.SSS")
    hdr.set("TELDEC",p['DEC'],"TCS declination: +-DD:MM:SS.SSS")
    hdr.set("RA",p['RA'],"Requested right Ascension: HH:MM:SS.SSS")
    hdr.set("DEC",p['DEC'],"Requested Declination: +- DD:MM:SS.SSS")

    """
        |a b| * |x| = |a x + b y|
        |c d|   |y|   |c x + d y|


    if p['NORTH'] == 'UP' or p['NORTH'] == 'DOWN':
        
        if p['NORTH'] == 'UP' and p['EAST'] == 'LEFT' :
            p['CD1_1'], p['CD1_2'] = np.cos(0.) , -np.sin(0.)
            p['CD2_1'], p['CD2_2'] = np.sin(0.) , np.cos(0.)
        elif p['NORTH'] == 'UP' and p['EAST'] == 'RIGHT' :
            p['CD1_1'], p['CD1_2'] = - np.cos(0.) , -np.sin(0.)
            p['CD2_1'], p['CD2_2'] = - np.sin(0.) , np.cos(0.)
        elif p['NORTH'] == 'DOWN' and p['EAST'] == 'LEFT' :
            p['CD1_1'], p['CD1_2'] = np.cos(0.) , +np.sin(0.)
            p['CD2_1'], p['CD2_2'] = np.sin(0.) , -np.cos(0.)
        elif p['NORTH'] == 'DOWN' and p['EAST'] == 'RIGHT' :
            p['CD1_1'], p['CD1_2'] = -np.cos(0.) , +np.sin(0.)
            p['CD2_1'], p['CD2_2'] = -np.sin(0.) , -np.cos(0.)

    elif p['NORTH'] == 'RIGHT' or p['NORTH'] == 'LEFT':
    
        if p['NORTH'] == 'RIGHT' and p['EAST'] == 'UP' :
            p['CD1_1'], p['CD1_2'] = np.cos(np.pi/2.) , -np.sin(np.pi/2.)
            p['CD2_1'], p['CD2_2'] = np.sin(np.pi/2.) , np.cos(np.pi/2.)
        elif p['NORTH'] == 'RIGHT' and p['EAST'] == 'DOWN' :
            p['CD1_1'], p['CD1_2'] = -np.cos(np.pi/2.) , -np.sin(np.pi/2.)
            p['CD2_1'], p['CD2_2'] = -np.sin(np.pi/2.) , np.cos(np.pi/2.)
        elif p['NORTH'] == 'LEFT' and p['EAST'] == 'UP' :
            p['CD1_1'], p['CD1_2'] = np.cos(-np.pi/2.) , +np.sin(-np.pi/2.)
            p['CD2_1'], p['CD2_2'] = np.sin(-np.pi/2.) , -np.cos(-np.pi/2.)
        elif p['NORTH'] == 'LEFT' and p['EAST'] == 'DOWN' :
            p['CD1_1'], p['CD1_2'] = np.cos(-np.pi/2.) , -np.sin(-np.pi/2.)
            p['CD2_1'], p['CD2_2'] = np.sin(-np.pi/2.) , np.cos(-np.pi/2.)
    """
    
    """
    hdr.set("WCSNAMEP",'PHYSICAL', "WCS name")
    
    hdr.set("CTYPE1",p['CTYPE1'],"Algorithm type for axis 1")
    hdr.set("CTYPE2",p['CTYPE2'],"Algorithm type for axis 2")

    hdr.set("CDELT1",p['CDELT1'],"RA pixel step (deg)")
    hdr.set("CDELT2",p['CDELT2'],"DEC pixel step (deg)")

    hdr.set("CRUNIT1",p['CRUNIT1'],"Unit of right ascension co-ordinates")
    hdr.set("CRUNIT2",p['CRUNIT2'],"Unit of declination co-ordinates")

    hdr.set("CRPIX1",p['CRPIX1'],"[pixel] Reference pixel along axis 1")
    hdr.set("CRPIX2",p['CRPIX2'],"[pixel] Reference pixel along axis 2")
    hdr.set("CRVAL1",radeg,"[deg] Right ascension at the reference pixel")
    hdr.set("CRVAL2",decdeg,"[deg] Declination at the reference pixel")

    hdr.set("CD1_1",p['CD1_1'],"Transformation matrix element")
    hdr.set("CD1_2",p['CD1_2'],"Transformation matrix element")
    hdr.set("CD2_1",p['CD2_1'],"Transformation matrix element")
    hdr.set("CD2_2",p['CD2_2'],"Transformation matrix element")
    """
    
    # set time zone
    timeZone = TimeDelta(p['TIMEZONE']*u.hour,scale='tai')
    hdr.set("TIMEZONE",p['TIMEZONE'],"time zone in hours e.g. BRT = UT - 3 hr")
    
    # set obstime
    obstime_str = ''
    if "DATE-OBS" in hdr.keys() :
        obstime_str = hdr["DATE-OBS"]
    elif "DATE" in hdr.keys() :
        obstime_str = hdr["DATE"]
    else :
        print("WARNING: could not get time keyword")
        
    if "," in obstime_str :
        obstime_str = obstime_str.replace(",",".")
        
    obstime = Time(obstime_str, format='isot', scale='utc', location=opd_location)
    
    if p['TIMETYPE'] == "LT" :
        obstime = obstime - timeZone
        
    jd = obstime.jd
    mjd = obstime.mjd
    
    # Set light travel time for source observed
    ltt_bary = obstime.light_travel_time(source)
    bjd = obstime.tdb.jd + ltt_bary
    
    #### HJD
    ltt_helio = obstime.light_travel_time(source, 'heliocentric') ### para o HJD
    hjd = obstime.tdb.jd + ltt_helio

    hdr.set("DATE",obstime.isot,"UT date file creation ISOT")
    hdr.set("DATE-OBS",obstime.isot,"UT date at start of exposure ISOT")
    hdr.set("UTDATE",obstime.isot,"UT date at start of exposure ISOT")
    hdr.set("LTDATE",(obstime+timeZone).isot,"LT date at start of exposure ISOT")
    hdr.set("JD",jd,"Julian date at start of exposure")
    hdr.set("JD-OBS",jd,"Julian date at start of exposure")
    hdr.set("MJD",mjd,"Modified Julian date at start of exposure")
    hdr.set("MJD-OBS",mjd,"Modified Julian date at start of exposure")
    hdr.set("BJD",bjd.value,"Barycentric Julian date at start of exposure")
    hdr.set("HJD",hjd.value,"Heliocentric Julian date at start of exposure")

    #sidereal = obstime.sidereal_time('apparent')
    #hdr.set("ST",sidereal,"Sidereal time")
    #hdr.set("SD",sidereal,"Sidereal time")

    # calculate exptime:
    exptime = convert_str2float_with_comma(hdr["EXPTIME"])
    hdr.set("EXPTIME",exptime,"Exposure time (s)")

    # calculate airmass
    airmass = source.transform_to(AltAz(obstime=obstime,location=opd_location)).secz
    
    hdr.set("AIRMASS",airmass.value,"Airmass at start of observation")

    #w = WCS(hdr)
    #sky = w.pixel_to_world(30, 40)
    #print(sky)

    return hdr


def fix_image_orientation (p, img_data, hdr) :

    p['ROTATED'] = 0
    p['FLIPPED'] = "None"

    if p['NORTH'] == 'UP' or p['NORTH'] == 'DOWN':
        if p['NORTH'] == 'UP' and p['EAST'] == 'LEFT' :
            pass
        elif p['NORTH'] == 'UP' and p['EAST'] == 'RIGHT' :
            img_data = np.fliplr(img_data)
            p['FLIPPED'] = "HORIZONTAL"
        elif p['NORTH'] == 'DOWN' and p['EAST'] == 'LEFT' :
            img_data = np.flipud(img_data)
            p['FLIPPED'] = "VERTICAL"
        elif p['NORTH'] == 'DOWN' and p['EAST'] == 'RIGHT' :
            img_data = np.rot90(img_data, k=2)
            p['ROTATED'] = 180
    elif p['NORTH'] == 'RIGHT' or p['NORTH'] == 'LEFT':
        if p['NORTH'] == 'RIGHT' and p['EAST'] == 'UP' :
            img_data = np.rot90(img_data, k=3)
            p['ROTATED'] = 270
        elif p['NORTH'] == 'RIGHT' and p['EAST'] == 'DOWN' :
            img_data = np.rot90(img_data)
            img_data = np.flipud(img_data)
            p['ROTATED'] = 90
            p['FLIPPED'] = "VERTICAL"
        elif p['NORTH'] == 'LEFT' and p['EAST'] == 'UP' :
            img_data = np.rot90(img_data)
            img_data = np.fliplr(img_data)
            p['ROTATED'] = 90
            p['FLIPPED'] = "HORIZONTAL"
        elif p['NORTH'] == 'LEFT' and p['EAST'] == 'DOWN' :
            img_data = np.rot90(img_data)
            p['ROTATED'] = 90

    hdr.set("ROTATED",p['ROTATED'],"Clockwise rotation applied to orig. img. [DEG]")
    hdr.set("FLIPPED",p['FLIPPED'],"Flip applied to original image (after rotation)")

    return img_data, p, hdr


parser = OptionParser()
parser.add_option("-r", "--datadir", dest="datadir", help="data directory",type='string',default="./")
parser.add_option("-i", "--input", dest="input", help="input data pattern",type='string',default="*.fits")
parser.add_option("-o", "--object", dest="object", help="Object name",type='string',default="")
parser.add_option("-y", "--obstype", dest="obstype", help="OBSTYPE",type='string',default="")
parser.add_option("-f", "--filter", dest="filter", help="Filter",type='string',default="")
parser.add_option("-R", "--ra", dest="ra", help="Right Ascension (HH:MM:SS.SSSSS)",type='string',default="")
parser.add_option("-D", "--dec", dest="dec", help="Declination (DD:MM:SS.SSSSS)",type='string',default="")
parser.add_option("-b", "--observer", dest="observer", help="Observer(s)",type='string',default="")

parser.add_option("-t", "--telescope", dest="telescope", help="Telescope",type='string',default="")
parser.add_option("-m", "--instrument", dest="instrument", help="Instrument",type='string',default="")
parser.add_option("-a", "--acqsys", dest="acqsys", help="Acquisition system",type='string',default="")
parser.add_option("-g", "--gain", dest="gain", help="CCD gain [e-/ADU]",type='float',default=0.)
parser.add_option("-n", "--readnoise", dest="readnoise", help="Readout noise [e-]",type='float',default=0.)
parser.add_option("-N", "--north", dest="north", help="North direction [UP, DOWN, RIGHT, LEFT]",type='string',default="")
parser.add_option("-E", "--east", dest="east", help="East direction [UP, DOWN, RIGHT, LEFT]",type='string',default="")

parser.add_option("-T", "--timekey", dest="timekey", help="Time key",type='string',default="")
parser.add_option("-Y", "--timetype", dest="timetype", help="Time type [UT or LT]",type='string',default="")
parser.add_option("-z", "--timezone", dest="timezone", help="Time zone",type='int',default=9999)

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

if options.object != "" :
    p['OBJECT'] = options.object
if options.obstype != "" :
    p['OBSTYPE'] = options.obstype

if options.filter != "" :
    p['FILTER'] = options.filter
  
if options.ra != "" :
    p['RA'] = options.ra
if options.ra != "" :
    p['DEC'] = options.dec

if options.observer != "" :
    p['OBSERVER'] = options.observer

if options.telescope != "" :
    p['TELESCOP'] = options.telescope
if options.instrument != "" :
    p['INSTRUME'] = options.instrument
if options.acqsys != "" :
    p['ACQSYS'] = options.acqsys
if options.gain != 0. :
    p['GAIN'] = options.gain
if options.readnoise != 0. :
    p['RDNOISE'] = options.readnoise
if options.north != "" :
    p['NORTH'] = options.north
if options.east != "" :
    p['EAST'] = options.east


if options.timekey != "" :
    p['TIMEKEY'] = options.timekey
if options.timetype != "" :
    p['TIMETYPE'] = options.timetype
if options.timezone != 9999 :
    p['TIMEZONE'] = options.timezone

for i in range(len(inputdata)) :

    print("Converting image {}/{} -> {}".format(i+1,len(inputdata),inputdata[i]))

    with fits.open(inputdata[i], mode='update') as hdu:
    #with fits.open(inputdata[i]) as hdu:

        # get main header
        header = hdu[0].header
        
        # set keywords related to observatory, instruments, configuration etc.
        header = set_siteobs_keys(p, header)
        
        # set keywords related to time and coordinates
        header = set_timecoords_keys(p, header)

        # get image data
        img_data = hdu[0].data
        
        # check if image data is in a cube format and make it 2D
        if len(np.shape(img_data)) == 3 :
            img_data = img_data[0]
        
        if len(np.shape(img_data)) != 2 :
            print("ERROR: image {} has shape {}, exiting ...".format(inputdata[i], np.shape(img_data)))
            exit()

        # fix image orientation to make sure North is up and East is left.
        img_data, p, header = fix_image_orientation(p, img_data, header)

        # set header and data to original hdu and update file
        hdu[0].header = header
        hdu[0].data = img_data
        hdu.flush()  # changes are written back to original fits


os.chdir(currdir)

