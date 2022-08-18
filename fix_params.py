# -----------------------------------------------------------------------------
#   define OBS FIX parameters
# -----------------------------------------------------------------------------
def load_obsfix_parameters(obstype='OBJECT') :
    
    #initialize parameters dictionary
    p = {}
    
    # SWITCHES
    p['UPDATE_FILE_NAME'] = False
    p['GET_SIMBAD_COORDS'] = False
    p['RECENTER_ON_TARGET'] = False
    p['TARGET_SUB_REGION'] = [0,1024,0,1024]
    p['GETCOORDSFROMHEADER'] = False

    # VARIABLES
    if obstype == 'OBJECT' :
        p['OBSTYPE'] = 'OBJECT'
        p['OBJECT'] = 'AU Mic'
        p['RA'] = '20:45:09.5324974119'
        p['DEC'] = '-31:20:27.237889841'
        p['CRPIX1'] = 435  # x pixel of target position
        p['CRPIX2'] = 941  # y pixel of target position
        p['CDELT1'] = 0.34 / (60 * 60)
        p['CDELT2'] = 0.34 / (60 * 60)

        p['SCIENCE'] = True

    elif obstype == 'FLAT' or obstype == 'ZERO':
        p['OBSTYPE'] = obstype
        p['OBJECT'] = obstype
        p['SCIENCE'] = False
        
    else :
        print("OBSTYPE = {} not recognized. Use OBJECT, FLAT or ZERO; exiting ...".format(obstype))
        exit()
        
    p['FILTER'] = 'I'
    p['OBSERVER'] = 'Eder, Laerte, Ted, Denis, Nelio'
    p['TELESCOP'] = 'BC0.60m' # 'BC0.60m', 'ZE0.60m', 'PE1.60m'
    
    p['INSTRUME'] = 'CAM2+IxonUltra9915'
    p['ACQSYS'] = 'OPDAcquisition'
    #p['ACQSYS'] = 'S4ACS'
    p['GAIN'] = 3.3
    p['RDNOISE'] = 6.57

    p['NORTH'] = 'RIGHT' # UP, DOWN, RIGHT, LEFT
    p['EAST'] = 'DOWN' # UP, DOWN, RIGHT, LEFT
    
    # TIME KEYWORD
    p['TIMEKEY'] = 'DATE'
    p['TIMETYPE'] = 'UT' # UT or LT
    p['TIMEZONE'] = -3
    
    # CONSTANTS
    p['CRAZYVALUE'] = -9999
    p['OBSALT'] = 1864.0
    p['OBSLONG'] = -45.5825
    p['OBSLAT'] = -22.53444444444445
    p['EQUINOX'] = 2000.0

    p['CTYPE1'] = 'RA---TAN'
    p['CTYPE2'] = 'DEC--TAN'

    p['CRUNIT1'] = 'deg     '
    p['CRUNIT2'] = 'deg     '
    
    p['CD1_1'], p['CD1_2'] = 1.0 , 0.0
    p['CD2_1'], p['CD2_2'] = 0.0 , 1.0

    #### KEYWORDS #####
    
    p['IMPORTANT_KEYS'] = [ 'OBSTYPE',
                        'OBJECT',
#                        'OBSERVER',
                        'FILTER',
                        'TELESCOP',
                        'INSTRUME',
                        'ACQSYS',
                        'EXPTIME',
#                        'OBSALT',
#                        'OBSLONG',
#                        'OBSLAT',
#                        'EQUINOX',
                        'DATE-OBS',
#                        'UTDATE',
                        'TIMEZONE',
#                        'LTDATE',
#                        'JD-OBS',
#                        'MJD-OBS',
                        'BJD',

#                        'RA',
                        'RA_DEG',
#                        'DEC',
                        'DEC_DEG',
#                        'EPOCH',
#                        'HA',
#                        'HA_DEG',
#                        'ST',
                        'AIRMASS',
#                        'CTYPE1',
#                        'CTYPE2',
#                        'CRUNIT1',
#                        'CRUNIT2',
#                        'CRPIX1',
#                        'CRPIX2',
#                        'CRVAL1',
#                        'CRVAL2',
#                        'CDELT1',
#                        'CDELT2',
#                        'CD1_1',
#                        'CD1_2',
#                        'CD2_1',
#                        'CD2_2',
#                        'PV2_1',
#                        'PV2_2',
#                        'PV2_3',
#                        'CCDTEMP',
                        'GAIN',
                        'RDNOISE',
#                        'CAMFOC',
#                        'TEMPEXT',
#                        'PRESSURE',
#                        'HUMIDITY'
                    ]
    return p

'''
CTYPE1  = 'RA---ZPN'           / Algorithm type for axis 1
CTYPE2  = 'DEC--ZPN'           / Algorithm type for axis 2
CRPIX1  =        2.9950000E+03 / [pixel] Reference pixel along axis 1
CRPIX2  =       -9.7296002E+02 / [pixel] Reference pixel along axis 2
CRVAL1  =        1.5847932E+02 / [deg] Right ascension at the reference pixel
CRVAL2  =        5.5045403E+01 / [deg] Declination at the reference pixel
CRUNIT1 = 'deg     '           / Unit of right ascension co-ordinates
CRUNIT2 = 'deg     '           / Unit of declination co-ordinates
CD1_1   =       -8.6063899E-08 / Transformation matrix element
CD1_2   =       -1.1131840E-04 / Transformation matrix element
CD2_1   =        1.1132305E-04 / Transformation matrix element
CD2_2   =       -1.8366745E-07 / Transformation matrix element
PV2_1   =             1.00E+00 / Pol.coeff. for pixel -> celestial coord
PV2_2   =         0.000000E+00 / Pol.coeff. for pixel -> celestial coord
PV2_3   =            -5.00E+01 / Pol.coeff. for pixel -> celestial coord
'''
