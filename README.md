# OPD FITS `fixtool`
Toolkit for correcting FITS files obtained from the Pico dos Dias observatory (OPD)

# Tool `fixrawdata.py`

The main tool is `fixrawdata.py`, which can be run from a command line directly from a shell terminal. 

The main functionality of this tool is to edit the FITS header of a given list of FITS files to include the necessary keywords defined in the `fix_params.py` file. The values of these keywords are also defined in the `fix_params.py`, but for a subset of these keywords it is possible to override their values with input values from the command line. The tool also rotates and flips the image to ensure the orientation is set to North up and East to the left. The following input arguments are available:

```
"-r", "--datadir", help="Data directory", type='string'
"-i", "--input", help="Input data pattern", type='string'

"-o", "--object", help="Object name", type='string'
"-y", "--obstype", help="Observation type [ZERO, FLAT, OBJECT]", type='string'
"-f", "--filter", help="Filter [U, B, V, R, I]", type='string'
"-b", "--observer", help="Observer(s)", type='string'
"-t", "--telescope", help="Telescope ['BC0.60m', 'ZE0.60m', 'PE1.60m']", type='string'
"-m", "--instrument", help="Instrument", type='string'
"-a", "--acqsys", help="Acquisition system", type='string'

"-g", "--gain", help="CCD gain [e-/ADU]", type='float'
"-n", "--readnoise", help="Readout noise [e-]", type='float'

"-R", "--ra", help="Right Ascension (HH:MM:SS.SSSSS)", type='string'
"-D", "--dec", help="Declination (DD:MM:SS.SSSSS)", type='string'

"-N", "--north", help="North direction [UP, DOWN, RIGHT, LEFT]", type='string'
"-E", "--east", help="East direction [UP, DOWN, RIGHT, LEFT]", type='string'

"-T", "--timekey", help="Time key", type='string'
"-Y", "--timetype", help="Time type [UT or LT]", type='string'
"-z", "--timezone", help="Time zone", type='int'
```

Below is a simple usage example using the default settings defined in the `fix_params.py` file:

```
python $PATH/fixrawdata.py --input=*.fits
```

Below are some usage examples where some keywords are defined on the command line:

```
ZE_22jul17_E -> python $PATH/fixrawdata.py --input=aumic_i_*.fits --telescope='ZE0.60m' --instrument='CAM4+Ixon4269' --acqsys="OPDAcquisition" --gain=3.8 --readnoise=8.2 --north='RIGHT' --east='DOWN' --timekey="DATE-OBS" --timetype="UT" --timezone=-3 -v

ZE_22jul17_W -> python $PATH/fixrawdata.py --input=aumic_i_*.fits --telescope='ZE0.60m' --instrument='CAM4+Ixon4269' --acqsys="OPDAcquisition" --gain=3.8 --readnoise=8.2 --north='LEFT' --east='UP' --timekey="DATE-OBS" --timetype="UT" --timezone=-3 -v
```

BC0.60m -> 22jul05:
```
Science: python $PATH/fixrawdata.py --input=*_c1_AUMIC_I.fits --telescope='BC0.60m' --instrument='CAM2+IxonUltra9915' --acqsys="GEI+ACS" --gain=3.3 --readnoise=6.57 --north='RIGHT' --east='DOWN' --timekey="DATE" --timetype="UT" --timezone=-3 -v

Bias: python $PATH/fixrawdata.py --input=*_c1_zero.fits --object="ZERO" --obstype="ZERO" --filter="None" --ra="00:00:00.00" --dec="00:00:00.00" --telescope='BC0.60m' --instrument='CAM2+IxonUltra9915' --acqsys="GEI+ACS" --gain=3.3 --readnoise=6.57 --north='RIGHT' --east='DOWN' --timekey="DATE" --timetype="UT" --timezone=-3 -v

Flat: python $PATH/fixrawdata.py --input=*_c1_flat_I.fits --object="FLAT" --obstype="FLAT" --filter="I" --ra="00:00:00.00" --dec="00:00:00.00" --telescope='BC0.60m' --instrument='CAM2+IxonUltra9915' --acqsys="GEI+ACS" --gain=3.3 --readnoise=6.57 --north='RIGHT' --east='DOWN' --timekey="DATE" --timetype="UT" --timezone=-3 -v
```

*****************************************
# WARNING: `fixrawdata.py` updates FITS files, so changes made cannot be undone!!!
*****************************************


# Tool `checkkeys.py`

The tool `checkkeys.py` can be used to print the important keys of a list of input FITS files. The important keys are defined as a list of keywords saved in the dictionary entry `p['IMPORTANT_KEYS']` in the `fix_params.py` file. One can uncomment/comment keywords in this list to include/exclude important keys to be printed. 

Below is a simple usage example of the tool `checkkeys.py`:

```
python $PATH/checkkeys.py --input=*.fits
```

