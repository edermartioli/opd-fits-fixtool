# OPD FITS `fixtool`
Toolkit for correcting FITS files obtained from the Pico dos Dias observatory (OPD)

The main tool is `fixrawdata.py`, which can be run from a command line directly from a shell terminal. 

The main functionality of this tool is to edit the FITS header of a given list of FITS files to include the necessary keywords defined in the `fix_params.py` file. The values of these keywords are also defined in the `fix_params.py`, but for a subset of these keywords it is possible to override their values with input values from the command line. The tool also rotates and flips the image to ensure the orientation is set to North up and East to the left. The following input arguments are available:

```
"-r", "--datadir", help="Data directory", type='string'
"-i", "--input", help="Input data pattern", type='string'
"-o", "--object", help="Object name", type='string'
"-y", "--obstype", help="OBSTYPE", type='string'
"-f", "--filter", help="Filter", type='string'
"-R", "--ra", help="Right Ascension (HH:MM:SS.SSSSS)", type='string'
"-D", "--dec", help="Declination (DD:MM:SS.SSSSS)", type='string'
"-b", "--observer", help="Observer(s)", type='string'

"-t", "--telescope", help="Telescope", type='string'
"-m", "--instrument", help="Instrument", type='string'
"-a", "--acqsys", help="Acquisition system", type='string'
"-g", "--gain", help="CCD gain [e-/ADU]", type='float'
"-n", "--readnoise", help="Readout noise [e-]", type='float'
"-N", "--north", help="North direction [UP, DOWN, RIGHT, LEFT]", type='string'
"-E", "--east", help="East direction [UP, DOWN, RIGHT, LEFT]", type='string'

"-T", "--timekey", help="Time key", type='string'
"-Y", "--timetype", help="Time type [UT or LT]", type='string'
"-z", "--timezone", help="Time zone", type='int'
```
