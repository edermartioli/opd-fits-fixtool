# OPD FITS `fixtool`
Toolkit for correcting FITS files obtained from the Pico dos Dias observatory (OPD)

The main tool is `fixrawdata.py`, which can be run from a command line directly from a shell terminal. 

The main functionality of this tool is to edit the FITS header of a given list of FITS files to include the necessary keywords defined in the `fix_params.py` file. The values of these keywords are also defined in the `fix_params.py`, but for a subset of these keywords it is possible to override their values with input values from the command line. The tool also rotates and flips the image to ensure the orientation is set to North up and East to the left. The following input arguments are available:

```
"-r", "--datadir", dest="datadir", help="data directory",type='string'
"-i", "--input", dest="input", help="input data pattern",type='string'
"-o", "--object", dest="object", help="Object name",type='string'
"-y", "--obstype", dest="obstype", help="OBSTYPE",type='string'
"-f", "--filter", dest="filter", help="Filter",type='string'
"-R", "--ra", dest="ra", help="Right Ascension (HH:MM:SS.SSSSS)",type='string'
"-D", "--dec", dest="dec", help="Declination (DD:MM:SS.SSSSS)",type='string'
"-b", "--observer", dest="observer", help="Observer(s)",type='string'

"-t", "--telescope", dest="telescope", help="Telescope",type='string'
"-m", "--instrument", dest="instrument", help="Instrument",type='string'
"-a", "--acqsys", dest="acqsys", help="Acquisition system",type='string'
"-g", "--gain", dest="gain", help="CCD gain [e-/ADU]",type='float'
"-n", "--readnoise", dest="readnoise", help="Readout noise [e-]",type='float'
"-N", "--north", dest="north", help="North direction [UP, DOWN, RIGHT, LEFT]",type='string'
"-E", "--east", dest="east", help="East direction [UP, DOWN, RIGHT, LEFT]",type='string'

"-T", "--timekey", dest="timekey", help="Time key",type='string'
"-Y", "--timetype", dest="timetype", help="Time type [UT or LT]",type='string'
"-z", "--timezone", dest="timezone", help="Time zone",type='int'
```
