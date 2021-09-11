![image](https://github.com/patternizer/ncc-stripes/blob/main/title_frame.png)

# ncc-stripes

Python code to merge geological, instrumental and climate model land surface air temperature anomalies and visualise them in the form of the [Climate Stripes](https://showyourstripes.info/) by lead scientist: Professor Ed Hawkins and uses data from PAGES2k, HadCRUT5 and FaIR v1.6.3.
  
## Contents

* `plot_ncc_stripes.py` - main script to be run with Python 3.8.11+

The first step is to clone the latest ncc-stripes code and step into the check out directory: 

    $ git clone https://github.com/patternizer/ncc-stripes.git
    $ cd ncc-stripes
    
### Using Standard Python 

The code should run with the [standard CPython](https://www.python.org/downloads/) installation and was tested 
in a conda virtual environment running a 64-bit version of Python 3.8.11+.

ncc-stripes can be run from sources directly, once the following data requirements are present in the DATA/ directory:

* `GloSAT.prelim01_reqSD_alternativegrid-178101-201912.timeseries.txt` - instrumental temperatures

Run with:

    $ python plot_ncc_stripes.py
        
## License

The code is distributed under terms and conditions of the Attribution 4.0 International (CC BY 4.0) license: 
https://creativecommons.org/licenses/by/4.0/ and is designed to visualise land surface air temperature data 
and provide an implementation akin to https://showyourstripes.info/ - graphics and lead scientist: Professor Ed Hawkins.

## Contact information

* [Michael Taylor](https://patternizer.github.io)


