![image](https://github.com/patternizer/ncc-stripes/blob/main/title_frame.png)

# ncc-stripes

Python code to merge geological, instrumental and climate model land surface air temperature anomalies and visualise them in the form of the [Climate Stripes](https://showyourstripes.info/) by lead scientist: Professor Ed Hawkins and uses data from PAGES2k, HadCRUT5 and FaIR v1.6.3.
  
## Contents

* `plot_ncc_stripes_futures.py` - climate stripes 500-2200 CE with 5 climate model forecasts from FaIR: SSP1-1.9, SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5
* `plot_ncc_stripes_pleistocene.py` - climate stripes 2.58 Myr ( < 2015 ) to 2200 CE for user-choice of SSP
* `plot_ncc_stripes_multiaxis.py` - climate stripes 65.5 Myr ( < 2015 ) to 2200 CE with epochs plotted in segments for user-choice of SSP
* `plot_ncc_stripes_rebinning.py` - climate stripes, bars and timeseries 65.5 Myr ( < 2015 ) to 2200 CE with log-linear re-binning implementation for user-choice of SSP

The first step is to clone the latest ncc-stripes code and step into the check out directory: 

    $ git clone https://github.com/patternizer/ncc-stripes.git
    $ cd ncc-stripes
    
### Using Standard Python 

The code should run with the [standard CPython](https://www.python.org/downloads/) installation and was tested 
in a conda virtual environment running a 64-bit version of Python 3.8.11+.

ncc-stripes can be run from sources directly, once the following data requirements are present in the DATA/ directory:

* `paleo_data_compilation.xls` - multi-proxy temperature estimate 65.5 Myr ( < 2015 ) to 1 CE ( source: [Greg Fergus](https://gergs.net/all_palaeotemps/), see also the discussion by Gavin Schmidt at [RealClimate](https://www.realclimate.org/index.php/archives/2014/03/can-we-make-better-graphs-of-global-temperature-history/) )
* `PAGES2k.txt` - 7000-member ensemble median global mean temperature reconstruction 1-1850 CE ( source: [PAGES2k Consortium](https://figshare.com/articles/dataset/Reconstruction_ensemble_median_and_95_range/8143094) )
* `HadCRUT5.csv` - instrumental record 1850-2020 CE ( source: [UKMO Hadley Centre](https://www.metoffice.gov.uk/hadobs/hadcrut5/) )
* `variability_realisation0.txt` - 1350-1546 proxy for climate model internal variability ( source: Professor Tim Osborn )
* `variability_realisation1.txt` - 1550-1746 proxy for climate model internal variability ( source: Professor Tim Osborn )
* `fair_rcp3pd.csv` - FaIR v1.6.3 climate model forecast for 2021-2200 CE for RCP3pd using an ensemble median constrained by HadCRUT5 observations ( source: [FaIR v1.6.3](https://github.com/OMS-NetZero/FAIR) )
* `fair_rcp45.csv` - FaIR v1.6.3 climate model forecast for 2021-2200 CE for RCP4.5using an ensemble median constrained by HadCRUT5 observations ( source: [FaIR v1.6.3](https://github.com/OMS-NetZero/FAIR) )
* `fair_rcp6.csv` - FaIR v1.6.3 climate model forecast for 2021-2200 CE for RCP6 using an ensemble median constrained by HadCRUT5 observations ( source: [FaIR v1.6.3](https://github.com/OMS-NetZero/FAIR) )
* `fair_rcp85.csv` - FaIR v1.6.3 climate model forecast for 2021-2200 CE for RCP8.5 using an ensemble median constrained by HadCRUT5 observations ( source: [FaIR v1.6.3](https://github.com/OMS-NetZero/FAIR) )
* `fair_ssp119.csv` - FaIR v1.6.3 climate model forecast for 2021-2200 CE for SSP1-1.9 using an ensemble median constrained by HadCRUT5 observations ( source: [FaIR v1.6.3](https://github.com/OMS-NetZero/FAIR) )
* `fair_ssp126.csv` - FaIR v1.6.3 climate model forecast for 2021-2200 CE for SSP1-2.6 using an ensemble median constrained by HadCRUT5 observations ( source: [FaIR v1.6.3](https://github.com/OMS-NetZero/FAIR) )
* `fair_ssp246.csv` - FaIR v1.6.3 climate model forecast for 2021-2200 CE for SSP2-4.5 using an ensemble median constrained by HadCRUT5 observations ( source: [FaIR v1.6.3](https://github.com/OMS-NetZero/FAIR) )
* `fair_ssp370.csv` - FaIR v1.6.3 climate model forecast for 2021-2200 CE for SSP3-7.0 using an ensemble median constrained by HadCRUT5 observations ( source: [FaIR v1.6.3](https://github.com/OMS-NetZero/FAIR) )
* `fair_ssp585.csv` - FaIR v1.6.3 climate model forecast for 2021-2200 CE for SSP5-8.5 using an ensemble median constrained by HadCRUT5 observations ( source: [FaIR v1.6.3](https://github.com/OMS-NetZero/FAIR) )

Run with for example:

    $ python plot_ncc_stripes_futures.py
    $ python plot_ncc_stripes_pleistocene.py
    $ python plot_ncc_stripes_multiaxis.py
    $ python plot_ncc_stripes_rebinning.py
        
## License

The code is distributed under terms and conditions of the Attribution 4.0 International (CC BY 4.0) license: 
https://creativecommons.org/licenses/by/4.0/ and is designed to visualise land surface air temperature data 
and provide an implementation akin to https://showyourstripes.info/ by Professor Ed Hawkins.

## Contact information

* [Michael Taylor](https://patternizer.github.io)


