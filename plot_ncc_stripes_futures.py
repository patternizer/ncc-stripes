#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: ncc_stripes_futures.py
# ( plots PAGES2k + HadCRUT5 + FaIR to visualise futures: no re-binning )
#------------------------------------------------------------------------------
# Version 0.3
# 30 September, 2021
# Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
# michael DOT a DOT taylor AT uea DOT ac DOT uk
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# IMPORT PYTHON LIBRARIES
#------------------------------------------------------------------------------

# Dataframe libraries:
import numpy as np
import pandas as pd
import xarray as xr

# Datetime libraries:
from datetime import datetime
import nc_time_axis
import cftime
from cftime import num2date, DatetimeNoLeap

# Plotting libraries:
import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt; plt.close('all')
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.cm import ScalarMappable
from matplotlib import rcParams
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

# Silence library version notifications
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# SETTINGS: 
#------------------------------------------------------------------------------

fontsize = 10
nsmooth = 2 # years
cbar_max = 6.0
t_start = 500 # >= 0

use_dark_theme = True
use_smoothing = True
use_overlay_axis = True
use_overlay_timeseries = True
use_overlay_colorbar = True

#projectionstr = 'RCP3pd'
#projectionstr = 'RCP45'
#projectionstr = 'RCP6'
#projectionstr = 'RCP85'
#projectionstr = 'SSP119'
#projectionstr = 'SSP126'
#projectionstr = 'SSP245'
#projectionstr = 'SSP370'
#projectionstr = 'SSP585'

projectionstr1 = 'SSP119'
projectionstr2 = 'SSP126'
projectionstr3 = 'SSP245'
projectionstr4 = 'SSP370'
projectionstr5 = 'SSP585'

baselinestr = 'baseline_1851_1900'
#baselinestr = 'baseline_1961_1990'
#baselinestr = 'baseline_1971_2000'

pathstr = 'DATA/'
pages2kstr = 'PAGES2k.txt'
hadcrut5str = 'HadCRUT5.csv'
fairstr1 = 'fair' + '_' + projectionstr1.lower() + '.csv' 
fairstr2 = 'fair' + '_' + projectionstr2.lower() + '.csv' 
fairstr3 = 'fair' + '_' + projectionstr3.lower() + '.csv' 
fairstr4 = 'fair' + '_' + projectionstr4.lower() + '.csv' 
fairstr5 = 'fair' + '_' + projectionstr5.lower() + '.csv' 
lovarstr = 'variability_realisation0.txt' # via Tim
hivarstr = 'variability_realisation1.txt' # via Tim
paleostr = 'paleo_data_compilation.xls'
pages2k_file = pathstr + pages2kstr 
hadcrut5_file = pathstr + hadcrut5str 
fair_file1 = pathstr + fairstr1 
fair_file2 = pathstr + fairstr2 
fair_file3 = pathstr + fairstr3 
fair_file4 = pathstr + fairstr4 
fair_file5 = pathstr + fairstr5 
lo_var_file = pathstr + lovarstr 
hi_var_file = pathstr + hivarstr 
paleo_file = pathstr + paleostr 

ipcc_rgb_txtfile = np.loadtxt("DATA/temp_div.txt") # IPCC AR6 temp div colormap file
cmap = mcolors.LinearSegmentedColormap.from_list('colormap', ipcc_rgb_txtfile) # ipcc_colormap
#cmap = plt.cm.get_cmap('RdBu_r')
#cmap = plt.cm.get_cmap('bwr')
                
#------------------------------------------------------------------------------
# DARK THEME
#------------------------------------------------------------------------------

if use_dark_theme == True:
    
    matplotlib.rcParams['text.usetex'] = True
#    rcParams['font.family'] = ['DejaVu Sans']
#    rcParams['font.sans-serif'] = ['Avant Garde']
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Avant Garde', 'Lucida Grande', 'Verdana', 'DejaVu Sans' ]
    plt.rc('text',color='white')
    plt.rc('lines',color='white')
    plt.rc('patch',edgecolor='white')
    plt.rc('grid',color='lightgray')
    plt.rc('xtick',color='white')
    plt.rc('ytick',color='white')
    plt.rc('axes',labelcolor='white')
    plt.rc('axes',facecolor='black')
    plt.rc('axes',edgecolor='lightgray')
    plt.rc('figure',facecolor='black')
    plt.rc('figure',edgecolor='black')
    plt.rc('savefig',edgecolor='black')
    plt.rc('savefig',facecolor='black')
    
else:

    matplotlib.rcParams['text.usetex'] = True
#    rcParams['font.family'] = ['DejaVu Sans']
#    rcParams['font.sans-serif'] = ['Avant Garde']
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Avant Garde', 'Lucida Grande', 'Verdana', 'DejaVu Sans' ]
    plt.rc('text',color='black')
    plt.rc('lines',color='black')
    plt.rc('patch',edgecolor='black')
    plt.rc('grid',color='lightgray')
    plt.rc('xtick',color='black')
    plt.rc('ytick',color='black')
    plt.rc('axes',labelcolor='black')
    plt.rc('axes',facecolor='white')    
    plt.rc('axes',edgecolor='black')
    plt.rc('figure',facecolor='white')
    plt.rc('figure',edgecolor='white')
    plt.rc('savefig',edgecolor='white')
    plt.rc('savefig',facecolor='white')

# Calculate current time

now = datetime.now()
currentdy = str(now.day).zfill(2)
currentmn = str(now.month).zfill(2)
currentyr = str(now.year)
titletime = str(currentdy) + '/' + currentmn + '/' + currentyr

#------------------------------------------------------------------------------
# LOAD: PAGES2k (via Ed Hawkins with thanks) --> df_pages2k
# NB: convert time to year.decimal
#------------------------------------------------------------------------------

# FORMAT:
# Year CE | raw instrumental target data | reconstruction ensemble 50th | 2.5th | 97.5th percentiles | 
# 31-year butterworth filtered instrumental target data | 31-year butterworth filtered reconstruction 50th | 
# 2.5th | 97.5th percentiles

nheader = 5
f = open(pages2k_file)
lines = f.readlines()
years = [] # [0001,2000]
obs = []
for i in range(nheader,len(lines)):
        words = lines[i].split()   
        year = words[0].zfill(4)
        val = (len(words)-1)*[None]            
        for j in range(len(val)):                                
            try: val[j] = float(words[j+1])                
            except: 
                pass                                 
        years.append(year)                                     
        obs.append(val)            
f.close()    
obs = np.array(obs)
t_pages2k = xr.cftime_range(start=years[0], periods=len(years), freq='A', calendar='gregorian')[0:1849]
ts_pages2k_instr = pd.to_numeric(obs[:,1][0:1849], errors='coerce')
ts_pages2k_recon = pd.to_numeric(obs[:,5][0:1849], errors='coerce')
ts_pages2k = np.append(ts_pages2k_recon[0:-36],ts_pages2k_instr[-36:],axis=None)
df_pages2k = pd.DataFrame()
df_pages2k['t_pages2k'] = t_pages2k.year.astype(float)
df_pages2k['ts_pages2k'] = ts_pages2k

#------------------------------------------------------------------------------
# LOAD: HadCRUT5 (via Tim Osborn and UKMO with thanks) --> df_hadcrut5
# NB: convert time to year.decimal
#------------------------------------------------------------------------------

hadcrut5 = pd.read_csv(hadcrut5_file)
t_hadcrut5_monthly = xr.cftime_range(start='1850', periods=len(hadcrut5), freq='MS', calendar='noleap')
ts_hadcrut5_monthly = hadcrut5['Anomaly (deg C)'].values
df_hadcrut5 = pd.DataFrame()
df_hadcrut5['t_hadcrut5'] = t_hadcrut5_monthly.year.astype(float) + t_hadcrut5_monthly.month.astype(float)/12.0
df_hadcrut5['ts_hadcrut5'] = ts_hadcrut5_monthly
years = np.unique(t_hadcrut5_monthly.year)
yearly = []
SD = []
for yyyy in years:
    year_data = df_hadcrut5[np.floor(df_hadcrut5['t_hadcrut5']).astype('int') == yyyy]['ts_hadcrut5']
    yearly_mean = np.nanmean(year_data)
    yearly_SD = np.nanstd(year_data)
    yearly.append(yearly_mean)
    SD.append(yearly_SD)
df_hadcrut5_yearly = pd.DataFrame()
df_hadcrut5_yearly['t_hadcrut5'] = years.astype('float')
df_hadcrut5_yearly['ts_hadcrut5'] = yearly
df_hadcrut5_yearly['ts_hadcrut5_SD'] = SD
df_hadcrut5_yearly = df_hadcrut5_yearly[df_hadcrut5_yearly.t_hadcrut5 <= 2020]

#-----------------------------------------------------------------------------
# LOAD: FaIR v1.6.3 projections (constrained by HadCRUT5-analysis) --> df_fair1, df_fair2 ...
# NB: convert time to year.decimal
#-----------------------------------------------------------------------------

fair1 = pd.read_csv(fair_file1)
fair2 = pd.read_csv(fair_file2)
fair3 = pd.read_csv(fair_file3)
fair4 = pd.read_csv(fair_file4)
fair5 = pd.read_csv(fair_file5)
df_fair1 = pd.DataFrame()
df_fair1['t_fair1'] = fair1.Year.values.astype('float')
df_fair1['ts_fair1'] = fair1.Global.values
df_fair2 = pd.DataFrame()
df_fair2['t_fair2'] = fair2.Year.values.astype('float')
df_fair2['ts_fair2'] = fair2.Global.values
df_fair3 = pd.DataFrame()
df_fair3['t_fair3'] = fair3.Year.values.astype('float')
df_fair3['ts_fair3'] = fair3.Global.values
df_fair4 = pd.DataFrame()
df_fair4['t_fair4'] = fair4.Year.values.astype('float')
df_fair4['ts_fair4'] = fair4.Global.values
df_fair5 = pd.DataFrame()
df_fair5['t_fair5'] = fair5.Year.values.astype('float')
df_fair5['ts_fair5'] = fair5.Global.values

#------------------------------------------------------------------------------
# LOAD: internal variability for FaIR v1.6.4 projections calculated by Tim from the instrumental record.
#------------------------------------------------------------------------------

nheader = 0
f_lo = open(lo_var_file)
f_hi = open(hi_var_file)
lines_lo = f_lo.readlines()
lines_hi = f_hi.readlines()
years = []
obs_lo = []
obs_hi = []
for i in range(nheader,180):
        words_lo = lines_lo[i].split()   
        words_hi = lines_hi[i].split()   
        year = int(words_lo[0].zfill(4)) + 671 # 1350 --> 2021 offset
        val_lo = (len(words_lo)-1)*[None]            
        val_hi = (len(words_hi)-1)*[None]            
        for j in range(len(val_lo)):                                
            try: 
                val_lo[j] = float(words_lo[j+1])                
                val_hi[j] = float(words_hi[j+1])                
            except: 
                pass                                 
        years.append(year)                                     
        obs_lo.append(val_lo)            
        obs_hi.append(val_hi)            
f_lo.close()    
f_hi.close()    
obs_lo = np.array(obs_lo).ravel()
obs_hi = np.array(obs_hi).ravel()
df_variability = pd.DataFrame({'lo_var':obs_lo, 'hi_var':obs_hi}, index=years)

if (projectionstr1 == 'SSP119') | (projectionstr1 == 'SSP126') | (projectionstr1 == 'SSP245') | (projectionstr1 == 'RCP3pd') | (projectionstr1 == 'RCP45'):
    df_fair1['ts_fair1'] = fair1.Global.values + df_variability.lo_var.values
elif (projectionstr1 == 'SSP370') | (projectionstr1 == 'SSP585') | (projectionstr1 == 'RCP6') | (projectionstr1 == 'RCP85'):
    df_fair1['ts_fair1'] = fair1.Global.values + df_variability.hi_var.values
if (projectionstr2 == 'SSP119') | (projectionstr2 == 'SSP126') | (projectionstr2 == 'SSP245') | (projectionstr2 == 'RCP3pd') | (projectionstr2 == 'RCP45'):
    df_fair2['ts_fair2'] = fair2.Global.values + df_variability.lo_var.values
elif (projectionstr2 == 'SSP370') | (projectionstr2 == 'SSP585') | (projectionstr2 == 'RCP6') | (projectionstr2 == 'RCP85'):
    df_fair2['ts_fair2'] = fair2.Global.values + df_variability.hi_var.values
if (projectionstr3 == 'SSP119') | (projectionstr3 == 'SSP126') | (projectionstr3 == 'SSP245') | (projectionstr3 == 'RCP3pd') | (projectionstr3 == 'RCP45'):
    df_fair3['ts_fair3'] = fair3.Global.values + df_variability.lo_var.values
elif (projectionstr3 == 'SSP370') | (projectionstr3 == 'SSP585') | (projectionstr3 == 'RCP6') | (projectionstr3 == 'RCP85'):
    df_fair3['ts_fair3'] = fair3.Global.values + df_variability.hi_var.values
if (projectionstr4 == 'SSP119') | (projectionstr4 == 'SSP126') | (projectionstr4 == 'SSP245') | (projectionstr4 == 'RCP3pd') | (projectionstr4 == 'RCP45'):
    df_fair4['ts_fair4'] = fair4.Global.values + df_variability.lo_var.values
elif (projectionstr4 == 'SSP370') | (projectionstr4 == 'SSP585') | (projectionstr4 == 'RCP6') | (projectionstr4 == 'RCP85'):
    df_fair4['ts_fair4'] = fair4.Global.values + df_variability.hi_var.values
if (projectionstr5 == 'SSP119') | (projectionstr5 == 'SSP126') | (projectionstr5 == 'SSP245') | (projectionstr5 == 'RCP3pd') | (projectionstr5 == 'RCP45'):
    df_fair5['ts_fair5'] = fair5.Global.values + df_variability.lo_var.values
elif (projectionstr5 == 'SSP370') | (projectionstr5 == 'SSP585') | (projectionstr5 == 'RCP6') | (projectionstr5 == 'RCP85'):
    df_fair5['ts_fair5'] = fair5.Global.values + df_variability.hi_var.values

#------------------------------------------------------------------------------
# LOAD: geological anomalies: 65.5229 Myr ( before 2015 )
# NB: paleo_file is from the "data_compilation" sheet in All_palaeotemps.xlsx  
# made available at: hhttp://gergs.net/?attachment_id=4310
#------------------------------------------------------------------------------

xl = pd.ExcelFile(paleo_file)
df_xl = xl.parse('Sheet1',header=2)

# FORMAT:
# Royer et al (2004)							                                                        Friedrich et al (2012) & Hansen et al (2013)				Zachos et al (2008) & Hansen et al (2013)					Lisiecki and Raymo (2005) & Hansen et al (2013)				EPICA Dome C, Antarctica  (x 0.5)		NGRIP, Greenland & Johnsen et al (1989) (x 0.5)				Marcott et al (2013)			Berkeley Earth land-ocean				    IPCC AR5 RCP8.5		
# Age My	 Royer / Veizer (x 2.0)	 Royer / Veizer - CO₂ from proxies (x 2.0)	Low	High	Axis	[]	Age My	Age ky before 2015	δ18O	Tdo	Ts	T anomaly		Age My	Age ky before 2015	δ18O	Tdo	Ts	T anomaly		Age My	Age ky before 2015	δ18O	Tdo	Ts	T anomaly		Age ky before 2015	T	T global		Age ky before 2015	δ18O	Ts	T anomaly	T global		Age ky before 2015	T	1σ		Decade	Age ky before 2015	T average		Year	Age ky before 2015	T

t_epica = df_xl.iloc[:,28] * -1.0e3 + 2015.0
ts_epica = df_xl.iloc[:,30]
t_lisiecki = df_xl.iloc[:,22] * -1.0e3 + 2015.0
ts_lisiecki = df_xl.iloc[:,26]
t_zachos = df_xl.iloc[:,15] * -1.0e3 + 2015.0
ts_zachos = df_xl.iloc[:,19]

# CONCATENATE: epochs and store in dataframe

t_paleo = np.array( list(t_epica) + list(t_lisiecki) + list(t_zachos) ).ravel().astype(int)
ts_paleo = np.array( list(ts_epica) + list(ts_lisiecki) + list(ts_zachos) ).ravel()
df_paleo = pd.DataFrame()
df_paleo['t_paleo'] = t_paleo
df_paleo['ts_paleo'] = ts_paleo

# TRIM: to 0 AD

df_paleo = df_paleo[ df_paleo.t_paleo < 1 ].dropna()

#------------------------------------------------------------------------------
# COMPUTE: climatological monthly mean (Ed's Climate Stripes: 1971-2000) from HadCRUT5 (monthly)
#------------------------------------------------------------------------------

mu_1851_1900 = np.nanmean( df_hadcrut5[(df_hadcrut5['t_hadcrut5']>=1851) & (df_hadcrut5['t_hadcrut5']<=1900)]['ts_hadcrut5'] ) # -0.507873106
mu_1961_1990 = np.nanmean( df_hadcrut5[(df_hadcrut5['t_hadcrut5']>=1961) & (df_hadcrut5['t_hadcrut5']<=1990)]['ts_hadcrut5'] ) #  0.005547222
mu_1971_2000 = np.nanmean( df_hadcrut5[(df_hadcrut5['t_hadcrut5']>=1971) & (df_hadcrut5['t_hadcrut5']<=2000)]['ts_hadcrut5'] ) #  0.176816667

if baselinestr == 'baseline_1851_1900':
    mu = mu_1851_1900
    baseline_start = 1851
    baseline_end = 1900
elif baselinestr == 'baseline_1961_1990':
    mu = mu_1961_1990
    baseline_start = 1961
    baseline_end = 1990
else:    
    mu = mu_1971_2000
    baseline_start = 1971
    baseline_end = 2000
    
baseline_midpoint = baseline_start+int(np.floor(baseline_end-baseline_start)/2)
cbarstr = r'Anomaly, $^{\circ}$C (from ' + str(baseline_start) + '-' + str(baseline_end) +')'

#------------------------------------------------------------------------------
# ALIGN: FaIR projections relative to chosen baseline
# NB: FaIR projections are calculated relative to 1880-2015 emission trends and so we
# need to subtract off the difference between the selected baseline and the 1851-1900 baseline
#------------------------------------------------------------------------------

df_fair1 = df_fair1 - ( mu - mu_1851_1900 )
df_fair2 = df_fair2 - ( mu - mu_1851_1900 )
df_fair3 = df_fair3 - ( mu - mu_1851_1900 )
df_fair4 = df_fair4 - ( mu - mu_1851_1900 )
df_fair5 = df_fair5 - ( mu - mu_1851_1900 )

#------------------------------------------------------------------------------
# MERGE: dataframes
# NB: PALEO + PAGES2k + HadCRUT5 + FaIR
#------------------------------------------------------------------------------

t = (np.floor(df_paleo.t_paleo).append(np.floor(df_pages2k.t_pages2k).append(np.floor(df_hadcrut5_yearly.t_hadcrut5))).append(np.floor(df_fair1.t_fair1)))
ts1 = (df_paleo.ts_paleo.append(df_pages2k.ts_pages2k.append(df_hadcrut5_yearly.ts_hadcrut5)).append(df_fair1.ts_fair1))
ts2 = (df_paleo.ts_paleo.append(df_pages2k.ts_pages2k.append(df_hadcrut5_yearly.ts_hadcrut5)).append(df_fair2.ts_fair2))
ts3 = (df_paleo.ts_paleo.append(df_pages2k.ts_pages2k.append(df_hadcrut5_yearly.ts_hadcrut5)).append(df_fair3.ts_fair3))
ts4 = (df_paleo.ts_paleo.append(df_pages2k.ts_pages2k.append(df_hadcrut5_yearly.ts_hadcrut5)).append(df_fair4.ts_fair4))
ts5 = (df_paleo.ts_paleo.append(df_pages2k.ts_pages2k.append(df_hadcrut5_yearly.ts_hadcrut5)).append(df_fair5.ts_fair5))

df = pd.DataFrame()
df['Year'] = t
df['Global1'] = ts1 - mu
df['Global2'] = ts2 - mu
df['Global3'] = ts3 - mu
df['Global4'] = ts4 - mu
df['Global5'] = ts5 - mu
df_sorted = df.sort_values(by=['Year']).dropna()
df = df_sorted.copy().reset_index(drop=True)

x = (df['Year']).values
y1 = np.array( df['Global1'] )
y2 = np.array( df['Global2'] )
y3 = np.array( df['Global3'] )
y4 = np.array( df['Global4'] )
y5 = np.array( df['Global5'] )
z1 = np.array(len(y1)*[ 1.0 ])
z2 = np.array(len(y2)*[ 1.0 ])
z3 = np.array(len(y3)*[ 1.0 ])
z4 = np.array(len(y4)*[ 1.0 ])
z5 = np.array(len(y5)*[ 1.0 ])

#------------------------------------------------------------------------------
# SMOOTH: rolling average ( nsmooth )  --> x, y, z
#------------------------------------------------------------------------------

if use_smoothing == True:
    x = x.astype(int)
    y1 = pd.Series( y1 ).rolling(nsmooth,center=True).mean().values
    y2 = pd.Series( y2 ).rolling(nsmooth,center=True).mean().values
    y3 = pd.Series( y3 ).rolling(nsmooth,center=True).mean().values
    y4 = pd.Series( y4 ).rolling(nsmooth,center=True).mean().values
    y5 = pd.Series( y5 ).rolling(nsmooth,center=True).mean().values
    mask = np.isfinite(y1)
    x = x[mask]
    y1 = y1[mask]
    y2 = y2[mask]
    y3 = y3[mask]
    y4 = y4[mask]
    y5 = y5[mask]
    z1 = z1[mask]    
    z2 = z2[mask]    
    z3 = z3[mask]    
    z4 = z4[mask]    
    z5 = z5[mask]    
    smoothstr = str(nsmooth).zfill(2) + 'yr-smooth' + '-'       
else:
    smoothstr = ''

#------------------------------------------------------------------------------
# RESCALE: colormap to max = cbar_max ( provide )
#------------------------------------------------------------------------------

y_norm1 = ( y1-y1.min() ) / ( y1.max() - y1.min() )
y_norm2 = ( y2-y2.min() ) / ( y2.max() - y2.min() )
y_norm3 = ( y3-y3.min() ) / ( y3.max() - y3.min() )
y_norm4 = ( y4-y4.min() ) / ( y4.max() - y4.min() )
y_norm5 = ( y5-y5.min() ) / ( y5.max() - y5.min() )

def rescale_colormap(cbar_max, y):
    
    colorscalefactor = cbar_max / df_hadcrut5.ts_hadcrut5.max()
    y_min = df_hadcrut5.ts_hadcrut5.min() * colorscalefactor
    y_max = df_hadcrut5.ts_hadcrut5.max() * colorscalefactor
    y_norm = (y - y_min) / (y_max - y_min) 
    y_norm_raw = ( y-y.min() ) / ( y.max() - y.min() )
    maxval = y_max
    minval = y_min
    colors = cmap( y_norm ) 
    norm = mcolors.TwoSlopeNorm( vmin=minval, vcenter=0.0, vmax=maxval) 
    sm = ScalarMappable( cmap=cmap, norm=norm )
    sm.set_array([])

    return colors, norm, sm

colors1, norm1, sm1 = rescale_colormap(cbar_max, y1)
colors2, norm2, sm2 = rescale_colormap(cbar_max, y2)
colors3, norm3, sm3 = rescale_colormap(cbar_max, y3)
colors4, norm4, sm4 = rescale_colormap(cbar_max, y4)
colors5, norm5, sm5 = rescale_colormap(cbar_max, y5)
    
#------------------------------------------------------------------------------
# PLOTS
#------------------------------------------------------------------------------

Y = y1[(x>=t_start)&(x<=2020)] # from t_start to 2020 CE 
Y1 = y1[ x>=2021 ] 		# SSP1-1.9
Y2 = y2[ x>=2021 ] 		# SSP1-2.6
Y3 = y3[ x>=2021 ] 		# SSP2-4.5
Y4 = y4[ x>=2021 ] 		# SSP3-7.0
Y5 = y5[ x>=2021 ] 		# SSP5-8.5

colors, norm, sm = rescale_colormap(cbar_max, Y)

Y_norm = ( Y-Y.min() ) / ( Y.max() - Y.min() )    
Y_norm1 = ( Y1-Y1.min() ) / ( Y1.max() - Y1.min() )    
Y_norm2 = ( Y2-Y2.min() ) / ( Y2.max() - Y2.min() )    
Y_norm3 = ( Y3-Y3.min() ) / ( Y3.max() - Y3.min() )    
Y_norm4 = ( Y4-Y4.min() ) / ( Y4.max() - Y4.min() )    
Y_norm5 = ( Y5-Y5.min() ) / ( Y5.max() - Y5.min() )    

figstr = 'climate-stripes' + '-' + 'projections' + '.png'
titlestr = 'Global mean anomaly since 500 CE: ' + 'SSP projections 2021-2200 CE'

ncolumns = int(np.floor((2200-t_start)/180))
if ncolumns < 2: ncolumns = 2
gs = gridspec.GridSpec(nrows=5, ncols=ncolumns, left=0, right=0.9, top=0.9, bottom=0.1, wspace=0.0, hspace=0.0)
fig = plt.figure(constrained_layout=True, figsize=(15,10) )
ax0 = plt.subplot(gs[:, 0:ncolumns-1]); ax0.axis('off') # row 0, col 0
ax1 = plt.subplot(gs[0, ncolumns-1]); ax1.axis('off') 	# row 0, col 1
ax2 = plt.subplot(gs[1, ncolumns-1]); ax2.axis('off') 	# row 1, col 1
ax3 = plt.subplot(gs[2, ncolumns-1]); ax3.axis('off') 	# row 2, col 1
ax4 = plt.subplot(gs[3, ncolumns-1]); ax4.axis('off') 	# row 3, col 1
ax5 = plt.subplot(gs[4, ncolumns-1]); ax5.axis('off') 	# row 4, col 1
ax0.bar( x[(x>=t_start)&(x<=2020)], z1[(x>=t_start)&(x<=2020)], color=colors, width=1.0 )
ax1.bar( x[x>=2021], z5[x>=2021], color=colors5[x>=2021], width=1.0 )
ax2.bar( x[x>=2021], z4[x>=2021], color=colors4[x>=2021], width=1.0 )
ax3.bar( x[x>=2021], z3[x>=2021], color=colors3[x>=2021], width=1.0 )
ax4.bar( x[x>=2021], z2[x>=2021], color=colors2[x>=2021], width=1.0 )
ax5.bar( x[x>=2021], z1[x>=2021], color=colors1[x>=2021], width=1.0 )
ax0.autoscale(tight=True)
ax1.autoscale(tight=True)
ax2.autoscale(tight=True)
ax3.autoscale(tight=True)
ax4.autoscale(tight=True)
ax5.autoscale(tight=True)
if use_overlay_axis == True: 
    ax0.axis('on')
    ax1.axis('off')
    ax2.axis('off')
    ax3.axis('off')
    ax4.axis('off')
    ax5.axis('on')
       
    def adjust_spines(ax, spines):
        for loc, spine in ax.spines.items():
            if loc in spines:
                spine.set_position(('outward', 5))
            else:
                spine.set_color('none')  
        if 'left' in spines:
            ax.yaxis.set_ticks_position('left')
        else:
            ax.yaxis.set_ticks([])
        if 'bottom' in spines:
            ax.xaxis.set_ticks_position('bottom')
        else:
            ax.xaxis.set_ticks([]) 
                        
    xlabels = map(str,np.arange(500,2100,500))
    plt.xticks( np.arange(500,2100,500), xlabels, rotation=0)
    plt.xlim([500,2020])        
    ax0.tick_params('both',length=7.5,width=2,which='major')             
    adjust_spines(ax0, ['bottom'])            
    ax0.spines['bottom'].set_linewidth(2)
    ax0.spines['bottom'].set_color('white')
    
    xlabels = map(str,np.arange(2050,2250,50))
    plt.xticks(np.arange(2050,2250,50),xlabels,rotation=0)
    plt.xlim([2020,2200])    
    ax5.tick_params('both',length=7.5,width=2,which='major')             
    adjust_spines(ax5, ['bottom'])            
    ax5.spines['bottom'].set_linewidth(2)
    ax5.spines['bottom'].set_color('white')
        
if use_overlay_colorbar == True:
    cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
    cbar = fig.colorbar(sm, cax=cbar_ax, shrink=0.5, extend='both')
    cbar.set_label( cbarstr, rotation=270, labelpad=25, fontsize=14 )
if use_overlay_timeseries == True: 
    ax0.plot( x[(x>=t_start)&(x<=2020)], 1*Y_norm/5., color='grey', ls='-', lw=1 )
    ax0.axvline(x=2021.0, ls='dashed', lw=2, color='black')
    ax0.text( x[x==t_start + int((2020-t_start)/2.0)], 0.5, str(t_start) + '-2020 CE', color='black', weight='bold', fontsize=20)
    ax1.plot( x[x>=2021], Y_norm5, color='black', ls='-', lw=1 )
    ax2.plot( x[x>=2021], Y_norm4*y4[-1]/y5[-1], color='black', ls='-', lw=1 )
    ax3.plot( x[x>=2021], Y_norm3*y3[-1]/y5[-1], color='black', ls='-', lw=1 )
    ax4.plot( x[x>=2021], Y_norm2*y2[-1]/y5[-1], color='black', ls='-', lw=1 )
    ax5.plot( x[x>=2021], Y_norm1*y1[-1]/y5[-1], color='black', ls='-', lw=1 )
    ax1.text(2030, 0.9, 'SSP5-8.5', color='black', weight='bold', fontsize=10)
    ax2.text(2030, 0.9, 'SSP3-7.0', color='black', weight='bold', fontsize=10)
    ax3.text(2030, 0.9, 'SSP2-4.5', color='black', weight='bold', fontsize=10)
    ax4.text(2030, 0.9, 'SSP1-2.6', color='black', weight='bold', fontsize=10)
    ax5.text(2030, 0.9, 'SSP1-1.9', color='black', weight='bold', fontsize=10)    
#plt.tight_layout()
plt.savefig( figstr, dpi=300 )
plt.close(fig)

#------------------------------------------------------------------------------
print('** END')

