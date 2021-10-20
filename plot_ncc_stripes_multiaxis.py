#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: plot_ncc_stripes_multiaxis.py
#------------------------------------------------------------------------------
# Version 0.2
# 11 October, 2021
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

# Statistics libraries:
from scipy import stats

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

fontsize = 14
nsmooth = 2 # years
cbar_max = 6.0
barwidthfraction = 1.0

use_hires_norwich = True

use_dark_theme = True
use_smoothing = False
use_overlay_axis = True
use_overlay_timeseries = True
use_overlay_colorbar = True

plot_forecast_variability = True
plot_color_mapping = True
plot_climate_timeseries = True
plot_climate_bars = True
plot_climate_stripes = True
 
#projectionstr = 'RCP3pd'
#projectionstr = 'RCP45'
#projectionstr = 'RCP6'
#projectionstr = 'RCP85'

projectionstr = 'SSP119'
#projectionstr = 'SSP126'
#projectionstr = 'SSP245'
#projectionstr = 'SSP370'
#projectionstr = 'SSP585'
 
#baselinestr = 'baseline_1851_1900'
baselinestr = 'baseline_1961_1990'
#baselinestr = 'baseline_1971_2000'

titlestr = 'Global mean anomaly, 65 Myr ( < 2015) - 2200 CE: ' + projectionstr
#titlestr = 'Global mean anomaly, 2.58 Myr ( < 2015) - 2200 CE: ' + projectionstr
 
pathstr = 'DATA/'
pages2kstr = 'PAGES2k.txt'
hadcrut5str = 'HadCRUT5.csv'
fairstr = 'fair' + '_' + projectionstr.lower() + '.csv' 
lovarstr = 'variability_realisation0.txt' 
hivarstr = 'variability_realisation1.txt' 
fairstr = 'fair' + '_' + projectionstr.lower() + '.csv' 
paleostr = 'paleo_data_compilation.xls'
pages2k_file = pathstr + pages2kstr 
hadcrut5_file = pathstr + hadcrut5str 
fair_file = pathstr + fairstr 
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
    
    matplotlib.rcParams['text.usetex'] = False
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

    matplotlib.rcParams['text.usetex'] = False
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
# METHODS
#------------------------------------------------------------------------------    

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
        
#-----------------------------------------------------------------------------
# LOAD: PAGES2k (via Ed Hawkins with thanks) --> df_pages2k
# NB: convert time to year.decimal
#-----------------------------------------------------------------------------

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

t_pages2k = xr.cftime_range(start=years[0], periods=len(years), freq='A', calendar='gregorian')[0:1850]
ts_pages2k_instr = pd.to_numeric(obs[:,1][0:1850], errors='coerce')
ts_pages2k_recon = pd.to_numeric(obs[:,5][0:1850], errors='coerce')
#ts_pages2k = np.append(ts_pages2k_recon[0:-36],ts_pages2k_instr[-36:],axis=None)
ts_pages2k = ts_pages2k_recon
df_pages2k = pd.DataFrame()
df_pages2k['t_pages2k'] = t_pages2k.year.astype(float)
df_pages2k['ts_pages2k'] = ts_pages2k

#-----------------------------------------------------------------------------
# LOAD: HadCRUT5 (via Tim Osborn and UKMO with thanks) --> df_hadcrut5
# NB: convert time to year.decimal
#-----------------------------------------------------------------------------

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
df_hadcrut5_yearly = df_hadcrut5_yearly[ (df_hadcrut5_yearly.t_hadcrut5 >= 1851) & (df_hadcrut5_yearly.t_hadcrut5 <= 2020) ]

#-----------------------------------------------------------------------------
# LOAD: FaIR v1.6.3 projections (constrained by HadCRUT5-analysis) --> df_fair
# NB: convert time to year.decimal
#-----------------------------------------------------------------------------

fair = pd.read_csv(fair_file,index_col=0)
df_fair = pd.DataFrame()
df_fair['t_fair'] = fair.Year.values.astype('float')

#-----------------------------------------------------------------------------
# LOAD: internal variability for FaIR v1.6.4 projections calculated by Tim from the instrumental record.
#-----------------------------------------------------------------------------

nheader = 0
f_lo = open(lo_var_file)
f_hi = open(hi_var_file)
lines_lo = f_lo.readlines()
lines_hi = f_hi.readlines()
years = []
obs_lo = []
obs_hi = []
#for i in range(nheader,len(lineslo)):
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

if (projectionstr == 'SSP119') | (projectionstr == 'SSP126') | (projectionstr == 'SSP245') | (projectionstr == 'RCP3pd') | (projectionstr == 'RCP45'):
    df_fair['ts_fair'] = fair.Global.values + df_variability.lo_var.values
elif (projectionstr == 'SSP370') | (projectionstr == 'SSP585') | (projectionstr == 'RCP6') | (projectionstr == 'RCP85'):
    df_fair['ts_fair'] = fair.Global.values + df_variability.hi_var.values
        
#-----------------------------------------------------------------------------
# LOAD: geological anomalies: 65.5229 Myr ( before 2015 )
# NB: paleo_file has copy-paste of "data_compilation" sheet values from All_palaeotemps.xlsx  
#     
#-----------------------------------------------------------------------------
 
#import xlrd
#workbook = xlrd.open_workbook(paleo_file)
#worksheet = workbook.sheet_by_index(0) # first sheet in workbook
#ncols = worksheet.utter_max_cols
#nrows = worksheet.utter_max_rows

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

df_epica = pd.DataFrame()
df_epica['t'] = t_epica
df_epica['ts'] = ts_epica
df_epica_sorted = df_epica.sort_values(by=['t']).dropna().reset_index(drop=True)
df_epica = df_epica_sorted[ (df_epica_sorted.t <= 0) ]

df_lisiecki = pd.DataFrame()
df_lisiecki['t'] = t_lisiecki
df_lisiecki['ts'] = ts_lisiecki
df_lisiecki_sorted = df_lisiecki.sort_values(by=['t']).dropna().reset_index(drop=True)
df_lisiecki = df_lisiecki_sorted[ (df_lisiecki_sorted.t <= 0) ]

df_zachos = pd.DataFrame()
df_zachos['t'] = t_zachos
df_zachos['ts'] = ts_zachos
df_zachos_sorted = df_zachos.sort_values(by=['t']).dropna().reset_index(drop=True)
df_zachos = df_zachos_sorted[ (df_zachos_sorted.t <= 0) ]

# CONCATENATE: epochs and store in dataframe

t_paleo = np.array( list(df_epica.t) + list(df_lisiecki.t) + list(df_zachos.t) ).ravel().astype(int)
ts_paleo = np.array( list(df_epica.ts) + list(df_lisiecki.ts) + list(df_zachos.ts) ).ravel()
df_paleo = pd.DataFrame()
df_paleo['t_paleo'] = t_paleo
df_paleo['ts_paleo'] = ts_paleo
df_paleo_sorted = df_paleo.sort_values(by=['t_paleo']).dropna().reset_index(drop=True)
df_paleo = df_paleo_sorted[ (df_paleo_sorted.t_paleo <= 0) ]

#------------------------------------------------------------------------------
# COMPUTE: baseline yearly means from HadCRUT5 ( monthly )
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
cbarstr = r'Anomaly, $^{\circ}$C ( from ' + str(baseline_start) + '-' + str(baseline_end) +' )'

#------------------------------------------------------------------------------
# ALIGN: FaIR projections relative to chosen baseline
# NB: FaIR projections are calculated relative to 1880-2015 emission trends and so we
# need to subtract off the difference between the selected baseline and the 1851-1900 baseline
#------------------------------------------------------------------------------

df_fair['ts_fair'] = df_fair['ts_fair'] - ( mu - mu_1851_1900 )

#------------------------------------------------------------------------------
# MERGE: dataframes
# NB: PALEO + PAGES2k + HadCRUT5 + FaIR
#------------------------------------------------------------------------------

t = (np.floor(df_paleo.t_paleo).append(np.floor(df_pages2k.t_pages2k).append(np.floor(df_hadcrut5_yearly.t_hadcrut5))).append(np.floor(df_fair.t_fair)))
ts = (df_paleo.ts_paleo.append(df_pages2k.ts_pages2k.append(df_hadcrut5_yearly.ts_hadcrut5)).append(df_fair.ts_fair))

df = pd.DataFrame()
df['Year'] = t
df['Global'] = ts - mu
df_sorted = df.sort_values(by=['Year']).dropna().reset_index(drop=True)
df = df_sorted

df.to_csv('ncc_data' + '_' + baselinestr + '_' + projectionstr + '.csv')

x_raw = df.Year
y_raw = df.Global

#------------------------------------------------------------------------------
# BINNED STATISTICS
#------------------------------------------------------------------------------

# PANEL 1: -66 Myr (K-Pg) to -2.58 Myr (Pleistocene) ---------------------------------

p1unit = 1e6
p1start = -66e6
p1end = -3e6
p1n = int( (p1end - p1start) / p1unit )

y1 = df[ df.Year < -2.0e6 ]['Global']
x1 = df[ df.Year < -2.0e6 ]['Year']
s1 = stats.binned_statistic( x1, y1, 'mean', bins=np.linspace( p1start, p1end, p1n + 1 ) ) 
p1x = np.linspace( p1start + p1unit/2., p1end - p1unit/2., p1n ) # [-65.5, -2.5] Myr --> 64 bars
p1y = s1.statistic

# PANEL 2: -2.58 Myr (Pleistocene) to -478 Kyr (Anglian Glacial) --------------

p2unit = 1e3
p2start = -2.6e6
p2end = -480e3
p2n = int( (p2end - p2start) / p2unit )

y2 = df[ (df.Year >= -2.6e6) & (df.Year < -480.0e3) ]['Global']
x2 = df[ (df.Year >= -2.6e6) & (df.Year < -480.0e3) ]['Year']
s2 = stats.binned_statistic( x2, y2, 'mean', bins=np.linspace( p2start, p2end, p2n + 1 ) )
p2x = np.linspace( p2start + p2unit/2., p2end - p2unit/2., p2n ) # [-2.55, -0.45] Myr --> 22 bars
p2y = s2.statistic
 
# PANEL 3: -478 Kyr (Anglian Glacial) to 11700 (Holocene) ---------------------

p3unit = 1e3
p3start = -480000
p3end = -12000
p3n = int( (p3end - p3start) / p3unit )

y3 = df[ (df.Year >= -500.0e3) & (df.Year < -12000) ]['Global']
x3 = df[ (df.Year >= -500.0e3) & (df.Year < -12000) ]['Year']
s3 = stats.binned_statistic( x3, y3, 'mean', bins=np.linspace( p3start, p3end, p3n + 1 ) )
p3x = np.linspace( p3start + p3unit/2, p3end - p3unit/2., p3n ) # [-475, -15] Kyr --> 47 bars
p3y = s3.statistic

# PANEL 4: 11700 (Holocene) to 500 CE (Norwich) -------------------------------

p4unit = 1e3
p4start = -12e3
p4end = 0
p4n = int( (p4end - p4start) / p4unit )

y4 = df[ (df.Year >= -11700) & (df.Year <= 0) ]['Global']
x4 = df[ (df.Year >= -11700) & (df.Year <= 0)]['Year']
s4 = stats.binned_statistic( x4, y4, 'mean', bins=np.linspace( p4start, p4end, p4n + 1) )
p4x = np.linspace( p4start + p4unit/2., p4end - p4unit/2., p4n ) # [-11500, +500] Kyr --> 13 bars
p4y = s4.statistic

# PANEL 5: 1 to 1850 CE ( PAGES2k ) -------------------------------------------

p5unit = 1e0
#p5start = 500
p5start = 10
p5end = 1850
p5n = int( (p5end - p5start) / p5unit )
y5 = df[ (df.Year >= 1) & (df.Year <= 1850) ]['Global']
x5 = df[ (df.Year >= 1) & (df.Year <= 1850) ]['Year']
s5 = stats.binned_statistic( x5, y5, 'mean', bins=np.linspace( p5start, p5end, p5n + 1) )
p5x = np.linspace( p5start + p5unit/2., p5end - p5unit/2., p5n ) # [1, 1850] --> 1850 bars
p5y = s5.statistic

# PANEL 6: 1850 to 2020 (HadCRUT5) --------------------------------------------

p6unit = 1e0 # decadal
p6start = 1851
p6end = 2020
p6n = int( (p6end - p6start) / p6unit )

y6 = df[ (df.Year >= 1851) & (df.Year <= 2020) ]['Global']
x6 = df[ (df.Year >= 1851) & (df.Year <= 2020) ]['Year']
s6 = stats.binned_statistic( x6, y6, 'mean', bins=np.linspace( p6start, p6end, p6n + 1) )
p6x = np.linspace( p6start + p6unit/2., p6end - p6unit/2., p6n ) # [+1851, 2020] --> 170 bars
p6y = s6.statistic

# PANEL 7: 2021 to 2200 CE (FaIR) ---------------------------------------------

p7unit = 1e0 # decadal
p7start = 2021
p7end = 2200
p7n = int( (p7end - p7start) / p7unit )

y7 = df[ df.Year >= 2021 ]['Global']
x7 = df[ df.Year >= 2021 ]['Year']
s7 = stats.binned_statistic( x7, y7, 'mean', bins=np.linspace( p7start, p7end, p7n + 1) )
p7x = np.linspace( p7start + p7unit/2., p7end - p7unit/2., p7n ) # [+2021, 2200] --> 180 bars
p7y = s7.statistic

#mask1 = np.isfinite(p1y); p1x=p1x[mask1]; p1y=p1y[mask1]; p1n=mask1.sum()
#mask2 = np.isfinite(p2y); p2x=p2x[mask2]; p2y=p2y[mask2]; p2n=mask2.sum()
#mask3 = np.isfinite(p3y); p3x=p3x[mask3]; p3y=p3y[mask3]; p3n=mask3.sum()
#mask4 = np.isfinite(p4y); p4x=p4x[mask4]; p4y=p4y[mask4]; p4n=mask4.sum()
#mask5 = np.isfinite(p5y); p5x=p5x[mask5]; p5y=p5y[mask5]; p5n=mask5.sum()
#mask6 = np.isfinite(p6y); p6x=p6x[mask6]; p6y=p6y[mask6]; p6n=mask6.sum()
#mask7 = np.isfinite(p7y); p7x=p7x[mask7]; p7y=p6y[mask7]; p7n=mask7.sum()

print('n1=', p1n)
print('n2=', p2n)
print('n3=', p3n)
print('n4=', p4n)
print('n5=', p5n)
print('n6=', p6n)
print('n7=', p7n)

x = np.array( list(p1x) + list(p2x) + list(p3x) + list(p4x) + list(p5x) + list(p6x) + list(p7x) )
y = np.array( list(p1y) + list(p2y) + list(p3y) + list(p4y) + list(p5y) + list(p6y) + list(p7y) )
#x = df['Year']
#y = df['Global']
z = np.array(len(y)*[1.0])

dg = pd.DataFrame()
dg['Year'] = x
dg['Global'] = y
dg.to_csv('ncc_data_binned' + '_' + baselinestr + '_' + projectionstr + '.csv')

#plt.plot( x_raw, y_raw)
#plt.plot( x, y)
#plt.xlim( x_raw[0], -2.58e6)
#plt.xlim( -2.58e6,-12000)
#plt.xlim( -12000, 0)
#plt.xlim( 1, 1850)
#plt.xlim( 1851, 2020)
#plt.xlim( 2021, 2200)
#plt.ylim(-1,8)

#------------------------------------------------------------------------------
# RESCALE: colormap to max = cbar_max ( provide )
#------------------------------------------------------------------------------

def rescale_colormap(cbar_max):
    
    colorscalefactor = cbar_max / df_hadcrut5.ts_hadcrut5.max()
#    y_min = df_hadcrut5.ts_hadcrut5.min() * colorscalefactor
#    y_max = df_hadcrut5.ts_hadcrut5.max() * colorscalefactor
    y_min = -cbar_max
    y_max = cbar_max
    y_norm = (y - y_min) / (y_max - y_min) 
    colors = cmap( y_norm ) 
    norm = mcolors.TwoSlopeNorm( vmin=y_min, vcenter=0.0, vmax=y_max) 
    sm = ScalarMappable( cmap=cmap, norm=norm )
    sm.set_array([])

    return colors, norm, sm

colors, norm, sm = rescale_colormap(cbar_max)

#==============================================================================
# PLOTS
#==============================================================================
                           
if plot_climate_bars == True:
 
    # NEW: climate line! 
    # plt.scatter(x,z,s=2,marker='s',color=colors)
    
    figstr = 'climate-bars' + '-' + projectionstr + '-' + str(nsmooth).zfill(2) + 'yr-smooth' + baselinestr + '-' + 'full_multiaxis' + '.png'
 
    gs = gridspec.GridSpec(nrows=1, ncols=7, left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.1, hspace=0.1)
    fig = plt.figure(constrained_layout=True, figsize=(15,10) )
        
    ax0 = fig.add_subplot(gs[0, 0]); ax0.axis('on') # row 0, col 0
    ax1 = fig.add_subplot(gs[0, 1]); ax1.axis('on') # row 0, col 1
    ax2 = fig.add_subplot(gs[0, 2]); ax2.axis('on') # row 0, col 2
    ax3 = fig.add_subplot(gs[0, 3]); ax3.axis('on') # row 0, col 3
    ax4 = fig.add_subplot(gs[0, 4]); ax4.axis('on') # row 0, col 4
    ax5 = fig.add_subplot(gs[0, 5]); ax5.axis('on') # row 0, col 5
    ax6 = fig.add_subplot(gs[0, 6]); ax6.axis('on') # row 0, col 6
    
    ax0.bar( x[0:p1n]/1e6, y[0:p1n], color=colors[0:p1n], width=barwidthfraction )
    ax1.bar( x[p1n:p1n+p2n]/1e3, y[p1n:p1n+p2n], color=colors[p1n:p1n+p2n], width=barwidthfraction )
    ax2.bar( x[p1n+p2n:p1n+p2n+p3n]/1e3, y[p1n+p2n:p1n+p2n+p3n], color=colors[p1n+p2n:p1n+p2n+p3n], width=barwidthfraction )
    ax3.bar( x[p1n+p2n+p3n:p1n+p2n+p3n+p4n]/1e3, y[p1n+p2n+p3n:p1n+p2n+p3n+p4n], color=colors[p1n+p2n+p3n:p1n+p2n+p3n+p4n], width=barwidthfraction )
    ax4.bar( x[p1n+p2n+p3n+p4n:p1n+p2n+p3n+p4n+p5n], y[p1n+p2n+p3n+p4n:p1n+p2n+p3n+p4n+p5n], color=colors[p1n+p2n+p3n+p4n:p1n+p2n+p3n+p4n+p5n], width=barwidthfraction )
    ax5.bar( x[p1n+p2n+p3n+p4n+p5n:p1n+p2n+p3n+p4n+p5n+p6n], y[p1n+p2n+p3n+p4n+p5n:p1n+p2n+p3n+p4n+p5n+p6n], color=colors[p1n+p2n+p3n+p4n+p5n:p1n+p2n+p3n+p4n+p5n+p6n], width=barwidthfraction )
    ax6.bar( x[p1n+p2n+p3n+p4n+p5n+p6n:], y[p1n+p2n+p3n+p4n+p5n+p6n:], color=colors[p1n+p2n+p3n+p4n+p5n+p6n:], width=barwidthfraction )

    ax0.set_xlabel('Millions of years BCE')
    ax1.set_xlabel('<---------------------')
    ax2.set_xlabel('Thousands of years BCE')
    ax3.set_xlabel('--------------------->')
    ax4.set_xlabel('<---------------------')
    ax5.set_xlabel('Year CE')
    ax6.set_xlabel('--------------------->')

    ax0.autoscale(tight=True)
    ax1.autoscale(tight=True)
    ax2.autoscale(tight=True)
    ax2.autoscale(tight=True)
    ax3.autoscale(tight=True)
    ax4.autoscale(tight=True)
    ax5.autoscale(tight=True)
    ax6.autoscale(tight=True)

    ax0.set_ylim([-5,15])
    ax1.set_ylim([-5,15])
    ax2.set_ylim([-5,15])
    ax3.set_ylim([-5,15])
    ax4.set_ylim([-5,15])
    ax5.set_ylim([-5,15])
    ax6.set_ylim([-5,15])
    
    ax1.yaxis.set_ticks([]) 
    ax2.yaxis.set_ticks([]) 
    ax3.yaxis.set_ticks([]) 
    ax4.yaxis.set_ticks([]) 
    ax5.yaxis.set_ticks([]) 
    ax6.yaxis.set_ticks([]) 

    ax0.xaxis.set_ticks( np.linspace(-60,-10,6) ); ax0.xaxis.set_ticklabels( np.linspace(-60,-10,6).astype(int) )
    ax1.xaxis.set_ticks( np.linspace(-2000,-1000,2) ); ax1.xaxis.set_ticklabels( np.linspace(-2000,-1000,2).astype(int) )
    ax2.xaxis.set_ticks( np.linspace(-500,-100,3) ); ax2.xaxis.set_ticklabels( np.linspace(-500,-100,3).astype(int) )
    ax3.xaxis.set_ticks( np.linspace(-10,-2,5) ); ax3.xaxis.set_ticklabels( np.linspace(-10,-2,5).astype(int) )
    ax4.xaxis.set_ticks( np.linspace(0,1500,4) ); ax4.xaxis.set_ticklabels( np.linspace(0,1500,4).astype(int) )
    ax5.xaxis.set_ticks( np.linspace(1850,1950,2) ); ax5.xaxis.set_ticklabels( np.linspace(1850,1950,2).astype(int) )
    ax6.xaxis.set_ticks( np.linspace(2000,2200,3) ); ax6.xaxis.set_ticklabels( np.linspace(2000,2200,3).astype(int) )

    ax0.spines[['right','top']].set_color('none')
    ax1.spines[['right','top']].set_color('none')
    ax2.spines[['right','top']].set_color('none')
    ax3.spines[['right','top']].set_color('none')
    ax4.spines[['right','top']].set_color('none')
    ax5.spines[['right','top']].set_color('none')
    ax6.spines[['right','top']].set_color('none')
  
    ax0.tick_params('both',length=7.5,width=2,which='major',color='lightgrey')             
    adjust_spines(ax0, ['bottom'])            
    ax0.spines['bottom'].set_linewidth(2)
    ax0.spines['bottom'].set_color('lightgrey')
    ax1.tick_params('both',length=7.5,width=2,which='major',color='grey')             
    adjust_spines(ax1, ['bottom'])            
    ax1.spines['bottom'].set_linewidth(2)
    ax1.spines['bottom'].set_color('grey')
    ax2.tick_params('both',length=7.5,width=2,which='major',color='lightgrey')             
    adjust_spines(ax2, ['bottom'])            
    ax2.spines['bottom'].set_linewidth(2)
    ax2.spines['bottom'].set_color('lightgrey')
    ax3.tick_params('both',length=7.5,width=2,which='major',color='grey')             
    adjust_spines(ax3, ['bottom'])            
    ax3.spines['bottom'].set_linewidth(2)
    ax3.spines['bottom'].set_color('grey')
    ax4.tick_params('both',length=7.5,width=2,which='major',color='lightgrey')             
    adjust_spines(ax4, ['bottom'])            
    ax4.spines['bottom'].set_linewidth(2)
    ax4.spines['bottom'].set_color('lightgrey')
    ax5.tick_params('both',length=7.5,width=2,which='major',color='grey')             
    adjust_spines(ax5, ['bottom'])            
    ax5.spines['bottom'].set_linewidth(2)
    ax5.spines['bottom'].set_color('grey')
    ax6.tick_params('both',length=7.5,width=2,which='major',color='lightgrey')             
    adjust_spines(ax6, ['bottom'])            
    ax6.spines['bottom'].set_linewidth(2)
    ax6.spines['bottom'].set_color('lightgrey')
                        
    if use_overlay_colorbar == True:
        cbar_ax = fig.add_axes([0.92, 0.25, 0.015, 0.5])
        cbar = fig.colorbar(sm, cax=cbar_ax, shrink=0.5, extend='both')
        cbar.set_label( cbarstr, rotation=270, labelpad=25, fontsize=fontsize )    
        
    fig.suptitle( titlestr, fontsize=fontsize )          
    plt.tick_params(labelsize=fontsize)    
    plt.tight_layout()
    plt.savefig( figstr, dpi=300 )
    plt.close(fig)

#------------------------------------------------------------------------------
print('** END')
