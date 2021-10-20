#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: plot_ncc_stripes_pleistocene.py
#------------------------------------------------------------------------------
# Version 0.1
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

fontsize = 10
nsmooth = 2 # years
cbar_max = 4.0
barwidthfraction = 1.0
t_start = -2580000
t_end = 2200

use_timemask = True
use_logarithm = False
use_log10_scale = False
use_hires_norwich = True

use_dark_theme = True
use_smoothing = True
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

#projectionstr = 'SSP119'
projectionstr = 'SSP126'
#projectionstr = 'SSP245'
#projectionstr = 'SSP370'
#projectionstr = 'SSP585'
 
baselinestr = 'baseline_1851_1900'
#baselinestr = 'baseline_1961_1990'
#baselinestr = 'baseline_1971_2000'

# titlestr = r'Global mean anomaly, 65 Myr ( $\le$ 2015) - 2200 CE: ' + projectionstr
titlestr = 'Global mean anomaly, 2.58 Myr ( < 2015) - 2200 CE: ' + projectionstr
 
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

t_pages2k = xr.cftime_range(start=years[0], periods=len(years), freq='A', calendar='gregorian')[0:1849]
ts_pages2k_instr = pd.to_numeric(obs[:,1][0:1849], errors='coerce')
ts_pages2k_recon = pd.to_numeric(obs[:,5][0:1849], errors='coerce')
ts_pages2k = np.append(ts_pages2k_recon[0:-36],ts_pages2k_instr[-36:],axis=None)
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
df_hadcrut5_yearly = df_hadcrut5_yearly[df_hadcrut5_yearly.t_hadcrut5 <= 2020]

#-----------------------------------------------------------------------------
# LOAD: FaIR v1.6.3 projections (constrained by HadCRUT5-analysis) --> df_fair
# NB: convert time to year.decimal
#-----------------------------------------------------------------------------

fair = pd.read_csv(fair_file)
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

# CONCATENATE: epochs and store in dataframe

t_paleo = np.array( list(t_epica) + list(t_lisiecki) + list(t_zachos) ).ravel().astype(int)
ts_paleo = np.array( list(ts_epica) + list(ts_lisiecki) + list(ts_zachos) ).ravel()
df_paleo = pd.DataFrame()
df_paleo['t_paleo'] = t_paleo
df_paleo['ts_paleo'] = ts_paleo

# TRIM: to +1 CE

df_paleo = df_paleo[ df_paleo.t_paleo < 1 ].dropna()

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

df_fair = df_fair - ( mu - mu_1851_1900 )

#------------------------------------------------------------------------------
# MERGE: dataframes
# NB: PALEO + PAGES2k + HadCRUT5 + FaIR
#------------------------------------------------------------------------------

t = (np.floor(df_paleo.t_paleo).append(np.floor(df_pages2k.t_pages2k).append(np.floor(df_hadcrut5_yearly.t_hadcrut5))).append(np.floor(df_fair.t_fair)))
ts = (df_paleo.ts_paleo.append(df_pages2k.ts_pages2k.append(df_hadcrut5_yearly.ts_hadcrut5)).append(df_fair.ts_fair))

df = pd.DataFrame()
df['Year'] = t
df['Global'] = ts - mu
df_sorted = df.sort_values(by=['Year']).dropna()
df = df_sorted.copy().reset_index(drop=True)

#------------------------------------------------------------------------------
# BINNED STATISTICS
#------------------------------------------------------------------------------

# PANEL 1: -66 Myr (K-Pg) to -2.58 Myr (Pleistocene) ---------------------------------

p1unit = 1e6
p1start = -66000000
p1end = -2000000
p1n = int( (p1end - p1start) / p1unit )

y1 = df[ df.Year < -2.0e6 ]['Global']
x1 = df[ df.Year < -2.0e6 ]['Year']
s1 = stats.binned_statistic( x1, y1, 'mean', bins=np.linspace( p1start, p1end, p1n + 1 ) ) 
p1x = np.linspace( p1start + p1unit/2., p1end - p1unit/2., p1n ) # [-65.5, -2.5] Myr --> 64 bars
p1y = s1.statistic

# PANEL 2: -2.58 Myr (Pleistocene) to -478 Kyr (Anglian Glacial) --------------

p2unit = 1e3
p2start = -2000000
p2end = -500000
p2n = int( (p2end - p2start) / p2unit )

y2 = df[ (df.Year >= -2.6e6) & (df.Year < -450.0e3) ]['Global']
x2 = df[ (df.Year >= -2.6e6) & (df.Year < -450.0e3) ]['Year']
s2 = stats.binned_statistic( x2, y2, 'mean', bins=np.linspace( p2start, p2end, p2n + 1 ) )
p2x = np.linspace( p2start + p2unit/2., p2end - p2unit/2., p2n ) # [-2.55, -0.45] Myr --> 22 bars
p2y = s2.statistic
 
# PANEL 3: -478 Kyr (Anglian Glacial) to 11700 (Holocene) ---------------------

p3unit = 1e3
p3start = -480000
p3end = -10000
p3n = int( (p3end - p3start) / p3unit )

y3 = df[ (df.Year >= -500.0e3) & (df.Year < -12000) ]['Global']
x3 = df[ (df.Year >= -500.0e3) & (df.Year < -12000) ]['Year']
s3 = stats.binned_statistic( x3, y3, 'mean', bins=np.linspace( p3start, p3end, p3n + 1 ) )
p3x = np.linspace( p3start + p3unit/2, p3end - p3unit/2., p3n ) # [-475, -15] Kyr --> 47 bars
p3y = s3.statistic

# PANEL 4: 11700 (Holocene) to 500 CE (Norwich) -------------------------------

p4unit = 1e2
p4start = -15000
#p4end = 1000
p4end = 50
p4n = int( (p4end - p4start) / p4unit )

#y4 = df[ (df.Year >= -20000) & (df.Year < 1000) ]['Global']
#x4 = df[ (df.Year >= -20000) & (df.Year < 1000)]['Year']
y4 = df[ (df.Year >= -20000) & (df.Year < 0) ]['Global']
x4 = df[ (df.Year >= -20000) & (df.Year < 0)]['Year']
s4 = stats.binned_statistic( x4, y4, 'mean', bins=np.linspace( p4start, p4end, p4n + 1) )
p4x = np.linspace( p4start + p4unit/2., p4end - p4unit/2., p4n ) # [-11500, +500] Kyr --> 13 bars
p4y = s4.statistic

# PANEL 5: 500 (Norwich) to 1850 CE (HadCRUT5) --------------------------------

if use_hires_norwich == True:
    p5unit = 1e1
else:
    p5unit = 1e2
#p5start = 500
p5start = 0
p5end = 1855
p5n = int( (p5end - p5start) / p5unit )

#y5 = df[ (df.Year >= 400) & (df.Year < 1900) ]['Global']
#x5 = df[ (df.Year >= 400) & (df.Year < 1900) ]['Year']
y5 = df[ (df.Year >= 0) & (df.Year < 1900) ]['Global']
x5 = df[ (df.Year >= 0) & (df.Year < 1900) ]['Year']
s5 = stats.binned_statistic( x5, y5, 'mean', bins=np.linspace( p5start, p5end, p5n + 1) )
p5x = np.linspace( p5start + p5unit/2., p5end - p5unit/2., p5n ) # [+550, 1850] --> 14 bars
p5y = s5.statistic

# PANEL 6: 1850 (HadCRUT5) to 2200 CE (FaIR) ----------------------------------

p6unit = 1e0 # decadal
p6start = 1850
p6end = 2200
p6n = int( (p6end - p6start) / p6unit )

y6 = df[ df.Year >= 1800 ]['Global']
x6 = df[ df.Year >= 1800 ]['Year']
s6 = stats.binned_statistic( x6, y6, 'mean', bins=np.linspace( p6start, p6end, p6n + 1) )
p6x = np.linspace( p6start + p6unit/2., p6end - p6unit/2., p6n ) # [+1855, 2195] --> 35 bars
p6y = s6.statistic

print('n1=', p1n)
print('n2=', p2n)
print('n3=', p3n)
print('n4=', p4n)
print('n5=', p5n)
print('n6=', p6n)

x = np.array( list(p1x) + list(p2x) + list(p3x) + list(p4x) + list(p5x) + list(p6x) )
y = np.array( list(p1y) + list(p2y) + list(p3y) + list(p4y) + list(p5y) + list(p6y) )

#x = df['Year']
#y = df['Global']
z = np.array(len(y)*[1.0])

#------------------------------------------------------------------------------
# SMOOTH: rolling average ( nsmooth )  --> x, y, z
#------------------------------------------------------------------------------

if use_smoothing == True:
    x = x.astype(int)
    y = pd.Series( y ).rolling(nsmooth,center=True).mean().values
    mask = np.isfinite(y)
    x = x[mask]
    y = y[mask]
    z = z[mask]    
    smoothstr = str(nsmooth).zfill(2) + 'yr-smooth' + '-'       
else:
    smoothstr = ''

#------------------------------------------------------------------------------
# ( EXPERIMENTAL ) log10 ( analytic continuation at ymin ) -----------------------
#------------------------------------------------------------------------------

if use_logarithm == True: 
    
    y = np.log10( y - np.nanmin(y) )        
    mask_inf = np.isinf(y)
    y[mask_inf] = np.nanmin(y[~mask_inf])
    yinterp = (pd.Series(y).rolling(nsmooth, center=True).mean()).values  # re-smooth to re-fill at singular value(s)
    y = yinterp    
    cbarstr = 'log₁₀ (Anomaly-min)'
    logstr = '-' + 'log10'
else:
    logstr = ''
    
#------------------------------------------------------------------------------
# TRIM: to [t_start, t_end]
#------------------------------------------------------------------------------

if use_timemask == True:
    
    timemask = (x>=t_start) & (x<=t_end)
    x = x[timemask]
    y = y[timemask]
    z = z[timemask]

#------------------------------------------------------------------------------
# RESCALE: colormap to max = cbar_max ( provide )
#------------------------------------------------------------------------------

y_norm_raw = ( y-np.nanmin(y) ) / ( np.nanmax(y) - np.nanmin(y) )

def rescale_colormap(cbar_max):
    
    colorscalefactor = cbar_max / df_hadcrut5.ts_hadcrut5.max()
    y_min = df_hadcrut5.ts_hadcrut5.min() * colorscalefactor
    y_max = df_hadcrut5.ts_hadcrut5.max() * colorscalefactor
    y_norm = (y - y_min) / (y_max - y_min) 
    maxval = y_max
    minval = y_min
    colors = cmap( y_norm ) 
    norm = mcolors.TwoSlopeNorm( vmin=minval, vcenter=0.0, vmax=maxval) 
    sm = ScalarMappable( cmap=cmap, norm=norm )
    sm.set_array([])

    return colors, norm, sm

colors, norm, sm = rescale_colormap(cbar_max)

#==============================================================================
# PLOTS
#==============================================================================
                            
# PLOT (1): climate timeseries ------------------------------------------------

if plot_climate_timeseries == True:

    figstr = 'climate-timeseries' + '-' + projectionstr + '-' + str(nsmooth).zfill(2) + 'yr-smooth' + logstr + '-' + baselinestr + '-' + 'pleistocene' + '.png'
    
    fig, ax = plt.subplots( figsize=(15,5) ); ax.axis('off')
    if use_log10_scale == True:                        
        # X = np.diff( np.log10( np.arange(len(x)) ))
        # Y = [ X[i+1]/X[i] for i in range(len(x)-2) ]
        plt.plot( np.array( [0.3,1] + list( np.arange( 2, len(x) ) ) ), y, ls='-', lw=0.5, color='grey', zorder=0 )
        plt.scatter( np.array( [0.3,1] + list( np.arange( 2, len(x) ) ) ), y, c=colors, s=10, cmap=cmap, norm=norm, zorder=1 )
        plt.plot( [0.3, len(x)], [0,0], ls='dashed', lw=0.5, color='white' )             
        ax.set_xscale("log")
    else:
        plt.plot( np.arange( len(x) ), y, ls='-', lw=0.5, color='grey', zorder=0 )
        plt.scatter( np.arange( len(x) ), y, c=colors, s=10, cmap=cmap, norm=norm, zorder=1 )
        plt.plot( [0, len(x)], [0,0], ls='dashed', lw=0.5, color='white' )               
    if use_overlay_axis == True: 
        ax.axis('on')
        ax=plt.gca(); ax.get_xaxis().set_visible(False)                
    if use_overlay_colorbar == True:    
        cbar = plt.colorbar( sm, shrink=0.5, extend='both' )
        cbar.set_label( cbarstr, rotation=270, labelpad=25, fontsize=fontsize )
    fig.suptitle( titlestr, fontsize=fontsize )       
    plt.tick_params(labelsize=fontsize)    
    plt.tight_layout()
    plt.savefig( figstr, dpi=300 )
    plt.close(fig)

# PLOT (2): climate bars ------------------------------------------------------

if plot_climate_bars == True:

    figstr = 'climate-bars' + '-' + projectionstr + '-' + str(nsmooth).zfill(2) + 'yr-smooth' + logstr + '-' + baselinestr + '-' + 'pleistocene' + '.png'
    
    fig, ax = plt.subplots( figsize=(15,5) ); ax.axis('off')
    if use_log10_scale == True:                
        plt.bar( np.array( [0.3,1] + list(np.arange(2, len(x))) ) + 2, y, color=colors, width=barwidthfraction )       
        ax.set_xscale("log")
    else:
        plt.bar( np.arange( len(x) ), y, color=colors, width=barwidthfraction )       
    if use_overlay_axis == True: 
        ax.axis('on')
        ax=plt.gca(); ax.get_xaxis().set_visible(False)                    
    if use_overlay_colorbar == True:    
        cbar = plt.colorbar( sm, shrink=0.5, extend='both' )
        cbar.set_label( cbarstr, rotation=270, labelpad=25, fontsize=fontsize ) 
    fig.suptitle( titlestr, fontsize=fontsize )          
    plt.tick_params(labelsize=fontsize)    
    plt.tight_layout()
    plt.savefig( figstr, dpi=300 )
    plt.close(fig)

# PLOT (3): climate stripes ---------------------------------------------------
       
if plot_climate_stripes == True:

    figstr = 'climate-stripes' + '-' + projectionstr + '-' + str(nsmooth).zfill(2) + 'yr-smooth' + logstr + '-' + baselinestr + '-' + 'pleistocene' + '.png'
        
    fig, ax = plt.subplots( figsize=(15,5) ); ax.axis('off')
    if use_log10_scale == True:                
        plt.bar( np.array( [0.3,1] + list(np.arange(2, len(x))) ) + 2, z, color=colors, width=barwidthfraction )       
        ax.set_xscale("log")
    else:
        plt.bar( np.arange( len(x) ), z, color=colors, width=barwidthfraction )   
    plt.ylim(0,1)        
    if use_overlay_axis == True: 
        ax.axis('on')
        ax=plt.gca(); ax.get_xaxis().set_visible(False)
    if use_overlay_timeseries == True: 
        if use_log10_scale == True:                
            plt.plot( np.array( [0.3/2,1] + list(np.arange(2, len(x))) ) + 2, y_norm_raw, color='black', ls='-', lw=1 )
        else:
            plt.plot( np.arange( len(x) ), y_norm_raw, color='black', ls='-', lw=1 )            
    if use_overlay_colorbar == True:
        cbar = plt.colorbar( sm, shrink=0.5, extend='both' )
        cbar.set_label( cbarstr, rotation=270, labelpad=25, fontsize=fontsize )
    fig.suptitle( titlestr, fontsize=fontsize )          
    plt.tick_params(labelsize=fontsize)    
    plt.tight_layout()
    plt.savefig( figstr, dpi=300 )
    plt.close(fig)

# CHECK: colour mapping ( norm ) ----------------------------------------------

if plot_color_mapping == True:
        
    figstr = 'color-mapping-variables' + '-' + projectionstr + '-' + str(nsmooth).zfill(2) + 'yr-smooth' + logstr + '-' + baselinestr + '-' + 'pleistocene' + '.png'
    titlestr = 'Colour mapping ( using IPCC AR6 divergent map )'
    
    fig, ax = plt.subplots( figsize=(15,5) )
    plt.plot( x, y_norm_raw, label='y_norm')
    plt.plot( x, y, label='y')
    plt.axhline(y=1, ls='dashed', lw=1)
    plt.axhline(y=0, ls='dashed', lw=1)
    plt.axvline(x=baseline_midpoint, ls='dashed', lw=1)
    plt.title( titlestr, fontsize=fontsize )        
    plt.legend(loc='lower left', ncol=2, markerscale=3, facecolor='black', framealpha=0.5, fontsize=fontsize)       
    plt.tick_params(labelsize=fontsize)    
    plt.tight_layout()
    plt.savefig( figstr, dpi=300 )
    plt.close(fig)

# CHECK: internal variability traces ------------------------------------------

if plot_forecast_variability == True:
        
    figstr = 'internal-variability.png'
    titlestr = 'Internal variability model ( low: 1350-1549, high: 1550-1749 )'
    
    fig, ax = plt.subplots( figsize=(15,5) )
    plt.plot( df_variability.lo_var.values, label='Low variability ( 1350-1549 )')
    plt.plot( df_variability.hi_var.values, label='High variability ( 1550-1749 )')
    plt.legend()
    plt.title( titlestr, fontsize=fontsize )        
    plt.legend(loc='lower left', ncol=2, markerscale=3, facecolor='black', framealpha=0.5, fontsize=fontsize)       
    plt.tick_params(labelsize=fontsize)    
    plt.tight_layout()
    plt.savefig( figstr, dpi=300 )
    plt.close(fig)

#------------------------------------------------------------------------------
print('** END')

#colors = ['crimson', 'dodgerblue', 'teal', 'limegreen', 'gold']

# ADD: source data annotations

#    datastr1 = r'$\bf{Climate\ model\ forecast}$' + '( 2021-2200 AD) : FaIR v1.6.4 / https://fair.readthedocs.io'        
#    datastr2 = r'$\bf{Instrumental\ Record}$' + '( 1850-2020 AD) : HadCRUT5.0.1 / https://www.metoffice.gov.uk/hadobs/hadcrut5'        
#    datastr3 = r'$\bf{Proxy\ Record}$' + '( 1-2000 AD) : PAGES2k / doi:10.6084/m9.figshare.8143094.v3'        
#    datastr4 = r'$\bf{Geological\ Record\ 1}$' + ': EPICA Dome C Antarctica (x0.5) / doi:10.1594/PANGAEA.683655'
#    datastr5 = r'$\bf{Geological\ Record\ 2}$' + ': Lisiecki & Raymo, 2005a,b / doi:10.1029/2005PA001164'   
#    datastr6 = r'$\bf{Geological\ Record\ 3}$' + ': Zachos et al, 2008 / doi:10.1038/nature06588'        
#    authorstr = r'$\bf{Graphic}$' + ': Michael Taylor, CRU / UEA' + ' -- ' + titletime
#    fig.suptitle(titlestr, fontsize=fontsize, color='white', fontweight='bold')        
#    plt.annotate(datastr1, xy=(500,350), xycoords='figure pixels', color='white', fontsize=fontsize) 
#    plt.annotate(datastr2, xy=(500,300), xycoords='figure pixels', color='white', fontsize=fontsize) 
#    plt.annotate(datastr3, xy=(500,250), xycoords='figure pixels', color='white', fontsize=fontsize) 
#    plt.annotate(datastr4, xy=(500,200), xycoords='figure pixels', color='white', fontsize=fontsize) 
#    plt.annotate(datastr5, xy=(500,150), xycoords='figure pixels', color='white', fontsize=fontsize) 
#    plt.annotate(datastr6, xy=(500,100), xycoords='figure pixels', color='white', fontsize=fontsize) 
#    plt.annotate(authorstr, xy=(500,30), xycoords='figure pixels', color='white', fontsize=fontsize, bbox=dict(boxstyle="square, pad=0.3", fc='black', edgecolor='white', linewidth=0.2))    
#    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)    

