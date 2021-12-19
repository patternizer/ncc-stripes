#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: plot_ncc_stripes_panel-3a.py
#------------------------------------------------------------------------------
# Version 0.1
# 2 December, 2021
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
cbar_max = 6.0
barwidthfraction = 1.0
t_start = 1400
t_end = 1850

use_timemask = True
use_logarithm = False
use_log10_scale = False
use_hires_norwich = False
use_data_cmax = False

use_dark_theme = True
use_smoothing = True
use_overlay_axis = True
use_overlay_timeseries = True
use_overlay_colorbar = True

plot_forecast_variability = True
plot_color_mapping = True
plot_climate_timeseries = False
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

titlestr = 'Global mean anomaly, 1400 CE - 1850 CE'
 
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

# PANEL 3a: 1400 - 1850 CE ----------------------------------------------------

p7unit = 1e0
p7start = 1400
p7end = 1850
p7n = int( (p7end - p7start) / p7unit )

y7 = df[ (df.Year >= 1400) & (df.Year < 1850) ]['Global']
x7 = df[ (df.Year >= 1400) & (df.Year < 1850) ]['Year']
s7 = stats.binned_statistic( x7, y7, 'mean', bins=np.linspace( p7start, p7end, p7n + 1) )
p7x = np.linspace( p7start + p7unit/2., p7end - p7unit/2., p7n ) 
p7y = s7.statistic

x = np.array( list(p7x) )
y = np.array( list(p7y) )
z = np.array(len(y)*[1.0])

print( 'p7n=', p7n )

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

if use_data_cmax == True:
    
    cbar_max = np.nanmax(y)

else:
    
    cbar_max = cbar_max

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

    figstr = 'climate-timeseries' + '-' + baselinestr + '-' + 'Panel-3a' + '.png'
    
    fig, ax = plt.subplots( figsize=(15,5) ); ax.axis('off')
    if use_log10_scale == True:                        
        plt.plot( np.array( [0.3,1] + list( np.arange( 2, len(x) ) ) ), y, ls='-', lw=0.5, color='grey', zorder=0 )
        plt.scatter( np.array( [0.3,1] + list( np.arange( 2, len(x) ) ) ), y, c=colors, s=10, cmap=cmap, norm=norm, zorder=1 )
        plt.plot( [0.3, len(x)], [0,0], ls='dashed', lw=0.5, color='white' )             
        ax.set_xscale("log")
    else:
#        plt.plot( np.arange( len(x) ), y, ls='-', lw=0.5, color='grey', zorder=0 )
#        plt.scatter( np.arange( len(x) ), y, c=colors, s=10, cmap=cmap, norm=norm, zorder=1 )
#        plt.plot( [0, len(x)], [0,0], ls='dashed', lw=0.5, color='white' )               
        plt.plot( x, y, ls='-', lw=0.5, color='grey', zorder=0 )
        plt.scatter( x, y, c=colors, s=10, cmap=cmap, norm=norm, zorder=1 )
        plt.plot( [x[0], x[-1]], [0,0], ls='dashed', lw=0.5, color='white' )               
    if use_overlay_axis == True: 
        ax.axis('on')
        ax=plt.gca(); ax.get_xaxis().set_visible(True)                
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

    figstr = 'climate-bars' + '-' + baselinestr + '-' + 'Panel-3a' + '.png'
    
    fig, ax = plt.subplots( figsize=(15,5) ); ax.axis('off')
    if use_log10_scale == True:                
        plt.bar( np.array( [0.3,1] + list(np.arange(2, len(x))) ) + 2, y, color=colors, width=barwidthfraction )       
        ax.set_xscale("log")
    else:
        plt.bar( np.arange( len(x) ), y, color=colors, width=barwidthfraction )       
    ax.axis('off')
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

    figstr = 'climate-stripes' + '-' + baselinestr + '-' + 'Panel-3a' + '.png'
        
    fig, ax = plt.subplots( figsize=(15,5) ); ax.axis('off')
    if use_log10_scale == True:                
        plt.bar( np.array( [0.3,1] + list(np.arange(2, len(x))) ) + 2, z, color=colors, width=barwidthfraction )       
        ax.set_xscale("log")
    else:
        plt.bar( np.arange( len(x) ), z, color=colors, width=barwidthfraction )   
    plt.ylim(0,1)        
    ax.axis('off')        
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

# PLOT (4): climate line ---------------------------------------------------
       
if plot_climate_stripes == True:

    figstr = 'climate-line' + '-' + baselinestr + '-' + 'Panel 3a' + '.png'
        
    fig, ax = plt.subplots( figsize=(15,5) ); ax.axis('off')
    plt.ylim(0,1)        
    ax.axis('off')        
    plt.plot( np.arange( len(x) ), y_norm_raw, color='white', ls='-', lw=1 )            
    fig.suptitle( titlestr, fontsize=fontsize )          
    plt.tick_params(labelsize=fontsize)    
    plt.tight_layout()
    plt.savefig( figstr, dpi=300 )
    plt.close(fig)
#------------------------------------------------------------------------------
print('** END')
