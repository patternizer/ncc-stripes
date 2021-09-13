#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: plot-glosat-stripes.py
#------------------------------------------------------------------------------
# Version 0.4
# 14 September, 2021
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
from matplotlib import cm
from matplotlib.cm import ScalarMappable
from matplotlib import rcParams
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import cmocean as cmo

# Silence library version notifications
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
#------------------------------------------------------------------------------

#----------------------------------------------------------------------------
# DARK BACKGROUND THEME
#----------------------------------------------------------------------------
matplotlib.rcParams['text.usetex'] = False
rcParams['font.family'] = ['DejaVu Sans']
rcParams['font.sans-serif'] = ['Avant Garde']
plt.rc('text',color='white')
plt.rc('lines',color='white')
plt.rc('patch',edgecolor='white')
plt.rc('grid',color='lightgray')
plt.rc('xtick',color='white')
plt.rc('ytick',color='white')
plt.rc('axes',edgecolor='lightgray')
plt.rc('axes',facecolor='black')
plt.rc('axes',labelcolor='white')
plt.rc('figure',facecolor='black')
plt.rc('figure',edgecolor='black')
plt.rc('savefig',edgecolor='black')
plt.rc('savefig',facecolor='black')

# Calculate current time

now = datetime.now()
currentmn = str(now.month)
if now.day == 1:
    currentdy = str(cal.monthrange(now.year,now.month-1)[1])
    currentmn = str(now.month-1)
else:
    currentdy = str(now.day-1)
if int(currentdy) < 10:
    currentdy = '0' + currentdy    
currentyr = str(now.year)
if int(currentmn) < 10:
    currentmn = '0' + currentmn
titletime = str(currentdy) + '/' + currentmn + '/' + currentyr

#------------------------------------------------------------------------------
# SETTINGS: 
#------------------------------------------------------------------------------

nsmooth = 10
use_smoothing = True
use_logarithm = False
use_only_hadcrut5 = False

use_overlay_timeseries = True
use_overlay_baseline = True
use_overlay_colorbar = True
use_overlay_axis = False
#use_overlay_annotations = True

projectionstr = 'RCP3pd'
#projectionstr = 'RCP45'
#projectionstr = 'RCP6'
#projectionstr = 'RCP85'

#baseline = 'baseline_1851_1900'
baseline = 'baseline_1961_1990'
#baseline = 'baseline_1971_2000'
 
fontsize = 20
pathstr = 'DATA/'
pages2kstr = 'PAGES2k.txt'
hadcrut5str = 'HadCRUT5.csv'
fairstr = 'fair' + '_' + projectionstr.lower() + '.csv' 
pages2k_file = pathstr + pages2kstr 
hadcrut5_file = pathstr + hadcrut5str 
fair_file = pathstr + fairstr 
        
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

t_pages2k = xr.cftime_range(start=years[0], periods=len(years), freq='A', calendar='gregorian')[499:1849]
ts_pages2k_instr = pd.to_numeric(obs[:,1][499:1849], errors='coerce')
ts_pages2k_recon = pd.to_numeric(obs[:,5][499:1849], errors='coerce')
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
for yyyy in years:
    year_data = df_hadcrut5[np.floor(df_hadcrut5['t_hadcrut5']).astype('int') == yyyy]['ts_hadcrut5']
    yearly_mean = np.nanmean(year_data)
    yearly.append(yearly_mean)
df_hadcrut5_yearly = pd.DataFrame()
df_hadcrut5_yearly['t_hadcrut5'] = years.astype('float')
df_hadcrut5_yearly['ts_hadcrut5'] = yearly
df_hadcrut5_yearly = df_hadcrut5_yearly[df_hadcrut5_yearly.t_hadcrut5 <= 2020]

#-----------------------------------------------------------------------------
# LOAD: FaIR v1.6.3 projections (constrained by HadCRUT5-analysis) --> df_fair
# NB: convert time to year.decimal
#-----------------------------------------------------------------------------

fair = pd.read_csv(fair_file)
df_fair = pd.DataFrame()
df_fair['t_fair'] = fair.Year.values.astype('float')
df_fair['ts_fair'] = fair.Global.values

#------------------------------------------------------------------------------
# COMPUTE: climatological monthly mean (Ed's Climate Stripes: 1971-2000) from HadCRUT5 (monthly)
#------------------------------------------------------------------------------

mu_1851_1900 = np.nanmean( df_hadcrut5[(df_hadcrut5['t_hadcrut5']>=1851) & (df_hadcrut5['t_hadcrut5']<=1900)]['ts_hadcrut5'] ) # -0.507873106
mu_1961_1990 = np.nanmean( df_hadcrut5[(df_hadcrut5['t_hadcrut5']>=1961) & (df_hadcrut5['t_hadcrut5']<=1990)]['ts_hadcrut5'] ) #  0.005547222
mu_1971_2000 = np.nanmean( df_hadcrut5[(df_hadcrut5['t_hadcrut5']>=1971) & (df_hadcrut5['t_hadcrut5']<=2000)]['ts_hadcrut5'] ) #  0.176816667

if baseline == 'baseline_1851_1900':
    mu = mu_1851_1900
    baseline_start = 1851
    baseline_end = 1900
elif baseline == 'baseline_1961_1990':
    mu = mu_1961_1990
    baseline_start = 1961
    baseline_end = 1990
else:    
    mu = mu_1971_2000
    baseline_start = 1971
    baseline_end = 2000
cbarstr = 'Anomaly (from ' + str(baseline_start) + '-' + str(baseline_end) +')'

#------------------------------------------------------------------------------
# MERGE: annual dataframes
#------------------------------------------------------------------------------

# PAGES2k + HadCRUT5:
#t = np.floor(df_pages2k.t_pages2k).append(np.floor(df_hadcrut5_yearly.t_hadcrut5)).reset_index(drop=True)
#ts = df_pages2k.ts_pages2k.append(df_hadcrut5_yearly.ts_hadcrut5).reset_index(drop=True)

# PAGES2k + HadCRUT5 + FaIR:
t = (np.floor(df_pages2k.t_pages2k).append(np.floor(df_hadcrut5_yearly.t_hadcrut5))).append(np.floor(df_fair.t_fair)).reset_index(drop=True)
ts = (df_pages2k.ts_pages2k.append(df_hadcrut5_yearly.ts_hadcrut5)).append(df_fair.ts_fair).reset_index(drop=True)

df = pd.DataFrame()
df['Year'] = t
df['Global'] = ts - mu

# COMPUTE: standard deviation of the annual average anomalies (1901-2000)
sigma = np.nanstd( df[(df['Year']>=1901) & (df['Year']<=2000)]['Global'] )

if use_only_hadcrut5 == True:
    x = (df[ (df['Year']>=1850) & (df['Year']<=2020) ]['Year']).values
    y = np.array( df[ (df['Year']>=1850) & (df['Year']<=2020) ]['Global'])
else:
    x = (df['Year']).values
    if use_smoothing == True:
        y = pd.Series(np.array( df['Global']) ).rolling(nsmooth,center=True).mean().values
    else:
        y = np.array( df['Global'] )
z = np.array(len(y)*[1.0])

dg = pd.DataFrame({'x':x,'y':y,'z':z})

if use_logarithm == True: 
    
    # EXPERIMENTAL: log10 ( with analytic continuation where y = ymin) 
    y_min = np.nanmin(y)
    mask = (y < 0.95 * y_min)
    y = np.log10( y-y_min )    
    y[mask] = np.min(y[ (x>1500) & (x<1800) ]) 				        # local epoch expectation minimum
    yinterp = (pd.Series(y).rolling(25,center=True).mean()).values	# re-smooth using in-fill at singular value(s)
    y = yinterp    
    cbarstr = 'log₁₀ (Anomaly-min)'

# DEFINE: colormap ------------------------------------------------------------

mask = np.isfinite(y)
y_min = y[mask].min()
y_max = y[mask].max()    
y_ptp = y[mask].ptp()    
y_norm = ((y[mask] - y_min ) / y_ptp)       # [0,1]
#y_norm = [i / max(y) for i in y]           # [ymin/ymax,1]

x = x[mask]
y = y[mask]
z = z[mask]

#cmap = plt.cm.get_cmap('RdBu_r')
cmap = plt.cm.get_cmap('bwr')

# 'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 
# 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 
# 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 
# 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 
# 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 
# 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 
# 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 
# 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 
# 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 
# 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 
# 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 
# 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 
# 'cmo.algae', 'cmo.algae_r', 'cmo.amp', 'cmo.amp_r', 'cmo.balance', 'cmo.balance_r', 
# 'cmo.curl', 'cmo.curl_r', 'cmo.deep', 'cmo.deep_r', 'cmo.delta', 
# 'cmo.delta_r', 'cmo.dense', 'cmo.dense_r', 'cmo.diff', 'cmo.diff_r', 
# 'cmo.gray', 'cmo.gray_r', 'cmo.haline', 'cmo.haline_r', 'cmo.ice', 
# 'cmo.ice_r', 'cmo.matter', 'cmo.matter_r', 'cmo.oxy', 'cmo.oxy_r', 'cmo.phase', 
# 'cmo.phase_r', 'cmo.rain', 'cmo.rain_r', 'cmo.solar', 'cmo.solar_r', 
# 'cmo.speed', 'cmo.speed_r', 'cmo.tarn', 'cmo.tarn_r', 'cmo.tempo', 
# 'cmo.tempo_r', 'cmo.thermal', 'cmo.thermal_r', 'cmo.topo', 
# 'cmo.topo_r', 'cmo.turbid', 'cmo.turbid_r', 
# 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 
# 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 
# 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 
# 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 
# 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 
# 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 
# 'hot_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 
# 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 
# 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 
# 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 
# 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 
# 'terrain', 'terrain_r', 'turbo', 'turbo_r', 'twilight', 'twilight_r', 
# 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'winter', 'winter_r'

#maxval = +2.6 * sigma
#minval = -2.6 * sigma
#maxval = +0.75
#minval = -0.75
maxval = np.nanmax( y )
minval = np.nanmin( y )
midval = minval + (maxval-minval)/2
if use_only_hadcrut5 == True:
    colors = cmap( midval + y_norm )
else:
    colors = cmap( np.abs(np.nanmin(y)) + y )
norm = mcolors.TwoSlopeNorm(vmin=minval, vcenter=0., vmax=maxval)
sm = ScalarMappable( cmap=cmap, norm=norm )
sm.set_array([])

#plt.plot(t,ts,'grey')
#plt.plot(x,y,'red')
#plt.plot(x,y_norm,'blue')
#plt.scatter( x, y, c=colors, cmap=cmap, norm=norm)
#plt.axhline(y=midval)
#plt.xlim(1750,2100)

#------------------------------------------------------------------------------
# PLOTS
#------------------------------------------------------------------------------

# PLOT (1): timeseries

figstr = 'climate-timeseries' + '-' + projectionstr + '-' + str(nsmooth).zfill(2) + 'yr-smooth' + '.png'
titlestr = 'Mean annual anomaly (global: ' + projectionstr + '): 500-2200 AD'

fig, ax = plt.subplots( figsize=(15,10) ); ax.axis('off')
plt.scatter( x, y, c=colors, cmap=cmap, norm=norm, zorder=1 )
if use_overlay_timeseries == True: plt.plot( x, y, color='grey', ls='-', lw=0.5, zorder=0 )    
if use_overlay_baseline == True: plt.fill_betweenx( [ np.nanmin(y), np.nanmax(y) ], baseline_start, baseline_end, facecolor='grey', alpha=0.2, zorder=0 )        
if use_overlay_axis == True: ax.axis('on')
if use_overlay_colorbar == True:    
    cbar = plt.colorbar( sm, shrink=0.5, extend='both' )
    cbar.set_label( cbarstr, rotation=270, labelpad=25, fontsize=14 )
if use_only_hadcrut5 == True:
    plt.plot( [1850, 2020], [0,0], lw=0.5, color='white' )
    plt.text( 1847, 0.02, '1850', weight='bold' )    
    plt.text( 2017, 0.02, '2020', weight='bold' )        
    titlestr = 'Mean annual anomaly (global): 1850-2020 AD'
    figstr = 'climate-timeseries-1850-2020.png'
else:
    plt.plot( [500,2200], [0,0], lw=0.5, color='white' )
    plt.text( 485, 0.02, '500', weight='bold' )    
    plt.text( 970, 0.02, '1000', weight='bold' )    
    plt.text( 1470, 0.02, '1500', weight='bold' )    
    plt.text( 1820, 0.02, '1850', weight='bold' )    
    plt.text( 1990, 0.02, '2020', weight='bold' )    
    plt.text( 2170, 0.02, '2200', weight='bold' )    
plt.title( titlestr, fontsize=fontsize )        
plt.savefig( figstr, dpi=300 )
plt.close(fig)

# PLOT (2): climate bars ( bar chart )

figstr = 'climate-bars' + '-' + projectionstr + '-' + str(nsmooth).zfill(2) + 'yr-smooth' + '.png'
titlestr = 'Mean annual anomaly (global: ' + projectionstr + '): 500-2200 AD'

fig, ax = plt.subplots( figsize=(15,10) ); ax.axis('off')
plt.bar( x, y, color=colors, width=1.0 )
if use_overlay_baseline == True: plt.fill_betweenx( [ np.nanmin(y), np.nanmax(y) ], baseline_start, baseline_end, facecolor='grey', alpha=0.2, zorder=0 )        
if use_overlay_axis == True: ax.axis('on')
if use_overlay_colorbar == True:    
    cbar = plt.colorbar( sm, shrink=0.5, extend='both' )
    cbar.set_label( cbarstr, rotation=270, labelpad=25, fontsize=14 )
if use_only_hadcrut5 == True:
    if use_overlay_timeseries == True: plt.step( list(x)+[x[-1]], list(y)+[y[-1]], where='mid', color='grey', ls='-', lw=1 )
    titlestr = 'Mean annual anomaly (global): 1850-2020 AD'
    figstr = 'climate-bars-1850-2020.png'
else:
    if use_overlay_timeseries == True: plt.plot( x, y, color='grey', ls='-', lw=1 )
    plt.text( 470, 0.02, '500', weight='bold' )    
    plt.text( 960, 0.02, '1000', weight='bold' )    
    plt.text( 1460, 0.02, '1500', weight='bold' )    
    plt.text( 1810, 0.02, '1850', weight='bold' )    
    plt.text( 1980, 0.02, '2020', weight='bold' )    
    plt.text( 2160, 0.02, '2200', weight='bold' )    
plt.title( titlestr, fontsize=fontsize )        
plt.savefig( figstr, dpi=300 )
plt.close(fig)

# PLOT (3): climate stripes

figstr = 'climate-stripes' + '-' + projectionstr + '-' + str(nsmooth).zfill(2) + 'yr-smooth' + '.png'
titlestr = 'Mean annual anomaly (global: ' + projectionstr + '): 500-2200 AD'

fig, ax = plt.subplots( figsize=(15,10) ); ax.axis('off')
plt.bar( x, z, color=colors, width=1.0, zorder=0 )
if use_overlay_axis == True: ax.axis('on')
if use_overlay_colorbar == True:
    cbar = plt.colorbar( sm, shrink=0.5, extend='both' )
    cbar.set_label( cbarstr, rotation=270, labelpad=25, fontsize=14 )
if use_only_hadcrut5 == True:
    if use_overlay_timeseries == True: plt.step( list(x)+[x[-1]], list(y_norm)+[y_norm[-1]], where='mid', color='black', ls='-', lw=1 )
    titlestr = 'Mean annual anomaly (global): 1850-2020 AD'
    figstr = 'climate-stripes-1850-2020.png'
else:
    if use_overlay_timeseries == True: plt.plot( x, y_norm, color='black', ls='-', lw=1 )
    plt.text( 470, -0.02, '500', weight='bold' )    
    plt.text( 960, -0.02, '1000', weight='bold' )    
    plt.text( 1460, -0.02, '1500', weight='bold' )    
    plt.text( 1810, -0.02, '1850', weight='bold' )    
    plt.text( 1980, -0.02, '2020', weight='bold' )    
    plt.text( 2160, -0.02, '2200', weight='bold' )    
plt.title( titlestr, fontsize=fontsize )        
plt.savefig( figstr, dpi=300 )
plt.close(fig)

#------------------------------------------------------------------------------
print('** END')
