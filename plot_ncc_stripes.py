#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: plot-glosat-stripes.py
#------------------------------------------------------------------------------
# Version 0.2
# 10 September, 2021
# Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# IMPORT PYTHON LIBRARIES
#------------------------------------------------------------------------------
import numpy as np
import pandas as pd
# Plotting libraries:
import matplotlib
import matplotlib.pyplot as plt; plt.close('all')
import matplotlib.cm as cm
from matplotlib import colors as mcol
from matplotlib.cm import ScalarMappable
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

import cmocean as cmo


# Silence library version notifications
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# SETTINGS: 
#------------------------------------------------------------------------------

fontsize = 20
use_horizontal_colorbar = False
plot_hemispheres = False

#------------------------------------------------------------------------------
# I/O: GloSAT.prelim01_reqSD_alternativegrid-178101-201912.timeseries.txt
#------------------------------------------------------------------------------

# LOAD: temperature anomaly (from 1961-1990) CSV into pandas dataframe --> df
pathstr = 'DATA/'
filenamestr = 'GloSAT.prelim01_reqSD_alternativegrid-178101-201912.timeseries.txt'
filename_txt = pathstr + filenamestr

headings = ['Year', 'Month', 'Global', 'NH', 'SH']
df = pd.DataFrame(columns = headings)
datain = pd.read_csv(filename_txt, delim_whitespace=True, header=4)
for i in range(len(df.columns)):
    df[df.columns[i]] = datain.values[:,i]    

# TRIM: 1850+
df = df[df['Year']>=1850]

# CONVERT: year and month columns to integer and add day=15 column
df['Day'] = 15
df[['Year','Month','Day']] = df[['Year','Month','Day']].applymap(np.int64) 

# CONVERT: (year, month, day) --> datetime 
df['Date'] = pd.to_datetime(df[['Year', 'Month', 'Day']])

# REPLACE: fillval with NaN
fillval = -9.99900
df.replace(fillval, np.NaN, inplace=True) 

# SAVE: data munged version to CSV
df.to_csv('data_munged.csv', sep=',', index=False, header=True, encoding='utf-8')

# COMPUTE: annual means dataframe --> da
da = pd.DataFrame()
da['Year'] = df.groupby(df['Year'].dropna()).mean()['Year']
da['Global_Annual_Mean'] = df.groupby(df['Year'].dropna()).mean()['Global'].values
da['NH_Annual_Mean'] = df.groupby(df['Year'].dropna()).mean()['NH'].values
da['SH_Annual_Mean'] = df.groupby(df['Year'].dropna()).mean()['SH'].values

if plot_hemispheres == True:

    # PLOT: Mean annual anomaly
    
    fig, ax = plt.subplots(figsize=(15,10))
    plt.bar(da['Year'], da['NH_Annual_Mean'], label='NH', color='pink')
    plt.bar(da['Year'], da['SH_Annual_Mean'], label='SH', color='lightblue')
    plt.plot(da['Year'], da['Global_Annual_Mean'], label='Global', color='black', linewidth=2)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.ylabel('Mean annual anomaly (from 1961-1990)', fontsize=fontsize)
    plt.legend(fontsize=fontsize)
    plt.title(filename_txt, fontsize=fontsize)
    plt.savefig(filenamestr + '_anomly_annual_mean.png')
    plt.close(fig)

#------------------------------------------------------------------------------
# climate stripes 
#------------------------------------------------------------------------------

# da = da[da['Year']>1978] # ERA-5

# COMPUTE: 1961-1990 monthly mean (Ed's Climate Stripes: 1971-2000)
mu_1851_1900 = np.nanmean( df[(df['Year']>1850) & (df['Year']<1901)]['Global'] ) # -0.507873106
mu_1961_1990 = np.nanmean( df[(df['Year']>1960) & (df['Year']<1991)]['Global'] ) #  0.005547222
mu_1971_2000 = np.nanmean( df[(df['Year']>1970) & (df['Year']<2001)]['Global'] ) #  0.176816667

mu = mu_1971_2000

# COMPUTE: standard deviation of the annual average anomalies (1901-2000)
sigma = np.nanstd( da[(da['Year']>1900) & (da['Year']<2001)]['Global_Annual_Mean'] )

x = da[da['Year']>1850]['Year']
y = np.array(da[da['Year']>1850]['Global_Annual_Mean'] - mu)
z = len(y)*[1.0]


mask = np.isfinite(y)
y_min = y[mask].min()    
y_max = y[mask].max()    
y_ptp = y[mask].ptp()    
y_norm = ((y - y_min) / y_ptp)    
    
# DEFINE: colormap

cmap = plt.cm.get_cmap('RdBu_r')

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

maxval = y_max
minval = -y_max

#colors = cmap( y/maxval + 0.5)
#sm = ScalarMappable(cmap=cmap, norm=plt.Normalize(np.nanmin(y),np.nanmax(y)))
#sm.set_array([])

colors = cmap( y_norm )
sm = ScalarMappable( cmap=cmap, norm=plt.Normalize(minval,maxval) )

# PLOT: mean annual anomaly (1900-2019) as climate stripe bars

fig, ax = plt.subplots(figsize=(15,10))
ax.axis('off')
plt.bar(x, y, color=colors, width=1.0)
if use_horizontal_colorbar == True:
    cbar = plt.colorbar(sm, shrink=0.5, orientation='horizontal')
    cbar.set_label('Anomaly (from 1971-2000)', rotation=0, labelpad=25, fontsize=fontsize)
else:
    cbar = plt.colorbar(sm, shrink=0.5)
    cbar.set_label('Anomaly (from 1971-2000)', rotation=270, labelpad=25, fontsize=fontsize)
plt.title('Mean annual anomaly: Global', fontsize=fontsize)
plt.savefig(filenamestr + '_anomaly-bars.png')
plt.close(fig)

# PLOT: mean annual anomaly (1900-2019) as climate stripes

fig, ax = plt.subplots(figsize=(15,10))
ax.axis('off')
plt.bar(x, z, color=colors, width=1.0)
if use_horizontal_colorbar == True:
    cbar = plt.colorbar(sm, shrink=0.5, orientation='horizontal')
    cbar.set_label('Mean annual anomaly (from 1971-2000)', rotation=0, labelpad=25, fontsize=fontsize)
else:
    cbar = plt.colorbar(sm, shrink=0.5)
    cbar.set_label('Anomaly (from 1971-2000)', rotation=270, labelpad=25, fontsize=fontsize)
plt.title('Mean annual anomaly: Global', fontsize=fontsize)
plt.savefig(filenamestr + '_anomaly-stripes.png')
plt.close(fig)

#------------------------------------------------------------------------------
print('** END')
