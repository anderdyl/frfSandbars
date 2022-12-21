import numpy as np
import matplotlib.pyplot as plt
import rasterio
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import datetime as DT
import os
import cv2
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib.patches as mpatches
import pickle
from datetime import datetime
from datetime import timedelta

# The following two files must downloaded from NAIP before running this file
# m_3607551_sw_18_1_20141005_20141118.jp2
# m_3607551_sw_18_1_20141005_20141118.jp2

fig = plt.figure(figsize = (10, 8))
jp21 = 'm_3607551_sw_18_1_20141005_20141118.jp2'
jp22 = 'm_3607550_se_18_1_20141005_20141118.jp2'

with rasterio.open(jp21) as f:
    im1r = f.read(1)
    im1g = f.read(2)
    im1b = f.read(3)
with rasterio.open(jp22) as f:
    im2r = f.read(1)
    im2g = f.read(2)
    im2b = f.read(3)

dataIM1r = np.array(im1r, dtype=im1r[0].dtype)
dataIM1g = np.array(im1g, dtype=im1g[0].dtype)
dataIM1b = np.array(im1b, dtype=im1b[0].dtype)
dataIM2r = np.array(im2r, dtype=im2r[0].dtype)
dataIM2g = np.array(im2g, dtype=im2g[0].dtype)
dataIM2b = np.array(im2b, dtype=im2b[0].dtype)

mergedr = np.hstack((dataIM2r,dataIM1r))
mergedg = np.hstack((dataIM2g,dataIM1g))
mergedb = np.hstack((dataIM2b,dataIM1b))
n, m = np.shape(mergedb)

rgbArray = np.zeros((n, m, 3), 'uint8')
rgbArray[..., 0] = mergedr  # imR[:, :, 0] * 255
rgbArray[..., 1] = mergedg  # imG[:, :, 0] * 255
rgbArray[..., 2] = mergedb  # imB[:, :, 0] * 255
imGrayscale = cv2.cvtColor(rgbArray, cv2.COLOR_BGR2GRAY)
ax2 = plt.subplot2grid((3,4),(0,0),rowspan=1,colspan=2)
ax2.imshow(rgbArray)

ax2.set_xlim([5500,8580])
ax2.set_ylim([2000,0])
ax2.set_xticks([])
ax2.set_yticks([])
ax2.plot([7040,8290],[1350,980],'r--')
ax2.plot([6950,8200],[1125,755],'r--')
ax2.plot([6950,7040],[1125,1350],'r--')
ax2.plot([8200,8290],[755,980],'r--')

ax2.add_artist(ScaleBar(dx=0.75,units='m'))

x_tail = 8290
y_tail = 1750
x_head = 8290
y_head = 1350
dx = x_head - x_tail
dy = y_head - y_tail
arrow = mpatches.FancyArrowPatch((x_tail, y_tail), (x_head, y_head),
                                 mutation_scale=16,color="k")
ax2.add_patch(arrow)
ax2.text(-0.05, 1.03, 'a.', transform=ax2.transAxes, size=14, weight='bold')

dbfile = open('sandbarsSouthernTransect_referencedMHHW.pickle', 'rb')
data = pickle.load(dbfile)
dbfile.close()
alllines = data['alllines']
xinterp = data['xinterp']
time = data['time']
ax3 = plt.subplot2grid((3,4),(0,2),rowspan=1,colspan=2)

# # Plotting all of the lines
t1 = 0
t2 = -1
ax3.set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(alllines)))))
for i in range(len(alllines)):
    ax3.plot(xinterp, alllines[i,:], color=[0.65,0.65,0.65],linewidth=0.5)# label=time[i])
ax3.plot(xinterp,np.mean(alllines,axis=0),'k')
ax3.plot(xinterp,np.mean(alllines,axis=0)+np.std(alllines,axis=0),'k--')
ax3.plot(xinterp,np.mean(alllines,axis=0)-np.std(alllines,axis=0),'k--')
ax3.set_xlim([0, 500])
ax3.set_ylim([-7.5, 1])
ax3.set_ylabel('Elevation (m, MHHW)')
ax3.set_xlabel('Cross-shore (m)')
ax3.text(-0.05, 1.03, 'b.', transform=ax3.transAxes, size=14, weight='bold')
plt.set_cmap('RdBu_r')
tg, xg = np.meshgrid(time, xinterp)
ax4 = plt.subplot2grid((3,4),(1,0),rowspan=1,colspan=4)
plt4 = ax4.pcolor(tg,xg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
cbaxes = inset_axes(ax4, width="30%", height="3%", loc=1)
cb1 = plt.colorbar(plt4, cax=cbaxes, orientation='horizontal')
cbaxes.set_xlabel('Elevation Anomaly (m)')
ax4.set_ylabel('Cross-shore (m)')
ax4.set_xlim([DT.datetime(1981,1,1),DT.datetime(2020,12,1)])
ax4.text(-0.05, 1.03, 'c.', transform=ax4.transAxes, size=14, weight='bold')
plt.style.use('default')
plt.set_cmap('RdBu_r')

wavedir = '/users/dylananderson/Documents/data/WIS_ST63218/'
files = os.listdir(wavedir)
files.sort()
files_path = [os.path.join(os.path.abspath(wavedir), x) for x in files]
wis = Dataset(files_path[0])

def getWIS(file):
    waves = Dataset(file)

    waveHs = waves.variables['waveHs'][:]
    waveTp = waves.variables['waveTp'][:]
    waveMeanDirection = waves.variables['waveMeanDirection'][:]
    waveTm = waves.variables['waveTm'][:]
    waveTm1 = waves.variables['waveTm1'][:]
    waveTm2 = waves.variables['waveTm2'][:]
    waveHsWindsea = waves.variables['waveHsWindsea'][:]
    waveTmWindsea = waves.variables['waveTmWindsea'][:]
    waveMeanDirectionWindsea = waves.variables['waveMeanDirectionWindsea'][:]
    waveSpreadWindsea = waves.variables['waveSpreadWindsea'][:]
    timeW = waves.variables['time'][:]
    waveTpSwell = waves.variables['waveTpSwell'][:]
    waveHsSwell = waves.variables['waveHsSwell'][:]
    waveMeanDirectionSwell = waves.variables['waveMeanDirectionSwell'][:]
    waveSpreadSwell = waves.variables['waveSpreadSwell'][:]
    output = dict()
    output['waveHs'] = waveHs
    output['waveTp'] = waveTp
    output['waveMeanDirection'] = waveMeanDirection
    output['waveTm'] = waveTm
    output['waveTm1'] = waveTm1
    output['waveTm2'] = waveTm2
    output['waveTpSwell'] = waveTpSwell
    output['waveHsSwell'] = waveHsSwell
    output['waveMeanDirectionSwell'] = waveMeanDirectionSwell
    output['waveSpreadSwell'] = waveSpreadSwell
    output['waveHsWindsea'] = waveHsWindsea
    output['waveTpWindsea'] = waveTmWindsea
    output['waveMeanDirectionWindsea'] = waveMeanDirectionWindsea
    output['waveSpreadWindsea'] = waveSpreadWindsea
    output['t'] = timeW
    return output

def datenum_to_datetime(datenum):
    """
    Convert Matlab datenum into Python datetime.
    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    days = datenum % 1
    hours = days % 1 * 24
    minutes = hours % 1 * 60
    seconds = minutes % 1 * 60
    return datetime.fromordinal(int(datenum)) \
           + timedelta(days=int(days)) \
           + timedelta(hours=int(hours)) \
           + timedelta(minutes=int(minutes)) \
           + timedelta(seconds=round(seconds)) \
           - timedelta(days=366)

Hs = []
Tp = []
Dm = []
hsSwell = []
tpSwell = []
dmSwell = []
hsWindsea = []
tpWindsea = []
dmWindsea = []

timeWave = []
for i in files_path:
    waves = getWIS(i)
    Hs = np.append(Hs,waves['waveHs'])
    Tp = np.append(Tp,waves['waveTp'])
    Dm = np.append(Dm,waves['waveMeanDirection'])
    hsSwell = np.append(hsSwell,waves['waveHsSwell'])
    tpSwell = np.append(tpSwell,waves['waveTpSwell'])
    dmSwell = np.append(dmSwell,waves['waveMeanDirectionSwell'])
    hsWindsea = np.append(hsWindsea,waves['waveHsWindsea'])
    tpWindsea = np.append(tpWindsea,waves['waveTpWindsea'])
    dmWindsea = np.append(dmWindsea,waves['waveMeanDirectionWindsea'])
    timeWave = np.append(timeWave,waves['t'].flatten())


def getArray(file):
    waves = Dataset(file)
    waveHs = waves.variables['waveHs'][:]
    waveTp = waves.variables['waveTp'][:]
    waveMeanDirection = waves.variables['waveMeanDirection'][:]
    timeW = waves.variables['time'][:]
    output = dict()
    output['waveHs'] = waveHs
    output['waveTp'] = waveTp
    output['waveMeanDirection'] = waveMeanDirection
    output['t'] = timeW
    return output

wavedir26 = '/users/dylananderson/Documents/data/26mArrayPlus17/'

# Need to sort the files to ensure correct temporal order...
files = os.listdir(wavedir26)
files.sort()
files_path = [os.path.join(os.path.abspath(wavedir26), x) for x in files[1:]]
array26m = Dataset(files_path[0])
Hs26m = []
Tp26m = []
Dm26m = []
timeWave26m = []
for i in files_path:
    waves26m = getArray(i)
    Hs26m = np.append(Hs26m,waves26m['waveHs'][0:-1:2])
    Tp26m = np.append(Tp26m,waves26m['waveTp'][0:-1:2])
    Dm26m = np.append(Dm26m,waves26m['waveMeanDirection'][0:-1:2])
    timeWave26m = np.append(timeWave26m,waves26m['t'][0:-1:2])

ind = np.where((Hs26m > 0))
hs26m = Hs26m[ind]
tp26m = Tp26m[ind]
dm26m = Dm26m[ind]
t26m = timeWave26m[ind]
tWave26m = [DT.datetime.fromtimestamp(x) for x in t26m]

hsCombined = np.append(Hs,hs26m)
#hsSmooth = moving_average(hsCombined,3)
tpCombined = np.append(Tp,tp26m)
dmCombined = np.append(Dm,dm26m)
tWave = [DT.datetime.fromtimestamp(x) for x in timeWave]
tC = np.append(np.array(tWave),tWave26m)

badDirs = np.where((dmCombined > 360))
dmCombined[badDirs] = dmCombined[badDirs]*np.nan
badtp = np.where((tpCombined < 1))
tpCombined[badtp] = tpCombined[badtp]*np.nan

waveNorm = dmCombined - 72
neg = np.where((waveNorm > 180))
waveNorm[neg[0]] = waveNorm[neg[0]]-360
offpos = np.where((waveNorm>90))
offneg = np.where((waveNorm<-90))
waveNorm[offpos[0]] = waveNorm[offpos[0]]*0
waveNorm[offneg[0]] = waveNorm[offneg[0]]*0

Lo = (9.81/(2*np.pi)) * np.square(tpCombined)
Ir = 0.122/(np.sqrt((hsCombined/Lo)))
HoverL = hsCombined/Lo
lwpC = 1025*np.square(hsCombined)*tpCombined*(9.81/(64*np.pi))*np.cos(waveNorm*(np.pi/180))*np.sin(waveNorm*(np.pi/180))
weC = np.square(hsCombined)*tpCombined
ws = (2.65-1)*9.81*np.square(0.00015)/(18*0.000001)
fV = hsCombined/(ws*tpCombined)

def moving_average(a, n=21):
   ret = np.cumsum(a, dtype=float)
   ret[n:] = ret[n:] - ret[:-n]
   return ret[n - 1:] / n

monthlyfV = moving_average(fV,24*30)

ax5 = plt.subplot2grid((3,4),(2,0),rowspan=1,colspan=4)
ax5.plot(tC,fV,color=[0.65,0.65,0.65],linewidth=0.5)
ax5.plot(tC[(15*24):-(15*24)+1],monthlyfV,'k')
ax5.set_xlim([DT.datetime(1981,1,1),DT.datetime(2020,12,1)])
ax5.set_ylabel('$\Omega$')
ax5.set_xlabel('Time')
ax5.text(-0.05, 1.03, 'd.', transform=ax5.transAxes, size=14, weight='bold')
plt.tight_layout()


