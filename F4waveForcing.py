import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.interpolate import interp1d
import os
from netCDF4 import Dataset
from itertools import groupby
import datetime as DT
import pickle
import copy
import sandbarUtils


def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.
    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="->", color=color),
        size=size
    )

def moving_average(a, n=21):
   ret = np.cumsum(a, dtype=float)
   ret[n:] = ret[n:] - ret[:-n]
   return ret[n - 1:] / n



dbfile2 = open('ceofsSouthernTransect.pickle', 'rb')
data2 = pickle.load(dbfile2)
dbfile2.close()

timeS = data2['time']
alllinesS = data2['alllines']
xinterpS = data2['xinterp']
SS = data2['S']
RS = data2['Rt']
thetaS = data2['thetaRadians']
theta2S = data2['thetaDegrees']
phitS = data2['phiRadian']
phit2S = data2['phiDegrees']
totalVS = data2['totalV']
lambaS = data2['lamda']
percentVS = data2['percentV']


time3 = timeS[1:-1]



thetasMode1 = np.arange(-180,360,5)
thetasMode2 = np.arange(-180,360,5)
meshTheta1, meshTheta2 = np.meshgrid(thetasMode1,thetasMode2)

thetasMode1 = np.arange(-180,360,5)
thetasMode2 = np.arange(-180,360,5)
outerBarX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterBarX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))

for tt in range(len(thetasMode1)):
   eofPred = SS[:, 0] * np.nanmean(RS[:, 0]) * np.cos(thetasMode1[tt] * np.pi / 180 - thetaS[:, 0])
   for yy in range(len(thetasMode2)):
      eofPred2 = SS[:, 1] * np.nanmean(RS[:, 1]) * np.cos(thetasMode2[yy] * np.pi / 180 - thetaS[:, 1])
      combined = np.real(np.mean(alllinesS,axis=0)+eofPred+eofPred2)
      xFRFbar, xFRFtrough = sandbarUtils.findSandBarAndTrough1D(xinterpS, combined, plotFname=None, smoothLengthScale=5, profileTrend=np.mean(alllinesS,axis=0))
      if xFRFbar is not None:
         if len(xFRFbar) > 1:
            outerBarX[yy,tt] = np.max(xFRFbar)
         else:
            innerOuterBarX[yy,tt] = xFRFbar

doublebar = copy.deepcopy(outerBarX)
doubleIndex = np.where((doublebar > 0))
doublebar[doubleIndex] = 0.5

singlebar = copy.deepcopy(innerOuterBarX)
singleIndex = np.where((singlebar > 0))
doublebar[singleIndex] = 0.75
wavedir = '/media/dylananderson/Elements/WIS_ST63218/'

# Need to sort the files to ensure correct temporal order...
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

wavedir26 = '/media/dylananderson/Elements/26mArrayPlus17/'
# Need to sort the files to ensure correct temporal order...
files = os.listdir(wavedir26)
files.sort()
files_path = [os.path.join(os.path.abspath(wavedir26), x) for x in files]
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
hsSmooth = moving_average(hsCombined,3)
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


hsOrdinaltime = np.array([i.toordinal() for i in tC])
time3Ordinal = np.array([i.toordinal() for i in time3])
avgMonthBeforeHs = np.nan * np.ones((len(time3),))
avgMonthBeforeWE = np.nan * np.ones((len(time3),))
avgMonthBeforeFV = np.nan * np.ones((len(time3),))
avgMonthBeforeLWP = np.nan * np.ones((len(time3),))

for ii in range(len(time3)):
    # beforeIndex = np.where((tC > DT.datetime(time5[ii].year,time5[ii].month,time5[ii].day)) & (tC < time5[ii]))
    beforeIndex = np.where((hsOrdinaltime > (time3Ordinal[ii]-30)) & (hsOrdinaltime < (time3Ordinal[ii])))

    tempHs = hsCombined[beforeIndex]
    avgMonthBeforeHs[ii] = np.nanmean(tempHs)
    tempWE = weC[beforeIndex]
    avgMonthBeforeWE[ii] = np.nanmean(tempWE)
    tempLWP = lwpC[beforeIndex]
    avgMonthBeforeLWP[ii] = np.nanmean(np.abs(tempLWP))
    tempFV = fV[beforeIndex]
    avgMonthBeforeFV[ii] = np.nanmean(tempFV)


lowBreak = 6
medLowBreaka = 6
medLowBreakb = 7.6
medHighBreaka = 7.6
medHighBreakb = 9.1
highBreaka = 9.1
highBreakb = 17
highesta = 18

binClassHs = np.nan * np.ones((len(time3),))
for ii in range(len(time3)):
    fvTemp = avgMonthBeforeFV[ii]
    if fvTemp > medLowBreaka and fvTemp < medLowBreakb:
        binClassHs[ii] = 1
    elif fvTemp >= medHighBreaka and fvTemp < medHighBreakb:
        binClassHs[ii] = 2
    elif fvTemp >= highBreaka and fvTemp < highBreakb:
        binClassHs[ii] = 3
    elif fvTemp >= highesta:
        binClassHs[ii] = 4
    elif fvTemp < lowBreak:
        binClassHs[ii] = 0

grouped_L = [(k, sum(1 for i in g)) for k,g in groupby(binClassHs)]

def repeatingNumbers(numList):
    i = 0
    output = []
    while i < len(numList) - 1:
        n = numList[i]
        startIndex = i
        while i < len(numList) - 1 and numList[i] == numList[i + 1]:
            i = i + 1

        endIndex = i

        output.append([n, startIndex, endIndex, (endIndex-startIndex + 1)])
        print("{0} >> {1}".format(n, [startIndex, endIndex]))
        i = i + 1


    return output

groupedOutput = repeatingNumbers(binClassHs)
lowWaves = []
medLowWaves = []
medHighWaves = []
bigWaves = []
highestWaves = []

for tt in range(len(groupedOutput)):
    if groupedOutput[tt][3] > 1:
        if groupedOutput[tt][0] == 0:
            lowWaves.append([groupedOutput[tt][1],groupedOutput[tt][2]])
        if groupedOutput[tt][0] == 1:
            medLowWaves.append([groupedOutput[tt][1], groupedOutput[tt][2]])
        if groupedOutput[tt][0] == 2:
            medHighWaves.append([groupedOutput[tt][1],groupedOutput[tt][2]])
        if groupedOutput[tt][0] == 3:
            bigWaves.append([groupedOutput[tt][1], groupedOutput[tt][2]])


fig = plt.figure(figsize=(10,10))

# ax100 = plt.subplot2grid((1,3),(0,0),rowspan=1,colspan=1)
ax100 = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)

con1 = ax100.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
coloredlines = cm.rainbow(np.linspace(0,1,len(lowWaves)))

lineSlopes = []
lineSlopeTimes = []
for ff in range(len(lowWaves)):
    si = lowWaves[ff][0]
    ei = lowWaves[ff][1]
    # si = medLowWaves[ff][0]
    # ei = medLowWaves[ff][1]
    # si = medHighWaves[ff][0]
    # ei = medHighWaves[ff][1]
    # si = bigWaves[ff][0]
    # ei = bigWaves[ff][1]
    #indexedWaves = avgMonthBeforeFV[si-1:ei+2]
    indexedTime = time3[si-2:ei]

    indexedOrdinalTime = [i.toordinal() for i in indexedTime]
    indexedTimeMonth = [i.month for i in indexedTime]
    mode1prop = mode1Phase3WA[si-2:ei]
    mode2prop = mode2Phase3WA[si-2:ei]
    mode1propWrap = np.unwrap((mode1prop + np.pi) % (2 * np.pi) - np.pi)
    mode2propWrap = np.unwrap((mode2prop + np.pi) % (2 * np.pi) - np.pi)
    x1 = np.sort(mode1propWrap*180/np.pi)
    y1 = mode2propWrap*180/np.pi
    mags = mode1Mag[si:ei+2]

    y1 = np.sort(y1)
    rise = np.abs(y1[-1]-y1[0])
    run = np.abs(x1[-1]-x1[0])

    slope = rise/run

    colormap = cm.coolwarm
    # colormap = cm.jet
    #normalize = mcolors.Normalize(vmin=0.33, vmax=3)
    normalize = mcolors.LogNorm(vmin=0.5, vmax=3)

    s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
    colorOfLine = colormap(normalize(slope))

    if x1[-1] < 5 and x1[-1] > 0:
        x1 = x1+30
    if y1[-1] < -20:
        y1 = y1 + 360
    newx = np.arange(x1[0],x1[-1],1)
    f = interp1d(x1,y1,kind='linear')
    newy = f(newx)

    if indexedTime[0].year < 2005 or indexedTime[0].year > 2008:
        if indexedTime[0].year < 2021:
            if np.nanmean(mags) > 0.4:

                if x1[-1] < 60 and y1[0]>-100:

                    lines1 = ax100.plot(newx,newy,'-',linewidth=2,color=colorOfLine,label='{}'.format(indexedTime[0].year))
                    lineSlopeTimes.append(indexedTime[0])
                    lineSlopes.append(slope)
                    add_arrow(lines1[0], color='k', size=25)
                elif x1[-1] > 60:
                    if x1[-1] > 200 and y1[-1]>200:
                        lines1 = ax100.plot(newx-360,newy,'-',linewidth=2,color=colorOfLine,label='{}'.format(indexedTime[0].year))
                        add_arrow(lines1[0], color='k', size=25)
                        lineSlopeTimes.append(indexedTime[0])
                        lineSlopes.append(slope)
                    else:
                        if indexedTime[0].year < 2014 and indexedTime[0].year > 1981:
                            if x1[0] < 73 and y1[0] < -65:
                                print('skipping')
                            else:
                                lines1 = ax100.plot(newx,newy,'-',linewidth=2,color=colorOfLine,label='{}'.format(indexedTime[0].year))
                                lineSlopeTimes.append(indexedTime[0])
                                lineSlopes.append(slope)
                            add_arrow(lines1[0], color='k', size=25)

ax100.set_xlim([-180, 280])
ax100.set_ylim([-160, 320])
ax100.set_xlabel(r'Offshore Propagation $\Longrightarrow$',fontsize=12)
ax100.set_ylabel(r'Onshore Propagation $\Longrightarrow$',fontsize=12)
ax100.set_title('Non-dimensional Fall Velocity < 6')
# plt.legend()
plt.show()




