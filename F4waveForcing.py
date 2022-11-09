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

def weightedMovingAverage(a,b,n=3):
    cut = np.floor(n/2)
    index = np.arange(int(cut),int(len(a)-cut),1)
    output = np.nan * np.ones((len(index),))
    counter = 0
    for ff in index:
        subset = a[int(ff-cut):int(ff+cut+1)]
        weights = (b[int(ff-cut):int(ff+cut+1)])
        output[counter] = np.average(subset,weights=weights)
        counter = counter+1
    return output



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



thetasMode1 = np.arange(0,450,2)
thetasMode2 = np.arange(-50,450,2)
meshTheta1, meshTheta2 = np.meshgrid(thetasMode1,thetasMode2)

thetasMode1 = np.arange(0,450,2)
thetasMode2 = np.arange(-50,450,2)
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
doublebar[singleIndex] = 0.25
wavedir = '/Users/dylananderson/Documents/data/WIS_ST63218/'

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

wavedir26 = '/Users/dylananderson/Documents/data/26mArrayPlus17/'
# Need to sort the files to ensure correct temporal order...
files = os.listdir(wavedir26)
files.sort()
files_path = [os.path.join(os.path.abspath(wavedir26), x) for x in files if x.endswith(".nc")]
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


mode1Phase3WA = weightedMovingAverage(np.unwrap(phit2S[:,0]*(np.pi/180)),RS[:,0],n=3)
mode2Phase3WA = weightedMovingAverage(np.unwrap(phit2S[:,1]*(np.pi/180)),RS[:,1],n=3)

mode1Mag3WA = weightedMovingAverage(RS[:,0],RS[:,0],n=3)
mode2Mag3WA = weightedMovingAverage(RS[:,1],RS[:,1],n=3)

mode1Wrap3WA = (mode1Phase3WA + np.pi) % (2 * np.pi) - np.pi
mode2Wrap3WA = (mode2Phase3WA + np.pi) % (2 * np.pi) - np.pi

mode1Mag = moving_average(RS[:,0],n=3)
mode2Mag = moving_average(RS[:,1],n=3)

fig = plt.figure(figsize=(14,8))
ax100 = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1)
con1 = ax100.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
coloredlines = cm.rainbow(np.linspace(0,1,len(lowWaves)))
lineSlopes = []
lineSlopeTimes = []
for ff in range(len(lowWaves)):
    si = lowWaves[ff][0]
    ei = lowWaves[ff][1]
    indexedTime = time3[si-2:ei]
    print(si)
    print(ei)
    print(indexedTime)
    indexedOrdinalTime = [i.toordinal() for i in indexedTime]
    indexedTimeMonth = [i.month for i in indexedTime]
    mode1prop = mode1Phase3WA[si-2:ei]
    mode2prop = mode2Phase3WA[si-2:ei]
    # mode1prop = mode1Phase3WA[si-1:ei+1]
    # mode2prop = mode2Phase3WA[si-1:ei+1]
    mode1propWrap = np.unwrap((mode1prop + np.pi) % (2 * np.pi) - np.pi)
    mode2propWrap = np.unwrap((mode2prop + np.pi) % (2 * np.pi) - np.pi)
    x1 = np.sort(mode1propWrap*180/np.pi)
    y1 = mode2propWrap*180/np.pi
    mags1 = mode1Mag[si:ei+2]
    mags2 = mode2Mag[si:ei+2]

    if len(x1) > 8:
        x1a = x1[0:6]
        x1b = x1[6:]
        y1a = np.sort(y1[0:6])
        y1b = np.sort(y1[6:])
        risea = np.abs(y1a[-1] - y1a[0])
        runa = np.abs(x1a[-1] - x1a[0])
        slopea = risea / runa

        riseb = np.abs(y1b[-1] - y1b[0])
        runb = np.abs(x1b[-1] - x1b[0])
        slopeb = riseb / runb

        colormap = cm.coolwarm
        # colormap = cm.jet
        # normalize = mcolors.Normalize(vmin=0.33, vmax=3)
        normalize = mcolors.LogNorm(vmin=0.66, vmax=3)

        s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
        colorOfLinea = colormap(normalize(slopea))
        colorOfLineb = colormap(normalize(slopeb))
        if y1a[-1] < -50:
            y1a = y1a + 360
        if y1a[0] < 20 and x1a[-1] < -30:
            # y1 = y1 + 360
            x1a = x1a + 360
        if y1a[0] > 120 and x1a[0] < -50:
            x1a = x1a + 360
        newxa = np.arange(x1a[0],x1a[-1],1)
        fa = interp1d(x1a,y1a,kind='linear')
        newya = fa(newxa)
        lines1a = ax100.plot(newxa, newya, '-', linewidth=2, color=colorOfLinea, label='{}'.format(indexedTime[0].year))
        add_arrow(lines1a[0], color='k', size=25)
        if y1b[-1] < -50:
            y1b = y1b + 360
        if y1b[0] < 20 and x1b[-1] < -30:
            # y1 = y1 + 360
            x1b = x1b + 360
        if y1b[0] > 120 and x1b[0] < -50:
            x1b = x1b + 360
        if x1b[0] < 0:
            x1b = x1b + 360
        newxb = np.arange(x1b[0],x1b[-1],1)
        fb = interp1d(x1b,y1b,kind='linear')
        newyb = fb(newxb)
        lines1b = ax100.plot(newxb, newyb, '-', linewidth=2, color=colorOfLineb, label='{}'.format(indexedTime[0].year))
        add_arrow(lines1b[0], color='k', size=25)
    else:
        y1 = np.sort(y1)
        rise = np.abs(y1[-1]-y1[0])
        run = np.abs(x1[-1]-x1[0])

        slope = rise/run

        colormap = cm.coolwarm
        # colormap = cm.jet
        #normalize = mcolors.Normalize(vmin=0.33, vmax=3)
        normalize = mcolors.LogNorm(vmin=0.66, vmax=3)

        s_map = cm.ScalarMappable(norm=normalize, cmap=colormap)
        colorOfLine = colormap(normalize(slope))

        # if x1[-1] < 5 and x1[-1] > 0:
        #     x1 = x1+30
        if y1[0] < -50:
            y1 = y1 + 360
        if y1[0] < 0 and x1[-1] < -30:
            # y1 = y1 + 360
            x1 = x1 + 360
        if y1[0] > 120 and x1[0] < -50:
            x1 = x1 + 360
        if x1[0] < 0:
            x1 = x1 + 360
        newx = np.arange(x1[0],x1[-1],1)
        f = interp1d(x1,y1,kind='linear')
        newy = f(newx)
        if indexedTime[0].year < 2005 or indexedTime[0].year > 2008:# and indexedTime[0].year < 2021:
            if np.mean(mags1) < 0.23:
                #lines1 = ax100.plot(newx,newy,'-',linewidth=2,color='k',label='{}'.format(indexedTime[0].year))
                print(np.mean(mags1))
                print(np.mean(mags2))

            # elif np.mean(mags2) < 0.23:
            #     lines1 = ax100.plot(newx,newy,'-',linewidth=2,color='k',label='{}'.format(indexedTime[0].year))
            #     print(np.mean(mags1))
            #     print(np.mean(mags2))
            else:
                lines1 = ax100.plot(newx, newy, '-', linewidth=2, color=colorOfLine, label='{}'.format(indexedTime[0].year))
                add_arrow(lines1[0], color='k', size=25)
                lineSlopes.append(slope)
                lineSlopeTimes.append(indexedTime[0])
                print(indexedTime[0])

    #
    # if indexedTime[0].year < 2005 or indexedTime[0].year > 2008:
    #     if indexedTime[0].year < 2021:
    #         if np.nanmean(mags) > 0.4:
    #
    #             if x1[-1] < 60 and y1[0]>-100:
    #
    #                 lines1 = ax100.plot(newx,newy,'-',linewidth=2,color=colorOfLine,label='{}'.format(indexedTime[0].year))
    #                 lineSlopeTimes.append(indexedTime[0])
    #                 lineSlopes.append(slope)
    #                 add_arrow(lines1[0], color='k', size=25)
    #             elif x1[-1] > 60:
    #                 if x1[-1] > 200 and y1[-1]>200:
    #                     lines1 = ax100.plot(newx-360,newy,'-',linewidth=2,color=colorOfLine,label='{}'.format(indexedTime[0].year))
    #                     add_arrow(lines1[0], color='k', size=25)
    #                     lineSlopeTimes.append(indexedTime[0])
    #                     lineSlopes.append(slope)
    #                 else:
    #                     if indexedTime[0].year < 2014 and indexedTime[0].year > 1981:
    #                         if x1[0] < 73 and y1[0] < -65:
    #                             print('skipping')
    #                         else:
    #                             lines1 = ax100.plot(newx,newy,'-',linewidth=2,color=colorOfLine,label='{}'.format(indexedTime[0].year))
    #                             lineSlopeTimes.append(indexedTime[0])
    #                             lineSlopes.append(slope)
    #                         add_arrow(lines1[0], color='k', size=25)

ax100.set_xlim([0, 440])
ax100.set_ylim([-50, 380])
ax100.set_xlabel(r'Offshore Propagation $\Longrightarrow$',fontsize=10)
ax100.set_ylabel(r'Onshore Propagation $\Longrightarrow$',fontsize=10)
# ax100.set_title('Non-dimensional Fall Velocity < 6')
ax100.set_title('Migrations during mild waves')
ax100.text(-0.05, 1.03, 'a.', transform=ax100.transAxes, size=14, weight='bold')

# plt.legend()
#



timeSubset= timeS#[1:]
lineSubset = []
threeMonthWE = []
threeMonthLWP = []
threeMonthFV = []
for hh in range(len(lineSlopeTimes)):
    tempInd = np.where(lineSlopeTimes[hh] == timeSubset)
    tempYear = timeSubset[tempInd[0][0]].year
    tempMonth = timeSubset[tempInd[0][0]].month
    tempDay = timeSubset[tempInd[0][0]].day

    if tempDay > 28:
        tempDay = 28
    if tempMonth > 4:
        allWaves4 = np.where((tC < DT.datetime(tempYear,tempMonth-3,tempDay)) & (tC > DT.datetime(tempYear,tempMonth-4,tempDay)))
        allWaves1 = np.where((tC < DT.datetime(tempYear,tempMonth-2,tempDay)) & (tC > DT.datetime(tempYear,tempMonth-3,tempDay)))
        allWaves2 = np.where((tC < DT.datetime(tempYear,tempMonth-1,tempDay)) & (tC > DT.datetime(tempYear,tempMonth-2,tempDay)))
        allWaves3 = np.where((tC < timeSubset[tempInd[0][0]]) & (tC > DT.datetime(tempYear,tempMonth-1,tempDay)))
    elif tempMonth == 4:
        allWaves4 = np.where((tC < DT.datetime(tempYear,tempMonth-3,tempDay)) & (tC > DT.datetime(tempYear-1,12,tempDay)))
        allWaves1 = np.where((tC < DT.datetime(tempYear,tempMonth-2,tempDay)) & (tC > DT.datetime(tempYear,tempMonth-3,tempDay)))
        allWaves2 = np.where((tC < DT.datetime(tempYear,tempMonth-1,tempDay)) & (tC > DT.datetime(tempYear,tempMonth-2,tempDay)))
        allWaves3 = np.where((tC < timeSubset[tempInd[0][0]]) & (tC > DT.datetime(tempYear,tempMonth-1,tempDay)))
    elif tempMonth == 3:
        allWaves4 = np.where((tC < DT.datetime(tempYear-1,12,tempDay)) & (tC > DT.datetime(tempYear-1,11,tempDay)))
        allWaves1 = np.where((tC < DT.datetime(tempYear-1,11,tempDay)) & (tC > DT.datetime(tempYear-1,10,tempDay)))
        allWaves2 = np.where((tC < DT.datetime(tempYear,tempMonth-1,tempDay)) & (tC > DT.datetime(tempYear,tempMonth-2,tempDay)))
        allWaves3 = np.where((tC < timeSubset[tempInd[0][0]]) & (tC > DT.datetime(tempYear,tempMonth-1,tempDay)))
    elif tempMonth == 2:
        allWaves2 = np.where((tC < DT.datetime(tempYear-1,10,tempDay)) & (tC > DT.datetime(tempYear-1,9,tempDay)))
        allWaves4 = np.where((tC < DT.datetime(tempYear-1,12,tempDay)) & (tC > DT.datetime(tempYear-1,11,tempDay)))
        allWaves1 = np.where((tC < DT.datetime(tempYear-1,11,tempDay)) & (tC > DT.datetime(tempYear-1,10,tempDay)))
        allWaves3 = np.where((tC < timeSubset[tempInd[0][0]]) & (tC > DT.datetime(tempYear,tempMonth-1,tempDay)))
    elif tempMonth == 1:
        allWaves2 = np.where((tC < DT.datetime(tempYear-1,10,tempDay)) & (tC > DT.datetime(tempYear-1,9,tempDay)))
        allWaves4 = np.where((tC < DT.datetime(tempYear-1,12,tempDay)) & (tC > DT.datetime(tempYear-1,11,tempDay)))
        allWaves1 = np.where((tC < DT.datetime(tempYear-1,11,tempDay)) & (tC > DT.datetime(tempYear-1,10,tempDay)))
        allWaves3 = np.where((tC < timeSubset[tempInd[0][0]]) & (tC > DT.datetime(tempYear-1,12,tempDay)))

    # if tempMonth > 2 and tempMonth < 10:
    #     print(lineSlopeTimes[hh])
    #     print('{}/{}/{}'.format(tempYear, tempMonth, tempDay))
    threeMonthWE.append([np.nanmean(weC[allWaves4]),np.nanmean(weC[allWaves1]),np.nanmean(weC[allWaves2]),np.nanmean(weC[allWaves3])])
    threeMonthLWP.append([np.nanmean(np.abs(lwpC[allWaves4])),np.nanmean(np.abs(lwpC[allWaves1])),np.nanmean(np.abs(lwpC[allWaves2])),np.nanmean(np.abs(lwpC[allWaves3]))])
    threeMonthFV.append([np.nanmean(fV[allWaves4]),np.nanmean(fV[allWaves1]),np.nanmean(fV[allWaves2]),np.nanmean(fV[allWaves3])])
    lineSubset.append(lineSlopes[hh])



lowSlopes = np.where((np.array(lineSubset) < 1.2))
highSlopes = np.where((np.array(lineSubset) > 1.2))


plotTime = [1,2,3,4]
ax1 = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
seasonalMean = np.mean(np.array(threeMonthWE)[highSlopes[0],:],axis=0)
seasonalStd = np.std(np.array(threeMonthWE)[highSlopes[0],:],axis=0)
ax1.plot(plotTime,seasonalMean,'o-',label='2 Bars Maintained',color='r')
ax1.fill_between(plotTime, seasonalMean - seasonalStd, seasonalMean + seasonalStd, color='r', alpha=0.2)

seasonalMean2 = np.mean(np.array(threeMonthWE)[lowSlopes[0],:],axis=0)
seasonalStd2 = np.std(np.array(threeMonthWE)[lowSlopes[0],:],axis=0)
ax1.plot(plotTime,seasonalMean2,'o-',label='New Bar Created',color='blue')
ax1.fill_between(plotTime, seasonalMean2 - seasonalStd2, seasonalMean2 + seasonalStd2, color='blue', alpha=0.2)

ax1.set_xticks([plotTime[0],plotTime[1],plotTime[2],plotTime[3]])
ax1.set_xticklabels(['4-months prior','3-months prior','2-months prior','1-month prior',])
ax1.legend()
ax1.set_title('Wave energy preceeding arrows in (a)')
ax1.text(-0.05, 1.03, 'b.', transform=ax1.transAxes, size=14, weight='bold')
ax1.set_ylabel('Wave Energy ')

import matplotlib.image as mpimg
ax2 = plt.subplot2grid((2,3),(0,1),rowspan=2,colspan=2)
ax2.imshow(mpimg.imread('conceptualFlowchart.png'))
ax2.axes.get_xaxis().set_visible(False)
ax2.axes.get_yaxis().set_visible(False)
ax2.text(-0.05, 1.02, 'c.', transform=ax2.transAxes, size=14, weight='bold')
ax2.set_title('Conceptual Flow Map for Interannual Sandbar Evolutions',fontsize=18)
ax2.spines["bottom"].set_linewidth(0)
ax2.spines["top"].set_linewidth(0)
ax2.spines["left"].set_linewidth(0)
ax2.spines["right"].set_linewidth(0)

plt.tight_layout()