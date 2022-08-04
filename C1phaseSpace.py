
import datetime as DT
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from glob import glob
import os
from scipy.interpolate import interp1d
import scipy.io
from matplotlib import gridspec
import pickle



import datetime as DT
import numpy as np
import peakutils, os
from matplotlib import pyplot as plt
from scipy import signal




dbfile = open('sandbarsSouthernTransect_referencedMHHW.pickle', 'rb')
data = pickle.load(dbfile)
dbfile.close()
alllines = data['alllines'][:,0:-1]
xinterp = data['xinterp'][0:-1]
time = data['time']#[0:-1]


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


dt = time[1:]-time[0:-1]

days = [t.days for t in dt]#[item for sublist in m for item in sublist]
days2 = np.array(days)


def P2R(radii, angles):
    return radii * np.exp(1j*angles)

def R2P(x):
    return np.abs(x), np.angle(x)


timeind = np.arange(0,len(time))
alllinesSubset = alllines[timeind, :]

eofPredog = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred2og = np.nan * np.ones((np.shape(alllinesSubset)))

for timestep in range(len(timeind)):
    mode = 0
    eofPredog[timestep, :] = RS[timestep, mode]* SS[:, mode] * np.cos(phitS[timestep, mode] - thetaS[:, mode])
    mode = 1
    eofPred2og[timestep, :] = RS[timestep, mode]* SS[:, mode] * np.cos(phitS[timestep, mode] - thetaS[:, mode])




def findSandBarAndTrough1D(xFRF, profile, plotFname=None, **kwargs):
    """ Finds multiple bars on a single profile, will also plot QA/QC plot.  The process begins by finding a trend (see
    keyword args for details) removing it from the data.  Then smooths the profile and trend and finds peaks on that
    detrended line.  It removes portions of the line that are above mean water level (Default=0), that are towards the
    offshore portion of the domain (default = 10%), and sandbars that are within x number of cells (default = 150).  It
    will return the indices of the bars and troughs of interest.

    Assumptions:
        profiles are shoreward to seaward, left to right

    Args:
        xFRF: cross-shore coordinates to describe profile values
        profile: bathymetry values (positive up)
        plotFname: if not None, will generate a QA/QC plot that (default=None)

    Keyword Args
        'waterLevel': finds only bars below this value (Default = 0)
        'trendOrder': removes trend of order (default = 3)
        'deepWaterPercentile': remove sandbars found in the last x percent of the domain (default = 0.9)
        'minSandBarSeparation':  separation in peaks to be found in cells(Default=150 Cells)
        'meanProfile': input a trend to remove data from
        'deepWaterCoordCutoff': a cutoff to say sandbars can't exist beyond this cross-shore coordinate (value must in
                xFRF coordinates)
        'verbose': if on, will print more information to screen (default = False)

    Returns:
        Peak (list): cross-shore location of sandbar crest positions
        Trough (list): cross-shore location of trough positions

    """
    #print('got this far')

    verbose = kwargs.get('verbose', False)
    waterLevel = kwargs.get('waterLevel', 0)
    polyFit = kwargs.get('trendOrder', 3)
    deepwaterPercentile = kwargs.get('deepWaterPercentile', None)
    deepwaterCoordinate = kwargs.get('deepWaterCoordCutoff', 450)
    minSandBarSeparation = kwargs.get('minSandBarSeparation', 5)
    smoothLine = kwargs.get('smoothLengthScale', 3)
    smoothBackground = kwargs.get('smoothBackLengthScale', 5)
    profileTrend = kwargs.get('profileTrend',
                              peakutils.baseline(profile, deg=polyFit))  # find the profiles general trend
    ################################################################
    # start working on data
    assert np.size(profile) == np.size(profileTrend) == np.size(xFRF), 'ProfileTrend must be same size as input profile data, and xFRF'
    #profile_smooth = sb.running_mean(profile, smoothLine)                                                 # first smooth the data to clean out noise
    #xFRF_smooth = sb.running_mean(xFRF, smoothLine)                                                       # smooth the  cross-shore coordinates
    #profileTrend_smooth = sb.running_mean(profileTrend, smoothLine)
    filterDegree = 2   # int(np.ceil(smoothLine/4)*2-1)  # always round to odd number (half of smoothline)
    smoothLine = int(np.ceil(smoothLine/2)*2+1)    # always round to odd number (up from smoothline)
    smoothBackground = int(np.ceil(smoothBackground/2)*2+1)
    profile_smooth = signal.savgol_filter(profile, smoothLine, filterDegree)
    xFRF_smooth = signal.savgol_filter(xFRF, smoothLine, filterDegree)
    profileTrend_smooth = signal.savgol_filter(profileTrend, smoothBackground, filterDegree)
    ### check profile Trend to make sure
    findPeaksOnThisLine = profile_smooth - profileTrend_smooth                                            # remove trend
    findPeaksOnThisLine[profile_smooth >= waterLevel] = np.nan                                            # find only the one's that are below the water level
    if not np.isnan(findPeaksOnThisLine).all():                                                           # if the whole thing is nans' don't do anything
        ###################################
        # 1. find sandbar first cut peaks #
        ###################################
        peakIdx = peakutils.indexes(findPeaksOnThisLine[~np.isnan(findPeaksOnThisLine)], min_dist=minSandBarSeparation)
        # peakIdx = peakutils.indexes(findPeaksOnThisLine[~np.isnan(findPeaksOnThisLine)], min_dist=minSandBarSeparation,
        #                    thres_abs=True, thres=0.25)
        peakIdx = np.argmin(np.isnan(findPeaksOnThisLine)) + peakIdx               # SHIFT back to the original baseline (with nans)
        #if deepwaterPercentile is not None:
        #    peakIdx = peakIdx[peakIdx < len(profile_smooth)*deepwaterPercentile]       # remove any peaks found oceanward of 90% of line
        #else:
        #    peakIdx = peakIdx[xFRF_smooth[peakIdx] < deepwaterCoordinate]

        peakIdx = peakIdx[::-1]  # flip peaks to move offshore to onshore
        ############################################
        # 1a. refine peaks to point of inflection  #
        ############################################
        # make sure now that each peak is actually a local maximum by looking at the slope at each peakIDX and finding
        # the nearest zero slope towards shore
        peakIdxNew, troughIdx = [], []
        shorelineIDX = np.nanargmin(np.abs(profile_smooth))  # identify shoreline on smoothed profile
        for pp, peak in enumerate(peakIdx):
            # find point most shoreward of point with slope greater than zero (add additional point to find point before)
            dElevation = np.diff(profile_smooth[shorelineIDX:peak])                   # take the derivative of the smoothed profile
            if (dElevation > -0.0).any():                                    # are any of the slopes positive
                idxDiff = np.argwhere(dElevation > -0.0).squeeze()   # find all values that have positive slope
                idxMax = np.max(idxDiff) + shorelineIDX              # find max cross-shore location add shoreline location
                peakIdxNew.append(idxMax + 1)                     # add one to find before point of inflection
        peakIdxNew = np.unique(peakIdxNew)
        #######################
        # 2. now find troughs #
        #######################
        for peak in peakIdxNew:
            #if profile_smooth[np.argmin(profile_smooth[:peak]).squeeze()] <  profile_smooth[peak]:      # check that its shallower than the sand bar
            troughIdx.append(np.argmin(profile_smooth[:peak.squeeze()]))
        ########################################################################################
        # 3. check to see if peaks are really peaks, find local maximum between bar and trough #
        ########################################################################################
        # if np.size(peakIdx) == np.size(troughIdx):                     # we have same number of peaks and troughs
        #     for pp, peak in enumerate(peakIdx):
        #         peakIdx[pp] = troughIdx[pp] + np.argmax(profile_smooth[troughIdx[pp]:peak])
        # else:                                                          # we found a peak, but no troughs to the sand bar
        #     # now reiterate the find method, but start from found peak move towards shore
        #     for peak in peakIdx:
        #         findPeaksOnThisLine[range(len(findPeaksOnThisLine)) > peak] = np.nan
        #####################################
        # Last: plot for QA/QC              #
        #####################################
        plt.ioff()                           # turn off plot visible
        if plotFname is not None:            # now plot if interested
            if verbose: print("plotting {}".format(plotFname))
            plt.figure(figsize=(8, 5))
            try:
                plt.suptitle(DT.datetime.strptime(plotFname.split('_')[-1].split('.')[0], '%Y%m%d'))
            except ValueError: # happens when there's not a date in the filename
                plt.suptitle(os.path.basename(plotFname))
            plt.plot(xFRF, profile, 'C1.', label='Raw')
            plt.plot(xFRF_smooth, profile_smooth, 'c.', ms=2, label='smoothed')
            plt.plot(xFRF_smooth, findPeaksOnThisLine, label='Find peaks on this line')
            plt.plot(xFRF_smooth, profileTrend_smooth, label='Trend')
            plt.plot([0, len(profile)], [0,0], 'k--')
            if np.size(peakIdx) > 0:
                plt.plot(xFRF_smooth[peakIdx], profile_smooth[peakIdx], 'r.', ms=5, label='Inital bar Location')
                plt.plot(xFRF_smooth[peakIdx], findPeaksOnThisLine[peakIdx], 'r.', ms=5)
            if np.size(peakIdxNew) >0:
                plt.plot(xFRF_smooth[peakIdxNew], profile_smooth[peakIdxNew], 'rd', ms=7, label='Refined Bar location')
            if np.size(troughIdx) >0:
                plt.plot(xFRF_smooth[troughIdx], profile_smooth[troughIdx], 'bo', ms=6, label='Trough Location')
            plt.legend(loc='upper right')
            plt.ylabel('Elevation NAVD88 [m]')
            plt.xlabel('cross-shore location [m]')
            plt.savefig(plotFname)
            plt.close()
        if np.size(peakIdxNew) > 0 and np.size(troughIdx) > 0:
            return xFRF_smooth[peakIdxNew], xFRF_smooth[troughIdx]
        elif np.size(peakIdxNew) > 0 and np.size(troughIdx) == 0:
            return xFRF_smooth[peakIdxNew], None
        elif np.size(peakIdxNew) == 0 and np.size(troughIdx) > 0:
            return None, xFRF_smooth[troughIdx]
        else:
            return None, None
    else:
        return None, None





thetasMode1 = np.arange(0,1200,5)
thetasMode2 = np.arange(0,1200,5)

eofHypo = np.nan * np.ones((len(thetasMode1),len(xinterpS)))
eofHypo2 = np.nan * np.ones((len(thetasMode1),len(xinterpS)))

for tt in range(len(thetasMode1)):
   eofHypo[tt,:] = SS[:, 0] * np.nanmean(RS[:, 0]) * np.cos(thetasMode1[tt] * np.pi / 180 - thetaS[:, 0]+180)
   eofHypo2[tt,:] = SS[:, 1] * np.nanmean(RS[:, 1]) * np.cos(thetasMode2[tt] * np.pi / 180 - thetaS[:, 1])


thetag, xintg = np.meshgrid(thetasMode1, xinterpS)









# thetasMode1 = np.arange(-180,180,2)
# thetasMode2 = np.arange(-100,260,2)

thetasMode1 = np.arange(0,360,2)
thetasMode2 = np.arange(-45,315,2)
innerBarX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
outerBarX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerBarZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
outerBarZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterBarX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterBarZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))

outerTroughZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
outerTroughX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerTroughZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerTroughX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterTroughZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterTroughX = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))

innerBarDZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
outerBarDZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))
innerOuterBarDZ = np.nan*np.ones((len(thetasMode2),len(thetasMode1)))

# innerBar = np.nan * np.ones((alllines.shape[0],))
# outerBar = np.nan * np.ones((alllines.shape[0],))
#
# zOuterBar = np.nan * np.ones((alllines.shape[0],))
# zInnerBar = np.nan * np.ones((alllines.shape[0],))
#
# innerTrough = np.nan * np.ones((alllines.shape[0],))
# outerTrough = np.nan * np.ones((alllines.shape[0],))
# zOuterTrough = np.nan * np.ones((alllines.shape[0],))
# zInnerTrough = np.nan * np.ones((alllines.shape[0],))
#

for tt in range(len(thetasMode1)):
   eofPred = SS[:, 0] * np.nanmean(RS[:, 0]) * np.cos(thetasMode1[tt] * np.pi / 180 - thetaS[:, 0])
   for yy in range(len(thetasMode2)):
      eofPred2 = SS[:, 1] * np.nanmean(RS[:, 1]) * np.cos(thetasMode2[yy] * np.pi / 180 - thetaS[:, 1])
      combined = np.real(np.mean(alllinesS,axis=0)+eofPred+eofPred2)
      fname = "/home/dylananderson/projects/duckGeomorph/sandBarTool/Mode1_{}_Mode2_{}.png".format(thetasMode1[tt],thetasMode2[yy])
      xFRFbar, xFRFtrough = findSandBarAndTrough1D(xinterpS, combined, plotFname=None, smoothLengthScale=5, profileTrend=np.mean(alllines,axis=0))
      if xFRFbar is not None:
            #for sandbarX in xFRFbar:
            #    ax[1].plot(sandbarX, time[tt], 'ro', label='bar')

         if len(xFRFbar) > 1:
            outerBarX[yy,tt] = np.max(xFRFbar)
            zOuterInd = np.where((np.round(2*np.max(xFRFbar))/2 == xinterpS))
            outerBarZ[yy,tt] = combined[zOuterInd[0]]
            outerTroughZ[yy,tt] = np.max(xFRFtrough)
            outerTroughZInd = np.where((np.round(2*np.max(xFRFtrough))/2 == xinterpS))
            outerTroughZ[yy,tt] = combined[outerTroughZInd[0]]

            innerBarX[yy,tt] = np.min(xFRFbar)
            zInnerInd = np.where((np.round(2*np.min(xFRFbar))/2 == xinterpS))
            innerBarZ[yy,tt] = combined[zInnerInd[0]]
            innerTroughZ[yy,tt] = np.max(xFRFtrough)
            innerTroughZInd = np.where((np.round(2*np.min(xFRFtrough))/2 == xinterpS))
            innerTroughZ[yy,tt] = combined[innerTroughZInd[0]]

            outerBarDZ[yy,tt] = outerBarZ[yy,tt]-outerTroughZ[yy,tt]
            innerBarDZ[yy,tt] = innerBarZ[yy,tt]-innerTroughZ[yy,tt]

            # innerOuterBarX[yy,tt] = np.max(xFRFbar)
            # zOuterInd = np.where((np.round(2*np.max(xFRFbar))/2 == xinterpS))
            # innerOuterBarZ[yy,tt] = combined[zOuterInd[0]]

         else:
            # innerBarX[yy,tt] = xFRFbar
            # zInnerInd = np.where((np.round(2*xFRFbar)/2 == xinterpS))
            # innerBarZ[yy,tt] = combined[zInnerInd[0]]
            #
            # outerBarX[yy,tt] = xFRFbar
            # zOuterInd = np.where((np.round(2*xFRFbar)/2 == xinterpS))
            # outerBarZ[yy,tt] = combined[zOuterInd[0]]

            innerOuterBarX[yy,tt] = xFRFbar
            zOuterInd = np.where((np.round(2*np.max(xFRFbar))/2 == xinterpS))
            innerOuterBarZ[yy,tt] = combined[zOuterInd[0]]
            innerOuterTroughZ[yy, tt] = xFRFtrough
            innerOuterTroughZInd = np.where((np.round(2 * np.max(xFRFtrough)) / 2 == xinterpS))
            innerOuterTroughZ[yy, tt] = combined[innerOuterTroughZInd[0]]
            innerOuterBarDZ[yy,tt] = innerOuterBarZ[yy,tt]-innerOuterTroughZ[yy,tt]


      # if xFRFtrough is not None:
        #     #for troughX in xFRFtrough:
        #     #    ax[1].plot(troughX, time[tt], 'bd', label='trough')
        #     if len(xFRFtrough) > 1:
        #         outerTrough[tt] = np.max(xFRFtrough)
        #         zOuterInd = np.where((np.round(2*np.max(xFRFtrough))/2 == bathyX))
        #         zOuterTrough[tt] = bathy[zOuterInd[0]]
        #
        #         innerTrough[tt] = np.min(xFRFtrough)
        #         zInnerInd = np.where((np.round(2*np.min(xFRFtrough))/2 == bathyX))
        #         zInnerTrough[tt] = bathy[zInnerInd[0]]
        #     else:
        #         outerTrough[tt] = xFRFtrough
        #         zOuterInd = np.where((np.round(2*xFRFtrough)/2 == bathyX))
        #         zOuterTrough[tt] = bathy[zOuterInd[0]]
        #barPatch = mpatches.Patch(color='red', label='sandbar')
        #troughPatch = mpatches.Patch(color='blue', label='trough')

meshTheta1, meshTheta2 = np.meshgrid(thetasMode1,thetasMode2)

# fig2 = plt.figure(figsize=(12,10))
# ax1 = plt.subplot2grid((2,2),(0,0),rowspan=1,colspan=1)
# p1 = ax1.pcolor(meshTheta1,meshTheta2,innerBarX,vmin=30,vmax=300)
# cb1 = plt.colorbar(p1,ax=ax1)
# cb1.set_label('xFRF')
# ax2 = plt.subplot2grid((2,2),(0,1),rowspan=1,colspan=1)
# p2 = ax2.pcolor(meshTheta1,meshTheta2,outerBarX,vmin=30,vmax=300)
# ax2.set_xlabel('Mode 1 Phase')
# ax1.set_xlabel('Mode 1 Phase')
# ax1.set_ylabel('Mode 2 Phase')
# ax1.set_title('Inner Bar Cross-shore Location')
# ax2.set_title('Outer Bar Cross-shore Location')
# cb2 = plt.colorbar(p2,ax=ax2)
# cb2.set_label('xFRF')
#
#
# ax3 = plt.subplot2grid((2,2),(1,0),rowspan=1,colspan=1)
# p3 = ax3.pcolor(meshTheta1,meshTheta2,innerBarDZ)#,vmin=-6,vmax=-1.5)
# cb3 = plt.colorbar(p3,ax=ax3)
# cb3.set_label('Depth (m)')
#
# ax4 = plt.subplot2grid((2,2),(1,1),rowspan=1,colspan=1)
# p4 = ax4.pcolor(meshTheta1,meshTheta2,outerBarDZ)#,vmin=-6,vmax=-1.5)
# cb4 = plt.colorbar(p4,ax=ax4)
# cb4.set_label('Depth (m)')
#
# ax4.set_xlabel('Mode 1 Phase')
# ax3.set_xlabel('Mode 1 Phase')
# ax3.set_ylabel('Mode 2 Phase')
# ax3.set_title('Inner Bar Depth')
# ax4.set_title('Outer Bar Depth')
# plt.tight_layout()
#
# plt.show()

# index1 = np.where((timeS > DT.datetime(1984,9,1)) & (timeS < DT.datetime(1988,9,1)))
# indexedTime = timeS[index1]
# indexedOrdinal = [i.toordinal() for i in indexedTime]
# plt.figure()
# mode1prop = np.unwrap(phitS[index1[0],0])
# mode2prop = np.unwrap(phitS[index1[0],1])
# plt.plot(mode1prop,mode2prop)




fig2 = plt.figure(figsize=(14,8))
ax1 = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1)
p1 = ax1.pcolor(meshTheta1,meshTheta2,innerBarX,vmin=30,vmax=300)
cb1 = plt.colorbar(p1,ax=ax1)
cb1.set_label('xFRF')
ax2 = plt.subplot2grid((2,3),(0,1),rowspan=1,colspan=1)
p2 = ax2.pcolor(meshTheta1,meshTheta2,outerBarX,vmin=30,vmax=300)
ax2.set_xlabel('Mode 1 Phase')
ax1.set_xlabel('Mode 1 Phase')
# ax3.set_xlabel('Mode 1 Phase')
ax1.set_ylabel('Mode 2 Phase')
ax1.set_title('Inner Bar Cross-shore Location')
ax2.set_title('Outer Bar Cross-shore Location')

cb2 = plt.colorbar(p2,ax=ax2)
cb2.set_label('xFRF')
ax1c = plt.subplot2grid((2,3),(0,2),rowspan=1,colspan=1)
p1c = ax1c.pcolor(meshTheta1,meshTheta2,innerOuterBarX,vmin=30,vmax=300)
cb1c = plt.colorbar(p1c,ax=ax1c)
cb1c.set_label('xFRF')
ax1c.set_title('Single Bar Cross-shore Location')

ax3 = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
p3 = ax3.pcolor(meshTheta1,meshTheta2,innerBarDZ,vmin=0,vmax=1)
cb3 = plt.colorbar(p3,ax=ax3)
cb3.set_label('Depth (m)')

ax4 = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
p4 = ax4.pcolor(meshTheta1,meshTheta2,outerBarDZ,vmin=0,vmax=1)
cb4 = plt.colorbar(p4,ax=ax4)
cb4.set_label('Depth (m)')

ax4c = plt.subplot2grid((2,3),(1,2),rowspan=1,colspan=1)
p4c = ax4c.pcolor(meshTheta1,meshTheta2,innerOuterBarDZ,vmin=0,vmax=1)
cb4c = plt.colorbar(p4c,ax=ax4c)
cb4c.set_label('Depth (m)')

ax4c.set_xlabel('Mode 1 Phase')
ax4.set_xlabel('Mode 1 Phase')
ax3.set_xlabel('Mode 1 Phase')
ax3.set_ylabel('Mode 2 Phase')
ax3.set_title('Inner Bar Magnitude')
ax4.set_title('Outer Bar Magnitude')
ax4c.set_title('Single Bar Magnitude')

plt.tight_layout()





plt.figure(figsize=(14,10))
ax101 = plt.subplot2grid((1,2),(0,0),rowspan=1,colspan=1)
p101 = ax101.pcolor(meshTheta1,meshTheta2,innerBarZ-outerBarZ)#vmin=0,vmax=1)
cb101 = plt.colorbar(p101,ax=ax101)
ax101.set_xlabel(r'Mode 2 Phase (offshore $\longrightarrow$)')
ax101.set_ylabel(r'Mode 2 Phase (onshore $\longrightarrow$)')
cb101.set_label('Difference in Bar Crest Depth')
ax101.set_title('Onshore Propagation Separates the Crest''s Depths')
# ax101.set_xlim([-180,180])
# ax101.set_ylim([-100,260])

ax102 = plt.subplot2grid((1,2),(0,1),rowspan=1,colspan=1)
p102 = ax102.pcolor(meshTheta1,meshTheta2,outerBarX-innerBarX)#vmin=0,vmax=1)
cb102 = plt.colorbar(p102,ax=ax102)
ax102.set_xlabel(r'Mode 2 Phase (offshore $\longrightarrow$)')
ax102.set_ylabel(r'Mode 2 Phase (onshore $\longrightarrow$)')
cb102.set_label('Difference in Bar Crest Cross-shore Location')
ax102.set_title('Offshore Propagation Converges the Crest''s Horizontal Distance')
# ax102.set_xlim([-180,180])
# ax102.set_ylim([-100,260])

plt.show()



import copy
doublebar = copy.deepcopy(outerBarX)
doubleIndex = np.where((doublebar > 0))
doublebar[doubleIndex] = 0.5

singlebar = copy.deepcopy(innerOuterBarX)
singleIndex = np.where((singlebar > 0))
doublebar[singleIndex] = 0.25









# plt.style.use('dark_background')
fig = plt.figure(figsize=(10,11))




ax1 = plt.subplot2grid((5,3),(3,0),rowspan=1,colspan=1)
p1 = ax1.pcolor(meshTheta1,meshTheta2,innerBarX,vmin=30,vmax=300,cmap='plasma')
# ax1.set_xlim([-180,180])
# ax1.set_ylim([-100,260])

# cb1 = plt.colorbar(p1,ax=ax1)
# cb1.set_label('xFRF')
# ax1.text(50,200,'Location')

ax2 = plt.subplot2grid((5,3),(3,1),rowspan=1,colspan=1)
p2 = ax2.pcolor(meshTheta1,meshTheta2,outerBarX,vmin=30,vmax=300,cmap='plasma')
# ax2.set_xlim([-180,180])
# ax2.set_ylim([-100,260])

# ax2.text(50,200,'Location')

# ax2.set_xlabel('Mode 1 Phase')
# ax1.set_xlabel('Mode 1 Phase')
# ax3.set_xlabel('Mode 1 Phase')
# ax1.set_ylabel('Mode 2 Phase')
ax1.set_ylabel(r'Onshore (deg) $\Longrightarrow$')

ax1.set_title('Inner Bar',fontsize=12)
ax2.set_title('Outer Bar',fontsize=12)

# cb2 = plt.colorbar(p2,ax=ax2)
# cb2.set_label('xFRF')
ax1c = plt.subplot2grid((5,3),(3,2),rowspan=1,colspan=1)
p1c = ax1c.pcolor(meshTheta1,meshTheta2,innerOuterBarX,vmin=30,vmax=300,cmap='plasma')
# ax1c.set_xlim([-180,180])
# ax1c.set_ylim([-100,260])


# cb1c = plt.colorbar(p1c,ax=ax1c)
# cb1c.set_label('xFRF')
ax1c.set_title('Single Bar',fontsize=12)

ax3 = plt.subplot2grid((5,3),(4,0),rowspan=1,colspan=1)
p3 = ax3.pcolor(meshTheta1,meshTheta2,innerBarDZ,vmin=0,vmax=1)
# ax3.set_xlim([-180,180])
# ax3.set_ylim([-100,260])

# cb3 = plt.colorbar(p3,ax=ax3)
# cb3.set_label('Depth (m)')
# ax3.text(50,200,'Magnitude')

ax4 = plt.subplot2grid((5,3),(4,1),rowspan=1,colspan=1)
p4 = ax4.pcolor(meshTheta1,meshTheta2,outerBarDZ,vmin=0,vmax=1)
# ax4.set_xlim([-180,180])
# ax4.set_ylim([-100,260])

# cb4 = plt.colorbar(p4,ax=ax4)
# cb4.set_label('Depth (m)')
# ax4.text(50,200,'Magnitude')

ax4c = plt.subplot2grid((5,3),(4,2),rowspan=1,colspan=1)
p4c = ax4c.pcolor(meshTheta1,meshTheta2,innerOuterBarDZ,vmin=0,vmax=1)
# ax4c.set_xlim([-180,180])
# ax4c.set_ylim([-100,260])

ax4c.set_xlabel(r'Offshore (deg) $\Longrightarrow$')
ax4.set_xlabel(r'Offshore (deg) $\Longrightarrow$')
ax3.set_xlabel(r'Offshore (deg) $\Longrightarrow$')
ax3.set_ylabel(r'Onshore (deg) $\Longrightarrow$')

plt.subplots_adjust(right=0.87)
cbar_ax = fig.add_axes([0.91, 0.28, 0.02, 0.13])
cb1c = plt.colorbar(p1c,cax=cbar_ax)
cb1c.set_label('Dist. from shore (m)')

cbar_ax2 = fig.add_axes([0.91, 0.11, 0.02, 0.13])
cb1c2 = plt.colorbar(p4c,cax=cbar_ax2)
cb1c2.set_label('Bar Magnitude (m)')


conceptbar_ax = fig.add_axes([0.55, 0.55, 0.4, 0.4])
#conceptbar_ax = plt.subplot2grid((5,3),(0,1),rowspan=3,colspan=2)
con = conceptbar_ax.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
conceptbar_ax.set_xlabel(r'Offshore (deg) $\Longrightarrow$')
conceptbar_ax.set_ylabel(r'Onshore (deg) $\Longrightarrow$')
conceptbar_ax.set_title('Morphologic Trends')


conceptbar_ax.arrow(x=250,y=-20,dx=0,dy=50,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.text(262,20,'Shallow Single')
conceptbar_ax.text(262,0, 'Bar Growth')

conceptbar_ax.arrow(x=80,y=100,dx=0,dy=135,width=4,facecolor='orange',edgecolor='orange')

conceptbar_ax.arrow(x=110,y=70,dx=100,dy=0,width=4,facecolor='blue',edgecolor='blue')
conceptbar_ax.text(120,80,'Inner Bar Growth')
conceptbar_ax.text(120,50, 'Outer Bar Decay')

conceptbar_ax.arrow(x=110,y=220,dx=100,dy=0,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.text(120,250,'Double Bar to')
conceptbar_ax.text(120,230, 'Double Terrace')

conceptbar_ax.arrow(x=80,y=0,dx=0,dy=40,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.arrow(x=80,y=0,dx=40,dy=0,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.text(90,30,'New Inner Bar')
conceptbar_ax.text(90,10, 'Formation')

conceptbar_ax.arrow(x=200,y=100,dx=0,dy=40,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.arrow(x=200,y=100,dx=40,dy=0,width=4,facecolor='green',edgecolor='green')
conceptbar_ax.text(210,130,'Outer Bar')
conceptbar_ax.text(210,110, 'Death')

conceptbar_ax.arrow(x=290,y=60,dx=0,dy=40,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.arrow(x=290,y=60,dx=40,dy=0,width=4,facecolor='black',edgecolor='black')

conceptbar_ax.arrow(x=250,y=230,dx=0,dy=40,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.arrow(x=250,y=230,dx=40,dy=0,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.text(262,265,'Single Bar from')
conceptbar_ax.text(262,245, 'Inner Terrace')

# conceptbar_ax.set_xlim([-180,180])
# conceptbar_ax.set_ylim([-100,260])
# plt.tight_layout()
plt.show()