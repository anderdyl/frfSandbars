#   Lets identify a transect to the north to check on how much nourishment evolution was caught

import datetime as DT
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import os
from scipy.interpolate import interp1d



geomorphdir = '/media/dylananderson/Elements/filteredFRF_GeomorphUpdate/'


files = os.listdir(geomorphdir)

files.sort()

files_path = [os.path.abspath(geomorphdir) for x in os.listdir(geomorphdir)]

def getBathy(file, lower, upper):
    bathy = Dataset(file)

    xs_bathy = bathy.variables['xFRF'][:]
    ys_bathy = bathy.variables['yFRF'][:]
    zs_bathy = bathy.variables['elevation'][:]
    ts_bathy = bathy.variables['time'][:]
    pr_bathy = bathy.variables['profileNumber'][:]

    zs_bathy = np.ma.masked_where((pr_bathy > upper), zs_bathy)
    ys_bathy = np.ma.masked_where((pr_bathy > upper), ys_bathy)
    xs_bathy = np.ma.masked_where((pr_bathy > upper), xs_bathy)
    pr_bathy = np.ma.masked_where((pr_bathy > upper), pr_bathy)
    ts_bathy = np.ma.masked_where((pr_bathy > upper), ts_bathy)

    zs_bathy = np.ma.masked_where((pr_bathy < lower), zs_bathy)
    ys_bathy = np.ma.masked_where((pr_bathy < lower), ys_bathy)
    xs_bathy = np.ma.masked_where((pr_bathy < lower), xs_bathy)
    pr_bathy = np.ma.masked_where((pr_bathy < lower), pr_bathy)
    ts_bathy = np.ma.masked_where((pr_bathy < lower), ts_bathy)

    output = dict()
    output['x'] = xs_bathy
    output['y'] = ys_bathy
    output['z'] = zs_bathy
    output['pr'] = pr_bathy
    output['t'] = ts_bathy

    return output


def interpBathy(xSub, zSub, x):

    f = interp1d(xSub, zSub, kind='linear', bounds_error=False)
    newz = f(x)
    newBathy = dict()
    newBathy['x'] = x
    newBathy['z'] = newz
    return newBathy

def interpBathy2(xSub,zSub,x):
    from Loess import Loess

    loess = Loess(xSub, zSub)
    zloess = np.nan * np.ones(len(x),)
    count = 0
    for xxx in x:
        zloess[count] = loess.estimate(xxx, window=5, use_matrix=False, degree=1)
        count = count+1
    newBathy = dict()
    newBathy['x'] = x
    newBathy['z'] = zloess
    return newBathy


subset = files.copy()

colormap = plt.cm.gist_ncar


labels = []

xinterp = np.arange(90, 650, 5)
higherResX = np.arange(90, 650, 1)
xinterpRef = np.arange(0,515,5)
referenceHeight = 0.45
deeperThan = 580
surveyType = dict()
noneOnshore = np.empty((len(xinterp),))
noneOffshore = np.empty((len(xinterp),))

count = 0
count2 = 0
count3 = 0
count4 = 0
count5 = 0
worstcount = 0
badCounter = 0
badCounterShort = 0
badCounterLong = 0
#fig = plt.figure(figsize=(10,10))
for i in range(len(subset)):

    file_params = subset[i].split('_')

    ## ### Southern Lines
    data1 = getBathy(os.path.join(geomorphdir, subset[i]), lower=-10, upper=20)
    data2 = getBathy(os.path.join(geomorphdir, subset[i]), lower=20, upper=60)
    data3 = getBathy(os.path.join(geomorphdir, subset[i]), lower=60, upper=100)
    data4 = getBathy(os.path.join(geomorphdir, subset[i]), lower=65430, upper=65460)
    data5 = getBathy(os.path.join(geomorphdir, subset[i]), lower=65480, upper=65520)

    temp = subset[i].split('_')
    if temp[1] == 'geomorphology':
        temp2 = temp[-1].split('.')
        surveydate = DT.datetime.strptime(temp2[0], '%Y%m%d')
    else:
        surveydate = DT.datetime.strptime(temp[1], '%Y%m%d')
    elevs = data1['z']
    cross = data1['x']
    crossind = np.argsort(data1['x'])
    crossS = cross[crossind]
    elevsS = elevs[crossind]

    elevs2 = data2['z']
    cross2 = data2['x']
    crossind2 = np.argsort(data2['x'])
    crossS2 = cross2[crossind2]
    elevsS2 = elevs2[crossind2]

    elevs3 = data3['z']
    cross3 = data3['x']
    crossind3 = np.argsort(data3['x'])
    crossS3 = cross3[crossind3]
    elevsS3 = elevs3[crossind3]

    elevs4 = data4['z']
    cross4 = data4['x']
    crossind4 = np.argsort(data4['x'])
    crossS4 = cross4[crossind4]
    elevsS4 = elevs4[crossind4]

    elevs5 = data5['z']
    cross5 = data5['x']
    crossind5 = np.argsort(data5['x'])
    crossS5 = cross5[crossind5]
    elevsS5 = elevs5[crossind5]

    xSub = np.ma.MaskedArray.filled(crossS, np.nan)
    zSub = np.ma.MaskedArray.filled(elevsS, np.nan)
    xSub2 = np.ma.MaskedArray.filled(crossS2, np.nan)
    zSub2 = np.ma.MaskedArray.filled(elevsS2, np.nan)
    xSub3 = np.ma.MaskedArray.filled(crossS3, np.nan)
    zSub3 = np.ma.MaskedArray.filled(elevsS3, np.nan)
    xSub4 = np.ma.MaskedArray.filled(crossS4, np.nan)
    zSub4 = np.ma.MaskedArray.filled(elevsS4, np.nan)
    xSub5 = np.ma.MaskedArray.filled(crossS5, np.nan)
    zSub5 = np.ma.MaskedArray.filled(elevsS5, np.nan)

    realValues = ~np.isnan(xSub)
    xSubNew = xSub[~np.isnan(xSub)]
    zSubNew = zSub[~np.isnan(xSub)]

    realValues2 = ~np.isnan(xSub2)
    xSubNew2 = xSub2[~np.isnan(xSub2)]
    zSubNew2 = zSub2[~np.isnan(xSub2)]

    realValues3 = ~np.isnan(xSub3)
    xSubNew3 = xSub3[~np.isnan(xSub3)]
    zSubNew3 = zSub3[~np.isnan(xSub3)]

    realValues4 = ~np.isnan(xSub4)
    xSubNew4 = xSub4[~np.isnan(xSub4)]
    zSubNew4 = zSub4[~np.isnan(xSub4)]

    realValues5 = ~np.isnan(xSub5)
    xSubNew5 = xSub5[~np.isnan(xSub5)]
    zSubNew5 = zSub5[~np.isnan(xSub5)]

    temp = np.hstack((np.hstack((np.hstack((realValues,realValues2)),realValues3)),realValues4))
    tempRealValues = np.hstack((temp,realValues5))
    tempXSub = np.hstack((np.hstack((np.hstack((xSubNew,xSubNew2)),xSubNew3)),xSubNew4))
    tempXSubNew = np.hstack((tempXSub,xSubNew5))


    if any(tempRealValues):
        if np.nanmax(tempXSubNew) > deeperThan:
            if np.nanmin(tempXSubNew) < 115:
                if len(tempXSubNew) > 6:

                    if any(realValues):
                        newdata = interpBathy(xSubNew, zSubNew, xinterp)
                    else:
                        newdata = []
                    if any(realValues2):
                        newdata2 = interpBathy(xSubNew2,zSubNew2, xinterp)
                    else:
                        newdata2 = []
                    if any(realValues3):
                        newdata3 = interpBathy(xSubNew3, zSubNew3, xinterp)
                    else:
                        newdata3 = []
                    if any(realValues4):
                        newdata4 = interpBathy(xSubNew4, zSubNew4, xinterp)
                    else:
                        newdata4 = []
                    if any(realValues5):
                        newdata5 = interpBathy(xSubNew5, zSubNew5, xinterp)
                    else:
                        newdata5 = []

                    if count == 0:
                        # Do we have data at all 5 transects?
                        if any(realValues) & any(realValues2) & any(realValues3) & any(realValues4) & any(realValues5):
                            avgProf = np.nanmean(
                                (newdata['z'], newdata2['z'], newdata3['z'], newdata4['z'], newdata5['z']), axis=0)
                            diff1 = np.abs(newdata['z'] - avgProf)
                            diff2 = np.abs(newdata2['z'] - avgProf)
                            diff3 = np.abs(newdata3['z'] - avgProf)
                            diff4 = np.abs(newdata4['z'] - avgProf)
                            diff5 = np.abs(newdata5['z'] - avgProf)
                            tempDiff = np.nanmax(
                                np.vstack((np.vstack((np.vstack((np.vstack((diff1, diff2)), diff3)), diff4)), diff5)),
                                axis=0)
                        else:
                            # What about data at only 4 of the 5?
                            if any(realValues) & any(realValues2) & any(realValues3) & any(realValues4):
                                avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata3['z'], newdata4['z']),
                                                     axis=0)
                                diff1 = np.abs(newdata['z'] - avgProf)
                                diff2 = np.abs(newdata2['z'] - avgProf)
                                diff3 = np.abs(newdata3['z'] - avgProf)
                                diff4 = np.abs(newdata4['z'] - avgProf)
                                tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((diff1, diff2)), diff3)), diff4)),
                                                     axis=0)

                            elif any(realValues) & any(realValues3) & any(realValues4) & any(realValues5):
                                avgProf = np.nanmean((newdata['z'], newdata3['z'], newdata4['z'], newdata5['z']),
                                                     axis=0)
                                diff1 = np.abs(newdata['z'] - avgProf)
                                diff5 = np.abs(newdata5['z'] - avgProf)
                                diff3 = np.abs(newdata3['z'] - avgProf)
                                diff4 = np.abs(newdata4['z'] - avgProf)
                                tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((diff1, diff5)), diff3)), diff4)),
                                                     axis=0)
                            elif any(realValues) & any(realValues2) & any(realValues4) & any(realValues5):
                                avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata4['z'], newdata5['z']),
                                                     axis=0)
                                diff1 = np.abs(newdata['z'] - avgProf)
                                diff2 = np.abs(newdata2['z'] - avgProf)
                                diff5 = np.abs(newdata5['z'] - avgProf)
                                diff4 = np.abs(newdata4['z'] - avgProf)
                                tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((diff1, diff2)), diff5)), diff4)),
                                                     axis=0)
                            elif any(realValues) & any(realValues2) & any(realValues3) & any(realValues5):
                                avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata3['z'], newdata5['z']),
                                                     axis=0)
                                diff1 = np.abs(newdata['z'] - avgProf)
                                diff2 = np.abs(newdata2['z'] - avgProf)
                                diff3 = np.abs(newdata3['z'] - avgProf)
                                diff5 = np.abs(newdata5['z'] - avgProf)
                                tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((diff1, diff2)), diff3)), diff5)),
                                                     axis=0)
                            elif any(realValues2) & any(realValues3) & any(realValues4) & any(realValues5):
                                avgProf = np.nanmean((newdata2['z'], newdata3['z'], newdata4['z'], newdata5['z']),
                                                     axis=0)
                                diff5 = np.abs(newdata5['z'] - avgProf)
                                diff2 = np.abs(newdata2['z'] - avgProf)
                                diff3 = np.abs(newdata3['z'] - avgProf)
                                diff4 = np.abs(newdata4['z'] - avgProf)
                                tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((diff5, diff2)), diff3)), diff4)),
                                                     axis=0)
                            else:
                                # What about if data is at only 3 of the 5?
                                if any(realValues) & any(realValues2) & any(realValues3):
                                    avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata3['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff2)), diff3)), axis=0)
                                elif any(realValues) & any(realValues2) & any(realValues4):
                                    avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata4['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff2)), diff4)), axis=0)
                                elif any(realValues) & any(realValues2) & any(realValues5):
                                    avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata5['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff5 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff2)), diff5)), axis=0)
                                elif any(realValues) & any(realValues3) & any(realValues4):
                                    avgProf = np.nanmean((newdata['z'], newdata3['z'], newdata4['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff4)), diff3)), axis=0)
                                elif any(realValues) & any(realValues3) & any(realValues5):
                                    avgProf = np.nanmean((newdata['z'], newdata3['z'], newdata5['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff5 = np.abs(newdata5['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff5)), diff3)), axis=0)
                                elif any(realValues) & any(realValues4) & any(realValues5):
                                    avgProf = np.nanmean((newdata['z'], newdata4['z'], newdata5['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    diff5 = np.abs(newdata5['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff4)), diff5)), axis=0)
                                elif any(realValues2) & any(realValues3) & any(realValues4):
                                    avgProf = np.nanmean((newdata2['z'], newdata3['z'], newdata4['z']), axis=0)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff4, diff2)), diff3)), axis=0)
                                elif any(realValues2) & any(realValues3) & any(realValues5):
                                    avgProf = np.nanmean((newdata2['z'], newdata3['z'], newdata5['z']), axis=0)
                                    diff5 = np.abs(newdata5['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff5, diff2)), diff3)), axis=0)
                                elif any(realValues2) & any(realValues4) & any(realValues5):
                                    avgProf = np.nanmean((newdata2['z'], newdata4['z'], newdata5['z']), axis=0)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff5 = np.abs(newdata5['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff4, diff2)), diff5)), axis=0)
                                elif any(realValues3) & any(realValues4) & any(realValues5):
                                    avgProf = np.nanmean((newdata3['z'], newdata4['z'], newdata5['z']), axis=0)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    diff5 = np.abs(newdata5['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff4, diff5)), diff3)), axis=0)


                                else:
                                    if any(realValues) & any(realValues2):
                                        avgProf = np.nanmean((newdata['z'], newdata2['z']), axis=0)
                                        diff1 = np.abs(newdata['z'] - avgProf)
                                        diff2 = np.abs(newdata2['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff1, diff2)), axis=0)
                                    elif any(realValues) & any(realValues3):
                                        avgProf = np.nanmean((newdata['z'], newdata3['z']), axis=0)
                                        diff1 = np.abs(newdata['z'] - avgProf)
                                        diff3 = np.abs(newdata3['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff1, diff3)), axis=0)
                                    elif any(realValues) & any(realValues4):
                                        avgProf = np.nanmean((newdata['z'], newdata4['z']), axis=0)
                                        diff1 = np.abs(newdata['z'] - avgProf)
                                        diff4 = np.abs(newdata4['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff1, diff4)), axis=0)
                                    elif any(realValues) & any(realValues5):
                                        avgProf = np.nanmean((newdata['z'], newdata5['z']), axis=0)
                                        diff1 = np.abs(newdata['z'] - avgProf)
                                        diff5 = np.abs(newdata5['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff1, diff5)), axis=0)
                                    elif any(realValues2) & any(realValues5):
                                        avgProf = np.nanmean((newdata2['z'], newdata5['z']), axis=0)
                                        diff5 = np.abs(newdata5['z'] - avgProf)
                                        diff2 = np.abs(newdata2['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff5, diff2)), axis=0)
                                    elif any(realValues2) & any(realValues4):
                                        avgProf = np.nanmean((newdata2['z'], newdata4['z']), axis=0)
                                        diff4 = np.abs(newdata4['z'] - avgProf)
                                        diff2 = np.abs(newdata2['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff4, diff2)), axis=0)
                                    elif any(realValues2) & any(realValues3):
                                        avgProf = np.nanmean((newdata2['z'], newdata3['z']), axis=0)
                                        diff3 = np.abs(newdata3['z'] - avgProf)
                                        diff2 = np.abs(newdata2['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff3, diff2)), axis=0)
                                    elif any(realValues3) & any(realValues4):
                                        avgProf = np.nanmean((newdata3['z'], newdata4['z']), axis=0)
                                        diff3 = np.abs(newdata3['z'] - avgProf)
                                        diff4 = np.abs(newdata4['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff3, diff4)), axis=0)
                                    elif any(realValues3) & any(realValues5):
                                        avgProf = np.nanmean((newdata3['z'], newdata5['z']), axis=0)
                                        diff3 = np.abs(newdata3['z'] - avgProf)
                                        diff5 = np.abs(newdata5['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff3, diff5)), axis=0)
                                    elif any(realValues4) & any(realValues5):
                                        avgProf = np.nanmean((newdata4['z'], newdata5['z']), axis=0)
                                        diff4 = np.abs(newdata4['z'] - avgProf)
                                        diff5 = np.abs(newdata5['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff4, diff5)), axis=0)
                                    else:
                                        if any(realValues):
                                            avgProf = newdata['z']
                                            tempDiff = 0 * avgProf
                                            # avgProf = np.nanmean((newdata['z']), axis=0)
                                        elif any(realValues2):
                                            avgProf = newdata2['z']
                                            tempDiff = 0 * avgProf
                                            # avgProf = np.nanmean((newdata['z'], newdata2['z']), axis=0)
                                        elif any(realValues3):
                                            avgProf = newdata3['z']
                                            tempDiff = 0 * avgProf
                                            # avgProf = np.nanmean((newdata['z'], newdata3['z']), axis=0)
                                        elif any(realValues4):
                                            avgProf = newdata4['z']
                                            tempDiff = 0 * avgProf
                                            # avgProf = np.nanmean((newdata['z'], newdata4['z']), axis=0)
                                        elif any(realValues5):
                                            avgProf = newdata5['z']
                                            tempDiff = 0 * avgProf

                        highResInterp = interpBathy(xinterp, avgProf, higherResX)
                        newZ = highResInterp['z']-referenceHeight
                        absz1 = np.abs(newZ)
                        minIndex1 = np.where((np.nanmin(absz1)==absz1))

                        refProf = newZ[minIndex1[0][0]:(int(minIndex1[0][0]+(len(higherResX)))):5]
                        refProf = refProf[0:len(xinterpRef)]

                        alllines = refProf
                        if file_params[1] == 'geomorphology':
                            surveyType = 'GPS'
                        else:
                            surveyType = file_params[6]
                        time = surveydate
                        maxDiff = tempDiff
                        count = count+1
                    else:
                        # # if count == 0:
                        # Do we have data at all 5 transects?
                        if any(realValues) & any(realValues2) & any(realValues3) & any(realValues4) & any(realValues5):
                            avgProf = np.nanmean(
                                (newdata['z'], newdata2['z'], newdata3['z'], newdata4['z'], newdata5['z']), axis=0)
                            diff1 = np.abs(newdata['z']-avgProf)
                            diff2 = np.abs(newdata2['z']-avgProf)
                            diff3 = np.abs(newdata3['z']-avgProf)
                            diff4 = np.abs(newdata4['z']-avgProf)
                            diff5 = np.abs(newdata5['z']-avgProf)
                            tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((np.vstack((diff1,diff2)),diff3)),diff4)),diff5)),axis=0)

                        else:
                            # What about data at only 4 of the 5?
                            if any(realValues) & any(realValues2) & any(realValues3) & any(realValues4):
                                avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata3['z'], newdata4['z']),
                                                     axis=0)
                                diff1 = np.abs(newdata['z'] - avgProf)
                                diff2 = np.abs(newdata2['z'] - avgProf)
                                diff3 = np.abs(newdata3['z'] - avgProf)
                                diff4 = np.abs(newdata4['z'] - avgProf)
                                tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((diff1, diff2)), diff3)), diff4)), axis=0)

                            elif any(realValues) & any(realValues3) & any(realValues4) & any(realValues5):
                                avgProf = np.nanmean((newdata['z'], newdata3['z'], newdata4['z'], newdata5['z']),
                                                     axis=0)
                                diff1 = np.abs(newdata['z'] - avgProf)
                                diff5 = np.abs(newdata5['z'] - avgProf)
                                diff3 = np.abs(newdata3['z'] - avgProf)
                                diff4 = np.abs(newdata4['z'] - avgProf)
                                tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((diff1, diff5)), diff3)), diff4)), axis=0)
                            elif any(realValues) & any(realValues2) & any(realValues4) & any(realValues5):
                                avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata4['z'], newdata5['z']),
                                                     axis=0)
                                diff1 = np.abs(newdata['z'] - avgProf)
                                diff2 = np.abs(newdata2['z'] - avgProf)
                                diff5 = np.abs(newdata5['z'] - avgProf)
                                diff4 = np.abs(newdata4['z'] - avgProf)
                                tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((diff1, diff2)), diff5)), diff4)), axis=0)
                            elif any(realValues) & any(realValues2) & any(realValues3) & any(realValues5):
                                avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata3['z'], newdata5['z']),
                                                     axis=0)
                                diff1 = np.abs(newdata['z'] - avgProf)
                                diff2 = np.abs(newdata2['z'] - avgProf)
                                diff3 = np.abs(newdata3['z'] - avgProf)
                                diff5 = np.abs(newdata5['z'] - avgProf)
                                tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((diff1, diff2)), diff3)), diff5)), axis=0)
                            elif any(realValues2) & any(realValues3) & any(realValues4) & any(realValues5):
                                avgProf = np.nanmean((newdata2['z'], newdata3['z'], newdata4['z'], newdata5['z']),
                                                     axis=0)
                                diff5 = np.abs(newdata5['z'] - avgProf)
                                diff2 = np.abs(newdata2['z'] - avgProf)
                                diff3 = np.abs(newdata3['z'] - avgProf)
                                diff4 = np.abs(newdata4['z'] - avgProf)
                                tempDiff = np.nanmax(np.vstack((np.vstack((np.vstack((diff5, diff2)), diff3)), diff4)), axis=0)
                            else:
                                # What about if data is at only 3 of the 5?
                                if any(realValues) & any(realValues2) & any(realValues3):
                                    avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata3['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff2)), diff3)), axis=0)
                                elif any(realValues) & any(realValues2) & any(realValues4):
                                    avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata4['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff2)), diff4)), axis=0)
                                elif any(realValues) & any(realValues2) & any(realValues5):
                                    avgProf = np.nanmean((newdata['z'], newdata2['z'], newdata5['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff5 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff2)), diff5)), axis=0)
                                elif any(realValues) & any(realValues3) & any(realValues4):
                                    avgProf = np.nanmean((newdata['z'], newdata3['z'], newdata4['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff4)), diff3)), axis=0)
                                elif any(realValues) & any(realValues3) & any(realValues5):
                                    avgProf = np.nanmean((newdata['z'], newdata3['z'], newdata5['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff5 = np.abs(newdata5['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff5)), diff3)), axis=0)
                                elif any(realValues) & any(realValues4) & any(realValues5):
                                    avgProf = np.nanmean((newdata['z'], newdata4['z'], newdata5['z']), axis=0)
                                    diff1 = np.abs(newdata['z'] - avgProf)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    diff5 = np.abs(newdata5['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff1, diff4)), diff5)), axis=0)
                                elif any(realValues2) & any(realValues3) & any(realValues4):
                                    avgProf = np.nanmean((newdata2['z'], newdata3['z'], newdata4['z']), axis=0)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff4, diff2)), diff3)), axis=0)
                                elif any(realValues2) & any(realValues3) & any(realValues5):
                                    avgProf = np.nanmean((newdata2['z'], newdata3['z'], newdata5['z']), axis=0)
                                    diff5 = np.abs(newdata5['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff5, diff2)), diff3)), axis=0)
                                elif any(realValues2) & any(realValues4) & any(realValues5):
                                    avgProf = np.nanmean((newdata2['z'], newdata4['z'], newdata5['z']), axis=0)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    diff2 = np.abs(newdata2['z'] - avgProf)
                                    diff5 = np.abs(newdata5['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff4, diff2)), diff5)), axis=0)
                                elif any(realValues3) & any(realValues4) & any(realValues5):
                                    avgProf = np.nanmean((newdata3['z'], newdata4['z'], newdata5['z']), axis=0)
                                    diff4 = np.abs(newdata4['z'] - avgProf)
                                    diff5 = np.abs(newdata5['z'] - avgProf)
                                    diff3 = np.abs(newdata3['z'] - avgProf)
                                    tempDiff = np.nanmax(np.vstack((np.vstack((diff4, diff5)), diff3)), axis=0)


                                else:
                                    if any(realValues) & any(realValues2):
                                        avgProf = np.nanmean((newdata['z'], newdata2['z']), axis=0)
                                        diff1 = np.abs(newdata['z'] - avgProf)
                                        diff2 = np.abs(newdata2['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff1, diff2)), axis=0)
                                    elif any(realValues) & any(realValues3):
                                        avgProf = np.nanmean((newdata['z'], newdata3['z']), axis=0)
                                        diff1 = np.abs(newdata['z'] - avgProf)
                                        diff3 = np.abs(newdata3['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff1, diff3)), axis=0)
                                    elif any(realValues) & any(realValues4):
                                        avgProf = np.nanmean((newdata['z'], newdata4['z']), axis=0)
                                        diff1 = np.abs(newdata['z'] - avgProf)
                                        diff4 = np.abs(newdata4['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff1, diff4)), axis=0)
                                    elif any(realValues) & any(realValues5):
                                        avgProf = np.nanmean((newdata['z'], newdata5['z']), axis=0)
                                        diff1 = np.abs(newdata['z'] - avgProf)
                                        diff5 = np.abs(newdata5['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff1, diff5)), axis=0)
                                    elif any(realValues2) & any(realValues5):
                                        avgProf = np.nanmean((newdata2['z'], newdata5['z']), axis=0)
                                        diff5 = np.abs(newdata5['z'] - avgProf)
                                        diff2 = np.abs(newdata2['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff5, diff2)), axis=0)
                                    elif any(realValues2) & any(realValues4):
                                        avgProf = np.nanmean((newdata2['z'], newdata4['z']), axis=0)
                                        diff4 = np.abs(newdata4['z'] - avgProf)
                                        diff2 = np.abs(newdata2['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff4, diff2)), axis=0)
                                    elif any(realValues2) & any(realValues3):
                                        avgProf = np.nanmean((newdata2['z'], newdata3['z']), axis=0)
                                        diff3 = np.abs(newdata3['z'] - avgProf)
                                        diff2 = np.abs(newdata2['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff3, diff2)), axis=0)
                                    elif any(realValues3) & any(realValues4):
                                        avgProf = np.nanmean((newdata3['z'], newdata4['z']), axis=0)
                                        diff3 = np.abs(newdata3['z'] - avgProf)
                                        diff4 = np.abs(newdata4['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff3, diff4)), axis=0)
                                    elif any(realValues3) & any(realValues5):
                                        avgProf = np.nanmean((newdata3['z'], newdata5['z']), axis=0)
                                        diff3 = np.abs(newdata3['z'] - avgProf)
                                        diff5 = np.abs(newdata5['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff3, diff5)), axis=0)
                                    elif any(realValues4) & any(realValues5):
                                        avgProf = np.nanmean((newdata4['z'], newdata5['z']), axis=0)
                                        diff4 = np.abs(newdata4['z'] - avgProf)
                                        diff5 = np.abs(newdata5['z'] - avgProf)
                                        tempDiff = np.nanmax(np.vstack((diff4, diff5)), axis=0)
                                    else:
                                        if any(realValues):
                                            avgProf = newdata['z']
                                            tempDiff = 0*avgProf
                                            # avgProf = np.nanmean((newdata['z']), axis=0)
                                        elif any(realValues2):
                                            avgProf = newdata2['z']
                                            tempDiff = 0*avgProf
                                            # avgProf = np.nanmean((newdata['z'], newdata2['z']), axis=0)
                                        elif any(realValues3):
                                            avgProf = newdata3['z']
                                            tempDiff = 0*avgProf
                                            # avgProf = np.nanmean((newdata['z'], newdata3['z']), axis=0)
                                        elif any(realValues4):
                                            avgProf = newdata4['z']
                                            tempDiff = 0*avgProf
                                            # avgProf = np.nanmean((newdata['z'], newdata4['z']), axis=0)
                                        elif any(realValues5):
                                            avgProf = newdata5['z']
                                            tempDiff = 0*avgProf

                        highResInterp = interpBathy(xinterp, avgProf, higherResX)
                        newZ = highResInterp['z']-referenceHeight
                        absz1 = np.abs(newZ)
                        minIndex1 = np.where((np.nanmin(absz1)==absz1))

                        refProf = newZ[minIndex1[0][0]:(int(minIndex1[0][0]+(len(higherResX)))):5]
                        refProf = refProf[0:len(xinterpRef)]

                        alllines = np.vstack((alllines, refProf))
                            #alllines = np.vstack((alllines, avgProf))
                        if file_params[1] == 'geomorphology':
                            surveyType = np.append(surveyType,'GPS')
                        else:
                            surveyType = np.append(surveyType,file_params[6])
                        time = np.append(time,surveydate)
                        maxDiff = np.vstack((maxDiff,tempDiff))

                        count = count+1

                else:
                    print('Data is less than 10 points long at line {}'.format(i))
                    newdata = []
            else:
                print('No data onshore of 125 meters at line {}'.format(i))
                print('Most onshore point at {}'.format(np.nanmin(tempXSubNew)))
                badCounterLong = badCounterLong + 1

        else:
            print('No data deeper than 600 meters at line {}'.format(i))
            print('Most offshore point at {}'.format(np.nanmax(tempXSubNew)))
            badCounterShort = badCounterShort + 1

    else:
        print('Survey with no data at this line for survey {}, count of {}'.format(i,count))
        badCounter = badCounter+1




volumes = np.empty(len(alllines,))
from numpy import trapz
for i in range(len(alllines)):
    volumes[i] = trapz(alllines[i,:]+10,dx=5)


# # Plotting all of the lines
t1 = 0
t2 = -1
subsetOfAllines = alllines[t1:t2]
fig, ax = plt.subplots(2,1)
ax[0].set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, len(subsetOfAllines)))))
for i in range(len(subsetOfAllines)):
    ax[0].plot(xinterpRef, subsetOfAllines[i,:], label=time[i])

ax[0].set_xlim([0, 550])

ax[0].set_ylim([-8, 4])
ax[0].set_title('Cross-shore profile variability')

plt.set_cmap('RdBu_r')

tg, xg = np.meshgrid(time, xinterpRef)
plt0 = ax[1].pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
fig.colorbar(plt0, ax=ax[1], orientation='horizontal')
ax[1].set_ylim([time[t1], time[t2]])
ax[1].set_title('Surveys (dev.)')

plt.style.use('default')
plt.set_cmap('RdBu_r')

figs = plt.figure()
ax = plt.subplot2grid((1,1),(0,0),rowspan=1,colspan=1)
pl = ax.pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
ax.set_ylim([time[t1], time[t2]])
ax.set_title('Surveys (dev.)')
ax.set_xlabel('Cross-shore (m)')
cb = figs.colorbar(pl, ax=ax, orientation='horizontal')
cb.set_label('Deviation from mean (m)')
plt.tight_layout()

plt.figure()
groups = np.unique(surveyType)
for name in groups:
    matched_indexes = []
    i = 0
    length = len(surveyType)

    while i < length:
        if name == surveyType[i]:
            matched_indexes.append(i)
        i += 1

    plt.scatter(time[matched_indexes],volumes[matched_indexes],marker='o',label=name)
plt.legend()
plt.ylabel('Volumes (m^3/m alongshore)')

tg, xg = np.meshgrid(time, xinterp)



morphoPickle = 'sandbarsSouthernTransect_referencedMHHW.pickle'
output = {}
output['time'] = time
output['alllines'] = alllines
output['xinterp'] = xinterpRef

import pickle
with open(morphoPickle,'wb') as f:
    pickle.dump(output, f)


