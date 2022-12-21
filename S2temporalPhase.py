import numpy as np
import matplotlib.pyplot as plt
import pickle

dbfile = open('sandbarsSouthernTransect_referencedMHHW.pickle', 'rb')
data = pickle.load(dbfile)
dbfile.close()
alllines = data['alllines']
xinterp = data['xinterp']
time = data['time']

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


def moving_average(a, n=5):
   ret = np.cumsum(a, dtype=float)
   ret[n:] = ret[n:] - ret[:-n]
   return ret[n - 1:] / n

def weightedMovingAverage(a,b,n=5):
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

mode1Phase = weightedMovingAverage(np.unwrap(phit2S[:,0]*(np.pi/180)),RS[:,0],n=5)
mode2Phase = weightedMovingAverage(np.unwrap(phit2S[:,1]*(np.pi/180)),RS[:,1],n=5)
mode3Phase = weightedMovingAverage(np.unwrap(phit2S[:,2]*(np.pi/180)),RS[:,2],n=5)

mode1Mag = moving_average(RS[:,0],n=5)
mode2Mag = moving_average(RS[:,1],n=5)
mode3Mag = moving_average(RS[:,2],n=5)

mode1Wrap = (mode1Phase + np.pi) % (2 * np.pi) - np.pi
mode2Wrap = (mode2Phase + np.pi) % (2 * np.pi) - np.pi
mode3Wrap = (mode3Phase + np.pi) % (2 * np.pi) - np.pi

index = np.where((mode1Wrap<0))
mode1Wrap[index[0]] = mode1Wrap[index[0]]+np.pi*2
mode = 0

fig = plt.figure(figsize=(10,5))
ax1 = plt.subplot2grid((3,2),(0,0),rowspan=1,colspan=2)
xg,tg = np.meshgrid(xinterpS,timeS)
plt0 = ax1.pcolor(tg,xg,(alllinesS-np.mean(alllinesS, axis=0)), vmin=-1.8, vmax=1.8,cmap='bwr')
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax1)
cax4 = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(plt0, cax=cax4, orientation='vertical');
cb.set_label('$\~z(x,t) (m)$')

ax = plt.subplot2grid((3,2),(1,0),rowspan=1,colspan=2,sharex=ax1)
# temp = RS[:,mode]*10
# temp1 = phit2S[:,mode]
temp = mode1Mag*10
temp1 = mode1Wrap
p3 = ax.scatter(timeS[4:-4],temp1,temp,mode1Mag,cmap='Reds')
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax)
cax3 = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(p3, cax=cax3, orientation='vertical');
cb.set_label('$R_{1}(t)$')

mode = 1
ax2 = plt.subplot2grid((3,2),(2,0),rowspan=1,colspan=2, sharex=ax1)
temp = mode2Mag*10
temp1 = mode2Wrap
p2 = ax2.scatter(timeS[4:-4],temp1,temp,mode2Mag,cmap='Blues')

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(p2, cax=cax, orientation='vertical');
cb.set_label('$R_{2}(t)$')

ax.set_xlim([timeS[0], timeS[-1]])
ax1.set_xlim([timeS[0], timeS[-1]])
ax2.set_xlim([timeS[0], timeS[-1]])
ax1.set_ylabel('Cross-shore (m)')
ax.set_yticks([1,3,5])
ax.set_yticklabels(['-2','0','2'])
ax.set_ylabel('$\Omega(t)$ (radians)')
ax2.set_ylabel('$\Omega(t)$ (radians)')
ax2.text(-0.05, 1.02, 'c.', transform=ax2.transAxes, size=14, weight='bold')
ax.text(-0.05, 1.02, 'b.', transform=ax.transAxes, size=14, weight='bold')
ax1.text(-0.05, 1.02, 'a.', transform=ax1.transAxes, size=14, weight='bold')
ax2.set_xlabel('Time')
plt.show()
