import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la
import numpy.linalg as npla
from scipy.signal import hilbert
from scipy.interpolate import interp1d
import pickle


dbfile = open('sandbarsSouthernTransect_referencedMHHW.pickle', 'rb')

data = pickle.load(dbfile)
dbfile.close()
alllinesOG = data['alllines']
xinterpOG = data['xinterp']
time = data['time']

def interpBathy(xSub, zSub, x):
    f = interp1d(xSub, zSub, kind='linear', bounds_error=False)
    newz = f(x)
    newBathy = dict()
    newBathy['x'] = x
    newBathy['z'] = newz
    return newBathy

xinterp = np.arange(0,510,5)
alllines = np.nan*np.ones((len(time),len(xinterp)))

for hh in range(len(time)):
    lowerRes = interpBathy(xinterpOG, alllinesOG[hh,:], xinterp)
    alllines[hh,:] = lowerRes['z']


nanmean = np.nanmean(alllines,axis=0)
finder = np.where(np.isnan(alllines))
alllines[finder] = nanmean[finder[1]]
demean = alllines - np.mean(alllines,axis=0)

data = (hilbert(demean.T))

data = data.T
c = np.matmul(np.conj(data).T,data)/np.shape(data)[0]




lamda, loadings = la.eigh(c)

lamda2, loadings2 = npla.eig(c)

ind = np.argsort(lamda[::-1])

lamda[::-1].sort()

loadings = loadings[:,ind]

pcs = np.dot(data, loadings)# / np.sqrt(lamda)
loadings = loadings# * np.sqrt(lamda)
pcsreal = np.real(pcs[:,0:200])
pcsimag = np.imag(pcs[:,0:200])
eofreal = np.real(loadings[:,0:200])
eofimag = np.imag(loadings[:,0:200])
S = np.power(loadings*np.conj(loadings),0.5) * np.sqrt(lamda)

theta = np.arctan2(eofimag,eofreal)
theta2 = theta*180/np.pi

Rt = np.power(pcs*np.conj(pcs),0.5) / np.sqrt(lamda)

phit = np.arctan2(pcsimag,pcsreal)
phit2 = phit*180/np.pi


mode = 0

fig, ax = plt.subplots(2,2)

ax[0,0].plot(xinterp, S[:,mode],'o')
ax[0,0].set_ylabel('Spatial Magnitude (m)')
ax[0,0].set_xlabel('Cross-shore (m)')
ax[1,0].plot(xinterp, theta2[:,mode],'o')
ax[1,0].set_ylabel('Spatial Phase (deg)')
ax[1,0].set_xlabel('Cross-shore (m)')
ax[0,1].plot(time,Rt[:,mode],'o')
ax[0,1].set_ylabel('Temporal Magnitude (m)')
ax[0,1].set_xlabel('Time')
ax[1,1].plot(time,phit2[:,mode],'o')
ax[1,1].set_ylabel('Temporal Phase (deg)')
ax[1,1].set_xlabel('Time')



PC1 = Rt[:, mode]*np.sin(phit[:, mode]) + Rt[:, mode]*np.cos(phit[:, mode])


totalV = np.sum(lamda)
percentV = lamda / totalV

ztemp = 0*np.ones(len(xinterp),)
timestep = 200
for mode in range(2):
    # ztemp = ztemp + Rt[timestep,mode]*np.sin(phit[timestep,mode]) * S[:,mode]*np.sin(theta[:,mode]) + Rt[timestep,mode]*np.cos(phit[timestep,mode]) * S[:,mode]*np.cos(theta[:,mode])
    ztemp = ztemp + Rt[timestep,mode]*S[:,mode]*np.cos(phit[timestep,mode] - theta[:,mode])


def P2R(radii, angles):
    return radii * np.exp(1j*angles)

def R2P(x):
    return np.abs(x), np.angle(x)



timeind = np.arange(0,len(time))
#timeind = np.arange(3, 50)
RtSubset = Rt[timeind, :]
phitSubset = phit[timeind, :]
phit2Subset = phit2[timeind, :]
timeSubset = time[timeind]
alllinesSubset = alllines[timeind, :]

eofPred = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred2 = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred3 = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred4 = np.nan * np.ones((np.shape(alllinesSubset)))

for timestep in range(len(timeind)):
    mode = 0
    eofPred[timestep, :] = RtSubset[timestep, mode]* S[:, mode] * np.cos(phitSubset[timestep, mode] - theta[:, mode])
    mode = 1
    eofPred2[timestep, :] = RtSubset[timestep, mode]* S[:, mode] * np.cos(phitSubset[timestep, mode] - theta[:, mode])

    mode = 2
    eofPred3[timestep, :] = RtSubset[timestep, mode]* S[:, mode] * np.cos(phitSubset[timestep, mode] - theta[:, mode])

    mode = 3
    eofPred4[timestep, :] = RtSubset[timestep, mode]* S[:, mode] * np.cos(phitSubset[timestep, mode] - theta[:, mode])



t1 = 0
t2 = -1
fig, ax = plt.subplots(1,5)
plt.set_cmap('RdBu_r')

tg, xg = np.meshgrid(time, xinterp)
plt0 = ax[0].pcolor(xg,tg,(alllines-np.mean(alllines, axis=0)).T, vmin=-1.8, vmax=1.8)
fig.colorbar(plt0, ax=ax[0], orientation='horizontal')
ax[0].set_ylim([time[t1], time[t2]])
ax[0].set_title('Surveys (dev.)')

plt1 = ax[1].pcolor(xg,tg,eofPred.T, vmin=-.75, vmax=0.75)
ax[1].set_ylim([time[t1], time[t2]])
fig.colorbar(plt1, ax=ax[1], orientation='horizontal')
ax[1].set_title('CEOF1 {:.2f}'.format(percentV[0]))
ax[1].get_yaxis().set_ticks([])

plt2 = ax[2].pcolor(xg,tg,eofPred2.T, vmin=-.85, vmax=.85)
ax[2].set_ylim([time[t1], time[t2]])
fig.colorbar(plt2, ax=ax[2], orientation='horizontal')
ax[2].set_title('CEOF2 {:.2f}'.format(percentV[1]))
ax[2].get_yaxis().set_ticks([])

plt3 = ax[3].pcolor(xg,tg,eofPred3.T, vmin=-.65, vmax=.65)
ax[3].set_ylim([time[t1], time[t2]])
ax[3].set_title('CEOF3 {:.2f}'.format(percentV[2]))
ax[3].get_yaxis().set_ticks([])

fig.colorbar(plt3, ax=ax[3], orientation='horizontal')

plt4 = ax[4].pcolor(xg,tg,eofPred4.T, vmin=-.65, vmax=.65)
ax[4].set_ylim([time[t1], time[t2]])
ax[4].set_title('CEOF4 {:.2f}'.format(percentV[3]))
ax[4].get_yaxis().set_ticks([])
fig.colorbar(plt4, ax=ax[4], orientation='horizontal')

plt.tight_layout(pad=0.5)

plt.show()


totalV = np.sum(lamda)
percentV = lamda / totalV
import datetime
def datetime2matlabdn(dt):
   ord = dt.toordinal()
   mdn = dt + datetime.timedelta(days = 366)
   frac = (dt-datetime.datetime(dt.year,dt.month,dt.day,0,0,0)).seconds / (24.0 * 60.0 * 60.0)
   return mdn.toordinal() + frac

matlabTime = np.zeros((len(time),))
for qq in range(len(time)):
    matlabTime[qq] = datetime2matlabdn(time[qq])


import scipy.io
output = dict()
output['alllines'] = alllines
output['xinterp'] = xinterp
output['time'] = matlabTime
output['S'] = S
output['Rt'] = Rt
output['thetaRadians'] = theta
output['thetaDegrees'] = theta2
output['phiRadian'] = phit
output['phiDegrees'] = phit2
output['totalV'] = totalV
output['lamda'] = lamda
output['percentV'] = percentV
scipy.io.savemat('ceofsSouthernTransect.mat',output)


morphoPickle = 'ceofsSouthernTransect.pickle'
output = {}
output['time'] = time
output['alllines'] = alllines
output['xinterp'] = xinterp
output['S'] = S
output['Rt'] = Rt
output['thetaRadians'] = theta
output['thetaDegrees'] = theta2
output['phiRadian'] = phit
output['phiDegrees'] = phit2
output['totalV'] = totalV
output['lamda'] = lamda
output['percentV'] = percentV


import pickle
with open(morphoPickle,'wb') as f:
    pickle.dump(output, f)






