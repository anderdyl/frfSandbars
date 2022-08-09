import pickle
import numpy as np
from matplotlib import pyplot as plt

def P2R(radii, angles):
   return radii * np.exp(1j * angles)

def R2P(x):
   return np.abs(x), np.angle(x)

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
days = [t.days for t in dt]
days2 = np.array(days)

n = 50

timeind = np.arange(0,len(time))
alllinesSubset = alllines[timeind, :]

eofPredog = np.nan * np.ones((np.shape(alllinesSubset)))
eofPred2og = np.nan * np.ones((np.shape(alllinesSubset)))

for timestep in range(len(timeind)):
    mode = 0
    eofPredog[timestep, :] = RS[timestep, mode]* SS[:, mode] * np.cos(phitS[timestep, mode] - thetaS[:, mode])
    mode = 1
    eofPred2og[timestep, :] = RS[timestep, mode]* SS[:, mode] * np.cos(phitS[timestep, mode] - thetaS[:, mode])

thetasMode1 = np.arange(0,1200,5)
thetasMode2 = np.arange(0,1200,5)
eofHypo = np.nan * np.ones((len(thetasMode1),len(xinterpS)))
eofHypo2 = np.nan * np.ones((len(thetasMode1),len(xinterpS)))

for tt in range(len(thetasMode1)):
   eofHypo[tt,:] = SS[:, 0] * np.nanmean(RS[:, 0]) * np.cos(thetasMode1[tt] * np.pi / 180 - thetaS[:, 0]+180)
   eofHypo2[tt,:] = SS[:, 1] * np.nanmean(RS[:, 1]) * np.cos(thetasMode2[tt] * np.pi / 180 - thetaS[:, 1])

thetag, xintg = np.meshgrid(thetasMode1, xinterpS)

t1 = 0
t2 = -1
fig = plt.figure(figsize=(14,10))
tg, xg = np.meshgrid(time, xinterp)
ax2 = plt.subplot2grid((7,4),(0,0),rowspan=7,colspan=1)
plt1 = ax2.pcolor(xg,tg,eofPredog.T, vmin=-.75, vmax=0.75,cmap='bwr')
ax2.set_ylim([time[t1], time[t2]])
ax2.set_title('Offshore CEOF (40%)',fontsize=12)
ax2.set_xlabel('Cross-shore (m)')
ax2.text(-0.1, 1.01, 'a.', transform=ax2.transAxes, size=14, weight='bold')
ax3 = plt.subplot2grid((7,4),(0,1),rowspan=7,colspan=1)
plt2 = ax3.pcolor(xg,tg,eofPred2og.T, vmin=-.85, vmax=.85,cmap='bwr')
ax3.set_ylim([time[t1], time[t2]])
ax3.set_title('Onshore CEOF (25%)',fontsize=12)
ax3.set_xlabel('Cross-shore (m)')
ax3.text(-0.1, 1.01, 'b.', transform=ax3.transAxes, size=14, weight='bold')

ax1 = plt.subplot2grid((7,4),(0,2),rowspan=3,colspan=1)
ax1b = ax1.twinx()
plt0 = ax1.pcolor(xintg,thetag,eofHypo.T, vmin=-0.75, vmax=0.75, cmap='bwr')
ax1.plot(xintg[:,0],360*np.ones((len(xintg[:,0],))),'--',color='orange')
ax1.plot(xintg[:,0],180*np.ones((len(xintg[:,0],))),'--',color='green')
plt0b = ax1b.plot(range(1200),np.ones(1200))
ax1b.cla()
ax1b.set_ylim(0,810)
ax1.set_ylim(0,810)
ax1b.set_yticks([0,98,196,294,392,490,589,687,785])
ax1b.set_yticklabels(['0','','1 yr','','2 yrs','','3 yrs','','4 yrs'])
ax1b.set_ylabel('Average Time Elapsed',fontsize=12)
ax1.set_title(r'Offshore Migration $\Longrightarrow$',fontsize=12)
ax1.set_ylabel('Phase (deg)',fontsize=12)
ax1.set_xlim(0,510)
ax1.text(-0.05, 1.045, 'c.', transform=ax1.transAxes, size=14, weight='bold')
ax22 = plt.subplot2grid((7,4),(3,2),rowspan=3,colspan=1)
ax22b = ax22.twinx()
plt12 = ax22.pcolor(xintg,thetag,eofHypo2.T, vmin=-0.75, vmax=0.75, cmap='bwr')
ax22.plot(xintg[:,0],20*np.ones((len(xintg[:,0],))),'--',color='orange')
ax22.plot(xintg[:,0],240*np.ones((len(xintg[:,0],))),'--',color='green')
plt22b = ax22b.plot(range(1200),np.ones(1200))
ax22b.cla()
ax22b.set_yticks([0,132,263,394,525,656,788,920,1052])
ax22b.set_yticklabels(['0','','1 yr','','2 yrs','','3 yrs','','4 yrs'])
ax22b.set_ylabel('Average Time Elapsed',fontsize=12)
ax22.set_ylabel('Phase (deg)',fontsize=12)
ax22.set_title(r'Onshore Migration $\Longleftarrow$',fontsize=12)
ax22.set_xlabel('Cross-shore',fontsize=12)
ax22.set_ylim([0,1100])
ax22b.set_ylim(0,1100)
ax22.set_xlim([0,510])
cbar_ax = fig.add_axes([0.58, 0.1, 0.12, 0.03])
cb1c = plt.colorbar(plt12, cax=cbar_ax, orientation='horizontal')
cb1c.set_label('Difference from Mean (m)')
ax22.text(-0.05, 1.03, 'd.', transform=ax22.transAxes, size=14, weight='bold')
ax31 = plt.subplot2grid((7,4),(0,3),rowspan=3,colspan=1)
ax31.plot(xintg[:,0],eofHypo2.T[:,0]+eofHypo.T[:,4]+np.mean(alllinesS,axis=0),color='orange')
ax31.set_ylabel('Depth (m)')
ax31.set_title('Example Profile Reconstruction',fontsize=12)
ax31.text(-0.18, 1.03, 'e.', transform=ax31.transAxes, size=14, weight='bold')
ax32 = plt.subplot2grid((7,4),(3,3),rowspan=3,colspan=1)
ax32.plot(xintg[:,0],eofHypo2.T[:,44]+eofHypo.T[:,36]+np.mean(alllinesS,axis=0),color='green')
ax32.set_title('Example Profile Reconstruction',fontsize=12)
ax32.set_ylabel('Depth (m)')
ax32.set_xlabel('Cross-shore (m)')
ax32.text(-0.18, 1.03, 'f.', transform=ax32.transAxes, size=14, weight='bold')
plt.tight_layout()
plt.show()

