import matplotlib.pyplot as plt
import pickle
import numpy as np
import sandbarUtils as sf

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

thetasMode1 = np.arange(0,1200,5)
thetasMode2 = np.arange(0,1200,5)

eofHypo = np.nan * np.ones((len(thetasMode1),len(xinterpS)))
eofHypo2 = np.nan * np.ones((len(thetasMode1),len(xinterpS)))

for tt in range(len(thetasMode1)):
   eofHypo[tt,:] = SS[:, 0] * np.nanmean(RS[:, 0]) * np.cos(thetasMode1[tt] * np.pi / 180 - thetaS[:, 0]+180)
   eofHypo2[tt,:] = SS[:, 1] * np.nanmean(RS[:, 1]) * np.cos(thetasMode2[tt] * np.pi / 180 - thetaS[:, 1])

thetag, xintg = np.meshgrid(thetasMode1, xinterpS)

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


for tt in range(len(thetasMode1)):
   eofPred = SS[:, 0] * np.nanmean(RS[:, 0]) * np.cos(thetasMode1[tt] * np.pi / 180 - thetaS[:, 0])
   for yy in range(len(thetasMode2)):
      eofPred2 = SS[:, 1] * np.nanmean(RS[:, 1]) * np.cos(thetasMode2[yy] * np.pi / 180 - thetaS[:, 1])
      combined = np.real(np.mean(alllinesS,axis=0)+eofPred+eofPred2)
      xFRFbar, xFRFtrough = sf.findSandBarAndTrough1D(xinterpS, combined, plotFname=None, smoothLengthScale=5, profileTrend=np.mean(alllines,axis=0))
      if xFRFbar is not None:

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

         else:

            innerOuterBarX[yy,tt] = xFRFbar
            zOuterInd = np.where((np.round(2*np.max(xFRFbar))/2 == xinterpS))
            innerOuterBarZ[yy,tt] = combined[zOuterInd[0]]
            innerOuterTroughZ[yy, tt] = xFRFtrough
            innerOuterTroughZInd = np.where((np.round(2 * np.max(xFRFtrough)) / 2 == xinterpS))
            innerOuterTroughZ[yy, tt] = combined[innerOuterTroughZInd[0]]
            innerOuterBarDZ[yy,tt] = innerOuterBarZ[yy,tt]-innerOuterTroughZ[yy,tt]

meshTheta1, meshTheta2 = np.meshgrid(thetasMode1,thetasMode2)

import copy
doublebar = copy.deepcopy(outerBarX)
doubleIndex = np.where((doublebar > 0))
doublebar[doubleIndex] = 0.5

singlebar = copy.deepcopy(innerOuterBarX)
singleIndex = np.where((singlebar > 0))
doublebar[singleIndex] = 0.25

fig = plt.figure(figsize=(10,11))
ax1 = plt.subplot2grid((5,3),(3,0),rowspan=1,colspan=1)
p1 = ax1.pcolor(meshTheta1,meshTheta2,innerBarX,vmin=30,vmax=300,cmap='plasma')
ax1.text(0.04, .85, 'e.', transform=ax1.transAxes, size=14, weight='bold')

ax2 = plt.subplot2grid((5,3),(3,1),rowspan=1,colspan=1)
p2 = ax2.pcolor(meshTheta1,meshTheta2,outerBarX,vmin=30,vmax=300,cmap='plasma')
ax2.text(0.04, .85, 'f.', transform=ax2.transAxes, size=14, weight='bold')

ax1.set_ylabel(r'Onshore (deg) $\Longrightarrow$')
ax1.set_title('Inner Bar',fontsize=12)
ax2.set_title('Outer Bar',fontsize=12)
ax1c = plt.subplot2grid((5,3),(3,2),rowspan=1,colspan=1)
p1c = ax1c.pcolor(meshTheta1,meshTheta2,innerOuterBarX,vmin=30,vmax=300,cmap='plasma')
ax1c.text(0.04, .85, 'g.', transform=ax1c.transAxes, size=14, weight='bold',color='white')

ax1c.set_title('Single Bar',fontsize=12)
ax3 = plt.subplot2grid((5,3),(4,0),rowspan=1,colspan=1)
p3 = ax3.pcolor(meshTheta1,meshTheta2,innerBarDZ,vmin=0,vmax=1)
ax3.text(0.04, .85, 'h.', transform=ax3.transAxes, size=14, weight='bold')

ax4 = plt.subplot2grid((5,3),(4,1),rowspan=1,colspan=1)
p4 = ax4.pcolor(meshTheta1,meshTheta2,outerBarDZ,vmin=0,vmax=1)
ax4.text(0.04, .85, 'i.', transform=ax4.transAxes, size=14, weight='bold')

ax4c = plt.subplot2grid((5,3),(4,2),rowspan=1,colspan=1)
p4c = ax4c.pcolor(meshTheta1,meshTheta2,innerOuterBarDZ,vmin=0,vmax=1)
ax4c.text(0.04, .85, 'j.', transform=ax4c.transAxes, size=14, weight='bold',color='white')

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

xMesh,yMesh = np.meshgrid(np.arange(0,1,0.1),np.arange(0,1,0.1))
label1 = np.zeros((np.shape(xMesh)))
shade1_ax = fig.add_axes([0.46,0.45,0.15,0.02])
shade1_ax.pcolor(xMesh,yMesh,label1,vmin=0,vmax=1,cmap='Greys')
shade1_ax.set_xticklabels('')
shade1_ax.set_yticklabels('')
shade1_ax.set_yticks([])
shade1_ax.set_xticks([])
shade1_ax.text(0.05,0.2,'Terrace (no bar)',weight='bold')

label2 = 0.25 * np.ones((np.shape(xMesh)))
shade2_ax = fig.add_axes([0.64,0.45,0.15,0.02])
shade2_ax.pcolor(xMesh,yMesh,label2,vmin=0,vmax=1,cmap='Greys')
shade2_ax.set_xticklabels('')
shade2_ax.set_yticklabels('')
shade2_ax.set_xticks([])
shade2_ax.set_yticks([])
shade2_ax.text(0.2,0.2,'Single Bar',weight='bold')

label3 = 0.5 * np.ones((np.shape(xMesh)))
shade3_ax = fig.add_axes([0.82,0.45,0.15,0.02])
shade3_ax.pcolor(xMesh,yMesh,label3,vmin=0,vmax=1,cmap='Greys')
shade3_ax.set_xticklabels('')
shade3_ax.set_yticklabels('')
shade3_ax.set_xticks([])
shade3_ax.set_yticks([])
shade3_ax.text(0.2,0.2,'Double Bar',weight='bold')


conceptbar_ax = fig.add_axes([0.48, 0.52, 0.48, 0.45])
con = conceptbar_ax.pcolor(meshTheta1,meshTheta2,doublebar,vmin=0,vmax=1,cmap='Greys')
conceptbar_ax.set_xlabel(r'Offshore (deg) $\Longrightarrow$')
conceptbar_ax.set_ylabel(r'Onshore (deg) $\Longrightarrow$')
conceptbar_ax.set_title('Morphologic Trends')
conceptbar_ax.arrow(x=250,y=-20,dx=0,dy=50,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.text(262,15,'Shallow Single')#,weight='bold')
conceptbar_ax.text(262,0, 'Bar Growth')#,weight='bold')
conceptbar_ax.arrow(x=80,y=100,dx=0,dy=135,width=4,facecolor='orange',edgecolor='orange')
conceptbar_ax.text(52,120, 'Inner Bar Decay',rotation=90)
conceptbar_ax.text(65,120, 'Outer Bar Maintained', rotation=90)
conceptbar_ax.text(88,120, 'Faster Inner Bar Migration', rotation=90)

conceptbar_ax.arrow(x=110,y=70,dx=100,dy=0,width=4,facecolor='blue',edgecolor='blue')
conceptbar_ax.text(120,80,'Inner Bar Growth')#,weight='bold')
conceptbar_ax.text(120,52, 'Outer Bar Decay')#,weight='bold')
conceptbar_ax.arrow(x=110,y=220,dx=100,dy=0,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.text(120,245,'Double Bar to')#,weight='bold')
conceptbar_ax.text(120,230, 'Double Terrace')#,weight='bold')
conceptbar_ax.arrow(x=80,y=0,dx=0,dy=40,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.arrow(x=80,y=0,dx=40,dy=0,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.text(90,25,'New Inner Bar')#,weight='bold')
conceptbar_ax.text(90,10, 'Formation')#,weight='bold')
conceptbar_ax.arrow(x=200,y=100,dx=0,dy=40,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.arrow(x=200,y=100,dx=40,dy=0,width=4,facecolor='green',edgecolor='green')
conceptbar_ax.text(210,125,'Outer Bar')#,weight='bold')
conceptbar_ax.text(210,110, 'Death')#,weight='bold')
conceptbar_ax.arrow(x=290,y=60,dx=0,dy=40,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.arrow(x=290,y=60,dx=40,dy=0,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.arrow(x=250,y=230,dx=0,dy=40,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.arrow(x=250,y=230,dx=40,dy=0,width=4,facecolor='black',edgecolor='black')
conceptbar_ax.text(262,260,'Single Bar from')#,weight='bold')
conceptbar_ax.text(262,245, 'Inner Terrace')#,weight='bold')
conceptbar_ax.text(-0.05, 1.02, 'd.', transform=conceptbar_ax.transAxes, size=14, weight='bold')

ax31 = fig.add_axes([0.13,0.82,0.27,0.145])
ax31.set_ylabel('Depth (m)')
ax31.set_xticklabels('')
thetaMode2 =60
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])
thetaMode1 = np.arange(100,210,10)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   if ii == 0:
        ax31.plot(xinterp,combined,color=colors[ii,:],label='Initial Profile')
   elif ii == len(thetaMode1)-1:
       ax31.plot(xinterp, combined, color=colors[ii, :],label='Final Profile')
   else:
       ax31.plot(xinterp, combined, color=colors[ii, :])
ax31.spines["bottom"].set_linewidth(3)
ax31.spines["top"].set_linewidth(3)
ax31.spines["left"].set_linewidth(3)
ax31.spines["right"].set_linewidth(3)
ax31.spines["bottom"].set_color('blue')
ax31.spines["top"].set_color('blue')
ax31.spines["left"].set_color('blue')
ax31.spines["right"].set_color('blue')
ax31.legend()
ax31.text(-0.18, 1.02, 'a.', transform=ax31.transAxes, size=14, weight='bold')

ax32 = fig.add_axes([0.13,0.66,0.27,0.145])
ax32.set_ylabel('Depth (m)')
ax32.set_xticklabels('')

thetaMode2 =75
eofPred2 = SS[:,0]*np.nanmean(RS[:,0]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,0])
thetaMode1 = np.arange(120,270,10)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,1] * np.nanmean(RS[:, 1]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 1])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax32.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))
ax32.spines["bottom"].set_linewidth(3)
ax32.spines["top"].set_linewidth(3)
ax32.spines["left"].set_linewidth(3)
ax32.spines["right"].set_linewidth(3)
ax32.spines["bottom"].set_color('orange')
ax32.spines["top"].set_color('orange')
ax32.spines["left"].set_color('orange')
ax32.spines["right"].set_color('orange')
ax32.text(-0.18, 1.02, 'b.', transform=ax32.transAxes, size=14, weight='bold')

ax33 = fig.add_axes([0.13,0.5,0.27,0.145])
ax33.set_ylabel('Depth (m)')
ax33.set_xlabel('Cross-shore (m)')
thetaMode2 =70
eofPred2 = SS[:,1]*np.nanmean(RS[:,1]) * np.cos(thetaMode2*np.pi/180 - thetaS[:,1])
thetaMode1 = np.arange(190,290,10)
import matplotlib.cm as cm
colors = cm.rainbow(np.linspace(0, 1, len(thetaMode1)))
for ii in range(len(thetaMode1)):
   eofPred = SS[:,0] * np.nanmean(RS[:, 0]) * np.cos(thetaMode1[ii]*np.pi/180 - thetaS[:, 0])
   combined = np.mean(alllinesS,axis=0)+eofPred+eofPred2
   ax33.plot(xinterp,combined,color=colors[ii,:],label=r'$\phi_1$ = {}'.format(thetaMode1[ii]))

ax33.spines["bottom"].set_linewidth(3)
ax33.spines["top"].set_linewidth(3)
ax33.spines["left"].set_linewidth(3)
ax33.spines["right"].set_linewidth(3)
ax33.spines["bottom"].set_color('green')
ax33.spines["top"].set_color('green')
ax33.spines["left"].set_color('green')
ax33.spines["right"].set_color('green')
ax33.text(-0.18, 1.02, 'c.', transform=ax33.transAxes, size=14, weight='bold')

plt.show()
