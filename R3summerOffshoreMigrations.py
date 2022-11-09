import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import datetime as DT
import os
import matplotlib.patches as mpatches
import pickle
from datetime import datetime
from datetime import timedelta


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




# not 12-15
# not 44-50
# not 97-100
# not 133-136
# not 150-153
# not 168-171
# not 185-190
# kind of but not 209-212
# not 224-232
# not 270-278
# not 306-209
# not 324-327
# not 357-358
# not 431-441
# not 458-461
# not 473-475

# maybe 72-74
# maybe 118-120
# maybe 201-203
# maybe 253-257
# maybe 409-412
# maybe 520-522

import matplotlib.cm as cm
plt.figure()

p1 = plt.subplot2grid((2,3),(0,0),rowspan=1,colspan=1)
si = 72
ei = 75
colors = cm.hsv(np.linspace(0, 1, (ei-si)+1))
for ff in range((ei-si)):
    qq = ff + si
    p1.plot(xinterp,alllines[qq,:],color=colors[ff,:],label=time[qq-1])
p1.legend()

p2 = plt.subplot2grid((2,3),(0,1),rowspan=1,colspan=1)
si = 118
ei = 121
colors = cm.hsv(np.linspace(0, 1, (ei-si)+1))
for ff in range((ei-si)):
    qq = ff + si
    p2.plot(xinterp,alllines[qq,:],color=colors[ff,:],label=time[qq-1])
p2.legend()

p3 = plt.subplot2grid((2,3),(0,2),rowspan=1,colspan=1)
si = 204
ei = 208
colors = cm.hsv(np.linspace(0, 1, (ei-si)+1))
for ff in range((ei-si)):
    qq = ff + si
    p3.plot(xinterp,alllines[qq,:],color=colors[ff,:],label=time[qq])
p3.legend()

p4 = plt.subplot2grid((2,3),(1,0),rowspan=1,colspan=1)
si = 253
ei = 258
colors = cm.hsv(np.linspace(0, 1, (ei-si)+1))
for ff in range((ei-si)):
    qq = ff + si
    p4.plot(xinterp,alllines[qq,:],color=colors[ff,:],label=time[qq-1])
p4.legend()

p5 = plt.subplot2grid((2,3),(1,1),rowspan=1,colspan=1)
si = 409
ei = 413
colors = cm.hsv(np.linspace(0, 1, (ei-si)+1))
for ff in range((ei-si)):
    qq = ff + si
    p5.plot(xinterp,alllines[qq,:],color=colors[ff,:],label=time[qq-1])
p5.legend()

p6 = plt.subplot2grid((2,3),(1,2),rowspan=1,colspan=1)
si = 596
ei = 600
colors = cm.hsv(np.linspace(0, 1, (ei-si)+1))
for ff in range((ei-si)):
    qq = ff + si
    p6.plot(xinterp,alllines[qq,:],color=colors[ff,:],label=time[qq-1])
p6.legend()
#
# p7 = plt.subplot2grid((3,3),(2,0),rowspan=1,colspan=1)
# si = 458
# ei = 462
# colors = cm.hsv(np.linspace(0, 1, (ei-si)+1))
# for ff in range((ei-si)):
#     qq = ff + si
#     p7.plot(xinterp,alllines[qq,:],color=colors[ff,:],label=qq)
# p7.legend()
#
# p8 = plt.subplot2grid((3,3),(2,1),rowspan=1,colspan=1)
# si = 473
# ei = 476
# colors = cm.hsv(np.linspace(0, 1, (ei-si)+1))
# for ff in range((ei-si)):
#     qq = ff + si
#     p8.plot(xinterp,alllines[qq,:],color=colors[ff,:],label=qq)
# p8.legend()
#
# p9 = plt.subplot2grid((3,3),(2,2),rowspan=1,colspan=1)
# si = 520
# ei = 523
# colors = cm.hsv(np.linspace(0, 1, (ei-si)+1))
# for ff in range((ei-si)):
#     qq = ff + si
#     p9.plot(xinterp,alllines[qq,:],color=colors[ff,:],label=qq)
# p9.legend()



