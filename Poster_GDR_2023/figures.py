"""
* Author : Hugo Jacquet (01/06/2023)

* Produce figures for poster to be exposed at :
	Grand Défi Théoriques pour les sciences du climat. 5,6,7 juin à l'Institut Henri Poincarré, Paris 5e
	
* Input :
	Simulation Aghulas current, average fields (in time, over last 120s)
* Output :
	- wind profiles over 1 box in cold SST and over 1 box in warm SST
	- momentum budget over the same boxes
	- SST, instanteous U at surface, SAR backscatter from Sentinel1 
	- largescale wind (=filtered) with boundary layer height in background
	
* python3 figures.py
"""


import sys
new_path = '/home/jacqhugo/scripts/simu_alex/'

if new_path not in sys.path:
    sys.path.append(new_path)

import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter
from Analysis_functions import interp_mass,interp_mass1D,calcul_EC,moyenne,ABLH
import xarray as xr
# plot related
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib import ticker
# //related
from joblib import Parallel
from joblib import delayed as jb_delayed


path = '/home/jacqhugo/scripts/simu_alex/MEAN_T/NO_DELETE/'
filename = 'MEAN_FIELDS'
filename_mom = 'Momentum_budget'

ds = xr.open_dataset(path+filename+'.nc')
ds_sst =  xr.open_dataset('/home/jacqhugo/WORKDIR/Alex_suite/AGHHR.1.T07h.004.nc')
ds_mom = xr.open_dataset(path+filename_mom+'.nc')

edge = 1 	# MNH HALO
dpi = 300
res = 50	# Horizontal resolution (m)

Z = ds.level[edge:-edge].values
X = ds.ni[edge:-edge].values
Y = ds.nj[edge:-edge].values

box_i = 2*np.array( [[195,255],[180, 240]]) # definition of the boxes
box_j = 2*np.array([[30,80], [360,410]])
boxnames = ['0: cold (=5)','1: warm (=4)']
nboxes = len(box_i)
	
#wind profile
if False:
	# dataset used : temporal average
	figsize = (3,40)
	fig, ax = plt.subplots(nboxes,1,figsize = figsize,constrained_layout=True,dpi=dpi) 
	for k in range(nboxes):
		ax[k].plot(5,0,color='b',label='V') # tricks to have proper legend
		ax[k].plot(ds.UT[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z,color='k',label='U')
		ax2 = ax[k].twiny()
		ax2.plot(ds.VT[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z,color='b',label='V')
		ax2.set_xlim([-0.6,0.6])
		ax[k].set_xlim([3,7.5])
		ax2.vlines(0,0,Z[-1],color='grey',linestyle='--') # mark the 0 for V
		ax[k].set_ylim([0,1000])
		ax[k].set_ylabel('Altitude (m)')
	ax[0].set_title('Wind speed (m.$s^{-1}$)')
	ax[0].text(6.2,900,'0: cold',fontsize=8)
	ax[1].text(6.2,900,'1: warm',fontsize=8)
	ax[0].legend(loc='upper left',fontsize=8)
	plt.show()
	fig.savefig('profiles_UV.png')
	
# budgets moments sur boxes
if False:
	# dataset used : temporal average
	figsize = (6,5)
	fig, ax = plt.subplots(2,nboxes,figsize = figsize,constrained_layout=True,dpi=dpi)
	fig.suptitle(r'Momentum budgets (m.s$^{-2}$)')
	for k in range(nboxes):
		h=1
		# U budget
		ax[0,k].plot(ds_mom.mUadvz[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='k',label='advz',linestyle='dotted')
		ax[0,k].plot(ds_mom.mUadvx[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='k',label='advx')
		ax[0,k].plot(ds_mom.mUadvy[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='k',label='advy',linestyle='--')
		ax[0,k].plot(ds_mom.mUpress[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='b',label='press')
		ax[0,k].plot(ds_mom.mUcorio[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='r',label='corio')
		ax[0,k].plot(ds_mom.mUreynolds[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='g',label='reynolds')
		ax[0,k].plot((ds_mom.mUadvz[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))
				+ds_mom.mUadvx[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))
				+ds_mom.mUadvy[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))
				+ds_mom.mUpress[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))
				+ds_mom.mUcorio[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))
				+ds_mom.mUreynolds[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))),
				Z/h,color='aqua',label='tendency',linestyle='--')
		ax[0,k].set_ylim([0,1000])
		# V budget
		ax[1,k].plot(ds_mom.mVadvz[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='k',label='advz',linestyle='dotted')
		ax[1,k].plot(ds_mom.mVadvx[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='k',label='advx')
		ax[1,k].plot(ds_mom.mVadvy[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='k',label='advy',linestyle='--')
		ax[1,k].plot(ds_mom.mVpress[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='b',label='press')
		ax[1,k].plot(ds_mom.mVcorio[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='r',label='corio')
		ax[1,k].plot(ds_mom.mVreynolds[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2)),Z/h,color='g',label='reynolds')
		ax[1,k].plot((ds_mom.mVadvz[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))
				+ds_mom.mVadvx[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))
				+ds_mom.mVadvy[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))
				+ds_mom.mVpress[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))
				+ds_mom.mVcorio[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))
				+ds_mom.mVreynolds[edge:-edge,box_j[k][0]:box_j[k][1],box_i[k][0]:box_i[k][1]].mean(axis=(1,2))),
				Z/h,color='aqua',label='tendency',linestyle='--')
		ax[1,k].set_ylim([0,1000])
		ax[1,k].set_xlim([-0.5*10**(-3),1*10**(-3)])
		ax[0,k].set_xlim([-0.5*10**(-3),1*10**(-3)])
		ax[1,k].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		ax[1,k].xaxis.major.formatter._useMathText = True
		ax[0,k].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
		ax[0,k].xaxis.major.formatter._useMathText = True
	ax[0,0].set_title('0: cold',loc='left')
	ax[0,1].set_title('1: warm',loc='left')
	ax[0,0].set_title('h = 452 m',loc='right',fontsize=9)
	ax[0,1].set_title('h = 677 m',loc='right',fontsize=9)
	ax[0,0].set_ylabel('U', rotation=0, size='large')
	ax[1,0].set_ylabel('V', rotation=0, size='large')
	ax[0,0].tick_params(axis='both',labelbottom=False,right=True,direction='out')
	ax[0,1].tick_params(axis='both',labelleft=False,labelbottom=False)
	ax[1,1].tick_params(axis='both',labelleft=False,top=True,direction='out')
	ax[1,0].tick_params(axis='both',right='True',top=True,direction='out')
	ax[0,0].legend(fontsize=8)
	fig.savefig('UV_budget_temporalmean.png')
	plt.show()
	
# SST, U surface, SAR	
if False:
	# data used : instanteous field for SST and U, pre processed Sentinel1 data for SAR image	
	figsize = (8,3)
	dpi=200
	fig, ax = plt.subplots(1,3,figsize = figsize,constrained_layout=True,dpi=dpi) # U,V moment
	
	# SST
	levels = np.arange(21,25,0.2)
	s = ax[0].contourf(X/1000,Y/1000,ds_sst.SST[edge:-edge,edge:-edge]-273,cmap = 'bwr',levels=15)
	plt.colorbar(s,ax=ax[0])
	ax[0].set_ylabel('Y (km)',fontsize=8)
	ax[0].set_xlabel('X (km)',fontsize=8)
	ax[0].set_title(r'(a) SST analysis ODYSSEA (°C)',fontsize=9)
	ax[0].set_aspect('equal')
	ax[0].set_ylim([0,43])
	
	# U at surface
	s = ax[1].pcolormesh(X/1000,Y/1000,ds_sst.UT[0,2,edge:-edge,edge:-edge],cmap = 'Greys_r',vmin=3,vmax=6)
	plt.colorbar(s,ax=ax[1])
	ax[1].set_xlabel('X (km)',fontsize=8)
	ax[1].set_title(r'(b) Surface U wind (m.s$^{-1}$)',fontsize=9)
	ax[1].set_aspect('equal')
	ax[1].set_ylim([0,43])
	ax[1].tick_params(axis='both',labelleft=False)
	
	# SAR, partir d'ici c'est recup du code de alex
	y = np.arange(-36.5, -35, 1.5/673)
	x = np.arange(23.5, 26, 2.5/910)
	X,Y = np.meshgrid(x,y)
	X_km = x * 111.3 * (- np.cos(-35.5))
	Y_km = y * 111.3
	sar = np.load('/home/jacqhugo/Satellite/SAR_cap_aiguilles/SAR_odl.npy')[::-1,:]
	xsar, ysar = np.meshgrid(X_km[300:458] - X_km[300], Y_km[300:475] - Y_km[300])
	im = ax[2].pcolormesh(xsar, ysar,sar[300:475, 300:458]/np.amax(sar[300:475, 300:450], axis=(0,1)), cmap='Greys_r')
	cbar = plt.colorbar(im, ax=ax[2])
	ax[2].set_ylim([0,43])
	ax[2].set_title('(c) Normalized SAR backscatter',fontsize=9)
	ax[2].set_aspect('equal')
	ax[2].set_xlabel('X (km)',fontsize=8)
	ax[2].tick_params(axis='both',labelleft=False)
	
	plt.show()
	fig.savefig('SST_Usurface_SAR.png')
	
# largescale wind with zi background	
if False:
	# Large scale wind is the temporal average wind minus mean wind over the domain and then filtered spatially
	# zi is computed as the height at which the potential temperature is 
	#	equal to the integrated potential temperature of lower heights plus 0,25K
	#
	# Parallel computing is done on Y dimension.
	#
	# -----------------------------------------------
	cmap = 'rainbow' # for zi background
	saving = False
	z_surface = 2 	# indice height of surface : 2=5m 10=30m 20=120m 30=230m 50=450m 70=720m 
	figsize = (5,5)
	factorK = 50	# =sigma for gaussian filter
	Karman = 0.4
	# -----------------------------------------------
	# Filtering
	WindU_SS = ds.UT[edge:-edge,edge:-edge,edge:-edge].values
	WindV_SS = ds.VT[edge:-edge,edge:-edge,edge:-edge].values
	mean_U = ds.UT[edge:-edge,edge:-edge,edge:-edge].mean(axis=(1,2)).values
	mean_V = ds.VT[edge:-edge,edge:-edge,edge:-edge].mean(axis=(1,2)).values
	jump = 100 #factorK//2, this is used to plot only some arrows
	sigma = factorK # here is the scale of filtered field (in indices)
	temp1 = WindU_SS[z_surface]-mean_U[z_surface]
	temp2 = WindV_SS[z_surface]-mean_V[z_surface]
	WindU_LS = gaussian_filter(temp1,[sigma,sigma])[jump//2::jump,jump//2::jump]
	WindV_LS = gaussian_filter(temp2,[sigma,sigma])[jump//2::jump,jump//2::jump]
	X_LS = X[jump//2::jump]
	Y_LS = Y[jump//2::jump]
	
	# Parallel processing for ABL height
	def process(y):
		L_zi_x = np.zeros(X.shape[0])
		for x in range(len(X)):
			L_zi_x[x] = ABLH('under0.25KTHTv',Z,ds.THTv[edge:-edge,y,x].values,percent=0.1,u_star=0.)
		return L_zi_x
	YSHAPE =Y.shape[0] # # 5
	Zi = np.zeros((YSHAPE,X.shape[0]))
	print('Start //')
	results = Parallel(n_jobs=24,backend='multiprocessing')(jb_delayed(process)(y) for y in range(YSHAPE)) # threading
	for y in range(YSHAPE):
		Zi[y,:] = results[y]
	print('End //')	
	
	add='_s'+str(sigma)
	print('height :',Z[z_surface])	
	print('U shape, small scale to large scale :',WindU_SS.shape,WindU_LS.shape)
	
	# PLOT -----------
	fig2, ax2 = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi) # large scale field
	# 	Back ground is ABL height
	s = ax2.pcolormesh(X/1000,Y/1000,Zi,cmap=cmap)
	plt.colorbar(s,ax=ax2)
	# 	Supperimposed is 3 SST contours (22.5,23,23.5) °C
	levels = [22.5,23,23.5]
	CS = ax2.contour(X/1000.,Y[300:]/1000.,ds_sst.SST[edge+300:-edge,
					edge:-edge]-273,
					levels=levels,
					colors='grey',linewidths=0.5)
				# here the 300 is to show only interesting lines, cutting artifacts
	fmt={k:str(k)+'°C' for k in CS.levels}
	ax2.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=12)
	# 	horizontal large scale wind vector field
	Q = ax2.quiver(X_LS/1000,Y_LS/1000,WindU_LS,WindV_LS,angles='xy',pivot='mid',headlength=4.5)
	ax2.quiverkey(Q, 0.75, 0.03, 1, r'$\Delta Wind$=$1 m/s$', labelpos='E',coordinates='figure',angle=0) # Reference arrow
	ax2.set_xlabel('X (km)')
	ax2.set_ylabel('Y (km)')
	ax2.set_title(r'$\vec U - \vec Umean $, ABL height (m) at z='+str(np.around(Z[z_surface],2))+'m')
	ax2.set_aspect('equal')
	if saving:
		# plot rectangle where boxes are located
		for i in range(nboxes):
			ax2.add_patch(plt.Rectangle((X[box_i[i,0]]/1000,Y[box_j[i,0]]/1000),
							X[box_i[i,1]]/1000-X[box_i[i,0]]/1000,
								Y[box_j[i,1]]/1000-Y[box_j[i,0]]/1000,
				linewidth=1,edgecolor='k',facecolor='none'))
			ax2.text(X[box_i[i,0]]/1000, Y[box_j[i,0]]/1000, str(i))
			ax2.add_patch(plt.Rectangle((X[box_i[i,0]]/1000,Y[box_j[i,0]]/1000),X[box_i[i,1]]/1000-X[box_i[i,0]]/1000,Y[box_j[i,1]]/1000-Y[box_j[i,0]]/1000,
				linewidth=1,edgecolor='k',facecolor='none'))		
		fig2.savefig('largescale_wind_z'+str(z_surface)+add+'.png')	
	plt.show()	

# coupe en X fixé de WT avec ABLH et en dessous SST et thtv0	
if True:
	ni_coupe = 400 # 400 100
	leftnj = 0
	rightnj = 899
	
	figsize = (7,3)
	saving = True
	method_ablh = 'under0.25KTHTv'
	
	fig, ax = plt.subplots(2,1,figsize = figsize,constrained_layout=True,dpi=dpi,gridspec_kw={'height_ratios': [5,2]})
	
	WindU_SS = ds.UT[edge:-edge,edge:-edge,edge:-edge].values
	WindV_SS = ds.VT[edge:-edge,edge:-edge,edge:-edge].values
	WindW_SS = ds.WT[edge:-edge,edge:-edge,edge:-edge].values
	THTV = ( ds.THT[edge:-edge,edge:-edge,edge:-edge].values * 
					(1+0.61*ds.RVT[edge:-edge,edge:-edge,edge:-edge].values) )		
	zi = np.zeros(Y.shape[0])
	for y in range(len(Y)):
		zi[y] = ABLH(method=method_ablh,Z=Z,col=THTV[:,y,ni_coupe],percent=0.1,u_star=0.)
	x = Y
	left,right = leftnj,rightnj
	ventU = WindU_SS[:,leftnj:rightnj,ni_coupe]
	ventV = WindV_SS[:,leftnj:rightnj,ni_coupe]
	ventW = WindW_SS[:,leftnj:rightnj,ni_coupe]
	surfSST = ds_sst['SST'][edge+leftnj:edge+rightnj,edge+ni_coupe]-273.15

	
	sc = ax[0].pcolormesh(Y[left:right]/1000,Z,ventW,cmap='bwr',vmin=-0.25,vmax=0.25)
	ax[0].plot(Y[left:right]/1000,gaussian_filter(zi[left:right],10),color='grey',linestyle = '--',linewidth=1) #gaussian_filter(zi[left:right],10,5)
	cb = plt.colorbar(sc,ax=ax[0])
	tick_locator = ticker.MaxNLocator(nbins=5)
	cb.locator = tick_locator
	cb.update_ticks()
	ax[0].set_ylabel('W (m/s)')
	ax[0].tick_params('both',labelbottom=False)
	
	ax[1].plot(Y[left:right]/1000,surfSST,color='k',label='SST (°C)')
	ax[1].plot(Y[left:right]/1000,gaussian_filter(THTV[2,leftnj:rightnj,ni_coupe],5)-273,color='grey',label=r'$\theta_{0}$ (°C)')
	ax[1].legend(fontsize=7)
	for axes in ax.flatten():
		axes.set_ylim([0,1000])
	ax[1].set_ylim([21.5,24.5])
	ax[1].set_xlim([Y[left]/1000,Y[right]/1000])
	
	print('Coupe à X=',np.round(res*ni_coupe/1000,2),'km')
	n_coupe = ni_coupe
	plt.show()
	if saving:
		fig.savefig('UVWthtv_X'+str(n_coupe)+'.png')
	
	
	
	
	
	
	
	
	
	
	
