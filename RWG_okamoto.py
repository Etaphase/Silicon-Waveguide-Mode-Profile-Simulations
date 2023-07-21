import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# OKAMOTO --> US
# x --> y (size is 2a = w)
# y --> z (size is 2d = h)
# z --> x
# E^x --> TE 
# E^y --> TM
  
pi = np.pi
WL = 1.55
k0 = 2*pi/WL
h = 0.22

n0 = 1.52
ep0 = n0**2
n1 = 3.48
ep1 = n1**2
Dep = ep1 - ep0

num_modes = 5
num_vals = 50
#w_vals = np.linspace(0.1,3.0,num_vals) # [4] = 0.5
w_vals = np.arange(0.1,3.0,0.02) # [6] = 0.22, [20] = 0.5
num_vals = len(w_vals)
itest = 45
w0 = w_vals[itest]

kyTE_list = np.zeros([num_modes,num_modes,num_vals]) # p,q,w
kzTE_list = np.zeros([num_modes,num_modes,num_vals]) # p,q,w
kyTM_list = np.zeros([num_modes,num_modes,num_vals]) # p,q,w
kzTM_list = np.zeros([num_modes,num_modes,num_vals]) # p,q,w
betaTE = np.zeros([num_modes,num_modes,num_vals]) # p,q,w
betaTM = np.zeros([num_modes,num_modes,num_vals])

for p in range(1,num_modes+1):
    for q in range(1,num_modes+1):
        for wi in range(num_vals):
            w = w_vals[wi]
            # TE
            try:
                kyTE_temp = fsolve(lambda ky: ky*w/2 - (p-1)*(pi/2) - np.arctan(ep1*np.sqrt(Dep*k0**2 - ky**2)/(ep0*ky)),k0)
                kyTE_list[p-1,q-1,wi] = kyTE_temp
                kzTE_temp = fsolve(lambda kz: kz*h/2 - (q-1)*(pi/2) - np.arctan(np.sqrt(Dep*k0**2 - kz**2)/kz),k0)
                kzTE_list[p-1,q-1,wi] = kzTE_temp
                betaTE[p-1,q-1,wi] = np.sqrt(ep1*k0**2 - kyTE_temp**2 - kzTE_temp**2)
            except:
                pass
            
            # TM
            try:
                kyTM_temp = fsolve(lambda ky: ky*w/2 - (p-1)*(pi/2) - np.arctan(np.sqrt(Dep*k0**2 - ky**2)/ky),k0)
                kyTM_list[p-1,q-1,wi] = kyTM_temp
                kzTM_temp = fsolve(lambda kz: kz*h/2 - (q-1)*(pi/2) - np.arctan(ep1*np.sqrt(Dep*k0**2 - kz**2)/(ep0*kz)),k0)
                kzTM_list[p-1,q-1,wi] = kzTM_temp
                betaTM[p-1,q-1,wi] = np.sqrt(ep1*k0**2 - kyTM_temp**2 - kzTM_temp**2)
            except:
                pass

plt.figure(dpi=150)
plt.plot(w_vals,betaTE[0,0,:]/k0,label='TE$_{1,1}$')
plt.plot(w_vals,betaTE[1,0,:]/k0,label='TE$_{2,1}$')
plt.plot(w_vals,betaTE[0,1,:]/k0,label='TE$_{1,2}$')
plt.plot(w_vals,betaTE[1,1,:]/k0,label='TE$_{2,2}$')
plt.xlabel('Waveguide width, $w$ [$\mu$m]')
plt.ylabel('Effective index, $n_{eff}$')
plt.title('TE modes (h = ' + str(h) + ' $\mu$m)')
plt.ylim([0,4])
plt.axvline(x=w0,ymin=0,ymax=4,ls='--',color='k',lw=0.5)
plt.legend(loc='best')
plt.show()

plt.figure(dpi=150)
plt.plot(w_vals,betaTM[0,0,:]/k0,label='TM$_{1,1}$')
plt.plot(w_vals,betaTM[1,0,:]/k0,label='TM$_{2,1}$')
plt.plot(w_vals,betaTM[0,1,:]/k0,label='TM$_{1,2}$')
plt.plot(w_vals,betaTM[1,1,:]/k0,label='TM$_{2,2}$')
plt.xlabel('Waveguide width, $w$ [$\mu$m]')
plt.ylabel('Effective index, $n_{eff}$')
plt.title('TM modes (h = ' + str(h) + ' $\mu$m)')
plt.ylim([0,4])
plt.axvline(x=w0,ymin=0,ymax=4,ls='--',color='k',lw=0.5)
plt.legend(loc='best')
plt.show()

#%% HIGHER ORDER MODES
# plt.figure(dpi=150)
# plt.plot(w_vals,betaTE[1,2,:]/k0,label='TE$_{2,3}$')
# plt.plot(w_vals,betaTE[2,1,:]/k0,label='TE$_{3,2}$')
# plt.plot(w_vals,betaTE[2,2,:]/k0,label='TE$_{3,3}$')
# plt.xlabel('Waveguide width, $w$ [$\mu$m]')
# plt.ylabel('Effective index, $n_{eff}$')
# plt.title('TE modes (h = ' + str(h) + ' $\mu$m)')
# plt.xlim([0.1,5.0])
# plt.ylim([0,4])
# plt.axvline(x=0.5,ymin=0,ymax=4,ls='--',color='k',lw=0.5)
# plt.legend(loc='best')
# plt.show()

# plt.figure(dpi=150)
# plt.plot(w_vals,betaTM[1,2,:]/k0,label='TM$_{2,3}$')
# plt.plot(w_vals,betaTM[2,1,:]/k0,label='TM$_{3,2}$')
# plt.plot(w_vals,betaTM[2,2,:]/k0,label='TM$_{3,3}$')
# plt.xlabel('Waveguide width, $w$ [$\mu$m]')
# plt.ylabel('Effective index, $n_{eff}$')
# plt.title('TM modes (h = ' + str(h) + ' $\mu$m)')
# plt.xlim([0.1,5.0])
# plt.ylim([0,4])
# plt.axvline(x=0.5,ymin=0,ymax=4,ls='--',color='k',lw=0.5)
# plt.legend(loc='best')
# plt.show()

#%% HOW MANY MODES?
print(betaTE[:,:,itest])
valsTE = np.unique(betaTE[:,:,itest])
print('num TE modes = ' + str(len(valsTE[~np.isnan(valsTE)])))

print(betaTM[:,:,itest])
valsTM = np.unique(betaTM[:,:,itest])
print('num TM modes = ' + str(len(valsTM[~np.isnan(valsTM)])))

#%% PLOT THE FIELDS
yz_lim = 3.0
y = np.linspace(-yz_lim,yz_lim,500)
z = np.linspace(-yz_lim,yz_lim,500)
yg,zg = np.meshgrid(y,z)

WG_boundaries = [
    ((yg < w0/2) & (yg > -w0/2)) & ((zg < h/2) & (zg > -h/2)),
    (yg > w0/2) & ((zg < h/2) & (zg > -h/2)),
    ((yg < w0/2) & (yg > -w0/2)) & (zg > h/2),
    (yg < -w0/2) & ((zg < h/2) & (zg > -h/2)),
    ((yg < w0/2) & (yg > -w0/2)) & (zg < -h/2)]

def Hz_TE(p0,q0):
    kyTE0 = kyTE_list[p0-1,q0-1,itest]
    kzTE0 = kzTE_list[p0-1,q0-1,itest]
    gammay0 = np.sqrt(Dep*k0**2 - kyTE0**2)
    gammaz0 = np.sqrt(Dep*k0**2 - kzTE0**2)
    return np.select(WG_boundaries,
        [np.cos(kyTE0*yg-(p0-1)*(pi/2)) * np.cos(kzTE0*zg-(q0-1)*(pi/2)),
        np.cos(kyTE0*w0/2-(p0-1)*(pi/2)) * np.exp(-gammay0*(yg-w0/2)) * np.cos(kzTE0*zg-(q0-1)*(pi/2)),
        np.cos(kyTE0*yg-(p0-1)*(pi/2)) * np.exp(-gammaz0*(zg-h/2)) * np.cos(kzTE0*h/2-(q0-1)*(pi/2)),
        np.cos(kyTE0*w0/2-(p0-1)*(pi/2)) * np.exp(gammay0*(yg+w0/2)) * np.cos(kzTE0*zg-(q0-1)*(pi/2)),
        np.cos(kyTE0*yg-(p0-1)*(pi/2)) * np.exp(gammaz0*(zg+h/2)) * np.cos(kzTE0*h/2-(q0-1)*(pi/2))])

def Hy_TM(p0,q0):
    kyTM0 = kyTM_list[p0-1,q0-1,itest]
    kzTM0 = kzTM_list[p0-1,q0-1,itest]
    gammay0 = np.sqrt(Dep*k0**2 - kyTM0**2)
    gammaz0 = np.sqrt(Dep*k0**2 - kzTM0**2)
    return np.select(WG_boundaries,
        [np.cos(kyTM0*yg-(p0-1)*(pi/2)) * np.cos(kzTM0*zg-(q0-1)*(pi/2)),
        np.cos(kyTM0*w0/2-(p0-1)*(pi/2)) * np.exp(-gammay0*(yg-w0/2)) * np.cos(kzTM0*zg-(q0-1)*(pi/2)),
        np.cos(kyTM0*yg-(p0-1)*(pi/2)) * np.exp(-gammaz0*(zg-h/2)) * np.cos(kzTM0*h/2-(q0-1)*(pi/2)),
        np.cos(kyTM0*w0/2-(p0-1)*(pi/2)) * np.exp(gammay0*(yg+w0/2)) * np.cos(kzTM0*zg-(q0-1)*(pi/2)),
        np.cos(kyTM0*yg-(p0-1)*(pi/2)) * np.exp(gammaz0*(zg+h/2)) * np.cos(kzTM0*h/2-(q0-1)*(pi/2))])

fig, ax = plt.subplots(dpi=150)
cax = ax.pcolormesh(y,z,Hz_TE(1,1),cmap='jet',shading='nearest',vmin=-1,vmax=1)
ax.set_title('Map of $E_y$ for the $TE_{1,1}$ mode')
ax.set_xlabel('y [$\mu$m]')
ax.set_ylabel('z [$\mu$m]')
ax.hlines(y=h/2,xmin=-w0/2,xmax=w0/2,linewidth=0.5,color='k')
ax.hlines(y=-h/2,xmin=-w0/2,xmax=w0/2,linewidth=0.5,color='k')
ax.vlines(x=w0/2,ymin=-h/2,ymax=h/2,linewidth=0.5,color='k')
ax.vlines(x=-w0/2,ymin=-h/2,ymax=h/2,linewidth=0.5,color='k')
fig.colorbar(cax)
plt.xlim([-1,1])
plt.ylim([-1,1])
plt.show()

fig, ax = plt.subplots(dpi=150)
cax = ax.pcolormesh(y,z,Hy_TM(1,2),cmap='jet',shading='nearest',vmin=-1,vmax=1)
ax.set_title('Map of $H_y$ for the $TM_{1,1}$ mode')
ax.set_xlabel('y [$\mu$m]')
ax.set_ylabel('z [$\mu$m]')
ax.hlines(y=h/2,xmin=-w0/2,xmax=w0/2,linewidth=0.5,color='k')
ax.hlines(y=-h/2,xmin=-w0/2,xmax=w0/2,linewidth=0.5,color='k')
ax.vlines(x=w0/2,ymin=-h/2,ymax=h/2,linewidth=0.5,color='k')
ax.vlines(x=-w0/2,ymin=-h/2,ymax=h/2,linewidth=0.5,color='k')
fig.colorbar(cax)
plt.xlim([-1,1])
plt.ylim([-1,1])
plt.show()

# why do regions 4 and 5 need added minus signs while sometimes they don't?
# solve in corner regions too (should be exp decay in both directions)
# solve for other field components using eqn 2.41, 2.43