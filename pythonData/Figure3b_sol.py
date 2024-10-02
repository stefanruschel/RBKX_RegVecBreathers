import numpy as np
import pylab as pl
from pydelay import dde23
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.misc import toimage
from scipy.signal import find_peaks
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.transforms as mtransforms
from colorsys import hls_to_rgb
import matplotlib.colors as mcolors

mpl.rcParams['axes.autolimit_mode'] = 'round_numbers'
mpl.rcParams['axes.xmargin'] = 0
mpl.rcParams['axes.ymargin'] = 0.05

def continuousPhase(complexTimeSeries):
    phi=np.angle(complexTimeSeries)
    for i in range(len(phi)-1):
        if phi[i+1]>phi[i]+np.pi:
            phi[i+1:]=phi[i+1:]-2*np.pi
        if phi[i+1]<phi[i]-np.pi:
            phi[i+1:]=phi[i+1:]+2*np.pi
    return phi
  
#-----------------------------------------------------------------------------
#	 simulation parameters 		
#-----------------------------------------------------------------------------
FileName=   'Figure3b'
epsa    =   0.00
epsp    =   0.355

# set system parameters 
# -----------------------------------------------------------------------
alpha   =   2.
beta    =   0.5
gammaG  =   0.01
gammaQ  =   gammaG
gammaS  =   10.*gammaG # can effectively be set to infty
A       =   2.5
B       =   2.0
a       =   10.*gammaG/gammaQ
cp      =   0.2         # coupling strength usually 0.2
tau     =   2000.
 
initPert=   10. # otherwise 100; very large initial conditions 1000 allow for "single" pulses


# set integration parameters 
# -----------------------------------------------------------------------
tstart 	=   2000.*tau
tfinal 	=   2060.*tau 
RelTol	=   9.9e-8
AbsTol	=   9.9e-8
sample	=   0.1	# solution sample resolution

#-----------------------------------------------------------------------------
# defintition of right hand side, solve
#-----------------------------------------------------------------------------
eqns 	= 	{ 
			'Ep:c': '1./2*(( (1.0+ii*alpha)*Gp -(1.0+ii*beta)*Qp -1.0)*Ep - (epsa+ii*epsp)*Em + cp*Ep(t-tau))',
			'Em:c': '1./2*(( (1.0+ii*alpha)*Gm - (1.0+ii*beta)*Qm -1.0)*Em  - (epsa+ii*epsp)*Ep + cp*Em(t-tau))',
			'Gp'   : 'gammaG*(A -  Gp*(1 + pow(abs(Ep),2)   )) - gammaS*(Gp-Gm)',
            'Qp'   : 'gammaQ*(B -  Qp*(1 + a*pow(abs(Ep),2) )) - gammaS*(Qp-Qm)',
            'Gm'   : 'gammaG*(A -  Gm*(1 + pow(abs(Em),2)   )) - gammaS*(Gm-Gp)',
            'Qm'   : 'gammaQ*(B -  Qm*(1 + a*pow(abs(Em),2) )) - gammaS*(Qm-Qp)'
		}

params = 	{ 
            'a'      :   a,
			'tau'    :   tau,
			'alpha'  :   alpha,
			'beta'   :   beta,
			'epsa'   :   epsa,
			'epsp'   :   epsp,
			'gammaS' :   gammaS,
			'gammaG' :   gammaG,
			'A'      :   A,
			'gammaQ' :   gammaQ,
			'B'      :   B,
			'cp'     :   cp
		} 

dde 		= 	dde23(eqns=eqns, params=params)
thist 		= 	np.linspace(0, tau, 1000)
Ephist 		= 	1.0*initPert*np.power(np.sin(1./2*2.*np.pi/tau*thist),20)*np.exp(0.j*np.pi/tau*thist) #(0.1+0.1j)*np.ones(len(thist))
Emhist 		= 	1.0*initPert*np.power(np.sin(1./2*2.*np.pi/tau*thist),20)*np.exp(2.j*np.pi/tau*thist+1.j*np.pi)
#Emhist 		= 	1.0*initPert*np.power(np.sin(1./2*2.*np.pi/tau*thist),20)*np.exp(3.j*np.pi/tau*thist)
Gphist 		= 	A*np.ones(len(thist))
Qphist 		= 	B*np.ones(len(thist))
Gmhist 		= 	A*np.ones(len(thist))
Qmhist 		= 	B*np.ones(len(thist))
dic 		= 	{'t' : thist, 'Ep': Ephist, 'Em': Emhist, 'Gp': Gphist, 'Qp': Qphist, 'Gm': Gmhist, 'Qm': Qmhist}
dde.set_sim_params(tfinal=tfinal)
dde.hist_from_arrays(dic)
dde.run()

spl 	= 	dde.sample(tstart, tfinal, sample)

t = dde.sol['t'] # get sol
Ep = dde.sol['Ep']
Em = dde.sol['Em']
Gp = dde.sol['Gp']
Gm = dde.sol['Gm']
Qp = dde.sol['Qp']
Qm = dde.sol['Qm']

t	=	spl['t']/tau
Ep	=	spl['Ep']
Em	=	spl['Em']
Ix  =   np.power(np.abs(Ep+Em),2)/2
Iy  =   np.power(np.abs(Ep-Em),2)/2
I   =   np.power(np.abs(Ep),2)+np.power(np.abs(Em),2)

peaks, _ = find_peaks(I, prominence=1, distance=5000)
peaksx, _ = find_peaks(Ix, prominence=1, distance=5000)
peaksy, _ = find_peaks(Iy, prominence=1, distance=5000)
Ipeaks = I[peaks]
tpeaks = t[peaks]

Ixavg = np.average(Ix[peaksx])
Iyavg = np.average(Iy[peaksy])

pl.figure(0)
pl.plot(t, Ix, '-b', label=r'$I_x$')
pl.plot(t, Ixavg+Iyavg -Iy, '-r', label=r'$C-I_y$')
pl.plot(tpeaks, Ipeaks, '*k', label=r'$I$')
pl.legend(loc="upper right")
pl.xlabel(r'$t/\tau$')
plt.axis('tight')
pl.xlim((tstart/tau,tfinal/tau))

pl.figure(1)
pl.plot(t, np.power(np.abs(Ep),2), '-b', label=r'$I_p$')
pl.plot(t, Ixavg+Iyavg -np.power(np.abs(Em),2), '-r', label=r'$C-I_m$')
pl.plot(tpeaks, Ipeaks, '*k', label=r'$I$')
pl.legend(loc="upper right")
pl.xlabel(r'$t/\tau$')
plt.axis('tight')
pl.xlim((tstart/tau,tfinal/tau))

pl.ion()
pl.show()
yon = raw_input('Hit y to save data. Hit anything to discard.')

#-----------------------------------------------------------------------------
# save 
#-----------------------------------------------------------------------------

if yon =='y':
    pl.savefig(FileName+".png")
    # save profile
    newfile = open(FileName+"_sol.txt", "w")
    for i in range(len(t)):
        # write line to output file
        newfile.write(str(t[i])+' '+str(np.real(Ep[i]))+' '+str(np.imag(Ep[i]))+' '+str(np.real(Em[i]))+' '+str(np.imag(Em[i]))+' '+str(Gp[i])+' '+str(Gm[i])+' '+str(Qp[i])+' '+str(Qm[i]))
        newfile.write("\n")
    newfile.close()
    # save parameters



