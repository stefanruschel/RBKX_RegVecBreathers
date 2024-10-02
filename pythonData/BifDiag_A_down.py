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

#-----------------------------------------------------------------------------
#	set simulation parameters 		
#-----------------------------------------------------------------------------
epsp    =   0.355
epsa    =   0.
tau     =   2000. #2000

sweapres =  5*40+1 # should be divisible by 2
maxPeaksnum = 100

# EVERYTHING ELSE REMAINS FIXED

# set system parameters 
# -----------------------------------------------------------------------
FileName=   'MaxPeaks_A'+str(epsp)+'_tau'+str(tau)+'_down'

# sweap
A_vals =  np.linspace(2.5, 2.1,sweapres)
maxPeaks  =  np.nan*np.zeros((sweapres, maxPeaksnum+1))
maxPeaksX  =  np.nan*np.zeros((sweapres, maxPeaksnum+1))

# more or less fixed 
alpha   =   2.
beta    =   0.5
gammaG  =   0.01
gammaQ  =   gammaG
gammaS  =   10.*gammaG # can effectively be set to infty
A       =   2.5
B       =   2.0
a       =   10.*gammaG/gammaQ
cp      =   0.2         # coupling strength usually 0.2
 
initPert=   100. # otherwise 100; very large initial conditions 1000 allow for "single" pulses

# set integration parameters 
# -----------------------------------------------------------------------

# redefined later
tstart 	=   1000*tau # 1000
tfinal 	=   1100*tau # 1100
RelTol	=   9.9e-8
AbsTol	=   9.9e-8
sample	=   .1	# solution sample resolution

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
			'A'      :   A_vals[0],
			'gammaQ' :   gammaQ,
			'B'      :   B,
			'cp'     :   cp
		} 

thist 		= 	np.linspace(0, tau, 1000)
Ephist 		= 	1.0*initPert*np.power(np.sin(1./2*2.*np.pi/tau*thist),20)*np.exp(2.j*np.pi/tau*thist) 
Emhist 		= 	1.0*initPert*np.power(np.sin(1./2*2.*np.pi/tau*thist),20)*np.exp(2.j*np.pi/tau*thist+1j*np.pi)
Gphist 		= 	A*np.ones(len(thist))
Qphist 		= 	B*np.ones(len(thist))
Gmhist 		= 	A*np.ones(len(thist))
Qmhist 		= 	B*np.ones(len(thist))
dic 		= 	{'t' : thist, 'Ep': Ephist, 'Em': Emhist, 'Gp': Gphist, 'Qp': Qphist, 'Gm': Gmhist, 'Qm': Qmhist}

# precompute suitable initial segment
dde = 	dde23(eqns=eqns, params=params)
dde.set_sim_params(tfinal=tfinal)
dde.hist_from_arrays(dic)
dde.run()
    
spl 	= 	dde.sample(tstart, tfinal, sample)
Ep = spl['Ep']
Em = spl['Em']
Gp = spl['Gp']
Gm = spl['Gm']
Qp = spl['Qp']
Qm = spl['Qm']

# new initial condition for next run 
if np.max(np.abs(Ep))>1.:
    locs        =   np.arange(-int(np.floor(tau/sample)),0)
    thist 		= 	np.linspace(0, tau, len(locs))
    Ephist 		= 	1*(np.real(Ep[locs])+1.j*np.imag(Ep[locs]))
    Emhist 		= 	1*(np.real(Em[locs])+1.j*np.imag(Em[locs]))
    Gphist 		= 	Gp[locs]
    Qphist 		= 	Qp[locs]
    Gmhist 		= 	Gm[locs]
    Qmhist 		= 	Qm[locs]
    dic 		= 	{'t' : thist, 'Ep': Ephist, 'Em': Emhist, 'Gp': Gphist, 'Qp': Qphist, 'Gm': Gmhist, 'Qm': Qmhist}


for i in range(len(A_vals)):
    params['A'] = A_vals[i]
    print i, params['A']
    dde = 	dde23(eqns=eqns, params=params)
    dde.set_sim_params(tfinal=tfinal, AbsTol=9.9e-5, RelTol=9.9e-5) #9.9e-5
    dde.hist_from_arrays(dic)
    dde.run()
    
    spl 	= 	dde.sample(tstart, tfinal, sample)
    t = spl['t']/tau # get sol
    Ep = spl['Ep']
    Em = spl['Em']
    Gp = spl['Gp']
    Gm = spl['Gm']
    Qp = spl['Qp']
    Qm = spl['Qm']
    
    # new initial condition for next run 
    locs        =   np.arange(-int(np.floor(tau/sample))-1,-1)
    thist 		= 	np.linspace(0, tau, len(locs))
    Ephist 		= 	1*(np.real(Ep[locs])+1.j*np.imag(Ep[locs]))
    Emhist 		= 	1*(np.real(Em[locs])+1.j*np.imag(Em[locs]))
    Gphist 		= 	Gp[locs]
    Qphist 		= 	Qp[locs]
    Gmhist 		= 	Gm[locs]
    Qmhist 		= 	Qm[locs]
    dic 		= 	{'t' : thist, 'Ep': Ephist, 'Em': Emhist, 'Gp': Gphist, 'Qp': Qphist, 'Gm': Gmhist, 'Qm': Qmhist}
    
    # Find peaks of signal, save and plot
    Ix = np.power(np.abs(Ep+Em),2)/2
    Iy = np.power(np.abs(Ep-Em),2)/2
    I = Iy # note this CHANGE!!! ###########################################################################
    
    # changed to total intensity
    peaks, _ = find_peaks(I, prominence=1, distance=tau*sample/2)
    peaksx, _ = find_peaks(Ix, prominence=1, distance=tau*sample/2)
    #peaks, _ = find_peaks(I, prominence=1)
    
    maxPeaks[i,0]=A_vals[i]
    maxPeaksX[i,0]=A_vals[i]
    
    if np.max(I)<1.:
        peaks = [0]
        I[0] = 0.
        
    if np.max(Ix)<1.:
        peaksx = [0]
        Ix[0] = 0.    
    
    if (np.var(I[peaks]) < 0.01):
        peaks = peaks[0]
    numPeaks=np.size(peaks)
    if numPeaks<=maxPeaksnum:
        maxPeaks[i,1:numPeaks+1]=I[peaks]
    else:
        maxPeaks[i,1:]=I[peaks[0:maxPeaksnum]]
        
    if (np.var(Ix[peaksx]) < 0.01):
        peaksx = peaksx[0]
    numPeaksx=np.size(peaksx)
    if numPeaksx<=maxPeaksnum:
        maxPeaksX[i,1:numPeaksx+1]=Ix[peaksx]
    else:
        maxPeaksX[i,1:]=Ix[peaksx[0:maxPeaksnum]]
        
    # plot Peaks versus A
    pl.figure(1)
    pl.plot(A_vals[i]+0*peaks, I[peaks], 'o', color='tab:orange', alpha=1., markerfacecolor='none')
    pl.plot(A_vals[i]+0*peaksx, Ix[peaksx], 'x', color='tab:blue', alpha=1., markerfacecolor='none')        
    
pl.xlabel(r'$A$')
pl.title(r'$\epsilon_a=$'+str(epsa)+'$, \epsilon_p=$'+str(epsp))
pl.axis('tight')

pl.ion()
pl.show()
yon = raw_input('Type <y> if you want to save figue and data under '+FileName+'.text/.eps?')

#-----------------------------------------------------------------------------
# save 
#-----------------------------------------------------------------------------

if yon =='y':
    pl.savefig(FileName+".png")
    # save profile
    newfile = open(FileName+".txt", "w")
    for i in range(sweapres):
        # write line to output file
        newfile.write(str(maxPeaks[i,0])+' ')
        for l in range(len(maxPeaks[i,:])-1):
            newfile.write(str(maxPeaks[i,l+1])+' ')
        newfile.write("\n")
    newfile.close()
    newfile = open(FileName+"X.txt", "w")
    for i in range(sweapres):
        # write line to output file
        newfile.write(str(maxPeaksX[i,0])+' ')
        for l in range(len(maxPeaksX[i,:])-1):
            newfile.write(str(maxPeaksX[i,l+1])+' ')
        newfile.write("\n")
    newfile.close()




