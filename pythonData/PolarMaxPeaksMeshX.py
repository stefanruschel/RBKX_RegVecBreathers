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
sweapres =  8*5 # must be divisible by 8 to fit golden ratio for figure
maxPeaksnum = 100

        # EVERYTHING ELSE REMAINS FIXED

# set system parameters 
# -----------------------------------------------------------------------
FileName=   'MaxPeaksMeshX'

# sweap
epsa_vals =  np.concatenate((np.linspace(-0.035, 0.015, sweapres), np.linspace(0.015, -0.035,sweapres)))
epsp_vals =  np.linspace(0.5, 0.1, 5*sweapres/8)
freqPeaks =  np.nan*np.zeros((2*sweapres, 5*sweapres/8))
varPeaks  =  np.nan*np.zeros((2*sweapres, 5*sweapres/8))
epsaMesh  =  np.nan*np.zeros((2*sweapres, 5*sweapres/8))
epspMesh  =  np.nan*np.zeros((2*sweapres, 5*sweapres/8))

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
tau     =   2000.#2000.
 
initPert=   100. # otherwise 100; very large initial conditions 1000 allow for "single" pulses

# set integration parameters 
# -----------------------------------------------------------------------
tstart 	=   2000*tau
tfinal 	=   tstart+maxPeaksnum*tau 
RelTol	=   9.9e-4
AbsTol	=   9.9e-4
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
			'epsa'   :   0.,
			'epsp'   :   0.3,
			'gammaS' :   gammaS,
			'gammaG' :   gammaG,
			'A'      :   A,
			'gammaQ' :   gammaQ,
			'B'      :   B,
			'cp'     :   cp
		} 

thist 		= 	np.linspace(0, tau, 1000)
Ephist 		= 	1.0*initPert*np.power(np.sin(1./2*2.*np.pi/tau*thist),20)*np.exp(2.j*np.pi/tau*thist) 
Emhist 		= 	1.0*initPert*np.power(np.sin(1./2*2.*np.pi/tau*thist),20)*np.exp(2.j*np.pi/tau*thist+1.j*np.pi)
Gphist 		= 	A*np.ones(len(thist))
Qphist 		= 	B*np.ones(len(thist))
Gmhist 		= 	A*np.ones(len(thist))
Qmhist 		= 	B*np.ones(len(thist))
dic 		= 	{'t' : thist, 'Ep': Ephist, 'Em': Emhist, 'Gp': Gphist, 'Qp': Qphist, 'Gm': Gmhist, 'Qm': Qmhist}

# precompute suitable initial segment
dde = dde23(eqns=eqns, params=params)
dde.set_sim_params(tfinal=tfinal)
dde.hist_from_arrays(dic)
dde.run()

spl= 	dde.sample(tstart, tfinal, sample)
Ep = spl['Ep']
Em = spl['Em']
Gp = spl['Gp']
Gm = spl['Gm']
Qp = spl['Qp']
Qm = spl['Qm']

# new initial condition for next run 
locs       =   np.arange(-int(np.floor(tau/sample)),0)
thist 	= 	np.linspace(0, tau, len(locs))
Ephist 	= 	1*(np.real(Ep[locs])+1.j*np.imag(Ep[locs]))
Emhist 	= 	1*(np.real(Em[locs])+1.j*np.imag(Em[locs]))
Gphist 	= 	Gp[locs]
Qphist 	= 	Qp[locs]
Gmhist 	= 	Gm[locs]
Qmhist 	= 	Qm[locs]
dic 	= 	{'t' : thist, 'Ep': Ephist, 'Em': Emhist, 'Gp': Gphist, 'Qp': Qphist, 'Gm': Gmhist, 'Qm': Qmhist}

for j in range(len(epsp_vals)):
    for i in range(len(epsa_vals)):
        params['epsa'] = epsa_vals[i]
        params['epsp'] = epsp_vals[j]
        dde = 	dde23(eqns=eqns, params=params)
        dde.set_sim_params(tfinal=tfinal, AbsTol=AbsTol, RelTol=RelTol)
        dde.hist_from_arrays(dic)
        dde.run()
        
        spl 	= 	dde.sample(tstart, tfinal, sample)
        t = spl['t'] # get sol
        Ep = spl['Ep']
        Em = spl['Em']
        Gp = spl['Gp']
        Gm = spl['Gm']
        Qp = spl['Qp']
        Qm = spl['Qm']

        # new initial condition for next run
        if (max(abs(Ep))>0.1) and (max(abs(Ep))<1000):
            locs        =   np.arange(-int(np.floor(tau/sample)),0)
            thist 		= 	np.linspace(0, tau, len(locs))
            Ephist 		= 	1*(np.real(Ep[locs])+1.j*np.imag(Ep[locs]))
            Emhist 		= 	1*(np.real(Em[locs])+1.j*np.imag(Em[locs]))
            Gphist 		= 	Gp[locs]
            Qphist 		= 	Qp[locs]
            Gmhist 		= 	Gm[locs]
            Qmhist 		= 	Qm[locs]
            dic 		= 	{'t' : thist, 'Ep': Ephist, 'Em': Emhist, 'Gp': Gphist, 'Qp': Qphist, 'Gm': Gmhist, 'Qm': Qmhist}
        else:
            thist 		= 	np.linspace(0, tau, 1000)
            Ephist 		= 	1.0*initPert*np.power(np.sin(1./2*2.*np.pi/tau*thist),20)*np.exp(2.j*np.pi/tau*thist) 
            Emhist 		= 	1.0*initPert*np.power(np.sin(1./2*2.*np.pi/tau*thist),20)*np.exp(2.j*np.pi/tau*thist+1.j*np.pi)
            Gphist 		= 	A*np.ones(len(thist))
            Qphist 		= 	B*np.ones(len(thist))
            Gmhist 		= 	A*np.ones(len(thist))
            Qmhist 		= 	B*np.ones(len(thist))
            dic 		= 	{'t' : thist, 'Ep': Ephist, 'Em': Emhist, 'Gp': Gphist, 'Qp': Qphist, 'Gm': Gmhist, 'Qm': Qmhist}

        # Find peaks of total intensity signal
        '''
        Ip = np.power(np.abs(Ep),2)/np.sqrt(2)
        Im = np.power(np.abs(Em),2)/np.sqrt(2)
        I = Ip + Im
        '''
        
        # Find peaks of x-polarized intensity signal
        I = np.power(np.abs(Ep+Em),2)/2
        
        peaks, _ = find_peaks(I, prominence=1, distance=5000)
        Ipeaks = I[peaks]
        tpeaks = t[peaks]
        
        '''
        spect = np.fft.fft(I)
        freq = np.fft.fftfreq(t.shape[-1])
        spect = spect[freq>0.0001]
        freq = freq[freq>0.0001]
        '''
        
        if len(peaks) > 0:
            varPeaks[i,j]=1.-np.amin(Ipeaks)/np.amax(Ipeaks)
        #freqPeaks[i,j]= 2*np.pi/freq[np.argmax(spect)]/tau   
        epsaMesh[i,j]=epsa_vals[i]
        epspMesh[i,j]=epsp_vals[j]
        print i, j, params['epsa'], params['epsp'], varPeaks[i,j], freqPeaks[i,j]
        
        # checking spectra 
        '''
        #pl.ion()
        fig=pl.figure(0)
        pl.plot(2*np.pi/freq/tau,spect)
        pl.xlim(0,15)
        pl.draw()
        pl.waitforbuttonpress(0) # this will wait for indefinite time
        plt.close(fig)
        '''
        
yon = 'y'#raw_input('Type <y> if initial peak was computed')

#-----------------------------------------------------------------------------
# save 
#-----------------------------------------------------------------------------

if yon =='y':
    #pl.savefig(FileName+".png")
    # save profile
    newfile = open(FileName+".txt", "w")
    for i in range(len(varPeaks[:,0])):
        for l in range(len(varPeaks[0,:])):
            newfile.write(str(varPeaks[i,l])+' ')
        newfile.write("\n")
    newfile.close()
    '''
    newfile = open(FileName+"_freq.txt", "w")
    for i in range(len(freqPeaks[:,0])):
        for l in range(len(freqPeaks[0,:])):
            newfile.write(str(freqPeaks[i,l])+' ')
        newfile.write("\n")
    newfile.close()
    '''
    newfile = open(FileName+"_epsa"+".txt", "w")
    for i in range(len(epsaMesh[:,0])):
        for l in range(len(epsaMesh[0,:])):
            newfile.write(str(epsaMesh[i,l])+' ')
        newfile.write("\n")
    newfile.close()
    
    newfile = open(FileName+"_epsp"+".txt", "w")
    for i in range(len(epspMesh[:,0])):
        for l in range(len(epspMesh[0,:])):
            newfile.write(str(epspMesh[i,l])+' ')
        newfile.write("\n")
    newfile.close()





