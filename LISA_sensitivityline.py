import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

c = 3.*1e8
L = 2.*1e9

def sn_acc(x):
    return (9.*1e-30/(2.*np.pi*x)**4.)*(1.+(1e-4)/x) 

def sn_sn(x):
    return 2.22*1e-23 

def sn_omn(x):
    return 2.65*1e-23

def sn(x):
    return (20./3.)*( (4.*sn_acc(x)+sn_sn(x)+sn_omn(x))/(L**2.) )*( 1.+( f/(0.41*c/(2.*L)) )**2. )

def h(x):
    return np.sqrt(x*sn(x))

freq = np.arange(-5,0.01,0.01)
freq = 10.**freq
strain = []
for f in freq:
    strain.append(h(f))

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size = '18',weight='bold')
plt.rc('axes', axisbelow=True)
plt.rcParams["axes.grid"] = True
plt.rcParams["xtick.major.pad"] = 10.
plt.rcParams["font.family"]='serif'
plt.rcParams["font.serif"][0] = 'times'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['xtick.major.width'] = 0.8
plt.rcParams['ytick.major.width'] = 0.8
plt.figure(figsize=(8, 8))
plt.plot(freq, strain, c='k', linestyle='-', linewidth = 2, label='Dayal et al. 2018',zorder=1,alpha =0.75)

def s_i(f):
    return 5.76*1e-48*( 1+(f_a/f)**2 )

def s_ii(f):
    return 3.6*1e-41

def r(f):
    return 1+(f/f_b)**2

def ksn(f):
    return ((10./3.)*( s_i(f)/(2.*np.pi*f)**4 + s_ii(f) )*r(f)+s_c(f))

def s_c(f):
    return A/2. * np.e**(-(f/f_1)**alpha) * f**(-7./3.) * ( 1+np.tanh((f_knee-f)/f_2) )

def kh(f):
    return np.sqrt(f*ksn(f))

A = 1.28*1e-44
f_a = 0.4*1e-3
f_b = 25*1e-3
alpha = 1.63

a_1 = -0.224
b_1 = -2.704
a_k = -0.361
b_k = -2.378
t_obs = 4
f_1 = 10.**( a_1*np.log10(t_obs)+b_1 )
f_2 = 10.**(-3.318)
f_knee = 10.**( a_k*np.log10(t_obs)+b_k )

freq = np.arange(-5,0.01,0.01) 
freq = 10.**freq
strain = []
for f in freq:
    strain.append(kh(f))

plt.plot(freq, strain, c='maroon', linestyle='-', linewidth = 3, label='Bonetti et al. 2020',zorder=1,alpha =0.75)
plt.xlabel('Frequency [Hz]', fontsize='22',weight = 'bold')
plt.ylabel('Characteristic Strain', fontsize='22')
plt.yscale('log')
plt.xscale('log')
plt.xlim(1e-5,1)
plt.ylim(1e-22,1e-16)
plt.legend(frameon=True, loc='upper right')
plt.savefig('sensitivity.png',format='png',dpi=400)
plt.close()
