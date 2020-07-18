import matplotlib
matplotlib.use('Agg')
import numpy as np
import emcee
import h5py
import matplotlib.pyplot as plt
import pickle
#preferrebly on the same order as the nMComb below
nM1 = 25
nM2 = 25
#shared bin numbers for m1 and m2, n(n+1)/2 steps rather than n2
#takes care of overestimation by factor 2
nMComb = 25
nMT = 5
ndim = nMT
nwak = 30
nsteps = 500
freq = np.arange(-4.0,-0.99,0.01)
freq = 10.**(freq)
zBinsSize = 0.5
error = 0.02

###defining the strain functions
c = 3.064e-7 #lightspeed in Mpc/yr  WA: c in Mpc/yr
G = 4.49E-33 #gravitational constant in Mpc, Mo, yr  WA: Gravitational constant in Megaparsec^3/solar mass/year^2

#functions to calculate min and max frequency of inspiral binary (1psc to ISCO)
def fmin(mt,z):
    y = 100./(2.*np.pi*3.154*(1.+z))*np.sqrt(G*(mt))
    return y
def fmax(mt,z):
    y = 1./(6.*np.pi*np.sqrt(6.)*(1.+z)*3.154e+7)*c**3./(G*(mt))
    return y

#defining the inner intgral
def integralM(mt,f,NBinned):
    #bn1 = int((np.log10(m1)-M1logBins[0])/M1BinSize)
    #bn2 = int((np.log10(m2)-M2logBins[0])/M2BinSize)
    bnt = int((np.log10(mt)-MTlogBins[0])/MTBinSize)
    summ = 0.
    for i in range(len(zBins)-1):
        if '%s' %(bnt) in NBinned.keys() and fmin(mt,zBins[i]) < f < fmax(mt,zBins[i]):
            summ = summ + (zBins[i+1]-zBins[i])*NBinned['%s'%(bnt)]/( (1.+(zBins[i]+zBins[i+1])/2.)**(1./3.) * zLength)
    return summ

#defining the chirp mass
def mChirp(m1,m2):
    return (m1*m2)**(3./5.)/(m1+m2)**(1./5.)

#defining the h(f) for each frequency
def h(f,NBinned):
    sum = 0.
    #i = 0
    for x in range(len(MComblogBins)):
        m1log = MComblogBins[x]
        m1 = 10.**(m1log+MCombBinSize/2.)
        for y in range(x,len(MComblogBins)):
            m2log = MComblogBins[y]
            m2 = 10.**(m2log+MCombBinSize/2.)
            totalMass = m1+m2
            sum = sum + integralM(totalMass,f,NBinned)*mChirp(m1,m2)**(5./3.)
            #print('%s-%s'%(i,sum))
        #i = i+1
    y = ( 4.*G**(5./3.)/(3.*np.pi**(1./3.)*c**2.) )**(0.5)*np.sqrt(sum)*f**(-2./3.) * (3.17e-8)**(2./3.)
    return f,y

###read the data and save it on simplified total mass bins + marginalize over z
f = open('Klein16_PopIII.dat','r')
lines = f.readlines()
f.close()

Z = []
M1 = []
M2 = []
N = {}
#maxFreq = {}
#minFreq = {}
for line in lines:
    #redshifts
    z = float(line.split(' ')[1])
    Z.append(z)

    #Mass of each blackhole
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])
    M1.append(m1)
    M2.append(m2)

    #calculate min and mass frequency for that specific total mass
    #maxFreq['%s-%s-%s'%(m1,m2,z)] = fmax(m1,m2,z)
    #minFreq['%s-%s-%s'%(m1,m2,z)] = fmin(m1,m2,z)

    #reading the comoving densities
    if '%s-%s-%s' %(m1,m2,z) in N.keys():
        N['%s-%s-%s'%(m1,m2,z)] = N['%s-%s-%s'%(m1,m2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        N['%s-%s-%s'%(m1,m2,z)] = float(line.split(' ')[-1].split('\n')[0])

#redshift bins with extrapolation for the last one
zBins = np.arange(0.0,20.+zBinsSize,zBinsSize)
zLength = zBins[-1]-zBins[0]

#unique m1,m2
M1 = np.unique(M1)
M2 = np.unique(M2)

#binning the M1 and M2
M1[-1] = M1[-1]+1.
M2[-1] = M2[-1]+1.

#combined bin production
MCombMin = min(M1[0],M2[0])
MCombMax = max(M1[-1],M2[-1])
MCombMin = np.log10(MCombMin)
MCombMax = np.log10(MCombMax)

Mtotal = [M1[0]+M2[0],M1[-1]+M2[-1]]
Mtotallog = np.log10(Mtotal)

M1log = np.log10(M1)
M2log = np.log10(M2)

M1BinSize = (M1log[-1]-M1log[0])/float(nM1)
M1logBins = []
for i in range(nM1+1):
    M1logBins.append(M1log[0]+M1BinSize*i)

M2BinSize = (M2log[-1]-M2log[0])/float(nM2)
M2logBins = []
for i in range(nM2+1):
    M2logBins.append(M2log[0]+M2BinSize*i)

MCombBinSize = (MCombMax-MCombMin)/float(nMComb)
MComblogBins = []
for i in range(nMComb+1):
    MComblogBins.append(MCombMin+MCombBinSize*i)

MTBinSize = (Mtotallog[-1]-Mtotallog[0])/nMT
MTlogBins = []
MTlogs = []
for i in range(nMT+1):
    MTlogBins.append(Mtotallog[0]+MTBinSize*i)

f = open('Klein16_PopIII.dat','r')
lines = f.readlines()
f.close()

NBinned = {}
#maxFreq = {}
#minFreq = {}
for line in lines:
    #redshifts
    z = float(line.split(' ')[1])

    #Mass of each blackhole
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])

    bn1 = int((np.log10(m1)-MComblogBins[0])/MCombBinSize)
    bn2 = int((np.log10(m2)-MComblogBins[0])/MCombBinSize)

    mt = m1 + m2
    bnt = int((np.log10(mt)-MTlogBins[0])/MTBinSize)

    #reading the comoving densities
    if '%s'%(bnt) in NBinned.keys():
        NBinned['%s'%(bnt)] = NBinned['%s'%(bnt)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        NBinned['%s'%(bnt)] = float(line.split(' ')[-1].split('\n')[0])

#saving the dictionary for tomorrow morning
pickle_out = open("NBinned.pickle","wb")
pickle.dump(NBinned, pickle_out)
pickle_out.close()

lines = []
for f in freq:
    x,y = h(f,NBinned)
    lines.append('%s %s\n'%(x,y))

f = open('output/strain_simplified.txt','w')
f.writelines(lines)
f.close()

mt = []
for key in NBinned.keys():
    bnt = int(key)
    mt.append(10.**((MTlogBins[bnt]+MTlogBins[bnt+1])/2.))
mt = np.sort(mt)
value = []
for m in mt:
    bnt = int((np.log10(m)-MTlogBins[0])/MTBinSize)
    value.append(NBinned['%s'%bnt])

mt = np.log10(mt)
value = np.log10(value)
plt.plot(mt,value)
plt.savefig('output/mPlot.eps',format='eps',dpi=1000)
plt.close()

#load the dictionary
pickle_in = open("NBinned.pickle","rb")
NBinned = pickle.load(pickle_in)

#real strain(f): observational data here (noisless)
m = open('output/strain_simplified.txt','r')
lines = m.readlines()
m.close()

data = np.zeros(len(freq))
freq = []
counter = 0
for line in lines:
    freq.append(float(line.split(' ')[0]))

    #noise
    workStrain = np.log10(float(line.split(' ')[1]))
    workStrain = np.random.normal(workStrain,-workStrain*error)
    #workStrain = workStrain + workStrain*error
    #strain[count] = 10.**(workStrain)

    data[counter] = 10.**(workStrain)
    #!data[counter] = float(line.split(' ')[1])
    counter = counter + 1

#this is in a sense my whole mcmc function 
def log_prob(par,truth):
    #counter = 0
    for w in par:
        if w < -6.0 or w > 1.:return np.NINF
        #if w > 3.0: return numpy.NINF
        #if original[counter] != 0.:
        #    if w>np.log10(original[counter])+5. or w<np.log10(original[counter])-5.:return np.NINF
        #counter = counter + 1
        #print(counter,original[counter-1])

    MBinned = {}
    for i in range(nMT):
        MBinned['%s'%(i)] = 10.**(par[i])
        #MBinned['%s-%s'%(i,j)] = 10.**(par[i*nM1+j])

    strain = np.zeros(len(freq))
    counter = 0
    for f in freq:
        x,y = h(f,MBinned)
        strain[counter] = y
        counter = counter + 1

    diff = (strain - truth) * 10.**20.
    #print(diff)
    #print(-0.5 * np.dot(diff, diff))
    return -0.5 * np.dot(diff, diff)

#the correct answer
original = np.zeros(nMT)
for i in range(nMT):
    if '%s'%(i) in NBinned.keys():
        original[i] = NBinned['%s'%(i)]

#starting value all random
par0 = np.random.rand(nwak,ndim)

filename = "chain.h5"
backend = emcee.backends.HDFBackend(filename)
backend.reset(nwak, ndim)

from multiprocessing import Pool
#multiprocessing.cpu_count()

with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwak,ndim,log_prob,args=[data],pool=pool,backend=backend)
    pos, prob, stat = sampler.run_mcmc(par0,nsteps,progress=True)

np.save('pos.npy', pos)

#state = sampler.run_mcmc(par0,100)

samples = sampler.get_chain(flat=True)

walkerValues = []
error = []
for i in range(nMT):
    sample = samples[:,i]
    #sample = np.sort(sample)
    peak = sample[int(nwak*50./100.)]
    down = sample[int(nwak*16./100.)]
    up = sample[int(nwak*84./100.)]
    diff = np.abs(up-down)
    error.append(diff)
    walkerValues.append(peak)

plt.errorbar(mt,walkerValues,yerr=error,marker='.',linestyle="None",markerfacecolor="red")
plt.xlabel('total mass($M_{sun}$)')
plt.ylabel('number density ($Mpc^{-3}$)')
plt.savefig('totalmas.pdf')

plt.close()
for i in range(nMT):
    plt.hist((samples[:,i]),nwak,color="k",histtype="step")
    if original[i] != 0.: plt.axvline(x=np.log10(original[i]))
    #plt.xscale('log')
    plt.savefig('output/%s_parameter.pdf'%(i))
    plt.close()                               
