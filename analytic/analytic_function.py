import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import glob
import emcee

fCut = -2.60
thresh = 6.80
nM1 = 10
nM2 = 10
nZ = 10
nSteps = 10000
skip = 15
endSkip = 6
nMT = nM1 #only for plotting purposes
ndCutoff = 0
popIII_directory = '../1snrThreshold/klein16_popIII/Klein16_PopIII.dat' #give the address to the popIII directory here

#constants
c = 3.064e-7 #lightspeed in Mpc/yr  WA: c in Mpc/yr
G = 4.49E-33 #gravitational constant in Mpc, Mo, yr  WA: Gravitational constant in Megaparsec^3/solar mass/year^2

#functions to calculate min and max frequency of inspiral binary (1psc to ISCO)
def fmin(mt,z):
    y = 100./(2.*np.pi*3.154*(1.+z))*np.sqrt(G*(mt))
    return y
def fmax(mt,z):
    y = 1./(6.*np.pi*np.sqrt(6.)*(1.+z)*3.154e+7)*c**3./(G*(mt))
    return y

#Binning Setup
f = open(popIII_directory,'r')
lines = f.readlines()
f.close()

Z = []
M1 = []
M2 = []
N = {}
for line in lines:
    #redshifts
    z = float(line.split(' ')[1])
    Z.append(z)

    #Mass of each blackhole
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])
    M1.append(m1)
    M2.append(m2)

    #reading the comoving densities
    if '%s-%s-%s' %(m1,m2,z) in N.keys():
        N['%s-%s-%s'%(m1,m2,z)] = N['%s-%s-%s'%(m1,m2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        N['%s-%s-%s'%(m1,m2,z)] = float(line.split(' ')[-1].split('\n')[0])

#redshift bins with extrapolation for the last one
zBins = Z
zBins.append(19.34)
zBins = np.sort(zBins)
zBins = np.unique(zBins)
zBinSize = (zBins[-1]-zBins[0])/float(nZ)
zNewBins = []
for i in range(nZ+1):
    zNewBins.append(zBins[0]+zBinSize*i)
zBins = zNewBins

#binning the M1 and M2
M1.append(50.0)
M1.append(4e+10)
M2.append(50.0)
M2.append(4e+10)
M1 = np.unique(M1)
M2 = np.unique(M2)
M1 = np.sort(M1)
M2 = np.sort(M2)

M1log = np.log10(M1)
M2log = np.log10(M2)

M1BinSize = (M1log[-1]-M1log[0])/nM1
M1logBins = []
for i in range(nM1+1):
    M1logBins.append(M1log[0]+M1BinSize*i)

M2BinSize = (M2log[-1]-M2log[0])/nM2
M2logBins = []
for i in range(nM2+1):
    M2logBins.append(M2log[0]+M2BinSize*i)

MT = []
MT.append(100.0)
MT.append(8e+10)
MTlog = np.log10(MT)
MTBinSize = (MTlog[-1]-MTlog[0])/nMT
MTlogBins = []
for i in range(nMT+1):
    MTlogBins.append(MTlog[0]+MTBinSize*i)

#reading the precalculated snr for each bin
with open('snr_%1.1f_%1.1f_%1.1f.pickle'%(nM1,nM2,nZ), 'rb') as f:
    snrDic = pkl.load(f)

#reading the data into bins
f = open(popIII_directory,'r')
lines = f.readlines()
f.close()

NBinned = {}
for line in lines:
    #redshifts
    z = float(line.split(' ')[1])
    bnZ = int((z-zBins[0])/zBinSize)

    #Mass of each blackhole
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])

    #finding the bin numbers for m1 and m2
    bn1 = int((np.log10(m1)-M1logBins[0])/M1BinSize)
    m1 = 10.**((M1logBins[bn1]+M1logBins[bn1+1])/2.)

    bn2 = int((np.log10(m2)-M2logBins[0])/M2BinSize)
    m2 = 10.**((M2logBins[bn2]+M2logBins[bn2+1])/2.)

    snr = snrDic['%s-%s-%1.8f' %(min(bn1,bn2),max(bn1,bn2),zBins[bnZ])]

    if snr < thresh:
        #reading the comoving densities
        if '%s-%s-%s'%(min(bn1,bn2),max(bn1,bn2),bnZ) in NBinned.keys():
            NBinned['%s-%s-%s'%(min(bn1,bn2),max(bn1,bn2),bnZ)] = NBinned['%s-%s-%s'%(min(bn1,bn2),max(bn1,bn2),bnZ)] + float(line.split(' ')[-1].split('\n')[0])
        else:
            NBinned['%s-%s-%s'%(min(bn1,bn2),max(bn1,bn2),bnZ)] = float(line.split(' ')[-1].split('\n')[0])

#updating the previous functions to be dependent on the input dictionary
def integralM(m1,m2,f,NBinned,FCUT):
    mt = m1+m2
    bn1 = int((np.log10(m1)-M1logBins[0])/M1BinSize)
    bn2 = int((np.log10(m2)-M2logBins[0])/M2BinSize)
    summ = 0.
    for i in range(len(zBins)-1):
        if '%s-%s-%s' %(bn1,bn2,i) in NBinned.keys() and fmax(mt,zBins[i])*(10.**FCUT) < f < fmax(mt,zBins[i]):
            summ = summ + NBinned['%s-%s-%s' %(bn1,bn2,i)]/( (1.+(zBins[i]+zBins[i+1])/2.)**(1./3.) )
    return summ

#defining the chirp mass
def mChirp(m1,m2):
    return (m1*m2)**(3./5.)/(m1+m2)**(1./5.)

#defining the h(f) for each frequency
def h(f,NBinned,FCUT):
    sum = 0.
    count = 0
    for m1log in M1logBins[:nM1]:
        m1 = 10.**(m1log+M1BinSize/2.)
        for m2log in M2logBins[count:nM2]:
            m2 = 10.**(m2log+M2BinSize/2.)
            sum = sum + integralM(m1,m2,f,NBinned,FCUT)*mChirp(m1,m2)**(5./3.)
        count = count + 1
    y = ( 4.*G**(5./3.)/(3.*np.pi**(1./3.)*c**2.) )**(0.5)*np.sqrt(sum)*f**(-2./3.) * (3.17e-8)**(2./3.)
    return f,y*2

freq = np.arange(-5.,-1.,0.05)
freq = 10.**(freq)

lines = []
for f in freq:
    x,y = h(f,NBinned,fCut)
    lines.append('%s %s\n'%(x,y))
