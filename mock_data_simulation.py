import matplotlib
matplotlib.use('Agg')
import numpy as np
import os
import h5py
import operator
import pickle
import random
import matplotlib.pyplot as plt
from astropy.cosmology import Planck15 as cosmo
import numpy
import PhenomA as pa
import LISA as li
import WaveformTools as wt
import os
import pickle as pickle

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}

matplotlib.rc('font', **font)

#reading the mass and redshifts from simulations and meanwhile calucalting the snr in filtering for snr>7.0
f = open('../0template/Klein16_PopIII.dat','r')
lines = f.readlines()
f.close()

#as always start with binning the mass and redshift spaces
Z = []
M1 = []
M2 = []
detections = {}
for line in lines:
    z = float(line.split(' ')[1])
    Z.append(z)
    m1 = float(line.split(' ')[2])
    M1.append(m1)
    m2 = float(line.split(' ')[3])
    M2.append(m2)

    #binary = wt.Binary(m1, m2, z=z)
    #binary.T_merge = T_merge
    #binary.SetFreqBounds(lisa)

    #freqs, X_char = binary.CalcStrain(lisa)
    #snr = binary.CalcSNR(freqs, X_char, lisa)   
    
    if '%s-%s-%s'%(m1,m2,z) in detections.keys():
        detections['%s-%s-%s'%(m1,m2,z)] = detections['%s-%s-%s'%(m1,m2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s-%s-%s'%(m1,m2,z)] = float(line.split(' ')[-1].split('\n')[0])

#zBins = Z (too much)
#zBins.append(20.00)
#zBins = numpy.sort(zBins)
#zBins = numpy.unique(zBins)
zBinsSize = 2.
zBins = numpy.arange(0.0,20.+zBinsSize,zBinsSize)

#unique m1,m2
M1 = numpy.unique(M1)
M2 = numpy.unique(M2)

#binning the M1 and M2
M1[-1] = M1[-1]+1.
M2[-1] = M2[-1]+1.

Mtotal = [M1[0]+M2[0],M1[-1]+M2[-1]]
Mtotallog = np.log10(Mtotal)

M1log = numpy.log10(M1)
M2log = numpy.log10(M2)

M1BinSize = (M1log[-1]-M1log[0])/20.
M1logBins = []
for i in range(21):
    M1logBins.append(M1log[0]+M1BinSize*i)

M2BinSize = (M2log[-1]-M2log[0])/20.
M2logBins = []
for i in range(21):
    M2logBins.append(M2log[0]+M2BinSize*i)

MTBinSize = (Mtotallog[-1]-Mtotallog[0])/8.
MTlogBins = []
MTlogs = []
for i in range(11):
    MTlogBins.append(Mtotallog[0]+MTBinSize*i)

###binning finished, now start actually reading the data
#!!!!! if it is enough to use the previously saved dictionary go to line 348
f = open('../0template/Klein16_PopIII.dat','r')
lines = f.readlines()
f.close()

simType = 1
detections = {}
for line in lines:
    z = float(line.split(' ')[1])
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])
    T_merge = 1.*li.YEAR

    binary = wt.Binary(m1 * pa.TSUN, m2 * pa.TSUN, z=z)
    binary.T_merge = T_merge
    binary.SetFreqBounds(lisa)

    freqs, X_char = binary.CalcStrain(lisa)
    snr = binary.CalcSNR(freqs, X_char, lisa)   
    
    #bn1 = int((numpy.log10(m1)-M1logBins[0])/M1BinSize)
    #bn2 = int((numpy.log10(m2)-M2logBins[0])/M2BinSize)
    bnZ = int((z-zBins[0])/zBinsSize)

    mt = m1 + m2
    bnt = int((np.log10(mt)-MTlogBins[0])/MTBinSize)

    if snr > 7.0:
        if '%s-%s-%s'%(simType,bnt,bnZ) in detections.keys():
            detections['%s-%s-%s'%(simType,bnt,bnZ)] = detections['%s-%s-%s'%(simType,bnt,bnZ)] + float(line.split(' ')[-1].split('\n')[0])
        else:
            detections['%s-%s-%s'%(simType,bnt,bnZ)] = float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s-%s-%s'%(simType,bnt,bnZ)] = 0.0

###
f = open('../0template/Klein16_Q3delays.dat','r')
lines = f.readlines()
f.close()

simType = 2
for line in lines:
    z = float(line.split(' ')[1])
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])
    T_merge = 1.*li.YEAR

    binary = wt.Binary(m1 * pa.TSUN, m2 * pa.TSUN, z=z)
    binary.T_merge = T_merge
    binary.SetFreqBounds(lisa)

    freqs, X_char = binary.CalcStrain(lisa)
    snr = binary.CalcSNR(freqs, X_char, lisa)   
    
    #bn1 = int((numpy.log10(m1)-M1logBins[0])/M1BinSize)
    #bn2 = int((numpy.log10(m2)-M2logBins[0])/M2BinSize)
    bnZ = int((z-zBins[0])/zBinsSize)

    mt = m1 + m2
    bnt = int((np.log10(mt)-MTlogBins[0])/MTBinSize)

    if snr > 7.0:
        if '%s-%s-%s'%(simType,bnt,bnZ) in detections.keys():
            detections['%s-%s-%s'%(simType,bnt,bnZ)] = detections['%s-%s-%s'%(simType,bnt,bnZ)] + float(line.split(' ')[-1].split('\n')[0])
        else:
            detections['%s-%s-%s'%(simType,bnt,bnZ)] = float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s-%s-%s'%(simType,bnt,bnZ)] = 0.0

###
f = open('../0template/Klein16_Q3nodelays.dat','r')
lines = f.readlines()
f.close()

simType = 3
for line in lines:
    z = float(line.split(' ')[1])
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])
    T_merge = 1.*li.YEAR

    binary = wt.Binary(m1 * pa.TSUN, m2 * pa.TSUN, z=z)
    binary.T_merge = T_merge
    binary.SetFreqBounds(lisa)

    freqs, X_char = binary.CalcStrain(lisa)
    snr = binary.CalcSNR(freqs, X_char, lisa)   
    
    #bn1 = int((numpy.log10(m1)-M1logBins[0])/M1BinSize)
    #bn2 = int((numpy.log10(m2)-M2logBins[0])/M2BinSize)
    bnZ = int((z-zBins[0])/zBinsSize)

    mt = m1 + m2
    bnt = int((np.log10(mt)-MTlogBins[0])/MTBinSize)

    if snr > 7.0:
        if '%s-%s-%s'%(simType,bnt,bnZ) in detections.keys():
            detections['%s-%s-%s'%(simType,bnt,bnZ)] = detections['%s-%s-%s'%(simType,bnt,bnZ)] + float(line.split(' ')[-1].split('\n')[0])
        else:
            detections['%s-%s-%s'%(simType,bnt,bnZ)] = float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s-%s-%s'%(simType,bnt,bnZ)] = 0.0

#just to save some time in the future lets save the dictionary!
pickle_out = open("outputs/rick.pickle","wb")
pickle.dump(detections, pickle_out)
pickle_out.close()

#load the dictionary
pickle_in = open("outputs/rick.pickle","rb")
detections = pickle.load(pickle_in)

#setting the bins with no merger equal to zero
for simType in range(1,4):
    for m in range(len(MTlogBins)-1):
        for z in range(len(zBins)-1):
            if '%s-%s-%s'%(simType,m,z) not in detections.keys():
                detections['%s-%s-%s'%(simType,m,z)] = 0.

###finished reading the simulations, calculate the eventRates per redshfit bin per year
def eventRate(N,z1,z2):
    return (1./(z2-z1))*4.*numpy.pi*N*(3.064e-7)*cosmo.comoving_distance((z1+z2)/2.).value**2.

#calculating the observable event rate for each bin - this is per redshift per mass bin per unit time
eventRates = {}
for simType in range(1,4):
    keyNumber = 0
    for m in range(len(MTlogBins)-1):
        for z in range(len(zBins)-1):
            if '%s-%s-%s'%(simType,m,z) in detections.keys():
                eventRates['%s-%s'%(simType,keyNumber)] = eventRate(detections['%s-%s-%s'%(simType,m,z)],zBins[z],zBins[z+1])
            else: 
                eventRates['%s-%s'%(simType,keyNumber)] = 0.
            keyNumber = keyNumber + 1

#calculating the total observable events per year by summing over mass and redshift bins
totalEventRate = {}
for simType in range(1,4):
    keyNumber = 0
    for m in range(len(MTlogBins)-1):
        for z in range(len(zBins)-1):
            if '%s'%simType not in totalEventRate.keys():
                totalEventRate['%s'%simType] = eventRates['%s-%s'%(simType,keyNumber)]*zBinsSize
            else:
                totalEventRate['%s'%simType] = totalEventRate['%s'%simType] + eventRates['%s-%s'%(simType,keyNumber)]*zBinsSize
            keyNumber = keyNumber + 1

#setting up this commulative event to know be able to calculate the number of events comming from each mass/redshit bin
commulativeKeys = {}
for simType in range(1,4):
    for i in range(1+int(len(eventRates.keys())/3)):
        if i == 0: 
            commulativeKeys['%s-%s'%(simType,i)] = 0.
        else:
            commulativeKeys['%s-%s'%(simType,i)] = commulativeKeys['%s-%s'%(simType,i-1)] + eventRates['%s-%s'%(simType,i-1)]*zBinsSize

#simulating the observations - sum of three random draws from our populations
simType = 2
observations = numpy.zeros(int(len(eventRates.keys())/3))
#for o in range(int(3.*totalEventRate['%s'%simType])):
for o in range(1000000):
    seed = numpy.random.rand(1)[0]*totalEventRate['%s'%simType]
    for i in range(int(len(eventRates.keys())/3)-1):
        if commulativeKeys['%s-%s'%(simType,i)]<seed<commulativeKeys['%s-%s'%(simType,i+1)]:
            observations[i] = observations[i] + 1

#define the distribution chi squared
def chiDist(observationF,expectedF):
    y = (observationF - expectedF)**2./expectedF
    return y

simType = 2
simTypeExpObs = numpy.zeros(int(len(eventRates.keys())/3))
for i in range(int(len(eventRates.keys())/3)):
    simTypeExpObs[i] = eventRates['%s-%s'%(simType,i)]*zBinsSize/totalEventRate['%s'%simType]

simTypeComm = numpy.zeros(1+int(len(eventRates.keys())/3))
for i in range(1+int(len(eventRates.keys())/3)):
    simTypeComm[i] = commulativeKeys['%s-%s'%(simType,i)]

for i in range(int(len(eventRates.keys())/3)):
    print('%s ----- %s'%(simTypeExpObs[i],observations[i]))

#calculating the probability distribution from the distribution of individually resolved mergers
#note that normalization is being performed
distProb = {}
sumProbes = 0.0
for simType in range(1,4):
    sumChi = 0.
    for i in range(int(len(eventRates.keys())/3)):
        e = eventRates['%s-%s'%(simType,i)]*zBinsSize/totalEventRate['%s'%simType]
        if e == 0.:
            e = 0.00000000000000001 
        #o = observations[i]/float(int(3.*totalEventRate['%s'%simType]))
        o = observations[i]/1000000.
        chi = chiDist(o,e)
        #print(chi)
        sumChi = sumChi + chi
    probability = prob(sumChi)
    distProb['%s'%simType] = probability
    sumProbes = sumProbes + probability

#normalizing the probaility
nDistProbes = {}
plotProbes = []
for tn in range(1,4):
    nDistProbes['%s'%tn] = distProb['%s'%tn]/sumProbes
    plotProbes.append(nDistProbes['%s'%tn])

#saving the normalized probability
lines = []
lines.append('probability distribution from the distribution only')
count = 1
for tn in templateNumber:
    lines.append('%s-%s'%(tn,nDistProbes['%s'%count]))
    count = count + 1 

f = open('outputs/LDistribution.txt','w')
f.writelines(lines)
f.close()

#plotting the probabilities
templateNumber = ['popIII','Q3delay','Q3noDelay']
plt.bar(templateNumber,plotProbes)
plt.xlabel('simulation')
plt.ylabel('probability')
plt.savefig('plots/LDistribution.eps',format='eps',dpi=1000)
plt.close()

#multiply the strain and distribution
finalProb = {}
summ = 0.
for simType in range(1,4):
    finalProb['%s'%simType] = nDistProbes['%s'%simType]*nProbes['%s'%simType]
    summ = summ + finalProb['%s'%simType]

#normalizing
nFinalProbes = {}
final = []
init = []
for tn in range(1,4):
    nFinalProbes['%s'%tn] = finalProb['%s'%tn]/summ
    final.append(nFinalProbes['%s'%tn])
    init.append(nProbes['%s'%tn])

#plotting
N = 3
ind = np.arange(N)  
width = 0.37       
fig = plt.figure()
ax = fig.add_subplot(111)

yvals = init
rects1 = ax.bar(ind, yvals, width, color='r', label='likelihood (strain)')
zvals = final
rects2 = ax.bar(ind+width, zvals, width, color='b', label ='posterior')

ax.set_ylabel('probability')
ax.set_xticks(ind+width/2)
ax.set_xticklabels( ('popIII', 'Q3delay', 'Q3noDelay') )
ax.legend()

def autolabel(rects):
    for rect in rects:
        h = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*h, '%d'%int(h),
                ha='center', va='bottom')

plt.savefig('plots/LFinal.eps',format='eps',dpi=1000)
plt.close()

#####
#stop
###are we going to use the number of detections?(yes/!no(for now)!) 

#the final prior on number of detections, chi squared
def chiNumber(nObserved,nExpect):
    y = ((nObserved-nExpect)**2.)/(2.)
    return y

#calculating the the probability distribution coming from the number of detections 
numberPrior = []
for simType in range(1,4):
    chi = chiNumber(int(totalEventRate['%s'%1]),totalEventRate['%s'%simType])
    probability = prob(chi)
    numberPrior.append(probability)
