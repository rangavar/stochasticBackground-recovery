#the output dictionary named 'usable_mock' gives the redshift, m1, and m2 for a drawn population of mergers in one year
#usable_mock is aranged as usable_mock[key] = [z][m1][m2] where keys are mergers IDs, which are integers going from 0 to total_number_of_mergers_in_one_year
import matplotlib
matplotlib.use('Agg')
import pickle
import numpy as np
from astropy.cosmology import Planck15 as cosmo
import matplotlib.pyplot as plot

#adding the extrapolation ratio: defined based on table 1 in https://arxiv.org/pdf/1511.05581.pdf
extrapolation_rate = 1.89 #1.89 for popIII, 1.72 for Q3_d, 1.98 for Q3_nod

period = 3 #in years
deltaZ = 0.1 #size of the z bins
NM1 = 100 #number of M1 mass bins
NM2 = 100 #number of M2 mass bins
NMT = 100 #number of total-mass bins 

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 14}

#reading the mass and redshifts from simulations and meanwhile calucalting the snr in filtering for snr>7.0
#f = open('/Users/nikos/work/DATA/Klein_et_al_2016/Klein16_PopIII.dat','r') #give the directory of popII or any other catalogue here
f = open('data/Klein16_PopIII.dat','r') #give the directory of popII or any other catalogue here
lines = f.readlines()
f.close()

#calculate the eventRates per redshfit bin per year, multiplied by the extrapolation_rate
def eventRate(N,z1,z2):
    return period * extrapolation_rate * (1./(z2-z1))*4.*np.pi*N*(3.064e-7)*cosmo.comoving_distance((z1+z2)/2.).value**2.

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

    if '%s-%s-%s'%(m1,m2,z) in detections.keys():
        detections['%s-%s-%s'%(m1,m2,z)] = detections['%s-%s-%s'%(m1,m2,z)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s-%s-%s'%(m1,m2,z)] = float(line.split(' ')[-1].split('\n')[0])

#binning z
zBinsSize = deltaZ
zBins = np.arange(0.0,20.+zBinsSize,zBinsSize)

#solving the python precison problem
new_zBins = np.zeros(len(zBins))
for i in range(len(zBins)):
    new_zBins[i] = '%2.2f'%zBins[i]
zBins = new_zBins

#didnt sort masses cause already sorted
M1 = np.unique(M1)
M2 = np.unique(M2)

#binning the M1 and M2
M1[-1] = M1[-1]+1.
M2[-1] = M2[-1]+1.

Mtotal = [M1[0]+M2[0],M1[-1]+M2[-1]]
Mtotallog = np.log10(Mtotal)

M1log = np.log10(M1)
M2log = np.log10(M2)

M1BinSize = (M1log[-1]-M1log[0])/NM1
M1logBins = []
for i in range(NM1+1):
    M1logBins.append(M1log[0]+M1BinSize*i)

M2BinSize = (M2log[-1]-M2log[0])/NM2
M2logBins = []
for i in range(NM2+1):
    M2logBins.append(M2log[0]+M2BinSize*i)

MTBinSize = (Mtotallog[-1]-Mtotallog[0])/NMT
MTlogBins = []
MTlogs = []
for i in range(NMT+1):
    MTlogBins.append(Mtotallog[0]+MTBinSize*i)

#binning finished, now reading the data into bins
###popIII
#f = open('/Users/nikos/work/DATA/Klein_et_al_2016/Klein16_PopIII.dat','r') #give the directory of popII or any other catalogue here
f = open('data/Klein16_PopIII.dat','r') #give the directory of popII or any other catalogue here
lines = f.readlines()
f.close()

detections = {}
for line in lines:
    z = float(line.split(' ')[1])
    bnZ = int((z-zBins[0])/zBinsSize)
    if '%s'%(bnZ) in detections.keys():
        detections['%s'%(bnZ)] = detections['%s'%(bnZ)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        detections['%s'%(bnZ)] = float(line.split(' ')[-1].split('\n')[0])

#event rate
eventRates = {}
for j in detections.keys():
    i = int(j) #control
    eventRates[zBins[i]] = eventRate(detections['%s'%i],zBins[i],zBins[i+1])*deltaZ

#adding the zero bins as well! 
for j in range(len(zBins)):
    if zBins[j] not in eventRates.keys():
        eventRates[zBins[j]] = 0.

#summing the event rates over the bins to get the total number of events
totalNEvenets = 0
for key in eventRates.keys():
    totalNEvenets = totalNEvenets + eventRates[key]

#now we want the distribution of M1 and M2 (or M_total) in each redshift bin
massDistribution = {} #keys = 'zbin-m1bin-m2bin
for line in lines:
    z = float(line.split(' ')[1])
    m1 = float(line.split(' ')[2])
    m2 = float(line.split(' ')[3])

    bnZ = int((z-zBins[0])/zBinsSize)
    bn1 = int((np.log10(m1)-M1logBins[0])/M1BinSize)
    bn2 = int((np.log10(m2)-M2logBins[0])/M2BinSize)

    mt = m1 + m2
    bnt = int((np.log10(mt)-MTlogBins[0])/MTBinSize)

    if '%s-%s-%s'%(bnZ,bn1,bn2) in massDistribution.keys():
        massDistribution['%s-%s-%s'%(bnZ,bn1,bn2)] = massDistribution['%s-%s-%s'%(bnZ,bn1,bn2)] + float(line.split(' ')[-1].split('\n')[0])
    else:
        massDistribution['%s-%s-%s'%(bnZ,bn1,bn2)] = float(line.split(' ')[-1].split('\n')[0])

#calculating the event rate for each of the new redshift_mass_mass bins
massDistribution_EventRate = {}
for key in massDistribution.keys():
    i = int(key.split('-')[0])
    massDistribution_EventRate[key] = eventRate(massDistribution[key],zBins[i],zBins[i+1])*deltaZ

#adding the missing bins (bins with zero event rate)
#making the mass distribution one dimensional for future use (now it is two dimensional: M1-M2)
oneD_massDistribution_EventRate = {}
for c1 in range(NM1):
    for c2 in range(NM2):
        for z in range(len(zBins)-1):
            if '%s-%s-%s' %(z,c1,c2) not in massDistribution_EventRate.keys():
                massDistribution_EventRate['%s-%s-%s' %(z,c1,c2)] = 0
            cT = NM2*c1 + c2 #recal this for conversion in the future
            oneD_massDistribution_EventRate['%s-%s' %(z,cT)] = massDistribution_EventRate['%s-%s-%s' %(z,c1,c2)]

#Cumulative oneD_massDistribution_EventRate to draw random realizations from in the future
cum_oneD_massDistribution_eventRate = {}
for z in range(len(zBins)-1):
    cum_oneD_massDistribution_eventRate['%s-%s'%(z,0)] = 0
    for cT in range(NM1*NM2):
        cum_oneD_massDistribution_eventRate['%s-%s'%(z,cT+1)] = cum_oneD_massDistribution_eventRate['%s-%s'%(z,cT)] + oneD_massDistribution_EventRate['%s-%s'%(z,cT)]

#Cumulative redshift distribution of the events
cum_zDistribution_eventRate = {}
cum_zDistribution_eventRate[0] = 0
for z in range(len(zBins)-1):
    cum_zDistribution_eventRate[z+1] = eventRates[zBins[z]] + cum_zDistribution_eventRate[z]

#drawing a realization from data
mock = {}
for merger in range(int(totalNEvenets)+1):
    #drawing a random redshift
    z_drawn = np.random.random()*totalNEvenets
    
    z_bin = 10000000 #defining redshift bin
    #finding the bin number for drawn redshift
    for cz in range(len(cum_zDistribution_eventRate)-1):
        if cum_zDistribution_eventRate[cz] < z_drawn < cum_zDistribution_eventRate[cz+1]:
            z_bin = cz
    
    #number of events in this redshift bin
    z_events = eventRates[zBins[z_bin]]
    
    #drawing M1 and M2 by drawing from the one_dimensional bin number/ and its one-D bin
    oneDM_drawn = np.random.random()*z_events
    oneDM_bin = 10000000
    for cm in range(NM1*NM2):
        if cum_oneD_massDistribution_eventRate['%s-%s'%(z_bin,cm)] < oneDM_drawn < cum_oneD_massDistribution_eventRate['%s-%s'%(z_bin,cm+1)]:
            oneDM_bin = cm

    #conversion to M1 and M2 bins
    M2_bin = np.mod(oneDM_bin,NM2)
    M1_bin = int(oneDM_bin/NM1)

    mock[merger] = np.zeros(3) #z,M1,M2 - all bin numbers
    mock[merger][0] = int(z_bin)
    mock[merger][1] = int(M1_bin)
    mock[merger][2] = int(M2_bin)

#now converting the results to a more usable format
#note that this is only for one year, and to get the second year you can re-run the code and add the results to the pre-saved dictionary of the first year
usable_mock = {}
for key in mock.keys():
    redshift_bin = int(mock[key][0])
    redshift = 0.5*(zBins[redshift_bin]+zBins[redshift_bin+1])

    mass1_bin = int(mock[key][1])
    mass1 = 10.**( 0.5*(M1logBins[mass1_bin]+M1logBins[mass1_bin+1]) )
    
    mass2_bin = int(mock[key][2])
    mass2 = 10.**( 0.5*(M2logBins[mass2_bin]+M2logBins[mass2_bin+1]) )

    usable_mock[key] = np.zeros(3)
    usable_mock[key][0] = redshift #you are right, this was missing!
    usable_mock[key][1] = mass1
    usable_mock[key][2] = mass2

f = open("popIII_1yr.pkl","wb")
pickle.dump(usable_mock,f)
f.close()
