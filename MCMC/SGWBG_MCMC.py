fCut = -2.60
thresh = 6.80
nM1 = 10
nM2 = 10
nZ = 10
nSteps = 10000
skip = 15 #adjusting the frequency range used
endSkip = 6 #same as above
nMT = nM1 #only for plotting purposes
ndCutoff = 1e-10 #minimum number density possible
scatterValue = 0.1 #scatter of the stochastic signal

#defining the prior function
def logPrior(vec):
    for v in vec:
        if v > 0:
            return -1
        if v < -10:
            return -1
        #FCUT prior
        if vec[-1] > fCut+0.5:
            return -1
        if vec[-1] < fCut-0.5:
            return -1
    return 0.

#defining the log-prob function
def logProb(sPrd,sObs,scatter):
    sPrd = sPrd + 1e-40
    sObs = sObs + 1e-40
    sPrd = -1.0*np.log10(sPrd)
    sObs = -1.0*np.log10(sObs)
    y = np.square(sPrd-sObs)/np.square(scatter)
    return y

#defining the total log-proabability function
def totalProb(vec,sObs):
    #checking against the prior
    if logPrior(vec) == -1:
        return np.NINF
    #reading the vector into a dictionary
    vNBinned = {}
    for i in range(len(vec)-1):
        bn1, bn2, bnz = findBin(nonzero_list[i])
        vNBinned['%s-%s-%s'%(bn1,bn2,bnz)] = 10.**vec[i]
    FCUT = vec[-1]
    #calculating the hc for the given set of data & calculating the chi-squared at each step
    scatter = scatterValue
    LOGPROB = 0.
    for fn in range(len(fArray)):
        x,sPrd = h(fArray[fn],vNBinned,FCUT)
        LOGPROB = LOGPROB + logProb(sPrd,sObs[fn],scatter)
    #calculating the (log) probaility
    LOGPROB = -0.5 * LOGPROB
    return LOGPROB

#non-zero list: list of all the bn1,bn2,bnz pairs (in the given order) which have non-zero number density after snr threshold
nonzero_list = []
nonzero_nd = []
nonzero_key = []
count = 0
for bn1 in range(nM1):
    for bn2 in range(nM2):
        for bnz in range(nZ):
            if '%s-%s-%s' %(bn1,bn2,bnz) in NBinned.keys():
                if NBinned['%s-%s-%s' %(bn1,bn2,bnz)] > ndCutoff:
                    nonzero_list.append(count)
                    nonzero_nd.append(NBinned['%s-%s-%s' %(bn1,bn2,bnz)])
                    nonzero_key.append('%s-%s-%s' %(bn1,bn2,bnz))
            count = count + 1

#function to go from count to bn1,bn2,bnz
def findBin(number):
    bnz = np.mod(number,nZ)
    bn2 = np.mod(int(number/nZ),nM2)
    bn1 = int(int(number/nZ)/nM2)
    return bn1, bn2, bnz

#taking the true number densities and producing a fake guess around it
nDim = len(nonzero_nd)+1
nWalkers = nDim*2
nonzero_nd.append(10**(fCut))
p0 = [np.log10(nonzero_nd)+np.append([1*np.random.randn(nDim-1)]+np.ones(nDim-1),np.random.normal(0,0.1)) for i in range(nWalkers)]
#test for if the findBin is working fine
for i in range(len(nonzero_list)):
    print(nonzero_key[i],findBin(nonzero_list[i]))
print('\n\ntotal number of dimensions = %s'%len(nonzero_list))

#setting up the MCMC
filename = 'chain.h5'
backend = emcee.backends.HDFBackend(filename)
backend.reset(nWalkers,nDim)

with Pool() as pool:
    sampler = emcee.EnsembleSampler(nWalkers,nDim,totalProb,args=[hArray],pool=pool,backend=backend)
    pos, prob, stat = sampler.run_mcmc(p0,nSteps,progress=True)
with open('sampler_%s_%s_%s.pickle'%(nM1,nM2,nZ), 'wb') as f:
    pkl.dump(sampler, f)
