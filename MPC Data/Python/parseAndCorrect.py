from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
import csv

mpcData = np.loadtxt('combinedMPCdata.csv', dtype='str', delimiter=',')

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

''' Bias correction for amplitude, following Binzel and Sauter, 1992
'''
def correctBias(variation):
    # single value cases
    if is_number(variation[0]):      # if a single number is given, keep it, unless >1.0
        num = float(variation[0])
        if num > 1.0:
            return wellSampledEstimate(num)
        else:
            return num
    elif is_number(variation[1]) and not is_number(variation[2]):       # if only a minimum exists
        num = float(variation[1])
        if num > 1.0:
            return wellSampledEstimate(num)
        else:
            return num        
    elif not is_number(variation[1]) and is_number(variation[2]):         # if only a maximum exists
        num = float(variation[2])
        if num > 1.0:
            return wellSampledEstimate(num)
        else:
            return num        

    # both min and max values exist
    elif is_number(variation[1]) and is_number(variation[2]):
        min = float(variation[1])
        max = float(variation[2])
        if max-min <=0.2:
            return (min+max)/2.0
        else:
            return wellSampledEstimate(max)
    
    else:       # nothing parsable was found
        return ''
            
''' Amplitude estimate for a well-sampled variation
    Follows Binzel and Sauter, 1992, p231-232, step 3, equations 2 and 3
    Assumes a = 1, then finds b and c to get A
'''
def wellSampledEstimate(max):
    a = 1.0
    b = math.e**(max/2.5)
    c = b/math.sqrt(2)
    Phi = math.pi/3.0       # 60 radians
    A = 2.5*math.log(a/b)-1.25*math.log((a**2*math.cos(Phi)**2+c**2*math.sin(Phi)**2)/(b**2*math.cos(Phi)**2+c**2*math.sin(Phi)**2))
    if A < 0:       # TODO - is this valid?
        return 0-A
    else:
        return A

parsedData = []
for i in range(1,len(mpcData)):
    row = []
    row = row + mpcData[i][:4].tolist()
    # check if the period column has modifiers it
    if ">" in mpcData[i][4]:        # period given is a minimum
        row = row + ['',mpcData[i][4][2:-1],'']
    elif "<" in mpcData[i][4]:      # period given is a maximum
        row = row + ['','',mpcData[i][4][2:-1]]
    else:
        row = row + [mpcData[i][4],'','']

    # check if the variation column has modifiers in it
    if ">" in mpcData[i][5]:        # variation given is a minimum
        row = row + ['',mpcData[i][5][2:-1],'']
    elif "<" in mpcData[i][5]:      # variation given is a maximum
        row = row + ['','',mpcData[i][5][2:-1]]
    elif "-" in mpcData[i][5]:      # variation given is a range
        split = mpcData[i][5].split('-')
        row = row + ['',split[0],split[1]]
    else:
        row = row + [mpcData[i][5],'','']

    # remove extra quotation marks
    for j in range(len(row)):
        row[j] = row[j].replace('\"','')

        
    row = row + [str(correctBias(row[-3:]))]

    row = row + mpcData[i][6:].tolist()     # add EMOID and H back in

    # remove extra quotation marks
    for j in range(len(row)):
        row[j] = row[j].replace('\"','')

    parsedData.append(row)

# turn corrected variation into an array
correctedVar = []
emoid = []
for i in range(len(parsedData)):
    if is_number(parsedData[i][10]) and is_number(parsedData[i][11]):
        correctedVar = correctedVar + [float(parsedData[i][10])]
        emoid = emoid + [float(parsedData[i][11])]


########################## PLOTTING PORTION ############

##bins = np.linspace(0.0,max(emoid),10000)
bin = np.linspace(0.0,0.05)
print bin
# only considering up to EMOID 0.05
binnedEmoid = np.digitize(emoid,bin)

varBins1 = []
varBins2 = []
varBins3 = []
varBins4 = []
# TODO - restructure this
##varBins4 = []
##varBins = [[] for i in range(max(binnedEmoid))]
for i in range(len(binnedEmoid)):
##    print 'binned ' + str(binnedEmoid[i])
##    print 'correctedVar ' + str(correctedVar[i])
##    varBins[binnedEmoid[i]] = varBins[binnedEmoid[i]] + [correctedVar[i]]
##    print varBins[binnedEmoid[i]]
    if binnedEmoid[i] == 1:
        varBins1 = varBins1 + [correctedVar[i]]
    elif binnedEmoid[i] == 2:
        varBins2 = varBins2 + [correctedVar[i]]
    elif binnedEmoid[i] == 3:
        varBins3 = varBins3 + [correctedVar[i]]
    elif binnedEmoid[i] == 4:
        varBins4 = varBins4 + [correctedVar[i]]
##    else:
##        print 'extra bin!'

# plot bars
plt.close('all')
# TODO - figure out how to use weights
weights = [[np.ones_like(varBins1)/len(varBins1)], [np.ones_like(varBins2)/len(varBins2)], [np.ones_like(varBins3)/len(varBins3)], [np.ones_like(varBins4)/len(varBins4)]]
barLabel = ['MOID <' + str(bin[1]), 'MOID ' + str(bin[1]) + '-' + str(bin[2]), 'MOID ' + str(bin[2]) + '-' + str(bin[3]), 'MOID ' + str(bin[3]) + '-' + str(bin[4])]
n, bins, patches = plt.hist([varBins1, varBins2, varBins3, varBins4], label = barLabel)

##n, bins, patches = plt.hist(varBins1, weights = weights)

##n, bins, patches = plt.hist([varBins1, varBins2, varBins3, varBins3], label = ['MOID <0.1', 'MOID 0.1-0.2', 'MOID 0.2-0.3', 'MOID 0.3-0.4'])
#plt.hist(varBins1, normed = True, stacked = True)  # TODO- is this norming correct? It doesn't seem to sum to 1...

# add guassian line
(mu1, sigma1) = norm.fit(varBins1)
(mu2, sigma2) = norm.fit(varBins2)
(mu3, sigma3) = norm.fit(varBins3)
(mu4, sigma4) = norm.fit(varBins4)
spacing = np.linspace(0,1.8)
y1 = mlab.normpdf(spacing, mu1, sigma1)
y2 = mlab.normpdf(spacing, mu2, sigma2)
y3 = mlab.normpdf(spacing, mu3, sigma3)
y4 = mlab.normpdf(spacing, mu4, sigma4)

##l1 = plt.plot(spacing, y1, 'b')
##l2 = plt.plot(spacing, y2, 'g')
##l3 = plt.plot(spacing, y3, 'r')
##l4 = plt.plot(spacing, y4, 'c')
print 'mu1 = ' + str(mu1) + ' sigma1 = ' + str(sigma1) + ' n1 = ' + str(len(varBins1))
print 'mu2 = ' + str(mu2) + ' sigma2 = ' + str(sigma2) + ' n2 = ' + str(len(varBins2))
print 'mu3 = ' + str(mu3) + ' sigma3 = ' + str(sigma3) + ' n3 = ' + str(len(varBins3))
print 'mu4 = ' + str(mu4) + ' sigma4 = ' + str(sigma4) + ' n4 = ' + str(len(varBins4))

plt.legend()
plt.xlabel('Corrected Light Curve Amplitude (mags)')
plt.ylabel('Number of Objects')
plt.show()

# save this into a text file
out = csv.writer(open("parsedMPCdata.csv","w"), delimiter=',',quoting=csv.QUOTE_ALL, lineterminator = '\n')
out.writerow(['number', 'name', 'year', 'letters', 'period', 'minPeriod', 'maxPeriod', 'variation', 'minVar', 'maxVar', 'correctedVar', 'EMOID', 'H'])
for i in range(len(parsedData)):
    out.writerow(parsedData[i])
print 'text saved'


