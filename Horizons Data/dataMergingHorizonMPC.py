import numpy as np
import matplotlib.pyplot as plt
import math
import csv

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def checkLists(lcList, orbList, combinedList):
    for i in range(len(lcList)):
        for j in range(len(orbList)):
            # check number and name
            if lcList[i][0] == orbList[j][0] and lcList[i][0] == orbList[j][0]:
                combinedList.append(lcList[i]+orbList[j][2:])
                break
##            elif lcList[i][2] == orbList[j][2] and lcList[i][3] == orbList[j][3]:
##                combinedList.append(lcList[i]+orbList[j][4:])
##                break
    return combinedList

LCparams = np.loadtxt('mpcLightcurveParameters.csv', dtype='str', delimiter=',')
horizons = np.loadtxt('horizons_data.csv', dtype='str', delimiter=',')

lcList = []

# split off designation components from the light curve database
for i in range(len(LCparams)):
    desig = str.split(LCparams[i][0])
    period = LCparams[i][2]
    variation = LCparams[i][3]
    number = ''
    name = ''
    year = ''
    letters = ''
    foundLetters = False
    for j in range(len(desig)):
        # check for an object number (in parentheses)
        if '(' in desig[j] and ')' in desig[j]:
            number = desig[j]
        # check for a year
        elif is_number(desig[j]):
            year = desig[j]
            letters = desig[j+1]
            foundLetters = True
        elif foundLetters == False and '(' not in desig[j]:
            name = desig[j]
    object = [number, name, year, letters, period, variation]
    lcList.append(object)

print 'done parsing LCparams'

combinedList = []

orbList = []
# split off designation components from the orbital paramters database (horizons)
for i in range(1,len(horizons)):
    number = ''
    name = ''
    year = ''
    letters = ''
    if is_number(horizons[i][1]):
        number = horizons[i][1]
        name = horizons[i][2]
        emoid = horizons[i][9]        
    object = [number, name, emoid]
    orbList.append(object)

combinedList = checkLists(lcList, orbList, combinedList)

print 'all lists added'
################

period = []
emoid = []
for i in range(len(combinedList)):
    if is_number(combinedList[i][4]) and is_number(combinedList[i][6]):
        period.append(float(combinedList[i][4]))
        emoid.append(float(combinedList[i][6]))
print "lists combined"
print period
print emoid

plt.close('all')
plt.figure()
plt.plot(emoid,period,'x')
period = np.array(period)
emoid = np.array(emoid)
plt.xlabel('MOID (AU)')
plt.ylabel('Period (h)')
plt.title('Period vs. MOID')
plt.ylim([0, 75])

plt.figure()
plt.hist(emoid, bins = 25)

averagePeriods = []
binnedEmoids = []
for i in range(int(math.ceil(max(emoid))/0.05)):
    sectionStart = i*0.05
    sectionEnd = (i+1)*0.05
    sectionPeriodList = []
    for j in range(len(emoid)):
        if emoid[j] >= sectionStart and emoid[j] < sectionEnd:
            sectionPeriodList.append(period[j])
    if len(sectionPeriodList)>0:
        binnedEmoids.append((sectionStart+sectionEnd)/2.0)
        averagePeriods.append(sum(sectionPeriodList)/len(sectionPeriodList))
    
plt.figure()
plt.plot(binnedEmoids,averagePeriods,'x')

# write to CSV
out = csv.writer(open("combinedMPCdata.csv","w"), delimiter=',',quoting=csv.QUOTE_ALL, lineterminator = '\n')

out.writerow(['number', 'name', 'year', 'letters', 'period', 'variation', 'EMOID', 'H'])
for i in range(len(combinedList)):
    out.writerow(combinedList[i])
print 'text saved'

plt.show()
