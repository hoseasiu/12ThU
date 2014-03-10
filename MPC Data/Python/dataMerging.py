import numpy as np
import matplotlib.pyplot as plt
import math

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
                combinedList.append(lcList[i]+orbList[j][4:5])
                break
            elif lcList[i][2] == orbList[j][2] and lcList[i][3] == orbList[j][3]:
                combinedList.append(lcList[i]+orbList[j][4:5])
                break
    return combinedList

LCparams = np.loadtxt('mpcLightcurveParameters.csv', dtype='str', delimiter=',')
amors = np.loadtxt('MPC_amors.csv', dtype='str', delimiter=',')
apollos = np.loadtxt('MPC_apollos.csv', dtype='str', delimiter=',')
atens = np.loadtxt('MPC_atens.csv', dtype='str', delimiter=',')

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
# split off designation components from the orbital paramters database (Amors)
for i in range(len(amors)):
    number = ''
    name = ''
    year = ''
    letters = ''
    if amors[i][0] != '' and not is_number(amors[i][0]):
        numName = amors[i][0].split()
        number = numName[0]
        name = numName[1]
    year = amors[i][1]
    letters = amors[i][2]
    emoid = amors[i][5]
    H = amors[i][6]
    object = [number, name, year, letters, emoid, H]
    orbList.append(object)

combinedList = checkLists(lcList, orbList, combinedList)
print 'done parsing Amors'

orbList = []
for i in range(len(apollos)):
    number = ''
    name = ''
    year = ''
    letters = ''
    if apollos[i][0] != '' and not is_number(apollos[i][0]):
        numName = apollos[i][0].split()
        number = numName[0]
        name = numName[1]
    year = apollos[i][1].split()[0]
    letters = apollos[i][1].split()[1]
    emoid = apollos[i][4]
    H = apollos[i][5]
    object = [number, name, year, letters, emoid, H]
    orbList.append(object)

combinedList = checkLists(lcList, orbList, combinedList)
print 'done parsing Apollos'

orbList = []
# split off designation components from the orbital paramters database (Amors)
for i in range(len(atens)):
    number = ''
    name = ''
    year = ''
    letters = ''
    if atens[i][0] != '' and not is_number(atens[i][0]):
        numName = atens[i][0].split()
        number = numName[0]
        name = numName[1]
    year = atens[i][1]
    letters = atens[i][2]
    emoid = atens[i][5]
    H = atens[i][6]
    object = [number, name, year, letters, emoid, H]
    orbList.append(object)

combinedList = checkLists(lcList, orbList, combinedList)
print 'done parsing Atens'

print 'all lists added'
################

period = []
emoid = []
for i in range(len(combinedList)):
    if is_number(combinedList[i][4]) and is_number(combinedList[i][6]):
        period.append(float(combinedList[i][4]))
        emoid.append(float(combinedList[i][6]))

plt.close('all')
plt.figure()
plt.plot(emoid,period,'x')
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
##np.savetxt('combinedData.txt', combinedList, delimiter = ',')
##print 'text saved'

plt.show()
