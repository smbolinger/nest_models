# sudo vim -o file1 file2 [open 2 files] 
# BLAH
# NOTE 5/16/25 - The percent bias responds more like I would expect when I use
#                the actual calculated DSR, not the assigned DSR (0.93 or 0.95)
#                BUT I still don't know why the calculated DSR is consistently low.

from datetime import datetime
import decimal
from decimal import Decimal
from itertools import product
import itertools
import numpy as np 
from os.path import exists
import os
from pathlib import Path
from scipy import optimize
import scipy.stats as stats
import sys

# -----------------------------------------------------------------------------
#  SETTINGS 
# -----------------------------------------------------------------------------
rng = np.random.default_rng(seed=102891) 
args = sys.argv

debug = False
if len(args) > 1:
    debug = args[1] == "debugTrue"

# -----------------------------------------------------------------------------
#  HELPER FUNCTIONS
# -----------------------------------------------------------------------------
# from https://stackoverflow.com/questions/13852700/create-file-but-if-name-exists-add-number
def uniquify(path):
    filename, extension = os.path.splitext(path)
    counter = 1
    # print("filename, extension:", filename, extension)

    while os.path.exists(path):
        path = filename + " (" + str(counter) + ")" + extension
        counter += 1

    return path

# -----------------------------------------------------------------------------
def searchSorted2(a, b):
    #out = np.zeros(a.shape)
    out = np.zeros((a.shape[0], len(b)))
    for i in range(len(a)):
        #out[i] = np.searchsorted(a[i], b[i])
        #print("sorted search of\n", b, "within\n", a[i])
        # if debug: print(">> sorted search of", b, "within", a[i])
        out[i] = np.searchsorted(a[i], b)
        if debug: print("sorted search of\n", b, "within\n", a[i], ":\n", out, out.shape)
        #print("index positions:", out, out.shape)
        # shouldn't the output have the shape of b?
    #print(">> index positions:\n", out, out.shape)
    return(out)

# -----------------------------------------------------------------------------
#   OPTIMIZER PARAMETERS: 
# -----------------------------------------------------------------------------
nruns    = 1 # number of times to run  the optimizer 
#nreps   = 1000 # number of replicates of simulated nest data
nreps   = 500
# nreps    = 5
# nreps    = 1
makeNest = True
num_out  = 21 # number zof output parameters
#dir_name = datetime.now().strftime("%m%d%Y")
# -----------------------------------------------------------------------------
#   NEST MODEL PARAMETERS:
# -----------------------------------------------------------------------------
# NOTE the main problem is that my fate-masking variable (storm activity) also 
#      leads to certain nest fates
# NOTE 2: how many varying params is a reasonable number?
breedingDays   = [150] 
discProb       = [0.7] 
# discProb       = [1] 
numNests       = [500,1000]
# numNests       = [500]
# numNests       = [10]
uncertaintyDays = [1]
print(
        "STATIC PARAMS: breeding season length:", breedingDays,
        "; discovery probability:", discProb
        # "; number of nests:", numNests
        # "how long nest fate is discoverable (in days):", uncertaintyDays
        )
# only really need to be in lists if they are part of the combinations below
# probSurv       = [0.95 ]   # daily prob of survival
probSurv       = [0.93, 0.95 ]   # daily prob of survival
probMortFlood  = [0.1] # 10% of failed nests - not of all nests
SprobSurv      = [0.2] # daily survival prob during storms - like intensity?
SprobMortFlood = [1.0] # all failed nests during storms fail due to flooding
stormDur       = [1,2,3]
# NOTE maybe duration should be in hours so I can test more values?
# stormDur       = [1,3]
stormFreq      = [1,2,3]
# stormFreq      = [1,3]
obsFreq        = [3,5,7]
# obsFreq        = [3]
hatchTime      = [16,20,24,28]
# hatchTime      = [16, 28]
# whether or not nests that end during storms are marked "unknown" (vs. flooded)
#assignUnknown  = [0,1] 
assignUnknown  = [0] 
#probCorrect     = 0.8
#fateCuesPresent = [0.7,0.8,0.9] 
fateCuesPresent = [0.8] 
# 80% chance that field cues about nest fate are present/observed
# based on around 80% success rate in identifying nest fate (from camera study)
# this is essentially uncertainty?
numMC           = 0 # number of nests misclassified
paramsList      = list(
    product(
        #    0        1         2          3         4          5          
        numNests, probSurv, SprobSurv, stormDur, stormFreq, hatchTime, 
        #   6          7              8           9           10
        # obsFreq, probMortFlood, breedingDays, discProb, fateCuesPresent
        obsFreq, probMortFlood, breedingDays, discProb 
        )
        )
# NOTE fate cues prob will now be decided based on obs_int (11-14)
# NOTE the above probably only needs to be for the params that i'm testing?
#print(">> params list:", paramsList, "length:", len(paramsList))
paramsArray = np.array(paramsList) # don't want prob surv to be an integer!
print("NUMBER OF PARAM SETS:", len(paramsList))
nrows   = len(paramsList)*nreps*nruns
# -----------------------------------------------------------------------------
#   IMPORT REAL DATA
#   weekly storm probability and weekly nest initiation probability:
# -----------------------------------------------------------------------------
init= np.genfromtxt(
        #fname="/mnt/c/Users/Sarah/Dropbox/nest_models/storm_init3.csv",
        fname="C:/Users/Sarah/Dropbox/nest_models/storm_init3.csv",
        dtype=float,
        delimiter=",",
        skip_header=1,
        usecols=2
        )
initProb = init / np.sum(init)
stormProb = np.genfromtxt(
        #fname="/mnt/c/Users/Sarah/Dropbox/nest_models/storm_init3.csv",
        fname="C:/Users/Sarah/Dropbox/nest_models/storm_init3.csv",
        dtype=float,
        delimiter=",",
        skip_header=1,
        usecols=3 # 4th column 
        )
storm_weeks2 = np.arange(14,29,1)
weekStart = (storm_weeks2 * 7) - 90
weekStart = weekStart.astype(int)
if debug:
    print(">> storm prob, by week:",stormProb, len(stormProb),"sum=",np.sum(stormProb))
    print(">> initiation prob, by week:",initProb,len(initProb),"sum=",np.sum(initProb))
    print(">> start day for each week:", weekStart, len(weekStart))
# -----------------------------------------------------------------------------
#   SAVE FILES 
# -----------------------------------------------------------------------------
# name for unique directory to hold all output:
dirName    = datetime.today().strftime('%m%d%Y_%H%M%S') 
print(">> save directory name:", dirName)
todaysDate = datetime.today().strftime("%Y%m%d")
likeFile   = Path(uniquify(Path.home() / 
                           'C://Users/Sarah/Dropbox/nest_models/py_output' / 
                           dirName / 
                           'ml_values.csv'))
likeFile.parent.mkdir(parents=True, exist_ok=True)
print(">> likelihood file path:", likeFile)
column_names     = np.array([
    'rep_ID', 'mark_s', 'psurv_est', 'ppred_est', 'pflood_est', 'stormsurv_est', 
    'stormpred_est', 'stormflood_est', 'storm_dur', 'storm_freq', 'psurv_real', 
    'stormsurv_real','pflood_real', 'stormflood_real', 'hatch_time','num_nests', 
    'obs_int', 'num_discovered','num_excluded', 'exception'
    ])
colnames = ', '.join([str(x) for x in column_names]) # needs to be string
# -----------------------------------------------------------------------------
#   FUNCTIONS
# -----------------------------------------------------------------------------
# Some are very small and specific (e.g. logistic function); others are 
# quite involved.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def stormGen(frq, dur):
    out = rng.choice(a=weekStart, size=frq, replace=False, p=stormProb)
    dr = np.arange(0, dur, 1)
    stormDays = [out + x for x in dr]
    stormDays = np.array(stormDays).flatten()
    print(">> storm days:", stormDays)
    return(stormDays)

# -----------------------------------------------------------------------------
# This function remaps values from R^2 into the lower left triangle located 
# within the unit square.
def triangle(x0, y0):
    if y0 > x0:
        ret = triangle(y0, x0)
        return ret[1], ret[0]

    r0 = np.sqrt( x0**2 + y0**2)
    m  = 1.0
    if y0 != x0:
        m = y0/x0

    theta = np.arctan(m)
    r3    = r0 * 1.0/(1.0 + m)
    x3    = r3 * np.cos(theta)
    y3    = r3 * np.sin(theta)
    return x3, y3

# -----------------------------------------------------------------------------
# This is just the logistic function
# Trying out type hints (PEP 484) to keep output from overflowing
#def logistic(x)->np.float128:
def logistic(x)->np.longdouble:
    #return 1.0/( 1.0 + math.exp(-x) )
    return 1.0/( 1.0 + np.exp(-x) )

# -----------------------------------------------------------------------------
# This function computes intersection of 2 arrays more quickly than intersect1d
#    > possible observations = intersection of observable & survey days

#@profile
def in1d_sorted(A,B): 
    idx = np.searchsorted(B, A)
    idx[idx==len(B)] = 0
    return A[B[idx] == A]

# -----------------------------------------------------------------------------
# This function creates the list of survey days by taking a random start date 
# from the first 5 breeding days and creating a range with step size determined
# by observation frequency. Then remove storm days.
def mk_surveys(stormDays, obsFreq, breedingDays):
    # first day of each week because the initiation probability is weekly 
    # the upper value should not be == to the total number of season days 
    # because then nests end after season is over 

    start       = rng.integers(1, high=5) # random day of 1st survey from 1st 5 breeding days          
    end         = start + breedingDays
    surveyDays  = np.arange(start, end, step=obsFreq)
    stormSurvey = np.isin(surveyDays, stormDays) 
    surveyDays  = surveyDays[np.isin(surveyDays, stormDays) == False] # keep only values that aren't in storm_days 
    print(">> all survey days, minus storms:\n", surveyDays, len(surveyDays)) 

    return(surveyDays)

# -----------------------------------------------------------------------------
def survey_int(surveyDays):
    surveyInts = np.array( [surveyDays[n] - surveyDays[n-1] for n in range(1, len(surveyDays)-1) ] )
    #print(">> interval between current survey and previous survey:\n", surveyInts, len(surveyInts))

    return(surveyInts)

# -----------------------------------------------------------------------------
#   NEST DATA COLUMNS AND PARAM LIST: 
# -----------------------------------------------------------------------------
# 0) ID number              | 4) flooded? (T/F)        | 8) j (last active)  
# 1) initiation date        | 5) surveys til discovery | 9) k (last checked) 
# 2) end date               | 6) discovered? (T/F)     | 10) assigned fate
# 3) hatched? (T/F)         | 7) i (first found)       | 11) num obs intervals 
# 12) final int storm (T/F) | 13) num other storms 

#paramsList      = list(product(
#          0           1         2          3         4          5
#        numNests, probSurv, SprobSurv, stormDur, stormFreq, hatchTime, 
#           6           7             8            9          10
#        obsFreq, probMortFlood, breedingDays, discProb, assignUnknown
#        ))
# -----------------------------------------------------------------------------
def randArgs():
    # Choose random initial values for the optimizer
    # These will be log-transformed before going through the likelihood function
    s     = rng.uniform(-10.0, 10.0)       
    mp    = rng.uniform(-10.0, 10.0)
    ss    = rng.uniform(-10.0, 10.0)
    mps   = rng.uniform(-10.0, 10.0)
    srand = rng.uniform(-10.0, 10.0) # should the MARK and matrix MLE start @ same value?
    z = np.array([s, mp, ss, mps, srand])

    return(z)

# -----------------------------------------------------------------------------
def ansTransform(ans):
        
    s0   = ans.x[0]         # Series of transformations of optimizer output.
    mp0  = ans.x[1]         # These make sure the output is between 0 and 1, 
    ss0  = ans.x[2]         # and that the three fate probabilities sum to 1.
    mps0 = ans.x[3]

    s1   = logistic(s0)
    mp1  = logistic(mp0)
    ss1  = logistic(ss0)
    mps1 = logistic(mps0)

    ret2 = triangle(s1, mp1)
    s2   = ret2[0]
    mp2  = ret2[1]
    mf2  = 1.0 - s2 - mp2

    ret3 = triangle(ss1, mps1)
    ss2  = ret3[0]
    mps2 = ret3[1]
    mfs2 = 1.0 - ss2 - mps2
    
    #ansTransformed = np.array([s2, mp2, mf2, ss2, mps2, mfs2], dtype=np.float128)
    ansTransformed = np.array([s2, mp2, mf2, ss2, mps2, mfs2], dtype=np.longdouble)
    #print(">> results as an array:\n", ansTransformed)
    #print(">> results (s, mp, mf, ss, mps, mfs, ex):\n", s2, mp2, mf2, ss2, mps2, mfs2, ex)
    return(ansTransformed)

# -----------------------------------------------------------------------------
def mk_nests(params, init, stormDays, surveyDays): 

    # 1. Unpack necessary parameters
    # NOTE about the params at the beginning of the script:
    # some have only 1 member, but they are still treated as arrays, not scalars

    hatchTime = int(params[5]) 
    stormDur  = int(params[3]) 
    stormFrq  = int(params[4]) 
    obsFreq   = int(params[6]) 
    numNests  = int(params[0]) 
    discProb  = params[9]
    pSurv     = params[1]       # daily survival probability
    pfMort    = params[2]       # prob of surviving (not flooding) during storm
    nCol      = 16              # number of output columns
    breedingDays  = params[8]   # number of days in breeding season
    fateCuesPresent = 0.6 if obsFreq > 5 else 0.66 if obsFreq == 5 else 0.75
    if debug: print(
        "observation frequency:", obsFreq, 
        "prob of correct fate:", fateCuesPresent
        )

    #@#print(">> discovery probability:", discProb)
    #@#print(">> hatch time:", hatchTime, "| storm duration:", stormDur, "| storm frequency:", stormFrq)
    #@#print(">> number of nests:", numNests, "| frequency of observations:", obsFreq)
    #@#print(">> probability of survival:", pSurv, "and failure due to flooding:", pfMort)
    #@#print(">> total days in the season:", breedingDays)

    # 2. Create lists of storm days and survey days, plus the intervals between surveys

    surveyInts = survey_int(surveyDays)
    nestData = np.zeros(shape=(numNests, nCol), dtype=int) 
    nestData[:,0] = np.arange(1,numNests+1) # nest ID numbers 

    # 3. Create a list of nest initiation dates (one for each nest) and make sure you have the correct number

    initWeek = rng.choice(a=weekStart, size=numNests, p=init)  # random starting weeks; len(a) must equal len(p)
    initiation = initWeek + rng.integers(7)                    # add a random number from 1 to 6 (?) 
    nestData[:,1] = initiation                                 # record to a column of the data array
    #print(">> initiation weeks:\n", initWeek, "\n>> initiation dates:\n", initiation, len(initiation)) 

    # 4. Decide how long each nest is active
    # >> use a negative binomial distribution - distribution of number of 
    #    failures until success 
    #     >> in this case, "success" is actually the nest failing, so use 
    #        1-pSurv (the failure probability) 
    #     >> gives you number of days until the nest fails (survival)
    #     >> if survival > incubation time, then the nest hatches 
    # >> then use survival to calculate end dates for each nest 
    #    (end = initiation + survival)

    survival = rng.negative_binomial(n=1, p=(1-pSurv), size=numNests) 
    survival = survival - 1 # but since the last trial is the success, need to subtract 1
    if debug: print(">> survival in days:\n", survival, len(survival)) 
    
    ## >> set values > incubation time to = incubation time (nest hatched): 
    ##      (need to because you are summing the survival time)
    survival[survival > hatchTime] = hatchTime
    nestEnd = initiation + survival # add num days survived to init date to get end date   
    if debug: print(">> end dates:\n", nestEnd, len(nestEnd)) 
    nestData[:,2] = nestEnd 
    hatched = survival >= hatchTime # the hatched nests survived for >= hatchTime days 
    ## NOTE THIS IS NOT THE TRUE HATCHED NUMBER; DOESN'T TAKE STORMS INTO ACCOUNT
    if debug: print("real hatch proportion:", hatched.sum()/numNests)
    nestData[:,15] = survival # number of days nest survived
    # NOTE Remember that int() only works for single values 

    # ---- FAILED NESTS --------------------------------------------------------

    # 5. Decide cause of failure for failed nests:
    # >> create a vector of probabilities, one for each failed nest, to decide
    #    whether nest flooded
    # >> probability value is compared to the prob of failure due to flooding
    #          on a regular day: prob of mortality = 1-DSR 
    #                            conditional prob of flooding = 0.05
    #          on a storm day: prob of mortality = 0.9 
    #                          conditional prob of flooding = 1

    pflood = rng.uniform(low=0, high=1, size=numNests) 
    # need to check whether this is the correct distribution 
    # NOTE: still needs to be conditional on nest having failed already...  
    # NOTE np.concatenate joins existing axes, while np.stack creates new ones
    nestPeriod = np.stack((initiation, nestEnd))
    nestPeriod = np.transpose(nestPeriod) # an array of start,end pairs 
    if debug: print(
        ">> start and end of nesting period:\n", 
        nestPeriod, 
        nestPeriod.shape
        )
    
# -----------------------------------------------------------------------------
    # >> stormNestIndex searches for storm days w/in active period of each nest
    #     >> returns index where storm day would be within the active interval: 
    #             0 = before init; 2 = after end; 1 = within interval
    #     >> fate cues should become harder to interpret after storms
    stormNestIndex = searchSorted2(nestPeriod, stormDays)
    if debug: print("where were storms in active period?", stormNestIndex)
    # if index == 1, then storm Day is within the period interval: 
    stormNest = np.any(stormNestIndex == 1, axis=1) 
    if debug: print(
        "which nests were active during a storm? How many?", 
        stormNest, 
        stormNest.sum()
        )
    numStorms = np.sum(stormNestIndex==1, axis=1) # axis=1 means summing over rows?
    # NOTE I *think* this is actually number of storm intervals, which is what 
    # we want for the likelihood function. so that would be good...
    if debug: print(">> number of storm intervals during nesting period:\n", numStorms)
    stormNest = stormNest.astype(int) # of those nests, tiny fraction do survive
# -----------------------------------------------------------------------------
    flooded = np.where(pflood>pfMort, 1, 0) # if pflood>pfMort, flooded=1, else flooded=0 
    # print("flooded", flooded)
    # since it's 1 and 0, can use arithmetic: 
    # print("flooded & storm:", stormNest + flooded)
    floodedAll = stormNest + flooded > 1 # both need to be true 
    #@#print("number flooded:", flooded.sum())
    # print(">> number flooded and during storm:\n", floodedAll.sum()) 
    # print(">> flooded & during storm:", floodedAll)
    #print(">> number of nests flooded:\n", sum(floodedAll==True), "of", len(floodedAll)) 
    #@#print(">> hatched as integers:\n", hatched.astype(int)) 
    #print(">> inverse of hatched:\n", (~hatched).astype(int)) 
    #floodedAndFailed = (floodedAll + (~hatched).astype(int)) > 1 
    # flooded == true and hatched == False; see NOTE above 
    # where flooded == true (1) and ~hatched == true(1), result will be > 1; else marked false
    # where you put the parentheses in ~hatched.astype(int) is important!  
    #nestData[:,4] = floodedAndFailed.astype(int) 
    nestData[:,4] = floodedAll.astype(int) 
    #print(">> did nest fail due to flooding:\n", floodedAndFailed, len(floodedAndFailed)) 
    #print(">> number of nests failed due to flooding:\n", sum(floodedAndFailed==True), len(floodedAndFailed)) 
    trueFate = np.empty(numNests) 
    #print("length where flooded = True", len(trueFate[floodedAll==True]))
    trueFate.fill(1) # nests that didn't flood or hatch were depredated 
    #trueFate[flooded==True] = 2 
    #trueFate[floodedAll==True] = 2 
    trueFate[hatched    == True] = 0 # was nest discovered?  
    trueFate[floodedAll == True] = 2 
    # trueFate[floodedAll ] = 2 
    # hatched fate is assigned last, so hatched is taking precedence over flooded
    # so maybe that's why the true DSR is always higher than 0.93, but what about true DSR of discovered?
    fates = [np.sum(trueFate==x) for x in range(3)]
    # print(">>>>> true final nest fates:\n", trueFate)

    # ---- TRUE DSR ------------------------------------------------------------

    # Calculate proportion of nests hatched and use to calculate true DSR
    #   daily mortality = num failed / total exposure days
    #     (num failed =  total-num hatched) 
    #     (total exposure days = add together survival periods)
    #   DSR = 1 - daily mortality
    trueHatch = trueFate==0 # true/false did nest hatch (after storms accounted for)?
    nestData[:,3] = trueHatch.astype(int)
    trueDSR2 = 1 - ( (numNests - trueHatch.sum()) / survival.sum() ) 
    #@#print(">>>> total exposure days:", survival.sum())
    # print(">>>>> and true DSR, calculated correctly:", trueDSR2)

    # ---- NEST DISCOVERY & OBSERVATION ----------------------------------------

    # NOTE this name is a little misleading - it's actually survey days til discovery
    daysTilDiscovery = rng.negative_binomial(n=1, p=discProb, size=numNests) # see above for explanation of p 
    nestData[:,5]    = daysTilDiscovery 
    # print(">> survey days until discovery:\n", daysTilDiscovery, len(daysTilDiscovery)) 

    # what are the first and last *possible* survey dates for the nest?  
    # this finds the index of where initiation would be inserted in surveyDays 
    position = np.searchsorted(surveyDays, initiation) 
    #@#print(">> position of initiation date in survey day list:\n", position, len(position)) 
    firstSurvey = surveyDays[position] 
    # NOTE last possible survey will be after nest ends (either hatch or fail) 
    position2 = np.searchsorted(surveyDays, nestEnd)
    #@#print(">> position of end date in survey day list:\n", position2, len(position2)) 
    lastSurvey = surveyDays[position2] 

    possSurveyDates = np.stack((firstSurvey, lastSurvey)) # double parens so numpy knows it's not two separate arguments
    possSurveyDates = np.transpose(possSurveyDates)
    # print(
            # ">> start of nest, end of nest, first possible survey, last possible survey:\n",
            # np.concatenate((nestPeriod,possSurveyDates), axis=1)
            # )
    totalSurvey = position2-position + 1# this gives the difference in index values from surveyDays, AKA number of surveys
    # number of surveys when nest is active (when it could be discovered)
    #print(">> total possible surveys to observe nest (position2 - position):\n", position2-position) 
    # print(">> total possible surveys to observe nest (position2 - position):\n", totalSurvey) 
    #totalReal  = totalSurvey - daysTilDiscovery + 1 # totalSurvey - (daysTilDiscovery - 1)
    totalReal  = totalSurvey - daysTilDiscovery  # totalSurvey - (daysTilDiscovery - 1)
    # print(">> actual number of surveys/observations for each nest:\n", totalReal, len(totalReal))
    #nestData[:,11] = totalReal-1 # number of observation intervals is number observations - 1 
    #nestData[:,11] = totalReal - 1 # number of observation intervals
    nestData[:,11] = position2-position # number of observation intervals
    # print("number of obs intervals:", nestData[:,11])
    # need to decide if it's number of observations or observation intervals

    discovered = daysTilDiscovery < (position2-position) 
    # print(">> nest discovered if surveys til discovery < total possible surveys (surveys while active):\n", discovered)
    # print(">> proportion of nests discovered:", sum(discovered)/numNests, "vs. expected proportion:", discProb)
    nestData[:,6] = discovered.astype(int) # convert to numeric for numpy array  

    firstFound = np.zeros(numNests) 
    firstFound[discovered==True] = surveyDays[position+daysTilDiscovery][discovered==True] 
    #print(">> nest first found:\n", firstFound, len(firstFound)) 

    lastActive = np.zeros(numNests)
    lastActive[discovered==True] = surveyDays[position2][discovered==True] 
    #print( ">> nest last active:\n", lastActive, len(lastActive)) 
    
    # Last checked will be one survey after last active, unless hatch == True
    lastChecked = np.zeros(numNests) 
    lastChecked[discovered==True] = surveyDays[position2+1][discovered==True] 
    #@#print(">> nest last checked, w/o hatch:\n", lastChecked, len(lastChecked)) 
    #lastChecked[hatched==True] = lastActive[hatched==True] # take hatch date into account
    lastChecked[trueHatch==True] = lastActive[trueHatch==True] # take hatch date into account
    #print(">> nest last checked:\n", lastChecked, len(lastChecked)) 

    nestData[:,7]  = firstFound 
    nestData[:,8]  = lastActive 
    nestData[:,9]  = lastChecked 

    # ---- STORMS DURING OBSERVATION PERIOD -------------------------------------------------------------------------

    obsPeriod      = lastChecked - firstFound 
    # print(">> length of observation period:\n", obsPeriod, len(obsPeriod)) # is the actual number of days in the period greater than expectedsPeriod) 
    # numStorms = (obsPeriod/obsFreq - totalReal) / obsFreq # number extra days in period divided by normal interval length 
    # NOTE already calculated numStorms earlier in function
    # print(">> number of storm intervals during nest active:", numStorms )
    # I changed how surveyInts was defined; now it's current minus previous, not next minus current
    # stormInterval = surveyInts > obsFreq
    #@#print(">> was there a storm in this interval?\n", stormInterval)
    stormIntFinal = surveyInts[position2] > obsFreq  # was obs interval longer than usual? (== there was a storm)
    #@#print(">> was there a storm in the final interval for a given nest?\n", stormIntFinal, len(stormIntFinal)) # basically ask is final interval > 3 
    # needs to be positionEnd -1 because of how intervals are calculated 
    nestData[:,12] = stormIntFinal.astype(int) #nestData[:,13] = numStorms - as.integer(stormFinal) 
    #@#print(">> stormIntFinal in integer form:\n", nestData[:,12], stormIntFinal.shape)
    nestData[:,13] = numStorms 

    assignedFate = np.zeros(numNests) # if there was no storm in the final interval, correct fate is assigned 
    assignedFate.fill(7) # default is unknown; fill with known fates if field cues allow
    # print(">> assigned fate array & its shape:\n", assignedFate, assignedFate.shape)

    #@#print("days from end date to final check:", lastChecked - nestEnd)
    #y = 3
    #assignedFate[lastChecked-nestEnd < y] = trueFate[lastChecked-nestEnd < y]

    fateCuesProb = rng.uniform(low=0, high=1, size=numNests)
    # print(">> random probs to compare to fateCuesPresent:\n", fateCuesProb, fateCuesProb.shape)
    assignedFate[fateCuesProb < fateCuesPresent] = trueFate[fateCuesProb < fateCuesPresent]
    # print("assigned fates:", assignedFate)
    # fate cues prob should be affecting all nest fates equally, not just failures.
    # print(">> proportion of nests assigned hatch fate:", np.sum((assignedFate==0)[discovered==True])/(sum(discovered==True)),"vs period survival:", pSurv**hatchTime)

    aFates = [np.sum((assignedFate == x)[discovered==True]) for x in range(4)]
    # this proportion needs to be out of nests discovered AND assigned
    aFatesProp = [np.sum((assignedFate == x)[discovered==True])/(np.sum(discovered==True)) for x in range(4)]
    tFates = [np.sum((trueFate == x)[discovered==True]) for x in range(4)]
    tFatesProp = [np.sum((trueFate == x)[discovered==True])/(np.sum(discovered==True)) for x in range(4)]
    # print(
    #         ">> assigned fate (hatched, depredated, flooded, unknown):", 
    #         # ">> assigned fate proportions (hatched, depredated, flooded, unknown):", 
    #         aFatesProp[0:3], 
    #         (np.sum(discovered==True) - np.sum(aFates)) / (np.sum(discovered==True)),
    #         # (np.sum(aFates==7) )/ (np.sum(discovered==True)),
    #         "\n\n>> proportions of known (assigned) fates (H, D, F):",
    #         aFates[0:3]/np.sum(aFates),
    #         "\n>> vs. actual proportions for discovered only (H, D, F):",
    #         tFatesProp[0:3],
    #         # np.sum(discovered==True)- np.sum(tFates)
    #         )
    nestData[:,10] = assignedFate
    expDays = survival.sum()
    # print(">> exposure days:", expDays)
    exp = np.zeros(numNests)
    exp.fill(expDays)
    nestData[:,14] = exp
 
    numDisc = discovered.sum()
    # print(">> nest discovered?", discovered)
    numDiscH = trueHatch[discovered==True].sum()
    # print(">> discovered & hatched:", discovered[trueHatch==True])
    survDisc = survival[discovered==True].sum()
    # print(">> calculate exposure days for discovered nests by summing this list:", survival[discovered==True])
    # print(">> number discovered:", numDisc, ", number discovered that hatched:", numDiscH, ", and exposure days (discovered nests):", survDisc)
    #trueDSR_disc = 1 - ((discovered.sum() - hatched.sum()) / survival.sum()) 
    trueDSR_disc = 1 - ((numDisc - numDiscH) / survDisc) # num failed / total exposure days
    # print(">> true DSR of discovered nests only:", trueDSR_disc) 
    # NOTE issue may be that not enough flooded nests are discovered, not that too many hatched nests are

    nests = np.stack((discovered, trueFate, firstFound, lastActive, lastChecked, totalReal, numStorms))
    nests = np.transpose(nests)
    # print(">> discovered?, fate, i, j, k, num surveys, num storm intervals:\n", nests)
    # print(
    #     ">> discovered?, fate, i, j, k, num surveys, num storm intervals:\n", 
    #     nests[:5,:], 
    #     "\n ... \n",
    #     nests[-5:,:]
    #     )

    #@#print(">> nest data:\n----id--ini-end-hch-fld-std-dsc-i--j--k-fate-nobs-sfin-nstm\n", nestData)

    # >>>>>>>>>>> SAVE NEST DATA AS CSV (may be more memory-intensive) >>>>>>>>>>>>>>
    #now       = datetime.now().strftime("%H%M%S")
    #nests  = Path.home()/'/mnt/c/Users/sarah/Dropbox/nest_models/py_output'/dirName/('nests'+now+'.csv')
    #nests.parent.mkdir(parents=True, exist_ok=True)

    #csvHead = "ID, initiation, end, hatched, flooded, surveys_til_discovery, discovered,\
    #           i, j, k, a_fate, num_obs, storm_during_final, num_storms" 

               #num_stm_survey_excl_final 
    #np.savetxt("nestdata.csv", nestData, fmt="%d", delimiter=",", header=csvHead) 
    #np.savetxt(nests, nestData, fmt="%d", delimiter=",", header=csvHead) 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    return(nestData) 

# -----------------------------------------------------------------------------
#   MAYFIELD
# -----------------------------------------------------------------------------
# Mayfield's original estimator was defined as: 
#                DSR = 1 - (# failed nests / # exposure days)
# so if DSR = 1 - daily mortality, then:
# daily mortality = # failed nests / # exposure days

# NOTE: Johnson (1979) provided a mathematical derivation that allowed the 
#       calculation of variance for the estimate.
# > for a single day:
#    > probability of survival is s    
#    > probability of failure is (1-s)
# > for interval of length k days:
#    > prob of survival is s**k 
#    > prob of failure is s**(1/2k-1)(1-s)
# > ex. - prob of a nest surviving three days and failing on the fourth is:
#         s*s*s*(1-s) 
#     > this assumes that a failed nest survived half (minus a day)
#       of interval and then failed
#
# Johnson's modified ML estimator: 
#        mortality = (f1 + sum(ft)) / (h1 + sum(t*ht) + f1 + 0.5 sum(t*ft)) 
# > created by differentiating the log-likelihood equation and setting to 
#   zero (maximizing)
# > f1 and h1 represent an interval between visits of one day, which is not 
#   used in our studies
#   > so we end up with: sum(ft) / (sum(t*ht) + 0.5*sum(t*ft)) 
#          where t = interval length, and 
#          f and h represent number of failures and hatches, respectively
# -----------------------------------------------------------------------------
def mayfield(ndata):
#    I am assuming the nest data that is input has already been filtered to only discovered nests w/ known fate
#    dat = ndata[
    hatched = np.sum(ndata[:,3])
    failed = len(ndata) - hatched 
    # print(">> Calculate Mayfield estimator.")
    # print(">> number of nests hatched:", hatched, "and failed", failed)
    mayf = failed / (hatched + 0.5*failed)
    # print(">> Mayfield estimator of daily mortality (1-DSR) =", mayf) 

    return(mayf)

# -----------------------------------------------------------------------------
#   PROGRAM MARK 
# -----------------------------------------------------------------------------
# The model used in Program MARK is based on Dinsmore (2002) -  
#        allows for variance in DSR & use of covariates

# This function is based on info in 'Program MARK: A Gentle Introduction' 

# It also has a wrapper function that transforms the initial optimizer values
# using the logistic function.
# This way, the optimizer can work over the range of -infinity:infinity, but
# the values fed to the function are between 0 and 1 (probabilities)

def prog_mark(s, ndata):

    # 1. grab the data for the input for MARK i
    # -- first, grab only discovered nests, then grab only needed rows:
    # print("all nests:", len(ndata)) 
    # # function is called only on discovered nests
    disc = ndata[np.where(ndata[:,6]==1)]
    # print("all nests:", len(ndata), "discovered nests only:", len(disc))
    # (nest ID, first found, last active, last checked, assigned fate)
    # inp[0] = ID | inp[1] = i | inp[2] = j | inp[3] = k | inp[4] = fate `
    # inp = ndata[:,np.r_[0,7:11]] # doesn't include index 11
    inp = disc[:,np.r_[0,7:11]] # doesn't include index 11
    #@#print("inp (ID, i, j, k, fate:)\n",inp)
    nocc = ndata[0,10] # load number of occasions from data
    l    = len(inp)
    #print("l=", l, "| s=", s, "| nocc=", nocc)
    #print("----------------------------")

    # 2. extract rows where j minus i does not equal zero (nest wasn't only 
    #    observed as active for one day)
    # > The model requires all nests to have at least two observations while 
    #   active
    inp = inp[(inp[:,2]-inp[:,1]) != 0] # access all rows; 
    #                                     create mask, index inp using it
    #@#print("----------------------------")
    #@#print(inp)
    #@#print("----------------------------")
    #@#print("length=", len(inp), "| s=", s)

    # 3. create vectors to store:
    # > a) the probability values for each nest
    # > b) the degrees of freedom for each value
    allp   = np.array(range(1,len(inp)), 
                      dtype=np.longdouble) # all nest probabilities 
    alldof = np.array(range(1,len(inp)), 
                      dtype=np.double) # all degrees of freedom
    #print("length of array:", len(allp))

    #4. fill the vectors
    for n in range(len(inp)-1): # want n to be the row NUMBER
        #print("nest ID=", n, "| i=", inp[n,1], "| j=", inp[n,2], "| k=", inp[n,3])

        # for the basic case where psurv is constant across all nests and times:
        #   - count the total number of alive days
        #   - count the number of days in the final interval (for failed nests)
        alive_days = inp[n,2] - inp[n,1] # interval from first found to last active = all nests are KNOWN to be alive
        alive_days = alive_days - 1 # since this is essentially 1-day intervals, need 1 fewer than total number

        final_int  = inp[n,3] - inp[n,2] 
        final_int  = final_int - 1
        # for failed nests, failure occurs in interval from last active to last checked 
        # for hatched nests, last active == last checked (assumption is that hatched nests are found on hatch day)
        #print("days nest was alive:",alive_days,"& final int:",final_int)
        # NOTE need nests to be alive for at least one interval

        # Probability equation: 
        # > daily probability of survival (DSR) raised to the power of intervals nest was known to be alive
        # > for hatched nests, that's it
        # > for failed nests, exact failure date can be unknown
        # > > but we know nest wasn't alive for the entire final interval
        # > > so add in probability of NOT surviving one interval (1-DSR)
        # > > > EX: if probability of surviving from day 1-3 is s1*s2*s3, then
        # > > >     probability of failure sometime during days 4-6 is 1-s4*s5*s6
        # > hatched nests also have one extra degree of freedom (dof)
        
        # FAILED NESTS will have final_int (last checked-last active) greater than zero
        if final_int > 0:
            p   = (s**alive_days)*(1-(s**final_int)) 
            dof = alive_days
            #print(">> nest", inp[n,0], "failed. likelihood=", p)
            # but if final_int==0 (survived), this will become zero
        # HATCHED NESTS will have final_int == 0
        else:
            p   = s**alive_days
            dof = alive_days + 1
            #print(">> nest", inp[n,0], "hatched. likelihood=", p)
            #print("probability=", p)

        allp[n]   = p # NOTE this line is throwing the Deprecation Warning
        #allp[n,0]   = p # NOTE this line is throwing the Deprecation Warning
        # apparently warning means that n is a 2d array, so need to index it
        # never mind, apparently it's one dimensional 
        alldof[n] = dof
        #    once we have all the probabilities then...?
        #    where do exposure days come into it?
        #    might just be alive_days + final_int
        #    technically, raise prob to power of frequency
        #    we are assuming each occurs only once
    #print(">> all nest cell probabilities:\n", allp)
    #print(">> all degrees of freedom:\n", alldof)
    #lnp  = np.log(allp) # vector of all log-transformed probabilities
    lnp  = -np.log(allp) # vector of all log-transformed probabilities
    #print("log of all nest cell probabilities:", lnp)
    #print(">> negative log likelihood of each nest cell probability:", lnp)
    #lnSum = lnp.sum()
    #NLL = -1*lnp.sum()
    NLL = lnp.sum()
    #print(">> sum to get negative log likelihood of the data:", NLL)
    return(NLL)

# -----------------------------------------------------------------------------
def mark_wrapper(srn, ndata):
    # This function calls the program MARK function when given a random starting value (srn) and some nest data (ndata)
    # > values given to the optimizer are transformed before being given to the MARK function
    # > > this allows a larger range of values for the optimizer to work over without overflow
    # > > but the values given to the function are still between 0 and 1, as is required
    # > Create vector to store the log-transformed values, then fill:
    s = np.ones(numNests, dtype=np.longdouble)
    s = logistic(srn)
    # NOTE is this multiple random starting values (for each nest) or one random starting value?
    #@#print("logistic of random starting value for program MARK:", s, s.dtype)
    # the logistic function tends to overflow if it's a normal float; make it np.float128
    ret = prog_mark(s, ndata)
    #@#print("ret=", ret)
    return ret

# ---------------------------------------------------------------------------------------
 # this function does the matrix multiplication for a SINGLE interval of length intElt days 
 # during observaton, nest state is assessed on each visit to form an observation history 
 # the function calculates the negative log likelihood of one interval from the observation history 
 # these can then be multiplied together to get the overall likelihood of the observation history 
 # the function takes 4 arguments: 
 #    > intElt - length in days of the observation interval being assessed 
 #    > initial state (stateI) - state of the nest at the beginning of this interval 
 # for analysis, these nest states are coded as a 1-dimensional matrix (vector): 
 # 
 #               [ 1 0 0 ] (alive)       [ 0 1 0 ] (failed - predation)        [ 0 0 1 ] (failed - flooding) 
 #
 # the formula used is from Etterson et al. (2007) 
 #    > in this case, the nest started the interval alive and ended it alive as well 
 #    > daily nest probabilities: s - survival; mp - mortality from predation; mf - mortality from flooding 
 #    > since these are daily probabilities, we raise the base matrix to the power of the number of days in the interval # 
 #                                          _         _  intElt           _   _ 
 #                 [ 1 0 0 ]               |  s  0  0  |                 |  1  | 
 #                                 *       |  mp 1  0  |            *    |  0  | 
 #                                         |_ mf 0  1 _|                 |_ 0 _|  
 #                              
 #           {  transpose(stateI)  *  useM, raised to intElt power  *    stateF  } 

# -----------------------------------------------------------------------------
#   THE LIKELIHOOD FUNCTION
# -----------------------------------------------------------------------------

# This function computes the overall likelyhood of the data given the model parameter estimates.
#
# The model parameters are expected to be received in the following order:
# - a_s   = probability of survival during non-storm days
# - a_mp  = conditional probability of predation given failure during non-storm days
# - a_mf  = conditional probability of flooding given failure during non-storm days
# - a_ss  = probability of survival during storm days
# - a_mfs = conditional probability of predation given failure during storm days
# - a_mps = conxditional probability of predation given failure during storm days
# - sM    = for program MARK?

#def like_old(a_s, a_mp, a_mf, a_ss, a_mfs, a_mps, nestData, stormDays, surveyDays, obs_int):
def like_old(argL, obsFreq, nestData, surveyDays, stormDays ):
  
    nCol = 18
    a_s, a_mp, a_mf, a_ss, a_mps, a_mfs, sM = argL
    obs_int = obsFreq
    likeData = np.zeros(shape=(numNests, nCol), dtype=np.longdouble) 

    stillAlive  = np.array([1, 0, 0])
    mortFlood   = np.array([0, 1, 0])
    mortPred    = np.array([0, 0, 1])

    # > starting matrix, from etterson 2007:
    startMatrix = np.array([[a_s,0,0], [a_mf,1,0], [a_mp,0,1]]) 
    # > use this matrix for storm weeks:
    stormMatrix = np.array([[a_ss,0,0], [a_mfs,1,0], [a_mps,0,1]]) 
        # how is this matrix actually being incorporated during the analysis?

    logLike = Decimal(0.0)          # initialize the overall likelihood counter
    #logLike = float(logLike)
    rowNum = 0
    for row in nestData:
    # columns for the likelihood comparison datarframe:
    #         num obs, log lik of nest, log lik period 1, period 2, etc.  
        
    # FOR EACH NEST -------------------------------------------------------------------------------------

        nest    = row         # choose one nest (multiple output components from mk_obs())
        # print('obs_int check: ', obs_int)

        disc    = nest[7].astype(int)    # first found
        endObs  = nest[9].astype(int)    # last observed
        hatched = nest[3]
        flooded = nest[4]

        if flooded == True & hatched == False:
            fate = 2
        elif flooded == False & hatched == False:
            fate = 3
        else:
            fate = 1
        # fate    = nest[7].astype(int)    # assigned fate

        if np.isnan(disc):
            ###print("this nest was not discovered but made it through")
            continue
        
        num     = len(np.arange(disc, endObs, obs_int)) + 1 # number o observations
        # print('#############################################################')
        # print('nest =', nest[0], 'row number =', rowNum, 'number of obs =', num)

        likeData[rowNum, 0:3]  = np.array([nest[0], fate, num] )

        # print("num=", num)
        obsDays = in1d_sorted(
            # (np.linspace(disc, endObs, num=num)), surveyDays)
            np.linspace(disc, endObs, num=num),
            surveyDays)
        # print("obs days for nest:", obsDays)
        obsPairs = np.fromiter(
            itertools.pairwise(obsDays), 
            dtype=np.dtype((int,2))
            ) # do elements of numpy arrays have to be floats?
        # print("date pairs in observation period:", obsPairs) 

        # > make a list of intervals between each pair of observations 
        #   (necessary for likelihood function)
        intList = obsPairs[:,1] - obsPairs[:,0]
        # print("interval list:", intList)
        
        # > start off with all intervals = alive
        obs     = [stillAlive for _ in range(len(obsPairs)+1)] 

        # > change the last obs if nest failed:
        if fate == 2:
            obs[-1] = mortFlood
        elif fate == 3:
            obs[-1] = mortPred

        # print("fate, obs = ", fate, " , ", obs) # check that last entry in obs corresponds to fate

        # if hatch, leave as is?

        # place this likelihood counter inside the for loop so it resets 
        # with each nest:
        logLikelihood = Decimal(0.0)   
        #logLikelihood = float(logLikelihood)

        # likeData[0, rowNum] = nest[0]
        obsNum = 0
        for i in range(len(obs)-1):
        # FOR EACH OBSERVATION OF THIS NEST ---------------------------------------------------------------------

            # print("observation number:", obsNum)
            
            intElt  = (intList[i-1]).astype(int)  # access the (i-1)th element of intList,
                                    # which is the interval from the (i-1)th
                                    # to the ith observation

            #stateF  = obs[i]
            stateF  = obs[i+1] 
            stateI  = obs[i]
            # print("stateF:",stateF)
            TstateI = np.transpose(stateI)
            # print("TstateI:", TstateI)

            # if any(d in storm_days for d in range(i-1, i)):
            if any(d in stormDays for d in range(i-1, i)):
                # if any of the days in the current observation interval (range) is also in storm days, use storm matrix
                # print("using storm matrix")
                lDay = np.dot(stateF, np.linalg.matrix_power(stormMatrix, intElt))
                # this is the dot product of the current state of the nest and the storm matrix ^ interval length
           # look into using @ instead of nest dot calls 
            else:
                # print("using normal matrix")
                lDay = np.dot(stateF, np.linalg.matrix_power(startMatrix, intElt))

            lPer = np.dot(lDay, TstateI)
            # print("likelihood for this interval:", lPer)

            logL = Decimal(- np.log(lPer))
            # print("negative log likelihood of this interval:", logL)

            #logL = float(logL)

            logLikelihood = logLikelihood + logL # add in the likelihood for this one observation
            # print("log likelihood of nest observation history:", logLikelihood)
            colNum = obsNum + 4 
            likeData[rowNum, colNum] = logL
            obsNum = obsNum + 1

        likeData[rowNum,3] = logLikelihood
        logLike = logLike + logLikelihood        # add in the likelihood for the observation history of this nest
        rowNum  = rowNum + 1
        # print("increment row number:", rowNum)
        # print("overall log likelihood so far:", logLike)
        

    # print(
    #     "nest num, fate, num obs, lik(obs hist), lik(each obs int) ... \n",
    #     likeData[:10,:] # print first ten rows
    # )
    return(logLike) 

# -----------------------------------------------------------------------------
# try to keep these in numpy:
def like(perfectInfo, hatchTime, argL, numNests, obsFreq, nestData, surveyDays):
    # perfectInfo == 0 or 1 to tell you whether you know all nest fates or not
# NEST DATA COLUMNS: 
# 0) ID number              | 4) flooded? (T/F)        | 8) j (last active)  
# 1) initiation date        | 5) surveys til discovery | 9) k (last checked) 
# 2) end date               | 6) discovered? (T/F)     | 10) assigned fate
# 3) hatched? (T/F)         | 7) i (first found)       | 11) num obs intervals 
# 12) final int storm (T/F) | 13) num other storms 
    # ---------------------------------------------------------------------------------------------------
    # 1. Unpack:
    #    a. Initial values for optimizer:
    a_s, a_mp, a_mf, a_ss, a_mps, a_mfs, sM = argL

    # b. Observation history values from nest data:
    #@#print(">> nest data going into likelihood function:\n", nestData, nestData.shape)
    flooded = nestData[:,4]
    hatched = nestData[:,3] # actual hatched nests (after storms accounted for)
    ff      = nestData[:,7]
    la      = nestData[:,8]
    lc      = nestData[:,9]
    numInt  = nestData[:,11]
    # print("number of obs intervals:", numInt)

    # obsDays  = int(la - ff) # doesn't work bc it's an array
    obsDays  = la - ff 
    finalInt = lc - la
    
    #ff, la, lc = nestData[:,8:10] # doesn't work bc they aren't lists
    # print(
    #     "i, j, k:", ff, la, lc, 
    #     # "observation days:", obsDays, type(obsDays),
    #     "final interval length:", finalInt)
    # print(">>>>> likelihood: hatched:", hatched.sum())
    # if perfectInfo == 0:
        # numIntTotal = nestData[:,11]
        #print(">> number of observation intervals for each nest:\n", numIntTotal)
    # else:
        # numIntTotal = hatchTime # what was this for???
        # numIntTotal = 
        #print(">> number of observation intervals for each nest, frequency = 1 day:\n", numIntTotal)
    # NOTE this SHOULD be number of obs (I think) but it's giving number of 
    #      obs intervals, so I'll roll with it

    # ---------------------------------------------------------------------------------------------------
    # 2. Initialize the overall likelihood counter; Decimal gives more precision
    logLike = Decimal(0.0)         
    # ---------------------------------------------------------------------------------------------------
    # 3. Define state vectors (1x3 matrices) - all the possible nest states 
    #       [1 0 0] (alive)   [0 1 0] (fail-predation)   [0 0 1] (fail-flood) 
    stillAlive = np.array([1,0,0]) 
    mortFlood  = np.array([0,1,0])
    mortPred   = np.array([0,0,1])
    # ---------------------------------------------------------------------------------------------------
    # 4. Create arrays to hold state vectors for all nests:
    #    a. state of nest on date nest was first found (stateFF)
    #    b. state of nest on date nest was last checked (stateLC) - this is the fate as observed
    # Could also calculate bassed on nest fate value
    # FOR THE INITIAL STATE (stateFF), just one vector (see notebook) - later in code
    stateEnd    = np.empty((numNests, 3))     # state at end of normal interval
    stateLC     = np.empty((numNests, 3)) 
     # > use broadcasting - fill doesn't work with arrays as the fill value:
    stateEnd[:] = stillAlive
    stateLC[:]  = mortPred   # default is still depredation

    stateLC[flooded==True] = mortFlood  # flooded status gets flooded state vector
    stateLC[hatched==True] = stillAlive # hatched nests stay alive the entire time

    # print(">> state at the end of a normal interval:\n",stateEnd) 
    # print(">> state on final nest observation:\n", stateLC)
    # ---------------------------------------------------------------------------------------------------
    # 5. Compose the matrix equation for one observation interval.
    #        The formula used is from Etterson et al. (2007) 
    #    For this, you need the nest state at the beginning and end of the interval, plus interval length
    #    > intElt - length in days of the observation interval being assessed 
    #    > initial state (stateI) - state of the nest at the beginning of this interval 
    #    > stateF - state of the nes at the end of this interval
    #    There is a transition matrix that is multiplied for each day in the interval 
    #    > in this case, the nest started the interval alive and ended it alive as well 
    #    > daily nest probabilities: s - survival; mp - mortality from predation; mf - mortality from flooding 
    #    > these are daily probabilities, so raise transition matrix to the power of number of days in interval  
    #
    #                                     _         _  intElt           _   _ 
    #              [ 1 0 0 ]             |  s  0  0  |                 |  1  | 
    #                             *      |  mp 1  0  |            *    |  0  | 
    #                                    |_ mf 0  1 _|                 |_ 0 _|  
    #                              
    #      {  transpose(stateI) * trMatrix, raised to intElt power * stateF } 
    #
    # Then, you can multiply this equation times number of intervals (numIntTotal)
    #    Single in  interval --> all intervals --> likelihood

    # in the following code, we calculate all matrix multiplications for all nests, and then turn them on/off
    # based on whether we need them. Would it be faster to only do the multiplications we need? That would require indexing
    # maybe calculate the one for normal and storm intervals (not final interval) once, and the final intervals for each nest
    # (based on what the actual fate is)
    # > all nests start interval as alive - don't need stateFF
    # > transition matrix, from etterson 2007:
    trMatrix = np.array([[a_s,0,0], [a_mf,1,0], [a_mp,0,1]]) 
    # print(">> transition matrix\n:", trMatrix, trMatrix.shape)
    TstateI = np.transpose(stillAlive) 
    ##print(">> transpose of initial state vector:", TstateI, TstateI.shape)
    #                                             this will depend on how many storms/when they are
    # NOTE need the transpose of each row, not entire matrix. 
    pwr = np.linalg.matrix_power(trMatrix, obsFreq) # raise the matrix to the power of the number of days in obs int
    pwrStm = np.linalg.matrix_power(trMatrix, obsFreq*2) # power equation for storm intervals (longer obs int)
    # print(">> transition matrix raised to the obs int power:\n", pwr, pwr.shape)

    normalInt = stateEnd@pwr@TstateI
    # print(
    #     ">> likelihood of one normal interval:\n",
    #     normalInt,
    #     normalInt.shape,
    #     normalInt.dtype)

    # The final interval is one of these two (ends in final state):
    normalFinal = stateLC@pwr@TstateI
    # print("final interval:", normalFinal, "and -log likelihood:", -np.log(normalFinal))
    stormFinal  = stateLC@pwrStm@TstateI

    logLik   = np.ones(numNests, dtype=np.longdouble) # this should give it enough precision & avoid errors
    logLik      = logLik * np.log(normalInt) * -1 # dtype changes to float64 unless you multiply it by itself
    # print("-log likelihood of 1 interval:", logLik)
    
    # logLikFin    = np.ones(numNests, dtype=np.longdouble)
    logLikFin = np.empty(numNests, dtype=np.longdouble)
    # logLikFin.fill(-np.log(normalFinal))
    logLikFin= -np.log(normalFinal)
    logLikFin[hatched == True] = 0
    logLikFin[flooded == True] = -np.log(stormFinal[flooded==True])
    # logLikFin[]    = logLikFin * (-np.log(normalFinal))
    # logLikFinStm = logLikFinStm * (-np.log(stormFinal))
    # print(">> log likelihood final interval, updated with storms/hatch:\n", logLikFin)

    stormDuringFin = nestData[:,12] # was there a storm during the final interval?
    logLikelihood  = (logLik*numInt) + (logLikFin)
    # print("log likelihood of each nest history:", logLikelihood)

    # for x in range(numNests):
        # print(">> likelihood equation: (",logLik[x],"*",numIntNorm[x],")+(",logLikStm[x],"*",numIntStm[x],")+(",logLikFin[x],"**(1 -",stormDuringFin[x],")+(", logLikFinStm[x],"**",stormDuringFin[x])
    #    print(
    #        f">> likelihood equation for nest {x}: " 
    #     #    f"{numInt[x]:.0f} * {logLik[x]:.5f} + "
    #        f"({numInt[x]:.0f} * {logLik[x]:.5f}) + {logLikFin[x]:.5f} ="
    #     #    f"{logLikFinStm[x]:.5f} * (1-{stormDuringFin[x]:.0f}) + " 
    #     #    f"{logLikFinStm[x]:.2f} * {stormDuringFin[x]:.2f} = "
    #        f"{logLikelihood[x]:.2f}")


    logLike        = np.sum(logLikelihood)
    # print("overall log likelihood:", logLike)
    return(logLike)


# -----------------------------------------------------------------------------
#   THE LIKELIHOOD WRAPPER FUNCTION
# -----------------------------------------------------------------------------
# The values are log-transformed before running them thru the likelihood 
# function, so the values given to the optimizer are the untransformed values,
# meaning the optimizer output will also be untransformed.
# >> Therefore, need to transform the output as well.
def like_smd( 
        x, perfectInfo, hatchTime, nestData, obsFreq, 
        stormDays, surveyDays, whichRet):

    # unpack the initial values:
    s0   = x[0]
    mp0  = x[1]
    ss0  = x[2]
    mps0 = x[3]
    sM   = x[4]
    #@#print("initial values:", s0, mp0, ss0, mps0, sM)

    # transform the initial values so all are between 0 and 1:
    s1   = logistic(s0)
    mp1  = logistic(mp0)
    ss1  = logistic(ss0)
    mps1 = logistic(mps0)
    #@#print("logistic-transformed initial values:", s1, mp1, ss1, mps1)

    # further transform so they remain in lower left triangle:
    tri1 = triangle(s1, mp1)
    tri2 = triangle(ss1, mps1)
    s2   = tri1[0]
    mp2  = tri1[1]
    ss2  = tri2[0]
    mps2 = tri2[1]
    #@#print("triangle-transformed initial values:", s2, mp2, ss2, mps2)

    # compute the conditional probability of mortality due to flooding:
    mf2  = 1.0 - s2 - mp2
    mfs2 = 1.0 - ss2 - mps2

    numNests = nestData.shape[0]
    #@#print(">> number of nests:", numNests)

    # call the likelihood function:
    argL = np.array([s2,mp2,mf2,ss2,mps2,mfs2, sM])
    #ret = like(argL, ndata, obs, storm, survey)
    #ret = like(argL, nestData, obsFreq, stormDays, surveyDays)
    if whichRet == 1:
        ret = like(perfectInfo, hatchTime, argL, numNests, 
                   obsFreq,  nestData, surveyDays)
        # print('like_smd(): Msg : ret = ', ret)
    
    # else:
    elif whichRet == 2:
        ret = like_old(argL, obsFreq, nestData, surveyDays, stormDays)
        # print('like_smd(): Msg : using old function; ret = ', ret)
    
    else:
        print(' argument whichRet is invalid ')
    
    # rets = np.array([ret, ret2])

    return(ret)

# -----------------------------------------------------------------------------
#   CREATE NEST DATA AND RUN THE OPTIMIZER 
# -----------------------------------------------------------------------------
# need to loop through the param combinations
# within the loop, need to unpack the params and run the optimizer

#def run_params(paramsList, dirName):
# NOTE need to create storms after params are chosen
# NOTE is there even a reason to have this in a function?

with open(likeFile, "wb") as f:

    parID     = 0

    for i in range(0, len(paramsList)): # for each set of params
        #par = paramsList[i]
        # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print(
            "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ "
            "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ")
        print(">>>>>>>>> param set number:", parID)
        # paramsArray is a 2d array; loop through the rows 
        # each row is a set of parameters we are trying
        par        = paramsArray[i] 
        # in an array, all values have the same numpy dtype (float in this case) 
        # after selecting the row, unpack the params & change dtype as needed:
        numNests   = par[0].astype(int)
        pSurv      = par[1]
        pSurvStorm = par[2]
        freq       = par[4].astype(int)
        dur        = par[3].astype(int)
        hatchTime  = par[5].astype(int)
        obsFreq    = par[6].astype(int)
        pMFlood    = par[7]
        brDays     = par[8]
        pDisc      = par[9]
        fateCues   = 0.6 if obsFreq > 5 else 0.66 if obsFreq == 5 else 0.75
        # fateCues   = par[10]
        #survival   = par[15]
        # Generate random list of storm days based on real weekly probabilities
        #  > Then create a list of survey days 
        # NOTE do you want new storms/survey days for each replicate 
        #      or each parameter set?
        stormDays  = stormGen(freq, dur)
        surveyDays = mk_surveys(stormDays, obsFreq, brDays)
        repID      = 0  # keep trackof replicates
        numOut     = 21 # number of output params

        # Create a file to store the nest data as you go: nestfile
        now      = datetime.now().strftime("%H%M%S")
        nestfile = Path(uniquify(Path.home()/
                                'C://Users/sarah/Dropbox/nest_models/py_output'/
                                dirName/
                                ('nests'+now+'.npy')))
        nestfile.parent.mkdir(parents=True, exist_ok=True)
        columns = [
                'rep_ID','ID','initiation','end_date','true_fate','first_found',
                'last_active','last_observed','ass_fate','n_obs','storm_true'
                ]
        colnames = ', '.join([str(x) for x in columns])
        print(">> nest data file address:", nestfile)
        print(
            ">>> nest params in this set:", 
            pSurv, 
            pSurvStorm, 
            dur, 
            freq, 
            hatchTime, 
            obsFreq, 
            fateCues
            )

        with open(nestfile, "wb") as n:
            # append data from each repl to nestfile as you loop through them
            for r in range(nreps):
                # if you write likelihoods to a file as you go, shouldn't need 
                # an array to store them -
                # this should save memory, at least in this case?
                print(
                    ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
                    ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                print(">>>>>>>>>>>> replicate ID:", repID)
                likeVal =  np.zeros(shape=(numOut), dtype=np.longdouble)
                try:
                    nestData = mk_nests(par, initProb, stormDays, surveyDays)
                    np.save(n, nestData) # make sure this is correct kind of save
                except IndexError as error:
                    print(
                        ">> IndexError in nest data:", 
                        error,
                        ". Go to next replicate")
                    continue
                # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                hatched = nestData[:,3] # true hatched (storms accounted for)
                #hatchProp = sum(hatched)/numNests
                expDays = nestData[0,14]
                trueDSR = 1 - ((numNests - hatched.sum()) / expDays)
                print(
                    "> ALL NESTS - hatched:",
                    hatched.sum(),
                    "failed:", 
                    numNests - hatched.sum(), 
                    "and exposure days:", 
                    expDays
                    )
                # NOTE rearrange columns to make more sense? 
                # Keep only discovered nests, then count them:
                nestData   = nestData[np.where(nestData[:,6]==1)] 
                discovered = nestData.shape[0]   
                hatched2   = nestData[:,3]
                print(
                    ">> nests with only one observation while active:",
                    np.where(nestData[:,7] == nestData[:,8])
                    ) 
                exclude  = ((nestData[:,10] == 7) | 
                            (nestData[:,7]==nestData[:,8]))                         
                # if there's only one observation, firstFound will == lastActive
                excluded = sum(exclude) # exclude = boolean array; sum = num True
                analyzed = discovered - excluded
                failed   = nestData.shape[0] - sum(nestData[:,3])
                expDisc  = sum(nestData[:,15])
                trueDSR_disc = 1 - ((nestData.shape[0] - sum(nestData[:,3])) / 
                                      sum(nestData[:,15]))
                nestData    = nestData[~(exclude),:]    # remove excluded nests 
                print(
                    "> DISCOVERED NESTS - total | analyzed:", 
                    discovered, "|", analyzed,
                    # "excluded from analysis:", excluded,
                    "failed:", failed, "|",
                    nestData.shape[0] - sum(nestData[:,3]),
                    "exposure days:", expDisc, "|", sum(nestData[:,15])
                    )
                trueDSR_analysis = 1 - ((nestData.shape[0] - sum(nestData[:,3])) / 
                                          sum(nestData[:,15]) )
                z = randArgs()
                # NOTE: optimizer for matrix function takes an array (z) but 
                #       optimizer for MARK takes one param (srand) which just 
                #       happens to be in z
                srand = z[4]
                perfectInfo = 0
                ex = np.longdouble("0.0")

                # -------------------------------------------------------------
                # Run the optimizer on the likelihood function you choose:
                # -------------------------------------------------------------

                # Run the optimizer with messages for exceptions
                # -------------------------------------------------------------
                
                # ------- LIKE_SMD --------------------------------------------
                # print("main.py: Msg: Running optimizer") 
                try:
                    ans = optimize.minimize(
                            like_smd, z, 
                            args=(
                                perfectInfo, 
                                hatchTime, 
                                nestData, 
                                obsFreq, 
                                stormDays, 
                                surveyDays, 
                                1
                                ),
                            method='Nelder-Mead'
                            # method='SLSQP'
                            )
                    ex = 0.0
                except decimal.InvalidOperation as error2:
                    ex=100.0
                    print(
                        ">> Error: invalid operation in decimal:", 
                        error2, 
                        "Go to next replicate."
                        )
                    continue
                except OverflowError as error3:
                    ex=200.0
                    print(
                        ">> Error: overflow error:", 
                        error3, 
                        "Go to next replicate."
                        )
                    continue
                #
                # #@#print("Success?", ans.success, ans.message)

                # # ------- LIKE_SMD_OLD ---------------------------------------
                # print("main.py: Msg: Running optimizer on old function")
                # try:
                #     ans2 = optimize.minimize(
                #             like_smd, z, 
                #             #args=(nestData, obsInt, stormDays, surveyDays),
                #             args=(
                #                 perfectInfo, 
                #                 hatchTime, 
                #                 nestData, 
                #                 obsFreq, 
                #                 stormDays, 
                #                 surveyDays,
                #                 2
                #                 ),
                #             method='Nelder-Mead'
                #             # method='SLSQP'
                #             )
                #     ex = 0.0
                # except decimal.InvalidOperation as error2:
                #     ex=100.0
                #     print(
                #         ">> Error: invalid operation in decimal:", 
                #         error2, 
                #         "Go to next replicate."
                #         )
                #     continue
                # except OverflowError as error3:
                #     ex=200.0
                #     print(
                #         ">> Error: overflow error:", 
                #         error3, 
                #         "Go to next replicate."
                #         )
                #     continue

                # print("> answer:", ans.x)
                # print("> answer (old function):", ans2.x)
                res = ansTransform(ans)
                # res2 = ansTransform(ans2)
                # print("> res:\n", res)
                # print("> res (old function):\n", res2)
                print(
                    ">> true DSR of all nests:", 
                    trueDSR, 
                    "discovered nests:",
                    trueDSR_disc,
                    "and nests used in analysis:", 
                    trueDSR_analysis
                    )

                #OPTIMIZER: MARK function
                markAns  = optimize.minimize(
                        mark_wrapper, srand,
                        args=(nestData), 
                        method='Nelder-Mead'
                        )
                # did the optimizer converge?
                #@#print(">> success?", markAns.success, markAns.message) 
                # this MLE seems unlikely to need exceptions, but who knows...
                #NOTE markAns is an "OptimizeResult" object - does not match 
                #     the other objects in like_val; need to extract "x"
                if markAns.success == True:
                    # Transform the MARK optimizer output so that it is 
                    # between 0 and 1
                    mark_s = logistic(markAns.x[0])
                    # answer.x is a list itself - need to index
                    # print("> logistic of MARK answer:", mark_s)
                else:
                    mark_s = 10001

                # print(">>> MARK:", mark_s, "new:", res[0], "old:", res2[0])

                # Compile the optimizer output for each replicate with other important info 
                # including params and the number of nests actually used in the analysis
                # can get a proportion because we know how many nests were generated initially
                #firstPart = np.array([repID, mark_s])
                #secondPart = np.array([dur, freq, pSurv, pSurvStorm, pMFlood, hatchTime, numNests, obsFreq, discovered, excluded, ex])
                #like_val   = np.concatenate((firstPart, res, secondPart))

                # -------------------------------------------------------------
                # Decide which likelihood function output to use for results:
                # -------------------------------------------------------------
                # s = res2[0]
                s2,mp2,mf2,ss2,mps2,mfs2 = res # unpack like() function output
                # s2,mp2,mf2,ss2,mps2,mfs2 = res2 # unpack like_old() function output
                like_val = [
                        repID,mark_s,s2,mp2,mf2,ss2,mps2,mfs2,dur,freq,
                        trueDSR,trueDSR_analysis,pSurv, pSurvStorm,pMFlood,
                        # hatchTime,numNests,obsFreq,discovered,excluded,ex
                        hatchTime,numNests,obsFreq,discovered,excluded
                        ]
                # is there a reason this wasn't a list?
                print(">> like_val:\n", like_val)
                #, "lengths:", [len(x) for x in like_val])
                #like_val = np.array(like_val, dtype=np.float128)
                like_val = np.array(like_val, dtype=np.longdouble)
                
                #like_val   = np.array([
                #    repID,mark_s,s2,mp2,mf2,ss2,mps2,mfs2,dur,freq,pSurv,pSurvStorm, 
                #    pMFlood,hatchTime,numNests,obsFreq,discovered,excluded,ex
                #    ], dtype=np.float128) # NOTE problem could be multiple types in same array?

                    #repID,mark_s,s2,mp2,mf2,ss2,mps2,mfs2,stormDur,stormFreq,probSurv,SprobSurv, 
                    #probMortFlood,SprobMortFlood,hatchTime,numNests,obs_int,discovered,excluded,
                    #ex,stormMat
                #####################
                # this is being saved w/o the delimiter
                # stack exchange says to make it a list of only one item
                #np.savetxt(f, like_val, delimiter=",")
                column_names     = np.array([
                    'rep_ID', 'mark_s', 'psurv_est', 'ppred_est', 'pflood_est', 
                    'stormsurv_est', 'stormpred_est', 'stormflood_est', 'storm_dur', 
                    'storm_freq', 'psurv_real', 'psurv_found', 'psurv_given',
                    'stormsurv_given','pflood_given', 'hatch_time','num_nests',
                    # 'obs_int', 'num_discovered','num_excluded', 'exception'
                    'obs_int', 'num_discovered','num_excluded'
                    ])
                # header values need to be stored to np.ndarray for np.savetxt; actually, needs to be a string
                colnames = ', '.join([str(x) for x in column_names])
                #if repID == 0:
                #print(">> saving likelihood values")
                if parID == 0 | repID == 0: # only the first line gets the header
                    np.savetxt(f, [like_val], delimiter=",", header=colnames)
                    #np.savetxt(f, like_val, delimiter=",", header=colnames)
                    # print(">> saving likelihood values")
                else:
                    np.savetxt(f, [like_val], delimiter=",")
                    #np.savetxt(f, like_val, delimiter=",")
                    # print(">> saving likelihood values")
                # except this adds the header before every line
                # could just remove them later in R

                repID = repID + 1
                # print(">> repID increased:", repID)

        parID = parID +1
        # print(">> param set ID increased:", parID)
        

