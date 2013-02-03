import sys
import math
import numpy as np
import random
from scipy.stats import poisson

def weighted_choice(probs):
    random_pos = random.random() * sum(probs)
    current_pos = 0.0
    for i,p in enumerate(probs):
        current_pos += p
        if random_pos < current_pos:
            return i
    return None

class BikeNetwork:
    def __init__(self, ridesPerDay, nStations, loud=True):
        self.G = nx.DiGraph()
        self.times = np.zeros(shape=(nStations,nStations))
        self.stationList = []
        self.bikeList = []
        self.amenityList = []

        self.nRides = 0
        self.nStations = 0
        self.nDockFail = []
        self.nBikeFail = []
        self.nOrigin = []
        self.nDest = []

        self.balancerList = []
        self.maxBalancers = 1
        self.bTime = 0
        self.nBalances = 0
        self.dropOffs = []
        self.pickUps = []

        # standard timestep: 1 minute
        self.dt = 1
        self.time = 0

        self.ridesPerDay = ridesPerDay

        # ride likelihoods will be weighted by gaussian with
        # std. dev. given by tWeight
        self.tWeight = 20

        # Diagnostic variable: if True, we'll print out a load of statements
        # during a simulation.
        self.loud = loud

        # Which balance algorithm to use?
        self.algo = 0

        ## Read in trip-generation probabilities
        filename = "attr_simple.txt"
        data=np.genfromtxt(filename,unpack=True, skip_header=1)
        self.timeList = np.array(data[0])
        self.resO = np.array(data[1])
        self.resD = np.array(data[2])
        self.empO = np.array(data[3])
        self.empD = np.array(data[4])
        self.svcO = np.array(data[5])
        self.svcD = np.array(data[6])
        self.attract = zip(self.resO, self.resD, self.empO, self.empD, self.svcO, self.svcD)

    def smallReport(self):
        if (self.loud):
            print "t = %d" % (self.time)
            print "\tRiders = %d" % (len(self.bikeList))
            for i in range(self.nStations):
                stn = self.stationList[i]
                print "\tStation %d has %d bikes and %d docks." % (i, stn.bikes, stn.docks)

    def bigReport(self):
        if (self.loud):
            print "t = %d:" % (self.time)
            print "\tRiders = %d" % (len(self.bikeList))
            for i in range(self.nStations):
                stn = self.stationList[i]
                print "\tStation %d has %d bikes and %d docks." % (i, stn.bikes, stn.docks)
                
            print "\tTotal rides: %d" % self.nRides
            
            for i in range(self.nStations):
                print "\t\tStation %d: %d origins, %d destinations, %d dock fails, %d bike fails." % (i, self.nOrigin[i], self.nDest[i], self.nDockFail[i], self.nBikeFail[i])

            print "\tTotal rebalances: %d" % self.nBalances
            print "\tTotal dock fails: %d" % sum(self.nDockFail)
            print "\tTotal bike fails: %d" % sum(self.nBikeFail)
            print "\tTotal balancer hours/day: %f" % (24.0*float(self.bTime)/self.time)

        ## Calculate total capital costs per day
        self.cES = 30564
        self.cS = 39301
        self.cM = 48039
        self.cL = 56776
        self.cD = 709+307
        
        self.capitalCost = 0
        for i in range(self.nStations):
            nDocks = self.stationList[i].docks 
            if (nDocks>= 7 and nDocks < 11):
                self.capitalCost += self.cES + (nDocks-7)*self.cD
            elif (nDocks<15):
                self.capitalCost += self.cS + (nDocks-11)*self.cD
            elif (nDocks<19):
                self.capitalCost += self.cM + (nDocks-15)*self.cD
            elif (nDocks<=23):
                self.capitalCost += self.cL + (nDocks-19)*self.cD

        self.balancerCost = 25000
        self.capitalCost += self.maxBalancers*self.balancerCost
        self.COC = 0.03
        self.yrs = 8
        self.capitalDC = self.capitalCost * self.COC / ( 1 - (1+self.COC)**(-self.yrs) ) / 365

        ## Calculate total labor costs
        ## $20/hour, in dollars per minute:
        self.laborUnitCost = 0.3333
        self.laborDC = self.laborUnitCost*self.bTime / (self.time/1440)

        ## Calculate bike/dock fail costs
        ## Assume each failure is worth $10 in customer outrage
        self.failCost = 20.0
        self.failDC = self.failCost * (sum(self.nDockFail) + sum(self.nBikeFail)) / (self.time/1440)

        if (self.loud):
            print ""
            print "\tDaily capital cost: $%.02f" % (self.capitalDC)
            print "\tDaily labor cost: $%.02f" % (self.laborDC)
            print "\tDaily failure cost: $%.02f" % (self.failDC)

        return [self.capitalDC, self.laborDC, self.failDC]

    def addStation(self,amenity, bikes, docks):
        self.stationList.append(Station(amenity=amenity, bikes=bikes, docks=docks))
        self.nStations+=1
        self.amenityList.append(amenity)
        if (self.loud):
            print "New station with %d population, %d employment, %d services, %d bikes, %d docks!" % (amenity[0], amenity[1], amenity[2], bikes, docks)
        self.nBikeFail.append(0)
        self.nDockFail.append(0)
        self.nOrigin.append(0)
        self.nDest.append(0)

    def addBalancer(self,instr):
        tTot = [instr[0][2]]
        ## Change times to running totals as opposed to segment lengths
        for i, el in enumerate(instr):
            if (i >0):
                tTot.append(el[2] + tTot[i-1])
            instr[i][2] = tTot[i]

        self.balancerList.append(Balancer(instr))
        self.nBalances += 1
        if (self.loud):
            print "\tNew balancer:"
            for el in instr:
                if (el[1]<0):
                    print "\t\tPick up %d at station %d in %d minutes." % (-el[1], el[0], el[2])
                elif (el[1]>0):
                    print "\t\tDrop off %d at station %d in %d minutes." % (el[1], el[0], el[2])

    def addRide(self, o, d):
        if (self.stationList[o].bikes == 0):
            ## If the origin station has no bikes, fail out.
            self.nBikeFail[o] += 1
            if (self.loud):
                print "\tTime=%d: bike fail at station %d!" % (self.time,o)

        else:
            ## Get travel time between orig and dest
            tTot = self.times[o][d]

            ## Create new biker with destination d
            self.bikeList.append(Biker(tTot,d))
            self.nRides += 1
            
            ## Remove a bike from origin station
            self.stationList[o].bikes -= 1
            if (self.loud):
                print "Time=%d: new ride started at station %d with destination %d." % (self.time, o, d)
                print "%d bikes, %d docks remaining at station %d." % (self.stationList[o].bikes,self.stationList[o].docks, o)

            ## Record the ride
            self.nOrigin[o]+=1
            self.nDest[d]+=1
            


    def advance(self):
        # Generate n new riders according to poisson distribution
        # with default frequency.

        ## 1440 minutes in a day:
        day = int(self.time) / 1440

        ## Attractiveness is split into 30-minute blocks:
        halfHour = int(int(self.time) - day*1440) / 30
        totAttr = sum(self.attract[halfHour])

        ## Attractiveness table is normalized to give about 0.6 rides per day.
        ## Renormalize to the desired rides per day:
        totAttr *= self.ridesPerDay/0.59
        n = poisson.rvs(totAttr*self.dt)
        
        for i in range(n):
            # generate each rider
            o, d = self.generateRide(halfHour)
            self.addRide(o,d)
            
        # move all the bikes by dt
        for i, bike in enumerate(self.bikeList):
            self.bikeList[i].advance(self.dt)
            #if (self.loud):
                #print "Bike %d now %d minutes from station %d." % (i, self.bikeList[i].tTot-self.bikeList[i].time, self.bikeList[i].dest)



        for i, bal in enumerate(self.balancerList):
            self.balancerList[i].advance(self.dt)

        ## check for completions (bikes and rebalancers)
        self.collectAll()

        self.rebalance()

        self.time += self.dt

    def rebalance(self):
        ## Generate list of current occupancies and "needs" (distance
        ## from half full).
        occupList = []
        myDtype = [('station',int), ('occup', float), ('docks', int), ('need', int), ('bikes', int), ('free',int)]
        for i in range(self.nStations):
            stn = self.stationList[i]
            occup = float(stn.bikes)/stn.docks
            need = int((0.5-occup)*stn.docks)
            occupList.append((i, occup, stn.docks, need, stn.bikes, stn.docks-stn.bikes))
        
        if (self.algo==0):
            if (len(self.balancerList) < self.maxBalancers):
                ## If we have a rebalancer available:

                ## Indices which will reference highest-priority pickups/drops
                pickUp = -1
                dropOff = -1
                
                ## First scan for dock fails, to find highest priority pick-ups:

                ## These lines turn occupList into a numpy array, and then
                ## sort by free docks (asc), breaking ties on size (desc).
                occupArr = np.array(occupList, dtype=myDtype)
                ## sort by size, asc
                sizeSortA = np.sort(occupArr, order=['docks'])
                ## reverse
                sizeSortD = sizeSortA[::-1]
                ## sort by bikes, asc
                freeSortA = np.sort(sizeSortD, order=['free'])

                for el in freeSortA:
                    if (el['free'] > 1):
                        ## Only looking for stations with 0 or 1 free docks.
                        break
                    k = el['station']
                    if (k not in self.pickUps):
                        ## Not already scheduled
                        pickUp = k
                        break

                ## sort by occup, desc
                occupSortA = np.sort(sizeSortD, order=['occup'])
                occupSortD = occupSortA[::-1]

                for el in occupSortD:
                    k = el['station']
                    if (pickUp != -1):
                        ## Already got one to pick up.
                        break
                    if (el['occup'] < 0.8):
                        ## No need to rebalance.
                        break
                    if (k not in self.pickUps):
                        ## Not already scheduled
                        pickUp = k
                        break

                ## Now scan for bike fails, to find highest priority drop-offs:

                ## Sort by bikes (asc).
                bikeSortA = np.sort(sizeSortD, order=['bikes'])

                for el in bikeSortA:
                    if (el['bikes'] > 1):
                        ## Only looking for stations with 0 or 1 bikes.
                        break
                    k = el['station']
                    if (k not in self.dropOffs):
                        ## Not already scheduled
                        dropOff = k
                        break

                ## sort by occup, asc
                for el in occupSortA:
                    k = el['station']
                    if (dropOff != -1):
                        ## Already got one to pick up
                        break
                    if (el['occup'] > 0.2):
                        ## No need to rebalance.
                        break
                    if (k not in self.dropOffs):
                        ## Not already scheduled
                        dropOff = k
                        break

                if (pickUp != -1 and dropOff != -1):
                    PUstn = self.stationList[pickUp]
                    DOstn = self.stationList[dropOff]
                    
                    ## Leave at least 40% of docks full:
                    freeBikes = int(PUstn.bikes-0.4*PUstn.docks)
                    
                    ## Leave at least 40% of docks free:
                    needBikes = int(0.6*DOstn.docks - DOstn.bikes)
                    transfer = min(freeBikes,needBikes)

                    ## Create list of instructions: pickup first then dropoff
                    newInstrList = [[pickUp, -transfer, 10],
                                    [dropOff, transfer, 5+transfer]]

                    ## Send new balancer into the field
                    self.addBalancer(newInstrList)

                    ## Keep track of which stations are scheduled
                    self.pickUps.append(pickUp)
                    self.dropOffs.append(dropOff)
                    
                elif (pickUp != -1 and dropOff == -1):
                    PUstn = self.stationList[pickUp]
                    
                    ## Haven't got a drop-off station yet
                    ## Drop off at two stations with most free docks
                    freeSortA = np.sort(sizeSortD, order=['free'])
                    freeSortD = freeSortA[::-1]

                    dropOff1 = freeSortD[0]['station']
                    dropOff2 = freeSortD[1]['station']
                    
                    DOstn1 = self.stationList[dropOff1]
                    DOstn2 = self.stationList[dropOff2]

                    ## Leave at least 40% of docks full:
                    freeBikes = int(PUstn.bikes-0.4*PUstn.docks)

                    ## Split the bikes between the two stations
                    howMany = int(freeBikes/2)

                    ## Generate instruction list
                    newInstrList = [[pickUp, -howMany*2, 10],
                                    [dropOff1, howMany, 5+howMany*2],
                                    [dropOff2, howMany, 5+howMany]]

                    ## Send new balancer into the field
                    self.addBalancer(newInstrList)

                    ## Keep track of schedules
                    self.pickUps.append(pickUp)
                    self.dropOffs.append(dropOff1)
                    self.dropOffs.append(dropOff2)

                elif (pickUp == -1 and dropOff != -1):
                    DOstn = self.stationList[dropOff]
                    
                    ## Haven't got a pick-up station yet
                    ## Pick up from two stations with most free bikes
                    bikeSortA = np.sort(sizeSortD, order=['bikes'])
                    bikeSortD = bikeSortA[::-1]

                    pickUp1 = bikeSortD[0]['station']
                    pickUp2 = bikeSortD[1]['station']
                    
                    PUstn1 = self.stationList[pickUp1]
                    PUstn2 = self.stationList[pickUp2]

                    ## Leave at least 40% of docks free:
                    needBikes = int(0.6*DOstn.docks - DOstn.bikes)

                    ## Split the bikes between the two stations
                    howMany = int(needBikes/2)

                    ## Generate instruction list
                    newInstrList = [[pickUp1, -howMany, 10],
                                    [pickUp2, -howMany, 5+howMany],
                                    [dropOff, howMany*2, 5+howMany]]

                    ## Send new balancer into the field
                    self.addBalancer(newInstrList)

                    ## Keep track of schedules
                    self.dropOffs.append(dropOff)
                    self.pickUps.append(pickUp1)
                    self.pickUps.append(pickUp2)


    def generateRide(self, halfHour):
        ## Choose origin with probability weighted by "origin" amenity values.
        ## Amenity = pop*pop-based attractiveness
        temp = np.array(self.amenityList)
        pop = temp[:,0]
        emp = temp[:,1]
        svc = temp[:,2]
        totAttrO = pop*self.resO[halfHour] + emp*self.empO[halfHour] + svc*self.svcO[halfHour]
        #print totAttrO
        orig = weighted_choice(totAttrO)

        ## Generate probability list using "destination" attractiveness,
        ## as well as a Gaussian weighting on distance.
        probList = []
        for i in range(self.nStations):
            if (i == orig):
                ## Don't generate trip to origin station.
                probList.append(0)
                continue

            time = self.times[orig][i]
            amenity = self.amenityList[i][0]*self.resD[halfHour] + self.amenityList[i][1]*self.empD[halfHour] + self.amenityList[i][2]*self.svcD[halfHour]
            prob = amenity*math.exp(-time**2.0/(2*self.tWeight**2.0))

            probList.append(prob)

        dest = weighted_choice(probList)

        return orig, dest

    def collectAll(self):
        # "Cycle" (hah) through the bikes, and see if they've reached
        # their destinations.
        killList = []
        for i, bike in enumerate(self.bikeList):
            if (self.bikeList[i].reached_dest()):
                d = self.bikeList[i].dest
                dStn = self.stationList[d]
                if (dStn.bikes < dStn.docks):
                    killList.append(i)
                    dStn.bikes+=1
                    if (self.loud):
                        print "\tTime=%d: bike %d reached destination %d!" % (self.time,i, d)
                        print "\tStation %d now has %d bikes, %d docks." % (d, dStn.bikes, dStn.docks)
                else:
                    if (self.loud):
                        print "\tTime=%d: dock fail for bike %d at destination %d!" % (self.time, i, d)
                    self.nDockFail[d] += 1
                    
                    ## Generate new destination: nearest station.
                    ## Find the nearest station.
                    minTime = 10000000
                    which = 0
                    for j in range(self.nStations):
                        time = times[d][j]
                        if (j != d and time < minTime):
                            minTime=time
                            which = j

                    # Assign new destination
                    self.bikeList[i].dest = which
                    self.bikeList[i].tTot = minTime
                    if (self.loud):
                        print "\tNew destination for bike %d is station %d in %d minutes." % (i, which, minTime)

        # Delete all bikes which successfully reached destination.
        killList2 = np.array(killList)
        for i in range(len(killList2)):
            del self.bikeList[killList2[i]]
            killList2 = killList2-1

        # Now cycle through the rebalancers and do the same thing.
        killList = []
        for i, bal in enumerate(self.balancerList):
            if (self.balancerList[i].reached_dest()):
                # get the destination
                d = bal.instr[0][0]
                dStn = self.stationList[d]
                killList.append(i)

                ## This if statement is for pick-ups.
                if (bal.instr[0][1] < 0):
                    if (self.loud): print "\tRebalancer %d picking up from station %d" % (i, d)
                    toGet = -bal.instr[0][1]
                    freeBikes = dStn.bikes - 1

                    ## Remove from pickup list.
                    self.pickUps.remove(d)

                    ## Are there enough bikes for pickup?
                    if (freeBikes >= toGet):
                        ## Good: now pick up toGet bikes.
                        if (self.loud): print "\tPicking up %d bikes" % (toGet)
                        self.balancerList[i].bikes += toGet
                        dStn.bikes -= toGet

                    elif (freeBikes < toGet):
                        ## Bad: only pick up freeBikes bikes.
                        if (self.loud): print "\tWanted to pick up %d, but only %d bikes free." % (toGet, freeBikes)
                        self.balancerList[i].bikes += freeBikes
                        dStn.bikes -= freeBikes

                    if (self.loud): print "\tStation now has %d bikes, %d docks." % (dStn.bikes,dStn.docks)                        

                ## This if statement is for drop-offs.
                elif (bal.instr[0][1] > 0):
                    if (self.loud): print "\tRebalancer %d dropping off at station %d" % (i, d)
                    toDrop = bal.instr[0][1]
                    balBikes = bal.bikes
                    freeSpaces = dStn.docks - dStn.bikes - 1

                    ## Remove from dropoff list.
                    self.dropOffs.remove(d)

                    ## Are there enough bikes to drop?
                    if (balBikes < toDrop):
                        ## Bad news: not enough bikes to follow orders.
                        if (self.loud): print "\tOops: meant to drop off %d bikes but only have %d." % (toDrop, balBikes)
                        toDrop = balBikes

                    ## Are there enough free spaces?
                    if (freeSpaces >= toDrop):
                        ## Good: drop the bikes.
                        if (self.loud): print "\tDropping off %d bikes" % (toDrop)
                        self.balancerList[i].bikes -= toDrop
                        dStn.bikes += toDrop
                        
                    ## Not enough free space?
                    elif (freeSpaces < toDrop):
                        ## Uh-oh. Not enough spaces.
                        ## Drop only freeSpaces bikes and then generate new destination.
                        if (self.loud): print "\tWanted to drop off %d bikes, but only %d docks free." % (toDrop, freeSpaces)
                        self.balancerList[i].bikes -= freeSpaces
                        dStn.bikes += freeSpaces

                        ## Then generate new destination.
                        ## Does the rebalancer have other dropoff destinations thereafter?
                        if (len(bal.instr) > 1):
                            ## The assumption is that if there are more stations in the
                            ## instruction sequence, at least one must be a dropoff.
                            ## Otherwise, we would end up with bikes on the rebalancer
                            ## at the end of the sequence: bad.
                        
                            ## Add the extra bikes to the next dropoff.
                            for k in range(1,len(bal.instr)):
                                if (bal.instr[k][1] > 0):
                                    self.balancerList[i].instr[k][1] += self.balancerList[i].bikes
                                    if (self.loud): print "\tAdded %d bikes to next dropoff." % (self.balancerList[i].bikes)
                                    break
                        else:
                            ## If no, generate new destination: station with most free docks.
                            maxDocksFree = 0
                            which = 0
                            for k in range(self.nStations):
                                stn = self.stationList[k]
                                if (stn.docks - stn.bikes > maxDocksFree):
                                    maxDocksFree = stn.docks-stn.bikes
                                    which = k
                        
                            ## Travel time for rebalancer is always 5 minutes drive-time + 1 per bike
                            travelTime = 5 + self.balancerList[i].bikes
                        
                            ## Generate a new instruction to drop off remaining bikes:
                            newInstr = [which, self.balancerList[i].bikes, bal.time + travelTime]
                            if (self.loud): print "\tNew destination: drop off remaining %d bikes at station %d in %d minutes." % (bal.bikes, which, travelTime)
                            ## Append it to balancer's instruction sequence:
                            self.balancerList[i].instr.append(newInstr)
                            self.dropOffs.append(which)

        ## Pop all completed instructions, and if a rebalancer has
        ## no more (and no bikes), kill it
        killList2 = np.array(killList)
        for i in range(len(killList2)):
            if (len(self.balancerList[killList2[i]].instr) > 1):
                ## There are remaining instructions: just pop, don't kill.
                self.balancerList[killList2[i]].instr.pop(0)
            else:
                ## No remaining instructions: kill.
                if (self.balancerList[killList2[i]].bikes > 0):
                    print "..."
                    print "ERROR ERROR: rebalancer has bikes but no more instructions!"
                    print "Lost %d bike(s)!" % (self.balancerList[killList2[i]].bikes)
                    print "..."
                self.bTime += self.balancerList[killList2[i]].time
                del self.balancerList[killList2[i]]
                killList2 = killList2-1

class Biker:
    def __init__(self, tTot, dest):
        # tTot is the total time of the given path
        self.tTot = tTot

        # dest is the destination node
        self.dest = dest

        self.time = 0

    def advance(self,dt):
        self.time += dt

    def reached_dest(self):
        if (self.time >= self.tTot):
            # Time's up: we're at the destination station.
            return True

class Balancer:
    def __init__(self, instr):
        # bikes is the current number of bikes being hauled
        self.bikes = 0
        
        # instrStack consists of a list of lists. Each element
        # of this list is of the form [station, # bikes to drop off,
        # time]. A negative dropoff number represents pickups.

        self.instr = instr

        # The stack of instructions will be inflexibly followed,
        # unless the number of bikes to be picked up/ dropped off
        # can't be, due to full or empty stations upon arrival.

        self.time = 0

    def advance(self,dt):
        self.time += dt

    def reached_dest(self):
        if (self.time >= self.instr[0][2]):
            # Time's up: we're at the destination station.
            return True

class Station:
    def __init__(self, amenity, bikes, docks):
        self.amenity = amenity
        self.bikes = bikes
        self.docks = docks


if __name__=='__main__':
    random.seed(10)
    np.random.seed(5)

    ## Initialize network

    ## Constants for linear model relating crow-flies distance and
    ## elevation change (in meters) to travel time (in minutes):
    c1 = 0.0059
    c2 = 0.0328

    ## Get the station data:
    filename = "pilot_prop.txt"
    data=np.genfromtxt(filename,unpack=True, skip_header=1)
    lat = np.array(data[1])*math.pi/180
    lon = np.array(data[2])*math.pi/180
    height = np.array(data[3])
    pop = np.array(data[4])
    emp = np.array(data[5])
    svc = np.array(data[6]+data[7]+data[8]+data[9])
    attr = np.array(data[4]) + np.array(data[5])
    stations = map(int, data[10])
    stateInit = stations

    bikeDocks = [[4,8],[7,13],[10,18],[13,22]]
    
    nStations = len(lat)

    nSims = 10
    costList = []
    for i in range(nSims):
        print "sim %d" % (i+1)
        net = BikeNetwork(ridesPerDay=194, nStations=nStations, loud=False)
        for i in range(nStations):
            net.addStation(amenity=[int(pop[i]),int(emp[i]),int(svc[i])], bikes=int(bikeDocks[stateInit[i]][0]), docks=int(bikeDocks[stateInit[i]][1]))
    
        ## Calculate times between stations according to simple linear model:
        times = np.zeros(shape=(nStations,nStations))
        for i in range(nStations):
            for j in range(nStations):
                deltaEl = height[j]-height[i]
                
                ## Distance in meters between lat-long pair:
                dist = 1609*2*3950*math.asin(math.sqrt(math.sin(0.5*(lat[i]-lat[j]))**2+math.cos(lat[i])*math.cos(lat[j])*math.sin(0.5*(lon[i]-lon[j]))**2))
                
                times[i][j] = c1*dist + c2*deltaEl
            #print "i = %d, max_time = %f" % (i, max(times[i]))
        net.times = times
    
        ## Run the actual simulations.
        while (net.time < 14400):
            #net.smallReport()
            net.advance()
    
        cost = net.bigReport()
        costList.append(cost)
        print cost
        
    costArr = np.array(costList)
    capEDCs = costArr[:,0]
    labEDCs = costArr[:,1]
    failEDCs = costArr[:,2]
    totEDCs = capEDCs + labEDCs + failEDCs
    print "average capital cost = %f, stddev = %f" % (np.average(capEDCs), np.std(capEDCs))
    print "average labor cost = %f, stddev = %f" % (np.average(labEDCs), np.std(labEDCs))
    print "average failure cost = %f, stddev = %f" % (np.average(failEDCs), np.std(failEDCs))
    print "average total cost = %f, stddev = %f" % (np.average(totEDCs),np.std(totEDCs))
