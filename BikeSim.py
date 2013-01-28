import networkx as nx
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
    def __init__(self):
        self.G = nx.DiGraph()
        self.bikeList = []
        self.amenityList = []

        self.nRides = 0
        self.nStations = 0
        self.nDockFail = 0
        self.nBikeFail = 0

        # standard timestep: 1 minute
        self.dt = 1
        self.time = 0

        # start at 100 rides per day (expressed on per-minute basis)
        self.frequency = 0.069

        # ride likelihoods will be weighted by gaussian with
        # std. dev. given by tWeight
        self.tWeight = 20

        # Diagnostic variable: if True, we'll print out a load of statements
        # during a simulation.
        self.loud = True

        # Which balance algorithm to use?
        self.algo = 0

        self.addStation(amenity=1000, bikes=10, docks=20)
        self.addStation(amenity=2000, bikes=10, docks=20)
        self.addStation(amenity=3000, bikes=10, docks=20)
        self.G.add_edge(0,1,time=5)
        self.G.add_edge(1,0,time=5)
        self.G.add_edge(0,2,time=10)
        self.G.add_edge(2,0,time=15) 
        self.G.add_edge(1,2,time=5)
        self.G.add_edge(2,1,time=10)

        self.kioskCost = 17000
        self.bikeCost = 1000
        self.dockCost = 1500
        self.capitalCost = 0
        for i in range(self.nStations):
            stn = self.G.node[i]
            bikes = stn['bikes']
            docks = stn['docks']
            cost = self.kioskCost + bikes*self.bikeCost + docks*self.dockCost
            self.capitalCost += cost
        self.COC = 0.03
        self.yrs = 8
        self.capitalEDC = self.capitalCost * self.COC / ( 1 - (1+self.COC)**(-self.yrs) ) / 365

    def smallReport(self):
        print "t = %d:" % (self.time)
        print "\tRiders = %d" % (len(self.bikeList))
        for i in range(self.nStations):
            stn = self.G.node[i]
            print "\tStation %d has %d bikes and %d docks." % (i, stn['bikes'], stn['docks'])

    def bigReport(self):
        print "t = %d:" % (self.time)
        print "\tRiders = %d" % (len(self.bikeList))
        for i in range(self.nStations):
            stn = self.G.node[i]
            print "\tStation %d has %d bikes and %d docks." % (i, stn['bikes'], stn['docks'])
        print "\tTotal dock fails: %d" % self.nDockFail
        print "\tTotal bike fails: %d" % self.nBikeFail
        print "\tCapital EDC: $%.02f" % self.capitalEDC

    def addStation(self,amenity, bikes, docks):
        self.G.add_node(self.nStations, amenity=amenity, bikes=bikes, docks=docks)
        self.nStations+=1
        self.amenityList.append(amenity)
        if (self.loud):
            print "New station with %d amenity, %d bikes, %d docks!" % (amenity, bikes, docks)

    def addRide(self, o, d):
        if (self.G.node[o]['bikes'] == 0):
            # if the origin station has no bikes, fail out
            self.nBikeFail += 1
            if (self.loud):
                print "\tTime=%d: bike fail at station %d!" % (self.time,o)

        else:
            # get travel time between orig and dest
            tTot = self.G[o][d]['time']
            
            # create new biker with destination d
            self.bikeList.append(Biker(tTot,d))
            self.nRides += 1
            
            # remove a bike from origin station
            self.G.node[o]['bikes'] -= 1
            if (self.loud):
                print "Time=%d: new ride started at station %d with destination %d." % (self.time, o, d)
                print "%d bikes, %d docks remaining at station %d." % (self.G.node[o]['bikes'],self.G.node[o]['docks'], o)

    def advance(self):
        # Generate n new riders according to poisson distribution
        # with default frequency.
        n = poisson.rvs(self.frequency*self.dt)
        
        for i in range(n):
            # generate each rider
            o, d = self.generateRide()
            self.addRide(o,d)
            
        # move all the bikes by dt
        for i, bike in enumerate(self.bikeList):
            self.bikeList[i].advance(self.dt)
            #if (self.loud):
                #print "Bike %d now %d minutes from station %d." % (i, self.bikeList[i].tTot-self.bikeList[i].time, self.bikeList[i].dest)

        # check for completions
        self.collectAll()

        self.rebalance()

        self.time += self.dt

    def rebalance(self):
        occup_list = []
        if (self.algo==0):
            for i in range(self.nStations):
                stn = self.G.node[i]
                occup = stn['bikes']/stn['docks']
        return

    def generateRide(self):
        # Choose origin with probability
        # weighted by amenity values.
        orig = weighted_choice(self.amenityList)

        # distance_list = list of distances to all other stations
        probList = []
        for i in range(self.nStations):
            if (self.G.has_edge(orig,i)):
                time = self.G[orig][i]['time']
                amenity = self.G.node[i]['amenity']
                prob = amenity*math.exp(-time**2.0/(2*self.tWeight**2.0))
            else:
                prob = 0
            probList.append(prob)

        dest = weighted_choice(probList)

        return orig, dest

    def collectAll(self):
        killList = []
        for i, bike in enumerate(self.bikeList):
            if (self.bikeList[i].reached_dest()):
                d = self.bikeList[i].dest
                dStn = self.G.node[d]
                if (dStn['bikes'] < dStn['docks']):
                    killList.append(i)
                    dStn['bikes']+=1
                    if (self.loud):
                        print "Time=%d: bike %d reached destination %d!" % (self.time,i, d)
                        print "Station %d now has %d bikes, %d docks." % (d, dStn['bikes'], dStn['docks'])
                else:
                    if (self.loud):
                        print "\tTime=%d: dock fail for bike %d at destination %d!" % (self.time, i, d)
                    self.nDockFail += 1
                    
                    # Generate new destination: nearest
                    # station.

                    # list of times to other stations
                    timeDict = nx.shortest_path_length(self.G,source=d)

                    # Find the shortest
                    minTime = 1000000
                    which = 0
                    for j in range(self.nStations):
                        time = timeDict[j]
                        if (time > 0 and time < minTime):
                            minTime=time
                            which = j

                    # Assign new destination
                    self.bikeList[i].dest = which
                    self.bikeList[i].tTot = minTime
                    if (self.loud):
                        print "New destination is station %d in %d minutes." % (which, minTime)

        # Delete all bikes which successfully reached destination.
        killList2 = np.array(killList)
        for i in range(len(killList2)):
            del self.bikeList[killList2[i]]
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
        
class Station:
    def __init__(self, amenity):
        # Amenity is (something like) the number of 
        # residential/commercial uses within the geographic
        # range of the station.

        # It's used to generate the number of rides at each station
        # (or more accurately the probability of any given station
        # being the start or endpoint of a ride)
        
        self.amenity = amenity

class Path:
    def __init__(self, tAB, tBA):
        # tAB and tBA are the times required to go from station A
        # to B and B to A respectively. Typically these will be
        # symmetric but not if (e.g.) a hill is involved.

        # Thinking of times as measured in minutes right now.
        
        self.tAB = tAB
        self.tBA = tBA


if __name__=='__main__':
    random.seed(3)
    np.random.seed(3)
    net = BikeNetwork()
    while (net.time < 1100):
        net.smallReport()
        net.advance()

    net.bigReport()
