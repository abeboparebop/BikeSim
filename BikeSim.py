import networkx as nx
import sys
import math
import numpy as np
from scipy.stats import poisson

if __name=='__main__':

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
        self.addStation(amenity=1000, bikes=10, docks=20)
        self.addStation(amenity=2000, bikes=10, docks=20)
        self.addStation(amenity=3000, bikes=10, docks=20)
        self.G.add_edge(1,2,tAB=5,tBA=5)
        self.G.add_edge(1,3,tAB=10,tBA=15)
        self.G.add_edge(2,3,tAB=5,tBA=10)
        self.bikeList = []
        self.amenityList = []
        self.nRides = 0
        self.nStations = 0
        self.nDockFail = 0
        self.nBikeFail = 0

        # standard timestep: 1 minute
        self.dt = 1

        # start at 100 rides per day (expressed on per-minute basis)
        self.frequency = 0.069

        # ride likelihoods will be weighted by gaussian with
        # std. dev. given by tWeight
        self.tWeight = 20

    def addStation(self,amenity, bikes, docks):
        self.G.add_node(self.nStations, 'amenity'=amenity, 'bikes'=bikes, 'docks'=docks)
        self.nStations+=1
        self.amenityList.append(amenity)

    def addRide(self, o, d):
        if (self.G[o]['bikes'] == 0):
            # if the origin station has no bikes, fail out
            print "bike fail!"
            self.nBikeFail += 1
            
        else:
            # get travel time between orig and dest
            tTot = self.G[o][d]['time']
            
            # create new biker with destination d
            self.bikeList.append(Biker(tTot,d))
            
            # remove a bike from origin station
            self.G[o]['bikes'] -= 1

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

        # check for completions
        self.collectAll()

    def generateRide(self):
        # Choose origin with probability
        # weighted by amenity values.
        orig = weighted_choice(amenityList)

        # distance_list = list of distances to all other stations
        probList = []
        for i in range(self.nStations):
            if (self.G.has_edge(orig,i)):
                time = self.G[orig][i]['time']
                amenity = self.G[i]['amenity']
                prob = amenity*math.exp(-time**2.0/(2*self.tWeight**2.0))
            else:
                prob = 0
            probList.append(prob)

        dest = weighted_choice(probList)

        return orig, dest

    def collectAll(self):
        killList = []
        for i, bike in enumerate(self.bikeList):
            if (self.bikeList[i].reached_dest):
                d = self.bikeList[i].dest
                if (self.G[d]['bikes'] < self.G[d]['docks']):
                    killList.append(i)
                    self.G[d]['bikes']+=1
                    print "a bike reached destination %d!" % d
                else:
                    print "dock fail!"
                    self.nDockFail += 1
                    
                    # Generate new destination: nearest
                    # station.

                    # list of times to other stations
                    timeDict = nx.shortest_path_length(G,source=d)

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
        print

    def reached_dest(self):
        if (self.time > self.tTot):
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
