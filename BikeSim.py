import networkx as nx
import sys
import math
import numpy as np
from scipy.stats import poisson

if __name=='__main__':
    
class BikeNetwork:
    def __init__(self):
        self.G = nx.Graph()
        self.G.add_node(amenity=1000, bikes=10, docks=20)
        self.G.add_node(amenity=2000, bikes=10, docks=20)
        self.G.add_node(amenity=3000, bikes=10, docks=20)
        self.G.add_edge(1,2,tAB=5,tBA=5)
        self.G.add_edge(1,3,tAB=10,tBA=15)
        self.G.add_edge(2,3,tAB=5,tBA=10)
        self.bikeList = []
        self.amenityList = []
        self.nRides = 0
        self.nDockFail = 0
        self.nBikeFail = 0

        # standard timestep: 1 minute
        self.dt = 1

        # start at 100 rides per day (expressed on per-minute basis)
        self.frequency = 0.069

    def addStation(self,amenity, bikes, docks):
        self.G.add_node

    def advance(self):
        # Generate new riders according to poisson distribution
        # with some frequency.
        n = poisson.rvs(self.frequency*self.dt)
        for i in range(n):
            self.generateRide()
        
        # move all the bikes by dt
        for i, bike in enumerate(self.bikeList):
            self.bikeList[i].advance(self.dt)

        # check for completions
        self.collectAll

    def generateRide(self):
        # Choose origin and destination with probabilities
        # weighted by amenity values.
        

    def collectAll(self):
        killList = []
        for i, bike in enumerate(self.bikeList):
            if (self.bikeList[i].reached_dest):
                dest = self.bikeList[i].dest
                if (self.G[dest].bikes < self.G[dest].docks):
                    killList.append(i)
                    self.G[dest].bikes+=1
                else:
                    print "dock fail!"
                    self.nDockFail += 1
                    # Eventually, generate new destination: nearest
                    # station.

                    # Right now, just kill the ride.
                    killList.append(i)

        killList2 = np.array(killList)
        for i in range(len(killList2)):
            del self.bikeList[killList2[i]]
            killList2 = killList2-1
                
class Biker:
    def __init__(self, t_tot, dest):
        # t_tot is the total time of the given path
        self.t_tot = t_tot

        # dest is the destination node
        self.dest = dest

        self.time = 0

    def advance(self,dt):
        self.time += dt
        print

    def reached_dest(self):
        if (self.time > self.t_tot):
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
