from copy import deepcopy
from Clusterpy.core.toolboxes.cluster.componentsAlg.areacl import AreaCl
from Clusterpy.core.toolboxes.cluster.componentsAlg.selectionTypeFunctions import selectionTypeDispatcher
from os import getpid
from time import time
import numpy as np


#     rm = RegionMakerNodes(am, pRegions,
#                     initialSolution = initialSolution,
#                    distanceType = distanceType,
#                    distanceStat = distanceStat,
#                    selectionType = selectionType,
#                    objectiveFunctionType = objectiveFunctionType,
#                    Dij = Dij,
#                    Cio = Cio)
class RegionMakerNodes:
    """
    This class deals with a large amount of methods required during both the
    construction and local search phases. This class takes the area instances and
    coordinate them during the solution process. It also send information to
    Memory when needed.
    """
    def __init__(self, am, pRegions=2, initialSolution=[],
                 seedSelection = "kmeans",
                 distanceType = "EuclideanSquared",
                 distanceStat = "Functional",
                 selectionType = "Minimum",
                 alpha = 0.2,
                 numRegionsType = "Exogenous",
                 objectiveFunctionType = "Functional",
                 threshold = 0.0,
                 weightsDistanceStat = [],
                 weightsObjectiveFunctionType = [],
                 indexDataStat = [],
                 indexDataOF = [],
                 Dij = {},
                 Cio = {}):
        
        """
        @type am: AreaManager
        @param am: Area manager object.

        @type pRegions: integer
        @keyword pRegions: Number of regions in scheme

        @type seeds: list
        @keyword seeds: List of area IDs for initial seeds.

        @type distanceType: string
        @keyword distanceType: Type of distance to be used, by default "EuclideanSquared"

        @type distanceStat: string
        @keyword distanceStat: Type of conversion used for summarizing distance, by defaults "Average"

        @type selectionType: string
        @keyword selectionType: Type of selection criterion for construction phase, by defaults "Minimum"

        @type alpha: float.
        @keyword alpha: float equal or between the interval [0,1]; for GRASP selection only.

        @type numRegionsType: string
        @keyword numRegionsType: Type of constructive method (Exogenous, EndogenousThreshold,
        EndogenousRange), by default "Exogenous"

        @type objectiveFunctionType: string
        @keyword objectiveFunctionType: Method to calculate the objective function, by default "Total"

        @type threshold: float
        @keyword threshold: Minimum population threshold to be satisfied for each region

        @type weightsDistanceStat: list
        @keyword weightsDistanceStat:

        @type weightsObjectiveFunctionStat: list
        @keyword weightsObjectiveFunctionStat:

        @type indexDataStat = list
        @keyword indexDataStat:

        @type indexDataOf = list
        @keyword indexDataOf:
        """
        self.am = am
        self.areas = deepcopy(am.areas)
        self.distanceType = distanceType
        self.distanceStat = distanceStat
        self.weightsDistanceStat = weightsDistanceStat
        self.indexDataStat = indexDataStat
        self.weightsObjectiveFunctionType = weightsObjectiveFunctionType
        self.indexDataOF = indexDataOF
        self.selectionType = selectionType
        self.objectiveFunctionType = objectiveFunctionType
        self.n = len(self.areas)
        self.unassignedAreas = self.areas.keys()
        self.assignedAreas = []
        self.area2Region = {}
        self.region2Area = {}
        self.potentialRegions4Area = {}
        self.intraBorderingAreas = {}
        self.candidateInfo = {}
        self.externalNeighs = set()
        self.alpha = alpha
        self.numRegionsType = numRegionsType
        self.neighSolutions = {(0,0): 99999999}
        self.regionMoves = set()
        self.distances = Dij
        #self.distances = {}
        self.NRegion = []
        self.N = 0
        self.data = {}
        self.objInfo = -1
        self.assignAreasNoNeighs()
        self.seeds = []
        self.potentialNodes = [k for k, v in self.am.y.items() if int(v[0]) == 1.0]
        self.dynimizingNodesDistance = Cio

        #  PREDEFINED NUMBER OF REGIONS

        seeds = []
        regions2createKeys = []
        emptyList = []
        c = 0
        lenUnassAreas = len(self.unassignedAreas)
        s = 0
        i = 0
        lseeds = 0
        if numRegionsType == "Exogenous":
            if not initialSolution:
                self.pRegions = pRegions
                seeds = self.kmeansInit(Cio) # Aca hay que revisar porque en el caso que un nodo dinamizador este muy cerca de tpodos los demas, igual lo selecicona asi no sea de primero
                self.seeds = seeds
                self.setSeeds(seeds) #Listo
                c = 0
                while lenUnassAreas > 0:
                    self.constructRegions(Dij = Dij, Cio = Cio, seeds = seeds)
                    lenUnassAreas = len(self.unassignedAreas)
                    c += 1
                print( 'pasa constructRegions')
                print (seeds)
                self.objInfo = self.getObj(Dij, Cio, seeds)
                print( 'pasa selfobjinfo')
            else:
                uniqueInitSolution = set(initialSolution)
                self.pRegions = len(uniqueInitSolution)
                seeds = []
                for s in uniqueInitSolution:
                    seeds.append(initialSolution.index(s))
                self.setSeeds(seeds)
                regions2create = {}
                c = 0

                for i in initialSolution:
                    regions2create.setdefault(i, []).append(c)
                    c += 1
                c = 0
                regions2createKeys = regions2create.keys()
                for i in regions2createKeys:
                    self.unassignedAreas = regions2create[i][1:]
                    lenUnassAreas = len(self.unassignedAreas)
                    while lenUnassAreas > 0:
                        self.constructRegions(filteredCandidates=self.unassignedAreas,
                                              filteredReg=i, Dij = Dij, Cio = Cio)
                        lenUnassAreas = len(self.unassignedAreas)
                        c += 1
                self.objInfo = self.getObj(Dij, Cio, seeds)

                #  NUMBER OF REGIONS IS ENDOGENOUS WITH A THRESHOLD VALUE

        if self.numRegionsType == "EndogenousThreshold":
            self.constructionStage = "growing"
            try:
                self.areas[self.areas.keys()[0]].thresholdVar
            except:
                self.extractThresholdVar()
            self.regionalThreshold = threshold
            c = 0
            self.feasibleRegions = {}
            self.regionValue = {}
            seeds = []
            for aID in self.areas:
                if self.areas[aID].thresholdVar >= self.regionalThreshold:
                    seed = aID
                    seeds = seeds + [seed]
                    self.regionValue[c] = self.areas[seed].thresholdVar
                    self.feasibleRegions[c] = [seed]
                    self.removeRegionAsCandidate()
                    c += 1
            self.setSeeds(seeds)
            while len(self.unassignedAreas) != 0:
                np.random.shuffle(self.unassignedAreas)
                vals = []
                for index in self.unassignedAreas:
                    vals += [self.areas[index].thresholdVar]
                seed = self.unassignedAreas[0]
                self.setSeeds([seed], c)
                self.regionValue[c] = self.areas[seed].thresholdVar
                if self.regionValue[c] >= self.regionalThreshold:
                    self.feasibleRegions[c] = [seed]
                    self.removeRegionAsCandidate()
                    c += 1
                else:
                    feasibleThreshold = 1
                    while self.regionValue[c] < self.regionalThreshold:
                        self.addedArea = -1
                        try:
                            self.constructRegions()
                            self.regionValue[c] += self.areas[self.addedArea].thresholdVar
                        except:
                            feasibleThreshold = 0
                            break
                    if feasibleThreshold == 1:
                        self.feasibleRegions[c] = self.region2Area[c]
                        self.removeRegionAsCandidate()
                    c += 1

        #  NUMBER OF REGIONS IS ENDOGENOUS WITH A RANGE VALUE

        if self.numRegionsType == "EndogenousRange":
            self.constructionStage = "growing"  #  there are two values for constructionStage: "growing" and "enclaves"
            try:
                self.areas[self.areas.keys()[0]].thresholdVar
            except:
                self.extractThresholdVar()
            self.regionalThreshold = threshold
            c = 0
            self.feasibleRegions = {}
            while len(self.unassignedAreas) != 0:

                #  select seed

                np.random.shuffle(self.unassignedAreas)
                seed = self.unassignedAreas[0]
                self.setSeeds([seed],c)

                #  regionRange contains the current range per region
                #  regionalThreshold is the predefined threshold value

                self.regionRange = {}
                maxValue = self.areas[seed].thresholdVar
                minValue = self.areas[seed].thresholdVar
                currentRange = maxValue - minValue
                self.regionRange[c] = currentRange

                # grow region if possible

                stop = 0
                while stop == 0:
                    upplim = maxValue + self.regionalThreshold - currentRange
                    lowlim = minValue - self.regionalThreshold + currentRange
                    feasibleNeigh = 0
                    toRemove = []
                    for ext in self.externalNeighs:
                        if self.areas[ext].thresholdVar <= upplim and self.areas[ext].thresholdVar >= lowlim:
                            feasibleNeigh = 1
                        if self.areas[ext].thresholdVar > upplim or self.areas[ext].thresholdVar < lowlim:
                            toRemove.append(ext)
                    self.toRemove = toRemove
                    if feasibleNeigh == 0:
                        stop = 1
                    if feasibleNeigh == 1:
                        try:
                            self.constructRegions()
                            if self.areas[self.addedArea].thresholdVar > maxValue:
                                maxValue = self.areas[self.addedArea].thresholdVar
                            if self.areas[self.addedArea].thresholdVar < minValue:
                                minValue = self.areas[self.addedArea].thresholdVar
                            currentRange = maxValue - minValue
                            self.regionRange[c] = currentRange
                        except:
                            stop = 1
                self.feasibleRegions[c] = self.region2Area[c]
                self.removeRegionAsCandidate()
                c += 1
        self.getIntraBorderingAreas()
    
    def kmeansInit(self, Cio):
        cachedDistances = {}
        ys = Cio
        n = int(np.sum(list(self.am.y.values())))
        #print( 'n', n) SE EJECUTA
        #n = 5
        #print 'pasa'
        #for k,v in self.am.y.iteritems():
        #    print v[0]
        #    if int(v[0]) == 1.0:
        #        print v[0] 
        potential = [k for k,v in self.am.y.items() if int(v[0]) == 1.0]
        
        #print 'Potential Nodes', potential
        distances = np.ones(n)
        #print 'distances', distances
        total = sum(distances)
        #print 'TOTAL', total
        #print 'AQUI'
        probabilities = list(map(lambda x: x / float(total), distances))
        seeds = []
        localDistanceType = self.distanceType  # ESTO CMABIARA
        
        returnDistance2Area = AreaCl.returnDistance2Area
        np.random.seed(int(time() * getpid()) % 4294967295)
        for k in range(self.pRegions):
            
            random = np.random.uniform(0, 1)
            find = False
            acum = 0
            cont = 0
            while not find:
                inf = acum
                sup = acum + probabilities[cont]
                if inf <= random <= sup:
                    find = True
                    seeds += [potential[cont]]
                    #print 'SSSSSeeds', seeds
                    #print 'SelfAMAreas Old', self.am.areas
                    selfAmAreas = { key: self.am.areas[key] for key in potential }
                    #print 'SelfAMAreas new', selfAmAreas
                    #selfAmAreas = self.am.areas #Estas

                    for area in selfAmAreas:
                        currentArea = selfAmAreas[area] # Esta
                        tempMap = []
                        
                        for x in seeds:
                            if x < area:
                                k = (x, area)
                            elif x > area:
                                k = (area, x)
                            else:
                                k = (0,0)
                            cached = cachedDistances.get(k, -1)
                            
                            if cached < 0:
                                '''
                                newDist = returnDistance2Area(currentArea,
                                                               selfAmAreas[x],
                                                               distanceType = localDistanceType)
                                '''
                                #print 'Current Area', area
                                #print 'Other Area', x
                                newDist = Cio[area, x]
                                #print 'newDist', newDist
                                tempMap.append(newDist)
                                cachedDistances[k] = newDist

                            else:
                                #print 'cachedddddddd'
                                tempMap.append(cached)
                                
                        #print 'Seeds', seeds
                        distancei = min(tempMap)
                        #print 'dis'
                        #print 'INDEX', potential.index(area)
                        distances[potential.index(area)] = distancei
                    total = sum(distances)
                    probabilities = list(map(lambda x: x / float(total), distances))
                else:
                    cont += 1
                    acum = sup
        del cachedDistances
        #print 'Seeds', seeds
        #print 'pasa por aca tambien'
        return seeds

    def constructRegions(self, filteredCandidates=-99, filteredReg=-99, Dij = {}, Cio = {}, seeds = []):
        """
        Construct potential regions per area
        """
        _d_stat = self.distanceStat
        _wd_stat = self.weightsDistanceStat
        _ida_stat = self.indexDataStat
        _fun_am_d2r = self.am.getDistance2Region

        lastRegion = 0
        #print 'Potential Regions 4 Area', self.potentialRegions4Area.keys()
        #print 'Potential Regions 4 Area', self.potentialRegions4Area.values()
        for areaID in self.potentialRegions4Area.keys():
            if len(self.areas)<areaID:
                print ('Problemas con el tamaño', len(self.areas), areaID)
            a = self.areas[areaID]
            regionIDs = list(self.potentialRegions4Area[areaID])
            for region in regionIDs:
                if (self.numRegionsType != "Exogenous" and
                    self.constructionStage == "growing"
                    and region in self.feasibleRegions):
                    #  Once a region reaches the threshold, the grow is
                    #  rejected until the assignation of enclaves
                    continue
                else:
                    if filteredCandidates == -99:
                        if (areaID not in self.newExternal and
                            region != self.changedRegion):
                            lastRegion = region
                            #print "me metí acá"
                            #print meRompo
                            pass
                        else:
                            _reg_dist = 0.0
                            if self.selectionType != "FullRandom":
                                _reg_dist = _fun_am_d2r(self.areas[areaID],
                                                        self.region2Area[region],
                                                        distanceStat = _d_stat,
                                                        weights = _wd_stat,
                                                        indexData = _ida_stat,
                                                        Dij = Dij,
                                                        Cio = Cio,
                                                        seeds = seeds)
                            self.candidateInfo[(areaID, region)] = _reg_dist

                    elif (filteredCandidates != -99 and
                          areaID in filteredCandidates and
                          region == filteredReg):
                        _reg_dist = _fun_am_d2r(self.areas[areaID],
                                                self.region2Area[region],
                                                distanceStat = _d_stat,
                                                weights = _wd_stat,
                                                indexData = _ida_stat)
                        self.candidateInfo[(areaID, region)] = _reg_dist

        if len(self.candidateInfo) == 0:
            self.changedRegion = lastRegion
        if self.numRegionsType == "EndogenousRange":
            self.filterCandidate(self.toRemove)
        selectionTypeDispatcher[self.selectionType](self)

    def assignAreasNoNeighs(self):
        """
        Assign to the region "-1" for the areas without neighbours
        """
        noNeighs = list(self.am.noNeighs)
        nr = -1
        for areaID in noNeighs:
            self.area2Region[areaID] = nr
            try:
                aid = self.unassignedAreas.remove(areaID)
            except:
                pass
            self.assignedAreas.append(areaID)
            setAssigned = set(self.assignedAreas)
        nr = nr - 1

    def setSeeds(self, seeds, c=0):
        '''
        Sets the initial seeds for clustering
        '''

        if self.numRegionsType == "Exogenous" and len(seeds) <= self.pRegions:
            idx = range(self.n)
            didx = list((set(idx) - set(seeds)) - self.am.noNeighs)
            np.random.shuffle(didx)
            self.seeds = seeds + didx[0:(self.pRegions - len(seeds))]
        else:
            self.seeds = seeds
            for seed in self.seeds:
                self.NRegion += [0]
                self.assignSeeds(seed, c)
                c += 1

    def extractThresholdVar(self):
        """
        Separate aggregation variables (data) from the variable selected
        to satisfy a threshold value (thresholdVar)
        """
        self.totalThresholdVar = 0.0
        for areaId in self.areas.keys():
            self.areas[areaId].thresholdVar = self.areas[areaId].data[-1]
            self.areas[areaId].data = self.areas[areaId].data[0: -1]
            self.totalThresholdVar += self.areas[areaId].thresholdVar

    def removeRegionAsCandidate(self):
        """
        Remove a region from candidates
        """
        for i in self.candidateInfo.keys():
          a, r = i
          if r in self.feasibleRegions:
            self.candidateInfo.pop(i)

    def recoverFromExtendedMemory(self, extendedMemory):
        """
        Recover a solution form the extended memory
        """
        self.objInfo = extendedMemory.objInfo
        self.area2Region = extendedMemory.area2Region
        self.region2Area = extendedMemory.region2Area
        self.intraBorderingAreas = extendedMemory.intraBorderingAreas
        self.seeds = extendedMemory.seeds