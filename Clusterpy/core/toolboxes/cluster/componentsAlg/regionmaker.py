from copy import deepcopy
import numpy as np

class RegionMaker:
    """
    This class deals with a large amount of methods required during both the
    construction and local search phases. This class takes the area instances and
    coordinate them during the solution process. It also send information to
    Memory when needed.
    """
    def __init__(self, am, pRegions=2, initialSolution=[],
                 seedSelection = "kmeans",
                 distanceType = "EuclideanSquared",
                 distanceStat = "Centroid",
                 selectionType = "Minimum",
                 alpha = 0.2,
                 numRegionsType = "Exogenous",
                 objectiveFunctionType = "SS",
                 threshold = 0.0,
                 weightsDistanceStat = [],
                 weightsObjectiveFunctionType = [],
                 indexDataStat = [],
                 indexDataOF = []):
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
        self.neighSolutions = {(0,0): 9999}
        self.regionMoves = set()
        self.distances = {}
        self.NRegion = []
        self.N = 0
        self.data = {}
        self.objInfo = -1
        self.assignAreasNoNeighs()

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
                seeds = self.kmeansInit()
                self.setSeeds(seeds)
                c = 0
                while lenUnassAreas > 0:
                    self.constructRegions()
                    lenUnassAreas = len(self.unassignedAreas)
                    c += 1
                self.objInfo = self.getObj()
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
                                              filteredReg=i)
                        lenUnassAreas = len(self.unassignedAreas)
                        c += 1
                self.objInfo = self.getObj()

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

    def returnRegions(self):
        """
        Return regions created
        """
        areasId = self.area2Region.keys()
        areasId = np.sort(areasId).tolist()
        return [self.area2Region[area] for area in areasId]