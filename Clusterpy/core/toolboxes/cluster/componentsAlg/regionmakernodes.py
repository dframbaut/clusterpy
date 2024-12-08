from copy import deepcopy
from Clusterpy.core.toolboxes.cluster.componentsAlg.areacl import AreaCl
from Clusterpy.core.toolboxes.cluster.componentsAlg.selectionTypeFunctions import selectionTypeDispatcher
from Clusterpy.core.toolboxes.cluster.componentsAlg.objFunctions import makeObjDict, objectiveFunctionTypeDispatcher
from os import getpid
from time import time
import numpy as np
from tqdm import tqdm 
import gc


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
                 Cio = {},
                 centroids = {}):
        
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

        self.unassignedAreas = list(self.areas.keys())

        self.assignedAreas = []
        self.area2Region = {}
        self.region2Area = {}
        # This is empty
        self.centroids = centroids
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
        self.potentialNodes = [k for k, v in self.am.y.items() if int(v[-1]) == 1.0]
        self.dynimizingNodesDistance = Cio

        #  PREDEFINED NUMBER OF REGIONS

        seeds = []
        regions2createKeys = []
        emptyList = []
        c = 0
        #print('unassignedAreas: ', self.unassignedAreas)
        lenUnassAreas = len(self.unassignedAreas)
        s = 0
        i = 0
        lseeds = 0
        if numRegionsType == "Exogenous":
            if not initialSolution:
                self.pRegions = pRegions # pRegions = 2
                seeds = self.kmeansInit(Cio) # Aca hay que revisar porque en el caso que un nodo dinamizador este muy cerca de todos los demas, igual lo selecciona asi no sea de primero
                self.seeds = seeds
                self.setSeeds(seeds) #Listo
                #print("Se ejecutó setSeeds")
                c = 0

                # ISSUES with lenUnassAreas
                pbar = tqdm(total=lenUnassAreas, desc="Constructing Regions", unit="area")
                # ====================================================================
                print('lenUnassAreas: ', lenUnassAreas)
                while lenUnassAreas > 0: #lenUnassAreas = 125
                    self.constructRegions(Dij = Dij, Cio = Cio, seeds = seeds)
                    lenUnassAreas = len(self.unassignedAreas)
                    pbar.update(1)
                    c += 1
                # ====================================================================
                pbar.close()
                print( 'pasa constructRegions')
                #print (seeds)

                # print('Arguments in colombia')
                # for attr_name in dir(self):
                #     if not attr_name.startswith("__"):
                #         try:
                #             attr_value = getattr(self, attr_name)
                #             print(f"Attribute {attr_name}: {attr_value}", '\n')
                #         except AttributeError:
                #             print(f"Attribute {attr_name} could not be accessed.")
    
                self.objInfo = self.getObj(Dij, Cio, seeds)
                #print( 'pasa selfobjinfo')
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
        self.getIntraBorderingAreas()
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

    def getIntraBorderingAreas(self):
        """
        Gets the intrabordering areas
        """
        self.intraBorderingAreas = {}
        if self.numRegionsType == "Exogenous":
            nr = range(self.pRegions)
        else:
            nr = self.feasibleRegions
        for regionID in nr:
            setNeighsNoRegion = set()
            try:
                areas2Eval = self.region2Area[regionID]
            except:
                areas2Eval = []

            for area in areas2Eval:
                setNeighsNoRegion.update(self.areas[area].neighs)
            setNeighsNoRegion.difference_update(areas2Eval)

            for neigh in setNeighsNoRegion:
                self.intraBorderingAreas.setdefault(neigh, set()).add(regionID)
    
    def kmeansInit(self, Cio):
        cachedDistances = {}
        ys = Cio
        #print('variable y: ', self.am.y)
        
        n = int(sum(value[-1] for value in self.am.y.values()))
        
        potential = [k for k,v in self.am.y.items() if int(v[-1]) == 1.0]
        #print("potencial", potential)

        #print ('Potential Nodes', potential)
        distances = np.ones(n) # n=125
        #print ('distances', distances)
        total = sum(distances) #total = 125
        #print( 'TOTAL', total)
        probabilities = list(map(lambda x: x / float(total), distances))
        seeds = []
        localDistanceType = self.distanceType  # ESTO CMABIARA
        
        returnDistance2Area = AreaCl.returnDistance2Area
        np.random.seed(int(time() * getpid()) % 4294967295)
        for k in range(self.pRegions): #pRegions = 2
            
            random = np.random.uniform(0, 1)
            find = False
            acum = 0
            cont = 0
            while not find: 
                inf = acum
                sup = acum + probabilities[cont] #0.08
                if inf <= random <= sup:
                    find = True
                    seeds += [potential[cont]]
                    #print ('SSSSSeeds', seeds) 
                    #print ('SelfAMAreas Old', self.am.areas)
                    selfAmAreas = { key: self.am.areas[key] for key in potential }
                    #print ('SelfAMAreas new', selfAmAreas)
                    #selfAmAreas = self.am.areas #Estas

                    for area in selfAmAreas.keys():
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
                                newDist = Cio.get(k, float('inf'))
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
        #print ('Potential Regions 4 Area', self.potentialRegions4Area.keys())
        #print ('Potential Regions 4 Area', self.potentialRegions4Area.values())
        
        for areaID in self.potentialRegions4Area.keys(): # loop of Neighs from seeds
            #print("areaID", areaID)
            #print("self.areas", self.areas)
            if len(self.areas)<areaID:
                print ('Problemas con el tamaño', len(self.areas), areaID)
            a = self.areas[areaID]
            regionIDs = list(self.potentialRegions4Area[areaID])
            #print('RegionIDs: ', regionIDs)
            
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

    def assignAreaStep1(self, areaID, regionID):
        """
        Assign an area to a region
        """
        a = self.areas[areaID]
        
        # print('Arguments in self.areas[areaID]')
        # for attr_name in dir(a):
        #     if not attr_name.startswith("__"):
        #         try:
        #             attr_value = getattr(a, attr_name)
        #             print(f"Attribute {attr_name}: {attr_value}", '\n')
        #         except AttributeError:
        #             print(f"Attribute {attr_name} could not be accessed.")
        
        neighs = a.neighs

        try:
            self.region2Area[regionID].append(areaID)
        except:
            self.region2Area[regionID] = [areaID]
        
        self.area2Region[areaID] = regionID

        # All 125 regions in the first loop
        #print('len of unassignedAreas: ', len(self.unassignedAreas))
        try:
            aid = self.unassignedAreas.remove(areaID)
        except Exception as e:
            print(f'Error removing {areaID} from unassignedAreas: {e}')
            pass

        # ====================================================================
        self.assignedAreas.append(areaID)
        setNeighs = set(neighs)
        setAssigned = set(self.assignedAreas)
        self.oldExternal = self.externalNeighs
        self.externalNeighs.update(setNeighs)
        self.externalNeighs.difference_update(setAssigned)
        self.newExternal = self.externalNeighs - self.oldExternal
        self.neighsMinusAssigned = setNeighs - setAssigned
        # ====================================================================
        
    def assignSeeds(self, areaID, regionID):
        """
        Assign an area to a region and updates potential regions for the neighs
        Parameters 
        """
        # Entro aquí
        self.assignAreaStep1(areaID, regionID)
        for neigh in self.neighsMinusAssigned:
            self.potentialRegions4Area.setdefault(neigh, set()).add(regionID)

        try:
            self.potentialRegions4Area.pop(areaID)
        except:
            pass
        self.changedRegion = 'null'
        self.newExternal = self.potentialRegions4Area.keys()

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
            except Exception as e:
                print(f'Could not remove {areaID} from unassignedAreas: {e}')
                pass
            self.assignedAreas.append(areaID)
            setAssigned = set(self.assignedAreas)
        nr = nr - 1

    def setSeeds(self, seeds, c=0):
        '''
        Sets the initial seeds for clustering
        '''

        if self.numRegionsType == "Exogenous" and len(seeds) <= self.pRegions:
            idx = list(range(self.n))
            #print('self.am.noNeighs', self.am.noNeighs, "\n")
            didx = list(set(idx) - (self.am.noNeighs - set(seeds))) 
            #print('didx: ', didx)

            np.random.shuffle(didx)
            #print('seeds in set seeds:', seeds + didx[0:(self.pRegions - len(seeds))])
            self.seeds = seeds + didx[0:(self.pRegions - len(seeds))]
        else:
            self.seeds = seeds

        for seed in self.seeds:
            # Si entra en el bucle
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

    def assignArea(self, areaID, regionID):
        """
        Assign an area to a region and updates potential regions for neighs
        """
        self.changedRegion = regionID
        self.addedArea = areaID
        self.assignAreaStep1(areaID, regionID)
        for neigh in self.neighsMinusAssigned:
            self.potentialRegions4Area.setdefault(neigh, set()).add(regionID)

        try:
            self.potentialRegions4Area.pop(areaID)
        except:
            pass
    
    def getObj(self, Dij = {}, Cio = {}, seeds = []):
        """
        Return the value of the objective function
        """
        #print 'getObj'
        
        if self.objInfo < 0:
            self.calcObj(Dij = Dij , Cio = Cio , seeds = seeds)
        return self.objInfo
    
    def calcObj(self, Dij = {}, Cio = {}, seeds = []):
        """
        Calculate the value of the objective function
        """
        #print 'calcObj'
        #print 'DIJ', Dij
        #print 'calcObj', self.getObjective(self.region2Area, Dij = Dij , Cio = Cio , seeds = seeds )
        self.objInfo, self.seeds = self.getObjective(self.region2Area, Dij = Dij , Cio = Cio , seeds = seeds )
        #self.objInfo = self.getObjective(self.region2Area, Dij = Dij , Cio = Cio , seeds = seeds )
    
    def getObjective(self, region2AreaDict, Dij = {}, Cio = {}, seeds = []):
        """
        Return the value of the objective function from regions2area dictionary

        This function acts as a proxy function since the idea behind the
        getObjective and getObjectiveFast is the same. When the non-fast
        approach is needed, this function will call getObjectiveFast with the
        extra parameter as None. This way the fast function will execute as the
        non-fast would have.
        """
        #print 'getObjective'
        #print 'DIJ', Dij
        #print 'OBJECTIVE', self.getObjectiveFast(region2AreaDict, Dij = Dij , Cio = Cio , seeds = seeds)
        return self.getObjectiveFast(region2AreaDict, modifiedRegions = None, Dij = Dij , Cio = Cio , seeds = seeds)

    def getObjectiveFast(self, region2AreaDict, modifiedRegions=[], Dij = {}, Cio = {}, seeds = []):
        """
        Return the value of the objective function from regions2area dictionary

        When this function gets called, the objectiveFunctionType property
        could be either a String representing the type of the function (the
        common case), or could be a list of function types, in which case it's
        necessary to iterate over all the functions.
        """
        distance = 0.0
        # equal to "Functional"
        _objFunType = self.objectiveFunctionType
        
        #tmp_obj= 9999999
        tmp_seeds = seeds
        seedsOut = seeds
        
        #print 'getObjectiveFast'
        #print 'DIJ', Dij
        if isinstance(_objFunType, "".__class__):
            if len(self.indexDataOF) == 0:
                indexData = range(len(self.areas[0].data))
            else:
                indexData = self.indexDataOF
            _fun = None
            
            if modifiedRegions == None:
                _fun = objectiveFunctionTypeDispatcher[_objFunType]
                #print 'pasa por la funcion 1'
                #distance = _fun(self, Dij = Dij , Cio = Cio , seeds = tmp_seeds)
                distance = _fun(self, region2AreaDict, Dij = Dij , Cio = Cio , seeds = tmp_seeds)
            

            else:              
                _fun = objectiveFunctionTypeDispatcher[_objFunType+'f']
                #print 'pasa por la funcion 2'
                #distance = _fun(self, Dij = Dij , Cio = Cio , seeds = tmp_seeds)
                distance=_fun(self, region2AreaDict, modifiedRegions, Dij = Dij , Cio = Cio , seeds = tmp_seeds)

            
        else:
            i = 0
            for oFT in _objFunType:
                if len(self.indexDataOF) == 0:
                    indexData = range(len(self.areas[0].data))
                else:
                    indexData = self.indexDataOF[i]

                if len(self.weightsObjectiveFunctionType) > 0:
                    _fun = objectiveFunctionTypeDispatcher[oFT]
                    distance += (self.weightsObjectiveFunctionType[i] *
                                 _fun(self, Dij = Dij , Cio = Cio , seeds = seeds))
                    i += 1
                else:
                    _fun = objectiveFunctionTypeDispatcher[oFT]
                    distance += _fun(self, Dij = Dij , Cio = Cio , seeds = seeds)
        
        #seedsOut = tmp_seeds
        #distance = tmp_obj 

        return distance, seedsOut
    
    def returnRegions(self):
        """
        Return regions created
        """
        areasId = list(self.area2Region.keys())
        #print('areasId content: ', areasId)
        areasId = np.sort(areasId).tolist()
        return [self.area2Region[area] for area in areasId]
    
    def checkFeasibility(self, regionID, areaID,
                         region2AreaDict = None):
        """
        Check feasibility from a change region (remove an area from a region)
        """
        #seeds = self.seed
        if areaID in self.seeds:
            feasible = 0
            #flag = 0
        else:            
            if not region2AreaDict:
                region2AreaDict = self.region2Area
            areas2Eval = list(region2AreaDict[regionID])
            a2r = set(region2AreaDict[regionID])
            aIDset = set([areaID])
            areas2Eval.remove(areaID)
            seedArea = areas2Eval[0]
            newRegion = set()
            newRegion.add(seedArea)
            newRegion.update(self.areas[seedArea].neighs)
            newRegion.intersection_update(areas2Eval)
            areas2Eval.remove(seedArea)
            flag = 1
            newAdded = newRegion.copy()
            newAdded.discard(seedArea)
            newNeighs = set()

            while flag:
                for area in newAdded:
                    newNeighs.update(self.areas[area].neighs)
                    areas2Eval.remove(area)
                newNeighs.intersection_update(a2r)
                newNeighs.difference_update(aIDset, newRegion, newAdded)
                newAdded = newNeighs.copy()
                newRegion.update(newAdded)

                ###
                # for area in newAdded:
                #     newNeighs = newNeighs | (((set(self.areas[area].neighs) & a2r) - aIDset) - newRegion)
                #     areas2Eval.remove(area)
                # newNeighs.difference_update(newAdded)
                # newAdded = newNeighs
                # newRegion.update(newAdded)
                ###

                if len(areas2Eval) == 0:
                    feasible = 1
                    flag = 0
                    break
                elif newNeighs == set() and len(areas2Eval) > 0:
                    feasible = 0
                    flag = 0
                    break

        return feasible
    
    def swapArea(self, area, newRegion, region2AreaDict, area2RegionDict):
        """
        Take an area from a region and give it to another
        """
        oldRegion = area2RegionDict[area]
        region2AreaDict[oldRegion].remove(area)
        region2AreaDict[newRegion].append(area)
        area2RegionDict[area] = newRegion
        if self.objectiveFunctionType == "GWalt":
            a = self.areas[area]
            lendata = len(a.data)
            adata0 = a.data[0]
            self.NRegion[newRegion] += adata0
            self.NRegion[oldRegion] -= adata0
            for index in range(1, lendata):
                self.data[newRegion][index - 1] += a.data[index] * adata0
                self.data[oldRegion][index - 1] -= a.data[index] * adata0
        if self.numRegionsType == "EndogenousThreshold":
            self.regionValue[newRegion] += self.areas[area].thresholdVar
            self.regionValue[oldRegion] -= self.areas[area].thresholdVar

    def updateTabuList(self, newValue, aList, endInd):
        """
        Add a new value to the tabu list.
        """
        return [newValue] + aList[0:endInd-1]

    def tabuMove(self, tabuLength = 5, convTabu = 5, typeTabu="random", Dij = {}, Cio = {}, seeds = []):
        """
        Conduct a solution to the best posible with tabu search
        """
        is_exact_type = (typeTabu == "exact")
        is_rand_type = (typeTabu == "random")
        aspireOBJ = self.objInfo
        currentOBJ = self.objInfo
        aspireSeeds = deepcopy(self.seeds)
        aspireRegions = self.returnRegions()
        print('CurrentRegions: ', aspireRegions)
        currentRegions = aspireRegions
        currentSeeds = deepcopy(self.seeds)
        region2AreaAspire = deepcopy(self.region2Area)
        print('Region2Area: ', self.region2Area)
        area2RegionAspire = deepcopy(self.area2Region)
        bestAdmisable = 9999999.0
        tabuList = [0]*tabuLength
        cBreak = []
        
        print('Función a minimizar: ', aspireOBJ)
        #areas4sketch = len(deepcopy(self.areas))
        admitedAreas = 50
        #print 'AREAS', areas4sketch
       
        c = 1
        self.round = 0
        resList = []
        epsilon = 1e-10
        r=1
        while c <= convTabu:
           
            if is_exact_type:
                self.objDict = makeObjDict(self)
                #print('pasa makeobjDict')
                self.allCandidates(Dij = Dij, Cio = Cio, seeds = seeds)
            else:
                moves = self.allMoves()

            if ((is_exact_type and len(self.neighSolutions) == 0) or
                (is_rand_type and len(moves) == 0)):
                c += convTabu
            else:
                if is_exact_type:
                    sortedk = sorted(self.neighSolutions)
                    end = len(sortedk)
                else:
                    end = len(moves)
                run = 0

                while run < end:
                    #print ('----------------------------------------------------------------')
                    if is_exact_type:
                        move = sortedk[run]
                    #    print('Movimientos: ', move, '\n')
                        area, region, seedsOut = move
                        obj4Move = self.neighSolutions[move]
                        candidate = 1
                    else:
                        candidate = 0
                        region2AreaCopy = deepcopy(self.region2Area)
                        area2RegionCopy = deepcopy(self.area2Region)
                        seedsCopy = deepcopy(self.seeds)
                        while (candidate == 0 and len(moves) > 0):
                            move = moves[np.random.randint(0, len(moves))]
                            moves.remove(move)
                            area, region, seeds = move
                            run += 1
                            regionIn = self.area2Region[area]
                            _feasible = self.checkFeasibility(regionIn, area) 
                            if _feasible == 1:
                                
                                #print 'ANTES DEL SWAP'
                                #obj4Move = self.recalcObj(region2AreaCopy, Dij = Dij, Cio = Cio, seeds = seeds)
                                #print obj4Move
                                #print region2AreaCopy
                                #print'SELF'
                                #print self.region2Area
                                #print self.objInfo
                       
                                self.swapArea(area, region, region2AreaCopy,
                                              area2RegionCopy)

                                if self.numRegionsType == "Exogenous":
                                    obj4Move, seedsOut = self.recalcObj(region2AreaCopy, Dij = Dij, Cio = Cio, seeds = seeds)
                                    #obj4Move = self.recalcObj(region2AreaCopy, Dij = Dij, Cio = Cio, seeds = seeds)
                                    
                                    #print 'DESPUES DEL SWAP'
                                    #print obj4Move
                                    #print region2AreaCopy
                                    #print'SELF'
                                    #print self.region2Area
                                    #print self.objInfo
                                    
                                    #print "esto es lo que evaluo", area, region, obj4Move 
                                    candidate = 1
                                elif (self.numRegionsType == "EndogenousThreshold" and
                                      self.regionValue[region] >= self.regionalThreshold and
                                      self.regionValue[regionIn] >= self.regionalThreshold):
                                    obj4Move, seedsOut = self.recalcObj(region2AreaCopy, Dij = Dij, Cio = Cio, seeds = seeds)
                                    #obj4Move = self.recalcObj(region2AreaCopy, Dij = Dij, Cio = Cio, seeds = seeds)
                                    candidate = 1
                                self.swapArea(area, regionIn, region2AreaCopy,
                                              area2RegionCopy)

                    tabuCount = 0
                    if candidate == 0:
                        c += convTabu
                        continue

                    if move in tabuList:
                        if (aspireOBJ - obj4Move) > epsilon:
                            # print( "--------CASO 1----------")
                            # print( "SI mejora, ACEPTA movimiento")
                            oldRegion = self.area2Region[area]
                            tabuList = self.updateTabuList((area, oldRegion, seedsOut),
                                                           tabuList, tabuLength)                   
                            self.moveArea(area, region, Dij = Dij, Cio = Cio, seeds = seedsOut)
                            #tabuList = self.updateTabuList((area, oldRegion, seeds),
                            #                               tabuList, tabuLength)                   
                            #self.moveArea(area, region, Dij = Dij, Cio = Cio, seeds = seeds)
                            self.objInfo = obj4Move
                            self.seeds = seedsOut
                            aspireOBJ = obj4Move                         
                            currentOBJ = obj4Move                          
                            aspireRegions = self.returnRegions()
                            currentRegions = aspireRegions
                            currentSeeds = seedsOut
                            currentSeeds = seeds
                            region2AreaAspire = deepcopy(self.region2Area)
                            area2RegionAspire = deepcopy(self.area2Region)
                            aspireSeeds = deepcopy(self.seeds)
                            bestAdmisable = obj4Move
                            cBreak.append(c)
                            c = 1
                            run = end
                            resList.append([obj4Move, aspireOBJ])
                            
                            #print '------ verificacion1 --------'
                            #print 'OF', aspireOBJ, currentOBJ
                            #print 'Movement', area, region
                            #print '-------', range(25)
                            #print 'Aspire-', aspireRegions
                            #print 'Current', currentRegions
                            #print 'Seeds', currentSeeds
                            
                            # if areas4sketch < admitedAreas:
                            #     rows = 5
                            #     cols = 5
                            #     number_areas = rows * cols
                            #     numRegions = 5
                            #     #instance = [clusterpy.createGrid(rows, cols),clusterpy.createGrid(rows, cols)]
                            #     NDs = [aspireSeeds, currentSeeds]#[5,8,17] #nodos dinamizadores
                            #     regions = [aspireRegions, currentRegions] #[2,2,2,0,0,2,2,2,0,0,2,2,1,0,0,1,1,1,1,0,1,1,1,1,1] #solución
                            #     OF = [aspireOBJ, currentOBJ]#4300 #Función objetivo
                            #     it = r #número de la iteración que se grafica (este número se usa para nombrar la figura que se exporta al folder outputMaps)


                            #     #self.drawMaps(it, instance, regions, NDs, OF)
                            #     r +=1

                        else:
                            # print( "--------CASO 2----------")
                            # print ("NO mejora, RECHAZA movimiento")
                            run += 1
                            tabuCount += 1
                            tabuList = self.updateTabuList((-1, 0, seedsOut),
                                                           tabuList, tabuLength)
                            #tabuList = self.updateTabuList((-1, 0, seeds),
                            #                               tabuList, tabuLength)
                            
                            
                            #print 'OF', aspireOBJ, currentOBJ
                            #print 'Movement', area, region
                            #print '-------', range(25)
                            #print 'Aspire-', aspireRegions
                            #print 'Current', currentRegions
                            #print 'Seeds', currentSeeds
                            
                            
                            # if areas4sketch < admitedAreas:
                            #     rows = 5
                            #     cols = 5
                            #     number_areas = rows * cols
                            #     numRegions = 5
                            #     #instance = [clusterpy.createGrid(rows, cols),clusterpy.createGrid(rows, cols)]
                            #     NDs = [aspireSeeds, currentSeeds]#[5,8,17] #nodos dinamizadores
                            #     regions = [aspireRegions, currentRegions] #[2,2,2,0,0,2,2,2,0,0,2,2,1,0,0,1,1,1,1,0,1,1,1,1,1] #solución
                            #     OF = [aspireOBJ, currentOBJ]#4300 #Función objetivo
                            #     it = r #número de la iteración que se grafica (este número se usa para nombrar la figura que se exporta al folder outputMaps)


                            #     #self.drawMaps(it, instance, regions, NDs, OF)
                            #     r +=1
                            
                            
                            
                            if tabuCount == end:
                                c = convTabu
                    else:
                        oldRegion = self.area2Region[area]
                        tabuList = self.updateTabuList((area, oldRegion, seedsOut),
                                                       tabuList, tabuLength)
                        #tabuList = self.updateTabuList((area, oldRegion, seeds),
                                                       #tabuList, tabuLength)
                        self.moveArea(area, region, Dij = Dij, Cio = Cio, seeds = seedsOut)
                        #self.moveArea(area, region, Dij = Dij, Cio = Cio, seeds = seeds)
                        self.objInfo = obj4Move
                        #self.seeds = seedsOut
                        currentOBJ = obj4Move
                        currentSeeds = seedsOut 
                        currentSeeds = seeds                      
                        if (aspireOBJ - obj4Move) > epsilon:
                            # print ("--------CASO 3----------")
                            # print ("SI mejora, ACEPTA movimiento")
                            aspireOBJ = obj4Move
                            aspireRegions = self.returnRegions()
                            currentRegions = self.returnRegions()
                            self.seeds = seedsOut
                            region2AreaAspire = deepcopy(self.region2Area)
                            area2RegionAspire = deepcopy(self.area2Region)
                            aspireSeeds = deepcopy(self.seeds)
                            cBreak.append(c)
                            c = 1
                            #print "mejora aspiracional"
                        else:
                            # print ("--------CASO 4----------")
                            # print ("NO mejora, ACEPTA movimiento")
                            currentRegions = self.returnRegions()
                            c += 1
                        bestAdmisable = obj4Move
                        run = end
                        resList.append([obj4Move, aspireOBJ])
                        
                        #print '------ verificacion2 --------'
                        #print 'OF', aspireOBJ, currentOBJ, 
                        #print 'Movement', area, region
                        #print '-------', range(25)
                        #print 'Aspire-', aspireRegions
                        #print 'Current', currentRegions
                        #print 'Seeds', currentSeeds
                        
                        # if areas4sketch < admitedAreas:
                        #     rows = 5
                        #     cols = 5
                        #     number_areas = rows * cols
                        #     numRegions = 5
                        #     #instance = [clusterpy.createGrid(rows, cols),clusterpy.createGrid(rows, cols)]
                        #     NDs = [aspireSeeds, currentSeeds]#[5,8,17] #nodos dinamizadores
                        #     regions = [aspireRegions, currentRegions] #[2,2,2,0,0,2,2,2,0,0,2,2,1,0,0,1,1,1,1,0,1,1,1,1,1] #solución
                        #     OF = [aspireOBJ, currentOBJ]#4300 #Función objetivo
                        #     it = r #número de la iteración que se grafica (este número se usa para nombrar la figura que se exporta al folder outputMaps)


                        #     #self.drawMaps(it, instance, regions, NDs, OF)
                        #     r +=1

        self.objInfo = aspireOBJ
        self.regions = aspireRegions
        self.region2Area = deepcopy(region2AreaAspire)
        self.area2Region = deepcopy(area2RegionAspire)
        self.seeds = deepcopy(aspireSeeds)
        self.resList = resList
        self.cBreak = cBreak
        
        aspireSeeds =  None
        currentSeeds = None
        region2AreaAspire = None
        area2RegionAspire = None
        gc.collect()

    def recalcObj(self, region2AreaDict, modifiedRegions=[], Dij = {}, Cio = {}, seeds = []):
        """
        Re-calculate the value of the objective function
        """
        #print 'recalcObj'
        if hasattr(self, 'objDict'):
            #print 'recalcObj', self.getObjectiveFast( region2AreaDict, Dij = Dij , Cio = Cio , seeds = seeds,)
            obj, seedsOut = self.getObjectiveFast( region2AreaDict, modifiedRegions, Dij = Dij , Cio = Cio , seeds = seeds,)
            #obj = self.getObjectiveFast( region2AreaDict, modifiedRegions, Dij = Dij , Cio = Cio , seeds = seeds,)
            #print 'Tiene el atributo'
        else:
            #print 'recalcObj', self.getObjectiveFast( region2AreaDict, Dij = Dij , Cio = Cio , seeds = seeds,)
            obj, seedsOut = self.getObjective(region2AreaDict, Dij = Dij , Cio = Cio , seeds = seeds)
            #obj = self.getObjective(region2AreaDict, Dij = Dij , Cio = Cio , seeds = seeds)
            #print 'Sin atributo'
            
        #print 'recalcObj', obj, seedsOut
        #print '+++++++++', region2AreaDict
        #print '*********', seedsOut
        return obj, seedsOut
    
    def allMoves(self):
        """
        Select all posible moves.
        """
        moves = []
        seeds = deepcopy(self.seeds)
        for area in self.intraBorderingAreas:
            regionIn = self.area2Region[area]
            regions4Move = list(self.intraBorderingAreas[area])
            if len(self.region2Area[regionIn]) > 1:
                for region in regions4Move:
                    moves.append((area, region, seeds))
        return moves

    def moveArea(self, areaID, regionID, Dij = {}, Cio = {}, seeds = []):
        """
        Move an area to a region
        """
        oldRegion = self.area2Region[areaID]
        self.region2Area[oldRegion].remove(areaID)
        self.region2Area[regionID].append(areaID)
        self.area2Region[areaID] = regionID
        a = self.areas[areaID]
        toUpdate = [areaID] + a.neighs
        if self.objectiveFunctionType == "GWalt":
            self.NRegion[regionID] += a.data[0]
            self.NRegion[oldRegion] -= a.data[0]
        if self.numRegionsType == "EndogenousThreshold":
            self.regionValue[regionID] += self.areas[areaID].thresholdVar
            self.regionValue[oldRegion] -= self.areas[areaID].thresholdVar
        try:
            for index in range(1, len(a.data)):
                self.data[regionID][index - 1] += a.data[index] * a.data[0]
            #for index in range(1, len(a.data)):
                self.data[oldRegion][index - 1] -= a.data[index] *a.data[0]
        except:
            pass
        for area in toUpdate:
            if area not in seeds:
                regionIn = self.area2Region[area]
                areasIdsIn = self.region2Area[regionIn]
                areasInNow = [self.areas[aID] for aID in areasIdsIn]
                areasInRegion = set(areasIdsIn)
                aNeighs = set(self.areas[area].neighs)
                neighsInOther = aNeighs - areasInRegion # neighs of this area in other regions
                if len(neighsInOther) == 0 and area in self.intraBorderingAreas:
                    self.intraBorderingAreas.pop(area)
                else:
                    borderRegions = set([])
                    for neigh in neighsInOther:
                        borderRegions = borderRegions | set([self.area2Region[neigh]])
                    if area in self.intraBorderingAreas:
                        self.intraBorderingAreas.pop(area)
                    self.intraBorderingAreas[area] = borderRegions
                    #_tmp = set()
                    #for neig in self.areas[area].neighs:
                    #    _tmp.add(self.area2Region[neig])

                    #_tmp.discard(self.area2Region[area])
                    #self.intraBorderingAreas[area] = _tmp
        self.calcObj(Dij = Dij, Cio = Cio, seeds = seeds)

    def allCandidates(self, Dij = {}, Cio = {}, seeds = []):
        """
        Select neighboring solutions.
        """
        intraCopy = deepcopy(self.intraBorderingAreas)
        reg2AreaCp = deepcopy(self.region2Area)
        area2RegionCp = deepcopy(self.area2Region)
        seedsCopy = deepcopy(self.seeds)
        neighSolutions = {}

        for area in intraCopy.keys():
            regionIn = self.area2Region[area]
            regions4Move = list(self.intraBorderingAreas[area])
            if (len(self.region2Area[regionIn]) > 1):
                for region in regions4Move:
                    _feasible = self.checkFeasibility(regionIn, area)
                    if _feasible == 1:
                        self.swapArea(area, region, reg2AreaCp, area2RegionCp)
                        if self.numRegionsType == "Exogenous":
                            modifiedRegions = [region, regionIn]
                            obj, seedsOut = self.recalcObj(reg2AreaCp, Dij = Dij, Cio = Cio, seeds = seedsCopy)
                            #print seedsOut
                            neighSolutions[(area, region, tuple(seedsCopy))] = obj
                        elif (self.numRegionsType == "EndogenousThreshold" and
                              self.regionValue[region] >= self.regionalThreshold and
                              self.regionValue[regionIn] >= self.regionalThreshold):
                            obj, seedsOut = self.recalcObj(reg2AreaCp, Dij = Dij, Cio = Cio, seeds = seedsCopy)
                            neighSolutions[(area, region, tuple(seedsCopy))] = obj
                        self.swapArea(area, regionIn, reg2AreaCp, area2RegionCp)

        self.neighSolutions = neighSolutions