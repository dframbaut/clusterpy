import numpy as np
from Clusterpy.core.toolboxes.cluster.componentsAlg.distanceFunctions import distMethods

def getObjectiveFunctionFunctional(regionMaker, region2AreaDict, Dij = {}, Cio = {}, seeds = []):
    """
    Objective function computed with the distance between all the areas
    in a region (Clique).
    """
    ofuncval = 0.0
    for region in region2AreaDict.values():
        size = len(region)
        distmatrix = np.zeros((size, size))
        #print 'Pase objFunctions'
        #print 'Dij', Dij
        #print 'Cio', Cio
        #print 'Seeds en fo', seeds
        for iti in range(size):
            areaidi = region[iti]
            #areai = np.array(regionMaker.areas[areaid].data)
            for itj in range(size):
                areaidj = region[itj]
                if areaidi < areaidj:
                    
                    #areaj = np.array(regionMaker.areas[areaid].data)
                    #distmatrix[iti][itj] = np.linalg.norm(areai - areaj)
                    #print 'Pase Aqui'
                    #print 'iti', iti
                    #print 'itj', itj
                    #print 'araid', areaidi, areaidj
                    #print 'Dij', Dij[areaidi, areaidj]
                    distmatrix[iti][itj] = np.array(Dij[areaidi, areaidj])
                    #print 'Pase Aquiii'
                #elif areaidi > areaidj:
                    #print 'Pase Aqu'
                #    distmatrix[iti][itj] = np.array(Dij[areaidj, areaidi])
                    #print 'Pase Aquuuuuu'
                if areaidi in seeds:
                    #print 'Mori Aqui'
                    distmatrix[iti][itj] = distmatrix[iti][itj] + np.array(Cio[areaidj, areaidi])

        ofuncval += distmatrix.sum()
        
    #print 'FO', ofuncval
    return ofuncval

def getObjectiveFunctionFunctionalFast(regionMaker, region2AreaDict, modifiedRegions, Dij = {}, Cio = {}, seeds = []):
    """
    Objective function computed with the distance between all the areas
    in a region (Clique).
    """
    ofuncval = 0.0
    for region in region2AreaDict.values():
        size = len(region)
        distmatrix = np.zeros((size, size))
        #print 'Pase objFunctions'
        #print 'Dij', Dij
        #print 'Cio', Cio
        #print 'Seeds', seeds
        for iti in range(size):
            areaidi = region[iti]
            #areai = np.array(regionMaker.areas[areaid].data)
            for itj in range(size):
                areaidj = region[itj]
                if areaidi < areaidj:
                    
                    #areaj = np.array(regionMaker.areas[areaid].data)
                    #distmatrix[iti][itj] = np.linalg.norm(areai - areaj)
                    #print 'Pase Aqui'
                    distmatrix[iti][itj] = np.array(Dij[areaidi, areaidj])
                    #print 'Pase Aquiii'
                #elif areaidi > areaidj:
                    #print 'Pase Aqu'
                #    distmatrix[iti][itj] = np.array(Dij[areaidj, areaidi])
                    #print 'Pase Aquuuuuu'
                if areaidi in seeds:
                    #print 'Mori Aqui'
                    distmatrix[iti][itj] = distmatrix[iti][itj] + np.array(Cio[areaidj, areaidi])

        ofuncval += distmatrix.sum()
        
    #print 'FO', ofuncval
    #print 'FAST'
    return ofuncval

objectiveFunctionTypeDispatcher = {}
# objectiveFunctionTypeDispatcher["SS"] = getObjectiveFunctionSumSquares
# objectiveFunctionTypeDispatcher["SSf"] = getObjectiveFunctionSumSquaresFast
# objectiveFunctionTypeDispatcher["complete"] = getObjectiveFunctionClique
objectiveFunctionTypeDispatcher["Functional"] = getObjectiveFunctionFunctional
objectiveFunctionTypeDispatcher["Functionalf"] = getObjectiveFunctionFunctionalFast

def makeObjDict(regionMaker, indexData=[]):
    """
    constructs a dictionary with the objective function per region
    """
    objDict = {}

    if len(regionMaker.indexDataOF) == 0:
        indexData = range(len(regionMaker.areas[0].data))
    else:
        indexData = regionMaker.indexDataOF
    for region in regionMaker.region2Area.keys():
        objDict[region] = 0.0
        areasIdsIn = regionMaker.region2Area[region]
        areasInNow = [regionMaker.areas[aID] for aID in areasIdsIn]
        dataAvg = regionMaker.am.getDataAverage(areasIdsIn, indexData)
        c = 1
        for area in areasInNow:
            areaData = []
            for index in indexData:
                areaData += [area.data[index]]
            data = [areaData] + [dataAvg]
            areaDistance = distMethods[regionMaker.distanceType](data)
            dist = areaDistance[0][0]
            objDict[region] += dist
    return objDict