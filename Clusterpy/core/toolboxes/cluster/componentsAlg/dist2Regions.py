import numpy as np

def getDistance2RegionCentroid(areaManager, area, areaList, indexData=[]):
    """
    The distance from area "i" to the attribute centroid of region "k"
    """
    sumAttributes = np.zeros(len(area.data))
    if len(areaManager.areas[areaList[0]].data) - len(area.data) == 1:
        for aID in areaList:
            sumAttributes += np.array(areaManager.areas[aID].data[0: -1])
    else:
        for aID in areaList:
            sumAttributes += np.array(areaManager.areas[aID].data)
    centroidRegion = sumAttributes/len(areaList)
    regionDistance = sum((np.array(area.data) - centroidRegion) ** 2)
    return regionDistance

def getDistance2RegionFunctional(areaManager, area, areaList, indexData=[], Dij = {}, Cio = {}, seeds = []):
    """
    The distance from area "i" to the attribute centroid of region "k"
    """
    #print 'Dij', Dij
    #print 'Cio', Cio
    #print 'Area Data', len(area.data)
    #print area.id
    sumAttributes = np.zeros(len(area.data))
    #print 'ACA SI LLEGA'
    #print 'Seeds', seeds
    if len(areaManager.areas[areaList[0]].data) - len(area.data) == 1:
        #print meRompo
        for aID in areaList:
            sumAttributes += np.array(areaManager.areas[aID].data[0: -1])
    else:
        for aID in areaList:
            #print 'ACA SI LLEGA 1'
            #print aID
            if area.id < aID:               
                #sumAttributes += np.array(areaManager.areas[aID].data)
                sumAttributes += np.array(Dij[area.id, aID])
            elif area.id > aID:
                #print 'ACA SI LLEGA 2'
                sumAttributes += np.array(Dij[aID, area.id])
            if aID in seeds:
                sumAttributes += np.array(Cio[area.id, aID])
        #print 'Pase dist 2 Regions'
                
    
    
    functionalRegion = sumAttributes #/len(areaList)
    regionDistance = functionalRegion #sum((np.array(area.data) - functionalRegion) ** 2)
    return regionDistance

distanceStatDispatcher = {}
distanceStatDispatcher["Centroid"] = getDistance2RegionCentroid
distanceStatDispatcher["Functional"] = getDistance2RegionFunctional