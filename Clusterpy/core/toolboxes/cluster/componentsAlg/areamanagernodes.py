from Clusterpy.core.toolboxes.cluster.componentsAlg.areacl import AreaCl
from Clusterpy.core.toolboxes.cluster.componentsAlg.dist2Regions import distanceStatDispatcher

class AreaManagerNodes:
    """
    This class contains operations at areal level, including the generation of
    instances of areas, a wide range of area2area and area2region distance
    functions.
    """
    #AreaManagerNodes(w, y, Dij, Cio, distanceType)
    def __init__(self, w, y, distanceType="EuclideanSquared", variance="false", Dij={}, Cij={}):
        """
        @type w: dictionary
        @param w: With B{key} = area Id, and B{value} = list with Ids of neighbours of
        each area.

        @type y: dictionary
        @param y: With B{key} = area Id, and B{value} = list with attribute
        values.

        @type distanceType: string
        @keyword distanceType: Function to calculate the distance between areas. Default value I{distanceType = 'EuclideanSquared'}.

        @type variance: boolean
        @keyword variance: Boolean indicating if the data have variance matrix. Default value I{variance = 'false'}.
        """
        self.y = y
        self.areas = {}
        self.noNeighs = set([])
        self.variance = variance
        self.distanceType = distanceType
        self.createAreas(w, y, Dij)
        self.distanceStatDispatcher = distanceStatDispatcher

    
    def createAreas(self, w, y, Dij):
        """
        Creates instances of areas based on a sparse weights matrix (w) and a
        data array (y).
        """
        n = len(self.y)
        #self.distances = {}
        self.distances = Dij
        noNeighs = []
        for key in range(n):
            data = y[key]
            try:
                neighbours = w[key]
            except:
                neighbours = {}
                w[key] = {}
            if len(w[key]) == 0:
                self.noNeighs = self.noNeighs | set([key])
            #print 'data', data
            a = AreaCl(key, neighbours, data)
            #a = AreaCl(key, neighbours, data, self.variance)
            self.areas[key] = a
        if len(self.noNeighs) > 0:
            print("Disconnected areas neighs: ", list(self.noNeighs))


    def returnDistance2Area(self, area, otherArea):
        """
        Returns the distance between two areas
        """
        i = 0
        j = 0
        dist = 0.0
        i = area.id
        j = otherArea.id

        if i < j:
            dist = self.distances[(i, j)]
        elif i == j:
            dist = 0.0
        else:
            dist = self.distances[(j, i)]
        return dist

    def getDataAverage(self, areaList, dataIndex):
        """
        Returns the attribute centroid of a set of areas
        """
        dataAvg = len(dataIndex) * [0.0]
        for aID in areaList:
            i = 0
            for index in dataIndex:
                dataAvg[i] += self.areas[aID].data[index] /len(areaList)
                i += 1
        return dataAvg

    def getDistance2Region(self, area, areaList, distanceStat="Centroid", weights=[], indexData=[], Dij = {}, Cio = {}, seeds = []):
        """
        Returns the distance from an area to a region (defined as a list of
        area IDs)
        """
        if isinstance(distanceStat, str):
            #print 'Index Data Stat', len(indexData)
            if len(indexData) == 0:
                indexData = range(len(area.data))
            #print self.distanceStatDispatcher[distanceStat](self, area, areaList, indexData, Dij, Cio)
            return self.distanceStatDispatcher[distanceStat](self, area, areaList, indexData, Dij, Cio, seeds)
        else:
            distance = 0.0
            i = 0
            for dS in distanceStat:
                if len(indexData) == 0:
                    indexDataDS = range(len(area.data))
                else:
                    indexDataDS = indexData[i]
                if len(weights) > 0:
                    distance += weights[i]
                    self.distanceStatDispatcher[dS](self, area, areaList, indexDataDS)
                else:
                    distance += self.distanceStatDispatcher[dS](self, area, areaList, indexDataDS)
                i += 1
            return distance

    def getDistance2AreaMin(self, area, areaList):
        """
        Return the ID of the area whitin a region that is closest to an area
        outside the region
        """
        areaMin = -1;
        distanceMin = 1e300
        for aID in areaList:
            if self.distances[area.id, aID] < distanceMin:
                areaMin = aID
                distanceMin = self.distances[area.id, aID]
        return areaMin

    def checkFeasibility(self, solution):
        """
        Checks feasibility of a candidate solution
        """
        n = len(solution)
        regions = {}
        for i in range(n):
            try:
                regions[solution[i]] = regions[solution[i]] + [i]
            except:
                regions[solution[i]] = [i]
        feasible = 1
        r = len(regions)
        for i in range(r):
            if len(regions[i]) > 0:
                newRegion = set([regions[i][0]])
                areas2Eval = set([regions[i][0]])

                while(len(areas2Eval) > 0):
                    area = areas2Eval.pop()
                    areaNeighs = (set(self.areas[area].neighs) & set(regions[i]))
                    areas2Eval = areas2Eval | (areaNeighs - newRegion)
                    newRegion = newRegion | areaNeighs
                if set(regions[i]) -newRegion != set([]):
                    feasible = 0
                    break
        return feasible